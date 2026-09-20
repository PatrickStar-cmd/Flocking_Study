function analysis = analyzeEffectiveInputSweep(inputData)
%ANALYZEEFFECTIVEINPUTSWEEP Estimate baseline-equivalent success thresholds.

if nargin < 1
    inputData = fullfile(fileparts(mfilename('fullpath')), ...
        'results', 'effective-input-quick.csv');
end
if istable(inputData)
    runs = inputData;
else
    runs = readtable(inputData, 'TextType', 'string');
end

isBaseline = runs.policy == "all" & runs.strengthMode == "original";
if ~any(isBaseline)
    error('Input data does not contain all-neighbor original baselines.');
end

nValues = unique(runs.N);
calibrationRows = struct([]);
groupRows = struct([]);
thresholdRows = struct([]);
calibrationIndex = 0;
groupIndex = 0;
thresholdIndex = 0;

for n = nValues.'
    nMask = runs.N == n;
    base = runs(nMask & isBaseline, :);
    goThreshold = percentileLinear(base.meanGO, 5);
    loThreshold = percentileLinear(base.meanLO, 5);
    pairThreshold = percentileLinear(base.meanPairDistanceNorm, 95);
    goMedian = median(base.meanGO, 'omitnan');
    loMedian = median(base.meanLO, 'omitnan');
    pairMedian = median(base.meanPairDistanceNorm, 'omitnan');
    goScale = max(goMedian-goThreshold, eps);
    loScale = max(loMedian-loThreshold, eps);
    pairScale = max(pairThreshold-pairMedian, eps);
    baseScore = min([(base.meanGO-goThreshold)/goScale, ...
        (base.meanLO-loThreshold)/loScale, ...
        (pairThreshold-base.meanPairDistanceNorm)/pairScale], [], 2);
    scoreThreshold = lowerOrderPercentile(baseScore, 5);

    calibrationIndex = calibrationIndex + 1;
    calibrationRows(calibrationIndex).N = n; %#ok<AGROW>
    calibrationRows(calibrationIndex).goThreshold = goThreshold;
    calibrationRows(calibrationIndex).loThreshold = loThreshold;
    calibrationRows(calibrationIndex).pairDistanceThreshold = pairThreshold;
    calibrationRows(calibrationIndex).goMedian = goMedian;
    calibrationRows(calibrationIndex).loMedian = loMedian;
    calibrationRows(calibrationIndex).pairDistanceMedian = pairMedian;
    calibrationRows(calibrationIndex).scoreThreshold = scoreThreshold;

    score = min([(runs.meanGO(nMask)-goThreshold)/goScale, ...
        (runs.meanLO(nMask)-loThreshold)/loScale, ...
        (pairThreshold-runs.meanPairDistanceNorm(nMask))/pairScale], [], 2);
    runs.baselineEquivalent(nMask) = score >= scoreThreshold;
    runs.baselineEquivalenceScore(nMask) = score;

    policies = unique(runs.policy(nMask));
    modes = unique(runs.strengthMode(nMask));
    for policy = policies.'
        for mode = modes.'
            pmMask = nMask & runs.policy == policy & runs.strengthMode == mode;
            if ~any(pmMask)
                continue;
            end
            kValues = unique(runs.k(pmMask));
            qualifyingK = [];
            qualifyingWilsonK = [];
            currentGroupIndices = zeros(numel(kValues), 1);
            currentRates = zeros(numel(kValues), 1);
            currentTrials = zeros(numel(kValues), 1);
            kIndex = 0;
            for k = kValues.'
                kIndex = kIndex + 1;
                mask = pmMask & runs.k == k;
                successes = sum(runs.baselineEquivalent(mask));
                trials = sum(mask);
                rate = successes/trials;
                [lower, upper] = wilsonInterval(successes, trials, 1.96);
                groupIndex = groupIndex + 1;
                groupRows(groupIndex).N = n; %#ok<AGROW>
                groupRows(groupIndex).policy = policy;
                groupRows(groupIndex).strengthMode = mode;
                groupRows(groupIndex).k = k;
                groupRows(groupIndex).trials = trials;
                groupRows(groupIndex).successes = successes;
                groupRows(groupIndex).successRate = rate;
                groupRows(groupIndex).wilsonLower95 = lower;
                groupRows(groupIndex).wilsonUpper95 = upper;
                groupRows(groupIndex).isotonicSuccessRate = NaN;
                currentGroupIndices(kIndex) = groupIndex;
                currentRates(kIndex) = rate;
                currentTrials(kIndex) = trials;
                if rate >= 0.95
                    qualifyingK(end + 1) = k; %#ok<AGROW>
                end
                if lower >= 0.95
                    qualifyingWilsonK(end + 1) = k; %#ok<AGROW>
                end
            end
            isotonicRates = isotonicNondecreasing(currentRates, currentTrials);
            for idx = 1:numel(currentGroupIndices)
                groupRows(currentGroupIndices(idx)).isotonicSuccessRate = ...
                    isotonicRates(idx);
            end
            qualifyingIsotonicK = kValues(isotonicRates >= 0.95);
            thresholdIndex = thresholdIndex + 1;
            thresholdRows(thresholdIndex).N = n; %#ok<AGROW>
            thresholdRows(thresholdIndex).policy = policy;
            thresholdRows(thresholdIndex).strengthMode = mode;
            if isempty(qualifyingK)
                thresholdRows(thresholdIndex).kStarEmpirical = NaN;
            else
                thresholdRows(thresholdIndex).kStarEmpirical = min(qualifyingK);
            end
            if isempty(qualifyingIsotonicK)
                thresholdRows(thresholdIndex).kStarIsotonic = NaN;
            else
                thresholdRows(thresholdIndex).kStarIsotonic = ...
                    min(qualifyingIsotonicK);
            end
            if isempty(qualifyingWilsonK)
                thresholdRows(thresholdIndex).kStarWilsonLower95 = NaN;
            else
                thresholdRows(thresholdIndex).kStarWilsonLower95 = ...
                    min(qualifyingWilsonK);
            end
        end
    end
end

analysis.runs = runs;
analysis.calibration = struct2table(calibrationRows);
analysis.groupSuccess = struct2table(groupRows);
analysis.thresholds = struct2table(thresholdRows);

outputDirectory = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end
writetable(analysis.calibration, fullfile(outputDirectory, 'calibration.csv'));
writetable(analysis.groupSuccess, fullfile(outputDirectory, 'group-success.csv'));
writetable(analysis.thresholds, fullfile(outputDirectory, 'thresholds.csv'));
end

function fitted = isotonicNondecreasing(values, weights)
% Weighted pool-adjacent-violators estimate of a nondecreasing sequence.
blockValue = values(:);
blockWeight = weights(:);
blockStart = (1:numel(values)).';
blockEnd = blockStart;
nBlocks = numel(values);
idx = 1;
while idx < nBlocks
    if blockValue(idx) <= blockValue(idx + 1)
        idx = idx + 1;
        continue;
    end
    combinedWeight = blockWeight(idx) + blockWeight(idx + 1);
    blockValue(idx) = (blockValue(idx)*blockWeight(idx) + ...
        blockValue(idx + 1)*blockWeight(idx + 1))/combinedWeight;
    blockWeight(idx) = combinedWeight;
    blockEnd(idx) = blockEnd(idx + 1);
    blockValue(idx + 1) = [];
    blockWeight(idx + 1) = [];
    blockStart(idx + 1) = [];
    blockEnd(idx + 1) = [];
    nBlocks = nBlocks - 1;
    idx = max(1, idx - 1);
end
fitted = zeros(numel(values), 1);
for idx = 1:nBlocks
    fitted(blockStart(idx):blockEnd(idx)) = blockValue(idx);
end
end

function value = percentileLinear(values, percentile)
values = sort(values(~isnan(values)));
if isempty(values)
    value = NaN;
    return;
end
position = 1 + (numel(values) - 1)*percentile/100;
lowerIndex = floor(position);
upperIndex = ceil(position);
if lowerIndex == upperIndex
    value = values(lowerIndex);
else
    fraction = position-lowerIndex;
    value = values(lowerIndex)*(1-fraction) + values(upperIndex)*fraction;
end
end

function value = lowerOrderPercentile(values, percentile)
% Select an observed lower-tail value so at least 1-p of samples pass.
values = sort(values(~isnan(values)));
if isempty(values)
    value = NaN;
    return;
end
index = max(1, ceil(numel(values)*percentile/100));
value = values(index);
end

function [lower, upper] = wilsonInterval(successes, trials, z)
if trials == 0
    lower = NaN;
    upper = NaN;
    return;
end
rate = successes/trials;
denominator = 1 + z^2/trials;
center = (rate + z^2/(2*trials))/denominator;
halfWidth = z*sqrt(rate*(1-rate)/trials + z^2/(4*trials^2))/denominator;
lower = max(0, center-halfWidth);
upper = min(1, center+halfWidth);
end
