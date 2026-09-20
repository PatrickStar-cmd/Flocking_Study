function results = runEffectiveInputControls(seeds)
%RUNEFFECTIVEINPUTCONTROLS Paired long-run baseline and zero-input controls.
%
% The controls use the neural-field parameters from the authors' example
% (N=20, Ns=100, T=8000, burn-in=6000, h_t=0.24). Each seed is run with
% the full all-neighbor social input and with social input disabled. The
% same initial condition seed is used in each pair.

if nargin < 1
    seeds = 1:10;
end
seeds = seeds(:).';
validateattributes(seeds, {'numeric'}, {'vector','integer','positive'});

root = fileparts(mfilename('fullpath'));
outputDirectory = fullfile(root, 'results');
if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end
csvPath = fullfile(outputDirectory, 'effective-input-controls.csv');
matPath = fullfile(outputDirectory, 'effective-input-controls.mat');

if isfile(csvPath)
    results = readtable(csvPath, 'TextType', 'string');
    if ~ismember('meanLO', results.Properties.VariableNames)
        % Recompute stale controls created before LO was added.
        results = table();
    end
else
    results = table();
end

base = defaultEffectiveInputConfig();
base.NAgents = 20;
base.NNeurons = 100;
base.NSteps = 8000;
base.burnIn = 6000;
base.recordStride = 100;
base.selectionPolicy = 'all';
base.neighborCount = base.NAgents - 1;
base.strengthMode = 'original';

conditions = ["all-original", "zero-input"];
for seed = seeds
    for condition = conditions
        if isempty(results)
            already = false;
        else
            already = results.seed == seed & results.condition == condition;
        end
        if any(already)
            fprintf('Skipping existing control seed=%d condition=%s\n', ...
                seed, condition);
            continue;
        end

        cfg = base;
        cfg.seed = seed;
        if condition == "zero-input"
            cfg.selectionPolicy = 'random';
            cfg.neighborCount = 0;
        end

        fprintf('Running control seed=%d condition=%s ...\n', seed, condition);
        simulation = simulateEffectiveSocialInput(cfg);
        row = makeRow(simulation, condition);
        if isempty(results)
            results = struct2table(row);
        else
            results = [results; struct2table(row)]; %#ok<AGROW>
        end
        writetable(results, csvPath);
        save(matPath, 'results');
        fprintf('  GO=%.4f pair=%.4f nearest=%.4f\n', ...
            row.meanGO, row.meanPairDistanceNorm, ...
            row.meanNearestDistanceNorm);
    end
end

results = sortrows(results, {'seed','condition'});
writetable(results, csvPath);
save(matPath, 'results');

summary = summarizeControls(results);
writetable(summary, fullfile(outputDirectory, ...
    'effective-input-controls-summary.csv'));
fprintf('\nControl summary:\n');
disp(summary);
end

function row = makeRow(simulation, condition)
cfg = simulation.config;
s = simulation.summary;
row.condition = string(condition);
row.seed = cfg.seed;
row.N = cfg.NAgents;
row.nNeurons = cfg.NNeurons;
row.nSteps = cfg.NSteps;
row.burnIn = cfg.burnIn;
row.totalSocialAttraction = cfg.totalSocialAttraction;
row.meanGO = s.globalOrder;
row.meanLO = s.localOrder;
row.meanPairDistanceNorm = s.meanPairDistanceNorm;
row.meanNearestDistanceNorm = s.meanNearestDistanceNorm;
row.meanInputAmplitude = s.meanInputAmplitude;
row.meanActiveNeighbors = s.meanActiveNeighbors;
row.meanAngularCoverage = s.meanAngularCoverage;
row.meanInputConcentration = s.meanInputConcentration;
row.meanFieldPeak = s.meanFieldPeak;
row.meanBumpConcentration = s.meanBumpConcentration;
row.meanSpeed = s.meanSpeed;
end

function summary = summarizeControls(results)
conditions = unique(results.condition, 'stable');
rows = table('Size', [numel(conditions), 14], ...
    'VariableTypes', [{'string'}, repmat({'double'}, 1, 13)], ...
    'VariableNames', {'condition','trials','meanGO','medianGO','goP05', ...
    'goP95','meanLO','medianLO','loP05','loP95', ...
    'meanPairDistanceNorm','medianPairDistanceNorm','pairP05','pairP95'});
for idx = 1:numel(conditions)
    condition = conditions(idx);
    data = results(results.condition == condition, :);
    rows.condition(idx) = condition;
    rows.trials(idx) = height(data);
    rows.meanGO(idx) = mean(data.meanGO, 'omitnan');
    rows.medianGO(idx) = median(data.meanGO, 'omitnan');
    rows.goP05(idx) = percentileLinear(data.meanGO, 5);
    rows.goP95(idx) = percentileLinear(data.meanGO, 95);
    rows.meanLO(idx) = mean(data.meanLO, 'omitnan');
    rows.medianLO(idx) = median(data.meanLO, 'omitnan');
    rows.loP05(idx) = percentileLinear(data.meanLO, 5);
    rows.loP95(idx) = percentileLinear(data.meanLO, 95);
    rows.meanPairDistanceNorm(idx) = mean(data.meanPairDistanceNorm, 'omitnan');
    rows.medianPairDistanceNorm(idx) = median(data.meanPairDistanceNorm, 'omitnan');
    rows.pairP05(idx) = percentileLinear(data.meanPairDistanceNorm, 5);
    rows.pairP95(idx) = percentileLinear(data.meanPairDistanceNorm, 95);
end
summary = rows;
end

function value = percentileLinear(values, percentile)
values = sort(values(~isnan(values)));
if isempty(values)
    value = NaN;
    return;
end
position = 1 + (numel(values)-1)*percentile/100;
lowerIndex = floor(position);
upperIndex = ceil(position);
if lowerIndex == upperIndex
    value = values(lowerIndex);
else
    fraction = position-lowerIndex;
    value = values(lowerIndex)*(1-fraction) + values(upperIndex)*fraction;
end
end
