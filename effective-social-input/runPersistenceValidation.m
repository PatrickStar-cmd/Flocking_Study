function rows = runPersistenceValidation(seeds)
%RUNPERSISTENCEVALIDATION Independent persistence endpoint validation.
% Runs 30 fresh seeds and writes a separate, resumable result file.
if nargin < 1
    seeds = 11:40;
end
root = fileparts(mfilename('fullpath'));
outDir = fullfile(root, 'results', 'persistence-validation');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
csvPath = fullfile(outDir, 'persistence-validation.csv');
if isfile(csvPath)
    rows = readtable(csvPath, 'TextType', 'string');
else
    rows = table();
end
thresholds = [0.3413210142885347, 0.5297546055663604, 0.24277616843099487];
conditions = {
    'all', 'original', 19;
    'random', 'fixed-total', 7;
    'random', 'fixed-total', 15;
    'balanced', 'fixed-total', 15};
for seed = seeds
    if ~isempty(rows) && sum(rows.seed == seed) == size(conditions, 1)
        continue;
    end
    seedRows = struct('seed', {}, 'policy', {}, 'strengthMode', {}, 'k', {}, ...
        'evaluationSteps', {}, 'sustainedSuccessFraction', {}, ...
        'persistence90', {}, 'firstFailureStep', {}, 'longestFailureSteps', {}, ...
        'failureRuns', {}, 'recoveryRuns', {}, 'meanGO', {}, 'meanLO', {}, ...
        'meanPairDistanceNorm', {});
    for c = 1:size(conditions, 1)
        cfg = defaultEffectiveInputConfig();
        cfg.NAgents = 20; cfg.NNeurons = 100; cfg.NSteps = 8000;
        cfg.burnIn = 6000; cfg.recordStride = 1; cfg.seed = seed;
        cfg.selectionPolicy = conditions{c, 1};
        cfg.strengthMode = conditions{c, 2};
        cfg.neighborCount = conditions{c, 3};
        simulation = simulateEffectiveSocialInput(cfg);
        ts = simulation.timeSeries;
        evalIdx = (cfg.burnIn + 1):cfg.NSteps;
        success = ts.globalOrder(evalIdx) >= thresholds(1) & ...
            ts.localOrder(evalIdx) >= thresholds(2) & ...
            ts.meanPairDistanceNorm(evalIdx) <= thresholds(3);
        [longestFailure, failureCount, firstFailure, recoveryCount] = ...
            runLengths(success);
        row.seed = seed; row.policy = string(cfg.selectionPolicy);
        row.strengthMode = string(cfg.strengthMode); row.k = cfg.neighborCount;
        row.evaluationSteps = numel(success);
        row.sustainedSuccessFraction = mean(success);
        row.persistence90 = row.sustainedSuccessFraction >= 0.90;
        row.firstFailureStep = firstFailure;
        row.longestFailureSteps = longestFailure;
        row.failureRuns = failureCount; row.recoveryRuns = recoveryCount;
        row.meanGO = simulation.summary.globalOrder;
        row.meanLO = simulation.summary.localOrder;
        row.meanPairDistanceNorm = simulation.summary.meanPairDistanceNorm;
        seedRows(end + 1) = row; %#ok<AGROW>
    end
    seedTable = struct2table(seedRows);
    if isempty(rows), rows = seedTable; else, rows = [rows; seedTable]; end %#ok<AGROW>
    rows = sortrows(rows, {'seed','policy','strengthMode','k'});
    writetable(rows, csvPath);
end
metadata.thresholds = thresholds;
metadata.persistenceTarget = 0.90;
metadata.seeds = seeds;
metadata.conditions = conditions;
metadata.note = 'Strict thresholds frozen from independent-validation.json; persistence target is 90% of post-burn-in time points.';
save(fullfile(outDir, 'persistence-validation.mat'), 'rows', 'metadata');
end

function [longestFailure, failureCount, firstFailure, recoveryCount] = runLengths(success)
failure = ~success(:).';
edges = diff([false, failure, false]);
starts = find(edges == 1); ends = find(edges == -1) - 1;
failureCount = numel(starts);
longestFailure = 0;
if ~isempty(starts), longestFailure = max(ends - starts + 1); end
firstFailure = find(failure, 1, 'first');
if isempty(firstFailure), firstFailure = NaN; end
recoveryCount = sum(diff([false, success(:).']) == 1 & [false, failure(1:end-1)]);
end
