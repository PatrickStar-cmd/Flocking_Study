function rows = runPersistencePilot(seeds)
%RUNPERSISTENCEPILOT Measure time-resolved success without changing dynamics.
% Uses strict thresholds frozen by independent-validation.json.
if nargin < 1
    seeds = 1:10;
end
root = fileparts(mfilename('fullpath'));
outDir = fullfile(root, 'results', 'persistence-pilot');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
thresholds = [0.3413210142885347, 0.5297546055663604, 0.24277616843099487];
conditions = {
    'all', 'original', 19;
    'random', 'fixed-total', 7;
    'random', 'fixed-total', 15;
    'balanced', 'fixed-total', 15};
rows = struct([]);
idx = 0;
for seed = seeds
    for c = 1:size(conditions, 1)
        cfg = defaultEffectiveInputConfig();
        cfg.NAgents = 20;
        cfg.NNeurons = 100;
        cfg.NSteps = 8000;
        cfg.burnIn = 6000;
        cfg.recordStride = 1;
        cfg.seed = seed;
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
        idx = idx + 1;
        rows(idx).seed = seed; %#ok<AGROW>
        rows(idx).policy = string(cfg.selectionPolicy);
        rows(idx).strengthMode = string(cfg.strengthMode);
        rows(idx).k = cfg.neighborCount;
        rows(idx).evaluationSteps = numel(success);
        rows(idx).sustainedSuccessFraction = mean(success);
        rows(idx).firstFailureStep = firstFailure;
        rows(idx).longestFailureSteps = longestFailure;
        rows(idx).failureRuns = failureCount;
        rows(idx).recoveryRuns = recoveryCount;
        rows(idx).meanGO = simulation.summary.globalOrder;
        rows(idx).meanLO = simulation.summary.localOrder;
        rows(idx).meanPairDistanceNorm = simulation.summary.meanPairDistanceNorm;
    end
end
rows = struct2table(rows);
writetable(rows, fullfile(outDir, 'persistence-pilot.csv'));
metadata.thresholds = thresholds;
metadata.seeds = seeds;
metadata.conditions = conditions;
metadata.recordStride = 1;
metadata.note = 'Strict thresholds are frozen from independent-validation.json; no threshold is re-estimated here.';
save(fullfile(outDir, 'persistence-pilot.mat'), 'rows', 'metadata');
end

function [longestFailure, failureCount, firstFailure, recoveryCount] = runLengths(success)
failure = ~success(:).';
edges = diff([false, failure, false]);
starts = find(edges == 1);
ends = find(edges == -1) - 1;
failureCount = numel(starts);
if isempty(starts)
    longestFailure = 0;
else
    longestFailure = max(ends - starts + 1);
end
firstFailure = find(failure, 1, 'first');
if isempty(firstFailure)
    firstFailure = NaN;
end
recoveryCount = sum(diff([false, success(:).']) == 1 & ...
    [false, failure(1:end-1)]);
end
