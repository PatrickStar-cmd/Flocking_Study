function results = runEffectiveInputSweep(profile)
%RUNEFFECTIVEINPUTSWEEP Run a resumable phase-one neighbor-input sweep.

if nargin < 1
    profile = 'quick';
end
profile = char(profile);

if strcmp(profile, 'n20-long')
    results = runN20LongSweep();
    return;
end
if strcmp(profile, 'n20-final')
    results = runN20LongSweep(1:30, 8, 1:18, 'n20-final');
    return;
end

cfg = defaultEffectiveInputConfig();
switch profile
    case 'quick'
        nValues = 12;
        kTemplate = [0, 1, 3, 5, 7, 10];
        seeds = 1:3;
        cfg.NNeurons = 40;
        cfg.NSteps = 400;
        cfg.burnIn = 200;
        cfg.recordStride = 20;
    case 'pilot'
        nValues = 20;
        kTemplate = [0, 1, 2, 3, 4, 5, 7, 10, 15];
        seeds = 1:10;
        cfg.NNeurons = 60;
        cfg.NSteps = 2000;
        cfg.burnIn = 1200;
        cfg.recordStride = 50;
    case 'full'
        nValues = [20, 40, 80];
        kTemplate = [0, 1, 2, 3, 4, 5, 7, 10, 15];
        seeds = 1:30;
        cfg.NNeurons = 100;
        cfg.NSteps = 8000;
        cfg.burnIn = 6000;
        cfg.recordStride = 100;
    otherwise
        error('profile must be quick, pilot, n20-long, n20-final, or full.');
end

policies = {'random','nearest','balanced'};
strengthModes = {'original','fixed-total'};
outputDirectory = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end
csvPath = fullfile(outputDirectory, ['effective-input-', profile, '.csv']);
matPath = fullfile(outputDirectory, ['effective-input-', profile, '.mat']);

rows = struct([]);
runIndex = 0;
for n = nValues
    cfg.NAgents = n;
    kValues = unique([kTemplate(kTemplate < n - 1), n - 1]);
    for seed = seeds
        % One explicit all-neighbor baseline per seed.
        baselineCfg = cfg;
        baselineCfg.seed = seed;
        baselineCfg.neighborCount = n - 1;
        baselineCfg.selectionPolicy = 'all';
        baselineCfg.strengthMode = 'original';
        baseline = simulateEffectiveSocialInput(baselineCfg);
        row = makeRow(baseline, profile);
        runIndex = runIndex + 1;
        if runIndex == 1
            rows = row;
        else
            rows(runIndex) = row; %#ok<AGROW>
        end
        fprintf('[%s] N=%d k=%d policy=all mode=original seed=%d GO=%.3f\n', ...
            profile, n, n - 1, seed, baseline.summary.globalOrder);

        for k = kValues(kValues < n - 1)
            if k == 0
                zeroCfg = cfg;
                zeroCfg.seed = seed;
                zeroCfg.neighborCount = 0;
                zeroCfg.selectionPolicy = 'random';
                zeroCfg.strengthMode = 'original';
                zeroSimulation = simulateEffectiveSocialInput(zeroCfg);
            end
            for policyIndex = 1:numel(policies)
                for modeIndex = 1:numel(strengthModes)
                    runCfg = cfg;
                    runCfg.seed = seed;
                    runCfg.neighborCount = k;
                    runCfg.selectionPolicy = policies{policyIndex};
                    runCfg.strengthMode = strengthModes{modeIndex};
                    if k == 0
                        simulation = zeroSimulation;
                        simulation.config.selectionPolicy = policies{policyIndex};
                        simulation.config.strengthMode = strengthModes{modeIndex};
                    else
                        simulation = simulateEffectiveSocialInput(runCfg);
                    end
                    row = makeRow(simulation, profile);
                    runIndex = runIndex + 1;
                    rows(runIndex) = row; %#ok<AGROW>
                    fprintf(['[%s] N=%d k=%d policy=%s mode=%s ', ...
                        'seed=%d GO=%.3f\n'], profile, n, k, ...
                        policies{policyIndex}, strengthModes{modeIndex}, ...
                        seed, simulation.summary.globalOrder);
                end
            end
        end

        results = struct2table(rows);
        writetable(results, csvPath);
        save(matPath, 'results');
    end
end
end

function row = makeRow(simulation, profile)
cfg = simulation.config;
s = simulation.summary;
row.profile = string(profile);
row.N = cfg.NAgents;
row.nNeurons = cfg.NNeurons;
row.nSteps = cfg.NSteps;
row.seed = cfg.seed;
row.k = cfg.neighborCount;
row.policy = string(cfg.selectionPolicy);
row.strengthMode = string(cfg.strengthMode);
row.meanGO = s.globalOrder;
row.meanLO = s.localOrder;
row.meanPairDistanceNorm = s.meanPairDistanceNorm;
row.meanNearestDistanceNorm = s.meanNearestDistanceNorm;
row.meanInputAmplitude = s.meanInputAmplitude;
row.meanVisibleNeighbors = s.meanVisibleNeighbors;
row.meanActiveNeighbors = s.meanActiveNeighbors;
row.meanAngularCoverage = s.meanAngularCoverage;
row.meanInputConcentration = s.meanInputConcentration;
row.meanFieldPeak = s.meanFieldPeak;
row.meanBumpConcentration = s.meanBumpConcentration;
row.meanSpeed = s.meanSpeed;
row.socialUnionLccFraction = s.socialUnionLccFraction;
row.socialAlgebraicConnectivity = s.socialAlgebraicConnectivity;
end
