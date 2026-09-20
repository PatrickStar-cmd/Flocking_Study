function results = runN20LongSweep(seeds, workerCount, kValues, profile)
%RUNN20LONGSWEEP Parallel, resumable N=20 long-run threshold sweep.

if nargin < 1
    seeds = 1:10;
end
if nargin < 2
    workerCount = 8;
end
if nargin < 3
    kValues = [1, 2, 3, 4, 5, 7, 10, 15];
end
if nargin < 4
    profile = 'n20-long';
end
seeds = seeds(:).';
validateattributes(seeds, {'numeric'}, {'vector','integer','positive'});
validateattributes(workerCount, {'numeric'}, ...
    {'scalar','integer','positive'});

root = fileparts(mfilename('fullpath'));
outputDirectory = fullfile(root, 'results');
if ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end
csvPath = fullfile(outputDirectory, ['effective-input-', profile, '.csv']);
matPath = fullfile(outputDirectory, ['effective-input-', profile, '.mat']);
if isfile(csvPath)
    results = readtable(csvPath, 'TextType', 'string');
else
    results = table();
end

pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes', workerCount);
end

base = defaultEffectiveInputConfig();
base.NAgents = 20;
base.NNeurons = 100;
base.NSteps = 8000;
base.burnIn = 6000;
base.recordStride = 100;
policies = {'random','nearest','balanced'};
modes = {'original','fixed-total'};
expectedRowsPerSeed = 1 + 6*(numel(kValues) + 1);

for seed = seeds
    if ~isempty(results) && sum(results.seed == seed) == expectedRowsPerSeed
        fprintf('Skipping complete n20-long seed=%d.\n', seed);
        continue;
    end

    jobs = cell(1, 2 + numel(kValues)*numel(policies)*numel(modes));
    jobIndex = 1;
    baseline = base;
    baseline.seed = seed;
    baseline.selectionPolicy = 'all';
    baseline.neighborCount = baseline.NAgents - 1;
    baseline.strengthMode = 'original';
    jobs{jobIndex} = baseline;

    jobIndex = jobIndex + 1;
    zero = base;
    zero.seed = seed;
    zero.selectionPolicy = 'random';
    zero.neighborCount = 0;
    zero.strengthMode = 'original';
    jobs{jobIndex} = zero;

    for k = kValues
        for policyIndex = 1:numel(policies)
            for modeIndex = 1:numel(modes)
                jobIndex = jobIndex + 1;
                cfg = base;
                cfg.seed = seed;
                cfg.neighborCount = k;
                cfg.selectionPolicy = policies{policyIndex};
                cfg.strengthMode = modes{modeIndex};
                jobs{jobIndex} = cfg;
            end
        end
    end

    fprintf('Running n20-long seed=%d (%d simulations on %d workers) ...\n', ...
        seed, numel(jobs), workerCount);
    simulations = cell(size(jobs));
    parfor idx = 1:numel(jobs)
        simulations{idx} = simulateEffectiveSocialInput(jobs{idx});
    end

    seedRows = makeRow(simulations{1}, profile);
    zeroSimulation = simulations{2};
    for policyIndex = 1:numel(policies)
        for modeIndex = 1:numel(modes)
            zeroSimulation.config.selectionPolicy = policies{policyIndex};
            zeroSimulation.config.strengthMode = modes{modeIndex};
            seedRows(end + 1) = makeRow(zeroSimulation, profile); %#ok<AGROW>
        end
    end
    for idx = 3:numel(simulations)
        seedRows(end + 1) = makeRow(simulations{idx}, profile); %#ok<AGROW>
    end

    if ~isempty(results)
        results(results.seed == seed, :) = [];
    end
    seedTable = struct2table(seedRows);
    if isempty(results)
        results = seedTable;
    else
        results = [results; seedTable]; %#ok<AGROW>
    end
    results = sortrows(results, {'seed','k','policy','strengthMode'});
    writetable(results, csvPath);
    save(matPath, 'results');
    baseRow = seedTable(seedTable.policy == "all", :);
    fprintf('Completed seed=%d: baseline GO=%.3f LO=%.3f pair=%.3f.\n', ...
        seed, baseRow.meanGO, baseRow.meanLO, baseRow.meanPairDistanceNorm);
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
