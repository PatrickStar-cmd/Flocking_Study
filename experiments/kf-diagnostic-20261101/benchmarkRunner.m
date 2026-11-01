function benchReport = benchmarkRunner(seed)
% BENCHMARKRUNNER Measure wall-clock runtime for each simulation condition.
% Runs full-length NSteps=2200 simulations for representative conditions
% and provides runtime statistics and whole-batch time projection.

if nargin < 1 || isempty(seed)
    seed = 2026110101;
end

thisDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(thisDir, '..', '..', 'effective-social-input');
addpath(thisDir);
addpath(srcDir);

fprintf('========================================================================\n');
fprintf('STARTING BENCHMARK RUNNER (Full Horizon: NSteps=2200, N=20, M=100)\n');
fprintf('Seed: %d, MATLAB: %s, Platform: %s\n', seed, version, computer);
fprintf('========================================================================\n\n');

cfgBase = defaultEffectiveInputConfig();
cfgBase.NAgents = 20;
cfgBase.NNeurons = 100;
cfgBase.NSteps = 2200;
cfgBase.burnIn = 1000;
cfgBase.dt = 0.3;
cfgBase.speedScale = 0.05;
cfgBase.arenaSize = 1000;
cfgBase.totalSocialAttraction = 0.24;
cfgBase.receptiveFieldWidth = 0.4;
cfgBase.connectivityExponent = 0.5;
cfgBase.beta = 1000;
cfgBase.allocentric = true;
cfgBase.selectionPolicy = 'fixed-id';
cfgBase.strengthMode = 'original';
cfgBase.distanceWeightMode = 'equal';
cfgBase.scheduledOcclusion = true;
cfgBase.neighborCount = 15;
cfgBase.occlusionHiddenCount = 7;
cfgBase.occlusionHiddenIds = [];
cfgBase.occlusionStartStep = 1401;
cfgBase.bearingMeasurementStd = 0.01;
cfgBase.bearingObservationVariance = 1e-4;
cfgBase.bearingProcessNoise = 1e-4;
cfgBase.predictionMaxAge = 120;
cfgBase.predictionMaxVariance = 1.0;
cfgBase.turnExperiment = false;
cfgBase.recordStride = 1;
cfgBase.seed = seed;

conditions = {
    'fresh',   0,   'none';
    'none',   20,   'none';
    'hold',   20,   'hold';
    'kf-mean',20,   'kf-mean';
    'kalman', 20,   'kalman';
    'none',  120,   'none';
    'hold',  120,   'hold';
    'kf-mean',120,  'kf-mean';
    'kalman', 120,  'kalman'
};

benchReport = struct();
benchReport.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
benchReport.seed = seed;
benchReport.runs = struct();

totalBenchTime = 0;
times = zeros(size(conditions, 1), 1);

for idx = 1:size(conditions, 1)
    condName = conditions{idx, 1};
    dur = conditions{idx, 2};
    mode = conditions{idx, 3};

    cfg = cfgBase;
    cfg.occlusionDurationSteps = dur;
    cfg.predictionMode = mode;

    label = sprintf('%s_D%d', condName, dur);
    fprintf('Benchmarking condition %d/%d: %-15s (Duration=%d, Mode=%-7s)... ', ...
        idx, size(conditions, 1), label, dur, mode);

    tStart = tic;
    res = simulateEffectiveSocialInputInstrumented(cfg);
    elapsed = toc(tStart);
    times(idx) = elapsed;
    totalBenchTime = totalBenchTime + elapsed;

    stepsPerSec = cfg.NSteps / elapsed;
    fprintf('Elapsed: %6.2f s (%6.1f steps/s)\n', elapsed, stepsPerSec);

    runInfo = struct();
    runInfo.condition = label;
    runInfo.duration = dur;
    runInfo.mode = mode;
    runInfo.elapsedSeconds = elapsed;
    runInfo.stepsPerSecond = stepsPerSec;
    benchReport.runs.(matlab.lang.makeValidName(label)) = runInfo;
end

meanTimePerRun = mean(times);
projected45BatchHours = (meanTimePerRun * 45) / 3600;

benchReport.meanTimePerRun = meanTimePerRun;
benchReport.projectedBatchTimeHours = projected45BatchHours;
benchReport.totalBenchTimeSeconds = totalBenchTime;

fprintf('\n------------------------------------------------------------------------\n');
fprintf('BENCHMARK SUMMARY:\n');
fprintf('  Average runtime per run:   %6.2f seconds\n', meanTimePerRun);
fprintf('  Projected 45-run batch:    %6.2f minutes (%5.2f hours)\n', ...
    (meanTimePerRun * 45) / 60, projected45BatchHours);
fprintf('------------------------------------------------------------------------\n');

reportPath = fullfile(thisDir, 'benchmark_report.json');
fid = fopen(reportPath, 'w');
if fid ~= -1
    fwrite(fid, jsonencode(benchReport), 'char');
    fclose(fid);
end
fprintf('Saved benchmark report to: %s\n', reportPath);
end
