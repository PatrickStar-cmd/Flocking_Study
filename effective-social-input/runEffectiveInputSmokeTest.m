function runEffectiveInputSmokeTest()
%RUNEFFECTIVEINPUTSMOKETEST Fast invariants and end-to-end output checks.

root = fileparts(mfilename('fullpath'));
addpath(root);

cfg = defaultEffectiveInputConfig();
cfg.NAgents = 8;
cfg.NNeurons = 24;
cfg.NSteps = 40;
cfg.burnIn = 20;
cfg.recordStride = 5;
cfg.seed = 17;
cfg.neighborCount = 3;
cfg.selectionPolicy = 'nearest';
cfg.strengthMode = 'original';

first = simulateEffectiveSocialInput(cfg);
second = simulateEffectiveSocialInput(cfg);
assert(isequal(first.trajectory.x, second.trajectory.x), ...
    'A fixed seed must reproduce the same trajectory.');

expectedOriginal = cfg.neighborCount*cfg.totalSocialAttraction/cfg.NAgents;
assert(abs(first.summary.meanInputAmplitude-expectedOriginal) < 1e-12, ...
    'Original mode input amplitude does not match k*hTotal/N.');

fixedCfg = cfg;
fixedCfg.strengthMode = 'fixed-total';
fixed = simulateEffectiveSocialInput(fixedCfg);
expectedFixed = (cfg.NAgents - 1)*cfg.totalSocialAttraction/cfg.NAgents;
assert(abs(fixed.summary.meanInputAmplitude-expectedFixed) < 1e-12, ...
    'Fixed-total mode did not preserve the reference input amplitude.');

zeroCfg = cfg;
zeroCfg.neighborCount = 0;
zero = simulateEffectiveSocialInput(zeroCfg);
assert(zero.summary.meanInputAmplitude == 0, ...
    'Zero-neighbor runs must have zero social input amplitude.');
assert(zero.summary.meanActiveNeighbors == 0, ...
    'Zero-neighbor runs must report zero active neighbors.');

allCfg = cfg;
allCfg.selectionPolicy = 'all';
allCfg.neighborCount = cfg.NAgents - 1;
allRun = simulateEffectiveSocialInput(allCfg);
expectedAll = (cfg.NAgents - 1)*cfg.totalSocialAttraction/cfg.NAgents;
assert(abs(allRun.summary.meanInputAmplitude-expectedAll) < 1e-12, ...
    'All-neighbor run does not match the authors'' amplitude convention.');

visibilityCfg = cfg;
visibilityCfg.NAgents = 2;
visibilityCfg.NNeurons = 24;
visibilityCfg.NSteps = 10;
visibilityCfg.burnIn = 5;
visibilityCfg.neighborCount = 1;
visibilityCfg.selectionPolicy = 'all';
visibilityCfg.visibilityMode = 'bearing-fov';
visibilityCfg.fieldOfView = pi/2;
visibilityCfg.detectionProbability = 1;
visibilityCfg.initialXY = [0, -10; 0, 0];
visibilityCfg.initialHeadings = [0; 0];
visibilityCfg.initialNeuralState = zeros(visibilityCfg.NNeurons, 2);
visibilityRun = simulateEffectiveSocialInput(visibilityCfg);
assert(abs(visibilityRun.summary.meanVisibleNeighbors-0.5) < 1e-12, ...
    'The two-agent FOV case must contain one directed visible edge.');
assert(abs(visibilityRun.summary.meanActiveNeighbors-0.5) < 1e-12, ...
    'Active neighbors must match directed visibility without prediction.');

missedCfg = visibilityCfg;
missedCfg.detectionProbability = 0;
missedRun = simulateEffectiveSocialInput(missedCfg);
assert(missedRun.summary.meanVisibleNeighbors == 0, ...
    'Zero detection probability must suppress all fresh observations.');
assert(missedRun.summary.meanActiveNeighbors == 0, ...
    'Missed observations must not activate the field without prediction.');

verifyAuthorsBaseline(root);

quickResults = runEffectiveInputSweep('quick');
analysis = analyzeEffectiveInputSweep(quickResults);
assert(~isempty(analysis.calibration), 'Calibration output is empty.');
assert(~isempty(analysis.groupSuccess), 'Grouped success output is empty.');
assert(all(first.summary.localOrder >= 0 & first.summary.localOrder <= 1), ...
    'Local order must remain in the normalized [0, 1] range.');

fprintf('All effective-social-input smoke tests passed.\n');
fprintf('Quick sweep rows: %d\n', height(quickResults));
fprintf('Results directory: %s\n', fullfile(root, 'results'));
end

function verifyAuthorsBaseline(root)
% Verify every saved state against the unmodified authors' implementation.

authorsRoot = fullfile(fileparts(root), 'allocentric-flocking-original');
authorsFunction = fullfile(authorsRoot, ...
    'simulateCollectiveFlockingNeuralField.m');
assert(isfile(authorsFunction), ...
    'Authors'' baseline function was not found at %s.', authorsFunction);
addpath(authorsRoot);

cfg = defaultEffectiveInputConfig();
cfg.NAgents = 5;
cfg.NNeurons = 20;
cfg.NSteps = 7;
cfg.burnIn = 3;
cfg.recordStride = 1;
cfg.neighborCount = cfg.NAgents - 1;
cfg.selectionPolicy = 'all';
cfg.seed = 29;
cfg.initialXY = reshape(mod((1:2*cfg.NAgents)*137, ...
    cfg.arenaSize), 2, cfg.NAgents);

W = buildRingConnectivity(cfg.NNeurons, cfg.connectivityExponent);
rng(cfg.seed, 'twister');
[authorsX, authorsY, authorsU, authorsHeadings] = ...
    simulateCollectiveFlockingNeuralField(cfg.NAgents, cfg.NNeurons, ...
    cfg.NSteps, cfg.dt, cfg.speedScale, cfg.arenaSize, W, ...
    cfg.globalInhibition, cfg.beta, ...
    cfg.totalSocialAttraction/cfg.NAgents, ...
    cfg.receptiveFieldWidth, cfg.initialXY, cfg.allocentric);

candidate = simulateEffectiveSocialInput(cfg);
tolerance = 1e-12;
assert(max(abs(candidate.trajectory.x-authorsX), [], 'all') < tolerance, ...
    'Position x no longer matches the authors'' all-neighbor baseline.');
assert(max(abs(candidate.trajectory.y-authorsY), [], 'all') < tolerance, ...
    'Position y no longer matches the authors'' all-neighbor baseline.');
assert(max(abs(candidate.trajectory.headings-authorsHeadings), [], 'all') ...
    < tolerance, ...
    'Heading no longer matches the authors'' all-neighbor baseline.');
assert(max(abs(candidate.finalNeuralState-authorsU(:, :, end)), [], 'all') ...
    < tolerance, ...
    'Neural state no longer matches the authors'' all-neighbor baseline.');
end
