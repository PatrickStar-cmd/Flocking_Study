function cfg = task3DefaultConfig()
%TASK3DEFAULTCONFIG Frozen Task 3 defaults (task3_plan_v1 / config schema task3_design_v1).
%
% Time unit is the model time unit T (not seconds). Every field is either
% copied from the handoff design (handoff/task3_config_v1.json) or is a new
% Task 3 field documented in the audit report. Nothing here is a physical
% constant; thresholds are design choices.

%% Base model (Section 3.1 of the Task 3 plan)
cfg.NAgents = 20;
cfg.NNeurons = 100;
cfg.dt = 0.3;
cfg.arenaSize = 1000;
cfg.beta = 1000;
cfg.speedScale = 0.05;
cfg.totalSocialAttraction = 0.24;
cfg.receptiveFieldWidth = 0.4;
cfg.neighborCount = 15;
cfg.selectionPolicy = 'fixed-id';
cfg.strengthMode = 'original';
cfg.distanceWeightMode = 'equal';
cfg.visibilityMode = 'ideal';
cfg.allocentric = true;
cfg.globalInhibition = 0;
cfg.connectivityExponent = 0.5;
cfg.recordStride = 1;
cfg.burnIn = 0;

%% Phase structure (Section 4.1)
cfg.initializationSteps = 5000;   % source-state steps, fixed, never seed-selected
cfg.commonPrefixSteps = 200;      % shared fresh-observation warmup
cfg.turnPrefixSteps = 20;         % cue starts in the last 20 prefix steps
cfg.postPrimarySteps = 100;       % main post window
cfg.recoverySteps = 300;          % total post-occlusion horizon for recovery

%% Prediction lifetime control (Section 3.2)
cfg.predictionMaxAge = 120;
cfg.predictionMaxVariance = 1e12; % compatibility bound only, not a reliability gate
cfg.assertNoVarianceCutoff = true;

%% Recovery definition (Section 7.3). Thresholds are design choices used to
%% define recovery, not stability guarantees.
cfg.recoveryGOThreshold = 0.05;
cfg.recoveryDistanceThreshold = 0.005;
cfg.recoveryConsecutiveSteps = 10;
cfg.recoveryTimeDefinition = 'first step of the first forward-confirmed window after reacquisition';

%% Region thresholds (Section 9.1). Minimum meaningful differences; they may be
%% revised only during development and must never be revised on validation data.
cfg.deltaGO = 0.01;
cfg.deltaDistance = 0.005;
cfg.deltaTurnRad = 0.05;

%% Resampling (Section 9.1). Unit is the seed cluster.
cfg.bootstrapRepeats = 5000;
cfg.bootstrapSeed = 70001;
cfg.confidenceLevel = 0.95;
cfg.resampleUnit = 'seed_cluster';
cfg.calibrationLocalStdMax = 0.5;
cfg.calibrationVarianceRatioRange = [0.5 2];
cfg.calibrationCoverageRange = [0.9 0.99];
cfg.regionLabels = {'beneficial','equivalent','harmful','tradeoff','undetermined'};
% Task 3 removes the age gate entirely (the plan's Section 3.2 freezes
% predictionMaxAge=120 for historical comparability but requires that the
% whole main grid be shown never to reach it). Branch runs set a huge value so
% that no edge is ever dropped for age, and the harness asserts afterwards
% that zero age-drop events occurred. This is what keeps the inactivated-edge
% diagnostics free of selection bias.
cfg.branchPredictionMaxAge = 1e9;
% Longest branch horizon (occlusion + recovery) the run will need. The prefix
% pre-generates the shared sensor-noise realisation to cover it, so every
% branch reads the same noise at the same absolute step.
cfg.branchHorizonSteps = 0;

%% Measurement model
cfg.bearingMeasurementStd = 0.01;
cfg.bearingObservationVariance = 0.01^2;
cfg.RPolicy = 'matched_to_measurement_std_squared';

%% Filter parameters (frozen after calibration; see calibration stage)
cfg.bearingProcessNoise = 1e-4;              % CV q, candidate until calibrated
cfg.bearingGaussMarkovProcessNoise = 1e-5;   % GM q, candidate until calibrated
cfg.bearingVelocityTau = 10;                 % GM tau, candidate until calibrated
cfg.predictionMode = 'none';
% Input-kernel construction (Section 7.2). 'wrapped' is the periodic
% wrapped-normal kernel with historical unit-peak scale, which makes
% the covariance heat convolution mass-exact and needs no clip. 'truncated'
% reproduces the historical sampled-Gaussian-plus-clip path and is retained
% only so the legacy contract test still exercises it.
cfg.occlusionKernelMode = 'wrapped';

%% Occlusion
cfg.scheduledOcclusion = true;
cfg.occlusionStartStep = 1;
cfg.occlusionDurationSteps = 0;
cfg.occlusionHiddenCount = 0;
cfg.occlusionHiddenIds = [];

%% Scenario / cue (Section 5.1)
cfg.turnExperiment = false;
cfg.turnLeaderIds = [1 2 3];
cfg.turnCueStartStep = 1;         % 1-based step within this run
cfg.turnCueDurationSteps = 0;
cfg.turnCueRate = 0;              % rad per model time unit
cfg.turnCueAmplitude = 0;
cfg.turnCueBaseAngles = [];

%% Task 3 fork indexing (Section 4.1: exact, shared observation stream)
% The sensor-noise realisation is generated once by the prefix run and carried
% through the filter state, so branches read the same absolute step. The neural
% field carries NO per-step noise (the authors perturb it only at
% initialization), so no neural stream needs to be checkpointed.
cfg.absoluteStepOffset = 0;       % steps elapsed before this run (fork continuation)

%% Storage control (keep prefix runs and short branches cheap)
cfg.storeObservationRecords = true;
cfg.storeTrajectory = true;
cfg.storeTimeSeries = true;
cfg.storeMeanAdjacency = false;
cfg.storeInputDiagnostics = false;   % extra per-step input decomposition
cfg.decompositionSteps = [];         % absolute steps at which to snapshot decomposition
% Historical runners consume observationRecords as a struct array of full
% n-by-n matrices. Task 3 uses the compact edge record by default because the
% legacy shape costs roughly 70 MB per 535-step branch. Set true only for the
% historical entry points (runPaperKFExtendedSimulation, runTurnOcclusionPilot)
% and the legacy contract tests.
cfg.legacyObservationRecord = false;

%% Identity
cfg.profile = 'smoke';
cfg.runId = '';
cfg.seed = 1;
end
