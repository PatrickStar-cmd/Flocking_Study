function cfg = defaultEffectiveInputConfig()
%DEFAULTEFFECTIVEINPUTCONFIG Parameters matching the authors' demo.

cfg.NAgents = 20;
cfg.NNeurons = 100;
cfg.NSteps = 8000;
cfg.dt = 0.3;
cfg.speedScale = 0.05;
cfg.arenaSize = 1000;
cfg.globalInhibition = 0;
cfg.beta = 1000;
cfg.totalSocialAttraction = 0.24;
cfg.receptiveFieldWidth = 0.4;
cfg.connectivityExponent = 0.5;
cfg.allocentric = true;

cfg.neighborCount = cfg.NAgents - 1;
cfg.selectionPolicy = 'all';
cfg.strengthMode = 'original';
% Optional range-bearing extension, disabled by default.
cfg.distanceWeightMode = 'equal';
cfg.distanceWeightScaleFraction = 0.5;
cfg.scheduledOcclusion = false;
cfg.occlusionStartStep = 401;
cfg.occlusionDurationSteps = 100;
cfg.occlusionHiddenCount = 2;
cfg.occlusionHiddenIds = [];
cfg.predictionMode = 'none';
cfg.bearingProcessNoise = 1e-4;
cfg.bearingGaussMarkovProcessNoise = 1e-5;
cfg.bearingVelocityTau = 10;
cfg.bearingObservationVariance = 1e-4;
% Sensor noise is explicit; matched to R in new scheduled-observation runs.
cfg.bearingMeasurementStd = 0.01;
cfg.predictionMaxAge = 120;
cfg.predictionMaxVariance = 1;

% Optional controlled-turn intervention for paired occlusion experiments.
% Leaders receive a weak rotating allocentric cue; disabled by default.
cfg.turnExperiment = false;
cfg.turnLeaderIds = [];
cfg.turnCueStartStep = 1;
cfg.turnCueDurationSteps = 0;
cfg.turnCueRate = 0;
cfg.turnCueAmplitude = 0;
cfg.turnCueBaseAngles = [];

% Bearing-only observation model. The ideal mode reproduces the paper;
% bearing-fov adds a heading-relative field of view and optional missed
% detections, without giving the controller a range measurement.
cfg.visibilityMode = 'ideal';
cfg.fieldOfView = 2*pi;
cfg.detectionProbability = 1;
cfg.occlusionEnabled = false;
cfg.agentRadius = 0;

cfg.seed = 1;
cfg.initialXY = [];
cfg.initialHeadings = [];
cfg.initialNeuralState = [];
cfg.initialPreferredDirections = [];

cfg.burnIn = 6000;
cfg.recordStride = 10;
cfg.storeMeanAdjacency = true;
end
