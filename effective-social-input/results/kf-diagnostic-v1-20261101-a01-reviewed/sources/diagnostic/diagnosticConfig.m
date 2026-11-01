function cfg=diagnosticConfig(seed,profile)
% Frozen developmental configuration; never modifies repository defaults.
if nargin<1, seed=2026110101; end
if nargin<2, profile='p1'; end
validateattributes(seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
cfg=defaultEffectiveInputConfig();
cfg.seed=seed; cfg.NAgents=20; cfg.NNeurons=100;
cfg.NSteps=2200; cfg.burnIn=1000; cfg.recordStride=1;
cfg.dt=.3; cfg.speedScale=.05; cfg.arenaSize=1000;
cfg.totalSocialAttraction=.24; cfg.receptiveFieldWidth=.4;
cfg.beta=1000; cfg.connectivityExponent=.5; cfg.allocentric=true;
cfg.strengthMode='original'; cfg.distanceWeightMode='equal';
cfg.scheduledOcclusion=true; cfg.selectionPolicy='fixed-id';
cfg.neighborCount=15; cfg.occlusionHiddenCount=7; cfg.occlusionHiddenIds=[];
cfg.occlusionStartStep=1401; cfg.occlusionDurationSteps=0;
cfg.predictionMode='none'; cfg.bearingMeasurementStd=.01;
cfg.bearingObservationVariance=1e-4; cfg.bearingProcessNoise=1e-4;
cfg.predictionMaxAge=120; cfg.predictionMaxVariance=1;
cfg.turnExperiment=false;
if strcmp(profile,'fast')
    cfg.NSteps=90; cfg.burnIn=10; cfg.occlusionStartStep=31;
else
    assert(strcmp(profile,'p1'),'Unknown profile.');
end
end
