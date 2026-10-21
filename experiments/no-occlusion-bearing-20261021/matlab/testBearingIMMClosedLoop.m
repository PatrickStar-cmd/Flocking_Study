function audit=testBearingIMMClosedLoop()
cfg=defaultEffectiveInputConfig(); cfg.NSteps=180; cfg.burnIn=80;
cfg.seed=6299; cfg.neighborCount=15; cfg.selectionPolicy='random';
cfg.recordStride=1; cfg.storeTimeSeries=true; cfg.storeTrajectory=true;
cfg.pilotSampleStride=20; cfg.bearingProcessNoise=.1;
cfg.bearingMeasurementStd=.1; cfg.bearingObservationVariance=.01;
cfg.nonGaussianOutlierProbability=0; cfg.nonGaussianOutlierScale=8;
cfg.noisyBearingMode='raw'; cfg.noiseDistribution='gaussian';
cfg.immLowQ=1e-4; cfg.immHighQ=.1; cfg.immTransition=[.98,.02;.10,.90];
old=simulateNoisyBearingKFPilot(cfg); raw=simulateBearingIMMClosedLoopPilot(cfg);
assert(isequal(old.trajectory,raw.trajectory),'Gaussian raw reference changed');
assert(isequal(old.timeSeries.globalOrder,raw.timeSeries.globalOrder));
assert(isequal(old.timeSeries.usedBearingRMSE,raw.timeSeries.usedBearingRMSE));
cfg.nonGaussianOutlierProbability=.05; cfg.noiseDistribution='mixture';
old=simulateNonGaussianBearingKFPilot(cfg); raw=simulateBearingIMMClosedLoopPilot(cfg);
assert(isequal(old.trajectory,raw.trajectory),'Mixture raw reference changed');
cfg.noisyBearingMode='selected-kf'; kf=simulateBearingIMMClosedLoopPilot(cfg);
cfg.noisyBearingMode='imm'; imm=simulateBearingIMMClosedLoopPilot(cfg);
for r={kf,imm}
 a=r{1}; assert(isequal(raw.initialPosition,a.initialPosition));
 assert(isequal(raw.initialHeadings,a.initialHeadings));
 assert(isequal(raw.initialNeuralStateUsed,a.initialNeuralStateUsed));
 assert(isequal(raw.pilot.selectedHistory,a.pilot.selectedHistory));
 assert(isequal(raw.pilot.finalNoiseStreamState,a.pilot.finalNoiseStreamState));
 assert(isequal(raw.pilot.finalMixtureStreamState,a.pilot.finalMixtureStreamState));
 assert(isequal(raw.pilot.evalNoise,a.pilot.evalNoise));
 for i=1:numel(raw.pilot.samples)
  assert(isequal(raw.pilot.samples{i}.standardNoise,a.pilot.samples{i}.standardNoise));
  assert(isequal(raw.pilot.samples{i}.outlierMask,a.pilot.samples{i}.outlierMask));
 end
end
assert(all(isnan(imm.pilot.observedBearingHistory)==isnan(kf.pilot.observedBearingHistory),'all'));
% Truth isolation: augmented truth/labels are ignored by the pure controller.
c=cfg; c.NAgents=3; c.NNeurons=20; c.neighborCount=2;
selected=~eye(3); ids=uint8([2,3;1,3;1,2]); y=nan(3); y(selected)=6.27;
obs=struct('selected',selected,'orderedIds',ids,'bearing',y,'trueBearing',zeros(3), ...
 'outlierMask',true(3),'position',zeros(2,3));
alpha=repmat(linspace(0,2*pi,21),3,1); alpha=alpha(:,1:20);
[~,s,~]=bearingIMMClosedLoopControllerStep(obs,alpha,[],c);
obs.bearing(selected)=.02;
[a,sA,~]=bearingIMMClosedLoopControllerStep(obs,alpha,s,c);
obs.trueBearing=ones(3)*3; obs.outlierMask=false(3); obs.position=ones(2,3)*999;
[b,sB,~]=bearingIMMClosedLoopControllerStep(obs,alpha,s,c);
assert(isequal(a,b) && isequal(sA,sB));
% Identical components yield the selected KF field within roundoff.
c.immLowQ=.1; c.immHighQ=.1; c.noisyBearingMode='imm';
[u,~,~]=bearingIMMClosedLoopControllerStep(obs,alpha,[],c);
c.noisyBearingMode='selected-kf';
[v,~,~]=bearingIMMClosedLoopControllerStep(obs,alpha,[],c);
sameField=max(abs(u-v),[],'all'); assert(sameField<1e-13);
assert(all(imm.pilot.estimator.p11>=0));
audit=struct('status','PASS','gaussianRawTrajectoryExact',true, ...
 'mixtureRawTrajectoryExact',true,'pairedInitializationSelectionNoiseAndLabels',true, ...
 'truthPositionLabelIsolation',true,'sameModelFieldMaxDifference',sameField, ...
 'observationHistoriesSameFreshMask',true);
disp(jsonencode(audit));
end
