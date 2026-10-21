function audit = testNoisyBearingKFPilot()
cfg=defaultEffectiveInputConfig(); cfg.NSteps=160; cfg.burnIn=80;
cfg.seed=5199; cfg.neighborCount=15; cfg.selectionPolicy='random';
cfg.recordStride=1; cfg.storeTimeSeries=true; cfg.storeTrajectory=true;
cfg.noisyBearingMode='raw'; cfg.pilotSampleStride=20;
cfg.bearingMeasurementStd=0;
old=simulateEffectiveSocialInput(cfg); raw=simulateNoisyBearingKFPilot(cfg);
assert(isequal(old.trajectory,raw.trajectory),'Zero-noise dynamics differ from legacy baseline');
assert(isequal(old.timeSeries.globalOrder,raw.timeSeries.globalOrder));
assert(isequal(old.timeSeries.meanPairDistanceNorm,raw.timeSeries.meanPairDistanceNorm));
cfg.bearingMeasurementStd=.1; cfg.bearingObservationVariance=.1^2;
raw=simulateNoisyBearingKFPilot(cfg); cfg.noisyBearingMode='kf-posterior';
kf=simulateNoisyBearingKFPilot(cfg);
assert(isequal(raw.initialPosition,kf.initialPosition));
assert(isequal(raw.initialHeadings,kf.initialHeadings));
assert(isequal(raw.initialNeuralStateUsed,kf.initialNeuralStateUsed));
assert(isequal(raw.pilot.selectedHistory,kf.pilot.selectedHistory));
assert(isequal(raw.pilot.finalNoiseStreamState,kf.pilot.finalNoiseStreamState));
for i=1:numel(raw.pilot.samples)
    assert(isequal(raw.pilot.samples{i}.standardNoise,kf.pilot.samples{i}.standardNoise));
end
assert(max(abs(raw.timeSeries.observedBearingRMSE-raw.timeSeries.usedBearingRMSE))<1e-14);
% Pure-controller truth-isolation, circle wrapping and PSD contracts.
n=3; c=cfg; c.NAgents=n; c.NNeurons=20; c.neighborCount=2;
alpha=repmat(linspace(0,2*pi,21),n,1); alpha=alpha(:,1:20);
selected=~eye(n); ids=uint8([2,3;1,3;1,2]);
y=nan(n); y(selected)=2*pi-.01;
obs=struct('selected',selected,'orderedIds',ids,'bearing',y,'trueBearing',zeros(n));
[~,state,~]=noisyBearingKFControllerStep(obs,alpha,[],c);
obs.bearing(selected)=.01; obs.trueBearing=ones(n)*5;
[a,stateA,dA]=noisyBearingKFControllerStep(obs,alpha,state,c);
obs.trueBearing=ones(n)*.1;
[b,stateB,~]=noisyBearingKFControllerStep(obs,alpha,state,c);
assert(isequal(a,b) && isequal(stateA,stateB));
assert(max(abs(dA.innovation))<.1);
assert(all(stateA.p11>=0 & stateA.p22>=0));
assert(all(stateA.p11.*stateA.p22-stateA.p12.^2>=-1e-14));
audit=struct('status','PASS','zeroNoiseLegacyExact',true,'truthIsolation',true, ...
    'pairedInitialization',true,'pairedSelection',true,'pairedNoise',true, ...
    'wrappedInnovation',true,'covariancePSD',true);
disp(jsonencode(audit));
end
