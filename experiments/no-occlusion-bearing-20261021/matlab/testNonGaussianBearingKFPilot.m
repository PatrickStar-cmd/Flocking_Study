function audit=testNonGaussianBearingKFPilot()
baseAudit=testNoisyBearingKFPilot();
cfg=defaultEffectiveInputConfig();cfg.NSteps=160;cfg.burnIn=80;cfg.seed=5199;
cfg.neighborCount=15;cfg.selectionPolicy='random';cfg.recordStride=1;
cfg.storeTimeSeries=true;cfg.storeTrajectory=true;cfg.pilotSampleStride=20;
cfg.bearingMeasurementStd=.1;cfg.bearingObservationVariance=.1^2;
cfg.nonGaussianOutlierProbability=0;cfg.nonGaussianOutlierScale=8;
for mode={'raw','kf-posterior'}
    cfg.noisyBearingMode=mode{1};
    gaussian=simulateNoisyBearingKFPilot(cfg);
    reduced=simulateNonGaussianBearingKFPilot(cfg);
    assert(isequal(gaussian.trajectory,reduced.trajectory),'p=0 differs from Gaussian dynamics');
    assert(isequal(gaussian.timeSeries.globalOrder,reduced.timeSeries.globalOrder));
    assert(isequal(gaussian.timeSeries.usedBearingRMSE,reduced.timeSeries.usedBearingRMSE));
end
cfg.nonGaussianOutlierProbability=.05;cfg.noisyBearingMode='raw';
raw=simulateNonGaussianBearingKFPilot(cfg);cfg.noisyBearingMode='kf-posterior';
kf=simulateNonGaussianBearingKFPilot(cfg);
assert(isequal(raw.initialPosition,kf.initialPosition));
assert(isequal(raw.initialHeadings,kf.initialHeadings));
assert(isequal(raw.initialNeuralStateUsed,kf.initialNeuralStateUsed));
assert(isequal(raw.pilot.selectedHistory,kf.pilot.selectedHistory));
assert(isequal(raw.pilot.finalNoiseStreamState,kf.pilot.finalNoiseStreamState));
assert(isequal(raw.pilot.finalMixtureStreamState,kf.pilot.finalMixtureStreamState));
assert(isequal(raw.pilot.evalNoise,kf.pilot.evalNoise));
for i=1:numel(raw.pilot.samples)
    assert(isequal(raw.pilot.samples{i}.standardNoise,kf.pilot.samples{i}.standardNoise));
    assert(isequal(raw.pilot.samples{i}.outlierMask,kf.pilot.samples{i}.outlierMask));
end
[e,labels,~,~]=drawNonGaussianBearingNoise(struct(),cfg,[1000,1000]);
v=var(e(:),1);mu=mean(e(:));fraction=mean(labels(:));
kurtosis=mean(e(:).^4)/mean(e(:).^2)^2;
assert(abs(mu)<.01 && abs(v-1)<.03 && abs(fraction-.05)<.002 && kurtosis>10);
audit=struct('status','PASS','baseGaussianContracts',baseAudit, ...
    'pZeroExactGaussianReduction',true,'mixturePairing',true, ...
    'empiricalUnitMean',mu,'empiricalUnitVariance',v, ...
    'empiricalOutlierFraction',fraction,'empiricalKurtosis',kurtosis, ...
    'theoreticalKurtosis',3*(.95+.05*8^4)/(.95+.05*8^2)^2);
disp(jsonencode(audit));
end
