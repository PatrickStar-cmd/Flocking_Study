function testPredictionModelContract()
% Missing measurements must exactly follow the prior, not a hidden update.
dt=.3; q=1e-4; R=1e-4;
t=bearingTrackStep([],6.27,dt,q,R); t=bearingTrackStep(t,.02,dt,q,R);
F=[1 dt;0 1]; Q=q*[dt^3/3 dt^2/2;dt^2/2 dt];
expected=F*t.P*F'+Q; predicted=bearingTrackStep(t,[],dt,q,R);
assert(max(abs(predicted.P-expected),[],'all')<1e-12);
assert(~predicted.updated && isnan(predicted.nis) && predicted.age==dt);
recovered=bearingTrackStep(predicted,.03,dt,q,R);
assert(recovered.updated && recovered.age==0 && isfinite(recovered.nis));
cfg=defaultEffectiveInputConfig(); cfg.NAgents=4; cfg.NNeurons=20;
cfg.NSteps=10; cfg.neighborCount=2; cfg.occlusionStartStep=3;
cfg.occlusionDurationSteps=2; cfg.occlusionHiddenCount=1;
p=[0 10 20 30;0 1 2 3]; alpha=repmat((0:19)*2*pi/20,4,1);
rng(123); before=rng;
[~,state,~,~,a]=scheduledOcclusionInput(p,alpha,[],cfg,1);
after=rng; assert(isequal(before,after));
cfg.predictionMode='kalman'; [~,~,~,~,b]=scheduledOcclusionInput(p,alpha,[],cfg,1);
assert(isequaln(a.observedBearing,b.observedBearing));
[~,state,~,~,hidden]=scheduledOcclusionInput(p,alpha,state,cfg,3);
assert(isnan(hidden.observedBearing(1,2)));
[~,~,~,~,visible]=scheduledOcclusionInput(p,alpha,state,cfg,5);
assert(isfinite(visible.observedBearing(1,2)));
disp('PASS: missing-update contract, reacquisition age, shared observation noise, RNG isolation');
testOcclusionPrediction;
end
