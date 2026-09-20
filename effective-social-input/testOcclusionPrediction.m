function testOcclusionPrediction()
testDistanceWeight;
t=bearingTrackStep([],2*pi-.01,.3,1e-4,1e-4);
t=bearingTrackStep(t,.01,.3,1e-4,1e-4);
assert(abs(t.x(2))<.2);
p=t.P(1,1); t=bearingTrackStep(t,[],.3,1e-4,1e-4);
assert(t.P(1,1)>p && min(eig(t.P))>=-1e-12);
cfg=defaultEffectiveInputConfig(); cfg.NAgents=4; cfg.NNeurons=40;
cfg.neighborCount=2; cfg.selectionPolicy='fixed-id'; cfg.scheduledOcclusion=true;
cfg.occlusionStartStep=3; cfg.occlusionDurationSteps=3; cfg.occlusionHiddenCount=1;
cfg.distanceWeightMode='exponential'; cfg.predictionMode='kalman';
cfg.NSteps=12; cfg.burnIn=1;
pos=[0 10 20 30;0 1 2 3]; alpha=repmat((0:39)*2*pi/40,4,1);
[~,state]=scheduledOcclusionInput(pos,alpha,[],cfg,1);
[~,state]=scheduledOcclusionInput(pos,alpha,state,cfg,2);
[a,~,~,~,rec]=scheduledOcclusionInput(pos,alpha,state,cfg,3);
pos(:,2)=[800;800]; [b]=scheduledOcclusionInput(pos,alpha,state,cfg,3);
assert(isequal(a(:,1),b(:,1))); % Target 2 is hidden to focal 1.
assert(~rec.fresh(1,2) && rec.active(1,2));
for mode={'none','hold','kalman','kalman-gm'}
 cfg.predictionMode=mode{1}; r=simulateEffectiveSocialInput(cfg);
 assert(all(r.timeSeries.meanVisibleNeighbors(3:5)==1));
 assert(r.timeSeries.meanVisibleNeighbors(6)==2);
 if strcmp(mode{1},'none'), assert(all(r.timeSeries.meanActiveNeighbors(3:5)==1));
 else, assert(all(r.timeSeries.meanActiveNeighbors(3:5)==2)); end
 assert(all(isfinite(r.timeSeries.globalOrder)));
end
% No-event input is identical for all memory modes.
cfg.occlusionDurationSteps=0; cfg.predictionMode='none'; a=simulateEffectiveSocialInput(cfg);
cfg.predictionMode='kalman'; b=simulateEffectiveSocialInput(cfg);
assert(isequal(a.finalNeuralState,b.finalNeuralState));
obs=struct('selected',[false true;false false],'fresh',[false true;false false], ...
 'bearing',[NaN 0;NaN NaN],'weight',[NaN .5;NaN NaN]);
opt=struct('dt',.3,'q',1e-4,'gmQ',1e-5,'velocityTau',10,'R',1e-4,'mode','kalman', ...
 'maxAge',100,'maxVariance',10,'sigma',.4,'h0',.012);
angles=repmat((0:99)*2*pi/100,2,1);
[sharp,state]=occlusionInputStep(obs,angles,[],opt);
obs.fresh(:)=false; obs.bearing(:)=NaN; obs.weight(:)=NaN;
state.tracks{1,2}.P=diag([.4 .01]);
[wide,~,~]=occlusionInputStep(obs,angles,state,opt);
assert(abs(sum(sharp(:,1))-sum(wide(:,1)))<1e-10 && max(wide(:,1))<max(sharp(:,1)));
opt.maxVariance=.01; [removed]=occlusionInputStep(obs,angles,state,opt);
assert(all(removed(:)==0));
disp('PASS: wrapped innovation, PSD, hidden-state isolation, schedule, no-event equality');
disp('PASS: circular mass conservation, peak reduction, covariance cutoff');
end
