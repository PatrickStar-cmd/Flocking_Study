function testDistanceWeight()
cfg=defaultEffectiveInputConfig(); d=[0 100 500];
[w,n]=distanceWeight(d,cfg); assert(isequal(w,ones(size(d)))); assert(abs(sum(n)-1)<1e-12);
cfg.distanceWeightMode='exponential'; cfg.distanceWeightScaleFraction=.5;
[w,n]=distanceWeight(d,cfg); assert(w(1)>w(2) && w(2)>w(3)); assert(abs(sum(n)-1)<1e-12);
disp('PASS: equal and exponential distance weights');
cfg=defaultEffectiveInputConfig(); cfg.NSteps=80; cfg.burnIn=40;
cfg.selectionPolicy='random'; cfg.neighborCount=7;
for mode={'original','fixed-total'}
 cfg.strengthMode=mode{1};
 old=simulateEffectiveSocialInputBeforeDistance(cfg); new=simulateEffectiveSocialInput(cfg);
 assert(isequal(old.finalNeuralState,new.finalNeuralState));
 names=fieldnames(old.timeSeries);
 for j=1:numel(names), assert(isequal(old.timeSeries.(names{j}),new.timeSeries.(names{j}))); end
end
cfg.distanceWeightMode='exponential'; cfg.strengthMode='fixed-total';
r=simulateEffectiveSocialInput(cfg);
assert(max(abs(r.timeSeries.meanInputAmplitude-.228))<1e-12);
cfg.neighborCount=0; r=simulateEffectiveSocialInput(cfg);
assert(all(r.timeSeries.meanInputMass==0));
% First-step positions: torus distance from x=1 to x=999 is 2, not 998.
cfg.NAgents=3; cfg.neighborCount=2; cfg.selectionPolicy='all'; cfg.strengthMode='original';
cfg.NSteps=1; cfg.burnIn=0; cfg.initialXY=[1 999 11;0 0 0];
r=simulateEffectiveSocialInput(cfg);
expected=.24/3*mean([exp(-2/500)+exp(-10/500),exp(-2/500)+exp(-12/500),exp(-10/500)+exp(-12/500)]);
assert(abs(r.timeSeries.meanInputAmplitude-expected)<1e-12);
disp('PASS: exact equal-mode regression, fixed mass, zero neighbors, periodic distance');
end
