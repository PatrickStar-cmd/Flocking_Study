function result = replayLocalSocialInput(cfg,schedule)
%REPLAYLOCALSOCIALINPUT Replay one focal neural field with prescribed input.
% This is a deterministic local mechanism tool. It has no positions,
% neighbour selection, or hidden state; schedule supplies bearings/weights.

required={'NNeurons','NSteps','dt','beta','globalInhibition', ...
    'receptiveFieldWidth','connectivityExponent','initialNeuralState'};
for i=1:numel(required), assert(isfield(cfg,required{i}),'Missing cfg.%s.',required{i}); end
assert(isfield(schedule,'theta') && isfield(schedule,'weights'));
theta=schedule.theta; weights=schedule.weights;
assert(size(theta,1)==cfg.NSteps && isequal(size(theta),size(weights)));
assert(all(isfinite(theta(:))) && all(isfinite(weights(:))) && all(weights(:)>=0));
assert(size(theta,2)>=1);

M=cfg.NNeurons; T=cfg.NSteps;
alpha=linspace(0,2*pi,M+1).'; alpha=alpha(1:end-1);
W=buildRingConnectivity(M,cfg.connectivityExponent);
u=cfg.initialNeuralState(:);
assert(numel(u)==M,'initialNeuralState must have NNeurons entries.');
names={'inputMass','fieldPeak','concentration','secondConcentration', ...
    'angularEntropy','maxAngularGap','angularCoverage','fieldConcentration', ...
    'activeCount','totalWeight','z1','z2','decodedDirection','bumpConcentration'};
for i=1:numel(names), result.timeSeries.(names{i})=nan(T,1); end
result.input=zeros(M,T); result.neuralState=zeros(M,T); result.neuralOutput=zeros(M,T);
result.alpha=alpha; result.config=cfg; result.schedule=schedule;
cosAlpha=cos(alpha); sinAlpha=sin(alpha);

for t=1:T
    [field,m]=localSocialInputMetrics(theta(t,:),weights(t,:),alpha,cfg.receptiveFieldWidth);
    result.input(:,t)=field;
    result.timeSeries.inputMass(t)=m.inputMass;
    result.timeSeries.fieldPeak(t)=m.fieldPeak;
    result.timeSeries.concentration(t)=m.concentration;
    result.timeSeries.secondConcentration(t)=m.secondConcentration;
    result.timeSeries.angularEntropy(t)=m.angularEntropy;
    result.timeSeries.maxAngularGap(t)=m.maxAngularGap;
    result.timeSeries.angularCoverage(t)=m.angularCoverage;
    result.timeSeries.fieldConcentration(t)=m.fieldConcentration;
    result.timeSeries.activeCount(t)=m.activeCount;
    result.timeSeries.totalWeight(t)=m.totalWeight;
    result.timeSeries.z1(t)=m.z1;
    result.timeSeries.z2(t)=m.z2;

    prior=tanh(cfg.beta*u);
    u=u+cfg.dt*(-u+(W*prior)/M-cfg.globalInhibition+field);
    out=tanh(cfg.beta*u); out(out<0)=0;
    result.neuralState(:,t)=u; result.neuralOutput(:,t)=out;
    cx=sum(out.*cosAlpha); cy=sum(out.*sinAlpha);
    result.timeSeries.decodedDirection(t)=mod(atan2(cy,cx),2*pi);
    result.timeSeries.bumpConcentration(t)=hypot(cx,cy)/(sum(out)+eps);
end
result.finalNeuralState=u;
end
