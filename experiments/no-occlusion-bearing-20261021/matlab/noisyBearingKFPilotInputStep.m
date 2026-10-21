function [input,adj,m,state] = noisyBearingKFPilotInputStep(position,alpha,state,cfg,step)
% Simulator-side sensor and evaluator. Only sanitized observations reach KF.
n=cfg.NAgents; k=cfg.neighborCount; M=cfg.NNeurons;
if isempty(state)
    state.noiseStream=RandStream('mt19937ar','Seed',cfg.seed+15000);
    state.initialNoiseStreamState=state.noiseStream.State;
    state.cv=[];
    state.selectedHistory=zeros(n,k,cfg.NSteps,'uint8');
    state.samples={};
end
standardNoise=randn(state.noiseStream,n,n);
adj=false(n); observed=nan(n); truth=nan(n);
orderedIds=zeros(n,k,'uint8');
for i=1:n
    candidates=[1:i-1,i+1:n];
    if k>=n-1, pick=1:n-1; else, pick=randperm(n-1,k); end
    ids=candidates(pick);
    orderedIds(i,:)=uint8(ids); adj(i,ids)=true;
    delta=position(:,ids)-position(:,i);
    delta=delta-cfg.arenaSize*round(delta/cfg.arenaSize);
    truth(i,ids)=mod(atan2(delta(2,:),delta(1,:)),2*pi);
    observed(i,ids)=mod(truth(i,ids)+cfg.bearingMeasurementStd*standardNoise(i,ids),2*pi);
end
state.selectedHistory(:,:,step)=orderedIds;
obs=struct('selected',adj,'orderedIds',orderedIds,'bearing',observed);
[input,state.cv,d]=noisyBearingKFControllerStep(obs,alpha,state.cv,cfg);
h0=cfg.totalSocialAttraction/n;
if strcmp(cfg.strengthMode,'fixed-total'), h0=(n-1)*h0/k; end
truthInput=zeros(M,n);
m.visibleCount=repmat(n-1,n,1); m.activeCount=repmat(k,n,1);
m.amplitudeSum=repmat(k*h0,n,1); m.distanceWeight=ones(n,1);
m.angularCoverage=zeros(n,1); m.inputConcentration=zeros(n,1);
for i=1:n
    ids=double(orderedIds(i,:));
    delta=abs(alpha(i,:).'-truth(i,ids)); delta=min(delta,2*pi-delta);
    truthInput(:,i)=h0*sum(exp(-0.5*(delta.^2)/(cfg.receptiveFieldWidth^2)),2);
    angles=sort(observed(i,ids));
    m.angularCoverage(i)=1-max([diff(angles),2*pi-angles(end)+angles(1)])/(2*pi);
    m.inputConcentration(i)=abs(mean(exp(1i*d.usedBearing(i,ids))));
end
m.fieldPeak=max(input,[],1).'; m.inputMass=sum(input,1).'*2*pi/M;
eObs=atan2(sin(observed(adj)-truth(adj)),cos(observed(adj)-truth(adj)));
eUsed=atan2(sin(d.usedBearing(adj)-truth(adj)),cos(d.usedBearing(adj)-truth(adj)));
m.observedBearingRMSE=sqrt(mean(eObs.^2));
m.usedBearingRMSE=sqrt(mean(eUsed.^2));
m.inputTruthRMS=sqrt(mean((input-truthInput).^2,'all'));
m.meanKFVariance=mean(d.variance(adj));
if mod(step-1,cfg.pilotSampleStride)==0 || step==cfg.NSteps
    sample=struct('step',step,'linearIndex',find(adj),'orderedIds',orderedIds, ...
        'standardNoise',standardNoise(adj),'observedBearing',observed(adj), ...
        'trueBearing',truth(adj),'usedBearing',d.usedBearing(adj), ...
        'input',input,'truthInput',truthInput,'alpha',alpha);
    state.samples{end+1}=sample;
end
end
