function [input,adj,m,state]=oracleBearingDiagnosticInputStep(position,alpha,state,cfg,step)
% Oracle diagnostic only: exact selected bearings drive the field, no filter.
% Unused noise streams are consumed to preserve pairing with archived runs.
n=cfg.NAgents; k=cfg.neighborCount; M=cfg.NNeurons; ix=find(~eye(n));
if isempty(state)
 state.estimator=[]; state.selectedHistory=zeros(n,k,cfg.NSteps,'uint8');
 state.samples={}; state.linearEdges=ix;
 state.observedBearingHistory=nan(cfg.NSteps,numel(ix));
 state.usedBearingEvaluation=nan(cfg.NSteps-cfg.burnIn,numel(ix));
 state.highProbabilityEvaluation=nan(cfg.NSteps-cfg.burnIn,numel(ix));
 state.evalNoise=struct('count',0,'sum',0,'sum2',0,'sum4',0,'outlierCount',0,'maxAbs',0);
end
[standardNoise,outlierMask,standardNormal,state]=drawNonGaussianBearingNoise(state,cfg,[n,n]);
adj=false(n); observed=nan(n); truth=nan(n); orderedIds=zeros(n,k,'uint8');
for i=1:n
 candidates=[1:i-1,i+1:n];
 if k>=n-1, pick=1:n-1; else, pick=randperm(n-1,k); end
 ids=candidates(pick); orderedIds(i,:)=uint8(ids); adj(i,ids)=true;
 delta=position(:,ids)-position(:,i);
 delta=delta-cfg.arenaSize*round(delta/cfg.arenaSize);
 truth(i,ids)=mod(atan2(delta(2,:),delta(1,:)),2*pi);
 observed(i,ids)=mod(truth(i,ids)+cfg.bearingMeasurementStd*standardNoise(i,ids),2*pi);
end
state.selectedHistory(:,:,step)=orderedIds;
state.observedBearingHistory(step,:)=observed(ix);
obs=struct('selected',adj,'orderedIds',orderedIds,'bearing',observed);
assert(cfg.bearingMeasurementStd==0 && strcmp(cfg.noisyBearingMode,'raw'));
[input,d]=oracleOnlyField(obs,alpha,cfg); state.estimator=[];
if step>cfg.burnIn
 at=step-cfg.burnIn; state.usedBearingEvaluation(at,:)=d.usedBearing(ix);
 state.highProbabilityEvaluation(at,:)=d.immHighProbability(ix);
 e=standardNoise(adj); stats=state.evalNoise;
 stats.count=stats.count+numel(e); stats.sum=stats.sum+sum(e);
 stats.sum2=stats.sum2+sum(e.^2); stats.sum4=stats.sum4+sum(e.^4);
 stats.outlierCount=stats.outlierCount+sum(outlierMask(adj));
 stats.maxAbs=max(stats.maxAbs,max(abs(e))); state.evalNoise=stats;
end
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
m.meanIMMHighProbability=mean(d.immHighProbability(adj));
if mod(step-1,cfg.pilotSampleStride)==0 || step==cfg.NSteps
 sample=struct('step',step,'linearIndex',find(adj),'orderedIds',orderedIds, ...
  'standardNoise',standardNoise(adj),'outlierMask',outlierMask(adj), ...
  'standardNormal',standardNormal(adj),'observedBearing',observed(adj), ...
  'trueBearing',truth(adj),'usedBearing',d.usedBearing(adj), ...
  'highProbability',d.immHighProbability(adj), ...
  'input',input,'truthInput',truthInput,'alpha',alpha);
 state.samples{end+1}=sample;
end
end

function [input,d]=oracleOnlyField(obs,alpha,cfg)
% Diagnostic arm only. No KF/IMM is called or updated with exact bearings.
n=size(obs.selected,1); input=zeros(cfg.NNeurons,n);
h0=cfg.totalSocialAttraction/n;
if strcmp(cfg.strengthMode,'fixed-total'), h0=(n-1)*h0/cfg.neighborCount; end
for i=1:n
 ids=double(obs.orderedIds(i,:));
 delta=abs(alpha(i,:).'-obs.bearing(i,ids)); delta=min(delta,2*pi-delta);
 input(:,i)=h0*sum(exp(-0.5*(delta.^2)/(cfg.receptiveFieldWidth^2)),2);
end
d=struct('usedBearing',obs.bearing,'variance',nan(n),'immHighProbability',nan(n));
end
