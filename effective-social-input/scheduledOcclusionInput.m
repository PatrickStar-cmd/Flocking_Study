function [input,state,adj,m,record] = scheduledOcclusionInput(position,alpha,state,cfg,step)
% Simulator-side sensor. Hidden truth is not passed to the controller.
if isempty(state)
    stream=RandStream('mt19937ar','Seed',cfg.seed);
    sensorNoise=cfg.bearingMeasurementStd*randn(stream,cfg.NAgents,cfg.NAgents,cfg.NSteps);
else
    sensorNoise=state.sensorNoise;
end
n=cfg.NAgents; selected=false(n); fresh=false(n); bearings=nan(n); weights=nan(n); trueBearings=nan(n);
hidden=step>=cfg.occlusionStartStep && step<cfg.occlusionStartStep+cfg.occlusionDurationSteps;
for i=1:n
 ids=mod(i-1+(1:cfg.neighborCount),n)+1; selected(i,ids)=true;
 count=min(cfg.occlusionHiddenCount,numel(ids));
 visible=ids;
 if hidden
  if isempty(cfg.occlusionHiddenIds)
   visible=ids(count+1:end);
  else
   hiddenIds=cfg.occlusionHiddenIds(i,1:count);
   visible=ids(~ismember(ids,hiddenIds));
  end
 end
 fresh(i,visible)=true;
 dx=position(:,visible)-position(:,i); dx=dx-cfg.arenaSize*round(dx/cfg.arenaSize);
 bearings(i,visible)=mod(atan2(dx(2,:),dx(1,:))+sensorNoise(i,visible,step),2*pi);
 weights(i,visible)=distanceWeight(hypot(dx(1,:),dx(2,:)),cfg);
 dxAll=position(:,ids)-position(:,i); dxAll=dxAll-cfg.arenaSize*round(dxAll/cfg.arenaSize);
 trueBearings(i,ids)=mod(atan2(dxAll(2,:),dxAll(1,:)),2*pi);
end
obs=struct('selected',selected,'fresh',fresh,'bearing',bearings,'weight',weights);
opt=struct('dt',cfg.dt,'q',cfg.bearingProcessNoise, ...
 'gmQ',cfg.bearingGaussMarkovProcessNoise,'velocityTau',cfg.bearingVelocityTau, ...
 'R',cfg.bearingObservationVariance, ...
 'mode',cfg.predictionMode,'maxAge',cfg.predictionMaxAge,'maxVariance',cfg.predictionMaxVariance, ...
 'sigma',cfg.receptiveFieldWidth,'h0',cfg.totalSocialAttraction/n);
if ~isempty(state), state=rmfield(state,'sensorNoise'); end
[input,state,d]=occlusionInputStep(obs,alpha,state,opt); adj=d.active;
state.sensorNoise=sensorNoise;
m.visibleCount=sum(fresh,2); m.activeCount=sum(adj,2);
m.amplitudeSum=sum(adj.*fillmissing(state.weight,'constant',0),2)*opt.h0;
m.inputMass=sum(input,1)'*2*pi/cfg.NNeurons;
m.fieldPeak=max(input,[],1)';
m.distanceWeight=sum(adj.*fillmissing(state.weight,'constant',0),2)./max(1,m.activeCount);
% Diagnostics of the actual injected field, not hidden true bearings.
z=sum(input.*exp(1i*alpha'),1); m.inputConcentration=abs(z)'./max(eps,sum(input,1)');
m.angularCoverage=nan(n,1); % Deliberately unavailable: no true-state shortcut.
record=struct('fresh',fresh,'active',logical(adj),'variance',d.variance,'age',d.age,'weight',d.weight);
record.observedBearing=bearings;
record.injectedBearing=d.center;
record.trueBearing=trueBearings;
end
