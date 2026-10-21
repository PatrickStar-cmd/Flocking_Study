function [input,state,d]=bearingIMMClosedLoopControllerStep(obs,alpha,state,cfg)
% Pure bearing controller: no positions, ranges, truth or noise labels read.
if ~strcmp(cfg.noisyBearingMode,'imm')
 legacy=cfg;
 if strcmp(cfg.noisyBearingMode,'selected-kf'), legacy.noisyBearingMode='kf-posterior';
 else, assert(strcmp(cfg.noisyBearingMode,'raw')); legacy.noisyBearingMode='raw'; end
 [input,state,d]=noisyBearingKFControllerStep(obs,alpha,state,legacy);
 d.immHighProbability=nan(size(obs.bearing));
 return
end
n=size(obs.selected,1); M=size(alpha,2);
assert(all(isnan(obs.bearing(~obs.selected))));
ix=find(~eye(n)); fresh=obs.selected(ix); y=obs.bearing(ix);
[state,diagnostic]=bearingIMMTrackStep(state,y,fresh,cfg.dt, ...
 [cfg.immLowQ,cfg.immHighQ],cfg.bearingObservationVariance,cfg.immTransition);
used=obs.bearing; used(ix(fresh))=state.x1(fresh);
variance=nan(n); variance(ix(fresh))=state.p11(fresh);
probability=nan(n); probability(ix)=state.mu(:,2);
input=zeros(M,n); h0=cfg.totalSocialAttraction/n;
if strcmp(cfg.strengthMode,'fixed-total'), h0=(n-1)*h0/cfg.neighborCount; end
for i=1:n
 ids=double(obs.orderedIds(i,:));
 delta=abs(alpha(i,:).'-used(i,ids)); delta=min(delta,2*pi-delta);
 input(:,i)=h0*sum(exp(-0.5*(delta.^2)/(cfg.receptiveFieldWidth^2)),2);
end
d=struct('usedBearing',used,'variance',variance,'immHighProbability',probability, ...
 'modelDiagnostics',diagnostic);
end
