function [input,cv,d] = noisyBearingKFControllerStep(obs,alpha,cv,cfg)
% Controller contract: no positions/ranges/true bearings are read.
n=size(obs.selected,1); M=size(alpha,2);
assert(all(isnan(obs.bearing(~obs.selected))));
ix=find(~eye(n)); fresh=obs.selected(ix); y=obs.bearing(ix);
[cv,nis,innov]=task3BatchTrack(cv,y,fresh,cfg.dt, ...
    cfg.bearingProcessNoise,cfg.bearingObservationVariance,Inf);
used=obs.bearing;
if strcmp(cfg.noisyBearingMode,'kf-posterior')
    used(ix(fresh))=cv.x1(fresh);
end
input=zeros(M,n); variance=nan(n);
variance(ix(fresh))=cv.p11(fresh);
h0=cfg.totalSocialAttraction/n;
if strcmp(cfg.strengthMode,'fixed-total')
    h0=(n-1)*h0/cfg.neighborCount;
end
for i=1:n
    ids=double(obs.orderedIds(i,:));
    delta=abs(alpha(i,:).'-used(i,ids));
    delta=min(delta,2*pi-delta);
    input(:,i)=h0*sum(exp(-0.5*(delta.^2)/(cfg.receptiveFieldWidth^2)),2);
end
d=struct('usedBearing',used,'variance',variance,'nis',nis,'innovation',innov);
end
