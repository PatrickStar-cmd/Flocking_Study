function [input,state,ctrl]=controllerInputFromObservation(alpha,state,obs,cfg,t)
%CONTROLLERINPUTFROMOBSERVATION Public controller boundary has no truth input.
n=cfg.NAgents; m=cfg.NNeurons;
if isempty(state), state=initializeInputState(cfg,t-1); end
input=zeros(m,n); active=false(n,n); injected=nan(n,n); edgeWeight=zeros(n,n);
assert(all(isnan(obs.bearing(obs.selected & ~obs.fresh))));
assert(all(isnan(obs.weight(obs.selected & ~obs.fresh))));
for i=1:n
    ids=find(obs.selected(i,:)); fresh=obs.fresh(i,ids);
    for j=ids
        if obs.fresh(i,j)
            center=obs.bearing(i,j); w=obs.weight(i,j);
            assert(isfinite(center)&&isfinite(w)&&w>=0);
            state.lastBearing(i,j)=center; state.lastWeight(i,j)=w; state.age(i,j)=0;
        else
            state.age(i,j)=state.age(i,j)+cfg.dt;
            if ~strcmp(cfg.mode,'hold') || ~isfinite(state.lastBearing(i,j)), continue; end
            center=state.lastBearing(i,j); w=state.lastWeight(i,j);
        end
        if ~isfinite(w) || w<=0, continue; end
        if strcmp(cfg.mode,'drop-normalized') && any(~fresh)
            nFresh=sum(fresh);
            if nFresh==0, continue; end
            w=cfg.A0/nFresh;
        end
        if ~isfinite(w) || w<=0, continue; end
        delta=atan2(sin(alpha(i,:)-center),cos(alpha(i,:)-center));
        input(:,i)=input(:,i)+w*exp(-0.5*(delta/cfg.receptiveFieldWidth).^2).';
        active(i,j)=true; injected(i,j)=center; edgeWeight(i,j)=w;
    end
end
state.driverIndex=t;
A=sum(edgeWeight,2); z1=zeros(n,1); z2=z1; direction=nan(n,1);
directionValid=false(n,1); entropy=nan(n,1); entropyNorm=entropy; maxGap=entropy;
for i=1:n
    ids=find(active(i,:)); if isempty(ids), continue; end
    w=edgeWeight(i,ids); theta=injected(i,ids);
    z1(i)=sum(w.*exp(1i*theta))/A(i); z2(i)=sum(w.*exp(2i*theta))/A(i);
    directionValid(i)=abs(z1(i))>1e-8;
    if directionValid(i), direction(i)=angle(z1(i)); end
    p=w/A(i); entropy(i)=-sum(p.*log(p));
    entropyNorm(i)=entropy(i)/max(log(numel(p)),1);
    s=angularSetMetrics(theta); maxGap(i)=s.maxGap;
end
ctrl.active=active; ctrl.injectedBearing=injected; ctrl.edgeWeight=edgeWeight;
ctrl.lastBearing=state.lastBearing; ctrl.lastWeight=state.lastWeight; ctrl.age=state.age;
ctrl.k_selected=sum(obs.selected,2); ctrl.k_fresh=sum(obs.fresh,2); ctrl.k_active=sum(active,2);
ctrl.A=A; ctrl.Mass=sum(input,1).'*2*pi/m; ctrl.Peak=max(input,[],1).';
ctrl.z1_real=real(z1); ctrl.z1_imag=imag(z1); ctrl.z1_mod=abs(z1);
ctrl.z2_real=real(z2); ctrl.z2_imag=imag(z2); ctrl.z2_mod=abs(z2);
ctrl.direction=direction; ctrl.directionValid=directionValid;
ctrl.direction_defined=directionValid; ctrl.weightEntropy=entropy;
ctrl.weightEntropyNormalized=entropyNorm; ctrl.maxGap=maxGap; ctrl.no_input=A==0;
ctrl.normalization_unavailable=strcmp(cfg.mode,'drop-normalized') & ...
    any(obs.selected & ~obs.fresh,2) & sum(obs.fresh,2)==0;
ctrl.inputField=input;
end
