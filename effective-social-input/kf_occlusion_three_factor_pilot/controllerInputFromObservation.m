function [input,state,ctrl]=controllerInputFromObservation(alpha,state,obs,cfg,t)
% 公共控制器边界：不接受position/truth，隐藏数据必须为NaN。
n=cfg.NAgents; m=cfg.NNeurons;
if isempty(state), state=initializeInputState(cfg,t-1); end
hidden=obs.selected & ~obs.fresh;
assert(all(isnan(obs.bearing(hidden))) && all(isnan(obs.weight(hidden))));
state.tracks=bearingTrackBankStep(state.tracks,obs,cfg);
% 替代R仅后台预热，不进入主控制输入，且不依据群体结果选择R。
if t<=6200
    for r=cfg.sensitivityRs
        key=sprintf('R_%d',round(-log10(r))); c=cfg; c.R=r;
        state.sensitivityTracks.(key)=bearingTrackBankStep(state.sensitivityTracks.(key),obs,c);
    end
end
state.age(hidden)=state.age(hidden)+cfg.dt; state.age(obs.fresh)=0;
state.lastBearing(obs.fresh)=obs.bearing(obs.fresh); state.lastWeight(obs.fresh)=obs.weight(obs.fresh);
active=obs.fresh; injected=obs.bearing; edgeWeight=zeros(n); edgeWeight(obs.fresh)=obs.weight(obs.fresh);
switch cfg.mode
    case 'hold'
        valid=hidden & isfinite(state.lastBearing); active(valid)=true;
        injected(valid)=state.lastBearing(valid); edgeWeight(valid)=state.lastWeight(valid);
    case 'KF-mean'
        valid=hidden & state.tracks.initialized & isfinite(state.lastWeight); active(valid)=true;
        injected(valid)=state.tracks.bearing(valid); edgeWeight(valid)=state.lastWeight(valid);
    case 'drop-normalized'
        count=sum(obs.fresh,2); ww=cfg.A0./count;
        w=repmat(ww,1,n); edgeWeight(obs.fresh)=w(obs.fresh);
    case {'fresh','drop'}
    otherwise, error('KF:Mode','Unknown method %s',cfg.mode);
end
active=active & edgeWeight>0; edgeWeight(~active)=0; injected(~active)=nan;
input=zeros(m,n); ids=sort(fixedTargetIds(cfg),2);
% 每个体按实际ID升序累加，保持r3的浮点求和次序。
for j=1:cfg.k
    ix=sub2ind([n n],(1:n).',ids(:,j));
    center=injected(ix); center(~isfinite(center))=0;
    delta=atan2(sin(alpha.'-center.'),cos(alpha.'-center.'));
    input=input+exp(-0.5*(delta/cfg.receptiveFieldWidth).^2).*edgeWeight(ix).';
end
state.driverIndex=t;
ctrl=struct();
if ~cfg.collectDiagnostics, return; end
A=sum(edgeWeight,2); aa=A; aa(aa==0)=1;
b=injected; b(~isfinite(b))=0;
z1=sum(edgeWeight.*exp(1i*b),2)./aa; z2=sum(edgeWeight.*exp(2i*b),2)./aa;
directionValid=A>0 & abs(z1)>1e-8; direction=angle(z1); direction(~directionValid)=nan;
p=edgeWeight./aa; pp=p; pp(pp==0)=1; entropy=-sum(p.*log(pp),2); entropy(A==0)=nan;
maxGap=nan(n,1);
for i=1:n, a=angularSetMetrics(injected(i,active(i,:))); maxGap(i)=a.maxGap; end
ctrl.active=active; ctrl.injectedBearing=injected; ctrl.edgeWeight=edgeWeight;
ctrl.lastBearing=state.lastBearing; ctrl.lastWeight=state.lastWeight; ctrl.age=state.age;
ctrl.k_selected=sum(obs.selected,2); ctrl.k_fresh=sum(obs.fresh,2); ctrl.k_active=sum(active,2);
ctrl.A=A; ctrl.Mass=sum(input,1).'*2*pi/m; ctrl.Peak=max(input,[],1).';
ctrl.z1_real=real(z1); ctrl.z1_imag=imag(z1); ctrl.z1_mod=abs(z1);
ctrl.z2_real=real(z2); ctrl.z2_imag=imag(z2); ctrl.z2_mod=abs(z2);
ctrl.direction=direction; ctrl.directionValid=directionValid; ctrl.direction_defined=directionValid;
ctrl.weightEntropy=entropy; ctrl.weightEntropyNormalized=entropy./max(log(sum(active,2)),1);
ctrl.maxGap=maxGap; ctrl.no_input=A==0;
ctrl.normalization_unavailable=strcmp(cfg.mode,'drop-normalized') & any(hidden,2) & sum(obs.fresh,2)==0;
ctrl.predictionBearing=state.tracks.bearing; ctrl.predictionVariance=state.tracks.P11;
ctrl.predictionNIS=state.tracks.nis; ctrl.predictionAge=state.age;
ctrl.predictionVelocity=state.tracks.velocity; ctrl.predictionP12=state.tracks.P12;
ctrl.predictionP22=state.tracks.P22; ctrl.priorBearing=state.tracks.priorBearing;
ctrl.priorVariance=state.tracks.priorVariance;
end
