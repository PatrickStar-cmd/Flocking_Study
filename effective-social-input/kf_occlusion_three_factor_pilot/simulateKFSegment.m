function [seg,cp]=simulateKFSegment(cfg,start,startStep,endStep)
% 与r3相同的神经场、Euler推进、速度和周期边界。每步记录输入前/状态后。
pos=start.position; head=start.heading; u=start.u; alpha=start.alpha; state=start.state;
rng(start.rng); n=cfg.NAgents; m=cfg.NNeurons; L=cfg.arenaSize;
W=buildRingConnectivity(m,cfg.connectivityExponent); rs=max(startStep,cfg.recordFromStep);
nr=max(0,endStep-rs); D=struct();
for f={'step','timeStart','timeEnd','GO','LO','pairDistance','nearestDistance','driverIndex'}
    D.(f{1})=nan(nr,1);
end
D.phase=strings(nr,1);
agents={'k_selected','k_fresh','k_active','A','Mass','Peak','z1_real','z1_imag','z1_mod',...
    'z2_real','z2_imag','z2_mod','direction','directionValid','direction_defined',...
    'weightEntropy','weightEntropyNormalized','maxGap','no_input','normalization_unavailable'};
for f=[agents {'selfFreshRelativeL2','neuralActivityMass','bumpConcentration','rawDecodedAngle',...
        'decodedAnalysisAngle','decoded_direction_defined'}], D.(f{1})=nan(nr,n); end
edges={'edgeWeight','age','lastBearing','lastWeight','observedBearing','observedWeight',...
    'injectedBearing','trueBearing','trueDistance','predictionBearing','predictionVariance',...
    'predictionNIS','predictionAge','predictionError','predictionVelocity','predictionP12',...
    'predictionP22','priorBearing','priorVariance','priorError'};
for f=edges, D.(f{1})=nan(n,n,nr); end
for f={'selected','fresh','hidden','active'}, D.(f{1})=false(n,n,nr); end
D.inputField=zeros(m,n,nr); D.selfFreshField=zeros(m,n,nr);
frames.position=nan(2,n,nr+1); frames.heading=nan(n,nr+1); frames.u=nan(m,n,nr+1);
frames.alpha=alpha; frames.step=rs+(0:nr); cp=struct();
if rs==startStep
    frames.position(:,:,1)=pos; frames.heading(:,1)=head; frames.u(:,:,1)=u;
end
for t=startStep+1:endStep
    c=cfg; c.collectDiagnostics=t>rs;
    [ext,state,d]=fourBranchInputStep(pos,alpha,state,c,t);
    u=u+cfg.dt*(-u+(W*tanh(cfg.beta*u))/m-cfg.globalInhibition+ext);
    y=tanh(cfg.beta*u); y(y<0)=0;
    cx=sum(y.*cos(alpha).',1); cy=sum(y.*sin(alpha).',1);
    moving=abs(cx)+abs(cy)>=cfg.headingThreshold;
    head(moving)=mod(atan2(cy(moving),cx(moving)),2*pi).';
    vx=cfg.speedScale*cx; vy=cfg.speedScale*cy;
    pos=mod(pos+cfg.dt*[vx;vy],L);
    if ismember(t,cfg.checkpointSteps)
        cp.(['step_' num2str(t)])=makeKFCheckpoint(pos,head,u,alpha,state,cfg,t);
    end
    if t==rs, frames.position(:,:,1)=pos; frames.heading(:,1)=head; frames.u(:,:,1)=u; end
    if t<=rs, continue; end
    ix=t-rs; oracle=buildFreshOracle(d.trueBearing,d.targetIds,alpha,cfg);
    speed=hypot(vx,vy); ux=zeros(1,n); uy=ux; nz=speed>1e-12;
    ux(nz)=vx(nz)./speed(nz); uy(nz)=vy(nz)./speed(nz);
    [go,lo,pair,near]=movementMetrics(pos,ux,uy,L);
    D.step(ix)=t; D.timeStart(ix)=(t-1)*cfg.dt; D.timeEnd(ix)=t*cfg.dt;
    if t<cfg.occlusionStartStep, D.phase(ix)="common_prefix";
    elseif t<cfg.recoveryStartStep, D.phase(ix)="occlusion";
    else, D.phase(ix)="recovery"; end
    D.GO(ix)=go; D.LO(ix)=lo; D.pairDistance(ix)=pair/L; D.nearestDistance(ix)=near/L;
    D.driverIndex(ix)=state.driverIndex;
    den=sum(oracle.^2,1); rel=sqrt(sum((ext-oracle).^2,1)./den); rel(den==0)=nan;
    D.selfFreshRelativeL2(ix,:)=rel;
    for f=agents, D.(f{1})(ix,:)=d.(f{1}).'; end
    nd=decodeNeuralActivity(y,alpha);
    for f=fieldnames(nd).', D.(f{1})(ix,:)=nd.(f{1}).'; end
    for f=edges, D.(f{1})(:,:,ix)=d.(f{1}); end
    for f={'selected','fresh','hidden','active'}, D.(f{1})(:,:,ix)=d.(f{1}); end
    D.inputField(:,:,ix)=ext; D.selfFreshField(:,:,ix)=oracle;
    frames.position(:,:,ix+1)=pos; frames.heading(:,ix+1)=head; frames.u(:,:,ix+1)=u;
end
seg.config=cfg; seg.startStep=rs; seg.endStep=endStep; seg.data=D; seg.state=frames;
seg.checkpoints=cp; seg.finalPosition=pos; seg.finalHeading=head; seg.finalNeuralState=u;
seg.finalAlpha=alpha; seg.finalState=state; seg.finalRng=rng;
seg.finalCheckpoint=makeKFCheckpoint(pos,head,u,alpha,state,cfg,endStep);
end
function [go,lo,pair,near]=movementMetrics(pos,ux,uy,L)
n=size(pos,2); go=hypot(sum(ux),sum(uy))/n;
dx=pos(1,:).'-pos(1,:); dy=pos(2,:).'-pos(2,:);
dx=dx-L*round(dx/L); dy=dy-L*round(dy/L); dist=hypot(dx,dy); vals=zeros(n,1);
for i=1:n
    [~,ord]=sort(dist(i,:),'ascend'); ids=ord(1:min(6,n));
    vals(i)=hypot(sum(ux(ids)),sum(uy(ids)))/numel(ids);
end
lo=mean(vals); dist(1:n+1:end)=nan;
pair=mean(dist,'all','omitnan'); near=mean(min(dist,[],2,'omitnan'),'omitnan');
end
