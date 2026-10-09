function [segment,cp]=simulateProtocolSegment(cfg,startState,startStep,endStep)
%SIMULATEPROTOCOLSEGMENT Continue one deterministic branch from a checkpoint.
state=startState.state; pos=startState.position; head=startState.heading; u=startState.u; alpha=startState.alpha; rng(startState.rng);
n=cfg.NAgents; m=cfg.NNeurons; L=cfg.arenaSize; W=buildRingConnectivity(m,cfg.connectivityExponent); nSteps=endStep-startStep;
D=struct(); scalarNames={'step','timeStart','timeEnd','GO','LO','pairDistance','nearestDistance','selfFreshA'};
for q=1:numel(scalarNames), D.(scalarNames{q})=nan(nSteps,1); end
D.phase=strings(nSteps,1); D.selfFreshField=zeros(m,n,nSteps); D.selfFreshMass=zeros(nSteps,n); D.selfFreshPeak=zeros(nSteps,n); D.selfFreshError=zeros(nSteps,n);
perAgent={'k_selected','k_fresh','k_active','A','Mass','Peak','z1_real','z1_imag','z1_mod','z2_real','z2_imag','z2_mod','direction','directionValid','weightEntropy','maxGap','no_input','normalization_unavailable'};
for q=1:numel(perAgent), D.(perAgent{q})=nan(nSteps,n); end
logicalFields={'directionValid','no_input','normalization_unavailable'}; for q=1:numel(logicalFields), D.(logicalFields{q})=false(nSteps,n); end
D.inputField=zeros(m,n,nSteps); D.edgeWeight=zeros(n,n,nSteps); D.age=nan(n,n,nSteps);
D.selected=false(n,n,nSteps); D.fresh=false(n,n,nSteps); D.hidden=false(n,n,nSteps); D.active=false(n,n,nSteps);
D.observedBearing=nan(n,n,nSteps); D.injectedBearing=nan(n,n,nSteps); D.trueBearing=nan(n,n,nSteps);
stateFrames=struct('position',nan(2,n,nSteps+1),'heading',nan(n,nSteps+1),'step',startStep+(0:nSteps));
stateFrames.position(:,:,1)=pos; stateFrames.heading(:,1)=head; checkpoints=struct();
for local=1:nSteps
    t=startStep+local;
    [ext,state,diag]=fourBranchInputStep(pos,alpha,state,cfg,t);
    oracle=buildFreshOracle(diag.trueBearing,diag.targetIds,alpha,cfg);
    neural=tanh(cfg.beta*u); u=u+cfg.dt*(-u+(W*neural)/m+ext);
    y=tanh(cfg.beta*u); y(y<0)=0; cX=sum(y.*cos(alpha).',1); cY=sum(y.*sin(alpha).',1);
    moving=abs(cX)+abs(cY)>1e-12; head(moving)=mod(atan2(cY(moving),cX(moving)),2*pi).';
    velX=cfg.speedScale*cX; velY=cfg.speedScale*cY; pos=mod(pos+cfg.dt*[velX;velY],L);
    speed=hypot(velX,velY); ux=zeros(1,n); uy=zeros(1,n); nz=speed>1e-12; ux(nz)=velX(nz)./speed(nz); uy(nz)=velY(nz)./speed(nz);
    [go,lo]=orderMetrics(pos,ux,uy,L,5); [pair,near]=distanceMetrics(pos,L);
    D.step(local)=t; D.timeStart(local)=(t-1)*cfg.dt; D.timeEnd(local)=t*cfg.dt;
    if t<=6200, D.phase(local)="common_prefix"; elseif t<6235, D.phase(local)="occlusion"; else, D.phase(local)="recovery"; end
    D.GO(local)=go; D.LO(local)=lo; D.pairDistance(local)=pair/L; D.nearestDistance(local)=near/L; D.selfFreshA(local)=cfg.A0;
    D.selfFreshField(:,:,local)=oracle; D.selfFreshMass(local,:)=sum(oracle,1)*2*pi/m; D.selfFreshPeak(local,:)=max(oracle,[],1); D.selfFreshError(local,:)=max(abs(diag.inputField-oracle),[],1);
    for q=1:numel(perAgent), name=perAgent{q}; D.(name)(local,:)=diag.(name).'; end
    D.inputField(:,:,local)=diag.inputField; D.edgeWeight(:,:,local)=diag.edgeWeight; D.age(:,:,local)=diag.age;
    D.selected(:,:,local)=diag.selected; D.fresh(:,:,local)=diag.fresh; D.hidden(:,:,local)=diag.hidden; D.active(:,:,local)=diag.active;
    D.observedBearing(:,:,local)=diag.observedBearing; D.injectedBearing(:,:,local)=diag.injectedBearing; D.trueBearing(:,:,local)=diag.trueBearing;
    stateFrames.position(:,:,local+1)=pos; stateFrames.heading(:,local+1)=head;
    if t==6000 || t==6200, checkpoints.(['step_' num2str(t)])=makeCheckpoint(t,pos,head,u,alpha,state,cfg); end
end
segment=struct('config',cfg,'startStep',startStep,'endStep',endStep,'data',D,'state',stateFrames,'checkpoints',checkpoints,'finalPosition',pos,'finalHeading',head,'finalNeuralState',u,'finalAlpha',alpha,'finalState',state); cp=segment.checkpoints;
end
function oracle=buildFreshOracle(trueBearing,targetIds,alpha,cfg)
n=cfg.NAgents; m=cfg.NNeurons; oracle=zeros(m,n);
for i=1:n
 ids=targetIds(i,:); theta=trueBearing(i,ids);
 for q=1:numel(ids)
  delta=atan2(sin(alpha(i,:)-theta(q)),cos(alpha(i,:)-theta(q))); oracle(:,i)=oracle(:,i)+cfg.w0*exp(-0.5*(delta/cfg.receptiveFieldWidth).^2).';
 end
end
end
function cp=makeCheckpoint(t,pos,head,u,alpha,state,cfg)
cp=struct('step',t,'position',pos,'heading',head,'u',u,'alpha',alpha,'targetIds',targetIds(cfg),'mask',state.mask,'maskReady',state.maskReady,'maskDiagnostics',state.maskDiagnostics,'lastBearing',state.lastBearing,'lastWeight',state.lastWeight,'age',state.age,'rng',rng,'driverIndex',state.driverIndex,'config',cfg);
end
function ids=targetIds(cfg), ids=zeros(cfg.NAgents,cfg.k); for i=1:cfg.NAgents, ids(i,:)=mod((i-1)+(1:cfg.k),cfg.NAgents)+1; end, end
function [go,lo]=orderMetrics(pos,ux,uy,L,k)
n=size(pos,2); go=hypot(sum(ux),sum(uy))/n; vals=zeros(n,1); dx=pos(1,:).'-pos(1,:); dy=pos(2,:).'-pos(2,:); dx=dx-L*round(dx/L); dy=dy-L*round(dy/L); dist=hypot(dx,dy); dist(1:n+1:end)=0;
for i=1:n, [~,ord]=sort(dist(i,:),'ascend'); ids=ord(1:min(k+1,n)); vals(i)=hypot(sum(ux(ids)),sum(uy(ids)))/numel(ids); end; lo=mean(vals); end
function [pair,near]=distanceMetrics(pos,L)
n=size(pos,2); dx=pos(1,:).'-pos(1,:); dy=pos(2,:).'-pos(2,:); dx=dx-L*round(dx/L); dy=dy-L*round(dy/L); d=hypot(dx,dy); d(1:n+1:end)=nan; pair=mean(d,'all','omitnan'); near=mean(min(d,[],2,'omitnan'),'omitnan'); end

