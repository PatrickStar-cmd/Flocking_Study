function [segment,cp]=simulateProtocolSegment(cfg,startState,startStep,endStep)
%SIMULATEPROTOCOLSEGMENT Same core path for uninterrupted and resumed runs.
state=startState.state; pos=startState.position; head=startState.heading;
u=startState.u; alpha=startState.alpha; rng(startState.rng);
n=cfg.NAgents; m=cfg.NNeurons; L=cfg.arenaSize;
W=buildRingConnectivity(m,cfg.connectivityExponent);
recordStart=max(startStep,cfg.recordFromStep); nr=max(0,endStep-recordStart);
D=struct(); scalar={'step','localStep','timeStart','timeEnd','GO','LO','pairDistance','nearestDistance','driverIndex'};
for q=1:numel(scalar), D.(scalar{q})=nan(nr,1); end
D.phase=strings(nr,1);
agent={'k_selected','k_fresh','k_active','A','Mass','Peak','z1_real','z1_imag','z1_mod',...
    'z2_real','z2_imag','z2_mod','direction','directionValid','direction_defined',...
    'weightEntropy','weightEntropyNormalized','maxGap','no_input','normalization_unavailable'};
extra={'selfFreshA','selfFreshMass','selfFreshPeak','selfFreshError','selfFreshRelativeL2',...
    'neuralActivityMass','bumpConcentration','rawDecodedAngle','decodedAnalysisAngle','decoded_direction_defined',...
    'hiddenR1','hiddenR2','hiddenMaxGap','remainingR1','remainingR2','remainingMaxGap','allR1','allR2','allMaxGap'};
for name=[agent extra], D.(name{1})=nan(nr,n); end
edge={'edgeWeight','age','lastBearing','lastWeight','observedBearing','observedWeight','injectedBearing','trueBearing','trueDistance'};
for name=edge, D.(name{1})=nan(n,n,nr); end
for name={'selected','fresh','hidden','active'}, D.(name{1})=false(n,n,nr); end
D.inputField=zeros(m,n,nr); D.selfFreshField=zeros(m,n,nr);
frames.position=nan(2,n,nr+1); frames.heading=nan(n,nr+1);
frames.u=nan(m,n,nr+1); frames.alpha=alpha; frames.step=recordStart+(0:nr);
checkpoints=struct();
if recordStart==startStep
    frames.position(:,:,1)=pos; frames.heading(:,1)=head; frames.u(:,:,1)=u;
end
for t=startStep+1:endStep
    [ext,state,d]=fourBranchInputStep(pos,alpha,state,cfg,t);
    neural=tanh(cfg.beta*u);
    u=u+cfg.dt*(-u+(W*neural)/m-cfg.globalInhibition+ext);
    y=tanh(cfg.beta*u); y(y<0)=0;
    cx=sum(y.*cos(alpha).',1); cy=sum(y.*sin(alpha).',1);
    moving=abs(cx)+abs(cy)>=cfg.headingThreshold;
    head(moving)=mod(atan2(cy(moving),cx(moving)),2*pi).';
    vx=cfg.speedScale*cx; vy=cfg.speedScale*cy;
    pos=mod(pos+cfg.dt*[vx;vy],L);
    if ismember(t,cfg.checkpointSteps)
        c.step=t; c.position=pos; c.heading=head; c.u=u; c.alpha=alpha;
        c.targetIds=fixedTargetIds(cfg); c.mask=state.mask; c.maskReady=state.maskReady;
        c.maskDiagnostics=state.maskDiagnostics; c.lastBearing=state.lastBearing;
        c.lastWeight=state.lastWeight; c.age=state.age; c.rng=rng;
        c.driverIndex=state.driverIndex; c.config=cfg;
        checkpoints.(['step_' num2str(t)])=c;
    end
    if t==recordStart
        frames.position(:,:,1)=pos; frames.heading(:,1)=head; frames.u(:,:,1)=u;
    end
    if t<=recordStart, continue; end
    local=t-recordStart;
    oracle=buildFreshOracle(d.trueBearing,d.targetIds,alpha,cfg);
    speed=hypot(vx,vy); ux=zeros(1,n); uy=ux; nz=speed>1e-12;
    ux(nz)=vx(nz)./speed(nz); uy(nz)=vy(nz)./speed(nz);
    [go,lo]=orderMetrics(pos,ux,uy,L,5); [pair,near]=distanceMetrics(pos,L);
    D.step(local)=t; D.localStep(local)=t-6000;
    D.timeStart(local)=(t-1)*cfg.dt; D.timeEnd(local)=t*cfg.dt;
    if t<cfg.occlusionStartStep
        D.phase(local)="common_prefix";
    elseif t<cfg.occlusionStartStep+cfg.occlusionDurationSteps
        D.phase(local)="occlusion";
    else
        D.phase(local)="recovery";
    end
    D.GO(local)=go; D.LO(local)=lo; D.pairDistance(local)=pair/L;
    D.nearestDistance(local)=near/L; D.driverIndex(local)=state.driverIndex;
    D.selfFreshField(:,:,local)=oracle; D.selfFreshA(local,:)=cfg.k*cfg.w0;
    D.selfFreshMass(local,:)=sum(oracle,1)*2*pi/m; D.selfFreshPeak(local,:)=max(oracle,[],1);
    D.selfFreshError(local,:)=max(abs(ext-oracle),[],1);
    den=sum(oracle.^2,1); rel=sqrt(sum((ext-oracle).^2,1)./den); rel(den==0)=nan;
    D.selfFreshRelativeL2(local,:)=rel;
    for name=agent, D.(name{1})(local,:)=d.(name{1}).'; end
    nd=decodeNeuralActivity(y,alpha);
    for name=fieldnames(nd).', D.(name{1})(local,:)=nd.(name{1}).'; end
    for name=edge, D.(name{1})(:,:,local)=d.(name{1}); end
    for name={'selected','fresh','hidden','active'}, D.(name{1})(:,:,local)=d.(name{1}); end
    for i=1:n
        ids=d.targetIds(i,:); h=ismember(ids,state.mask(i,:));
        sets={ids(h),ids(~h),ids}; labels={'hidden','remaining','all'};
        for j=1:3
            a=angularSetMetrics(d.trueBearing(i,sets{j}));
            D.([labels{j} 'R1'])(local,i)=a.R1;
            D.([labels{j} 'R2'])(local,i)=a.R2;
            D.([labels{j} 'MaxGap'])(local,i)=a.maxGap;
        end
    end
    D.inputField(:,:,local)=ext;
    frames.position(:,:,local+1)=pos; frames.heading(:,local+1)=head; frames.u(:,:,local+1)=u;
end
segment.config=cfg; segment.startStep=recordStart; segment.endStep=endStep;
segment.data=D; segment.state=frames; segment.checkpoints=checkpoints;
segment.finalPosition=pos; segment.finalHeading=head; segment.finalNeuralState=u;
segment.finalAlpha=alpha; segment.finalState=state; segment.finalRng=rng;
cp=checkpoints;
end
function [go,lo]=orderMetrics(pos,ux,uy,L,k)
n=size(pos,2); go=hypot(sum(ux),sum(uy))/n; vals=zeros(n,1);
dx=pos(1,:).'-pos(1,:); dy=pos(2,:).'-pos(2,:);
dx=dx-L*round(dx/L); dy=dy-L*round(dy/L); dist=hypot(dx,dy);
dist(1:n+1:end)=0;
for i=1:n
    [~,ord]=sort(dist(i,:),'ascend'); ids=ord(1:min(k+1,n));
    vals(i)=hypot(sum(ux(ids)),sum(uy(ids)))/numel(ids);
end
lo=mean(vals);
end
function [pair,near]=distanceMetrics(pos,L)
n=size(pos,2); dx=pos(1,:).'-pos(1,:); dy=pos(2,:).'-pos(2,:);
dx=dx-L*round(dx/L); dy=dy-L*round(dy/L); d=hypot(dx,dy);
d(1:n+1:end)=nan; pair=mean(d,'all','omitnan'); near=mean(min(d,[],2,'omitnan'),'omitnan');
end
