% Independent read-only audit. Do not call B's file-writing runner/tests.
root=fileparts(mfilename('fullpath'));repo=fileparts(fileparts(root));
src='D:/Work/Study/UESTC资料存档/课程/通关项目/WHQ/sq1-fourbranch-protocol12-seed1-r2';
pkg=fullfile(src,'results','sq1-fourbranch-protocol12-seed1-r2-results-v3');
addpath(fullfile(repo,'effective-social-input'));addpath(src,'-begin');
out=struct('matlabVersion',version,'modes',[]); modes={'fresh','drop','drop-normalized','hold'};runs=cell(4,1);
for q=1:4,x=load(fullfile(pkg,[modes{q} '_seed1.mat']),'run');runs{q}=x.run;end
summary=readtable(fullfile(pkg,'summary.csv'),'TextType','string'); checks=cell(4,1);
for q=1:4
 r=runs{q}; d=r.segment.data;c=r.config;s=r.segment.state; ref=runs{1};cpEvent=r.checkpoint6200(1);
 a=struct('mode',modes{q},'steps',numel(d.step),'stateFrames',size(s.position,3),'hasCorrectStepRange',isequal(d.step,(6001:6600).'),...
 'commonFieldsExact',true,'commonStateExact',isequal(s.position(:,:,1:201),ref.segment.state.position(:,:,1:201))&&isequal(s.heading(:,1:201),ref.segment.state.heading(:,1:201)),...
 'checkpoint6000Exact',isequaln(r.checkpoint6000,ref.checkpoint6000));
 a.checkpoint6200ElementCount=numel(r.checkpoint6200);
 names=fieldnames(d);for z=1:numel(names)
  name=names{z};v=d.(name);vref=ref.segment.data.(name);
  if ndims(v)==3, equal=isequaln(v(:,:,1:200),vref(:,:,1:200));else,equal=isequaln(v(1:200,:),vref(1:200,:));end
  a.commonFieldsExact=a.commonFieldsExact&&equal;
 end
 a.Aevent=unique(d.A(201:234,:));a.Aafter=unique(d.A(235:end,:));a.eventFresh=unique(d.k_fresh(201:234,:));a.eventActive=unique(d.k_active(201:234,:));
 a.recoveryFresh=unique(d.k_fresh(235:end,:));a.selectedAlways8=all(d.k_selected(:)==8);a.hiddenObservedNaN=all(isnan(d.observedBearing(d.hidden)));
 a.phaseCorrect=all(d.phase(1:200)=="common_prefix")&&all(d.phase(201:234)=="occlusion")&&all(d.phase(235:end)=="recovery");
 a.frozenMaskExact=isequal(r.checkpoint6200.mask,ref.checkpoint6200.mask);
 a.ageFirst=[];a.ageLast=[];a.ageRecovery=[];a.holdCacheError=0;
 if q>1
  mask=d.hidden(:,:,201);v=d.age(:,:,201);a.ageFirst=unique(v(mask));v=d.age(:,:,234);a.ageLast=unique(v(mask));v=d.age(:,:,235);a.ageRecovery=unique(v(mask));
  if q==4
   for j=201:234,v=d.injectedBearing(:,:,j);cache=cpEvent.lastBearing;a.holdCacheError=max(a.holdCacheError,max(abs(v(mask)-cache(mask))));end
  end
 end
 a.maxARecomputeError=0;a.maxMomentError=0;a.maxWeightEntropyError=0;a.maxMassError=0;a.maxFieldError=0;a.maxOracleError=0;
 for j=1:600
  w=d.edgeWeight(:,:,j);th=d.injectedBearing(:,:,j);th(~d.active(:,:,j))=0;A=sum(w,2);
  z1=sum(w.*exp(1i*th),2)./A;z2=sum(w.*exp(2i*th),2)./A;
  a.maxARecomputeError=max(a.maxARecomputeError,max(abs(A-d.A(j,:).')));
  a.maxMomentError=max(a.maxMomentError,max(abs([real(z1)-d.z1_real(j,:).';imag(z1)-d.z1_imag(j,:).';real(z2)-d.z2_real(j,:).';imag(z2)-d.z2_imag(j,:).'])));
  p=w./A;p(p==0)=1; entropy=-sum((w./A).*log(p),2);a.maxWeightEntropyError=max(a.maxWeightEntropyError,max(abs(entropy-d.weightEntropy(j,:).')));
  a.maxMassError=max(a.maxMassError,max(abs(sum(d.inputField(:,:,j),1)*2*pi/100-d.Mass(j,:))));
 end
 for j=[1 200 201 234 235 334 600]
  field=zeros(100,20);oracle=field;alpha=r.checkpoint6000.alpha;
  for i=1:20
   for k=find(d.active(i,:,j))
    delta=atan2(sin(alpha(i,:)-d.injectedBearing(i,k,j)),cos(alpha(i,:)-d.injectedBearing(i,k,j)));
    field(:,i)=field(:,i)+d.edgeWeight(i,k,j)*exp(-.5*(delta/.4).^2).';
   end
   ids=cpEvent.targetIds(i,:);
   for k=ids,delta=atan2(sin(alpha(i,:)-d.trueBearing(i,k,j)),cos(alpha(i,:)-d.trueBearing(i,k,j)));oracle(:,i)=oracle(:,i)+c.w0*exp(-.5*(delta/.4).^2).';end
  end
  a.maxFieldError=max(a.maxFieldError,max(abs(field-d.inputField(:,:,j)),[],'all'));a.maxOracleError=max(a.maxOracleError,max(abs(oracle-d.selfFreshField(:,:,j)),[],'all'));
 end
 a.maxPairError=0;a.maxNearestError=0;a.maxGOFromMotionError=0;a.maxLOFromMotionError=0;
 for j=1:600
  pos=s.position(:,:,j+1);prev=s.position(:,:,j);L=c.arenaSize;
  dx=pos(1,:).'-pos(1,:);dy=pos(2,:).'-pos(2,:);dx=dx-L*round(dx/L);dy=dy-L*round(dy/L);dist=hypot(dx,dy);
  pair=mean(dist(triu(true(20),1)))/L;near=dist;near(1:21:end)=inf;near=mean(min(near,[],2))/L;
  a.maxPairError=max(a.maxPairError,abs(pair-d.pairDistance(j)));a.maxNearestError=max(a.maxNearestError,abs(near-d.nearestDistance(j)));
  delta=pos-prev;delta=delta-L*round(delta/L);v=delta/c.dt;speed=hypot(v(1,:),v(2,:));ux=zeros(1,20);uy=ux;nz=speed>1e-12;ux(nz)=v(1,nz)./speed(nz);uy(nz)=v(2,nz)./speed(nz);
  go=hypot(sum(ux),sum(uy))/20;lv=zeros(20,1);
  for i=1:20,[~,idx]=sort(dist(i,:));idx=idx(1:6);lv(i)=hypot(sum(ux(idx)),sum(uy(idx)))/6;end
  a.maxGOFromMotionError=max(a.maxGOFromMotionError,abs(go-d.GO(j)));a.maxLOFromMotionError=max(a.maxLOFromMotionError,abs(mean(lv)-d.LO(j)));
 end
 ix=235:334;a.primaryDistance=mean(d.pairDistance(ix)-ref.segment.data.pairDistance(ix));a.primaryGOLoss=mean(ref.segment.data.GO(ix)-d.GO(ix));a.primaryLOLoss=mean(ref.segment.data.LO(ix)-d.LO(ix));
 row=find(summary.mode==string(modes{q}));a.summaryMaxError=max(abs([mean(d.GO(ix))-summary.meanGO(row),mean(d.LO(ix))-summary.meanLO(row),mean(d.pairDistance(ix))-summary.meanPairDistanceNorm(row)]));
 a.missingFields=setdiff({'lastBearing','lastWeight','observedWeight','neuralMass','bumpConcentration','decoded_direction_defined'},names);checks{q}=a;
end
out.modes=checks;fprintf('Saved data and input/geometry recomputation complete.\n');
x=load(fullfile(pkg,'T0_results.mat'),'t0');tn=fieldnames(x.t0);out.T0=cell(numel(tn),1);
for q=1:numel(tn),z=x.t0.(tn{q});base=x.t0.fresh;out.T0{q}=struct('mode',tn{q},'allDataExact',isequaln(z.data,base.data),'stateExact',isequaln(z.state,base.state),'finalUExact',isequaln(z.finalNeuralState,base.finalNeuralState),'finalCacheExact',isequaln(z.finalState,base.finalState));end
clear x;
% Re-run only the 600-step fresh tail from B's checkpoint using the supplied source.
c=runs{1}.config;[tail,~]=resumeProtocolFromCheckpoint(runs{1}.checkpoint6000,c,6600);
out.resume=struct('allDataExact',isequaln(tail.data,runs{1}.segment.data),'allStateExact',isequaln(tail.state,runs{1}.segment.state),'finalUExact',isequaln(tail.finalNeuralState,runs{1}.segment.finalNeuralState),'finalCacheExact',isequaln(tail.finalState,runs{1}.segment.finalState));
out.resume.maxPositionDiff=max(abs(tail.state.position-runs{1}.segment.state.position),[],'all');out.resume.maxInputDiff=max(abs(tail.data.inputField-runs{1}.segment.data.inputField),[],'all');out.resume.maxUDiff=max(abs(tail.finalNeuralState-runs{1}.segment.finalNeuralState),[],'all');out.resume.maxGODiff=max(abs(tail.data.GO-runs{1}.segment.data.GO));
out.resume.firstInputDiff=max(abs(tail.data.inputField(:,:,1)-runs{1}.segment.data.inputField(:,:,1)),[],'all');out.resume.firstPositionDiff=max(abs(tail.state.position(:,:,2)-runs{1}.segment.state.position(:,:,2)),[],'all');
out.resume.numericFieldDiffs=struct();for z=fieldnames(tail.data).',name=z{1};if isnumeric(tail.data.(name)),aa=tail.data.(name);bb=runs{1}.segment.data.(name);ix=isfinite(aa)&isfinite(bb);out.resume.numericFieldDiffs.(name)=max(abs(aa(ix)-bb(ix)));end,end
% Initialization is compared independently to the preserved old core.
c=defaultFourBranchConfig();s=initializeCoreState(c);rng(1);pos=rand(2,20)*1000;head=zeros(20,1);u=zeros(100,20);for i=1:20,head(i)=2*pi*rand;u(:,i)=.1*randn(100,1);end
out.initializationExact=isequal(pos,s.position)&&isequal(head,s.heading)&&isequal(u,s.u);
% A short same-environment repeat isolates determinism from cross-version data comparisons.
[a1,~]=simulateProtocolSegment(c,s,0,3);[a2,~]=simulateProtocolSegment(c,s,0,3);out.sameEnvironmentRepeatExact=isequaln(a1,a2);
% Demonstrate ID tie handling with equal bearings around a cyclic wrap.
pos=zeros(2,20);pos(1,17)=-1;alpha=s.alpha;cfg=c;cfg.maskRule='max_hidden_R1';cfg.maskReadyStep=1;
rng(5678);before=rng;[~,maskstate,~]=fourBranchInputStep(pos,alpha,[],cfg,1);after=rng;out.maskGlobalRngUnchanged=isequal(before,after);
out.tieObserver=17;out.tieActual=maskstate.mask(17,:);out.tieExpected=[1 2 3 4];
% k=0 adapter/indicator behavior, without any trajectory integration.
cfg=c;cfg.k=0;cfg.neighborCount=0;cfg.h=0;cfg.maskReadyStep=1;
try,[field,~,zero]=fourBranchInputStep(s.position,s.alpha,[],cfg,1);out.zeroInput=struct('executes',true,'fieldZero',all(field(:)==0),'momentsAllNaN',all(isnan(zero.z1)),'AZero',all(zero.A==0));catch ME,out.zeroInput=struct('executes',false,'message',ME.message);end
% Verify phase labels actually follow a changed, already registered duration.
cfg=c;cfg.occlusionDurationSteps=10;cfg.recoveryStartStep=6211;
try,[~,~]=resumeProtocolFromCheckpoint(runs{1}.checkpoint6200,cfg,6212);out.eventCheckpointResume=struct('executes',true);catch ME,out.eventCheckpointResume=struct('executes',false,'message',ME.message);end
% Build a scalar start from the recorded event state to isolate the phase issue.
ev=runs{1}.checkpoint6200(1);st=struct();st.position=ev.position;st.heading=ev.heading;st.u=ev.u;st.alpha=ev.alpha;st.rng=ev.rng;st.step=ev.step;st.state=runs{1}.segment.finalState;
st.state.lastBearing=ev.lastBearing;st.state.lastWeight=ev.lastWeight;st.state.age=ev.age;st.state.mask=ev.mask;st.state.maskReady=ev.maskReady;st.state.driverIndex=ev.driverIndex;
[short,~]=simulateProtocolSegment(cfg,st,6200,6212);out.duration10RecoveryPhase=char(short.data.phase(11));
fid=fopen(fullfile(root,'mat_audit.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(out,PrettyPrint=true));fclose(fid);
fprintf('Audit completed. Not a full G1-G15 acceptance or a scientific experiment.\n');
