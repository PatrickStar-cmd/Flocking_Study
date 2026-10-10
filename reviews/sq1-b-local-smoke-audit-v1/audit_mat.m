% Read-only audit of B's saved data; all new outputs stay beside this script.
root=fileparts(mfilename('fullpath'));
src='D:/Work/Study/UESTC资料存档/课程/通关项目/paper/experiment/sq1_fourbranch_pilot_20260929';
pkg=fullfile(src,'results','sq1-fourbranch-pilot-v1-seed1-retry');
addpath(src,'-begin'); cleanup=onCleanup(@()rmpath(src));
fprintf('Loading cumulative hold MAT...\n'); s=load(fullfile(pkg,'hold_seed1.mat'),'runs'); runs=s.runs;
out=struct(); out.matlabVersion=version; out.savedCellCount=numel(runs); out.modes=cell(4,1);
for k=1:4
 r=runs{k}; D=r.data; first=r.records{6001}; last=r.records{6200};
 c=struct('mode',r.config.mode,'config',r.config,'recordCount',numel(r.records),'dataFields',{fieldnames(D)},'resultFields',{fieldnames(r)},...
 'allLOEqualsGO',isequaln(D.LO,D.GO),'phaseValues',unique(D.phase),'selected',unique(D.k_selected),...
 'eventFreshCount',unique(D.k_fresh(6001:end,:)),'eventActiveCount',unique(D.k_active(6001:end,:)),...
 'ageAllInfAt6001',all(isinf(first.age(:))),'ageAllInfAt6200',all(isinf(last.age(:))),...
 'prefixGOExact',isequal(D.GO(1:6000,:),runs{1}.data.GO(1:6000,:)),...
 'checkpoint6000Exact',isequal(r.checkpoints.step_6000.u,runs{1}.checkpoints.step_6000.u)&&isequal(r.checkpoints.step_6000.position,runs{1}.checkpoints.step_6000.position));
 cp=r.checkpoints.step_6200; pos=cp.position; L=r.config.arenaSize;
 dx=pos(1,:).'-pos(1,:);dy=pos(2,:).'-pos(2,:);dx=dx-L*round(dx/L);dy=dy-L*round(dy/L);dist=hypot(dx,dy);
 independent=mean(dist(triu(true(20),1)));c.finalPairDistance=independent;c.savedPairError=abs(independent-D.pairDistance(end,1));
 y=max(tanh(r.config.beta*cp.u),0);cx=sum(y.*cos(cp.alpha).',1);cy=sum(y.*sin(cp.alpha).',1);sp=hypot(cx,cy);ux=zeros(1,20);uy=ux; moving=sp*r.config.speedScale>1e-12;ux(moving)=cx(moving)./sp(moving);uy(moving)=cy(moving)./sp(moving);
 lv=zeros(20,1);for i=1:20,[~,order]=sort(dist(i,:));nb=order(1:6);lv(i)=hypot(sum(ux(nb)),sum(uy(nb)))/6;end
 c.independentFinalLO=mean(lv);c.savedFinalLO=D.LO(end,1);c.independentFinalGO=hypot(sum(ux),sum(uy))/20;
 c.eventMeanGO=mean(D.GO(6001:6200,1));c.eventMeanPairNormalized=mean(D.pairDistance(6001:6200,1))/L;
 c.eventPairedDistanceNormalized=mean(D.pairDistance(6001:6200,1)-runs{1}.data.pairDistance(6001:6200,1))/L;
 c.inputAmplitudeMean=mean(D.inputAmplitude(6001:6200,:),'all');
 c.holdCacheUnchanged=isequaln(first.injectedBearing(first.hidden),last.injectedBearing(last.hidden));
 c.savedRecordFields=fieldnames(first);
 out.modes{k}=c;
end
t=load(fullfile(pkg,'T0_results.mat'));out.T0=t;
% Short isolated calls expose amplitude, age and global-RNG behavior.
cfg=defaultFourBranchConfig(); cfg.occlusionStartStep=2;cfg.occlusionDurationSteps=2;
pos=runs{1}.checkpoints.step_6000.position; alpha=runs{1}.checkpoints.step_6000.alpha;
out.synthetic=cell(4,1);
for k=1:4
 cfg.mode=runs{k}.config.mode;rng(9876);before=rng;
 [~,state,~]=fourBranchInputStep(pos,alpha,[],cfg,1);after=rng;
 [field,state,d]=fourBranchInputStep(pos,alpha,state,cfg,2);
 % Fit field to independent unit-weight kernels to recover actual coefficient sum.
 ids=find(d.active(1,:));basis=zeros(100,1);
 for j=ids,delta=atan2(sin(alpha(1,:)-d.injectedBearing(1,j)),cos(alpha(1,:)-d.injectedBearing(1,j)));basis=basis+exp(-.5*(delta/.4).^2).';end
 coefficient=(basis'*field(:,1))/(basis'*basis);z1=mean(exp(1i*d.injectedBearing(1,ids)));z2=mean(exp(2i*d.injectedBearing(1,ids)));
 out.synthetic{k}=struct('mode',cfg.mode,'actualTotalWeight',coefficient*numel(ids),'globalRngChanged',~isequal(before,after),...
 'hiddenAges',d.age(d.hidden),'reportedZ1',d.z1(1),'reportedZ2',d.z2(1),'correctZ1Real',real(z1),'correctZ1Imag',imag(z1),'correctZ2Real',real(z2),'correctZ2Imag',imag(z2));
end
fid=fopen(fullfile(root,'mat_audit.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(out,PrettyPrint=true));fclose(fid);
fprintf('PASS audit completed: four saved modes read; findings written. This is not protocol acceptance.\n');
