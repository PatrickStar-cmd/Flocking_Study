% Independent calculations on the ORIGINAL delivered MAT, not C's regenerated files.
root=fileparts(mfilename('fullpath'));src=fullfile(root,'source');pkg=fullfile(src,'results','sq1-fourbranch-protocol12-seed1-r3-results-v4');
out=struct();modes={'fresh','drop','drop-normalized','hold'};runs=cell(4,1);checks=cell(4,1);
for q=1:4,x=load(fullfile(pkg,[modes{q} '_seed1.mat']),'run');runs{q}=x.run;end
tab=readtable(fullfile(pkg,'summary.csv'),'TextType','string');
for q=1:4
 r=runs{q};d=r.segment.data;s=r.segment.state;c=r.config;f=runs{1};
 a.mode=modes{q};a.steps=numel(d.step);a.frames=size(s.position,3);a.cp6000Scalar=isscalar(r.checkpoint6000);a.cp6200Scalar=isscalar(r.checkpoint6200);a.maskDiagnosticsCount=numel(r.checkpoint6200.maskDiagnostics);
 a.commonPositionExact=isequaln(s.position(:,:,1:201),f.segment.state.position(:,:,1:201));a.commonUExact=isequaln(s.u(:,:,1:201),f.segment.state.u(:,:,1:201));
 a.Aevent=unique(d.A(201:234,:));a.entropyEvent=unique(d.weightEntropy(201:234,:));a.phaseCorrect=all(d.phase(1:200)=="common_prefix")&&all(d.phase(201:234)=="occlusion")&&all(d.phase(235:end)=="recovery");
 a.maxPairError=0;a.maxGOError=0;a.maxLOError=0;a.maxFieldError=0;a.maxMomentError=0;a.maxSelfL2Error=0;a.maxSummaryError=0;a.maxHoldCacheError=0;
 for t=1:600
  p=s.position(:,:,t+1);pr=s.position(:,:,t);L=c.arenaSize;
  delta=p-pr;delta=delta-L*round(delta/L);v=delta/c.dt;speed=hypot(v(1,:),v(2,:));unit=zeros(2,20);nz=speed>1e-12;unit(:,nz)=v(:,nz)./speed(nz);
  dx=p(1,:).'-p(1,:);dy=p(2,:).'-p(2,:);dx=dx-L*round(dx/L);dy=dy-L*round(dy/L);dist=hypot(dx,dy);
  pair=mean(dist(triu(true(20),1)))/L;go=norm(sum(unit,2))/20;lv=zeros(20,1);
  for i=1:20,[~,idx]=sort(dist(i,:));lv(i)=norm(sum(unit(:,idx(1:6)),2))/6;end
  a.maxPairError=max(a.maxPairError,abs(pair-d.pairDistance(t)));a.maxGOError=max(a.maxGOError,abs(go-d.GO(t)));a.maxLOError=max(a.maxLOError,abs(mean(lv)-d.LO(t)));
  w=d.edgeWeight(:,:,t);theta=d.injectedBearing(:,:,t);theta(~d.active(:,:,t))=0;A=sum(w,2);z1=sum(w.*exp(1i*theta),2)./A;z2=sum(w.*exp(2i*theta),2)./A;
  a.maxMomentError=max(a.maxMomentError,max(abs([real(z1)-d.z1_real(t,:).';imag(z1)-d.z1_imag(t,:).';real(z2)-d.z2_real(t,:).';imag(z2)-d.z2_imag(t,:).'])));
  den=sum(d.selfFreshField(:,:,t).^2,1);l2=sqrt(sum((d.inputField(:,:,t)-d.selfFreshField(:,:,t)).^2,1)./den);a.maxSelfL2Error=max(a.maxSelfL2Error,max(abs(l2-d.selfFreshRelativeL2(t,:))));
  if ismember(t,[1 200 201 234 235 334 600])
   field=zeros(100,20);
   for i=1:20,for j=find(d.active(i,:,t)),delta=atan2(sin(s.alpha(i,:)-d.injectedBearing(i,j,t)),cos(s.alpha(i,:)-d.injectedBearing(i,j,t)));field(:,i)=field(:,i)+w(i,j)*exp(-delta.^2/(2*.4^2)).';end,end
   a.maxFieldError=max(a.maxFieldError,max(abs(field-d.inputField(:,:,t)),[],'all'));
  end
  if q==4 && t>=201 && t<=234
   mask=d.hidden(:,:,t);actual=d.injectedBearing(:,:,t);cache=r.checkpoint6200.lastBearing;a.maxHoldCacheError=max(a.maxHoldCacheError,max(abs(actual(mask)-cache(mask))));
  end
 end
 rows=find(tab.mode==string(modes{q}));
 for idx=rows.'
  ix=d.step>=tab.windowStart(idx)&d.step<=tab.windowEnd(idx);ref=f.segment.data;
  expected=[mean(d.GO(ix)),mean(d.LO(ix)),mean(d.pairDistance(ix)),mean(d.pairDistance(ix)-ref.pairDistance(ix)),mean(ref.GO(ix)-d.GO(ix)),mean(ref.LO(ix)-d.LO(ix))];
  saved=[tab.meanGO(idx),tab.meanLO(idx),tab.meanPairDistanceNorm(idx),tab.pairDistanceIncrease(idx),tab.GOLoss(idx),tab.LOLoss(idx)];a.maxSummaryError=max(a.maxSummaryError,max(abs(expected-saved)));
 end
 ix=235:334;a.primaryDistance=mean(d.pairDistance(ix)-f.segment.data.pairDistance(ix));checks{q}=a;
end
out.modes=checks;
x=load(fullfile(pkg,'reference_same_state_replay.mat'),'replay');rp=x.replay;out.sameState=cell(4,1);
for q=1:4
 cand=rp.candidates.(matlab.lang.makeValidName(modes{q}));fr=rp.candidates.fresh;
 rel=sqrt(sum((cand.field-fr.field).^2,1)./sum(fr.field.^2,1));rel=squeeze(rel).';
 out.sameState{q}=struct('mode',modes{q},'maxL2RecomputeError',max(abs(rel-cand.relativeL2),[],'all'),'commonPrefixExact',isequal(cand.field(:,:,1:200),fr.field(:,:,1:200)),'meanEventL2',mean(cand.relativeL2(201:234,:),'all'));
end
x=load(fullfile(pkg,'T0_results.mat'),'t0');base=x.t0.fresh;out.T0=cell(4,1);
for q=1:4,b=x.t0.(matlab.lang.makeValidName(modes{q}));out.T0{q}=struct('mode',modes{q},'dataExact',isequaln(base.data,b.data),'stateExact',isequaln(base.state,b.state),'cacheExact',isequaln(base.finalState,b.finalState),'rngExact',isequaln(base.finalRng,b.finalRng));end
fid=fopen(fullfile(root,'saved_mat_audit.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(out,PrettyPrint=true));fclose(fid);
fprintf('Original delivered MAT: independent 2400-step distances/GO/LO/moments, 560 agent-step fields, 16 summary rows and two replay types checked.\n');
