function screenBearingPrediction()
% Offline screening, not new swarm dynamics or per-step calibration.
root=fileparts(mfilename('fullpath')); src=fullfile(root,'results','time-series-validation-20260909');
out=fullfile(root,'results','bearing-screen-v1'); if ~exist(out,'dir'), mkdir(out); end
qs=[1e-8 1e-6 1e-4]; R=.01^2; horizons=[1 2 4]; starts=[41 61];
rows=cell(30*2*20*3*3,9); idx=0; rates=[];
for seed=101:130
 x=load(fullfile(src,sprintf('seed-%d-condition-1.mat',seed)),'result');
 p=x.result.trajectory; dt=p.time(2)-p.time(1); assert(abs(dt-30)<1e-10);
 rng(seed+9000,'twister');
 for i=1:20
  j=mod(i,20)+1;
  dx=p.x(j,:)-p.x(i,:); dy=p.y(j,:)-p.y(i,:);
  L=x.result.config.arenaSize; dx=dx-L*round(dx/L); dy=dy-L*round(dy/L);
  theta=atan2(dy,dx); obs=theta+.01*randn(size(theta));
  if seed<=115
   d=diff(theta(1:61)); rates=[rates abs(atan2(sin(d),cos(d)))/dt]; %#ok<AGROW>
  end
  for qi=1:3
   track=[];
   for ti=1:max(starts)
    track=bearingTrackStep(track,obs(ti),dt,qs(qi),R);
    if ~ismember(ti,starts), continue; end
    predicted=track;
    for h=1:max(horizons)
     predicted=bearingTrackStep(predicted,[],dt,qs(qi),R);
     if ~ismember(h,horizons), continue; end
     truth=theta(ti+h);
     ek=atan2(sin(predicted.x(1)-truth),cos(predicted.x(1)-truth));
     eh=atan2(sin(obs(ti)-truth),cos(obs(ti)-truth));
     idx=idx+1; rows(idx,:)={seed,i,ti,qs(qi),h*dt,ek^2,eh^2,predicted.P(1,1),abs(ek)<=1.96*sqrt(predicted.P(1,1))};
    end
   end
  end
 end
end
t=cell2table(rows,'VariableNames',{'seed','focal','startIndex','q','T','kalmanSquaredError','holdSquaredError','variance','covered'});
loss=zeros(1,3);
for i=1:3, loss(i)=mean(t.kalmanSquaredError(t.seed<=115 & t.q==qs(i))); end
[~,best]=min(loss); selectedQ=qs(best);
summary=cell(3,6);
for h=1:3
 mask=t.seed>115 & t.q==selectedQ & t.T==horizons(h)*30;
 v=t(mask,:); perSeed=zeros(15,1);
 for i=1:15
  s=v.seed==115+i; perSeed(i)=mean(v.kalmanSquaredError(s)-v.holdSquaredError(s));
 end
 summary(h,:)={horizons(h)*30,sqrt(mean(v.kalmanSquaredError)),sqrt(mean(v.holdSquaredError)),mean(v.covered),sum(perSeed<0),mean(v.variance)};
end
writetable(t,fullfile(out,'errors.csv'));
report=cell2table(summary,'VariableNames',{'T','kalmanRMSE','holdRMSE','nominal95Coverage','seedsImprovedOf15','meanVariance'});
writetable(report,fullfile(out,'summary.csv'));
meta.selectedQ=selectedQ; meta.candidates=qs; meta.trainingMSE=loss; meta.R=R;
meta.dt=30; meta.angularTimeScale=.4/median(rates); meta.trainSeeds=101:115; meta.evaluationSeeds=116:130;
meta.note='Historical coarse trajectories; cyclic 20 directed edges, two cutoffs per seed; edges are not independent replicates. Noise is synthetic. Not dt=.3 calibration.';
fid=fopen(fullfile(out,'metadata.json'),'w'); fprintf(fid,'%s',jsonencode(meta)); fclose(fid);
disp(meta); disp(report);
end
