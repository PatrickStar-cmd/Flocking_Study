function analyzeTimeSeriesValidation()
root=fileparts(mfilename('fullpath'));
src=fullfile(root,'results','time-series-validation-20260909');
out=fullfile(src,'analysis-v2'); if ~exist(out,'dir'), mkdir(out); end
meta=load(fullfile(src,'metadata.mat'));
assert(numel(dir(fullfile(src,'seed-*-condition-*.mat')))==180);
values=zeros(30,6,7); details=cell(180,13); index=0; maxError=0;
for si=1:30
 for c=1:6
  seed=100+si;
  x=load(fullfile(src,sprintf('seed-%03d-condition-%d.mat',seed,c)));
  assert(x.complete && x.result.config.seed==seed);
  cfg=meta.base; cfg.seed=seed; cfg.selectionPolicy=meta.conditions{c,1};
  cfg.strengthMode=meta.conditions{c,2}; cfg.neighborCount=meta.conditions{c,3};
  assert(isequaln(x.result.config,cfg) && isequal(x.thresholds,meta.thresholds));
  ts=x.result.timeSeries;
  names=fieldnames(ts);
  for j=1:numel(names), assert(numel(ts.(names{j}))==8000 && all(isfinite(ts.(names{j})))); end
  m=[ts.globalOrder(6001:8000),ts.localOrder(6001:8000),ts.meanPairDistanceNorm(6001:8000)];
  [b1,s1]=blockTimeMetrics(m,cfg.dt,100,x.thresholds);
  [b3,s3]=blockTimeMetrics(m,cfg.dt,300,x.thresholds);
  assert(isequal(b1,x.blocks100) && isequal(b3,x.blocks300));
  assert(isequaln(s1,x.summary100) && isequaln(s3,x.summary300));
  assert(sum(b1.nSteps)==2000 && sum(b3.nSteps)==2000 && b3.nSteps(end)==200);
  old=[x.result.summary.globalOrder,x.result.summary.localOrder,x.result.summary.meanPairDistanceNorm];
  maxError=max(maxError,max(abs(mean(m,1)-old)));
  values(si,c,:)=[old,s1.pointFraction,s1.durationWeightedPass,s3.durationWeightedPass,s1.longestFailureTime];
  index=index+1;
  details(index,:)={seed,c,old(1),old(2),old(3),s1.pointFraction,s1.durationWeightedPass,s3.durationWeightedPass,s1.longestFailureTime,s1.failureRuns,s1.recoveries,s1.initialFailureLeftCensored,s1.terminalFailureCensored};
 end
end
assert(maxError<1e-10);
writetable(cell2table(details,'VariableNames',{'seed','condition','GO','LO','distance','pointFraction','block100','block300','longestFailureTime','failureRuns','recoveries','leftCensored','rightCensored'}),fullfile(out,'runs.csv'));
rng(20260909,'twister'); samples=randi(30,30,5000);
metrics={'GO','LO','distance','pointFraction','block100','block300','longestFailureTime'};
groups=cell(42,6); idx=0;
for c=1:6
 for j=1:7
  v=values(:,c,j); boot=mean(v(samples),1); ci=prctile(boot,[2.5 97.5]);
  idx=idx+1; groups(idx,:)={c,metrics{j},mean(v),median(v),ci(1),ci(2)};
 end
end
writetable(cell2table(groups,'VariableNames',{'condition','metric','mean','median','lower95','upper95'}),fullfile(out,'groups.csv'));
% All five contrasts with baseline; plus prespecified mechanistic comparisons.
contrasts=[2 1;3 1;4 1;5 1;6 1;4 3;5 4;4 6];
paired=cell(size(contrasts,1)*7,6); idx=0;
for p=1:size(contrasts,1)
 for j=1:7
  d=values(:,contrasts(p,1),j)-values(:,contrasts(p,2),j);
  ci=prctile(mean(d(samples),1),[2.5 97.5]);
  idx=idx+1; paired(idx,:)={contrasts(p,1),contrasts(p,2),metrics{j},mean(d),ci(1),ci(2)};
 end
end
writetable(cell2table(paired,'VariableNames',{'conditionA','conditionB','metric','meanDifference','lower95','upper95'}),fullfile(out,'paired.csv'));
audit.files=180; audit.seeds=101:130; audit.maxSummaryError=maxError;
audit.bootstrapReplicates=5000; audit.bootstrapSeed=20260909;
audit.note='Pointwise exploratory percentile intervals; seed is resampling unit. Not adjusted for multiple comparisons; thresholds fixed.';
fid=fopen(fullfile(out,'audit.json'),'w'); fprintf(fid,'%s',jsonencode(audit)); fclose(fid);
disp(cell2table(groups,'VariableNames',{'condition','metric','mean','median','lower95','upper95'}));
end
