function analyzeDistanceWeightPilot()
root=fileparts(mfilename('fullpath'));
src=fullfile(root,'results','distance-weight-pilot-v1');
out=fullfile(src,'analysis'); if ~exist(out,'dir'), mkdir(out); end
meta=load(fullfile(src,'metadata.mat'));
assert(numel(dir(fullfile(src,'seed-*.mat')))==40);
metrics={'GO','LO','distance','amplitude','inputMass','distanceWeight'};
values=zeros(5,2,4,6); rows=cell(40,10); idx=0; err=0;
for s=1:5
 for ki=1:2
  k=meta.ks(ki);
  for c=1:4
   seed=meta.seeds(s);
   x=load(fullfile(src,sprintf('seed-%d-k-%d-condition-%d.mat',seed,k,c)));
   assert(x.complete);
   cfg=meta.base; cfg.seed=seed; cfg.selectionPolicy='random'; cfg.neighborCount=k;
   cfg.distanceWeightMode=meta.conditions{c,1}; cfg.strengthMode=meta.conditions{c,2};
   assert(isequaln(cfg,x.result.config));
   t=x.result.timeSeries; names=fieldnames(t);
   for j=1:numel(names)
    z=t.(names{j}); assert(numel(z)==8000 && all(isfinite(z)));
    err=max(err,abs(mean(z(6001:8000))-x.result.summary.(names{j})));
   end
   assert(all(t.meanActiveNeighbors==k));
   assert(all(t.meanDistanceWeight>=0 & t.meanDistanceWeight<=1));
   if c>=3, assert(max(abs(t.meanInputAmplitude-.228))<1e-12); end
   if c==1, assert(max(abs(t.meanInputAmplitude-.012*k))<1e-12); end
   q=[mean(t.globalOrder(6001:8000)),mean(t.localOrder(6001:8000)),mean(t.meanPairDistanceNorm(6001:8000)),mean(t.meanInputAmplitude(6001:8000)),mean(t.meanInputMass(6001:8000)),mean(t.meanDistanceWeight(6001:8000))];
   values(s,ki,c,:)=q;
   idx=idx+1; rows(idx,:)=[{seed,k,string(meta.conditions{c,1}),string(meta.conditions{c,2})},num2cell(q)];
  end
 end
end
assert(err<1e-12);
writetable(cell2table(rows,'VariableNames',[{'seed','k','weightMode','strengthMode'},metrics]),fullfile(out,'runs.csv'));
groups=cell(48,7); idx=0;
for ki=1:2
 for c=1:4
  for j=1:6
   v=values(:,ki,c,j); idx=idx+1;
   groups(idx,:)={meta.ks(ki),string(meta.conditions{c,1}),string(meta.conditions{c,2}),metrics{j},mean(v),min(v),max(v)};
  end
 end
end
writetable(cell2table(groups,'VariableNames',{'k','weightMode','strengthMode','metric','mean','min','max'}),fullfile(out,'groups.csv'));
% Exact enumeration of all 5^5 paired seed bootstrap resamples. Exploratory.
[a,b,c,d,e]=ndgrid(1:5); samples=[a(:),b(:),c(:),d(:),e(:)];
contrasts=[2 1;4 3;4 2]; pairs=cell(36,8); idx=0;
for ki=1:2
 for p=1:3
  for j=1:6
   delta=values(:,ki,contrasts(p,1),j)-values(:,ki,contrasts(p,2),j);
   ci=prctile(mean(delta(samples),2),[2.5 97.5]); idx=idx+1;
   pairs(idx,:)={meta.ks(ki),contrasts(p,1),contrasts(p,2),metrics{j},mean(delta),ci(1),ci(2),sum(delta>0)};
  end
 end
end
writetable(cell2table(pairs,'VariableNames',{'k','conditionA','conditionB','metric','meanDifference','bootstrapLower95','bootstrapUpper95','positivePairs'}),fullfile(out,'paired.csv'));
audit.files=40; audit.seeds=meta.seeds; audit.maxSummaryError=err;
audit.note='Five seeds per condition. Bootstrap intervals enumerate 3125 paired resamples and are exploratory, unadjusted, not reliability certification.';
fid=fopen(fullfile(out,'audit.json'),'w'); fprintf(fid,'%s',jsonencode(audit)); fclose(fid);
disp('PASS: 40 configurations, full time series, summaries, neighbor counts and amplitude invariants');
disp(cell2table(groups,'VariableNames',{'k','weightMode','strengthMode','metric','mean','min','max'}));
end
