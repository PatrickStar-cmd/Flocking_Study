function testBlockTimeMetrics()
thresholds=[.5 .5 .5]; good=repmat([1 1 0],2000,1);
[b,s]=blockTimeMetrics(good,.3,300,thresholds);
assert(height(b)==7 && b.nSteps(end)==200 && sum(b.nSteps)==2000);
assert(s.durationWeightedPass==1 && s.failureRuns==0);
bad=repmat([0 0 1],2000,1);
[~,s]=blockTimeMetrics(bad,.3,300,thresholds);
assert(s.durationWeightedPass==0 && s.longestFailureTime==600 && s.terminalFailureCensored);
tail=good; tail(1801:end,:)=bad(1801:end,:);
[~,s]=blockTimeMetrics(tail,.3,300,thresholds);
assert(abs(s.durationWeightedPass-.9)<1e-12 && s.terminalFailureCensored && s.recoveries==0);
short=good; short(501,:)=bad(501,:);
[~,s]=blockTimeMetrics(short,.3,100,thresholds);
assert(abs(s.pointFraction-.9995)<1e-12 && s.durationWeightedPass==1 && s.recoveries==1);
assert(~s.terminalFailureCensored && abs(s.longestFailureTime-.3)<1e-12);
assert(max(abs(s.weightedMeans-mean(short,1)))<1e-12);
cfg=defaultEffectiveInputConfig(); cfg.NSteps=40; cfg.burnIn=20; cfg.seed=1;
r=simulateEffectiveSocialInput(cfg); t=r.timeSeries;
m=[t.globalOrder(21:40),t.localOrder(21:40),t.meanPairDistanceNorm(21:40)];
[~,s]=blockTimeMetrics(m,cfg.dt,13,thresholds);
assert(max(abs(s.weightedMeans-[r.summary.globalOrder r.summary.localOrder r.summary.meanPairDistanceNorm]))<1e-12);
disp('PASS: tail weighting, success/failure, short failure, censoring and simulator summary regression');
end
