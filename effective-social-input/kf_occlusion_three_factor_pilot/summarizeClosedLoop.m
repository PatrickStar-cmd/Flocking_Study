function rows=summarizeClosedLoop(seg,reference,cfg,condition,runid)
% 所有群体偏离相对于同seed、同A0的fresh参考计算。
steps=seg.data.step; [ok,ri]=ismember(steps,reference.data.step); assert(all(ok));
d=seg.data; ref=reference.data; s=cfg.occlusionStartStep; r=cfg.recoveryStartStep;
bounds=[s r-1;r r+99;r r+299;s cfg.NSteps]; labels={'occlusion','recovery100','recovery300','event400'};
for w=1:4
    ix=steps>=bounds(w,1) & steps<=bounds(w,2); jj=ri(ix);
    row.run_id=runid; row.reference_id=sprintf('fresh_seed%d_A%d',cfg.seed,round(cfg.A0*1000));
    row.seed=cfg.seed; row.condition_id=condition.id; row.h=cfg.h; row.structure=condition.label;
    row.A0=cfg.A0; row.mode=cfg.mode; row.window=labels{w}; row.windowStart=bounds(w,1);
    row.windowEnd=bounds(w,2); row.nSteps=sum(ix);
    row.meanGO=mean(d.GO(ix)); row.meanLO=mean(d.LO(ix)); row.meanPairDistanceNorm=mean(d.pairDistance(ix));
    row.pairDistanceIncrease=mean(d.pairDistance(ix)-ref.pairDistance(jj));
    row.GOLoss=mean(ref.GO(jj)-d.GO(ix)); row.LOLoss=mean(ref.LO(jj)-d.LO(ix));
    row.meanActive=mean(d.k_active(ix,:),'all'); row.meanA=mean(d.A(ix,:),'all');
    row.meanMass=mean(d.Mass(ix,:),'all'); row.meanPeak=mean(d.Peak(ix,:),'all');
    row.selfFreshL2Mean=mean(d.selfFreshRelativeL2(ix,:),'all','omitnan');
    row.selfFreshL2Max=max(d.selfFreshRelativeL2(ix,:),[],'all');
    hidden=d.hidden(:,:,ix); e=d.predictionError(:,:,ix); v=d.predictionVariance(:,:,ix);
    row.hiddenBearingRMSE=sqrt(mean(e(hidden).^2,'omitnan'));
    row.hiddenCoverage95=mean(abs(e(hidden))<=1.96*sqrt(v(hidden)));
    row.hiddenMeanVariance=mean(v(hidden),'omitnan');
    rows(w)=row; %#ok<AGROW>
end
end
