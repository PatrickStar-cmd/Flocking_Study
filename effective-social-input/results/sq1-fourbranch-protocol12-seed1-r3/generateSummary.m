function tab=generateSummary(outDir)
%GENERATESUMMARY Window means and signed paired contrasts from saved steps.
modes={'fresh','drop','drop-normalized','hold'}; x=load(fullfile(outDir,'fresh_seed1.mat'),'run');
ref=x.run.segment.data; cfg=x.run.config;
s=cfg.occlusionStartStep; d=cfg.occlusionDurationSteps; r=s+d;
windowNames={'occlusion','recovery100','recovery300','event400'};
bounds=[s s+d-1;r r+99;r r+299;s 6600]; rows=struct([]); count=0;
for q=1:4
    x=load(fullfile(outDir,[modes{q} '_seed1.mat']),'run'); a=x.run.segment.data;
    for w=1:4
        ix=a.step>=bounds(w,1)&a.step<=bounds(w,2); count=count+1;
        row.run_id=[modes{q} '_seed1_random_d34']; row.reference_id='fresh_seed1';
        row.mode=modes{q}; row.seed=cfg.seed; row.mask_rule=cfg.maskRule; row.duration_steps=d;
        row.window=windowNames{w}; row.windowStart=bounds(w,1); row.windowEnd=bounds(w,2); row.nSteps=sum(ix);
        row.meanGO=mean(a.GO(ix)); row.meanLO=mean(a.LO(ix)); row.meanPairDistanceNorm=mean(a.pairDistance(ix));
        row.meanActive=mean(a.k_active(ix,:),'all'); row.meanA=mean(a.A(ix,:),'all');
        row.pairDistanceIncrease=mean(a.pairDistance(ix)-ref.pairDistance(ix));
        row.GOLoss=mean(ref.GO(ix)-a.GO(ix)); row.LOLoss=mean(ref.LO(ix)-a.LO(ix));
        row.selfFreshL2Mean=mean(a.selfFreshRelativeL2(ix,:),'all','omitnan');
        row.selfFreshL2Max=max(a.selfFreshRelativeL2(ix,:),[],'all');
        if count==1, rows=row; else, rows(count)=row; end
    end
end
tab=struct2table(rows);
for metric={'pairDistanceIncrease','GOLoss','LOLoss'}
    name=metric{1}; normMinusDrop=nan(height(tab),1); holdMinusNorm=normMinusDrop;
    for w=1:4
        b=tab.window==string(windowNames{w});
        normMinusDrop(b)=tab.(name)(b & tab.mode=="drop-normalized")-tab.(name)(b & tab.mode=="drop");
        holdMinusNorm(b)=tab.(name)(b & tab.mode=="hold")-tab.(name)(b & tab.mode=="drop-normalized");
    end
    tab.(['normalizedMinusDrop_' name])=normMinusDrop; tab.(['holdMinusNormalized_' name])=holdMinusNorm;
end
end
