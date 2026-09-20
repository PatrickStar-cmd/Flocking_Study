function runTimeSeriesValidation()
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','time-series-validation-20260909');
if ~exist(out,'dir'), mkdir(out); end
seeds=101:130;
conditions={'all','original',19; 'random','original',0; ...
 'random','fixed-total',7; 'random','fixed-total',15; ...
 'balanced','fixed-total',15; 'random','original',15};
thresholds=[0.3413210142885347 0.5297546055663604 0.24277616843099487];
base=defaultEffectiveInputConfig(); base.recordStride=100;
metaPath=fullfile(out,'metadata.mat');
if ~isfile(metaPath)
    matlabVersion=version;
    save(metaPath,'base','seeds','conditions','thresholds','matlabVersion');
end
for seed=seeds
 for c=1:size(conditions,1)
    path=fullfile(out,sprintf('seed-%03d-condition-%d.mat',seed,c));
    if isfile(path)
        previous=load(path,'complete');
        assert(isfield(previous,'complete') && previous.complete);
        fprintf('Skip seed %d condition %d\n',seed,c); continue;
    end
    cfg=base; cfg.seed=seed; cfg.selectionPolicy=conditions{c,1};
    cfg.strengthMode=conditions{c,2}; cfg.neighborCount=conditions{c,3};
    fprintf('Start seed %d condition %d (%s %s k=%d)\n',seed,c,cfg.selectionPolicy,cfg.strengthMode,cfg.neighborCount);
    result=simulateEffectiveSocialInput(cfg);
    ts=result.timeSeries; ix=(cfg.burnIn+1):cfg.NSteps;
    m=[ts.globalOrder(ix) ts.localOrder(ix) ts.meanPairDistanceNorm(ix)];
    [blocks100,summary100]=blockTimeMetrics(m,cfg.dt,100,thresholds);
    [blocks300,summary300]=blockTimeMetrics(m,cfg.dt,300,thresholds);
    assert(max(abs(summary300.weightedMeans-[result.summary.globalOrder result.summary.localOrder result.summary.meanPairDistanceNorm]))<1e-10);
    complete=true;
    tmp=[path '.partial'];
    save(tmp,'result','blocks100','blocks300','summary100','summary300','thresholds','complete','-v7');
    movefile(tmp,path);
    fprintf('Saved seed %d condition %d\n',seed,c);
 end
end
disp('COMPLETE: 180 simulations saved with full time series');
end
