function rows = runBlockMetricPilot(seeds)
%RUNBLOCKMETRICPILOT Save fixed-time block means for endpoint design.
if nargin < 1, seeds = 51:55; end
root=fileparts(mfilename('fullpath')); outDir=fullfile(root,'results','block-metric-pilot');
if ~exist(outDir,'dir'), mkdir(outDir); end
thresholds=[0.3413210142885347,0.5297546055663604,0.24277616843099487];
conditions={'all','original',19;'random','fixed-total',7;'random','fixed-total',15;'balanced','fixed-total',15};
rows=struct('seed',{},'policy',{},'strengthMode',{},'k',{},'blockSize',{},'nBlocks',{},'meanGO',{},'meanLO',{},'meanPairDistance',{},'passBlocks',{},'passFraction',{},'durationWeightedPass',{});
for seed=seeds
 for c=1:size(conditions,1)
  cfg=defaultEffectiveInputConfig(); cfg.NAgents=20; cfg.NNeurons=100; cfg.NSteps=8000; cfg.burnIn=6000; cfg.recordStride=1; cfg.seed=seed; cfg.selectionPolicy=conditions{c,1}; cfg.strengthMode=conditions{c,2}; cfg.neighborCount=conditions{c,3};
  s=simulateEffectiveSocialInput(cfg); ix=(cfg.burnIn+1):cfg.NSteps;
  for bs=[100 300]
   edges=[1:bs:numel(ix),numel(ix)+1];
   for b=1:numel(edges)-1
    q=ix(edges(b):edges(b+1)-1); g=mean(s.timeSeries.globalOrder(q)); l=mean(s.timeSeries.localOrder(q)); d=mean(s.timeSeries.meanPairDistanceNorm(q));
    row.seed=seed; row.policy=string(cfg.selectionPolicy); row.strengthMode=string(cfg.strengthMode); row.k=cfg.neighborCount; row.blockSize=bs; row.nBlocks=1; row.meanGO=g; row.meanLO=l; row.meanPairDistance=d; row.passBlocks=double(g>=thresholds(1)&&l>=thresholds(2)&&d<=thresholds(3)); row.passFraction=1; row.durationWeightedPass=row.passBlocks; rows(end+1)=row; %#ok<AGROW>
   end
  end
 end
end
rows=struct2table(rows); writetable(rows,fullfile(outDir,'block-metrics.csv')); save(fullfile(outDir,'block-metrics.mat'),'rows','thresholds','conditions');
end
