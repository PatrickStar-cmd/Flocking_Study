function rows = runPersistenceEndpointSensitivity(seeds)
%RUNPERSISTENCEENDPOINTSENSITIVITY Compare persistence endpoint definitions.
if nargin < 1, seeds = 41:50; end
root = fileparts(mfilename('fullpath'));
outDir = fullfile(root, 'results', 'persistence-endpoint-sensitivity');
if ~exist(outDir, 'dir'), mkdir(outDir); end
csvPath = fullfile(outDir, 'endpoint-sensitivity.csv');
thresholds = [0.3413210142885347, 0.5297546055663604, 0.24277616843099487];
conditions = {'all','original',19; 'random','fixed-total',7; ...
    'random','fixed-total',15; 'balanced','fixed-total',15};
rows = struct('seed',{},'policy',{},'strengthMode',{},'k',{}, ...
    'pointFraction',{},'block100_90',{},'block100_95',{}, ...
    'block300_90',{},'block300_95',{},'longestFailureSteps',{}, ...
    'summaryGO',{},'summaryLO',{},'summaryPairDistance',{});
for seed = seeds
    for c = 1:size(conditions,1)
        cfg = defaultEffectiveInputConfig();
        cfg.NAgents=20; cfg.NNeurons=100; cfg.NSteps=8000;
        cfg.burnIn=6000; cfg.recordStride=1; cfg.seed=seed;
        cfg.selectionPolicy=conditions{c,1}; cfg.strengthMode=conditions{c,2};
        cfg.neighborCount=conditions{c,3};
        result=simulateEffectiveSocialInput(cfg); ts=result.timeSeries;
        ix=(cfg.burnIn+1):cfg.NSteps;
        ok=ts.globalOrder(ix)>=thresholds(1) & ts.localOrder(ix)>=thresholds(2) & ...
            ts.meanPairDistanceNorm(ix)<=thresholds(3);
        row.seed=seed; row.policy=string(cfg.selectionPolicy); row.strengthMode=string(cfg.strengthMode); row.k=cfg.neighborCount;
        row.pointFraction=mean(ok); row.block100_90=blockEndpoint(ok,100,.90); row.block100_95=blockEndpoint(ok,100,.95);
        row.block300_90=blockEndpoint(ok,300,.90); row.block300_95=blockEndpoint(ok,300,.95);
        row.longestFailureSteps=longestFailure(ok); row.summaryGO=result.summary.globalOrder;
        row.summaryLO=result.summary.localOrder; row.summaryPairDistance=result.summary.meanPairDistanceNorm;
        rows(end+1)=row; %#ok<AGROW>
    end
end
rows=struct2table(rows); writetable(rows,csvPath); save(fullfile(outDir,'endpoint-sensitivity.mat'),'rows','thresholds','conditions');
metadata.seeds=seeds; metadata.thresholds=thresholds; metadata.conditions=conditions;
metadata.note='Endpoint sensitivity only; thresholds frozen from independent-validation.json.';
save(fullfile(outDir,'metadata.mat'),'metadata');
end

function pass=blockEndpoint(ok,blockSize,target)
nBlocks=floor(numel(ok)/blockSize); if nBlocks==0, pass=false; return; end
blockRates=zeros(nBlocks,1);
for i=1:nBlocks, ix=(i-1)*blockSize+(1:blockSize); blockRates(i)=mean(ok(ix)); end
pass=mean(blockRates>=target)>=target;
end

function longest=longestFailure(ok)
f=~ok(:).'; e=diff([false f false]); s=find(e==1); t=find(e==-1)-1;
if isempty(s), longest=0; else, longest=max(t-s+1); end
end
