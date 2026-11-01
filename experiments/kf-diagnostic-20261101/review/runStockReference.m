function runStockReference()
% Independent full-length reference: unmodified repository simulator only.
here=fileparts(mfilename('fullpath'));
root=fileparts(fileparts(fileparts(here)));
restoredefaultpath;
addpath(fullfile(root,'effective-social-input'));
cfg=defaultEffectiveInputConfig();
cfg.seed=2026110101;
cfg.NSteps=2200; cfg.burnIn=1000; cfg.recordStride=1;
cfg.scheduledOcclusion=true; cfg.selectionPolicy='fixed-id';
cfg.neighborCount=15; cfg.occlusionHiddenCount=7;
cfg.occlusionStartStep=1401; cfg.occlusionDurationSteps=120;
cfg.predictionMode='kalman';
path=fullfile(here,'stock_full_reference.mat');
assert(~isfile(path),'Refuse to overwrite stock reference.');
started=tic;
result=simulateEffectiveSocialInput(cfg);
seconds=toc(started);
save(path,'result','cfg','seconds','-v7');
occ=1401:1520;
expected=15*ones(120,1); expected(101:120)=8;
assert(isequal(result.timeSeries.meanActiveNeighbors(occ),expected));
assert(all(result.timeSeries.meanVisibleNeighbors(occ)==8));
record=struct('status','PASS','source','unmodified repository simulator', ...
    'seed',cfg.seed,'seconds',seconds,'N',20,'k',15,'hidden',7, ...
    'D_steps',120,'age_units','simulation time, dt=0.3', ...
    'firstPTruncationMissingStep',101,'truncatedHiddenEdgeSteps',2800, ...
    'eligibleHiddenEdgeSteps',16800,'truncationFraction',20/120, ...
    'GO_occ',mean(result.timeSeries.globalOrder(occ)), ...
    'pair_occ',mean(result.timeSeries.meanPairDistanceNorm(occ)), ...
    'matlabVersion',version);
jsonPath=fullfile(here,'stock_full_reference.json');
assert(~isfile(jsonPath),'Refuse to overwrite stock reference metadata.');
fid=fopen(jsonPath,'w');
assert(fid>=0); cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(record,PrettyPrint=true));
fprintf('STOCK_FULL_REFERENCE_PASS seconds=%.3f P_cutoff=101 fraction=%.9f\n',seconds,20/120);
end
