function rows=runOracleBearingDiagnostic(outDir,upstream)
% Add six oracle controls to an archived 72-run paired experiment.
assert(~isfolder(outDir),'Existing results are protected'); mkdir(outDir);
diary(fullfile(outDir,'run.log')); cleanup=onCleanup(@() diary('off'));
meta=jsondecode(fileread(fullfile(upstream,'metadata.json')));
complete=jsondecode(fileread(fullfile(upstream,'completion.json')));
assert(strcmp(complete.status,'COMPLETED') && complete.runs==72);
archived=fullfile(upstream,'source-snapshot'); addpath(archived,'-begin');
keys=fieldnames(meta.sourceHashes); snapshot=fullfile(outDir,'source-snapshot'); mkdir(snapshot);
for i=1:numel(keys)
 name=strrep(keys{i},'_m','.m'); path=fullfile(archived,name);
 assert(strcmp(sha256(path),meta.sourceHashes.(keys{i}))); copyfile(path,fullfile(snapshot,name));
end
newNames={'simulateOracleBearingDiagnostic.m','oracleBearingDiagnosticInputStep.m', ...
 'runOracleBearingDiagnostic.m'}; newHashes=struct();
for i=1:numel(newNames)
 path=which(newNames{i}); copyfile(path,fullfile(snapshot,newNames{i}));
 newHashes.(matlab.lang.makeValidName(newNames{i}))=sha256(path);
end
copyfile(fullfile(upstream,'metadata.json'),fullfile(outDir,'upstream-metadata.json'));
copyfile(fullfile(upstream,'summary.csv'),fullfile(outDir,'upstream-summary.csv'));
protected={'defaultEffectiveInputConfig.m','scheduledOcclusionInput.m', ...
 'simulateEffectiveSocialInput.m','task3BatchInputStep.m'}; guards=struct();
repoCode=fileparts(mfilename('fullpath'));
for i=1:numel(protected), guards.(matlab.lang.makeValidName(protected{i}))=sha256(fullfile(repoCode,protected{i})); end
cfg=meta.baseConfig; cfg.noiseDistribution='gaussian'; cfg.nonGaussianOutlierProbability=0;
cfg.nonGaussianOutlierScale=8; cfg.bearingMeasurementStd=0;
cfg.bearingObservationVariance=.01; cfg.noisyBearingMode='raw';
% R remains positive for schema compatibility; the oracle never calls a filter.
[~,branch]=system('git branch --show-current'); [~,commit]=system('git rev-parse HEAD');
info=struct('schema','oracle-bearing-diagnostic-v1','date','2026-10-21', ...
 'branch',strtrim(branch),'sourceCommit',strtrim(commit),'matlabVersion',version, ...
 'upstream',upstream,'baseConfig',cfg,'seeds',6201:6206,'newUniqueTrajectories',6, ...
 'reusedUniqueTrajectories',72,'conditionTableRows',96,'newSourceHashes',newHashes, ...
 'preexistingTrackedFileHashes',guards,'evidence','exploratory diagnostic extension of previously examined seeds', ...
 'primary','Oracle minus raw mean GO per condition; four paired seed t intervals and Holm tests', ...
 'exampleSeedChosenBeforeRun',6201,'oracleReuse','one identical zero-noise trajectory per seed reused across four noise conditions', ...
 'interpretation','oracle is a diagnostic reference, not a mathematical upper bound on GO; no occlusion added');
writeJSON(fullfile(outDir,'metadata.json'),info);
% Short exact regression and truth-path separation.
c=cfg; c.NSteps=140; c.burnIn=40; c.pilotSampleStride=20; c.seed=6298;
a=simulateOracleBearingDiagnostic(c); b=simulateBearingIMMClosedLoopPilot(c);
assert(isequal(a.trajectory,b.trajectory)); assert(isempty(a.pilot.estimator));
assert(all(a.timeSeries.usedBearingRMSE==0) && all(a.timeSeries.inputTruthRMS==0));
c.noiseDistribution='mixture'; c.nonGaussianOutlierProbability=.05;
d=simulateOracleBearingDiagnostic(c); assert(isequal(a.trajectory,d.trajectory));
audit=testBearingIMMClosedLoop();
audit.oracleEqualsZeroNoiseLegacy=true; audit.oracleHasNoFilterState=true;
audit.oracleIndependentOfUnusedNoiseDistribution=true;
% Full old raw replay checks software/environment compatibility before reuse.
old=load(fullfile(upstream,'gaussian-seed-6201-std-100-triple.mat'),'triple');
reference=simulateBearingIMMClosedLoopPilot(old.triple.raw.config);
assert(isequal(reference.trajectory,old.triple.raw.trajectory));
assert(isequal(reference.timeSeries.globalOrder,old.triple.raw.timeSeries.globalOrder));
audit.full8000StepArchivedRawExact=true;
writeJSON(fullfile(outDir,'software-test-audit.json'),audit);
clear a b c d old reference;
rows=struct([]);
for seed=6201:6206
 cfg.seed=seed; started=tic; oracle=simulateOracleBearingDiagnostic(cfg);
 old=load(fullfile(upstream,sprintf('gaussian-seed-%d-std-010-triple.mat',seed)),'triple');
 for f={'initialPosition','initialHeadings','initialNeuralStateUsed'}
  assert(isequal(oracle.(f{1}),old.triple.raw.(f{1})));
 end
 assert(isequal(oracle.pilot.selectedHistory,old.triple.raw.pilot.selectedHistory));
 assert(isequal(oracle.pilot.finalNoiseStreamState,old.triple.raw.pilot.finalNoiseStreamState));
 assert(isequal(oracle.pilot.finalMixtureStreamState,old.triple.raw.pilot.finalMixtureStreamState));
 assert(isempty(oracle.pilot.estimator));
 assert(all(oracle.timeSeries.usedBearingRMSE==0) && all(oracle.timeSeries.inputTruthRMS==0));
 s=oracle.summary;
 row=struct('seed',seed,'method',"oracle",'meanGO',s.globalOrder,'meanLO',s.localOrder, ...
 'meanPairDistanceNorm',s.meanPairDistanceNorm,'meanSpeed',s.meanSpeed, ...
 'usedBearingRMSE',0,'inputTruthRMS',0,'elapsedSeconds',toc(started));
 if isempty(rows), rows=row; else, rows(end+1)=row; end
 path=fullfile(outDir,sprintf('oracle-seed-%d.mat',seed)); assert(~isfile(path)); save(path,'oracle','-v7');
 writetable(struct2table(rows),fullfile(outDir,'progress-oracle-summary.csv'));
 fprintf('ORACLE %d/6 seed=%d GO=%.6f LO=%.6f dist/L=%.6f elapsed=%.1fs\n', ...
  seed-6200,seed,row.meanGO,row.meanLO,row.meanPairDistanceNorm,row.elapsedSeconds);
 clear oracle old;
end
writetable(struct2table(rows),fullfile(outDir,'oracle-summary.csv'));
for i=1:numel(protected)
 assert(strcmp(guards.(matlab.lang.makeValidName(protected{i})),sha256(fullfile(repoCode,protected{i}))));
end
writeJSON(fullfile(outDir,'completion.json'),struct('status','COMPLETED','newRuns',6, ...
 'reusedRuns',72,'preexistingTrackedFilesUnchanged',true));
end

function writeJSON(path,data)
assert(~isfile(path)); fid=fopen(path,'w'); assert(fid>=0); c=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(data,PrettyPrint=true));
end

function h=sha256(path)
fid=fopen(path,'rb'); assert(fid>=0); c=onCleanup(@() fclose(fid)); bytes=fread(fid,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
