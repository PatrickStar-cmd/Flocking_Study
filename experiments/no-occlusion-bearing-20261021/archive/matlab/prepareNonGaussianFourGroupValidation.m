function prepareNonGaussianFourGroupValidation(outDir,upstream)
% Freeze the new-seed study before workers produce any scientific outcome.
assert(~isfolder(outDir),'Existing results are protected'); mkdir(outDir);
diary(fullfile(outDir,'prepare.log')); c=onCleanup(@() diary('off'));
old=jsondecode(fileread(fullfile(upstream,'metadata.json')));
source=fullfile(upstream,'source-snapshot'); snapshot=fullfile(outDir,'source-snapshot'); mkdir(snapshot);
names=fieldnames(old.sourceHashes); hashes=struct();
for i=1:numel(names)
 name=strrep(names{i},'_m','.m'); file=fullfile(source,name);
 assert(strcmp(sha256(file),old.sourceHashes.(names{i})));
 copyfile(file,fullfile(snapshot,name)); hashes.(names{i})=sha256(file);
end
extra={'simulateOracleBearingDiagnostic.m','oracleBearingDiagnosticInputStep.m', ...
 'prepareNonGaussianFourGroupValidation.m','runNonGaussianFourGroupValidation.m'};
repoCode=fileparts(mfilename('fullpath'));
for i=1:numel(extra)
 file=fullfile(repoCode,extra{i}); copyfile(file,fullfile(snapshot,extra{i}));
 hashes.(matlab.lang.makeValidName(extra{i}))=sha256(file);
end
scripts={'analyzeNonGaussianFourGroupValidation.py','analyze_oracle_bearing_diagnostic_20261021.py', ...
 'recheck_bearing_imm_replay_pilot.py'};
scriptDir=fullfile(outDir,'analysis-scripts'); mkdir(scriptDir);
scriptHashes=struct();
for i=1:numel(scripts)
 file=fullfile(fileparts(repoCode),'handoff',scripts{i});
 copyfile(file,fullfile(scriptDir,scripts{i})); scriptHashes.(matlab.lang.makeValidName(scripts{i}))=sha256(file);
end
copyfile(fullfile(fileparts(repoCode),'handoff','non-gaussian-four-group-validation-20261021-plan.md'),fullfile(outDir,'pre-run-plan.md'));
copyfile(fullfile(upstream,'metadata.json'),fullfile(outDir,'upstream-metadata.json'));
copyfile(fullfile(upstream,'upstream-frozen-selection.json'),fullfile(outDir,'upstream-frozen-selection.json'));
protected={'defaultEffectiveInputConfig.m','scheduledOcclusionInput.m', ...
 'simulateEffectiveSocialInput.m','task3BatchInputStep.m'}; guards=struct();
for i=1:numel(protected), guards.(matlab.lang.makeValidName(protected{i}))=sha256(fullfile(repoCode,protected{i})); end
cfg=old.baseConfig; cfg.noiseDistribution='mixture'; cfg.nonGaussianOutlierProbability=.05;
cfg.nonGaussianOutlierScale=8; cfg.bearingMeasurementStd=.01; cfg.bearingObservationVariance=.0001;
cfg.noisyBearingMode='raw'; cfg.seed=73999;
assert(cfg.NAgents==20 && cfg.NSteps==8000 && cfg.burnIn==6000 && cfg.neighborCount==15);
assert(cfg.bearingProcessNoise==.1 && cfg.immLowQ==.0001 && cfg.immHighQ==.1);
[~,branch]=system('git branch --show-current'); [~,head]=system('git rev-parse HEAD');
meta=struct('schema','non-gaussian-four-group-validation-v1','date','2026-10-21', ...
 'branch',strtrim(branch),'sourceCommit',strtrim(head),'matlabVersion',version, ...
 'seeds',73001:73064,'nSeeds',64,'expectedRuns',256,'workerCount',2, ...
 'methods',{{'oracle','raw','selected-kf','imm'}},'baseConfig',cfg, ...
 'sourceHashes',hashes,'analysisScriptHashes',scriptHashes,'preexistingTrackedFileHashes',guards, ...
 'upstream',upstream,'evidence','frozen new-seed independent validation in a pilot-selected condition', ...
 'primaryEndpoint','mean GO at steps6001:8000', ...
 'contrasts',{{'oracle-raw','selected-kf-raw','imm-raw','imm-selected-kf','oracle-selected-kf','oracle-imm'}}, ...
 'familySize',6,'alphaFamily',.05,'practicalGOMargin',.05, ...
 'primaryInference','seed-paired t tests with Holm6; simultaneous 95% family Bonferroni t intervals; ordinary95 intervals also retained', ...
 'sensitivity','50000 paired seed bootstrap, Bonferroni percentile intervals; no sequential stopping or retuning', ...
 'bootstrapSeed',2026102173,'bootstrapRepeats',50000, ...
 'sampleSizeRationale','pilot worst paired SD0.303 inflated10percent; n64 approximately80percent power for GO difference0.15 at conservative alpha0.05/6; not powered to guarantee equivalence at0.05', ...
 'exampleSeeds',[73001,73016,73032,73064], ...
 'filterReplaySeeds',[73001,73064], ...
 'interpretation','oracle is a diagnostic reference, not a guaranteed GO upper bound; no added occlusion; means are not stability or safety');
writeJSON(fullfile(outDir,'metadata.json'),meta);
addpath(snapshot,'-begin');
audit=testBearingIMMClosedLoop();
c0=cfg; c0.NSteps=140; c0.burnIn=40; c0.pilotSampleStride=20; c0.bearingMeasurementStd=0;
a=simulateOracleBearingDiagnostic(c0); b=simulateBearingIMMClosedLoopPilot(c0);
assert(isequal(a.trajectory,b.trajectory)); assert(isempty(a.pilot.estimator));
assert(all(a.timeSeries.usedBearingRMSE==0) && all(a.timeSeries.inputTruthRMS==0));
% A full old trajectory establishes exact compatibility of the frozen program.
historical=load(fullfile(upstream,'mixture-seed-6201-std-010-triple.mat'),'triple');
r=simulateBearingIMMClosedLoopPilot(historical.triple.raw.config);
assert(isequal(r.trajectory,historical.triple.raw.trajectory));
assert(isequal(r.timeSeries.globalOrder,historical.triple.raw.timeSeries.globalOrder));
audit.oracleZeroNoiseExact=true; audit.oracleNoFilterState=true;
audit.fullArchived8000StepMixtureRawExact=true;
writeJSON(fullfile(outDir,'software-test-audit.json'),audit);
mkdir(fullfile(outDir,'worker-1')); mkdir(fullfile(outDir,'worker-2'));
writeJSON(fullfile(outDir,'ready.json'),struct('status','FROZEN_AND_TESTED_BEFORE_NEW_SEEDS','nSeeds',64,'runs',256));
fprintf('READY:64 independent new seeds,256 trajectories,2 disjoint workers.\n');
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
