function result=runBearingIMMReplayPilot(outDir,stage,gaussianRoot,mixtureRoot)
% Two sequential stages. Validation cannot begin without frozen selection.
assert(any(strcmp(stage,{'development','validation'})));
if strcmp(stage,'development')
 assert(~isfolder(outDir),'Existing output directory is protected'); mkdir(outDir);
else
 assert(isfile(fullfile(outDir,'frozen-selection.json')),'Development must finish first');
 assert(~isfolder(fullfile(outDir,'validation')),'Validation directory already exists');
end
stageDir=fullfile(outDir,stage); mkdir(stageDir);
diary(fullfile(stageDir,'run.log')); diaryCleanup=onCleanup(@() diary('off'));
P=[.98,.02;.10,.90]; grid=[.01,.1,1,10]; low=1e-4;
if strcmp(stage,'development')
 snapshot=fullfile(outDir,'source-snapshot'); mkdir(snapshot);
 names={'defaultEffectiveInputConfig.m','task3DefaultConfig.m','buildRingConnectivity.m', ...
  'simulateEffectiveSocialInput.m','task3BatchTrack.m','simulateNoisyBearingKFPilot.m', ...
  'noisyBearingKFPilotInputStep.m','noisyBearingKFControllerStep.m', ...
  'drawNonGaussianBearingNoise.m','simulateNonGaussianBearingKFPilot.m', ...
  'nonGaussianBearingKFPilotInputStep.m','testNoisyBearingKFPilot.m', ...
  'bearingIMMTrackStep.m','testBearingIMMTrack.m','buildBearingReplaySequence.m', ...
  'replayBearingModels.m','runBearingIMMReplayPilot.m'};
 hashes=struct();
 for i=1:numel(names)
  path=which(names{i}); assert(~isempty(path)); copyfile(path,fullfile(snapshot,names{i}));
  hashes.(matlab.lang.makeValidName(names{i}))=sha256(path);
 end
 protected={'defaultEffectiveInputConfig.m','scheduledOcclusionInput.m', ...
  'simulateEffectiveSocialInput.m','task3BatchInputStep.m'}; guards=struct();
 for i=1:numel(protected), guards.(matlab.lang.makeValidName(protected{i}))=sha256(which(protected{i})); end
 [~,branch]=system('git branch --show-current'); [~,commit]=system('git rev-parse HEAD');
 meta=struct('date','2026-10-14','schema','bearing-imm-replay-pilot-v1', ...
  'matlabVersion',version,'branch',strtrim(branch),'sourceCommit',strtrim(commit), ...
  'developmentSeeds',5101:5106,'validationSeeds',6101:6106, ...
  'noiseStd',[.01,.10],'distributions',{{'gaussian','mixture'}}, ...
  'highQGrid',grid,'lowQ',low,'transitionPerModelStep',P, ...
  'sourceHashes',hashes,'preexistingTrackedFileHashes',guards, ...
  'gaussianRoot',gaussianRoot,'mixtureRoot',mixtureRoot, ...
  'reference','raw-observation-controlled trajectory; identical observations for all offline estimators', ...
  'primary','IMM minus development-selected best single-Q KF RMSE in four conditions; seed paired t CI, Holm four comparisons', ...
  'evidence','development selection followed by new-seed small-sample test, not a closed-loop experiment');
 writeJSON(fullfile(outDir,'metadata.json'),meta);
 audit=testBearingIMMTrack(); writeJSON(fullfile(outDir,'software-test-audit.json'),audit);
 seeds=5101:5106; methods=method('raw','raw',0);
 methods(end+1)=method('original-KF','kf',low);
 for q=grid, methods(end+1)=method(sprintf('single-q-%g',q),'kf',q); end
 for q=grid, methods(end+1)=method(sprintf('IMM-q-%g',q),'imm',q); end
else
 frozen=jsondecode(fileread(fullfile(outDir,'frozen-selection.json')));
 meta=jsondecode(fileread(fullfile(outDir,'metadata.json')));
 names=fieldnames(meta.sourceHashes);
 for i=1:numel(names)
  path=fullfile(outDir,'source-snapshot',strrep(names{i},'_m','.m'));
  assert(strcmp(sha256(path),meta.sourceHashes.(names{i})),'Snapshot changed after development');
 end
 seeds=6101:6106; P=frozen.transitionPerModelStep;
 methods=[method('raw','raw',0),method('original-KF','kf',low), ...
  method('same-high-Q-KF','kf',frozen.immHighQ), ...
  method('selected-single-KF','kf',frozen.bestSingleQ), ...
  method('IMM','imm',frozen.immHighQ)];
end
allRows=struct([]); sensorAudits=cell(1,24); sourceHashes=struct(); idx=0;
distributions={'gaussian','mixture'};
for di=1:2
 distribution=distributions{di};
 for sd=[.01,.10]
  for seed=seeds
   idx=idx+1; started=tic;
   tag=sprintf('%s-seed-%d-std-%03d',distribution,seed,round(sd*1000));
   if strcmp(stage,'development')
    if di==1, root=gaussianRoot; else, root=mixtureRoot; end
    source=fullfile(root,sprintf('seed-%d-std-%03d-pair.mat',seed,round(sd*1000)));
    loaded=load(source,'pair'); r=loaded.pair.raw;
    sourceHashes.(matlab.lang.makeValidName(tag))=sha256(source);
   else
    cfg=frozen.baseConfig; cfg.seed=seed; cfg.bearingMeasurementStd=sd;
    cfg.bearingObservationVariance=sd^2; cfg.noisyBearingMode='raw';
    cfg.nonGaussianOutlierProbability=.05; cfg.nonGaussianOutlierScale=8;
    if di==1, r=simulateNoisyBearingKFPilot(cfg); else, r=simulateNonGaussianBearingKFPilot(cfg); end
    referenceFile=fullfile(stageDir,['reference-',tag,'.mat']);
    save(referenceFile,'r','-v7');
    sourceHashes.(matlab.lang.makeValidName(tag))=sha256(referenceFile);
   end
   [seq,sensorAudits{idx}]=buildBearingReplaySequence(r,distribution);
   baseConfig=r.config; clear r loaded;
   seqFile=fullfile(stageDir,['sequence-',tag,'.mat']); save(seqFile,'seq','-v7');
   [rows,o]=replayBearingModels(seq,methods,P);
   if idx==1, allRows=rows; else, allRows=[allRows,rows]; end %#ok<AGROW>
   if strcmp(stage,'validation')
    replayFile=fullfile(stageDir,['estimates-',tag,'.mat']); save(replayFile,'o','-v7');
   end
   fprintf('%s %d/24 %s elapsed %.1fs',stage,idx,tag,toc(started));
   for m=1:numel(rows), fprintf(' %s=%.6f',rows(m).method,rows(m).RMSE); end
   fprintf('\n');
   writetable(struct2table(allRows),fullfile(stageDir,'progress-summary.csv'));
   clear seq o;
  end
 end
end
writetable(struct2table(allRows),fullfile(stageDir,'summary.csv'));
assert(numel(unique(string({allRows.method})))==numel(methods));
writeJSON(fullfile(stageDir,'sensor-audit.json'),struct('status','PASS', ...
 'sequences',24,'audits',{sensorAudits},'referenceHashes',sourceHashes));
if strcmp(stage,'development')
 immScores=zeros(size(grid)); singleGrid=[low,grid]; singleScores=zeros(size(singleGrid));
 for i=1:numel(grid), immScores(i)=score(allRows,'imm',grid(i)); end
 for i=1:numel(singleGrid), singleScores(i)=score(allRows,'kf',singleGrid(i)); end
 [~,a]=min(immScores); [~,b]=min(singleScores);
 frozen=struct('status','FROZEN_BEFORE_VALIDATION','immHighQ',grid(a), ...
  'bestSingleQ',singleGrid(b),'lowQ',low,'transitionPerModelStep',P, ...
  'immHighQGrid',grid,'immScores',immScores,'singleQGrid',singleGrid, ...
  'singleScores',singleScores,'criterion','mean log(RMSE/raw RMSE), equal weight to each seed-condition sequence', ...
  'baseConfig',baseConfig,'developmentSequences',24, ...
  'developmentSeeds',5101:5106,'validationSeeds',6101:6106);
 writeJSON(fullfile(outDir,'frozen-selection.json'),frozen);
 fprintf('FROZEN IMM high q=%g; best single q=%g\n',frozen.immHighQ,frozen.bestSingleQ);
end
keys=fieldnames(meta.preexistingTrackedFileHashes);
for i=1:numel(keys)
 name=strrep(keys{i},'_m','.m');
 assert(strcmp(sha256(which(name)),meta.preexistingTrackedFileHashes.(keys{i})), ...
  'Preexisting tracked file changed during this stage');
end
result=struct('status','COMPLETED','stage',stage,'sequences',idx, ...
 'methods',numel(methods),'summaryRows',numel(allRows),'preexistingTrackedFilesUnchanged',true);
writeJSON(fullfile(stageDir,'completion.json'),result);
end

function m=method(name,type,q)
m=struct('name',name,'type',type,'q',q);
end

function s=score(rows,type,q)
if strcmp(type,'imm'), prefix='IMM-q-'; else, prefix='single-q-'; end
if strcmp(type,'kf') && q==1e-4, name='original-KF'; else, name=sprintf('%s%g',prefix,q); end
selected=rows(strcmp(string({rows.method}),name)); raw=rows(strcmp(string({rows.method}),'raw'));
assert(numel(selected)==24 && numel(raw)==24);
assert(isequal([selected.seed],[raw.seed]) && isequal([selected.measurementStd],[raw.measurementStd]));
s=mean(log([selected.RMSE]./[raw.RMSE]));
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
