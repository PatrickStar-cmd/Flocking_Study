function runNonGaussianFourGroupValidation(outDir,worker)
% Two independent single-thread MATLAB processes; disjoint seed/file sets.
assert(any(worker==[1,2]));
ready=jsondecode(fileread(fullfile(outDir,'ready.json'))); assert(ready.runs==256);
meta=jsondecode(fileread(fullfile(outDir,'metadata.json')));
snapshot=fullfile(outDir,'source-snapshot'); addpath(snapshot,'-begin');
keys=fieldnames(meta.sourceHashes);
for i=1:numel(keys)
 name=strrep(keys{i},'_m','.m'); assert(strcmp(sha256(fullfile(snapshot,name)),meta.sourceHashes.(keys{i})));
end
dest=fullfile(outDir,sprintf('worker-%d',worker));
assert(~isfile(fullfile(dest,'started.json')),'Worker results already exist; no silent resume');
writeJSON(fullfile(dest,'started.json'),struct('worker',worker,'startedUTC',char(datetime('now','TimeZone','UTC'))));
diary(fullfile(dest,'run.log')); cleanup=onCleanup(@() diary('off'));
seeds=meta.seeds(worker:2:end); methods={'oracle','raw','selected-kf','imm'};
rows=struct([]); runCount=0;
for seed=seeds(:).'
 proof=[]; seedRows=struct([]); fileHashes=struct();
 for mi=1:4
  method=methods{mi}; cfg=meta.baseConfig; cfg.seed=seed; started=tic;
  if mi==1
   cfg.noisyBearingMode='raw'; cfg.bearingMeasurementStd=0;
   result=simulateOracleBearingDiagnostic(cfg);
   assert(isempty(result.pilot.estimator));
   assert(all(result.timeSeries.usedBearingRMSE==0) && all(result.timeSeries.inputTruthRMS==0));
  else
   cfg.noisyBearingMode=method; result=simulateBearingIMMClosedLoopPilot(cfg);
  end
  current=pairProof(result);
  if mi==1, proof=current; else, assert(isequaln(proof,current),'Pairing contract failed'); end
  ev=cfg.burnIn+1:cfg.NSteps; s=result.summary;
  row=struct('seed',seed,'method',string(method),'meanGO',s.globalOrder,'meanLO',s.localOrder, ...
   'meanPairDistanceNorm',s.meanPairDistanceNorm,'meanSpeed',s.meanSpeed, ...
   'observedBearingRMSE',sqrt(mean(result.timeSeries.observedBearingRMSE(ev).^2)), ...
   'usedBearingRMSE',sqrt(mean(result.timeSeries.usedBearingRMSE(ev).^2)), ...
   'inputTruthRMS',sqrt(mean(result.timeSeries.inputTruthRMS(ev).^2)), ...
   'elapsedSeconds',toc(started));
  if isempty(seedRows), seedRows=row; else, seedRows(end+1)=row; end
  if isempty(rows), rows=row; else, rows(end+1)=row; end
  file=sprintf('seed-%d-%s.mat',seed,method); path=fullfile(dest,file);
  assert(~isfile(path)); save(path,'result','-v7');
  fileHashes.(matlab.lang.makeValidName(file))=sha256(path);
  runCount=runCount+1;
  fprintf('WORKER %d RUN %d/128 seed=%d method=%s elapsed=%.1fs\n',worker,runCount,seed,method,row.elapsedSeconds);
  clear result current;
 end
 csvPath=fullfile(dest,sprintf('seed-%d-summary.csv',seed)); assert(~isfile(csvPath));
 writetable(struct2table(seedRows),csvPath);
 writeJSON(fullfile(dest,sprintf('completed-seed-%d.json',seed)),struct('status','COMPLETED','seed',seed, ...
  'runs',4,'allInitializationSelectionNoiseAndLabelsPaired',true,'fileHashes',fileHashes));
 writetable(struct2table(rows),fullfile(dest,'progress-summary.csv'));
 clear proof;
end
writetable(struct2table(rows),fullfile(dest,'summary.csv'));
writeJSON(fullfile(dest,'completion.json'),struct('status','COMPLETED','worker',worker,'seeds',seeds,'runs',runCount));
end

function p=pairProof(r)
p=struct('initialPosition',r.initialPosition,'initialHeadings',r.initialHeadings, ...
 'initialNeuralState',r.initialNeuralStateUsed,'selectedHistory',r.pilot.selectedHistory, ...
 'finalNoiseStreamState',r.pilot.finalNoiseStreamState, ...
 'finalMixtureStreamState',r.pilot.finalMixtureStreamState,'evalNoise',r.pilot.evalNoise);
for i=1:numel(r.pilot.samples)
 s=r.pilot.samples{i}; p.samples{i}=struct('standardNoise',s.standardNoise,'outlierMask',s.outlierMask,'standardNormal',s.standardNormal);
end
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
