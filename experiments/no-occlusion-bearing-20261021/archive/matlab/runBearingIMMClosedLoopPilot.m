function rows=runBearingIMMClosedLoopPilot(outDir,frozenPath)
assert(~isfolder(outDir),'Existing results are protected'); mkdir(outDir);
diary(fullfile(outDir,'run.log')); cleanup=onCleanup(@() diary('off'));
frozen=jsondecode(fileread(frozenPath));
assert(strcmp(frozen.status,'FROZEN_BEFORE_VALIDATION'));
assert(frozen.immHighQ==.1 && frozen.bestSingleQ==.1 && frozen.lowQ==1e-4);
cfg=frozen.baseConfig; cfg.bearingProcessNoise=frozen.bestSingleQ;
cfg.immLowQ=frozen.lowQ; cfg.immHighQ=frozen.immHighQ;
cfg.immTransition=frozen.transitionPerModelStep;
seeds=6201:6206; noiseStd=[.01,.10]; distributions={'gaussian','mixture'};
modes={'raw','selected-kf','imm'}; labels={'raw','kf','imm'};
names={'defaultEffectiveInputConfig.m','task3DefaultConfig.m','buildRingConnectivity.m', ...
 'simulateEffectiveSocialInput.m','task3BatchTrack.m','noisyBearingKFControllerStep.m', ...
 'simulateNoisyBearingKFPilot.m','noisyBearingKFPilotInputStep.m', ...
 'simulateNonGaussianBearingKFPilot.m','nonGaussianBearingKFPilotInputStep.m', ...
 'drawNonGaussianBearingNoise.m','bearingIMMTrackStep.m', ...
 'bearingIMMClosedLoopControllerStep.m','bearingIMMClosedLoopInputStep.m', ...
 'simulateBearingIMMClosedLoopPilot.m','testBearingIMMClosedLoop.m', ...
 'runBearingIMMClosedLoopPilot.m'};
snapshot=fullfile(outDir,'source-snapshot'); mkdir(snapshot); hashes=struct();
for i=1:numel(names)
 path=which(names{i}); assert(~isempty(path)); copyfile(path,fullfile(snapshot,names{i}));
 hashes.(matlab.lang.makeValidName(names{i}))=sha256(path);
end
copyfile(frozenPath,fullfile(outDir,'upstream-frozen-selection.json'));
protected={'defaultEffectiveInputConfig.m','scheduledOcclusionInput.m', ...
 'simulateEffectiveSocialInput.m','task3BatchInputStep.m'}; guards=struct();
for i=1:numel(protected), guards.(matlab.lang.makeValidName(protected{i}))=sha256(which(protected{i})); end
[~,branch]=system('git branch --show-current'); [~,commit]=system('git rev-parse HEAD');
metadata=struct('schema','bearing-imm-closed-loop-pilot-v1','date','2026-10-14', ...
 'matlabVersion',version,'branch',strtrim(branch),'sourceCommit',strtrim(commit), ...
 'baseConfig',cfg,'seeds',seeds,'noiseStd',noiseStd,'distributions',{distributions}, ...
 'modes',{modes},'expectedRuns',72,'expectedTriples',24,'sourceHashes',hashes, ...
 'preexistingTrackedFileHashes',guards,'upstreamFrozenParameterSHA256',sha256(frozenPath), ...
 'primary','IMM minus selected single-Q KF mean actual GO, steps6001:8000; four seed-paired comparisons, Holm', ...
 'evidence','new-seed small-sample closed-loop test after upstream development selection', ...
 'pairing','same initialization, all-step selected edges, noise and labels; observations differ after geometry diverges', ...
 'timeUnit','model time','noiseStreamSeedOffset',15000,'mixtureStreamSeedOffset',25000);
writeJSON(fullfile(outDir,'metadata.json'),metadata);
audit=testBearingIMMClosedLoop(); writeJSON(fullfile(outDir,'software-test-audit.json'),audit);
rows=struct([]); run=0; tripleCount=0;
for di=1:2
 cfg.noiseDistribution=distributions{di};
 if di==1, cfg.nonGaussianOutlierProbability=0; else, cfg.nonGaussianOutlierProbability=.05; end
 cfg.nonGaussianOutlierScale=8;
 for sd=noiseStd
  cfg.bearingMeasurementStd=sd; cfg.bearingObservationVariance=sd^2;
  for seed=seeds
   cfg.seed=seed; triple=struct();
   for mi=1:3
    cfg.noisyBearingMode=modes{mi}; started=tic;
    r=simulateBearingIMMClosedLoopPilot(cfg); s=r.summary; eval=cfg.burnIn+1:cfg.NSteps;
    run=run+1;
    row=struct('seed',seed,'distribution',string(distributions{di}),'measurementStd',sd, ...
     'method',string(modes{mi}),'meanGO',s.globalOrder,'meanLO',s.localOrder, ...
     'meanPairDistanceNorm',s.meanPairDistanceNorm,'meanSpeed',s.meanSpeed, ...
     'observedBearingRMSE',sqrt(mean(r.timeSeries.observedBearingRMSE(eval).^2)), ...
     'usedBearingRMSE',sqrt(mean(r.timeSeries.usedBearingRMSE(eval).^2)), ...
     'inputTruthRMS',sqrt(mean(r.timeSeries.inputTruthRMS(eval).^2)), ...
     'elapsedSeconds',toc(started));
    if run==1, rows=row; else, rows(run)=row; end
    triple.(labels{mi})=r;
    fprintf('RUN %d/72 %s seed=%d std=%.2f %s GO=%.6f LO=%.6f dist/L=%.6f angleRMSE=%.6f elapsed=%.1fs\n', ...
     run,distributions{di},seed,sd,modes{mi},row.meanGO,row.meanLO,row.meanPairDistanceNorm,row.usedBearingRMSE,row.elapsedSeconds);
   end
   a=triple.raw;
   for mi=2:3
    b=triple.(labels{mi});
    assert(isequal(a.initialPosition,b.initialPosition));
    assert(isequal(a.initialHeadings,b.initialHeadings));
    assert(isequal(a.initialNeuralStateUsed,b.initialNeuralStateUsed));
    assert(isequal(a.pilot.selectedHistory,b.pilot.selectedHistory));
    assert(isequal(a.pilot.finalNoiseStreamState,b.pilot.finalNoiseStreamState));
    assert(isequal(a.pilot.finalMixtureStreamState,b.pilot.finalMixtureStreamState));
    assert(isequal(a.pilot.evalNoise,b.pilot.evalNoise));
    for i=1:numel(a.pilot.samples)
     assert(isequal(a.pilot.samples{i}.standardNoise,b.pilot.samples{i}.standardNoise));
     assert(isequal(a.pilot.samples{i}.outlierMask,b.pilot.samples{i}.outlierMask));
    end
   end
   file=fullfile(outDir,sprintf('%s-seed-%d-std-%03d-triple.mat',distributions{di},seed,round(sd*1000)));
   assert(~isfile(file)); save(file,'triple','-v7'); tripleCount=tripleCount+1;
   writetable(struct2table(rows),fullfile(outDir,'progress-summary.csv'));
   clear triple r a b;
  end
 end
end
writetable(struct2table(rows),fullfile(outDir,'summary.csv'));
keys=fieldnames(guards);
for i=1:numel(keys)
 assert(strcmp(sha256(which(strrep(keys{i},'_m','.m'))),guards.(keys{i})));
end
writeJSON(fullfile(outDir,'completion.json'),struct('status','COMPLETED','runs',run, ...
 'triples',tripleCount,'preexistingTrackedFilesUnchanged',true));
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
