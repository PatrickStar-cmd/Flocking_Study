function rows = runNoisyBearingKFPilot(outDir)
% Frozen exploratory paired experiment. Existing output directories are refused.
assert(nargin==1 && ~isfolder(outDir),'Output directory already exists or is missing');
mkdir(outDir); diary(fullfile(outDir,'run.log')); cleanup=onCleanup(@() diary('off'));
cfg=defaultEffectiveInputConfig(); cfg.neighborCount=15;
cfg.selectionPolicy='random'; cfg.strengthMode='original';
cfg.scheduledOcclusion=false; cfg.turnExperiment=false; cfg.openArena=false;
cfg.recordStride=1; cfg.storeTimeSeries=true; cfg.storeTrajectory=true;
cfg.pilotSampleStride=100; cfg.bearingProcessNoise=1e-4;
seeds=5101:5106; noiseStd=[.01,.10]; modes={'raw','kf-posterior'};
metadata=struct('schema','noisy-bearing-kf-pilot-v1','date','2026-10-14', ...
    'evidence','exploratory','matlabVersion',version,'baseConfig',cfg, ...
    'seeds',seeds,'noiseStd',noiseStd,'modes',{modes},'expectedRuns',24, ...
    'primary','mean actual GO in steps 6001:8000', ...
    'KFInput','same-step posterior mean after CV prediction and noisy update', ...
    'pairing','same initialization, neighbor selection, and noise perturbations; closed-loop bearings may diverge', ...
    'timeUnit','model time','noiseStreamSeedOffset',15000);
sourceNames={'simulateEffectiveSocialInput.m','defaultEffectiveInputConfig.m', ...
    'task3DefaultConfig.m','buildRingConnectivity.m','task3BatchTrack.m', ...
    'simulateNoisyBearingKFPilot.m','noisyBearingKFPilotInputStep.m', ...
    'noisyBearingKFControllerStep.m','testNoisyBearingKFPilot.m','runNoisyBearingKFPilot.m'};
sourceDir=fullfile(outDir,'source-snapshot'); mkdir(sourceDir);
hashes=struct();
for i=1:numel(sourceNames)
    path=which(sourceNames{i}); assert(~isempty(path));
    copyfile(path,fullfile(sourceDir,sourceNames{i}));
    key=matlab.lang.makeValidName(sourceNames{i}); hashes.(key)=sha256(path);
end
metadata.sourceHashes=hashes;
[~,metadata.branch]=system('git branch --show-current');
[~,metadata.sourceCommit]=system('git rev-parse HEAD');
writeJSON(fullfile(outDir,'metadata.json'),metadata);
audit=testNoisyBearingKFPilot(); writeJSON(fullfile(outDir,'software-test-audit.json'),audit);
rows=struct([]); runIndex=0;
for seed=seeds
 for sd=noiseStd
    pair=struct(); cfg.seed=seed; cfg.bearingMeasurementStd=sd;
    cfg.bearingObservationVariance=sd^2;
    for mi=1:2
        cfg.noisyBearingMode=modes{mi}; started=tic;
        result=simulateNoisyBearingKFPilot(cfg);
        s=result.summary;
        r=struct('seed',seed,'noiseStd',sd,'method',string(modes{mi}), ...
            'meanGO',s.globalOrder,'meanLO',s.localOrder, ...
            'meanPairDistanceNorm',s.meanPairDistanceNorm,'meanSpeed',s.meanSpeed, ...
            'observedBearingRMSE',sqrt(mean(result.timeSeries.observedBearingRMSE(6001:8000).^2)), ...
            'usedBearingRMSE',sqrt(mean(result.timeSeries.usedBearingRMSE(6001:8000).^2)), ...
            'inputTruthRMS',sqrt(mean(result.timeSeries.inputTruthRMS(6001:8000).^2)), ...
            'meanKFVariance',s.meanKFVariance,'elapsedSeconds',toc(started));
        runIndex=runIndex+1;
        if runIndex==1, rows=r; else, rows(runIndex)=r; end
        if mi==1, pair.raw=result; else, pair.kf=result; end
        fprintf('RUN %d/24 seed=%d std=%.2f mode=%s GO=%.6f LO=%.6f pair/L=%.6f angleRMSE=%.6f elapsed=%.2fs\n', ...
            runIndex,seed,sd,modes{mi},r.meanGO,r.meanLO,r.meanPairDistanceNorm,r.usedBearingRMSE,r.elapsedSeconds);
    end
    assert(isequal(pair.raw.initialPosition,pair.kf.initialPosition));
    assert(isequal(pair.raw.initialHeadings,pair.kf.initialHeadings));
    assert(isequal(pair.raw.initialNeuralStateUsed,pair.kf.initialNeuralStateUsed));
    assert(isequal(pair.raw.pilot.selectedHistory,pair.kf.pilot.selectedHistory));
    assert(isequal(pair.raw.pilot.finalNoiseStreamState,pair.kf.pilot.finalNoiseStreamState));
    for i=1:numel(pair.raw.pilot.samples)
        assert(isequal(pair.raw.pilot.samples{i}.standardNoise,pair.kf.pilot.samples{i}.standardNoise));
    end
    filename=fullfile(outDir,sprintf('seed-%d-std-%03d-pair.mat',seed,round(sd*1000)));
    assert(~isfile(filename)); save(filename,'pair','-v7');
 end
end
writetable(struct2table(rows),fullfile(outDir,'summary.csv'));
writeJSON(fullfile(outDir,'completion.json'),struct('status','COMPLETED','runs',runIndex,'pairs',12));
end

function writeJSON(path,data)
assert(~isfile(path)); fid=fopen(path,'w'); assert(fid>=0); c=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(data,PrettyPrint=true));
end

function h=sha256(path)
fid=fopen(path,'rb'); c=onCleanup(@() fclose(fid)); bytes=fread(fid,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
