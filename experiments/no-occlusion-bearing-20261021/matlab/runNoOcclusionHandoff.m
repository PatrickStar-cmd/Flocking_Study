function rows=runNoOcclusionHandoff(outDir,profile)
% Portable runner. Frozen dynamics/filters are unchanged; profiles describe
% reproduction or software checks, not a newly completed scientific study.
if nargin<2, profile='smoke'; end
here=fileparts(mfilename('fullpath')); root=fileparts(here);
addpath(here,'-begin');
assert(~isfolder(outDir) && ~isfile(outDir),'Refuse existing output');
cfg=jsondecode(fileread(fullfile(root,'configs','base-config.json')));
switch profile
 case 'smoke'
  seeds=6299; distributions={'gaussian','mixture'}; levels=[.01,.1];
  cfg.NSteps=80; cfg.burnIn=40; cfg.pilotSampleStride=20;
 case 'pilot'
  seeds=6201:6206; distributions={'gaussian','mixture'}; levels=[.01,.1];
 case 'validation64'
  seeds=73001:73064; distributions={'mixture'}; levels=.01;
 otherwise
  error('profile must be smoke, pilot, or validation64');
end
assert(~cfg.scheduledOcclusion && ~cfg.occlusionEnabled && ~cfg.turnExperiment && ~cfg.openArena);
mkdir(outDir); mkdir(fullfile(outDir,'runs'));
diary(fullfile(outDir,'run.log')); cleanup=onCleanup(@() diary('off'));
methods={'oracle','raw','selected-kf','imm'}; rows=struct([]); uniqueRuns=0;
metadata=struct('schema','no-occlusion-handoff-v1','profile',profile,'seeds',seeds, ...
 'baseConfig',cfg,'noiseStd',levels,'distributions',{distributions},'methods',{methods}, ...
 'oracleReuse','one exact-bearing trajectory per seed reused across noise conditions', ...
 'evidence','smoke is software only; pilot and validation64 reproduce previously used seeds', ...
 'frozenSourceManifest','provenance/frozen-source-manifest.json','matlabVersion',version);
writeJSON(fullfile(outDir,'metadata.json'),metadata);
for seed=seeds
 c=cfg; c.seed=seed; c.noisyBearingMode='raw'; c.bearingMeasurementStd=0;
 oracle=simulateOracleBearingDiagnostic(c); uniqueRuns=uniqueRuns+1;
 assert(isempty(oracle.pilot.estimator));
 assert(all(oracle.timeSeries.usedBearingRMSE==0));
 result=oracle; save(fullfile(outDir,'runs',sprintf('seed-%d-oracle.mat',seed)),'result','-v7');
 for di=1:numel(distributions)
  for sd=levels
   c=cfg; c.seed=seed; c.noiseDistribution=distributions{di};
   c.bearingMeasurementStd=sd; c.bearingObservationVariance=sd^2;
   c.nonGaussianOutlierProbability=.05*strcmp(c.noiseDistribution,'mixture');
   c.nonGaussianOutlierScale=8;
   for mi=1:4
    if mi==1
     result=oracle;
     matFile=sprintf('runs/seed-%d-oracle.mat',seed);
    else
     c.noisyBearingMode=methods{mi}; result=simulateBearingIMMClosedLoopPilot(c);
     uniqueRuns=uniqueRuns+1;
     matFile=sprintf('runs/%s-seed-%d-std-%03d-%s.mat',c.noiseDistribution,seed,round(sd*1000),methods{mi});
     save(fullfile(outDir,matFile),'result','-v7');
    end
    assert(isequal(oracle.initialPosition,result.initialPosition));
    assert(isequal(oracle.initialHeadings,result.initialHeadings));
    assert(isequal(oracle.initialNeuralStateUsed,result.initialNeuralStateUsed));
    assert(isequal(oracle.pilot.selectedHistory,result.pilot.selectedHistory));
    if mi==2
     rawNoise=result.pilot.evalNoise;
     rawNoiseState=result.pilot.finalNoiseStreamState;
     rawMixtureState=result.pilot.finalMixtureStreamState;
    elseif mi>2
     assert(isequal(rawNoise,result.pilot.evalNoise));
     assert(isequal(rawNoiseState,result.pilot.finalNoiseStreamState));
     assert(isequal(rawMixtureState,result.pilot.finalMixtureStreamState));
    end
    ev=cfg.burnIn+1:cfg.NSteps; s=result.summary;
    row=struct('seed',seed,'distribution',string(c.noiseDistribution),'measurementStd',sd, ...
     'method',string(methods{mi}),'meanGO',s.globalOrder,'meanLO',s.localOrder, ...
     'meanPairDistanceNorm',s.meanPairDistanceNorm,'meanSpeed',s.meanSpeed, ...
     'usedBearingRMSE',sqrt(mean(result.timeSeries.usedBearingRMSE(ev).^2)), ...
     'inputTruthRMS',sqrt(mean(result.timeSeries.inputTruthRMS(ev).^2)), ...
     'matFile',string(matFile));
    if isempty(rows), rows=row; else, rows(end+1)=row; end
   end
   writetable(struct2table(rows),fullfile(outDir,'progress-summary.csv'));
  end
 end
 clear oracle result rawNoise;
end
writetable(struct2table(rows),fullfile(outDir,'summary.csv'));
expectedUnique=numel(seeds)*(1+3*numel(levels)*numel(distributions));
assert(uniqueRuns==expectedUnique);
assert(numel(rows)==numel(seeds)*4*numel(levels)*numel(distributions));
writeJSON(fullfile(outDir,'completion.json'),struct('status','COMPLETED','uniqueRuns',uniqueRuns, ...
 'summaryRows',numel(rows),'pairedSeeds',numel(seeds),'noOcclusion',true));
fprintf('COMPLETED %s: %d unique trajectories, %d condition rows\n',profile,uniqueRuns,numel(rows));
end

function writeJSON(path,data)
assert(~isfile(path)); fid=fopen(path,'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(data,PrettyPrint=true));
end
