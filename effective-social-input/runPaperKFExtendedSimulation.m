function runPaperKFExtendedSimulation(profile)
%RUNPAPERKFEXTENDEDSIMULATION Paired KF compensation sweep under occlusion.
%
% Profiles:
%   smoke : 2 new seeds, 2 occlusion loads, 2 horizons
%   pilot : 6 new seeds, 3 occlusion loads, 5 horizons
%   full  : 30 new seeds, 3 occlusion loads, 6 horizons
%
% The full profile is intentionally stored under a new result directory. It
% does not overwrite the original paper-kf-simulation results.

if nargin < 1, profile = "smoke"; end
profile = string(profile);
[seeds,hiddenCounts,durationSteps,maxWorkers] = profileGrid(profile);

root = fileparts(mfilename('fullpath'));
out = fullfile(root,'results','paper-kf-extended-v3',char(profile));
ckptDir = fullfile(out,'checkpoints');
seedDir = fullfile(out,'seed-results');
if ~exist(out,'dir'), mkdir(out); end
if ~exist(ckptDir,'dir'), mkdir(ckptDir); end
if ~exist(seedDir,'dir'), mkdir(seedDir); end

modes = {'none','hold','mean','kalman','kalman-gm'};
warmupSteps = 200;
postSteps = 100;
recoverySteps = 300;
branchSteps = warmupSteps + max(durationSteps) + recoverySteps;

base = defaultEffectiveInputConfig();
base.NAgents = 20; base.NNeurons = 100;
base.NSteps = 5000; base.burnIn = 3000;
base.recordStride = 1; base.scheduledOcclusion = true;
base.selectionPolicy = 'fixed-id';
base.neighborCount = 15;
base.occlusionStartStep = warmupSteps + 1;
base.occlusionDurationSteps = 0;
base.predictionMode = 'none';
base.distanceWeightMode = 'equal';
base.strengthMode = 'original';
base.bearingMeasurementStd = 0.01;

pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes',min(maxWorkers,feature('numcores')));
end

seedTables = cell(numel(seeds),1);
seedCurves = cell(numel(seeds),1);

fprintf('profile=%s seeds=%d hidden=%s T=%s workers=%d\n', ...
    profile,numel(seeds),mat2str(hiddenCounts),mat2str(durationSteps*base.dt),maxWorkers);

parfor si = 1:numel(seeds)
    seed = seeds(si);
    cfg = base;
    cfg.seed = seed;
    ckptPath = fullfile(ckptDir,sprintf('checkpoint-%d.mat',seed));
    checkpoint = loadOrCreateCheckpoint(ckptPath,cfg);

    branch = cfg;
    branch.NSteps = branchSteps;
    branch.burnIn = 0;
    branch.initialXY = checkpoint.finalPosition;
    branch.initialHeadings = checkpoint.finalHeadings;
    branch.initialNeuralState = checkpoint.finalNeuralState;
    branch.initialPreferredDirections = checkpoint.finalPreferredDirections;
    branch.occlusionStartStep = warmupSteps + 1;
    branch.occlusionDurationSteps = 0;
    branch.predictionMode = 'none';
    control = simulateEffectiveSocialInput(branch);
    checkpointGO = mean(control.timeSeries.globalOrder(1:warmupSteps));

    nBranches = numel(hiddenCounts)*numel(durationSteps)*numel(modes);
    rows = cell(nBranches,25);
    curveSums = zeros(numel(hiddenCounts),numel(durationSteps),numel(modes), ...
        branchSteps,2);
    idx = 0;
    for hi = 1:numel(hiddenCounts)
        hiddenCount = hiddenCounts(hi);
        for ti = 1:numel(durationSteps)
            duration = durationSteps(ti);
            eventIx = (warmupSteps+1):(warmupSteps+duration);
            postIx = (warmupSteps+duration+1):(warmupSteps+duration+postSteps);
            for mi = 1:numel(modes)
                branch.occlusionHiddenCount = hiddenCount;
                branch.occlusionDurationSteps = duration;
                branch.predictionMode = modes{mi};
                result = simulateEffectiveSocialInput(branch);
                assertCommonPrefix(control,result,warmupSteps);

                [predErr,predErr90,predCount] = hiddenPredictionError(result);
                rec = recoveryTime(result.timeSeries.globalOrder, ...
                    control.timeSeries.globalOrder,postIx,branch.dt);
                eventInputMass = mean(result.timeSeries.meanInputMass(eventIx));
                eventActive = mean(result.timeSeries.meanActiveNeighbors(eventIx));
                eventFresh = mean(result.timeSeries.meanVisibleNeighbors(eventIx));
                eventConcentration = mean(result.timeSeries.meanInputConcentration(eventIx));
                eventGOLoss = mean(control.timeSeries.globalOrder(eventIx) - ...
                    result.timeSeries.globalOrder(eventIx));
                eventGODeviation = mean(abs(control.timeSeries.globalOrder(eventIx) - ...
                    result.timeSeries.globalOrder(eventIx)));
                postGOLoss = mean(control.timeSeries.globalOrder(postIx) - ...
                    result.timeSeries.globalOrder(postIx));
                postGODeviation = mean(abs(control.timeSeries.globalOrder(postIx) - ...
                    result.timeSeries.globalOrder(postIx)));
                eventPairIncrease = mean(result.timeSeries.meanPairDistanceNorm(eventIx) - ...
                    control.timeSeries.meanPairDistanceNorm(eventIx));
                eventPairDeviation = mean(abs(result.timeSeries.meanPairDistanceNorm(eventIx) - ...
                    control.timeSeries.meanPairDistanceNorm(eventIx)));
                postPairIncrease = mean(result.timeSeries.meanPairDistanceNorm(postIx) - ...
                    control.timeSeries.meanPairDistanceNorm(postIx));
                postPairDeviation = mean(abs(result.timeSeries.meanPairDistanceNorm(postIx) - ...
                    control.timeSeries.meanPairDistanceNorm(postIx)));
                eventGOPreservation = mean(result.timeSeries.globalOrder(eventIx)) / ...
                    mean(control.timeSeries.globalOrder(eventIx));
                eventMinGO = min(result.timeSeries.globalOrder(eventIx));

                idx = idx + 1;
                rows(idx,:) = {seed,hiddenCount,hiddenCount/branch.neighborCount, ...
                    duration*branch.dt,duration,modes{mi},checkpointGO,eventGOLoss, ...
                    eventGODeviation,postGOLoss,postGODeviation,eventGOPreservation, ...
                    eventMinGO,eventPairIncrease,eventPairDeviation,postPairIncrease, ...
                    postPairDeviation,eventInputMass,eventActive,eventFresh, ...
                    eventConcentration,predErr,predErr90,predCount,rec};

                curveSums(hi,ti,mi,:,1) = ...
                    reshape(squeeze(curveSums(hi,ti,mi,:,1)),1,[]) + ...
                    (control.timeSeries.globalOrder - result.timeSeries.globalOrder).';
                curveSums(hi,ti,mi,:,2) = ...
                    reshape(squeeze(curveSums(hi,ti,mi,:,2)),1,[]) + ...
                    (result.timeSeries.meanPairDistanceNorm - ...
                    control.timeSeries.meanPairDistanceNorm).';
            end
        end
    end

    names = {'seed','hiddenCount','hiddenFraction','T','durationSteps','mode', ...
        'checkpointGO','eventGOLoss','eventGODeviation','postGOLoss','postGODeviation', ...
        'eventGOPreservation','eventMinGO','eventPairDistanceIncrease', ...
        'eventPairDistanceDeviation','postPairDistanceIncrease', ...
        'postPairDistanceDeviation','eventInputMass','eventActiveNeighbors', ...
        'eventFreshNeighbors','eventInputConcentration', ...
        'meanPredictionError','p90PredictionError','predictionErrorCount','recoveryTime'};
    seedTable = cell2table(rows,'VariableNames',names);
    curveMeta = struct('durationSteps',durationSteps,'hiddenCounts',hiddenCounts, ...
        'modes',{modes},'branchSteps',branchSteps,'warmupSteps',warmupSteps, ...
        'postSteps',postSteps,'recoverySteps',recoverySteps);
    writeSeedResults(seedDir,seed,seedTable,curveSums,curveMeta);
    seedTables{si} = seedTable;
    seedCurves{si} = curveSums;
    fprintf('seed %d complete\n',seed);
end

R = vertcat(seedTables{:});
R = sortrows(R,{'hiddenCount','T','mode','seed'});
writetable(R,fullfile(out,'summary.csv'));

curveSums = zeros(size(seedCurves{1}));
for si = 1:numel(seedCurves)
    curveSums = curveSums + seedCurves{si};
end
curveMean = curveSums/numel(seeds);
curveMeta = struct('durationSteps',durationSteps,'hiddenCounts',hiddenCounts, ...
    'modes',{modes},'branchSteps',branchSteps,'warmupSteps',warmupSteps, ...
    'postSteps',postSteps,'recoverySteps',recoverySteps,'seeds',seeds, ...
    'dt',base.dt,'controls',{'curveMean(:,:,:,:,1)=mean(controlGO-resultGO); ' ... 
    'curveMean(:,:,:,:,2)=mean(resultPair-controlPair)'});
save(fullfile(out,'curves.mat'),'curveMean','curveMeta','-v7');

meta = struct('profile',char(profile),'N',base.NAgents,'M',base.NNeurons, ...
    'dt',base.dt,'speedScale',base.speedScale,'sigma_s',base.receptiveFieldWidth, ...
    'seeds',seeds,'hiddenCounts',hiddenCounts,'neighborCount',base.neighborCount, ...
    'durationSteps',durationSteps,'T',durationSteps*base.dt,'modes',{modes}, ...
    'warmupSteps',warmupSteps,'postSteps',postSteps,'recoverySteps',recoverySteps, ...
    'measurementStd',base.bearingMeasurementStd,'bearingProcessNoise',base.bearingProcessNoise, ...
    'bearingGaussMarkovProcessNoise',base.bearingGaussMarkovProcessNoise, ...
    'bearingVelocityTau',base.bearingVelocityTau, ...
    'bearingObservationVariance',base.bearingObservationVariance, ...
    'predictionMaxAge',base.predictionMaxAge,'predictionMaxVariance',base.predictionMaxVariance, ...
    'note','Paired checkpoint-fork exploratory simulation with new seeds; no threshold tuning.');
fid = fopen(fullfile(out,'metadata.json'),'w');
fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true));
fclose(fid);

fprintf('wrote %d rows to %s\n',height(R),fullfile(out,'summary.csv'));
end

function [seeds,hiddenCounts,durationSteps,maxWorkers] = profileGrid(profile)
switch profile
    case "smoke"
        seeds = 201:202;
        hiddenCounts = [7 11];
        durationSteps = [10 34];
        maxWorkers = 2;
    case "pilot"
        seeds = 201:206;
        hiddenCounts = [7 11 13];
        durationSteps = [10 20 34 60 100];
        maxWorkers = 6;
    case "full"
        seeds = 201:230;
        hiddenCounts = [7 11 13];
        durationSteps = [10 20 34 60 100 200];
        maxWorkers = 8;
    otherwise
        error('profile must be smoke, pilot, or full.');
end
end

function checkpoint = loadOrCreateCheckpoint(path,cfg)
if exist(path,'file')
    loaded = load(path,'checkpoint');
    checkpoint = loaded.checkpoint;
else
    checkpoint = simulateEffectiveSocialInput(cfg);
    data = struct('checkpoint',checkpoint);
    save(path,'-struct','data','-v7');
end
end

function writeSeedResults(seedDir,seed,seedTable,curveSums,curveMeta)
writetable(seedTable,fullfile(seedDir,sprintf('seed-%d-summary.csv',seed)));
data = struct('curveSums',curveSums,'curveMeta',curveMeta);
save(fullfile(seedDir,sprintf('seed-%d-curves.mat',seed)),'-struct','data','-v7');
end

function assertCommonPrefix(control,result,n)
fields = fieldnames(control.timeSeries);
for i = 1:numel(fields)
    assert(isequaln(control.timeSeries.(fields{i})(1:n),result.timeSeries.(fields{i})(1:n)), ...
        'Branch diverged before the occlusion event.');
end
end

function [meanError,p90Error,count] = hiddenPredictionError(result)
errors = [];
records = result.observationRecords;
for t = 1:numel(records)
    record = records{t};
    if isempty(record), continue; end
    mask = record.active & ~record.fresh;
    if ~any(mask(:)), continue; end
    injected = record.injectedBearing(mask);
    truth = record.trueBearing(mask);
    valid = isfinite(injected) & isfinite(truth);
    if any(valid)
        delta = atan2(sin(injected(valid)-truth(valid)),cos(injected(valid)-truth(valid)));
        errors = [errors; abs(delta(:))]; %#ok<AGROW>
    end
end
if isempty(errors)
    meanError = NaN; p90Error = NaN; count = 0;
else
    meanError = mean(errors);
    p90Error = sort(errors);
    p90Error = p90Error(max(1,ceil(0.9*numel(p90Error))));
    count = numel(errors);
end
end

function rec = recoveryTime(go,controlGO,postIx,dt)
tol = 0.05;
ok = abs(go(postIx)-controlGO(postIx)) <= tol;
j = find(conv(double(ok),ones(10,1),'same') >= 10,1);
if isempty(j)
    rec = NaN;
else
    rec = (j-1)*dt;
end
end
