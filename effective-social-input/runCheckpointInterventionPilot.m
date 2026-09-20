function runCheckpointInterventionPilot()
% Exploratory paired checkpoint intervention; no predictive model is used.
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','checkpoint-intervention-pilot-v1');
assert(~exist(out,'dir'),'Output directory exists; use a new name for reruns.');
mkdir(out);

seeds=136:137;
durationSteps=[10 34 100];
modes={'none','hold'};
warmupSteps=200;
recoverySteps=300;
branchSteps=warmupSteps+max(durationSteps)+recoverySteps;
eventStart=warmupSteps+1;
postSteps=100;

base=defaultEffectiveInputConfig();
base.NAgents=20;
base.NNeurons=100;
base.NSteps=6000;
base.burnIn=4000;
base.recordStride=1;
base.scheduledOcclusion=true;
base.selectionPolicy='fixed-id';
base.neighborCount=15;
base.occlusionHiddenCount=7;
base.occlusionStartStep=2;
base.occlusionDurationSteps=0;
base.predictionMode='none';
base.distanceWeightMode='equal';
base.strengthMode='original';

rows=cell(numel(seeds)*(1+numel(durationSteps)*numel(modes)),14);
idx=0;
for seed=seeds
    base.seed=seed;
    checkpoint=simulateEffectiveSocialInput(base);
    checkpointPath=fullfile(out,sprintf('checkpoint-%d.mat',seed));
    save(checkpointPath,'checkpoint','-v7');

    branch=base;
    branch.NSteps=branchSteps;
    branch.burnIn=0;
    branch.initialXY=checkpoint.finalPosition;
    branch.initialHeadings=checkpoint.finalHeadings;
    branch.initialNeuralState=checkpoint.finalNeuralState;
    branch.initialPreferredDirections=checkpoint.finalPreferredDirections;
    branch.occlusionStartStep=eventStart;
    branch.occlusionDurationSteps=0;
    control=simulateEffectiveSocialInput(branch);
    save(fullfile(out,sprintf('seed-%d-control.mat',seed)),'control','-v7');
    idx=idx+1;
    rows(idx,:)={seed,'control',0,0,mean(control.timeSeries.globalOrder),0,0, ...
        mean(control.timeSeries.meanPairDistanceNorm),0,0, ...
        mean(control.timeSeries.meanInputMass),0,NaN, ...
        mean(checkpoint.timeSeries.globalOrder(end-199:end))};

    for duration=durationSteps
        for mi=1:numel(modes)
            branch.predictionMode=modes{mi};
            branch.occlusionDurationSteps=duration;
            result=simulateEffectiveSocialInput(branch);
            assertCommonPrefix(control,result,warmupSteps);
            save(fullfile(out,sprintf('seed-%d-%s-T%g.mat',seed,modes{mi},duration*branch.dt)),'result','-v7');

            eventIx=eventStart:(eventStart+duration-1);
            postIx=(eventStart+duration):(eventStart+duration+postSteps-1);
            idx=idx+1;
            rows(idx,:)={seed,modes{mi},duration*branch.dt,duration, ...
                mean(result.timeSeries.globalOrder), ...
                mean(control.timeSeries.globalOrder(eventIx)-result.timeSeries.globalOrder(eventIx)), ...
                mean(control.timeSeries.globalOrder(postIx)-result.timeSeries.globalOrder(postIx)), ...
                mean(result.timeSeries.meanPairDistanceNorm), ...
                mean(result.timeSeries.meanPairDistanceNorm(eventIx)-control.timeSeries.meanPairDistanceNorm(eventIx)), ...
                mean(result.timeSeries.meanPairDistanceNorm(postIx)-control.timeSeries.meanPairDistanceNorm(postIx)), ...
                mean(result.timeSeries.meanInputMass), ...
                mean(result.timeSeries.meanInputMass(eventIx)-control.timeSeries.meanInputMass(eventIx)), ...
                mean(result.timeSeries.meanVisibleNeighbors(eventIx)), ...
                mean(checkpoint.timeSeries.globalOrder(end-199:end))};
        end
    end
    fprintf('Completed checkpoint intervention seed %d\n',seed);
end

resultTable=cell2table(rows,'VariableNames',{'seed','mode','T','durationSteps', ...
    'fullGO','eventGOLoss','postGOLoss','fullPairDistanceNorm', ...
    'eventPairDistanceIncrease','postPairDistanceIncrease','fullInputMass', ...
    'eventInputMassDifference','eventFreshNeighbors','checkpointGO'});
writetable(resultTable,fullfile(out,'results.csv'));

audit=struct('createdAt',char(datetime('now','TimeZone','Asia/Shanghai')), ...
    'purpose','Paired checkpoint causal pilot without prediction', ...
    'seeds',seeds,'seedStatus','Previously used; exploratory only', ...
    'N',20,'k',15,'hiddenCount',7,'distanceWeightMode','equal', ...
    'strengthMode','original','bearingMeasurementStd',base.bearingMeasurementStd, ...
    'checkpointSteps',base.NSteps,'warmupSteps',warmupSteps, ...
    'durationSteps',durationSteps,'T',durationSteps*base.dt, ...
    'recoverySteps',recoverySteps,'postSteps',postSteps, ...
    'branchSteps',branchSteps,'matlabVersion',version);
fid=fopen(fullfile(out,'audit.json'),'w');
fprintf(fid,'%s',jsonencode(audit,'PrettyPrint',true));
fclose(fid);
disp(resultTable);
end

function assertCommonPrefix(control,result,warmupSteps)
fields=fieldnames(control.timeSeries);
for i=1:numel(fields)
    assert(isequaln(control.timeSeries.(fields{i})(1:warmupSteps), ...
        result.timeSeries.(fields{i})(1:warmupSteps)));
end
assert(isequaln(control.trajectory.x(:,1:warmupSteps+1),result.trajectory.x(:,1:warmupSteps+1)));
assert(isequaln(control.trajectory.y(:,1:warmupSteps+1),result.trajectory.y(:,1:warmupSteps+1)));
assert(isequaln(control.trajectory.headings(:,1:warmupSteps+1),result.trajectory.headings(:,1:warmupSteps+1)));
end
