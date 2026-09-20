function runPaperKFSimulation()
%RUNPAPERKFSIMULATION Four-method paired intervention matching the draft PDF.
% Methods: fresh-drop, hold-last, KF-mean, KF-uncertainty.
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','paper-kf-simulation');
if ~exist(out,'dir'), mkdir(out); end
seeds=101:105; durationSteps=[10 34 100];
warmupSteps=200; recoverySteps=300; branchSteps=warmupSteps+max(durationSteps)+recoverySteps;
eventStart=warmupSteps+1; postSteps=100;
modes={'none','hold','mean','kalman'};
base=defaultEffectiveInputConfig();
base.NAgents=20; base.NNeurons=100; base.NSteps=5000; base.burnIn=3000;
base.recordStride=1; base.scheduledOcclusion=true; base.selectionPolicy='fixed-id';
base.neighborCount=15; base.occlusionHiddenCount=7; base.occlusionStartStep=2;
base.occlusionDurationSteps=0; base.predictionMode='none'; base.distanceWeightMode='equal';
base.strengthMode='original'; base.bearingMeasurementStd=0.01;
rows={}; idx=0;
for seed=seeds
  base.seed=seed;
  checkpoint=simulateEffectiveSocialInput(base);
  save(fullfile(out,sprintf('checkpoint-%d.mat',seed)),'checkpoint','-v7');
  branch=base; branch.NSteps=branchSteps; branch.burnIn=0;
  branch.initialXY=checkpoint.finalPosition; branch.initialHeadings=checkpoint.finalHeadings;
  branch.initialNeuralState=checkpoint.finalNeuralState; branch.initialPreferredDirections=checkpoint.finalPreferredDirections;
  branch.occlusionStartStep=eventStart;
  branch.occlusionDurationSteps=0; branch.predictionMode='none';
  control=simulateEffectiveSocialInput(branch);
  for Tsteps=durationSteps
    eventIx=eventStart:(eventStart+Tsteps-1);
    postIx=(eventStart+Tsteps):(eventStart+Tsteps+postSteps-1);
    for mi=1:numel(modes)
      branch.occlusionDurationSteps=Tsteps; branch.predictionMode=modes{mi};
      result=simulateEffectiveSocialInput(branch);
      assertCommonPrefixPaper(control,result,warmupSteps);
      save(fullfile(out,sprintf('seed-%d-%s-T%g.mat',seed,modes{mi},Tsteps*base.dt)),'result','-v7');
      idx=idx+1;
      rows(idx,:)={seed,modes{mi},Tsteps*base.dt,mean(result.timeSeries.globalOrder), ...
       mean(control.timeSeries.globalOrder(eventIx)-result.timeSeries.globalOrder(eventIx)), ...
       mean(control.timeSeries.meanPairDistanceNorm(postIx)-result.timeSeries.meanPairDistanceNorm(postIx)), ...
       mean(result.timeSeries.meanVisibleNeighbors(eventIx)),mean(result.timeSeries.meanActiveNeighbors(eventIx)), ...
       mean(result.timeSeries.meanInputMass(eventIx)),mean(result.timeSeries.meanInputConcentration(eventIx))};
    end
    fprintf('seed %d T %.1f complete\n',seed,Tsteps*base.dt);
  end
end
R=cell2table(rows,'VariableNames',{'seed','mode','T','fullGO','eventGOLoss','postPairDistanceIncrease','eventFreshNeighbors','eventActiveNeighbors','eventInputMass','eventInputConcentration'});
writetable(R,fullfile(out,'summary.csv'));
meta=struct('model','paper draft v0.1','N',20,'M',100,'dt',base.dt,'v_s',base.speedScale,'sigma_s',base.receptiveFieldWidth,'seeds',seeds,'T',durationSteps*base.dt,'modes',{modes},'hiddenCount',7,'k',15,'measurementStd',base.bearingMeasurementStd,'note','Exploratory paired closed-loop simulation; thresholds and formal success endpoint remain to be frozen.');
fid=fopen(fullfile(out,'metadata.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
end

function assertCommonPrefixPaper(control,result,n)
fields=fieldnames(control.timeSeries);
for i=1:numel(fields)
 assert(isequaln(control.timeSeries.(fields{i})(1:n),result.timeSeries.(fields{i})(1:n)));
end
end
