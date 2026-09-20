function runOcclusionIntegrationPilot()
testOcclusionPrediction;
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results','occlusion-integration-v1');
if ~exist(out,'dir'), mkdir(out); end
cfg=defaultEffectiveInputConfig(); cfg.NAgents=8; cfg.NNeurons=40;
cfg.NSteps=600; cfg.burnIn=200; cfg.recordStride=10;
cfg.scheduledOcclusion=true; cfg.selectionPolicy='fixed-id'; cfg.neighborCount=4;
cfg.occlusionHiddenCount=2; cfg.occlusionStartStep=301;
cfg.occlusionDurationSteps=100; cfg.distanceWeightMode='exponential';
rows=cell(6,6); idx=0;
for seed=136:137
 for mode={'none','hold','kalman'}
  cfg.seed=seed; cfg.predictionMode=mode{1};
  result=simulateEffectiveSocialInput(cfg);
  filename=fullfile(out,sprintf('seed-%d-%s.mat',seed,mode{1}));
  assert(~isfile(filename),'Output exists; use a new directory for reruns.');
  save(filename,'result','-v7');
  ix=301:400; ts=result.timeSeries; idx=idx+1;
  rows(idx,:)={seed,mode{1},mean(ts.globalOrder(ix)),mean(ts.meanInputMass(ix)),mean(ts.meanVisibleNeighbors(ix)),mean(ts.meanActiveNeighbors(ix))};
 end
end
tableOut=cell2table(rows,'VariableNames',{'seed','mode','occludedGO','inputMass','freshCount','activeCount'});
writetable(tableOut,fullfile(out,'integration.csv')); disp(tableOut);
end
