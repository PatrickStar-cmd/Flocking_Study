function result=diagnose_noisy_kf_rmse_20261014(root,outDir)
% Replay archived KF trajectories and exact sensor stream; no dynamics rerun.
assert(~isfolder(outDir)); mkdir(outDir);
source=fullfile(root,'source-snapshot'); addpath(source);
files=dir(fullfile(root,'seed-*-std-010-pair.mat'));
assert(numel(files)==6); rows=struct([]);
maxRMSEError=0; maxNoiseStateDifference=0; examples=cell(1,6);
for fi=1:numel(files)
 loaded=load(fullfile(files(fi).folder,files(fi).name)); r=loaded.pair.kf; cfg=r.config;
 n=cfg.NAgents; ix=find(~eye(n)); [recv,nei]=ind2sub([n n],ix);
 stream=RandStream('mt19937ar','Seed',cfg.seed+15000); tr=[];
 evalSteps=(cfg.burnIn+1):cfg.NSteps; E=cfg.NAgents*cfg.neighborCount;
 errors=zeros(numel(evalSteps),E); rawErrors=errors; gains=errors;
 priorErrors=errors; recentJump=false(size(errors)); absDelta=errors;
 previousRelative=[]; jumpHistory=false(n*n,11); traces=cell(numel(evalSteps),1);
 for step=1:cfg.NSteps
  pos=[r.trajectory.x(:,step).';r.trajectory.y(:,step).'];
  rx=pos(1,:)-pos(1,:).'; ry=pos(2,:)-pos(2,:).';
  rx=rx-cfg.arenaSize*round(rx/cfg.arenaSize);
  ry=ry-cfg.arenaSize*round(ry/cfg.arenaSize);
  truth=mod(atan2(ry,rx),2*pi);
  noise=cfg.bearingMeasurementStd*randn(stream,n,n);
  ids=double(r.pilot.selectedHistory(:,:,step)); selected=false(n);
  for i=1:n, selected(i,ids(i,:))=true; end
  fresh=selected(ix); y=nan(size(ix)); y(fresh)=mod(truth(ix(fresh))+noise(ix(fresh)),2*pi);
  if isempty(tr)
   old=false(size(ix)); predicted=zeros(size(ix));
  else
   old=tr.has; predicted=mod(tr.x1+cfg.dt*tr.x2,2*pi);
  end
  [tr,~,~]=task3BatchTrack(tr,y,fresh,cfg.dt,cfg.bearingProcessNoise,cfg.bearingObservationVariance,Inf);
  err=atan2(sin(tr.x1(fresh)-truth(ix(fresh))),cos(tr.x1(fresh)-truth(ix(fresh))));
  obsErr=atan2(sin(y(fresh)-truth(ix(fresh))),cos(y(fresh)-truth(ix(fresh))));
  archivedRMSE=r.timeSeries.usedBearingRMSE(step);
  maxRMSEError=max(maxRMSEError,abs(sqrt(mean(err.^2))-archivedRMSE));
  nowRelative=[rx(:),ry(:)]; jump=false(n*n,1);
  if ~isempty(previousRelative)
   jump=any(abs(nowRelative-previousRelative)>cfg.arenaSize/2,2);
  end
  jumpHistory=[jumpHistory(:,2:end),jump]; previousRelative=nowRelative;
  if step>cfg.burnIn
   ei=step-cfg.burnIn;
   errors(ei,:)=err; rawErrors(ei,:)=obsErr;
   gains(ei,:)=tr.p11(fresh)/cfg.bearingObservationVariance;
   priorErrors(ei,:)=atan2(sin(predicted(fresh)-truth(ix(fresh))),cos(predicted(fresh)-truth(ix(fresh))));
   recentJump(ei,:)=any(jumpHistory(ix(fresh),:),2);
   if step>1
    priorPos=[r.trajectory.x(:,step-1).';r.trajectory.y(:,step-1).'];
    dx=priorPos(1,:)-priorPos(1,:).';dy=priorPos(2,:)-priorPos(2,:).';
    dx=dx-cfg.arenaSize*round(dx/cfg.arenaSize);dy=dy-cfg.arenaSize*round(dy/cfg.arenaSize);
    oldTruth=mod(atan2(dy,dx),2*pi);
    absDelta(ei,:)=abs(atan2(sin(truth(ix(fresh))-oldTruth(ix(fresh))),cos(truth(ix(fresh))-oldTruth(ix(fresh)))));
   end
   traces{ei}=ix(fresh);
  end
 end
 assert(maxRMSEError<1e-12,'Archived posterior errors not exactly replayed');
 assert(isequal(stream.State,r.pilot.finalNoiseStreamState),'Noise stream mismatch');
 err=errors(:); obs=rawErrors(:); sq=err.^2; [sorted,order]=sort(sq,'descend');
 topOne=sum(sorted(1:ceil(numel(err)*.01)))/sum(sq);
 jump=recentJump(:); prior=priorErrors(:); gain=gains(:);
 row=struct('seed',cfg.seed,'measurementStd',cfg.bearingMeasurementStd, ...
  'count',numel(err),'rawRMSE',sqrt(mean(obs.^2)),'kfRMSE',sqrt(mean(sq)), ...
  'rawMedianAbs',median(abs(obs)),'kfMedianAbs',median(abs(err)), ...
  'kfP95Abs',quantile(abs(err),.95),'kfP99Abs',quantile(abs(err),.99), ...
  'kfBetterFraction',mean(abs(err)<abs(obs)), ...
  'topOnePercentSquaredErrorShare',topOne,'recentJumpCount',sum(jump), ...
  'recentJumpSquaredErrorShare',sum(sq(jump))/sum(sq), ...
  'noRecentJumpRMSE',sqrt(mean(sq(~jump))), ...
  'meanGain',mean(gain),'priorRMSE',sqrt(mean(prior.^2)), ...
  'predictionPartMSE',mean(((1-gain).*prior).^2), ...
  'measurementPartMSE',mean((gain.*obs).^2), ...
  'crossTerm',mean(2*(1-gain).*prior.*gain.*obs));
 if fi==1, rows=row; else, rows(fi)=row; end
 % Store one worst-error trace identity for each seed; no sample is removed.
 [at,edge]=ind2sub(size(errors),order(1)); linear=traces{at}(edge);
 [receiver,neighbor]=ind2sub([n n],linear);
 examples{fi}=struct('seed',cfg.seed,'step',evalSteps(at),'receiver',receiver, ...
  'neighbor',neighbor,'posteriorError',errors(at,edge),'observationError',rawErrors(at,edge), ...
  'recentJump',recentJump(at,edge),'currentTrueStepChange',absDelta(at,edge));
 fprintf('seed=%d rawRMSE=%.6f KF=%.6f median=%.6f better=%.3f top1%%MSE=%.3f jumpsMSE=%.3f\n', ...
  row.seed,row.rawRMSE,row.kfRMSE,row.kfMedianAbs,row.kfBetterFraction,topOne,row.recentJumpSquaredErrorShare);
 % Full arrays support a second independent count/statistic implementation.
 replayFile=fullfile(outDir,sprintf('seed-%d-replay-errors.mat',cfg.seed));
 save(replayFile,'errors','rawErrors','gains','priorErrors','recentJump','absDelta','-v7');
end
writetable(struct2table(rows),fullfile(outDir,'diagnostics.csv'));
result=struct('status','PASS_EXACT_FILTER_REPLAY','sourceDirectory',root, ...
 'interpretation','post-hoc exploratory diagnosis; no parameter tuning or dynamics rerun', ...
 'maxArchivedRMSEError',maxRMSEError,'rows',rows,'worstExamples',{examples}, ...
 'recentJumpDefinition','minimum-image relative displacement changed by more than L/2 now or in prior 10 steps; includes 11 instants');
fid=fopen(fullfile(outDir,'replay-audit.json'),'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
