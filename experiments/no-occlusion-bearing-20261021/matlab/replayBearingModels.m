function [rows,o]=replayBearingModels(seq,methods,transition)
% All estimators read the same fresh observations on the frozen raw trajectory.
cfg=seq.config; T=cfg.NSteps; E=size(seq.bearing,2); B=cfg.burnIn; M=numel(methods);
states=cell(1,M); o.estimate=nan(T-B,E,M);
o.immHighProbability=nan(T-B,E); o.immVariance=nan(T-B,E);
o.methodNames={methods.name}; o.methods=methods; weakMeans=0;
for t=1:T
 y=seq.bearing(t,:).'; fresh=seq.fresh(t,:).';
 for m=1:M
  method=methods(m);
  if strcmp(method.type,'raw')
   estimate=y;
  elseif strcmp(method.type,'kf')
   [states{m},~,~]=task3BatchTrack(states{m},y,fresh,cfg.dt,method.q, ...
    cfg.bearingObservationVariance,Inf); estimate=states{m}.x1;
  else
   [states{m},d]=bearingIMMTrackStep(states{m},y,fresh,cfg.dt, ...
    [1e-4,method.q],cfg.bearingObservationVariance,transition);
   estimate=states{m}.x1; weakMeans=weakMeans+d.weakCircularMeans;
   if t>B && M<=5
    o.immHighProbability(t-B,:)=states{m}.mu(:,2);
    o.immVariance(t-B,:)=states{m}.p11;
   end
  end
  if t>B, o.estimate(t-B,:,m)=estimate; end
 end
end
o.finalStates=states; o.weakCircularMeans=weakMeans;
mask=seq.fresh(B+1:end,:); truth=seq.truth(B+1:end,:);
delta=seq.trueDelta(B+1:end,:); jump=mask & delta>.5;
quiet=mask & delta<.02; recent=false(size(mask));
for lag=0:10
 source=max(1,B+1-lag):T-lag;
 recent=recent | seq.trueDelta(source,:)>.5;
end
recent=recent & mask;
[eventTime,eventEdge]=find(jump); eventN=numel(eventTime);
o.events=struct('relativeStep',eventTime,'edgeColumn',eventEdge, ...
 'periodicSwitch',seq.periodicSwitch(sub2ind(size(seq.truth),eventTime+B,eventEdge)), ...
 'durationSteps',zeros(eventN,M),'recovered',false(eventN,M));
rows=struct([]);
for m=1:M
 error=atan2(sin(o.estimate(:,:,m)-truth),cos(o.estimate(:,:,m)-truth));
 a=abs(error(mask)); sorted=sort(a); count=numel(a); sq=a.^2;
 top=sort(sq,'descend'); good=mask & abs(error)<=3*cfg.bearingMeasurementStd;
 for e=1:eventN
  t=eventTime(e); edge=eventEdge(e); finish=min(size(mask,1),t+30);
  g=good(t:finish,edge);
  at=find(g(1:end-2) & g(2:end-1) & g(3:end),1);
  if isempty(at), o.events.durationSteps(e,m)=finish-t;
  else, o.events.durationSteps(e,m)=at+1; o.events.recovered(e,m)=true; end
 end
 row=struct('seed',cfg.seed,'distribution',string(seq.distribution), ...
  'measurementStd',cfg.bearingMeasurementStd,'method',string(methods(m).name), ...
  'q',methods(m).q,'count',count,'RMSE',sqrt(mean(sq)), ...
  'medianAbs',median(a),'P95Abs',sorted(ceil(.95*count)), ...
  'P99Abs',sorted(ceil(.99*count)), ...
  'topOnePercentMSEShare',sum(top(1:ceil(.01*count)))/sum(sq), ...
  'quietCount',sum(quiet,'all'),'quietRMSE',sqrt(mean(error(quiet).^2)), ...
  'jumpCount',sum(jump,'all'),'jumpRMSE',sqrt(mean(error(jump).^2)), ...
  'afterJumpCount',sum(recent,'all'),'afterJumpRMSE',sqrt(mean(error(recent).^2)), ...
  'recoveryEventCount',eventN,'recoveredEventCount',sum(o.events.recovered(:,m)), ...
  'censoredEventCount',sum(~o.events.recovered(:,m)));
 if m==1, rows=row; else, rows(m)=row; end
end
o.parameters=struct('jumpThresholdRad',.5,'quietThresholdRad',.02, ...
 'recoveryThresholdSigma',3,'recoveryConsecutiveFreshSteps',3,'recoveryHorizonSteps',30, ...
 'recoveryDefinition','third of three consecutive model steps with fresh observation and abs error <= 3 sigma; censor at 30 steps or evaluation end');
end
