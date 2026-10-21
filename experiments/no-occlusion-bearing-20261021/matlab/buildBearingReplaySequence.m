function [seq,audit]=buildBearingReplaySequence(r,distribution)
% Simulator/evaluator-side reconstruction. The estimator receives only y/fresh.
cfg=r.config; n=cfg.NAgents; T=cfg.NSteps; ix=find(~eye(n)); E=numel(ix);
assert(cfg.recordStride==1 && size(r.trajectory.x,2)==T+1);
seq=struct('config',cfg,'distribution',distribution,'linearEdges',ix, ...
 'bearing',nan(T,E),'truth',zeros(T,E),'fresh',false(T,E), ...
 'outlier',false(T,E),'trueDelta',nan(T,E),'periodicSwitch',false(T,E));
noiseStream=RandStream('mt19937ar','Seed',cfg.seed+15000); mixtureState=[];
maxArchiveError=0; maxSampleError=0; lastRelative=[]; previousTruth=[];
sampleAt=1;
for t=1:T
 x=r.trajectory.x(:,t); y=r.trajectory.y(:,t);
 dx=x.'-x; dy=y.'-y;
 dx=dx-cfg.arenaSize*round(dx/cfg.arenaSize);
 dy=dy-cfg.arenaSize*round(dy/cfg.arenaSize);
 truth=mod(atan2(dy,dx),2*pi);
 if strcmp(distribution,'gaussian')
  z=randn(noiseStream,n,n); outlier=false(n);
 else
  [z,outlier,~,mixtureState]=drawNonGaussianBearingNoise(mixtureState,cfg,[n,n]);
 end
 ids=double(r.pilot.selectedHistory(:,:,t)); selected=false(n);
 for i=1:n, selected(i,ids(i,:))=true; end
 fresh=selected(ix); observed=nan(E,1);
 observed(fresh)=mod(truth(ix(fresh))+cfg.bearingMeasurementStd*z(ix(fresh)),2*pi);
 seq.truth(t,:)=truth(ix); seq.bearing(t,:)=observed;
 seq.fresh(t,:)=fresh; seq.outlier(t,:)=outlier(ix);
 relative=[dx(:),dy(:)];
 if t>1
  seq.trueDelta(t,:)=abs(atan2(sin(truth(ix)-previousTruth),cos(truth(ix)-previousTruth)));
  periodic=any(abs(relative-lastRelative)>cfg.arenaSize/2,2);
  seq.periodicSwitch(t,:)=periodic(ix);
 end
 previousTruth=truth(ix); lastRelative=relative;
 err=atan2(sin(observed(fresh)-truth(ix(fresh))),cos(observed(fresh)-truth(ix(fresh))));
 maxArchiveError=max(maxArchiveError,abs(sqrt(mean(err.^2))-r.timeSeries.observedBearingRMSE(t)));
 if sampleAt<=numel(r.pilot.samples) && r.pilot.samples{sampleAt}.step==t
  a=r.pilot.samples{sampleAt};
  assert(isequal(a.linearIndex,find(selected)));
  assert(isequal(a.standardNoise,z(selected)));
  es=atan2(sin(observed(fresh)-a.observedBearing),cos(observed(fresh)-a.observedBearing));
  maxSampleError=max(maxSampleError,max(abs(es))); sampleAt=sampleAt+1;
 end
end
if strcmp(distribution,'gaussian')
 assert(isequal(noiseStream.State,r.pilot.finalNoiseStreamState));
else
 assert(isequal(mixtureState.noiseStream.State,r.pilot.finalNoiseStreamState));
 assert(isequal(mixtureState.mixtureStream.State,r.pilot.finalMixtureStreamState));
end
assert(sampleAt==numel(r.pilot.samples)+1 && maxArchiveError<1e-12 && maxSampleError<1e-12);
assert(all(sum(seq.fresh,2)==cfg.NAgents*cfg.neighborCount));
audit=struct('status','PASS_EXACT_SENSOR_REPLAY','seed',cfg.seed, ...
 'distribution',distribution,'measurementStd',cfg.bearingMeasurementStd, ...
 'maxArchivedObservedRMSEDifference',maxArchiveError, ...
 'maxArchivedSampleAngleDifference',maxSampleError,'finalRandomStatesExact',true, ...
 'evaluationCount',sum(seq.fresh(cfg.burnIn+1:end,:),'all'));
end
