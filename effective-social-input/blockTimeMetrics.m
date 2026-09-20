function [blocks,summary] = blockTimeMetrics(metrics,dt,blockSize,thresholds)
% Each row is a time step; columns are GO, LO, normalized pair distance.
assert(size(metrics,2)==3 && ~isempty(metrics) && all(isfinite(metrics),'all'));
validateattributes(blockSize,{'numeric'},{'scalar','integer','positive'});
validateattributes(dt,{'numeric'},{'scalar','positive','finite'});
rows=struct('firstEvaluationStep',{},'lastEvaluationStep',{},'nSteps',{}, ...
    'duration',{},'GO',{},'LO',{},'pairDistance',{},'passes',{});
for first=1:blockSize:size(metrics,1)
    last=min(first+blockSize-1,size(metrics,1));
    v=mean(metrics(first:last,:),1);
    r.firstEvaluationStep=first; r.lastEvaluationStep=last;
    r.nSteps=last-first+1; r.duration=r.nSteps*dt;
    r.GO=v(1); r.LO=v(2); r.pairDistance=v(3);
    r.passes=v(1)>=thresholds(1) && v(2)>=thresholds(2) && v(3)<=thresholds(3);
    rows(end+1)=r; %#ok<AGROW>
end
blocks=struct2table(rows);
summary.durationWeightedPass=sum(blocks.duration.*blocks.passes)/sum(blocks.duration);
summary.weightedMeans=sum(blocks.nSteps.*[blocks.GO blocks.LO blocks.pairDistance],1)/sum(blocks.nSteps);
ok=metrics(:,1)>=thresholds(1) & metrics(:,2)>=thresholds(2) & metrics(:,3)<=thresholds(3);
summary.pointFraction=mean(ok);
e=diff([false;~ok;false]); starts=find(e==1); ends=find(e==-1)-1;
summary.failureRuns=numel(starts);
summary.longestFailureTime=0;
if ~isempty(starts), summary.longestFailureTime=max(ends-starts+1)*dt; end
summary.terminalFailureCensored=~ok(end);
summary.initialFailureLeftCensored=~ok(1);
summary.recoveries=sum(diff(ok)==1);
end
