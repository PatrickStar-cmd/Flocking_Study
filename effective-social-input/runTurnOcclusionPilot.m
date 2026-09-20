function runTurnOcclusionPilot(profile, seedOverride)
%RUNTURNOCCLUSIONPILOT Paired occlusion experiment with a controlled turn.
% A small leader set receives a weak rotating allocentric cue while all other
% agents must follow through the bearing-only social pathway. The control and
% every compensation method share the same checkpoint, noise stream, hidden
% edges, and turn cue. Metrics include relative degradation from control.
if nargin<1, profile="smoke"; end
profile=string(profile);
[seeds,hiddenCounts,durationSteps,maxWorkers]=grid(profile);
if nargin >= 2 && ~isempty(seedOverride)
    seeds = seedOverride(:).';
end
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','turn-occlusion-pilot-v1',char(profile));
if ~exist(out,'dir'), mkdir(out); end
ckdir=fullfile(out,'checkpoints'); if ~exist(ckdir,'dir'), mkdir(ckdir); end
seedDir=fullfile(out,'seed-results'); if ~exist(seedDir,'dir'), mkdir(seedDir); end

modes={'none','hold','mean','kalman','kalman-gm'};
warm=200; preTurn=20; post=100; recovery=300; leaders=1:3;
base=defaultEffectiveInputConfig();
base.NAgents=20; base.NNeurons=100; base.NSteps=5000; base.burnIn=3000;
base.recordStride=1; base.scheduledOcclusion=true; base.selectionPolicy='fixed-id';
base.neighborCount=15; base.occlusionStartStep=warm+1; base.predictionMode='none';
base.distanceWeightMode='equal'; base.strengthMode='original'; base.bearingMeasurementStd=.01;
 base.turnExperiment=true; base.turnLeaderIds=leaders;
base.turnCueRate=.025; base.turnCueAmplitude=.08;

seedTables=cell(numel(seeds),1);
useParallel = profile == "full" && (nargin < 2 || isempty(seedOverride));
if useParallel
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', min(maxWorkers, 4));
    end
    parfor si=1:numel(seeds)
        seedTables{si}=runOneTurnSeed(seeds(si),base,hiddenCounts,durationSteps,modes, ...
            ckdir,seedDir,warm,preTurn,post,recovery,leaders);
    end
else
    for si=1:numel(seeds)
        seedTables{si}=runOneTurnSeed(seeds(si),base,hiddenCounts,durationSteps,modes, ...
            ckdir,seedDir,warm,preTurn,post,recovery,leaders);
    end
end
R=sortrows(vertcat(seedTables{:}),{'hiddenCount','T','mode','seed'});
% Subset workers write per-seed CSVs only; the final aggregation is done after
% all workers finish so that concurrent MATLAB processes never overwrite the
% shared summary file.
if nargin < 2 || isempty(seedOverride)
    writetable(R,fullfile(out,'summary.csv'));
end
meta=struct('profile',char(profile),'seeds',seeds,'N',20,'M',100,'k',15,'leaders',leaders, ...
 'hiddenCounts',hiddenCounts,'T',durationSteps*base.dt,'dt',base.dt,'turnRate',base.turnCueRate, ...
 'turnCueAmplitude',base.turnCueAmplitude,'modes',{modes},'warmupSteps',warm,'postSteps',post, ...
 'preTurnSteps',preTurn,'turnStartBeforeOcclusion',preTurn*base.dt, ...
 'recoverySteps',recovery,'measurementStd',base.bearingMeasurementStd, ...
 'note','Controlled leader-turn intervention; relative degradation is absolute deviation from paired no-occlusion turn control normalized by control magnitude.');
if nargin < 2 || isempty(seedOverride)
    fid=fopen(fullfile(out,'metadata.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
end
end

function seedTable=runOneTurnSeed(seed,base,hiddenCounts,durationSteps,modes,ckdir,seedDir,warm,preTurn,post,recovery,leaders)
cfg=base; cfg.seed=seed;
cpPath=fullfile(ckdir,sprintf('checkpoint-%d.mat',seed));
if exist(cpPath,'file')
    z=load(cpPath);
    if isfield(z,'checkpoint')
        cp=z.checkpoint;
    else
        cfg.turnExperiment=false; cfg.turnCueAmplitude=0;
        cp=simulateEffectiveSocialInput(cfg); checkpoint=cp; save(cpPath,'checkpoint','-v7');
    end
else
    cfg.turnExperiment=false; cfg.turnCueAmplitude=0;
    cp=simulateEffectiveSocialInput(cfg); checkpoint=cp; save(cpPath,'checkpoint','-v7');
end
rows=cell(numel(hiddenCounts)*numel(durationSteps)*numel(modes),26); idx=0;
for hi=1:numel(hiddenCounts)
    h=hiddenCounts(hi); hiddenIds=makeHiddenIds(cfg.NAgents,cfg.neighborCount,h,leaders);
    for ti=1:numel(durationSteps)
        dur=durationSteps(ti); total=warm+dur+recovery;
        common=base; common.NSteps=total; common.burnIn=0; common.initialXY=cp.finalPosition;
        common.initialHeadings=cp.finalHeadings; common.initialNeuralState=cp.finalNeuralState;
        common.initialPreferredDirections=cp.finalPreferredDirections; common.seed=seed;
        common.occlusionStartStep=warm+1; common.occlusionHiddenCount=h;
        common.occlusionHiddenIds=hiddenIds; common.turnCueBaseAngles=cp.finalHeadings;
        common.turnCueStartStep=warm-preTurn+1;
        common.turnCueDurationSteps=preTurn+dur+post;
        controlCfg=common; controlCfg.occlusionDurationSteps=0; controlCfg.predictionMode='none';
        control=simulateEffectiveSocialInput(controlCfg);
        eventIx=(warm+1):(warm+dur); postIx=(warm+dur+1):(warm+dur+post);
        for mi=1:numel(modes)
            branch=common; branch.occlusionDurationSteps=dur; branch.predictionMode=modes{mi};
            result=simulateEffectiveSocialInput(branch); assertPrefix(control,result,warm);
            m=collectMetrics(control,result,eventIx,postIx,leaders,branch);
            idx=idx+1; rows(idx,:)={seed,h,h/15,dur*branch.dt,modes{mi}, ...
                m.eventGORelativeDeviation,m.postGORelativeDeviation,m.eventLORelativeDeviation, ...
                m.postPairRelativeDeviation,m.eventPairRelativeDeviation,m.eventInputMassRelativeDeviation, ...
                m.eventTurnTrackingError,m.postTurnTrackingError,m.eventFollowerLeaderError, ...
                m.controlEventGO,m.treatmentEventGO,m.controlPostPair,m.treatmentPostPair, ...
                m.eventGOLoss,m.postPairIncrease,m.eventActive,m.eventFresh,m.predictionError, ...
                m.p90PredictionError,m.recoveryTime,m.turnAngle};
        end
    end
end
names={'seed','hiddenCount','hiddenFraction','T','mode','eventGORelativeDeviation', ...
    'postGORelativeDeviation','eventLORelativeDeviation','postPairRelativeDeviation', ...
    'eventPairRelativeDeviation','eventInputMassRelativeDeviation','eventTurnTrackingError', ...
    'postTurnTrackingError','eventFollowerLeaderError','controlEventGO','treatmentEventGO', ...
    'controlPostPair','treatmentPostPair','eventGOLoss','postPairIncrease','eventActiveNeighbors', ...
    'eventFreshNeighbors','meanPredictionError','p90PredictionError','recoveryTime','turnAngle'};
seedTable=cell2table(rows(1:idx,:),'VariableNames',names);
writetable(seedTable,fullfile(seedDir,sprintf('seed-%d-summary.csv',seed)));
fprintf('turn seed %d complete\n',seed);
end

function ids=makeHiddenIds(n,k,h,leaders)
ids=zeros(n,h);
for i=1:n
 cand=mod(i-1+(1:k),n)+1;
 priority=cand(ismember(cand,leaders)); rest=cand(~ismember(cand,leaders));
 ordered=[priority rest]; ids(i,:)=ordered(1:h);
end
end

function m=collectMetrics(c,r,eventIx,postIx,leaders,cfg)
cg=c.timeSeries.globalOrder; rg=r.timeSeries.globalOrder;
cl=c.timeSeries.localOrder; rl=r.timeSeries.localOrder;
cp=c.timeSeries.meanPairDistanceNorm; rp=r.timeSeries.meanPairDistanceNorm;
cm=c.timeSeries.meanInputMass; rm=r.timeSeries.meanInputMass;
m.eventGORelativeDeviation=mean(abs(rg(eventIx)-cg(eventIx)))/max(eps,mean(abs(cg(eventIx))));
m.postGORelativeDeviation=mean(abs(rg(postIx)-cg(postIx)))/max(eps,mean(abs(cg(postIx))));
m.eventLORelativeDeviation=mean(abs(rl(eventIx)-cl(eventIx)))/max(eps,mean(abs(cl(eventIx))));
m.postPairRelativeDeviation=mean(abs(rp(postIx)-cp(postIx)))/max(eps,mean(abs(cp(postIx))));
m.eventPairRelativeDeviation=mean(abs(rp(eventIx)-cp(eventIx)))/max(eps,mean(abs(cp(eventIx))));
m.eventInputMassRelativeDeviation=mean(abs(rm(eventIx)-cm(eventIx)))/max(eps,mean(abs(cm(eventIx))));
[m.eventTurnTrackingError,m.postTurnTrackingError,m.eventFollowerLeaderError,m.turnAngle]=turnErrors(r,c,eventIx,postIx,leaders,cfg);
m.controlEventGO=mean(cg(eventIx)); m.treatmentEventGO=mean(rg(eventIx));
m.controlPostPair=mean(cp(postIx)); m.treatmentPostPair=mean(rp(postIx));
m.eventGOLoss=mean(cg(eventIx)-rg(eventIx)); m.postPairIncrease=mean(rp(postIx)-cp(postIx));
m.eventActive=mean(r.timeSeries.meanActiveNeighbors(eventIx)); m.eventFresh=mean(r.timeSeries.meanVisibleNeighbors(eventIx));
[m.predictionError,m.p90PredictionError]=predictionError(r);
m.recoveryTime=recoveryTime(rg,cg,postIx,cfg.dt);
end

function [e1,e2,ef,angle]=turnErrors(r,c,eventIx,postIx,leaders,cfg)
tr=r.trajectory.headings; tc=c.trajectory.headings; s=cfg.turnCueStartStep; rate=cfg.turnCueRate; dt=cfg.dt;
allIx=[eventIx postIx]+1; err=[]; ferr=[]; ectrl=[];
for step=allIx
 t=step-1; target=mod(cfg.turnCueBaseAngles(leaders)+rate*(t-s+1)*dt,2*pi);
 targetMean=atan2(imag(mean(exp(1i*target))),real(mean(exp(1i*target))));
 leader=atan2(imag(mean(exp(1i*tr(leaders,step)))),real(mean(exp(1i*tr(leaders,step)))));
 followers=setdiff(1:cfg.NAgents,leaders); follower=tr(followers,step);
 err(end+1)=abs(atan2(sin(leader-targetMean),cos(leader-targetMean))); %#ok<AGROW>
 ferr(end+1)=mean(abs(atan2(sin(follower-leader),cos(follower-leader)))); %#ok<AGROW>
 leaderC=atan2(imag(mean(exp(1i*tc(leaders,step)))),real(mean(exp(1i*tc(leaders,step)))));
 ectrl(end+1)=abs(atan2(sin(leaderC-targetMean),cos(leaderC-targetMean))); %#ok<AGROW>
end
n=numel(eventIx); e1=mean(err(1:n)); e2=mean(err(n+1:end)); ef=mean(ferr(1:n)); angle=rate*cfg.dt*numel(eventIx);
end

function [mu,p90]=predictionError(r)
e=[];
for t=1:numel(r.observationRecords)
    z=r.observationRecords{t};
    if isempty(z), continue; end
    mask=z.active & ~z.fresh;
    if ~any(mask(:)), continue; end
    d=atan2(sin(z.injectedBearing(mask)-z.trueBearing(mask)),cos(z.injectedBearing(mask)-z.trueBearing(mask)));
    e=[e;abs(d(:))]; %#ok<AGROW>
end
if isempty(e), mu=NaN; p90=NaN; else, mu=mean(e); q=sort(e); p90=q(max(1,ceil(.9*numel(q)))); end
end

function rec=recoveryTime(x,c,ix,dt)
ok=abs(x(ix)-c(ix))<=.05; j=find(conv(double(ok),ones(10,1),'same')>=10,1);
if isempty(j), rec=NaN; else, rec=(j-1)*dt; end
end

function assertPrefix(a,b,n)
f=fieldnames(a.timeSeries); for i=1:numel(f), assert(isequaln(a.timeSeries.(f{i})(1:n),b.timeSeries.(f{i})(1:n))); end
end

function [seeds,h,d,maxWorkers]=grid(p)
switch p
 case "smoke", seeds=201:202; h=11; d=[10 34]; maxWorkers=2;
 case "pilot", seeds=201:206; h=[7 11 13]; d=[10 20 34 60 100]; maxWorkers=6;
 case "full", seeds=201:230; h=[7 11 13]; d=[10 20 34 60 100 200]; maxWorkers=8;
 otherwise, error('profile must be smoke, pilot, or full');
end
end
