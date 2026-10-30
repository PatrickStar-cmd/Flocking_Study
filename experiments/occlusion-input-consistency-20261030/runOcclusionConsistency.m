function plan=runOcclusionConsistency(action,runId,seedOverride)
% Default only reports the plan. Explicit 'execute' starts model trajectories.
if nargin<1, action='plan'; end
root=fileparts(mfilename('fullpath')); cfg=jsondecode(fileread(fullfile(root,'config.v1.json')));
if strcmp(action,'model-smoke')
    cfg.primary=struct('profile','software_smoke','pre_steps',0,'seeds_first',2026103998,'seed_count',1);
    cfg.history_controls=cfg.history_controls([]); cfg.fresh_counts=[0;19]; cfg.post_steps=25;
    cfg.terminal_steps=10; cfg.sensitivity_terminal_steps=[5;20];
    cfg.trajectory_counts=struct('primary',2,'history_controls',0,'total',2); cfg.workers=1;
end
profiles=[cfg.primary; cfg.history_controls(:)]; [~,ix]=sort([profiles.pre_steps]); profiles=profiles(ix);
plan=struct('status','planned_not_run', ...
 'primaryTrajectories',cfg.trajectory_counts.primary,'controlTrajectories',cfg.trajectory_counts.history_controls, ...
 'S0',cfg.model.single_edge_weight*cfg.fresh_counts(:).','profiles',{profiles});
if strcmp(action,'plan'), disp(plan); return; end
assert(any(strcmp(action,{'execute','model-smoke'})),'Use plan, execute, or model-smoke.');
assert(nargin>=2 && ~isempty(regexp(runId,'^[A-Za-z0-9][A-Za-z0-9_-]*$','once')),'Provide a fresh run ID.');
if nargin<3, seedOverride=[]; end
out=fullfile(root,'results',runId); assert(~exist(out,'dir'),'A run ID is immutable; use a new ID.'); mkdir(out);
sourceDir=fullfile(out,'sources'); mkdir(sourceDir);
for name={'runOcclusionConsistency.m','movementConsistency.m','summarizeConsistency.m','config.v1.json','analyze.py','requirements.txt'}
    copyfile(fullfile(root,name{1}),fullfile(sourceDir,name{1}));
end
fid=fopen(fullfile(out,'config.json'),'w'); fprintf(fid,'%s',jsonencode(cfg,'PrettyPrint',true)); fclose(fid);
info=struct('matlab',version,'started',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),'runId',runId,'seedOverride',seedOverride, ...
 'evidence','scheduled fixed-protocol simulation; every seed and occlusion level retained');
fid=fopen(fullfile(out,'execution.json'),'w'); fprintf(fid,'%s',jsonencode(info,'PrettyPrint',true)); fclose(fid);
if isempty(gcp('nocreate')), parpool('Processes',cfg.workers); end
seedEnd=max(arrayfun(@(p)p.seeds_first+p.seed_count-1,profiles));
seeds=cfg.primary.seeds_first:seedEnd;
if ~isempty(seedOverride)
    assert(all(ismember(seedOverride,seeds)) && numel(unique(seedOverride))==numel(seedOverride));
    seeds=seedOverride(:).';
end
parfor si=1:numel(seeds)
    runSeed(seeds(si),cfg,profiles,out);
end
info.finished=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')); info.status='complete';
fid=fopen(fullfile(out,'execution.json'),'w'); fprintf(fid,'%s',jsonencode(info,'PrettyPrint',true)); fclose(fid);
if strcmp(action,'model-smoke'), verifyModelSmoke(root,out,cfg); end
end

function verifyModelSmoke(root,out,cfg)
% Future authorized model regression, never called by checkPreparation.
repo=fileparts(fileparts(root)); addpath(fullfile(repo,'effective-social-input'));
seed=cfg.primary.seeds_first; folder=fullfile(out,sprintf('seed-%d',seed));
first=load(fullfile(folder,'initial.mat')); x=first.initial;
full=load(fullfile(folder,'software_smoke-r19.mat')); zero=load(fullfile(folder,'software_smoke-r00.mat'));
old=defaultEffectiveInputConfig(); old.seed=seed; old.NSteps=cfg.post_steps; old.burnIn=0; old.recordStride=1;
old.initialXY=x.position; old.initialHeadings=x.heading; old.initialNeuralState=x.u; old.initialPreferredDirections=x.alpha.';
legacy=simulateEffectiveSocialInput(old);
neuralError=max(abs(legacy.finalNeuralState-full.finalState.u),[],'all');
positionError=max(abs(legacy.finalPosition-full.finalState.position),[],'all');
goError=max(abs(legacy.timeSeries.globalOrder-full.series(:,1)));
loError=max(abs(legacy.timeSeries.localOrder-full.series(:,2)));
assert(neuralError<1e-9 && positionError<1e-9 && goError<1e-9 && loError<1e-9);
assert(all(zero.meta.initialSocialInputPerAgent==0)&&all(zero.series(:,7)==0));
audit=struct('status','PASS','evidence','future model software smoke','neuralError',neuralError, ...
    'positionError',positionError,'GOError',goError,'LOError',loError,'zeroInputExact',true);
fid=fopen(fullfile(out,'model-smoke-checks.json'),'w');fprintf(fid,'%s',jsonencode(audit,'PrettyPrint',true));fclose(fid);
end

function runSeed(seed,cfg,profiles,out)
c=cfg.model; n=c.N; m=c.M; c.ids=mod((0:n-1).'+(1:c.full_neighbors),n)+1;
grid=(0:m-1)'*2*pi/m; d=abs(grid-grid.'); d=min(d,2*pi-d);
W=cos(pi*(d/pi).^c.connection_exponent); c.Wfft=fft(W(:,1));
rng(seed,'twister'); initial.position=rand(2,n)*c.L; initial.u=zeros(m,n); initial.heading=zeros(1,n);
for i=1:n, initial.heading(i)=2*pi*rand; initial.u(:,i)=.1*randn(m,1); end
initial.alpha=mod(grid+initial.heading,2*pi); initial.referenceAxis=initial.heading;
stream=RandStream('mt19937ar','Seed',seed+100000); order=zeros(n,c.full_neighbors);
for i=1:n, order(i,:)=randperm(stream,c.full_neighbors); end
initial.maskOrder=order; state=initial; previousStep=0; rows=struct([]);
seedDir=fullfile(out,sprintf('seed-%d',seed)); mkdir(seedDir);
save(fullfile(seedDir,'initial.mat'),'initial','-v7');
for p=1:numel(profiles)
    profile=profiles(p); if seed>=profile.seeds_first+profile.seed_count, continue; end
    assert(profile.pre_steps>=previousStep);
    for t=previousStep+1:profile.pre_steps
        I=socialField(state.position,state.alpha,true(n,c.full_neighbors),c);
        [state,~]=advance(state,I,c);
    end
    previousStep=profile.pre_steps; checkpoint=state; checkpoint.step=profile.pre_steps;
    checkpoint.maskOrder=order; checkpoint.rng=rng;
    save(fullfile(seedDir,[profile.profile '-checkpoint.mat']),'checkpoint','-v7');
    for r=cfg.fresh_counts(:).'
        visible=false(n,c.full_neighbors);
        for i=1:n, visible(i,order(i,1:r))=true; end
        x=checkpoint; z0=bearings(x.position,c); w=c.single_edge_weight*double(visible);
        % Manuscript chosen q*=1: scalar input is the edge-weight sum.
        initialInput=sum(w,2); delta=z0-x.referenceAxis.';
        vectorX=sum(w.*cos(delta),2); vectorY=sum(w.*sin(delta),2);
        field0=socialField(x.position,x.alpha,visible,c);
        initialMass=sum(field0,1).'*2*pi/m;
        series=zeros(cfg.post_steps,7); velocity=zeros(2,n,cfg.post_steps);
        position=zeros(2,n,cfg.post_steps+1); position(:,:,1)=x.position;
        for t=1:cfg.post_steps
            I=socialField(x.position,x.alpha,visible,c); [x,v]=advance(x,I,c);
            q=movementConsistency(x.position,v,c.L);
            series(t,:)=[q.GO q.LO q.C q.movingFraction q.meanSpeed q.pairDistance mean(sum(I,1)*2*pi/m)];
            velocity(:,:,t)=v; position(:,:,t+1)=x.position;
        end
        s=summarizeConsistency(series(:,1),series(:,2),cfg.terminal_steps);
        row=struct('seed',seed,'profile',profile.profile,'preSteps',profile.pre_steps,'remaining',r, ...
            'hidden',c.full_neighbors-r,'occlusionFraction',(c.full_neighbors-r)/c.full_neighbors, ...
            'initialSocialInput',mean(initialInput),'initialFieldMass',mean(initialMass), ...
            'initialVectorX',mean(vectorX),'initialVectorY',mean(vectorY), ...
            'finalConsistency',s.finalConsistency,'finalGO',s.finalGO,'finalLO',s.finalLO, ...
            'finalConsistencyQ10',s.finalConsistencyQ10,'finalJointHighFraction',s.finalJointHighFraction, ...
            'finalMovingFraction',mean(series(end-cfg.terminal_steps+1:end,4)), ...
            'finalMeanSpeed',mean(series(end-cfg.terminal_steps+1:end,5)), ...
            'finalPairDistance',mean(series(end-cfg.terminal_steps+1:end,6)));
        for window=cfg.sensitivity_terminal_steps(:).'
            ss=summarizeConsistency(series(:,1),series(:,2),window);
            row.(sprintf('C_last%d',window))=ss.finalConsistency;
        end
        rows(end+1)=row; %#ok<AGROW>
        meta=struct('seed',seed,'profile',profile.profile,'remaining',r,'visible',visible, ...
            'initialSocialInputPerAgent',initialInput,'initialVectorX',vectorX,'initialVectorY',vectorY, ...
            'initialFieldMassPerAgent',initialMass,'seriesColumns',{{'GO','LO','C','movingFraction','meanSpeed','pairDistanceNorm','fieldMass'}}, ...
            'measurementNoise',0,'postSteps',cfg.post_steps,'terminalSteps',cfg.terminal_steps);
        finalState=x;
        save(fullfile(seedDir,sprintf('%s-r%02d.mat',profile.profile,r)), ...
            'meta','series','velocity','position','finalState','-v7');
    end
end
writetable(struct2table(rows),fullfile(seedDir,'summary.csv'));
fprintf('seed %d complete (%d trajectories)\n',seed,numel(rows));
end

function z=bearings(position,c)
dx=reshape(position(1,c.ids),c.N,c.full_neighbors)-position(1,:).';
dy=reshape(position(2,c.ids),c.N,c.full_neighbors)-position(2,:).';
dx=dx-c.L*round(dx/c.L); dy=dy-c.L*round(dy/c.L); z=mod(atan2(dy,dx),2*pi);
end

function [I,z]=socialField(position,alpha,visible,c)
z=bearings(position,c); observed=z; observed(~visible)=nan;
% Only visible values enter the neural input. Missing values have zero weight.
center=observed; center(~visible)=0;
delta=mod(alpha-reshape(center,1,c.N,c.full_neighbors)+pi,2*pi)-pi;
I=c.single_edge_weight*sum(exp(-.5*(delta/c.kernel_sigma).^2).*reshape(visible,1,c.N,c.full_neighbors),3);
end

function [x,v]=advance(x,I,c)
x.u=x.u+c.dt*(-x.u+real(ifft(fft(tanh(c.beta*x.u),[],1).*c.Wfft,[],1))/c.M-c.global_inhibition+I);
y=max(0,tanh(c.beta*x.u)); v=c.speed_scale*[sum(y.*cos(x.alpha),1);sum(y.*sin(x.alpha),1)];
moving=abs(v(1,:))+abs(v(2,:))>=c.speed_scale*1e-9;
x.heading(moving)=mod(atan2(v(2,moving),v(1,moving)),2*pi);
x.position=mod(x.position+c.dt*v,c.L);
end
