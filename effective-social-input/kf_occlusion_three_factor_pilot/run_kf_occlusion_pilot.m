function outDir=run_kf_occlusion_pilot(outputName,profile)
% 一次性顺序执行，拒绝覆盖历史结果。正式结果只来自profile='pilot'。
root=fileparts(mfilename('fullpath')); addpath(root,'-begin');
if nargin<1, outputName=['pilot_' datestr(now,'yyyymmdd_HHMMSS')]; end
if nargin<2, profile='pilot'; end
outDir=fullfile(root,'results',outputName); assert(~exist(outDir,'dir'),'KF:Exists','Use a new output name');
mkdir(outDir); mkdir(fullfile(outDir,'logs')); mkdir(fullfile(outDir,'references')); mkdir(fullfile(outDir,'conditions'));
diary(fullfile(outDir,'logs','run.log')); cleanup=onCleanup(@()diary('off')); wall=tic;
cfg=defaultKFExperimentConfig(); conditions=pilotConditions();
if strcmp(profile,'smoke'), cfg.seeds=1; conditions=conditions(2); end
cfg.profile=profile; cfg.outputName=outputName;
writeJson(fullfile(outDir,'config.json'),cfg); writeJson(fullfile(outDir,'conditions.json'),conditions);
source=sourceManifest(root); writeJson(fullfile(outDir,'source_manifest.json'),source);
manifest=struct([]); closed=struct([]); same=struct([]); offline=struct([]); masks=struct([]);
try
    fprintf('Start %s; MATLAB %s; profile %s\n',datestr(now,31),version,profile);
    checks=run_software_checks(fullfile(outDir,'logs')); %#ok<NASGU>
    for seed=cfg.seeds
        levels=unique([conditions.A0]);
        for level=levels
            c=cfg; c.seed=seed; c.A0=level; c.w0=level/c.k; c.mode='fresh';
            fprintf('Reference seed%d A0 %.3f: 0 -> 6600\n',seed,level); timer=tic;
            fresh=simulateKFSegment(c,initializeCoreState(c),0,c.NSteps);
            cpbase=fresh.checkpoints.step_6200;
            refid=sprintf('fresh_seed%d_A%d',seed,round(level*1000));
            save(fullfile(outDir,'references',[refid '.mat']),'fresh','cpbase','-v7.3');
            entry=manifestEntry(refid,seed,level,0,'reference','fresh',['references/' refid '.mat'],toc(timer));
            manifest=[manifest entry]; %#ok<AGROW>
            writeJson(fullfile(outDir,'RUN_MANIFEST.json'),manifest);
            for ci=find([conditions.A0]==level)
                cond=conditions(ci); c.h=cond.h; c.maskRule=cond.maskRule;
                [mask,md]=makeProtocolMasks(cpbase.position,c,cpbase.inputState.lastBearing);
                cp=cpbase; cp.config=c; cp.inputState.mask=mask; cp.inputState.maskDiagnostics=md; cp.inputState.maskReady=true;
                key=sprintf('seed%d_%s',seed,cond.id); folder=fullfile(outDir,'conditions',key); mkdir(folder);
                save(fullfile(folder,'checkpoint6200.mat'),'cp','cond','-v7.3');
                for i=1:c.NAgents
                    rr.seed=seed; rr.condition_id=cond.id; rr.agent=i; rr.h=c.h; rr.structure=cond.label;
                    rr.A0=level; rr.hiddenIds=strtrim(sprintf('%d ',mask(i,:)));
                    rr.R1=md{i}.hidden.R1; rr.R2=md{i}.hidden.R2; rr.maxGap=md{i}.hidden.maxGap;
                    rr.maskBearingStep=6200; masks=[masks rr]; %#ok<AGROW>
                end
                fprintf('Offline + same-state %s\n',key);
                [diagnostic,dr]=run_kf_prediction_diagnostics(fresh,cp,c,cond);
                save(fullfile(folder,'prediction_diagnostics.mat'),'diagnostic','-v7.3'); offline=[offline dr]; %#ok<AGROW>
                [replay,sr]=run_same_state_input_experiment(fresh,cp,c,cond);
                save(fullfile(folder,'same_state_replay.mat'),'replay','-v7.3'); same=[same sr]; clear replay diagnostic;
                c.mode='fresh'; rows=summarizeClosedLoop(fresh,fresh,c,cond,refid); closed=[closed rows]; %#ok<AGROW>
                for mi=2:5
                    c.mode=c.modes{mi}; runid=[key '_' c.mode]; fprintf('Closed-loop %s\n',runid); timer=tic;
                    seg=run_closed_loop_experiment(cp,c);
                    save(fullfile(folder,[c.mode '.mat']),'seg','-v7.3');
                    rows=summarizeClosedLoop(seg,fresh,c,cond,runid); closed=[closed rows]; %#ok<AGROW>
                    entry=manifestEntry(runid,seed,level,c.h,cond.label,c.mode,['conditions/' key '/' c.mode '.mat'],toc(timer));
                    manifest=[manifest entry]; clear seg; %#ok<AGROW>
                    writeJson(fullfile(outDir,'RUN_MANIFEST.json'),manifest);
                end
                writetable(struct2table(closed),fullfile(outDir,'closed_loop_summary.csv'));
                writetable(struct2table(same),fullfile(outDir,'same_state_summary.csv'));
                writetable(struct2table(offline),fullfile(outDir,'prediction_diagnostics_summary.csv'));
                writetable(struct2table(masks),fullfile(outDir,'hidden_sets.csv'));
            end
            clear fresh cpbase;
        end
    end
    writetable(struct2table(manifest),fullfile(outDir,'RUN_MANIFEST.csv'));
    audit.status='COMPLETE'; audit.elapsedSeconds=toc(wall); audit.actualTrajectoryCount=numel(manifest);
    audit.expectedTrajectoryCount=numel(cfg.seeds)*(4*numel(conditions)+numel(unique([conditions.A0])));
    assert(audit.actualTrajectoryCount==audit.expectedTrajectoryCount);
    audit.profile=profile; audit.seeds=cfg.seeds; audit.workerCount=1; audit.MATLAB=version;
    audit.exclusions='none'; audit.analysisUnit='seed'; audit.parameterSelection=cfg.parameterSelection;
    audit.KF_uncertainty='not implemented in this pilot; covariance diagnostic only';
    audit.factorialInteractions='not tested: seven one-factor conditions only';
    writeJson(fullfile(outDir,'audit.json'),audit);
    fprintf('COMPLETE: %d trajectories; %.1f seconds\n',numel(manifest),toc(wall));
catch ME
    failure.status='FAILED'; failure.message=getReport(ME,'extended','hyperlinks','off'); failure.elapsedSeconds=toc(wall);
    failure.lastManifest=manifest; writeJson(fullfile(outDir,'failure.json'),failure);
    fprintf('%s\n',failure.message); rethrow(ME);
end
end
function entry=manifestEntry(id,seed,A0,h,structure,mode,file,elapsed)
entry=struct('run_id',id,'seed',seed,'A0',A0,'h',h,'structure',structure,'mode',mode,...
    'file',file,'elapsedSeconds',elapsed,'status','COMPLETE');
end
