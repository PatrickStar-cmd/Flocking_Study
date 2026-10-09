function outDir=run_protocol12_seed1_r3(outputName)
%RUN_PROTOCOL12_SEED1_R3 New output only; seed1 is software evidence.
root=fileparts(mfilename('fullpath')); addpath(root); addpath(fullfile(root,'tests')); addpath(fullfile(root,'sources','A'),'-end');
if nargin<1, outputName='sq1-fourbranch-protocol12-seed1-r3-results'; end
outDir=fullfile(root,'results',outputName); createOutputDirectory(outDir);
diary(fullfile(outDir,'logs','run.log')); cleanup=onCleanup(@()diary('off'));
started=datetime('now'); wall=tic;
try
    cfg=defaultFourBranchConfig(); writeJson(fullfile(outDir,'config.json'),cfg);
    source=buildSourceManifest(root,outDir);
    fprintf('MATLAB %s; seed1 N20 M100 k8 h4 random d34\n',version);
    c=cfg; c.recordFromStep=6000;
    fprintf('Uninterrupted 0 -> 6600 reference (recording 6001 -> 6600).\n');
    [continuous,~]=simulateProtocolSegment(c,initializeCoreState(c),0,6600);
    cp6000=continuous.checkpoints.step_6000; cp6200=continuous.checkpoints.step_6200;
    assert(isscalar(cp6000)&&isscalar(cp6200));
    assert(numel(cp6000.maskDiagnostics)==20 && numel(cp6200.maskDiagnostics)==20);
    save(fullfile(outDir,'checkpoints.mat'),'cp6000','cp6200','-v7.3');
    modes={'fresh','drop','drop-normalized','hold'}; branchTimes=zeros(4,1);
    for q=1:4
        c=cfg; c.mode=modes{q}; clock=tic;
        fprintf('Branch %s from cp6000.\n',modes{q});
        seg=resumeProtocolFromCheckpoint(cp6000,c,6600);
        run.mode=modes{q}; run.seed=cfg.seed; run.config=c;
        run.run_id=[modes{q} '_seed1_random_d34']; run.reference_id='fresh_seed1';
        run.checkpoint6000=cp6000; run.checkpoint6200=seg.checkpoints.step_6200;
        run.segment=seg; run.source=source;
        save(fullfile(outDir,[modes{q} '_seed1.mat']),'run','-v7.3'); branchTimes(q)=toc(clock);
        if q==1, fresh=run; end
    end
    save(fullfile(outDir,'continuous_reference.mat'),'continuous','-v7.3');
    fprintf('T=0 four branch test data.\n'); t0=struct();
    for q=1:4
        c=cfg; c.mode=modes{q}; c.occlusionDurationSteps=0; c.recoveryStartStep=c.occlusionStartStep;
        t0.(matlab.lang.makeValidName(modes{q}))=resumeProtocolFromCheckpoint(cp6000,c,6600);
    end
    save(fullfile(outDir,'T0_results.mat'),'t0','-v7.3');
    fprintf('Fresh-reference same-state replay.\n'); replay=buildReferenceReplay(fresh);
    save(fullfile(outDir,'reference_same_state_replay.mat'),'replay','-v7.3');
    writetable(generateSummary(outDir),fullfile(outDir,'summary.csv'));
    rows=struct([]);
    for q=1:4
        rows(q).run_id=[modes{q} '_seed1_random_d34']; rows(q).seed=1; rows(q).mode=modes{q};
        rows(q).mask_rule='random'; rows(q).duration_steps=34; rows(q).reference_id='fresh_seed1';
        rows(q).checkpoint_id='seed1_step6000'; rows(q).file_path=[modes{q} '_seed1.mat'];
        rows(q).status='COMPLETED'; rows(q).error_reason=''; rows(q).is_shared_reference=q==1;
    end
    writetable(struct2table(rows),fullfile(outDir,'RUN_MANIFEST.csv'));
    for q=1:4
        rows(q).run_id=[modes{q} '_seed1_T0']; rows(q).duration_steps=0;
        rows(q).file_path=['T0_results.mat:t0.' matlab.lang.makeValidName(modes{q})];
    end
    writetable(struct2table(rows),fullfile(outDir,'T0_TEST_MANIFEST.csv'));
    fprintf('Run G1-G15 assertions.\n'); report=run_tests(outDir);
    metadata.protocol_version='1.2'; metadata.status='GENERATED_PENDING_C_REVIEW';
    metadata.formalBatchAllowed=false; metadata.seed=1; metadata.protocolSHA=source.protocolSHA;
    metadata.configSHA=source.configSHA; metadata.shared_reference='fresh_seed1';
    metadata.failedTests={report(strcmp({report.status},'FAIL')).test_id};
    writeJson(fullfile(outDir,'metadata.json'),metadata);
    audit.matlabVersion=version; audit.workerCount=0; audit.execution='serial';
    audit.command=['restoredefaultpath; addpath(package_root); run_protocol12_seed1_r3(''' outputName ''')'];
    audit.startTime=char(started); audit.endTime=char(datetime('now')); audit.wallSeconds=toc(wall);
    audit.branchSeconds=branchTimes; audit.actualRoot=root; audit.overwriteProtection=true;
    writeJson(fullfile(outDir,'audit.json'),audit);
    fprintf('Finished: %d PASS, %d FAIL, %d NOT_RUN.\n',sum(strcmp({report.status},'PASS')),sum(strcmp({report.status},'FAIL')),sum(strcmp({report.status},'NOT_RUN')));
catch ME
    failure.status='FAILED'; failure.time=char(datetime('now')); failure.error= getReport(ME,'extended','hyperlinks','off');
    writeJson(fullfile(outDir,'failure_logs','runner_failure.json'),failure); rethrow(ME);
end
end
