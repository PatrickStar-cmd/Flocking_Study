function reviewInstrumentation()
% Codex-owned regression independent of AGY's P0 implementation.
here=fileparts(mfilename('fullpath')); expDir=fileparts(here);
root=fileparts(fileparts(expDir));
restoredefaultpath;
addpath(fullfile(root,'effective-social-input')); addpath(expDir);
assert(~isfile(fullfile(here,'instrumentation-regression.json')));
cfg=defaultEffectiveInputConfig;
cfg.seed=2026110101; cfg.NSteps=80; cfg.burnIn=10;
cfg.recordStride=1; cfg.scheduledOcclusion=true;
cfg.selectionPolicy='fixed-id'; cfg.neighborCount=15;
cfg.occlusionHiddenCount=7; cfg.occlusionStartStep=31;
cfg.occlusionDurationSteps=20;
modes={'none','hold','kf-mean','kalman'};
elapsed=zeros(1,numel(modes));
for i=1:numel(modes)
    cfg.predictionMode=modes{i};
    stock=simulateEffectiveSocialInput(cfg);
    timer=tic; instrumented=simulateEffectiveSocialInputInstrumented(cfg);
    elapsed(i)=toc(timer);
    compareStock(stock,instrumented);
    fprintf('PASS passive instrumentation short mode=%s seconds=%.3f\n',modes{i},elapsed(i));
end
reference=load(fullfile(here,'stock_full_reference.mat'));
timer=tic;
instrumented=simulateEffectiveSocialInputInstrumented(reference.cfg);
fullSeconds=toc(timer);
compareStock(reference.result,instrumented);
save(fullfile(here,'instrumented_full_reference.mat'),'instrumented','fullSeconds','-v7');
report=struct('status','PASS','modes',{modes},'shortSeconds',elapsed, ...
    'fullSteps',2200,'fullMode','kalman','fullSeconds',fullSeconds, ...
    'comparison','Exact isequaln for all stock result fields and observation fields');
fid=fopen(fullfile(here,'instrumentation-regression.json'),'w');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
fprintf('PASS full instrumented versus stock, exact; seconds=%.3f\n',fullSeconds);
end

function compareStock(a,b)
fields=fieldnames(a);
for f=1:numel(fields)
    key=fields{f};
    if strcmp(key,'observationRecords')
        assert(numel(a.(key))==numel(b.(key)));
        for t=1:numel(a.(key))
            original=a.(key){t}; added=b.(key){t};
            names=fieldnames(original);
            for k=1:numel(names)
                assert(isequaln(original.(names{k}),added.(names{k})), ...
                    'Observation mismatch at step %d field %s',t,names{k});
            end
        end
    else
        assert(isequaln(a.(key),b.(key)),'Stock field mismatch: %s',key);
    end
end
end
