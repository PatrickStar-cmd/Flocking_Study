function resultsDir=runReviewedBatch(seeds,suffix)
% Incremental P1 runner. Results are developmental, paired by seed.
if nargin<1, seeds=2026110101:2026110105; end
if nargin<2, suffix=char(datetime('now','Format','yyyyMMdd-HHmmss')); end
validateattributes(seeds,{'numeric'},{'vector','nonempty','integer','>=',0,'<=',2^32-1});
assert(numel(unique(seeds))==numel(seeds),'Repeated seeds are not independent runs.');
assert(ischar(suffix) && ~isempty(regexp(suffix,'^[A-Za-z0-9_-]+$','once')), ...
    'Suffix must be a plain identifier.');
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
src=fullfile(root,'effective-social-input'); addpath(src); addpath(here);
gatePath=fullfile(here,'p0-reviewed-p1-v2.json');
assert(isfile(gatePath),'Run reviewed P0 before P1.');
gate=jsondecode(fileread(gatePath)); assert(strcmp(gate.status,'PASS'));
resultsDir=fullfile(src,'results',['kf-diagnostic-v1-20261101-' suffix]);
assert(~isfolder(resultsDir) && ~isfile(resultsDir),'Refuse to overwrite run directory.');
mkdir(resultsDir); mkdir(fullfile(resultsDir,'raw')); mkdir(fullfile(resultsDir,'logs'));
mkdir(fullfile(resultsDir,'sources')); mkdir(fullfile(resultsDir,'sources','core'));
mkdir(fullfile(resultsDir,'sources','diagnostic'));
logPath=fullfile(resultsDir,'logs','batch.log');
manifestPath=fullfile(resultsDir,'manifest.json');
[code,commit]=system('git rev-parse HEAD'); assert(code==0); commit=strtrim(commit);
[code,branch]=system('git branch --show-current'); assert(code==0); branch=strtrim(branch);
[code,dirty]=system('git status --short'); assert(code==0);
assert(strcmp(commit,'f0b797220db5d6f2b6ed883e3838614e67cebbe5'));
cfgBase=diagnosticConfig(seeds(1),'p1');
diagnosticJson(fullfile(resultsDir,'config.json'),cfgBase);
copyfile(gatePath,fullfile(resultsDir,'p0-reviewed-p1-v2.json'));
copyfile(fullfile(here,'p0-reviewed-fast-v2.json'),fullfile(resultsDir,'p0-reviewed-fast-v2.json'));
sources=struct('relativePath',{},'sha256',{});
for folder={src,here}
    entries=dir(fullfile(folder{1},'*.m'));
    if strcmp(folder{1},src), kind='core'; else, kind='diagnostic'; end
    for i=1:numel(entries)
        input=fullfile(folder{1},entries(i).name);
        relative=fullfile('sources',kind,entries(i).name);
        target=fullfile(resultsDir,relative); copyfile(input,target);
        sources(end+1)=struct('relativePath',strrep(relative,'\','/'), ...
            'sha256',diagnosticSha256(target)); %#ok<AGROW>
    end
end
diagnosticJson(fullfile(resultsDir,'source_hashes.json'),sources);
ledger=table(seeds(:),repmat("exploratory-development",numel(seeds),1), ...
    'VariableNames',{'seed','evidence'});
writetable(ledger,fullfile(resultsDir,'seed_ledger.csv'));
audit=struct('branch',branch,'commit',commit,'dirty',dirty,'matlabVersion',version, ...
    'platform',computer,'seeds',seeds,'N',20,'k',15,'hiddenCount',7, ...
    'maskRule','index-fixed','evidence','exploratory-development', ...
    'p0Report','p0-reviewed-p1-v2.json','modelChanged',false, ...
    'instrumentation','passive copies; exact stock regressions in review/', ...
    'stateIndexing','uRecord(:,:,t) is state after t-1; checkpoint at pre-update t0', ...
    'windowPre',[1201 1400],'windowCommon',[1401 1920],'recoverySteps',400);
diagnosticJson(fullfile(resultsDir,'audit.json'),audit);
manifest=struct('status','IN_PROGRESS','seeds',seeds,'total_conditions',9*numel(seeds), ...
    'completed_count',0,'failed_count',0,'runs',struct());
diagnosticJson(manifestPath,manifest);
conditions={'fresh',0,'none';'none',20,'none';'hold',20,'hold'; ...
    'kf-mean',20,'kf-mean';'kalman',20,'kalman';'none',120,'none'; ...
    'hold',120,'hold';'kf-mean',120,'kf-mean';'kalman',120,'kalman'};
started=tic;
for seed=seeds(:).'
    reference=[];
    for c=1:size(conditions,1)
        label=conditions{c,1}; duration=conditions{c,2}; mode=conditions{c,3};
        id=sprintf('seed%d_mode_%s_D%d',seed,label,duration);
        key=matlab.lang.makeValidName(id);
        relative=['raw/' id '.mat']; path=fullfile(resultsDir,relative);
        assert(~isfile(path));
        entry=struct('run_id',id,'seed',seed,'mode',label,'duration',duration, ...
            'mat_file',relative,'status','RUNNING');
        manifest.runs.(key)=entry; diagnosticJson(manifestPath,manifest);
        logLine(logPath,sprintf('START %d/%d %s',manifest.completed_count+1,9*numel(seeds),id));
        try
            cfg=cfgBase; cfg.seed=seed; cfg.predictionMode=mode;
            cfg.occlusionDurationSteps=duration;
            timer=tic; result=simulateEffectiveSocialInputInstrumented(cfg); seconds=toc(timer);
            if isempty(reference)
                reference=result;
            else
                assertDiagnosticPrefix(reference,result,1400);
            end
            for t=1:2200
                hidden=t>=1401 && t<1401+duration;
                expectedFresh=15-7*hidden;
                expectedActive=15;
                if hidden && (strcmp(mode,'none') || ...
                        (any(strcmp(mode,{'kf-mean','kalman'})) && t>=1501))
                    expectedActive=8;
                end
                rec=result.observationRecords{t};
                assert(all(sum(rec.fresh,2)==expectedFresh));
                assert(all(sum(rec.active,2)==expectedActive));
                p11=rec.rawVariance; p22=rec.rawP22;
                p12=rec.rawP12; p21=rec.rawP21; selected=isfinite(rec.trueBearing);
                assert(all(abs(p12(selected)-p21(selected))<1e-12));
                assert(all(p11(selected).*p22(selected)-p12(selected).*p21(selected)>=-1e-12));
            end
            run_id=id; evidence='exploratory-development'; %#ok<NASGU>
            save(path,'result','cfg','run_id','commit','evidence','seconds','-v7');
            entry.status='COMPLETED'; entry.elapsed_seconds=seconds;
            entry.sha256=diagnosticSha256(path);
            pre=result.timeSeries.globalOrder(1201:1400);
            entry.pre_GO_mean=mean(pre);
            entry.pre_GO_first100=mean(pre(1:100)); entry.pre_GO_last100=mean(pre(101:200));
            manifest.completed_count=manifest.completed_count+1;
            manifest.runs.(key)=entry; diagnosticJson(manifestPath,manifest);
            logLine(logPath,sprintf('DONE %d/%d %s %.3fs', ...
                manifest.completed_count,9*numel(seeds),id,seconds));
            clear result;
        catch ME
            entry.status='FAILED'; entry.error=getReport(ME,'extended','hyperlinks','off');
            manifest.runs.(key)=entry; manifest.status='FAILED';
            manifest.failed_count=manifest.failed_count+1;
            diagnosticJson(manifestPath,manifest);
            logLine(logPath,entry.error); rethrow(ME);
        end
    end
    clear reference;
end
manifest.status='COMPLETED'; manifest.total_elapsed_seconds=toc(started);
diagnosticJson(manifestPath,manifest);
fprintf('PASS P1 batch %d unique trajectories: %s\n',manifest.completed_count,resultsDir);
end

function logLine(path,message)
line=sprintf('[%s] %s\n',char(datetime('now','TimeZone','Asia/Shanghai')),message);
fprintf('%s',line); fid=fopen(path,'a','n','UTF-8'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s',line);
end
