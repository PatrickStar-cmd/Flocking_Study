function report=run_tests(outDir)
%RUN_TESTS No gate is marked PASS unless its assertions executed successfully.
root=fileparts(mfilename('fullpath')); addpath(root); addpath(fullfile(root,'tests'));
addpath(fullfile(root,'sources','A'),'-end');
if nargin<1, outDir=fullfile(root,'results','sq1-fourbranch-protocol12-seed1-r3-results'); end
report=repmat(struct('test_id','','config',defaultFourBranchConfig(),'seed',1,'command','',...
    'status','NOT_RUN','maxDiff',nan,'count',0,'tolerance',nan,'logPath','','details','Not executed'),15,1);
for id=1:15
    report(id).test_id=['G' num2str(id)]; report(id).command=sprintf('protocolAssertions(%d, result_directory)',id);
    report(id).logPath=['logs/G' num2str(id) '.log'];
    diary('off'); diary(fullfile(outDir,report(id).logPath));
    fprintf('%s actual execution start: %s; MATLAB %s\n',report(id).test_id,char(datetime('now')),version);
    try
        r=protocolAssertions(id,outDir);
        report(id).maxDiff=r.maxDiff; report(id).count=r.count; report(id).tolerance=r.tolerance;
        report(id).details=r.details; report(id).status='PASS';
    catch ME
        report(id).status='FAIL'; report(id).details=getReport(ME,'extended','hyperlinks','off');
        fprintf('%s\n',report(id).details);
    end
    fprintf('%s %s at %s\n',report(id).test_id,report(id).status,char(datetime('now'))); diary('off');
    writeJson(fullfile(outDir,'tests.json'),report);
end
diary(fullfile(outDir,'logs','run.log'));
end
