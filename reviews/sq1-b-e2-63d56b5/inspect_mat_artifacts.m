% Read-only inspection of pinned MAT files. Does not invoke experiment runners.
root=fileparts(mfilename('fullpath'));
base=fullfile(root,'source','effective-social-input','run_B_E2');
files=dir(fullfile(base,'run_B_E2_results','raw','*.mat'));
checks=cell(numel(files),1);
for i=1:numel(files)
    d=load(fullfile(files(i).folder,files(i).name),'rawResult');
    r=d.rawResult;
    c=struct('file',files(i).name,'seed',r.config.seed, ...
        'NSteps',r.config.NSteps,'burnIn',r.config.burnIn, ...
        'recordStride',r.config.recordStride, ...
        'timeSeriesLength',numel(r.timeSeries.globalOrder), ...
        'GO_summary',r.summary.globalOrder, ...
        'GO_recomputed',mean(r.timeSeries.globalOrder(r.config.burnIn+1:end)), ...
        'timeSeriesFields',{fieldnames(r.timeSeries)}, ...
        'topLevelFields',{fieldnames(r)});
    checks{i}=c;
end
d=load(fullfile(base,'run_B_E2_results','B_E2_runs.mat'));
csvTable=readtable(fullfile(base,'run_B_E2_results','B_E2_runs.csv'));
report=struct(); report.matlabVersion=version;
report.formalMode=d.mode; report.formalConfig=d.base;
report.formalSeeds=d.seeds; report.formalRows=height(d.results);
report.formalCsvMatMaxGODifference=max(abs(d.results.meanGO-csvTable.meanGO));
report.rawChecks=checks;
report.maxRawGORecomputeError=max(cellfun(@(x)abs(x.GO_summary-x.GO_recomputed),checks));
try
    report.runnerCodeAnalyzer=checkcode(fullfile(base,'run_B_E2_experiments.m'),'-id');
catch ME
    report.runnerCodeAnalyzerError=ME.message;
end
fid=fopen(fullfile(root,'mat_inspection.json'),'w','n','UTF-8');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
fprintf('Inspected %d raw files; formal rows=%d; raw GO error=%.3g; CSV/MAT GO error=%.3g\n', ...
    numel(files),height(d.results),report.maxRawGORecomputeError,report.formalCsvMatMaxGODifference);
