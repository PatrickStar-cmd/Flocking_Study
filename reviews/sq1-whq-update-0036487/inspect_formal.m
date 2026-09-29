% Read saved results only; never invoke B's runner or modify source data.
root=fileparts(mfilename('fullpath'));
tab=readtable(fullfile(root,'source','run_B_1_03_summary_results','B_E2_runs.csv'),'TextType','string');
files=dir(fullfile(root,'raw','*.mat'));
assert(numel(files)==480 && height(tab)==480);
csvNames={'meanGO','meanLO','meanPairDistanceNorm','meanNearestDistanceNorm','meanInputAmplitude','meanActiveNeighbors','meanAngularCoverage','meanInputConcentration','meanFieldPeak','meanBumpConcentration','meanSpeed','meanInputMass'};
tsNames=csvNames; tsNames{1}='globalOrder'; tsNames{2}='localOrder';
maxMeanError=zeros(size(csvNames)); maxSummaryError=0; maxDistanceError=0;
maxAmplitudeError=0; maxCountError=0; checks=cell(numel(files),1);
for j=1:numel(files)
    d=load(fullfile(files(j).folder,files(j).name),'rawResult'); r=d.rawResult; c=r.config;
    assert(c.NSteps==8000 && c.burnIn==6000 && c.recordStride==100 && c.NAgents==20 && c.NNeurons==100);
    assert(c.dt==0.3 && ~c.scheduledOcclusion && ~c.occlusionEnabled && strcmp(c.visibilityMode,'ideal'));
    idx=find(tab.seed==c.seed & tab.k==c.neighborCount & tab.strengthMode==string(c.strengthMode) & tab.selectionPolicy==string(c.selectionPolicy));
    prefix=extractBefore(string(files(j).name),'_'+string(c.strengthMode)+'_');
    idx=idx(tab.experiment(idx)==prefix); assert(isscalar(idx));
    for m=1:numel(csvNames)
        v=r.timeSeries.(tsNames{m}); assert(numel(v)==8000 && all(isfinite(v)));
        value=mean(v(6001:8000));
        maxMeanError(m)=max(maxMeanError(m),abs(value-tab.(csvNames{m})(idx)));
        maxSummaryError=max(maxSummaryError,abs(value-r.summary.(tsNames{m})));
    end
    assert(isequal(size(r.trajectory.x),[20 81]));
    assert(max(abs(r.trajectory.time-(0:100:8000)*c.dt))<1e-10);
    % Independent periodic pair-distance calculation at every saved nonzero frame.
    for t=2:81
        xy=[r.trajectory.x(:,t),r.trajectory.y(:,t)];
        total=0; count=0;
        for a=1:19
            for b=a+1:20
                delta=abs(xy(a,:)-xy(b,:)); delta=min(delta,c.arenaSize-delta);
                total=total+sqrt(sum(delta.^2)); count=count+1;
            end
        end
        maxDistanceError=max(maxDistanceError,abs(total/count/c.arenaSize-r.timeSeries.meanPairDistanceNorm((t-1)*100)));
    end
    if strcmp(c.strengthMode,'original'), amplitude=c.neighborCount*c.totalSocialAttraction/c.NAgents;
    else, amplitude=(c.NAgents-1)*c.totalSocialAttraction/c.NAgents; end
    maxAmplitudeError=max(maxAmplitudeError,max(abs(r.timeSeries.meanInputAmplitude-amplitude)));
    maxCountError=max(maxCountError,max(abs(r.timeSeries.meanActiveNeighbors-c.neighborCount)));
    checks{j}=struct('file',files(j).name,'seed',c.seed,'steps',c.NSteps,'burnIn',c.burnIn,'trajectoryFrames',size(r.trajectory.x,2));
end
report=struct('matlabVersion',version,'rawCount',numel(files),'verifiedTimeSeriesValues',numel(files)*8000*numel(csvNames), ...
    'csvMetricNames',{csvNames},'maxCsvMeanErrors',maxMeanError,'maxSummaryError',maxSummaryError, ...
    'maxIndependentPairDistanceError',maxDistanceError,'distanceFramesChecked',480*80, ...
    'maxAmplitudeError',maxAmplitudeError,'maxActiveCountError',maxCountError, ...
    'timeSeriesFields',{fieldnames(r.timeSeries)},'topLevelFields',{fieldnames(r)},'rawChecks',{checks});
assert(max(maxMeanError)<1e-12 && maxSummaryError<1e-12 && maxDistanceError<1e-12 && maxAmplitudeError<1e-12 && maxCountError==0);
d=load(fullfile(root,'source','run_B_1_03_summary_results','B_E2_runs.mat'));
report.csvMatSummaryMaxGOError=max(abs(d.results.meanGO-tab.meanGO));
assert(report.csvMatSummaryMaxGOError<1e-12);
fid=fopen(fullfile(root,'mat_audit.json'),'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
fprintf('PASS: %d formal MAT files, 12 time-series means, %d independent distance frames; max CSV mean error %.3g, distance error %.3g.\n',numel(files),480*80,max(maxMeanError),maxDistanceError);
