% Read pinned inputs; do not invoke A's runner which overwrites its results.
root=fileparts(mfilename('fullpath')); src=fullfile(root,'source','effective-social-input');
addpath(src,'-begin'); clean=onCleanup(@()rmpath(src));
testLocalSocialInput();
d=load(fullfile(src,'results','sq1-local-v1','e1-results.mat'));
checks=cell(numel(d.cases),1); assert(numel(d.cases)==9);
for j=1:numel(d.cases)
    r=replayLocalSocialInput(d.cfg,d.cases{j}.schedule); previous=d.results{j};
    names=fieldnames(r.timeSeries); err=0;
    for k=1:numel(names)
        err=max(err,max(abs(r.timeSeries.(names{k})-previous.timeSeries.(names{k}))));
    end
    fieldErr=max(abs(r.input(:)-previous.input(:)));
    stateErr=max(abs(r.neuralState(:)-previous.neuralState(:)));
    checks{j}=struct('name',d.cases{j}.name,'seriesError',err,'inputError',fieldErr,'stateError',stateErr);
    assert(max([err fieldErr stateErr])<1e-12);
end
alpha=(0:99)*2*pi/100;
[~,concentrated]=localSocialInputMetrics([0 0 0 0],ones(1,4)/4,alpha,.4);
[~,spread]=localSocialInputMetrics((0:3)*pi/2,ones(1,4)/4,alpha,.4);
[~,opposed]=localSocialInputMetrics([0 pi],[.5 .5],alpha,.4);
emptySupported=true; emptyMessage='';
try, localSocialInputMetrics([],[],alpha,.4); catch ME, emptySupported=false; emptyMessage=ME.message; end
% Independent scalar kernel formula on a shifted agent-specific neural grid.
angles=[.1 2.2 6.2]; weights=[.03 .09 .108]; grid=mod(alpha+.73,2*pi);
[field,metrics]=localSocialInputMetrics(angles,weights,grid,.4);
independent=zeros(100,1);
for m=1:100
    for j=1:3
        delta=atan2(sin(grid(m)-angles(j)),cos(grid(m)-angles(j)));
        independent(m)=independent(m)+weights(j)*exp(-.5*(delta/.4)^2);
    end
end
formulaError=max(abs(independent-field)); assert(formulaError<1e-12);
fresh=d.cases{8}.schedule.theta; hold=d.cases{9}.schedule.theta;
report=struct('matlabVersion',version,'testLocalSocialInputPassed',true,'caseCount',9,'caseChecks',{checks}, ...
    'independentShiftedGridError',formulaError,'weightSum',metrics.totalWeight, ...
    'concentratedEntropy',concentrated.angularEntropy,'spreadEntropy',spread.angularEntropy, ...
    'concentratedR1',abs(concentrated.z1),'spreadR1',abs(spread.z1), ...
    'opposedR1',abs(opposed.z1),'opposedHasDirection',opposed.hasDirection, ...
    'emptyVectorsSupported',emptySupported,'emptyVectorError',emptyMessage, ...
    'freshHoldFirst39ExactlyEqual',isequal(fresh(1:39),hold(1:39)), ...
    'freshHoldFirstValueDifference',hold(1)-fresh(1), ...
    'holdEntireScheduleConstant',all(hold==fresh(40)), ...
    'dispersedBumpConcentration',d.results{6}.timeSeries.bumpConcentration(end), ...
    'dispersedDecodedDirection',d.results{6}.timeSeries.decodedDirection(end));
fid=fopen(fullfile(root,'local_recheck.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(report,PrettyPrint=true)); fclose(fid);
fprintf('PASS: A test, 9 saved replays, independent kernel error %.3g. Diagnostics: opposed hasDirection=%d; empty vectors=%d; common prefix=%d.\n',formulaError,opposed.hasDirection,emptySupported,report.freshHoldFirst39ExactlyEqual);
