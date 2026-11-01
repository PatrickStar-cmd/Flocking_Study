function report=runReviewedP0(profile,outPath)
% Reviewed eight-point software gate. No scientific superiority is tested here.
if nargin<1, profile='fast'; end
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
addpath(here); addpath(fullfile(root,'effective-social-input'));
if nargin<2, outPath=fullfile(here,['p0-reviewed-' profile '-v2.json']); end
assert(~isfile(outPath),'Refuse to overwrite P0 report.');
report=struct('profile',profile,'status','IN_PROGRESS','matlabVersion',version, ...
    'evidence','software test','timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
[code,commit]=system('git rev-parse HEAD'); assert(code==0);
[code,branch]=system('git branch --show-current'); assert(code==0);
[code,dirty]=system('git status --short'); assert(code==0);
report.commit=strtrim(commit); report.branch=strtrim(branch); report.dirty=dirty;
files={'simulateEffectiveSocialInput.m','scheduledOcclusionInput.m', ...
    'occlusionInputStep.m','bearingTrackStep.m','defaultEffectiveInputConfig.m'};
for i=1:numel(files)
    report.coreHashes.(matlab.lang.makeValidName(files{i}))= ...
        diagnosticSha256(fullfile(root,'effective-social-input',files{i}));
end
report.checks=struct();
report.checks.P0_1=struct('passed',true,'detail','Actual Git/environment/core SHA256 recorded');
diagnosticJson(outPath,report);
cfg=diagnosticConfig(2026110101,profile); cfg.occlusionDurationSteps=20;
cfg.predictionMode='kalman';
timer=tic;
a=simulateEffectiveSocialInputInstrumented(cfg);
b=simulateEffectiveSocialInputInstrumented(cfg);
assert(isequaln(a,b),'Repeated configuration must be exactly deterministic.');
stock=simulateEffectiveSocialInput(cfg);
stockFields=fieldnames(stock);
for k=1:numel(stockFields)
    key=stockFields{k};
    if strcmp(key,'observationRecords')
        for t=1:cfg.NSteps
            names=fieldnames(stock.observationRecords{t});
            for f=1:numel(names)
                assert(isequaln(stock.observationRecords{t}.(names{f}), ...
                    a.observationRecords{t}.(names{f})));
            end
        end
    else
        assert(isequaln(stock.(key),a.(key)),'Stock mismatch %s',key);
    end
end
clear b stock;
report.checks.P0_2=struct('passed',true,'detail','Exact repeat and all original result fields');
modes={'none','hold','kf-mean','kalman'};
cfg.occlusionDurationSteps=0; cfg.predictionMode='none';
fresh=simulateEffectiveSocialInputInstrumented(cfg);
for k=2:numel(modes)
    cfg.predictionMode=modes{k}; b=simulateEffectiveSocialInputInstrumented(cfg);
    fields={'trajectory','timeSeries','uRecord','positionsRecord','headingsRecord', ...
        'sensorNoise','preOcclusionCheckpoint','observationRecords','finalNeuralState'};
    for f=1:numel(fields), assert(isequaln(fresh.(fields{f}),b.(fields{f}))); end
end
report.checks.P0_3=struct('passed',true,'detail','D0 all four modes, full physical/neural/filter records');
cfg.occlusionDurationSteps=20;
minEigen=inf; maxAsym=0; maxMassError=0;
for k=1:numel(modes)
    cfg.predictionMode=modes{k}; b=simulateEffectiveSocialInputInstrumented(cfg);
    assertDiagnosticPrefix(fresh,b,cfg.occlusionStartStep-1);
    for t=1:cfg.NSteps
        rec=b.observationRecords{t};
        selected=isfinite(rec.trueBearing);
        assert(all(sum(selected,2)==15));
        hidden=t>=cfg.occlusionStartStep && t<cfg.occlusionStartStep+20;
        expectedFresh=15-7*hidden;
        assert(all(sum(rec.fresh,2)==expectedFresh));
        assert(all(isnan(rec.observedBearing(selected & ~rec.fresh))));
        if hidden && strcmp(modes{k},'none'), expectedActive=8; else, expectedActive=15; end
        assert(all(sum(rec.active,2)==expectedActive));
        expectedMass=expectedActive*.012*sqrt(2*pi)*.4;
        massError=abs(b.timeSeries.meanInputMass(t)-expectedMass);
        maxMassError=max(maxMassError,massError);
        assert(massError<1e-10,'Discrete Gaussian mass mismatch');
        p11=rec.rawVariance(selected); p22=rec.rawP22(selected);
        p12=rec.rawP12(selected); p21=rec.rawP21(selected);
        assert(all(isfinite([p11;p22;p12;p21])));
        asym=max(abs(p12-p21)); maxAsym=max(maxAsym,asym);
        eigenMin=.5*(p11+p22-sqrt((p11-p22).^2+4*p12.*p21));
        minEigen=min(minEigen,min(eigenMin));
        assert(asym<1e-12 && min(eigenMin)>=-1e-12);
    end
end
report.checks.P0_4=struct('passed',true,'detail','All four prefixes including last fresh state and complete tracker/sensor/RNG checkpoint');
% Existing independent contracts include hidden-position perturbation, covariance
% convolution mass, wrapped innovations, PSD and separate sensor random stream.
testPredictionModelContract;
angles=repmat((0:39)*2*pi/40,2,1);
obs=struct('selected',logical([0 1;0 0]),'fresh',logical([0 1;0 0]), ...
    'bearing',[NaN .1;NaN NaN],'weight',[NaN 1;NaN NaN]);
opt=struct('dt',.3,'q',1e-4,'gmQ',1e-5,'velocityTau',10,'R',1e-4, ...
    'mode','kalman','maxAge',120,'maxVariance',1,'sigma',.4,'h0',.012);
[~,state]=occlusionInputStepInstrumented(obs,angles,[],opt);
obs.fresh(:)=false; obs.bearing(:)=NaN; obs.weight(:)=NaN;
[input1,state1]=occlusionInputStepInstrumented(obs,angles,state,opt);
obs.trueBearing=123*ones(2); obs.truePosition=-999*ones(2);
[input2,state2]=occlusionInputStepInstrumented(obs,angles,state,opt);
assert(isequaln(input1,input2) && isequaln(state1,state2));
obs.bearing(1,2)=.5; caught=false;
try, occlusionInputStepInstrumented(obs,angles,state,opt);
catch ME, caught=strcmp(ME.identifier,'MATLAB:assertion:failed'); end
assert(caught,'Non-NaN hidden observation must fail its explicit assertion.');
report.checks.P0_5=struct('passed',true,'detail','Hidden truth poison ignored; non-NaN hidden observation rejected; stock hidden-position contract');
report.checks.P0_6=struct('passed',true,'maxMassError',maxMassError,'detail','Selected15 fresh8; expected mass per active edge and heat-convolution contract');
assert(abs(cfg.bearingObservationVariance-cfg.bearingMeasurementStd^2)<eps);
report.checks.P0_7=struct('passed',true,'minimumCovarianceEigenvalue',minEigen, ...
    'maximumCovarianceAsymmetry',maxAsym,'R',cfg.bearingObservationVariance);
% Exercise age-only, variance-only, joint and Hold-ignores-P classifications.
obs=rmfield(obs,{'trueBearing','truePosition'}); obs.bearing(:)=NaN;
state.age(1,2)=5; state.tracks{1,2}.P=eye(2)*10;
opt.maxAge=1; opt.maxVariance=1;
[~,~,d]=occlusionInputStepInstrumented(obs,angles,state,opt); assert(d.cutoffReason(1,2)==23);
opt.mode='hold'; [~,~,d]=occlusionInputStepInstrumented(obs,angles,state,opt);
assert(d.cutoffReason(1,2)==21,'Hold cannot have P cutoff');
opt.mode='kalman'; opt.maxAge=120;
[~,~,d]=occlusionInputStepInstrumented(obs,angles,state,opt); assert(d.cutoffReason(1,2)==22);
opt.maxAge=1; opt.maxVariance=1e9;
[~,~,d]=occlusionInputStepInstrumented(obs,angles,state,opt); assert(d.cutoffReason(1,2)==21);
if strcmp(profile,'p1')
    cfg.occlusionDurationSteps=120; cfg.predictionMode='kalman';
    long=simulateEffectiveSocialInputInstrumented(cfg); counts=zeros(120,1);
    for k=1:120
        rec=long.observationRecords{1400+k}; counts(k)=sum(rec.cutoffReason==22,'all');
        assert(~any(ismember(rec.cutoffReason,[21 23]),'all'));
    end
    assert(all(counts(1:100)==0) && all(counts(101:end)==140));
    report.longPcutoffEdgeSteps=sum(counts);
end
report.checks.P0_8=struct('passed',true,'detail','Age-only/P-only/joint/Hold-ignore-P explicit fixtures; full-profile checks D120 timing');
report.status='PASS'; report.overallPassed=true; report.seconds=toc(timer);
diagnosticJson(outPath,report);
fprintf('PASS REVIEWED P0 %s seconds=%.3f\n',profile,report.seconds);
end
