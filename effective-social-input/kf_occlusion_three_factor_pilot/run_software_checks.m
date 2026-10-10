function report=run_software_checks(outDir)
% 数值/因果接口回归。PASS只说明实现检查，不说明科学收益。
root=fileparts(mfilename('fullpath')); addpath(root);
if ~exist(outDir,'dir'), mkdir(outDir); end
cfg=defaultKFExperimentConfig(); cfg.NSteps=90; cfg.recordFromStep=0;
cfg.occlusionStartStep=41; cfg.occlusionDurationSteps=10; cfg.recoveryStartStep=51;
cfg.checkpointSteps=[40 45]; cfg.maskReadyStep=40;
report=struct('test',{},'status',{},'maxDifference',{},'detail',{});
names={'vectorized_KF_matches_scalar','circular_angle_and_PSD','T0_five_methods',...
    'h0_five_methods','checkpoint_resume','hidden_truth_interface','weight_and_recovery_contract',...
    'h8_no_fresh_boundary','A0_zero_boundary','legacy_r3_core_parity'};
for id=1:numel(names)
    diff=0; detail='Executed assertions';
    try
        switch id
            case 1
                c=cfg; c.NAgents=2; c.k=1; st=initializeInputState(c); tr=[];
                for t=1:120
                    o.selected=[false true;true false]; o.fresh=o.selected;
                    o.bearing=nan(2); o.weight=nan(2); obs=mod(6+0.02*t,2*pi);
                    if t>=60 && t<=90, o.fresh(1,2)=false; else, o.bearing(1,2)=obs; end
                    o.bearing(2,1)=0.7; o.weight(o.fresh)=c.w0;
                    st.tracks=bearingTrackBankStep(st.tracks,o,c);
                    if o.fresh(1,2), z=obs; else, z=[]; end
                    tr=scalarKF(tr,z,c); bank=st.tracks;
                    delta=max(abs([bank.bearing(1,2);bank.velocity(1,2)]-tr.x));
                    delta=max(delta,max(abs([bank.P11(1,2) bank.P12(1,2);bank.P12(1,2) bank.P22(1,2)]-tr.P),[],'all'));
                    diff=max(diff,delta); assert(diff<1e-11);
                end
            case 2
                c=cfg; c.NAgents=2; c.k=1; st=initializeInputState(c);
                for t=1:200
                    o.selected=[false true;true false]; o.fresh=o.selected;
                    o.bearing=[nan mod(pi-0.1+0.03*t,2*pi);0.7 nan]; o.weight=c.w0*double(o.selected);
                    st.tracks=bearingTrackBankStep(st.tracks,o,c); b=st.tracks;
                    assert(b.P11(1,2)>=0 && b.P22(1,2)>=0);
                    assert(b.P11(1,2)*b.P22(1,2)-b.P12(1,2)^2>=-1e-13);
                    if t>30, assert(b.nis(1,2)<1e-3); end
                end
            case {3,4}
                c=cfg; if id==3, c.occlusionDurationSteps=0; c.recoveryStartStep=41; else, c.h=0; end
                for k=1:5
                    c.mode=c.modes{k}; x=simulateKFSegment(c,initializeCoreState(c),0,90);
                    if k==1, ref=x; else, assert(isequaln(ref.state,x.state)); assert(isequaln(ref.data,x.data)); end
                end
            case 5
                c=cfg; c.mode='KF-mean'; a=simulateKFSegment(c,initializeCoreState(c),0,90);
                cp=a.checkpoints.step_45; b=simulateKFSegment(c,checkpointToStart(cp),45,90);
                assert(isequaln(a.data.inputField(:,:,46:end),b.data.inputField));
                assert(isequaln(a.state.position(:,:,46:end),b.state.position));
                assert(isequaln(a.state.u(:,:,46:end),b.state.u));
                assert(isequaln(a.finalState,b.finalState) && isequaln(a.finalRng,b.finalRng));
            case 6
                c=cfg; c.mode='KF-mean'; a=simulateKFSegment(c,initializeCoreState(c),0,40);
                st=a.finalState; [o,truth]=sensorObservation(a.finalPosition,st,c,41);
                [u1,s1]=controllerInputFromObservation(a.finalAlpha,st,o,c,41);
                truth.trueBearing(:)=1e99; %#ok<NASGU> evaluation truth is not an argument
                [u2,s2]=controllerInputFromObservation(a.finalAlpha,st,o,c,41);
                assert(isequaln(u1,u2) && isequaln(s1,s2));
                bad=o; bad.bearing(o.selected & ~o.fresh)=0;
                failed=false; try, controllerInputFromObservation(a.finalAlpha,st,bad,c,41); catch, failed=true; end
                assert(failed);
                source=fileread(fullfile(root,'controllerInputFromObservation.m'));
                assert(~contains(source,'trueBearing') && ~contains(source,'trueDistance'));
                detail='Controller has no truth/position argument; rejects finite hidden measurements';
            case 7
                for k=1:5
                    c=cfg; c.mode=c.modes{k}; a=simulateKFSegment(c,initializeCoreState(c),0,60); d=a.data;
                    oc=d.step>=41 & d.step<=50; rec=d.step==51;
                    if strcmp(c.mode,'drop'), expected=c.A0*(1-c.h/c.k); else, expected=c.A0; end
                    assert(max(abs(d.A(oc,:)-expected),[],'all')<1e-14);
                    assert(all(d.k_active(rec,:)==8)); assert(max(abs(d.A(rec,:)-c.A0))<1e-14);
                    assert(max(d.selfFreshRelativeL2(rec,:))<1e-14);
                    h=d.hidden(:,:,41:50); observed=d.observedBearing(:,:,41:50);
                    assert(all(isnan(observed(h))));
                end
            case 8
                for mode={'drop','drop-normalized','KF-mean'}
                    c=cfg; c.h=8; c.mode=mode{1}; a=simulateKFSegment(c,initializeCoreState(c),0,60);
                    if strcmp(mode{1},'KF-mean'), assert(all(a.data.k_active(41:50,:)==8,'all'));
                    else, assert(all(a.data.A(41:50,:)==0,'all')); end
                    assert(all(isfinite(a.state.position),'all'));
                end
            case 9
                c=cfg; c.A0=0; c.w0=0; c.mode='KF-mean'; a=simulateKFSegment(c,initializeCoreState(c),0,60);
                assert(all(a.data.A==0,'all')); assert(all(isnan(a.data.selfFreshRelativeL2),'all'));
            case 10
                originalPath=path; pathCleanup=onCleanup(@()path(originalPath));
                addpath(fullfile(root,'sources','r3'),'-begin'); c=cfg; c.mode='fresh';
                old=simulateProtocolSegment(c,initializeCoreState(c),0,90);
                rmpath(fullfile(root,'sources','r3')); addpath(root,'-begin');
                new=simulateKFSegment(c,initializeCoreState(c),0,90);
                diff=max(abs(old.state.position-new.state.position),[],'all');
                diff=max(diff,max(abs(old.state.u-new.state.u),[],'all'));
                assert(diff<1e-10); clear pathCleanup;
        end
        status='PASS';
    catch ME
        status='FAIL'; detail=getReport(ME,'extended','hyperlinks','off');
    end
    report(id)=struct('test',names{id},'status',status,'maxDifference',diff,'detail',detail);
    fprintf('%s: %s (max diff %.3g)\n',names{id},status,diff);
end
writetable(struct2table(report),fullfile(outDir,'software_checks.csv'));
writeJson(fullfile(outDir,'software_checks.json'),report);
assert(all(strcmp({report.status},'PASS')),'KF:SoftwareGate','Software checks failed; inspect preserved log');
end
function tr=scalarKF(tr,observation,c)
F=[1 c.dt;0 1]; Q=c.q*[c.dt^3/3 c.dt^2/2;c.dt^2/2 c.dt]; R=c.R;
if isempty(tr), tr.x=[observation;0]; tr.P=diag([R c.initialVelocityVariance]); return; end
tr.x=F*tr.x; tr.P=F*tr.P*F'+Q;
if ~isempty(observation)
    v=atan2(sin(observation-tr.x(1)),cos(observation-tr.x(1))); k=tr.P(:,1)/(tr.P(1,1)+R);
    tr.x=tr.x+k*v; A=eye(2)-k*[1 0]; tr.P=A*tr.P*A'+k*R*k';
end
tr.x(1)=mod(tr.x(1),2*pi); tr.P=(tr.P+tr.P')/2;
end
