function result=protocolAssertions(id,outDir)
%PROTOCOLASSERTIONS Actual protocol gates; all thresholds are set before runs.
cfg=defaultFourBranchConfig(); root=fileparts(fileparts(mfilename('fullpath')));
result.maxDiff=0; result.count=0; result.tolerance=0; result.details='';
modes={'fresh','drop','drop-normalized','hold'};
switch id
case 1
    % Equal all-neighbor topology is required for an honest original-core comparison.
    addpath(fullfile(root,'sources','legacy'),'-end');
    c=cfg; c.NAgents=5; c.NNeurons=24; c.k=4; c.h=0; c.w0=.228/4;
    c.occlusionStartStep=100; c.checkpointSteps=9; c.recordFromStep=0;
    maximum=0;
    for inhibition=[0 .02]
        c.globalInhibition=inhibition; s=initializeCoreState(c);
        rng(c.seed,'twister'); rand(2,c.NAgents);
        W=buildRingConnectivity(c.NNeurons,c.connectivityExponent);
        [x,y,u,h]=simulateCollectiveFlockingNeuralField(c.NAgents,c.NNeurons,9,c.dt,c.speedScale,c.arenaSize,W,inhibition,c.beta,c.w0,c.receptiveFieldWidth,s.position,true);
        a=simulateProtocolSegment(c,s,0,9);
        diffs=[max(abs(squeeze(a.state.position(1,:,:))-x),[],'all'),...
            max(abs(squeeze(a.state.position(2,:,:))-y),[],'all'),...
            max(abs(a.state.u-u),[],'all'),max(abs(a.state.heading-h),[],'all')];
        maximum=max(maximum,max(diffs));
    end
    result.maxDiff=maximum; result.tolerance=1e-12; assert(maximum<=result.tolerance);
    result.count=2*9; result.details='Original unmodified author core, equivalent all-neighbor topology, two inhibition values; per-step position/u/heading. Heading threshold restored to original 1e-9; no fixed-k versus old dynamic topology comparison.';
case 2
    targets=fixedTargetIds(cfg);
    for i=1:cfg.NAgents
        assert(numel(unique(targets(i,:)))==8 && ~ismember(i,targets(i,:)));
    end
    for q=1:4
        r=loadRun(outDir,modes{q}); a=r.segment.data;
        for t=1:600
            for i=1:20, assert(isequal(find(a.selected(i,:,t)),sort(targets(i,:)))); end
        end
    end
    result.count=48000; result.details='Every saved step/branch uses the same eight unique cyclic IDs and no self edge.';
case 3
    x=load(fullfile(outDir,'T0_results.mat'),'t0'); a=x.t0.fresh;
    for q=2:4
        b=x.t0.(matlab.lang.makeValidName(modes{q})); z=compareSegments(a,b,6001,6600);
        assert(z.mismatchCount==0); assert(~b.finalState.maskReady);
    end
    writetable(table(string(modes(2:4)).',zeros(3,1),zeros(3,1),'VariableNames',{'mode','maxDiff','mismatchCount'}),fullfile(outDir,'T0_maxdiff.csv'));
    result.count=3*600; result.details='Exact equality of all exported steps, 601 frames including u/alpha, final cache and RNG; d0 does not enable a mask.';
case 4
    rows=struct([]); count=0;
    for q=1:4
        r=loadRun(outDir,modes{q});
        if q==1
            x=load(fullfile(outDir,'continuous_reference.mat'),'continuous'); a=x.continuous;
        else
            c=r.config; c.recordFromStep=6000;
            a=simulateProtocolSegment(c,initializeCoreState(c),0,6600);
        end
        for step=[6000 6200]
            cp=a.checkpoints.(['step_' num2str(step)]);
            assert(isscalar(cp)&&numel(cp.maskDiagnostics)==20);
            b=resumeProtocolFromCheckpoint(cp,r.config,6600);
            z=compareSegments(a,b,step+1,6600);
            assert(z.mismatchCount==0);
            count=count+1; rows(count).mode=modes{q}; rows(count).checkpointStep=step;
            rows(count).maxDiff=z.maxDiff; rows(count).mismatchCount=z.mismatchCount;
            rows(count).fieldsCompared=z.fields; rows(count).maskDiagnosticsCount=numel(cp.maskDiagnostics);
        end
        z=compareSegments(a,r.segment,6001,6600); assert(z.mismatchCount==0);
        fprintf('G4 %s: uninterrupted6600/cp6000/cp6200/saved-run exact equality\n',modes{q});
    end
    writetable(struct2table(rows),fullfile(outDir,'resume_check.csv'));
    result.count=8; result.details='All four modes run uninterrupted from initialization to 6600 in this MATLAB session; both scalar checkpoints resume exactly; every data field, per-step u, static alpha, cache/age, final RNG and driverIndex are compared.';
case 5
    a=loadRun(outDir,'fresh');
    for q=2:4
        b=loadRun(outDir,modes{q}); z=compareSegments(a.segment,b.segment,6001,6200);
        assert(z.mismatchCount==0); assert(isequaln(a.checkpoint6200.mask,b.checkpoint6200.mask));
        assert(b.checkpoint6200.maskReady);
    end
    result.count=600; result.details='All saved common-prefix fields including u and cache are exactly equal; mask is frozen using step6200-start geometry.';
case 6
    expected=[.228 .114 .228 .228]; rows=struct([]);
    for q=1:4
        r=loadRun(outDir,modes{q}); ix=r.segment.data.phase=="occlusion";
        e=max(abs(r.segment.data.A(ix,:)-expected(q)),[],'all'); assert(e<1e-12);
        rows(q).mode=modes{q}; rows(q).expectedA=expected(q); rows(q).maxDiff=e;
    end
    writetable(struct2table(rows),fullfile(outDir,'weight_contract.csv'));
    count=0;
    for k=[0 1 4 8]
        for h=unique([0 floor(k/2) k])
            c=cfg; c.k=k; c.h=h; c.w0=.228/max(1,k);
            s=initializeCoreState(c); [mask,~]=makeProtocolMasks(s.position,c); assert(isequal(size(mask),[20 h]));
            for q=1:4
                c.mode=modes{q}; [o,cache]=syntheticObservation(c,zeros(1,k),repmat(c.w0,1,k));
                hidden=1:h; ids=find(o.selected(1,:));
                for i=1:20
                    ids=find(o.selected(i,:)); o.fresh(i,ids(hidden))=false;
                    o.bearing(i,ids(hidden))=nan; o.weight(i,ids(hidden))=nan;
                end
                [field,~,d]=controllerInputFromObservation(s.alpha,cache,o,c,6201);
                if strcmp(c.mode,'fresh') || strcmp(c.mode,'drop'), want=(k-h)*c.w0;
                elseif strcmp(c.mode,'drop-normalized'), want=c.A0*(k>h);
                else, want=k*c.w0; end
                assert(max(abs(d.A-want))<1e-12 && all(isfinite(field(:))));
                if want==0, assert(all(d.no_input) && all(d.k_active==0) && all(field(:)==0)); end
                if h==k && h>0 && strcmp(c.mode,'drop-normalized'), assert(all(d.normalization_unavailable)); end
                count=count+1;
            end
        end
    end
    c=cfg; c.mode='hold'; [o,~]=syntheticObservation(c,zeros(1,8),ones(1,8));
    o.fresh(:)=false; o.bearing(:)=nan; o.weight(:)=nan;
    [f,~,d]=controllerInputFromObservation(s.alpha,[],o,c,1); assert(all(f(:)==0)&&all(d.k_active==0));
    for q=1:4
        c.mode=modes{q}; [o,cache]=syntheticObservation(c,zeros(1,8),zeros(1,8));
        for i=1:20
            ids=find(o.selected(i,:)); o.fresh(i,ids(1:4))=false;
            o.bearing(i,ids(1:4))=nan; o.weight(i,ids(1:4))=nan;
        end
        [f,~,d]=controllerInputFromObservation(s.alpha,cache,o,c,1);
        assert(all(f(:)==0)&&all(d.no_input)&&all(d.k_active==0));
    end
    result.count=count+5; result.tolerance=1e-12; result.details='Main per-agent amplitude and k0/k1/k4/k8, h0/hk/intermediate, all missing, zero weights in all four modes with missing edges, and hold with no cache.';
case 7
    count=0;
    for duration=[0 1 10 34 100]
        c=cfg; c.NAgents=8; c.NNeurons=40; c.k=4; c.h=2; c.w0=c.A0/c.k;
        c.occlusionStartStep=4; c.maskReadyStep=3; c.occlusionDurationSteps=duration;
        c.checkpointSteps=3; c.mode='drop-normalized'; c.recordFromStep=0;
        a=simulateProtocolSegment(c,initializeCoreState(c),0,5+duration);
        assert(sum(a.data.phase=="occlusion")==duration);
        if duration>0, assert(all(a.data.k_fresh(4:3+duration,:)==c.k-c.h,'all')); end
        assert(a.data.phase(4+duration)=="recovery");
        assert(all(a.data.k_fresh(4+duration,:)==c.k));
        w=a.data.edgeWeight(:,:,4+duration); assert(all(w(w>0)==c.w0));
        if duration==0, assert(~a.finalState.maskReady); end
        count=count+1;
    end
    c=cfg; c.occlusionDurationSteps=10; c.mode='drop';
    x=load(fullfile(outDir,'checkpoints.mat'),'cp6000'); a=resumeProtocolFromCheckpoint(x.cp6000,c,6212);
    assert(a.data.phase(a.data.step==6211)=="recovery");
    result.count=count+1; result.details='d0/1/10/34/100 input and dynamic phase boundaries; explicit absolute6211 d10 recovery and w0 restoration.';
case 8
    r=loadRun(outDir,'hold'); a=r.segment.data; mask=r.checkpoint6200.mask;
    count=0;
    for i=1:20
        for j=mask(i,:)
            for t=6201:6234
                q=t-6000;
                assert(a.lastBearing(i,j,q)==r.checkpoint6200.lastBearing(i,j));
                assert(a.injectedBearing(i,j,q)==r.checkpoint6200.lastBearing(i,j));
                assert(abs(a.age(i,j,q)-(t-6200)*cfg.dt)<1e-12); count=count+1;
            end
            assert(a.age(i,j,235)==0); assert(a.fresh(i,j,235));
        end
    end
    result.count=count; result.tolerance=1e-12; result.details='Every hidden hold edge retains last legal step6200-start bearing, age0.3 through10.2, no expiry, recovery resets to0.';
case 9
    s=initializeCoreState(cfg); prior=rng; [mask,d]=makeProtocolMasks(s.position,cfg); assert(isequal(prior,rng));
    [again,~]=makeProtocolMasks(s.position,cfg); assert(isequal(mask,again));
    stream=RandStream('mt19937ar','Seed',100001);
    for i=1:20
        pick=randi(stream,70); assert(isequal(mask(i,:),d{i}.candidateIds(pick,:)));
        assert(size(unique(d{i}.candidateIds,'rows'),1)==70);
    end
    uniformStream=RandStream('mt19937ar','Seed',73991); histogram=accumarray(randi(uniformStream,70,[70000 1]),1,[70 1]);
    assert(all(abs(histogram-1000)<150));
    writeJson(fullfile(outDir,'logs','mask_uniformity.json'),struct('seed',73991,'draws',70000,'counts',histogram,'expected',1000,'predeclaredMaxDeviation',150));
    for rule={'max_hidden_R1','min_hidden_R1'}
        c=cfg; c.maskRule=rule{1}; [~,d]=makeProtocolMasks(s.position,c);
        for i=1:20
            if strcmp(c.maskRule,'max_hidden_R1'), best=max(d{i}.R1); else, best=min(d{i}.R1); end
            assert(abs(d{i}.R1(d{i}.chosenIndex)-best)<=1e-12);
        end
        b=rand(20,20)*2*pi; b(1,2)=pi-1e-9; b(1,3)=-pi+1e-9;
        [a,~]=makeProtocolMasks(s.position,c,b); [rot,~]=makeProtocolMasks(s.position,c,b+.773); assert(isequal(a,rot));
        [wrap,~]=makeProtocolMasks(s.position,c,mod(b+pi,2*pi)-pi); assert(isequal(a,wrap));
        [tie,~]=makeProtocolMasks(s.position,c,zeros(20)); assert(isequal(tie(17,:),[1 2 3 4]));
    end
    result.count=70000+20*3; result.tolerance=1e-12; result.details='70 unique lexicographic actual-ID tuples, independent reproducible MT19937 uniform candidate draws, fixed-seed uniformity check, extrema, actualID tie, rotation and +/-pi wrap invariance; global RNG unaffected by masks.';
case 10
    s=initializeCoreState(cfg); maximum=0;
    for q=2:4
        c=cfg; c.mode=modes{q}; cache=initializeInputState(c); [cache.mask,cache.maskDiagnostics]=makeProtocolMasks(s.position,c); cache.maskReady=true;
        [legal,~]=sensorObservation(s.position,cache,c,6200);
        [~,cache]=controllerInputFromObservation(s.alpha,cache,legal,c,6200);
        [o,truth]=sensorObservation(s.position,cache,c,6201);
        truthPerturbed=truth; truthPerturbed.trueBearing(o.selected&~o.fresh)=17;
        [f,a,da]=controllerInputFromObservation(s.alpha,cache,o,c,6201);
        [g,b,db]=controllerInputFromObservation(s.alpha,cache,o,c,6201);
        assert(~isequaln(truth,truthPerturbed)); assert(isequaln(f,g)&&isequaln(a,b)&&isequaln(da,db));
        maximum=max(maximum,max(abs(f-g),[],'all'));
    end
    ids=fixedTargetIds(cfg); oracle1=buildFreshOracle(truth.trueBearing,ids,s.alpha,cfg);
    oracle2=buildFreshOracle(truthPerturbed.trueBearing,ids,s.alpha,cfg); assert(max(abs(oracle1-oracle2),[],'all')>1e-4);
    result.maxDiff=maximum; result.count=3; result.details='Fixed legal observations/cache/alpha; hidden diagnostic truth changed without controller access. All nonfresh output/cache fields exactly unchanged; offline truth oracle changes. Public signature excludes position/truth.';
case 11
    c=cfg; c.k=4; c.h=0; s=initializeCoreState(c);
    thetaCases={zeros(1,4),[0 0 pi pi],[0 pi/2 pi 3*pi/2]};
    expectedR1=[1 0 0]; expectedR2=[1 1 0]; expectedGap=[2*pi pi pi/2];
    for q=1:3
        [o,cache]=syntheticObservation(c,thetaCases{q},ones(1,4));
        [~,~,d]=controllerInputFromObservation(s.alpha,cache,o,c,1);
        assert(abs(d.z1_mod(1)-expectedR1(q))<1e-12 && abs(d.z2_mod(1)-expectedR2(q))<1e-12);
        assert(abs(d.maxGap(1)-expectedGap(q))<1e-12); assert(d.directionValid(1)==(q==1));
        assert(abs(d.weightEntropy(1)-log(4))<1e-12);
    end
    c.k=8; [o,cache]=syntheticObservation(c,zeros(1,8),ones(1,8));
    [~,~,d]=controllerInputFromObservation(s.alpha,cache,o,c,1); assert(abs(d.weightEntropy(1)-log(8))<1e-12);
    c.k=0; c.h=0; [o,cache]=syntheticObservation(c,[],[]);
    [f,~,d]=controllerInputFromObservation(s.alpha,cache,o,c,1);
    assert(all(f(:)==0)&&all(d.A==0)&&all(d.z1_mod==0)&&all(d.z2_mod==0)&&all(d.no_input));
    assert(all(isnan(d.direction))&&all(isnan(d.maxGap))&&all(isnan(d.weightEntropy)));
    result.count=5; result.tolerance=1e-12; result.details='Single direction, opposing double peak, uniform angles and empty input; R1/R2/gap/direction flags; ln4/ln8 unnormalized label weight entropy.';
case 12
    x=load(fullfile(outDir,'reference_same_state_replay.mat'),'replay'); r=x.replay;
    maximum=0;
    for q=1:4
        f=r.candidates.(matlab.lang.makeValidName(modes{q})).field;
        maximum=max(maximum,max(abs(f(:,:,1:200)-r.candidates.fresh.field(:,:,1:200)),[],'all'));
    end
    assert(maximum==0);
    ix=r.step>=6201 & r.step<=6234;
    assert(max(abs(r.candidates.drop_normalized.A(ix,:)-r.candidates.fresh.A(ix,:)),[],'all')<1e-12);
    assert(max(abs(r.candidates.drop_normalized.field(:,:,ix)-r.candidates.fresh.field(:,:,ix)),[],'all')>1e-4);
    fresh=loadRun(outDir,'fresh'); c=cfg; c.h=0;
    f=struct(); f.config=c; f.segment=fresh.segment; f.checkpoint6000=fresh.checkpoint6000;
    f.checkpoint6200=fresh.checkpoint6200; z=buildReferenceReplay(f);
    for q=2:4, assert(isequaln(z.candidates.fresh.field,z.candidates.(matlab.lang.makeValidName(modes{q})).field)); end
    for q=1:4
        a=loadRun(outDir,modes{q});
        for j=[1 200 201 234 235 600]
            oracle=buildFreshOracle(a.segment.data.trueBearing(:,:,j),fixedTargetIds(cfg),a.segment.state.alpha,cfg);
            assert(isequal(oracle,a.segment.data.selfFreshField(:,:,j)));
        end
    end
    result.count=4*6+3; result.details='Two distinct oracle types, common-state prefix/h0 equality, normalized A equality with unequal field structure, and branch-own oracle reconstruction.';
case 13
    actual=readtable(fullfile(outDir,'summary.csv')); expected=generateSummary(outDir);
    maximum=0; vars=expected.Properties.VariableNames;
    for name=vars
        if isnumeric(expected.(name{1})), maximum=max(maximum,max(abs(actual.(name{1})-expected.(name{1})),[],'all')); end
    end
    assert(maximum<1e-12); manifest=readtable(fullfile(outDir,'RUN_MANIFEST.csv'));
    assert(height(manifest)==4&&numel(unique(manifest.run_id))==4&&all(string(manifest.status)=="COMPLETED"));
    source=jsondecode(fileread(fullfile(outDir,'source_manifest.json')));
    for row=source.sourceFiles.'
        assert(strcmp(row.SHA256,fileSHA256(fullfile(root,strrep(row.relativePath,'/',filesep)))));
    end
    for q=1:4
        info=whos('-file',fullfile(outDir,[modes{q} '_seed1.mat'])); assert(numel(info)==1 && strcmp(info.name,'run'));
        a=loadRun(outDir,modes{q}); assert(numel(a.segment.data.step)==600 && size(a.segment.state.position,3)==601);
        assert(all(isfield(a.segment.data,{'lastBearing','lastWeight','observedWeight','neuralActivityMass','decoded_direction_defined','selfFreshRelativeL2','hiddenR2'})));
        assert(all(isnan(a.segment.data.observedBearing(a.segment.data.hidden))));
    end
    % Independent periodic distance, velocity GO/LO and input formula reconstruction.
    for q=1:4
        a=loadRun(outDir,modes{q});
        for t=1:600
            p=a.segment.state.position(:,:,t+1); prev=a.segment.state.position(:,:,t);
            v=p-prev; v=v-cfg.arenaSize*round(v/cfg.arenaSize); speed=hypot(v(1,:),v(2,:));
            unit=zeros(2,20); valid=speed>cfg.dt*1e-12; unit(:,valid)=v(:,valid)./speed(valid);
            go=norm(sum(unit,2))/20; distances=zeros(20);
            for i=1:20
                for j=1:20
                    delta=p(:,j)-p(:,i); delta=delta-cfg.arenaSize*round(delta/cfg.arenaSize); distances(i,j)=norm(delta);
                end
            end
            pair=mean(distances(triu(true(20),1)))/cfg.arenaSize;
            lo=zeros(20,1);
            for i=1:20, [~,ord]=sort(distances(i,:)); lo(i)=norm(sum(unit(:,ord(1:6)),2))/6; end
            maximum=max([maximum abs(go-a.segment.data.GO(t)) abs(mean(lo)-a.segment.data.LO(t)) abs(pair-a.segment.data.pairDistance(t))]);
            assert(all(abs(sum(a.segment.data.edgeWeight(:,:,t),2)-a.segment.data.A(t,:).')<1e-12));
            if ismember(t,[1 200 201 234 235 600])
                for i=1:20
                    active=a.segment.data.active(i,:,t);
                    th=a.segment.data.injectedBearing(i,active,t);
                    weight=a.segment.data.edgeWeight(i,active,t);
                    [field,metrics]=localSocialInputMetrics(th,weight,a.segment.state.alpha(i,:),cfg.receptiveFieldWidth);
                    maximum=max([maximum max(abs(field-a.segment.data.inputField(:,i,t))),...
                        abs(metrics.inputMass-a.segment.data.Mass(t,i)),...
                        abs(real(metrics.z1)-a.segment.data.z1_real(t,i)),...
                        abs(imag(metrics.z2)-a.segment.data.z2_imag(t,i))]);
                end
            end
        end
    end
    assert(maximum<=1e-12);
    result.maxDiff=maximum; result.count=2400; result.tolerance=1e-12;
    result.details='Saved single-run MAT/keys/600steps/601frames, source hashes, missing fields and hidden NaN checked; CSV regenerated; all2400 saved steps independently recompute periodic pair distance, GO, nearest5+self LO and A; 480 agent-step cases independently reconstruct Gaussian field, Mass and moments.';
case 14
    assert(strcmp(which('simulateProtocolSegment'),fullfile(root,'simulateProtocolSegment.m')));
    before=fileSHA256(fullfile(outDir,'fresh_seed1.mat')); caught=false;
    try, run_protocol12_seed1_r3(lastFolder(outDir)); catch ME, caught=strcmp(ME.identifier,'SQ1:Exists'); end
    assert(caught&&strcmp(before,fileSHA256(fullfile(outDir,'fresh_seed1.mat'))));
    oldpath=path; cleanup=onCleanup(@()path(oldpath)); restoredefaultpath; addpath(root);
    c=defaultFourBranchConfig(); c.NAgents=8; c.NNeurons=40; c.k=4; c.h=2; c.w0=c.A0/4;
    c.occlusionStartStep=4; c.checkpointSteps=3; c.mode='drop';
    s=initializeCoreState(c); a=simulateProtocolSegment(c,s,0,8); assert(all(isfinite(a.state.u(:))));
    % Real deliberately invalid input failure is retained as negative-test evidence.
    c.maskRule='INVALID_TEST_RULE'; caught=false;
    try, makeProtocolMasks(s.position,c); catch ME
        caught=true; f=fopen(fullfile(outDir,'failure_logs','expected_invalid_mask.log'),'w');
        fprintf(f,'%s\n',getReport(ME,'extended','hyperlinks','off')); fclose(f);
    end
    assert(caught); result.count=3;
    result.details='Clean MATLAB search path runs short end-to-end; existing output rejected before writes, MAT hash unchanged; deliberate invalid mask failure retained in failure_logs (expected negative test).';
case 15
    addpath(fullfile(root,'sources','A'),'-end'); c=cfg; c.k=4; c.h=0; s=initializeCoreState(c); maximum=0;
    for shift=[0 .71 -2.83]
        grid=mod(s.alpha+shift,2*pi); theta=[0 pi-1e-10 -pi+1e-10 1.3]; weights=[.01 .03 .07 .04];
        [o,cache]=syntheticObservation(c,theta,weights); [f,~,d]=controllerInputFromObservation(grid,cache,o,c,1);
        for i=1:20
            [ref,m]=localSocialInputMetrics(theta,weights,grid(i,:),c.receptiveFieldWidth);
            scalar=zeros(100,1);
            for neuron=1:100
                for edge=1:4
                    delta=atan2(sin(grid(i,neuron)-theta(edge)),cos(grid(i,neuron)-theta(edge)));
                    scalar(neuron)=scalar(neuron)+weights(edge)*exp(-delta^2/(2*c.receptiveFieldWidth^2));
                end
            end
            maximum=max([maximum max(abs(f(:,i)-ref)) max(abs(f(:,i)-scalar)) abs(d.weightEntropy(i)-m.angularEntropy)]);
        end
    end
    one=c; one.k=1; [o,a]=syntheticObservation(one,.8,.1); [f,~,~]=controllerInputFromObservation(s.alpha,a,o,one,1);
    two=c; two.k=2; [o,a]=syntheticObservation(two,[.8 .8],[.05 .05]); [g,~,~]=controllerInputFromObservation(s.alpha,a,o,two,1);
    maximum=max(maximum,max(abs(f-g),[],'all')); W=buildRingConnectivity(100,.5);
    uf=s.u+c.dt*(-s.u+W*tanh(c.beta*s.u)/100+f);
    ug=s.u+c.dt*(-s.u+W*tanh(c.beta*s.u)/100+g); maximum=max(maximum,max(abs(uf-ug),[],'all'));
    y=ones(100,20); nd=decodeNeuralActivity(y,s.alpha);
    assert(all(~nd.decoded_direction_defined)&&all(isnan(nd.decodedAnalysisAngle))&&all(isfinite(nd.rawDecodedAngle)));
    nd0=decodeNeuralActivity(zeros(100,20),s.alpha); assert(all(~nd0.decoded_direction_defined)&&all(nd0.neuralActivityMass==0));
    assert(maximum<=1e-12); result.maxDiff=maximum; result.tolerance=1e-12; result.count=63;
    result.details='Pinned A@992e71a and independent scalar Gaussian formula on actual shifted alpha grids and +/-pi; label split field/neural equivalence; uniform and zero neural activity flags. Empty/zero/hidden/no-cache and ln4/ln8 cases are also asserted in G6/G10/G11, invoked by the same suite.';
end
fprintf('%s\nmaxDiff=%.17g count=%d tolerance=%.17g\n',result.details,result.maxDiff,result.count,result.tolerance);
end
function r=loadRun(outDir,mode)
x=load(fullfile(outDir,[mode '_seed1.mat']),'run'); r=x.run;
end
function [obs,state]=syntheticObservation(c,theta,weights)
n=c.NAgents; ids=fixedTargetIds(c); obs.selected=false(n); obs.fresh=false(n);
obs.bearing=nan(n); obs.weight=nan(n); state=initializeInputState(c);
for i=1:n
    obs.selected(i,ids(i,:))=true; obs.fresh(i,ids(i,:))=true;
    obs.bearing(i,ids(i,:))=theta; obs.weight(i,ids(i,:))=weights;
    state.lastBearing(i,ids(i,:))=theta; state.lastWeight(i,ids(i,:))=weights; state.age(i,ids(i,:))=0;
end
end
function name=lastFolder(path)
[~,name]=fileparts(path);
end
