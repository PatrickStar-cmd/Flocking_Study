function testLocalSocialInput()
%TESTLOCALSOCIALINPUT Deterministic contracts for the E1 local replay.
cfg=localTestConfig(); alpha=linspace(0,2*pi,cfg.NNeurons+1); alpha=alpha(1:end-1);
[field,m]=localSocialInputMetrics([0 pi/2],[0.25 0.75],alpha,cfg.receptiveFieldWidth);
assert(abs(m.inputMass-sum(field)*2*pi/cfg.NNeurons)<1e-12);
[~,zero]=localSocialInputMetrics(0,0,alpha,cfg.receptiveFieldWidth);
assert(zero.inputMass==0 && ~zero.hasDirection && zero.maxAngularGap==2*pi);

theta=[0 pi/2]; w=[0.25 0.75];
splitTheta=[0 0 pi/2 pi/2]; splitW=[.1 .15 .3 .45];
[f1,a]=localSocialInputMetrics(theta,w,alpha,cfg.receptiveFieldWidth);
[f2,b]=localSocialInputMetrics(splitTheta,splitW,alpha,cfg.receptiveFieldWidth);
assert(max(abs(f1-f2))<1e-12 && abs(a.inputMass-b.inputMass)<1e-12);

s.theta=repmat(theta,cfg.NSteps,1); s.weights=repmat(w,cfg.NSteps,1);
r1=replayLocalSocialInput(cfg,s); r2=replayLocalSocialInput(cfg,s);
assert(isequaln(r1.timeSeries,r2.timeSeries));
assert(max(abs(sum(r1.input,1)-r1.timeSeries.inputMass.'*cfg.NNeurons/(2*pi)))<1e-12);
split.theta=repmat(splitTheta,cfg.NSteps,1); split.weights=repmat(splitW,cfg.NSteps,1);
rs=replayLocalSocialInput(cfg,split);
assert(max(abs(r1.input(:)-rs.input(:)))<1e-12);
assert(max(abs(r1.neuralState(:)-rs.neuralState(:)))<1e-12);
disp('PASS testLocalSocialInput');
end

function cfg=localTestConfig()
cfg.NNeurons=64; cfg.NSteps=12; cfg.dt=.3; cfg.beta=1000;
cfg.globalInhibition=0; cfg.receptiveFieldWidth=.4;
cfg.connectivityExponent=.5; cfg.initialNeuralState=zeros(cfg.NNeurons,1);
end
