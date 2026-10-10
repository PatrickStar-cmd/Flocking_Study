function cfg = defaultFourBranchConfig()
%DEFAULTFOURBRANCHCONFIG Protocol 1.2 seed-1 smoke configuration.
cfg.NAgents=20; cfg.NNeurons=100; cfg.NSteps=6600; cfg.dt=0.3;
cfg.arenaSize=1000; cfg.speedScale=0.05; cfg.beta=1000; cfg.globalInhibition=0;
cfg.connectivityExponent=0.5; cfg.receptiveFieldWidth=0.4; cfg.allocentric=true;
cfg.k=8; cfg.neighborCount=8; cfg.h=4; cfg.A0=0.228; cfg.w0=cfg.A0/cfg.k;
cfg.checkpointSteps=[6000 6200]; cfg.maskReadyStep=6200;
cfg.occlusionStartStep=6201; cfg.occlusionDurationSteps=34; cfg.recoveryStartStep=6235;
cfg.seed=1; cfg.maskRule='random'; cfg.mode='fresh'; cfg.recordStride=1;
cfg.initialXY=[]; cfg.initialHeadings=[]; cfg.initialNeuralState=[];
cfg.recordFromStep=0; cfg.protocolVersion='1.2'; cfg.headingThreshold=1e-9;
end
