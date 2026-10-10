function cfg=defaultKFExperimentConfig()
% 三因素开发pilot。R是无噪声观测中的假定方差/数值正则，并非传感器噪声。
cfg=defaultFourBranchConfig();
cfg.protocolVersion='KF-3factor-pilot-1.0';
cfg.seeds=1:5; cfg.checkpointSteps=6200; cfg.recordFromStep=6000;
cfg.q=1e-6; cfg.R=1e-4; cfg.initialVelocityVariance=1;
cfg.sensitivityRs=[1e-6 1e-2];
cfg.KFModel='constant-angular-velocity'; cfg.observationNoiseStd=0;
cfg.parameterSelection='fixed prior exploratory candidate; no new outcome tuning';
cfg.saveFields=true; cfg.collectDiagnostics=true;
cfg.modes={'fresh','drop','drop-normalized','hold','KF-mean'};
cfg.exclusionRule='none; preserve failures'; cfg.workerCount=1;
end
