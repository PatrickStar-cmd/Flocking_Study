function seg=run_closed_loop_experiment(cp,cfg)
% 每种方法从完全相同的6200完整checkpoint分叉。
seg=simulateKFSegment(cfg,checkpointToStart(cp),cp.step,cfg.NSteps);
end
