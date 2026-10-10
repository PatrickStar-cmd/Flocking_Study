function s=checkpointToStart(cp)
% 完整续跑：神经状态、合法缓存、KF协方差、随机状态、驱动步数一起恢复。
s.position=cp.position; s.heading=cp.heading; s.u=cp.u; s.alpha=cp.alpha;
s.rng=cp.rng; s.step=cp.step; s.state=cp.inputState;
assert(s.state.driverIndex==cp.step && isfield(s.state,'tracks'));
end
