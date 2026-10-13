function s=checkpointToStart(cp)
assert(isscalar(cp),'SQ1:Checkpoint','Checkpoint must be scalar.');
s.position=cp.position; s.heading=cp.heading; s.u=cp.u; s.alpha=cp.alpha;
s.rng=cp.rng; s.step=cp.step;
s.state.lastBearing=cp.lastBearing; s.state.lastWeight=cp.lastWeight;
s.state.age=cp.age; s.state.mask=cp.mask; s.state.maskReady=cp.maskReady;
s.state.maskDiagnostics=cp.maskDiagnostics; s.state.driverIndex=cp.driverIndex;
assert(isscalar(s.state) && numel(s.state.maskDiagnostics)==cp.config.NAgents);
end
