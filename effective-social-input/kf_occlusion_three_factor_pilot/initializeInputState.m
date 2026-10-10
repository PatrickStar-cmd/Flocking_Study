function state=initializeInputState(cfg,t)
% 所有边的合法缓存。tracks仅从obs更新，不含隐藏真值。
if nargin<2, t=0; end
n=cfg.NAgents;
state.lastBearing=nan(n); state.lastWeight=nan(n); state.age=inf(n);
state.mask=zeros(n,cfg.h); state.maskReady=false; state.maskDiagnostics=cell(n,1);
state.driverIndex=t; state.tracks.initialized=false(n);
for f={'bearing','velocity','P11','P12','P22','nis','priorBearing','priorVariance'}
    state.tracks.(f{1})=nan(n);
end
for r=cfg.sensitivityRs
    state.sensitivityTracks.(sprintf('R_%d',round(-log10(r))))=state.tracks;
end
end
