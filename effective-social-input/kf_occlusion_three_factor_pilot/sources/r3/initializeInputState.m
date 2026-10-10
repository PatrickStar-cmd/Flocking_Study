function state=initializeInputState(cfg,t)
%INITIALIZEINPUTSTATE Empty legal observation cache; no truth is cached.
if nargin<2, t=0; end
n=cfg.NAgents;
state.lastBearing=nan(n,n); state.lastWeight=nan(n,n); state.age=inf(n,n);
state.mask=zeros(n,cfg.h); state.maskReady=false;
state.maskDiagnostics=cell(n,1); state.driverIndex=t;
end
