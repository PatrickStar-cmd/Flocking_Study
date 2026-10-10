function start=initializeCoreState(cfg)
%INITIALIZECORESTATE Match legacy per-agent random draw order.
rng(cfg.seed,'twister'); n=cfg.NAgents; m=cfg.NNeurons; L=cfg.arenaSize;
position=rand(2,n)*L; heading=zeros(n,1); u=zeros(m,n);
for i=1:n, heading(i)=2*pi*rand; u(:,i)=0.1*randn(m,1); end
alpha0=linspace(0,2*pi,m+1); alpha0=alpha0(1:end-1); alpha=mod(repmat(alpha0,n,1)+heading,2*pi);
start=struct('position',position,'heading',heading,'u',u,'alpha',alpha,'state',[], 'rng',rng,'step',0);
end
