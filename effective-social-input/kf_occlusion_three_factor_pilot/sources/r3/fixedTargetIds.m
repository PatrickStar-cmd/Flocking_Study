function ids=fixedTargetIds(cfg)
validateattributes(cfg.k,{'numeric'},{'scalar','integer','>=',0,'<=',cfg.NAgents-1});
ids=zeros(cfg.NAgents,cfg.k);
for i=1:cfg.NAgents, ids(i,:)=mod((i-1)+(1:cfg.k),cfg.NAgents)+1; end
end
