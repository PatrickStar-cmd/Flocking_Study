function field=buildFreshOracle(bearing,targetIds,alpha,cfg)
%BUILDFRESHORACLE Offline branch-own geometry; never used by the controller.
field=zeros(cfg.NNeurons,cfg.NAgents);
for i=1:cfg.NAgents
    for j=sort(targetIds(i,:))
        d=atan2(sin(alpha(i,:)-bearing(i,j)),cos(alpha(i,:)-bearing(i,j)));
        field(:,i)=field(:,i)+cfg.w0*exp(-0.5*(d/cfg.receptiveFieldWidth).^2).';
    end
end
end
