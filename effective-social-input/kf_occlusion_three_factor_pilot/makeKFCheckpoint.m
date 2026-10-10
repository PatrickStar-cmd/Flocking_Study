function cp=makeKFCheckpoint(pos,head,u,alpha,state,cfg,t)
cp.step=t; cp.position=pos; cp.heading=head; cp.u=u; cp.alpha=alpha;
cp.inputState=state; cp.rng=rng; cp.config=cfg; cp.targetIds=fixedTargetIds(cfg);
end
