function [segment,checkpoints]=resumeProtocolFromCheckpoint(cp,cfg,endStep)
%RESUMEPROTOCOLFROMCHECKPOINT Resume without initialization or index reset.
if nargin<3, endStep=cfg.NSteps; end
[segment,checkpoints]=simulateProtocolSegment(cfg,checkpointToStart(cp),cp.step,endStep);
end
