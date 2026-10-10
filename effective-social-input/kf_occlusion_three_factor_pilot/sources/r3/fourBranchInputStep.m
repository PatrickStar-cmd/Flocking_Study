function [input,state,diag]=fourBranchInputStep(position,alpha,state,cfg,t)
%FOURBRANCHINPUTSTEP Sensor, legal cache and separate offline truth channel.
if isempty(state), state=initializeInputState(cfg,t-1); end
if t==cfg.occlusionStartStep-1 && ~state.maskReady && cfg.occlusionDurationSteps>0
    [state.mask,state.maskDiagnostics]=makeProtocolMasks(position,cfg);
    state.maskReady=true;
end
[obs,truth,ids]=sensorObservation(position,state,cfg,t);
[input,state,diag]=controllerInputFromObservation(alpha,state,obs,cfg,t);
diag.selected=obs.selected; diag.fresh=obs.fresh; diag.hidden=obs.selected & ~obs.fresh;
diag.observedBearing=obs.bearing; diag.observedWeight=obs.weight;
diag.trueBearing=truth.trueBearing; diag.trueDistance=truth.trueDistance;
diag.targetIds=ids; diag.mask=state.mask; diag.driverIndex=t;
end
