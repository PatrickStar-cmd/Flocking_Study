function [input,state,d]=fourBranchInputStep(position,alpha,state,cfg,t)
% 传感器 -> 合法KF/缓存 -> 神经输入；最后才添加评价用真值。
if isempty(state), state=initializeInputState(cfg,t-1); end
if t==cfg.occlusionStartStep-1 && ~state.maskReady && cfg.occlusionDurationSteps>0
    [state.mask,state.maskDiagnostics]=makeProtocolMasks(position,cfg); state.maskReady=true;
end
[obs,truth,ids]=sensorObservation(position,state,cfg,t);
[input,state,d]=controllerInputFromObservation(alpha,state,obs,cfg,t);
if ~cfg.collectDiagnostics, return; end
d.selected=obs.selected; d.fresh=obs.fresh; d.hidden=obs.selected & ~obs.fresh;
d.observedBearing=obs.bearing; d.observedWeight=obs.weight;
d.trueBearing=truth.trueBearing; d.trueDistance=truth.trueDistance;
d.predictionError=atan2(sin(d.predictionBearing-truth.trueBearing),cos(d.predictionBearing-truth.trueBearing));
d.priorError=atan2(sin(d.priorBearing-truth.trueBearing),cos(d.priorBearing-truth.trueBearing));
d.targetIds=ids; d.mask=state.mask; d.driverIndex=t;
end
