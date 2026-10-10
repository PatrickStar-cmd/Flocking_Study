function [obs,truth,targetIds]=sensorObservation(position,state,cfg,t)
% 传感器和离线评价层读取几何；控制器只能收到obs，不能收到truth。
n=cfg.NAgents; targetIds=fixedTargetIds(cfg);
selected=false(n); ix=sub2ind([n n],repmat((1:n).',1,cfg.k),targetIds);
selected(ix)=true;
dx=position(1,:)-position(1,:).'; dy=position(2,:)-position(2,:).';
dx=dx-cfg.arenaSize*round(dx/cfg.arenaSize); dy=dy-cfg.arenaSize*round(dy/cfg.arenaSize);
truth.trueBearing=mod(atan2(dy,dx),2*pi); truth.trueDistance=hypot(dx,dy);
truth.trueBearing(~selected)=nan; truth.trueDistance(~selected)=nan;
fresh=selected;
event=t>=cfg.occlusionStartStep && t<cfg.occlusionStartStep+cfg.occlusionDurationSteps;
if event && ~strcmp(cfg.mode,'fresh') && state.maskReady
    if ~isempty(state.mask)
        hi=sub2ind([n n],repmat((1:n).',1,size(state.mask,2)),state.mask); fresh(hi)=false;
    end
end
obs.selected=selected; obs.fresh=fresh; obs.bearing=truth.trueBearing;
obs.bearing(~fresh)=nan; obs.weight=nan(n); obs.weight(fresh)=cfg.w0;
end
