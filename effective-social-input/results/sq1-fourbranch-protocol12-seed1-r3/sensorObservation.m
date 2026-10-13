function [obs,truth,targetIds]=sensorObservation(position,state,cfg,t)
%SENSOROBSERVATION Only this sensor/evaluation layer can read geometry.
n=cfg.NAgents; targetIds=fixedTargetIds(cfg);
obs.selected=false(n,n); obs.fresh=false(n,n);
obs.bearing=nan(n,n); obs.weight=nan(n,n);
truth.trueBearing=nan(n,n); truth.trueDistance=nan(n,n);
event=t>=cfg.occlusionStartStep && t<cfg.occlusionStartStep+cfg.occlusionDurationSteps;
for i=1:n
    ids=targetIds(i,:); obs.selected(i,ids)=true;
    dxy=position(:,ids)-position(:,i); dxy=dxy-cfg.arenaSize*round(dxy/cfg.arenaSize);
    b=mod(atan2(dxy(2,:),dxy(1,:)),2*pi);
    truth.trueBearing(i,ids)=b; truth.trueDistance(i,ids)=hypot(dxy(1,:),dxy(2,:));
    hidden=false(size(ids));
    if event && ~strcmp(cfg.mode,'fresh') && state.maskReady
        hidden=ismember(ids,state.mask(i,:));
    end
    visible=~hidden; obs.fresh(i,ids(visible))=true;
    obs.bearing(i,ids(visible))=b(visible); obs.weight(i,ids(visible))=cfg.w0;
end
end
