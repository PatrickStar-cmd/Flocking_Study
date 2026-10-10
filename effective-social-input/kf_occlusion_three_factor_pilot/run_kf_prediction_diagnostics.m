function [diagnostic,rows]=run_kf_prediction_diagnostics(fresh,cp,cfg,condition)
% fresh轨迹上进行纯离线预测：不推进群体，不从隐藏真值更新KF。
% 三个R均使用完整6200步相同合法历史。参数事前固定，不据结果挑R。
Rs=[1e-6 cfg.R 1e-2]; rows=struct([]); diagnostic=struct();
step=fresh.data.step; ix=find(step>=cfg.occlusionStartStep & step<cfg.recoveryStartStep);
ids=fixedTargetIds(cfg); selected=false(cfg.NAgents);
selected(sub2ind(size(selected),repmat((1:cfg.NAgents).',1,cfg.k),ids))=true;
hidden=false(size(selected)); hidden(sub2ind(size(hidden),repmat((1:cfg.NAgents).',1,cfg.h),cp.inputState.mask))=true;
for r=1:3
    c=cfg; c.R=Rs(r); st=cp.inputState;
    if Rs(r)~=cfg.R
        st.tracks=cp.inputState.sensitivityTracks.(sprintf('R_%d',round(-log10(c.R))));
    end
    error=nan(cfg.NAgents,cfg.NAgents,numel(ix)); variance=error; bearing=error; holdError=error;
    for jj=1:numel(ix)
        j=ix(jj); o.selected=selected; o.fresh=selected & ~hidden;
        o.bearing=fresh.data.trueBearing(:,:,j); truth=o.bearing; o.bearing(~o.fresh)=nan;
        o.weight=nan(size(selected)); o.weight(o.fresh)=c.w0;
        st.tracks=bearingTrackBankStep(st.tracks,o,c);
        e=atan2(sin(st.tracks.bearing-truth),cos(st.tracks.bearing-truth)); e(~hidden)=nan;
        error(:,:,jj)=e; v=st.tracks.P11; v(~hidden)=nan; variance(:,:,jj)=v;
        b=st.tracks.bearing; b(~hidden)=nan; bearing(:,:,jj)=b;
        he=atan2(sin(cp.inputState.lastBearing-truth),cos(cp.inputState.lastBearing-truth));
        he(~hidden)=nan; holdError(:,:,jj)=he;
    end
    % 恢复首步NIS使用预测创新，在更新前比较；不把无噪声NIS当传感器校准证明。
    j=find(step==cfg.recoveryStartStep); o.fresh=selected; o.bearing=fresh.data.trueBearing(:,:,j);
    o.weight=c.w0*double(selected); st.tracks=bearingTrackBankStep(st.tracks,o,c);
    nis=st.tracks.nis(hidden); e=error(isfinite(error)); v=variance(isfinite(variance)); he=holdError(isfinite(holdError));
    row.seed=c.seed; row.condition_id=condition.id; row.h=c.h; row.structure=condition.label;
    row.A0=c.A0; row.R=c.R; row.q=c.q; row.primaryR=c.R==cfg.R;
    row.warmup='full_6200_steps';
    row.nHiddenEdgeSteps=numel(e); row.rmseKF=sqrt(mean(e.^2)); row.rmseHold=sqrt(mean(he.^2));
    row.KFminusHoldRMSE=row.rmseKF-row.rmseHold;
    row.coverage95=mean(abs(e)<=1.96*sqrt(v)); row.meanVariance=mean(v);
    row.recoveryNISmean=mean(nis); row.recoveryNISmedian=median(nis);
    row.fullCircle95Fraction=mean(1.96*sqrt(v)>=pi);
    row.numericalFailures=sum(~isfinite(st.tracks.bearing(hidden)) | ~isfinite(st.tracks.P11(hidden)));
    if r==1, rows=row; else, rows(r)=row; end
    field=sprintf('R_%d',round(-log10(c.R)));
    diagnostic.(field)=struct('step',step(ix),'hidden',hidden,'error',error,'variance',variance,...
        'bearing',bearing,'holdError',holdError,'recoveryNIS',st.tracks.nis);
end
diagnostic.config=cfg; diagnostic.mask=cp.inputState.mask; diagnostic.rows=rows;
end
