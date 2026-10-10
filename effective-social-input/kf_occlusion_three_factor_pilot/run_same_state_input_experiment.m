function [replay,rows]=run_same_state_input_experiment(fresh,cp,cfg,condition)
% 所有方法使用同一fresh物理轨迹，只比较输入；候选不会推进运动。
idx=find(fresh.data.step>6200); step=fresh.data.step(idx); rows=struct([]);
replay.reference_id=sprintf('fresh_seed%d_A%d',cfg.seed,round(cfg.A0*1000));
replay.step=step; replay.mask=cp.inputState.mask; replay.config=cfg;
replay.referencePosition=fresh.state.position(:,:,idx); replay.alpha=fresh.state.alpha;
names={'A','Mass','Peak','z1_real','z1_imag','z1_mod','z2_real','z2_imag','z2_mod','direction','directionValid','k_active'};
for q=1:5
    c=cfg; c.mode=c.modes{q}; state=cp.inputState; cand=struct();
    cand.field=zeros(cfg.NNeurons,cfg.NAgents,numel(idx)); cand.relativeL2=nan(numel(idx),cfg.NAgents);
    for f=names, cand.(f{1})=nan(numel(idx),cfg.NAgents); end
    for f={'predictionVariance','predictionError','predictionAge','predictionNIS'}
        cand.(f{1})=nan(cfg.NAgents,cfg.NAgents,numel(idx));
    end
    for jj=1:numel(idx)
        j=idx(jj); [field,state,d]=fourBranchInputStep(fresh.state.position(:,:,j),fresh.state.alpha,state,c,step(jj));
        cand.field(:,:,jj)=field; reference=fresh.data.inputField(:,:,j);
        den=sum(reference.^2,1); rel=sqrt(sum((field-reference).^2,1)./den); rel(den==0)=nan;
        cand.relativeL2(jj,:)=rel;
        for f=names, cand.(f{1})(jj,:)=d.(f{1}).'; end
        for f={'predictionVariance','predictionError','predictionAge','predictionNIS'}
            cand.(f{1})(:,:,jj)=d.(f{1});
        end
    end
    replay.candidates.(matlab.lang.makeValidName(c.mode))=cand;
    occ=step>=cfg.occlusionStartStep & step<cfg.recoveryStartStep;
    row.seed=c.seed; row.condition_id=condition.id; row.h=c.h; row.structure=condition.label;
    row.A0=c.A0; row.mode=c.mode; row.window='occlusion'; row.nSteps=sum(occ);
    row.relativeL2Mean=mean(cand.relativeL2(occ,:),'all','omitnan');
    row.relativeL2Max=max(cand.relativeL2(occ,:),[],'all');
    for f={'A','Mass','Peak','z1_mod','z2_mod','k_active'}
        row.(['mean_' f{1}])=mean(cand.(f{1})(occ,:),'all','omitnan');
    end
    if q==1, rows=row; else, rows(q)=row; end
end
end
