function replay=buildReferenceReplay(fresh)
%BUILDREFERENCEREPLAY Four candidate fields on ONE fresh physical trajectory.
cfg=fresh.config; seg=fresh.segment; modes={'fresh','drop','drop-normalized','hold'};
replay.reference_id='fresh_seed1'; replay.step=seg.data.step;
replay.mask=fresh.checkpoint6200.mask; replay.alpha=seg.state.alpha;
replay.cacheAt6200=fresh.checkpoint6200.lastBearing;
for q=1:4
    c=cfg; c.mode=modes{q}; s=checkpointToStart(fresh.checkpoint6000); state=s.state;
    candidate.field=zeros(cfg.NNeurons,cfg.NAgents,numel(replay.step));
    candidate.relativeL2=nan(numel(replay.step),cfg.NAgents);
    names={'A','Mass','Peak','z1_real','z1_imag','z1_mod','z2_real','z2_imag','z2_mod','maxGap'};
    for name=names, candidate.(name{1})=nan(numel(replay.step),cfg.NAgents); end
    for j=1:numel(replay.step)
        [field,state,d]=fourBranchInputStep(seg.state.position(:,:,j),seg.state.alpha,state,c,replay.step(j));
        candidate.field(:,:,j)=field;
        den=sum(seg.data.inputField(:,:,j).^2,1);
        v=sqrt(sum((field-seg.data.inputField(:,:,j)).^2,1)./den); v(den==0)=nan;
        candidate.relativeL2(j,:)=v;
        for name=names, candidate.(name{1})(j,:)=d.(name{1}).'; end
    end
    replay.candidates.(matlab.lang.makeValidName(modes{q}))=candidate;
end
end
