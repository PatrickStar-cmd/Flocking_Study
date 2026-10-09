function finalizePackage(outDir,r2Dir)
%FINALIZEPACKAGE Post-run audit only; never changes trajectory/test data.
root=fileparts(mfilename('fullpath')); addpath(root); started=datetime('now');
rows=struct([]); modes={'fresh','drop','drop-normalized','hold'};
for q=1:4
    old=load(fullfile(r2Dir,[modes{q} '_seed1.mat']),'run');
    fresh=load(fullfile(outDir,[modes{q} '_seed1.mat']),'run');
    a=old.run.segment; b=fresh.run.segment;
    row.mode=modes{q}; row.r2_file_path=fullfile(r2Dir,[modes{q} '_seed1.mat']);
    row.r2_SHA256=fileSHA256(row.r2_file_path); row.r3_file_path=[modes{q} '_seed1.mat'];
    row.r3_SHA256=fileSHA256(fullfile(outDir,row.r3_file_path));
    row.cp6000_position_maxdiff=max(abs(old.run.checkpoint6000.position-fresh.run.checkpoint6000.position),[],'all');
    row.prefix_position_maxdiff=max(abs(a.state.position(:,:,1:201)-b.state.position(:,:,1:201)),[],'all');
    row.prefix_heading_maxdiff=max(abs(a.state.heading(:,1:201)-b.state.heading(:,1:201)),[],'all');
    row.prefix_input_maxdiff=max(abs(a.data.inputField(:,:,1:200)-b.data.inputField(:,:,1:200)),[],'all');
    row.full_position_maxdiff=max(abs(a.state.position-b.state.position),[],'all');
    displacement=a.state.position-b.state.position;
    displacement=displacement-b.config.arenaSize*round(displacement/b.config.arenaSize);
    row.full_position_periodic_maxdiff=max(abs(displacement),[],'all');
    row.full_input_maxdiff=max(abs(a.data.inputField-b.data.inputField),[],'all');
    row.final_u_maxdiff=max(abs(a.finalNeuralState-b.finalNeuralState),[],'all');
    row.mask_changed_observers=sum(any(sort(old.run.checkpoint6200(1).mask,2)~=sort(fresh.run.checkpoint6200.mask,2),2));
    row.note='r2 random candidate index order vs r3 sorted actual-ID tuples; r3 is regenerated from initial state, not reused r2 data';
    if q==1, rows=row; else, rows(q)=row; end
end
writetable(struct2table(rows),fullfile(outDir,'r2_r3_comparison.csv'));
x=load(fullfile(outDir,'reference_same_state_replay.mat'),'replay'); replay=x.replay;
ix=replay.step>=6201 & replay.step<=6234; rr=struct([]);
for q=1:4
    c=replay.candidates.(matlab.lang.makeValidName(modes{q}));
    row=struct(); row.mode=modes{q}; row.reference_id='fresh_seed1';
    row.window='occlusion6201_6234'; row.meanRelativeL2=mean(c.relativeL2(ix,:),'all','omitnan');
    row.maxRelativeL2=max(c.relativeL2(ix,:),[],'all'); row.meanA=mean(c.A(ix,:),'all');
    if q==1, rr=row; else, rr(q)=row; end
end
writetable(struct2table(rr),fullfile(outDir,'reference_replay_summary.csv'));
for name={'README.md','DATA_DICTIONARY.md','FIXES.md','load_pilot_example.m'}
    copyfile(fullfile(root,name{1}),fullfile(outDir,name{1}));
end
audit.command='finalizePackage(result_directory, original_r2_v3_directory)';
audit.scriptSHA256=fileSHA256(mfilename('fullpath')+".m"); audit.startTime=char(started);
audit.endTime=char(datetime('now')); audit.matlabVersion=version;
audit.scope='post-run r2 comparison, replay summary and document copies; trajectory data and test evidence not modified';
writeJson(fullfile(outDir,'postprocess_audit.json'),audit);
fprintf('Postprocess complete. r2/r3 fresh position diff: %.17g; changed mask observers: %d\n',rows(1).full_position_maxdiff,rows(1).mask_changed_observers);
end
