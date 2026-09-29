% Read-only comparison of every label against its frozen representative.
root=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(root));
out=fullfile(root,'alias_checks.json');
assert(~isfile(out),'Output exists; preserve it and use a new review directory.');
map=readtable(fullfile(root,'RAW_FILE_MAP.csv'),'TextType','string');
assert(height(map)==480);
ids=unique(map.canonical_run_id); assert(numel(ids)==240);
checks=cell(height(map),1); allPassed=true; nComparisons=0;
for gi=1:numel(ids)
    ix=find(map.canonical_run_id==ids(gi));
    repIx=ix(map.is_representative(ix)==1); assert(isscalar(repIx));
    d=load(fullfile(repo,map.raw_path(repIx)),'rawResult'); reference=d.rawResult;
    for j=ix'
        d=load(fullfile(repo,map.raw_path(j)),'rawResult'); r=d.rawResult;
        assert(r.config.seed==map.seed(j) && r.config.neighborCount==map.k(j));
        assert(string(r.config.strengthMode)==map.strengthMode(j));
        assert(string(r.config.selectionPolicy)==map.selectionPolicy(j));
        c=struct('canonical_run_id',char(ids(gi)),'raw_path',char(map.raw_path(j)), ...
            'representative_raw_path',char(map.raw_path(repIx)), ...
            'is_representative',j==repIx);
        % Config labels may differ only under the two declared model equivalences.
        [a,changedA]=normalizeConfig(r.config); [b,changedB]=normalizeConfig(reference.config);
        c.normalizedConfigExact=isequaln(a,b);
        c.configNormalization=unique([changedA changedB]);
        names=fieldnames(reference); names(strcmp(names,'config'))=[];
        c.payloadFieldNamesExact=isequal(sort(fieldnames(r)),sort(fieldnames(reference)));
        exact=true(size(names));
        for k=1:numel(names)
            exact(k)=isfield(r,names{k}) && isequaln(r.(names{k}),reference.(names{k}));
        end
        c.payloadFields=names; c.payloadFieldExact=exact;
        c.allPayloadExact=c.payloadFieldNamesExact && all(exact);
        c.n_agents=r.config.NAgents; c.n_neurons=r.config.NNeurons;
        c.n_steps=r.config.NSteps; c.burn_in=r.config.burnIn; c.dt=r.config.dt;
        c.time_series_length=numel(r.timeSeries.globalOrder);
        c.trajectory_frames=numel(r.trajectory.time); c.trajectory_stride_steps=r.config.recordStride;
        c.time_series_fields=fieldnames(r.timeSeries); c.trajectory_fields=fieldnames(r.trajectory);
        assert(c.n_steps==8000 && c.burn_in==6000 && c.time_series_length==8000 && c.trajectory_frames==81);
        for f=1:numel(c.time_series_fields)
            values=r.timeSeries.(c.time_series_fields{f});
            assert(numel(values)==8000 && all(isfinite(values)));
        end
        assert(max(abs(r.trajectory.time-(0:100:8000)*c.dt))<1e-10);
        checks{j}=c; allPassed=allPassed && c.allPayloadExact && c.normalizedConfigExact;
        nComparisons=nComparisons+(j~=repIx);
    end
end
report=struct('matlabVersion',version,'allPassed',allPassed, ...
    'canonicalRuns',numel(ids),'nonRepresentativeComparisons',nComparisons,'checks',{checks});
fid=fopen(out,'w','n','UTF-8'); assert(fid>0);
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true)); fclose(fid);
assert(allPassed,'Alias mismatch: inspect alias_checks.json before producing the canonical index.');
fprintf('PASS: 480 files, 240 canonical runs, %d alias comparisons; all saved payload fields exactly equal.\n',nComparisons);

function [c,changes]=normalizeConfig(c)
changes={};
if strcmp(c.selectionPolicy,'angle-dispersed')
    c.selectionPolicy='balanced'; changes{end+1}='angle-dispersed -> balanced';
end
if c.neighborCount==c.NAgents-1
    assert(strcmp(c.selectionPolicy,'random'));
    assert(ismember(c.strengthMode,{'original','fixed-total'}));
    c.strengthMode='full-equivalent'; changes{end+1}='k19 strengthMode -> full-equivalent';
end
end
