function checkPreparation()
% Static analyzer and constructed metric fixtures only: zero model steps.
root=fileparts(mfilename('fullpath')); addpath(root);
names={'runOcclusionConsistency.m','movementConsistency.m','summarizeConsistency.m','checkPreparation.m'};
diagnostics=struct(); parseErrors=0;
for i=1:numel(names)
    messages=checkcode(fullfile(root,names{i}),'-id');
    diagnostics.(erase(names{i},'.m'))=messages;
    for j=1:numel(messages)
        if strcmp(messages(j).id,'PARSE'), parseErrors=parseErrors+1; end
    end
end
assert(parseErrors==0);
L=1000; n=20; position=[10:10:200;zeros(1,n)]; v=[ones(1,n);zeros(1,n)];
a=movementConsistency(position,v,L); assert(abs(a.GO-1)<1e-14 && abs(a.LO-1)<1e-14 && abs(a.C-1)<1e-14);
stopped=movementConsistency(position,zeros(2,n),L); assert(stopped.C==0&&stopped.GO==0&&stopped.LO==0);
position=[10+(0:9) 400+(0:9);zeros(1,n)]; v=[ones(1,10) -ones(1,10);zeros(1,n)];
split=movementConsistency(position,v,L); assert(split.GO==0&&abs(split.LO-1)<1e-14&&split.C==0);
v=[ones(1,10) zeros(1,10);zeros(1,n)]; half=movementConsistency(position,v,L); assert(abs(half.GO-.5)<1e-14);
s=summarizeConsistency([1 0],[0 1],2); assert(s.finalConsistency==0 && s.finalGO==.5 && s.finalLO==.5);
% The default runner action only returns a plan and creates no trajectory.
plan=runOcclusionConsistency(); assert(plan.primaryTrajectories==10000 && plan.controlTrajectories==4000);
cfg=jsondecode(fileread(fullfile(root,'config.v1.json')));
assert(numel(cfg.fresh_counts)==20 && cfg.primary.seed_count==500);
assert(abs(plan.S0(end)-.228)<1e-14 && plan.S0(1)==0);
out=fullfile(root,'checks'); if ~exist(out,'dir'), mkdir(out); end
audit=struct('status','PASS','matlab',version,'model_steps_executed',0,'model_trajectories_executed',0, ...
 'static_parse_errors',parseErrors,'aligned',[a.GO a.LO a.C], ...
 'opposite_local_clusters',[split.GO split.LO split.C],'stationary',[stopped.GO stopped.LO stopped.C], ...
 'time_pairing_check',s.finalConsistency,'plan_only_check',true,'diagnostics',diagnostics);
fid=fopen(fullfile(out,'matlab-preparation.json'),'w'); fprintf(fid,'%s',jsonencode(audit,'PrettyPrint',true)); fclose(fid);
disp(audit);
end
