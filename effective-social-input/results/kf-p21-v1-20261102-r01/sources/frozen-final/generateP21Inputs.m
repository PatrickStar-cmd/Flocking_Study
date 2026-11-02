function out=generateP21Inputs(stage,seeds,suffix,runRoot)
% New Fresh references only, identical physical configuration to frozen P2.
assert(any(strcmp(stage,{'development','validation'})));
if nargin<3,suffix='r01';end
assert(ischar(suffix) && ~isempty(regexp(suffix,'^[A-Za-z0-9_-]+$','once')));
here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));
src=fullfile(root,'effective-social-input');p1=fullfile(root,'experiments','kf-diagnostic-20261101');
addpath(src);addpath(p1);
parent=fullfile(src,'results','kf-p21-v1-20261102-r01');
if nargin>=4,parent=char(runRoot);end
assert(isfile(fullfile(parent,'g0-r02','audit.json')),'G0 required.');
g=jsondecode(fileread(fullfile(parent,'g0-r02','audit.json')));assert(strcmp(g.status,'PASS'));
assert(isfile(fullfile(parent,'g0-r02','independent_labels.json')),'Independent G0 label audit required.');
assert(isfile(fullfile(parent,'g1','audit.json')),'G1 required.');
g=jsondecode(fileread(fullfile(parent,'g1','audit.json')));assert(strcmp(g.status,'PASS'));
validateattributes(seeds,{'numeric'},{'vector','integer','nonempty','>=',0,'<=',2^32-1});
assert(numel(unique(seeds))==numel(seeds));
if strcmp(stage,'development')
 assert(isequal(seeds(:).',2026111001:2026111012));
else
 frozen=fullfile(parent,'development-r01','FROZEN_DECISION.json');assert(isfile(frozen));
 expected=strtrim(fileread(fullfile(parent,'development-r01','FROZEN_DECISION.sha256')));
 assert(strcmp(diagnosticSha256(frozen),expected),'Frozen decision hash mismatch.');
 d=jsondecode(fileread(frozen));assert(~isempty(d.q_star),'No admissible q, validation prohibited.');
 assert(isequal(seeds(:),d.validation_seeds(:)),'Validation seeds must match frozen list.');
end
out=fullfile(parent,[stage '-' suffix]);assert(~isfolder(out),'Refuse overwrite.');
mkdir(out);mkdir(fullfile(out,'raw'));mkdir(fullfile(out,'raw','input'));mkdir(fullfile(out,'logs'));
diary(fullfile(out,'logs','generation.log'));cleanup=onCleanup(@()diary('off'));
inputs=struct('source',{},'sha256',{},'seed',{},'scenario',{},'bytes',{});
pairing=struct('seed',{},'prefix_equal',{},'sensor_noise_equal',{},'heading_change_rms',{});
manifest=struct('status','RUNNING','stage',stage,'seeds',seeds,'inputs',inputs);
diagnosticJson(fullfile(out,'input_manifest.json'),manifest);
try
for seed=seeds(:).'
 cfg=diagnosticConfig(seed,'p1');
 fprintf('START %s seed %d base\n',stage,seed);
 base=simulateEffectiveSocialInputInstrumented(cfg);
 relative=sprintf('raw/input/seed%d_base.mat',seed);
 compact=extract(base,seed);save(fullfile(out,relative),'-struct','compact','-v7');
 inputs(end+1)=entry(out,relative,seed,'base'); %#ok<AGROW>
 diagnosticJson(fullfile(out,sprintf('config_%d_base.json',seed)),cfg);
 cfg.turnExperiment=true;cfg.turnLeaderIds=1:20;cfg.turnCueStartStep=1401;
 cfg.turnCueDurationSteps=400;cfg.turnCueRate=.02;cfg.turnCueAmplitude=.06;
 cfg.turnCueBaseAngles=base.headingsRecord(:,1401);
 fprintf('START %s seed %d turn\n',stage,seed);
 turned=simulateEffectiveSocialInputInstrumented(cfg);
 assertDiagnosticPrefix(base,turned,1400);
 assert(isequal(base.sensorNoise,turned.sensorNoise));
 relative=sprintf('raw/input/seed%d_turn.mat',seed);
 compact=extract(turned,seed);save(fullfile(out,relative),'-struct','compact','-v7');
 inputs(end+1)=entry(out,relative,seed,'turn'); %#ok<AGROW>
 diagnosticJson(fullfile(out,sprintf('config_%d_turn.json',seed)),cfg);
 h=atan2(sin(base.headingsRecord(:,1402:1801)-turned.headingsRecord(:,1402:1801)),cos(base.headingsRecord(:,1402:1801)-turned.headingsRecord(:,1402:1801)));
 pairing(end+1)=struct('seed',seed,'prefix_equal',true,'sensor_noise_equal',true,'heading_change_rms',sqrt(mean(h(:).^2))); %#ok<AGROW>
 manifest.inputs=inputs;manifest.pairing=pairing;
 diagnosticJson(fullfile(out,'input_manifest.json'),manifest);
 fprintf('DONE %s seed %d (%d inputs)\n',stage,seed,numel(inputs));
 clear base turned compact;
end
manifest.status='COMPLETED';manifest.matlab_version=version;
[~,commit]=system('git rev-parse HEAD');manifest.source_commit=strtrim(commit);
manifest.seeds=seeds;diagnosticJson(fullfile(out,'input_manifest.json'),manifest);
writetable(struct2table(pairing),fullfile(out,'turn_pairing.csv'));
catch ME
 manifest.status='FAILED';manifest.error=getReport(ME,'extended','hyperlinks','off');
 diagnosticJson(fullfile(out,'input_manifest.json'),manifest);rethrow(ME);
end
end

function d=extract(r,seed)
i=repelem((1:20).',7);offset=repmat((1:7).',20,1);j=mod(i-1+offset,20)+1;
ix=sub2ind([20 20],i,j);truth=zeros(2200,140);observed=truth;
for t=1:2200
 rec=r.observationRecords{t};assert(all(rec.fresh(ix)));
 truth(t,:)=rec.trueBearing(ix).';observed(t,:)=rec.observedBearing(ix).';
end
d=struct('true_bearing',truth,'observed_bearing',observed,'observer_id',i,'target_id',j, ...
 'positions',r.positionsRecord,'headings',r.headingsRecord,'global_order',r.timeSeries.globalOrder, ...
 'mean_speed',r.timeSeries.meanSpeed,'seed',seed,'dt',.3,'arena_size',1000);
end

function e=entry(out,relative,seed,scenario)
p=fullfile(out,relative);info=dir(p);
e=struct('source',relative,'sha256',diagnosticSha256(p),'seed',seed,'scenario',scenario,'bytes',info.bytes);
end
