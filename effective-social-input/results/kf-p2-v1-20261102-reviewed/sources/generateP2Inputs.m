function out=generateP2Inputs(suffix)
% Fixed physical records; estimator replay never changes the reference motion.
if nargin<1, suffix='reviewed'; end
assert(ischar(suffix) && ~isempty(regexp(suffix,'^[A-Za-z0-9_-]+$','once')));
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
src=fullfile(root,'effective-social-input'); p1=fullfile(root,'experiments','kf-diagnostic-20261101');
addpath(src); addpath(p1); addpath(here);
old=fullfile(src,'results','kf-diagnostic-v1-20261101-a01-reviewed');
out=fullfile(src,'results',['kf-p2-v1-20261102-' suffix]);
assert(~isfolder(out),'Refuse to overwrite any previous output.');
mkdir(out); mkdir(fullfile(out,'raw')); mkdir(fullfile(out,'raw','input'));
mkdir(fullfile(out,'logs')); mkdir(fullfile(out,'sources'));
diary(fullfile(out,'logs','input-generation.log')); clean=onCleanup(@()diary('off'));
copyfile(fullfile(here,'PROTOCOL.md'),fullfile(out,'PROTOCOL.md'));
[code,commit]=system('git rev-parse HEAD'); assert(code==0);
[code,branch]=system('git branch --show-current'); assert(code==0);
design=struct('dt',.3,'R',1e-4,'q_cv',1e-4,'q_gm',1e-5,'tau',10, ...
 'max_age',120,'max_variance',1,'durations',[20 120],'recovery_steps',400, ...
 'pre_steps',200,'bootstrap_seed',2026110200,'bootstrap_replicates',10000, ...
 'seeds',2026110101:2026110105,'evidence','exploratory-development', ...
 'source_commit',strtrim(commit),'branch',strtrim(branch), ...
 'case_selection','truth-only RMS wrapped first differences, candidates1201:20:1681,120steps; ties earliest', ...
 'turn_start',1401,'turn_duration',400,'turn_amplitude',.06,'turn_rate',.02, ...
 'turn_leaders',1:20,'N',20,'k',15,'hidden_count',7,'mask_rule','index-fixed', ...
 'matlab_version',version,'state_index','observation and positions at step t are before dynamics update t');
diagnosticJson(fullfile(out,'protocol.json'),design);
oldManifest=jsondecode(fileread(fullfile(old,'manifest.json')));
cases=struct('case_id',{},'seed',{},'stratum',{},'source',{},'t0',{},'selection_score',{});
inputs=struct('source',{},'sha256',{},'bytes',{},'upstream_mat',{},'upstream_sha256',{});
pairing=struct('seed',{},'exact_prefix_steps',{},'noise_equal',{},'turn_heading_change_rms',{},'turn_bearing_change_rms',{});
try
 for seed=design.seeds
  oldId=sprintf('seed%d_mode_fresh_D0',seed);
  upstream=fullfile(old,'raw',[oldId '.mat']);
  upstreamSha=diagnosticSha256(upstream);
  assert(strcmp(upstreamSha,oldManifest.runs.(oldId).sha256),'P1 source SHA mismatch.');
  a=load(upstream,'result','cfg'); base=a.result; baseCfg=a.cfg; clear a;
  [input,ids]=extract(base,seed,baseCfg.dt);
  baseRelative=sprintf('raw/input/seed%d_base.mat',seed);
  saveCompact(fullfile(out,baseRelative),input);
  inputs(end+1)=inputEntry(out,baseRelative,upstream,upstreamSha); %#ok<AGROW>
  rates=atan2(sin(diff(input.true_bearing)),cos(diff(input.true_bearing)))/design.dt;
  candidates=1201:20:1681; scores=zeros(size(candidates));
  for j=1:numel(candidates)
   t=candidates(j); block=rates(t-1:t+118,:); scores(j)=sqrt(mean(block(:).^2));
  end
  [lo,jlo]=min(scores); [hi,jhi]=max(scores);
  cases(end+1)=makeCase(seed,'S',baseRelative,candidates(jlo),lo); %#ok<AGROW>
  cases(end+1)=makeCase(seed,'V',baseRelative,candidates(jhi),hi); %#ok<AGROW>
  selection=table(candidates(:),scores(:),'VariableNames',{'t0','truth_bearing_rate_rms'});
  writetable(selection,fullfile(out,sprintf('selection_seed%d.csv',seed)));
  diagnosticJson(fullfile(out,sprintf('config_seed%d_base.json',seed)),baseCfg);
  turnCfg=baseCfg; turnCfg.turnExperiment=true; turnCfg.turnLeaderIds=1:20;
  turnCfg.turnCueStartStep=1401; turnCfg.turnCueDurationSteps=400;
  turnCfg.turnCueRate=.02; turnCfg.turnCueAmplitude=.06;
  turnCfg.turnCueBaseAngles=base.headingsRecord(:,1401);
  fprintf('Generating matched turn record seed %d\n',seed);
  turned=simulateEffectiveSocialInputInstrumented(turnCfg);
  assertDiagnosticPrefix(base,turned,1400);
  assert(isequal(base.sensorNoise,turned.sensorNoise));
  h=atan2(sin(turned.headingsRecord(:,1402:1801)-base.headingsRecord(:,1402:1801)), ...
      cos(turned.headingsRecord(:,1402:1801)-base.headingsRecord(:,1402:1801)));
  [input,turnIds]=extract(turned,seed,turnCfg.dt); assert(isequal(ids,turnIds));
  truthBase=baseTruth(base,ids); delta=atan2(sin(input.true_bearing(1402:1801,:)-truthBase(1402:1801,:)), ...
      cos(input.true_bearing(1402:1801,:)-truthBase(1402:1801,:)));
  pairing(end+1)=struct('seed',seed,'exact_prefix_steps',1400,'noise_equal',true, ...
      'turn_heading_change_rms',sqrt(mean(h(:).^2)),'turn_bearing_change_rms',sqrt(mean(delta(:).^2))); %#ok<AGROW>
  turnRelative=sprintf('raw/input/seed%d_turn.mat',seed);
  saveCompact(fullfile(out,turnRelative),input);
  inputs(end+1)=inputEntry(out,turnRelative,upstream,upstreamSha); %#ok<AGROW>
  r=atan2(sin(diff(input.true_bearing)),cos(diff(input.true_bearing)))/design.dt;
  b=r(1400:1519,:);
  cases(end+1)=makeCase(seed,'T',turnRelative,1401,sqrt(mean(b(:).^2))); %#ok<AGROW>
  diagnosticJson(fullfile(out,sprintf('config_seed%d_turn.json',seed)),turnCfg);
  fprintf('PASS paired source seed %d, S t0=%d V t0=%d, headRMS=%g\n',seed,candidates(jlo),candidates(jhi),pairing(end).turn_heading_change_rms);
  clear input base turned truthBase;
 end
 design.cases=cases;
 diagnosticJson(fullfile(out,'design.json'),design);
 diagnosticJson(fullfile(out,'input_manifest.json'),struct('status','COMPLETED','count',numel(inputs),'inputs',inputs,'pairing',pairing));
 writetable(struct2table(cases),fullfile(out,'cases.csv'));
 writetable(struct2table(pairing),fullfile(out,'turn_pairing.csv'));
 names={'bearingTrackStep.m','bearingTrackStepGaussMarkov.m','defaultEffectiveInputConfig.m'};
 for k=1:numel(names), copyfile(fullfile(src,names{k}),fullfile(out,'sources',names{k})); end
 names={'diagnosticConfig.m','simulateEffectiveSocialInputInstrumented.m','scheduledOcclusionInputInstrumented.m','occlusionInputStepInstrumented.m','assertDiagnosticPrefix.m','diagnosticSha256.m','diagnosticJson.m'};
 for k=1:numel(names), copyfile(fullfile(p1,names{k}),fullfile(out,'sources',names{k})); end
 copyfile(fullfile(here,'generateP2Inputs.m'),fullfile(out,'sources','generateP2Inputs.m'));
 fprintf('PASS %d sources, %d fixed cases: %s\n',numel(inputs),numel(cases),out);
catch ME
 diagnosticJson(fullfile(out,'FAILED.json'),struct('status','FAILED','error',getReport(ME,'extended','hyperlinks','off')));
 rethrow(ME);
end
end

function [data,ids]=extract(result,seed,dt)
N=20; ids=zeros(140,2); e=0;
for i=1:N
 for offset=1:7
  e=e+1; ids(e,:)=[i mod(i-1+offset,N)+1];
 end
end
idx=sub2ind([N N],ids(:,1),ids(:,2)); T=numel(result.observationRecords);
truth=zeros(T,140); observed=zeros(T,140);
for t=1:T
 rec=result.observationRecords{t}; assert(all(rec.fresh(idx)));
 truth(t,:)=rec.trueBearing(idx).'; observed(t,:)=rec.observedBearing(idx).';
end
data=struct('true_bearing',truth,'observed_bearing',observed,'observer_id',ids(:,1), ...
 'target_id',ids(:,2),'positions',result.positionsRecord,'headings',result.headingsRecord, ...
 'global_order',result.timeSeries.globalOrder,'mean_speed',result.timeSeries.meanSpeed, ...
 'seed',seed,'dt',dt,'arena_size',1000,'state_index','before update t');
end

function saveCompact(path,data)
assert(~isfile(path)); save(path,'-struct','data','-v7');
end

function entry=inputEntry(out,relative,upstream,sha)
info=dir(fullfile(out,relative));
entry=struct('source',relative,'sha256',diagnosticSha256(fullfile(out,relative)), ...
 'bytes',info.bytes,'upstream_mat',upstream,'upstream_sha256',sha);
end

function c=makeCase(seed,stratum,source,t0,score)
c=struct('case_id',sprintf('%s_seed%d',stratum,seed),'seed',seed,'stratum',stratum, ...
 'source',source,'t0',t0,'selection_score',score);
end

function truth=baseTruth(result,ids)
idx=sub2ind([20 20],ids(:,1),ids(:,2)); truth=zeros(numel(result.observationRecords),140);
for t=1:size(truth,1), rec=result.observationRecords{t}; truth(t,:)=rec.trueBearing(idx).'; end
end
