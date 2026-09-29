function results = run_B_1_experiments(mode)
%RUN_B_1_EXPERIMENTS Run B's E2 experiments 1-4.
% Use RUN_B_1_EXPERIMENTS('pilot') for a small exploratory run.
if nargin < 1, mode = 'formal'; end
thisDir = fileparts(mfilename('fullpath'));
sourceRoot = 'C:\Users\WHQ\Studying\大三上学期\无人机群\Flocking_Study-main\Flocking_Study-main\effective-social-input';
addpath(sourceRoot); addpath(thisDir);
if strcmp(mode,'pilot'), seeds=1:2; nSteps=1200; burnIn=600; recordStride=20; saveRaw=true;
elseif strcmp(mode,'formal'), seeds=1:20; nSteps=8000; burnIn=6000; recordStride=100; saveRaw=true;
else, error('mode must be pilot or formal'); end
if strcmp(mode,'formal'), out=fullfile(thisDir,'results_formal'); else, out=fullfile(thisDir,'results_pilot'); end; raw=fullfile(out,'raw');
if ~isfolder(out), mkdir(out); end; if saveRaw && ~isfolder(raw), mkdir(raw); end
base=defaultEffectiveInputConfig(); base.NAgents=20; base.NNeurons=100; base.NSteps=nSteps; base.burnIn=burnIn; base.recordStride=recordStride;
base.distanceWeightMode='equal'; base.visibilityMode='ideal'; base.occlusionEnabled=false; base.scheduledOcclusion=false; base.predictionMode='none'; base.allocentric=true;
rows=repmat(emptyRow(),0,1);
rows=runGrid(rows,'E1_strength_policy',[2 8 19],{'original','fixed-total'},{'random'},base,seeds,saveRaw,raw);
rows=runGrid(rows,'E2_neighbor_count',[2 4 8 12 19],{'original','fixed-total'},{'random'},base,seeds,saveRaw,raw);
rows=runGrid(rows,'E3_angular_structure',8,{'fixed-total'},{'angle-concentrated','random','angle-dispersed'},base,seeds,saveRaw,raw);
rows=runGrid(rows,'E4_selection_strategy',8,{'fixed-total'},{'random','nearest','balanced','angle-concentrated','angle-dispersed'},base,seeds,saveRaw,raw);
results=struct2table(rows); writetable(results,fullfile(out,'B_E2_runs.csv')); save(fullfile(out,'B_E2_runs.mat'),'results','mode','base','seeds','-v7'); writeReport(out,mode); writeAudit(out,mode,base,seeds,raw,thisDir);
fprintf('Finished. Results are in %s\n',out);
end
function rows=runGrid(rows,experiment,kValues,modes,policies,base,seeds,saveRaw,raw)
for im=1:numel(modes), for ip=1:numel(policies), for ik=1:numel(kValues), for is=1:numel(seeds)
cfg=base; cfg.seed=seeds(is); cfg.neighborCount=kValues(ik); cfg.strengthMode=modes{im}; cfg.selectionPolicy=policies{ip};
sim=simulateBGroupE2(cfg); s=sim.summary; r=emptyRow(); r.experiment=string(experiment); r.strengthMode=string(cfg.strengthMode); r.selectionPolicy=string(cfg.selectionPolicy); r.k=cfg.neighborCount; r.seed=cfg.seed;
r.meanGO=s.globalOrder; r.meanLO=s.localOrder; r.meanPairDistanceNorm=s.meanPairDistanceNorm; r.meanNearestDistanceNorm=s.meanNearestDistanceNorm; r.meanInputAmplitude=s.meanInputAmplitude; r.meanActiveNeighbors=s.meanActiveNeighbors; r.meanAngularCoverage=s.meanAngularCoverage; r.meanInputConcentration=s.meanInputConcentration; r.meanFieldPeak=s.meanFieldPeak; r.meanBumpConcentration=s.meanBumpConcentration; r.meanSpeed=s.meanSpeed; r.meanInputMass=s.meanInputMass; r.socialUnionLccFraction=s.socialUnionLccFraction; r.socialAlgebraicConnectivity=s.socialAlgebraicConnectivity; rows(end+1,1)=r;
if saveRaw, rawResult=sim; fn=sprintf('%s_%s_%s_k%d_seed%d.mat',experiment,cfg.strengthMode,cfg.selectionPolicy,cfg.neighborCount,cfg.seed); save(fullfile(raw,fn),'rawResult','-v7'); end
end,end,end,end
end
function r=emptyRow()
r=struct('experiment',"",'strengthMode',"",'selectionPolicy',"",'k',NaN,'seed',NaN,'meanGO',NaN,'meanLO',NaN,'meanPairDistanceNorm',NaN,'meanNearestDistanceNorm',NaN,'meanInputAmplitude',NaN,'meanActiveNeighbors',NaN,'meanAngularCoverage',NaN,'meanInputConcentration',NaN,'meanFieldPeak',NaN,'meanBumpConcentration',NaN,'meanSpeed',NaN,'meanInputMass',NaN,'socialUnionLccFraction',NaN,'socialAlgebraicConnectivity',NaN);
end
function writeReport(out,mode)
fid=fopen(fullfile(out,'B_E2_RESULT_EXPLANATION.md'),'w'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# B-E2 实验 1–4 结果说明\n\n运行模式：%s。\n\n',mode);
fprintf(fid,'这是无遮挡探索性实验，不包含遮挡、KF、预测补偿或观测缺失。\n\n');
fprintf(fid,'- 实验1：比较 original 与 fixed-total 的输入强度。\n- 实验2：扫描 k=2,4,8,12,19。\n- 实验3：固定 k=8 和总输入，比较集中、随机、分散方位结构。\n- 实验4：固定 k=8 和总输入，比较五种选邻策略。\n\n');
fprintf(fid,'主要指标：meanGO、meanLO、meanPairDistanceNorm、meanInputAmplitude、meanAngularCoverage、meanInputConcentration、socialUnionLccFraction 和 socialAlgebraicConnectivity。\n\n');
fprintf(fid,'先检查实验1的输入幅值是否符合设计，再解释实验2–4的群体指标。最终结论必须使用独立种子和冻结后的统计方案。\n');
end
function writeAudit(out,mode,base,seeds,rawDir,thisDir)
audit.mode = mode;
audit.generatedAt = datestr(now,31);
audit.runner = fullfile(thisDir,'run_B_1_experiments.m');
audit.simulator = fullfile(thisDir,'simulateBGroupE2.m');
audit.config = base;
audit.seeds = seeds;
audit.rawFileCount = numel(dir(fullfile(rawDir,'*.mat')));
audit.matlabVersion = version;
audit.outputDirectory = out;
audit.command = sprintf("run_B_1_experiments('%s')",mode);
fid = fopen(fullfile(out,'B_E2_formal_audit.json'),'w');
assert(fid>0);
fwrite(fid,jsonencode(audit,'PrettyPrint',true),'char');
fclose(fid);
end


