function runSq1LocalE1()
%RUNSQ1LOCALE1 Generate the deterministic local mechanism replay (E1).
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','sq1-local-v1'); if ~exist(out,'dir'), mkdir(out); end
cfg=defaultLocalConfig();
cases=buildCases(cfg);
results=cell(size(cases)); rows=cell(0,12);
for i=1:numel(cases)
    results{i}=replayLocalSocialInput(cfg,cases{i}.schedule);
    r=results{i}; ts=r.timeSeries;
    rows(end+1,:)={cases{i}.name,mean(ts.inputMass),max(ts.fieldPeak), ...
        mean(ts.concentration),mean(ts.secondConcentration),mean(ts.maxAngularGap), ...
        mean(ts.angularEntropy),mean(ts.bumpConcentration), ...
        ts.decodedDirection(1),ts.decodedDirection(end), ...
        mean(ts.activeCount),mean(ts.totalWeight)}; %#ok<AGROW>
end
summary=cell2table(rows,'VariableNames',{'case','meanInputMass','maxFieldPeak', ...
    'meanZ1','meanZ2','meanMaxGap','meanAngularEntropy','meanBumpConcentration', ...
    'initialDecodedDirection','finalDecodedDirection','meanActiveCount','meanTotalWeight'});
writetable(summary,fullfile(out,'summary.csv'));
save(fullfile(out,'e1-results.mat'),'cfg','cases','results','summary','-v7');
writeReadme(out,cfg,cases);
disp(summary);
end

function cfg=defaultLocalConfig()
cfg.NNeurons=100; cfg.NSteps=120; cfg.dt=.3; cfg.beta=1000;
cfg.globalInhibition=0; cfg.receptiveFieldWidth=.4;
cfg.connectivityExponent=.5; cfg.initialNeuralState=zeros(cfg.NNeurons,1);
end

function cases=buildCases(cfg)
T=cfg.NSteps; t=(0:T-1).';
cases={};
theta=[0 pi/2 pi];
for scale=[.5 1 2]
    cases{end+1}=makeCase(sprintf('strength-scale-%.1f',scale), ...
        repmat(theta,T,1),repmat(scale*[.3 .4 .3],T,1)); %#ok<AGROW>
end
cases{end+1}=makeCase('structure-concentrated',repmat([0 0 0 0],T,1),repmat(.25,T,4));
cases{end+1}=makeCase('structure-bimodal',repmat([-pi/3 -pi/3 pi/3 pi/3],T,1),repmat(.25,T,4));
u=(0:7)*2*pi/8;
cases{end+1}=makeCase('structure-dispersed',repmat(u,T,1),repmat(1/8,T,8));
cases{end+1}=makeCase('label-split-equivalent',repmat([0 0 pi/2 pi/2],T,1),repmat([.15 .15 .35 .35],T,1));
fresh=mod(.15+.025*t,2*pi); hold=ones(T,1)*fresh(40);
cases{end+1}=makeCase('fresh-bearing',fresh,repmat(1,T,1));
cases{end+1}=makeCase('hold-last-bearing',hold,repmat(1,T,1));
end

function c=makeCase(name,theta,weights)
c.name=name; c.schedule.theta=theta; c.schedule.weights=weights;
end

function writeReadme(out,cfg,cases)
fid=fopen(fullfile(out,'README.md'),'w');
fprintf(fid,['# SQ1 local E1 replay v1\n\n' ...
    'Formal model definition: `../../SQ1_LOCAL_INPUT_MODEL.md`.\n\n' ...
    'Deterministic single-agent replay for input intensity and bearing structure. ' ...
    'No group positions, random seeds, KF, range, or occlusion are used.\n\n' ...
    '- Grid: NNeurons=%d, dt=%.3g, steps=%d, sigma=%.3g.\n' ...
    '- Input: `I(alpha)=sum_j w_j exp(-d(alpha,theta_j)^2/(2 sigma^2))`.\n' ...
    '- `inputMass` is the discrete angular integral; `z1/z2`, maximum gap and entropy describe bearing structure.\n' ...
    '- The label-split case tests identical field and neural response after splitting one weighted bearing into duplicate labels.\n' ...
    '- Fresh versus hold-last is a local time-pattern comparison, not a group performance or prediction claim.\n\n' ...
    'Cases: %d. Results are exploratory development evidence and are not a population threshold or stability guarantee.\n'], ...
    cfg.NNeurons,cfg.dt,cfg.NSteps,cfg.receptiveFieldWidth,numel(cases));
fclose(fid);
end
