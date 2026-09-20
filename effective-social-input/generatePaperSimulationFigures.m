function generatePaperSimulationFigures()
%GENERATEPAPERSIMULATIONFIGURES Build reproducible figures for the paper draft.
% Results are separated into mechanism evidence, analytical KF properties,
% offline prediction diagnostics, and exploratory group intervention results.
root = fileparts(mfilename('fullpath'));
out = fullfile(root,'results','paper-simulation-figures');
if ~exist(out,'dir'), mkdir(out); end

%% Fig. 1: input mechanism sweep (100 seeds, N=20)
t = readtable(fullfile(root,'results','effective-input-n20-100.csv'));
[~,uniqueRows] = unique(t(:,{'seed','policy','strengthMode','k'}),'rows','stable');
t = t(uniqueRows,:);
ks = 0:19;
panels = {'meanInputAmplitude','meanAngularCoverage','meanInputConcentration', ...
          'meanGO','meanLO','meanPairDistanceNorm'};
ylabs = {'Mean injected input amplitude','Angular coverage [0,1]', ...
         'Input concentration [0,1]','Global order (GO)','Local order (LO)', ...
         'Mean pair distance / L'};
styles = {'-o','-s','--^'};
labels = {'random / original','random / fixed-total','balanced / fixed-total'};
queries = {strcmp(string(t.policy),'random') & strcmp(string(t.strengthMode),'original'), ...
           strcmp(string(t.policy),'random') & strcmp(string(t.strengthMode),'fixed-total'), ...
           strcmp(string(t.policy),'balanced') & strcmp(string(t.strengthMode),'fixed-total')};
f = figure('Color','w','Position',[70 70 1500 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for pi = 1:numel(panels)
    nexttile; hold on;
    for qi = 1:numel(queries)
        means = nan(size(ks)); cis = nan(size(ks));
        for ii = 1:numel(ks)
            v = t.(panels{pi})(queries{qi} & t.k==ks(ii));
            v = v(isfinite(v));
            means(ii) = mean(v); cis(ii) = 1.96*std(v)/sqrt(max(1,numel(v)));
        end
        errorbar(ks,means,cis,styles{qi},'LineWidth',1.5,'MarkerSize',4);
    end
    xlabel('Selected neighbors k'); ylabel(ylabs{pi}); grid on; box off;
    if pi==1, legend(labels,'Location','best'); end
    if pi==3, ylim([0 1]); end
end
sgtitle('Allocentric flocking mechanism sweep (N=20, 100 seeds)');
exportgraphics(f,fullfile(out,'fig1_input_mechanism.png'),'Resolution',240); close(f);

%% Fig. 2: analytical covariance growth and Fourier information retention
dt = 0.3; q = 1e-4; R = 1e-4; P = diag([R 1]);
F = [1 dt;0 1]; Q = q*[dt^3/3 dt^2/2;dt^2/2 dt]; H = [1 0];
for ii=1:500
  Pm=F*P*F'+Q; K=Pm*H'/(H*Pm*H'+R);
  P=(eye(2)-K*H)*Pm*(eye(2)-K*H)'+K*R*K';
end
p11=P(1,1); p12=P(1,2); p22=P(2,2);
T = linspace(0,40,1001);
V = p11 + 2*T*p12 + T.^2*p22 + q*T.^3/3;
gamma = exp(-V/2);
f = figure('Color','w','Position',[100 100 1200 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(T,V,'LineWidth',2); hold on;
for x = [3 10.2 30], xline(x,'--','LineWidth',1); end
xlabel('Blackout duration T'); ylabel('Bearing variance V(T) [rad^2]');
title('Covariance growth during missing observations'); grid on; box off;
nexttile; plot(T,gamma,'LineWidth',2); hold on;
yline(exp(-0.5),'--','gamma=0.607','LineWidth',1.2);
for x = [3 10.2 30], xline(x,'--','LineWidth',1); end
xlabel('Blackout duration T'); ylabel('First-harmonic retention gamma(T)');
title('Direction information retained by convolution'); ylim([0 1.02]); grid on; box off;
sgtitle('Analytical prediction-interface properties');
exportgraphics(f,fullfile(out,'fig2_covariance_information.png'),'Resolution',240); close(f);

%% Fig. 3: offline bearing prediction diagnostic
s = readtable(fullfile(root,'results','bearing-damped-model-v1','summary.csv'));
f = figure('Color','w','Position',[100 100 1200 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
semilogx(s.T,s.holdRMSE,'-o','LineWidth',1.8,'MarkerSize',6); hold on;
semilogx(s.T,s.cvRMSE,'-s','LineWidth',1.8,'MarkerSize',6);
semilogx(s.T,s.dampedRMSE,'-^','LineWidth',1.8,'MarkerSize',6);
xlabel('Blackout duration T'); ylabel('Circular bearing RMSE [rad]');
title('Offline prediction error'); legend({'hold-last','constant-omega KF','damped-omega KF'},'Location','northwest'); grid on; box off;
nexttile;
semilogx(s.T,100*s.dampedCoverage95,'-o','LineWidth',1.8,'MarkerSize',6); hold on;
yline(95,'--','Nominal 95%','LineWidth',1.2);
xlabel('Blackout duration T'); ylabel('Empirical 95% coverage [%]'); ylim([0 100]);
title('Covariance calibration remains incomplete'); grid on; box off;
sgtitle('Bearing-only prediction: medium-window mean improvement, no blanket superiority');
annotation(f,'textbox',[.29 .005 .42 .035],'String','Development evaluation seeds 140--141; exploratory evidence','EdgeColor','none','HorizontalAlignment','center','FontAngle','italic');
exportgraphics(f,fullfile(out,'fig3_prediction_offline.png'),'Resolution',240); close(f);

%% Fig. 3b: combined analytical and empirical uncertainty diagnostics
f = figure('Color','w','Position',[70 70 1450 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(T,V,'LineWidth',2); hold on;
for x = [3 10.2 30], xline(x,'--','LineWidth',1); end
xlabel('Blackout duration T'); ylabel('Bearing variance V(T) [rad^2]');
title('Analytical covariance growth'); grid on; box off;
nexttile; plot(T,gamma,'LineWidth',2); hold on;
yline(exp(-0.5),'--','gamma=0.607','LineWidth',1.2);
for x = [3 10.2 30], xline(x,'--','LineWidth',1); end
xlabel('Blackout duration T'); ylabel('First-harmonic retention');
title('Direction information retained'); ylim([0 1.02]); grid on; box off;
nexttile;
semilogx(s.T,s.holdRMSE,'-o','LineWidth',1.8,'MarkerSize',6); hold on;
semilogx(s.T,s.cvRMSE,'-s','LineWidth',1.8,'MarkerSize',6);
semilogx(s.T,s.dampedRMSE,'-^','LineWidth',1.8,'MarkerSize',6);
xlabel('Blackout duration T'); ylabel('Circular bearing RMSE [rad]');
title('Offline mean prediction error');
legend({'hold-last','constant-omega KF','damped-omega KF'},'Location','northwest'); grid on; box off;
nexttile;
semilogx(s.T,100*s.dampedCoverage95,'-o','LineWidth',1.8,'MarkerSize',6); hold on;
yline(95,'--','Nominal 95%','LineWidth',1.2);
xlabel('Blackout duration T'); ylabel('Empirical 95% coverage [%]'); ylim([0 100]);
title('Covariance calibration remains incomplete'); grid on; box off;
sgtitle('Bearing prediction diagnostics: mean accuracy and uncertainty calibration');
exportgraphics(f,fullfile(out,'kf-uncertainty-diagnostics.png'),'Resolution',240); close(f);

%% Fig. 4: paired group occlusion pilot
p = readtable(fullfile(root,'results','checkpoint-intervention-pilot-v1','results.csv'));
p = p(p.T>0,:); durations = unique(p.T)'; modes = {'none','hold'};
f = figure('Color','w','Position',[80 80 1450 760]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
metricNames = {'eventGOLoss','postPairDistanceIncrease','eventInputMassDifference','fullInputMass'};
yNames = {'GO loss during blackout','Post-blackout pair-distance increase','Input-mass change during blackout','Mean input mass'};
for mi = 1:4
    nexttile; hold on;
    for ci = 1:numel(modes)
        means = nan(size(durations));
        for di = 1:numel(durations)
            z = p(strcmp(string(p.mode),modes{ci}) & p.T==durations(di),:);
            means(di) = mean(z.(metricNames{mi}),'omitnan');
        end
        semilogx(durations,means,'-o','LineWidth',1.8,'MarkerSize',6);
    end
    xlabel('Blackout duration T'); ylabel(yNames{mi}); grid on; box off;
    if mi==1, legend({'delete hidden inputs','hold last bearing'},'Location','best'); end
end
sgtitle('Checkpoint intervention: stale input and input deletion have different costs');
annotation(f,'textbox',[.30 .005 .40 .03],'String','Exploratory paired pilot (N=20, seeds 136--137); not a confirmation of KF benefit','EdgeColor','none','HorizontalAlignment','center','FontAngle','italic');
exportgraphics(f,fullfile(out,'fig4_occlusion_pilot.png'),'Resolution',240); close(f);

%% Fig. 5: representative trajectories from stored time series
files = {fullfile(root,'results','time-series-validation-20260909','seed-101-condition-1.mat'), ...
         fullfile(root,'results','time-series-validation-20260909','seed-101-condition-2.mat')};
names = {'all neighbors, original input','zero fresh input'};
f = figure('Color','w','Position',[100 100 1200 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for fi=1:2
    d = load(files{fi},'result'); r=d.result; tr=r.trajectory;
    nexttile; hold on;
    for a=1:r.config.NAgents, plot(tr.x(a,:),tr.y(a,:),'LineWidth',0.7); end
    scatter(tr.x(:,1),tr.y(:,1),18,'k','filled'); scatter(tr.x(:,end),tr.y(:,end),24,'r','filled');
    xlim([0 r.config.arenaSize]); ylim([0 r.config.arenaSize]); axis square;
    xlabel('x'); ylabel('y'); title(names{fi}); grid on; box on;
end
sgtitle('Representative allocentric trajectories (periodic domain; red = final state)');
exportgraphics(f,fullfile(out,'fig5_trajectories.png'),'Resolution',240); close(f);

fid=fopen(fullfile(out,'README.md'),'w');
 fprintf(fid,['# Paper simulation figures\n\nGenerated by `generatePaperSimulationFigures.m`.\n\n' ...
 'Fig. 1 uses the N=20, 100-seed mechanism sweep. Fig. 2 is the analytical covariance formula. ' ...
 'Fig. 3 uses the offline bearing-prediction development summary. Fig. 3b combines analytical covariance ' ...
 'and empirical prediction calibration. Fig. 4 uses the paired checkpoint intervention pilot. ' ...
 'Fig. 5 uses stored time series for seed 101. Prediction and intervention figures are exploratory and are not universal stability or safety guarantees.\n']);
fclose(fid);
fprintf('Generated figures in %s\n',out);
end
