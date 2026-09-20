function plotPaperKFExtendedResults(profile)
%PLOTPAPERKFEXTENDEDRESULTS Figures and grouped summaries for KF sweep v3.

if nargin < 1, profile = "pilot"; end
root = fileparts(mfilename('fullpath'));
out = fullfile(root,'results','paper-kf-extended-v3',char(profile));
summaryPath = fullfile(out,'summary.csv');
assert(exist(summaryPath,'file')==2,'Missing summary.csv for profile %s.',profile);

t = readtable(summaryPath);
t.mode = string(t.mode);
modes = ["none","hold","mean","kalman","kalman-gm"];
labels = ["fresh-drop","hold-last","KF-mean","KF-unc","KF-GM-unc"];
hiddenCounts = unique(t.hiddenCount)';
Ts = unique(t.T)';

metrics = ["eventGODeviation","postPairDistanceDeviation"];
metricLabels = ["Mean absolute GO deviation during occlusion", ...
    "Mean absolute pair-distance deviation after occlusion"];

f = figure('Color','w','Position',[60 60 1550 820]);
tiledlayout(2,numel(hiddenCounts),'TileSpacing','compact','Padding','compact');
for hi = 1:numel(hiddenCounts)
    for pi = 1:numel(metrics)
        nexttile;
        hold on;
        for mi = 1:numel(modes)
            means = nan(size(Ts)); ci = nan(size(Ts));
            for ti = 1:numel(Ts)
                x = t(t.hiddenCount==hiddenCounts(hi) & t.T==Ts(ti) & t.mode==modes(mi),:);
                [means(ti),ci(ti)] = meanCI(x.(metrics(pi)));
            end
            errorbar(Ts,means,ci,'-o','LineWidth',1.5,'MarkerSize',5);
        end
        xlabel('Blackout duration T'); ylabel(metricLabels(pi)); grid on; box off;
        title(sprintf('%d/%d neighbor inputs hidden',hiddenCounts(hi),15));
        if pi==1 && hi==1, legend(labels,'Location','best'); end
    end
end
sgtitle('KF compensation under intermittent bearing observations');
exportgraphics(f,fullfile(out,'kf-performance-summary.png'),'Resolution',240);
close(f);

f = figure('Color','w','Position',[60 60 1500 430]);
tiledlayout(1,numel(hiddenCounts),'TileSpacing','compact','Padding','compact');
for hi = 1:numel(hiddenCounts)
    nexttile; hold on;
    for mode = ["hold","kalman","kalman-gm"]
        means = nan(size(Ts)); ci = nan(size(Ts));
        for ti = 1:numel(Ts)
            x = t(t.hiddenCount==hiddenCounts(hi) & t.T==Ts(ti) & t.mode==mode,:);
            [means(ti),ci(ti)] = meanCI(x.meanPredictionError);
        end
        errorbar(Ts,means,ci,'-o','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('Blackout duration T'); ylabel('Hidden-edge bearing error [rad]');
    grid on; box off; title(sprintf('%d/%d hidden',hiddenCounts(hi),15));
    if hi==1, legend(["hold-last","KF-unc","KF-GM-unc"],'Location','best'); end
end
sgtitle('Prediction error on hidden directed edges');
exportgraphics(f,fullfile(out,'kf-prediction-error.png'),'Resolution',240);
close(f);

f = figure('Color','w','Position',[60 60 1550 760]);
tiledlayout(2,numel(hiddenCounts),'TileSpacing','compact','Padding','compact');
compareModes = ["mean","kalman","kalman-gm"];
compareLabels = ["KF-mean","KF-unc","KF-GM-unc"];
for pi = 1:numel(metrics)
    for hi = 1:numel(hiddenCounts)
        nexttile; hold on;
        values = nan(numel(Ts),numel(compareModes));
        for ci = 1:numel(Ts)
            baseline = t(t.hiddenCount==hiddenCounts(hi) & t.T==Ts(ci) & t.mode=="none",:);
            for mi = 1:numel(compareModes)
                treated = t(t.hiddenCount==hiddenCounts(hi) & t.T==Ts(ci) & ...
                    t.mode==compareModes(mi),:);
                [tf,loc] = ismember(treated.seed,baseline.seed);
                baseValue = baseline.(metrics(pi))(loc(tf));
                value = treated.(metrics(pi))(tf);
                values(ci,mi) = (mean(baseValue)-mean(value))/max(eps,abs(mean(baseValue)));
            end
        end
        bar(Ts,100*values);
        xlabel('Blackout duration T'); ylabel('Reduction vs fresh-drop [%]');
        grid on; box off; title(sprintf('%d/%d hidden',hiddenCounts(hi),15));
        if pi==1 && hi==1
            legend(compareLabels,'Location','best');
        end
    end
end
sgtitle('Paired reduction in occlusion-induced deviation relative to fresh-drop');
exportgraphics(f,fullfile(out,'kf-relative-improvement.png'),'Resolution',240);
close(f);

f = figure('Color','w','Position',[60 60 1550 780]);
tiledlayout(2,numel(hiddenCounts),'TileSpacing','compact','Padding','compact');
compareModes = ["mean","kalman","kalman-gm"];
compareLabels = ["KF-mean - hold","KF-unc - hold","KF-GM-unc - hold"];
for pi = 1:numel(metrics)
    for hi = 1:numel(hiddenCounts)
        nexttile; hold on;
        for mi = 1:numel(compareModes)
            means = nan(size(Ts)); ci = nan(size(Ts));
            for ti = 1:numel(Ts)
                holdRows = t(t.hiddenCount==hiddenCounts(hi) & t.T==Ts(ti) & ...
                    t.mode=="hold",:);
                methodRows = t(t.hiddenCount==hiddenCounts(hi) & t.T==Ts(ti) & ...
                    t.mode==compareModes(mi),:);
                [tf,loc] = ismember(methodRows.seed,holdRows.seed);
                differences = methodRows.(metrics(pi))(tf) - holdRows.(metrics(pi))(loc(tf));
                [means(ti),ci(ti)] = meanCI(differences);
            end
            errorbar(Ts,means,ci,'-o','LineWidth',1.5,'MarkerSize',5);
        end
        yline(0,'k--','LineWidth',1);
        xlabel('Blackout duration T');
        if pi==1
            ylabel('GO deviation difference');
        else
            ylabel('Pair-distance deviation difference');
        end
        title(sprintf('%d/%d hidden (negative = better than hold)',hiddenCounts(hi),15));
        grid on; box off;
        if pi==1 && hi==1
            legend(compareLabels,'Location','best');
        end
    end
end
sgtitle('Incremental Kalman effect relative to hold-last (negative is better)');
exportgraphics(f,fullfile(out,'kf-versus-hold.png'),'Resolution',240);
close(f);

analyzeExtendedResults(t,out,hiddenCounts,Ts,modes);
end

function analyzeExtendedResults(t,out,hiddenCounts,Ts,modes)
metrics = ["eventGOLoss","eventGODeviation","postGOLoss","postGODeviation", ...
    "eventPairDistanceIncrease","eventPairDistanceDeviation", ...
    "postPairDistanceIncrease","postPairDistanceDeviation","meanPredictionError", ...
    "p90PredictionError","recoveryTime"];
rows = cell(0,7);
idx = 0;
for h = hiddenCounts
    for T = Ts
        for mode = modes
            x = t(t.hiddenCount==h & t.T==T & t.mode==mode,:);
            for metric = metrics
                v = x.(metric);
                v = v(isfinite(v));
                if isempty(v), continue; end
                idx = idx + 1;
                [mu,ci] = meanCI(v);
                rows(idx,:) = {h,T,mode,metric,numel(v),mu,ci};
            end
        end
    end
end
summary = cell2table(rows,'VariableNames', ...
    {'hiddenCount','T','mode','metric','n','mean','ci95'});
writetable(summary,fullfile(out,'group-summary.csv'));
end

function [mu,ci] = meanCI(x)
x = x(:);
n = numel(x);
mu = mean(x);
if n < 2
    ci = NaN;
else
    ci = 1.96*std(x)/sqrt(n);
end
end
