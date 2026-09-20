function plotCheckpointInterventionPilot()
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','checkpoint-intervention-pilot-v1');
t=readtable(fullfile(out,'results.csv'));
t=t(t.T>0,:);
durations=unique(t.T)';
modes={'none','hold'};
colors=[0.8500 0.3250 0.0980;0 0.4470 0.7410];

f=figure('Color','w','Position',[100 100 1400 560]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for mi=1:2
    means=zeros(size(durations));
    for di=1:numel(durations)
        v=t(strcmp(t.mode,modes{mi}) & t.T==durations(di),:);
        means(di)=mean(v.eventGOLoss);
        plot(repmat(durations(di),height(v),1),v.eventGOLoss,'o', ...
            'Color',[colors(mi,:) .35],'MarkerSize',5,'HandleVisibility','off');
    end
    semilogx(durations,means,'-o','Color',colors(mi,:),'LineWidth',2,'MarkerSize',7);
end
yline(0,':','LineWidth',1);
xlabel('Occlusion duration T'); ylabel('GO loss during occlusion');
title('Order response is not monotone');
legend({'Delete hidden inputs','Hold last bearing'},'Location','northwest'); grid on; box off;

nexttile; hold on;
for mi=1:2
    means=zeros(size(durations));
    for di=1:numel(durations)
        v=t(strcmp(t.mode,modes{mi}) & t.T==durations(di),:);
        means(di)=mean(v.postPairDistanceIncrease)*1000;
        plot(repmat(durations(di),height(v),1),1000*v.postPairDistanceIncrease,'o', ...
            'Color',[colors(mi,:) .35],'MarkerSize',5,'HandleVisibility','off');
    end
    semilogx(durations,means,'-o','Color',colors(mi,:),'LineWidth',2,'MarkerSize',7);
end
yline(0,':','LineWidth',1);
xlabel('Occlusion duration T'); ylabel('Pair-distance increase after occlusion');
title('Dispersion accumulates with blackout duration');
legend({'Delete hidden inputs','Hold last bearing'},'Location','northwest'); grid on; box off;

annotation(f,'textbox',[.32 .005 .36 .04], ...
    'String','Exploratory paired pilot; two previously used seeds', ...
    'EdgeColor','none','HorizontalAlignment','center','FontAngle','italic');
exportgraphics(f,fullfile(out,'checkpoint-pilot.png'),'Resolution',220);
close(f);
end
