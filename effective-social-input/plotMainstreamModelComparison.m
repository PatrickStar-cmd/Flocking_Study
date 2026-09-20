function plotMainstreamModelComparison()
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results','mainstream-model-comparison');
t=readtable(fullfile(out,'summary.csv')); models=unique(string(t.model),'stable');
labels={'Proposed bearing-only KF-unc','Vicsek (heading)','Boids (position + heading)','Potential-consensus (position + heading)'};
metrics={'eventGO','eventGOLoss','postPairDistanceIncrease','minimumDistance','collisionFraction','controlEffort'};
ylabs={'GO during blackout','GO loss vs. no blackout','Pair-distance change vs. control','Minimum pair distance / L','Dangerous-encounter fraction (d<10)','Mean turning effort [rad/time]'};
f=figure('Color','w','Position',[60 60 1500 900]); tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for pi=1:numel(metrics)
 nexttile; means=zeros(numel(models),1); ci=means;
 for mi=1:numel(models)
  v=t.(metrics{pi})(string(t.model)==models(mi)); means(mi)=mean(v,'omitnan'); ci(mi)=1.96*std(v,'omitnan')/sqrt(sum(isfinite(v)));
 end
 bar(1:numel(models),means); hold on; errorbar(1:numel(models),means,ci,'.k','LineWidth',1.3);
 xticks(1:numel(models)); xticklabels({'Proposed','Vicsek','Boids','Potential-consensus'}); xtickangle(15);
 ylabel(ylabs{pi}); grid on; box off;
end
sgtitle('Intermittent-observation robustness versus representative mainstream flocking models');
annotation(f,'textbox',[.17 .003 .66 .035],'String','Same initial group state and T=30 blackout; mainstream baselines use richer observations and are not parameter-optimized','EdgeColor','none','HorizontalAlignment','center','FontAngle','italic');
exportgraphics(f,fullfile(out,'mainstream-comparison.png'),'Resolution',240); close(f);
end
