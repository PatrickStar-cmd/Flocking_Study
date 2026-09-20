function plotBearingDampedResults()
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','bearing-damped-model-v1');
s=readtable(fullfile(out,'summary.csv'));

f=figure('Color','w','Position',[100 100 1400 560]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
semilogx(s.T,s.holdRMSE,'-o','LineWidth',2,'MarkerSize',7); hold on;
semilogx(s.T,s.cvRMSE,'-s','LineWidth',2,'MarkerSize',7);
semilogx(s.T,s.dampedRMSE,'-^','LineWidth',2,'MarkerSize',7);
xlabel('Occlusion duration T');
ylabel('Circular bearing RMSE (rad)');
title('Development evaluation: mean forecast error');
legend({'Hold','Constant angular velocity','Damped angular velocity'},'Location','northwest');
grid on; box off;

nexttile;
semilogx(s.T,100*s.dampedCoverage95,'-o','LineWidth',2,'MarkerSize',7); hold on;
yline(95,'--','Nominal 95%','LineWidth',1.5);
xlabel('Occlusion duration T');
ylabel('Empirical interval coverage (%)');
title('Damped-model covariance consistency');
ylim([0 100]);
grid on; box off;

annotation(f,'textbox',[.31 .005 .38 .04], ...
    'String','Seeds 140-141; previously viewed development set', ...
    'EdgeColor','none','HorizontalAlignment','center','FontAngle','italic');
exportgraphics(f,fullfile(out,'prediction-comparison.png'),'Resolution',220);
close(f);
end
