function plotPaperKFSimulation()
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results','paper-kf-simulation');
t=readtable(fullfile(out,'summary.csv')); modes={'none','hold','mean','kalman'};
labels={'fresh-drop','hold-last','KF-mean','KF-uncertainty'};
f=figure('Color','w','Position',[70 70 1500 800]); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
names={'eventGOLoss','postPairDistanceIncrease','eventInputMass','eventInputConcentration'};
ylabels={'GO loss during blackout','Post-blackout pair-distance increase','Mean injected input mass','Input concentration'};
for pi=1:4
 nexttile; hold on;
 for mi=1:numel(modes)
   means=[]; lo=[]; hi=[]; Ts=unique(t.T)';
   for T=Ts
    x=t(strcmp(string(t.mode),modes{mi}) & t.T==T,:); v=x.(names{pi});
    means(end+1)=mean(v); lo(end+1)=1.96*std(v)/sqrt(height(x)); hi(end+1)=lo(end);
   end
   errorbar(Ts,means,lo,hi,'-o','LineWidth',1.6,'MarkerSize',5);
 end
 xlabel('Blackout duration T'); ylabel(ylabels{pi}); grid on; box off;
 if pi==1, legend(labels,'Location','best'); end
end
sgtitle('PDF-aligned paired comparison under intermittent bearing observations');
annotation(f,'textbox',[.25 .005 .5 .03],'String','N=20, M=100, dt=0.3, k=15, 7 hidden bearings; seeds 101--105; mean +/- 95% normal CI','EdgeColor','none','HorizontalAlignment','center','FontAngle','italic');
exportgraphics(f,fullfile(out,'paper-kf-comparison.png'),'Resolution',240); close(f);
end
