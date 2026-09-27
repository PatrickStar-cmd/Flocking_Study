function plotSq1LocalE1()
%PLOTSQ1LOCALE1 Plot compact input and local response diagnostics.
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','sq1-local-v1');
x=load(fullfile(out,'e1-results.mat'));
names=cellfun(@(c)c.name,x.cases,'UniformOutput',false);
pick=[1 4 5 6 8 9];
fig=figure('Visible','off','Color','w','Position',[100 100 1100 700]);
subplot(2,2,1); hold on;
for i=pick, plot(x.results{i}.timeSeries.inputMass,'LineWidth',1.2); end
xlabel('step'); ylabel('input mass'); grid on;
legend(names(pick),'Interpreter','none','Location','best');
subplot(2,2,2); hold on;
for i=pick, plot(x.results{i}.timeSeries.concentration,'LineWidth',1.2); end
xlabel('step'); ylabel('|z_1|'); ylim([0 1.05]); grid on;
subplot(2,2,3); hold on;
for i=pick, plot(x.results{i}.timeSeries.fieldPeak,'LineWidth',1.2); end
xlabel('step'); ylabel('field peak'); grid on;
subplot(2,2,4); hold on;
for i=[8 9], plot(x.results{i}.timeSeries.decodedDirection,'LineWidth',1.2); end
xlabel('step'); ylabel('decoded direction (rad)'); grid on;
legend(names([8 9]),'Interpreter','none','Location','best');
sgtitle('SQ1 local E1 deterministic replay');
saveas(fig,fullfile(out,'e1-mechanism.png')); close(fig);
end
