function plotTurnOcclusionPilot(profile)
%PLOTTURNOCCLUSIONPILOT Plot relative degradation for the controlled-turn test.
if nargin<1, profile="smoke"; end
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','turn-occlusion-pilot-v1',char(profile));
t=readtable(fullfile(out,'summary.csv')); t.mode=string(t.mode);
modes=["none","hold","mean","kalman","kalman-gm"];
labels=["fresh-drop","hold-last","KF-mean","KF-unc","KF-GM-unc"];
Ts=unique(t.T)'; h=unique(t.hiddenCount)';
metrics={"eventGORelativeDeviation","postPairRelativeDeviation","eventFollowerLeaderError"};
ylabels={"Relative GO deviation","Relative post-blackout pair-distance deviation", ...
 "Follower--leader angular error [rad]"};
f=figure('Color','w','Position',[60 60 1450 1050]); tiledlayout(3,numel(h),'TileSpacing','compact','Padding','compact');
for pi=1:3
 for hi=1:numel(h)
  nexttile; hold on;
  for mi=1:numel(modes)
   mu=nan(size(Ts)); ci=nan(size(Ts));
   for ti=1:numel(Ts)
    x=t(t.hiddenCount==h(hi)&t.T==Ts(ti)&t.mode==modes(mi),:);
    v=x.(metrics{pi}); mu(ti)=mean(v,'omitnan'); ci(ti)=1.96*std(v,'omitnan')/sqrt(max(1,sum(isfinite(v))));
   end
   errorbar(Ts,mu,ci,'-o','LineWidth',1.4,'MarkerSize',5);
  end
  xlabel('Blackout duration T'); ylabel(ylabels{pi}); grid on; box off;
  title(sprintf('%d/%d directed inputs hidden',h(hi),15));
  if pi==1&&hi==1, legend(labels,'Location','best'); end
 end
end
sgtitle('Controlled-turn occlusion: relative degradation from paired no-occlusion control');
annotation(f,'textbox',[.20 .005 .60 .03],'String','Turn begins 20 steps (6 model-time units) before blackout; leaders receive the same cue in control and all branches','EdgeColor','none','HorizontalAlignment','center','FontAngle','italic');
exportgraphics(f,fullfile(out,'relative-degradation.png'),'Resolution',240); close(f);

% Timeline diagram for explaining the causal intervention.
f=figure('Color','w','Position',[180 180 1100 300]); hold on;
pre=6; event=max(Ts); post=30;
patch([0 pre pre 0],[0 0 1 1],[.78 .88 .98],'EdgeColor','none');
patch([pre pre+event pre+event pre],[0 0 1 1],[.98 .78 .72],'EdgeColor','none');
patch([pre+event pre+event+post pre+event+post pre+event],[0 0 1 1],[.80 .93 .80],'EdgeColor','none');
xline(pre,'k--','LineWidth',1.3); xline(pre+event,'k--','LineWidth',1.3);
text(pre/2,.5,{'fresh observations +','leader turns'},'HorizontalAlignment','center');
text(pre+event/2,.5,{'blackout +','leader turns'},'HorizontalAlignment','center');
text(pre+event+post/2,.5,{'post-blackout','recovery'},'HorizontalAlignment','center');
xlim([0 pre+event+post]); ylim([0 1]); yticks([]); xlabel('Model time relative to turn start');
title('Causal timing of the controlled-turn intervention'); box on;
exportgraphics(f,fullfile(out,'intervention-timeline.png'),'Resolution',240); close(f);
end
