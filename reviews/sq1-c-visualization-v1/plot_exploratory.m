% Visualize existing summaries only; no simulator is called.
root=fileparts(mfilename('fullpath'));
paired=readtable(fullfile(root,'..','sq1-c-paired-analysis-v1','paired_effects.csv'),'TextType','string');
deltas=readtable(fullfile(root,'..','sq1-c-paired-analysis-v1','paired_seed_differences.csv'),'TextType','string');
runs=readtable(fullfile(root,'..','sq1-c-condition-ledger-v1','B_E2_CANONICAL_RUNS.csv'),'TextType','string');
assert(height(paired)==33 && height(runs)==240 && height(deltas)==660);
set(groot,'defaultAxesFontName','Microsoft YaHei','defaultTextFontName','Microsoft YaHei');
set(groot,'defaultAxesFontSize',11,'defaultFigureColor','w');
metrics=["meanGO","meanLO","meanPairDistanceNorm"];
titles=["GO：全局有序度","LO：局部有序度","平均对距 / L"];
colors=[0.12 .39 .68; .12 .55 .47; .75 .36 .16];
ids=["strength_k2","strength_k4","strength_k8","strength_k12"];
labels=["k = 2","k = 4","k = 8","k = 12"];
forest(root,paired,deltas,ids,labels,metrics,titles,colors, ...
    '01_strength_effects','输入制度：fixed-total 相对 original 的配对差异');
ids=["policy_concentrated","policy_balanced","policy_nearest"];
labels=["角度集中","角度分散 / balanced","最近邻"];
forest(root,paired,deltas,ids,labels,metrics,titles,colors, ...
    '03_policy_effects','选邻策略：相对 random 的配对差异（k = 8，固定总幅值）');

f=figure('Visible','off','Position',[50 50 1500 740]);
layout=tiledlayout(1,3,'Padding','loose','TileSpacing','loose');
layout.OuterPosition=[0 .16 1 .84];
ks=[2 4 8 12 19]; cfgs=["C04","C05","C06","C07","C08"];
for m=1:3
    ax=nexttile;hold(ax,'on');vals=zeros(20,5);
    for j=1:5
        q=runs(runs.canonical_config_id==cfgs(j),:);q=sortrows(q,'seed');
        assert(height(q)==20 && isequal(q.seed,(1:20)'));
        vals(:,j)=q.(metrics(m));
    end
    light=.80+.20*colors(m,:);
    plot(ax,ks,vals','-','Color',light,'LineWidth',.65);
    plot(ax,ks,mean(vals,1),'-o','Color',colors(m,:),'LineWidth',2.6, ...
        'MarkerFaceColor',colors(m,:),'MarkerSize',7);
    xticks(ax,ks);xlim(ax,[1 20]);grid(ax,'on');box(ax,'off');
    title(ax,titles(m));xlabel(ax,'邻居数 k');
    if m<3,ylim(ax,[0 1]);else,ylim(ax,[0 max(vals(:))*1.13]);end
end
sgtitle('固定总幅值下的邻居数响应：20个种子与均值','FontSize',19,'FontWeight','bold');
annotation(f,'textbox',[.05 .005 .90 .09],'String', ...
    {'浅线：同一编号种子的结果；深线：20个种子的均值。连线仅帮助阅读离散条件，不表示中间k已被测试。', ...
    'k < 19 使用动态随机选邻，k = 19 为全邻居；每条件总幅值 A = 0.228。探索性数据，不代表最优k或普适阈值。'}, ...
    'EdgeColor','none','FontSize',11,'HorizontalAlignment','center');
savefigures(f,root,'02_neighbor_response');
disp('PASS: exported three PNG/SVG figures from fixed summaries; no simulation run.');

function forest(root,paired,deltas,ids,labels,metrics,titles,colors,name,heading)
f=figure('Visible','off','Position',[50 50 1600 740]);
layout=tiledlayout(1,3,'Padding','loose','TileSpacing','loose');
layout.OuterPosition=[0 .16 1 .84];
for m=1:3
    ax=nexttile;hold(ax,'on');xline(ax,0,'--','Color',[.4 .4 .4]);
    for j=1:numel(ids)
        p=paired(paired.contrast_id==ids(j) & paired.metric==metrics(m),:);
        d=deltas(deltas.contrast_id==ids(j) & deltas.metric==metrics(m),:);d=sortrows(d,'seed');
        assert(height(p)==1 && height(d)==20);
        assert(abs(mean(d.difference)-p.mean_difference)<1e-12);
        jitter=linspace(-.16,.16,20)';
        scatter(ax,d.difference,j+jitter,20,[.76 .80 .84],'filled');
        plot(ax,[p.bootstrap_ci95_low p.bootstrap_ci95_high],[j j], ...
            '-','Color',colors(m,:),'LineWidth',3.5);
        scatter(ax,p.mean_difference,j,65,colors(m,:),'filled','MarkerEdgeColor','w');
    end
    yticks(ax,1:numel(ids));yticklabels(ax,labels);ylim(ax,[.45 numel(ids)+.55]);
    set(ax,'YDir','reverse');grid(ax,'on');box(ax,'off');title(ax,titles(m));
    xlabel(ax,'处理条件 − 参照条件');
end
sgtitle(heading,'FontSize',19,'FontWeight','bold');
annotation(f,'textbox',[.03 .005 .94 .09],'String', ...
    {'浅点：20个种子的配对差；彩色点/线：平均差及逐项95%配对bootstrap区间（20,000次，未做多重比较校正）。', ...
    'GO/LO正值表示更有序；对距负值仅表示更紧凑。已用种子的探索性结果，不是独立验证或安全保证。'}, ...
    'EdgeColor','none','FontSize',11,'HorizontalAlignment','center');
savefigures(f,root,name);
end

function savefigures(f,root,name)
exportgraphics(f,fullfile(root,[name '.png']),'Resolution',160);
set(f,'PaperPositionMode','auto');print(f,fullfile(root,[name '.svg']),'-dsvg');
close(f);
end
