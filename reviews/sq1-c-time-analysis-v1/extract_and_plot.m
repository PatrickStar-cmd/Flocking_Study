% Read-only exploratory time analysis; full-resolution curves, no simulation.
root=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(root));
assert(~isfile(fullfile(root,'run_windows.csv')),'Preserve existing outputs; use a new directory.');
idx=readtable(fullfile(repo,'reviews','sq1-c-raw-index-v1','CANONICAL_RAW_INDEX.csv'),'TextType','string');
old=readtable(fullfile(repo,'reviews','sq1-c-condition-ledger-v1','B_E2_CANONICAL_RUNS.csv'),'TextType','string');
effects=readtable(fullfile(repo,'reviews','sq1-c-paired-analysis-v1','paired_effects.csv'),'TextType','string');
cs=effects(effects.metric=="meanGO",:); assert(height(cs)==11);
fields={'globalOrder','localOrder','meanPairDistanceNorm'};
metrics={'meanGO','meanLO','meanPairDistanceNorm'};
labels={'全局有序性 GO','局部有序性 LO','归一化平均对距'};
windows=[6001 8000;6001 7000;7001 8000]; wn={'main','early','late'};
values=zeros(8000,20,12,3); seen=false(20,12); rows=cell(2160,9); nr=0; maxOldError=0;
for i=1:height(idx)
    c=str2double(extractAfter(idx.canonical_config_id(i),'C')); s=idx.seed(i);
    assert(~seen(s,c)); seen(s,c)=true;
    d=load(fullfile(repo,idx.raw_path(i)),'rawResult'); r=d.rawResult;
    assert(r.config.seed==s && r.config.NSteps==8000 && r.config.dt==0.3);
    oi=find(old.canonical_run_id==idx.canonical_run_id(i)); assert(isscalar(oi));
    for m=1:3
        v=r.timeSeries.(fields{m}); assert(numel(v)==8000 && all(isfinite(v)));
        values(:,s,c,m)=v;
        maxOldError=max(maxOldError,abs(mean(v(6001:8000))-old.(metrics{m})(oi)));
        for w=1:3
            seg=v(windows(w,1):windows(w,2)); nr=nr+1;
            rows(nr,:)={char(idx.canonical_config_id(i)),s,metrics{m},wn{w},windows(w,1),windows(w,2),mean(seg),std(seg),max(seg)-min(seg)};
        end
    end
end
assert(all(seen,'all') && nr==2160 && maxOldError<1e-12);
runTable=cell2table(rows,'VariableNames',{'config','seed','metric','window','start_step','end_step','mean','time_sd','time_range'});
writetable(runTable,fullfile(root,'run_windows.csv'));
set(groot,'defaultAxesFontName','Microsoft YaHei','defaultTextFontName','Microsoft YaHei');
steps=(1:8000)'; blue=[.08 .38 .7]; orange=[.85 .35 .08];
pairRows=cell(1980,8); np=0;
for j=1:height(cs)
    a=str2double(extractAfter(cs.treatment_config(j),'C')); b=str2double(extractAfter(cs.reference_config(j),'C'));
    f=figure('Visible','off','Color','w','Position',[30 30 1400 830]);
    layout=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
    layout.OuterPosition=[0 .11 1 .89];
    title(layout,sprintf('%02d  %s',j,cs.contrast(j)),'FontSize',17,'Interpreter','none');
    for m=1:3
        av=values(:,:,a,m); bv=values(:,:,b,m); dv=av-bv;
        ax=nexttile(layout,m); hold(ax,'on');
        h1=plot(ax,steps,mean(av,2),'Color',blue,'LineWidth',1.5);
        h2=plot(ax,steps,mean(bv,2),'Color',orange,'LineWidth',1.5);
        xline(ax,6000,'--','评价窗起点','Color',[.5 .2 .2],'LabelVerticalAlignment','bottom');
        xline(ax,7000,':','Color',[.5 .5 .5]); xlim(ax,[0 8000]);
        title(ax,labels{m}); ylabel(ax,'条件均值'); xlabel(ax,'模拟步（每步 0.3 时间单位）'); grid(ax,'on');
        if m<3, ylim(ax,[0 1]); end
        legend(ax,[h1 h2],{sprintf('处理 %s',cs.treatment_config(j)),sprintf('参照 %s',cs.reference_config(j))},'Location','best','FontSize',9);
        ax=nexttile(layout,m+3); hold(ax,'on');
        plot(ax,steps,dv,'Color',[.82 .82 .82],'LineWidth',.35);
        plot(ax,steps,mean(dv,2),'Color',blue,'LineWidth',1.7);
        yline(ax,0,'k--'); xline(ax,6000,'--','Color',[.5 .2 .2]); xline(ax,7000,':','Color',[.5 .5 .5]);
        xlim(ax,[0 8000]); ylabel(ax,'处理 − 参照'); xlabel(ax,'模拟步'); grid(ax,'on');
        if m<3, title(ax,'正值：处理更有序'); else, title(ax,'负值：处理更紧凑'); end
        for w=1:3
            ia=mean(av(windows(w,1):windows(w,2),:),1); ib=mean(bv(windows(w,1):windows(w,2),:),1);
            for s=1:20
                np=np+1; pairRows(np,:)={char(cs.contrast_id(j)),metrics{m},wn{w},s,ia(s),ib(s),ia(s)-ib(s),char(cs.contrast(j))};
            end
        end
    end
    annotation(f,'textbox',[.03 .015 .94 .065],'String', ...
        '上排：两条件各20个种子的均值。下排：灰线为各种子配对差，蓝线为平均差。所有曲线未平滑。\newline虚线6000后为原评价窗，点线7000分隔前后两半。曲线没有置信带；探索性结果，不代表持续成功或缺失恢复。', ...
        'EdgeColor','none','FontSize',11,'Interpreter','tex');
    stem=sprintf('%02d_%s',j,cs.contrast_id(j));
    exportgraphics(f,fullfile(root,[stem '.png']),'Resolution',140);
    print(f,fullfile(root,[stem '.svg']),'-dsvg'); close(f);
end
assert(np==1980); pairTable=cell2table(pairRows,'VariableNames',{'contrast_id','metric','window','seed','treatment','reference','difference','contrast'});
writetable(pairTable,fullfile(root,'paired_windows.csv'));
audit=struct('matlabVersion',version,'canonicalRuns',height(idx),'maxOldSummaryError',maxOldError,'runWindowRows',nr,'pairWindowRows',np,'figures',height(cs),'noSmoothing',true);
fid=fopen(fullfile(root,'extraction_audit.json'),'w'); fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); fclose(fid);
fprintf('PASS: 240 runs, %d run-window rows, %d paired rows, 11 figures; old-summary error %.3g.\n',nr,np,maxOldError);
