"""Frozen paired-seed analysis. No fitting/tuning and no exclusion on outcome."""
from pathlib import Path
import sys,json,csv
import numpy as np
import pandas as pd
import h5py
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
METRICS=['GO','GO_followers','follower_turn_error_deg','leader_turn_error_deg','fixed5_distance','pair_distance',
    'nearest_distance','minimum_distance','leader_follower_distance','S','S_abs_error','bearing_error_deg',
    'shadow_kf_rmse_rad','hidden_active_edges','gated_age_edges','gated_variance_edges']
MODES=['fresh','drop','drop-normalized','hold','linear-trend','kf-mean','kf-uncertainty']
COLORS={'fresh':'#202020','drop':'#d55e00','drop-normalized':'#cc79a7','hold':'#888888',
        'linear-trend':'#8c564b','kf-mean':'#0072b2','kf-uncertainty':'#009e73'}
CASES={'development':(20,18,200),'main':(20,18,200),'duration6':(20,6,200),'duration30':(20,30,200),
       'early':(20,18,140),'scale40':(40,18,200),'scale60':(60,18,200)}

def holm(p):
    p=np.asarray(p);order=np.argsort(p);ordered=np.maximum.accumulate((len(p)-np.arange(len(p)))*p[order]);out=np.empty_like(p)
    out[order]=np.minimum(1,ordered);return out

def analyze(out):
    out=Path(out);dest=out/'analysis';dest.mkdir(exist_ok=True);figdir=dest/'figures';figdir.mkdir(exist_ok=True)
    records=pd.read_csv(out/'RUN_MANIFEST.csv');assert len(records)==858 and (records.status=='COMPLETE').all()
    cache={};paths={};raw=[];curves=[];engineering=[]
    for r in records.itertuples():
        df=pd.read_csv(out/r.metrics);cache[(r.batch,r.seed,r.mode)]=df
        paths[(r.batch,r.seed,r.mode)]=out/r.data
    for case,(N,T,start) in CASES.items():
        seeds=[2026101401,2026101402,2026101403] if case=='development' else list(range(2026101501,2026101531))
        methods=[m for m in MODES if (case,seeds[0],m) in cache]
        if 'fresh' not in methods:methods=['fresh']+methods
        dur=round(T/.3);windows=dict(occlusion=(start,start+dur),recovery100=(start+dur,start+dur+100),terminal30=(770,800))
        for seed in seeds:
            refcase=case if case in ['development','main','scale40','scale60'] else 'main'
            fresh=cache[(refcase,seed,'fresh')].copy()
            if start!=200:
                # Reuse fresh physics while choosing the correct pre-occlusion neighbours.
                with h5py.File(paths[(refcase,seed,'fresh')]) as f:pos=f['states/xy'][:]
                delta=pos[start,:,None,:]-pos[start,None,:,:];delta-=1000*np.round(delta/1000)
                d=np.linalg.norm(delta,axis=2);np.fill_diagonal(d,np.inf);ids=np.argsort(d,axis=1)[:,:5]
                delta=pos[1:,:,None,:]-pos[1:,None,:,:];delta-=1000*np.round(delta/1000)
                distance=np.linalg.norm(delta,axis=3)
                fixed=distance[:,np.arange(N)[:,None],ids].mean(axis=(1,2));fixed[:start]=np.nan
                fresh['fixed5_distance']=fixed
            for mode in methods:
                df=fresh.copy() if mode=='fresh' else cache[(case,seed,mode)].copy()
                df['time']=(df.step+1-start)*.3
                for window,(a,b) in windows.items():
                    row=dict(case=case,N=N,T=T,start=start,seed=seed,mode=mode,window=window,nSteps=b-a,
                        referenceReused=mode=='fresh' and refcase!=case)
                    for metric in METRICS:
                        x=df[metric].to_numpy()[a:b];y=fresh[metric].to_numpy()[a:b]
                        row[metric]=float(np.nanmean(x)) if np.isfinite(x).any() else np.nan
                        row[metric+'_minus_fresh']=float(np.nanmean(x-y)) if np.isfinite(x-y).any() else np.nan
                        row[metric+'_abs_deviation_fresh']=float(np.nanmean(np.abs(x-y))) if np.isfinite(x-y).any() else np.nan
                    row['follower_angle_deg']=float(df.follower_turn_angle_deg.iloc[a:b].mean())
                    row['follower_resultant']=float(df.follower_resultant.iloc[a:b].mean())
                    row['minimum_distance_window_min']=float(df.minimum_distance.iloc[a:b].min())
                    raw.append(row)
                keep=['step','time','GO','GO_followers','LO','fixed5_distance','follower_turn_error_deg','S','S_abs_error',
                      'pair_distance','hidden_active_edges','gated_variance_edges']
                ts=df[keep].copy();ts['case']=case;ts['seed']=seed;ts['mode']=mode;curves.append(ts)
            if case in ['development','main']:
                drop=cache[(case,seed,'drop')];lf=fresh.leader_turn_error_deg.iloc[770:].mean();ff=fresh.follower_turn_error_deg.iloc[770:].mean()
                angle=fresh.follower_turn_angle_deg.iloc[770:].mean();resultant=fresh.follower_resultant.iloc[770:].min()
                event=(drop.follower_turn_error_deg-fresh.follower_turn_error_deg).iloc[start:start+dur].mean()
                stage1=lf<=10 and ff<=15 and angle>=30
                engineering.append(dict(case=case,seed=seed,leaderFinalError=lf,followerFinalError=ff,
                    followerFinalAngle=angle,followerMinResultant=resultant,taskGate=stage1,
                    dropEventErrorIncrease=event,dropGateEvaluated=bool(stage1),dropGate=bool(stage1 and event>=3)))
    raw=pd.DataFrame(raw);curve=pd.concat(curves,ignore_index=True)
    paired=[]
    for (case,seed,window),g in raw.groupby(['case','seed','window'],sort=False):
        s=g.set_index('mode')
        for method in ['kf-mean','kf-uncertainty','linear-trend']:
            if method not in s.index:continue
            baselines=['hold','fresh']
            if case in ['development','main','duration6','duration30']:baselines+=['drop','drop-normalized']
            if method.startswith('kf') and 'linear-trend' in s.index:baselines+=['linear-trend']
            if method=='kf-uncertainty':baselines+=['kf-mean']
            for baseline in baselines:
                if baseline==method or baseline not in s.index:continue
                for metric in METRICS:
                    paired.append(dict(case=case,N=int(g['N'].iloc[0]),T=int(g['T'].iloc[0]),seed=seed,window=window,
                        method=method,baseline=baseline,metric=metric,delta=s.loc[method,metric]-s.loc[baseline,metric],
                        methodValue=s.loc[method,metric],baselineValue=s.loc[baseline,metric]))
    paired=pd.DataFrame(paired);draws=np.random.default_rng(20261014).integers(0,30,(10000,30));aggregates=[]
    for (case,window,method,baseline,metric),g in paired.groupby(['case','window','method','baseline','metric'],sort=False):
        g=g.sort_values('seed');x=g.delta.to_numpy();good=np.isfinite(x);vals=x[good];n=len(vals)
        base=dict(case=case,window=window,method=method,baseline=baseline,metric=metric,nSeedsTotal=len(x),nFinite=n)
        if n==0:aggregates.append(base|dict(mean=np.nan,SD=np.nan,median=np.nan,minimum=np.nan,maximum=np.nan,tCI_low=np.nan,tCI_high=np.nan,bootstrap_low=np.nan,bootstrap_high=np.nan,pvalue=np.nan,nPositive=0,nNegative=0));continue
        mean=vals.mean();sd=vals.std(ddof=1) if n>1 else np.nan
        half=stats.t.ppf(.975,n-1)*sd/np.sqrt(n) if n>1 else np.nan
        p=2*stats.t.sf(abs(mean/(sd/np.sqrt(n))),n-1) if n>1 and sd>0 else (1. if mean==0 else np.nan)
        bl,bh=(np.quantile(vals[draws].mean(axis=1),[.025,.975]) if n==30 else [np.nan,np.nan])
        aggregates.append(base|dict(mean=mean,SD=sd,median=np.median(vals),minimum=vals.min(),maximum=vals.max(),
            tCI_low=mean-half if case!='development' else np.nan,tCI_high=mean+half if case!='development' else np.nan,
            bootstrap_low=bl,bootstrap_high=bh,pvalue=p if case!='development' else np.nan,
            nPositive=int((vals>1e-10).sum()),nNegative=int((vals< -1e-10).sum())))
    agg=pd.DataFrame(aggregates);agg['p_holm']=np.nan;agg['comparison_family']='descriptive'
    primary=(agg.case=='main')&(agg.window=='occlusion')&(agg.method=='kf-mean')&(agg.baseline=='hold')&(agg.metric=='GO')
    agg.loc[primary,'comparison_family']='primary'
    secondary=(agg.case=='main')&(agg.window=='occlusion')&(agg.metric=='GO')&(
        ((agg.method=='kf-mean')&agg.baseline.isin(['drop','drop-normalized','linear-trend']))|
        ((agg.method=='kf-uncertainty')&(agg.baseline=='kf-mean')))
    agg.loc[secondary,'p_holm']=holm(agg.loc[secondary,'pvalue']);agg.loc[secondary,'comparison_family']='main_secondary_4'
    expansion=agg.case.isin(['duration6','duration30','early','scale40','scale60'])&(agg.window=='occlusion')&(
        agg.method=='kf-mean')&(agg.baseline=='hold')&(agg.metric=='GO')
    agg.loc[expansion,'p_holm']=holm(agg.loc[expansion,'pvalue']);agg.loc[expansion,'comparison_family']='expansion_5'
    assert primary.sum()==1 and secondary.sum()==4 and expansion.sum()==5
    for name,df in [('window_metrics_by_seed',raw),('paired_differences_by_seed',paired),('paired_statistics',agg),
        ('time_series_all_conditions',curve),('engineering_diagnostics',pd.DataFrame(engineering))]:
        df.to_csv(dest/(name+'.csv'),index=False,encoding='utf-8-sig')
    np.save(dest/'bootstrap_seed_indices.npy',draws)
    plot_all(figdir,raw,paired,agg,curve,records,out)
    write_report(out,dest,raw,paired,agg,pd.DataFrame(engineering))
    print(agg[primary|expansion][['case','mean','tCI_low','tCI_high','pvalue','p_holm','nPositive','nNegative']].to_string(index=False),flush=True)
    return raw,paired,agg

def save(fig,path):
    fig.savefig(str(path)+'.png',dpi=170,bbox_inches='tight');fig.savefig(str(path)+'.pdf',bbox_inches='tight');plt.close(fig)

def plot_all(figdir,raw,paired,agg,curve,records,out,summary_only=False):
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    panels=[('GO','GO (all agents)'),('LO','Local direction order'),('fixed5_distance','Fixed-five distance (model length)'),
        ('follower_turn_error_deg','Follower target error (deg)'),('S','Social-input score S'),('S_abs_error','Own-state S absolute error')]
    fig,axes=plt.subplots(2,3,figsize=(14,7),constrained_layout=True)
    for ax,(metric,label) in zip(axes.flat,panels):
        for mode in MODES:
            g=curve[(curve.case=='main')&(curve['mode']==mode)&(curve.time>=-4)&(curve.time<=48)].groupby('time')[metric].agg(['mean','min','max'])
            ax.plot(g.index,g['mean'],color=COLORS[mode],label=mode,lw=1.3)
            ax.fill_between(g.index,g['min'],g['max'],color=COLORS[mode],alpha=.06)
        ax.axvspan(0,18,color='gray',alpha=.12);ax.axvline(18,color='gray',ls='--');ax.set_xlim(-4,48)
        ax.set_xlabel('Time relative to occlusion (step-end labels)');ax.set_ylabel(label)
    axes.flat[0].legend(fontsize=7)
    fig.suptitle('N=20; T=18; 30 evaluation seeds; 51 missing bearing edges\nMean lines, seed-range shading, no smoothing')
    save(fig,figdir/'main_six_panel')
    fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for ax,metric in zip(axes,['GO','GO_followers','follower_turn_error_deg']):
        g=paired[(paired.case=='main')&(paired.window=='occlusion')&(paired.method=='kf-mean')&(paired.baseline=='hold')&(paired.metric==metric)].sort_values('seed')
        ax.scatter(np.arange(1,31),g.delta,s=18,color=COLORS['kf-mean']);ax.axhline(0,color='gray',ls='--')
        ax.axhline(g.delta.mean(),color='black',lw=1);ax.set_xlabel('Evaluation seed index');ax.set_ylabel('KF mean minus Hold');ax.set_title(metric)
    save(fig,figdir/'main_paired_seed_effects')
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
    names=['main','duration6','duration30','early','scale40','scale60']
    for i,case in enumerate(names):
        g=paired[(paired.case==case)&(paired.window=='occlusion')&(paired.method=='kf-mean')&(paired.baseline=='hold')&(paired.metric=='GO')]
        s=agg[(agg.case==case)&(agg.window=='occlusion')&(agg.method=='kf-mean')&(agg.baseline=='hold')&(agg.metric=='GO')].iloc[0]
        ax.scatter(np.full(len(g),i)+np.linspace(-.12,.12,len(g)),g.delta,color='#0072b2',s=9,alpha=.35)
        ax.errorbar(i,s['mean'],yerr=[[s['mean']-s.tCI_low],[s.tCI_high-s['mean']]],color='black',fmt='D',capsize=4)
    ax.axhline(0,color='gray',ls='--');ax.set_xticks(range(6),['Main T18','T6','T30','Earlier onset','N40','N60'])
    ax.set_ylabel('Occlusion GO(KF mean) - GO(Hold)');save(fig,figdir/'robustness_effects')
    if summary_only:return
    # Complete tails for every seed and every case, not just a selected example.
    full=figdir/'all_seed_time_series';full.mkdir(exist_ok=True)
    for (case,seed),g in curve.groupby(['case','seed'],sort=False):
        fig,axes=plt.subplots(2,2,figsize=(10,6),constrained_layout=True)
        for ax,metric in zip(axes.flat,['GO','GO_followers','follower_turn_error_deg','fixed5_distance']):
            for mode in MODES:
                d=g[g['mode']==mode]
                if len(d):ax.plot(d.time,d[metric],label=mode,color=COLORS[mode],lw=1)
            T=CASES[case][1];ax.axvspan(0,T,color='gray',alpha=.1);ax.set_ylabel(metric);ax.set_xlabel('Relative time')
        axes.flat[0].legend(fontsize=6);fig.suptitle(f'{case}; seed {seed}; complete 800 steps')
        # Complete seed figures as PDF, lightweight preview for first predefined seed only.
        fig.savefig(full/f'{case}_seed{seed}.pdf',bbox_inches='tight')
        if seed in [2026101401,2026101501]:fig.savefig(full/f'{case}_seed{seed}.png',dpi=140,bbox_inches='tight')
        plt.close(fig)
    trajectories=figdir/'all_seed_trajectories';trajectories.mkdir(exist_ok=True)
    # All 858 trajectories; each case/seed page shares equal coordinate ranges.
    for (batch,seed),runs in records.groupby(['batch','seed'],sort=False):
        modes=list(runs['mode']);data=[];bounds=[];N=int(runs.N.iloc[0]);L=3 if N==20 else int(.15*N)
        for r in runs.itertuples():
            with h5py.File(out/r.data) as f:p=f['states/xy'][:]
            data.append((r,p));bounds.append(p[100:])
        low=np.concatenate(bounds).min(axis=(0,1));high=np.concatenate(bounds).max(axis=(0,1))
        fig,axes=plt.subplots(2,4,figsize=(15,8),constrained_layout=True)
        for ax,(r,p) in zip(axes.flat,data):
            for i in range(N):
                # Raw periodic coordinates retained; NaN breaks remove wrap-crossing strokes only.
                q=p[100:,i].copy();jumps=np.any(np.abs(np.diff(q,axis=0))>500,axis=1);q[np.r_[False,jumps]]=np.nan
                ax.plot(q[:,0],q[:,1],color='#d55e00' if i<L else '#0072b2',lw=.5,alpha=.6)
            for step,marker in [(100,'o'),(int(r.start),'s'),(int(r.start+r.durationSteps),'^')]:
                ax.scatter(p[step,:,0],p[step,:,1],s=8,marker=marker,color='black',alpha=.4)
            ax.set_xlim(low[0]-2,high[0]+2);ax.set_ylim(low[1]-2,high[1]+2);ax.set_aspect('equal');ax.set_title(r.mode,fontsize=8)
        for ax in axes.flat[len(data):]:ax.set_visible(False)
        fig.suptitle(f'{batch}; seed {seed}; red leaders / blue followers; circle task onset / square occlusion / triangle recovery')
        fig.savefig(trajectories/f'{batch}_seed{seed}.pdf',bbox_inches='tight');plt.close(fig)

def table(df):
    def fmt(x):return f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x)
    return '\n'.join(['|'+'|'.join(df.columns)+'|','|'+'|'.join(['---']*len(df.columns))+'|']+['|'+'|'.join(map(fmt,row))+'|' for row in df.itertuples(index=False,name=None)])

def write_report(out,dest,raw,paired,agg,engineering):
    sel=(agg.case=='main')&(agg.window=='occlusion')&(agg.method=='kf-mean')&(agg.baseline=='hold')
    core=agg[sel&(agg.metric=='GO')].iloc[0]
    verdict='支持该场景下KF-mean的平均GO改善' if core.tCI_low>0 else ('显示平均GO下降' if core.tCI_high<0 else '平均GO差的区间跨0，改善证据不足')
    metrics=agg[sel&agg.metric.isin(['GO','GO_followers','follower_turn_error_deg','fixed5_distance','pair_distance','S_abs_error'])]
    expansions=agg[(agg.window=='occlusion')&(agg.method=='kf-mean')&(agg.baseline=='hold')&(agg.metric=='GO')&(agg.case!='development')]
    main=raw[(raw.case=='main')&(raw.window=='occlusion')]
    means=main.groupby('mode')[['GO','GO_followers','follower_turn_error_deg','fixed5_distance','pair_distance','S_abs_error']].mean().reindex(MODES).reset_index()
    method_comparisons=agg[(agg.case=='main')&(agg.window=='occlusion')&(agg.metric=='GO')&(
        ((agg.method=='kf-mean')&agg.baseline.isin(['drop','drop-normalized','linear-trend']))|
        ((agg.method=='kf-uncertainty')&(agg.baseline=='kf-mean')))]
    text=['# 领航者方位缺测实验：实际结果','',
        '按v2方案完整执行858条科学运行：18条原开发复刻＋840条新评估。评估种子2026101501—2026101530，共30个独立种子；各方法/窗口/规模属于配对条件，不是858个独立随机样本。没有无转向机制对照。代码未据评估结果调参、换种子或改变窗口。', '',
        f'**主比较：{verdict}。遮挡60步内GO(KF-mean)−GO(Hold)均值为{core["mean"]:+.6f}，配对t的95%区间[{core.tCI_low:+.6f}, {core.tCI_high:+.6f}]；bootstrap区间[{core.bootstrap_low:+.6f}, {core.bootstrap_high:+.6f}]。30种子中{int(core.nPositive)}个正差、{int(core.nNegative)}个负差。**', '',
        'GO取值0—1，差+.01表示增加1个百分点，不代表相对提升1%。种子正差计数仅用1e-10数值容差，不是实际效益阈值。统计证据与实际幅度须分别判断。', '',
        '## 主场景：七方法的窗口均值','',table(means), '',
        '完整参考Fresh使用含噪当前观测；新鲜项所有方法都直接使用含噪观测。只有编号1—3接收余弦目标方向提示；跟随者没有目标指令。目标在遮挡初期达到45°后保持，各方法实际领航轨迹闭环演化，不是强制共用轨迹。', '',
        '## KF-mean相对Hold：效果量和区间','',table(metrics[['metric','mean','SD','tCI_low','tCI_high','bootstrap_low','bootstrap_high','nPositive','nNegative']]), '',
        'GO/GO_followers差越大越好；目标误差和输入误差差越小越好。距离差更小表示更紧凑，但不能据此保证安全或无分裂。不同指标不能互相替代。没有非劣效界限，不声称“其他性能完全无损”。', '',
        '## 较强预测基线与不确定性输入','',table(method_comparisons[['method','baseline','mean','tCI_low','tCI_high','p_holm']]),'',
        '简单Linear-trend已取得大部分相对Hold的GO收益；主场景KF-mean相对它的平均GO增益约0.000136，虽在本比较族下有统计证据，但实际幅度很小。KF-uncertainty相对KF-mean的平均GO差约0.000037，区间跨0，不能声称主场景下有额外GO优势。隐藏方位方差门控在主场景未触发，不能将收益归因于该门控。', '',
        '## 扩展条件：遮挡期KF-mean相对Hold','',table(expansions[['case','mean','tCI_low','tCI_high','pvalue','p_holm','nPositive','nNegative']]), '',
        '主比较未经多重校正；T6/T30/提前时机/N40/N60作为5项扩展比较族采用Holm，完整统计在paired_statistics.csv。时长复用Fresh，不增加独立参考样本；早期遮挡的冻结5邻居按第141步更新前重新选定；规模扩展同时改变布局与图分配，仅支持所定义的扩展场景。', '',
        '## 原协议工程诊断（不筛选种子）','',table(engineering), '',
        '工程阈值：Fresh终末30步领航目标误差≤10°、跟随目标误差≤15°、跟随者平均方向≥30°；通过后描述Drop遮挡期误差增量是否≥3°。所有种子均保留。小圆周合矢量或跨角边界须结合原始角度解释。', '',
        f'主场景30个评估种子中，Fresh任务工程检查通过{int(engineering[engineering.case=="main"].taskGate.sum())}个；Drop误差增量至少3°的检查通过{int(engineering[engineering.case=="main"].dropGate.sum())}个。绝大多数种子没有达到原工程规则所设的3°增量，说明不能将该场景一律描述为达到预设严重缺测影响；也没有据此删除样本。遮挡期跟随者对目标的绝对误差仍约47°，短期GO改善不能解读为已经准确完成目标方向跟随；终末窗口才用于评价最终响应。', '',
        '## 软件、复现与数据','',
        '- 原参考脚本哈希匹配；六个原方法与参考模型、T=0七方法、标量KF、隐藏真值隔离、趋势外推、独立积分验证均已执行。具体看checks/preflight。',
        '- independent检查逐运行从保存状态重算GO、跟随者GO、任务误差、对距、冻结邻居距离和社会输入S，并检查共享前缀、合法观测、权重与文件SHA256。最终状态以checks/independent/validation.json为准。',
        '- 原始运行runs/以HDF5无损保存完整状态、神经场、逐边估计和输入，每个运行另有800行CSV；inputs/按seed/N保存共享初态和完整噪声。',
        '- analysis/保存逐种子窗口统计、配对差、t区间、bootstrap及Holm结果。figures/保存六面板、种子效果图、扩展图；all_seed_time_series和all_seed_trajectories包含全部场景种子页面，不按结果挑示例。', '',
        '## 可用于论文的有边界表述','',
        f'在所定义的领航者目标方向引导与方位缺测场景中，KF-mean相对保持最后观测的遮挡期全群GO配对差为{core["mean"]:+.6f}（95%区间{core.tCI_low:+.6f}至{core.tCI_high:+.6f}）。该结果来自30组初态与观测噪声重复。是否同时改善跟随者任务响应，以及在不同遮挡长度、时机和群体规模下保持收益，应按上述实测表分别陈述。', '',
        '不研究“收益由转向机制造成”的因果结论。本轮是近似有序初态、固定图、人工方位缺测的抽象仿真；没有真实几何障碍、平台飞行、避碰或闭环稳定性证明。30种子仅覆盖规定的扰动分布。多实验数量不能替代论文的算法/任务贡献。']
    (ROOT/'RESULTS_说明.md').write_text('\n'.join(text),encoding='utf-8')
    (dest/'analysis_metadata.json').write_text(json.dumps(dict(primaryComparison='main/occlusion/GO/kf-mean-minus-hold',
        nEvaluationSeeds=30,bootstrapDraws=10000,bootstrapSeed=20261014,secondaryFamilySize=4,expansionFamilySize=5,
        referenceReuse='Fresh physics reused across T; early fixed-neighbour metric recomputed'),indent=2),encoding='utf-8')

if __name__=='__main__':analyze(ROOT/(sys.argv[1] if len(sys.argv)>1 else 'results/full_v2_run1'))
