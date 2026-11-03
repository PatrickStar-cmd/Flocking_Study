"""Additional requested figures, paired history contrasts and Chinese report."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'v2'))
import analyze_v2 as v2

def main(run_id):
    run=ROOT/'results'/run_id; source=run/'analysis-v2'
    df=pd.read_csv(source/'by_seed.csv'); agg=pd.read_csv(source/'condition_means.csv')
    bins=pd.read_csv(source/'input_bins.csv')
    out=run/'delivery-v1';out.mkdir(exist_ok=False)
    cfg=v2.base.CFG; primary=cfg['primary']['profile'];z=df[df.profile==primary]
    # Retain prespecified empty bins explicitly, with missing estimates.
    complete=[]
    for profile,g in bins.groupby('profile'):
        t=g.set_index('input_bin').reindex(range(-1,19));t['profile']=profile
        t[['cases','seeds']]=t[['cases','seeds']].fillna(0).astype(int)
        complete.append(t.reset_index())
    pd.concat(complete).to_csv(out/'all_input_bins.csv',index=False)
    plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'pdf.fonttype':42,'font.size':11})
    def save(fig,name):
        fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=180);fig.savefig(out/(name+'.pdf'));plt.close(fig)
    b=bins[bins.profile==primary].sort_values('input_mean')
    fig,ax=plt.subplots(figsize=(10,6))
    ax.scatter(z.initialSocialInput,z.finalConsistency,s=5,alpha=.045)
    ax.plot(b.input_mean,b.consistency_mean,'o-',color='#0072b2');ax.fill_between(b.input_mean,b.lo,b.hi,alpha=.25)
    ax.set(xlabel='遮挡后初始三因素社会输入 S₀',ylabel='最终运动一致性 C_final',ylim=(0,1.02),title='主实验：500个种子；固定分箱与种子整体bootstrap 95%区间')
    save(fig,'主图_社会输入与最终一致性')
    # One set of cluster draws applies to every bin, metric and window.
    cols=['finalGO','finalLO','finalConsistency','C_last500','C_last2000']; seed_ids=sorted(z.seed.unique())
    zz=z.copy();x=zz.initialSocialInput.to_numpy();zz['bin']=np.where(x==0,-1,np.minimum(np.floor(x/.012).astype(int),18))
    bb=sorted(zz.bin.unique());gg=zz.groupby(['seed','bin'])
    count=gg.size().unstack(fill_value=0).reindex(index=seed_ids,columns=bb,fill_value=0).to_numpy(float)
    sums=np.stack([gg[c].sum().unstack(fill_value=0).reindex(index=seed_ids,columns=bb,fill_value=0).to_numpy(float) for c in cols],axis=-1)
    rng=np.random.default_rng(2026103999); boot=np.full((10000,len(bb),len(cols)),np.nan)
    for start in range(0,10000,100):
        ii=rng.integers(len(seed_ids),size=(100,len(seed_ids)));den=count[ii].sum(axis=1)
        num=sums[ii].sum(axis=1);boot[start:start+100]=np.divide(num,den[:,:,None],out=np.full_like(num,np.nan),where=den[:,:,None]>0)
    metric_rows=[]
    for j,bi in enumerate(bb):
        for k,c in enumerate(cols):
            lo,hi=np.nanquantile(boot[:,j,k],[.025,.975]);metric_rows.append(dict(bin=bi,metric=c,mean=sums[:,j,k].sum()/count[:,j].sum(),lo=lo,hi=hi))
    mm=pd.DataFrame(metric_rows);mm.to_csv(out/'binned_metrics_and_windows.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    for ax,metrics in zip(axes,[cols[:3],['C_last500','finalConsistency','C_last2000']]):
        for c in metrics:
            t=mm[mm.metric==c].set_index('bin').loc[bb];xx=b.set_index('input_bin').loc[bb].input_mean
            ax.plot(xx,t['mean'],'o-',label=c);ax.fill_between(xx,t.lo,t.hi,alpha=.12)
        ax.set(xlabel='初始三因素社会输入 S₀',ylim=(0,1.02));ax.legend(fontsize=8);ax.grid(alpha=.2)
    axes[0].set_title('GO、LO及逐时取最小值的C');axes[1].set_title('最后500、1000、2000步窗口敏感性')
    save(fig,'指标与时间窗口对照')
    # Across-history contrasts use the same first 100 seeds and resampling indices.
    idx=rng.integers(100,size=(10000,100)); history=[]
    sub=df[df.seed<2026103101]
    ref=sub[sub.profile==primary].pivot(index='seed',columns='remaining',values='finalConsistency').sort_index()
    for profile in ['formation','long_history']:
        t=sub[sub.profile==profile].pivot(index='seed',columns='remaining',values='finalConsistency').loc[ref.index,ref.columns]
        delta=t.to_numpy()-ref.to_numpy();samples=delta[idx].mean(axis=1)
        for j,r in enumerate(ref.columns):
            lo,hi=np.quantile(samples[:,j],[.025,.975]);history.append(dict(profile=profile,remaining=int(r),difference=delta[:,j].mean(),lo=lo,hi=hi,seeds=100))
    pd.DataFrame(history).to_csv(out/'paired_history_differences.csv',index=False)
    # Equal-complexity descriptive regressions, grouped by initial seed, without tuning.
    order=np.array(seed_ids);rng.shuffle(order);folds=np.array_split(order,5);cv={}
    for axis in ['initialSocialInput','initialTotalStrength']:
        errors=[]
        for fold in folds:
            test=z.seed.isin(fold);x=z[axis].to_numpy()/.228;X=np.column_stack([np.ones(len(x)),x,x*x]);y=z.finalConsistency.to_numpy()
            coef=np.linalg.lstsq(X[~test],y[~test],rcond=None)[0];errors.extend((X[test]@coef-y[test])**2)
        cv[axis]=float(np.sqrt(np.mean(errors)))
    (out/'grouped_descriptive_rmse.json').write_text(json.dumps(cv,indent=2),encoding='utf-8')
    # Independent direct evaluation of S0 using stored visible edges and event heading.
    seed=2026103001;profile=primary;r=8
    cp=loadmat(run/f'seed-{seed}'/f'{profile}-checkpoint.mat',simplify_cells=True)['checkpoint']
    raw=loadmat(run/f'seed-{seed}'/f'{profile}-r{r:02d}.mat',simplify_cells=True)
    ids=(np.arange(20)[:,None]+np.arange(1,20)[None,:])%20
    d=cp['position'][:,ids]-cp['position'][:,:,None];d-=1000*np.rint(d/1000)
    theta=np.arctan2(d[1],d[0]);S=float(np.mean(np.sum(.012*raw['meta']['visible']*np.exp(np.cos(theta-np.asarray(cp['heading'])[:,None])-1),axis=1)))
    expected=df[(df.seed==seed)&(df.profile==profile)&(df.remaining==r)].initialSocialInput.iloc[0]
    assert abs(S-expected)<1e-12
    (out/'independent_input_check.json').write_text(json.dumps(dict(status='PASS',S0=S,table_S0=expected,error=abs(S-expected))),encoding='utf-8')
    # Describe supported ranges; no forced monotonic curve or selected universal cutoff.
    supported=b[(b.seeds>=30)&(b.input_bin>=0)].sort_values('input_mean')
    gaps=np.diff(supported.consistency_mean);j=int(np.argmax(gaps)) if len(gaps) else 0
    sensitive=supported.iloc[j:j+2];high=supported.tail(3)
    text=['# 遮挡后初始社会输入与最终运动一致性：实验报告.v1','',
        '采用v2三因素输入与原始无补偿神经场；主实验500种子，历史对照各100种子，共14000条轨迹。每条遮挡后6000步，主终点为末1000步逐时min(GO,LO)的均值。',
        '区间采用10000次独立种子整体重采样；同一种子下20个遮挡条件一起抽取。精确零输入单列，正输入按0.012固定分箱。', '',
        '## 主实验与平台、敏感区','',
        '以下数值来自实际运行。正输入箱仅在至少30个独立种子覆盖时用于区间描述；全部箱（含空箱）的数量见all_input_bins.csv。','',
        '| S₀箱内均值 | C_final | 95%区间 | 轨迹数 | 独立种子数 |','|---|---|---|---|---|']
    for row in b.itertuples():text.append(f'| {row.input_mean:.6f} | {row.consistency_mean:.4f} | [{row.lo:.4f}, {row.hi:.4f}] | {row.cases} | {row.seeds} |')
    text+=['',f'高输入端三个有充分种子覆盖的箱，其C均值范围为{high.consistency_mean.min():.4f}至{high.consistency_mean.max():.4f}；对应S₀均值范围{high.input_mean.min():.6f}至{high.input_mean.max():.6f}。平台判断结合这些差值与逐箱区间，报告实际高端起伏。',
        f'支持箱中最大相邻均值上升发生于S₀={sensitive.input_mean.min():.6f}至{sensitive.input_mean.max():.6f}。这是分箱描述的敏感候选区间；两端95%区间与箱宽0.012共同给出当前分辨率。本轮保留完整非单调曲线。', '',
        '## 三因素输入S与总强度A','',f'相同轨迹、相同五折种子分组、同为二次回归的探索性留组预测RMSE：S={cv["initialSocialInput"]:.5f}，A={cv["initialTotalStrength"]:.5f}。更小的RMSE表示在该固定描述模型下更好的留种子预测。完整散点和条件均值对照见analysis-v2/social_input_and_strength_control.png。本批单边强度固定，A与邻居数成正比；S进一步记录遮挡时相对方位。','',
        '## 初始历史','']
    for profile in ['short_history','formation','long_history']:
        t=agg[agg.profile==profile];zero=t[t.remaining==0].iloc[0];full=t[t.remaining==19].iloc[0]
        text.append(f'{profile}：零保留边C={zero.finalConsistency:.4f}（95%区间{zero.finalConsistency_lo:.4f}至{zero.finalConsistency_hi:.4f}）；19条边C={full.finalConsistency:.4f}（{full.finalConsistency_lo:.4f}至{full.finalConsistency_hi:.4f}）。')
    text+=['', '历史间比较采用共同前100种子，逐保留数量的配对差值及95%区间见paired_history_differences.csv。每种历史单独解释其输入—响应关系。','',
        '## 窗口与复核','', '500/1000/2000步窗口结果及区间见指标与时间窗口对照.png和binned_metrics_and_windows.csv；主终点始终为末1000步。GO、LO、移动比例、平均速度及周期对距全部保存在condition_means.csv。',
        f'独立S₀复算误差{abs(S-expected):.3g}。逐步位置与速度独立重算的GO、LO、C检查见analysis-v2/audit.json；包含三个历史各r=0、8、19的原始轨迹。','',
        '## 后续KF对比覆盖区间','',f'优先覆盖上述敏感候选S₀区间{ sensitive.input_mean.min():.6f}至{sensitive.input_mean.max():.6f}及左右各一个0.012宽的分箱，同时保留零输入、全输入及高输入端对照。按共同检查点与原遮挡掩码进行配对比较，历史条件分层报告。','',
        '本轮为预先指定参数的模拟数据与探索性平台/敏感区描述。论文整合另行安排。原模型目录与历史结果保持原样，结果与软件修复保存在YCwork待审查；公共回退使用revert。']
    (out/'实验报告.v1.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    (out/'completion.json').write_text(json.dumps(dict(status='complete',trajectories=len(df),seed_count=int(df.seed.nunique()),bootstrap_repeats=10000,report='实验报告.v1.md'),indent=2),encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',required=True);main(p.parse_args().run)
