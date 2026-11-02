"""Predeclared descriptive development/validation figures; no case selection."""
import json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot(parent):
    dev=parent/'development-r01'
    decision=json.loads((dev/'FROZEN_DECISION.json').read_text(encoding='utf-8'))
    dest=parent/'analysis/figures';assert not dest.exists();dest.mkdir(parents=True)
    seed=pd.read_csv(dev/'analysis/seed_metrics.csv');seed=seed[seed.scenario=='combined']
    sweep=pd.read_csv(dev/'analysis/q_sweep.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,group,D in zip(axes,['clean-fast','clean-slow'],[120,20]):
        f=seed[(seed.geometry==group)&(seed.D==D)]
        for q in sweep.q:
            y=f[f.method==f'CV_q{q:g}'].rmse
            ax.scatter(np.full(len(y),np.log10(q)),y,s=12,color='#0072b2',alpha=.5)
        means=[f[f.method==f'CV_q{q:g}'].rmse.mean() for q in sweep.q]
        ax.plot(np.log10(sweep.q),means,'o-',label='CV: seed means')
        ax.axhline(f[f.method=='Hold'].rmse.mean(),ls='--',color='black',label='Hold mean')
        ax.set_title(f'Development: {group}, D={D}; 12 seeds')
        ax.set_xlabel('log10 q');ax.set_ylabel('Mean-window circular RMSE (rad)');ax.legend(fontsize=8)
    fig.suptitle('Q1: same truth-labeled supports, all five Q candidates; points=physical seeds')
    fig.tight_layout();fig.savefig(dest/'Q1_q_grid.png',dpi=170);plt.close(fig)
    coverage=pd.read_csv(dev/'analysis/coverage_by_age.csv.gz')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    sub=coverage[(coverage.geometry=='clean-fast')&(coverage.D==120)]
    for method in sorted(sub.method.unique()):
        f=sub[sub.method==method].groupby(['seed','age_band']).coverage.mean().groupby('age_band').mean()
        axes[0].plot(f.index,f.values,'o-',label=method)
    axes[0].axhspan(.90,.98,color='gray',alpha=.12);axes[0].set_ylim(0,1.02)
    axes[0].set_ylabel('Approx. Gaussian coverage');axes[0].set_xlabel('Occlusion age band (steps)')
    bins=pd.read_csv(dev/'analysis/p_bins.csv.gz')
    sub=bins[(bins.geometry=='clean-fast')&(bins.D==120)]
    for method in sorted(sub.method.unique()):
        f=sub[sub.method==method].groupby('P_lo')[['mean_P','mean_error_squared']].mean()
        axes[1].loglog(f.mean_P,f.mean_error_squared,'o-',label=method)
    axes[1].plot([1e-6,10],[1e-6,10],'k--');axes[1].set_xlabel('Mean P (rad^2)');axes[1].set_ylabel('Mean squared error (rad^2)')
    axes[0].legend(fontsize=7);axes[1].legend(fontsize=7)
    fig.suptitle('Q2: development clean-fast; posterior geometry labels; 12 physical seeds\nDescriptive calibration, no independent edge/time inference')
    fig.tight_layout();fig.savefig(dest/'Q2_calibration.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    model='CV_q0.0001'
    for group in ['G-clean','G-switch','G-close','G-ambiguous']:
        f=seed[(seed.geometry==group)&(seed.D==120)&(seed.method==model)]
        x=['G-clean','G-switch','G-close','G-ambiguous'].index(group)
        axes[0].scatter(np.full(len(f),x),f.rmse,s=15,label=f'{group}: {f.seed.nunique()} seeds')
        axes[1].scatter(np.full(len(f),x),f.mean_window_p99,s=15)
    for ax in axes:
        ax.set_xticks(range(4),['clean','switch','close','ambiguous']);ax.set_xlabel('Posterior truth-based geometry class')
    axes[0].set_ylabel('Mean-window RMSE (rad)');axes[1].set_ylabel('Mean-window P99 absolute error (rad)');axes[0].legend(fontsize=7)
    fig.suptitle('G3: development default CV, D120; actual available seeds shown')
    fig.tight_layout();fig.savefig(dest/'G3_geometry_errors.png',dpi=170);plt.close(fig)
    nis=pd.read_csv(dev/'analysis/nis_tail_by_geometry.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    categories=['G-clean','switch_event','switch_near2','G-close','all']
    for ax,metric in zip(axes,['median','p99']):
        f=nis[nis.method=='CV_q0.0001'].groupby(['seed','geometry'])[metric].mean()
        for k,group in enumerate(categories):
            values=f.xs(group,level='geometry') if group in f.index.get_level_values('geometry') else []
            ax.scatter(np.full(len(values),k),values,s=15)
        ax.set_xticks(range(len(categories)),['clean','switch','switch +/-2','close','all'],rotation=15)
        ax.set_yscale('log');ax.set_ylabel(f'Fresh NIS {metric} (log)')
    fig.suptitle('G4: development default CV; full-fresh observation NIS, 12 seeds\nSame sources and noise; conditional descriptive tails, not Gaussian proof')
    fig.tight_layout();fig.savefig(dest/'G4_NIS_tails.png',dpi=170);plt.close(fig)
    rates=[]
    for p in sorted((dev/'labels').glob('*.npz')):
        with np.load(p) as x:rates.extend(x['rates'][x['labels']=='G-clean'].tolist())
    fig,ax=plt.subplots(figsize=(7,4));ax.hist(np.log10(np.maximum(rates,1e-12)),bins=70,color='#0072b2')
    for cut in decision['rate_cutpoints']:ax.axvline(np.log10(cut),ls='--',color='red')
    ax.set_xlabel('log10 clean continuous bearing-rate RMS (rad/model time)');ax.set_ylabel('Descriptive edge-window count')
    ax.set_title('G2: development truth-only terciles; 12 seeds, repeated windows')
    fig.tight_layout();fig.savefig(dest/'G2_clean_rate_thresholds.png',dpi=170);plt.close(fig)
    validation=parent/'validation-r01'
    fig,ax=plt.subplots(figsize=(8,4))
    if validation.exists() and (validation/'analysis/validation_contrasts.csv').exists():
        contrast=pd.read_csv(validation/'analysis/validation_contrasts.csv')
        f=contrast[(contrast.geometry=='clean-fast')&(contrast.D==120)&(contrast.scenario=='combined')&(contrast.method==f"CV_q{decision['q_star']:g}")]
        audit=json.loads((validation/'analysis/validation_audit.json').read_text(encoding='utf-8'))
        ax.scatter(f.seed.astype(str).str[-3:],f.delta_RMSE,s=18)
        mean=audit['primary']['mean'];lo=audit['primary']['ci_lower'];hi=audit['primary']['ci_upper']
        ax.axhline(mean,color='red',label=f'mean={mean:.4f}, seed-bootstrap95% [{lo:.4f},{hi:.4f}]')
        ax.axhline(0,color='black',ls='--');ax.tick_params(axis='x',rotation=90);ax.set_ylabel('Calibrated CV - Hold RMSE (rad)')
        ax.set_title(f"V1: independent validation, clean-fast D120, n={len(f)} seeds");ax.legend(fontsize=8)
    else:
        ax.axis('off');ax.text(.5,.5,'V1: NOT RUN\nDevelopment gate found no admissible global q.\nReserved validation seeds remain untouched.',ha='center',va='center',transform=ax.transAxes)
    fig.tight_layout();fig.savefig(dest/'V1_validation_status.png',dpi=170);plt.close(fig)
    print('PASS 6 overview figures; G0/G1 and raw age curves supplied separately')


if __name__=='__main__':plot(Path(sys.argv[1]))
