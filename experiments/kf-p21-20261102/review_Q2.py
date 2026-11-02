"""Correct only Q2 chronological plotting order, preserving its earlier version."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main(root):
    target=root/'review/figures/Q2_chronological.png';assert not target.exists()
    stage=root/'development-r01'
    coverage=pd.read_csv(stage/'analysis/coverage_by_age.csv.gz')
    coverage=coverage[(coverage.geometry=='clean-fast')&(coverage.D==120)]
    bins=pd.read_csv(stage/'analysis/p_bins.csv.gz');bins=bins[(bins.geometry=='clean-fast')&(bins.D==120)]
    order=['1-20','21-60','61-100','101-120'];fig,axes=plt.subplots(1,2,figsize=(11,4.7))
    for method in sorted(coverage.method.unique()):
        f=coverage[coverage.method==method].groupby(['seed','age_band']).coverage.mean().groupby('age_band').mean().reindex(order)
        axes[0].plot(range(4),f.values,'o-',label=method)
        b=bins[bins.method==method].groupby(['seed','P_lo'])[['mean_P','mean_error_squared']].mean().groupby('P_lo').mean()
        axes[1].loglog(b.mean_P,b.mean_error_squared,'o-',label=method)
    axes[0].set_xticks(range(4),order);axes[0].axhspan(.90,.98,color='gray',alpha=.12);axes[0].set_ylim(0,1.02)
    axes[0].set_xlabel('Occlusion age band (steps), chronological order');axes[0].set_ylabel('Approx. Gaussian coverage')
    axes[1].plot([1e-6,10],[1e-6,10],'k--');axes[1].set_xlabel('Mean P (rad^2)');axes[1].set_ylabel('Mean squared error (rad^2)')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Q2 FINAL: development clean-fast D120; 12 physical seeds\nDescriptive seed means; posterior geometry labels; shaded 90-98% is a reference',fontsize=10)
    fig.tight_layout(rect=(0,0,1,.88));fig.savefig(target,dpi=180,bbox_inches='tight');plt.close(fig)
    source=root/'P21_FINAL_REPORT.md';report=root/'P21_FINAL_REPORT_v2.md';assert not report.exists()
    report.write_text(source.read_text(encoding='utf-8')+'\n\n## 最终绘图复核补充\n\n正式Q2图为review/figures/Q2_chronological.png，按1–20、21–60、61–100、101–120步排序。初版analysis/figures/Q2_calibration.png按字符串排序，保留作QA历史，不能用于时间顺序解释。此次只修复图序并明确逐seed均值，不改变CSV、原始预测、冻结判据或q_star=null结论。正式入口见INDEX.md。\n',encoding='utf-8')
    print('PASS final chronological Q2 and additive report v2')


if __name__=='__main__':main(Path(sys.argv[1]))
