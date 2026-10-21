"""Readable presentation of already-frozen, independently verified statistics."""
import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, r'C:\Users\Administrator\AppData\Local\Temp\noisy-bearing-kf-pydeps')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

font_manager.fontManager.addfont(r'C:\Windows\Fonts\msyh.ttc')
plt.rcParams.update({'font.family': 'Microsoft YaHei', 'axes.unicode_minus': False, 'font.size': 11})


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run',type=Path,required=True); args=ap.parse_args()
    root=args.run.resolve(); stats=root/'independent-analysis/statistics.json'
    assert json.loads((root/'independent-analysis/audit.json').read_text(encoding='utf-8'))['status']=='PASS'
    rows=json.loads(stats.read_text(encoding='utf-8'))['comparisons']
    out=root/'presentation'; assert not out.exists(); out.mkdir()
    shutil.copyfile(__file__,out/'plot-source.py')
    names={'oracle':'真实方位','raw':'带噪直接输入','selected-kf':'KF','imm':'IMM'}
    fig,ax=plt.subplots(figsize=(12,6.8))
    fig.subplots_adjust(left=.24,right=.97,bottom=.23,top=.85)
    ax.axvspan(-.05,.05,color='#d8e7dc',alpha=.7); ax.axvline(0,color='grey',lw=1)
    for i,r in enumerate(rows):
        ax.errorbar(r['difference'],i,xerr=[[r['difference']-r['family95Low']],[r['family95High']-r['difference']]],fmt='o',color='#245a81',capsize=5,lw=2)
        ax.plot([r['bootstrapFamily95Low'],r['bootstrapFamily95High']],[i+.12,i+.12],color='#b16d22',lw=2)
    ax.set_yticks(range(6),[names[r['a']]+' − '+names[r['b']] for r in rows]); ax.invert_yaxis()
    ax.set_xlim(-.14,.105); ax.set_xticks([-.1,-.05,0,.05,.1]); ax.grid(axis='x',alpha=.2)
    ax.set_xlabel('评价窗平均 GO 的配对差值；正值表示前者更高',labelpad=10)
    fig.suptitle('64个新种子：四组对照的六对差值',fontsize=16,y=.98)
    fig.text(.61,.91,'已考虑6次比较的同时95%区间；全部区间均包含0',ha='center',fontsize=11)
    fig.legend(handles=[Line2D([0],[0],marker='o',color='#245a81',label='预定的同时 t 区间'),
                        Line2D([0],[0],color='#b16d22',label='配对种子 bootstrap 敏感性区间'),
                        Patch(facecolor='#d8e7dc',label='本轮操作等效带 ±0.05')],
               loc='lower center',bbox_to_anchor=(.52,.035),ncol=3,frameon=False,fontsize=10)
    fig.savefig(out/'six-contrast-intervals-clear.png',dpi=160); plt.close(fig)
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    (out/'metadata.json').write_text(json.dumps({'statisticsSHA256':sha(stats),'plotSourceSHA256':sha(__file__),
        'change':'Legend moved outside axes; no statistics changed; original figures retained.',
        'sourceAudit':'independent-analysis/audit.json'},indent=2)+'\n',encoding='utf-8')
    print('Presentation figure saved; verified statistics and original figures unchanged.')


if __name__=='__main__': main()
