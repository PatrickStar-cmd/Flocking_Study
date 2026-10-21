import sys,json,csv,hashlib
from pathlib import Path
sys.path.insert(0,r'C:\Users\Administrator\AppData\Local\Temp\noisy-bearing-kf-pydeps')
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(r'D:\大三上专业课\卓越计划\effective-social-input\results\noisy-bearing-kf-rmse-diagnosis-20261014-v1')
source=Path(r'D:\大三上专业课\卓越计划\effective-social-input\results\noisy-bearing-kf-pilot-20261014-v2')
audit=json.loads((root/'replay-audit.json').read_text())
assert audit['status']=='PASS_EXACT_FILTER_REPLAY'
checks=[]; hashes={}
for row in audit['rows']:
    p=root/f"seed-{row['seed']}-replay-errors.mat"; hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    d=loadmat(p); e=d['errors'].ravel(); raw=d['rawErrors'].ravel(); jump=d['recentJump'].ravel().astype(bool)
    sq=e*e; top=np.partition(sq,-6000)[-6000:].sum()/sq.sum()
    values={'rawRMSE':np.sqrt(np.mean(raw**2)),'kfRMSE':np.sqrt(np.mean(sq)),
      'rawMedianAbs':np.median(np.abs(raw)),'kfMedianAbs':np.median(np.abs(e)),
      'kfBetterFraction':np.mean(np.abs(e)<np.abs(raw)),
      'topOnePercentSquaredErrorShare':top,'recentJumpSquaredErrorShare':sq[jump].sum()/sq.sum(),
      'noRecentJumpRMSE':np.sqrt(sq[~jump].mean())}
    err=max(abs(float(v)-row[k]) for k,v in values.items()); assert err<1e-12
    checks.append({'seed':row['seed'],'count':len(e),'maxDifference':err,'values':{k:float(v) for k,v in values.items()}})
assert sum(c['count'] for c in checks)==3600000
example=next(x for x in audit['worstExamples'] if x['seed']==5104)
pair=loadmat(source/'seed-5104-std-010-pair.mat',squeeze_me=True,struct_as_record=False)['pair']
r=pair.kf; d=loadmat(root/'seed-5104-replay-errors.mat'); center=example['step']
steps=np.arange(center-25,center+31); recv=example['receiver']-1; nei=example['neighbor']-1
truth=[]; obsErr=[]; kfErr=[]; priorErr=[]; angles=[]
for step in steps:
    x=r.trajectory.x[:,step-1];y=r.trajectory.y[:,step-1]; L=r.config.arenaSize
    dx=x[nei]-x[recv];dy=y[nei]-y[recv];dx-=L*np.round(dx/L);dy-=L*np.round(dy/L)
    angles.append(np.arctan2(dy,dx))
    ids=np.asarray(r.pilot.selectedHistory[:,:,step-1],dtype=int)-1
    sel=np.zeros((20,20),dtype=bool);sel[np.arange(20)[:,None],ids]=True
    indices=np.flatnonzero(sel.ravel(order='F')); target=recv+20*nei
    col=np.flatnonzero(indices==target)
    if len(col):
      ei=step-6001; obsErr.append(abs(d['rawErrors'][ei,col[0]]));kfErr.append(abs(d['errors'][ei,col[0]]));priorErr.append(abs(d['priorErrors'][ei,col[0]]))
    else:obsErr.append(np.nan);kfErr.append(np.nan);priorErr.append(np.nan)
change=np.abs(np.angle(np.exp(1j*np.diff(angles))))
fig,axes=plt.subplots(2,1,figsize=(9,6),layout='constrained',sharex=True)
axes[0].plot(steps[1:],change,'o-',markersize=3,color='#555555');axes[0].set_ylabel('True bearing step change (rad)')
axes[1].plot(steps,obsErr,'o-',markersize=3,label='Current noisy observation',color='#376c96')
axes[1].plot(steps,kfErr,'o-',markersize=3,label='KF posterior',color='#c2513c')
axes[1].plot(steps,priorErr,'--',label='KF prior',color='#b29860',alpha=.8)
for ax in axes:ax.axvline(center,color='gray',linestyle=':');ax.grid(alpha=.2)
axes[1].set_ylabel('Absolute angular error (rad)');axes[1].set_xlabel('Model step');axes[1].legend()
fig.suptitle('Archived example: seed 5104, receiver 13 -> neighbor 2\nPeriodic minimum-image switch at step 6768; gaps denote unselected edges',fontsize=11)
fig.savefig(root/'worst-error-example.png',dpi=180);plt.close(fig)
other=next(x for x in audit['worstExamples'] if x['seed']==5102)
r2=loadmat(source/'seed-5102-std-010-pair.mat',squeeze_me=True,struct_as_record=False)['pair'].kf
distances=[]
for step in (other['step']-1,other['step']):
    rec=other['receiver']-1;nei=other['neighbor']-1
    dx=r2.trajectory.x[nei,step-1]-r2.trajectory.x[rec,step-1]
    dy=r2.trajectory.y[nei,step-1]-r2.trajectory.y[rec,step-1]
    L=r2.config.arenaSize;dx-=L*np.round(dx/L);dy-=L*np.round(dy/L)
    distances.append(float(np.hypot(dx,dy)))
verify={'status':'PASS_SECOND_IMPLEMENTATION','points':3600000,'checks':checks,'fileHashes':hashes,
        'nearPassExample':{**other,'previousAndCurrentDistance':distances},
        'evidence':'post-hoc mechanistic diagnosis; points not independent samples'}
(root/'independent-recheck.json').write_text(json.dumps(verify,indent=2),encoding='utf-8')
rows=audit['rows']
text=['# 为什么本轮低噪声KF的RMSE更大？','',
'Material Passport：academic-research-suite / experiment-agent；2026-10-14。Verification Status：ANALYZED。事后机制诊断，不是新确认性实验。','',
'**主要原因是尾部大误差：KF在大部分观测点上去噪有效，但常角速度预测不能及时适应突变方位，留下少量很大的后验误差。RMSE对平方误差敏感，因此总体反而变大。**','',
'诊断复用原实验KF组全部8000步轨迹、完整选邻序列和完全相同的传感器随机流，重放原参数滤波器，不积分动力学、不改变q/R。各步RMSE与归档一致；评价窗口每种子600000个边×步误差，共3600000个点。独立统计单位仍是6个种子，不能把360万点当作独立重复。','',
f"滤波重放与归档RMSE最大差 {audit['maxArchivedRMSEError']:.3g}；Python第二实现重算RMSE、中位数、优于观测占比、前1%平方误差占比和近期周期跳变占比，最大差 {max(c['maxDifference'] for c in checks):.3g}。",'',
'| seed | 原始观测RMSE | KF RMSE | KF绝对误差中位数 | KF更准点占比 | 最大1%误差占总平方误差 | 近期周期跳变占总平方误差 |',
'|---|---:|---:|---:|---:|---:|---:|']
for r in rows:text.append(f"| {r['seed']} | {r['rawRMSE']:.6f} | {r['kfRMSE']:.6f} | {r['kfMedianAbs']:.6f} | {100*r['kfBetterFraction']:.1f}% | {100*r['topOnePercentSquaredErrorShare']:.1f}% | {100*r['recentJumpSquaredErrorShare']:.1f}% |")
text+=['','原始观测绝对误差中位数约.00674 rad；KF约.00432–.00465 rad。KF有65.8%–69.4%的点更准，说明不能把较大RMSE解释为处处更差。每个seed最大的1%误差贡献88.1%–97.1%的总平方误差；不删除这些误差，它们是方法真实代价。','',
'## 突变来自哪里','',
'周期空间使用最近图像的相对位移。当某个分量越过±L/2，代表邻居的最近图像可能切换，感知方位可突然改变接近pi；这是几何代表的改变，不是个体撞墙或物理折返。代码虽对±pi/2pi的角度环绕作处理，但真实代表方位变化接近pi并非普通数值包裹问题，wrapped innovation不会消除这个变化。','',
'近期跳变定义为当前或此前10步（共11时刻）最短图像相对位移有分量改变超过L/2。在4个seed中，这些点贡献88.4%–93.5%总平方误差；另外两个seed为15.4%和49.6%，不能把所有误差全归因于周期边界。该时间窗是事后诊断选择，不作因果份额或总体阈值。','',
f"非周期跳变例：seed5102、第6958步，个体14对2的真实方位一步变化 {other['currentTrueStepChange']:.6f} rad；两时刻最近距离为 {distances[0]:.6g}/{distances[1]:.6g} 个模型距离单位。这是近距离运动时的快速几何变化，观测误差仅 {other['observationError']:.6f} rad，KF后验误差却为 {other['posteriorError']:.6f} rad。",'',
'## 为什么观测很准，KF仍保留大误差','',
'低噪声评价窗平均当前观测增益约.487，即后验保留约一半预测方向、引入约一半新观测。在局部小误差条件下，e_post≈(1−K)e_prior+K*measurement_noise。模型预测正确时能去噪；方向突变后，错误预测仍有权重，当前正确观测不能在一步内完全纠正。固定Q/R不根据残差幅度自动改变增益。','',
f"具体归档例：seed5104，第6768步，真实方位变化 {example['currentTrueStepChange']:.6f} rad，当前观测误差 {example['observationError']:.6f} rad，KF后验误差 {example['posteriorError']:.6f} rad。",'',
'![归档突变与误差](worst-error-example.png)','',
'## 解释边界','',
'这里诊断的是噪声std=.01、q=1e-4、随机k15、周期空间、二状态CV-KF。真实邻居方位受双方速度和距离共同影响，不等于邻居航向，不保证匀角速度。较高噪声.10的总体RMSE均值原本较小，本轮没有重放其全部轨迹。需要先在同一带噪观测序列上验证突变识别/滤波重置或模型改进，再讨论闭环是否有收益；本次没有实施或证明任何改进。','',
'新增诊断MATLAB/Python脚本、本目录重放误差数组与报告；原模型、原24条轨迹、原统计和原源码不变。分支exp/noisy-bearing-kf-pilot-v1，HEAD f8b73cf；未提交、推送、合并。回退只撤销新增诊断实现，诊断输出留作历史记录。']
(root/'RMSE原因诊断.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
print(json.dumps({'status':verify['status'],'maxRMSEReplayError':audit['maxArchivedRMSEError'],
 'maxIndependentDifference':max(c['maxDifference'] for c in checks),'nearPassExample':verify['nearPassExample']},ensure_ascii=False,indent=2))
