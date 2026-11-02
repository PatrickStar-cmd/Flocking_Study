"""Compile the four planned decisions without turning development into validation."""
import json,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd


def main(parent):
    target=parent/'P21_FINAL_REPORT.md';assert not target.exists()
    dev=parent/'development-r01';val=parent/'validation-r01'
    frozen=json.loads((dev/'FROZEN_DECISION.json').read_text(encoding='utf-8'))
    seed=pd.read_csv(dev/'analysis/seed_metrics.csv');seed=seed[seed.scenario=='combined']
    sweep=pd.read_csv(dev/'analysis/q_sweep.csv')
    nis=pd.read_csv(parent/'g0-r02/nis_mass_contribution.csv')
    synthetic=pd.read_csv(parent/'g1/synthetic_summary.csv')
    cv=nis[nis.method=='KF-CV'].set_index('geometry')
    fraction=float(cv.loc['switch_event','fraction_edge_times']+cv.loc['switch_near2','fraction_edge_times'])
    mass=float(cv.loc['switch_event','fraction_NIS_sum']+cv.loc['switch_near2','fraction_NIS_sum'])
    clean_nis=float(cv.loc['G-clean','sum_NIS']/cv.loc['G-clean','n_edge_times'])
    rows=[]
    for geometry in ['all','G-clean','G-switch','G-close','clean-slow','clean-mid','clean-fast']:
        for D in [20,120]:
            f=seed[(seed.geometry==geometry)&(seed.D==D)]
            h=f[f.method=='Hold'].set_index('seed');c=f[f.method=='CV_q0.0001'].set_index('seed');g=f[f.method=='GM_default'].set_index('seed')
            ids=h.index.intersection(c.index)
            rows.append(f"| {geometry} | {D} | {len(ids)} | {h.rmse.mean():.5f} | {c.rmse.mean():.5f} | {g.rmse.mean():.5f} | {(c.loc[ids].rmse-h.loc[ids].rmse).mean():+.5f} |")
    qrows=[]
    for _,r in sweep.iterrows():qrows.append(f"| {r.q:g} | {r.objective_fast_D120_RMSE:.5f} | {r.slow_D20_ratio_to_Hold:.3f} | {r.age_coverage_min:.3%}–{r.age_coverage_max:.3%} | {r.admissible} |")
    a=synthetic[(synthetic.config_id=='A1')&(synthetic.D==20)].groupby('method').circular_rmse.mean()
    if frozen['q_star'] is None:
        validation_section='开发阶段没有候选同时满足预先冻结的支持数量、slow短期保护与严重失配排除条件。因此q_star=null；按照运行前规则关闭验证门槛，未生成、未读取任何预留验证seed。**不能给出独立验证的ΔRMSE或CI，也没有V1/V2验证收益结果。** 停止并非软件失败，不得将训练集最小RMSE强行命名为可信q*。'
        next_decision='本轮不支持把调Q后的KF升级为正式方案或进入P3。先解释开发表中短期保护与条件覆盖的冲突；若另行提出几何/运动模型或确认性默认算法检验，必须新协议、新版本，不能放宽本轮门槛或回写结果。'
    else:
        audit=json.loads((val/'analysis/validation_audit.json').read_text(encoding='utf-8'))
        p=audit['primary']
        validation_section=f"冻结全局q*={frozen['q_star']:g}后才生成{frozen['validation_n']}个新物理seed，每个base/turn各一条。主比较clean-fast D120有{p['n_seeds']}个支持seed，平均配对差{p['mean']:+.5f}rad，95%seed-cluster bootstrap区间[{p['ci_lower']:+.5f},{p['ci_upper']:+.5f}]，相对改善{p['mean_relative_improvement']:.2%}。保护性slow D20均值比={audit['slow_D20_ratio']:.3f}。完整seed点和所有次要条件保留在validation_contrasts.csv/validation_summary.csv，不能用边或窗口放大n。"
        next_decision='进入P3的意见必须综合主效应是否达到0.02rad、区间、短期保护、条件覆盖/NIS和尾部反例；本轮仅归档开环验证，不执行P3，也不以总体覆盖率接近95%宣布可靠。'
    report=f'''# P2.1 几何解耦与KF校准实验报告

## Material Passport

计划：用户提供的P2.1 Protocol v1.0，原文件按原字节保存于experiments/kf-p21-20261102/PROTOCOL.md；运行前未指定细则见RUN_RULES.md。原物理/滤波基线Hao_L@81e0080。旧5seed只用于G0审计，新开发12seed为2026111001–2026111012，验证完全隔离。本轮仅开环；P3/P4未执行。

## 四个必须回答的判断

1. **几何贡献**：旧P2观测窗口中switch当步及±2步只占{fraction:.3%}边时刻，却贡献{mass:.2%}的CV NIS总量；G-clean逐点时刻平均CV NIS约{clean_nis:.3f}。这反驳“高总体NIS必然代表所有连续运动都失配”的解释。该统计含重叠窗口、5个复用物理seed，是条件性描述，不是独立因果干预。close、near-half-box、aliasing及其支持数量全部保留。
2. **连续快变化优势**：下表给出12个新开发seed在清洁continuous-fast中的默认CV/Hold对照。开发层的三分位切点为{frozen['rate_cutpoints'][0]:.8f}/{frozen['rate_cutpoints'][1]:.8f}rad/模型时间，仅从开发真值计算一次。它是相对运动分层，不是在线可用门控。
3. **Q校准**：五个候选全部报告，只允许一个全局q，禁止逐seed/边/D分别挑选。冻结状态为`{frozen['status']}`，q*=`{frozen['q_star']}`；细节与候选表见下。
4. **下一阶段**：{next_decision}

## 执行范围与前置验收

A门槛复核144个已登记P2文件SHA与全部270行RMSE/NIS。G0独立位置代数核对2814000条逐边时刻标签，全部窗口分类通过；G0-r02为正式版本，g0初版因图单位标签有误仅作为QA历史保留。G-clean/G-switch/G-close/G-ambiguous严格互斥；缺层输出NaN，不填0，不删困难样本。

G1使用20条固定随机流、8个解析配置、20/120步共320条件，三方法960行独立复算通过。D为无碰撞近距离掠过，E为连续位置下minimum-image切换，所有方法共享噪声和掩码。常方位A的20步平均RMSE：Hold={a['Hold']:.5f}rad，CV={a['KF-CV']:.5f}rad，GM={a['KF-GM']:.5f}rad；普通观测噪声在速度估计后被外推也能形成误差，不能只把所有短期劣势归于几何跳变。合成结果不与物理独立seed混合。

新开发24条物理参考轨迹全部Fresh，参数N20/神经元100/k15、1000周期边界、dt=.3、观测std=.01、R=1e-4沿用P2；转向刺激1401步起、全20个体、幅度.06、速率.02、400步，前缀完整状态/随机源匹配。全部25起点1201:20:1681、两种D、两类来源保留；25×24×2=1200个候选源窗口，7个不同方法配置，总8400个方法/窗口条件。独立推断单位始终只有12个开发物理seed。

## 开发集默认方法结果

每窗口先在符合几何标签的边和全部遮挡时刻求圆周RMSE，再在每seed内等权平均可用候选窗口，最后跨seed汇总。两个源的候选集合均预先定义；空层不制造0。单位rad。下表不是独立验证结果。

| 几何层 | D步 | 支持seed | Hold | 默认CV | 默认GM | CV−Hold |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

## Q网格、保护与冻结

主要开发目标clean-fast D120 seed均值RMSE最小。保护约束clean-slow D20的均值比不超过1.10；覆盖90%–98%是参考，不是严格高斯证明。RUN_RULES在任何新结果打开前将严重失配工作定义为主/保护年龄段跨seed平均覆盖<80%或>99.5%，非有限/非PSD同样不合格。支持门槛和近似等价1%规则也在运行前写定，不根据本轮结果放宽。

| CV q | clean-fast D120 RMSE | slow D20 / Hold | 年龄段覆盖范围 | 可接受 |
|---:|---:|---:|---|---|
{chr(10).join(qrows)}

GM仍为q=1e-5,tau10的原整体配置，不把不同Q参数化的数值直接解释成模型结构优越。门限age>120模型时间或P>1仅作离线部署诊断，主RMSE保留全部预测；首次截断、比例、连续截断长度、角速度及仅在清洁差分上的误差均已保存。

## 独立验证状态

{validation_section}

FROZEN_DECISION.json及SHA固定了Q、切点、源代码、主响应、关注差异.02rad、样本量公式与seed列表。任何后续改变须新版本；不能把本轮停止/负结果覆盖为新方法成功。

## 数据、数值与复核

全部5候选CV和GM与未修改P2实现及原MATLAB逐步对照：6模型×3边×2200步，最大方位差约1.06e-13rad，P差3.55e-15，NIS绝对差5.82e-11。10项P2.1软件测试通过；缺失/未来污染、真值隔离、MATLAB取整、PSD、普通坐标回绕不误报、140边形状和age单位均核验。共享前缀暖启动和压缩存储相对完整重跑通过逐位恢复测试。

raw/input保留物理位置与方位观测。raw/replay采用共享全前缀+float64位模式XOR，精确保留所有25窗口的遮挡/400步恢复预测，前200步来自共享fresh。engine.recover_branch可逐位还原，不是只存汇总或近似压缩。独立审计从位置重建标签、速率和开发切点，并从原始预测复算每个非空窗口RMSE。SHA、真实样本支持、所有候选/空层、NIS分位尾部、P分箱、年龄覆盖、极端误差及门限记录完整保存。

G0审计、G1机理、开发选择和可能的独立验证证据分开。KF能否被神经场利用：**本轮不回答**。KF能否提高GO：**本轮不回答**。不声称二维位置可观测性、闭环稳定性或安全保证。

## 文件保护、发布与回退

仅新增P2.1目录；原始allocentric-flocking-original、effective-social-input既有模型、P1/P2脚本和历史结果未改。新结果位于kf-p21-v1-20261102-r01，分阶段目录与QA记录均保留。当前任务分支{ subprocess.check_output(['git','branch','--show-current'],text=True).strip() }。发布前逐文件核验，远端写入遵照用户授权；本报告不把本地完成写成已上传。公共回退使用对应提交git revert，不强推、不改main。
'''
    target.write_text(report,encoding='utf-8')
    print('Wrote',target)


if __name__=='__main__':main(Path(sys.argv[1]))
