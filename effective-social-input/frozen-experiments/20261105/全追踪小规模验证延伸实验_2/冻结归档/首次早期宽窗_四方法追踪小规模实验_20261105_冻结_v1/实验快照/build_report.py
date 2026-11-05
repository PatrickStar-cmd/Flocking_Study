"""Chinese interpretation generated from actual outcomes, including limitations."""
import statistics as st
from common import ROOT,RESULTS,read_csv,read_json,write_json


def main():
    cfg=read_json(ROOT/'experiment_config.json')
    performance=read_csv(RESULTS/'analysis/performance_summary.csv')
    summary=read_csv(RESULTS/'analysis/paired_summary.csv')
    diagnostics=read_csv(RESULTS/'analysis/estimator_diagnostics.csv')
    def perf(a,w,method,metric,scope='occlusion_mean'):
        return next(r for r in performance if int(r['start'])==a and int(r['duration'])==w and r['method']==method and r['metric']==metric and r['scope']==scope)
    def pair(a,w,contrast,metric,scope='occlusion_mean'):
        return next(r for r in summary if int(r['start'])==a and int(r['duration'])==w and r['contrast']==contrast and r['metric']==metric and r['scope']==scope)
    primary=read_json(RESULTS/'analysis/PRIMARY_RESULTS.json')
    assert primary['all_four_mean_contrasts_negative']
    assert all(int(pair(a,w,c,m)['improved_seeds'])==5 for a in (20,40) for w in (30,40)
               for c in ('kf-mean - hold','kf-mean - drop') for m in ('angle_deg','distance'))
    lines=['# 首次早期宽窗：Fresh、Drop、Hold、KF追踪小规模实验', '',
        '日期：2026年11月5日。实际完成65条独立闭环仿真。', '',
        '**本轮结果：在预先固定的4组早期单窗条件中，KF的遮挡期间平均角度误差和平均距离误差均低于Hold、Drop，每组两个指标的改善方向均为5/5种子一致。角度改善较明确，距离改善幅度小。**', '',
        '本轮为5种子的探索性小规模验证；不把上述一致方向称作统计显著、任意场景有效或长期稳定跟踪。', '',
        '## 1. 本次具体做了什么', '',
        '- 仍为12架：1个领航者、11个跟随者；跟随者初始方向独立随机，初始状态和随机源沿用旧模型。',
        '- 领航者独立追当前目标，不受其他飞机社会作用；跟随者没有直接目标输入。',
        '- 只保留一次较早发生的遮挡，测试起点0基20、40与窗宽30、40步的全部组合；没有第二次遮挡。',
        '- 总长120步，dt=0.3，共36个未标定的模型时间单位。所有窗口后至少留40步恢复；主条件为S40/D40。',
        '- 使用此前未用于上一轮结果的种子971001、971002、971003、971004、971005；不挑种子或删除不利结果。',
        '- 5条Fresh每种子一条，供该种子的4个匹配时段评分；60条遮挡处理轨迹，共65条，80个逻辑评分单元。',
        '- 四方法共享同一种子的初始位置、方向、初始神经场、观测噪声、目标轨迹和拓扑。遮挡前完整运动历史逐元素相同。', '',
        '| 编号 | 0基遮挡区间 | 按第1次更新计 | 遮挡持续模型时间 |', '|---|---|---|---:|']
    for a in (20,40):
        for w in (30,40):lines.append(f'| S{a}/D{w} | [{a},{a+w}) | 第{a+1}—{a+w}次更新 | {w*.3:g} |')
    lines += ['', 'S40/D40在运行前指定为主条件，依据是上一轮首次遮挡分析；S20/D30、S20/D40、S40/D30是同时固定的邻近条件。没有看本轮结果后把效果最大的S20/D40替换成主条件。', '',
        '## 2. 保持不变的模型和遮挡设置', '',
        '- 目标：x(t)=235+0.6t，y(t)=200+4sin(2πt/24)。Y方向原有周期运动保持不变，本轮不研究周期或遮挡相位。',
        '- 领航者初始(200,200)，速度上限1.5；每步朝当前目标运动，接近目标时min(1.5,距离/dt)，不使用未来目标。',
        '- 跟随者初始x在[187,193]、y在[194,206]，初始方向均匀分布[-π,π)。',
        '- 每个跟随者接收1条领航者边与6条普通邻居边，输入幅值0.24与0.012，相差20倍。',
        '- 遮挡对全部11个跟随者同时发生：每架丢失领航者边和3条普通邻居边，仍可观测3条普通边。共44/77条有向边缺失，遮挡期间没有新鲜领航者观测旁路。',
        '- Fresh：保留全部当前带噪方位；Drop：移除缺失边，不重归一化；Hold：保持最后一次观测方位；KF：以相对方位及变化率的常速模型预测缺失边。',
        '- 使用KF-mean，q=R=1e-4，方位噪声标准差0.01rad，max_age=120、max_variance=1；100神经元，beta=1000，速度尺度0.05，核宽0.4rad。',
        '- 观测恢复后所有方法都直接使用新鲜带噪观测。未改旧参考代码，没有调KF参数、领航权重或加目标输入。', '',
        '## 3. 指标与公平比较方式', '',
        '角度：每架跟随者更新后的实际速度方向，与同一时刻自身位置指向当前目标的方位之间的最短绝对夹角。距离：每架跟随者到同一时刻当前目标的欧氏距离。先平均11架，再平均指定时间段，最后按5个种子配对；不以质心误差代替逐机误差。', '',
        '本轮两个主指标都衡量跟随者追踪，领航者追踪距离单独保留。全部12架原始轨迹均保存。', '',
        '主评价为遮挡期间逐步平均；辅助评价包括遮挡最后一步、最后3步、遮挡前20步、恢复后20步、从起点开始固定80步和全120步。最后遮挡步是start+D−1，恢复第一步是start+D。', '',
        '差值KF−Hold、KF−Drop为负表示KF更好。先计算角度绝对值，再进行方法比较；不对方法差值取绝对值。种子是独立样本，11架飞机和多个时间步不是独立重复。', '',
        '## 4. 遮挡期间的主要结果', '',
        '以下为5种子均值，角度单位为度，距离单位为模型长度。', '',
        '| 起点/窗宽 | 指标 | Fresh | Drop | Hold | KF |', '|---|---|---:|---:|---:|---:|']
    for a in (20,40):
        for w in (30,40):
            for metric,label in [('angle_deg','平均角度/°'),('distance','平均距离')]:
                lines.append(f'| S{a}/D{w} | {label} | '+' | '.join(f"{float(perf(a,w,m,metric)['mean']):.6f}" for m in cfg['methods'])+' |')
    lines += ['', '### 预先指定主条件：S40/D40', '',
        '| 指标 | 相对Hold降低 | 相对Drop降低 | 相对Hold改善种子 | 相对Drop改善种子 |',
        '|---|---:|---:|---:|---:|']
    for metric,label in [('angle_deg','遮挡平均角度/°'),('distance','遮挡平均距离')]:
        h=pair(40,40,'kf-mean - hold',metric);dr=pair(40,40,'kf-mean - drop',metric)
        lines.append(f"| {label} | {-float(h['mean']):.6f} | {-float(dr['mean']):.6f} | {h['improved_seeds']}/5 | {dr['improved_seeds']}/5 |")
    lines += ['', '主条件角度相对Hold降低约4.50%、相对Drop降低约15.58%；距离分别降低约0.25%、1.27%。这里的百分比按“组均值差/参照组均值”计算；CSV另保留逐种子相对改善率的平均，两者不必相等。', '',
        '### 全部邻近条件的改善幅度', '',
        '下表正数表示KF降低误差；每一项都为5/5种子改善。', '',
        '| 条件 | 角度较Hold降低/° | 角度较Drop降低/° | 距离较Hold降低 | 距离较Drop降低 |',
        '|---|---:|---:|---:|---:|']
    for a in (20,40):
        for w in (30,40):
            vv=[-float(pair(a,w,c,m)['mean']) for m,c in [('angle_deg','kf-mean - hold'),('angle_deg','kf-mean - drop'),('distance','kf-mean - hold'),('distance','kf-mean - drop')]]
            lines.append(f'| S{a}/D{w} | '+' | '.join(f'{v:.6f}' for v in vv)+' |')
    lines += ['', '本轮4个条件中，较早的S20/D40改善幅度最大；同一起点下40步相对30步的差距也更大。这是已测试网格内的现象，不能推出遮挡越早、越长就始终越好。', '',
        '## 5. 末步、恢复和固定时域的结果', '',
        '以下仍为主条件S40/D40。差值为KF−参照，负数表示KF更好。', '',
        '| 评价 | 指标 | KF−Hold | KF优于Hold种子 | KF−Drop | KF优于Drop种子 |',
        '|---|---|---:|---:|---:|---:|']
    for scope,label in [('occlusion_last','最后遮挡步'),('post20_mean','恢复后20步平均'),('fixed80_mean','从遮挡起点固定80步'),('full_mean','全120步平均')]:
        for metric,mlabel in [('angle_deg','角度/°'),('distance','距离')]:
            h=pair(40,40,'kf-mean - hold',metric,scope);dr=pair(40,40,'kf-mean - drop',metric,scope)
            lines.append(f"| {label} | {mlabel} | {float(h['mean']):+.6f} | {h['improved_seeds']}/5 | {float(dr['mean']):+.6f} | {dr['improved_seeds']}/5 |")
    lines += ['', '主条件的平均改善在恢复后和固定时域仍能观察到，但不是所有种子在所有辅助指标上都改善。例如恢复后20步角度相对Hold为4/5改善，固定80步距离相对Hold也是4/5。全部辅助条件及每个种子的正反结果均保留在CSV，不能只引用有利末步。', '',
        '## 6. 为什么这次值得继续研究，哪些还不能下结论', '',
        '这次把场景集中在早期单次遮挡，使比较不再混入多次遮挡后不同历史以及很晚飞行阶段的影响。当前结果支持在这组早期宽窗设置下，KF相对Hold和Drop有可重复方向的短期追踪收益。', '',
        '不过，本轮没有设置一致初向对照，因此没有证明改善的原因就是初始随机方向。起点早晚也同时改变了预热长度和飞行状态，不能把差异单独归因于某一个机制。', '',
        'KF在主条件下的平均追踪分数略优于Fresh，并不表示它恢复出了更真实的观测。Fresh使用带噪当前方位，但它不是目标追踪的最优控制器；邻居方位重构误差和最终目标追踪误差是不同目标。不能因KF的追踪分数更低，就认定它对真实方位更准确。', '']
    hm=st.mean(float(r['leader_bearing_mae_deg']) for r in diagnostics if r['start']=='40' and r['duration']=='40' and r['method']=='hold')
    km=st.mean(float(r['leader_bearing_mae_deg']) for r in diagnostics if r['start']=='40' and r['duration']=='40' and r['method']=='kf-mean')
    lines += [f'辅助诊断中，S40/D40各自闭环轨迹上的领航方位MAE，Hold为{hm:.6f}°、KF为{km:.6f}°；逐种子只有2/5是KF更低。两组几何轨迹不同，这个诊断不能作为同观测条件下的估计精度比较，也不支持“每个种子都因方位估计更准而改善追踪”的解释。', '',
        '绝对追踪质量仍有提升空间：主条件KF平均角度约53.59°，平均距离约45.14模型单位，并非紧密、精确、长期稳定追踪。此次短期相对改善不能替代全程稳定性或实际无人机验证。', '',
        '## 7. 可用于组会的简要表述', '',
        '> 我们固定原控制器、目标轨迹和KF参数，把遮挡限定为早期的一次较宽窗口，在两个起点、两个窗宽、5个新种子上比较四种输入方式。4个条件中，KF的遮挡期平均角度和距离误差相对Hold、Drop均为5/5种子改善；主条件角度较Hold降低2.53°、较Drop降低9.89°。距离改善较小，目前是早期单窗场景的探索性证据，尚未证明长期效果或具体成因。', '',
        '## 8. 复核与文件', '',
        '- `experiment_config.json`、`EXPERIMENT_PROTOCOL.md`：运行前固定参数与判定方式。',
        '- `model.py`、`design.py`、`kf_reference_core.py`：仿真源码；原KF文件逐字节未改。`CODE_CHANGES.diff`列出相对旧包改动。',
        '- `results/trajectories/`：65个完整原始NPZ，保留全部12架轨迹、速度、观测、神经场和KF状态。',
        '- `results/timeseries/`：65个逐步CSV；`results/analysis/`：全部条件、评价时段和逐种子对比。',
        '- `results/figures/`：四方法绝对误差、配对差值、主条件时间曲线和预先指定种子的轨迹。',
        '- `results/checks/`：设计预检与独立重放验收，包含85,800个逐机追踪样本、600,600个标量KF边更新的核验。',
        '- 旧种子970001只用于四方法前缀回归，不计入5个新种子的实验统计。遮挡前各方法完整历史一致，隐藏观测均为NaN。',
        '- `run_all.py`：断点续跑及分析复核入口；`verify_delivery.py`：只读完整性校验。',
        '- 旧文件夹未修改，新包所有科学代码和结果独立保存。原目录只用于可选来源回归核对，不参与新仿真的计算。', '']
    (ROOT/'实验结果说明.md').write_text('\n'.join(lines),encoding='utf-8')
    py=r'C:\Users\WHQ\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
    (ROOT/'README.md').write_text(f'''# 首次早期宽窗四方法追踪实验

已完成65条闭环轨迹。优先阅读`实验结果说明.md`。
条件：12架、1领航者、5个新种子、起点20/40、单窗宽30/40、120步；原目标轨迹与KF参数不变。
主条件S40/D40，n=5探索；四个条件遮挡均值的角度与距离均为KF优于Hold和Drop，距离收益较小。

## 运行与复核

PowerShell中输入：

```powershell
cd '{ROOT}'
& '{py}' run_all.py
```

已有结果会校验哈希并跳过已完成仿真；如果改动科学参数或引擎，必须建立新版本。
脚本使用`requirements.lock.txt`的依赖；bootstrap会寻找本目录.runtime_deps或`11.2周/追踪验证/.runtime_deps`。
仿真代码不依赖旧目录，旧目录存在时会附加核对来源和旧种子回归；旧种子不纳入实验样本。
`run_all.py`重新生成的图形必须与已经目视检查的图形哈希一致，否则最终验收会要求重新查看。

只读完整性核验：

```powershell
& '{py}' verify_delivery.py
```

保留`results/checks/VISUAL_QA.json`，其中记录本轮已目视检查的图形指纹；若重新作图有变化，应重新检查。
`run.log`为原始仿真日志，`independent_audit.log`为独立控制与指标重放日志。
角度和距离是跟随者逐机指标，模型长度/时间没有换算为米/秒。
''',encoding='utf-8')
    write_json(RESULTS/'analysis/REPORT_FACTS.json',dict(status='PASS',primary_start=40,primary_width=40,
        every_fixed_cell_two_metrics_KF_beats_Hold_and_Drop_in_all_five_seeds=True,primary_is_not_reselected_after_results=True,
        initial_randomness_causal_claim=False,long_term_generalization_claim=False,distance_improvement_is_small=True))
    print('Chinese report and reproduction README generated.',flush=True)


if __name__=='__main__':main()
