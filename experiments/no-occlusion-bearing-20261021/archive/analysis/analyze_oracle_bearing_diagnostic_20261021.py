"""Read-only upstream reuse, independent trajectory metrics, paired plots/report."""
import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, r'C:\Users\Administrator\AppData\Local\Temp\noisy-bearing-kf-pydeps')
import numpy as np
from scipy.io import loadmat
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont(r'C:\Windows\Fonts\msyh.ttc')
plt.rcParams.update({'font.family': 'Microsoft YaHei', 'axes.unicode_minus': False,
                     'font.size': 10, 'savefig.dpi': 150})
METHODS = ['oracle', 'raw', 'selected-kf', 'imm']
NAMES = {'oracle': '真实方位（诊断）', 'raw': '带噪直接输入', 'selected-kf': 'KF', 'imm': 'IMM'}
COLORS = ['#292929', '#dc7928', '#3677b5', '#29a081']
SEEDS = list(range(6201, 6207))
CONDITIONS = [(d, s) for d in ['gaussian', 'mixture'] for s in [.01, .1]]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''): h.update(block)
    return h.hexdigest()


def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def title(d, sd):
    return f'{"高斯" if d == "gaussian" else "非高斯混合"}，观测标准差 {sd:.2f} rad'


def metrics(r):
    """Reconstruct movement order from positions, not recorded headings."""
    c = r['config']; x, y = r['trajectory']['x'], r['trajectory']['y']
    n, T, B = int(c['NAgents']), int(c['NSteps']), int(c['burnIn'])
    L, dt = float(c['arenaSize']), float(c['dt'])
    assert (n, T, B, x.shape, y.shape) == (20, 8000, 6000, (20, 8001), (20, 8001))
    dx, dy = np.diff(x, axis=1), np.diff(y, axis=1)
    dx -= L * np.round(dx / L); dy -= L * np.round(dy / L)
    vx, vy = dx / dt, dy / dt; speed = np.hypot(vx, vy)
    ux = np.divide(vx, speed, out=np.zeros_like(vx), where=speed > 1e-12)
    uy = np.divide(vy, speed, out=np.zeros_like(vy), where=speed > 1e-12)
    go = np.hypot(ux.sum(0), uy.sum(0)) / n
    px, py = x[:, B+1:].T, y[:, B+1:].T
    rx, ry = px[:, :, None] - px[:, None, :], py[:, :, None] - py[:, None, :]
    rx -= L*np.round(rx/L); ry -= L*np.round(ry/L)
    dist = np.hypot(rx, ry); ii, jj = np.triu_indices(n, 1)
    pair = dist[:, ii, jj].mean(1) / L
    near = np.argsort(dist, axis=2)[:, :, :6]
    lx = np.take_along_axis(ux[:, B:].T[:, None, :], near, axis=2).sum(2) / 6
    ly = np.take_along_axis(uy[:, B:].T[:, None, :], near, axis=2).sum(2) / 6
    lo = np.hypot(lx, ly).mean(1)
    ts = r['timeSeries']
    errors = {'GO': float(np.max(np.abs(go - ts['globalOrder']))),
              'LO': float(np.max(np.abs(lo - ts['localOrder'][B:]))),
              'distance': float(np.max(np.abs(pair - ts['meanPairDistanceNorm'][B:]))),
              'speed': float(np.max(np.abs(speed.mean(0) - ts['meanSpeed'])))}
    assert max(errors.values()) < 1e-9, errors
    row = {'seed': int(c['seed']), 'meanGO': float(go[B:].mean()), 'meanLO': float(lo.mean()),
           'meanPairDistanceNorm': float(pair.mean()), 'meanSpeed': float(speed[:, B:].mean()),
           'usedBearingRMSE': float(np.sqrt(np.mean(ts['usedBearingRMSE'][B:]**2))),
           'inputTruthRMS': float(np.sqrt(np.mean(ts['inputTruthRMS'][B:]**2)))}
    for key, field in [('meanGO', 'globalOrder'), ('meanLO', 'localOrder'),
                       ('meanPairDistanceNorm', 'meanPairDistanceNorm'), ('meanSpeed', 'meanSpeed')]:
        assert abs(row[key] - r['summary'][field]) < 1e-9
    return row, errors, go


def oracle_truth_check(r):
    c, p = r['config'], r['pilot']; L = float(c['arenaSize']); n = 20
    x, y = r['trajectory']['x'][:, :-1].T, r['trajectory']['y'][:, :-1].T
    edges = np.asarray(p['linearEdges']).astype(int)-1; rec, nei = edges % n, edges // n
    dx, dy = x[:, nei]-x[:, rec], y[:, nei]-y[:, rec]
    dx -= L*np.round(dx/L); dy -= L*np.round(dy/L)
    truth = np.arctan2(dy, dx) % (2*np.pi)
    ids = np.asarray(p['selectedHistory']).astype(int).transpose(2, 0, 1)-1
    chosen = np.zeros((8000, n, n), bool)
    chosen[np.arange(8000)[:, None, None], np.arange(n)[None, :, None], ids] = True
    fresh = chosen[:, rec, nei]; obs = p['observedBearingHistory']; used = p['usedBearingEvaluation']
    assert np.array_equal(np.isfinite(obs), fresh)
    assert np.array_equal(np.isfinite(used), fresh[6000:])
    assert np.all(fresh.sum(1) == 300)
    delta = obs[fresh]-truth[fresh]
    err = float(np.max(np.abs(np.arctan2(np.sin(delta), np.cos(delta)))))
    assert err < 1e-12 and np.array_equal(obs[6000:], used, equal_nan=True)
    assert np.size(p['estimator']) == 0
    assert np.all(r['timeSeries']['inputTruthRMS'] == 0)
    for s in p['samples']:
        assert np.array_equal(s['input'], s['truthInput'])
    return {'maxTrueBearingDifference': err, 'selectedEdgeSteps': int(fresh[6000:].sum()),
            'allArchivedInputSamplesExact': True, 'noFilterState': True}


def pair_stat(rows, d, sd, baseline):
    key = {(r['seed'], r['method']): r for r in rows
           if r['distribution'] == d and r['measurementStd'] == sd}
    a = np.array([key[s, 'oracle']['meanGO'] for s in SEEDS])
    b = np.array([key[s, baseline]['meanGO'] for s in SEEDS]); delta = a-b
    mean, se = float(delta.mean()), float(delta.std(ddof=1)/np.sqrt(6))
    margin = float(student_t.ppf(.975, 5))*se
    p = float(2*student_t.sf(abs(mean/se), 5)) if se else (1. if mean == 0 else 0.)
    return {'distribution': d, 'measurementStd': sd, 'baseline': baseline, 'nSeeds': 6,
            'oracleMeanGO': float(a.mean()), 'baselineMeanGO': float(b.mean()),
            'oracleMinusBaseline': mean, 'CI95Low': mean-margin, 'CI95High': mean+margin,
            'pUnadjusted': p, 'oracleHigherSeeds': int(np.sum(delta > 0)), 'differences': delta.tolist()}


def plot_tracks(quartet, d, sd, seed, out, full=False):
    fig, axes = plt.subplots(2, 2, figsize=(10.7, 10.5), sharex=True, sharey=True)
    colors = plt.get_cmap('tab20')(np.arange(20)); start = 0 if full else 6000
    for ax, method in zip(axes.ravel(), METHODS):
        r = quartet[method]; x, y = r['trajectory']['x']/1000, r['trajectory']['y']/1000
        for agent in range(20):
            xx, yy = x[agent, start:].copy(), y[agent, start:].copy()
            jumps = np.flatnonzero((np.abs(np.diff(xx)) > .5) | (np.abs(np.diff(yy)) > .5))+1
            xx[jumps] = np.nan; yy[jumps] = np.nan
            ax.plot(xx, yy, color=colors[agent], lw=.75, alpha=.8)
            ax.scatter(x[agent, start], y[agent, start], s=16, facecolors='none', edgecolors=[colors[agent]], lw=.7)
            ax.scatter(x[agent, -1], y[agent, -1], s=22, color=colors[agent], marker='^')
        vx = x[:, -1]-x[:, -2]; vy = y[:, -1]-y[:, -2]
        vx -= np.round(vx); vy -= np.round(vy); mag = np.hypot(vx, vy)
        ax.quiver(x[:, -1], y[:, -1], vx/mag, vy/mag, color=colors,
                  angles='xy', scale_units='xy', scale=35, width=.003)
        go = np.mean(r['timeSeries']['globalOrder'][6000:])
        dd = np.mean(r['timeSeries']['meanPairDistanceNorm'][6000:])
        ax.set_title(f'{NAMES[method]}\n评价窗平均 GO={go:.3f}，对距/L={dd:.3f}')
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel='x / L', ylabel='y / L')
        ax.set_aspect('equal'); ax.grid(alpha=.15)
    window = '全时段 0–2400' if full else '评价时段 1800–2400'
    fig.suptitle(f'{title(d, sd)} · seed {seed}\n{window}（模型时间），每组20个体', fontsize=14)
    fig.text(.5, .015, '同色表示同一个体；空心圆为窗口起点，三角与箭头为终点和方向；周期边界跳线已断开。', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .035, 1, .93))
    suffix = 'full' if full else 'evaluation'
    fig.savefig(out/f'{d}-std-{round(sd*1000):03d}-seed-{seed}-{suffix}.png'); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run', type=Path, required=True); args = ap.parse_args()
    root = args.run.resolve(); meta = json.loads((root/'metadata.json').read_text(encoding='utf-8'))
    assert json.loads((root/'completion.json').read_text())['status'] == 'COMPLETED'
    upstream = Path(meta['upstream']); out = root/'analysis'; assert not out.exists(); out.mkdir()
    tracks = out/'trajectories'; tracks.mkdir()
    oldaudit = json.loads((upstream/'independent-analysis'/'audit.json').read_text(encoding='utf-8'))
    # Upstream audit key is inspected rather than silently accepting changed data.
    old_hashes = oldaudit.get('rawFileSHA256', oldaudit.get('fileHashes', oldaudit.get('rawFileHashes', oldaudit.get('sha256', {}))))
    assert old_hashes, f'Missing upstream raw-file hashes: {list(oldaudit)}'
    oracle = {}; oracle_rows = {}; checks = []; maxerr = {}; source_hashes = {}; all_rows = []; series = {}
    for seed in SEEDS:
        path = root/f'oracle-seed-{seed}.mat'; source_hashes[path.name] = sha(path)
        r = loadmat(path, simplify_cells=True)['oracle']; oracle[seed] = r
        row, errors, go = metrics(r); oracle_rows[seed] = row
        checks.append({'seed': seed, 'kind': 'oracle', 'metrics': errors, 'truth': oracle_truth_check(r)})
        for k, v in errors.items(): maxerr[k] = max(maxerr.get(k, 0.), v)
        series[seed, 'oracle'] = go
    with (upstream/'summary.csv').open(encoding='utf-8-sig', newline='') as f:
        archived = list(csv.DictReader(f))
    lookup = {(int(r['seed']), r['distribution'], float(r['measurementStd']), r['method']): r for r in archived}
    assert len(lookup) == len(archived) == 72
    for d, sd in CONDITIONS:
        for seed in SEEDS:
            path = upstream/f'{d}-seed-{seed}-std-{round(sd*1000):03d}-triple.mat'
            digest = sha(path); assert digest == old_hashes[path.name], path
            source_hashes[path.name] = digest
            triple = loadmat(path, simplify_cells=True)['triple']; quartet = {'oracle': oracle[seed]}
            orow = {**oracle_rows[seed], 'distribution': d, 'measurementStd': sd, 'method': 'oracle'}
            all_rows.append(orow)
            for label, method in [('raw', 'raw'), ('kf', 'selected-kf'), ('imm', 'imm')]:
                r = triple[label]; quartet[method] = r
                for field in ['initialPosition', 'initialHeadings', 'initialNeuralStateUsed']:
                    assert np.array_equal(oracle[seed][field], r[field])
                for field in ['selectedHistory', 'finalNoiseStreamState', 'finalMixtureStreamState']:
                    assert np.array_equal(oracle[seed]['pilot'][field], r['pilot'][field])
                for a, b in zip(oracle[seed]['pilot']['samples'], r['pilot']['samples']):
                    assert np.array_equal(a['standardNormal'], b['standardNormal'])
                row, errors, go = metrics(r); row.update(distribution=d, measurementStd=sd, method=method)
                old = lookup[seed, d, sd, method]
                for field in ['meanGO', 'meanLO', 'meanPairDistanceNorm', 'meanSpeed', 'usedBearingRMSE', 'inputTruthRMS']:
                    assert abs(row[field]-float(old[field])) < 1e-9
                all_rows.append(row); series[seed, d, sd, method] = go
                checks.append({'seed': seed, 'distribution': d, 'std': sd, 'method': method, 'metrics': errors})
                for k, v in errors.items(): maxerr[k] = max(maxerr.get(k, 0.), v)
            plot_tracks(quartet, d, sd, seed, tracks)
            plot_tracks(quartet, d, sd, seed, tracks, full=True)
            print(f'PAIRED {d} std={sd} seed={seed}: oracle GO={orow["meanGO"]:.6f}', flush=True)
            del triple, quartet, r
    assert len(all_rows) == len({(r['seed'], r['distribution'], r['measurementStd'], r['method']) for r in all_rows}) == 96
    write_csv(out/'independently-recomputed.csv', all_rows)
    primary = [pair_stat(all_rows, d, sd, 'raw') for d, sd in CONDITIONS]
    order = np.argsort([r['pUnadjusted'] for r in primary]); previous = 0.
    for rank, index in enumerate(order):
        adjusted = min(1., max(previous, (4-rank)*primary[index]['pUnadjusted']))
        primary[index]['pHolmFour'] = adjusted; previous = adjusted
    secondary = [pair_stat(all_rows, d, sd, m) for d, sd in CONDITIONS for m in ['selected-kf', 'imm']]
    (out/'paired-statistics.json').write_text(json.dumps({'primary': primary, 'secondaryExploratory': secondary}, ensure_ascii=False, indent=2), encoding='utf-8')
    write_csv(out/'primary-paired-statistics.csv', [{k: v for k, v in r.items() if k != 'differences'} for r in primary])
    means = []
    for d, sd in CONDITIONS:
        for m in METHODS:
            group = [r for r in all_rows if r['distribution'] == d and r['measurementStd'] == sd and r['method'] == m]
            means.append({'distribution': d, 'measurementStd': sd, 'method': m, 'nSeeds': 6,
                          **{key: float(np.mean([r[key] for r in group])) for key in
                             ['meanGO', 'meanLO', 'meanPairDistanceNorm', 'meanSpeed', 'usedBearingRMSE', 'inputTruthRMS']}})
    write_csv(out/'condition-means.csv', means)
    (out/'condition-means.json').write_text(json.dumps(means, ensure_ascii=False, indent=2), encoding='utf-8')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), sharey=True)
    for ax, (d, sd) in zip(axes.ravel(), CONDITIONS):
        for seed in SEEDS:
            values = [next(r['meanGO'] for r in all_rows if (r['seed'], r['distribution'], r['measurementStd'], r['method']) == (seed, d, sd, m)) for m in METHODS]
            ax.plot(range(4), values, '-o', color=plt.get_cmap('tab10')(seed-6201), alpha=.65, lw=1, ms=4, label=str(seed))
        ax.plot(range(4), [next(r['meanGO'] for r in means if (r['distribution'], r['measurementStd'], r['method']) == (d, sd, m)) for m in METHODS], 'kD-', lw=2.2, ms=6, label='6种子均值')
        ax.set_xticks(range(4), [NAMES[m] for m in METHODS], fontsize=9)
        ax.set_title(title(d, sd)); ax.set_ylim(0, 1.03); ax.grid(alpha=.2); ax.set_ylabel('评价窗口平均 GO')
    axes[0, 0].legend(fontsize=8, ncol=2)
    fig.suptitle('全部6种子的配对结果；同色连线为同一种子\n四面板复用相同6条真实方位轨迹，独立样本仍为每条件 n=6', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .92)); fig.savefig(out/'paired-go-all-seeds.png'); plt.close(fig)
    for mean_mode in [False, True]:
        fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True, sharey=True)
        time = np.arange(1, 8001)*.3
        for ax, (d, sd) in zip(axes, CONDITIONS):
            for method, color in zip(METHODS, COLORS):
                curves = [series[s, 'oracle'] if method == 'oracle' else series[s, d, sd, method] for s in SEEDS]
                values = np.mean(curves, axis=0) if mean_mode else curves[0]
                ax.plot(time, values, color=color, lw=1.2, alpha=.9, label=NAMES[method])
            ax.axvspan(1800, 2400, color='grey', alpha=.1); ax.set_title(title(d, sd), loc='left', fontsize=10)
            ax.set_ylabel('GO'); ax.set_ylim(0, 1.03); ax.grid(alpha=.15)
        axes[0].legend(ncol=4, fontsize=9); axes[-1].set_xlabel('模型时间（灰区为评价窗口）')
        fig.suptitle('6种子逐时刻均值（未平滑）' if mean_mode else '预先指定示例 seed6201：原始 GO 时间序列（未平滑）')
        fig.tight_layout(rect=(0, 0, 1, .96)); fig.savefig(out/('mean-go-timeseries.png' if mean_mode else 'seed6201-go-timeseries.png')); plt.close(fig)
    audit = {'status': 'PASS', 'newUniqueTrajectories': 6, 'reusedUniqueTrajectories': 72,
             'tableRows': 96, 'pairedQuartets': 24, 'independentSeedsPerCondition': 6,
             'maxMetricErrors': maxerr, 'checks': checks, 'rawFileSHA256': source_hashes,
             'upstreamDataUnchangedHashCheck': True, 'oraclePairingAllInitialStatesAllStepsSelections': True,
             'newFilterStateReplay': False, 'trajectoryPanels': 48, 'figureCount': 51,
             'exampleSeedPrespecified': 6201, 'evidence': meta['evidence']}
    (out/'audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    fallacies = {
        'Simpson': 'Four conditions retained separately; oracle is repeated across conditions, not pooled as 24 independent controls.',
        'Ecological': 'Inference unit is seed-level simulated group, not individual animal or edge.',
        'Berkson': 'All six previously specified seeds retained; limited fixed-parameter population.',
        'Collider': 'No adjustment on post-treatment GO, geometry or error.',
        'BaseRate': 'No diagnostic classification claim or success-rate endpoint.',
        'RegressionToMean': 'No extreme-seed selection; example seed fixed before oracle results.',
        'Survivorship': 'All 78 unique scientific trajectories included; software/regression runs are not independent scientific repeats.',
        'LookElsewhere': 'All four primary contrasts reported, Holm adjusted; eight other GO contrasts explicitly exploratory.',
        'ForkingPaths': 'Plan saved before experiment; prior viewed seeds make this an exploratory extension.',
        'CorrelationCausation': 'Controlled total input intervention does not identify neural denoising or error mediation.',
        'ReverseCausality': 'Closed-loop motion affects later bearings; lower RMSE is not assumed to cause higher GO.'}
    (out/'methodology-audit.json').write_text(json.dumps({
        'coverage': '11/11 checked', 'interpretation': 'CAUTION: n=6 paired t intervals assume approximately normal seed differences; no equivalence test',
        'fallacyChecks': fallacies}, ensure_ascii=False, indent=2), encoding='utf-8')
    lines = ['# 真实方位输入诊断实验报告', '',
             'Material Passport：academic-research-suite / experiment-agent；2026-10-21。', '',
             '证据等级：对已查看种子的探索性诊断扩展；不是新的独立确认性验证。', '',
             '新增6条真实方位轨迹，配对复用此前72条带噪直接输入/KF/IMM轨迹。四条件共24个四元组、96条条件记录，实际78条唯一轨迹；真实方位对照跨条件重复，每条件独立单位始终为6个种子。', '',
             'N20/M100、逐步随机k15、original等权强度.24、核宽.4、beta1000、speedScale .05、dt .3、周期L1000；8000步、评价6001:8000。无新增遮挡或领航。KF/IMM参数与上游冻结值一致，无重新选参。', '',
             'oracle在自身当前几何下使用所选邻居的真实方位，不调用或更新KF/IMM；实际估计器只接收带噪观测。四组共有初态、全部选邻和基础随机数流；反馈分叉后位置、真实方位和观测角度不同。', '',
             '## 全种子结果', '', '| 条件 | 真实方位 GO | 带噪直接 GO | KF GO | IMM GO |', '|---|---:|---:|---:|---:|']
    for d, sd in CONDITIONS:
        vals = [next(r['meanGO'] for r in means if (r['distribution'], r['measurementStd'], r['method']) == (d, sd, m)) for m in METHODS]
        lines.append(f'| {title(d, sd)} | '+ ' | '.join(f'{v:.6f}' for v in vals)+' |')
    lines += ['', '主要比较为真实方位−带噪直接输入的平均GO。普通单项paired t 95%区间，Holm校正4项p；区间跨0不能证明等效。', '',
              '| 条件 | 差值 | 95%区间 | Holm p | 真实方位更高种子数 |', '|---|---:|---|---:|---:|']
    for r in primary:
        lines.append(f'| {title(r["distribution"], r["measurementStd"])} | {r["oracleMinusBaseline"]:+.6f} | [{r["CI95Low"]:+.6f}, {r["CI95High"]:+.6f}] | {r["pHolmFour"]:.4f} | {r["oracleHigherSeeds"]}/6 |')
    lines += ['', '## 次要指标', '', '| 条件 | 方法 | LO | 平均对距/L | 方位RMSE rad |', '|---|---|---:|---:|---:|']
    for r in means:
        lines.append(f'| {title(r["distribution"], r["measurementStd"])} | {NAMES[r["method"]]} | {r["meanLO"]:.6f} | {r["meanPairDistanceNorm"]:.6f} | {r["usedBearingRMSE"]:.6f} |')
    lines += ['', '## 轨迹与时间序列', '',
              '预先固定seed6201为示例；不是按GO选择。全部6种子×4条件均保存评价窗口与全时段四方法并排轨迹，共48张；同一个体颜色一致，空心圆为起点，三角/箭头为终点/方向。周期边界跳线断开。灰区评价窗口与轨迹显示窗口一致。', '',
              '![全种子配对GO](paired-go-all-seeds.png)', '',
              '![高斯 .10，seed6201](trajectories/gaussian-std-100-seed-6201-evaluation.png)', '',
              '![非高斯 .10，seed6201](trajectories/mixture-std-100-seed-6201-evaluation.png)', '',
              '![原始GO时间序列](seed6201-go-timeseries.png)', '',
              '## 解释边界', '',
              '真实方位是零观测误差的诊断参考，不是GO数学上界；神经场和运动反馈可使噪声或时间滤波改变状态转换，因此不能将方法均值排序直接解释为神经环去噪机制。未做等效性/非劣性检验，当前宽区间不能证明“原模型已完全抗扰”或“KF永远无用”。未扫描连续遮挡、更多噪声强度或神经时间尺度，也不能回答中间是否存在滤波有效区间。', '',
              '平均GO衡量方向一致性；需同时看LO和对距。高GO不是持续稳定、凝聚或安全保证。所有统计以种子为单位，不能把时间点和邻居边当独立重复。oracle−KF/IMM的探索性配对结果见paired-statistics.json。n=6的t区间依赖种子差值近似正态，未证实该假设；结论谨慎。方法审计11/11项见methodology-audit.json。', '',
              '## 验证与文件保护', '',
              '短时oracle与零噪旧raw轨迹精确一致；oracle无滤波状态，真实方位及输入误差为0；未使用的噪声分布变化不改变oracle轨迹。完整8000步重跑6201 Gaussian .10 raw与上游精确一致。KF/IMM真值隔离契约通过。', '',
              f'独立从全部78条位置重建GO、LO、对距和速度；最大误差：{json.dumps(maxerr, ensure_ascii=False)}。新增6条对照全部方位独立重建、全部选邻配对；上游24个MAT的SHA256与原审计一致。本次没有重新重放所有KF/IMM内部状态，上游验证范围仍见上游审计。', '',
              '原始模型、所有旧结果、旧滤波器、论文和既有tracked修改只读。本轮只新增代码/结果。分支exp/oracle-bearing-diagnostic-v1，HEAD f8b73cf，无commit/push/merge。建议保存实验分支，代码审查后再集成；回退只撤销新增实现，保留全部结果。']
    (out/'真实方位诊断实验报告.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    (root/'README.md').write_text('# 真实方位诊断对照\n\n完整结果见 [实验报告](analysis/真实方位诊断实验报告.md)。本轮新增6条真实方位轨迹，复用上游72条轨迹；所有旧数据只读。\n\nmetadata.json、source-snapshot/、software-test-audit.json、oracle-summary.csv、completion.json记录参数与测试。analysis/记录独立复核、96条条件表、种子配对统计与51张图。\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'means': means, 'primary': primary}, ensure_ascii=False, indent=2))


if __name__ == '__main__': main()
