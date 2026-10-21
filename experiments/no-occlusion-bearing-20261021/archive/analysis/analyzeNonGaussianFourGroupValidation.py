"""Frozen 64-seed four-arm study: independent metrics, six contrasts and figures."""
import argparse
import csv
import gc
import importlib
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
import analyze_oracle_bearing_diagnostic_20261021 as reference
import recheck_bearing_imm_replay_pilot as matrix_reference

METHODS = ['oracle', 'raw', 'selected-kf', 'imm']
CONTRASTS = [('oracle', 'raw'), ('selected-kf', 'raw'), ('imm', 'raw'),
             ('imm', 'selected-kf'), ('oracle', 'selected-kf'), ('oracle', 'imm')]
NAMES = reference.NAMES
COLORS = reference.COLORS


def wrap(x): return np.arctan2(np.sin(x), np.cos(x))


def write_json(path, data):
    assert not path.exists(), path
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


def proof(r):
    p = r['pilot']
    return {'initialPosition': r['initialPosition'], 'initialHeadings': r['initialHeadings'],
            'initialNeuralStateUsed': r['initialNeuralStateUsed'],
            'selectedHistory': p['selectedHistory'], 'finalNoiseStreamState': p['finalNoiseStreamState'],
            'finalMixtureStreamState': p['finalMixtureStreamState'], 'evalNoise': p['evalNoise'],
            'samples': [{k: s[k] for k in ['standardNoise', 'standardNormal', 'outlierMask']} for s in p['samples']]}


def equal_proof(a, b):
    for k in ['initialPosition', 'initialHeadings', 'initialNeuralStateUsed', 'selectedHistory',
              'finalNoiseStreamState', 'finalMixtureStreamState']:
        assert np.array_equal(a[k], b[k]), k
    assert a['evalNoise'] == b['evalNoise']
    assert len(a['samples']) == len(b['samples']) == 81
    for x, y in zip(a['samples'], b['samples']):
        for k in x: assert np.array_equal(x[k], y[k]), k


def check_input(r, method, replay):
    c, p, ts = r['config'], r['pilot'], r['timeSeries']
    assert c['noiseDistribution'] == 'mixture' and c['nonGaussianOutlierProbability'] == .05
    assert c['nonGaussianOutlierScale'] == 8 and c['bearingObservationVariance'] == .0001
    assert c['bearingProcessNoise'] == .1 and c['immLowQ'] == .0001 and c['immHighQ'] == .1
    assert np.array_equal(c['immTransition'], [[.98, .02], [.10, .90]])
    sigma = float(c['bearingMeasurementStd']); assert sigma == (0. if method == 'oracle' else .01)
    x, y = r['trajectory']['x'][:, :-1].T, r['trajectory']['y'][:, :-1].T
    edges = np.asarray(p['linearEdges']).astype(int)-1; rec, nei = edges % 20, edges // 20
    dx, dy = x[:, nei]-x[:, rec], y[:, nei]-y[:, rec]
    dx -= 1000*np.round(dx/1000); dy -= 1000*np.round(dy/1000)
    truth = np.arctan2(dy, dx) % (2*np.pi)
    ids = np.asarray(p['selectedHistory']).astype(int).transpose(2, 0, 1)-1
    selected = np.zeros((8000, 20, 20), dtype=bool)
    selected[np.arange(8000)[:, None, None], np.arange(20)[None, :, None], ids] = True
    fresh = selected[:, rec, nei]; obs = p['observedBearingHistory']; used = p['usedBearingEvaluation']
    assert obs.shape == (8000, 380) and used.shape == (2000, 380)
    assert np.array_equal(np.isfinite(obs), fresh) and np.array_equal(np.isfinite(used), fresh[6000:])
    assert np.all(fresh.sum(1) == 300) and p['evalNoise']['count'] == 600000
    eo, eu = wrap(obs[6000:]-truth[6000:]), wrap(used-truth[6000:]); mask = fresh[6000:]
    obs_step = np.sqrt(np.nansum(eo**2, axis=1)/300)
    used_step = np.sqrt(np.nansum(eu**2, axis=1)/300)
    errors = {'observedAngleStepRMSE': float(np.max(np.abs(obs_step-ts['observedBearingRMSE'][6000:]))),
              'usedAngleStepRMSE': float(np.max(np.abs(used_step-ts['usedBearingRMSE'][6000:])))}
    if method in ['oracle', 'raw']: assert np.array_equal(obs[6000:], used, equal_nan=True)
    if method == 'oracle':
        assert np.max(np.abs(eu[mask])) < 1e-12 and np.size(p['estimator']) == 0
        assert np.all(ts['inputTruthRMS'] == 0)
    alpha = np.asarray(r['finalPreferredDirections']); kernel = float(c['receptiveFieldWidth'])
    h = float(c['totalSocialAttraction'])/20; max_field = 0.; max_noise = 0.; max_sample_rms = 0.
    for s in p['samples']:
        step = int(s['step'])-1; ix = np.asarray(s['linearIndex']).astype(int)-1
        assert np.array_equal(ix, np.flatnonzero(selected[step].ravel(order='F')))
        re, ne = ix % 20, ix // 20
        assert np.max(np.abs(wrap(s['trueBearing']-truth[step, np.searchsorted(edges, ix)]))) < 1e-12
        independent_noise = s['standardNormal']*np.where(s['outlierMask'], 8., 1.)/np.sqrt(4.15)
        max_noise = max(max_noise, float(np.max(np.abs(independent_noise-s['standardNoise']))))
        assert np.max(np.abs(wrap(s['observedBearing']-s['trueBearing']-sigma*independent_noise))) < 1e-12
        order = np.asarray(s['orderedIds']).astype(int)-1; fields = []
        for angles in [s['usedBearing'], s['trueBearing']]:
            centers = np.full((20, 20), np.nan); centers[re, ne] = angles
            field = np.zeros((100, 20))
            for receiver in range(20):
                delta = np.abs(alpha[receiver, :, None]-centers[receiver, order[receiver]][None, :])
                delta = np.minimum(delta, 2*np.pi-delta)
                field[:, receiver] = h*np.exp(-.5*delta**2/kernel**2).sum(1)
            fields.append(field)
        max_field = max(max_field, float(np.max(np.abs(fields[0]-s['input']))), float(np.max(np.abs(fields[1]-s['truthInput']))))
        max_sample_rms = max(max_sample_rms, abs(float(np.sqrt(np.mean((fields[0]-fields[1])**2)))-float(ts['inputTruthRMS'][step])))
        if step >= 6000:
            cols = np.searchsorted(edges, ix)
            assert np.max(np.abs(wrap(s['usedBearing']-used[step-6000, cols]))) < 1e-12
    errors.update(archivedInputField=max_field, archivedInputRMS=max_sample_rms, generatedNoise=max_noise)
    assert max(errors.values()) < 1e-9, errors
    filter_check = None
    if replay and method in ['selected-kf', 'imm']:
        state = None; max_angle = 0.; max_prob = 0.
        for step in range(8000):
            if method == 'selected-kf':
                state, _, _ = matrix_reference.kf_matrix(state, obs[step], fresh[step], .3, .1, .0001)
                estimate = state['x'][:, 0]
            else:
                state, xx, _ = matrix_reference.imm_matrix(state, obs[step], fresh[step], .3, .1, .0001, np.array([[.98, .02], [.10, .90]]))
                estimate = xx[:, 0]
                if step >= 6000:
                    max_prob = max(max_prob, float(np.max(np.abs(state['mu'][:, 1]-p['highProbabilityEvaluation'][step-6000]))))
            if step >= 6000:
                f = fresh[step]; max_angle = max(max_angle, float(np.max(np.abs(wrap(estimate[f]-used[step-6000, f])))))
        assert max_angle < 1e-8 and max_prob < 1e-8
        filter_check = {'seed': int(c['seed']), 'method': method, 'steps': 8000,
                        'maxAngleError': max_angle, 'maxProbabilityError': max_prob}
    return {'observedBearingRMSE': float(np.sqrt(np.mean(eo[mask]**2))),
            'usedBearingRMSE': float(np.sqrt(np.mean(eu[mask]**2))),
            'inputTruthRMS': float(np.sqrt(np.mean(ts['inputTruthRMS'][6000:]**2)))}, errors, filter_check


def classify(low, high, margin=.05):
    if low > margin: return '支持超过0.05的正差异'
    if high < -margin: return '支持超过0.05的负差异'
    if low >= -margin and high <= margin: return '支持本轮±0.05实际等效'
    if low > 0: return '支持正差异，未确认超过0.05'
    if high < 0: return '支持负差异，未确认超过0.05'
    return '方向未定；按区间报告可排除的差异'


def statistics(rows, meta):
    seeds = list(map(int, meta['seeds'])); n = len(seeds); df = n-1
    key = {(r['seed'], r['method']): r for r in rows}
    values = np.array([[key[s, m]['meanGO'] for m in METHODS] for s in seeds])
    delta = np.column_stack([values[:, METHODS.index(a)]-values[:, METHODS.index(b)] for a, b in CONTRASTS])
    means = delta.mean(0); se = delta.std(0, ddof=1)/np.sqrt(n)
    t95 = float(student_t.ppf(.975, df)); simultaneous = float(student_t.ppf(1-.05/12, df))
    rng = np.random.default_rng(int(meta['bootstrapSeed'])); repeats = int(meta['bootstrapRepeats'])
    boot = np.empty((repeats, 6))
    for start in range(0, repeats, 2000):
        stop = min(start+2000, repeats); ix = rng.integers(0, n, size=(stop-start, n))
        boot[start:stop] = delta[ix].mean(1)
    bootstrap = np.quantile(boot, [.05/12, 1-.05/12], axis=0)
    comparisons = []
    for i, (a, b) in enumerate(CONTRASTS):
        low, high = float(means[i]-simultaneous*se[i]), float(means[i]+simultaneous*se[i])
        r = {'contrast': f'{a}-{b}', 'a': a, 'b': b, 'nSeeds': n,
             'difference': float(means[i]), 'pairedSD': float(delta[:, i].std(ddof=1)),
             'CI95Low': float(means[i]-t95*se[i]), 'CI95High': float(means[i]+t95*se[i]),
             'family95Low': low, 'family95High': high,
             'pUnadjusted': float(2*student_t.sf(abs(means[i]/se[i]), df)),
             'bootstrapFamily95Low': float(bootstrap[0, i]), 'bootstrapFamily95High': float(bootstrap[1, i]),
             'positiveSeeds': int(np.sum(delta[:, i] > 0)), 'negativeSeeds': int(np.sum(delta[:, i] < 0)),
             'canExcludeGainAtLeast005': high < .05, 'canExcludeLossAtLeast005': low > -.05,
             'decision': classify(low, high)}
        comparisons.append(r)
    previous = 0.
    for rank, index in enumerate(np.argsort([r['pUnadjusted'] for r in comparisons])):
        adjusted = min(1., max(previous, (6-rank)*comparisons[index]['pUnadjusted']))
        comparisons[index]['pHolmSix'] = adjusted; previous = adjusted
    summaries = []
    for m in METHODS:
        group = [key[s, m] for s in seeds]; result = {'method': m, 'nSeeds': n}
        for field in ['meanGO', 'meanLO', 'meanPairDistanceNorm', 'meanSpeed', 'observedBearingRMSE', 'usedBearingRMSE', 'inputTruthRMS']:
            v = np.array([r[field] for r in group]); av = float(v.mean()); sem = float(v.std(ddof=1)/np.sqrt(n))
            result.update({field: av, field+'CI95Low': av-t95*sem, field+'CI95High': av+t95*sem})
        summaries.append(result)
    return comparisons, summaries, delta


def plot_summary(rows, contrasts, summaries, series, seeds, out):
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.axvspan(-.05, .05, color='#d8e7dc', alpha=.7, label='本轮操作等效带 ±0.05')
    ax.axvline(0, color='grey', lw=1)
    for i, r in enumerate(contrasts):
        ax.errorbar(r['difference'], i, xerr=[[r['difference']-r['family95Low']], [r['family95High']-r['difference']]], fmt='o', color='#245a81', capsize=5)
        ax.plot([r['bootstrapFamily95Low'], r['bootstrapFamily95High']], [i+.12, i+.12], color='#b16d22', lw=2)
    ax.set_yticks(range(6), [f'{NAMES[r["a"]]} − {NAMES[r["b"]]}' for r in contrasts]); ax.invert_yaxis()
    ax.set_xlabel('评价窗平均 GO 的配对差值；正值表示前者更高')
    ax.set_title('64个新种子：六对比较的同时95%区间\n蓝：Bonferroni t区间；橙：配对种子bootstrap敏感性区间')
    ax.grid(axis='x', alpha=.2); ax.legend(loc='lower right'); fig.tight_layout(); fig.savefig(out/'six-contrast-intervals.png'); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 6)); rng = np.random.default_rng(20261021)
    values = [[r['meanGO'] for r in rows if r['method'] == m] for m in METHODS]
    parts = ax.violinplot(values, positions=range(4), showmeans=False, showextrema=False)
    for p, color in zip(parts['bodies'], COLORS): p.set_facecolor(color); p.set_alpha(.2)
    for i, (v, color) in enumerate(zip(values, COLORS)):
        ax.scatter(i+rng.uniform(-.1, .1, len(v)), v, s=18, color=color, alpha=.65)
        a = summaries[i]; ax.errorbar(i, a['meanGO'], yerr=[[a['meanGO']-a['meanGOCI95Low']], [a['meanGOCI95High']-a['meanGO']]], fmt='D', color='black', capsize=5)
    ax.set_xticks(range(4), [NAMES[m] for m in METHODS]); ax.set_ylabel('每种子评价窗平均 GO'); ax.set_ylim(0, 1.02)
    ax.set_title('全部64个新种子；黑点与误差线为方法均值及普通95%区间\n方法间推断使用配对差值区间，不由这些误差线是否重叠判断')
    ax.grid(axis='y', alpha=.15); fig.tight_layout(); fig.savefig(out/'all-seed-go-distributions.png'); plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharex=True, sharey=True)
    key = {(r['seed'], r['method']): r for r in rows}
    for ax, (a, b) in zip(axes.ravel(), CONTRASTS):
        x = [key[s, b]['meanGO'] for s in seeds]; y = [key[s, a]['meanGO'] for s in seeds]
        ax.plot([0, 1], [0, 1], '--', color='grey'); ax.scatter(x, y, s=22, alpha=.7)
        ax.set(xlabel=NAMES[b]+' GO', ylabel=NAMES[a]+' GO', xlim=(0, 1), ylim=(0, 1)); ax.grid(alpha=.15)
    fig.suptitle('64个新种子的配对散点；对角线上方表示纵轴方法更高'); fig.tight_layout(rect=(0, 0, 1, .96)); fig.savefig(out/'paired-go-scatter.png'); plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 4.8)); time = np.arange(1, 8001)*.3
    for m, color in zip(METHODS, COLORS): ax.plot(time, np.mean([series[s, m] for s in seeds], axis=0), color=color, label=NAMES[m])
    ax.axvspan(1800, 2400, color='grey', alpha=.12); ax.set(xlabel='模型时间；灰区为评价窗', ylabel='逐时刻GO均值', ylim=(0, 1))
    ax.set_title('64个新种子的逐时刻均值，未平滑；不是持续稳定或成功比例'); ax.legend(ncol=4); ax.grid(alpha=.15)
    fig.tight_layout(); fig.savefig(out/'mean-go-timeseries.png'); plt.close(fig)


def report(root, out, meta, summaries, contrasts, audit):
    lines = ['# 低强度非高斯四组对照：64个新种子独立验证', '',
             'Material Passport：academic-research-suite / experiment-agent；2026-10-21。', '',
             '## 你的问题或任务', '', '在新种子中检验真实方位、带噪直接输入、KF与IMM的群体效果，判断前期差距能否重复。', '',
             '## 我的判断', '', '判断以以下六对同时区间为准；不把未显著当等效，不把估计误差下降当成群体收益。', '']
    for r in contrasts:
        lines.append(f'- {NAMES[r["a"]]}−{NAMES[r["b"]]}：{r["difference"]:+.4f}，同时95%区间[{r["family95Low"]:+.4f}, {r["family95High"]:+.4f}]；{r["decision"]}。')
    lines += ['', '## 我做了什么', '',
              '冻结64个新种子73001–73064，四臂256条轨迹；旧6201–6206不进入本轮统计。N20/M100、逐步随机k15、等权h_t=.24、核宽.4、beta1000、dt=.3、L1000、8000步、评价6001:8000；没有新增遮挡或领航。', '',
              '噪声为同方差95%小Gaussian/5%大Gaussian混合，标准差比8，总体σ=.01 rad。单KF q=.1；IMM低/高q=1e-4/.1、转移[[.98,.02],[.10,.90]]、R=.0001；不重新调参。真实方位仅进入oracle对照，实际估计器只用各自带噪观测。四组共享初态、全部选邻、标准扰动和标签，闭环后实际角度可不同。', '',
              '## 结果', '', '| 方法 | GO均值 | GO普通95%区间 | LO均值 | 平均对距/L | 方位RMSE rad |', '|---|---:|---|---:|---:|---:|']
    for r in summaries:
        lines.append(f'| {NAMES[r["method"]]} | {r["meanGO"]:.6f} | [{r["meanGOCI95Low"]:.6f}, {r["meanGOCI95High"]:.6f}] | {r["meanLO"]:.6f} | {r["meanPairDistanceNorm"]:.6f} | {r["usedBearingRMSE"]:.6f} |')
    lines += ['', '| 配对比较 | GO差值 | 六项同时95%区间 | Holm p | 正差值种子数 |', '|---|---:|---|---:|---:|']
    for r in contrasts:
        lines.append(f'| {NAMES[r["a"]]}−{NAMES[r["b"]]} | {r["difference"]:+.6f} | [{r["family95Low"]:+.6f}, {r["family95High"]:+.6f}] | {r["pHolmSix"]:.6g} | {r["positiveSeeds"]}/64 |')
    lines += ['', '同时区间按Bonferroni六项构造：每项覆盖99.1667%，全家族至少95%；与普通95%区间区别见contrasts.csv。±.05为本轮预定实际差异操作标准，不是物理或生物阈值。bootstrap同时百分位区间见contrasts.csv，作为敏感性分析，不替换预定主要检验。', '',
              '![六对差值与同时区间](six-contrast-intervals.png)', '',
              '![全部种子的GO](all-seed-go-distributions.png)', '',
              '![预先指定73001四组轨迹](trajectories/mixture-std-010-seed-73001-evaluation.png)', '',
              '## 原因分析', '',
              '本轮只能识别固定条件下的总效果，不能单独区分神经环自身抗扰、滤波滞后和闭环反馈机制。oracle是零观测误差参考，不是GO上界；无噪声输入也不保证更高有序度。平均GO、凝聚和持续稳定是不同问题。', '',
              '如果前6种子的差距未在新种子重复，应以本轮新种子推断为主，旧结果保留为探索性发现；不能用合并旧种子或调整参数制造显著性。条件由前期探索筛出，所以证据属于条件内独立验证，不能泛化到其他噪声、遮挡或参数。', '',
              '## 后续建议', '',
              '先按六对同时区间决定哪些差异得到支持、哪些大小的收益已可排除；不因结果不显著自动继续加种子或换滤波器。若群体收益与误差改善分离，后续再做固定状态下的时间结构或神经环机制干预。', '',
              '## 验证与交付', '',
              f'全部256条位置轨迹独立重建GO/LO/对距/速度及全部评价窗方位误差，最大误差为：{json.dumps(audit["maxErrors"], ensure_ascii=False)}。独立样本始终为64个种子，153600000个边×步误差不构成153600000次独立重复。', '',
              '全部归档输入场样本独立重建，未声称完整2000步输入场逐步第二实现重建；预定73001和73064的KF/IMM各8000步独立矩阵重放，共4条，未声称所有滤波器内部状态均被第二实现重放。全部初态/选邻/基础噪声/标签配对、文件哈希、计数、零噪oracle和旧轨迹回归通过。', '',
              '源码/参数/方案均在新种子开始前冻结。两个单线程MATLAB进程写不相交目录。原始模型、旧仿真器/估计器、历史结果、论文和已有tracked改动保持不变。分支exp/non-gaussian-four-group-validation-v1，HEAD f8b73cf；无commit/push/merge。建议保留实验分支，审查后再决定集成；回退只撤销新增代码，结果继续保留。']
    (out/'64新种子独立验证报告.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run', type=Path, required=True); args = ap.parse_args()
    root = args.run.resolve(); meta = json.loads((root/'metadata.json').read_text(encoding='utf-8'))
    seeds = list(map(int, meta['seeds'])); assert seeds == list(range(73001, 73065))
    for worker in [1, 2]:
        done = json.loads((root/f'worker-{worker}/completion.json').read_text(encoding='utf-8'))
        assert done['status'] == 'COMPLETED' and done['runs'] == 128
    for key, digest in meta['analysisScriptHashes'].items():
        name = key[:-3]+'.py' if key.endswith('_py') else key
        assert reference.sha(root/'analysis-scripts'/name) == digest
    out = root/'independent-analysis'; assert not out.exists(), 'Existing analysis protected'; out.mkdir()
    tracks = out/'trajectories'; tracks.mkdir(); rows = []; checks = []; replays = []; hashes = {}; series = {}; max_errors = {}
    for si, seed in enumerate(seeds):
        worker = 1 if seed % 2 else 2; folder = root/f'worker-{worker}'
        marker = json.loads((folder/f'completed-seed-{seed}.json').read_text(encoding='utf-8')); assert marker['runs'] == 4
        published = list(csv.DictReader((folder/f'seed-{seed}-summary.csv').open(encoding='utf-8-sig', newline='')))
        assert len(published) == 4 and {r['method'] for r in published} == set(METHODS)
        published = {r['method']: r for r in published}; quartet = {}; first = None
        for method in METHODS:
            file = folder/f'seed-{seed}-{method}.mat'; digest = reference.sha(file)
            key = re.sub('[^A-Za-z0-9_]', '_', file.name)
            assert marker['fileHashes'][key] == digest; hashes[str(file.relative_to(root))] = digest
            r = loadmat(file, simplify_cells=True)['result']; quartet[method] = r
            p = proof(r)
            if first is None: first = p
            else: equal_proof(first, p)
            row, error, go = reference.metrics(r); row['method'] = method
            angles, input_error, replay = check_input(r, method, seed in meta['filterReplaySeeds'])
            row.update(angles); error.update(input_error)
            for name, value in error.items(): max_errors[name] = max(max_errors.get(name, 0.), value)
            old = published[method]
            for name in ['meanGO', 'meanLO', 'meanPairDistanceNorm', 'meanSpeed', 'observedBearingRMSE', 'usedBearingRMSE', 'inputTruthRMS']:
                assert abs(row[name]-float(old[name])) < 1e-9, (seed, method, name)
            rows.append(row); series[seed, method] = go
            checks.append({'seed': seed, 'method': method, 'maxErrors': error, 'evaluationEdgeSteps': 600000})
            if replay: replays.append(replay)
        reference.plot_tracks(quartet, 'mixture', .01, seed, tracks)
        if seed in meta['exampleSeeds']: reference.plot_tracks(quartet, 'mixture', .01, seed, tracks, full=True)
        print(f'INDEPENDENT CHECK {si+1}/64: four runs, pairing, errors and trajectory figure PASS', flush=True)
        del quartet, first, r, p; gc.collect()
    assert len(rows) == len({(r['seed'], r['method']) for r in rows}) == 256
    assert len(replays) == 4 and len(hashes) == 256
    reference.write_csv(out/'independently-recomputed.csv', rows)
    contrasts, summaries, delta = statistics(rows, meta)
    reference.write_csv(out/'contrasts.csv', contrasts); reference.write_csv(out/'method-means.csv', summaries)
    write_json(out/'statistics.json', {'comparisons': contrasts, 'methodMeans': summaries,
                                      'simultaneousFamily': 'six contrasts, Bonferroni95% family t intervals',
                                      'practicalMargin': .05, 'bootstrapRepeats': 50000})
    reference.write_csv(out/'per-seed-differences.csv', [{'seed': s, **{f'{a}-{b}': float(delta[i, j]) for j, (a, b) in enumerate(CONTRASTS)}} for i, s in enumerate(seeds)])
    plot_summary(rows, contrasts, summaries, series, seeds, out)
    audit = {'status': 'PASS', 'runs': 256, 'seeds': seeds, 'pairedQuartets': 64,
             'independentUnit': 'seed', 'noExcludedRuns': True, 'maxErrors': max_errors,
             'allPairingChecksPassed': True, 'inputFieldsScope': 'all81 archived samples per run, not all2000 evaluation steps',
             'filterReplay': replays, 'fileSHA256': hashes, 'checks': checks,
             'evidence': meta['evidence'], 'figureCount': len(list(out.rglob('*.png')))}
    write_json(out/'audit.json', audit); report(root, out, meta, summaries, contrasts, audit)
    reference.write_csv(root/'summary.csv', rows)
    write_json(root/'completion.json', {'status': 'COMPLETED_AND_INDEPENDENTLY_VERIFIED', 'runs': 256, 'seeds': 64})
    (root/'README.md').write_text('# 64新种子四组独立验证\n\n[结果报告](independent-analysis/64新种子独立验证报告.md)\n\n冻结方案见pre-run-plan.md，参数和源码哈希见metadata.json；worker-1/2保存全部256条独立MAT轨迹与日志。independent-analysis保存全种子指标、六对同时区间、bootstrap敏感性、64组轨迹图及审计。\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'means': summaries, 'comparisons': contrasts}, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__': main()
