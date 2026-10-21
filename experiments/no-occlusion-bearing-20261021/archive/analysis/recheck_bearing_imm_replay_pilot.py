"""Independent matrix-form KF/IMM, geometry, errors and seed-level statistics."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, r'C:\Users\Administrator\AppData\Local\Temp\noisy-bearing-kf-pydeps')
import numpy as np
from scipy.io import loadmat
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def wrap(x):
    return np.arctan2(np.sin(x), np.cos(x))


def initial_state(n):
    return {'x': np.zeros((n, 2)), 'P': np.zeros((n, 2, 2)),
            'has': np.zeros(n, dtype=bool)}


def kf_matrix(s, y, fresh, dt, q, R):
    """Full matrix prediction and Joseph update, independently from MATLAB scalars."""
    n = len(y)
    if s is None:
        s = initial_state(n)
    old = s['has'].copy()
    init = fresh & ~old
    F = np.array([[1., dt], [0., 1.]])
    Q = float(q) * np.array([[dt**3 / 3, dt**2 / 2], [dt**2 / 2, dt]])
    xp = s['x'] @ F.T
    xp[:, 0] %= 2 * np.pi
    Pp = F @ s['P'] @ F.T + Q
    s['x'][old] = xp[old]
    s['P'][old] = Pp[old]
    S = Pp[:, 0, 0] + R
    innovation = np.full(n, np.nan)
    upd = fresh & old
    v = wrap(y[upd] - xp[upd, 0])
    K = Pp[upd, :, 0] / S[upd, None]
    s['x'][upd] = xp[upd] + K * v[:, None]
    s['x'][upd, 0] %= 2 * np.pi
    A = np.repeat(np.eye(2)[None, :, :], int(upd.sum()), axis=0)
    A[:, :, 0] -= K
    s['P'][upd] = A @ Pp[upd] @ A.transpose(0, 2, 1) + R * K[:, :, None] * K[:, None, :]
    innovation[upd] = v
    s['x'][init] = np.column_stack([y[init] % (2 * np.pi), np.zeros(int(init.sum()))])
    s['P'][init] = np.diag([R, 1.])
    s['has'][init] = True
    return s, innovation, S


def circular_mix(models, weights):
    angles = np.column_stack([m['x'][:, 0] for m in models])
    z = np.sum(weights * np.exp(1j * angles), axis=1)
    theta = np.angle(z) % (2 * np.pi)
    weak = np.abs(z) < 1e-10
    best = weights.argmax(axis=1)
    theta[weak] = angles[np.arange(len(theta))[weak], best[weak]]
    omega = sum(weights[:, i] * models[i]['x'][:, 1] for i in range(2))
    x = np.column_stack([theta, omega])
    P = np.zeros((len(theta), 2, 2))
    for i in range(2):
        delta = models[i]['x'] - x
        delta[:, 0] = wrap(delta[:, 0])
        P += weights[:, i, None, None] * (models[i]['P'] + delta[:, :, None] * delta[:, None, :])
    return x, P


def imm_matrix(s, y, fresh, dt, highq, R, transition):
    n = len(y)
    if s is None:
        s = {'models': [initial_state(n), initial_state(n)], 'mu': np.full((n, 2), .5)}
    old = s['models'][0]['has'].copy()
    init = fresh & ~old
    prior = s['mu'] @ transition
    mixed = []
    for j in range(2):
        w = s['mu'] * transition[:, j][None, :] / prior[:, j, None]
        x, P = circular_mix(s['models'], w)
        m = {k: v.copy() for k, v in s['models'][j].items()}
        m['x'][old], m['P'][old] = x[old], P[old]
        mixed.append(m)
    innovations, variances = [], []
    for j, q in enumerate((1e-4, highq)):
        m, v, S = kf_matrix(mixed[j], y, fresh, dt, q, R)
        s['models'][j] = m
        innovations.append(v)
        variances.append(S)
    s['mu'][old] = prior[old]
    upd = fresh & old
    innovation = np.column_stack(innovations)[upd]
    S = np.column_stack(variances)[upd]
    logw = np.log(prior[upd]) - .5 * (np.log(2 * np.pi * S) + innovation**2 / S)
    if upd.any():
        w = np.exp(logw - logw.max(axis=1, keepdims=True))
        s['mu'][upd] = w / w.sum(axis=1, keepdims=True)
    s['mu'][init] = .5
    x, P = circular_mix(s['models'], s['mu'])
    return s, x, P


def independent_replay(seq, output, frozen):
    bearing = seq['bearing']
    fresh = seq['fresh'].astype(bool)
    cfg = seq['config']
    B, T = int(cfg['burnIn']), int(cfg['NSteps'])
    methods = output['methods']
    states = [None] * len(methods)
    max_angle, max_prob, max_variance = 0., 0., 0.
    for t in range(T):
        y, f = bearing[t], fresh[t]
        for m, method in enumerate(methods):
            if method['type'] == 'raw':
                estimate = y
            elif method['type'] == 'kf':
                states[m], _, _ = kf_matrix(states[m], y, f, float(cfg['dt']), float(method['q']), float(cfg['bearingObservationVariance']))
                estimate = states[m]['x'][:, 0]
            else:
                states[m], x, P = imm_matrix(states[m], y, f, float(cfg['dt']), float(method['q']), float(cfg['bearingObservationVariance']), np.asarray(frozen['transitionPerModelStep']))
                estimate = x[:, 0]
                if t >= B:
                    max_prob = max(max_prob, float(np.max(np.abs(states[m]['mu'][:, 1] - output['immHighProbability'][t-B]))))
                    max_variance = max(max_variance, float(np.max(np.abs(P[:, 0, 0] - output['immVariance'][t-B]))))
            if t >= B:
                # Uninitialized unselected edges do not affect observations; all
                # initialized predictions are compared, not just sampled edges.
                valid = np.isfinite(estimate)
                error = wrap(estimate[valid] - output['estimate'][t-B, valid, m])
                max_angle = max(max_angle, float(np.max(np.abs(error))) if len(error) else 0.)
    assert max_angle < 1e-8 and max_prob < 1e-8 and max_variance < 1e-7, (max_angle, max_prob, max_variance)
    return {'maxAngleDifference': max_angle, 'maxProbabilityDifference': max_prob,
            'maxVarianceDifference': max_variance, 'steps': T}


def check_geometry(seq, reference):
    cfg = seq['config']
    n, T = int(cfg['NAgents']), int(cfg['NSteps'])
    x = reference['trajectory']['x'][:, :T].T
    y = reference['trajectory']['y'][:, :T].T
    dx, dy = x[:, None, :] - x[:, :, None], y[:, None, :] - y[:, :, None]
    L = float(cfg['arenaSize'])
    dx -= L * np.round(dx / L)
    dy -= L * np.round(dy / L)
    edges = np.asarray(seq['linearEdges']).astype(int) - 1
    rec, nei = edges % n, edges // n
    truth = np.arctan2(dy[:, rec, nei], dx[:, rec, nei]) % (2 * np.pi)
    maximum = float(np.max(np.abs(wrap(truth - seq['truth']))))
    assert maximum < 1e-12
    ids = np.asarray(reference['pilot']['selectedHistory']).astype(int).transpose(2, 0, 1) - 1
    selected = np.zeros((T, n, n), dtype=bool)
    selected[np.arange(T)[:, None, None], np.arange(n)[None, :, None], ids] = True
    fresh = selected[:, rec, nei]
    assert np.array_equal(fresh, seq['fresh'].astype(bool))
    delta = np.full_like(truth, np.nan)
    delta[1:] = np.abs(wrap(truth[1:] - truth[:-1]))
    assert np.nanmax(np.abs(delta - seq['trueDelta'])) < 1e-12
    switch = np.zeros_like(fresh)
    switch[1:] = (np.abs(np.diff(dx[:, rec, nei], axis=0)) > L / 2) | (np.abs(np.diff(dy[:, rec, nei], axis=0)) > L / 2)
    assert np.array_equal(switch, seq['periodicSwitch'].astype(bool))
    return maximum


def recompute_metrics(seq, output):
    cfg = seq['config']
    B = int(cfg['burnIn'])
    truth = seq['truth'][B:]
    mask = seq['fresh'][B:].astype(bool)
    delta = seq['trueDelta'][B:]
    quiet, jump = mask & (delta < .02), mask & (delta > .5)
    recent = np.zeros_like(mask)
    for lag in range(11):
        recent |= seq['trueDelta'][B-lag:len(seq['truth'])-lag] > .5
    recent &= mask
    edge, when = np.nonzero(jump.T)
    events = output['events']
    assert np.array_equal(np.atleast_1d(events['relativeStep']).astype(int), when + 1)
    assert np.array_equal(np.atleast_1d(events['edgeColumn']).astype(int), edge + 1)
    saved_duration = np.asarray(events['durationSteps']).reshape(len(when), len(output['methods']))
    saved_recovered = np.asarray(events['recovered']).reshape(len(when), len(output['methods'])).astype(bool)
    rows = []
    for m, method in enumerate(output['methods']):
        error = wrap(output['estimate'][:, :, m] - truth)
        a = np.abs(error[mask])
        ordered = np.sort(a)
        count = len(a)
        sq = a*a
        take = int(np.ceil(.01 * count))
        top = np.partition(sq, count-take)[-take:].sum() / sq.sum()
        good = mask & (np.abs(error) <= 3 * float(cfg['bearingMeasurementStd']))
        duration = np.minimum(30, len(mask)-1-when)
        recovered = np.zeros(len(when), dtype=bool)
        for lag in range(2, 31):
            valid = (~recovered) & (when + lag < len(mask))
            ids = np.flatnonzero(valid)
            ok = good[when[ids]+lag, edge[ids]] & good[when[ids]+lag-1, edge[ids]] & good[when[ids]+lag-2, edge[ids]]
            ids = ids[ok]
            duration[ids], recovered[ids] = lag, True
        assert np.array_equal(duration, saved_duration[:, m])
        assert np.array_equal(recovered, saved_recovered[:, m])
        row = {'seed': int(cfg['seed']), 'distribution': seq['distribution'],
               'measurementStd': float(cfg['bearingMeasurementStd']), 'method': method['name'],
               'q': float(method['q']), 'count': count, 'RMSE': float(np.sqrt(sq.mean())),
               'medianAbs': float(np.median(a)), 'P95Abs': float(ordered[int(np.ceil(.95*count))-1]),
               'P99Abs': float(ordered[int(np.ceil(.99*count))-1]), 'topOnePercentMSEShare': float(top),
               'quietCount': int(quiet.sum()), 'quietRMSE': float(np.sqrt(np.mean(error[quiet]**2))),
               'jumpCount': int(jump.sum()), 'jumpRMSE': float(np.sqrt(np.mean(error[jump]**2))),
               'afterJumpCount': int(recent.sum()), 'afterJumpRMSE': float(np.sqrt(np.mean(error[recent]**2))),
               'recoveryEventCount': len(when), 'recoveredEventCount': int(recovered.sum()),
               'censoredEventCount': int((~recovered).sum())}
        # Descriptive KM at the predeclared horizon, retaining endpoint censoring.
        survival = 1.
        for time in np.unique(duration):
            deaths = int(np.sum((duration == time) & recovered))
            at_risk = int(np.sum(duration >= time))
            survival *= 1 - deaths / at_risk
        row['KMRecoveredBy30Steps'] = 1-survival if len(when) else np.nan
        row['medianRecoveredDelayModelTime'] = float(np.median(duration[recovered])*cfg['dt']) if recovered.any() else np.nan
        rows.append(row)
    high = output['immHighProbability']
    diagnostics = {'meanHighProbabilityQuiet': float(np.mean(high[quiet])),
                   'meanHighProbabilityJump': float(np.mean(high[jump])) if jump.any() else None,
                   'meanHighProbabilityOutlier': float(np.mean(high[mask & seq['outlier'][B:].astype(bool)])) if seq['distribution']=='mixture' else None,
                   'weakCircularMeans': int(output['weakCircularMeans']),
                   'endpointCensoredEvents': int(np.sum((saved_duration[:, -1] < 30) & ~saved_recovered[:, -1]))}
    return rows, diagnostics


def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def paired_comparison(rows, distribution, sd, base):
    selected = [r for r in rows if r['distribution']==distribution and r['measurementStd']==sd]
    key = {(r['seed'],r['method']): r for r in selected}
    seeds = sorted({r['seed'] for r in selected})
    diff = np.array([key[s,'IMM']['RMSE']-key[s,base]['RMSE'] for s in seeds])
    mean, se = diff.mean(), diff.std(ddof=1)/np.sqrt(len(diff))
    interval = student_t.ppf(.975, len(diff)-1)*se
    p = float(2*student_t.sf(abs(mean/se), len(diff)-1)) if se else (0. if mean else 1.)
    return {'distribution': distribution, 'measurementStd': sd, 'baseline': base,
            'nSeeds': len(seeds), 'meanDifference': float(mean), 'CI95Low': float(mean-interval),
            'CI95High': float(mean+interval), 'pUnadjusted': p,
            'IMMImprovedSeeds': int(np.sum(diff<0)), 'differences': diff.tolist()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    out = root/'independent-analysis'
    assert not out.exists(), 'Existing independent analysis is protected'
    out.mkdir()
    frozen = json.loads((root/'frozen-selection.json').read_text(encoding='utf-8'))
    assert frozen['status']=='FROZEN_BEFORE_VALIDATION'
    dev = read_csv(root/'development'/'summary.csv')
    dev_by = {(int(float(r['seed'])),r['distribution'],float(r['measurementStd']),r['method']):r for r in dev}
    scores = {}
    for typ, candidates, prefix in [('imm', frozen['immHighQGrid'],'IMM-q-'),('single',frozen['singleQGrid'],'single-q-')]:
        values = []
        for q in candidates:
            name = 'original-KF' if typ=='single' and q==1e-4 else f'{prefix}{q:g}'
            logs=[]
            for seed in range(5101,5107):
                for dist in ['gaussian','mixture']:
                    for sd in [.01,.1]:
                        a=dev_by[seed,dist,sd,name]
                        b=dev_by[seed,dist,sd,'raw']
                        logs.append(np.log(float(a['RMSE'])/float(b['RMSE'])))
            values.append(float(np.mean(logs)))
        scores[typ]=values
    assert np.max(np.abs(np.asarray(scores['imm'])-frozen['immScores']))<1e-12
    assert np.max(np.abs(np.asarray(scores['single'])-frozen['singleScores']))<1e-12
    assert frozen['immHighQ']==frozen['immHighQGrid'][int(np.argmin(scores['imm']))]
    assert frozen['bestSingleQ']==frozen['singleQGrid'][int(np.argmin(scores['single']))]
    archived=read_csv(root/'validation'/'summary.csv')
    by_key={(int(float(r['seed'])),r['distribution'],float(r['measurementStd']),r['method']):r for r in archived}
    assert len(archived)==len(by_key)==120
    recomputed, audits, hashes=[],[],{}
    max_metric=0.
    files=sorted((root/'validation').glob('sequence-*.mat'))
    assert len(files)==24
    for i,path in enumerate(files,1):
        tag=path.stem.removeprefix('sequence-')
        seq=loadmat(path,simplify_cells=True)['seq']
        output_path=path.with_name('estimates-'+tag+'.mat')
        output=loadmat(output_path,simplify_cells=True)['o']
        reference_path=path.with_name('reference-'+tag+'.mat')
        reference=loadmat(reference_path,simplify_cells=True)['r']
        geometry=check_geometry(seq,reference)
        full_replay=independent_replay(seq,output,frozen)
        rows,diagnostics=recompute_metrics(seq,output)
        for row in rows:
            a=by_key[row['seed'],row['distribution'],row['measurementStd'],row['method']]
            for name in a:
                if name in ['distribution','method']:
                    assert a[name]==row[name]
                else:
                    discrepancy=abs(float(a[name])-row[name])
                    max_metric=max(max_metric,discrepancy)
                    assert discrepancy<1e-10, (tag,name,discrepancy)
        recomputed.extend(rows)
        audits.append({'tag':tag,'geometryMaxDifference':geometry,**full_replay,**diagnostics})
        for p in [path,output_path,reference_path]:
            hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        print(f'RECHECK {i}/24 {tag} angle={full_replay["maxAngleDifference"]:.3g} probability={full_replay["maxProbabilityDifference"]:.3g}',flush=True)
        del seq,output,reference
    assert {r['seed'] for r in recomputed}==set(range(6101,6107))
    assert all(r['count']==600000 for r in recomputed)
    write_csv(out/'independently-recomputed.csv',recomputed)
    primary=[paired_comparison(recomputed,d,sd,'selected-single-KF') for d in ['gaussian','mixture'] for sd in [.01,.1]]
    order=np.argsort([x['pUnadjusted'] for x in primary]); previous=0.
    for rank,index in enumerate(order):
        adjusted=min(1.,max(previous,(4-rank)*primary[index]['pUnadjusted']))
        primary[index]['pHolmFour']=adjusted; previous=adjusted
    secondary=[paired_comparison(recomputed,d,sd,b) for d in ['gaussian','mixture'] for sd in [.01,.1] for b in ['raw','original-KF','same-high-Q-KF']]
    (out/'paired-statistics.json').write_text(json.dumps({'primary':primary,'secondaryUnadjustedExploratory':secondary},indent=2),encoding='utf-8')
    write_csv(out/'primary-paired-statistics.csv',[{k:v for k,v in r.items() if k!='differences'} for r in primary])
    methods=['raw','original-KF','same-high-Q-KF','selected-single-KF','IMM']
    labels=['Raw','Original KF','Same high Q','Selected KF','IMM']
    colors=['#5c7386','#a75645','#bd9458','#658464','#4269a6']
    for metric,filename,title in [('RMSE','paired-rmse.png','Angular RMSE: identical observations, six new seeds'),('P99Abs','paired-p99.png','99th percentile absolute error: identical observations')]:
        fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
        for ax,(dist,sd) in zip(axes.flat,[(d,s) for d in ['gaussian','mixture'] for s in [.01,.1]]):
            group=[r for r in recomputed if r['distribution']==dist and r['measurementStd']==sd]
            key={(r['seed'],r['method']):r for r in group}
            values=np.array([[key[seed,m][metric] for m in methods] for seed in range(6101,6107)])
            for value in values: ax.plot(range(5),value,color='#b3bac3',alpha=.6,lw=.8,marker='.',ms=3)
            means=values.mean(axis=0)
            for x,mean,color in zip(range(5),means,colors):ax.plot(x,mean,'o',color=color,ms=7,zorder=4)
            ax.set_xticks(range(5),labels,rotation=18)
            ax.set_yscale('log');ax.set_ylabel(metric+' (rad)');ax.set_title(f'{dist}, sigma={sd:.2f}')
            ax.grid(axis='y',alpha=.2)
        fig.suptitle(title+'\nGray lines: paired seeds; colored points: seed means. No closed-loop IMM run.',fontsize=12)
        fig.savefig(out/filename,dpi=170);plt.close(fig)
    # Illustrative failure selected by a predeclared baseline, not IMM success.
    seq=loadmat(root/'validation'/'sequence-gaussian-seed-6101-std-010.mat',simplify_cells=True)['seq']
    o=loadmat(root/'validation'/'estimates-gaussian-seed-6101-std-010.mat',simplify_cells=True)['o']
    truth=seq['truth'][6000:];error=np.abs(wrap(o['estimate']-truth[:,:,None]));mask=seq['fresh'][6000:].astype(bool)
    target=np.where(mask,error[:,:,1],-1.)
    at,edge=np.unravel_index(np.argmax(target),target.shape)
    start,end=max(0,at-15),min(len(truth),at+36);steps=np.arange(start,end)+6001
    fig,axes=plt.subplots(3,1,figsize=(10,7),layout='constrained',sharex=True)
    axes[0].plot(steps,seq['trueDelta'][6000+start:6000+end,edge],color='#555555');axes[0].set_ylabel('True step change\n(rad)')
    for m in [0,1,3,4]:
        a=error[start:end,edge,m].copy();a[~mask[start:end,edge]]=np.nan
        axes[1].plot(steps,a,'o-',label=labels[m],color=colors[m],ms=3,lw=1)
    axes[1].set_ylabel('Absolute error\n(rad)');axes[1].legend(ncol=4,fontsize=9)
    axes[2].plot(steps,o['immHighProbability'][start:end,edge],color=colors[-1]);axes[2].set_ylabel('IMM high-Q\nprobability');axes[2].set_ylim(-.03,1.03)
    for ax in axes:ax.axvline(at+6001,color='gray',ls=':');ax.grid(alpha=.2)
    axes[2].set_xlabel('Model step')
    linear=int(np.asarray(seq['linearEdges'])[edge])-1
    fig.suptitle(f'Illustrative worst original-KF error: Gaussian sigma=.01, seed 6101\nreceiver {linear%20+1} -> neighbor {linear//20+1}; gaps are unselected edges',fontsize=11)
    fig.savefig(out/'baseline-worst-event-example.png',dpi=170);plt.close(fig)
    audit={'status':'PASS_INDEPENDENT_MATRIX_REPLAY_AND_STATISTICS','sequences':24,'summaryRows':120,
           'newSeeds':list(range(6101,6107)),'selectedPointsPerSequence':600000,
           'independentUnit':'seed; edges, events and time points are not independent repetitions',
           'fullFilterReplaySteps':24*8000,'maxMetricDifference':max_metric,
           'maxAngleDifference':max(x['maxAngleDifference'] for x in audits),
           'maxProbabilityDifference':max(x['maxProbabilityDifference'] for x in audits),
           'maxVarianceDifference':max(x['maxVarianceDifference'] for x in audits),
           'developmentSelectionRecomputed':True,'checks':audits,'fileHashes':hashes}
    (out/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    aggregate=[]
    for dist in ['gaussian','mixture']:
        for sd in [.01,.1]:
            for m in methods:
                group=[r for r in recomputed if r['distribution']==dist and r['measurementStd']==sd and r['method']==m]
                aggregate.append({'distribution':dist,'measurementStd':sd,'method':m,'nSeeds':len(group),
                                  **{k:float(np.mean([r[k] for r in group])) for k in ['RMSE','medianAbs','P95Abs','P99Abs','quietRMSE','jumpRMSE','afterJumpRMSE','KMRecoveredBy30Steps','medianRecoveredDelayModelTime']}})
    write_csv(out/'condition-means.csv',aggregate)
    (out/'condition-means.json').write_text(json.dumps(aggregate,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k not in ['checks','fileHashes']},indent=2),flush=True)
    print(json.dumps({'frozenIMMHighQ':frozen['immHighQ'],'frozenBestSingleQ':frozen['bestSingleQ'],'primary':primary,'conditionMeans':aggregate},indent=2),flush=True)


if __name__=='__main__':
    main()
