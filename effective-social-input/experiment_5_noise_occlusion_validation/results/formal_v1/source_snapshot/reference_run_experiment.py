"""Development-only informed-turn experiment; all output directories are exclusive.

The Python implementation is checked against the repository's MATLAB field,
bearingTrackStep and occlusionInputStep before scientific runs. Hidden truth is
confined to the sensor and offline metrics, never supplied to Memory.input().
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
import platform
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SEEDS = [2026101401, 2026101402, 2026101403]
MODES = ['fresh', 'drop', 'drop-normalized', 'hold', 'kf-mean', 'kf-uncertainty']


def wrap(x):
    return (x + np.pi) % (2*np.pi) - np.pi


def ring(m):
    angles = np.arange(m)*2*np.pi/m
    return np.cos(np.pi*np.sqrt(np.abs(wrap(angles[:, None]-angles[None, :]))/np.pi))


@dataclass
class Config:
    n: int = 20
    m: int = 100
    k: int = 15
    dt: float = .3
    arena: float = 1000.
    beta: float = 1000.
    speed_scale: float = .05
    h0: float = .012
    sigma: float = .4
    q: float = 1e-4
    R: float = 1e-4
    noise_std: float = .01
    max_age: float = 120.
    max_variance: float = 1.
    warm_steps: int = 100
    lead_steps: int = 20
    total_steps: int = 500
    cue_amplitude: float = .08
    cue_shape: str = 'gaussian'
    turn_rate: float = .025
    turn_angle: float = np.pi/4
    geometry_scale: float = 1.


class Memory:
    """Vectorized equivalent of the existing scalar CV-KF and input interface."""
    def __init__(self, selected, cfg):
        self.selected = selected.copy()
        self.cfg = cfg
        self.x = np.zeros((*selected.shape, 2))
        self.P = np.zeros((*selected.shape, 2, 2))
        self.initialized = np.zeros_like(selected)
        self.last = np.full(selected.shape, np.nan)
        self.age = np.full(selected.shape, np.inf)

    def input(self, observed, fresh, alpha, heading, mode):
        # Only the sanitized observation matrix enters this function.
        assert not np.any(fresh & ~self.selected)
        assert np.all(np.isfinite(observed[fresh]))
        assert np.all(np.isnan(observed[~fresh]))
        c = self.cfg
        F = np.array([[1., c.dt], [0., 1.]])
        Q = c.q*np.array([[c.dt**3/3, c.dt**2/2], [c.dt**2/2, c.dt]])
        old = self.initialized.copy()
        new = fresh & ~old
        self.x[old] = self.x[old] @ F.T
        self.P[old] = F @ self.P[old] @ F.T + Q
        self.age[old] += c.dt
        updated = old & fresh
        innovation = wrap(observed[updated]-self.x[updated, 0])
        prior = self.P[updated].copy()
        K = prior[:, :, 0]/(prior[:, 0, 0]+c.R)[:, None]
        self.x[updated] += K*innovation[:, None]
        A = np.broadcast_to(np.eye(2), prior.shape).copy()
        A[:, :, 0] -= K
        self.P[updated] = A @ prior @ A.transpose(0, 2, 1) + c.R*K[:, :, None]*K[:, None, :]
        self.x[new, 0] = observed[new]
        self.x[new, 1] = 0
        self.P[new] = np.diag([c.R, 1.])
        self.initialized |= new
        self.x[:, :, 0] %= 2*np.pi
        self.last[fresh] = observed[fresh]
        self.age[fresh] = 0
        active = self.selected & self.initialized
        missing = active & ~fresh
        if mode in ('drop', 'drop-normalized', 'fresh'):
            active &= fresh
        else:
            active &= self.age <= c.max_age
            if mode.startswith('kf'):
                active &= self.P[:, :, 0, 0] <= c.max_variance
        center = np.where(fresh, observed, self.last if mode == 'hold' else self.x[:, :, 0])
        center = np.where(active, center, 0.)
        variance = np.where(active & ~fresh, self.P[:, :, 0, 0], 0.) if mode == 'kf-uncertainty' else np.zeros_like(center)
        weights = active.astype(float)*c.h0
        if mode == 'drop-normalized':
            weights *= self.selected.sum(1)[:, None]/np.maximum(1, active.sum(1))[:, None]
        delta = wrap(alpha[:, None, :]-center[:, :, None])
        kernels = np.exp(-.5*(delta/c.sigma)**2)
        convolve = variance > 0
        freq = np.fft.fftfreq(alpha.shape[1])*alpha.shape[1]
        if np.any(convolve):
            kernels[convolve] = np.maximum(0., np.fft.ifft(np.fft.fft(kernels[convolve], axis=-1)*np.exp(-.5*variance[convolve, None]*freq**2), axis=-1).real)
        field = np.sum(weights[:, :, None]*kernels, axis=1).T
        factor = np.exp(np.cos(center-heading[:, None])-1.)
        if np.any(convolve):
            # E[exp(cos(theta-heading)-1)] under the same wrapped-normal
            # posterior. This is an offline scalar diagnostic, not a control.
            grid = np.arange(c.m)*2*np.pi/c.m
            coef = (np.fft.fft(np.exp(np.cos(grid)-1.))/c.m).real
            harmonics = np.arange(1, min(24, c.m//2))
            d = (center-heading[:, None])[convolve]
            factor[convolve] = coef[0]+2*np.sum(coef[harmonics]*np.exp(-.5*variance[convolve, None]*harmonics**2)*np.cos(d[:, None]*harmonics), axis=1)
        S = np.sum(weights*factor, axis=1)
        A_sum = weights.sum(1)
        details = dict(active=active, center=center, variance=variance, weights=weights,
                       S=S, A=A_sum, C=np.divide(S, A_sum, out=np.zeros_like(S), where=A_sum>0),
                       S_memory=np.sum(weights*factor*missing, axis=1),
                       prediction_variance=self.P[:, :, 0, 0].copy())
        return field, details


def topology(c):
    selected = np.zeros((c.n, c.n), bool)
    followers = np.arange(3, c.n)
    for i in range(3):
        ids = np.r_[np.setdiff1d(np.arange(3), i), followers[:c.k-2]]
        selected[i, ids] = True
    for f, i in enumerate(followers):
        selected[i, :3] = True
        selected[i, followers[(f+np.arange(1, c.k-2)) % len(followers)]] = True
    hidden = selected.copy()
    hidden[:3] = False
    hidden[:, 3:] = False
    assert np.all(selected.sum(1) == c.k) and not np.any(np.diag(selected))
    assert hidden.sum() == 51
    # All possible source-to-follower edges are masked; the follower ring
    # remains connected, preventing a fresh relay path around the mask.
    assert np.all(hidden[3:, :3])
    return selected, hidden


def initialize(seed, c):
    # Controlled ordered initial state, not a selected successful endpoint.
    # PCG64 is explicitly distinct from historical MATLAB mt19937ar seeds.
    rng = np.random.Generator(np.random.PCG64(seed))
    xy = np.empty((c.n, 2))
    xy[:3] = [[492, 492], [492, 500], [492, 508]]
    j = np.arange(c.n-3)
    xy[3:, 0] = 480-8*(j//4)
    xy[3:, 1] = 500+8*(j%4-1.5)
    xy = 500+c.geometry_scale*(xy-500)+rng.normal(0, .4, xy.shape)
    heading = rng.normal(0, .03, c.n)
    alpha = (np.arange(c.m)[None, :]*2*np.pi/c.m+heading[:, None]) % (2*np.pi)
    u = .1*np.cos(alpha.T-heading[None, :])+rng.normal(0, .001, (c.m, c.n))
    noise = rng.normal(0, c.noise_std, (c.total_steps, c.n, c.n))
    return xy, heading, alpha, u, noise


def geometry(xy, c):
    r = xy[None, :, :]-xy[:, None, :]
    r -= c.arena*np.round(r/c.arena)
    return np.arctan2(r[:, :, 1], r[:, :, 0]) % (2*np.pi), np.linalg.norm(r, axis=-1)


def simulate(seed, c, mode, duration, turning=True, initial=None, selected=None, hidden=None):
    if initial is None:
        initial = initialize(seed, c)
    xy, heading, alpha, u, noise = [v.copy() for v in initial]
    if selected is None:
        selected, hidden = topology(c)
    memory = Memory(selected, c)
    W = ring(c.m)
    start = c.warm_steps+c.lead_steps
    trajectory = np.empty((c.total_steps+1, c.n, 2))
    headings = np.empty((c.total_steps+1, c.n))
    velocity = np.empty((c.total_steps, c.n, 2))
    trajectory[0], headings[0] = xy, heading
    rows = []
    fixed_nearest = None
    checkpoint = {}
    for step in range(c.total_steps):
        theta, distance = geometry(xy, c)
        if step == start:
            checkpoint = {f'checkpoint_{key}':value.copy() for key,value in
                          dict(xy=xy,heading=heading,alpha=alpha,u=u,
                               filter_x=memory.x,filter_P=memory.P,
                               filter_initialized=memory.initialized,
                               last_bearing=memory.last,age=memory.age).items()}
            d0 = distance.copy()
            np.fill_diagonal(d0, np.inf)
            fixed_nearest = np.argsort(d0, axis=1)[:, :5]
        blackout = mode != 'fresh' and start <= step < start+duration
        fresh = selected & ~(hidden if blackout else np.zeros_like(hidden))
        observed = np.where(fresh, (theta+noise[step]) % (2*np.pi), np.nan)
        field, diagnostic = memory.input(observed, fresh, alpha, heading, mode)
        # External task cue is solely a leader input. Its schedule never
        # enters the follower controller or KF.
        elapsed = max(0, (step-c.warm_steps+1)*c.dt)
        command = min(c.turn_angle, c.turn_rate*elapsed) if turning else 0.
        if c.cue_shape == 'cosine':
            # Explicit exogenous leader-only steering: signed cosine excites
            # the desired half-ring and inhibits the opposite half-ring.
            # Follower social input and field dynamics are unchanged.
            cue = c.cue_amplitude*np.cos(alpha[:3]-command)
        else:
            cue = c.cue_amplitude*np.exp(-.5*(wrap(alpha[:3]-command)/c.sigma)**2)
        field[:, :3] += cue.T
        S_oracle = c.h0*np.sum(selected*np.exp(np.cos(theta-heading[:, None])-1.), axis=1)
        missing_S = c.h0*np.sum(hidden*np.exp(np.cos(theta-heading[:, None])-1.), axis=1)
        prediction_mask = hidden & blackout
        err = wrap(diagnostic['center']-theta)
        pred_active = prediction_mask & diagnostic['active']
        # All hidden edges also have a shadow KF; score on all edge-time
        # points to avoid a covariance-gating selection bias.
        shadow_error = wrap(memory.x[:, :, 0]-theta)
        hidden_err = shadow_error[prediction_mask]
        coverage = np.abs(hidden_err) <= 1.96*np.sqrt(memory.P[:, :, 0, 0][prediction_mask])
        u += c.dt*(-u + W@np.tanh(c.beta*u)/c.m + field)
        positive = np.maximum(np.tanh(c.beta*u), 0.)
        vx = c.speed_scale*np.sum(positive*np.cos(alpha.T), axis=0)
        vy = c.speed_scale*np.sum(positive*np.sin(alpha.T), axis=0)
        v = np.c_[vx, vy]
        speed = np.linalg.norm(v, axis=1)
        moving = speed > 1e-12
        direction = np.divide(v, speed[:, None], out=np.zeros_like(v), where=moving[:, None])
        heading[moving] = np.arctan2(vy[moving], vx[moving]) % (2*np.pi)
        xy = (xy+c.dt*v) % c.arena
        trajectory[step+1], headings[step+1], velocity[step] = xy, heading, v
        _, distance_after = geometry(xy, c)
        nearest = np.argsort(distance_after, axis=1)[:, :6]
        local = np.linalg.norm(direction[nearest].sum(axis=1), axis=1)/6
        d_no_self = distance_after.copy()
        np.fill_diagonal(d_no_self, np.inf)
        frozen_distance = float(np.mean(distance_after[np.arange(c.n)[:, None], fixed_nearest])) if fixed_nearest is not None else np.nan
        E = np.abs(wrap(heading-command))
        leader_dir = np.angle(np.mean(np.exp(1j*heading[:3])))
        row = dict(step=step, time=(step+1-start)*c.dt, command_deg=np.degrees(command),
                   blackout=int(blackout), GO=np.linalg.norm(direction.sum(0))/c.n,
                   LO=local.mean(), fixed5_distance=frozen_distance,
                   nearest_distance=np.min(d_no_self, axis=1).mean(),
                   pair_distance=distance_after[np.triu_indices(c.n, 1)].mean(),
                   leader_follower_distance=distance_after[3:, :3].mean(),
                   minimum_distance=np.min(d_no_self),
                   follower_turn_error_deg=np.degrees(E[3:]).mean(),
                   leader_turn_error_deg=np.degrees(E[:3]).mean(),
                   follower_leader_error_deg=np.degrees(np.abs(wrap(heading[3:]-leader_dir))).mean(),
                   follower_turn_angle_deg=np.degrees(np.angle(np.mean(np.exp(1j*heading[3:])))),
                   leader_turn_angle_deg=np.degrees(leader_dir),
                   mean_speed=speed.mean(), moving_fraction=moving.mean(),
                   S=diagnostic['S'][3:].mean(), A=diagnostic['A'][3:].mean(), C=diagnostic['C'][3:].mean(),
                   S_memory=diagnostic['S_memory'][3:].mean(), S_oracle=S_oracle[3:].mean(),
                   S_abs_error=np.abs(diagnostic['S'][3:]-S_oracle[3:]).mean(),
                   S_hidden_oracle=missing_S[3:].mean(),
                   fresh_edges=int(fresh.sum()), active_edges=int(diagnostic['active'].sum()),
                   hidden_active_edges=int(pred_active.sum()),
                   bearing_error_deg=np.degrees(np.abs(err[pred_active])).mean() if pred_active.any() else np.nan,
                   shadow_kf_rmse_rad=np.sqrt(np.mean(hidden_err**2)) if hidden_err.size else np.nan,
                   shadow_kf_coverage=coverage.mean() if coverage.size else np.nan,
                   hidden_bearing_change_deg=np.degrees(np.abs(wrap(theta-hidden_last))[hidden]).mean() if step >= start else np.nan)
        if step == start-1:
            hidden_last = theta.copy()
        rows.append(row)
    return dict(rows=rows, xy=trajectory, headings=headings, velocity=velocity, u=u, alpha=alpha,
                final_x=memory.x, final_P=memory.P, noise=noise, selected=selected, hidden=hidden,
                **checkpoint)


def write_csv(path, rows):
    with path.open('x', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def save_run(directory, seed, mode, T, turning, result):
    name = f'seed-{seed}_{"turn" if turning else "straight"}_{mode}_T-{T:g}'
    write_csv(directory/(name+'.csv'), result['rows'])
    np.savez_compressed(directory/(name+'.npz'), **{k:v for k,v in result.items() if k != 'rows'})
    return name


def assert_prefix(a, b, steps):
    for key in ['xy', 'headings']:
        assert np.array_equal(a[key][:steps+1], b[key][:steps+1]), key
    assert np.array_equal(a['velocity'][:steps], b['velocity'][:steps])
    for key in [k for k in a if k.startswith('checkpoint_')]:
        assert np.array_equal(a[key],b[key],equal_nan=True), key


def run_stage(out, c, gate_only=False):
    out.mkdir(parents=True, exist_ok=False)
    data = out/'runs'
    data.mkdir()
    metadata = dict(config=asdict(c), seeds=SEEDS, modes=MODES, horizons=[6,10.2,18],
                    evidence='exploratory development; all three seeds retained', rng='NumPy PCG64',
                    python=platform.python_version(), numpy=np.__version__,
                    branch=subprocess.check_output(['git','branch','--show-current'], cwd=ROOT, text=True).strip(),
                    base_commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
                    leader_ids_1based=[1,2,3], hidden_edges=51, total_edges=300,
                    gate_rule=dict(fresh_final_follower_error_deg=15, fresh_final_leader_error_deg=10,
                                   fresh_follower_turn_min_deg=30, drop_event_error_increase_deg=3),
                    note='Gate cutoffs are declared engineering diagnostics, not scientific thresholds. Initial ordered state is constructed. No seeds excluded. Fresh output reused across horizons.')
    source_paths=[HERE/'run_experiment.py', ROOT/'simulateEffectiveSocialInput.m',
                  ROOT/'bearingTrackStep.m',ROOT/'occlusionInputStep.m',ROOT/'buildRingConnectivity.m']
    metadata['source_sha256_at_start']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    (out/'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    records = []
    for seed in SEEDS:
        print('RUN',out.name,'seed',seed, flush=True)
        fresh = simulate(seed,c,'fresh',0)
        fresh_name = save_run(data,seed,'fresh',0,True,fresh)
        straight = simulate(seed,c,'fresh',0,False)
        straight_name = save_run(data,seed,'fresh',0,False,straight)
        for T in ([18.] if gate_only else [6.,10.2,18.]):
            duration = round(T/c.dt)
            modes = ['drop'] if gate_only else MODES[1:]
            for mode in modes:
                result = simulate(seed,c,mode,duration)
                assert_prefix(fresh,result,c.warm_steps+c.lead_steps)
                name = save_run(data,seed,mode,T,True,result)
                records.append(dict(seed=seed, turning=1, mode=mode,T=T,name=name,control=fresh_name))
            result = simulate(seed,c,'drop',duration,False)
            assert_prefix(straight,result,c.warm_steps+c.lead_steps)
            name = save_run(data,seed,'drop',T,False,result)
            records.append(dict(seed=seed,turning=0,mode='drop',T=T,name=name,control=straight_name))
        # Full no-event checks for all strategies, not merely shared prefix.
        for mode in MODES[1:]:
            zero = simulate(seed,c,mode,0)
            assert_prefix(fresh,zero,c.total_steps)
            assert np.array_equal(fresh['u'],zero['u'])
        print('PASS prefix and T=0 checks',seed,flush=True)
    write_csv(out/'manifest.csv',records)
    return records


def export_reference_request(out):
    out.mkdir(parents=True, exist_ok=False)
    c = Config(total_steps=12,warm_steps=100,lead_steps=20,cue_amplitude=0,noise_std=0)
    xy,heading,alpha,u,noise = initialize(SEEDS[0],c)
    selected=np.zeros((20,20),bool)
    for i in range(20): selected[i,(i+np.arange(1,16))%20]=True
    hidden=np.zeros_like(selected)
    result=simulate(SEEDS[0],c,'fresh',0,initial=(xy,heading,alpha,u,noise),selected=selected,hidden=hidden)
    req=dict(xy=xy.T.tolist(),heading=heading.tolist(),alpha=alpha.tolist(),u=u.tolist())
    (out/'request.json').write_text(json.dumps(req),encoding='utf-8')
    np.savez_compressed(out/'python-reference.npz',xy=result['xy'],heading=result['headings'],u=result['u'])


def validate_reference(out):
    ref=json.loads((out/'matlab-reference.json').read_text(encoding='utf-8'))
    py=np.load(out/'python-reference.npz')
    errors={
        'field_final_state_max_abs':float(np.max(np.abs(np.array(ref['u'])-py['u']))),
        'trajectory_max_abs':float(np.max(np.abs(np.stack([ref['x'],ref['y']],axis=-1).transpose(1,0,2)-py['xy']))),
        'heading_circular_max_abs':float(np.max(np.abs(wrap(np.array(ref['heading']).T-py['heading']))))}
    for k,v in errors.items(): assert v<1e-8,(k,v)
    # Compare the existing scalar KF on wrap crossing plus 8 missing samples.
    c=Config(n=2,m=40)
    selected=np.array([[False,True],[False,False]])
    alpha=np.tile(np.arange(40)*2*np.pi/40,(2,1))
    samples=ref['samples']
    for mode,matlab_mode in [('hold','hold'),('kf-mean','mean'),('kf-uncertainty','kalman')]:
        mem=Memory(selected,c)
        max_field=max_x=max_P=0.
        for ix,y in enumerate(samples):
            fresh=np.zeros((2,2),bool)
            obs=np.full((2,2),np.nan)
            if y is not None:
                fresh[0,1]=True;obs[0,1]=y
            field,d=mem.input(obs,fresh,alpha,np.zeros(2),mode)
            expected=ref['filters'][matlab_mode][ix]
            max_field=max(max_field,float(np.max(np.abs(field-np.array(expected['field'])))))
            max_x=max(max_x,float(np.max(np.abs(wrap(mem.x[0,1]-np.array(expected['x']))))))
            max_P=max(max_P,float(np.max(np.abs(mem.P[0,1]-np.array(expected['P'])))))
        errors[mode]=dict(field=max_field,x=max_x,P=max_P)
        assert max(max_field,max_x,max_P)<1e-10,errors[mode]
    (out/'validation.json').write_text(json.dumps(dict(passed=True,errors=errors),indent=2),encoding='utf-8')
    print(json.dumps(errors,indent=2))


def contracts():
    c=Config(total_steps=8,warm_steps=2,lead_steps=2)
    selected,hidden=topology(c)
    initial=initialize(SEEDS[0],c)
    # Sanitization: changing hidden truth must have no controller effect.
    mem=Memory(selected,c); xy,head,alpha,_,_=initial
    theta,_=geometry(xy,c)
    mem.input(np.where(selected,theta,np.nan),selected,alpha,head,'kf-uncertainty')
    a,b=copy.deepcopy(mem),copy.deepcopy(mem)
    fresh=selected & ~hidden
    truth_changed=theta.copy();truth_changed[hidden]=wrap(truth_changed[hidden]+2.)
    ia,da=a.input(np.where(fresh,theta,np.nan),fresh,alpha,head,'kf-uncertainty')
    ib,db=b.input(np.where(fresh,truth_changed,np.nan),fresh,alpha,head,'kf-uncertainty')
    assert np.array_equal(ia,ib) and np.array_equal(a.x,b.x)
    assert da['active'].sum()==300 and fresh.sum()==249
    # Strong uncertainty converges to a constant input while preserving mass.
    a.P[hidden,0,0]=.4
    sharp,_=b.input(np.where(fresh,theta,np.nan),fresh,alpha,head,'kf-mean')
    wide,_=a.input(np.where(fresh,theta,np.nan),fresh,alpha,head,'kf-uncertainty')
    assert np.max(np.abs(sharp.sum(0)-wide.sum(0)))<1e-8
    # Exact independent f=exp(cos-1) quadrature for posterior S.
    x=.7;V=.4;g=np.arange(100)*2*np.pi/100
    coeff=np.fft.fft(np.exp(np.cos(g)-1.))/100
    hs=np.arange(1,24)
    expected=coeff[0].real+2*np.sum(coeff[hs].real*np.exp(-.5*V*hs**2)*np.cos(hs*x))
    eps=np.linspace(-8*np.sqrt(V),8*np.sqrt(V),100001)
    quadrature=np.trapezoid(np.exp(np.cos(x+eps)-1.)*np.exp(-eps**2/(2*V))/np.sqrt(2*np.pi*V),eps)
    assert abs(expected-quadrature)<1e-10
    # Missing-before-initialization, variance/age cutoffs and reacquisition.
    empty=Memory(selected,c)
    zero,dd=empty.input(np.full_like(theta,np.nan),np.zeros_like(selected),alpha,head,'kf-mean')
    assert not zero.any() and not dd['active'].any()
    cut=copy.deepcopy(mem); cut.P[hidden,0,0]=2.
    zero,dd=cut.input(np.where(fresh,theta,np.nan),fresh,alpha,head,'kf-mean')
    assert not dd['active'][hidden].any()
    _,dd=cut.input(np.where(selected,theta,np.nan),selected,alpha,head,'kf-mean')
    assert dd['active'].sum()==300
    stale=copy.deepcopy(mem);stale.age[hidden]=121.
    _,dd=stale.input(np.where(fresh,theta,np.nan),fresh,alpha,head,'hold')
    assert not dd['active'][hidden].any()
    print('PASS: 51-edge cut; 249 fresh edges; hidden truth isolation; kernel mass; posterior S quadrature; uninitialized/cutoff/reacquisition',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('action',choices=['reference-request','reference-check','contracts','run'])
    p.add_argument('--out',type=Path)
    p.add_argument('--amplitude',type=float,default=.08)
    p.add_argument('--cue-shape',choices=['gaussian','cosine'],default='gaussian')
    p.add_argument('--lead-steps',type=int,default=20)
    p.add_argument('--steps',type=int,default=500)
    p.add_argument('--scale',type=float,default=1.)
    p.add_argument('--gate-only',action='store_true')
    args=p.parse_args()
    if args.action=='reference-request': export_reference_request(args.out)
    elif args.action=='reference-check': validate_reference(args.out)
    elif args.action=='contracts': contracts()
    else:
        contracts()
        run_stage(args.out,Config(cue_amplitude=args.amplitude,cue_shape=args.cue_shape,lead_steps=args.lead_steps,total_steps=args.steps,geometry_scale=args.scale),args.gate_only)
