"""Read the prespecified pilot design; no outcome-dependent parameters."""
import bootstrap
import numpy as np
from common import ROOT, read_json
from kf_reference_core import Config

SETTINGS = read_json(ROOT / 'experiment_config.json')
N = SETTINGS['n']
LEADER = 0
FOLLOWERS = np.arange(1, N)
STEPS = SETTINGS['steps']
DT = SETTINGS['dt']
SEEDS = tuple(SETTINGS['pilot_seeds'])
PROFILES = tuple(p['name'] for p in SETTINGS['target_profiles'])
DURATIONS = tuple(SETTINGS['durations_steps'])
STARTS = tuple(SETTINGS['occlusion_starts_zero_based'])
METHODS = tuple(SETTINGS['methods'])


def config():
    c, kf = SETTINGS['controller'], SETTINGS['kf']
    return Config(n=N, m=c['neurons'], k=7, dt=DT, arena=c['arena'],
                  beta=c['beta'], speed_scale=c['speed_scale'], h0=.012,
                  sigma=c['sigma_rad'], q=kf['q'], R=kf['R'],
                  noise_std=SETTINGS['observation_noise_std_rad'],
                  max_age=kf['max_age'], max_variance=kf['max_variance'],
                  total_steps=STEPS, cue_amplitude=0., cue_shape='cosine')


def target(profile, times=None):
    if times is None:
        times = np.arange(STEPS + 1) * DT
    times = np.asarray(times, dtype=float)
    p = next(p for p in SETTINGS['target_profiles'] if p['name'] == profile)
    if p['type'] == 'static':
        return np.broadcast_to(np.asarray(p['position'], dtype=float), times.shape + (2,)).copy()
    return np.stack((p['x0'] + p['vx'] * times,
                     p['y0'] + p['amplitude_y'] * np.sin(2 * np.pi * times / p['period'] + p['phase'])), axis=-1)


def topology():
    selected = np.zeros((N, N), dtype=bool)
    hidden = selected.copy()
    for rank, agent in enumerate(FOLLOWERS):
        peers = FOLLOWERS[(rank + np.arange(1, 7)) % len(FOLLOWERS)]
        selected[agent, LEADER] = True
        selected[agent, peers] = True
        hidden[agent, LEADER] = True
        hidden[agent, peers[:3]] = True
    return selected, hidden


def edge_scales():
    scales = np.ones((N, N))
    scales[FOLLOWERS, LEADER] = 20.
    return scales


def initial(seed):
    rng = np.random.default_rng(seed)
    xy = np.empty((N, 2))
    xy[0] = (200., 200.)
    xy[FOLLOWERS] = np.c_[rng.uniform(187., 193., len(FOLLOWERS)),
                          rng.uniform(194., 206., len(FOLLOWERS))]
    headings = rng.normal(0., .03, N)
    headings[0] = 0.
    headings[FOLLOWERS] = np.random.default_rng(seed + 1104000).uniform(-np.pi, np.pi, len(FOLLOWERS))
    alpha = (np.arange(100)[None, :] * 2 * np.pi / 100 + headings[:, None]) % (2 * np.pi)
    u = .1 * np.cos(alpha.T - headings[None, :]) + rng.normal(0., .001, (100, N))
    noise = rng.normal(0., .01, (STEPS, N, N))
    return xy, headings, alpha, u, noise


def windows(duration, start):
    steps=np.arange(STEPS)
    return (steps>=start)&(steps<start+duration)


def jobs():
    return [(s,m,a,w) for s in SEEDS for m,a,w in [('fresh',0,0)]+[(m,a,w) for a in STARTS for w in DURATIONS for m in METHODS if m!='fresh']]


def run_id(seed, method, start, duration):
    return f'{seed}_S{start:03d}_D{duration:03d}_{method}'
