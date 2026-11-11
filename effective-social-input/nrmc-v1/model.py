"""Allocentric neural field, matching simulateEffectiveSocialInput.m.

Agent-major arrays are used here. The original implementation remains read-only.
"""
import numpy as np


def wrap(x):
    return (x + np.pi) % (2*np.pi) - np.pi


def min_image(x, length):
    q = x / length
    return x - length * np.sign(q) * np.floor(np.abs(q) + .5)


class Field:
    def __init__(self, cfg):
        self.c = cfg
        self.n, self.m = cfg['N'], cfg['M']
        self.ids = (np.arange(self.n)[:, None] + np.arange(1, self.n)[None]) % self.n
        d = np.arange(self.m) * 2*np.pi/self.m
        d = np.minimum(d, 2*np.pi-d)
        self.wfft = np.fft.fft(np.cos(np.pi*(d/np.pi)**cfg['connection_exponent']))
        self.freq = np.fft.fftfreq(self.m, 1/self.m)

    def initial(self, seed):
        r = np.random.default_rng(np.random.SeedSequence([seed, 0, self.n]))
        pos = r.random((2, self.n)).T * self.c['L']
        heading = np.zeros(self.n); u = np.zeros((self.n, self.m))
        for i in range(self.n):
            heading[i] = 2*np.pi*r.random()
            u[i] = .1*r.standard_normal(self.m)
        alpha = (2*np.pi*np.arange(self.m)[None]/self.m + heading[:, None]) % (2*np.pi)
        return dict(position=pos, u=u, heading=heading, alpha=alpha)

    def bearings(self, state):
        d = min_image(state['position'][self.ids]-state['position'][:, None], self.c['L'])
        return np.arctan2(d[..., 1], d[..., 0])

    def kernels(self, alpha, centers, variance=None):
        # centers: (..., N, K); result (..., N, K, M).
        d = wrap(alpha[:, None, :] - centers[..., None])
        k = np.exp(-.5*(d/self.c['kernel_sigma'])**2)
        if variance is not None:
            k = np.fft.ifft(np.fft.fft(k, axis=-1) *
                            np.exp(-.5*variance[..., None]*self.freq**2), axis=-1).real
            k = np.maximum(k, 0)
        return k

    def observed_input(self, alpha, y, fresh):
        if not np.isfinite(y[fresh]).all() or not np.isnan(y[~fresh]).all():
            raise ValueError('Observed bearings must be finite; hidden bearings must be NaN')
        return self.c['edge_weight'] * (self.kernels(alpha, np.where(fresh, y, 0)) *
                                        fresh[..., None]).sum(axis=-2)

    def base(self, u, visible_input, external=None):
        recurrence = np.fft.ifft(np.fft.fft(np.tanh(self.c['beta']*u), axis=-1)*self.wfft,
                                 axis=-1).real/self.m
        b = u + self.c['dt']*(-u+recurrence-self.c['global_inhibition']+visible_input)
        return b if external is None else b+self.c['dt']*external

    def response(self, base, inputs, alpha):
        # inputs may contain a candidate/sample axis.
        multi = inputs.ndim == 3
        a = alpha[:, None, :] if multi else alpha
        b = base[:, None, :] if multi else base
        out = np.maximum(0, np.tanh(self.c['beta']*(b+self.c['dt']*inputs)))
        return self.c['speed_scale']*np.stack([(out*np.cos(a)).sum(axis=-1),
                                              (out*np.sin(a)).sum(axis=-1)], axis=-1)

    def advance(self, state, total_input, external=None):
        state['u'] = self.base(state['u'], total_input, external)
        v = self.response(state['u'], np.zeros_like(total_input), state['alpha'])
        moving = np.abs(v).sum(axis=-1) >= self.c['speed_scale']*1e-9
        state['heading'][moving] = np.arctan2(v[moving, 1], v[moving, 0]) % (2*np.pi)
        state['position'] = (state['position']+self.c['dt']*v) % self.c['L']
        return v

    def metrics(self, state, v):
        speed = np.linalg.norm(v, axis=1)
        unit = np.divide(v, speed[:, None], out=np.zeros_like(v), where=speed[:, None] > 1e-12)
        d = min_image(state['position'][:, None]-state['position'][None], self.c['L'])
        distance = np.linalg.norm(d, axis=-1)
        np.fill_diagonal(distance, np.inf)
        ids = np.argsort(distance, axis=1, kind='stable')[:, :5]
        go = np.linalg.norm(unit.mean(axis=0))
        lo = np.linalg.norm(unit+unit[ids].sum(axis=1), axis=1).mean()/6
        return np.array([go, lo, min(go, lo), speed.mean()])


def clone(state):
    return {k: v.copy() for k, v in state.items()}
