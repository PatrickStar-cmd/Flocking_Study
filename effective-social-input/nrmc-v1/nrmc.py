"""NRMC-v1: causal bearing filter bank and one-step response selection."""
import copy
import numpy as np
from scipy.special import logsumexp
from model import wrap


def log_wrapped_normal(e, variance):
    v = np.asarray(variance)
    out = np.empty_like(e)
    narrow = v < 1
    a, b = e[narrow], v[narrow]
    out[narrow] = logsumexp(-.5*(a[:, None]+2*np.pi*np.arange(-1, 2))**2/b[:, None], axis=1) - .5*np.log(2*np.pi*b)
    a, b = e[~narrow], v[~narrow]
    f = np.arange(1, 13)
    out[~narrow] = np.log((1+2*(np.exp(-.5*b[:, None]*f**2)*np.cos(a[:, None]*f)).sum(axis=1))/(2*np.pi))
    return out


class Bank:
    def __init__(self, shape, cfg, dt):
        self.cfg = cfg; self.R = max(cfg['R'], 1e-8)
        rho = np.exp(-dt/cfg['tau']); tau = cfg['tau']
        self.F = np.array([[[1, 0], [0, 0]], [[1, dt], [0, 1]],
                           [[1, tau*(1-rho)], [0, rho]]], float)
        self.Q = np.zeros((3, 2, 2))
        self.Q[0, 0, 0] = cfg['q_H']*dt
        self.Q[1] = cfg['q_C']*np.array([[dt**3/3, dt**2/2], [dt**2/2, dt]])
        self.Q[2] = cfg['q_D']*np.array([[tau*tau*dt-tau**3*(2*(1-rho)-(1-rho*rho)/2), tau*tau*(1-rho)**2/2],
                                        [tau*tau*(1-rho)**2/2, tau*(1-rho*rho)/2]])
        self.x = np.zeros((3, *shape, 2)); self.P = np.zeros((3, *shape, 2, 2))
        self.pi = np.full((3, *shape), 1/3)
        self.seen = np.zeros(shape, bool); self.last = np.full(shape, np.nan)

    def copy(self):
        return copy.deepcopy(self)

    def update(self, y, fresh):
        if not np.isnan(y[~fresh]).all() or not np.isfinite(y[fresh]).all():
            raise ValueError('Unsanitized observation')
        old = self.seen.copy(); new = fresh & ~old
        self.x = np.einsum('mab,mijb->mija', self.F, self.x)
        self.P = np.einsum('mab,mijbc,mdc->mijad', self.F, self.P, self.F)+self.Q[:, None, None]
        upd = fresh & old
        if upd.any():
            e = wrap(y[upd][None]-self.x[:, upd, 0])
            p = self.P[:, upd].copy()
            S = p[..., 0, 0]+self.R
            prior = (1-self.cfg['epsilon'])*self.pi[:, upd]+self.cfg['epsilon']/3
            lp = np.log(prior)+log_wrapped_normal(e, S)
            self.pi[:, upd] = np.exp(lp-logsumexp(lp, axis=0, keepdims=True))
            gain = p[..., :, 0]/S[..., None]
            self.x[:, upd] += gain*e[..., None]
            A = np.broadcast_to(np.eye(2), p.shape).copy(); A[..., :, 0] -= gain
            self.P[:, upd] = A@p@np.swapaxes(A, -1, -2)+self.R*gain[..., :, None]*gain[..., None, :]
        if new.any():
            self.x[:, new, 0] = y[new]
            self.x[:, new, 1] = 0
            self.P[:, new] = 0
            self.P[:, new, 0, 0] = self.R
            self.P[1:, new, 1, 1] = 1
        self.x[..., 0] = wrap(self.x[..., 0])
        self.P = .5*(self.P+np.swapaxes(self.P, -1, -2))
        self.last[fresh] = y[fresh]; self.seen |= fresh


def candidates(fields, gamma):
    # gamma-major ordering implements the smallest-gamma tie rule.
    parts = [np.zeros_like(fields[:, 0])]
    labels = [(-1, 0.)]
    for g in gamma[1:]:
        for b in range(fields.shape[1]):
            parts.append(g*fields[:, b]); labels.append((b, float(g)))
    return np.stack(parts, axis=1), np.asarray(labels)


def choose(field, base, alpha, inputs, scenarios):
    target = field.response(base, scenarios, alpha).mean(axis=1)
    predicted = field.response(base, inputs, alpha)
    loss = ((predicted-target[:, None])**2).sum(axis=-1)
    index = np.argmin(loss, axis=1)
    return inputs[np.arange(field.n), index], index, loss


def compensate(method, bank, y, fresh, u, alpha, field, cfg, random_key, external=None):
    """No positions, true hidden bearings, GO, or LO enter this interface."""
    visible = field.observed_input(alpha, y, fresh)
    hidden = ~fresh & bank.seen
    diag = dict(gamma=0., component=-1., selected_risk=0., mix_risk=0.)
    if method in ('drop', 'fresh') or not hidden.any():
        return visible, diag
    if method == 'hold':
        k = field.kernels(alpha, np.nan_to_num(bank.last))
        return visible+field.c['edge_weight']*(k*hidden[..., None]).sum(axis=1), diag
    if method in ('cv-kf', 'd-kf'):
        m = 1 if method == 'cv-kf' else 2
        k = field.kernels(alpha, bank.x[m, ..., 0], bank.P[m, ..., 0, 0])
        return visible+field.c['edge_weight']*(k*hidden[..., None]).sum(axis=1), diag
    if method not in ('nrmc', 'nrmc-mix', 'nrmc-d'):
        raise ValueError(method)
    k = field.kernels(alpha, bank.x[..., 0], bank.P[..., 0, 0])
    single = field.c['edge_weight']*(k*hidden[None, ..., None]).sum(axis=2)
    mix = field.c['edge_weight']*(k*bank.pi[..., None]*hidden[None, ..., None]).sum(axis=(0, 2))
    if method == 'nrmc-mix':
        return visible+mix, diag
    fields = np.concatenate([single.transpose(1, 0, 2), mix[:, None]], axis=1)
    if method == 'nrmc-d':
        fields = fields[:, 2:3]
    actions, labels = candidates(fields, cfg['gamma'])
    r = np.random.default_rng(np.random.SeedSequence(random_key))
    shape = (field.n, cfg['samples'], field.n-1)
    categorical = r.random(shape)
    cum = np.cumsum(bank.pi, axis=0).transpose(1, 2, 0)
    component = (categorical[..., None] > cum[:, None, :, :2]).sum(axis=-1)
    mu = bank.x[..., 0].transpose(1, 2, 0)
    var = bank.P[..., 0, 0].transpose(1, 2, 0)
    ii = np.arange(field.n)[:, None, None]; jj = np.arange(field.n-1)[None, None, :]
    theta = mu[ii, jj, component]+np.sqrt(var[ii, jj, component])*r.standard_normal(shape)
    scenarios = np.zeros((field.n, cfg['samples'], field.m))
    # A per-edge loop limits peak memory while preserving exact Gaussian kernels.
    for j in range(field.n-1):
        d = wrap(alpha[:, None]-theta[:, :, j, None])
        scenarios += hidden[:, j, None, None]*np.exp(-.5*(d/field.c['kernel_sigma'])**2)
    scenarios *= field.c['edge_weight']
    base = field.base(u, visible, external)
    selected, ix, loss = choose(field, base, alpha, actions, scenarios)
    active = hidden.any(axis=1)
    diag = dict(gamma=float(labels[ix[active], 1].mean()),
                component=float(labels[ix[active], 0].mean()),
                selected_risk=float(loss[np.arange(field.n), ix][active].mean()),
                mix_risk=float(loss[:, -1][active].mean()))
    return visible+selected, diag
