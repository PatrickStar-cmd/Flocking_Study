"""Target-tracking and repeated-short-occlusion experiment model.

The bearing KF and neural-field implementation come from the validated
reference_run_experiment module. This experiment fixes N=30 and L=6,
uses a moving target, and applies three short, separated occlusion windows.
"""
from __future__ import annotations

import numpy as np

import reference_run_experiment as ref


N = 30
LEADER_COUNT = 6
K_NEIGHBORS = 15
HIDDEN_PER_FOLLOWER = 9
OCCLUSION_WINDOWS = ((180, 200), (300, 320), (420, 440))
TOTAL_STEPS = 600
NOISE_LEVELS = (0.01, 0.035)
METHODS = ("drop", "hold", "kf-mean")
CONDITIONS = (("tracking", True),)
MASKS = ("leader_loss",)


def config(noise_std: float = 0.01, total_steps: int = TOTAL_STEPS):
    return ref.Config(
        n=N,
        k=K_NEIGHBORS,
        lead_steps=20,
        total_steps=total_steps,
        cue_amplitude=0.32,
        cue_shape="cosine",
        noise_std=noise_std,
        R=noise_std ** 2,
    )


def role_ids(seed: int):
    permutation = np.random.default_rng(seed + 1_406_117).permutation(N)
    return permutation[:LEADER_COUNT], permutation[LEADER_COUNT:]


def initialize(seed: int, cfg):
    """Generate common initial state and a standard-normal noise tensor scaled by cfg."""
    base_cfg = config(noise_std=1.0, total_steps=cfg.total_steps)
    xy, heading, alpha, u, standard_noise = ref.initialize(seed, base_cfg)
    return xy, heading, alpha, u, standard_noise * cfg.noise_std


def _is_strongly_connected(adjacency: np.ndarray) -> bool:
    n = len(adjacency)

    def reached(matrix):
        seen = {0}
        stack = [0]
        while stack:
            node = stack.pop()
            for nxt in np.flatnonzero(matrix[node]):
                nxt = int(nxt)
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return len(seen) == n

    return reached(adjacency) and reached(adjacency.T)


def topology(cfg, seed: int, mask_type: str):
    """Build a fixed-degree graph and a matched nine-edge mask."""
    if mask_type not in MASKS:
        raise ValueError(f"Unknown mask type: {mask_type}")
    leaders, followers = role_ids(seed)
    follower_pos = {int(agent): pos for pos, agent in enumerate(followers)}
    selected = np.zeros((cfg.n, cfg.n), dtype=bool)

    follower_links = cfg.k - len(leaders)
    for agent in followers:
        pos = follower_pos[int(agent)]
        selected[agent, leaders] = True
        ring_ids = followers[(pos + np.arange(1, follower_links + 1)) % len(followers)]
        selected[agent, ring_ids] = True

    for rank, agent in enumerate(leaders):
        selected[agent, leaders[leaders != agent]] = True
        needed = cfg.k - (len(leaders) - 1)
        follower_ids = followers[(rank + np.arange(needed)) % len(followers)]
        selected[agent, follower_ids] = True

    hidden = np.zeros_like(selected)
    for agent in followers:
        candidates = np.flatnonzero(selected[agent] & np.isin(np.arange(cfg.n), followers))
        rng = np.random.default_rng(seed * 31 + int(agent) * 101 + 4_019)
        ordered = rng.permutation(candidates)
        if mask_type == "leader_loss":
            hidden[agent, leaders] = True
            hidden[agent, ordered[: HIDDEN_PER_FOLLOWER - len(leaders)]] = True
        else:
            hidden[agent, ordered[:HIDDEN_PER_FOLLOWER]] = True

    assert np.all(selected.sum(axis=1) == cfg.k)
    assert not np.any(np.diag(selected))
    assert np.all(hidden.sum(axis=1)[followers] == HIDDEN_PER_FOLLOWER)
    assert not np.any(hidden[leaders])
    if mask_type == "leader_loss":
        assert np.all(hidden[np.ix_(followers, leaders)])
    else:
        assert not np.any(hidden[np.ix_(followers, leaders)])
    assert _is_strongly_connected(selected)
    return selected, hidden, leaders, followers


def blackout_at(step: int) -> bool:
    return any(start <= step < end for start, end in OCCLUSION_WINDOWS)


def target_position(step: int, tracking: bool) -> np.ndarray:
    """A moving target: straight baseline or a finite smooth lateral route."""
    u = min(1.0, max(0.0, step / float(TOTAL_STEPS - 1)))
    x = 620.0 + 220.0 * u
    y = 500.0 + (60.0 * np.sin(2.0 * np.pi * u) if tracking else 0.0)
    return np.array([x, y], dtype=float)


class ForecastMemory(ref.Memory):
    """Ten-sample circular linear-trend baseline retained from experiment_3."""

    def __init__(self, selected, cfg):
        super().__init__(selected, cfg)
        self.history = np.full((10, *selected.shape), np.nan)
        self.count = np.zeros(selected.shape, dtype=int)
        self.previous_fresh = np.zeros_like(selected)
        self.slope = np.zeros(selected.shape)

    def input(self, observed, fresh, alpha, heading, mode):
        if mode != "linear-trend":
            return super().input(observed, fresh, alpha, heading, mode)

        restarted = fresh & ~self.previous_fresh
        self.history[:, restarted] = np.nan
        self.count[restarted] = 0
        self.history[1:, fresh] = self.history[:-1, fresh]
        self.history[0, fresh] = observed[fresh]
        self.count[fresh] = np.minimum(self.count[fresh] + 1, 10)

        dif = ref.wrap(self.history[:-1] - self.history[1:])
        unwrapped = np.concatenate(
            [self.history[:1], self.history[:1] - np.cumsum(dif, axis=0)], axis=0
        )
        valid = np.arange(10)[:, None, None] < self.count[None]
        y = np.where(valid, unwrapped, 0.0)
        t = -np.arange(10)[:, None, None] * self.cfg.dt
        count = np.maximum(self.count, 1)
        st = np.sum(t * valid, axis=0)
        sy = np.sum(y, axis=0)
        numerator = np.sum(t * y, axis=0) - st * sy / count
        denominator = np.sum(t * t * valid, axis=0) - st * st / count
        self.slope = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
        self.previous_fresh = fresh.copy()

        _, details = super().input(observed, fresh, alpha, heading, "hold")
        active = details["active"]
        safe_age = np.where(self.initialized, self.age, 0.0)
        center = np.where(fresh, observed, (self.last + safe_age * self.slope) % (2 * np.pi))
        center = np.where(active, center, 0.0)
        weights = details["weights"]
        kernels = np.exp(-0.5 * (ref.wrap(alpha[:, None, :] - center[:, :, None]) / self.cfg.sigma) ** 2)
        field = np.sum(weights[:, :, None] * kernels, axis=1).T
        factor = np.exp(np.cos(center - heading[:, None]) - 1.0)
        social_input = np.sum(weights * factor, axis=1)
        weight_sum = np.sum(weights, axis=1)
        details.update(
            center=center,
            S=social_input,
            A=weight_sum,
            C=np.divide(social_input, weight_sum, out=np.zeros_like(social_input), where=weight_sum > 0),
            S_memory=np.sum(weights * factor * (active & ~fresh), axis=1),
        )
        return field, details


def simulate(seed: int, mode: str, tracking: bool, mask_type: str, initial=None, cfg=None):
    cfg = config() if cfg is None else cfg
    if initial is None:
        initial = initialize(seed, cfg)
    xy, heading, alpha, u, noise = [np.array(value, copy=True) for value in initial]
    selected, hidden, leaders, followers = topology(cfg, seed, mask_type)
    memory = ForecastMemory(selected, cfg)
    ring = ref.ring(cfg.m)
    rows = []
    fixed_nearest = None
    theta_before_window = None

    for step in range(cfg.total_steps):
        theta, distance = ref.geometry(xy, cfg)
        if step == OCCLUSION_WINDOWS[0][0]:
            d0 = distance.copy()
            np.fill_diagonal(d0, np.inf)
            fixed_nearest = np.argsort(d0, axis=1)[:, :5]
        if any(step == start - 1 for start, _ in OCCLUSION_WINDOWS):
            theta_before_window = theta.copy()

        blackout = mode != "fresh" and blackout_at(step)
        fresh = selected & ~(hidden if blackout else np.zeros_like(hidden))
        observed = np.where(fresh, (theta + noise[step]) % (2 * np.pi), np.nan)
        field, details = memory.input(observed, fresh, alpha, heading, mode)

        target = target_position(step, tracking)
        target_heading = np.arctan2(target[1] - xy[:, 1], target[0] - xy[:, 0]) % (2 * np.pi)
        cue = cfg.cue_amplitude * np.cos(alpha[leaders] - target_heading[leaders, None])
        field[:, leaders] += cue.T

        oracle = cfg.h0 * np.sum(selected * np.exp(np.cos(theta - heading[:, None]) - 1.0), axis=1)
        prediction_error = ref.wrap(details["center"] - theta)
        prediction_mask = hidden & blackout & details["active"]

        u += cfg.dt * (-u + ring @ np.tanh(cfg.beta * u) / cfg.m + field)
        positive = np.maximum(np.tanh(cfg.beta * u), 0.0)
        velocity = cfg.speed_scale * np.c_[
            np.sum(positive * np.cos(alpha.T), axis=0),
            np.sum(positive * np.sin(alpha.T), axis=0),
        ]
        speed = np.linalg.norm(velocity, axis=1)
        moving = speed > 1e-12
        direction = np.divide(velocity, speed[:, None], out=np.zeros_like(velocity), where=moving[:, None])
        heading[moving] = np.arctan2(velocity[moving, 1], velocity[moving, 0]) % (2 * np.pi)
        xy = (xy + cfg.dt * velocity) % cfg.arena

        _, after = ref.geometry(xy, cfg)
        no_self = after.copy()
        np.fill_diagonal(no_self, np.inf)
        leader_direction = np.angle(np.mean(np.exp(1j * heading[leaders])))
        follower_direction = np.angle(np.mean(np.exp(1j * heading[followers])))
        target_error = np.degrees(np.abs(ref.wrap(heading - target_heading)))
        pair_indices = np.triu_indices(cfg.n, 1)
        rows.append(dict(
            step=step,
            n=cfg.n,
            leaders=len(leaders),
            noise_std=float(cfg.noise_std),
            R=float(cfg.R),
            target_x=float(target[0]),
            target_y=float(target[1]),
            target_heading_mean_deg=float(np.degrees(np.angle(np.mean(np.exp(1j * target_heading))))),
            blackout=int(blackout),
            blackout_window=int(next((i + 1 for i, (start, end) in enumerate(OCCLUSION_WINDOWS) if start <= step < end), 0)),
            GO_all=float(np.linalg.norm(direction.sum(axis=0)) / cfg.n),
            GO_followers=float(np.linalg.norm(direction[followers].sum(axis=0)) / len(followers)),
            GO_leaders=float(np.linalg.norm(direction[leaders].sum(axis=0)) / len(leaders)),
            follower_target_error_deg=float(target_error[followers].mean()),
            leader_target_error_deg=float(target_error[leaders].mean()),
            follower_leader_error_deg=float(np.degrees(np.abs(ref.wrap(heading[followers] - leader_direction))).mean()),
            follower_mean_direction_deg=float(np.degrees(follower_direction)),
            leader_mean_direction_deg=float(np.degrees(leader_direction)),
            pair_distance=float(after[pair_indices].mean()),
            fixed5_distance=float(after[np.arange(cfg.n)[:, None], fixed_nearest].mean()),
            minimum_distance=float(no_self.min()),
            S_abs_error=float(np.abs(details["S"][followers] - oracle[followers]).mean()),
            fresh_edges=int(fresh.sum()),
            hidden_edges=int((hidden & blackout).sum()),
            hidden_active_edges=int(prediction_mask.sum()),
            hidden_prediction_error_deg=(float(np.degrees(np.abs(prediction_error[prediction_mask])).mean()) if prediction_mask.any() else np.nan),
            hidden_bearing_change_deg=(float(np.degrees(np.abs(ref.wrap(theta - theta_before_window))[hidden]).mean()) if theta_before_window is not None and blackout else np.nan),
        ))

    return rows

