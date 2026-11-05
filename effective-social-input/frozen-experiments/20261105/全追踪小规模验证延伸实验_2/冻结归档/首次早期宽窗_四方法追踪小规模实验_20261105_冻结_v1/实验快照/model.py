"""Closed-loop followers; the leader alone knows the current target position."""
import bootstrap
import numpy as np
import design as d
import kf_reference_core as ref


def simulate(seed, profile, method, duration, start):
    c = d.config()
    xy, heading, alpha, u, noise = d.initial(seed)
    initial_xy, initial_u = xy.copy(), u.copy()
    goals = d.target(profile)
    selected, hidden = d.topology()
    scales = d.edge_scales()
    memory = ref.Memory(selected, c, scales)
    ring = ref.ring(c.m)
    blocked = d.windows(duration, start) if method != 'fresh' else np.zeros(d.STEPS, dtype=bool)
    positions = np.empty((d.STEPS + 1, d.N, 2))
    headings = np.empty((d.STEPS + 1, d.N))
    velocities = np.empty((d.STEPS, d.N, 2))
    observations = np.full((d.STEPS, d.N, d.N), np.nan)
    centers = np.empty_like(observations)
    active = np.empty(observations.shape, bool)
    fresh_log = np.empty_like(active)
    filter_state = np.empty(observations.shape + (2,))
    filter_covariance = np.empty(observations.shape + (2, 2))
    input_mass = np.empty((d.STEPS, d.N))
    effective_social_input = np.empty_like(input_mass)
    positions[0], headings[0] = xy, heading
    for step in range(d.STEPS):
        theta, _ = ref.geometry(xy, c)
        available = selected & ~hidden if blocked[step] else selected.copy()
        obs = np.where(available, (theta + noise[step]) % (2 * np.pi), np.nan)
        # No target, position, hidden truth or future data enters Memory.input.
        field, diag = memory.input(obs, available, alpha, heading, method)
        uf = u[:, d.FOLLOWERS]
        uf = uf + c.dt * (-uf + ring @ np.tanh(c.beta * uf) / c.m + field[:, d.FOLLOWERS])
        u[:, d.FOLLOWERS] = uf
        positive = np.maximum(np.tanh(c.beta * uf), 0.)
        vf = c.speed_scale * np.c_[np.sum(positive * np.cos(alpha[d.FOLLOWERS].T), 0),
                                   np.sum(positive * np.sin(alpha[d.FOLLOWERS].T), 0)]
        # Leader uses the current target only, with no look-ahead or social term.
        remaining = goals[step] - xy[d.LEADER]
        length = np.linalg.norm(remaining)
        vl = (remaining / length * min(1.5, length / c.dt)
              if length > 1e-12 else np.zeros(2))
        velocity = np.empty((d.N, 2))
        velocity[d.FOLLOWERS], velocity[d.LEADER] = vf, vl
        speeds = np.linalg.norm(velocity, axis=1)
        moving = speeds > 1e-12
        heading[moving] = np.arctan2(velocity[moving, 1], velocity[moving, 0])
        xy = xy + c.dt * velocity
        positions[step + 1], headings[step + 1], velocities[step] = xy, heading, velocity
        observations[step], centers[step] = obs, diag['center']
        active[step], fresh_log[step] = diag['active'], available
        filter_state[step], filter_covariance[step] = memory.x, memory.P
        input_mass[step], effective_social_input[step] = diag['A'], diag['S']
    assert np.isfinite(positions).all() and np.isfinite(velocities).all()
    assert np.all(positions >= 0.) and np.all(positions < c.arena), 'Boundary crossing'
    return dict(seed=np.array(seed), profile=np.array(profile), method=np.array(method),
                duration=np.array(duration), occlusion_start=np.array(start), positions=positions, headings=headings,
                velocities=velocities, goals=goals, alpha=alpha, initial_xy=initial_xy,
                initial_u=initial_u, final_u=u, noise=noise, selected=selected, hidden=hidden,
                blocked=blocked, observations=observations, centers=centers, active=active,
                fresh=fresh_log, edge_scales=scales, filter_state=filter_state,
                filter_covariance=filter_covariance, input_mass=input_mass,
                effective_social_input=effective_social_input)
