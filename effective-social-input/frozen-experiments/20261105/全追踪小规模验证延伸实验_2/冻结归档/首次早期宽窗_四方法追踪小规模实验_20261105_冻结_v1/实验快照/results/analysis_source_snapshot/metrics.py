"""Post-update, time-aligned follower tracking scores; lower is better."""
import bootstrap
import numpy as np
import design as d


def series(run):
    delta = run['goals'][1:, None, :] - run['positions'][1:, d.FOLLOWERS]
    distances = np.linalg.norm(delta, axis=-1)
    velocity = run['velocities'][:, d.FOLLOWERS]
    speeds = np.linalg.norm(velocity, axis=-1)
    desired = np.arctan2(delta[..., 1], delta[..., 0])
    headings = np.arctan2(velocity[..., 1], velocity[..., 0])
    err = headings - desired
    valid = (distances > 1e-9) & (speeds > 1e-12)
    angle = np.where(valid, np.degrees(np.abs(np.arctan2(np.sin(err), np.cos(err)))), np.nan)
    count = valid.sum(1)
    mean_angle = np.divide(np.nansum(angle, axis=1), count,
                           out=np.full(d.STEPS, np.nan), where=count > 0)
    unit = np.divide(velocity, speeds[..., None], out=np.zeros_like(velocity), where=speeds[..., None] > 1e-12)
    centroid = run['positions'][1:, d.FOLLOWERS].mean(1)
    leader_delta = run['goals'][1:] - run['positions'][1:, 0]
    return dict(angle_deg=mean_angle, distance=distances.mean(1),
                angle_valid_fraction=valid.mean(1), stationary_fraction=(speeds <= 1e-12).mean(1),
                leader_target_distance=np.linalg.norm(leader_delta, axis=1),
                mean_follower_speed=speeds.mean(1),
                follower_GO=np.linalg.norm(unit.sum(1), axis=1) / len(d.FOLLOWERS),
                follower_spread=np.linalg.norm(run['positions'][1:, d.FOLLOWERS] - centroid[:, None], axis=-1).mean(1),
                fresh_edges=run['fresh'].sum((1, 2)), active_edges=run['active'].sum((1, 2)),
                mean_social_mass=run['input_mass'][:, d.FOLLOWERS].mean(1),
                mean_effective_social_input=run['effective_social_input'][:, d.FOLLOWERS].mean(1))
