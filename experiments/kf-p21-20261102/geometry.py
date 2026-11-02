"""experiments/kf-p21-20261102/geometry.py

Geometric decoupling and classification module for P2.1.
Implements exact minimum-image calculations, branch jump detection,
window-level and time-step-level classification, and continuous angular rate estimation.

Standalone geometry library: strictly isolated from all filter / estimation logic.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple, Union
import numpy as np

# Window classification constants
LABEL_CLEAN: str = "G-clean"
LABEL_SWITCH: str = "G-switch"
LABEL_CLOSE: str = "G-close"
LABEL_AMBIGUOUS: str = "G-ambiguous"
WINDOW_LABELS: Tuple[str, ...] = (
    LABEL_CLEAN,
    LABEL_SWITCH,
    LABEL_CLOSE,
    LABEL_AMBIGUOUS,
)

# Time-step level fresh observation classification constants
FRESH_SWITCH_EVENT: str = "switch_event"
FRESH_SWITCH_NEAR2: str = "switch_near2"
FRESH_CLOSE: str = "G-close"
FRESH_CLEAN: str = "G-clean"
FRESH_AMBIGUOUS: str = "G-ambiguous"
FRESH_LABELS: Tuple[str, ...] = (
    FRESH_SWITCH_EVENT,
    FRESH_SWITCH_NEAR2,
    FRESH_CLOSE,
    FRESH_CLEAN,
    FRESH_AMBIGUOUS,
)


def round_matlab(x: Union[float, int, np.ndarray]) -> Union[float, np.ndarray]:
    """MATLAB-compatible round function: sign(x) * floor(abs(x) + 0.5).

    Rounds half-way cases away from zero (e.g. 0.5 -> 1.0, -0.5 -> -1.0, 1.5 -> 2.0).
    """
    arr = np.asarray(x, dtype=np.float64)
    # np.sign(0.0) returns 0.0
    res = np.sign(arr) * np.floor(np.abs(arr) + 0.5)
    if np.isscalar(x):
        return float(res.item())
    return res


def wrap(angle: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Wrap angle differences to [-pi, pi) using atan2(sin, cos)."""
    return np.arctan2(np.sin(angle), np.cos(angle))


def geometry(
    source: Dict[str, Any],
    *,
    arena_size: float = 1000.0,
    dt: float = 0.3,
) -> Dict[str, np.ndarray]:
    """Compute exact minimum-image relative vectors, branch jumps, and geometry metrics.

    Parameters
    ----------
    source : dict
        Loaded from scipy.io.loadmat(simplify_cells=True). Required keys:
          - 'positions': [2, N, T] finite coordinates in [0, arena_size)
          - 'true_bearing': [T, E] ground truth bearings in radians
          - 'observer_id': [E] 1-based observer agent indices
          - 'target_id': [E] 1-based target agent indices
        Optional keys:
          - 'arena_size': float, default 1000.0
          - 'dt': float, default 0.3
    arena_size : float
        Periodic domain edge length (default 1000.0).
    dt : float
        Sampling step interval in model time (default 0.3).

    Returns
    -------
    dict with keys:
      - 'distance': [T, E] float64 Euclidean distance ||d_min||
      - 'image_switch': [T, E] bool flag for nearest-periodic-image branch jump (> arena_size/2)
      - 'near_range': [T, E] bool flag for close proximity (r < 5.0)
      - 'near_half_box': [T, E] bool flag for proximity to periodic boundary choice (margin < 10.0)
      - 'near_aliasing': [T, E] bool flag for sampling aliasing (|wrap(Delta theta)| >= 0.9*pi)
      - 'angular_rate_clean': [T, E] float64 continuous angular rate in rad/model time (NaN on invalid/first step)
      - 'd_min': [T, E, 2] float64 minimum-image relative displacement vector
      - 'image_index': [T, E, 2] float64 integer image offset k such that d_min = d_raw - L*k
    """
    if "positions" not in source:
        raise KeyError("source dict must contain 'positions'")
    if "true_bearing" not in source:
        raise KeyError("source dict must contain 'true_bearing'")
    if "observer_id" not in source or "target_id" not in source:
        raise KeyError("source dict must contain 'observer_id' and 'target_id'")

    pos = np.asarray(source["positions"], dtype=np.float64)
    if pos.ndim != 3 or pos.shape[0] != 2:
        raise ValueError(f"positions must have shape [2, N, T], got {pos.shape}")
    if not np.all(np.isfinite(pos)):
        raise ValueError("positions array contains non-finite values")

    _, N, T = pos.shape
    if T < 1:
        raise ValueError("Trajectory length T must be >= 1")

    tb = np.asarray(source["true_bearing"], dtype=np.float64)
    if tb.ndim != 2 or tb.shape[0] != T:
        raise ValueError(f"true_bearing must have shape [T, E] with T={T}, got {tb.shape}")
    E = tb.shape[1]

    obs_id = np.asarray(source["observer_id"]).reshape(-1)
    tgt_id = np.asarray(source["target_id"]).reshape(-1)
    if len(obs_id) != E or len(tgt_id) != E:
        raise ValueError(f"observer_id and target_id must have length E={E}")

    obs_idx = obs_id.astype(int) - 1
    tgt_idx = tgt_id.astype(int) - 1
    if np.any(obs_idx < 0) or np.any(obs_idx >= N):
        raise ValueError(f"observer_id indices must be in 1-based range [1, {N}]")
    if np.any(tgt_idx < 0) or np.any(tgt_idx >= N):
        raise ValueError(f"target_id indices must be in 1-based range [1, {N}]")

    # Use explicit argument if provided, else check source dict
    L = float(source.get("arena_size", arena_size))
    step_dt = float(source.get("dt", dt))
    if L <= 0.0:
        raise ValueError(f"arena_size must be positive, got {L}")
    if step_dt <= 0.0:
        raise ValueError(f"dt must be positive, got {step_dt}")

    # Extract positions per edge:
    # pos is [2, N, T]. pos[:, obs_idx, :] is [2, E, T]. Transpose to [T, E, 2].
    p_obs = pos[:, obs_idx, :].transpose(2, 1, 0)  # [T, E, 2]
    p_tgt = pos[:, tgt_idx, :].transpose(2, 1, 0)  # [T, E, 2]

    # Raw relative displacement
    d_raw = p_tgt - p_obs  # [T, E, 2]

    # Minimum image calculation with exact MATLAB rounding rule
    image_index = round_matlab(d_raw / L)  # [T, E, 2]
    d_min = d_raw - L * image_index  # [T, E, 2]

    # Distance
    distance = np.linalg.norm(d_min, axis=2)  # [T, E]

    # Nearest periodic image branch jump:
    # Defined on transitions t-1 -> t: |d_min,k(t) - d_min,k(t-1)| > L / 2
    # First step (t=0) has no derivative -> False.
    image_switch = np.zeros((T, E), dtype=bool)
    if T > 1:
        diff_dmin = np.abs(d_min[1:] - d_min[:-1])  # [T-1, E, 2]
        image_switch[1:] = np.any(diff_dmin > (L / 2.0), axis=2)

    # Diagnostic geometric flags
    near_range = distance < 5.0  # [T, E] bool

    # near_half_box: min_k ||d_min,k| - L/2| < 10.0
    half_box_margin = np.abs(np.abs(d_min) - (L / 2.0))  # [T, E, 2]
    min_margin = np.min(half_box_margin, axis=2)  # [T, E]
    near_half_box = min_margin < 10.0  # [T, E] bool

    # near_aliasing: |wrap(Delta theta)| >= 0.9 * pi
    # First step (t=0) has no derivative -> False.
    near_aliasing = np.zeros((T, E), dtype=bool)
    if T > 1:
        diff_tb = wrap(tb[1:] - tb[:-1])  # [T-1, E]
        near_aliasing[1:] = np.abs(diff_tb) >= (0.9 * np.pi)

    # Continuous clean angular rate:
    # wrap(Delta theta) / dt, valid only when:
    # 1) t >= 1 (first step no derivative)
    # 2) no image_switch at transition (image_switch[t] is False)
    # 3) neither adjacent point is in near_range (near_range[t] and near_range[t-1] are both False)
    angular_rate_clean = np.full((T, E), np.nan, dtype=np.float64)
    if T > 1:
        valid_diff = (~image_switch[1:]) & (~near_range[1:]) & (~near_range[:-1])
        angular_rate_clean[1:] = np.where(valid_diff, diff_tb / step_dt, np.nan)

    return {
        "distance": distance,
        "image_switch": image_switch,
        "near_range": near_range,
        "near_half_box": near_half_box,
        "near_aliasing": near_aliasing,
        "angular_rate_clean": angular_rate_clean,
        "d_min": d_min,
        "image_index": image_index,
    }


def classify_window(
    g: Dict[str, np.ndarray],
    t0: int,
    D: int,
    pre_steps: int = 200,
) -> Tuple[np.ndarray, np.ndarray]:
    """Classify edges for an occlusion window with strict RUN_RULES precedence.

    Parameters
    ----------
    g : dict
        Geometry dictionary returned by geometry().
    t0 : int
        1-based simulation step where occlusion begins (e.g. 1201, 1401, 1501).
    D : int
        Occlusion duration in steps (e.g. 20, 120).
    pre_steps : int
        Number of visible historical steps preceding occlusion (default 200).

    Returns
    -------
    labels : [E] ndarray of str
        Mutually exclusive labels: 'G-clean', 'G-switch', 'G-close', 'G-ambiguous'.
    rates : [E] ndarray of float64
        RMS of valid clean angular rates during the occlusion period (NaN if no valid steps).
    """
    T, E = g["distance"].shape

    # 1-based t0 converts to 0-based index
    occ_start = int(t0) - 1
    occ_end = occ_start + int(D)
    hist_start = occ_start - int(pre_steps)

    labels = np.full(E, LABEL_AMBIGUOUS, dtype=object)
    rates = np.full(E, np.nan, dtype=np.float64)

    # Check window boundary validity
    if hist_start < 0 or occ_end > T:
        return labels, rates

    occ_sl = slice(occ_start, occ_end)
    hist_occ_sl = slice(hist_start, occ_end)

    switches_occ = g["image_switch"][occ_sl]  # [D, E]
    near_occ = g["near_range"][occ_sl]  # [D, E]
    half_box_occ = g["near_half_box"][occ_sl]  # [D, E]

    switches_hist_occ = g["image_switch"][hist_occ_sl]  # [pre_steps + D, E]
    aliasing_hist_occ = g["near_aliasing"][hist_occ_sl]  # [pre_steps + D, E]

    clean_rates_occ = g["angular_rate_clean"][occ_sl]  # [D, E]

    for e in range(E):
        # 1. G-switch: at least one branch switch during occlusion
        if np.any(switches_occ[:, e]):
            labels[e] = LABEL_SWITCH
        # 2. G-close: no switch during occlusion, but distance < 5 during occlusion
        elif np.any(near_occ[:, e]):
            labels[e] = LABEL_CLOSE
        # 3. G-clean:
        #    - no switch across full history (200 steps) + occlusion (D steps)
        #    - no near-aliasing across full history + occlusion
        #    - distance >= 5 during occlusion (already guaranteed by elif above)
        #    - margin to half-box boundary >= 10 during occlusion (no near_half_box in occ)
        elif (
            (not np.any(switches_hist_occ[:, e]))
            and (not np.any(aliasing_hist_occ[:, e]))
            and (not np.any(half_box_occ[:, e]))
        ):
            labels[e] = LABEL_CLEAN
        else:
            labels[e] = LABEL_AMBIGUOUS

        # Compute RMS of valid clean angular rates during occlusion
        rates_e = clean_rates_occ[:, e]
        fin = np.isfinite(rates_e)
        if np.any(fin):
            rates[e] = float(np.sqrt(np.mean(rates_e[fin] ** 2)))
        else:
            rates[e] = np.nan

    return labels, rates


def fresh_labels(g: Dict[str, np.ndarray]) -> np.ndarray:
    """Classify each time step and edge into mutually exclusive categories.

    Strict priority order per step (t, e):
      1. 'switch_event': step has an image_switch
      2. 'switch_near2': within +/- 2 steps of an image_switch on that edge
      3. 'G-close': distance < 5.0 (near_range)
      4. 'G-clean': not near half-box boundary and not near aliasing
      5. 'G-ambiguous': all remaining steps

    Parameters
    ----------
    g : dict
        Geometry dictionary returned by geometry().

    Returns
    -------
    labels : [T, E] ndarray of object (strings)
        Mutually exclusive time-step classification.
    """
    T, E = g["distance"].shape
    labels = np.full((T, E), FRESH_AMBIGUOUS, dtype=object)

    switches = g["image_switch"]  # [T, E] bool
    near_rng = g["near_range"]  # [T, E] bool
    half_box = g["near_half_box"]  # [T, E] bool
    aliasing = g["near_aliasing"]  # [T, E] bool

    for e in range(E):
        sw_e = switches[:, e]
        sw_indices = np.where(sw_e)[0]

        # Construct +/- 2 mask around any switch event
        near2_mask = np.zeros(T, dtype=bool)
        for idx in sw_indices:
            start = max(0, idx - 2)
            end = min(T, idx + 3)
            near2_mask[start:end] = True

        near2_only = near2_mask & (~sw_e)

        # Mutually exclusive priority assignment
        p1 = sw_e
        p2 = near2_only
        p3 = (~near2_mask) & near_rng[:, e]
        p4 = (~near2_mask) & (~near_rng[:, e]) & (~half_box[:, e]) & (~aliasing[:, e])
        p5 = (~p1) & (~p2) & (~p3) & (~p4)

        labels[p1, e] = FRESH_SWITCH_EVENT
        labels[p2, e] = FRESH_SWITCH_NEAR2
        labels[p3, e] = FRESH_CLOSE
        labels[p4, e] = FRESH_CLEAN
        labels[p5, e] = FRESH_AMBIGUOUS

    return labels
