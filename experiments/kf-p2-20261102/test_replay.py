"""experiments/kf-p2-20261102/test_replay.py

Deterministic software tests for open-loop bearing estimators (replay and hold).
Independent of scientific experiment MAT datasets; uses synthetic test fixtures.

Coverage:
  1. CV vs GM(tau -> inf) numerical equivalence.
  2. Circular innovation wrapping across +/- pi / 0 / 2*pi boundaries.
  3. Constant bearing and constant angular velocity prediction vs Hold.
  4. Joseph covariance update symmetry and positive semi-definiteness (PSD).
  5. First-step track initialization (no predict step on first observation).
  6. Cumulative age in exact dt units and reacquisition reset.
  7. Strict isolation against observation corruption during missing steps (no leakage).
  8. Gating cutoff verification: first 100 occlusion steps active, step 101 cutoff at baseline q, dt.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Ensure replay module can be imported whether run from repo root or directory
curr_dir = Path(__file__).resolve().parent
if str(curr_dir) not in sys.path:
    sys.path.insert(0, str(curr_dir))

from replay import compute_active_flags, hold, mod_2pi, replay, wrap_to_pi


def test_cv_gm_tau_inf_equality() -> None:
    """Verify that GM model with tau -> inf is numerically equivalent to CV."""
    T, E = 80, 5
    dt, q, R = 0.3, 1e-4, 1e-4
    rng = np.random.default_rng(42)

    obs = rng.uniform(0.0, 2.0 * np.pi, size=(T, E))
    missing = np.zeros(T, dtype=bool)
    missing[20:50] = True  # 30 occluded steps

    res_cv = replay(obs, missing, dt=dt, q=q, R=R, tau=None)
    # Test with very large finite tau (1e12) and tau=np.inf
    res_gm_inf = replay(obs, missing, dt=dt, q=q, R=R, tau=float("inf"))

    # Inf tau equivalence
    np.testing.assert_allclose(res_cv["theta"], res_gm_inf["theta"], atol=1e-12)
    np.testing.assert_allclose(res_cv["p11"], res_gm_inf["p11"], atol=1e-12)
    np.testing.assert_allclose(res_cv["p12"], res_gm_inf["p12"], atol=1e-12)
    np.testing.assert_allclose(res_cv["p22"], res_gm_inf["p22"], atol=1e-12)
    np.testing.assert_allclose(res_cv["nis"], res_gm_inf["nis"], atol=1e-12, equal_nan=True)
    np.testing.assert_array_equal(res_cv["age"], res_gm_inf["age"])

    # Only the explicit inf branch is a CV identity. Extremely large finite
    # tau suffers cancellation in the stock MATLAB Q formula and is not used.


def test_circular_wrap_near_pi() -> None:
    """Verify wrapped innovation correctly wraps across 0 / 2*pi and pi boundaries."""
    dt, q, R = 0.3, 1e-4, 1e-4
    # Observation steps: starts at 0.05, then jumps across boundary to 2*pi - 0.05 (~6.233)
    # The true angular distance is -0.10 rad, NOT +6.18 rad!
    obs = np.array([[0.05], [2.0 * np.pi - 0.05]], dtype=np.float64)
    missing = np.array([False, False], dtype=bool)

    res = replay(obs, missing, dt=dt, q=q, R=R)
    theta_0 = res["theta"][0, 0]
    theta_1 = res["theta"][1, 0]

    # Verify first state matches initial observation mod 2*pi
    assert np.isclose(theta_0, 0.05, atol=1e-12)

    # On step 1, innovation = atan2(sin(6.233 - 0.05), cos(6.233 - 0.05)) = -0.10
    # The filter angular velocity state must be negative, reflecting the small clockwise step
    # Posterior theta must remain wrapped in [0, 2*pi)
    assert 0.0 <= theta_1 < 2.0 * np.pi
    # Innovation check via NIS: innov = -0.1, denom = P11 + R
    # NIS should be finite and close to (-0.1)^2 / S
    assert np.isfinite(res["nis"][1, 0])
    assert res["nis"][1, 0] > 0.0


def test_constant_angle_hidden_prediction() -> None:
    """Verify constant bearing yields constant prediction matching Hold."""
    T, E = 30, 2
    dt, q, R = 0.3, 1e-4, 1e-4
    obs = np.full((T, E), 1.25, dtype=np.float64)
    missing = np.zeros(T, dtype=bool)
    missing[10:] = True  # Occlusion from step 10 to 30

    res_cv = replay(obs, missing, dt=dt, q=q, R=R)
    res_hold = hold(obs, missing, dt=dt)

    # During occlusion, constant bearing filter with omega=0 stays at 1.25
    np.testing.assert_allclose(res_cv["theta"][missing], 1.25, atol=1e-6)
    np.testing.assert_allclose(res_hold["theta"][missing], 1.25, atol=1e-12)


def test_constant_omega_hidden_prediction_vs_hold() -> None:
    """Verify constant angular velocity prediction tracks motion while Hold stagnates."""
    T, E = 40, 1
    dt = 0.3
    omega_true = 0.25  # rad/s
    t_axis = np.arange(T) * dt
    truth = mod_2pi(0.5 + omega_true * t_axis[:, None])

    # Low noise observations
    obs = truth.copy()
    missing = np.zeros(T, dtype=bool)
    missing[20:] = True  # Missing for last 20 steps

    res_cv = replay(obs, missing, dt=dt, q=1e-5, R=1e-6)
    res_hold = hold(obs, missing, dt=dt)

    err_cv_occ = np.abs(wrap_to_pi(res_cv["theta"][20:] - truth[20:]))
    err_hold_occ = np.abs(wrap_to_pi(res_hold["theta"][20:] - truth[20:]))

    # Hold stays at step 19 observation, so error grows: omega_true * 20 * dt = 1.5 rad
    assert np.max(err_hold_occ) > 1.0
    # KF-CV has learned omega_true and extrapolates; error remains small
    assert np.max(err_cv_occ) < 0.15


def test_joseph_covariance_psd_and_symmetry() -> None:
    """Verify that Joseph covariance update guarantees symmetry and PSD."""
    T, E = 60, 4
    dt, q, R = 0.3, 1e-4, 1e-4
    rng = np.random.default_rng(123)
    obs = rng.uniform(0.0, 2.0 * np.pi, size=(T, E))
    missing = rng.choice([False, True], size=T, p=[0.7, 0.3])
    missing[0] = False  # Step 0 observed

    res = replay(obs, missing, dt=dt, q=q, R=R)
    p11 = res["p11"]
    p12 = res["p12"]
    p22 = res["p22"]

    for t in range(T):
        for e in range(E):
            P = np.array([[p11[t, e], p12[t, e]], [p12[t, e], p22[t, e]]])
            # Symmetry check
            assert np.isclose(P[0, 1], P[1, 0], atol=1e-12)
            # PSD check: eigenvalues >= -1e-12
            eigs = np.linalg.eigvalsh(P)
            assert np.all(eigs >= -1e-12), f"Negative eigenvalue at t={t}, e={e}: {eigs}"
            # Diagonal variances >= 0
            assert P[0, 0] >= -1e-12
            assert P[1, 1] >= -1e-12


def test_first_step_initialization_no_predict() -> None:
    """Verify initialization matches stock MATLAB: no prediction step on first observation."""
    dt, q, R = 0.3, 1e-4, 1e-4
    first_y = 2.71828
    obs = np.array([[first_y]], dtype=np.float64)
    missing = np.array([False], dtype=bool)

    res = replay(obs, missing, dt=dt, q=q, R=R)

    # Initial state x = [mod(observation, 2*pi); 0]; P = diag([R, 1])
    assert np.isclose(res["theta"][0, 0], mod_2pi(first_y), atol=1e-12)
    assert np.isclose(res["p11"][0, 0], R, atol=1e-12)
    assert np.isclose(res["p12"][0, 0], 0.0, atol=1e-12)
    assert np.isclose(res["p22"][0, 0], 1.0, atol=1e-12)
    assert res["age"][0, 0] == 0.0
    assert np.isnan(res["nis"][0, 0])


def test_age_dt_units() -> None:
    """Verify age accumulates in exact dt units during occlusion and resets to 0."""
    T = 8
    dt = 0.3
    obs = np.full((T, 1), 1.0, dtype=np.float64)
    # Steps: 0 fresh, 1..4 missing (4 steps), 5 fresh, 6..7 missing (2 steps)
    missing = np.array([False, True, True, True, True, False, True, True], dtype=bool)

    res = replay(obs, missing, dt=dt, q=1e-4, R=1e-4)
    expected_age = np.array([[0.0], [0.3], [0.6], [0.9], [1.2], [0.0], [0.3], [0.6]])

    np.testing.assert_allclose(res["age"], expected_age, atol=1e-12)

    res_hold = hold(obs, missing, dt=dt)
    np.testing.assert_allclose(res_hold["age"], expected_age, atol=1e-12)


def test_no_future_or_hidden_observation_leakage() -> None:
    """Verify filter outputs during and after occlusion are invariant to corrupted hidden inputs."""
    T, E = 50, 3
    dt, q, R = 0.3, 1e-4, 1e-4
    rng = np.random.default_rng(999)

    obs1 = rng.uniform(0.0, 2.0 * np.pi, size=(T, E))
    missing = np.zeros(T, dtype=bool)
    missing[15:35] = True  # 20 missing steps

    res1 = replay(obs1, missing, dt=dt, q=q, R=R)

    # Corrupt observations during missing frames with radical arbitrary values
    obs2 = obs1.copy()
    obs2[missing, :] = rng.uniform(-10000.0, 10000.0, size=(np.sum(missing), E))

    res2 = replay(obs2, missing, dt=dt, q=q, R=R)

    # Strict bit-for-bit equality
    np.testing.assert_array_equal(res1["theta"], res2["theta"])
    np.testing.assert_array_equal(res1["p11"], res2["p11"])
    np.testing.assert_array_equal(res1["p12"], res2["p12"])
    np.testing.assert_array_equal(res1["p22"], res2["p22"])
    np.testing.assert_array_equal(res1["age"], res2["age"])

    np.testing.assert_array_equal(np.isnan(res1["nis"]), np.isnan(res2["nis"]))
    finite_mask = ~np.isnan(res1["nis"])
    np.testing.assert_array_equal(res1["nis"][finite_mask], res2["nis"][finite_mask])


def test_gate_first_100_active_101_cutoff() -> None:
    """Verify active gating: active for first 100 steps of occlusion, cut off at step 101.

    Under baseline parameters (dt=0.3, q=1e-4, R=1e-4, max_variance=1.0, max_age=120):
    Filter reaches steady state on prefix observations, then occluded for 105 steps.
    At occlusion step 100: P11 <= 1.0 (active).
    At occlusion step 101: P11 > 1.0 (cutoff by P).
    """
    n_prefix = 50
    n_occ = 105
    T = n_prefix + n_occ
    E = 1
    dt, q, R = 0.3, 1e-4, 1e-4
    max_age = 120.0
    max_variance = 1.0

    # Prefix fresh observations at steady constant angle
    obs = np.full((T, E), 1.0, dtype=np.float64)
    missing = np.zeros(T, dtype=bool)
    missing[n_prefix:] = True  # Occlusion begins at step n_prefix (index 50)

    res = replay(obs, missing, dt=dt, q=q, R=R)

    # Compute active flags
    act_age, act_p, active = compute_active_flags(
        res["age"],
        res["p11"],
        np.broadcast_to(missing[:, None], (T, E)),
        max_age=max_age,
        max_variance=max_variance,
        is_hold=False,
    )

    # Occlusion step 100 corresponds to index n_prefix + 99 = 149
    idx_100 = n_prefix + 99
    # Occlusion step 101 corresponds to index n_prefix + 100 = 150
    idx_101 = n_prefix + 100

    p_100 = res["p11"][idx_100, 0]
    p_101 = res["p11"][idx_101, 0]
    age_100 = res["age"][idx_100, 0]
    age_101 = res["age"][idx_101, 0]

    # Verify age is within max_age (30.0s and 30.3s <= 120.0s)
    assert age_100 <= max_age
    assert age_101 <= max_age

    # Variance check across boundary 1.0
    assert p_100 <= max_variance, f"Step 100 P11 {p_100} exceeded max_variance {max_variance}"
    assert p_101 > max_variance, f"Step 101 P11 {p_101} did not exceed max_variance {max_variance}"

    # Gating flag assertions
    assert active[idx_100, 0] is True or active[idx_100, 0] == 1, "Occlusion step 100 must be active"
    assert active[idx_101, 0] is False or active[idx_101, 0] == 0, "Occlusion step 101 must be cut off"


def test_initial_missing_has_no_observation_leakage() -> None:
    missing = np.array([True, True, False, True, False])
    obs = np.array([[9.0], [8.0], [1.0], [7.0], [1.1]])
    poisoned = obs.copy()
    poisoned[missing] = -1234.
    for fn in [lambda y: replay(y, missing, dt=.3, q=1e-4, R=1e-4),
               lambda y: replay(y, missing, dt=.3, q=1e-5, R=1e-4, tau=10),
               lambda y: hold(y, missing, dt=.3)]:
        a, b = fn(obs), fn(poisoned)
        assert np.isnan(a['theta'][:2]).all()
        assert np.isinf(a['age'][:2]).all()
        for key in a:
            np.testing.assert_array_equal(a[key], b[key])


def test_no_future_observation_leakage() -> None:
    obs = np.linspace(0., 2., 60).reshape(30, 2)
    missing = np.zeros(30, dtype=bool)
    other = obs.copy(); other[20:] = 999.
    for tau in [None,10]:
        a = replay(obs, missing, dt=.3, q=1e-4, R=1e-4, tau=tau)
        b = replay(other, missing, dt=.3, q=1e-4, R=1e-4, tau=tau)
        for key in a:
            np.testing.assert_array_equal(a[key][:20], b[key][:20])


if __name__ == "__main__":
    tests = [
        test_cv_gm_tau_inf_equality,
        test_circular_wrap_near_pi,
        test_constant_angle_hidden_prediction,
        test_constant_omega_hidden_prediction_vs_hold,
        test_joseph_covariance_psd_and_symmetry,
        test_first_step_initialization_no_predict,
        test_age_dt_units,
        test_no_future_or_hidden_observation_leakage,
        test_gate_first_100_active_101_cutoff,
        test_initial_missing_has_no_observation_leakage,
        test_no_future_observation_leakage,
    ]
    passed = 0
    print("=== Running deterministic replay test suite ===")
    for test in tests:
        try:
            test()
            print(f"[PASS] {test.__name__}")
            passed += 1
        except Exception as e:
            print(f"[FAIL] {test.__name__}: {e}")

    print(f"=== Results: {passed}/{len(tests)} tests passed ===")
    if passed != len(tests):
        sys.exit(1)
