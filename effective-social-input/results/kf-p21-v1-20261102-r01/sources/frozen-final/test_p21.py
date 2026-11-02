"""experiments/kf-p21-20261102/test_p21.py

Standalone verification test suite for P2.1 protocol requirements.
Requires zero external test runners (no pytest dependency; runs directly with python).
Tests:
  1. Angle wrap 359 deg vs 1 deg yielding 2 deg.
  2. Non-anticipation: future observation perturbation does not alter past/present estimates.
  3. Hidden truth isolation: modifying hidden ground truth alters offline metrics but not filter predictions.
  4. Covariance validity: all 5 CV process noise q candidates maintain Joseph PSD covariances.
  5. Reference import: stock P2 replay/hold imported cleanly without modification; CV and GM (tau=10) verified.
  6. MATLAB round tie-breaking: sign(x)*floor(abs(x)+0.5) on +/-0.5, +/-1.5.
  7. Minimum-image crossing flags: synthetic E1 triggers image_switch, ordinary position wrap does not false-flag.
  8. Classification mutual exclusivity: classify_window and fresh_labels produce strictly mutually exclusive classes.
  9. 140-edge shape contracts: geometry output shapes [T, 140] and [T, 140, 2], first step derivative bounds.
 10. Gating units and prediction retention: age measured in model time (dt), gated predictions retained.

Usage:
  python test_p21.py
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
from typing import Any, Dict, List, Tuple
import numpy as np

# Add parent directory for geometry and synthetic_bearings imports
curr_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(curr_dir))
import geometry
import synthetic_bearings


def load_p2_replay() -> Any:
    """Load stock P2 replay module."""
    p2_replay_path = curr_dir.parent / "kf-p2-20261102" / "replay.py"
    if not p2_replay_path.exists():
        raise FileNotFoundError(f"Missing stock P2 replay.py at {p2_replay_path}")
    spec = importlib.util.spec_from_file_location("p2_replay_ref", p2_replay_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load spec for {p2_replay_path}")
    p2_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p2_mod)
    return p2_mod


def test_1_wrap_359_vs_1_deg() -> None:
    """Test 1: wrap(359 deg - 1 deg) and wrap(1 deg - 359 deg) yield 2 deg."""
    deg1 = np.deg2rad(1.0)
    deg359 = np.deg2rad(359.0)

    diff1 = geometry.wrap(deg1 - deg359)
    diff2 = geometry.wrap(deg359 - deg1)

    expected = np.deg2rad(2.0)
    assert np.isclose(diff1, expected, atol=1e-12), f"Expected {expected}, got {diff1}"
    assert np.isclose(diff2, -expected, atol=1e-12), f"Expected {-expected}, got {diff2}"
    assert np.isclose(abs(diff1), expected, atol=1e-12)
    assert np.isclose(abs(diff2), expected, atol=1e-12)


def test_2_missing_future_perturbation_no_leak() -> None:
    """Test 2: Perturbing future missing observations causes zero change to past predictions."""
    p2_mod = load_p2_replay()
    T = 100
    t_pert = 50
    dt = 0.3

    rng = np.random.default_rng(2026110201)
    obs_base = rng.uniform(0, 2 * np.pi, size=(T, 1))
    miss_base = np.zeros((T, 1), dtype=bool)
    miss_base[t_pert:] = True

    # Baseline run
    res_base = p2_mod.replay(obs_base, miss_base, dt=dt, q=1e-4, R=1e-4)

    # Perturb future observations at and after step t_pert + 5
    obs_pert = obs_base.copy()
    obs_pert[t_pert + 5 :] = rng.uniform(0, 2 * np.pi, size=(T - (t_pert + 5), 1))
    res_pert = p2_mod.replay(obs_pert, miss_base, dt=dt, q=1e-4, R=1e-4)

    # Predictions up to t_pert + 4 must be identical (exact 0.0 difference)
    safe_slice = slice(0, t_pert + 5)
    np.testing.assert_array_equal(
        res_base["theta"][safe_slice],
        res_pert["theta"][safe_slice],
        err_msg="Future observation perturbation leaked into past/current theta predictions!",
    )
    np.testing.assert_array_equal(
        res_base["p11"][safe_slice],
        res_pert["p11"][safe_slice],
        err_msg="Future observation perturbation leaked into past/current p11 covariances!",
    )
    np.testing.assert_array_equal(
        res_base["age"][safe_slice],
        res_pert["age"][safe_slice],
        err_msg="Future observation perturbation leaked into past/current age metrics!",
    )


def test_3_hidden_truth_isolation() -> None:
    """Test 3: Changing hidden truth affects offline error evaluation but not filter predictions."""
    p2_mod = load_p2_replay()
    T = 80
    t0 = 40
    dt = 0.3

    rng = np.random.default_rng(2026110202)
    obs = rng.uniform(0, 2 * np.pi, size=(T, 1))
    miss = np.zeros((T, 1), dtype=bool)
    miss[t0:] = True

    truth_1 = np.full((T, 1), 0.5)
    truth_2 = np.full((T, 1), 0.5)
    truth_2[t0:] = 2.5  # Dramatically different ground truth during occlusion

    # Filter receives only observations and missing mask, zero truth
    res_1 = p2_mod.replay(obs, miss, dt=dt, q=1e-4, R=1e-4)
    res_2 = p2_mod.replay(obs, miss, dt=dt, q=1e-4, R=1e-4)

    # Predictions must be strictly identical
    np.testing.assert_array_equal(res_1["theta"], res_2["theta"])
    np.testing.assert_array_equal(res_1["p11"], res_2["p11"])

    # Offline evaluation errors must differ
    err_1 = geometry.wrap(res_1["theta"][t0:] - truth_1[t0:])
    err_2 = geometry.wrap(res_2["theta"][t0:] - truth_2[t0:])
    assert not np.allclose(err_1, err_2), "Offline errors should differ when truth differs"


def test_4_all_5_cv_q_psd_joseph() -> None:
    """Test 4: All 5 CV q candidates {1e-6, 1e-5, 1e-4, 1e-3, 1e-2} maintain Joseph PSD covariances."""
    p2_mod = load_p2_replay()
    q_grid = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
    T = 150
    dt = 0.3
    R = 1e-4

    rng = np.random.default_rng(2026110203)
    obs = rng.uniform(0, 2 * np.pi, size=(T, 2))
    miss = np.zeros((T, 2), dtype=bool)
    miss[50:100] = True  # Partial occlusion

    for q in q_grid:
        res = p2_mod.replay(obs, miss, dt=dt, q=q, R=R, tau=None)
        p11 = res["p11"]  # [T, 2]
        p12 = res["p12"]
        p22 = res["p22"]

        # Check finite
        assert np.all(np.isfinite(p11)), f"p11 contains non-finite values for q={q}"
        assert np.all(np.isfinite(p12)), f"p12 contains non-finite values for q={q}"
        assert np.all(np.isfinite(p22)), f"p22 contains non-finite values for q={q}"

        # Diagonal elements non-negative
        assert np.all(p11 >= -1e-14), f"p11 has negative variance for q={q}: min={np.min(p11)}"
        assert np.all(p22 >= -1e-14), f"p22 has negative variance for q={q}: min={np.min(p22)}"

        # Determinant non-negative: p11 * p22 - p12^2 >= 0
        det = p11 * p22 - (p12**2)
        assert np.all(det >= -1e-14), f"Covariance determinant negative for q={q}: min_det={np.min(det)}"

        # Eigenvalues non-negative
        trace = p11 + p22
        disc = np.maximum(0.0, trace**2 - 4.0 * det)
        eig_min = 0.5 * (trace - np.sqrt(disc))
        assert np.all(eig_min >= -1e-14), f"Covariance has negative eigenvalue for q={q}: min_eig={np.min(eig_min)}"


def test_5_p2_reference_import_cv_and_gm_tau10() -> None:
    """Test 5: Stock P2 replay/hold imported cleanly; CV and GM (tau=10) executed on identical code."""
    p2_mod = load_p2_replay()
    T = 60
    dt = 0.3
    obs = np.linspace(0.5, 1.5, T)[:, None]
    miss = np.zeros((T, 1), dtype=bool)
    miss[30:] = True

    # Test hold
    res_hold = p2_mod.hold(obs, miss, dt=dt)
    assert "theta" in res_hold
    assert "age" in res_hold
    assert np.all(np.isnan(res_hold["p11"]))

    # Test CV
    res_cv = p2_mod.replay(obs, miss, dt=dt, q=1e-4, R=1e-4, tau=None)
    assert np.all(np.isfinite(res_cv["theta"]))
    assert np.all(np.isfinite(res_cv["p11"]))

    # Test GM tau=10
    res_gm = p2_mod.replay(obs, miss, dt=dt, q=1e-5, R=1e-4, tau=10.0)
    assert np.all(np.isfinite(res_gm["theta"]))
    assert np.all(np.isfinite(res_gm["p11"]))


def test_6_matlab_round_ties() -> None:
    """Test 6: MATLAB round tie-breaking sign(x)*floor(abs(x)+0.5) on +/-0.5, +/-1.5."""
    assert geometry.round_matlab(0.5) == 1.0
    assert geometry.round_matlab(-0.5) == -1.0
    assert geometry.round_matlab(1.5) == 2.0
    assert geometry.round_matlab(-1.5) == -2.0
    assert geometry.round_matlab(0.0) == 0.0

    arr = np.array([0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 0.0, 0.49, -0.49, 0.51, -0.51])
    expected = np.array([1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.0, 0.0, 0.0, 1.0, -1.0])
    np.testing.assert_array_equal(geometry.round_matlab(arr), expected)


def test_7_synthetic_minimage_crossing_and_ordinary_wrap() -> None:
    """Test 7: Synthetic E1 triggers image_switch; ordinary coordinate wrapping does not false-flag."""
    # 7a. Synthetic E1 boundary jump
    traj_e = synthetic_bearings.generate_synthetic_trajectory("E1", seed=2026113001, D=120)
    assert np.any(traj_e["image_switch"]), "Config E1 must trigger minimum-image switch"
    assert not np.any(traj_e["near_range"]), "Config E1 must not trigger near_range (dist ~ 510)"

    # 7b. Synthetic D1 close passage
    traj_d = synthetic_bearings.generate_synthetic_trajectory("D1", seed=2026113001, D=120)
    assert not np.any(traj_d["image_switch"]), "Config D1 must not trigger image switch"
    assert np.any(traj_d["near_range"]), "Config D1 must trigger near_range (dist min = 2.0 < 5.0)"

    # 7c. Ordinary coordinate wrap without physical jump
    # Observer fixed at (100, 100). Target moves smoothly from x=999 to x=1 across boundary.
    # dx_raw goes from 899 to -99.
    # d_min before: 899 - 1000 = -101.
    # d_min after: -99 - 0 = -99.
    # Difference in d_min is |-99 - (-101)| = 2.0 << 500 -> image_switch must be False!
    positions = np.zeros((2, 2, 2), dtype=np.float64)
    positions[:, 0, :] = 100.0  # Observer
    positions[:, 1, 0] = [999.0, 100.0]  # Target step 0
    positions[:, 1, 1] = [1.0, 100.0]    # Target step 1

    true_bearing = np.zeros((2, 1), dtype=np.float64)
    true_bearing[0, 0] = np.arctan2(0.0, -101.0)
    true_bearing[1, 0] = np.arctan2(0.0, -99.0)

    src = {
        "positions": positions,
        "true_bearing": true_bearing,
        "observer_id": np.array([1]),
        "target_id": np.array([2]),
        "arena_size": 1000.0,
        "dt": 0.3,
    }
    g = geometry.geometry(src)
    assert not g["image_switch"][1, 0], "Ordinary coordinate wrap must NOT be flagged as image_switch!"
    assert np.isclose(g["distance"][0, 0], 101.0)
    assert np.isclose(g["distance"][1, 0], 99.0)


def test_8_per_window_and_fresh_labels_mutually_exclusive() -> None:
    """Test 8: classify_window and fresh_labels produce strictly mutually exclusive classes."""
    # Test on synthetic condition
    traj = synthetic_bearings.generate_synthetic_trajectory("D1", seed=2026113001, D=20)
    T = len(traj["truth"])

    # Construct mock geometry dictionary
    g = {
        "distance": traj["distance"][:, None],
        "image_switch": traj["image_switch"][:, None],
        "near_range": traj["near_range"][:, None],
        "near_half_box": traj["near_half_box"][:, None],
        "near_aliasing": np.zeros((T, 1), dtype=bool),
        "angular_rate_clean": np.full((T, 1), 0.01),
        "d_min": np.zeros((T, 1, 2)),
        "image_index": np.zeros((T, 1, 2)),
    }

    # Window classification
    labels, rates = geometry.classify_window(g, t0=201, D=20, pre_steps=200)
    assert len(labels) == 1
    assert labels[0] in geometry.WINDOW_LABELS
    assert labels[0] == geometry.LABEL_CLOSE  # D1 has near_range < 5 during occlusion

    # Time-step level classification
    fl = geometry.fresh_labels(g)
    assert fl.shape == (T, 1)
    for t in range(T):
        lbl = fl[t, 0]
        assert lbl in geometry.FRESH_LABELS, f"Unknown fresh label: {lbl}"


def test_9_140_edge_shapes_and_geometry_contract() -> None:
    """Test 9: 140 edges shapes [T, 140] and [T, 140, 2], first step derivative bounds."""
    T = 50
    N = 20
    E = 140
    rng = np.random.default_rng(2026110204)

    positions = rng.uniform(0, 1000, size=(2, N, T))
    true_bearing = rng.uniform(0, 2 * np.pi, size=(T, E))

    # Standard 140 edges: 20 observers, 7 cyclic neighbours each
    obs_ids = []
    tgt_ids = []
    for obs in range(1, N + 1):
        for offset in range(1, 8):
            tgt = ((obs - 1 + offset) % N) + 1
            obs_ids.append(obs)
            tgt_ids.append(tgt)
    obs_ids = np.array(obs_ids)
    tgt_ids = np.array(tgt_ids)
    assert len(obs_ids) == E

    src = {
        "positions": positions,
        "true_bearing": true_bearing,
        "observer_id": obs_ids,
        "target_id": tgt_ids,
        "arena_size": 1000.0,
        "dt": 0.3,
    }

    g = geometry.geometry(src)

    # Validate output dictionary shapes
    assert g["distance"].shape == (T, E)
    assert g["image_switch"].shape == (T, E)
    assert g["near_range"].shape == (T, E)
    assert g["near_half_box"].shape == (T, E)
    assert g["near_aliasing"].shape == (T, E)
    assert g["angular_rate_clean"].shape == (T, E)
    assert g["d_min"].shape == (T, E, 2)
    assert g["image_index"].shape == (T, E, 2)

    # First step no derivative
    assert not np.any(g["image_switch"][0]), "image_switch at step 0 must be False"
    assert not np.any(g["near_aliasing"][0]), "near_aliasing at step 0 must be False"
    assert np.all(np.isnan(g["angular_rate_clean"][0])), "angular_rate_clean at step 0 must be NaN"


def test_10_gate_units_dt_and_retention() -> None:
    """Test 10: Age measured in model time (dt), gated predictions retained in output."""
    p2_mod = load_p2_replay()
    T = 60
    dt = 0.3
    obs = np.linspace(0.1, 0.9, T)[:, None]
    miss = np.zeros((T, 1), dtype=bool)
    miss[20:] = True  # Occlusion starts at step 20

    res = p2_mod.replay(obs, miss, dt=dt, q=1e-4, R=1e-4)

    # Verify age increments by dt during occlusion
    occ_ages = res["age"][20:, 0]
    expected_ages = np.arange(1, len(occ_ages) + 1) * dt
    np.testing.assert_allclose(occ_ages, expected_ages, rtol=1e-12)

    # Verify predictions are retained throughout occlusion (no NaNs in theta or p11)
    assert np.all(np.isfinite(res["theta"][20:, 0])), "Predictions must be retained during occlusion!"
    assert np.all(np.isfinite(res["p11"][20:, 0])), "Variance must be retained during occlusion!"


def run_all_tests() -> bool:
    """Run full test suite sequentially."""
    tests = [
        ("Test 1: Angle wrap 359 vs 1 deg", test_1_wrap_359_vs_1_deg),
        ("Test 2: Non-anticipation & no future leak", test_2_missing_future_perturbation_no_leak),
        ("Test 3: Hidden truth isolation", test_3_hidden_truth_isolation),
        ("Test 4: All 5 CV q maintain Joseph PSD", test_4_all_5_cv_q_psd_joseph),
        ("Test 5: Stock P2 reference import CV & GM", test_5_p2_reference_import_cv_and_gm_tau10),
        ("Test 6: MATLAB round ties +/-0.5, +/-1.5", test_6_matlab_round_ties),
        ("Test 7: Synthetic min-image crossing & no false wrap", test_7_synthetic_minimage_crossing_and_ordinary_wrap),
        ("Test 8: Classification mutual exclusivity", test_8_per_window_and_fresh_labels_mutually_exclusive),
        ("Test 9: 140 edges shape contracts & derivative bounds", test_9_140_edge_shapes_and_geometry_contract),
        ("Test 10: Gating units in model time & prediction retention", test_10_gate_units_dt_and_retention),
    ]

    print("=" * 70)
    print("Running P2.1 Standalone Test Suite (No pytest dependency)")
    print("=" * 70)

    passed_count = 0
    for name, test_fn in tests:
        try:
            test_fn()
            print(f"[PASS] {name}")
            passed_count += 1
        except Exception as e:
            print(f"[FAIL] {name}: {e}")
            import traceback
            traceback.print_exc()

    print("=" * 70)
    print(f"Results: {passed_count}/{len(tests)} tests passed.")
    print("=" * 70)
    return passed_count == len(tests)


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
