"""Codex independent checks; does not import the experiment's analysis code."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[3]
REVIEW = Path(__file__).resolve().parent


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args])


def hashes():
    paths = git("ls-files", "-z").decode("utf-8").split("\0")
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
            for p in paths if p}


def covariance_prediction():
    dt, q, r = 0.3, 1e-4, 1e-4
    f = np.array([[1.0, dt], [0.0, 1.0]])
    qmat = q * np.array([[dt**3 / 3, dt**2 / 2], [dt**2 / 2, dt]])
    p = np.diag([r, 1.0])
    for _ in range(2, 1401):
        p = f @ p @ f.T + qmat
        k = p[:, 0] / (p[0, 0] + r)
        a = np.eye(2) - np.outer(k, [1.0, 0.0])
        p = a @ p @ a.T + np.outer(k, k) * r
        p = (p + p.T) / 2
    last_fresh = p.tolist()
    missing_variance = []
    for _ in range(120):
        p = f @ p @ f.T + qmat
        p = (p + p.T) / 2
        missing_variance.append(float(p[0, 0]))
    truncated = [i + 1 for i, v in enumerate(missing_variance) if v > 1.0]
    return {
        "method": "Independent numpy Joseph/Riccati covariance recurrence",
        "last_fresh_P": last_fresh,
        "first_P_cutoff_missing_step": min(truncated),
        "D120_P_cutoff_steps": len(truncated),
        "D120_hidden_edge_steps": 20 * 7 * 120,
        "D120_expected_P_cutoff_edge_steps": 20 * 7 * len(truncated),
        "D120_age_time": 120 * dt,
        "age_threshold_time": 120,
        "D120_expected_age_cutoff_edge_steps": 0,
        "P_by_missing_step": missing_variance,
        "status": "analytic software expectation; MATLAB comparison pending",
    }


def verify_geometry(path, loaded=None):
    """Reconstruct geometry and velocity metrics without importing AGY code."""
    content = loadmat(path, simplify_cells=True) if loaded is None else loaded
    result = content["result"]
    cfg = result["config"]
    traj = result["trajectory"]
    series = result["timeSeries"]
    position = np.stack([traj["x"], traj["y"]], axis=-1).transpose(1, 0, 2)
    steps, n = int(cfg["NSteps"]), int(cfg["NAgents"])
    length, dt = float(cfg["arenaSize"]), float(cfg["dt"])
    assert position.shape == (steps + 1, n, 2)
    assert int(cfg["recordStride"]) == 1
    displacement = np.diff(position, axis=0)
    displacement -= length * np.round(displacement / length)
    speed = np.linalg.norm(displacement, axis=2) / dt
    unit = np.divide(displacement, dt * speed[..., None],
                     out=np.zeros_like(displacement), where=speed[..., None] > 1e-12)
    global_order = np.linalg.norm(unit.sum(axis=1), axis=1) / n
    difference = position[1:, :, None, :] - position[1:, None, :, :]
    difference -= length * np.round(difference / length)
    distance = np.linalg.norm(difference, axis=-1)
    i, j = np.triu_indices(n, 1)
    pair = distance[:, i, j].mean(axis=1) / length
    nearest_distance = distance.copy()
    nearest_distance[:, np.arange(n), np.arange(n)] = np.inf
    nearest = nearest_distance.min(axis=2).mean(axis=1) / length
    order = np.argsort(distance, axis=2, kind="stable")[:, :, :min(6, n)]
    neighbors = unit[np.arange(steps)[:, None, None], order]
    local_order = np.linalg.norm(neighbors.sum(axis=2), axis=2).mean(axis=1) / min(6, n)
    reconstructed = {
        "globalOrder": global_order, "localOrder": local_order,
        "meanPairDistanceNorm": pair, "meanNearestDistanceNorm": nearest,
        "meanSpeed": speed.mean(axis=1),
    }
    errors = {key: float(np.max(np.abs(value - np.asarray(series[key]))))
              for key, value in reconstructed.items()}
    passed = all(value < 1e-9 for value in errors.values())
    return {"file": str(Path(path).resolve()), "N": n, "steps": steps,
            "seed": int(cfg["seed"]), "mode": cfg["predictionMode"],
            "max_absolute_errors": errors, "pass": passed}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["baseline", "verify", "covariance", "geometry"])
    parser.add_argument("mat_files", nargs="*")
    args = parser.parse_args()
    baseline_path = REVIEW / "protected_baseline.json"
    if args.mode == "baseline":
        data = {"commit": git("rev-parse", "HEAD").decode().strip(),
                "branch": git("branch", "--show-current").decode().strip(),
                "tracked_sha256": hashes()}
        with baseline_path.open("x", encoding="utf-8") as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
        print(json.dumps({"tracked_files": len(data["tracked_sha256"]),
                          "commit": data["commit"]}))
    elif args.mode == "verify":
        baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
        now = hashes()
        changed = [p for p, h in baseline["tracked_sha256"].items() if now.get(p) != h]
        print(json.dumps({"tracked_files_checked": len(baseline["tracked_sha256"]),
                          "changed": changed, "pass": not changed}, ensure_ascii=False))
        if changed:
            raise SystemExit(1)
    elif args.mode == "covariance":
        data = covariance_prediction()
        path = REVIEW / "independent_covariance_expected.json"
        with path.open("x", encoding="utf-8") as stream:
            json.dump(data, stream, indent=2)
        print(json.dumps({k: v for k, v in data.items() if k != "P_by_missing_step"}))
    else:
        assert args.mat_files, "Supply one or more MAT files."
        records = [verify_geometry(path) for path in args.mat_files]
        print(json.dumps(records, ensure_ascii=False, indent=2))
        if not all(record["pass"] for record in records):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
