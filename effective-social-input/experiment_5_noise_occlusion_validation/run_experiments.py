"""Run the small paired noise-by-occlusion validation experiment."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from threadpoolctl import threadpool_limits

import experiment_model as model

ROOT = Path(__file__).resolve().parent
SEEDS = list(range(2026105001, 2026105021))
TOTAL_STEPS = model.TOTAL_STEPS


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tasks():
    rows = []
    for seed in SEEDS:
        leaders, _ = model.role_ids(seed)
        for noise_std in model.NOISE_LEVELS:
            noise_label = f"noise{noise_std:.3f}".replace(".", "p")
            for mode in model.METHODS:
                run_id = f"N{model.N}_L{model.LEADER_COUNT}_tracking_leader_loss_{noise_label}_{mode}_seed{seed}"
                rows.append(dict(
                    run_id=run_id, seed=seed, n=model.N, leaders=len(leaders),
                    followers=model.N-len(leaders), leader_ids_0based=json.dumps(leaders.tolist()),
                    noise_std=noise_std, R=noise_std ** 2, condition="tracking",
                    tracking=1, mask_type="leader_loss", mode=mode,
                    data=f"runs/{run_id}.csv",
                ))
            run_id = f"N{model.N}_L{model.LEADER_COUNT}_tracking_reference_{noise_label}_seed{seed}"
            rows.append(dict(
                run_id=run_id, seed=seed, n=model.N, leaders=len(leaders),
                followers=model.N-len(leaders), leader_ids_0based=json.dumps(leaders.tolist()),
                noise_std=noise_std, R=noise_std ** 2, condition="tracking", tracking=1,
                mask_type="reference", mode="fresh", data=f"runs/{run_id}.csv",
            ))
    assert len(rows) == len(SEEDS) * len(model.NOISE_LEVELS) * (len(model.METHODS) + 1)
    return rows


def write_json(path: Path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def preflight():
    seed = SEEDS[0]
    cfg_low = model.config(noise_std=model.NOISE_LEVELS[0], total_steps=500)
    cfg_mid = model.config(noise_std=model.NOISE_LEVELS[1], total_steps=500)
    initial_low = model.initialize(seed, cfg_low)
    initial_mid = model.initialize(seed, cfg_mid)
    assert np.allclose(initial_low[0], initial_mid[0])
    assert np.allclose(initial_low[1], initial_mid[1])
    assert np.allclose(initial_low[2], initial_mid[2])
    assert np.allclose(initial_mid[4] / model.NOISE_LEVELS[1], initial_low[4] / model.NOISE_LEVELS[0])
    selected, hidden, leaders, followers = model.topology(cfg_low, seed, "leader_loss")
    assert selected.sum() == model.N * model.K_NEIGHBORS
    assert hidden.sum() == model.HIDDEN_PER_FOLLOWER * len(followers)
    assert np.all(hidden.sum(axis=1)[followers] == model.HIDDEN_PER_FOLLOWER)
    assert model.config(model.NOISE_LEVELS[0]).R == model.NOISE_LEVELS[0] ** 2
    assert model.config(model.NOISE_LEVELS[1]).R == model.NOISE_LEVELS[1] ** 2
    prefix_metrics = ["GO_all", "GO_followers", "pair_distance"]
    for noise_std, initial, cfg in ((model.NOISE_LEVELS[0], initial_low, cfg_low), (model.NOISE_LEVELS[1], initial_mid, cfg_mid)):
        hold = model.simulate(seed, "hold", True, "leader_loss", initial=initial, cfg=cfg)
        kf = model.simulate(seed, "kf-mean", True, "leader_loss", initial=initial, cfg=cfg)
        fresh = model.simulate(seed, "fresh", True, "leader_loss", initial=initial, cfg=cfg)
        a = np.array([[row[key] for key in prefix_metrics] for row in hold[:model.OCCLUSION_WINDOWS[0][0]]])
        for candidate in (kf, fresh):
            b = np.array([[row[key] for key in prefix_metrics] for row in candidate[:model.OCCLUSION_WINDOWS[0][0]]])
            assert np.allclose(a, b, rtol=0, atol=1e-12)
    return dict(
        status="PASS", tested_seed=seed, expected_runs=len(tasks()), n=model.N,
        leaders=model.LEADER_COUNT, hidden_edges_per_follower=model.HIDDEN_PER_FOLLOWER,
        noise_levels=list(model.NOISE_LEVELS), R_values=[s ** 2 for s in model.NOISE_LEVELS],
        topology_edges=int(selected.sum()), hidden_edges=int(hidden.sum()),
        checks=[
            "common initial state across noise levels",
            "common standard-normal observation noise scaled by sigma",
            "R equals sigma squared at each noise level",
            "three 20-step windows with 100-step gaps",
            "nine hidden edges per follower",
            "no method differences before first blackout",
        ],
    )


def worker(row, output_root):
    output_root = Path(output_root)
    input_path = output_root / "inputs" / f"seed{row['seed']}.npz"
    with np.load(input_path) as data:
        initial = tuple(data[key] for key in ("xy", "heading", "alpha", "u", "standard_noise"))
    cfg = model.config(noise_std=float(row["noise_std"]))
    initial = initial[:4] + (initial[4] * cfg.noise_std,)
    with threadpool_limits(limits=1):
        rows = model.simulate(int(row["seed"]), row["mode"], True, "leader_loss", initial=initial, cfg=cfg)
    output = output_root / row["data"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    return row | dict(status="COMPLETE", rows=len(rows), bytes=output.stat().st_size, sha256=sha256(output))


def prepare_output(out: Path, manifest):
    out.mkdir(parents=True, exist_ok=True)
    for name in ("runs", "inputs", "analysis", "figures", "source_snapshot", "checks"):
        (out / name).mkdir(exist_ok=True)
    planned_path = out / "planned_manifest.json"
    planned_ids = [row["run_id"] for row in manifest]
    if planned_path.exists():
        saved = json.loads(planned_path.read_text(encoding="utf-8"))
        if [row["run_id"] for row in saved] != planned_ids:
            raise RuntimeError("Existing output uses a different run manifest; refusing to mix runs.")
    else:
        write_json(planned_path, manifest)
    metadata_path = out / "metadata.json"
    if not metadata_path.exists():
        write_json(metadata_path, dict(
            status="PREPARED", experiment="Moderate observation noise and repeated short occlusion validation",
            n=model.N, leaders=model.LEADER_COUNT, followers=model.N-model.LEADER_COUNT,
            k_neighbors=model.K_NEIGHBORS, hidden_edges_per_follower=model.HIDDEN_PER_FOLLOWER,
            occlusion_windows=[list(window) for window in model.OCCLUSION_WINDOWS], total_steps=TOTAL_STEPS,
            noise_levels=list(model.NOISE_LEVELS), R_values=[s ** 2 for s in model.NOISE_LEVELS],
            q=model.config().q, methods=list(model.METHODS), conditions=["tracking"], masks=["leader_loss"],
            seeds=SEEDS, planned_runs=len(manifest), primary_metric="follower_target_error_deg",
            paired_design="Same initial state and standard-normal observation noise tensor within each seed; noise levels differ only by scale",
            rng="NumPy Generator(PCG64)", python=platform.python_version(), numpy=np.__version__,
        ))
    write_json(out / "leader_ids_by_seed.json", {str(seed): model.role_ids(seed)[0].tolist() for seed in SEEDS})
    for seed in SEEDS:
        path = out / "inputs" / f"seed{seed}.npz"
        if not path.exists():
            initial = model.initialize(seed, model.config(noise_std=1.0))
            np.savez_compressed(path, xy=initial[0], heading=initial[1], alpha=initial[2], u=initial[3], standard_noise=initial[4])
    for source in ROOT.glob("*.py"):
        shutil.copy2(source, out / "source_snapshot" / source.name)
    write_json(out / "source_hashes.json", {p.name: sha256(p) for p in sorted((out / "source_snapshot").glob("*.py"))})


def run(args):
    manifest = tasks()
    out = ROOT / "results" / args.name
    prepare_output(out, manifest)
    write_json(out / "preflight.json", preflight())
    if (out / "audit.json").exists() and json.loads((out / "audit.json").read_text(encoding="utf-8")).get("status") == "COMPLETE":
        print(f"Already complete: {out}", flush=True)
        return
    complete, pending, failures = [], [], []
    for row in manifest:
        path = out / row["data"]
        if not path.exists():
            pending.append(row)
        else:
            with path.open("r", encoding="utf-8-sig", newline="") as handle:
                rows = sum(1 for _ in csv.reader(handle)) - 1
            if rows != TOTAL_STEPS:
                raise RuntimeError(f"Existing partial trace must be inspected manually: {path}")
            complete.append(row | dict(status="COMPLETE", rows=rows, bytes=path.stat().st_size, sha256=sha256(path)))
    with (out / "run.log").open("a", encoding="utf-8", buffering=1) as log:
        log.write(f"START pending={len(pending)} reused={len(complete)}\n")
        if pending:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = {pool.submit(worker, row, str(out)): row for row in pending}
                for future in as_completed(futures):
                    row = futures[future]
                    try:
                        complete.append(future.result())
                    except Exception:
                        failure = row | dict(status="FAILED", error=traceback.format_exc())
                        failures.append(failure); log.write(json.dumps(failure, ensure_ascii=False) + "\n")
                    finished = len(complete) + len(failures)
                    if finished % 20 == 0 or finished == len(manifest):
                        message = f"{finished}/{len(manifest)} finished; failures={len(failures)}"
                        print(message, flush=True); log.write(message + "\n")
    if failures:
        write_json(out / "failures.json", failures)
        raise RuntimeError(f"{len(failures)} runs failed; completed files are preserved for resume.")
    complete.sort(key=lambda row: row["run_id"])
    with (out / "RUN_MANIFEST.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(complete[0])); writer.writeheader(); writer.writerows(complete)
    write_json(out / "audit.json", dict(status="COMPLETE", expected_runs=len(manifest), completed_runs=len(complete), failed_runs=0, independent_seed_count=len(SEEDS), preflight="PASS", analysis="PENDING"))
    metadata_path = out / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8")); metadata.update(status="SIMULATION_COMPLETE", completed_runs=len(complete)); metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"SIMULATION COMPLETE: {len(complete)} runs at {out}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--name", default="formal_v1"); parser.add_argument("--workers", type=int, default=6); run(parser.parse_args())
