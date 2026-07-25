r"""Checkpointed 100k-network GM4 scale-up (resumable).

Extends the GM4 transfer experiment (``run_gm4_scaled``) to very large N so the
predictive-null and Fiedler-gradient-localizer conclusions can be read off a
*scaling curve* with shrinking confidence intervals rather than a single point.

Design notes
------------
- Same per-network distribution as ``gm4._scaled_net`` (seed ``70000 + k``,
  ``n_nodes = 500``, radius sweep ``[0.055, 0.095]``), so the first 2000 networks
  reproduce the earlier ``run_gm4_scaled(2000)`` run exactly — this is a clean
  superset scale-up, not a new experiment.
- The localizer origin ranks + spatial-null draws are computed **inline** while
  the network graph is still in memory, so the aggregation never rebuilds a
  network (rebuilding 35k inducible graphs would roughly double wall-clock).
- Work is sharded (default 5000 networks/shard) and each shard is written
  atomically to ``outputs/scaled100k/shard_XXXX.json``. Re-running skips shards
  already on disk, so the job survives container restarts and resumes.
- ``aggregate()`` reads all shards and reports predict + localize at a ladder of
  N values (the scaling curve). The heavy nested-CV/bootstrap predict test runs
  only at the ladder points, not per shard.
"""
from __future__ import annotations

import json
import os
from typing import Dict, List

import numpy as np

from asb.config import Config
from asb.experiments.gm1 import COMPETITOR_COLS, _evaluate_pair
from asb.experiments.gm2 import localizer_fields
from asb.experiments.gm3 import gm3_extra_features
from asb.experiments.gm4 import FROZEN_FHN, _perm_verdict, _scaled_net
from asb.experiments.realcohort import FROZEN_SFI
from asb.features import subject_features
from asb.transfer.fhn import induce_fhn

_KEYS = ("grad_phi2", "perron", "combined", "wdegree", "fibrosis")
_SEVEN = ("total", "region_max", "region_mean", "region_std", "top1", "top2", "top3")


def _process_big(k: int) -> Dict[str, object]:
    """One network: FHN label + subject/subspace features + inline localizer ranks."""
    G, radius, burden = _scaled_net(k)
    label = induce_fhn(G, FROZEN_FHN)
    feats = subject_features(G, Config(seed=0, sfi=FROZEN_SFI), np.random.default_rng(k))
    feats.update(gm3_extra_features(G, k_dim=2, include_exact=False))
    rec: Dict[str, object] = {
        "seed": int(k), "radius": float(radius), "burden": float(burden),
        "inducible": bool(label.inducible),
        "features": {kk: float(v) for kk, v in feats.items()},
    }
    if label.inducible and label.reentry_origin is not None:
        origin = int(label.reentry_origin)
        fields = localizer_fields(G)
        n = int(G.n_nodes)
        rng = np.random.default_rng(900000 + k)
        rand = rng.integers(0, n, size=min(200, n))
        loc = {}
        for key in _KEYS:
            fld = np.asarray(fields[key], float)
            loc[key] = {
                "rank": float(np.mean(fld >= fld[origin])),
                "perm": [float(np.mean(fld >= fld[int(v)])) for v in rand],
            }
        rec["origin"] = origin
        rec["loc"] = loc
    return rec


def run_100k(n_total: int = 100_000, *, shard: int = 5000, n_jobs: int = 4,
             out: str = "outputs/scaled100k") -> None:
    """Build+label ``n_total`` networks in resumable shards."""
    os.makedirs(out, exist_ok=True)
    n_shards = (n_total + shard - 1) // shard
    for si in range(n_shards):
        path = os.path.join(out, f"shard_{si:04d}.json")
        if os.path.exists(path) and os.path.getsize(path) > 2:
            continue
        lo, hi = si * shard, min((si + 1) * shard, n_total)
        if n_jobs > 1:
            from multiprocessing import Pool
            with Pool(n_jobs) as pool:
                recs = pool.map(_process_big, range(lo, hi))
        else:
            recs = [_process_big(k) for k in range(lo, hi)]
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(recs, fh)
        os.replace(tmp, path)
        n_ind = sum(int(r["inducible"]) for r in recs)
        print(f"[shard {si + 1}/{n_shards}] N={hi} (+{len(recs)}, {n_ind} unstable)",
              flush=True)
    print(f"ALL_SHARDS_DONE (n_total={n_total})", flush=True)


def _load(out: str, limit: int | None = None) -> List[dict]:
    recs: List[dict] = []
    for fn in sorted(os.listdir(out)):
        if not fn.startswith("shard_") or not fn.endswith(".json"):
            continue
        with open(os.path.join(out, fn)) as fh:
            recs.extend(json.load(fh))
        if limit is not None and len(recs) >= limit:
            return recs[:limit]
    return recs


def _predict_at(recs: List[dict], *, n_boot: int, n_splits: int, seed: int) -> dict:
    import pandas as pd
    y = np.array([int(r["inducible"]) for r in recs], int)
    groups = [f"snet_{r['seed']}" for r in recs]
    rows = [r["features"] for r in recs]
    X = pd.DataFrame(rows).reindex(sorted({k for row in rows for k in row}), axis=1).fillna(0.0)
    variants = {"single_vector": [f"sfi_{s}" for s in _SEVEN if f"sfi_{s}" in X.columns],
                "subspace": sorted(c for c in X.columns if c.startswith("sfisub_"))}
    cfg = {"n_splits": n_splits, "seed": seed, "n_boot": n_boot}
    out = {}
    if min(int(y.sum()), len(y) - int(y.sum())) >= 2:
        for vname, cols in variants.items():
            out[f"competitors_vs_+{vname}"] = _evaluate_pair(
                X, y, groups, COMPETITOR_COLS, cols, cfg)
    return out


def _localize_at(recs: List[dict], seed: int) -> dict:
    rng = np.random.default_rng(seed)
    subj = [r for r in recs if r.get("loc")]
    out = {"n_localized": len(subj), "fields": {}}
    for key in _KEYS:
        ranks = np.array([r["loc"][key]["rank"] for r in subj], float)
        perms = [np.asarray(r["loc"][key]["perm"], float) for r in subj]
        if ranks.size == 0:
            continue
        mean_rank, perm_p = _perm_verdict(ranks, perms, rng)
        out["fields"][key] = {"mean_origin_rank": float(mean_rank), "perm_p": float(perm_p),
                              "verdict": "KEEP" if (np.isfinite(perm_p) and perm_p < 0.05)
                              else "DELETE"}
    return out


def aggregate(out: str = "outputs/scaled100k", *, ladder=(2000, 5000, 10000, 25000, 50000,
              100000), n_boot: int = 5000, n_splits: int = 5, seed: int = 0,
              metrics_path: str = "results/gm4_100k_metrics.json") -> dict:
    """Read all shards; report predict + localize on a scaling ladder."""
    allrecs = _load(out)
    have = len(allrecs)
    curve = []
    for N in ladder:
        if N > have:
            break
        sub = allrecs[:N]
        n_ind = sum(int(r["inducible"]) for r in sub)
        pred = _predict_at(sub, n_boot=n_boot, n_splits=n_splits, seed=seed)
        loc = _localize_at(sub, seed)
        curve.append({"N": N, "n_unstable": n_ind, "predict": pred, "localize": loc})
        gp = loc["fields"].get("grad_phi2", {})
        print(f"[N={N:>6}] unstable={n_ind:>6} "
              f"grad_phi2 rank={gp.get('mean_origin_rank', float('nan')):.3f} "
              f"p={gp.get('perm_p', float('nan')):.4f} {gp.get('verdict', '')}", flush=True)
    result = {"n_built": have, "ladder": curve,
              "distribution": {"n_nodes": 500, "radius_range": [0.055, 0.095],
                               "seed_base": 70000}}
    os.makedirs(os.path.dirname(metrics_path), exist_ok=True)
    with open(metrics_path, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
    return result
