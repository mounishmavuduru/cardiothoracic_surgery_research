r"""Does COARSE fibrosis explain the UW cohort's low inducibility? The last candidate.

Two candidate causes of the UW/Boyle cohort's anomalous rate have already been tested and
neither survives as the explanation:

  * sealed atrial orifices (elemTag 164 mapped to half-conducting tissue) -- repairing them
    moves 6/82 to 8/82, not significant;
  * a degenerate fibre field -- destroying the fibre field on Roney, which has a real one,
    leaves the rate EXACTLY unchanged at 20/62 (McNemar p = 1.0). Refuted.

The remaining structural difference is the fibrosis representation. After the cap elements
are dropped, UW carries only tags 111 and 115, i.e. a BINARY fibrotic/not field averaged onto
vertices. Roney carries a continuous LGE intensity ratio with 276-1071 distinct levels. Graded
border zones -- partially coupled tissue between healthy and dense scar -- are exactly where
slow conduction and unidirectional block arise, so losing them is a mechanistically plausible
reason for reentry to stop forming.

Arms (Roney anatomy throughout, everything else frozen)::

    Q_A_continuous       fibrosis as shipped                      <- control
    Q_B_binary_half      fibrosis := 1[f > 0.5]                   coarse AND lower burden
    Q_C_binary_matched   fibrosis := 1[f > t], t chosen per       coarse, burden PRESERVED
                         subject so the mean burden is unchanged

Q_C is the one that matters: it destroys only the GRADATION while holding mean fibrosis
burden fixed, so a drop there cannot be dismissed as "less fibrosis". Q_B is included to show
what confounding burden with gradation would have looked like.

Labels are simulator verdicts and are never joined to the clinical outcome column.

Usage:  python scripts/roney_fibrosis_quantization.py [--limit 62] [--jobs 4]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple

import numpy as np

ARMS = ("Q_A_continuous", "Q_B_binary_half", "Q_C_binary_matched")


def _quantize(fib: np.ndarray, arm: str) -> np.ndarray:
    """Apply the arm's fibrosis transform. Q_C preserves the mean burden exactly."""
    if arm == "Q_A_continuous":
        return fib
    if arm == "Q_B_binary_half":
        return (fib > 0.5).astype(float)
    if arm == "Q_C_binary_matched":
        m = float(np.mean(fib))
        if m <= 0.0:
            return np.zeros_like(fib)
        if m >= 1.0:
            return np.ones_like(fib)
        # threshold at the (1 - m) quantile so P(f > t) == m
        t = float(np.quantile(fib, 1.0 - m))
        return (fib > t).astype(float)
    raise ValueError(arm)


def _run(job: Tuple[str, str]) -> Dict[str, object]:
    arm, path = job
    warnings.simplefilter("ignore")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

    from asb.experiments.realcohort import (COARSEN_NODES, FROZEN_BURST_CLS,
                                            FROZEN_MONO, _subject_seed)
    from asb.labels.monodomain import induce_monodomain
    from asb.substrate.roney import coarsen_mesh, load_roney_mesh

    name = os.path.basename(path)
    t0 = time.time()
    mesh = load_roney_mesh(path)
    coarse = coarsen_mesh(mesh, COARSEN_NODES)

    before = float(np.mean(coarse.fibrosis))
    n_levels_before = int(np.unique(np.round(coarse.fibrosis, 4)).size)
    coarse.fibrosis = _quantize(coarse.fibrosis, arm)
    after = float(np.mean(coarse.fibrosis))
    n_levels_after = int(np.unique(np.round(coarse.fibrosis, 4)).size)

    lab = induce_monodomain(coarse, FROZEN_MONO, np.random.default_rng(_subject_seed(name)),
                            burst_cls=FROZEN_BURST_CLS)
    return {
        "arm": arm, "subject": name, "n_nodes": int(coarse.n_points),
        "inducible": bool(lab.inducible),
        "sustained_ms": float(lab.meta.get("best_sustained_ms", 0.0)),
        "fibrosis_mean_before": before, "fibrosis_mean_after": after,
        "levels_before": n_levels_before, "levels_after": n_levels_after,
        "seconds": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=62)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", default="results/roney_fibrosis_quantization.json")
    a = ap.parse_args()

    paths = sorted(glob.glob("data/roney/Mesh_*.vtk"), key=lambda p: (os.path.getsize(p), p))
    if a.limit:
        paths = paths[: a.limit]
    jobs = [(arm, p) for arm in ARMS for p in paths]
    print(f"[quantization] {len(paths)} meshes x {len(ARMS)} arms = {len(jobs)} runs, "
          f"{a.jobs} workers", flush=True)

    rows: List[Dict[str, object]] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, r in enumerate(ex.map(_run, jobs), 1):
            rows.append(r)
            if i % 20 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  elapsed {(time.time()-t0)/60:.1f} min", flush=True)

    summary = {}
    print(f"\n{'arm':<22} {'inducible':>12} {'rate':>8} {'burden':>8} {'levels':>8}")
    for arm in ARMS:
        sel = [r for r in rows if r["arm"] == arm]
        k = sum(bool(r["inducible"]) for r in sel)
        summary[arm] = {
            "n": len(sel), "n_inducible": k, "inducible_fraction": k / max(len(sel), 1),
            "mean_burden": float(np.mean([r["fibrosis_mean_after"] for r in sel])),
            "mean_levels": float(np.mean([r["levels_after"] for r in sel])),
        }
        s = summary[arm]
        print(f"{arm:<22} {k:>5}/{len(sel):<6} {s['inducible_fraction']:>7.1%} "
              f"{s['mean_burden']:>8.3f} {s['mean_levels']:>8.1f}")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({
            "description": "does destroying fibrosis GRADATION (not burden) suppress reentry?",
            "control_arm": "Q_A_continuous",
            "key_arm": "Q_C_binary_matched",
            "uw_reference": {"n": 82, "n_inducible": 8,
                             "note": "corrected substrate; UW fibrosis is binary at cell level"},
            "summary": summary, "rows": rows,
        }, fh, indent=2)
    print(f"\nwrote {a.out}   total {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
