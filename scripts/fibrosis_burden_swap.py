r"""Does fibrosis BURDEN explain the UW cohort's low inducibility? Swap it both ways.

Three substrate hypotheses are dead (sealed orifices, degenerate fibres, coarse fibrosis).
The geometry audit then found something none of them covered: the two cohorts do not carry
comparable amounts of fibrosis. Coarsened to the frozen 2000 nodes, mean fibrosis is 0.349
on Roney and 0.244 on UW -- UW has about 30 % less. The frozen protocol's own calibration
reports Spearman(fibrosis, sustained reentry) = 0.75, so burden drives the label, and a
cohort with less of it should be less inducible. That is not an anomaly; it is the expected
behaviour of a working labeller.

This tests it in both directions, which is what makes it decisive rather than suggestive:

    R_A_roney_native      Roney, untouched                    <- control, must give 20/62
    R_B_roney_at_uw       Roney, fibrosis rescaled to UW's mean (0.244)
    U_A_uw_native         UW, untouched                       <- control, must give 8/82
    U_B_uw_at_roney       UW, fibrosis rescaled to Roney's mean (0.349)

Rescaling is multiplicative and clipped to [0, 1], so the spatial pattern and the gradation
are preserved and only the level moves. If burden is the explanation, R_B falls toward the
UW rate and U_B rises toward the Roney rate. If neither moves, burden is excluded too and
the anomaly survives a fourth control.

Labels are simulator verdicts and are never joined to the clinical outcome column.

Usage:  python scripts/fibrosis_burden_swap.py [--roney-limit 62] [--jobs 4]
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

#: Cohort mean fibrosis at the frozen coarsening, from results/cohort_geometry_audit.json.
RONEY_MEAN = 0.349
UW_MEAN = 0.244

ARMS = {
    "R_A_roney_native": ("roney", None),
    "R_B_roney_at_uw": ("roney", UW_MEAN),
    "U_A_uw_native": ("uw", None),
    "U_B_uw_at_roney": ("uw", RONEY_MEAN),
}


def _rescale(fib: np.ndarray, target_mean: float) -> np.ndarray:
    """Scale fibrosis multiplicatively to a target mean, preserving pattern and gradation."""
    cur = float(np.mean(fib))
    if cur <= 0 or target_mean is None:
        return fib
    return np.clip(fib * (target_mean / cur), 0.0, 1.0)


def _run(job: Tuple[str, str]) -> Dict[str, object]:
    arm, path = job
    warnings.simplefilter("ignore")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

    from asb.experiments.realcohort import (COARSEN_NODES, FROZEN_BURST_CLS,
                                            FROZEN_MONO, _subject_seed)
    from asb.labels.monodomain import induce_monodomain
    from asb.substrate.roney import coarsen_mesh, load_roney_mesh
    from asb.substrate.uw_boyle import load_uw_mesh

    cohort, target = ARMS[arm]
    name = ("uw_" if cohort == "uw" else "") + os.path.basename(path)
    t0 = time.time()
    mesh = load_uw_mesh(path) if cohort == "uw" else load_roney_mesh(path)
    coarse = coarsen_mesh(mesh, COARSEN_NODES)

    before = float(np.mean(coarse.fibrosis))
    if target is not None:
        coarse.fibrosis = _rescale(coarse.fibrosis, target)
    after = float(np.mean(coarse.fibrosis))

    lab = induce_monodomain(coarse, FROZEN_MONO, np.random.default_rng(_subject_seed(name)),
                            burst_cls=FROZEN_BURST_CLS)
    return {"arm": arm, "cohort": cohort, "subject": name,
            "inducible": bool(lab.inducible),
            "sustained_ms": float(lab.meta.get("best_sustained_ms", 0.0)),
            "fibrosis_before": before, "fibrosis_after": after,
            "seconds": round(time.time() - t0, 1)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--roney-limit", type=int, default=62)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", default="results/fibrosis_burden_swap.json")
    a = ap.parse_args()

    rp = sorted(glob.glob("data/roney/Mesh_*.vtk"),
                key=lambda p: (os.path.getsize(p), p))[: a.roney_limit]
    up = sorted(glob.glob("data/uw_boyle/*_PreAbl_rev1.vtk"))

    jobs: List[Tuple[str, str]] = []
    for arm, (cohort, _) in ARMS.items():
        jobs += [(arm, p) for p in (rp if cohort == "roney" else up)]
    print(f"[burden-swap] {len(jobs)} runs ({len(rp)} Roney x2, {len(up)} UW x2), "
          f"{a.jobs} workers", flush=True)

    rows: List[Dict[str, object]] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, r in enumerate(ex.map(_run, jobs), 1):
            rows.append(r)
            if i % 25 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  elapsed {(time.time()-t0)/60:.1f} min", flush=True)

    summary = {}
    print(f"\n{'arm':<20} {'inducible':>12} {'rate':>8} {'burden before':>14} {'after':>8}")
    for arm in ARMS:
        sel = [r for r in rows if r["arm"] == arm]
        k = sum(bool(r["inducible"]) for r in sel)
        summary[arm] = {
            "n": len(sel), "n_inducible": k, "inducible_fraction": k / max(len(sel), 1),
            "mean_burden_before": float(np.mean([r["fibrosis_before"] for r in sel])),
            "mean_burden_after": float(np.mean([r["fibrosis_after"] for r in sel])),
        }
        s = summary[arm]
        print(f"{arm:<20} {k:>5}/{len(sel):<6} {s['inducible_fraction']:>7.1%} "
              f"{s['mean_burden_before']:>14.3f} {s['mean_burden_after']:>8.3f}")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({"description": "swap fibrosis burden between cohorts, both directions",
                   "roney_mean": RONEY_MEAN, "uw_mean": UW_MEAN,
                   "controls": {"R_A_roney_native": "must reproduce 20/62",
                                "U_A_uw_native": "must reproduce 8/82"},
                   "summary": summary, "rows": rows}, fh, indent=2)
    print(f"\nwrote {a.out}   total {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
