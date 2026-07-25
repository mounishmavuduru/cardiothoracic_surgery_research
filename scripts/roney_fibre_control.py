r"""The clean causal test: strip the fibres from a substrate that HAS them.

The UW 2x2 (``scripts/uw_substrate_ablation.py``) varies two defects on a cohort where
both are present. This script runs the complementary and much cleaner experiment: take
the Roney cohort, which ships genuinely varying fibres, and destroy ONLY the fibre field,
holding anatomy, fibrosis and every solver parameter fixed.

Measured contrast between the two cohorts' released substrates:

    field                       Roney (Zenodo 5801337)      UW/Boyle (Dryad kkwh70sg0)
    mean fibre directional      0.945 - 0.979               0.000 (exactly)
      spread
    distinct fibre directions   ~one per vertex             ONE, mesh-wide
    fibrosis                    continuous IIR,             categorical elemTag,
                                276-1071 distinct levels    3 levels

Arms::

    R_A_real_fibres       untouched Roney            <- control; must reproduce 20/62
    R_B_constant_fibres   fibres := (1, 0, 0)        <- UW's degenerate field, same anatomy

If R_B collapses toward the UW cohort's 7%, fibre degeneracy alone is sufficient to
explain the anomaly. If R_B stays near 32%, fibres are not the driver and the sealed
atrial openings (UW tag 164) carry the blame instead.

Labels here are simulator verdicts. They are never joined to the clinical outcome column.

Usage:  python scripts/roney_fibre_control.py [--limit 62] [--jobs 7]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple

import numpy as np

ARMS = {
    "R_A_real_fibres": "real",
    "R_B_constant_fibres": "constant",
}


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
    if ARMS[arm] == "constant":
        coarse.fibres = np.tile(np.array([1.0, 0.0, 0.0]), (coarse.points.shape[0], 1))

    spread = float(np.mean(np.linalg.norm(coarse.fibres - coarse.fibres.mean(axis=0), axis=1)))
    lab = induce_monodomain(coarse, FROZEN_MONO, np.random.default_rng(_subject_seed(name)),
                            burst_cls=FROZEN_BURST_CLS)
    return {
        "arm": arm, "subject": name, "n_nodes": int(coarse.n_points),
        "inducible": bool(lab.inducible),
        "sustained_ms": float(lab.meta.get("best_sustained_ms", 0.0)),
        "fibrosis_mean": float(coarse.fibrosis.mean()),
        "fibre_spread": spread, "seconds": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=62)
    ap.add_argument("--jobs", type=int, default=7)
    ap.add_argument("--out", default="results/roney_fibre_control.json")
    a = ap.parse_args()

    import glob
    # Same deterministic ordering the frozen pipeline uses: (filesize, path).
    paths = sorted(glob.glob("data/roney/Mesh_*.vtk"), key=lambda p: (os.path.getsize(p), p))
    if a.limit:
        paths = paths[: a.limit]
    jobs = [(arm, p) for arm in ARMS for p in paths]
    print(f"[roney-control] {len(paths)} meshes x {len(ARMS)} arms = {len(jobs)} runs, "
          f"{a.jobs} workers", flush=True)

    rows: List[Dict[str, object]] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, r in enumerate(ex.map(_run, jobs), 1):
            rows.append(r)
            if i % 20 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  elapsed {(time.time()-t0)/60:.1f} min", flush=True)

    summary = {}
    print(f"\n{'arm':<24} {'inducible':>12} {'rate':>8} {'fib_mean':>9} {'fib_spread':>11}")
    for arm in ARMS:
        sel = [r for r in rows if r["arm"] == arm]
        k = sum(bool(r["inducible"]) for r in sel)
        summary[arm] = {
            "n": len(sel), "n_inducible": k, "inducible_fraction": k / max(len(sel), 1),
            "mean_fibrosis": float(np.mean([r["fibrosis_mean"] for r in sel])),
            "mean_fibre_spread": float(np.mean([r["fibre_spread"] for r in sel])),
        }
        s = summary[arm]
        print(f"{arm:<24} {k:>5}/{len(sel):<6} {s['inducible_fraction']:>7.1%} "
              f"{s['mean_fibrosis']:>9.3f} {s['mean_fibre_spread']:>11.3f}")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({
            "description": "destroy ONLY the fibre field on a cohort that has one",
            "control_arm": "R_A_real_fibres",
            "recorded_frozen_result": {"n": 62, "n_inducible": 20},
            "uw_reference": {"n": 82, "n_inducible": 6, "note": "degenerate fibres AND sealed orifices"},
            "summary": summary, "rows": rows,
        }, fh, indent=2)
    print(f"\nwrote {a.out}   total {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
