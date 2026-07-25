r"""Why is the UW/Boyle cohort only ~7% inducible? A 2x2 substrate ablation.

The same frozen monodomain Mitchell-Schaeffer protocol gives 20/62 (32.3%) inducible
on the Roney cohort and 6/82 (7.3%) on UW/Boyle. Two defects in the released UW
substrate could each explain that, and this script separates them by brute force
rather than by argument.

Factor 1 -- OPENINGS. Element tag 164 is not myocardium: it is the set of caps sealing
the pulmonary veins and the mitral valve. Census over all 82 meshes: four to six large
components, each a topological disc (five in 75 meshes, four in five, six in two -- the
usual pulmonary-vein variation); borders fibrotic tag 115 on only ~0.08% of its incident
edges against 14-20% by chance; unchanged by ablation; removing it opens the surface by a
median of 7 boundary loops. The frozen loader mapped it to
fibrosis 0.5, i.e. half-conducting tissue, which seals the atrium's orifices so activation
crosses them instead of circling them -- and reentry anchored on those orifices is a
principal AF mechanism.

Factor 2 -- FIBRES. All 164 released meshes carry a ``fiber`` array that is a constant
(1, 0, 0); measured directional spread is exactly 0.0. Conduction anisotropy is applied
relative to the LOCAL fibre direction, so a single global direction degenerates into a
fixed coordinate bias and removes the fibre heterogeneity that seeds unidirectional block.

Arms::

    A  caps kept     + constant fibres   <- the frozen pipeline; must reproduce 6/82
    B  caps dropped  + constant fibres
    C  caps kept     + varying fibres
    D  caps dropped  + varying fibres

The "varying" fibre field is a PROBE, not a correction: it is the global +z direction
projected onto each vertex tangent plane. It is not anatomical and no result here should
be read as a measurement of these patients. Its only job is to test whether inducibility
is sensitive to fibre heterogeneity at all.

Everything else -- coarsening, seeds, membrane parameters, the induction battery -- is the
frozen configuration, untouched. Labels here are simulator verdicts and are never joined
to the clinical outcome column.

Usage:  python scripts/uw_substrate_ablation.py [--limit N] [--jobs K]
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

ARMS: Dict[str, Dict[str, object]] = {
    "A_caps_kept_constant_fibres":    {"drop": False, "fibres": "asis"},
    "B_caps_dropped_constant_fibres": {"drop": True,  "fibres": "asis"},
    "C_caps_kept_varying_fibres":     {"drop": False, "fibres": "tangential"},
    "D_caps_dropped_varying_fibres":  {"drop": True,  "fibres": "tangential"},
}


def _vertex_normals(points: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Area-weighted vertex normals."""
    p0, p1, p2 = points[faces[:, 0]], points[faces[:, 1]], points[faces[:, 2]]
    fn = np.cross(p1 - p0, p2 - p0)
    vn = np.zeros_like(points)
    for k in range(3):
        np.add.at(vn, faces[:, k], fn)
    nrm = np.linalg.norm(vn, axis=1, keepdims=True)
    return vn / np.maximum(nrm, 1e-12)


def _tangential_fibres(points: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Probe field: global +z projected onto each vertex tangent plane, unit length.

    Spatially varying because the surface normal varies. Not anatomical.
    """
    nrm = _vertex_normals(points, faces)
    d = np.array([0.0, 0.0, 1.0])
    t = d - (nrm @ d)[:, None] * nrm
    ln = np.linalg.norm(t, axis=1, keepdims=True)
    # Where +z is parallel to the normal the projection vanishes; use +x there.
    fallback = np.array([1.0, 0.0, 0.0]) - (nrm @ np.array([1.0, 0.0, 0.0]))[:, None] * nrm
    t = np.where(ln > 1e-6, t, fallback)
    return t / np.maximum(np.linalg.norm(t, axis=1, keepdims=True), 1e-12)


def _run(job: Tuple[str, str]) -> Dict[str, object]:
    """Worker: one (arm, mesh) cell of the design."""
    arm, path = job
    warnings.simplefilter("ignore")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

    from asb.experiments.realcohort import (COARSEN_NODES, FROZEN_BURST_CLS,
                                            FROZEN_MONO, _subject_seed)
    from asb.labels.monodomain import induce_monodomain
    from asb.substrate.roney import coarsen_mesh
    from asb.substrate.uw_boyle import DROP_TAGS, load_uw_mesh

    spec = ARMS[arm]
    name = "uw_" + os.path.basename(path)
    t0 = time.time()
    mesh = load_uw_mesh(path, drop_tags=(DROP_TAGS if spec["drop"] else ()))
    coarse = coarsen_mesh(mesh, COARSEN_NODES)
    if spec["fibres"] == "tangential":
        coarse.fibres = _tangential_fibres(coarse.points, coarse.faces)

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
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--out", default="results/uw_substrate_ablation.json")
    a = ap.parse_args()

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    from asb.substrate.uw_boyle import uw_preablation_paths

    paths = uw_preablation_paths("data/uw_boyle")
    if a.limit:
        paths = paths[: a.limit]
    jobs = [(arm, p) for arm in ARMS for p in paths]
    print(f"[ablation] {len(paths)} meshes x {len(ARMS)} arms = {len(jobs)} runs, "
          f"{a.jobs} workers", flush=True)

    rows: List[Dict[str, object]] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, r in enumerate(ex.map(_run, jobs), 1):
            rows.append(r)
            if i % 20 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  elapsed {(time.time()-t0)/60:.1f} min", flush=True)

    summary = {}
    print(f"\n{'arm':<34} {'inducible':>12} {'rate':>8} {'fib_mean':>9} {'fib_spread':>11}")
    for arm in ARMS:
        sel = [r for r in rows if r["arm"] == arm]
        k = sum(bool(r["inducible"]) for r in sel)
        summary[arm] = {
            "n": len(sel), "n_inducible": k, "inducible_fraction": k / max(len(sel), 1),
            "mean_fibrosis": float(np.mean([r["fibrosis_mean"] for r in sel])),
            "mean_fibre_spread": float(np.mean([r["fibre_spread"] for r in sel])),
        }
        s = summary[arm]
        print(f"{arm:<34} {k:>5}/{len(sel):<6} {s['inducible_fraction']:>7.1%} "
              f"{s['mean_fibrosis']:>9.3f} {s['mean_fibre_spread']:>11.3f}")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({
            "description": "2x2 substrate ablation isolating why the UW cohort is ~7% inducible",
            "reference_roney": {"n": 62, "n_inducible": 20, "inducible_fraction": 20 / 62},
            "frozen_arm": "A_caps_kept_constant_fibres",
            "recorded_frozen_result": {"n": 82, "n_inducible": 6},
            "caveat": "the 'varying' fibre field is a non-anatomical probe, not a correction",
            "summary": summary, "rows": rows,
        }, fh, indent=2)
    print(f"\nwrote {a.out}   total {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
