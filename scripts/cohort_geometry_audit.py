r"""Why is the UW cohort so rarely inducible? Compare the two cohorts' GEOMETRY.

Three substrate hypotheses have been tested and refuted: sealed atrial orifices
(+2/82, not significant), a degenerate fibre field (no effect on the rate at all,
McNemar p = 1.0) and coarse fibrosis (the wrong direction -- binarising fibrosis
roughly doubles inducibility). What has NOT been checked is whether the two cohorts
present the solver with comparable tissue in the first place.

Two candidates, both measurable without running a single simulation:

1. FIBROSIS BURDEN. The frozen protocol's own calibration reports Spearman(fibrosis,
   sustained reentry) = 0.75, so burden drives the label. If UW simply carries less
   fibrosis than Roney, a lower rate is expected and is not an anomaly at all.

2. EFFECTIVE RESOLUTION. Both cohorts are coarsened to a FIXED 2000 nodes. If their
   surface areas differ, equal node counts mean unequal element sizes, hence unequal
   numerical conduction. That matters here more than it usually would: the six-tier
   convergence study found the labeller does not converge in the tested window
   (rates 58.3, 16.7, 25.0, 20.8, 25.0, 37.5 % from 1500 to 48000 nodes), so a
   resolution difference between cohorts maps directly onto a rate difference.

Reports both, per cohort, with the statistics needed to decide whether either is
large enough to matter. Pure geometry -- no solver, no labels.

Usage:  python scripts/cohort_geometry_audit.py [--limit 40]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
warnings.simplefilter("ignore")

from asb.experiments.realcohort import COARSEN_NODES  # noqa: E402
from asb.substrate.roney import coarsen_mesh, load_roney_mesh  # noqa: E402
from asb.substrate.uw_boyle import load_uw_mesh  # noqa: E402


def surface_area_mm2(points: np.ndarray, faces: np.ndarray) -> float:
    p0, p1, p2 = points[faces[:, 0]], points[faces[:, 1]], points[faces[:, 2]]
    return float(0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1).sum())


def mean_edge_mm(points: np.ndarray, faces: np.ndarray) -> float:
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    d = np.linalg.norm(points[e[:, 0]] - points[e[:, 1]], axis=1)
    return float(d.mean())


def describe(name: str, meshes) -> dict:
    rows = []
    for label, mesh in meshes:
        coarse = coarsen_mesh(mesh, COARSEN_NODES)
        rows.append({
            "subject": label,
            "area_mm2": surface_area_mm2(np.asarray(coarse.points), np.asarray(coarse.faces)),
            "n_nodes": int(coarse.n_points),
            "mean_edge_mm": mean_edge_mm(np.asarray(coarse.points), np.asarray(coarse.faces)),
            "fibrosis_mean": float(np.mean(coarse.fibrosis)),
            "fibrosis_p90": float(np.percentile(coarse.fibrosis, 90)),
            "frac_above_half": float(np.mean(coarse.fibrosis > 0.5)),
        })
    arr = {k: np.array([r[k] for r in rows], dtype=float)
           for k in rows[0] if k != "subject"}
    summary = {k: {"mean": float(v.mean()), "sd": float(v.std()),
                   "min": float(v.min()), "max": float(v.max())}
               for k, v in arr.items()}
    print(f"\n=== {name}  (n={len(rows)}) ===")
    for k in ("area_mm2", "n_nodes", "mean_edge_mm", "fibrosis_mean",
              "fibrosis_p90", "frac_above_half"):
        s = summary[k]
        print(f"  {k:<16} mean {s['mean']:>9.3f}  sd {s['sd']:>8.3f}  "
              f"range {s['min']:.3f} .. {s['max']:.3f}")
    return {"summary": summary, "rows": rows}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--out", default="results/cohort_geometry_audit.json")
    a = ap.parse_args()

    rp = sorted(glob.glob("data/roney/Mesh_*.vtk"),
                key=lambda p: (os.path.getsize(p), p))[: a.limit]
    up = sorted(glob.glob("data/uw_boyle/*_PreAbl_rev1.vtk"))[: a.limit]
    print(f"loading {len(rp)} Roney and {len(up)} UW meshes "
          f"(coarsened to {COARSEN_NODES} nodes, the frozen setting)")

    roney = describe("Roney", [(os.path.basename(p), load_roney_mesh(p)) for p in rp])
    uw = describe("UW/Boyle", [(os.path.basename(p), load_uw_mesh(p)) for p in up])

    print("\n=== ratio, UW relative to Roney ===")
    for k in ("area_mm2", "mean_edge_mm", "fibrosis_mean", "frac_above_half"):
        r = roney["summary"][k]["mean"]
        u = uw["summary"][k]["mean"]
        print(f"  {k:<16} Roney {r:>9.3f}   UW {u:>9.3f}   ratio {u / r if r else float('nan'):.3f}")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({"description": "cohort geometry and fibrosis-burden comparison",
                   "coarsen_nodes": COARSEN_NODES,
                   "roney": roney, "uw": uw}, fh, indent=2)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
