r"""Why does openCARP return zero on the UW cohort? Check propagation itself.

openCARP labels 4/82 on UW against 22/100 on Roney, and raising UW's fibrosis burden to
Roney's level does not recover the rate -- unlike the monodomain solver, where the burden
swap moves it in both directions. That leaves a real discrepancy, and the first thing to
rule out is the most basic: that the wave propagates on UW anatomy at all.

There is a specific mechanism to suspect. Both cohorts are coarsened to a fixed 2000 nodes,
but UW meshes are ~18 % larger in surface area, so their elements are ~7 % longer. A
monodomain discretisation has a maximum element size above which a wavefront cannot be
resolved and propagation fails outright rather than merely slowing. If that threshold sits
between the two cohorts' element sizes, openCARP would report no reentry on UW for a
numerical reason that has nothing to do with the tissue -- and the label would be an
artefact, not a measurement.

Reports, per mesh, the peak fraction of tissue depolarised and the number of activation
wavefronts, alongside the mesh's mean edge length. Propagation failure and genuine
non-inducibility look completely different here: the first never activates more than the
stimulus footprint, the second activates widely and then dies.

Usage:  python scripts/opencarp_uw_propagation.py [--n 6]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import tempfile
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
warnings.simplefilter("ignore")

from asb.experiments.realcohort import COARSEN_NODES, _subject_seed  # noqa: E402
from asb.labels.opencarp_run import (OpenCARPConfig, export_for_opencarp,  # noqa: E402
                                     read_igb)
from asb.substrate.roney import coarsen_mesh, load_roney_mesh  # noqa: E402
from asb.substrate.uw_boyle import load_uw_mesh  # noqa: E402


def mean_edge_mm(points: np.ndarray, faces: np.ndarray) -> float:
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    return float(np.linalg.norm(points[e[:, 0]] - points[e[:, 1]], axis=1).mean())


def probe(cohort: str, path: str, cfg: OpenCARPConfig) -> dict:
    from asb.labels.opencarp_run import build_par, find_opencarp, _to_wsl_path
    import subprocess

    name = ("uw_" if cohort == "uw" else "") + os.path.basename(path)
    mesh = load_uw_mesh(path) if cohort == "uw" else load_roney_mesh(path)
    coarse = coarsen_mesh(mesh, COARSEN_NODES)
    edge = mean_edge_mm(np.asarray(coarse.points), np.asarray(coarse.faces))

    wd = tempfile.mkdtemp(prefix="asb_prop_")
    ex = export_for_opencarp(coarse, wd, cfg)
    pts_um = np.asarray(coarse.points, dtype=float) * 1.0e3
    rng = np.random.default_rng(_subject_seed(name))
    sites = rng.choice(pts_um.shape[0], size=min(cfg.n_pacing_sites, pts_um.shape[0]),
                       replace=False)
    par = os.path.join(wd, "run.par")
    with open(par, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(build_par(ex, cfg, [pts_um[s] for s in sites]))

    binary = find_opencarp()
    sim = os.path.join(wd, "out")
    lin = binary[4:]
    lib = os.path.dirname(os.path.dirname(lin)) + "/lib"
    cmd = (f'export LD_LIBRARY_PATH="{lib}:$LD_LIBRARY_PATH"; "{lin}" '
           f'+F "{_to_wsl_path(par)}" -meshname "{_to_wsl_path(ex["basename"])}" '
           f'-simID "{_to_wsl_path(sim)}"')
    subprocess.run(["wsl", "-e", "bash", "-lc", cmd], capture_output=True, text=True,
                   timeout=1800)

    igb = os.path.join(sim, "vm.igb")
    out = {"cohort": cohort, "subject": name, "n_nodes": int(coarse.n_points),
           "mean_edge_mm": edge, "fibrosis_mean": float(np.mean(coarse.fibrosis))}
    if not os.path.isfile(igb):
        out["error"] = "no vm.igb"
        return out
    vm = read_igb(igb)
    above = vm > 0.5
    frac = above.mean(axis=1)
    out.update({"peak_depolarised_fraction": float(frac.max()),
                "mean_depolarised_fraction": float(frac.mean()),
                "n_nodes_ever_activated": int(above.any(axis=0).sum()),
                "frac_nodes_ever_activated": float(above.any(axis=0).mean())})
    import shutil
    shutil.rmtree(wd, ignore_errors=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--out", default="results/opencarp_uw_propagation.json")
    a = ap.parse_args()

    cfg = OpenCARPConfig()
    rp = sorted(glob.glob("data/roney/Mesh_*.vtk"),
                key=lambda p: (os.path.getsize(p), p))[: a.n]
    up = sorted(glob.glob("data/uw_boyle/*_PreAbl_rev1.vtk"))[: a.n]

    rows = []
    print(f"{'cohort':<7} {'subject':<26} {'edge mm':>8} {'peak depol':>11} "
          f"{'ever activated':>15}")
    for cohort, paths in (("roney", rp), ("uw", up)):
        for p in paths:
            r = probe(cohort, p, cfg)
            rows.append(r)
            if "error" in r:
                print(f"{cohort:<7} {r['subject'][:26]:<26} {'':>8} {r['error']}")
                continue
            print(f"{cohort:<7} {r['subject'][:26]:<26} {r['mean_edge_mm']:>8.2f} "
                  f"{r['peak_depolarised_fraction']:>11.3f} "
                  f"{r['frac_nodes_ever_activated']:>15.3f}")

    ok = [r for r in rows if "error" not in r]
    for cohort in ("roney", "uw"):
        sel = [r for r in ok if r["cohort"] == cohort]
        if not sel:
            continue
        print(f"\n{cohort}: mean edge {np.mean([r['mean_edge_mm'] for r in sel]):.2f} mm, "
              f"peak depolarised {np.mean([r['peak_depolarised_fraction'] for r in sel]):.3f}, "
              f"nodes ever activated "
              f"{np.mean([r['frac_nodes_ever_activated'] for r in sel]):.3f}")
    u = [r for r in ok if r["cohort"] == "uw"]
    if u:
        act = np.mean([r["frac_nodes_ever_activated"] for r in u])
        print("\nVERDICT:", "propagation FAILS on UW -- the zero is a numerical artefact"
              if act < 0.25 else
              "UW propagates normally -- the zero is not a propagation failure")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({"description": "does openCARP propagate on UW anatomy at all?",
                   "rows": rows}, fh, indent=2)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
