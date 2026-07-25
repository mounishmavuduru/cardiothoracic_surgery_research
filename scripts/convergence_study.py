"""Large resolution-convergence study: do inducibility verdicts converge as the mesh is refined?

Usage: python scripts/convergence_study.py [roney|uw]   (default roney — the informative cohort;
UW is only 9.8% inducible (8/82, corrected substrate) so it has too few borderline cases to
stress-test convergence).

The n=8 pilot showed only 50% verdict agreement across {1500,2000,3000} on Roney. This pushes toward
FINER resolution on more subjects to answer: do verdicts STABILIZE as resolution rises (convergence),
or keep drifting? Resumable — one JSON per (subject,resolution) under outputs/convergence_<cohort>/,
skipped if present, so it survives container restarts.

Scope (honest): SPATIAL convergence at the labeller's fixed timestep. Full space+time convergence
against a validated solver is the openCARP Claude Science block; this narrows, not closes. Real hearts
are finite (~144 total); the knob is nodes/mesh, capped by CPU. 'A million' is not applicable.
"""
import glob
import json
import os
import sys

import numpy as np

from asb.experiments.realcohort import FROZEN_BURST_CLS, FROZEN_MONO, _subject_seed
from asb.labels.monodomain import induce_monodomain
from asb.substrate.roney import coarsen_mesh, load_roney_mesh
from asb.substrate.uw_boyle import load_uw_mesh, uw_preablation_paths

RES = [1500, 3000, 6000, 12000, 24000, 48000]
N_SUBJECTS = 24
COHORT = sys.argv[1] if len(sys.argv) > 1 else "roney"
OUTDIR = f"outputs/convergence_{COHORT}"


def _cohort_paths():
    if COHORT == "roney":
        return sorted(glob.glob("data/roney/Mesh_*.vtk"), key=os.path.getsize)[:N_SUBJECTS]
    return sorted(uw_preablation_paths(), key=os.path.getsize)[:N_SUBJECTS]


def _load_and_name(path):
    if COHORT == "roney":
        return load_roney_mesh(path), os.path.basename(path)          # matches main Roney seeds
    return load_uw_mesh(path), "uw_" + os.path.basename(path)          # matches UW cohort seeds


def _key(path, res):
    return f"{os.path.basename(path).replace('.vtk','')}__{res}"


def _run_one(path, res):
    import time
    t0 = time.time()
    mesh, name = _load_and_name(path)
    coarse = coarsen_mesh(mesh, res)
    label = induce_monodomain(coarse, FROZEN_MONO, np.random.default_rng(_subject_seed(name)),
                              burst_cls=FROZEN_BURST_CLS)
    return {"subject": os.path.basename(path), "cohort": COHORT, "target_res": res,
            "n_nodes": int(coarse.n_points), "inducible": bool(label.inducible),
            "sustained_ms": float(label.meta.get("best_sustained_ms", 0.0)),
            "seconds": round(time.time() - t0, 1)}


def _work(args):
    p, r = args
    try:
        rec = _run_one(p, r)
    except Exception as e:  # noqa: BLE001
        rec = {"subject": os.path.basename(p), "target_res": r, "error": str(e)}
    with open(os.path.join(OUTDIR, _key(p, r) + ".json"), "w") as fh:
        json.dump(rec, fh)
    return rec


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    pairs = [(p, r) for r in RES for p in _cohort_paths()]           # coarse-to-fine
    todo = [(p, r) for (p, r) in pairs
            if not os.path.exists(os.path.join(OUTDIR, _key(p, r) + ".json"))]
    print(f"[convergence:{COHORT}] {len(todo)} (subject,res) pairs "
          f"(res={RES}, {N_SUBJECTS} subjects) -> {OUTDIR}", flush=True)
    from multiprocessing import Pool
    with Pool(4) as pool:
        for rec in pool.imap_unordered(_work, todo):
            if "error" in rec:
                print(f"  ERR {rec['subject']} res={rec['target_res']}: {rec['error'][:60]}", flush=True)
            else:
                print(f"  {rec['subject']} n={rec['n_nodes']} res={rec['target_res']} "
                      f"inducible={rec['inducible']} ({rec['seconds']}s)", flush=True)
    print("CONVERGENCE_BUILD_DONE", flush=True)


if __name__ == "__main__":
    main()
