"""Large resolution-convergence study: do inducibility verdicts converge as the mesh is refined?

The n=8 pilot showed only 50% verdict agreement across {1500,2000,3000}. This pushes toward FINER
resolution (up to ~native) on more subjects to answer the real question: do verdicts STABILIZE as
resolution increases (convergence), or keep drifting? Resumable — one JSON per (subject,resolution)
under outputs/convergence/, skipped if present, so it survives container restarts.

Scope (honest): this is SPATIAL convergence at the labeller's fixed timestep. Full space+time
convergence against a validated solver is the openCARP Claude Science block; this narrows, not closes.
Real hearts are finite (~144 total); the knob here is nodes/mesh, capped by CPU (~native ~75k for the
0.5mm UW meshes). 'A million' is not applicable — there is no million-heart dataset.
"""
import glob
import json
import os

import numpy as np

from asb.experiments.realcohort import FROZEN_BURST_CLS, FROZEN_MONO, _subject_seed
from asb.labels.monodomain import induce_monodomain
from asb.substrate.roney import coarsen_mesh
from asb.substrate.uw_boyle import load_uw_mesh, uw_preablation_paths

RES = [1500, 3000, 6000, 12000, 24000, 48000]
N_SUBJECTS = 24
OUTDIR = "outputs/convergence"


def _pairs():
    paths = uw_preablation_paths()
    paths = sorted(paths, key=lambda p: os.path.getsize(p))[:N_SUBJECTS]
    # coarse-to-fine so early results arrive fast and the trend builds up
    return [(p, r) for r in RES for p in paths]


def _key(path, res):
    return f"{os.path.basename(path).replace('.vtk','')}__{res}"


def _run_one(path, res):
    name = "uw_" + os.path.basename(path)
    import time
    t0 = time.time()
    mesh = load_uw_mesh(path)
    coarse = coarsen_mesh(mesh, res)
    label = induce_monodomain(coarse, FROZEN_MONO, np.random.default_rng(_subject_seed(name)),
                              burst_cls=FROZEN_BURST_CLS)
    return {"subject": os.path.basename(path), "target_res": res, "n_nodes": int(coarse.n_points),
            "inducible": bool(label.inducible),
            "sustained_ms": float(label.meta.get("best_sustained_ms", 0.0)),
            "seconds": round(time.time() - t0, 1)}


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    todo = [(p, r) for (p, r) in _pairs()
            if not os.path.exists(os.path.join(OUTDIR, _key(p, r) + ".json"))]
    print(f"[convergence] {len(todo)} (subject,res) pairs to run "
          f"(res={RES}, {N_SUBJECTS} subjects)", flush=True)
    from multiprocessing import Pool

    def _work(args):
        p, r = args
        try:
            rec = _run_one(p, r)
        except Exception as e:  # noqa: BLE001
            rec = {"subject": os.path.basename(p), "target_res": r, "error": str(e)}
        with open(os.path.join(OUTDIR, _key(p, r) + ".json"), "w") as fh:
            json.dump(rec, fh)
        return rec

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
