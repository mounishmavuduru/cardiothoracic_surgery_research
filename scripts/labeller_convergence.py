"""Interim label-fidelity bridge: is the monodomain inducibility verdict resolution-stable?

Runs the frozen monodomain labeller on real meshes coarsened to several node counts and checks
whether the inducible/not verdict agrees across resolutions. Agreement => the label is not a
coarse-mesh artifact. This does NOT replace openCARP cross-validation (a different simulator);
it is internal-robustness evidence that narrows, not closes, the fidelity gap.
"""
import glob
import json
import os
from collections import defaultdict

import numpy as np

from asb.experiments.realcohort import FROZEN_BURST_CLS, FROZEN_MONO, _subject_seed
from asb.labels.monodomain import induce_monodomain
from asb.substrate.roney import coarsen_mesh, load_roney_mesh

RES = [1500, 2000, 3000]


def _one(args):
    path, target = args
    name = os.path.basename(path)
    mesh = load_roney_mesh(path)
    coarse = coarsen_mesh(mesh, target)
    label = induce_monodomain(coarse, FROZEN_MONO, np.random.default_rng(_subject_seed(name)),
                              burst_cls=FROZEN_BURST_CLS)
    return {"subject": name, "target": target, "n_nodes": int(coarse.n_points),
            "inducible": bool(label.inducible),
            "sustained_ms": float(label.meta.get("best_sustained_ms", 0.0))}


if __name__ == "__main__":
    from multiprocessing import Pool
    paths = sorted(glob.glob("data/roney/Mesh_*.vtk"), key=lambda p: (os.path.getsize(p), p))[:8]
    if not paths:
        print("no Roney meshes on disk; skipping"); raise SystemExit(0)
    jobs = [(p, r) for p in paths for r in RES]
    with Pool(3) as pool:
        recs = pool.map(_one, jobs)
    bysub = defaultdict(dict)
    for r in recs:
        bysub[r["subject"]][r["target"]] = r["inducible"]
    consistent = sum(1 for d in bysub.values() if len(set(d.values())) == 1)
    rate = {res: float(np.mean([r["inducible"] for r in recs if r["target"] == res])) for res in RES}
    out = {"n_subjects": len(bysub), "resolutions": RES,
           "verdict_consistent_subjects": consistent, "total_subjects": len(bysub),
           "verdict_agreement_frac": consistent / max(1, len(bysub)),
           "inducibility_rate_by_resolution": {str(k): v for k, v in rate.items()},
           "records": recs}
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/labeller_convergence.json", "w", encoding="utf-8"), indent=2)
    print(f"verdict consistent across {RES}: {consistent}/{len(bysub)} subjects "
          f"({consistent / max(1,len(bysub)):.0%})", flush=True)
    print("inducibility rate by resolution:", {k: round(v, 2) for k, v in rate.items()}, flush=True)
    print("CONVERGENCE_DONE", flush=True)
