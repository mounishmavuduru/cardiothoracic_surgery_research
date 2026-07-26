r"""Is openCARP's zero on the UW cohort a solver failure, or the burden result repeating?

Gate attempt 3 passed on Roney (10/30 = 33.3 %, Spearman +0.481) and returned 0/30 on
UW. Read carelessly that is half a calibration failure. Read against what is already
established it is the expected answer: `results/fibrosis_burden_swap.json` showed on the
monodomain solver that fibrosis burden, not any substrate defect, sets the inducibility
rate, and that Roney rescaled to UW's burden collapses to 1/62 (1.6 %). UW carries ~30 %
less fibrosis than Roney (0.244 against 0.349), so a near-zero rate is what a correctly
calibrated solver should produce there.

Asserting that would be rationalising. This tests it: rescale UW's fibrosis up to Roney's
burden, exactly as the monodomain burden swap did, and run openCARP. If the rate rises
toward the Roney figure, UW's zero is a property of the cohort's fibrosis and not of the
solver, and the same conclusion has now been reached independently by two solvers. If it
stays at zero, openCARP genuinely fails on UW anatomy and must be reported as such.

Note which way the falsification runs: this can only support the burden explanation by
MOVING the rate. A null result here counts against it.

Usage:  python scripts/opencarp_uw_burden_check.py [--n 30] [--jobs 2]
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

#: Cohort mean fibrosis at the frozen coarsening (results/cohort_geometry_audit.json).
RONEY_MEAN = 0.349

ARMS = {"native": None, "at_roney_burden": RONEY_MEAN}


def _run(job: Tuple[str, str]) -> Dict[str, object]:
    arm, path = job
    warnings.simplefilter("ignore")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    import tempfile
    from asb.experiments.realcohort import COARSEN_NODES, _subject_seed
    from asb.labels.opencarp_run import OpenCARPConfig, label_with_opencarp
    from asb.substrate.roney import coarsen_mesh
    from asb.substrate.uw_boyle import load_uw_mesh

    name = "uw_" + os.path.basename(path)
    coarse = coarsen_mesh(load_uw_mesh(path), COARSEN_NODES)
    before = float(np.mean(coarse.fibrosis))
    target = ARMS[arm]
    if target is not None and before > 0:
        coarse.fibrosis = np.clip(coarse.fibrosis * (target / before), 0.0, 1.0)
    after = float(np.mean(coarse.fibrosis))

    wd = tempfile.mkdtemp(prefix="asb_uwb_")
    try:
        lab = label_with_opencarp(coarse, OpenCARPConfig(),
                                  np.random.default_rng(_subject_seed(name)),
                                  workdir=wd, keep=False, timeout_s=1800)
        return {"arm": arm, "subject": name, "inducible": bool(lab.inducible),
                "burden_before": before, "burden_after": after}
    except Exception as exc:                                    # noqa: BLE001
        return {"arm": arm, "subject": name, "error": str(exc)[:200]}
    finally:
        import shutil
        shutil.rmtree(wd, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--out", default="results/opencarp_uw_burden_check.json")
    a = ap.parse_args()

    paths = sorted(glob.glob("data/uw_boyle/*_PreAbl_rev1.vtk"))[: a.n]
    jobs = [(arm, p) for arm in ARMS for p in paths]
    print(f"[uw-burden] {len(paths)} UW meshes x {len(ARMS)} arms = {len(jobs)} openCARP runs",
          flush=True)

    rows: List[Dict[str, object]] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, r in enumerate(ex.map(_run, jobs), 1):
            rows.append(r)
            if i % 10 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  elapsed {(time.time()-t0)/60:.1f} min", flush=True)

    ok = [r for r in rows if "error" not in r]
    summary = {}
    print(f"\n{'arm':<18} {'inducible':>12} {'rate':>8} {'burden':>9}")
    for arm in ARMS:
        sel = [r for r in ok if r["arm"] == arm]
        if not sel:
            continue
        k = sum(bool(r["inducible"]) for r in sel)
        summary[arm] = {"n": len(sel), "n_inducible": k,
                        "inducible_fraction": k / len(sel),
                        "mean_burden": float(np.mean([r["burden_after"] for r in sel]))}
        s = summary[arm]
        print(f"{arm:<18} {k:>5}/{len(sel):<6} {s['inducible_fraction']:>7.1%} "
              f"{s['mean_burden']:>9.3f}")

    nat = summary.get("native", {}).get("inducible_fraction", 0.0)
    up = summary.get("at_roney_burden", {}).get("inducible_fraction", 0.0)
    print(f"\nRoney reference under the same solver: 10/30 = 33.3%")
    print("VERDICT:", "burden explains UW's low rate on openCARP too"
          if up > nat + 0.05 else
          "rescaling does NOT recover the rate -- openCARP may genuinely fail on UW anatomy")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({"description": "does raising UW fibrosis to Roney's burden recover "
                                  "openCARP inducibility?",
                   "roney_reference": {"n": 30, "n_inducible": 10, "rate": 10 / 30},
                   "summary": summary, "rows": rows}, fh, indent=2)
    print(f"\nwrote {a.out}   total {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
