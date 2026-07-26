r"""Re-pass the pre-registered calibration gate on openCARP, and measure concordance.

`docs/PRE_REGISTRATION.md` section 7.1 records openCARP as the intended ground truth and the
monodomain Mitchell-Schaeffer solver as a substitution forced by its unavailability. openCARP
now runs (`scripts/opencarp_validate.sh`), so the substitution can be retired -- but only
after the same calibration gate section 7.4 applied to the substitute is re-passed on it:

    (a) inducible fraction in the 10-40 % band, and
    (b) a POSITIVE rank correlation between fibrosis burden and the verdict.

A labeller that fails either is not usable, whichever solver produced it.

The script also reports something worth having in its own right: per-subject CONCORDANCE
between the two solvers on identical meshes. The monodomain labels are the ones every
in-silico result in the manuscript rests on, and their agreement with a standard
reaction-diffusion solver is a direct measure of how much that substitution cost. Given the
already-documented per-subject instability of the monodomain labeller -- destroying a fibre
field flips 10/62 verdicts while leaving the rate unchanged -- concordance is expected to be
imperfect, and how imperfect is the number that matters.

Both solvers see the same coarsened mesh, the same fibrosis field and the same frozen
parameters. Labels are simulator verdicts and are never joined to the clinical outcomes.

Usage:  python scripts/opencarp_calibrate.py [--roney 30] [--uw 30] [--jobs 2]
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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
warnings.simplefilter("ignore")


def _run(job: Tuple[str, str]) -> Dict[str, object]:
    cohort, path = job
    warnings.simplefilter("ignore")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    import tempfile
    from asb.experiments.realcohort import COARSEN_NODES, _subject_seed
    from asb.labels.opencarp_run import OpenCARPConfig, label_with_opencarp
    from asb.substrate.roney import coarsen_mesh, load_roney_mesh
    from asb.substrate.uw_boyle import load_uw_mesh

    name = ("uw_" if cohort == "uw" else "") + os.path.basename(path)
    mesh = load_uw_mesh(path) if cohort == "uw" else load_roney_mesh(path)
    coarse = coarsen_mesh(mesh, COARSEN_NODES)
    wd = tempfile.mkdtemp(prefix="asb_cal_")
    t0 = time.time()
    try:
        lab = label_with_opencarp(coarse, OpenCARPConfig(),
                                  np.random.default_rng(_subject_seed(name)),
                                  workdir=wd, keep=False, timeout_s=1800)
        out = {"cohort": cohort, "subject": name, "inducible": bool(lab.inducible),
               "fib_burden": float(np.mean(coarse.fibrosis)),
               "seconds": round(time.time() - t0, 1), **{k: v for k, v in lab.meta.items()
                                                         if k not in ("solver",)}}
    except Exception as exc:                                   # noqa: BLE001
        out = {"cohort": cohort, "subject": name, "error": str(exc)[:300]}
    finally:
        import shutil
        shutil.rmtree(wd, ignore_errors=True)
    return out


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Rank correlation without a scipy dependency in the worker path."""
    if x.size < 3 or len(set(y.tolist())) < 2:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    d = float(np.sqrt((rx ** 2).sum() * (ry ** 2).sum()))
    return float((rx * ry).sum() / d) if d > 0 else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--roney", type=int, default=30)
    ap.add_argument("--uw", type=int, default=30)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--out", default="results/opencarp_calibration.json")
    a = ap.parse_args()

    import glob
    rp = sorted(glob.glob("data/roney/Mesh_*.vtk"),
                key=lambda p: (os.path.getsize(p), p))[: a.roney]
    up = sorted(glob.glob("data/uw_boyle/*_PreAbl_rev1.vtk"))[: a.uw]
    jobs = [("roney", p) for p in rp] + [("uw", p) for p in up]
    print(f"[calibrate] {len(jobs)} openCARP runs, {a.jobs} workers", flush=True)

    rows: List[Dict[str, object]] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, r in enumerate(ex.map(_run, jobs), 1):
            rows.append(r)
            if i % 10 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  elapsed {(time.time()-t0)/60:.1f} min", flush=True)

    errs = [r for r in rows if "error" in r]
    if errs:
        print(f"\n{len(errs)} runs errored; first: {errs[0]['error'][:200]}")
    ok = [r for r in rows if "error" not in r]

    # ---- gate ----------------------------------------------------------------
    gate: Dict[str, object] = {}
    print(f"\n{'cohort':<8} {'inducible':>12} {'rate':>8} {'Spearman(fib,label)':>21} {'gate':>8}")
    for cohort in ("roney", "uw"):
        sel = [r for r in ok if r["cohort"] == cohort]
        if not sel:
            continue
        y = np.array([1.0 if r["inducible"] else 0.0 for r in sel])
        f = np.array([r["fib_burden"] for r in sel])
        rate = float(y.mean())
        rho = _spearman(f, y)
        passed = bool(0.10 <= rate <= 0.40 and (rho > 0 if rho == rho else False))
        gate[cohort] = {"n": len(sel), "n_inducible": int(y.sum()),
                        "inducible_fraction": rate, "spearman_fibrosis": rho,
                        "gate_passed": passed}
        print(f"{cohort:<8} {int(y.sum()):>5}/{len(sel):<6} {rate:>7.1%} "
              f"{rho:>21.3f} {'PASS' if passed else 'FAIL':>8}")

    # ---- concordance with the monodomain labels the manuscript rests on -------
    conc: Dict[str, object] = {}
    try:
        from asb.experiments.realcohort import cohort_records
        from asb.experiments.uw_cohort import uw_cohort_records
        mono = {}
        for rec in cohort_records("data/roney", limit=None, n_jobs=1,
                                  cache_dir="outputs", verbose=False):
            mono[rec.subject] = bool(rec.inducible)
        for rec in uw_cohort_records("data/uw_boyle", n_jobs=1, cache_dir="outputs",
                                     verbose=False):
            mono[rec.subject] = bool(rec.inducible)
        pairs = [(r, mono[r["subject"]]) for r in ok if r["subject"] in mono]
        if pairs:
            agree = sum(1 for r, m in pairs if bool(r["inducible"]) == m)
            both = sum(1 for r, m in pairs if r["inducible"] and m)
            only_c = sum(1 for r, m in pairs if r["inducible"] and not m)
            only_m = sum(1 for r, m in pairs if (not r["inducible"]) and m)
            conc = {"n_compared": len(pairs), "agreement": agree / len(pairs),
                    "both_inducible": both, "only_opencarp": only_c,
                    "only_monodomain": only_m}
            print(f"\nconcordance with monodomain on {len(pairs)} shared subjects: "
                  f"{agree}/{len(pairs)} = {agree/len(pairs):.1%}")
            print(f"  both inducible {both} | openCARP only {only_c} | monodomain only {only_m}")
        else:
            print("\n(no cached monodomain labels matched; concordance skipped)")
    except Exception as exc:                                    # noqa: BLE001
        print(f"\n(concordance unavailable: {str(exc)[:160]})")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({"description": "openCARP calibration gate and monodomain concordance",
                   "gate_band": [0.10, 0.40], "gate": gate,
                   "concordance_with_monodomain": conc,
                   "n_errors": len(errs), "rows": rows}, fh, indent=2)
    print(f"\nwrote {a.out}   total {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
