r"""Does the fibrosis bin count change the openCARP verdict? It must not.

openCARP assigns conductivity and membrane parameters through discrete ``gregion`` /
``imp_region`` element tags, so the continuous fibrosis field has to be binned before the
solver can see it. That is a discretisation this project has already shown to be dangerous:
`results/roney_fibrosis_quantization.json` found that binarising fibrosis roughly DOUBLES
inducibility on the monodomain solver (20/62 -> 40/62 at matched burden). A coarse bin count
would therefore import exactly that artefact and silently inflate every openCARP label.

The default is 20 bins. This sweeps 5 / 10 / 20 / 40 on a fixed subset and reports whether
the verdict is stable. Stability is a precondition for using these labels at all, not a
nice-to-have: if the rate moves materially with the bin count, the labels are measuring the
discretisation rather than the tissue, and the sweep says which bin count is safe.

Reports the inducible fraction per bin count and the per-subject verdict agreement against
the finest setting, which is the stricter test -- an unchanged aggregate rate can still hide
wholesale relabelling, as the fibre control demonstrated elsewhere in this project.

Usage:  python scripts/opencarp_bin_sensitivity.py [--n 20] [--jobs 2]
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

BIN_COUNTS = (5, 10, 20, 40)


def _run(job: Tuple[int, str]) -> Dict[str, object]:
    n_bins, path = job
    warnings.simplefilter("ignore")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    import tempfile
    from asb.experiments.realcohort import COARSEN_NODES, _subject_seed
    from asb.labels.opencarp_run import OpenCARPConfig, label_with_opencarp
    from asb.substrate.roney import coarsen_mesh, load_roney_mesh

    name = os.path.basename(path)
    coarse = coarsen_mesh(load_roney_mesh(path), COARSEN_NODES)
    wd = tempfile.mkdtemp(prefix="asb_bins_")
    try:
        lab = label_with_opencarp(coarse, OpenCARPConfig(n_bins=n_bins),
                                  np.random.default_rng(_subject_seed(name)),
                                  workdir=wd, keep=False, timeout_s=1800)
        return {"n_bins": n_bins, "subject": name, "inducible": bool(lab.inducible),
                "fib_burden": float(np.mean(coarse.fibrosis))}
    except Exception as exc:                                    # noqa: BLE001
        return {"n_bins": n_bins, "subject": name, "error": str(exc)[:200]}
    finally:
        import shutil
        shutil.rmtree(wd, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20, help="meshes per bin count")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--out", default="results/opencarp_bin_sensitivity.json")
    a = ap.parse_args()

    paths = sorted(glob.glob("data/roney/Mesh_*.vtk"),
                   key=lambda p: (os.path.getsize(p), p))[: a.n]
    jobs = [(b, p) for b in BIN_COUNTS for p in paths]
    print(f"[bin-sweep] {len(paths)} meshes x {len(BIN_COUNTS)} bin counts = {len(jobs)} runs",
          flush=True)

    rows: List[Dict[str, object]] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, r in enumerate(ex.map(_run, jobs), 1):
            rows.append(r)
            if i % 10 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  elapsed {(time.time()-t0)/60:.1f} min", flush=True)

    ok = [r for r in rows if "error" not in r]
    finest = {r["subject"]: bool(r["inducible"])
              for r in ok if r["n_bins"] == max(BIN_COUNTS)}

    summary = {}
    print(f"\n{'n_bins':>7} {'inducible':>12} {'rate':>8} {'agreement vs finest':>21}")
    for b in BIN_COUNTS:
        sel = [r for r in ok if r["n_bins"] == b]
        if not sel:
            continue
        k = sum(bool(r["inducible"]) for r in sel)
        shared = [r for r in sel if r["subject"] in finest]
        agree = sum(1 for r in shared if bool(r["inducible"]) == finest[r["subject"]])
        summary[b] = {"n": len(sel), "n_inducible": k,
                      "inducible_fraction": k / max(len(sel), 1),
                      "agreement_vs_finest": agree / max(len(shared), 1)}
        print(f"{b:>7} {k:>5}/{len(sel):<6} {k/max(len(sel),1):>7.1%} "
              f"{agree}/{len(shared)} = {agree/max(len(shared),1):>9.1%}")

    rates = [summary[b]["inducible_fraction"] for b in summary]
    spread = (max(rates) - min(rates)) if rates else 0.0
    stable = spread <= 0.10
    print(f"\nraw spread across ALL bin counts: {spread:.1%}  -> "
          f"{'STABLE (<=10 pp)' if stable else 'UNSTABLE (>10 pp) by the raw test'}")

    # The raw spread is a crude proxy and is reported above exactly as it was
    # pre-specified, pass or fail. But the question a discretisation study actually asks
    # is whether the answer CONVERGES under refinement, and a coarse-limit outlier does
    # not bear on that -- the mesh-convergence study has the same shape, with its 1500-node
    # tier far off the rest. So the refinement-limit reading is reported alongside, never
    # instead: it is only meaningful if the finer settings agree with each other.
    refined = [b for b in sorted(summary) if b >= 10]
    conv = False
    if len(refined) >= 2:
        r_rates = [summary[b]["inducible_fraction"] for b in refined]
        r_spread = max(r_rates) - min(r_rates)
        agree = min(summary[b]["agreement_vs_finest"] for b in refined)
        conv = r_spread <= 0.05 and agree >= 0.95
        print(f"refinement limit (bins >= 10): spread {r_spread:.1%}, "
              f"worst per-subject agreement {agree:.1%}  -> "
              f"{'CONVERGED' if conv else 'NOT CONVERGED'}")
        print(f"operating point is {20} bins, which is "
              f"{'inside' if conv and 20 in refined else 'NOT inside'} the converged regime")

    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump({"description": "openCARP verdict sensitivity to the fibrosis bin count",
                   "bin_counts": list(BIN_COUNTS), "rate_spread": spread,
                   "stable_raw_spread_test": bool(stable),
                   "converged_above_10_bins": bool(conv),
                   "operating_bins": 20,
                   "summary": {str(k): v for k, v in summary.items()},
                   "rows": rows}, fh, indent=2)
    print(f"\nwrote {a.out}   total {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
