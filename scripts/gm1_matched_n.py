r"""Like-for-like GM1: the corrected substrate at the ORIGINAL cohort sizes.

The full recompute (``results/gm1_expanded_metrics.json``) changed two things at once
relative to the previously reported table: the substrate was corrected (UW cap elements
dropped) AND the Roney cohort grew from 62 to 100 meshes, because all 100 are now on disk
where previously only 62 were. A reader comparing the two tables cannot tell which change
moved the numbers.

This script removes that confound by rerunning the identical analysis with Roney capped at
its original 62 meshes -- the same subset, since ``cohort_records`` orders deterministically
by (filesize, path) -- giving Roney 62 + UW 82 = 144, exactly the n of the reported table.
Any difference from the reported numbers is then attributable to the substrate fix alone.

Per-subject records are cached by ``_config_tag()``, so this reuses the labelling work
already done and costs only the evaluation.

Usage:  python scripts/gm1_matched_n.py [--roney-limit 62]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from asb.experiments.realcohort import cohort_records  # noqa: E402
from asb.experiments.uw_cohort import _gm1_on, uw_cohort_records  # noqa: E402

#: The previously reported table, computed on the SEALED-atrium substrate with Roney n=62.
#: Kept here so the comparison is explicit rather than left to the reader.
PREVIOUS = {
    "roney":    {"n": 62,  "lr": -0.030, "gbt": +0.036, "p_lr": 0.64, "p_gbt": 0.26},
    "uw":       {"n": 82,  "lr": -0.057, "gbt": +0.104, "p_lr": 0.29, "p_gbt": 0.66},
    "combined": {"n": 144, "lr": -0.004, "gbt": +0.051, "p_lr": 0.92, "p_gbt": 0.155},
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--roney-limit", type=int, default=62)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--out", default="results/gm1_matched_n.json")
    a = ap.parse_args()

    roney = cohort_records("data/roney", limit=a.roney_limit, n_jobs=2, cache_dir="outputs")
    uw = uw_cohort_records("data/uw_boyle", n_jobs=2, cache_dir="outputs")
    sets = {"roney": roney, "uw": uw, "combined": roney + uw}

    cfg = {"n_splits": 5, "seed": 0, "n_boot": a.n_boot}
    cohorts = {name: _gm1_on(recs, cfg) for name, recs in sets.items()}

    print(f"\n{'cohort':<10} {'n':>4} {'ind':>7}   competitors_vs_+SFI dAUC (p)")
    print(f"{'':<10} {'':>4} {'':>7}   {'lr':>18}  {'gbt':>18}   [previous, sealed substrate]")
    for name, c in cohorts.items():
        cv = c.get("competitors_vs_+SFI")
        prev = PREVIOUS[name]
        if not cv:
            print(f"{name:<10} {c['n_subjects']:>4} {c['n_inducible']:>7}   (not evaluable)")
            continue
        lr, gbt = cv["lr"], cv["gbt"]
        print(f"{name:<10} {c['n_subjects']:>4} {c['n_inducible']:>3}/{c['n_subjects']:<3} "
              f"{lr['grouped_delta_auc']:>+11.3f} ({lr['delong_p']:.3f}) "
              f"{gbt['grouped_delta_auc']:>+11.3f} ({gbt['delong_p']:.3f})"
              f"   [was {prev['lr']:+.3f} / {prev['gbt']:+.3f}]")

    out = {
        "description": ("GM1 on the CORRECTED substrate at the ORIGINAL cohort sizes, so the "
                        "substrate fix is isolated from the Roney cohort growing 62 -> 100"),
        "roney_limit": a.roney_limit,
        "previous_reported_sealed_substrate": PREVIOUS,
        "cohorts": cohorts,
        "config": {"n_boot": a.n_boot, "n_splits": 5, "seed": 0},
        "label_source": "monodomain_ms",
    }
    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
