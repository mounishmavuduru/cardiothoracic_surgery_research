r"""The primary endpoint re-scored against openCARP's labels instead of the stand-in's.

The manuscript's Conclusion claimed the openCARP agreement "bounds how much of the
null could be an artefact of the stand-in".  Concordance on the *label* does not
bound the *incremental AUC*, so that claim needed either withdrawal or the actual
experiment.  This is the actual experiment.

Nothing about the features changes: SFI and every competitor are closed-form given
the mesh and never see a label.  Only the outcome column is swapped, from the frozen
monodomain Mitchell--Schaeffer verdict to openCARP's verdict on the same subject, and
the identical GM1 evaluation is re-run.  If the null survives an independent
reaction--diffusion solver's labels, it is not an artefact of the stand-in.

The openCARP arm is less powered by construction (26 positives against 42), so a
wider interval here is expected and is not evidence of anything; what matters is
whether the point estimate moves and whether any arm clears the pre-registered gate.

Writes results/gm1_opencarp_labels.json.

Usage:  python scripts/gm1_opencarp_labels.py
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from asb.experiments.realcohort import cohort_records  # noqa: E402
from asb.experiments.uw_cohort import _gm1_on, uw_cohort_records  # noqa: E402

OPENCARP = "results/opencarp_full_cohort.json"
GATE = 0.05


def opencarp_labels() -> dict:
    with open(OPENCARP, encoding="utf-8") as fh:
        d = json.load(fh)
    return {r["subject"]: bool(r["inducible"]) for r in d["rows"]}


def relabel(records, labels):
    """Swap the outcome column, leaving every feature untouched."""
    out, missing = [], []
    for rec in records:
        if rec.subject not in labels:
            missing.append(rec.subject)
            continue
        out.append(dataclasses.replace(rec, inducible=labels[rec.subject]))
    if missing:
        raise SystemExit(f"{len(missing)} subjects have no openCARP verdict, e.g. {missing[:3]}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--out", default="results/gm1_opencarp_labels.json")
    a = ap.parse_args()

    labels = opencarp_labels()
    roney = relabel(cohort_records("data/roney", limit=None, n_jobs=2, cache_dir="outputs"), labels)
    uw = relabel(uw_cohort_records("data/uw_boyle", n_jobs=2, cache_dir="outputs"), labels)
    sets = {"roney": roney, "uw": uw, "combined": roney + uw}

    cfg = {"n_splits": 5, "seed": 0, "n_boot": a.n_boot}
    cohorts = {name: _gm1_on(recs, cfg) for name, recs in sets.items()}

    print(f"{'cohort':<10} {'n':>4} {'ind':>8}   competitors_vs_+SFI dAUC (DeLong p)")
    for name, c in cohorts.items():
        cv = c.get("competitors_vs_+SFI")
        if not cv:
            print(f"{name:<10} {c['n_subjects']:>4} {c['n_inducible']:>3}/{c['n_subjects']:<4}"
                  f"   (not evaluable: too few positives)")
            continue
        lr, gbt = cv["lr"], cv["gbt"]
        print(f"{name:<10} {c['n_subjects']:>4} {c['n_inducible']:>3}/{c['n_subjects']:<4} "
              f"lr {lr['grouped_delta_auc']:>+7.3f} ({lr['delong_p']:.3f})   "
              f"gbt {gbt['grouped_delta_auc']:>+7.3f} ({gbt['delong_p']:.3f})")

    cleared = []
    for name, c in cohorts.items():
        for fam in ("competitors_vs_+SFI", "fibrosis_vs_+SFI"):
            if fam not in c:
                continue
            for clf in ("lr", "gbt"):
                r = c[fam][clf]
                if r["grouped_delta_auc"] >= GATE and r["delong_p"] < 0.05:
                    cleared.append(f"{name}/{fam}/{clf}")

    out = {
        "description": (
            "GM1 primary endpoint re-scored against openCARP verdicts instead of the "
            "monodomain Mitchell-Schaeffer stand-in. Features are identical and label-free; "
            "only the outcome column differs. Tests whether the predictive null is an "
            "artefact of the stand-in labeller."
        ),
        "label_source": "opencarp",
        "gate": GATE,
        "cohorts": cohorts,
        "arms_clearing_gate_and_significant": cleared,
        "null_survives_independent_solver": len(cleared) == 0,
        "config": {"n_boot": a.n_boot, "n_splits": 5, "seed": 0},
    }
    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"\narms clearing the {GATE} gate at p<0.05: {cleared if cleared else 'none'}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
