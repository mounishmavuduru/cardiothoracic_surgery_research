"""Per-cohort openCARP/monodomain concordance, chance-corrected.

The manuscript reported a single pooled agreement of 156/182 (85.7%).  That
figure has two problems a referee will raise immediately: it pools the UW
cohort, on which the openCARP calibration gate is *not* met, and on an outcome
this imbalanced raw agreement is dominated by shared negatives.  This script
splits the concordance by cohort and reports Cohen's kappa and an exact McNemar
test alongside the raw figure, so the paper can say what the agreement is worth
rather than only how large it is.

Writes results/opencarp_concordance_detail.json.
"""
from __future__ import annotations

import json
import math
import os
from typing import Dict, List, Tuple

OPENCARP = "results/opencarp_full_cohort.json"
MONO_CACHE = [
    "outputs/cohort_records_223f1d7c52.json",
    "outputs/uw_cohort_records_223f1d7c52.json",
]
OUT = "results/opencarp_concordance_detail.json"


def load_monodomain() -> Dict[str, bool]:
    """subject -> inducible, from the cached records the manuscript's labels use."""
    labels: Dict[str, bool] = {}
    for path in MONO_CACHE:
        if not os.path.exists(path):
            raise SystemExit(f"missing monodomain cache: {path}")
        with open(path, encoding="utf-8") as fh:
            payload = json.load(fh)
        recs = payload["records"] if isinstance(payload, dict) and "records" in payload else payload
        for rec in recs:
            labels[rec["subject"]] = bool(rec["inducible"])
    return labels


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar on the discordant pairs (binomial, p=0.5)."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2 ** n)
    return float(min(1.0, 2.0 * tail))


def cohen_kappa(a: int, b: int, c: int, d: int) -> float:
    """a=both pos, b=openCARP only, c=monodomain only, d=both neg."""
    n = a + b + c + d
    if n == 0:
        return float("nan")
    po = (a + d) / n
    p_yes = ((a + b) / n) * ((a + c) / n)
    p_no = ((c + d) / n) * ((b + d) / n)
    pe = p_yes + p_no
    if abs(1.0 - pe) < 1e-12:
        return float("nan")
    return float((po - pe) / (1.0 - pe))


def summarise(pairs: List[Tuple[bool, bool]]) -> Dict[str, object]:
    """pairs = [(opencarp_inducible, monodomain_inducible), ...]"""
    a = sum(1 for o, m in pairs if o and m)
    b = sum(1 for o, m in pairs if o and not m)
    c = sum(1 for o, m in pairs if (not o) and m)
    d = sum(1 for o, m in pairs if (not o) and (not m))
    n = len(pairs)
    return {
        "n": n,
        "both_inducible": a,
        "only_opencarp": b,
        "only_monodomain": c,
        "both_not_inducible": d,
        "n_agree": a + d,
        "raw_agreement": (a + d) / n if n else float("nan"),
        "agreement_on_positives": a / (a + b + c) if (a + b + c) else float("nan"),
        "cohens_kappa": cohen_kappa(a, b, c, d),
        "mcnemar_exact_p": mcnemar_exact(b, c),
        "opencarp_rate": (a + b) / n if n else float("nan"),
        "monodomain_rate": (a + c) / n if n else float("nan"),
    }


def main() -> None:
    with open(OPENCARP, encoding="utf-8") as fh:
        carp = json.load(fh)
    mono = load_monodomain()

    by_cohort: Dict[str, List[Tuple[bool, bool]]] = {}
    missing = []
    for row in carp["rows"]:
        subj = row["subject"]
        if subj not in mono:
            missing.append(subj)
            continue
        by_cohort.setdefault(row["cohort"], []).append((bool(row["inducible"]), mono[subj]))

    if missing:
        raise SystemExit(f"{len(missing)} openCARP subjects have no monodomain label, e.g. {missing[:3]}")

    pooled = [p for v in by_cohort.values() for p in v]
    out = {
        "description": (
            "openCARP vs monodomain concordance, split by cohort and chance-corrected. "
            "The pooled figure spans the UW cohort, on which the openCARP calibration "
            "gate is not met."
        ),
        "pooled": summarise(pooled),
        "by_cohort": {k: summarise(v) for k, v in sorted(by_cohort.items())},
        "gate_passed": {k: carp["gate"][k]["gate_passed"] for k in carp["gate"]},
    }

    os.makedirs("results", exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)

    def show(label: str, s: Dict[str, object]) -> None:
        print(
            f"{label:<10} n={s['n']:<4} agree={s['n_agree']}/{s['n']} "
            f"({s['raw_agreement']:.1%})  kappa={s['cohens_kappa']:+.3f}  "
            f"McNemar p={s['mcnemar_exact_p']:.4f}  "
            f"carp={s['opencarp_rate']:.1%} mono={s['monodomain_rate']:.1%}"
        )

    print("openCARP / monodomain concordance\n")
    show("pooled", out["pooled"])
    for k, s in out["by_cohort"].items():
        show(k, s)
    print(f"\ncalibration gate passed: {out['gate_passed']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
