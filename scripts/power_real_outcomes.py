"""Power for the real clinical-outcome endpoint (pre-registration section 8.3).

Reproduces the power table frozen in ``docs/PRE_REGISTRATION.md`` section 8.3 before any
feature-outcome association was computed. Only the *marginal* event counts from the UW
recurrence file are used (48 events / 34 non-events); marginals are needed to size a study
and are not a feature-outcome association, so running this pre-analysis does not constitute
peeking at the endpoint.

Model: two nested classifier scores (baseline, baseline+SFI) under a binormal ROC, sharing a
fraction ``r`` of their noise -- nested models are strongly correlated in practice, and the
paired DeLong test gains most of its sensitivity from exactly that correlation. Power is the
empirical rejection rate of the repo's Sun-Xu fast-midrank DeLong test.

Usage:  python scripts/power_real_outcomes.py
"""
from __future__ import annotations

import json
import os
from typing import Dict

import numpy as np
from scipy.stats import norm

from asb.evaluation import delong_test

OUT = "results/power_real_outcomes.json"

#: Marginals from the UW/Boyle 2-year recurrence file (see pre-registration section 8.2).
N_EVENTS = 48
N_NONEVENTS = 34

NSIM = 4000
ALPHA = 0.05
SEED = 20260724

BASE_AUCS = (0.65, 0.70, 0.75)
DELTAS = (0.03, 0.05, 0.07, 0.10, 0.15)
CORRELATIONS = (0.80, 0.50)


def _mu_for_auc(auc: float) -> float:
    """Binormal separation giving a target AUC with unit-variance scores."""
    return np.sqrt(2.0) * norm.ppf(auc)


def power(auc_base: float, delta: float, r: float, rng: np.random.Generator) -> float:
    """Empirical power of the paired DeLong test for an incremental ``delta`` in AUC."""
    y = np.concatenate([np.ones(N_EVENTS), np.zeros(N_NONEVENTS)])
    n = y.size
    mu_base, mu_aug = _mu_for_auc(auc_base), _mu_for_auc(auc_base + delta)
    shared, indep = np.sqrt(r), np.sqrt(1.0 - r)
    hits = 0
    for _ in range(NSIM):
        z = rng.standard_normal(n)
        base = mu_base * y + shared * z + indep * rng.standard_normal(n)
        aug = mu_aug * y + shared * z + indep * rng.standard_normal(n)
        res = delong_test(y, aug, base)
        p = res.get("p_value", res.get("p"))
        if p is not None and p < ALPHA:
            hits += 1
    return hits / NSIM


def required_n(delta: float, r: float, rng: np.random.Generator,
               target: float = 0.80, event_fraction: float = None) -> int:
    """Smallest total cohort size reaching ``target`` power, holding the event rate fixed.

    Answers the question a reader will ask about an underpowered null: not "was this
    study big enough" (it was not) but "how big would it have to be". The event
    fraction is held at the observed 48/82 so the answer is about THIS population.
    """
    if event_fraction is None:
        event_fraction = N_EVENTS / (N_EVENTS + N_NONEVENTS)
    lo, hi = 40, 4000
    best = hi
    while lo <= hi:
        mid = (lo + hi) // 2
        n1 = max(int(round(mid * event_fraction)), 2)
        n0 = max(mid - n1, 2)
        pw = _power_at(n1, n0, 0.70, delta, r, rng)
        if pw >= target:
            best = mid
            hi = mid - 1
        else:
            lo = mid + 1
    return best


def _power_at(n1: int, n0: int, auc_base: float, delta: float, r: float,
              rng: np.random.Generator, nsim: int = 800) -> float:
    """Power at an arbitrary (n1, n0); a lighter-weight sibling of :func:`power`."""
    y = np.concatenate([np.ones(n1), np.zeros(n0)])
    n = y.size
    mu_b, mu_a = _mu_for_auc(auc_base), _mu_for_auc(auc_base + delta)
    shared, indep = np.sqrt(r), np.sqrt(1.0 - r)
    hits = 0
    for _ in range(nsim):
        z = rng.standard_normal(n)
        base = mu_b * y + shared * z + indep * rng.standard_normal(n)
        aug = mu_a * y + shared * z + indep * rng.standard_normal(n)
        res = delong_test(y, aug, base)
        p = res.get("p_value", res.get("p"))
        if p is not None and p < ALPHA:
            hits += 1
    return hits / nsim


def main() -> None:
    rng = np.random.default_rng(SEED)
    print(f"n_events={N_EVENTS}  n_nonevents={N_NONEVENTS}  total={N_EVENTS + N_NONEVENTS}"
          f"  alpha={ALPHA}  nsim={NSIM}")
    print("power = P[DeLong two-sided p < alpha] for baseline+SFI vs baseline\n")

    grid: Dict[str, Dict[str, Dict[str, float]]] = {}
    for r in CORRELATIONS:
        print(f"--- shared-noise correlation r={r:.2f} ---")
        print(f"{'base AUC':>9} " + " ".join(f"{'d=+%.2f' % d:>8}" for d in DELTAS))
        grid[f"{r:.2f}"] = {}
        for auc_base in BASE_AUCS:
            row = [power(auc_base, d, r, rng) for d in DELTAS]
            grid[f"{r:.2f}"][f"{auc_base:.2f}"] = {f"{d:.2f}": v for d, v in zip(DELTAS, row)}
            print(f"{auc_base:>9.2f} " + " ".join(f"{v:>8.2f}" for v in row))
        print()

    print("cohort size needed for 80% power, holding the observed event rate "
          f"({N_EVENTS}/{N_EVENTS + N_NONEVENTS} = "
          f"{N_EVENTS / (N_EVENTS + N_NONEVENTS):.1%}) fixed:")
    print(f"{'target dAUC':>12} {'r=0.80':>10} {'r=0.50':>10}")
    required: Dict[str, Dict[str, int]] = {}
    for d in (0.05, 0.07, 0.10):
        row = [required_n(d, r, rng) for r in (0.80, 0.50)]
        required[f"{d:.2f}"] = {"r_0.80": int(row[0]), "r_0.50": int(row[1])}
        print(f"{d:>12.2f} {row[0]:>10,} {row[1]:>10,}")
    print(f"\n(this study has n = {N_EVENTS + N_NONEVENTS})")

    # Aggregate statistics only. The per-patient outcome column is restricted-use
    # human-subject data (pre-registration section 8.1) and is never written here.
    out = {
        "description": (
            "Power for the pre-registered clinical endpoint (pre-registration section 8.3), "
            "computed from the marginal event counts alone. No feature-outcome association "
            "is computed, inspected, or stored. Aggregate statistics only; the per-patient "
            "outcome column is restricted-use and is not written to this file."
        ),
        "n_events": N_EVENTS,
        "n_nonevents": N_NONEVENTS,
        "n_total": N_EVENTS + N_NONEVENTS,
        "event_rate": N_EVENTS / (N_EVENTS + N_NONEVENTS),
        "alpha": ALPHA,
        "n_simulations": NSIM,
        "seed": SEED,
        "power_grid": grid,
        "n_required_for_80pct_power": required,
        "power_at_prereg_gate_0.05": {
            "r_0.80_base_0.70": grid["0.80"]["0.70"]["0.05"],
            "r_0.50_base_0.70": grid["0.50"]["0.70"]["0.05"],
        },
        "minimum_detectable_dauc": [0.10, 0.15],
    }
    os.makedirs("results", exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
