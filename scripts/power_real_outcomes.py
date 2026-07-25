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

import numpy as np
from scipy.stats import norm

from asb.evaluation import delong_test

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


def main() -> None:
    rng = np.random.default_rng(SEED)
    print(f"n_events={N_EVENTS}  n_nonevents={N_NONEVENTS}  total={N_EVENTS + N_NONEVENTS}"
          f"  alpha={ALPHA}  nsim={NSIM}")
    print("power = P[DeLong two-sided p < alpha] for baseline+SFI vs baseline\n")
    for r in CORRELATIONS:
        print(f"--- shared-noise correlation r={r:.2f} ---")
        print(f"{'base AUC':>9} " + " ".join(f"{'d=+%.2f' % d:>8}" for d in DELTAS))
        for auc_base in BASE_AUCS:
            row = [power(auc_base, d, r, rng) for d in DELTAS]
            print(f"{auc_base:>9.2f} " + " ".join(f"{v:>8.2f}" for v in row))
        print()


if __name__ == "__main__":
    main()
