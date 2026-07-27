r"""Compute and store every derived statistic the manuscript quotes.

`scripts/verify_manuscript_numbers.py` checks that each number in the paper traces to a
stored result. Several quoted values are legitimately *derived* rather than measured --
confidence intervals on the convergence tiers, the dispersion of the validity radius, a
concordance count that the openCARP result file stores only as a fraction, and the largest
per-subject change in sustained-reentry duration under the fibre control. Asserting those
in the text without storing them makes them untraceable, which is the defect the checker
exists to find, so they are computed here from the primary results and written out.

Nothing here re-runs a simulation; every input is an existing file under results/.

Usage:  python scripts/derived_statistics.py
"""
from __future__ import annotations

import json
import os
import re

import numpy as np
from scipy.stats import beta


def clopper_pearson(k: int, n: int, alpha: float = 0.05):
    """Exact binomial interval. Degenerate ends are handled explicitly."""
    lo = float(beta.ppf(alpha / 2, k, n - k + 1)) if k > 0 else 0.0
    hi = float(beta.ppf(1 - alpha / 2, k + 1, n - k)) if k < n else 1.0
    return lo, hi


def convergence_intervals() -> dict:
    """Exact binomial intervals for each mesh-resolution tier."""
    d = json.load(open("results/convergence_summary_roney.json", encoding="utf-8"))
    tiers = {}
    # tolerate either a list of tier records or a rate-keyed mapping
    # The file stores {"inducibility_rate": {"1500": 0.583, ...}} with a shared n_subjects.
    recs = d.get("tiers") or d.get("rows") or []
    rates = d.get("inducibility_rate") or d.get("rates")
    if not recs and isinstance(rates, dict):
        n_all = int(d.get("n_subjects", 24))
        recs = []
        for k, v in rates.items():
            # each entry is either a bare rate or {"rate": r, "n": n}
            if isinstance(v, dict):
                recs.append({"nodes": int(k), "inducible_fraction": float(v["rate"]),
                             "n": int(v.get("n", n_all))})
            else:
                recs.append({"nodes": int(k), "inducible_fraction": float(v), "n": n_all})
    for r in recs:
        n = int(r.get("n") or r.get("n_subjects") or 24)
        frac = r.get("inducible_fraction")
        if frac is None:
            continue
        k = int(round(float(frac) * n))
        lo, hi = clopper_pearson(k, n)
        tiers[str(r.get("nodes") or r.get("n_nodes"))] = {
            "n": n, "n_inducible": k, "rate": float(frac),
            "ci95_low": lo, "ci95_high": hi,
        }
    return tiers


def rho_dispersion() -> dict:
    """Median, IQR and range of the validity radius across the scanned meshes."""
    txt = open("results/validity_radius_scan.txt", encoding="utf-8").read()
    vals = [float(m) for m in re.findall(r"\s(\d+\.\d)\s+False", txt)]
    if not vals:
        vals = [float(x) for x in re.findall(r"\s(\d{3,6}\.\d)\s", txt)]
    a = np.array(vals, dtype=float)
    return {"n_meshes": int(a.size), "median": float(np.median(a)),
            "iqr_low": float(np.percentile(a, 25)), "iqr_high": float(np.percentile(a, 75)),
            "min": float(a.min()), "max": float(a.max()),
            "rho_star": 3.0, "min_over_rho_star": float(a.min() / 3.0),
            "all_exceed_rho_star": bool((a > 3.0).all())}


def opencarp_concordance_count() -> dict:
    """The file stores an agreement fraction; the paper quotes the count."""
    d = json.load(open("results/opencarp_full_cohort.json", encoding="utf-8"))
    c = d.get("concordance_with_monodomain", {})
    n = int(c.get("n_compared", 0))
    agree = int(round(float(c.get("agreement", 0.0)) * n))
    return {"n_compared": n, "n_agree": agree, "agreement": c.get("agreement"),
            "both_inducible": c.get("both_inducible"),
            "only_opencarp": c.get("only_opencarp"),
            "only_monodomain": c.get("only_monodomain")}


def fibre_control_dynamics() -> dict:
    """Largest per-subject change in sustained duration when fibres are destroyed."""
    d = json.load(open("results/roney_fibre_control.json", encoding="utf-8"))
    a = {r["subject"]: r["sustained_ms"] for r in d["rows"] if r["arm"] == "R_A_real_fibres"}
    b = {r["subject"]: r["sustained_ms"] for r in d["rows"] if r["arm"] == "R_B_constant_fibres"}
    shared = sorted(set(a) & set(b))
    diff = np.array([abs(a[s] - b[s]) for s in shared], dtype=float)
    return {"n_subjects": len(shared),
            "n_changed": int((diff > 0).sum()),
            "max_abs_change_ms": float(diff.max()) if diff.size else 0.0,
            "median_abs_change_ms": float(np.median(diff)) if diff.size else 0.0}


def main() -> None:
    out = {}
    for name, fn in (("convergence_binomial_intervals", convergence_intervals),
                     ("validity_radius_dispersion", rho_dispersion),
                     ("opencarp_concordance_count", opencarp_concordance_count),
                     ("fibre_control_dynamics", fibre_control_dynamics)):
        try:
            out[name] = fn()
            print(f"{name}: ok")
        except Exception as exc:                                # noqa: BLE001
            out[name] = {"error": str(exc)[:200]}
            print(f"{name}: FAILED -- {str(exc)[:120]}")

    ci = out.get("convergence_binomial_intervals", {})
    if ci and "error" not in ci:
        print("\nconvergence tiers, exact binomial 95% intervals:")
        for nodes, v in sorted(ci.items(), key=lambda kv: int(kv[0])):
            print(f"  {nodes:>6} nodes  {v['n_inducible']:>2}/{v['n']}  {v['rate']*100:5.1f}%"
                  f"  [{v['ci95_low']*100:4.1f}, {v['ci95_high']*100:4.1f}]")
    r = out.get("validity_radius_dispersion", {})
    if r and "error" not in r:
        print(f"\nrho over {r['n_meshes']} meshes: median {r['median']:.0f}, "
              f"IQR [{r['iqr_low']:.0f}, {r['iqr_high']:.0f}], "
              f"range [{r['min']:.0f}, {r['max']:.0f}], "
              f"min is {r['min_over_rho_star']:.0f}x rho*")
    c = out.get("opencarp_concordance_count", {})
    if c and "error" not in c:
        print(f"concordance: {c['n_agree']}/{c['n_compared']}")
    f = out.get("fibre_control_dynamics", {})
    if f and "error" not in f:
        print(f"fibre control: {f['n_changed']}/{f['n_subjects']} subjects changed, "
              f"max |change| {f['max_abs_change_ms']:.1f} ms")

    os.makedirs("results", exist_ok=True)
    with open("results/derived_statistics.json", "w", encoding="utf-8") as fh:
        json.dump({"description": "statistics quoted in the manuscript that are derived "
                                  "from, rather than stored in, the primary result files",
                   **out}, fh, indent=2)
    print("\nwrote results/derived_statistics.json")


if __name__ == "__main__":
    main()
