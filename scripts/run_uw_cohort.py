"""Label the UW/Boyle cohort and run the expanded (cross-cohort) GM1."""
import time

from asb.experiments.uw_cohort import run_gm1_expanded

if __name__ == "__main__":
    t0 = time.time()
    m = run_gm1_expanded(n_jobs=2, n_boot=10000)
    print(f"=== UW EXPANDED GM1 DONE in {time.time() - t0:.0f}s ===", flush=True)
    for name, c in m["cohorts"].items():
        line = (f"[{name:9s}] n={c['n_subjects']:>3} inducible={c.get('n_inducible', 0):>3} "
                f"({c.get('inducible_fraction', 0):.2f})")
        cv = c.get("competitors_vs_+SFI")
        if cv:
            line += (f"  competitors_vs_+SFI: grouped dAUC={cv['grouped_delta_auc']:+.3f} "
                     f"DeLong_p={cv['delong_p']:.4f} {'MET' if cv['endpoint_met'] else 'no'}")
        print(line, flush=True)
    print("=== UW_DONE ===", flush=True)
