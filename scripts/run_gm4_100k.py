"""Launch/resume the 100k-network GM4 scale-up, then aggregate the scaling curve."""
import time

from asb.experiments.gm4_scale import aggregate, run_100k

if __name__ == "__main__":
    t0 = time.time()
    run_100k(n_total=100_000, shard=5000, n_jobs=4, out="outputs/scaled100k")
    print(f"=== BUILD DONE in {time.time() - t0:.0f}s ===", flush=True)
    aggregate(out="outputs/scaled100k",
              ladder=(2000, 5000, 10000, 25000, 50000, 100000),
              n_boot=5000, metrics_path="results/gm4_100k_metrics.json")
    print("=== AGGREGATE_DONE ===", flush=True)
