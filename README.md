# AtrialSpectralBench

**A spectral-fragility benchmark for surgery-triggered atrial fibrillation — entirely in silico.**

AtrialSpectralBench (`asb`) studies where a fibrotic atrium becomes reentry-inducible
after diffuse perioperative stress, using a single closed-form spectral-perturbation
result as an interpretable add-on feature inside an otherwise established
cohort → feature → inducibility-classifier pipeline.

## Status & headline results (honest outcome)

The pre-registered primary endpoint (SFI adds ΔAUC ≥ 0.05 with DeLong p<0.05) is **not met** —
and per `docs/PRE_REGISTRATION.md` §1, reporting this null is an *accepted success of the protocol*.
The full writeup is **`docs/paper/PREPRINT.md`**; every number traces to `results/*.json`.

- **Predictive null**, replicated on two real cohorts (Roney n=62, UW/Boyle n=82) and up to
  **100,000** synthetic excitable networks: SFI adds no practically meaningful signal (ΔAUC ≈ +0.003;
  significant only at N≥10⁴, an effect needing **~4,150 patients** to detect — clinically undetectable).
- **Validity radius (why it fails):** the first-order SFI is valid only while `ρ = ‖ΔL‖/(λ₃−λ₂) = O(1)`;
  real atria have **ρ ≈ 2422** (~10³× past), so the linear biomarker provably cannot work.
- **ρ scaling law:** the spectral gap collapses `∝ N⁻¹·⁰⁵` on the real atrial surface (Weyl −1.0), so ρ
  grows lawfully with resolution — the failure is structural, not incidental.
- **Transfer:** the Fiedler-gradient localizer `|∇φ₂|` locates the instability origin above a spatial
  null across **three** dynamical media (cardiac, FitzHugh–Nagumo, Kuramoto) — though not uniquely
  better than substrate imaging, for a reason made precise.
- **Falsification protocol:** naive methodology would have called SFI a working biomarker; our
  pre-registered guards correctly reject it.

**Honest limitations:** labels are a monodomain simulator (not clinical, not the deferred openCARP
gold standard); see `docs/paper/PREPRINT.md` §4 and the openCARP validation in `ROADMAP.md` §3.1.
Read order: `docs/PRE_REGISTRATION.md` → `docs/SFI_THEORY.md` → `docs/paper/PREPRINT.md` →
`notebooks/lab_notebook.md`. Reproduce via `docs/paper/REPRODUCE.md`.

## The protected novel seed

The load-bearing contribution is the closed-form sensitivity of the atrial algebraic
connectivity (Fiedler value λ₂) to edge uncoupling:

> **dλ₂/dw_ij = (φ_i − φ_j)²**   (φ = Fiedler vector, φ₂)

Because the weighted graph Laplacian factorizes as
`L = Σ_(i,j) w_ij (e_i − e_j)(e_i − e_j)ᵀ`, the first-order sensitivity of any simple
eigenvalue λ with unit eigenvector φ is `dλ/dw_ij = φᵀ(∂L/∂w_ij)φ = (φ_i − φ_j)²`.

Aggregated under a **stochastic, fibrosis-weighted, diffuse perioperative-uncoupling
field** Δw (NOT literal surgical cuts), this defines a per-region **Spectral Fragility
Index (SFI)** — the expected drop in algebraic connectivity per unit perioperative
uncoupling:

```
SFI(R) = E[Δλ₂] ≈ Σ_(i,j)∈R E[Δw_ij] (φ_i − φ_j)²
                   + Σ_{k≠2} (ψ_kᵀ ΔL φ₂)² / (λ₂ − λ_k)
```

Reentry-initiation hotspots are hypothesized where `|∇φ₂|` co-localizes with the
Perron (dominant-activation-mode) eigenvector. A spectral-projector / subspace
generalization keeps SFI well-defined when λ₂ is near-degenerate.

## In-silico, no-patient framing

This project is **100% in silico**: zero patients, zero IRB (only a documented
human-data *exemption* determination for public de-identified geometries). Ground-truth
inducibility labels come from an EP simulator's verdict — **openCARP** when available, and
an in-repo `mock_ep` excitable-media stand-in otherwise.

The `mock_ep` labels are a **development convenience and are never clinical POAF**. Every
public surface that emits or consumes them says so. Claims are scoped to "a cheap spectral
surrogate for a specific EP-simulator's inducibility verdict" — there is no FDA/VICTRE/ISCT
language, and nothing here is validated against real post-operative AF.

## Real data & openCARP are deferred / network-gated

The build environment has **no network access** (zenodo.org / opencarp.org / GitHub are
blocked by policy). Accordingly:

- Real-data loaders (Rodero Zenodo 4506930, Roney Zenodo 5801337, Fibre Atlas Zenodo
  3764917) and the openCARP runner are written as clean, documented interfaces that raise
  an informative `DataUnavailableError` / `OpenCARPUnavailableError` if the data or binary
  is absent. **They are not executed here.**
- The heavy openCARP inducibility sweep and final publication figures are a **deferred**
  phase; until then the synthetic cohort plus `mock_ep` labels exercise the full pipeline
  end to end, so it is one config-run away from real labels.

## Quickstart

```bash
# 1. Create and activate a virtual environment (Python 3.11, CPU-only)
python3 -m venv .venv
source .venv/bin/activate

# 2. Install the package (editable) with dashboard + dev extras
pip install -e ".[dashboard,dev]"

# 3. Run the analytic + property test suite (hard CI gates)
pytest -q

# 4. Run the end-to-end pipeline on the synthetic cohort
asb run --config configs/default.yaml
#   -> outputs/results_report.md, outputs/metrics.json, and the five figures

# 5. Launch the interactive dashboard (needs the [dashboard] extra)
streamlit run src/asb/dashboard/app.py
```

## Determinism & reproducibility

Every stochastic function takes an explicit `rng: np.random.Generator` or `seed: int` —
no global `np.random` state, no wall-clock seeds. All experiment knobs live in a single
YAML config (`configs/default.yaml`), emitted from `asb.config.Config`.

## License

MIT — see [LICENSE](LICENSE).
