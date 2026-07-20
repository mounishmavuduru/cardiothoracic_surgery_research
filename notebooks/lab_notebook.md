# AtrialSpectralBench — Lab Notebook

> **Why this file exists.** This is a dated research notebook. Science-fair and
> peer judges routinely inspect the lab notebook to verify that the work was done
> incrementally, honestly, and reproducibly — that hypotheses were written down
> *before* results, that dead ends and negative results were recorded, and that
> every claimed figure traces back to a seeded, config-driven run. Keep it append-only:
> add new dated entries, never silently rewrite old ones. If an earlier entry was wrong,
> add a later entry that corrects it and say so.

## How to keep this notebook

- **One entry per working session**, newest at the bottom (or top — pick one and be
  consistent). Always start with a real date (`YYYY-MM-DD`).
- Record: what you set out to do, what you actually did, the exact command / config /
  seed, what you observed (numbers, not vibes), what surprised you, and the next step.
- Log **negative and null results** with equal weight — a clean null
  ("SFI ≡ re-encoded fibrosis") is an accepted, reportable outcome for this project.
- Reference the commit hash and the config file for anything reproducible.
- Never invent data or backfill a date. If you did nothing on a day, there is no entry.

---

## Entry template (copy for each new session)

### YYYY-MM-DD — <short title>

- **Goal:** what question or task this session targets.
- **Setup:** branch / commit hash, config file, seed(s), environment notes.
- **Did:** the concrete steps taken (commands run, files changed, experiments launched).
- **Observed:** measured results — metrics, AUC/ΔAUC, plots produced, test pass/fail,
  error messages. Paste the numbers, not a summary.
- **Interpretation:** what the numbers mean; whether they support or falsify the
  pre-registered hypothesis; caveats and threats to validity.
- **Surprises / dead ends:** anything unexpected, including things that did not work.
- **Next:** the single most important next action.

---

## Entries

### YYYY-MM-DD — Phase 0: scaffold & environment

- **Goal:** stand up the repository skeleton so every later phase is one command away.
- **Setup:** branch `claude/third-idea-project-plan-*`; Python 3.11 CPU-only venv at
  `.venv`; config `configs/default.yaml` emitted from `asb.config.Config`.
- **Did:** created the venv, editable-installed the package with `[dashboard,dev]`
  extras, verified `import asb, asb.types, asb.config` and the default-config emit,
  wrote LICENSE / README / Makefile / CI / Dockerfile / this notebook / the
  pre-registration document.
- **Observed:** _fill in — package version, numpy/scipy/sklearn versions, whether the
  backbone import and config emit succeeded._
- **Interpretation:** _fill in._
- **Surprises / dead ends:** _fill in._
- **Next:** implement the spectral engine (Phase 1) and its analytic CI gates
  (`λ_k = 2 − 2cos(kπ/N)` for path/ring/grid).

---

### 2026-07-20 — STEP 1 verify engine + STEP 0 orientation; E0/E1 approach decision

- **Goal:** fresh-container engine verification (pytest + synthetic end-to-end), then
  scout the real-data + solver landscape to fix the E0/E1 approach before
  pre-registration and before generating/looking at any inducibility label.
- **Setup:** branch `claude/third-idea-project-plan-rmerjr`, commit at HEAD `deaa982`
  (identical to the `...c2c2bc` branch). Python 3.11.15 CPU-only venv; editable install
  of `[dashboard,dev]`. Config `configs/default.yaml` (seed 0). Network now open to
  zenodo.org / opencarp.org.
- **Did:**
  1. Read EXPERIMENTAL_PLAN, IDEA3_PROJECT_PLAN, PRE_REGISTRATION, BUILD_SPEC and all
     `src/asb/` compute modules. Engine matches the frozen contract.
  2. `pytest -q` → **89 passed in 15.35 s**.
  3. `asb run --config configs/default.yaml` → outputs/ populated (metrics.json, report,
     5 PNGs). **Synthetic plumbing baseline (mock_ep labels — NOT a result):** 48
     subjects / 8 shape families, 16 inducible (33.3 %, both classes present).
     fibrosis: grouped AUC 0.859 / naive 0.861. fibrosis+SFI: grouped 0.920 / naive
     0.910. DeLong ΔAUC = **+0.0605**, z = 1.94, **p = 0.053**. Colocalization observed
     0.285 vs null-5th-pct 0.146 → **not passed** (p = 0.20). These come from the
     `mock_ep` eikonal surrogate and are a PLUMBING CHECK ONLY.
  4. Scouted the three Zenodo archives + solver install:
     - **Roney LA virtual cohort 5801337** — 100 VTK meshes, 3.51 GB, **no separate
       label files**. Downloaded the smallest (`Mesh_83815278.vtk`, 21.9 MB, 97 161
       pts / 193 438 tris). Embedded arrays: **`UAC1`,`UAC2`** (∈[0,1] universal atrial
       coords), **`IIR`** (Image-Intensity-Ratio, LGE fibrosis proxy: min 0.41 / max
       1.59 / mean 1.08; IIR>1.2 → 19.4 % fibrotic), **`fiber_endo`,`fiber_epi`**
       (per-cell fibre vectors). Coordinates are in microns. This is REAL patient-derived
       LA anatomy with fibrosis + UAC + fibres — a large upgrade over synthetic icospheres.
     - **Rodero 4506930** — 27.5 GB (won't fit the 29 GB free disk); whole-heart, not
       LA-with-fields. Skipped in favour of Roney LA meshes.
     - **atrialmtk 10139306** — single 872 MB zip (toolkit; bundles openCARP/meshtool).
     - **openCARP** — not available via apt / conda / pip; no binary on PATH. A source
       build (PETSc etc.) on 4 CPUs is out of scope for this ephemeral CPU session, and
       the user explicitly DEFERS the full openCARP sweep to Claude Science.
- **Observed / decision (logged before any label is generated):**
  - The published Roney inducibility/reentry labels are **not** in the archive as files,
    so "reuse published Roney labels" is not directly possible from 5801337. The mesh
    arrays give anatomy+substrate only.
  - **E1 ground-truth labeler for this session = a genuine CPU monodomain
    Mitchell–Schaeffer reaction–diffusion solve on the (coarsened) real Roney meshes**,
    verified (planar CV, APD restitution, spiral induction, mesh/Δt convergence) and
    calibrated to literature CV/ERP. This is exactly the "descoped openCARP: monodomain +
    phenomenological Mitchell–Schaeffer on coarsened meshes" of plan §8.4 — genuine
    nonlinear excitable-media reentry, driven by wavelength/source–sink/fibrosis and
    **independent of λ₂** (so not circular with SFI). It is labelled `source='monodomain_ms'`
    and is NEVER called openCARP or clinical POAF. openCARP itself stays the wired,
    one-command **deferred** target (`labels/opencarp.py`) for the Claude-Science full
    sweep; I will attempt an opportunistic openCARP cross-check if it proves installable.
  - **E0** in this session = calculation-verification I can do on CPU: solver
    verification above + λ₂/label mesh-&-Δt convergence (ASME V&V-40). The openCARP
    Niederer reproduction is wired but deferred (solver absent).
- **Interpretation:** engine is green and faithful to the contract; the synthetic ΔAUC
  hint (+0.06 at p≈0.05) is encouraging plumbing but scientifically inert (mock labels).
  The real science now runs on real Roney anatomy with monodomain-MS labels.
- **Surprises / dead ends:** Rodero too big for disk; openCARP not trivially installable;
  Roney labels not archived as files — all resolved by the decision above.
- **Next:** freeze the pre-registration (label-source reality made explicit + Δw
  calibration frozen) and git-tag it BEFORE generating any label; then E1 (real loader +
  coarsening + monodomain-MS labeler + inducible-fraction gate) → GM1 head-to-head.

---

### 2026-07-20 — E1 build: real loader, monodomain-MS solver, and inducible-fraction calibration (GATE PASSED)

- **Goal:** build the real-anatomy substrate + a genuine nonlinear inducibility labeler
  and calibrate the inducible fraction into the pre-registered 10–40 % band.
- **Setup:** branch `claude/third-idea-project-plan-rmerjr`; new modules
  `src/asb/substrate/roney.py` (VTK loader for UAC/IIR/fibre + field-preserving
  vertex-clustering coarsening) and `src/asb/labels/monodomain.py` (monodomain
  Mitchell–Schaeffer solver, implicit-diffusion cotangent Laplacian). 30 smallest Roney
  meshes (Zenodo 5801337) downloaded, coarsened to ~2 000 nodes.
- **Did:**
  1. Verified the solver physics (E0 calculation-verification, CPU): single-cell **APD90
     = 218 ms healthy / 121 ms fibrotic** (τ_close 110 → 55 ms); planar-wave **CV = 0.43–
     1.24 m/s** across d0, i.e. squarely in the 0.3–1.2 m/s human-LA band (Heida 2021).
  2. Fixed an explicit-diffusion CFL blow-up (V→1.5 on coarse meshes) by switching to
     **implicit backward-Euler diffusion** (factorize (I−dt·A) once); peak V now ≤ 0.95.
  3. Found the half-field cross-field inducer anchors *geometric* reentry independent of
     substrate (induces in healthy too) → **dropped from the label battery**. Adopted
     **burst pacing** (substrate-discriminating: rapid pacing between healthy and fibrotic
     ERP breaks waves only where refractoriness is short/heterogeneous).
  4. Defined inducibility as **self-sustained supra-threshold activity ≥ reentry_min_ms
     after pacing ends** (a normal paced response transits + repolarizes well inside it).
  5. Confirmed **monotonicity**: on one mesh, scaling fibrosis 0→0.17→0.34→0.44 gave
     sustained reentry 224→246→599→1000 ms.
- **Observed (24–30-mesh parallel calibration):**
  - First pass (d0 0.12, 2 CL × 2 sites, min 600 ms): fraction **0.75** — too high (AF
    cohort is high-fibrosis), Spearman(fib,sustained)=0.67.
  - **Frozen pass** (d0 **0.20** ≈ healthy CV 0.87 m/s along-fibre, burst CL **150**,
    **2 sites**, n_burst 6, reentry_min **650 ms**, coarsen 2000): **inducible fraction
    = 0.33 (10/30)** — IN the 10–40 % gate, near ~26 %. **Spearman(fibrosis_burden,
    sustained) = 0.748.** sustained_ms is cleanly **bimodal** (non-inducible ≤ 584 ms;
    inducible ≥ 854 ms) so the label is robust to the exact threshold (0.33 across
    600–850 ms).
  - Headroom for SFI is real: highest-burden mesh (fib 0.77) is **non-inducible** while a
    0.43-burden mesh is inducible — inducibility depends on fibrosis *pattern/geometry*,
    not burden alone.
- **Interpretation:** E1 gate **PASSED** — a genuine nonlinear monodomain-MS ground truth,
  physiologically calibrated, substrate-driven, in the target prevalence band. Ready to
  freeze the pre-registration and run GM1. Label framing: the primary endpoint is
  "inducibility classification"; the acute perioperative Δw stress is embedded in the SFI
  feature (fibrosis-weighted expected λ₂ drop), the label is the simulator's inducibility
  verdict on the fibrotic substrate. Post-stress relabelling (apply the Δw field, then
  simulate) is a wired extension.
- **Surprises / dead ends:** explicit-diffusion CFL blow-up on coarse meshes; the
  half-field cross-field protocol's substrate-independent geometric reentry. Both resolved.
- **Next:** freeze pre-registration (frozen monodomain protocol + Δw = 0.36 from
  f_CV = 0.20 → f_D = 1−0.8² and MS/IIR constants) + git-tag; scale the cohort; run GM1.

---

### 2026-07-20 — GM1 ★ PRIMARY head-to-head — NULL on the pre-registered endpoint (reportable), with the validity-radius mechanism

- **Goal:** the pre-registered primary test — does per-region SFI add grouped ΔAUC ≥ 0.05
  with DeLong p < 0.05 over the real competitor set, on real Roney anatomy + monodomain-MS labels?
- **Setup:** commit at `prereg-frozen-20260720` (5f5a616) frozen BEFORE this analysis.
  N = **62** Roney patients, coarsened ~2000 nodes; **20 inducible (0.32)**, both classes;
  patient-held-out GroupKFold (each patient its own group); classifiers LR (nested
  GroupKFold C-selection) + GBT; DeLong + **10 000-resample** bootstrap. An independent
  code audit (subagent) confirmed no leakage / DeLong sign-flip / bootstrap error /
  spec_-column contamination / λ₂-circularity in the label path → the number is trustworthy.
- **Observed (grouped = patient-held-out; naive alongside):**

  | comparison | clf | grouped base | grouped +SFI | grouped **ΔAUC** | **DeLong p** | naive base | naive +SFI | boot ΔAUC [95% CI] |
  |---|---|---|---|---|---|---|---|---|
  | competitors (6) +SFI | LR | 0.832 | 0.802 | **−0.030** | 0.643 | 0.849 | 0.774 | −0.031 [−0.167, +0.086] |
  | competitors (6) +SFI | GBT | 0.820 | 0.856 | **+0.036** | 0.257 | 0.804 | 0.836 | +0.035 [−0.029, +0.096] |
  | fibrosis-het (3) +SFI | LR | 0.781 | 0.798 | +0.017 | 0.594 | 0.818 | 0.777 | +0.017 [−0.045, +0.082] |
  | fibrosis-het (3) +SFI | GBT | 0.792 | 0.861 | **+0.069** | 0.064 | 0.776 | 0.837 | +0.069 [−0.004, +0.144] |

- **Interpretation — NULL on the primary endpoint (an accepted, pre-registered outcome).**
  Against the full 6-competitor set (which already includes λ₂-alone) SFI adds nothing
  significant: −0.030 (LR) / +0.036 (GBT, p=0.26). The strongest signal is GBT adding SFI
  to the fibrosis-heterogeneity-only baseline (ΔAUC +0.069, bootstrap one-sided p(Δ≤0)=0.030)
  but DeLong p=0.064 misses <0.05 and it does not clear the *full* competitor set. This is
  the pre-registered "SFI ≡ re-encoded fibrosis" null. **Not tuned; reported as-is.**
- **The mechanism (this is the GM3 "honest core", arriving early and quantified):** the
  λ₂–λ₃ **spectral gap is minuscule on real LA anatomy** (cohort gap min 2.2e-4, median
  2.6e-3, max 5.4e-3). The frozen perioperative stress has ‖ΔL‖ ≈ 3.2–3.8, so the
  first-order validity ratio **‖ΔL‖/(λ₃−λ₂) has median ≈ 2422 and is > the safety bound on
  0/12 meshes checked** (`results/validity_radius_scan.txt`). The single-Fiedler-vector SFI
  is therefore applied **thousands of times outside its perturbative regime** — λ₂ is
  near-degenerate everywhere, the Fiedler vector is ill-conditioned, so the linear SFI
  feature is noisy and collapses to fibrosis. This is *exactly* the pre-registered pivot:
  "when/why does spectral fragility collapse to fibrosis?" — answered with a hard number.
- **Surprises:** the gap is even smaller than expected on real anatomy → the validity
  radius is violated universally, not marginally. GBT (nonlinear) extracts a little more
  from SFI than LR, hinting SFI's residual value is in nonlinear interactions.
- **Next (paused for PI decision):** the pre-registered fallback makes **GM3 the lead** —
  formalize the validity radius (Weyl / Davis–Kahan) and test whether the **subspace
  (projector) SFI**, which stays well-defined under near-degeneracy, recovers the signal the
  single-vector SFI loses. Also available: exact-Δλ₂ / MC SFI features (valid at large ‖ΔL‖).
  Artifacts: `results/gm1_metrics.json`, `results/gm1_report.md`, `results/validity_radius_scan.txt`.
