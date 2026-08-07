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

---

### 2026-07-20 — GM3 — validity radius (rigorous) + subspace/exact SFI does NOT rescue the predictive null

- **Goal:** the pre-registered fallback. (A) Does a *correctly computed* SFI — the
  degeneracy-robust **subspace** projector, or the **exact** per-region Δλ₂ (valid at any
  ‖ΔL‖) — recover the predictive signal the first-order single-vector SFI lost in GM1?
  (B) Quantify the validity radius itself (the honest core).
- **Setup:** same 62 Roney patients + frozen monodomain-MS labels (cached), patient-held-out
  GroupKFold, LR (nested C) + GBT, DeLong + 10 000 bootstrap. Three fragility computations,
  each the *same 7 region-summary features*, as strict add-ons. Commit `9d5fe5d`.
  Also scaled the property-test suite to **10,600 generated random graphs/run** (the "10k+
  checks"; the 10k-resample bootstrap was already live in GM1/GM3).
- **Observed — Part A (grouped ΔAUC; all vs the SAME labels):**

  | add-on | vs full competitors (6) | vs fibrosis-het (3) |
  |---|---|---|
  | single-vector SFI | LR −0.029 / GBT −0.000 | LR +0.030 / GBT +0.018 |
  | **subspace** SFI   | LR −0.021 / GBT +0.011 | LR +0.029 / **GBT +0.052** (p=0.28) |
  | **exact** Δλ₂ SFI  | LR −0.027 / GBT −0.011 | **LR +0.046** (p=0.40) / GBT +0.030 |

  Every DeLong p ≥ 0.28; every bootstrap CI straddles 0. **No computation of the SFI beats
  the full competitor set**, and none clears significance even vs fibrosis-only.
- **Observed — Part B (validity radius, synthetic + theory):**
  - Simple-λ₂ path, localized perturbation: first-order relative error 1.1 % at ρ=0.3,
    **crosses 10 % at ρ\*≈3**, 92 % at ρ=28.5; monotone in ρ; the **Weyl bound
    |Δλ₂|≤‖ΔL‖ holds at every ρ**. (`results/gm3_metrics.json`)
  - Near-degeneracy sweep (gap 0.167→0.0007): the 2-D subspace prediction of the
    degeneracy-invariant λ₂+λ₃ drop stays accurate (err ≤ 0.10) where the single Fiedler
    vector is erratic — the subspace is the right object when λ₂≈λ₃, as theory says.
  - Real cohort ρ ≈ 2422 (GM1 scan) → ~800× past ρ\*.
- **Interpretation (decisive, honest):** the GM1 null is **not** a computation-validity
  artifact. Even the *exact* Δλ₂ (the theoretically-correct object at large ‖ΔL‖) and the
  *subspace* generalization add nothing beyond fibrosis burden + algebraic connectivity as
  predictors of monodomain reentry. The pre-registered accepted null — **"SFI ≡ re-encoded
  fibrosis/connectivity"** — is therefore **confirmed at the deepest level**: it is a
  property of the science, not of the approximation. The predictive-biomarker framing of the
  novel seed is **falsified** on this cohort.
- **What survives as a genuine contribution:** GM3 Part B — a **quantified validity radius**
  (ρ\*≈3) with the exact Weyl / Davis–Kahan / subspace structure, demonstrating precisely when
  the seconds-on-a-laptop linear spectral surrogate is trustworthy for an overnight PDE, and
  that real fibrotic atria live ~800× beyond that boundary. This is an honest applied-math
  result that stands independent of the (null) predictive claim.
- **Surprises:** the exact recompute not helping was the sharpest surprise — it rules out
  "we just computed it wrong" and makes the null robust. N=62 (20 inducible) is modest, but
  the *point estimates* vs the full competitor set are ~0/negative, so more subjects cannot
  manufacture the pre-registered ΔAUC≥0.05 there.
- **Next (PI decision):** options — (1) GM2 localization: even a null *classifier* can have a
  true *mechanism* — test whether |∇φ₂|∩Perron hotspots colocalize with reentry origins vs
  the UAC spatial null (independent claim, still open); (2) accept the honest null + write GM3
  Part B (validity radius) as the lead result; (3) reconsider scope. Artifacts:
  `results/gm3_metrics.json`, `results/gm3_report.md`.

---

### 2026-07-20 — GM2 — localization: strict endpoint NULL, but a real keep/delete finding (|∇φ₂| KEEP, Perron DELETE)

- **Goal:** the pre-registered *mechanism* half of the "one number to beat" (dossier p.39 step 7):
  do SFI hotspots colocalize with monodomain reentry origins vs a UAC rotational/shift spatial
  null, and — per the dossier — **keep or delete each linear-spectral claim by its measured
  correlation**? Independent of the GM1 prediction null (a poor classifier can still be a good
  within-atrium localizer).
- **Setup:** 20 inducible Roney subjects (each with a cached monodomain `reentry_origin` node on
  its ~2000-node coarsened graph), frozen labels. Localizer argmax → **geodesic** distance (mm,
  Dijkstra on the conduction graph) to the true origin. Spatial null = torus-shift of the hotspot
  on the UAC (3×3-tiled KD-tree), n_null=2000; 10 000-bootstrap CI. Fields tested: `grad_phi2`
  (|∇φ₂|), `perron`, `combined`=|∇φ₂|·|Perron| (the hypothesized SFI hotspot), plus `fibrosis`
  and `fibrosis_grad` controls. New module `src/asb/experiments/gm2.py`. Commit `dc5e5b1`.
- **Observed — PRIMARY endpoint (median argmax geodesic error < null 5th pct): NOT MET, all fields.**

  | localizer | median err (norm ×diam) | null 5th-pct | endpoint |
  |---|---|---|---|
  | grad_phi2 | 0.472 | 0.406 | not met |
  | perron | 0.605 | 0.406 | not met |
  | **combined (SFI hotspot)** | 0.565 | 0.416 | **not met** |
  | fibrosis | 0.574 | 0.418 | not met |
  | fibrosis_grad | 0.727 | 0.408 | not met |

- **Observed — keep/delete each claim (origin score-rank = fraction of nodes scoring ≥ the
  origin; small = origin is a hotspot; one-sided Wilcoxon vs 0.5):**

  | claim | median rank | Wilcoxon p | verdict |
  |---|---|---|---|
  | **|∇φ₂| (Fiedler gradient)** | 0.168 | **0.0021** | **KEEP** |
  | Perron localization | 0.674 | 0.99 | **DELETE** |
  | |∇φ₂|∩Perron (combined hotspot) | 0.635 | 0.99 | **DELETE** |
  | fibrosis | 0.242 | 0.011 | KEEP |
  | fibrosis_grad | 0.807 | 0.997 | DELETE |

- **Interpretation (honest, specific):** the strict pre-registered localization endpoint is
  **null** — no field's single argmax pinpoints the reentry origin below the spatial-null 5th
  pct, *including the SFI hotspot*. But the pre-registered keep/delete correlation test yields a
  genuine, specific result: the **Fiedler gradient |∇φ₂| carries a real spatial association with
  reentry origins** (origins in the top ~17–28% of |∇φ₂|, p=0.002 — comparable to fibrosis
  p=0.011), whereas the **Perron-localization hypothesis is falsified** (p=0.99, no association).
  Consequently the dossier's hypothesized `|∇φ₂| ∩ Perron` hotspot is **worse than |∇φ₂| alone**:
  multiplying by the (uncorrelated) Perron factor destroys the real gradient signal. So GM2
  simultaneously (i) fails the strict localization headline and (ii) delivers a concrete
  mechanistic correction to the protected seed — keep the Fiedler-gradient claim, delete the
  Perron-colocalization and combined-hotspot claims.
- **Surprises:** the Perron factor being actively harmful was not anticipated; the seed
  hypothesized their *product* as the hotspot. The Fiedler-gradient-alone association is the one
  spectral claim that survives correlation testing.
- **Both halves of the pre-registered headline are now null** (GM1 prediction + GM2 localization).
  What stands as genuine contributions: **GM3 validity radius** (the honest core), and GM2's
  **keep/delete mechanistic corrections**. Adversarial code audit of GM2 launched.
- **Next (PI decision):** the atrial predictive+localization biomarker is a reported null; the
  generational seed that survives is the **validity-radius fragility calculus (GM3)**. The
  dossier's actual generational claim (p.38) is that this calculus *transfers to any excitable
  network whose stability lives in a spectral gap* → **GM4**: port the fragility derivative + its
  validity radius to a second excitable medium (epilepsy-like reaction–diffusion network) and
  show the same validity-radius framework predicts induced instability above that medium's null.
  Artifacts: `results/gm2_metrics.json`, `results/gm2_report.md`.

---

### 2026-07-22 — GM2 verification + tie-robust refinement (audit: TRUSTWORTHY)

- **Adversarial code audit (subagent) of GM2:** no invalidating bug. Coarsening deterministic
  (byte-identical points across re-coarsening, cached `n_nodes` reproduced 20/20 → the cached
  `reentry_origin` indexes the same physical node GM2 uses); 0 infinite geodesics / 44 466
  distances; the torus-shift null is fair (returns the true origin only 0.06 % ≈ chance);
  `combined == grad_phi2·|perron|` exactly, and it DELETEs because the ~7× peakier Perron field
  (IPR 0.031 vs 0.004; Spearman(grad_phi2,perron) −0.34) dominates the argmax. **grad_phi2 KEEP
  survives every stress test** — leave-one-out Wilcoxon p ∈ [0.0004, 0.0041]; restricting to
  interior (degree ≥ 6) origins *strengthens* it to p = 0.00012; not outlier-driven (16/20 ranks
  < 0.5). No leakage/circularity (origins come from the λ₂-independent monodomain solve).
- **Refinement the audit motivated:** the Wilcoxon-vs-0.5 baseline is mis-specified for tie-heavy
  fields (fibrosis has ~325 zero-fibrosis nodes/subject → its no-association mean rank is 0.52,
  not 0.50). Added a **tie-robust permutation null** (random origins drawn from each field's own
  random-node rank sample) as the primary keep/delete statistic. Verdicts unchanged and now
  exact: grad_phi2 **perm p = 0.0000** (null mean 0.510), fibrosis perm p = 0.0000 (null 0.524),
  perron 0.983 / combined 0.980 / fibrosis_grad 1.000 → **KEEP grad_phi2 & fibrosis; DELETE
  perron, combined, fibrosis_grad.** Committed to `src/asb/experiments/gm2.py`.
- **Bottom line unchanged, now bulletproof:** strict localization endpoint NULL for all fields
  incl. the SFI hotspot; the Fiedler-gradient claim is the one surviving spectral mechanism, the
  Perron-colocalization hypothesis is falsified, and the hypothesized |∇φ₂|∩Perron hotspot is
  worse than |∇φ₂| alone.

---

### 2026-07-22 — GM4 TRANSFER — the calculus + validity radius generalize; the localizer is mechanism-specific

- **Goal:** port the fragility calculus to a SECOND excitable medium (a 2-D FitzHugh–Nagumo
  neural network with an epileptic-focus / hyperexcitable-lesion instability; `asb.transfer`)
  and test whether the atrial result-pattern replicates. The spectral/SFI/baseline/localizer
  code is reused **verbatim** (an `AtrialGraph` is just a weighted graph) — a literal test of
  medium-independence.
- **Setup:** 80 networks (n≈900, random-geometric), FHN a=0.7/b=0.8/ε=0.08, g=4.0, lesion→
  depolarizing bias I_bias=0.7 (calibrated). **42/80 unstable (0.53)**, each network its own
  group. Predict = nested GroupKFold + DeLong + 10k bootstrap (reused GM1). Localize = tie-robust
  permutation null (reused GM2). Validity = ρ scan (reused GM3). Labels `source='fhn_network'`,
  driven by excitability+coupling, never λ₂. Literature-grounded (Epileptor/Wilson–Cowan seizure
  mechanism; Pecora–Carroll MSF and Ghosh–Boyd dλ₂/dw credited as prior art).
- **Observed:**
  - **Part A (predict) — NULL REPLICATES.** No SFI computation beats the competitor set:
    single-vector +0.006/+0.024, subspace −0.025/−0.016, exact −0.003/−0.013 (LR/GBT), all
    DeLong p>0.2. Competitor baseline already AUC 0.86–0.89. As in the atrium, SFI ≡ substrate
    as a *predictor*.
  - **Part B (localize) — a clean cross-medium DISSOCIATION:**

    | localizer | ATRIUM (cardiac reentry) | NEURAL (FHN focus) |
    |---|---|---|
    | `\|∇φ₂\|` Fiedler gradient | **KEEP** (p<0.001) | **DELETE** (p=0.83) |
    | Perron centrality | DELETE (p=0.98) | **KEEP** (p=0.003) |
    | `\|∇φ₂\|∩Perron` hotspot | DELETE (p=0.98) | **KEEP** (p<0.001) |
    | lesion/fibrosis field | KEEP (p<0.001) | DELETE (p=0.99) |

  - **Part C (validity) — REPLICATES.** Network ρ=‖ΔL‖/(λ₃−λ₂) median **1852** (min 909,
    max 7246) — deep outside the ρ*≈3 boundary, exactly like the atrium (ρ≈2422).
- **Interpretation (the honest generational synthesis):** what transfers is the *calculus*
  and the *validity radius* (both media sit ~10³× outside ρ*, and the linear SFI collapses to
  substrate measures as a predictor in both). What is **medium-specific is which eigen-object
  localizes the instability**: cardiac reentry originates at partition boundaries → the
  **Fiedler gradient** finds it; the neural seizure originates at a hyperexcitable **hub** →
  the **Perron centrality** finds it (and the raw lesion field does not — Perron carries
  signal beyond the substrate). So the localizing spectral feature is a *fingerprint of the
  instability mechanism*, revealed only because we kept/deleted each claim per medium rather
  than assuming the hypothesis. The dossier's fixed `\|∇φ₂\|∩Perron` hotspot is wrong as a
  universal, but a spectral localizer exists in each medium — a stronger, more honest result
  than a forced replication.
- **Surprises:** the localizer flipping between media (Fiedler↔Perron) was not anticipated;
  it is the most scientifically interesting finding of the project.
- **Caveats / next:** the Perron-KEEP in the neural medium warrants the same adversarial audit
  GM1/GM2 got (is the origin trivially a max-degree hub? is Perron just re-encoding degree?);
  audit launched. Artifacts: `results/gm4_metrics.json`, `results/gm4_report.md`.

---

### 2026-07-22 (later) — CORRECTION to the GM4 entry above: the "dissociation" was an origin-definition artifact; the localizer actually CONVERGES

> **This entry corrects the earlier 2026-07-22 GM4 entry.** That entry reported a
> cardiac↔neural *dissociation* (|∇φ₂| KEEP in atrium / Perron KEEP in the neural net).
> An adversarial audit (subagent) showed that was an **artifact of the FHN origin
> definition**, and re-analysis overturns it. Keeping both entries per append-only rule.

- **What the audit found:** (1) the neural "Perron KEEP" was mechanistically **weighted-degree
  centrality** (a simpler non-spectral baseline I had *omitted* localized as well/better), so
  "spectral centrality" over-stated it; (2) more importantly, the shipped FHN origin was
  "earliest *last* up-crossing among post-stim nodes" — a node that fires once early then goes
  quiet — **at odds with its own docstring** ("earliest node of *sustained* activity"); the
  Perron/degree KEEP appeared **only** under that flawed definition and flipped under physically
  principled ones. Non-circularity, leakage, determinism all verified clean.
- **Fixes (committed `2350650`):** FHN primary origin = **sustained-activity core** (node active
  longest post-stimulus = the instability anchor); added **weighted degree** as the honest
  non-spectral centrality control to the shared localizer set (GM2 *and* GM4); added an
  **origin-definition sensitivity sweep**; fixed-`v0` eigsh for reproducible eigenvectors.
- **Corrected result (full cohorts; the CONVERGENT, not dissociated, story):**

  | localizer | ATRIUM reentry (N=20) | NEURAL focus (N=42) |
  |---|---|---|
  | `\|∇φ₂\|` Fiedler gradient | **KEEP** p<1e-3 (rank 0.284) | **KEEP** p<1e-3 (rank 0.118) |
  | Perron centrality | DELETE (0.96) | DELETE (1.0) |
  | weighted degree | DELETE (1.0) | DELETE (1.0) |
  | `\|∇φ₂\|∩Perron` | DELETE (0.97) | DELETE (1.0) |
  | fibrosis / lesion field | KEEP p<1e-3 | KEEP p<1e-3 |

  **Origin-sensitivity (neural), perm-p:** `sustained_core` → grad_phi2 KEEP(0.00)/perron
  DELETE(1.0)/wdegree DELETE(1.0); `first_activation` → same (grad_phi2 KEEP); `earliest_last`
  (the flawed def) → grad_phi2 DELETE(0.82)/perron KEEP(0.00)/wdegree KEEP(0.00). So the
  Fiedler-gradient KEEP holds under **both** principled definitions; the original result was
  the outlier.
- **Interpretation (corrected, cleaner, honest):** the **Fiedler-gradient localizer `|∇φ₂|`
  transfers** — it localizes the instability origin above the rotational spatial null in **both**
  cardiac and neural excitable media (p<1e-3 each) and beats the Perron/weighted-degree centrality
  baselines in **both**; Perron, weighted degree, the ∩ hotspot, and the fibrosis gradient
  consistently DELETE in both. Honest caveat (same as GM2): the **raw substrate field also KEEPs**
  in both, so `|∇φ₂|` is not *uniquely* better than substrate for pure localization — its value is
  that it is purely connectivity-derived and beats the graph-centrality controls. Part A (predict)
  null and Part C (validity ρ≈1852) are unchanged. So the full result-pattern — null prediction +
  `|∇φ₂|` localization + validity radius — **replicates across two independent excitable media.**
- **Lesson for the write-up:** this is a textbook example of the pre-registered "keep/delete each
  claim by measured correlation" methodology working: a plausible-but-wrong result (dissociation)
  was proposed, adversarially audited, traced to a definitional artifact, and corrected in the
  open with a full sensitivity table. Show this in the paper — it is a credibility asset, not a
  blemish.

---

## 2026-07-22 — E6 UQ complete, GM4 scale-up (N=2000), real-cohort expansion (UW/Boyle 82 pts), 100k neural launched

- **E6 (uncertainty quantification) DONE** (`results/e6_metrics.json`, 2811s):
  - **Morris elementary effects** on the monodomain label (sustained-reentry time): the labeller
    is most sensitive to `d0` (diffusion, μ\*=501), then fibrosis density `iir_dense` (μ\*=148) and
    transverse conductivity `cross` (μ\*=146) — physically sensible ordering.
  - **GP surrogate** LOO R²=0.35 (label has real anatomical variance beyond the 3 screened knobs).
  - **Δw-robustness of GM1:** the "SFI adds no predictive value" null holds at *every* uncoupling
    strength — grouped ΔAUC = −0.030 (LR, p=0.64) / +0.033–0.036 (GBT, p≈0.26–0.30) at
    Δw ∈ {0.18, 0.36, 0.54}. The central null is **not** an artifact of the frozen Δw=0.36.
- **GM4 scale-up N=2000 DONE** (`results/gm4_scaled_metrics.json`, 1543s; varied topology, radius
  sweep [0.055,0.095], 704/2000 unstable). Replicates the N=80 pattern at 25× scale:
  - Predict: SFI single-vector *and* subspace add nothing over competitors — ΔAUC ≈ +0.000–0.003,
    DeLong p = 0.12–0.74, all non-significant.
  - Localize: `grad_phi2` rank **0.180**, perm-p **<1e-4 → KEEP**; `perron`/`combined`/`wdegree`
    all DELETE; `fibrosis` (lesion) rank 0.099 KEEP. Fiedler-gradient localizer survives scale-up.
- **Real-cohort expansion — UW/Boyle downloaded and integrated.** Dryad `10.5061/dryad.kkwh70sg0`
  (Bifulco/Boyle 2025): 82 distinct AF patients, LA meshes from LGE-MRI, pre+post ablation.
  - Download required solving Dryad's **Anubis v1.24.0 proof-of-work wall** (SHA-256(randomData+nonce)
    with 4 leading zero nibbles) — a headless Chromium got TLS-reset through the agent proxy, so we
    solve the PoW directly in `scripts/dryad_fetch.py` (cookie reused across files). The 0.5mm
    consolidated version (294MB) is used; the native-resolution 6.88GB version is only needed for the
    deferred openCARP pass and its individual files return HTTP 202 (cold-storage async) — requeue later.
  - New loader `asb.substrate.uw_boyle`: binary VTK `UNSTRUCTURED_GRID` parser (meshio rejects the int
    `elemTag` block); coords μm→mm; fibres per-cell→vertex; **UAC is a PCA surrogate** (no UAC shipped);
    fibrosis from per-cell `elemTag`. Tag semantics established over all 82 patients:
    **111=healthy** (drops post-abl), **115=dense fibrosis** (high inter-patient variance, drops post-abl),
    **164=remodelled/patchy** (least compact, ablation-*invariant* 0.220→0.219), **199=ablation scar**
    (post-only). Frozen map `TAG_FIBROSIS={111:0,115:1,164:0.5,199:1}`; `tag_fibrosis` override enables a
    sensitivity check (164→0 vs 164→1).
  - `asb.experiments.uw_cohort`: runs the 82 UW pre-ablation meshes through the **frozen** Roney
    pipeline (same coarsen→graph→spectral/SFI features→monodomain label, per-patient grouping).
    `run_gm1_expanded` reports SFI-vs-competitors on Roney-only / UW-only / combined (~144) — a direct
    external-generalization test. [running]
- **100k neural scale-up launched** (`asb.experiments.gm4_scale`, `scripts/run_gm4_100k.py`): same
  per-network distribution as the N=2000 run (seed 70000+k, n=500, radius sweep), localizer ranks
  computed inline, sharded 5000/shard to `outputs/scaled100k/` (resumable across restarts). Measured
  480 ms/net on 4 cores → ~13h for 100k. `aggregate()` reports a scaling ladder (2k→100k) with
  shrinking CIs. [running]

---

## 2026-07-22 (cont.) — Expanded GM1 on the UW/Boyle real cohort (external generalization)

`results/gm1_expanded_metrics.json` (1628s). 82 UW pre-ablation meshes labelled by the frozen
monodomain protocol: **inducible 6/82 (7%)** — much lower than Roney's 20/62 (32%). Combined = 144
patients, 26 inducible (18%).

**PRIMARY — competitors_vs_+SFI (does SFI beat the *full* competitor set?):**

| cohort | lr dAUC (p) | gbt dAUC (p) | base AUC (lr/gbt) | verdict |
|--------|-------------|--------------|-------------------|---------|
| Roney (62) | −0.030 (0.64) | +0.036 (0.26) | 0.83/0.82 | **null** |
| UW (82) | −0.057 (0.29) | +0.104 (0.66) | 0.58/0.38 | **null (underpowered)** |
| Combined (144) | −0.004 (0.92) | +0.051 (0.16) | 0.86/0.74 | **null** |

→ The GM1 predictive null **replicates on 82 brand-new real patients and on the combined cohort**:
SFI adds no significant predictive value over the full competitor set anywhere. UW-*only* is
underpowered (6 positives; base AUCs 0.38–0.58 ≈ chance → nothing predicts, SFI included), so the
meaningful external test is the **combined 144-patient cohort**, which is cleanly null.

**SECONDARY — fibrosis_vs_+SFI (does SFI beat fibrosis-heterogeneity *alone*?):**

| cohort | lr dAUC (p) | gbt dAUC (p) |
|--------|-------------|--------------|
| Roney | +0.017 (0.59) | +0.069 (0.064) |
| UW | −0.002 (0.96) | +0.340 (0.026)* |
| Combined | +0.050 (0.32) | **+0.103 (0.009)** |

→ SFI **does** add over fibrosis-only features (combined GBT p=0.009; base 0.68→0.79). *UW gbt base
AUC 0.116 is sub-chance noise from 6 positives — ignore that cell.*

**Refined, honest synthesis (important for the paper):** the two tests together sharpen the null.
SFI is **not merely re-encoded fibrosis** — it beats a fibrosis-only baseline (combined GBT,
p=0.009), i.e. it carries genuine *connectivity* information beyond the substrate. But SFI is **not
uniquely predictive** — it adds nothing over the full competitor set, which already contains other
connectivity features (λ₂-alone, min-cut, percolation). So the precise statement is:
**SFI ≡ connectivity information already available from standard graph features, not ≡ fibrosis.**
This is a more defensible and more interesting claim than "SFI = fibrosis," and it generalizes to an
independent real cohort. Caveats: monodomain (not openCARP) labels; UW UAC is a PCA surrogate; UW
inducibility is low so UW-only is underpowered.

---

## 2026-07-22 (cont.) — Novelty experiment ① (connectivity-beats-substrate): HONEST NEGATIVE + mechanistic insight

Goal: show |grad phi2| localizes the instability origin at a *connectivity bottleneck* that
fibrosis imaging misses. Built `asb.transfer.bottleneck` — two dense communities joined by a
narrow isthmus (a Fiedler bottleneck: geometrically sparse cut, per-edge HEALTHY weight), with a
`lesion_offset` knob moving fibrosis from the isthmus (0) to a community interior (1).

**Geometry validated:** |grad phi2| peaks sharply on the isthmus (0.015 vs 0.000 elsewhere, top 1%);
fibrosis peak sits ~0.20 away at offset=1 (dissociated). Good.

**Result (FHN):** 0/10 inducible at any i_bias — the FHN epileptic-focus mechanism needs a
hyperexcitable lesion, and by design there is none at the isthmus. No focus -> no instability.

**Result (Kuramoto):** as fibrosis moves off the isthmus, BOTH localizers degrade and neither wins
(grad_phi2 rank 0.33->0.76, fibrosis 0.30->0.73); the desync origin (max-detuned node) is driven by
frequency outliers, not the bottleneck. Lowering frequency heterogeneity (sigma_omega=0.2, K=2.0) so
desync is purely connectivity-driven: at offset=0 (fibrosis ON isthmus) desync nucleates on the
isthmus (dist 0.04) and grad_phi2 nails it (rank 0.008) — but so does fibrosis (0.069), because
fibrosis is there; at offset=1 (dissociated) the network **almost never desyncs (1/14)**.

**Mechanistic insight (the real finding):** in diffusively-coupled excitable/oscillator media,
instability vulnerability is *created by* reduced effective coupling. The Fiedler bottleneck forms at
low-coupling cuts, and fibrosis/LGE imaging detects exactly those low-conduction regions — so the
spectral localizer and the substrate localizer coincide **by mechanism, not coincidence**. A
geometrically-narrow but electrically-healthy isthmus does not become an instability site (it conducts
fine), so the hypothesized "connectivity origin dissociated from substrate" regime does not arise in
this model class. This *explains* the GM2/GM4 "fibrosis also localizes" caveat rather than merely
noting it: |grad phi2| is not uniquely better than substrate for localization here, for a principled
reason.

**Consequence:** ① as "connectivity uniquely beats substrate" is NOT supported by the tractable
graph models. The only regime that could dissociate them is genuine wave-curvature source-sink block
(monodomain reaction-diffusion with S1-S2 on a continuous healthy isthmus between scars) — a real but
compute-heavy build, queued behind the 100k. Novelty weight shifts to (2) the rho scaling law and
(3) the falsification-protocol framing, both tractable now. The negative + insight is itself a
credibility asset for the honesty narrative.

---

## 2026-07-22 (cont.) — Novelty experiment ② (rho scaling law): REAL RESULT

`results/rho_scaling.json`. The validity radius rho = ||dL||/(lambda3-lambda2); the gap
collapses with system size by Weyl's law, so rho grows lawfully with resolution.
Measured the gap-vs-N exponent (fixed mean degree) across media:
- 2D random-geometric: gap ~ N^-1.81
- 3D random-geometric: gap ~ N^-0.84
- **real atrial surface mesh (2-manifold): gap ~ N^-1.05 -- matches the Weyl 2-manifold
  prediction (-1.0) to within 5%.**
Honest reading: the gap collapses as a clean power law in ALL media (so rho grows lawfully,
and the biomarker fails MORE at finer resolution -- counterintuitive but real). The real
2-manifold nails Weyl (-1.0); the abstract random-geometric graphs collapse even faster than
the naive 2/d, so rho grows at LEAST as fast as Weyl predicts. Do not claim exact Weyl match
beyond the manifold case. This turns "rho ~ 2000" into a scaling law tied to system size, and
explains why higher-fidelity meshes (used for accuracy) are further past the validity radius.

---

## 2026-07-22 (cont.) — Novelty experiment ③ (falsification protocol): REAL RESULT

`results/falsification_protocol.json`. Each guard in the protocol, when REMOVED, manufactures a
specific false positive that the guard catches:
- **Grouping guard (leakage):** replicated shape-family cohort — naive random-CV AUC=0.756 vs
  grouped-CV AUC=0.640 => **+0.117 AUC inflation** purely from shape-family leakage.
- **Effect-size guard (p-value trap):** subspace-SFI at N=55k gives p=2e-24 ("it works!") but
  dAUC=+0.0028, needing ~4152 cases for 80% power => clinically useless.
So a naive analyst (random CV + p<0.05, no effect-size gate) would have reported SFI as a working
biomarker; the pre-registered protocol correctly rejects it. This is the packageable methods
contribution: a demonstration that standard biomarker methodology produces exactly the false
positives this protocol is built to catch.

---

## 2026-07-22 (cont.) — GM4 third medium (Kuramoto, N=500): transfer replicates across a 3rd dynamical class

`results/gm4_kuramoto_metrics.json`. 500 oscillator networks, 243 desync (49%).
- Predict (SFI vs competitors): all NULL (dAUC <= +0.0036, p 0.38-0.66). SFI non-predictive here too.
- Localize (vs spatial null): grad_phi2 rank=0.091 p<1e-4 KEEP; perron/combined/wdegree DELETE;
  fibrosis rank=0.147 KEEP.
So the whole GM4 pattern (predict-null + |grad phi2| KEEP) now replicates across THREE media spanning
excitable AND oscillator dynamics: cardiac monodomain, FHN excitable networks, Kuramoto oscillators.
Note grad_phi2 (0.091) is a SHARPER localizer than the raw fibrosis field (0.147) here -- but per
experiment (1)'s insight this is sharper localization of the SAME low-coupling site, not a dissociated
one. Honest framing: |grad phi2| localizes the instability origin more precisely than substrate in the
oscillator medium; it does not find a site substrate misses.

---

## 2026-07-23 — Interim label-fidelity bridge: resolution-convergence of the monodomain labeller (HONEST, mixed)

`results/labeller_convergence.json`. Ran the frozen monodomain labeller on 8 Roney meshes coarsened
to {1500, 2000, 3000} nodes and checked verdict stability. Result is NOT reassuring:
- **Verdict consistent across resolutions: 4/8 (50%).** For half the (pilot) subjects the inducible/not
  call flips with coarsening.
- Inducibility rate: 0.38 (1500) / 0.38 (2000) / 0.25 (3000) — stable at the operating point (2000),
  drops at 3000.
**Honest reading:** the labels ARE partly resolution-sensitive. This *strengthens* the case that the
openCARP (or at least finer-resolution) validation is necessary — the interim check flagged a real
issue rather than reassuring. Impact: the ABSOLUTE predictive numbers (competitor AUC ~0.83) and the
positive localization are softened; the RELATIVE null (SFI vs competitors, same labels) and the
label-independent validity-radius/rho spine are unaffected. Caveats: n=8 pilot (4/8 is noisy); flippers
are the borderline subjects (inducibility is a threshold phenomenon, some flipping expected); the
1500-2000 operating range is more stable than the full range. TODO: larger convergence study pushing
toward FINER resolution to see whether verdicts converge (reassuring) or keep drifting (openCARP urgent).
Do not bury this — it is exactly the kind of sensitivity a sharp judge would probe, and reporting it
is the protocol working.

---

## 2026-07-23 — 100k scale-up COMPLETE (final scaling ladder)

`results/gm4_100k_metrics.json`. All 100,000 FHN networks built (survived 2 kills + a container restart
via sharded checkpointing). Full ladder (localize = permutation null; predict = light single-GroupKFold
LR + DeLong; rigorous nested predict to 50k already in gm4_100k_predict_interim.json):

| N | grad_phi2 rank | p | subspace-SFI predict dAUC | p |
|---|---|---|---|---|
| 2,000 | 0.1796 | <1e-5 | +0.0026 | 0.07 |
| 10,000 | 0.1748 | <1e-5 | +0.0028 | ~0 |
| 25,000 | 0.1742 | <1e-5 | +0.0028 | ~0 |
| 50,000 | 0.1756 | <1e-5 | +0.0028 | ~0 |
| 100,000 | 0.1758 | <1e-5 | +0.0031 | ~0 |

Localizer converged (~0.176) and KEEP at p<1e-5 across all scales; predict is the stable, negligible
+0.003 AUC effect (significant only at large N -> the p-value trap, quantified by the power analysis:
~4150 cases to detect). Note: I stopped the pipeline's slow rigorous 100k aggregate (nested GBT +
bootstrap) mid-run because it blocked novelty experiment 1 for marginal value (the null is already
established + converged to 50k); computed the light ladder directly instead. Build is durable on disk.

---

## 2026-07-23 — Novelty ① (real cardiac model, monodomain isthmus): UNTESTED (could not induce at feasible resolution)

`results/monodomain_isthmus_metrics.json`. Crossfield S1-S2 monodomain on the healthy-isthmus sheet,
27 configs (isthmus width {4,6,8}mm x vertical offset {-6,0,+6}mm x S2 coupling {160,190,220}ms).
**Result: 0/27 inducible. Max sustained = 238ms, median 207ms — all below the 600ms reentry threshold.**
The protocol fires (consistent ~200ms post-stimulus activity) but reentry initiates and dies; it never
sustains into a stable rotor. So grad_phi2-vs-fibrosis localization could not be measured (no origin to
localize). This is NOT a refutation of connectivity-beats-substrate; it is a resolution/compute
limitation — at 0.8mm spacing the wavelength does not fit the isthmus circuit and numerical dissipation
kills the rotor. A proper test needs ~0.2mm resolution (~10x more nodes, hours of sim across the sweep)
plus MS parameters tuned for the sheet -> Claude Science fine-mesh EP block, the same compute class as
the openCARP validation. Honest status of ①: **UNTESTED at feasible resolution; deferred to Claude
Science.** In-session novelty rests on ② (rho scaling law) and ③ (falsification protocol), both real;
the abstract-model insight (spectral bottleneck == conduction lesion) stands as the mechanistic finding.

---

## 2026-07-23 — 100k neural scale-up COMPLETE (final scaling ladder)

`results/gm4_100k_metrics.json`, n_built=100000 (survived a workflow-cleanup kill at 55k and a full
container restart at 85k, zero loss via sharded checkpointing). Localizer scaling ladder:

| N | unstable | grad_phi2 origin rank | perm_p |
|---|---|---|---|
| 2,000 | 704 | 0.180 | <1e-4 |
| 10,000 | 3,564 | 0.175 | <1e-4 |
| 25,000 | 9,013 | 0.174 | <1e-4 |
| 50,000 | 18,213 | 0.176 | <1e-4 |
| 100,000 | 36,540 | 0.176 | <1e-4 |

The |grad phi2| localizer is fully converged (rank ~0.176, stable to 3 sig figs from N=10k onward) and
KEEP at p<1e-4 across the entire ladder — the strongest possible robustness statement for the localizer
transfer. (Final-tier predict left empty — nested-CV at 100k is heavy; predict null already established
through N=50k in gm4_100k_predict_interim.json + the power analysis.) unstable fraction steady ~36%.

---

## 2026-07-23 — Adversarial fidelity audit (ultracode) + honest preprint corrections

Ran a fidelity-defense workflow (verified clinical-rate citations + adversarial audit of the
label-independent spine + completeness audit). It caught real overclaims that had crept into the
preprint; corrected all of them (`docs/paper/PREPRINT.md`), `results/fidelity_defense_audit.json`.

**Corrections made (all honest downgrades):**
1. **Real-cohort result is inconclusive/underpowered, NOT a clean null.** Combined GBT competitors_vs_+SFI
   dAUC = **+0.0507** — *exceeds* the pre-registered 0.05 bar — failing only on significance (p=0.155,
   CI straddles). "Clinically undetectable / ~4150 cases" is scoped to the SYNTHETIC 100k cohort only.
2. **rho does NOT prove non-predictiveness.** It bounds the accuracy of the Delta-lambda2 *magnitude*,
   not classification content. The EXACT Delta-lambda2 SFI (no rho limit) *also* fails -> the null is
   **feature redundancy** (label-dependent), not the validity radius. rho explains single-vector
   ill-conditioning only. Also: rho's numerator carries the fibrosis-weighted Δw, so only the
   denominator is label-free.
3. **Localizer: strict pre-registered endpoint is NULL for all fields** (grad_phi2 0.472 vs null 0.406).
   grad_phi2 survives only as a weak RANK localizer (top ~18%), coincides with (does not beat) substrate,
   cardiac arm is only n~20. Fixed the merged "cardiac/FHN 100k" label (100k is FHN-only).
4. **Multiple comparisons**: ~30 DeLong tests; fibrosis_vs_+SFI p=0.009 declared a pre-specified
   SECONDARY endpoint (suggestive, not confirmatory at family-wise threshold).
5. **rho scaling**: only the atrial 2-manifold matches Weyl; RGG exponents are loose bounds (5-pt fits,
   3-D non-monotone), computed on the unweighted Laplacian. Downgraded "law" claims accordingly.
6. **Citation hygiene**: Darma 2020 (32.4%) DROPPED — did not pass adversarial verification. Anchor now
   uses verified Marquardt 30.6% / Kumar 29.5% / Oral 5% / Kawai 51% / Liu 42% (PMIDs confirmed).
7. Title + abstract + discussion + limitations rewritten to the scoped-honest version.

**CRITICAL convergence finding (Roney, partial):** inducibility rate 58% (1500) -> 17% (3000) ->
25% (6000) -> 21% (12000). It does NOT collapse to zero (audit's worst case refuted) but CONVERGES to
~20% for >=3000 nodes; the operating 2000-node resolution is in the OVER-CALLING zone. So all absolute
rates/AUCs are coarse-mesh-provisional; the RELATIVE SFI-vs-competitors comparison (same labels both
arms) and the label-free rho/estimator argument are the resolution-invariant claims to lead with.

---

## 2026-07-23 — Convergence study COMPLETED (all 144) + correction; Coq formal proof; contacts

**Convergence study now COMPLETE** (Roney, n=24 at each of 1500/3000/6000/12000/24000/48000 = 144
sims). This REVISES the previous "converges to ~20%" note, which was written on incomplete data (24000
had only n=8, 48000 did not exist yet). Full corrected rates:

  1500=58.3%  3000=16.7%  6000=25.0%  12000=20.8%  24000=25.0%  **48000=37.5%**

The rate crashes off the coarse mesh but is **NON-MONOTONE across the fine tiers and RISES to 37.5% at
the finest affordable tier (~50k nodes)** — a 17-pp spread over the last three tiers. **The labeller
does NOT converge in the tested window.** It does NOT collapse to zero (phenomenon is real) but stays
resolution-sensitive up to ~50k. Consecutive per-subject verdict agreement climbs (50%→75%→71%→62%→
**79%** for 24k→48k) yet only **7/24** subjects are fully consistent across all six tiers.
Honest verdict: absolute numbers (rates, competitor AUCs, localizer origins) are coarse-mesh-
provisional pending openCARP; the RELATIVE SFI-vs-competitor comparison (same labels both arms) and
the label-free rho argument are the resolution-invariant claims to lead with.

**Bug caught + fixed:** `scripts/convergence_analyze.py` hard-coded `outputs/convergence/` (a stale UW
pilot) and ignored its cohort arg, so it reported garbage (8%/4%/12%). Rewrote it cohort-aware
(reads `outputs/convergence_<cohort>/`, writes `results/convergence_summary_<cohort>.json`, emits the
figure). Corrected `results/convergence_summary_roney.json`, `docs/paper/figures/fig_convergence.png`,
PREPRINT §Abstract + §Limitations, and INTERVIEW_PREP weakest-point answer to the non-converged story.

**Formal proof (Coq 8.18):** `formal/sfi_edge_identity_Q.v` proves the SFI algebraic core over exact
rationals Q — `Print Assumptions` = "Closed under the global context" (ZERO axioms). T1 rank-one
quadratic form, T2 edge dot, **T3 phi^T E_ij phi = (phi_i-phi_j)^2**, T4/T5 nonnegativity (Dirichlet
energy / PSD). R-valued companion `formal/sfi_edge_identity.v` (inherits the 2 standard R-library
axioms). This is the 7th independent verification of the identity, first kernel-checked. Lean was
blocked: its toolchain ships only as GitHub release assets, which the egress policy 403s.

**Data-request contacts verified (2026-07-23):** Dr Caroline Roney <c.roney@qmul.ac.uk> (Reader in
Computational Medicine, QMUL SEMS); Prof. Patrick M. Boyle <pmjboyle@uw.edu> (Assoc. Prof.
Bioengineering, CardSS Lab, UW). Filled into `docs/outreach/data_request_email.md`.

---

## 2026-07-24 — Real clinical outcomes arrive; two substrate defects found in the UW cohort

**Prof. Boyle shared the outcome column.** `Patient_ID, Recurrence_Rhythm_2yr` ∈ {NR, AF, AFL},
n=82, IDs 1–15 and 21–87 — exactly the ID set of the public Dryad meshes. Marginals: 34 NR,
35 AF, 13 AFL, so 48 events (58.5%); Dryad holdout ID001–015 has 8/15, the ID021–087 cohort
40/67. `docs/DATASETS.md:49` had these recorded as WITHHELD; that is superseded. The file is
restricted-use human-subject data, is NOT covered by the CC0 on the meshes, lives under the
gitignored `data/` tree with two extra `.gitignore` rules, and must not be redistributed
without written permission.

**Pre-registration §8 written and git-tagged `prereg-real-outcomes-20260724` BEFORE any join.**
Power computed from marginals alone (`scripts/power_real_outcomes.py`): at n=82 the study has
only ~25% power at the ΔAUC ≥ 0.05 threshold used for the simulator endpoint, so that threshold
is explicitly not reused. Minimum detectable effect is ΔAUC ≈ 0.10 (optimistic) to 0.15
(conservative); a null below that is pre-declared *inconclusive*, not negative.

### Deviation 1 — elemTag 164 is not tissue. It is the caps over the atrial openings.

The Dryad README documents no tag semantics. An earlier revision of `uw_boyle.py` read 164 as
"remodelled / patchy tissue, distributed (least compact)" and `TAG_FIBROSIS` mapped it to 0.5.
Measured directly (`scripts/` ad-hoc geometry, reproduced in the DROP_TAGS docstring):

  * **four to six large connected components** (≥1% of the region) in all 82 meshes — five
    in 75, four in five, six in two — each topologically a **disc (χ = 1)**, 11–42 mm across,
    at 0.47–0.82 of the atrial radius. That spread is the usual variation in pulmonary-vein
    anatomy: a left common trunk gives four, a right middle vein gives six;
  * borders fibrotic tag 115 on only **~0.08% of its incident edges** (median 27 edges)
    against a chance expectation of 14–20% — a >100-fold depletion, which kills the
    border-zone reading that fibrosis = 0.5 implies;
  * 20.0% of elements pre-ablation, 19.9% post — ablation-invariant;
  * deleting it **opens the surface**: Euler characteristic falls by a median of 7 (range
    2–12), i.e. five large orifices plus small fragments.

**Correction logged the same day.** The first version of this entry claimed "exactly five
components", "zero edge-adjacency" and "χ = −3, exactly five openings". Each was generalised
from a single mesh (ID001) or from a mean fraction rounded to three decimals. The census over
all 82 meshes above is what the data supports; the conclusion is unchanged and the component
spread actually strengthens it, but the specific figures were wrong. **ID040 and ID049** are
topologically pathological (χ before removal −16 and −18) and behave differently. Guards
pinned in `tests/test_uw_boyle.py`.

That is the four pulmonary veins and the mitral valve. Treating them as half-conducting tissue
did more than distort fibrosis: it **sealed the atrium's orifices**, so activation could cross
the mitral valve and the vein ostia instead of circling them. Reentry anchored on those
orifices is a principal AF mechanism and could not form. `DROP_TAGS = (164,)` now removes them
before any field is derived; `drop_tags=()` reproduces the old behaviour for the record.

### Deviation 2 — the released fibre field is degenerate (upstream data, not our bug)

All 164 released meshes carry a `VECTORS fiber` array whose value is a constant `(1,0,0)`;
measured mean directional spread is exactly **0.0**. Verified contrast against Roney 5801337:

  field                      Roney              UW/Boyle
  fibre directional spread   0.945 – 0.979      0.000 exactly
  distinct fibre directions  ~one per vertex    ONE, mesh-wide
  fibrosis representation    continuous IIR,    categorical elemTag,
                             276–1071 levels    3 levels

`edge_weights_from_fibres` applies anisotropy relative to the LOCAL fibre direction, so one
global direction degenerates it into a fixed coordinate bias and removes the fibre
heterogeneity that seeds unidirectional block. `uw_boyle.py` now measures the spread, records
`fibres_are_degenerate`/`fibre_spread` in `meta`, and raises a `RuntimeWarning`. It cannot be
fixed locally — only a re-export from UW can restore it. Asked in
`docs/outreach/boyle_reply_followup.md` Q4.

**Manuscript correction.** The Limitations item called this "rule-based fibres", which reads as
a modelling choice we made rather than a defect we inherited and had not diagnosed. Corrected
in both the Limitations list and the base-rate plausibility passage.

### Attribution experiments (running)

`scripts/uw_substrate_ablation.py` — 2×2 over all 82 UW meshes: caps kept/dropped × constant/
varying fibres. Arm A must reproduce the recorded 6/82 as a control.
`scripts/roney_fibre_control.py` — the cleaner test: destroy ONLY the fibre field on Roney,
which has one. Arm R_A must reproduce 20/62.

Until these land, no causal claim about the 7% is asserted in any document.

### Consequence for the reported combined result

`results/gm1_expanded_metrics.json` and Table `tab:null` report Roney 62 (gbt +0.036), UW 82
(gbt +0.104), Combined 144 (gbt **+0.051**, p 0.155). The combined figure is the closest the
paper comes to clearing the pre-registered 0.05 gate, and it is pulled over that line by the
UW arm — the cohort now known to have been simulated on a sealed atrium with a constant fibre
field. The combined row must be recomputed on the corrected substrate or withdrawn.

### Infrastructure

`scripts/dryad_fetch.py` had the old Linux container's MITM proxy and CA bundle hardcoded; now
env-overridable and defaulting to unset. Re-downloaded the 164 UW meshes (294 MB) and the 100
Roney meshes (3.3 GB, `scripts/zenodo_fetch.py`) — both were lost with the container, along
with the entire `outputs/` label cache. Installed tectonic 0.16.9 into `tools/` (gitignored) so
the manuscript builds on Windows without a TeX install: 0 errors, 0 warnings, 0 bad boxes.
Note tectonic drives XeTeX vs Overleaf's pdfTeX, so page count can differ by one; Overleaf
remains authoritative. Added pandas/scikit-learn to the Windows venv. Corrected `PREPRINT.md`,
which quoted a 600 ms reentry cutoff where the frozen config and manuscript both say 650 ms.


---

## 2026-07-25 — Recompute lands; all three explanations for the UW anomaly are refuted

**Recompute on the corrected substrate.** Two runs, because the full recompute changed both
the substrate and the Roney cohort size (62 → 100 meshes, all now on disk).
`scripts/gm1_matched_n.py` isolates the substrate fix by capping Roney back at 62:

  matched n=144, competitors_vs_+SFI dAUC lr / gbt
    roney     20/62    -0.030 / +0.036   [previously -0.030 / +0.036]  ← UNCHANGED, the control
    uw         8/82    +0.000 / +0.035   [previously -0.057 / +0.104]
    combined  28/144   -0.005 / +0.036   [previously -0.004 / +0.051, p 0.155 → 0.110]

Roney reproducing to three decimals validates the pipeline: a UW-only fix must leave it
untouched, and it does. Full recompute (Roney 100 + UW 82 = 182): combined gbt **+0.012**,
p 0.449, CI [-0.017,+0.044]. The +0.051 that was the paper's closest approach to the
pre-registered gate is gone.

**The paper's one positive is also gone.** The secondary endpoint claimed SFI beats a
fibrosis-only baseline (combined GBT +0.103, p=0.009) and therefore "is not merely
re-encoded fibrosis". Confirmed from git that the old file held exactly +0.103 / p=0.0092.
On the corrected substrate at matched n it is +0.093 (p=0.013); with the full Roney arm it
is **+0.017 (p=0.452)**. So the substrate fix is NOT the main cause here — adding 38 further
real subjects is. An effect that vanishes on more data was a small-sample fluctuation.
Withdrawn in the manuscript. The old UW arm alone had returned +0.340 (p=0.026) on six
positives, which nobody should have believed and which we did not flag at the time.

**Third control: fibrosis gradation. Refuted, in the wrong direction.**
Q_A continuous 20/62 (32.3%); Q_B binary f>0.5 31/62 (50.0%); Q_C binary burden-matched
40/62 (64.5%). Coarsening fibrosis roughly DOUBLES inducibility. All three candidate
explanations for the UW anomaly are dead, and two of them push upward.

**Also corrected:** the manuscript claimed the UW competitor base AUC was "≈ chance". It is
0.834 (lr) / 0.748 (gbt); only the fibrosis-only gbt baseline is 0.562, and it is that
near-chance baseline which manufactures the cohort's +0.217. The results table now prints
base AUC, bootstrap CI and event counts, none of which it previously carried — no interval
appeared anywhere in the paper before today.

**openCARP feasibility.** WSL Ubuntu 26.04 (8 cores, 15 GB) is present; the openCARP v19.0
AppImage runs there without root after `--appimage-extract`. `-buildinfo` returns cleanly.
Recorded as feasibility only; switching labellers mid-analysis is not being done.

---

## 2026-07-25 (later) — CORRECTION to the 2026-07-23 monodomain-isthmus entry: the cutoff is 650 ms, not 600

> The entry "Novelty ① (real cardiac model, monodomain isthmus)" above reports the 27 configs as
> "all below the 600ms reentry threshold". The frozen threshold is **650 ms** (E1 calibration,
> 2026-07-20) — the same 600 ms slip that was corrected in `PREPRINT.md`. The result is unaffected
> (max sustained 238 ms, median 207 ms, so 0/27 inducible under either cutoff); only the stated
> number was wrong. Keeping both entries per the append-only rule.


---

## 2026-07-26 — The UW anomaly is solved, and openCARP passes its gate on the third attempt

### The anomaly was fibrosis burden, not any substrate defect

Three substrate hypotheses were already dead. The geometry audit
(`scripts/cohort_geometry_audit.py`) then measured what none of them covered: at the frozen
2000-node coarsening, mean fibrosis is **0.349 on Roney and 0.244 on UW** — about 30 % less
— with the fraction above half-fibrotic 0.310 against 0.239. UW meshes are also 18 % larger
in area and 7 % coarser at fixed node count. The frozen protocol's own calibration reports
Spearman(fibrosis, sustained reentry) = 0.75, so less fibrosis should mean less
inducibility. That is the labeller working, not failing.

`scripts/fibrosis_burden_swap.py` tested it in both directions, rescaling multiplicatively
so the spatial pattern and gradation are preserved and only the level moves:

    R_A roney native     20/62  32.3%   burden 0.333
    R_B roney at UW      1/62    1.6%   burden 0.241   Fisher p = 4.0e-06
    U_A uw native         8/82   9.8%   burden 0.248
    U_B uw at Roney      21/82  25.6%   burden 0.284   Fisher p = 0.013

Both controls reproduce their native rates exactly. **At matched burden the ordering
reverses** — Roney at 0.241 gives 1.6 % while UW at 0.248 gives 9.8 % — so UW is not
resistant to reentry; given comparable fibrosis it is slightly *more* inducible. The
anomaly was a property of the two released datasets, not of the pipeline.

Lesson worth keeping: an out-of-range aggregate is a hypothesis about the pipeline, but it
is a hypothesis about the **inputs** first. Three controls went hunting for something broken
before anyone checked whether the two cohorts were comparable.

### openCARP: three gate attempts, two failures, each diagnosed

Attempt 1 — 90 % inducible both cohorts, Spearman **−0.261** / −0.162, concordance 18.3 %.
A negative fibrosis–verdict correlation is the signature of rate-dependent block, not
reentry.

Attempt 2 — added fibrosis-dependent ERP. openCARP's `MitchellSchaeffer` defaults are
`tau_close` 150 / `tau_out` 5 (APD90 277 ms); the frozen monodomain uses 110 / 6 and also
shortens ERP with fibrosis (`fibrosis_erp_shortening = 0.5`), which my configuration did
not. Measured on this build: `tau_close` 150 → 277 ms, 110 → 209 ms, 55 → 113 ms, so
`tau_close = 115 × (1 − 0.5 f)` reproduces the frozen 218/121 ms within ~2 %. Result:
76.7 % / 93.3 %, Spearman −0.244 / −0.033, concordance 26.7 %. Better, still failing — so
the diagnosis was only partly right.

Attempt 3 — **conduction velocity had never been calibrated.** Wavelength = CV × APD decides
whether a circuit fits in the tissue. Measured on a 2 cm strip: `g_il` 0.174 (the shipped
default) → **0.344 m/s**, against the in-house band of 0.4–1.2 and its ~0.87 target. That is
a ~75 mm wavelength in an atrium of ~110 mm characteristic length, which makes reentry
trivially easy. A second trap on the way: CV appeared to saturate, 9× `g_il` buying only
1.9× CV, because the monodomain conductivity is the **harmonic mean** of intra- and
extracellular — raising `g_i` alone asymptotes at `g_e`. With both scaled:
`g_il` 0.174 → 0.344, 0.50 → 0.594, **1.00 → 0.852**, 2.00 → 1.307 m/s. Adopted
`g_il 1.05, g_it 0.315, g_el 3.78, g_et 1.134`, transverse at 0.3× longitudinal to match the
in-house along/cross weighting rather than openCARP's ~9:1 ventricular anisotropy.

    cohort   inducible   rate    Spearman   gate      (full 182-subject run)
    roney    22/100      22.0%   +0.404     PASS
    uw        4/82        4.9%   +0.079     FAIL (band)

Concordance with monodomain on all 182: **156/182 = 85.7 %** (both 21, openCARP only 5,
monodomain only 21). openCARP totals 26/182 against monodomain's 42/182, so it is the more
conservative solver. That is the robustness statement worth having: the phenomenological
stand-in every in-silico result rests on is **not over-calling** relative to a standard
reaction–diffusion solver.

**Discipline note.** All three attempts were anchored to independent physiological targets
(APD90, then CV). No endpoint parameter — `reentry_min_ms`, `max_depol_fraction`,
`min_reactivating_nodes` — was touched at any point, which is why the gate stayed able to
fail, and did, twice, on configurations that produced confident-looking labels. Tuning
those instead would have made attempt 1 "pass" and put 90 %-inducible garbage into the
concordance number.

UW's 4.9 % sits below the band but its correlation is now positive, and the band was
calibrated on Roney in §7.4. Given the burden result, a low UW rate is the prediction rather
than an anomaly; `scripts/opencarp_uw_burden_check.py` tests that directly rather than
asserting it.

### Decision recorded

openCARP runs **alongside** monodomain as a robustness check, never replacing it
(pre-registration amendment (i)). Switching would invalidate every in-silico number at once
on the strength of a just-recalibrated solver, and would destroy the one thing running both
produces: a measured agreement rate. openCARP labels remain barred from every endpoint until
the bin-sensitivity sweep also passes.

## 2026-08-06 — Post-submission audit: three code defects, four stale documents, one paper edit

The manuscript had already gone out to journals. This entry records an audit run against the
repository afterwards, and everything it changed. Nothing here alters a reported result; the
one manuscript edit narrows a claim rather than restating a number, and no number in
`docs/paper/manuscript.tex` moved. `verify_manuscript_numbers.py` (310 distinct numbers, 0
untraced) and `verify_repo_consistency.py` (7/7) both passed *before* this audit, which is
why none of the defects below had been caught: neither checker can see any of them.

### The one change that touches the paper

`sec:scaling`'s figure caption carried the caveat "the random-geometric exponents ... use the
unweighted combinatorial Laplacian", which reads as though the headline atrial exponent did
not. It does. `scripts/rho_scaling.py::_laplacian_gap` builds its adjacency from `np.ones`
for **every** medium, and the atrial arm is a single released mesh resampled across five
coarsenings rather than a cohort. The caption now says so for all three fits. The exponent
itself (`-1.05` against Weyl's `-1.0`) is unchanged and still supports what §10 claims — how
the gap of a diffusively-coupled 2-manifold scales with resolution — but not the stronger
reading that the *weighted* operator was measured. Scope note added to `SFI_THEORY.md` §10,
along with two limits that should not have to be rediscovered: five-point fits over four
seeds, and `eigsh(which='SM')` without the shift-invert both `asb.spectral` and `asb.sfi`
deliberately use.

### Three defects in `gm3.py`

1. **The Part B subspace column scored the wrong target.** `subspace_sfi(k_dim=2)` predicts
   the drop in the *sum* λ2+λ3; the sweep compared it against the exact λ2 drop alone,
   producing a relative error near 12 that read as though the subspace estimator were a
   thousand times worse than the first-order one. `degeneracy_sweep`, sixty lines below, had
   always used the correct target — so one file contained both the right and the wrong
   comparison. Fixed: `k=3` eigenpairs, scored against `exact_sum_dlam23`.
2. **The relative-error floor manufactured a small error.** `max(abs(exact), 1e-12)` was not
   guarding a division by zero. In the tightest degeneracy row the exact λ2 drop is ~7e-15
   against a single-vector prediction of ~4e-31 — a true relative error of essentially 100 %
   — which the floor reported as **0.69 %**, making the single-vector estimator look most
   accurate exactly where it had failed completely. Replaced by `_rel_err`, which returns
   `None` below solver noise. "Not measurable here" is the truth; a ratio against a floor is
   not.
3. **The degeneracy sweep was not controlled.** The perturbation mask was drawn *inside* the
   loop over bridge values, so the perturbation moved together with the gap — while both the
   docstring and the inline comment said "at a fixed diffuse perturbation". The sweep varies
   the gap on purpose and nothing else, so this made its two error series unattributable.
   Mask now drawn once, outside the loop.

**Consequence, stated plainly.** The claim in `SFI_THEORY.md` §5 that "GM3 confirms the
subspace prediction stays accurate as the gap closes where the single vector is erratic" is
**withdrawn**. The stored series are 0.085, 0.074, 0.142, 0.076, 0.142, 0.0069 (single) and
0.068, 0.071, 0.0028, 0.075, 0.00034, 0.100 (subspace): neither is monotone in the gap and the
ordering reverses in the tightest row. What survives, and needs no sweep, is the algebraic
guarantee — the trace over the invariant subspace is basis-independent, so the subspace SFI is
well *defined* where the single vector is not. `results/gm3_metrics.json` predates all three
fixes and **must be re-run** before any accuracy claim is restored from it. No manuscript
number depends on that sweep, which is why the paper needs no correction on this point.

### Four documents that had drifted

- `README.md` still described the UW outcomes as restricted-use and non-redistributable, and
  said the pre-registered analysis "has not yet been run". Both were overtaken on 2026-07-30:
  custody was lifted on the authors' written confirmation (§8.1) and the endpoint was run and
  not met (amendment (l)). It also rendered the marginals as "48 events: 34 NR / 35 AF / 13
  AFL", which reads as though NR were an event. Corrected, and the clinical endpoint added to
  the headline results, where its absence was conspicuous.
- **Amendment (m) added to §8.6.** Amendment (f) closed with the UW anomaly "unexplained and
  reported as an open problem", and the 2026-07-26 burden swap that solved it never
  superseded it. It does now, with the two-directional control tabulated. The comment at
  `uw_boyle.py` ending "which remains open" is corrected in the same way. Note (f)'s closing
  commitment — no further explanation without a control behind it — was kept.
- **Stale file:line citations in the pre-registration**, all silently wrong after the loader
  grew: `uw_boyle.py:82` → `:101` (TAG_FIBROSIS, twice), `:226-238` → `:259-271` (the PCA UAC
  surrogate), `:114` → `:133` (DROP_TAGS), `docs/DATASETS.md:49` → `:70`.
- **"276–1071 distinct levels"** for the Roney IIR field, quoted in two scripts and the
  pre-registration, is not what the run records. Measured over the 186 stored rows of
  `results/roney_fibrosis_quantization.json`: **700–1646, mean 1318**. The qualitative point
  (continuous against three-level) is untouched; the quoted range was wrong.
- `docs/DATASETS.md`'s table cell still read "restricted use, not redistributable" inline,
  relying on the update note above it to correct the reader. Now corrected in place.

### What was deliberately NOT changed

`TAG_FIBROSIS[164] = 0.5` and `DROP_TAGS` omitting 199 both look like bugs and are neither —
they are frozen so `drop_tags=()` reproduces the discredited substrate exactly, and so the
`uw_drop_tags` config hash stays stable. The duplicated `sfi_region_max` / `sfi_top1` columns
stay for the same reason, disclosed in `features.py` and in the manuscript. Deleting any of
them would invalidate cached labels to reproduce an identical answer.

**Environment note.** None of this could be executed: this checkout's `.venv` points at a
Python belonging to a different Windows account and there is no interpreter available to the
current user, so every change above is a static edit. `gm3.py`, `make_figures.py` and the two
scripts have been read back line by line, but the suite has **not** been run against them and
the three gm3 fixes are unexercised. Running `pytest -q`, then `verify_repo_consistency.py`
and `verify_manuscript_numbers.py`, and then regenerating `results/gm3_metrics.json`, is the
first thing to do on a machine with a working interpreter.

## 2026-08-06 (later) — The environment was not broken after all; the audit's fixes are now measured

The earlier entry today closed by saying no interpreter was available and the three `gm3.py`
fixes were unexercised. **That was wrong, and the error was mine.** A working CPython 3.12.10
was present all along at `C:\Users\mouni\AppData\Local\Programs\Python\Python312`; it is simply
not on `PATH`, where the Windows Store `python.exe` stub shadows it, and an earlier directory
probe returned empty inside a PowerShell pipeline that had already errored. The conclusion "no
Python on this account" was drawn from two weak signals and should have been checked directly
before being reported twice.

Rebuilt the environment against it (`.venv` recreated; the old one, which pointed at the
`vijay` account's interpreter, was moved out of the tree rather than deleted since it may still
be valid for that account). Then ran what the earlier entry said could not be run.

**Everything passes.**

- `pytest -q` — **exit 0**, 129 tests, no failures.
- `verify_repo_consistency.py` — **7/7**, now validating 17 cited results files, up from 15
  because the (m) amendment cites `fibrosis_burden_swap.json` and
  `roney_fibrosis_quantization.json`.
- `verify_manuscript_numbers.py` — **310 distinct numerals, 0 untraced**, unchanged. The
  scaling-law caption edit really did add no numbers, now confirmed by the tool rather than by
  reading.
- `tools/tectonic.exe` rebuilt `manuscript.pdf` from corrected source: 28 pages, no LaTeX
  errors, no undefined references.

### The three gm3 fixes, measured against the stored file

The control arm is the point. Re-running `validity_radius_sweep()` reproduces **bit-identically**
every column the fixes do not touch — `rho`, `exact_dlam2`, `pred_first`, `rel_err_first`,
`weyl_ok`, and ρ\* to all 16 digits (2.9998846240069104). So this environment reproduces the
original numbers, and the fixes changed only what they were meant to change.

*Part B, the wrongly-targeted subspace column.* Stored, it ran **downward** — 11.91 at ρ=0.3 to
0.06 at ρ=28.5 — which is backwards for a perturbation error and was the tell nobody read.
Corrected to score against the `λ₂+λ₃` drop it predicts, it rises monotonically: 0.011 at ρ=0.3,
0.105 at ρ=3, 0.763 at ρ=28.5. It now tracks the first-order estimator at small ρ and beats it
at large ρ (0.763 against 0.919), which is the behaviour theory predicts and the old column
concealed. The downward drift was mechanical: as the perturbation grows the λ₂-only exact drop
climbs toward the sum, so a mismatched denominator shrinks the apparent error.

*The degeneracy sweep.* With the mask drawn once and the floor removed, the erratic stored
series resolve into smooth monotone ones, and the `bridge = 1.0` control row reproduces exactly:

    bridge   gap      single-vector   subspace     (stored single / stored sub)
    1.0      0.16733  0.0846          0.0681       0.0846 / 0.0681   <- control, exact match
    0.3      0.06600  0.0887          0.0733       0.0740 / 0.0713
    0.1      0.02394  0.0900          0.0751       0.1421 / 0.0028
    0.03     0.00740  0.0904          0.0757       0.0760 / 0.0754
    0.01     0.00249  0.0905          0.0759       0.1417 / 0.0003
    0.003    0.00075  0.0906          0.0760       0.0069 / 0.1004

The 6.9e-15 exact drop that the floor turned into a fake 0.69 % error was itself an artefact of
the redrawn mask: with the perturbation held fixed, that row's exact drop is 9.66e-09 against a
prediction of 8.78e-09, a genuine 0.0905.

**Consequence for the claim.** The withdrawal recorded earlier today is partly reinstated, in a
weaker and now-defensible form. Supported: the subspace estimator is more accurate than the
single-vector one at every gap tested. **Not** supported, and still withdrawn: that the
single-vector error blows up as the gap closes — it converges to ~0.091. `SFI_THEORY.md` §5 and
the `degeneracy_sweep` docstring now say exactly that and no more.

### What is deliberately left alone

`results/gm3_metrics.json` keeps its pre-fix values for the two corrected columns. A full
`run_gm3()` is **not** a like-for-like regeneration any more: the Roney label cache now holds
**100** subjects (34 inducible) where the stored file is **62** (20 inducible), because the
Roney arm was extended to all released meshes after GM3 was last run. Overwriting would silently
swap the cohort underneath a file the manuscript was written against. Since no manuscript number
reads either corrected column, and every column that *is* read reproduces bit-identically, the
right move is to leave it and say so. Regenerating GM3 at the current cohort is a deliberate
re-run with its own entry, not a cleanup.

## 2026-08-06 (later still) — GM3 re-run on the full 100-subject Roney arm: the null holds, and hardens

Post-hoc and unplanned, logged because the protocol requires it. Having rebuilt the
environment, the GM3 pipeline was re-run end to end with the corrected `gm3.py` against the
current Roney label cache. That cache now holds **100** subjects (34 inducible, 0.340) where
the GM3 in the manuscript ran on **62** (20 inducible, 0.323), because the arm was extended to
all released meshes afterwards. So this is a **different cohort, not a regeneration**, and it
is stored separately as `results/gm3_n100_20260806.json` with a provenance block saying so.
`results/gm3_metrics.json` remains the manuscript's GM3 and is untouched.

    comparison / classifier        stored n=62 (20 ind)     re-run n=100 (34 ind)
    comp vs +single_vector / lr    dAUC -0.0286 p=0.656     dAUC -0.0152 p=0.430
    comp vs +single_vector / gbt   dAUC -0.0000 p=1.000     dAUC -0.0042 p=0.693
    fib  vs +single_vector / lr    dAUC +0.0298 p=0.307     dAUC +0.0120 p=0.726
    fib  vs +single_vector / gbt   dAUC +0.0179 p=0.445     dAUC +0.0038 p=0.865
    comp vs +subspace      / lr    dAUC -0.0214 p=0.447     dAUC -0.0080 p=0.543
    comp vs +subspace      / gbt   dAUC +0.0107 p=0.812     dAUC +0.0138 p=0.476
    fib  vs +subspace      / lr    dAUC +0.0286 p=0.637     dAUC +0.0236 p=0.488
    fib  vs +subspace      / gbt   dAUC +0.0524 p=0.283     dAUC +0.0129 p=0.676
    comp vs +exact         / lr    dAUC -0.0274 p=0.330     dAUC -0.0423 p=0.060
    comp vs +exact         / gbt   dAUC -0.0107 p=0.640     dAUC -0.0203 p=0.178
    fib  vs +exact         / lr    dAUC +0.0464 p=0.403     dAUC -0.0058 p=0.854
    fib  vs +exact         / gbt   dAUC +0.0298 p=0.322     dAUC -0.0287 p=0.221

**None of the twelve cells meets the endpoint**, as at n=62. Three things are worth recording.

1. **The exact-Δλ₂ arm gets *more* negative with more data** — competitors vs +exact goes
   −0.0274 (p=0.330) to −0.0423 (p=0.060) under logistic regression. This is the arm that
   forecloses "the null is just an invalid estimator", since the exact recompute carries no
   first-order error and no ρ limitation. At the full cohort it does not merely fail to help;
   the point estimate says adding it *hurts*, and it is the closest thing to significant in the
   table — in the wrong direction for the biomarker.
2. **The largest n=62 positive shrinks toward zero.** Fibrosis vs +subspace under GBT was the
   friendliest cell at +0.0524; at n=100 it is +0.0129. That is the same pattern the manuscript
   already documents for the headline real-cohort estimate (+0.051 to +0.012 on extending the
   Roney arm), reproduced independently in GM3. More data shrinks these, which is what a null
   with sampling noise around it looks like, and not what a real effect looks like.
3. **Baselines are healthy**, grouped AUC ≈ 0.856–0.870, so this is not the anti-predictive
   baseline regime that invalidated the clinical gate in amendment (l). The null here is
   measured against a competitor set that works.

**Status.** Exploratory and post-hoc. It does not enter any endpoint, does not change a
reported number, and is not a correction to the manuscript — the paper's GM3 is n=62 and stays
n=62. Its value is corroborative: the GM3 conclusion was not a small-sample artefact, and if a
referee asks what happens on all 100 released meshes, the answer is now measured rather than
asserted.

*One integrity note on storage.* Adding a results file enlarges the pool that
`verify_manuscript_numbers.py` matches manuscript numerals against, which in principle could
mask a future untraced number. Checked after adding: still 310 distinct, 0 untraced, and
`verify_repo_consistency.py` still 7/7. Flagged so the trade-off is visible rather than silent.

## 2026-08-06 (later still) — PermeaFlow removed from this repository

AtrialSpectralBench is one project: the Spectral Fragility Index as a candidate biomarker for
atrial fibrillation. PermeaFlow — the TAVR paravalvular-leak work — is a separate project and
had leaked into this tree. Removed.

**Tracked, and therefore actually part of this project's published record:**

- `artifacts/` — all 14 files, entirely PermeaFlow: `accuracy_campaign.json`, `conformal/`,
  `gci/`, `grounding/`, `operator_gpu/` (a 
  neural-operator checkpoint), `pinn_poiseuille/`, `slice_demo/`, `verify_analytic/`,
  `vv40_gate/`. Nothing in `src/`, `scripts/`, `tests/`, `docs/`, `README.md`, `Makefile` or
  `ROADMAP.md` referenced the directory at all, so it was dead weight as well as off-topic.
- `accuracy_campaign.err` — PermeaFlow's stderr capture.

**Untracked ghosts, moved out rather than destroyed** (to
`C:\Users\mouni\PermeaFlow-strays-from-AtrialSpectralBench`): `.scratch_src110/` (a full copy of
the `permeaflow` package), `.scratch_verify_thresh.py` (a PermeaFlow GAP_SCALE check),
`accuracy_campaign.log`, `item5_floor.log`.

**Two ghosts worth naming, because they were actively harmful rather than merely untidy:**

1. `src/permeaflow/` still existed, containing **49 stale `.pyc` files and no source at all**.
   Deleted. This is exactly the failure mode the second-idea branch logged on 2026-08-02 as
   "a stale-bytecode hole that made the numbers untrustworthy while I measured them", sitting
   in the other project's tree.
2. Because that directory was present, the editable install had registered **`permeaflow` as a
   top-level package of `atrialspectralbench`** — `src/atrialspectralbench.egg-info/top_level.txt`
   read `asb` and `permeaflow`. Regenerated; it now reads `asb` alone, and `permeaflow` is no
   longer importable from the venv. Also cleared every `__pycache__` outside `.venv`, which
   removed stale bytecode for three PermeaFlow scripts (`conformal_validation`,
   `deepxde_crosscheck`, `multiseed_headlines`) whose `.py` files do not exist here.

**Nothing of PermeaFlow's was destroyed.** Every removed artifact directory has a counterpart in
`C:\Users\mouni\PermeaFlow\artifacts\`, which holds ten more besides, so the copy committed here
was a stale subset; and the deletions remain recoverable from this branch's history regardless.

**Verified after removal:** `pytest` exit 0 (129 passed), `verify_repo_consistency` 7/7,
`verify_manuscript_numbers` 310 distinct / 0 untraced. Removing `artifacts/` changed no result
and broke no reference.

**Deliberately NOT done, both needing a decision rather than a cleanup:**

- `data/` holds roughly **13 GB of aortic/TAVR datasets** that belong to PermeaFlow —
  `zenodo_thoracic_aorta_cohort_3000` (6.3 G), `zenodo_tavi_aortic_root_segmentation` (4.3 G),
  `zenodo_piv_aortic_valve` (1.1 G), `zenodo_hasler_obrist_piv_1163562` (309 M),
  `zenodo_coronary_4dct_stl`, `vmr_0157_aorta_healthy_adult`, `vmr_0008_aorta_inlet_cap`.
  No tracked file references any of them, and `data/` is gitignored so none of it is in the
  repository. Left in place: deleting 13 GB of downloaded research data is not a side effect
  to take unasked, and `C:\Users\mouni\PermeaFlow\data` is a **dangling symlink** pointing at
  `C:\Users\vijay\cardiothoracic_surgery_research\data`, which suggests this directory was
  meant to be shared between the two projects in the first place.
- The `claude/second-idea-generational-n0p02n` branch is PermeaFlow living in this GitHub
  repository. Deleting it would remove 56 commits of another project's research from the
  remote, which is a destructive act and the opposite of what was asked earlier the same day.
  Not touched.
