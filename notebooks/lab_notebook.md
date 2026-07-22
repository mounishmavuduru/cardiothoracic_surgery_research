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
