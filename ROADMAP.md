# ROADMAP — from here to a finished, defensible research paper

This is the complete plan for finishing **AtrialSpectralBench**: what is done, what
remains for **Claude Code** (CPU-feasible, do now), what is deferred to **Claude
Science** (GPU/heavy, wired as one-command entry points), how the **paper** is
assembled, and how to **rehearse the interview**. It supersedes nothing in
`EXPERIMENTAL_PLAN.md` / `docs/PRE_REGISTRATION.md` — it operationalizes them given
the results now in hand.

---

## 0. The honest scientific status (read this first)

The pre-registered predictive hypothesis was **tested rigorously and falsified**, and
that is a *valid, pre-registered, reportable outcome* — not a failure. What survives is
real and defensible. The paper's spine is therefore **not** "a new POAF biomarker"; it
is **"exactly when a closed-form spectral-fragility surrogate is trustworthy for
nonlinear excitable-media instability, and what genuinely transfers."**

| Generational-maker | Pre-registered claim | Result | Status |
|---|---|---|---|
| **GM1 — predict** | SFI beats competitors on inducibility (ΔAUC≥0.05, DeLong p<0.05) | grouped ΔAUC ≈ 0/negative vs full competitor set; exact & subspace SFI do **not** rescue it | **NULL (robust)** |
| **GM2 — localize** | `\|∇φ₂\|∩Perron` finds reentry origin < spatial-null 5th pct | strict endpoint **null**; but **KEEP `\|∇φ₂\|`** (perm p<1e-3), **DELETE Perron & the ∩ hotspot** | **partial positive** |
| **GM3 — validity radius** | quantify where the linear surrogate holds | ρ\*≈3 (first-order 10% error), Weyl bound holds, subspace handles degeneracy; **real atria at ρ≈2422** | **HOLDS (lead result)** |
| **GM4 — transfer** | fragility calculus generalizes to a 2nd excitable medium | full pattern replicates in an FHN neural net: predictive **null** + `\|∇φ₂\|` **localizer KEEP** (p<1e-3, beats Perron/degree) + validity radius (ρ≈1852) | **HOLDS (convergent)** |

**Three genuine contributions** the write-up defends:
1. **The validity radius (GM3)** — a quantified, theory-matched boundary (Weyl /
   Davis–Kahan) for when the seconds-on-a-laptop derivative is a faithful surrogate for
   an overnight PDE, plus the empirical demonstration that real fibrotic atria sit
   ~800× outside it. This is the dossier's stated "honest core."
2. **A measured mechanism correction (GM2)** — the Fiedler gradient `|∇φ₂|` genuinely
   localizes reentry origins, but the hypothesized Perron co-localization is *falsified*
   and actively hurts. We kept/deleted each claim by measured correlation, as promised.
3. **Cross-medium generality (GM4)** — the *pattern* replicates in an independent
   excitable medium (a FitzHugh–Nagumo neural network with an epileptic-focus
   instability): the linear SFI is a null *predictor*, the **Fiedler gradient `|∇φ₂|`
   localizes** the instability origin above the spatial null in *both* media (p<1e-3,
   beating Perron/weighted-degree centrality controls), and the validity radius is the
   same (ρ≈1852 ≫ ρ*≈3). Establishes the calculus as graph-universal, not
   cardiac-specific. *(An initial GM4 run suggested a cardiac↔neural localizer
   dissociation; adversarial audit traced it to an origin-definition artifact and the
   corrected analysis shows convergence — the correction is logged in the notebook and is
   itself a credibility asset.)*

**Honesty guardrails locked in** (never violate in the write-up): labels are simulator
verdicts, never clinical POAF; no FDA/VICTRE/in-silico-clinical-trial language; the
`dλ₂/dw=(φ_i−φ_j)²` identity is credited as classical (Ghosh–Boyd 2006) and the
λ₂↔excitable-stability link as prior art (Pecora–Carroll 1998; Bomela 2020) — our
contribution is the *spatial-localization use* and its *transfer*, not the formula.

---

## 1. What is DONE (Claude Code, this environment)

- **Engine** (`src/asb/`): spectral, SFI (incl. subspace + exact-Δλ₂), baselines,
  features, evaluation (DeLong, grouped/naive GroupKFold, bootstrap, colocalization),
  figures, pipeline, dashboard. **All frozen to `BUILD_SPEC.md`.**
- **Verification / tests:** analytic spectral gates (path/ring/grid `λ_k=2−2cos`),
  derivative-identity finite-difference, subspace fallback, Perron positivity — as
  **property-based tests generating 10,600 random graphs/run** plus ~100 example tests.
  All green.
- **E1 real substrate:** real Roney LA cohort loader (Zenodo 5801337; UAC + IIR-fibrosis
  + fibre arrays), field-preserving coarsening. 62 patients labelled.
- **Ground-truth labeller:** a genuine **monodomain Mitchell–Schaeffer** reaction–diffusion
  solver (cotangent-Laplacian diffusion), calibrated (inducible ≈0.32, Spearman(fibrosis,
  reentry)=0.75), `source='monodomain_ms'`.
- **GM1** (`experiments/gm1.py`), **GM2** (`gm2.py`), **GM3** (`gm3.py`), **GM4**
  (`gm4.py` + `asb/transfer/`) — all config-driven, cached, one-command.
- **Pre-registration** frozen + git-tagged **before** any label was inspected; dated lab
  notebook kept append-only.
- **Δw calibration** frozen to literature perioperative CV slowing (f_CV≈0.20 → f_D=1−0.8²
  ≈ **0.36**), recorded in the pre-registration.

---

## 2. What REMAINS for Claude Code (CPU-feasible — do these next, in order)

These need **no** GPU and can all run in this kind of session. Ordered by leverage.

### 2.1 Finish + verify GM4 (in progress)
- Complete the 80-network FHN run; **adversarially verify** (subagent) for circularity
  (label must not touch λ₂), leakage (each network its own group), and that the
  "generality" claim is stated honestly. Record verdict in the notebook.
- If Part B `|∇φ₂|` KEEPs in the neural medium too → the transfer headline is real.

### 2.2 E6 — Uncertainty quantification (descoped, CPU)
- **Morris elementary-effects screening** over the 4 pre-registered knobs: conduction
  scale `d0` (CV), fibrosis/IIR threshold, fibre anisotropy (`along/cross`), and the Δw
  magnitude. Outcome metrics: inducible fraction, GM1 grouped ΔAUC, GM2 `|∇φ₂|` perm-p,
  GM3 ρ. Report which inputs the **conclusions** (not just λ₂) are robust to. ~dozens of
  monodomain runs → hours on CPU; cache aggressively.
- **Gaussian-process surrogate** over the same 4-D box trained on the Morris runs, to
  interpolate the inducible-fraction and ΔAUC response surfaces. `scikit-learn`
  `GaussianProcessRegressor`. Report main effects + a screening-level sensitivity bar.
- **Δw sensitivity** specifically: re-run GM1/GM3 at Δw ∈ {0.18, 0.36, 0.54} to show the
  null and the validity-radius conclusion are not artifacts of the frozen Δw.
- New module: `src/asb/experiments/e6_uq.py` (Morris + GP), one-command.

### 2.3 E7 — Credibility case (CPU, mostly analysis + writing)
- **Literature range-checks** (`experiments/e7_credibility.py`): confirm our numbers sit
  in published ranges — inducible fraction vs pooled post-CABG POAF ~26% (PMC10823463);
  monodomain planar CV in 0.3–1.2 m/s (a solver-verification gate we can compute now via
  `measure_planar_cv`); AUC magnitudes vs Bieging (0.64→0.78) / Varela (0.68–0.71).
- **Anatomy-realism check:** compare Roney LA shape/UAC statistics to the public 2018 LA
  Segmentation Challenge geometries (download a handful; overlap of size/sphericity
  distributions). Anatomy realism **only**, never clinical validation.
- **Independent-geometry generalization:** re-run the spectral/SFI feature extraction on
  unrelated open geometries (Vascular Model Repository .vtp) to show the pipeline isn't
  overfit to one shape model (report feature ranges, not a new AUC).
- **Convergence (calculation verification):** we already coarsen to ~2000 nodes; add a
  mesh-refinement + Δt-reduction study showing λ₂ and the inducibility label are stable
  (ASME V&V-40) — `experiments/e0_convergence.py` exists; extend + run on 3–5 meshes.

### 2.4 Figures (publication-quality, matplotlib) — `experiments/make_figures.py`
The dossier's five + our result figures:
1. φ₂-colored atrium (real Roney mesh).
2. `|∇φ₂|` hotspot vs monodomain reentry origin (the GM2 KEEP), with the spatial-null.
3. Eigenvalue spectrum + the λ₂–λ₃ gap distribution across the cohort.
4. ROC panel: competitors vs competitors+SFI (grouped), both classifiers (the GM1 null,
   shown honestly with CIs).
5. **The validity-radius curve** (GM3): first-order/subspace/exact error vs ρ, with ρ\*
   and the real-cohort ρ marked — this is the money figure.
6. **GM4 transfer panel:** the FHN network colored by `|∇φ₂|` with the instability origin,
   + a side-by-side "atrium vs neural network" keep/delete comparison.
7. E6 sensitivity bars.

### 2.5 Deliverables / writing scaffolds (Claude Code can draft)
- **SFI theory note** (`docs/SFI_THEORY.md`): the full perturbation derivation
  (`L=Σw(e_i−e_j)(e_i−e_j)ᵀ` ⇒ `dλ₂/dw=(φ_i−φ_j)²`), the second-order resolvent term, Weyl
  vs Davis–Kahan, the subspace projector SFI, and the validity radius — with the
  Ghosh–Boyd / Pecora–Carroll prior-art credits stated up front.
- **`results_report.md`**: auto-assembled from all `*_metrics.json` (a small script that
  renders every GM's table + the honest narrative).
- **Dashboard** polish: add the real-Roney + FHN-network views and the validity-radius
  slider.
- **README** update: the pre-registered hypothesis, the honest outcome, reproduction
  commands.

---

## 3. What is DEFERRED to Claude Science (GPU/heavy — wired, not run here)

Each is already (or should be) a **one-command config-driven entry point** so Claude
Science runs it with a single call. None changes the code's logic — only scale/solver.

### 3.1 Real openCARP labels (replaces `monodomain_ms`)  ★ highest priority
- **Reproduce the Niederer N-version benchmark** (PMID 21969679) in openCARP — the E0
  gate. Interface already stubbed in `labels/opencarp.py::reproduce_niederer_benchmark`.
- **Full monodomain inducibility sweep**: openCARP + Mitchell-Schaeffer, S1–S2/burst from
  multiple sites, pre/post-stress, stochastic realizations, over the **whole** Roney (and
  optionally Rodero-derived) cohort. Emit `source='opencarp'` labels via
  `labels/opencarp.py::run_opencarp`.
- **Re-run GM1/GM2/GM3 on openCARP labels** (same `run_gm*` entry points, just point the
  cohort at openCARP labels): confirm the null / the `|∇φ₂|` localization / the validity
  radius hold with the *real* solver. **This is the single most important deferred step**
  — it upgrades every result from "monodomain-MS verdict" to "openCARP verdict."
- **Label-concordance cross-check**: agreement of `monodomain_ms` vs `opencarp` on a
  subset (Cohen's κ) — quantifies how good the CPU surrogate was.

### 3.2 Full UQ at scale
- Full **Morris + Sobol/GP** over the 4-D parameter box with openCARP labels (thousands of
  PDE evals) — the compute-heavy version of §2.2.

### 3.3 Bayesian hierarchical validation
- Hierarchical model of ΔAUC / localization across shape families with credible
  intervals (partial pooling), to formalize "no incremental value" with a posterior.

### 3.4 Scale the cohort + GM4 generality
- Grow to **hundreds of distinct real anatomies** — verified candidate datasets are catalogued
  in **`docs/DATASETS.md`**; add the **UW/Boyle Dryad set** (`10.5061/dryad.kkwh70sg0`, ~82
  distinct real patients with LGE fibrosis, `.vtk`) first (nearly triples distinct real N), then
  the Nagel bi-atrial SSM. Tag patient = group for grouped CV (variants of one shape ≈ 1 group).
- GM4 at scale: more networks and **more topologies** (small-world Watts–Strogatz,
  scale-free) + a second neural model (Wilson–Cowan / Epileptor node) to show the transfer
  is not FHN-specific.

### 3.5 Publication-scale figures
- Re-render Figures 1–7 at journal resolution with the full cohort.

---

## 4. The paper — structure and how results map to claims

Target: an open preprint (arXiv/bioRxiv) + the ISEF board. Suggested structure:

1. **Title / abstract** — "When is a closed-form spectral-fragility derivative a
   trustworthy surrogate for nonlinear excitable-media instability? A pre-registered
   in-silico study across cardiac and neural media." Lead with the validity radius +
   the honest null.
2. **Intro** — the white space (spectral perturbation as the object of study), the
   prior art credited honestly (Zahid min-cut, Sun GFT, Falkenberg percolation, FibMap;
   Ghosh–Boyd, Pecora–Carroll).
3. **Methods** — substrate (real Roney LA + coarsening), the SFI math (+ theory note),
   the monodomain-MS labeller (+ openCARP deferred), the competitor set, the
   leakage-controlled evaluation, the pre-registration + Δw calibration.
4. **Results** —
   - E0/E1 verification + calibration (inducible fraction, CV range, convergence);
   - **GM1 null** (grouped & naive, both classifiers, DeLong, bootstrap) — reported as a
     pre-registered accepted outcome;
   - **GM2** — the `|∇φ₂|` KEEP, Perron/hotspot DELETE, with the permutation null;
   - **GM3** — the validity radius (the lead figure);
   - **GM4** — the transfer result;
   - E6 UQ robustness.
5. **Discussion** — why fragility collapses to fibrosis/connectivity as a *predictor* but
   `|∇φ₂|` still *localizes*; the validity radius as the reusable contribution; what the
   openCARP re-run (Claude Science) would confirm.
6. **Limitations** — simulator verdict not clinical POAF; monodomain-MS not openCARP
   (deferred); N=62; single fibre model; Δw calibration uncertainty (with the E6
   sensitivity).
7. **Reproducibility** — repo, seeds, configs, commit hashes, dated notebook, one-command
   `run_gm*`.

**The one honest sentence** (replaces the falsified headline): *"A closed-form per-edge
eigenvalue-sensitivity map is a faithful surrogate for the exact connectivity drop only
inside a validity radius ‖ΔL‖ ≲ 3(λ₃−λ₂); real fibrotic atria sit ~800× outside it, so the
linear Spectral Fragility Index does not beat fibrosis-heterogeneity as an inducibility
predictor — yet the Fiedler gradient still localizes reentry origins above a spatial null,
and the same behavior replicates in a neural excitable medium."*

---

## 5. Interview / poster rehearsal (25 of 100 ISEF points — budget the last month)

- **Own every symbol**: derive `dλ₂/dw=(φ_i−φ_j)²` on the whiteboard; explain Weyl
  (eigenvalue) vs Davis–Kahan (subspace) and why the subspace SFI exists.
- **Rehearse the null**: a strong candidate *defends a negative result*. "We pre-registered
  it, we did not tune, and here is the mechanism (ρ≈2422)." Judges reward this.
- **Anticipate the 7 killer objections** (thin novelty, circularity, linear-vs-nonlinear,
  over-claimed ground truth, compute-infeasible UQ, perturbation misuse, near-degeneracy +
  leakage) — each is pre-answered *in the build*; know where.
- **Two-audience poster**: a math panel (perturbation theory + validity radius) and a
  bio/CBIO panel (atrial substrate + the transfer). Figures 5 (validity radius) and 6
  (transfer) are the centerpieces.
- **Dated notebook + preprint** ready to show provenance.

---

## 6. Concrete ordered checklist

- [x] Engine, tests (10.6k property graphs), pre-registration + tag, E1 substrate,
      monodomain-MS labeller, GM1, GM2, GM3.
- [ ] **GM4 finish + verify** (running).
- [ ] E6 UQ (Morris + GP + Δw sensitivity) — `experiments/e6_uq.py`.
- [ ] E7 credibility (literature ranges, 2018 LA-Challenge realism, VMR generalization,
      convergence study) — `experiments/e7_credibility.py`.
- [ ] Figures 1–7 — `experiments/make_figures.py`.
- [ ] `docs/SFI_THEORY.md` + auto `results_report.md` + README + dashboard polish.
- [ ] Draft the preprint + poster.
- [ ] **[Claude Science]** openCARP: Niederer benchmark → full inducibility sweep →
      re-run GM1/GM2/GM3 on `opencarp` labels → concordance vs monodomain-MS.
- [ ] **[Claude Science]** full UQ, Bayesian hierarchical validation, cohort scale-up,
      GM4 multi-topology/multi-model, journal-resolution figures.

*Keep committing to `claude/third-idea-project-plan-rmerjr`; keep the lab notebook
append-only; never tune to an endpoint; report every null.*
