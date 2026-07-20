# Idea 3 — Spectral Graph Theory for Post-Operative AF: Project Plan

**Project name (working):** AtrialSpectralBench
**Category framing:** Computational Biology & Bioinformatics (CBIO), with Mathematics (MATH) as the strong alternate because the load-bearing novelty is a closed-form spectral-perturbation result.
**Constraint honored throughout:** 100% in-silico. Zero patients, zero IRB (only a documented SRC human-data *exemption* determination for the public de-identified geometries, filed before any experiment).

---

## 1. The one sentence

The closed-form sensitivity of atrial algebraic connectivity to edge uncoupling,

> **dλ₂/dw_ij = (φ_i − φ_j)²**   (φ = Fiedler vector)

defines a per-region **Spectral Fragility Index (SFI)** — the expected drop in atrial algebraic connectivity per unit perioperative uncoupling — that predicts where a fibrotic atrium becomes reentry-inducible after diffuse perioperative stress, validated entirely against open cardiac meshes + free openCARP simulations.

## 2. The protected novel seed (survives every de-scope)

The **only** thing that must survive is:

- The closed-form Fiedler-value sensitivity `dλ₂/dw_ij = (φ_i − φ_j)²`,
- aggregated under a **stochastic, fibrosis-weighted, diffuse perioperative-uncoupling field** Δw (NOT literal surgical cuts),
- defined as a per-region SFI,
- with reentry-initiation hotspots hypothesized where `|∇φ₂|` co-localizes with the Perron (dominant-activation-mode) eigenvector,
- and a **spectral-projector / subspace** generalization that stays well-defined when λ₂ is near-degenerate.

Everything else (openCARP labels, ML classifier, cohort, dashboard, UQ) is a known, published pipeline. The contribution is **one new interpretable feature inside an established pipeline** — position it that way, never as a new paradigm.

## 3. Why the mechanism is honest

CABG does **not** incise the left atrium. POAF after CABG is driven by diffuse inflammatory edema, transient ischemia, autonomic APD/refractoriness heterogeneity, and atrial stretch. So we model those as a distribution of **small, diffuse edge-weight reductions** — which is *exactly* the small-‖ΔL‖ regime where first-order eigenvalue perturbation theory is valid. Physiology and math align instead of fighting. The literal "incision = edge deletion" idea is reserved only for surgeries where it is physiologically true (Maze, atriotomy, transplant), and for those we recompute λ₂ **exactly** rather than using the linear formula.

---

## 4. The killer objections and the required defenses

A hostile ISEF/EP/CFD judge will raise these. Each must be pre-answered **in the build**, not just in talk.

| # | Objection | Required defense (build it in) |
|---|-----------|-------------------------------|
| 1 | **Thin novelty** — the openCARP-cohort → feature → inducibility-classifier pipeline is already published (Roney 2022; 2025 fibrosis-heterogeneity inducibility paper). "You renamed min-cut" (Zahid/Trayanova 2016 is the spectral cousin of min-cut via Cheeger). | Cite those papers **up front**. Run SFI as an **add-on/ablation** to their published fibrosis-heterogeneity feature set on the **same** openCARP labels. Pre-register the falsifiable hypothesis that SFI adds incremental AUC / localization beyond fibrosis entropy + patch-size. |
| 2 | **Circularity** — edge weights are reduced in fibrotic nodes AND the openCARP label is driven by that same fibrosis field, so SFI may just re-encode fibrosis burden. | Real competitor baselines (below), DeLong test on shape-family-held-out folds. A clean null ("SFI ≡ re-encoded fibrosis") is an **accepted, reportable** outcome — say so out loud. |
| 3 | **Linear-spectral vs nonlinear-excitable-media mismatch** — λ₂/Fiedler/Perron are linear diffusion/centrality quantities; AF reentry is governed by wavelength (CV×ERP), source-sink mismatch, unidirectional block. "|∇φ₂| ∩ Perron predicts reentry origin" is asserted, not proven. "Spectral radius = reentry" has no mechanistic basis. | Turn the claim into a **falsifiable spatial test**: colocalize `|∇φ₂| ∩ Perron` hotspots with openCARP reentry origins vs a **UAC rotational/shift spatial-null**, report an overlap statistic with CI. Keep or **delete** each linear-spectral claim by its measured correlation. Drop "spectral radius = reentry" unless it empirically correlates. |
| 4 | **Over-claimed ground truth** — no real POAF patient anywhere; the label is a monodomain sim of synthetic anatomy under an invented uncoupling field. Calling it a "POAF predictor" / "in-silico clinical trial per FDA/VICTRE" is unsupported. | Rescope claims to "a cheap spectral surrogate for a **specific EP-simulator's inducibility verdict**." **Delete all FDA/VICTRE/ISCT language.** Calibrate Δw to literature-measured perioperative CV/ERP changes; report sensitivity to that choice. State explicitly it is unvalidated against real POAF. Use the 2018 LA Segmentation Challenge geometries **only** as an anatomy-realism check, never as clinical validation. |
| 5 | **Compute-infeasible UQ** — full Sobol over 4 params × openCARP cost × cohort is thousands of GPU-hours. | Descope UQ to **Morris elementary-effects screening + a Gaussian-process/PCE surrogate** trained on a few dozen openCARP runs. Benchmark real wall-clock on the target free GPU **before** fixing cohort size. Reuse published Roney-cohort labels where possible. |
| 6 | **Perturbation-theory misuse** — first/second-order expansion is valid only for ‖ΔL‖ small vs the spectral gap; Maze/discrete cuts break it. "Bounded by Weyl" is the wrong tool for the eigenvector term. | State the **validity radius ‖ΔL‖ ≪ λ₃−λ₂ explicitly** and enforce it as a runtime guard. Use **Weyl only for eigenvalue bounds**; use **Davis–Kahan / Bauer–Fike for the eigenvector/subspace term**. Recompute λ₂ exactly for any discrete cut. |
| 7 | **Near-degenerate λ₂ + synthetic-N leakage** — small λ₂−λ₃ gap makes the single Fiedler vector ill-defined; PCA-resampled anatomies are near-duplicates that leak train↔test and inflate AUC. | Auto-switch to the **spectral-projector (subspace) SFI** when the gap is small; report the λ₂−λ₃ gap distribution across the cohort. **GroupKFold by shape family**, assign each PCA-resampled anatomy to its parent's group; report **both grouped and naive AUC** to expose leakage. |

---

## 5. The applied-math core (own every symbol in the interview)

- Weighted Laplacian `L = D − W`, with W = fibre-anisotropic conduction couplings (fast along-fibre, slow cross-fibre), reduced in fibrotic/edematous nodes.
- Because `L = Σ_(i,j) w_ij (e_i − e_j)(e_i − e_j)ᵀ`, the first-order sensitivity of any simple eigenvalue λ with unit eigenvector φ is `dλ/dw_ij = φᵀ(∂L/∂w_ij)φ = (φ_i − φ_j)²`.
- Applied to λ₂: **`dλ₂/dw_ij = (φ₂,i − φ₂,j)²`** — a per-edge fragility.
- Region aggregate under a stress distribution:
  `SFI(R) = E[Δλ₂] ≈ Σ_(i,j)∈R E[Δw_ij](φ_i − φ_j)² + Σ_{k≠2} (ψ_kᵀ ΔL φ₂)²/(λ₂ − λ_k)`
  — the second term is the second-order (large-edit) resolvent correction, bounded by Weyl's inequalities for the eigenvalue part.
- **Cheeger's inequality** links λ₂ to conductance, tying SFI to substrate isolation without collapsing to a single min-cut.
- **Perron–Frobenius** gives the dominant eigenpair of the nonnegative adjacency; the Perron vector's inverse-participation-ratio flags where activation localizes.
- Hotspots predicted where `|∇φ₂|` and Perron localization overlap.

---

## 6. The number to beat (falsifiable headline)

> On shape-family-held-out (GroupKFold) folds, adding the per-region SFI to the published fibrosis-heterogeneity feature set raises openCARP-inducibility classification AUC by **ΔAUC ≥ 0.05 with DeLong p < 0.05**, while the same SFI hotspot map localizes the simulator's reentry-initiation site with **median geodesic error below the rotational spatial-null's 5th percentile**.

---

## 7. Verified open datasets & free tools (no hospital)

| Resource | Role | Access |
|----------|------|--------|
| Rodero et al. **1000 synthetic four-chamber meshes** | Anchor anatomical cohort (PCA-SSM, no real person) | Zenodo **4506930** (CC) |
| Roney et al. **LA virtual-cohort meshes + fibrosis + UAC + fibre fields** | Weighted-graph substrate + heterogeneity + reproducible line placement | Zenodo **5801337**, DOI 10.5281/zenodo.5801337; paper Circ AE 2022 DOI 10.1161/circep.121.010253 |
| Roney et al. **Human Atrial Fibre Atlas** (7 ex-vivo DT-MRI) | Anisotropy for physiologically correct edge weights; CARP `.pts/.elem/.lon` | Zenodo **3764917**; Ann Biomed Eng 2020 DOI 10.1007/s10439-020-02525-w |
| **atrialmtk** (pcmlab) | Turns SSM anatomies into bilayer/volumetric meshes with fibres/regions; bundles openCARP + meshtool | GitHub pcmlab/atrialmtk (GPL-3.0); Zenodo 10139306 |
| **openCARP** | Monodomain/eikonal EP solver → synthetic ground-truth inducibility labels; ships Niederer N-version benchmark | opencarp.org (free); benchmark PMID 21969679 |
| **2018 LA Segmentation Challenge** (Cardiac Atlas Project) | Real, ethics-cleared, public LA geometries — anatomy-realism check only | atriaseg2018.cardiacatlas.org |
| **Vascular Model Repository + SimVascular** | Independent geometries for generalization tests; optional low-WSS CFD arm | vascularmodel.org / simvascular.github.io (BSD) |

---

## 8. The one-season methodology (as pre-registered, descoped to feasible)

1. **Substrate (laptop, free).** Public 2018 LA Segmentation geometries + openCARP example meshes as realism anchors. Build ~8–12 base atria; morph each into a few statistical-shape variants → ~40–80 graphs. Paint fibrosis from LGE proxy + parametric patchy fields. **Tag every variant with its base shape family for grouped CV.**
2. **Spectral engine (seconds).** Weighted Laplacian, conductance reduced in fibrotic nodes (calibrated to literature CV ranges). Compute λ₂, φ₂, λ₃ via sparse `eigsh`.
3. **Protected seed — SFI.** `dλ₂/dw_ij=(φ₂_i−φ₂_j)²`; diffuse fibrosis-weighted stochastic uncoupling field Δw calibrated to literature perioperative CV/ERP changes; aggregate to per-region SFI (analytic expectation + Monte Carlo). Enforce validity radius `‖ΔL‖ ≪ λ₃−λ₂`; auto-switch to subspace SFI when the gap is small; recompute λ₂ exactly for any discrete-cut comparison.
4. **Ground truth (descoped openCARP).** Monodomain + phenomenological Mitchell-Schaeffer ionics on coarsened meshes; S1–S2/burst from a handful of pacing sites, pre- and post-uncoupling. Extract binary inducible/not labels + reentry-origin coordinates. Explicitly a **simulator verdict, not clinical POAF**.
5. **Real competitors (kill the strawman).** Fibrosis burden, spatial entropy/patch-size, deterministic min-cut, percolation threshold, **λ₂-alone**.
6. **Primary test.** Nested GroupKFold (shape-family-held-out) inducibility classifier; ΔAUC of adding SFI via **DeLong**; report grouped AND naive AUC.
7. **Mechanism test (falsifiable).** Colocalize `|∇φ₂| ∩ Perron` hotspots with openCARP reentry origins vs a UAC rotational/shift spatial-null (with CI); keep/delete each claim by measured correlation.
8. **UQ (descoped).** Morris elementary-effects + a GP surrogate on the uncoupling-field parameters over a few dozen runs — not full Sobol.

---

## 9. Validation without a single patient

- **Code verification:** reproduce the Niederer N-version cardiac-EP benchmark (ships with openCARP) — verifies the engine that makes your labels.
- **Analytic graph verification (hard CI gate):** path/ring/grid Laplacian eigenvalues `λ_k = 2 − 2cos(kπ/N)`; known single-edge-cut → analytic Δλ₂. Validates the spectral pipeline independent of any data.
- **Held-out synthetic ground truth:** openCARP inducibility label is the reference; ROC-AUC, sensitivity/specificity on never-seen virtual subjects with proper CV.
- **Convergence / mesh-independence:** λ₂ and predicted labels stable under mesh refinement + time-step reduction (ASME V&V-40 calculation verification).
- **Literature range-checking:** POAF incidence ~26% pooled after CABG; LA conduction velocities ~0.3–1.2 m/s; predictive power comparable to real shape-based predictors (Bieging LAA-shape+CHADS-VASc AUC 0.640→0.778; Varela LA-shape AUC ~0.68–0.71).
- **Anatomical realism check:** SSM synthetic LA shape statistics overlap the public 2018 Atrial Segmentation geometries.
- **Independent-geometry generalization:** re-run on unrelated open geometries (VMR) to show it isn't overfit to one shape model.

---

## 10. Deliverables & figures

Open **AtrialSpectralBench** repo + dataset; SFI theory note with the perturbation derivation; interactive dashboard; figures: (1) φ₂-colored atria, (2) SFI hotspot maps vs openCARP reentry origins, (3) eigenvalue spectra, (4) ROC vs clinical baseline vs full-EP gold standard, (5) Sobol/Morris sensitivity bars.

---

## 11. Division of labor — Claude Code vs Claude Science

### ⚙ Claude Code BUILDS the engine (all software, unit-tested — no heavy compute)
- Scaffold Python project (NumPy/SciPy-sparse, NetworkX/igraph, scikit-learn, optional PyTorch Geometric).
- Graph construction: weighted adjacency W from mesh geodesic connectivity + fibre anisotropy + fibrosis reduction; degree matrix D.
- Laplacian variants: `L = D − W`, symmetric-normalized `L_sym`, random-walk normalized.
- Spectrum: `scipy.sparse.linalg.eigsh` → Fiedler value/vector; spectral gap; near-zero-eigenvalue count.
- Perron–Frobenius: dominant eigenvector + spectral radius.
- Spectral-feature extractor: spectral gap, eigenvalue distribution, spectral entropy, heat-kernel signatures, Cheeger-constant estimate.
- **The SFI module** (the novel seed): closed-form derivative, stochastic Δw field, analytic + Monte-Carlo aggregation, validity-radius guard, subspace/projector fallback.
- openCARP **runner code** (config, coarsening, pacing protocol, label + reentry-origin extractor) — built and unit-tested on tiny/mock meshes.
- Competitor baselines (fibrosis burden, entropy, min-cut, percolation, λ₂-alone).
- Evaluation stack: nested GroupKFold, DeLong, colocalization vs spatial-null, ROC/AUC, calibration, bootstrap CIs, permutation tests, SHAP.
- Property tests (Laplacian PSD + zero row-sum; Fiedler on path/cycle; Perron positivity; analytic-spectrum CI gates).
- Streamlit/Plotly dashboard; config-driven reproducible seeded repo; Dockerfile; CI.

### 🔬 Claude Science RUNS the experiments (heavy compute + stats + figures) — **DEFERRED**
- Execute the openCARP monodomain inducibility sweep to generate real ground-truth labels.
- Run the spectral robustness sweeps; edge/node-removal Fiedler sweeps.
- Statistical validation: bootstrap + permutation tests, Bayesian hierarchical comparison across graph families with credible intervals.
- Morris/GP UQ over the uncoupling-field parameters.
- Generate the publication figures (Fiedler heatmaps, eigenvalue spectra, ROC panels, partition visualizations) and `results_report.md`.
- Literature synthesis grounding the findings.

---

## 12. Sequenced first steps — front-load ALL Claude Code work

> **Governing decision:** everything Claude Code can build and unit-test **without** running the expensive openCARP sweep comes first. The engine is written so that the moment Claude Science access arrives, it is one config-run away from generating labels and figures. Analytic/mock tests stand in for real labels until then.

**Phase 0 — Scaffold & environment (Claude Code)**
- `pyproject.toml`/conda env with pinned NumPy, SciPy, NetworkX/igraph, scikit-learn, PyVista, matplotlib; seed control; Dockerfile; GitHub Actions CI; dated-notebook stub; README with the pre-registered hypothesis written down.
- File the SRC / human-data **exemption determination** paperwork now (before any experiment) — the 2018 LA Challenge + public Zenodo data need a documented exemption, not full IRB.

**Phase 1 — Spectral engine + analytic CI gates (Claude Code, NO data needed)**
- Laplacian builders (`L`, `L_sym`, random-walk).
- `eigsh` spectrum → λ₂, φ₂, λ₃, spectral gap, Perron pair.
- **Hard CI gate:** reproduce `λ_k = 2 − 2cos(kπ/N)` for path/ring/grid graphs; known single-cut Δλ₂. *(This is the week-1–2 credibility gate; it needs no dataset and no simulator.)*

**Phase 2 — The protected novel seed: SFI (Claude Code, NO data needed)**
- Implement `dλ₂/dw_ij = (φ₂_i − φ₂_j)²`.
- Stochastic fibrosis-weighted diffuse Δw field; analytic-expectation + Monte-Carlo per-region aggregation.
- Validity-radius guard (`‖ΔL‖ ≪ λ₃−λ₂`); subspace/projector SFI fallback for near-degenerate λ₂; exact-λ₂ recompute path for discrete cuts.
- Unit tests on small synthetic graphs where the answer is analytically known.

**Phase 3 — Data pipeline + one real atrium (Claude Code)**
- Download Rodero 1000-mesh cohort (Zenodo 4506930); run **atrialmtk** to extract one LA surface graph with fibres + Universal Atrial Coordinates.
- Assemble weighted graph: edge weights from fibre anisotropy, reduced in fibrotic nodes sampled from Zenodo 5801337.
- Compute λ₂, φ₂, and the SFI hotspot map **on that one atrium** — this is the entire novel seed running in seconds on a laptop. *(Milestone: the first real figure.)*

**Phase 4 — openCARP integration code (Claude Code, built NOT run at scale)**
- Runner: mesh coarsening, Mitchell-Schaeffer monodomain config, S1–S2/burst pacing protocol, inducibility + reentry-origin label extractor.
- Validate the runner on a **tiny** example mesh so the code path is proven; the full sweep is deferred to Claude Science.
- Wire the Niederer N-version benchmark reproduction as a runnable-but-not-yet-run verification target.

**Phase 5 — Competitor baselines + cohort generator (Claude Code)**
- Fibrosis burden, spatial entropy/patch-size, deterministic min-cut, percolation threshold, λ₂-alone.
- Cohort builder: 8–12 base atria → statistical-shape variants → ~40–80 graphs, each tagged with its **base shape family** for grouped CV.

**Phase 6 — Evaluation + ML stack (Claude Code, tested on mock labels)**
- Nested GroupKFold (shape-family-held-out), DeLong ΔAUC, colocalization vs UAC rotational spatial-null, ROC/calibration/bootstrap CIs, SHAP.
- Run the whole pipeline end-to-end on **mock/placeholder labels** so it is proven and one real-label run away from a result.

**Phase 7 — Dashboard, reproducibility, property tests (Claude Code)**
- Streamlit/Plotly dashboard: atrial graph colored by φ₂ + per-subject SFI report.
- Property-test suite (Laplacian PSD, zero row-sum, Perron positivity, analytic-spectrum gates) as CI.
- Freeze configs/seeds; write the reproducible run instructions.

**➡ HANDOFF POINT — everything above needs no Claude Science.**

**Phase 8 — Claude Science (DEFERRED until access):** run the openCARP inducibility sweep for real labels → run the primary GroupKFold+DeLong test → mechanism colocalization test → Morris/GP UQ → generate final figures + `results_report.md` → Bayesian hierarchical validation.

---

## 13. The universal sequence (every idea)

**Verify → Reproduce → Innovate → Validate → Quantify uncertainty → Write → Rehearse the interview.**

- Reproduce a published baseline BEFORE innovating (proves the pipeline, reveals the real gap).
- File the SRC/human-data exemption BEFORE any experiment.
- Keep a dated notebook (judges inspect it).
- Budget the last month for the interview — it is 25 of 100 points; Creativity+Impact and Presentation together are 55%.
