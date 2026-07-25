# The Spectral Fragility Index — theory note

The load-bearing mathematics of AtrialSpectralBench, with the validity boundary made
explicit and the prior art credited up front. Nothing here is asserted without its
proof or its citation; the empirical companion is GM3 (`results/gm3_metrics.json`).

## 0. Credits (this is not new mathematics — the *use* is the contribution)

- The identity `∂λ₂/∂w_ij = (φ_i − φ_j)²` is a **classical** eigenvalue-perturbation
  result and is the engine of Fiedler-vector edge-editing heuristics — **Ghosh & Boyd,
  "Growing well-connected graphs," IEEE CDC 2006**; algebraic connectivity itself is
  **Fiedler, Czech. Math. J. 1973**.
- That `λ₂` governs the **stability of coupled excitable/oscillator networks** is the
  **master-stability-function** result — **Pecora & Carroll, PRL 1998** — with the
  synchronizability eigenratio `λ_N/λ₂` from **Barahona & Pecora, PRL 2002**; the Fiedler
  value is even used as a real-time seizure detector (**Bomela et al., Sci. Rep. 2020**).

**What is ours** is (i) aggregating the per-edge derivative under a *stochastic,
fibrosis-weighted, diffuse* perioperative-uncoupling field into a per-region index; (ii)
testing it as a leakage-controlled *add-on* against real competitors and a spatial null;
(iii) mapping its **validity radius** empirically and matching it to theory; and (iv)
showing which parts *transfer* to a second excitable medium.

## 1. The weighted conduction graph

Nodes = mesh vertices; edges = triangle edges; edge weight `w_ij > 0` = a fibre-anisotropic
conduction conductance, attenuated by local fibrosis. The combinatorial Laplacian is

```
L = D − W ,   D = diag(Σ_j w_ij) .
```

`L` is symmetric positive-semidefinite with a single zero eigenvalue (constant vector) for
a connected graph. Its eigenpairs `0 = λ₁ < λ₂ ≤ λ₃ ≤ …` with eigenvectors `φ₁, φ₂, …`;
`λ₂` is the **algebraic connectivity** and `φ₂` the **Fiedler vector**.

## 2. The closed-form fragility (first order)

`L` is a sum of rank-one edge terms:

```
L = Σ_{(i,j)∈E} w_ij (e_i − e_j)(e_i − e_j)ᵀ  ⇒  ∂L/∂w_ij = (e_i − e_j)(e_i − e_j)ᵀ .
```

For a **simple** eigenvalue `λ` with unit eigenvector `φ`, first-order perturbation theory
(Hellmann–Feynman) gives

```
∂λ/∂w_ij = φᵀ (∂L/∂w_ij) φ = (φ_i − φ_j)² .
```

Applied to the Fiedler pair:

> **`∂λ₂/∂w_ij = (φ₂,ᵢ − φ₂,ⱼ)²`**  — a non-negative **per-edge fragility**.

Reducing a conductance `w_ij` by a small `Δw_ij ≥ 0` (uncoupling) therefore **lowers**
`λ₂` by `Δw_ij (φ₂,ᵢ − φ₂,ⱼ)²` to first order. (`asb.sfi.edge_fragility`.)

> **Machine-verified** (`scripts/verify_sfi_identity.py`, sympy + mpmath at 40 digits): the
> first-order edge identity `φᵀ(e_i−e_j)(e_i−e_j)ᵀφ = (φ_i−φ_j)²` and the Laplacian edge
> structure `∂L/∂w_ij = (e_i−e_j)(e_i−e_j)ᵀ` are proved symbolically; on a concrete 5-node graph
> the finite-difference `∂λ₂/∂w_ij` matches `(φ₂,ᵢ−φ₂,ⱼ)²` to `|err| ≈ 2×10⁻¹⁷` and the §4
> second-order resolvent term to `≈ 5×10⁻¹⁷`. The derivation is not asserted; it is checked.
> **Six independent guarantees** (adding `scripts/verify_math_rigor.py`): (P4) the identity holds
> *exactly* with rational weights (zero floating point); (P5) an **interval-arithmetic certified
> enclosure** proves the true eigenvalue is within `~4×10⁻⁴⁰` of the computed `λ₂`; (P6) `λ₂`/`φ₂`
> agree across **scipy + networkx + a hand-rolled inverse iteration** (spread `2.5×10⁻¹⁵`), ruling
> out an implementation bug in any single solver.
>
> **Kernel-checked in a proof assistant** (`formal/sfi_edge_identity_Q.v`, Coq 8.18): the algebraic
> core is re-proved inside Coq and checked by its small trusted kernel — a level of assurance
> above a computer-algebra system. Theorems: **T1** the rank-one quadratic form `φᵀ(vvᵀ)φ = (v·φ)²`;
> **T2** the edge vector `(e_i−e_j)·φ = φ_i−φ_j`; **T3** the SFI edge sensitivity `φᵀE_ijφ = (φ_i−φ_j)²`;
> **T4/T5** each SFI term and the aggregate are non-negative (the SFI is a genuine Dirichlet energy;
> the Laplacian form is PSD). The proof is over the **exact rationals `Q`**, so `Print Assumptions`
> reports *"Closed under the global context"* — **zero axioms** (the R-valued companion
> `formal/sfi_edge_identity.v` is identical mathematically but inherits Coq's two standard
> real-number axioms). The analytic Hellmann–Feynman step (differentiability of a simple eigenvalue)
> is the classical cited result (Ghosh–Boyd 2006); it is *not* re-proved here and is not claimed to be.
> Lean/Mathlib was attempted but its toolchain binaries are GitHub release assets blocked by the
> compute environment's egress policy; Coq provides the equivalent kernel-checked guarantee.

## 3. Per-region aggregation → the SFI

Model diffuse perioperative stress as a distribution of small, non-negative, fibrosis-
weighted fractional reductions `Δw` (mean `E[Δw_ij] = frac_ij · w_ij`, `frac_ij` growing with
local fibrosis; calibrated to literature perioperative CV slowing, **Δw ≈ 0.36**, frozen in
`docs/PRE_REGISTRATION.md`). The **per-region Spectral Fragility Index** is the expected
`λ₂` drop contributed by region `R`:

```
SFI(R) = E[Δλ₂ from R] ≈ Σ_{(i,j)∈R} E[Δw_ij] (φ₂,ᵢ − φ₂,ⱼ)²
                          + Σ_{k≠2} (φ_kᵀ ΔL φ₂)² / (λ₂ − λ_k)  +  O(‖ΔL‖³).
```

The first term is `asb.sfi.sfi_region` (analytic expectation); `asb.sfi.sfi_monte_carlo`
recomputes `λ₂` exactly per sampled field as the ground-truth check. The second term is the
**second-order resolvent correction** — the source of the validity radius below.

## 4. The validity radius (the honest core)

The first-order form is trustworthy only while the perturbation is small **relative to the
spectral gap** `λ₃ − λ₂`. Two different tools bound two different things (a common error is
to use Weyl for both):

- **Weyl's inequality** bounds the **eigenvalue** move: `|Δλ₂| ≤ ‖ΔL‖`. (Verified to hold at
  every ρ in GM3.)
- **Davis–Kahan / Bauer–Fike** bound the **eigenvector/subspace** move, which scales like
  `‖ΔL‖ / (λ₃ − λ₂)`. This — not Weyl — governs when `φ₂` (hence the single-vector SFI)
  degrades.

Define `ρ = ‖ΔL‖ / (λ₃ − λ₂)`. The relative error of the first-order Δλ₂ prediction grows
`~ ρ`. **Empirically (GM3): the first-order error crosses 10% at ρ\* ≈ 3**, is 1% at ρ = 0.3,
and 92% at ρ = 28 — monotone in ρ, with Weyl holding throughout. `asb.sfi.validity_radius_ok`
enforces `‖ΔL‖ ≤ safety · (λ₃ − λ₂)` as a runtime guard.

**On real excitable media the gap is tiny**, so ρ is enormous: atrial cohort median
ρ ≈ 2422 — ~10³× past ρ\*. The single-vector linear SFI is thus provably outside its validity
radius on real anatomy: by Davis–Kahan `φ₂` is ill-conditioned and the first-order Δλ₂ *magnitude*
is unreliable. **This does NOT by itself explain the predictive null.** ρ bounds estimator accuracy,
not the classification content of a monotone ranking feature; and the **exact** Δλ₂ SFI (no ρ
limitation, §6) is *also* non-predictive (GM1), so the null is a *feature-redundancy* result —
SFI's information is already carried by the standard connectivity features — which is label-dependent
and separate from ρ. Two honest caveats on ρ itself: its numerator `‖ΔL‖` carries the
fibrosis-weighted Δw (so only the denominator `λ₃−λ₂` is label-free), and ρ\*≈3 is a path-graph
constant, so "10³× past" is order-of-magnitude. The label-free content of this section is precisely:
*the single-vector estimator is numerically invalid on real tissue* — a real, portable a-priori
screen — not "the biomarker provably cannot predict."

## 5. Near-degeneracy → the subspace (projector) SFI

When `λ₂ ≈ λ₃` the individual `φ₂` is an arbitrary rotation within the low invariant subspace,
so `edge_fragility` is ill-conditioned. The **trace over the k-dimensional invariant subspace**
`P = Σ_c φ_c φ_cᵀ` is basis-independent and well-defined:

```
s_ijᵀ P s_ij = Σ_c (φ_c,ᵢ − φ_c,ⱼ)²  ,   s_ij = e_i − e_j ,
```

the first-order contribution to the drop in the **sum** of the clustered eigenvalues
(`asb.sfi.subspace_sfi`). GM3 confirms the subspace prediction stays accurate as the gap
closes where the single vector is erratic.

## 6. Discrete cuts → exact recompute

For genuinely discrete edits (Maze lesion, atriotomy) `‖ΔL‖` is large and the expansion is
abandoned: `λ₂` is recomputed **exactly** after the rank-one edit
(`asb.sfi.exact_delta_lambda2`). This is the honest boundary of the linear theory and the
reference for the derivative-identity finite-difference test.

## 7. The localizer

Reentry/instability **initiation** was hypothesized at `|∇φ₂| ∩ Perron` (`asb.sfi.hotspot_map`,
with `|∇φ₂|` the edge-averaged Fiedler gradient and the Perron vector the dominant nonnegative
adjacency eigenpair). GM2/GM4 **kept `|∇φ₂|`** as a weak **rank** localizer (tie-robust
permutation p < 1e-3 in both cardiac and neural media) while the strict rotational-spatial-null
endpoint is **null for every field**, and **deleted the Perron/∩ and centrality claims**
(they do not correlate). The keep/delete was decided by measured correlation, never assumed.

## 8. Supporting classical facts (used, not claimed as new)

- **Cheeger's inequality** `λ₂/2 ≤ h ≤ √(2λ₂)` ties `λ₂` to conductance (substrate isolation)
  without collapsing to a single deterministic min-cut. (`asb.spectral.cheeger_estimate`.)
- **Perron–Frobenius** gives the nonnegative dominant adjacency eigenpair; its
  inverse-participation-ratio flags activation localization. (`asb.spectral.perron`.)
- Analytic gates: path/ring/grid Laplacian eigenvalues `λ_k = 2 − 2cos(kπ/N)` and the
  known single-edge-cut Δλ₂ — hard CI checks that validate the pipeline independent of data.

## 9. The validity-radius criterion (stated, with derivation)

Section 4 reported the validity radius empirically (ρ\* ≈ 3). Here is the statement made
explicit. **Read §9.4 first: the mathematics below is entirely classical — this is not a new
theorem, and must never be presented as one.**

### 9.1 Setup
`L = D − W` symmetric PSD with `0 = λ₁ < λ₂ ≤ λ₃ ≤ …`, Fiedler pair `(λ₂, φ₂)`, `φ₂` simple,
unit-norm. `ΔL` a symmetric perturbation (our fibrosis-weighted uncoupling field is itself a
Laplacian, so `ΔL` is symmetric PSD). The linear SFI uses `Δλ₂⁽¹⁾ = φ₂ᵀ ΔL φ₂ =
Σ_{(i,j)} Δw_ij (φ₂,ᵢ − φ₂,ⱼ)²`.

### 9.2 Second-order truth (Rayleigh–Schrödinger)
For a simple eigenvalue,

```
Δλ₂ = φ₂ᵀΔL φ₂  +  Σ_{k≠2} |φ₂ᵀ ΔL φ_k|² / (λ₂ − λ_k)  +  O(‖ΔL‖³).
```

The `k = 3` term dominates (smallest denominator); it is negative and bounded by
`|φ₂ᵀΔLφ₃|²/(λ₃−λ₂) ≤ ‖ΔL‖²/(λ₃−λ₂)` (unit vectors ⇒ `|φ₂ᵀΔLφ_k| ≤ ‖ΔL‖`).

### 9.3 The criterion
Writing the leading correction `S₂` and `Δλ₂⁽¹⁾ ≍ c·‖ΔL‖` (`c ∈ (0,1]` = alignment of `ΔL`
with the `φ₂` direction), the **relative error of the first-order SFI** is

```
|Δλ₂ − Δλ₂⁽¹⁾| / |Δλ₂⁽¹⁾|  ≈  |S₂|/|Δλ₂⁽¹⁾|  ≲  (1/c)·‖ΔL‖/(λ₃−λ₂)  =  (1/c)·ρ ,
      ρ := ‖ΔL‖ / (λ₃ − λ₂).
```

The relative error is **`Θ(ρ)`** — linear in perturbation-size / spectral-gap. GM3 measures
exactly this: 1 % at ρ = 0.3, ≈10 % at ρ\* ≈ 3, 92 % at ρ = 28. The constant `c` is
problem-dependent, which is **why ρ\* is calibrated empirically, not asserted**.

**Eigenvector corollary (Davis–Kahan).** `sin∠(φ₂, φ₂′) ≤ √2·‖ΔL‖/(λ₃−λ₂) = √2·ρ`. For
`ρ ≳ 1`, `φ₂` rotates into the `φ₃` direction and every `φ₂`-derived quantity (edge fragility,
`|∇φ₂|`, single-vector SFI) is ill-posed — precisely the regime where §5's basis-independent
subspace SFI is the only well-defined object.

**In one sentence:** *a Fiedler-based spectral perturbation biomarker is trustworthy only while
`ρ = ‖ΔL‖/(λ₃−λ₂) = O(1)`.* Real excitable media have a tiny gap, so ρ is enormous —
**atrial-cohort median ρ ≈ 2422, neural-cohort median ρ ≈ 1852**, ~10³× past ρ\*. This proves
the single-vector estimator is numerically invalid on real tissue; it does **not** explain the
GM1 predictive null, which is separately a feature-redundancy result (§4 — the exact Δλ₂ SFI,
immune to ρ, is also non-predictive). It *is* a **screening test any spectral fragility
biomarker can be subjected to** in any domain (seizure-focus localization, connectomics, power grids): compute ρ; if `ρ ≫ 1`,
distrust the linear spectral biomarker regardless of in-sample correlation.

### 9.4 Novelty scope (state this honestly, always)
The mathematics in §9.1–9.3 — Weyl's inequality, Davis–Kahan, second-order Rayleigh–Schrödinger
perturbation, and `ρ` governing when perturbation theory is valid — is **classical numerical
linear algebra (1970s and earlier). Nothing here is a new theorem.** Any applied mathematician
knows perturbation theory fails once the perturbation is comparable to the eigenvalue gap. The
**contribution is not mathematical**; it is: (i) recognizing that this exact ratio is the right
*a priori* validity check for a spectral fragility biomarker; (ii) *measuring* that real atrial
and neural media violate it by ~10³; and (iii) using it to *mechanistically explain* a negative
result instead of merely reporting a null. That is a methodological / framing contribution.
Presenting §9 as "a new theorem" would be false and would (rightly) collapse under expert
questioning; presenting it as "applying a classical criterion as a biomarker screening test,
and quantifying the violation on real tissue" is defensible and still genuinely useful.

## 10. The ρ scaling law (why the failure is not incidental)

§9 is a per-instance criterion. It has a *scaling* consequence that is the honest quantitative
contribution. The gap `λ₃ − λ₂` of a diffusively-coupled medium collapses with system size by
**Weyl's law**: for a `d`-dimensional domain the low Laplacian eigenvalues scale as
`λ_k ~ (k/N)^{2/d}`, so the gap shrinks with `N` and therefore

```
ρ = ‖ΔL‖ / (λ₃ − λ₂)  grows lawfully with resolution.
```

**Measured** (`results/rho_scaling.json`, gap-vs-`N` at fixed mean degree):

| medium | fitted gap exponent | reference |
|---|---|---|
| real atrial surface (2-manifold) | `N^{-1.05}` | Weyl 2-manifold `−1.0` (matched to 5 %) |
| 2-D random-geometric | `N^{-1.81}` | (faster than naive `2/d`) |
| 3-D random-geometric | `N^{-0.84}` | (faster than naive `2/d`) |

So the biomarker does not fail by accident on one dataset — **it fails *predictably*, and *more* as the
medium is meshed finer** (the opposite of "more resolution helps"). The real 2-manifold matches Weyl's
`−1.0`; the abstract graphs collapse even faster, so `ρ` grows *at least* as fast as Weyl predicts.

**Honest scope.** Weyl's law is classical. What is ours is (i) tying it to *biomarker validity* rather
than PDE discretization error, (ii) measuring the gap-scaling exponent across independent media, and
(iii) the counterintuitive corollary that higher-fidelity meshes are *further* past the validity radius.
Only the real manifold matches the Weyl exponent exactly; we do not claim more.

## 11. The falsification protocol (what actually caught the false positive)

The guards are not decoration; each one, removed, manufactures a specific false positive
(`results/falsification_protocol.json`):

- **Shape-family GroupKFold** — remove it and shape-family **leakage inflates AUC by +0.117**
  (naive random-CV 0.756 vs grouped 0.640 on a controlled replicated cohort).
- **Effect-size gate** — remove it and the **p-value trap** fires: subspace-SFI is `p ≈ 2×10⁻²⁴` at
  large `N`, which *looks* like a hit, but `ΔAUC = +0.003` needs `~4150` cases to detect at 80 % power
  (clinical cohorts are ~200–1000) — statistically significant, practically useless.
- **Spatial null (torus rotation) + keep/delete by measured correlation** — remove it and
  spatially-smooth centrality fields "localize" the origin spuriously.

**The claim.** A naive analyst (random CV + `p<0.05`, no effect-size gate, no spatial null) would have
reported SFI as a *working* reentry biomarker. The pre-registered protocol correctly rejects it. This
is the packageable methods contribution — a worked demonstration that standard biomarker methodology
produces exactly the false positives this protocol is built to catch, on a candidate plausible enough
to fool it.
