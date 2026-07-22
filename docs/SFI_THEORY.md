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
ρ ≈ 2422, neural-network cohort median ρ ≈ 1852 — both ~10³× past ρ\*. The linear SFI is thus
provably outside its validity radius on real anatomy, which is *why* it collapses to a
substrate re-encoding as a predictor (GM1 null). This is the project's central, honest result.

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
adjacency eigenpair). GM2/GM4 **kept `|∇φ₂|`** (localizes the origin above a spatial null in
both cardiac and neural media, p < 1e-3) and **deleted the Perron/∩ and centrality claims**
(they do not correlate). The keep/delete was decided by measured correlation, never assumed.

## 8. Supporting classical facts (used, not claimed as new)

- **Cheeger's inequality** `λ₂/2 ≤ h ≤ √(2λ₂)` ties `λ₂` to conductance (substrate isolation)
  without collapsing to a single deterministic min-cut. (`asb.spectral.cheeger_estimate`.)
- **Perron–Frobenius** gives the nonnegative dominant adjacency eigenpair; its
  inverse-participation-ratio flags activation localization. (`asb.spectral.perron`.)
- Analytic gates: path/ring/grid Laplacian eigenvalues `λ_k = 2 − 2cos(kπ/N)` and the
  known single-edge-cut Δλ₂ — hard CI checks that validate the pipeline independent of data.
