# Formal (proof-assistant) verification of the SFI algebraic core

Kernel-checked proofs of the Spectral Fragility Index edge-sensitivity identity, complementing the
computer-algebra / numeric checks in `scripts/verify_sfi_identity.py` and
`scripts/verify_math_rigor.py`.

## Files
- **`sfi_edge_identity_Q.v`** — the primary result, over the exact rationals `Q`. **Axiom-free**:
  every `Print Assumptions` reports *"Closed under the global context"*. This is the formal analogue
  of the exact-rational `sympy.Rational` check (P4).
- **`sfi_edge_identity.v`** — the same theorems over the reals `R`. Mathematically identical, but
  inherits Coq's two standard real-number-library axioms (classical Dedekind-real completeness and
  functional extensionality). Kept to show the statements hold over ℝ as well.

## Theorems (both files)
| Name | Statement | Meaning |
|------|-----------|---------|
| `T1_rank_one_quadratic_form` | `φᵀ(vvᵀ)φ = (v·φ)²` | quadratic form of a rank-one edge matrix |
| `dot_edge` (T2) | `(e_i − e_j)·φ = φ_i − φ_j` | the edge vector selects the potential difference |
| `T3_sfi_edge_sensitivity` | `φᵀE_ijφ = (φ_i − φ_j)²` | **the SFI edge term** (`E_ij = ∂L/∂w_ij`) |
| `T4_sfi_term_nonneg` | `0≤w ⇒ 0 ≤ w(φ_i−φ_j)²` | each SFI contribution is non-negative |
| `T5_aggregate_nonneg` | `Σ nonneg terms ≥ 0` | SFI is a Dirichlet energy; the Laplacian form is PSD |

Combined with the **classical** Hellmann–Feynman theorem
`∂λ₂/∂w_ij = φ₂ᵀ(∂L/∂w_ij)φ₂` (Ghosh–Boyd 2006 — an analysis result we cite, not re-prove), T3 is
exactly the closed form `∂λ₂/∂w_ij = (φ₂,ᵢ − φ₂,ⱼ)²`.

## Reproduce
```
sudo apt-get install -y coq        # Coq 8.18 (any 8.x with QArith/Reals works)
coqc formal/sfi_edge_identity_Q.v  # prints "Closed under the global context" x4  -> zero axioms
coqc formal/sfi_edge_identity.v    # prints the two standard R-library axioms
```

## Why not Lean?
Lean/Mathlib was attempted. Lean's toolchain is distributed only as GitHub *release assets*, which
the compute environment's egress policy blocks with a hard 403 (release-asset host), and building
Lean 4 from source is infeasible here. Coq provides the same kernel-checked guarantee, so the formal
verification is done in Coq.
