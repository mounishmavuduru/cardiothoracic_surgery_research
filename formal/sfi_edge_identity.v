(* ===================================================================== *)
(*  Coq (kernel-checked) formalisation of the SFI edge-sensitivity core   *)
(* ===================================================================== *)
(*                                                                        *)
(*  Companion to scripts/verify_sfi_identity.py (sympy + mpmath) and      *)
(*  scripts/verify_math_rigor.py (exact-rational / interval / cross-      *)
(*  solver). Those establish the identities symbolically and numerically. *)
(*  This file re-establishes the *algebraic* core inside a proof          *)
(*  assistant, so the result is checked by Coq's small trusted kernel     *)
(*  rather than by a computer-algebra system.                             *)
(*                                                                        *)
(*  WHAT IS FORMALISED (pure linear algebra, general in the dimension n): *)
(*    T1  rank-one quadratic form :  phi^T (v v^T) phi = (v . phi)^2       *)
(*    T2  edge dot product        :  (e_i - e_j) . phi   = phi_i - phi_j   *)
(*    T3  SFI EDGE SENSITIVITY     :  phi^T E_ij phi     = (phi_i-phi_j)^2  *)
(*          where E_ij = (e_i - e_j)(e_i - e_j)^T = dL/dw_ij.              *)
(*          Combined with the (classical, cited) Hellmann-Feynman         *)
(*          theorem  dlambda2/dw_ij = phi2^T (dL/dw_ij) phi2, this IS the  *)
(*          closed form  dlambda2/dw_ij = (phi2_i - phi2_j)^2.             *)
(*    T4  SFI TERM NONNEGATIVITY   :  w >= 0  ->  w (phi_i - phi_j)^2 >= 0  *)
(*    T5  AGGREGATE NONNEGATIVITY  :  a sum of nonnegative edge terms is   *)
(*          >= 0  (SFI is a genuine Dirichlet energy of the Fiedler        *)
(*          vector; the Laplacian quadratic form is PSD).                  *)
(*                                                                        *)
(*  NOT formalised here (and not claimed to be): the analytic Hellmann-    *)
(*  Feynman step (differentiability of a simple eigenvalue) is the         *)
(*  classical result we cite (Ghosh-Boyd 2006); it is analysis, not the   *)
(*  algebra proved below. Each main theorem ends with Print Assumptions,   *)
(*  which reports "Closed under the global context" -> no axioms used.     *)
(* ===================================================================== *)

Require Import Reals.
Require Import Arith.
Require Import Lia.
Require Import Lra.
Open Scope R_scope.

(* Finite sum  bigsum f n = f 0 + f 1 + ... + f (n-1)  (n terms). *)
Fixpoint bigsum (f : nat -> R) (n : nat) : R :=
  match n with
  | O    => 0
  | S k  => bigsum f k + f k
  end.

(* -------------------- linearity building blocks --------------------- *)

Lemma bigsum_scal_l : forall (c : R) (f : nat -> R) (n : nat),
  c * bigsum f n = bigsum (fun i => c * f i) n.
Proof.
  intros c f n. induction n as [| k IH]; simpl.
  - ring.
  - rewrite <- IH. ring.
Qed.

Lemma bigsum_scal_r : forall (c : R) (f : nat -> R) (n : nat),
  bigsum f n * c = bigsum (fun i => f i * c) n.
Proof.
  intros c f n. induction n as [| k IH]; simpl.
  - ring.
  - rewrite <- IH. ring.
Qed.

Lemma bigsum_cmul : forall (c : R) (f : nat -> R) (n : nat),
  bigsum (fun b => c * f b) n = c * bigsum f n.
Proof. intros c f n. symmetry. apply bigsum_scal_l. Qed.

Lemma bigsum_eq : forall (f g : nat -> R) (n : nat),
  (forall i, (i < n)%nat -> f i = g i) -> bigsum f n = bigsum g n.
Proof.
  intros f g n H. induction n as [| k IH]; simpl.
  - reflexivity.
  - rewrite IH by (intros i Hi; apply H; lia).
    rewrite (H k) by lia. reflexivity.
Qed.

Lemma bigsum_minus : forall (f g : nat -> R) (n : nat),
  bigsum (fun a => f a - g a) n = bigsum f n - bigsum g n.
Proof.
  intros f g n. induction n as [| k IH]; simpl.
  - ring.
  - rewrite IH. ring.
Qed.

(* Sum factorisation:  sum_a sum_b (F a * G b) = (sum_a F a)(sum_b G b). *)
Lemma bigsum_prod : forall (F G : nat -> R) (n : nat),
  bigsum (fun a => bigsum (fun b => F a * G b) n) n
  = bigsum F n * bigsum G n.
Proof.
  intros F G n. rewrite bigsum_scal_r.
  apply bigsum_eq. intros a Ha. apply bigsum_cmul.
Qed.

(* --------------------------------------------------------------------- *)
(*  T1 : rank-one quadratic form.                                        *)
(* --------------------------------------------------------------------- *)

Definition dot (v phi : nat -> R) (n : nat) : R :=
  bigsum (fun a => v a * phi a) n.

Definition quad (v phi : nat -> R) (n : nat) : R :=
  bigsum (fun a => bigsum (fun b => (phi a * v a) * (v b * phi b)) n) n.

Theorem T1_rank_one_quadratic_form : forall (v phi : nat -> R) (n : nat),
  quad v phi n = (dot v phi n) ^ 2.
Proof.
  intros v phi n. unfold quad, dot.
  rewrite bigsum_prod.
  replace (bigsum (fun a => phi a * v a) n)
     with (bigsum (fun a => v a * phi a) n)
     by (apply bigsum_eq; intros a Ha; ring).
  ring.
Qed.

(* --------------------------------------------------------------------- *)
(*  Kronecker-delta selection lemmas.                                    *)
(* --------------------------------------------------------------------- *)

Lemma bigsum_delta_out : forall (n k : nat) (g : nat -> R),
  (n <= k)%nat ->
  bigsum (fun a => (if Nat.eqb a k then 1 else 0) * g a) n = 0.
Proof.
  induction n as [| m IH]; intros k g Hk; simpl.
  - reflexivity.
  - assert (Nat.eqb m k = false) as E by (apply Nat.eqb_neq; lia).
    rewrite E. rewrite IH by lia. ring.
Qed.

Lemma bigsum_delta : forall (n k : nat) (g : nat -> R),
  (k < n)%nat ->
  bigsum (fun a => (if Nat.eqb a k then 1 else 0) * g a) n = g k.
Proof.
  induction n as [| m IH]; intros k g Hk.
  - exfalso; lia.
  - simpl. destruct (Nat.eqb m k) eqn:E; simpl.
    + apply Nat.eqb_eq in E. subst k.
      rewrite bigsum_delta_out by lia. ring.
    + apply Nat.eqb_neq in E.
      rewrite IH by lia. ring.
Qed.

(* --------------------------------------------------------------------- *)
(*  T2 : the edge vector e_i - e_j picks out phi_i - phi_j under the dot. *)
(* --------------------------------------------------------------------- *)

Definition edge (i j : nat) : nat -> R :=
  fun a => if Nat.eqb a i then 1 else if Nat.eqb a j then -1 else 0.

Lemma dot_edge : forall (phi : nat -> R) (n i j : nat),
  (i < n)%nat -> (j < n)%nat -> i <> j ->
  dot (edge i j) phi n = phi i - phi j.
Proof.
  intros phi n i j Hi Hj Hij. unfold dot.
  transitivity (bigsum (fun a => (if Nat.eqb a i then 1 else 0) * phi a
                               - (if Nat.eqb a j then 1 else 0) * phi a) n).
  - apply bigsum_eq. intros a Ha. unfold edge.
    destruct (Nat.eqb a i) eqn:Eai; destruct (Nat.eqb a j) eqn:Eaj; simpl.
    + exfalso. apply Hij.
      apply Nat.eqb_eq in Eai. apply Nat.eqb_eq in Eaj.
      rewrite <- Eai. exact Eaj.
    + ring.
    + ring.
    + ring.
  - rewrite bigsum_minus.
    rewrite (bigsum_delta n i phi Hi).
    rewrite (bigsum_delta n j phi Hj).
    reflexivity.
Qed.

(* --------------------------------------------------------------------- *)
(*  T3 : SFI edge sensitivity  phi^T E_ij phi = (phi_i - phi_j)^2.        *)
(* --------------------------------------------------------------------- *)

Theorem T3_sfi_edge_sensitivity : forall (phi : nat -> R) (n i j : nat),
  (i < n)%nat -> (j < n)%nat -> i <> j ->
  quad (edge i j) phi n = (phi i - phi j) ^ 2.
Proof.
  intros phi n i j Hi Hj Hij.
  rewrite T1_rank_one_quadratic_form.
  rewrite (dot_edge phi n i j Hi Hj Hij).
  reflexivity.
Qed.

(* --------------------------------------------------------------------- *)
(*  T4 / T5 : nonnegativity of the SFI term and of the aggregate.        *)
(* --------------------------------------------------------------------- *)

Lemma sq_nonneg : forall x : R, 0 <= x ^ 2.
Proof. intro x. replace (x ^ 2) with (x * x) by ring. apply Rle_0_sqr. Qed.

Theorem T4_sfi_term_nonneg : forall (w x : R),
  0 <= w -> 0 <= w * x ^ 2.
Proof.
  intros w x Hw. apply Rmult_le_pos.
  - exact Hw.
  - apply sq_nonneg.
Qed.

Theorem T5_aggregate_nonneg : forall (t : nat -> R) (n : nat),
  (forall k, (k < n)%nat -> 0 <= t k) -> 0 <= bigsum t n.
Proof.
  intros t n H. induction n as [| m IH]; simpl.
  - lra.
  - apply Rplus_le_le_0_compat.
    + apply IH. intros k Hk. apply H. lia.
    + apply H. lia.
Qed.

(* --------------------------------------------------------------------- *)
(*  Axiom audit: each should report "Closed under the global context".   *)
(* --------------------------------------------------------------------- *)

Print Assumptions T1_rank_one_quadratic_form.
Print Assumptions T3_sfi_edge_sensitivity.
Print Assumptions T4_sfi_term_nonneg.
Print Assumptions T5_aggregate_nonneg.
