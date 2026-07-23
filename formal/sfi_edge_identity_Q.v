(* ===================================================================== *)
(*  AXIOM-FREE Coq formalisation of the SFI edge-sensitivity core over Q  *)
(* ===================================================================== *)
(*                                                                        *)
(*  The companion file sfi_edge_identity.v proves the same identities     *)
(*  over the reals R, but Coq's real-number library is built on two       *)
(*  standard axioms (classical Dedekind-real completeness + functional    *)
(*  extensionality), so `Print Assumptions` there lists them. This file   *)
(*  re-proves the algebraic core over the *rationals* Q -- exact          *)
(*  arithmetic, constructed without axioms -- so `Print Assumptions`      *)
(*  reports "Closed under the global context": ZERO axioms. This is the   *)
(*  formal analogue of the exact-rational sympy.Rational check in         *)
(*  scripts/verify_math_rigor.py (P4).                                    *)
(*                                                                        *)
(*  Q carries a setoid equality (==), so identities are stated with ==.   *)
(*                                                                        *)
(*    T1  rank-one quadratic form :  phi^T (v v^T) phi == (v.phi)^2        *)
(*    T2  edge dot product        :  (e_i - e_j).phi    == phi_i - phi_j   *)
(*    T3  SFI EDGE SENSITIVITY     :  phi^T E_ij phi     == (phi_i-phi_j)^2 *)
(*    T4  SFI term nonnegativity   :  0<=w -> 0 <= w*(phi_i-phi_j)^2       *)
(*    T5  aggregate nonnegativity  :  sum of nonnegative edge terms >= 0   *)
(* ===================================================================== *)

Require Import QArith.
Require Import Arith.
Require Import Lia.
Open Scope Q_scope.

Fixpoint bigsum (f : nat -> Q) (n : nat) : Q :=
  match n with
  | O    => 0
  | S k  => bigsum f k + f k
  end.

(* -------------------- linearity building blocks --------------------- *)

Lemma bigsum_scal_l : forall (c : Q) (f : nat -> Q) (n : nat),
  c * bigsum f n == bigsum (fun i => c * f i) n.
Proof.
  intros c f n. induction n as [| k IH]; simpl.
  - ring.
  - rewrite <- IH. ring.
Qed.

Lemma bigsum_scal_r : forall (c : Q) (f : nat -> Q) (n : nat),
  bigsum f n * c == bigsum (fun i => f i * c) n.
Proof.
  intros c f n. induction n as [| k IH]; simpl.
  - ring.
  - rewrite <- IH. ring.
Qed.

Lemma bigsum_cmul : forall (c : Q) (f : nat -> Q) (n : nat),
  bigsum (fun b => c * f b) n == c * bigsum f n.
Proof. intros c f n. symmetry. apply bigsum_scal_l. Qed.

Lemma bigsum_eq : forall (f g : nat -> Q) (n : nat),
  (forall i, (i < n)%nat -> f i == g i) -> bigsum f n == bigsum g n.
Proof.
  intros f g n H. induction n as [| k IH]; simpl.
  - reflexivity.
  - rewrite IH by (intros i Hi; apply H; lia).
    rewrite (H k) by lia. reflexivity.
Qed.

Lemma bigsum_minus : forall (f g : nat -> Q) (n : nat),
  bigsum (fun a => f a - g a) n == bigsum f n - bigsum g n.
Proof.
  intros f g n. induction n as [| k IH]; simpl.
  - ring.
  - rewrite IH. ring.
Qed.

Lemma bigsum_prod : forall (F G : nat -> Q) (n : nat),
  bigsum (fun a => bigsum (fun b => F a * G b) n) n == bigsum F n * bigsum G n.
Proof.
  intros F G n. rewrite bigsum_scal_r.
  apply bigsum_eq. intros a Ha. apply bigsum_cmul.
Qed.

(* -------------------------- T1 : rank-one form ---------------------- *)

Definition dot (v phi : nat -> Q) (n : nat) : Q :=
  bigsum (fun a => v a * phi a) n.

Definition quad (v phi : nat -> Q) (n : nat) : Q :=
  bigsum (fun a => bigsum (fun b => (phi a * v a) * (v b * phi b)) n) n.

Theorem T1_rank_one_quadratic_form : forall (v phi : nat -> Q) (n : nat),
  quad v phi n == dot v phi n * dot v phi n.
Proof.
  intros v phi n. unfold quad, dot.
  rewrite bigsum_prod.
  assert (H : bigsum (fun a => phi a * v a) n == bigsum (fun a => v a * phi a) n)
    by (apply bigsum_eq; intros a Ha; ring).
  rewrite H. reflexivity.
Qed.

(* ------------------- Kronecker-delta selection ---------------------- *)

Lemma bigsum_delta_out : forall (n k : nat) (g : nat -> Q),
  (n <= k)%nat ->
  bigsum (fun a => (if Nat.eqb a k then 1 else 0) * g a) n == 0.
Proof.
  induction n as [| m IH]; intros k g Hk; simpl.
  - reflexivity.
  - assert (E : Nat.eqb m k = false) by (apply Nat.eqb_neq; lia).
    rewrite E.
    replace (if false then 1 else 0) with (0:Q) by reflexivity.
    rewrite IH by lia. ring.
Qed.

Lemma bigsum_delta : forall (n k : nat) (g : nat -> Q),
  (k < n)%nat ->
  bigsum (fun a => (if Nat.eqb a k then 1 else 0) * g a) n == g k.
Proof.
  induction n as [| m IH]; intros k g Hk.
  - exfalso; lia.
  - simpl. destruct (Nat.eqb m k) eqn:E.
    + apply Nat.eqb_eq in E. subst k.
      replace (if true then 1 else 0) with (1:Q) by reflexivity.
      rewrite bigsum_delta_out by lia. ring.
    + apply Nat.eqb_neq in E.
      replace (if false then 1 else 0) with (0:Q) by reflexivity.
      rewrite IH by lia. ring.
Qed.

(* ---------------------------- T2 : edge dot ------------------------- *)

Definition edge (i j : nat) : nat -> Q :=
  fun a => if Nat.eqb a i then 1 else if Nat.eqb a j then -1 else 0.

Lemma dot_edge : forall (phi : nat -> Q) (n i j : nat),
  (i < n)%nat -> (j < n)%nat -> i <> j ->
  dot (edge i j) phi n == phi i - phi j.
Proof.
  intros phi n i j Hi Hj Hij. unfold dot.
  transitivity (bigsum (fun a => (if Nat.eqb a i then 1 else 0) * phi a
                               - (if Nat.eqb a j then 1 else 0) * phi a) n).
  - apply bigsum_eq. intros a Ha. unfold edge.
    destruct (Nat.eqb a i) eqn:Eai; destruct (Nat.eqb a j) eqn:Eaj.
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

(* ------------------------ T3 : SFI edge sensitivity ----------------- *)

Theorem T3_sfi_edge_sensitivity : forall (phi : nat -> Q) (n i j : nat),
  (i < n)%nat -> (j < n)%nat -> i <> j ->
  quad (edge i j) phi n == (phi i - phi j) * (phi i - phi j).
Proof.
  intros phi n i j Hi Hj Hij.
  rewrite T1_rank_one_quadratic_form.
  rewrite (dot_edge phi n i j Hi Hj Hij).
  reflexivity.
Qed.

(* --------------------- T4 / T5 : nonnegativity ---------------------- *)

Lemma Qopp_nonneg : forall x : Q, x <= 0 -> 0 <= - x.
Proof.
  intros x H. setoid_replace 0 with (- 0) by ring.
  apply Qopp_le_compat. exact H.
Qed.

Lemma Qsqr_nonneg : forall x : Q, 0 <= x * x.
Proof.
  intro x. destruct (Qlt_le_dec x 0) as [Hlt | Hge].
  - setoid_replace (x * x) with (- x * - x) by ring.
    apply Qmult_le_0_compat; apply Qopp_nonneg; apply Qlt_le_weak; exact Hlt.
  - apply Qmult_le_0_compat; exact Hge.
Qed.

Theorem T4_sfi_term_nonneg : forall (w x : Q),
  0 <= w -> 0 <= w * (x * x).
Proof.
  intros w x Hw. apply Qmult_le_0_compat.
  - exact Hw.
  - apply Qsqr_nonneg.
Qed.

Lemma Qplus_nonneg : forall a b : Q, 0 <= a -> 0 <= b -> 0 <= a + b.
Proof.
  intros a b Ha Hb. setoid_replace 0 with (0 + 0) by ring.
  apply Qplus_le_compat; assumption.
Qed.

Theorem T5_aggregate_nonneg : forall (t : nat -> Q) (n : nat),
  (forall k, (k < n)%nat -> 0 <= t k) -> 0 <= bigsum t n.
Proof.
  intros t n H. induction n as [| m IH]; simpl.
  - apply Qle_refl.
  - apply Qplus_nonneg.
    + apply IH. intros k Hk. apply H. lia.
    + apply H. lia.
Qed.

(* --------------- axiom audit: expect "Closed under the global context" *)
Print Assumptions T1_rank_one_quadratic_form.
Print Assumptions T3_sfi_edge_sensitivity.
Print Assumptions T4_sfi_term_nonneg.
Print Assumptions T5_aggregate_nonneg.
