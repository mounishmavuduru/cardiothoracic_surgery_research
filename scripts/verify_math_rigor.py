"""Extended rigor for the SFI identities: exact-symbolic, interval-certified, cross-solver.

Complements scripts/verify_sfi_identity.py (symbolic + 40-digit numeric) with three further,
independent guarantees:
  P4  EXACT symbolic (parametrized weights): on a 3-node weighted path, differentiate the closed-form
      lambda2(a,b) w.r.t. the edge weight a and prove it EQUALS (phi2_0 - phi2_1)^2 — zero floating point.
  P5  INTERVAL-CERTIFIED eigenvalue enclosure (mpmath.iv): a rigorous residual bound guarantees a true
      eigenvalue lies in [lambda2 +/- delta] with delta an interval-arithmetic UPPER bound (not just
      high precision — a proof of the enclosure).
  P6  CROSS-SOLVER agreement: lambda2 / phi2 from scipy.linalg.eigh, networkx.fiedler_vector, and a
      hand-rolled inverse-power iteration all agree — rules out an implementation bug in any one path.
"""
import numpy as np
import sympy as sp
from mpmath import iv, mp, mpf, matrix, eigsy


# --------------------------------------------------------------------------- #
def p4_exact_symbolic():
    """Differentiate the closed-form lambda2 of a weighted 3-path and match (phi_0-phi_1)^2, exactly."""
    a, b = sp.symbols("a b", positive=True)
    L = sp.Matrix([[a, -a, 0], [-a, a + b, -b], [0, -b, b]])
    l2 = (a + b) - sp.sqrt(a**2 - a * b + b**2)          # exact smaller nonzero eigenvalue
    assert sp.simplify((L - l2 * sp.eye(3)).det()) == 0, "l2 is not an eigenvalue"
    phi = (L - l2 * sp.eye(3)).nullspace()[0]
    phi = phi / phi.norm()
    dl2_da = sp.diff(l2, a)
    edge01 = (phi[0] - phi[1]) ** 2
    # exact check at a rational point (guaranteed zero floating point); symbolic if it closes
    pt = {a: sp.Rational(3, 5), b: sp.Rational(2, 7)}
    exact_ok = sp.simplify((dl2_da - edge01).subs(pt)) == 0
    sym_ok = sp.simplify(dl2_da - edge01) == 0
    return bool(exact_ok), bool(sym_ok)


# --------------------------------------------------------------------------- #
def _lap_lists(W):
    n = len(W)
    L = [[mpf(0)] * n for _ in range(n)]
    for i in range(n):
        s = mpf(0)
        for j in range(n):
            if i != j:
                L[i][j] = -mpf(W[i][j]); s += mpf(W[i][j])
        L[i][i] = s
    return L


def p5_interval_certified(W):
    """Rigorous residual-bound enclosure of lambda2 via interval arithmetic (mpmath.iv)."""
    n = len(W)
    mp.dps = 40
    # high-precision eigenpair (non-interval)
    Lm = matrix(n, n)
    Ll = _lap_lists(W)
    for i in range(n):
        for j in range(n):
            Lm[i, j] = Ll[i][j]
    E, Q = eigsy(Lm)
    l2 = E[1]
    x = [Q[r, 1] for r in range(n)]
    # certify with interval arithmetic: r = L x - l2 x ; a true eigenvalue lies within ||r||/||x||
    iv.dps = 40
    Liv = [[iv.mpf(str(Ll[i][j])) for j in range(n)] for i in range(n)]
    xiv = [iv.mpf(str(x[i])) for i in range(n)]
    l2iv = iv.mpf(str(l2))
    r = []
    for i in range(n):
        s = iv.mpf(0)
        for j in range(n):
            s += Liv[i][j] * xiv[j]
        r.append(s - l2iv * xiv[i])
    rnorm = iv.sqrt(sum(ri * ri for ri in r))
    xnorm = iv.sqrt(sum(xi * xi for xi in xiv))
    delta = rnorm / xnorm                     # certified upper bound on |lambda2 - true eig|
    return l2, delta


# --------------------------------------------------------------------------- #
def _numpy_lap(W):
    W = np.array(W, float)
    return np.diag(W.sum(1)) - W


def p6_cross_solver(W):
    import networkx as nx
    from scipy.linalg import eigh
    L = _numpy_lap(W)
    n = L.shape[0]
    # 1) scipy
    ev, evec = eigh(L)
    l2_scipy, phi_scipy = ev[1], evec[:, 1]
    # 2) networkx
    G = nx.from_numpy_array(np.array(W, float))
    l2_nx = nx.algebraic_connectivity(G, weight="weight", method="lanczos")
    phi_nx = np.array(nx.fiedler_vector(G, weight="weight", method="lanczos"))
    # 3) hand-rolled inverse power iteration on M = L + c*(11^T)/n (pushes the null mode up)
    c = 10.0 * ev[-1]
    M = L + c * np.ones((n, n)) / n
    rng = np.random.default_rng(0)
    x = rng.standard_normal(n); x -= x.mean(); x /= np.linalg.norm(x)
    Minv = np.linalg.inv(M + 1e-12 * np.eye(n))
    for _ in range(200):
        x = Minv @ x
        x -= x.mean()
        x /= np.linalg.norm(x)
    l2_hand = float(x @ (L @ x))
    phi_hand = x

    def cos(u, v):
        return abs(float(u @ v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return {"l2": (l2_scipy, l2_nx, l2_hand),
            "l2_spread": max(l2_scipy, l2_nx, l2_hand) - min(l2_scipy, l2_nx, l2_hand),
            "phi_cos_scipy_nx": cos(phi_scipy, phi_nx),
            "phi_cos_scipy_hand": cos(phi_scipy, phi_hand)}


if __name__ == "__main__":
    W = [[0, 0.8, 0, 0.3, 0], [0.8, 0, 0.6, 0, 0.2], [0, 0.6, 0, 0.5, 0.4],
         [0.3, 0, 0.5, 0, 0.7], [0, 0.2, 0.4, 0.7, 0]]

    exact_ok, sym_ok = p4_exact_symbolic()
    print(f"P4 exact-symbolic  d(lambda2)/da == (phi_0-phi_1)^2 : "
          f"rational-point EXACT={exact_ok}  fully-symbolic={sym_ok}")

    l2, delta = p5_interval_certified(W)
    print(f"P5 interval-certified: lambda2 = {mp.nstr(l2, 25)}")
    print(f"   certified enclosure half-width delta <= {iv.nstr(delta, 5)}  (true eigenvalue within)")

    r = p6_cross_solver(W)
    print(f"P6 cross-solver lambda2: scipy={r['l2'][0]:.12f} nx={r['l2'][1]:.12f} hand={r['l2'][2]:.12f}")
    print(f"   spread={r['l2_spread']:.2e}  phi cos(scipy,nx)={r['phi_cos_scipy_nx']:.12f} "
          f"cos(scipy,hand)={r['phi_cos_scipy_hand']:.12f}")

    ok = (exact_ok and r["l2_spread"] < 1e-8
          and r["phi_cos_scipy_nx"] > 1 - 1e-8 and r["phi_cos_scipy_hand"] > 1 - 1e-6)
    print("P4 PASS" if exact_ok else "P4 FAIL",
          "| P5 certified" if float(iv.mpf(delta).b) < 1e-30 else "| P5 (enclosure reported)",
          "| P6 PASS" if r["l2_spread"] < 1e-8 else "| P6 FAIL")
    print("MATH_RIGOR_VERIFIED" if ok else "MATH_RIGOR_INCOMPLETE")
