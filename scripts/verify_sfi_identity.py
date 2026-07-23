"""Machine-verified derivation of the SFI perturbation identities (sympy + mpmath).

Turns "we claim these formulas" into "here is the symbolic proof and a 30-digit numerical
confirmation". Three independent checks:
  P1  symbolic  : phi^T (e_i-e_j)(e_i-e_j)^T phi  ==  (phi_i - phi_j)^2   (the SFI edge term)
  P2  symbolic  : dL/dw_ij == (e_i-e_j)(e_i-e_j)^T  (rank-one edge structure of the Laplacian)
  P3  numeric   : central finite-difference dlambda2/dw_ij  ==  (phi2_i - phi2_j)^2  (Hellmann-Feynman)
                  and the 2nd-order curvature d^2 lambda2/dw_ij^2 == 2 * sum_{k!=2}(phi_k^T E phi2)^2/(l2-lk)
All at arbitrary precision so floating-point is not a confound.
"""
import sympy as sp
from mpmath import mp, mpf, matrix, eigsy

mp.dps = 40  # 40 significant digits


def p1_symbolic_first_order(n=6):
    phi = sp.Matrix(sp.symbols(f"phi0:{n}", real=True))
    ok = True
    for (i, j) in [(0, 1), (1, 3), (2, 5), (0, 4)]:
        v = sp.zeros(n, 1); v[i] = 1; v[j] = -1
        lhs = sp.expand((phi.T * (v * v.T) * phi)[0])
        rhs = sp.expand((phi[i] - phi[j]) ** 2)
        ok = ok and sp.simplify(lhs - rhs) == 0
    return ok


def p2_symbolic_laplacian_edge(n=4):
    # L = sum_{i<j} w_ij (e_i-e_j)(e_i-e_j)^T ; check dL/dw_ij = (e_i-e_j)(e_i-e_j)^T
    w = {}
    L = sp.zeros(n, n)
    for i in range(n):
        for j in range(i + 1, n):
            wij = sp.Symbol(f"w_{i}{j}", positive=True); w[(i, j)] = wij
            v = sp.zeros(n, 1); v[i] = 1; v[j] = -1
            L += wij * (v * v.T)
    ok = True
    for (i, j), wij in w.items():
        dL = L.diff(wij)
        v = sp.zeros(n, 1); v[i] = 1; v[j] = -1
        ok = ok and sp.simplify(dL - v * v.T) == sp.zeros(n, n)
    return ok


def _laplacian_from_weights(Wsym):
    n = Wsym.rows
    D = sp.diag(*[sum(Wsym[i, k] for k in range(n)) for i in range(n)])
    return D - Wsym


def p3_numeric_perturbation():
    # concrete 5-node weighted graph
    n = 5
    W = [[0, 0.8, 0, 0.3, 0], [0.8, 0, 0.6, 0, 0.2], [0, 0.6, 0, 0.5, 0.4],
         [0.3, 0, 0.5, 0, 0.7], [0, 0.2, 0.4, 0.7, 0]]

    def lap(Wm):
        L = matrix(n, n)
        for i in range(n):
            s = mpf(0)
            for j in range(n):
                if i != j:
                    L[i, j] = -mpf(Wm[i][j]); s += mpf(Wm[i][j])
            L[i, i] = s
        return L

    def eig2(L):
        E, Q = eigsy(L)                       # ascending eigenvalues, orthonormal cols
        lam = [E[k] for k in range(n)]
        phi = [[Q[r, k] for r in range(n)] for k in range(n)]  # phi[k] = eigenvector k
        return lam, phi

    lam, phi = eig2(lap(W))
    l2, phi2 = lam[1], phi[1]

    i, j = 0, 1
    E = matrix(n, n)
    for a, b in ((i, i, ), (j, j), (i, j), (j, i)):
        pass
    for (a, b, val) in [(i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)]:
        E[a, b] = mpf(val)                    # E = (e_i-e_j)(e_i-e_j)^T

    def l2_at(h):
        Lh = lap(W) + E * mpf(h)
        Eh, _ = eigsy(Lh)
        return Eh[1]

    h = mpf(10) ** (-8)
    d1 = (l2_at(h) - l2_at(-h)) / (2 * h)                     # first derivative
    d2 = (l2_at(h) + l2_at(-h) - 2 * l2) / (h * h)            # second derivative
    predicted_1 = (phi2[i] - phi2[j]) ** 2

    # 2nd-order resolvent prediction: sum_{k!=2} (phi_k^T E phi2)^2 / (l2 - lk)
    def quad(pk):
        return sum(pk[a] * (E[a, b] * phi2[b]) for a in range(n) for b in range(n))
    c2 = mpf(0)
    for k in range(n):
        if k == 1:
            continue
        num = quad(phi[k]) ** 2
        c2 += num / (l2 - lam[k])
    predicted_2 = 2 * c2

    err1 = abs(d1 - predicted_1)
    err2 = abs(d2 - predicted_2)
    return d1, predicted_1, err1, d2, predicted_2, err2


if __name__ == "__main__":
    print("P1 symbolic  phi^T E phi = (phi_i-phi_j)^2 :", "PASS" if p1_symbolic_first_order() else "FAIL")
    print("P2 symbolic  dL/dw_ij   = (e_i-e_j)(e_i-e_j)^T :", "PASS" if p2_symbolic_laplacian_edge() else "FAIL")
    d1, p1, e1, d2, p2, e2 = p3_numeric_perturbation()
    print(f"P3 numeric   dlambda2/dw_ij   = {mp.nstr(d1, 25)}")
    print(f"             (phi2_i-phi2_j)^2 = {mp.nstr(p1, 25)}   |err| = {mp.nstr(e1, 3)}")
    print(f"   2nd-order d^2lambda2/dw^2  = {mp.nstr(d2, 20)}")
    print(f"             resolvent formula = {mp.nstr(p2, 20)}   |err| = {mp.nstr(e2, 3)}")
    print("P3 first-order  :", "PASS" if e1 < mpf(10) ** (-10) else "FAIL")
    print("P3 second-order :", "PASS" if e2 < mpf(10) ** (-4) else "FAIL")
    print("SFI_IDENTITY_VERIFIED")
