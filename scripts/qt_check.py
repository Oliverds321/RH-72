#!/usr/bin/env python3
"""qt_check.py — numerical anchors for the companion notes transform bound, mirrored tail
lemma, F2 pair repair).  Rule 1/11: every constant transcribed from the tree gets a check.

(A) ramp-derivative constants: for the ACTUAL rhoTwo (theta = g(x)g(1-x), g = exp(-1/x)),
    verify ||theta^(k-1)||_1 / int(theta) <= B A^k k^{2k}  (the content of the record's
    proved ramp lemma ||phi^(k)||_1 <= 2Bw(A/w)^k k^{2k}), A = 36/e, B = 2e^8, k = 1..8.
(B) QT.a envelope: |phihat(r)| <= C_env e^{|Im|L/2} exp(-(2/e)sqrt(w r/A)),
    C_env = e^2 max(2Bw, L), on the real line (Im = 0), actual rhoTwo taper, L=30, w=2.
(C) row-sum closed bound: sum_{k<d} exp(-c sqrt(w(D+kh)/A))
      <= exp(-c t_D)(1 + (L/2pi)(2A/w)(t_D/c + 1/c^2)),  t_D = sqrt(wD/A), c = 4/e.
(D) closing-condition ladder: supplied exponent (4/e or 2/e)sqrt(w D0/A) vs needed
    L/2 (+ 2logQ) + logP at D0 = T/log(cal L), w = c2 log(cal L), for the c2 ladder,
    at the critical T = L^2 and the adopted T = L^{2+eps}.
(E) F2 pair repair toy: on random Hermitian blocks, verify the pair perturbation
    inequality  4trG - frobSq G - [4B_tr + 2 B_F sqrt(frobSq G) + B_F^2] <= 4trA - frobSq A
    for the direct sum, and compare the pair loss vs the scalar-B loss (round-3 F2 reading).
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 40
E = mp.e
A_G = mp.mpf(36) / E          # record's Gevrey A
B_G = 2 * mp.e**8             # record's Gevrey B

print("=== (A) ramp-derivative constants of rhoTwo vs the record's (A,B) = (36/e, 2e^8) ===")
import sympy as sp
x = sp.symbols('x')
theta_sym = sp.exp(-1/x) * sp.exp(-1/(1-x))
int_theta = mp.quad(lambda v: mp.e**(-1/v) * mp.e**(-1/(1-v)), [0, 0.5, 1])
print(f"    int_0^1 theta = {mp.nstr(int_theta, 8)}   (1/int = {mp.nstr(1/int_theta, 6)}; "
      f"B = 2e^8 = {mp.nstr(B_G, 6)} conservative by {mp.nstr(B_G*int_theta/2, 4)}x... "
      f"[B vs 1/int: ratio {mp.nstr(B_G/(1/int_theta), 4)}])")
d = theta_sym
for k in range(1, 9):
    # ||rho2^(k)||_1 = ||theta^(k-1)||_1 / int(theta)
    f = sp.lambdify(x, sp.Abs(d), 'mpmath')
    l1 = mp.quad(lambda v: f(v), mp.linspace(0.001, 0.999, 41)) / int_theta
    bound = B_G * A_G**k * mp.mpf(k)**(2*k)
    ok = "OK " if l1 <= bound else "FAIL"
    print(f"    k={k}: ||rho2^(k)||_1 = {mp.nstr(l1, 4):>12}   bound B A^k k^2k = {mp.nstr(bound, 4):>12}   "
          f"ratio {mp.nstr(l1/bound, 3):>10}  {ok}")
    d = sp.diff(d, x)

print("=== (B) QT.a envelope on the real line: actual rhoTwo taper, L = 30, w = 2 ===")
L_t, w_t = mp.mpf(30), mp.mpf(2)
C_env = E**2 * max(2 * B_G * w_t, L_t)
theta_mp = lambda v: mp.e**(-1/v) * mp.e**(-1/(1-v)) if 0 < v < 1 else mp.mpf(0)
def phihat_abs(r):
    """|phihat(r)| = (2/(r int)) |int_0^1 theta(v) sin(r(L/2 - w v)) dv|  (one part, exact)."""
    r = mp.mpf(r)
    nseg = max(8, int(3 * r * w_t / mp.pi))
    val = mp.quad(lambda v: theta_mp(v) * mp.sin(r * (L_t/2 - w_t * v)), mp.linspace(0, 1, nseg + 1))
    return abs(2 * val / (r * int_theta))
print("    r        |phihat(r)|      envelope         log-slack (nats)")
for r in (5, 10, 20, 40, 80, 160, 320):
    ph = phihat_abs(r)
    env = C_env * mp.e**(-(2/E) * mp.sqrt(w_t * r / A_G))
    slack = mp.log(env / ph) if ph > 0 else mp.inf
    ok = "OK " if ph <= env else "FAIL"
    print(f"    {r:>4}   {mp.nstr(ph, 6):>14}   {mp.nstr(env, 6):>14}   {mp.nstr(slack, 4):>8}  {ok}")

print("=== (C) row-sum closed bound, c = 4/e ===")
c4 = 4 / E
def row_lhs(D, L, w, d):
    h = 2 * mp.pi / L
    return sum(mp.e**(-c4 * mp.sqrt(w * (D + k * h) / A_G)) for k in range(d))
def row_rhs(D, L, w):
    tD = mp.sqrt(w * D / A_G)
    return mp.e**(-c4 * tD) * (1 + (L / (2 * mp.pi)) * (2 * A_G / w) * (tD / c4 + 1 / c4**2))
for (D, L, w, d) in [(10, 30, 2, 200), (50, 30, 2, 200), (100, 120, 8, 500)]:
    lhs, rhs = row_lhs(D, L, w, d), row_rhs(D, L, w)
    ok = "OK " if lhs <= rhs else "FAIL"
    print(f"    D={D:>4} L={L:>4} w={w}: lhs = {mp.nstr(lhs, 5):>12}  rhs = {mp.nstr(rhs, 5):>12}  "
          f"ratio {mp.nstr(lhs/rhs, 3)}  {ok}")

print("=== (D) closing-condition ladder:  supplied vs needed exponent (nats) ===")
print("    supplied(2f) = (4/e)sqrt(c2 T/A);  needed(eta/L) = L/2 + logP + log(L/eta),")
print("    needed(Q^-2) = L/2 + 2logQ + logP;  single-factor supplies half.  eta = 0.01.")
lam = mp.mpf('1.1931581210')   # dyadic lambda*
for eps_exp in (mp.mpf(2), mp.mpf('2.2')):
    print(f"    --- T = (logQ)^{mp.nstr(eps_exp,3)},  D0 = T/log(calL), w = c2 log(calL) ---")
    for logQ in (50, 230.3, 1000, 10000):
        calL = mp.mpf(logQ)                     # calL = logQ (1+o(1))
        T = calL**eps_exp
        Lb = lam * calL
        loglog = mp.log(calL)
        for c2 in (3, 45, 170):
            wD0 = c2 * loglog * (T / loglog)    # = c2 T
            sup2 = (4 / E) * mp.sqrt(wD0 / A_G)
            # polynomial prefactor P = A0 * C_env^2 * 2calL * S(D0)^2 / L, A0 ~ 1;
            w_ = c2 * loglog
            Cenv = E**2 * max(2 * B_G * w_, Lb)
            tD0 = mp.sqrt(wD0 / A_G)
            S = 1 + (Lb / (2 * mp.pi)) * (2 * A_G / w_) * (tD0 / c4 + 1 / c4**2)
            logP = mp.log(2 * calL * Cenv**2 * S**2 / Lb)
            need_eta = Lb / 2 + logP + mp.log(Lb / mp.mpf('0.01'))
            need_Q2 = Lb / 2 + 2 * calL + logP
            verdict = ("2f:eta" + ("+" if sup2 >= need_eta else "-")
                       + " 2f:Q2" + ("+" if sup2 >= need_Q2 else "-")
                       + " 1f:eta" + ("+" if sup2 / 2 >= need_eta else "-")
                       + " 1f:Q2" + ("+" if sup2 / 2 >= need_Q2 else "-"))
            print(f"      logQ={logQ:>7} c2={c2:>4}: supplied={mp.nstr(sup2, 6):>10} "
                  f"needed eta/L={mp.nstr(need_eta, 6):>10} Q^-2={mp.nstr(need_Q2, 6):>10}   {verdict}")

print("=== (E) F2 pair repair: direct-sum toy, 60 blocks of dim 8 ===")
rng = np.random.default_rng(17)
M, n = 60, 8
lossG_true = 0.0; trG = 0.0; frG = 0.0
b_tr_sum = 0.0; b_sq_sum = 0.0; b_max = 0.0
lhs_pair_terms = []
for i in range(M):
    Gh = rng.standard_normal((n, n)); Gh = (Gh + Gh.T) / 2 + 3 * np.eye(n)
    Eh = rng.standard_normal((n, n)) * 0.05; Eh = (Eh + Eh.T) / 2
    Ah = Gh - Eh
    trG += np.trace(Gh); frG += np.sum(Gh * Gh)
    lossG_true += (4 * np.trace(Ah) - np.sum(Ah * Ah)) - (4 * np.trace(Gh) - np.sum(Gh * Gh))
    tn = np.abs(np.linalg.eigvalsh(Eh)).sum()          # trace norm per block
    b_tr_sum += tn; b_sq_sum += tn**2; b_max = max(b_max, tn)
B_tr = b_tr_sum                       # >= |tr E_fam|
B_F = np.sqrt(b_sq_sum)               # >= ||E_fam||_F  (since ||E||_F <= ||E||_1 per block)
pair_loss = 4 * B_tr + 2 * B_F * np.sqrt(frG) + B_F**2
scalar_B = b_tr_sum                   # round-3 reading: single scalar B = |F| theta0/(aL)-type
scalar_loss = scalar_B * (4 + 2 * np.sqrt(frG) + scalar_B)
print(f"    true |loss| = {abs(lossG_true):.3f}")
print(f"    pair-lemma guaranteed loss  = {pair_loss:.3f}   (holds: {abs(lossG_true) <= pair_loss})")
print(f"    scalar-B guaranteed loss    = {scalar_loss:.3f}   (inflation vs pair: "
      f"{scalar_loss / pair_loss:.1f}x; grows like sqrt(#blocks) in the dominant term)")
