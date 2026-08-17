#!/usr/bin/env python3
"""q7_check.py — numerical anchor for Lemma Q7 (rule 1/11: every object gets an anchor identity).

Checks, in the record's paperFT convention  h_f(z) = INT f(u) e^{izu} du :

  (1) the Parseval identity   Mform[u1,u2] := INT_{IxI} Phi(t-t')^2 u1(t) u2(t') dt dt'
                                           =  INT g(s) F1(s) conj(F2(s)) ds
      with Phi = paperFT(phi^2), g = phi^2 * phi^2 (convolution = autocorr for even phi^2),
      F_j(s) = INT_I u_j(t) e^{its} dt   -- NO 2*pi prefactor.
  (2) g >= 0 pointwise (the gv_nonneg surrogate).
  (3) the coefficient display: for nu = P_chi (Legendre mod 5, X = 200),
      F(s) = -(1/2pi) SUM_n Lambda(n) n^{-1/2} [chi(n) D_T(s - log n) + chi(n) D_T(s + log n)]
      with D_T(v) = (e^{2iTv} - e^{iTv})/(iv), matches direct quadrature of INT_I P e^{its} dt.
  (4) diagonal bookkeeping: INT g(s) |D_T(s - log n)|^2 ds  ~  2*pi*T*g(log n)  (rel. err O(1/T)),
      the constant behind the (T/pi) family-diagonal display.

Toy scales only (T = 6, L0 = 4): this validates CONVENTIONS AND STRUCTURE, not asymptotics.
"""
import numpy as np

# ---------- taper and derived objects ----------
L0 = 4.0                       # support of phi: [-L0/2, L0/2]
def phi2(x):                   # phi^2, even, C^1, compactly supported
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    m = np.abs(x) <= L0 / 2
    out[m] = np.cos(np.pi * x[m] / L0) ** 4
    return out

NX = 4001
xs = np.linspace(-L0 / 2, L0 / 2, NX)
wx = np.full(NX, xs[1] - xs[0]); wx[0] *= 0.5; wx[-1] *= 0.5   # trapezoid
p2 = phi2(xs)

def Phi(tau):                  # paperFT(phi^2)(tau) = INT phi^2(x) e^{i tau x} dx  (real, even)
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    return (p2 * wx) @ np.cos(np.outer(xs, tau))

def g_of(s):                   # (phi^2 * phi^2)(s), support [-L0, L0]
    s = np.atleast_1d(np.asarray(s, dtype=float))
    return np.array([np.sum(p2 * phi2(v - xs) * wx) for v in s])

# ---------- window ----------
T = 6.0
I0, I1 = T, 2 * T
NT = 1200
ts = np.linspace(I0, I1, NT)
wt = np.full(NT, ts[1] - ts[0]); wt[0] *= 0.5; wt[-1] *= 0.5

NS = 6001
ss = np.linspace(-L0, L0, NS)                 # supp g
ws = np.full(NS, ss[1] - ss[0]); ws[0] *= 0.5; ws[-1] *= 0.5
gs = g_of(ss)

_h = ts[1] - ts[0]
_uniq = Phi(np.arange(-(NT - 1), NT) * _h) ** 2          # Phi^2 on the 2NT-1 unique diffs
_K = np.empty((NT, NT))
for i in range(NT):
    _K[i, :] = _uniq[(NT - 1) + (i - np.arange(NT))]

def mform(nu1, nu2):
    """INT_{IxI} Phi(t - t')^2 nu1(t) nu2(t') dt dt'  by direct 2-D quadrature."""
    return float((nu1 * wt) @ _K @ (nu2 * wt))

def fhat(nu, s_arr):
    """F(s) = INT_I nu(t) e^{i t s} dt."""
    return (nu * wt) @ np.exp(1j * np.outer(ts, s_arr))

def parseval_rhs(nu1, nu2):
    F1 = fhat(nu1, ss); F2 = fhat(nu2, ss)
    return float(np.real(np.sum(gs * F1 * np.conj(F2) * ws)))

# ---------- (1) Parseval identity on random smooth densities ----------
rng = np.random.default_rng(23)
print("(1) Parseval identity  Mform = INT g |F|^2  (no 2pi):")
for trial in range(3):
    c = rng.standard_normal(5); f = rng.uniform(0.2, 2.0, 5)
    nu = sum(c[k] * np.cos(f[k] * ts) for k in range(5)) + 1.5
    lhs = mform(nu, nu); rhs = parseval_rhs(nu, nu)
    print(f"    trial {trial}: LHS={lhs:.12e}  RHS={rhs:.12e}  rel={abs(lhs-rhs)/abs(lhs):.2e}")

# ---------- (2) g >= 0 ----------
print(f"(2) min g on grid = {gs.min():.3e}  (>= 0 required)")

# ---------- (3) coefficient display for nu = P_chi, Legendre mod 5 ----------
def legendre5(n):
    r = n % 5
    return {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}[r]

X = 200
ns, lams, chis = [], [], []
for n in range(2, X + 1):
    m, p = n, None
    for q in range(2, n + 1):
        if n % q == 0:
            p = q; break
    k, m = 0, n
    while m % p == 0:
        m //= p; k += 1
    if m == 1:                      # n = p^k -> Lambda(n) = log p
        ns.append(n); lams.append(np.log(p)); chis.append(legendre5(n))
ns = np.array(ns); lams = np.array(lams); chis = np.array(chis, dtype=float)

P = -(1 / np.pi) * np.array([np.sum(lams * chis * ns ** -0.5 * np.cos(t * np.log(ns))) for t in ts])

def D_T(v):
    v = np.asarray(v, dtype=complex)
    out = np.where(np.abs(v) < 1e-12, T + 0j, (np.exp(2j * I0 * v) - np.exp(1j * I0 * v)) / (1j * v))
    return out

F_direct = fhat(P, ss)
F_closed = np.zeros(NS, dtype=complex)
for (n, lam, ch) in zip(ns, lams, chis):
    F_closed += -(1 / (2 * np.pi)) * lam * n ** -0.5 * ch * (D_T(ss - np.log(n)) + D_T(ss + np.log(n)))
err3 = np.max(np.abs(F_direct - F_closed)) / np.max(np.abs(F_direct))
print(f"(3) coefficient display (a_n(u), D_T weights): max rel err = {err3:.2e}")
lhsP = mform(P, P); rhsP = parseval_rhs(P, P)
print(f"    and Parseval on P_chi: LHS={lhsP:.10e}  RHS={rhsP:.10e}  rel={abs(lhsP-rhsP)/abs(lhsP):.2e}")

# ---------- (4) diagonal bookkeeping constant + the smearing term ----------
print("(4) INT g |D_T(s - log n)|^2 ds  vs  2*pi*T*g(log n)  [deviation = the D_T-smearing")
print("    term of Q7.iii, relative O((1+|g'/g|)/T); must shrink ~1/T]:")
def DT_at(v, Ta):
    v = np.asarray(v, dtype=complex)
    return np.where(np.abs(v) < 1e-12, Ta + 0j,
                    (np.exp(2j * Ta * v) - np.exp(1j * Ta * v)) / (1j * v))
for n in (2, 7):
    devs = []
    for Ta in (6.0, 24.0, 96.0):
        v = float(np.sum(gs * np.abs(DT_at(ss - np.log(n), Ta)) ** 2 * ws))
        tgt = 2 * np.pi * Ta * g_of(np.log(n))[0]
        devs.append(abs(v / tgt - 1))
    print(f"    n={n}: rel dev at T=6/24/96 = {devs[0]:.3e} / {devs[1]:.3e} / {devs[2]:.3e}"
          f"   (ratios {devs[0]/devs[1]:.1f}x, {devs[1]/devs[2]:.1f}x per 4x in T)")
