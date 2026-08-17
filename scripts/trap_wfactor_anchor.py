#!/usr/bin/env python3
"""trap_wfactor_anchor.py -- anchors for two quoted-but-previously-unanchored figures.

(1) THE 0.3693 TRAP (Corollary 3 / §12.3): the value of the variational problem with the
    in-zone coefficient doubled (in-zone weight 2|α|, out-zone 2C at C = π⁴/18, i.e. the
    blanket-positivity misreading the corollary warns against). Exact-cell Nyström on the
    kernel K(y) = 2·min-region + 2C·outer, support [-b,b] optimised, v ≥ 0 checked a
    posteriori.
(2) THE DECLINED q/φ(q) WEIGHT (§6 remark): ⟨q/φ(q)⟩_{φ*} at Q = 2×10⁶, the implied
    C_w, and P(C_w) via the gated free-boundary solver (fb.py) -> the +1.29/+1.01 figures.
"""
import numpy as np
from math import pi, log
import fb

print("=== (1) the 0.3693 trap: in-zone 2|a|, out-zone 2C (C = 2pi^4/27, the dyadic constant of Corollary 3) ===")
C = 2*pi**4/27
def solve(b, n):
    h = 2.0*b/n; t = -b + h*(np.arange(n)+0.5)
    D = np.abs(t[:,None]-t[None,:])
    K = np.where(D <= 1.0, 2.0*D, 2.0*C*D)  # 2|y| in-zone; 2C|y| out-zone
    A = np.eye(n)*h + K*h*h
    one = np.ones(n)*h
    x = np.linalg.solve(A, one); v = x/(one@x)
    return v@A@v, v
best = None
for b in np.arange(0.30, 1.01, 0.01):
    B, v = solve(b, 900)
    if v.min() >= -1e-9 and (best is None or B < best[0]):
        best = (B, b)
B_c, b_c = best
# refine at the best support
Bs = {}
for n in (800, 1600, 3200):
    B, v = solve(b_c, n); Bs[n] = B
Bext = Bs[3200] + (Bs[3200]-Bs[1600])/3
print(f"  best support b = {b_c:.2f} (lam = {2*b_c:.2f});  P = 2 - B:")
for n in (800,1600,3200): print(f"    n={n}: P = {2-Bs[n]:.6f}")
print(f"  Richardson: P = {2-Bext:.6f}   (quoted trap value 0.3693; review value 0.369289)")
print(f"  min v/max v at n=3200: {v.min()/v.max():+.2e}  (v >= 0 non-binding check)")

print("=== (2) the q/phi(q) weight, declined in §6 ===")
Q = 2*10**6
phi = np.arange(Q+1); is_p = np.ones(Q+1, bool); is_p[:2] = False
for p in range(2, int(Q**0.5)+1):
    if is_p[p]: is_p[p*p::p] = False
phi = np.arange(Q+1, dtype=np.float64)
for p in range(2, Q+1):
    if is_p[p]: phi[p::p] *= (1 - 1/p)
# phi*(q) = number of primitive characters mod q = sum_{d|q} mu(d) phi(q/d); use phi* via known identity:
# phi*(q) = q * prod_{p||q}(1-2/p) * prod_{p^2|q}(1-1/p)^2  -- compute via factorization sieve
smallest = np.zeros(Q+1, int)
for p in range(2, Q+1):
    if is_p[p]: smallest[p::p] = np.where(smallest[p::p]==0, p, smallest[p::p])
phistar = np.zeros(Q+1)
phistar[1] = 1
for q in range(3, Q+1):
    x, val = q, 1.0
    while x > 1:
        p = smallest[x]; a = 0
        while x % p == 0: x //= p; a += 1
        val *= (p-2) if a == 1 else (p**a)*(1-1/p)**2
    phistar[q] = val
w = phistar[3:] * (np.arange(3, Q+1)/phi[3:])
num, den = w.sum(), phistar[3:].sum()
ratio = num/den
Cw = (pi**4/18)/ratio
print(f"  <q/phi(q)>_phi* at Q=2e6 = {ratio:.6f}   (paper: 1.2965)")
print(f"  C_w = (pi^4/18)/ratio = {Cw:.4f}          (paper: ~4.174)")
r = fb.payoff_a(Cw)
r2 = fb.payoff_a((2*pi**4/27)/ratio)
print(f"  P(C_w)      = {r['P']:.6f}  -> +{100*(r['P']-0.7212835668):.2f} pp vs q<=Q   (paper: ~+1.29)")
print(f"  P(C_w,dyad) = {r2['P']:.6f} -> +{100*(r2['P']-0.7099167448):.2f} pp vs dyadic (paper: ~+1.01)")
