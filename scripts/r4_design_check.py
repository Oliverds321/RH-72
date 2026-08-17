#!/usr/bin/env python3
"""r4_design_check.py — reconcile round-4 F1/F2 (buffer domination, effective range) and
validate the rebalanced design.

Budget model = the round-4 adversarial's (their Sec 8, read off Assembly.lean):
  L1 zone breakpoint = s * (log T - log 2pi)/LL      (s = |dP/ddelta| at the family's C)
  L2 zone edge       = s * K_dp * log(LL)/LL         (delta' = K_dp loglogQ/logQ; K_dp = 3
                                                      per R4; this session's the companion notes
                                                      K >= 2 -- take 3, conservative)
  L3 ramp            = 4 w / L
  L4 buffer          = 6 D0 / T
  L5 ends            = 4 LL log(LL) / T
subject to: two-factor Gevrey closing  w*D0 >= c_close * LL^2  (c_close = A e^2 lam^2/64,
            pair-repair eta/L target), and admissibility 8w <= L.

(A) reproduce the R4 table for the SHIPPED design (D0 = T/logLL, w = 169 logLL, T = LL^2).
(B) the rebalanced concrete design (T, D0, w) = (logQ^{3+eps}, LL^2 logLL, 3 logLL):
    admissibility, closing margin, budget, P_eff across Q; rate behaviour.
(C) best-achievable P_eff over (w, D0, eps) — reproduce R4's ~0.62 at Q=1e100 and the
    ~10^1064 threshold for 0.698.
"""
import numpy as np
from math import log, pi, sqrt, e

A_G = 36 / e
P_QQ, LAM_QQ, S_QQ = 0.7212835668, 1.2507321515, 0.5073   # q<=Q family (R4-remeasured slope)
K_DP = 3.0

def budget(logQ, T, D0, w, lam=LAM_QQ, s=S_QQ):
    LL = logQ + log(T) - log(2 * pi)
    L = lam * LL
    L1 = s * (log(T) - log(2 * pi)) / LL
    L2 = s * K_DP * log(LL) / LL
    L3 = 4 * w / L
    L4 = 6 * D0 / T
    L5 = 4 * LL * log(LL) / T
    admiss = 8 * w <= L
    c_close = A_G * e ** 2 * lam ** 2 / 64
    closes = w * D0 >= c_close * LL ** 2
    return LL, L1 + L2 + L3 + L4 + L5, (L1, L2, L3, L4, L5), admiss, closes

print("=== (A) shipped design (D0 = T/logLL, w = 169 logLL, T = LL^2) — R4 F1/F2 reproduction ===")
for lq in (230.3, 2302.6, 23026.0):        # Q = 1e100, 1e1000, 1e10000
    LL0 = lq                                # first pass; iterate LL = logQ + logT - log2pi
    for _ in range(5):
        T = LL0 ** 2
        LL0 = lq + log(T) - log(2 * pi)
    D0, w = T / log(LL0), 169 * log(LL0)
    LL, tot, parts, adm, cl = budget(lq, T, D0, w)
    print(f"  Q=1e{lq/2.302585:.0f}: 8w<=L? {'yes' if adm else 'NO':>3}  "
          f"L1..L5 = {parts[0]:.3f}/{parts[1]:.3f}/{parts[2]:.3f}/{parts[3]:.3f}/{parts[4]:.4f}"
          f"  total = {tot:.3f}  P_eff = {P_QQ - tot:+.3f}")
print("  (R4 quotes: at 1e100 ramp 12.97, buffer 1.095, total 14.22, P_eff -13.51; wrange NO)")

print("=== (B) rebalanced concrete design: T = (logQ)^(3+eps), D0 = LL^2 logLL, w = 3 logLL ===")
for eps in (0.0, 0.5):
    print(f"  --- eps = {eps} ---")
    for lq in (57.6, 230.3, 691, 2302.6, 23026.0):   # 1e25, 1e100, 1e300, 1e1000, 1e10000
        LL0 = lq
        for _ in range(6):
            T = lq ** (3 + eps)
            LL0 = lq + log(T) - log(2 * pi)
        LL = LL0
        D0, w = LL ** 2 * log(LL), 3 * log(LL)
        LLc, tot, parts, adm, cl = budget(lq, T, D0, w)
        print(f"    Q=1e{lq/2.302585:>6.0f}: adm {'y' if adm else 'N'} closes {'y' if cl else 'N'}"
              f"  L1..L5 = {parts[0]:.4f}/{parts[1]:.4f}/{parts[2]:.4f}/{parts[3]:.4f}/{parts[4]:.5f}"
              f"  total = {tot:.4f}  P_eff = {P_QQ - tot:+.4f}")
    # rate check: total should be O(logLL/LL)
print("  rate check at eps=0: total*LL/log(LL) across Q (should stabilise = O(loglogQ/logQ) class):")
vals = []
for lq in (230.3, 2302.6, 23026.0, 230260.0):
    for _ in range(6):
        T = lq ** 3
        LL = lq + log(T) - log(2 * pi)
    D0, w = LL ** 2 * log(LL), 3 * log(LL)
    _, tot, _, _, _ = budget(lq, T, D0, w)
    vals.append(tot * LL / log(LL))
print("   ", " / ".join(f"{v:.2f}" for v in vals))

print("=== (C) best-achievable P_eff over (w, D0, eps) — R4 reproduction ===")
def best_Peff(lq):
    best = -99
    for eps in np.arange(0.2, 3.2, 0.05):
        T = lq ** (2 + eps)
        LL = lq + log(T) - log(2 * pi)
        L = LAM_QQ * LL
        c_close = A_G * e ** 2 * LAM_QQ ** 2 / 64
        # optimal D0 given w*D0 = c_close LL^2: minimize 4w/L + 6D0/T
        # w = c LL^2/D0 -> f(D0) = 4 c LL^2/(L D0) + 6 D0/T -> D0* = sqrt(2 c LL^2 T/(3 L))
        D0 = sqrt(2 * c_close * LL ** 2 * T / (3 * L))
        w = c_close * LL ** 2 / D0
        if 8 * w > L: w = L / 8; D0 = c_close * LL ** 2 / w
        if D0 > T / 3: continue
        _, tot, _, adm, cl = budget(lq, T, D0, w)
        if adm and cl: best = max(best, P_QQ - tot)
    return best
for lq, tag in ((230.3, "1e100"), (2302.6, "1e1000"), (2450.0, "1e1064")):
    print(f"    Q={tag}: best P_eff = {best_Peff(lq):+.4f}")
print("  (R4 quotes: ~0.62 best at 1e100 [dyadic P; ours is q<=Q so slightly higher is fine];")
print("   0.698 reachable near 1e1064)")
