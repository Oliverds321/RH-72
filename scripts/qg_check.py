#!/usr/bin/env python3
"""qg_check.py — anchors for the companion notes g: assembly).

(1) Nystrom QP for the payoff functional B(v) = psi(0) + int_{|a|<=1}|a|psi + C int_{1<|a|<=lam}|a|psi,
    psi = v*v, int v = 1, supp v in [-b, b]:  B(v) = v^T M v with
    M = h*P_mirror + h^2*W, W_ij = w_C(t_i + t_j), and the RUN.md trap handled
    (kernel value (1+C)/2 on |t_i+t_j| = 1).  CALIBRATION GATES (rule 14):
      flat v on [-1/2,1/2] -> B = 4/3;  MT at lam = 1 -> 0.6725007 (Richardson);
      C = pi^4/18 at b = lam*/2 -> 0.7212836 (Richardson).
(2) THE RAMP STUDY (item g.5): multiply the optimum v* by a smooth edge ramp of relative
    width rho, renormalise, measure Delta P(rho).  Claim: because the free boundary has
    v*(b) = 0 (v* vanishes LINEARLY at the edge), the ramp cost is SECOND order:
    Delta P ~ rho^2, not the budgeted O(rho) = O(w/L).  Measure the exponent.
"""
import numpy as np

C_QQ = np.pi ** 4 / 18
LAM_STAR, B_STAR, P_STAR = 1.2507321515, 0.6253660758, 0.7212835668

def solve_B(b, C, n):
    """min v^T M v  s.t.  h*sum(v) = 1, supp [-b,b]; returns (B, v, t, h)."""
    t = np.linspace(-b, b, n)
    h = t[1] - t[0]
    lam = 2 * b
    S = t[:, None] + t[None, :]
    A = np.abs(S)
    W = np.where(A <= 1, A, np.where(A <= lam + 1e-12, C * A, 0.0))
    W[np.abs(A - 1) < h / 2] = (1 + C) / 2          # the RUN.md averaged-kernel trap
    M = h ** 2 * W
    M[np.arange(n), n - 1 - np.arange(n)] += h      # psi(0) = h * sum v_i v_{mirror(i)}
    M = 0.5 * (M + M.T)
    one = np.ones(n) * h
    x = np.linalg.solve(M, one)
    v = x / (one @ x)
    B = float(v @ M @ v)
    return B, v, t, h

def richardson(b, C, n1=2001, n2=4001):
    B1, *_ = solve_B(b, C, n1)
    B2, *_ = solve_B(b, C, n2)
    return 2 * B2 - B1

print("=== (1) calibration gates (rule 14) ===")
# flat-v gate: B(flat on [-1/2,1/2]) with C arbitrary (support edge exactly 1/2)
n = 2001
t = np.linspace(-0.5, 0.5, n); h = t[1] - t[0]
v = np.ones(n) / (n * h)
S = np.abs(t[:, None] + t[None, :])
W = np.where(S <= 1, S, 0.0)
M = h ** 2 * W
M[np.arange(n), n - 1 - np.arange(n)] += h
Bflat = float(v @ (0.5 * (M + M.T)) @ v)
print(f"    flat-v gate: B = {Bflat:.6f}  (target 4/3 = 1.333333)")
BMT = richardson(0.5, C_QQ)      # lam = 1: kernel never sees the C zone
print(f"    MT gate (lam = 1): P = {2 - BMT:.7f}  (target 0.6725007)")
BQQ = richardson(B_STAR, C_QQ)
print(f"    C = pi^4/18 at b* : P = {2 - BQQ:.7f}  (target {P_STAR:.7f})")
BDY = richardson(0.5965790605, 2 * np.pi ** 4 / 27)
print(f"    dyadic gate       : P = {2 - BDY:.7f}  (target 0.7099167)")

print("=== (2) ramp study at C = pi^4/18, b = b* ===")
n = 4001
B0, v0, t, h = solve_B(B_STAR, C_QQ, n)
print(f"    unramped P (n={n}) = {2 - B0:.7f};  min v*/max v* = {v0.min()/v0.max():.2e} "
      f"(v >= 0 {'OK' if v0.min() > -1e-12 else 'VIOLATED'})")
# edge behaviour: v* should vanish ~linearly at |t| = b (free boundary)
edge = np.abs(t) > B_STAR * 0.98
slope_fit = np.polyfit(B_STAR - np.abs(t[edge]), v0[edge], 1)
print(f"    edge fit v(t) ~ {slope_fit[0]:.3f}*(b-|t|) + {slope_fit[1]:.2e}  (linear vanish)")
print("    rho        P(rho)        DeltaP        DeltaP/rho^2")
rows = []
for rho in (0.02, 0.04, 0.08, 0.16):
    ramp = np.ones(n)
    m = np.abs(t) > B_STAR * (1 - rho)
    x = (B_STAR - np.abs(t[m])) / (B_STAR * rho)      # 0 at edge, 1 at inner
    ramp[m] = np.sin(0.5 * np.pi * np.clip(x, 0, 1)) ** 2
    vr = v0 * ramp
    vr = vr / (h * vr.sum())
    M = None
    # rebuild M for B evaluation
    S = np.abs(t[:, None] + t[None, :]); lam = 2 * B_STAR
    W = np.where(S <= 1, S, np.where(S <= lam + 1e-12, C_QQ * S, 0.0))
    W[np.abs(S - 1) < h / 2] = (1 + C_QQ) / 2
    M = h ** 2 * W
    M[np.arange(n), n - 1 - np.arange(n)] += h
    M = 0.5 * (M + M.T)
    Br = float(vr @ M @ vr)
    dP = (2 - B0) - (2 - Br)
    rows.append((rho, dP))
    print(f"    {rho:5.2f}   {2 - Br:.7f}   {dP:+.3e}    {dP / rho**2:+.4f}")
ex = np.polyfit(np.log([r for r, _ in rows]), np.log([abs(d) for _, d in rows]), 1)[0]
print(f"    measured scaling exponent of DeltaP in rho: {ex:.2f}  (2 = second order, the")
print(f"    free-boundary prediction; 1 would be the budgeted first-order O(w/L))")
