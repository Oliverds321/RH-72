"""Round-7 independent free-boundary solver for the C-penalised variational problem.

Same derivation as my round-3 implementation (kept as a separate file so this round's
numbers are produced by code in this round's repro/):

    B(v) = int v^2 + iint W_{C,a}(t-s) v(t)v(s),   W_{C,a}(y) = |y| (|y| <= a), C|y| (|y| > a)
    minimise over v >= 0, int v = 1, supp v subset [-b, b];  lam = 2b;  P = 2 - min B.

Rescale t -> t/a so the kernel break sits at 1:
    minimise (1/a) int w^2 + a iint W_C(t-s) w w  over int w = 1, w >= 0
EL twice-differentiated:
    w'' + a^2 [ 2w + (C-1)(w(t-1)+w(t+1)+w'(t-1)+w'(t+1)) ] = 0
  bulk |t| <= 1-b :  w = A cos(a sqrt2 t)
  edge 1-b <= |t| <= b, characteristic (rho^2 + 2a^2)^2 = a^4 (C-1)^2 (1 - rho^2)
  free boundary: C^1 match at t0 = 1-b, and w(b) = 0.   B = mu.
lam in the ORIGINAL variable is 2*a*b.

Gates (mandatory, run by every script that imports this):
  analytic MT value 0.672500703679412 ; flat-v 4/3 ; 0.7212835668 ; 0.7099167448.
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

SQ2 = math.sqrt(2.0)


def basis_a(C, a):
    K2 = (a * a * (C - 1.0)) ** 2
    p = 4 * a * a + K2
    qq = 4 * a ** 4 - K2
    s = math.sqrt(p * p - 4 * qq)
    out = []
    for z in ((-p + s) / 2.0, (-p - s) / 2.0):
        if z > 0:
            r = math.sqrt(z)
            beta = -(z + 2 * a * a) / (a * a * (C - 1.0) * (1.0 + r))
            s0 = math.log(-beta) / (2 * r)
            out.append((lambda s, r=r, s0=s0: math.sinh(r * (s - s0)),
                        lambda s, r=r, s0=s0: r * math.cosh(r * (s - s0))))
        else:
            wv = math.sqrt(-z)
            beta = -(-wv * wv + 2 * a * a) / (a * a * (C - 1.0) * (1.0 + 1j * wv))
            th = np.angle(beta)
            out.append((lambda s, wv=wv, th=th: math.cos(wv * s - th / 2),
                        lambda s, wv=wv, th=th: -wv * math.sin(wv * s - th / 2)))
    return out


def profile_a(C, a, b):
    (F1, F1p), (F2, F2p) = basis_a(C, a)
    t0 = 1.0 - b
    s0 = t0 - 0.5
    k = a * SQ2
    M = np.array([[F1(s0), F2(s0)], [F1p(s0), F2p(s0)]])
    rhs = np.array([math.cos(k * t0), -k * math.sin(k * t0)])
    c1, c2 = np.linalg.solve(M, rhs)
    def w(t):
        t = abs(t)
        return math.cos(k * t) if t <= t0 else c1 * F1(t - 0.5) + c2 * F2(t - 0.5)
    return w, w(b)


def payoff_a(C, a=1.0):
    """free-boundary payoff with orthogonality edge at |alpha| <= a."""
    g = lambda b: profile_a(C, a, b)[1]
    bs = np.linspace(0.5001, 1.9, 3000)
    vals = np.array([g(x) for x in bs])
    sg = np.sign(vals)
    idx = np.where(sg[:-1] * sg[1:] < 0)[0]
    if len(idx) == 0:
        raise RuntimeError("no free boundary root at C=%g a=%g" % (C, a))
    b = brentq(g, bs[idx[0]], bs[idx[0] + 1], xtol=1e-15, rtol=8.9e-16)
    w, res = profile_a(C, a, b)
    t0 = 1.0 - b
    WC = lambda y: abs(y) if abs(y) <= 1 else C * abs(y)
    Iv = 2 * (quad(w, 0, t0, limit=200)[0] + quad(w, t0, b, limit=200)[0])
    pts = [0, t0, b] if b <= 1 else [0, t0, 1, b]
    Iw = sum(quad(lambda s: WC(s) * w(s), pts[i], pts[i + 1], limit=200)[0]
             for i in range(len(pts) - 1))
    mu = ((1 / a) * w(0) + a * 2 * Iw) / Iv
    return dict(C=C, a=a, b=b, lam=2 * a * b, B=mu, P=2 - mu, res=res)


def payoff_grid(C, lam, n=4001):
    """constrained Nystrom solve at fixed lam (kernel jump averaged at |t-s| = 1)."""
    b = lam / 2.0
    t = np.linspace(-b, b, n)
    h = t[1] - t[0]
    D = np.abs(t[:, None] - t[None, :])
    W = np.where(np.abs(D - 1) < 1e-12, (1.0 + C) / 2.0 * D, np.where(D < 1, D, C * D))
    A = W + np.eye(n) / h
    one = np.ones(n)
    y = np.linalg.solve(A, one)
    x = y / (one @ y)
    if x.min() < 0:
        act = np.ones(n, bool)
        for _ in range(300):
            idx = np.where(act)[0]
            As = A[np.ix_(idx, idx)]
            os_ = np.ones(len(idx))
            ys = np.linalg.solve(As, os_)
            xs = ys / (os_ @ ys)
            if xs.min() >= -1e-14:
                break
            act[idx[xs < 0]] = False
        x = np.zeros(n)
        x[np.where(act)[0]] = xs
    B = x @ A @ x / (x.sum() ** 2)
    return B, 2 - B


def flat_v_payoff(C, lam):
    """B for the flat admissible v = 1/lam on [-lam/2, lam/2] (gate: C=1, lam=1 -> 2-2/3)."""
    b = lam / 2.0
    v = 1.0 / lam
    Iv2 = v * v * lam
    # iint W_C(t-s) v v = v^2 * iint W_C over the square
    def inner(t):
        pts = sorted(set([-b, b] + [x for x in (t - 1, t + 1) if -b < x < b]))
        return sum(quad(lambda u: (abs(t - u) if abs(t - u) <= 1 else C * abs(t - u)),
                        pts[i], pts[i + 1], limit=200)[0] for i in range(len(pts) - 1))
    Idd = v * v * quad(inner, -b, b, limit=200)[0]
    return Iv2 + Idd


def gates(verbose=True):
    """the four mandatory calibration gates; returns dict of residuals."""
    from mpmath import mp, mpf, sqrt as msqrt, sin as msin, cos as mcos
    mp.dps = 30
    r2 = msqrt(2)
    muMT = (msin(1 / r2) / r2 + mcos(1 / r2)) / (r2 * msin(1 / r2))
    P_MT = float(2 - muMT)
    B1, P1 = payoff_grid(1.0, 1.0, n=4001)
    B2, P2 = payoff_grid(1.0, 1.0, n=8001)
    rich = 2 * P2 - P1
    flat = flat_v_payoff(1.0, 1.0)
    p_qq = payoff_a(math.pi ** 4 / 18)
    p_dy = payoff_a(2 * math.pi ** 4 / 27)
    out = dict(MT_analytic=P_MT, MT_target=0.672500703679412,
               MT_nystrom_richardson=rich, flat_B=flat, flat_target=4.0 / 3.0,
               P_qq=p_qq['P'], P_qq_target=0.7212835668,
               lam_qq=p_qq['lam'], P_dyad=p_dy['P'], P_dyad_target=0.7099167448,
               lam_dyad=p_dy['lam'], res_qq=p_qq['res'])
    if verbose:
        print("=== MANDATORY GATES ===")
        print("  analytic MT              P = %.15f   (target 0.672500703679412, d=%.1e)"
              % (P_MT, P_MT - 0.672500703679412))
        print("  Nystrom C=1,lam=1 Rich.  P = %.9f    (vs analytic MT, d=%.1e)"
              % (rich, rich - P_MT))
        print("  flat-v C=1,lam=1         B = %.12f   (target 4/3, d=%.1e)"
              % (flat, flat - 4.0 / 3.0))
        print("  free boundary C=pi^4/18  P = %.10f  lam*=%.10f  (target 0.7212835668, d=%.1e)"
              % (p_qq['P'], p_qq['lam'], p_qq['P'] - 0.7212835668))
        print("  free boundary C=2pi^4/27 P = %.10f  lam*=%.10f  (target 0.7099167448, d=%.1e)"
              % (p_dy['P'], p_dy['lam'], p_dy['P'] - 0.7099167448))
        print("  |v(b*)| residual at C=pi^4/18: %.1e" % abs(p_qq['res']))
    return out


if __name__ == "__main__":
    gates()
