"""budget.py -- independent implementation of the design budget (round 7).

Rows (the author's own row set; see DRAFT 10.3, scripts/r4_design_check.py, round-5 p2_budget):
  zone     : curvature-honest  P(C,1) - P(C,a)   with  a = (1-d')(1 - l/LL)
             (linear model, the paper's: s*(1-a) with s = 0.5073)
  ramp     : c_ramp * w/L      c_ramp = 6 (the paper's conservative r2 = 6w/L; 4 = the
                               sharper retired reading, kept selectable for the comparison)
  buffer   : 6 D0/T            priced at the SHARP zero density, not the proved A0 (paper 10.3)
  ends     : 6 LL log LL / T
  minor    : L6..L11 of round 5
constraints: 8w <= L ([eq:wrange]); 10L <= D0 <= T/3; T >= 10 LL log LL;
             honest Gevrey closing (4/e) sqrt(wD0/A) >= L/2 + log(prefactor) + log(L/eta).
"""
import math
from math import log, pi, sqrt, e, exp
import numpy as np
from scipy.interpolate import CubicSpline
import fb

C_QQ = pi ** 4 / 18
A_G = 36 / e
B_G = 2 * e ** 8
c4 = 4 / e
S_LIN = 0.5073
REC = 0.672500703679412
DPDLOGC = (0.7212835668 - 0.6980745436) / log(2.0)

_spl = {}


def curve(C=C_QQ):
    """cubic splines a -> P(C,a), a -> lam*(C,a) on a = 1-delta in [0.60, 1]."""
    if C not in _spl:
        ag = np.concatenate([np.linspace(0.60, 0.94, 35), np.linspace(0.945, 1.0, 12)])
        Ps, Ls = [], []
        for a in ag:
            r = fb.payoff_a(C, a)
            Ps.append(r['P']); Ls.append(r['lam'])
        _spl[C] = (CubicSpline(ag, np.array(Ps)), CubicSpline(ag, np.array(Ls)),
                   fb.payoff_a(C)['P'])
    return _spl[C]


def closing_margin(w, D0, L, LL, T, eta=1.0):
    """(4/e)sqrt(wD0/A) - [L/2 + log prefactor + log(L/eta)]; the LEMMA_QT.b(5) prefactor."""
    C_env = e ** 2 * max(2 * B_G * w, L)
    tD = sqrt(w * D0 / A_G)
    S = 1 + (L / (2 * pi)) * (2 * A_G / w) * (tD / c4 + 1 / c4 ** 2)
    logpre = 2 * log(C_env) + log(2 * (LL + log(4 * T))) + 2 * log(S) - log(L)
    return c4 * tD - (L / 2 + logpre + log(L / eta))


def opt_ramp_buffer(L, LL, T, c_ramp=4.0, nw=500):
    """minimise c_ramp*w/L + 6*D0/T over the feasible set; returns (w, D0, cost) or None."""
    D0max = T / 3.0
    D0min = 10.0 * L
    if D0min >= D0max:
        return None
    def feas(w, D0):
        # w >= 1 is a regime hypothesis of the reused cap-free ends lemmas
        # (NOTE_QR QR.1: "8 <= L, 1 <= w, 4 <= c_rho, T-floors") -- enforced here.
        return (w >= 1.0) and (8 * w <= L * (1 + 1e-12)) and (D0min <= D0 <= D0max * (1 + 1e-12)) \
            and (closing_margin(w, D0, L, LL, T) >= 0.0)
    best = None
    for w in np.exp(np.linspace(log(1.0), log(L / 8), nw)):
        if not feas(w, D0max):
            continue
        lo, hi = D0min, D0max
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if feas(w, mid):
                hi = mid
            else:
                lo = mid
        cost = c_ramp * w / L + 6 * hi / T
        if best is None or cost < best[2]:
            best = (w, hi, cost)
    return best


def design(logQ, r=3.0, K=3.0, c_ramp=4.0, curvature=True, minor=True,
           lam_reopt=True, C=C_QQ, nw=500):
    Pspl, Lspl, P0 = curve(C)
    T = logQ ** r
    l = log(T) - log(2 * pi)
    LL = logQ + l
    if T < 10 * LL * log(LL):
        return None
    dp = K * log(LL) / LL
    a = (1 - dp) * (1 - l / LL)
    lam = float(Lspl(a)) if lam_reopt else float(Lspl(1.0))
    L = lam * LL
    zone = (P0 - float(Pspl(a))) if curvature else S_LIN * (1 - a)
    got = opt_ramp_buffer(L, LL, T, c_ramp=c_ramp, nw=nw)
    if got is None:
        return None
    w, D0, rb = got
    rows = dict(zone=zone, ramp=c_ramp * w / L, buffer=6 * D0 / T,
                ends=6 * LL * log(LL) / T)
    if minor:
        rows['minor'] = (S_LIN * 0.5 / LL
                         + 1.3 * 7.7 * sqrt(C_QQ * log(T * LL) / (T * LL))
                         + DPDLOGC * 2.763953 * sqrt(log(L) / (T * L))
                         + 0.35 * C_QQ * 2.763953 * sqrt(log(L) / (T * L))  # L12: IN-ZONE cross (Lemma 4.5): C*rho_U * in-zone sensitivity (0.35, conservative vs the review-computed ~0.32)
                         + DPDLOGC * 0.2 / T
                         + 3 * log(LL) ** 2 / L ** 2
                         + DPDLOGC * exp((lam - 2) * logQ))
    tot = sum(rows.values())
    return dict(logQ=logQ, r=r, K=K, T=T, LL=LL, l=l, dp=dp, a=a, delta=1 - a, lam=lam, L=L,
                w=w, D0=D0, wD0=w * D0, rows=rows, total=tot, P_eff=P0 - tot, P0=P0,
                margin=closing_margin(w, D0, L, LL, T), c_ramp=c_ramp,
                rate_const=tot * LL / log(LL))


def threshold(target, lo=20.0, hi=1e9, **kw):
    """smallest log Q with P_eff > target (bisection; P_eff is increasing in log Q)."""
    d = design(hi, **kw)
    if d is None or d['P_eff'] <= target:
        return None
    for _ in range(80):
        mid = sqrt(lo * hi)
        d = design(mid, **kw)
        if d is not None and d['P_eff'] > target:
            hi = mid
        else:
            lo = mid
    return hi
