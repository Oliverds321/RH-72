#!/usr/bin/env python3
"""ext_checks.py -- external-review checks B1.1 and B1.5.
B1.1: is qg's v>=0 violation a discretisation artifact? Grid ladder on the unconstrained
      Nystrom solve at C=pi^4/18, b=b*: record min(v)/max(v) vs n, fit the decay rate.
B1.5: closing-condition margin sensitivity: require margin >= 0 / 1 nat / 0.1*L and
      recompute P_eff at 1e25/1e100/1e300."""
import numpy as np
from math import pi, log
import fb, budget_q as bd

print("=== B1.1: v>=0 ladder (C=pi^4/18, b*=0.6253660758, unconstrained solve) ===")
C, b = pi**4/18, 0.6253660758
prev = None
for n in (500, 1000, 2000, 4000, 8000):
    h = 2*b/n; t = -b + h*(np.arange(n)+0.5)
    D = np.abs(t[:,None]-t[None,:])
    K = np.where(D <= 1.0, D, C*D)
    A = np.eye(n)*h + K*h*h
    one = np.ones(n)*h
    x = np.linalg.solve(A, one); v = x/(one@x)
    r = v.min()/v.max()
    rate = "" if prev is None else f"  ratio vs prev: {r/prev:+.3f}"
    print(f"  n={n:5d}  min v/max v = {r:+.3e}{rate}")
    prev = r

print("=== B1.5: closing-margin sensitivity (design (3.5,3), 6w/L) ===")
orig = bd.closing_margin
lq = lambda x: x*log(10)
for label, mfun in (("margin >= 0 (printed table)", lambda L: 0.0),
                    ("margin >= 1 nat", lambda L: 1.0),
                    ("margin >= 0.1*L", lambda L: 0.1*L)):
    bd.closing_margin = lambda w, D0, L, LL, T, eta=1.0, _m=mfun: orig(w, D0, L, LL, T, eta) - _m(L)
    out = []
    for q in (25, 100, 300):
        d = bd.design(lq(q), r=3.5, K=3.0, c_ramp=6.0, curvature=True, minor=True)
        out.append(f"1e{q}: {d['P_eff']:+.4f}" if d else f"1e{q}: INFEASIBLE")
    print(f"  {label:28s} " + "   ".join(out))
bd.closing_margin = orig
