#!/usr/bin/env python3
"""qf_check.py — anchors for NOTE_QF (item f: q-uniform EF constants + ParamsQ).

(1) LGrowth anchor: ||L(s,chi)|| <= q*|s|/Re(s) on Re s > 0 (the tree's explicit q-uniform
    growth bound, LGrowth.lean:210) — verified numerically for actual primitive characters
    (L via Hurwitz zeta), a grid of s.
(2) A0 anchor on the real dataset: for all 60 primitive characters (moduli 5..19) in
    family_zeros.pkl, measure max over unit windows of N_chi(t, t+1] / log(q(|t|+3)) —
    the target local-count shape N_chi(t,t+1] <= A0 log(q(|t|+3)) with ABSOLUTE A0.
    (Diagnostics only, never evidence — round-1 rule.)
"""
import sys, os, pickle, numpy as np, mpmath as mp
mp.mp.dps = 25
_default = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "family_zeros.pkl")
if not os.path.exists(_default):
    _default = "data/family_zeros.pkl"
DATASET = sys.argv[1] if len(sys.argv) > 1 else _default

def prim_chars(q):
    def order(x):
        k, y = 1, x % q
        while y != 1: y = y * x % q; k += 1
        return k
    g = next(u for u in range(2, q) if order(u) == q - 1)
    dlog = {}; x = 1
    for k in range(q - 1): dlog[x] = k; x = x * g % q
    return [(j, np.array([0] + [np.exp(2j * np.pi * jj * dlog[a] / (q - 1))
                                for a in range(1, q)]))
            for jj in range(1, q - 1) for j in [jj]]

def Lchi(s, q, tbl):
    return mp.mpf(q) ** (-s) * sum(complex(tbl[a]) * mp.zeta(s, mp.mpf(a) / q)
                                   for a in range(1, q) if tbl[a] != 0)

print("=== (1) LGrowth bound ||L(s,chi)|| <= q|s|/Re s  (tree: LGrowth.lean:210) ===")
worst = 0.0
for q in (5, 7, 11):
    for (j, tbl) in prim_chars(q)[:2]:
        for s in (mp.mpc(0.5, 2), mp.mpc(0.5, 14.1), mp.mpc(1.0, 30), mp.mpc(2.0, 5), mp.mpc(0.25, 8)):
            v = abs(Lchi(s, q, tbl))
            bd = q * abs(s) / s.real
            worst = max(worst, float(v / bd))
print(f"    worst ratio over q in {{5,7,11}}, 2 chars each, 5 s-values: {worst:.4f}  "
      f"({'OK' if worst <= 1 else 'FAIL'})")

print("=== (2) A0 on the real zero dataset (60 primitive chars, moduli 5..19) ===")
with open(DATASET, 'rb') as f:
    Z = pickle.load(f)
amax, arg = 0.0, None
for (q, j), rec in Z.items():
    zs = np.asarray(rec['zeros'] if isinstance(rec, dict) and 'zeros' in rec else rec[0])
    for t0 in np.arange(0.0, max(zs) - 1, 0.5):
        cnt = int(np.sum((zs > t0) & (zs <= t0 + 1)))
        val = cnt / np.log(q * (abs(t0) + 3))
        if val > amax: amax, arg = val, (q, j, t0, cnt)
print(f"    max N_chi(t,t+1]/log(q(|t|+3)) = {amax:.3f}  at (q,j,t,count) = {arg}")
print(f"    => measured A0 on the dataset = {amax:.2f}; the classical statement carries an")
print(f"       absolute A0 (Landau-type); [diagnostic only]")
