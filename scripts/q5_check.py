#!/usr/bin/env python3
"""q5_check.py — numerical anchors for the companion notes character-sum arithmetic).

  (1) EXACT aggregated identity (integer precision), Q = 60:
        sum_{q<=Q} sum*_{chi mod q} chi(n) conj(chi(m))
          = sum_{f | (n-m), f<=Q, (f,nm)=1} phi(f) * M_{nm}(Q/f),
      M_{nm}(y) := sum_{r<=y, (r,nm)=1} mu(r).      [q=1 included: contributes 1]
  (2) crude bound |S| <= Q tau(|n-m|)  (no PNT).
  (3) the zone sum S(Y) = sum_{n != m <= Y} Lambda Lambda (nm)^{-1/2} tau(|n-m|):
      measured growth exponent j in S(Y) ~ Y (log Y)^j  (crude proof gives j <= 3).
  (4) delta'-calibration: relative error Q^{-delta'} L^{j-2} at delta' = K loglog Q/log Q:
      K = 1 is NOT enough at the proved j = 3; K = 2 is safe.  (This is the write-through's
      sharpening of RAMANUJAN_NOTE — its own Sec 5 item 2 warned exactly here.)
"""
import numpy as np
from sympy import totient, mobius, divisors, gcd, factorint
from itertools import product

# ---------- full character groups (adapted from work/ramanujan_check.py, validated there) ----
def _pp_chars(p, a):
    q = p ** a
    units = [u for u in range(1, q) if u % p != 0] if p != 2 else [u for u in range(1, q, 2)]
    if p == 2 and a == 1: return [{1: 1.0 + 0j}]
    if p == 2 and a >= 3:
        o5 = 2 ** (a - 2); exps = {}
        for e1 in range(2):
            for e5 in range(o5):
                u = (pow(q - 1, e1, q) * pow(5, e5, q)) % q
                exps[u] = (e1, e5)
        return [{u: np.exp(2j * np.pi * (j1 * e1 / 2 + j5 * e5 / o5))
                 for u, (e1, e5) in exps.items()} for j1 in range(2) for j5 in range(o5)]
    n_u = len(units)
    def order(x):
        k, y = 1, x % q
        while y != 1: y = y * x % q; k += 1
        return k
    g = next(u for u in units if order(u) == n_u)
    dlog = {}; x = 1
    for k in range(n_u): dlog[x] = k; x = x * g % q
    return [{u: np.exp(2j * np.pi * j * dlog[u] / n_u) for u in units} for j in range(n_u)]

def char_group(q):
    fac = factorint(q)
    comp = [_pp_chars(p, a) for p, a in fac.items()]
    mods = [p ** a for p, a in fac.items()]
    chars = []
    for combo in product(*comp):
        tbl = np.zeros(q, complex)
        for n in range(1, q):
            if gcd(n, q) != 1: continue
            v = 1.0 + 0j
            for c, m in zip(combo, mods):
                v *= c[n % m] if m > 1 else 1.0
            tbl[n] = v
        chars.append(tbl)
    return chars

def conductor(q, tbl):
    for f in sorted(divisors(q)):
        if all(abs(tbl[a % q] - 1) < 1e-9 for a in range(1, q + 1, f) if gcd(a, q) == 1):
            return f
    return q

Q = 60
PRIM = {q: [t for t in char_group(q) if conductor(q, t) == q] for q in range(3, Q + 1)}
# anchor identity for the builder itself (rule 11): family count
tot = sum(len(v) for v in PRIM.values()) + 1
print(f"builder anchor: sum phi*(q), q<=60 = {tot}  vs (18/pi^4) 60^2 = {18/np.pi**4*3600:.1f}")

def M_coprime(y, k):
    return sum(int(mobius(r)) for r in range(1, int(y) + 1) if gcd(r, k) == 1)

print("=== (1) exact aggregated identity + (2) crude bound, Q = 60 ===")
allok = True
for (n, m) in [(2, 3), (2, 7), (3, 13), (5, 29), (11, 13), (7, 37), (4, 9), (25, 27), (2, 32)]:
    S = 1.0 + 0j                                     # q = 1 term
    for q in range(3, Q + 1):
        if gcd(n * m, q) != 1: continue
        S += sum(t[n % q] * np.conj(t[m % q]) for t in PRIM[q])
    rhs = sum(int(totient(f)) * M_coprime(Q // f, n * m)
              for f in divisors(abs(n - m)) if f <= Q and gcd(f, n * m) == 1)
    crude = Q * len(divisors(abs(n - m)))
    ok = abs(S - rhs) < 1e-6 and abs(S) <= crude
    allok &= ok
    print(f"    (n,m)=({n:2d},{m:2d}): S = {S.real:9.3f}  exact-formula = {rhs:5d}  "
          f"crude Q*tau = {crude:4d}  {'OK' if ok else 'FAIL'}")
print("    identity + bound:", "ALL OK" if allok else "FAILURES")

print("=== (3) zone sum S(Y) growth ===")
def zone_sum(Y):
    lam = np.zeros(Y + 1)
    for p in range(2, Y + 1):
        if all(p % d for d in range(2, int(p ** 0.5) + 1)):
            pk = p
            while pk <= Y: lam[pk] = np.log(p); pk *= p
    idx = np.nonzero(lam)[0]
    # tau table
    taut = np.zeros(Y + 1, dtype=int)
    for d in range(1, Y + 1): taut[d::d] += 1
    s = 0.0
    for i, n in enumerate(idx):
        for m in idx[i + 1:]:
            s += 2 * lam[n] * lam[m] / np.sqrt(n * m) * taut[m - n]
    return s
prev = None
print("    Y        S(Y)       S/Y      S/(Y logY)   S/(Y logY^2)   local exponent j")
for Y in (250, 500, 1000, 2000, 4000):
    s = zone_sum(Y)
    lY = np.log(Y)
    if prev:
        j = (np.log(s / Y) - np.log(prev[1] / prev[0])) / (np.log(lY) - np.log(np.log(prev[0])))
    print(f"    {Y:5d}  {s:9.1f}  {s/Y:8.3f}   {s/(Y*lY):8.4f}     {s/(Y*lY**2):8.5f}"
          + (f"        {j:5.2f}" if prev else ""))
    prev = (Y, s)

print("=== (4) delta'-calibration: relative = exp(-K loglogQ) * L^(j-2), proved j = 3 ===")
print("    logQ     K=1        K=1.5      K=2")
for logQ in (50, 230.3, 1000, 10000):
    ll = np.log(logQ)
    row = [np.exp(-K * ll) * logQ ** 1 for K in (1, 1.5, 2)]   # L^{j-2} = L at j=3
    print(f"    {logQ:7.1f}  {row[0]:9.3f}  {row[1]:9.3f}  {row[2]:9.4f}")
print("    (K=1: Theta(1), does NOT vanish; K=1.5 slow; K=2: -> 0 like 1/logQ.  At the")
print("     measured j from (3) the margins improve; the LEMMA adopts K = 2 as proved-safe.)")
