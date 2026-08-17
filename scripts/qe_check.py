#!/usr/bin/env python3
"""qe_check.py — numerical anchor for the companion notes ends mirror / B_fam bound).

On the REAL family of all 60 primitive characters with prime moduli {5,7,11,13,17,19}
(the project's dataset family, labeling chi_j(a) = e^{2 pi i j ind_g(a)/(q-1)}):

  (1) digamma envelope: |Re psi(1/4 + kappa/2 + i tau/2)| <= log(2 + |tau|/2) + 4.3
      (the archimedean ingredient of B_fam), verified on a tau grid.
  (2) B_fam bound: sum_chi nu_{X,chi}(tau)^2  <=  2 sum_chi mubound_q(tau)^2
        + 2 (1/pi^2)(X + Q^2 - 1) sum_{n<=X} Lambda(n)^2/n,
      verified pointwise at tau in the window, near field, and far field.
  (3) the LS step at fixed tau: sum_chi |sum_n Lambda(n) n^{-1/2-i tau} chi(n)|^2
        <= (X + Q^2 - 1) sum_n Lambda(n)^2/n   (tau-twisted coefficients, unchanged l2 mass).
  (4) far-field shape: B_fam(tau) grows like the log+ term, ratio to (QL)(1 + log+(|tau|/4T))
      bounded.
"""
import numpy as np
from scipy.special import digamma

PRIMES = [5, 7, 11, 13, 17, 19]
Q = 19
T = 10.0
X = 100
LL = np.log(Q * T / (2 * np.pi))          # cal L = log(QT/2pi)

def prim_chars(q):
    """primitive characters mod prime q via smallest primitive root (dataset labeling)."""
    def order(x):
        k, y = 1, x % q
        while y != 1: y = y * x % q; k += 1
        return k
    g = next(u for u in range(2, q) if order(u) == q - 1)
    dlog = {}; x = 1
    for k in range(q - 1): dlog[x] = k; x = x * g % q
    out = []
    for j in range(1, q - 1):                       # j = 0 is principal (imprimitive)
        tbl = np.zeros(q, complex)
        for a in range(1, q): tbl[a] = np.exp(2j * np.pi * j * dlog[a] / (q - 1))
        kappa = 0 if abs(tbl[(q - 1) % q] - 1) < 1e-9 else 1   # parity from chi(-1)
        out.append((q, j, kappa, tbl))
    return out

FAM = [c for q in PRIMES for c in prim_chars(q)]
print(f"family: {len(FAM)} primitive characters, moduli {PRIMES}  (expect 60)")

# Lambda(n) for n <= X
lam = np.zeros(X + 1)
for p in range(2, X + 1):
    if all(p % d for d in range(2, p)):
        pk = p
        while pk <= X: lam[pk] = np.log(p); pk *= p
ns = np.arange(1, X + 1)
sumL2 = np.sum(lam[1:] ** 2 / ns)

def mu_q(q, kappa, tau):
    return (1 / (2 * np.pi)) * (np.log(q / np.pi)
        + np.real(digamma(0.25 + kappa / 2 + 0.5j * tau)))

def P_chi(tbl, q, tau):
    vals = tbl[ns % q]
    return -(1 / np.pi) * np.real(np.sum(lam[1:] * vals * ns ** (-0.5 - 1j * tau)))

print("=== (1) digamma envelope |Re psi(1/4 + k/2 + i t/2)| <= log(2+|t|/2) + 4.3 ===")
taus_env = np.concatenate([np.linspace(0, 50, 300), np.geomspace(50, 5000, 100)])
worst = 0.0
for kap in (0, 1):
    v = np.abs(np.real(digamma(0.25 + kap / 2 + 0.5j * taus_env)))
    env = np.log(2 + taus_env / 2) + 4.3
    worst = max(worst, np.max(v / env))
print(f"    worst ratio over both parities, tau in [0, 5000]: {worst:.4f}  "
      f"({'OK' if worst <= 1 else 'FAIL'})")

print("=== (2) B_fam bound pointwise ===")
print("    tau      sum nu^2      bound        ratio")
for tau in (0.0, 10.0, 15.0, 20.0, 40.0, 100.0, 1000.0):
    s = sum((mu_q(q, kap, tau) + P_chi(tbl, q, tau)) ** 2 for (q, j, kap, tbl) in FAM)
    bmu = sum(((1 / (2 * np.pi)) * (np.log(q / np.pi) + np.log(2 + tau / 2) + 4.3)) ** 2
              for (q, j, kap, tbl) in FAM)
    bound = 2 * bmu + 2 * (1 / np.pi ** 2) * (X + Q ** 2 - 1) * sumL2
    print(f"    {tau:7.1f}  {s:11.3f}  {bound:11.3f}   {s / bound:.4f}  "
          f"{'OK' if s <= bound else 'FAIL'}")

print("=== (3) LS at fixed tau, tau-twisted coefficients ===")
budget = (X + Q ** 2 - 1) * sumL2
for tau in (0.0, 12.5, 333.0):
    s = 0.0
    for (q, j, kap, tbl) in FAM:
        vals = tbl[ns % q]
        s += abs(np.sum(lam[1:] * ns ** (-0.5 - 1j * tau) * vals)) ** 2
    print(f"    tau={tau:7.1f}:  family sum = {s:10.3f}   budget (X+Q^2-1)||a||^2 = "
          f"{budget:10.3f}   ratio {s / budget:.4f}  {'OK' if s <= budget else 'FAIL'}")

print("=== (4) far-field shape: B_fam(tau) / [ c Q LL (1 + log+(|tau|/4T)) ] ===")
cQL = Q * LL
for tau in (15.0, 40.0, 200.0, 2000.0, 2e4):
    s = np.sqrt(sum((mu_q(q, kap, tau) + P_chi(tbl, q, tau)) ** 2
                    for (q, j, kap, tbl) in FAM))
    shape = cQL * (1 + max(np.log(tau / (4 * T)), 0))
    print(f"    tau={tau:9.1f}:  B_fam = {s:9.3f}   Q*LL*(1+log+) = {shape:9.3f}   "
          f"ratio {s / shape:.4f}")
