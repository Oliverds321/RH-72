#!/usr/bin/env python3
"""gate3_so16.py -- THIRD calibration gate: Sono [So16]'s GRH benchmark kernel
F(y) = min(|y|,1) at lambda = 2 (supp v in [-1,1]).  Published value M = 0.93228262.
B(v) = int v^2 + iint F(t-s) v(t)v(s);  minimise over int v = 1 (sign-unconstrained
linear solve; positivity of the optimiser checked a posteriori).  Richardson in n.
"""
import numpy as np
def solve(n):
    # midpoint cells on [-1,1]
    h = 2.0/n; t = -1.0 + h*(np.arange(n)+0.5)
    D = np.abs(t[:,None]-t[None,:])
    K = np.minimum(D,1.0)*h*h          # kernel term, cell-integrated (midpoint)
    A = np.eye(n)*h + K                # int v^2 -> h * v_i^2 (midpoint)
    one = np.ones(n)*h
    # minimise v^T A v  s.t.  one^T v = 1  ->  A v = mu * one/ (2) ; scale-free:
    x = np.linalg.solve(A, one)
    v = x/(one@x)
    B = v@A@v
    return B, v
Bs = {}
for n in (400, 800, 1600, 3200):
    B, v = solve(n)
    Bs[n] = B
    print(f"  n={n:5d}  B = {B:.10f}  P = {2-B:.10f}  min v/max v = {v.min()/v.max():+.2e}")
# Richardson: convergence is second order in h (measured ratio ~4), so B_ext = B_n + (B_n - B_{n/2})/3
r1 = Bs[3200] + (Bs[3200]-Bs[1600])/3; r2 = Bs[1600] + (Bs[1600]-Bs[800])/3
print(f"  Richardson (2nd order):  P = {2-r1:.10f}  (prev {2-r2:.10f})")
print(f"  target [So16] M = 0.93228262   delta = {abs((2-r1)-0.93228262):.2e}")
