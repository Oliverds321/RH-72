"""
High-precision C-penalised plateau payoff, corrected.

Kernel W_C(x) = |x| (|x|<=1), C|x| (|x|>1)  -- DISCONTINUOUS at |x|=1.
W_C'' = 2 delta(x) + (C-1)[delta(x-1)+delta(x+1) + delta'(x-1) - delta'(x+1)].

EL on support [-b,b] (1/2 < b < 1), v even, v + W_C*v = mu on supp:
  middle |t| < t0=1-b:  v'' + 2v = 0            -> v = A cos(sqrt2 t)
  edge   t0 < t < b:    with s = t-1/2 in (-d,d), d = b-1/2, w(s)=v(1/2+s):
      w''(s) + 2 w(s) + (C-1)[ w(-s) - w'(-s) ] = 0
  characteristic: rho in roots of (rho^2+2)^2 = (C-1)^2 (1-rho^2),
  i.e. z=rho^2:  z^2 + (4+D^2) z + (4-D^2) = 0,  D = C-1.
  Each z-root gives ONE real solution  w_z(s) = Re[ e^{rho s} + beta e^{-rho s} ],
  beta = -(z+2)/(D(1+rho))   (the pairing beta(rho)beta(-rho)=1 collapses the pair).
Matching: v in C^1 at t0 (2 conds); free boundary v(b)=0 (1 cond) -> root-find b.
mu = v(0) + 2 int_0^b s v(s) ds  (|s|<b<1 so C never enters here);  B = mu / int v.
P = 2 - B.

Independent check: Nystrom trapezoid with kink-aligned grid and AVERAGED kernel
value (1+C)/2*|x| on the discontinuity |t_i-t_j|=1, KKT solve at fixed b.
"""
import numpy as np
from mpmath import mp, mpf, mpc, cos, sin, sqrt, exp, quad, findroot, pi, re, arg, log

mp.dps = 40
SQ2 = sqrt(mpf(2))
MT_exact = mpf(3)/2 - (1/SQ2)*(cos(1/SQ2)/sin(1/SQ2))

# ---------- solver B: closed-form zones ----------
def zone_basis(C):
    C = mpf(C); D = C - 1
    disc = sqrt((4+D**2)**2 - 4*(4-D**2))
    z1 = (-(4+D**2) + disc)/2
    z2 = (-(4+D**2) - disc)/2
    basis = []
    for z in (z1, z2):
        rho = sqrt(mpc(z))              # complex-safe
        beta = -(z+2)/(D*(1+rho))
        def w(s, rho=rho, beta=beta):   # value
            return re(exp(rho*s) + beta*exp(-rho*s))
        def wp(s, rho=rho, beta=beta):  # derivative
            return re(rho*exp(rho*s) - rho*beta*exp(-rho*s))
        basis.append((w, wp, rho, beta))
    return basis

def check_basis(C, d=mpf('0.12')):
    """verify each basis fn satisfies the functional ODE at random points"""
    C = mpf(C); D = C-1
    errs = []
    for (w, wp, rho, beta) in zone_basis(C):
        h = mpf(10)**-10
        for s in (mpf('0.03'), mpf('-0.07'), mpf('0.11')):
            wpp = (w(s+h) - 2*w(s) + w(s-h))/h**2
            res = wpp + 2*w(s) + D*(w(-s) - wp(-s))
            errs.append(abs(res))
    return max(errs)

def solveB_fixed(C, b):
    """EL solution at fixed support [-b,b]; A=1. Returns dict.

    Junction conditions at t0 = 1-b (s0 = -d):
      value:  c1 w1(-d) + c2 w2(-d) = cos(sqrt2 t0)
      slope:  v jumps by -v(b) at the support edge, so (W*v)' jumps by
              (C-1) v(-b) at t0 and v' must jump DOWN by (C-1)v(-b):
              v'(t0+) = v'(t0-) - (C-1) v(b)   [v(b) = w(d) by evenness]
      =>  c1 [w1p(-d) + (C-1) w1(d)] + c2 [w2p(-d) + (C-1) w2(d)]
            = -sqrt2 sin(sqrt2 t0).
    At the free boundary (v(b)=0) this reduces to plain C^1 matching."""
    C = mpf(C); b = mpf(b); d = b - mpf(1)/2; t0 = 1 - b
    (w1, w1p, *_), (w2, w2p, *_) = zone_basis(C)
    M = mp.matrix([[w1(-d), w2(-d)],
                   [w1p(-d) + (C-1)*w1(d), w2p(-d) + (C-1)*w2(d)]])
    r = mp.matrix([cos(SQ2*t0), -SQ2*sin(SQ2*t0)])
    c = mp.lu_solve(M, r)
    c1, c2 = c[0], c[1]
    wfun  = lambda s: c1*w1(s) + c2*w2(s)
    def v(t):
        t = mpf(t); at = abs(t)
        if at <= t0: return cos(SQ2*at)
        return wfun(at - mpf(1)/2)
    edge = wfun(d)
    I1 = quad(lambda s: s*v(s), [0, t0, b])
    mu = v(0) + 2*I1
    Iv = 2*quad(v, [0, t0, b])
    B = mu/Iv
    return dict(B=B, mu=mu, Iv=Iv, edge=edge, v=v, b=b)

def plateau_exact(C, seed=None):
    C = mpf(C)
    f = lambda b: solveB_fixed(C, b)['edge']
    # bracket the sign change of the edge value, then bisect + refine
    lo, hi = mpf('0.502'), mpf('0.98')
    flo = f(lo)
    step = mpf('0.02'); bprev, fprev = lo, flo
    bcur = lo + step
    while bcur < hi:
        fcur = f(bcur)
        if fprev*fcur < 0:
            break
        bprev, fprev = bcur, fcur
        bcur += step
    else:
        raise RuntimeError(f"no bracket for C={C}")
    b = findroot(f, (bprev, bcur), solver='anderson')
    sol = solveB_fixed(C, b)
    return sol

# ---------- solver A: Nystrom with averaged discontinuity ----------
def kernel_avg(x, C):
    ax = np.abs(x)
    out = np.where(ax <= 1.0, ax, C*ax)
    on = np.isclose(ax, 1.0, rtol=0, atol=1e-12)
    out = np.where(on, (1.0+C)/2.0, out)
    return out

def B_nystrom(C, b, K):
    h = 1.0/K
    J = int(round(b*K)); assert abs(J*h-b) < 1e-12
    t = np.linspace(-b, b, 2*J+1)
    w = np.full(2*J+1, h); w[0] = w[-1] = h/2
    M = np.diag(w) + np.outer(w, w)*kernel_avg(t[:,None]-t[None,:], C)
    a = w
    x = np.linalg.solve(0.5*(M+M.T), a)
    mu = 1.0/(a@x)
    return mu, mu*x, t

def kkt_outside(C, sol, npts=30, reach=1.3):
    C = mpf(C); b = sol['b']; v = sol['v']; Iv = sol['Iv']; mun = sol['mu']/Iv
    def h(t):
        t = mpf(t)
        def integrand(s):
            x = abs(t-s); wgt = x if x <= 1 else C*x
            return wgt*v(s)
        pts = sorted(set([float(-b), float(b)] +
                         ([float(t-1)] if -b < t-1 < b else []) +
                         ([float(t+1)] if -b < t+1 < b else [])))
        return quad(integrand, pts)/Iv
    worst = None
    for tt in np.linspace(float(b)+1e-5, float(b)+reach, npts):
        val = h(tt) - mun
        if worst is None or val < worst[1]: worst = (tt, val)
    return worst, float(mun)

if __name__ == "__main__":
    print("basis ODE residual (C=pi^4/18):", mp.nstr(check_basis(pi**4/18), 3))
    print("basis ODE residual (C=2):      ", mp.nstr(check_basis(mpf(2)), 3))

    print("\n=== fixed-b cross-validation, C = pi^4/18 ===")
    C = pi**4/18; Cf = float(C)
    for b in (0.55, 0.6, 0.62):
        BB = solveB_fixed(C, mpf(b))['B']
        row = [f"exact B={mp.nstr(BB,12)}"]
        vals = []
        for K in (100, 200, 400, 800):
            BA, _, _ = B_nystrom(Cf, b, K)
            vals.append(BA)
            row.append(f"K{K}:{BA-float(BB):+.2e}")
        # Richardson
        rich = (4*vals[-1]-vals[-2])/3
        row.append(f"Rich:{rich-float(BB):+.2e}")
        print(f"  b={b}: " + "  ".join(row))

    print("\n=== plateau values (solver B, free boundary) ===")
    targets = [("2", mpf(2)), ("3.3", mpf('3.3')), ("4.4", mpf('4.4')),
               ("pi^4/18", pi**4/18), ("2pi^4/27", 2*pi**4/27),
               ("8", mpf(8)), ("2*pi^4/18", pi**4/9), ("3*pi^4/18", pi**4/6),
               ("20", mpf(20)), ("40", mpf(40))]
    results = {}
    for name, C in targets:
        sol = plateau_exact(C)
        P = 2 - sol['B']
        results[name] = sol
        print(f"  C={mp.nstr(C,10):>12} ({name:>9}):  b*={mp.nstr(sol['b'],10)}  "
              f"lam*={mp.nstr(2*sol['b'],10)}  P={mp.nstr(P,12)}")

    print("\n=== reviewer comparison ===")
    print("  C=2:      mine", mp.nstr(2-results['2']['B'],8), " FW(n=601) 0.7914495, doc 0.7910")
    print("  pi^4/18:  mine", mp.nstr(2-results['pi^4/18']['B'],8), " adversarial 0.7213760/0.7222214, novelty 0.7221")
    print("  2pi^4/27: mine", mp.nstr(2-results['2pi^4/27']['B'],8), " adversarial 0.7099166, novelty 0.7108")
    print("  C=8:      mine", mp.nstr(2-results['8']['B'],8), " doc 0.7067, adversarial 0.7064667")
    print("  C=20:     mine", mp.nstr(2-results['20']['B'],8), " doc 0.6867, adversarial 0.6866161")
    print("  C=40:     mine", mp.nstr(2-results['40']['B'],8), " adversarial 0.6790873")

    print("\n=== Nystrom plateau check, C = pi^4/18 (scan b on grid, Richardson) ===")
    for K in (200, 400, 800):
        best = None
        prev = None
        for J in range(K//2, int(0.9*K)):
            b = J/K
            B, v, t = B_nystrom(Cf, b, K)
            if v.min() < -1e-12:
                break
            prev = (b, B)
        print(f"  K={K}: largest feasible b={prev[0]:.4f}  P={2-prev[1]:.8f}")

    print("\n=== KKT outside-support check ===")
    for name in ("pi^4/18", "2"):
        worst, mun = kkt_outside([c for n,c in targets if n==name][0], results[name], npts=25)
        print(f"  C={name}: min over (b, b+1.3] of h(t)-mu = {float(worst[1]):+.3e} at t={worst[0]:.3f}  (mu={mun:.6f})")
