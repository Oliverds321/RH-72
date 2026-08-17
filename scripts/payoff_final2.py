"""Final headline numbers; direct functional evaluation safeguard in float precision."""
from mpmath import mp, mpf, cos, sin, sqrt, quad, pi, nstr
import numpy as np
from scipy import integrate
import payoff_hp2 as P2

mp.dps = 40
MT = P2.MT_exact

def direct_B_float(C, sol):
    """B(v) evaluated directly with scipy on the closed-form v (float precision)."""
    Cf = float(C); b = float(sol['b']); t0 = 1-b
    vv = sol['v']
    v = lambda t: float(vv(t))
    pts = [-b, -t0, 0.0, t0, b]
    I2 = integrate.quad(lambda t: v(t)**2, -b, b, points=pts, limit=200)[0]
    def inner(t):
        f = lambda s: (abs(t-s) if abs(t-s) <= 1 else Cf*abs(t-s))*v(s)
        sp = sorted(set(pts + ([t-1] if -b < t-1 < b else [])
                            + ([t+1] if -b < t+1 < b else []) + [t]))
        return integrate.quad(f, -b, b, points=[p for p in sp if -b <= p <= b], limit=200)[0]
    Icross = integrate.quad(lambda t: v(t)*inner(t), -b, b, points=pts, limit=100)[0]
    return (I2 + Icross)/float(sol['Iv'])**2

print("MT anchor       :", nstr(MT, 20))
for name, C in [("pi^4/18", pi**4/18), ("2pi^4/27", 2*pi**4/27),
                ("2*pi^4/18", pi**4/9), ("3*pi^4/18", pi**4/6)]:
    sol = P2.plateau_exact(C)
    Pv = 2 - sol['B']
    Bd = direct_B_float(C, sol)
    print(f"\nC = {nstr(C,18)}  ({name})")
    print("  lambda* =", nstr(2*sol['b'], 18))
    print("  P (EL)  =", nstr(Pv, 18))
    print("  P (dir) =", f"{2-Bd:.12f}", "   |EL-dir| =", f"{abs(float(sol['B'])-Bd):.2e}")
    print("  gain    = +%s pp" % nstr(100*(Pv-MT), 8))

for C in (40, 100):   # C -> infinity approach to MT from above
    sol = P2.plateau_exact(mpf(C))
    print(f"\nC={C}: P = {nstr(2-sol['B'],12)}  (-> MT from above)  lam*={nstr(2*sol['b'],8)}")
