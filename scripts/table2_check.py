#!/usr/bin/env python3
"""table2_check.py — THE script that computes DRAFT §10.4's Table 2 (v0.5).

Provenance: the free-boundary solver (fb.py) and budget model (budget_q.py) are the
round-7/8 adversarial review's reproduction code, shipped in the round-8 review bundle,
adopted here after verification against the paper's own figures. One change from the
reviewer's model: the ends row is charged at OUR budget's coefficient 6·ℒ·logℒ/T
(DRAFT §8 / Lemma 8.2), not the reviewer's 4 (immaterial at 2 d.p. everywhere).
The four mandatory calibration gates run at import (fb.gates()).

Outputs:
 (1) Table 2 AS PRINTED in v0.5: design (r,K) = (3.5,3), ramp row charged at the
     conservative r₂ = 6w/L of §10.3, curvature-honest zone cost (measured P(δ)
     secant at the design offset), all minor rows carried.
 (2) The sharper 4w/L variant quoted in §10.4's follow-on sentence.
 (3) Thresholds for both readings; the quantified gap between them.
 (4) w* at the design of record (the §10.3 parenthetical) and the rate constant.
"""
from math import log
import fb, budget_q as bd

lq = lambda n: n * log(10)
REC = bd.REC
print("=== gates ===")
fb.gates()

def table(c_ramp, label):
    print(f"\n=== Table 2 at ramp {label} (design (3.5,3), curvature, minor rows in) ===")
    hdr = "%9s %8s %7s %7s %8s %8s %8s %8s %9s %9s | %8s %8s"
    print(hdr % ("Q","LL","lam*","delta","zone","ramp","buffer","ends","minor","TOTAL","P_eff","w*"))
    out = {}
    for n in (25, 100, 300, 1000, 10000):
        d = bd.design(lq(n), r=3.5, K=3.0, c_ramp=c_ramp, curvature=True, minor=True)
        rr = d['rows']; out[n] = d
        print(hdr % ("1e%d"%n, "%.1f"%d['LL'], "%.4f"%d['lam'], "%.4f"%d['delta'],
              "%.4f"%rr['zone'], "%.4f"%rr['ramp'], "%.4f"%rr['buffer'], "%.5f"%rr['ends'],
              "%.5f"%rr['minor'], "%.4f"%d['total'], "%+.4f"%d['P_eff'], "%.2f"%d['w']))
    for tgt, name in ((0.0,"non-vacuous"), (0.5,"P_eff>0.5"), (REC,"P_eff>record"), (2/3.,"P_eff>2/3")):
        x = bd.threshold(tgt, r=3.5, K=3.0, c_ramp=c_ramp, curvature=True, minor=True)
        print("   %-14s from 1e%.0f" % (name, x/log(10)))
    return out

t6 = table(6.0, "6w/L  (r2 = 6w/L, THE PRINTED TABLE)")
t4 = table(4.0, "4w/L  (the sharper reading)")

print("\n=== the gap, quantified (for §10.4's follow-on sentence) ===")
for n in (25, 100, 300, 1000):
    print("   1e%-5d  6w/L %+.4f   4w/L %+.4f   gain %+.4f" %
          (n, t6[n]['P_eff'], t4[n]['P_eff'], t4[n]['P_eff'] - t6[n]['P_eff']))
for tgt, name in ((0.0,"non-vacuous"), (0.5,">0.5"), (REC,">record")):
    a = bd.threshold(tgt, r=3.5, K=3.0, c_ramp=6.0, curvature=True, minor=True)
    b = bd.threshold(tgt, r=3.5, K=3.0, c_ramp=4.0, curvature=True, minor=True)
    print("   %-12s 6w/L 1e%-5.0f 4w/L 1e%-5.0f  (%.2f orders)" %
          (name, a/log(10), b/log(10), (a-b)/log(10)))

print("\n=== w* and the rate constant at the design of record ===")
d = bd.design(lq(100), r=3.5, K=3.0, c_ramp=6.0, curvature=True, minor=True)
print("   w* at 1e100, (3.5,3), 6w/L objective : %.2f   (D0 = %.3e, wD0/LL^2 = %.2f)" %
      (d['w'], d['D0'], d['wD0']/d['LL']**2))
d3 = bd.design(lq(100), r=3.0, K=3.0, c_ramp=4.0, curvature=True, minor=True)
print("   [historical: w* at r=3, 4w/L objective: %.2f — the region the old '~2.3' was quoted from]" % d3['w'])
for cr in (4.0, 6.0):
    dd = bd.design(2302600.0, r=3.5, K=3.0, c_ramp=cr, curvature=True, minor=True)
    print("   rate const total*LL/logLL at logQ=2.3e6, ramp %.0fw/L : %.3f  [s(r+K) = %.3f]" %
          (cr, dd['rate_const'], bd.S_LIN*6.5))
