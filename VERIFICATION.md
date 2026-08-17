# Verifying the result, top to bottom

Everything numerical in the paper traces to a script in `scripts/` with its reference
output in `logs/`. Requirements: Python 3.10+, `numpy scipy mpmath`.

## 1. The headline constant (the variational layer, §§10–11 of the paper)

| script | what it establishes | reference log |
|---|---|---|
| `payoff_hp2.py` | **P = 0.7212835668…** to ten digits, by closed-form Euler–Lagrange zones cross-checked against an independent Nyström discretisation, with the free-boundary and KKT (complementarity) conditions verified — plus the dyadic and parity-family constants | `log_payoff.txt` |
| `payoff_final2.py` | direct functional evaluation of the same constants in independent float arithmetic | (same log) |
| `gate3_so16.py` | calibration gate on a *third party's published number*: the functional, run on Sono's GRH kernel min(\|α\|,1) at λ = 2, reproduces his published M = 0.93228262 to ten digits | `log_gate3.txt` |
| `qg_check.py` | a deliberately coarse independent grid anchor (4–5 digits; its printed v ≥ 0 flag is a discretisation artifact, absent in the fine solves — see paper §1.6) | `log_qg.txt` |

Three calibration gates must pass in any reimplementation: the Montgomery–Taylor value
0.672500703679412 (closed form), the flat-profile value 4/3, and Sono's 0.93228262.

## 2. The finite-Q table (paper §10.4, Table 2)

| script | what it establishes | reference log |
|---|---|---|
| `table2_check.py` (uses `fb.py`, `budget_q.py`) | the printed effective-proportion table at the design of record, under the conservative ramp accounting, with all budget rows carried; also the sharper alternative reading and the thresholds | `log_table2.txt` |
| `trap_wfactor_anchor.py` | two quoted figures anchored: the 0.3693 "trap" value that Corollary 3 warns against (in-zone coefficient doubled), and the §6 declined-weight figures (⟨q/φ(q)⟩ = 1.2965, C_w = 4.174, +1.30/+1.01) | `log_trap_wfactor.txt` |

## 3. Lemma-level numerical anchors (diagnostics, not evidence — they check conventions and transcription)

| script | paper section | what it checks | log |
|---|---|---|---|
| `q7_check.py` | §4 | the Parseval identity in the paper's Fourier convention (to 2×10⁻¹⁴) | `log_q7.txt` |
| `q5_check.py` | §5 | the aggregated character-sum identity, exact to integer precision, 9/9 cases | `log_q5.txt` |
| `qt_check.py` | §7 | the Gevrey taper transform bound and the tail's closing ladder | `log_qt.txt` |
| `qe_check.py` | §8 | the family mean-square envelope on real L-function data | `log_qe.txt` |
| `qf_check.py` | §9 | the q-uniform growth bound and the local-count shape on the 60-character dataset in `data/` | `log_qf.txt` |
| `r4_design_check.py` | §10 | the budget model across designs (historical: includes the superseded designs for comparison) | `log_r4_design.txt` |

## 4. What cannot be verified from this repository alone

The six imported statements H1–H6 (paper §2.1) are theorems of
[zeta-23-lean](https://github.com/anthropics/zeta-23-lean); verify them there (`lake build`
+ its `#print axioms` comparator). The detailed prose derivations behind the paper's
sections (the companion working notes) and the full multi-round review record are packaged
separately and available on request; they will be added to this repository.
