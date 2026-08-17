# 72.128% of the Generalized Riemann Hypothesis, on average, unconditionally

> The **Generalized Riemann Hypothesis** predicts that every zero of every Dirichlet
> L-function lies on the critical line. This work proves — unconditionally, with no
> hypothesis — that at least **72.128%** of them do, and are simple besides, averaged over
> *every* primitive character of *every* modulus q ≤ Q. That is the largest fraction of
> GRH ever established for this family, by 8 to 15 points over the standing records.
>
> Derived by an AI system and put through repeated adversarial review before you ever saw
> it. **This repository is the verification core:
> the paper, and every script, log, and dataset behind every number in it — enough to check
> the result top to bottom.**

*(To be precise about what is and is not claimed: this does **not** prove the Riemann
Hypothesis or GRH, and it makes no claim about any individual L-function — a zero-density
subfamily could in principle carry all its zeros off the line. It proves the stated
proportion **on family average**, which is the form in which the field's strongest
unconditional knowledge of GRH exists.)*

---

## The result

**Theorem 1.** For the family of all primitive Dirichlet characters of modulus q ≤ Q, in
ordinate windows [T, 2T] with T = (log Q)^{3+ε}: unconditionally, at least

**P = 0.7212835668… − O(log log Q / log Q)**

of the zeros are simple and lie on the critical line.

| statement | family | value | previous record |
|---|---|---|---|
| Theorem 1 | all primitive, q ≤ Q, equal weight | **0.72128** | — (no prior result for this family) |
| Corollary 2 | dyadic conductor q ∈ (Q/2, Q] | **0.70992** | 0.56 — Conrey–Iwaniec–Soundararajan (2013), at r ≥ 6, up to window shape and weight |
| Corollary 3 | even primitive, dyadic | **0.69194** | 0.6044 — Sono (2021) |
| Corollary 3′ | **odd** primitive, dyadic | **0.69194** | — (no prior result at all) |

For calibration: the strongest *published* per-character unconditional proportion is Wu's
0.4074 (2019); the strongest per-character result of any kind is 0.67250, from the preprint
this work builds on. The engine is a one-sided use of the classical multiplicative large
sieve — which is lossless precisely when the family is larger than the polynomial — in
place of the GRH pair-correlation asymptotic, with a C-penalised variational problem
converting the sieve constant into a proportion.

## The foundation

This work builds on **[R]** — the August 2026 Anthropic preprint whose 0.67250
per-character theorem is the strongest result of its kind — by importing exactly **six
statements (H1–H6, §2.1 of the paper)**, each a `sorry`-free theorem in [R]'s public
Lean 4 artifact ([zeta-23-lean](https://github.com/anthropics/zeta-23-lean)), whose axiom
audit reports only the three standard axioms. [R] is, at this writing, unposted and
unrefereed; the paper says so up front, and §2.1 isolates every import precisely so the
dependency can be judged on its own terms — a reader who prefers to treat the result as
conditional on H1–H6 can do exactly that, and the Lean formalization (below) will make the
entire chain mechanical. Everything *beyond* those six statements — the family theorem,
the sieve mechanism, the two-zone corollaries — is this work's own, and is what the review
record below has been testing.

## Reproduce it

```
pip install numpy scipy mpmath
python scripts/payoff_hp2.py      # the headline constant to ten digits, with KKT checks
python scripts/gate3_so16.py      # the third calibration gate: Sono's published 0.93228262
python scripts/table2_check.py    # the finite-Q orientation table, both accountings
```

Every number quoted in the paper traces to a script in [`scripts/`](scripts/) and a log
shipped beside it; the 5,230-zero dataset used by the diagnostic anchors is in
[`data/`](data/).

## Authorship and AI disclosure

**Author: Oliver D'Souza.** The mathematics, numerics, and exposition were produced
substantially by Claude (Anthropic) under the author's direction; verification was by
process — the adversarial review record above, the calibration gates, and the Lean track —
rather than by the author's own line-by-line reading, and the paper's disclosure statement
(§1.6) says exactly that in arXiv's own vocabulary. No AI system is listed as an author.
The foundation [R] is itself an AI-generated work. This repository exists so that nobody
has to take anyone's word — human or machine — for any of it.

## Verification

The work was checked through multiple rounds of independent, adversarial AI review —
including sessions given no project history at all — with every review and point-by-point
response preserved in a companion archive available to reviewers and collaborators. The
construction has survived every round; the headline constants have been reproduced by
seven independent implementations; and the variational layer reproduces two third-party
*published* constants — Montgomery–Taylor's 0.672500703679 and Sono's 0.93228262 — as
mandatory calibration gates that every shipped solver must pass.

## Roadmap

- **Now:** external review continues; new reviews are folded into the record as they land.
- **Next:** arXiv (math.NT) posting — the endorsement process is under way.
- **In progress:** the Lean 4 formalization — a skeleton freezing every lemma statement,
  filled in parallel tracks, held to the same standard [R] set: `sorry`-free,
  axiom-audited, machine-checked end to end. It will be added to this repository as it
  develops, along with the full review archive.

## Reviews welcome

Reviews, comments, error reports, and collaboration are all welcome — via issues here or
by contact with the author.

---

*License: Apache-2.0. Cite via `CITATION.cff`. The paper: [`paper/`](paper/).*
