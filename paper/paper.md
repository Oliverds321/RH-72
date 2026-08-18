# Simple zeros on the critical line for the primitive Dirichlet family: 72.128% of the Generalized Riemann Hypothesis, on average, unconditionally

**Preprint — 2026-08-18 (v0.8).**

**Author: Oliver D'Souza.**

**Abstract.** We prove that, unconditionally, at least 0.72128 of the zeros of Dirichlet
L-functions L(s, χ), averaged over all primitive characters of modulus q ≤ Q in the ordinate
window [T, 2T] with T = (log Q)^{3+ε}, are simple and lie on the critical line, with an
explicit rate of convergence. On the corresponding families this improves the unconditional
records — 56% over all primitive characters at q ≍ Q (Conrey–Iwaniec–Soundararajan [CIS2],
at their window floor T ≥ (log Q)⁶, which our statement reaches at r ≥ 6, and up to window
shape: [CIS2] counts |γ| ≤ T where we count dyadic blocks) and 0.6044 over
the even primitive characters (Sono [So21]) — by 14.99 and 8.75 points; both partners carry
a smooth conductor weight where our count is sharp and unweighted. Per character, the
strongest PUBLISHED unconditional proportion is Wu's 0.4074 [Wu19]; the 0.672500703679 of
the August 2026 preprint [R], on which our architecture builds and which is at this writing
unposted to any preprint server and unrefereed, is likewise per character, in a range of T
relative to q disjoint from ours. The q ≤ Q family of Theorem 1, in which every primitive
character of every modulus up to Q is counted with equal weight, has no counterpart in the
literature. Our results consume six statements imported from [R], isolated in §2.1, each a
sorry-free theorem in [R]'s public Lean 4 artifact — at [R]'s own parameter class; their
use at THIS paper's parameters (bandwidth λ > 1) rests on the file-level audit of §9, and
the re-parameterisation that would certify it mechanically is an open obligation of the
formalization track. The input is one-sided: the multiplicative large
sieve — lossless precisely when the family is larger than the polynomial, which holds here
at every bandwidth below 2 — replaces the GRH pair-correlation asymptotic at the cost of an
explicit constant, and a C-penalised variational problem converts that constant into a
proportion. No asymptotic large sieve and no prime-correlation asymptotic is used, and no
unproved conjecture is assumed: beyond classical prime number theory, every input is one of
the six imported statements of §2.1 — proved theorems of [R], consumed here as premises.
The per-character architecture (the rank–trace certificate and the Weil-form grid Gram) is
imported from the August 2026 preprint [R]; everything family-level is new.


---

## 1. Introduction

### 1.1 The result

Let 𝔉_Q denote the set of primitive Dirichlet characters of modulus q ≤ Q, and for χ ∈ 𝔉_Q
let N_χ(T, 2T) count the nontrivial zeros ρ = β + iγ of L(s, χ) with T < γ ≤ 2T, and
N^s_{0,χ}(T, 2T) those that are simple and on the critical line. The architecture we build on
does not assume the zeros lie on the line, nor that γ is real for β = 1/2 — handling possible
off-line zeros unconditionally is the point.

**Theorem 1.** Fix ε > 0. There are Q₀(ε) and an effectively computable c(ε) such that for
Q ≥ Q₀(ε) and T = T(Q) = (log Q)^{3+ε}, unconditionally,

>   Σ_{χ∈𝔉_Q} N^s_{0,χ}(T, 2T) ≥ (P − c(ε)·(log log Q)/(log Q)) · Σ_{χ∈𝔉_Q} N_χ(T, 2T),

where **P = 0.7212835668…** is the value of an explicit variational problem (§10) at
C = π⁴/18. (The same statement holds at T = (log Q)^{r+ε} for every r ≥ 3, with c = c(r, ε)
effectively computable; a computed — not proved — value of the construction's budget constant
is in the remark of §10.3, outside this statement.)

**Corollary 2 (dyadic).** The same holds for the dyadic subfamily q ∈ (Q/2, Q] with
**P_dyad = 0.7099167448…** (C = 2π⁴/27, λ* = 1.1931581210), for every T = (log Q)^{r+ε},
r ≥ 3. The comparison partner in the literature is Conrey–Iwaniec–Soundararajan [CIS2] (all
primitive characters, q ≍ Q, no parity restriction, weighted Ψ(q/Q)/φ(q) — the weight
largely cancels in a ratio of two identically weighted sums, though we do not claim it
exactly does): their theorem gives 14/25 = 56%
simple-and-on-line (58.65% announced at optimised parameters), in the window |γ| ≤ T with
floor T ≥ (log Q)⁶. At r ≥ 6 our statement enters [CIS2]'s window class (asymptotic budget
constant ≈ 4.6 rather than ≈ 3.3; at finite Q the r = 6 design is in fact STRONGER, since
the buffer and ramp rows fall faster in T than the zone row rises — §10.4), and there the
comparison is like-for-like up to window shape ([CIS2] counts all |γ| ≤ T; we count the
dyadic block [T, 2T] — the two statements are never literally identical, since this method
needs T ≥ (log Q)³ and cannot cover the low blocks): **+14.99 points** over 14/25 (+12.34
over the announced optimum). Theorem 1's q ≤ Q family has no comparison partner in the
literature.

**Corollary 3 (even primitive, dyadic — the Sono-comparable statement).** Sono's family
[So21] is the EVEN primitive characters of dyadic conductor. The even subfamily carries
exactly half the character mass — |𝔉_even| = ½|𝔉| + O(Q), provably: #even-prim(q) =
(φ*(q) + S(q))/2 with S(q) := Σ_{χ prim mod q} χ(−1) multiplicative (S(p) = −1 for odd p,
S(p^a) = 0 for odd p and a ≥ 2; S(2) = 0, S(4) = −1, S(2^a) = 0 for a ≥ 3), so S is
supported on the squarefree-odd q and 4·(squarefree-odd) and Σ_{q≤Q} S(q) = O(Q). The
subfamily passage is a TWO-ZONE argument, and the two zones use different mechanisms —
positivity alone, applied in both zones, would double the in-zone coefficient as well and
prove only 0.3693 (§11 at in-zone weight 2|α|), which is NOT this corollary. IN-ZONE
(n, m ≤ Y = Q^{1−δ′}), the even-subfamily main term is computed exactly by parity
projection, Σ_{χ even} = ½Σ*_χ(1 + χ(−1)): the 1-half is half the full-family form, in step
with the halved denominator, and the χ(−1)-half pairs n ≡ −m (mod q) and is governed by the
aggregated identity in its n + m form, Lemma 5.2′ (stated in §5) — in-zone
n + m ≤ 2Q^{1−δ′} ≪ Q, so Lemma 5.3′'s divisor bound applies and it is negligible at the
same zone-edge calibration. The in-zone coefficient of the variational problem is therefore
UNCHANGED. OUT-ZONE (the C-penalised region), where n reaches X = Q^λ and n + m outruns
every modulus — exactly where the projection route fails by Q^{λ−1} (§12.3) — the sieve
input passes to the subfamily by positivity of its summands: the budget N + Q² − 1 is
unchanged, the denominator halves, and the OUT-ZONE constant doubles, C → 2C. The error
rows (ends, cross, tail) are one-sided or absolute-valued, so the subfamily is bounded by
the full family at the cost of a factor ≤ 2 against the halved denominator — each is priced
with an order-of-magnitude margin in the budget, and doubling them moves nothing. Finally
the denominator is a ZERO count, not the character count just proved: 𝒩_even = ½𝒩 +
O(QTℒ), by partial summation of N_χ(T, 2T) = (T/2π)(log(qT/2π) + 2 log 2 − 1) + O(log qT)
against
Σ_{q≤x}S(q) = O(x) — relative error O(1/Q). (Parity uniformity holds where it is touched
at all: μ_q's parity dependence is exactly π Re cot(π(¼ + iτ/2)) = O(e^{−πτ}) by the
reflection formula — beyond floating point at τ ≍ T.) The argument with out-zone constant
2C gives, for even primitive χ with q ∈ (Q/2, Q]:
**P_even,dyad = 0.6919434301…** (C = 4π⁴/27, λ* = 1.1015998422) — **+8.75 points over the
unconditional family record 0.6044 [So21]**. (Even primitive over q ≤ Q: 0.6980745436 at
C = π⁴/9.) The identical argument with the projector ½Σ*(1 − χ(−1)) gives the ODD
primitive family at the same constants: **P_odd,dyad = 0.6919434301** (q ∈ (Q/2, Q];
0.6980745436 over q ≤ Q) — a class in which the literature contains no result at all
([So21] is stated for even characters only, in its (1.8)). Remark: a parity-refined large
sieve at budget ½(N + Q² − 1) would restore C in the out-zone, but the n ≡ −m orthogonality
route — which IS what the corollary uses in-zone — fails OUT-ZONE at λ > 1 by Q^{λ−1}, the
same mechanism as the cross term of §4, and no parity-restricted sieve with the ½ constant
appears to be available; mirroring Sono's smooth weight W ≤ 1 on [1, 2] would further
inflate C by 3/(2∫W(x)x dx).

Both statements are family averages: they say nothing about any individual L(s, χ), and a
zero-density subfamily could in principle carry all its zeros off the line. This is the same
statement class as [CIS2] and [So21] ([So16] is GRH-conditional and counts low-lying zeros
weighted by |Φ̂(iγ)|² for a fixed Φ, so its count concentrates at |γ| ≪ 1/log Q — a
different class on both grounds; the harmonic 1/φ(q) weight is NOT a class separator, since
[CIS2] itself carries Ψ(q/Q)/φ(q) — see Corollary 2's statement).

### 1.2 What was known

Per character, the strongest PUBLISHED unconditional proportion of simple-and-on-line zeros
of Dirichlet L-functions is Wu's **0.4074** [Wu19] (with 0.4172 on the line), for every
primitive χ mod q with log q = o(log T); for the family of all primitive characters of a
single FIXED modulus, Dickinson [Di24] gives 38.2% on the line, at T ≫ q^ε; for ζ alone the corresponding constants are
Pratt–Robles–Zaharescu–Zeindler's 0.407511 and 0.417293 [PRZZ] — the latter is the "more
than five-twelfths" quoted as 41.7% in the recent expository literature [GS26]. Higher
still, but not published: the preprint [R], dated 11 August 2026 and revised in place
13 August 2026, proves for ζ and — its **Theorem B** — for every FIXED primitive χ that at
least 3/2 − (1/√2)cot(1/√2) = **0.672500703679…** of the zeros in dyadic windows [T, 2T]
are simple and on the critical line, at full bandwidth λ = 1 with the Montgomery–Taylor
window and T ≥ T₀ depending on the character ([R]'s Remark 6.1: no λ < 1 improves the
constants). It is Lean-formalized (the artifact retains the lettering of the 10 August
version, in which the Dirichlet case is Theorem E); it is not posted to a preprint server,
and we are aware of no journal submission, refereeing, erratum or withdrawal. The constant
itself is not new: it appears in [BGST2, Thm 2] under a narrow-vertical-box hypothesis, and
[R]'s contribution is removing that hypothesis. We make no claim about [R]'s reception.
**Our construction is not a corollary of [R]'s Theorem B**: that theorem is per character,
with T₀ and the error constants depending on q — [R] states the Dirichlet case with O_q(·)
errors and claims no uniformity — and it supplies no family statement at T = polylog(Q);
the uniformity in q is the work, and it is where every difficulty of this paper lives. The
same disjointness applies to [Wu19], whose hypothesis log q = o(log T) is the exact reverse
of our range, in which log T ≍ r log log Q = o(log q) for almost every member of the family.
We record for orientation that our family average 0.72128 exceeds both per-character
constants (by 4.88 and 31.39 points), but the statements are of different kinds — a
per-character bound implies the family average in its own range, not conversely — and we
make no claim about any individual L(s, χ).

On family average the unconditional records are, by class: over the EVEN primitive
characters of dyadic conductor, Sono's 0.6044 simple-and-on-line (0.6107 on-line) [So21];
over ALL primitive characters at q ≍ Q with the harmonic weight Ψ(q/Q)/φ(q),
Conrey–Iwaniec–Soundararajan's 14/25 = 56% (58.65% announced at numerically optimised
parameters) [CIS2]; and, for the weaker quantity "simple, critical line not required" over
that same weighted family, Wu's 0.60261, or 0.66433 under GRH [Wu16]. Under GRH the
low-lying benchmark is [So16]'s 93.22% — a different quantity again (§1.1). On the
corresponding families, Corollary 2 exceeds [CIS2] by **+14.99 points** (at r ≥ 6, inside
[CIS2]'s (log Q)⁶ window floor and up to window shape and the harmonic weight — see the
corollary) and Corollary 3 exceeds [So21] by **+8.75 points**. Since simple-and-on-line
implies simple, Corollary 2 also gives ≥ 0.70992 simple, exceeding [Wu16]'s GRH figure
unconditionally by +4.56 points, modulo his harmonic weight. (Each corollary pays its own
family's constant; §12.2–12.3.) On the procedural point: dependence on unrefereed work is
not unusual in this literature — the asymptotic large sieve [CIS1] is listed by its first
author among unpublished papers and has remained so since 2011, and [CIS2], [CLLR], [So16]
(via [CLLR]), [CL] and [Wu16] all rest on it, [Wu16] citing it by its arXiv number as
unpublished. §2.1 isolates every statement this paper imports from [R], so the dependency
is checkable.

### 1.3 The idea, and the one-sided large sieve

The certificate of [R] consumes, for each χ, a second moment of a zero-counting linear form
— a family pair-correlation-type quantity — **with a single sign: an upper bound suffices.**
Under GRH the relevant object F(α) is known asymptotically for the Dirichlet family: F(α) =
min(|α|, 1) on |α| < 2 — for the PRIMITIVE family this is Chandee–Lee–Liu–Radziwiłł
[CLLR, Thm 2] (their Corollary 3 is the character-sum consequence); the antecedent, over
ALL characters rather than the primitive ones, is Özlük [Öz] — the primitivity is exactly
the defect [CLLR] name and fix. Our Corollary 2 produces a one-sided unconditional
counterpart of that statement's CONTENT for the dyadic family, where the normalisation
matches exactly — for the q ≤ Q family of Theorem 1 the consumed object is the large-sieve
majorant of §12.1, not [CLLR]'s F_Φ.
Unconditionally, in the t-aspect, the trivial bound loses a factor log T ([BGST, Thm 1] —
the unconditional t-aspect Montgomery theorem — and [Go], survey).
The observation this paper runs on — we are not aware of it being recorded, and we state it
as a mechanism rather than a theorem — is that **on family average the loss is a constant**:
the multiplicative large sieve

>   Σ_{q≤Q} Σ*_{χ mod q} |Σ_{n≤N} a_n χ(n)|² ≤ (N + Q² − 1)·Σ|a_n|²

is lossless precisely when the polynomial is shorter than the family, N ≪ Q² — which holds
in the q-aspect for every bandwidth λ < 2 and never in the t-aspect. The price of
one-sidedness is the constant C = Q²/|𝔉_Q| = π⁴/18 + o(1) — the "smaller than Q² by a
constant factor" of [CIS1] — paid on the off-diagonal range 1 < |α| ≤ λ where GRH-truth
would be F = 1. The certificate turns that overshoot into a proportion via a C-penalised
variational problem (§10) whose value exceeds the record's constant for every finite C.

**Remark (the size of the margin).** Corollary 2 exceeds [CIS2] by fifteen points using a
sieve inequality available since 1974, and the margin is large for this literature. The
structural explanation is short: the certificate the sieve plugs into is days old, so the
combination could not have been tried before. We do not ask the reader to weigh
explanations — a result of this profile can be wrong at any margin, and the repository is
built so that checking, not argument, settles it.

### 1.4 Why q ≤ Q, against the field's convention

Every family result in this literature normalises q ≍ Q (we verified fifteen papers without
exception — [Öz], [CIS1], [CIS2], [CLLR], [So16], [So21], [Wu16], [CL], [DPR], [FM], [Mi],
[ILS], [GZ], [BP], [LM]; the restriction lives in the smooth weight's support, and [CIS1]
states the convention verbatim: "the conductor q is restricted by a smooth function Ψ of
compact support in ℝ₊, so q ≍ Q" — [CIS2] and [Wu16] inherit Ψ from it). The convention is
inherited from machinery: the asymptotic large sieve
[CIS1] — the engine under all the asymptotic results — is stated for q ≍ Q. **This paper
does not use the asymptotic large sieve.** (The polylog-T window, by contrast, is NOT part
of the novelty claim: [So21] already operates from the weaker floor (log Q)², below our
(log Q)^{3+ε}. What no prior result covers is the FAMILY — every primitive character of
every modulus q ≤ Q, at equal weight.) The classical sieve's native and sharp range is
exactly q ≤ Q; restricting to the dyadic family costs (the budget N + Q² − 1 is charged
against the largest modulus while the dyadic family supplies only 3/4 of the characters),
which is precisely why C_dyad = 2π⁴/27 > π⁴/18. Our family choice is a design argument, not
a liability. The price is a normalisation discipline (§12): the consumed object is a
large-sieve majorant of a family second moment at the common scale ℒ = log(QT/2π) — it is
NOT the form factor F_Φ of [CLLR], which is correctly normalised only at q ≍ Q; the
conductor-average ⟨log q⟩_{φ*} = log Q − 1/2 + O(1/Q^{…}) is stated with its constant, in the
tradition of Fiorilli–Miller [FM] and Miller [Mi, (2.8) and Rem. 3.2], who compute resp.
recommend exactly this substitution.

### 1.5 Honesty about effectivity

Theorem 1 is an asymptotic statement, as is every prior result in this line; unusually, we
have computed our construction's finite-Q content, and we state it exactly once, in a remark
labelled "for orientation" (§10.4). Nothing effective is claimed anywhere else in the paper,
and no effective figure appears in this introduction or in the abstract by design.

### 1.6 Ancillary computations and disclosure

The numerical constants in this paper (the variational values of §10–11 and the orientation
figures of §10.4) were computed by scripts included as ancillary files, whose identities and
conventions are cross-checked against toy-scale exact computations and a dataset of 5,230
zeros of the 60 primitive L-functions of moduli 5–19; the variational constants are
reproduced to ten digits by three independent implementations (the high-precision
closed-form/Nyström pair is shipped as payoff_hp2.py; the coarse grid anchor qg_check.py
reaches 4–5 digits, and its printed v ≥ 0 flag is a discretisation artifact of its own
edge-sampling scheme: on an exact-cell discretisation min v is POSITIVE at every tested
resolution and decays like the grid spacing — v vanishing linearly at the free boundary —
with the ladder n = 500…8000 shipped in ext_checks.py/log_ext_checks.txt; the closed-form
solve is positive with strict KKT complementarity). The
orientation figures of §10.4 are computed by the shipped table2_check.py. These checks
validate transcription and conventions; they are not part of any proof.

**Disclosure of generative-AI use.** This work was produced substantially by a generative
AI system — Claude (Anthropic) — working under the direction of the named author: the
derivation and optimisation of the construction of §§4–12, the numerical work of §§10–11
and the ancillary scripts, and the drafting of the exposition are all substantially
AI-generated, and this is disclosed in accordance with arXiv's policy on generative AI
language tools. No AI system is listed as an author. Verification of the mathematical
content has been by process rather than by the author's own line-by-line reading: eight
rounds of independent adversarial review (themselves AI-conducted, on blinded bundles,
including one fully independent session with no project history), numerical calibration
gates that reproduce two published third-party constants (§11), reproducible ancillary
scripts for every quoted figure, and the planned Lean formalization of §13, which is the
intended final arbiter. The named author takes full responsibility for the entire contents
of the paper, irrespective of how any part of it was generated. The preprint [R], from
which the six hypotheses H1–H6 of §2.1 are imported, is itself an AI-generated work; §2.1
isolates every statement taken from it, and its public Lean 4 artifact is the intended
means of checking those statements independently of this paper. [At LaTeX conversion this
disclosure moves to an unnumbered statement after the acknowledgements.]

---

## 2. Notation, standing conventions, and imported hypotheses

### 2.1 What is imported from [R], as numbered inputs

Everything this paper uses from the unrefereed preprint [R] is one of the following six
statements, each a sorry-free theorem in [R]'s Lean artifact (file references in the
ancillary notes): **(H1)** the rank–trace certificate: for the direct sum of Hermitian
blocks, #{simple, on-line} ≥ 4 tr Ĝ − ‖Ĝ‖²_F − 2𝒩 − 3N_{II} − (perturbation loss), linear
in (tr, ‖·‖²_F); **(H2)** the Weil-form grid Gram identity Gz = Gp per primitive χ with
density ν_{X,χ} = μ_q + P_{X,χ} and no pole term; **(H3)** the pointwise positivity
g = φ²⋆φ² ≥ 0; **(H4)** the Gevrey-2 profile ϱ₂ with derivative constants (A, B) = (36/e,
2e⁸) and the ramp lemma ‖φ^{(k)}‖₁ ≤ 2Bw(A/w)^k k^{2k}; **(H5)** the ν-generic ends
majorants (the cap-free L-intrinsic forms) with their explicit constants; **(H6)** the
q-uniform local count N_χ(t, t+1] ≤ A₀·log(q(|t|+3)) with absolute A₀. The λ ≤ 1 audit:
[R]'s T-aspect development frequently assumes λ ≤ 1 (equivalently X ≤ T, which fails here at
EVERY λ since X = (QT/2π)^λ ≫ T); the audit of §9 (file-level record in the companion working notes) shows that no PROOF among H1–H6's consumes it — H1 is bandwidth-free
by its signature, H3–H6 are cap-free by construction, and the cap-gated material ([R]'s
prime side, ends instantiation, payoff algebra, calE analytics) is replaced, not imported —
with one Lean-structural caveat stated in §9: [R]'s Params.Valid bundles λ ≤ 1 (and w ≥ 1),
so formalization-track reuse requires the ParamsQ weakening first. This is a LIVE
obligation, named as a risk: the artifact certifies H1–H6 at [R]'s parameter class
(λ ≤ 1, w ≥ 1, D₀ = √T); at this paper's parameters their validity rests on the §9 audit —
"no proof consumes the gated fields" — until the ParamsQ re-parameterisation compiles. The
w ≥ 1 clamp in §10.3's design is a case of the artifact's validity class already shaping
the design. One further artifact
fact, for completeness: none of H1–H6 routes through [R]'s PairCeiling/ modules — the one
place the artifact carries an extra-Lean hypothesis (EnclOK: 256 interval-arithmetic
enclosures computed outside Lean); PairCeiling is imported only by the artifact's top-level
umbrella and its axiom comparator, not by any module H1–H6 lives in (verified by import
grep).

### 2.2 Notation

Primitive χ mod q, q ≤ Q; κ = κ(χ) ∈ {0, 1} the parity. ℒ := log(QT/2π); bandwidth
parameter λ ∈ (0, 2); L := λℒ; X := e^L = (QT/2π)^λ. Window I := [T, 2T], T = (log Q)^{3+ε}.
Grid τ_k := T + kh, h := 2π/L, d := ⌊LT/2π⌋ (≍ ℒ^{1+r+ε} — ℒ^{4.5} at the design of
record: polylog). paperFT h_f(z) :=
∫ f(u)e^{izu}du. Taper: the Gevrey-2 profile ϱ₂ of [R] (their P1–P4: normalized primitive of
e^{−1/x}e^{−1/(1−x)}, derivative constants (A, B) = (36/e, 2e⁸)), φ(u) := ϱ₂((L/2 − |u|)/w),
ramp width w and buffer D₀ — the values w = 3 log ℒ, D₀ = ℒ² log ℒ displayed in §7 are ONE
admissible choice (the REFERENCE design), at which every closing margin is proven with
(log ℒ)²-room; the design of record for the budget is §10.3's jointly-optimised BALANCED
design (w = Θ(1), D₀ balanced),
whose closing condition is met with equality by construction at the corner design — §7.2;
NOT a fortiori from the reference design's: the balanced design trades all its margin for budget. Φ := h_{φ²}; g := φ²⋆φ² ≥ 0;
I′ := (T − D₀, 2T + D₀]. Zone edge **δ′ := 3 log log Q/log Q**. Densities: μ_q(τ) :=
(1/2π)[log(q/π) + Re ψ(¼ + κ/2 + iτ/2)]; P_{X,χ}(τ) := −(1/π) Re Σ_{n≤X} Λ(n)χ(n)n^{−1/2−iτ};
ν_{X,χ} := μ_q + P_{X,χ} — **no pole term** (L(s,χ) entire for primitive χ ≠ χ₀). 𝒩 :=
Σ_χ N_χ(T, 2T) ≍ Q²Tℒ. D_T(v) := ∫_I e^{iτv}dτ. |𝔉_Q| = (18/π⁴)Q²(1 + o(1)); C :=
Q²/|𝔉_Q| → π⁴/18 (dyadic: 2π⁴/27).

All parameter SHAPES are forced (the values within the admissible region are the design's
to choose — §10.3): T ≥ (log Q)³ by the buffer/rate trade (§10.3); D₀ ≫ L by the tail (§7);
w by the Gevrey closing condition; δ′ by the zone-edge calibration (§5); the common taper
(one coefficient vector for the whole family) by the large sieve itself.

---

## 3. The certificate

Per χ, the grid Gram G(χ)_{kl} := ∫ φ̂(τ − τ_k)φ̂(τ − τ_l) ν_{X,χ}(τ)dτ; by the q-uniform
explicit formula (§9) this equals the zero-side matrix Σ_ρ m_ρ h_φ(γ_ρ − τ_k)·conj(h_φ(…));
hat units Ĝ := G/(aL²), a := L⁻¹∫φ². The certificate, displayed in full (H1; the −2𝒩 and
−3N_{II} window terms are part of the inequality, not corrections):

>   Σ_χ N^s_{0,χ}(T, 2T) ≥ 4 tr Ĝ_fam − ‖Ĝ_fam‖²_F − 2𝒩 − 3·N_{II,fam}
>                           − [4B_tr + 2B_F·‖Ĝ_fam‖_F + B_F²],

linear in (tr, ‖·‖²_F), so the direct sum over 𝔉_Q aggregates identically — with exactly
one exception, the tail perturbation, which aggregates by the pair mechanism of §7.3.

**Proposition 3.1 (quantitative pair moment certificate).** Suppose for a given (Q, T):
(i) the display above; (ii) 𝒩(1 − r₁) ≤ tr Ĝ_fam and ‖Ĝ_fam‖²_F ≤ (κ_C + r₂)𝒩, where
κ_C is the Frobenius main-term constant of §§4–6 (NOT the parity κ(χ) of §2.2); (iii)
N_{II,fam} ≤ r₃𝒩; (iv) B_tr ≤ r₄𝒩 and B_F ≤ r₅√𝒩. Then
>   Σ_χ N^s_{0,χ} ≥ [2 − κ_C − (4r₁ + r₂ + 3r₃ + 4r₄ + 2r₅√(κ_C+r₂) + r₅²)]·𝒩,
so the rate of Theorem 1 is carried by the explicit budget of §10.3 (whose rows are exactly
the rᵢ), not lost to an ε-limit. *Proof:* the moment algebra of [R]'s certificate with the
scalar B split into (B_tr, B_F) — trace linearity, the Frobenius triangle inequality, and
√frobSq(Ê) ≤ B_F; three displays. ∎ (The moment inputs are §§4–9.)

The prime side must deliver: the trace ≈ 𝒩 (§9), and frobSq ≤ (κ_C + o(1))·𝒩 with
κ_C = the bracket the variational problem penalises (§§4–6). Everything else is error, and
every error class is named in the ledger (§10.2).

---

## 4. The two-zone analysis in the dual variable

This is the paper's central lemma (full derivations in the companion working notes; conventions anchor-verified to 2×10⁻¹⁴).

**Lemma 4.1 (Parseval).** For real u₁, u₂ ∈ L¹(I): ∬_{I×I} Φ(τ−τ′)² u₁(τ)u₂(τ′) dτdτ′ =
∫_ℝ g(s) F₁(s) conj(F₂(s)) ds, F_j(s) := ∫_I u_j(τ)e^{iτs}dτ — no 2π factor in the paperFT
convention. In particular the PP block 𝓜[P_χ, P_χ] = ∫ g(s)|F_χ(s)|² ds ≥ 0.

**Lemma 4.2 (positivity).** g = φ²⋆φ² ≥ 0 pointwise ([R] gv_nonneg). Consequently every
s-zone 0 ≤ ∫_U g|F|² ≤ ∫_ℝ g|F|² separately: the α-zones of the certificate are zones of
s = αℒ, restrictions of this integral. (What this structure excludes is an ℓ¹-bound on one
half of a split character sum — the (‖a‖₁ + √C‖a‖₂)² failure mode; n-range splits with an
ℓ²/sieve bound on each part are used, and priced, in Lemma 4.4 and §5.) (The pointwise positivity — not the
automatic ψ̂ ≥ 0 — is the load-bearing statement; it is forced by v = φ² ≥ 0, which is
forced by the Gram realization.)

**Lemma 4.3 (family consumption at constant 1).** With F_χ(s) = A_χ(s) + B_χ(s),
A_χ(s) = Σ_n a_n(s)χ(n), B_χ(s) = Σ_n b_n(s)χ̄(n), a_n(s) = −(1/2π)Λ(n)n^{−1/2}D_T(s − log n)
(character-independent weights; no root number, no Gauss sum, no log q), the large sieve
applies POINTWISE in s to each squared half; the χ/χ̄ CROSS term is discharged by
Cauchy–Schwarz in χ pointwise in s — the field's own tool for a χ(nm) bilinear form
([Va, Lemma 6]; [IK, Ch. 7]) — followed by band separation: b_m(s) = conj(a_m(−s))
(exactly: conj D_T(v) = D_T(−v)), and every pair (n, m) is s-separated by log(nm) ≥ log 4
BECAUSE Λ(1) = 0 forces n, m ≥ 2. Split the s-integral at 0; on s ≥ 0 the b-half carries no
T-peak (|D_T(s + log m)| ≤ 2/log m), so **sup_{s≥0}** ‖b(s)‖₂² ≤ (log L + M + O(1))/π²; the
s < 0 half is identical under the mirror s ↦ −s, since g is even and ‖b(s)‖₂ = ‖a(−s)‖₂
(on s < 0 the roles of the two halves exchange). The cross term then costs a multiplicative
(1 + ρ) with **ρ ≤ √(24/π)·√((log L + O(1))/(TL))** (the closed form of the constant
2.7640) — at design scales O(ℒ^{−(1+r)/2}√log ℒ), a δP ≈ 2×10⁻⁵ at Q = 10²⁵ falling to
≈ 3×10⁻⁶ at 10¹⁰⁰ (this is the OUT-zone pricing; the in-zone charge is Lemma 4.5's, ~50×
larger and still sub-minor), well below the smallest budget row. (A factor-√2
bookkeeping in the constant — whether the mirror doubles the half-integral or reproduces
it — is left explicitly conservative: even at √(48/π) the term does not register.) **Remark (a new obligation, and where it comes from).**
CLLR, CIS, Sono and BGST never meet this cross term — each arranges it away (a one-sided
explicit formula whose dual prime sum vanishes; χ×χ̄ bilinear forms by definition; σ ↔ 1−σ
symmetrization). A Gram/Weil realization requires a REAL density, so both halves are present
by construction: the cross term is an obligation specific to this architecture, and it is
discharged by the standard tool for exactly this object. Since g ≥ 0, the budget
X + Q² − 1 = Q²(1 + O(Q^{−δ})) then factors out of the s-integral with no further loss:

>   Σ_χ 𝓜[P_χ, P_χ] ≤ C·|𝔉_Q|·(T/π)·Σ_{n≤X} (Λ(n)²/n)·g(log n)·(1 + o(1)),

i.e. family total ≤ C × family diagonal — and the right side matches the diagonal main term
of [R] (their prop:PP), of which ONLY the constant (T/π) transfers (their error term O(L²X)
and their λ ≤ 1 hypothesis do not; the agreement is a consistency check, not an import). No
Gallagher lemma, no hybrid sieve: the hybrid budget Q²T + N is needed only when N ≫ Q², a
regime this design never enters. Named error terms: the D_T-smearing, which is O(1/T) **in
mass-weighted aggregate, zone by zone at O(log(TL)/(TΔ)) on a zone of width Δ** (the
pointwise version is false near the zone edge and is not used); and the zone-boundary term,
which we state as a lemma:

**Lemma 4.4 (zone boundary).** Restricting the s-integral to U = {|s| ≤ s₀} does not
restrict n. Split a = a′ + a″ at n = Y = e^{s₀} and apply CS-in-χ + Lemma 6.1 (never an ℓ¹
bound) to the a″-part: with R := ∫_U g‖a″‖₂² ≤ (‖g‖_∞/π²)·ℒ·log(Tℒ)·(1+o(1)) (the D_T tails
carry mass 4/y beyond distance y) and P := ∫_U g‖a′‖₂² = (T/2π)∫₀^{s₀} u g(u) du·(1+o(1)),
the relative inflation of the in-zone form is at most 2C√(R/P) = O(√(C log(Tℒ)/(Tℒ))) —
≈ 5×10⁻⁴ at Q = 10¹⁰⁰ at the design of record (the 1.5×10⁻³ formerly quoted was the r = 3
figure); budget row L₇. ∎

**Lemma 4.5 (in-zone cross term).** The χ/χ̄ cross term of the in-zone form — the third
piece of Σ_χ ∫_U g|F_χ|² = Σ|A_χ|² + Σ|B_χ|² + 2Re Σ A_χ conj(B_χ) — is NOT covered by
Lemma 5.2's orthogonality (the correction of §4 withdrew the character-sum route for the
cross term), and bounding it by CS-in-χ + Lemma 6.1 carries the sieve's C even in-zone. Two
bounds, either sufficient. (i) Conservative, and the one the budget charges: the in-zone
cross term inflates the in-zone form by at most C·ρ_U, with ρ_U the cross-term ratio of
this section restricted to U — ≈ 5×10⁻⁴ at Q = 10¹⁰⁰, a δP ≈ 1.6×10⁻⁴ after the in-zone
sensitivity; charged as minor row L₁₂. (The out-zone pricing of the cross term, row ~3×10⁻⁶
at 10¹⁰⁰, does NOT cover this — the in-zone charge is ~50× larger and was previously
unpriced.) (ii) Sharper, and the structural reason the term is harmless: in-zone the a′/a″
split of Lemma 4.4 restricts both halves to n, m ≤ Y = Q^{1−δ′}, where the ℓ¹ × divisor
route that the correction withdrew GLOBALLY is valid LOCALLY — the failing case is
a λ > 1 phenomenon (n reaching X = Q^λ), and in-zone the effective bandwidth is below 1:
the aggregated character sum of the cross pairs is bounded by Q·τ(·) (the divisor average
of Lemma 5.3), giving cross ≪ Q·‖a′‖₁‖b′‖₁·polylog ≍ Q^{2−δ′}·polylog against the family
diagonal ≍ Q²Tℒ — negligible at the same zone-edge calibration as Lemma 5.3, with the a″
tails absorbed by Lemma 4.4's row. ∎ (Corollary 3's in-zone parity projection consumes
Lemmas 4.4–4.5 verbatim: its "in-zone coefficient 1" claim rests on them.)

The in-zone breakpoint is |s| ≤ (1 − δ′)·log Q, not (1 − δ′)ℒ (the difference is budget
row L₁). (Full proofs in the companion working notes, which replace the withdrawn
divisor-bound treatment of the cross term.) ∎

---

## 5. The character-sum arithmetic

(Full derivations in the companion working notes; the aggregated identity anchor-verified to integer precision, 9/9 cases.)

**Lemma 5.1.** For (nm, q) = 1: Σ*_{χ mod q} χ(n)χ̄(m) = Σ_{f|q, f|(n−m)} μ(q/f)φ(f) — the
standard primitive-character orthogonality relation, in the form displayed in [CIS1, §3].

**Lemma 5.2 (aggregated, exact).** For n ≠ m: Σ_{q≤Q}Σ*_χ χ(n)χ̄(m) =
Σ_{f|(n−m), f≤Q, (f,nm)=1} φ(f)·M_{nm}(Q/f), with M_{nm}(y) the coprimality-restricted
Mertens function. (Moduli with (q, nm) > 1 contribute zero since χ vanishes off coprimality.)
Two-line consequence: |Σ_{q≤Q}Σ*_χ χ(n)χ̄(m)| ≤ Q·τ(|n−m|), and ≤ Q·τ(n−1) for the linear
sums. **No PNT is used**; M(y) = o(y) sharpens constants only ([CIS1]'s own Δ′ + Δ″ split is
this optimization).

**Lemma 5.3 (zone bound).** S(Y) := Σ_{n≠m≤Y} Λ(n)Λ(m)(nm)^{−1/2}τ(|n−m|) ≤ C·Y(log Y)³
unconditionally (divisor average — no pointwise τ ≤ Y^ε, whose ε is NOT o(δ′) and would
poison the zone edge). Hence the low–low zone (n, m ≤ Y = Q^{1−δ′}) family off-diagonal is
≪ Q^{−δ′}ℒ·(family diagonal) — and the calibration δ′ = K·log log Q/log Q needs **K ≥ 2**:
at K = 1 the relative error is Θ(1), not o(1), and K = 2 is the proved minimum (the companion working notes); the design adopts K = 3, a conservative choice, not a forced one (the
budget constant s(r + K) prices each unit of K at 0.5073). The payoff cost 0.5073·δ′ sits
inside Theorem 1's rate.

**Lemmas 5.2′/5.3′ (the n + m variants).** Lemma 5.1 applied at the congruence −n ≡ m
(mod f) and aggregated exactly as in Lemma 5.2 gives the identical identity with the
divisor sum over f | (n + m); and since n, m ≥ 2 forces n + m ≥ 4 with no excluded diagonal
(n = −m is impossible for positive n, m), Lemma 5.3's divisor average applies verbatim with
τ(n + m) — the n + m case is, if anything, easier. These are the forms Corollary 3's
in-zone parity projection consumes (§1.1, §12.3).

No case split between the dyadic and full families, and no subfamily sieve, occurs anywhere
in this section — Lemmas 5.2/5.2′ are uniform over q ≤ Q (the parity SUBFAMILY appears only
in Corollary 3's projector, which consumes them at full-family level).

---

## 6. The large sieve

**Lemma 6.1** (Gallagher [Ga67]; Montgomery–Vaughan [MV73], [MV74] — the additive
N − 1 + δ⁻¹ input is [MV74]'s Hilbert inequality, the multiplicative deduction is classical;
see [IK, Thm 7.13] or [Va, Thm 4]; in the exact form quoted, [PT25, (1.1)]).
Σ_{q≤Q}Σ*_{χ mod q}|Σ_{n≤N}a_nχ(n)|² ≤ (N + Q² − 1)Σ_n|a_n|².

This is the paper's only sieve input, consumed in Lemma 4.3 (pointwise in s), Lemma 8.1
(pointwise in τ), and nowhere else. The multiplicative form is not in the Lean tree of [R],
but its additive engine is: the Montgomery–Vaughan generalized Hilbert inequality is proved
in the artifact (Zeta23/MV/), so the from-scratch build is the multiplicative deduction
(duality + Gauss sums + primitive decomposition) on top of an in-tree input — still the
largest single build if this paper's result is to sit beside the record's standard of
evidence (§13).

**Remark (the q/φ(q) weight, declined on the record).** The sharp classical form carries a
weight: Σ_{q≤Q}(q/φ(q))Σ*_χ|·|² ≤ (N + Q²)Σ|a_n|² [Va, Thm 4; IK, Thm 7.13; CIS1, (1.4)].
Dropping q/φ(q) ≥ 1, as Lemma 6.1 does, is valid but not free: under the (q/φ(q))φ*-measure
the family constant would fall to C_w ≈ 4.174 (⟨q/φ(q)⟩_{φ*} = 1.2965, measured at
Q = 2×10⁶), worth ≈ +1.29 points at q ≤ Q and ≈ +1.01 dyadic. We decline it deliberately:
the weighted object is a q/φ(q)-weighted family average — a different statement class from
the sharp-count corollaries whose like-for-like comparisons (+14.99, +8.75) this paper is
built to make — and Lemma 5.2's exact identity, §12.2's conductor averages and the shell
arithmetic are all stated for the unweighted count. The option is recorded; the price of
the cleaner statement is ≈ 1.3 points.

---

## 7. The tail

(Full derivations in the companion working notes; anchors: ramp constants k = 1…8; envelope on the true ϱ₂; the closing ladder; the
60-block pair-vs-scalar toy.)

**7.1 The Gevrey transform bound.** For any GevreyProfile(2, A, B) taper with 2w ≤ L, all
z ∈ ℂ: ‖φ̂(z)‖ ≤ e²·max(2Bw, L)·e^{|Im z|·L/2}·exp(−(2/e)√(w|z|/A)) — from [R]'s proved ramp
lemma ‖φ^{(k)}‖₁ ≤ 2Bw(A/w)^k k^{2k} by k-fold parts and optimization (k = ⌊√r/e⌋; the
integer floor absorbed by the small-r branch). The general limitation is the Ingham/
Paley–Wiener logarithmic-integral criterion [In34, PW, Hö, Thm 1.3.5]: a Q-power of decay at
distance D₀ with bounded ramp fails for EVERY compactly supported taper; for the Gevrey-2
class the constraint is wD₀ ≳ L², met by the reference design with wD₀ = 3ℒ²(log ℒ)² and
by the design of record (§10.3) with wD₀/ℒ² = 8.51 at Q = 10²⁵, 3.80 at 10¹⁰⁰, 2.87 at
10³⁰⁰, decreasing to the asymptote ≈ 2.4 — the (1 + o(1)) is large at headline scales, and
table2_check.py prints the actual values. (A
near-optimal non-quasianalytic class would relax L² to L·polylog; we reuse [R]'s profile and
widen the buffer instead. The record's own hard-coded choice D₀ = √T is not merely
suboptimal here: forcing it makes the ramp row Θ(ℒ^{1−r/2}), which breaks the theorem's rate
class for every r < 4 — so the D₀ generalisation, a field of the ParamsQ re-parameterisation
in the formalization plan, is required, not cosmetic.)

**7.2 The mirrored tail lemma.** Each squared Gram-column entry carries BOTH taper factors at
the same argument ([R]'s own |γ − τ_k|⁻⁴ structure — the two-factor convention is the faithful
mirror), giving suppression exp(−(4/e)√(w·dist/A)) per entry against the off-line
amplification e^{L/2} = X^{1/2} (|Im γ_ρ| ≤ 1/2 against φ̂ of exponential type L/2 — the
unavoidable unconditional price). Row sums by exponential telescoping, window sums by the
q-uniform local count (§9); the closing condition (4/e)√(wD₀/A) ≥ L/2 + log(prefactor) +
log(L/η) holds with Θ(ℒ log ℒ) to spare at the reference design; the design of record
sits ON the constraint by construction — the optimiser drives D₀ to the boundary, so the
margin is 0 nats at every Q (table2_check.py), which is admissible: the condition is an
inequality with the j-sum prefactor and the log(L/η) slack already inside it, and equality
delivers exactly the target θ₀ (the former "+3% at 10¹⁰⁰" pricing referred to the interior
optimum that the w ≥ 1 clamp replaced, and is withdrawn). Sensitivity, computed
(ext_checks.py): REQUIRING a margin of 1 nat moves P_eff by ≤ 0.003 at 10²⁵ and ≤ 10⁻⁴
from 10¹⁰⁰ on; even a forced margin of 0.1·L (≈ 28 nats at 10¹⁰⁰) costs only
0.021/0.002/0.0004 at 10²⁵/10¹⁰⁰/10³⁰⁰ — the zero-margin design is a boundary convenience,
not a knife edge. Either way ‖Ẽ_χ‖ ≤ θ₀ → 0 faster
than any required rate, uniformly over the family. The window j-sum is collected here at crude explicit constants;
the honest collection (per-window decrement 0.350/ℒ at the reference design,
Σ_j = 2.86ℒ·(1 + o(1)); ≈ 0.29/ℒ at the design of record, same order) is in the companion
working notes; the crude form exceeds it by the factor 0.199ℒ —
conservative, never wrong. 

**7.3 Aggregation without a Q-power.** Over the direct sum: the operator norm is a MAX (the
Weyl consumer pays nothing); the trace is additive, |tr Ê_fam| ≤ |𝔉|θ₀/(aL) =: B_tr; the
Frobenius norm is √-additive, ‖Ê_fam‖_F ≤ √|𝔉|·θ₀/(aL) =: B_F. The pair perturbation lemma
(4 tr Ĝ − ‖Ĝ‖²_F − [4B_tr + 2B_F‖Ĝ‖_F + B_F²] ≤ 4 tr Â − ‖Â‖²_F; three lines, generalizing
[R]'s scalar version, which is the case B_tr = B_F) feeds Proposition 3.1's (iv) with
per-block θ₀ → 0 only: **no negative power of Q is demanded anywhere.** (The verbatim scalar
interface would demand θ₀ = o(Q^{−1}·polylog); the pair costs nothing and removes it.)

---

## 8. The ends

(Full derivations in the companion working notes; anchors on the real 60-character family.)

**Lemma 8.1 (family mean square).** B_fam(τ)² := Σ_χ ν_{X,χ}(τ)² ≤ c₁²Q²(ℒ + log⁺(|τ|/4T)
+ C₀)², with explicit absolute (c₁, C₀): the μ-part by digamma envelopes; the P-part by
Lemma 6.1 at fixed τ (the τ-twist leaves ‖a‖₂ invariant — the prime part never touches the
log⁺ growth). **No X-dependence** — the standard large-sieve substitution of a family mean
value for a pointwise bound [Mo71, Ch. 12; IK, Ch. 10; Ga67], applied at the (B, ν)-generic
seam of [R]'s ends layer, where the T-aspect's pointwise B = l + 4√X (the λ ≤ 1 wall) is
replaced by B_fam ≍ Qℒ.

**Lemma 8.2 (bilinear ends bound).** Σ_χ (𝓔₁ + 𝓔₂)(ν_χ) ≤ ∬ W(τ, τ′) B_fam(τ)B_fam(τ′),
with W the ν-free majorant kernels of [R] — quoting the CAP-FREE L-intrinsic constants
(|𝓔₁| ≤ 4(90 + 32cϱ²)·L³B²·(L + l); |𝓔₂| ≤ C₂′·L³B²·(1 + log L)(L + l)), because the capped
forms' hypothesis L ≤ 2l — innocuous in the T-aspect, where L = λl — FAILS here by a factor
10–60. Family relative order: **Θ(ℒ log ℒ/T)** (settling a flagged question from the companion working notes in
favour of the larger figure: the discrepancy was the capped-vs-cap-free shape, not
bookkeeping); budget row 6ℒ log ℒ/T, 0.4% of the budget — nothing downstream moves. The hat-normalization bookkeeping
reproduces the record's λ ≤ 1 ends wall exactly in the T-aspect (consistency anchor).
The bilinear write-out with pinned constants (c₁ = 0.41, C₀ = 6) is in the companion working notes; its
Lean-level instantiation is deferred to the formalization.

---

## 9. The archimedean and explicit-formula layer, uniformly in q

(Full derivations in the companion working notes.) The Weil explicit formula Gz(χ) = Gp(χ) with density ν_{X,χ} and no pole term;
RvM-χ with absolute constant; the local count N_χ(t, t+1] ≤ A₀·log(q(|t|+3)) with ABSOLUTE
A₀. Every analytic constant in the chain is already explicit-and-linear in q in [R]'s own
development (‖L(s,χ)‖ ≤ q‖s‖/Re s; ‖L(w,χ)‖ ≤ 8q(|w.im|+3)); the uniformity is re-packaging,
not new analysis. The buffer count NII_χ ≤ 3A₀D₀·log-scale gives the L₄ budget row.
**The λ ≤ 1 / X ≤ T audit (full file-level table in the companion working notes):** [R]'s Params.Valid bundles lam_le_one, so the reuse obligation is
the ParamsQ weakening (budgeted in the formalization plan), and Tail/ is NOT "reused verbatim"
(theta0_le hard-codes D₀² = T) — both corrected from QR.4's first version; the moment core,
the ends _L variants, the Taper machinery and the uniform local count are cap-free as
audited. The LocalCountChi/MainTermChi re-runs with explicit constants: the companion working notes (the
local count is [R]'s own uniform theorem); the ParamsQ class is formalization-track.

---

## 10. Assembly, the budget, and the proof of Theorem 1

**10.1** Proposition 3.1's proof is the elementary moment algebra of [R]'s certificate with
the pair split (the companion working notes). **Remark (deferred to the formalization):** its Lean-adjacent ε-filter
threading. **10.2 The ledger.** Every explicit-formula error term is classified [P]/[R]/[CS]/
[LS]/[X]/[M] (19 error classes in 17 ledger rows; no term with
τ-mass ≍ T outside the positive form); each row states its normalization (coefficient/tilde/
hat). The aggregation discharge: exactly three components are bounded only in family sum
(the PP off-diagonal, the μP cross, the ends), each consumed by a certified family mechanism
(§§4, 5, 8); no per-block-then-multiply step exists (tabulated in the companion working notes).
**10.3 The budget (the jointly-optimised design of record).** Principal terms: zone breakpoint
s·(log T)/ℒ and zone edge s·δ′ (s = 0.5073, the measured small-δ payoff SECANT at
C = π⁴/18 — the secant at δ ≈ 0.005; the δ → 0 slope is 0.5042, and Table 2 uses the larger
secant at the actual design offset — these two pin the rate class); ramp: the diagonal sandwich (L−2w)³/6 gives relative deficit 1 − (1 − 2w/L)³ ≈ 6w/L in the
frobSq PP-diagonal main (the sandwich is [R]'s, PrimeSideB.lean:384/:386, cited with the
term enumeration in the companion working notes; the
coefficient expansion 1 − (1 − 2w/L)³ = 6w/L − 12(w/L)² + 8(w/L)³ is displayed here, so
6w/L is sharp and charging it in full over-charges for every w > 0) — charged
conservatively in full as r₂ = 6w/L, sign restated (the tapered bracket under-weights the
|α| ≤ 1 mass); the previously quoted 4w/L ([R]'s calE shape) and any "4r₁" reading are
retired; buffer 6D₀/T — priced at the SHARP zero density (N_II at the true T/2π-per-χ
scale), not at the proved absolute A₀ of §9, which is larger: charging the proved constant
would move the finite-Q thresholds, never the asymptotics, and §10.4 states which
convention its table uses; ends (§8) — plus the minor rows (conductor spread O(1/ℒ);
Lemmas 4.4–4.5; the cross term). **The ramp and buffer leave the leading order entirely under the
balanced choice**: minimising the ramp-plus-buffer rows (6w/L + 6D₀/T) subject to the honest Gevrey closing condition
(4/e)√(wD₀/A) ≥ L/2 + log(prefactor) + log(L/η) gives w* = Θ(ℒ^{(3−r)/2}) — Θ(1) at r = 3,
slowly decreasing for r > 3 — CLAMPED at the regime floor w ≥ 1 carried by the reused
lemmas ([R]'s Params.Valid field one_le_w; side conditions in the companion working notes), which binds at
the design of record from Q ≈ 10⁵⁰; so the design sets w = max(1, w*) (= 1 at Q = 10¹⁰⁰;
the previously quoted ≈ 2.3 was the r = 3 optimum under the retired 4w/L objective, and is
withdrawn), D₀ balanced accordingly, every side condition satisfied ([eq:wrange], the QT.a
optimisation, wD₀ ≫ e²A), and ramp + buffer = O(1/ℒ) — the clamp contributes 6/(λℒ), the
same order with constant 6/λ*. One accounting question is stated openly for review: [R]'s
aggregate calE carries a trace-side w/L, which is NOT a separate charge in our accounting —
the hat normalisation a := L⁻¹∫φ² is computed for the TAPERED φ, so the taper's first-order
mass deficit cancels in tr Ĝ against a𝒩, and the trace-side residue is the buffer/tail
rows already charged; if a reviewer rules otherwise, the fully conservative charge is
r₂ + 4r₁ = 10w/L, costing ≤ 0.012 at Q ≥ 10¹⁰⁰ with the rate class and constant untouched.
Total = O(log log Q/log Q), pinned by the two zone terms. **Remark (the constant, computed,
not proved).** The asymptotic budget constant is s·(r + K) at T = (log Q)^{r+ε} and
δ′ = K log log Q/log Q: **≈ 3.3 at the design (r, K) = (3.5, 3)** (2.54 at (r, K) = (3, 2), K = 2 being the
minimum §5 actually proves; ≈ 4.6 at (6, 3), the choice that additionally sits inside
[CIS2]'s window floor (log Q)⁶) — a computed property of the jointly-optimised construction,
NOT in the theorem, whose implied constant is effectively computable. (With the w ≥ 1 clamp
the measured constant at log Q = 2.3×10⁶ is 3.56, converging to s(r + K) from above: the
clamp's 6/λℒ sits a log ℒ below the leading term.) **Distinct objects:**
this ramp row (the taper's diagonal-mass cost, first order in w/L) and §11's ramp cost (the
variational profile's edge perturbation, third order — provably Θ(ρ³)) are different
quantities; both appear because the taper acts once on the window and once at the free
boundary.
**10.4 Table 2 — for orientation only** [the only effective statements in the paper; they
appear in this remark and nowhere else — not in the abstract, not in the introduction]:
Design (r, K) = (3.5, 3) — T = (log Q)^{3.5}, δ′ = 3 log log Q/log Q; zone cost from the
measured P(δ) secant at the design offset; the ramp row charged at §10.3's conservative
r₂ = 6w/L with the w ≥ 1 regime floor enforced; the buffer row at the sharp zero density
(§10.3 — the one row NOT charged at its proved constant); ALL rows carried including the
minor ones, Lemma 4.5's in-zone cross row now among them (computed by the shipped
table2_check.py; log_table2.txt):
**P_eff = +0.20 at 10²⁵, +0.61 at 10¹⁰⁰, +0.68 at 10³⁰⁰, +0.71 at 10¹⁰⁰⁰, → 0.72128;
non-vacuous from ≈ 10²⁰; > 0.5 from ≈ 10⁵¹; exceeding the per-character record from
≈ 10²³⁶.** Under the sharper 4w/L reading every figure improves — by +0.044 at 10²⁵ and by
≤ +0.007 from 10¹⁰⁰ on — and the thresholds fall to ≈ 10¹⁸ / 10⁴⁸ / 10²²⁰. Thresholds below
≈ 10³⁰ are extrapolations of an asymptotic model and should be quoted as such, if at all.
We have not attempted to optimize the implied
constant; an effective version of Theorem 1 would require explicit forms of the inputs of
§§5–9 and is a separate undertaking. (For calibration: no prior FAMILY-AVERAGE result in this line — CIS, CLLR, Sono — states
any effective range; the computed threshold here is believed to be, if anything, smaller
than those constructions would yield.) **10.5** Combining §§4–9 into
Proposition 3.1's hypotheses at λ = λ* proves Theorem 1 at each fixed ε > 0, with the rate
carried by the budget rows rᵢ — no ε-limit is taken, and Q₀(ε) is the largest of §§5–9's
thresholds; the dyadic corollary runs the same argument at C = 2π⁴/27 (its finite-Q
arithmetic from the measured dyadic P(δ) secant at the design offset — 0.64 at the 10¹⁰⁰
offset; the small-δ secant 0.5485 is not usable at finite Q).

## 11. The variational layer

The framework is Sono's [So16] (the RKHS/ODE class; sign-dual: Chirre–Gonçalves–de Laat);
the C-penalised free-boundary solution is a calculation within it: minimize B(v) = ψ(0) +
∫_{|α|≤1}|α|ψ + C∫_{1<|α|≤λ}|α|ψ, ψ = v⋆v, over v ≥ 0, ∫v = 1, supp ⊆ [−λ/2, λ/2];
P = 2 − min B. Optimal v: cos(√2t) in the bulk; edge zones from the quartic (ρ² + 2)² =
(C−1)²(1 − ρ²) with the reflection collapse; free boundary v(b*) = 0. Values (ten digits;
three calibration gates: MT = 0.672500703679412; flat-v 4/3; and [So16]'s PUBLISHED GRH
benchmark reproduced on the strictly better kernel F = min(|α|, 1) at λ = 2 —
0.9322826239 against Sono's 0.93228262, gate3_so16.py — a gate that lands on someone
else's published number. One honest caveat, stated because a referee will notice it: all
three gates are special cases — fixed kernel or a C-limit — so they certify the SOLVER
against outside ground truth; the finite-C free-boundary regime the headline uses is
checked separately, by the closed-form/Nyström agreement to twelve digits and the strict
KKT complementarity of payoff_hp2.py):
π⁴/18 → 0.7212835668 at λ* = 1.2507321515; 2π⁴/27 → 0.7099167448 at λ* = 1.1931581210;
gain positive for EVERY finite C (no break-even); λ = 1⁺ slope ≈ 0.68 (weakly C-dependent:
0.6847 at C = 1, 0.6819 at C = 20). Remark (an
INADMISSIBLE limit, for orientation only): at C = 1 — unattainable, since C ≥ π⁴/18 — the
plateau is exactly 2 − π/(2√2) = 0.8892792655 at λ = π/√2 = 2.2214, which also exceeds the
sieve's own range λ < 2; the constrained λ = 2 ceiling is 0.88836541. Neither is comparable
to [So16]'s GRH 93.22% (M = 0.93228262), which optimises the strictly better kernel
F = min(|α|, 1) — the third calibration gate above. Two documented solver traps (the kernel jump at |α| = 1; the
constrained-window junction) — any new implementation must pass the gates. The taper's ramp does not disturb the optimum, and provably: v* vanishes linearly at the
free boundary, so a ramp of relative width ρ gives ‖δv‖₁ = O(ρ²) and ‖δv‖₂² = O(ρ³); B is
stationary at v* for mass-preserving perturbations, and of the second variation the kernel
part pairs two ℓ¹-masses (O(ρ⁴)) while ψ(0) = ∫v² is an L²-norm — so **ΔP = Θ(ρ³), with the
∫v² term as the reason** (stable across three ramp shapes and three
grids; measured exponent 3.03–3.09). The free-boundary condition itself protects the profile.

## 12. Normalisation of the q ≤ Q family

**12.1** The consumed object is a large-sieve majorant of the family second moment of a
prime-supported Dirichlet polynomial at the common scale ℒ — defined by the construction,
single-scale by fiat of the common taper (itself forced by Lemma 6.1's single coefficient
vector). It is not [CLLR]'s F_Φ. **12.2** The per-character scales that DO survive (μ_q, the
trace, the denominators N_χ) enter as bounded powers of c_q = log q/ℒ under φ*-weights, and
the certificate consumes a ratio of sums, never an average of ratios: ⟨c²⟩/⟨c⟩ and 1/⟨c⟩ are
1 + O(log log Q/log Q) (driven by log T/ℒ; the conductor-spread contribution alone is
O(1/ℒ) — the same order as the rate, and absorbed by it). The conductor averages, CHARACTER-COUNT-weighted (φ*-weighted), by partial
summation against Σ_{q≤x}φ*(q) = (18/π⁴)x² + O(x log x): ⟨log q⟩ = log Q − 1/2 over q ≤ Q,
and log Q + (log 2)/3 − 1/2 = log Q − 0.26895 over (Q/2, Q] — for contrast, the UNWEIGHTED
averages are log Q − 1 and log Q + log 2 − 1 = log Q − 0.30685; the latter is the quantity
Fiorilli–Miller compute, cited here as the precedent for computing (rather than discarding)
the conductor-averaging correction, with [Mi, Rem. 3.2] for the family-wide-scale
recommendation ([Mi, (2.8)] itself is stated for subsets of [N, 2N] and does not cover
q ≤ Q — the two-line derivation above is what does). Every subfamily this paper states pays
its own constant, itemised: Corollary 2 (dyadic: C = 2π⁴/27, conductor average
log Q − 0.26895) and Corollary 3 (even dyadic: C = 4π⁴/27 via sieve positivity, the parity
restriction being exactly a subfamily localisation — §12.3).
**12.3** The Sono comparison is Corollary 3's, and the parity restriction is priced, not
waved: the passage is the two-zone hybrid of Corollary 3 — parity projection IN-ZONE
(½Σ*(1 + χ(−1)), the n ≡ −m orthogonality by Lemma 5.2′, valid there because
n + m ≤ 2Q^{1−δ′} ≪ Q), positivity OUT-ZONE at C → 2C (the projection route fails there at
λ > 1 by Q^{λ−1} — the identical mechanism as the §4 cross term, with no band separation
available since both halves share a band — and no parity-restricted sieve at budget
½(N + Q² − 1) appears in the literature). The mechanism generalises exactly as far as its
TWO obligations do: a subfamily of positive relative density α gets the argument at
out-zone constant C/α — above the λ = 1 floor 0.672500703679, since §11's gain is positive
for every finite C — only if (i) its conductor distribution matches the full family's to
leading order (else §12.2 must be redone, as Corollary 2 does with its own C and its own
⟨log q⟩), AND (ii) it admits an in-zone orthogonality projector expressing its sum through
full-family instances of Lemmas 5.2/5.2′, as the parity classes do. A general density-α
subfamily has no such projector, and positivity in-zone gives coefficient 1/α — the 0.3693
trap named in Corollary 3: selecting within each modulus the half of the characters with
the largest in-zone mass satisfies (i) exactly and fails the floor, so (ii) is not
removable. The proved cases are exactly three: even and odd (the parity projectors), and
dyadic (conductor-restricted, with §12.2 redone).
Residual mismatches with [So21]: the smooth weight (inflates C by 3/(2∫Wx dx) if mirrored —
stated, not paid, since our sharp cutoff is the stronger statement) and nothing else:
[T, 2T] matches [So21] and [R] exactly. ([CIS2]'s window is |γ| ≤ T with floor (log Q)⁶ —
our T = (log Q)^{3+ε} sits BELOW that floor; at r ≥ 6 the theorem enters [CIS2]'s window
class too, at budget constant ≈ 4.6.)

## 13. Formalization outlook

A Lean 4 formalization is planned and scoped: most of §§3–10 mirrors [R]'s tree at
file-level pointers recorded in the companion working notes; the from-scratch builds are
Lemma 6.1 (the sharp multiplicative large sieve) and the re-parameterisation of [R]'s
parameter structure (λ > 1, free D₀). It is not claimed as done; when complete it will be
held to [R]'s own standard (sorry-free, axiom-audited).

## References

[BGST] S. A. C. Baluyot, D. A. Goldston, A. I. Suriajaya, C. L. Turnage-Butterbaugh, *An
unconditional Montgomery theorem for pair correlation of zeros of the Riemann
zeta-function*, Acta Arith. 214 (2024), 357–376; arXiv:2306.04799 ·
[BGST2] iid., *Pair correlation of zeros of the Riemann zeta function I: proportions of
simple zeros and critical zeros*, arXiv:2501.14545 (unpublished; its Thm 2 carries
0.672500703679 conditionally, under a narrow-vertical-box hypothesis — [R]'s contribution
is removing it) ·
[BP] S. Baluyot, K. Pratt, *Dirichlet L-functions of quadratic characters of prime
conductor at the central point*, J. Eur. Math. Soc. 24 (2022), 369–460 ·
[CGdL] A. Chirre, F. Gonçalves, D. de Laat, *Pair correlation estimates for the zeros of
the zeta function via semidefinite programming*, Adv. Math. 361 (2020), art. 106926 ·
[CIS1] J. B. Conrey, H. Iwaniec, K. Soundararajan, *The asymptotic large sieve*,
arXiv:1105.1176 (2011; unpublished — so listed on Conrey's own publication page; its (1.4)
carries the q/φ(q) weight, and its introduction the convention sentence quoted in §1.4) ·
[CIS2] iid., *Critical zeros of Dirichlet L-functions*, J. reine angew. Math. 681 (2013),
175–198; arXiv:1105.1177 ·
[CIS3] iid., *The mean square of the product of a Dirichlet L-function and a Dirichlet
polynomial*, Funct. Approx. Comment. Math. 61 (2019), 147–177; arXiv:1808.02879 (published;
the mean-value engine under [So21]) ·
[CL] V. Chandee, Y. Lee, *n-level density of the low lying zeros of primitive Dirichlet
L-functions*, Adv. Math. 369 (2020), art. 107185; arXiv:1706.02848 ·
[CLLR] V. Chandee, Y. Lee, S.-C. Liu, M. Radziwiłł, *Simple zeros of primitive Dirichlet
L-functions and the asymptotic large sieve*, Q. J. Math. 65 (2014), 63–87; arXiv:1211.6725 ·
[Di24] M. Dickinson, *Zeros of Dirichlet L-functions near the critical line*, Mathematika
(2024), doi:10.1112/mtk.12239; arXiv:2211.06264 (38.2% on the line for the primitive
family of a single fixed modulus, T ≫ q^ε) ·
[DPR] S. Drappeau, K. Pratt, M. Radziwiłł, *One-level density estimates for Dirichlet
L-functions with extended support*, Algebra & Number Theory 17 (2023), 805–829 ·
[FM] D. Fiorilli, S. J. Miller, *Surpassing the Ratios Conjecture in the 1-level density of
Dirichlet L-functions*, Algebra & Number Theory 9 (2015), 13–52; arXiv:1111.3896 ·
[Ga67] P. X. Gallagher, *The large sieve*, Mathematika 14 (1967), 14–20 ·
[Go] D. A. Goldston, *Notes on pair correlation of zeros and prime numbers*, in *Recent
Perspectives in Random Matrix Theory and Number Theory*, LMS Lecture Note Ser. 322,
Cambridge Univ. Press (2005), 79–110; arXiv:math/0412313 ·
[GS26] D. A. Goldston, A. I. Suriajaya, *Zeta Zeros on the Critical Line*, arXiv:2511.20059
(expository, v2 5 Feb 2026; its 41.7%/40.7% are [PRZZ]'s constants, for ζ) ·
[GZ] P. Gao, L. Zhao, *One level density of low-lying zeros of families of L-functions*,
Compositio Math. 147 (2011), 1–18 ·
[Hö] L. Hörmander, *The Analysis of Linear Partial Differential Operators I*, Grundlehren
256, Springer (1983), Thm 1.3.5 ·
[IK] H. Iwaniec, E. Kowalski, *Analytic Number Theory*, AMS Colloq. 53 (2004); Thm 7.13
(the multiplicative large sieve), Ch. 7, Ch. 10 ·
[ILS] H. Iwaniec, W. Luo, P. Sarnak, *Low lying zeros of families of L-functions*, Publ.
Math. IHÉS 91 (2000), 55–131 ·
[In34] A. E. Ingham, *A note on Fourier transforms*, J. London Math. Soc. 9 (1934), 29–32 ·
[LM] J. Levinson, S. J. Miller, *The n-level densities of low-lying zeros of quadratic
Dirichlet L-functions*, Acta Arith. 161 (2013), 145–182; arXiv:1208.0930 ·
[Mi] S. J. Miller, *Density functions for families of Dirichlet characters* (unpublished
note, 1999) ·
[Mo71] H. L. Montgomery, *Topics in Multiplicative Number Theory*, LNM 227 (1971) ·
[MV73] H. L. Montgomery, R. C. Vaughan, *The large sieve*, Mathematika 20 (1973), 119–134 ·
[MV74] iid., *Hilbert's inequality*, J. London Math. Soc. (2) 8 (1974), 73–82 ·
[Öz] A. E. Özlük, *On the q-analogue of the pair correlation conjecture*, J. Number Theory
59 (1996), 319–351 ·
[PRZZ] K. Pratt, N. Robles, A. Zaharescu, D. Zeindler, *More than five-twelfths of the
zeros of ζ are on the critical line*, Res. Math. Sci. 7 (2020), art. 20; arXiv:1802.10521 ·
[PT25] A. Pascadi, J. Thorner, *Large sieves for GL_n and applications*, arXiv:2508.14888
(2025; cited for its (1.1), the sharp multiplicative large sieve) ·
[PW] R. Paley, N. Wiener, *Fourier Transforms in the Complex Domain*, AMS Colloq. 19
(1934), Thm XII ·
[R] "Claude" (Anthropic), *More Than Two Thirds of the Zeros of the Riemann Zeta Function
Are Simple and On the Critical Line*, preprint, dated 11 August 2026 (released 10 August
2026; revised in place 13 August 2026, per the publisher's changelog; not posted to arXiv;
no journal submission, refereeing, DOI or erratum known to us); current PDF:
www-cdn.anthropic.com/95c246936988e43127bc6b2ceb7077c1dad2d68e.pdf (the 10 August version,
titled "…Lie on the Critical Line" and lettered A–E, remains live at
…/564f962e60643842f5fcb4a17c9dbc8f608f1c37.pdf); Lean 4 artifact:
github.com/anthropics/zeta-23-lean (Lean v4.33.0-rc2; Mathlib 51e6992e; Apache-2.0;
sorry-free; #print axioms reports only propext, Classical.choice, Quot.sound; the artifact
retains the A–E lettering, in which the Dirichlet case is Theorem E) ·
[So16] K. Sono, *A Note on Simple Zeros of Primitive Dirichlet L-Functions*, Bull. Aust.
Math. Soc. 93 (2016), 19–30 (the GRH low-lying benchmark, M = 0.93228262) ·
[So21] K. Sono, *Zeros of Dirichlet L-functions on the critical line*, J. Number Theory
(2025; PII S0022314X25000319); arXiv:2105.07422 (2021). The arXiv page carries no journal
reference, which misled earlier drafts of this paper; the even-primitive restriction is in
its (1.8); window floor (log Q)² ·
[Va] R. C. Vaughan, *The Large Sieve* (Penn State lecture notes,
personal.science.psu.edu/rcv4/LargeSieve.pdf), Lemma 6, Thm 4 ·
[Wu16] X. Wu, *Distinct Zeros and Simple Zeros for the Family of Dirichlet L-Functions*,
Q. J. Math. 67 (2016), 757–779; arXiv:1206.1679 (unconditionally 60.261% simple — critical
line NOT required — and 66.433% under GRH, over the Ψ(q/Q)/φ(q)-weighted primitive family;
cites [CIS1] by arXiv number as unpublished) ·
[Wu19] X. Wu, *The twisted mean square and critical zeros of Dirichlet L-functions*,
Math. Z. 293 (2019), 825–865; arXiv:1802.09704 (per character, log q = o(log T): 41.72%
on the line, 40.74% simple-and-on-line — the strongest PUBLISHED per-character
unconditional proportions).
