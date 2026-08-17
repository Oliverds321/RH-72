# Provenance

**Author: Oliver D'Souza.** The mathematics, numerics, and exposition of this work were
produced substantially by a generative AI system — Claude (Anthropic) — working under the
author's direction, across a set of isolated AI sessions coordinated by the author. No AI
system is listed as an author. The paper's §1.6 carries the formal disclosure in arXiv's
vocabulary; the author takes full responsibility for the entire contents.

**How it was checked.** Verification was by process rather than by any single reader:

- **Blind adversarial review, eight rounds.** Each round, the complete working bundle was
  delivered to reviewer sessions that shared no context or history with the authoring
  session — the isolation was enforced physically, by the author carrying archives between
  sessions by hand. Reviewers were tasked to attack, not to approve; two additional
  reviews were performed by a fully independent session given a deliberately neutral
  bundle with no suggested reading order. Every review, verdict, and point-by-point
  response is preserved verbatim in the project record.
- **Numerical anchors and calibration gates** for every quoted constant (see
  `VERIFICATION.md`), including gates that reproduce two independently published
  third-party constants.
- **A Lean 4 formalization**, in progress, held to the same standard as the foundation
  work's artifact (sorry-free, axiom-audited). This is the intended final arbiter.

**What the process caught.** Across the rounds, review found and fixed two proof-level
defects and a larger number of presentation and citation errors — all before any public
release. The construction itself has not been successfully attacked in any round.

**The foundation.** Six statements (paper §2.1, H1–H6) are imported from [R], an
AI-generated Anthropic preprint that is at this writing unposted and unrefereed but
machine-verified in its public Lean artifact. The paper isolates the imports so the
dependency can be judged cleanly.

**The full record.** The detailed working notes (complete derivations behind each paper
section) and the complete review archive are packaged as companion documentation —
available to reviewers and collaborators on request, and intended for addition to this
repository. Nothing in the companion record contradicts anything here; it is the same
project at working depth.

**Timeline.** The originating prompt, the construction, the paper, and the eight review
rounds all occurred within a single day (2026-08-17). This is stated as fact, not as
methodology advice.
