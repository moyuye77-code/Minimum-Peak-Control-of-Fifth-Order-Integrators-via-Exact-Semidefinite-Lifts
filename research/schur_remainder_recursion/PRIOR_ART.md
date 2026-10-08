# Attribution and next priority question

2026-09-28. This is a bounded audit, not exhaustive priority clearance.

1. The exponential moment transform and the Schur lower boundary are classical
   Markov-moment machinery; see the prior frozen general_congruence_gate
   source audit. Block elimination, completing a square and translating
   polynomial coefficients are classical algebra. The all-order lemma here
   gives a specific restricted certificate reduction; it does not rename
   those tools as a new optimization framework.

2. J. Nie, First Order Conditions for Semidefinite Representations of Convex
   Sets Defined by Rational or Singular Polynomials, arXiv:0806.4721v1:
   https://arxiv.org/html/0806.4721 . Reread definitions in Section 2.1,
   Hessian-certificate discussion in Section 2.3.1, rational definitions
   (3.7)--(3.8), the adjacent denominator/domain statements and Section 3.3.1
   on epigraphs. Separate positive multipliers and domain-weighted SOS are
   prior art. An arbitrary multiplier need not stay separate after a
   reference-dependent translation: our determinant restriction is explicit.
   Our lift projection and strict-domain enforcement are proved directly;
   a theorem invocation is not used to insert missing zero denominators.

3. J. Gouveia and T. Netzer, Positive Polynomials and Projections of
   Spectrahedra, arXiv:0911.2750v2:
   https://arxiv.org/html/0911.2750 . Reread Section 3.1, Proposition 3.1,
   Lemma 3.2/proof, Proposition 3.3/proof and Corollary 3.4/proof. The closure
   of a shadow is a shadow, using affine polarity twice after handling strict
   feasibility/affine hull. This supplies the final K5 existence step, not a
   compact block-size or solver-performance guarantee. The HTML generated
   title-page date differs from its arXiv version stamp; no new publication
   date is inferred from it.

4. Direct searches included integrator/fifth-order/finite semidefinite lift,
   bounded-density moment SDP, Markov Schur translation, zonoid shadows and
   rational q-module/perspective constructions. They did not identify an
   explicitly matched finite affine K5/non-SOCP theorem or the same remainder
   certificate. Many hits are asymptotic moment hierarchies, nonlinear moment
   tests, or unrelated uses of fifth relaxation order. Those are not evidence
   for or against priority of this exact claim. No all-order novelty claim is
   made; the closest specific prior theorem still needs a focused audit.

The strongest candidate consequence is now a control reach-set category
separation, finite SDP versus no finite SOCP for K5, not simply a smaller
fourth-order implementation. The general lower bound is previous work in
this project based on existing extension-degree theory. Neither the separate
ingredients nor the current search justify claiming a completed TAC paper.
An all-order certificate remains unsolved, while the user's publication goal
does not automatically require settling every order.
