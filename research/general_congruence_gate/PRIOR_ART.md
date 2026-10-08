# Source scope for the polynomial-congruence gate

2026-09-28. No claim of exhaustive search or cleared priority.

## Exponential transform and reconstruction

L. Gosse and O. Runborg, *Existence, uniqueness and a constructive solution
algorithm for a class of finite Markov moment problems*, arXiv:0809.3714.
Source: https://arxiv.org/pdf/0809.3714 . A local unchanged copy is at
../paper/sources/Gosse_Runborg_0809.3714.pdf and its hash is in the archive.

Read PDF pages 9-10 (Lemmas 1-3, exponential recurrence), 18-22 (Theorem 3,
Remark 7, residue Lemma 4, Proposition 4, Theorem 4 and its proof; Corollary 1
only through the part on page 22). Pages 9 and 18 visually inspected. The
PDF has 25 pages; this is not a full-paper proof audit. Its first-page arXiv
stamp and generated title-page date differ; neither is treated as a verified
journal publication date. The text's m_k includes branch power sums in the
later algebraic problem. Our y_j are integral moments, so that m_k=k*y_(k-1);
our transform additionally uses the negative exponent and 1-exp(-Y).

The source already supplies the exponential/Toeplitz machinery and residue
interpretation. THEOREM 3 additionally assumes a singular extension; it must
not be quoted as uniqueness of every density with an arbitrary even prefix.
The all-order interior boundary derivation in THEORY.md is an explicit
classical support/Schur argument, not a newly claimed Markov theorem.

## What is required to turn a PMI into an affine lift

J. Nie, *Polynomial Matrix Inequality and Semidefinite Representation*,
arXiv:0908.0364v2 HTML: https://arxiv.org/html/0908.0364 . Re-read Section 2
definition of matrix concavity and Theorem 2.2/proof; Section 4 definition
(4.2), Lemma 4.1, Theorem 4.2/proof, Corollary 4.3 and its scope. Matrix
concavity is a hypothesis of these sufficient conditions, not a consequence
of convexity of the set defined by the PMI. Our obstruction prevents applying
this hypothesis to the specified polynomial-congruence family. It says
nothing negative about Nie's theorem or every other possible representation.
Closedness, domain/interiority and denominator conditions cannot be omitted.

C. Scheiderer, *Spectrahedral shadows*, arXiv:1612.07048 HTML:
https://arxiv.org/html/1612.07048 . Read Theorem 3.4/proof, Remark 3.5,
Remark 3.7 and Corollary 3.8/proof, and the adjacent necessary-factorization
argument. A fixed finite-dimensional square space is essential. These
criteria do not mean that pointwise rational convexity, or separately
obtaining some positivity certificate for each support, proves an affine
lift. No application of the obstruction examples in Section 4 to K_d has
been established here.

## Search and interpretation

Queries included polynomial L1 norms and spectrahedra, bounded-density
semidefinite representations, integrator reach-set shadows, Markov moment
matrix concavity, and zonoid shadows. They located the known direct-control
classification and moment/lifting references, not a matched published
polynomial-congruence no-go theorem or a usable all-order affine lift.
This limited negative search is NOT priority clearance. Search-result crawl
dates were not treated as publication dates. Unread theses and unrelated
results were not used to support a claim.

The new result's role is to close a proof-template search class. It does not
decide the general SDP or K4 SOC question, and it is not promoted to an
independent main contribution. The next mathematical test is the uniform
first-order certificate for the rational Schur boundary, including prefix
domains and singular closure; standard rational lifting itself is prior art.

Read-only PDF inspection initially hit a console encoding error; it was
repeated with UTF-8. Web screenshots failed, so the original PDF was rendered
locally and visually checked. No source PDF was edited or re-exported.
