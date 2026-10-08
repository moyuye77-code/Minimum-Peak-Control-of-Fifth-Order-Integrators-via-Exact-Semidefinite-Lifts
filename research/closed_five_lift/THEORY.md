# A direct finite affine lift of the complete K5, including singular inputs

2026-09-28. This implements the existing category-separation core; it is not
a second innovation or a completed journal submission. The earlier frozen
existence and sparse-space proofs remain valid. The construction here avoids
an unimplemented twice-polar closure by proving a particular smaller pencil
already has the complete closed projection.

## 1. Explicit prescription

Let y=(m,z,q,r,v), x=(m,z,q,r), bar y=(1,1/2,1/3,1/4,1/5)-y. Use the
same Delta, N, p=m Delta, and polynomial p*L5 as the frozen parent proofs.
Require x in the existing seven-block affine K4 lift, with its six auxiliaries.

For each of x and bar x introduce a separate linear functional ell on the
span of the following matrix entries and the six output polynomials
(p,pm,pz,pq,pr,pL5):

| weight | feature vector |
|---|---|
| 1 | (N,Delta) |
| 1 | (p,Delta) |
| Delta | (q,m,m²,z) |
| Delta | (mz,z,m,m²) |
| Delta | (m³,m²,m) |
| p | (z,m,1) |
| p | (m²,m,1) |

For every row (g,b), impose ell(g*b*b') PSD. Also require

    ell(p)=1,  ell(p*x_i)=x_i,  ell(p*L5)<=v,

with barred outputs and upper value 1/5-v for the second copy. All entries
are affine in functional coordinates, not nonlinear evaluations at x.

The exact coefficient matrix of these polynomials has row rank 19. Its first
six selected basis elements are precisely (p,pm,pz,pq,pr,pL5). RREF provides
an exact rational basis and expresses all matrix entries in it. This is only
linear elimination, not an approximation or a numerical rank threshold.
Eliminating the first five functional coordinates leaves 14 auxiliaries per
side. Together with K4 this gives 6+14+14=34 auxiliary scalars,
five output coordinates, 21 PSD blocks and two scalar inequalities.
The PSD orders are K4's (3,3,3,2,2,2,2), followed by two copies of
(2,2,4,4,3,3,3). The maximum order is four.

**Theorem.** The projection of this finite affine model is exactly K5.
No strictly positive denominator constraint or tolerance is part of the model.
The following two proofs, not numerical tests, establish both inclusions.

## 2. Sufficiency, including every zero denominator

First suppose x belongs to int(K4). Both p(x) and p(bar x) are positive.
Every square in the parent's full joint identity still belongs to one of the
seven retained spaces. Take its reference equal to the output x and apply
ell. Normalization and the output equations cancel the first-order term;
all reference weights are positive. Thus ell(pL5)>=L5(x). The barred copy
gives the complementary upper bound. The classical interior lower/upper
description therefore proves y in K5.

The deleted blocks, unit_domain and delta_domain, were used only to infer
strict positivity of the output denominator without any prefix constraint.
They are not needed in the preceding argument, because int(K4) supplies it.

Now let any tuple in the proposed model have possibly boundary prefix.
The density rho=1/2 gives a feasible central tuple by point evaluation and
the existing atomic K4 witness. Mix the ENTIRE feasible lifted tuple with
this central tuple, with weight epsilon>0 on the latter. Every pencil and
scalar inequality remains feasible by affinity and convexity. Its prefix is
in int(K4), since the center is interior. The preceding argument puts its
output y_epsilon in K5. Compactness of K5 and y_epsilon->y prove y in K5.

This reasoning does not assume that a projected spectrahedron is closed.
It establishes that every point in this particular projection lies in the
known closed target. Nor does it replace the necessity proof below.

## 3. Necessity: finite witnesses, not infinite auxiliary limits

For a physical point y, first take rho_epsilon=(1-epsilon)rho+epsilon/2 and
its moments y_epsilon. For epsilon>0 the point-evaluation functional

    ell_epsilon(f)=f(x_epsilon)/p(x_epsilon)

satisfies all seven PSD blocks, output equalities and the lower inequality;
the barred copy has the same property. The numerator is always a polynomial.
We now prove that the 19 SELECTED functional coordinates have finite limits.
The raw monomial evaluations need not have finite limits and are not used.

If m>0 and Delta>0, nothing is singular. Suppose m>0 and Delta=0.
The physical unshifted exponential Hankel condition, after the established
invertible congruence and the first Schur elimination, has leading entry
Delta/m and off-diagonal N/m. Its positive semidefiniteness implies N=0.
This is a necessary classical moment condition, not the new lift theorem.

Put kappa=q-z²/m-m³/12=Delta/m. It is concave for m>0, since z²/m and m³
are convex there. The central kappa is strictly positive. Thus

    kappa(x_epsilon)>=epsilon*kappa(center),
    Delta(x_epsilon)>=c*epsilon for sufficiently small epsilon,
    N(x_epsilon)=O(epsilon), m_epsilon bounded away from zero.

All entries of the first two normalized matrices are combinations of
N²/p, N*Delta/p=N/m, Delta²/p=Delta/m, p²/p=p and p*Delta/p=Delta.
They have finite limits; N²/p tends to zero. The three Delta-weighted
matrices reduce to b*b'/m, hence have finite limits. The two p-weighted
matrices reduce to b*b'. Finally pL5/p=H2+N²/p has a finite limit, and
the other outputs are 1 and x_i. These are all selected basis candidates.

If m=0, nonnegative density implies x=0. On the chosen central path
m,z,q,r are O(epsilon), Delta and N are O(epsilon²), and p is of order
epsilon³ with positive leading coefficient. Thus the first two normalized
blocks are bounded (indeed vanish); every feature in a Delta-weighted block
is O(epsilon), so b*b'/m is bounded; the p-weighted blocks and L5 are bounded
as well. The output moment v is zero. This covers the only zero-mass prefix.
The same arguments apply to the complement, including full density rho=1.

Every selected basis element is an actual retained matrix entry or an output
polynomial. Its normalized evaluation is a bounded rational function of
epsilon, so its singularity is removable. Define its finite functional value
by that limit. Linear polynomial relations commute with limits, so the 19
values define ONE consistent functional on the whole required space. Each
resulting matrix is a finite PSD limit; all affine equations persist. Taking
limits in ell_epsilon(pL5)<=v_epsilon proves the scalar inequality too.
The existing K4 atomic witness supplies finite prefix auxiliaries.

This produces an actual finite feasible tuple for y, not merely a sequence
with divergent auxiliaries. `witness()` computes these removable limits in
exact rational arithmetic for rational physical moments. It is not a generic
membership oracle for arbitrary supplied prefixes.

## 4. What changes and what does not

The model directly realizes the bound 3<=sxdeg(K5)<=4 already proved in the
previous turn. The lower bound is unchanged. Unlike the earlier nine-block
open model, this one includes singular points with finite coordinates, and
the 21-block/34-auxiliary count now refers to FULL K5. Minimality, all-order
SDP existence and a universal numerical advantage are not claimed.

Removing the two blocks without the K4 prefix, without retaining all common
polynomial relations, or without the finite-limit necessity proof would not
establish this theorem. It is not a general rule that dropping domain blocks
automatically closes a rational lift.

The construction continues to use classical Markov moments, domain-weighted
first-order certificates and sparse polynomial spaces. Their general methods
are prior art. Direct K5 priority and scientific importance remain under audit.
The present support-value checks are static numerical diagnostics, not new
closed-loop simulations, real-time guarantees or a completed paper.
