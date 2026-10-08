# Sparse certificate spaces and a block-order bound for the closed K5

2026-09-28. A refinement of the same category-separation construction, not an
additional supposed innovation or a completed control experiment. The frozen
parent proof remains unchanged. The new upper bound below requires an explicit
block-preserving closure argument; it does not merely assume that taking
closure keeps a displayed open model or its auxiliary solutions.

## 1. Nine small blocks replace the total-degree spaces

Use the notation and exact eight-term joint identity (5) of
`../schur_remainder_recursion/THEORY.md`: x=(m,z,q,r), Delta=mq-z²-m⁴/12,
N=mr-zq-m³z/6, p=m Delta. Let ell be one linear functional on the finite
polynomial space spanned by all products listed below and p, p*x_i, p*L5.
Impose ell(g*b*b') PSD for each row:

| name | g | feature vector b | block order |
|---|---|---|---|
| unit_A | 1 | (N, Delta) | 2 |
| unit_h | 1 | (p, Delta) | 2 |
| unit_domain | 1 | (m, Delta) | 2 |
| delta_B | Delta | (q, m, m², z) | 4 |
| delta_D | Delta | (mz, z, m, m²) | 4 |
| delta_h | Delta | (m³, m², m) | 3 |
| p_C | p | (z, m, 1) | 3 |
| p_h | p | (m², m, 1) | 3 |
| delta_domain | Delta | (1, m) | 2 |

Also impose ell(p)=1, output x_i=ell(p*x_i), and v>=ell(p*L5).
Nonlinear feature polynomials are evaluated by the common linear functional;
they are NOT evaluated at the output. All matrix entries are affine in ell's
scalar coordinates. Using unrelated block variables without the shared
polynomial linear relations would be a different, unproved relaxation.

The implementation expands products into a common monomial coefficient vector
and constructs every entry and output as an exact rational linear expression.
It does not make a nonlinear substitution of output coordinates into an LMI.

### Why the smaller model remains exact

Every squared polynomial in the parent's identity lies in the stated space:
terms 1--8 respectively use unit_A, delta_B, p_C, delta_D, unit_h, p_h,
p_h, delta_h. Their coordinates depend only polynomially on the reference
u=(M,Z,Q,R). `square_coordinates()` gives and checks all eight reconstructions.

The reverse proof must also force its output into the strict domain. It uses
unit_domain on (m,Delta) to get ell(Delta²)>0 from ell(p)=1; delta_domain on
(1,m) to get M=ell(pm)>0; delta_B on (m,z) for the perspective inequality;
delta_h on (m²,m) and p_h on (m,1) for the cubic remainder inequality.
Those are exactly the ingredients of parent identity (8), giving
Delta(output)>0 without assuming a representing measure.

Applying the parent identity at the reference equal to the output then proves
ell(p*L5)>=L5(output). Each term is nonnegative by one of the nine PSD blocks
and the positive reference generators. Conversely ell(f)=f(x)/p(x) supplies
rank-one PSD witnesses for every x with m,Delta>0. Therefore the projection is
exactly the same open E5, not just an outer approximation.

The block orders are now at most four. No block-order minimality is claimed.
The open representation still does not accept zero denominators.

## 2. Closure preserves a block-order ceiling: full proof

This is the blockwise version of the standard affine-polar argument, not a
new general closure theorem. See Gouveia--Netzer, Section 3.1, Proposition 3.1
and Corollary 3.4: https://arxiv.org/html/0911.2750 . The regularity repair is
the same product-face compression already written in the project's frozen
normal-cone lower-bound proof.

**Lemma.** If nonempty convex C has a finite affine lift over real PSD blocks
of order at most r>=1, then so does its closure.

Write C={Lw+c : Aw=b, w in Q}, where Q is a finite product of such blocks.
Free scalar variables can be differences of nonnegative scalars. For each
block take the span of the ranges of all feasible block matrices. Finitely
many feasible tuples span these subspaces; their average is positive definite
on every retained subspace. Compress each block to that subspace and discard
zero blocks. The same C now has a strictly feasible lift with block orders
no larger than r. Equalities and the affine output map are retained.

Define the cone of affine functions nonnegative on C:

    A(C) = {(a,z): a+z'x>=0 for every x in C}.

For (a,z) in this cone, the primal infimum of z'Lw over Aw=b,w in Q is finite
below. Strict primal feasibility gives an attained dual optimum, even when
the primal infimum is not attained. Consequently

    A(C) = {(a,z): exists eta,
              L'*z-A'*eta in Q,
              a+z'*c+b'*eta >= 0}.                 (1)

Conversely (1) implies nonnegativity by weak duality. Thus A(C) has a lift
with the same block ceiling and one extra scalar inequality. This works for
unbounded and nonclosed C; no interchange of closure and projection is used.

Apply the same argument, including face compression if needed, to the
nonempty convex set A(C). It proves that A(A(C)) has the same block ceiling.
Separation of a point from a closed convex set gives

    closure(C) = {x : (0,1,x) belongs to A(A(C))}.   (2)

Here the first coordinate 0 is the constant of an affine functional on
(a,z), and (1,x) are its linear coefficients. Formula (2) asks exactly that
a+z'x>=0 for all affine functions nonnegative on C. It is an affine slice
and does not enlarge the PSD blocks. This proves the lemma.

The statement controls the maximum block order, not the number of blocks,
the number of auxiliary scalars, their conditioning or an efficient algorithm
for finding all required face-compression bases.

## 3. Consequence for the same K5 core

Use the parent's interior K4 homothety, which adds only order-two blocks to
the existing maximum-order-three K4 lift. Intersect it with two copies of the
present nine-block open lift for the lower and complementary upper bound.
This gives exactly the K5 band over int(K4), with maximum block order four.
Its closure is full K5 by mixing densities with 1/2. Applying the lemma proves

    3 <= sxdeg(K5) <= 4.

The lower bound is the previous general theorem. Thus finite SDP versus no
finite SOCP now has a controlled block-order upper bound, including every
singular boundary point. It remains unknown whether three or four is minimal.

This does NOT claim that simply substituting a singular point into the nine
open blocks works, or that the final closed model has only those nine blocks.
Its actual twice-polar coefficient implementation is not constructed here.
No optimizer or new closed-loop simulation was run for this result.

The earlier 126-order open model and its frozen report remain valid but are
superseded for the block ceiling by this explicit sparse-space proof plus
the blockwise closure lemma. The new conclusion is not inferred from test
counts or numerical success. Priority and paper significance remain under
audit; the generic affine-polar operation and sparse-basis selection are not
claimed as separately original methods.
