# Why polynomial triangular convexification does not extend the K4 proof

2026-09-28. This is a proof-template obstruction, not nonexistence of an SDP
lift. It applies to an entire class of coordinate-dependent polynomial
congruences, not only weighted-homogeneous ones. The rational alternative
at the end is classical moment geometry and an unresolved construction
route, not a claimed new all-order SDP theorem.

## 1. Classical coefficients and the class under test

For y=(y0,...,y_(d-1)), define

    Y(w)=sum_j y_j w^(j+1),
    b_j=[w^(j+1)](1-exp(-Y(w))).

These exponential moment coefficients are classical; the sign/scaling
translation to Gosse--Runborg is recorded in PRIOR_ART.md. If e0=1 and
exp(-Y)=sum e_k w^k, then

    n e_n = -sum_(k=1)^n k y_(k-1) e_(n-k),  b_j=-e_(j+1).

For n>=3 let H_n(y)=[b_(i+j)]_(i,j=0,...,n-1), using at least 2n-1 original
moments. Consider any lower-triangular polynomial matrix T(y) with constant
nonzero diagonal. Entries may depend on all the original moments and have
arbitrary finite polynomial degree. The proposed proof strategy would make
G(y)=T(y)H_n(y)T(y)' matrix-concave and then use a first-order matrix-lifting
certificate. Constant congruences preserve matrix concavity, so diagonal
normalization reduces to unit lower triangular T without loss.

**Proposition.** No such T makes G matrix-concave on a nonempty open set
with y0>0. In particular it cannot be matrix-concave throughout the interior
of the bounded-density moment body K_(2n-1).

The proposition concerns these unshifted Hankel blocks and polynomial
triangular congruences. It does NOT cover arbitrary rational or nontriangular
changes, different matrix inequalities, or all affine lifts of K_d. In
particular it is not an SDP nonrepresentability theorem.

## 2. A zero-curvature row forces affine off-diagonal entries

For a twice-differentiable symmetric matrix function G, matrix concavity
on an open convex neighborhood requires -D_h^2 G(y) PSD for every direction h.
A PSD matrix with a zero diagonal entry has its entire corresponding row
zero. Consequently, if G_00 is affine, every G_0j must have zero Hessian.
For polynomial entries this is a polynomial identity, so each G_0j is
affine as a polynomial, not merely at a sampled point. The same argument
applies on an open subset; neighborhoods are enough for the contradiction.

Write the first five moments as (m,z,q,r,v). The coefficients are

    b0=m,
    b1=z-m^2/2,
    b2=q-m*z+m^3/6,
    b3=r-m*q-z^2/2+m^2*z/2-m^4/24,
    b4=v-m*r-q*z+m^2*q/2+m*z^2/2-m^3*z/6+m^5/120.

Only the first three rows/columns matter. Let T10=P, T20=Q, T21=R be
arbitrary polynomials in all available moments. Since G00=m,

    G01 = z-m^2/2+m*P

must be affine. Divisibility by m gives P=m/2+c for a real constant c.
Subtract c times row zero by a CONSTANT congruence; this preserves matrix
concavity, polynomial triangularity and unit diagonal. Hence assume
P=m/2 and G01=z.

Next G02=b2+R*b1+Q*m is affine. At m=0 it equals q+z*R(0,z,q,r,v,...).
Polynomial divisibility by z implies R|_(m=0) is a constant. Subtract that
constant times row one, and then a constant times row zero to remove a
possible linear m term. These are again constant congruences. We obtain
G02=q, and R=m*F for some arbitrary polynomial F. The equation G02=q
then forces

    Q=z-m^2/6-F*(z-m^2/2).

Thus every candidate has, up to constant congruences, the normalized block

    T = [[1,0,0],
         [m/2,1,0],
         [z-m^2/6-F*(z-m^2/2),m*F,1]].                 (1)

This derivation does NOT assume weighted homogeneity or bound the degree
of F. Restricting a search to a single constant parameter is unnecessary.

## 3. The remaining polynomial freedom is forced away

Direct multiplication gives

    G11=q-m^3/12,
    Delta=m*q-z^2-m^4/12,
    G12=r-m^2*z/6+(F-1/2)*Delta.                       (2)

For every direction with h_m=0, D_h^2 G11=0. The zero-row condition therefore
forces G12 to be affine in all the remaining moment variables, with m held
fixed. This identity holds polynomially. Regard the polynomials as elements
of R[m][z,q,r,v,...]. Delta has total degree two in the bracketed variables.
If F-1/2 were nonzero, its product with Delta would have degree at least two
(degrees of nonzero products add over the integral domain R[m]). The other
terms in (2) are affine in those variables and cannot cancel that highest
degree. Therefore F=1/2 identically.

Every possible polynomial candidate is consequently forced to

    G* = [[m, z, q],
          [z, q-m^3/12, r-m^2*z/6],
          [q, r-m^2*z/6, v-m^2*q/12-m*z^2/4+m^5/720]]. (3)

The Hessian of its last diagonal entry in (m,q) is

    [[m^3/36-q/6, -m/6],
     [-m/6,       0]],

whose determinant is -m^2/36<0 when m>0. It is indefinite, contradicting
even scalar concavity of G*22. This proves the proposition for n=3.
For n>3, the leading three-by-three block of T H_n T' depends only on the
leading three-by-three block of T and the first five moments. Extra moments
may occur in P,Q,R, but were already allowed in the polynomial argument.
A principal submatrix of a matrix-concave function must be matrix-concave,
so the same contradiction applies to every n>=3.

K_d really has open feasible neighborhoods: at rho=1/2 the moment map of
polynomial perturbations h of degree <=d-1 is the invertible Hilbert Gram
matrix. Every sufficiently small moment perturbation has bounded h and
keeps 0<1/2+h<1. Thus the open-set reasoning is legitimate in the original
physical moment body, not just in infeasible ambient directions.

## 4. An explicit physical Jensen violation

Use the densities

    rho_plus(t)=3/4-3t/8,
    rho_minus(t)=1/4+3t/8,  0<=t<=1.

Both lie strictly between zero and one. Their first five moments are

    y_plus =(9/16,1/4,5/32,9/80,7/80),
    y_minus=(7/16,1/4,17/96,11/80,9/80),
    y_center=(1/2,1/4,1/6,1/8,1/10)=(y_plus+y_minus)/2.

For (3),

    [(G*(y_plus)+G*(y_minus))/2-G*(y_center)]_22
       =43/6291456 >0.

The opposite inequality is necessary for matrix concavity. Feasibility is
certified by the actual densities, not by the transformed PMI under test.
The same densities provide moment vectors of every higher order. No numerical
optimizer or floating eigenvalue is used in this obstruction.

## 5. What remains available: classical rational interior boundaries

The failed polynomial matrix template must not be conflated with failure
of rational scalar boundary methods. Here is a self-contained derivation
of the classical all-order Schur route, to fix the next research question.

Let x be an interior point of K_(d-1), and let L_d(x) be the minimum last
moment integral t^(d-1)*rho with this prefix. Weak-star compactness attains
the minimum. The value is finite and convex on its prefix domain, so it
has a subgradient lambda at every interior x. Optimality is then pointwise
minimization of the monic polynomial

    p(t)=t^(d-1)-sum_(j=0)^(d-2) lambda_j*t^j.

The minimizing density is uniquely rho=1_{p<0}, up to null sets. If it had
k<d-1 genuine interior switches, a degree-k polynomial with those switches
and the appropriate sign would expose its prefix in K_(d-1), contradicting
interiority. Thus p has exactly d-1 distinct roots in (0,1).

If d=2n+1, rho is the indicator of n disjoint interior intervals. If
d=2n+2, it is the indicator of n+1 intervals, the first starting at zero
and all the right endpoints strictly below one. For intervals [a_i,b_i],
direct integration gives

    1-exp(-Y(w))=1-prod_i (1-b_i*w)/(1-a_i*w)
                =w*sum_i c_i/(1-a_i*w),
    c_i=-prod_j(a_i-b_j)/prod_(j!=i)(a_i-a_j)>0.        (4)

Positivity follows from interlacing: each factor ratio with j!=i is positive
and -(a_i-b_i)>0. Hence b_j=sum_i c_i*a_i^j. This is the classical residue/
Vandermonde representation of the exponential transform, not a new theorem.

Define s=0 for d=2n+1 and s=1 for d=2n+2. The matrix

    B=[b_(s+i+j)]_(i,j=0,...,n) = [[A,b],[b',b_(d-1)]]

has rank n at the minimizing density. Its leading A is positive definite:
there are exactly n distinct effective atoms with positive weights; in the
even case multiplication by the left endpoint removes the atom at zero.
All entries except the last are determined by x. Since
b_(d-1)=y_(d-1)+P_d(x), the Schur equality yields

    L_d(x)=b' A(x)^(-1) b-P_d(x).                       (5)

Thus the interior lower boundary is a rational convex function for every
order, and the upper boundary is

    U_d(x)=1/d-L_d((1,1/2,...,1/(d-1))-x).

Every last moment between these extrema is realized by a convex mixture.
The full K_d is the closure of this interior-prefix band: mixing any feasible
density with 1/2 produces interior prefixes and converges back to it.
This closure statement does not permit substituting singular A into (5).

At d=5 an equivalent expression from (3) is

    A=[[m,z],[z,q-m^3/12]], B=[q,r-m^2*z/6]',
    L5=m^2*q/12+m*z^2/4-m^5/720+B' A^(-1) B.

Standard rational Schur congruences can diagonalize (3), so the polynomial
obstruction does not prohibit that alternative. HOWEVER, neither rational
convexity nor a nonlinear PSD condition supplies a finite affine SDP lift.
To use Nie's q-module route one still needs a uniform, domain-weighted
sum-of-squares identity for the denominator-cleared first-order remainder
and correct treatment of the prefix domain and its singular boundary.
No such all-order certificate, induction, or all-order SDP exists in this
module. An arbitrary positivity theorem with parameter-dependent degrees
or denominators does not fill that gap.

## 6. Decision and limits

Stop searching polynomial lower-triangular congruences for matrix concavity
of the unshifted blocks: the entire class is excluded, at arbitrary degree.
Do not relabel another coefficient choice or a higher Hankel order as a
new attempt at that route. The canonical scalar rational boundary remains
the next precise construction test; demand a uniform identity or a genuine
obstruction, not a fifth-order numeric fit alone. Literature audit remains
necessary before counting any result as a main contribution.

The proposition is useful exclusion evidence, not an extra Q1 paper or a
substitute for a substantial representation theorem. The earlier K4 lift,
general block-size lower bound, and all unfavorable computation records
remain unchanged. New optimization calls: zero. New closed loops: zero.
