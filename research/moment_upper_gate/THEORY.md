# Exact upper-representation gate: what works and what does not

2026-09-28. These results resolve two construction checks left open in
GENERAL_MOMENT_LIFT_GATE.md. They do not establish an exact SDP lift for general
K_d, and are not being counted as separate new paper contributions. The
all-order lower bound in moment_lift_obstruction remains unchanged.

## 1. A finite domination hierarchy is never the required exact lift

Write ell_j=1/(j+1). For R>=1 define the standard finite moment relaxation

    H_R(y) = [y_(i+j)]_(i,j=0,...,R),
    G_(R-1)(y) = [y_(i+j+1)-y_(i+j+2)]_(i,j=0,...,R-1),
    D_R = {y in R^(2R+1): H_R(y), H_R(ell-y),
                          G_(R-1)(y), G_(R-1)(ell-y) >=PSD 0}.

Every density 0<=rho<=1 gives such a y. This follows by integrating f^2 and
t(1-t)f^2 against rho and 1-rho. This is a standard moment/localizing-matrix
domination construction, not a new hierarchy [L].

**Non-exactness check.** For every finite R>=1, and every 2<=d<=2R+1,
the first d coordinates of D_R strictly contain K_d.

Here is an explicit rational certificate. Set w=1/(R+1)^2 and take

    y=(w,0,0,...,0).                                            (A)

H_R(y)=w e_0 e_0^T and G(y)=0. The complement localizing matrix is a positive
definite Gram matrix for t(1-t) on [0,1]. The remaining condition is

    H_R(ell)-w e_0 e_0^T >=PSD 0.

It holds because every degree-R polynomial f satisfies

    f(0)^2 <= (R+1)^2 integral_0^1 f(t)^2 dt.                    (B)

For completeness, let P_j(t)=(1/j!) (d/dt)^j[t^j(t-1)^j]. Integration by parts
shows that these shifted Legendre polynomials are orthogonal, P_j(0)=(-1)^j,
and integral P_j^2=1/(2j+1). Expanding f=sum a_j P_j and applying
Cauchy--Schwarz gives (B), since sum_(j=0)^R(2j+1)=(R+1)^2. Equivalently the
first diagonal entry of the inverse Hilbert matrix of size R+1 is (R+1)^2;
this is classical Hilbert/Christoffel-kernel algebra, not a new inequality.

No bounded density has the first two moments (w,0): a nonnegative density with
integral t rho=0 vanishes almost everywhere on (0,1], hence has zero mass.
Alternatively K_2 requires y_1>=y_0^2/2. The linear objective p(t)=w-t gives

    value at (A) = w^2,
    max over K_d = integral_0^w (w-t) dt = w^2/2.

The resulting relaxation-gap lower bound is 1/[2(R+1)^4]. The direction here
depends on R; this is not asserted to be the exact relaxation gap or its sharp
asymptotic rate, and is not normalized to unit Euclidean coefficient norm.

There is also a strict gap for one *fixed* objective p(t)=1/2-t at every R.
Let s=1/2, h_j=s^(j+1)/(j+1), and B=H_R(ell-h), a positive definite Gram
matrix on [s,1]. Set

    epsilon = 1/(e_0^T B^(-1) e_0) > 0,
    y=h+epsilon e_0.                                           (C)

All four matrices in D_R remain PSD: H_R(h) increases by a rank-one PSD
matrix; its complement is B-epsilon e_0e_0^T, PSD by Cauchy--Schwarz in the
B inner product; the two G matrices do not change. The true support is 1/8,
whereas (C) attains 1/8+epsilon/2. For rational s every entry is rational.
This construction does not assert that (C) is optimal in D_R.

A qualitative reason extends to every sign-changing polynomial objective of
degree at most 2R:
the true maximizer rho=1_{p>0} gives positive definite H and G for both rho
and 1-rho whenever their supports contain open subintervals. It is therefore
an interior point of every finite D_R, so a nonzero linear objective cannot
be maximized there. This is a strict finite-truncation gap, not a failure of
the full infinite hierarchy, whose characterization in [L] remains valid.

**Consequence for our research:** increasing the order of this particular
hierarchy cannot supply a finite exact matching upper representation, even
for K_2, which already has an exact SOC description. Failure of a hierarchy
is logically different from nonexistence of any SDP lift.

## 2. Completion of the four-moment constructive boundary

Let (m,z,q) be in K_3. For 0<m<1 introduce

    D=z-m^2/2, Q=q-m^3/3, E=m(1-m)-D.

The classic three-moment conditions are equivalent to

    0<=D<=m(1-m),
    m D + D^2/m <= Q <= (m+1)D - D^2/(1-m).                  (P)

If D>0, put b=Q/D-m and

    Delta=(m-b)^2+4D,
    a=(m+b-sqrt(Delta))/2, c=(m+b+sqrt(Delta))/2.

Condition (P) implies D/m<=b<=1-D/(1-m). Thus mb>=D and
(1-m)(1-b)>=D; also 0<=m,b<=1. The roots of

    (t-m)(t-b)-D=0

therefore satisfy 0<=a<=b<=c<=1. The middle ordering follows from
sqrt(Delta)>=|m-b|; the endpoint ordering follows by comparing Delta to
(m+b)^2 and (2-m-b)^2. This closes the ordering gap in the earlier exploratory
calculation, including a=0 or c=1.

The density rho_- = 1_[0,a] + 1_[b,c] has moments m,z,q. Indeed a+c=m+b,
ac=mb-D imply

    a+c-b=m,
    (a^2+c^2-b^2)/2=m^2/2+D=z,
    (a^3+c^3-b^3)/3=m^3/3+D(m+b)=q.

For any other density with these three moments, use

    W(t)=(t-a)(t-b)(t-c).

On rho_-=1, W<=0 and rho-rho_-<=0; on its complement W>=0 and
rho-rho_->=0. Hence integral W(rho-rho_-)>=0. All terms of degree <=2
cancel, so the fourth moment r=integral t^3 rho satisfies r>=L, attained by
rho_-, where direct expansion gives

    L(m,z,q)=Q^2/D - m Q + m^2 D + D^2/2 + m^4/4.           (L4)

This equals the earlier nonlinear Markov expression in the original gate.
If D=0, z=m^2/2 forces rho=1_[0,m] almost everywhere: apply the same sign
argument to t-m to see equality is possible only for that density. Thus
q=m^3/3 and r=m^4/4. Define L accordingly. The cases m=0 and m=1 give the
zero and full densities and are handled separately, with L=0 and L=1/4.

Apply the lower construction to the complement prefix

    (mbar,zbar,qbar)=(1-m,1/2-z,1/3-q)

and complement its density to obtain the upper boundary

    U(m,z,q)=1/4-L(mbar,zbar,qbar).                            (U4)

Both extremizers have the specified prefix. The existence of a K_3 density
and the two necessary bounds ensure L<=U. Convex mixtures of the two
extremizers attain every r in [L,U]. Thus (P), (L4), and (U4), with these
endpoint definitions, characterize K_4 exactly and construct a bounded
density for every point. This is a classical moment extremal construction,
not a claim of a new truncated-moment theorem or a general affine SDP lift.

## 3. A precise positive result: fixed-mass K_4 fibers have SOC lifts

Fix m in (0,1) as DATA, not an optimization variable. Define D,Q as above;
define Dbar=zbar-mbar^2/2 and Qbar=qbar-mbar^3/3. These quantities are affine
in the optimization variables z,q. With four auxiliaries A,B,Abar,Bbar, impose

    D,Dbar>=0,
    [[m,D],[D,Q-mD]] >=PSD 0,
    [[mbar,Dbar],[Dbar,Qbar-mbar Dbar]] >=PSD 0,
    [[D,Q],[Q,A]] >=PSD 0, [[1,D],[D,2B]] >=PSD 0,
    [[Dbar,Qbar],[Qbar,Abar]] >=PSD 0,
    [[1,Dbar],[Dbar,2Bbar]] >=PSD 0,
    r >= A+B-mQ+m^2 D+m^4/4,
    1/4-r >= Abar+Bbar-mbar Qbar+mbar^2 Dbar+mbar^4/4.       (S4)

These are six order-two PSD constraints, hence six 3D SOC constraints, plus
linear inequalities. The first two give the K_3 prefix. The remaining four
give the closed perspectives in (L4) and (U4). Conversely choose
A=Q^2/D, B=D^2/2, and complementary values, using A=0 when D=Q=0.
Then the bound inequalities are exactly the fourth-moment bounds. Degenerate
m=0 or 1 fibers are single points and need only linear constraints.

This is an exact affine lift for each fixed m. It does NOT provide a joint
affine lift for variable m: terms mQ, m^2D, m^4 and definitions D,Q would
then become nonlinear. Taking the union over m is not a valid finite conic
construction. The full K_4 SOC question remains unresolved by this work.

For a fourth-order integrator, m is the normalized total input, hence fixing
it is equivalent to fixing the highest-derivative state increment. An exact
measurement imposing this condition can simplify that particular fiber; it
does not simplify every general or noisy observation history. No new control
performance or communication optimality claim follows from (S4).

## 4. Literature audit and research decision

[L] J. B. Lasserre, *The moment problem with bounded density*,
https://arxiv.org/pdf/math/0607463 (v2, 2007). Read Theorem 3, its proof and
the subsequent finite-matrix interpretation on pp.5--7. It characterizes a
full moment sequence by all orders, not a uniformly finite exact lift for
every truncated vector. The present check is consistent with that result.

[GR] L. Gosse and O. Runborg, *Resolution of the finite Markov moment problem*,
https://arxiv.org/pdf/0910.5791. Read the problem definition, four-step
generalized-eigenvalue reconstruction and Theorem 1 statement. Also inspected
the exponential-transform Lemma 2 and the positivity/interlacing discussion
of https://arxiv.org/pdf/0809.3714; the complete 25-page proof was not audited.
These are direct antecedents for reconstructing interval endpoints, so that
operation and the nonlinear transformed moment matrices are not claimed new.
They are not automatically original-coordinate affine PSD lifts.

Classical Hilbert and orthogonal-polynomial machinery underlies (A)--(C).
W. Kahan's notes https://people.eecs.berkeley.edu/~wkahan/MathH110/HilbMats.pdf
were inspected at the inverse/roundoff discussion (p.5), not read in full.
The checks here use rational arithmetic, not floating-point eigenvalues of
ill-conditioned Hilbert matrices. No novelty is claimed for their inversion.

The earlier Averkov-based lower bound still has no located exact matching
statement in the targeted searches. That is not a complete citation audit.
The present results do not upgrade its priority status and do not establish
that the general SDP question is currently open throughout the literature.

**Decision:** eliminate finite truncations of the ordinary domination
hierarchy as a route to an exact matching upper lift. Retain the explicit K_4
boundary and fixed-m fibers as construction tools. Neither is promoted to a
standalone Q1 core. The unresolved central issue is the variable-m coupling
and general-order representation, not more heuristics or extra experiments.
