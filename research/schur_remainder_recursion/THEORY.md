# Domain-weighted Schur remainder reduction and the first non-SOC SDP case

2026-09-28. Internally proved research results, not cleared priority or a
completed journal paper. Classical exponential moments, rational moment
lifting, and shadow-closure operations retain their original attribution.
An all-order certificate induction is NOT proved below.

## 1. General odd-order algebra: eliminate one coordinate and one determinant

For n>=1 write x=(y0,...,y_(2n-1)), m=y0, t=y_(2n-1), and w for the
remaining 2n-1 coordinates. Define the classical exponential coefficients

    b_j=[z^(j+1)](1-exp(-sum_k y_k z^(k+1))).

The rational interior lower boundary F_n=L_(2n+1) is the Schur boundary of
the unshifted Hankel block. Write

    A=[b_(i+j)]_(i,j=0,...,n-1) = [[A0,a],[a',c]],
    b=[b_n,...,b_(2n-1)]' = [v,t+d]',
    b_(2n)=y_(2n)-m*t+e(w).

Here A,a,c,v,d,e depend on w alone. Put D_j=det [b_(i+k)]_(i,k=0,...,j-1),
D_0=1, and work on the open determinant domain D_1,...,D_n>0, so A>0.
The needed formulas, with empty-block terms zero for n=1, are

    s=c-a'A0^(-1)a=D_n/D_(n-1),
    eta=a'A0^(-1)v-d-m*s/2,
    H_n(w)=(v+(m/2)a)'A0^(-1)(v+(m/2)a)-m*d-m^2*c/4-e.

Block Gaussian elimination and completing one square give the identity

    F_n(w,t)=H_n(w)+(t-eta(w))^2/s(w).                 (1)

Thus H_n has denominator dividing D_(n-1), whereas F_n has denominator
dividing D_n. This is an identity for every n, not a pattern inferred from
low-order runs. The implementation calibrates it at n=1,2,3. Importantly,
H_n is not in general F_(n-1): (1) is not already a closed induction.

## 2. Translation makes the square stationary at any reference

The time-translation map on raw integral moments is the linear map

    (T_c y)_j=sum_(k=0)^j binom(j,k)(-c)^(j-k)y_k.

Its action on exponential coefficients is the same binomial map. Indeed
substituting z/(1+c z) in the generating series gives precisely that rule.
The leading Hankel blocks transform by unit lower-triangular binomial
congruences, so every D_j is unchanged. This is algebraic; T_c need not keep
the physical support interval [0,1] fixed. The larger determinant domain is
translation invariant, and that is the domain used for this reduction.

Applying determinant invariance also to the next Hankel block gives

    F_n(T_c x)=F_n(x)+sum_(k=0)^(2n-1) binom(2n,k)(-c)^(2n-k)y_k. (2)

For fixed c the difference is affine in x. Its coefficient on t is -2n c.
Since derivative_t F_n=2(t-eta)/s, at a reference u=(w0,t0) choose

    c_u=(t0-eta(w0))/(n*s(w0)).

Then derivative_t F_n(T_(c_u)u)=0. If x'=T_(c_u)x and u'=T_(c_u)u,
the square in (1) has zero value AND zero gradient at u'. Affine terms
and fixed linear coordinate changes preserve first-order remainders, hence

    R_F_n(x,u)=(t'-eta(w'))^2/s(w')+R_H_n(w',w0').     (3)

When differentiating in x, c_u is held fixed. Forgetting this would give an
incorrect remainder. Its eventual dependence on u is rational with only
leading determinant denominators.

### A restricted certificate equivalence, not an arbitrary-lift equivalence

Consider joint first-order SOS certificates whose domain weights are products
of the leading D_j at x and u, and whose separate clearing multipliers p(x),
q(u) are monomials in the D_j. Exponents and square degrees may depend on n
but are fixed for all points and reference points at that n.

**Lemma.** F_n has a certificate in this class if and only if H_n does.

For the forward construction from H_n to F_n, substitute the rational
translation in its certificate. All D_j weights and p(x) remain invariant,
so no mixed x,u denominator enters p. Clear the finitely many u-only
denominators with even powers of the D_j(u). Add the square in (3), using
1/s=D_(n-1)/D_n and clearing the determinant denominators of eta. Multiplying
by further determinant powers preserves the preordering: even powers enter
squares, odd powers enter the generator products. Thus the resulting square
space is finite and joint-polynomial, with separate determinant multipliers.

Conversely, substitute t=eta(w), t0=eta(w0) into a certificate for F_n.
The t-gradient at the reference is zero, so its remainder becomes R_H_n.
The D_j do not depend on t. Clearing the finitely many denominators of eta
by even determinant powers again gives a certificate of the specified type.

The restriction to invariant determinant multipliers is substantive. For an
arbitrary p_H, p_H(T_(c_u)x) could depend jointly on x and u and need not
admit separate clearing. The lemma therefore does not decide all rational
certificates or all affine SDP lifts. It supplies a precise uniform reduction
target; H_n remains an unsolved family beyond the base cases below.

## 3. The n=2 reduced certificate closes exactly

For x=(m,z,q,r), let Delta=mq-z^2-m^4/12 and N=mr-zq-m^3z/6. Then

    L5=H2(m,z,q)+N^2/(m Delta),
    H2=q^2/m+m^2q/12+mz^2/4-m^5/720.

Use u=(M,Z,Q,R), h=m-M, Delta0=Delta(u), N0=N(u), and kappa=Delta/m.
An exact six-term certificate for the reduced remainder is

    R_H2 = [q-mQ/M+mMh/12]^2/m
           +[Mz+(m-2M)Z]^2/(6M)
           +[(2m-M)z-mZ]^2/(12m)
           +h^2*kappa(x)/12+h^2*kappa(u)/6
           +h^4*(3M+2m)/360.                         (4)

It is nonnegative whenever m,M,Delta,Delta0>0. It follows by completing the
q square after subtracting h^2*kappa(x)/12+h^2*kappa(u)/6. The remaining
z,Z quadratic splits into the second and third squares in (4), and the
remaining scalar is exactly h^4*(3M+2m)/360. All identities are checked
symbolically in independent variables, not by fitting an SOS Gram matrix.

For this case c_u=N0/(2 Delta0). Set

    A=Delta0*N-N0*Delta,
    B=Delta0*(Mq-mQ+mM^2h/12)-N0*(Mz-mZ),
    C=Delta0*(Mz+(m-2M)Z)-N0*M*h,
    D=Delta0*((2m-M)z-mZ)-N0*m*h,
    p(x)=m*Delta, p(u)=M*Delta0.

Equations (3)--(4), after clearing separate denominators, give

    p(x)p(u)^2 R_L5(x,u)
      = (MA)^2 + Delta*B^2
        +(m Delta)*M*C^2/6 + Delta*(MD)^2/12
        +(Delta M Delta0 h)^2/12
        +(m Delta)*(M Delta0)*(Delta0 h)^2/6
        +(m Delta)*M*(M Delta0 h^2)^2/120
        +Delta*(m M Delta0 h^2)^2/180.                (5)

This is a finite domain-weighted certificate, with generators 1,m,Delta,
m Delta at each endpoint. The maximum degree in x is 10. The negative
unweighted witness from the prior turn is not contradicted: (5) changes the
clearing multiplier and includes domain weights. The initial development
expansion omitted the M^2 factor in the first term; exact cancellation
detected the error, which was corrected before any archive or theorem claim.

## 4. An explicit finite affine lift of the open one-sided epigraph

Let Z be a linear functional on polynomials in four formal coordinates
(m,z,q,r) of total degree at most 10. There are 1001 unreduced coefficients.
Let V_j denote all monomials of degree <=j. Impose

    Z(V5 V5')>=0,          order 126,
    Z(m V4 V4')>=0,        order 70,
    Z(Delta V3 V3')>=0,    order 35,
    Z(p V2 V2')>=0,        order 15,
    Z(p)=1,
    output x_i=Z(p*x_i),
    output v>=Z(p*L5).                               (6)

p*L5 is a polynomial of degree <=10. All constraints in (6) are affine
equalities/inequalities and finite affine PSD blocks; they are not nonlinear
constraints on the output. This is intentionally unreduced, not a size or
speed claim. Its exact projection is

    E5={ (x,v): m>0, Delta(x)>0, v>=L5(x) }.          (7)

For any point of (7), the functional Z(f)=f(x)/p(x) makes every matrix a
positive multiple of an outer product and satisfies all constraints.

For the converse, first prove that a feasible functional's output is inside
the domain, without assuming a representing measure. Write its output as
(M,Z0,Q,R). The Delta-weighted block on (1,m) is

    [[Z(Delta),1],[1,M]],

so M>0. The unweighted block on (m,Delta) has off-diagonal Z(p)=1, so
Z(Delta^2)>0. The Delta-weighted block on (m,z) gives
Z(Delta*z^2)>=Z0^2/M. Finally, with Lambda(f)=Z(p f),

    Lambda(m^3)-M^3
      =Z(Delta*[m(m-M)]^2)+2M*Z(p*(m-M)^2)>=0.

Normalization and the output equations make the Taylor linear terms vanish.
Using p*(q-z^2/m-m^3/12)=Delta^2 gives the exact identity

    Q-Z0^2/M-M^3/12
      =Z(Delta^2)+[Z(Delta*z^2)-Z0^2/M]
         +[Lambda(m^3)-M^3]/12 > 0.                 (8)

Thus Delta(output)>0. Apply Z to (5), taking u to be that output. Every
square is controlled by the displayed blocks, and all reference generators
are now positive. The linear remainder vanishes, proving
Z(p L5)>=L5(output). This proves both directions of (7).

The rational lift does NOT contain zero denominators. In particular no
finite auxiliary can represent the zero prefix in (6). A point-evaluation
functional along density epsilon/2 has Z(1)=1/p(x_epsilon) tending to
infinity. Singular outputs require a separate closure construction.

## 5. Full K5, including its degenerate boundary

The previously proved finite affine K4 lift can be used to represent
int K4 by a standard central homothety. Fix its center c=(1/2,1/4,1/6,1/8).
Write a K4 lift as A0+sum y_i A_i+sum z_j B_j>=0, including all its blocks.
Use epsilon+tau=1 and epsilon,tau>0, encoded by two 2x2 blocks
[[epsilon,1],[1,a]]>=0 and [[tau,1],[1,b]]>=0. Require

    tau*A0+sum (x_i-epsilon*c_i) A_i+sum w_j B_j>=0.

Its projection is precisely x=epsilon*c+tau*y with y in K4, hence int K4.
Every interior x has this expression for sufficiently small epsilon>0;
conversely mixing any K4 point with its interior center stays interior.
All entries above are affine. No multiplication of unknown x and epsilon
is present, and no zero-tau homogenization issue is used.

Intersect this interior-prefix representation with (6) for (x,v) and a
second copy for (bar x,1/5-v), where

    bar x=(1,1/2,1/3,1/4)-x.

The classical lower/upper boundary description makes this intersection
exactly the K5 band above int K4. Both required determinants are positive
there. The full K5 is its closure: mix any feasible density with 1/2 and
let the mixing weight tend to zero. Conversely K5 is compact and closed.

By Gouveia--Netzer, Corollary 3.4, the closure of a spectrahedral shadow is
again a spectrahedral shadow. This separately established operation proves
that full K5 has a finite affine SDP lift. It does not permit plugging a zero
denominator into (6), nor does it preserve the stated open-model block sizes
by assertion. We do not give or claim a 126-block size bound for full K5.

Combining this construction with the previous proved general lower bound
sxdeg(K_d)>=ceil(d/2) yields an internally proved category separation:

    K5 has a finite affine SDP lift, but no finite SOCP lift.

This consequence concerns any finite product of Lorentz cones, not just a
particular solver syntax. It is not a new lower-bound proof, a minimal SDP
size theorem, an all-order SDP theorem, or cleared originality.

## 6. What would complete the general mechanism, and what is not required

To complete this construction route one needs determinant-module certificates
for the reduced family H_n for every n, not just Schur identities at more
orders. n=3 already has a 15-term numerator and denominator proportional to
D2; its expression is archived as a calibration, not an SOS certificate.
The reduction lemma does not identify H_n with a previously settled family.
Even-order bodies are affine projections of the next odd-order body, so an
all-odd-order construction would suffice for all K_d, but it is unproved.

The user's publication goal is a substantial, original, coherent contribution,
not an artificially imposed demand to solve every order. The K5 category
separation is therefore a concrete candidate to audit for independent value;
it is not counted as a second paper or added to a list of small innovations.
Its priority, significance, practical representation cost and control evidence
remain to be assessed before choosing it as the manuscript's central result.
The current paper is not declared complete. No new optimization or closed-loop
experiment is claimed here, and older unfavorable performance data remain.
