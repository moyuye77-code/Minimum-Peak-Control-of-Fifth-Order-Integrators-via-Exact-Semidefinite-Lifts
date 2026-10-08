# A joint affine SDP lift of the four bounded-density moments

2026-09-28. This is an internally proved construction, not a claim of priority,
an SOC representation, or a completed journal paper. It replaces the previous
"variable-m SDP unresolved" status for K4, but not the all-order SDP question.

## 1. Statement and explicit matrices

Let

    K4 = {(m,z,q,r) = integral_0^1 (1,t,t^2,t^3) rho(t) dt : 0<=rho<=1}.

**Theorem.** A point (m,z,q,r) belongs to K4 if and only if there exist six
real auxiliary variables (a,b,c,d,e,f) for which all seven matrices below are
positive semidefinite. Every entry is affine in the ten variables.

    H = [[1,m,a], [m,a,b], [a,b,c]],

    J_sigma = [[1, m, z+sigma*a/6],
               [m, a, d+sigma*b/6],
               [z+sigma*a/6, d+sigma*b/6, f+sigma*e/3+c/36]], sigma=+1,-1,

    T_minus = [[z-a/2, d-b/2], [d-b/2, e-c/2]],
    T_plus  = [[m-z-a/2, a-d-b/2], [a-d-b/2, b-e-c/2]],

    G_minus = [[z-a/2, q-d/2-b/12],
               [q-d/2-b/12, r-f/2-e/4]],

    B_plus = -a/4+b/12-d/2+m/2-q+z/2,
    C_plus = -a/8-d/2+e/4-f/2+m/4-r+3*z/4,
    G_plus = [[m-z-a/2, B_plus], [B_plus, C_plus]].

Thus K4 has an affine lift with three PSD blocks of order three, four PSD
blocks of order two, and six auxiliary scalar variables. No additional
prefix constraints or constraints fixing m are needed. In particular,

    2 <= sxdeg(K4) <= 3.

The lower bound is the earlier general result, not a consequence of the
present construction. Whether the lower or upper endpoint is exact remains
unresolved here: this does not exclude a different SOC lift of K4.

## 2. Equivalent nonlinear boundary conditions, including singular points

Define

    D = z-m^2/2, E = m-z-m^2/2,
    B = q-m*z/2-m^3/12, F0=z^2/2+m^2*z/4,
    G(m,z,q,r) = [[D,B],[B,r-F0]].

Use barred coordinates

    (mbar,zbar,qbar,rbar)=(1-m,1/2-z,1/3-q,1/4-r).

We first establish the exact, but NONLINEAR, description

    K4 = {y : G(y)>=PSD 0, G(ybar)>=PSD 0}.                 (1)

The lower boundary already proved in moment_upper_gate/THEORY.md equals

    L(m,z,q) = F0 + B^2/D,  D>0.                           (2)

The equality to that previous expression follows by expansion. Set

    qlo=z^2/m+m^3/12,
    qhi=1/3-(1/2-z)^2/(1-m)-(1-m)^3/12.

For 0<m<1 and D,E>0, direct algebra gives

    qhi-qlo = D*E/[m*(1-m)],
    L(y_prefix)+L(ybar_prefix)-1/4
       = m*(1-m)*(q-qlo)*(q-qhi)/(D*E).                    (3)

If both matrices in (1) are PSD, their leading entries imply D,E>=0, hence
m^2<=m and 0<=m<=1. In the strictly positive case, their Schur complements
and (3) force qlo<=q<=qhi. The first three moments are therefore in K3, and
the two bounds (2) are exactly the necessary and sufficient K4 bounds proved
in the preceding module. Conversely, every K4 point obeys these matrices.

There is no division-by-zero gap in (1). For 0<m<1 and D=0, a PSD matrix
with a zero leading diagonal entry must have B=0. Thus z=m^2/2 and
q=m^3/3, the prefix density 1_[0,m]. Its complement has positive E; the
other bound forces r<=m^4/4, while the first forces r>=m^4/4. The case E=0
is symmetric. If m=0, D,E>=0 force z=0; the two matrices then force q=r=0.
The case m=1 follows by complementation and gives (1,1/2,1/3,1/4).

## 3. The identity that prevents the affine lift from relaxing K4

For any real parameter g, put

    F_g(s,t) = t^2/2+s^2*t/4
               -g*(s^3/12+s*t/2)+g^2*(s^2/8-t/4).

Let a reference point be (M,Z), write u=s-M, v=t-Z, and
D(s,t)=t-s^2/2, D_ref=Z-M^2/2. The following polynomial identity holds
for all real s,t,M,Z,g, without any genericity assumption:

    F_g(s,t)-F_g(M,Z)-grad F_g(M,Z) dot (u,v)
      = 1/2*[v+(M-g)*u/2+u^2/6]^2
          +D(s,t)*u^2/12+D_ref*u^2/6+u^4/36.              (4)

It can be verified by expanding coefficients; exact.py checks the identity
symbolically. If D(s,t),D_ref>=0, the right side is nonnegative. More
importantly, the square factors occupy small fixed polynomial spaces:

    u, u^2 belong to span{1,s,s^2};
    v+(M-g)*u/2+u^2/6 belongs to span{1,s,t+s^2/6}.

This permits a finite positive-functional argument, not just a convexity
assertion. Ordinary convexity alone would not establish any finite SDP lift.

## 4. Sufficiency of the seven affine LMIs

Given a feasible tuple, define a LINEAR FUNCTIONAL Lambda on the span of

    1,s,t,s^2,s^3,s^4,s*t,s^2*t,t^2

by assigning its values

    1,m,z,a,b,c,d,e,f,

respectively. Lambda is not assumed to possess a representing measure.
These are coordinate-polynomial auxiliaries, not extra time moments of rho.

The matrix H certifies Lambda[h^2]>=0 for h in span{1,s,s^2}.
J_plus certifies the same on span{1,s,t+s^2/6}.
T_minus certifies Lambda[D(s,t)*h^2]>=0 for h in span{1,s}.

First, H implies a>=m^2. The leading entries of T_minus and T_plus imply
z>=a/2 and m-z>=a/2, so a<=m. Therefore 0<=m<=1 and

    D_ref=z-m^2/2>=0, E_ref=m-z-m^2/2>=0.

Apply Lambda to (4) with (M,Z)=(m,z). The linear Taylor term vanishes.
All four terms on the right are nonnegative by H, J_plus, T_minus and
D_ref>=0. Hence the following Jensen inequality holds for every real g:

    Lambda(F_g) >= F_g(m,z).                               (5)

The quadratic form of G_minus evaluated at (-g/2,1) says

    r >= g*q+Lambda(F_g) >= g*q+F_g(m,z), for every g.       (6)

For D_ref>0, taking the supremum over g gives (2). If D_ref=0,
finiteness of (6) forces B=0 and r>=F0. Thus (6), together with
D_ref>=0, is equivalent to G(m,z,q,r)>=PSD 0, including its singular face.

Now use sbar=1-s, tbar=1/2-t. Polynomial spaces transform as follows:

    span{1,sbar,sbar^2} = span{1,s,s^2},
    span{1,sbar,tbar+sbar^2/6} = span{1,s,t-s^2/6},
    D(sbar,tbar) = s-t-s^2/2.

Consequently H, J_minus and T_plus prove (5) in the barred coordinates.
Their coordinate moments are

    abar=1-2m+a, bbar=1-3m+3a-b,
    cbar=1-4m+6a-4b+c,
    dbar=1/2-z-m/2+d,
    ebar=1/2-m+a/2-z+2d-e,
    fbar=1/4-z+f.

G_plus is exactly G_minus's affine formula with these barred variables.
Repeating (6) therefore establishes G(ybar)>=PSD 0. Equation (1) now
proves that the projected point belongs to K4. This proof uses neither
Slater's condition nor a closure-of-projection argument.

## 5. Necessity and the uncompressed comparison

For any y in K4 choose the atomic coordinate functional

    (a,b,c,d,e,f)=(m^2,m^3,m^4,m*z,m^2*z,z^2).

H and J_sigma are outer products. T_minus=D*[1,m]*[1,m]^T and
T_plus=E*[1,m]*[1,m]^T. The final two matrices equal the exact nonlinear
matrices in (1). Every LMI holds. This proves the reverse inclusion.

An easier but larger variant replaces H,J_plus,J_minus by

    M4 = [[1,m,z,a],[m,a,d,b],[z,d,f,e],[a,b,e,c]] >=PSD 0.

It certifies every square in span{1,s,t,s^2}. The three compressed Gram
matrices are congruences of M4. The proof above shows that only those three
restricted spaces are needed; no equivalence of arbitrary Gram cones is
claimed. Both variants have the same K4 projection and six auxiliaries.

## 6. Control consequence and limits

For a fourth-order scalar integrator with bounded measurable input |u|<=J,
fixed J,T>0, a known initial state, and u(T*t)=J*(2*rho(t)-1), the input
increment in the component k integrals below the input, k=0,1,2,3, is

    J*T^(k+1)/k! * [2*sum_(j=0)^k (-1)^j binomial(k,j)*y_j - 1/(k+1)].

This is an invertible affine map, so the full one-segment reachable set
has the stated exact finite SDP lift. It imposes no input time mesh and
does not fix the total input. Independent input-box chains follow by products.
Multi-segment dynamics, linear observations and initial sets with known
SDP lifts can be joined using standard affine/product/intersection operations.
Those closure operations are not claimed new. Finite observed knots can be
joined without adding a time discretization of the unknown input.

There is no assertion of exact Euclidean input-ball coupling, nonlinear
quadrotor dynamics, inter-knot state safety, communication optimality, or
floating-point certification. Exact representability is a mathematical
statement; a numerical solver may return inaccurate or infeasible results,
especially on singular faces. The present numerical audit is not a closed
loop and does not establish runtime advantage over analytic support oracles.

## 7. Relation to prior tools

The use of finite moment/Gram liftings with first-order nonnegativity
certificates is established methodology, notably Nie's rational-function
and matrix-concavity work. The explicit identity (4), its sparse polynomial
spaces, the resulting seven affine matrices, and the boundary-complete K4
application are the specific candidate result here, not a new generic
semidefinite-representability principle. See PRIOR_ART.md for the reading
scope and unresolved priority audit. Existing tools do not by themselves
either establish or invalidate novelty of this particular application.
