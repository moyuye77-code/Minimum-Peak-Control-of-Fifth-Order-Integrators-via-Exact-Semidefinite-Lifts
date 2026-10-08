# Saturated-input normals and a block-size obstruction

2026-09-28. Research proposition with a written proof; not an independently
reviewed theorem, a novelty clearance, or a completed control paper.

## 1. Object and exact scope

For a positive integer d, use the original moment coordinates

    K_d = { y : y_j = integral_0^1 t^j rho(t) dt, j=0,...,d-1,
                  rho measurable, 0 <= rho <= 1 almost everywhere }.

The semidefinite extension degree sxdeg(C) is the smallest positive integer r
such that C has an exact affine lift over finitely many real PSD blocks of
order at most r. Set it to infinity if no finite PSD lift exists. Arbitrarily
many auxiliary variables and blocks are allowed; a nonlinear change of the
output coordinates is not an affine lift. Scalar inequalities are order-one
blocks, and unrestricted variables can be differences of nonnegative scalars.

**Proposition.** For every d >= 1,

    sxdeg(K_d) >= ceil(d/2).

Consequently K_d has no exact lift over any finite product of Lorentz cones
when d >= 5. This rules out *all* finite SOCP lifts in that range, not just a
particular DCP expression or the old eight-cone construction. It does not rule
out a finite SDP lift with larger blocks, approximations, root-based support
evaluation, or nonlinear formulations.

The proof below connects classical conic duality to an existing polynomial
extension-degree lower bound. We do not claim either classical ingredient as
new. Whether this consequence for bounded-density moment bodies is new remains
under audit.

## 2. Basic geometry and saturated-input normal

K_d is convex and compact: the order interval [0,1] in L-infinity is weak-star
compact, and its finitely many integrals against L-one monomials are continuous
in that topology. Define

    v_j = 1/(j+1),  c = v/2,  C_d = K_d-c,
    p_lambda(t) = sum_{j=0}^{d-1} lambda_j t^j.

Pointwise maximization over rho gives

    h_C_d(lambda) = (1/2) integral_0^1 |p_lambda(t)| dt.

For nonzero lambda the polynomial is not identically zero, hence this support
is positive. Its minimum over the Euclidean unit sphere is positive by
continuity and compactness. Thus C_d contains a ball about zero; it is a
full-dimensional convex body. In particular v belongs to K_d and corresponds
to rho=1.

With the outward-normal convention

    N_K_d(v) = { lambda : lambda dot (y-v) <= 0 for every y in K_d },

we have the exact coefficient-space identity

    N_K_d(v) = P_{d-1}([0,1]),
    P_m([0,1]) = { lambda : p_lambda(t) >= 0 for every t in [0,1] }.

Indeed lambda dot (y-v) = integral p_lambda(rho-1). Nonnegative p_lambda makes
this nonpositive. Conversely, a negative value of a continuous polynomial
gives a negative interval of positive length (also when the value occurs at an
endpoint). Set rho=0 there and rho=1 elsewhere to obtain a strict positive
value. This proves both directions, including the zero polynomial.

An equivalent check uses the polar face

    F = {lambda in C_d polar : lambda dot c = 1}
      = {lambda : p_lambda >= 0 on [0,1], integral p_lambda = 2}.

The identity follows from integral |p| <= 2 and integral p = 2. Its conic hull
is P_{d-1}([0,1]). The main proof does not rely on exchanging closure with a
projection or on an unproved tangent-cone limit.

## 3. Transfer lemma, with the missing regularity issue handled

**Lemma.** If a nonempty convex set C has a finite affine PSD lift with blocks
of order at most r, then for any v in C its normal cone N_C(v) has such a lift
with block orders at most max(r,1).

Write the initial representation as

    C = {L X + c0 : A X = b, X in Q},
    Q = product_i S_+^{r_i},  r_i <= r.

Do not assume the initial lift has Slater points. For each block, let W_i be
the span of the ranges of all feasible X_i. Finitely many feasible matrices
suffice to span these finite-dimensional W_i. Average all the corresponding
feasible tuples. The averaged tuple has range W_i in every block: the kernel
of a sum of PSD matrices is the intersection of their kernels. Every feasible
tuple lies in the product face supported on these W_i. Compress to orthonormal
bases of W_i, discarding zero-dimensional blocks. The same set C is now
represented with PSD blocks of no larger order and a strictly feasible tuple.
Below A,L,Q denote these compressed data.

For any lambda in N_C(v), the primal conic maximization

    maximize <L*lambda, X> subject to A X=b, X in Q

has finite value lambda dot (v-c0), attained by any preimage of v. Strict primal
feasibility gives an attained conic dual of the same value. Since Q is
self-dual, this yields the representation

    N_C(v) = {lambda : there exists eta with
                      A*eta - L*lambda in Q,
                      b dot eta = lambda dot (v-c0)}.                 (1)

The reverse implication follows directly from weak duality: for every feasible
X, lambda dot (LX+c0-v) = -<A*eta-L*lambda,X> <= 0. Thus (1) is exact. Its
constraints are linear plus the same PSD block sizes; free lambda and eta can
be encoded with order-one blocks if needed. No assertion about attainment of
the *old unreduced* support dual is used. This is compatible with earlier
boundary examples where an unreduced dual failed to attain its infimum.

This is a standard duality consequence, written explicitly for the present
application. The proper-lift/Slater argument is also the mechanism in
Gouveia--Parrilo--Thomas [GPT], Theorem 2.4 and the following discussion of lifts
to faces. Their Proposition 2.8 records related polar/face closure properties.

## 4. Imported polynomial lower bound and its application

Let k=floor((d-1)/2). For k>=1, consider the cone P_{2k}([0,1]). It is closed
and convex as an intersection of evaluation halfspaces, pointed because p and
-p can both be nonnegative only if p=0, and full-dimensional because the
constant polynomial one stays positive under small coefficient perturbations.
The interval [0,1] has nonempty interior. Averkov's
Theorem 2 [A, arXiv v2 numbering] implies sxdeg(P_{2k}([0,1]))>k: choose
arbitrarily large finite sets S of distinct points in (0,1), and, for every
k-element subset T of S, take

    f_T(t) = product_{a in T} (t-a)^2.

This polynomial has degree 2k, is nonnegative on the interval, vanishes at T,
and is strictly positive at S outside T. These are exactly the theorem's
hypotheses. This witness construction is classical too; see Remark 29 in [A].
The arbitrarily-large-S quantifier is essential. No finite computation below
establishes the universal lower bound; that part relies on [A].

If d-1=2k+1, intersect P_{d-1} with the linear condition that the coefficient
of t^{2k+1} is zero. This slice is P_{2k} and cannot require larger lift blocks
than P_{d-1}. For d-1=2k there is no slicing step. Section 2 and the transfer
lemma therefore give

    sxdeg(K_d) >= sxdeg(P_{d-1}) >= sxdeg(P_{2k}) >= k+1 = ceil(d/2).

For d=1,2, the displayed bound is just the definitional lower bound one.

Any finite-dimensional Lorentz cone has a lift over finitely many 3D Lorentz
cones: introduce successive partial Euclidean norms, bounding each pairwise
norm by the next auxiliary variable. The forward implication follows from
monotonicity of the norm and nonnegativity of the auxiliaries; the reverse
uses the exact partial norms. The 3D cone is linearly isomorphic to S_+^2 via

    (a,b,c) -> [[a+b,c],[c,a-b]].

One- and two-dimensional Lorentz cones are polyhedral. Hence any finite SOCP
lift would imply sxdeg(K_d)<=2, contrary to the bound for d>=5. The exclusion
allows cones of arbitrary dimension, not only 3D cones chosen in the old code.

## 5. Control consequence and low-order boundary

For a scalar dth-order integrator at a fixed T>0, fixed initial state x0, and
|u|<=J with J>0, write rho(t)=(u(Tt)/J+1)/2. In derivative-reversed coordinates,
the input contribution to state x_{d-k}(T), k=0,...,d-1, is

    J T^{k+1}/k! * [2 sum_{j=0}^k (-1)^j binom(k,j) y_j - 1/(k+1)].

This is a triangular affine map of y with nonzero diagonal
2 J T^{k+1}(-1)^k/k!. Adding the deterministic initial-state flow preserves
invertibility. Thus the exact full-state reach set has the same extension
degree as K_d. Arbitrary positive horizons and nondegenerate scalar input
intervals therefore inherit the lower bound.

Independent input-box chains inherit at least the bound of their largest-order
component by projecting onto it. This does not establish the same bound for
every partial observation, every conditional history, a coupled input ball,
or a nonlinear drone model. An observation can collapse the set to a point.
A universal exact SOCP estimator covering all histories and including the
unobserved full-state case at d>=5 would, however, contradict this proposition.

The resulting current classification is deliberately incomplete:

| order d | current justified conclusion |
|---|---|
| 1 | polytope; order-one lift |
| 2,3 | finite SOC lifts available from the classic moment conditions |
| 4 | this argument gives only block-size >=2; SOC existence unresolved here |
| >=5 | no finite SOC lift; any finite PSD lift needs a block >=ceil(d/2) |

K_2 is not polyhedral: maximize the polynomial s-t, for 0<s<1, to expose the
distinct point (s,s^2/2), using rho=1_[0,s]. Thus it has infinitely many exposed
points. Projection K_3 -> K_2 excludes a polyhedral K_3 lift too. Their minimal
block sizes are exactly two, but this elementary observation is not a new main
contribution. The new lower bound does not decide if a larger-block SDP exists
or whether the bound is sharp in higher order.

## 6. Source and novelty audit

[A] G. Averkov, *Optimal Size of Linear Matrix Inequalities in Semidefinite
Approaches to Polynomial Optimization*, SIAM J. Applied Algebra and Geometry
3(1), 128--151 (2019), DOI 10.1137/18M1201342.
https://epubs.siam.org/doi/10.1137/18M1201342
https://arxiv.org/html/1806.08656v2
Read: definitions, Theorem 2, Section 3.3, Lemmas 23--25 and Theorem 2 proof,
Remark 29. The theorem numbering here is v2, not asserted to match the journal.

[GPT] J. Gouveia, P. A. Parrilo, R. R. Thomas, *Lifts of Convex Sets and Cone
Factorizations*, Mathematics of Operations Research 38(2), 248--264 (2013),
DOI 10.1287/moor.1120.0575. https://arxiv.org/pdf/1111.3164
Read: Theorem 2.4 proof, face/properness discussion, Corollary 2.6 and
Proposition 2.8. The normal-cone formula (1) is derived here by that standard
duality mechanism; it is not claimed as a new general lift theorem.

[HH] S. Haddad, A. Halder, *The Curious Case of Integrator Reach Sets, Part I:
Basic Theory*. https://arxiv.org/html/2102.11423
Rechecked Section IV-E: it distinguishes spectrahedra from their projections
and discusses uncertainty about the latter. Our obstruction concerns small
PSD blocks, so it does not resolve the general projected-SDP question. An old
paper's open-question statement does not establish its current open status.

Targeted searches on 2026-09-28 covered integrator reach sets plus extension
degree/normal cones, bounded-density or L-moment bodies plus SOC/spectrahedral
lifts, and zonoids plus Averkov. No directly matching statement was located in
these searches. That is not an exhaustive citation audit and not proof of
priority. [A] is the closest proof-engine prior; [HH] is a direct control-object
prior. Generic modern spectrahedral state-estimation operations remain prior
art, not new contributions of this proposition.

**Decision:** retain as a structural research candidate, not a cleared Q1/TAC
core. Unlike adding a fourth-moment formula, it constrains every system order
and every finite SOC lift in the stated class. Its short derivation from strong
existing results may still limit standalone significance. A complete original
representation boundary, a sharp constructive counterpart, or a substantial
certified-computation consequence must be established before a paper is
declared ready. Neither tests nor extra communication heuristics fill that gap.

## 7. What the software checks, and what it cannot check

The small sidecar uses exact rational arithmetic to check witness polynomials,
their zero/positive incidence pattern, polar-face normalization, sign-changing
support examples, and the integrator affine-map identity against symbolic
integration. It does not numerically decide nonrepresentability. It runs no
optimizer and no new closed loop; previous solver records remain unchanged.

## 8. Later audit of the entire ordinary-normal route at order four

The proof in ../four_soc_gate/THEORY.md identifies the normal cone at any
boundary density with k genuine interior switches as the linear image
W*P_{d-1-k}([0,1]), where W has exactly those switches and the density's sign
orientation. This handles every boundary point, including the saturated
ones above. For d=4 every such cone has an SOC lift, so selecting a different
ordinary normal cannot strengthen this argument to a block-order-three
lower bound. It does not prove that K4 itself has an SOC lift and does not
rule out global extension-complexity methods. The general lower bound in
this document is unchanged; the sharper K4 question remains unresolved.
