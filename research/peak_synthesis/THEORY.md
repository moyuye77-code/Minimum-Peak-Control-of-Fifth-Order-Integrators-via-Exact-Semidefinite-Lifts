# Exact peak-input synthesis as a consequence of the K4 lift

2026-09-28. This is an application of the candidate representation in
`../moment_four_lift/THEORY.md`, not a new minimum-norm control problem or
a separate main contribution. Perfect-spline/minimum-maximum-derivative
problems have classical antecedents; see PRIOR_ART.md. Mathematical exactness
below does not imply exact floating-point solver output.

## 1. Problem and scope

There are A independent fourth-order chains, with known initial states,

    x_a=(p_a,v_a,a_a,j_a),   x'_a=(v_a,a_a,j_a,u_a),   a=1,...,A.

Fix T>0 and finitely many rows (t_i,w_i,b_i,epsilon_i), where
0<t_i<=T and epsilon_i>=0. A row can mix derivatives and different axes.
The problem is

    minimize J = ess sup_{t in [0,T]} max_a |u_a(t)|,
    subject to |w_i^T x(t_i)-b_i| <= epsilon_i for every i.

Inputs range over all essentially bounded measurable functions. The initial
state is a parameter, not an optimization variable. Times are fixed. This is
a componentwise peak bound, not a Euclidean ball coupling simultaneous inputs.
The rows constrain only the listed times. No inter-waypoint state, collision,
thrust, attitude, disturbance-rejection or closed-loop guarantee is asserted.
The fourth derivative of position is the input; an L-infinity objective is
different from the usual integral-squared minimum-snap objective.

## 2. Homogenized moment cone, including its apex

Recall K4={integral_0^1 (1,s,s^2,s^3)rho(s)ds: 0<=rho<=1}.
Write the seven affine pencils of the previous theorem as M_k(y,zeta),
where y=(m,z,q,r), zeta=(a,b,c,d,e,f). The only nonzero constant entries
are the upper-left 1 in H and in each J_sigma.

Define Mhat_k(J,Y,Z) by changing precisely these three constants from 1
to J and replacing y,zeta with Y,Z. All entries remain affine. Then

    C = {(J,Y): J>=0, exists Z, Mhat_k(J,Y,Z)>=0 for all k}
      = {(J,J*y): J>=0, y in K4}.                         (1)

For J>0 this follows by division: Mhat_k=J*M_k(Y/J,Z/J).
It is essential to check J=0 instead of assuming a homogenized projection
has no spurious recession points. At J=0, use the fact that a zero diagonal
entry of a PSD matrix forces its entire row and column to vanish:

1. H forces m=a=0, then b=0; its last diagonal gives c>=0.
2. J_sigma forces z=0. The top row of T_minus now forces d=0.
3. T_minus and T_plus give e-c/2>=0 and -e-c/2>=0, hence c=e=0.
4. The last diagonal of J_sigma gives f>=0. G_minus forces q=0.
5. The last diagonals of G_minus and G_plus give r-f/2>=0 and
   -r-f/2>=0, hence f=r=0.

Thus all moments and auxiliaries vanish at the apex; conversely the all-zero
tuple is feasible. This proves (1) without a Slater or positive-peak assumption.

## 3. Finite SDP equal to the continuous synthesis problem

Partition [0,T] only at the prescribed row times, obtaining N segments [l,r].
For each axis and segment introduce four Y moments and six auxiliaries with
the cone in (1), sharing the same scalar J across all segments and axes.
For h=r-l and s=(t-l)/h, write

    v_a(s)=(u_a(l+h*s)+J)/2,
    Y_j=integral_0^1 s^j v_a(s)ds, j=0,...,3.

Then 0<=v_a<=J if and only if the homogenized cone is satisfied, in the
sense of existence of a representing input. Its signed input moments are

    integral_0^1 s^j u_a(l+h*s)ds = 2Y_j-J/(j+1).         (2)

Remove the known zero-input contribution exactly. Put

    beta_i = b_i-w_i^T exp(F*t_i)x(0),
    k_ia(t) = sum_{d=0}^3 w_{i,4a+d}(t_i-t)^(3-d)/(3-d)!
              for t<=t_i, and 0 otherwise.

Indices here use d=0 for position; the matrix F is the block chain dynamics.
On each preceding segment expand

    h*k_ia(l+h*s) = sum_{j=0}^3 c_iaj^(l,r) s^j.

The exact finite row is

    |sum_{a,segments,j} c_iaj^(l,r)[2Y_aj^(l,r)-J/(j+1)]-beta_i|
       <= epsilon_i.                                   (3)

Minimize J subject to J>=0, all seven PSD blocks per axis/segment, and (3).
There are 10AN+1 scalar variables before solver canonicalization, with
3AN order-three and 4AN order-two PSD blocks. There is no control time mesh.

Proof of equivalence: every feasible measurable input gives the moment
variables and feasible auxiliaries in (1), and convolution gives (3).
Conversely (1) supplies a bounded measurable input on each axis/segment.
Concatenate these inputs; convolution defines a continuous state trajectory
with the fixed initial state and satisfies all rows. Values at a finite
number of segment endpoints do not affect essential bounds or integrals.
At J=0 the apex proof gives the unique zero input almost everywhere.
Thus both feasible peak sets and optimal values agree. No assertion about
numerical feasibility follows solely from a solver status.

If the continuous problem is feasible, the infimum is attained: restrict
inputs to a bounded minimizing sequence, use weak-star compactness of its
L-infinity ball, and pass the finitely many L1-kernel integrals to the limit.
The finite SDP has a corresponding feasible optimizer by (1).

## 4. Continuous dual and independently verifiable bounds

For arbitrary real lambda, put

    k_a(t)=sum_i lambda_i k_ia(t),
    I(lambda)=sum_a integral_0^T |k_a(t)|dt,
    B(lambda)=sum_i lambda_i beta_i-sum_i epsilon_i |lambda_i|.

Every feasible input of peak J satisfies B(lambda)<=J*I(lambda).
Consequently a rigorous upper bound Ihat>=I gives the lower bound

    L=max(0,B/Ihat), when Ihat>0.                         (4)

If Ihat=0 and B<=0 take L=0. Ihat=0 with B>0 proves inconsistency;
the present pilot consists only of feasible problems and does not implement
a general infeasibility-certificate interface.

The standard continuous norm dual is

    sup_lambda B(lambda), subject to I(lambda)<=1.       (5)

For completeness, equality of its supremum and the primal value follows
in finite-dimensional measurement space. Let Z be the compact symmetric
image of the unit input ball under the kernel map, and E the compact box
of errors. The extended value f(b)=inf{J>=0: b in JZ+E} is proper convex
and closed: its sublevels are cZ+E, hence compact. Its conjugate is
support_E(lambda) if support_Z(lambda)<=1 and +infinity otherwise.
Here support_Z=I and support_E=sum epsilon_i|lambda_i|. Biconjugacy gives
(5). This statement does not require claiming finite dual attainment or
convergence of the numerical SLSQP implementation. The certificate (4)
uses only weak duality and remains sound for unsuccessful solver proposals.

On each row-time segment, k_a is cubic. Both methods obtain numerical
proposals without approximating the control by a fixed time grid: the SDP
uses (1), while the baseline integrates cubic absolute values using real
roots and analytic antiderivatives. Both then use the same certificate path.

For a cubic with Bernstein coefficients b_0,...,b_3 on [l,r],

    integral_l^r |p(t)|dt <= (r-l)/4 * sum_j |b_j|.       (6)

Subdivision tightens (6), but its validity is independent of refinement
termination. The producer uses power-to-Bernstein conversion; the independent
checker recomputes the four coefficients from endpoint values and slopes.
It also checks kernel reconstruction and complete interval coverage.

A primal certificate is a rational piecewise-constant input. Its exact
peak U and exact waypoint values are verified with rational arithmetic and
closed-form convolution, independent of floating solvers. Together these
certify L<=J_star<=U. The stopping condition is

    U-L <= 10^(-6) * max(1,U).                            (7)

Both proposals use the identical inner-LP recovery: initial quarter-segment
mesh, their own cubic roots, up to four iterations and 257 knots. Refined
meshes are inner feasible families, not the exact outer SDP. Floating
controls are rationalized, selected active rows solved exactly, and ALL
rows rechecked. U is recomputed after correction. This recovery is not
proved complete, particularly for degenerate exact coupled measurements.
Failures are retained, not interpreted as infeasibility of the continuous
problem. The tested pipeline is proposal plus shared recovery, not bare SDP.

## 5. Analytic calibration, not a new theorem

For one chain on [0,T], zero initial state, terminal state (d,0,0,0),

    J_star=384|d|/T^4.                                  (8)

At T=d=1, use amplitude 384 and alternating signs +,-,+,- with cuts
(1-1/sqrt(2))/2, 1/2, (1+1/sqrt(2))/2. Exact integration gives (1,0,0,0).
For the four terminal rows, lambda=(384,-192,40,-4) gives the kernel

    p(t)=-4(2t-1)(8t^2-8t+1).

Its signs are those of the proposed input and integral |p|=1, giving a
matching dual lower bound 384. Time/amplitude scaling proves (8), with
the zero case immediate. The benchmark optimum is used only to audit
results; it is not passed to either solver or recovery routine.

## 6. What this establishes and what it does not

This gives a precise synthesis consequence of the candidate K4 representation:
an exact, fixed-size-per-segment affine SDP for unrestricted continuous inputs
and finite linear waypoint constraints. Homogenization, concatenation, duality,
perfect splines and numerical certification are supporting tools, not five
additional innovations. Priority of the concrete lift remains unestablished.
Neither this theorem nor the small deterministic pilot establishes a universal
speed advantage, arbitrary numerical reliability, or a safe flight controller.
