# Fifth-order constrained synthesis: validation of the same representation

2026-09-28. This applies the complete K5 model to a classical minimum-L-infinity
input problem. Neither that objective, multi-segment composition, continuous
duality, nor certificate recovery is claimed as an additional innovation.

## 1. Continuous problem and exact affine SDP

There are A independent fifth-order integrator chains, with known initial
states, and measurable scalar inputs u_a. Minimize

    J = ess sup_(t in [0,T]) max_a |u_a(t)|

subject to |w_i' x(t_i)-b_i|<=epsilon_i at finitely many fixed times.
Rows may couple different axes and derivatives. This is a componentwise input
bound; it is not an instantaneous Euclidean ball, a quadrotor attitude model,
an obstacle constraint, or continuous-time state safety.

Homogenize every constant of the complete closed K5 pencil by J. In each
functional use (J,Y0,Y1,Y2,Y3,14 auxiliary entries); in the upper copy use
(J,J-Y0,J/2-Y1,J/3-Y2,J/4-Y3,14 entries). The prefix has the established
homogeneous K4 pencil with six auxiliaries. Replace the scalar inequalities
by Y4>=ell_minus(pF) and J/5-Y4>=ell_plus(pF). All entries remain affine.

For J>0 division by J gives exactly Y/J in K5. The apex must be checked:
the previously proved K4 apex forces Y0=...=Y3=0 when J=0. In the retained
seven-block spaces the polynomial identity

    pF = N² + Delta*q² + p²/12 + Delta*(mz)²/3 + Delta*m^6/180

expresses ell(pF) as a positive combination of five PSD diagonal entries.
It is therefore nonnegative even at J=0. The lower and upper scalar rows
then imply Y4>=0 and -Y4>=0. Thus the only output at the apex is Y=0;
the all-zero auxiliaries realize it. Auxiliary uniqueness is not needed.
Consequently the homogenized projection is exactly {(J,Jy):J>=0,y in K5}.

Partition time only at the prescribed row times. On a segment [l,r] put
h=r-l and v_a(s)=(u_a(l+hs)+J)/2. Its five moments satisfy

    Y_j=integral_0^1 s^j v_a(s)ds,
    integral_0^1 s^j u_a(l+hs)ds=2Y_j-J/(j+1), j=0,...,4.

Remove the known ballistic initial-state contribution beta_i=b_i-w_i'Phi(t_i)x0.
The exact continuous kernel is

    k_ia(t)=sum_(d=0)^4 w_(i,5a+d)*(t_i-t)^(4-d)/(4-d)!, t<=t_i,

and zero afterward. Expand h*k_ia(l+hs) in s. Using its five coefficients
against the signed moments above gives the exact finite affine waypoint row.
The SDP has 39*A*N+1 variables, 21*A*N PSD blocks of order at most four,
two additional scalar inequalities per axis-segment, and the physical rows.

Every feasible measurable input supplies such moments. Conversely the exact
homogenized cone supplies an input on each axis-segment; concatenation and
convolution satisfy all rows with the same peak. The zero-peak case was
proved separately. Hence the feasible peak values and infima agree, without
input-time discretization. A bounded minimizing sequence has a weak-star
convergent subsequence, so feasible continuous problems attain their optimum.

## 2. Strong continuous baseline and exact certificates

For any multiplier lambda define

    B(lambda)=lambda'beta-sum_i epsilon_i*|lambda_i|,
    I(lambda)=sum_a integral_0^T |sum_i lambda_i*k_ia(t)|dt.

The standard continuous dual maximizes B subject to I<=1. The finite
measurement image of the unit input ball is compact and convex. As in the
previous fourth-order derivation, the closed gauge with measurement error
box has this conjugate, proving equality of suprema without assuming dual
attainment or numerical convergence. Its numerical implementation computes
quartic roots and polynomial antiderivatives; it is not a time-grid outer LP.

Certificates require only weak duality. Rational lambda and an exact upper
bound Ihat>=I give L=max(0,B/Ihat). Each quartic on [l,r] has Bernstein
coefficients b0,...,b4, and Ihat is the sum of (r-l)*sum|bj|/5 after rational
subdivision. This follows from the nonnegative Bernstein partition of unity.
The producer uses power-to-Bernstein conversion. The independent checker
uses endpoint values and derivatives, including

    b2=p(l)+(r-l)*p'(l)/2+(r-l)^2*p''(l)/12,

and verifies every kernel, leaf partition and rational arithmetic identity.

A piecewise constant rational input gives a feasible upper bound U after
every requested row and its exact tolerance is checked. The checker evaluates
convolution directly; the task generator uses sequential state transitions.
Neither checker imports the optimization producer or floating libraries.

Both proposal methods feed only their multipliers into the same bounded
recovery process: add rounded numerical roots to a rational mesh, solve an
inner feasible-trajectory LP, correct selected row equations rationally,
and recheck ALL rows. Numerical rank selection may fail; exact checking then
rejects that trajectory rather than silently relaxing data. LP duals can
improve the final lower bound. Initial and final dual certificates are saved
separately; the shared recovery's success is not attributed entirely to SDP.
The SDP primal moments are not directly converted into the reported input.

## 3. Input-defined, physically scaled stopping rule

For each row let Ui be its quartic Bernstein upper bound for
sum_a integral_0^(t_i)|k_ia|, using the single full interval [0,t_i]. Set

    S=max_i (|beta_i|-epsilon_i)_+/Ui,

omitting zero kernels with zero numerator. If all numerators vanish, zero
input is feasible and define S=1 for reporting. Otherwise 0<S<=J_opt follows
from each row's own weak bound. The independent checker recomputes S from
the task data, not from a proposed solution or known optimum.

The protocol requires U-L<=1e-5*S. For nonzero cases this is at least as
strict as a 1e-5 relative-to-optimum requirement. A mere valid bracket does
not count as a converged result. Cases with scales 1e-4 and 1e3 are not
silently evaluated using a common absolute tolerance of 1e-5.

## 4. Scope of the experiment

The frozen protocol specifies 15 deterministic tasks and both methods, with
all statuses, warnings, recovery failures and nonconverged brackets retained.
Analytic anchor values are used only for subsequent checking. Their sign-
quartic generating controls are not passed to either optimization or recovery.
Other targets use undisclosed piecewise-constant controls solely to ensure
feasibility; these inputs are not claimed optimal. Scaled copies are not iid
samples. This is constrained open-loop planning, not closed-loop simulation.

The new evidence tests whether the representation can enter a verified
control computation. It cannot, by itself, establish priority, major practical
advantage, stability, recursive feasibility, intersample safety or acceptance
in any journal. Existing negative fourth-order comparisons remain valid.
