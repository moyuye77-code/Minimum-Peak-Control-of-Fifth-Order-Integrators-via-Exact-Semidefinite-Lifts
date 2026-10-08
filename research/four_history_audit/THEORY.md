# Conditional four-integrator certificates: correctness, not a new main theory

2026-09-28. This module tests the computational consequences of the joint K4
lift. The new main candidate remains that lift and the separate representation
lower bound, not the standard duality, Bernstein bounds or recovery below.

## Model and exact finite SDP

Let x=(p,v,a,j), x'=(v,a,j,u), |u(s)|<=J almost everywhere. The initial state
lies in the specified closed box P=[c-r,c+r]. Arrived observation i is

    |w_i^T x(t_i)-y_i| <= epsilon_i,
    0<=t_i<=received_i<=now<=T.

All data in the pilot are rational. No future or unreceived observation can
enter either solver or the independent verifier. The quantity sought is
max w^T x(T) over every bounded measurable input and compatible initial state.

Put a K4 block on each interval between distinct observation timestamps,
0 and T, use the exact affine state dynamics, and impose the observation
strips. Necessity follows by taking the four moments of an actual input on
each interval. For sufficiency, realize each feasible K4 point as a density
on its interval using the established boundary/mixing construction. Concatenate
the inputs. The affine dynamics and observations are then satisfied exactly.
This establishes an exact finite SDP, including zero-error observations and
empty histories. It does not establish that floating-point infeasibility
statuses are valid mathematical certificates.

Joining lifted sets by affine equations and intersections is standard prior
machinery, not an additional novel framework. Continuous state constraints
between observation knots are NOT imposed or certified here.

## Continuous dual used as the strong comparator

Let b(t,w)_k=sum_(r=0)^k w_r*t^(k-r)/(k-r)! and

    k(t,w;s)=sum_(r=0)^3 w_r*(t-s)^(3-r)/(3-r)! for s<=t,
              0 otherwise.

For arbitrary real lambda_i, every consistent trajectory obeys

    w^T x(T) <= U(lambda)
      = h_P(b(T,w)-sum lambda_i b(t_i,w_i))
        +sum(lambda_i*y_i+epsilon_i*|lambda_i|)
        +J integral_0^T |k(T,w;s)-sum lambda_i k(t_i,w_i;s)| ds.   (1)

This follows by substituting the dynamics, adding observation residuals and
maximizing independently over the initial box and |u|<=J. It is a standard
support/Lagrange bound, not a new information-theoretic law. Under nonempty
feasibility, the infimum of (1) is the support value: the trajectory parameter
set P times the input ball is compact in the product Euclidean/weak-star
topology, the finitely many observation maps are continuous, and the
convex-concave minimax interchange applies to the observation Lagrangian.
There need not be a finite minimizing multiplier on singular histories.

On each timestamp interval the residual kernel is cubic. The comparator
minimizes (1) with SLSQP, splitting signed multipliers into nonnegative parts
and representing uncertain-initial-box absolute values by epigraph variables.
Cubic zeros are found numerically; each sign interval is integrated by its
polynomial antiderivative, with an analytic gradient. There is no uniform
time mesh in this comparator. Its optimizer success flag and floating-point
root integrals are proposals only, not certificates or claims of optimality.

## Independent upper and lower certificates

For any rational multiplier vector, subdivide the residual cubic into
intervals [l,r]. If its degree-three Bernstein coefficients are B0,...,B3,

    integral_l^r |p(s)| ds <= (r-l)/4 * sum |Bi|.               (2)

The Bernstein basis is nonnegative and each basis function integrates to
(r-l)/4, which proves (2). Rational interval endpoints and coefficients make
the whole upper bound rational. All leaves must partition [0,T], must not
straddle observation timestamps, and must contain the correct residual
polynomial. Stopping subdivision at a depth cap preserves validity, although
it need not preserve the requested tightness.

A rational, piecewise-constant input plus rational initial state supplies a
lower bound only after exact checks of the initial box, input bounds and
EVERY observation. An inner LP helps propose such an input. Its mesh does
not restrict the unknown continuous input in the outer problem. Numerical
inner infeasibility is not a certificate that the continuous history is
inconsistent, even when the floating-point solver labels it infeasible.

The independent check.py uses the convolution formula for trajectory states,
where the producer uses sequential state propagation. It recomputes cubic
Bernstein coefficients with endpoint Hermite formulas rather than the
producer's power-basis conversion. It imports no optimizer, NumPy, or
certificate producer. Accepted gaps are rational differences U-L; no PSD
eigenvalue tolerance is used to declare an input physically feasible.

The first pilot used one recovery LP on a mesh from each method's own
multiplier roots. Many flat or almost flat residual kernels do not provide
sufficient switching locations. A shared follow-up recovery adds rational
inner knots and repeats the LP, obtaining new multipliers and an exact
checked input. Both methods receive the same rules and iteration/cell caps,
not the other's result. All additional LP calls and original failures are
retained. The resulting pipeline is HYBRID; its success is not attributed
to a bare SDP or bare SLSQP call.

Some recovery proposals use slightly narrowed noisy strips to help rational
rounding. The accepted trajectory is checked against the original strips,
and upper certificates always use the original noise bounds. No zero-error
strip is narrowed. This is numerical certificate engineering, not a claimed
new control contribution or a general convergence theorem.

## Singular-face certificate and the stress family

If the initial state is known exactly, a zero-error terminal observation
attains the maximum or minimum obtained by a constant input u=+J or -J,
and its nonzero polynomial kernel has one sign on [0,T], the input is forced
to that constant almost everywhere. A nonzero polynomial has finitely many
zeros; equality in the integral bound forces saturation elsewhere. Checking
the resulting trajectory against all other observations and the target
gives an exact singleton support certificate. Positivity in the implemented
check is sufficient Bernstein-coefficient positivity, not an asserted
complete polynomial-positivity algorithm.

This is the classical equality case of a support inequality / exposed-face
reduction. It covers a useful boundary case, not every degenerate history.
The separate pure-arithmetic four_history_saturation_check.py verifies it.

For T=J=1, x(0)=0, terminal position observation

    |p(1)-1/24|<=epsilon, epsilon=ell^4/12, 0<=ell<=1,

the exact maximum of -j(1) is -1+2*ell. To see this, let d=1-u in [0,2].
The position deficit is integral (1-s)^3*d(s)/6. For a fixed mass of d,
the cheapest placement is at the right endpoint. The density d=2 on the
last interval of length ell costs ell^4/12 and attains mass 2*ell.
It corresponds to u=+1 before 1-ell and -1 afterwards. For ell>0 the
support multiplier lambda=-6/ell^3 attains the dual value. At ell=0,
the feasible input is identically +1, but finite multipliers need not
attain the infimum; the exposed-face certificate gives value -1.

The analytic values and multipliers are used only by independent tests
and stress-case ground truth, not supplied to either numerical proposal.
Recognizing an exactly saturated observation uses the observation itself,
not the case label or saved answer.

## Scope of the experiment

Four deterministic signal/noise/prior configurations, three time scales
(T=1/2,1,2), and five linear objectives give 60 regular queries. There are
12 arrived mixed observations per query, with receipt times explicitly
after sampling. Five analytic boundary queries complete the 65-case pilot.
The time-scaled copies are not independent random samples. All regular
truths have |u|<=4/5 and strictly feasible noisy observations.

This is a static conditional-estimation benchmark, not a closed-loop swarm
simulation, request policy, packet-loss study, or proof of collision safety.
Cold single-run timing and post-hoc recovery timings cannot establish a
general speed claim. Rational-certificate completeness for arbitrary
zero-error/algebraic histories is not asserted.
