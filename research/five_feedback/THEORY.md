# Scope of the fifth-order feedback validation

2026-09-28. Uses the frozen exact fifth-order synthesis formulation. This is
application evidence for the same representation theorem, not a new MPC,
delay compensation, safety filter, or stability theorem.

The plant is a box-input fifth-order integrator, each axis with state
(position,velocity,acceleration,jerk,snap). Exact state is measured at sample
times. At each available request, solve the same three-second minimum-peak
terminal-rest problem. Scaling a terminal row by (5-d)!/H^(5-d) does not
change its zero equality. A valid rational certificate with peak at most five
is eligible even if its optimality gap is wider than the synthesis target.

Two execution modes deliberately have different information-to-actuation
timelines. Ideal-instant planning ignores all computation time in the plant;
it is a reference condition only. In the one-period mode the queued old
command is executed while a new plan is computed. Predict the next sampling
state using the current exact state and that known, quantized queue. Solve
from the prediction, without access to the coming external disturbance.

An eligible result completing by the next boundary replaces the queue at
that boundary, never earlier. A late result is discarded. Its worker remains
busy until the measured completion time: every sampling request before then
is skipped. Thus an overrun is not assumed to cancel for free. After completion,
new work starts only at the next sampling time. The old plan continues until
its horizon ends and then commands zero. Initial lack of a plan also commands
zero. Neither expired-plan zero nor the stale plan is a certified safe backup.

Measurements are offline perf-counter durations, used as service times in a
discrete-event replay. They are not worst-case execution times, hard-real-time
tests, or reproducible timing guarantees. JSON serialization, instrumentation
and final archive replay are outside the planning timer. State prediction,
fresh construction, candidate solve, common recovery and the final independent
plan check are inside. The simulator does not move physical hardware.

With disturbance d_a constant in one sample and actuator quantum q, the
executed-versus-intended one-step derivative error is bounded exactly by

    (|d_a|+q/2)*delta^(5-j)/(5-j)!, j=0,...,4.

For a delayed request the prediction already uses the quantized old command,
so its next-state discrepancy is bounded by |d_a|*delta^(5-j)/(5-j)!.
These follow directly from the integrator convolution and are local error
bounds only. The nominal terminal equality of a new plan does not transfer
to a disturbed or delayed plant trajectory. No recursive feasibility or
convergence proof follows from these bounds or successful episodes.

An independent replay imports only the rational synthesis checker, standard
fractions/math/statistics, and no producer or optimizer. It reconstructs the
queue, absolute-time prefix, missed-deadline/busy decisions, problem data,
disturbance and state chain. Producer propagation is sequential; checker
propagation is direct convolution. Both compute the exact time integral of
position squared on constant-input intervals using independently arranged
polynomial arithmetic. RMS is over all axes and the full six seconds;
sampled-position RMS is also retained and is not substituted for it.
