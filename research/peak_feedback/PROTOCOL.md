# Feedback validation protocol (specified before the first archived run)

2026-09-28. Tests the existing peak-input synthesis primitive in feedback;
no new controller innovation, robust stability or collision-safety theorem.
Frozen peak_synthesis sources are imported unchanged.

Four deterministic scenarios, two methods, 20 replans per episode: 8 closed
loops and 160 planning queries. Scenarios are one-axis nominal, two-axis
nominal, two-axis with a finite input-disturbance pulse, and two-axis with
persistent alternating disturbance. Initial position is 0.1 for axis one,
-0.08 for axis two; all other initial derivatives vanish. The pulse is
(0.05,-0.04) during steps 4--7 inclusive; persistent disturbances are
period-eight signs with amplitudes 1/30 and 1/40, phase offset two on axis two.
They are not random draws or independent reliability trials.

At each 0.25-second step, the planner receives the exact current state and
asks for zero position/velocity/acceleration/jerk after a fresh two-second
horizon. It receives no future disturbance. Each method replans from its own
resulting state, using the same model, information and solver/recovery budget.
An independently verified feasible plan with peak <=5 is eligible; a wider
optimality gap is recorded but does not invalidate actuation feasibility.
If no eligible certificate is obtained, explicitly command zero for the next
step and count a fallback. This fallback is not claimed safe or terminally
feasible. Do not clip an infeasible plan and call it a certified solution.

Execute the entire plan prefix, including its internal switches, up to 0.25
seconds. Actuator levels are rounded exactly to the nearest multiple of 1e-9
(ties to even), and bounded by 5. Apply the external constant disturbance
over that step. The executed input is therefore different from the certified
nominal plan; do not attribute its exact terminal equality to the plant.
Rational propagation avoids integration-step error and prevents uncontrolled
denominator growth from repeated unquantized rational recovery.

The producer uses sequential polynomial state propagation. The replay checker
uses convolution and independently reconstructs command clipping, disturbances,
waypoint targets, accepted-plan prefixes and the entire feedback state chain.
The one-step deviation from the planned prefix is checked exactly against
(|disturbance_a|+quantum/2)*period^(4-d)/(4-d)! for derivative d. This is a
local execution-error bound, not a global disturbance/stability guarantee.

Report every proposal warning, recovery failure, planning time, fallback,
certificate gap, command peak, sampled-position RMS and final component
errors. No all-time path safety, communication cost, flight validity or
real-time hardware deadline is tested. Timing is offline and fixed method
order, not randomized runtime evidence. Nonunique minimizers can generate
different feedback trajectories even for two exact formulations.

Run commands from E:/dig/flock:

    python -m pytest research/peak_feedback/test_feedback.py -q
    python -m research.peak_feedback.experiment --create
    python -m research.peak_feedback.experiment

Creation refuses to overwrite an archive. Subsequent engineering must preserve
this first run and its source manifest, including unfavorable results.
