# Fixed constrained fifth-order pilot

2026-09-28. Freeze this file and all implementation sources before dispatch.

15 cases, in model.cases() order: two zero-peak ballistic tasks (one/two axes);
total-input anchors with peak 1e-4 and 1e3; three five-terminal-state quartic
anchors (unit time, half time with peak 3/2, and endpoint-near switches);
eight waypoint tasks from axes 1/2, seeds 1/2, intermediate tolerance 0/1/200.
The latter include quarter-time position/velocity rows, coupled across axes,
and exact full terminal states. Keep all cases; no selection after outcomes.

Methods: the exact homogeneous closed-K5 SDP and the continuous quartic-root
dual. Exactly one proposal run per method/case, 30 proposal calls total.
Even cases run SDP first, odd cases dual first. Each proposal constructs its
own optimization problem. Common coefficient templates may be cached; no
warm optimization, cross-method warm seed, solution sharing or retry is used.
One run per pair cannot support statistical speed claims.

SDP: CLARABEL, abs/rel gap and feasibility tolerance 1e-10, max_iter=300,
max_threads=1, warm_start=False. Continuous dual: SLSQP, analytic objective
gradient and root-based norm gradient, maxiter=600, ftol=1e-11. Raw physical
input amplitudes are preserved. Row normalization is fixed in the task
definition and shared by both methods, not tuned after outcomes.

Shared certificate stage: rational multipliers limited to denominator 1e9;
four recovery-LP iterations maximum; quarter-segment initial mesh; root
locations rounded to denominator 1e10; at most 257 mesh nodes. LP tolerances
1e-10, one thread, exact rational row correction afterward. Bernstein norm
upper enclosure tolerance 1e-12 times relative interval length, depth cap 32.
Record every LP call and warning. Do not rerun with a larger budget in this
archive when a case fails to reach the requested gap.

Convergence is U-L<=1e-5*S, where S is the task-only weak row bound described
in THEORY.md. Save valid but nonconverged certificates distinctly from failures
to recover an exact primal input. Initial-dual and final-dual bounds are both
preserved for attribution. Known peak values and generating input profiles
must not enter either proposal or recovery API.

Record full proposal data, final rational trajectory, dual leaf certificate,
attempt logs, initial bound, timings including recovery and final checking,
and versions. Write each trial exclusively as soon as it finishes. Independent
read-only replay checks every certificate and all hashes without optimization.
No new closed-loop episode is part of this pilot, and no safety claim follows.
