# Predetermined fifth-order feedback and service-time protocol

2026-09-28. Freeze this document and every implementation source before the
first archived optimization. Preserve existing synthesis and feedback archives.

Four deterministic scenarios: one-axis nominal, two-axis nominal, two-axis
finite pulse, two-axis persistent alternating external input. Initial positions
are 1/10 and -2/25, and all other derivatives are zero. The pulse is
(1/20,-1/25) on slots 3--5. Persistent signs have period four slots, plus on
the first two; the second axis is shifted one slot, amplitudes 1/30 and 1/40.
Disturbance values are requested only after planning and are not planner input.

Each episode has twelve half-second slots (six seconds), horizon three seconds,
command cap five and ties-to-even actuator quantum 1e-9. All five terminal
states per axis target zero. No intermediate obstacle or flight model is added.
The cap constrains commanded input, not commanded-plus-external disturbance.

For each scenario test ideal-instant and one-period-deadline modes with both
the closed-K5 SDP and continuous quartic-root dual: 16 episodes, 192 slots,
at most 192 planning queries. Alternate method order across scenario/mode
pairs. No warm starts, cross-method solution sharing or retries. Both methods
receive their own resulting exact/predicted state, identical public targets,
and the same frozen synthesis proposal and recovery budgets. A wide but valid
and peak-eligible certificate may be actuated; retain its nonconverged flag.

In deadline mode the plan is intended to start exactly one sample after the
request. Include queue prediction and plan checking in measured service time.
Use the integer nanosecond comparison elapsed<=500000000. While planning, the
old queue commands the plant. Overdue results are discarded, but the worker
stays busy for its full service time; skip later samples occurring before
completion. Do not cancel an overrun for free or apply its result retroactively.
If a new result is invalid or over cap, keep the old queue. If there is none,
or it has expired, command zero without claiming safety or terminal feasibility.

This is offline measured-service-time replay. Timing is one machine/run,
not a latency distribution or WCET. Freeze all thresholds now; do not increase
the deadline or select only successful episodes after observing results.
The ideal mode records period exceedances but deliberately ignores them in
plant propagation; do not label that mode real-time.

Save the manifest before all formal solves, each of 192 slot records
immediately, and each complete episode separately. Save all raw proposals,
certificates, warnings, elapsed nanoseconds, queue identities, predictions,
actual actions and rational state transitions. Independently replay the
entire archive, without solving, after completion. Report all outcomes:
valid/accurate/eligible/accepted plans, late discards, busy skips, no-plan and
expired-plan intervals, exact continuous-position RMS, sampled RMS, final
errors in all five derivatives, peak command, service time and optimization
counts. No tolerance failures or unstable-looking trajectories are omitted.

The primary research claim remains finite SDP versus no finite SOCP for K5.
These tests assess its use and limitations; they are not separate innovations,
random reliability samples, novelty clearance or a journal acceptance test.
