# Theory-aligned diagnostics, 2026-10-02

This internal protocol is fixed before the first production solve. It is not
an external preregistration. Existing experiments and sources are not changed.
Outputs use an exclusive new directory; no retry, solver switching, deleted
case, or outcome-dependent tolerance change is permitted in a completed run.

## A. Singular limits and membership diagnostics

Six indicator densities: empty, full, [0,2/5], [3/5,1], [1/4,3/4], and
[0,1/4] union [1/2,3/4]. Mix each with density 1/2 at epsilon in
{1, 1/10, 1/100, 1/10000, 1/10^6, 1/10^8, 1/10^10, 0}: 48 cases.
The repeats at epsilon=1 are intentional, not independent random samples.

Construct each selected functional coordinate from the manuscript's explicit
19-polynomial basis along the exact rational mixing path. At epsilon=0 use
the removable limit. Check both seven-block copies and the seven prefix
blocks using exact rational PSD pivots and scalar/output equalities. Record
the coordinate norm, distance to the limiting coordinates, and denominator.
Raw 1/p is recorded only for positive p: its divergence is not divergence of
the selected coordinates. The saved witnesses are checked against matrix
entries parsed from the manuscript, not merely the solver's returned status.

For each point separately solve the fixed-output lift, minimizing the infinity
norm of its 34 auxiliaries (a diagnostic objective, not a control algorithm).
No exact witness is passed as an optimizer seed. Record all outputs, residuals,
statuses and warnings. A floating diagnostic pass means output error, negative
PSD eigenvalue and negative scalar slack are each <=1e-7, and an optimal or
optimal_inaccurate status. It is NOT an exact feasibility certificate.

Additional membership controls use zero/full/[1/4,3/4] densities and their
exact supporting polynomials -1, 1, -(t-1/4)(t-3/4). For eta=1e-3,1e-6,1e-9,
test the physical inside mixture and y+eta*c/(c.c). The latter violates the
original moment body's support inequality by exactly eta. These 18 tests
explicitly record numerical acceptance of an exactly outside point; no
floating result is treated as a counterexample to an exact theorem.

## B. Model/recovery attribution

Re-solve all 15 existing fifth-order synthesis tasks with one fresh SDP
proposal each. Do not supply known peaks or generating controls to any
optimizer. Recompute waypoint residuals directly by integrating the chain
kernels against the raw returned homogeneous moments; recompute PSD residuals.
The raw diagnostic threshold is 1e-7*max(1,abs(J)) for cone residuals and
1e-7 after row scaling by max(1,abs(adjusted target),row tolerance).
This is a floating diagnostic, not an exact physical-moment certificate.

Run the unchanged recovery twice: with SDP multipliers and with zero multipliers,
alternating the order by task index. Both have four LP iterations, 257 nodes,
and the original exact gap target U-L<=1e-5*S. Save initial L0, final L and U,
every recovery attempt and certificate. Recompute recovered controls' moments
by exact integration and record their distance to the raw SDP moment vector;
different moments are not necessarily an error because optima can be nonunique.
If U>0, check the recovered physical moments' complete lift with exact witnesses.
If U=0, verify zero controls and zero homogeneous moments directly.

The zero-seed route can learn multipliers from its own LPs. It is not a method
without dual information, nor is it a replacement for the archived dual baseline.
Do not ascribe all certificate success to the SDP or claim exclusivity.

## Common settings and limits

CLARABEL: abs/rel gap and feasibility tolerances 1e-10, max_iter=300,
max_threads=1, warm_start=False. One solve per case; exception/status/warning
records remain in the output. No new closed-loop or speed comparison is made.
Record package versions, source hashes, all cases and per-record hashes. A
read-only checker re-evaluates exact witnesses, raw residuals and certificates.
Finite cases do not establish universal boundary completeness, non-SOC
representability, recovery completeness, stability, or real-time feasibility.
