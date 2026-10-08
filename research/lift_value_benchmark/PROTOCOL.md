# Matched formulation-value benchmark (declared before solving)

2026-09-28. The preceding SOC gate did not establish minimality. This tests
computational value of the same exact representation, not a new algorithm,
an originality clearance, a real-time claim, or an extra main contribution.

Methods on exactly the same physical problems:

1. Seven-block homogeneous K4 SDP (three order-3, four order-2; six auxiliaries).
2. Exact generic polynomial-coordinate SDP (orders 6/3/3/2/2; twelve auxiliaries).
3. Continuous polynomial-root support dual, solved by SLSQP.

The generic construction is imported algebraically from the audited exact
coordinate lift, not the larger rational lift or an ordinary time-moment
truncation. Homogenization replaces affine constants by J. For J>0 divide
by J; at J=0 the seven-block restrictions force the output to zero. Extra
generic auxiliary recession directions do not change that output conclusion.

Fixed problem suite (14 cases):

- Scale: 1/2 axes x 2/4/8 waypoint segments, T=1, amplitude=1 (6 cases).
- Time conditioning: the 2-axis 4-segment pattern with T=1/4 and T=4 (2).
- Near saturation: 1 axis, 4 segments, input +1 until 1-epsilon then -1,
  epsilon in {1/10,1/100,1/1000,1/10000} (4); exact optimum J=1.
- Amplitude conditioning: 2 axes, 4 segments, amplitudes 1/10000 and 1000 (2).

General cases use a withheld rational three-cell-per-segment generating
input. Intermediate position/velocity rows couple axes; terminal rows fix
the full state. Intermediate tolerances are zero in this suite. Exact
generation witnesses only verify that these are feasible problems; neither
their controls nor their peak is supplied to any solver or recovery method.
Noisy observations were studied in the older pilot, not silently included here.

Three repetitions per case and method, 126 planning queries. Seed 20260928
fixes the case permutation and a cyclic method-order offset. Each method
occupies each position exactly once for each case. These are repeated timings
of deterministic, partly scaled cases, not 126 independent random trials.
All problems are rebuilt cold; no warm starts or retained model factorizations.
Library import and one-time symbolic coefficient-template generation are
outside per-query timings and reported separately. No external competing
benchmark jobs are launched. Single-thread numerical settings are requested.

SDPs: Clarabel, absolute/relative gap and feasibility tolerances 1e-10,
300 iterations, max_threads=1. Both use the same coefficient assembly for
physical constraints. Root dual: SLSQP, 600 iterations, ftol=1e-11, same
root-integral/gradient routine as the frozen pilot. Failures/warnings remain.

All proposals use the identical frozen trajectory-recovery and independent
rational certificate checker: at most 4 LP calls, at most 257 mesh knots.
Only multipliers seed recovery. The SDP primal variables and true generating
controls do not enter recovery. Retain lower-only evidence if no primal is
recovered, and retain valid-but-not-converged certificates.

For amplitude A, use recovery tolerance 1e-6*min(1,A), under the frozen
criterion gap <= tolerance*max(1,U). Also report gap/max(A,U) and the suite
criterion gap <= 1e-6*max(A,U). The recovery stopping rule is no looser than
the suite criterion: small-amplitude cases cannot pass on an absolute 1e-6
gap alone; large-amplitude cases may receive a stricter stopping tolerance.
The same choice is used for all three methods on a case.

Report model construction, canonicalization, solve-call wall time, proposal
total, recovery, final independent check, and full pipeline time separately.
SDP dimensions refer to the actual canonicalized model. Primary comparison is
per-case median full cost among cases where both methods succeed in all
three repetitions, with all other cases listed rather than excluded from
reliability totals. Secondary stage timings cannot substitute for end-to-end
value. No universal speed or reliability claim follows from this finite suite.

Every completed trial is saved immediately in an exclusive file. Interrupted
runs resume only missing trials after source/manifest checks; never restart
a live process or overwrite completed trials. Final verification rechecks all
rational certificates, calibration values and source hashes without solving.
