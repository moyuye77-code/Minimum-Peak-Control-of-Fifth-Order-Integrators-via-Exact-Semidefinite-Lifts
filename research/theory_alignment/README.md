# Theory-aligned supplementary experiments

Completed run: `results/run-20261002-r2/`. Date: 2026-10-02.
Internal protocol fixed before solving; not external preregistration.
Historical sources, records, tasks and timing results are unchanged.

## What is tested

| Claim | New evidence | Limit |
|---|---|---|
| Finite selected coordinates at singular boundaries | 48 exact rational witnesses, parsed manuscript matrices, six paths down to epsilon=0 | Specific paths, not the universal analytic proof |
| Numerical membership versus exact membership | 9 inside and 9 exactly separated outside points | Floating residual acceptance is not exact feasibility |
| Raw homogeneous SDP respects original task rows numerically | 15 fresh proposals, direct binomial convolution and cone residuals | Does not certify exact physical raw moments |
| Contribution of proposal multipliers to shared recovery | Paired SDP-seeded and zero-seeded recovery on all 15 tasks | Fixed-budget comparison, not exclusive SDP certification |
| Recovered physical input maps into the complete lift; zero apex | 53 nonzero-peak segment witnesses plus 2 zero tasks per route | Forward direction only, not inversion of arbitrary SDP moments |

## Completed results

- Boundary paths: 48/48 exact witness checks, 48/48 floating diagnostic passes,
  **48/48 precision warnings**. Selected functional coordinates have infinity
  norm at most one. Maximum distance to the endpoint at epsilon=1e-10 is
  less than 1.5e-10. Values and limits are stored exactly as rational strings.
- Inside/outside tests: inside diagnostic acceptance 9/9. All six outside
  points at support violation 1e-3 or 1e-6 are reported infeasible. **All three
  outside points at violation 1e-9 pass the 1e-7 numerical diagnostic** and
  return optimal_inaccurate. This is an intentional tolerance-boundary stress
  test, not a claim that an exact SDP lift contains those points.
- Raw synthesis: 15/15 meet declared scaled residual checks; 13 warnings.
  Maximum scaled original-row violation is 2.1849068110313397e-9.
- SDP-seeded recovery: 15/15 valid brackets, 14/15 target gaps, 19 recovery LPs.
  Initial L0 suffices with the final U in all 14 successes; no lower-bound
  improvement in recovery. The unsuccessful task remains waypoints_a2_s1_e0.
- Zero-seeded recovery: 15/15 valid brackets, 4/15 target gaps, 48 LPs.
  Its LPs improve the initial lower bound in 13 tasks; it is not a dual-free
  algorithm. Neither route receives known optima or generating controls.
- Exact recovered moments: 53 nonzero-peak segment witness checks and two
  zero-apex checks for each route. The maximum raw-to-recovered moment
  discrepancy for SDP-seeded recovery is 0.16249219485332497; matching moments
  is not imposed, so there is no direct raw-moment inversion claim.
- Completed-run cost: **81 fresh SDP calls plus 67 recovery LP calls**.
  No fresh dual benchmark, timing-superiority test or closed-loop experiment.

The separate first attempt in `results/run-20261002/` stopped when serializing
a NumPy boolean in the first synthesis record. It retains all 66 membership
records, its manifest, a partial synthesis record, and original local sources.
The correction only casts the boolean and pre-serializes JSON. The second
run repeats the identical protocol; neither thresholds nor cases were selected
after seeing outcomes. The first attempt incurred 67 additional SDP calls
and 2 recovery LP calls. Do not aggregate it with the completed run as another
independent statistical trial. See its `INTERRUPTED.md`.

## Reproduce or verify

Use the recorded numerical dependencies from the repository's
`requirements-numerical.txt`. Versions were NumPy 2.3.5, SciPy 1.17.1,
SymPy 1.14.0, CVXPY 1.7.5, Clarabel 0.11.1, Python 3.13 on Windows.
This run used already installed dependencies from the original research
workspace; a clean install was not performed and no package was installed.

From the repository root, read-only replay (no optimizer calls):

```text
python -B -m research.theory_alignment.verify --output research/theory_alignment/results/run-20261002-r2
python -B -m unittest research.theory_alignment.test_alignment -v
```

Fresh re-solving must use a NEW output directory:

```text
python -B -m research.theory_alignment.experiment --output NEW_EMPTY_RUN_DIRECTORY
```

The directory must not exist. No resume/overwrite option is provided. All
completed-run JSON records and the source/protocol manifest have SHA256 hashes;
the checker rejects changes. It also disables optimizer calls during replay.
The coefficient matrices are parsed from the current manuscript and checked
against exact witnesses. The recorded initial manuscript hash predates the
addition of these experimental results; the mathematical matrix entries did
not change. These checks establish reproducibility of the recorded calculations,
not the universal analytic claims or novelty.
