# Minimum-peak fourth-order-chain pilot

This is a control-synthesis consequence of the candidate exact K4 SDP, not a
new peak-control problem or a completed paper. See THEORY.md and PRIOR_ART.md.
The Chinese result summary is `../PEAK_SYNTHESIS_RESULTS.md`.

## Model and fairness

Known initial state; fixed finite waypoint times and linear rows (including
cross-axis rows); minimize common componentwise input peak. The outer SDP
does not discretize the input. The continuous-root dual baseline also has
no control time mesh. Both use the same rational inner-trajectory recovery
and independent exact upper/lower-bound checker. No analytic optimum,
manufacturing input or other method's multipliers enters either solver.

The initial state and zero-input response are eliminated in exact arithmetic
before conversion to solver floats. No waypoint tolerance is enlarged during
recovery. Exact recomputed input peaks include all rational corrections.

## Prespecified pilot and limits

18 deterministic problems: 6 scaled analytic rest-to-rest cases; 8 waypoint
cases from 1/2 axes, 2 manufactured inputs and exact/noisy intermediate rows;
2 ballistic zero-peak cases; 2 required-total-input scales (1e-4 and 1000).
Time/amplitude copies and noisy/exact copies are not independent random trials.
Inputs used to construct the waypoint data are withheld from the solvers.

SDP: Clarabel, absolute/relative gap and feasibility tolerance 1e-10, 300
iterations. Root-dual: SLSQP, 600 iterations, ftol 1e-11. Recovery: at most
4 LP calls and 257 knots; root-guided refinement from that method's own
multipliers. Its initial quarter-segment mesh is an inner witness family,
not the claim of exact finite outer representation. Bernstein upper bounds
use rational subdivision and remain valid at the depth cap.

Certificate target: U-L <= 1e-6*max(1,U). Times include cold proposal, recovery
and independent checking, measured once in fixed method order. No speed
ranking beyond this run, flight experiment or closed-loop experiment is claimed.

## Reproduction

Run from `E:/dig/flock` with the installed Python environment:

    python -m research.peak_synthesis.experiment
    python -m pytest research/peak_synthesis/test_synthesis.py -q

The first command replays all exact certificates, source hashes and known
calibration optima without re-solving the problems. `--create` was used once
to create the immutable archive; it refuses to overwrite an existing file.
Preserve the Python sources and `results/verification.json`. Any follow-up
engineering or ablation should use a separately versioned module/archive.

The archive stores every proposal status, warning, timing, recovery attempt,
certificate and library version. `check.py` imports only Python rational
arithmetic and elementary math, not the producer or an optimization library.
Full regression output is `results/regression.xml` when available.

Archive SHA256:

    e8bdfa4eaeca997f7388f50899ad2bacf638c0995909576f61dbf1cd22dc75da

This pilot does not replace or repair the separate conditional-estimation
archive; its one near-saturated SDP recovery failure remains preserved.
