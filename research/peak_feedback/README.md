# Feedback use of the frozen peak-input synthesis primitive

See [pre-run protocol](PROTOCOL.md) and [complete results and limitations](../PEAK_FEEDBACK_RESULTS.md).
This directory contains the first archived deterministic feedback experiment;
it is not a new controller or a real-time/stability/safety certificate.

- `protocol.py`: predetermined scenarios, state-dependent planning queries,
  quantized prefix actuation and rational plant propagation.
- `experiment.py`: two frozen proposal/recovery pipelines, exclusive archive
  creation, source manifest and replay entry point.
- `check.py`: separate exact-convolution replay of the causal state chain and
  nominal certificates, command rounding, disturbances and local error bounds.
- `test_feedback.py`: twelve unit checks made before the first archive run.
- `results/verification.json`: eight episodes, 160 queries, 337 optimizer calls.
- `results/regression.xml`: full research regression, not numerical experiment data.

Replay from E:/dig/flock:

    python -m research.peak_feedback.experiment

The `--create` option refuses to overwrite the archive. Keep all Python sources
and frozen parent dependencies unchanged; new experiments require a new sidecar.
Documentation updates do not change the archived algorithms.

All 160 queries produced eligible plans and target-gap certificates, but 19/80
SDP planning times exceeded the simulated 0.25-second period. Computation delay
was not modeled in the plant. This negative result is part of the experiment,
not grounds to delete cases or call the loop real-time. Both methods have
nonzero final errors; only commands, not commands plus exogenous disturbance,
are bounded by the declared actuator cap.
