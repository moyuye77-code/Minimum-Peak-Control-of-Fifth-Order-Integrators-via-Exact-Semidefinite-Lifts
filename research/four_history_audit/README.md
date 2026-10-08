# Four-integrator conditional estimation and certificate audit

The exact K4 model is reused unchanged from moment_four_lift. Read THEORY.md
for the conditional model, certificates, comparator and honest scope.

From E:/dig/flock:

    python -m pytest research/four_history_audit/test_history.py research/test_four_history_repair.py -q
    python -m research.four_history_audit.experiment
    python -m research.four_history_repair

The last two commands recheck frozen archives without rerunning optimizers.
The initial archive retains all 130 method/query outcomes and their failures.
The recovery archive links it by SHA256 and retains every additional attempt.
Creation, with --create, refuses to overwrite an existing result.

The code at ../four_history_repair.py was added after the first archive;
the first-stage source files and archive remain unchanged. Recovery is
supporting implementation work, not an additional claimed innovation.
The remaining near-saturated SDP failure is deliberately not deleted.
