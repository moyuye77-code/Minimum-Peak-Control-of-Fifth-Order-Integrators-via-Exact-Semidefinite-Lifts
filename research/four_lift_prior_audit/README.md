# Four-lift direct-prior audit

Run from E:/dig/flock:

    python -m research.four_lift_prior_audit.verify
    python -m pytest research/four_lift_prior_audit/test_audit.py -q

The first command reconstructs exact formulas and checks the immutable
results/verification.json and source hashes. `--create` is only for initial
creation; it refuses to overwrite an existing archive. It never optimizes.

THEORY.md supplies the mathematical implication and literature reading
scope. Exact coefficient checks verify calculations, not novelty. The
archive's 76 atomic witnesses check implementation and singular cases; they
do not prove the universal equality of feasible sets by sampling.

The original four-state solver, conditional-estimation code, and their
positive and negative archives remain unchanged. No hypothesis is relaxed.
