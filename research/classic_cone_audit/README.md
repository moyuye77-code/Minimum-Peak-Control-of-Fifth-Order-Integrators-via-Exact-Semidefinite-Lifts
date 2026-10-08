# Standard-conic-compiler audit

This package checks which parts of the existing eight-SOC moment lift are
automatic applications of standard convex atoms. It does not optimize a query,
change the legacy solver, or run a new closed loop.

- `THEORY.md`: exact bidirectional elimination, endpoints included, source scope,
  and the decision to demote the eight-cone construction to supporting machinery.
- `compiler.py`: high-level classic inequalities and actual primitive maps from
  the project's existing CVXPY 1.7.5 installation.
- `check.py`: expected integer affine maps and model dimensions, independently
  stated rather than inferred from the producer's counts.
- `results/verification.json`: frozen records for three primitive maps and five
  history sizes, including hashes of participating code and relevant old sources.

Commands from the project root:

    python -m research.classic_cone_audit.verify
    python -m research.classic_cone_audit.verify --recompile
    python -m pytest research/classic_cone_audit -q

The archive is created once with `--create`, which refuses to replace an existing
record. No optimizer calls occur even with `--recompile`. The local exact proof,
not equality of finite test outputs, establishes mathematical equivalence.
