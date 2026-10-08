# Reproducibility archive: guide

Article: Minimum-Peak Control of Fifth-Order Integrators via Exact Semidefinite Lifts

Intended journal: Applied Mathematics & Optimization (no acceptance or publication claim)

Author and corresponding author: Zihan Yang

Affiliation: School of Information and Communication Engineering, University of
Electronic Science and Technology of China, Chengdu, China

Correspondence: 768216265@qq.com

Archive file: ESM_1.zip, version 2026-10-09. The printable entry guide is `README.pdf`.

Repository: https://github.com/moyuye77-code/Minimum-Peak-Control-of-Fifth-Order-Integrators-via-Exact-Semidefinite-Lifts

Versioned release: https://github.com/moyuye77-code/Minimum-Peak-Control-of-Fifth-Order-Integrators-via-Exact-Semidefinite-Lifts/releases/tag/reproducibility-2026-10-09

The materials are hosted in this public repository, not in a journal-uploaded
supplement. Legacy archive and metadata filenames are retained for compatibility.

This guide accompanies Minimum-Peak Control of Fifth-Order Integrators via
Exact Semidefinite Lifts. Paths below are relative to the extracted
reproducibility package. Run commands from that directory.

## 1. Read-only replay

The standard-library entry point verifies the package inventory, manuscript
structure, reference metadata and exact fifth-order certificates:

```text
python -B -S research/paper/review_checks.py
```

It does not require NumPy, SciPy, SymPy, CVXPY or a solver. Do not use Python's
`-O` flag or `PYTHONOPTIMIZE`: the frozen checkers contain assertions.
The wrapper refuses optimized execution and checks the inventory again after
the replay.

To replay the fifth-order synthesis and feedback records alone:

```text
python -B -S research/paper/reproduce.py
```

| Recorded study | Expected result |
|---|---|
| Fifth-order synthesis | 30 valid certificates; accuracy 14/15 SDP and 12/15 dual |
| Jointly accurate synthesis tasks | 12 pairs; SDP faster in 0; median cost ratio about 1.97250724 |
| Fifth-order feedback | 16 episodes, 192 slots, 158 planning queries |
| Ideal mode | 48/48 plans activated for each method |
| Half-second deadline, SDP | 17 queries, 0 activated, 17 late, 31 busy skips |
| Half-second deadline, dual | 45 queries, 42 activated, 3 late, 3 busy skips |

Feasibility and primal/dual certificate checks use rational arithmetic.
Floating tolerances apply only to displayed summary quantities, such as RMS,
not to exact feasibility, schedule decisions or gap classification.

## 2. Symbolic and numerical checks

With the packages in `requirements-numerical.txt` available:

```text
python -B research/paper/review_checks.py --numerical
```

This mode reconstructs the polynomial identities and boundary witnesses,
replays the supporting studies, and runs the regression tests. It blocks
CVXPY solve and SciPy optimization entry points during replay. It does not
generate new solver outputs or timing measurements.

The manuscript coefficient checker can also be run separately:

```text
python -B research/paper/check_lift_coefficients.py
python -B -m unittest research.paper.test_reader_entrypoints -v
python -B -m unittest research.paper.test_review_clarifications -v
```

The checker parses the seven additional pencils in Section 3.2.1 and checks
all 67 displayed matrix entries, the 19-coordinate basis, the joint remainder
identity and the zero-peak identity. The regression tests include deliberately
altered coefficients, a missing matrix, a changed certificate and an incorrectly
activated late plan; each alteration must be rejected.

These calculations do not replace the analytic boundary-existence proofs or
the universal non-SOC lower-bound argument.

The additional symbolic script and its recorded outputs are under
`research/math_derivation_audit/20261002/received/` and `rerun/`.
Its 49 named algebra checks are compared with the recorded rerun. That script
transcribes an earlier manuscript rather than parsing the current source;
the current-source checks above provide the connection to the displayed formulas.
If executing the additional script, copy it to a new directory first: it writes
`symbolic_audit_results.json` alongside itself.

## 3. Claim-to-evidence map

| Manuscript component | Supporting files |
|---|---|
| Four-moment prefix, Section 3.1 | `research/moment_four_lift/` |
| Fifth-order weighted identity | `research/schur_remainder_recursion/` |
| Finite boundary witnesses and 21-block model | `research/closed_five_lift/` |
| Necessary block bound, Section 3.3 | Written proof and cited prior results; `research/moment_lift_obstruction/` checks finite examples only |
| Minimum-peak equivalence and zero apex | Section 3.4 and `research/five_synthesis/THEORY.md` |
| Fifth-order synthesis, Table 1 | `research/five_synthesis/results/` |
| Singular-limit and membership diagnostics | `research/theory_alignment/results/run-20261002-r2/` |
| Fifth-order recovery attribution, Table 2 | `research/theory_alignment/results/run-20261002-r2/synthesis/` |
| Fourth-order matched comparison, Table 3 | `research/lift_value_benchmark/` |
| Fourth-order recovery ablation, Table 4 | `research/peak_synthesis/` and `research/peak_recovery_results/` |
| Fifth-order feedback, Tables 5-6 | `research/five_feedback/results/` |
| Fourth-order instantaneous feedback, Table 7 | `research/peak_feedback/` |

The bibliography checker compares eleven formal reference entries against saved
Crossref responses in `research/paper/reference_metadata/`:

```text
python -B -S research/paper/check_references.py
```

This checks bibliographic fields, not whether each citation establishes a
mathematical claim.

## 4. Re-solving versus replaying

The completed supplementary run contains 48 boundary-path solves, 18
inside/outside diagnostics, 15 synthesis proposals and 67 recovery LPs:
81 SDP calls in total. Its separate interrupted first attempt is retained.
See `research/theory_alignment/README.md` and its `PROTOCOL.md` for the
precise task definitions, warnings, tolerances and limitations.

The following commands require the packages in `requirements-numerical.txt`,
including CVXPY and Clarabel, even though replay does not call an optimizer.
Install these dependencies before running the commands; the verification
modules import the numerical model definitions.

```text
python -B -m research.theory_alignment.verify --output research/theory_alignment/results/run-20261002-r2
python -B -m unittest research.theory_alignment.test_alignment -v
```

The completed study includes three exactly outside points accepted at floating
tolerance and one SDP-seeded recovery target failure. Neither is removed.

For a new supplementary solve, use a new, nonexistent output directory:

```text
python -B -m research.theory_alignment.experiment --output NEW_EMPTY_RUN_DIRECTORY
```

Do not run creation or writing flags on the distributed frozen records.
The historical experiment interfaces use write-once creation modes; new
experiments require separate working directories.

Recorded numerical versions are Python 3.13, NumPy 2.3.5, SciPy 1.17.1,
SymPy 1.14.0, CVXPY 1.7.5 and Clarabel 0.11.1 on Windows.
`requirements-numerical.txt` records the numerical package versions;
it is not a transitive lockfile. Libraries are not included.

Stored service times can be replayed as scheduling inputs. Fresh solve times
depend on hardware and software, so the deadline results concern the recorded
implementation, not all possible implementations. Hashes check consistency
within this distribution, not an external timestamp.

## 5. Document source

The source `research/paper/representability.tex` contains the text, proofs,
tables and bibliography. Keep it with `sn-jnl.cls` and `sn-mathphys-num.bst`.
There are no external figure inputs and no BibTeX pass is needed.
The structure checker validates source layout and references; it does not
compile the PDF or establish mathematical correctness.
