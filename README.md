# Reproducibility archive

Article: Minimum-Peak Control of Fifth-Order Integrators via Exact Semidefinite Lifts

Intended journal: Applied Mathematics & Optimization (no acceptance or publication claim)

Author and corresponding author: Zihan Yang

Affiliation: School of Information and Communication Engineering, University of
Electronic Science and Technology of China, Chengdu, China

Correspondence: 768216265@qq.com

Archive file: ESM_1.zip. Version: 2026-10-09.

Repository: https://github.com/moyuye77-code/Minimum-Peak-Control-of-Fifth-Order-Integrators-via-Exact-Semidefinite-Lifts

Versioned release: https://github.com/moyuye77-code/Minimum-Peak-Control-of-Fifth-Order-Integrators-via-Exact-Semidefinite-Lifts/releases/tag/reproducibility-2026-10-09

These supporting materials are distributed through the public repository, not
as a journal-uploaded supplement. The legacy filename `ESM_1.zip` and metadata
path `research/paper/ONLINE_RESOURCE_1.json` are retained for compatibility;
they do not assert that a journal hosts or has published this archive.

Caption: Source code, numerical records, rational verification scripts
and a reproducibility guide for the exact semidefinite formulation and its
computational studies, including warnings, failed cases and boundary diagnostics.

Start with `README.pdf` for a printable package guide. The inventory
`bundle-manifest.json` records the article and author identification both at
package level and for each listed file. This metadata is kept separate from
the frozen scientific files so that their contents and hashes remain intact.

This is a scientific subset of the research archive for Minimum-Peak Control
of Fifth-Order Integrators via Exact Semidefinite Lifts. The manuscript is at
`output/pdf/representability-review.pdf`; its source and official template
files are under `research/paper/`.

## Start without a solver installation

Extract the ZIP to a NEW directory. From that directory:

```text
python -B -S research/paper/review_checks.py
```

This verifies every file against the bundled SHA256 inventory, checks manuscript
structure and registered reference metadata, and replays historical fifth-order
synthesis certificates and feedback records with rational arithmetic.
No NumPy, SciPy, SymPy, CVXPY or solver is required for this command.
Expected fifth-order results: 30 valid certificates, accuracy 14/15 versus 12/15,
12 jointly accurate pairs with SDP slower in all of them, 16 feedback episodes,
192 slots, and zero SDP activations under the recorded half-second deadline.

## Extended checks with existing numerical libraries

With the versions in `requirements-numerical.txt` available:

```text
python -B research/paper/review_checks.py --numerical
```

The extended checks reconstruct the current manuscript's matrix identities,
boundary witnesses, fourth-order comparisons and recovery ablation, and the
October 2 theory-aligned experiment. They run the targeted regression tests.
The wrapper blocks CVXPY solve and SciPy optimization entry points during
replay. These commands do not regenerate experimental solver outputs or timings.
Run without `-O` or `PYTHONOPTIMIZE` because frozen scientific checkers use assertions.
The wrapper refuses optimized mode and rechecks the file inventory at the end.

The supplementary study includes 48 exact boundary witnesses, 18 membership
diagnostics and 15 recovery-attribution tasks. All three outside points at
violation 1e-9 were numerically accepted; 48 boundary solves warned about
precision. These limitations and the interrupted first attempt are retained.
Consult `research/theory_alignment/README.md` and the protocol for details.

## What these checks do not establish

This is evidence replay and algebra checking, not a formal proof,
universal novelty search, new timing benchmark or clean-install test.
Numerical libraries are not included. Historical solver source files are unchanged;
some contain original local path hints or old research claims. The current
manuscript and its explicit limitations take precedence.

The bibliographic metadata checker compares the eleven references with the saved
DOI-registry responses; it does not assess the mathematical relevance of a
citation. The stored symbolic-check results are cross-checked, not silently
rerun over the originals.

The inventory checks consistency within this distribution, not an independently
timestamped origin. It deliberately omits itself from its file list. No full-text
third-party papers, credentials, Git history, installed packages, private author
worksheets or unrelated early UAV experiments are included. Do not run creation
flags such as `--create`, `--run` or `--write` on the frozen records.
