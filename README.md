# Minimum-Peak Control of Fifth-Order Integrators via Exact Semidefinite Lifts

Reproducibility materials accompanying the manuscript by **Zihan Yang**,
University of Electronic Science and Technology of China.
Contact: **768216265@qq.com**.

## Download

- [Versioned release: 2026-10-09](https://github.com/moyuye77-code/Minimum-Peak-Control-of-Fifth-Order-Integrators-via-Exact-Semidefinite-Lifts/releases/tag/reproducibility-2026-10-09)
- [Download ESM_1.zip](https://github.com/moyuye77-code/Minimum-Peak-Control-of-Fifth-Order-Integrators-via-Exact-Semidefinite-Lifts/releases/download/reproducibility-2026-10-09/ESM_1.zip)
- [Archive checksum](SHA256SUMS.txt)

The same archive is also available as `ESM_1.zip` in this repository. Extract
it to a new directory and read `README.pdf` or `README.md` inside it.
The historical filename is retained for compatibility: the archive is hosted
here, not supplied as a journal-uploaded "Online Resource 1".

The archive contains the matching manuscript PDF and LaTeX source, mathematical
verification scripts, experimental code, frozen numerical records, rational
certificates, and a detailed claim-to-file guide. Its internal manifest lists
726 payload files with sizes, SHA256 checksums, and article identification.
Warnings, failed cases, and an interrupted experiment are retained.

## Verify the extracted archive

From the **extracted archive directory**, run:

```text
python -B -S research/paper/review_checks.py
```

This standard-library-only command checks the inventory, manuscript structure,
reference metadata, and stored rational certificates.

For the extended algebra and regression checks, use the dependencies in
`requirements-numerical.txt`, then run:

```text
python -B research/paper/review_checks.py --numerical
```

These are read-only replays, not fresh optimization experiments or timing
benchmarks. Do not use `-O`, `PYTHONOPTIMIZE`, or experiment-creation flags on
the distributed records. See `research/paper/REPRODUCIBILITY.md` for details.

## Scope

The manuscript presents an explicit exact finite semidefinite lift for the
fifth-order bounded-density moment body and its minimum-peak synthesis
formulation. The construction has 34 auxiliary scalars and 21 PSD blocks of
order at most four. It does not establish real-time superiority or determine
the exact fifth-order semidefinite extension degree (currently bounded by 3
and 4). The recorded SDP route is slower on all 12 jointly accurate comparison
tasks and activates no plans under the recorded half-second deadline.

Public availability is not journal acceptance, peer-review certification, or a
formal proof. The manuscript retains its actual AI-use disclosure. Existing
third-party notices in the archive remain applicable.
