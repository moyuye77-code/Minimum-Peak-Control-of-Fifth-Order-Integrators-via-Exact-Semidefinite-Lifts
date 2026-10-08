# Repository files and versioned archive

The code, frozen data, manuscript snapshot, and reproducibility guides are
available directly in this repository. They are the exact unpacked contents
of **ESM_1.zip, version 2026-10-09**.

- Browse `research/` for code, data, and detailed guides.
- Read `README.pdf` for the printable guide.
- Read `output/pdf/representability-review.pdf` for the matching manuscript.
- [Download the same version from Releases](https://github.com/moyuye77-code/Minimum-Peak-Control-of-Fifth-Order-Integrators-via-Exact-Semidefinite-Lifts/releases/tag/reproducibility-2026-10-09).
- [Browse the versioned file snapshot](https://github.com/moyuye77-code/Minimum-Peak-Control-of-Fifth-Order-Integrators-via-Exact-Semidefinite-Lifts/tree/materials-2026-10-09).

If you clone this repository, run the commands in `README.md` from the repository
root; there is no need to extract the ZIP again. If you only download the ZIP,
run them from its extracted root instead. For example:

```text
python -B -S research/paper/review_checks.py
```

The 726 files listed in `bundle-manifest.json` match the released ZIP byte for
byte. `.gitattributes` prevents Git from changing their line endings. This
layout note, Git attributes, the outer ZIP, and its outer checksum are repository
distribution files outside that internal payload inventory.

The original release tag and ZIP are preserved; `materials-2026-10-09` identifies
the additional commit that exposes the identical unpacked files for browsing.
Checksums do not constitute independent proof verification or journal acceptance.
