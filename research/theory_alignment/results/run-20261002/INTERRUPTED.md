# Interrupted implementation run, retained without overwriting

2026-10-02. All 48 boundary and 18 membership records were saved. The first
synthesis case was solved and recovered, but serializing its diagnostic failed
because a NumPy boolean was not JSON serializable. The partial synthesis/00.json
is NOT a complete scientific record. There is no final summary or completed
verification for this directory. Original local sources are in interrupted-source.

Correction: cast that diagnostic to Python bool and serialize before opening
the output file; add JSON serialization to the regression test. No mathematical
formula, solver setting, case, target threshold or recovery budget was changed.

The complete protocol is rerun once in run-20261002-r2. All original files are
retained, including this failed record. The first attempt's 67 SDP calls and
the two first-case recovery calls are additional to the completed run counts;
these are not silently included in a successful-run timing comparison.
