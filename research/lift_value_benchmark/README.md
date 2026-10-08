# Matched value test of the four-moment lift

2026-09-28. Completed, frozen, not novelty-cleared. See
[the full Chinese result report](../LIFT_VALUE_BENCHMARK_RESULTS.md),
[the pre-solve protocol](PROTOCOL.md), and [baseline/scope proof](THEORY.md).

- 14 deterministic problems, three methods, three order-balanced repetitions:
  126 queries. Repeats are timings, not independent problem samples.
- Valid exact primal/dual certificates: 42/42 for each method.
- Precision target: seven-block 24/42 (8/14 problems), generic polynomial
  24/42 (8/14), root dual 27/42 (9/14). Each problem's three outcomes agree.
- Seven versus generic: six faster out of seven jointly successful problems;
  median paired full-cost ratio 0.9247003. Seven versus root dual: two faster
  out of six; median ratio 1.3753754. No universal performance ranking.
- 417 optimizer calls: 126 outer plus 291 recovery LPs. No new closed loops.
- All trial results, warning/status fields and failed precision targets are
  preserved. Protocol, code, suite and schedule hashes are recorded before
  solving; completed trials are exclusive files and never overwritten.
- Independent replay completed with the same hash. Full research regression:
  1087 passed in 239.20 seconds, zero failures/errors/skips; 14
  temporary-directory cleanup warnings. See results/regression.xml.

Final archive SHA256:
`0ec4b032e4676b0e66eb84d950764ffadbf7ee6fe886d806fcf28b4c8c22049b`.

From the repository root, replay without solving:

```powershell
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python313\python.exe' -m research.lift_value_benchmark.experiment
```

`--run` refuses to overwrite a complete archive. Preserve the frozen sources
and protocol; any subsequent normalized implementation or recovery study
must be a distinct experiment, not a replacement of this one.

This experiment does not establish minimum block size, prior-art clearance,
real-time performance, stability, collision safety or publication readiness.
