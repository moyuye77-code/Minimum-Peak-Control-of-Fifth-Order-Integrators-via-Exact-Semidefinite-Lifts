# General-order polynomial-congruence gate

2026-09-28. See [the result report](../GENERAL_CONGRUENCE_GATE_RESULTS.md),
[the proof](THEORY.md), and [source scope](PRIOR_ART.md).

For every unshifted exponential Hankel block H_n with n>=3, no polynomial
lower-triangular congruence with constant nonzero diagonal can make the
matrix function concave on an open region with positive mass. This holds
at arbitrary polynomial degree and allows dependence on all moment variables.
It is NOT a theorem excluding all SDP lifts, rational transformations, or
nontriangular constructions. The proof eliminates a search class; it is not
promoted to a separate paper contribution.

A physical exact midpoint witness gives a Jensen defect 43/6291456. The
all-order rational interior-boundary Schur formula is documented separately
as classical machinery, not an affine lift. Its uniform first-order
positivity certificate is still unproved here.

21 targeted tests passed. Zero optimizers and zero new closed-loop runs.
The full research regression passed 1108 tests in 238.34 seconds, with
zero failures/errors/skips and 15 temporary-directory cleanup warnings.
The archive was written once and independently replayed with SHA256
`d7ca6297931c152e4a59abeda6781fbc9738735afb05a129ebb0d551d3beeb6d`.

Replay from E:\dig\flock without solving:

```powershell
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python313\python.exe' -m research.general_congruence_gate.verify
```

`--write` uses exclusive creation and refuses to overwrite evidence. Frozen
source includes THEORY.md, PRIOR_ART.md and the original inspected source
PDF. Later work must not silently change this archive or previous experiments.
