# Domain-weighted Schur remainder reduction

The all-order block/translation argument reduces a restricted determinant-module
certificate question to a lower-denominator rational family. It is not a closed
all-order induction. The first nontrivial base has an exact weighted certificate,
an explicit open-epigraph affine SDP lift, and a separate full K5 closure proof.
With the existing lower bound this gives finite SDP versus no finite SOCP for K5.
Priority, significance and numerical usefulness have not been cleared.

See [THEORY.md](THEORY.md) and [PRIOR_ART.md](PRIOR_ART.md).
27 targeted tests passed in 4.68 seconds. The full research suite passed 1153
tests in 242.00 seconds (0 failures/errors/skips; 17 temporary-directory cleanup
warnings). No new optimization or closed-loop experiments. Sources and
dependencies are frozen in verification.json; independent replay passed.

Archive SHA256:
`56284ca0969752f1cb1cd12ac616e6e5b85e6519cb752605b38017a18c308e66`.
See [the Chinese result report](../SCHUR_REMAINDER_RECURSION_RESULTS.md) and
[the scoped direct comparison](DIRECT_COMPARISON.md) for current adjudication.

Replay, without optimizing:

```powershell
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python313\python.exe' -m research.schur_remainder_recursion.verify
```

The --write option creates the archive exclusively and refuses to overwrite it.
Do not edit frozen inputs after archiving. The large open-lift blocks are not
claimed to be minimal or to bound the final closed K5 representation.
