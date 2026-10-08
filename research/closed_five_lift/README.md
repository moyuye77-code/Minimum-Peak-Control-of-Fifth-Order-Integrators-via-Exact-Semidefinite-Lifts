# Direct closed fifth-order lift

The full K5, including singular boundaries, has the implemented affine lift
with 39 variables (five outputs and 34 auxiliaries), 21 PSD blocks of maximum
order four, and two scalar inequalities. See [THEORY.md](THEORY.md).

25 targeted tests passed in 4.07 seconds. Full research regression: 1192
passed in 244.52 seconds, no failures/errors/skips; 19 temporary-directory
cleanup warnings. The frozen 18-case support audit passes its declared
1e-7 diagnostic threshold in every case, but 16 have solver precision warnings.
No new closed-loop experiment. See [the full report](../CLOSED_FIVE_LIFT_RESULTS.md).

Read-only replay of exact witnesses and saved numerical residuals:

```powershell
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python313\python.exe' -m research.closed_five_lift.verify
```

Archive SHA256:
`ad064aa5b9c2e06955f0bc95136fee5d8786177d348babe50747fa9f4c1ae705`.

Source/protocol hashes were frozen before the 18 optimizer calls. Each trial
and archive uses exclusive creation. Do not edit frozen source or overwrite
results. README and external reports are not part of the frozen input set.
Numerical success is not a proof of exact returned-state feasibility or novelty.
