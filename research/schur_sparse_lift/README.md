# Sparse spaces for the same fifth-order category separation

See [THEORY.md](THEORY.md) for the exact open lift and the block-preserving
closure proof. They give **3 <= sxdeg(K5) <= 4** for the full closed body.
The open lift has 50 shared monomial coefficients and nine PSD blocks of
orders 2,2,2,4,4,3,3,3,2. Its affine coefficient model is built exactly in
algebra.py. The final closed model's twice-polar coefficient implementation,
block count and auxiliary count are not supplied.

14 targeted tests passed in 1.47 seconds. Source/dependency archive creation
and independent replay passed. Full research regression: 1167 passed in
244.71 seconds, no failures/errors/skips, 18 temporary-directory cleanup
warnings. No optimizer or new closed loop.

Frozen archive SHA256:
`a26d9bf78e9995871055df2c5237270d8a881f3cd59b3cf54453bc052486860e`.

Read-only replay:

```powershell
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python313\python.exe' -m research.schur_sparse_lift.verify
```

Do not edit frozen inputs or overwrite verification.json. Follow-up reports
and README are outside that frozen source set. The generic sparse-space and
affine-polar methods remain prior art; this is a refinement of the one K5
candidate, not another paper or a cleared originality claim.
