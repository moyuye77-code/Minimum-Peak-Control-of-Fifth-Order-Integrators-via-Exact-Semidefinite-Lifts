# Certified fifth-order constrained synthesis

This frozen pilot applies the complete K5 representation to a classical
continuous minimum-peak-input problem. It is not a new control objective,
closed-loop validation, or an additional claimed innovation.

Read `THEORY.md`, `PROTOCOL.md`, and `../FIVE_SYNTHESIS_RESULTS.md`.
The authoritative data are the pre-solve `results/manifest.json`, 30 individual
`results/trials/*.json` records and `results/verification.json`.

Read-only replay from the repository root:

```powershell
python -m research.five_synthesis.experiment
```

Do not use `--create` on the existing archive. Do not edit the frozen Python
sources, THEORY or PROTOCOL. A changed experiment belongs in a new archive.

Both methods produced 15/15 valid rational brackets; 14/15 SDP and 12/15
continuous-dual runs met the input-scaled gap target. Thirteen SDP precision
warnings and three non-success dual statuses remain recorded. SDP was slower
on all 12 commonly converged tasks. The shared recovery reads multipliers,
not SDP primal moments or known generating trajectories. Initial and final
lower bounds are separate to retain attribution.
