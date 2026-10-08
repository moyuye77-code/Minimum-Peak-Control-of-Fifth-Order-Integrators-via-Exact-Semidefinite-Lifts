# Fixed initial numerical audit of the direct closed K5 model

2026-09-28. Written before any optimization call in this module.

This is a support-function implementation gate, not a closed-loop or a
publication-level performance study. It does not replace independent paper
priority/importance review or the user's requested reliable simulations.

Use nine specified root families, each with both leading signs, in this order:
empty; (1/2); (1/4,3/4); (1/5,1/2,4/5); (1/5,2/5,3/5,4/5);
(1/100,1/3,2/3,99/100); (1/2,1/2); (1/4,1/4,3/4,3/4);
(49/100,51/100). Pad coefficients to degree four and normalize by their
Euclidean norm. Exactly 18 optimization calls, one fresh non-warm solve per
case. No retry, rescaling, solver switching or case removal in this archive.

The exact optimum is independently the integral of the positive part of the
given polynomial. Known rational roots partition [0,1]; rational integration
supplies both exact support and a feasible bang-bang density. Also evaluate
a separate floating-point polynomial-root oracle from coefficients alone.
The known-root oracle is not a timing competitor; record the numerical root
oracle separately, without claiming a general timing advantage from 18 cases.

CLARABEL settings: abs/rel gap and feasibility tolerances 1e-10, max_iter=300,
max_threads=1, warm_start=False. Record construction and canonicalization once,
and per-case solve-call timing, status, iterations, warnings, returned point,
all 34 auxiliary values, minimum eigenvalue of all 21 pencils, two scalar
slacks and signed support error. Preserve every failure and precision warning.

Diagnostic pass: finite outputs; normalized absolute support error<=1e-7;
all pencil eigenvalues and scalar slacks>=-1e-7; returned objective agrees with
the returned point within 1e-9. This is a numerical acceptance threshold,
not an exact rational certificate for the solver-returned state. Separately
check the exact analytic maximizing point with finite rational lift witnesses.

Write each trial exclusively immediately after its run. Freeze the protocol
and source hashes before dispatch. Read-only verification must recompute
targets and residuals from saved outputs without invoking an optimizer.
