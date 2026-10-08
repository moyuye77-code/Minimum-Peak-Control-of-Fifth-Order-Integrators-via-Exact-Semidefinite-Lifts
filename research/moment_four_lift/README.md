# Joint four-moment exact SDP construction

The internally proved result is in [THEORY.md](THEORY.md); priority limits
are in [PRIOR_ART.md](PRIOR_ART.md). This is not an SOC theorem, a completed
paper, or a communication algorithm.

From E:/dig/flock, with the existing Python runtime:

    python -m pytest research/moment_four_lift/test_lift.py -q
    python -m research.moment_four_lift.verify --create
    python -m research.moment_four_lift.verify --create --numeric
    python -m research.moment_four_lift.verify
    python -m research.moment_four_lift.verify --numeric

Creation refuses to overwrite evidence. The exact archive is rebuilt from
rational witnesses and symbolic polynomial identities; source digests are
checked. Numerical results preserve statuses and warnings from both models.
They compare normalized support values with exact rational integration for
46 polynomial directions, including repeated roots and saturated inputs.
Numerical errors/eigenvalue residuals are diagnostics, not certified bounds.
No closed-loop simulation is performed by this package.
