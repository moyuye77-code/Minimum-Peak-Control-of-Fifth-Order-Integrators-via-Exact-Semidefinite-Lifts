# Moment-lift block obstruction

`THEORY.md` states and proves a lower bound `sxdeg(K_d) >= ceil(d/2)` using
standard conic duality and Averkov's polynomial extension-degree theorem.
It excludes finite exact SOC lifts at order five and above, but not general
SDP lifts. Novelty/priority and standalone paper significance are not cleared.

`exact.py`, `verify.py`, and the tests check finite exact algebra witnesses and
the state/moment map. These checks are not a numerical proof of impossibility.
No optimizer or new closed-loop experiment is run. Older archives are untouched.

From the project root:

    python -m pytest research/moment_lift_obstruction -q
    python -m research.moment_lift_obstruction.verify --create
    python -m research.moment_lift_obstruction.verify

`--create` refuses to overwrite an existing archive. The record contains exact
rational coefficients and evaluations for ten orders, plus source hashes.
