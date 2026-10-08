# Upper-representation gate

Two focused checks, not two new claimed paper contributions:

- Every finite truncation of the ordinary moment-domination hierarchy has a
  strict gap, including at a fixed linear objective; exact rational witnesses.
- The four-moment boundary has a constructive closed characterization, and
  each fixed-mass fiber admits six 3D SOC constraints. Full variable-mass K4
  and general-order exact affine SDP lifts remain unresolved here.

Read THEORY.md for proofs, classical-source attribution, and limitations.
No optimization solve or new closed-loop simulation is performed.

    python -m pytest research/moment_upper_gate -q
    python -m research.moment_upper_gate.verify --create
    python -m research.moment_upper_gate.verify

Archive creation refuses to overwrite existing evidence.
