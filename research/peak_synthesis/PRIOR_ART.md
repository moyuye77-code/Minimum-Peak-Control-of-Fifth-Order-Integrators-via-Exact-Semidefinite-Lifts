# Scope audit: minimum-peak waypoint synthesis

Checked 2026-09-28. This note records what was actually accessible. No claim
of exhaustive priority clearance follows from a keyword search.

## Classical problem, not a new problem statement

D. E. McClure, *Perfect spline solutions of L-infinity extremal problems by
control methods*, Journal of Approximation Theory 15(3), 226-242, 1975,
[publisher / DOI](https://www.sciencedirect.com/science/article/pii/0021904575901057),
DOI 10.1016/0021-9045(75)90105-7.

The publisher-indexed abstract describes minimizing a maximum norm of a
linear differential operator under interpolation restrictions and perfect-spline
solutions. This is direct prior scope for the minimum-peak derivative problem.
The full publisher page/PDF was not acquired (fetch failed/403); no numbered
theorem or proof is claimed read. It does not establish that our seven-block
K4 lift was already published, but rules out presenting the problem itself,
or a generic bang-bang/perfect-spline principle, as new.

## Recent SDP motion planning is relevant but solves another problem

B. P. Graesdal, A. Amice, P. A. Parrilo and R. Tedrake,
*Semidefinite Relaxations for Collision-Free Motion Planning*,
[arXiv:2606.14063v1](https://arxiv.org/html/2606.14063v1), 2026-06-12.

Read abstract, introduction/related work, Sections III-A--D and IV-A including
Lemma 2 and its proof. These formulate fixed-degree polynomial trajectories,
integral-squared derivative costs and continuous spherical-obstacle avoidance;
fixed-curve avoidance has an SOS representation, while synthesis is nonconvex
and uses an SDP relaxation. Our unrestricted-input, peak-norm, finite-waypoint
convex problem has no such obstacle guarantee. Thus this is an important scope
comparison, not a matching baseline or evidence we outperform collision-free
planners. Later theorems and experiments were not fully audited here.

## Candidate claim and non-claims

The candidate is still the specific exact affine K4 representation from
`../moment_four_lift/THEORY.md`, with its prior algebra audit in
`../four_lift_prior_audit/THEORY.md`. This module shows a control-synthesis
consequence by homogenizing the same representation, including peak zero.

Do not claim novelty for minimum maximum derivative, waypoint interpolation,
independent-axis concatenation, ordinary convex duality, analytic cubic roots,
Bernstein bounds, rational witness checking, or a generic SDP application.
Do not infer duplication solely because established tools derive a result.
The direct published status of this particular lift remains unresolved.
The general PSD-block lower bound is part of the same representability theme,
not an unrelated extra module added to reach a journal label.

## Remaining research gates

1. Obtain and compare direct bounded-density-moment / integrator-reachability
   representations at the theorem and coordinate-transformation level.
2. Test the recovery routine on its own: both numerical routes currently
   benefit from a common adaptive LP. Attribute any gain to the correct part.
3. Broaden conditioning/scaling and coupled-waypoint cases before a reliability
   or complexity claim. The present 18 problems are deterministic and related.
4. If using a flight/control narrative, add a clearly scoped closed-loop study
   and matched information/actuation constraints. Open-loop waypoint feasibility
   cannot support all-time safety claims.
5. Integrate the one defensible representation result and consequences into a
   manuscript only with explicit attribution and unresolved priority noted.
