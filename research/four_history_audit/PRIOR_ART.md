# Added primary-source check and contribution boundaries

2026-09-28. This extends, rather than replaces, the priority caveats in
../moment_four_lift/PRIOR_ART.md.

Henrion and Lasserre, *Graph Recovery from Incomplete Moment Information*,
author manuscript dated November 5, 2020:
https://homepages.laas.fr/henrion/papers/momgraphcompletion.pdf .
Read the introduction, Theorems 4.1/4.2 statements, their hierarchy
definitions (4.2)/(4.4), and Section 5's L-moment application including
the indicator-density simplification. The full proofs in Section 8 were
not audited. The stated recovery is asymptotic as the hierarchy order
grows, followed by graph reconstruction. Those read results do not themselves
give the explicit finite affine K4 lift used here. This is a scope distinction,
not a proof that no equivalent formulation exists elsewhere.

The newly implemented observer uses standard affine history assembly and
continuous Lagrange duality. Exact Bernstein bounds, rational primal checking,
and saturation equality cases are verification tools, not new main claims.
The experimental contribution is evidence about the candidate's computational
use and limits; adding these tools does not create additional publishable
innovations by count.

The current priority audit remains incomplete. No new claim is made that
the old reachability question remains open throughout the literature.
Specific equivalence to Markov moment constructions and later SDP
representability work is not resolved by this limited comparison.
