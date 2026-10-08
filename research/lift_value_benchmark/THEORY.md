# Scope of the matched formulation benchmark

2026-09-28. This note explains the baseline and interpretation; it does not
introduce a new optimization method. The pre-solve choices are frozen in
PROTOCOL.md and results/manifest.json. This note was written after the run.

## Exact generic homogeneous baseline

Let Lambda act on the fifteen coordinate monomials s^i t^j, i+j<=4, with
Lambda(1)=1, Lambda(s)=m and Lambda(t)=z. The remaining twelve functional
values are auxiliary scalars; q,r are separate output coordinates. The
audited polynomial-coordinate lift has blocks

    Lambda(V2 V2'), Lambda(D V1 V1'), Lambda(E V1 V1'),
    Lambda(G), Lambda(Gbar),

where V2=(1,s,t,s^2,st,t^2), V1=(1,s,t), D=t-s^2/2,
E=s-t-s^2/2. Their orders are 6,3,3,2,2. This is NOT the non-exact finite
truncation of a time-moment domination hierarchy. Its output projection is
exactly K4, including boundary points; see ../four_lift_prior_audit/THEORY.md.

Replace Lambda(1)=1 by Lambda(1)=J and replace all other output/auxiliary
scalars by their homogeneous values. For J>0, dividing all entries by J
reduces to the known exact affine lift. For J=0, restricting the generic
blocks to the seven small Gram/localizing spaces gives precisely the
homogeneous seven-block system. Its apex argument forces all four output
coordinates, and its six retained auxiliaries, to zero. The six additional
generic auxiliaries need not all vanish, but cannot create nonzero output
directions. Conversely J=0 and all variables zero is feasible. Thus both
project onto the same closed cone {(J,J*y): J>=0, y in K4}.

The exact polynomial restriction identities are tested with symbolic J;
this is not a numerical argument based on near-zero solver output. The
rest of the synthesis assembly is identical for the two SDPs. With A axes
and N segments, their actual scalar-variable counts are 1+10AN versus
1+16AN. In the largest case A=2,N=8, canonicalization produced 161 versus
257 variables and 517 versus 661 matrix rows, with 36 equality rows in
both. These dimensions alone do not imply a wall-time advantage.

## Withheld generation and a known-optimum family

The general suite fixes exact waypoint data using rational piecewise-constant
inputs. Only the Problem object is passed to any proposal or recovery method;
the generating mesh, controls and peak are not solver inputs. Independent
convolution checks verify the generation, not optimality of those controls.

For the saturated family, T=1 and u=1 before a=1-epsilon, u=-1 afterwards.
The prescribed terminal acceleration and jerk fix the pairing

    integral_0^1 (a-t) u(t) dt = integral_0^1 |a-t| dt > 0.

Any competing input with peak J must have this pairing <=
J*integral_0^1 |a-t| dt. Consequently J>=1, and the generating input attains
1. The two terminal quantities suffice for this lower bound; additional
waypoints cannot invalidate it. This calibration is classical support
duality, not an innovation. Every archived certificate encloses J*=1.

## What the experiment measures

All three methods supply multipliers to the same frozen recovery process;
the SDP primal moments and objectives are not used. Its certificate combines
an exact feasible control with a rational continuous dual lower bound. The
protocol caps recovery at four LPs and 257 knots, so success means success of
that particular hybrid pipeline within this budget, not of every possible
implementation of a mathematical formulation.

For generating amplitude A, report g=(U-L)/max(A,U), and accept g<=1e-6.
The frozen recovery API uses tolerance=1e-6*min(1,A) with its old
gap<=tolerance*max(1,U) stopping rule. Since

    min(1,A)*max(1,U) <= max(A,U),

that rule is never looser. At A=1000, generic and root-dual outputs meet the
suite criterion but not the stricter recovery stopping rule. Both status
fields are archived. The same tolerances apply to all methods per case.
The small-amplitude case does not receive an absolute-1e-6 shortcut.

Cold-model timings separate assembly, canonicalization, solve call,
recovery and independent final checking. Solve-call time includes the
modeling interface's cached solve dispatch; it is not a solver-kernel-only
measurement. Componentwise medians do not add to the median total.
One-time imports and the 0.09637-second symbolic template initialization
are excluded from per-query totals. No persistent/warm-start deployment
comparison was performed. Raw physical rows, not dimensionless normalized
problems, are used by all methods; conditioning failures do not show that
standard rescaling could not resolve them.

## Post-run failure localization, without new solves

All 291 recovery LPs returned status zero and every attempted rational
primal correction passed. All 126 final certificates are valid. The failures
are failures to close the gap, not failure to obtain a feasible trajectory.
Some runs exhaust four LPs; others encounter the next-mesh knot cap earlier.
Returned best trajectories may come from an earlier iteration, so their
knot counts alone must not be used to diagnose the stopping cause.

For the SDP failures with epsilon=1/100,1/1000,1/10000, the lower bounds
are already within 3e-9 of the known optimum 1. The excess certified upper
bound dominates the gap. This locates the failure in primal recovery within
the budget, not in a disproof of the exact outer lift. It does not establish
which recovery change would solve the issue without further experiments.
For two axes/eight segments, both bounds and recovery remain relevant;
the root-dual outer solve reaches its iteration limit. The generic warning
labels or a root solver's successful exit cannot substitute for certificates.

There is no broad accuracy/runtime dominance in this suite. The seven-block
compression has a modest matched cost benefit over the generic polynomial
model, and complementary success cases, but this alone does not establish
the significance or priority required of the intended research paper.
