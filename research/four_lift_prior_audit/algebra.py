"""Derive the candidate from classical transforms and generic moment matrices.

All checks are symbolic or rational. This is not a novelty detector: a matching
algebraic pipeline does not show that its K4 application was published before.
"""
from fractions import Fraction
import sympy as sp

from research.moment_four_lift.exact import matrices as candidate_matrices

s, t, q, r = sp.symbols('s t q r')
S, T, alpha, beta, tau = sp.symbols('S T alpha beta tau')
V2 = (sp.S.One, s, t, s*s, s*t, t*t)
V1 = V2[:3]
EXPONENTS = tuple((degree-j, j) for degree in range(5)
                  for j in range(degree+1))


def exponential_coefficients():
    """Coefficient extraction, not a hard-coded four-moment inequality."""
    w = sp.Symbol('w')
    p = s*w+t*w*w+q*w**3+r*w**4
    # Terms beyond p^4 have no effect on coefficients through w^4.
    series = sp.expand(sum((-1)**(k+1)*p**k/sp.factorial(k)
                           for k in range(1, 5)))
    return tuple(sp.expand(series.coeff(w, j)) for j in range(1, 5))


def classical_and_congruent():
    aa = exponential_coefficients()
    hankel = sp.Matrix([[aa[1], aa[2]], [aa[2], aa[3]]])
    change = sp.Matrix([[1, 0], [s/2, 1]])
    return hankel, change, (change*hankel*change.T).applyfunc(sp.expand)


def curvature_identities():
    _, _, G = classical_and_congruent()
    xi = sp.Matrix([alpha, beta])
    phi = sp.expand(-(xi.T*G*xi)[0])
    D = t-s*s/2
    D0 = T-S*S/2
    vec = sp.Matrix([alpha+beta*s/2, beta])
    factored = vec*vec.T+sp.diag(beta*beta*D/2, 0)
    hessian_residual = (sp.hessian(phi, (s, t))-factored).applyfunc(sp.expand)
    u, v = s-S, t-T
    ref = {s:S, t:T}
    rem = sp.expand(phi-phi.subs(ref)-sum(
        sp.diff(phi, z).subs(ref)*(z-ref[z]) for z in (s, t)))
    homogeneous_sos = (
        (beta*(v+S*u/2+u*u/6)+alpha*u)**2/2
        + beta*beta*(D*u*u/12+D0*u*u/6+u**4/36))
    direction = sp.Matrix([u, v])
    along = factored.subs({s:S+tau*u, t:T+tau*v}, simultaneous=True)
    integrated = sp.integrate(sp.expand(
        (1-tau)*(direction.T*along*direction)[0]), (tau, 0, 1))
    domain_along = D.subs({s:S+tau*u, t:T+tau*v}, simultaneous=True)
    return dict(
        hessian=list(hessian_residual),
        domain_segment=sp.expand(domain_along-
            ((1-tau)*D0+tau*D+tau*(1-tau)*u*u/2)),
        integral_remainder=sp.expand(rem-integrated),
        homogeneous_remainder=sp.expand(rem-homogeneous_sos),
        zero_beta=sp.expand(homogeneous_sos.subs(beta, 0)-alpha*alpha*u*u/2),
    )


def moment_data():
    yy = {ij:sp.Symbol('Y'+str(ij[0])+str(ij[1])) for ij in EXPONENTS}
    yy[(0,0)] = sp.S.One
    yy[(1,0)] = s
    yy[(0,1)] = t
    return yy


def linearize(poly, yy):
    """q,r occur affinely; only coordinate variables s,t are lifted."""
    return sp.expand(sum(coef*yy[ij] for ij, coef in
                         sp.Poly(sp.expand(poly), s, t).terms()))


def generic_matrices(yy=None):
    yy = moment_data() if yy is None else yy
    gram = sp.Matrix([[linearize(a*b, yy) for b in V2] for a in V2])
    domains = (t-s*s/2, s-t-s*s/2)
    localizers = [sp.Matrix([[linearize(d*a*b, yy) for b in V1]
                            for a in V1]) for d in domains]
    _, _, G = classical_and_congruent()
    barred = G.subs({s:1-s, t:sp.Rational(1,2)-t,
                    q:sp.Rational(1,3)-q, r:sp.Rational(1,4)-r},
                   simultaneous=True)
    ends = [mat.applyfunc(lambda x:linearize(x, yy)) for mat in (G, barred)]
    return [gram, *localizers, *ends]


def compression_residuals():
    yy = moment_data()
    generic = generic_matrices(yy)
    aa = tuple(yy[ij] for ij in ((2,0),(3,0),(4,0),(1,1),(2,1),(0,2)))
    candidate = list(map(sp.Matrix, candidate_matrices((s,t,q,r), aa)))
    changes = [sp.Matrix([[1,0,0,0,0,0], [0,1,0,0,0,0],
                           [0,0,0,1,0,0]])]
    for sign in (1,-1):
        changes.append(sp.Matrix([[1,0,0,0,0,0], [0,1,0,0,0,0],
                                  [0,0,1,sp.Rational(sign,6),0,0]]))
    derived = [C*generic[0]*C.T for C in changes]
    select = sp.Matrix([[1,0,0],[0,1,0]])
    derived.extend(select*mat*select.T for mat in generic[1:3])
    derived.extend(generic[3:])
    return [(a-b).applyfunc(sp.expand) for a,b in zip(candidate, derived)]


def atomic_matrices(point):
    m, z, second, third = map(sp.Rational, point)
    yy = {ij:m**ij[0]*z**ij[1] for ij in EXPONENTS}
    return [mat.subs({q:second, r:third}) for mat in generic_matrices(yy)]


def formula_records():
    original, change, transformed = classical_and_congruent()
    yy = moment_data()
    return dict(
        exponential_coefficients=list(map(str, exponential_coefficients())),
        hankel=[[str(x) for x in row] for row in original.tolist()],
        change=[[str(x) for x in row] for row in change.tolist()],
        transformed=[[str(x) for x in row] for row in transformed.tolist()],
        basis=list(map(str,V2)),
        moment_variables={str(ij):str(v) for ij,v in yy.items()},
        generic_matrices=[[[str(x) for x in row] for row in mat.tolist()]
                          for mat in generic_matrices(yy)],
    )
