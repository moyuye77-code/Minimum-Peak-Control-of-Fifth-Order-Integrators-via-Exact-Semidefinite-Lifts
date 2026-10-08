import pytest
import sympy as s
from .algebra import (
    identities, certificate_squares, fifth_schur, h0, X, U, AT, DELTA,
    DELTA0, P, P0, m, z, q, r, M, Z, Q, R, c, rank_reduction,
    exponential_coefficients,
)


@pytest.mark.parametrize('identity', (
    'schur_split', 'domain_invariance', 'numerator_translation',
    'boundary_translation', 'reduced_certificate', 'full_joint_certificate',
))
def test_exact_joint_identity(identity):
    assert identities()[identity] == 0


def test_certificate_basis_and_degree():
    for coefficient, gx, gu, f in certificate_squares():
        assert coefficient > 0
        assert gx in (1, m, DELTA, m*DELTA)
        assert gu in (1, M, DELTA0, M*DELTA0)
        assert not s.denom(f).free_symbols
        assert s.Poly(gx*f*f, X).total_degree() <= 10
    assert [s.binomial(4+d,d) for d in (5,4,3,2)] == [126,70,35,15]
    assert s.binomial(4+10,10) == 1001
    assert s.Poly(s.cancel(P*fifth_schur()[2]), X).total_degree() <= 10


@pytest.mark.parametrize('n', (1,2,3))
def test_all_order_block_lemma_calibrations(n):
    data = rank_reduction(n)
    y = data['variables']
    assert s.cancel(data['F']-data['reduced']-(y[-1]-data['eta'])**2/data['sigma']) == 0
    denominator = s.fraction(data['reduced'])[1]
    previous = data['determinants'][-2] if n > 1 else 1
    assert not s.cancel(denominator/previous).free_symbols
    moved = tuple(sum(s.binomial(j,k)*(-c)**(j-k)*y[k] for k in range(j+1))
                  for j in range(2*n))
    at = dict(zip(y,moved))
    for det in data['determinants']:
        assert s.expand(det.xreplace(at)-det) == 0
    b = exponential_coefficients(y)
    for j, bj in enumerate(exponential_coefficients(moved)):
        assert s.expand(bj-sum(s.binomial(j,k)*(-c)**(j-k)*b[k] for k in range(j+1))) == 0


def test_reduction_agrees_with_previous_explicit_base():
    data = rank_reduction(2)
    assert s.cancel(data['reduced'].xreplace(dict(zip(data['variables'],X)))-h0()) == 0


def density_moments(a,b,count=5):
    return tuple(s.Rational(a)/(j+1)+s.Rational(b)/(j+2) for j in range(count))


@pytest.mark.parametrize('fraction', (s.Rational(1,7),s.Rational(1,2),s.Rational(6,7)))
def test_rational_moment_lift_on_exact_two_point_functionals(fraction):
    points = (density_moments(s.Rational(1,3),s.Rational(1,5))[:4],
              density_moments(s.Rational(2,3),s.Rational(-1,7))[:4])
    weights = (fraction,1-fraction)
    def functional(poly):
        return sum(w*s.cancel(poly/P).subs(dict(zip(X,point)))
                   for point,w in zip(points,weights))
    assert functional(P) == 1
    output = tuple(functional(P*t) for t in X)
    mo,zo,qo,ro = output
    assert mo > 0
    # The 1-block on (m,Delta) proves Z(Delta^2)>0 since Z(m Delta)=1.
    mass_gap = functional(DELTA**2)
    assert mass_gap > 0
    perspective_gap = functional(DELTA*z*z)-zo*zo/mo
    cubic_gap = functional(P*m**3)-mo**3
    assert perspective_gap >= 0 and cubic_gap >= 0
    assert qo-zo*zo/mo-mo**3/12 == mass_gap+perspective_gap+cubic_gap/12
    assert DELTA.subs(dict(zip(X,output))) > 0
    L = fifth_schur()[2]
    value = functional(s.cancel(P*L))
    assert value >= L.subs(dict(zip(X,output)))
    for point in points:
        at = dict(zip(X,point))
        assert all(g.subs(at)>0 for g in (m,DELTA,m*DELTA))


@pytest.mark.parametrize('epsilon', (s.Rational(1,2),s.Rational(1,10),s.Rational(1,100)))
@pytest.mark.parametrize('bands', ((), ((0,1),), ((s.Rational(1,4),s.Rational(3,4)),),
                                  ((0,s.Rational(1,4)),(s.Rational(1,2),s.Rational(3,4)))))
def test_singular_physical_boundary_approached_not_substituted(epsilon,bands):
    original = tuple(sum((b**(j+1)-a**(j+1))/s.Integer(j+1) for a,b in bands)
                     for j in range(5))
    center = tuple(s.Rational(1,2*(j+1)) for j in range(5))
    mixed = tuple((1-epsilon)*a+epsilon*b for a,b in zip(original,center))
    comp = tuple(s.Rational(1,j+1)-mixed[j] for j in range(5))
    L = fifth_schur()[2]
    for values in (mixed,comp):
        at = dict(zip(X,values[:4]))
        assert P.subs(at)>0
        assert L.subs(at)<=values[4]
    assert all(a-b == epsilon*(c0-b) for a,b,c0 in zip(mixed,original,center))


def test_no_finite_point_evaluation_at_zero_mass():
    assert P.subs(dict(zip(X,(0,0,0,0)))) == 0
    # Normalize Z(P)=1 only at positive denominators, then take output closure.
    epsilon = s.Symbol('epsilon', positive=True)
    point = tuple(epsilon/(2*(j+1)) for j in range(4))
    reciprocal = 1/P.subs(dict(zip(X,point)))
    assert s.limit(reciprocal,epsilon,0,dir='+') == s.oo
