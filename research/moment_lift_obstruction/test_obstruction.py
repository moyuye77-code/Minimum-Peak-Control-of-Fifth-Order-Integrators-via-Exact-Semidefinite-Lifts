from copy import deepcopy
from fractions import Fraction as F
from math import comb, factorial
import pytest
import sympy as sp
from .exact import integral, value, moments, squared_root_polynomial, input_state, witness_record
from .verify import build, audit


@pytest.mark.parametrize('d', range(3, 13))
def test_arbitrary_order_witness_examples(d):
    record = witness_record(d)
    k = (d-1)//2
    assert len(record['witnesses']) == comb(k+3, k)
    assert record['lower_bound'] == (d+1)//2
    t = sp.Symbol('t')
    for row in record['witnesses']:
        roots = list(map(sp.Rational, row['roots']))
        symbolic = sp.prod((t-a)**2 for a in roots).expand()
        assert [sp.Poly(symbolic, t).nth(j) for j in range(2*k+1)] == list(map(sp.Rational, row['coefficients']))
        assert sp.integrate(symbolic, (t, 0, 1)) == sp.Rational(row['integral'])


@pytest.mark.parametrize('d', [1, 2, 3, 4, 5, 8, 12])
def test_moment_to_integrator_map_against_symbolic_integration(d):
    pieces = ((F(0), F(1, 5), F(1)), (F(1, 3), F(2, 3), F(2, 7)),
              (F(3, 4), F(1), F(1)))
    horizon, bound = F(7, 3), F(5, 2)
    y = moments(pieces, d)
    state = input_state(y, horizon, bound)
    t = sp.Symbol('t')
    for k in range(d):
        kernel = (1-t)**k
        density_part = sum(sp.Rational(level.numerator, level.denominator)*
                           sp.integrate(kernel, (t, sp.Rational(a.numerator, a.denominator),
                                                     sp.Rational(b.numerator, b.denominator)))
                           for a, b, level in pieces)
        exact = sp.Rational(bound.numerator, bound.denominator)*sp.Rational(horizon.numerator, horizon.denominator)**(k+1)/factorial(k)*(
            2*density_part-sp.integrate(kernel, (t, 0, 1)))
        assert exact == sp.Rational(state[k].numerator, state[k].denominator)


def test_negative_polynomial_region_invalidates_saturated_normal():
    # p(t)=t-1/2. Removing the negative half strictly improves the functional.
    coefficients = (F(-1, 2), F(1))
    y = moments(((F(1, 2), F(1), F(1)),), 2)
    full = (F(1), F(1, 2))
    gain = sum(c*(a-b) for c, a, b in zip(coefficients, y, full))
    assert gain == F(1, 8) == -integral(coefficients, F(0), F(1, 2))
    centered_support = (-integral(coefficients, F(0), F(1, 2))+
                        integral(coefficients, F(1, 2), F(1)))/2
    assert centered_support == F(1, 8)


def test_endpoint_negativity_is_not_lost():
    p = (F(-1, 4), F(1))
    assert value(p, F(0)) < 0
    assert -integral(p, F(0), F(1, 4)) == F(1, 32)


def test_polynomial_gram_identity_for_two_roots():
    t = sp.Symbol('t')
    a, b = sp.Rational(1, 3), sp.Rational(2, 3)
    z = sp.Matrix([1, t, t*t])
    v = sp.Matrix([a*b, -a-b, 1])
    assert sp.expand((z.T*(v*v.T)*z)[0] - (t-a)**2*(t-b)**2) == 0
    assert list(map(sp.Rational, squared_root_polynomial((F(1, 3), F(2, 3))))) == [sp.Poly((t-a)**2*(t-b)**2, t).nth(j) for j in range(5)]


def test_bad_archive_and_false_experiment_count_rejected():
    data = build()
    assert audit(data)['orders'] == 10
    bad = deepcopy(data)
    bad['orders'][0]['witnesses'][0]['evaluations'][0] = '999'
    with pytest.raises(AssertionError):
        audit(bad)
    bad = deepcopy(data)
    bad['new_closed_loop_runs'] = 1
    with pytest.raises(AssertionError):
        audit(bad)


def test_degenerate_horizon_and_invalid_density_rejected():
    with pytest.raises(ValueError):
        input_state((F(1),), horizon=F(0))
    with pytest.raises(ValueError):
        moments(((F(0), F(1), F(2)),), 3)
    with pytest.raises(ValueError):
        witness_record(2)
