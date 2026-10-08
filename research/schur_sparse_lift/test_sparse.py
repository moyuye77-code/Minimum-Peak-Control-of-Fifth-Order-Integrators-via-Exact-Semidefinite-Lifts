import pytest
import sympy as s
from .algebra import (
    blocks, square_coordinates, certificate_squares, affine_model, records,
    X, U, m, z, q, r, M, DELTA, P,
)


@pytest.mark.parametrize('index', range(8))
def test_all_joint_squares_have_fixed_sparse_spaces(index):
    _, weight, _, root = certificate_squares()[index]
    name, coordinates = square_coordinates()[index]
    actual_weight, basis = {a:(b,c) for a,b,c in blocks()}[name]
    assert s.expand(weight-actual_weight) == 0
    assert s.expand(root-s.Add(*(a*b for a,b in zip(coordinates,basis)))) == 0
    assert all(not (set(X)&c.free_symbols) for c in coordinates)


def test_affine_coefficient_model_round_trip():
    model = affine_model()
    monomials = tuple(s.prod(x**power for x,power in zip(X,exponent))
                      for exponent in model['support'])
    at = dict(zip(model['coefficients'],monomials))
    for name,matrix in model['matrices'].items():
        recovered = matrix.xreplace(at)
        assert recovered.applyfunc(s.expand) == model['polynomial_matrices'][name]
        assert matrix == matrix.T
        assert all(s.Poly(t,model['coefficients']).total_degree()<=1 for t in matrix)
    for expression,poly in zip(model['outputs'],model['output_polynomials']):
        assert s.expand(expression.xreplace(at)-poly) == 0
    assert model['outputs'][0] in model['coefficients'] or s.Poly(
        model['outputs'][0],model['coefficients']).total_degree() == 1


def test_strict_domain_proof_has_all_required_subspaces():
    named = {name:(weight,basis) for name,weight,basis in blocks()}
    assert named['unit_domain'] == (1,(m,DELTA))
    assert named['delta_domain'] == (DELTA,(1,m))
    # (m,z) handles the perspective; (m^2,m) handles the cubic remainder.
    assert m in named['delta_B'][1] and z in named['delta_B'][1]
    assert m*m in named['delta_h'][1] and m in named['delta_h'][1]
    assert m in named['p_h'][1] and 1 in named['p_h'][1]
    assert s.expand(P*(q-z*z/m-m**3/12)-DELTA**2) == 0
    # After applying ell, use ell(P)=1 and ell(P*m)=M.
    rhs = DELTA*(m*(m-M))**2+2*M*P*(m-M)**2
    assert s.expand(rhs-(P*m**3-3*M*M*P*m+2*M**3*P)) == 0


@pytest.mark.parametrize('rho', (s.Rational(1,4),s.Rational(1,2),s.Rational(3,4)))
def test_point_evaluation_witness_all_sparse_blocks(rho):
    point = tuple(rho/s.Integer(i+1) for i in range(4))
    at = dict(zip(X,point))
    denominator = P.subs(at)
    assert denominator > 0
    model = affine_model()
    evaluation = {var:s.prod(a**power for a,power in zip(point,exponent))/denominator
                  for var,exponent in zip(model['coefficients'],model['support'])}
    assert model['outputs'][0].subs(evaluation) == 1
    assert tuple(expr.subs(evaluation) for expr in model['outputs'][1:5]) == point
    for name,weight,basis in blocks():
        vector = s.Matrix(basis).subs(at)
        factor = weight.subs(at)/denominator
        assert factor > 0
        assert model['matrices'][name].subs(evaluation) == factor*vector*vector.T


def test_exact_record_and_block_ceiling():
    data = records()
    assert data['square_residuals'] == ['0']*8
    assert len(data['blocks']) == 9
    assert data['maximum_block_order'] == 4
    assert data['optimization_calls'] == 0
    assert data['closed_model_coefficients_constructed'] is False
