"""Exact affine coefficient model; no optimization or numerical SOS fitting."""
from functools import lru_cache
import sympy as s
from research.schur_remainder_recursion.algebra import (
    X, U, m, z, q, r, M, Z, Q, R, DELTA, DELTA0, N, N0, P, P0,
    certificate_squares, fifth_schur,
)


def blocks():
    # name, weight, polynomial feature vector. All share ONE linear functional.
    return (
        ('unit_A', s.S.One, (N, DELTA)),
        ('unit_h', s.S.One, (P, DELTA)),
        ('unit_domain', s.S.One, (m, DELTA)),
        ('delta_B', DELTA, (q, m, m*m, z)),
        ('delta_D', DELTA, (m*z, z, m, m*m)),
        ('delta_h', DELTA, (m**3, m*m, m)),
        ('p_C', P, (z, m, s.S.One)),
        ('p_h', P, (m*m, m, s.S.One)),
        ('delta_domain', DELTA, (s.S.One, m)),
    )


def square_coordinates():
    return (
        ('unit_A', (M*DELTA0, -M*N0)),
        ('delta_B', (M*DELTA0, -Q*DELTA0-M**3*DELTA0/12+N0*Z,
                     M*M*DELTA0/12, -M*N0)),
        ('p_C', (M*DELTA0, Z*DELTA0-M*N0,
                 -2*M*Z*DELTA0+M*M*N0)),
        ('delta_D', (2*M*DELTA0, -M*M*DELTA0,
                     -M*Z*DELTA0+M*M*N0, -M*N0)),
        ('unit_h', (M*DELTA0, -M*M*DELTA0)),
        ('p_h', (s.S.Zero, DELTA0, -M*DELTA0)),
        ('p_h', (M*DELTA0, -2*M*M*DELTA0, M**3*DELTA0)),
        ('delta_h', (M*DELTA0, -2*M*M*DELTA0, M**3*DELTA0)),
    )


def polynomial(poly):
    return s.Poly(s.cancel(poly), X)


@lru_cache(maxsize=1)
def affine_model():
    matrix_polys = {
        name: s.Matrix(len(basis), len(basis),
                       lambda i, j: s.expand(weight*basis[i]*basis[j]))
        for name, weight, basis in blocks()
    }
    output_polys = (P,) + tuple(P*t for t in X) + (s.cancel(P*fifth_schur()[2]),)
    polys = [p for mat in matrix_polys.values() for p in mat] + list(output_polys)
    support = tuple(sorted({exponent for p in polys for exponent in polynomial(p).monoms()}))
    coefficients = s.symbols('ell0:'+str(len(support)))
    lookup = dict(zip(support, coefficients))

    def linear(poly):
        return s.Add(*(coef*lookup[exponent] for exponent, coef in polynomial(poly).terms()))

    matrices = {name: mat.applyfunc(linear) for name, mat in matrix_polys.items()}
    outputs = tuple(linear(p) for p in output_polys)
    return dict(support=support, coefficients=coefficients, polynomial_matrices=matrix_polys,
                matrices=matrices, output_polynomials=output_polys, outputs=outputs)


def records():
    data = affine_model()
    named = {name: (weight, basis) for name, weight, basis in blocks()}
    residuals = []
    for term, (name, coords) in zip(certificate_squares(), square_coordinates()):
        coefficient, gx, gu, root = term
        weight, basis = named[name]
        assert s.expand(gx-weight) == 0
        assert coefficient > 0
        assert all(not (set(X) & value.free_symbols) for value in coords)
        residuals.append(s.expand(root-s.Add(*(a*b for a,b in zip(coords,basis)))))
    assert residuals == [0]*8
    for matrix in data['matrices'].values():
        assert matrix == matrix.T
        assert all(s.Poly(entry, data['coefficients']).total_degree() <= 1 for entry in matrix)
    assert all(s.Poly(entry, data['coefficients']).total_degree() <= 1 for entry in data['outputs'])
    return dict(schema=1, square_residuals=list(map(str,residuals)),
                blocks=[dict(name=name, weight=str(weight), basis=list(map(str,basis)),
                             order=len(basis)) for name,weight,basis in blocks()],
                scalar_functional_coordinates=len(data['support']),
                maximum_block_order=max(len(basis) for _,_,basis in blocks()),
                monomial_support=[list(item) for item in data['support']],
                full_K5_extension_degree_lower=3, full_K5_extension_degree_upper=4,
                closed_model_coefficients_constructed=False,
                priority_cleared=False, optimization_calls=0, new_closed_loops=0)
