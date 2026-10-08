"""Exact polynomial-space reduction and finite singular-boundary witnesses."""
from functools import lru_cache
import sympy as s
from research.schur_sparse_lift.algebra import blocks as original_blocks
from research.schur_remainder_recursion.algebra import X, m, z, q, r, P, DELTA, N, fifth_schur

CENTER = tuple(s.Rational(1, 2*(j+1)) for j in range(5))


def blocks():
    return tuple(item for item in original_blocks()
                 if item[0] not in ('unit_domain', 'delta_domain'))


@lru_cache(maxsize=1)
def model():
    polys = [s.expand(P), *(s.expand(P*x) for x in X), s.expand(s.cancel(P*fifth_schur()[2]))]
    entries = []
    for name, weight, basis in blocks():
        locations = []
        for i in range(len(basis)):
            for j in range(i, len(basis)):
                locations.append((i,j,len(polys)))
                polys.append(s.expand(weight*basis[i]*basis[j]))
        entries.append((name, len(basis), locations))
    support = tuple(sorted({a for p in polys for a in s.Poly(p, X).monoms()}))
    W = s.Matrix([[s.Poly(p,X).coeff_monomial(a) for a in support] for p in polys])
    _, row_indices = W.T.rref()
    B = W[list(row_indices), :]
    _, columns = B.rref()
    coordinates = W[:,list(columns)]*B[:,list(columns)].inv()
    assert coordinates*B == W
    theta = s.symbols('theta0:'+str(B.rows))
    expressions = coordinates*s.Matrix(theta)
    pencils = {}
    for name, size, locations in entries:
        mat = s.zeros(size)
        for i,j,k in locations:
            mat[i,j] = mat[j,i] = expressions[k]
        pencils[name] = mat
    return dict(basis=tuple(polys[k] for k in row_indices), theta=theta,
                matrices=pencils, outputs=tuple(expressions[:6]),
                coordinates=coordinates, source_polynomials=tuple(polys),
                row_indices=row_indices, monomials=support)


def evaluate(values):
    data = model()
    at = dict(zip(data['theta'],values))
    return ({name:mat.xreplace(at) for name,mat in data['matrices'].items()},
            tuple(p.xreplace(at) for p in data['outputs']))


@lru_cache(maxsize=128)
def witness(prefix):
    """Rational finite functional, including zero denominator via a regular path.

    Caller supplies a physical prefix. This function is not a feasibility oracle.
    Every chosen basis polynomial is an actual retained entry or output, so its
    evaluation quotient has a removable singularity on the proved path.
    """
    prefix = tuple(map(s.Rational,prefix))
    at = dict(zip(X,prefix))
    denominator = P.subs(at)
    if denominator > 0:
        values = tuple(s.cancel(f/P).subs(at) for f in model()['basis'])
        return values, 0
    if denominator != 0:
        raise ValueError('negative denominator is not a physical prefix')
    e = s.Symbol('eps', positive=True)
    path = dict(zip(X,((1-e)*a+e*b for a,b in zip(prefix,CENTER))))
    denominator_path = s.Poly(s.expand(P.subs(path)),e)
    valuation = min(power[0] for power,coef in denominator_path.terms() if coef != 0)
    values = []
    for f in model()['basis']:
        quotient = s.cancel(f.subs(path)/denominator_path.as_expr())
        value = quotient.subs(e,0)
        if not value.is_Rational:
            raise AssertionError(('nonremovable coordinate', f, quotient, value))
        values.append(value)
    return tuple(values), valuation


def interval_point(intervals):
    return tuple(sum((s.Rational(b)**(j+1)-s.Rational(a)**(j+1))/s.Integer(j+1)
                     for a,b in intervals) for j in range(5))


def complement(point):
    return tuple(s.Rational(1,j+1)-v for j,v in enumerate(point))


def boundary_cases():
    return (
        ('zero', interval_point(())),
        ('full', interval_point(((0,1),))),
        ('left_interval', interval_point(((0,s.Rational(2,5)),))),
        ('right_interval', interval_point(((s.Rational(3,5),1),))),
        ('interior_interval', interval_point(((s.Rational(1,4),s.Rational(3,4)),))),
        ('two_intervals', interval_point(((0,s.Rational(1,4)),(s.Rational(1,2),s.Rational(3,4))))),
        ('four_switches', interval_point(((0,s.Rational(1,5)),(s.Rational(2,5),s.Rational(3,5)),(s.Rational(4,5),1)))),
        ('constant_half', CENTER),
        ('near_zero', tuple(v/100 for v in CENTER)),
        ('near_full', complement(tuple(v/100 for v in CENTER))),
    )
