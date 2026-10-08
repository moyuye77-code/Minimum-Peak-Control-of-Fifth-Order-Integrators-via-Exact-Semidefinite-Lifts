"""Rational matrix witnesses and symbolic identities for the joint lift.

The six auxiliary entries describe a positive linear functional on selected
polynomials in the COORDINATES (m,z), not higher time moments of a density.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
import sympy as sp

from research.moment_upper_gate.exact import psd_pivots, fourth_bounds


def atomic_aux(m, z):
    m, z = F(m), F(z)
    return (m*m, m**3, m**4, m*z, m*m*z, z*z)


def matrices(point, aux, compressed=True):
    """Only affine entries; works for Fraction, sympy, or CVXPY scalars."""
    m, z, q, r = point
    a, b, c, d, e, f = aux
    if compressed:
        result = [[[1, m, a], [m, a, b], [a, b, c]]]
        for sign in (1, -1):
            result.append([[1, m, z+sign*a/6],
                           [m, a, d+sign*b/6],
                           [z+sign*a/6, d+sign*b/6, f+sign*e/3+c/36]])
    else:
        result = [[[1, m, z, a], [m, a, d, b],
                   [z, d, f, e], [a, b, e, c]]]
    result.extend([
        [[z-a/2, d-b/2], [d-b/2, e-c/2]],
        [[m-z-a/2, a-d-b/2], [a-d-b/2, b-e-c/2]],
        [[z-a/2, q-d/2-b/12], [q-d/2-b/12, r-f/2-e/4]],
        [[m-z-a/2, -a/4+b/12-d/2+m/2-q+z/2],
         [-a/4+b/12-d/2+m/2-q+z/2,
          -a/8-d/2+e/4-f/2+m/4-r+3*z/4]],
    ])
    return result


def interval_point(intervals):
    return tuple(sum((F(hi)**(j+1)-F(lo)**(j+1))/F(j+1)
                     for lo, hi in intervals) for j in range(4))


def polynomial_direction(roots, sign=1):
    coeff = [F(sign)]
    for root in roots:
        out = [F(0)]*(len(coeff)+1)
        for j, value in enumerate(coeff):
            out[j] -= F(root)*value
            out[j+1] += value
        coeff = out
    return tuple(coeff+[F(0)]*(4-len(coeff)))


def direction_record(roots, sign=1):
    coeff = polynomial_direction(roots, sign)
    knots = sorted({F(0), F(1), *(F(x) for x in roots if 0 < x < 1)})
    intervals = []
    for lo, hi in zip(knots, knots[1:]):
        mid = (lo+hi)/2
        if sum(c*mid**j for j, c in enumerate(coeff)) > 0:
            intervals.append((lo, hi))
    point = interval_point(intervals)
    true_support = sum(c*y for c, y in zip(coeff, point))
    return dict(roots=list(map(str, roots)), sign=sign,
                direction=list(map(str, coeff)),
                intervals=[[str(lo), str(hi)] for lo, hi in intervals],
                point=list(map(str, point)), support=str(true_support))


def directions():
    grid = (F(1,4), F(1,2), F(3,4))
    roots = [()]
    for degree in (1, 2, 3):
        roots.extend(combinations_with_replacement(grid, degree))
    roots.extend([(F(0), F(1)), (F(-1), F(1,3), F(2)),
                  (F(1,100), F(1,2), F(99,100))])
    return [direction_record(rr, sign) for rr in roots for sign in (1, -1)]


def point_record(point, label):
    point = tuple(map(F, point))
    lo, hi = fourth_bounds(*point[:3])
    assert lo <= point[3] <= hi
    pivots = [[str(x) for x in psd_pivots(matrix)]
              for matrix in matrices(point, atomic_aux(*point[:2]))]
    return dict(label=label, point=list(map(str, point)),
                auxiliary=list(map(str, atomic_aux(*point[:2]))), pivots=pivots)


def witness_records():
    records = [point_record(row['point'], f'direction_{i}')
               for i, row in enumerate(directions())]
    # Both singular faces, interior densities, and lower/upper fourth bounds.
    for m in (F(0), F(1,100), F(1,4), F(1,2), F(3,4), F(99,100), F(1)):
        records.append(point_record(interval_point([(0,m)]), f'prefix_{m}'))
        records.append(point_record(interval_point([(1-m,1)]), f'suffix_{m}'))
        records.append(point_record(tuple(m/F(j+1) for j in range(4)), f'uniform_{m}'))
    for m,z,q in ((F(1,2),F(1,4),F(1,6)),
                  (F(1,4),F(1,8),F(1,12)),
                  (F(3,4),F(3,8),F(1,4))):
        lo, hi = fourth_bounds(m,z,q)
        for name,r in (('lo',lo),('mid',(lo+hi)/2),('hi',hi)):
            records.append(point_record((m,z,q,r), f'fiber_{m}_{name}'))
    return records


def identity_residuals():
    m,z,q,M,Z,g = sp.symbols('m z q M Z g')
    D = z-m*m/2
    D0 = Z-M*M/2
    Fg = z*z/2+m*m*z/4-g*(m**3/12+m*z/2)+g*g*(m*m/8-z/4)
    ref = {m:M,z:Z}
    remainder = Fg-Fg.subs(ref)-sum(sp.diff(Fg,v).subs(ref)*(v-ref[v]) for v in (m,z))
    u,v = m-M,z-Z
    sos = (v+(M-g)*u/2+u*u/6)**2/2 + D*u*u/12+D0*u*u/6+u**4/36
    L = z*z/2+m*m*z/4+(q-m*z/2-m**3/12)**2/D
    old = m*q+z*z/2-m*m*z/2+m**4/24+(q-m*z+m**3/6)**2/D
    E = m-z-m*m/2
    lo = z*z/m+m**3/12
    hi = sp.Rational(1,3)-(sp.Rational(1,2)-z)**2/(1-m)-(1-m)**3/12
    Lbar = L.subs({m:1-m,z:sp.Rational(1,2)-z,q:sp.Rational(1,3)-q}, simultaneous=True)
    return {name:str(sp.factor(value)) for name,value in dict(
        bregman=remainder-sos,
        boundary=L-old,
        prefix_gap=hi-lo-D*E/(m*(1-m)),
        fourth_gap=L+Lbar-sp.Rational(1,4)-m*(1-m)*(q-lo)*(q-hi)/(D*E),
        conjugate=(g*q+Fg).subs(g,2*(q-m*z/2-m**3/12)/D)-L,
    ).items()}
