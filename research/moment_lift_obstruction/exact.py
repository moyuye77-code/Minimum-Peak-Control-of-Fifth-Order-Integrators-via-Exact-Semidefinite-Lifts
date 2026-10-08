"""Small rational polynomial primitives, not a nonrepresentability oracle."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial


def multiply(a, b):
    result = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return tuple(result)


def value(coefficients, t):
    result = F(0)
    for coefficient in reversed(coefficients):
        result = result*t + coefficient
    return result


def integral(coefficients, left=F(0), right=F(1)):
    return sum((c*(right**(j+1)-left**(j+1))/F(j+1)
                for j, c in enumerate(coefficients)), F(0))


def squared_root_polynomial(roots):
    coefficients = (F(1),)
    for root in roots:
        coefficients = multiply(coefficients, (root*root, -2*root, F(1)))
    return coefficients


def witness_record(d):
    if d < 3:
        raise ValueError('The nontrivial polynomial witness starts at order 3')
    k = (d-1)//2
    grid = tuple(F(i, k+4) for i in range(1, k+4))
    witnesses = []
    for roots in combinations(grid, k):
        coefficients = squared_root_polynomial(roots)
        area = integral(coefficients)
        assert len(coefficients)-1 == 2*k <= d-1 and area > 0
        evaluations = [value(coefficients, t) for t in grid]
        assert all((v == 0 if t in roots else v > 0)
                   for t, v in zip(grid, evaluations))
        normalized = tuple(2*c/area for c in coefficients)
        assert integral(normalized) == 2
        # Dot product with the centered saturated point equals one.
        assert sum((c/F(2*(j+1)) for j, c in enumerate(normalized)), F(0)) == 1
        witnesses.append(dict(roots=list(map(str, roots)),
                              coefficients=list(map(str, coefficients)),
                              integral=str(area), evaluations=list(map(str, evaluations))))
    return dict(order=d, k=k, lower_bound=k+1, grid=list(map(str, grid)),
                witnesses=witnesses)


def moments(pieces, d):
    """Integrate disjoint density pieces (left,right,level) in [0,1]."""
    if d < 1:
        raise ValueError('Positive moment dimension required')
    previous = F(0)
    result = [F(0)]*d
    for left, right, level in pieces:
        if not previous <= left < right <= 1 or not 0 <= level <= 1:
            raise ValueError('Ordered disjoint bounded-density pieces required')
        for j in range(d):
            result[j] += level*(right**(j+1)-left**(j+1))/F(j+1)
        previous = right
    return tuple(result)


def input_state(y, horizon=F(1), bound=F(1)):
    """Return input contributions in reverse derivative order, initial flow excluded."""
    if horizon <= 0 or bound <= 0:
        raise ValueError('Nondegenerate horizon and input bound required')
    return tuple(bound*horizon**(k+1)/F(factorial(k)) *
                 (2*sum(((-1)**j*comb(k, j)*y[j] for j in range(k+1)), F(0))
                  - F(1, k+1)) for k in range(len(y)))
