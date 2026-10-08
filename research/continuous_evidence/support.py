"""Certified continuous-time support, using inner LPs only as proposals.

Unknown jerk is arbitrary measurable, not restricted to the LP mesh. Every
accepted bound has an exact rational trajectory/dual-Bernstein certificate.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial
import warnings
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import qr
import sympy as sp


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


@dataclass(frozen=True)
class Observation:
    time: F
    order: int  # 0 position, 1 velocity, 2 actual acceleration, NOT a nominal input
    value: F
    error: F
    received: F

    def record(self):
        return dict(time=str(self.time), order=self.order, value=str(self.value),
                    error=str(self.error), received=str(self.received))


def initial_row(t, weights):
    return [sum((weights[r]*t**(k-r)/factorial(k-r) for r in range(k+1)), F(0))
            for k in range(3)]


def kernel(t, weights, s):
    if s > t:
        return F(0)
    return sum((weights[r]*(t-s)**(2-r)/factorial(2-r) for r in range(3)), F(0))


def integrated_kernel(t, weights, left, right):
    right = min(right, t)
    if left >= right:
        return F(0)
    return sum((weights[r]*((t-left)**(3-r)-(t-right)**(3-r))/factorial(3-r)
                for r in range(3)), F(0))


def row(t, weights, mesh):
    return initial_row(t, weights)+[integrated_kernel(t, weights, l, r)
                                    for l, r in zip(mesh, mesh[1:])]


def basis(order):
    return tuple(F(k == order) for k in range(3))


def bernstein(coeff, left, right):
    """Power polynomial c0+c1*s+c2*s^2, exact degree-2 Bernstein coefficients."""
    a, b, c = coeff
    h = right-left
    p = a+b*left+c*left*left
    return p, p+h*(b+2*c*left)/2, a+b*right+c*right*right


def polynomial(t, weights):
    p, v, a = weights
    return [p*t*t/2+v*t+a, -p*t-v, p/2]


def dual_certificate(center, radii, jerk, observations, target, weights, lam,
                     integral_tolerance=F(1, 10**10)):
    """Any rational multipliers give a rigorous upper support bound."""
    c = initial_row(target, weights)
    value = F(0)
    knots = sorted({F(0), target} | {o.time for o in observations})
    for o, l in zip(observations, lam):
        c = [x-l*y for x, y in zip(c, initial_row(o.time, basis(o.order)))]
        value += l*o.value+abs(l)*o.error
    value += dot(c, center)+dot(list(map(abs, c)), radii)
    leaves = []
    refine_scores = []
    def enclose(poly, l, r, depth=0):
        bs = bernstein(poly, l, r)
        upper = (r-l)*sum(map(abs, bs), F(0))/3
        lower = abs((r-l)*sum(bs, F(0))/3)
        if upper-lower <= integral_tolerance*(r-l)/target or depth >= 28:
            leaves.append(dict(left=str(l), right=str(r), polynomial=list(map(str, poly)),
                               bernstein=list(map(str, bs)), upper=str(upper)))
            return upper
        mid = (l+r)/2
        return enclose(poly, l, mid, depth+1)+enclose(poly, mid, r, depth+1)
    integral = F(0)
    for l, r in zip(knots, knots[1:]):
        mid = (l+r)/2
        poly = polynomial(target, weights)
        for o, weight in zip(observations, lam):
            if mid < o.time:
                poly = [x-weight*y for x, y in zip(poly, polynomial(o.time, basis(o.order)))]
        integral += enclose(poly, l, r)
    value += jerk*integral
    return dict(multipliers=list(map(str, lam)), residual_initial=list(map(str, c)),
                leaves=leaves, integral_upper=str(integral), upper=str(value))


def recover_vertex(result, A, b, bounds):
    """Round bound variables; solve a small independent active measurement system."""
    n = len(bounds)
    x = [None]*n
    free = []
    for i, ((lo, hi), numeric) in enumerate(zip(bounds, result.x)):
        if abs(numeric-float(lo)) < 1e-7:
            x[i] = lo
        elif abs(numeric-float(hi)) < 1e-7:
            x[i] = hi
        else:
            free.append(i)
    if free:
        active = [i for i, slack in enumerate(result.ineqlin.residual) if abs(slack) < 1e-7]
        M = np.array([[float(A[i][j]) for j in free] for i in active])
        assert len(active) >= len(free)
        _, R, pivot = qr(M.T, pivoting=True, mode='economic')
        selected = [active[i] for i in pivot[:len(free)]]
        exact_A = sp.Matrix([[sp.Rational(A[i][j]) for j in free] for i in selected])
        exact_b = sp.Matrix([sp.Rational(b[i]-sum((A[i][j]*v for j, v in enumerate(x)
                                                if v is not None), F(0))) for i in selected])
        solved = exact_A.inv()*exact_b
        for j, value in zip(free, solved):
            x[j] = F(int(value.p), int(value.q))
    assert all(lo <= v <= hi for v, (lo, hi) in zip(x, bounds))
    assert all(dot(a, x) <= rhs for a, rhs in zip(A, b))
    return x


def solve(center, radii, jerk, observations, now, target, weights=(F(1), F(0), F(0)),
          tolerance=F(1, 10**6), max_iterations=14):
    assert 0 < target and now <= target and jerk > 0
    assert len(center) == len(radii) == len(weights) == 3 and min(radii) >= 0
    assert all(0 <= o.time <= o.received <= now and o.error >= 0 and o.order in (0, 1, 2)
               for o in observations), 'Undelivered/future evidence must not enter a causal bound'
    mesh = sorted({F(0), target} | {o.time for o in observations})
    logs = []
    for iteration in range(max_iterations):
        objective = row(target, weights, mesh)
        A, b = [], []
        for o in observations:
            a = row(o.time, basis(o.order), mesh)
            A += [a, [-v for v in a]]
            b += [o.value+o.error, -o.value+o.error]
        bounds = [(c-r, c+r) for c, r in zip(center, radii)]+[(-jerk, jerk)]*(len(mesh)-1)
        with warnings.catch_warnings(record=True) as notices:
            warnings.simplefilter('always')
            result = linprog(-np.array(list(map(float, objective))),
                             A_ub=np.array([[float(v) for v in a] for a in A]) if A else None,
                             b_ub=np.array(list(map(float, b))) if b else None,
                             bounds=[tuple(map(float, bound)) for bound in bounds],
                             method='highs', options=dict(threads=1))
        if result.status == 2:
            # An infeasible inner mesh DOES NOT prove the continuous history inconsistent.
            logs.append(dict(status='inner_mesh_infeasible', cells=len(mesh)-1))
            mesh = sorted(set(mesh) | {(l+r)/2 for l, r in zip(mesh, mesh[1:])})
            continue
        assert result.status == 0, result.message
        values = recover_vertex(result, A, b, bounds)
        lam = [F(float(-result.ineqlin.marginals[2*i]+result.ineqlin.marginals[2*i+1])).limit_denominator(10**9)
               for i in range(len(observations))]
        dual = dual_certificate(center, radii, jerk, observations, target, weights, lam)
        lower = dot(objective, values)
        gap = F(dual['upper'])-lower
        assert gap >= 0
        logs.append(dict(status='certified', cells=len(mesh)-1, gap=str(gap),
                         solver_notices=[str(n.message) for n in notices]))
        record = dict(center=list(map(str, center)), radii=list(map(str, radii)), jerk=str(jerk),
                      observations=[o.record() for o in observations], now=str(now), target=str(target),
                      weights=list(map(str, weights)), mesh=list(map(str, mesh)),
                      trajectory=list(map(str, values)), lower=str(lower), dual=dual,
                      gap=str(gap), tolerance=str(tolerance), converged=gap <= tolerance, iterations=logs)
        if gap <= tolerance:
            return record
        # Refine only cells where the continuous residual changes sign. No jerk
        # discretization is used in the dual certificate, so the bound remains valid.
        scores = []
        for l, r in zip(mesh, mesh[1:]):
            poly = polynomial(target, weights)
            for o, weight in zip(observations, lam):
                if (l+r)/2 < o.time:
                    poly = [x-weight*y for x, y in zip(poly, polynomial(o.time, basis(o.order)))]
            bs = bernstein(poly, l, r)
            score = (r-l)*(sum(map(abs, bs), F(0))-abs(sum(bs, F(0))))/3
            if score > 0:
                scores.append((score, (l+r)/2))
        if not scores:
            return record  # Preserve a certified but nonconverged result, never relabel it.
        mesh = sorted(set(mesh) | {mid for _, mid in sorted(scores, reverse=True)[:6]})
    if logs and logs[-1]['status'] == 'certified':
        return record
    raise RuntimeError('No feasible continuous witness recovered; inner infeasibility is inconclusive')
