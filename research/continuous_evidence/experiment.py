"""Continuous measurement histories and one-hold planar separation tests.

This is NOT a recurrent distributed request controller or a full flocking trial.
Payload alternatives are fixed in advance; no reply is used before its arrival.
"""
from fractions import Fraction as F
from itertools import combinations
import warnings
import numpy as np
from scipy.optimize import linprog
from .support import Observation, solve, dot


def truth(seed, t):
    p, v, a = F(0), F(0), F(seed-1, 20)
    stamps = [F(0), F(3, 10), F(11, 20), F(4, 5), F(7, 5), F(2)]
    jerks = [F(1, 4), F(-1, 3), F(1, 5), F(-1, 4), F(1, 3)]
    if seed % 2:
        jerks = [-x for x in jerks]
    for l, r, j in zip(stamps, stamps[1:], jerks):
        h = max(F(0), min(t, r)-l)
        p, v, a = p+h*v+h*h*a/2+h**3*j/6, v+h*a+h*h*j/2, a+h*j
    return p, v, a


def observations(seed, payload=None):
    data = []
    for k in range(6):
        t = F(k, 5)
        for order, error in [(0, F(1, 25)), (1, F(1, 25))]:
            noise = error*F(((k+seed+order) % 5)-2, 3)
            data.append(Observation(t, order, truth(seed, t)[order]+noise, error, t))
    if payload is not None:
        stamp = F(4, 5)
        # A more precise timestamped sample arrives 0.2 seconds after sampling.
        data.append(Observation(stamp, payload, truth(seed, stamp)[payload], F(1, 1000), F(1)))
    return data


def derivative_lp(samples, times, noise, jerk, a_radius=F(1)):
    """Haimovich et al. (7)-(9), augmented by the SAME initial v/a prior."""
    N = len(times); n = 2*N
    A, b = [], []
    def add(values, bound):
        row = [F(0)]*n
        for i, v in values.items(): row[i] = F(v)
        A.append(row); b.append(F(bound))
    for k in range(1, N):
        h = times[k]-times[k-1]
        for sign in (-1, 1):
            add({N+k: sign, N+k-1: -sign}, jerk*h)
            add({k-1: sign, k: -sign, N+k: sign*h}, jerk*h*h/2)
    for k, sample in enumerate(samples):
        add({k: 1}, sample+noise); add({k: -1}, -sample+noise)
    for sign in (-1, 1):
        add({0: sign}, 0); add({N: sign}, a_radius)
    output = []
    for sign in (-1, 1):
        c = [F(0)]*n; c[-1] = F(sign)
        with warnings.catch_warnings(record=True):
            result = linprog(-np.array(list(map(float, c))),
                             A_ub=np.array([[float(v) for v in row] for row in A]),
                             b_ub=np.array(list(map(float, b))), bounds=[(None, None)]*n,
                             method='highs', options=dict(threads=1))
        assert result.status == 0
        primal = [F(float(x)).limit_denominator(10**8) for x in result.x]
        dual = [F(float(-x)).limit_denominator(10**8) for x in result.ineqlin.marginals]
        output.append(dict(A=[list(map(str, row)) for row in A], b=list(map(str, b)),
                           primal=list(map(str, primal)), dual=list(map(str, dual)),
                           objective=list(map(str, c)), value=str(dot(c, primal)), sign=sign))
    return output


def projection(constraints, nominal):
    A = [tuple(map(F, c['normal'])) for c in constraints]+[(F(1), F(0)), (F(-1), F(0)),
                                                         (F(0), F(1)), (F(0), F(-1))]
    b = [F(c['rhs']) for c in constraints]+[F(1)]*4
    feasible = lambda x: all(dot(a, x) <= rhs for a, rhs in zip(A, b))
    candidates = []
    if feasible(nominal):
        return nominal, [F(0)]*len(A)
    for i, (a, rhs) in enumerate(zip(A, b)):
        weight = (dot(a, nominal)-rhs)/dot(a, a)
        x = tuple(v-weight*w for v, w in zip(nominal, a))
        if weight >= 0 and feasible(x):
            mu = [F(0)]*len(A); mu[i] = weight
            candidates.append((sum((x[j]-nominal[j])**2 for j in range(2)), x, mu))
    for i, j in combinations(range(len(A)), 2):
        a, d = A[i], A[j]
        det = a[0]*d[1]-a[1]*d[0]
        if not det: continue
        x = ((b[i]*d[1]-a[1]*b[j])/det, (a[0]*b[j]-b[i]*d[0])/det)
        v = [nominal[k]-x[k] for k in range(2)]
        m1 = (v[0]*d[1]-d[0]*v[1])/det
        m2 = (a[0]*v[1]-v[0]*a[1])/det
        if m1 >= 0 and m2 >= 0 and feasible(x):
            mu = [F(0)]*len(A); mu[i] = m1; mu[j] = m2
            candidates.append((sum((x[j]-nominal[j])**2 for j in range(2)), x, mu))
    if not candidates:
        return None
    _, x, mu = min(candidates, key=lambda item: item[0])
    return x, mu


def run():
    supports, baselines, comparisons, plans = [], [], [], []
    def insert(label, center, radii, jerk, obs, target, weights):
        record = solve(center, radii, jerk, obs, F(1), target, weights, tolerance=F(1, 100000))
        assert record['converged'], (label, record['gap'])
        index = len(supports)
        supports.append(dict(label=label, certificate=record))
        return index
    # Same-data continuous derivative comparison, including an extremal-curvature history.
    times = [F(i, 4) for i in range(5)]
    for seed in range(4):
        noise = F(1, 100)
        samples = [(t*t/2 if seed == 0 else truth(seed, t)[1]) for t in times]
        obs = [Observation(t, 1, v, noise, t) for t, v in zip(times, samples)]
        ids = [insert('derivative_'+str(seed)+'_'+str(sign), [F(0)]*3, [F(0), F(0), F(1)],
                      F(1), obs, F(1), (F(0), F(0), F(sign))) for sign in (-1, 1)]
        lp = derivative_lp(samples, times, noise, F(1))
        for item in lp: baselines.append(dict(seed=seed, certificate=item))
        new_width = sum((F(supports[i]['certificate']['dual']['upper']) for i in ids), F(0))
        old_width = sum((F(item['value']) for item in lp), F(0))
        comparisons.append(dict(seed=seed, support_indices=ids, certified_width=str(new_width),
                                lp_width=str(old_width)))
        print('derivative', seed, float(new_width), float(old_width), flush=True)
    # Three neighbors, continuous radial motion, one controlled planar ego.
    normals = [(F(1), F(0)), (F(-3, 5), F(4, 5)), (F(-3, 5), F(-4, 5))]
    for payload in (None, 0, 1, 2):
        constraints = []
        starts = []
        for seed, normal in enumerate(normals):
            obs = observations(seed, payload)
            offset = F(4, 5)-truth(seed, F(1))[0]
            start = insert('start_'+str(payload)+'_'+str(seed),
                           [F(0)]*3, [F(1, 50), F(3, 100), F(1, 2)], F(1, 2),
                           obs, F(1), (F(-1), F(0), F(0)))
            starts.append(dict(normal=list(map(str, normal)), offset=str(offset), support_index=start))
            for h in (F(1, 4), F(1, 2)):
                index = insert('plan_'+str(payload)+'_'+str(seed)+'_'+str(h),
                               [F(0)]*3, [F(1, 50), F(3, 100), F(1, 2)], F(1, 2),
                               obs, F(1)+h, (F(-1), F(0), F(0)))
                upper = F(supports[index]['certificate']['dual']['upper'])
                # For t<=1.5: |a_neighbor|<=.5+.5*1.5=1.25. Ego |u_axis|<=1.
                M = F(5, 4)+sum(map(abs, normal), F(0))
                buffer = M*F(1, 4)**2/8
                rhs = 2*(offset-upper-F(3, 5)-buffer)/h**2
                constraints.append(dict(normal=list(map(str, normal)), rhs=str(rhs),
                                        support_index=index, offset=str(offset), distance='3/5',
                                        buffer=str(buffer), acceleration_bound=str(M), step='1/4'))
        nominal = (F(1), F(2, 5))
        result = projection(constraints, nominal)
        assert result is not None, ('one-hold QP infeasible', payload)
        control, mu = result
        cost = sum(((control[i]-nominal[i])**2 for i in range(2)), F(0))/2
        plans.append(dict(payload=payload, fixed_reply_packets=0 if payload is None else 3,
                          nominal=list(map(str, nominal)), control=list(map(str, control)),
                          multipliers=list(map(str, mu)), constraints=constraints, starts=starts,
                          intervention=str(cost)))
        print('plan', payload, list(map(float, control)), float(cost), flush=True)
    truths = []
    for seed in range(3):
        signs = -1 if seed % 2 else 1
        truths.append(dict(initial=['0', '0', str(F(seed-1, 20))],
                           mesh=['0', '3/10', '11/20', '4/5', '7/5', '2'],
                           jerks=list(map(str, [signs*j for j in [F(1, 4), F(-1, 3), F(1, 5), F(-1, 4), F(1, 3)]]))))
    return dict(supports=supports, derivative_lp=baselines, comparisons=comparisons,
                plans=plans, physical_truths=truths)
