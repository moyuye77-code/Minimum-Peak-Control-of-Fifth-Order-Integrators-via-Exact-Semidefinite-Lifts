"""Independent rational certificate checker (no optimizer/producer imports)."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from hashlib import sha256
import json

ROOT = Path(__file__).resolve().parent


def dot(a, b):
    assert len(a) == len(b)
    return sum((x*y for x, y in zip(a, b)), F(0))


def state_at(values, mesh, time):
    p, v, a = values[:3]
    for left, right, jerk in zip(mesh, mesh[1:], values[3:]):
        h = max(F(0), min(time, right)-left)
        p, v, a = p+h*v+h*h*a/2+h**3*jerk/6, v+h*a+h*h*jerk/2, a+h*jerk
    return p, v, a


def initial_row(t, weights):
    return [weights[0], weights[0]*t+weights[1],
            weights[0]*t*t/2+weights[1]*t+weights[2]]


def poly(t, weights):
    return [weights[0]*t*t/2+weights[1]*t+weights[2],
            -weights[0]*t-weights[1], weights[0]/2]


def check_support(record):
    center, radii = list(map(F, record['center'])), list(map(F, record['radii']))
    J, now, T = F(record['jerk']), F(record['now']), F(record['target'])
    weights = list(map(F, record['weights']))
    mesh, values = list(map(F, record['mesh'])), list(map(F, record['trajectory']))
    assert len(center) == len(radii) == len(weights) == 3
    assert J > 0 and 0 <= now <= T and T > 0 and min(radii) >= 0
    assert mesh[0] == 0 and mesh[-1] == T and all(l < r for l, r in zip(mesh, mesh[1:]))
    assert len(values) == len(mesh)+2
    assert all(abs(x-c) <= r for x, c, r in zip(values[:3], center, radii))
    assert all(abs(j) <= J for j in values[3:])
    observations = record['observations']
    for o in observations:
        stamp, receive, eps = F(o['time']), F(o['received']), F(o['error'])
        assert 0 <= stamp <= receive <= now and o['order'] in (0, 1, 2) and eps >= 0
        assert abs(state_at(values, mesh, stamp)[o['order']]-F(o['value'])) <= eps
    primal = dot(weights, state_at(values, mesh, T))
    assert primal == F(record['lower'])
    dual = record['dual']
    lam = list(map(F, dual['multipliers']))
    assert len(lam) == len(observations)
    c = initial_row(T, weights)
    bound = F(0)
    for o, l in zip(observations, lam):
        e = [F(i == o['order']) for i in range(3)]
        c = [a-l*b for a, b in zip(c, initial_row(F(o['time']), e))]
        bound += l*F(o['value'])+abs(l)*F(o['error'])
    assert c == list(map(F, dual['residual_initial']))
    bound += dot(c, center)+dot(list(map(abs, c)), radii)
    total, last = F(0), F(0)
    for leaf in dual['leaves']:
        left, right = F(leaf['left']), F(leaf['right'])
        assert left == last and left < right <= T
        assert not any(left < F(o['time']) < right for o in observations)
        expected = poly(T, weights)
        for o, l in zip(observations, lam):
            if right <= F(o['time']):
                e = [F(i == o['order']) for i in range(3)]
                expected = [x-l*y for x, y in zip(expected, poly(F(o['time']), e))]
        assert expected == list(map(F, leaf['polynomial']))
        a, b, c2 = expected
        pleft = a+b*left+c2*left*left
        pright = a+b*right+c2*right*right
        bs = [pleft, pleft+(right-left)*(b+2*c2*left)/2, pright]
        assert bs == list(map(F, leaf['bernstein']))
        upper = (right-left)*sum(map(abs, bs), F(0))/3
        assert upper == F(leaf['upper'])
        total += upper; last = right
    assert last == T and total == F(dual['integral_upper'])
    bound += J*total
    assert bound == F(dual['upper']) and bound >= primal
    gap = bound-primal
    assert gap == F(record['gap'])
    assert record['converged'] == (gap <= F(record['tolerance']))
    return dict(gap=float(gap), primal_intervals=len(mesh)-1, dual_intervals=len(dual['leaves']))


def check_lp(data):
    A = [list(map(F, row)) for row in data['A']]
    b, x, mu, c = [list(map(F, data[k])) for k in ('b', 'primal', 'dual', 'objective')]
    assert len(A) == len(b) == len(mu) and len(x) == len(c) and min(mu) >= 0
    assert all(dot(a, x) <= rhs for a, rhs in zip(A, b))
    assert all(sum((row[j]*w for row, w in zip(A, mu)), F(0)) == c[j] for j in range(len(c)))
    assert dot(c, x) == dot(b, mu) == F(data['value'])


def audit(data, check_sources=True):
    if check_sources:
        files = sorted(ROOT.glob('*.py'))+[ROOT/'THEORY.md', ROOT/'PRIOR_ART.md']
        assert data['source_manifest'] == {p.name: sha256(p.read_bytes()).hexdigest() for p in files}
    totals = []
    assert len(data['supports']) == 44 and len(data['derivative_lp']) == 8
    assert len(data['comparisons']) == 4 and {c['seed'] for c in data['comparisons']} == set(range(4))
    assert len(data['plans']) == 4 and [p['payload'] for p in data['plans']] == [None, 0, 1, 2]
    assert len(data['physical_truths']) == 3
    for entry in data['supports']:
        totals.append(check_support(entry['certificate']))
        assert entry['certificate']['converged']
    for baseline in data['derivative_lp']:
        check_lp(baseline['certificate'])
    for comparison in data['comparisons']:
        new = sum((F(data['supports'][i]['certificate']['dual']['upper'])
                   for i in comparison['support_indices']), F(0))
        old = sum((F(r['certificate']['value']) for r in data['derivative_lp']
                   if r['seed'] == comparison['seed']), F(0))
        assert new == F(comparison['certified_width']) and old == F(comparison['lp_width'])
        # Reconstruct the published Taylor LP, with the explicitly declared same initial prior.
        observations = data['supports'][comparison['support_indices'][0]]['certificate']['observations']
        assert len(comparison['support_indices']) == 2
        for index, sign in zip(comparison['support_indices'], (-1, 1)):
            r = data['supports'][index]['certificate']
            assert r['observations'] == observations
            assert list(map(F, r['center'])) == [0, 0, 0] and list(map(F, r['radii'])) == [0, 0, 1]
            assert list(map(F, r['weights'])) == [0, 0, sign]
            assert F(r['jerk']) == F(r['now']) == F(r['target']) == 1
        assert len(observations) == 5 and all(o['order'] == 1 and F(o['error']) == F(1, 100) for o in observations)
        times = [F(o['time']) for o in observations]
        assert times == [F(k, 4) for k in range(5)]
        A, b = [], []
        def add(entries, rhs):
            row = [F(0)]*10
            for j, value in entries.items(): row[j] = F(value)
            A.append(row); b.append(F(rhs))
        for k in range(1, 5):
            for sign in (-1, 1):
                add({5+k: sign, 4+k: -sign}, F(1, 4))
                add({k-1: sign, k: -sign, 5+k: F(sign, 4)}, F(1, 32))
        for k, o in enumerate(observations):
            add({k: 1}, F(o['value'])+F(1, 100))
            add({k: -1}, -F(o['value'])+F(1, 100))
        for sign in (-1, 1):
            add({0: sign}, 0); add({5: sign}, 1)
        group = [r['certificate'] for r in data['derivative_lp'] if r['seed'] == comparison['seed']]
        assert len(group) == 2 and {r['sign'] for r in group} == {-1, 1}
        for r in group:
            assert A == [list(map(F, row)) for row in r['A']] and b == list(map(F, r['b']))
            assert list(map(F, r['objective'])) == [F(0)]*9+[F(r['sign'])]
    for plan in data.get('plans', []):
        u = list(map(F, plan['control']))
        assert len(u) == 2 and max(map(abs, u)) <= 1
        assert plan['fixed_reply_packets'] == (0 if plan['payload'] is None else 3)
        assert len(plan['starts']) == 3 and len(plan['constraints']) == 6
        assert [list(map(F, s['normal'])) for s in plan['starts']] == [[1, 0], [F(-3, 5), F(4, 5)], [F(-3, 5), F(-4, 5)]]
        for neighbor, start in enumerate(plan['starts']):
            normal = list(map(F, start['normal']))
            assert dot(normal, normal) == 1
            r = data['supports'][start['support_index']]['certificate']
            assert F(r['now']) == F(r['target']) == 1
            assert list(map(F, r['weights'])) == [-1, 0, 0]
            assert list(map(F, r['center'])) == [0, 0, 0]
            assert list(map(F, r['radii'])) == [F(1, 50), F(3, 100), F(1, 2)]
            assert F(r['jerk']) == F(1, 2)
            physical = data['physical_truths'][neighbor]
            mesh = list(map(F, physical['mesh']))
            values = list(map(F, physical['initial']+physical['jerks']))
            assert all(abs(j) <= F(r['jerk']) for j in values[3:])
            assert all(abs(v-F(c)) <= F(rad) for v, c, rad in zip(values[:3], r['center'], r['radii']))
            assert F(start['offset'])+state_at(values, mesh, F(1))[0] == F(4, 5)
            for o in r['observations']:
                assert abs(state_at(values, mesh, F(o['time']))[o['order']]-F(o['value'])) <= F(o['error'])
            base = r['observations'][:12]
            assert len(base) == 12
            assert [(F(o['time']), o['order'], F(o['error']), F(o['received'])) for o in base] == [
                (F(k, 5), order, F(1, 25), F(k, 5)) for k in range(6) for order in (0, 1)]
            assert len(r['observations']) == (12 if plan['payload'] is None else 13)
            if plan['payload'] is not None:
                packet = r['observations'][-1]
                assert packet['order'] == plan['payload'] and F(packet['time']) == F(4, 5)
                assert F(packet['received']) == 1 and F(packet['error']) == F(1, 1000)
            M = F(5, 4)+sum(map(abs, normal), F(0))
            assert F(start['offset'])-F(r['dual']['upper']) >= F(3, 5)+M*F(1, 4)**2/8
        for check in plan['constraints']:
            assert dot(list(map(F, check['normal'])), u) <= F(check['rhs'])
            bound = data['supports'][check['support_index']]['certificate']
            assert list(map(F, bound['weights'])) == [F(-1), F(0), F(0)]
            horizon = F(bound['target'])-F(bound['now'])
            # n'(p_j - p_i) >= d at samples; continuous interpolation buffer below.
            low = F(check['offset'])-F(bound['dual']['upper'])
            assert low-F(check['distance'])-F(check['buffer']) == F(check['rhs'])*horizon*horizon/2
            assert F(check['buffer']) == F(check['acceleration_bound'])*F(check['step'])**2/8
            assert F(check['acceleration_bound']) == F(5, 4)+sum(map(abs, map(F, check['normal'])), F(0))
            matching = [s for s in plan['starts'] if s['normal'] == check['normal']]
            assert len(matching) == 1 and matching[0]['offset'] == check['offset']
            initial = data['supports'][matching[0]['support_index']]['certificate']
            assert bound['observations'] == initial['observations']
            assert bound['center'] == initial['center'] and bound['radii'] == initial['radii']
            assert F(bound['jerk']) == F(1, 2) and F(bound['now']) == 1
        for start in plan['starts']:
            targets = sorted(F(data['supports'][item['support_index']]['certificate']['target'])
                             for item in plan['constraints'] if item['normal'] == start['normal'])
            assert targets == [F(5, 4), F(3, 2)]
        # Verify exact Euclidean projection via KKT multipliers over ALL halfplanes.
        nominal = list(map(F, plan['nominal']))
        grad = [u[i]-nominal[i] for i in range(2)]
        A = [list(map(F, item['normal'])) for item in plan['constraints']]
        b = [F(item['rhs']) for item in plan['constraints']]
        A += [[F(1), F(0)], [F(-1), F(0)], [F(0), F(1)], [F(0), F(-1)]]
        b += [F(1)]*4
        mu = list(map(F, plan['multipliers']))
        assert len(mu) == len(A) and min(mu) >= 0
        for a, rhs, w in zip(A, b, mu):
            assert dot(a, u) <= rhs and w*(rhs-dot(a, u)) == 0
            grad = [g+w*v for g, v in zip(grad, a)]
        assert grad == [0, 0]
        assert sum(((u[i]-nominal[i])**2 for i in range(2)), F(0))/2 == F(plan['intervention'])
    return dict(status='passed', continuous_support_certificates=len(totals),
                exact_lp_baseline_certificates=len(data['derivative_lp']),
                exact_plan_projection_certificates=len(data.get('plans', [])),
                max_support_gap=max(x['gap'] for x in totals),
                primal_jerk_intervals=sum(x['primal_intervals'] for x in totals),
                dual_polynomial_intervals=sum(x['dual_intervals'] for x in totals))


def main():
    path = ROOT/'results'/'verification.json'
    result = audit(json.loads(path.read_text(encoding='utf-8')))
    result['input_sha256'] = sha256(path.read_bytes()).hexdigest()
    (path.parent/'independent_audit.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
