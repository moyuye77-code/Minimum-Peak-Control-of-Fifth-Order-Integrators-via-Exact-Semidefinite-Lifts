"""Fresh, non-overwriting experiments tied to the manuscript's exact claims."""
import argparse
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
from math import comb, factorial
import hashlib
import json
import platform
import sys
import time
import warnings

import numpy as np
import sympy as sp
from research.classic_cone_audit.compiler import cp
from research.closed_five_lift.solver import constraints
from research.moment_four_lift.exact import matrices as prefix_matrices, atomic_aux
from research.moment_upper_gate.exact import psd_pivots
from research.paper.check_lift_coefficients import polynomial_data, parse_pencils
from research.five_synthesis.model import cases, knots, homogeneous
from research.five_synthesis.proposals import sdp
from research.five_synthesis.certify import certify
from research.five_synthesis.check import check, check_lower, reference_scale, rhs_for

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MANUSCRIPT = ROOT / 'research/paper/representability.tex'
CENTER = tuple(F(1, 2*(j+1)) for j in range(5))
EPSILONS = tuple(map(F, ('1', '1/10', '1/100', '1/10000', '1/1000000',
                              '1/100000000', '1/10000000000', '0')))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path, data):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def complement(point):
    return tuple(F(1, j+1)-v for j, v in enumerate(point))


def interval_point(intervals):
    return tuple(sum((F(b)**(j+1)-F(a)**(j+1) for a, b in intervals), F(0))/(j+1)
                 for j in range(5))


def endpoints():
    return (('zero', interval_point(())), ('full', interval_point(((0, 1),))),
            ('left', interval_point(((0, F(2, 5)),))),
            ('right', interval_point(((F(3, 5), 1),))),
            ('central', interval_point(((F(1, 4), F(3, 4)),))),
            ('two_intervals', interval_point(((0, F(1, 4)), (F(1, 2), F(3, 4))))))


@lru_cache(None)
def path_functions(endpoint):
    x, _, _, p, _, basis, _ = polynomial_data()
    e = sp.Symbol('epsilon')
    substitution = dict(zip(x, ((1-e)*sp.Rational(a)+e*sp.Rational(b)
                               for a, b in zip(endpoint, CENTER))))
    denominator = sp.expand(p.subs(substitution))
    values = tuple(sp.cancel(b.subs(substitution)/denominator) for b in basis)
    limits = tuple(F(v.subs(e, 0)) for v in values)
    return e, denominator, values, limits


def path_record(name, endpoint, epsilon):
    point = tuple((1-epsilon)*a+epsilon*b for a, b in zip(endpoint, CENTER))
    sides = []
    for end in (endpoint, complement(endpoint)):
        e, denominator, functions, limits = path_functions(end)
        values = tuple(F(v.subs(e, sp.Rational(epsilon))) for v in functions)
        den = F(denominator.subs(e, sp.Rational(epsilon)))
        sides.append(dict(values=list(map(str, values)), limit=list(map(str, limits)),
                          denominator=str(den), raw_inverse=None if not den else str(1/den),
                          norm=str(max(map(abs, values))),
                          limit_error=str(max(abs(a-b) for a, b in zip(values, limits)))))
    return dict(name=name, epsilon=str(epsilon), point=list(map(str, point)), sides=sides)


@lru_cache(None)
def manuscript_pencils():
    return parse_pencils(MANUSCRIPT.read_text(encoding='utf-8'))


def check_witness(point, values_by_side):
    theta, matrices = manuscript_pencils()
    for mat in prefix_matrices(tuple(point[:4]), atomic_aux(point[0], point[1])):
        psd_pivots(mat)
    for side, values in zip((point, complement(point)), values_by_side):
        assert tuple(values[:5]) == (F(1), *side[:4])
        assert values[5] <= side[4]
        at = dict(zip(theta, map(sp.Rational, values)))
        for matrix in matrices:
            psd_pivots(matrix.xreplace(at).tolist())
    return True


def direct_witness(point):
    x, _, _, p, _, basis, _ = polynomial_data()
    out = []
    for side in (point, complement(point)):
        at = dict(zip(x, map(sp.Rational, side[:4])))
        den = p.xreplace(at)
        if den > 0:
            out.append(tuple(F(b.xreplace(at)/den) for b in basis))
        else:
            assert den == 0
            out.append(path_functions(side)[3])
    return out


def numeric_residuals(point, candidate):
    if candidate is None:
        return dict(diagnostic_pass=False)
    y, a, lo, hi = map(np.asarray, (candidate['y'], candidate['prefix'],
                                  candidate['lower'], candidate['upper']))
    _, pencils, slacks = constraints(y, a, lo, hi)
    eig = min(float(np.linalg.eigvalsh(g.value).min()) for _, g in pencils)
    slack = min(map(float, slacks))
    error = float(np.max(np.abs(y-np.array(list(map(float, point))))))
    return dict(minimum_eigenvalue=eig, minimum_scalar_slack=slack,
                output_error=error, diagnostic_pass=(eig >= -1e-7 and slack >= -1e-7
                                                     and error <= 1e-7))


def solve_membership(point):
    y = cp.Variable(5)
    a, lo, hi = cp.Variable(6), cp.Variable(14), cp.Variable(14)
    cons, _, _ = constraints(y, a, lo, hi)
    problem = cp.Problem(cp.Minimize(cp.norm_inf(cp.hstack([a, lo, hi]))),
                         cons+[y == np.array(list(map(float, point)))])
    failure = None
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            problem.solve(solver='CLARABEL', tol_gap_abs=1e-10, tol_gap_rel=1e-10,
                          tol_feas=1e-10, max_iter=300, max_threads=1, warm_start=False)
        except cp.error.SolverError as exc:
            failure = str(exc)
    candidate = None
    if all(v.value is not None and np.isfinite(v.value).all() for v in (y, a, lo, hi)):
        candidate = dict(y=y.value.tolist(), prefix=a.value.tolist(),
                         lower=lo.value.tolist(), upper=hi.value.tolist())
    residuals = numeric_residuals(point, candidate)
    residuals['diagnostic_pass'] &= problem.status in ('optimal', 'optimal_inaccurate')
    return dict(status=problem.status, failure=failure, warnings=[str(w.message) for w in caught],
                candidate=candidate, residuals=residuals, elapsed_seconds=time.perf_counter()-started)


def membership_controls():
    for name, point in endpoints():
        if name not in ('zero', 'full', 'central'):
            continue
        c = {'zero': (F(-1), F(0), F(0), F(0), F(0)),
             'full': (F(1), F(0), F(0), F(0), F(0)),
             'central': (F(-3, 16), F(1), F(-1), F(0), F(0))}[name]
        support = sum((a*b for a, b in zip(c, point)), F(0))
        norm2 = sum((a*a for a in c), F(0))
        for eta in map(F, ('1/1000', '1/1000000', '1/1000000000')):
            for kind in ('inside', 'outside'):
                target = (tuple((1-eta)*a+eta*b for a, b in zip(point, CENTER)) if kind == 'inside'
                          else tuple(a+eta*b/norm2 for a, b in zip(point, c)))
                margin = sum((a*b for a, b in zip(c, target)), F(0))-support
                assert (margin == eta) if kind == 'outside' else (margin < 0)
                yield dict(name=name, kind=kind, eta=str(eta), point=list(map(str, target)),
                           direction=list(map(str, c)), support=str(support), exact_violation=str(margin))


def raw_synthesis(case, proposal):
    if not all(k in proposal for k in ('value', 'moments', 'auxiliaries')):
        return dict(diagnostic_pass=False)
    J = proposal['value']
    Y, aux = np.array(proposal['moments']), np.array(proposal['auxiliaries'])
    times = knots(case); n = len(times)-1
    cone_eig = []; cone_slacks = [J]
    for i in range(case.axes*n):
        mats, ss = homogeneous(J, list(Y[i]), list(aux[i]))
        cone_eig.extend(float(np.linalg.eigvalsh(np.asarray(g, dtype=float)).min()) for g in mats)
        cone_slacks.extend(map(float, ss))
    errors = []
    # Direct binomial convolution, without the producer's local_kernel function.
    for row in case.waypoints:
        value = 0.
        for axis in range(case.axes):
            for k, (left, right) in enumerate(zip(times, times[1:])):
                if right > row.time:
                    continue
                width = right-left
                for derivative in range(5):
                    degree = 4-derivative
                    for j in range(degree+1):
                        coef = (row.weights[axis*5+derivative]*comb(degree, j)*
                                (row.time-left)**(degree-j)*(-width)**j*width/factorial(degree))
                        value += float(coef)*(2*Y[axis*n+k, j]-J/(j+1))
        target = float(rhs_for(row.record(), list(case.initial), case.horizon))
        errors.append(max(0., abs(value-target)-float(row.error))/max(1., abs(target), float(row.error)))
    eig, slack, error = min(cone_eig), min(cone_slacks), max(errors)
    threshold = 1e-7*max(1., abs(J))
    return dict(minimum_eigenvalue=eig, minimum_scalar_slack=slack, maximum_scaled_row_error=error,
                cone_tolerance=threshold, diagnostic_pass=(eig >= -threshold and slack >= -threshold
                                                           and error <= 1e-7))


def recovered_moments(case, certificate):
    mesh = list(map(F, certificate['mesh'])); controls = list(map(F, certificate['controls']))
    peak = F(certificate['upper']); count = len(mesh)-1; times = knots(case); out = []
    for axis in range(case.axes):
        for left, right in zip(times, times[1:]):
            width = right-left
            row = []
            for j in range(5):
                integral = sum((controls[axis*count+k]*
                                (((min(right, b)-left)/width)**(j+1)-
                                 ((max(left, a)-left)/width)**(j+1))/(j+1)
                                for k, (a, b) in enumerate(zip(mesh, mesh[1:]))
                                if max(left, a) < min(right, b)), F(0))
                row.append((integral+peak/F(j+1))/2)
            out.append(row)
    return out


def attribution(case, proposal, recovery):
    initial = check_lower(case.record(), recovery['initial_dual'])
    certificate = recovery['certificate']
    out = dict(initial_lower=str(initial), certified=certificate.get('certified', False),
               converged=certificate.get('converged', False))
    if not out['certified']:
        return out
    lower, upper = check(certificate); scale = reference_scale(case.record())
    moments = recovered_moments(case, certificate)
    point_checks = []
    if upper == 0:
        assert all(v == 0 for row in moments for v in row)
    else:
        for row in moments:
            point = tuple(v/upper for v in row)
            # Independent manuscript-basis construction; any physical prefix is allowed.
            point_checks.append(check_witness(point, direct_witness(point)))
    out.update(lower=str(lower), upper=str(upper), lower_improvement=str(lower-initial),
               initial_gap_scaled=str((upper-initial)/scale), final_gap_scaled=str((upper-lower)/scale),
               initial_gap_pass=upper-initial <= F(1,100000)*scale,
               homogeneous_moments=[[str(v) for v in row] for row in moments],
               recovered_lift_checks=len(point_checks), recovered_zero_apex=(upper == 0))
    if 'moments' in proposal:
        discrepancy=max(abs(float(v)-proposal['moments'][i][j])
                        for i,row in enumerate(moments) for j,v in enumerate(row))
        out['raw_to_recovered_moment_distance'] = discrepancy
    return out


def source_hashes():
    paths = list(HERE.glob('*.py'))+[HERE/'PROTOCOL.md']
    for directory in ('closed_five_lift','five_synthesis','schur_sparse_lift',
                      'schur_remainder_recursion','moment_four_lift','moment_upper_gate',
                      'peak_synthesis','four_history_audit','classic_cone_audit','general_congruence_gate'):
        paths.extend((ROOT/'research'/directory).glob('*.py'))
    paths.append(ROOT/'research/paper/check_lift_coefficients.py')
    return {str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sorted(set(paths))}


def summarize(folder):
    paths = [json.loads(p.read_text()) for p in sorted((folder/'boundary').glob('*.json'))]
    controls = [json.loads(p.read_text()) for p in sorted((folder/'membership').glob('*.json'))]
    synthesis = [json.loads(p.read_text()) for p in sorted((folder/'synthesis').glob('*.json'))]
    out = dict(boundary=dict(cases=len(paths), exact_passes=sum(r['exact_pass'] for r in paths),
               numerical_passes=sum(r['numeric']['residuals']['diagnostic_pass'] for r in paths),
               warnings=sum(bool(r['numeric']['warnings']) for r in paths),
               max_selected_norm=max(float(F(v['norm'])) for r in paths for v in r['sides']),
               errors_at_smallest_epsilon={r['name']:max(float(F(v['limit_error'])) for v in r['sides'])
                                         for r in paths if r['epsilon']=='1/10000000000'}),
               membership={}, synthesis={})
    for eta in ('1/1000','1/1000000','1/1000000000'):
        chosen=[r for r in controls if r['eta']==eta]
        out['membership'][eta]={kind:dict(cases=sum(r['kind']==kind for r in chosen),
            numeric_accepts=sum(r['numeric']['residuals']['diagnostic_pass'] for r in chosen if r['kind']==kind))
            for kind in ('inside','outside')}
    out['synthesis']['raw_diagnostic_passes']=sum(r['raw']['diagnostic_pass'] for r in synthesis)
    for route in ('sdp_seed','zero_seed'):
        aa=[r['attribution'][route] for r in synthesis]
        out['synthesis'][route]=dict(cases=len(aa),certified=sum(r['certified'] for r in aa),
            converged=sum(r['converged'] for r in aa),
            initial_gap_passes=sum(r.get('initial_gap_pass',False) for r in aa),
            improved_lower=sum(F(r.get('lower_improvement','0'))>0 for r in aa),
            recovered_lift_checks=sum(r.get('recovered_lift_checks',0) for r in aa),
            zero_apex_checks=sum(r.get('recovered_zero_apex',False) for r in aa))
    out['fresh_sdp_calls']=len(paths)+len(controls)+len(synthesis)
    out['recovery_lp_calls']=sum(v['repair_lp_calls'] for r in synthesis for v in r['recovery'].values())
    return out


def create(folder):
    import scipy, clarabel
    folder.mkdir(parents=True, exist_ok=False)
    for sub in ('boundary','membership','synthesis'):
        (folder/sub).mkdir()
    write_new(folder/'manifest.json',dict(protocol_sha256=sha(HERE/'PROTOCOL.md'),sources=source_hashes(),
        manuscript_sha256_at_start=sha(MANUSCRIPT),python=sys.version,platform=platform.platform(),
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,
                      cvxpy=cp.__version__,clarabel=clarabel.__version__),
        epsilons=list(map(str,EPSILONS)),endpoints=[(name,list(map(str,p))) for name,p in endpoints()],
        synthesis_tasks=[case.record() for case,_ in cases()]))
    for i,(name,end) in enumerate(endpoints()):
        for j,epsilon in enumerate(EPSILONS):
            record=path_record(name,end,epsilon);point=tuple(map(F,record['point']))
            record['exact_pass']=check_witness(point,[tuple(map(F,v['values'])) for v in record['sides']])
            record['numeric']=solve_membership(point)
            write_new(folder/'boundary'/f'{i:02d}_{j:02d}.json',record)
            print('boundary',name,epsilon,record['numeric']['status'],record['numeric']['residuals']['diagnostic_pass'],flush=True)
    for i,record in enumerate(membership_controls()):
        record['numeric']=solve_membership(tuple(map(F,record['point'])))
        write_new(folder/'membership'/f'{i:02d}.json',record)
        print('membership',record['name'],record['kind'],record['eta'],record['numeric']['status'],flush=True)
    for i,(case,known) in enumerate(cases()):
        proposal=sdp(case);recovery={};diagnostics={}
        order=('sdp_seed','zero_seed') if i%2==0 else ('zero_seed','sdp_seed')
        for route in order:
            recovery[route]=certify(case,proposal if route=='sdp_seed' else {'lam':[0.]*len(case.waypoints)})
            diagnostics[route]=attribution(case,proposal,recovery[route])
            if known is not None and diagnostics[route]['certified']:
                assert F(diagnostics[route]['lower'])<=known<=F(diagnostics[route]['upper'])
        record=dict(case=case.record(),known_peak=None if known is None else str(known),proposal=proposal,
                    raw=raw_synthesis(case,proposal),recovery=recovery,attribution=diagnostics)
        write_new(folder/'synthesis'/f'{i:02d}.json',record)
        print('synthesis',case.label,record['raw']['diagnostic_pass'],
              {k:v['converged'] for k,v in diagnostics.items()},flush=True)
    summary=summarize(folder)
    write_new(folder/'summary.json',summary)
    write_new(folder/'hashes.json',{str(p.relative_to(folder)).replace('\\','/'):sha(p)
              for p in sorted(folder.rglob('*.json')) if p.name!='hashes.json'})
    print(json.dumps(summary,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();create(args.output.resolve())


if __name__=='__main__':
    main()
