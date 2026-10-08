"""Shared rational lower-witness / continuous upper-bound procedure."""
from fractions import Fraction as F
import time
import warnings
import numpy as np
from scipy.optimize import linprog
from research.continuous_evidence.support import recover_vertex
from .model import dot, trajectory_row, segment_polynomials, dual_certificate
from .proposals import real_roots


def certify(case,proposal,tolerance=F(1,10**6)):
    case.validate()
    started=time.perf_counter()
    if 'lam' not in proposal or not all(np.isfinite(v) for v in proposal['lam']):
        return dict(certified=False,reason='no_finite_multiplier',repair_lp_calls=0,
                    elapsed_seconds=time.perf_counter()-started)
    lam=[F(float(v)).limit_denominator(10**9) for v in proposal['lam']]
    upper=dual_certificate(case,lam)
    mesh={F(0),case.target}|{o.time for o in case.observations}
    for l,r,p in segment_polynomials(case,lam):
        for root in real_roots(p,float(l),float(r)):
            rounded=F(round(root*10**8),10**8)
            if l<rounded<r: mesh.add(rounded)
    mesh=sorted(mesh)
    objective=trajectory_row(case.target,case.weights,mesh)
    A,b=[],[]
    for o in case.observations:
        row=trajectory_row(o.time,o.weights,mesh)
        A.extend((row,[-v for v in row])); b.extend((o.value+o.error,-o.value+o.error))
    bounds=[(c-r,c+r) for c,r in zip(case.center,case.radii)]+[(-case.bound,case.bound)]*(len(mesh)-1)
    with warnings.catch_warnings(record=True) as notices:
        warnings.simplefilter('always')
        result=linprog(-np.array(list(map(float,objective))),A_ub=np.array([[float(v) for v in row] for row in A]),
                       b_ub=np.array(list(map(float,b))),bounds=[tuple(map(float,x)) for x in bounds],
                       method='highs',options=dict(threads=1,primal_feasibility_tolerance=1e-10,
                                                  dual_feasibility_tolerance=1e-10))
    base=dict(case=case.record(),dual=upper,mesh=list(map(str,mesh)),repair_lp_calls=1,
              lp_status=int(result.status),lp_message=result.message,
              notices=[str(w.message) for w in notices],tolerance=str(tolerance))
    if result.status!=0:
        return dict(base,certified=False,reason='inner_lp_failed_not_continuous_infeasibility',
                    elapsed_seconds=time.perf_counter()-started)
    try:
        values=recover_vertex(result,A,b,bounds)
    except (AssertionError,ValueError,ZeroDivisionError) as exc:
        return dict(base,certified=False,reason='rational_primal_recovery_failed',detail=str(exc),
                    elapsed_seconds=time.perf_counter()-started)
    lower=dot(objective,values); gap=F(upper['upper'])-lower
    assert gap>=0, 'A valid primal/dual pair cannot have negative gap'
    return dict(base,certified=True,trajectory=list(map(str,values)),lower=str(lower),
                gap=str(gap),converged=gap<=tolerance,elapsed_seconds=time.perf_counter()-started)
