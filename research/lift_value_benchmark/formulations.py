"""The same physical optimization rows, two exact homogeneous K4 lifts."""
from functools import lru_cache
from fractions import Fraction as F
import time
import warnings
import numpy as np
import sympy as sp
from scipy.optimize import minimize
from research.classic_cone_audit.compiler import cp
from research.four_lift_prior_audit import algebra as generic
from research.four_history_audit.proposals import numeric_integral
from research.peak_synthesis.model import homogeneous,knots,adjusted,local_kernel,kernel

GENERIC_EXPONENTS=tuple(e for e in generic.EXPONENTS if sum(e)>=2)


@lru_cache(maxsize=1)
def templates():
    yy=generic.moment_data()
    variables=(generic.s,generic.t,generic.q,generic.r)+tuple(yy[e] for e in GENERIC_EXPONENTS)
    zero={v:0 for v in variables};out=[]
    for block in generic.generic_matrices(yy):
        rows=[]
        for row in block.tolist():
            converted=[]
            for expression in row:
                assert sp.Poly(expression,*variables).total_degree()<=1
                constant=F(expression.subs(zero))
                terms=tuple((j,F(sp.diff(expression,v))) for j,v in enumerate(variables)
                            if sp.diff(expression,v)!=0)
                converted.append((constant,terms))
            rows.append(converted)
        out.append(rows)
    return tuple(out)


def generic_homogeneous(J,y,aux):
    values=list(y)+list(aux)
    assert len(values)==16
    return [[[constant*J+sum((c*values[j] for j,c in terms),0)
              for constant,terms in row] for row in block] for block in templates()]


def build_sdp(case,method):
    assert method in ('seven_block','generic_polynomial')
    case.validate();times=knots(case);n=len(times)-1
    J=cp.Variable(name='peak');y=cp.Variable((case.axes*n,4))
    aux=cp.Variable((case.axes*n,6 if method=='seven_block' else 12))
    cons=[J>=0];pencils=[]
    for i in range(case.axes*n):
        raw=homogeneous(J,list(y[i,:]),list(aux[i,:])) if method=='seven_block' else generic_homogeneous(J,list(y[i,:]),list(aux[i,:]))
        blocks=[cp.bmat(b) for b in raw];pencils.extend(blocks);cons.extend(b>>0 for b in blocks)
    handles=[]
    for o,b in zip(case.waypoints,adjusted(case)):
        assert o.error==0,'This exact-waypoint suite does not exercise noisy rows'
        value=0
        for axis in range(case.axes):
            for k,(left,right) in enumerate(zip(times,times[1:])):
                kk=local_kernel(o,axis,left,right)
                value+=sum(float(kk[j])*(2*y[axis*n+k,j]-J/(j+1)) for j in range(4))
        row=value==float(b);handles.append(row);cons.append(row)
    problem=cp.Problem(cp.Minimize(J),cons);assert problem.is_dcp()
    return problem,J,y,aux,pencils,handles


def sdp(case,method):
    start=time.perf_counter();problem,J,y,aux,pencils,handles=build_sdp(case,method)
    built=time.perf_counter()
    data,_,_=problem.get_problem_data('CLARABEL');canonical=time.perf_counter()
    dims=data['dims'];failure=None;value=None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            value=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,
                                tol_feas=1e-10,max_iter=300,max_threads=1,warm_start=False)
        except cp.error.SolverError as exc:failure=str(exc)
    solved=time.perf_counter()
    out=dict(method=method,status='solver_error' if failure else problem.status,failure=failure,
             warnings=[str(w.message) for w in caught],optimization_calls=1,
             construction_seconds=built-start,canonicalization_seconds=canonical-built,
             solve_call_seconds=solved-canonical,proposal_seconds=solved-start,
             variables=int(data['A'].shape[1]),matrix_rows=int(data['A'].shape[0]),
             equalities=int(dims.zero),PSD_blocks=list(map(int,dims.psd)),
             SOC_blocks=list(map(int,dims.soc)),
             iterations=None if problem.solver_stats is None else problem.solver_stats.num_iters)
    if value is not None and np.isfinite(value) and all(h.dual_value is not None for h in handles):
        out.update(value=float(value),lam=[float(-h.dual_value) for h in handles],
                   minimum_eigenvalue=min(float(np.linalg.eigvalsh(g.value).min()) for g in pencils),
                   moments=y.value.tolist(),auxiliaries=aux.value.tolist())
    out['proposal_seconds']=time.perf_counter()-start
    return out


def continuous_dual(case):
    start=time.perf_counter();case.validate();times=knots(case);n=len(case.waypoints)
    assert all(o.error==0 for o in case.waypoints)
    rhs=np.array(list(map(float,adjusted(case))));pieces=[]
    for axis in range(case.axes):
        for l,r in zip(times,times[1:]):
            ks=np.array([[float(v) for v in kernel(o,axis)] if r<=o.time else [0.]*4
                         for o in case.waypoints])
            pieces.append((float(l),float(r),ks))
    def normgrad(lam):
        value=0.;grad=np.zeros(n)
        for l,r,ks in pieces:
            v,g=numeric_integral(lam@ks,l,r,ks);value+=v;grad+=g
        return value,grad
    norm=normgrad(rhs)[0];z0=rhs/(2*norm) if norm>0 else np.zeros(n)
    cons=[dict(type='ineq',fun=lambda z:1-normgrad(z)[0],jac=lambda z:-normgrad(z)[1])]
    built=time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result=minimize(lambda z:-rhs@z,z0,jac=lambda z:-rhs,method='SLSQP',
                        constraints=cons,options=dict(maxiter=600,ftol=1e-11))
    solved=time.perf_counter()
    finite=bool(np.isfinite(result.fun) and np.all(np.isfinite(result.x)))
    return dict(method='continuous_root_dual',status=str(result.message),success=bool(result.success),
                value=float(-result.fun) if finite else None,
                lam=result.x.tolist() if finite else [],
                nonfinite_output=None if finite else {'fun':str(result.fun),'x':list(map(str,result.x))},
                iterations=int(result.nit),
                numeric_dual_norm=float(normgrad(result.x)[0]) if finite else None,optimization_calls=1,
                warnings=[str(w.message) for w in caught],construction_seconds=built-start,
                canonicalization_seconds=0.,solve_call_seconds=solved-built,
                proposal_seconds=time.perf_counter()-start,variables=n)
