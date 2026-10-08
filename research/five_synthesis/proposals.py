"""Two continuous optimization proposals; neither supplies a primal certificate."""
import time
import warnings
import numpy as np
from scipy.optimize import minimize,LinearConstraint
from research.classic_cone_audit.compiler import cp
from research.four_history_audit.proposals import real_roots
from .model import homogeneous,knots,adjusted,local_kernel,kernel


def numeric_integral(p,left,right,other):
    cuts=[left]+real_roots(p,left,right)+[right]
    total=0.;grad=np.zeros(len(other));powers=np.arange(1,6)
    for lo,hi in zip(cuts,cuts[1:]):
        sign=float(np.sign(np.polynomial.polynomial.polyval((lo+hi)/2,p)))
        primitive=(hi**powers-lo**powers)/powers
        total+=sign*float(np.dot(p,primitive));grad+=sign*(other@primitive)
    return total,grad


def build_sdp(case):
    case.validate();times=knots(case);n=len(times)-1
    J=cp.Variable(name='peak');Y=cp.Variable((case.axes*n,5));aux=cp.Variable((case.axes*n,34))
    cons=[J>=0];pencils=[];slacks=[]
    for i in range(case.axes*n):
        raw,ss=homogeneous(J,list(Y[i,:]),list(aux[i,:]))
        gs=[cp.bmat(g) for g in raw];pencils.extend(gs);slacks.extend(ss)
        cons.extend(g>>0 for g in gs);cons.extend(v>=0 for v in ss)
    handles=[]
    for row,b in zip(case.waypoints,adjusted(case)):
        value=0
        for axis in range(case.axes):
            for k,(left,right) in enumerate(zip(times,times[1:])):
                kk=local_kernel(row,axis,left,right)
                value+=sum(float(kk[j])*(2*Y[axis*n+k,j]-J/(j+1)) for j in range(5))
        if row.error==0:
            con=value==float(b);cons.append(con);handles.append((con,None))
        else:
            up=value<=float(b+row.error);lo=value>=float(b-row.error)
            cons.extend((up,lo));handles.append((up,lo))
    problem=cp.Problem(cp.Minimize(J),cons);assert problem.is_dcp()
    return problem,J,Y,aux,pencils,slacks,handles


def sdp(case):
    started=time.perf_counter();problem,J,Y,aux,pencils,slacks,handles=build_sdp(case)
    built=time.perf_counter();data,_,_=problem.get_problem_data('CLARABEL');compiled=time.perf_counter()
    failure=None;value=None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            value=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,
                                tol_feas=1e-10,max_iter=300,max_threads=1,warm_start=False)
        except cp.error.SolverError as exc:failure=str(exc)
    solved=time.perf_counter()
    out=dict(method='closed_five_sdp',status='solver_error' if failure else problem.status,failure=failure,
             warnings=[str(w.message) for w in caught],optimization_calls=1,
             construction_seconds=built-started,canonicalization_seconds=compiled-built,
             solve_call_seconds=solved-compiled,elapsed_seconds=solved-started,
             variables=int(data['A'].shape[1]),PSD_blocks=list(map(int,data['dims'].psd)))
    if value is not None and np.isfinite(value) and Y.value is not None and aux.value is not None and np.all(
        np.isfinite(Y.value)) and np.all(np.isfinite(aux.value)) and all(
        up.dual_value is not None and (lo is None or lo.dual_value is not None) for up,lo in handles):
        out.update(value=float(value),lam=[float(-up.dual_value+(lo.dual_value if lo is not None else 0)) for up,lo in handles],
                   moments=Y.value.tolist(),auxiliaries=aux.value.tolist(),
                   minimum_eigenvalue=min(float(np.linalg.eigvalsh(g.value).min()) for g in pencils),
                   minimum_scalar_slack=min(float(v.value) for v in slacks))
    return out


def continuous_dual(case):
    case.validate();started=time.perf_counter();times=knots(case);n=len(case.waypoints)
    rhs=np.array(list(map(float,adjusted(case))));eps=np.array([float(row.error) for row in case.waypoints])
    noisy=np.flatnonzero(eps>0);size=n+len(noisy);pieces=[]
    for axis in range(case.axes):
        for left,right in zip(times,times[1:]):
            ks=np.array([[float(v) for v in kernel(row,axis)] if right<=row.time else [0.]*5 for row in case.waypoints])
            pieces.append((float(left),float(right),ks))
    def normgrad(lam):
        value=0.;gradient=np.zeros(n)
        for left,right,ks in pieces:
            v,g=numeric_integral(lam@ks,left,right,ks);value+=v;gradient+=g
        return value,gradient
    objective=lambda z:-rhs@z[:n]+eps[noisy]@z[n:]
    gradient=np.concatenate((-rhs,eps[noisy]))
    cons=[dict(type='ineq',fun=lambda z:1-normgrad(z[:n])[0],
               jac=lambda z:np.concatenate((-normgrad(z[:n])[1],np.zeros(len(noisy)))))]
    if len(noisy):
        rows=[]
        for j,i in enumerate(noisy):
            for sign in (-1,1):
                row=np.zeros(size);row[i]=sign;row[n+j]=1;rows.append(row)
        cons.append(LinearConstraint(np.array(rows),0,np.inf))
    z0=np.zeros(size);normal=normgrad(rhs)[0]
    if normal>0:z0[:n]=rhs/(2*normal);z0[n:]=np.abs(z0[:n][noisy])
    built=time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result=minimize(objective,z0,jac=lambda z:gradient,method='SLSQP',constraints=cons,
                        bounds=[(None,None)]*n+[(0,None)]*len(noisy),options=dict(maxiter=600,ftol=1e-11))
    solved=time.perf_counter();finite=bool(np.isfinite(result.fun) and np.all(np.isfinite(result.x)))
    return dict(method='continuous_root_dual',status=str(result.message),success=bool(result.success),
                value=float(-result.fun) if finite else None,lam=result.x[:n].tolist() if finite else [],
                iterations=int(result.nit),numeric_dual_norm=float(normgrad(result.x[:n])[0]) if finite else None,
                optimization_calls=1,warnings=[str(w.message) for w in caught],
                construction_seconds=built-started,canonicalization_seconds=0.,solve_call_seconds=solved-built,
                elapsed_seconds=solved-started)
