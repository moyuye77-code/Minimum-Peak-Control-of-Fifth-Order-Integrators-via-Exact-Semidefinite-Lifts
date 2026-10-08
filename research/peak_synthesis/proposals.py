import time
import warnings
import numpy as np
from scipy.optimize import minimize,LinearConstraint
from research.classic_cone_audit.compiler import cp
from research.four_history_audit.proposals import numeric_integral
from .model import homogeneous,knots,adjusted,local_kernel,kernel


def build_sdp(case):
    case.validate(); times=knots(case); n=len(times)-1
    J=cp.Variable(name='input_peak'); y=cp.Variable((case.axes*n,4)); aux=cp.Variable((case.axes*n,6))
    cons=[J>=0]; pencils=[]
    for i in range(case.axes*n):
        blocks=[cp.bmat(m) for m in homogeneous(J,y[i,:],aux[i,:])]
        pencils.extend(blocks); cons.extend(g>>0 for g in blocks)
    rhs=adjusted(case); handles=[]
    for o,b in zip(case.waypoints,rhs):
        value=0
        for a in range(case.axes):
            for k,(l,r) in enumerate(zip(times,times[1:])):
                kk=local_kernel(o,a,l,r)
                value+=sum(float(kk[j])*(2*y[a*n+k,j]-J/(j+1)) for j in range(4))
        if o.error==0:
            c=value==float(b); cons.append(c); handles.append((c,None))
        else:
            upper=value<=float(b+o.error); lower=value>=float(b-o.error)
            cons.extend((upper,lower)); handles.append((upper,lower))
    problem=cp.Problem(cp.Minimize(J),cons)
    assert problem.is_dcp()
    return problem,J,y,aux,pencils,handles


def sdp(case):
    start=time.perf_counter()
    problem,J,y,aux,pencils,handles=build_sdp(case)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            value=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,
                                tol_feas=1e-10,max_iter=300,warm_start=False)
            failure=None
        except cp.error.SolverError as exc:
            value=None; failure=str(exc)
    out=dict(method='homogeneous_sdp',status='solver_error' if failure else problem.status,
             failure=failure,warnings=[str(w.message) for w in caught],optimization_calls=1,
             elapsed_seconds=time.perf_counter()-start,segments=len(knots(case))-1)
    if value is not None and np.isfinite(value) and all(h.dual_value is not None for h,_ in handles):
        out.update(value=float(value),lam=[float(-h.dual_value+(lo.dual_value if lo is not None else 0))
                                          for h,lo in handles],
                   minimum_eigenvalue=min(float(np.linalg.eigvalsh(g.value).min()) for g in pencils),
                   moments=y.value.tolist(),auxiliaries=aux.value.tolist())
    return out


def continuous_dual(case):
    case.validate(); start=time.perf_counter(); times=knots(case)
    n=len(case.waypoints); rhs=np.array(list(map(float,adjusted(case))))
    eps=np.array([float(o.error) for o in case.waypoints]); noisy=np.flatnonzero(eps>0)
    size=n+len(noisy); pieces=[]
    for a in range(case.axes):
        for l,r in zip(times,times[1:]):
            kernels=np.array([[float(v) for v in kernel(o,a)] if r<=o.time else [0.]*4
                              for o in case.waypoints])
            pieces.append((float(l),float(r),kernels))
    def normgrad(lam):
        value=0.; gradient=np.zeros(n)
        for l,r,ks in pieces:
            v,g=numeric_integral(lam@ks,l,r,ks)
            value+=v; gradient+=g
        return value,gradient
    def objective(z):
        return -rhs@z[:n]+eps[noisy]@z[n:]
    gradient=np.concatenate((-rhs,eps[noisy]))
    def constraint(z): return 1-normgrad(z[:n])[0]
    def jac(z): return np.concatenate((-normgrad(z[:n])[1],np.zeros(len(noisy))))
    cons=[dict(type='ineq',fun=constraint,jac=jac)]
    if len(noisy):
        rows=[]
        for j,i in enumerate(noisy):
            for sign in (-1,1):
                row=np.zeros(size); row[i]=sign; row[n+j]=1; rows.append(row)
        cons.append(LinearConstraint(np.array(rows),0,np.inf))
    z0=np.zeros(size)
    nn=normgrad(rhs)[0]
    if nn>0:
        z0[:n]=rhs/(2*nn); z0[n:]=np.abs(z0[:n][noisy])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result=minimize(objective,z0,jac=lambda z:gradient,method='SLSQP',
                        constraints=cons,bounds=[(None,None)]*n+[(0,None)]*len(noisy),
                        options=dict(maxiter=600,ftol=1e-11))
    return dict(method='continuous_root_dual',status=str(result.message),success=bool(result.success),
                value=float(-result.fun),lam=result.x[:n].tolist(),iterations=int(result.nit),
                numeric_dual_norm=float(normgrad(result.x[:n])[0]),optimization_calls=1,
                warnings=[str(w.message) for w in caught],elapsed_seconds=time.perf_counter()-start)
