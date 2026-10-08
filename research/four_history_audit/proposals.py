"""Both numerical proposals solve the continuous problem, not a time-grid outer set."""
from fractions import Fraction as F
from math import factorial, comb
import time
import warnings
import numpy as np
from scipy.optimize import minimize, LinearConstraint
from research.classic_cone_audit.compiler import cp
from research.moment_four_lift.exact import matrices
from .model import initial_row, kernel_poly


def sdp(case):
    case.validate()
    started=time.perf_counter()
    knots=sorted({F(0),case.target}|{o.time for o in case.observations})
    n=len(knots)-1
    X=cp.Variable((n+1,4))
    M=cp.Variable((n,4)); A=cp.Variable((n,6))
    cons=[]; gram=[]
    for i,(c,r) in enumerate(zip(case.center,case.radii)):
        if r==0: cons.append(X[0,i]==float(c))
        else: cons.extend([X[0,i]<=float(c+r),X[0,i]>=float(c-r)])
    for k,(l,r) in enumerate(zip(knots,knots[1:])):
        pencils=[cp.bmat(mat) for mat in matrices(M[k,:],A[k,:])]
        gram.extend(pencils); cons.extend(mat>>0 for mat in pencils)
        h=float(r-l)
        for i in range(4):
            degree=3-i
            predicted=sum(X[k,j]*h**(j-i)/factorial(j-i) for j in range(i,4))
            predicted+=float(case.bound)*h**(degree+1)/factorial(degree)*(
                2*sum((-1)**j*comb(degree,j)*M[k,j] for j in range(degree+1))-1/(degree+1))
            cons.append(X[k+1,i]==predicted)
    handles=[]
    for o in case.observations:
        value=np.array(list(map(float,o.weights)))@X[knots.index(o.time),:]
        if o.error==0:
            con=value==float(o.value); cons.append(con); handles.append((con,None))
        else:
            up=value<=float(o.value+o.error); lo=value>=float(o.value-o.error)
            cons.extend((up,lo)); handles.append((up,lo))
    problem=cp.Problem(cp.Maximize(np.array(list(map(float,case.weights)))@X[-1,:]),cons)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            value=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,
                                tol_feas=1e-10,max_iter=300,warm_start=False)
            failure=None
        except cp.error.SolverError as exc:
            value=None; failure=str(exc)
    report=dict(method='joint_sdp',status='solver_error' if failure else problem.status,
                failure=failure,elapsed_seconds=time.perf_counter()-started,
                warnings=[str(w.message) for w in caught],segments=n,
                optimization_calls=1,variables=4*(n+1)+10*n,PSD_blocks=7*n)
    if failure is None and value is not None and np.isfinite(value) and all(up.dual_value is not None for up,lo in handles):
        report.update(value=float(value),lam=[float(up.dual_value-(lo.dual_value if lo is not None else 0)) for up,lo in handles],
                      states=X.value.tolist(),moments=M.value.tolist(),
                      minimum_eigenvalue=min(float(np.linalg.eigvalsh(g.value).min()) for g in gram))
    return report


def real_roots(p,l,r):
    p=np.trim_zeros(np.array(p,dtype=float),'b')
    if len(p)<=1: return []
    roots=np.polynomial.polynomial.polyroots(p)
    return sorted({float(z.real) for z in roots if abs(z.imag)<1e-8 and l<z.real<r})


def numeric_integral(p,l,r,other=None):
    """Analytic antiderivatives on floating-point root intervals; NOT certified."""
    knots=[l]+real_roots(p,l,r)+[r]
    total=0.0
    integrals=np.zeros(0 if other is None else len(other))
    anti=np.array(p)/np.arange(1,5)
    for left,right in zip(knots,knots[1:]):
        middle=(left+right)/2
        sign=float(np.sign(np.polynomial.polynomial.polyval(middle,p)))
        powers=(right**np.arange(1,5)-left**np.arange(1,5))
        total+=sign*float(anti@powers)
        if other is not None: integrals+=sign*(other@(powers/np.arange(1,5)))
    return total,integrals


def root_dual(case):
    case.validate()
    started=time.perf_counter()
    obs=case.observations; n=len(obs)
    base=np.array(list(map(float,initial_row(case.target,case.weights))))
    initial=np.array([[float(v) for v in initial_row(o.time,o.weights)] for o in obs]).T
    center=np.array(list(map(float,case.center))); radii=np.array(list(map(float,case.radii)))
    values=np.array([float(o.value) for o in obs]); eps=np.array([float(o.error) for o in obs])
    target=np.array(list(map(float,kernel_poly(case.target,case.weights))))
    knots=sorted({F(0),case.target}|{o.time for o in obs})
    pieces=[]
    for l,r in zip(knots,knots[1:]):
        kernels=np.array([[float(v) for v in kernel_poly(o.time,o.weights)] if r<=o.time else [0.]*4 for o in obs])
        pieces.append((float(l),float(r),kernels))
    uncertain=np.flatnonzero(radii>0); size=2*n+len(uncertain)
    calls=0
    def fun(x):
        nonlocal calls
        calls+=1
        lam=x[:n]-x[n:2*n]
        res=base-initial@lam
        val=float(values@lam+eps@(x[:n]+x[n:2*n])+center@res+radii[uncertain]@x[2*n:])
        grad_lam=values-initial.T@center
        for l,r,kernels in pieces:
            integral,gradient=numeric_integral(target-lam@kernels,l,r,kernels)
            val+=float(case.bound)*integral
            grad_lam-=float(case.bound)*gradient
        grad=np.concatenate((grad_lam+eps,-grad_lam+eps,radii[uncertain]))
        return val,grad
    constraints=[]
    if len(uncertain):
        rows=[]; rhs=[]
        for j,idx in enumerate(uncertain):
            for sign in (1,-1):
                row=np.zeros(size)
                row[:n]=sign*initial[idx]; row[n:2*n]=-sign*initial[idx]; row[2*n+j]=1
                rows.append(row); rhs.append(sign*base[idx])
        constraints=[LinearConstraint(np.array(rows),np.array(rhs),np.inf)]
    x0=np.zeros(size); x0[2*n:]=np.abs(base[uncertain])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result=minimize(fun,x0,method='SLSQP',jac=True,bounds=[(0,None)]*size,
                        constraints=constraints,options=dict(maxiter=400,ftol=1e-11))
    return dict(method='continuous_root_dual',status=str(result.message),success=bool(result.success),
                elapsed_seconds=time.perf_counter()-started,iterations=int(result.nit),
                function_evaluations=calls,value=float(result.fun),
                lam=(result.x[:n]-result.x[n:2*n]).tolist(),optimization_calls=1,
                warnings=[str(w.message) for w in caught],segments=len(pieces))
