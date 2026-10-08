"""Shared inner witnesses and continuous dual certificates, not a new algorithm claim."""
from fractions import Fraction as F
from math import comb
import time
import warnings
import numpy as np
import sympy as sp
from scipy.linalg import qr
from scipy.optimize import linprog
from sympy.matrices.exceptions import NonInvertibleMatrixError
from research.four_history_audit.proposals import real_roots
from .model import knots,adjusted,dot,kernel_segments,input_row


def bernstein(p,l,r):
    shifted=[sum((p[j]*comb(j,k)*l**(j-k)*(r-l)**k for j in range(k,4)),F(0)) for k in range(4)]
    return [sum((shifted[k]*F(comb(i,k),comb(3,k)) for k in range(i+1)),F(0)) for i in range(4)]


def dual_bound(case,lam):
    assert len(lam)==len(case.waypoints)
    rhs=adjusted(case)
    numerator=dot(lam,rhs)-sum((abs(v)*o.error for v,o in zip(lam,case.waypoints)),F(0))
    leaves=[]
    def enclose(axis,p,l,r,depth=0):
        bs=bernstein(p,l,r)
        up=(r-l)*sum(map(abs,bs),F(0))/4
        lo=abs((r-l)*sum(bs,F(0))/4)
        if up-lo<=F(1,10**12)*(r-l)/case.horizon or depth>=32:
            leaves.append([axis,str(l),str(r),*map(str,p)])
            return up
        mid=(l+r)/2
        return enclose(axis,p,l,mid,depth+1)+enclose(axis,p,mid,r,depth+1)
    norm=sum((enclose(a,p,l,r) for a,l,r,p in kernel_segments(case,lam)),F(0))
    if norm==0:
        assert numerator<=0, 'This multiplier proves infeasibility, outside the feasible pilot'
        lower=F(0)
    else: lower=max(F(0),numerator/norm)
    return dict(lam=list(map(str,lam)),leaves=leaves,numerator=str(numerator),
                norm_upper=str(norm),lower=str(lower))


def rational_multipliers(values,n):
    if len(values)!=n or not all(np.isfinite(v) for v in values): return [F(0)]*n
    scale=max(map(abs,values),default=0)
    return [F(0)]*n if scale==0 else [F(float(v/scale)).limit_denominator(10**9) for v in values]


def recover_controls(case,mesh,numeric):
    """Eliminate fixed initial state, then solve a small exact active-row system.

    Peak is recomputed after correction; no bound or waypoint is relaxed.
    This may reject a floating proposal and is not a completeness assertion.
    """
    A=[input_row(case,o,mesh) for o in case.waypoints]; b=adjusted(case)
    active=[]; targets=[]
    for row,rhs,o in zip(A,b,case.waypoints):
        value=float(np.array(list(map(float,row)))@numeric)
        if o.error==0:
            active.append(row); targets.append(rhs)
        elif abs(value-float(rhs))>=float(o.error)-1e-8:
            active.append(row); targets.append(rhs+(o.error if value>=float(rhs) else -o.error))
    x=[F(round(float(v)*10**12),10**12) for v in numeric]
    if active:
        mat=np.array([[float(v) for v in row] for row in active])
        scales=np.maximum(np.linalg.norm(mat,axis=1),1e-300)
        norm=mat/scales[:,None]
        _,R,rp=qr(norm.T,pivoting=True,mode='economic')
        diagonal=np.abs(np.diag(R)); rank=int(np.count_nonzero(diagonal>1e-11))
        selected=list(map(int,rp[:rank]))
        _,_,cp=qr(norm[selected,:],pivoting=True,mode='economic')
        pivots=list(map(int,cp[:rank])); free=[j for j in range(len(x)) if j not in pivots]
        exact_A=sp.Matrix([[sp.Rational(active[i][j]) for j in pivots] for i in selected])
        exact_b=sp.Matrix([sp.Rational(targets[i]-sum((active[i][j]*x[j] for j in free),F(0)))
                           for i in selected])
        solved=exact_A.inv()*exact_b
        for j,v in zip(pivots,solved): x[j]=F(int(v.p),int(v.q))
    assert all(abs(dot(row,x)-rhs)<=o.error for row,rhs,o in zip(A,b,case.waypoints))
    return x


def inner_lp(case,mesh):
    n=case.axes*(len(mesh)-1); Aeq=[]; beq=[]; Aub=[]; bub=[]; handles=[]
    for o,b in zip(case.waypoints,adjusted(case)):
        row=[F(0)]+input_row(case,o,mesh)
        scale=max([abs(b),o.error,F(1,10**12)]+list(map(abs,row)))
        numeric=[float(v/scale) for v in row]
        if o.error==0:
            handles.append(('eq',len(Aeq),scale)); Aeq.append(numeric); beq.append(float(b/scale))
        else:
            handles.append(('ineq',len(Aub),scale)); Aub.extend([numeric,[-v for v in numeric]])
            bub.extend([float((b+o.error)/scale),float((-b+o.error)/scale)])
    for j in range(n):
        for sign in (1,-1):
            row=[0.]*(n+1); row[0]=-1.; row[j+1]=sign; Aub.append(row); bub.append(0.)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result=linprog([1.]+[0.]*n,A_ub=np.array(Aub),b_ub=np.array(bub),
                       A_eq=np.array(Aeq) if Aeq else None,b_eq=np.array(beq) if beq else None,
                       bounds=[(0,None)]+[(None,None)]*n,method='highs',
                       options=dict(threads=1,primal_feasibility_tolerance=1e-10,
                                    dual_feasibility_tolerance=1e-10))
    lam=[]
    if result.status==0:
        for kind,i,scale in handles:
            val=result.eqlin.marginals[i] if kind=='eq' else result.ineqlin.marginals[i]-result.ineqlin.marginals[i+1]
            lam.append(float(val)/float(scale))
    return result,lam,[str(w.message) for w in caught]


def certify(case,proposal,max_iterations=4,tolerance=F(1,10**6)):
    started=time.perf_counter()
    lam=rational_multipliers(proposal.get('lam',[]),len(case.waypoints))
    bestdual=dual_bound(case,lam); best=None; attempts=[]
    points=knots(case); mesh=set(points)
    for l,r in zip(points,points[1:]): mesh.update(l+(r-l)*F(j,4) for j in range(1,4))
    for iteration in range(max_iterations):
        for axis,l,r,p in kernel_segments(case,lam):
            for root in real_roots(p,float(l),float(r)):
                point=F(round(root*10**10),10**10)
                if l<point<r: mesh.add(point)
        mesh=sorted(mesh)
        if len(mesh)>257: break
        result,values,notices=inner_lp(case,mesh)
        log=dict(iteration=iteration,cells=len(mesh)-1,status=int(result.status),
                 message=result.message,warnings=notices)
        attempts.append(log)
        if result.status==0:
            lam=rational_multipliers(values,len(case.waypoints))
            dual=dual_bound(case,lam)
            if F(dual['lower'])>F(bestdual['lower']): bestdual=dual
            try:
                controls=recover_controls(case,mesh,result.x[1:])
                peak=max(map(abs,controls))
                if best is None or peak<F(best['upper']):
                    best=dict(case=case.record(),mesh=list(map(str,mesh)),controls=list(map(str,controls)),
                              upper=str(peak),certified=True,tolerance=str(tolerance))
                log['primal_recovered']=True
            except (AssertionError,ValueError,ZeroDivisionError,NonInvertibleMatrixError) as exc:
                log['primal_recovered']=False; log['failure']=type(exc).__name__+': '+str(exc)
            if best is not None:
                gap=F(best['upper'])-F(bestdual['lower']); assert gap>=0
                best.update(dual=bestdual,lower=bestdual['lower'],gap=str(gap),
                            converged=gap<=tolerance*max(F(1),F(best['upper'])))
                from .check import check
                check(best)
                log['gap']=str(gap)
                if best['converged']: break
        # This is only an inner feasible trajectory family, not the exact outer model.
        mesh=set(mesh)|{(l+r)/2 for l,r in zip(mesh,mesh[1:])}
    if best is None:
        best=dict(case=case.record(),certified=False,reason='no_exact_primal',dual=bestdual)
    return dict(certificate=best,attempts=attempts,repair_lp_calls=len(attempts),
                elapsed_seconds=time.perf_counter()-started)
