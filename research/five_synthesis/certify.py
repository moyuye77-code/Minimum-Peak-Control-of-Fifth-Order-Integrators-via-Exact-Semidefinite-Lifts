"""Shared standard recovery, generalized to quartic kernels; no novel-method claim."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import time
import warnings
import numpy as np
import sympy as s
from scipy.linalg import qr
from scipy.optimize import linprog
from sympy.matrices.exceptions import NonInvertibleMatrixError
from research.four_history_audit.proposals import real_roots
from research.peak_synthesis.certify import rational_multipliers
from .model import knots,adjusted,dot,kernel,kernel_segments,input_row


def bernstein(poly,left,right):
    shifted=[sum((poly[j]*comb(j,k)*left**(j-k)*(right-left)**k for j in range(k,5)),F(0)) for k in range(5)]
    return [sum((shifted[k]*F(comb(i,k),comb(4,k)) for k in range(i+1)),F(0)) for i in range(5)]


@lru_cache(maxsize=128)
def reference_scale(case):
    candidates=[]
    for row,b in zip(case.waypoints,adjusted(case)):
        norm=sum((row.time*sum(map(abs,bernstein(kernel(row,a),F(0),row.time)),F(0))/5
                  for a in range(case.axes)),F(0))
        numerator=max(F(0),abs(b)-row.error)
        if norm:candidates.append(numerator/norm)
        else:assert numerator==0,'infeasible zero measurement kernel'
    return max(candidates,default=F(0)) or F(1)


def dual_bound(case,lam):
    assert len(lam)==len(case.waypoints)
    numerator=dot(lam,adjusted(case))-sum((abs(v)*row.error for v,row in zip(lam,case.waypoints)),F(0))
    leaves=[]
    def enclose(axis,poly,left,right,depth=0):
        bs=bernstein(poly,left,right)
        up=(right-left)*sum(map(abs,bs),F(0))/5
        lo=abs((right-left)*sum(bs,F(0))/5)
        if up-lo<=F(1,10**12)*(right-left)/case.horizon or depth>=32:
            leaves.append([axis,str(left),str(right),*map(str,poly)])
            return up
        middle=(left+right)/2
        return enclose(axis,poly,left,middle,depth+1)+enclose(axis,poly,middle,right,depth+1)
    norm=sum((enclose(a,p,left,right) for a,left,right,p in kernel_segments(case,lam)),F(0))
    if norm==0:assert numerator<=0;lower=F(0)
    else:lower=max(F(0),numerator/norm)
    return dict(lam=list(map(str,lam)),leaves=leaves,numerator=str(numerator),norm_upper=str(norm),lower=str(lower))


def recover_controls(case,mesh,numeric):
    A=[input_row(case,row,mesh) for row in case.waypoints];b=adjusted(case)
    active=[];targets=[]
    for row,rhs,obs in zip(A,b,case.waypoints):
        value=float(np.array(list(map(float,row)))@numeric)
        if obs.error==0:active.append(row);targets.append(rhs)
        elif abs(value-float(rhs))>=float(obs.error)-1e-8:
            active.append(row);targets.append(rhs+(obs.error if value>=float(rhs) else -obs.error))
    controls=[F(round(float(v)*10**12),10**12) for v in numeric]
    if active:
        mat=np.array([[float(v) for v in row] for row in active])
        scales=np.maximum(np.linalg.norm(mat,axis=1),1e-300);normal=mat/scales[:,None]
        _,R,perm=qr(normal.T,pivoting=True,mode='economic')
        rank=int(np.count_nonzero(np.abs(np.diag(R))>1e-11));selected=list(map(int,perm[:rank]))
        _,_,columns=qr(normal[selected,:],pivoting=True,mode='economic')
        pivots=list(map(int,columns[:rank]));free=[j for j in range(len(controls)) if j not in pivots]
        exact_A=s.Matrix([[s.Rational(active[i][j]) for j in pivots] for i in selected])
        exact_b=s.Matrix([s.Rational(targets[i]-sum((active[i][j]*controls[j] for j in free),F(0))) for i in selected])
        solved=exact_A.inv()*exact_b
        for j,value in zip(pivots,solved):controls[j]=F(int(value.p),int(value.q))
    assert all(abs(dot(row,controls)-rhs)<=obs.error for row,rhs,obs in zip(A,b,case.waypoints))
    return controls


def inner_lp(case,mesh):
    n=case.axes*(len(mesh)-1);Aeq=[];beq=[];Aub=[];bub=[];handles=[]
    for obs,b in zip(case.waypoints,adjusted(case)):
        row=[F(0)]+input_row(case,obs,mesh)
        scale=max([abs(b),obs.error,F(1,10**12)]+list(map(abs,row)))
        numeric=[float(v/scale) for v in row]
        if obs.error==0:handles.append(('eq',len(Aeq),scale));Aeq.append(numeric);beq.append(float(b/scale))
        else:
            handles.append(('ineq',len(Aub),scale));Aub.extend([numeric,[-v for v in numeric]])
            bub.extend([float((b+obs.error)/scale),float((-b+obs.error)/scale)])
    for j in range(n):
        for sign in (1,-1):
            row=[0.]*(n+1);row[0]=-1.;row[j+1]=sign;Aub.append(row);bub.append(0.)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result=linprog([1.]+[0.]*n,A_ub=np.array(Aub),b_ub=np.array(bub),
                       A_eq=np.array(Aeq) if Aeq else None,b_eq=np.array(beq) if beq else None,
                       bounds=[(0,None)]+[(None,None)]*n,method='highs',
                       options=dict(threads=1,primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10))
    lam=[]
    if result.status==0:
        for kind,i,scale in handles:
            value=result.eqlin.marginals[i] if kind=='eq' else result.ineqlin.marginals[i]-result.ineqlin.marginals[i+1]
            lam.append(float(value)/float(scale))
    return result,lam,[str(w.message) for w in caught]


def certify(case,proposal,max_iterations=4,tolerance=F(1,100000)):
    started=time.perf_counter();scale=reference_scale(case)
    lam=rational_multipliers(proposal.get('lam',[]),len(case.waypoints))
    bestdual=dual_bound(case,lam);initial_dual=bestdual
    best=None;attempts=[];points=knots(case);mesh=set(points)
    for left,right in zip(points,points[1:]):mesh.update(left+(right-left)*F(j,4) for j in range(1,4))
    for iteration in range(max_iterations):
        for axis,left,right,p in kernel_segments(case,lam):
            for root in real_roots(p,float(left),float(right)):
                value=F(round(root*10**10),10**10)
                if left<value<right:mesh.add(value)
        mesh=sorted(mesh)
        if len(mesh)>257:break
        result,values,notices=inner_lp(case,mesh)
        log=dict(iteration=iteration,cells=len(mesh)-1,status=int(result.status),message=result.message,warnings=notices)
        attempts.append(log)
        if result.status==0:
            lam=rational_multipliers(values,len(case.waypoints));dual=dual_bound(case,lam)
            if F(dual['lower'])>F(bestdual['lower']):bestdual=dual
            try:
                controls=recover_controls(case,mesh,result.x[1:]);peak=max(map(abs,controls))
                if best is None or peak<F(best['upper']):
                    best=dict(case=case.record(),mesh=list(map(str,mesh)),controls=list(map(str,controls)),
                              upper=str(peak),certified=True,tolerance=str(tolerance),scale=str(scale))
                log['primal_recovered']=True
            except (AssertionError,ValueError,ZeroDivisionError,NonInvertibleMatrixError) as exc:
                log['primal_recovered']=False;log['failure']=type(exc).__name__+': '+str(exc)
            if best is not None:
                gap=F(best['upper'])-F(bestdual['lower']);assert gap>=0
                best.update(dual=bestdual,lower=bestdual['lower'],gap=str(gap),converged=gap<=tolerance*scale)
                from .check import check
                check(best);log['gap']=str(gap)
                if best['converged']:break
        mesh=set(mesh)|{(left+right)/2 for left,right in zip(mesh,mesh[1:])}
    if best is None:best=dict(case=case.record(),certified=False,reason='no_exact_primal',dual=bestdual)
    return dict(certificate=best,initial_dual=initial_dual,attempts=attempts,
                repair_lp_calls=len(attempts),elapsed_seconds=time.perf_counter()-started)
