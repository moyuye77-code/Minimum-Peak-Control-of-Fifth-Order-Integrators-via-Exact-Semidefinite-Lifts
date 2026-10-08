"""Independent rational checker: no producer, array library, or optimizer imports."""
from fractions import Fraction as F
from math import factorial,comb


def validated(case):
    assert case['order']==5
    T=F(case['horizon']);initial=list(map(F,case['initial']));rows=case['waypoints']
    assert T>0 and initial and len(initial)%5==0 and rows
    for row in rows:
        assert 0<F(row['time'])<=T and F(row['error'])>=0 and len(row['weights'])==len(initial)
    return T,initial,rows


def evolve(initial,mesh,controls,time):
    axes=len(initial)//5;n=len(mesh)-1;out=[]
    assert len(controls)==axes*n
    for a in range(axes):
        for i in range(5):
            value=sum((initial[5*a+j]*time**(j-i)/factorial(j-i) for j in range(i,5)),F(0))
            for k,(left,right) in enumerate(zip(mesh,mesh[1:])):
                rr=min(right,time)
                if left<rr:value+=controls[a*n+k]*((time-left)**(5-i)-(time-rr)**(5-i))/factorial(5-i)
            out.append(value)
    return out


def endpoint_bernstein(p,left,right):
    value=lambda t:sum((p[j]*t**j for j in range(5)),F(0))
    deriv=lambda t:sum((j*p[j]*t**(j-1) for j in range(1,5)),F(0))
    second=lambda t:sum((j*(j-1)*p[j]*t**(j-2) for j in range(2,5)),F(0))
    h=right-left
    return [value(left),value(left)+h*deriv(left)/4,
            value(left)+h*deriv(left)/2+h*h*second(left)/12,
            value(right)-h*deriv(right)/4,value(right)]


def row_poly(obs,axis):
    time=F(obs['time']);out=[F(0)]*5
    for i in range(5):
        weight=F(obs['weights'][5*axis+i]);degree=4-i
        for j in range(degree+1):out[j]+=weight*comb(degree,j)*time**(degree-j)*(-1)**j/factorial(degree)
    return out


def rhs_for(obs,initial,T):
    zero=evolve(initial,[F(0),T],[F(0)]*(len(initial)//5),F(obs['time']))
    return F(obs['value'])-sum((F(w)*x for w,x in zip(obs['weights'],zero)),F(0))


def reference_scale(case):
    T,initial,rows=validated(case);axes=len(initial)//5;candidates=[]
    for row in rows:
        time=F(row['time'])
        norm=sum((time*sum(map(abs,endpoint_bernstein(row_poly(row,a),F(0),time)),F(0))/5
                  for a in range(axes)),F(0))
        numerator=max(F(0),abs(rhs_for(row,initial,T))-F(row['error']))
        if norm:candidates.append(numerator/norm)
        else:assert numerator==0
    return max(candidates,default=F(0)) or F(1)


def check_lower(case,dual):
    T,initial,rows=validated(case);axes=len(initial)//5
    lam=list(map(F,dual['lam']));assert len(lam)==len(rows)
    numerator=sum((v*rhs_for(row,initial,T)-abs(v)*F(row['error']) for row,v in zip(rows,lam)),F(0))
    last=[F(0)]*axes;total=F(0)
    for leaf in dual['leaves']:
        assert len(leaf)==8
        axis=leaf[0];assert type(axis) is int and 0<=axis<axes
        left,right,*p=map(F,leaf[1:])
        assert left==last[axis] and left<right<=T
        assert not any(left<F(row['time'])<right for row in rows)
        wanted=[F(0)]*5
        for row,value in zip(rows,lam):
            if right<=F(row['time']):wanted=[a+value*b for a,b in zip(wanted,row_poly(row,axis))]
        assert p==wanted
        bs=endpoint_bernstein(p,left,right)
        total+=(right-left)*sum(map(abs,bs),F(0))/5;last[axis]=right
    assert last==[T]*axes and total==F(dual['norm_upper']) and numerator==F(dual['numerator'])
    if total==0:assert numerator<=0;lower=F(0)
    else:lower=max(F(0),numerator/total)
    assert lower==F(dual['lower'])
    return lower


def check(record):
    assert record['certified']
    case=record['case'];T,initial,rows=validated(case)
    mesh=list(map(F,record['mesh']));controls=list(map(F,record['controls']))
    assert mesh[0]==0 and mesh[-1]==T and all(a<b for a,b in zip(mesh,mesh[1:]))
    assert len(controls)==len(initial)//5*(len(mesh)-1)
    upper=max(map(abs,controls));assert upper==F(record['upper'])
    for row in rows:
        point=evolve(initial,mesh,controls,F(row['time']))
        value=sum((F(w)*x for w,x in zip(row['weights'],point)),F(0))
        assert abs(value-F(row['value']))<=F(row['error']),'Waypoint mismatch'
    lower=check_lower(case,record['dual']);gap=upper-lower
    assert gap>=0 and lower==F(record['lower']) and gap==F(record['gap'])
    tolerance=F(record['tolerance']);scale=reference_scale(case)
    assert tolerance>0 and scale==F(record['scale'])
    assert record['converged']==(gap<=tolerance*scale)
    return lower,upper
