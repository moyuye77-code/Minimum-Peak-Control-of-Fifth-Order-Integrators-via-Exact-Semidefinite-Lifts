"""Standalone exact checker. No producer imports, NumPy, SymPy, or optimization."""
from fractions import Fraction as F
from math import factorial,comb


def validated(case):
    T=F(case['horizon']); initial=list(map(F,case['initial'])); obs=case['waypoints']
    assert T>0 and initial and len(initial)%4==0 and obs
    for o in obs:
        assert 0<F(o['time'])<=T and F(o['error'])>=0
        assert len(o['weights'])==len(initial)
    return T,initial,obs


def evolve(initial,mesh,controls,time):
    axes=len(initial)//4; n=len(mesh)-1; out=[]
    assert len(controls)==axes*n
    for a in range(axes):
        for i in range(4):
            v=sum((initial[4*a+j]*time**(j-i)/factorial(j-i) for j in range(i,4)),F(0))
            for k,(l,r) in enumerate(zip(mesh,mesh[1:])):
                rr=min(r,time)
                if l<rr: v+=controls[a*n+k]*((time-l)**(4-i)-(time-rr)**(4-i))/factorial(4-i)
            out.append(v)
    return out


def check_lower(case,dual):
    T,initial,obs=validated(case); axes=len(initial)//4
    lam=list(map(F,dual['lam'])); assert len(lam)==len(obs)
    numerator=F(0)
    for o,v in zip(obs,lam):
        weights=list(map(F,o['weights'])); time=F(o['time'])
        zero=evolve(initial,[F(0),T],[F(0)]*axes,time)
        b=F(o['value'])-sum((w*x for w,x in zip(weights,zero)),F(0))
        numerator+=v*b-abs(v)*F(o['error'])
    last=[F(0)]*axes; total=F(0)
    for leaf in dual['leaves']:
        assert len(leaf)==7
        axis=leaf[0]; assert isinstance(axis,int) and 0<=axis<axes
        left,right,*p=map(F,leaf[1:])
        assert left==last[axis] and left<right<=T
        assert not any(left<F(o['time'])<right for o in obs)
        wanted=[F(0)]*4
        for o,v in zip(obs,lam):
            time=F(o['time'])
            if right<=time:
                for i in range(4):
                    w=F(o['weights'][4*axis+i]); deg=3-i
                    for j in range(deg+1):
                        wanted[j]+=v*w*comb(deg,j)*time**(deg-j)*(-1)**j/factorial(deg)
        assert p==wanted
        value=lambda x:sum((p[j]*x**j for j in range(4)),F(0))
        derivative=lambda x:p[1]+2*p[2]*x+3*p[3]*x*x
        h=right-left
        bs=[value(left),value(left)+h*derivative(left)/3,
            value(right)-h*derivative(right)/3,value(right)]
        total+=h*sum(map(abs,bs),F(0))/4
        last[axis]=right
    assert last==[T]*axes
    assert total==F(dual['norm_upper']) and numerator==F(dual['numerator'])
    if total==0:
        assert numerator<=0; lower=F(0)
    else: lower=max(F(0),numerator/total)
    assert lower==F(dual['lower'])
    return lower


def check(record):
    assert record['certified']
    case=record['case']; T,initial,obs=validated(case)
    mesh=list(map(F,record['mesh'])); controls=list(map(F,record['controls']))
    assert mesh[0]==0 and mesh[-1]==T and all(l<r for l,r in zip(mesh,mesh[1:]))
    assert len(controls)==len(initial)//4*(len(mesh)-1)
    upper=max(map(abs,controls)); assert upper==F(record['upper'])
    for o in obs:
        x=evolve(initial,mesh,controls,F(o['time']))
        value=sum((F(w)*v for w,v in zip(o['weights'],x)),F(0))
        assert abs(value-F(o['value']))<=F(o['error']), 'Waypoint mismatch'
    lower=check_lower(case,record['dual']); gap=upper-lower
    assert gap>=0 and lower==F(record['lower']) and gap==F(record['gap'])
    tolerance=F(record['tolerance']); assert tolerance>0
    assert record['converged']==(gap<=tolerance*max(F(1),upper))
    return lower,upper
