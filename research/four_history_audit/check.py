"""Independent certificate replay. Imports no producer, array library or optimizer."""
from fractions import Fraction as F
from math import factorial,comb


def dot(a,b):
    assert len(a)==len(b)
    return sum((x*y for x,y in zip(a,b)),F(0))


def evolve(initial,mesh,controls,t):
    out=[]
    for i in range(4):
        value=sum((initial[j]*t**(j-i)/factorial(j-i) for j in range(i,4)),F(0))
        for l,r,u in zip(mesh,mesh[1:],controls):
            right=min(r,t)
            if l<right: value+=u*((t-l)**(4-i)-(t-right)**(4-i))/factorial(4-i)
        out.append(value)
    return out


def initial_row(t,w):
    return [sum((w[i]*t**(k-i)/factorial(k-i) for i in range(k+1)),F(0)) for k in range(4)]


def polynomial(t,w):
    # Independent direct expansion of (t-s)^j.
    coeff=[F(0)]*4
    for i,weight in enumerate(w):
        degree=3-i
        for j in range(degree+1):
            coeff[j]+=weight*comb(degree,j)*t**(degree-j)*(-1)**j/factorial(degree)
    return coeff


def check_upper(case,dual):
    center,radii,w=(list(map(F,case[k])) for k in ('center','radii','weights'))
    T,now,J=F(case['target']),F(case['now']),F(case['bound'])
    assert len(center)==len(radii)==len(w)==4 and min(radii)>=0
    assert J>0 and 0<=now<=T and T>0
    obs=case['observations']; lam=list(map(F,dual['lam']))
    assert len(lam)==len(obs)
    residual=initial_row(T,w); upper=F(0)
    for o,v in zip(obs,lam):
        t,received,eps=F(o['time']),F(o['received']),F(o['error'])
        ow=list(map(F,o['weights']))
        assert len(ow)==4 and 0<=t<=received<=now and eps>=0
        residual=[a-v*b for a,b in zip(residual,initial_row(t,ow))]
        upper+=v*F(o['value'])+abs(v)*eps
    upper+=dot(center,residual)+dot(radii,list(map(abs,residual)))
    last=F(0); total=F(0)
    for leaf in dual['leaves']:
        assert len(leaf)==6
        left,right,*p=map(F,leaf)
        assert left==last and left<right<=T
        assert not any(left<F(o['time'])<right for o in obs)
        expected=polynomial(T,w)
        for o,v in zip(obs,lam):
            if right<=F(o['time']):
                expected=[a-v*b for a,b in zip(expected,polynomial(F(o['time']),list(map(F,o['weights']))))]
        assert expected==p
        h=right-left
        # Hermite endpoint formula for cubic Bernstein coefficients.
        val=lambda s:sum((p[j]*s**j for j in range(4)),F(0))
        der=lambda s:p[1]+2*p[2]*s+3*p[3]*s*s
        bs=(val(left),val(left)+h*der(left)/3,val(right)-h*der(right)/3,val(right))
        total+=h*sum(map(abs,bs),F(0))/4
        last=right
    assert last==T
    upper+=J*total
    assert upper==F(dual['upper'])
    return upper


def check(record):
    assert record['certified']
    case=record['case']
    upper=check_upper(case,record['dual'])
    center,radii=(list(map(F,case[k])) for k in ('center','radii'))
    mesh=list(map(F,record['mesh'])); x=list(map(F,record['trajectory']))
    T,J=F(case['target']),F(case['bound'])
    assert mesh[0]==0 and mesh[-1]==T and all(l<r for l,r in zip(mesh,mesh[1:]))
    assert len(x)==len(mesh)+3
    assert all(abs(v-c)<=r for v,c,r in zip(x[:4],center,radii))
    assert all(abs(u)<=J for u in x[4:])
    for o in case['observations']:
        value=dot(list(map(F,o['weights'])),evolve(x[:4],mesh,x[4:],F(o['time'])))
        assert abs(value-F(o['value']))<=F(o['error'])
    lower=dot(list(map(F,case['weights'])),evolve(x[:4],mesh,x[4:],T))
    assert lower==F(record['lower'])
    gap=upper-lower
    assert gap>=0 and gap==F(record['gap'])
    assert record['converged']==(gap<=F(record['tolerance']))
    if case['truth_support'] is not None:
        assert lower<=F(case['truth_support'])<=upper
    return gap
