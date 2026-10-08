from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial, comb


def dot(a,b):
    assert len(a)==len(b)
    return sum((x*y for x,y in zip(a,b)),F(0))


@dataclass(frozen=True)
class Observation:
    time: F
    weights: tuple
    value: F
    error: F
    received: F

    def record(self):
        return dict(time=str(self.time),weights=list(map(str,self.weights)),
                    value=str(self.value),error=str(self.error),received=str(self.received))


@dataclass(frozen=True)
class Case:
    label: str
    center: tuple
    radii: tuple
    bound: F
    observations: tuple
    now: F
    target: F
    weights: tuple
    truth_support: object = None

    def validate(self):
        assert len(self.center)==len(self.radii)==len(self.weights)==4
        assert self.bound>0 and min(self.radii)>=0 and 0<=self.now<=self.target and self.target>0
        assert all(len(o.weights)==4 and 0<=o.time<=o.received<=self.now and o.error>=0
                   for o in self.observations), 'Unarrived evidence or invalid observation'

    def record(self):
        self.validate()
        return dict(label=self.label,center=list(map(str,self.center)),radii=list(map(str,self.radii)),
                    bound=str(self.bound),observations=[o.record() for o in self.observations],
                    now=str(self.now),target=str(self.target),weights=list(map(str,self.weights)),
                    truth_support=None if self.truth_support is None else str(self.truth_support))


def from_record(r):
    return Case(r['label'],tuple(map(F,r['center'])),tuple(map(F,r['radii'])),F(r['bound']),
                tuple(Observation(F(o['time']),tuple(map(F,o['weights'])),F(o['value']),
                                  F(o['error']),F(o['received'])) for o in r['observations']),
                F(r['now']),F(r['target']),tuple(map(F,r['weights'])),
                None if r['truth_support'] is None else F(r['truth_support']))


def initial_row(t,w):
    return [sum((w[r]*t**(k-r)/factorial(k-r) for r in range(k+1)),F(0)) for k in range(4)]


def kernel_poly(t,w):
    return [sum((w[r]*F(comb(3-r,j))*t**(3-r-j)*(-1)**j/factorial(3-r)
                 for r in range(4-j)),F(0)) for j in range(4)]


def integrated_kernel(t,w,l,r):
    r=min(t,r)
    if l>=r: return F(0)
    return sum((w[k]*((t-l)**(4-k)-(t-r)**(4-k))/factorial(4-k) for k in range(4)),F(0))


def trajectory_row(t,w,mesh):
    return initial_row(t,w)+[integrated_kernel(t,w,l,r) for l,r in zip(mesh,mesh[1:])]


def state(initial,mesh,inputs,t):
    x=list(initial)
    for l,r,u in zip(mesh,mesh[1:],inputs):
        h=max(F(0),min(t,r)-l)
        x=[sum((x[k]*h**(k-i)/factorial(k-i) for k in range(i,4)),F(0))
           +u*h**(4-i)/factorial(4-i) for i in range(4)]
    return x


def segment_polynomials(case,lam):
    knots=sorted({F(0),case.target}|{o.time for o in case.observations})
    for l,r in zip(knots,knots[1:]):
        p=kernel_poly(case.target,case.weights)
        for o,v in zip(case.observations,lam):
            if r<=o.time:
                p=[a-v*b for a,b in zip(p,kernel_poly(o.time,o.weights))]
        yield l,r,p


def bernstein(p,l,r):
    # Power in s -> power in v=(s-l)/(r-l) -> degree-three Bernstein.
    shifted=[sum((p[j]*comb(j,k)*l**(j-k)*(r-l)**k for j in range(k,4)),F(0)) for k in range(4)]
    return [sum((shifted[k]*F(comb(i,k),comb(3,k)) for k in range(i+1)),F(0)) for i in range(4)]


def dual_certificate(case,lam,tolerance=F(1,10**10)):
    case.validate()
    assert len(lam)==len(case.observations)
    residual=initial_row(case.target,case.weights)
    base=F(0)
    for o,v in zip(case.observations,lam):
        residual=[a-v*b for a,b in zip(residual,initial_row(o.time,o.weights))]
        base+=v*o.value+abs(v)*o.error
    base+=dot(residual,case.center)+dot(list(map(abs,residual)),case.radii)
    leaves=[]
    def recurse(p,l,r,depth=0):
        bs=bernstein(p,l,r)
        up=(r-l)*sum(map(abs,bs),F(0))/4
        lo=abs((r-l)*sum(bs,F(0))/4)
        if up-lo<=tolerance*(r-l)/case.target or depth==30:
            leaves.append([str(l),str(r),*map(str,p)])
            return up
        middle=(l+r)/2
        return recurse(p,l,middle,depth+1)+recurse(p,middle,r,depth+1)
    integral=sum((recurse(p,l,r) for l,r,p in segment_polynomials(case,lam)),F(0))
    return dict(lam=list(map(str,lam)),leaves=leaves,upper=str(base+case.bound*integral))


def cases():
    """Prespecified deterministic pilot, not a random-population performance study."""
    out=[]
    for T in (F(1,2),F(1),F(2)):
        scales=tuple(T**(4-i)/factorial(4-i) for i in range(4))
        for seed in range(4):
            mesh=[T*t for t in (F(0),F(1,7),F(2,5),F(5,8),F(6,7),F(1))]
            inputs=[F(((seed+2)*(j+3))%9-4,5) for j in range(5)]
            initial=(F(0),)*4
            obs=[]
            for k,tfrac in enumerate((F(1,5),F(2,5),F(3,5),F(4,5))):
                t=T*tfrac
                xx=state(initial,mesh,inputs,t)
                for channel in range(3):
                    if channel<2:
                        w=tuple(F(i==channel)/scales[i] for i in range(4))
                    else:
                        w=(F(0),F(0),F(1,2)/scales[2],F((-1)**seed,2)/scales[3])
                    eps=F(1,100) if seed<2 else F(1,1000)
                    offset=eps*F(((k+seed+channel)%5)-2,4)
                    obs.append(Observation(t,w,dot(w,xx)+offset,eps,t+T/20))
            radii=tuple(scale/50 if seed%2 else F(0) for scale in scales)
            for axis in range(5):
                w=tuple(F((-1)**seed if i==axis else 0)/scales[i] for i in range(4)) if axis<4 else tuple(F((-1)**i,4)/scales[i] for i in range(4))
                out.append(Case(f'history_T{T}_s{seed}_axis{axis}',initial,radii,F(1),
                                tuple(obs),9*T/10,T,w))
    for ell in (F(0),F(1,100),F(1,10),F(1,2),F(1)):
        obs=(Observation(F(1),(F(1),F(0),F(0),F(0)),F(1,24),ell**4/12,F(1)),)
        out.append(Case(f'boundary_ell{ell}',(F(0),)*4,(F(0),)*4,F(1),obs,F(1),F(1),
                        (F(0),F(0),F(0),F(-1)),-1+2*ell))
    return out
