from dataclasses import dataclass
from fractions import Fraction as F
from math import comb,factorial
from research.moment_four_lift.exact import matrices


def dot(a,b):
    assert len(a)==len(b)
    return sum((x*y for x,y in zip(a,b)),F(0))


@dataclass(frozen=True)
class Waypoint:
    time:F
    weights:tuple
    value:F
    error:F=F(0)

    def record(self):
        return dict(time=str(self.time),weights=list(map(str,self.weights)),
                    value=str(self.value),error=str(self.error))


@dataclass(frozen=True)
class Problem:
    label:str
    horizon:F
    initial:tuple
    waypoints:tuple

    @property
    def axes(self): return len(self.initial)//4

    def validate(self):
        assert self.horizon>0 and len(self.initial)>0 and len(self.initial)%4==0
        assert len(self.waypoints)>0
        for o in self.waypoints:
            assert len(o.weights)==len(self.initial) and 0<o.time<=self.horizon
            assert o.error>=0

    def record(self):
        self.validate()
        return dict(label=self.label,horizon=str(self.horizon),initial=list(map(str,self.initial)),
                    waypoints=[o.record() for o in self.waypoints])


def from_record(c):
    return Problem(c['label'],F(c['horizon']),tuple(map(F,c['initial'])),
                   tuple(Waypoint(F(o['time']),tuple(map(F,o['weights'])),F(o['value']),F(o['error']))
                         for o in c['waypoints']))


def homogeneous(J,y,aux):
    """Only the three constant Gram top-left entries become J. No division."""
    out=matrices(y,aux)
    for block in out[:3]: block[0][0]=J
    return out


def ballistic(initial,time):
    return [sum((initial[4*a+j]*time**(j-i)/factorial(j-i) for j in range(i,4)),F(0))
            for a in range(len(initial)//4) for i in range(4)]


def adjusted(problem):
    return [o.value-dot(o.weights,ballistic(problem.initial,o.time)) for o in problem.waypoints]


def knots(problem):
    return sorted({F(0),problem.horizon}|{o.time for o in problem.waypoints})


def kernel(o,axis):
    w=o.weights[4*axis:4*axis+4]
    return [sum((w[i]*comb(3-i,j)*o.time**(3-i-j)*(-1)**j/factorial(3-i)
                 for i in range(4-j)),F(0)) for j in range(4)]


def kernel_segments(problem,lam):
    points=knots(problem)
    for axis in range(problem.axes):
        for l,r in zip(points,points[1:]):
            p=[F(0)]*4
            for o,v in zip(problem.waypoints,lam):
                if r<=o.time: p=[x+v*y for x,y in zip(p,kernel(o,axis))]
            yield axis,l,r,p


def local_kernel(o,axis,l,r):
    if r>o.time: return [F(0)]*4
    p=kernel(o,axis); h=r-l
    return [h**(j+1)*sum((p[k]*comb(k,j)*l**(k-j) for k in range(j,4)),F(0))
            for j in range(4)]


def input_row(problem,o,mesh):
    row=[]
    for axis in range(problem.axes):
        for l,r in zip(mesh,mesh[1:]):
            rr=min(r,o.time)
            row.append(F(0) if rr<=l else sum((
                o.weights[4*axis+i]*((o.time-l)**(4-i)-(o.time-rr)**(4-i))/factorial(4-i)
                for i in range(4)),F(0)))
    return row


def state(initial,mesh,controls,time):
    # Sequential dynamics; independent checker uses convolution instead.
    axes=len(initial)//4; cells=len(mesh)-1; x=list(initial)
    assert len(controls)==axes*cells
    for k,(l,r) in enumerate(zip(mesh,mesh[1:])):
        h=max(F(0),min(time,r)-l)
        x=[sum((x[4*a+j]*h**(j-i)/factorial(j-i) for j in range(i,4)),F(0))
           +controls[a*cells+k]*h**(4-i)/factorial(4-i)
           for a in range(axes) for i in range(4)]
    return x


def cases():
    """18 deterministic pilot problems. Time-scaled copies are not iid trials."""
    out=[]
    for T in (F(1,2),F(1),F(2)):
        for J in (F(1),F(3,2)):
            terminal=(J*T**4/384,F(0),F(0),F(0))
            obs=[]
            for i in range(4):
                w=tuple(F(j==i)*factorial(4-i)/T**(4-i) for j in range(4))
                obs.append(Waypoint(T,w,dot(w,terminal)))
            out.append((Problem(f'rest_T{T}_J{J}',T,(F(0),)*4,tuple(obs)),J))
    for axes in (1,2):
        for seed in (1,2):
            for error in (F(0),F(1,200)):
                T=F(1); initial=tuple(F((i+seed)%5-2,10) for i in range(4*axes))
                mesh=list(map(F,(0,F(1,7),F(2,5),F(5,8),F(6,7),1)))
                controls=[F(((seed+a+2)*(j+2))%11-5,4) for a in range(axes) for j in range(5)]
                obs=[]
                for t0 in (F(1,4),F(1,2),F(3,4),F(1)):
                    x=state(initial,mesh,controls,t0)
                    # Terminal derivatives and intermediate position/velocity.
                    for i in range(4 if t0==1 else 2):
                        for a in range(axes):
                            w=[F(0)]*(4*axes); w[4*a+i]=F(factorial(4-i))
                            if axes==2 and t0<1:
                                w[4*(1-a)+i]=F((-1)**(a+seed),2)*factorial(4-i)
                            obs.append(Waypoint(t0,tuple(w),dot(w,x),error if t0<1 else F(0)))
                out.append((Problem(f'waypoints_a{axes}_s{seed}_e{error}',T,initial,tuple(obs)),None))
    for axes in (1,2):
        initial=tuple(F(i-2,5) for i in range(4*axes)); obs=[]
        for time in (F(1,2),F(1)):
            x=ballistic(initial,time)
            for i in range(4*axes):
                w=tuple(F(j==i) for j in range(4*axes))
                obs.append(Waypoint(time,w,x[i]))
        out.append((Problem(f'zero_peak_a{axes}',F(1),initial,tuple(obs)),F(0)))
    for J in (F(1,10000),F(1000)):
        # Exact peak lower bound follows from a required total input integral.
        o=Waypoint(F(1),(F(0),F(0),F(0),F(1)),J)
        out.append((Problem(f'total_input_J{J}',F(1),(F(0),)*4,(o,)),J))
    assert len(out)==18
    return out
