from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
import sympy as s
from research.closed_five_lift.algebra import model as closed_model
from research.peak_synthesis.model import homogeneous as prefix_homogeneous

ORDER=5


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
        return dict(time=str(self.time),weights=list(map(str,self.weights)),value=str(self.value),error=str(self.error))


@dataclass(frozen=True)
class Problem:
    label:str
    horizon:F
    initial:tuple
    waypoints:tuple
    @property
    def axes(self):return len(self.initial)//ORDER
    def validate(self):
        assert self.horizon>0 and self.initial and len(self.initial)%ORDER==0 and self.waypoints
        for row in self.waypoints:
            assert len(row.weights)==len(self.initial) and 0<row.time<=self.horizon and row.error>=0
    def record(self):
        self.validate()
        return dict(order=ORDER,label=self.label,horizon=str(self.horizon),initial=list(map(str,self.initial)),
                    waypoints=[row.record() for row in self.waypoints])


def from_record(record):
    assert record['order']==ORDER
    return Problem(record['label'],F(record['horizon']),tuple(map(F,record['initial'])),
                   tuple(Waypoint(F(row['time']),tuple(map(F,row['weights'])),F(row['value']),F(row['error']))
                         for row in record['waypoints']))


@lru_cache(maxsize=1)
def exact_templates():
    data=closed_model();theta=data['theta']
    return tuple(tuple(tuple(tuple((k,F(expression.diff(var))) for k,var in enumerate(theta)
                                   if expression.diff(var)!=0) for expression in row)
                       for row in mat.tolist()) for mat in data['matrices'].values())


def homogeneous(J,Y,aux):
    assert len(Y)==5 and len(aux)==34
    out=prefix_homogeneous(J,Y[:4],aux[:6])
    for values in ([J,*Y[:4],*aux[6:20]],
                   [J,*(J*F(1,i+1)-Y[i] for i in range(4)),*aux[20:34]]):
        for matrix in exact_templates():
            out.append([[sum((coef*values[k] for k,coef in terms),0) for terms in row] for row in matrix])
    return out,[Y[4]-aux[6],J*F(1,5)-Y[4]-aux[20]]


def ballistic(initial,time):
    return [sum((initial[5*a+j]*time**(j-i)/factorial(j-i) for j in range(i,5)),F(0))
            for a in range(len(initial)//5) for i in range(5)]


def adjusted(case):return [row.value-dot(row.weights,ballistic(case.initial,row.time)) for row in case.waypoints]
def knots(case):return sorted({F(0),case.horizon}|{row.time for row in case.waypoints})


def kernel(row,axis):
    weights=row.weights[5*axis:5*axis+5]
    return [sum((weights[i]*comb(4-i,j)*row.time**(4-i-j)*(-1)**j/factorial(4-i)
                 for i in range(5-j)),F(0)) for j in range(5)]


def kernel_segments(case,lam):
    points=knots(case)
    for axis in range(case.axes):
        for left,right in zip(points,points[1:]):
            p=[F(0)]*5
            for row,value in zip(case.waypoints,lam):
                if right<=row.time:p=[a+value*b for a,b in zip(p,kernel(row,axis))]
            yield axis,left,right,p


def local_kernel(row,axis,left,right):
    if right>row.time:return [F(0)]*5
    p=kernel(row,axis);h=right-left
    return [h**(j+1)*sum((p[k]*comb(k,j)*left**(k-j) for k in range(j,5)),F(0)) for j in range(5)]


def input_row(case,row,mesh):
    out=[]
    for axis in range(case.axes):
        for left,right in zip(mesh,mesh[1:]):
            rr=min(right,row.time)
            out.append(F(0) if rr<=left else sum((row.weights[5*axis+i]*
                       ((row.time-left)**(5-i)-(row.time-rr)**(5-i))/factorial(5-i)
                       for i in range(5)),F(0)))
    return out


def state(initial,mesh,controls,time):
    axes=len(initial)//5;cells=len(mesh)-1;x=list(initial)
    assert len(controls)==axes*cells
    for k,(left,right) in enumerate(zip(mesh,mesh[1:])):
        h=max(F(0),min(time,right)-left)
        x=[sum((x[5*a+j]*h**(j-i)/factorial(j-i) for j in range(i,5)),F(0))
           +controls[a*cells+k]*h**(5-i)/factorial(5-i) for a in range(axes) for i in range(5)]
    return x


def cases():
    """15 fixed problems; generating controls never leave this function."""
    out=[]
    for axes in (1,2):
        initial=tuple(F(i-2,5) for i in range(5*axes));rows=[]
        for t in (F(1,2),F(1)):
            target=ballistic(initial,t)
            for i in range(5*axes):
                rows.append(Waypoint(t,tuple(F(j==i) for j in range(5*axes)),target[i]))
        out.append((Problem(f'zero_a{axes}',F(1),initial,tuple(rows)),F(0)))
    for peak in (F(1,10000),F(1000)):
        row=Waypoint(F(1),(F(0),)*4+(F(1),),peak)
        out.append((Problem(f'total_J{peak}',F(1),(F(0),)*5,(row,)),peak))
    anchors=((F(1),F(1),(F(1,5),F(2,5),F(3,5),F(4,5))),
             (F(1,2),F(3,2),(F(1,5),F(2,5),F(3,5),F(4,5))),
             (F(1),F(1),(F(1,100),F(1,3),F(2,3),F(99,100))))
    for number,(T,peak,roots) in enumerate(anchors):
        mesh=(F(0),*(T*x for x in roots),T)
        controls=tuple(peak*((-1)**i) for i in range(5))
        target=state((F(0),)*5,mesh,controls,T);rows=[]
        for i in range(5):
            weights=tuple(F(j==i)*factorial(5-i)/T**(5-i) for j in range(5))
            rows.append(Waypoint(T,weights,dot(weights,target)))
        out.append((Problem(f'anchor_{number}',T,(F(0),)*5,tuple(rows)),peak))
    for axes in (1,2):
        for seed in (1,2):
            for error in (F(0),F(1,200)):
                initial=tuple(F((i+seed)%5-2,10) for i in range(5*axes))
                mesh=(F(0),F(1,7),F(2,5),F(5,8),F(6,7),F(1))
                controls=[F(((seed+a+2)*(j+2))%11-5,4) for a in range(axes) for j in range(5)]
                rows=[]
                for t in (F(1,4),F(1,2),F(3,4),F(1)):
                    target=state(initial,mesh,controls,t)
                    for i in range(5 if t==1 else 2):
                        for a in range(axes):
                            weights=[F(0)]*(5*axes);weights[5*a+i]=F(factorial(5-i))
                            if axes==2 and t<1:weights[5*(1-a)+i]=F((-1)**(a+seed),2)*factorial(5-i)
                            rows.append(Waypoint(t,tuple(weights),dot(weights,target),error if t<1 else F(0)))
                out.append((Problem(f'waypoints_a{axes}_s{seed}_e{error}',F(1),initial,tuple(rows)),None))
    assert len(out)==15
    return out
