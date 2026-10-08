from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial
from research.five_synthesis.model import Problem,Waypoint,state

HORIZON=F(3)
PERIOD=F(1,2)
STEPS=12
CAP=F(5)
QUANTUM=F(1,10**9)
PERIOD_NS=500_000_000
MODES=('ideal_instant','one_period_deadline')


@dataclass(frozen=True)
class Scenario:
    label:str
    initial:tuple
    disturbance:str

    @property
    def axes(self):return len(self.initial)//5

    def record(self):
        return dict(label=self.label,initial=list(map(str,self.initial)),disturbance=self.disturbance,
                    horizon=str(HORIZON),period=str(PERIOD),steps=STEPS,cap=str(CAP),
                    quantum=str(QUANTUM),exact_state=True,order=5)


def scenarios():
    one=(F(1,10),F(0),F(0),F(0),F(0))
    two=one+(F(-2,25),F(0),F(0),F(0),F(0))
    return (Scenario('single_nominal',one,'none'),Scenario('pair_nominal',two,'none'),
            Scenario('pair_pulse',two,'pulse'),Scenario('pair_persistent',two,'persistent'))


def disturbance(scenario,step):
    assert 0<=step<STEPS
    if scenario.disturbance=='none':return (F(0),)*scenario.axes
    if scenario.disturbance=='pulse':
        return (F(1,20),F(-1,25)) if 3<=step<6 else (F(0),F(0))
    assert scenario.disturbance=='persistent'
    return tuple((1 if (step+a)%4<2 else -1)*F(1,30+10*a) for a in range(scenario.axes))


def planning_problem(initial):
    """Only the measured/predicted state and public rest target enter the planner."""
    rows=[]
    for i in range(len(initial)):
        weights=[F(0)]*len(initial)
        weights[i]=F(factorial(5-i%5))/HORIZON**(5-i%5)
        rows.append(Waypoint(HORIZON,tuple(weights),F(0)))
    return Problem('feedback_five_to_rest',HORIZON,tuple(initial),tuple(rows))


def quantize(value):return max(-CAP,min(CAP,round(value/QUANTUM)*QUANTUM))


def window(active,start,axes):
    """Absolute-time remaining plan, followed by zero when its horizon expires."""
    end=start+PERIOD
    if active is None:return [F(0),PERIOD],[F(0)]*axes
    origin=F(active['start']);plan=active['certificate']
    full=[origin+F(t) for t in plan['mesh']];u=list(map(F,plan['controls']))
    n=len(full)-1;cuts=sorted({start,end}|{t for t in full if start<t<end})
    assert start>=origin and len(u)==axes*n
    intended=[]
    for a in range(axes):
        for left,right in zip(cuts,cuts[1:]):
            index=next((k for k,(l,r) in enumerate(zip(full,full[1:])) if l<=left and right<=r),None)
            intended.append(F(0) if index is None else u[a*n+index])
    return [t-start for t in cuts],intended


def predict(current,mesh,intended):
    return state(current,mesh,list(map(quantize,intended)),PERIOD)


def execute(current,mesh,intended,external):
    axes=len(current)//5;n=len(mesh)-1
    command=list(map(quantize,intended))
    actual=[u+external[a] for a in range(axes) for u in command[a*n:(a+1)*n]]
    next_state=state(current,mesh,actual,PERIOD)
    nominal=state(current,mesh,intended,PERIOD)
    integral=F(0);point=list(current)
    for k,(left,right) in enumerate(zip(mesh,mesh[1:])):
        h=right-left;controls=[actual[a*n+k] for a in range(axes)]
        for a in range(axes):
            coeff=[point[5*a+j]/factorial(j) for j in range(5)]+[controls[a]/factorial(5)]
            integral+=sum((x*y*h**(i+j+1)/(i+j+1) for i,x in enumerate(coeff)
                           for j,y in enumerate(coeff)),F(0))
        point=state(point,[F(0),h],controls,h)
    assert point==next_state
    for a in range(axes):
        for d in range(5):
            bound=(abs(external[a])+QUANTUM/2)*PERIOD**(5-d)/factorial(5-d)
            assert abs(next_state[5*a+d]-nominal[5*a+d])<=bound
    return dict(mesh=list(map(str,mesh)),intended=list(map(str,intended)),command=list(map(str,command)),
                actual=list(map(str,actual)),disturbance=list(map(str,external)),
                next_state=list(map(str,next_state)),nominal_state=list(map(str,nominal)),
                position_squared_integral=str(integral))
