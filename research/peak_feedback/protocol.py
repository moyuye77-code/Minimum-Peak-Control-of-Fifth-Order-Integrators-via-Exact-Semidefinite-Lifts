from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial
from research.peak_synthesis.model import Problem,Waypoint,state

HORIZON=F(2)
PERIOD=F(1,4)
STEPS=20
INPUT_CAP=F(5)
QUANTUM=F(1,10**9)


@dataclass(frozen=True)
class Scenario:
    label:str
    initial:tuple
    disturbance:str

    @property
    def axes(self): return len(self.initial)//4

    def record(self):
        return dict(label=self.label,initial=list(map(str,self.initial)),disturbance=self.disturbance,
                    horizon=str(HORIZON),period=str(PERIOD),steps=STEPS,
                    input_cap=str(INPUT_CAP),actuator_quantum=str(QUANTUM),
                    exact_state_measurement=True)


def scenarios():
    one=(F(1,10),F(0),F(0),F(0))
    two=one+(F(-2,25),F(0),F(0),F(0))
    return (Scenario('single_nominal',one,'none'),
            Scenario('pair_nominal',two,'none'),
            Scenario('pair_pulse',two,'pulse'),
            Scenario('pair_persistent',two,'persistent'))


def disturbance(scenario,step):
    assert 0<=step<STEPS
    if scenario.disturbance=='none': return (F(0),)*scenario.axes
    if scenario.disturbance=='pulse':
        return (F(1,20),F(-1,25)) if 4<=step<8 else (F(0),F(0))
    assert scenario.disturbance=='persistent'
    return tuple((F(1) if (step+2*a)%8<4 else F(-1))*F(1,30+10*a)
                 for a in range(scenario.axes))


def planning_problem(current):
    """Uses only the current exact state and public constant target/horizon."""
    rows=[]
    for a in range(len(current)//4):
        for d in range(4):
            w=[F(0)]*len(current); w[4*a+d]=F(factorial(4-d))/HORIZON**(4-d)
            rows.append(Waypoint(HORIZON,tuple(w),F(0)))
    return Problem('feedback_to_rest',HORIZON,tuple(current),tuple(rows))


def quantize(value):
    # Exact ties-to-even quantization; the cap is an integer multiple of QUANTUM.
    return max(-INPUT_CAP,min(INPUT_CAP,round(value/QUANTUM)*QUANTUM))


def execute(current,certificate,external):
    """Exact plant propagation of the clipped prefix and quantized actuator levels."""
    axes=len(current)//4
    if certificate is None:
        mesh=[F(0),PERIOD]; intended=[F(0)]*axes
    else:
        oldmesh=list(map(F,certificate['mesh'])); controls=list(map(F,certificate['controls']))
        cells=len(oldmesh)-1
        mesh=[t for t in oldmesh if t<PERIOD]+[PERIOD]
        kept=len(mesh)-1
        intended=[controls[a*cells+k] for a in range(axes) for k in range(kept)]
    cells=len(mesh)-1
    commanded=list(map(quantize,intended))
    actual=[u+external[a] for a in range(axes) for u in commanded[a*cells:(a+1)*cells]]
    next_state=state(current,mesh,actual,PERIOD)
    planned=state(current,mesh,intended,PERIOD)
    for a in range(axes):
        for d in range(4):
            bound=(abs(external[a])+QUANTUM/2)*PERIOD**(4-d)/factorial(4-d)
            assert abs(next_state[4*a+d]-planned[4*a+d])<=bound
    return dict(mesh=list(map(str,mesh)),intended=list(map(str,intended)),
                commanded=list(map(str,commanded)),actual=list(map(str,actual)),
                disturbance=list(map(str,external)),next_state=list(map(str,next_state)),
                planned_prefix_state=list(map(str,planned)))
