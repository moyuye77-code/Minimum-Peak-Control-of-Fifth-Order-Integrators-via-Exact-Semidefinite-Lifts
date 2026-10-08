from fractions import Fraction as F
from math import factorial
from random import Random
from research.peak_synthesis.model import Problem,Waypoint,state,dot
from research.peak_synthesis.check import evolve

METHODS=('seven_block','generic_polynomial','continuous_root_dual')
SEED=20260928
REPETITIONS=3


def make_case(axes,segments,T=F(1),amplitude=F(1),epsilon=None,label=None):
    if epsilon is None:
        mesh=tuple(T*F(k,3*segments) for k in range(3*segments+1))
        controls=tuple(amplitude*F(((a+2)*(k+3)+k*k)%13-6,6)
                       for a in range(axes) for k in range(3*segments))
    else:
        assert axes==1 and 0<epsilon<1
        mesh=(F(0),T*(1-epsilon),T)
        controls=(amplitude,-amplitude)
    initial=(F(0),)*(4*axes)
    obs=[]
    for n in range(1,segments+1):
        time=T*F(n,segments); x=state(initial,mesh,controls,time)
        assert x==evolve(initial,mesh,controls,time)
        for derivative in range(4 if n==segments else 2):
            for axis in range(axes):
                w=[F(0)]*(4*axes);w[4*axis+derivative]=F(factorial(4-derivative))
                if axes>1 and n<segments:
                    w[4*((axis+1)%axes)+derivative]=F((-1)**(axis+n),3)*factorial(4-derivative)
                obs.append(Waypoint(time,tuple(w),dot(w,x)))
    case=Problem(label,T,initial,tuple(obs));case.validate()
    witness={'mesh':list(map(str,mesh)),'controls':list(map(str,controls)),
             'peak':str(max(map(abs,controls)))}
    return case,witness


def suite():
    out=[]
    for axes in (1,2):
        for segments in (2,4,8):
            label=f'scale_a{axes}_n{segments}'
            c,w=make_case(axes,segments,label=label)
            out.append(dict(case=c,witness=w,group='scale',amplitude=F(1),known=None))
    for T in (F(1,4),F(4)):
        c,w=make_case(2,4,T=T,label=f'time_T{T}')
        out.append(dict(case=c,witness=w,group='time',amplitude=F(1),known=None))
    for eps in (F(1,10),F(1,100),F(1,1000),F(1,10000)):
        c,w=make_case(1,4,epsilon=eps,label=f'saturation_e{eps}')
        out.append(dict(case=c,witness=w,group='saturation',amplitude=F(1),known=F(1)))
    for A in (F(1,10000),F(1000)):
        c,w=make_case(2,4,amplitude=A,label=f'amplitude_A{A}')
        out.append(dict(case=c,witness=w,group='amplitude',amplitude=A,known=None))
    assert len(out)==14
    return out


def schedule():
    rng=Random(SEED);ids=list(range(14));rng.shuffle(ids);out=[]
    for case_id in ids:
        offset=rng.randrange(3)
        for repetition in range(REPETITIONS):
            for position in range(3):
                method=METHODS[(offset+repetition+position)%3]
                out.append(dict(case_id=case_id,repetition=repetition,
                                position=position,method=method))
    return out


def record(item):
    return dict(case=item['case'].record(),witness=item['witness'],group=item['group'],
                amplitude=str(item['amplitude']),
                known=None if item['known'] is None else str(item['known']))
