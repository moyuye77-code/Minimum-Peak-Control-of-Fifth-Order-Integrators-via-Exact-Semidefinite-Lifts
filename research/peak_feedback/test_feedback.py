from copy import deepcopy
from fractions import Fraction as F
from math import factorial
import pytest
from research.peak_synthesis.check import evolve
from .protocol import (scenarios,planning_problem,quantize,execute,disturbance,HORIZON,
                       PERIOD,QUANTUM,INPUT_CAP)


def fake_plan():
    return dict(mesh=['0','1/10','1/5','2'],
                controls=['1/3','-1/7','1/11','-1/4','2/7','-3/11'])


def test_planner_only_uses_present_state_and_fixed_public_parameters():
    initial=scenarios()[1].initial; case=planning_problem(initial)
    assert case.initial==initial and case.horizon==HORIZON and len(case.waypoints)==8
    assert all(o.time==HORIZON and o.value==0 and o.error==0 for o in case.waypoints)


@pytest.mark.parametrize('numerator',[-10000000001,-10000000000,-1,0,1,10000000000,10000000001])
def test_exact_quantization_and_cap(numerator):
    x=F(numerator,2*10**9); y=quantize(x)
    assert abs(y)<=INPUT_CAP and (y/QUANTUM).denominator==1
    if abs(x)<=INPUT_CAP: assert abs(x-y)<=QUANTUM/2


def test_prefix_cuts_partial_cell_and_plant_matches_independent_convolution():
    x=scenarios()[1].initial; plan=fake_plan(); noise=(F(1,20),F(-1,25))
    act=execute(x,plan,noise)
    assert list(map(F,act['mesh']))==[F(0),F(1,10),F(1,5),PERIOD]
    assert list(map(F,act['next_state']))==evolve(x,list(map(F,act['mesh'])),list(map(F,act['actual'])),PERIOD)
    assert plan==fake_plan()


def test_zero_fallback_is_explicit_and_not_a_feasible_plan_claim():
    x=scenarios()[1].initial; act=execute(x,None,(F(0),F(0)))
    assert act['commanded']==['0','0'] and act['mesh']==['0',str(PERIOD)]
    assert list(map(F,act['next_state']))==list(x)


def test_exact_local_execution_error_bound_with_disturbance():
    x=scenarios()[1].initial; w=(F(1,20),F(-1,25)); act=execute(x,fake_plan(),w)
    actual=list(map(F,act['next_state'])); planned=list(map(F,act['planned_prefix_state']))
    for a in range(2):
        for d in range(4):
            assert abs(actual[4*a+d]-planned[4*a+d]) <= (abs(w[a])+QUANTUM/2)*PERIOD**(4-d)/factorial(4-d)


def test_prespecified_scenarios_disturbances_are_not_random_trials():
    ss=scenarios(); assert len(ss)==4
    assert all(disturbance(ss[0],k)==(F(0),) for k in range(20))
    assert disturbance(ss[2],3)==(F(0),F(0))
    assert disturbance(ss[2],4)==(F(1,20),F(-1,25))
    assert disturbance(ss[2],8)==(F(0),F(0))
    assert all(abs(w)<=F(1,20) for s in ss for k in range(20) for w in disturbance(s,k))
