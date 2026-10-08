from copy import deepcopy
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import ast
import pytest
from research.five_synthesis.certify import dual_bound,reference_scale,certify
from research.five_synthesis.check import evolve
from .protocol import (Scenario,scenarios,planning_problem,window,predict,execute,quantize,
                       HORIZON,PERIOD,CAP,QUANTUM,disturbance)
from .check import extract,replay
from . import experiment


def zero_recovery(case,proposal):
    dual=dual_bound(case,[F(0)]*len(case.waypoints))
    certificate=dict(case=case.record(),mesh=['0','3'],controls=['0']*case.axes,
                     upper='0',lower='0',gap='0',dual=dual,certified=True,converged=True,
                     tolerance='1/100000',scale=str(reference_scale(case)))
    return dict(certificate=certificate,initial_dual=dual,repair_lp_calls=0,attempts=[],elapsed_seconds=0.)


def stub_episode(monkeypatch,elapsed=100_000_000,mode='one_period_deadline'):
    def proposal(case):
        return dict(method='closed_five_sdp',status='optimal',optimization_calls=0,warnings=[])
    monkeypatch.setitem(experiment.METHODS,'closed_five_sdp',proposal)
    monkeypatch.setattr(experiment,'certify',zero_recovery)
    values=iter([v for k in range(12) for v in (k*10**10,k*10**10+elapsed)])
    return experiment.episode(Scenario('zero_test',(F(0),)*5,'none'),mode,'closed_five_sdp',clock=lambda:next(values))


def test_terminal_tasks_use_all_five_states_without_future_disturbance():
    assert len(scenarios())==4
    for scenario in scenarios():
        task=planning_problem(scenario.initial)
        assert task.initial==scenario.initial and task.horizon==HORIZON
        for i,row in enumerate(task.waypoints):
            expected=[F(0)]*len(scenario.initial);expected[i]=F(factorial(5-i%5))/HORIZON**(5-i%5)
            assert row.weights==tuple(expected) and row.value==row.error==0
        assert all(len(disturbance(scenario,k))==scenario.axes for k in range(12))


@pytest.mark.parametrize('offset',(F(0),F(1,2),F(5,2),F(3),F(7,2)))
def test_absolute_queue_switches_and_expiry(offset):
    plan=dict(mesh=['0','1/7','1/3','5/2','3'],controls=['1','-2','3','-4','-1','2','-3','4'])
    active=dict(source=2,start='1',certificate=plan)
    mesh,controls=window(active,F(1)+offset,2)
    assert (mesh,controls)==extract(plan,F(1),F(1)+offset,PERIOD,2)
    assert mesh[0]==0 and mesh[-1]==PERIOD
    if offset>=3:assert not any(controls)


def test_quantization_and_exact_disturbed_execution_integral():
    assert quantize(QUANTUM/2)==0 and quantize(3*QUANTUM/2)==2*QUANTUM
    assert quantize(CAP+QUANTUM)==CAP
    current=[F(j-3,11) for j in range(10)]
    mesh=[F(0),F(1,7),F(1,3),PERIOD]
    intended=[F(j-2,13) for j in range(6)];external=[F(1,20),F(-1,25)]
    act=execute(current,mesh,intended,external)
    assert list(map(F,act['next_state']))==evolve(current,mesh,list(map(F,act['actual'])),PERIOD)
    prediction=predict(current,mesh,intended)
    for a in range(2):
        for d in range(5):
            assert F(act['next_state'][5*a+d])-prediction[5*a+d]==external[a]*PERIOD**(5-d)/factorial(5-d)
    # Constant position two yields integral 4*period exactly.
    constant=execute([F(2),F(0),F(0),F(0),F(0)],[F(0),PERIOD],[F(0)],[F(0)])
    assert F(constant['position_squared_integral'])==4*PERIOD


def test_on_time_plan_never_actuates_early(monkeypatch):
    episode=stub_episode(monkeypatch);summary=replay(episode)
    assert episode['steps'][0]['execution_source'] is None
    assert episode['steps'][1]['execution_source']==0
    assert summary['queries']==summary['accepted']==12 and summary['no_plan_slots']==1


def test_late_result_occupies_worker_and_is_not_reused(monkeypatch):
    episode=stub_episode(monkeypatch,elapsed=1_100_000_000)
    summary=replay(episode)
    assert summary['queries']==4 and summary['busy_skips']==8 and summary['accepted']==0
    assert summary['no_plan_slots']==12 and summary['over_period_queries']==4
    assert [row['decision'] for row in episode['steps'][:3]]==['late_discard','busy_skip','busy_skip']


def test_deadline_equality_is_accepted(monkeypatch):
    summary=replay(stub_episode(monkeypatch,elapsed=500_000_000))
    assert summary['accepted']==12 and summary['busy_skips']==summary['over_period_queries']==0


def test_ideal_reference_records_but_does_not_hide_exceedances(monkeypatch):
    episode=stub_episode(monkeypatch,elapsed=1_100_000_000,mode='ideal_instant')
    summary=replay(episode)
    assert summary['queries']==summary['accepted']==summary['over_period_queries']==12
    assert summary['busy_skips']==0 and episode['steps'][0]['execution_source']==0


@pytest.mark.parametrize('kind',('early','accepted','elapsed','busy','prediction','target','disturbance','state','integral','source','final'))
def test_independent_checker_rejects_schedule_and_plant_tampering(monkeypatch,kind):
    episode=deepcopy(stub_episode(monkeypatch));row=episode['steps'][0]
    if kind=='early':row['execution_source']=0
    elif kind=='accepted':row['request']['accepted']=False
    elif kind=='elapsed':row['request']['elapsed_ns']=600_000_000
    elif kind=='busy':row['busy_until_after_ns']=0
    elif kind=='prediction':row['request']['problem']['initial'][0]='1'
    elif kind=='target':row['request']['problem']['waypoints'][0]['value']='1'
    elif kind=='disturbance':row['execution']['disturbance'][0]='1'
    elif kind=='state':row['execution']['next_state'][0]='1'
    elif kind=='integral':row['execution']['position_squared_integral']='1'
    elif kind=='source':row['active_after']=1
    else:episode['final_state'][0]='1'
    with pytest.raises(AssertionError):replay(episode)


def test_checker_has_no_producer_or_optimizer_imports():
    tree=ast.parse(Path(__file__).with_name('check.py').read_text(encoding='utf-8'))
    imports=[node.module for node in ast.walk(tree) if isinstance(node,ast.ImportFrom)]
    assert set(imports)=={'fractions','math','statistics','research.five_synthesis.check'}
    assert not any(isinstance(node,ast.Import) for node in ast.walk(tree))


def test_frozen_schedule_has_every_combination_once_and_alternates():
    entries=experiment.schedule();assert len(entries)==len(set(entries))==16
    for i in range(0,16,2):
        assert entries[i][:2]==entries[i+1][:2]
        assert entries[i][2]!=entries[i+1][2]


def test_actual_common_recovery_on_terminal_zero_case():
    task=planning_problem((F(0),)*5)
    result=certify(task,{'lam':[0.]*5})
    assert result['certificate']['certified'] and result['certificate']['converged']
    assert F(result['certificate']['upper'])==0


def nonzero_episode(monkeypatch,peak=F(1),late_after_first=False):
    # A single constant-input plan reaches rest; later candidate recoveries fail.
    at_start=[peak*(-1)**(5-d)*HORIZON**(5-d)/factorial(5-d) for d in range(5)]
    initial=tuple(sum((at_start[j]*(-PERIOD)**(j-i)/factorial(j-i) for j in range(i,5)),F(0)) for i in range(5))
    count=0
    def recover(case,proposal):
        nonlocal count
        dual=dual_bound(case,[F(0)]*5)
        if count==0:
            assert list(case.initial)==at_start
            plan=dict(case=case.record(),mesh=['0','3'],controls=[str(peak)],upper=str(peak),
                      lower='0',gap=str(peak),dual=dual,certified=True,converged=False,
                      tolerance='1/100000',scale=str(reference_scale(case)))
        else:plan=dict(case=case.record(),certified=False,dual=dual,reason='test_injected_failure')
        count+=1
        return dict(certificate=plan,initial_dual=dual,attempts=[],repair_lp_calls=0,elapsed_seconds=0.)
    monkeypatch.setitem(experiment.METHODS,'closed_five_sdp',lambda case:
                        dict(method='closed_five_sdp',status='optimal',optimization_calls=0,warnings=[]))
    monkeypatch.setattr(experiment,'certify',recover)
    values=iter([v for k in range(12) for v in (k*10**10,k*10**10+
                  (1_100_000_000 if late_after_first and k>0 else 100_000_000))])
    return experiment.episode(Scenario('nonzero_test',initial,'none'),'one_period_deadline',
                              'closed_five_sdp',clock=lambda:next(values))


def test_wide_valid_plan_then_failed_replans_continue_queue_and_expire(monkeypatch):
    episode=nonzero_episode(monkeypatch);summary=replay(episode)
    assert summary['accepted']==summary['accepted_wide']==1
    assert summary['ineligible']==11 and summary['no_plan_slots']==1 and summary['expired_plan_slots']==5
    assert summary['input_peak']==1 and all(F(x)==0 for x in episode['final_state'])
    assert F(episode['steps'][1]['execution']['intended'][0])==1


def test_overcap_feasible_plan_is_not_clipped_and_called_eligible(monkeypatch):
    episode=nonzero_episode(monkeypatch,peak=F(6));summary=replay(episode)
    assert summary['accepted']==0 and summary['ineligible']==12 and summary['input_peak']==0


def test_busy_overrun_preserves_existing_nonzero_queue(monkeypatch):
    episode=nonzero_episode(monkeypatch,late_after_first=True);summary=replay(episode)
    assert summary['accepted']==1 and summary['busy_skips']>0
    assert episode['steps'][2]['decision']=='busy_skip'
    assert episode['steps'][2]['execution_source']==0 and F(episode['steps'][2]['execution']['command'][0])==1
