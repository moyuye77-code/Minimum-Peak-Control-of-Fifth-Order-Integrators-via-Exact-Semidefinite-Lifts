from fractions import Fraction as F
import pytest
from research import peak_recovery_ablation as ab
from research.peak_synthesis.model import cases,Problem,Waypoint,knots
from research.peak_synthesis.check import check


def simple():
    return Problem('standalone',F(1),(F(0),)*4,
                   (Waypoint(F(1),(F(0),F(0),F(0),F(1)),F(1)),))


@pytest.mark.parametrize('iterations',[1,4])
def test_wrapper_only_supplies_zero_seed_and_iteration_budget(monkeypatch,iterations):
    c=simple(); calls=[]
    def fake(case,proposal,**kw):
        calls.append((case,proposal,kw)); return 'sentinel'
    monkeypatch.setattr(ab,'certify',fake)
    assert ab.run(c,iterations)=='sentinel'
    assert calls==[(c,dict(method='zero_seed',lam=[0.],optimization_calls=0),
                    dict(max_iterations=iterations))]


@pytest.mark.parametrize('iterations',[1,4])
def test_actual_zero_seed_certificate_and_initial_mesh(iterations):
    case=simple(); result=ab.run(case,iterations); rec=result['certificate']
    assert check(rec)==(F(1),F(1)) and rec['converged']
    assert result['attempts'][0]['cells']==4*(len(knots(case))-1)
    assert result['repair_lp_calls']==1


def test_all_prespecified_cases_get_no_information_seed():
    for case,_ in cases():
        assert ab.seed(case)['lam']==[0.]*len(case.waypoints)
        assert set(ab.seed(case))=={'method','lam','optimization_calls'}


def test_budget_outside_preregistered_arms_rejected():
    with pytest.raises(AssertionError): ab.run(simple(),2)
