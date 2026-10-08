from fractions import Fraction as F
from collections import Counter
import sympy as sp
import pytest
from . import cases,formulations as f
from research.four_lift_prior_audit import algebra as g
from research.peak_synthesis.model import homogeneous,knots,local_kernel
from research.peak_synthesis.check import evolve


def test_suite_generation_and_schedule_are_fixed():
    data=cases.suite();assert len(data)==14
    assert Counter(x['group'] for x in data)==dict(scale=6,time=2,saturation=4,amplitude=2)
    schedule=cases.schedule();assert len(schedule)==126 and schedule==cases.schedule()
    assert len({(r['case_id'],r['repetition'],r['method']) for r in schedule})==126
    for case_id in range(14):
        for method in cases.METHODS:
            assert sorted(r['position'] for r in schedule if r['case_id']==case_id and r['method']==method)==[0,1,2]


def test_withheld_witnesses_satisfy_all_exact_rows():
    for item in cases.suite():
        case=item['case'];w=item['witness'];mesh=list(map(F,w['mesh']));u=list(map(F,w['controls']))
        assert max(map(abs,u))<=item['amplitude']
        for o in case.waypoints:
            x=evolve(case.initial,mesh,u,o.time)
            assert sum(a*b for a,b in zip(o.weights,x))==o.value


def test_generic_coefficients_reconstruct_original_affine_model():
    yy=g.moment_data();variables=(g.s,g.t,g.q,g.r)+tuple(yy[e] for e in f.GENERIC_EXPONENTS)
    built=f.generic_homogeneous(sp.S.One,variables[:4],variables[4:])
    for a,b in zip(built,g.generic_matrices(yy)):
        assert (sp.Matrix(a)-b).applyfunc(sp.expand)==sp.zeros(b.rows)


def test_homogeneous_generic_restricts_to_same_seven_blocks():
    J=sp.Symbol('J');y=sp.symbols('m z q r');aux=sp.symbols('u0:12')
    blocks=list(map(sp.Matrix,f.generic_homogeneous(J,y,aux)))
    aa=dict(zip(f.GENERIC_EXPONENTS,aux))
    small=list(map(sp.Matrix,homogeneous(J,y,[aa[e] for e in ((2,0),(3,0),(4,0),(1,1),(2,1),(0,2))])))
    changes=[sp.Matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,1,0,0]])]
    for sign in (1,-1):
        changes.append(sp.Matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,sp.Rational(sign,6),0,0]]))
    select=sp.Matrix([[1,0,0],[0,1,0]])
    derived=[C*blocks[0]*C.T for C in changes]+[select*b*select.T for b in blocks[1:3]]+blocks[3:]
    for a,b in zip(small,derived):assert (a-b).applyfunc(sp.expand)==sp.zeros(a.rows)


@pytest.mark.parametrize('method,per_segment,sizes',[
    ('seven_block',10,[3,3,3,2,2,2,2]),('generic_polynomial',16,[6,3,3,2,2])])
def test_actual_compiled_sizes(method,per_segment,sizes):
    case=cases.suite()[0]['case'];problem,*_=f.build_sdp(case,method)
    data,_,_=problem.get_problem_data('CLARABEL')
    n=(len(knots(case))-1)*case.axes
    assert data['A'].shape[1]==1+per_segment*n
    assert list(data['dims'].psd)==sizes*n and not data['dims'].soc


def test_saturation_optimum_is_one_by_exact_linear_support():
    t=sp.Symbol('t')
    for item in cases.suite():
        if item['group']!='saturation':continue
        a=sp.Rational(item['witness']['mesh'][1]);p=a-t
        integral=sp.integrate(p,(t,0,a))-sp.integrate(p,(t,a,1))
        assert integral>0 and item['known']==1
        # The last two terminal rows weight acceleration by 2 and jerk by 1.
        obs=item['case'].waypoints
        pairing=sp.Rational(obs[-2].value)/2-(1-a)*sp.Rational(obs[-1].value)
        assert pairing==integral


@pytest.mark.parametrize('A',[F(1,10000),F(1),F(1000)])
def test_recovery_threshold_is_not_looser_than_suite_target(A):
    for U in (A/2,A,2*A):
        assert min(F(1),A)*max(F(1),U)<=max(A,U)
