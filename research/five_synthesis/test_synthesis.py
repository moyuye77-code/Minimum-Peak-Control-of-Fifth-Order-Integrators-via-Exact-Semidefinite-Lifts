from copy import deepcopy
from fractions import Fraction as F
from math import factorial
import ast
from pathlib import Path
import pytest
import sympy as s
from .model import (Problem,Waypoint,cases,homogeneous,state,ballistic,adjusted,
                    input_row,dot,local_kernel,kernel,knots,from_record)
from .certify import bernstein,dual_bound,reference_scale,certify
from .check import check,check_lower,evolve,endpoint_bernstein,reference_scale as checked_scale
from .proposals import build_sdp
from research.closed_five_lift.algebra import witness,complement,model as closed_model
from research.moment_four_lift.exact import atomic_aux
from research.moment_upper_gate.exact import psd_pivots


def simple_case():
    return Problem('simple',F(1),(F(0),)*5,(Waypoint(F(1),(F(0),)*4+(F(1),),F(1)),))


def simple_certificate():
    case=simple_case();dual=dual_bound(case,[F(1)])
    return dict(case=case.record(),mesh=['0','1'],controls=['1'],upper='1',lower='1',gap='0',
                dual=dual,tolerance='1/100000',scale='1',converged=True,certified=True)


def test_homogeneous_pencils_and_slacks_are_exactly_linear():
    J=s.Symbol('J');Y=s.symbols('Y0:5');aux=s.symbols('a0:34')
    left,slacks=homogeneous(J,Y,aux)
    right,ss=homogeneous(s.S.One,[v/J for v in Y],[v/J for v in aux])
    for a,b in zip(left,right):
        assert (s.Matrix(a)-J*s.Matrix(b)).applyfunc(s.cancel)==s.zeros(len(a))
        assert all(s.Poly(v,J,*Y,*aux).total_degree()<=1 for row in a for v in row)
    assert all(s.cancel(a-J*b)==0 for a,b in zip(slacks,ss))


def test_last_coordinate_apex_argument_has_a_nonnegative_diagonal_identity():
    data=closed_model();gs=data['matrices']
    nonnegative=gs['unit_A'][0,0]+gs['delta_B'][0,0]+gs['unit_h'][0,0]/12
    nonnegative+=gs['delta_D'][0,0]/3+gs['delta_h'][0,0]/180
    assert s.expand(nonnegative-data['outputs'][5])==0
    zero,slacks=homogeneous(F(0),[F(0)]*5,[F(0)]*34)
    assert all(v==0 for g in zero for row in g for v in row) and slacks==[0,0]


@pytest.mark.parametrize('peak',(F(0),F(1,10000),F(2),F(1000)))
@pytest.mark.parametrize('point',((s.Rational(1,2),s.Rational(1,4),s.Rational(1,6),s.Rational(1,8),s.Rational(1,10)),
                                 (s.S.One,s.Rational(1,2),s.Rational(1,3),s.Rational(1,4),s.Rational(1,5))))
def test_finite_scaled_witness_including_apex(peak,point):
    lo,_=witness(point[:4]);hi,_=witness(complement(point)[:4])
    aux=list(atomic_aux(*point[:2]))+list(map(F,lo[5:]))+list(map(F,hi[5:]))
    mats,slacks=homogeneous(peak,[peak*F(v) for v in point],[peak*v for v in aux])
    for mat in mats:psd_pivots(mat)
    assert min(slacks)>=0


@pytest.mark.parametrize('coefficients',((F(1),F(-2),F(3),F(-4),F(5)),(F(0),F(0),F(0),F(0),F(1)),(F(-1),)*5))
def test_quartic_bernstein_independent_formula(coefficients):
    left,right=F(1,7),F(4,5)
    bs=bernstein(coefficients,left,right)
    assert bs==endpoint_bernstein(coefficients,left,right)
    t=s.Symbol('t')
    expansion=sum(s.Rational(v)*s.binomial(4,j)*t**j*(1-t)**(4-j) for j,v in enumerate(bs))
    wanted=sum(s.Rational(v)*(s.Rational(left)+s.Rational(right-left)*t)**j for j,v in enumerate(coefficients))
    assert s.expand(expansion-wanted)==0


def test_local_kernel_and_two_independent_state_integrations():
    initial=tuple(F(j-4,13) for j in range(10));mesh=(F(0),F(1,4),F(2,3),F(1))
    controls=tuple(F(j-3,7) for j in range(6));time=F(4,5)
    assert state(initial,mesh,controls,time)==evolve(initial,mesh,controls,time)
    row=Waypoint(time,tuple(F((-1)**j,j+1) for j in range(10)),F(0))
    case=Problem('mixed',F(1),initial,(row,))
    assert dot(row.weights,state(initial,mesh,controls,time))==dot(input_row(case,row,mesh),controls)-adjusted(case)[0]
    t=s.Symbol('t');left,right=F(1,7),F(3,5)
    kk=local_kernel(row,1,left,right)
    expected=(right-left)*sum(s.Rational(w)*(s.Rational(time-left)-s.Rational(right-left)*t)**(4-i)/factorial(4-i)
                             for i,w in enumerate(row.weights[5:]))
    assert s.expand(expected-sum(s.Rational(v)*t**j for j,v in enumerate(kk)))==0


def test_known_anchors_have_exact_supporting_multipliers():
    t=s.Symbol('t')
    roots_by_label={
        'anchor_0':(F(1,5),F(2,5),F(3,5),F(4,5)),
        'anchor_1':(F(1,5),F(2,5),F(3,5),F(4,5)),
        'anchor_2':(F(1,100),F(1,3),F(2,3),F(99,100))}
    for case,known in cases():
        if case.label not in roots_by_label:continue
        roots=tuple(case.horizon*r for r in roots_by_label[case.label])
        polynomial=s.Poly(s.prod(t-s.Rational(r) for r in roots),t)
        wanted=s.Matrix([polynomial.nth(j) for j in range(5)])
        K=s.Matrix([[s.Rational(kernel(row,0)[j]) for row in case.waypoints] for j in range(5)])
        lam=K.inv()*wanted
        mesh=(F(0),*roots,case.horizon)
        norm=sum((-1)**j*s.integrate(polynomial.as_expr(),(t,s.Rational(a),s.Rational(b)))
                 for j,(a,b) in enumerate(zip(mesh,mesh[1:])))
        assert sum(v*s.Rational(rhs) for v,rhs in zip(lam,adjusted(case)))==s.Rational(known)*norm


def test_simple_certificate_and_input_only_reference_scale():
    record=simple_certificate();assert check(record)==(F(1),F(1))
    for case,known in cases():
        case.validate();assert from_record(case.record())==case
        assert reference_scale(case)==checked_scale(case.record())
        if known is not None and known>0:assert reference_scale(case)<=known


@pytest.mark.parametrize('kind',('control','lower','kernel','leaf','waypoint','initial','peak','flag','scale','order'))
def test_checker_rejects_tampering(kind):
    record=deepcopy(simple_certificate())
    if kind=='control':record['controls'][0]='0'
    elif kind=='lower':record['dual']['lower']='2'
    elif kind=='kernel':record['dual']['leaves'][0][-1]='2'
    elif kind=='leaf':record['dual']['leaves']=[]
    elif kind=='waypoint':record['case']['waypoints'][0]['value']='2'
    elif kind=='initial':record['case']['initial'][4]='1'
    elif kind=='peak':record['upper']='0'
    elif kind=='flag':record['converged']=False
    elif kind=='scale':record['scale']='2'
    else:record['case']['order']=4
    with pytest.raises(AssertionError):check(record)


def test_checker_does_not_import_producer_or_numerical_libraries():
    source=Path(__file__).with_name('check.py').read_text(encoding='utf-8')
    tree=ast.parse(source)
    imported=[node.module for node in ast.walk(tree) if isinstance(node,ast.ImportFrom)]
    assert set(imported)=={'fractions','math'}
    assert not any(isinstance(node,ast.Import) for node in ast.walk(tree))


def test_actual_recovery_consumes_only_multipliers_not_hidden_primal_or_known_solution():
    case=simple_case()
    a=certify(case,{'lam':[1.]},max_iterations=1)['certificate']
    b=certify(case,{'lam':[1.],'known_peak':999,'controls':['bad'],'moments':None},max_iterations=1)['certificate']
    assert a==b and a['certified'] and a['converged']


def test_joint_model_compilation_dimensions():
    allcases=cases();assert len(allcases)==15
    for case,_ in (allcases[0],allcases[4],allcases[-1]):
        problem,*_=build_sdp(case);data,_,_=problem.get_problem_data('CLARABEL')
        n=case.axes*(len(knots(case))-1)
        assert data['A'].shape[1]==39*n+1
        assert list(data['dims'].psd)==([3,3,3,2,2,2,2]+[2,2,4,4,3,3,3]*2)*n
