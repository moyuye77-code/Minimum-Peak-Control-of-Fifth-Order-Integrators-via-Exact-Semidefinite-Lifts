from copy import deepcopy
from fractions import Fraction as F
from math import factorial
import pytest
import sympy as sp
from .model import (Problem,Waypoint,cases,homogeneous,input_row,state,dot,adjusted,
                    local_kernel,kernel,knots,from_record)
from .proposals import build_sdp,sdp,continuous_dual
from .certify import dual_bound,certify
from .check import check,check_lower,evolve
from research.moment_four_lift.exact import matrices,atomic_aux,psd_pivots


def simple_case():
    return Problem('simple',F(1),(F(0),)*4,
                   (Waypoint(F(1),(F(0),F(0),F(0),F(1)),F(1)),))


def simple_certificate():
    c=simple_case(); dual=dual_bound(c,[F(1)])
    return dict(case=c.record(),mesh=['0','1'],controls=['1'],upper='1',lower='1',
                gap='0',dual=dual,tolerance='1/1000000',converged=True,certified=True)


def test_homogenization_exact_scaling_and_affineness():
    J=sp.Symbol('J'); y=sp.symbols('m z q r'); aux=sp.symbols('a b c d e f')
    left=homogeneous(J,y,aux)
    right=matrices([v/J for v in y],[v/J for v in aux])
    for a,b in zip(left,right):
        assert (sp.Matrix(a)-J*sp.Matrix(b)).applyfunc(sp.simplify)==sp.zeros(len(a))
        for row in a:
            for value in row: assert sp.Poly(value,J,*y,*aux).total_degree()<=1


@pytest.mark.parametrize('J',[F(0),F(1,10000),F(2),F(1000)])
def test_scaled_atomic_witness_including_zero(J):
    y=(F(1,2),F(1,4),F(1,6),F(1,8)); aux=atomic_aux(*y[:2])
    for matrix in homogeneous(J,[J*v for v in y],[J*v for v in aux]): psd_pivots(matrix)


def test_rest_to_rest_analytic_value_without_float_roots():
    t=sp.Symbol('t'); a=(1-1/sp.sqrt(2))/2; b=sp.Rational(1,2); c=(1+1/sp.sqrt(2))/2
    mesh=(sp.S(0),a,b,c,sp.S(1)); signs=(1,-1,1,-1)
    state=[sp.simplify(sum(sign*384*sp.integrate((1-t)**(3-i)/sp.factorial(3-i),(t,l,r))
                          for l,r,sign in zip(mesh,mesh[1:],signs))) for i in range(4)]
    assert state==[1,0,0,0]
    p=sum(v*(1-t)**(3-i)/sp.factorial(3-i) for i,v in enumerate((384,-192,40,-4)))
    assert sp.expand(p+4*(2*t-1)*(8*t*t-8*t+1))==0
    norm=sp.simplify(sum(sign*sp.integrate(p,(t,l,r)) for l,r,sign in zip(mesh,mesh[1:],signs)))
    assert norm==1


def test_local_kernel_matches_direct_symbolic_integral():
    t=sp.Symbol('s'); o=Waypoint(F(1),(F(1),F(-2),F(3),F(-4)),F(0))
    l,r=F(1,5),F(3,4); kk=local_kernel(o,0,l,r)
    expression=(sp.Rational(r-l)*sum(sp.Rational(w)*(1-sp.Rational(l)-sp.Rational(r-l)*t)**(3-i)
                                    /sp.factorial(3-i) for i,w in enumerate(o.weights)))
    assert sp.expand(expression-sum(sp.Rational(k)*t**i for i,k in enumerate(kk)))==0


def test_independent_dynamics_and_mixed_rows():
    initial=tuple(F(i,11) for i in range(8)); mesh=list(map(F,(0,F(1,4),F(2,3),1)))
    controls=tuple(F(i-2,3) for i in range(6)); time=F(4,5)
    a=state(initial,mesh,controls,time); b=evolve(initial,mesh,controls,time)
    assert a==b
    o=Waypoint(time,tuple(F((-1)**i,i+1) for i in range(8)),F(0))
    case=Problem('mixed',F(1),initial,(o,))
    assert dot(o.weights,a)==dot(input_row(case,o,mesh),controls)-adjusted(case)[0]


def test_simple_exact_certificate_and_scale_invariant_dual():
    rec=simple_certificate(); assert check(rec)==(F(1),F(1))
    for scale in (F(1,100),F(2),F(1000)):
        assert check_lower(rec['case'],dual_bound(simple_case(),[scale]))==1


@pytest.mark.parametrize('kind',['control','lower','kernel','leaf','waypoint','initial','peak','flag'])
def test_checker_rejects_tampering(kind):
    rec=deepcopy(simple_certificate())
    if kind=='control': rec['controls'][0]='0'
    elif kind=='lower': rec['dual']['lower']='2'
    elif kind=='kernel': rec['dual']['leaves'][0][3]='2'
    elif kind=='leaf': rec['dual']['leaves']=[]
    elif kind=='waypoint': rec['case']['waypoints'][0]['value']='2'
    elif kind=='initial': rec['case']['initial'][3]='1'
    elif kind=='peak': rec['upper']='0'
    else: rec['converged']=False
    with pytest.raises(AssertionError): check(rec)


def test_declared_cases_compile_and_roundtrip_without_solving():
    allcases=cases(); assert len(allcases)==18
    for case,known in allcases:
        assert from_record(case.record())==case
        case.validate()
    for case,_ in (allcases[0],allcases[10],allcases[14]):
        p,*_=build_sdp(case); data,_,_=p.get_problem_data('CLARABEL')
        dims=data['dims']; n=case.axes*(len(knots(case))-1)
        assert list(dims.psd)==[3,3,3,2,2,2,2]*n
        assert data['A'].shape[1]==10*n+1


@pytest.mark.parametrize('method',[sdp,continuous_dual])
def test_two_proposal_smoke_certificates(method):
    case=simple_case(); proposed=method(case); result=certify(case,proposed)
    cert=result['certificate']; assert cert['certified'] and cert['converged']
    lo,up=check(cert); assert lo<=1<=up
