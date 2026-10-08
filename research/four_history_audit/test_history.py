from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import pytest
import sympy as sp
from .model import (cases,kernel_poly,integrated_kernel,bernstein,state,trajectory_row,
                    dot,dual_certificate,Observation)
from .proposals import sdp,root_dual,numeric_integral
from .certify import certify
from .check import check,check_upper,evolve


def test_data_and_causality():
    data=cases(); assert len(data)==65 and len({c.label for c in data})==65
    for case in data: case.validate()
    c=data[0]
    late=replace(c.observations[0],received=c.now+1)
    with pytest.raises(AssertionError):
        replace(c,observations=(late,)).validate()


def test_cubic_integral_and_bernstein_independently():
    t=sp.Symbol('t')
    w=(F(2),F(-3),F(4),F(5))
    p=kernel_poly(F(7,5),w)
    poly=sum(sp.Rational(v)*t**j for j,v in enumerate(p))
    l,r=F(1,7),F(3,4)
    assert sp.integrate(poly,(t,sp.Rational(l),sp.Rational(r)))==sp.Rational(integrated_kernel(F(7,5),w,l,r))
    u=sp.Symbol('u')
    bs=bernstein(p,l,r)
    bern=sum(sp.Rational(bs[j])*sp.binomial(3,j)*u**j*(1-u)**(3-j) for j in range(4))
    assert sp.expand(poly.subs(t,sp.Rational(l)+sp.Rational(r-l)*u)-bern)==0


def test_trajectory_two_independent_implementations():
    x=(F(1,3),F(-2),F(1),F(-1,2))
    mesh=[F(0),F(1,3),F(3,4),F(1)]
    u=[F(-1),F(1,2),F(1)]
    w=(F(1),F(2),F(-3),F(4))
    for t in (F(0),F(1,5),F(1,3),F(4,5),F(1)):
        assert state(x,mesh,u,t)==evolve(x,mesh,u,t)
        assert dot(trajectory_row(t,w,mesh),list(x)+u)==dot(w,state(x,mesh,u,t))


@pytest.mark.parametrize('ell',[F(1,100),F(1,10),F(1,2),F(1)])
def test_analytic_extremal_dual_and_trajectory(ell):
    c=next(c for c in cases() if c.label==f'boundary_ell{ell}')
    lam=[-6/ell**3]
    dual=dual_certificate(c,lam)
    upper=check_upper(c.record(),dual)
    assert 0<=upper-c.truth_support<=F(1,10**9)
    mesh=[F(0),1-ell,F(1)] if ell<1 else [F(0),F(1)]
    controls=[F(1),F(-1)] if ell<1 else [F(-1)]
    xx=evolve([F(0)]*4,mesh,controls,F(1))
    assert -xx[-1]==c.truth_support
    assert xx[0]==F(1,24)-ell**4/12


def test_numeric_cubic_integral_repeated_roots():
    for p in ([F(-1,8),F(3,4),F(-3,2),F(1)],
              [F(1,4),F(-1),F(1),F(0)],[F(0)]*4):
        val,_=numeric_integral(p,0.,1.)
        c=replace(cases()[0],observations=(),radii=(F(0),)*4,center=(F(0),)*4,
                  target=F(1),now=F(1),weights=(F(-6)*p[3],2*p[2]+6*p[3],
                      -p[1]-2*p[2]-3*p[3],sum(p)))
        assert kernel_poly(c.target,c.weights)==list(p)
        up=check_upper(c.record(),dual_certificate(c,[]))
        assert abs(val-float(up))<1e-8


@pytest.fixture(scope='module')
def boundary_records():
    case=next(c for c in cases() if c.label=='boundary_ell1/2')
    records=[]
    for solver in (sdp,root_dual):
        proposal=solver(case)
        certified=certify(case,proposal)
        assert certified['certified']
        check(certified)
        records.append(certified)
    return records


def test_both_solvers_and_common_verifier(boundary_records):
    assert all(r['converged'] for r in boundary_records)


@pytest.mark.parametrize('corruption',['trajectory','upper','arrival','polynomial','gap'])
def test_independent_checker_rejects_corruption(boundary_records,corruption):
    bad=deepcopy(boundary_records[0])
    if corruption=='trajectory': bad['trajectory'][-1]='2'
    elif corruption=='upper': bad['dual']['upper']='-100'
    elif corruption=='arrival': bad['case']['observations'][0]['received']='2'
    elif corruption=='polynomial': bad['dual']['leaves'][0][-1]='1234'
    elif corruption=='gap': bad['gap']='0'
    with pytest.raises(AssertionError): check(bad)
