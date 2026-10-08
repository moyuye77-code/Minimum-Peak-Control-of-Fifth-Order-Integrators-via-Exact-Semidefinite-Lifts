import sympy as s
import pytest
from . import algebra as g


@pytest.mark.parametrize('count',[2,4,5,7])
def test_recurrence_against_independent_truncated_exponential(count):
    y=s.symbols('y0:'+str(count));w=s.Symbol('w')
    P=sum(value*w**(j+1) for j,value in enumerate(y))
    # Truncate after each multiplication rather than using the same recurrence.
    power=s.S.One;result=s.S.Zero
    for k in range(1,count+1):
        power=s.Poly(s.expand(power*P),w)
        power=sum(coef*w**ex[0] for ex,coef in power.terms() if ex[0]<=count)
        result+=(-1)**(k+1)*power/s.factorial(k)
    assert tuple(s.expand(result).coeff(w,j+1) for j in range(count))==g.exponential_coefficients(y)


def test_first_row_conditions_and_forced_form():
    _,G=g.raw_family();m,z,q,r,v=g.Y
    assert s.diff(G[0,1],m,2)==2*g.c-1
    equations=[s.Poly(G[0,2]-q,m,z,q).coeff_monomial(ex)
               for ex in (m*z,m**3)]
    assert s.solve(equations,(g.a,g.b))=={g.a:1-g.f,g.b:g.f/2-s.Rational(1,6)}
    C,F=g.normalized_family()
    assert list(F[0,:])==[m,z,q]
    assert F[1,1]==q-m**3/12
    assert s.expand(F[1,2]-r+m*m*z/6-(g.f-s.Rational(1,2))*(m*q-z*z-m**4/12))==0


@pytest.mark.parametrize('F',[s.Rational(1,2),g.m,g.z,g.q,s.Rational(1,2)+g.z*g.q])
def test_arbitrary_polynomial_family_identity(F):
    _,G=g.normalized_family(F);m,z,q,r,v=g.Y
    assert s.expand(G[1,2]-r+m*m*z/6-(F-s.Rational(1,2))*(m*q-z*z-m**4/12))==0
    hz=s.hessian(G[1,2],(z,q,r,v))
    assert (hz==s.zeros(4))==(F==s.Rational(1,2))


def test_local_indefinite_hessian_and_midpoint_are_exact():
    G=g.forced_matrix();m,z,q,r,v=g.Y
    h=s.hessian(G[2,2],(m,q))
    assert s.factor(h.det())==-m*m/36
    w=g.midpoint_witness()
    assert all((a+b)/2==c for a,b,c in zip(w['plus'],w['minus'],w['center']))
    assert w['defect'][2,2]==s.Rational(43,6291456)
    for intercept,slope in [(s.Rational(3,4),s.Rational(-3,8)),
                             (s.Rational(1,4),s.Rational(3,8))]:
        assert all(0<intercept+slope*t<1 for t in (0,1))


@pytest.mark.parametrize('order',[3,4,5])
def test_leading_block_persists_at_larger_orders(order):
    yy=g.Y+s.symbols('extra0:'+str(2*order-6)) if order>3 else g.Y
    H=g.hankel(yy,order)
    C=s.eye(order);C[:3,:3]=g.normalized_family()[0]
    for i in range(3,order):
        for j in range(i):C[i,j]=s.Symbol('c'+str(i)+str(j))
    assert (C*H*C.T)[:3,:3].applyfunc(s.expand)==g.normalized_family()[1]


def test_interval_transform_and_schur_calibrations():
    data=g.records()
    assert [c['d'] for c in data['interval_calibrations']]==[5,6]
    assert data['optimization_calls']==data['new_closed_loop_runs']==0
    assert data['general_SDP_decided'] is False


def test_schur_singularity_is_not_silently_divided():
    # A single initial interval is on the prefix boundary, not in its interior.
    yy=tuple(s.Rational(1,2)**(j+1)/(j+1) for j in range(4))
    with pytest.raises(AssertionError):g.boundary_prediction(yy)


@pytest.mark.parametrize('d',[2,3,4])
def test_general_schur_formula_agrees_with_old_low_order_boundaries(d):
    m,z,q,r,v=g.Y
    D=z-m*m/2;B=q-m*z/2-m**3/12
    expected={2:m*m/2,3:z*z/m+m**3/12,4:z*z/2+m*m*z/4+B*B/D}
    assert s.factor(g.boundary_prediction(g.Y[:d-1])-expected[d])==0


@pytest.mark.parametrize('d',[7,8])
def test_larger_interior_interval_calibrations(d):
    bands=[(s.Rational(1,10),s.Rational(2,10)),
           (s.Rational(4,10),s.Rational(5,10)),
           (s.Rational(7,10),s.Rational(8,10))]
    if d==8:bands=[(0,s.Rational(1,20))]+bands
    yy,aa,ww=g.interval_data(bands,d)
    assert all(w>0 for w in ww)
    assert g.boundary_prediction(yy[:-1])==yy[-1]
