from copy import deepcopy
from fractions import Fraction as F
import pytest
import sympy as sp
from .exact import (fourth_bounds, lower_endpoints, fiber_matrices, psd_pivots,
                    hierarchy_record, prefix_records, lower_fourth)
from .verify import ORDERS, build, audit


@pytest.mark.parametrize('R',ORDERS)
def test_exact_hierarchy_witnesses(R):
    pure=hierarchy_record(R)
    fixed=hierarchy_record(R,True)
    assert F(pure['gap_lower_bound'])==F(1,2*(R+1)**4)
    assert F(fixed['true_support'])==F(1,8)
    # Independent orthogonal-polynomial evaluation of the inverse-Gram scalar.
    kernel=sum((2*j+1)*sp.legendre(j,-3)**2 for j in range(R+1))
    assert F(fixed['mass_added'])==F(1,2*int(kernel))
    assert F(fixed['gap_lower_bound'])>0


def test_all_prefix_endpoints_and_integrals():
    t=sp.Symbol('t')
    for row in prefix_records():
        m,z,q=map(F,(row['m'],row['z'],row['q']))
        for mm,zz,qq in ((m,z,q),(1-m,F(1,2)-z,F(1,3)-q)):
            a,b,c=lower_endpoints(mm,zz,qq)
            assert bool(a>=0) and bool(b>=a) and bool(c>=b) and bool(c<=1)
            computed=[sp.simplify(sp.integrate(t**j,(t,0,a))+sp.integrate(t**j,(t,b,c))) for j in range(4)]
            assert computed==[sp.Rational(x) for x in (mm,zz,qq,lower_fourth(mm,zz,qq))]
        lo,hi=fourth_bounds(m,z,q)
        if 0<m<1:
            for r in sorted({lo,(lo+hi)/2,hi}):
                matrices,slacks=fiber_matrices(m,z,q,r)
                assert len(matrices)==6 and min(slacks)>=0
                for matrix in matrices:
                    psd_pivots(matrix)


def test_symbolic_prefix_and_fourth_identities():
    m,D,q=sp.symbols('m D q')
    z=D+m*m/2
    upper=sp.Rational(1,3)-(sp.Rational(1,2)-z)**2/(1-m)-(1-m)**3/12
    assert sp.factor(upper-(m**3/3+(m+1)*D-D**2/(1-m)))==0
    Q=q-m**3/3
    markov=m*q+z*z/2-m*m*z/2+m**4/24+(q-m*z+m**3/6)**2/D
    assert sp.factor(markov-(Q*Q/D-m*Q+m*m*D+D*D/2+m**4/4))==0


def test_psd_checker_rejects_zero_pivot_offdiagonal():
    with pytest.raises(AssertionError):
        psd_pivots([[0,1],[1,1]])
    with pytest.raises(AssertionError):
        psd_pivots([[1,2],[2,1]])
    assert psd_pivots([[1,1],[1,1]])==[1,0]


def test_invalid_fourth_moment_and_endpoint_prefix():
    lo,hi=fourth_bounds(F(1,2),F(1,4),F(1,6))
    assert (lo,hi)==(F(15,128),F(17,128))
    with pytest.raises(ValueError):
        fiber_matrices(F(1,2),F(1,4),F(1,6),lo-F(1,100))
    with pytest.raises(ValueError):
        lower_fourth(F(0),F(0),F(1))


def test_false_general_lift_and_corrupt_certificate_rejected():
    data=build()
    audit(data)
    bad=deepcopy(data)
    bad['full_K4_SOC_lift_established']=True
    with pytest.raises(AssertionError):
        audit(bad)
    bad=deepcopy(data)
    bad['hierarchy'][0]['gap_lower_bound']='0'
    with pytest.raises(AssertionError):
        audit(bad)
