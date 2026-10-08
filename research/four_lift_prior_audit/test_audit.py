from copy import deepcopy
from fractions import Fraction as F
import pytest
import sympy as sp
from .algebra import (s,t,q,r,S,T,alpha,beta, exponential_coefficients,
                      classical_and_congruent,curvature_identities,generic_matrices,
                      compression_residuals,atomic_matrices,moment_data)
from .verify import build,audit
from research.moment_four_lift.exact import psd_pivots


def test_exponential_coefficients_independently_by_derivatives():
    w=sp.Symbol('w')
    expression=1-sp.exp(-(s*w+t*w*w+q*w**3+r*w**4))
    for degree,a in enumerate(exponential_coefficients(),1):
        assert sp.expand(sp.diff(expression,w,degree).subs(w,0)/sp.factorial(degree)-a)==0


def test_congruence_is_invertible_on_all_boundaries():
    H,C,G=classical_and_congruent()
    assert C.det()==1
    assert sp.expand(H.det()-G.det())==0
    assert G==sp.Matrix([[t-s*s/2,q-s*t/2-s**3/12],
                          [q-s*t/2-s**3/12,r-t*t/2-s*s*t/4]])


def test_hessian_integration_and_homogeneous_certificate():
    ids=curvature_identities()
    assert ids.pop('hessian')==[0]*4
    assert all(v==0 for v in ids.values())


def test_global_sos_convexity_hypothesis_is_false():
    _,_,G=classical_and_congruent()
    # xi=(0,1), s=0,t=-1: a negative Hessian direction. Domain is essential.
    H=sp.hessian(-G[1,1],(s,t)).subs({s:0,t:-1})
    assert H[0,0]==-sp.Rational(1,2)


def test_strict_domain_interior_hypothesis_is_false():
    # Density 1_[0,1/2] lies in K4, but on D=t-s^2/2=0.
    point=(F(1,2),F(1,8),F(1,24),F(1,64))
    assert point[1]-point[0]**2/2==0
    for matrix in atomic_matrices(point):
        psd_pivots(matrix.tolist())


def test_generic_coefficients_are_affine_with_expected_sizes():
    variables=list(moment_data().values())[1:]+[q,r]
    mats=generic_matrices()
    assert [m.rows for m in mats]==[6,3,3,2,2]
    assert len(set(variables))==16
    for mat in mats:
        for entry in mat:
            assert sp.Poly(entry,*variables).total_degree()<=1


def test_all_seven_compressed_blocks_are_congruences_or_identical():
    for matrix in compression_residuals():
        assert matrix==sp.zeros(matrix.rows,matrix.cols)


def test_generic_lift_permits_nondirac_coordinate_functionals():
    a=(F(1,4),F(1,8),F(1,12),F(1,16))
    b=(F(3,4),F(3,8),F(1,4),F(3,16))
    ma,mb=atomic_matrices(a),atomic_matrices(b)
    for left,right in zip(ma,mb):
        psd_pivots(((left+right)/2).tolist())
    assert (a[0]**2+b[0]**2)/2 != ((a[0]+b[0])/2)**2


@pytest.mark.parametrize('field',['novelty_cleared','arbitrary_auxiliary_completion_proved',
                                   'published_identical_K4_result_located'])
def test_archive_rejects_unearned_claims(field):
    data=build()
    data[field]=True
    with pytest.raises(AssertionError):
        audit(data)


def test_archive_rejects_modified_transform():
    data=build()
    data['formulas']['exponential_coefficients'][3]='r'
    with pytest.raises(AssertionError):
        audit(data)
