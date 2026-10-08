from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import pytest
import sympy as sp
from .exact import (matrices, atomic_aux, identity_residuals, witness_records,
                    directions, interval_point, psd_pivots, fourth_bounds)
from .verify import build_exact, audit
from .solver import build, compile_record


def test_exact_identities():
    assert set(identity_residuals().values())=={'0'}


def test_all_witnesses():
    assert len(witness_records())==76
    for row in witness_records():
        for matrix in matrices(tuple(map(F,row['point'])), tuple(map(F,row['auxiliary'])),False):
            psd_pivots(matrix)


@pytest.mark.parametrize('compressed,blocks',[(True,[3,3,3,2,2,2,2]),(False,[4,2,2,2,2])])
def test_compilation_is_affine_and_has_no_hidden_domain(compressed,blocks):
    rec=compile_record(compressed)
    assert rec['psd_dimensions']==blocks
    assert rec['variables']==10
    assert rec['soc_dimensions']==[]
    assert rec['equalities']==rec['nonnegative_rows']==0


def test_pencils_are_affine():
    yy=sp.symbols('m z q r')
    aa=sp.symbols('a b c d e f')
    for matrix in matrices(yy,aa):
        for row in matrix:
            for entry in row:
                assert sp.Poly(entry,*(yy+aa)).total_degree()<=1


def test_complement_pencils_match_direct_substitution():
    m,z,q,r,a,b,c,d,e,f=sp.symbols('m z q r a b c d e f')
    original=matrices((m,z,q,r),(a,b,c,d,e,f))
    complement=(1-m,sp.Rational(1,2)-z,sp.Rational(1,3)-q,sp.Rational(1,4)-r)
    abar=(1-2*m+a,1-3*m+3*a-b,1-4*m+6*a-4*b+c,
          sp.Rational(1,2)-z-m/2+d,
          sp.Rational(1,2)-m+a/2-z+2*d-e,
          sp.Rational(1,4)-z+f)
    direct=sp.Matrix(matrices(complement,abar)[5])
    assert (direct-sp.Matrix(original[6])).applyfunc(sp.expand)==sp.zeros(2)
    # [1,1-m] is a basis change from [1,m].
    C=sp.Matrix([[1,0],[1,-1]])
    directT=sp.Matrix(matrices(complement,abar)[3])
    assert (directT-C*sp.Matrix(original[4])*C.T).applyfunc(sp.expand)==sp.zeros(2)


def test_gram_compression_is_congruence_of_full_matrix():
    m,z,q,r,a,b,c,d,e,f=sp.symbols('m z q r a b c d e f')
    yy=(m,z,q,r); aa=(a,b,c,d,e,f)
    full=sp.Matrix(matrices(yy,aa,False)[0])
    compact=matrices(yy,aa)
    changes=[sp.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1]])]
    for sign in (1,-1):
        changes.append(sp.Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,sp.Rational(sign,6)]]))
    for mat,C in zip(compact,changes):
        assert (sp.Matrix(mat)-C*full*C.T).applyfunc(sp.expand)==sp.zeros(3)


def test_independent_exact_support_by_antiderivative():
    t=sp.Symbol('t')
    for row in directions():
        coeff=list(map(sp.Rational,row['direction']))
        polynomial=sum(c*t**j for j,c in enumerate(coeff))
        integrated=sum(sp.integrate(polynomial,(t,sp.Rational(lo),sp.Rational(hi))) for lo,hi in row['intervals'])
        assert integrated==sp.Rational(row['support'])


@pytest.mark.parametrize('point',[(F(1,2),F(1,4),F(1,6),F(15,128)-F(1,100)),
                                  (F(1,2),F(1,4),F(1,6),F(17,128)+F(1,100)),
                                  (F(0),F(0),F(0),F(1,10)),
                                  (F(1),F(1,2),F(1,3),F(1,5))])
def test_atomic_invalid_fourth_moment_is_rejected(point):
    # This checks a witness, not infeasibility of all non-atomic lifts.
    with pytest.raises(AssertionError):
        for matrix in matrices(point,atomic_aux(*point[:2])):
            psd_pivots(matrix)


def test_archive_rejects_strengthened_claim_or_modified_witness():
    data=build_exact()
    audit(data)
    bad=deepcopy(data); bad['full_K4_SOC_lift_proved']=True
    with pytest.raises(AssertionError): audit(bad)
    bad=deepcopy(data); bad['witnesses'][0]['point'][3]='123'
    with pytest.raises(AssertionError): audit(bad)


def test_nonlinear_pair_matches_previous_characterization_on_grid():
    grids=[map(F,values) for values in (
        ('-1','0','1/4','1/2','1','5/4'),
        ('-1/8','0','1/8','1/4','1/2','5/8'),
        ('-1/10','0','1/8','1/6','1/4','1/3','1/2'),
        ('-1/10','0','1/10','1/8','1/4','2/5'))]
    checked=0
    for point in product(*grids):
        try:
            lo,hi=fourth_bounds(*point[:3])
            original=lo<=point[3]<=hi
        except ValueError:
            original=False
        try:
            for matrix in matrices(point,atomic_aux(*point[:2]))[-2:]:
                psd_pivots(matrix)
            pair=True
        except AssertionError:
            pair=False
        assert pair==original
        checked+=1
    assert checked==1512


def test_nonatomic_feasible_lift_and_projected_point():
    x=interval_point([(F(0),F(1,4)),(F(1,2),F(3,4))])
    y=interval_point([(F(1,5),F(4,5))])
    weight=F(2,5)
    point=tuple(weight*v+(1-weight)*w for v,w in zip(x,y))
    aux=tuple(weight*v+(1-weight)*w for v,w in zip(atomic_aux(*x[:2]),atomic_aux(*y[:2])))
    assert aux!=atomic_aux(*point[:2])
    for matrix in matrices(point,aux):
        psd_pivots(matrix)
    lo,hi=fourth_bounds(*point[:3])
    assert lo<=point[3]<=hi
