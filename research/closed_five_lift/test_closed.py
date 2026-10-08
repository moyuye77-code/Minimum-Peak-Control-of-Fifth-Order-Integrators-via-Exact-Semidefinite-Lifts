import pytest
import sympy as s
from fractions import Fraction
from .algebra import model, blocks, evaluate, witness, boundary_cases, complement, X, P, CENTER
from research.moment_upper_gate.exact import psd_pivots
from research.moment_four_lift.exact import matrices, atomic_aux


def test_exact_affine_space_reconstruction():
    data = model()
    substituted = data['coordinates']*s.Matrix(data['basis'])
    assert all(s.expand(a-b)==0 for a,b in zip(substituted,data['source_polynomials']))
    assert len(blocks()) == 7
    assert max(len(v) for _,_,v in blocks()) == 4
    assert all(s.Poly(f,data['theta']).total_degree()<=1 for mat in data['matrices'].values() for f in mat)


@pytest.mark.parametrize('name,point', boundary_cases())
@pytest.mark.parametrize('barred', (False,True))
def test_finite_exact_primal_witness_on_singular_and_regular_points(name,point,barred):
    point = complement(point) if barred else point
    values, valuation = witness(point[:4])
    assert all(v.is_Rational for v in values)
    mats, outputs = evaluate(values)
    assert outputs[0] == 1
    assert outputs[1:5] == point[:4]
    assert outputs[5] <= point[4]
    assert valuation in (0,1,3)
    for mat in mats.values():
        psd_pivots([[Fraction(v) for v in row] for row in mat.tolist()])
    prefix = tuple(map(Fraction,point[:4]))
    for mat in matrices(prefix,atomic_aux(prefix[0],prefix[1])):
        psd_pivots(mat)


def test_zero_prefix_is_actually_feasible_not_an_open_lift_substitution():
    values,valuation = witness((0,0,0,0))
    mats,outputs = evaluate(values)
    assert valuation == 3 and outputs == (1,0,0,0,0,0)
    assert mats['p_C'][-1,-1] == 1
    assert mats['p_h'][-1,-1] == 1
    assert all(f == 0 for name,mat in mats.items() if name not in ('p_C','p_h') for f in mat)


def test_removed_domain_blocks_are_not_silently_retained():
    names = {name for name,_,_ in blocks()}
    assert 'unit_domain' not in names and 'delta_domain' not in names
    assert P.subs(dict(zip(X,(0,0,0,0)))) == 0
    # Positive-mass central point supplies the common feasible point in the proof.
    assert P.subs(dict(zip(X,CENTER[:4]))) > 0


def test_compiled_model_has_only_affine_psd_constraints():
    from .solver import compile_record
    data=compile_record()
    assert data['variables']==39 and data['equalities']==0 and data['nonnegative_rows']==2
    assert data['PSD_blocks']==[3,3,3,2,2,2,2]+[2,2,4,4,3,3,3]*2
    assert data['SOC_blocks']==[] and data['dcp']


def test_prespecified_exact_support_oracle():
    from .cases import cases
    from .experiment import root_oracle
    import numpy as np
    rows=cases()
    assert len(rows)==18
    for case in rows:
        coefs=np.array([float(Fraction(v)) for v in case['coefficients']])
        support,_=root_oracle(coefs)
        assert abs(support-float(Fraction(case['support'])))<1e-10
        point=tuple(s.Rational(v) for v in case['point'])
        assert sum(s.Rational(a)*b for a,b in zip(case['coefficients'],point))==s.Rational(case['support'])
