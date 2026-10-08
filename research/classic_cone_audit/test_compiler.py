from copy import deepcopy
from fractions import Fraction as F
import json
import pytest
import sympy as sp
from .compiler import cp, compile_record, primitive_records
from .check import check_primitives, check_record, exact_point_check
from .verify import ROOT, audit, SIZES


@pytest.mark.parametrize('n', SIZES)
def test_automatic_model_uses_only_standard_soc(n, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('This audit must not solve an optimization problem')
    monkeypatch.setattr(cp.Problem, 'solve', forbidden)
    check_record(compile_record(n))


def test_actual_canonicalizer_affine_maps():
    check_primitives(primitive_records())


def test_markov_transform_identity():
    m, z, q = sp.symbols('m z q')
    a0, a1, a2 = m, z-m*m/2, q-m*z+m**3/6
    assert sp.expand(a0*a2-a1*a1-(m*q-z*z-m**4/12)) == 0


@pytest.mark.parametrize('m', [F(0), F(1, 4), F(1, 2), F(3, 4), F(1)])
def test_exact_lifts_including_closed_endpoints(m):
    for z in sorted({m*m/2, m/2, m-m*m/2}):
        lo = z*z/m+m**3/12 if m else F(0)
        hi = F(1, 3)-(F(1, 2)-z)**2/(1-m)-(1-m)**3/12 if m != 1 else F(1, 3)
        for q in sorted({lo, (lo+hi)/2, hi}):
            assert exact_point_check(m, z, q)


def test_false_cone_count_rejected():
    record = compile_record(1)
    record['soc_dimensions'].pop()
    with pytest.raises(AssertionError):
        check_record(record)


def test_changed_power_map_rejected():
    data = primitive_records()
    data['cube']['soc_affine_rows'][1][2][-1] = 1
    with pytest.raises(AssertionError):
        check_primitives(data)


def test_archive():
    data = json.loads((ROOT/'results'/'verification.json').read_text(encoding='utf-8'))
    assert audit(data)['compilations'] == 5
    bad = deepcopy(data)
    bad['new_closed_loop_runs'] = 1
    with pytest.raises(AssertionError):
        audit(bad)
