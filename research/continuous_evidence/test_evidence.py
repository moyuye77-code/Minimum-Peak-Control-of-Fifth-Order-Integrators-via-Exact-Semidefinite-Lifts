from fractions import Fraction as F
from pathlib import Path
import copy
import json
import pytest
from .audit import audit, check_support, check_lp
from .support import Observation, solve

ROOT = Path(__file__).resolve().parent


@pytest.fixture(scope='module')
def data():
    return json.loads((ROOT/'results'/'verification.json').read_text(encoding='utf-8'))


def test_complete_archive(data):
    assert audit(data)['status'] == 'passed'


@pytest.mark.parametrize('index', range(44))
def test_support_certificate(data, index):
    result = check_support(data['supports'][index]['certificate'])
    assert result['gap'] <= 1e-5


@pytest.mark.parametrize('index', range(8))
def test_exact_baseline_lp(data, index):
    check_lp(data['derivative_lp'][index]['certificate'])


def test_no_evidence_closed_form():
    r = solve([F(0)]*3, [F(0)]*3, F(1), [], F(0), F(1))
    check_support(r)
    assert F(r['lower']) == F(r['dual']['upper']) == F(1, 6)


def test_coarse_inner_mesh_is_not_a_continuous_upper_bound():
    w = (F(1), F(-1, 4), F(0))
    r = solve([F(0)]*3, [F(0)]*3, F(1), [], F(0), F(1), w, max_iterations=1)
    check_support(r)
    assert F(r['lower']) == F(1, 24)
    assert F(r['dual']['upper']) >= F(1, 16)
    assert not r['converged']
    refined = solve([F(0)]*3, [F(0)]*3, F(1), [], F(0), F(1), w)
    check_support(refined)
    assert refined['converged'] and F(refined['lower']) == F(1, 16)


def test_undelivered_packet_rejected():
    packet = Observation(F(1, 2), 2, F(0), F(1, 10), F(2))
    with pytest.raises(AssertionError, match='Undelivered'):
        solve([F(0)]*3, [F(1)]*3, F(1), [packet], F(1), F(1))


def test_reversed_timestamp_rejected():
    packet = Observation(F(1), 0, F(0), F(1, 10), F(1, 2))
    with pytest.raises(AssertionError):
        solve([F(0)]*3, [F(1)]*3, F(1), [packet], F(1), F(1))


def test_infeasible_inner_mesh_not_a_physical_infeasibility_verdict():
    packets = [Observation(F(0), 0, F(2), F(0), F(0))]
    with pytest.raises(RuntimeError, match='inconclusive'):
        solve([F(0)]*3, [F(0)]*3, F(1), packets, F(0), F(1), max_iterations=2)


@pytest.mark.parametrize('tamper', ['initial', 'gap', 'coverage', 'polynomial', 'measurement', 'multiplier'])
def test_corrupt_certificate_rejected(data, tamper):
    r = copy.deepcopy(data['supports'][0]['certificate'])
    if tamper == 'initial': r['trajectory'][0] = '999'
    elif tamper == 'gap': r['gap'] = '-1'
    elif tamper == 'coverage': r['dual']['leaves'].pop()
    elif tamper == 'polynomial': r['dual']['leaves'][0]['polynomial'][0] = '999'
    elif tamper == 'measurement': r['observations'][0]['value'] = '999'
    else: r['dual']['multipliers'][0] = '999'
    with pytest.raises(AssertionError):
        check_support(r)


def test_conditional_improvement_not_claimed_universal(data):
    gaps = [F(x['lp_width'])-F(x['certified_width']) for x in data['comparisons']]
    assert gaps[0] > F(3, 100)
    assert gaps[1:] == [0, 0, 0]


def test_preselected_payload_comparison(data):
    assert [p['payload'] for p in data['plans']] == [None, 0, 1, 2]
    assert [p['fixed_reply_packets'] for p in data['plans']] == [0, 3, 3, 3]
    values = [F(p['intervention']) for p in data['plans']]
    assert values[0] > values[1] > values[2] > values[3]


def test_start_of_first_hold_cannot_be_omitted(data):
    bad = copy.deepcopy(data)
    bad['plans'][0]['starts'][0]['offset'] = '-100'
    with pytest.raises(AssertionError): audit(bad)


def test_manifest_completeness(data):
    bad = copy.deepcopy(data)
    bad['source_manifest'].pop('support.py')
    with pytest.raises(AssertionError): audit(bad)


def test_missing_plan_start_rejected(data):
    bad = copy.deepcopy(data)
    bad['plans'][0]['starts'].pop()
    with pytest.raises(AssertionError): audit(bad)


def test_missing_unimproved_history_rejected(data):
    bad = copy.deepcopy(data)
    bad['comparisons'].pop()
    with pytest.raises(AssertionError): audit(bad)
