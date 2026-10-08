"""Read-only replay: saved exact data, manuscript matrices and input certificates."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
from unittest.mock import patch
from .experiment import (cp, sha, source_hashes, endpoints, EPSILONS, path_record,
                         check_witness, numeric_residuals, membership_controls,
                         raw_synthesis, attribution, summarize, cases)


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def verify(folder):
    manifest = read(folder/'manifest.json')
    assert manifest['sources'] == source_hashes(), 'Source mismatch'
    hashes = read(folder/'hashes.json')
    assert set(hashes) == {str(p.relative_to(folder)).replace('\\','/')
                          for p in folder.rglob('*.json') if p.name != 'hashes.json'}
    for name, value in hashes.items():
        assert sha(folder/name) == value, name
    assert len(list((folder/'boundary').glob('*.json'))) == 48
    assert len(list((folder/'membership').glob('*.json'))) == 18
    assert len(list((folder/'synthesis').glob('*.json'))) == 15
    for i,(name,end) in enumerate(endpoints()):
        for j,epsilon in enumerate(EPSILONS):
            row = read(folder/'boundary'/f'{i:02d}_{j:02d}.json')
            expected = path_record(name,end,epsilon)
            assert {k:row[k] for k in expected} == expected
            point = tuple(map(F,row['point']))
            assert row['exact_pass'] == check_witness(point,[tuple(map(F,v['values'])) for v in row['sides']])
            check_numeric(point,row['numeric'])
    for i,expected in enumerate(membership_controls()):
        row=read(folder/'membership'/f'{i:02d}.json')
        assert {k:row[k] for k in expected} == expected
        check_numeric(tuple(map(F,row['point'])),row['numeric'])
    for i,(case,known) in enumerate(cases()):
        row=read(folder/'synthesis'/f'{i:02d}.json')
        assert row['case']==case.record()
        assert row['known_peak']==(None if known is None else str(known))
        assert row['raw']==raw_synthesis(case,row['proposal'])
        for route in ('sdp_seed','zero_seed'):
            rec=row['recovery'][route]
            assert rec['repair_lp_calls']==len(rec['attempts'])
            expected=attribution(case,row['proposal'],rec)
            assert row['attribution'][route]==expected
            if known is not None and expected['certified']:
                assert F(expected['lower'])<=known<=F(expected['upper'])
    summary=summarize(folder)
    assert read(folder/'summary.json')==summary
    return summary


def check_numeric(point,row):
    got=numeric_residuals(point,row['candidate'])
    got['diagnostic_pass'] &= row['status'] in ('optimal','optimal_inaccurate')
    assert got==row['residuals']


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with patch.object(cp.Problem,'solve',side_effect=AssertionError('Replay must not call a solver')):
        result=verify(args.output.resolve())
    print(json.dumps(dict(summary=result,new_optimizer_calls=0,verified=True),indent=2))


if __name__=='__main__':
    main()
