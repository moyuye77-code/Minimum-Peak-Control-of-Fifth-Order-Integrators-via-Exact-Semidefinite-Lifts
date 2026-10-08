import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp
from .algebra import (curvature_identities, compression_residuals, formula_records,
                      generic_matrices, atomic_matrices)
from research.moment_four_lift.exact import witness_records, psd_pivots

ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT.parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    paths = list(ROOT.glob('*.py'))
    paths += [RESEARCH/'moment_four_lift'/'exact.py',
              RESEARCH/'moment_upper_gate'/'exact.py']
    return {str(p.relative_to(RESEARCH)):digest(p) for p in sorted(paths)}


def build():
    identities = curvature_identities()
    identities = {k:([str(x) for x in v] if isinstance(v,list) else str(v))
                  for k,v in identities.items()}
    assert identities['hessian'] == ['0']*4
    assert all(v=='0' for k,v in identities.items() if k!='hessian')
    compressed = compression_residuals()
    assert all(mat==sp.zeros(mat.rows,mat.cols) for mat in compressed)
    witnesses = []
    for row in witness_records():
        pivots = [[str(v) for v in psd_pivots(mat.tolist())]
                  for mat in atomic_matrices(row['point'])]
        witnesses.append(dict(label=row['label'], point=row['point'], pivots=pivots))
    return dict(schema=1, sources=sources(), formulas=formula_records(),
                identities=identities, compression_entries=sum(m.rows*m.cols for m in compressed),
                generic_block_sizes=[m.rows for m in generic_matrices()],
                generic_auxiliaries=12, candidate_auxiliaries=6,
                exact_atomic_witnesses=witnesses,
                generic_and_candidate_have_same_K4_projection=True,
                arbitrary_auxiliary_completion_proved=False,
                published_identical_K4_result_located=False,
                novelty_cleared=False, new_optimizer_calls=0, new_closed_loop_runs=0)


def audit(data):
    assert data==build(), 'Evidence, claim, or source hash mismatch'
    return {key:data[key] for key in (
        'compression_entries','generic_block_sizes','generic_auxiliaries',
        'candidate_auxiliaries','generic_and_candidate_have_same_K4_projection',
        'arbitrary_auxiliary_completion_proved','published_identical_K4_result_located',
        'novelty_cleared','new_optimizer_calls','new_closed_loop_runs')}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--create',action='store_true')
    args=parser.parse_args()
    path=ROOT/'results'/'verification.json'
    if args.create:
        if path.exists():
            raise FileExistsError('Preserve archived evidence')
        data=build()
        audit(data)
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x',encoding='utf-8') as stream:
            json.dump(data,stream,indent=2)
            stream.write('\n')
    data=json.loads(path.read_text(encoding='utf-8'))
    print(json.dumps(audit(data),indent=2))
    print('witnesses',len(data['exact_atomic_witnesses']))
    print('archive_sha256',digest(path))


if __name__=='__main__':
    main()
