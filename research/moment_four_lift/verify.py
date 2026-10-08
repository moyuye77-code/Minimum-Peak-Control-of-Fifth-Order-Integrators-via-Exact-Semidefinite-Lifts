import argparse
import hashlib
import json
from pathlib import Path
from .exact import identity_residuals, witness_records, directions

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    return {p.name:digest(p) for p in sorted(ROOT.glob('*.py'))}


def build_exact():
    residuals=identity_residuals()
    assert set(residuals.values())=={'0'}
    return dict(schema=1, sources=sources(), identities=residuals,
                witnesses=witness_records(), directions=directions(),
                full_K4_SDP_lift_proved=True, full_K4_SOC_lift_proved=False,
                max_PSD_block=3, auxiliary_variables=6,
                new_closed_loop_runs=0)


def audit(data):
    assert data==build_exact(), 'Exact reconstruction or source digest mismatch'
    return dict(witnesses=len(data['witnesses']), identities=len(data['identities']),
                support_directions=len(data['directions']), max_PSD_block=3,
                full_K4_SDP_lift_proved=True, full_K4_SOC_lift_proved=False,
                new_closed_loop_runs=0)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--create', action='store_true')
    parser.add_argument('--numeric', action='store_true')
    args=parser.parse_args()
    filename='numerical.json' if args.numeric else 'verification.json'
    path=ROOT/'results'/filename
    if args.create:
        if path.exists():
            raise FileExistsError('Existing evidence is preserved')
        if args.numeric:
            from .solver import run_supports
            data=dict(schema=1, sources=sources(), variants=[run_supports(True),run_supports(False)])
        else:
            data=build_exact()
            audit(data)
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x',encoding='utf-8') as stream:
            json.dump(data,stream,indent=2)
            stream.write('\n')
    data=json.loads(path.read_text(encoding='utf-8'))
    if args.numeric:
        assert data['sources']==sources()
        for variant in data['variants']:
            cases=variant['cases']
            finite=[r for r in cases if r['signed_error'] is not None]
            print(json.dumps(dict(compilation=variant['compilation'],
                                  cases=len(cases), statuses={s:sum(r['status']==s for r in cases) for s in sorted({r['status'] for r in cases})},
                                  finite_results=len(finite),
                                  max_abs_error=max((abs(r['signed_error']) for r in finite),default=None),
                                  min_eigenvalue=min((min(r['minimum_eigenvalues']) for r in finite),default=None),
                                  warnings=sum(len(r['warnings']) for r in cases)),indent=2))
    else:
        print(json.dumps(audit(data),indent=2))
    print('archive_sha256',digest(path))


if __name__=='__main__':
    main()
