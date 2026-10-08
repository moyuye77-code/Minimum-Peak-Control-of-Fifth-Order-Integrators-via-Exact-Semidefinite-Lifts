"""Reproduce a small compiler audit, not a solver or closed-loop benchmark."""
import argparse
import hashlib
import json
from pathlib import Path
from .compiler import cp, RESEARCH, compile_record, primitive_records
from .check import check_record, check_primitives

ROOT = Path(__file__).resolve().parent
SIZES = (1, 2, 4, 8, 16)
DEPENDENCIES = [
    '.deps/sls_baseline/cvxpy/reductions/dcp2cone/canonicalizers/power_canon.py',
    '.deps/sls_baseline/cvxpy/reductions/dcp2cone/canonicalizers/quad_over_lin_canon.py',
    '.deps/sls_baseline/cvxpy/utilities/power_tools.py',
    '.deps/sls_baseline/cvxpy/atoms/elementwise/power.py',
    'moment_information/moments.py', 'moment_soc/solver.py']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    paths = list(ROOT.glob('*.py'))+[RESEARCH/rel for rel in DEPENDENCIES]
    return {str(path.relative_to(RESEARCH)).replace('\\', '/'): digest(path)
            for path in sorted(paths)}


def build():
    return dict(schema=1, cvxpy_version=cp.__version__, primitives=primitive_records(),
                compilations=[compile_record(n) for n in SIZES],
                new_optimizer_calls=0, new_closed_loop_runs=0, sources=sources())


def audit(data, recompile=False):
    assert data['schema'] == 1 and data['cvxpy_version'] == cp.__version__ == '1.7.5'
    assert data['new_optimizer_calls'] == data['new_closed_loop_runs'] == 0
    assert data['sources'] == sources()
    check_primitives(data['primitives'])
    assert [row['segments'] for row in data['compilations']] == list(SIZES)
    for row in data['compilations']:
        check_record(row)
    if recompile:
        assert data['primitives'] == primitive_records()
        assert data['compilations'] == [compile_record(n) for n in SIZES]
    return dict(compilations=5, primitive_mappings=3, soc_per_segment=8,
                manual_variables='12N+3', automatic_variables='14N+3',
                new_optimizer_calls=0, new_closed_loop_runs=0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--create', action='store_true')
    parser.add_argument('--recompile', action='store_true')
    args = parser.parse_args()
    path = ROOT/'results'/'verification.json'
    if args.create:
        if path.exists():
            raise FileExistsError('Existing archive is preserved')
        data = build()
        audit(data)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    data = json.loads(path.read_text(encoding='utf-8'))
    print(json.dumps(audit(data, recompile=args.recompile), indent=2))
    print('archive_sha256', digest(path))


if __name__ == '__main__':
    main()
