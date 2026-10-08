"""Freeze finite algebra examples; the general theorem is in THEORY.md."""
import argparse
import hashlib
import json
from pathlib import Path
from .exact import witness_record

ROOT = Path(__file__).resolve().parent
ORDERS = tuple(range(3, 13))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    return {p.name: digest(p) for p in sorted(ROOT.glob('*.py'))}


def build():
    return dict(schema=1, orders=[witness_record(d) for d in ORDERS],
                universal_proof='THEORY.md plus the cited Averkov theorem; not finite sampling',
                new_optimizer_calls=0, new_closed_loop_runs=0, sources=sources())


def audit(data):
    assert data == build(), 'Archive differs from exact reconstruction or source hashes'
    rows = data['orders']
    return dict(orders=len(rows), witness_polynomials=sum(len(r['witnesses']) for r in rows),
                exact_evaluations=sum(len(r['grid'])*len(r['witnesses']) for r in rows),
                new_optimizer_calls=0, new_closed_loop_runs=0,
                finite_checks_are_not_the_universal_proof=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    path = ROOT/'results'/'verification.json'
    if args.create:
        if path.exists():
            raise FileExistsError('Existing archive is preserved')
        data = build()
        audit(data)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('x', encoding='utf-8') as stream:
            json.dump(data, stream, indent=2)
            stream.write('\n')
    data = json.loads(path.read_text(encoding='utf-8'))
    print(json.dumps(audit(data), indent=2))
    print('archive_sha256', digest(path))


if __name__ == '__main__':
    main()
