"""Exclusive archive creation and read-only replay of the sparse lift evidence."""
import argparse
import hashlib
import json
from pathlib import Path
from .algebra import records

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'results'/'verification.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    paths = sorted(ROOT.glob('*.py'))+[ROOT/'THEORY.md']
    dependencies = (
        'schur_remainder_recursion/algebra.py',
        'schur_remainder_recursion/THEORY.md',
        'schur_remainder_recursion/results/verification.json',
        'moment_four_lift/THEORY.md',
        'moment_lift_obstruction/THEORY.md',
    )
    paths.extend(ROOT.parent/item for item in dependencies)
    return dict(evidence=records(), sources={str(p.relative_to(ROOT.parent)):digest(p) for p in paths})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    data = build()
    if args.write:
        OUTPUT.parent.mkdir(exist_ok=True)
        with OUTPUT.open('x', encoding='utf-8') as stream:
            json.dump(data, stream, indent=2)
            stream.write('\n')
    assert json.loads(OUTPUT.read_text(encoding='utf-8')) == data
    evidence = data['evidence']
    print('sparse certificate replay passed; 3 <= sxdeg(K5) <= 4 uses written closure proof')
    print('functional_coordinates', evidence['scalar_functional_coordinates'])
    print('open_blocks', [item['order'] for item in evidence['blocks']])
    print('archive_sha256', digest(OUTPUT))


if __name__ == '__main__':
    main()
