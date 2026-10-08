"""Freeze and replay exact construction evidence, with no solver calls."""
import argparse
import hashlib
import json
from pathlib import Path
from .algebra import records

ROOT=Path(__file__).resolve().parent
OUTPUT=ROOT/'results'/'verification.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    dependencies = ('general_congruence_gate/algebra.py',
                    'general_congruence_gate/THEORY.md',
                    'moment_four_lift/THEORY.md',
                    'moment_lift_obstruction/THEORY.md')
    paths = sorted(ROOT.glob('*.py'))+[ROOT/'THEORY.md',ROOT/'PRIOR_ART.md']
    paths += [ROOT.parent/d for d in dependencies]
    return dict(evidence=records(),
                sources={str(p.relative_to(ROOT.parent)):digest(p) for p in paths})


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    data=build()
    if args.write:
        OUTPUT.parent.mkdir(exist_ok=True)
        with OUTPUT.open('x',encoding='utf-8') as f:
            json.dump(data,f,indent=2)
            f.write('\n')
    assert json.loads(OUTPUT.read_text(encoding='utf-8'))==data
    print('exact weighted certificate replay passed; no general induction or priority claim')
    print('archive_sha256',digest(OUTPUT))


if __name__=='__main__':
    main()
