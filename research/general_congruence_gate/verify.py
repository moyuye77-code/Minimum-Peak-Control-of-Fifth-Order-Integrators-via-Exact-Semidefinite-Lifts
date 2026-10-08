"""Freeze/replay symbolic and rational evidence; never run an optimizer."""
import argparse
import hashlib
import json
from pathlib import Path
from .algebra import records

ROOT=Path(__file__).resolve().parent
OUTPUT=ROOT/'results'/'verification.json'
SOURCE=ROOT.parent/'paper'/'sources'/'Gosse_Runborg_0809.3714.pdf'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    paths=sorted(ROOT.glob('*.py'))+[ROOT/'THEORY.md',ROOT/'PRIOR_ART.md',SOURCE]
    return dict(evidence=records(),sources={str(p.relative_to(ROOT.parent)):digest(p) for p in paths})


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    data=build()
    if args.write:
        OUTPUT.parent.mkdir(exist_ok=True)
        with OUTPUT.open('x',encoding='utf-8') as f:json.dump(data,f,indent=2);f.write('\n')
    assert json.loads(OUTPUT.read_text(encoding='utf-8'))==data
    print('symbolic/rational replay passed; general SDP representation not decided')
    print('archive_sha256',digest(OUTPUT))


if __name__=='__main__':main()
