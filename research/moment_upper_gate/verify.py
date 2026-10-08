import argparse
import hashlib
import json
from pathlib import Path
from .exact import hierarchy_record, prefix_records

ROOT = Path(__file__).resolve().parent
ORDERS = (1,2,3,4,6,8)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    return {p.name:digest(p) for p in sorted(ROOT.glob('*.py'))}


def build():
    return dict(schema=1, hierarchy=[hierarchy_record(R,fixed) for R in ORDERS for fixed in (False,True)],
                prefixes=prefix_records(), sources=sources(),
                new_optimizer_calls=0, new_closed_loop_runs=0,
                full_K4_SOC_lift_established=False)


def audit(data):
    assert data == build(), 'Exact reconstruction or source hash mismatch'
    return dict(hierarchy_witnesses=len(data['hierarchy']), prefix_records=len(data['prefixes']),
                new_optimizer_calls=0,new_closed_loop_runs=0, full_K4_SOC_lift_established=False)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--create',action='store_true')
    args=parser.parse_args()
    path=ROOT/'results'/'verification.json'
    if args.create:
        if path.exists():
            raise FileExistsError('Existing archive is preserved')
        data=build()
        audit(data)
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x',encoding='utf-8') as stream:
            json.dump(data,stream,indent=2)
            stream.write('\n')
    print(json.dumps(audit(json.loads(path.read_text(encoding='utf-8'))),indent=2))
    print('archive_sha256',digest(path))


if __name__ == '__main__':
    main()
