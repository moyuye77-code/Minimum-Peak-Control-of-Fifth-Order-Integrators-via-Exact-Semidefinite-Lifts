"""Zero-seed recovery controls. Frozen synthesis sources are imported unchanged."""
import argparse
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
from statistics import median
import time
from research.peak_synthesis.model import cases,knots
from research.peak_synthesis.certify import certify
from research.peak_synthesis.check import check,check_lower
from research.peak_synthesis.experiment import sources as parent_sources

ROOT=Path(__file__).resolve().parent
OUTPUT=ROOT/'peak_recovery_results'/'verification.json'
PARENT=ROOT/'peak_synthesis'/'results'/'verification.json'
ARMS=(('quarter_mesh_only',1),('zero_seed_adaptive',4))


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    out=parent_sources()
    for name in ('peak_recovery_ablation.py','test_peak_recovery_ablation.py'):
        out[name]=digest(ROOT/name)
    return out


def seed(case):
    # No archive, known optimum, moments, previous control or proposal enters.
    return dict(method='zero_seed',lam=[0.]*len(case.waypoints),optimization_calls=0)


def run(case,iterations):
    assert iterations in (1,4)
    return certify(case,seed(case),max_iterations=iterations)


def audit(data):
    assert data['sources']==sources() and data['parent_sha256']==digest(PARENT)
    parent=json.loads(PARENT.read_text(encoding='utf-8'))
    assert parent['sources']==parent_sources()
    expected=cases(); assert len(data['records'])==len(expected)*len(ARMS)
    assert data['new_closed_loop_runs']==0
    for i,(case,known) in enumerate(expected):
        originals=parent['records'][2*i:2*i+2]
        for j,(name,budget) in enumerate(ARMS):
            row=data['records'][2*i+j]; recovery=row['recovery']; rec=recovery['certificate']
            assert row['case']==rec['case']==case.record()
            assert row['arm']==name and row['iterations']==budget and row['seed']==seed(case)
            assert row['known_peak']==(str(known) if known is not None else None)
            assert 1<=len(recovery['attempts'])==recovery['repair_lp_calls']<=budget
            assert recovery['attempts'][0]['cells']==4*(len(knots(case))-1)
            if rec.get('certified'):
                low,up=check(rec)
                if known is not None: assert low<=known<=up
                # Independently valid intervals for one optimum must intersect.
                for old in originals:
                    assert old['case']==case.record()
                    old_rec=old['recovery']['certificate']
                    if old_rec.get('certified'):
                        a,b=check(old_rec); assert max(a,low)<=min(b,up)
            else: check_lower(row['case'],rec['dual'])
    out={}
    for name,_ in ARMS:
        rows=[r for r in data['records'] if r['arm']==name]
        valid=[r for r in rows if r['recovery']['certificate'].get('certified')]
        out[name]=dict(problems=len(rows),certified=len(valid),
            converged=sum(r['recovery']['certificate'].get('converged',False) for r in valid),
            lp_calls=sum(r['recovery']['repair_lp_calls'] for r in rows),
            max_relative_gap=max((float(F(r['recovery']['certificate']['gap'])/
                                       max(F(1),F(r['recovery']['certificate']['upper'])))
                                  for r in valid),default=None),
            median_seconds=median(r['elapsed_seconds'] for r in rows),
            timing_scope='single run fixed arm order; not comparable randomized timings')
    out['new_closed_loop_runs']=0
    return out


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--create',action='store_true')
    args=parser.parse_args()
    if args.create:
        if OUTPUT.exists(): raise FileExistsError('Preserve ablation archive')
        rows=[]
        for case,known in cases():
            for name,budget in ARMS:
                started=time.perf_counter(); recovery=run(case,budget)
                rec=recovery['certificate']
                if rec.get('certified'): check(rec)
                rows.append(dict(case=case.record(),arm=name,iterations=budget,seed=seed(case),
                    known_peak=str(known) if known is not None else None,
                    recovery=recovery,elapsed_seconds=time.perf_counter()-started))
                print(case.label,name,rec.get('certified'),rec.get('converged'),
                      recovery['repair_lp_calls'],float(F(rec['gap'])) if 'gap' in rec else None,flush=True)
        data=dict(schema=1,sources=sources(),parent_sha256=digest(PARENT),records=rows,
                  new_closed_loop_runs=0)
        audit(data); OUTPUT.parent.mkdir(parents=True,exist_ok=True)
        with OUTPUT.open('x',encoding='utf-8') as stream:
            json.dump(data,stream,indent=2); stream.write('\n')
    data=json.loads(OUTPUT.read_text(encoding='utf-8'))
    print(json.dumps(audit(data),indent=2)); print('archive_sha256',digest(OUTPUT))


if __name__=='__main__': main()
