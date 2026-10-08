import argparse
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
from statistics import median
import time
from .model import cases
from .proposals import sdp,continuous_dual
from .certify import certify
from .check import check,check_lower

ROOT=Path(__file__).resolve().parent
RESEARCH=ROOT.parent


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    paths=list(ROOT.glob('*.py'))
    paths += [RESEARCH/name for name in (
        'moment_four_lift/exact.py','moment_upper_gate/exact.py',
        'four_history_audit/proposals.py','four_history_audit/model.py',
        'classic_cone_audit/compiler.py')]
    return {str(p.relative_to(RESEARCH)):digest(p) for p in sorted(paths)}


def summary(data):
    out={}
    for name in ('homogeneous_sdp','continuous_root_dual'):
        rows=[r for r in data['records'] if r['proposal']['method']==name]
        certified=[r for r in rows if r['recovery']['certificate'].get('certified')]
        statuses=sorted({r['proposal']['status'] for r in rows})
        out[name]=dict(problems=len(rows),certified=len(certified),
            converged=sum(r['recovery']['certificate'].get('converged',False) for r in certified),
            proposal_statuses={s:sum(r['proposal']['status']==s for r in rows) for s in statuses},
            proposal_calls=sum(r['proposal']['optimization_calls'] for r in rows),
            recovery_lp_calls=sum(r['recovery']['repair_lp_calls'] for r in rows),
            max_relative_gap=max((float(F(r['recovery']['certificate']['gap'])/
                                  max(F(1),F(r['recovery']['certificate']['upper'])))
                                  for r in certified),default=None))
    pairs=[]
    for a,b in zip(data['records'][::2],data['records'][1::2]):
        if a['recovery']['certificate'].get('converged') and b['recovery']['certificate'].get('converged'):
            pairs.append((a['elapsed_seconds'],b['elapsed_seconds']))
    out['paired_timing']=dict(pairs=len(pairs),sdp_faster=sum(a<b for a,b in pairs),
        sdp_median=median(a for a,b in pairs) if pairs else None,
        dual_median=median(b for a,b in pairs) if pairs else None,
        scope='one cold run; full proposal+shared recovery+independent check; no random order repetitions')
    out['new_closed_loop_runs']=0
    return out


def audit(data):
    assert data['sources']==sources()
    expected=cases(); assert len(data['records'])==2*len(expected)
    for i,(case,known) in enumerate(expected):
        for j,method in enumerate(('homogeneous_sdp','continuous_root_dual')):
            row=data['records'][2*i+j]
            assert row['case']==case.record() and row['proposal']['method']==method
            assert row['known_peak']==(str(known) if known is not None else None)
            rec=row['recovery']['certificate']; assert rec['case']==case.record()
            assert row['recovery']['repair_lp_calls']==len(row['recovery']['attempts'])
            if rec.get('certified'):
                lower,upper=check(rec)
                if known is not None: assert lower<=known<=upper
            else: check_lower(row['case'],rec['dual'])
    return summary(data)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--create',action='store_true')
    args=parser.parse_args(); path=ROOT/'results'/'verification.json'
    if args.create:
        if path.exists(): raise FileExistsError('Preserve frozen results')
        rows=[]
        for case,known in cases():
            for proposal_function in (sdp,continuous_dual):
                started=time.perf_counter(); proposal=proposal_function(case)
                recovery=certify(case,proposal); certificate=recovery['certificate']
                if certificate.get('certified'): check(certificate)
                rows.append(dict(case=case.record(),known_peak=str(known) if known is not None else None,
                                 proposal=proposal,recovery=recovery,elapsed_seconds=time.perf_counter()-started))
                print(case.label,proposal['method'],proposal['status'],certificate.get('certified'),
                      certificate.get('converged'),recovery['repair_lp_calls'],
                      float(F(certificate['gap'])) if 'gap' in certificate else None,flush=True)
        import numpy,scipy,sympy
        from research.classic_cone_audit.compiler import cp
        import clarabel
        data=dict(schema=1,sources=sources(),records=rows,
                  versions=dict(numpy=numpy.__version__,scipy=scipy.__version__,sympy=sympy.__version__,
                                cvxpy=cp.__version__,clarabel=clarabel.__version__),new_closed_loop_runs=0)
        audit(data); path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x',encoding='utf-8') as stream:
            json.dump(data,stream,indent=2); stream.write('\n')
    data=json.loads(path.read_text(encoding='utf-8'))
    print(json.dumps(audit(data),indent=2)); print('archive_sha256',digest(path))


if __name__=='__main__': main()
