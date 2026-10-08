"""Run the fixed pilot and preserve every success, warning and failed check."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
from .model import cases
from .proposals import sdp,root_dual
from .certify import certify
from .check import check,check_upper

ROOT=Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    own={str(p.relative_to(ROOT.parent)):digest(p) for p in sorted(ROOT.glob('*.py'))}
    for folder,name in (('moment_four_lift','exact.py'),('continuous_evidence','support.py')):
        p=ROOT.parent/folder/name
        own[str(p.relative_to(ROOT.parent))]=digest(p)
    return own


def audit(data):
    assert data['sources']==sources(), 'Source digest mismatch'
    expected={case.label:case.record() for case in cases()}
    assert len(data['records'])==2*len(expected)
    assert {(r['case']['label'],r['method']) for r in data['records']}=={
        (label,method) for label in expected for method in ('joint_sdp','continuous_root_dual')}
    for record in data['records']:
        assert record['case']==expected[record['case']['label']]
        certificate=record.get('certificate',{})
        if certificate.get('certified'):
            assert certificate['case']==record['case']
            check(certificate)
        elif 'dual' in certificate:
            check_upper(record['case'],certificate['dual'])
    return summarize(data)


def summarize(data):
    summary={}
    for method in ('joint_sdp','continuous_root_dual'):
        rows=[r for r in data['records'] if r['method']==method]
        valid=[r for r in rows if r.get('certificate',{}).get('certified')]
        conv=[r for r in valid if r['certificate']['converged']]
        failures={}
        for r in rows:
            if not r.get('certificate',{}).get('certified'):
                reason=r.get('certificate',{}).get('reason',r.get('failure','unknown'))
                failures[reason]=failures.get(reason,0)+1
        summary[method]=dict(cases=len(rows),certified=len(valid),converged=len(conv),failures=failures,
            proposal_statuses={status:sum(r.get('proposal',{}).get('status')==status for r in rows)
                               for status in sorted({r.get('proposal',{}).get('status','missing') for r in rows})},
            maximum_certified_gap=max((float(F(r['certificate']['gap'])) for r in valid),default=None),
            proposal_calls=sum(r.get('proposal',{}).get('optimization_calls',0) for r in rows),
            repair_lp_calls=sum(r.get('certificate',{}).get('repair_lp_calls',0) for r in rows))
    pairs=[]
    for a,b in zip(data['records'][::2],data['records'][1::2]):
        assert a['case']['label']==b['case']['label']
        ca,cb=a.get('certificate',{}),b.get('certificate',{})
        if ca.get('certified') and cb.get('certified'):
            assert max(F(ca['lower']),F(cb['lower']))<=min(F(ca['dual']['upper']),F(cb['dual']['upper']))
        if ca.get('converged') and cb.get('converged'):
            pairs.append((a['elapsed_seconds'],b['elapsed_seconds']))
    if pairs:
        from statistics import median
        summary['paired_converged_timing']=dict(pairs=len(pairs),SDP_faster=sum(a<b for a,b in pairs),
            SDP_median=median(a for a,b in pairs),dual_median=median(b for a,b in pairs),
            scope='single cold run per prespecified case; full proposal+certificate+independent check')
    summary['new_closed_loop_runs']=0
    return summary


def run(selected=None):
    records=[]
    for case in (cases() if selected is None else selected):
        for method,proposer in (('joint_sdp',sdp),('continuous_root_dual',root_dual)):
            started=time.perf_counter()
            record=dict(case=case.record(),method=method)
            try:
                proposal=proposer(case); record['proposal']=proposal
                certificate=certify(case,proposal); record['certificate']=certificate
                if certificate.get('certified'): check(certificate)
                elif 'dual' in certificate: check_upper(case.record(),certificate['dual'])
            except Exception as exc:
                # Persist the failure, including checker failures; never turn one into success.
                record['failure']=f'{type(exc).__name__}: {exc}'
                if 'certificate' in record:
                    record['certificate']['certified']=False
                    record['certificate']['reason']='exception_or_independent_check_failed'
            record['elapsed_seconds']=time.perf_counter()-started
            records.append(record)
            c=record.get('certificate',{})
            print(case.label,method,record.get('failure',c.get('reason','checked')),
                  c.get('converged',False),float(F(c['gap'])) if 'gap' in c else None,flush=True)
    return dict(schema=1,sources=sources(),records=records,new_closed_loop_runs=0)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--create',action='store_true')
    parser.add_argument('--name',default='verification.json')
    args=parser.parse_args()
    assert Path(args.name).name==args.name and args.name.endswith('.json')
    path=ROOT/'results'/args.name
    if args.create:
        if path.exists(): raise FileExistsError('Preserve existing experiment')
        data=run()
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x',encoding='utf-8') as stream:
            json.dump(data,stream,indent=2); stream.write('\n')
    data=json.loads(path.read_text(encoding='utf-8'))
    print(json.dumps(audit(data),indent=2))
    print('archive_sha256',digest(path))


if __name__=='__main__': main()
