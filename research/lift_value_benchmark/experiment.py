"""Immutable per-query records; serial, resumable matched benchmark."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import platform
from statistics import median
import time

from .cases import suite,schedule,record,METHODS
from .formulations import sdp,continuous_dual,templates
from research.peak_synthesis.certify import certify
from research.peak_synthesis.check import check,check_lower,evolve

ROOT=Path(__file__).resolve().parent
RESEARCH=ROOT.parent
RESULTS=ROOT/'results'
MANIFEST=RESULTS/'manifest.json'
OUTPUT=RESULTS/'verification.json'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    paths=list(ROOT.glob('*.py'))+[ROOT/'PROTOCOL.md']
    paths += [RESEARCH/name for name in ('four_lift_prior_audit/algebra.py',
        'moment_four_lift/exact.py','moment_upper_gate/exact.py',
        'classic_cone_audit/compiler.py','four_history_audit/proposals.py',
        'four_history_audit/model.py','peak_synthesis/model.py',
        'peak_synthesis/certify.py','peak_synthesis/check.py')]
    return {str(p.relative_to(RESEARCH)):digest(p) for p in sorted(paths)}


def write_new(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def audit_manifest(data):
    assert data['sources']==sources(),'Do not resume after source/protocol changes'
    expected=[record(item) for item in suite()]
    assert data['cases']==expected and data['schedule']==schedule()
    for item in expected:
        case=item['case'];w=item['witness'];initial=tuple(map(F,case['initial']))
        mesh=tuple(map(F,w['mesh']));controls=tuple(map(F,w['controls']))
        assert mesh[0]==0 and mesh[-1]==F(case['horizon'])
        assert F(w['peak'])==max(map(abs,controls))
        for o in case['waypoints']:
            x=evolve(initial,mesh,controls,F(o['time']))
            assert sum(F(a)*b for a,b in zip(o['weights'],x))==F(o['value'])


def assess(item,recovery):
    c=recovery['certificate'];A=F(item['amplitude'])
    if c.get('certified'):
        lower,upper=check(c)
        if item['known'] is not None:assert lower<=F(item['known'])<=upper
        gap=upper-lower;scaled=gap/max(A,upper)
        return dict(valid=True,legacy_converged=c['converged'],
                    suite_converged=scaled<=F(1,10**6),scaled_gap=str(scaled),
                    lower=str(lower),upper=str(upper))
    check_lower(item['case'],c['dual'])
    return dict(valid=False,legacy_converged=False,suite_converged=False)


def run_trial(item,entry):
    case=item['case'];method=entry['method'];started=time.perf_counter()
    proposal=continuous_dual(case) if method=='continuous_root_dual' else sdp(case,method)
    tolerance=F(1,10**6)*min(F(1),item['amplitude'])
    recstart=time.perf_counter();recovery=certify(case,proposal,tolerance=tolerance)
    recovered=time.perf_counter();assessment=assess(record(item),recovery)
    checked=time.perf_counter()
    return dict(schedule=entry,case=case.record(),proposal=proposal,recovery=recovery,
                assessment=assessment,requested_tolerance=str(tolerance),
                recovery_wall_seconds=recovered-recstart,check_seconds=checked-recovered,
                full_seconds=checked-started)


def audit_trial(item,entry,row):
    assert row['schedule']==entry and row['case']==item['case']
    assert row['proposal']['method']==entry['method']
    assert row['requested_tolerance']==str(F(1,10**6)*min(F(1),F(item['amplitude'])))
    recovery=row['recovery'];assert recovery['certificate']['case']==item['case']
    assert recovery['repair_lp_calls']==len(recovery['attempts'])<=4
    if recovery['certificate'].get('certified'):
        assert recovery['certificate']['tolerance']==row['requested_tolerance']
        assert len(recovery['certificate']['mesh'])<=257
    assert row['assessment']==assess(item,recovery)
    assert row['proposal']['optimization_calls']==1
    for key in ('construction_seconds','canonicalization_seconds','solve_call_seconds','proposal_seconds'):
        assert row['proposal'][key]>=0
    assert row['full_seconds']>=row['proposal']['proposal_seconds']
    return row['assessment']


def summary(rows):
    out={}
    for method in METHODS:
        rr=[r for r in rows if r['schedule']['method']==method]
        valid=[r for r in rr if r['assessment']['valid']]
        statuses=sorted({r['proposal']['status'] for r in rr})
        out[method]=dict(queries=len(rr),valid=len(valid),
            suite_converged=sum(r['assessment']['suite_converged'] for r in rr),
            legacy_converged=sum(r['assessment']['legacy_converged'] for r in rr),
            warnings=sum(len(r['proposal']['warnings']) for r in rr),
            statuses={s:sum(r['proposal']['status']==s for r in rr) for s in statuses},
            recovery_lp_calls=sum(r['recovery']['repair_lp_calls'] for r in rr),
            max_scaled_gap=max((float(F(r['assessment']['scaled_gap'])) for r in valid),default=None),
            medians_seconds={k:median(r['proposal'][k] for r in rr) for k in
                ('construction_seconds','canonicalization_seconds','solve_call_seconds','proposal_seconds')},
            recovery_median=median(r['recovery_wall_seconds'] for r in rr),
            full_median=median(r['full_seconds'] for r in rr))
    paired=[]
    for case_id in range(14):
        group={m:[r for r in rows if r['schedule']['case_id']==case_id and r['schedule']['method']==m]
               for m in METHODS}
        data={'case_id':case_id,'case':next(iter(group.values()))[0]['case']['label'],
              'methods':{m:{'suite_converged':sum(r['assessment']['suite_converged'] for r in rr),
                    'full_median':median(r['full_seconds'] for r in rr),
                    'solve_median':median(r['proposal']['solve_call_seconds'] for r in rr)} for m,rr in group.items()}}
        paired.append(data)
    out['cases']=paired
    comparisons={}
    for other in METHODS[1:]:
        common=[p for p in paired if p['methods']['seven_block']['suite_converged']==3
                and p['methods'][other]['suite_converged']==3]
        ratios=[p['methods']['seven_block']['full_median']/p['methods'][other]['full_median'] for p in common]
        comparisons[other]=dict(common_cases=len(common),seven_faster=sum(x<1 for x in ratios),
                               median_full_cost_ratio=median(ratios) if ratios else None)
    out['paired_comparison']=comparisons
    out['optimizer_calls']=len(rows)+sum(r['recovery']['repair_lp_calls'] for r in rows)
    out['new_closed_loop_runs']=0
    return out


def audit(data):
    manifest=read(MANIFEST);audit_manifest(manifest)
    assert data['manifest_sha256']==digest(MANIFEST)
    assert len(data['trial_files'])==len(manifest['schedule'])==126
    rows=[]
    for i,entry in enumerate(manifest['schedule']):
        path=RESULTS/'trials'/f'{i:03}.json'
        assert data['trial_files'][path.name]==digest(path)
        row=read(path);audit_trial(manifest['cases'][entry['case_id']],entry,row);rows.append(row)
    assert data['summary']==summary(rows)
    return data['summary']


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',action='store_true');args=p.parse_args()
    if args.run:
        if OUTPUT.exists():raise FileExistsError('Completed archive exists; replay without --run')
        if not MANIFEST.exists():
            import numpy,scipy,sympy,clarabel
            from research.classic_cone_audit.compiler import cp
            started=time.perf_counter();templates();initialization=time.perf_counter()-started
            data=dict(schema=1,sources=sources(),cases=[record(x) for x in suite()],schedule=schedule(),
                versions=dict(numpy=numpy.__version__,scipy=scipy.__version__,sympy=sympy.__version__,
                              cvxpy=cp.__version__,clarabel=clarabel.__version__,python=platform.python_version()),
                platform=platform.platform(),processor=platform.processor(),
                thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},
                template_initialization_seconds=initialization)
            audit_manifest(data);write_new(MANIFEST,data)
        manifest=read(MANIFEST);audit_manifest(manifest);items=suite();templates()
        rows=[]
        for i,entry in enumerate(manifest['schedule']):
            path=RESULTS/'trials'/f'{i:03}.json'
            if path.exists():
                row=read(path);audit_trial(manifest['cases'][entry['case_id']],entry,row)
                print('replayed',i,flush=True)
            else:
                item=items[entry['case_id']]
                print('starting',i,item['case'].label,entry['method'],flush=True)
                row=run_trial(item,entry);audit_trial(manifest['cases'][entry['case_id']],entry,row)
                write_new(path,row)
                print('completed',i,row['assessment'],round(row['full_seconds'],4),flush=True)
            rows.append(row)
        data=dict(schema=1,manifest_sha256=digest(MANIFEST),
                  trial_files={p.name:digest(p) for p in sorted((RESULTS/'trials').glob('*.json'))},
                  summary=summary(rows),novelty_cleared=False,new_closed_loop_runs=0)
        audit(data);write_new(OUTPUT,data)
    data=read(OUTPUT);print(json.dumps(audit(data),indent=2));print('archive_sha256',digest(OUTPUT))


if __name__=='__main__':main()
