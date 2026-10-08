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
RESULTS=ROOT/'results'
METHODS={'closed_five_sdp':sdp,'continuous_root_dual':continuous_dual}


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    paths=sorted(ROOT.glob('*.py'))+[ROOT/'THEORY.md',ROOT/'PROTOCOL.md']
    dependencies=('closed_five_lift/algebra.py','closed_five_lift/THEORY.md','closed_five_lift/results/verification.json',
                  'schur_sparse_lift/algebra.py','schur_remainder_recursion/algebra.py',
                  'general_congruence_gate/algebra.py','moment_four_lift/exact.py','moment_upper_gate/exact.py',
                  'peak_synthesis/model.py','peak_synthesis/certify.py','four_history_audit/proposals.py',
                  'four_history_audit/model.py','classic_cone_audit/compiler.py')
    paths.extend(ROOT.parent/name for name in dependencies)
    return {str(path.relative_to(ROOT.parent)):digest(path) for path in paths}


def write_new(path,data):
    with path.open('x',encoding='utf-8') as stream:
        json.dump(data,stream,indent=2,allow_nan=False);stream.write('\n')


def schedule():
    names=tuple(METHODS)
    return [(index,name) for index in range(len(cases())) for name in (names if index%2==0 else names[::-1])]


def summary(rows):
    out={}
    for name in METHODS:
        chosen=[r for r in rows if r['proposal']['method']==name]
        good=[r for r in chosen if r['recovery']['certificate'].get('certified')]
        converged=[r for r in good if r['recovery']['certificate']['converged']]
        out[name]=dict(cases=len(chosen),certified=len(good),converged=len(converged),
                       warning_cases=sum(bool(r['proposal']['warnings']) for r in chosen),
                       recovery_lp_calls=sum(r['recovery']['repair_lp_calls'] for r in chosen),
                       nonconverged=[r['case']['label'] for r in chosen if not r['recovery']['certificate'].get('converged')],
                       maximum_scaled_gap=max((float(F(r['recovery']['certificate']['gap'])/F(r['recovery']['certificate']['scale']))
                                               for r in good),default=None))
    pairs=[]
    for index in range(len(cases())):
        pair={r['proposal']['method']:r for r in rows if r['index']==index}
        if all(r['recovery']['certificate'].get('converged') for r in pair.values()):
            a=pair['closed_five_sdp']['elapsed_seconds'];b=pair['continuous_root_dual']['elapsed_seconds'];pairs.append((a,b))
    out['paired_timing']=dict(common_converged=len(pairs),sdp_faster=sum(a<b for a,b in pairs),
                             median_sdp_to_dual_ratio=median(a/b for a,b in pairs) if pairs else None,
                             scope='one alternating-order run per case; full proposal+recovery+check; not iid')
    out['optimization_calls']=sum(r['proposal']['optimization_calls']+r['recovery']['repair_lp_calls'] for r in rows)
    out['new_closed_loops']=0
    return out


def audit():
    manifest=json.loads((RESULTS/'manifest.json').read_text(encoding='utf-8'))
    assert manifest['sources']==sources()
    expected=[dict(case=c.record(),known_peak=None if k is None else str(k)) for c,k in cases()]
    assert manifest['cases']==expected and manifest['schedule']==[list(x) for x in schedule()]
    paths=sorted((RESULTS/'trials').glob('*.json'));assert len(paths)==30
    rows=[]
    for number,(index,name) in enumerate(schedule()):
        path=RESULTS/'trials'/f'{number:02d}.json';row=json.loads(path.read_text(encoding='utf-8'))
        case,known=cases()[index]
        assert row['index']==index and row['case']==case.record() and row['proposal']['method']==name
        assert row['known_peak']==(None if known is None else str(known))
        rec=row['recovery'];assert rec['repair_lp_calls']==len(rec['attempts'])
        initial=check_lower(row['case'],rec['initial_dual'])
        certificate=rec['certificate'];assert certificate['case']==case.record()
        if certificate.get('certified'):
            lower,upper=check(certificate)
            if known is not None:assert lower<=known<=upper and initial<=known
        else:check_lower(row['case'],certificate['dual'])
        rows.append(row)
    return dict(schema=1,summary=summary(rows),sources=sources(),
                trial_hashes={path.name:digest(path) for path in paths},
                manifest_sha256=digest(RESULTS/'manifest.json'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--create',action='store_true');args=parser.parse_args()
    output=RESULTS/'verification.json'
    if args.create:
        RESULTS.mkdir(exist_ok=True);(RESULTS/'trials').mkdir(exist_ok=True)
        import numpy,scipy,sympy,clarabel
        from research.classic_cone_audit.compiler import cp
        write_new(RESULTS/'manifest.json',dict(schema=1,sources=sources(),
                    cases=[dict(case=c.record(),known_peak=None if k is None else str(k)) for c,k in cases()],
                    schedule=schedule(),versions=dict(numpy=numpy.__version__,scipy=scipy.__version__,
                    sympy=sympy.__version__,cvxpy=cp.__version__,clarabel=clarabel.__version__)))
        for number,(index,name) in enumerate(schedule()):
            case,known=cases()[index];started=time.perf_counter()
            proposal=METHODS[name](case);recovery=certify(case,proposal);certificate=recovery['certificate']
            if certificate.get('certified'):check(certificate)
            else:check_lower(case.record(),certificate['dual'])
            row=dict(index=index,case=case.record(),known_peak=None if known is None else str(known),
                     proposal=proposal,recovery=recovery,elapsed_seconds=time.perf_counter()-started)
            write_new(RESULTS/'trials'/f'{number:02d}.json',row)
            print(number,case.label,name,proposal['status'],certificate.get('certified'),
                  certificate.get('converged'),recovery['repair_lp_calls'],flush=True)
        write_new(output,audit())
    verified=audit();assert json.loads(output.read_text(encoding='utf-8'))==verified
    print(json.dumps(verified['summary'],indent=2));print('archive_sha256',digest(output))


if __name__=='__main__':main()
