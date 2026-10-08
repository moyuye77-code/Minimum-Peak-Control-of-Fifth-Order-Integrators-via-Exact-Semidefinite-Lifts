"""Read-only, standard-library replay of the paper's fifth-order evidence.

Run with `python -B -S research/paper/reproduce.py`. Never use -O: the frozen
independent checkers use assertions. This entry point cannot launch solvers.
"""
from pathlib import Path
from fractions import Fraction as F
from statistics import median
import hashlib
import json
import math
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
RESEARCH = ROOT/'research'
ANCHORS = {
    'five_synthesis': 'e459fc4b0c75923d6a82a379462f9425287ba5589bf2eade9dc1fc3ce315ed36',
    'five_feedback': '5aa4b5aaac5845230ad65f1754967b59ce5ac33d6e28f5b5208d96d918347c3d',
}
METHODS = ('closed_five_sdp','continuous_root_dual')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def inside(root, relative):
    # Manifests were written on Windows; normalize separators, not file contents.
    path = (root/relative.replace('\\','/')).resolve()
    require(path.is_relative_to(root.resolve()), 'path outside the archive')
    return path


def verified_archive(name):
    folder = RESEARCH/name/'results'
    summary = folder/'verification.json'
    require(digest(summary)==ANCHORS[name], name+': verification anchor mismatch')
    saved = read(summary)
    require(digest(folder/'manifest.json')==saved['manifest_sha256'], name+': manifest mismatch')
    manifest = read(folder/'manifest.json')
    require(saved['sources']==manifest['sources'], name+': inconsistent sources')
    for relative,expected in manifest['sources'].items():
        require(digest(inside(RESEARCH,relative))==expected, 'source changed: '+relative)
    key = 'trial_hashes' if name=='five_synthesis' else 'record_hashes'
    base = folder/'trials' if name=='five_synthesis' else folder
    for relative,expected in saved[key].items():
        require(digest(inside(base,relative))==expected, 'record changed: '+relative)
    return folder,manifest,saved


def same_summary(actual, expected):
    """Exact discrete/rational fields; tolerate only final floating display rounding."""
    if isinstance(expected,dict):
        require(isinstance(actual,dict) and actual.keys()==expected.keys(),'summary keys')
        for key in expected:same_summary(actual[key],expected[key])
    elif isinstance(expected,list):
        require(isinstance(actual,list) and len(actual)==len(expected),'summary list')
        for a,b in zip(actual,expected):same_summary(a,b)
    elif isinstance(expected,float):
        require(isinstance(actual,(float,int)) and math.isfinite(actual)
                and math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-14),'floating summary mismatch')
    else:
        require(type(actual) is type(expected) and actual==expected,'exact summary mismatch')


def load_checkers():
    require(__debug__, 'Do not run replay with -O or PYTHONOPTIMIZE.')
    sys.path.insert(0,str(ROOT))
    from research.five_synthesis.check import check,check_lower
    from research.five_feedback.check import replay
    return check,check_lower,replay


def synthesis(folder, manifest, saved, check, check_lower):
    require(len(manifest['cases'])==15 and len(manifest['schedule'])==30,'pilot size')
    require({p.name for p in (folder/'trials').glob('*.json')}==set(saved['trial_hashes']),'pilot record set')
    wanted = [[i,method] for i in range(15) for method in (METHODS if i%2==0 else METHODS[::-1])]
    require(manifest['schedule']==wanted,'pilot method schedule')
    rows=[]
    for number,(index,method) in enumerate(wanted):
        row=read(folder/'trials'/f'{number:02d}.json');spec=manifest['cases'][index]
        require(row['index']==index and row['case']==spec['case'],'pilot input')
        require(row['known_peak']==spec['known_peak'] and row['proposal']['method']==method,'pilot label')
        rec=row['recovery'];cert=rec['certificate']
        require(cert['case']==spec['case'] and len(rec['attempts'])==rec['repair_lp_calls'],'pilot certificate input')
        initial=check_lower(spec['case'],rec['initial_dual'])
        if cert.get('certified'):
            lo,hi=check(cert)
            if spec['known_peak'] is not None:
                known=F(spec['known_peak']);require(lo<=known<=hi and initial<=known,'known anchor excluded')
        else:
            check_lower(spec['case'],cert['dual'])
        rows.append(row)
    result={}
    for method in METHODS:
        selected=[r for r in rows if r['proposal']['method']==method]
        valid=[r for r in selected if r['recovery']['certificate'].get('certified')]
        accurate=[r for r in valid if r['recovery']['certificate']['converged']]
        result[method]=dict(cases=len(selected),certified=len(valid),converged=len(accurate),
            warning_cases=sum(bool(r['proposal']['warnings']) for r in selected),
            recovery_lp_calls=sum(r['recovery']['repair_lp_calls'] for r in selected),
            nonconverged=[r['case']['label'] for r in selected if not r['recovery']['certificate'].get('converged')],
            maximum_scaled_gap=max((float(F(r['recovery']['certificate']['gap'])/F(r['recovery']['certificate']['scale'])) for r in valid),default=None))
    ratios=[]
    for i in range(15):
        pair={r['proposal']['method']:r for r in rows if r['index']==i}
        if all(r['recovery']['certificate'].get('converged') for r in pair.values()):
            ratios.append(pair[METHODS[0]]['elapsed_seconds']/pair[METHODS[1]]['elapsed_seconds'])
    result['paired_timing']=dict(common_converged=len(ratios),sdp_faster=sum(v<1 for v in ratios),
        median_sdp_to_dual_ratio=median(ratios),
        scope='one alternating-order run per case; full proposal+recovery+check; not iid')
    result['optimization_calls']=sum(r['proposal']['optimization_calls']+r['recovery']['repair_lp_calls'] for r in rows)
    result['new_closed_loops']=0
    same_summary(result,saved['summary'])
    return result


def feedback(folder,manifest,saved,replay):
    require(len(manifest['scenarios'])==4,'feedback scenarios')
    wanted=[]
    for i in range(4):
        for j,mode in enumerate(('ideal_instant','one_period_deadline')):
            wanted.extend([i,mode,m] for m in (METHODS if (i+j)%2==0 else METHODS[::-1]))
    require(manifest['schedule']==wanted,'feedback schedule')
    require(len(list((folder/'episodes').glob('*.json')))==16,'episode count')
    require(len(list((folder/'steps').glob('*.json')))==192,'slot count')
    summaries=[]
    for number,(index,mode,method) in enumerate(wanted):
        episode=read(folder/'episodes'/f'{number:02d}.json')
        require(episode['scenario']==manifest['scenarios'][index],'episode input')
        require((episode['mode'],episode['method'])==(mode,method),'episode identity')
        require(len(episode['steps'])==12,'episode slots')
        for step,row in enumerate(episode['steps']):
            require(read(folder/'steps'/f'{number:02d}_{step:02d}.json')==row,'slot copy mismatch')
        values=dict(scenario=episode['scenario']['label'],mode=mode,method=method,**replay(episode))
        same_summary(values,saved['episodes'][number]);summaries.append(values)
    require(sum(v['slots'] for v in summaries)==saved['control_slots']==192,'total slots')
    require(sum(v['queries'] for v in summaries)==saved['planning_queries']==158,'total queries')
    require(sum(v['optimization_calls'] for v in summaries)==saved['optimization_calls']==410,'call count')
    aggregates=[]
    for mode in ('ideal_instant','one_period_deadline'):
        for method in METHODS:
            group=[v for v in summaries if (v['mode'],v['method'])==(mode,method)]
            aggregates.append(dict(mode=mode,method=method,**{key:sum(v[key] for v in group)
                for key in ('queries','certified','converged','accepted','over_period_queries','busy_skips')}))
    return dict(episodes=16,slots=192,queries=158,groups=aggregates)


def main():
    require(__debug__, 'Do not run replay with -O or PYTHONOPTIMIZE.')
    # Validate distributed checkers before importing them.
    pilot=verified_archive('five_synthesis');closed=verified_archive('five_feedback')
    check,check_lower,replay=load_checkers()
    report=dict(synthesis=synthesis(*pilot,check,check_lower),feedback=feedback(*closed,replay))
    forbidden={'numpy','scipy','sympy','cvxpy','clarabel'}
    require(not(forbidden&set(sys.modules)),'unexpected third-party numeric import')
    require(not any(name.endswith(('.proposals','.certify','.experiment','.solver')) for name in sys.modules),'producer imported')
    report.update(new_optimizer_calls=0,new_closed_loop_runs=0,
                  stored_time_replay_not_new_timing=True,anchors=ANCHORS)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
