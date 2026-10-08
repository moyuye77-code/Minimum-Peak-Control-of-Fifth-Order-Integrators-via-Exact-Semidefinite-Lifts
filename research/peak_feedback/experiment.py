import argparse
import hashlib
import json
from pathlib import Path
import time
from fractions import Fraction as F
from research.peak_synthesis.proposals import sdp,continuous_dual
from research.peak_synthesis.certify import certify
from research.peak_synthesis.check import check as check_plan
from research.peak_synthesis.experiment import sources as parent_sources
from .protocol import scenarios,planning_problem,disturbance,execute,INPUT_CAP,STEPS
from .check import replay

ROOT=Path(__file__).resolve().parent
OUTPUT=ROOT/'results'/'verification.json'
METHODS=(sdp,continuous_dual)


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    out=parent_sources()
    out.update({str(p.relative_to(ROOT.parent)):digest(p) for p in sorted(ROOT.glob('*.py'))})
    return out


def episode(scenario,method):
    current=scenario.initial; rows=[]
    for step in range(STEPS):
        problem=planning_problem(current); started=time.perf_counter()
        proposed=method(problem); recovery=certify(problem,proposed)
        certificate=recovery['certificate']; accepted=False
        if certificate.get('certified'):
            _,upper=check_plan(certificate); accepted=upper<=INPUT_CAP
        elapsed=time.perf_counter()-started
        # The disturbance is obtained only after planning, never passed to the planner.
        act=execute(current,certificate if accepted else None,disturbance(scenario,step))
        rows.append(dict(index=step,state=list(map(str,current)),problem=problem.record(),
                         proposal=proposed,recovery=recovery,accepted=accepted,
                         planning_seconds=elapsed,execution=act))
        current=tuple(map(F,act['next_state']))
        print(scenario.label,method.__name__,step,proposed['status'],accepted,
              certificate.get('converged'),float(max(abs(current[4*a]) for a in range(scenario.axes))),flush=True)
    rec=dict(scenario=scenario.record(),method=rows[0]['proposal']['method'],
             steps=rows,final_state=list(map(str,current)))
    replay(rec)
    return rec


def audit(data):
    assert data['sources']==sources() and len(data['episodes'])==8
    summaries=[]
    for index,scenario in enumerate(scenarios()):
        for j,name in enumerate(('homogeneous_sdp','continuous_root_dual')):
            rec=data['episodes'][2*index+j]
            assert rec['scenario']==scenario.record() and rec['method']==name
            assert all(row['proposal']['method']==name for row in rec['steps'])
            summaries.append(dict(scenario=scenario.label,method=name,**replay(rec)))
    return dict(episodes=summaries,new_closed_loop_runs=8,
                planning_queries=sum(s['steps'] for s in summaries),
                archived_optimization_calls=sum(s['optimization_calls'] for s in summaries),
                scope='four deterministic scenarios times two methods; not independent stochastic trials')


def main():
    p=argparse.ArgumentParser(); p.add_argument('--create',action='store_true'); args=p.parse_args()
    if args.create:
        if OUTPUT.exists(): raise FileExistsError('Preserve feedback archive')
        episodes=[episode(s,m) for s in scenarios() for m in METHODS]
        import numpy,scipy,sympy,clarabel
        from research.classic_cone_audit.compiler import cp
        data=dict(schema=1,sources=sources(),episodes=episodes,
                  versions=dict(numpy=numpy.__version__,scipy=scipy.__version__,sympy=sympy.__version__,
                                cvxpy=cp.__version__,clarabel=clarabel.__version__))
        audit(data); OUTPUT.parent.mkdir(parents=True,exist_ok=True)
        with OUTPUT.open('x',encoding='utf-8') as stream:
            json.dump(data,stream,indent=2); stream.write('\n')
    data=json.loads(OUTPUT.read_text(encoding='utf-8'))
    print(json.dumps(audit(data),indent=2)); print('archive_sha256',digest(OUTPUT))


if __name__=='__main__': main()
