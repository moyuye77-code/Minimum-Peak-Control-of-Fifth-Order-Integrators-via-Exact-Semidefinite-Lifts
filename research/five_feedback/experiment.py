import argparse
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
import time
from research.five_synthesis.proposals import sdp,continuous_dual
from research.five_synthesis.certify import certify
from research.five_synthesis.check import check as check_plan
from research.five_synthesis.experiment import sources as parent_sources
from .protocol import (scenarios,planning_problem,disturbance,window,predict,execute,
                       CAP,STEPS,PERIOD,PERIOD_NS,MODES)
from .check import replay

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
METHODS={'closed_five_sdp':sdp,'continuous_root_dual':continuous_dual}


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    out=parent_sources()
    paths=sorted(ROOT.glob('*.py'))+[ROOT/'PROTOCOL.md',ROOT/'THEORY.md',ROOT.parent/'five_synthesis/results/verification.json']
    out.update({str(p.relative_to(ROOT.parent)):digest(p) for p in paths})
    return out


def write_new(path,data):
    with path.open('x',encoding='utf-8') as stream:
        json.dump(data,stream,indent=2,allow_nan=False);stream.write('\n')


def schedule():
    out=[];names=tuple(METHODS)
    for i in range(4):
        for j,mode in enumerate(MODES):
            out.extend((i,mode,name) for name in (names if (i+j)%2==0 else names[::-1]))
    return out


def episode(scenario,mode,name,persist=None,clock=time.perf_counter_ns):
    current=scenario.initial;active=None;busy_until=0;rows=[]
    for step in range(STEPS):
        now=step*PERIOD;now_ns=step*PERIOD_NS
        before=None if active is None else active['source'];busy_before=busy_until;req=None;accept=False
        if mode=='one_period_deadline' and now_ns<busy_until:decision='busy_skip'
        else:
            started=clock()
            if mode=='one_period_deadline':
                mesh,intended=window(active,now,scenario.axes);initial=predict(current,mesh,intended)
            else:initial=current
            problem=planning_problem(initial);proposal=METHODS[name](problem);recovery=certify(problem,proposal)
            plan=recovery['certificate'];eligible=False
            if plan.get('certified'):
                _,upper=check_plan(plan);eligible=upper<=CAP
            elapsed=clock()-started;assert elapsed>=0
            accept=eligible and (mode=='ideal_instant' or elapsed<=PERIOD_NS)
            decision='accepted' if accept else ('late_discard' if mode=='one_period_deadline' and elapsed>PERIOD_NS else 'ineligible')
            req=dict(problem=problem.record(),proposal=proposal,recovery=recovery,elapsed_ns=elapsed,
                     eligible=eligible,accepted=accept)
            if mode=='one_period_deadline':busy_until=now_ns+elapsed
            elif accept:active=dict(source=step,start=str(now),certificate=plan)
        mesh,intended=window(active,now,scenario.axes)
        execution_source=None if active is None else active['source']
        # Current external input is first requested after the planner has finished.
        act=execute(current,mesh,intended,disturbance(scenario,step))
        if mode=='one_period_deadline' and accept:
            active=dict(source=step,start=str(now+PERIOD),certificate=req['recovery']['certificate'])
        row=dict(index=step,state=list(map(str,current)),active_before=before,
                 busy_until_before_ns=busy_before,request=req,decision=decision,
                 execution_source=execution_source,execution=act,
                 active_after=None if active is None else active['source'],busy_until_after_ns=busy_until)
        if persist:persist(step,row)
        rows.append(row);current=tuple(map(F,act['next_state']))
        print(scenario.label,mode,name,step,decision,
              None if req is None else req['recovery']['certificate'].get('converged'),flush=True)
    rec=dict(scenario=scenario.record(),mode=mode,method=name,steps=rows,final_state=list(map(str,current)))
    replay(rec)
    return rec


def audit():
    manifest=json.loads((RESULTS/'manifest.json').read_text(encoding='utf-8'))
    assert manifest['sources']==sources() and manifest['scenarios']==[s.record() for s in scenarios()]
    assert manifest['schedule']==[list(v) for v in schedule()]
    summaries=[];hashes={}
    for number,(i,mode,name) in enumerate(schedule()):
        path=RESULTS/'episodes'/f'{number:02d}.json'
        episode_data=json.loads(path.read_text(encoding='utf-8'))
        assert episode_data['scenario']==scenarios()[i].record()
        assert (episode_data['mode'],episode_data['method'])==(mode,name)
        for k,row in enumerate(episode_data['steps']):
            raw=RESULTS/'steps'/f'{number:02d}_{k:02d}.json'
            assert json.loads(raw.read_text(encoding='utf-8'))==row
            hashes[str(raw.relative_to(RESULTS))]=digest(raw)
        hashes[str(path.relative_to(RESULTS))]=digest(path)
        summaries.append(dict(scenario=scenarios()[i].label,mode=mode,method=name,**replay(episode_data)))
    assert len(list((RESULTS/'episodes').glob('*.json')))==16 and len(list((RESULTS/'steps').glob('*.json')))==192
    return dict(schema=1,episodes=summaries,new_closed_loop_runs=16,control_slots=192,
                planning_queries=sum(s['queries'] for s in summaries),
                optimization_calls=sum(s['optimization_calls'] for s in summaries),
                sources=sources(),record_hashes=hashes,manifest_sha256=digest(RESULTS/'manifest.json'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--create',action='store_true');args=parser.parse_args()
    output=RESULTS/'verification.json'
    if args.create:
        for directory in (RESULTS,RESULTS/'episodes',RESULTS/'steps'):directory.mkdir(exist_ok=True)
        import numpy,scipy,sympy,clarabel
        from research.classic_cone_audit.compiler import cp
        write_new(RESULTS/'manifest.json',dict(schema=1,sources=sources(),scenarios=[s.record() for s in scenarios()],
                    schedule=schedule(),versions=dict(numpy=numpy.__version__,scipy=scipy.__version__,
                    sympy=sympy.__version__,cvxpy=cp.__version__,clarabel=clarabel.__version__)))
        for number,(i,mode,name) in enumerate(schedule()):
            save=lambda k,row:write_new(RESULTS/'steps'/f'{number:02d}_{k:02d}.json',row)
            rec=episode(scenarios()[i],mode,name,persist=save)
            write_new(RESULTS/'episodes'/f'{number:02d}.json',rec)
        write_new(output,audit())
    result=audit();assert result==json.loads(output.read_text(encoding='utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k not in ('sources','record_hashes')},indent=2))
    print('archive_sha256',digest(output))


if __name__=='__main__':main()
