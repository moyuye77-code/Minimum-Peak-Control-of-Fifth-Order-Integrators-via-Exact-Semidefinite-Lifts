"""Independent schedule and plant replay; imports no producer or solver."""
from fractions import Fraction as F
from math import factorial,sqrt
from statistics import median
from research.five_synthesis.check import check as check_plan,check_lower,evolve


def extract(plan,origin,start,delta,axes):
    """Midpoint lookup independent of the producer's interval containment lookup."""
    if plan is None:return [F(0),delta],[F(0)]*axes
    mesh=list(map(F,plan['mesh']));levels=list(map(F,plan['controls']));n=len(mesh)-1
    cuts=sorted({F(0),delta}|{origin+t-start for t in mesh if start<origin+t<start+delta})
    controls=[]
    for a in range(axes):
        for left,right in zip(cuts,cuts[1:]):
            time=start+(left+right)/2-origin
            indices=[k for k in range(n) if mesh[k]<=time<mesh[k+1]]
            assert len(indices)<=1
            controls.append(levels[a*n+indices[0]] if indices else F(0))
    return cuts,controls


def replay(episode):
    cfg=episode['scenario'];delta=F(cfg['period']);horizon=F(cfg['horizon'])
    cap=F(cfg['cap']);quantum=F(cfg['quantum']);mode=episode['mode']
    assert (delta,horizon,cap,quantum,cfg['steps'],cfg['order'])==(F(1,2),F(3),F(5),F(1,10**9),12,5)
    assert mode in ('ideal_instant','one_period_deadline') and cfg['exact_state']
    current=list(map(F,cfg['initial']));axes=len(current)//5;assert axes in (1,2)
    assert len(episode['steps'])==12
    active=None;origin=F(0);source=None;busy_until=0
    queried=accepted=certified=converged=late=skips=invalid=wide=0
    no_plan=expired=optimization_calls=warnings=non_success=0
    integral=F(0);peak=F(0);times=[];samples=[current];initial_sufficient=0
    for index,row in enumerate(episode['steps']):
        start=index*delta;now_ns=index*500_000_000
        assert row['index']==index and list(map(F,row['state']))==current
        assert row['active_before']==source and row['busy_until_before_ns']==busy_until
        mesh,intended=extract(active,origin,start,delta,axes)
        command=[max(-cap,min(cap,round(u/quantum)*quantum)) for u in intended]
        expected_prediction=evolve(current,mesh,command,delta) if mode=='one_period_deadline' else current
        req=row['request'];eligible=False;accept=False
        if mode=='one_period_deadline' and now_ns<busy_until:
            assert req is None and row['decision']=='busy_skip';skips+=1
        else:
            assert req is not None;queried+=1
            case=req['problem'];plan=req['recovery']['certificate'];proposal=req['proposal']
            assert case['order']==5 and case['label']=='feedback_five_to_rest'
            assert F(case['horizon'])==horizon and list(map(F,case['initial']))==expected_prediction
            assert len(case['waypoints'])==axes*5
            for i,obs in enumerate(case['waypoints']):
                weights=[F(0)]*(5*axes);weights[i]=F(factorial(5-i%5))/horizon**(5-i%5)
                assert list(map(F,obs['weights']))==weights
                assert F(obs['time'])==horizon and F(obs['value'])==F(obs['error'])==0
            assert proposal['method']==episode['method'] and plan['case']==case
            assert len(req['recovery']['attempts'])==req['recovery']['repair_lp_calls']
            initial=check_lower(case,req['recovery']['initial_dual'])
            if plan.get('certified'):
                lo,upper=check_plan(plan);certified+=1;converged+=bool(plan['converged'])
                initial_sufficient+=upper-initial<=F(plan['tolerance'])*F(plan['scale'])
                eligible=upper<=cap
            else:check_lower(case,plan['dual'])
            elapsed=req['elapsed_ns'];assert type(elapsed) is int and elapsed>=0
            times.append(elapsed/1e9);on_time=elapsed<=500_000_000
            late+=not on_time
            accept=eligible and (mode=='ideal_instant' or on_time)
            wanted='accepted' if accept else ('late_discard' if mode=='one_period_deadline' and not on_time else 'ineligible')
            assert row['decision']==wanted and req['eligible']==eligible and req['accepted']==accept
            accepted+=accept;invalid+=not eligible;wide+=accept and not plan['converged']
            optimization_calls+=proposal['optimization_calls']+req['recovery']['repair_lp_calls']
            warnings+=bool(proposal.get('warnings'))
            if 'success' in proposal:non_success+=not proposal['success']
            else:non_success+=proposal['status'] not in ('optimal','optimal_inaccurate')
            if mode=='one_period_deadline':busy_until=now_ns+elapsed
            elif accept:
                active=plan;origin=start;source=index
                mesh,intended=extract(active,origin,start,delta,axes)
        executed_source=source
        assert row['execution_source']==executed_source
        no_plan+=active is None
        expired+=active is not None and start>=origin+horizon
        act=row['execution'];assert list(map(F,act['mesh']))==mesh
        assert list(map(F,act['intended']))==intended
        n=len(mesh)-1;command=list(map(F,act['command']));actual=list(map(F,act['actual']))
        if cfg['disturbance']=='none':external=[F(0)]*axes
        elif cfg['disturbance']=='pulse':external=[F(1,20),F(-1,25)] if 3<=index<6 else [F(0)]*axes
        else:
            assert cfg['disturbance']=='persistent'
            external=[(1 if (index+a)%4<2 else -1)*F(1,30+10*a) for a in range(axes)]
        assert list(map(F,act['disturbance']))==external
        assert len(command)==len(actual)==len(intended)==axes*n
        for k,(u,v) in enumerate(zip(intended,command)):
            assert v==max(-cap,min(cap,round(u/quantum)*quantum)) and abs(v)<=cap
            assert abs(u-v)<=quantum/2 and actual[k]==v+external[k//n]
        next_state=evolve(current,mesh,actual,delta);nominal=evolve(current,mesh,intended,delta)
        assert list(map(F,act['next_state']))==next_state and list(map(F,act['nominal_state']))==nominal
        for a in range(axes):
            for d in range(5):
                assert abs(next_state[5*a+d]-nominal[5*a+d])<=(abs(external[a])+quantum/2)*delta**(5-d)/factorial(5-d)
                if mode=='one_period_deadline':
                    assert abs(next_state[5*a+d]-expected_prediction[5*a+d])<=abs(external[a])*delta**(5-d)/factorial(5-d)
        local_integral=F(0)
        for k,(left,right) in enumerate(zip(mesh,mesh[1:])):
            at_left=evolve(current,mesh,actual,left)
            for a in range(axes):
                coeff=[at_left[5*a+j]/factorial(j) for j in range(5)]+[actual[a*n+k]/factorial(5)]
                squared=[sum((coeff[i]*coeff[d-i] for i in range(6) if 0<=d-i<6),F(0)) for d in range(11)]
                local_integral+=sum((c*(right-left)**(d+1)/(d+1) for d,c in enumerate(squared)),F(0))
        assert local_integral==F(act['position_squared_integral']) and local_integral>=0
        integral+=local_integral;peak=max([peak]+list(map(abs,command)))
        if mode=='one_period_deadline' and accept:
            active=req['recovery']['certificate'];origin=start+delta;source=index
        assert row['active_after']==source and row['busy_until_after_ns']==busy_until
        current=next_state;samples.append(current)
    assert list(map(F,episode['final_state']))==current
    return dict(slots=12,queries=queried,accepted=accepted,certified=certified,converged=converged,
                accepted_wide=wide,ineligible=invalid,over_period_queries=late,busy_skips=skips,
                no_plan_slots=no_plan,expired_plan_slots=expired,optimization_calls=optimization_calls,
                proposal_warning_queries=warnings,non_success_queries=non_success,
                initial_lower_sufficient=initial_sufficient,input_peak=float(peak),
                continuous_position_rms=sqrt(float(integral/(12*delta*axes))),
                sampled_position_rms=sqrt(float(sum((x[5*a]**2 for x in samples for a in range(axes)),F(0))/(13*axes))),
                final_max_abs_by_derivative=[float(max(abs(current[5*a+d]) for a in range(axes))) for d in range(5)],
                median_planning_seconds=median(times) if times else None,
                maximum_planning_seconds=max(times,default=0))
