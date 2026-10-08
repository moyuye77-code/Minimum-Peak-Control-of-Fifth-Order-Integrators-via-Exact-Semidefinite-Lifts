"""Independent feedback replay, with rational plant dynamics from the old checker."""
from fractions import Fraction as F
from math import factorial,sqrt
from statistics import median
from research.peak_synthesis.check import check as check_plan,evolve


def replay(episode):
    cfg=episode['scenario']; h=F(cfg['horizon']); delta=F(cfg['period'])
    cap=F(cfg['input_cap']); quantum=F(cfg['actuator_quantum'])
    assert (h,delta,cap,quantum,cfg['steps'])==(F(2),F(1,4),F(5),F(1,10**9),20)
    current=list(map(F,cfg['initial'])); axes=len(current)//4
    assert axes in (1,2) and cfg['exact_state_measurement']
    assert len(episode['steps'])==cfg['steps']
    samples=[current]; input_peak=F(0); accepted=converged=fallbacks=0; calls=0
    for index,row in enumerate(episode['steps']):
        assert row['index']==index and list(map(F,row['state']))==current
        plan=row.get('recovery',{}).get('certificate',{})
        case=row['problem']; assert list(map(F,case['initial']))==current
        assert F(case['horizon'])==h and len(case['waypoints'])==4*axes
        for i,o in enumerate(case['waypoints']):
            expected=[F(0)]*(4*axes); expected[i]=F(factorial(4-i%4))/h**(4-i%4)
            assert list(map(F,o['weights']))==expected
            assert F(o['time'])==h and F(o['value'])==0 and F(o['error'])==0
        eligible=False
        if plan.get('certified'):
            assert plan['case']==case
            low,upper=check_plan(plan); eligible=upper<=cap
            converged+=bool(plan['converged'])
        assert row['accepted']==eligible
        accepted+=eligible; fallbacks+=not eligible
        calls+=row.get('proposal',{}).get('optimization_calls',0)
        calls+=row.get('recovery',{}).get('repair_lp_calls',0)
        act=row['execution']; mesh=list(map(F,act['mesh']))
        assert mesh[0]==0 and mesh[-1]==delta and all(l<r for l,r in zip(mesh,mesh[1:]))
        m=len(mesh)-1
        intended=list(map(F,act['intended'])); commanded=list(map(F,act['commanded']))
        actual=list(map(F,act['actual'])); external=list(map(F,act['disturbance']))
        assert len(intended)==len(commanded)==len(actual)==axes*m and len(external)==axes
        if eligible:
            fullmesh=list(map(F,plan['mesh'])); controls=list(map(F,plan['controls']))
            fullcells=len(fullmesh)-1
            assert mesh==[t for t in fullmesh if t<delta]+[delta]
            assert intended==[controls[a*fullcells+k] for a in range(axes) for k in range(m)]
        else:
            assert mesh==[F(0),delta] and intended==[F(0)]*axes
        for k,(u,v) in enumerate(zip(intended,commanded)):
            assert v==max(-cap,min(cap,round(u/quantum)*quantum))
            assert abs(u-v)<=quantum/2 and abs(v)<=cap
            assert actual[k]==v+external[k//m]
        if cfg['disturbance']=='none': expected_external=[F(0)]*axes
        elif cfg['disturbance']=='pulse':
            expected_external=[F(1,20),F(-1,25)] if 4<=index<8 else [F(0)]*axes
        else:
            assert cfg['disturbance']=='persistent'
            expected_external=[(1 if (index+2*a)%8<4 else -1)*F(1,30+10*a) for a in range(axes)]
        assert external==expected_external
        planned=evolve(current,mesh,intended,delta)
        next_state=evolve(current,mesh,actual,delta)
        assert planned==list(map(F,act['planned_prefix_state']))
        assert next_state==list(map(F,act['next_state']))
        for a in range(axes):
            for d in range(4):
                local_bound=(abs(external[a])+quantum/2)*delta**(4-d)/factorial(4-d)
                assert abs(next_state[4*a+d]-planned[4*a+d])<=local_bound
        input_peak=max([input_peak]+list(map(abs,commanded)))
        current=next_state; samples.append(current)
    assert list(map(F,episode['final_state']))==current
    # RMS is of sampled positions, not a continuous-time integral metric.
    rms=sqrt(float(sum((x[4*a]**2 for x in samples for a in range(axes)),F(0))/(len(samples)*axes)))
    return dict(steps=len(episode['steps']),accepted=accepted,converged=converged,fallbacks=fallbacks,
        optimization_calls=calls,input_peak=float(input_peak),position_rms_samples=rms,
        final_max_abs_by_derivative=[float(max(abs(current[4*a+d]) for a in range(axes))) for d in range(4)],
        median_planning_seconds=median(row['planning_seconds'] for row in episode['steps']))
