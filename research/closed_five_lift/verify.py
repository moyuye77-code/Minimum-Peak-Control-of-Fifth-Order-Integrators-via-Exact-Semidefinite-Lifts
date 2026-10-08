"""Read-only source, exact witness, and saved floating-point residual audit."""
import argparse
import hashlib
import json
from fractions import Fraction as F
import numpy as np
import sympy as s
from .algebra import model,witness,evaluate,complement,boundary_cases
from .solver import pencil
from .cases import cases
from .experiment import ROOT,RESULTS,hashes,exclusive,root_oracle
from research.moment_four_lift.exact import matrices,atomic_aux
from research.moment_upper_gate.exact import psd_pivots


def exact_record(point):
    result=[]
    for side in (point,complement(point)):
        values,valuation=witness(tuple(side[:4]))
        mats,outputs=evaluate(values)
        assert outputs[0]==1 and outputs[1:5]==tuple(side[:4]) and outputs[5]<=side[4]
        pivots={name:list(map(str,psd_pivots(mat.tolist()))) for name,mat in mats.items()}
        result.append(dict(valuation=valuation,values=list(map(str,values)),
                           lower_slack=str(side[4]-outputs[5]),pivots=pivots))
    prefix=tuple(map(F,point[:4]))
    for mat in matrices(prefix,atomic_aux(prefix[0],prefix[1])):psd_pivots(mat)
    return result


def verify_numeric(row):
    case=row['case'];raw=np.array([float(F(v)) for v in case['coefficients']])
    norm=np.linalg.norm(raw);target=float(F(case['support']))/norm
    assert abs(row['normalization']-norm)<=1e-15 and abs(row['exact_normalized_support']-target)<=1e-15
    oracle,_=root_oracle(raw/norm)
    assert abs(oracle-row['oracle_value'])<=1e-13 and abs(oracle-target)<=1e-10
    if not row['available']:
        assert row['diagnostic_pass'] is False
        return
    y=np.array(row['point']);a4=np.array(row['aux4']);lo=np.array(row['auxlo']);hi=np.array(row['auxhi'])
    mats=list(matrices(y[:4],a4))
    mats.extend(mat for _,mat in pencil([1,*y[:4],*lo]))
    mats.extend(mat for _,mat in pencil([1,*(1/(i+1)-y[i] for i in range(4)),*hi]))
    eigen=min(float(np.linalg.eigvalsh(np.asarray(mat,dtype=float)).min()) for mat in mats)
    slacks=[y[4]-lo[0],.2-y[4]-hi[0]]
    error=row['value']-target
    assert abs(eigen-row['minimum_eigenvalue'])<=1e-12
    assert np.max(np.abs(np.array(slacks)-row['scalar_slacks']))<=1e-14
    assert abs(error-row['signed_error'])<=1e-14
    accepted=bool(abs(error)<=1e-7 and eigen>=-1e-7 and min(slacks)>=-1e-7
                  and abs(np.dot(raw/norm,y)-row['value'])<=1e-9)
    assert accepted==row['diagnostic_pass']


def build_record():
    manifest=json.loads((RESULTS/'manifest.json').read_text(encoding='utf-8'))
    assert manifest['sources']==hashes() and manifest['cases']==cases()
    rows=[]
    for case in cases():
        row=json.loads((RESULTS/'trials'/f"{case['index']:02d}.json").read_text(encoding='utf-8'))
        assert row['case']==case
        verify_numeric(row);rows.append(row)
    exact=[dict(name=name,point=list(map(str,point)),sides=exact_record(point)) for name,point in boundary_cases()]
    for case in cases():
        exact.append(dict(name=f"support_{case['index']}",point=case['point'],
                          sides=exact_record(tuple(s.Rational(v) for v in case['point']))))
    return dict(schema=1,exact=exact,functional_dimension=len(model()['basis']),
                numerical_cases=len(rows),diagnostic_passes=sum(r['diagnostic_pass'] for r in rows),
                warning_cases=sum(bool(r['warnings']) for r in rows),
                failure_cases=[r['case']['index'] for r in rows if not r['diagnostic_pass']],
                maximum_absolute_support_error=max((abs(r['signed_error']) for r in rows if r['available']),default=None),
                optimization_calls=sum(r['optimization_calls'] for r in rows),new_closed_loops=0,
                sources=hashes(),trial_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                               for p in sorted((RESULTS/'trials').glob('*.json'))})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    record=build_record();path=RESULTS/'verification.json'
    if args.write:exclusive(path,record)
    assert json.loads(path.read_text(encoding='utf-8'))==record
    print('closed K5 exact witnesses and all saved numerical residuals replayed; no optimizer called')
    print({k:record[k] for k in ('functional_dimension','numerical_cases','diagnostic_passes','warning_cases',
                                'failure_cases','maximum_absolute_support_error','optimization_calls','new_closed_loops')})
    print('archive_sha256',hashlib.sha256(path.read_bytes()).hexdigest())


if __name__=='__main__':main()
