"""Run the prespecified support gate once, keeping failures and all witnesses."""
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
import time
import warnings
import numpy as np
from .solver import build, compile_record
from .cases import cases

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'


def hashes():
    paths=sorted(ROOT.glob('*.py'))+[ROOT/'THEORY.md',ROOT/'PROTOCOL.md']
    dependencies=('schur_sparse_lift/algebra.py','schur_sparse_lift/THEORY.md',
                  'schur_sparse_lift/results/verification.json',
                  'schur_remainder_recursion/algebra.py','schur_remainder_recursion/THEORY.md',
                  'moment_four_lift/exact.py','moment_four_lift/THEORY.md',
                  'moment_lift_obstruction/THEORY.md','moment_upper_gate/exact.py')
    paths.extend(ROOT.parent/name for name in dependencies)
    return {str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def exclusive(path,data):
    with path.open('x',encoding='utf-8') as stream:
        json.dump(data,stream,indent=2,allow_nan=False)
        stream.write('\n')


def root_oracle(coef):
    coefficients=np.trim_zeros(np.asarray(coef,dtype=float),'b')
    roots=np.polynomial.polynomial.polyroots(coefficients) if len(coefficients)>1 else []
    cuts=sorted([0.,1.]+[float(z.real) for z in roots if abs(z.imag)<1e-7 and 0<z.real<1])
    point=np.zeros(5)
    for left,right in zip(cuts,cuts[1:]):
        if np.polynomial.polynomial.polyval((left+right)/2,coefficients)>0:
            point+=np.array([(right**(j+1)-left**(j+1))/(j+1) for j in range(5)])
    return float(np.dot(coef,point)),point.tolist()


def main():
    RESULTS.mkdir(exist_ok=True)
    trials=RESULTS/'trials';trials.mkdir(exist_ok=True)
    # Exclusive manifest creation prevents accidental second runs or protocol edits.
    exclusive(RESULTS/'manifest.json',dict(schema=1,sources=hashes(),cases=cases(),
                                         optimization_calls_planned=18,new_closed_loops=0))
    start=time.perf_counter()
    problem,y,a4,alo,ahi,direction,pencils,slacks=build()
    built=time.perf_counter()
    data,_,_=problem.get_problem_data('CLARABEL')
    compiled=time.perf_counter()
    exclusive(RESULTS/'compilation.json',dict(compilation=compile_record(),
                  construction_seconds=built-start,canonicalization_seconds=compiled-built))
    for case in cases():
        raw=np.array([float(F(v)) for v in case['coefficients']])
        norm=float(np.linalg.norm(raw));direction.value=raw/norm
        oracle_start=time.perf_counter()
        oracle_value,oracle_point=root_oracle(raw/norm)
        oracle_elapsed=time.perf_counter()-oracle_start
        failure=None;value=None;started=time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                value=problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,
                                    tol_feas=1e-10,max_iter=300,max_threads=1,warm_start=False)
            except Exception as exc:
                failure=type(exc).__name__+': '+str(exc)
        elapsed=time.perf_counter()-started
        available=(failure is None and value is not None and np.isfinite(value)
                   and all(v.value is not None and np.all(np.isfinite(v.value)) for v in (y,a4,alo,ahi)))
        row=dict(case=case,normalization=norm,exact_normalized_support=float(F(case['support']))/norm,
                 status='solver_error' if failure else problem.status,failure=failure,
                 warnings=[str(w.message) for w in caught],solve_call_seconds=elapsed,
                 iterations=None if problem.solver_stats is None else problem.solver_stats.num_iters,
                 oracle_value=oracle_value,oracle_point=oracle_point,oracle_seconds=oracle_elapsed,
                 optimization_calls=1,new_closed_loops=0,available=bool(available))
        if available:
            row.update(value=float(value),point=y.value.tolist(),aux4=a4.value.tolist(),
                       auxlo=alo.value.tolist(),auxhi=ahi.value.tolist(),
                       minimum_eigenvalue=min(float(np.linalg.eigvalsh(g.value).min()) for _,g in pencils),
                       scalar_slacks=[float(g.value) for g in slacks],
                       signed_error=float(value)-row['exact_normalized_support'])
            row['diagnostic_pass']=bool(abs(row['signed_error'])<=1e-7
                and row['minimum_eigenvalue']>=-1e-7 and min(row['scalar_slacks'])>=-1e-7
                and abs(np.dot(direction.value,y.value)-value)<=1e-9)
        else:
            row['diagnostic_pass']=False
        exclusive(trials/f"{case['index']:02d}.json",row)
        print(case['index'],row['status'],row['diagnostic_pass'],row.get('signed_error'),flush=True)


if __name__=='__main__':
    main()
