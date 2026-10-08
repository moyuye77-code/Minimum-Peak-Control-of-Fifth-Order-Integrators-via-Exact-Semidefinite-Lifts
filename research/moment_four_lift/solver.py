"""Numerical audit only: floating-point support comparisons are not certificates."""
import time
import warnings
from fractions import Fraction as F
import numpy as np
from research.classic_cone_audit.compiler import cp
from .exact import matrices, directions


def build(compressed=True):
    y = cp.Variable(4, name='four_moments')
    aux = cp.Variable(6, name='coordinate_functional')
    objective = cp.Parameter(4, name='support_direction')
    pencils = [cp.bmat(mat) for mat in matrices(y, aux, compressed)]
    constraints = [mat >> 0 for mat in pencils]
    problem = cp.Problem(cp.Maximize(objective@y), constraints)
    assert problem.is_dcp()
    return problem, y, aux, objective, pencils


def compile_record(compressed=True):
    problem, *_ = build(compressed)
    data, _, _ = problem.get_problem_data('CLARABEL')
    dims = data['dims']
    return dict(compressed=compressed, dcp=problem.is_dcp(),
                variables=int(data['A'].shape[1]),
                psd_dimensions=list(map(int,dims.psd)),
                soc_dimensions=list(map(int,dims.soc)),
                equalities=int(dims.zero), nonnegative_rows=int(dims.nonneg))


def run_supports(compressed=True):
    problem, y, aux, objective, pencils = build(compressed)
    results = []
    for row in directions():
        raw = np.array([float(F(x)) for x in row['direction']])
        scale = np.linalg.norm(raw)
        objective.value = raw/scale
        started = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            failure = None
            try:
                value = problem.solve(solver='CLARABEL', warm_start=False,
                                      tol_gap_abs=1e-10, tol_feas=1e-10,
                                      tol_gap_rel=1e-10, max_iter=300)
            except cp.error.SolverError as exc:
                value, failure = None, str(exc)
        target = float(F(row['support']))/scale
        available = failure is None and value is not None and np.isfinite(value) and y.value is not None
        eig = [float(np.linalg.eigvalsh(mat.value).min()) for mat in pencils] if available else []
        results.append(dict(direction=row['direction'], exact_support=row['support'],
                            normalization=float(scale), true_normalized_support=target,
                            status='solver_error' if failure else problem.status,
                            solver_error=failure, value=float(value) if available else None,
                            signed_error=float(value-target) if available else None,
                            point=y.value.tolist() if available else None,
                            auxiliary=aux.value.tolist() if available else None,
                            minimum_eigenvalues=eig,
                            elapsed_seconds=time.perf_counter()-started,
                            warnings=[str(w.message) for w in caught]))
    return dict(compilation=compile_record(compressed), cases=results,
                solver='CLARABEL', numeric_only=True, closed_loop_runs=0)
