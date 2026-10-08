"""Executable affine closed-K5 pencil; floating point is not a proof."""
from functools import lru_cache
import numpy as np
from research.classic_cone_audit.compiler import cp
from research.moment_four_lift.exact import matrices as four_matrices
from .algebra import model


@lru_cache(maxsize=1)
def templates():
    data = model()
    out = []
    for name,matrix in data['matrices'].items():
        array = np.zeros((matrix.rows,matrix.cols,len(data['theta'])))
        for i in range(matrix.rows):
            for j in range(matrix.cols):
                for k,t in enumerate(data['theta']):
                    array[i,j,k] = float(matrix[i,j].diff(t))
        out.append((name,array))
    return tuple(out)


def pencil(values):
    return [(name, [[sum(float(c)*values[k] for k,c in enumerate(row) if c)
                    for row in matrix_row] for matrix_row in array])
            for name,array in templates()]


def constraints(y,aux4,auxlo,auxhi):
    raw = [(f'prefix_{i}',g) for i,g in enumerate(four_matrices(list(y[:4]),list(aux4)))]
    lo = [1,*list(y[:4]),*list(auxlo)]
    hi = [1,*(1/(i+1)-y[i] for i in range(4)),*list(auxhi)]
    raw.extend(('lower_'+name,g) for name,g in pencil(lo))
    raw.extend(('upper_'+name,g) for name,g in pencil(hi))
    pencils = [(name,cp.bmat(g)) for name,g in raw]
    slacks = [y[4]-auxlo[0], 0.2-y[4]-auxhi[0]]
    return [g>>0 for _,g in pencils]+[g>=0 for g in slacks], pencils, slacks


def build():
    y = cp.Variable(5,name='five_moments')
    aux4 = cp.Variable(6,name='prefix_aux')
    auxlo = cp.Variable(14,name='lower_functional_aux')
    auxhi = cp.Variable(14,name='upper_functional_aux')
    direction = cp.Parameter(5,name='support_direction')
    cons,pencils,slacks = constraints(y,aux4,auxlo,auxhi)
    problem = cp.Problem(cp.Maximize(direction@y),cons)
    assert problem.is_dcp()
    return problem,y,aux4,auxlo,auxhi,direction,pencils,slacks


def compile_record():
    problem,*_ = build()
    data,_,_ = problem.get_problem_data('CLARABEL')
    dims = data['dims']
    return dict(variables=int(data['A'].shape[1]),equalities=int(dims.zero),
                nonnegative_rows=int(dims.nonneg),PSD_blocks=list(map(int,dims.psd)),
                SOC_blocks=list(map(int,dims.soc)),dcp=problem.is_dcp())
