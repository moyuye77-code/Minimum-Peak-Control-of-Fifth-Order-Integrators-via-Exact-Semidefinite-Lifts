"""Classic moment inequalities only; no hand-written cones and no solve call."""
from fractions import Fraction as F
from pathlib import Path
import sys
import numpy as np

RESEARCH = Path(__file__).resolve().parents[1]
for folder in (RESEARCH/'.deps'/'sls_baseline', RESEARCH/'.deps'/'moment_soc'):
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))
import cvxpy as cp
from cvxpy.reductions.dcp2cone.canonicalizers.power_canon import power_canon
from cvxpy.reductions.dcp2cone.canonicalizers.quad_over_lin_canon import quad_over_lin_canon


def classic_constraints(m, z, q):
    return [m >= 0, m <= 1,
            cp.square(m) <= 2*z, cp.square(m) <= 2*(m-z),
            cp.quad_over_lin(z, m)+cp.power(m, 3)/12 <= q,
            cp.quad_over_lin(0.5-z, 1-m)+cp.power(1-m, 3)/12 <= 1/3-q]


def build(segments):
    assert isinstance(segments, int) and segments >= 1
    # Nonuniform, declared observation knots. No mesh between these knots.
    times = [F(k*k, segments*segments) for k in range(segments+1)]
    states = cp.Variable((segments+1, 3), name='node_states')
    moments = cp.Variable((segments, 3), name='segment_moments')
    constraints = [states[0] <= 0, states[0] >= 0]
    for k in range(segments):
        m, z, q = (moments[k, i] for i in range(3))
        p, v, a = (states[k, i] for i in range(3))
        h = float(times[k+1]-times[k])
        constraints += classic_constraints(m, z, q)
        terminal = cp.hstack([p+h*v+h*h*a/2+h**3*(m-2*z+q-1/6),
                              v+h*a+h*h*(2*(m-z)-0.5), a+h*(2*m-1)])
        constraints.append(states[k+1] == terminal)
        observed = states[k+1, 0]+((-1)**k/3)*states[k+1, 1]+((k+1)/(segments+1))*states[k+1, 2]
        constraints += [observed <= 1, observed >= -1]
    objective = states[-1, 0]-states[-1, 1]/3+states[-1, 2]/2
    return cp.Problem(cp.Maximize(objective), constraints), times


def compile_record(segments):
    problem, times = build(segments)
    assert problem.is_dcp()
    data, chain, _ = problem.get_problem_data('CLARABEL')
    dims = data['dims']
    return dict(segments=segments, observation_knots=list(map(str, times)),
                dcp=True, variables=int(data['A'].shape[1]),
                equalities=int(dims.zero), nonnegative_rows=int(dims.nonneg),
                soc_dimensions=list(map(int, dims.soc)), psd_dimensions=list(map(int, dims.psd)),
                exponential_cones=int(dims.exp), power_exponents=list(map(str, dims.p3d)),
                canonicalizers=[type(item).__name__ for item in chain.reductions],
                new_optimizer_calls=0)


def affine_soc_rows(cones, variables):
    """Read actual primitive canonicalizer output on an affine basis.

    All tested primitive coefficients are integers exactly representable as
    binary floats. This is NOT rational recovery of an arbitrary numeric model.
    """
    def evaluate():
        values = []
        for cone in cones:
            assert isinstance(cone, cp.constraints.second_order.SOC)
            row = np.concatenate([np.asarray(arg.value).ravel() for arg in cone.args])
            assert len(row) == 3 and all(float(x).is_integer() for x in row)
            values.append([int(x) for x in row])
        return values

    for var in variables:
        var.value = np.zeros(var.shape)
    constant = evaluate()
    columns = []
    for var in variables:
        var.value = np.ones(var.shape)
        value = evaluate()
        columns.append([[b-a for a, b in zip(c, v)] for c, v in zip(constant, value)])
        var.value = np.zeros(var.shape)
    return [[[constant[k][j]]+[col[k][j] for col in columns] for j in range(3)]
            for k in range(len(cones))]


def primitive_records():
    x = cp.Variable(name='x')
    square = cp.power(x, 2)
    square_t, square_cones = power_canon(square, [x])
    cube = cp.power(x, 3)
    cube_t, cube_cones = power_canon(cube, [x])
    extras = {v.id: v for c in cube_cones for v in c.variables() if v.id not in (x.id, cube_t.id)}
    assert len(extras) == 1 and cube.p_rational == 3 and cube.approx_error == 0
    m, z = cp.Variable(name='m'), cp.Variable(name='z')
    ratio_t, ratio_cones = quad_over_lin_canon(cp.quad_over_lin(z, m), [z, m])
    return dict(
        square=dict(variable_order=['x', 't'], exponent=str(square.p_rational),
                    soc_affine_rows=affine_soc_rows(square_cones, [x, square_t])),
        cube=dict(variable_order=['x', 't', 'v'], exponent=str(cube.p_rational),
                  approximation_error=str(cube.approx_error),
                  weights=list(map(str, cube.w)),
                  soc_affine_rows=affine_soc_rows(cube_cones, [x, cube_t, next(iter(extras.values()))])),
        ratio=dict(variable_order=['m', 'z', 't'],
                   soc_affine_rows=affine_soc_rows(ratio_cones, [m, z, ratio_t])))
