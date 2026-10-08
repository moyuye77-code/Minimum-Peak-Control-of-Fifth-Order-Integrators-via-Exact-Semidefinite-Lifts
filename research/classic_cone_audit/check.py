"""Closed-form checks of compiler records, independent of its construction."""
from fractions import Fraction as F


def check_primitives(data):
    assert data['square'] == dict(variable_order=['x', 't'], exponent='2',
        soc_affine_rows=[[[1, 0, 1], [-1, 0, 1], [0, 2, 0]]])
    assert data['cube'] == dict(variable_order=['x', 't', 'v'], exponent='3',
        approximation_error='0.0', weights=['1/3', '2/3'],
        soc_affine_rows=[[[1, 0, 0, 1], [1, 0, 0, -1], [0, 2, 0, 0]],
                         [[0, 1, 1, 0], [0, -1, 1, 0], [0, 0, 0, 2]]])
    assert data['ratio'] == dict(variable_order=['m', 'z', 't'],
        soc_affine_rows=[[[0, 1, 0, 1], [0, 1, 0, -1], [0, 0, 2, 0]]])


def check_record(data):
    n = data['segments']
    assert n >= 1 and data['dcp'] and data['new_optimizer_calls'] == 0
    assert list(map(F, data['observation_knots'])) == [F(k*k, n*n) for k in range(n+1)]
    assert data['variables'] == 14*n+3
    assert data['equalities'] == 3*n and data['nonnegative_rows'] == 8*n+6
    assert data['soc_dimensions'] == [3]*(8*n)
    assert data['psd_dimensions'] == data['power_exponents'] == []
    assert data['exponential_cones'] == 0
    assert data['canonicalizers'] == ['FlipObjective', 'Dcp2Cone', 'CvxAttr2Constr',
                                      'ConeMatrixStuffing', 'CLARABEL']


def exact_point_check(m, z, q):
    """Witness lift sanity check; universal equivalence is in THEORY.md."""
    assert 0 <= m <= 1 and m*m/2 <= z <= m-m*m/2
    s0 = z*z/m if m else F(0)
    s1 = (F(1, 2)-z)**2/(1-m) if m != 1 else F(0)
    assert s0+m**3/12 <= q <= F(1, 3)-s1-(1-m)**3/12
    # Compiler's two square epigraph variables are eliminated in the manual lift.
    r0 = r1 = m*m
    assert r0 <= 2*z and r1 <= 2*(m-z)
    cones = [(r0, F(1), m), (r1, F(1), m),
             (m, s0, z), (m*m, F(1), m), (m, m**3, m*m),
             (1-m, s1, F(1, 2)-z), ((1-m)**2, F(1), 1-m),
             (1-m, (1-m)**3, (1-m)**2)]
    assert all(x >= 0 and y >= 0 and x*y-w*w >= 0 for x, y, w in cones)
    return True
