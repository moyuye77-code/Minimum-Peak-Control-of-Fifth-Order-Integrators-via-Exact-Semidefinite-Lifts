from fractions import Fraction as F
import sympy as sp


def q_bounds(m, z):
    if not F(0) <= m <= 1 or not m*m/2 <= z <= m-m*m/2:
        raise ValueError('Invalid two-moment prefix')
    if m == 0:
        return F(0), F(0)
    if m == 1:
        return F(1, 3), F(1, 3)
    return z*z/m+m**3/12, F(1, 3)-(F(1, 2)-z)**2/(1-m)-(1-m)**3/12


def lower_fourth(m, z, q):
    low, high = q_bounds(m, z)
    if not low <= q <= high:
        raise ValueError('Invalid three-moment prefix')
    D, Q = z-m*m/2, q-m**3/3
    if D == 0:
        assert Q == 0
        return m**4/4
    return Q*Q/D-m*Q+m*m*D+D*D/2+m**4/4


def fourth_bounds(m, z, q):
    low = lower_fourth(m, z, q)
    high = F(1, 4)-lower_fourth(1-m, F(1, 2)-z, F(1, 3)-q)
    assert low <= high
    return low, high


def lower_endpoints(m, z, q):
    lower_fourth(m, z, q)
    D, Q = z-m*m/2, q-m**3/3
    if not D:
        return sp.S(0), sp.S(0), sp.Rational(m.numerator, m.denominator)
    b = Q/D-m
    delta = (m-b)**2+4*D
    return ((sp.Rational(m+b)-sp.sqrt(sp.Rational(delta)))/2,
            sp.Rational(b),
            (sp.Rational(m+b)+sp.sqrt(sp.Rational(delta)))/2)


def fiber_matrices(m, z, q, r):
    if not 0 < m < 1:
        raise ValueError('Endpoint fibers use point constraints')
    low, high = fourth_bounds(m, z, q)
    if not low <= r <= high:
        raise ValueError('Fourth moment outside fiber')
    matrices, slacks = [], []
    for mass, first, second, third in ((m,z,q,r), (1-m,F(1,2)-z,F(1,3)-q,F(1,4)-r)):
        D, Q = first-mass**2/2, second-mass**3/3
        A, B = (Q*Q/D if D else F(0)), D*D/2
        matrices.extend([[[mass,D],[D,Q-mass*D]], [[D,Q],[Q,A]], [[F(1),D],[D,2*B]]])
        slacks.append(third-(A+B-mass*Q+mass**2*D+mass**4/4))
    return matrices, slacks


def psd_pivots(matrix):
    """Exact Schur-complement PSD test, including zero pivots."""
    work = [[F(x) for x in row] for row in matrix]
    n = len(work)
    assert all(len(row) == n for row in work)
    assert all(work[i][j] == work[j][i] for i in range(n) for j in range(n))
    pivots = []
    for k in range(n):
        pivot = work[k][k]
        assert pivot >= 0, 'Negative exact pivot'
        pivots.append(pivot)
        if pivot == 0:
            assert all(work[k][j] == 0 for j in range(k+1,n)), 'Nonzero row at zero pivot'
        else:
            for i in range(k+1,n):
                for j in range(k+1,n):
                    work[i][j] -= work[i][k]*work[k][j]/pivot
    return pivots


def domination_matrices(y):
    R = (len(y)-1)//2
    assert len(y) == 2*R+1
    complement = [F(1,j+1)-x for j,x in enumerate(y)]
    matrices = []
    for values in (y, complement):
        matrices.append([[values[i+j] for j in range(R+1)] for i in range(R+1)])
        matrices.append([[values[i+j+1]-values[i+j+2] for j in range(R)] for i in range(R)])
    return matrices


def hierarchy_record(R, fixed_direction=False):
    if R < 1:
        raise ValueError('Positive relaxation order required')
    if fixed_direction:
        s = F(1,2)
        base = [s**(j+1)/F(j+1) for j in range(2*R+1)]
        B = sp.Matrix([[sp.Rational(F(1,i+j+1)-base[i+j]) for j in range(R+1)] for i in range(R+1)])
        epsilon = F(1)/F(B.inv()[0,0])
        y = base[:]
        y[0] += epsilon
        true_support = s*s/2
        relaxed_value = s*y[0]-y[1]
    else:
        epsilon = F(1,(R+1)**2)
        s = epsilon
        y = [epsilon]+[F(0)]*(2*R)
        true_support, relaxed_value = s*s/2, s*s
    matrices = domination_matrices(y)
    pivots = [psd_pivots(matrix) for matrix in matrices]
    assert relaxed_value > true_support
    assert y[1] < y[0]*y[0]/2
    return dict(order=R, fixed_direction=fixed_direction, mass_added=str(epsilon),
                y=list(map(str,y)), exact_pivots=[list(map(str,p)) for p in pivots],
                objective_constant=str(s), true_support=str(true_support),
                feasible_relaxation_value=str(relaxed_value),
                gap_lower_bound=str(relaxed_value-true_support))


def prefix_records():
    records = []
    for m in (F(0), F(1,4), F(1,2), F(3,4), F(1)):
        for z in sorted({m*m/2,m/2,m-m*m/2}):
            low, high = q_bounds(m,z)
            for q in sorted({low,(low+high)/2,high}):
                lower, upper = fourth_bounds(m,z,q)
                records.append(dict(m=str(m),z=str(z),q=str(q),lower=str(lower),upper=str(upper)))
    return records
