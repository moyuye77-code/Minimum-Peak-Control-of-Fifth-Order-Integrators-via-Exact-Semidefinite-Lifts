"""Exact identities only. No numerical SOS fits or optimization."""
from functools import lru_cache
import sympy as s
from research.general_congruence_gate.algebra import (
    fifth_schur, exponential_coefficients, boundary_prediction, m, z, q, r,
)

X = (m, z, q, r)
M, Z, Q, R, c = s.symbols('M Z Q R c')
U = (M, Z, Q, R)
AT = dict(zip(X, U))
DELTA = m*q-z*z-m**4/12
DELTA0 = DELTA.xreplace(AT)
N = m*r-z*q-m**3*z/6
N0 = N.xreplace(AT)
P = m*DELTA
P0 = M*DELTA0


def remainder(f, variables, reference):
    at = dict(zip(variables, reference))
    return f-f.xreplace(at)-sum(s.diff(f, t).xreplace(at)*(t-u)
                               for t, u in zip(variables, reference))


def h0():
    return q*q/m+m*m*q/12+m*z*z/4-m**5/720


def h0_certificate():
    h = m-M
    return ((q-m*Q/M+m*M*h/12)**2/m
            +(M*z+(m-2*M)*Z)**2/(6*M)
            +((2*m-M)*z-m*Z)**2/(12*m)
            +h*h*DELTA/(12*m)+h*h*DELTA0/(6*M)
            +h**4*(3*M+2*m)/360)


def translated():
    return (m, z-c*m, q-2*c*z+c*c*m,
            r-3*c*q+3*c*c*z-c**3*m)


@lru_cache(maxsize=3)
def rank_reduction(n):
    if n < 1:
        raise ValueError('positive Hankel order required')
    y = s.symbols('y0:'+str(2*n))
    last = s.Symbol('last')
    b = exponential_coefficients(y+(last,))
    A = s.Matrix(n, n, lambda i,j: b[i+j])
    v = s.Matrix(b[n:2*n-1])
    a = A[:n-1, n-1]
    cc = A[n-1,n-1]
    dd = s.expand(b[2*n-1]-y[-1])
    ee = s.expand(b[2*n]-last+y[0]*y[-1])
    if n == 1:
        sigma = cc
        eta = -dd-y[0]*sigma/2
        reduced = -y[0]*dd-y[0]**2*cc/4-ee
    else:
        inverse = A[:n-1,:n-1].inv()
        sigma = cc-(a.T*inverse*a)[0]
        eta = (a.T*inverse*v)[0]-dd-y[0]*sigma/2
        w = v+y[0]*a/2
        reduced = (w.T*inverse*w)[0]-y[0]*dd-y[0]**2*cc/4-ee
    F = boundary_prediction(y)
    sigma, eta, reduced = map(s.factor, (sigma, eta, reduced))
    assert s.cancel(F-reduced-(y[-1]-eta)**2/sigma) == 0
    return dict(n=n, variables=y, F=F, sigma=sigma, eta=eta, reduced=reduced,
                determinants=tuple(s.factor(A[:j,:j].det()) for j in range(1,n+1)))


def certificate_squares():
    h = m-M
    A = DELTA0*N-N0*DELTA
    B = DELTA0*(M*q-m*Q+m*M*M*h/12)-N0*(M*z-m*Z)
    C = DELTA0*(M*z+(m-2*M)*Z)-N0*M*h
    D = DELTA0*((2*m-M)*z-m*Z)-N0*m*h
    # coefficient, generator at x, generator at u, polynomial square root
    return (
        (s.S.One, s.S.One, s.S.One, M*A),
        (s.S.One, DELTA, s.S.One, B),
        (s.Rational(1,6), m*DELTA, M, C),
        (s.Rational(1,12), DELTA, s.S.One, M*D),
        (s.Rational(1,12), s.S.One, s.S.One, DELTA*M*DELTA0*h),
        (s.Rational(1,6), m*DELTA, M*DELTA0, DELTA0*h),
        (s.Rational(1,120), m*DELTA, M, M*DELTA0*h*h),
        (s.Rational(1,180), DELTA, s.S.One, m*M*DELTA0*h*h),
    )


@lru_cache(maxsize=1)
def identities():
    L = fifth_schur()[2]
    translated_L = L.xreplace(dict(zip(X, translated())))
    residuals = {
        'schur_split': s.cancel(L-h0()-N*N/(m*DELTA)),
        'domain_invariance': s.expand(DELTA.xreplace(dict(zip(X, translated())))-DELTA),
        'numerator_translation': s.expand(N.xreplace(dict(zip(X, translated())))-N+2*c*DELTA),
        'boundary_translation': s.cancel(translated_L-L+4*c*r-6*c*c*q+4*c**3*z-c**4*m),
        'reduced_certificate': s.cancel(remainder(h0(), X[:3], U[:3])-h0_certificate()),
    }
    lhs = s.cancel(P*P0**2*remainder(L, X, U))
    rhs = sum(coef*gx*gu*f*f for coef, gx, gu, f in certificate_squares())
    residuals['full_joint_certificate'] = s.expand(lhs-rhs)
    return residuals


def records():
    checks = identities()
    assert all(v == 0 for v in checks.values()), [k for k,v in checks.items() if v != 0]
    terms = []
    for coef, gx, gu, f in certificate_squares():
        degree_x = s.Poly(gx*f*f, X).total_degree()
        terms.append(dict(coefficient=str(coef), x_generator=str(gx),
                          u_generator=str(gu), square=str(f), degree_x=degree_x))
    pL = s.cancel(P*fifth_schur()[2])
    assert not s.denom(pL).free_symbols
    reductions = []
    for n in range(1,4):
        data = rank_reduction(n)
        reductions.append({key:str(data[key]) for key in ('n','sigma','eta','reduced')})
    return dict(schema=1, exact_residuals={k:str(v) for k,v in checks.items()},
                denominator=str(P), terms=terms, pL_degree=s.Poly(pL, X).total_degree(),
                rank_reduction_calibrations=reductions,
                general_induction_proved=False, priority_cleared=False,
                optimization_calls=0, new_closed_loops=0)


if __name__ == '__main__':
    for key, value in identities().items():
        print(key, value == 0, flush=True)
    print('maximum certificate degree in x', max(t['degree_x'] for t in records()['terms']))
