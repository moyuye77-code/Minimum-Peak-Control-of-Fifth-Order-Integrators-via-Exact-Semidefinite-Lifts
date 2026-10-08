#!/usr/bin/env python3
"""Independent exact-algebra audit of representability-review.pdf.

Dependencies: Python 3.10+, SymPy. Run: python verify_derivations.py
This script is transcribed from the PDF, not imported from the author's archive.
It verifies algebraic claims only; it is not a formal proof of the analytic theorems.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as S

R = S.Rational
results: dict[str, object] = {"sympy_version": S.__version__, "checks": []}

def check(name: str, expression: S.Expr) -> None:
    residual = S.cancel(expression)
    ok = residual == 0
    results["checks"].append({"name": name, "passed": bool(ok)})
    print(('PASS' if ok else 'FAIL') + ': ' + name, flush=True)
    if not ok:
        raise AssertionError(f'{name}: residual = {S.factor(residual)}')

# Four-moment boundary and first-order identities.
s,t,w0 = S.symbols('s t w0')
M,Z,g = S.symbols('M Z g')
u,w = s-M,t-Z
D=t-s**2/2
D0=Z-M**2/2
Fg=t**2/2+s**2*t/4-g*(s**3/12+s*t/2)+g**2*(s**2/8-t/4)
sub2={s:M,t:Z}
remainder=Fg-Fg.xreplace(sub2)-S.diff(Fg,s).xreplace(sub2)*u-S.diff(Fg,t).xreplace(sub2)*w
rhs11=(w+(M-g)*u/2+u**2/6)**2/2+D*u**2/12+D0*u**2/6+u**4/36
check('Eq. (11), polynomial first-order identity', remainder-rhs11)

V=S.symbols('V')
B=w0-s*t/2-s**3/12
F0=t**2/2+s**2*t/4
P=D*F0+B**2
sub3={s:M,t:Z,w0:V}
B0=B.xreplace(sub3)
P0=P.xreplace(sub3)
# Polynomial-only evaluation of D(x)D(X)^2 times the rational Taylor remainder.
lhs14=P*D0**2-D*P0*D0-D*sum((S.diff(P,a).xreplace(sub3)*D0-P0*S.diff(D,a).xreplace(sub3))*(a-sub3[a]) for a in (s,t,w0))
rhs14=(D0*B-B0*D)**2+D*(D0*(w+M*u/2+u**2/6)-B0*u)**2/2+(D*D0*u)**2/12+D*D0*(D0*u)**2/6+D*(D0*u**2)**2/36
check('Eq. (14), rational first-order identity after clearing denominators', S.expand(lhs14-rhs14))

m,z,q,r,v=S.symbols('m z q r v')
y=(m,z,q,r)
D4=z-m**2/2
E4=m-z-m**2/2
B4=q-m*z/2-m**3/12
F04=z**2/2+m**2*z/4
L4=F04+B4**2/D4
bar={m:1-m,z:R(1,2)-z,q:R(1,3)-q,r:R(1,4)-r}
qmin=z**2/m+m**3/12
qmax=R(1,3)-(R(1,2)-z)**2/(1-m)-(1-m)**3/12
check('Eq. (6), q-plus minus q-minus', qmax-qmin-D4*E4/(m*(1-m)))
check('Eq. (6), complementary fourth-boundary factorization', L4+L4.xreplace(bar)-R(1,4)-m*(1-m)*(q-qmin)*(q-qmax)/(D4*E4))

# K4 lift, exact polynomial witnesses and complementary transformation.
a,b,c,d,e,f=S.symbols('a b c d e f')
H4=S.Matrix([[1,m,a],[m,a,b],[a,b,c]])
J4=lambda sig:S.Matrix([[1,m,z+sig*a/6],[m,a,d+sig*b/6],[z+sig*a/6,d+sig*b/6,f+sig*e/3+c/36]])
Tm=S.Matrix([[z-a/2,d-b/2],[d-b/2,e-c/2]])
Tp=S.Matrix([[m-z-a/2,a-d-b/2],[a-d-b/2,b-e-c/2]])
Qm=S.Matrix([[z-a/2,q-d/2-b/12],[q-d/2-b/12,r-f/2-e/4]])
Bp=-a/4+b/12-d/2+m/2-q+z/2
Cp=-a/8-d/2+e/4-f/2+m/4-r+3*z/4
Qp=S.Matrix([[m-z-a/2,Bp],[Bp,Cp]])
atomic={a:m**2,b:m**3,c:m**4,d:m*z,e:m**2*z,f:z**2}
G4=S.Matrix([[D4,B4],[B4,r-F04]])
expected4=[S.Matrix([1,m,m**2])*S.Matrix([1,m,m**2]).T]
expected4 += [S.Matrix([1,m,z+sig*m**2/6])*S.Matrix([1,m,z+sig*m**2/6]).T for sig in (1,-1)]
expected4 += [weight*S.Matrix([1,m])*S.Matrix([1,m]).T for weight in (D4,E4)]
expected4 += [G4,G4.xreplace(bar)]
actual4=[H4,J4(1),J4(-1),Tm,Tp,Qm,Qp]
for i,(mat,expect) in enumerate(zip(actual4,expected4),start=1):
    check(f'Theorem 2 atomic witness, block {i}', sum(S.expand(e0)**2 for e0 in (mat.xreplace(atomic)-expect)))
baraux={**bar,a:1-2*m+a,b:1-3*m+3*a-b,c:1-4*m+6*a-4*b+c,d:R(1,2)-z-m/2+d,e:R(1,2)-m+a/2-z+2*d-e,f:R(1,4)-z+f}
check('Theorem 2, complemented Q-minus is Q-plus', sum(S.expand(e0)**2 for e0 in Qm.xreplace(baraux)-Qp))

# Exponential-transform Hankel congruences.
wseries=S.symbols('wseries')
poly=m*wseries+z*wseries**2+q*wseries**3+r*wseries**4+v*wseries**5
expcoeff=[S.expand(sum((-1)**(j+1)*poly**j/S.factorial(j) for j in range(1,6))).coeff(wseries,k) for k in range(1,6)]
Hankel=S.Matrix(3,3,lambda i,j:expcoeff[i+j])
Tr=S.Matrix([[1,0,0],[m/2,1,0],[z/2+m**2/12,m/2,1]])
Gt=S.Matrix([[m,z,q],[z,q-m**3/12,r-m**2*z/6],[q,r-m**2*z/6,v-m**2*q/12-m*z**2/4+m**5/720]])
check('Lemma 4, all 9 Hankel congruence entries', sum(S.expand(e0)**2 for e0 in Tr*Hankel*Tr.T-Gt))
Tr2=S.Matrix([[1,0],[m/2,1]])
HH2=S.Matrix([[expcoeff[1],expcoeff[2]],[expcoeff[2],expcoeff[3]]])
check('Remark 3, all 4 Hankel congruence entries', sum(S.expand(e0)**2 for e0 in Tr2*HH2*Tr2.T-G4))

# Fifth-order construction.
Delta=m*q-z**2-m**4/12
N=m*r-z*q-m**3*z/6
p=m*Delta
HH=q**2/m+m**2*q/12+m*z**2/4-m**5/720
numer=S.cancel(p*HH+N**2)
Q,Rr=S.symbols('Q R')
ref={m:M,z:Z,q:Q,r:Rr}
Delta0=Delta.xreplace(ref)
N0=N.xreplace(ref)
p0=p.xreplace(ref)
num0=numer.xreplace(ref)
h=m-M
Astar=Delta0*N-N0*Delta
Bstar=Delta0*(M*q-m*Q+m*M**2*h/12)-N0*(M*z-m*Z)
Cstar=Delta0*(M*z+(m-2*M)*Z)-N0*M*h
Dstar=Delta0*((2*m-M)*z-m*Z)-N0*m*h
lhs16=numer*p0**2-p*num0*p0-p*sum((S.diff(numer,a).xreplace(ref)*p0-num0*S.diff(p,a).xreplace(ref))*(a-ref[a]) for a in y)
rhs16=(M*Astar)**2+Delta*Bstar**2+p*M*Cstar**2/6+Delta*(M*Dstar)**2/12+(Delta*M*Delta0*h)**2/12+p*(M*Delta0)*(Delta0*h)**2/6+p*M*(M*Delta0*h**2)**2/120+Delta*(m*M*Delta0*h**2)**2/180
check('Eq. (16), all coefficients with 8 independent indeterminates', S.expand(lhs16-rhs16))
rhs20=N**2+Delta*q**2+p**2/12+Delta*(m*z)**2/3+Delta*m**6/180
check('Eq. (20), homogeneous apex identity', numer-rhs20)
check('Lemma 4, Schur complement determinant', S.det(Gt)-Delta*(v-(HH+N**2/p)))

basis=[p,p*m,p*z,p*q,p*r,numer,N**2,N*Delta,Delta**2,p**2,p*Delta,Delta*q**2,Delta*m**2*q,Delta*q*z,Delta*m**3,Delta*m**4,Delta*m**2*z,Delta*m**2*z**2,p*z**2]
weights=[S.Integer(1),S.Integer(1),Delta,Delta,Delta,p,p]
features=[(N,Delta),(p,Delta),(q,m,m**2,z),(m*z,z,m,m**2),(m**3,m**2,m),(z,m,1),(m**2,m,1)]
weighted=[weight*S.Matrix(feat)*S.Matrix(feat).T for weight,feat in zip(weights,features)]
theta=S.symbols('theta0:19')
g1=theta[3]-theta[8]-theta[15]/12
g2=6*(theta[4]-theta[13]-theta[7])
g3=12*(theta[12]-theta[10]-theta[18])
g4=180*(theta[5]-theta[6]-theta[11])-15*theta[9]-60*theta[17]
pencils=[S.Matrix([[theta[6],theta[7]],[theta[7],theta[8]]]),S.Matrix([[theta[9],theta[10]],[theta[10],theta[8]]]),S.Matrix([[theta[11],theta[3],theta[12],theta[13]],[theta[3],theta[1],theta[14],theta[2]],[theta[12],theta[14],theta[15],theta[16]],[theta[13],theta[2],theta[16],g1]]),S.Matrix([[theta[17],theta[18],theta[16],g2],[theta[18],g1,theta[2],theta[16]],[theta[16],theta[2],theta[1],theta[14]],[g2,theta[16],theta[14],theta[15]]]),S.Matrix([[g4,g3,theta[15]],[g3,theta[15],theta[14]],[theta[15],theta[14],theta[1]]]),S.Matrix([[theta[18],theta[16],theta[2]],[theta[16],theta[14],theta[1]],[theta[2],theta[1],theta[0]]]),S.Matrix([[g3,theta[15],theta[14]],[theta[15],theta[14],theta[1]],[theta[14],theta[1],theta[0]]])]
subs_basis=dict(zip(theta,basis))
entry_count=0
for bi,(pencil,wp) in enumerate(zip(pencils,weighted),start=1):
    differences=pencil.xreplace(subs_basis)-wp
    for e0 in differences:
        assert S.expand(e0)==0, (bi,e0)
        entry_count+=1
    check(f'Section 3.2.1 affine coefficient interface block {bi}', sum(S.expand(e0)**2 for e0 in differences))
results['matrix_entries_checked']=entry_count
basis_polys=[S.Poly(e0,*y) for e0 in basis]
monomials=sorted(set().union(*(set(pp.monoms()) for pp in basis_polys)))
CM=S.Matrix([[pp.coeff_monomial(mon) for pp in basis_polys] for mon in monomials])
rank=CM.rank()
assert rank==19
results['basis_rank']=rank
results['coefficient_matrix_shape']=list(CM.shape)
print(f'PASS: coefficient space rank = {rank}; explicit entries = {entry_count}',flush=True)
# Each square in (16) is explicitly represented in the declared feature space.
coeffs=[S.Matrix([M*Delta0,-M*N0]),S.Matrix([Delta0*M,-Delta0*Q-Delta0*M**3/12+N0*Z,Delta0*M**2/12,-N0*M]),S.Matrix([Delta0*M,Delta0*Z-N0*M,-2*Delta0*M*Z+N0*M**2]),S.Matrix([2*M*Delta0,-M**2*Delta0,M*(N0*M-Delta0*Z),-M*N0]),S.Matrix([M*Delta0,-M**2*Delta0]),S.Matrix([0,Delta0,-M*Delta0]),S.Matrix([M*Delta0,-2*M**2*Delta0,M**3*Delta0]),S.Matrix([M*Delta0,-2*M**2*Delta0,M**3*Delta0])]
feature_ids=[0,2,5,3,1,5,6,4]
squares=[M*Astar,Bstar,Cstar,M*Dstar,Delta*M*Delta0*h,Delta0*h,M*Delta0*h**2,m*M*Delta0*h**2]
for i,(co,fi,sq) in enumerate(zip(coeffs,feature_ids,squares),start=1):
    check(f'Eq. (16) square {i} lies in retained feature space', (co.T*S.Matrix(features[fi]))[0]-sq)

# Singularity analysis: exact boundary witnesses, not floating samples.
eps,alpha=S.symbols('eps alpha',positive=True)
center={m:R(1,2),z:R(1,4),q:R(1,6),r:R(1,8)}
kappa=q-z**2/m-m**3/12
check('Boundary: kappa(center) = 1/32',kappa.xreplace(center)-R(1,32))
Hess=S.hessian(kappa,(m,z))
check('Boundary: determinant of (m,z) Hessian of kappa = 1', Hess.det()-1)
interval={m:m,z:alpha*m+m**2/2,q:alpha**2*m+alpha*m**2+m**3/3,r:alpha**3*m+R(3,2)*alpha**2*m**2+alpha*m**3+m**4/4}
interval_v=alpha**4*m+2*alpha**3*m**2+2*alpha**2*m**3+alpha*m**4+m**5/5
check('Boundary: single-interval Delta vanishes',Delta.xreplace(interval))
check('Boundary: single-interval N vanishes',N.xreplace(interval))
check('Boundary: single-interval H equals physical fifth moment',HH.xreplace(interval)-interval_v)
zero_path={key:eps*value for key,value in center.items()}
check('Boundary zero path: Delta exact order 2',Delta.xreplace(zero_path)-eps**2*(4-eps**2)/192)
check('Boundary zero path: N equals Delta',N.xreplace(zero_path)-Delta.xreplace(zero_path))
check('Boundary zero path: p exact order 3',p.xreplace(zero_path)-eps**3*(4-eps**2)/384)
limits_zero=[S.limit(S.cancel(el.xreplace(zero_path)/p.xreplace(zero_path)),eps,0,dir='+') for el in basis]
assert limits_zero==[S.Integer(1)]+[S.Integer(0)]*18
results['zero_path_basis_limits']=[str(x) for x in limits_zero]
print('PASS: all 19 selected functional values have exact zero-density limits',flush=True)
full={m:1,z:R(1,2),q:R(1,3),r:R(1,4)}
full_path={a:(1-eps)*full[a]+eps*center[a] for a in y}
limits_full=[S.limit(S.cancel(el.xreplace(full_path)/p.xreplace(full_path)),eps,0,dir='+') for el in basis]
assert all(el.is_finite for el in limits_full)
results['full_path_basis_limits']=[str(x) for x in limits_full]
print('PASS: all 19 selected functional values have finite full-density limits',flush=True)
# Generic single-interval witness: explicit limits, implied by N^2/p -> 0.
interval_theta=[1,m,interval[z],interval[q],interval[r],interval_v,0,0,0,0,0,interval[q]**2/m,m*interval[q],interval[q]*interval[z]/m,m**2,m**3,m*interval[z],m*interval[z]**2,interval[z]**2]
# Normalized first two blocks vanish; five other blocks are direct outer products.
for idx in range(7):
    if idx<2:
        expected=S.zeros(len(features[idx]))
    else:
        bb=S.Matrix(features[idx]).xreplace(interval)
        expected=bb*bb.T/(m if weights[idx]==Delta else 1)
    actual=pencils[idx].xreplace(dict(zip(theta,interval_theta)))
    check(f'Boundary symbolic interval witness block {idx+1}',sum(S.cancel(el)**2 for el in actual-expected))

# Bernstein coefficient and support-anchor calculations.
l,hseg=S.symbols('l hseg')
a0,a1,a2,a3,a4=S.symbols('a0 a1 a2 a3 a4')
tt=S.symbols('tt')
quartic=a0+a1*tt+a2*tt**2+a3*tt**3+a4*tt**4
power_local=S.Poly(S.expand(quartic.subs(tt,l+hseg*tt)),tt)
b2=sum(power_local.nth(k)*S.binomial(2,k)/S.binomial(4,k) for k in range(3))
check('Section 3.5 quartic Bernstein middle coefficient',b2-(quartic.subs(tt,l)+hseg*S.diff(quartic,tt).subs(tt,l)/2+hseg**2*S.diff(quartic,tt,2).subs(tt,l)/12))
roots=[(1-1/S.sqrt(2))/2,R(1,2),(1+1/S.sqrt(2))/2]
knots=[0,*roots,1]
anchor=[S.simplify(sum((-1)**i*S.integrate(384*(1-tt)**(3-j)/S.factorial(3-j),(tt,knots[i],knots[i+1])) for i in range(4))) for j in range(4)]
assert anchor==[1,0,0,0]
results['fourth_order_rest_to_rest_anchor']=list(map(str,anchor))
print('PASS: fourth-order rest-to-rest control peak 384 reaches (1,0,0,0)',flush=True)
kernel=384*(1-tt)**3/6-192*(1-tt)**2/2+40*(1-tt)-4
absintegral=S.simplify(sum((-1)**i*S.integrate(kernel,(tt,knots[i],knots[i+1])) for i in range(4)))
check('Section 4.4 stated dual absolute kernel integral equals 1',absintegral-1)

results['all_checks_passed']=all(row['passed'] for row in results['checks'])
results['number_of_named_checks']=len(results['checks'])
out=Path(__file__).with_name('symbolic_audit_results.json')
out.write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
print(f'Wrote {out}; all checks passed.',flush=True)
