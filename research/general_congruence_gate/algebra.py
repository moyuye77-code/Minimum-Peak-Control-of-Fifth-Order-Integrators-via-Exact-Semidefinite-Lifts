"""Exact algebra only; no SDP/nonrepresentability inference from samples."""
from functools import lru_cache
import sympy as s

m,z,q,r,v=s.symbols('m z q r v')
Y=(m,z,q,r,v)
c,a,b,f=s.symbols('c a b f')


def exponential_coefficients(moments):
    """b_j=[w^(j+1)] (1-exp(-sum y_j w^(j+1)))."""
    ee=[s.S.One]
    for n in range(1,len(moments)+1):
        ee.append(s.expand(-sum(k*moments[k-1]*ee[n-k]
                                for k in range(1,n+1))/n))
    return tuple(-x for x in ee[1:])


def hankel(moments,order,shift=0):
    assert len(moments)>=2*order-1+shift
    aa=exponential_coefficients(moments)
    return s.Matrix(order,order,lambda i,j:aa[i+j+shift])


@lru_cache(maxsize=1)
def raw_family():
    C=s.Matrix([[1,0,0],[c*m,1,0],[a*z+b*m*m,f*m,1]])
    return C,(C*hankel(Y,3)*C.T).applyfunc(s.expand)


def normalized_family(F=f):
    """F may be ANY polynomial, not only a weighted-homogeneous constant."""
    C=s.Matrix([[1,0,0],[m/2,1,0],
                [z-m*m/6-F*(z-m*m/2),m*F,1]])
    return C,(C*hankel(Y,3)*C.T).applyfunc(s.expand)


@lru_cache(maxsize=1)
def forced_matrix():
    return normalized_family(s.Rational(1,2))[1]


def linear_density_moments(intercept,slope,count=5):
    return tuple(s.Rational(intercept)/(j+1)+s.Rational(slope)/(j+2)
                 for j in range(count))


def midpoint_witness():
    plus=linear_density_moments(s.Rational(3,4),s.Rational(-3,8))
    minus=linear_density_moments(s.Rational(1,4),s.Rational(3,8))
    center=linear_density_moments(s.Rational(1,2),0)
    G=forced_matrix()
    at=lambda point:G.subs(dict(zip(Y,point)))
    defect=((at(plus)+at(minus))/2-at(center)).applyfunc(s.factor)
    return dict(plus=plus,minus=minus,center=center,defect=defect)


def interval_data(intervals,count):
    intervals=tuple((s.Rational(l),s.Rational(rr)) for l,rr in intervals)
    assert all(l<rr for l,rr in intervals)
    assert all(intervals[i][1]<intervals[i+1][0] for i in range(len(intervals)-1))
    moments=tuple(sum((rr**(j+1)-l**(j+1))/(j+1) for l,rr in intervals)
                  for j in range(count))
    lefts=[l for l,rr in intervals]
    weights=tuple(s.factor(-s.prod(l-rr for _,rr in intervals)/
                          s.prod(l-other for k,other in enumerate(lefts) if k!=i))
                  for i,l in enumerate(lefts))
    return moments,tuple(lefts),weights


def boundary_prediction(prefix):
    """Classical interior Schur formula. Caller must establish its domain."""
    d=len(prefix)+1;shift=1-d%2;n=(d-1)//2
    highest=s.Symbol('highest')
    H=hankel(tuple(prefix)+(highest,),n+1,shift)
    A=H[:n,:n];B=H[:n,n]
    offset=s.expand(H[n,n]-highest)
    if n==0:return s.factor(-offset)
    assert A.det()!=0,'The interior Schur formula must not divide by zero'
    return s.factor((B.T*A.inv()*B)[0]-offset)


def fifth_schur():
    G=forced_matrix();A=G[:2,:2];B=G[:2,2]
    L=s.factor(v-G[2,2]+(B.T*A.inv()*B)[0])
    return A,B,L


def records():
    C,G=raw_family();N,F=normalized_family();forced=forced_matrix()
    delta=m*q-z*z-m**4/12
    assert s.expand(F[1,2]-(r-m*m*z/6+(f-s.Rational(1,2))*delta))==0
    hessian=s.hessian(forced[2,2],(m,q))
    witness=midpoint_witness()
    assert witness['defect'][2,2]==s.Rational(43,6291456)>0
    assert s.factor(hessian.det())==-m*m/36
    A,B,L=fifth_schur()
    assert s.factor(L-boundary_prediction(Y[:4]))==0
    intervals=(((s.Rational(1,10),s.Rational(3,10)),
                (s.Rational(6,10),s.Rational(8,10))),
               ((s.Rational(0),s.Rational(1,10)),
                (s.Rational(3,10),s.Rational(5,10)),
                (s.Rational(7,10),s.Rational(9,10))))
    calibrations=[]
    for ii,(bands,d) in enumerate(zip(intervals,(5,6))):
        yy,nodes,weights=interval_data(bands,d)
        assert all(w>0 for w in weights)
        coeff=exponential_coefficients(yy)
        assert all(s.factor(coef-sum(w*x**j for w,x in zip(weights,nodes)))==0
                   for j,coef in enumerate(coeff))
        assert boundary_prediction(yy[:-1])==yy[-1]
        calibrations.append(dict(d=d,moments=list(map(str,yy)),
                                 nodes=list(map(str,nodes)),weights=list(map(str,weights))))
    mat=lambda mm:[[str(x) for x in row] for row in mm.tolist()]
    return dict(schema=1,coefficients=list(map(str,exponential_coefficients(Y))),
                generic_congruence=mat(C),generic_matrix=mat(G),
                normalized_congruence=mat(N),normalized_matrix=mat(F),
                forced_matrix=mat(forced),mixed_hessian=mat(hessian),
                mixed_hessian_determinant=str(s.factor(hessian.det())),
                midpoint={k:mat(val) if isinstance(val,s.MatrixBase) else list(map(str,val))
                          for k,val in witness.items()},
                fifth_schur=str(L),interval_calibrations=calibrations,
                scope='polynomial lower-triangular congruences with constant nonzero diagonal only',
                general_SDP_decided=False,novelty_cleared=False,
                optimization_calls=0,new_closed_loop_runs=0)
