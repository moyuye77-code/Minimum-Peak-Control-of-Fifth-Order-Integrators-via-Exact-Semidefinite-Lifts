"""Check manuscript pencils independently. Requires SymPy, never an optimizer.

No imports from research constructions. Checks algebra/transcription, not
the analytic boundary or nonrepresentability proof.
"""
from pathlib import Path
import json
import re
import sympy as s


def polynomial_data():
    x = s.symbols('m z q r')
    m, z, q, r = x
    delta = m*q-z*z-m**4/12
    n = m*r-z*q-m**3*z/6
    p = m*delta
    h = q*q/m+m*m*q/12+m*z*z/4-m**5/720
    f = h+n*n/p
    basis = (p,p*m,p*z,p*q,p*r,s.cancel(p*f),n*n,n*delta,delta*delta,
             p*p,p*delta,delta*q*q,delta*m*m*q,delta*q*z,delta*m**3,
             delta*m**4,delta*m*m*z,delta*m*m*z*z,p*z*z)
    features = ((1,(n,delta)),(1,(p,delta)),(delta,(q,m,m*m,z)),
                (delta,(m*z,z,m,m*m)),(delta,(m**3,m*m,m)),
                (p,(z,m,1)),(p,(m*m,m,1)))
    return x,delta,n,p,f,basis,features


def parse_pencils(source):
    section = re.split(r'\\(?:sub)*section(?:\*)?\{',
                       source.split(r'\label{app:coefficients}', 1)[1], maxsplit=1)[0]
    theta = s.symbols('theta0:19')
    names = {str(v):v for v in theta}

    def expression(raw):
        raw = re.sub(r'\\theta_\{(\d+)\}',r'theta\1',raw)
        raw = re.sub(r'\\gamma_\{(\d+)\}',r'gamma\1',raw).strip()
        tokens = re.findall(r'theta\d+|gamma\d+|\d+|[()+*/-]',raw)
        assert ''.join(tokens) == re.sub(r'\s+','',raw), raw
        assert all(t in names for t in tokens if t.startswith(('theta','gamma')))
        raw = re.sub(r'(\d)(\()',r'\1*\2',raw)
        raw = re.sub(r'(\d)(theta)',r'\1*\2',raw)
        return s.sympify(raw,locals=names)

    for number,raw in re.findall(r'^\\gamma_\{(\d+)\}&=(.*)$',section,re.M):
        names['gamma'+number] = expression(raw.rstrip('\\,. '))
    assert set(names) == {str(t) for t in theta}|{'gamma'+str(i) for i in range(1,5)}
    matches = re.findall(r'M_\{(\d+)\}\(\\theta\)&=\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}',section,re.S)
    assert [number for number,_ in matches] == list(map(str,range(1,8)))
    mats = [s.Matrix([[expression(cell) for cell in row.strip().split('&')]
                      for row in raw.strip().split(r'\\')]) for _,raw in matches]
    return theta,mats


def check_pencils(source):
    x,delta,n,p,f,basis,features = polynomial_data()
    theta,mats = parse_pencils(source)
    assert len(basis) == len(theta) == 19
    for mat,(weight,vector) in zip(mats,features):
        assert mat == mat.T and mat.rows == len(vector)
        assert all(s.Poly(e,theta).total_degree() <= 1 for e in mat)
        expanded = mat.xreplace(dict(zip(theta,basis)))
        desired = weight*s.Matrix(vector)*s.Matrix(vector).T
        assert all(s.expand(e) == 0 for e in expanded-desired), 'incorrect pencil coefficient'
    support = sorted({a for b in basis for a in s.Poly(b,*x).monoms()})
    coefficients = s.Matrix([[s.Poly(b,*x).coeff_monomial(a) for a in support] for b in basis])
    assert coefficients.rank() == 19
    m,z,q,r = x
    assert s.cancel(p*f-(n*n+delta*q*q+p*p/12+delta*(m*z)**2/3+delta*m**6/180)) == 0
    return dict(functional_rank=19,pencil_orders=[v.rows for v in mats],
                matrix_entries_checked=sum(v.rows**2 for v in mats),apex_identity=True)


def check_remainder():
    x,delta,n,p,f,_,_ = polynomial_data()
    m,z,q,r = x
    u = s.symbols('M Z Q R');M,Z,Q,R = u
    replace = dict(zip(x,u));d0=delta.xreplace(replace);n0=n.xreplace(replace)
    pu=p.xreplace(replace);fu=f.xreplace(replace);h=m-M
    a=d0*n-n0*delta
    b=d0*(M*q-m*Q+m*M*M*h/12)-n0*(M*z-m*Z)
    c=d0*(M*z+(m-2*M)*Z)-n0*M*h
    d=d0*((2*m-M)*z-m*Z)-n0*m*h
    rhs=(M*a)**2+delta*b*b+p*M*c*c/6+delta*(M*d)**2/12
    rhs+=(delta*M*d0*h)**2/12+p*(M*d0)*(d0*h)**2/6
    rhs+=p*M*(M*d0*h*h)**2/120+delta*(m*M*d0*h*h)**2/180
    remainder=f-fu-sum(s.diff(f,xi).xreplace(replace)*(xi-ui) for xi,ui in zip(x,u))
    numerator=s.together(p*pu*pu*remainder-rhs).as_numer_denom()[0]
    assert s.Poly(s.expand(numerator),*x,*u).is_zero
    return True


def main():
    source=Path(__file__).with_name('representability.tex').read_text(encoding='utf-8')
    report=check_pencils(source)
    report.update(joint_remainder_identity=check_remainder(),
                  optimizer_calls=0,boundary_proof_automatically_verified=False,
                  nonrepresentability_proof_automatically_verified=False)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
