"""Fixed-coframe geometric comparisons, separate from gravitational modes."""
import sympy as s
from .homogeneous import curvature,shape_residual
NAMES={'T3':s.zeros(3),'Nil':s.diag(0,0,-1),'Sol':s.diag(0,1,-1),'S3':2*s.eye(3)}
def rows():
    return [dict(name=name,n=n,rank=n.rank(),constant_closed_dimension=3-n.rank(),R=curvature(n,s.eye(3)),X=shape_residual(n)) for name,n in NAMES.items()]
def paths():
    z=s.Symbol('s',real=True);nil=NAMES['Nil'];sol=s.Matrix([[0,0,0],[0,0,1],[0,1,0]])
    specs=[('torus_shape',s.zeros(3),(s.Rational(1,2),-s.Rational(1,2),0)),('nil_flat',nil,(s.Rational(1,2),-s.Rational(1,2),0)),('nil_vertical',nil,(-s.Rational(1,2),-s.Rational(1,2),1)),('sol_base',sol,(s.Rational(1,3),-s.Rational(1,6),-s.Rational(1,6))),('sol_isometry',sol,(0,1,-1))]
    out=[]
    for name,n,v in specs:
        E=s.diag(*[s.exp(a*z) for a in v]);g=E.T*E;R=s.simplify(curvature(n,g));m=s.simplify(E*n*E.T/E.det())
        W=s.trace(m)*s.eye(3)/2-m
        out.append(dict(name=name,exponents=v,metric=g,orthonormal_structure=m,lc_rotation_coefficients=W,R=R,difference=s.simplify(R-R.subs(z,0)),derivative=s.diff(R,z),series=s.series(R-R.subs(z,0),z,0,5).removeO()))
    return out
def spectra():
    ell=s.Symbol('l',integer=True,positive=True);k=s.Symbol('k',integer=True,nonzero=True);j=s.Symbol('j',integer=True,nonnegative=True);q=s.Symbol('q',positive=True)
    m,n=s.symbols('m n',integer=True)
    return dict(sphere_rough=(ell*(ell+2)-2)/q**2,sphere_multiplicity=2*(ell-1)*(ell+3),
                nil_toroidal=(m*m+n*n)/q**2,nil_heisenberg=(k*k/(4*s.pi**2)+s.Abs(k)*(2*j+1)/(2*s.pi))/q**2,
                nil_signed_k_multiplicity=s.Abs(k))
