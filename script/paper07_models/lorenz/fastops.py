"""Exact coefficient arithmetic in local K=k exp(-xi), M=l exp(xi)."""
from optools import *
from sympy.polys.matrices import DomainMatrix
oldw,oldk,oldl=w,k,l
k,l=s.symbols('K M',real=True,nonzero=True);w=s.S(1)
domain=s.QQ.frac_field(t,F,k,l)
TT,FF,KK,LL=[domain.from_sympy(v) for v in (t,F,k,l)]
def dm(M):
    M=M.subs({oldw:1,oldk:k,oldl:l},simultaneous=True)
    return DomainMatrix([[domain.from_sympy(x) for x in row] for row in M.tolist()],M.shape,domain).to_sparse()
def red(M):return dm(M).to_Matrix()
def deriv(x,i):
    return x.field.new(x.numer.diff(i)*x.denom-x.numer*x.denom.diff(i),x.denom**2)
def _dx(A):return A.applyfunc(lambda x:-KK*deriv(x,2)+LL*deriv(x,3),domain)
def dx(M):return _dx(dm(M)).to_Matrix()
def dt(M):
    return dm(M).applyfunc(lambda x:deriv(x,0)+(2*TT-1/(2*TT))*FF*deriv(x,1),domain).to_Matrix()
def trim(a):
    a=[red(M) for M in a]
    while len(a)>1 and a[-1]==s.zeros(*a[-1].shape):a.pop()
    return a
def add(*ops):
    shape=ops[0][0].shape;z=dm(s.zeros(*shape))
    out=[]
    for j in range(max(map(len,ops))):
        q=z
        for a in ops:
            if j<len(a):q=q+dm(a[j])
        out.append(q.to_Matrix())
    return trim(out)
def compose(a,b):
    aa=list(map(dm,a));bb=list(map(dm,b))
    out=[dm(s.zeros(a[0].rows,b[0].cols)) for _ in range(len(a)+len(b)-1)]
    for i,A in enumerate(aa):
        for j,B in enumerate(bb):
            D=B
            for r in range(i+1):
                out[i+j-r]=out[i+j-r]+(A*D)*domain.convert(comb(i,r))
                D=_dx(D)
    return trim([M.to_Matrix() for M in out])
def adjoint(a):
    out=[dm(s.zeros(a[0].cols,a[0].rows)) for _ in a]
    for j,A in enumerate(a):
        D=dm(A.H)
        for r in range(j+1):
            out[j-r]=out[j-r]+D*domain.convert((-1)**j*comb(j,r))
            D=_dx(D)
    return trim([M.to_Matrix() for M in out])
def timed(a):return trim([dt(M) for M in a])
def iszero(a):return all(M.is_zero_matrix for M in map(red,a))
def scale(c,a):
    return trim([(dm(M)*domain.from_sympy(s.sympify(c))).to_Matrix() for M in a])
def extract(a,rows,cols):return trim([M.extract(rows,cols) for M in a])
def jet_dx(M,nvars):
    out=dx(M)
    assert M[:,-nvars:]==s.zeros(M.rows,nvars),'jet overflow'
    out[:,nvars:]+=M[:,:-nvars]
    return red(out)
def fromjets(M,nvars):
    assert M.cols%nvars==0
    return trim([M[:,j:j+nvars] for j in range(0,M.cols,nvars)])
