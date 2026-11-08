"""Matrix differential operators with coefficients on the left of d/dxi."""
from common import *
import sympy as s
from math import comb
t,F,w=s.symbols('t F w',positive=True)
k,l=s.symbols('k l',real=True,nonzero=True)
def red(M):return M.applyfunc(lambda q:s.cancel(s.expand(q)))
def dx(M):return M.diff(w)*w
def dt(M):return red(M.diff(t)+(2*t-1/(2*t))*F*M.diff(F))
def trim(a):
    a=[red(M) for M in a]
    while len(a)>1 and a[-1]==s.zeros(*a[-1].shape):a.pop()
    return a
def zero(m,n):return [s.zeros(m,n)]
def add(*ops):
    n=max(map(len,ops));shape=ops[0][0].shape
    return trim([sum((a[j] for a in ops if j<len(a)),s.zeros(*shape)) for j in range(n)])
def scale(c,a):return trim([c*M for M in a])
def compose(a,b):
    out=[s.zeros(a[0].rows,b[0].cols) for _ in range(len(a)+len(b)-1)]
    for i,A in enumerate(a):
        for j,B in enumerate(b):
            D=B
            for r in range(i+1):
                out[i+j-r]+=comb(i,r)*A*D
                D=dx(D)
    return trim(out)
def adjoint(a):
    out=[s.zeros(a[0].cols,a[0].rows) for _ in range(len(a))]
    for j,A in enumerate(a):
        D=A.H
        for r in range(j+1):
            out[j-r]+=(-1)**j*comb(j,r)*D
            D=dx(D)
    return trim(out)
def timed(a):return trim([dt(M) for M in a])
def iszero(a):return all(M==s.zeros(*M.shape) for M in trim(a))
def extract(a,rows,cols):return trim([M.extract(rows,cols) for M in a])
def asjets(a,nvars,njet):
    assert a[0].cols==nvars and len(a)<=njet
    return s.Matrix.hstack(*(a[j] if j<len(a) else s.zeros(a[0].rows,nvars) for j in range(njet)))
def fromjets(M,nvars):
    assert M.cols%nvars==0
    return trim([M[:,j:j+nvars] for j in range(0,M.cols,nvars)])
def jet_dx(M,nvars):
    out=dx(M)
    assert M[:,-nvars:]==s.zeros(M.rows,nvars),'jet overflow'
    out[:,nvars:]+=M[:,:-nvars]
    return red(out)
def parse_ops(raw):
    env={str(v):v for v in (t,F,w,k,l)}
    return [s.Matrix([[s.sympify(v,locals=env) for v in row] for row in A]) for A in raw]
