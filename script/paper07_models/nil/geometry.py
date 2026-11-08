"""Exact frame and coordinate geometry. nabla_a E_b=Gamma[c,a,b] E_c."""
import sympy as s
from itertools import product
PAIRS=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
def simp(x):return s.factor(s.cancel(s.expand(x)))
def symmetric(v,dual=False):
    M=s.zeros(3)
    for (i,j),x in zip(PAIRS,v):M[i,j]=M[j,i]=x/(2 if dual and i!=j else 1)
    return M
def constants(n):return [[[ -sum(n[i,l]*s.LeviCivita(l,j,k) for l in range(3)) for k in range(3)] for j in range(3)] for i in range(3)]
def connection(g,C,D,reduce=simp):
    d=g.rows;inv=g.inv();G={}
    for c,a,b in product(range(d),repeat=3):
        G[c,a,b]=reduce(sum(inv[c,e]*(D(a,g[b,e])+D(b,g[a,e])-D(e,g[a,b])-sum(C[f][b][e]*g[a,f]-C[f][e][a]*g[b,f]-C[f][a][b]*g[e,f] for f in range(d))) for e in range(d))/2)
    return G
def ricci(g,C,D,reduce=simp):
    d=g.rows;G=connection(g,C,D,reduce)
    R=s.zeros(d)
    for c,b in product(range(d),repeat=2):
        R[c,b]=reduce(sum(D(a,G[a,b,c])-D(b,G[a,a,c])+sum(G[e,b,c]*G[a,a,e]-G[e,a,c]*G[a,b,e]-C[e][a][b]*G[a,e,c] for e in range(d)) for a in range(d)))
    return R,G
def linear_ricci(g,h,C,D0,D1,reduce=simp):
    d=g.rows;inv=g.inv();di=-inv*h*inv
    G0=connection(g,C,D0,reduce);G1={}
    def kos(M,D,e,a,b):
        return D(a,M[b,e])+D(b,M[a,e])-D(e,M[a,b])-sum(C[f][b][e]*M[a,f]-C[f][e][a]*M[b,f]-C[f][a][b]*M[e,f] for f in range(d))
    for c,a,b in product(range(d),repeat=3):
        G1[c,a,b]=reduce(sum(inv[c,e]*kos(h,D1,e,a,b)+di[c,e]*kos(g,D0,e,a,b) for e in range(d))/2)
    R=s.zeros(d)
    for c,b in product(range(d),repeat=2):
        R[c,b]=reduce(sum(D1(a,G1[a,b,c])-D1(b,G1[a,a,c])+sum(G1[e,b,c]*G0[a,a,e]+G0[e,b,c]*G1[a,a,e]-G1[e,a,c]*G0[a,b,e]-G0[e,a,c]*G1[a,b,e]-C[e][a][b]*G1[a,e,c] for e in range(d)) for a in range(d)))
    return R,G1
def backgrounds():
    u=s.Symbol('u',positive=True);t=s.Symbol('t',positive=True);z=s.Symbol('tau',real=True)
    A=1+u**4;F=s.exp(t*t)/s.sqrt(t)
    return {
        'Nil':dict(time=u,N=3*s.sqrt(A),g=s.diag(A,A,16*u*u/(9*A)),n=s.diag(0,0,-1),interval=['1','2'],parameter='b=1; literature time=u^3'),
        'Sol':dict(time=t,N=s.sqrt(F),g=s.diag(F,t,t),n=s.Matrix([[0,0,0],[0,0,1],[0,1,0]]),interval=['1','2'],parameter='C=1; fixed base length log((3+sqrt(5))/2)'),
        'S3':dict(time=z,N=s.sqrt((1+z*z)/(1-z*z)),g=s.diag(4*(1+z*z),4*(1+z*z),16*(1-z*z)/(1+z*z)),n=2*s.eye(3),interval=['1/4','1/2'],parameter='ell=1,m=0; sigma=-omega/2')}
def adm_data(bg):
    t,N,g=bg['time'],bg['N'],bg['g'];inv=g.inv()
    # Keep positive branches explicit for the chosen time intervals.
    vol=s.prod(s.sqrt(g[i,i]) for i in range(3))
    if bg['n']==2*s.eye(3):
        z=t;vol=16*s.sqrt((1-z*z)*(1+z*z))
    K=s.diff(g,t)/(2*N);tr=s.trace(inv*K)
    pi=(vol*(inv*K*inv-tr*inv)).applyfunc(s.simplify)
    return dict(volume=vol,K=K,traceK=s.simplify(tr),pi=pi)
