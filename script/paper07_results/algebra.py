"""Geometric torsion, Cartan forms, and the complete connection equation.

Signature (-,+,+,+); epsilon_0123=+1; all antisymmetric indices are summed.
The 24 K coordinates are K_ab,c with a<b and both Lorentz indices lowered.
"""
from itertools import combinations,product,permutations
import sympy as s
SIGN=(-1,1,1,1)
def clean(f):return {k:s.factor(v) for k,v in f.items() if s.expand(v)!=0}
def add(*fs):
    out={}
    for f in fs:
        for k,v in f.items():out[k]=out.get(k,0)+v
    return clean(out)
def scale(a,f):return clean({k:a*v for k,v in f.items()})
def wedge(f,g):
    out={}
    for k,a in f.items():
        for l,b in g.items():
            kl=k+l
            if len(set(kl))<len(kl):continue
            sign=(-1)**sum(kl[i]>kl[j] for i in range(len(kl)) for j in range(i+1,len(kl)))
            idx=tuple(sorted(kl));out[idx]=out.get(idx,0)+sign*a*b
    return clean(out)
def exterior(f,de,derivative):
    out={}
    for idx,v in f.items():
        out=add(out,wedge({(0,):derivative(v)},{idx:1}))
        for pos,a in enumerate(idx):out=add(out,scale((-1)**pos*v,wedge(wedge({idx[:pos]:1},de[a]),{idx[pos+1:]:1})))
    return out
def contortion():
    slots=[(a,b,c) for a,b in combinations(range(4),2) for c in range(4)]
    y=s.symbols('K0:24',real=True);kd=dict(zip(slots,y))
    def K(a,b,c):return 0 if a==b else (1 if a<b else -1)*kd[min(a,b),max(a,b),c]
    rk=s.expand(sum(SIGN[a]*SIGN[b]*SIGN[e]*(K(a,e,a)*K(e,b,b)-K(a,e,b)*K(e,b,a)) for a,b,e in product(range(4),repeat=3)))
    torsion=[SIGN[a]*(K(a,d,c)-K(a,c,d)) for a in range(4) for c,d in combinations(range(4),2)]
    M=s.hessian(rk,y);A=s.Matrix(torsion).jacobian(y)
    return dict(coordinates=y,scalar=rk,hessian=M,torsion_map=A,
                current=[s.expand(2*SIGN[a]*sum(SIGN[b]*K(a,b,b) for b in range(4))) for a in range(4)])
def class_a_forms():
    q=s.Symbol('q',positive=True);H,HD,eta,etaD,V,VD=s.symbols('H HD eta etaD V VD',real=True)
    n=s.symbols('n1:4',real=True);D=lambda x:sum(s.diff(x,a)*b for a,b in [(q,H*q),(H,HD),(eta,etaD),(V,VD)])
    e=[{(i,):s.S.One} for i in range(4)];cyc=[(1,2,3),(2,3,1),(3,1,2)]
    de=[{}]+[add(scale(H,wedge(e[0],e[i])),scale(n[i-1]/q,wedge(e[j],e[k]))) for i,j,k in cyc]
    om=[[{} for _ in range(4)] for _ in range(4)]
    for i in range(1,4):om[0][i]=om[i][0]=scale(H-V/3,e[i])
    for i,j,k in permutations(range(1,4)):om[j][k]=scale(s.LeviCivita(i-1,j-1,k-1)*(sum(n)/2-n[i-1]-eta)/q,e[i])
    T=[add(de[a],*[wedge(om[a][b],e[b]) for b in range(4)]) for a in range(4)]
    R=[[add(exterior(om[a][b],de,D),*[wedge(om[a][c],om[c][b]) for c in range(4)]) for b in range(4)] for a in range(4)]
    def comp(f,a,b):return 0 if a==b else (1 if a<b else -1)*f.get(tuple(sorted((a,b))),0)
    scalar=s.factor(sum(SIGN[b]*comp(R[a][b],a,b) for a,b in product(range(4),repeat=2)))
    Q=add(*[scale(SIGN[a],wedge(e[a],T[a])) for a in range(4)])
    tt=add(*[scale(SIGN[a],wedge(T[a],T[a])) for a in range(4)])
    eer=add(*[scale(SIGN[a],wedge(wedge(e[a],e[b]),R[a][b])) for a,b in product(range(4),repeat=2)])
    eps=add(*[scale(s.LeviCivita(a,b,c,d)*SIGN[d],wedge(wedge(e[a],e[b]),R[c][d])) for a,b,c,d in permutations(range(4))])
    return dict(q=q,H=H,HD=HD,eta=eta,etaD=etaD,V=V,VD=VD,n=n,R=scalar,Q=Q,NY=exterior(Q,de,D),TT=tt,EER=eer,epsilon_EER=eps,T=T)
