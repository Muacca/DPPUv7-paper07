"""Independent coordinate Ricci expansion about Minkowski space, through order two."""
import sympy as s
from itertools import product
ETA=s.diag(-1,1,1,1)
def ricci_linear(H,coords):
    G={}
    for a,b,c in product(range(4),repeat=3):
        G[a,b,c]=sum(ETA[a,d]*(s.diff(H[d,c],coords[b])+s.diff(H[d,b],coords[c])-s.diff(H[b,c],coords[d]))/2 for d in range(4))
    return s.Matrix(4,4,lambda b,c:sum(s.diff(G[a,b,c],coords[a])-s.diff(G[a,b,a],coords[c]) for a in range(4)))
def ricci_quadratic(H,coords):
    G1={};G2={};raised=ETA*H*ETA
    low={(d,b,c):(s.diff(H[d,c],coords[b])+s.diff(H[d,b],coords[c])-s.diff(H[b,c],coords[d]))/2 for d,b,c in product(range(4),repeat=3)}
    for a,b,c in product(range(4),repeat=3):
        G1[a,b,c]=sum(ETA[a,d]*low[d,b,c] for d in range(4))
        G2[a,b,c]=-sum(raised[a,d]*low[d,b,c] for d in range(4))
    R=s.Matrix(4,4,lambda b,c:sum(s.diff(G2[a,b,c],coords[a])-s.diff(G2[a,b,a],coords[c]) for a in range(4))+sum(G1[a,a,d]*G1[d,b,c]-G1[a,c,d]*G1[d,b,a] for a,d in product(range(4),repeat=2)))
    return R.applyfunc(s.expand)
def einstein_quadratic(H,coords):
    R1=ricci_linear(H,coords);R2=ricci_quadratic(H,coords)
    trace1=s.trace(ETA*R1)
    trace2=s.trace(ETA*R2)-s.trace(ETA*H*ETA*R1)
    return (R2-ETA*trace2/2-H*trace1/2).applyfunc(s.expand)
