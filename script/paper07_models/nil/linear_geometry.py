"""Linearized Ricci tensor in a noncoordinate frame with structure constants (used by operators.py)."""
from geometry import *

def red(M):return M.applyfunc(lambda x:s.cancel(s.expand(x)))

def ricci_linear(g,h,C,D0,D1,reduce=red):
    n=g.rows;inv=g.inv();sample=next(iter(h.values()));zero=s.zeros(*sample.shape)
    G0=connection(g,C,D0,reduce=simp)
    di={(a,b):-sum((inv[a,i]*h[i,j]*inv[j,b] for i,j in product(range(n),repeat=2)),zero) for a,b in product(range(n),repeat=2)}
    def kos0(e,a,b):return D0(a,g[b,e])+D0(b,g[a,e])-D0(e,g[a,b])-sum(C[f][b][e]*g[a,f]-C[f][e][a]*g[b,f]-C[f][a][b]*g[e,f] for f in range(n))
    def kos1(e,a,b):return D1(a,h[b,e])+D1(b,h[a,e])-D1(e,h[a,b])-sum((C[f][b][e]*h[a,f]-C[f][e][a]*h[b,f]-C[f][a][b]*h[e,f] for f in range(n)),zero)
    G1={(c,a,b):reduce(sum((inv[c,e]*kos1(e,a,b)+di[c,e]*kos0(e,a,b) for e in range(n)),zero)/2) for c,a,b in product(range(n),repeat=3)}
    R={}
    for c,b in product(range(n),repeat=2):
        R[c,b]=reduce(sum((D1(a,G1[a,b,c])-D1(b,G1[a,a,c])+sum((G1[e,b,c]*G0[a,a,e]+G0[e,b,c]*G1[a,a,e]-G1[e,a,c]*G0[a,b,e]-G0[e,a,c]*G1[a,b,e]-C[e][a][b]*G1[a,e,c] for e in range(n)),zero) for a in range(n)),zero))
    return R,G1
