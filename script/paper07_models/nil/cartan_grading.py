"""Cartan (connection-elimination) map in the central Nil grading: invertibility, rotation covariance and helicity multiplicities. Writes cartan_grading.json."""
from common import *
from geometry import *
from itertools import combinations
sys.path.insert(0,str(REPO/'script'))
from paper07_results import algebra

if __name__=='__main__':
    ck=Checks();R=s.zeros(4);R[1,2]=1;R[2,1]=-1
    slots=[(a,b,c) for a,b in combinations(range(4),2) for c in range(4)]
    w=list(s.symbols('w0:24'));wd=dict(zip(slots,w))
    val=lambda a,b,c:0 if a==b else (1 if a<b else -1)*wd[min(a,b),max(a,b),c]
    rw=s.Matrix([sum(R[a,e]*val(e,b,c)+R[b,e]*val(a,e,c)+R[c,e]*val(a,b,e) for e in range(4)) for a,b,c in slots]).jacobian(w)
    t=list(s.symbols('t0:24'));ts=[(a,b,c) for a in range(4) for b,c in combinations(range(4),2)];td=dict(zip(ts,t))
    tv=lambda a,b,c:0 if b==c else (1 if b<c else -1)*td[a,min(b,c),max(b,c)]
    rt=s.Matrix([sum(R[a,e]*tv(e,b,c)+R[b,e]*tv(a,e,c)+R[c,e]*tv(a,b,e) for e in range(4)) for a,b,c in ts]).jacobian(t)
    ref=algebra.contortion();T=ref['torsion_map'];M=s.hessian(ref['scalar'],ref['coordinates'])
    ck('Cartan algebra map is invertible',T.det()!=0)
    ck('Cartan map intertwines all tensor rotation components',rt*T==T*rw)
    ck('contortion quadratic form rotation invariant',rw.T*M+M*rw==s.zeros(24))
    eigen=rw.eigenvals();expected={-2*s.I:2,-s.I:6,0:8,s.I:6,2*s.I:2}
    ck('connection helicity multiplicities',eigen==expected)
    counts={str(L):sum(mult for spin,mult in [(-2,2),(-1,6),(0,8),(1,6),(2,2)] if L-spin>=0) for L in range(-2,3)}
    save('cartan_grading.json',dict(checks=ck.items,connection_rotation=rw,torsion_rotation=rt,Cartan_map=T,spin_multiplicity={str(k/s.I):v for k,v in eigen.items()},connection_components_per_positive_m_r_block=counts,argument='The pointwise Cartan map preserves tensor spin and oscillator index and is invertible. It therefore restricts to an invertible map in each finite total-grade block, including low-end blocks. The derivative source is equivariant by the same ladder/tensor grading as the metric gradient. Eliminating auxiliary contortion adds no propagating connection pair.',endpoint_scope='The exact complete NY current e^a wedge T_a and the EH contortion current are retained as boundary representatives. Their difference shifts the pulled one-form by delta B and leaves the two-form unchanged. No fresh explicit block endpoint-current matrix or full unreduced Dirac chain is claimed in this first implementation.'))
    print(len(ck.items),'Cartan grading checks',flush=True)
