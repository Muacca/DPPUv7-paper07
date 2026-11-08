"""Exact ADM operators with noncommuting Wigner derivative matrices."""
from common import *
from geometry import *
from harmonics import generators

def red(M):return M.applyfunc(lambda x:s.cancel(s.expand(x)))
def bdata():
    t=s.Symbol('tau',real=True);A=1+t*t;B=1-t*t
    g=s.diag(4*A,4*A,16*B/A);inv=g.inv()
    nv=16*A;nov=1/(16*B);von=16*B
    pi=red(von/2*(inv*s.diff(g,t)*inv-s.trace(inv*s.diff(g,t))*inv))
    return t,g,pi,nv,nov

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

def build(j,ck):
    t,g,pi,nv,nov=bdata();inv=g.inv();D,realS,labels=generators(j);d=len(labels);nq=6*d
    h={};zero=s.zeros(d,nq)
    for a,(i,l) in enumerate(PAIRS):
        M=s.zeros(d,nq);M[:,a*d:(a+1)*d]=s.eye(d);h[i,l]=h[l,i]=M
    C=constants(2*s.eye(3));tr=lambda M:sum((inv[i,l]*M[l,i] for i,l in product(range(3),repeat=2)),zero)
    R0,_=ricci(g,C,lambda a,f:0);scalar0=s.trace(inv*R0)
    print('j',j,'spatial variation',flush=True)
    R1,_=ricci_linear(g,h,C,lambda a,f:0,lambda a,f:D[a]*f)
    scalar1=red(tr(R1)-sum(((inv*R0*inv)[i,l]*h[i,l] for i,l in product(range(3),repeat=2)),zero));trh=tr(h)
    Ein0=inv*R0*inv-scalar0*inv/2
    forces=[]
    for i,l in PAIRS:
        de=sum((inv[i,a]*R1[a,b]*inv[b,l]-(inv*R0*inv)[b,l]*inv[i,a]*h[a,b]-(inv*R0*inv)[i,a]*h[a,b]*inv[b,l] for a,b in product(range(3),repeat=2)),zero)
        de+=-scalar1*inv[i,l]/2+scalar0*sum((inv[i,a]*h[a,b]*inv[b,l] for a,b in product(range(3),repeat=2)),zero)/2+trh*Ein0[i,l]/2
        forces.append(red(nv*(1 if i==l else 2)*de))
    V=s.Matrix.vstack(*forces);ck('potential Hermitian',red(V-V.H)==s.zeros(nq))
    hs=s.symbols('h0:6');ps=s.symbols('p0:6');hm=symmetric(hs);pm=symmetric(ps,True)
    S0=pi*g;S1=pm*g+pi*hm;S2=pm*hm
    T0=s.trace(S0*S0)-s.trace(S0)**2/2
    T1=2*s.trace(S0*S1)-s.trace(S0)*s.trace(S1)
    T2=s.trace(S1*S1)+2*s.trace(S0*S2)-s.trace(S1)**2/2-s.trace(S0)*s.trace(S2)
    hh=s.trace(inv*hm)
    Hkin=red(s.hessian(nov*(T2-hh*T1/2+(hh**2/8+s.trace(inv*hm*inv*hm)/4)*T0),hs+ps))
    H=s.Matrix(s.kronecker_product(Hkin,s.eye(d)));H[:nq,:nq]+=V
    scalar_kin=s.Matrix([nov*(T1-hh*T0/2)]).jacobian(hs+ps)
    CH=s.Matrix(s.kronecker_product(red(scalar_kin),s.eye(d)));CH[:,:nq]-=nv*(scalar1+trh*scalar0/2)
    hfull={key:M.row_join(s.zeros(d,nq)) for key,M in h.items()};pfull={}
    for a,(i,l) in enumerate(PAIRS):
        M=s.zeros(d,2*nq);M[:,nq+a*d:nq+(a+1)*d]=s.eye(d)/(1 if i==l else 2);pfull[i,l]=pfull[l,i]=M
    zero2=s.zeros(d,2*nq);mom=[]
    for a in range(3):
        M=D[a]*sum((pi[i,l]*hfull[i,l] for i,l in product(range(3),repeat=2)),zero2)
        M-=2*sum((D[i]*(pi[i,l]*hfull[a,l]+pfull[i,l]*g[a,l]) for i,l in product(range(3),repeat=2)),zero2)
        M-=sum((pi[i,l]*(C[b][a][i]*hfull[b,l]+C[b][a][l]*hfull[i,b])+pfull[i,l]*(C[b][a][i]*g[b,l]+C[b][a][l]*g[i,b]) for i,l,b in product(range(3),repeat=3)),zero2)
        mom.append(red(M))
    CM=red(s.Matrix.vstack(CH,*mom));H=red(H)
    J=s.zeros(2*nq);J[:nq,nq:]=s.eye(nq);J[nq:,:nq]=-s.eye(nq)
    ck('linear constraints first class',red(CM*J*CM.H)==s.zeros(4*d))
    ck('full Hamiltonian Hermitian',red(H-H.H)==s.zeros(2*nq))
    return dict(j=str(j),time=t,labels=labels,g=g,pi=pi,Nvolume=nv,N_over_volume=nov,D=D,reality=realS,H=H,C=CM,J=J,spatial_Ricci_variation=[R1[a,b] for a,b in PAIRS])

if __name__=='__main__':
    j=s.Rational(sys.argv[1]);ck=Checks();r=build(j,ck)
    save('operators_j'+str(j).replace('/','_')+'.json',dict(checks=ck.items,result=r))
    print('saved',j,flush=True)
