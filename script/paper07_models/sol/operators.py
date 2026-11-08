"""Exact differential operators of a general Sol dual orbit: Hamiltonian and the four constraints, formal self-adjointness and constraint propagation. Writes operators_exact.json."""
from optools import *
from geometry import PAIRS,symmetric,constants,ricci
from linear_geometry import ricci_linear
from itertools import product

def background():
    return s.diag(F,t,t),s.diag(-1/F,-t-1/(4*t),-t-1/(4*t)),t*F,1/t
def build(ck):
    g,pi,nv,nov=background();inv=g.inv();nvars=6;njet=5
    c=constants(s.Matrix([[0,0,0],[0,0,1],[0,1,0]]))
    def derivative(a,M):
        return jet_dx(M,M.cols//njet) if a==0 else red((s.I*k/w if a==1 else s.I*l*w)*M)
    zeroq=s.zeros(1,nvars*njet);h={}
    for q,(i,j) in enumerate(PAIRS):
        M=zeroq.copy();M[q]=1;h[i,j]=h[j,i]=M
    print('Spatial Ricci variation',flush=True)
    R0,_=ricci(g,c,lambda a,f:0)
    R1,_=ricci_linear(g,h,c,lambda a,f:0,derivative)
    scalar0=s.trace(inv*R0)
    tr=lambda M:sum((inv[i,j]*M[j,i] for i,j in product(range(3),repeat=2)),zeroq)
    scalar1=red(tr(R1)-sum(((inv*R0*inv)[i,j]*h[i,j] for i,j in product(range(3),repeat=2)),zeroq));trh=tr(h)
    Ein0=inv*R0*inv-scalar0*inv/2
    forces=[]
    for i,j in PAIRS:
        de=sum((inv[i,a]*R1[a,b]*inv[b,j]-(inv*R0*inv)[b,j]*inv[i,a]*h[a,b]-(inv*R0*inv)[i,a]*h[a,b]*inv[b,j] for a,b in product(range(3),repeat=2)),zeroq)
        de+=-scalar1*inv[i,j]/2+scalar0*sum((inv[i,a]*h[a,b]*inv[b,j] for a,b in product(range(3),repeat=2)),zeroq)/2+trh*Ein0[i,j]/2
        forces.append(red(nv*(1 if i==j else 2)*de))
    V=fromjets(s.Matrix.vstack(*forces),6)
    hs=s.symbols('h0:6');ps=s.symbols('p0:6');hm=symmetric(hs);pm=symmetric(ps,True)
    S0=pi*g;S1=pm*g+pi*hm;S2=pm*hm
    T0=s.trace(S0*S0)-s.trace(S0)**2/2
    T1=2*s.trace(S0*S1)-s.trace(S0)*s.trace(S1)
    T2=s.trace(S1*S1)+2*s.trace(S0*S2)-s.trace(S1)**2/2-s.trace(S0)*s.trace(S2)
    hh=s.trace(inv*hm)
    Hkin=red(s.hessian(nov*(T2-hh*T1/2+(hh**2/8+s.trace(inv*hm*inv*hm)/4)*T0),hs+ps))
    H=[Hkin]+[s.zeros(12) for _ in range(len(V)-1)]
    for j,M in enumerate(V):H[j][:6,:6]+=M
    scalar_kin=red(s.Matrix([nov*(T1-hh*T0/2)]).jacobian(hs+ps))
    CH=asjets([scalar_kin],12,njet)
    Rop=fromjets(nv*(scalar1+trh*scalar0/2),6)
    for j,M in enumerate(Rop):CH[:,j*12:j*12+6]-=M
    hfull={};pfull={};zz=s.zeros(1,12*njet)
    for q,(i,j) in enumerate(PAIRS):
        M=zz.copy();M[q]=1;hfull[i,j]=hfull[j,i]=M
        M=zz.copy();M[6+q]=s.S(1)/(1 if i==j else 2);pfull[i,j]=pfull[j,i]=M
    mom=[]
    for a in range(3):
        M=derivative(a,sum((pi[i,j]*hfull[i,j] for i,j in product(range(3),repeat=2)),zz))
        M-=2*sum((derivative(i,pi[i,j]*hfull[a,j]+pfull[i,j]*g[a,j]) for i,j in product(range(3),repeat=2)),zz)
        M-=sum((pi[i,j]*(c[b][a][i]*hfull[b,j]+c[b][a][j]*hfull[i,b])+pfull[i,j]*(c[b][a][i]*g[b,j]+c[b][a][j]*g[i,b]) for i,j,b in product(range(3),repeat=3)),zz)
        mom.append(red(M))
    C=fromjets(s.Matrix.vstack(CH,*mom),12);H=trim(H)
    ck('all coefficients exact with no Float atoms',not any(M.has(s.Float) for M in H+C))
    J=s.zeros(12);J[:6,6:]=s.eye(6);J[6:,:6]=-s.eye(6)
    ck('full Hamiltonian formally self-adjoint',iszero(add(H,scale(-1,adjoint(H)))))
    ck('all four constraints first class as differential operators',iszero(compose(compose(C,[J]),adjoint(C))))
    B=[s.zeros(4),s.zeros(4)]
    B[0][0,0]=t-1/(4*t);B[1][0,1]=1;B[0][0,2]=F*s.I*k/(t*w);B[0][0,3]=F*s.I*l*w/t
    ck('time-dependent full constraint propagation',iszero(add(timed(C),compose(compose(C,[J]),H),scale(-1,compose(B,C)))))
    gauge=compose([J],adjoint(C))
    ck('gauge intertwiner with time-dependent multipliers',iszero(add(timed(gauge),scale(-1,compose(compose([J],H),gauge)),scale(-1,compose(gauge,adjoint(B))))))
    ck('Hamiltonian order two',len(H)==3)
    ck('constraints order two',len(C)==3)
    ck('fiber-zero sector is distinct',all(not M.has(k,l) for M in [A.subs({k:0,l:0}) for A in C]))
    return dict(H=H,C=C,J=J,gauge=gauge,propagation=B,metric=g,pi=pi,Nvolume=nv,N_over_volume=nov,spatial_R=scalar0)

if __name__=='__main__':
    ck=Checks();data=build(ck)
    save('operators_exact.json',dict(checks=ck.items,data=data,convention='F=exp(t^2)/sqrt(t), w=exp(xi); left coefficients multiply powers of d/dxi. k,l are nonzero real components of a dual-lattice representative including 2pi. Variables (h11,h22,h33,h12,h13,h23,p11,p22,p33,p12,p13,p23), offdiagonal canonical p=2pi. Lapse is deltaN/N.',scope='Exact full infinite orbit operators, without spatial truncation. Same orbit is closed under evolution, constraints and gauge. No finite scalar mode closure is assumed.'))
    print(len(ck.items),'full orbit operator checks',flush=True)
