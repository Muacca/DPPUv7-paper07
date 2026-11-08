"""Graded Hamiltonian and lapse/shift constraint operators of the central-nonzero Nil sector on a guarded oscillator range. Writes operators_g<size>.json."""
from common import *
from linear_geometry import *
from basis import *

def bdata():
    u=s.Symbol('u',positive=True);A=1+u**4
    g=s.diag(A,A,16*u*u/(9*A));pi=s.diag(-4/(9*A),-4/(9*A),-u*u)
    return u,g,pi,4*u*A,9/(4*u)

def build(size,ck):
    t,g,pi,nv,nov=bdata();inv=g.inv();b=s.Symbol('b',positive=True);D=oscillator(size,b);d=size;nq=6*d
    h={};zero=s.zeros(d,nq)
    for a,(i,l) in enumerate(PAIRS):
        M=s.zeros(d,nq);M[:,a*d:(a+1)*d]=s.eye(d);h[i,l]=h[l,i]=M
    C=constants(s.diag(0,0,-1));tr=lambda M:sum((inv[i,l]*M[l,i] for i,l in product(range(3),repeat=2)),zero)
    R0,_=ricci(g,C,lambda a,f:0);scalar0=s.trace(inv*R0)
    print('guard size',size,'spatial variation',flush=True)
    R1,_=ricci_linear(g,h,C,lambda a,f:0,lambda a,f:D[a]*f)
    scalar1=red(tr(R1)-sum(((inv*R0*inv)[i,l]*h[i,l] for i,l in product(range(3),repeat=2)),zero));trh=tr(h)
    Ein0=inv*R0*inv-scalar0*inv/2
    forces=[]
    for i,l in PAIRS:
        de=sum((inv[i,a]*R1[a,b]*inv[b,l]-(inv*R0*inv)[b,l]*inv[i,a]*h[a,b]-(inv*R0*inv)[i,a]*h[a,b]*inv[b,l] for a,b in product(range(3),repeat=2)),zero)
        de+=-scalar1*inv[i,l]/2+scalar0*sum((inv[i,a]*h[a,b]*inv[b,l] for a,b in product(range(3),repeat=2)),zero)/2+trh*Ein0[i,l]/2
        forces.append(red(nv*(1 if i==l else 2)*de))
    V=s.Matrix.vstack(*forces)
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
    defect=red(H-H.H)
    return dict(time=t,frequency_sqrt=b,g=g,pi=pi,Nvolume=nv,N_over_volume=nov,D=D,H=H,C=CM,J=J,guard_Hermiticity_defect=defect)

def project(data,size,ck):
    U,V=helicity();Q=s.kronecker_product(U,s.eye(size));T=s.diag(Q,Q.inv().H)
    row=s.kronecker_product(s.diag(s.ones(1),V.H),s.eye(size))
    H=red(T.H*data['H']*T);C=red(row*data['C']*T)
    qgrade=s.diag(*[spin+n for spin in SPINS for n in range(size)])
    grade=s.diag(qgrade,qgrade);cgrade=s.diag(*[spin+n for spin in [0]+VSPINS for n in range(size)])
    interior=[r*size+n for r in range(12) for n in range(size-2)]
    grade_defect=red(grade*H-H*grade)
    ck('tensor Hamiltonian grading exact on guard interior',grade_defect[:,interior]==s.zeros(12*size,len(interior)))
    ck('all lapse and shift constraints respect grading',red(cgrade*C-C*grade)==s.zeros(4*size,12*size))
    blocks={}
    for L in range(-2,3):
        qi,ci,metric,mult=selections(L,size);ni=len(qi);zi=qi+[6*size+i for i in qi]
        otherz=[a for a in range(12*size) if a not in zi];otherc=[a for a in range(4*size) if a not in ci]
        h=H.extract(zi,zi);c=C.extract(ci,zi);J=s.zeros(2*ni);J[:ni,ni:]=s.eye(ni);J[ni:,:ni]=-s.eye(ni)
        ck(f'L={L} Hamiltonian has no discarded component',H.extract(otherz,zi)==s.zeros(len(otherz),2*ni))
        ck(f'L={L} constraints have no discarded target',C.extract(otherc,zi)==s.zeros(len(otherc),2*ni))
        ck(f'L={L} gauge map has no discarded component',C.extract(ci,otherz)==s.zeros(len(ci),len(otherz)))
        ck(f'L={L} constraints first class',red(c*J*c.H)==s.zeros(len(ci)))
        ck(f'L={L} projected Hamiltonian Hermitian',h==h.H)
        ck(f'L={L} guard defect annihilates physical block',red(data['guard_Hermiticity_defect']*T[:,zi])==s.zeros(12*size,2*ni))
        blocks[str(L)]=dict(L=L,metric_slots=metric,multiplier_slots=mult,H=h,C=c,J=J,original_embedding=T[:,zi],multiplier_embedding=row.H[:,ci],metric_dimension_complex=ni,constraint_rows=len(ci))
    return blocks

if __name__=='__main__':
    size=int(sys.argv[1]) if len(sys.argv)>1 else 7;ck=Checks();data=build(size,ck);blocks=project(data,size,ck)
    save('operators_g'+str(size)+'.json',dict(checks=ck.items,guard_size=size,data=data,blocks=blocks,scope='Positive central frequency mu=b^2, b>0. All oscillator paths up to two derivatives retained before projection. L=-2,-1,0,1,2 are exact closed blocks, not numerical oscillator truncations.'))
    print(len(ck.items),'operator and closure checks',flush=True)
