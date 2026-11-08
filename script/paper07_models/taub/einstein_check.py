"""Independent four-dimensional Einstein tensor on full ADM evolution.

All metric components, lapse and shift, and their time jets remain explicit.
The ADM evolution must give G_ij=0 and G_0mu proportional to the constraints.
"""
from reduce_block import *

if __name__=='__main__':
    j=s.Rational(sys.argv[1]);t,data=read_operators(j);H,C,J=data['H'],data['C'],data['J'];D,_,labels=generators(j);d=len(labels)
    nq=6*d;nz=12*d;n=24*d;ck=Checks();_,g,pi,nv,nov=bdata();N2=(1+t*t)/(1-t*t)
    evo=s.zeros(n);evo[:nz,:nz]=J*H;evo[:nz,nz:nz+4*d]=J*C.H
    evo[nz:nz+4*d,nz+4*d:nz+8*d]=s.eye(4*d);evo[nz+4*d:nz+8*d,nz+8*d:]=s.eye(4*d)
    h={};zero=s.zeros(d,n)
    for a,b in product(range(4),repeat=2):h[a,b]=zero.copy()
    for a,(i,l) in enumerate(PAIRS):
        M=zero.copy();M[:,a*d:(a+1)*d]=s.eye(d);h[i+1,l+1]=h[l+1,i+1]=M
    h[0,0][:,nz:nz+d]=-2*N2*s.eye(d)
    for i in range(3):
        M=zero.copy();M[:,nz+(i+1)*d:nz+(i+2)*d]=g[i,i]*s.eye(d);h[0,i+1]=h[i+1,0]=M
    g4=s.diag(-N2,*g.diagonal());C4=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)];cs=constants(2*s.eye(3))
    for a,b,c in product(range(3),repeat=3):C4[a+1][b+1][c+1]=cs[a][b][c]
    D0=lambda a,f:s.diff(f,t) if a==0 else 0
    D1=lambda a,M:red(s.diff(M,t)+M*evo) if a==0 else D[a-1]*M
    print('j',j,'four dimensional Einstein tensor',flush=True)
    ric,_=ricci_linear(g4,h,C4,D0,D1);scalar=red(sum((g4.inv()[a,b]*ric[b,a] for a,b in product(range(4),repeat=2)),zero))
    ein={(a,b):red(ric[a,b]-g4[a,b]*scalar/2) for a,b in product(range(4),repeat=2)}
    for i,l in PAIRS:ck(f'spatial Einstein {i}{l} on full ADM evolution',ein[i+1,l+1]==zero)
    padded=C.row_join(s.zeros(4*d,n-nz))
    # H uses the fractional lapse alpha, so C_H=N times the ADM scalar.
    ck('Einstein 00 equals Hamiltonian constraint',red(ein[0,0]+padded[:d,:]*N2/(2*nv))==zero)
    for a in range(3):ck(f'Einstein 0{a+1} equals momentum constraint',red(ein[0,a+1]+padded[(a+1)*d:(a+2)*d,:]*nov/2)==zero)
    save('einstein_j'+str(j).replace('/','_')+'.json',dict(checks=ck.items,j=str(j),Einstein_constraint_relation='G00=-N^2 C_H/(2 N sqrt(g)); G0i=-N C_i/(2 sqrt(g)); Gij=0 on the full linear ADM evolution, before setting constraints to zero. All lapse/shift time jets retained and cancel appropriately.',scope='Exact rational-time matrix identities for the complete right-index and metric block; the independent spectator left multiplicity is unchanged.'))
    print('Einstein checks passed',len(ck.items),flush=True)
