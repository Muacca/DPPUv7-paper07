"""Independent four-dimensional check of the central-nonzero Nil reconstruction at grade L=2: all ten Einstein components. Writes einstein_check.json."""
from reduce_blocks import *

if __name__=='__main__':
    L=2;size=7;u,b,raw,a=load_block(L,size);H,C,J=a['H'],a['C'],a['J'];P=a['original_embedding'];M=a['multiplier_embedding'];D=oscillator(size,b)
    nq=H.rows//2;nz=2*nq;nr=C.rows;nt=nz+3*nr;ck=Checks();_,g,pi,nv,nov=bdata();N2=9*(1+u**4)
    evo=s.zeros(nt);evo[:nz,:nz]=J*H;evo[:nz,nz:nz+nr]=J*C.H
    evo[nz:nz+nr,nz+nr:nz+2*nr]=s.eye(nr);evo[nz+nr:nz+2*nr,nz+2*nr:]=s.eye(nr)
    z=s.zeros(size,nt);h={(a0,b0):z.copy() for a0,b0 in product(range(4),repeat=2)}
    for r,(i,k) in enumerate(PAIRS):
        hh=z.copy();hh[:,:nz]=P[r*size:(r+1)*size,:];h[i+1,k+1]=h[k+1,i+1]=hh
    h[0,0][:,nz:nz+nr]=-2*N2*M[:size,:]
    for i in range(3):
        hh=z.copy();hh[:,nz:nz+nr]=g[i,i]*M[(i+1)*size:(i+2)*size,:];h[0,i+1]=h[i+1,0]=hh
    g4=s.diag(-N2,*g.diagonal());cs=constants(s.diag(0,0,-1));C4=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for a0,b0,c0 in product(range(3),repeat=3):C4[a0+1][b0+1][c0+1]=cs[a0][b0][c0]
    D0=lambda i,f:s.diff(f,u) if i==0 else 0
    D1=lambda i,f:red(s.diff(f,u)+f*evo) if i==0 else D[i-1]*f
    print('Independent full Einstein check L=2',flush=True)
    R,_=ricci_linear(g4,h,C4,D0,D1);scalar=red(sum((g4.inv()[i,k]*R[k,i] for i,k in product(range(4),repeat=2)),z))
    G={(i,k):red(R[i,k]-g4[i,k]*scalar/2) for i,k in product(range(4),repeat=2)}
    for i,k in PAIRS:ck('Einstein spatial '+str((i,k)),G[i+1,k+1]==z)
    CC=red(M*C).row_join(s.zeros(4*size,3*nr))
    ck('Hamiltonian Einstein identity',red(G[0,0]+CC[:size,:]*N2/(2*nv))==z)
    for i in range(3):ck('momentum Einstein identity '+str(i),red(G[0,i+1]+CC[(i+1)*size:(i+2)*size,:]*nov/2)==z)
    save('einstein_check.json',dict(checks=ck.items,L=L,frequency='mu=b^2>0 arbitrary',scope='All ten four-dimensional Einstein components, with every tensor component in the closed block and all lapse/shift time jets. Spatial Einstein equations vanish on ADM evolution; temporal components equal the original constraints.'))
    print(len(ck.items),'four-dimensional checks',flush=True)
