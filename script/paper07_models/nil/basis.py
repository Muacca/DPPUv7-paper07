"""Helicity tensor basis, oscillator functions and grade selections for the central-nonzero Nil sector. Writes basis.json."""
from common import *
from geometry import *

SPINS=[2,1,0,0,-1,-2]
VSPINS=[1,0,-1]

def helicity():
    vp=s.Matrix([1,s.I,0])/s.sqrt(2);vm=s.conjugate(vp);v0=s.Matrix([0,0,1])
    forms=[vp*vp.T,(vp*v0.T+v0*vp.T)/s.sqrt(2),(vp*vm.T+vm*vp.T)/s.sqrt(2),v0*v0.T,(vm*v0.T+v0*vm.T)/s.sqrt(2),vm*vm.T]
    U=s.Matrix.hstack(*[s.Matrix([h[i,j] for i,j in PAIRS]) for h in forms])
    return U,s.Matrix.hstack(vp,v0,vm)

def oscillator(size,b):
    a=s.zeros(size)
    for n in range(1,size):a[n-1,n]=s.sqrt(n)
    return [b*(a-a.T)/s.sqrt(2),s.I*b*(a+a.T)/s.sqrt(2),s.I*b*b*s.eye(size)]

def selections(L,size):
    metric=[(r,L-spin) for r,spin in enumerate(SPINS) if 0<=L-spin<size]
    mult=([(0,L)] if 0<=L<size else [])+[(1+r,L-spin) for r,spin in enumerate(VSPINS) if 0<=L-spin<size]
    return [r*size+n for r,n in metric],[r*size+n for r,n in mult],metric,mult

if __name__=='__main__':
    ck=Checks();b=s.Symbol('b',positive=True);U,V=helicity();W=s.diag(1,1,1,2,2,2)
    ck('tensor Frobenius normalization',s.simplify(U.H*W*U)==s.eye(6))
    ck('vector normalization',s.simplify(V.H*V)==s.eye(3))
    T=s.diag(U,U.inv().H);J=s.zeros(12);J[:6,6:]=s.eye(6);J[6:,:6]=-s.eye(6)
    ck('tensor and dual momentum canonical',s.simplify(T.H*J*T)==J)
    rot=s.Matrix([[0,1,0],[-1,0,0],[0,0,0]]) # components of L_R on frame covectors
    ck('vector weights',s.simplify(rot*V-V*s.diag(*[s.I*x for x in VSPINS]))==s.zeros(3))
    # Transform the tensor basis by the same frame action on both indices.
    tensor_rot=s.zeros(6)
    q=s.symbols('q0:6');hm=symmetric(q);lh=rot*hm+hm*rot.T
    tensor_rot=s.Matrix([lh[i,j] for i,j in PAIRS]).jacobian(q)
    ck('tensor weights',s.simplify(tensor_rot*U-U*s.diag(*[s.I*x for x in SPINS]))==s.zeros(6))
    D=oscillator(8,b);number=s.diag(*range(8))
    ck('oscillator rotation relation X',number*D[0]-D[0]*number==s.I*D[1])
    ck('oscillator rotation relation Y',number*D[1]-D[1]*number==-s.I*D[0])
    defect=D[0]*D[1]-D[1]*D[0]-D[2]
    ck('Heisenberg relation on guard interior',defect[:7,:7]==s.zeros(7))
    ck('truncation defect explicitly nonzero at upper edge',defect[7,7]==-8*s.I*b*b)
    x,mu=s.symbols('x mu',real=True);P=s.Symbol('P',positive=True);m,r,k=s.symbols('m r k',integer=True)
    # P=2pi, mu=m/P; reindex k -> k+1 in x deck transformation.
    ck.zero('deck oscillator argument reindex',x+P+(r+k*m)/(m/P)-(x+(r+(k+1)*m)/(m/P)))
    ck.zero('deck y phase reindex',m+r+k*m-(r+(k+1)*m))
    ck.zero('conjugation preserves oscillator argument',(-r-k*m)/(-mu)-(r+k*m)/mu)
    # A local horizontal rotation does not descend through the fixed lattice.
    xx,yy,zz=s.symbols('x y z',real=True);R=s.Matrix([-yy,xx,(xx*xx-yy*yy)/2]);coords=s.Matrix([xx,yy,zz]);deck=s.Matrix([xx+P,yy,zz+P*yy])
    descent=R.subs(dict(zip(coords,deck)),simultaneous=True)-deck.jacobian(coords)*R
    ck('grading is not a global rotation Killing field',descent!=s.zeros(3,1))
    blocks={str(L):dict(metric_slots=selections(L,8)[2],multiplier_slots=selections(L,8)[3]) for L in range(-2,4)}
    g1,g2=s.symbols('g1 g2',positive=True)
    scalar=-(D[0]*D[0]/g1+D[1]*D[1]/g2)
    ck.zero('nonLRS scalar ladder coupling',scalar[2,0]+b*b*s.sqrt(2)*(1/g1-1/g2)/2)
    save('basis.json',dict(checks=ck.items,tensor_transform=U,vector_transform=V,spins=SPINS,vector_spins=VSPINS,blocks=blocks,Heisenberg_upper_edge_defect=defect,rotation_deck_residual=descent,nonLRS_n0_to_n2=scalar[2,0],scalar_basis='P=2pi, mu=m/P !=0; Phi_(m,r,n)=sqrt(P) exp(i mu z) sum_(k in Z) exp(i(r+k m)y) f_n(x+(r+k m)/mu), where f_n is the real unit-L2 Hermite function of frequency abs(mu). r=0,...,abs(m)-1. Haar measure is normalized by P^4.',real_field='conjugate(Phi_m,r,n)=Phi_-m,(-r mod abs(m)),n. Tensor conjugation also reverses spin. Use L=n+sign(m)*spin, invariant under this conjugation.',closure_argument='[iN,X]=-Y and [iN,Y]=X. The diagonal tensor action plus iN makes gradient, symmetrized gradient, divergence, curvature variation and ADM maps equivariant when g and pi have equal horizontal axes. All tensor spins are bounded, so a fixed L has finitely many slots. This grading is an algebraic operator on each Schrödinger representation, not a descended rotation isometry.',boundary='The displayed finite oscillator matrices are only a guard implementation. A second-order operator acting on n<=L+2 reaches at most L+4; retain that entire range before projecting. Never impose a finite-dimensional Heisenberg commutator at the guard edge.',nonLRS_boundary='If g1!=g2, the scalar operator contains nonzero n<->n+2 links for arbitrarily high n. The present finite grading no longer closes; this is not a theorem excluding every possible alternative basis.'))
    print(len(ck.items),'basis checks',flush=True)
