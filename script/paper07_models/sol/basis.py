"""Sol fiber lattice, fixed monodromy and the dual-lattice orbit basis on the unfolded line. Writes basis_v2.json."""
from optools import *
import numpy as np
from scipy.integrate import quad
from scipy.special import eval_hermite,factorial

if __name__=='__main__':
    ck=Checks();phi=(1+s.sqrt(5))/2;lam=(3+s.sqrt(5))/2;L=s.log(lam);c=s.sqrt(1+phi**2)
    A=s.Matrix([[2,1],[1,1]]);S=s.Matrix([[1,-phi],[phi,1]])/c;D=s.diag(1/lam,lam)
    ck('oriented area-one fiber lattice',s.simplify(S.det())==1)
    ck('orthogonal eigen-coordinate choice',s.simplify(S*S.T)==s.eye(2))
    ck('fixed monodromy intertwines lattice',s.simplify(S*A-D*S)==s.zeros(2))
    a,b=s.symbols('a b',integer=True);q=2*s.pi*S*s.Matrix([a,b])
    ck('dual lattice phases integral',s.simplify(S.T*q-2*s.pi*s.Matrix([a,b]))==s.zeros(2,1))
    Q=a*a-a*b-b*b
    ck('arithmetic orbit invariant',s.expand(Q.subs({a:2*a+b,b:a+b},simultaneous=True)-Q)==0)
    ck('product invariant and normalization',s.simplify(q[0]*q[1]-4*s.pi**2*Q/s.sqrt(5))==0)
    xx,yy,zz=s.symbols('x y z',real=True);j=s.Symbol('j',integer=True)
    ck('orbit sum deck reindex',s.simplify((xx+L)+j*L-(xx+(j+1)*L))==0)
    ck('fiber wavevector deck reindex',s.simplify(D.T*s.Matrix([k*lam**(-j),l*lam**j])-s.Matrix([k*lam**(-j-1),l*lam**(j+1)]))==s.zeros(2,1))
    ck('X2 orbit intertwiner',s.simplify(s.exp(-xx)*k*lam**(-j)-k*s.exp(-(xx+j*L)))==0)
    ck('X3 orbit intertwiner',s.simplify(s.exp(xx)*l*lam**j-l*s.exp(xx+j*L))==0)
    d1=[s.zeros(1),s.ones(1)];d2=[s.Matrix([[s.I*k/w]])];d3=[s.Matrix([[s.I*l*w]])]
    ck('Sol bracket X1 X2',iszero(add(compose(d1,d2),scale(-1,compose(d2,d1)),d2)))
    ck('Sol bracket X1 X3',iszero(add(compose(d1,d3),scale(-1,compose(d3,d1)),scale(-1,d3))))
    ck('Sol bracket X2 X3',iszero(add(compose(d2,d3),scale(-1,compose(d3,d2)))))
    for i,d in enumerate((d1,d2,d3)):ck('anti-adjoint generator '+str(i),iszero(add(adjoint(d),d)))
    CK=len(ck.items)
    pn=float(phi);ln=float(L);sn=np.array(S.evalf(),float);qn=2*np.pi*sn@np.array([1.,0.])
    def herm(n,x):
        return np.pi**(-.25)*eval_hermite(n,x)*np.exp(-x*x/2)/np.sqrt(2.**n*factorial(n))
    def mode(n,x,y,z,cut=12,sign=1):
        jj=np.arange(-cut,cut+1);weights=herm(n,x+jj*ln)
        return np.sqrt(ln)*np.sum(weights*np.exp(1j*sign*(qn[0]*np.exp(-jj*ln)*y+qn[1]*np.exp(jj*ln)*z)))
    errs=[]
    for n in (0,1,3,6):
        x,y,z=.23,.19,-.12
        norm=quad(lambda z:herm(n,z)**2,-np.inf,np.inf,epsabs=1e-12)[0]
        ck('Hermite norm n='+str(n),abs(norm-1)<1e-11)
        v=mode(n,x,y,z)
        errors=[abs(mode(n,x+ln,y/float(lam),z*float(lam))-v),abs(mode(n,x,y,z,sign=-1)-np.conj(v))]
        for col in range(2):errors.append(abs(mode(n,x,y+sn[0,col],z+sn[1,col])-v))
        errors.append(abs(mode(n,x,y,z,cut=14)-v));errs+=errors
        ck('quotient and conjugacy controls n='+str(n),max(errors)<2e-10)
    coeffs=[]
    for n in range(9):
        value=quad(lambda x:herm(n,x)*np.exp(x)*herm(0,x),-12,12,epsabs=1e-12)[0]
        expected=np.exp(.25)*(1/np.sqrt(2))**n/np.sqrt(factorial(n))
        ck('nonzero exponential coupling n='+str(n),abs(value-expected)<1e-11 and value>0)
        coeffs.append(value)
    save('basis_v2.json',dict(checks=ck.items,exact_checks=CK,numerical_checks=len(ck.items)-CK,lattice=S,A=A,lambda_=lam,base_length=L,representative=[1,0],dual_wavevector=q.subs({a:1,b:0}),orbit_invariant=Q,arithmetic_product=4*s.pi**2*Q/s.sqrt(5),operators=[d1,d2,d3],max_quotient_error=max(errs),exp_couplings_n0=coeffs,
        basis='U_q f = sqrt(L) sum_j f(x+j L) exp(i[(lambda^-j k)y+(lambda^j l)z]); f is a function on R. With normalized quotient Haar measure dx dy dz/L, U_q is isometric from L2(R). Fixed Hermite functions give an orthonormal basis for this orbit.',
        domain='Smooth invariant core: intersection over all integers a,b>=0 of {f: exp(a abs(xi)) derivative^b f in L2(R)}. Compactly supported smooth functions and Hermite functions lie in it. No periodic or quasiperiodic condition is imposed on f on R.',
        conjugacy='U_-q conjugate(f)=conjugate(U_q f). Positive-eigenvalue monodromy never identifies q and -q for q!=0. One complete real sector uses both entire orbits.',
        closure='All invariant frame tensors descend unchanged. D1=d/dxi, D2=i k exp(-xi), D3=i l exp(xi) act within the same orbit. All homogeneous-coefficient gravity operators preserve an orbit but generally couple infinitely many fixed Hermite functions.',
        arithmetic='For nonzero integer (a,b), a^2-a*b-b^2 is nonzero since sqrt(5) is irrational. Both eigen-direction wavevector components are nonzero. Representative changes correspond to translation of f by an integer multiple of L.',
        no_finite_Hermite_closure='exp(xi) psi_0 has coefficient exp(1/4) 2^(-n/2)/sqrt(n!) for every n>=0. A finite scalar Hermite cutoff is not an exact gravity block. This does not classify every possible tensor-dependent basis.',
        reference='Bolsinov-Dullin-Veselov, arXiv:math-ph/0503046v1, section 5, equation 19 and Theorem 2; convention and lattice derived here.'))
    print(len(ck.items),'basis checks',flush=True)
