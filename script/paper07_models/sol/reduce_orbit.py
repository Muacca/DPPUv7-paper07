"""Differential right inverse, gauge and canonical embedding of a general Sol orbit. Writes reduction.json."""
from operators import *
from fastops import *

def load_data():
    a=json.loads((OUT/'operators_exact.json').read_text(encoding='utf8'))['data']
    U=s.diag(1,1,1,s.I,s.I,1);T=s.diag(U,U);L=s.diag(1,1,-s.I,-s.I)
    return trim([(T.H*M*T).applyfunc(s.expand) for M in parse_ops(a['H'])]),trim([(L*M*T).applyfunc(s.expand) for M in parse_ops(a['C'])]),parse_ops([a['J']])[0]

def lifts(ck,C):
    K=k/w;M=l*w;V=K*K+M*M;d=2-F*V/(2*t);A=F*(4*t*t-1)/(2*t)
    rows=s.eye(6).row_join(s.zeros(6,12));f,b,r0,r1,r2,r3=[rows[i,:] for i in range(6)]
    deriv=lambda x:jet_dx(x,6)
    aa=red((-deriv(b)-(K*K-M*M)*r0+A*(K*K-M*M)*f-V*r1/(2*t)-F*V*deriv(f)/t-K*r2/t+M*r3/t)/d)
    p22=red((r0-A*f)/2+r1/(4*t)+F*deriv(f)/(2*t)-F*aa/(4*t))
    p33=red((r0-A*f)/2-r1/(4*t)-F*deriv(f)/(2*t)+F*aa/(4*t))
    p12=red((aa+b)/(2*K));p13=red((aa-b)/(2*M))
    p23=red((-r2/t-deriv(p12)-p12-2*K*p22)/M)
    U=fromjets(s.Matrix.vstack(f,p22,p33,p12,p13,p23),6)
    E=extract(U,range(6),[0,1]);R=extract(U,range(6),range(2,6))
    T=s.zeros(6,2);T[0,0]=1;T[3,1]=K;T[4,1]=-M
    Cp=extract(C,range(4),range(6,12));Cq=extract(C,range(4),range(6))
    ck('explicit momentum kernel',iszero(compose(Cp,E)))
    ck('explicit differential right inverse',iszero(add(compose(Cp,R),[-s.eye(4)])))
    ck('physical coordinates paired with momentum amplitudes',iszero(add(compose([T.H],E),[-s.eye(2)])))
    ck('particular momentum is orthogonal to physical configurations',iszero(compose([T.H],R)))
    ck('complete differential canonical splitting',iszero(add(compose(E,[T.H]),compose(R,Cp),[-s.eye(6)])))
    P=scale(-1,compose(compose(R,Cq),[T]))
    lift=[]
    for j in range(max(len(E),len(P))):
        A0=s.zeros(12,4)
        if j==0:A0[:6,:2]=T
        if j<len(P):A0[6:,:2]=P[j]
        if j<len(E):A0[6:,2:]=E[j]
        lift.append(A0)
    gauge=[]
    for M0 in adjoint(R):
        G=s.zeros(4,12);G[:,:6]=M0;gauge.append(G)
    return trim(lift),gauge,dict(momentum_kernel=E,momentum_right_inverse=R,configuration_lift=T,particular_momentum=P,denominator=d)

if __name__=='__main__':
    class VerboseChecks(Checks):
        def __call__(self,name,ok):
            super().__call__(name,ok);print('PASS',name,flush=True)
    ck=VerboseChecks();H,C,J=load_data()
    print('Constructing exact differential constraint reduction',flush=True)
    E,G,parts=lifts(ck,C)
    Jr=s.zeros(4);Jr[:2,2:]=s.eye(2);Jr[2:,:2]=-s.eye(2)
    ck('all constraints eliminated',iszero(compose(C,E)))
    ck('differential gauge fixed',iszero(compose(G,E)))
    ck('gauge-constraint matrix is identity',iszero(add(compose(compose(G,[J]),adjoint(C)),[-s.eye(4)])))
    ck('reduced symplectic form canonical',iszero(add(compose(compose(adjoint(E),[-J]),E),[Jr])))
    corr=compose(compose(adjoint(E),[-J]),timed(E))
    ck('time-dependent pullback correction evaluated',iszero(corr))
    ck('gauge rows commute',iszero(compose(compose(G,[J]),adjoint(G))))
    ck('reduced kinetic correction is self-adjoint',iszero(add(corr,scale(-1,adjoint(corr)))))
    ck('uniform denominator bound negative',bool((2-4*s.pi**2*s.E/s.sqrt(5)).is_negative))
    ck('F/t increasing for t>=1',(2*t-3/(2*t)).subs(t,1)>0 and s.diff(2*t-3/(2*t),t).is_positive)
    ck('no Floats in reduction',not any(M.has(s.Float) for M in E+G))
    save('reduction.json',dict(checks=ck.items,embedding=E,gauge=G,reduced_J=Jr,time_correction=corr,parts=parts,
        reduced_H_factorization='Hr=E^dagger H_work E-E^dagger Omega Edot; the second term is exactly zero in this gauge. This factorization is the exact saved differential operator, not a finite-mode matrix.',
        multiplier_factorization='lambda=-(Gdot E+G J H_work E)y in working constraint phases. The gauge-constraint operator G J C^dagger equals the identity.',
        reconstruction_proof=[
            'Let R=J H E+J C^dagger lambda_op-Edot-E Jr Hr.',
            'C R=0 follows from Cdot+C J H=B C, C E=0, and C J C^dagger=0.',
            'G R=0 follows from G J C^dagger=I, the displayed multiplier, and d(G E)/dt=0.',
            'E^dagger Omega R=0 follows from the definition of Hr and E^dagger Omega E=Omega_r.',
            'The verified complete momentum splitting E_p T^dagger+R_p C_p=I and its adjoint imply ker(C) intersect ker(G)=range(E).',
            'Thus R=E v and Omega_r v=0, implying R=0. This proves every full evolution component without expanding the large factorized residual.',
        ],
        computation_scope='Full C/H differential identities, explicit right inverse, canonical pullback and time correction were checked coefficient by coefficient. A direct expanded residual calculation was interrupted; the exact reconstruction conclusion uses the displayed algebraic proof and the verified complete splitting. No expanded Hr coefficient array is claimed in the adopted artifact.',
        variables='Working Cartesian h,p have phase U=diag(1,1,1,i,i,1); physical original h=U h_work and p=U p_work. Constraints C_work=diag(1,1,-i,-i) C_original diag(U,U). In this working basis Q=(h11,(h12/K-h13/M)/2), Pi=(p11,K*p12-M*p13), with h22=h33=h23=0 and h12/K+h13/M=0. K=k0 exp(-xi), M=l0 exp(xi); saved rational coefficients use local K,M. General physical Q before gauge fixing is adjoint(momentum_kernel) h_work.',
        denominator_certificate='d=2-F*(K^2+M^2)/(2t). On every nonzero dual lattice orbit |k*l|>=4*pi^2/sqrt(5), F/t>=e for 1<=t<=2, hence d<=2-4*pi^2*e/sqrt(5)<0 uniformly in xi.',
        scope='Two complex canonical field pairs on the whole real line per complex orbit, and the conjugate orbit for a real field. Infinite-dimensional field reduction, not two total pairs in a finite block. Differential changes of variables use integration by parts on the specified invariant domain; time correction and full lapse/shift reconstruction retained.'))
    print(len(ck.items),'exact reduction checks',flush=True)
