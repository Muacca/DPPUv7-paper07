"""Taub KID controls: Taub-parameter tangents, time displacement as pure linear gauge, and the left and axial charges. Writes kid_controls.json."""
from reduce_block import *

def load_reduction(j):
    t=s.Symbol('tau',real=True);env={'tau':t}
    raw=json.loads((OUT/('reduction_j'+str(j).replace('/','_')+'.json')).read_text(encoding='utf8'))
    mat=lambda key:s.Matrix([[s.sympify(v,locals=env) for v in row] for row in raw[key]])
    return raw,{key:mat(key) for key in ('embedding','physical_selector','gauge','independent_constraints','reduced_J','reduced_H')}

def axial(j):
    D,_,labels=generators(j);d=len(labels);C=constants(2*s.eye(3));q=s.symbols('q0:6');h=symmetric(q)
    lie=s.Matrix([-sum(C[l][2][i]*h[l,k]+C[l][2][k]*h[i,l] for l in range(3)) for i,k in PAIRS])
    S=lie.jacobian(q)
    return s.Matrix(s.kronecker_product(s.eye(6),D[2])+s.kronecker_product(S,s.eye(d))),S

if __name__=='__main__':
    ck=Checks();low={}
    for j in [s.S(0),s.Rational(1,2),s.S(1),s.Rational(3,2)]:
        t,a=read_operators(j);C,H,J=a['C'],a['H'],a['J'];d=int(2*j+1);nq=6*d;field=Field(j,t);norm,mul=field.norm,field.mul
        kernel=C.subs(t,s.Rational(1,3)).H.nullspace()
        for n,v in enumerate(kernel):
            ck(f'j={j} KID kernel valid for all time {n}',norm(C.H*v)==s.zeros(12*d,1))
            ck(f'j={j} KID has zero lapse {n}',v[:d,0]==s.zeros(d,1))
        B,_=axial(j);T=s.diag(B,-B.H)
        ck(f'j={j} axial action symplectic',T.H*J+J*T==s.zeros(12*d))
        ck(f'j={j} full Hamiltonian axial invariant',norm(mul(T.H,H)+mul(H,T))==s.zeros(12*d))
        raw,r=load_reduction(j);E,sel,Hr,Jr=r['embedding'],r['physical_selector'],r['reduced_H'],r['reduced_J']
        # Q=Re(p^dag B h), with the real-field pairing over the full left multiplicity.
        A=s.zeros(2*nq);A[nq:,:nq]=B/2;A[:nq,nq:]=B.H/2
        basis=s.Matrix([[s.sympify(v,locals={'tau':t}) for v in row] for row in raw.get('constant_canonical_basis_transform',s.eye(2*nq).tolist())])
        Ew=red(basis.inv()*E);Aw=red(basis.H*A*basis)
        rf=Field(None if j==s.Rational(3,2) else j,t)
        phase=s.I if j==s.Rational(3,2) else 1
        Aq=red(Aw/phase);Qwork=rf.mul(Ew.H,Aq,Ew);Q=phase*Qwork
        ck(f'j={j} axial charge conserved after reduction',rf.norm(s.diff(Qwork,t)+rf.mul((Jr*Hr).H,Qwork)+rf.mul(Qwork,Jr,Hr))==s.zeros(*Q.shape))
        # Spatial gauge directions do not change the quadratic KID form on C=0.
        if j==s.Rational(3,2):
            weights=s.diag(1,s.sqrt(3),s.sqrt(3),1);rs=s.kronecker_product(s.diag(1,-s.I,1,-s.I),weights.inv())
            Cwork=red(rs*C*basis)
            valid=rf.mul(Ew.H,Aq,J,Cwork.H)==s.zeros(2*int(raw['linear_reduced_pairs_per_left_copy']),4*d)
        else:valid=mul(E.H,A,J,C.H)==s.zeros(2*int(raw['linear_reduced_pairs_per_left_copy']),4*d)
        ck(f'j={j} axial charge annihilates pure gauge on constrained data',valid)
        low[str(j)]=dict(rank=raw['linear_constraint_rank_per_left_copy'],adjoint_kernel=kernel,missing_per_left_copy=len(kernel),missing_full_real=d*len(kernel),linear_reduced_pairs_full_real=raw['linear_reduced_pairs_full_real_j'],axial_generator=B,axial_charge_matrix=Q)
        print('KID controls',j,flush=True)
    # Parameter variations of the exact Taub family are j=0 and are not gauge.
    t,a=read_operators(0);_,r=load_reduction(0);ell,m=s.symbols('ell m',real=True);A=1+t*t;V=4-4*t*t-8*m*t/ell
    g=ell**2*s.diag(4*A,4*A,4*V/A);inv=g.inv();pi=red(2*ell**2*V*(inv*s.diff(g,t)*inv-s.trace(inv*s.diff(g,t))*inv))
    zz=s.Matrix([g[i,k] for i,k in PAIRS]+[pi[i,k]*(1 if i==k else 2) for i,k in PAIRS])
    family=red(zz.jacobian([ell,m]).subs({ell:1,m:0}));C=a['C'];J=a['J'];G=r['gauge'];Cr=r['independent_constraints']
    ck('two exact Taub parameter tangents obey linear constraints',red(C*family)==s.zeros(4,2))
    project=s.eye(12)-J*Cr.H*(G*J*Cr.H).inv()*G
    reduced_family=red(r['physical_selector']*project*family)
    pairing=red(reduced_family.T*(-r['reduced_J'])*reduced_family)
    ck('Taub parameter tangents span one nondegenerate linear pair',pairing.det()!=0)
    bg=zz.subs({ell:1,m:0});ck('time displacement is pure linear gauge',red(project*s.diff(bg,t))==s.zeros(12,1))
    qs=s.Matrix(s.symbols('h11 h22 h33 h12 h13 h23'));ps=s.Matrix(s.symbols('p11 p22 p33 p12 p13 p23'));_,B0=axial(0)
    charge0=s.expand((ps.T*B0*qs)[0])
    save('kid_controls.json',dict(checks=ck.items,blocks=low,Taub_parameter_tangents=family,reduced_parameter_tangents=reduced_family,parameter_symplectic_pairing=pairing,homogeneous_axial_charge=charge0,three_SU2_charges='For the complete real fixed-j coefficient set, Q_a=Re sum_{component,n,m,l} conjugate(p_component,m,n) L_a[m,l] h_component,l,n, where L_a=-transpose(D_a). The invariant frame has zero Lie derivative under these spatial isometries. Their conservation follows from identical independent left copies and the canonical pairing; reality is imposed across m,n.',interpretation='j=0 carries the axial U(1) KID, j=1 carries the three left SU(2) KIDs after restoring the left multiplicity. The 3 and 21 pair counts are linearized quotients at a singular background, not smooth nonlinear physical phase-space dimensions. All four quadratic momentum maps must vanish for a tangent to an exact family. No sufficiency theorem or explicit nonlinear lift is established.'))
    print(len(ck.items),'KID and low label checks',flush=True)
