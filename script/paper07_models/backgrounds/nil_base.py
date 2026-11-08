"""All six metric components on the m_center=0, nonzero base-wave block.

The potential Hessian is computed from the variation of the spatial Einstein
tensor, independently of the Gowdy comparison used only at the end.
"""
from common import *
from geometry import *

def build(ck):
    bg=backgrounds()['Nil'];u=bg['time'];N=bg['N'];g=bg['g'];ad=adm_data(bg);vol=ad['volume'];pi=ad['pi'];inv=g.inv()
    k=s.Symbol('k',positive=True);hs=list(s.symbols('h11 h22 h33 h12 h13 h23',real=True));ps=list(s.symbols('p11 p22 p33 p12 p13 p23',real=True))
    h=symmetric(hs);p=symmetric(ps,True);z=s.Matrix(hs+ps);C=constants(bg['n'])
    def linear_simp(f):
        f=s.expand(f)
        return s.Add(*[simp(f.coeff(v))*v for v in hs]) if any(f.has(v) for v in hs) else simp(f)
    print('full spatial Einstein variation',flush=True)
    R0,_=ricci(g,C,lambda a,f:0)
    print('spatial background complete',flush=True)
    R1,_=linear_ricci(g,h,C,lambda a,f:0,lambda a,f:s.I*k*f if a==0 else 0,reduce=linear_simp)
    print('spatial linear Ricci complete',flush=True)
    scalar0=s.trace(inv*R0);scalar1=linear_simp(s.trace(inv*R1-inv*h*inv*R0));trh=s.trace(inv*h)
    E0=inv*R0*inv-scalar0*inv/2
    dE=inv*R1*inv-inv*h*inv*R0*inv-inv*R0*inv*h*inv-scalar1*inv/2+scalar0*inv*h*inv/2+trh*E0/2
    force=s.Matrix([N*vol*dE[i,j]*(1 if i==j else 2) for i,j in PAIRS])
    V=force.jacobian(hs).applyfunc(simp)
    print('potential Hessian complete',flush=True)
    ck('spatial Hessian self-adjoint on the complete wave block',(V-V.subs(k,-k).T).applyfunc(simp)==s.zeros(6))
    S0=pi*g;S1=p*g+pi*h;S2=p*h
    T0=s.trace(S0*S0)-s.trace(S0)**2/2
    T1=2*s.trace(S0*S1)-s.trace(S0)*s.trace(S1)
    T2=s.trace(S1*S1)+2*s.trace(S0*S2)-s.trace(S1)**2/2-s.trace(S0)*s.trace(S2)
    kin2=N/vol*(T2-trh*T1/2+(trh*trh/8+s.trace(inv*h*inv*h)/4)*T0)
    HH=s.hessian(kin2,z).applyfunc(simp);HH[:6,:6]+=V
    print('kinetic Hessian complete',flush=True)
    HC=N/vol*(T1-trh*T0/2)-N*vol*(scalar1+trh*scalar0/2)
    moms=[]
    for a in range(3):
        deriv=s.I*k if a==0 else 0
        f=deriv*s.trace(pi*h)-2*s.I*k*sum(pi[0,j]*h[a,j]+p[0,j]*g[a,j] for j in range(3))
        f-=sum(pi[i,j]*(C[l][a][i]*h[l,j]+C[l][a][j]*h[i,l])+p[i,j]*(C[l][a][i]*g[l,j]+C[l][a][j]*g[i,l]) for i,j,l in product(range(3),repeat=3))
        moms.append(s.expand(f))
    CM=s.Matrix([HC]+moms).jacobian(z).applyfunc(simp)
    print('constraint rows complete',flush=True)
    J=s.zeros(12);J[:6,6:]=s.eye(6);J[6:,:6]=-s.eye(6);Omega=-J
    ck('four linear constraints first class',(CM*J*CM.subs(k,-k).T).applyfunc(simp)==s.zeros(4))
    # Certify rank and time-dependent propagation without fixing a gauge.
    piv=[6,9,10,7];Cp=CM[:,piv];det=simp(Cp.det())
    ck('constraint rank four for k>0,u>0',det!=0)
    free=[j for j in range(12) if j not in piv]
    null=s.zeros(12,8)
    for a,j in enumerate(free):null[j,a]=1
    solvep=(-Cp.inv()*CM[:,free]).applyfunc(simp)
    for a,j in enumerate(piv):null[j,:]=solvep[a,:]
    ck('constraint solution matrix',(CM*null).applyfunc(simp)==s.zeros(4,8))
    prop=(s.diff(CM,u)+CM*J*HH).applyfunc(simp)
    ck('explicit time derivative included in constraint propagation',(prop*null).applyfunc(simp)==s.zeros(4,8))
    closure=(prop[:,piv]*Cp.inv()).applyfunc(simp)
    ck('closed propagation identity',(prop-closure*CM).applyfunc(simp)==s.zeros(4,12))
    print('constraints propagate; reducing the full block',flush=True)
    A=g[0,0];B=g[2,2];pa0=2*A*pi[0,0]
    gauges=s.Matrix([hs[3],hs[4],hs[1]/A+hs[2]/B,2*A*ps[0]+2*pi[0,0]*hs[0]])
    GM=gauges.jacobian(z)
    faddeev=(GM*J*CM.subs(k,-k).T).applyfunc(simp)
    fd=simp(faddeev.det());ck('four gauge conditions admissible',fd!=0)
    Pdef=hs[2]/B;Qdef=hs[5]/B
    PPdef=-A*ps[1]+B*ps[2]+(A*pi[1,1]+B*pi[2,2])*Pdef
    PQdef=B*ps[5]+2*B*pi[1,1]*Qdef
    physical=s.Matrix([Pdef,Qdef,PPdef,PQdef]).jacobian(z)
    matrix=CM.col_join(GM).col_join(physical)
    target=s.zeros(12,4);target[8:,:]=s.eye(4)
    E=(matrix.inv()*target).applyfunc(simp)
    ck('all constraints eliminated',(CM*E).applyfunc(simp)==s.zeros(4))
    ck('all gauge conditions satisfied',(GM*E).applyfunc(simp)==s.zeros(4))
    ck('physical coordinates normalized',(physical*E).applyfunc(simp)==s.eye(4))
    Jr=s.zeros(4);Jr[:2,2:]=s.eye(2);Jr[2:,:2]=-s.eye(2)
    ck('pulled symplectic form canonical',(E.subs(k,-k).T*Omega*E).applyfunc(simp)==-Jr)
    time_term=(E.subs(k,-k).T*Omega*s.diff(E,u)).applyfunc(simp)
    HR=(E.subs(k,-k).T*HH*E-time_term).applyfunc(simp)
    RR=4*u/3;ee=16*u*u/(9*A*A)
    expected=s.zeros(4)
    expected[0,0]=3*RR*(k*k+2*ee);expected[1,1]=3*RR*ee*k*k
    expected[0,1]=-6*s.I*RR*ee*k;expected[1,0]=6*s.I*RR*ee*k
    expected[2,2]=3/RR;expected[3,3]=3/(RR*ee)
    ck('full ADM reduction equals independently reduced twisted Gowdy Hamiltonian',(HR-expected).applyfunc(simp)==s.zeros(4))
    multipliers=(-faddeev.inv()*(s.diff(GM,u)*E+GM*J*HH*E)).applyfunc(simp)
    dynamic=(J*HH*E+J*CM.subs(k,-k).T*multipliers-s.diff(E,u)-E*Jr*HR).applyfunc(simp)
    ck('all full Hamilton equations reconstructed',dynamic==s.zeros(12,4))
    ck('reduced Hamiltonian Hermitian',(HR-HR.subs(k,-k).T).applyfunc(simp)==s.zeros(4))
    return dict(checks=ck.items,background=bg,canonical_background=ad,variables=z,spatial_Ricci_variation=R1,spatial_scalar_variation=scalar1,quadratic_Hessian=HH,constraints=CM,constraint_minor=det,propagation_matrix=closure,gauge_rows=GM,gauge_determinant=fd,physical_rows=physical,embedding=E,time_dependent_symplectic_correction=time_term,reduced_Hessian=HR,multiplier_reconstruction=multipliers,area_R=RR,exp_2P0=ee,physical_pairs_both_real_quadratures=4,propagating_field_types=2,scope='central index=0; any nonzero integer base wavevector (n1,n2), k=sqrt(n1^2+n2^2); LRS background allows algebraic frame alignment without rotating the fixed lattice',verdict='NIL_FULL_METRIC_NONZERO_BASE_BLOCK_REDUCED')

if __name__=='__main__':
    ck=Checks();result=build(ck);save('nil_base.json',result);print(len(ck.items),'exact matrix identities',flush=True)
