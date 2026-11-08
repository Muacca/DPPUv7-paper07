"""Lapse/shift-retaining Gowdy ADM reduction and independent 4D linear Ricci."""
from common import *
from geometry import *

if __name__=='__main__':
    ck=Checks();t,x=s.symbols('t x',positive=True);k=s.Symbol('k',positive=True)
    af,rf,pf,qf=[s.Function(n)(x) for n in ('a','R','P','Q')]
    g=s.Matrix([[s.exp(2*af),0,0],[0,rf*s.exp(pf),rf*s.exp(pf)*qf],[0,rf*s.exp(pf)*qf,rf*(s.exp(-pf)+s.exp(pf)*qf*qf)]])
    C0=[[[s.S(0) for _ in range(3)] for _ in range(3)] for _ in range(3)]
    print('spatial scalar from full three-metric',flush=True)
    R3,_=ricci(g,C0,lambda a,f:s.diff(f,x) if a==0 else 0)
    expected=s.exp(-2*af)*(-2*s.diff(rf,x,2)/rf+s.diff(rf,x)**2/(2*rf**2)+2*s.diff(af,x)*s.diff(rf,x)/rf-(s.diff(pf,x)**2+s.exp(2*pf)*s.diff(qf,x)**2)/2)
    ck.zero('unrestricted Gowdy spatial scalar',s.trace(g.inv()*R3)-expected)
    a,R,P,Q,N=s.symbols('a R P Q N',positive=True)
    va,vR,vP,vQ=s.symbols('va vR vP vQ',real=True)
    gs=g.subs({af:a,rf:R,pf:P,qf:Q})
    dg=sum((s.diff(gs,u)*v for u,v in zip([a,R,P,Q],[va,vR,vP,vQ])),s.zeros(3))
    H=gs.inv()*dg/(2*N)
    kin=s.exp(a)*R*N*(s.trace(H*H)-s.trace(H)**2)
    Lkin=s.exp(a)/N*(-2*va*vR-vR*vR/(2*R)+R*(vP*vP+s.exp(2*P)*vQ*vQ)/2)
    ck.zero('full ADM kinetic contraction',kin-Lkin)
    pa,pR,pP,pQ=s.symbols('pa pR pP pQ',real=True)
    velocities={va:N*s.exp(-a)*(-pR/2+pa/(4*R)),vR:-N*s.exp(-a)*pa/2,vP:N*s.exp(-a)*pP/R,vQ:N*s.exp(-a-2*P)*pQ/R}
    for v,m in zip([va,vR,vP,vQ],[pa,pR,pP,pQ]):ck.zero('Legendre '+str(m),s.diff(Lkin,v).subs(velocities)-m)
    kineticH=s.exp(-a)*(-pa*pR/2+pa*pa/(8*R)+(pP*pP+s.exp(-2*P)*pQ*pQ)/(2*R))
    ck.zero('Hamiltonian kinetic Legendre',sum(m*velocities[v] for v,m in zip([va,vR,vP,vQ],[pa,pR,pP,pQ]))-Lkin.subs(velocities)-N*kineticH)
    # Linear constraints, one complex k mode (the -k partner is retained).
    A,rho,p,r,da,dR,dp,dr=s.symbols('A rho p r delta_pa delta_pR delta_pP delta_pr')
    zz=s.Matrix([A,rho,p,r,da,dR,dp,dr]);J=s.zeros(8);J[:4,4:]=s.eye(4);J[4:,:4]=-s.eye(4)
    scalar=dR+(t-1/(4*t))*da+(2-1/(2*t*t)-2*k*k)*rho+2*t*s.I*k*p
    momentum=-2*s.I*k*A+s.I*k*(-2*t-1/(2*t))*rho+2*dp-s.I*k*da
    CC=s.Matrix([scalar,momentum]).jacobian(zz)
    ck('linear scalar/momentum first class', (CC*J*CC.subs(k,-k).T).applyfunc(s.simplify)==s.zeros(2))
    GG=s.Matrix([rho,da]).jacobian(zz);gfp=(GG*J*CC.subs(k,-k).T).applyfunc(s.simplify)
    ck('nonzero-mode gauge fixing invertible',s.simplify(gfp.det())!=0)
    sub={rho:0,da:0,A:dp/(s.I*k),dR:-2*t*s.I*k*p}
    ck.zero('solve scalar constraint',scalar.subs(sub));ck.zero('solve momentum constraint',momentum.subs(sub))
    # Canonical reduction: (rho,delta_pa)=0; delta_pa dA + delta_pR drho=0.
    phys=s.Matrix([p,r,dp,dr]);embedding=zz.subs(sub).jacobian(phys)
    omega=-J
    reducedOmega=(embedding.T*omega*embedding).applyfunc(s.simplify)
    target=s.zeros(4);target[:2,2:]=-s.eye(2);target[2:,:2]=s.eye(2)
    ck('pulled reduced two-form canonical',reducedOmega==target)
    pd,rd=s.symbols('pd rd',real=True)
    L2=t*(pd*pd+rd*rd-k*k*p*p-(k*k+4)*r*r)/2
    H2=(dp*dp+dr*dr)/(2*t)+t*(k*k*p*p+(k*k+4)*r*r)/2
    ck.zero('reduced Legendre transform',dp*pd+dr*rd-L2-H2.subs({dp:t*pd,dr:t*rd}) if False else (dp*pd+dr*rd-L2).subs({pd:dp/t,rd:dr/t})-H2)
    # Verify the apparent cross term from Q=e^{-2x} r is an x derivative.
    rx=s.Symbol('rx');ck.zero('twist cross term is spatial derivative',-(rx-2*r)**2-(-rx*rx-4*r*r)-4*r*rx)
    print('full four-dimensional linearized Ricci',flush=True)
    bg=backgrounds()['Sol'];tt=bg['time'];t=tt
    fp=s.Function('p')(t);fr=s.Function('r')(t);fl=s.Function('ell')(t)
    g4=s.diag(-bg['N']**2,*bg['g'].diagonal());h=s.zeros(4);F=bg['g'][0,0]
    h[0,0]=-F*fl/2;h[1,1]=F*fl/2;h[2,2]=t*fp;h[3,3]=-t*fp;h[2,3]=h[3,2]=t*fr
    C=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    cs=constants(bg['n'])
    for i,j,l in product(range(3),repeat=3):C[i+1][j+1][l+1]=cs[i][j][l]
    D0=lambda a,f:s.diff(f,t) if a==0 else 0
    D1=lambda a,f:s.diff(f,t) if a==0 else (s.I*k*f if a==1 else 0)
    Ric,G1=linear_ricci(g4,h,C,D0,D1)
    rules={s.diff(fp,t,2):-s.diff(fp,t)/t-k*k*fp,s.diff(fr,t,2):-s.diff(fr,t)/t-(k*k+4)*fr,s.diff(fl,t,2):4*s.I*k*(fp+t*s.diff(fp,t)),s.diff(fl,t):4*s.I*k*t*fp,fl:4*t*s.diff(fp,t)/(s.I*k)}
    residual=Ric.subs(rules,simultaneous=True).applyfunc(s.simplify)
    for i in range(4):
        for j in range(i,4):ck.zero(f'Sol full linear Ricci {i}{j}',residual[i,j])
    ck.zero('ell time and space constraints compatible',s.diff(4*t*s.diff(fp,t)/(s.I*k),t).subs(rules)-4*s.I*k*t*fp)
    save('sol_gowdy.json',dict(checks=ck.items,spatial_scalar=expected,ADM_kinetic=Lkin,canonical_H_kinetic=kineticH,linear_constraints=[scalar,momentum],gauge_fixing=['delta R=0','delta pa=0; nonzero modes only'],gauge_bracket=gfp,gauge_determinant=s.factor(gfp.det()),constraint_solution=sub,reduced_action=L2,reduced_H=H2,reduced_Omega=reducedOmega,linear_Ricci=Ric,Ricci_residual=residual,mode_equations=['p_tt+p_t/t+k^2 p=0','r_tt+r_t/t+(k^2+4)r=0'],mode_scope='k=2*pi*n/L, n!=0; one real cosine or sine quadrature has two pairs; both real quadratures have four',normalization='EH coefficient and reference spatial volume factored out; areal time t',verdict='SOL_RESTRICTED_GOWDY_TWO_FIELDS_REDUCED'))
    print(len(ck.items),'symbolic checks',flush=True)
