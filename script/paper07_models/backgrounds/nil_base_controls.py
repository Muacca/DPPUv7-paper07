"""Independent controls of the Nil base-wave reduction: full reconstruction of lapse and shifts, the linearized Ricci tensor on the reduced dynamics and the two coupled wave equations. Writes nil_base_controls.json."""
from common import *
from geometry import *

if __name__=='__main__':
    ck=Checks();bg=backgrounds()['Nil'];u=bg['time'];k=s.Symbol('k',positive=True)
    raw=json.loads((OUT/'nil_base.json').read_text(encoding='utf8'));env={'u':u,'k':k}
    mat=lambda key:s.Matrix([[s.sympify(v,locals=env) for v in row] for row in raw[key]])
    E=mat('embedding');HR=mat('reduced_Hessian');lapse_shift=mat('multiplier_reconstruction')
    y=s.Matrix(s.symbols('P Q PiP PiQ',real=True));J=s.zeros(4);J[:2,2:]=s.eye(2);J[2:,:2]=-s.eye(2)
    evolution=J*HR;dt=lambda f:s.diff(f,u)+sum(s.diff(f,y[i])*(evolution*y)[i] for i in range(4))
    def li(f):
        f=s.expand(f)
        return s.Add(*[simp(f.coeff(v))*v for v in y]) if any(f.has(v) for v in y) else simp(f)
    z=E*y;mult=lapse_shift*y
    ck('fiber shifts vanish after full reconstruction',li(mult[2])==0 and li(mult[3])==0)
    ck('base shift vanishes for nonzero mode',li(mult[1])==0)
    ck.zero('lapse agrees with conformal-base reconstruction',li(mult[0]-z[0]/(2*bg['g'][0,0])))
    g4=s.diag(-bg['N']**2,*bg['g'].diagonal());h4=s.zeros(4)
    h4[0,0]=-2*bg['N']**2*mult[0];h4[1:,1:]=symmetric(list(z[:6,0]))
    for i in range(3):h4[0,i+1]=h4[i+1,0]=bg['g'][i,i]*mult[i+1]
    C=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    cc=constants(bg['n'])
    for i,j,l in product(range(3),repeat=3):C[i+1][j+1][l+1]=cc[i][j][l]
    print('independent full four-dimensional Ricci reconstruction',flush=True)
    D0=lambda a,f:s.diff(f,u) if a==0 else 0
    D1=lambda a,f:dt(f) if a==0 else (s.I*k*f if a==1 else 0)
    ric,_=linear_ricci(g4,h4,C,D0,D1,reduce=li)
    for a in range(4):
        for b in range(a,4):ck.zero('full linearized Ricci '+str(a)+str(b),ric[a,b])
    RR=4*u/3;ee=16*u*u/(9*(1+u**4)**2);P0prime=1/u-4*u**3/(1+u**4)
    dy=evolution*y;ddy=dy.applyfunc(dt)
    ck.zero('first coupled wave equation',li(ddy[0]+dy[0]/u+9*(k*k+2*ee)*y[0]-18*s.I*k*ee*y[1]))
    ck.zero('second coupled wave equation',li(ddy[1]+(1/u+2*P0prime)*dy[1]+9*k*k*y[1]+18*s.I*k*y[0]))
    a,b=s.symbols('n1 n2',real=True);norm=s.sqrt(a*a+b*b)
    rot=s.Matrix([[a/norm,b/norm,0],[-b/norm,a/norm,0],[0,0,1]])
    ck('algebraic horizontal frame alignment orthogonal',(rot*rot.T).applyfunc(s.simplify)==s.eye(3))
    ck('alignment preserves Nil structure',(rot*bg['n']*rot.T-bg['n']).applyfunc(s.simplify)==s.zeros(3))
    ck('alignment maps the given covector to k direction',(rot*s.Matrix([a,b,0])-s.Matrix([norm,0,0])).applyfunc(s.simplify)==s.zeros(3,1))
    save('nil_base_controls.json',dict(checks=ck.items,reconstructed_metric=h4,reconstructed_lapse_shift=mult,Ricci_on_reduced_dynamics=ric,field_evolution=evolution,physical_variables=y,KID_necessary_condition='The only global KID is Z=partial_z. For the entire central-invariant sector L_Z g=L_Z pi=0, so the exact integrated momentum map integral pi^{ij} L_Z g_ij and its quadratic Taub coefficient vanish identically. No general sufficiency theorem or explicit nonlinear lift is claimed.',alignment=rot,lattice_statement='n1,n2 remain the original integer lattice labels. This is a tensor/frame basis change for each covector, not a quotient rotation or identification of distinct modes.',verdict='NIL_FOUR_DIMENSIONAL_AND_GLOBAL_CONTROLS_PASSED'))
    print(len(ck.items),'independent controls',flush=True)
