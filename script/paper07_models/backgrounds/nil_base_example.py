"""Finite-time real parity block; compare Lagrangian and Hamiltonian evolution."""
from common import *
import sympy as s
import numpy as np
from scipy.integrate import solve_ivp

if __name__=='__main__':
    ck=Checks();u=s.Symbol('u',positive=True);k=s.Symbol('k',positive=True)
    R=4*u/3;e=16*u*u/(9*(1+u**4)**2)
    V=s.Matrix([[k*k+2*e,-2*e*k],[-2*e*k,e*k*k]])
    ck.zero('potential determinant',V.det()-e*k*k*(k*k-2*e))
    ck.zero('e monotone derivative',s.diff(e,u)-32*u*(1-3*u**4)/(9*(1+u**4)**3))
    ck.zero('e at left endpoint',e.subs(u,1)-s.Rational(4,9))
    # u=1+v proves the derivative has negative numerator throughout u>=1.
    v=s.Symbol('v',nonnegative=True)
    ck('negative derivative numerator on u>=1',all(c<0 for c in s.Poly(s.expand((1-3*u**4).subs(u,1+v)),v).all_coeffs()))
    ck('lattice lower bound for k^2-2e',1-2*s.Rational(4,9)>0)
    raw=json.loads((OUT/'nil_base.json').read_text(encoding='utf8'))
    parse=lambda key:s.Matrix([[s.sympify(a,locals={'u':u,'k':k}) for a in row] for row in raw[key]])
    HR=parse('reduced_Hessian');E=parse('embedding');C=parse('constraints')
    J=s.zeros(4);J[:2,2:]=s.eye(2);J[2:,:2]=-s.eye(2)
    # cos for P, sin for Q: complex coefficient map (p,-i*q,Pi_p,-i*Pi_q).
    parity=s.diag(1,-s.I,1,-s.I)
    Hr=(s.conjugate(parity).T*HR*parity).applyfunc(s.simplify)
    ck('real parity Hamiltonian',all(s.im(a)==0 for a in Hr))
    ck('real parity potential matches action',(Hr[:2,:2]-3*R*V).applyfunc(s.simplify)==s.zeros(2))
    kval=1.0;ts=np.linspace(1,2,201)
    y0=np.array([.2,.07,-.1,.03]) # p,p_u,q,q_u
    en=lambda x:16*x*x/(9*(1+x**4)**2)
    rn=lambda x:4*x/3
    def lag_rhs(x,y):
        ee=en(x);p0prime=1/x-4*x**3/(1+x**4)
        return [y[1],-y[1]/x-9*(kval*kval+2*ee)*y[0]+18*kval*ee*y[2],
                y[3],-(1/x+2*p0prime)*y[3]-9*kval*kval*y[2]+18*kval*y[0]]
    hf=s.lambdify((u,k),J*Hr,'numpy')
    can0=[y0[0],y0[2],rn(1)*y0[1]/3,rn(1)*en(1)*y0[3]/3]
    runs=[];solutions=[]
    for tol in (1e-10,1e-12):
        lag=solve_ivp(lag_rhs,[1,2],y0,t_eval=ts,method='DOP853',rtol=tol,atol=tol/100)
        ham=solve_ivp(lambda x,y:np.asarray(hf(x,kval),dtype=float)@y,[1,2],can0,t_eval=ts,method='DOP853',rtol=tol,atol=tol/100)
        ck('both integrators successful '+str(tol),lag.success and ham.success)
        canonical=np.array([lag.y[0],lag.y[2],rn(ts)*lag.y[1]/3,rn(ts)*en(ts)*lag.y[3]/3])
        err=float(np.max(np.abs(canonical-ham.y)))
        ck('independent formulation agreement '+str(tol),err<100*tol)
        runs.append(dict(rtol=tol,max_formulation_difference=err,lagrangian_evaluations=lag.nfev,hamiltonian_evaluations=ham.nfev))
        solutions.append(canonical)
    ck('refinement improves formulation agreement',runs[1]['max_formulation_difference']<runs[0]['max_formulation_difference'])
    refinement=float(np.max(np.abs(solutions[1]-solutions[0])))
    ck('refinement difference bounded',refinement<1e-8)
    ef=s.lambdify((u,k),E,'numpy');cf=s.lambdify((u,k),C,'numpy')
    residual=[]
    for x,y in zip(ts,solutions[-1].T):
        z=np.asarray(ef(x,kval),dtype=complex)@np.diag([1,-1j,1,-1j])@y
        residual.append(float(np.max(np.abs(np.asarray(cf(x,kval),dtype=complex)@z))))
    ck('all four reconstructed constraints',max(residual)<1e-11)
    save('nil_base_example.json',dict(checks=ck.items,mode=dict(central_index=0,base_labels=[1,0],k=kval,real_basis='P=sqrt(2)*p*cos(x), Q=sqrt(2)*q*sin(x); same parity for conjugate momenta'),interval=[1,2],time=ts.tolist(),canonical_solution=solutions[-1].tolist(),runs=runs,refinement_difference=refinement,max_constraint_residual=max(residual),positivity='For u in [1,2], 0<e<=4/9 and every nonzero lattice wave has k^2>=1. Thus k^2-2e>=1/9, so both kinetic and potential quadratic forms are positive. The Hamiltonian depends explicitly on time; this is not a conserved energy or a long-time stability theorem.',scope='linear finite-interval example; no nonlinear lift',verdict='NIL_FINITE_INTERVAL_EXAMPLE_VERIFIED'))
    print(runs,'constraint residual',max(residual),flush=True)
