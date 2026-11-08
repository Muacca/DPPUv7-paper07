"""Finite-time central Nil example (m=1, r=0, L=2): reduced versus full ADM evolution, reconstruction of constraints and gauge, and the central KID condition. Writes finite_example.json."""
from reduce_blocks import *
import numpy as np
from scipy.integrate import solve_ivp

if __name__=='__main__':
    ck=Checks();L=2;u,b,raw,a=load_block(L);allred=json.loads((OUT/'reduction.json').read_text(encoding='utf8'));r=allred['blocks'][str(L)];env={'u':u,'b':b}
    mat=lambda key:parse_matrix(r[key],env)
    E,G,Hr,Jr=map(mat,['embedding','gauge','reduced_H','reduced_J']);H,C,J=a['H'],a['C'],a['J'];nq=H.rows//2;nf=Hr.rows//2
    mu=1/(2*np.pi);bval=np.sqrt(mu);ts=np.linspace(1,2,201)
    fn=lambda M:s.lambdify((u,b),M,'numpy',cse=True)
    ef,hf,cf,gf,rf=[fn(M) for M in [E,H,C,G,Jr*Hr]];jn=np.array(J,dtype=complex)
    y0=np.zeros(2*nf,complex);y0[nf]=.03;y0[nf+1]=-.02;z0=np.asarray(ef(1,bval),dtype=complex)@y0
    def red_rhs(t,y):return np.asarray(rf(t,bval),dtype=complex)@y
    def full_rhs(t,z):
        HH=np.asarray(hf(t,bval),dtype=complex);CC=np.asarray(cf(t,bval),dtype=complex);GG=np.asarray(gf(t,bval),dtype=complex)
        dz=jn@HH@z;lag=np.linalg.solve(GG@jn@CC.conj().T,-GG@dz)
        return dz+jn@CC.conj().T@lag
    results=[];sols=[]
    for tol in [1e-10,1e-12]:
        rr=solve_ivp(red_rhs,[1,2],y0,t_eval=ts,rtol=tol,atol=tol/100,method='DOP853')
        ff=solve_ivp(full_rhs,[1,2],z0,t_eval=ts,rtol=tol,atol=tol/100,method='DOP853')
        ck('both finite-time evolutions '+str(tol),rr.success and ff.success)
        zz=np.array([np.asarray(ef(t,bval),dtype=complex)@y for t,y in zip(ts,rr.y.T)]).T
        error=float(np.max(np.abs(zz-ff.y)));ck('reduced/full evolution agreement '+str(tol),error<200*tol)
        results.append(dict(rtol=tol,max_difference=error,reduced_evaluations=rr.nfev,full_evaluations=ff.nfev));sols.append(rr.y)
    ck('refinement improves agreement',results[1]['max_difference']<results[0]['max_difference'])
    constraints=[];gauges=[];charges=[]
    for t,z in zip(ts,zz.T):
        constraints.append(float(np.max(np.abs(np.asarray(cf(t,bval),dtype=complex)@z))))
        gauges.append(float(np.max(np.abs(np.asarray(gf(t,bval),dtype=complex)@z))))
        charges.append(float(-2*mu*np.imag(np.vdot(z[nq:],z[:nq]))))
    ck('all four constraints reconstructed',max(constraints)<1e-10)
    ck('gauge reconstructed',max(gauges)<1e-12)
    ck('central KID necessary condition maintained',max(abs(q) for q in charges)<1e-10)
    ck('nonzero physical metric develops',float(np.max(np.abs(zz[:nq])))>1e-5)
    counter=np.zeros(2*nf,complex);counter[0]=1;counter[nf]=1j;zbad=np.asarray(ef(1,bval),dtype=complex)@counter
    qbad=float(-2*mu*np.imag(np.vdot(zbad[nq:],zbad[:nq])))
    ck('counterexample still satisfies linear constraints',np.max(np.abs(np.asarray(cf(1,bval),dtype=complex)@zbad))<1e-10)
    ck('counterexample has nonzero central charge',abs(qbad-2*mu)<1e-12 and qbad>0)
    save('finite_example.json',dict(checks=ck.items,L=L,central_index=1,lattice_copy=0,central_frequency=mu,interval=[1,2],time=ts.tolist(),initial_reduced_data_real=y0.real.tolist(),runs=results,reduced_solution_real=sols[-1].real.tolist(),reduced_solution_imag=sols[-1].imag.tolist(),max_constraint_residual=max(constraints),max_gauge_residual=max(gauges),central_charge=charges,max_central_charge=max(abs(q) for q in charges),counterexample=dict(reduced_q=[1,0],reduced_p=['i','0'],charge=qbad,expected='2 mu'),real_field='Include conjugate m=-1,r=0 coefficients and conjugate tensor helicities; this doubles the signed central momentum contribution rather than cancelling it. The reported charge includes both sectors.',scope='Finite-time linear wave satisfying the central KID necessary condition, plus a linearly admissible but quadratically obstructed datum. No nonlinear solution family or long-time stability conclusion.'))
    print(results,'constraints',max(constraints),'Q',max(abs(q) for q in charges),flush=True)
