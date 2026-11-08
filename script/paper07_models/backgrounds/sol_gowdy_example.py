"""Finite-time Sol Gowdy example: numerical evolution of the two reduced fields compared with the Bessel solutions, and the linear periodic reconstruction condition. Writes sol_gowdy_example.json."""
from common import *
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import jv

if __name__=='__main__':
    ck=Checks();L=float(np.log((3+np.sqrt(5))/2));k=2*np.pi/L;w=np.sqrt(k*k+4)
    ts=np.linspace(1,2,201)
    exact=np.array([jv(0,k*ts),-k*jv(1,k*ts),.5*jv(0,w*ts),-.5*w*jv(1,w*ts)])
    def rhs(t,y):return [y[1],-y[1]/t-k*k*y[0],y[3],-y[3]/t-w*w*y[2]]
    runs=[]
    for tol in (1e-10,1e-12):
        sol=solve_ivp(rhs,[1,2],exact[:,0],t_eval=ts,method='DOP853',rtol=tol,atol=tol/100)
        ck('integration success '+str(tol),sol.success)
        err=float(np.max(np.abs(sol.y-exact)))
        ck('Bessel independent error '+str(tol),err<tol*50)
        runs.append(dict(rtol=tol,max_error=err,steps_evaluations=sol.nfev))
    ck('refinement improves accuracy',runs[1]['max_error']<runs[0]['max_error'])
    # Real unit-mean-square waves. ell is fixed by both linear constraints.
    xs=np.linspace(0,L,257)[:-1]
    periodic=[]
    for ti,ui,dui,vi,dvi in zip(ts,exact[0],exact[1],exact[2],exact[3]):
        p=np.sqrt(2)*ui*np.cos(k*xs);pt=np.sqrt(2)*dui*np.cos(k*xs)
        r=np.sqrt(2)*vi*np.cos(k*xs);rt=np.sqrt(2)*dvi*np.cos(k*xs)
        px=-np.sqrt(2)*k*ui*np.sin(k*xs);rx=-np.sqrt(2)*k*vi*np.sin(k*xs)
        # This is the FIRST-order mean condition, not the nonlinear charge.
        periodic.append(float(abs(np.mean(2*pt))))
    ck('linear periodic reconstruction condition',max(periodic)<1e-12)
    # Bessel J0 equation checked from J0'=-J1 and J1'=(J0-J2)/2.
    rp=k*k*(jv(2,k*ts)-jv(0,k*ts))/2+exact[1]/ts+k*k*exact[0]
    rr=.5*w*w*(jv(2,w*ts)-jv(0,w*ts))/2+exact[3]/ts+w*w*exact[2]
    ck('p Bessel equation residual',np.max(np.abs(rp))<1e-12)
    ck('r Bessel equation residual',np.max(np.abs(rr))<1e-12)
    save('sol_gowdy_example.json',dict(checks=ck.items,L=L,k=k,areal_squared_coefficients=[k*k,k*k+4],interval=[1,2],runs=runs,time=ts.tolist(),solution=exact.tolist(),linear_periodic_residual=max(periodic),scope='linear two-field standing-wave example; no nonlinear family or long-time stability claimed',verdict='SOL_FINITE_INTERVAL_EXAMPLE_VERIFIED'))
    print(runs,flush=True)
