"""Finite-time Taub example at j=3/2: reduced versus full ADM evolution and its refinement. Writes finite_example.json."""
from kid_controls import *
import numpy as np
from scipy.integrate import solve_ivp

if __name__=='__main__':
    j=s.Rational(3,2);t,a=read_operators(j);raw,r=load_reduction(j);ck=Checks();d=4;nq=24;nf=8
    E,Hr,Jr=r['embedding'],r['reduced_H'],r['reduced_J'];C,H,J=a['C'],a['H'],a['J'];G,Cr=r['gauge'],r['independent_constraints']
    lamb=lambda M:s.lambdify(t,M,'numpy',cse=True)
    ef,hf,cf,gf,crf,rf=[lamb(M) for M in (E,H,C,G,Cr,Jr*Hr)];Jn=np.array(J,dtype=complex)
    y0=np.zeros(16);y0[8]=.03;y0[13]=-.02;ts=np.linspace(.25,.5,201)
    z0=np.array(ef(.25),dtype=complex)@y0
    ck('nonzero reduced initial data',np.linalg.norm(y0)>0)
    ck('initial metric perturbation zero',np.max(np.abs(z0[:nq]))==0)
    def reduced_rhs(x,y):return np.asarray(rf(x),dtype=complex)@y
    def full_rhs(x,z):
        HH=np.asarray(hf(x),dtype=complex);CC=np.asarray(crf(x),dtype=complex);GG=np.asarray(gf(x),dtype=complex)
        dz=Jn@HH@z;mu=np.linalg.solve(GG@Jn@CC.conj().T,-GG@dz)
        return dz+Jn@CC.conj().T@mu
    runs=[];solutions=[]
    for tol in (1e-10,1e-12):
        rr=solve_ivp(reduced_rhs,[.25,.5],y0.astype(complex),t_eval=ts,rtol=tol,atol=tol/100,method='DOP853')
        ff=solve_ivp(full_rhs,[.25,.5],z0,t_eval=ts,rtol=tol,atol=tol/100,method='DOP853')
        ck('both evolutions successful '+str(tol),rr.success and ff.success)
        reconstructed=np.array([np.asarray(ef(x),dtype=complex)@y for x,y in zip(ts,rr.y.T)]).T
        error=float(np.max(np.abs(reconstructed-ff.y)))
        ck('reduced/full ADM agreement '+str(tol),error<tol*200)
        runs.append(dict(rtol=tol,max_difference=error,reduced_evaluations=rr.nfev,full_evaluations=ff.nfev));solutions.append(rr.y)
    ck('refinement improves cross formulation agreement',runs[-1]['max_difference']<runs[0]['max_difference'])
    ck('nontrivial metric develops',np.max(np.abs(reconstructed[:nq]))>1e-5)
    D,S,labels=generators(j);left=[np.array(-v.T,dtype=complex) for v in D];B,_=axial(j);Bn=np.array(B,dtype=complex)
    # Populate m=+j and its conjugate m=-j partner; all other left indices vanish.
    partner=np.zeros((d,d),complex)
    for n,label in enumerate(labels):partner[d-1-n,n]=(-1)**int(j-label)
    conjugate_all=np.kron(np.eye(12),partner)
    constraints=[];gauges=[];charges=[];reality=[]
    for x,z in zip(ts,reconstructed.T):
        zp=conjugate_all@np.conjugate(z)
        for zz in (z,zp):
            constraints.append(float(np.max(np.abs(np.asarray(cf(x),dtype=complex)@zz))))
            gauges.append(float(np.max(np.abs(np.asarray(gf(x),dtype=complex)@zz))))
        h=np.zeros((6,d,d),complex);p=np.zeros_like(h)
        h[:,d-1,:]=z[:nq].reshape(6,d);p[:,d-1,:]=z[nq:].reshape(6,d)
        h[:,0,:]=zp[:nq].reshape(6,d);p[:,0,:]=zp[nq:].reshape(6,d)
        for m0,ml in enumerate(labels):
            for n0,nl in enumerate(labels):reality.append(float(np.max(np.abs(np.conjugate(h[:,m0,n0])-(-1)**int(ml-nl)*h[:,d-1-m0,d-1-n0]))))
        qleft=[float(np.real(sum(np.vdot(p[v],L@h[v]) for v in range(6)))) for L in left]
        qaxis=float(np.real(np.vdot(z[nq:],Bn@z[:nq])+np.vdot(zp[nq:],Bn@zp[:nq])))
        charges.append(qleft+[qaxis])
    ck('all reconstructed constraints',max(constraints)<1e-10)
    ck('gauge including conjugate partner',max(gauges)<1e-12)
    ck('full SU2 real field condition',max(reality)<1e-12)
    ck('all four KID quadratic necessary charges',np.max(np.abs(charges))<1e-10)
    save('finite_example.json',dict(checks=ck.items,j='3/2',interval=['1/4','1/2'],time=ts.tolist(),initial_reduced_data=y0.tolist(),mode_statement='Full normalized Wigner j=3/2 block. Example excites left m=+3/2 and its real-field conjugate m=-3/2; all right indices and metric components produced by the reduction are retained.',runs=runs,reduced_solution_real=solutions[-1].real.tolist(),reduced_solution_imag=solutions[-1].imag.tolist(),max_constraint_residual=max(constraints),max_gauge_residual=max(gauges),max_reality_residual=max(reality),KID_charges=charges,max_KID_charge=float(np.max(np.abs(charges))),refinement_difference=float(np.max(np.abs(solutions[-1]-solutions[0]))),scope='Nontrivial smooth linear finite-time solution with all four quadratic KID necessary conditions. No nonlinear lift or horizon/long-time stability claim.'))
    print(runs,'max constraint',max(constraints),'max KID',float(np.max(np.abs(charges))),flush=True)
