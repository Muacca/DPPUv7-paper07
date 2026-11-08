"""Continuous-time residuals of one Lorenz trajectory (N=1535, W=5.25): analytic first and second derivatives of the accepted seventh-degree DOP853 dense polynomial at twenty off-stage times, substituted into the residual audit. Writes harmonic_dense_residual.json."""
from residuals import *
from scipy.integrate._ivp.rk import DOP853,Dop853DenseOutput
import inspect,scipy

N=1535;width=5.25;tol=3e-14
g=Sine(N,width);raw=json.loads((OUT/'harmonic_operators.json').read_text())
fn=evaluator(*local_ops([raw['B0'],raw['Bt'],raw['Bx']]))
def rhs(t,y):
    H,V=y.reshape(2,N,10);dx,ddx=g.derivatives(H);b0,bt,bx=fn(t,g.K,g.M)
    acc=ddx+np.einsum('qab,qb->qa',b0,H)+np.einsum('qab,qb->qa',bt,V)+np.einsum('qab,qb->qa',bx,dx)
    return np.array([V,acc]).ravel()
targets=list(1.025+.05*np.arange(20));records=[]
class DenseDOP(DOP853):
    def _dense_output_impl(self):
        pp=super()._dense_output_impl()
        while targets and pp.t_old<=targets[0]<=pp.t:
            t=targets.pop(0);x=(t-pp.t_old)/pp.h
            y=np.zeros_like(pp.y_old);d=np.zeros_like(y);dd=np.zeros_like(y)
            for i,f in enumerate(reversed(pp.F)):
                y=y+f
                if i%2==0:
                    dd=dd*x+2*d;d=d*x+y;y=y*x
                else:
                    dd=dd*(1-x)-2*d;d=d*(1-x)-y;y=y*(1-x)
            y+=pp.y_old;d/=pp.h;dd/=pp.h**2
            HV=y.reshape(2,N,10);Hd=d.reshape(2,N,10)[0];Hdd=dd.reshape(2,N,10)[0]
            records.append(dict(t=float(t),H=HV[0],Ht=Hd,Htt=Hdd,step_width=float(pp.h),
                kinematic_L2=float(np.sqrt(g.length/(N+1)*np.sum((Hd-HV[1])**2))),
                polynomial_value_difference=float(np.max(np.abs(y-pp(t))))))
        return pp
times=np.sort(np.r_[np.linspace(1,2,101),targets])
sol=solve_ivp(rhs,[1,2],g.initial().ravel(),method=DenseDOP,t_eval=times,rtol=tol,atol=tol/100)
assert sol.success and len(records)==20
a=Audit(g,Q=2048);rows=[]
for r in records:
    y,row=a.evaluate(r['t'],np.array([r['H'],r['Ht']]),actual=r['Htt'])
    row.update(kinematic_L2=r['kinematic_L2'],polynomial_value_difference=r['polynomial_value_difference'],step_width=r['step_width'])
    rows.append(row);print('Actual dense',r['t'],'wave',row['wave_residual_L2'],'rho_red',row['rho_red'],flush=True)
save('harmonic_dense_residual.json',dict(N=N,width=width,rtol=tol,success=sol.success,nfev=sol.nfev,results=rows,
     source=inspect.getsource(Dop853DenseOutput),scipy_version=scipy.__version__,
     method='First and second analytic derivatives of the actual accepted seventh-degree DOP853 dense polynomial at twenty off-stage times. H, Ht, Htt enter the independent four-dimensional Ricci operator and the reduced fourth-order system (SG-05); no generator substitution for the time derivatives.'))
