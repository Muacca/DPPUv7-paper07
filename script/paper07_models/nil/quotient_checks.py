"""Numerical controls of the Nil quotient basis: deck reindexing, reality, ladder relations and Haar normalization. Writes quotient_checks.json."""
from common import *
import numpy as np
from scipy.special import eval_hermite,gammaln
from scipy.integrate import quad

P=2*np.pi
def hermite(n,x,mu):
    if n<0:return np.zeros_like(x,dtype=float)
    z=np.sqrt(abs(mu))*np.asarray(x)
    return (abs(mu)/np.pi)**.25*np.exp(-z*z/2-(n*np.log(2)+gammaln(n+1))/2)*eval_hermite(n,z)

def mode(m,r,n,x,y,z,cut=10,derivative=None):
    mu=m/P;k=np.arange(-cut,cut+1);q=r+k*m;xi=x+q/mu;f=hermite(n,xi,mu)
    if derivative=='X':f=np.sqrt(abs(mu)/2)*(np.sqrt(n)*hermite(n-1,xi,mu)-np.sqrt(n+1)*hermite(n+1,xi,mu))
    if derivative=='Y':f=1j*mu*xi*f
    if derivative=='Z':f=1j*mu*f
    return np.sqrt(P)*np.exp(1j*mu*z)*np.sum(np.exp(1j*q*y)*f)

if __name__=='__main__':
    ck=Checks();errors={};normal=[]
    for m in [-2,-1,1,2]:
        for r in range(abs(m)):
            for n in [0,1,2,4]:
                mu=m/P;key=f'{m},{r},{n}';ee=[]
                for x,y,z in [(.17,.31,.7),(.8*P,.4*P,.6*P*P)]:
                    f=mode(m,r,n,x,y,z)
                    ee.extend([abs(mode(m,r,n,x+P,y,z+P*y)-f),abs(mode(m,r,n,x,y+P,z)-f),abs(mode(m,r,n,x,y,z+P*P)-f),abs(np.conj(f)-mode(-m,(-r)%abs(m),n,x,y,z)),abs(f-mode(m,r,n,x,y,z,cut=14))])
                    xp=mode(m,r,n,x,y,z,derivative='X');yp=mode(m,r,n,x,y,z,derivative='Y')
                    if m>0:
                        ee.append(abs(xp+1j*yp+np.sqrt(2*mu*(n+1))*mode(m,r,n+1,x,y,z)))
                        ee.append(abs(xp-1j*yp-(np.sqrt(2*mu*n)*mode(m,r,n-1,x,y,z) if n else 0)))
                errors[key]=float(max(ee));ck('deck/reality/ladder and tail '+key,max(ee)<1e-11)
                value=quad(lambda x:sum(hermite(n,x+(r+k*m)/mu,mu)**2 for k in range(-10,11)),0,P,epsabs=1e-12,epsrel=1e-12)[0]
                normal.append(float(abs(value-1)));ck('normalized Haar mode '+key,abs(value-1)<1e-11)
    save('quotient_checks.json',dict(checks=ck.items,max_mode_identity_error=max(errors.values()),max_normalization_error=max(normal),mode_errors=errors,scope='Numerical control of the analytic deck/reindex, conjugation, ladder and normalization derivations. Summation cut 10 compared with 14. Not a numerical proof of all oscillator levels or all central indices.'))
    print(len(ck.items),'quotient numerical controls',max(errors.values()),max(normal),flush=True)
