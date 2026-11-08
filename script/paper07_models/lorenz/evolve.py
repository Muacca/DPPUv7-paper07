"""Sine discretization of the independently derived Lorenz-gauge wave system."""
from numerics import *
from scipy.fft import dst,idst,dct
import time

class Sine:
    def __init__(self,N,width):
        phi=(1+np.sqrt(5))/2;c=np.sqrt(1+phi*phi)
        self.k0=2*np.pi/c;self.l0=phi*self.k0;self.center=.5*np.log(self.k0/self.l0)
        self.N=N;self.width=width;self.a=self.center-width;self.length=2*width
        self.x=self.a+np.arange(1,N+1)*self.length/(N+1)
        self.K=self.k0*np.exp(-self.x);self.M=self.l0*np.exp(self.x)
        self.k=np.arange(1,N+1)*np.pi/self.length
    def derivatives(self,H):
        c=dst(H,type=1,axis=0,norm='ortho')
        second=idst(-self.k[:,None]**2*c,type=1,axis=0,norm='ortho')
        ext=np.zeros((self.N+2,H.shape[1]));ext[1:-1]=self.k[:,None]*c
        first=dct(ext,type=1,axis=0)[1:-1]*np.sqrt(2/(self.N+1))/2
        return first,second
    def interpolators(self,x,order):
        angle=(x-self.a)[:,None]*self.k[None,:]
        return [np.sqrt(2/(self.N+1))*np.sin(angle+j*np.pi/2)*self.k[None,:]**j for j in range(order+1)]
    def initial(self):
        raw=json.loads((OUT/'initial_maps.json').read_text())
        ih=evaluator(*local_ops(raw['initial_h']));iv=evaluator(*local_ops(raw['initial_v']))
        original=np.array(json.loads((OUT/'inputs/wave_initial.json').read_text())['initial']).reshape(4,32)
        r=(self.x-self.center)*np.sqrt(8)
        P=[p*8**(.25+j/2)*np.exp(-r*r/2)[:,None] for j,p in enumerate(hermite_polynomials(r,32,4))]
        jets=[p@original.T for p in P]
        apply=lambda op:sum(np.einsum('qab,qb->qa',a,jets[j]) for j,a in enumerate(op))
        return np.array([apply(ih(1.,self.K,self.M)),apply(iv(1.,self.K,self.M))])

if __name__=='__main__':
    N=int(sys.argv[1]);width=float(sys.argv[2]);name=sys.argv[3]
    tol=float(sys.argv[4]) if len(sys.argv)>4 else 1e-10
    g=Sine(N,width);raw=json.loads((OUT/'harmonic_operators.json').read_text())
    fn=evaluator(*local_ops([raw['B0'],raw['Bt'],raw['Bx']]))
    init=g.initial();started=time.monotonic();reached=[.95]
    def rhs(t,y):
        if t>=reached[0]+.1:reached[0]=t;print(name,round(t,4),'seconds',round(time.monotonic()-started,1),flush=True)
        H,V=y.reshape(2,N,10);dx,ddx=g.derivatives(H);b0,bt,bx=fn(t,g.K,g.M)
        acc=ddx+np.einsum('qab,qb->qa',b0,H)+np.einsum('qab,qb->qa',bt,V)+np.einsum('qab,qb->qa',bx,dx)
        return np.array([V,acc]).ravel()
    sol=solve_ivp(rhs,[1,2],init.ravel(),t_eval=np.linspace(1,2,101),rtol=tol,atol=tol/100,method='DOP853')
    save(name,dict(N=N,width=width,time=sol.t.tolist(),solution=sol.y.T.reshape(-1,2,N,10).tolist(),
                   success=sol.success,nfev=sol.nfev,rtol=tol,atol=tol/100,seconds=time.monotonic()-started,
                   basis='Time-fixed sine basis on finite window, homogeneous Dirichlet metric at endpoints. Orbit coordinate is not made periodic.',
                   initial='Same saved reduced initial profile; exact geometric ADM/Lorenz initial map, pointwise sampled.'))
    print('Done',name,sol.success,sol.nfev,flush=True)
