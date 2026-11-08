"""Residual audit of a saved Lorenz trajectory on the interior: wave, Lorenz, Einstein and ADM residuals and the normalized residual rho_red of the reduced fourth-order system. Writes audit_<trajectory>.json."""
from evolve import *
from scipy.special import roots_legendre
from math import comb

def apply(ops,jets):return sum(np.einsum('qab,qb->qa',a,jets[j]) for j,a in enumerate(ops))
def norm(v,w):return float(np.sqrt(np.sum(w[:,None]*v*v)))
def compile_ops(raw):return evaluator(*local_ops(raw))

class Audit:
    def __init__(self,g,core=4.,Q=384):
        self.g=g;self.core=core
        r,w=roots_legendre(Q);self.x=g.center+core*r;self.w=core*w
        self.K=g.k0*np.exp(-self.x);self.M=g.l0*np.exp(self.x)
        self.P=g.interpolators(self.x,7)
        base=json.loads((OUT/'harmonic_operators.json').read_text())
        maps=json.loads((OUT/'initial_maps.json').read_text())
        jets=json.loads((OUT/'coefficient_jets.json').read_text())
        self.wave=compile_ops([base['B0'],base['Bt'],base['Bx']])
        self.lor=compile_ops([base['lorenz']]);self.lort=compile_ops([jets['lorenz_t']])
        self.ricci=compile_ops([base['ricci']])
        self.yh=[compile_ops(a) for a in jets['reduced_h']]
        self.yv=[compile_ops(a) for a in jets['reduced_v']]
        self.yht=compile_ops(jets['reduced_h_t']);self.yvt=compile_ops(jets['reduced_v_t'])
        self.ch=compile_ops(maps['constraint_h']);self.cv=compile_ops(maps['constraint_v'])
        self.hr=compile_ops(json.loads((OUT/'inputs/expanded_Hr.json').read_text())['coefficients'])
    def evaluate(self,t,HV,full=True,cutoff=None,actual=None):
        H,V=HV;g=self.g
        hc=dst(H,type=1,axis=0,norm='ortho');vc=dst(V,type=1,axis=0,norm='ortho')
        # A fixed projection of the computed trajectory. Its derivative is
        # projected from the full trajectory, not recomputed after filtering.
        if cutoff is not None:hc[cutoff:]=0;vc[cutoff:]=0
        hj=[p@hc for p in self.P];vj=[p@vc for p in self.P[:7]]
        yh=[f(t,self.K,self.M) for f in self.yh];yv=[f(t,self.K,self.M) for f in self.yv]
        y=apply(yh[0],hj)+apply(yv[0],vj)
        if not full:return y
        dx,ddx=g.derivatives(H);b0,bt,bx=self.wave(t,g.K,g.M)
        acceleration=ddx+np.einsum('qab,qb->qa',b0,H)+np.einsum('qab,qb->qa',bt,V)+np.einsum('qab,qb->qa',bx,dx)
        if actual is not None:acceleration=actual
        ac=dst(acceleration,type=1,axis=0,norm='ortho');aj=[p@ac for p in self.P[:3]]
        if cutoff is not None:
            ac[cutoff:]=0;aj=[p@ac for p in self.P[:3]]
        b0,bt,bx=self.wave(t,self.K,self.M)
        cont=hj[2]+np.einsum('qab,qb->qa',b0,hj[0])+np.einsum('qab,qb->qa',bt,vj[0])+np.einsum('qab,qb->qa',bx,hj[1])
        wres=aj[0]-cont
        timejet=np.c_[hj[0],vj[0],hj[1],aj[0],vj[1],hj[2]]
        Fg=np.einsum('qab,qb->qa',self.lor(t,self.K,self.M)[0],timejet)
        R=np.einsum('qab,qb->qa',self.ricci(t,self.K,self.M)[0],timejet)
        bgF=np.exp(t*t)/np.sqrt(t);scalar=-R[:,0]/bgF+R[:,4]/bgF+(R[:,5]+R[:,6])/t
        Einstein=R.copy()
        for q,d in [(0,-bgF),(4,bgF),(5,t),(6,t)]:Einstein[:,q]-=d*scalar/2
        constraint=apply(self.ch(t,self.K,self.M),hj)+apply(self.cv(t,self.K,self.M),vj)
        yy=[]
        for d in range(5):
            val=sum(comb(d,r)*(apply(yh[r],hj[d-r:])+apply(yv[r],vj[d-r:])) for r in range(d+1))
            yy.append(val)
        ry=apply(self.hr(t,self.K,self.M),yy);rhs=np.c_[ry[:,2:],-ry[:,:2]]
        dy=apply(self.yht(t,self.K,self.M),hj)+apply(yh[0],vj)+apply(self.yvt(t,self.K,self.M),vj)+apply(yv[0],aj)
        err=norm(dy-rhs,self.w)
        return y,dict(t=float(t),wave_residual_L2=norm(wres,self.w),lorenz_L2=norm(Fg,self.w),Einstein_L2=norm(Einstein,self.w),
                     ADM_constraint_L2=norm(constraint,self.w),reduced_strong_L2=err,rho_red=err/max(1,norm(rhs,self.w)),
                     reduced_rhs_L2=norm(rhs,self.w),reduced_L2=norm(y,self.w),metric_L2=norm(hj[0],self.w))

if __name__=='__main__':
    name=sys.argv[1];Q=int(sys.argv[2]) if len(sys.argv)>2 else 384
    wave=json.loads((OUT/name).read_text());g=Sine(wave['N'],wave['width']);a=Audit(g,Q=Q)
    fields=[];rows=[]
    for i,t in enumerate(wave['time']):
        value=a.evaluate(t,np.array(wave['solution'][i]),full=i%5==0)
        if i%5==0:
            y,row=value;rows.append(row)
            print(name,t,'Lorenz',row['lorenz_L2'],'rho_red',row['rho_red'],flush=True)
        else:y=value
        fields.append(y.tolist())
    suffix='_Q'+str(Q) if len(sys.argv)>2 else ''
    save('audit_'+name.removesuffix('.json')+suffix+'.json',dict(input=name,core_half_width=4.,quadrature=Q,time=wave['time'],
        results=rows,reduced_fields=fields,scope='Interior |xi-center|<=4. For W>=5.25 and 1<=t<=2, continuum characteristic speed 1 leaves at least 0.25 margin to the artificial Dirichlet boundary. Finite-window convergence is separately compared; no global error theorem is asserted.'))
