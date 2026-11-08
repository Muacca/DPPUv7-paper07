"""Space, window and time comparisons of the reconstructed reduced variables on the common interior Gauss-Legendre quadrature. Writes harmonic_comparisons.json."""
from evolve import *
from scipy.special import roots_legendre

class Fields:
    def __init__(self,name,x):
        self.raw=json.loads((OUT/name).read_text());self.g=Sine(self.raw['N'],self.raw['width'])
        self.solution=np.array(self.raw['solution']);self.P=self.g.interpolators(x,3)
        self.K=self.g.k0*np.exp(-x);self.M=self.g.l0*np.exp(x)
        m=json.loads((OUT/'initial_maps.json').read_text())
        self.yh=evaluator(*local_ops(m['reduced_h']));self.yv=evaluator(*local_ops(m['reduced_v']))
    def fields(self,index):
        t=self.raw['time'][index];H,V=self.solution[index]
        h=[p@dst(H,type=1,axis=0,norm='ortho') for p in self.P]
        v=[p@dst(V,type=1,axis=0,norm='ortho') for p in self.P[:3]]
        y=sum(np.einsum('qab,qb->qa',a,h[j]) for j,a in enumerate(self.yh(t,self.K,self.M)))+sum(np.einsum('qab,qb->qa',a,v[j]) for j,a in enumerate(self.yv(t,self.K,self.M)))
        return h[0],v[0],y

def compare(a,b,Q=2048):
    g=Sine(1,5.25);r,w=roots_legendre(Q);x=g.center+4*r;w=4*w
    aa=Fields(a,x);bb=Fields(b,x);errors=[]
    for i,t in enumerate(aa.raw['time']):
        av=aa.fields(i);bv=bb.fields(i)
        errors.append(dict(t=t,metric=float(np.sqrt(np.sum(w[:,None]*(av[0]-bv[0])**2))),
                           velocity=float(np.sqrt(np.sum(w[:,None]*(av[1]-bv[1])**2))),
                           reduced=float(np.sqrt(np.sum(w[:,None]*(av[2]-bv[2])**2)))))
    return dict(a=a,b=b,core_half_width=4.,quadrature=Q,errors=errors,
                maxima={key:max(x[key] for x in errors) for key in ['metric','velocity','reduced']})

if __name__=='__main__':
    cases=[
      ('harmonic_N1023_W525_fine.json','harmonic_N1535_W525_fine.json'),
      ('harmonic_N1535_W525_fine.json','harmonic_N1663_W56875_fine.json'),
      ('harmonic_N1535_W525.json','harmonic_N1535_W525_fine.json')]
    rows=[]
    for a,b in cases:
        row=compare(a,b);rows.append(row);print(a,b,row['maxima'],flush=True)
    save('harmonic_comparisons.json',rows)
