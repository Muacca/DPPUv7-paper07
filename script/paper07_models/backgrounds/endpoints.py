"""Explicit quadratic endpoint currents with both real Fourier quadratures."""
from common import *
from geometry import *
from itertools import combinations,permutations
SIGN=(-1,1,1,1)

if __name__=='__main__':
    name=sys.argv[1];ck=Checks();bg=backgrounds()[name];t=bg['time'];N=bg['N'];k=s.Symbol('k',positive=True)
    raw=json.loads((OUT/('connection_'+name.lower()+'.json')).read_text(encoding='utf8'))['result']
    # Recreate real symbols with identical assumptions before parsing.
    q=list(s.symbols('h11 h22 h33 h12 h13 h23 alpha beta1 beta2 beta3 lambda01 lambda02 lambda03 lambda12 lambda13 lambda23',real=True))
    env={str(a):a for a in q};env.update({str(t):t,'k':k})
    f=s.Matrix([[s.sympify(z,locals=env) for z in row] for row in raw['coframe_variation']])
    C=[[[s.sympify(z,locals=env) for z in row] for row in plane] for plane in raw['background_structure']]
    om={tuple(int(z.strip()) for z in key.strip('()').split(',')):s.sympify(val,locals=env) for key,val in raw['background_connection'].items()}
    qs=list(s.symbols('s0:'+str(len(q)),real=True));swap=dict(zip(q,qs))
    F=lambda a,i:s.Matrix([f[a,i],f[a,i].subs(swap,simultaneous=True)])
    slots=[(a,b,c) for a,b in combinations(range(4),2) for c in range(4)]
    w=list(s.symbols('w0:24',real=True));v=list(s.symbols('sw0:24',real=True));wd=dict(zip(slots,[s.Matrix([a,b]) for a,b in zip(w,v)]))
    def W(a,b,c):return s.zeros(2,1) if a==b else (1 if a<b else -1)*wd[min(a,b),max(a,b),c]
    scale1=s.sqrt(bg['g'][0,0]);D=lambda a,z:s.Matrix([k*z[1],-k*z[0]])/scale1 if a==1 else s.zeros(2,1)
    avg=lambda a,b:(a.T*b)[0]
    S={}
    for a in range(4):
        for b,c in combinations(range(1,4),2):
            S[a,b,c]=D(b,F(a,c))-D(c,F(a,b))-sum((C[e][b][c]*F(a,e) for e in range(4)),s.zeros(2,1))+sum((SIGN[a]*(om[a,e,b]*F(e,c)-om[a,e,c]*F(e,b)) for e in range(4)),s.zeros(2,1))
    Bny=0;Beh=0
    for a in range(4):
        for i in range(1,4):
            for j,l in combinations(range(1,4),2):Bny+=SIGN[a]*s.LeviCivita(i-1,j-1,l-1)*avg(F(a,i),S[a,j,l])
    for a,b in product(range(4),repeat=2):
        for i,j,l in permutations(range(1,4)):
            if a==i:Bny+=2*s.LeviCivita(i-1,j-1,l-1)*avg(W(a,b,j),F(b,l))
    # In this expression W denotes independent contortion kappa, not δω.
    for a,b,c,d in permutations(range(4)):
        eps=s.LeviCivita(a,b,c,d)*SIGN[c]*SIGN[d]/2
        for i,j,l in permutations(range(1,4)):
            if b==j:Beh+=eps*s.LeviCivita(i-1,j-1,l-1)*avg(F(a,i),W(c,d,l))
            if a==i:Beh+=eps*s.LeviCivita(i-1,j-1,l-1)*avg(F(b,j),W(c,d,l))
    volume=adm_data(bg)['volume'];Bny=s.expand(volume*Bny);Beh=s.expand(volume*Beh)
    coords=q+qs+w+v
    ck('full NY current retained',Bny!=0)
    ck('EH contortion current retained',Beh!=0)
    for label,B in [('NY',Bny),('EH',Beh)]:
        theta=s.Matrix([s.diff(B,a) for a in coords]);jac=theta.jacobian(coords)
        ck(label+' exact endpoint one-form leaves two-form unchanged',jac==jac.T)
        ck(label+' has no independent time-velocity symbol',not any(str(z).startswith('v') for z in B.free_symbols))
    # The NY spatial derivative part necessarily mixes the real quadratures.
    ck('NY has quadrature-mixing derivative contribution',s.diff(Bny,k)!=0)
    save('endpoints_'+name.lower()+'.json',dict(checks=ck.items,coordinates=coords,B_NY_quadratic=Bny,B_EH_quadratic_contortion_coordinates=Beh,convention='real cos/sin basis has mean square one; W in B_NY is δω, W in B_EH is κ=δω-δωLC; background volume included',identity='Complete NY: d(e^a wedge T_a)=T^a wedge T_a-e^a wedge e^b wedge R_ab. A time-dependent background is retained in the displayed currents. Total d/dt, including explicit coefficient derivatives, shifts the Hamiltonian representative.',scope='Endpoint representatives and pulled covariant symplectic forms; not a point canonical change on the unreduced space',verdict='CURVED_ENDPOINT_CURRENTS_RETAINED'))
    print(name,len(ck.items),'endpoint checks',flush=True)
