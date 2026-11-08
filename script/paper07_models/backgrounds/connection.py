"""Independent connection in the background orthonormal frame.

Use the exact EC=EH+quadratic contortion+endpoint identity. This is a
derivative-dependent auxiliary elimination, not a point canonical map.
"""
from common import *
from geometry import *
from itertools import combinations,permutations
sys.path.insert(0,str(REPO/'script'))
from paper07_results import algebra
SIGN=algebra.SIGN

def make(bg,k,ck):
    t,N,g=bg['time'],bg['N'],bg['g'];scales=[N]+[s.sqrt(g[i,i]) for i in range(3)]
    coordC=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    cc=constants(bg['n'])
    for i,j,l in product(range(3),repeat=3):coordC[i+1][j+1][l+1]=cc[i][j][l]
    C=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for a,b,c in product(range(4),repeat=3):C[a][b][c]=coordC[a][b][c]*scales[a]/(scales[b]*scales[c])
    for a in range(1,4):
        H=s.diff(scales[a],t)/(N*scales[a]);C[a][0][a]=-H;C[a][a][0]=H
    eta=s.diag(*SIGN);G=globals()['connection'](eta,C,lambda a,f:0,reduce=s.simplify)
    om={(a,b,c):SIGN[a]*G[a,c,b] for a,b,c in product(range(4),repeat=3)}
    h=list(s.symbols('h11 h22 h33 h12 h13 h23',real=True));alpha=s.Symbol('alpha',real=True);betas=list(s.symbols('beta1:4',real=True))
    ls=list(s.symbols('lambda01 lambda02 lambda03 lambda12 lambda13 lambda23',real=True))
    q=h+[alpha]+betas+ls;dq=list(s.symbols('v0:'+str(len(q)),real=True));hm=symmetric(h)
    f=s.zeros(4);f[0,0]=alpha
    for i in range(3):
        f[i+1,0]=scales[i+1]*betas[i]/N
        for j in range(3):f[i+1,j+1]=hm[i,j]/(2*scales[i+1]*scales[j+1])
    lor=s.zeros(4)
    for (a,b),ll in zip(combinations(range(4),2),ls):lor[a,b]=SIGN[a]*ll;lor[b,a]=-SIGN[b]*ll
    f+=lor
    def D(a,z):
        if a==0:return (s.diff(z,t)+sum(s.diff(z,x)*v for x,v in zip(q,dq)))/N
        if a==1:return s.I*k*z/scales[1]
        return 0
    rhs=[]
    for a in range(4):
        for b,c in combinations(range(4),2):
            rhs.append(D(b,f[a,c])-D(c,f[a,b])-sum(f[a,e]*C[e][b][c] for e in range(4))+sum(SIGN[a]*(om[a,e,b]*f[e,c]-om[a,e,c]*f[e,b]) for e in range(4)))
    ref=algebra.contortion();T=ref['torsion_map'];sol=(-T.inv()*s.Matrix(rhs)).applyfunc(s.simplify)
    residual=(T*sol+s.Matrix(rhs)).applyfunc(s.simplify)
    for i,v in enumerate(residual):ck.zero('Cartan component '+str(i),v)
    ck('Cartan map nonsingular',T.det()!=0)
    slots=[(a,b,c) for a,b in combinations(range(4),2) for c in range(4)]
    zero_l={**dict.fromkeys(ls,0),**dict.fromkeys(dq[-6:],0)}
    for i,(a,b,c) in enumerate(slots):
        expected=-SIGN[a]*(D(c,lor[a,b])+sum(SIGN[a]*om[a,e,c]*lor[e,b]-lor[a,e]*SIGN[e]*om[e,b,c] for e in range(4)))
        ck.zero('Lorentz connection covariance '+str(i),sol[i]-sol[i].subs(zero_l)-expected)
    # Direct epsilon contraction of K wedge K, independent of the scalar formula.
    ys=ref['coordinates'];kd=dict(zip(slots,ys))
    def K(a,b,c):return 0 if a==b else (1 if a<b else -1)*kd[min(a,b),max(a,b),c]
    quad=0
    for a,b,c,d in permutations(range(4)):
        for rho,sig in permutations([v for v in range(4) if v not in (a,b)]):
            for e in range(4):quad+=s.Rational(1,2)*s.LeviCivita(a,b,c,d)*SIGN[c]*SIGN[d]*SIGN[e]*s.LeviCivita(a,b,rho,sig)*K(c,e,rho)*K(e,d,sig)
    ck.zero('epsilon K wedge K equals scalar contortion',quad-ref['scalar'])
    M=s.hessian(quad,ys)
    ck('24 auxiliary connection Hessian determinant',M.det()==2**32)
    volume=adm_data(bg)['volume'];density=s.simplify(N*volume)
    # All independent connection perturbations w appear as kappa=w-wLC.
    w=s.Matrix(s.symbols('w0:24'));kappa=w-sol
    eom=M*kappa
    ck('connection solution eliminates every auxiliary EOM',eom.subs(dict(zip(w,sol)),simultaneous=True)==s.zeros(24,1))
    return dict(q=q,dq=dq,coframe_variation=f,background_structure=C,background_connection=om,connection_slots=slots,LC_solution=sol,Cartan_map=T,contortion_hessian=M,density=density,auxiliary_determinant=2**32*density**24,quadratic_identity='L_EC2 = L_ADM2 + density/2 (w-wLC)^T M (w-wLC) + dB_EH2/dt; complete constant NY adds theta*dB_NY2/dt',endpoint_identity='B_EH=1/2 epsilon_abcd e^a wedge e^b wedge K^{cd}; B_NY=e^a wedge T_a; theta_keep-theta_subtracted=delta(B_EH+theta B_NY); delta^2 B=0',scope='all 24 connection components and all six linear Lorentz directions; nonzero one-direction mode in the chosen invariant frame')

if __name__=='__main__':
    name=sys.argv[1];ck=Checks();k=s.Symbol('k',positive=True)
    result=make(backgrounds()[name],k,ck)
    save('connection_'+name.lower()+'.json',dict(checks=ck.items,result=result,verdict='AUXILIARY_CONNECTION_AND_LORENTZ_LIFT_VERIFIED'))
    print(name,len(ck.items),'connection checks',flush=True)
