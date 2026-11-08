"""Time-gauge tetrad + all 24 independent connection components, one real Fourier block."""
from common import *
import sympy as s
from itertools import product,permutations,combinations
sys.path.insert(0,str(REPO/'script'))
from paper07_results import torus,algebra
SIGN=(-1,1,1,1)
def build(k,include_lorentz=False):
    hs=s.symbols('h11 h22 h33 h12 h13 h23',real=True)
    vs=s.symbols('dh11 dh22 dh33 dh12 dh13 dh23',real=True)
    bs=s.symbols('b1:4',real=True);alpha=s.Symbol('lapse',real=True)
    q=[*hs,*bs,alpha];v=[*vs,*s.symbols('db1 db2 db3 da',real=True)]
    hm=torus.symmetric(hs);F={(a,m):s.zeros(2,1) for a,m in product(range(4),repeat=2)}
    F[0,0]=s.Matrix([alpha,0])
    for i in range(3):
        F[i+1,0]=s.Matrix([0,bs[i]])
        for j in range(3):F[i+1,j+1]=s.Matrix([hm[i,j]/2,0])
    if include_lorentz:
        ls=s.symbols('lambda01 lambda02 lambda03 lambda12 lambda13 lambda23',real=True)
        dl=s.symbols('dlambda01 dlambda02 dlambda03 dlambda12 dlambda13 dlambda23',real=True)
        q.extend(ls);v.extend(dl)
        for (a,b),value in zip(combinations(range(4),2),ls):
            lf=s.Matrix([0,value]) if a==0 else s.Matrix([value,0])
            F[a,b]+=SIGN[a]*lf;F[b,a]-=SIGN[b]*lf
    slots=[(a,b,c) for a,b in combinations(range(4),2) for c in range(4)]
    ws=s.symbols('w0:24',real=True);wd=dict(zip(slots,ws))
    def field(a,b,c):
        if a==b:return s.zeros(2,1)
        val=(1 if a<b else -1)*wd[min(a,b),max(a,b),c]
        boost=(a==0 or b==0)
        cosine=(boost and c!=0) or (not boost and c==0)
        return s.Matrix([val,0]) if cosine else s.Matrix([0,val])
    def D(f,r):
        if r==0:return f.subs(dict(zip(q,v)),simultaneous=True)
        if r==3:return s.Matrix([k*f[1],-k*f[0]])
        return s.zeros(2,1)
    avg=lambda f,g:(f.T*g)[0]
    cross=0;quad=0;boundary=0;raw=0
    dws=s.symbols('dw0:24',real=True)
    for a,b,c,d in permutations(range(4)):
        ec=s.LeviCivita(a,b,c,d)*SIGN[c]*SIGN[d]
        for mu,r,sigma in permutations([i for i in range(4) if i!=b]):
            e=ec*s.LeviCivita(mu,b,r,sigma)
            wf=field(c,d,sigma)
            cross-=e*avg(D(F[a,mu],r),wf)
            dw=wf.subs(dict(zip(ws,dws)),simultaneous=True) if r==0 else (s.Matrix([k*wf[1],-k*wf[0]]) if r==3 else s.zeros(2,1))
            raw+=e*avg(F[a,mu],dw)
            if r==0:boundary+=e*avg(F[a,mu],wf)
        for rho,sigma in permutations([i for i in range(4) if i not in (a,b)]):
            for eidx in range(4):quad+=s.Rational(1,2)*ec*SIGN[eidx]*s.LeviCivita(a,b,rho,sigma)*avg(field(c,eidx,rho),field(eidx,d,sigma))
    return dict(q=q,v=v,ws=list(ws),dws=list(dws),slots=slots,F=F,field=field,cross=s.expand(cross),quad=s.expand(quad),boundary=s.expand(boundary),raw=s.expand(raw+quad),L=s.expand(cross+quad))

def matrix_closure(L,q,v):
    n=len(q);W=s.hessian(L,v)
    assert W==s.zeros(n),'This first-order closure requires zero velocity Hessian'
    theta=s.Matrix([s.diff(L,d) for d in v]);B=theta.jacobian(q)
    assert theta==B*s.Matrix(q)
    A=s.hessian(L.subs(dict.fromkeys(v,0)),q)
    P=(-B).row_join(s.eye(n));J=s.zeros(2*n);J[:n,n:]=s.eye(n);J[n:,:n]=-s.eye(n)
    HH=s.zeros(2*n);HH[:n,:n]=-A
    Phi=P;gens=[n]
    for _ in range(8):
        M=Phi*J*P.T
        null=M.T.nullspace()
        rows=s.Matrix.vstack(*[(u.T*Phi*J*HH) for u in null]) if null else s.zeros(0,2*n)
        joint=Phi.col_join(rows);rref,piv=joint.rref();new=rref[:len(piv),:]
        if len(piv)==Phi.rows:break
        gens.append(len(piv)-Phi.rows);Phi=new
    else:raise RuntimeError('closure limit')
    C=Phi*J*Phi.T;sc=C.rank();fc=Phi.rows-sc
    return dict(n=n,primary=n,generations=gens,total=Phi.rows,second_class=sc,first_class=fc,pairs=s.Rational(2*n-sc-2*fc,2),primary_rank=(P*J*P.T).rank(),theta=theta,primary_rows=P,constraint_rows=Phi,poisson=C,H=-s.Rational(1,2)*(s.Matrix(q).T*A*s.Matrix(q))[0])

if __name__=='__main__':
    ck=Checks();k=s.Symbol('k',positive=True)
    print('Constructing first-order EC quadratic action',flush=True)
    z=build(k);L=z['L'];w=z['ws'];q=z['q'];v=z['v']
    M=s.hessian(L,w);C=s.Matrix([s.diff(L,a).subs(dict.fromkeys(w,0)) for a in w])
    ck('24-component connection Hessian nondegenerate',M.det()==2**32)
    ref=algebra.contortion();ck('independent epsilon contraction equals contortion quadratic form',s.expand(z['quad']-ref['scalar'].subs(dict(zip(ref['coordinates'],w))))==0)
    solution=-M.inv()*C
    reduced=s.expand(L.subs(dict(zip(w,solution)),simultaneous=True))
    expected=torus.mode_action(k)['L']
    ck('eliminated first-order action equals full quadratic ADM action',s.expand(reduced-expected)==0)
    dtB=sum(s.diff(z['boundary'],a)*b for a,b in zip(q+w,v+z['dws']))
    ck('raw epsilon action differs by the explicit time endpoint derivative',s.expand(z['raw']-L-dtB)==0)
    # Independent torsion equations from d f^a + omega^a_b wedge e0^b.
    def deriv(f,mu):
        if mu==0:return f.subs(dict(zip(q,v)),simultaneous=True)
        if mu==3:return s.Matrix([k*f[1],-k*f[0]])
        return s.zeros(2,1)
    torsion=[]
    for a in range(4):
        for mu,nu in combinations(range(4),2):
            T=deriv(z['F'][a,nu],mu)-deriv(z['F'][a,mu],nu)+SIGN[a]*(z['field'](a,nu,mu)-z['field'](a,mu,nu))
            residual=T.subs(dict(zip(w,solution)),simultaneous=True).applyfunc(s.expand)
            ck(f'Cartan torsion {a}{mu}{nu}',residual==s.zeros(2,1));torsion.append(residual)
    print('24 Cartan equations agree; closing all primary constraints',flush=True)
    closure=matrix_closure(L,q+w,v+list(s.symbols('vdw0:24')))
    ck('nonzero-mode first-order physical pairs',closure['pairs']==2)
    ck('first-order closure FC8 SC48',closure['first_class']==8 and closure['second_class']==48)
    metric=torus.mode_closure(k)
    ck('same metric pair count',metric['physical_pairs']==closure['pairs'])
    pulltheta=closure['theta'][:10,0].subs(dict(zip(w,solution)),simultaneous=True)
    metricp=s.Matrix([s.diff(reduced,d) for d in v])
    ck('canonical endpoint one-form pulls back to ADM momenta',s.simplify(pulltheta-metricp)==s.zeros(10,1))
    A,B,dA,dB=s.symbols('A B dA dB',real=True)
    sub=dict(zip(q,[A,-A,0,B,0,0,0,0,0,0]));sub.update(dict(zip(v,[dA,-dA,0,dB,0,0,0,0,0,0])))
    TT=s.expand(reduced.subs(sub,simultaneous=True))
    ck('TT kinetic and gradient coefficients',s.expand(TT-(dA*dA+dB*dB-k*k*(A*A+B*B))/2)==0)
    thetaTT=s.Matrix([pulltheta[0]-pulltheta[1],pulltheta[3]]).subs(sub,simultaneous=True).applyfunc(s.simplify)
    ck('TT reduced symplectic form canonical',thetaTT==s.Matrix([dA,dB]))
    save('connection_modes.json',dict(checks=ck.items,normalization='q0=N0=c=1; unit-mean-square real cosine/sine block; k>0 symbolic',q=serial(q),v=serial(v),connection_slots=z['slots'],connection_variables=serial(w),connection_hessian=serial(M),connection_solution=serial(solution),first_order_action=str(L),metric_action=str(reduced),endpoint_EH=str(z['boundary']),pulled_theta=serial(pulltheta),TT_action=str(TT),TT_theta=serial(thetaTT),closure=serial(closure),metric_closure={key:serial(metric[key]) for key in ['velocity_rank','primary_count','generations','first_class','second_class','physical_pairs']},verdict='FIRST_ORDER_AND_METRIC_NONZERO_MODE_REDUCTIONS_IDENTICAL'))
    print(len(ck.items),'Connection checks; 34 configurations,',closure['first_class'],'FC,',closure['second_class'],'SC,',closure['pairs'],'pairs',flush=True)
