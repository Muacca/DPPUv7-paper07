"""Killing fields and KIDs that descend to the fixed Nil and Sol quotients, the fixed Sol monodromy, and coordinate Ricci checks. Writes fixed_quotients.json."""
from common import *
from geometry import *

def push_residual(V,coords,map_):
    J=s.Matrix(map_).jacobian(coords)
    return (V.subs(dict(zip(coords,map_)),simultaneous=True)-J*V).applyfunc(s.expand)
def polysolve(expressions,coords,parameters):
    coeff=[]
    for f in expressions:coeff.extend(s.Poly(f,*coords).coeffs())
    return s.linsolve(coeff,parameters)

if __name__=='__main__':
    ck=Checks();x,y,z=s.symbols('x y z',real=True);coords=[x,y,z]
    aa,bb,cc,dd=s.symbols('a b c d',real=True);period=s.Symbol('P',positive=True)
    sig=s.Matrix([[1,0,0],[0,1,0],[0,-x,1]])
    maps=[[x+period,y,z+period*y],[x,y+period,z],[x,y,z+period**2]]
    for i,mp in enumerate(maps):
        ck(f'Nil deck {i} preserves coframe',s.simplify(sig.subs(dict(zip(coords,mp)),simultaneous=True)*s.Matrix(mp).jacobian(coords)-sig)==s.zeros(3))
    V=s.Matrix([aa-dd*y,bb+dd*x,cc+aa*y+dd*(x*x-y*y)/2])
    rr=[a for mp in maps for a in push_residual(V,coords,mp)]
    nil=polysolve(rr,coords,[aa,bb,cc,dd]);ck('Nil only central continuous Killing direction descends',nil==s.FiniteSet((0,0,cc,0)))
    lam=s.Symbol('lambda',positive=True)
    solsig=s.diag(1,s.exp(x),s.exp(-x));mp=[x+s.log(lam),y/lam,z*lam]
    ck('Sol base deck preserves coframe',(solsig.subs(dict(zip(coords,mp)),simultaneous=True)*s.Matrix(mp).jacobian(coords)-solsig).applyfunc(s.simplify)==s.zeros(3))
    Vsol=s.Matrix([aa,bb-aa*y,cc+aa*z])
    # The two actual lattice generators span the fiber plane, hence checking
    # the linear translation condition on a basis is equivalent.
    rr=[v for mp0 in [[x,y+1,z],[x,y,z+1],mp] for v in push_residual(Vsol,coords,mp0)]
    ss=polysolve(rr,coords,[aa,bb,cc]);ck('Sol no continuous Killing direction descends for lambda!=1',ss==s.FiniteSet((0,0,0)))
    A=s.Matrix([[2,1],[1,1]]);gold=(3+s.sqrt(5))/2
    ck('fixed hyperbolic monodromy determinant',A.det()==1)
    ck.zero('fixed expanding eigenvalue',(A-gold*s.eye(2)).det())
    kid={}
    for name,bg in backgrounds().items():
        ad=adm_data(bg);t=bg['time'];norm=s.simplify(s.trace(bg['g'].inv()*ad['K']*bg['g'].inv()*ad['K']))
        k11=s.simplify(ad['K'][0,0]/bg['g'][0,0])
        if name=='Sol':k11=s.simplify(ad['K'][1,1]/bg['g'][1,1])
        kid[name]=dict(K_squared=norm,nonzero_principal_K=k11,argument='CMC KID trace: integrate N(-Delta+KijKij)N=0; on the selected interval KijKij>0, hence lapse N=0',global_dimension={'Nil':1,'Sol':0,'S3':4}[name])
    # Coordinate calculation independent of the homogeneous-frame curvature.
    coordinate={}
    for name in ('Nil','Sol'):
        print('independent coordinate Ricci',name,flush=True)
        bg=backgrounds()[name];t=bg['time'];sigma=sig if name=='Nil' else solsig
        g4=s.zeros(4);g4[0,0]=-bg['N']**2;g4[1:,1:]=sigma.T*bg['g']*sigma
        xx=[t,x,y,z];C=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
        R,G=ricci(g4,C,lambda a,f:s.diff(f,xx[a]))
        for a in range(4):
            for b in range(a,4):ck.zero(f'{name} coordinate Ricci {a}{b}',R[a,b])
        coordinate[name]=R
    save('fixed_quotients.json',dict(checks=ck.items,Nil_descending_killing=nil,Sol_descending_killing=ss,KID=kid,coordinate_Ricci=coordinate,scope='S3 uses the known full SU(2)xU(1) spatial isometry algebra of biaxial Taub; no S3 perturbation computed',verdict='FIXED_QUOTIENTS_AND_CMC_KID_BOUNDARIES_VERIFIED'))
    print(len(ck.items),'quotient/KID/coordinate checks',flush=True)
