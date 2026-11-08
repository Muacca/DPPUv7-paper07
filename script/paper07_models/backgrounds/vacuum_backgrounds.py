"""The three selected vacuum backgrounds on Nil, Sol and the Taub three-sphere: four-dimensional Ricci flatness, ADM data and homogeneous constraint ranks. Writes vacuum_backgrounds.json."""
from common import *
from geometry import *
from itertools import combinations

def homogeneous_rank(bg,pi):
    x=s.symbols('g11 g22 g33 g12 g13 g23');p=s.symbols('p11 p22 p33 p12 p13 p23')
    g=symmetric(x);P=symmetric(p,True);n=bg['n']
    curvature=s.trace(n*g)**2/2-s.trace(n*g*n*g)
    F=s.trace(P*g*P*g)-s.trace(P*g)**2/2-curvature
    As=[s.Matrix(3,3,lambda i,j:sum(n[i,z]*s.LeviCivita(z,k,j) for z in range(3))) for k in range(3)]
    Js=[s.trace(P*(a.T*g+g*a)) for a in As]
    sub={**dict(zip(x,[bg['g'][i,j] for i,j in PAIRS])),**dict(zip(p,[pi[i,j]*(1 if i==j else 2) for i,j in PAIRS]))}
    D=s.Matrix(Js).jacobian(x+p).subs(sub).applyfunc(s.simplify)
    H=s.Matrix([F]+Js).jacobian(x+p).subs(sub).applyfunc(s.simplify)
    r=D.rank();rr=H.rank()
    def witness(M,r):
        for rows in combinations(range(M.rows),r):
            for cols in combinations(range(M.cols),r):
                det=s.factor(M.extract(rows,cols).det())
                if det!=0:return dict(rows=rows,cols=cols,minor=det)
    return dict(momentum_rank=r,constraint_rank=rr,momentum_minor=witness(D,r),constraint_minor=witness(H,rr),F=s.simplify(F.subs(sub)),J=[s.simplify(j.subs(sub)) for j in Js],momentum_jacobian=D,constraint_jacobian=H)

if __name__=='__main__':
    ck=Checks();out={}
    for name,bg in backgrounds().items():
        print('frame Ricci',name,flush=True)
        t=bg['time'];g4=s.diag(-bg['N']**2,*bg['g'].diagonal());c=constants(bg['n'])
        C=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
        for i,j,k in product(range(3),repeat=3):C[i+1][j+1][k+1]=c[i][j][k]
        D=lambda a,f:s.diff(f,t) if a==0 else 0
        R,G=ricci(g4,C,D)
        for i in range(4):
            for j in range(i,4):ck.zero(f'{name} Ricci {i}{j}',R[i,j])
        ck(f'{name} Ricci symmetric',R==R.T)
        ad=adm_data(bg);K=ad['K'];inv=bg['g'].inv()
        R3=(s.trace(bg['n']*bg['g'])**2/2-s.trace(bg['n']*bg['g']*bg['n']*bg['g']))/bg['g'].det()
        ck.zero(f'{name} independent ADM Hamiltonian',R3+ad['traceK']**2-s.trace(inv*K*inv*K))
        hr=homogeneous_rank(bg,ad['pi'])
        ck.zero(f'{name} canonical Hamiltonian',hr['F'])
        for j,v in enumerate(hr['J']):ck.zero(f'{name} homogeneous momentum {j}',v)
        ck(f'{name} ranks', (hr['momentum_rank'],hr['constraint_rank'])==({'Nil':(2,3),'Sol':(3,4),'S3':(2,3)}[name]))
        # Coordinate-frame torsion and metric compatibility, all components.
        for a,b,cidx in product(range(4),repeat=3):
            ck.zero(f'{name} torsion {a}{b}{cidx}',G[cidx,a,b]-G[cidx,b,a]-C[cidx][a][b])
            ck.zero(f'{name} metric compatibility {a}{b}{cidx}',D(a,g4[b,cidx])-sum(G[e,a,b]*g4[e,cidx]+G[e,a,cidx]*g4[b,e] for e in range(4)))
        out[name]=dict(background=bg,ADM=ad,Ricci=R,spatial_scalar=s.factor(R3),homogeneous=hr,connection=G)
    save('vacuum_backgrounds.json',dict(checks=ck.items,results=out,verdict='THREE_VACUUM_BACKGROUNDS_VERIFIED'))
    print(len(ck.items),'background checks',flush=True)
