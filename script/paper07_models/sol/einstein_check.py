"""Independent four-dimensional Einstein identities for the Sol orbit operators: spatial equations and the Hamiltonian and momentum constraints. Writes einstein_check.json."""
from operators import *

def load_data():
    a=json.loads((OUT/'operators_exact.json').read_text(encoding='utf8'))['data']
    return parse_ops(a['H']),parse_ops(a['C']),parse_ops([a['J']])[0]

if __name__=='__main__':
    ck=Checks();H,C,J=load_data();g,pi,nv,nov=background()
    nt=24;njet=6;ga=compose([J],adjoint(C));A=compose([J],H)
    evo=[]
    for order in range(max(len(A),len(ga))):
        M=s.zeros(nt)
        if order<len(A):M[:12,:12]=A[order]
        if order<len(ga):M[:12,12:16]=ga[order]
        if order==0:M[12:16,16:20]=s.eye(4);M[16:20,20:24]=s.eye(4)
        evo.append(M)
    zz=s.zeros(1,nt*njet)
    h={(a,b):zz.copy() for a,b in product(range(4),repeat=2)}
    for q,(a,b) in enumerate(PAIRS):
        M=zz.copy();M[0,q]=1;h[a+1,b+1]=h[b+1,a+1]=M
    h[0,0][0,12]=-2*F
    for a in range(3):
        M=zz.copy();M[0,13+a]=g[a,a];h[0,a+1]=h[a+1,0]=M
    g4=s.diag(-F,F,t,t);cs=constants(s.Matrix([[0,0,0],[0,0,1],[0,1,0]]))
    c4=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for a,b,c in product(range(3),repeat=3):c4[a+1][b+1][c+1]=cs[a][b][c]
    def D0(a,f):return s.cancel(s.diff(f,t)+(2*t-1/(2*t))*F*s.diff(f,F)) if a==0 else 0
    def D1(a,M):
        if a==0:return red(dt(M)+asjets(compose(fromjets(M,nt),evo),nt,njet))
        if a==1:return jet_dx(M,nt)
        return red((s.I*k/w if a==2 else s.I*l*w)*M)
    R0,_=ricci(g4,c4,D0)
    ck('background Ricci in rational F representation',R0==s.zeros(4))
    print('Four-dimensional full-orbit Ricci variation',flush=True)
    R,_=ricci_linear(g4,h,c4,D0,D1)
    scalar=red(sum((g4.inv()[a,b]*R[b,a] for a,b in product(range(4),repeat=2)),zz))
    G={(a,b):red(R[a,b]-g4[a,b]*scalar/2) for a,b in product(range(4),repeat=2)}
    for a,b in PAIRS:ck('spatial Einstein '+str((a,b)),G[a+1,b+1]==zz)
    CC=[]
    for M in C:
        B=s.zeros(4,nt);B[:,:12]=M;CC.append(B)
    CJ=asjets(CC,nt,njet)
    ck('Einstein Hamiltonian constraint identity',red(G[0,0]+CJ[0,:]/(2*t))==zz)
    for a in range(3):ck('Einstein momentum constraint identity '+str(a),red(G[0,a+1]+CJ[a+1,:]/(2*t))==zz)
    save('einstein_check.json',dict(checks=ck.items,scope='Exact differential identities for all xi and arbitrary nonzero real k,l, before spatial discretization. All metric components and lapse/shift time jets included. ADM evolution implies G_ij=0 and G_00=-C_H/(2t), G_0i=-C_i/(2t); all ten components vanish on the full constraints.'))
    print(len(ck.items),'four-dimensional checks',flush=True)
