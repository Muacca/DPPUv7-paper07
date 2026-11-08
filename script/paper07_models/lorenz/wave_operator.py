"""Independent four-dimensional linear Ricci operator and Lorenz hyperbolic reduction."""
from common import *
from geometry import connection,ricci,constants,PAIRS
from linear_geometry import ricci_linear
from numerics import local_ops
import sympy as s
from itertools import product

t,F=s.symbols('t F',positive=True);K,M=s.symbols('K M',real=True,nonzero=True)
tpairs=[(0,0),(0,1),(0,2),(0,3),(1,1),(2,2),(3,3),(1,2),(1,3),(2,3)]
phase=[1,1,s.I,s.I,1,1,1,s.I,s.I,1]
jets=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
zero=s.zeros(1,60)
def red(A):return A.applyfunc(lambda q:s.cancel(s.expand(q)))
def dtime(q):return s.diff(q,t)+(2*t-1/(2*t))*F*s.diff(q,F)
def dspace(q):return -K*s.diff(q,K)+M*s.diff(q,M)
def D0(a,q):return s.cancel(dtime(q)) if a==0 else s.S(0)
def D1(a,A):
    if a in (2,3):return red((s.I*K if a==2 else s.I*M)*A)
    out=A.applyfunc(dtime if a==0 else dspace)
    for j,(dt,dx) in enumerate(jets):
        target=(dt+(a==0),dx+(a==1))
        if target in jets:
            n=jets.index(target);out[:,10*n:10*n+10]+=A[:,10*j:10*j+10]
        elif A[:,10*j:10*j+10]!=s.zeros(A.rows,10):raise ValueError('jet overflow')
    return red(out)
h={(a,b):zero.copy() for a,b in product(range(4),repeat=2)}
for q,(a,b) in enumerate(tpairs):
    row=zero.copy();row[q]=phase[q];h[a,b]=h[b,a]=row
g=s.diag(-F,F,t,t);iv=g.inv()
c3=constants(s.Matrix([[0,0,0],[0,0,1],[0,1,0]]))
c=[[[s.S(0) for k in range(4)] for j in range(4)] for i in range(4)]
for i,j,k in product(range(3),repeat=3):c[i+1][j+1][k+1]=c3[i][j][k]
ck=Checks();R0,G=ricci(g,c,D0);ck('background four-dimensional vacuum',R0==s.zeros(4))
print('Independent linear Ricci',flush=True)
R,_=ricci_linear(g,h,c,D0,D1,reduce=red)
trace=sum((iv[a,a]*h[a,a] for a in range(4)),zero)
lor=[]
for b in range(4):
    row=sum((iv[a,a]*(D1(a,h[a,b])-sum((G[d,a,a]*h[d,b]+G[d,a,b]*h[a,d] for d in range(4)),zero)) for a in range(4)),zero)-D1(b,trace)/2
    lor.append(red(row))
wave=[];ric=[]
for q,(a,b) in enumerate(tpairs):
    sym=(D1(a,lor[b])+D1(b,lor[a])-sum(((G[d,a,b]+G[d,b,a])*lor[d] for d in range(4)),zero))/2
    wave.append(red((R[a,b]-sym)/phase[q]));ric.append(red(R[a,b]/phase[q]))
W=s.Matrix.vstack(*wave);RR=s.Matrix.vstack(*ric)
ck('Lorenz reduced time principal coefficient',W[:,30:40]==s.eye(10)/(2*F))
ck('Lorenz reduced space principal coefficient',W[:,50:60]==-s.eye(10)/(2*F))
ck('no mixed time-space principal derivative',W[:,40:50]==s.zeros(10))
B0=red(-2*F*W[:,:10]);Bt=red(-2*F*W[:,10:20]);Bx=red(-2*F*W[:,20:30])
Lf=red(s.diag(1,1,-s.I,-s.I)*s.Matrix.vstack(*lor))
ck('real working hyperbolic coefficients',not any(a.has(s.I) for a in [B0,Bt,Bx,Lf,RR]))
print('Geometric ADM recovery',flush=True)
g3=s.diag(F,t,t);gi=g3.inv();G3=connection(g3,c3,lambda a,q:0)
hh=s.Matrix(3,3,lambda a,b:h[a+1,b+1])
# Matrix-of-row objects is handled explicitly to avoid SymPy block coercion.
hs={(a,b):h[a+1,b+1] for a,b in product(range(3),repeat=2)}
shift=[h[0,a+1] for a in range(3)];lapse=-h[0,0]/(2*F)
W0=s.diag(F*(t-1/(4*t)),s.Rational(1,2),s.Rational(1,2))
trW=s.trace(gi*W0)
trh=sum((gi[a,a]*hs[a,a] for a in range(3)),zero)
dW={}
for a,b in product(range(3),repeat=2):
    shiftder=D1(a+1,shift[b])+D1(b+1,shift[a])-sum(((G3[d,a,b]+G3[d,b,a])*shift[d] for d in range(3)),zero)
    dW[a,b]=red((D1(0,hs[a,b])-shiftder)/2-lapse*W0[a,b])
dtr=sum((gi[a,a]*dW[a,a]-gi[a,a]**2*W0[a,a]*hs[a,a] for a in range(3)),zero)
zrows=[h[a+1,b+1]/phase[tpairs.index((a+1,b+1))] for a,b in PAIRS]
for a,b in PAIRS:
    pi=t*(trh*(gi[a,a]*W0[a,b]*gi[b,b]-trW*gi[a,b])/2+
          gi[a,a]*dW[a,b]*gi[b,b]-
          gi[a,a]*hs[a,b]*gi[b,b]*(W0[a,a]*gi[a,a]+W0[b,b]*gi[b,b])-
          dtr*gi[a,b]+trW*gi[a,a]*hs[a,b]*gi[b,b])
    pi=red(pi*(1 if a==b else 2)/phase[tpairs.index((a+1,b+1))])
    zrows.append(pi)
Z=red(s.Matrix.vstack(*zrows))
ck('ADM recovery first derivatives only',Z[:,30:]==s.zeros(12,30))
ck('ADM recovery real',not Z.has(s.I))
save('harmonic_operators.json',dict(checks=ck.items,pairs=tpairs,phases=phase,jets=jets,B0=B0,Bt=Bt,Bx=Bx,
       lorenz=Lf,ricci=RR,adm=Z,metric=g,
       convention='Working h phases (1,1,i,i,1,1,1,i,i,1); Htt=Hxx+Bx Hx+Bt Ht+B0 H. Lorenz F_b=nabla^a h_ab-(1/2)nabla_b trace(h). Built by direct Ricci variation minus symmetrized covariant derivative of F, independently of the reduced Hamiltonian.',
       boundary='No boundary imposed in the symbolic operator. Finite-window numerical Dirichlet approximation and an interior domain-of-dependence check are recorded separately.'))
print('Saved hyperbolic operators',flush=True)
