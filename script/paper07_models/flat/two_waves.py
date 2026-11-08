"""Two nonparallel standing waves: all sourced Fourier sectors through second order."""
from common import *
import sympy as s
from functools import lru_cache
from perturbation import ETA,ricci_linear,einstein_quadratic
ck=Checks();t,x,y,z=s.symbols('t x y z',real=True)
a,b,h=s.symbols('a b h',real=True);coords=[t,x,y,z]
def clean(f):return s.simplify(s.expand(f.rewrite(s.exp)).subs(h**2,(a*a+b*b)/24))
f=a*s.cos(x)*s.cos(t);g=b*s.cos(y)*s.cos(t)
H=s.diag(-6*h*t,2*h*t+g,2*h*t+f,2*h*t-f-g)
barH=H-ETA*s.trace(ETA*H)/2
ck('first-order Ricci zero',ricci_linear(H,coords).applyfunc(s.simplify)==s.zeros(4))
ck('first-order harmonic gauge',all(s.simplify(sum(ETA[m,m]*s.diff(barH[m,n],coords[m]) for m in range(4)))==0 for n in range(4)))
print('Constructing quadratic coordinate Einstein source',flush=True)
source=einstein_quadratic(H,coords)
X,Y,U=s.symbols('X Y U',nonzero=True)
def exponent(node):
    v=s.expand(node.args[0]/s.I)
    nx,ny,nt=[s.diff(v,q) for q in (x,y,t)]
    assert all(n.is_Integer for n in [nx,ny,nt]) and s.expand(v-nx*x-ny*y-nt*t)==0,node
    return X**nx*Y**ny*U**nt
def modes(expr):
    expr=s.expand(expr.rewrite(s.exp)).replace(lambda q:q.func==s.exp,exponent)
    out={}
    for term in s.Add.make_args(s.expand(expr)):
        powers=term.as_powers_dict();nx=int(powers.get(X,0));ny=int(powers.get(Y,0));nt=int(powers.get(U,0));m=int(powers.get(t,0))
        coef=s.cancel(term/(X**nx*Y**ny*U**nt*t**m))
        key=(nx,ny);out.setdefault(key,{})
        out[key][nt,m]=out[key].get((nt,m),0)+coef
    return {key:{nm:s.expand(v).subs(h**2,(a*a+b*b)/24) for nm,v in terms.items() if s.expand(v).subs(h**2,(a*a+b*b)/24)!=0} for key,terms in out.items()}
entries={(i,j):modes(source[i,j]) for i in range(4) for j in range(4)}
ks=sorted({k for bank in entries.values() for k,terms in bank.items() if terms})
def src(i,j,k):return sum(c*t**m*s.exp(s.I*w*t) for (w,m),c in entries[i,j].get(k,{}).items())
def CS(lam):
    if lam==0:return s.S.One,t
    w=s.sqrt(lam);return (s.exp(s.I*w*t)+s.exp(-s.I*w*t))/2,(s.exp(s.I*w*t)-s.exp(-s.I*w*t))/(2*s.I*w)
@lru_cache(None)
def response(lam,w,m):
    resonant=lam==w*w
    start=2 if resonant and w==0 else (1 if resonant else 0)
    degree=m+start
    coeff=s.symbols('v0:'+str(degree-start+1))
    P=sum(v*t**j for j,v in zip(range(start,degree+1),coeff))
    poly=s.Poly(s.diff(P,t,2)+2*s.I*w*s.diff(P,t)+(lam-w*w)*P-t**m,t)
    sol=s.solve(poly.all_coeffs(),coeff,dict=True)[0];P=P.subs(sol)
    part=P*s.exp(s.I*w*t);C,S=CS(lam)
    ans=s.expand(part-part.subs(t,0)*C-s.diff(part,t).subs(t,0)*S)
    assert clean(s.diff(ans,t,2)+lam*ans-t**m*s.exp(s.I*w*t))==0
    assert clean(ans.subs(t,0))==0 and clean(s.diff(ans,t).subs(t,0))==0
    return ans
results={};psi_terms={}
print('Solving sourced wave equations in',len(ks),'spatial Fourier sectors',flush=True)
for k in ks:
    lam=k[0]**2+k[1]**2;G=s.Matrix(4,4,lambda i,j:src(i,j,k));kv=[k[0],k[1],0]
    for nu in range(4):ck(f'source conservation k={k} nu={nu}',clean(-s.diff(G[0,nu],t)+s.I*sum(kv[i]*G[i+1,nu] for i in range(3)))==0)
    if lam==0:
        ck('integrated energy source vanishes',clean(G[0,0])==0)
        ck('three integrated momentum sources vanish',all(clean(G[0,i])==0 for i in range(1,4)))
        psi=s.S.Zero
    else:psi=clean(2*G[0,0].subs(t,0)/lam)
    psi_terms[str(k)]=str(psi)
    ck(f'zero initial momentum source k={k}',all(clean(G[0,i].subs(t,0))==0 for i in range(1,4)))
    C,S=CS(lam);B=s.zeros(4)
    for i in range(4):
        for j in range(4):
            val=sum(-2*c*response(lam,w,m) for (w,m),c in entries[i,j].get(k,{}).items())
            if i==j and i>0:val+=psi*C
            if i==0 and j>0:val+=s.I*kv[j-1]*psi*S
            if j==0 and i>0:val+=s.I*kv[i-1]*psi*S
            B[i,j]=s.expand(val)
    for i in range(4):
        for j in range(i,4):ck(f'second-order reduced Einstein equation k={k} {i}{j}',clean(s.diff(B[i,j],t,2)+lam*B[i,j]+2*G[i,j])==0)
    for nu in range(4):ck(f'second-order harmonic gauge k={k} nu={nu}',clean(-s.diff(B[0,nu],t)+s.I*sum(kv[i]*B[i+1,nu] for i in range(3)))==0)
    Q=B-ETA*s.trace(ETA*B)/2
    results[str(k)]=dict(source=serial(G),trace_reversed_second=serial(B),metric_second=serial(Q))
for k in ks:
    nk=(-k[0],-k[1])
    ck(f'reality pairing {k}',all(clean(s.conjugate(s.sympify(results[str(k)]['metric_second'][i][j],locals={'t':t,'a':a,'b':b,'h':h}))-s.sympify(results[str(nk)]['metric_second'][i][j],locals={'t':t,'a':a,'b':b,'h':h}))==0 for i in range(4) for j in range(i,4)))
save('two_waves.json',dict(checks=ck.items,parameters='h^2=(a^2+b^2)/24; q0=1; periods 2pi',metric_first=serial(H),second_order_wavevectors=ks,initial_psi_modes=psi_terms,modes=results,response_basis_count=response.cache_info().currsize,verdict='NONPARALLEL_TWO_WAVE_SECOND_ORDER_SOLUTION_CONSTRUCTED',scope='all quadratic generated sectors retained; exact symbolic residuals; no finite-amplitude or long-time claim'))
print(len(ck.items),'two-wave checks;',len(ks),'Fourier sectors; local equations and harmonic gauge hold through order two',flush=True)
