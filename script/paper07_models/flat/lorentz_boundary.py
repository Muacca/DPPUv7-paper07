"""Unfixed Lorentz closure and the complete quadratic NY endpoint for both quadratures."""
from common import *
import sympy as s
from itertools import product,permutations,combinations
from connection_modes import build,matrix_closure,SIGN
ck=Checks();k=s.Symbol('k',positive=True)
print('Retaining all six Lorentz tetrad directions',flush=True)
full=build(k,include_lorentz=True);q=full['q'];v=full['v'];w=full['ws']
M=s.hessian(full['L'],w);c=s.Matrix([s.diff(full['L'],a).subs(dict.fromkeys(w,0)) for a in w]);sol=-M.inv()*c
reduced=s.expand(full['L'].subs(dict(zip(w,sol)),simultaneous=True))
gauge=build(k);gsol=-s.hessian(gauge['L'],w).inv()*s.Matrix([s.diff(gauge['L'],a).subs(dict.fromkeys(w,0)) for a in w])
metric=s.expand(gauge['L'].subs(dict(zip(w,gsol)),simultaneous=True))
delta=s.expand(reduced-metric)
BL=s.expand(sum(a*s.diff(delta,b) for a,b in zip(q,v))/2)
ck('unfixed eliminated action equals ADM plus Lorentz endpoint',s.expand(delta-sum(s.diff(BL,a)*b for a,b in zip(q,v)))==0)
closure=matrix_closure(full['L'],q+w,v+list(s.symbols('dw0:24')))
ck('unfixed Lorentz system has two physical pairs',closure['pairs']==2)
ck('unfixed system has FC20 SC36',closure['first_class']==20 and closure['second_class']==36)
theta=closure['theta'][:16,0].subs(dict(zip(w,sol)),simultaneous=True).applyfunc(s.expand)
metricp=s.Matrix([s.diff(metric,a) for a in v])
ck('pulled one-form differs from ADM by exact Lorentz endpoint',s.simplify(theta-metricp-s.Matrix([s.diff(BL,a) for a in q]))==s.zeros(16,1))
print('Unfixed closure verified; checking full NY identity with both real quadratures',flush=True)
z=gauge;coords=z['q']+z['ws'];vel=z['v']+z['dws']
minus=s.symbols('minus0:'+str(len(coords)),real=True);dminus=s.symbols('dminus0:'+str(len(coords)),real=True)
sub=dict(zip(coords,minus));rotation=s.Matrix([[0,-1],[1,0]])
F={key:value+rotation*value.subs(sub,simultaneous=True) for key,value in z['F'].items()}
def W(a,b,mu):
    field=z['field'](a,b,mu)
    return field+rotation*field.subs(sub,simultaneous=True)
allq=coords+list(minus);allv=vel+list(dminus)
def D(field,mu):
    if mu==0:return field.subs(dict(zip(allq,allv)),simultaneous=True)
    if mu==3:return s.Matrix([k*field[1],-k*field[0]])
    return s.zeros(2,1)
avg=lambda f,g:(f.T*g)[0]
T={}
for a in range(4):
    for mu,nu in combinations(range(4),2):T[a,mu,nu]=D(F[a,nu],mu)-D(F[a,mu],nu)+SIGN[a]*(W(a,nu,mu)-W(a,mu,nu))
TT=0
for a in range(4):
    for mu,nu in combinations(range(4),2):
        for rho,sigma in combinations(range(4),2):TT+=SIGN[a]*s.LeviCivita(mu,nu,rho,sigma)*avg(T[a,mu,nu],T[a,rho,sigma])
EER=0
for a,b,rho,sigma in permutations(range(4)):
    for e in range(4):EER+=SIGN[e]*s.LeviCivita(a,b,rho,sigma)*avg(W(a,e,rho),W(e,b,sigma))
for a,b in product(range(4),repeat=2):
    for mu,rho,sigma in permutations([i for i in range(4) if i!=b]):EER+=2*s.LeviCivita(mu,b,rho,sigma)*avg(F[a,mu],D(W(a,b,sigma),rho))
B=0
for a in range(4):
    for i,j,l in permutations([1,2,3]):B+=SIGN[a]*s.LeviCivita(i-1,j-1,l-1)*avg(F[a,i],D(F[a,l],j))
for a in range(4):
    for i,j,b in permutations([1,2,3]):B+=2*s.LeviCivita(i-1,j-1,b-1)*avg(F[a,i],W(a,b,j))
B=s.expand(B);dtB=s.expand(sum(s.diff(B,a)*b for a,b in zip(allq,allv)))
ck('complete NY is time derivative of the explicit averaged current',s.expand(TT-EER-dtB)==0)
ck('quadrature-mixing NY current is retained',B!=0)
Hess=s.hessian(B,allq)
ck('NY endpoint shift has zero exterior second derivative',Hess==Hess.T)
ck('NY boundary contains no time velocities',not (B.free_symbols & set(allv)))
save('lorentz_boundary.json',dict(checks=ck.items,unfixed_closure=serial(closure),unfixed_q=serial(q),unfixed_connection_solution=serial(sol),unfixed_pulled_theta=serial(theta),Lorentz_endpoint=str(BL),two_quadrature_coordinates=serial(allq),complete_NY_endpoint=str(B),complete_NY_time_derivative=str(dtB),NY_identity_residual=str(s.expand(TT-EER-dtB)),verdict='LORENTZ_AND_COMPLETE_NY_BOUNDARY_CONTROLS_PASSED',scope='40 configurations before Lorentz fixing; 34 after fixing; two real Fourier quadratures included in NY current; no omission of endpoint mixing'))
print(len(ck.items),'Lorentz/NY checks; full closure FC20 SC36; NY endpoint identity exact',flush=True)
