"""Lorentz covariance of the solved connection and an undivided zero-mode control."""
from common import *
import sympy as s
from itertools import combinations
from connection_modes import build,matrix_closure,SIGN
ck=Checks();k=s.Symbol('k',positive=True);f=build(k,True);g=build(k)
def solve(z):
    w=z['ws'];return -s.hessian(z['L'],w).inv()*s.Matrix([s.diff(z['L'],a).subs(dict.fromkeys(w,0)) for a in w])
fs=solve(f);gs=solve(g)
ld=dict(zip(combinations(range(4),2),f['q'][10:]))
def lam(a,b):
    if a==b:return s.zeros(2,1)
    value=(1 if a<b else -1)*ld[min(a,b),max(a,b)]
    return s.Matrix([0,value]) if (a==0 or b==0) else s.Matrix([value,0])
def D(field,mu):
    if mu==0:return field.subs(dict(zip(f['q'],f['v'])),simultaneous=True)
    if mu==3:return s.Matrix([k*field[1],-k*field[0]])
    return s.zeros(2,1)
for a,b,c in f['slots']:
    wfull=f['field'](a,b,c).subs(dict(zip(f['ws'],fs)),simultaneous=True)
    wfixed=g['field'](a,b,c).subs(dict(zip(g['ws'],gs)),simultaneous=True)
    ck(f'solved connection Lorentz transformation {a}{b}{c}',(wfull-wfixed+D(lam(a,b),c)).applyfunc(s.expand)==s.zeros(2,1))
# k=0 is reconstructed before any symbolic-k divisions in the closure.
controls={}
for include in [False,True]:
    z=build(s.S.Zero,include);cl=matrix_closure(z['L'],z['q']+z['ws'],z['v']+list(s.symbols('dw0:24')))
    ck('zero-mode has six pairs '+str(include),cl['pairs']==6)
    controls['unfixed' if include else 'fixed']={key:serial(cl[key]) for key in ['n','primary','generations','total','second_class','first_class','pairs']}
save('connection_controls.json',dict(checks=ck.items,zero_mode_controls=controls,scope='zero-mode is only a quadratic control; no claim of nonlinear regularity at the static torus'))
print(len(ck.items),'controls; full Lorentz covariance; undivided k=0 yields six quadratic pairs')
