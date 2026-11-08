"""Construct and verify a periodic standing wave plus mean expansion to second order."""
from common import *
import sympy as s
from perturbation import ETA,ricci_linear,ricci_quadratic
ck=Checks();t,x,y,z=s.symbols('t x y z',real=True);a,h,eps=s.symbols('a h eps',real=True)
coords=[t,x,y,z]
def clean(f):
    f=s.expand(f).subs(h**2,a**2/24)
    return s.trigsimp(s.expand_trig(s.expand(f)),method='fu')
A1=B1=h*t;P1=a*s.cos(x)*s.cos(t)/2
A2=h*h*t*t/2+a*a*(1-s.cos(2*t)-s.cos(2*x))/32
B2=-h*h*t*t-a*a*s.cos(2*x)*s.cos(2*t)/64
P2=h*a*(s.sin(t)-t*s.cos(t))*s.cos(x)/2
# ds^2 = exp(2A)(-dt^2+dx^2) + exp(2B+2P)dy^2 + exp(2B-2P)dz^2.
H=s.diag(-2*A1,2*A1,2*(B1+P1),2*(B1-P1))
Q=s.diag(-(2*A2+2*A1*A1),2*A2+2*A1*A1,2*(B2+P2)+2*(B1+P1)**2,2*(B2-P2)+2*(B1-P1)**2)
R1=ricci_linear(H,coords).applyfunc(clean)
R2=(ricci_linear(Q,coords)+ricci_quadratic(H,coords)).applyfunc(clean)
for i in range(4):
    for j in range(i,4):
        ck(f'linear spacetime Ricci {i}{j}',R1[i,j]==0)
        ck(f'second-order spacetime Ricci {i}{j}',R2[i,j]==0)
# Independent reduced warped-product equations (all coefficients in eps^2).
warped={
 'area':s.diff(B2,t,2)-s.diff(B2,x,2)+2*h*h,
 'polarization':s.diff(P2,t,2)-s.diff(P2,x,2)+2*h*s.diff(P1,t),
 'tt':-s.diff(A2,t,2)+s.diff(A2,x,2)-2*s.diff(B2,t,2)-2*s.diff(P1,t)**2,
 'xx':s.diff(A2,t,2)-s.diff(A2,x,2)-2*s.diff(B2,x,2)+2*h*h-2*s.diff(P1,x)**2,
 'tx':-2*s.diff(B2,t,x)-2*s.diff(P1,t)*s.diff(P1,x)
}
for name,res in warped.items():ck('warped-product '+name,clean(res)==0)
ck('periodic x metric',all(clean(f.subs(x,x+2*s.pi)-f)==0 for f in [A2,B2,P2]))
ck('has zero initial displacement and velocity',P2.subs(t,0)==0 and s.diff(P2,t).subs(t,0)==0)
ck('linear metric has intended shear-free homogeneous part',clean((H[1,1]+H[2,2]+H[3,3])-6*h*t)==0)
f=2*P1
E=s.simplify(s.integrate(2*(s.diff(f,t)**2+s.diff(f,x)**2),(x,0,2*s.pi))/(2*s.pi))
charge=s.simplify(s.integrate(2*s.diff(f,t)*s.diff(f,x),(x,0,2*s.pi))/(2*s.pi))
ck('wave energy',E==a*a);ck('translation charge',charge==0)
save('standing_wave.json',dict(checks=ck.items,parameters='q0=1; h^2=a^2/24; both signs of h; spatial period 2pi',ansatz='A=eps A1+eps^2 A2; B=eps B1+eps^2 B2; P=eps P1+eps^2 P2',coefficients=serial(dict(A1=A1,B1=B1,P1=P1,A2=A2,B2=B2,P2=P2)),metric_first=serial(H),metric_second=serial(Q),ricci_first=serial(R1),ricci_second=serial(R2),energy=str(E),translation=str(charge),verdict='SECOND_ORDER_LOCAL_VACUUM_SOLUTION_CONSTRUCTED',scope='periodic standing wave; fixed finite time interval; not an exact finite-amplitude or long-time solution'))
print(len(ck.items),'standing-wave checks; all spacetime Ricci coefficients through order two vanish')
