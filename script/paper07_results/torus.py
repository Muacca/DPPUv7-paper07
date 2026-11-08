"""Flat-torus quadratic ADM action and second-order compact constraints.

For mode counts set q0=N0=c=1. The real Fourier basis has unit mean square.
Metric coefficients multiply cos, shift coefficients multiply sin, k is along z.
"""
import sympy as s
from .homogeneous import symmetric
from .constraints import quadratic_dirac
def spatial_density(dh):
    """Second coefficient of sqrt(g) R, modulo a spatial divergence.

Computed from the two Christoffel products, not a spectral polynomial.
Only the derivative in direction 2 is nonzero.
"""
    d=lambda i,j,k:dh[i,j] if k==2 else 0
    G=lambda a,i,j:(d(a,j,i)+d(a,i,j)-d(i,j,a))/2
    return s.expand(sum(G(k,i,l)*G(l,i,k)-G(k,i,i)*G(l,k,l) for i in range(3) for k in range(3) for l in range(3)))
def mode_action(k):
    h=s.symbols('h11 h22 h33 h12 h13 h23',real=True);dh=s.symbols('dh11 dh22 dh33 dh12 dh13 dh23',real=True)
    B=s.symbols('b1:4',real=True);a=s.Symbol('lapse',real=True)
    hm=symmetric(h);v=symmetric(dh)
    for i in range(3):v[i,2]-=k*B[i];v[2,i]-=k*B[i]
    kinetic=(s.trace(v*v)-s.trace(v)**2)/4
    potential=spatial_density(-k*hm)
    L=s.expand(kinetic+potential+a*k*k*(h[0]+h[1]))
    return dict(L=L,h=h,dh=dh,shift=B,lapse=a)
def mode_closure(k):
    m=mode_action(k)
    return quadratic_dirac(m['L'],[*m['h'],*m['shift'],m['lapse']],[*m['dh'],*s.symbols('db1 db2 db3 da')])
def transverse_action():
    k=s.Symbol('k',positive=True);m=mode_action(k);A,B,dA,dB=s.symbols('A B dA dB',real=True)
    sub=dict(zip(m['h'],[A,-A,0,B,0,0]));sub.update(dict(zip(m['dh'],[dA,-dA,0,dB,0,0])));sub.update({b:0 for b in m['shift']});sub[m['lapse']]=0
    return s.factor(m['L'].subs(sub)),s.hessian(m['L'].subs(sub),(dA,dB))
def compact_benchmarks():
    x,t=s.symbols('x t',real=True);amp,k=s.symbols('a k',positive=True)
    # One polarization has norm squared 2; spatial average cos^2=sin^2=1/2.
    energy_standing=s.simplify(2*amp**2*k**2*(s.sin(k*t)**2+s.cos(k*t)**2)/2)
    energy_travelling=2*amp**2*k**2
    momentum_travelling=-amp**2*k**2
    # Six h^2 = tr Sigma^2 + (1/4)<dot gamma^2 + (grad gamma)^2>.
    return dict(standing_energy=energy_standing,standing_momentum=0,travelling_energy=energy_travelling,travelling_momentum=momentum_travelling,
                standing_minimum_expansion_squared=s.factor(energy_standing/24))
