"""Independent controls of the Sol Gowdy reduction: background and linear constraints, the quadratic reduced Hamiltonian and the exact expansion of the periodicity charge. Writes sol_gowdy_controls.json."""
from common import *
import sympy as s

if __name__=='__main__':
    ck=Checks();t,k=s.symbols('t k',positive=True);ep=s.Symbol('epsilon');I=s.I
    A,rho,p,r,da,dR,dp,dr,px,rx=s.symbols('A rho p r delta_pa delta_pR delta_pP delta_pr px rx')
    R=t+ep*rho;pa=-2+ep*da;pR=-2*t-1/(2*t)+ep*dR;pP=ep*dp
    Dx=-pa*pR/2+pa*pa/(8*R)+(pP*pP+s.exp(-2*ep*p)*(ep*dr)**2)/(2*R)-2*ep*k*k*rho-(ep*I*k*rho)**2/(2*R)-2*(ep*I*k*A)*(ep*I*k*rho)+R*((2+ep*I*k*p)**2+s.exp(2*ep*p)*(ep*(I*k-2)*r)**2)/2
    C=pa*ep*I*k*A+pR*ep*I*k*rho+pP*(2+ep*I*k*p)+ep*dr*ep*(I*k-2)*r-ep*I*k*da
    ck.zero('full scalar background',Dx.subs(ep,0))
    expected=dR+(t-1/(4*t))*da+(2-1/(2*t*t)-2*k*k)*rho+2*t*I*k*p
    ck.zero('linear scalar independently differentiated',s.diff(Dx,ep).subs(ep,0)-expected)
    ck.zero('linear momentum independently differentiated',s.diff(C,ep).subs(ep,0)-(-2*I*k*A+I*k*(-2*t-1/(2*t))*rho+2*dp-I*k*da))
    Hred=1/(2*t)+((ep*dp)**2+s.exp(-2*ep*p)*(ep*dr)**2)/(2*t)+t*((2+ep*px)**2+s.exp(2*ep*p)*(ep*(rx-2*r))**2)/2
    H2=s.diff(Hred,ep,2).subs(ep,0)/2
    expectedH=(dp*dp+dr*dr)/(2*t)+t*(px*px+rx*rx+4*r*r)/2
    ck.zero('quadratic reduced H plus spatial endpoint',H2-expectedH+2*t*r*rx)
    # Exact periodicity constraint before expanding: mean(PiP P_x+PiQ Q_x)=0.
    charge=s.expand(ep*dp*(2+ep*px)+ep*dr*ep*(rx-2*r))
    ck.zero('first order charge',charge.coeff(ep,1)-2*dp)
    ck.zero('second order source includes twist term',charge.coeff(ep,2)-(dp*px+dr*(rx-2*r)))
    save('sol_gowdy_controls.json',dict(checks=ck.items,full_scalar=Dx,full_momentum=C,quadratic_reduced_H_before_spatial_boundary=H2,exact_charge_expansion=charge,explanation='Second-order mean PiP is available to enforce the periodic constraint; setting all second-order homogeneous corrections to zero is an extra restriction. Sol has no global KID, so this is not the flat-torus four-charge obstruction.'))
    print(len(ck.items),'reduction controls',flush=True)
