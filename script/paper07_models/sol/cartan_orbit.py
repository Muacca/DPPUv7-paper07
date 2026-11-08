"""Cartan (connection-elimination) map on a Sol orbit: the full 24-component torsion map and the contortion quadratic form. Writes cartan_orbit.json."""
from common import *
import sympy as s
sys.path.insert(0,str(REPO/'script'))
from paper07_results import algebra

if __name__=='__main__':
    ck=Checks();data=algebra.contortion();T=data['torsion_map'];v=data['coordinates']
    M=s.hessian(data['scalar'],v)
    ck('full 24-component torsion map invertible',T.rows==24 and T.det()!=0)
    ck('full contortion quadratic form nondegenerate',M.rows==24 and M.det()!=0)
    ck('pointwise Cartan inverse exact',T*T.inv()==s.eye(24))
    ck('Cartan map has no spatial or mode-dependent coefficients',not T.free_symbols)
    save('cartan_orbit.json',dict(checks=ck.items,torsion_map=T,inverse=T.inv(),quadratic_hessian=M,
        orbit_argument='The invariant coframe and all tensor products descend to the fixed Sol quotient. The constant full Cartan map, and the linear Levi-Civita derivative source made from D1,D2,D3, act within the same infinite orbit. Thus auxiliary connection elimination preserves the orbit and adds no propagating connection field pair.',
        endpoint_scope='Retain the full NY current e^a wedge T_a and EH contortion current. On the full torsion-free connection solution the complete NY current vanishes identically; endpoint differences shift the pulled one-form by an exact variation and do not change its two-form. No new explicit full quadratic endpoint-current matrix or unreduced Dirac chain is computed in this package.',
        source='The public paper07_results/algebra.py implements the full contortion form; background Cartan and endpoint conventions are those of Appendices A-B.'))
    print(len(ck.items),'Cartan orbit checks',flush=True)
