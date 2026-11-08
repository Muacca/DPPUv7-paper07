"""Highest-derivative coefficient of the reduced fourth-order operator: rank one and nilpotent after multiplication by J_r. Writes principal_check.json."""
from common import *
from fastops import *
from numerics import local_ops
raw=json.loads((OUT/'expanded_Hr.json').read_text())
H,args=local_ops(raw['coefficients']);A=H[-1]
ck=Checks()
ck('fourth derivative coefficient nonzero',A[3,3]!=0)
ck('fourth derivative coefficient rank one by exact outer product',red(A*A[3,3]-A[:,3]*A[3,:])==s.zeros(4))
J=s.zeros(4);J[:2,2:]=s.eye(2);J[2:,:2]=-s.eye(2)
ck('leading generator coefficient square zero',red(J*A*J*A)==s.zeros(4))
save('principal_check.json',dict(checks=ck.items,order=4,leading_H_rank=1,leading_JH_nilpotence_index=2,
       scope='Exact local highest-derivative coefficients. This alone does not determine full hyperbolicity, dispersion, wellposedness, or instability.'))
print(ck.items,flush=True)
