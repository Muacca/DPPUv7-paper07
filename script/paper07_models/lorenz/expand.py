"""Exact expansion of the reduced fourth-order operator H_r = E^dagger H E from the canonical reduction input. Writes expanded_Hr.json."""
from common import *
from fastops import *
from numerics import local_ops
raw=json.loads((OUT/'inputs/reduction_core.json').read_text())
E,args=local_ops(raw['embedding']);H,_=local_ops(raw['H_work'])
print('Compose H E',flush=True)
HE=compose(H,E)
print('Compose E adjoint H E',flush=True)
Hr=compose(adjoint(E),HE)
save('expanded_Hr.json',dict(coefficients=Hr,order=len(Hr)-1,method='Exact rational composition over QQ(t,F,K,M); same E and H as the canonical reduction input.'))
print('Saved order',len(Hr)-1,flush=True)
