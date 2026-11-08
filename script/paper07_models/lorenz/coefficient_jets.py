"""Spatial and time derivatives (jets) of the maps to the reduced variables and of the Lorenz constraint operator, used by the residual audit. Writes coefficient_jets.json."""
from common import *
from fastops import *
from numerics import local_ops
raw=json.loads((OUT/'initial_maps.json').read_text())
result={}
for key in ['reduced_h','reduced_v']:
    aa,args=local_ops(raw[key]);der=[aa]
    for i in range(4):der.append([dx(a) for a in der[-1]])
    result[key]=der;result[key+'_t']=[dt(a) for a in aa]
base=json.loads((OUT/'harmonic_operators.json').read_text())
L=local_ops([base['lorenz']])[0][0]
result['lorenz_t']=dt(L)
save('coefficient_jets.json',result)
print('Saved coefficient jets for independent reduced residual',flush=True)
