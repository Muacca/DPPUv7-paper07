"""Maps from the constrained initial ADM data to the Lorenz system, and the projection that removes gauge directions before comparison with the reduced variables. Writes initial_maps.json."""
from common import *
from numerics import local_ops
from fastops import *
def mat(raw):return local_ops([raw])[0][0]
data=json.loads((OUT/'harmonic_operators.json').read_text())
base=json.loads((OUT/'inputs/reduction_core.json').read_text())
E,args=local_ops(base['embedding']);H,_=local_ops(base['H_work']);C,_=local_ops(base['C_work'])
J=mat(base['J']);Jr=mat(base['reduced_J']);ck=Checks()
B0,Bt,Bx=[mat(data[n]) for n in ['B0','Bt','Bx']]
Lor=mat(data['lorenz']);L0,Lt,Lx=Lor[:,:10],Lor[:,10:20],Lor[:,20:30]
Z=mat(data['adm']);Zh=[Z[:,:10],Z[:,20:30]];Zv=[Z[:,10:20]]
I=s.zeros(10,12);I[4:,:6]=s.eye(6)
A=compose([J],H);ck('initial ADM metric velocity algebraic',all(a[:6,:]==s.zeros(6,12) for a in A[1:]))
V0=s.zeros(10,12);V0[4:,:]=A[0][:6,:]
V1=s.zeros(10,12);inv=Lt[:,:4].inv()
V0[:4,:]=red(-inv*(L0*I+Lt*V0));V1[:4,:]=red(-inv*Lx*I)
V=[V0,V1]
ck('initial Lorenz constraints exactly zero',iszero(add(compose([L0,Lx],[I]),compose([Lt],V))))
ck('same initial ADM metric and momentum',iszero(add(compose(Zh,[I]),compose(Zv,V),[-s.eye(12)])))
# Differentiate Lorenz F along the independent wave equation.
Fd_h=[red(dt(L0)+Lt*B0),red(dt(Lx)+Lt*Bx),Lt]
Fd_v=[red(L0+dt(Lt)+Lt*Bt),Lx]
Fdt_init=add(compose(Fd_h,[I]),compose(Fd_v,V))
ck('Lorenz time derivative zero on constrained initial data',iszero(compose(Fdt_init,E)))
projection=compose([-Jr],compose(adjoint(E),[J]))
ck('reduced projection left-inverts embedding',iszero(add(compose(projection,E),[-s.eye(4)])))
ck('reduced projection kills all gauge generators',iszero(compose(projection,compose([J],adjoint(C)))))
Yh=compose(projection,Zh);Yv=compose(projection,Zv)
initial_h=compose([I],E);initial_v=compose(V,E)
Ydh=add(timed(Yh),compose(Yv,[B0,Bx,s.eye(10)]))
Ydv=add(Yh,timed(Yv),compose(Yv,[Bt]))
ch=compose(C,Zh);cv=compose(C,Zv)
save('initial_maps.json',dict(checks=ck.items,initial_h=initial_h,initial_v=initial_v,initial_h_from_ADM=[I],initial_v_from_ADM=V,
       adm_h=Zh,adm_v=Zv,reduced_projection=projection,reduced_h=Yh,reduced_v=Yv,reduced_dt_h=Ydh,reduced_dt_v=Ydv,
       constraint_h=ch,constraint_v=cv,lorenz_dt_h=Fd_h,lorenz_dt_v=Fd_v,
       initial_gauge='H00=H0i=0 at t=1. Spatial H and Hdot reproduce the unchanged constrained ADM initial data. Hdot00,Hdot0i solve Lorenz F=0. Independent hyperbolic equations imply Fdot=0 on that initial data, checked exactly.',
       reduced_map='y=-Jr E^dagger J z, valid on the constraint surface and annihilating the full linearized gauge generators. Changes of gauge are not compared by equating raw four-dimensional metric components.'))
print(len(ck.items),'exact mapping checks',flush=True)
