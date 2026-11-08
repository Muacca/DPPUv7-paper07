"""Canonical reduction core of the Sol orbit (constraints eliminated, canonical pullback), the input of the Lorenz comparison. Writes reduction_core.json."""
from reduce_orbit import *

if __name__=='__main__':
    ck=Checks();H,C,J=load_data();E,G,parts=lifts(ck,C)
    Jr=s.zeros(4);Jr[:2,2:]=s.eye(2);Jr[2:,:2]=-s.eye(2)
    ck('constraints eliminated',iszero(compose(C,E)))
    ck('canonical pullback',iszero(add(compose(compose(adjoint(E),[-J]),E),[Jr])))
    save('reduction_core.json',dict(checks=ck.items,embedding=E,gauge=G,parts=parts,reduced_J=Jr,H_work=H,C_work=C,J=J,
        coefficient_variables='t,F,K,M with F=exp(t^2)/sqrt(t), K=k0 exp(-xi), M=l0 exp(xi)',
        phase='original h,p=diag(1,1,1,i,i,1) times working h,p; original constraints=diag(1,1,i,i) times working constraints',
        reduced_H_factorization='E^dagger H_work E, since the verified time correction vanishes; all adjoints are differential adjoints on the specified smooth orbit domain.'))
    print(len(ck.items),'core reduction checks',flush=True)
