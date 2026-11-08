"""Wigner harmonics on the Taub three-sphere: generators and their commutation relations. Writes harmonics.json."""
from common import *
import sympy as s

def generators(j):
    j=s.Rational(j);d=int(2*j+1);labels=[-j+n for n in range(d)]
    plus=s.zeros(d)
    for n,m in enumerate(labels[:-1]):plus[n+1,n]=s.sqrt((j-m)*(j+m+1))
    minus=plus.T
    J=[(plus+minus)/2,(plus-minus)/(2*s.I),s.diag(*labels)]
    D=[2*s.I*a for a in J]
    S=s.Matrix(d,d,lambda a,b:(-1)**int(j-labels[a]) if b==d-1-a else 0)
    return D,S,labels

if __name__=='__main__':
    ck=Checks();results={}
    for j in [0,s.Rational(1,2),1,s.Rational(3,2)]:
        D,S,labels=generators(j);d=len(labels)
        for a in range(3):
            ck(f'j={j} derivative skew adjoint {a}',D[a].H==-D[a])
            ck(f'j={j} real field intertwiner {a}',s.conjugate(D[a])*S==S*D[a])
            for b in range(3):ck(f'j={j} structure {a}{b}',D[a]*D[b]-D[b]*D[a]==-2*sum((s.LeviCivita(a,b,c)*D[c] for c in range(3)),s.zeros(d)))
        ck(f'j={j} Casimir',sum((a*a for a in D),s.zeros(d))==-4*j*(j+1)*s.eye(d))
        ck(f'j={j} real structure square',S*S==(-1)**int(2*j)*s.eye(d))
        results[str(j)]=dict(generators=D,conjugation_intertwiner=S,labels=labels,scalar_harmonics=d*d,full_real_coefficient_dimension=d*d)
    save('harmonics.json',dict(checks=ck.items,results=results,conventions='sigma=-omega/2, [E_a,E_b]=-2 epsilon_abc E_c; E acts on a fixed left Wigner index through D_a=2i J_a on the right coefficient vector. Y_mn=sqrt(2j+1) D^j_mn, normalized Haar measure. Y_mn*=(-1)^(m-n)Y_-m,-n. A full fixed-j scalar coefficient set has (2j+1)^2 real dimensions; half-integer j is retained.',closure='Background frame coefficients depend only on time. Every derivative, constraint and gauge map is a finite tensor expression in these D_a and frame structure constants, so all right indices and tensor components at fixed j and left index close. The full left multiplicity is restored with the stated reality condition.'))
    print(len(ck.items),'harmonic checks',flush=True)
