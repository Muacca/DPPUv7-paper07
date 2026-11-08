"""Linear gauge reduction of the central-nonzero Nil blocks, real pair counts and the central KID condition. Writes reduction.json."""
from operators import *
from itertools import combinations
from rational_field import RationalField

def parse_matrix(x,env):return s.Matrix([[s.sympify(v,locals=env) for v in row] for row in x]) if x else s.zeros(0,0)
def load_block(L,size=7):
    u=s.Symbol('u',positive=True);b=s.Symbol('b',positive=True);env={'u':u,'b':b}
    raw=json.loads((OUT/f'operators_g{size}.json').read_text(encoding='utf8'));block=raw['blocks'][str(L)]
    nq=block['metric_dimension_complex'];nc=block['constraint_rows']
    mats={key:parse_matrix(block[key],env) for key in ('H','C','J','original_embedding','multiplier_embedding')}
    if nc==0:mats['C']=s.zeros(0,2*nq);mats['multiplier_embedding']=s.zeros(4*size,0)
    return u,b,block,mats

def positive_polynomial(expr,u,b):
    v=s.Symbol('v',nonnegative=True)
    poly=s.Poly(s.expand(expr.subs(u,v+1)),v,b);normalized=s.Poly(poly.as_expr()/poly.LC(),v,b)
    coeff=normalized.coeffs();at_left=s.Poly(normalized.as_expr().subs(v,0),b)
    return all(c.is_nonnegative for c in coeff) and at_left.as_expr()!=0,dict(shifted_normalized=normalized.as_expr(),coefficients=coeff,left_boundary=at_left.as_expr())

def reduce(L,ck):
    u,b,raw,a=load_block(L);C,H,J=a['C'],a['H'],a['J'];nq=H.rows//2;r=C.rows
    weight=lambda n:s.sqrt(2**n*s.factorial(n))
    qw=[weight(n)*(s.sqrt(2) if q==2 else 1)*(s.I if SPINS[q]%2 else 1) for q,n in raw['metric_slots']]
    mw=[weight(n)*(s.sqrt(2) if q in [1,3] else 1)*(s.I if q in [1,3] else 1)*(s.I if q else 1) for q,n in raw['multiplier_slots']]
    Q=s.diag(*qw);T=s.diag(Q,Q.inv().H);M=s.diag(*mw)
    H=red(T.H*H*T);C=red(M.H*C*T)
    ck(f'L={L} constant basis change canonical',red(T.H*J*T)==J)
    ck(f'L={L} working matrices real rational',all(x.is_rational_function(u,b) and not x.has(s.I) for x in list(H)+list(C)))
    field=RationalField(u,b);norm,mul=field.norm,field.mul
    ck(f'L={L} full constraint row rank at control point',C.subs({u:s.Rational(3,2),b:1}).rank()==r)
    free=list(range(nq));piv=[];minor=s.S(1);cert={};pinv=s.zeros(0)
    if r:
        for cols in combinations(range(nq),r):
            pivot=C[:,[nq+i for i in cols]];det=s.factor(pivot.det())
            if det==0:continue
            num,den=s.fraction(det);ok1,c1=positive_polynomial(num,u,b);ok2,c2=positive_polynomial(den,u,b)
            if ok1 and ok2:piv=list(cols);minor=det;cert=dict(numerator=c1,denominator=c2);pinv=field.inv(pivot);break
        ck(f'L={L} a coordinate gauge valid for all u>=1 and b>0',len(piv)==r)
    free=[i for i in range(nq) if i not in piv];nf=len(free);E=s.zeros(2*nq,2*nf);G=s.zeros(r,2*nq)
    for i,k in enumerate(free):E[k,i]=1;E[nq+k,nf+i]=1
    if r:
        slave=-mul(pinv,C,E)
        for i,k in enumerate(piv):E[nq+k,:]=slave[i,:];G[i,k]=1
    Jr=s.zeros(2*nf);Jr[:nf,nf:]=s.eye(nf);Jr[nf:,:nf]=-s.eye(nf)
    ck(f'L={L} constraints eliminated',mul(C,E)==s.zeros(r,2*nf))
    ck(f'L={L} canonical pullback',mul(E.H,-J,E)==-Jr)
    correction=mul(E.H,-J,s.diff(E,u));ck(f'L={L} time correction evaluated',correction==s.zeros(2*nf))
    Hr=norm(mul(E.H,H,E)-correction)
    F=G*J*C.H
    mult=-mul(field.inv(F),G,J,H,E) if r else s.zeros(0,2*nf)
    ck(f'L={L} all full equations reconstructed',norm(mul(J,H,E)+mul(J,C.H,mult)-s.diff(E,u)-mul(E,Jr,Hr))==s.zeros(2*nq,2*nf))
    # A basis for ker(C) before fixing a gauge.
    remaining=[i for i in range(2*nq) if i not in [nq+k for k in piv]];null=s.zeros(2*nq,len(remaining))
    for i,k in enumerate(remaining):null[k,i]=1
    if r:
        slave=-mul(pinv,C,null)
        for i,k in enumerate(piv):null[nq+k,:]=slave[i,:]
    ck(f'L={L} time-dependent constraint propagation',mul(norm(s.diff(C,u)+mul(C,J,H)),null)==s.zeros(r,len(remaining)))
    ck(f'L={L} reduced Hamiltonian Hermitian',Hr==Hr.H)
    # Central KID moment map per positive-frequency representative and its conjugate.
    qform=s.zeros(2*nq);qform[:nq,nq:]=-s.I*b*b*s.eye(nq);qform[nq:,:nq]=s.I*b*b*s.eye(nq)
    charge=s.I*mul(E.H,qform/s.I,E)
    ck(f'L={L} central charge canonical',charge==s.I*b*b*(-Jr))
    ck(f'L={L} central charge conserved',norm(s.diff(charge/s.I,u)+mul((Jr*Hr).H,charge/s.I)+mul(charge/s.I,Jr,Hr))==s.zeros(2*nf))
    return dict(L=L,metric_slots=raw['metric_slots'],multiplier_slots=raw['multiplier_slots'],metric_dimension_complex=nq,constraint_rank_complex=r,reduced_pairs_complex=nf,reduced_pairs_real_per_fixed_positive_m_and_r=2*nf,momentum_minor=minor,minor_basis='constant canonical working basis',nonzero_certificate=cert,gauge_indices=piv,free_indices=free,embedding=T*E,gauge=G*T.inv(),reduced_J=Jr,reduced_H=Hr,time_dependent_correction=correction,multiplier_reconstruction=M*mult,central_charge_matrix=charge,constant_canonical_transform=T,multiplier_basis_transform=M)

if __name__=='__main__':
    ck=Checks();a=json.loads((OUT/'operators_g7.json').read_text(encoding='utf8'));bb=json.loads((OUT/'operators_g8.json').read_text(encoding='utf8'))
    env={'u':s.Symbol('u',positive=True),'b':s.Symbol('b',positive=True)}
    for L in range(-2,3):
        for key in ['H','C','J']:
            x=parse_matrix(a['blocks'][str(L)][key],env);y=parse_matrix(bb['blocks'][str(L)][key],env)
            ck(f'L={L} guard-size independent {key}',red(x-y)==s.zeros(*x.shape))
    results={}
    for L in range(-2,3):
        print('Reducing L',L,flush=True);results[str(L)]=reduce(L,ck)
    save('reduction.json',dict(checks=ck.items,blocks=results,frequency='b=sqrt(mu)>0; physical mu=m/(2pi), m positive integer. All identities and gauge certificates hold for every b>0, hence for every permitted m.',rank_argument='For the LRS compact Nil background the only global KID is central Z, which has central index zero. Closure makes a null vector of an individual block constraint adjoint a genuine smooth KID. There is none in m!=0, so each complete block has full row rank. This extends the rank count to all L>=2; symbolic matrices below explicitly check L=-2,-1,0,1,2.',scope='Linear gauge reduction and necessary central KID condition. Real pairs include the conjugate negative central sector, and must be multiplied by abs(m) for the lattice multiplicity.'))
    print(len(ck.items),'reduction and guard checks',flush=True)
