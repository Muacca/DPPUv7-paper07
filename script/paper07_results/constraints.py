"""Dirac closure for constant-coefficient quadratic finite-dimensional actions.

This helper deliberately rejects nonlinear constraints. Their singular strata
must be studied in the original constraint, not by linear coefficient ranks.
"""
from itertools import combinations
import sympy as s
from .homogeneous import bracket
def quadratic_dirac(L,q,v):
    q=list(q);v=list(v);p=list(s.symbols('p0:'+str(len(q))));z=q+p
    try:degree=s.Poly(L,*(q+v)).total_degree()
    except s.PolynomialError as exc:raise ValueError('action must be polynomial in coordinates and velocities') from exc
    if degree>2:raise ValueError('action must be quadratic in coordinates and velocities')
    W=s.hessian(L,v)
    if any(x.free_symbols for x in W):raise ValueError('velocity Hessian must be constant')
    rank=W.rank();idx=next((a for a in combinations(range(len(q)),rank) if W.extract(a,a).det()!=0),())
    other=[i for i in range(len(q)) if i not in idx]
    zero={v[i]:0 for i in other}
    sol=dict(zero)
    if idx:
        b=s.Matrix([p[i]-s.diff(L,v[i]).subs({w:0 for w in v}) for i in idx])
        sol.update(dict(zip([v[i] for i in idx],W.extract(idx,idx).inv()*b)))
    phi=[s.expand(p[j]-s.diff(L,v[j]).subs(sol,simultaneous=True)) for j in other]
    H=s.expand((sum(p[i]*v[i] for i in range(len(q)))-L).subs(sol,simultaneous=True))
    constraints=[];generations=[]
    def row(f):
        if any(s.diff(f,a,b)!=0 for a in z for b in z):raise ValueError('nonlinear constraint: '+str(f))
        return [s.diff(f,a) for a in z]+[f.subs({a:0 for a in z})]
    def insert(f):
        f=s.expand(f)
        if f==0:return False
        if not (f.free_symbols & set(z)):raise ValueError('inconsistent constant constraint')
        rows=[row(g) for g in constraints]
        if s.Matrix(rows+[row(f)]).rank()> (s.Matrix(rows).rank() if rows else 0):constraints.append(f);return True
        return False
    for f in phi:insert(f)
    generations.append(len(constraints))
    for _ in range(12):
        # Total Hamiltonian: only ORIGINAL primary multipliers are arbitrary.
        # Giving secondary constraints independent multipliers would miss chains.
        consistency=s.Matrix([[bracket(f,g,q,p) for g in phi] for f in constraints])
        additions=[]
        for null in consistency.T.nullspace():
            additions.append(s.expand(sum(null[i]*bracket(f,H,q,p) for i,f in enumerate(constraints))))
        num=sum(insert(f) for f in additions)
        if num==0:break
        generations.append(num)
    else:raise RuntimeError('closure did not terminate')
    C=s.Matrix([[bracket(f,g,q,p) for g in constraints] for f in constraints]);r=C.rank();nf=len(constraints)-r
    return dict(velocity_rank=rank,primary_count=len(phi),generations=generations,constraints=constraints,
                poisson=C,second_class=r,first_class=nf,physical_pairs=s.Rational(2*len(q)-2*nf-r,2),H=H,q=q,p=p)
def scalar_routes(k):
    q,N,a,b,dq=s.symbols('q N a b dq',real=True)
    # c=1, k=1 is unit round; k=0 is the flat common-scale family.
    L=s.expand(6*N*q**3*b**2-12*q*q*dq*b+6*N*q**3*(2*k*a/q-a*a))
    eliminated=s.factor(L.subs({a:k/q,b:dq/(N*q)}))
    eps=s.Symbol('eps');u,n,A,B,du=s.symbols('u n A B du',real=True)
    repl={q:1+eps*u,N:1+eps*n,a:k+eps*A,b:eps*B,dq:eps*du}
    expanded=s.expand(L.subs(repl,simultaneous=True));L1=expanded.coeff(eps,1);L2=expanded.coeff(eps,2)
    reduced2=s.series(eliminated.subs(repl,simultaneous=True),eps,0,3).removeO().expand().coeff(eps,2)
    full=quadratic_dirac(L2,[u,n,A,B],[du,*s.symbols('dn dA dB')])
    reduced=quadratic_dirac(reduced2,[u,n],[du,s.Symbol('dn')])
    return dict(exact=L,eliminated=eliminated,linear=L1,quadratic=L2,reduced_quadratic=reduced2,full=full,reduced=reduced)
def legacy_trace():
    """The v1 l=0 trace block at N=q=1; a historical algebraic control."""
    u,n,a,b,db,du=s.symbols('s n a b db ds',real=True)
    L=-8*u*db+12*b*b-12*a*a-4*a*u-12*n*a-4*n*u
    ibp=s.expand(L+8*u*db+8*du*b)
    reduced=s.factor(ibp.subs({a:-u/6-n/2,b:-du/3}))
    full=quadratic_dirac(L,[u,n,a,b],[du,*s.symbols('dn da'),db])
    metric=quadratic_dirac(reduced,[u,n],[du,s.Symbol('dn')])
    # Reproduce the v1 erroneous rule phi -> {phi,Hc}=0, omitting multipliers.
    p=s.symbols('P0:4');q=[u,n,a,b];phi=[p[0],p[1],p[2],p[3]+8*u];H=-L.subs(db,0)
    candidates=phi+[bracket(f,H,q,p) for f in phi];rows=[];kept=[]
    for f in candidates:
        row=[s.diff(f,x) for x in q+list(p)]
        if s.Matrix(rows+[row]).rank()>(s.Matrix(rows).rank() if rows else 0):rows.append(row);kept.append(f)
    C=s.Matrix([[bracket(f,g,q,p) for g in kept] for f in kept]);r=C.rank();nf=len(kept)-r
    wrong=s.Rational(8-2*nf-r,2)
    return dict(first_order=L,integrated=ibp,reduced=reduced,full=full,metric=metric,naive_pairs=wrong)
