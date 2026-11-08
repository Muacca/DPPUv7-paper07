"""Linear reduction of one Taub Wigner block j: constraint ranks, reduced pairs, embedding and multipliers. Writes reduction_j<j>.json."""
from operators import *
from itertools import combinations
from exact_field import Field,ExpressionField

def read_operators(j):
    raw=json.loads((OUT/('operators_j'+str(j).replace('/','_')+'.json')).read_text(encoding='utf8'))['result']
    t=s.Symbol('tau',real=True);env={'tau':t}
    mat=lambda key:s.Matrix([[s.sympify(v,locals=env) for v in row] for row in raw[key]])
    return t,{key:mat(key) for key in ('H','C','J','g','pi','reality')}

def reduce_block(j,ck):
    t,data=read_operators(j);C,H,J=data['C'],data['H'],data['J'];d=int(2*j+1);nq=6*d
    transform=s.eye(2*nq);rowscale=s.eye(4*d)
    if j==s.Rational(3,2):
        # Constant canonical rescaling removes algebraic coefficients during elimination.
        weights=s.diag(1,s.sqrt(3),s.sqrt(3),1)
        Qscale=s.kronecker_product(s.diag(1,1,1,s.I,1,s.I),weights)
        transform=s.diag(Qscale,Qscale.inv().H);rowscale=s.kronecker_product(s.diag(1,-s.I,1,-s.I),weights.inv())
        C=red(rowscale*C*transform);H=red(transform.H*H*transform)
        ck('constant harmonic rescaling is canonical',transform.H*J*transform==J)
        ck('working matrices real rational over tau',not C.has(s.sqrt(3),s.I) and not H.has(s.sqrt(3),s.I))
    field=Field(None,t) if j==s.Rational(3,2) else ExpressionField();norm,mul=field.norm,field.mul
    at=s.Rational(1,3);C0=C.subs(t,at);rows=list(C0.T.rref()[1]);r=len(rows);Cr=C[rows,:];Cp=Cr[:,nq:]
    print('j',j,'constraint rank',r,'selecting gauge',flush=True)
    cp0=Cp.subs(t,at);piv=None
    if r==4*d:
        for comp in combinations(range(6),4):
            cols=[a*d+n for a in comp for n in range(d)]
            if cp0[:,cols].rank()==r:piv=cols;break
    if piv is None:piv=list(cp0.rref()[1])
    ck('momentum pivot rank equals complete constraint rank',len(piv)==r)
    free=[n for n in range(nq) if n not in piv];nf=len(free)
    pivot=Cp[:,piv];minor=s.factor(pivot.det(method='domain-ge'))
    print('j',j,'minor',minor,flush=True)
    pinv=field.inv(pivot)
    print('j',j,'pivot inverted',flush=True)
    if r==4*d:
        ck('all constraint rows represented',mul(pivot,pinv)==s.eye(r) and Cr==C)
    else:
        ck('all constraint rows represented',norm(C-mul(C[:,[nq+p for p in piv]],pinv,Cr))==s.zeros(*C.shape))
    # Constant configuration-coordinate gauge makes the reduced coordinates canonical.
    E=s.zeros(2*nq,2*nf)
    for a,n in enumerate(free):E[n,a]=1;E[nq+n,nf+a]=1
    slave=-mul(pinv,Cr,E)
    print('j',j,'embedding solved',flush=True)
    for a,n in enumerate(piv):E[nq+n,:]=slave[a,:]
    G=s.zeros(r,2*nq)
    for a,n in enumerate(piv):G[a,n]=1
    ck('constraints vanish on embedding',mul(C,E)==s.zeros(4*d,2*nf))
    ck('gauge vanishes on embedding',G*E==s.zeros(r,2*nf))
    Jr=s.zeros(2*nf);Jr[:nf,nf:]=s.eye(nf);Jr[nf:,:nf]=-s.eye(nf);Om=-J
    ck('reduced two form canonical',mul(E.H,Om,E)==-Jr)
    correction=mul(E.H,Om,s.diff(E,t));ck('time correction vanishes in this constant coordinate gauge',correction==s.zeros(2*nf))
    Hr=norm(mul(E.H,H,E)-correction)
    print('j',j,'reduced Hamiltonian obtained',flush=True)
    ck('reduced Hamiltonian Hermitian',norm(Hr-Hr.H)==s.zeros(2*nf))
    F=mul(G,J,Cr.H);multipliers=-mul(field.inv(F),G,J,H,E)
    print('j',j,'multipliers obtained',flush=True)
    ck('all original Hamilton equations reconstructed',norm(mul(J,H,E)+mul(J,Cr.H,multipliers)-s.diff(E,t)-mul(E,Jr,Hr))==s.zeros(2*nq,2*nf))
    # An independent null basis verifies constraint propagation before gauge fixing.
    nonp=[n for n in range(2*nq) if n not in [nq+p for p in piv]]
    null=s.zeros(2*nq,len(nonp))
    for a,n in enumerate(nonp):null[n,a]=1
    zsl=-mul(pinv,Cr,null)
    for a,n in enumerate(piv):null[nq+n,:]=zsl[a,:]
    prop=norm(s.diff(C,t)+mul(C,J,H))
    ck('time dependent full constraint propagation',mul(prop,null)==s.zeros(4*d,len(nonp)))
    # The minor has no zeros or poles on the selected real closed interval.
    num,den=s.fraction(minor)
    # Strip nonzero numerical algebraic constants; exact polynomial root counts.
    def realpoly(expr):
        pol=s.Poly(expr,t);return s.Poly(pol.as_expr()/pol.LC(),t)
    zero_count=realpoly(num).count_roots(s.Rational(1,4),s.Rational(1,2))
    pole_count=realpoly(den).count_roots(s.Rational(1,4),s.Rational(1,2))
    ck('gauge minor nonzero throughout interval',zero_count==0 and pole_count==0)
    ck('same rank at both interval endpoints',C.subs(t,s.Rational(1,4)).rank()==r and C.subs(t,s.Rational(1,2)).rank()==r)
    select=s.zeros(2*nf,2*nq)
    for a,n in enumerate(free):select[a,n]=1;select[nf+a,nq+n]=1
    Eoriginal=transform*E;select_original=select*transform.inv();Goriginal=G*transform.inv()
    row_sub=rowscale.extract(rows,rows)
    Cr_original=row_sub.inv()*Cr*transform.inv();mult_original=row_sub.H*multipliers
    return dict(checks=ck.items,j=str(j),right_dimension=d,left_multiplicity=d,full_real_metric_dimension=6*d*d,linear_constraint_rank_per_left_copy=r,linear_reduced_pairs_per_left_copy=nf,linear_reduced_pairs_full_real_j=d*nf,constraint_rows=rows,gauge_coordinate_indices=piv,free_coordinate_indices=free,minor=minor,minor_basis='working constant canonical rescaling; original normalized Wigner basis restored in embedding and multiplier matrices',minor_roots_in_interval=int(zero_count),minor_poles_in_interval=int(pole_count),embedding=Eoriginal,physical_selector=select_original,gauge=Goriginal,independent_constraints=Cr_original,reduced_J=Jr,reduced_H=Hr,time_dependent_correction=correction,multipliers_for_independent_rows=mult_original,constant_canonical_basis_transform=transform,scope='Linearized constraint and gauge quotient at the chosen symmetric background. For j=0 and j=1, missing KID rows mean these dimensions are not dimensions of a smooth nonlinear physical phase space; quadratic KID charges must also be imposed.')

if __name__=='__main__':
    j=s.Rational(sys.argv[1]);ck=Checks();r=reduce_block(j,ck)
    save('reduction_j'+str(j).replace('/','_')+'.json',r);print('reduced',j,len(ck.items),flush=True)
