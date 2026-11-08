"""Full Wigner block Cartan solution and complete endpoint representatives."""
from operators import *
from itertools import combinations,permutations
sys.path.insert(0,str(REPO/'script'))
from paper07_results import algebra
SIGN=algebra.SIGN

if __name__=='__main__':
    j=s.Rational(3,2);D,_,labels=generators(j);d=len(labels);ck=Checks()
    a,c,N=s.symbols('a c N',positive=True);Ha,Hc=s.symbols('H_a H_c',real=True);sc=[N,a,a,c];Hs=[0,Ha,Ha,Hc]
    C=[[[s.S(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for i,b,k in product(range(1,4),repeat=3):C[i][b][k]=-2*s.LeviCivita(i-1,b-1,k-1)*sc[i]/(sc[b]*sc[k])
    for i in range(1,4):C[i][0][i]=-Hs[i];C[i][i][0]=Hs[i]
    eta=s.diag(*SIGN);Gamma=connection(eta,C,lambda i,f:0)
    om={(a0,b0,i):SIGN[a0]*Gamma[a0,i,b0] for a0,b0,i in product(range(4),repeat=3)}
    nq=16*d;nv=32*d;zero=s.zeros(d,nv)
    inject=[]
    for r in range(16):
        M=zero.copy();M[:,r*d:(r+1)*d]=s.eye(d);inject.append(M)
    f={(a0,b0):zero.copy() for a0,b0 in product(range(4),repeat=2)}
    f[0,0]=inject[6]
    for i in range(3):f[i+1,0]=inject[7+i]
    for r,(i,k) in enumerate(PAIRS):f[i+1,k+1]=f[k+1,i+1]=inject[r]/2
    lor={(a0,b0):zero.copy() for a0,b0 in product(range(4),repeat=2)}
    for r,(i,k) in enumerate(combinations(range(4),2)):
        lor[i,k]=SIGN[i]*inject[10+r];lor[k,i]=-SIGN[k]*inject[10+r]
    f={key:val+lor[key] for key,val in f.items()}
    def deriv(i,M):
        if i:return D[i-1]*M/sc[i]
        return s.zeros(d,nq).row_join(M[:,:nq])/N
    rhs=[]
    for i in range(4):
        for b,k in combinations(range(4),2):
            rhs.append(deriv(b,f[i,k])-deriv(k,f[i,b])-sum((C[e][b][k]*f[i,e] for e in range(4)),zero)+sum((SIGN[i]*(om[i,e,b]*f[e,k]-om[i,e,k]*f[e,b]) for e in range(4)),zero))
    ref=algebra.contortion();Tm=s.Matrix(s.kronecker_product(ref['torsion_map'],s.eye(d)))
    W=red(-s.kronecker_product(ref['torsion_map'].inv(),s.eye(d))*s.Matrix.vstack(*rhs))
    ck('all 96 Cartan rows solved',red(Tm*W+s.Matrix.vstack(*rhs))==s.zeros(24*d,nv))
    ck('Cartan block rank 96',ref['torsion_map'].det()!=0)
    slots=[(i,k,b) for i,k in combinations(range(4),2) for b in range(4)]
    lorcols=list(range(10*d,16*d))+list(range(nq+10*d,nv))
    for r,(i,k,b) in enumerate(slots):
        expected=-SIGN[i]*(deriv(b,lor[i,k])+sum((SIGN[i]*om[i,e,b]*lor[e,k]-lor[i,e]*SIGN[e]*om[e,k,b] for e in range(4)),zero))
        ck('Lorentz covariance slot '+str(r),red(W[r*d:(r+1)*d,lorcols]-expected[:,lorcols])==s.zeros(d,12*d))
    M=ref['matrix'] if 'matrix' in ref else s.hessian(ref['scalar'],ref['coordinates'])
    ck('24-component auxiliary Hessian determinant',M.det()==2**32)
    # Spatially integrated endpoint matrices on q,w (not qdot).
    size=(16+24)*d;F={key:val[:,:nq].row_join(s.zeros(d,24*d)) for key,val in f.items()};z=s.zeros(d,size)
    ws={}
    for r,key in enumerate(slots):
        Z=z.copy();Z[:,nq+r*d:nq+(r+1)*d]=s.eye(d);ws[key]=Z
    def wc(i,k,b):return z if i==k else (1 if i<k else -1)*ws[min(i,k),max(i,k),b]
    ds=lambda i,Z:D[i-1]*Z/sc[i]
    TT={}
    for i in range(4):
        for b,k in combinations(range(1,4),2):TT[i,b,k]=ds(b,F[i,k])-ds(k,F[i,b])-sum((C[e][b][k]*F[i,e] for e in range(4)),z)+sum((SIGN[i]*(om[i,e,b]*F[e,k]-om[i,e,k]*F[e,b]) for e in range(4)),z)
    Bny=s.zeros(size);Beh=s.zeros(size)
    for i in range(4):
        for b in range(1,4):
            for k,l in combinations(range(1,4),2):Bny+=SIGN[i]*s.LeviCivita(b-1,k-1,l-1)*(F[i,b].H*TT[i,k,l])
    for i,k in product(range(4),repeat=2):
        for b,l,m in permutations(range(1,4)):
            if i==b:Bny+=2*s.LeviCivita(b-1,l-1,m-1)*(wc(i,k,l).H*F[k,m])
    for i,k,l,m in permutations(range(4)):
        eps=s.LeviCivita(i,k,l,m)*SIGN[l]*SIGN[m]/2
        for b,v,w in permutations(range(1,4)):
            if k==v:Beh+=eps*s.LeviCivita(b-1,v-1,w-1)*(F[i,b].H*wc(l,m,w))
            if i==b:Beh+=eps*s.LeviCivita(b-1,v-1,w-1)*(F[k,v].H*wc(l,m,w))
    Bny=red(a*a*c*(Bny+Bny.H)/2);Beh=red(a*a*c*(Beh+Beh.H)/2)
    ck('complete NY endpoint nonzero',Bny!=s.zeros(size));ck('contortion endpoint nonzero',Beh!=s.zeros(size))
    ck('both real endpoint quadratic forms Hermitian',Bny==Bny.H and Beh==Beh.H)
    pull=s.zeros(size,nv);pull[:nq,:nq]=s.eye(nq);pull[nq:,:]=W
    first_pull=red(pull.H*Bny*pull)
    # A quadratic endpoint pullback requires the mean second-order connection.
    # This auxiliary equation does not construct second-order gravitational data.
    wfirst={key:W[r*d:(r+1)*d,:] for r,key in enumerate(slots)}
    wf=lambda i,k,b:s.zeros(d,nv) if i==k else (1 if i<k else -1)*wfirst[min(i,k),max(i,k),b]
    mean=lambda U,V:(U.H*V+V.H*U)/2
    source2=[]
    for i in range(4):
        for b,k in combinations(range(4),2):
            source2.append(red(sum((SIGN[i]*(mean(wf(i,e,b),f[e,k])-mean(wf(i,e,k),f[e,b])) for e in range(4)),s.zeros(nv))))
    ti=ref['torsion_map'].inv()
    wsecond=[red(-sum((ti[r,v]*source2[v] for v in range(24)),s.zeros(nv))) for r in range(24)]
    for r in range(24):ck('mean second-order Cartan slot '+str(r),red(sum((ref['torsion_map'][r,v]*wsecond[v] for v in range(24)),s.zeros(nv))+source2[r])==s.zeros(nv))
    w2=dict(zip(slots,wsecond))
    get2=lambda i,k,b:(1 if i<k else -1)*w2[min(i,k),max(i,k),b]
    second_endpoint=red(a*a*c*sum((s.LeviCivita(i-1,b-1,k-1)*get2(i,k,b) for i,b,k in permutations(range(1,4))),s.zeros(nv)))
    ck('complete second-order NY pullback vanishes',red(first_pull+second_endpoint)==s.zeros(nv))
    # Sparse storage avoids thousands of literal zero entries.
    sparse=lambda mat:dict(shape=mat.shape,entries={str((i,k)):mat[i,k] for i in range(mat.rows) for k in range(mat.cols) if mat[i,k]!=0})
    save('connection_endpoints.json',dict(checks=ck.items,j=str(j),right_dimension=d,coframe_components=16,connection_components_per_left_copy=24*d,Cartan_solution=sparse(W),first_order_only_NY_pullback=sparse(first_pull),mean_second_order_NY_endpoint=sparse(second_endpoint),NY_endpoint=sparse(Bny),EH_contortion_endpoint=sparse(Beh),parameters='a=2 sqrt(1+tau^2), c=4 sqrt((1-tau^2)/(1+tau^2)), N=sqrt((1+tau^2)/(1-tau^2)); H_a=a_dot/(N a), H_c=c_dot/(N c). H_a,H_c are algebraically independent real expansion parameters in the identity checks; no positivity assumption is needed for them.',coframe_convention='q=(orthonormal metric six, fractional lapse, orthonormal shift three, lower Lorentz six), then its coordinate time derivatives. Endpoint matrices act on q and w. In NY, w=delta omega; in EH, w=independent contortion kappa. Hermitian forms are summed over all left indices with the Wigner real-field condition.',auxiliary_identity='L_EC2=L_ADM2+(N sqrt(g))/2 sum_m kappa_m^dag (M tensor I_d) kappa_m + dB_EH2/dtau. det(M)=2^32; complete constant NY adds theta dB_NY2/dtau. Endpoints shift the one-form by delta B and leave the pulled two-form unchanged.',scope='Explicit all-component Cartan and Lorentz checks, auxiliary action identity and endpoint representatives. Not a fresh unreduced 40-component Dirac chain, not a point canonical map, and not a separate bulk TT/EER re-contraction.'))
    print(len(ck.items),'connection and endpoint checks',flush=True)
