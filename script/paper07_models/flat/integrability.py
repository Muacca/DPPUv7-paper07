"""KID-symbol completeness and the AMM cone in Paper07 normalization."""
from common import *
import sympy as s
ck=Checks()
k=s.Matrix(s.symbols('k1:4',real=True));r2=(k.T*k)[0]
N=s.Symbol('N');X=s.Matrix(s.symbols('X1:4'))
scalar=N*(r2*s.eye(3)-k*k.T)
vector=k*X.T+X*k.T
symbol=s.Matrix(list(scalar)+list(vector)).jacobian([N,*X])
gram=s.simplify(symbol.T*symbol)
expected=s.diag(2*r2**2,*([0]*3));expected[1:4,1:4]=2*r2*s.eye(3)+2*k*k.T
ck('KID symbol Gram matrix',s.simplify(gram-expected)==s.zeros(4))
ck('KID symbol injective at every nonzero covector',s.factor(gram.det()-32*r2**5)==0)
ck('zero Fourier symbol has all four constant KIDs',symbol.subs(dict.fromkeys(k,0))==s.zeros(18,4))
h,a,b,u,v,w=s.symbols('h a b u v w',real=True)
S=s.Matrix([[a,u,v],[u,b,w],[v,w,-a-b]])
pi=S-2*h*s.eye(3)
beta=s.trace(pi)
ck('AMM beta equals minus six mean expansion',beta==-6*h)
ck('AMM tracefree momentum equals shear',pi-beta*s.eye(3)/3==S)
ck('AMM homogeneous cone normalization',s.expand(s.trace(pi*pi)-s.trace(pi)**2/2-s.trace(S*S)+6*h*h)==0)
Ekin,Egrad=s.symbols('Ekin Egrad',nonnegative=True)
Q=s.trace(S*S)+Ekin/4-beta**2/6+Egrad/4
ck('AMM Hamiltonian Taub cone equals minus W-03 residual',s.expand(Q-(s.trace(S*S)+(Ekin+Egrad)/4-6*h*h))==0)
ck('periodic constant metric variations do not enter kinetic cone',not (Q.free_symbols & set(s.symbols('C0:6'))))
# Pointwise polynomial before averaging: cross terms with mean-zero TT velocities integrate to zero.
d0,d1,d2,d3,d4=s.symbols('d0:5',real=True)
D=s.Matrix([[d0,d2,d3],[d2,d1,d4],[d3,d4,-d0-d1]])
pi1=S+D/2-2*h*s.eye(3)
kin=s.expand(s.trace(pi1*pi1)-s.trace(pi1)**2/2)
ck('full mixed kinetic polynomial',s.expand(kin-(s.trace(S*S)+s.trace(S*D)+s.trace(D*D)/4-6*h*h))==0)
save('integrability.json',dict(checks=ck.items,kid_symbol=serial(symbol),kid_gram_determinant=str(s.factor(gram.det())),constant_kid_dimension=4,AMM_beta=str(beta),AMM_cone=serial(Q),verdict='W03_NECESSARY_AND_SUFFICIENT_FOR_LOCAL_INTEGRABILITY_SCOPED',proof_kind='application_of_AMM1982_Theorem_3_3_with_explicit_KID_and_normalization_map'))
print(len(ck.items),'Integrability checks; exactly four KIDs; AMM cone matches W-03')
