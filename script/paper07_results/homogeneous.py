"""Full symmetric homogeneous metrics on a fixed class-A coframe.

d sigma^i=(1/2)n_ij epsilon_jkl sigma^k wedge sigma^l.
c=V_ref/(2 kappa^2). Off-diagonal canonical momenta equal twice pi^ij.
F=c sqrt(det g) H is a positive rescaling of the scalar constraint.
"""
import sympy as s
PAIRS=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
def symmetric(v,dual=False):
    M=s.zeros(3)
    for (i,j),x in zip(PAIRS,v):M[i,j]=M[j,i]=x/(2 if dual and i!=j else 1)
    return M
def bracket(f,h,x,p):return s.expand(sum(s.diff(f,a)*s.diff(h,b)-s.diff(f,b)*s.diff(h,a) for a,b in zip(x,p)))
def system():
    x=s.symbols('g11 g22 g33 g12 g13 g23',real=True);p=s.symbols('p11 p22 p33 p12 p13 p23',real=True)
    g=symmetric(x);pi=symmetric(p,True);n=s.diag(*s.symbols('n1:4',real=True));c=s.Symbol('c',positive=True)
    P=s.expand(s.trace(n*g)**2/2-s.trace(n*g*n*g));T=s.expand(s.trace(pi*g*pi*g)-s.trace(pi*g)**2/2)
    A=[s.Matrix(3,3,lambda i,j:sum(n[i,z]*s.LeviCivita(z,k,j) for z in range(3))) for k in range(3)]
    J=[s.expand(s.trace(pi*(a.T*g+g*a))) for a in A]
    F=s.expand(T-c*c*P)
    return dict(x=x,p=p,g=g,pi=pi,n=n,c=c,P=P,T=T,F=F,J=J,A=A)
def representative(ns):
    z=system();g0=s.diag(1,2,3);n=s.diag(*ns);P=s.trace(n*g0)**2/2-s.trace(n*g0*n*g0)
    if P==0: mixed=s.diag(1,1,0)
    elif P<0: mixed=s.sqrt(-2*P/3)*s.eye(3)
    else:mixed=s.diag(s.sqrt(2*P),0,0)
    pi=mixed*g0.inv();subs=dict(zip(z['x'],[1,2,3,0,0,0]));subs.update(dict(zip(z['p'],[pi[0,0],pi[1,1],pi[2,2],0,0,0])))
    subs.update(dict(zip(z['n'].diagonal(),ns)));subs[z['c']]=1
    M=s.Matrix(z['J']).jacobian(z['x']+z['p']).subs(subs)
    C=s.Matrix([z['F']]+z['J']).jacobian(z['x']+z['p']).subs(subs)
    return dict(substitutions=subs,momentum_rank=M.rank(),constraint_rank=C.rank(),F=s.simplify(z['F'].subs(subs)),J=[s.simplify(j.subs(subs)) for j in z['J']],M=M,C=C)
def curvature(n,g):
    return s.factor((s.trace(n*g)**2/2-s.trace(n*g*n*g))/g.det())
def shape_residual(n):
    chi=(s.trace(n)**2/2-s.trace(n*n))/6
    return s.simplify(2*n*n-s.trace(n)*n+4*chi*s.eye(3))
