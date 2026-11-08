"""Numerical helpers: parsing exact operator coefficients, evaluating them on grids, and normalized Hermite functions."""
from common import *
import sympy as s
import numpy as np
from scipy.special import roots_hermite
from scipy.integrate import solve_ivp
from scipy.fft import dct

def local_ops(raw):
    t,F=s.symbols('t F',positive=True);K,M=s.symbols('K M',real=True,nonzero=True)
    env={'t':t,'F':F,'K':K,'M':M}
    return [s.Matrix([[s.sympify(x,locals=env) for x in row] for row in A]) for A in raw],(t,F,K,M)

def evaluator(ops,variables):
    positions=[];expressions=[]
    for j,A in enumerate(ops):
        for a in range(A.rows):
            for b in range(A.cols):
                if A[a,b]!=0:positions.append((j,a,b));expressions.append(A[a,b])
    fn=s.lambdify(variables,expressions,'numpy',cse=True)
    def evaluate(t,K,M):
        vals=fn(t,np.exp(t*t)/np.sqrt(t),K,M)
        out=[np.zeros((len(K),A.rows,A.cols)) for A in ops]
        for (j,a,b),v in zip(positions,vals):out[j][:,a,b]=v
        return out
    return evaluate

def hermite_polynomials(nodes,N,order):
    size=N+order+1
    P=np.zeros((len(nodes),size));P[:,0]=np.pi**(-.25)
    if size>1:P[:,1]=np.sqrt(2)*nodes*P[:,0]
    for n in range(1,size-1):P[:,n+1]=np.sqrt(2/(n+1))*nodes*P[:,n]-np.sqrt(n/(n+1))*P[:,n-1]
    D=np.zeros((size,size))
    for n in range(size):
        if n:D[n-1,n]=np.sqrt(n/2)
        if n+1<size:D[n+1,n]=-np.sqrt((n+1)/2)
    derivs=[];B=np.eye(size)
    for j in range(order+1):
        derivs.append(P@B[:,:N]);B=D@B
    return derivs
