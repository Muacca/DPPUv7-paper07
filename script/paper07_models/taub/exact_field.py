"""Exact matrix arithmetic helpers (expression and rational fields), used by reduce_block.py."""
import sympy as s
from sympy.polys.matrices import DomainMatrix

class ExpressionField:
    def norm(self,M):return M.applyfunc(lambda x:s.cancel(s.expand(x)))
    def mul(self,*matrices):
        v=matrices[0]
        for M in matrices[1:]:v=v*M
        return self.norm(v)
    def inv(self,M):return self.norm(M.inv(method='DM'))

class Field:
    def __init__(self,j,t):
        if j is None:
            self.domain=s.QQ.frac_field(t)
            return
        ext=[s.I]+([s.sqrt(2)] if j==1 else [s.sqrt(3)] if j==s.Rational(3,2) else [])
        self.domain=s.QQ.algebraic_field(*ext).frac_field(t)
    def dm(self,M):
        return DomainMatrix([[self.domain.from_sympy(x) for x in row] for row in M.tolist()],M.shape,self.domain).to_sparse()
    def norm(self,M):return self.dm(M).to_Matrix()
    def mul(self,*matrices):
        v=self.dm(matrices[0])
        for M in matrices[1:]:v=v*self.dm(M)
        return v.to_Matrix()
    def inv(self,M):return self.dm(M).inv().to_Matrix()
