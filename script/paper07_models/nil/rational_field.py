"""Exact matrix arithmetic over a rational function field (SymPy DomainMatrix), used by reduce_blocks.py."""
import sympy as s
from sympy.polys.matrices import DomainMatrix
class RationalField:
    def __init__(self,*variables):self.domain=s.QQ.frac_field(*variables)
    def dm(self,M):return DomainMatrix([[self.domain.from_sympy(x) for x in row] for row in M.tolist()],M.shape,self.domain).to_sparse()
    def norm(self,M):return self.dm(M).to_Matrix()
    def mul(self,*matrices):
        d=self.dm(matrices[0])
        for M in matrices[1:]:d=d*self.dm(M)
        return d.to_Matrix()
    def inv(self,M):return self.dm(M).inv().to_Matrix()
