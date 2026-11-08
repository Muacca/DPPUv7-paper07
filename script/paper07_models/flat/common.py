"""Output conventions for the mathematical examples in this directory."""
from pathlib import Path
import json,os,sys
WORK=Path(__file__).resolve().parent
REPO=WORK.parents[2]
if 'PAPER07_OUTPUT' not in os.environ:
    raise RuntimeError('Use python script/perturbations.py recompute --model ... --output NEW_DIRECTORY')
OUT=Path(os.environ['PAPER07_OUTPUT']).resolve()
def write(p,t):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(t.rstrip()+'\n')
def serial(x):
    import sympy as s
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    if isinstance(x,s.MatrixBase):return serial(x.tolist())
    if isinstance(x,s.Basic):return str(x)
    return x
def save(name,obj):write(OUT/name,json.dumps(serial(obj),ensure_ascii=False,indent=2))
class Checks:
    def __init__(self):self.items=[]
    def __call__(self,name,ok):
        self.items.append(dict(name=name,passed=bool(ok)))
        if not ok:raise AssertionError(name)
    def zero(self,name,expr):
        import sympy as s
        v=s.simplify(s.expand(expr))
        if v!=0:
            # Factor radicands so positive-time perfect squares are recognized.
            # No forced power rules or numerical tolerance are used.
            v=v.replace(lambda z:z.is_Pow and z.exp.is_Rational and not z.exp.is_Integer,
                        lambda z:s.Pow(s.factor(z.base),z.exp))
            v=s.simplify(v)
        self(name,v==0)
