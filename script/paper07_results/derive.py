"""Exact producers and meaningful algebraic checks for the paper."""
import json,hashlib
from pathlib import Path
import sympy as s
from . import algebra,homogeneous,geometric,torus,constraints
def serial(x):
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    if isinstance(x,s.MatrixBase):return [[serial(v) for v in row] for row in x.tolist()]
    if isinstance(x,s.Basic):return str(x)
    return x
def closure_summary(r):return {k:r[k] for k in ['velocity_rank','primary_count','generations','first_class','second_class','physical_pairs','constraints','poisson','H']}
def derive():
    checks=[]
    def check(name,condition):
        ok=bool(condition);checks.append(dict(name=name,passed=ok))
        if not ok:raise AssertionError(name)
    c=algebra.contortion();check('full contortion Hessian invertible',c['hessian'].rank()==24);check('torsion coordinates invertible',c['torsion_map'].rank()==24)
    check('scalar is exact quadratic',s.expand(c['scalar']-(s.Matrix(c['coordinates']).T*c['hessian']*s.Matrix(c['coordinates']))[0]/2)==0)
    f=algebra.class_a_forms();volume=(0,1,2,3)
    check('complete NY exterior identity',algebra.add(f['TT'],algebra.scale(-1,f['EER']),algebra.scale(-1,f['NY']))=={})
    check('epsilon action equals twice scalar',s.expand(f['epsilon_EER'].get(volume,0)-2*f['R'])==0)
    check('boundary Q geometric sign',s.simplify(f['Q'][(1,2,3)]-6*f['eta']/f['q'])==0)
    scalar_lc=f['R'].subs({f['eta']:0,f['V']:0,f['VD']:0})
    delta=s.factor(f['R']-scalar_lc)
    check('vector divergence and torsion square',s.expand(delta+6*f['eta']**2/f['q']**2-s.Rational(2,3)*f['V']**2+2*f['VD']+6*f['H']*f['V'])==0)
    action=dict(scalar=f['R'],scalar_difference=delta,Q=f['Q'],NY=f['NY'],TT=f['TT'],EER=f['EER'],epsilon_to_scalar_ratio=2,contortion_scalar=c['scalar'],contortion_coordinates=c['coordinates'],contortion_hessian=c['hessian'],contortion_determinant=c['hessian'].det(),torsion_map_determinant=c['torsion_map'].det(),divergence_current=c['current'])
    z=homogeneous.system();x=z['x'];p=z['p']
    for i in range(3):check(f'momentum {i} preserves scalar',homogeneous.bracket(z['J'][i],z['F'],x,p)==0)
    for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:check(f'momentum algebra {i}{j}',s.expand(homogeneous.bracket(z['J'][i],z['J'][j],x,p)+z['n'][k,k]*z['J'][k])==0)
    regular=[]
    for name,ns,rank in [('T3',(0,0,0),0),('Nil',(0,0,-1),2),('Sol',(0,1,-1),3),('S3',(2,2,2),3)]:
        r=homogeneous.representative(ns);check(name+' constraint shell',r['F']==0 and all(j==0 for j in r['J']));check(name+' regular Jacobian',r['momentum_rank']==rank and r['constraint_rank']==rank+1)
        rows=list(range(rank+1));cols=[8] if rank==0 else [6,11,10,9][:rank+1]
        minor=s.factor(r['C'].extract(rows,cols).det());check(name+' regular nonzero minor',minor!=0)
        regular.append(dict(name=name,momentum_rank=rank,metric_configuration_dof=5-rank,substitutions=r['substitutions'],constraint_jacobian=r['C'],minor_rows=rows,minor_columns=cols,minor_determinant=minor))
    homogeneous_data=dict(metric=z['g'],canonical_momentum=z['pi'],x=x,p=p,P=z['P'],T=z['T'],F=z['F'],J=z['J'],regular_representatives=regular)
    routes={}
    for k in [0,1]:
        r=constraints.scalar_routes(k);check('scalar route agreement '+str(k),r['full']['physical_pairs']==r['reduced']['physical_pairs'])
        routes[str(k)]={key:r[key] for key in ['exact','eliminated','linear','quadratic','reduced_quadratic']};routes[str(k)].update(full=closure_summary(r['full']),reduced=closure_summary(r['reduced']))
    old=constraints.legacy_trace();check('legacy correction resolves discrepancy',old['full']['physical_pairs']==old['metric']['physical_pairs']==1 and old['naive_pairs']==0)
    routes['legacy_trace']={key:old[key] for key in ['first_order','integrated','reduced','naive_pairs']};routes['legacy_trace'].update(full=closure_summary(old['full']),reduced_closure=closure_summary(old['metric']))
    modes={}
    for k in [0,1,2]:
        r=torus.mode_closure(k);check('Fourier mode count '+str(k),r['physical_pairs']==(6 if k==0 else 2));modes[str(k)]=closure_summary(r)
    k=s.Symbol('k',positive=True);r=torus.mode_closure(k);check('symbolic nonzero Fourier closure',r['physical_pairs']==2 and r['first_class']==8 and r['second_class']==0);modes['nonzero_symbolic']=closure_summary(r)
    Lt,Kt=torus.transverse_action();check('TT kinetic positive',Kt==s.eye(2));modes.update(TT_action=Lt,TT_kinetic=Kt,compact_examples=torus.compact_benchmarks())
    geom=geometric.rows();paths=geometric.paths()
    for row in geom:check('tracefree Ricci '+row['name'],s.trace(row['X'])==0)
    for row in paths:
        check('volume preserving '+row['name'],s.simplify(row['metric'].det())==1)
        m=row['orthonormal_structure'];W=row['lc_rotation_coefficients']
        residual=[s.simplify(sum(m[i,z]*s.LeviCivita(z,j,k)+s.LeviCivita(z,i,k)*W[z,j]-s.LeviCivita(z,i,j)*W[z,k] for z in range(3))) for i in range(3) for j in range(3) for k in range(j+1,3)]
        check('finite-path LC torsion '+row['name'],all(a==0 for a in residual))
    nu=s.Symbol('nu',nonzero=True,real=True);rot=s.diag(-nu,-nu,0)
    check('rotating flat frame curvature',homogeneous.curvature(rot,s.eye(3))==0);check('rotating flat frame Ricci',homogeneous.shape_residual(rot)==s.zeros(3))
    geometry_data=dict(representatives=geom,paths=paths,rotating_frame=dict(n=rot,rank=2,constant_closed_dimension=1,torus_b1=3),spectra=geometric.spectra())
    sample=[-0.8+i*0.02 for i in range(81)];zz=s.Symbol('s',real=True)
    figure_data=dict(counts=[dict(name=row['name'],closed=row['constant_closed_dimension'],regular_pairs=regular[i]['metric_configuration_dof']) for i,row in enumerate(geom)],
        paths=[dict(name=row['name'],s=sample,R=[float(row['R'].subs(zz,a)) for a in sample]) for row in paths],
        sphere=[dict(l=l,eigenvalue=l*(l+2)-2,multiplicity=2*(l-1)*(l+3)) for l in range(2,9)],
        nil=[dict(k=k,full=float(k*k/(4*s.pi**2)+k/(2*s.pi)),horizontal=float(k/(2*s.pi))) for k in range(1,9)],
        charges=[dict(energy=e,sigma=[-2+i*0.02 for i in range(201)],h=[float(s.sqrt(((-2+i*0.02)**2+e/4)/6)) for i in range(201)]) for e in [0,1,4]])
    return dict(action=action,homogeneous=homogeneous_data,routes=routes,torus=modes,geometry=geometry_data,figures=figure_data,checks=checks)
def recompute(output):
    output=Path(output)
    if output.exists():raise FileExistsError('Choose a new output directory')
    data=derive();output.mkdir(parents=True)
    hashes={}
    for name,obj in data.items():
        raw=(json.dumps(serial(obj),ensure_ascii=False,indent=2,sort_keys=True)+'\n').encode('utf-8');(output/(name+'.json')).write_bytes(raw);hashes[name+'.json']=hashlib.sha256(raw).hexdigest()
    (output/'manifest.json').write_text(json.dumps(dict(schema=1,description='Derived formulas; geometric spectra are distinct from physical mode spectra',sha256=hashes),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    return len(data['checks'])
