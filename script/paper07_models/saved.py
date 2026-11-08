"""Select reproducible publication data from a complete model recomputation."""
from pathlib import Path
import hashlib,json,shutil

def publication_scope(dense,comparisons):
    """Theorem scope and numerical evidence are separate publication metadata."""
    return dict(
        flat_integrability='Known theorem applied under compact boundaryless CMC and sufficient regularity hypotheses; Appendix D.',
        second_order='Explicit finite-order corrections; distinct from the theorem proving local exact-family existence.',
        curved_reduction='Selected backgrounds and sectors in Appendices F-G; no universal curved physical spectrum.',
        sol_lorenz_interior_supported=True,sol_reduced_fourth_order_convergence_established=False,
        whole_line_error_guarantee=False,all_curved_backgrounds_nonlinear_lifting_sufficiency_established=False,
        nil_nonlinear_lifting_sufficiency_established=False,taub_nonlinear_lifting_sufficiency_established=False,
        sol_local_integrability=dict(status='PROVED_SCOPED',
            background='BF-01 on the fixed compact Sol quotient; fixed t0 in [1,2].',
            data='Real metric and momentum-density variations in H^s x H^(s-1), fixed integer s >= 6, satisfying all linearized vacuum constraints.',
            conclusion='Tangency to an exact vacuum constraint-data family, local Cauchy development and LC lift to source-free EC with constant complete NY.',
            basis='Fischer-Marsden-Moncrief (1980), Theorem 2.1, with the function-space and elliptic arguments of SL-01 through SL-05.',
            modes='Nonzero-base Gowdy reconstruction; E(t0)y for y in D^4 with original tensor phases and conjugate orbit restored; continuous SH-03 Hermite initial profile; infinite sums obeying SL-07.',
            boundary='Small family parameter and local time. No fixed CMC value, gauge, Gowdy ansatz or single-orbit closure is imposed. Arbitrary reduced L2 data, numerical trajectory certification, SG-05 well-posedness and nonlinear stability are not established.'),
        max_reduced_comparisons=[r['maxima']['reduced'] for r in comparisons],
        max_dense_residual={k:max(r[k] for r in dense['results']) for k in ['wave_residual_L2','Einstein_L2','lorenz_L2','ADM_constraint_L2','rho_red']},
        limitations=['Finite-time interior numerical evidence is not a whole-line theorem.','Convergence of the reduced fourth-order system (SG-05) under the Lorenz evolution is not established.','Finite low-label checks are not an all-j or non-LRS classification.','Sol local integrability concerns admissible full initial data, not a closed realization of SG-05 or fixed-sector nonlinear families.'])

def export(source,output):
    import numpy as np
    from scipy.special import roots_legendre
    source=Path(source);output=Path(output)
    if output.exists():raise FileExistsError('Choose a new export directory')
    output.mkdir(parents=True)
    selected=[]
    for group in ['flat','backgrounds','nil','sol','taub','lorenz']:
        for p in sorted((source/group).glob('*.json')):
            if p.name=='manifest.json' or p.name.startswith('operators_g'):continue
            if group=='lorenz' and (p.name.startswith('audit_') or p.name.startswith('harmonic_N')):continue
            dest=output/group/p.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
            selected.append(dest.relative_to(output).as_posix())
    def write(name,obj):
        dest=output/name;dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    # Keep the numerical diagnostics without shipping four redundant full trajectories.
    for p in sorted((source/'lorenz').glob('audit_*.json')):
        obj=json.loads(p.read_text());obj.pop('reduced_fields')
        write('lorenz/'+p.name,obj)
    raw=json.loads((source/'lorenz/audit_harmonic_N1535_W525_fine_Q2048.json').read_text())
    x,w=roots_legendre(raw['quadrature']);ix=np.r_[np.arange(0,len(x),4),len(x)-1]
    fields=np.asarray(raw['reduced_fields']);ti=[0,50,100]
    write('figures.json',dict(description='Saved samples of reconstructed reduced variables on the Sol interior; no evolution is performed by the renderer.',
        orbit=[1,0],x_relative=(4*x[ix]).tolist(),time=[raw['time'][i] for i in ti],
        Q1=fields[ti][:,ix,0].tolist(),Pi2=fields[ti][:,ix,3].tolist(),
        source='lorenz/audit_harmonic_N1535_W525_fine_Q2048.json',
        scope='One orbit, one initial profile, finite time and interior only; not a physical spectrum.'))
    dense=json.loads((source/'lorenz/harmonic_dense_residual.json').read_text())
    comparisons=json.loads((source/'lorenz/harmonic_comparisons.json').read_text())
    write('scope.json',publication_scope(dense,comparisons))
    hashes={p.relative_to(output).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.rglob('*.json'))}
    write('manifest.json',dict(schema=1,edition='v2',sha256=hashes,
        omitted='Intermediate Nil oscillator buffers and full Sol trajectories are generated by recompute; only model coefficients, controls, finite examples, diagnostics and plot samples are distributed.'))
    return len(hashes)

def render(data,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data=Path(data);output=Path(output)
    if output.exists():raise FileExistsError('Choose a new figure directory')
    manifest=json.loads((data/'manifest.json').read_text())
    for name,digest in manifest['sha256'].items():
        if hashlib.sha256((data/name).read_bytes()).hexdigest()!=digest:raise ValueError('Hash mismatch: '+name)
    z=json.loads((data/'figures.json').read_text());output.mkdir(parents=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.16,'savefig.dpi':180})
    fig,axs=plt.subplots(1,2,figsize=(8.3,3.3),layout='constrained')
    colors=['#235b8f','#d0782a','#56836b']
    for ax,key,label in zip(axs,['Q1','Pi2'],[r'$Q_1$',r'$\Pi_2$']):
        for t,row,color in zip(z['time'],z[key],colors):ax.plot(z['x_relative'],row,label=f't = {t:g}',color=color)
        ax.set_xlabel(r'$\xi-\xi_0$');ax.set_ylabel(label);ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2));ax.legend(frameon=False)
    fig.suptitle('Sol orbit (1, 0): finite-time interior example',fontsize=12)
    name='fig05_sol_interior'
    fig.savefig(output/(name+'.png'),bbox_inches='tight',metadata={'Software':'Matplotlib'})
    fig.savefig(output/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None,'Creator':'Paper07 supplement'})
    plt.close(fig)
    return [name]
