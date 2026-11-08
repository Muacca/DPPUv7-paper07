"""Render supplied numeric data; this command does not run symbolic derivations."""
from pathlib import Path
import json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
COLORS=['#235b8f','#d0782a','#56836b','#785b8c']
def render(data,output):
    data=Path(data);output=Path(output)
    if output.exists():raise FileExistsError('Choose a new figure directory')
    manifest=json.loads((data/'manifest.json').read_text(encoding='utf-8'))
    for name,digest in manifest['sha256'].items():
        if hashlib.sha256((data/name).read_bytes()).hexdigest()!=digest:raise ValueError('Data hash mismatch: '+name)
    z=json.loads((data/'figures.json').read_text(encoding='utf-8'));output.mkdir(parents=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.16,'savefig.dpi':180})
    names=[]
    def save(fig,name):
        fig.savefig(output/(name+'.png'),bbox_inches='tight',metadata={'Software':'Matplotlib'})
        fig.savefig(output/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None,'Creator':'Paper07 supplement'})
        plt.close(fig);names.append(name)
    fig,axs=plt.subplots(1,2,figsize=(8.3,3.1),layout='constrained');labels=[r['name'] for r in z['counts']]
    for ax,key,title in zip(axs,['closed','regular_pairs'],['Fixed-coframe closed one-forms','Regular homogeneous metric pairs']):
        vals=[r[key] for r in z['counts']];ax.bar(labels,vals,color=COLORS);ax.set_ylim(0,6);ax.set_yticks(range(7));ax.set_title(title,fontsize=11)
        for i,v in enumerate(vals):ax.text(i,v+.12,str(v),ha='center')
    save(fig,'fig01_distinct_counts')
    fig,ax=plt.subplots(figsize=(6.6,3.8),layout='constrained')
    for i,row in enumerate(z['charges']):
        ax.plot(row['sigma'],row['h'],color=COLORS[i],label=f"E = {row['energy']}");ax.plot(row['sigma'],-np.array(row['h']),color=COLORS[i])
    ax.set_xlabel(r'Shear amplitude $\sigma$ (one unit-norm direction)');ax.set_ylabel(r'Mean expansion $h$');ax.legend(frameon=False,ncol=3,loc='upper center');ax.set_ylim(-1.2,1.35);ax.set_title(r'Compact energy condition: $6h^2=\sigma^2+E/4$',fontsize=11)
    save(fig,'fig02_compact_constraint')
    fig,axs=plt.subplots(1,2,figsize=(8.3,3.2),layout='constrained')
    for row in z['paths']:
        if row['name']=='nil_flat':axs[0].plot(row['s'],row['R'],label='Nil horizontal shape',color=COLORS[0])
        if row['name']=='nil_vertical':axs[0].plot(row['s'],row['R'],label='Nil vertical shape',color=COLORS[1])
        if row['name']=='sol_base':axs[1].plot(row['s'],row['R'],label='Sol base shape',color=COLORS[1])
        if row['name']=='sol_isometry':axs[1].plot(row['s'],row['R'],label='Sol isometric path',color=COLORS[0])
    for ax in axs:ax.set_xlabel('Finite path parameter s');ax.set_ylabel(r'$q_0^2\,{}^{(3)}R$');ax.legend(frameon=False,fontsize=9)
    save(fig,'fig03_finite_paths')
    fig,axs=plt.subplots(1,2,figsize=(8.3,3.2),layout='constrained')
    axs[0].plot([r['l'] for r in z['sphere']],[r['eigenvalue'] for r in z['sphere']],'o-',color=COLORS[0]);axs[0].set_xlabel('Tensor harmonic label l');axs[0].set_ylabel(r'$q^2\lambda_{\rm rough}$');axs[0].set_title('Sphere: spatial TT rough Laplacian',fontsize=10)
    for key,label,color in [('full','Full scalar eigenvalue',COLORS[0]),('horizontal','Horizontal contribution',COLORS[1])]:axs[1].plot([r['k'] for r in z['nil']],[r[key] for r in z['nil']],'o-',label=label,color=color)
    axs[1].set_xlabel('Central label k (j = 0)');axs[1].set_ylabel(r'$q^2\lambda_{\rm scalar}$');axs[1].legend(frameon=False,fontsize=8);axs[1].set_title('Nil: fixed-lattice scalar Laplacian',fontsize=10)
    save(fig,'fig04_geometric_operators')
    return names
