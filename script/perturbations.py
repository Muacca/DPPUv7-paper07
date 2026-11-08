"""Recompute the paper's specified perturbation models into a new directory."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
MODELS=Path(__file__).resolve().parent/'paper07_models'
COMMANDS={
 'flat':[['integrability.py'],['standing_wave.py'],['two_waves.py'],['connection_modes.py'],['lorentz_boundary.py'],['connection_controls.py']],
 'backgrounds':[['vacuum_backgrounds.py'],['fixed_quotients.py'],['sol_gowdy.py'],['sol_gowdy_controls.py'],['sol_gowdy_example.py'],['connection.py','Sol'],['endpoints.py','Sol'],['nil_base.py'],['nil_base_controls.py'],['connection.py','Nil'],['endpoints.py','Nil'],['nil_base_example.py']],
 'taub':[['harmonics.py']]+[cmd for j in ['0','1/2','1','3/2'] for cmd in [['operators.py',j],['reduce_block.py',j]]]+[['einstein_check.py','3/2'],['connection_block.py'],['kid_controls.py'],['finite_example.py']],
 'nil':[['basis.py'],['quotient_checks.py'],['operators.py','7'],['operators.py','8'],['reduce_blocks.py'],['einstein_check.py'],['cartan_grading.py'],['finite_example.py']],
 'sol':[['basis.py'],['operators.py'],['einstein_check.py'],['reduce_orbit.py'],['reduction_core.py'],['cartan_orbit.py']],
 'lorenz':[['expand.py'],['principal.py'],['wave_operator.py'],['initial_maps.py'],['coefficient_jets.py']]+[
   cmd for n,w,name,tol,q in [(1023,5.25,'harmonic_N1023_W525_fine.json','3e-14',1536),(1535,5.25,'harmonic_N1535_W525_fine.json','3e-14',2048),(1663,5.6875,'harmonic_N1663_W56875_fine.json','3e-14',2304),(1535,5.25,'harmonic_N1535_W525.json','1e-12',2048)]
   for cmd in [['evolve.py',str(n),str(w),name,tol],['residuals.py',name,str(q)]]]+[['compare.py'],['dense_residual.py']]
}

def run_model(model,out,canonical=None):
    out.mkdir(parents=True)
    if model=='lorenz':
        inputs=out/'inputs';inputs.mkdir()
        source=canonical or ROOT/'data/perturbations/sol/reduction_core.json'
        if not source.is_file():raise FileNotFoundError('Recompute sol first and pass --canonical-input PATH_TO_reduction_core.json')
        shutil.copyfile(source,inputs/'reduction_core.json')
        data=ROOT/'data/initial';manifest=json.loads((data/'manifest.json').read_text())
        assert hashlib.sha256((data/'sol_profile.json').read_bytes()).hexdigest()==manifest['sha256']['sol_profile.json']
        shutil.copyfile(data/'sol_profile.json',inputs/'wave_initial.json')
    env=os.environ.copy();env.update(PAPER07_OUTPUT=str(out.resolve()),SYMPY_GROUND_TYPES='python',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONNOUSERSITE='1',PYTHONHASHSEED='1')
    env.pop('PYTHONPATH',None)
    for args in COMMANDS[model]:
        print(model,*args,flush=True)
        completed=subprocess.run([sys.executable,'-B',str(MODELS/model/args[0]),*args[1:]],env=env,cwd=MODELS/model)
        if completed.returncode:raise RuntimeError('Producer failed: '+' '.join(args))
        # The reduced fourth-order system (SG-05) is an explicit input to its residual audit.
        if model=='lorenz' and args[0]=='expand.py':shutil.copyfile(out/'expanded_Hr.json',out/'inputs/expanded_Hr.json')
    hashes={p.relative_to(out).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*.json')}
    (out/'manifest.json').write_text(json.dumps(dict(model=model,commands=COMMANDS[model],sha256=hashes),indent=2)+'\n',encoding='utf8')

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('recompute');q.add_argument('--model',choices=['all',*COMMANDS],required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--canonical-input',type=Path)
    q=sub.add_parser('export');q.add_argument('--input',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('render');q.add_argument('--data',type=Path,default=ROOT/'data/perturbations');q.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve()
    if out.exists():raise FileExistsError('Choose a new output directory; saved results are never overwritten')
    if a.command=='render':
        from paper07_models.saved import render
        print('Rendered:',', '.join(render(a.data,out)));return
    if a.command=='export':
        from paper07_models.saved import export
        print('Exported files:',export(a.input,out));return
    if a.model=='all':
        out.mkdir()
        for model in COMMANDS:run_model(model,out/model,out/'sol/reduction_core.json' if model=='lorenz' else None)
    else:run_model(a.model,out,a.canonical_input.resolve() if a.canonical_input else None)
if __name__=='__main__':main()
