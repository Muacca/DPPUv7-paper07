"""Paper07 formula recomputation and saved-data figure rendering."""
import argparse
from pathlib import Path
def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('recompute');r.add_argument('--output',type=Path,required=True)
    r=sub.add_parser('render');r.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'data/revised');r.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.command=='recompute':
        from paper07_results.derive import recompute
        print('Exact checks passed:',recompute(args.output))
    else:
        from paper07_results.render import render
        print('Rendered:',', '.join(render(args.data,args.output)))
if __name__=='__main__':main()
