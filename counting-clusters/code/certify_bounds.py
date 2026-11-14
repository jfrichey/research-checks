"""Reproduce certified strip slopes and bounds, with compilation in /tmp."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import tempfile


def slope_interval(d):
    estimate=F(int(d['expectation_numerator']),int(d['Q'])*d['D'])
    error=F(d['reward_bound']*int(d['residual_numerator']),int(d['Q']))
    return estimate-error,estimate+error


def outward(x,places,up=False):
    scale=10**places
    z=x*scale
    n=-((-z.numerator)//z.denominator) if up else z.numerator//z.denominator
    return f'{n//scale}.{n%scale:0{places}d}'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-width',type=int,default=10)
    parser.add_argument('--output-dir',type=Path,
                        help='Optional destination for JSON and dyadic certificates.')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    known=json.loads((root/'data/exact-strips.json').read_text())
    results={}
    if args.output_dir:
        (args.output_dir/'certificates').mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='counting-clusters-') as tmp:
        exe=Path(tmp)/'certify_strips'
        subprocess.run(['g++','-O3','-std=c++17',str(root/'code/certify_strips.cpp'),
                        '-o',str(exe)],check=True)
        print('k  certified lower bound  certified upper bound',flush=True)
        for k in range(1,args.max_width+1):
            for wired in (False,True):
                kind='wired' if wired else 'free'
                cmd=[str(exe),str(k),str(int(wired))]
                if args.output_dir:
                    cmd.append(str(args.output_dir/'certificates'/f'certificate-{kind}-{k}.txt'))
                cp=subprocess.run(cmd,capture_output=True,text=True,check=True)
                d=json.loads(cp.stdout)
                lo,hi=slope_interval(d)
                d['interval_lower']=str(lo)
                d['interval_upper']=str(hi)
                key=f'square_{kind}_{k}'
                if key in known:
                    exact=F(known[key]['slope'])
                    assert lo<=exact<=hi,(key,lo,exact,hi)
                results[key]=d
            if k>=2:
                lo=slope_interval(results[f'square_free_{k}'])[0]-slope_interval(results[f'square_free_{k-1}'])[1]
                hi=slope_interval(results[f'square_wired_{k}'])[1]-slope_interval(results[f'square_wired_{k-1}'])[0]
                print(k,outward(lo,8),outward(hi,8,up=True),flush=True)
            if args.output_dir:
                (args.output_dir/'certified-strips.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    main()
