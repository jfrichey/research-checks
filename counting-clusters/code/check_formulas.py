"""Certify the small-strip closed forms for every positive length.

The augmented transfer matrix includes the expected cumulative cluster count.
An exact Krylov annihilation check plus Cayley-Hamilton proves the recurrence
for all lengths, after which finitely many initial values identify the formula.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from check_strips import chain, finite


def certify(k, wired, triangular, formula, q=None):
    states, transition, rewards, den = chain(k, wired, triangular)
    size = len(states)
    dimension = size + 1

    def apply(v):
        return [sum(F(count,den)*v[j] for j,count in row.items())
                + F(rewards[i],den)*v[-1]
                for i,row in enumerate(transition)] + [v[-1]]

    def shifted(v, scalar):
        return [a-scalar*b for a,b in zip(apply(v),v)]

    h = [F(0)]*size + [F(1)]
    roots = [F(1),F(1)] + ([] if q is None else [q])
    for root in roots:
        h = shifted(h,root)
    # Initial row vector is (unit mass at reset state, initial cluster count).
    for m in range(1,dimension+1):
        h = apply(h)
        assert h[0] + int(wired)*h[-1] == 0, (k,wired,triangular,m)
    # Cayley-Hamilton extends these zero residuals to all m >= 1.
    degree = len(roots)
    values = finite((states,transition,rewards,den),degree,wired)
    assert values == [formula(m) for m in range(1,degree+1)]
    return {'random_width':k,'wired':wired,'triangular':triangular,
            'states':size,'augmented_dimension':dimension,
            'annihilator_roots':list(map(str,roots)),
            'krylov_checks':dimension,'initial_values':list(map(str,values)),
            'verified_for_all_positive_lengths':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    specs=[
        ('square_free_1',1,False,False,lambda m:F(m+1,4),None),
        ('square_free_2',2,False,False,lambda m:F(5*m+7,16),None),
        ('square_free_3',3,False,False,lambda m:(1183*m+1945+F(8)**(2-m))/3136,F(1,8)),
        ('square_wired_1',1,True,False,lambda m:F(1),None),
        ('square_wired_2',2,True,False,lambda m:(3*m+41+F(4)**(1-m))/36,F(1,4)),
        ('triangular_free_1',1,False,True,lambda m:F(m+1,4),None),
        ('triangular_free_2',2,False,True,lambda m:F(m+2,4),None),
        ('triangular_free_3',3,False,True,lambda m:(203*m+573+F(8)**(2-m))/784,F(1,8)),
        ('triangular_wired_1',1,True,True,lambda m:F(1),None),
        ('triangular_wired_2',2,True,True,lambda m:(3*m+83+F(4)**(2-m))/72,F(1,4)),
    ]
    output={key:certify(k,wired,triangular,formula,q)
            for key,k,wired,triangular,formula,q in specs}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(output,indent=2)+'\n')
    print(f'Certified {len(output)} closed forms for every positive length.')


if __name__=='__main__':
    main()
