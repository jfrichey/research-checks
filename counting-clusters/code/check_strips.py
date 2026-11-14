"""Independent site-percolation strip calculations; standard library only."""
from collections import Counter
from fractions import Fraction as F
import json
import time
import argparse
from pathlib import Path


def step(state, mask, triangular=False, periodic=False):
    w = len(state)
    parent = list(range(2*w))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    def join(i, j):
        a, b = root(i), root(j)
        if a == b:
            return 0
        parent[a] = b
        return 1
    seen = {}
    for i, a in enumerate(state):
        if a:
            if a in seen:
                join(i, seen[a])
            else:
                seen[a] = i
    delta = mask.bit_count()
    for i in range(w):
        if mask >> i & 1:
            if state[i]:
                delta -= join(i, w+i)
            if i and (mask >> (i-1) & 1):
                delta -= join(w+i-1, w+i)
            if triangular and i and state[i-1]:
                delta -= join(i-1, w+i)
    if periodic and w > 2 and (mask & 1) and (mask >> (w-1) & 1):
        delta -= join(w, 2*w-1)
    labels = {}
    out = []
    for i in range(w):
        if mask >> i & 1:
            r = root(w+i)
            if r not in labels:
                labels[r] = len(labels)+1
            out.append(labels[r])
        else:
            out.append(0)
    return tuple(out), delta


def chain(k, wired=False, triangular=False, periodic=False):
    w = k+int(wired)
    masks = [(x<<1)|1 for x in range(1<<k)] if wired else list(range(1<<k))
    initial = (1,)+(0,)*k if wired else (0,)*k
    states, index, transitions, rewards = [initial], {initial:0}, [], []
    for state in states:
        row, reward = Counter(), 0
        for mask in masks:
            new, delta = step(state, mask, triangular, periodic)
            if new not in index:
                index[new] = len(states)
                states.append(new)
            row[index[new]] += 1
            reward += delta
        transitions.append(row)
        rewards.append(reward)
    return states, transitions, rewards, len(masks)


def solve_exact(a, b):
    """Gaussian elimination with exact rational arithmetic."""
    n = len(b)
    a = [[F(v) for v in row]+[F(rhs)] for row,rhs in zip(a,b)]
    for j in range(n):
        pivot = next(i for i in range(j,n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        for i in range(j+1,n):
            if a[i][j]:
                q = a[i][j]/a[j][j]
                a[i][j] = F(0)
                for h in range(j+1,n+1):
                    a[i][h] -= q*a[j][h]
    x = [F(0)]*n
    for i in range(n-1,-1,-1):
        x[i] = (a[i][-1]-sum(a[i][j]*x[j] for j in range(i+1,n)))/a[i][i]
    return x


def slope_exact(data):
    states, trans, rewards, den = data
    n = len(states)
    a = [[0]*n for _ in range(n)]
    for i,row in enumerate(trans):
        a[i][i] -= den
        for j,v in row.items():
            a[j][i] += v
    a[-1] = [1]*n
    pi = solve_exact(a,[0]*(n-1)+[1])
    assert sum(pi) == 1 and all(x >= 0 for x in pi)
    for j in range(n):
        assert sum(pi[i]*trans[i].get(j,0) for i in range(n)) == den*pi[j]
    return sum(p*r for p,r in zip(pi,rewards))/den


def finite(data, m, wired=False):
    states, trans, rewards, den = data
    counts = [0]*len(states)
    counts[0] = 1
    total, denom = int(wired), 1
    result = []
    for _ in range(m):
        total = total*den+sum(c*r for c,r in zip(counts,rewards))
        new = [0]*len(states)
        for i,c in enumerate(counts):
            for j,v in trans[i].items():
                new[j] += c*v
        counts = new
        denom *= den
        result.append(F(total,denom))
    return result


def brute(m,k,wired=False,triangular=False,torus=False):
    w=k+int(wired)
    total=0
    for mask in range(1<<(m*k)):
        alive={(y,x) for y in range(m) for x in range(w)
               if (wired and x==0) or (mask >> (y*k+x-int(wired)) & 1)}
        clusters=0
        while alive:
            clusters+=1
            stack=[alive.pop()]
            while stack:
                y,x=stack.pop()
                nb=[(y+1,x),(y-1,x),(y,x+1),(y,x-1)]
                if triangular:
                    nb.extend([(y+1,x+1),(y-1,x-1)])
                if torus:
                    nb=[(a%m,b%w) for a,b in nb]
                for q in nb:
                    if q in alive:
                        alive.remove(q)
                        stack.append(q)
        total+=clusters
    return F(total,1<<(m*k))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-width',type=int,default=5,
                        help='Largest square strip width; width 7 takes a few minutes.')
    parser.add_argument('--output',type=Path,help='Optional exact result JSON file.')
    args=parser.parse_args()
    output={}
    for triangular in (False,True):
        label='triangular' if triangular else 'square'
        for wired in (False,True):
            kind='wired' if wired else 'free'
            for k in range(1,args.max_width+1 if not triangular else 4):
                start=time.time()
                data=chain(k,wired,triangular)
                values=finite(data,6,wired)
                if k<=3:
                    for m in range(1,4):
                        assert values[m-1]==brute(m,k,wired,triangular)
                slope=slope_exact(data)
                key=f'{label}_{kind}_{k}'
                output[key]={'states':len(data[0]),'slope':str(slope),'decimal':float(slope),
                             'finite_1_to_6':list(map(str,values))}
                print(key,output[key],f'{time.time()-start:.2f}s',flush=True)
                if args.output:
                    args.output.parent.mkdir(parents=True,exist_ok=True)
                    args.output.write_text(json.dumps(output,indent=2)+'\n')
    for k in (3,4):
        slope=slope_exact(chain(k,periodic=True))
        print('square_cylinder',k,str(slope),float(slope),flush=True)
