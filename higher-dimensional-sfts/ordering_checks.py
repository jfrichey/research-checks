"""Exact finite-box tests for single-pattern overlap ordering.

Research code, developed by ChatGPT from Jacob Richey's research direction.
Integer calculations are checked against direct configuration enumeration on
small boxes. Finite-volume results do not certify infinite-volume ordering.
Only the Python standard library is required for the exact calculations.
"""

from collections import Counter
from itertools import product
from math import comb
import argparse
import json
import random


def rectangle(height, width):
    return tuple(product(range(height), range(width)))


def displacements(shape):
    return tuple(sorted({(x-a, y-b) for x,y in shape for a,b in shape
                         if (x-a, y-b) > (0, 0)}))


def agreement(shape, word):
    letters = dict(zip(shape, word))
    result = 0
    for k, (dx, dy) in enumerate(displacements(shape)):
        if all(letters[x,y] == letters[x-dx,y-dy] for x,y in shape
               if (x-dx,y-dy) in letters):
            result |= 1 << k
    return result


def classes(shape, q=2):
    groups = {}
    for word in product(range(q), repeat=len(shape)):
        mask = agreement(shape, word)
        groups.setdefault(mask, []).append(word)
    return groups


def finite_moments(shape, anchors, agree, q=2):
    """Return unnormalised binomial moments on the union of all translates."""
    sites = sorted({(x+a,y+b) for x,y in shape for a,b in anchors})
    indices = {z:i for i,z in enumerate(sites)}
    supports = [sum(1 << indices[x+a,y+b] for x,y in shape) for a,b in anchors]
    disp = {d:i for i,d in enumerate(displacements(shape))}
    incompatible = []
    for j, (x,y) in enumerate(anchors):
        bad = 0
        for i, (a,b) in enumerate(anchors[:j]):
            delta = (x-a,y-b)
            delta = max(delta, (-delta[0],-delta[1]))
            k = disp.get(delta)
            if k is not None and not (agree >> k) & 1:
                bad |= 1 << i
        incompatible.append(bad)
    powers = [q**k for k in range(len(sites)+1)]
    moments = [0]*(len(anchors)+1)

    def visit(next_index, selected, union, count):
        moments[count] += powers[len(sites)-union.bit_count()]
        for j in range(next_index,len(anchors)):
            if not selected & incompatible[j]:
                visit(j+1, selected | (1 << j), union | supports[j], count+1)
    visit(0,0,0,0)
    return moments


def distribution(moments):
    return [sum((-1)**(k-r)*comb(k,r)*moments[k]
                for k in range(r,len(moments))) for r in range(len(moments))]


def direct_distribution(shape, word, anchors, q=2):
    sites = sorted({(x+a,y+b) for x,y in shape for a,b in anchors})
    indices = {z:i for i,z in enumerate(sites)}
    copies = [[indices[x+a,y+b] for x,y in shape] for a,b in anchors]
    counts = [0]*(len(anchors)+1)
    for x in product(range(q), repeat=len(sites)):
        n = sum(all(x[i] == letter for i,letter in zip(copy,word)) for copy in copies)
        counts[n] += 1
    return counts


def hinge_gaps(u,v):
    """Integer multiples of E[(N_v-k)+]-E[(N_u-k)+]."""
    return [sum((r-k)*(v[r]-u[r]) for r in range(k+1,len(u)))
            for k in range(len(u)-1)]


def word_string(word, width):
    return ';'.join(''.join(map(str,word[i:i+width]))
                    for i in range(0,len(word),width))


def exact_scan(height,width,anchor_height,anchor_width,q=2):
    shape = rectangle(height,width)
    anchors = rectangle(anchor_height,anchor_width)
    groups = classes(shape,q)
    laws = {g:distribution(finite_moments(shape,anchors,g,q)) for g in groups}
    for law in laws.values():
        assert all(x >= 0 for x in law)
        assert sum(law) == q**((height+anchor_height-1)*(width+anchor_width-1))
    violations = []
    pairs = 0
    for u,ug in groups.items():
        for v,vg in groups.items():
            if u != v and u & v == u:
                pairs += 1
                gaps = hinge_gaps(laws[u],laws[v])
                if min(gaps) < 0:
                    violations.append(dict(u=word_string(ug[0],width),
                        v=word_string(vg[0],width),mask_u=u,mask_v=v,
                        hinge_gaps=gaps,law_u=laws[u],law_v=laws[v]))
    return dict(shape=[height,width],anchor_shape=[anchor_height,anchor_width],
                q=q,agreement_classes=len(groups),inclusion_pairs=pairs,
                convex_order_failures=violations)


def irregular_scan(seed=20261020, trials=500):
    rng = random.Random(seed)
    total_pairs = 0
    for trial in range(trials):
        shape = tuple(sorted(rng.sample(rectangle(3,3),rng.randrange(3,9))))
        anchors = tuple(sorted(rng.sample(rectangle(4,4),rng.randrange(5,13))))
        groups = classes(shape)
        laws = {g:distribution(finite_moments(shape,anchors,g)) for g in groups}
        for u,ug in groups.items():
            for v,vg in groups.items():
                if u != v and u & v == u:
                    total_pairs += 1
                    gaps = hinge_gaps(laws[u],laws[v])
                    if min(gaps) < 0:
                        return dict(seed=seed,trial=trial,shape=shape,anchors=anchors,
                                    u=ug[0],v=vg[0],hinge_gaps=gaps,
                                    law_u=laws[u],law_v=laws[v])
    return dict(seed=seed,trials=trials,inclusion_pairs=total_pairs,failures=[])


def smoke():
    shape = rectangle(2,2)
    anchors = rectangle(2,2)
    for g,words in classes(shape).items():
        assert distribution(finite_moments(shape,anchors,g)) == direct_distribution(shape,words[0],anchors)
    shape = rectangle(2,3)
    anchors = rectangle(2,2)
    for g,words in classes(shape).items():
        assert distribution(finite_moments(shape,anchors,g)) == direct_distribution(shape,words[0],anchors)
    print('Exact marked-moment versus direct configuration checks passed.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--shape', nargs=2, type=int, default=[2,3])
    parser.add_argument('--anchors', nargs=2, type=int, default=[3,3])
    parser.add_argument('--q', type=int, default=2)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--irregular',type=int,default=0)
    args = parser.parse_args()
    if args.smoke:
        smoke()
    result = irregular_scan(trials=args.irregular) if args.irregular else exact_scan(*args.shape,*args.anchors,args.q)
    print(json.dumps(result,indent=2), flush=True)


def binary_representatives(height,width):
    """One representative per agreement set, using bitwise rectangle overlaps."""
    shape=rectangle(height,width)
    specs=[]
    for dx,dy in displacements(shape):
        lower=sum(1 << (r*width+c) for r,c in shape
                  if 0<=r+dx<height and 0<=c+dy<width)
        specs.append((dx*width+dy,lower))
    groups={}
    for bits in range(1 << (height*width)):
        mask=0
        for i,(shift,lower) in enumerate(specs):
            if ((bits ^ (bits >> shift)) & lower)==0:mask |= 1 << i
        if mask not in groups:
            groups[mask]=tuple((bits >> j)&1 for j in range(height*width))
    return groups
