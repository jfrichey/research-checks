"""Exact finite decorated-drain inventories; see drain-inventory.md.

Standard library only. Count coefficients include unused sleep coins along
the completed drains. A positive census is not a proof for all inventories.
Research implementation written by Codex; no independent user validation.
"""

from verify_drained_paths import fixed_path_occupation, execute
from fractions import Fraction as F
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import product
from math import comb, factorial
from random import Random
import json

@lru_cache(None)
def binomial_poly(counts):
    return {ks: prod_int(comb(m,k) for m,k in zip(counts,ks)) for ks in product(*(range(m+1) for m in counts))}

def prod_int(values):
    a=1
    for v in values:a*=v
    return a

def add_counts(poly,other):
    out=defaultdict(int)
    for a,c in poly.items():
        for b,d in other.items():out[tuple(x+y for x,y in zip(a,b))]+=c*d
    return dict(out)

@lru_cache(None)
def occupation_poly(paths,n,target=0):
    @lru_cache(None)
    def rec(positions,queue,residents):
        if not queue:
            if residents[target]<0:return {}
            unused=Counter(z for path,i in zip(paths,positions) for z in path[i+1:-1])
            return binomial_poly(tuple(unused[z] for z in range(n)))
        label=queue[0];z=paths[label][positions[label]]
        pending=list(queue[1:]);indices=list(positions);rest=list(residents)
        def emit(label,indices,pending):
            indices[label]+=1
            if paths[label][indices[label]]!=-1:pending.append(label)
        if rest[z]>=0:emit(rest[z],indices,pending);rest[z]=-1
        sleeping=list(rest);sleeping[z]=label
        yes=rec(tuple(indices),tuple(pending),tuple(sleeping))
        emit(label,indices,pending)
        no=rec(tuple(indices),tuple(pending),tuple(rest))
        out=defaultdict(int,no)
        for key,count in yes.items():
            k=list(key);k[z]+=1;out[tuple(k)]+=count
        return dict(out)
    value=rec((0,)*len(paths),tuple(range(len(paths))),(-1,)*n)
    m=sum(len(path)-1 for path in paths)
    assert F(sum(value.values()),2**m)==fixed_path_occupation(paths,(F(1,2),)*n,target)
    return value

class TooLarge(Exception):pass

def inventory_poly(starts,stock,split,n,limit=5000):
    left=dict(stock);count=0;poly=defaultdict(int)
    def rec(label,path,finished):
        nonlocal count
        z=path[-1]
        if z==-1:
            paths=finished+(tuple(path),)
            if label+1==len(starts):
                if any(left.values()):return
                count+=1
                if count>limit:raise TooLarge
                a,b=paths[:split],paths[split:]
                na=Counter(z for path in a for z in path[:-1]);nb=Counter(z for path in b for z in path[:-1])
                plus=add_counts(occupation_poly(b,n),binomial_poly(tuple(na[z] for z in range(n))))
                minus=add_counts(occupation_poly(a,n),binomial_poly(tuple(nb[z] for z in range(n))))
                for key,value in plus.items():poly[key]+=value
                for key,value in minus.items():poly[key]-=value
                return
            rec(label+1,[starts[label+1]],paths);return
        for (site,d),remaining in list(left.items()):
            if site==z and remaining:
                left[site,d]-=1;rec(label,path+[d],finished);left[site,d]+=1
    rec(0,[starts[0]],())
    return count,{key:value for key,value in poly.items() if value}

def determinant(matrix):
    """Exact determinant, including the empty determinant 1."""
    a=[list(map(F,row)) for row in matrix];value=F(1);n=len(a)
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[j],a[pivot]=a[pivot],a[j];value=-value
        entry=a[j][j];value*=entry
        for i in range(j+1,n):
            factor=a[i][j]/entry
            for k in range(j+1,n):a[i][k]-=factor*a[j][k]
    assert value.denominator==1
    return value.numerator

def path_family_count(stock,n):
    """Last-exit determinant count; starting multiplicities are div(stock)."""
    out=[0]*n;incoming=[0]*n
    for (z,w),a in stock.items():
        assert a>=0
        out[z]+=a
        if w>=0:incoming[w]+=a
    if any(out[z]<incoming[z] for z in range(n)):return 0
    active=[z for z in range(n) if out[z]];idx={z:i for i,z in enumerate(active)}
    L=[[int(z==w)*out[z] for w in active] for z in active]
    for (z,w),a in stock.items():
        if a and w>=0:
            if w not in idx:return 0
            L[idx[z]][idx[w]]-=a
    value=F(determinant(L))
    for z in active:value*=factorial(out[z]-1)
    for a in stock.values():value/=factorial(a)
    assert value.denominator==1,(stock,value)
    return value.numerator


def forest_inventory_poly(starts,stock,split,n,target=0,guard=500000):
    """Count the same pool by ordinary arrivals and last-exit forest closure.

    The guard is a computational limit, never a history truncation. A case
    exceeding it raises TooLarge and yields no claimed coefficient result.
    """
    edges=tuple(sorted(stock));initial=tuple(stock[e] for e in edges)
    divergence=Counter()
    for (z,w),a in stock.items():
        divergence[z]+=a
        if w>=0:divergence[w]-=a
    assert all(divergence[z]==Counter(starts)[z] for z in range(n))
    @lru_cache(None)
    def drain_count(counts):
        return path_family_count(dict(zip(edges,counts)),n)
    @lru_cache(None)
    def solve(counts,lo,hi,residents):
        nonlocal visited
        visited += 1
        if visited > guard: raise TooLarge
        if not lo and not hi:
            delta=int(residents[n+target]!=-2)-int(residents[target]!=-2)
            if not delta:return {}
            ways=drain_count(counts)
            if not ways:return {}
            unused=[0]*n
            for (z,w),a in zip(edges,counts):unused[z]+=a
            return {key:delta*ways*value for key,value in binomial_poly(tuple(unused)).items()}
        c=0 if lo else 1;queues=[lo,hi];z=queues[c][0];rest=list(queues[c][1:]);r=list(residents)
        old=r[c*n+z];r[c*n+z]=-2
        if old>=0:rest.append(old)
        result=defaultdict(int)
        for j,((site,d),a) in enumerate(zip(edges,counts)):
            if site!=z or not a:continue
            after=list(counts);after[j]-=1
            for bit in (0,1):
                r2=list(r);pending=list(rest)
                if bit:r2[c*n+z]=d
                elif d>=0:pending.append(d)
                qq=list(queues);qq[c]=tuple(pending)
                term=solve(tuple(after),qq[0],qq[1],tuple(r2))
                for key,value in term.items():
                    k=list(key);k[z]+=bit;result[tuple(k)]+=value
        return {key:value for key,value in result.items() if value}
    visited=0
    poly=solve(initial,tuple(starts[:split]),tuple(starts[split:]),(-2,)*(2*n))
    return drain_count(initial),poly,solve.cache_info().currsize

def forest_audit():
    rng=Random(8813)
    for trial in range(24):
        n=3;starts=[0,2,0,2,1];paths=[]
        for source in starts:
            path=[source]
            for _ in range(rng.randrange(0,3)):
                d=rng.choice([0,1,2,-1]);path.append(d)
                if d==-1:break
            if path[-1]!=-1:path.append(-1)
            paths.append(path)
        stock=Counter((a,b) for path in paths for a,b in zip(path,path[1:]))
        expected=inventory_poly(starts,stock,2,n)
        actual=forest_inventory_poly(starts,stock,2,n)
        assert expected==actual[:2],(trial,expected,actual)
    examples=[]
    for kind in range(8):
        n=4;common=[0,2] if kind%2 else [1,2];extra=1 if kind%2 else 3
        starts=common+common+[extra]
        stock=Counter((z,-1) for z in starts)
        cycles=([(0,1),(1,0)],[(0,2),(2,0)],[(1,2),(2,1)],[(2,3),(3,2)])
        for i in range(2+kind//2):stock.update(cycles[i%4])
        ways,poly,states=forest_inventory_poly(starts,stock,2,n)
        assert all(value>=0 for value in poly.values()),(kind,stock,poly)
        examples.append(dict(kind=kind,edges=sum(stock.values()),bundles=ways,states=states))
    return dict(check='last_exit_forest_enumeration',independent_comparisons=24,
                interacting_cycle_examples=examples,all_coefficients_nonnegative=True,
                status='Exact counting method; finite positivity evidence only')


def forest_refinement_obstruction():
    """Fixing the pooled last-exit forest is a false stronger assertion."""
    counts=defaultdict(Counter);word_counts=Counter()
    for row0 in ((1,-1,-1),(-1,1,-1),(-1,-1,1)):
        for row1 in ((2,-1),(-1,2)):
            if row0[-1]==1 and row1[-1]==2:continue
            rows=(row0,row1,(0,));forest=tuple(row[-1] for row in rows)
            word_counts[forest]+=1
            split=[0,0,0];z=0
            while z!=-1:
                d=rows[z][split[z]];split[z]+=1;z=d
            coordinates=[(site,index) for site,row in enumerate(rows)
                         for index in range(len(row))]
            arrows={((0,site,index+1) if index<split[site]
                     else (1,site,index+1-split[site])): rows[site][index]
                    for site,index in coordinates}
            for values in product((0,1),repeat=6):
                bits={((0,site,index+1) if index<split[site]
                       else (1,site,index+1-split[site])): bit
                      for (site,index),bit in zip(coordinates,values)}
                data=execute(bits,[[0],[0,1]],arrows=arrows)
                k=tuple(sum(bit for (site,index),bit in zip(coordinates,values)
                            if site==z) for z in range(3))
                counts[forest][k]+=int(0 in data['supports'][1])-int(0 in data['supports'][0])
    actual={forest:{key:value for key,value in poly.items() if value}
            for forest,poly in counts.items()}
    expected={(-1,-1,0):{(1,0,0):-1,(2,0,0):-1,(1,2,0):1,(2,2,0):1},
              (-1,2,0):{(1,0,0):1,(2,0,0):1,(1,1,0):1,(2,1,0):1},
              (1,-1,0):{}}
    assert actual==expected
    assert word_counts=={(-1,-1,0):2,(-1,2,0):2,(1,-1,0):1}
    stock={(0,-1):2,(0,1):1,(1,-1):1,(1,2):1,(2,0):1}
    ways,poly=inventory_poly([0,0,1],stock,1,3)
    total=Counter()
    for part in actual.values():total.update(part)
    assert ways==5 and poly=={key:value for key,value in total.items() if value}
    assert all(value>=0 for value in poly.values())
    return dict(check='fixed_last_exit_forest_obstruction',bundles=5,
                forest_word_counts={str(key):value for key,value in word_counts.items()},
                negative_coefficients_in_forest=[-1,-1],full_pool_nonnegative=True,
                verified=True,status='Only fixed-forest refinement refuted')


if __name__=='__main__':
    stock={(0,2):2,(2,0):2,(0,-1):3,(1,0):1}
    bundles, coefficients = inventory_poly([0,0,1],stock,1,3)
    assert bundles == 6
    assert coefficients == {(1,0,2): 1, (2,0,2): 3, (3,0,2): 3, (4,0,2): 1}
    print(json.dumps(dict(check='six_bundle_cycle_redistribution', bundles=bundles,
                         coefficients={str(key): value for key, value in coefficients.items()},
                         verified=True)), flush=True)
    rng=Random(8813);tested=too_long=too_many=0
    for trial in range(240):
        n=4;common=[0,2] if trial%2 else [1,2];extra=1 if trial%2 else 3;starts=common+common+[extra]
        paths=[]
        for source in starts:
            path=[source]
            for _ in range(rng.randrange(0,3)):
                d=rng.choice(list(range(n))+[-1]);path.append(d)
                if d==-1:break
            if path[-1]!=-1:path.append(-1)
            paths.append(path)
        stock=Counter((a,b) for path in paths for a,b in zip(path,path[1:]))
        if sum(stock.values())>11:
            too_long+=1
            continue
        try:count,poly=inventory_poly(starts,stock,2,n)
        except TooLarge:
            too_many+=1
            continue
        tested+=1
        negatives=[(key,value) for key,value in poly.items() if value<0]
        if negatives:
            raise AssertionError((trial, starts, count, negatives, dict(stock), paths))
    print(json.dumps(dict(check='decorated_drain_inventory_census', seed=8813,
                         generated_cases=240, exact_pools_checked=tested,
                         excluded_above_11_edges=too_long, excluded_above_5000_bundles=too_many,
                         all_coefficients_nonnegative=True,
                         status='Finite exact evidence; general inventory inequality open')), flush=True)

    print(json.dumps(forest_audit()),flush=True)
    print(json.dumps(forest_refinement_obstruction()),flush=True)
