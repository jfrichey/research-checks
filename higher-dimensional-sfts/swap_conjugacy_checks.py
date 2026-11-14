"""Exact checks for dimension-independent simultaneous pattern swaps.

Research code developed by ChatGPT from Jacob Richey's research direction.
The all-dimensional proof is in entropy-ordering.md. Local tests verify a
nontrivial binary 2x3 example; graph enumeration concerns this specified
swap criterion and alphabet permutations, not all possible conjugacies.
"""
import argparse
from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path

EVIDENCE = Path(__file__).with_name('swap-conjugacy-evidence.json')


def geometry(shape):
    sites=tuple(product(*(range(n) for n in shape))); index={z:i for i,z in enumerate(sites)}
    differences=sorted({tuple(x-y for x,y in zip(a,b)) for a in sites for b in sites if a>b})
    pairs=[[(i,index[t]) for i,s in enumerate(sites)
            if (t:=tuple(x-y for x,y in zip(s,d))) in index] for d in differences]
    return sites,differences,pairs


def cross(pairs,u,v):
    return sum(1<<k for k,edges in enumerate(pairs) if all(u[i]==v[j] for i,j in edges))


def partition(word):
    labels={}
    return tuple(labels.setdefault(c,len(labels)) for c in word)


def graph_audit(shape,q):
    sites,ds,pairs=geometry(shape);groups=defaultdict(list)
    for word in product(range(q),repeat=len(sites)):
        groups[cross(pairs,word,word)].append(word)
    records=[]
    for g,words in sorted(groups.items()):
        parent=list(range(len(words)));swaps=0;permutations=0
        def root(i):
            while parent[i]!=i:
                parent[i]=parent[parent[i]];i=parent[i]
            return i
        for i,u in enumerate(words):
            for j in range(i+1,len(words)):
                v=words[j];relabel=partition(u)==partition(v)
                swap=cross(pairs,u,v)==g and cross(pairs,v,u)==g
                if swap:
                    # Independent coordinate-level check of overlap protection.
                    changed={k for k in range(len(sites)) if u[k]!=v[k]}
                    for k,edges in enumerate(pairs):
                        if g>>k&1:
                            assert all(a not in changed and b not in changed for a,b in edges)
                    swaps+=1
                if relabel:permutations+=1
                if swap or relabel:parent[root(i)]=root(j)
        components=sorted(Counter(root(i) for i in range(len(words))).values())
        records.append(dict(g=g,patterns=len(words),swap_edges=swaps,
                            relabel_edges=permutations,components=components))
    return dict(shape=list(shape),q=q,groups=records)


def local_example():
    # Flip the top middle bit of 000/011 or 010/011.
    u=(0,0,0,0,1,1);v=(0,1,0,0,1,1)
    sites,ds,pairs=geometry((2,3))
    masks=[cross(pairs,a,b) for a in (u,v) for b in (u,v)]
    assert len(set(masks))==1
    positions=list(product(range(3),range(-2,3)))
    def value(x,r,c):
        marker=(x[r,c-1]==0 and x[r,c+1]==0 and x[r+1,c-1]==0
                and x[r+1,c]==1 and x[r+1,c+1]==1)
        return x[r,c]^marker
    checked=0
    for bits in product(range(2),repeat=len(positions)):
        x=dict(zip(positions,bits))
        y={(r,c):value(x,r,c) for r in range(2) for c in range(-1,2)}
        assert value(y,0,0)==x[0,0]
        before=tuple(x[r,c] for r in range(2) for c in range(-1,2))
        after=tuple(y[r,c] for r in range(2) for c in range(-1,2))
        assert (after==u)==(before==v)
        assert (after==v)==(before==u)
        checked+=1
    return dict(shape=[2,3],u=list(u),v=list(v),
                positive_agreements=[list(d) for i,d in enumerate(ds) if masks[0]>>i&1],
                exact_local_contexts=checked)


def failed_naive_replacement():
    # Insert the 3x3 patch into the constant-one background.
    u={(0,0):0,(0,1):0,(1,0):0,(1,1):1}
    v={(0,0):1,(0,1):0,(1,0):0,(1,1):0}
    x={(r,c):int((r,c)==(0,0)) for r in range(3) for c in range(3)}
    def occurrences(patch,word):
        return {(a,b) for a in range(-2,4) for b in range(-2,4)
                if all(patch.get((a+r,b+c),1)==value for (r,c),value in word.items())}
    assert not occurrences(x,u) and occurrences(x,v)=={(0,0)}
    y=dict(x);y.update(u)
    assert occurrences(y,v)=={(1,1)}
    sites,ds,pairs=geometry((2,2))
    assert cross(pairs,tuple(u.values()),tuple(u.values()))==cross(pairs,tuple(v.values()),tuple(v.values()))
    return dict(u=[[0,0],[0,1]],v=[[1,0],[0,0]],
                input_patch=[[1,0,0],[0,0,0],[0,0,0]],exterior=1,
                old_v_anchors=[[0,0]],new_v_anchors=[[1,1]])


def audit():
    systems=[graph_audit(shape,q) for shape,q in
             [((2,2),2),((2,3),2),((2,4),2),((3,3),2),((2,2),3)]]
    return dict(status='Exact local and finite graph checks; all-dimensional conjugacy uses the working proof.',
                local_example=local_example(),naive_replacement_failure=failed_naive_replacement(),
                systems=systems)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--rebuild',action='store_true')
    args=parser.parse_args();result=audit()
    if args.rebuild:EVIDENCE.write_text(json.dumps(result,indent=2)+'\n')
    else:assert result==json.loads(EVIDENCE.read_text())
    print('SWAP CONJUGACY CHECKS VERIFIED',result['local_example']['exact_local_contexts'],flush=True)
    for r in result['systems']:
        print(r['shape'],r['q'],'classes',len(r['groups']),
              'disconnected',sum(len(g['components'])>1 for g in r['groups']),flush=True)


if __name__=='__main__':main()
