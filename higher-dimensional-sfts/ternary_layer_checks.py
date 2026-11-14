"""Independent replay of shared-letter ternary cover certificates.

Research code developed by ChatGPT from Jacob Richey's specifications.
Positive witnesses are literal patterns and line directions. Negative cases
exhaust forced-component colorings in their natural order, independently of
the heuristic search order. These are proof-route obstructions, not entropy
counterexamples. The analytic shared-letter theorem is in entropy-ordering.md.
"""
import argparse
from collections import Counter
from functools import lru_cache
import json
from math import gcd
from pathlib import Path

from ordering_checks import agreement, rectangle, displacements
import alphabet_sweep_checks as alphabet

ROOT = Path(__file__).parent
EVIDENCE = ROOT/'ternary-layer-evidence.json'


def all_words(sites, g):
    """Natural-order restricted-growth colorings with literal final checks."""
    n = len(sites); parent = list(range(n)); index = {z:i for i,z in enumerate(sites)}
    ds = displacements(sites)
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for k, (dx, dy) in enumerate(ds):
        if not g >> k & 1:
            continue
        for i, (x, y) in enumerate(sites):
            if (x-dx,y-dy) in index:
                parent[root(i)] = root(index[x-dx,y-dy])
    ids = {}; labels = []
    for i in range(n):
        r = root(i)
        if r not in ids: ids[r] = len(ids)
        labels.append(ids[r])
    c = len(ids); at = [[] for _ in range(c)]
    for k, (dx, dy) in enumerate(ds):
        if g >> k & 1: continue
        edges = {(labels[i], labels[index[x-dx,y-dy]]) for i, (x,y) in enumerate(sites)
                 if (x-dx,y-dy) in index and labels[i] != labels[index[x-dx,y-dy]]}
        assert edges
        at[max(max(a,b) for a,b in edges)].append(edges)
    values = [0]*c
    def visit(i, top):
        if i == c:
            word = tuple(values[k] for k in labels)
            assert agreement(sites, word) == g
            yield word; return
        for color in range(min(top+2,3)):
            values[i] = color
            if all(any(values[a] != values[b] for a,b in edges) for edges in at[i]):
                yield from visit(i+1,max(top,color))
    yield from visit(0,-1)


def audit(data):
    source = json.loads(alphabet.EVIDENCE.read_text())
    expected = {tuple(r['shape']): {(g,h) for g,h,*_ in r['generators']}
                for r in source['critical_lines'] if r['q'] == 3}
    statistics = []; all_negative = []
    for record in data['systems']:
        dims = tuple(record['shape']); sites = rectangle(*dims); ds = displacements(sites)
        @lru_cache(None)
        def profile(text):
            word = tuple(map(int,text)); values = dict(zip(sites,word))
            assert set(word) <= {0,1,2}
            mismatches = [tuple((x,y) for x,y in sites if (x-dx,y-dy) in values
                          and values[x,y] != values[x-dx,y-dy]) for dx,dy in ds]
            return agreement(sites,word), mismatches
        positives = {(row[0],row[1]) for row in record['certificates']}
        negatives = {(row[0],row[1]) for row in record['failures']}
        assert not positives & negatives
        assert len(positives) == len(record['certificates'])
        assert positives | negatives == expected[dims]
        normals = Counter()
        for g,h,u,v,a,b,cu,cv in record['certificates']:
            assert g != h and g & h == g and gcd(a,b) == 1
            assert tuple(c == str(cu) for c in u) == tuple(c == str(cv) for c in v)
            gu, before = profile(u); gv, after = profile(v)
            assert (gu,gv) == (g,h)
            projection = lambda z:a*z[0]+b*z[1]
            assert all(not y or x and min(map(projection,y)) <= max(map(projection,x))
                       for x,y in zip(before,after))
            normals['axis' if a == 0 or b == 0 else 'oblique'] += 1
        for negative in record['negative_shared_search']:
            g,h = negative['g'],negative['h']; counts=[]; masks=[]
            for k in (g,h):
                count=0; classes=set()
                for word in all_words(sites,k):
                    count += 1
                    classes.update(sum(1<<i for i,c in enumerate(word) if c == letter)
                                   for letter in set(word))
                counts.append(count);masks.append(classes)
            assert counts == [negative['source_count'],negative['target_count']]
            assert not masks[0] & masks[1]
            assert negative['shared_pairs'] == 0
            all_negative.append(dict(shape=list(dims),g=g,h=h,representative_counts=counts))
        statistics.append(dict(shape=list(dims),shared_covers=len(positives),
                               no_shared_letter_covers=len(negatives),**normals))
        print('TERNARY LAYERS REPLAYED',dims,len(positives),len(negatives),dict(normals),flush=True)
    assert sum(r['shared_covers'] for r in statistics) == 30756
    assert len(all_negative) == 32
    return dict(systems=statistics,negative_certificates=all_negative)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record-audit',action='store_true')
    args=parser.parse_args(); data=json.loads(EVIDENCE.read_text()); result=audit(data)
    if args.record_audit:
        data['status']='Literal shared-letter line certificates and independent exhaustive negative audits; entropy implications use the analytic proofs.'
        data['independent_audit']=result
        EVIDENCE.write_text(json.dumps(data,indent=2)+'\n')
    else:
        assert result == data['independent_audit']
    print('TERNARY LAYER CERTIFICATES VERIFIED',flush=True)


if __name__ == '__main__':
    main()
