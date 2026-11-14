"""Exact audits of DI's theorem at off-target coin-count endpoints.

The proof is in drain-inventory.md. One audit compares the positive coefficient
expansion with exhaustive decorated-path enumeration; the other compares the
IDLA formula with independent raw-ARW absorption at sleep parameters 0 and 1.
No positive result here asserts the general mixed-coin inventory inequality.
Research code written by Codex, without independent user validation.
"""

from verify_response_renewal import stopped, inverse, multiply
from collections import defaultdict
from fractions import Fraction as F
import json
from verify_drain_inventory import path_family_count, inventory_poly, TooLarge
from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb
from random import Random

def corner_coefficients(common,extra,stock,n,T,target=0):
    edges=tuple(sorted(stock));initial=tuple(stock[e] for e in edges)
    T=frozenset(T);R=set(range(n))-T
    mv=sum(a for (z,w),a in stock.items() if z==target)
    zero=(0,)*(mv+1)
    def add(a,b):return tuple(x+y for x,y in zip(a,b))
    def choices(counts,z,allowed=None):
        for j,((site,d),a) in enumerate(zip(edges,counts)):
            if site!=z or not a or (allowed is not None and d not in allowed):continue
            rest=list(counts);rest[j]-=1
            yield d,tuple(rest)
    @lru_cache(None)
    def drain(counts):return path_family_count(dict(zip(edges,counts)),n)
    @lru_cache(None)
    def loops(counts,z,ell,hit,allowed,need):
        out=zero
        if z==target and (need<0 or hit):
            ways=drain(counts)
            if ways:
                b=mv-ell
                assert b>=1
                out=tuple(ways*comb(b-1,k-1) if 1<=k<=b else 0 for k in range(mv+1))
        for d,after in choices(counts,z,allowed):
            out=add(out,loops(after,d,ell+int(d==target),hit or d==need,allowed,need))
        return out
    def gain(counts,C,y):
        if y==-1:return zero
        D=tuple(sorted(C+(y,)))
        if target not in C and y!=target:return zero
        need=-1 if y==target else y
        allowed=frozenset(R|set(D))
        return loops(counts,target,0,False,allowed,need)
    @lru_cache(None)
    def idla(counts,C,label,z):
        settled=z==-1 or (z in T and z not in C)
        if settled:
            if label==len(common):return gain(counts,C,z)
            D=C if z==-1 else tuple(sorted(C+(z,)))
            label+=1
            root=common[label] if label<len(common) else extra
            return idla(counts,D,label,root)
        out=zero
        for d,after in choices(counts,z):out=add(out,idla(after,C,label,d))
        return out
    return idla(initial,(),0,common[0] if common else extra)

def coefficient_audit():
    cases=[([0],1,{(0,2):2,(2,0):2,(0,-1):3,(1,0):1},3),
           ([0],1,{(0,-1):2,(0,1):1,(1,-1):1,(1,2):1,(2,0):1},3)]
    rng=Random(42319)
    for trial in range(24):
        n=3;common=[0,2] if trial%2 else [1];extra=1 if trial%2 else 2
        starts=common+common+[extra];paths=[]
        for source in starts:
            path=[source]
            for _ in range(rng.randrange(0,3)):
                d=rng.choice([-1,0,1,2]);path.append(d)
                if d==-1:break
            if path[-1]!=-1:path.append(-1)
            paths.append(path)
        stock=Counter((a,b) for path in paths for a,b in zip(path,path[1:]))
        cases.append((common,extra,stock,n))
    checked=0
    for i,(common,extra,stock,n) in enumerate(cases):
        ways,poly=inventory_poly(common+common+[extra],stock,len(common),n)
        m=Counter()
        for (z,w),a in stock.items():m[z]+=a
        for bits in product((0,1),repeat=n-1):
            T={0}|{z+1 for z,b in enumerate(bits) if b}
            actual=corner_coefficients(common,extra,stock,n,T)
            expected=tuple(poly.get((k,)+tuple(m[z+1]*b for z,b in enumerate(bits)),0) for k in range(m[0]+1))
            assert actual==expected,(i,T,stock,actual,expected)
            checked+=1
    return dict(check='positive_corner_coefficients',exact_slices=checked,pools=len(cases),all_equal=True)

K=((F(0),F(1,4),F(1,4),F(0),F(1,2)),
   (F(0),F(0),F(1,2),F(1,2),F(0)),
   (F(0),F(0),F(1,2),F(1,2),F(0)),
   (F(1,2),F(0),F(0),F(1,2),F(0)))

def boundary_step(K,p,s):
    n=len(K);z=next(z for z,a in enumerate(s) if a>0)
    sleep=p[z] if s[z]==1 else F(0)
    out=defaultdict(F)
    if sleep:
        t=list(s);t[z]=-1;out[tuple(t)]+=sleep
    for w,q in enumerate(K[z]):
        if not q or sleep==1:continue
        t=list(s);t[z]-=1
        if w<n:t[w]=2 if t[w]==-1 else t[w]+1
        out[tuple(t)]+=(1-sleep)*q
    return out

def corner_payoff(K,T,v,p):
    n=len(K);T=frozenset(T);R=set(range(n))-T
    @lru_cache(None)
    def exit_law(C,root):
        available=sorted(T-set(C))+[n]
        if root in available:return {root:F(1)}
        inside=sorted(R|set(C))
        G=inverse([[F(z==w)-K[z][w] for w in inside] for z in inside])
        H=multiply(G,[[K[z][w] for w in available] for z in inside])
        return dict(zip(available,H[inside.index(root)]))
    @lru_cache(None)
    def payoff(C):
        if v not in C:return F(0)
        inside=sorted(R|(set(C)-{v}))
        H=multiply(inverse([[F(z==w)-K[z][w] for w in inside] for z in inside]),
                   [[K[z][v]] for z in inside])
        ret=K[v][v]+sum(K[v][z]*H[i][0] for i,z in enumerate(inside))
        return p/(1-(1-p)*ret)
    def expectation(eta):
        law={():F(1)}
        for z,count in enumerate(eta):
            for _ in range(count):
                nxt=defaultdict(F)
                for C,weight in law.items():
                    for y,q in exit_law(C,z).items():
                        D=C if y==n else tuple(sorted(C+(y,)))
                        nxt[D]+=weight*q
                law=nxt
        return sum(q*payoff(C) for C,q in law.items())
    return expectation,payoff

def probability_audit():
    initial=[bits for bits in product((0,1),repeat=4) if sum(bits)<=3]
    initial += [(2,1,0,0),(0,2,1,0)]
    checked=0
    for mask in product((0,1),repeat=3):
        p=(F(2,5),)+tuple(map(F,mask));T={0}|{z+1 for z,b in enumerate(mask) if b}
        predict,payoff=corner_payoff(K,T,0,p[0])
        laws,states=stopped(initial,lambda s:boundary_step(K,p,s),lambda s:all(a<=0 for a in s))
        for eta in initial:
            exact=sum(q for B,q in laws[eta].items() if B[0]==-1)
            assert predict(eta)==exact,(mask,eta,predict(eta),exact)
            checked+=1
        subsets=[tuple(z for z,b in zip(sorted(T),bits) if b) for bits in product((0,1),repeat=len(T))]
        assert all(payoff(C)<=payoff(D) for C in subsets for D in subsets if set(C)<=set(D))

    return dict(check='corner_idla_stabilization',exact_comparisons=checked,all_equal=True,all_payoffs_increasing=True)

def self_loop_elevation(a,b,r,loops):
    """Positive Bernstein coefficients after inserting a fixed self-loop count."""
    m=r+loops;out=[0]*(m+1)
    for i in range(loops+1):
        j=loops-i
        pa=comb(r+i-1,i)
        pb=comb(r+i-2,i) if r>1 else int(i==0)
        free=r+i-1
        for u in range(free+1):
            out[u+1]+=a*pa*comb(free,u)
            out[u]+=b*pb*comb(free,u)
    return tuple(out)

def self_loop_audit():
    rng=Random(90141);checked=0
    for trial in range(32):
        n=3;common=[0,1] if trial%2 else [0];extra=2;starts=common+common+[extra]
        # Vertex 1 can have initial/self visits and exits to 0 or 2, but no external incoming edge.
        paths=[]
        for source in starts:
            path=[source]
            for _ in range(rng.randrange(0,3)):
                d=rng.choice([-1,0,2]);path.append(d)
                if d==-1:break
            if path[-1]!=-1:path.append(-1)
            paths.append(path)
        stock=Counter((a,b) for path in paths for a,b in zip(path,path[1:]))
        if not sum(a for (z,w),a in stock.items() if z==1):continue
        ways,base=inventory_poly(starts,stock,len(common),n)
        r=sum(a for (z,w),a in stock.items() if z==1);m0=sum(a for (z,w),a in stock.items() if z==0);m2=sum(a for (z,w),a in stock.items() if z==2)
        for loops in range(3):
            full=Counter(stock);full[1,1]+=loops
            _,poly=inventory_poly(starts,full,len(common),n)
            for k0 in range(m0+1):
                for bit in (0,1):
                    a=base.get((k0,r,m2*bit),0);b=base.get((k0,0,m2*bit),0)
                    expect=self_loop_elevation(a,b,r,loops)
                    actual=tuple(poly.get((k0,k,m2*bit),0) for k in range(r+loops+1))
                    assert actual==expect,(trial,loops,stock,k0,bit,actual,expect)
                    checked+=1
    source_checks=checked
    rng=Random(426171);checked=passed=excluded=0
    allowed={0:[0,-1],1:[0,-1],2:[0,1,3,-1],3:[0,1,2,-1]}
    for trial in range(16):
        starts=[0,2,0,2,3];stock=Counter()
        for source in starts:
            z=source
            for _ in range(rng.randrange(1,3)):
                d=rng.choice(allowed[z]);stock[z,d]+=1;z=d
                if z==-1:break
            if z!=-1:stock[z,-1]+=1
        m=Counter()
        for (z,w),a in stock.items():m[z]+=a
        r=m[1]
        if not r:continue
        try:_,base=inventory_poly(starts,stock,2,4)
        except TooLarge:excluded+=1;continue
        for loops in (1,2):
            full=Counter(stock);full[1,1]+=loops
            try:_,poly=inventory_poly(starts,full,2,4)
            except TooLarge:excluded+=1;continue
            for k0 in range(m[0]+1):
                for b2,b3 in product((0,1),repeat=2):
                    a=base.get((k0,r,b2*m[2],b3*m[3]),0)
                    b=base.get((k0,0,b2*m[2],b3*m[3]),0)
                    expect=self_loop_elevation(a,b,r,loops)
                    actual=tuple(poly.get((k0,k,b2*m[2],b3*m[3]),0) for k in range(r+loops+1))
                    assert expect==actual,(trial,stock,loops,k0,b2,b3,expect,actual)
                    checked+=1
            passed+=1
    return dict(check='mixed_sites_without_external_returns',source_only_checks=source_checks,
                incoming_checks=checked,incoming_pools=passed,excluded=excluded,all_equal=True)

if __name__=='__main__':
    print(json.dumps(coefficient_audit()),flush=True)
    print(json.dumps(probability_audit()),flush=True)
    print(json.dumps(self_loop_audit()),flush=True)
