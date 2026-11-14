"""Exact audit of the new avoidance-defect basis, using Fraction arithmetic.

Rebuilds the signed renewal matrix from raw ARW transitions, eliminating
sleepers, performing IDLA, and removing occupied sources. No saved numerical
matrix is used. The chosen four-site directed example has p=1/2 everywhere.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json


def eye(n):
    return [[F(i==j) for j in range(n)] for i in range(n)]


def inverse(M):
    n=len(M);a=[list(r)+b for r,b in zip(M,eye(n))]
    for i in range(n):
        j=next(j for j in range(i,n) if a[j][i])
        a[i],a[j]=a[j],a[i];q=a[i][i];a[i]=[x/q for x in a[i]]
        for j in range(n):
            if j!=i and a[j][i]:
                q=a[j][i];a[j]=[x-q*y for x,y in zip(a[j],a[i])]
    return [r[n:] for r in a]


def multiply(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def solve_rows(rows, rhs):
    """Sparse rational elimination, followed by back substitution."""
    rows=[dict(r) for r in rows];rhs=[list(b) for b in rhs]
    count=len(rows);width=len(rhs[0]) if count else 0
    for i in range(count):
        pivot=rows[i][i]
        assert pivot>0
        tail={k:q/pivot for k,q in rows[i].items() if k>i and q}
        rows[i]=tail;rhs[i]=[q/pivot for q in rhs[i]]
        for j in range(i+1,count):
            q=rows[j].pop(i,F(0))
            if not q:continue
            for k,a in tail.items():
                rows[j][k]=rows[j].get(k,F(0))-q*a
                if not rows[j][k]:rows[j].pop(k)
            rhs[j]=[a-q*b for a,b in zip(rhs[j],rhs[i])]
    values=[[F(0)]*width for _ in rows]
    for i in range(count-1,-1,-1):
        values[i]=[rhs[i][v]-sum(q*values[j][v] for j,q in rows[i].items())
                   for v in range(width)]
    return values


def stopped(initial, step, stop):
    initial=list(dict.fromkeys(initial))
    states=[s for s in initial if not stop(s)];idx={s:i for i,s in enumerate(states)}
    boundary={s:i for i,s in enumerate(s for s in initial if stop(s))}
    rows=[];exits=[]
    for s in states:
        row=defaultdict(F);row[idx[s]]=F(1);out=defaultdict(F)
        for t,q in step(s).items():
            if not q:continue
            if stop(t):
                if t not in boundary:boundary[t]=len(boundary)
                out[boundary[t]]+=q
            else:
                if t not in idx:idx[t]=len(states);states.append(t)
                row[idx[t]]-=q
        rows.append(row);exits.append(out)
    rhs=[[b.get(j,F(0)) for j in range(len(boundary))] for b in exits]
    sol=solve_rows(rows,rhs)
    ans={}
    for s in initial:
        if stop(s):ans[s]={s:F(1)}
        else:ans[s]={b:sol[idx[s]][j] for b,j in boundary.items() if sol[idx[s]][j]}
    return ans,len(states)


def arw_step(K,p,s,frozen=None):
    n=len(K);z=next(z for z,a in enumerate(s) if a>0 and z!=frozen)
    out=defaultdict(F);t=list(s)
    if t[z]==1:t[z]=-1
    out[tuple(t)]+=p[z]
    for w,q in enumerate(K[z]):
        if not q:continue
        t=list(s);t[z]-=1
        if w<n:t[w]=2 if t[w]==-1 else t[w]+1
        out[tuple(t)]+=(1-p[z])*q
    return out


def exact_renewal(K,p,N):
    n=len(K)
    allactive=[a for a in product(range(N+1),repeat=n) if 0<sum(a)<=N]
    binary=[a for a in allactive if max(a)<=1]
    edges=[(a,x) for a in binary for x in range(n) if not a[x]]
    ei={e:i for i,e in enumerate(edges)}
    def idla_step(a):
        z=next(z for z,k in enumerate(a) if k>1);out=defaultdict(F)
        for w,q in enumerate(K[z]):
            b=list(a);b[z]-=1
            if w<n:b[w]+=1
            out[tuple(b)]+=q
        return out
    idla,_=stopped(allactive,idla_step,lambda a:max(a)<=1)
    @lru_cache(None)
    def occupied_source(a,x):
        if not a[x]:return {(a,x):F(1)}
        A=[z for z,k in enumerate(a) if k];B=[z for z,k in enumerate(a) if not k]
        M=[[F(z==w)-K[z][w] for w in A] for z in A]
        R=[[K[z][w] for w in B] for z in A]
        h=multiply(inverse(M),R)
        return {(a,w):h[A.index(x)][j] for j,w in enumerate(B) if h[A.index(x)][j]}
    @lru_cache(None)
    def allactive_response(a,x):
        out=defaultdict(F)
        for b,q in idla[a].items():
            for edge,t in occupied_source(b,x).items():out[edge]+=q*t
        return dict(out)
    @lru_cache(None)
    def expand(c,x):
        if min(c)>=0:return allactive_response(c,x)
        z=next(z for z,k in enumerate(c) if k==-1)
        base=list(c);base[z]=0
        active=base[:];active[z]=1
        out=defaultdict(F)
        for e,q in expand(tuple(active),x).items():out[e]+=q/p[z]
        for w,q in enumerate(K[z]):
            if not q:continue
            b=base[:]
            if w<n:b[w]=2 if b[w]==-1 else b[w]+1
            for e,t in expand(tuple(b),x).items():out[e]-=(1-p[z])*q*t/p[z]
        return dict(out)
    traps={};counts=[]
    for y in range(n):
        starts=[a for a in binary if not a[y]]
        law,count=stopped(starts,lambda c:arw_step(K,p,c,y),
                         lambda c:all(k<=0 for z,k in enumerate(c) if z!=y))
        traps[y]=law;counts.append(count)
    M=[[F(0)]*len(edges) for _ in edges];H={}
    for a in binary:
        k=[p[z] if not a[z] else F(0) for z in range(n)]
        IQ=[[F(z==w)-(1-k[z])*K[z][w] for w in range(n)] for z in range(n)]
        Green=inverse(IQ)
        H[a]=[[Green[z][w]*k[w] for w in range(n)] for z in range(n)]
    for i,(a,x) in enumerate(edges):
        for y in range(n):
            if a[y] or not H[a][x][y]:continue
            for c,q in traps[y][a].items():
                if c[y]<=0:continue
                for e,t in expand(c,y).items():M[i][ei[e]]+=H[a][x][y]*q*t
    return edges,M,H,counts


def main():
    K=[[F(0),F(0),F(0),F(1,2),F(1,2)],
       [F(1,4),F(0),F(1,4),F(1,4),F(1,4)],
       [F(0),F(0),F(0),F(1,2),F(1,2)],
       [F(0),F(0),F(1,2),F(0),F(1,2)]]
    p=[F(1,2)]*4;v=2
    edges,M,H,counts=exact_renewal(K,p,2)
    subsets=[a for a in product([0,1],repeat=4) if sum(a)<=2 and not a[v]]
    si={a:i for i,a in enumerate(subsets)}
    G=[[F(0)]*len(subsets) for _ in edges]
    for i,(a,x) in enumerate(edges):
        if a in si:G[i][si[a]]=H[a][x][v]
    IM=[[F(i==j)-M[i][j] for j in range(len(edges))] for i in range(len(edges))]
    L=multiply(inverse(IM),G)
    row=L[edges.index(((1,1,0,0),3))]
    expected={(0,0,0,1):-F(1,1560),(0,1,0,1):F(1,26),
              (1,0,0,1):F(7,240),(1,1,0,0):F(2,15)}
    assert {s:q for s,q in zip(subsets,row) if q}==expected
    # For common set {0,3}, force the particle at 0 to jump. It exits
    # immediately with probability 1/2. If it reaches 3, move the surplus:
    # with probability 1/2 it exits, or otherwise hits the target 2.
    # Therefore J_0 a({0,3}) = (3/4) a({3}).
    defect=[F(0)]*len(subsets)
    defect[si[(1,0,0,1)]]=F(1)
    defect[si[(0,0,0,1)]]=-F(3,4)
    alpha=F(1,1170)
    remainder=[q-alpha*d for q,d in zip(row,defect)]
    assert all(q>=0 for q in remainder)
    assert remainder[si[(1,0,0,1)]]==F(53,1872)
    # Independent direct raw-instruction absorption checks the real response
    # and the actual avoidance values used in the identity.
    lo=(1,1,0,0);hi=(1,1,0,1)
    laws,rawcount=stopped([lo,hi],lambda s:arw_step(K,p,s),lambda s:all(c<=0 for c in s))
    values={s:sum(q*(t[v]==-1) for t,q in law.items()) for s,law in laws.items()}
    avoid,avoidcount=stopped(subsets,lambda s:arw_step(K,p,s),
                            lambda s:s[v]>0 or all(c<=0 for c in s))
    a={s:sum(q for t,q in law.items() if t[v]<=0) for s,law in avoid.items()}
    actual=values[hi]-values[lo]
    assert sum(q*a[s] for s,q in zip(subsets,row))==actual
    assert sum(q*a[s] for s,q in zip(subsets,remainder))+alpha*sum(q*a[s] for s,q in zip(subsets,defect))==actual
    out=dict(K=[[str(q) for q in r] for r in K],p=[str(q) for q in p],
             common=lo,extra=3,target=v,edges=len(edges),trap_state_counts=counts,
             exact_coefficient_identity=[dict(A=s,coefficient=str(q)) for s,q in zip(subsets,row) if q],
             defect=dict(A=[0,3],forced_source=0,value='a({0,3}) - 3/4 a({3})',weight=str(alpha)),
             nonnegative_remainder=[dict(A=s,coefficient=str(q)) for s,q in zip(subsets,remainder) if q],
             actual_avoidance=[dict(A=s,value=str(a[s])) for s in subsets],
             lower_occupation=str(values[lo]),upper_occupation=str(values[hi]),
             actual_increment=str(actual),raw_state_count=rawcount,avoidance_state_count=avoidcount,
             verified=True)
    with open('avoidance_exact_verification.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
