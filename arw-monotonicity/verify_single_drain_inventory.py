"""Exact checks of the stronger one-family inventory candidate in DI.

Lower and upper share common complete decorated paths; the upper also uses
one extra path. The aggregate inventory can vary their common geometry.
General coefficient positivity remains open. Standard library only;
research code written by Codex, without independent user validation.
"""

from verify_drain_inventory import (occupation_poly,add_counts,binomial_poly,
                                   path_family_count,TooLarge)
from verify_inventory_corners import corner_coefficients, stopped, boundary_step
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import product
from random import Random
import gc,json
from math import comb,factorial
from fractions import Fraction as F

def single_inventory(common,extra,stock,n,limit=10000,target=0):
    starts=common+[extra];left=dict(stock);count=0;poly=defaultdict(int)
    def rec(label,path,finished):
        nonlocal count
        z=path[-1]
        if z==-1:
            paths=finished+(tuple(path),)
            if label+1==len(starts):
                if any(left.values()):return
                count+=1
                if count>limit:raise TooLarge
                plus=occupation_poly(paths,n,target)
                nc=Counter(paths[-1][:-1])
                minus=add_counts(occupation_poly(paths[:-1],n,target),binomial_poly(tuple(nc[z] for z in range(n))))
                for key,value in plus.items():poly[key]+=value
                for key,value in minus.items():poly[key]-=value
                return
            rec(label+1,[starts[label+1]],paths);return
        for (site,d),a in list(left.items()):
            if site!=z or not a:continue
            left[site,d]-=1;rec(label,path+[d],finished);left[site,d]+=1
    rec(0,[starts[0]],())
    return count,{key:value for key,value in poly.items() if value}


def single_dp(common,extra,stock,n,target=0,guard=150000,coefficient_guard=750000):
    edges=tuple(sorted(stock));initial=tuple(stock[e] for e in edges);visited=0;cached_coefficients=0
    @lru_cache(None)
    def closure(counts):
        ways=path_family_count(dict(zip(edges,counts)),n)
        if not ways:return {}
        m=[0]*n
        for (z,w),a in zip(edges,counts):m[z]+=a
        return {key:ways*value for key,value in binomial_poly(tuple(m)).items()}
    @lru_cache(None)
    def rec(counts,queue,residents):
        nonlocal visited,cached_coefficients
        visited+=1
        if visited>guard:raise TooLarge
        if not queue:return closure(counts) if residents[target]!=-2 else {}
        z=queue[0];pending=list(queue[1:]);r=list(residents);old=r[z];r[z]=-2
        if old>=0:pending.append(old)
        out=defaultdict(int)
        for j,((site,d),a) in enumerate(zip(edges,counts)):
            if site!=z or not a:continue
            after=list(counts);after[j]-=1
            for bit in (0,1):
                rr=list(r);qq=list(pending)
                if bit:rr[z]=d
                elif d>=0:qq.append(d)
                for key,value in rec(tuple(after),tuple(qq),tuple(rr)).items():
                    kk=list(key);kk[z]+=bit;out[tuple(kk)]+=value
        cached_coefficients+=len(out)
        if cached_coefficients>coefficient_guard:raise TooLarge
        return {key:value for key,value in out.items() if value}
    lower=rec(initial,tuple(common),(-2,)*n)
    upper=rec(initial,tuple(common+[extra]),(-2,)*n)
    poly=defaultdict(int,upper)
    for key,value in lower.items():poly[key]-=value
    return path_family_count(stock,n),{key:value for key,value in poly.items() if value},visited


base_terms=[((1,1,1,1),(1,1,1,0),1,1),
            ((2,1,0,2),(1,1,1,0),1,1),
            ((1,2,2,0),(1,1,2,0),2,1),
            ((2,2,1,1),(1,1,2,0),2,2),
            ((3,2,0,2),(1,1,2,0),2,1)]
extra_terms=[((0,0,1,1),(1,0,0,0),1,1),
             ((1,0,0,2),(1,0,0,0),1,1),
             ((0,1,2,0),(1,0,1,0),2,1),
             ((1,1,1,1),(1,0,1,0),2,2),
             ((2,1,0,2),(1,0,1,0),2,1)]

def positive_formula(J,common):
    inv=(J.get((0,1),0),J.get((1,0),0),J.get((0,-1),0),J.get((1,-1),0))
    m0=inv[0]+inv[2];m1=inv[1]+inv[3];out=defaultdict(int)
    terms=[(term,True) for term in base_terms]
    if common==1:terms += [(term,False) for term in extra_terms]
    for (e,coins,h,mult),second_den in terms:
        left=tuple(a-b for a,b in zip(inv,e))
        if left[2:]!=(0,0) or left[0]!=left[1] or left[0]<0:continue
        total=left[0]
        for i in range(total+1):
            for j in range(total-i+1):
                k=total-i-j
                if k and not second_den:continue
                for a in range(k+1):
                    for b in range(k-a+1):
                        c=k-a-b
                        coeff=mult*comb(h+i-1,i)*factorial(k)//(factorial(a)*factorial(b)*factorial(c))
                        o0,z0,o1,z1=coins
                        o0+=b;z0+=j+a+c;o1+=c;z1+=j+a+b
                        pad0=m0-o0-z0;pad1=m1-o1-z1
                        assert min(pad0,pad1)>=0,(inv,e,i,j,k,coins,pad0,pad1)
                        for b0 in range(pad0+1):
                            for b1 in range(pad1+1):
                                out[o0+b0,o1+b1]+=coeff*comb(pad0,b0)*comb(pad1,b1)
    return dict(out)


def insert_self_loops(poly,m,site,loops):
    out=defaultdict(int);r=m[site]
    for key,value in poly.items():
        k=key[site]
        for i in range(loops+1):
            j=loops-i
            left=comb(k+i-1,i) if k else int(i==0)
            right=comb(r+j-1,j) if r else int(j==0)
            for u in range(i+1):
                kk=list(key);kk[site]+=u
                out[tuple(kk)]+=value*left*right*comb(i,u)
    return {key:value for key,value in out.items() if value}


def two_site_audit():
    count=loop_checks=0
    for r in range(6):
        for inv in ((r,r,1,1),(r,r+1,2,0),(r+1,r,0,2)):
            J=dict(zip(((0,1),(1,0),(0,-1),(1,-1)),inv))
            for common in (0,1):
                _,actual=single_inventory([common],1-common,J,2)
                predicted=positive_formula(J,common)
                assert actual==predicted,(r,inv,common,actual,predicted)
                count+=1
    for r in (1,2):
        J={(0,1):r,(1,0):r,(0,-1):1,(1,-1):1}
        for common in (0,1):
            base=positive_formula(J,common)
            for l0,l1 in ((1,0),(0,1),(1,1),(2,1)):
                full=dict(J);full[0,0]=l0;full[1,1]=l1
                predicted=insert_self_loops(base,(r+1,r+1),0,l0)
                predicted=insert_self_loops(predicted,(r+1+l0,r+1),1,l1)
                _,actual=single_inventory([common],1-common,full,2)
                assert actual==predicted,(r,common,l0,l1,actual,predicted)
                loop_checks+=1
    return dict(check='two_site_positive_generating_formula',exact_pools=count,
                self_loop_insertions=loop_checks,all_equal=True)


three_cycle_edges=Counter({(0,1):1,(1,2):1,(2,0):1})
three_cycle_terms={0:[({(1,-1):1},(1,1,0,1),False),
          ({(1,-1):1},(1,1,1,0),True),
          ({(1,2):1,(2,-1):1},(1,1,1,0),True),
          ({(1,2):1,(2,0):1},(1,1,2,0),True)],
       1:[({(0,-1):1},(0,1,1,1),False),
          ({(0,1):1,(1,-1):1},(0,1,1,1),False),
          ({(0,1):1,(1,2):1},(1,1,1,1),False),
          ({(0,-1):1},(1,0,1,1),True),
          ({(0,1):1},(2,0,1,1),True)]}

def positive_cycle_formula(J,target):
    m=Counter()
    for (z,w),a in J.items():m[z]+=a
    out=defaultdict(int)
    for prefix,coins,second_den in three_cycle_terms[target]:
        used=three_cycle_edges+Counter(prefix)
        rem={edge:J.get(edge,0)-used[edge] for edge in set(J)|set(used)}
        if min(rem.values())<0:continue
        bound=min(rem[e] for e in three_cycle_edges)
        for j in range(bound+1):
            for k in range(bound-j+1 if second_den else 1):
                rest={e:a-(j+k)*three_cycle_edges[e] for e,a in rem.items()}
                ways=path_family_count(rest,3)
                if not ways:continue
                for a in range(k+1):
                    for b in range(k-a+1):
                        c=k-a-b;weight=ways*factorial(k)//(factorial(a)*factorial(b)*factorial(c))
                        o0,z0,o1,z1=coins
                        o0+=b;z0+=j+a+c;o1+=c;z1+=j+a+b
                        f0=m[0]-o0-z0;f1=m[1]-o1-z1
                        assert min(f0,f1)>=0
                        for u in range(f0+1):
                            for v in range(f1+1):
                                out[o0+u,o1+v,m[2]]+=weight*comb(f0,u)*comb(f1,v)
    return dict(out)

def three_cycle_audit():
    rng=Random(8928)
    for trial in range(16):
        u,b,g,p0,p1=[rng.choice([F(1,4),F(1,2),F(3,4)]) for _ in range(5)]
        K=((F(0),u,F(0),1-u),(F(0),F(0),b,1-b),(g,F(0),F(0),1-g))
        p=(p0,p1,F(1));q0=1-p0;q1=1-p1;t=u*b*g;D=1-q0*q1*t;D2=1-(1-p0*p1)*t
        low0=p0+q0*u*q1*b*b*g*p0/D
        low1=p1+p1*q0*q0*q1*u*u*b*b*g/D
        up0=p0/D+p0*p1*q0*t*(1-q1*b*g)/(D*D2)
        up1=p1/D+p0*p1*q1*t*(1-q0*u)/(D*D2)
        laws,_=stopped(((1,1,0),(1,1,1)),lambda s:boundary_step(K,p,s),lambda s:all(a<=0 for a in s))
        for eta,expected in (((1,1,0),(low0,low1)),((1,1,1),(up0,up1))):
            actual=tuple(sum(weight for state,weight in laws[eta].items() if state[v]==-1) for v in (0,1))
            assert actual==expected,(trial,eta,actual,expected)

    checked=0
    for a in product(range(4),repeat=3):
        sink=tuple(1+a[z-1]-a[z] for z in range(3))
        if min(sink)<0:continue
        J=Counter({(z,(z+1)%3):a[z] for z in range(3)})
        J.update({(z,-1):sink[z] for z in range(3)})
        if not path_family_count(J,3):continue
        m2=sum(c for (z,w),c in J.items() if z==2)
        for target in (0,1):
            _,all_coins=single_inventory([0,1],2,J,3,target=target)
            actual={key:value for key,value in all_coins.items() if key[2]==m2}
            predicted=positive_cycle_formula(J,target)
            assert predicted==actual,(a,sink,target,predicted,actual)
            checked+=1
    return dict(check='permanent_extra_on_directed_three_cycle',raw_parameter_cases=16,
                exact_inventory_target_comparisons=checked,all_equal=True)


if __name__=='__main__':
    J={(0,2):2,(2,0):2,(0,-1):2,(1,0):1}
    ways,c=single_inventory([0],1,J,3)
    assert ways==3 and c=={(1,0,2):1,(2,0,2):2,(3,0,2):1}
    rng=Random(9162);corner_slices=0
    for trial in range(32):
        n=4;common=[0,2] if trial%2 else [1,2];extra=1 if trial%2 else 3;J=Counter()
        for source in common+[extra]:
            z=source
            for _ in range(rng.randrange(0,3)):
                d=rng.choice([-1,0,1,2,3]);J[z,d]+=1;z=d
                if z==-1:break
            if z!=-1:J[z,-1]+=1
        expected=single_inventory(common,extra,J,n)
        actual=single_dp(common,extra,J,n)
        assert expected==actual[:2],(trial,expected,actual)
        assert all(value>=0 for value in actual[1].values()),(trial,J,actual)
        m=Counter()
        for (z,w),value in J.items():m[z]+=value
        for bits in product((0,1),repeat=3):
            T={0}|{z+1 for z,bit in enumerate(bits) if bit}
            predicted=corner_coefficients(common,extra,J,n,T)
            values=tuple(actual[1].get((k,)+tuple(m[z+1]*bit for z,bit in enumerate(bits)),0)
                         for k in range(m[0]+1))
            assert predicted==values,(trial,T,J,predicted,values)
            corner_slices+=1
    print(json.dumps(dict(check='single_family_enumeration',exact_comparisons=32,
                          endpoint_formula_slices=corner_slices,all_equal=True,
                          three_bundle_identity=True,all_coefficients_nonnegative=True)),flush=True)
    checked=[];excluded=[]
    for kind in (0,1,2,6):
        n=4;common=[0,2] if kind<6 else [0,1,2];extra=1 if kind<6 else 3
        J=Counter((z,-1) for z in common+[extra])
        for z in range(n):J[z,(z+1)%n]+=1
        for repeat in range(1+kind%4):
            a=repeat%n;b=(a+2)%n;J[a,b]+=1;J[b,a]+=1
        if kind%2:
            for z in range(n):J[z,(z+1)%n]+=1
        try:ways,c,states=single_dp(common,extra,J,n)
        except TooLarge:
            excluded.append(kind);gc.collect();continue
        assert all(value>=0 for value in c.values()),(kind,J,c)
        checked.append(dict(kind=kind,edges=sum(J.values()),families=ways,states=states))
        gc.collect()
    print(json.dumps(dict(check='single_family_interacting_cycles',completed=checked,
                          excluded_at_computational_guards=excluded,
                          status='Finite evidence; general one-family inequality open')),flush=True)
    print(json.dumps(two_site_audit()),flush=True)
    print(json.dumps(three_cycle_audit()),flush=True)
