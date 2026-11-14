"""Exact excursion-block and graphical-dual audits; standard library only.

AI-written research code. General proper-target coefficient positivity fails.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import verify_response_renewal as rn
from verify_creator_potentials import determinant, subsets


def full_target_excursion_dual_audit():
    """Physical pre-loss activity process and its positive subset dual."""
    from collections import defaultdict
    from math import prod
    entries=balance=values=0
    for case in rn.CASES:
        K=[[F(q) for q in row] for row in case["K"]]
        p=list(map(F,case["p"]))
        n=len(p);V=(1<<n)-1;masks=range(1<<n)
        excursions={}
        for i in range(n):
            start=(i,0)
            def walk_step(state):
                z,mask=state
                return {(j,mask | (1<<j) if j<n else mask):q
                        for j,q in enumerate(K[z]) if q}
            laws,_=rn.stopped([start],walk_step,
                              lambda state:state[0]==n or state[0]==i and state[1]!=0)
            excursions[i]=laws[start]
        dual={};phase={}
        fail=("lost",)
        for A in masks:
            row=defaultdict(F)
            for i in range(n):
                if not A & (1<<i):
                    row[A]+=F(1,n)
                    continue
                row[A & ~(1<<i)]+=p[i]/n
                for (end,visited),mass in excursions[i].items():
                    row[fail if end==n else A | visited]+=(1-p[i])*mass/n
            phase[A]=dict(row)
        determinantC={}
        for C in masks:
            sites=[i for i in range(n) if C & (1<<i)]
            determinantC[C]=determinant([[F(i==j)-K[i][j] for j in sites] for i in sites])
        for H in masks:
            row=defaultdict(F)
            for i in range(n):
                row[H & ~(1<<i)]+=p[i]/n
                if H & (1<<i):
                    row[H]+=(1-p[i])/n
                    continue
                e=sum(q for (end,visited),q in excursions[i].items()
                      if end==n or visited & H)
                assert e == determinantC[V ^ H]/determinantC[V ^ (H | (1<<i))]
                row[H | (1<<i)]+=(1-p[i])*e/n
                row[H]+=(1-p[i])*(1-e)/n
            dual[H]=dict(row)
        weight={H:prod((1-p[i])/p[i] for i in range(n) if H & (1<<i))
                    /determinantC[V ^ H] for H in masks}
        total=sum(weight.values())
        stationary={H:q/total for H,q in weight.items()}
        for H in masks:
            for J in masks:
                assert stationary[H]*dual[H].get(J,0) == stationary[J]*dual[J].get(H,0)
                balance+=1
        for A in masks:
            for H in masks:
                lhs=sum(q*F(not B & H) for B,q in phase[A].items() if B!=fail)
                rhs=sum(q*F(not A & J) for J,q in dual[H].items())
                assert lhs == rhs,(A,H,lhs,rhs)
                entries+=1
        killed,_=rn.stopped(masks,lambda A:phase[A],lambda A:A==0 or A==fail)
        starts={S:tuple(-1 if S & (1<<i) else 1 for i in range(n)) for S in masks}
        raw,_=rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                         lambda state:all(q<=0 for q in state))
        for S,start in starts.items():
            value=sum(q for H,q in stationary.items() if H & S == H)
            assert value == killed[V ^ S].get(0,F(0))
            assert value == raw[start].get((-1,)*n,F(0))
            values+=1
    # With a vacancy, the analogous activation-field sign can fail.
    K=[[F(0),F(0),F(1),F(0)],[F(0),F(0),F(1),F(0)],
       [F(1,2),F(0),F(0),F(1,2)]]
    p=[F(1,2)]*3
    starts=[(1,1,0),(1,-1,0)]
    raw,_=rn.stopped(starts,lambda state:rn.arw_step(K,p,state),
                     lambda state:all(q<=0 for q in state))
    values_partial=[sum(q for end,q in raw[start].items() if end[0]==-1) for start in starts]
    assert values_partial == [F(43,70),F(4,7)]
    return dict(excursion_duality_entries=entries, detailed_balance_entries=balance,
                full_target_values=values, partial_base_activation_coefficient="-3/70")



def proper_target_field_obstruction():
    """A negative full-input coefficient on an irreducible three-site graph."""
    K=[[F(0),F(0),F(1,2),F(1,2)],
       [F(1,2),F(0),F(0),F(1,2)],
       [F(15,16),F(3,64),F(0),F(1,64)]]
    p=[F(1,1001)]*3
    starts={(): (1,1,1), (1,): (1,-1,1),
            (2,): (1,1,-1), (1,2): (1,-1,-1)}
    raw,_=rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                     lambda state:all(q<=0 for q in state))
    moments={C:sum(q for end,q in raw[state].items() if end[0]==-1)
             for C,state in starts.items()}
    coefficient=moments[1,2]-moments[1,]-moments[2,]+moments[()]
    assert coefficient==F(-4885292907192310219717609728515625,
                          21076002824735985616104392981206548077528)<0
    return dict(full_base_negative_coefficient=str(coefficient),
                target=0,marks=(1,2))


def binary_excursion_audit():
    """All binary mixed starts, including vacancies, against raw ARW."""
    checked=blocks=0
    for case in rn.CASES:
        K=[[F(q) for q in row] for row in case["K"]]
        p=list(map(F,case["p"]))
        n=len(p)
        cache={}
        def step(state):
            nonlocal blocks
            i=next(i for i,a in enumerate(state) if a>0)
            B=sum(1<<j for j,a in enumerate(state) if a)
            key=(B,i)
            if key not in cache:
                start=(i,0)
                def walk_step(st):
                    z,mask=st
                    return {(j,mask | (1<<j) if j<n else mask):q
                            for j,q in enumerate(K[z]) if q}
                walks,_=rn.stopped([start],walk_step,
                    lambda st:st[0]==n or st[1]!=0 and
                    (st[0]==i or not B & (1<<st[0])))
                cache[key]=walks[start]
                blocks+=1
            row=defaultdict(F)
            asleep=list(state);asleep[i]=-1
            row[tuple(asleep)]+=p[i]
            for (j,visited),mass in cache[key].items():
                end=list(state);end[i]=0
                for z in range(n):
                    if visited & (1<<z):end[z]=1
                row[tuple(end)]+=(1-p[i])*mass
            assert sum(row.values()) == 1
            return dict(row)
        starts=list(product((-1,0,1),repeat=n))
        macro,_=rn.stopped(starts,step,lambda state:all(q<=0 for q in state))
        raw,_=rn.stopped(starts,lambda state:rn.arw_step(K,p,state),
                         lambda state:all(q<=0 for q in state))
        for start in starts:
            terminals=set(macro[start]) | set(raw[start])
            assert all(macro[start].get(end,F(0)) == raw[start].get(end,F(0))
                       for end in terminals)
            checked+=len(terminals)
    return dict(binary_excursion_terminal_entries=checked, excursion_blocks=blocks)


def directed_cycle_field_audit():
    """Circulation renewal and exact reduction to one-hole coefficients."""
    from math import prod
    checked=renewals=0
    for a in ([F(1,2),F(2,3),F(3,4)],[F(1),F(1),F(1,2)]):
        n=3
        p=[F(1,3),F(2,5),F(3,7)]
        K=[[F(0)]*(n+1) for _ in range(n)]
        for i in range(n):K[i][(i+1)%n]=a[i];K[i][n]=1-a[i]
        V=frozenset(range(n))
        sets=[frozenset(C) for C in subsets(range(n))]
        starts={(R,S):tuple(0 if i not in R else -1 if i in S else 1 for i in range(n))
                for R in sets for S in sets if S <= R}
        raw,_=rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                         lambda state:all(q<=0 for q in state))
        c=prod(a)
        for T in sets:
            moment={(R,S):sum(q for end,q in raw[start].items() if all(end[i]==-1 for i in T))
                    for (R,S),start in starts.items()}
            theta={(R,D):sum((-1)**len(D-S)*moment[R,S] for S in sets if S <= D)
                   for R in sets for D in sets if D <= R}
            mu=moment[V,frozenset()]
            for i in range(n):
                r=(1-p[i])/p[i]
                assert theta[V,frozenset((i,))] == r*(1-c)*(mu-moment[V-{i},frozenset()])
                checked+=1
                for D in sets:
                    if not D or i in D:continue
                    alpha=F(1);z=i
                    while z not in D:
                        alpha*=a[z];z=(z+1)%n
                    assert theta[V,D | {i}] == r*(theta[V,D]-(1-alpha)*theta[V-{i},D])
                    checked+=1
            rhs=prod(p)
            for i in range(n):
                E=frozenset(range(i))
                prefix=F(1)
                visited=set()
                for k in range(n):
                    z=(i+k)%n
                    if k:visited.add(z)
                    beta=prefix*(1-a[z])
                    rhs+=(1-p[i])*prod(p[j] for j in E)*beta*moment[V-{i},E-visited]
                    prefix*=a[z]
            assert (1-c*(1-prod(p)))*mu == rhs
            renewals+=1
    return dict(cycle_marked_field_entries=checked, cycle_first_loss_renewals=renewals)


def full_field_current_audit():
    """Exact edge currents and normalized potential from raw one-hole laws."""
    from math import prod
    edges=squares=reconstructions=0
    for case in rn.CASES:
        K=[[F(q) for q in row] for row in case["K"]]
        p=list(map(F,case["p"]))
        n=len(p);V=frozenset(range(n))
        sets=[frozenset(C) for C in subsets(range(n))]
        supports=[V]+[V-{i} for i in range(n)]
        starts={(R,S):tuple(0 if i not in R else -1 if i in S else 1 for i in range(n))
                for R in supports for S in sets if S <= R}
        raw,_=rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                         lambda state:all(q<=0 for q in state))
        minors={C:determinant([[F(i==j)-K[i][j] for j in sorted(C)] for i in sorted(C)])
                for C in sets}
        weight={H:prod((1-p[i])/p[i] for i in H)/minors[V-H] for H in sets}
        total=sum(weight.values());pi={H:q/total for H,q in weight.items()}
        h={};e={}
        for H in sets:
            U=sorted(V-H)
            G=rn.inverse([[F(i==j)-K[i][j] for j in U] for i in U])
            for i,row in zip(U,G):
                h[H,i]=sum(q*K[j][n] for j,q in zip(U,row))
                e[H,i]=minors[V-H]/minors[V-(H | {i})]
        for T in sets:
            moments={(R,S):sum(q for end,q in raw[start].items()
                              if all(end[i]==-1 for i in T))
                     for (R,S),start in starts.items()}
            coeff={(R,H):sum((-1)**len(H-S)*moments[R,S] for S in sets if S <= H)
                   for R in supports for H in sets if H <= R}
            gradient={}
            for H in sets:
                for i in V-H:
                    j=(1-p[i])*e[H,i]*h[H,i]*coeff[V-{i},H]
                    current=(1-p[i])*e[H,i]*coeff[V,H]-p[i]*coeff[V,H | {i}]
                    assert current == j
                    gradient[H,i]=h[H,i]*coeff[V-{i},H]/pi[H]
                    assert coeff[V,H]/pi[H]-coeff[V,H | {i}]/pi[H | {i}] == gradient[H,i]
                    edges+=1
            potential={frozenset():F(0)}
            for H in sets:
                if H:
                    i=max(H);B=H-{i}
                    potential[H]=potential[B]+gradient[B,i]
                for i in V-H:
                    for j in V-H-{i}:
                        if j<i:continue
                        assert gradient[H,i]+gradient[H | {i},j] == gradient[H,j]+gradient[H | {j},i]
                        squares+=1
            mean=sum(pi[H]*potential[H] for H in sets)
            for H in sets:
                assert coeff[V,H] == pi[H]*(1+mean-potential[H])
                reconstructions+=1
    return dict(full_field_edge_currents=edges, zero_circulation_squares=squares,
                potential_reconstruction_entries=reconstructions)


if __name__ == "__main__":
    print(proper_target_field_obstruction(),flush=True)
    print(full_field_current_audit(), flush=True)
    print(directed_cycle_field_audit(), flush=True)
    print(binary_excursion_audit(), flush=True)
    print(full_target_excursion_dual_audit(), flush=True)
