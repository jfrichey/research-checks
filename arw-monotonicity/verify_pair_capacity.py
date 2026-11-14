"""Exact finite audits of the N2 pair-to-triple capacity condition.

Standard library only. Raw ARW absorption, all Hall cuts, explicit rational
flows, and pair renewal/source/transience identities on two specified
four-site kernels. No parameter search.
These finite checks do not prove the arbitrary-graph condition.
One fixed mixed-contact obstruction records why all-active input matters.
AI-written research code; see response-renewal.md Section 26.
"""
from collections import deque
from fractions import Fraction as F
from itertools import combinations, product
import verify_response_renewal as rn


def inclusion_flow(pairs, triples, debts, supplies):
    """Return a rational allocation meeting every pair deficit."""
    source=('source',); sink=('sink',)
    residual={}; neighbors={}
    def edge(a,b,c):
        residual[a,b]=c; residual[b,a]=F(0)
        neighbors.setdefault(a,[]).append(b)
        neighbors.setdefault(b,[]).append(a)
    total=sum(debts.values(),F(0))
    for B in triples:
        edge(source,('triple',B),supplies[B])
        for C in pairs:
            if C<=B:
                edge(('triple',B),('pair',C),total)
    for C in pairs:
        edge(('pair',C),sink,debts[C])
    value=F(0)
    while value<total:
        parent={source:None}; queue=deque([source])
        while queue and sink not in parent:
            a=queue.popleft()
            for b in neighbors.get(a,[]):
                if b not in parent and residual[a,b]>0:
                    parent[b]=a; queue.append(b)
        assert sink in parent, 'infeasible pair capacity'
        path=[]; b=sink
        while parent[b] is not None:
            a=parent[b]; path.append((a,b)); b=a
        delta=min(residual[e] for e in path)
        for a,b in path:
            residual[a,b]-=delta; residual[b,a]+=delta
        value+=delta
    flow={(B,C):total-residual[('triple',B),('pair',C)]
          for B in triples for C in pairs if C<=B}
    assert all(w>=0 for w in flow.values())
    assert all(sum((w for (B,D),w in flow.items() if D==C),F(0))==debts[C]
               for C in pairs)
    assert all(sum((w for (D,C),w in flow.items() if D==B),F(0))<=supplies[B]
               for B in triples)
    return flow


def killed_pair_laws(K,p,pairs):
    """Solve each killed domain once, for all its two-root inputs.

    On exact pair outputs no particle was lost, so these entries are also
    the original-domain no-visit subprobabilities. In particular the A-to-A
    entry with x killed is kappa; no reciprocal diagonal formula is used.
    """
    n=len(p); laws={}; count=0
    for z in range(n):
        domain=[i for i in range(n) if i!=z]
        killed_K=[[K[i][j] for j in domain]+[K[i][n]+K[i][z]]
                  for i in domain]
        killed_p=[p[i] for i in domain]
        assert all(sum(row)==1 for row in killed_K)
        starts={A:tuple(int(i in A) for i in domain)
                for A in pairs if z not in A}
        raw,states=rn.stopped(
            starts.values(),lambda state:rn.arw_step(killed_K,killed_p,state),
            lambda state:all(q<=0 for q in state))
        count+=states
        laws[z]={A:{frozenset(domain[i] for i,q in enumerate(end) if q==-1):weight
                    for end,weight in raw[start].items()}
                 for A,start in starts.items()}
        assert all(sum(law.values())==1 for law in laws[z].values())
    return laws,count


def pair_renewal_audit(A,x,p,pairs,triples,lo,hi,W,killed,h):
    """Check both choices of missing root against the same absorption laws."""
    r=[(1-q)/q for q in p]
    ell=[lo.get(C,F(0)) for C in pairs]
    upper=[hi.get(C,F(0)) for C in pairs]
    tau=[hi.get(B,F(0)) for B in triples]
    kappa=killed[x][A].get(A,F(0))
    a,b=sorted(A)
    assert kappa>=p[a]*p[b]>0,(A,x,kappa)
    minimum=None
    for policy,choose in [('smallest',min),('largest',max)]:
        selected=[x if C==A else choose(A-C) for C in pairs]
        assert all(z not in C and z in A|{x} and (C==A or z in A)
                   for C,z in zip(pairs,selected))
        P=[[r[z]*W[D,z].get(C,F(0)) if z in D else F(0)
            for D in pairs] for C,z in zip(pairs,selected)]
        Q=[[r[z]*W[B,z].get(C,F(0)) if z in B else F(0)
            for B in triples] for C,z in zip(pairs,selected)]
        assert all(q>=0 for row in P+Q for q in row)
        assert all(P[i][i]==0 for i in range(len(pairs)))
        for i,(C,z) in enumerate(zip(pairs,selected)):
            Pell=sum((q*mass for q,mass in zip(P[i],ell)),F(0))
            Pu=sum((q*mass for q,mass in zip(P[i],upper)),F(0))
            Qtau=sum((q*mass for q,mass in zip(Q[i],tau)),F(0))
            assert Pell==ell[i]-(kappa if C==A else 0),(A,x,policy,C,'lower')
            assert Pu+Qtau==upper[i],(A,x,policy,C,'upper')

            # h is the pair layer of the uniform mixture of all pair-root laws.
            # Its exact source is independently available from the killed laws.
            Ph=sum((q*h[D] for q,D in zip(P[i],pairs)),F(0))
            source=sum((killed[z][D].get(C,F(0)) for D in pairs if z not in D),
                       F(0))/len(pairs)
            defect=h[C]-Ph
            assert defect==source,(A,x,policy,C,'mixture source')
            c,d=sorted(C)
            assert defect>=p[c]*p[d]/len(pairs)>0,(A,x,policy,C,'transience')
            assert 0<h[C]<=1
            gap=defect/h[C]
            assert 0<gap<=1
            minimum=gap if minimum is None else min(minimum,gap)
    return minimum


def audit(K,p,archived=None):
    n=len(p); vertices=frozenset(range(n))
    pairs=list(map(frozenset,combinations(range(n),2)))
    triples=list(map(frozenset,combinations(range(n),3)))
    starts={A:tuple(int(z in A) for z in range(n)) for A in pairs+triples}
    raw,count=rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                         lambda state:all(q<=0 for q in state))
    laws={A:{frozenset(z for z,q in enumerate(end) if q==-1):weight
             for end,weight in raw[start].items()} for A,start in starts.items()}
    if archived is not None:
        subsets=[frozenset(i for i,b in enumerate(bits) if b)
                 for bits in archived['subsets']]
        saved={A:{B:F(q) for B,q in zip(subsets,row) if F(q)}
               for A,row in zip(subsets,archived['stabilization_matrix'])}
        assert all(laws[A]==saved[A] for A in starts)
    # These continuations and killed laws are shared by every root comparison
    # and both policies; the original-domain raw laws above are not recomputed.
    _,_,_,raw_W=rn.positive_operator(K,p,3)
    W={(frozenset(B),z):{frozenset(C):q for C,q in law.items()}
       for (B,z),law in raw_W.items()}
    killed,killed_count=killed_pair_laws(K,p,pairs)
    h={C:sum((laws[A].get(C,F(0)) for A in pairs),F(0))/len(pairs)
       for C in pairs}
    families=[D for size in range(1,len(pairs)+1)
              for D in combinations(pairs,size)]
    cuts=0; comparisons=0; minimum=None; minimum_gap=None
    for A in pairs:
        for x in sorted(vertices-A):
            lo=laws[A]; hi=laws[A|{x}]
            assert sum(lo.values())==sum(hi.values())==1
            assert all(q>=0 and len(C)<=2 for C,q in lo.items())
            assert all(q>=0 and len(C)<=3 for C,q in hi.items())
            debts={C:max(F(0),lo.get(C,F(0))-hi.get(C,F(0))) for C in pairs}
            supplies={B:hi.get(B,F(0)) for B in triples}
            for D in families:
                demand=sum((debts[C] for C in D),F(0))
                supply=sum((supplies[B] for B in triples
                            if any(C<=B for C in D)),F(0))
                slack=supply-demand
                assert slack>=0,(A,x,D,slack)
                minimum=slack if minimum is None else min(minimum,slack)
                cuts+=1
            inclusion_flow(pairs,triples,debts,supplies)
            gap=pair_renewal_audit(A,x,p,pairs,triples,lo,hi,W,killed,h)
            minimum_gap=gap if minimum_gap is None else min(minimum_gap,gap)
            comparisons+=1
    return dict(comparisons=comparisons,hall_cuts=cuts,
                raw_transient_states=count,minimum_slack=str(minimum),
                renewal_policies=['smallest','largest'],
                renewal_systems=2*comparisons,kappa_checks=comparisons,
                killed_transient_states=killed_count,
                minimum_contraction_gap=str(minimum_gap))


def mixed_contact_obstruction():
    """A mixed-state failure, not a counterexample to all-active Hall or N2."""
    K=[[F(0),F(1,2),F(0),F(1,2)],
       [F(1,2),F(0),F(1,2),F(0)],
       [F(0),F(1,2),F(0),F(1,2)]]
    p=[F(1,2)]*3
    lower=(1,-1,0); upper=(2,-1,0)
    raw,count=rn.stopped([lower,upper],lambda state:rn.arw_step(K,p,state),
                         lambda state:all(q<=0 for q in state))
    # Joint occupation of {0,1}: the upper event includes the full triple.
    values={start:sum((q for end,q in law.items() if end[0]==end[1]==-1),F(0))
            for start,law in raw.items()}
    assert values[lower]==F(45,82)
    assert values[upper]==F(45323,99876)
    increment=values[upper]-values[lower]
    assert increment==F(-9487,99876)<0
    return dict(lower=str(values[lower]),upper=str(values[upper]),
                increment=str(increment),raw_transient_states=count,
                scope='mixed inputs only')


def union_restoration_obstruction():
    """Stopping at either of two pair events destroys fixed-root confluence."""
    from verify_stack_exploration import stabilize
    from verify_root_restoration import counts

    base={(z,i):(1,-1) for z in range(3) for i in range(1,20)}
    base.update({(0,1):(0,2),(1,1):(0,2),
                 (2,1):(1,-1),(2,2):(0,-1)})
    def event(B):
        return 0 in B and bool(B & {1,2})
    assert stabilize(base,[0,1])[0]==set()
    outcomes=[]
    for order in [(0,), (1,0)]:
        book=dict(base)
        for z in order:
            result=stabilize(book,[0,1])
            assert not event(result[0]) and z not in result[0]
            index=max(i for x,i in result[1] if x==z)
            assert book[z,index][0]==0
            book[z,index]=1,book[z,index][1]
        result=stabilize(book,[0,1])
        assert event(result[0])
        outcomes.append((sorted(result[0]),dict(counts(result))))
    assert outcomes==[([0,2],{1:1}),([0,1],{})]
    return dict(initial=[0,1],pair_family=[[0,1],[0,2]],
                restoration_orders=[[0],[1,0]],outcomes=outcomes,
                scope='event-stopped restoration only; no Hall counterexample')


def killed_sleep_insertions(K,p,v):
    """Raw no-v-visit addition and signed sleep insertion on V without v."""
    n=len(p)
    sites=[z for z in range(n) if z!=v]
    supports=[frozenset(z for j,z in enumerate(sites) if mask>>j&1)
              for mask in range(1<<len(sites))]
    pos={C:i for i,C in enumerate(supports)}
    size=len(supports)
    identity=[[F(i==j) for j in range(size)] for i in range(size)]
    starts={}
    for C in supports:
        for z in sites:
            state=[-int(y in C) for y in range(n)]
            state[z]=2 if z in C else 1
            starts[C,z]=tuple(state)
    raw,count=rn.stopped(starts.values(),lambda s:rn.arw_step(K,p,s),
                         lambda s:s[v]>0 or all(q<=0 for q in s))
    A={z:[[raw[starts[C,z]].get(tuple(-int(y in D) for y in range(n)),F(0))
           for D in supports] for C in supports] for z in sites}
    J={z:[[K[z][n]*identity[i][j]+sum(K[z][y]*A[y][i][j] for y in sites)
           for j in range(size)] for i in range(size)] for z in range(n)}
    B={z:[[(A[z][i][j]-(1-p[z])*J[z][i][j])/p[z]
           for j in range(size)] for i in range(size)] for z in sites}
    for z in sites:
        square=rn.multiply(J[z],J[z])
        r=(1-p[z])/p[z]
        for C in supports:
            i=pos[C]
            if z not in C:
                assert B[z][i]==[F(D==C|{z}) for D in supports]
            else:
                expected=[(1-r)*J[z][i][j]+r*square[pos[C-{z}]][j]
                          for j in range(size)]
                assert B[z][i]==expected
            if (1-p[z])*(1-K[z][z])<=p[z]:
                assert min(B[z][i])>=0 and sum(B[z][i])<=1
        for y in sites:
            assert rn.multiply(B[z],B[y])==rn.multiply(B[y],B[z])
            assert rn.multiply(B[z],A[y])==rn.multiply(A[y],B[z])
    return supports,pos,A,J,B,count


def high_sleep_budget_audit():
    """Exact killed-kernel budget and N2 decomposition on fixed saved kernels.

    Positivity of each B_z needs only p_z>=1/2, even when other p's are
    smaller. The scalar decomposition is tested when every off-target
    probability meets the bound. No parameter sweep is used.
    """
    cases=[]
    for case_number in [2,3,4]:
        case=rn.CASES[case_number]
        K=[[F(q) for q in row] for row in case['K']]
        p=list(map(F,case['p'])); n=len(p)
        full=[frozenset(z for z,q in enumerate(bits) if q)
              for bits in case['subsets']]
        full_pos={C:i for i,C in enumerate(full)}
        P=[[F(q) for q in row] for row in case['stabilization_matrix']]
        _,_,_,W=rn.positive_operator(K,p,3)
        budgets=decompositions=bulk_decompositions=raw_states=0
        for v in range(n):
            supports,pos,A,J,B,count=killed_sleep_insertions(K,p,v)
            raw_states+=count
            sites=[z for z in range(n) if z!=v]
            nu=J[v][pos[frozenset()]]
            assert all(not q for C,q in zip(supports,nu) if len(C)>1)
            w={}
            for C in supports:
                w[C]=sum(J[v][pos[C]])
                row=nu[:]
                for z in sorted(C):
                    row=rn.multiply([row],B[z])[0]
                assert sum(row)==w[C]
                if len(C)<=2:
                    assert w[C]==sum(W[tuple(sorted(C|{v})),v].values(),F(0))
            for a,b in combinations(sites,2):
                da=[1-sum(row) for row in B[a]]
                db=[1-sum(row) for row in B[b]]
                rhs=(nu[pos[frozenset({a})]]*
                     sum(q*d for q,d in zip(B[a][pos[frozenset({a})]],db))+
                     nu[pos[frozenset({b})]]*
                     sum(q*d for q,d in zip(B[b][pos[frozenset({b})]],da)))
                slack=w[frozenset({a})]+w[frozenset({b})]-w[frozenset()]-w[frozenset({a,b})]
                assert slack==rhs
                if min(p[a],p[b])>=F(1,2):
                    assert slack>=0
                    budgets+=1
            if any(p[z]<F(1,2) for z in sites):
                continue
            r=(1-p[v])/p[v]
            scale=1/(1+r*w[frozenset()])
            weights={C|{v}:r*scale*(w[frozenset()]-w[C])
                     for C in w if C}
            def no_hit(roots):
                if v in roots:
                    return F(0)
                row=[F(not C) for C in supports]
                for z in sorted(roots):
                    row=rn.multiply([row],A[z])[0]
                return sum(row)
            pairs=[C for C in full if len(C)==2 and v in C]
            larger=[C for C in full if len(C)>=3 and v in C]
            for roots in [C for C in full if len(C)==2]:
                batches=[frozenset({x}) for x in range(n) if x not in roots]
                if n-len(roots)>1:
                    batches.append(frozenset(range(n))-roots)
                for batch in batches:
                    upper=roots|batch
                    lo=P[full_pos[roots]]; hi=P[full_pos[upper]]
                    terms=[scale*(no_hit(roots)-no_hit(upper))]
                    for C in pairs:
                        gain=sum(q for D,q in zip(full,hi) if C<=D)-lo[full_pos[C]]
                        assert gain>=0
                        terms.append(weights[C]*gain)
                    for D in larger:
                        budget=weights[D]-sum(weights[C] for C in pairs if C<=D)
                        assert budget>=0
                        terms.append(budget*hi[full_pos[D]])
                    increment=sum(q-h for C,q,h in zip(full,hi,lo) if v in C)
                    assert min(terms)>=0 and increment==sum(terms)
                    if len(batch)==1:
                        decompositions+=1
                    else:
                        bulk_decompositions+=1
        cases.append(dict(case=case_number,raw_transient_states=raw_states,
                          two_background_budgets=budgets,N2_decompositions=decompositions,
                          bulk_decompositions=bulk_decompositions))

    # A low-sleep failure of the coefficient budget, not of N2 itself.
    K=[[F(0),F(1,2),F(0),F(1,2)],
       [F(0),F(0),F(1),F(0)],
       [F(1),F(0),F(0),F(0)]]
    p=[F(1,2),F(1,4),F(1,2)]
    supports,pos,_,J,_,_=killed_sleep_insertions(K,p,0)
    values=[sum(J[0][pos[frozenset(C)]]) for C in [(),(1,),(2,),(1,2)]]
    assert values==[F(13,16),F(9,16),F(5,8),F(1,2)]
    assert values[1]+values[2]-values[0]-values[3]==F(-1,8)
    return dict(high_sleep_audits=cases,low_sleep_w_values=list(map(str,values)),
                low_sleep_budget='-1/8',scope='coefficient failure only at low sleep')


def quantitative_pair_audit():
    """Entrance-mass bound with its exact remainder and one scoped failure."""
    from verify_creator_potentials import idla_addition_matrices
    from verify_feedback_cone import binary_stabilization
    from verify_creation_cone import simplex

    totals=[]
    for number in [2,3,4]:
        case=rn.CASES[number]
        K=[[F(q) for q in row] for row in case['K']]
        p=list(map(F,case['p'])); n=len(p)
        supports=[frozenset(i for i,q in enumerate(bits) if q)
                  for bits in case['subsets']]
        pos={C:i for i,C in enumerate(supports)}
        P=[[F(q) for q in row] for row in case['stabilization_matrix']]
        D=idla_addition_matrices(K,supports)
        plus=[rn.multiply(D[x],P) for x in range(n)]
        pairs=[C for C in supports if len(C)==2]
        beta={}
        for C in pairs:
            boundary=[B for B in supports if C<=B and len(B)<=3]
            Q=[[P[pos[E]][pos[B]] for B in boundary] for E in boundary]
            m=dict(zip(boundary,map(sum,rn.inverse(Q))))
            beta[C]={B:m[B]/m[C] for B in boundary if B!=C}
            assert all(0<q<=1 for q in beta[C].values())
        checks=0
        for roots in pairs:
            for x in range(n):
                lo=P[pos[roots]]; hi=plus[x][pos[roots]]
                for C in pairs:
                    alpha=lo[pos[C]]/P[pos[C]][pos[C]]
                    component=plus[x][pos[C]]
                    residual=[q-alpha*h for q,h in zip(hi,component)]
                    assert min(residual)>=0
                    deficit=lo[pos[C]]-hi[pos[C]]
                    exact=sum(b*alpha*component[pos[B]] for B,b in beta[C].items())-residual[pos[C]]
                    bound=sum(b*hi[pos[B]] for B,b in beta[C].items())
                    assert deficit==exact<=bound
                    checks+=1
        totals.append(dict(case=number,exact_remainder_checks=checks))

    K=[[F(0),F(0),F(1,2),F(1,2)],
       [F(1,2),F(0),F(0),F(1,2)],
       [F(15,16),F(3,64),F(0),F(1,64)]]
    p=[F(1,1001)]*3
    supports,P,E,D=binary_stabilization(K,p)
    pos={C:i for i,C in enumerate(supports)}
    _,_,_,W=rn.positive_operator(K,p,3)
    full=frozenset(range(3)); target=0
    w=lambda C:sum(W[tuple(sorted(C)),target].values(),F(0))
    base=w({target}); budget=base-w(full); rare_beta={}
    for C in [frozenset({0,1}),frozenset({0,2})]:
        boundary=[C,full]
        Q=[[P[pos[E]][pos[B]] for B in boundary] for E in boundary]
        masses=list(map(sum,rn.inverse(Q)))
        rare_beta[C]=masses[1]/masses[0]
        budget-=(base-w(C))*masses[1]/masses[0]
    assert budget==F(-266112409230694467411733683,
                     56456560891831241242814096384)

    # Retaining only the triple part of the minorization remainder also
    # fails on C_2. This fixed exact feasible witness is not a root law.
    source=2
    addition=rn.multiply(rn.multiply(E,D[source]),P)
    plus=rn.multiply(D[source],P)
    states=[C for C in supports if len(C)<=2]
    constraints=[]
    for z in range(3):
        for C in states:
            if z in C:
                continue
            constraints.append([
                F(B==C) if z not in B else
                -(1-p[z])/p[z]*W[tuple(sorted(B)),z].get(tuple(sorted(C)),F(0))
                for B in states])
    g=[]
    for C in states:
        value=(base-w(full))*addition[pos[C]][pos[full]]
        if C in rare_beta:
            value-=(base-w(C))*rare_beta[C]*plus[pos[C]][pos[full]]/P[pos[C]][pos[C]]
        g.append(value)
    lp=[[-q for q in row] for row in constraints]+[[F(1)]*len(states)]
    rhs=[F(0)]*len(constraints)+[F(1)]
    value,mu,_,_=simplex(lp,rhs,[-q for q in g])
    assert min(mu)>=0 and sum(mu)<=1
    assert all(sum(q*h for q,h in zip(row,mu))>=0 for row in constraints)
    assert value==-sum(q*h for q,h in zip(g,mu))==F(
        8811470822597717846819589644771241224919,
        11046683866931891518105167894089536955986428385262336)
    assert value>0
    occupation_response=sum(weight*(sum(addition[pos[C]][pos[B]]
                            for B in supports if target in B)-int(target in C))
                            for C,weight in zip(states,mu))
    assert occupation_response>=0
    return dict(quantitative_pair_audits=totals,
                rare_sleep_weighted_triangle_slack=str(budget),
                source_weighted_cone_budget_slack=str(-value),
                cone_witness_occupation_response=str(occupation_response),
                cone_witness=[(sorted(C),str(q)) for C,q in zip(states,mu) if q],
                scope='pointwise combined budget fails; individual pair bound holds')


def released_resident_audit():
    """Defects stay nonnegative even when killed sleep insertion is signed."""
    K=[[F(0),F(1,2),F(0),F(1,2)],
       [F(0),F(0),F(1),F(0)],
       [F(1),F(0),F(0),F(0)]]
    loop_K=[row[:] for row in K]
    loop_K[1][1],loop_K[1][2]=F(3,4),F(1,4)
    audits=[]
    for label,kernel,p,expected in [
            ('nonterminal_low',K,[F(1,2),F(1,4),F(1,2)],F(-1,8)),
            ('terminal_low',K,[F(1,2),F(2,3),F(1,4)],F(1,24)),
            ('effective_high',loop_K,[F(1,2),F(1,4),F(1,2)],F(1,28))]:
        supports,pos,_,J,B,_=killed_sleep_insertions(kernel,p,0)
        checked=0
        for z in (1,2):
            square=rn.multiply(J[z],J[z])
            r=(1-p[z])/p[z]
            for C in supports:
                defect=1-sum(B[z][pos[C]])
                if z not in C:
                    assert defect==0
                    continue
                j1=sum(J[z][pos[C]])
                j2=sum(square[pos[C-{z}]])
                assert 0<=j2<=j1<=1
                assert defect==1-j1+r*(j1-j2)>=0
                checked+=1
        w=lambda C:sum(J[0][pos[frozenset(C)]])
        slack=w({1})+w({2})-w(set())-w({1,2})
        assert slack==expected
        if label=='nonterminal_low':
            assert 1-sum(B[1][pos[frozenset({1})]])==2
        elif label=='terminal_low':
            # Site 2 has low sleep but cannot reach 1 before the target.
            d1=[1-sum(row) for row in B[1]]
            assert sum(q*d for q,d in zip(B[2][pos[frozenset({2})]],d1))==0
        else:
            assert p[1]/(1-(1-p[1])*kernel[1][1])==F(4,7)
            assert all(min(row)>=0 and sum(row)<=1 for row in B[1])
        audits.append(dict(case=label,p=list(map(str,p)),resident_release_inequalities=checked,
                           two_background_budget=str(slack)))
    return dict(released_resident_audits=audits)


def resident_release_transfer_audit():
    """A one-site odds change is a geometric number of base-law releases."""
    from verify_feedback_cone import binary_stabilization
    case=rn.CASES[4]
    K=[[F(q) for q in row] for row in case['K']]
    p=[F(1,2)]*4; n=4; site=1
    supports,P,E,D=binary_stabilization(K,p)
    pos={C:i for i,C in enumerate(supports)}; size=len(supports)
    A=[rn.multiply(rn.multiply(E,mat),P) for mat in D]
    J=[[K[site][n]*F(i==j)+sum(K[site][y]*A[y][i][j] for y in range(n))
        for j in range(size)] for i in range(size)]
    release=[J[pos[C-{site}]] if site in C else [F(C==B) for B in supports]
             for C in supports]
    assert all(min(row)>=0 and sum(row)==1 for row in release)
    assert all(not weight for C,row in zip(supports,release)
               for B,weight in zip(supports,row) if len(B)>len(C))
    changed=p[:]; changed[site]=F(1,4)
    _,new_P,new_E,_=binary_stabilization(K,changed)
    t=p[site]*((1-changed[site])/changed[site]-(1-p[site])/p[site])
    matrix=[[F(i==j)+t*(F(i==j)-release[i][j])
             for j in range(size)] for i in range(size)]
    assert rn.multiply(new_E,P)==matrix
    transfer=rn.inverse(matrix)
    assert transfer==rn.multiply(E,new_P)
    assert all(min(row)>=0 and sum(row)==1 for row in transfer)
    assert t/(1+t)==F(1,2)
    return dict(resident_release_transfer_entries=size*size,site=site,
                old_p=str(p[site]),new_p=str(changed[site]),
                geometric_continue_probability=str(t/(1+t)))


def actual_source_budget_audit():
    """Fixed cases passing a source budget that fails on other actual inputs."""
    from verify_feedback_cone import binary_stabilization
    archived=rn.CASES[4]
    cases=[('RN12_cycle',
            [[F(31,32) if j==(i+1)%4 else F(0) for j in range(4)]+[F(1,32)]
             for i in range(4)],[F(1,2),F(1,2),F(1,8),F(1,8)]),
           ('archived_four_rare',[[F(q) for q in row] for row in archived['K']],
            [F(1,1001)]*4),
           ('RN22_rare_three',
            [[F(0),F(0),F(1,2),F(1,2)],
             [F(1,2),F(0),F(0),F(1,2)],
             [F(15,16),F(3,64),F(0),F(1,64)]],[F(1,1001)]*3)]
    audits=[]
    for name,K,p in cases:
        n=len(p)
        supports,P,_,D=binary_stabilization(K,p)
        pos={C:i for i,C in enumerate(supports)}
        plus=[rn.multiply(mat,P) for mat in D]
        pairs=[C for C in supports if len(C)==2]
        triples=[B for B in supports if len(B)==3]
        mass={}
        for C in pairs:
            boundary=[B for B in supports if C<=B and len(B)<=3]
            Q=[[P[pos[E]][pos[B]] for B in boundary] for E in boundary]
            mass[C]=dict(zip(boundary,map(sum,rn.inverse(Q))))
        _,_,_,W=rn.positive_operator(K,p,3)
        checks=0; minimum=None
        for v in range(n):
            def w(C):
                return sum(W[tuple(sorted(C)),v].values(),F(0))
            base=w({v})
            weights={C:base-w(C) for C in pairs+triples if v in C}
            for roots in pairs:
                for x in range(n):
                    if x in roots:
                        continue
                    lo=P[pos[roots]]; upper=plus[x][pos[roots]]
                    for B in triples:
                        if v not in B:
                            continue
                        slack=(weights[B]*upper[pos[B]]-
                               sum(weights[C]*lo[pos[C]]*mass[C][B]*
                                   plus[x][pos[C]][pos[B]]
                                   for C in pairs if v in C and C<=B))
                        assert slack>=0,(name,roots,x,v,B,slack)
                        minimum=slack if minimum is None else min(minimum,slack)
                        checks+=1
        audits.append(dict(case=name,source_dependent_triple_budgets=checks,
                           minimum_slack=str(minimum)))
    assert sum(row['source_dependent_triple_budgets'] for row in audits)==297
    return dict(actual_source_budget_audits=audits,scope='finite evidence only')


def rare_source_suppression_audit():
    """Uniform-odds source suppression on the fixed RN22 rare-three kernel.

    The general first-contact identity is
      B_y-A_y = L_y(Insert_y-E_y A_y),
    where E_y embeds supports avoiding y. On creation laws it follows by
    holding the extra active/sleeping resident until first contact; the
    invertible binary creation matrix then gives the operator identity.
    Thus f(y)=mu_(R+delta_y)(C) is harmonic on R and subharmonic off
    C without R, with f(sink)=0 when |C|=|R|+1. Stopping the ordinary
    walk at that boundary or the sink gives Harnack, including sources
    already in R.

    For A={0,a}, x=b, C={0,b}, Harnack bounds alpha by h_a(b), and the
    source-to-full ratio by h_b(a). Here there are exactly three sites,
    so the latter ratio equals kappa=h_b(a) by surplus IDLA. Both products
    h_a(b)h_b(a) equal 3/136. Small positive quadratic cores then prove
    the source-weighted triple budget for every common odds r>0. The
    exact-pair residual is retained separately. No alpha gap polynomial,
    large occupation certificate, SymPy dependency, or parameter sweep.
    """
    from verify_feedback_cone import binary_stabilization

    K=[[F(0),F(0),F(1,2),F(1,2)],
       [F(1,2),F(0),F(0),F(1,2)],
       [F(15,16),F(3,64),F(0),F(1,64)]]
    p=[F(1,1001)]*3
    odds=F(1000)
    full=frozenset(range(3))
    supports,P,E,D=binary_stabilization(K,p)
    pos={C:i for i,C in enumerate(supports)}
    size=len(supports)
    plus=[rn.multiply(mat,P) for mat in D]
    additions=[rn.multiply(E,mat) for mat in plus]
    _,_,_,W=rn.positive_operator(K,p,3)

    def harmonic_extension(boundary,values,sink=F(0)):
        remaining=sorted(full-boundary)
        inverse=rn.inverse([[F(i==j)-K[i][j] for j in remaining]
                            for i in remaining])
        rhs=[K[i][3]*sink+sum(K[i][j]*values[j] for j in boundary)
             for i in remaining]
        result=dict(values)
        result.update({i:sum(q*b for q,b in zip(row,rhs))
                       for i,row in zip(remaining,inverse)})
        return result

    hit={j:harmonic_extension(frozenset({j}),{j:F(1)}) for j in (1,2)}
    assert hit[2][1]==F(1,4) and hit[1][2]==F(3,34)
    lam=hit[2][1]*hit[1][2]
    assert lam==F(3,136)

    # Independent raw W rows audit the full first-contact matrix identity.
    # Insert_y is only applied on the smaller space of supports avoiding y.
    deletion={}
    operator_entries=0
    for y in range(3):
        avoided=[C for C in supports if y not in C]
        J=[[K[y][3]*F(i==j)+sum(K[y][z]*additions[z][i][j] for z in range(3))
            for j in range(size)] for i in range(size)]
        B=[[(additions[y][i][j]-(1-p[y])*J[i][j])/p[y]
            for j in range(size)] for i in range(size)]
        L=[[F(C==T) if y not in C else
            -odds*W[tuple(sorted(C)),y].get(tuple(sorted(T)),F(0))
            for T in avoided] for C in supports]
        insert_minus_add=[[F(T==C|{y})-additions[y][pos[C]][pos[T]]
                           for T in supports] for C in avoided]
        lhs=[[B[i][j]-additions[y][i][j] for j in range(size)]
             for i in range(size)]
        assert lhs==rn.multiply(L,insert_minus_add)
        assert lhs==[[odds*(additions[y][i][j]-J[i][j]) for j in range(size)]
                      for i in range(size)]
        assert rn.multiply(additions[y],additions[y])==rn.multiply(additions[y],J)
        no_visit=rn.multiply(P,L)
        assert all(q>=0 for row in no_visit for q in row)
        assert all(not any(no_visit[pos[R]]) for R in supports if y in R)
        deletion[y]=L
        operator_entries+=size*size
    assert operator_entries==192

    # R+delta_y includes a duplicate when y is in R; it is not set union.
    harnack_checks=0
    exact_support_checks=0
    for R in supports:
        for C in supports:
            if len(C)==len(R)+1:
                boundary=C-R
                f={y:plus[y][pos[R]][pos[C]] for y in range(3)}
                assert P[pos[R]][pos[C]]==0 and boundary
                bound=harmonic_extension(boundary,{y:f[y] for y in boundary})
                for y in range(3):
                    mean=sum(K[y][z]*f[z] for z in range(3))
                    if y in R:
                        assert f[y]==mean
                    if y not in boundary:
                        assert f[y]<=mean
                    assert f[y]<=bound[y]
                    harnack_checks+=1
            if C<=R:
                for y in range(3):
                    assert plus[y][pos[R]][pos[C]]<=P[pos[R]][pos[C]]
                    exact_support_checks+=1
    assert harnack_checks==45 and exact_support_checks==81

    # Polynomials below have ascending rational coefficients. Their common
    # denominator is (1+r)^2; no polynomial fitting or numeric sampling.
    def add(*polynomials):
        return tuple(sum((poly[i] if i<len(poly) else F(0)) for poly in polynomials)
                     for i in range(max(map(len,polynomials))))

    def scale(poly,factor):
        return tuple(F(factor)*q for q in poly)

    def multiply(left,right):
        out=[F(0)]*(len(left)+len(right)-1)
        for i,a in enumerate(left):
            for j,b in enumerate(right):
                out[i+j]+=a*b
        return tuple(out)

    def evaluate(poly,r):
        return sum(q*r**i for i,q in enumerate(poly))

    one_plus_r=(F(1),F(1))
    rpoly=(F(0),F(1))
    # c=p+q/2 is survival of one active at site 1 before hitting target 0.
    c_numerator=(F(1),F(1,2))
    r_times_c=multiply(rpoly,c_numerator)
    r_times_denominator=multiply(rpoly,one_plus_r)
    # F_ab: start active 2, with initial sleepers a at 1 and b at 2.
    # Successful emissions from 2 go to sink (1/64) or site 1 (3/64).
    f00=add(one_plus_r,scale(add(r_times_denominator,scale(r_times_c,3)),F(1,64)))
    f10=add(one_plus_r,scale(add(r_times_denominator,scale(r_times_c,F(3,2))),F(1,64)))
    f01=add(scale(add(one_plus_r,scale(c_numerator,3)),F(1,64)),
            scale(add(scale(r_times_denominator,2),scale(r_times_c,21)),F(1,8192)))
    f11=add(scale(add(one_plus_r,scale(c_numerator,F(3,2))),F(1,64)),
            scale(add(scale(r_times_denominator,4),scale(r_times_c,21)),F(1,16384)))
    pair01=frozenset({0,1}); pair02=frozenset({0,2})
    # u_C=w_0({0})-w_0(C); the forced first jump reaches 2 with probability 1/2.
    U={pair01:scale(add(f00,scale(f10,-1)),32768),
       pair02:scale(add(f00,scale(f01,-1)),32768),
       full:scale(add(f00,scale(f11,-1)),32768)}
    assert U[pair01]==(0,768,384)
    assert U[pair02]==(30720,33444,1230)
    assert U[full]==(31488,33870,1251)
    positive_cores={}
    for x,A,C in [(2,pair01,pair02),(1,pair02,pair01)]:
        core=add(U[full],scale(U[A],-1),scale(U[C],-lam))
        expected=(scale(multiply((2,1),(2176,71)),F(3,17)) if x==1 else
                  scale((698368,733590,19037),F(3,68)))
        assert core==expected and all(q>0 for q in core)
        positive_cores[x]=core

    # One exact raw solve audits both original roots and repeated-source
    # components, separately from the binary-inverse construction above.
    starts=[(1,1,0),(1,0,1),(1,1,1),(1,0,2),(1,2,0)]
    raw,raw_states=rn.stopped(starts,lambda state:rn.arw_step(K,p,state),
                              lambda state:all(q<=0 for q in state))
    assert raw_states==28
    terminal={C:tuple(-int(i in C) for i in range(3)) for C in supports}
    raw_rows={start:[raw[start].get(terminal[C],F(0)) for C in supports]
              for start in starts}
    upper=P[pos[full]]
    assert raw_rows[1,1,1]==upper
    tau=upper[pos[full]]
    w=lambda C:sum(W[tuple(sorted(C)),0].values(),F(0))
    w0=w(frozenset({0}))
    u={C:w0-w(C) for C in (pair01,pair02,full)}
    assert all(u[C]==evaluate(U[C],odds)/(65536*(1+odds)**2) for C in U)
    beta={C:(P[pos[C]][pos[C]]-upper[pos[C]])/tau for C in (pair01,pair02)}
    assert all(0<b<=1 for b in beta.values())
    negative_bracket=u[full]-sum(u[C]*beta[C] for C in beta)
    assert negative_bracket==F(-266112409230694467411733683,
                               56456560891831241242814096384)<0
    scale_v=odds/(1+odds*w0)
    cases=[]
    for x,a,A,C,duplicate in [(2,1,pair01,pair02,(1,0,2)),
                              (1,2,pair02,pair01,(1,2,0))]:
        lower=P[pos[A]]
        assert raw_rows[tuple(int(i in A) for i in range(3))]==lower
        component=plus[x][pos[C]]
        assert raw_rows[duplicate]==component
        self_mass=P[pos[C]][pos[C]]
        alpha=lower[pos[C]]/self_mass
        bar=hit[x][a]
        kappa=hit[a][x]
        assert 0<=alpha<=bar and bar*kappa==lam
        ratio=component[pos[full]]/tau
        assert ratio<=kappa and ratio==kappa  # Equality uses exactly three sites.
        assert component==[(1-kappa)*q+kappa*t for q,t in zip(P[pos[C]],upper)]
        sigma=[q-alpha*t for q,t in zip(lower,P[pos[C]])]
        residual=rn.multiply([sigma],additions[x])[0]
        assert min(sigma)>=0 and min(residual)>=0
        assert residual==[q-alpha*t for q,t in zip(upper,component)]
        assert not any(rn.multiply([lower],deletion[0])[0])
        exact_pair=residual[pos[C]]
        triple_piece=beta[C]*residual[pos[full]]
        assert residual[pos[full]]==(1-alpha*kappa)*tau
        assert exact_pair==(1-alpha)*self_mass-triple_piece
        assert exact_pair+triple_piece==(1-alpha)*self_mass
        core=evaluate(positive_cores[x],odds)/(65536*(1+odds)**2)
        terms=[core,u[A]*(1-beta[A]),u[C]*kappa*(bar-alpha),
               u[C]*alpha*kappa*(1-beta[C])]
        triple_budget=u[full]-u[A]*beta[A]-u[C]*beta[C]*alpha*kappa
        assert min(terms)>=0 and core>0 and triple_budget==sum(terms)>0
        gain=sum(q-h for B,q,h in zip(supports,upper,lower) if 0 in B)
        assert gain==scale_v*(u[C]*exact_pair+tau*triple_budget)>0
        assert gain==scale_v*(u[C]*(exact_pair+triple_piece)+tau*negative_bracket)
        cases.append(dict(roots=sorted(A),source=x,other_pair=sorted(C),
                          alpha=str(alpha),harnack_bar=str(bar),kappa=str(kappa),
                          source_weight=str(alpha*ratio),source_weight_bound=str(lam),
                          exact_pair_residual=str(exact_pair),
                          weighted_triple_residual=str(triple_piece),
                          triple_budget_divided_by_tau=str(triple_budget),
                          occupation_gain=str(gain)))
    return dict(kernel=[[str(q) for q in row] for row in K],p=list(map(str,p)),
                common_odds=str(odds),first_contact_operator_entries=operator_entries,
                source_harnack_checks=harnack_checks,exact_support_decrease_checks=exact_support_checks,
                raw_transient_states=raw_states,raw_law_entries=len(starts)*size,
                h_1_hits_2=str(hit[2][1]),h_2_hits_1=str(hit[1][2]),lambda_bound=str(lam),
                uniform_odds_domain='r>0',
                u_numerators_over_65536_one_plus_r_squared={
                    ','.join(map(str,sorted(C))):list(map(str,poly)) for C,poly in U.items()},
                positive_core_coefficients_ascending={x:list(map(str,poly))
                                                      for x,poly in positive_cores.items()},
                negative_unweighted_source_bracket=str(negative_bracket),comparisons=cases,
                scope='two root laws on this fixed kernel for every common odds; not whole-cone closure')



def cycle_weighted_budget_obstruction():
    """A cycle obstruction, and two source-weighted repairs, in RN26.

    Standalone fragment for verify_pair_capacity.py; no main invocation.
    Standard library only, exact rational raw absorption and polynomial
    identities, with no parameter search or generated certificate files.

    On 0->1 with probability 1/4, 1->2 surely, 2->0 with probability
    1/2, uniform sleep p, target 0, the beta-weighted triangle fails
    throughout 0<p<=1/10 and its normalized slack tends to -1/64.
    The common positive factor r_0 S_00 is omitted from all budgets.
    This is not a counterexample to occupation monotonicity.

    For the same cycle, the two source-weighted bounds with root sets
    01 and 02 hold for every 0<p<1: even replacing both beta factors
    by one gives positive quadratic slacks, with cross weight 1/8.
    """
    from fractions import Fraction as F
    import verify_response_renewal as rn

    K = [[F(0), F(1, 4), F(0), F(3, 4)],
         [F(0), F(0), F(1), F(0)],
         [F(1, 2), F(0), F(0), F(1, 2)]]
    p0 = F(1, 10)
    sleeps = [p0] * 3
    V = frozenset(range(3))
    pairs = [V - {i} for i in range(3)]
    C01, C02 = frozenset((0, 1)), frozenset((0, 2))
    starts = {C: tuple(int(i in C) for i in range(3)) for C in pairs + [V]}
    # The two repeated-source inputs independently audit the IDLA factor.
    repeated = {C01: (1, 0, 2), C02: (1, 2, 0)}
    raw, state_count = rn.stopped(
        list(starts.values()) + list(repeated.values()),
        lambda s: rn.arw_step(K, sleeps, s),
        lambda s: all(a <= 0 for a in s))

    def mass(start, C):
        return raw[start].get(tuple(-int(i in C) for i in range(3)), F(0))

    full = mass(starts[V], V)
    assert full == F(8, 7001)
    entrance, beta = {}, {}
    for C in (C01, C02):
        # In the boundary order C,V, this is the raw Q_C of RN12/16.
        Q = [[mass(starts[E], B) for B in (C, V)] for E in (C, V)]
        mC, mV = map(sum, rn.inverse(Q))
        assert mC == F(56223271, 633520)
        assert 1 <= mV <= mC
        entrance[C] = (mC, mV)
        beta[C] = mV / mC
        assert beta[C] == (Q[0][0] - Q[1][0]) / full
    assert beta[C01] == F(7618127, 8649734)
    assert beta[C02] == F(50085190, 56223271)

    # Independent RN25 pair renewal, specialized to uniform-p cycles.
    jump = [row[:3] for row in K]
    jump2 = rn.multiply(jump, jump)
    q0, c = 1 - p0, F(1, 8)
    R = [[q0 * jump[i][j] + p0 * q0 * jump2[i][j]
          for j in range(3)] for i in range(3)]
    G = rn.inverse([[F(i == j) - R[i][j] for j in range(3)] for i in range(3)])
    z = [q0 / p0 * (1 - sum(R[i]) - p0**2 * c) for i in range(3)]
    Z = (1 - c) / p0**3 + c
    assert full == 1 / Z
    for i in range(3):
        assert mass(starts[V], pairs[i]) == sum(G[i][j] * z[j] for j in range(3)) / Z
        for j in range(3):
            assert mass(starts[pairs[j]], pairs[i]) == p0**2 * G[i][j]
    for i in (1, 2):
        assert beta[pairs[i]] == Z * p0**2 * G[i][i] - sum(G[i][j] * z[j] for j in range(3))

    # Raw killed forced-jump laws, without using the pair renewal.
    backgrounds = [frozenset(), frozenset((1,)), frozenset((2,)), frozenset((1, 2))]
    after_jump = {D: (0, 2 if 1 in D else 1, -int(2 in D)) for D in backgrounds}
    killed, killed_count = rn.stopped(
        after_jump.values(), lambda s: rn.arw_step(K, sleeps, s),
        lambda s: s[0] > 0 or all(a <= 0 for a in s))
    w = {D: F(3, 4) + F(1, 4) * sum(weight for end, weight in killed[start].items()
                                     if end[0] == 0)
         for D, start in after_jump.items()}
    assert [w[D] for D in backgrounds] == [F(719, 800), F(1321, 1600),
                                           F(1339, 1600), F(2521, 3200)]
    weights = {C: w[frozenset()] - w[C - {0}] for C in (C01, C02, V)}
    slack = weights[V] - sum(weights[C] * beta[C] for C in (C01, C02))
    assert slack == F(-772388791, 89957233600)

    # Small univariate identities: coefficients are in increasing powers of p.
    def trim(a):
        a = list(a)
        while len(a) > 1 and a[-1] == 0:
            a.pop()
        return tuple(a)

    def add(*polys):
        out = [F(0)] * max(map(len, polys))
        for a in polys:
            for i, x in enumerate(a):
                out[i] += x
        return trim(out)

    def scale(a, x):
        return trim(F(x) * y for y in a)

    def mul(*polys):
        out = (F(1),)
        for a in polys:
            new = [F(0)] * (len(out) + len(a) - 1)
            for i, x in enumerate(out):
                for j, y in enumerate(a):
                    new[i + j] += x * y
            out = trim(new)
        return out

    def evaluate(a, x):
        return sum(y * x**i for i, y in enumerate(a))

    one, p, q = (F(1),), (F(0), F(1)), (F(1), F(-1))
    D = add(one, scale(mul(q, q, q), -c), scale(mul(p, q, q), -3*c),
            scale(mul(p, p, p, q, q, q), -c*c))
    A = add(one, scale(mul(p, q, q), -c))
    U = add(q, scale(mul(p, p, q, q), c))
    # (I-qK-pqK^2)(A I+U K+q K^2)=D I, using K^3=c I.
    assert add(A, scale(mul(q, q), -c), scale(mul(p, q, U), -c)) == D
    assert add(U, scale(mul(q, A), -1), scale(mul(p, q, q), -c)) == (F(0),)
    assert add(q, scale(mul(q, U), -1), scale(mul(p, q, A), -1)) == (F(0),)
    assert D == tuple(F(x, 64) for x in (56, 0, 24, -17, 3, -3, 1))
    # This decomposition proves D>0 throughout [0,1].
    assert scale(D, 64) == add((F(56),), mul(p, p, add((F(7),), scale(q, 17))),
                              scale(mul(p, p, p, p, q), 3), mul(p, p, p, p, p, p))
    N = {}
    for i, ai, next_ai in ((2, F(1, 2), F(1, 4)), (1, F(1), F(1, 2))):
        row_sum = add(A, scale(U, ai), scale(q, ai * next_ai))
        numerator = add(mul(add((1-c,), scale(mul(p, p, p), c)), A),
                        scale(mul(q, D), -1), scale(mul(p, p, q, row_sum), c))
        assert numerator[0] == 0
        N[i] = trim(numerator[1:])  # Exact division by p; beta_i=N_i/D.
        assert evaluate(N[i], p0) / evaluate(D, p0) == beta[pairs[i]]
    assert N[2] == tuple(F(x, 128) for x in (98, 6, 46, -25, 7, -5, 1))
    assert N[1] == tuple(F(x, 64) for x in (49, 10, 9, -5, 2, -1))

    # Direct two-buffer queue formulas, valid for every uniform p.
    t = scale(add(one, p), F(1, 2))
    w00 = add((F(3, 4),), scale(add(p, mul(q, t)), F(1, 4)))
    w10 = add((F(3, 4),), scale(mul(t, t), F(1, 4)))
    w01 = add((F(3, 4),), scale(add(p, scale(mul(q, t), F(1, 2))), F(1, 4)))
    w11 = add((F(3, 4),), scale(mul(t, t), F(1, 8)))
    a01, a02, aV = [add(w00, scale(other, -1)) for other in (w10, w01, w11)]
    assert a01 == tuple(F(x, 16) for x in (1, 2, -3))
    assert a02 == tuple(F(x, 16) for x in (1, 0, -1))
    assert aV == tuple(F(x, 32) for x in (3, 6, -5))
    assert [evaluate(a, p0) for a in (a01, a02, aV)] == [weights[C] for C in (C01, C02, V)]
    numerator = add(mul(aV, D), scale(mul(a01, N[2]), -1), scale(mul(a02, N[1]), -1))
    assert numerator == tuple(F(x, 1024) for x in (-14, 57, 54, 37, -9, 1, 2, 1, -1))
    assert evaluate(numerator, p0) / evaluate(D, p0) == slack
    limit = numerator[0] / D[0]
    assert limit == F(-1, 64)
    assert N[1][0] / D[0] == N[2][0] / D[0] == F(7, 8)
    # Drop negative nonconstant terms and bound positive ones at p=1/10.
    # Since D>0, this proves failure on the entire interval 0<p<=1/10.
    interval_upper = numerator[0] + sum(max(x, F(0)) * p0**i
                                      for i, x in enumerate(numerator) if i)
    assert interval_upper == F(-77229879, 10240000000) < 0

    # The stronger source-weighted bounds need only the proved beta<=1.
    source01 = add(aV, scale(a01, -1), scale(a02, -c))
    source02 = add(aV, scale(a02, -1), scale(a01, -c))
    assert source01 == tuple(F(x, 128) for x in (3, 8, 5))
    assert source02 == tuple(F(x, 128) for x in (3, 22, -9))
    assert scale(source02, 128) == add((F(3),), scale(p, 13), scale(mul(p, q), 9))
    # Actual cross-source factors are c*U/A and c*q/A. Both are <=c:
    # A-U=p[1-c*q^2*(1+p)], A-q=p[1-c*q^2], and A>0.
    assert add(A, scale(U, -1)) == mul(p, add(one, scale(mul(q, q, add(one, p)), -c)))
    assert add(A, scale(q, -1)) == mul(p, add(one, scale(mul(q, q), -c)))
    walk = rn.inverse([[F(i == j) - jump[i][j] for j in range(3)] for i in range(3)])
    hit = lambda i, j: walk[i][j] / walk[j][j]
    assert hit(1, 2) == 1 and hit(2, 1) == c
    source_checks = []
    for roots, other, x, missing, gamma_poly, lower_poly in (
            (C01, C02, 2, 1, U, source01), (C02, C01, 1, 2, q, source02)):
        assert mass(repeated[roots], V) == hit(x, missing) * full
        gamma = (mass(starts[roots], other) / mass(starts[other], other)
                 * mass(repeated[roots], V) / full)
        assert gamma == c * evaluate(gamma_poly, p0) / evaluate(A, p0)
        assert 0 <= gamma <= c
        source_slack = weights[V] - weights[roots] * beta[roots] - gamma * weights[other] * beta[other]
        assert source_slack >= evaluate(lower_poly, p0) > 0
        source_checks.append(dict(roots=sorted(roots), extra_source=x,
                                  actual_cross_weight=str(gamma), exact_slack=str(source_slack)))

    return dict(kernel=[[str(x) for x in row] for row in K], uniform_p=str(p0),
                target=0, full_mass=str(full), raw_transient_states=state_count,
                killed_transient_states=killed_count, pair_renewal_entries=12,
                w_values=[str(w[D]) for D in backgrounds],
                beta_01=str(beta[C01]), beta_02=str(beta[C02]),
                normalized_weighted_triangle_slack=str(slack),
                slack_numerator_coefficients=[str(x) for x in numerator],
                common_denominator_coefficients=[str(x) for x in D],
                rare_sleep_limit=str(limit), failure_interval='0<p<=1/10',
                numerator_upper_on_interval=str(interval_upper),
                source_weighted_uniform_p_range='0<p<1',
                source_01_lower_bound='(3+8*p+5*p^2)/128',
                source_02_lower_bound='(3+22*p-9*p^2)/128',
                source_checks=source_checks,
                scope='weighted triangle fails; two specified actual-source budgets hold; no general N2 claim')


def source_maximum_principle_audit():
    """First-contact operator identity and ordinary-walk source bounds."""
    from verify_feedback_cone import binary_stabilization
    case=rn.CASES[4]
    K=[[F(q) for q in row] for row in case['K']]
    p=list(map(F,case['p'])); n=len(p)
    supports,P,E,D=binary_stabilization(K,p)
    pos={C:i for i,C in enumerate(supports)}; size=len(supports)
    A=[rn.multiply(rn.multiply(E,mat),P) for mat in D]
    plus=[rn.multiply(mat,P) for mat in D]
    _,_,_,W=rn.positive_operator(K,p,n)
    residual_columns=[]; operator_checks=0
    for y in range(n):
        columns=[C for C in supports if y not in C]
        r=(1-p[y])/p[y]
        L=[[F(B==C) if y not in B else
            -r*W[tuple(sorted(B)),y].get(tuple(sorted(C)),F(0))
            for C in columns] for B in supports]
        contact=[[F(B==C|{y})-A[y][pos[C]][pos[B]] for B in supports]
                 for C in columns]
        J=[[K[y][n]*F(i==j)+sum(K[y][z]*A[z][i][j] for z in range(n))
            for j in range(size)] for i in range(size)]
        assert rn.multiply(L,contact)==[
            [r*(A[y][i][j]-J[i][j]) for j in range(size)] for i in range(size)]
        residual_columns.append([
            sum(L[i][k]*A[y][pos[C]][pos[frozenset()]] for k,C in enumerate(columns))/r
            for i in range(size)])
        operator_checks+=size*size
    green=rn.inverse([[F(i==j)-K[i][j] for j in range(n)] for i in range(n)])
    for x in range(n):
        for C in supports:
            left=F(not C)-A[x][pos[C]][pos[frozenset()]]
            right=sum(green[x][y]*residual_columns[y][pos[C]] for y in range(n))
            assert left==right
    h={}
    for j in range(n):
        other=[i for i in range(n) if i!=j]
        inv=rn.inverse([[F(i==z)-K[i][z] for z in other] for i in other])
        values=rn.multiply(inv,[[K[i][j]] for i in other])
        h[j]={i:values[k][0] for k,i in enumerate(other)}; h[j][j]=F(1)
    harnack=exact_support=strict=0
    for roots in supports:
        for x in range(n):
            for C in supports:
                if C<=roots:
                    assert plus[x][pos[roots]][pos[C]]<=P[pos[roots]][pos[C]]
                    exact_support+=1
        for j in range(n):
            if j in roots:
                continue
            C=roots|{j}; diagonal=P[pos[C]][pos[C]]
            assert diagonal>0
            for i in range(n):
                ratio=plus[i][pos[roots]][pos[C]]/diagonal
                assert ratio<=h[j][i]
                if i not in C and h[j][i]>0:
                    assert ratio<h[j][i]
                    strict+=1
                harnack+=1
    # Ordinary target-capacity IDLA: parking is allowed only in C.
    from functools import lru_cache
    @lru_cache(None)
    def hit(missing,protected=None):
        targets=sorted(missing); other=[i for i in range(n) if i not in missing]
        jump=[[(1-p[i] if protected is not None and i not in protected else 1)*q
               for q in K[i]] for i in range(n)]
        inv=rn.inverse([[F(i==z)-jump[i][z] for z in other] for i in other])
        values=rn.multiply(inv,[[jump[i][z] for z in targets] for i in other]) if other else []
        rows={i:dict(zip(targets,values[k])) for k,i in enumerate(other)}
        rows.update({i:{z:F(i==z) for z in targets} for i in targets})
        return rows
    def fills(order,target,killed=False):
        @lru_cache(None)
        def step(k,parked):
            if k==len(order):
                return F(parked==target)
            missing=target-parked
            return sum(weight*step(k+1,parked|{j})
                       for j,weight in hit(missing,target if killed else None)[order[k]].items())
        return step(0,frozenset())
    idla_checks=duplicate_checks=0
    maximum_lifetime=max(map(sum,green))
    def sandwich(order,target,ratio):
        upper=fills(order,target); lower=fills(order,target,True)
        assert upper==fills(tuple(reversed(order)),target)
        assert lower==fills(tuple(reversed(order)),target,True)
        assert lower<=ratio<=upper
        epsilon=max((p[y] for y in range(n) if y not in target),default=F(0))
        assert upper-lower<=len(target)*epsilon*maximum_lifetime
    for roots in supports:
        for target in supports:
            if len(roots)!=len(target):
                continue
            sandwich(tuple(sorted(roots)),target,
                     P[pos[roots]][pos[target]]/P[pos[target]][pos[target]])
            idla_checks+=1
    for roots in supports:
        for x in roots:
            for target in supports:
                if len(target)!=len(roots)+1:
                    continue
                sandwich(tuple(sorted(roots))+(x,),target,
                         plus[x][pos[roots]][pos[target]]/P[pos[target]][pos[target]])
                duplicate_checks+=1
    return dict(first_contact_operator_entries=operator_checks,
                nonempty_green_entries=n*size,source_harnack_bounds=harnack,
                strict_replacement_bounds=strict,exact_support_decreases=exact_support,
                target_capacity_IDLA_sandwiches=idla_checks,
                duplicate_source_sandwiches=duplicate_checks)


def rare_sleep_expansion_audit(K=None):
    """Raw zero-sleep Taylor coefficients, independent of the IDLA formulas."""
    if K is None:
        case=rn.CASES[4]
        K=[[F(q) for q in row] for row in case['K']]
        n=len(K)
        K=[[q/2 for q in row[:n]]+[(1+row[n])/2] for row in K]
    n=len(K)
    supports=[frozenset(i for i,b in enumerate(bits) if b) for bits in product((0,1),repeat=n)]
    starts={A:tuple(int(i in A) for i in range(n)) for A in supports}
    states=[s for s in starts.values() if any(q>0 for q in s)]
    pos={s:i for i,s in enumerate(states)}
    for s in states:
        for t,q in rn.arw_step(K,[F(1,2)]*n,s).items():
            if q and any(a>0 for a in t) and t not in pos:
                pos[t]=len(states);states.append(t)
    zero=[F(0)]*n
    rows=[]; derivatives=[]; exits=[]; dexits=[]
    for s in states:
        base=rn.arw_step(K,zero,s)
        z=next(i for i,a in enumerate(s) if a>0)
        unit=zero[:];unit[z]=F(1)
        changed=rn.arw_step(K,unit,s)
        derivative={t:changed.get(t,F(0))-base.get(t,F(0)) for t in set(base)|set(changed)}
        rows.append({pos[s]:F(1)})
        rhs=[F(0)]*n; drhs=[F(0)]*n; drow={}
        for t,q in base.items():
            if t in pos:
                rows[-1][pos[t]]=rows[-1].get(pos[t],F(0))-q
            else:
                for v in range(n):rhs[v]+=q*F(t[v]==-1)
        for t,q in derivative.items():
            if t in pos:
                drow[pos[t]]=q
            else:
                for v in range(n):drhs[v]+=q*F(t[v]==-1)
        derivatives.append((z,drow));exits.append(rhs);dexits.append(drhs)
    q0=rn.old.solve_rows(rows,exits)
    rhs1=[[F(0)]*(n*n) for s in states]
    for i,(z,drow) in enumerate(derivatives):
        for v in range(n):
            rhs1[i][z*n+v]=dexits[i][v]+sum(weight*q0[j][v] for j,weight in drow.items())
    q1=rn.old.solve_rows(rows,rhs1)
    terms=list(combinations(range(n),2))+[(i,i) for i in range(n)]
    rhs2=[[F(0)]*(len(terms)*n) for s in states]
    for i,(z,drow) in enumerate(derivatives):
        for k,(a,b) in enumerate(terms):
            for v in range(n):
                if z==a:rhs2[i][k*n+v]+=sum(weight*q1[j][b*n+v] for j,weight in drow.items())
                if a!=b and z==b:rhs2[i][k*n+v]+=sum(weight*q1[j][a*n+v] for j,weight in drow.items())
    q2=rn.old.solve_rows(rows,rhs2)
    green=rn.inverse([[F(i==j)-K[i][j] for j in range(n)] for i in range(n)])
    h=lambda i,j:green[i][j]/green[j][j]
    assert all(0<h(i,j)<1 for i in range(n) for j in range(n) if i!=j)
    def coefficient(table,A,col):
        return table[pos[starts[A]]][col] if A else F(0)
    first=second=strict=0
    for A in supports:
        for x in range(n):
            if x in A:continue
            upper=A|{x}
            for v in range(n):
                value=coefficient(q1,upper,v*n+v)-coefficient(q1,A,v*n+v)
                expected=green[v][v]*h(x,v)
                for a in A:expected*=1-h(a,v)
                assert value==expected,(A,x,v,value,expected)
                first+=1
                if v not in A:
                    assert value>0;strict+=1
                    continue
                for c in range(n):
                    a,b=sorted((v,c));k=terms.index((a,b))
                    value=coefficient(q2,upper,k*n+v)-coefficient(q2,A,k*n+v)
                    expected=F(0)
                    if c not in A:
                        expected=green[v][c]*green[c][v]*(1-h(c,v))*h(x,c)
                        for y in A-{v}:expected*=1-h(y,c)
                        assert expected>0;strict+=1
                    assert value==expected,(A,x,v,c,value,expected)
                    second+=1
    # At zero outside sleep, trace onto a fixed two-site target. This
    # checks the whole-input IDLA identity, not only equal-count inputs.
    from functools import lru_cache
    target=frozenset({0,n-1}); target_sites=sorted(target)
    face_p=[F(i+1,n+1) if i in target else F(0) for i in range(n)]
    raw,face_states=rn.stopped(starts.values(),lambda s:rn.arw_step(K,face_p,s),
                              lambda s:all(a<=0 for a in s))
    end=tuple(-int(i in target) for i in range(n))
    diagonal=raw[starts[target]].get(end,F(0))
    trace_IminusK=rn.inverse([[green[i][j] for j in target_sites] for i in target_sites])
    det=(trace_IminusK[0][0]*trace_IminusK[1][1]
         -trace_IminusK[0][1]*trace_IminusK[1][0])
    odds=[(1-face_p[i])/face_p[i] for i in target_sites]
    Z=1+det*(odds[0]*odds[1]+odds[1]/trace_IminusK[0][0]
             +odds[0]/trace_IminusK[1][1])
    assert diagonal==1/Z
    @lru_cache(None)
    def hitting(missing):
        boundary=sorted(missing); other=[i for i in range(n) if i not in missing]
        inv=rn.inverse([[F(i==j)-K[i][j] for j in other] for i in other])
        values=rn.multiply(inv,[[K[i][j] for j in boundary] for i in other])
        rows={i:dict(zip(boundary,values[k])) for k,i in enumerate(other)}
        rows.update({i:{j:F(i==j) for j in boundary} for i in boundary})
        return rows
    def fill(order):
        @lru_cache(None)
        def step(k,parked):
            if parked==target:return F(1)
            if k==len(order):return F(0)
            hit=hitting(target-parked)[order[k]]
            return ((1-sum(hit.values()))*step(k+1,parked)
                    +sum(weight*step(k+1,parked|{j}) for j,weight in hit.items()))
        return step(0,frozenset())
    for A in supports:
        assert raw[starts[A]].get(end,F(0))/diagonal==fill(tuple(sorted(A)))
        assert fill(tuple(sorted(A)))==fill(tuple(sorted(A,reverse=True)))
        if target<=A:
            assert raw[starts[A]]==raw[starts[target]]
    return dict(raw_states=len(states),first_derivatives=first,
                second_derivatives=second,strict_coefficients=strict,
                zero_outside_trace_states=face_states,trace_IDLA_checks=len(supports),
                trace_diagonal=str(diagonal))


def source_switch_target_audit():
    from fractions import Fraction as F
    from itertools import combinations
    import verify_response_renewal as rn
    from verify_feedback_cone import binary_stabilization

    K = [[F(0), F(1, 2), F(0), F(0), F(1, 2)],
         [F(0), F(0), F(1, 3), F(1, 3), F(1, 3)],
         [F(1, 3), F(0), F(0), F(1, 3), F(1, 3)],
         [F(1, 2), F(0), F(0), F(0), F(1, 2)]]
    p, v = (F(2, 3), F(1, 2), F(9, 10), F(1, 16)), 1
    n = len(p)
    supports, P, E, D = binary_stabilization(K, p)
    m, pos = len(supports), {C: i for i, C in enumerate(supports)}
    pairs = [frozenset(C) for C in combinations(range(n), 2)]
    mu = {R: P[pos[R]] for R in pairs}
    f = [F(v in C) for C in supports]
    dot = lambda a, b: sum((x*y for x, y in zip(a, b)), F(0))
    A = [rn.multiply(rn.multiply(E, D[x]), P) for x in range(n)]
    J = [[[K[x][n]*F(i == j)+sum(K[x][y]*A[y][i][j]
             for y in range(n)) for j in range(m)] for i in range(m)]
         for x in range(n)]
    B = [[[(A[x][i][j]-(1-p[x])*J[x][i][j])/p[x]
           for j in range(m)] for i in range(m)] for x in range(n)]
    Af = [[dot(row, f) for row in M] for M in A]
    g = [[dot(B[x][i], f)-f[i] for i in range(m)] for x in range(n)]
    states = [(R, x) for R in pairs for x in range(n)]
    response = {(R, x): dot(mu[R], Af[x])-dot(mu[R], f)
                for R, x in states}
    boundary = [(R, x) for R, x in states
                if x not in R and x != v and p[x] < F(1, 2)]
    interior = [s for s in states if s not in boundary]
    ipos, bpos = ({s: i for i, s in enumerate(interior)},
                  {s: i for i, s in enumerate(boundary)})
    Q = [[F(0)]*len(interior) for _ in interior]
    H = [[F(0)]*len(boundary) for _ in interior]
    reward = [F(0)]*len(interior)
    _, _, _, W = rn.positive_operator(K, p, 2)
    common = switches = outside = target_rows = no_visit = 0
    for row, (R, x) in enumerate(interior):
        edges = {}
        def add(state, weight):
            edges[state] = edges.get(state, F(0))+weight
        if x in R:
            for y in range(n):
                add((R, y), K[x][y])
            common += 1
        elif x == v:
            r = (1-p[v])/p[v]
            for C in supports:
                if v in C or len(C) > 2:
                    continue
                N = mu[R][pos[C]]-r*sum(
                    mu[R][j]*W[tuple(sorted(T)), v].get(tuple(sorted(C)), F(0))
                    for j, T in enumerate(supports) if v in T and len(T) <= 2)
                assert N >= 0
                no_visit += 1
                reward[row] += N*(1-Af[v][pos[C]])/r
            for y in range(n):
                add((R, y), K[v][y])
            target_rows += 1
        else:
            C = frozenset({v, x})
            assert p[x] >= F(1, 2)
            assert all(min(b) >= 0 and sum(b) == 1 for b in B[x])
            alpha = mu[R][pos[C]]/mu[C][pos[C]]
            sigma = [a-alpha*b for a, b in zip(mu[R], mu[C])]
            assert 0 <= alpha <= 1 and min(sigma) >= 0
            assert sigma[pos[C]] == 0 and sum(sigma) == 1-alpha
            assert all(g[x][j] >= 0 for j, T in enumerate(supports)
                       if len(T) <= 2 and T != C)
            h = dot(sigma, g[x])
            assert 0 <= h <= 1-alpha
            reward[row] = p[x]*h
            for y in range(n):
                add((R, y), (1-p[x])*K[x][y])
            add((C, x), p[x]*alpha)
            switches += 1
            outside += int(v not in R)
        assert reward[row] >= 0 and all(w >= 0 for w in edges.values())
        assert sum(edges.values()) <= 1
        assert response[R, x] == reward[row]+sum(
            w*response[s] for s, w in edges.items())
        for s, w in edges.items():
            if s in ipos:
                Q[row][ipos[s]] += w
            else:
                H[row][bpos[s]] += w

    G = rn.inverse([[F(i == j)-Q[i][j] for j in range(len(interior))]
                    for i in range(len(interior))])
    assert all(w >= 0 for row in G for w in row)
    weights = rn.multiply(G, H)
    paid = [dot(row, reward) for row in G]
    for i, s in enumerate(interior):
        assert min(weights[i]) >= 0 and sum(weights[i]) <= 1
        assert response[s] == paid[i]+sum(
            w*response[b] for w, b in zip(weights[i], boundary))
    assert (common, switches, outside, target_rows) == (12, 6, 2, 3)
    assert (len(interior), len(boundary), no_visit) == (21, 3, 21)
    return dict(root_pairs=len(pairs), common_root_equations=common,
                high_source_switches=switches, target_outside_switches=outside,
                positive_target_reward_equations=target_rows,
                no_visit_residual_coordinates=no_visit,
                interior_states=len(interior), low_extra_boundary_states=len(boundary),
                exact_resolvent_rows=len(interior))


def cycle_actual_source_budget_obstruction():
    """Exact RN26 (127i) failure on a single-leak three-cycle.

    Integration fragment; standard library plus the maintained raw solver.
    No main invocation, parameter search, or generated certificate files.
    Vertices (v,a,b)=(0,1,2): v->a/s equally, a->b, b->v.
    Sleeps (1/2,t,1/2), original A={v,b}, extra a, other C={v,a}.
    The normalized source budget fails for 0<t<=1/50; the exact pair
    remainder restores a strictly positive occupation response.
    """
    from fractions import Fraction as F
    from itertools import permutations
    import verify_response_renewal as rn

    # Tiny univariate polynomial arithmetic, coefficients in ascending order.
    def trim(a):
        a = list(map(F, a))
        while len(a) > 1 and not a[-1]:
            a.pop()
        return tuple(a)

    def add(*args):
        out = [F(0)] * max(map(len, args))
        for a in args:
            for i, x in enumerate(a):
                out[i] += x
        return trim(out)

    def scale(a, x):
        return trim(F(x)*y for y in a)

    def mul(*args):
        out = (F(1),)
        for a in args:
            new = [F(0)] * (len(out)+len(a)-1)
            for i, x in enumerate(out):
                for j, y in enumerate(a):
                    new[i+j] += x*y
            out = trim(new)
        return out

    def ev(a, x):
        return sum(y*x**i for i, y in enumerate(a))

    zero, one, t, q = (F(0),), (F(1),), (F(0), F(1)), (F(1), F(-1))
    D = (F(40), F(15), F(1))
    pp, qq, edge = [(F(1,2),), t, (F(1,2),)], [(F(1,2),), q, (F(1,2),)], [F(1,2), F(1), F(1)]
    # H=diag(p) Rcal diag(p)^(-1). Its entries follow RN25 (115).
    H = [[zero for j in range(3)] for i in range(3)]
    for i in range(3):
        j, k = (i+1)%3, (i+2)%3
        H[i][j] = scale(qq[i], edge[i])
        H[i][k] = scale(mul(qq[i], pp[j]), edge[i]*edge[j])
    M = [[add(one if i == j else zero, scale(H[i][j], -1))
          for j in range(3)] for i in range(3)]
    det = add(*(scale(mul(*(M[i][perm[i]] for i in range(3))),
                     (-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3)))
                for perm in permutations(range(3))))
    adj = [[zero for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            rows, cols = [k for k in range(3) if k != j], [k for k in range(3) if k != i]
            adj[i][j] = scale(add(mul(M[rows[0]][cols[0]], M[rows[1]][cols[1]]),
                                  scale(mul(M[rows[0]][cols[1]], M[rows[1]][cols[0]]), -1)),
                              (-1)**(i+j))
    assert det == scale(D, F(1,64))
    assert all(c > 0 for c in D)  # Invertibility for t>=0.
    for i in range(3):
        for j in range(3):
            assert add(*(mul(M[i][k], adj[k][j]) for k in range(3))) == (det if i == j else zero)
    # zhat=diag(p)z. Pair pi(R_i)=tau*(adj zhat)_i/(p_i det).
    zhat = [(F(1,4),), scale(q, F(1,8)), (F(1,8),)]
    Az = [add(*(mul(adj[i][j], zhat[j]) for j in range(3))) for i in range(3)]
    assert Az[1] == scale(q, F(1,2))
    assert Az[2] == (F(9,32), F(-1,32))
    # Z_V*p_0*p_1*p_2=(4+t)/8; beta_i=(Z_V*prod(p)*adj_ii-Az_i)/(p_i det).
    zp = (F(1,2), F(1,8))
    beta_num = {1: (F(36), F(-1)), 2: (F(20), F(26), F(2))}
    for i in (1,2):
        numerator = add(mul(zp, adj[i][i]), scale(Az[i], -1))
        assert mul(numerator, D) == mul(beta_num[i], pp[i], det)
    # alpha=mu_A(C)/mu_C(C)=adj[2,1]/adj[2,2]=2/(7+t).
    assert mul(adj[2][1], (F(7), F(1))) == scale(adj[2][2], 2)

    # Direct two-buffer no-return calculation (target jump has sink mass 1/2).
    w0, wA, wC, wV = (F(3,4), F(1,4)), (F(1,2), F(1,2)), (F(1,2), F(1,4)), (F(1,2),)
    uA, uC, uV = [add(w0, scale(w, -1)) for w in (wA,wC,wV)]
    assert (uA,uC,uV) == ((F(1,4), F(-1,4)), (F(1,4),), (F(1,4), F(1,4)))
    numerator = (F(-12), F(596), F(193), F(22), F(1))
    denominator = scale(mul((F(7), F(1)), D), 4)
    assert add(mul(uV, denominator),
               scale(mul(uA, beta_num[1], (F(7), F(1))), -4),
               scale(mul(uC, beta_num[2]), -8)) == numerator
    assert numerator[0]/denominator[0] == F(-3,280)
    assert all(i*c > 0 for i,c in enumerate(numerator) if i)
    assert ev(numerator,F(1,50)) == F(-16399,6250000) < 0
    # Positive exact-pair repair: rho=(1-alpha)pi(C), tau=2t/(4+t).
    repair_num = mul((F(9), F(-1)), (F(5), F(1)))
    positive = (F(24), F(84), F(15), F(1))
    assert add(numerator, scale(repair_num, 4)) == mul((F(7), F(1)), positive)
    assert all(c > 0 for c in positive)

    # Independent raw physical absorption at one strictly positive parameter.
    t0 = F(1,100)
    K = [[F(0),F(1,2),F(0),F(1,2)],
         [F(0),F(0),F(1),F(0)], [F(1),F(0),F(0),F(0)]]
    p = [F(1,2),t0,F(1,2)]
    V, A, C = frozenset(range(3)), frozenset((0,2)), frozenset((0,1))
    pairs = [V-{i} for i in range(3)]
    starts = {E: tuple(int(i in E) for i in range(3)) for E in pairs+[V]}
    duplicate = (1,2,0)
    raw, raw_count = rn.stopped(list(starts.values())+[duplicate],
                               lambda s: rn.arw_step(K,p,s), lambda s: all(a<=0 for a in s))
    assert raw[duplicate] == raw[starts[V]]
    def mass(start, E):
        return raw[start].get(tuple(-int(i in E) for i in range(3)), F(0))
    tau = mass(starts[V],V)
    assert tau == F(2,401) == 2*t0/(4+t0)
    beta, entrance = {}, {}
    for E in (A,C):
        boundary = (E,V)
        Q = [[mass(starts[J],L) for L in boundary] for J in boundary]
        mE, mV = map(sum,rn.inverse(Q))
        entrance[E], beta[E] = (mE,mV), mV/mE
        assert beta[E] == (Q[0][0]-Q[1][0])/tau
    assert beta[A] == F(359900,401501) == ev(beta_num[1],t0)/ev(D,t0)
    assert beta[C] == F(202602,401501) == ev(beta_num[2],t0)/ev(D,t0)
    assert entrance[A][1] == F(3599,1598)
    assert entrance[C][1] == F(101301,1402)

    # Original RN25 coefficients, independently evaluated without H/adj.
    r = [(1-x)/x for x in p]
    Z = [1+r[j]+r[k]+r[j]*r[k] for i in range(3) for j,k in [sorted(V-{i})]]
    R, z = [], []
    for i in range(3):
        row = [F(0)]*3
        for j in V-{i}:
            k = next(iter(V-{i,j}))
            row[j] = r[i]/Z[i]*(K[i][j]*(1+r[k])+K[i][k]*K[k][j])
        R.append(row)
        j,k = sorted(V-{i})
        z.append(r[i]/Z[i]*(F(1,2)+(Z[i]-1)*K[i][3]
                 +r[k]*K[i][j]*K[j][3]+r[j]*K[i][k]*K[k][3]))
    G = rn.inverse([[F(i==j)-R[i][j] for j in range(3)] for i in range(3)])
    for i in range(3):
        assert mass(starts[V],pairs[i]) == tau*sum(G[i][j]*z[j] for j in range(3))
        for j in range(3):
            assert mass(starts[pairs[j]],pairs[i]) == G[i][j]/Z[j]
            assert p[i]*R[i][j]/p[j] == ev(H[i][j],t0)
            assert p[i]*G[i][j]/p[j] == ev(adj[i][j],t0)/ev(det,t0)
        assert p[i]*z[i] == ev(zhat[i],t0)

    backgrounds = [frozenset(), frozenset((2,)), frozenset((1,)), frozenset((1,2))]
    after = {E: (0,2 if 1 in E else 1,-int(2 in E)) for E in backgrounds}
    killed, killed_count = rn.stopped(after.values(), lambda s: rn.arw_step(K,p,s),
                                     lambda s: s[0]>0 or all(a<=0 for a in s))
    w = {E: F(1,2)+F(1,2)*sum(prob for end,prob in killed[start].items() if end[0]==0)
         for E,start in after.items()}
    assert [w[E] for E in backgrounds] == [ev(a,t0) for a in (w0,wA,wC,wV)]
    u = {E: w[frozenset()]-w[E-{0}] for E in (A,C,V)}
    alpha = mass(starts[A],C)/mass(starts[C],C)
    kappa = mass(duplicate,V)/tau
    gamma = alpha*kappa
    assert gamma == alpha == F(200,701) == 2/(t0+7)
    assert kappa == 1
    slack = u[V]-u[A]*beta[A]-gamma*u[C]*beta[C]
    assert slack == ev(numerator,t0)/ev(denominator,t0) == F(-602067799,112580880400)
    assert tau*slack == F(-602067799,22572466520200)
    rho = mass(starts[V],C)-alpha*mass(duplicate,C)
    assert rho == (1-alpha)*mass(starts[V],C) == F(360319200,112862332601) > 0
    credit = u[C]*rho/tau
    assert credit == ev(repair_num,t0)/((t0+7)*ev(D,t0))
    repaired = slack+credit
    assert repaired == ev(positive,t0)/(4*ev(D,t0)) == F(24841501,160600400) > 0
    singleton_sigma = mass(starts[A],frozenset({0}))-alpha*mass(starts[C],frozenset({0}))
    assert singleton_sigma == F(80000,491401)
    direct_sleep_credit = u[C]*t0*singleton_sigma/tau
    assert direct_sleep_credit == F(40100,491401)
    assert slack+direct_sleep_credit == F(6018026512901,78919197160400) > 0
    assert t0*singleton_sigma <= rho
    Sv = p[0]*rn.inverse([[F(i==j)-(1-p[i])*K[i][j] for j in range(3)] for i in range(3)])[0][0]
    assert r[0]*Sv == 4/(t0+7)
    gain = sum(prob for end,prob in raw[starts[V]].items() if end[0]==-1)-sum(
        prob for end,prob in raw[starts[A]].items() if end[0]==-1)
    assert gain == r[0]*Sv*tau*repaired == F(49683002,112862332601)
    return dict(sleep_t=str(t0), roots=sorted(A), extra_source=1,
                raw_transient_states=raw_count, killed_transient_states=killed_count,
                pair_renewal_entries=12, beta_A=str(beta[A]), beta_C=str(beta[C]),
                actual_source_ratio=str(gamma), full_mass=str(tau),
                normalized_slack=str(slack), slack_without_rvSvv=str(tau*slack),
                slack_numerator=[str(x) for x in numerator],
                failure_interval='0<t<=1/50', boundary_limit='-3/280',
                positive_repair_numerator=[str(x) for x in positive],
                repaired_normalized_slack=str(repaired), occupation_gain=str(gain),
                direct_sleep_credit=str(direct_sleep_credit),
                scope='Actual (127i) fails; retained exact-pair remainder repairs this comparison.')


def low_source_feeder_budget_audit():
    """One four-site graph, two raw priorities, and a univariate certificate."""
    from collections import defaultdict
    from fractions import Fraction as F
    from math import comb
    import verify_response_renewal as rn
    from verify_feedback_cone import binary_stabilization

    # Vertex order v,a,b,z; the final column is the sink.  No arrow enters z.
    K = [[F(0), F(1, 2), F(0), F(0), F(1, 2)],
         [F(0), F(0), F(1), F(0), F(0)],
         [F(1, 2), F(0), F(0), F(0), F(1, 2)],
         [F(0), F(1), F(0), F(0), F(0)]]
    t = F(1, 16)
    p = [F(1, 2), t, F(1, 2), F(1, 2)]
    v, a, z = 0, 1, 3
    A, R, C, B = map(frozenset, ((0, 3), (0, 1), (0, 2), (0, 1, 2)))
    single = frozenset((v,))
    supports, P, E, D = binary_stabilization(K, p)
    pos = {S: i for i, S in enumerate(supports)}
    plus = rn.multiply(D[a], P)
    active_add = rn.multiply(E, plus)
    assert all(min(row) >= 0 and sum(row) == 1 for row in active_add)
    pairs = [S for S in supports if v in S and len(S) == 2]
    triples = [S for S in supports if v in S and len(S) == 3]
    boundary = {S: [T for T in supports if S <= T and len(T) <= 3]
                for S in pairs}
    initial = {S: tuple(int(i in S) for i in range(4))
               for S in set().union(*map(set, boundary.values()), {single})}
    extra = {S: tuple(int(i in S) + int(i == a) for i in range(4))
             for S in pairs}
    terminal = {S: tuple(-int(i in S) for i in range(4)) for S in supports}

    # Independent transition implementation, highest active vertex first.
    def reverse_step(state):
        x = max(i for i, count in enumerate(state) if count > 0)
        out = defaultdict(F)
        asleep = list(state)
        if asleep[x] == 1:
            asleep[x] = -1
        out[tuple(asleep)] += p[x]
        for y, probability in enumerate(K[x]):
            if not probability:
                continue
            end = list(state)
            end[x] -= 1
            if y < 4:
                end[y] = 2 if end[y] == -1 else end[y] + 1
            out[tuple(end)] += (1-p[x])*probability
        assert sum(out.values()) == 1
        return dict(out)

    starts = set(initial.values()) | set(extra.values())
    stop = lambda state: all(count <= 0 for count in state)
    ordinary_step = lambda state: rn.arw_step(K, p, state)
    raw, ordinary_count = rn.stopped(starts, ordinary_step, stop)
    reversed_raw, reverse_count = rn.stopped(starts, reverse_step, stop)
    for start in starts:
        assert raw[start] == reversed_raw[start]
    for S, start in initial.items():
        assert [raw[start].get(terminal[T], F(0)) for T in supports] == P[pos[S]]
    for S, start in extra.items():
        assert [raw[start].get(terminal[T], F(0)) for T in supports] == plus[pos[S]]

    def occupation(row):
        return sum(row[pos[S]] for S in supports if v in S)

    # Audit both feeder branches as full sublaws, not only one triple mass.
    pi, lower, upper = P[pos[B]], P[pos[A]], plus[pos[A]]
    for S in supports:
        if z not in S:
            assert lower[pos[S]] == P[pos[R]][pos[S]]/2
            assert upper[pos[S]] == pi[pos[S]]/2
            assert lower[pos[S | {z}]] == P[pos[single]][pos[S]]/2
            assert upper[pos[S | {z}]] == P[pos[R]][pos[S]]/2
    assert plus[pos[R]] == plus[pos[C]] == pi
    full_gain = occupation(upper)-occupation(lower)
    assert full_gain == (occupation(pi)-occupation(P[pos[single]]))/2
    assert full_gain == F(7922477, 499968237) > 0

    # Full four-site entrance matrices include the triples containing z.
    # Their B-row agrees with the closed three-site (S,B) calculation.
    mass, beta = {}, {}
    for S in pairs:
        Q = [[raw[initial[U]].get(terminal[T], F(0)) for T in boundary[S]]
             for U in boundary[S]]
        mass[S] = dict(zip(boundary[S], map(sum, rn.inverse(Q))))
        beta[S] = {T: mass[S][T]*P[pos[S]][pos[S]] for T in boundary[S]}
        assert all(F(0) < x <= 1 for x in beta[S].values())
        if S <= B:
            small = [[P[pos[U]][pos[T]] for T in (S, B)] for U in (S, B)]
            assert list(map(sum, rn.inverse(small)))[1] == mass[S][B]
    assert mass[R][B] == F(81717, 3856)
    assert mass[C][B] == F(201, 68)

    # Forced v jump, killed upon any return to v.  All seven backgrounds
    # needed for the complete four-site decomposition are audited.
    w_backgrounds = [single] + pairs + triples
    first = {}
    for S in w_backgrounds:
        law = defaultdict(F)
        for y, probability in enumerate(K[v]):
            if not probability:
                continue
            state = [-int(i in S and i != v) for i in range(4)]
            if y < 4:
                state[y] = 2 if state[y] == -1 else 1
            law[tuple(state)] += probability
        first[S] = law
    killed_starts = set().union(*(set(law) for law in first.values()))
    killed_stop = lambda state: state[v] > 0 or stop(state)
    killed, killed_count = rn.stopped(killed_starts, ordinary_step, killed_stop)
    killed_reverse, killed_reverse_count = rn.stopped(
        killed_starts, reverse_step, killed_stop)
    for start in killed_starts:
        # A hit stops immediately; the unfinished configuration at that hit
        # depends on priority.  Only the surviving stable sublaw is invariant.
        assert {end: q for end, q in killed[start].items() if end[v] == 0 and q} == {
            end: q for end, q in killed_reverse[start].items() if end[v] == 0 and q}
    w = {S: sum(probability*sum(q for end, q in killed[start].items()
                               if end[v] == 0)
                for start, probability in law.items()) for S, law in first.items()}
    assert [w[S] for S in (single, R, C, B)] == [
        (7+t)/8, (11+3*t)/16, (11+5*t)/16, (19+3*t)/32]
    for S in w_backgrounds:
        if z in S:
            assert w[S] == w[S-{z}]
    weight = {S: w[single]-w[S] for S in pairs+triples}
    assert [weight[S] for S in (R, C, B)] == [(3-t)/16, 3*(1-t)/16, (9+t)/32]
    common = 8/(15+t)  # r_v*S_vv; here r_v=1.
    assert common == P[pos[single]][pos[single]] == F(128, 241)

    q = t*t+31*t+208
    rho = P[pos[R]][pos[C]]/P[pos[C]][pos[C]]
    assert rho == (1-t)*(t+16)/(16-t) == F(257, 272)
    assert P[pos[R]][pos[R]] == 8*t*(t+15)/q
    assert P[pos[R]][pos[C]] == 4*(1-t)*(t+16)/q
    assert P[pos[C]][pos[C]] == 4*(16-t)/q
    assert pi[pos[B]] == 4*t/(12+t) == F(4, 193)
    assert pi[pos[R]] == 2*t*(17-t)*(t+24)/((12+t)*q)
    assert pi[pos[C]] == 32*(1-t)*(t+24)/((12+t)*q)
    assert beta[R][B] == (5*t*t+115*t+312)/(2*q) == F(27239, 35830)
    assert beta[C][B] == (7*t+188)/q == F(3216, 3583)
    budget = {}
    for T in triples:
        budget[T] = weight[T]*upper[pos[T]]-sum(
            weight[S]*lower[pos[S]]*mass[S][T]*plus[pos[S]][pos[T]]
            for S in pairs if S <= T)
    normalized = budget[B]/upper[pos[B]]
    assert normalized == weight[B]-weight[R]*beta[R][B]-weight[C]*rho*beta[C][B]
    assert normalized == F(-53071, 9745760)
    assert budget[B] == F(-53071, 940465840)
    assert common*budget[B] == F(-424568, 14165766715) < 0

    # Keep the discarded exact-pair remainders.  Check both the source-aware
    # decomposition and (127f) for the complete four-site target response.
    pair_remainder, epsilon = {}, {}
    for S in pairs:
        alpha = lower[pos[S]]/P[pos[S]][pos[S]]
        sigma = [u-alpha*c for u, c in zip(lower, P[pos[S]])]
        assert min(sigma) >= 0 and sigma[pos[S]] == 0
        after = rn.multiply([sigma], active_add)[0]
        assert after == [u-alpha*c for u, c in zip(upper, plus[pos[S]])]
        assert min(after) >= 0
        pair_remainder[S] = after[pos[S]]
        epsilon[S] = after[pos[S]]+sum(
            beta[S][T]*after[pos[T]] for T in triples if S <= T)
    assert pair_remainder[R] == pair_remainder[A] == 0
    assert pair_remainder[C] == (1-rho)*pi[pos[C]]/2
    restored = weight[C]*pair_remainder[C]/upper[pos[B]]
    corrected = normalized+restored
    assert restored == F(259875, 1949152)
    assert corrected == (2*t*t+43*t+51)/(2*q) == F(2291, 17915) > 0
    core_gain = occupation(pi)-occupation(P[pos[R]])
    assert core_gain == common*pi[pos[B]]*corrected == F(1172992, 833280395)
    assert full_gain == common*(sum(budget.values())+sum(
        weight[S]*pair_remainder[S] for S in pairs))
    assert full_gain == common*(sum(weight[S]*epsilon[S] for S in pairs)+sum(
        upper[pos[T]]*(weight[T]-sum(weight[S]*beta[S][T]
                                     for S in pairs if S <= T)) for T in triples))

    # A smaller, generally available repair: tau=sigma_C A_a has tau L_a=0,
    # hence tau(C)>=r_a W_a(B,C)tau(B).  Different bad triples v,a,b use
    # different exact pairs v,b, so these credits are never reused.
    # Here W_a(B,C)=5/16, independent of t since a is a killing boundary.
    a_start = (-1, 0, 2, 0)
    a_stop = lambda state: state[a] > 0 or stop(state)
    a_killed, _ = rn.stopped([a_start], ordinary_step, a_stop)
    a_killed_reverse, _ = rn.stopped([a_start], reverse_step, a_stop)
    assert {end: x for end, x in a_killed[a_start].items() if end[a] == 0 and x} == {
        end: x for end, x in a_killed_reverse[a_start].items() if end[a] == 0 and x}
    w_a_self = a_killed[a_start][terminal[C]]
    assert w_a_self == F(5, 16)
    tau_B = upper[pos[B]]-(lower[pos[C]]/P[pos[C]][pos[C]])*plus[pos[C]][pos[B]]
    self_pair_credit = (1-t)/t*w_a_self*tau_B
    assert pair_remainder[C] >= self_pair_credit >= 0
    local_repair = weight[C]*self_pair_credit/upper[pos[B]]
    assert local_repair == F(50625, 1114112)
    assert normalized+local_repair == F(798257467, 19959316480) > 0

    # Integer polynomial arithmetic, coefficients in ascending order.
    # This checks the entire t-family; no interpolation or numerical screen.
    def trim(a):
        a = list(a)
        while len(a) > 1 and a[-1] == 0:
            a.pop()
        return tuple(a)

    def add(*polys):
        out = [0]*max(map(len, polys))
        for poly in polys:
            for i, x in enumerate(poly):
                out[i] += x
        return trim(out)

    def scale(poly, c):
        return trim(c*x for x in poly)

    def mul(*polys):
        out = (1,)
        for poly in polys:
            new = [0]*(len(out)+len(poly)-1)
            for i, x in enumerate(out):
                for j, y in enumerate(poly):
                    new[i+j] += x*y
            out = trim(new)
        return out

    def det2(a, b, c, d):
        return add(mul(a, d), scale(mul(b, c), -1))

    # M=diag(4,8t,8)(I-H), with H the compact three-site pair renewal.
    M = [[(4,), (0, -2), (0, -1)],
         [(-1, 1), (0, 8), (-4, 4)],
         [(-2,), (0, -1), (8,)]]
    adj = [[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rows, cols = [k for k in range(3) if k != j], [k for k in range(3) if k != i]
            adj[i][j] = scale(det2(M[rows[0]][cols[0]], M[rows[0]][cols[1]],
                                    M[rows[1]][cols[0]], M[rows[1]][cols[1]]), (-1)**(i+j))
    Q, var, one_minus, sixteen_minus = (208, 31, 1), (0, 1), (1, -1), (16, -1)
    for i in range(3):
        for j in range(3):
            assert add(*(mul(M[i][k], adj[k][j]) for k in range(3))) == (
                mul(var, Q) if i == j else (0,))
    assert scale(adj[2][2], 4) == scale(mul(var, (15, 1)), 8)
    assert scale(adj[1][2], 4) == scale(mul(one_minus, (16, 1)), 4)
    assert scale(adj[1][1], 2) == scale(sixteen_minus, 4)
    rhs = [(4, 1), (5, -5), (10,)]  # Twice diag(4,8t,8) times renewal source.
    full_pair = [add(*(mul(adj[i][j], rhs[j]) for j in range(3))) for i in range(3)]
    assert full_pair[2] == mul(var, (17, -1), (24, 1))
    assert full_pair[1] == scale(mul(one_minus, (24, 1)), 16)
    beta_R_num, beta_C_num = (312, 115, 5), (188, 7)
    assert add(scale(mul((12, 1), (15, 1)), 4),
               scale(mul((17, -1), (24, 1)), -1)) == beta_R_num
    assert add(mul((12, 1), sixteen_minus),
               scale(mul(one_minus, (24, 1)), -8)) == mul(var, beta_C_num)
    N = (192, -2539, 794, 110, 3)
    raw_numerator = add(mul((9, 1), sixteen_minus, Q),
                        scale(mul((3, -1), beta_R_num, sixteen_minus), -1),
                        scale(mul(one_minus, one_minus, (16, 1), beta_C_num), -6))
    assert raw_numerator == scale(N, -16)
    # E(t)=-N/[2(16-t)Q].  Q>0, 16-t>0 on 0<t<=1/2.
    # On 0<t<=1/16: N>=192-2539/16=533/16>0.
    assert N[0]+N[1]*F(1, 16) == F(533, 16) > 0
    assert min(N[2:]) > 0
    # The retained pair remainder cancels the problematic denominator:
    # -N+3(1-t)^2(t+14)(t+24)=(16-t)(2t^2+43t+51).
    assert add(scale(N, -1), scale(mul(one_minus, one_minus, (14, 1), (24, 1)), 3)) == mul(
        sixteen_minus, (51, 43, 2))

    # The smaller self-triple credit also suffices throughout 0<t<=1/2.
    # Its numerator over 256(16-t)Q has positive Bernstein coefficients
    # after t=u/2, 0<=u<=1.  Verify the conversion, not sampled values.
    local_num = add(scale(N, -128), scale(mul(one_minus, one_minus, (14, 1), Q), 15))
    assert local_num == (19104, 247262, -76537, -5785, 261, 15)
    bernstein = [sum(local_num[j]*F(1, 2)**j*F(comb(k, j), comb(5, j))
                     for j in range(k+1)) for k in range(6)]
    assert bernstein == [F(19104), F(219151, 5), F(2665719, 40),
                         F(6997601, 80), F(8499381, 80), F(3932621, 32)]
    back = [sum(bernstein[k]*comb(5, k)*comb(5-k, j-k)*(-1)**(j-k)
                for k in range(j+1)) for j in range(6)]
    assert back == [x*F(1, 2)**j for j, x in enumerate(local_num)]
    assert min(bernstein) > 0

    return dict(
        kernel=[[str(x) for x in row] for row in K], sleep=[str(x) for x in p],
        roots=sorted(A), extra=a, target=v, bad_triple=sorted(B),
        raw_transient_states=ordinary_count, reversed_transient_states=reverse_count,
        killed_transient_states=killed_count, killed_reversed_states=killed_reverse_count,
        four_site_budget_slack=str(common*budget[B]),
        budget_slack_without_common_factor=str(budget[B]),
        normalized_budget_slack=str(normalized), full_target_gain=str(full_gain),
        normalized_pair_repair=str(restored), normalized_repaired_budget=str(corrected),
        self_triple_pair_credit=str(local_repair),
        normalized_self_triple_repaired_budget=str(normalized+local_repair),
        self_triple_repair_Bernstein_on_0_half=[str(x) for x in bernstein],
        core_target_gain=str(core_gain), bad_numerator_coefficients=list(N),
        positive_repair_numerator_coefficients=[51, 43, 2],
        failure_interval='0<p_a<=1/16', positive_renewal_determinant='(t^2+31*t+208)/256',
        scope='(127i) fails for actual low-extra input; target gain and exact pair repair are positive')


def cheap_pair_repair_obstruction():
    """The two local credit bounds do not cover every individual triple."""
    from verify_feedback_cone import binary_stabilization
    K=[[F(0),F(1,2),F(0),F(0),F(1,2)],
       [F(0),F(0),F(1),F(0),F(0)],
       [F(1),F(0),F(0),F(0),F(0)],
       [F(0),F(1),F(0),F(0),F(0)]]
    p=[F(1,2),F(1,16),F(1,2),F(1,2)];v,a=0,1
    A,R,C,B=map(frozenset,[(0,3),(0,1),(0,2),(0,1,2)])
    supports,P,E,D=binary_stabilization(K,p);pos={S:i for i,S in enumerate(supports)}
    plus=rn.multiply(D[a],P);active=rn.multiply(E,plus)
    starts={S:tuple(int(i in S) for i in range(4)) for S in [A,R,C,B]}
    extra={S:tuple(int(i in S)+int(i==a) for i in range(4)) for S in [A,R,C]}
    raw,count=rn.stopped(list(starts.values())+list(extra.values()),
                         lambda s:rn.arw_step(K,p,s),lambda s:all(x<=0 for x in s))
    terminal={S:tuple(-int(i in S) for i in range(4)) for S in supports}
    for S,start in starts.items():
        assert [raw[start].get(terminal[T],F(0)) for T in supports]==P[pos[S]]
    for S,start in extra.items():
        assert [raw[start].get(terminal[T],F(0)) for T in supports]==plus[pos[S]]
    _,_,_,W=rn.positive_operator(K,p,3)
    w=lambda S:sum(W[tuple(sorted(S)),v].values(),F(0))
    u={S:w(frozenset({v}))-w(S) for S in [R,C,B]}
    assert [u[S] for S in [R,C,B]]==[F(1,4),F(15,64),F(17,64)]
    alpha={};beta={};sigma={};tau={}
    lower=P[pos[A]];upper=plus[pos[A]];mass=upper[pos[B]]
    assert mass==F(1,65)
    for S in [R,C]:
        boundary=[T for T in supports if S<=T and len(T)<=3]
        inverse=rn.inverse([[P[pos[U]][pos[T]] for T in boundary] for U in boundary])
        beta[S]=sum(inverse[boundary.index(B)])*P[pos[S]][pos[S]]
        alpha[S]=lower[pos[S]]/P[pos[S]][pos[S]]
        sigma[S]=[x-alpha[S]*y for x,y in zip(lower,P[pos[S]])]
        tau[S]=rn.multiply([sigma[S]],active)[0]
        assert min(sigma[S])>=0 and min(tau[S])>=0
    assert [beta[S] for S in [R,C]]==[F(5538,10481),F(9200,10481)]
    slack=u[B]-sum(u[S]*beta[S]*alpha[S]*plus[pos[S]][pos[B]]/mass for S in [R,C])
    assert slack==F(-83033,1331087)
    assert sigma[R][pos[frozenset({v})]]==tau[R][pos[R]]==0
    wa=W[tuple(sorted(B)),a].get(tuple(sorted(C)),F(0));assert wa==F(1,8)
    cheap=u[C]*15*wa*tau[C][pos[B]]/mass
    full=u[C]*tau[C][pos[C]]/mass
    assert cheap==F(21825,1040384)
    assert slack+cheap==F(-451458511,10904264704)<0
    assert full==F(174600,1331087)
    assert slack+full==F(721,10481)>0
    ac=frozenset({1,2})
    assert tau[C][pos[C]]==15*(tau[C][pos[R]]/2+tau[C][pos[ac]]/4+tau[C][pos[B]]/8)
    assert sigma[C][pos[frozenset({0})]]==sigma[C][pos[frozenset({2})]]==0
    gain=sum(upper[pos[S]]-lower[pos[S]] for S in supports if v in S)
    assert gain==F(2167809,76982945)>0
    return dict(raw_transient_states=count,normalized_source_budget=str(slack),
                normalized_cheap_credit=str(cheap),cheap_repaired_slack=str(slack+cheap),
                full_pair_repaired_slack=str(slack+full),occupation_gain=str(gain))


def main():
    case=next(c for c in rn.CASES if len(c['K'])==4)
    K=[[F(q) for q in row] for row in case['K']]
    p=list(map(F,case['p']))
    print('archived_four_site',audit(K,p,case),flush=True)
    K=[[0,F(1,2),0,0,F(1,2)],
       [0,0,F(1,3),F(1,3),F(1,3)],
       [F(1,3),0,0,F(1,3),F(1,3)],
       [F(1,2),0,0,0,F(1,2)]]
    K=[[F(q) for q in row] for row in K]
    p=[F(2,3),F(1,2),F(9,10),F(9,10)]
    result=audit(K,p)
    assert F(result['minimum_slack'])==F(432635184999258001002,
                                        131921836276734018522239)
    print('heterogeneous_RN6',result,flush=True)
    print('mixed_contact_obstruction',mixed_contact_obstruction(),flush=True)
    print('union_restoration_obstruction',union_restoration_obstruction(),flush=True)
    print('high_sleep_budget_audit',high_sleep_budget_audit(),flush=True)
    print('quantitative_pair_audit',quantitative_pair_audit(),flush=True)
    print('released_resident_audit',released_resident_audit(),flush=True)
    print('resident_release_transfer_audit',resident_release_transfer_audit(),flush=True)
    print('actual_source_budget_audit',actual_source_budget_audit(),flush=True)
    print('cheap_pair_repair_obstruction',cheap_pair_repair_obstruction(),flush=True)
    print('low_source_feeder_budget_audit',low_source_feeder_budget_audit(),flush=True)
    print('cycle_actual_source_budget_obstruction',cycle_actual_source_budget_obstruction(),flush=True)
    print('rare_source_suppression_audit',rare_source_suppression_audit(),flush=True)
    print('cycle_weighted_budget_obstruction',cycle_weighted_budget_obstruction(),flush=True)
    print('source_maximum_principle_audit',source_maximum_principle_audit(),flush=True)
    print('source_switch_target_audit',source_switch_target_audit(),flush=True)
    print('rare_sleep_expansion_audit',rare_sleep_expansion_audit(),flush=True)
    path=[[F(0)]*5 for _ in range(4)]
    for i in range(4):
        for j in (i-1,i+1):
            path[i][j if 0<=j<4 else 4]+=F(1,2)
    print('rare_sleep_wired_path_audit',rare_sleep_expansion_audit(path),flush=True)


if __name__=='__main__':
    main()

