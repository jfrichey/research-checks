#!/usr/bin/env python3
"""Exact audits of directed-cycle queue reversal; standard library only.

RN (130),(131) are checked on all point-mass input/test-function bases.
RN (132) is checked on every nonempty Boolean down-set indicator,
which generates the cone of nonnegative decreasing densities, and on a
strictly positive nonproduct decreasing density. All arithmetic is rational.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import comb, prod


def bern(bits, pars):
    return prod((p if b else 1-p for b, p in zip(bits, pars)), start=F(1))


def binomial(n, p):
    return [(h, F(comb(n, h))*p**h*(1-p)**(n-h))
            for h in range(n+1)]


def pass_law(b, M, p, a):
    """Direct buffer processing and routing, with no change of measure."""
    row = {(M, ()): F(1)}
    for i in range(len(b)):
        out = defaultdict(F)
        for (incoming, z), w in row.items():
            choices = [(b[i], F(1))] if incoming == 0 else [(0, 1-p[i]), (1, p[i])]
            for bit, probability in choices:
                L = incoming+b[i]-bit
                for h, routing in binomial(L, a[i]):
                    out[h, z+(bit,)] += w*probability*routing
        row = dict(out)
    ans = {(z, h): w for (h, z), w in row.items() if w}
    assert sum(ans.values()) == 1
    return ans


def killed_busy_binomial(b, z, M, tilted, h):
    """Expectation in the positive-return formula; kill at an idle buffer."""
    row = {M: F(1)}
    for i in range(len(b)):
        out = defaultdict(F)
        for incoming, w in row.items():
            if incoming == 0:
                continue
            L = incoming+b[i]-z[i]
            if i == len(b)-1:
                out[L] += w
            else:
                for k, routing in binomial(L, tilted[i]):
                    out[k] += w*routing
        row = dict(out)
    return sum((w*comb(L, h) for L, w in row.items() if L >= h), F(0))


def aux_joint(M, p, beta, tilted):
    """Generate the auxiliary (Z, reconstructed B) joint law directly."""
    row = {(M, (), ()): F(1)}
    for i in range(len(p)):
        out = defaultdict(F)
        for (incoming, zs, bs), w in row.items():
            for z in (0, 1):
                wz = beta[i] if z else 1-beta[i]
                old = [(z, F(1))] if incoming == 0 else [(0, 1-p[i]), (1, p[i])]
                for b, wb in old:
                    L = incoming+b-z
                    routes = [(0, F(1))] if i == len(p)-1 else binomial(L, tilted[i])
                    for nxt, wr in routes:
                        out[nxt, zs+(z,), bs+(b,)] += w*wz*wb*wr
        row = dict(out)
    ans = defaultdict(F)
    for (_, z, b), w in row.items():
        ans[z, b] += w
    assert sum(ans.values()) == 1
    return ans


def suffix_audit(a):
    p = (F(1,3), F(2,5), F(3,7))
    m = len(p)
    states = list(product((0,1), repeat=m))
    gamma = tuple(1-prod(a[i:]) for i in range(m))
    beta = tuple(p[i]/(p[i]+(1-p[i])*gamma[i]) for i in range(m))
    tilted = tuple(a[i]*gamma[i+1]/gamma[i] for i in range(m-1))
    assert all(0 <= x <= 1 for x in tilted)
    wb = {b: bern(b, beta) for b in states}
    wp = {b: bern(b, p) for b in states}
    downsets = []
    for mask in range(1, 1 << len(states)):
        D = {b for k, b in enumerate(states) if mask & (1 << k)}
        if all(c in D for b in D for c in states if all(x <= y for x,y in zip(c,b))):
            downsets.append(D)
    assert len(downsets) == 19
    raw_densities = [{b: F(b in D) for b in states} for D in downsets]
    raw_densities.append({b: F(1+3*(b[0]==0)+2*(b[1]==b[2]==0)
                               +(b[0]==b[2]==0)) for b in states})
    densities = []
    for f in raw_densities:
        total = sum(wb[b]*f[b] for b in states)
        densities.append({b: f[b]/total for b in states})
    counts = dict(positive_identity=0, zero_identity=0, positive_density_covers=0,
                  zero_target_budget=0, product_reset=0, no_contact_product=0)
    # Independently audit (128): every site has its own initial active resident.
    active = {(0, ()): F(1)}
    for i in range(m):
        nxt = defaultdict(F)
        for (incoming, zs), w in active.items():
            for z in (0,1):
                for h, wr in binomial(incoming+1-z, a[i]):
                    nxt[h, zs+(z,)] += w*(p[i] if z else 1-p[i])*wr
        active = dict(nxt)
    for z in states:
        expected = prod((p[i] if z[i] else (1-p[i])*gamma[i]
                         for i in range(m)), start=F(1))
        assert active.get((0,z), F(0)) == expected
        counts['no_contact_product'] += 1
    for M in range(5):
        laws = {b: pass_law(b, M, p, a) for b in states}
        aux = aux_joint(M, p, beta, tilted)
        for b in states:
            for z in states:
                lhs = wb[b]*laws[b].get((z,0), F(0))
                rhs = gamma[0]**M*aux[z,b]
                assert lhs == rhs, ('zero identity', a, M, b, z, lhs, rhs)
                counts['zero_identity'] += 1
            for h in range(1, M+m+1):
                density = {}
                for z in states:
                    lhs = wb[b]*laws[b].get((z,h), F(0))
                    rhs = wb[z]*gamma[0]**M*(a[-1]/gamma[-1])**h*wp[b] \
                          *killed_busy_binomial(b,z,M,tilted,h)
                    assert lhs == rhs, ('positive identity', a, M, b, z, h, lhs, rhs)
                    counts['positive_identity'] += 1
                    density[z] = laws[b].get((z,h), F(0))/wb[z]
                for z in states:
                    for i in range(m):
                        if z[i] == 0:
                            zz = list(z); zz[i] = 1
                            assert density[z] >= density[tuple(zz)]
                            counts['positive_density_covers'] += 1
        for z in states:
            assert sum(wb[b]*laws[b].get((z,0), F(0)) for b in states) == gamma[0]**M*wb[z]
            counts['product_reset'] += 1
        for f in densities:
            zero = {z: sum(wb[b]*f[b]*laws[b].get((z,0), F(0)) for b in states)
                    for z in states}
            for T in states:
                probability = sum(w for z,w in zero.items() if all(not t or z[i] for i,t in enumerate(T)))
                budget = prod((beta[i] for i,t in enumerate(T) if t), start=F(1))*sum(zero.values())
                assert probability <= budget, ('target budget', a,M,f,T,probability,budget)
                counts['zero_target_budget'] += 1
    return dict(p=list(map(str,p)), a=list(map(str,a)), gamma=list(map(str,gamma)),
                beta=list(map(str,beta)), tilted=list(map(str,tilted)),
                incoming_batches=list(range(5)), nonempty_downsets=19,
                extra_positive_nonproduct_density=True, checks=counts)


def physical_audit(a):
    import verify_response_renewal as rn
    from verify_creator_potentials import wake_field_coefficients
    n = 4
    p = (F(1,3), F(2,5), F(3,7), F(4,9))
    single_leak = all(x == 1 for x in a[:-1])
    c = prod(a)
    K = [[F(0)]*(n+1) for _ in range(n)]
    for i in range(n):
        K[i][(i+1)%n] = a[i]
        K[i][n] = 1-a[i]
    starts = list(product((-1,0,1), repeat=n))
    # Raw one-instruction ARW includes ineffective sleeps at piles.
    laws, transient = rn.stopped(starts, lambda s: rn.arw_step(K,p,s),
                                 lambda s: all(x <= 0 for x in s))
    assert all(sum(row.values()) == 1 and all(w >= 0 for w in row.values())
               for row in laws.values())
    sets = [frozenset(i for i in range(n) if mask & (1<<i)) for mask in range(1<<n)]
    V = sets[-1]
    def moment(start, T):
        return sum((w for end,w in laws[start].items() if all(end[i] == -1 for i in T)), F(0))
    mixed = {S: tuple(-1 if i in S else 1 for i in range(n)) for S in sets}
    active = {R: tuple(int(i in R) for i in range(n)) for R in sets}
    theta = {(T,C): sum((-1)**len(C-S)*moment(mixed[S],T) for S in sets if S <= C)
             for T in sets for C in sets}
    P = [[sum((w for end,w in laws[active[R]].items()
               if frozenset(i for i,x in enumerate(end) if x == -1) == B), F(0))
          for B in sets] for R in sets]
    algebra, _, _ = wake_field_coefficients(K,p,sets,P)
    counts = defaultdict(int)
    for T in sets:
        assert sum(theta[T,C] for C in sets) == 1
        counts['field_row_normalization'] += 1
        for C in sets:
            odds = prod(((1-p[i])/p[i] for i in C), start=F(1))
            assert theta[T,C] == odds*algebra[C][sets.index(T)]
            counts['raw_vs_operator_field'] += 1
            if single_leak:
                assert theta[T,C] >= 0
                counts['nonnegative_field_coefficients'] += 1
            if C and single_leak:
                j = max(C)
                factor = prod(((1-p[i])/p[i] for i in C if i != j), start=F(1))
                assert theta[T,C] == factor*theta[T,frozenset((j,))]
                counts['marked_reduction'] += 1
        for j in range(n):
            gain = moment(active[V],T)-moment(active[V-{j}],T)
            assert theta[T,frozenset((j,))] == (1-p[j])/p[j]*(1-c)*gain
            assert gain >= 0
            counts['singleton_full_minus_one_identity'] += 1
    for j in range(n):
        other = sorted(V-{j})
        gamma = {}
        for i in other:
            route = F(1); k = i
            while k != j:
                route *= a[k]; k = (k+1)%n
            gamma[i] = 1-route
        normal = prod((p[i]+(1-p[i])*gamma[i] for i in other), start=F(1))
        beta = {i: p[i]/(p[i]+(1-p[i])*gamma[i]) for i in other}
        weights = {}
        for B in sets:
            if j not in B:
                weights[B] = prod((beta[i] if i in B else 1-beta[i] for i in other), start=F(1))
        start = active[V-{j}]
        stopped, _ = rn.stopped([start], lambda s: rn.arw_step(K,p,s),
                                lambda s: s[j] > 0 or all(x <= 0 for x in s))
        for B,w in weights.items():
            end = tuple(-int(i in B) for i in range(n))
            assert stopped[start].get(end,F(0)) == normal*w
            counts['no_contact_product_entries'] += 1
        for T in sets:
            U = sum(w*moment(tuple(1 if i==j else -int(i in B) for i in range(n)),T)
                    for B,w in weights.items())
            budget = prod((beta[i] for i in T if i != j), start=F(1))
            assert theta[T,frozenset((j,))] == normal*(budget-U)
            assert U <= budget
            counts['physical_product_theta_identity_and_budget'] += 1
    return dict(p=list(map(str,p)), a=list(map(str,a)), raw_binary_initial_states=len(starts),
                raw_transient_states=transient, checks=dict(counts),
                target_3_singletons=[str(theta[frozenset((3,)),frozenset((j,))]) for j in range(n)],
                coefficient_mask_order='rows target T, columns marked C; bit i is vertex i',
                coefficients=[[str(theta[T,C]) for C in sets] for T in sets])



def branching_queue_audit():
    """DAG zero-return identity and a positive-return density obstruction."""
    import verify_response_renewal as rn
    p=(F(99,100),)*2
    routes=((F(9,10),F(1,20),F(1,20)),(F(0),F(1,2),F(1,2)))
    h=(F(1,2),)*2
    beta=tuple(p[i]/(p[i]+(1-p[i])*h[i]) for i in range(2))
    bits=list(product((0,1),repeat=2))
    mu={b:bern(b,beta) for b in bits}
    def route(L,probs):
        a,r,s=probs
        return [(k,H,comb(L,k)*comb(L-k,H)*a**k*r**H*s**(L-k-H))
                for k in range(L+1) for H in range(L-k+1)]
    def batch(b,M):
        row={(M,0,()):F(1)}
        for i in range(2):
            out=defaultdict(F)
            for (incoming,H,zs),mass in row.items():
                if incoming==0:
                    out[0,H,zs+(b[i],)]+=mass
                    continue
                for z in (0,1):
                    for nxt,J,w in route(incoming+b[i]-z,routes[i]):
                        out[nxt,H+J,zs+(z,)]+=mass*(p[i] if z else 1-p[i])*w
            row=dict(out)
        assert all(nxt==0 for (nxt,_,_),w in row.items() if w)
        result={(z,H):w for (_,H,z),w in row.items() if w}
        assert sum(result.values())==1
        return result
    def step(state):
        config,H=state
        i=next(i for i,x in enumerate(config) if x>0)
        out=defaultdict(F)
        end=list(config)
        if end[i]==1:
            end[i]=-1
        out[tuple(end),H]+=p[i]
        for dest,w in enumerate(routes[i]):
            if not w:
                continue
            end=list(config);end[i]-=1
            if dest==0:
                end[i+1]=2 if end[i+1]==-1 else end[i+1]+1
            out[tuple(end),H+int(dest==1)]+=(1-p[i])*w
        return dict(out)
    starts={(M,b):((M+b[0] if M else -b[0],-b[1]),0)
            for M in range(4) for b in bits}
    raw,_=rn.stopped(starts.values(),step,
                     lambda state:all(x<=0 for x in state[0]))
    downsets=[{b for i,b in enumerate(bits) if mask>>i&1}
              for mask in range(1,1<<len(bits))]
    downsets=[D for D in downsets if all(c in D for b in D for c in bits
                                        if all(x<=y for x,y in zip(c,b)))]
    assert len(downsets)==5
    entries=budgets=0
    for M in range(4):
        laws={b:batch(b,M) for b in bits}
        aux=aux_joint(M,p,beta,(F(9,10),))
        for b in bits:
            physical={(tuple(int(x==-1) for x in end),H):w
                      for (end,H),w in raw[starts[M,b]].items()}
            assert physical==laws[b]
            for z in bits:
                assert mu[b]*laws[b].get((z,0),0)==h[0]**M*aux[z,b]
                entries+=1
        for D in downsets:
            total=sum(mu[b] for b in D)
            zero={z:sum((mu[b]*laws[b].get((z,0),0) for b in D),F(0))/total
                  for z in bits}
            for T in bits:
                value=sum(w for z,w in zero.items()
                          if all(not t or z[i] for i,t in enumerate(T)))
                bound=prod(beta[i] for i,t in enumerate(T) if t)*sum(zero.values())
                assert value<=bound
                budgets+=1
        if M==1:
            output=defaultdict(F)
            for b in bits:
                if not b[1]:
                    for key,w in laws[b].items():
                        output[key]+=mu[b]*w/(1-beta[1])
            density=[output[z,1]/mu[z] for z in ((0,0),(1,0))]
            assert density==[F(2014477,500000),F(2147409,200000)]
            increase=density[1]-density[0]
            assert increase==F(6708091,1000000)>0
            assert sum(w for (z,H),w in output.items() if H==1)==F(2334727,39800000)
            # The failure persists for a strictly positive decreasing density.
            epsilon=F(1,100)
            positive=defaultdict(F)
            for b in bits:
                f=epsilon+(1-epsilon)*F(not b[1])/(1-beta[1])
                for key,w in laws[b].items():
                    positive[key]+=mu[b]*f*w
            assert positive[(1,0),1]/mu[1,0]>positive[(0,0),1]/mu[0,0]
            # Ordinary stochastic domination is too weak even for H=0.
            below={(0,0):F(1,39601),(0,1):F(198,199),
                   (1,0):F(198,39601),(1,1):F(0)}
            assert sum(below.values())==1
            for mask in range(1,1<<len(bits)):
                U={b for i,b in enumerate(bits) if mask>>i&1}
                if all(c in U for b in U for c in bits
                       if all(x<=y for x,y in zip(b,c))):
                    assert sum(below[b] for b in U)<=sum(mu[b] for b in U)
            zero={z:sum(below[b]*laws[b].get((z,0),0) for b in bits) for z in bits}
            excess=sum(w for z,w in zero.items() if z[1])-beta[1]*sum(zero.values())
            assert excess==F(368899839,78805990000)>0
            first={z:sum(mu[b]*laws[b].get((z,3),0) for b in bits) for z in bits}
            assert first[0,0]==F(970299,79202000000)
            assert all(not w for z,w in first.items() if z!=(0,0))
        if M==3:
            second=[laws[0,0].get((z,1),F(0)) for z in bits]
            assert second==[F(n,80000000) for n in (3297,395604,788040,38811960)]
            conditional={z:w/(sum(second)*mu[z]) for z,w in zip(bits,second)}
            reachable_increase=conditional[1,0]-conditional[0,0]
            assert reachable_increase==F(27047483,39998901)>0
    return dict(dag_zero_return_entries=entries,dag_target_budgets=budgets,
                positive_return_density_increase=str(increase),
                stochastic_order_zero_return_budget_excess=str(excess),
                reachable_two_pass_density_increase=str(reachable_increase))



def partial_root_queue_audit():
    """One exact partial-root law, its positive tilt and zero-pass composition."""
    import verify_response_renewal as rn
    p=(F(1,2),)*4
    K=[[F(0)]*5 for _ in range(4)]
    for i in range(4):
        K[i][(i+1)%4]=F(1,2); K[i][4]=F(1,2)
    h=(F(7,8),F(3,4),F(1,2))
    beta=(F(8,15),F(4,7),F(2,3))
    t=tuple((F(1,2)+F(1,2)*x)/x for x in h)
    bits=list(product((0,1),repeat=3))
    starts={M:(0,1+M,0,1) for M in range(3)}
    raw,_=rn.stopped(starts.values(),lambda state:rn.arw_step(K,p,state),
                     lambda state:state[0]>0 or all(x<=0 for x in state))
    nu={M:{tuple(int(x==-1) for x in end[1:]):w
           for end,w in raw[start].items() if end[0]<=0}
        for M,start in starts.items()}
    assert sum(nu[0].values())==F(93,128)
    assert nu[0][0,0,0]==F(11,128)
    assert nu[0][1,0,0]==F(1,8)
    assert sum(prod(h[i] for i,z in enumerate(B) if z)*w
               for B,w in nu[0].items())==h[0]*h[2]
    norm=sum(nu[0].values())
    density={B:nu[0].get(B,0)/(norm*bern(B,beta)) for B in bits}
    assert density[0,0,0]==F(55,31)<F(70,31)==density[1,0,0]
    marginal=sum(w for B,w in nu[0].items() if B[0])/norm
    assert marginal==F(16,31)<beta[0]

    # Independent visited-set tilt: all-active path input under (Khat,beta).
    Khat=[[F(0),F(3,7),F(0),F(4,7)],
          [F(0),F(0),F(1,3),F(2,3)],
          [F(0),F(0),F(0),F(1)]]
    initial=((1,0,1),5)
    def tilted_step(state):
        config,visited=state
        return {(end,visited | sum(1<<i for i,x in enumerate(end) if x>0)):w
                for end,w in rn.arw_step(Khat,beta,config).items()}
    tilted,_=rn.stopped([initial],tilted_step,
                        lambda state:all(x<=0 for x in state[0]))
    weighted=defaultdict(F)
    for (end,visited),w in tilted[initial].items():
        B=tuple(int(x==-1) for x in end)
        weighted[B]+=h[0]*h[2]*w*prod(t[i] for i in range(3) if visited>>i&1)
    assert dict(weighted)==nu[0]
    compositions=budgets=prefixes=0
    for M in range(3):
        composed=defaultdict(F)
        for B,w in nu[0].items():
            for (Z,H),q in pass_law(B,M,(F(1,2),)*3,(F(1,2),)*3).items():
                if H==0:
                    composed[Z]+=w*q
        assert dict(composed)==nu[M]
        compositions+=1
        total=sum(composed.values())
        for T in bits:
            value=sum(w for B,w in composed.items()
                      if all(not z or B[i] for i,z in enumerate(T)))
            assert value<=total*prod(beta[i] for i,z in enumerate(T) if z)
            budgets+=1
        for i in range(3):
            for prefix in product((0,1),repeat=i):
                denominator=sum(w for B,w in composed.items() if B[:i]==prefix)
                numerator=sum(w for B,w in composed.items() if B[:i]==prefix and B[i])
                assert numerator<=beta[i]*denominator
                prefixes+=1
    # The actual partial-root marginal cannot replace the product budget.
    K2=[[F(0)]*5 for _ in range(4)]
    for i,a in enumerate((F(1),F(1),F(1),F(2,3))):
        K2[i][(i+1)%4]=a; K2[i][4]=1-a
    p2=(F(1,3),F(2,5),F(3,7),F(4,9))
    eta=(1,0,1,0)
    active={i:tuple(eta[z]+int(z==i) for z in range(4)) for i in range(4)}
    mixed={i:tuple(-1 if z==i else eta[z] for z in range(4)) for i in (1,3)}
    starts2=[eta]+list(active.values())+list(mixed.values())
    raw2,_=rn.stopped(starts2,lambda state:rn.arw_step(K2,p2,state),
                      lambda state:all(x<=0 for x in state))
    def moment(start):
        return sum(w for end,w in raw2[start].items() if end[2]==-1)
    theta={i:moment(mixed[i])-moment(active[i]) if i in mixed else F(0)
           for i in range(4)}
    assert theta[1]==F(-15889009571721264,1490750968359917015)<0
    gain={i:moment(active[i])-moment(eta) for i in range(4)}
    for i in range(4):
        assert theta[i]==(1-p2[i])/p2[i]*(gain[i]-sum(K2[i][j]*gain[j] for j in range(4)))
    return dict(partial_root_tilt=True,zero_pass_compositions=compositions,
                joint_budgets=budgets,prefix_bounds=prefixes,
                nondecreasing_density_pair=('55/31','70/31'),
                actual_marginal='16/31',product_bound='8/15',
                negative_partial_coefficient=str(theta[1]),source_poisson_entries=4)


if __name__ == '__main__':
    print(partial_root_queue_audit(),flush=True)
    print(branching_queue_audit(),flush=True)
    for name,a in [('single_leak',(F(1),F(1),F(2,3))),
                   ('varied_leaks',(F(3,4),F(2,3),F(3,5)))]:
        print(name, suffix_audit(a), flush=True)
    for name,a in [('single_leak',(F(1),F(1),F(1),F(2,3))),
                   ('varied_leaks',(F(3,4),F(2,3),F(3,5),F(4,5)))]:
        result=physical_audit(a)
        print(name, {k:v for k,v in result.items() if k != 'coefficients'}, flush=True)
