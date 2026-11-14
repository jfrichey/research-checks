"""Exact audit of the three-site pair renewal and target Schur complement.

Run with python3 -B. Uses existing read-only ARW solvers; writes no files.
Three specified heterogeneous kernels, not a numerical screen.
AI-written; identities are audited separately from the general sign certificate.
"""
from fractions import Fraction as F
from math import prod
from pathlib import Path
import sys

import verify_response_renewal as rn
import verify_creator_potentials as cp

V = frozenset(range(3))
SUBS = [frozenset(C) for C in cp.subsets(V)]


def remove_self_loops(K, p):
    return (
        [[F(0) if i == j else K[i][j] / (1-K[i][i])
          for j in range(4)] for i in range(3)],
        [p[i] / (1-(1-p[i])*K[i][i]) for i in range(3)],
    )


def audit(K, p):
    r = [(1-a)/a for a in p]
    k = [row[3] for row in K]
    def minor(C):
        return cp.determinant([[F(i == j)-K[i][j] for j in sorted(C)]
                               for i in sorted(C)])
    det = minor(V)
    d = [minor(V-{i}) for i in range(3)]
    Z = [1+d[i]*(r[j]+r[l]+r[j]*r[l])
         for i in range(3) for j, l in [sorted(V-{i})]]
    ZV = det*sum(prod(r[i] for i in V-C)/minor(C) for C in SUBS)
    H, z = [], []
    for i in range(3):
        row = [F(0)]*3
        for j in V-{i}:
            l = next(iter(V-{i, j}))
            row[j] = r[i]/Z[i]*(K[i][j]*(1+d[i]*r[l])+K[i][l]*K[l][j])
        H.append(row)
        j, l = sorted(V-{i})
        z.append(r[i]/Z[i]*(det/d[i]+(Z[i]-1)*k[i]
                 +d[i]*(r[l]*K[i][j]*k[j]+r[j]*K[i][l]*k[l])))
    IminusH = [[F(i == j)-H[i][j] for j in range(3)] for i in range(3)]
    G = rn.inverse(IminusH)
    beta = [sum(a*b for a, b in zip(row, z))/ZV for row in G]

    starts = {B: tuple(int(i in B) for i in range(3)) for B in SUBS}
    laws, _ = rn.stopped(starts.values(), lambda s: rn.arw_step(K, p, s),
                         lambda s: all(a <= 0 for a in s))
    def pair(B, i):
        return laws[starts[B]].get(tuple(-int(j != i) for j in range(3)), F(0))
    assert beta == [pair(V, i) for i in range(3)]
    assert laws[starts[V]][(-1, -1, -1)] == 1/ZV
    for x in range(3):
        assert [G[i][x]/Z[x] for i in range(3)] == [pair(V-{x}, i) for i in range(3)]

    # Independent killed raw absorption for each last-arrival continuation.
    _, _, _, WW = rn.positive_operator(K, p, 3)
    for i in range(3):
        Ri=tuple(sorted(V-{i}))
        for j in V-{i}:
            assert H[i][j] == r[i]*WW[tuple(sorted(V-{j})),i].get(Ri,F(0))
        assert z[i] == r[i]*WW[(0,1,2),i].get(Ri,F(0))
    assert all(q>=0 for row in G for q in row)
    results = []
    for v in range(3):
        U = sorted(V-{v})
        j, l = U
        inv = rn.inverse([[F(a == b)-K[a][b] for b in U] for a in U])
        h = dict(zip(U, [sum(a*k[b] for a, b in zip(row, U)) for row in inv]))
        invs = rn.inverse([[F(a == b)-(1-p[a])*K[a][b] for b in U] for a in U])
        u = dict(zip(U, [sum(a*(p[b]+(1-p[b])*k[b]) for a, b in zip(row, U))
                         for row in invs]))
        localsets = [frozenset(C) for C in cp.subsets(range(2))]
        th = {C: F(0) for C in localsets}
        if min(h.values()) > 0:
            for D in localsets:
                weight = prod(1-h[a] if b in D else h[a] for b, a in enumerate(U))
                row = cp.two_site_wake_field(K[j][l]*h[l]/h[j], K[l][j]*h[j]/h[l],
                                            [p[j], p[l]], D)
                for C, a in row.items():
                    th[C] += weight*a
        else:
            # A two-site killed block remains well-defined when its Doob potential vanishes.
            mixed={C:tuple(0 if i==v else -1 if U.index(i) in C else 1
                           for i in range(3)) for C in localsets}
            killed,_=rn.stopped(mixed.values(),lambda state:rn.arw_step(K,p,state),
                                lambda state:state[v]>0 or all(q<=0 for q in state))
            moment={C:sum(q for end,q in killed[start].items() if end[v]==0)
                    for C,start in mixed.items()}
            th={C:sum((-1)**len(C-D)*moment[D] for D in localsets if D<=C)
                for C in localsets}
        theta = {frozenset(U[b] for b in C): a for C, a in th.items()}
        t0 = theta[frozenset()]
        w = {frozenset({v}): k[v]+sum(K[v][a]*u[a] for a in U)}
        for a, b in [(j, l), (l, j)]:
            w[frozenset({v, a})] = (
                k[v]+K[v][b]*(t0+theta[frozenset({a})])
                +K[v][a]*(k[a]*u[a]+K[a][b]*t0))
        w[V] = (k[v]+t0*sum(K[v][a]*h[a] for a in U)
                +sum(K[v][a]*theta[frozenset({b})]*k[a]
                     for a, b in [(j, l), (l, j)]))
        for B, val in w.items():
            assert val == sum(WW[tuple(sorted(B)), v].values()), (v, B, val)
        sv = cp.first_sleep(K, p)[v][v]
        coeff = [F(0) if i == v else sv*r[v]*(w[frozenset({v})]-w[V-{i}])
                 for i in range(3)]
        bstar = sv*r[v]*(w[frozenset({v})]-w[V])
        assert min(coeff+[bstar]) >= 0
        for x in V-{v}:
            gain = ((bstar+sum(coeff[i]*sum(G[i][j]*z[j] for j in range(3))
                              for i in range(3)))/ZV
                    -sum(coeff[i]*G[i][x] for i in range(3))/Z[x])
            actual = (sum(q for s, q in laws[starts[V]].items() if s[v] == -1)
                      -sum(q for s, q in laws[starts[V-{x}]].items() if s[v] == -1))
            assert gain == actual
            J = sorted(V-{x})
            M = rn.inverse([[F(a == b)-H[a][b] for b in J] for a in J])
            Mr = [sum(row[t]*H[J[t]][x] for t in range(2)) for row in M]
            Mz = [sum(row[t]*z[J[t]] for t in range(2)) for row in M]
            tau = 1-sum(H[x][a]*q for a, q in zip(J, Mr))
            A = coeff[x]+sum(coeff[a]*q for a, q in zip(J, Mr))
            B = bstar+sum(coeff[a]*q for a, q in zip(J, Mz))
            zbar = z[x]+sum(H[x][a]*q for a, q in zip(J, Mz))
            debt = ZV/Z[x]-zbar
            assert tau > 0
            assert gain == (tau*B-A*debt)/(ZV*tau)
            results.append(dict(x=x, v=v, gain=str(gain), debt_positive=debt > 0,
                                pair_decreases=[i for i in range(3)
                                                if beta[i] < G[i][x]/Z[x]]))
    return dict(pair_system_det=str(cp.determinant(IminusH)), comparisons=results)


if __name__ == '__main__':
    case = rn.CASES[3]
    original_K = [list(map(F, row)) for row in case['K']]
    original_p = list(map(F, case['p']))
    cases = [remove_self_loops(original_K, original_p),
             ([[F(0), F(1, 3), F(1, 2), F(1, 6)],
               [F(1, 5), F(0), F(2, 5), F(2, 5)],
               [F(1, 4), F(1, 2), F(0), F(1, 4)]],
              [F(2, 5), F(3, 7), F(4, 9)])]
    cases.append(([[F(0),F(1),F(0),F(0)],
                   [F(0),F(0),F(1),F(0)],
                   [F(1,2),F(0),F(0),F(1,2)]],
                  [F(1,3),F(2,5),F(3,7)]))
    for index, (K, p) in enumerate(cases):
        print(dict(case=index, **audit(K, p)), flush=True)
