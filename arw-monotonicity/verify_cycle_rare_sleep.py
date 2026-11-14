"""Exact fixed-ray Taylor audit, using only the maintained raw ARW solver.

Recorded verification status, 2026-10-16. The agent reports both fixed
checks passed; the general proof awaits independent audit (RN section 26).
No automatic run and no parameter search. Sites 0,...,n-1 follow the
linear order after the leak; the last site jumps to 0 with probability c.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import prod
import sys
sys.dont_write_bytecode = True
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_response_renewal as rn


def fixed_ray_coefficients(n, c, rates, degree):
    K = [[F(0)]*(n+1) for _ in range(n)]
    for i in range(n-1):
        K[i][i+1] = F(1)
    K[-1][0], K[-1][-1] = F(c), 1-F(c)
    supports = [frozenset(i for i,b in enumerate(bits) if b)
                for bits in product((0,1), repeat=n)]
    starts = {A: tuple(int(i in A) for i in range(n)) for A in supports}
    states = [s for s in starts.values() if any(a>0 for a in s)]
    pos = {s:i for i,s in enumerate(states)}
    for s in states:
        for end,w in rn.arw_step(K,[F(1,2)]*n,s).items():
            if w and any(a>0 for a in end) and end not in pos:
                pos[end] = len(states)
                states.append(end)
    rows, changes, exits, dexits = [], [], [], []
    for s in states:
        base = rn.arw_step(K,[F(0)]*n,s)
        z = next(i for i,a in enumerate(s) if a>0)
        unit = [F(0)]*n
        unit[z] = F(1)
        changed = rn.arw_step(K,unit,s)
        delta = {end:F(rates[z])*(changed.get(end,F(0))-base.get(end,F(0)))
                 for end in set(base)|set(changed)}
        row, rhs, drhs, drow = {pos[s]:F(1)}, [F(0)]*n, [F(0)]*n, {}
        for end,w in base.items():
            if end in pos:
                row[pos[end]] = row.get(pos[end],F(0))-w
            else:
                for v in range(n):
                    rhs[v] += w*F(end[v]==-1)
        for end,w in delta.items():
            if end in pos:
                drow[pos[end]] = w
            else:
                for v in range(n):
                    drhs[v] += w*F(end[v]==-1)
        rows.append(row)
        changes.append(drow)
        exits.append(rhs)
        dexits.append(drhs)
    previous = rn.old.solve_rows(rows,exits)
    output = {A:[[F(0)]*n for _ in range(degree+1)] for A in supports}
    for k in range(degree+1):
        if k:
            rhs = [[sum(w*previous[j][v] for j,w in drow.items())
                    +(dexits[i][v] if k==1 else 0) for v in range(n)]
                   for i,drow in enumerate(changes)]
            previous = rn.old.solve_rows(rows,rhs)
        for A in supports:
            if A:
                output[A][k] = previous[pos[starts[A]]][:]
    return output, len(states)


def differences(table):
    n = len(next(iter(table.values()))[0])
    answer = {}
    for A in table:
        for x in range(n):
            if x in A:
                continue
            for v in range(n):
                coeff = tuple(upper[v]-lower[v]
                              for upper,lower in zip(table[A|{x}],table[A]))
                answer[A,x,v] = coeff
    return answer


def leaky_target_leading(n, c, rates, A, x):
    """Agent-derived leading term for target n-1, p_i=epsilon*rates[i].

    Returns (degree, coefficient). This is a fixed positive rate ray;
    no uniform claim over all heterogeneous approaches to zero is made.
    """
    A, c, rates = frozenset(A), F(c), tuple(map(F,rates))
    v, d, k = n-1, 1-c, len(A)
    assert 0<c<1 and x not in A and all(a>0 for a in rates)
    r = next(i for i in range(n) if i not in A)
    h = lambda j: F(1) if x<=j else c
    tails = [sum(rates[j:v]) for j in range(v)]

    def no_visit_coefficient():
        roots = sorted(A)
        return rates[v]/d*sum(prod(rates[j] for j in C)
            for C in combinations(range(v),k)
            if all(a<=j for a,j in zip(roots,C)))

    if v not in A and k-r<=1:
        return k+1, no_visit_coefficient()
    remaining = A-set(range(r))-{v}
    b = min(remaining) if remaining else v
    assert b>r
    if r==0:
        correction = c*d**(k-3)*rates[v]*sum(
            h(j)*rates[j]*tails[j] for j in range(b))
    else:
        correction = c*c*d**(k-r-3)*rates[v]*sum(
            h(C[-1])*prod(rates[j] for j in C)*tails[C[0]]
            for C in combinations(range(b),r+1))
    if v not in A and k-r==2:
        correction += no_visit_coefficient()
    return r+3, correction


def before_next_root_leading(n, c, rates, A, x, v):
    """General-target leading term when v precedes the next old root.

    If the hypothesis is absent, return None, not a claimed zero response.
    """
    A, c, rates = frozenset(A), F(c), tuple(map(F,rates))
    d, k = 1-c, len(A)
    r = next(i for i in range(n) if i not in A)
    remaining = A-set(range(r))
    b = min(remaining) if remaining else n
    if v>=b:
        return None
    h = lambda j: F(1) if x<=j else c
    coeff = d**(k-r-1)*sum(
        h(C[-1])*prod(rates[j] for j in C)*(1 if v==C[-1] else c)
        for C in combinations(range(b),r+1) if v in C)
    return r+1, coeff


def single_leak_rare_leading_audit():
    """Two fixed identity audits, not a parameter search.

    The first also certifies all binary comparisons on the half-leak
    four-cycle along the common-p ray. The second checks the stated
    formulas with unequal rates and a different leak probability.
    """
    output = []
    for c,rates in ((F(1,2),(1,1,1,1)), (F(2,3),(2,3,5,7))):
        table,states = fixed_ray_coefficients(4,c,rates,5)
        dif = differences(table)
        leading = {key:next(((k,a) for k,a in enumerate(coeff) if a),None)
                   for key,coeff in dif.items()}
        leaky_checks = prefix_checks = 0
        for (A,x,v),actual in leading.items():
            if v==3:
                predicted = leaky_target_leading(4,c,rates,A,x)
                assert actual==predicted,(A,x,v,actual,predicted)
                assert predicted[1]>0
                leaky_checks += 1
            predicted = before_next_root_leading(4,c,rates,A,x,v)
            if predicted is not None:
                assert actual==predicted,(A,x,v,actual,predicted)
                assert predicted[1]>0
                prefix_checks += 1
        if rates==(1,1,1,1):
            assert len(leading)==128
            assert all(value is not None and value[1]>0 for value in leading.values())
            assert max(value[0] for value in leading.values())==5
            pairs = [value for (A,x,v),value in leading.items() if len(A)==2]
            assert len(pairs)==48 and max(value[0] for value in pairs)==4
            assert leading[frozenset((0,3)),1,3] == (4,F(8))
            assert leading[frozenset((0,3)),2,3] == (4,F(13,2))
        output.append(dict(c=str(c),rates=rates,raw_transient_states=states,
                           leaky_target_formula_checks=leaky_checks,
                           before_next_root_formula_checks=prefix_checks,
                           first_nonzero_positive=sum(bool(a and a[1]>0) for a in leading.values()),
                           maximum_observed_degree=max(a[0] for a in leading.values() if a)))
    return output
