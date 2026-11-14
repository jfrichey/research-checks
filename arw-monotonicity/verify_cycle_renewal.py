"""Exact audits of the arbitrary-cycle renewal identities in drain-inventory.md.

The proofs cover all cycle lengths. These rational checks compare the explicit
matrices with independently enumerated homogeneous complete-block inverses.
Research code by Codex; no independent user-validation claim.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import prod
from random import Random
import json

from verify_inventory_corners import inverse, multiply


def arc(i,k,n):
    out = []
    while i!=k:
        out.append(i)
        i = (i+1)%n
    return out


def determinant(matrix):
    matrix = [list(row) for row in matrix]
    answer = F(1)
    for i in range(len(matrix)):
        pivot = next(j for j in range(i,len(matrix)) if matrix[j][i])
        if pivot!=i:
            matrix[i],matrix[pivot] = matrix[pivot],matrix[i]
            answer = -answer
        value = matrix[i][i]
        answer *= value
        for j in range(i+1,len(matrix)):
            coefficient = matrix[j][i]/value
            for k in range(i+1,len(matrix)):
                matrix[j][k] -= coefficient*matrix[i][k]
    return answer


def at_least(q,minimum):
    return sum(prod(q[i] if bit else 1-q[i] for i,bit in enumerate(bits))
               for bits in product((0,1),repeat=len(q)) if sum(bits)>=minimum)


def homogeneous_row(A,p,r):
    """No sink jumps: enumerate the signed initial jumps and deterministic IDLA."""
    n = len(p)
    out = {}
    for size in range(len(A)+1):
        for selected in combinations(A,size):
            state = [int(z in A and z not in selected) for z in range(n)]
            weight = F((-1)**size)*prod(1-p[z] for z in selected)
            for z in selected:
                state[(z+1)%n] += 1
                weight *= r[z]
            while max(state)>1:
                z = next(z for z,count in enumerate(state) if count>1)
                state[z] -= 1
                state[(z+1)%n] += 1
                weight *= r[z]
            B = tuple(z for z,count in enumerate(state) if count)
            out[B] = out.get(B,F(0))+weight
    return out


def vacancy_audit(n,p,r):
    q = tuple(1-z for z in p)
    t,P = prod(r),prod(p)
    M = [[F(k==i)-(q[i]*prod(r[z] for z in arc(i,k,n))
                    *prod(p[z] for z in arc((i+1)%n,k,n)) if i!=k else 0)
          for i in range(n)] for k in range(n)]
    sets = [tuple(z for z in range(n) if z!=k) for k in range(n)]
    rows = [homogeneous_row(A,p,r) for A in sets]
    assert M == [[row.get(B,F(0)) for B in sets] for row in rows]
    D = 1-sum(P**(ell-1)*at_least(q,ell+1)*t**ell for ell in range(1,n))
    assert D==determinant(M)>0
    W = [[p[i]*r[i] if k==(i+1)%n else F(0) for i in range(n)] for k in range(n)]
    IminusW = [[F(k==i)-W[k][i] for i in range(n)] for k in range(n)]
    b = [t*q[i]*prod(p[j] for j in range(n) if j!=i) for i in range(n)]
    assert multiply(IminusW,M) == [
        [(1+b[i] if k==i else 0)-(r[i] if k==(i+1)%n else 0)
         for i in range(n)] for k in range(n)]
    inv = inverse(M)
    for k in range(n):
        for i in range(n):
            if k==i:
                predicted = 1+q[i]*sum(P**(ell-1)*at_least(
                    [q[j] for j in range(n) if j!=i],ell)*t**ell
                    for ell in range(1,n))/D
            else:
                path = arc(i,k,n)
                predicted = q[i]*prod(r[z] for z in path)*prod(
                    1+b[z] for z in range(n) if z not in path+[k])/D
            assert inv[k][i]==predicted,(n,k,i)
    return n*n


def pair_audit(n,p,r):
    states = list(combinations(range(n),2))
    index = {A:i for i,A in enumerate(states)}
    rows = [homogeneous_row(A,p,r) for A in states]
    matrix = [[row.get(B,F(0)) for B in states] for row in rows]
    L = [[F(0)]*len(states) for _ in states]
    T = [[F(0)]*len(states) for _ in states]
    for row,(i,j) in enumerate(states):
        prefix = F(1)
        for k in range(i,j):
            L[row][index[k,j]] += prefix
            prefix *= (1-p[k])*r[k]
        dest = (j+1)%n
        first = (1-p[j])*r[j]
        if dest==i:
            first *= r[i]
            dest = (i+1)%n
        T[row][index[tuple(sorted((i,dest)))]] += first
        T[row][index[tuple(sorted((j,(j+1)%n)))]] += p[j]*r[j]*prefix
        assert sum(T[row])<1
    expected = [[F(i==j)-T[i][j] for j in range(len(states))]
                for i in range(len(states))]
    assert multiply(L,matrix)==expected
    assert determinant(L)==1
    assert all(value>=0 for row in multiply(inverse(expected),L) for value in row)
    return len(states)**2


if __name__ == '__main__':
    rng = Random(88421)
    vacancies = pairs = 0
    for n in range(2,9):
        p = tuple(rng.choice([F(1,3),F(1,2),F(3,4)]) for _ in range(n))
        r = tuple(rng.choice([F(1,2),F(2,3),F(1)]) for _ in range(n))
        vacancies += vacancy_audit(n,p,r)
        pairs += pair_audit(n,p,r)
    print(json.dumps(dict(check='arbitrary_cycle_renewal_identities',
                         cycle_lengths=list(range(2,9)),vacancy_inverse_entries=vacancies,
                         pair_matrix_entries=pairs,all_equal=True)),flush=True)
