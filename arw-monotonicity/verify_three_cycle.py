"""Exact certificate for all sleep parameters on a directed three-cycle.

Integer polynomial certificate plus independent rational raw-ARW checks.
See drain-inventory.md. This is not a proof for arbitrary directed graphs.
Research implementation by Codex; no independent user-validation claim.
"""

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product, permutations
from math import comb, prod as number_product
from random import Random
import hashlib
import json

from verify_three_site_permanent import Poly
from verify_inventory_corners import stopped, boundary_step

# Poly has eight coordinates; the last two remain zero in this certificate.
p = tuple(Poly.variable(i) for i in range(3))
a = tuple(Poly.variable(i+3) for i in range(3))
q = tuple(1-x for x in p)
ONE = Poly.constant(1)
ZERO = Poly.constant(0)


def prod(values):
    result = ONE
    for value in values:
        result = result*value
    return result


def det(matrix):
    n = len(matrix)
    result = ZERO
    for perm in permutations(range(n)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1,n))
        result += sign*prod(matrix[i][perm[i]] for i in range(n))
    return result


@lru_cache(None)
def idla(state):
    if max(state) <= 1:
        return {tuple(i for i,c in enumerate(state) if c): ONE}
    z = next(i for i,c in enumerate(state) if c > 1)
    out = {}
    for dest,weight in (((z+1)%3,a[z]),(-1,1-a[z])):
        next_state = list(state)
        next_state[z] -= 1
        if dest >= 0:
            next_state[dest] += 1
        for B,c in idla(tuple(next_state)).items():
            out[B] = out.get(B,ZERO)+weight*c
    return out


def full_inverse_row():
    out = {}
    for size in range(4):
        for selected in combinations(range(3),size):
            for choices in product((0,1),repeat=size):
                state = [int(z not in selected) for z in range(3)]
                weight = (-1)**size*prod(q[z] for z in selected)
                for z,bit in zip(selected,choices):
                    weight *= a[z] if bit else 1-a[z]
                    if bit:
                        state[(z+1)%3] += 1
                for B,c in idla(tuple(state)).items():
                    out[B] = out.get(B,ZERO)+weight*c
    return out


def occupation_numerators():
    t = prod(a)
    D1 = 1-prod(q)*t
    Q1 = (p[0],p[0]*q[1]*a[1]*q[2]*a[2],p[0]*q[2]*a[2])
    matrix = [[Poly.constant(int(i==j)) for j in range(3)] for i in range(3)]
    rhs = []
    for i in range(3):
        j = (i+1)%3
        matrix[i][j] -= q[i]*p[j]*a[i]*a[j]
        matrix[i][(i-1)%3] -= q[j]*a[j]
        rhs.append(p[i]*p[j]*int(0 in (i,j))*D1
                   +q[j]*(1-a[j])*Q1[i]
                   +q[i]*p[j]*a[i]*(1-a[j])*Q1[j]
                   +q[i]*p[j]*(1-a[i])*int(j==0)*D1)
    D2 = det(matrix)
    expected = (1-t*(prod(q)+sum((p[i]*q[(i+1)%3]*q[(i+2)%3]
                                 for i in range(3)),ZERO))
                -t*t*prod(p)*prod(q))
    assert D2.terms == expected.terms
    pair = []
    for i in range(3):
        cramer = [[rhs[row] if col==i else matrix[row][col]
                   for col in range(3)] for row in range(3)]
        pair.append(det(cramer))
    row = full_inverse_row()
    D3 = 1-(1-prod(p))*t
    assert row[0,1,2].terms == D3.terms
    numerator = prod(p)*D1*D2
    pair_index = {(0,1):0,(1,2):1,(0,2):2}
    for B,value in row.items():
        if len(B)==1:
            numerator -= value*Q1[B[0]]*D2
        elif len(B)==2:
            numerator -= value*pair[pair_index[B]]
    increments = {x:numerator-D3*pair[(x+1)%3] for x in range(3)}
    return (D1,D2,D3),Q1,pair,numerator,increments


def divide_one_minus(poly,index):
    rows = defaultdict(dict)
    for key,value in poly.terms.items():
        rest = key[:index]+key[index+1:]
        rows[rest][key[index]] = value
    quotient = {}
    for rest,coefficients in rows.items():
        running = 0
        degree = max(coefficients)
        for i in range(degree):
            running += coefficients.get(i,0)
            if running:
                quotient[rest[:index]+(i,)+rest[index:]] = running
        assert coefficients.get(degree,0) == -running
    out = Poly(quotient)
    assert (out*(1-Poly.variable(index))).terms == poly.terms
    return out


def certificate(increments):
    results = []
    for x,numerator in increments.items():
        # Extract p0, and for x != 0 also q0*p_x*a0*a1*a2.
        powers = [1,0,0,0,0,0,0,0]
        if x:
            powers[x] += 1
            powers[3:6] = [1,1,1]
        terms = {}
        for key,value in numerator.terms.items():
            reduced = tuple(i-j for i,j in zip(key,powers))
            assert min(reduced) >= 0
            terms[reduced] = value
        residual = Poly(terms)
        if x:
            residual = divide_one_minus(residual,0)
        degrees = tuple(max(key[i] for key in residual.terms) for i in range(6))
        coefficients = defaultdict(int)
        for key,value in residual.terms.items():
            for output in product(*(range(i,d+1) for i,d in zip(key,degrees))):
                factor = 1
                for d,i,j in zip(degrees,key,output):
                    factor *= comb(d-i,j-i)
                coefficients[output] += value*factor
        positive = {key:value for key,value in coefficients.items() if value}
        assert all(value > 0 for value in positive.values())
        # Check the basis transformation at independent rational arguments.
        values = tuple(F(i+2,2*i+7) for i in range(6))
        expanded = sum(F(value)*number_product(
            z**k*(1-z)**(d-k) for z,k,d in zip(values,key,degrees))
            for key,value in positive.items())
        assert expanded == residual.evaluate(values+(F(0),F(0)))
        digest = hashlib.sha256(json.dumps(sorted(positive.items()),
                                separators=(',',':')).encode()).hexdigest()
        results.append(dict(extra=x,residual_degrees=degrees,
                            positive_coefficients=len(positive),
                            negative_coefficients=0,coefficient_sha256=digest))
    assert [r['positive_coefficients'] for r in results] == [7027,2756,982]
    return dict(check='directed_three_cycle_integer_certificates',
                certificates=results,verified=True)


def occupation_audit(data):
    denominators,singles,pairs,upper,_ = data
    rng = Random(718012)
    for trial in range(24):
        ps = tuple(rng.choice([F(1,4),F(1,2),F(2,3),F(4,5)]) for _ in range(3))
        aa = tuple(rng.choice([F(0),F(1,3),F(3,5),F(1)]) for _ in range(3))
        if all(value==1 for value in aa):
            aa = (F(1,2),)+aa[1:]
        values = ps+aa+(F(0),F(0))
        D1,D2,D3 = [poly.evaluate(values) for poly in denominators]
        assert min(D1,D2,D3) > 0
        expected = {tuple(int(z==i) for z in range(3)):poly.evaluate(values)/D1
                    for i,poly in enumerate(singles)}
        expected.update({tuple(int(z in (i,(i+1)%3)) for z in range(3)):
                         poly.evaluate(values)/(D1*D2) for i,poly in enumerate(pairs)})
        expected[1,1,1] = upper.evaluate(values)/(D1*D2*D3)
        K = tuple(tuple(aa[i] if z==(i+1)%3 else 1-aa[i] if z==3 else F(0)
                        for z in range(4)) for i in range(3))
        laws,_ = stopped(expected,lambda state:boundary_step(K,ps,state),
                         lambda state:all(c<=0 for c in state))
        for eta,answer in expected.items():
            actual = sum(weight for state,weight in laws[eta].items() if state[0]==-1)
            assert answer == actual,(trial,eta,answer,actual)
    return dict(check='directed_three_cycle_raw_occupations',
                rational_cases=24,occupations=168,all_equal=True)


if __name__ == '__main__':
    data = occupation_numerators()
    print(json.dumps(certificate(data[-1])),flush=True)
    print(json.dumps(occupation_audit(data)),flush=True)
