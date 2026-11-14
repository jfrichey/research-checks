"""Unpadded three-cycle inventory positivity: an infinite-family certificate.

All nonempty nested binary inputs and arbitrary cycle counts.
The polynomial tail has degree at most two in the cycle count; its forward
differences give a finite certificate for every remaining count. The separate
small-inventory audit compares with ordinary-buffer/forest enumeration.
Standard library only. Research code by Codex, without independent review.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product, permutations
from math import comb
import hashlib
import json

from verify_three_site_permanent import Poly
from verify_three_cycle import det, prod, ONE, ZERO
from verify_single_drain_inventory import single_dp

p = tuple(Poly.variable(i) for i in range(3))
r = tuple(Poly.variable(i+3) for i in range(3))
# A numerator has total sink degree equal to its number of upper roots.
# Thus d2 can be set to one and restored from that homogeneous degree.
d = (Poly.variable(6),Poly.variable(7),ONE)
q = tuple(1-value for value in p)


def full_row():
    @lru_cache(None)
    def idla(state):
        if max(state) <= 1:
            return {tuple(i for i,count in enumerate(state) if count):ONE}
        z = next(i for i,count in enumerate(state) if count>1)
        out = {}
        for dest,weight in (((z+1)%3,r[z]),(-1,d[z])):
            nxt = list(state)
            nxt[z] -= 1
            if dest>=0:
                nxt[dest] += 1
            for B,value in idla(tuple(nxt)).items():
                out[B] = out.get(B,ZERO)+weight*value
        return out

    out = {}
    for size in range(4):
        for selected in combinations(range(3),size):
            for bits in product((0,1),repeat=size):
                state = [int(z not in selected) for z in range(3)]
                weight = (-1)**size*prod(q[z] for z in selected)
                for z,bit in zip(selected,bits):
                    weight *= r[z] if bit else d[z]
                    if bit:
                        state[(z+1)%3] += 1
                for B,value in idla(tuple(state)).items():
                    out[B] = out.get(B,ZERO)+weight*value
    return out


def inventory_numerators():
    t = prod(r)
    c = tuple(d[i]+r[i]*d[(i+1)%3]+r[i]*r[(i+1)%3]*d[(i+2)%3]
              for i in range(3))
    D1 = 1-prod(q)*t
    single = (p[0]*c[0],p[0]*c[0]*q[1]*r[1]*q[2]*r[2],p[0]*c[0]*q[2]*r[2])
    M = [[Poly.constant(int(i==j)) for j in range(3)] for i in range(3)]
    rhs = []
    for i in range(3):
        j = (i+1)%3
        M[i][j] -= q[i]*p[j]*r[i]*r[j]
        M[i][(i-1)%3] -= q[j]*r[j]
        rhs.append(p[i]*p[j]*c[i]*c[j]*int(0 in (i,j))*D1
                   +(1-t)*(q[j]*d[j]*single[i]
                            +q[i]*p[j]*r[i]*d[j]*single[j]
                            +q[i]*p[j]*d[i]*c[j]*int(j==0)*D1))
    D2 = det(M)
    pairs = []
    for i in range(3):
        matrix = [[rhs[row] if col==i else M[row][col]
                   for col in range(3)] for row in range(3)]
        pairs.append(det(matrix))
    row = full_row()
    D3 = 1-(1-prod(p))*t
    assert row[0,1,2].terms == D3.terms
    upper = prod(p)*prod(c)*D1*D2
    for B,value in row.items():
        if len(B)==1:
            upper -= value*single[B[0]]*D2*(1-t)*(1-t)
        elif len(B)==2:
            upper -= value*pairs[{(0,1):0,(1,2):1,(0,2):2}[B]]*(1-t)
    out = {(tuple(z for z in range(3) if z!=x),x):
           upper-c[x]*D3*pairs[(x+1)%3] for x in range(3)}
    for common,extra in permutations(range(3),2):
        i = {(0,1):0,(1,2):1,(0,2):2}[tuple(sorted((common,extra)))]
        out[(common,),extra] = pairs[i]-c[extra]*single[common]*D2
    return out


def split_patterns(numerators):
    tables = {}
    for (common,x),poly in numerators.items():
        particles = len(common)+1
        roots = tuple(int(i in common or i==x) for i in range(3))
        by_exit = defaultdict(dict)
        for key,value in poly.terms.items():
            exits = (key[6],key[7],particles-key[6]-key[7])
            assert min(exits)>=0
            assert all(key[3+i]+exits[i]-key[3+(i-1)%3]==roots[i] for i in range(3))
            by_exit[exits][key] = value
        assert len(by_exit)==comb(particles+2,2)
        table = {}
        for exits,terms in sorted(by_exit.items()):
            base = tuple(min(key[3+i] for key in terms) for i in range(3))
            coefficients = defaultdict(dict)
            for key,value in terms.items():
                powers = tuple(key[3+i]-base[i] for i in range(3))
                assert len(set(powers))==1
                coefficients[powers[0]][key[:3]+(0,)*5] = value
            table[exits] = base,{j:Poly(terms) for j,terms in coefficients.items()}
        tables[common,x] = table
    return tables


def bernstein(poly,degrees):
    out = defaultdict(int)
    for key,value in poly.terms.items():
        assert all(i<=degree for i,degree in zip(key,degrees)),(key,degrees)
        assert not any(key[3:])
        for output in product(*(range(i,degree+1) for i,degree in zip(key,degrees))):
            factor = value
            for i,j,degree in zip(key,output,degrees):
                factor *= comb(degree-i,j-i)
            out[output] += factor
    return {key:value for key,value in out.items() if value}


def drain_coefficient(coefficients,ell,particles):
    return sum((value*comb(ell-j+particles-1,particles-1)
                for j,value in coefficients.items() if j<=ell),ZERO)


def certificate(tables):
    records = []
    checks = positive = zeros = 0
    A3 = 1-prod(p)
    for (common,x),table in tables.items():
        particles = len(common)+1
        for exits,(base,coefficients) in table.items():
            L = max(coefficients)
            def degree(ell):
                return tuple(a+b+ell for a,b in zip(base,exits))
            G = lambda ell:drain_coefficient(coefficients,ell,particles)
            prefix = ZERO
            for ell in range(L+1):
                prefix = G(ell)+(A3*prefix if particles==3 and x==0 else ZERO)
                array = bernstein(prefix,degree(ell))
                assert all(value>0 for value in array.values()),(x,exits,ell)
                checks += 1
                positive += len(array)
                zeros += (degree(ell)[0]+1)*(degree(ell)[1]+1)*(degree(ell)[2]+1)-len(array)
                records.append((common,x,exits,'prefix',ell,sorted(array.items())))
            tails = tuple(sum(((-1)**(order-j)*comb(order,j)*G(L+j)
                               for j in range(order+1)),ZERO)
                          for order in range(particles))
            for order,tail in enumerate(tails):
                array = bernstein(tail,degree(L))
                assert all(value>0 for value in array.values()),(x,exits,'tail',order)
                checks += 1
                positive += len(array)
                zeros += (degree(L)[0]+1)*(degree(L)[1]+1)*(degree(L)[2]+1)-len(array)
                records.append((common,x,exits,'tail',order,sorted(array.items())))
    digest = hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest()
    return dict(check='unpadded_three_cycle_infinite_inventory_certificate',
                binary_comparisons=9,comparison_exit_patterns=66,bernstein_arrays=checks,
                positive_coefficients=positive,zero_coefficients=zeros,
                negative_coefficients=0,coefficient_sha256=digest,
                all_cycle_counts_certified=True)


def full_coefficients(coefficients,maximum,particles):
    result = [drain_coefficient(coefficients,ell,particles) for ell in range(maximum+1)]
    A1 = prod(q)
    A2 = A1+sum((p[i]*q[(i+1)%3]*q[(i+2)%3] for i in range(3)),ZERO)
    B2 = prod(p)*A1
    factors = [(A1,ZERO),(A2,B2)]
    if particles==3:
        factors.append((1-prod(p),ZERO))
    for linear,quadratic in factors:
        out = []
        for ell,value in enumerate(result):
            out.append(value+(linear*out[ell-1] if ell else ZERO)
                       +(quadratic*out[ell-2] if ell>1 else ZERO))
        result = out
    return result


def enumeration_audit(tables):
    checked = states = 0
    for (common,x),table in tables.items():
        particles = len(common)+1
        for exits,(base,coefficients) in table.items():
            start = min(base)
            minimal = tuple(a-start for a in base)
            values = full_coefficients(coefficients,2-start,particles)
            for loops in range(3):
                forward = tuple(a+loops for a in minimal)
                stock = {(i,(i+1)%3):count for i,count in enumerate(forward) if count}
                stock.update({(i,-1):count for i,count in enumerate(exits) if count})
                degree = tuple(a+b for a,b in zip(forward,exits))
                expected = {} if loops<start else bernstein(values[loops-start],degree)
                _,actual,count = single_dp(list(common),x,stock,3)
                assert expected == actual,(x,exits,loops,expected,actual)
                checked += 1
                states += count
    return dict(check='independent_three_cycle_inventory_enumeration',
                exact_pools=checked,visited_states=states,all_equal=True)


if __name__ == '__main__':
    tables = split_patterns(inventory_numerators())
    print(json.dumps(certificate(tables)),flush=True)
    print(json.dumps(enumeration_audit(tables)),flush=True)
