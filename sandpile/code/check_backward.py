"""Exact rational check of the three-site multi-toppling obstruction."""
from fractions import Fraction as F
from itertools import product

gap = (F(2, 5),) * 3
potential_bound = (F(6, 5), F(8, 5), F(6, 5))

def feasible(second_middle):
    # Endpoint second releases are at least 1/2, so both endpoint counts <= 1.
    # A third middle release exceeds its potential bound in both examples.
    cumulative = ((0, F(9, 10)), (0, F(3, 5), F(3, 5)+second_middle), (0, F(9, 10)))
    good = []
    for k in product(range(2), range(3), range(2)):
        if k == (0, 0, 0):
            continue
        c = tuple(cumulative[i][k[i]] for i in range(3))
        if any(c[i] >= potential_bound[i] for i in range(3)):
            continue
        bc = (c[0]-c[1]/2, c[1]-(c[0]+c[2])/2, c[2]-c[1]/2)
        if all(bc[i] < gap[i] for i in range(3)):
            good.append((k, bc))
    return good

first = feasible(F(3, 5))
second = feasible(F(19, 20))
assert first == [((1, 2, 1), (F(3, 10),)*3)]
assert not second
print('Second middle release 3/5: only inverse vector is (1,2,1), with B C = (3/10,3/10,3/10).')
print('Second middle release 19/20: no nonzero inverse vector. First-card tests agree in both cases.')
