#!/usr/bin/env python3
"""At most two predeclared exact actual-queue discriminators, no sweep.

Question: do the centered kappa-source coefficients in report.md stay nonnegative?
Each run solves one 15- or 21-state rational stationary system, checks all
its equations and active binomial marginal, and compares the coefficient
formula with independently aggregated service-flux weights. No simulation.
The optional particle-Palm check reuses that same law for marked-source
equations, exact backward lifetimes, and two failed coupling shortcuts.
"""
import argparse
import hashlib
import json
import resource
import time
from fractions import Fraction as F
from math import comb, prod
from pathlib import Path

resource.setrlimit(resource.RLIMIT_AS, (1900 * 1024**2,) * 2)
resource.setrlimit(resource.RLIMIT_CPU, (45, 45))
started = time.monotonic()
parser = argparse.ArgumentParser()
parser.add_argument('case', choices=['n4-half', 'n5-fourfifths'])
parser.add_argument('--output', type=Path, help='Create a new JSON certificate; otherwise print it')
parser.add_argument('--particle-palm', action='store_true',
                    help='Also audit the one/two-particle marked resolvents in this same law')
args = parser.parse_args()
n, p = {'n4-half': (4, F(1, 2)), 'n5-fourfifths': (5, F(4, 5))}[args.case]
alpha, q = F(1), 1-p
states = [(i, r) for i in range(n+1) for r in range(n-i+1)]
index = {v: j for j, v in enumerate(states)}
dimension = len(states)
Q = [[F(0) for _ in states] for _ in states]
for col, (i, r) in enumerate(states):
    transitions = []
    if i < n:
        transitions.append(((i+1, max(r-1, 0)), p*(n-i)))
    if i:
        transitions.append(((i-1, r+1), q*i))
    if r:
        transitions.append(((i, r-1), alpha))
    for dest, rate in transitions:
        Q[index[dest]][col] += rate
        Q[col][col] -= rate
assert all(sum(Q[row][col] for row in range(dimension)) == 0
           for col in range(dimension))

system = [row[:] + [F(0)] for row in Q[:-1]]
system.append([F(1)]*dimension + [F(1)])
for col in range(dimension):
    pivot = next(row for row in range(col, dimension) if system[row][col])
    system[col], system[pivot] = system[pivot], system[col]
    scale = system[col][col]
    system[col] = [v/scale for v in system[col]]
    for row in range(dimension):
        if row != col and system[row][col]:
            scale = system[row][col]
            system[row] = [v-scale*w for v, w in zip(system[row], system[col])]
pi = dict(zip(states, (row[-1] for row in system)))
assert sum(pi.values()) == 1 and min(pi.values()) > 0
assert all(sum(Q[row][col]*pi[state] for col, state in enumerate(states)) == 0
           for row in range(dimension))
assert all(sum(pi[i, r] for r in range(n-i+1)) == comb(n, i)*p**i*q**(n-i)
           for i in range(n+1))

kappa = {(s, j): p*(n-s)*pi[j, s-j]
          -alpha*sum(pi[u, s+1-u] for u in range(j+1))
         for s in range(1, n) for j in range(s)}
assert min(kappa.values()) >= 0
omega = {(i, r): kappa[i+r+1, i]
         for i in range(n-1) for r in range(n-1-i)}
for (i, r), value in omega.items():
    h = n-2-i-r
    lhs = ((1+alpha+p*(n-1-i)+q*i)*value
           -p*(n-i)*omega.get((i-1, r+1), F(0))
           -q*(i+1)*omega.get((i+1, r-1), F(0))
           -alpha*omega.get((i, r+1), F(0)))
    rhs = (p*q*(h+1)*pi[i, r+1]+p*alpha*pi[i, r+2]
           +p*p*(h+1)*pi.get((i-1, r+2), F(0)))
    assert lhs == rhs and rhs > 0
rows = []
for k in range(2, n):
    m = n-k
    outside = [(i, r) for i in range(m+1) for r in range(m-i+1)]
    # Direct hypergeometric projection, retaining threshold and service flux.
    N = [{state: F(0) for state in outside} for _ in range(k+1)]
    w = [{state: F(0) for state in outside} for _ in range(k+1)]
    for (I, R), mass in pi.items():
        H = n-I-R
        for a in range(min(k, I)+1):
            for b in range(min(k-a, R)+1):
                c = k-a-b
                if c > H:
                    continue
                weight = mass * F(comb(I,a)*comb(R,b)*comb(H,c), comb(n,k))
                t, state = a+b, (I-a,R-b)
                N[t][state] += weight
                if R:
                    w[t][state] += weight*F(b,R)
    assert sum(sum(v.values()) for v in N) == 1
    Z = {ell: sum(sum(N[t].values()) for t in range(ell,k+1))
         for ell in range(2,k+1)}
    for s in range(1,m+1):
        for j in range(s):
            weights = []
            for t in range(1,k+1):
                lifted = F(comb(n-s-t,k-t),comb(n,k))*sum(
                    comb(J,j)*comb(s+t-J-2,s-j-1)*kappa[s+t-1,J]
                    for J in range(j,j+t))
                direct = p*(m-s+1)*w[t][j,s-j-1] - alpha*sum(
                    w[t][i,s-i] for i in range(j+1))
                assert lifted == direct and lifted >= 0
                weights.append(lifted)
            total = sum(weights)
            for ell in range(2,k+1):
                tail = sum(weights[ell-1:])
                centered = tail-Z[ell]*total
                rows.append({
                    'k':k, 'm':m, 's':s, 'j':j, 'ell':ell,
                    'a_weights':[str(v) for v in weights],
                    'kappa_weighted_tail':str(tail/total) if total else None,
                    'ordinary_block_tail':str(Z[ell]),
                    'centered_coefficient':str(centered),
                    'sign':(centered>0)-(centered<0),
                })

def particle_palm_check():
    """Check actual factorially weighted sources, not a sign conjecture."""
    def falling(value, degree):
        return prod(range(value-degree+1, value+1))

    def solve_exact(matrix, rhs):
        augmented = [row[:] + [value] for row, value in zip(matrix, rhs)]
        size = len(rhs)
        for col in range(size):
            pivot = next(row for row in range(col, size) if augmented[row][col])
            augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
            scale = augmented[col][col]
            augmented[col] = [value/scale for value in augmented[col]]
            for row in range(size):
                if row != col and augmented[row][col]:
                    scale = augmented[row][col]
                    augmented[row] = [v-scale*w for v, w in zip(augmented[row], augmented[col])]
        answer = [row[-1] for row in augmented]
        assert all(sum(v*w for v, w in zip(row, answer)) == value
                   for row, value in zip(matrix, rhs))
        return answer

    records = []
    for d in (1, 2):
        marked = [(i, r, a) for i in range(n-d+1)
                  for r in range(n-d-i+1) for a in range(d+1)]
        where = {state: j for j, state in enumerate(marked)}
        size = len(marked)
        generator = [[F(0)]*size for _ in marked]
        nu, source = [], []
        for col, (i, r, a) in enumerate(marked):
            h, b = n-d-i-r, d-a
            R = r+b
            transitions = []
            if i:
                transitions.append(((i-1, r+1, a), q*i))
            if a:
                transitions.append(((i, r, a-1), q*a))
            if R:
                if r:
                    transitions.extend([((i+1, r-1, a), p*r*(1+F(h,R))),
                                        ((i, r-1, a), alpha*F(r,R))])
                if b:
                    transitions.append(((i, r, a+1), p*b*(1+F(h,R))))
                generator[col][col] -= alpha*F(b,R)
            elif h:
                transitions.append(((i+1, 0, d), p*h))
            for dest, rate in transitions:
                generator[where[dest]][col] += rate
                generator[col][col] -= rate
            nu.append(comb(d,a)*falling(i+a,a)*falling(r+b,b)*pi[i+a,r+b])
            source.append(d*p*(n-i-d+1)*falling(i+d-1,d-1)*pi[i+d-1,0]
                          if r == 0 and a == d else F(0))
        assert all(-sum(generator[row][col]*nu[col] for col in range(size)) == source[row]
                   for row in range(size))
        assert sum(nu) == sum(falling(i+r,d)*weight for (i,r),weight in pi.items() if i+r>=d)
        conditional = [sum(nu[j] for j,(i,r,a) in enumerate(marked) if i+r==s)/sum(nu)
                       for s in range(n-d+1)]
        raw = [sum(falling(i+r,d)*weight for (i,r),weight in pi.items() if i+r==s+d)
               for s in range(n-d+1)]
        assert conditional == [weight/sum(raw) for weight in raw]
        record = {'marks':d, 'states':size, 'stationary_source_residual_pass':True,
                  'factorial_mass_and_conditional_rank_pass':True,
                  'conditional_remaining_rank':list(map(str,conditional)),
                  'all_active_source':[str(source[where[i,0,d]]) for i in range(n-d+1)]}
        backward = [[-generator[col][row] for col in range(size)] for row in range(size)]
        lifetime = solve_exact(backward, [F(1)]*size)
        if n == 4 and p == F(1,2) and alpha == 1:
            selected = [(1,0,1),(1,1,1)] if d == 1 else [(0,0,2),(0,1,2)]
            record['selected_lifetimes'] = {','.join(map(str,state)):str(lifetime[where[state]])
                                            for state in selected}
            record['transformed_common_service_rate'] = str(
                lifetime[where[selected[0]]]/lifetime[where[selected[1]]])
        if d == 2:
            full_time = solve_exact(backward, [F(i+r==n-d) for i,r,a in marked])
            top = where[n-d,0,d]
            unconditioned_full = sum(F(comb(i+r,n-d),comb(n,n-d))*weight
                                     for (i,r),weight in pi.items() if i+r>=n-d)
            largest = full_time[top]/lifetime[top]
            averaged = sum(v*w for v,w in zip(source,full_time))/sum(v*w for v,w in zip(source,lifetime))
            assert averaged == conditional[-1]
            record.update({'largest_entrance_lifetime':str(lifetime[top]),
                           'largest_entrance_full_time':str(full_time[top]),
                           'largest_entrance_normalized_full':str(largest),
                           'unconditioned_outside_full':str(unconditioned_full),
                           'largest_entrance_excess':str(largest-unconditioned_full),
                           'actual_averaged_full':str(averaged),
                           'actual_averaged_gap':str(unconditioned_full-averaged)})
            if n == 4 and p == F(1,2) and alpha == 1:
                assert largest == F(65040,100103) > unconditioned_full
                assert averaged == F(74880,193481) < unconditioned_full
        records.append(record)
    if n == 4 and p == F(1,2) and alpha == 1:
        difference = F(records[1]['transformed_common_service_rate'])-F(records[0]['transformed_common_service_rate'])
        assert difference == F(-877418425986147,105857052385440136) < 0
        records[1]['two_minus_one_mark_service_rate'] = str(difference)
    return records


marked_checks = particle_palm_check() if args.particle_palm else None
result = {
    'scope':'One exact actual-clock law; stronger centered coefficient test only',
    'n':n,'p':str(p),'alpha':str(alpha),
    'ordinary_loopless_lambda':str(F(n+1,n)*p/q),
    'count_state_dimension':dimension,
    'stationary_equations_and_binomial_active_marginal_pass':True,
    'positive_kappa_lift_vs_direct_service_flux_pass':True,
    'kappa_killed_queue_resolvent_pass':True,
    'kappa_reduced_state_count':len(omega),
    'stationary_counts':{f'{i},{r}':str(pi[i,r]) for i,r in states},
    'kappa':{f'{s},{j}':str(v) for (s,j),v in kappa.items()},
    'centered_rows':rows,
    'checked_coefficient_count':len(rows),
    'negative_count':sum(row['sign']<0 for row in rows),
    'all_centered_coefficients_nonnegative':all(row['sign']>=0 for row in rows),
    'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'elapsed_seconds':time.monotonic()-started,
    'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
    'limits':{'cpu_seconds':45,'address_space_mib':1900},
}
if marked_checks is not None:
    result['particle_palm_checks'] = marked_checks
if args.output is None:
    print(json.dumps(result, indent=2))
else:
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in [
        'n', 'p', 'count_state_dimension', 'checked_coefficient_count',
        'negative_count', 'elapsed_seconds', 'peak_rss_mib']}, indent=2))
    print('output='+str(args.output))
