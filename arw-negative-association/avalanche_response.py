"""Three declared exact discriminators, not a parameter sweep.

Candidate: reverse avalanche loss decreases stochastically with terminal
rank for homogeneous complete ARW. This would imply two consecutive losses
are negatively quadrant dependent via forward loss monotonicity. n=4,p=1/4
already defeats TP2 of the addition kernel, and n=6,p=1/2 is beyond the
known complete-graph full-NA range. A failure parks this sufficient route;
a pass warrants a structural argument, not more examples.

The ordinary three-site sink-end path at rate one separately tests whether
conditioning on each exact loss >=2 lowers the terminal law. Its failure
rules out a universal concavity shortcut for the loss PGF, not NA itself.
"""
from fractions import Fraction as F
from math import comb
import json
import resource
import time
import argparse
import hashlib
from pathlib import Path

resource.setrlimit(resource.RLIMIT_AS, (1900*1024**2,)*2)
resource.setrlimit(resource.RLIMIT_CPU, (45,45))

def solve(a, b):
    n, k = len(a), len(b[0])
    m = [row[:] + rhs[:] for row,rhs in zip(a,b)]
    for c in range(n):
        v = next(r for r in range(c,n) if m[r][c])
        m[c],m[v] = m[v],m[c]
        pivot=m[c][c]
        m[c]=[v/pivot for v in m[c]]
        for r in range(n):
            if r!=c and m[r][c]:
                z=m[r][c]
                m[r]=[x-z*y for x,y in zip(m[r],m[c])]
    result=[row[n:] for row in m]
    assert all(sum(a[i][j]*result[j][h] for j in range(n))==b[i][h]
               for i in range(n) for h in range(k))
    return result

def case(n,p,alpha):
    q=1-p
    unit=lambda i:[F(j==i) for j in range(n+1)]
    R=[[[F(j==0) for j in range(n+1)]]]
    for b in range(1,n+2):
        d=min(b,n+1)
        a=[[F(0) for _ in range(d)] for _ in range(d)]
        rhs=[[alpha*R[b-1][y][j] for j in range(n+1)] for y in range(d)]
        for y in range(d):
            birth,death=p*(n-y),q*y
            a[y][y]=alpha+birth+death
            if y:
                a[y][y-1]=-death
            if y+1<d:
                a[y][y+1]=-birth
            elif b<=n:
                rhs[y][b]+=birth
        rows=solve(a,rhs)
        if b<=n:
            rows.append(unit(b))
        assert all(sum(row)==1 and min(row)>=0 for row in rows)
        R.append(rows)
    # H is hub addition. Actual uniform-site addition has an initial local
    # refresh. Its rank kernel is T=((n+alpha)H-alpha I)/n.
    H=[[R[j+1][j][i] for j in range(n+1)] for i in range(n+1)]
    T=[[(n+alpha)/n*H[i][j]-alpha/n*F(i==j) for j in range(n+1)]
       for i in range(n+1)]
    for j in range(n+1):
        # Independent local-refresh construction of the same column.
        explicit=[p*F(n-j,n)*F(i==j+1)
                  +(p*F(j,n)+q*F(n-j,n))*H[i][j]
                  +(q*F(j,n)*sum(H[i][k]*H[k][j-1] for k in range(n+1)) if j else 0)
                  for i in range(n+1)]
        assert explicit==[T[i][j] for i in range(n+1)]
        assert sum(T[i][j] for i in range(n+1))==1
        assert min(T[i][j] for i in range(n+1))>=0
    r=[sum(F(comb(n,y))*p**y*q**(n-y)*R[n][y][i] for y in range(n+1))
       for i in range(n+1)]
    assert sum(r)==1 and min(r)>0
    assert all(sum(T[i][j]*r[j] for j in range(n+1))==r[i] for i in range(n+1))
    reverse=[[sum(r[k]*T[j][k] for k in range(n+1) if 1+k-j==l)/r[j]
              for l in range(n+2)] for j in range(n+1)]
    gaps=[(sum(reverse[j+1][:l+1])-sum(reverse[j][:l+1]),j,l)
          for j in range(n) for l in range(n+1)]
    losses=[[F(0) for _ in range(n+2)] for _ in range(n+2)]
    for i in range(n+1):
        for j in range(n+1):
            for k in range(n+1):
                mass=r[i]*T[j][i]*T[k][j]
                if mass:
                    losses[1+i-j][1+j-k]+=mass
    marg=[sum(row) for row in losses]
    assert sum(marg)==1 and sum(i*v for i,v in enumerate(marg))==1
    assert [sum(losses[i][j] for i in range(n+2)) for j in range(n+2)]==marg
    tails=[(sum(losses[i][j] for i in range(a,n+2) for j in range(b,n+2))
            -sum(marg[a:])*sum(marg[b:]),a,b)
           for a in range(1,n+2) for b in range(1,n+2)]
    exact_loss_gaps=[]
    for ell in range(2,n+2):
        for threshold in range(1,n+1):
            cov=sum(r[j]*reverse[j][ell] for j in range(threshold,n+1))-sum(r[threshold:])*marg[ell]
            exact_loss_gaps.append((cov,threshold,ell))
    def G(z): return sum(v*z**i for i,v in enumerate(marg))
    def G2(z,w): return sum(losses[i][j]*z**i*w**j for i in range(n+2) for j in range(n+2))
    t=u=F(1,4); s=1-t-u
    HH=(G(1-t)*G(1-u)+G2(s,s)-G2(1-t,s)-G2(1-u,s))/(t*u)
    return {'n':n,'p':str(p),'alpha':str(alpha),'rank_law':list(map(str,r)),
            'addition_columns_normalized':True,'independent_refresh_identity':True,
            'stationarity_from_saturated_sampler':True,
            'reverse_loss_distributions':[[str(v) for v in row] for row in reverse],
            'minimum_reverse_cdf_gap':list(map(str,min(gaps))),
            'reverse_order_failures':[[str(x),j,l] for x,j,l in gaps if x<0],
            'maximum_serial_loss_tail_covariance':list(map(str,max(tails))),
            'serial_NQD_failures':[[str(x),a,b] for x,a,b in tails if x>0],
            'exact_loss_rank_domination_failures':[[str(x),a,b] for x,a,b in exact_loss_gaps if x>0],
            'H_at_t_u_quarter':str(HH)}

def path_concavity_obstruction():
    from exact import absorption, graph_kernel, stable, add, increasing_events, sp
    n=3
    kernel=graph_kernel([(1,3),(0,2),(1,)])
    starts=[(1,)*n]+[add(stable(mask,n),v) for v in range(n) for mask in range(1<<n)]
    laws,count=absorption(starts,kernel,(sp.Rational(1,2),)*n)
    other,count2=absorption(starts,kernel,(sp.Rational(1,2),)*n,reverse=True)
    assert laws==other
    mu=laws[0,:]
    failures=[]
    checks=0
    rank_reverse_failures=[]
    serial_failures=[]
    for v in range(n):
        addition=laws[1+v*8:1+(v+1)*8,:]
        assert mu*addition==mu
        loss_joint=[[sum(mu[x]*addition[x,y] for x in range(8)
                         if 1+x.bit_count()-y.bit_count()==ell)
                     for y in range(8)] for ell in range(n+2)]
        for ell in range(2,n+2):
            mass=loss_joint[ell]
            pl=sum(mass)
            for event in increasing_events(tuple(range(n))):
                cov=sum(mass[y] for y in event)-pl*sum(mu[y] for y in event)
                checks+=1
                if cov>0:
                    failures.append({'covariance':str(cov),'entrance':v,'loss':ell,
                                     'event_true_masks':list(event)})
        ranks=[sum(mu[y] for y in range(8) if y.bit_count()==j) for j in range(n+1)]
        conditional=[[sum(loss_joint[ell][y] for y in range(8) if y.bit_count()==j)/ranks[j]
                      for ell in range(n+2)] for j in range(n+1)]
        for j in range(n):
            for ell in range(n+1):
                gap=sum(conditional[j+1][:ell+1])-sum(conditional[j][:ell+1])
                if gap<0:
                    rank_reverse_failures.append({'entrance':v,'lower_rank':j,
                                                  'loss_threshold':ell,'cdf_gap':str(gap)})
        pair=[[sp.S.Zero for _ in range(n+2)] for _ in range(n+2)]
        for x in range(8):
            for y in range(8):
                for z in range(8):
                    mass=mu[x]*addition[x,y]*addition[y,z]
                    if mass:
                        pair[1+x.bit_count()-y.bit_count()][1+y.bit_count()-z.bit_count()]+=mass
        marginal=[sum(row) for row in pair]
        assert sum(marginal)==1 and sum(ell*w for ell,w in enumerate(marginal))==1
        assert [sum(pair[i][j] for i in range(n+2)) for j in range(n+2)]==marginal
        for a in range(1,n+2):
            for b in range(1,n+2):
                cov=sum(pair[i][j] for i in range(a,n+2) for j in range(b,n+2))-sum(marginal[a:])*sum(marginal[b:])
                if cov>0:
                    serial_failures.append({'entrance':v,'thresholds':[a,b],'covariance':str(cov)})
    assert failures==[{'covariance':'12293539433/18807104329000','entrance':2,'loss':2,
                      'event_true_masks':[1,2,3,4,5,6,7]}]
    return {'graph':'sink--0--1--2','common_sleep_rate':'1','pure_entrances':'all',
            'stationary_law':list(map(str,mu)),
            'transient_counts':[count,count2],
            'absorption_equations_both_priorities_and_all_addition_stationarity':True,
            'configuration_exact_loss_checks':checks,
            'configuration_exact_loss_domination_failures':failures,
            'rank_reverse_loss_order_failures':rank_reverse_failures,
            'serial_loss_NQD_failures':serial_failures}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    start=time.monotonic()
    out={'cases':[case(4,F(1,4),F(1)),case(6,F(1,2),F(1))],
         'path_concavity_obstruction':path_concavity_obstruction(),
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'elapsed_seconds':time.monotonic()-start,
         'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
    if args.output:
        args.output.write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'complete_reverse_order_failures':[c['reverse_order_failures'] for c in out['cases']],
                          'path':out['path_concavity_obstruction'],
                          'elapsed_seconds':out['elapsed_seconds'],
                          'peak_rss_mib':out['peak_rss_mib']},indent=2))
    else:
        print(json.dumps(out,indent=2))
