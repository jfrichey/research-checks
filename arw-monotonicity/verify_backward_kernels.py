"""Exact audits of backward creator kernels and target-cycle certificates.

Standard library only. General sign claims remain open; exact finite
transport certificates are not a proof for arbitrary graphs. AI-written.
"""
from itertools import product
import verify_creation_cone as cc
F = cc.F

def feasible(A,b):
 m=len(b); n=len(A[0])
 rows=[[F(a) for a in row]+[F(i==j) for j in range(m)]+[F(b[i])] for i,row in enumerate(A)]
 for i in range(m):
  if rows[i][-1]<0:
   rows[i]=[-a for a in rows[i]];rows[i][n+i]=1
 obj=[F(0)]*n+[F(1)]*m+[F(0)]
 for row in rows:obj=[a-b for a,b in zip(obj,row)]
 tab=rows+[obj];basis=list(range(n,n+m))
 for it in range(20000):
  j=next((j for j,a in enumerate(tab[-1][:-1]) if a<0),None)
  if j is None:
   sol=[F(0)]*n
   for i,k in enumerate(basis):
    if k<n:sol[k]=tab[i][-1]
   if tab[-1][-1]==0:
    assert all(sum(a*x for a,x in zip(row,sol))==rhs for row,rhs in zip(A,b))
   return tab[-1][-1],sol,it
  opts=[i for i in range(m) if tab[i][j]>0];assert opts
  i=min(opts,key=lambda i:(tab[i][-1]/tab[i][j],basis[i]))
  pivot=tab[i][j];tab[i]=[a/pivot for a in tab[i]]
  for k in range(m+1):
   if k!=i and tab[k][j]:
    fac=tab[k][j];tab[k]=[a-fac*b for a,b in zip(tab[k],tab[i])]
  basis[i]=j
 raise RuntimeError('pivot guard')


def root_hitting_audit(case):
 K=[[F(q) for q in row] for row in case['K']];p=list(map(F,case['p']));n=len(p)
 S=[tuple(i for i,b in enumerate(bits) if b) for bits in case['subsets']]
 P=[[F(q) for q in row] for row in case['stabilization_matrix']]
 pi=dict(zip(S,P[S.index(tuple(range(n)))]))
 assert all(w>0 for w in pi.values())
 _,_,_,W=cc.rn.positive_operator(K,p,n)
 R={}
 for C in S:
  for z in range(n):
   if z in C:continue
   R[C,z]={B:(1-p[z])/p[z]*pi[B]*W[B,z].get(C,F(0))/pi[C] for B in S if z in B}
   assert sum(R[C,z].values())==1
   assert all(not w or len(B)>=len(C) for B,w in R[C,z].items())
   assert R[C,z][tuple(sorted(C+(z,)))]>0
 checks=0
 for A in S:
  k=len(A);trans=[C for C in S if len(C)<=k and C!=A]
  for choice in [min,max]:
   Q=[];b=[]
   for C in trans:
    z=choice(set(A)-set(C));row=R[C,z]
    Q.append([F(C==D)-row.get(D,F(0)) for D in trans]);b.append(row.get(A,F(0)))
   inv=cc.rn.inverse(Q)
   H={C:sum(w*q for w,q in zip(row,b)) for C,row in zip(trans,inv)};H[A]=F(1)
   scale=P[S.index(A)][S.index(A)]/pi[A]
   for C in S:
    assert P[S.index(A)][S.index(C)]==pi[C]*scale*H.get(C,F(0)),(A,C)
    checks+=1
 minor=0
 for A in S:
  for D in S:
   ratio=P[S.index(D)][S.index(A)]/P[S.index(A)][S.index(A)]
   for C in S:
    assert P[S.index(D)][S.index(C)]>=ratio*P[S.index(A)][S.index(C)]
    minor+=1
 return {'n':n,'reverse_rows':len(R),'root_hitting_entries':checks,'minorization_entries':minor,'verified':True}

def root_exit_audit(case):
    """Audit all boundary columns, not just hitting the original root support."""
    K = [[F(q) for q in row] for row in case["K"]]
    p = list(map(F, case["p"]))
    n = len(p)
    S = [tuple(i for i, b in enumerate(bits) if b) for bits in case["subsets"]]
    P = [[F(q) for q in row] for row in case["stabilization_matrix"]]
    pi = dict(zip(S, P[S.index(tuple(range(n)))]))
    _, _, _, W = cc.rn.positive_operator(K, p, n)
    reverse = {(C, z): {B: (1-p[z])/p[z] * pi[B] *
                         W[B, z].get(C, F(0)) / pi[C]
                         for B in S if z in B}
               for C in S for z in range(n) if z not in C}
    checked = 0
    for roots in S:
        boundary = [B for B in S if set(roots) <= set(B)]
        Q = [[P[S.index(D)][S.index(B)] for B in boundary] for D in boundary]
        gamma = cc.rn.multiply(cc.rn.inverse(Q), [P[S.index(D)] for D in boundary])
        assert all(q >= 0 for row in gamma for q in row)
        H = {(C, B): pi[B] * gamma[j][i] / pi[C]
             for i, C in enumerate(S) for j, B in enumerate(boundary)}
        for C in S:
            assert sum(H[C, B] for B in boundary) == 1
            for B in boundary:
                if C in boundary:
                    assert H[C, B] == F(C == B)
                if len(C) > len(B):
                    assert H[C, B] == 0
                for z in set(roots) - set(C):
                    assert H[C, B] == sum(q * H[D, B] for D, q in reverse[C, z].items())
                    checked += 1
    return dict(n=n, root_exit_harmonic_entries=checked, all_boundaries_verified=True)


def entrance_bias_obstruction():
    """The positive entrance components need not each favor even a common root."""
    n = 4
    K = [[F(0)] * 5 for _ in range(n)]
    for z in range(n):
        K[z][(z+1) % n] = F(31, 32)
        K[z][n] = F(1, 32)
    p = [F(1, 2), F(1, 2), F(1, 8), F(1, 8)]
    S = [tuple(i for i, b in enumerate(bits) if b)
         for bits in product((0, 1), repeat=n) if sum(bits) <= 3]
    starts = [tuple(int(z in A) for z in range(n)) for A in S]
    laws, _ = cc.rn.stopped(starts, lambda s: cc.rn.arw_step(K, p, s),
                           lambda s: all(a <= 0 for a in s))
    P = [[laws[st].get(tuple(-int(z in B) for z in range(n)), F(0))
          for B in S] for st in starts]
    roots, marked, target = (0, 2), (0, 2, 3), 0
    boundary = [B for B in S if set(roots) <= set(B)]
    Q = [[P[S.index(D)][S.index(B)] for B in boundary] for D in boundary]
    gamma = cc.rn.multiply(cc.rn.inverse(Q), [P[S.index(D)] for D in boundary])
    assert all(q >= 0 for row in gamma for q in row)
    masses = [sum(row) for row in gamma]
    baseline = sum(q for C, q in zip(S, P[S.index(roots)]) if target in C)
    gains = [sum(q for C, q in zip(S, row) if target in C) / mass - baseline
             for row, mass in zip(gamma, masses)]
    diff = gains[boundary.index(marked)]
    assert diff == -F(139677858771237304857855520692679096496067445789542055936,
                      20105157511516783257569470544841023362289343683299598529833)
    increments = {}
    for x in (1, 3):
        upper = tuple(sorted(roots + (x,)))
        weights = [P[S.index(upper)][S.index(B)] * mass
                   for B, mass in zip(boundary, masses)]
        assert sum(weights) == 1
        increment = sum(q for C, q in zip(S, P[S.index(upper)]) if target in C) - baseline
        assert increment == sum(w * g for w, g in zip(weights, gains))
        assert increment > 0
        increments[x] = str(increment)
    return dict(entrance_target_gain=str(diff), actual_increments_positive=sorted(increments))


def source_green_obstruction():
    """The one-extra boundary matrix need not be a killed-chain Green matrix."""
    supports = [[1, 3, 5], [2, 5], [0, 3, 5], [4, 5], [0, 5]]
    K = [[F(0)] * 6 for _ in range(5)]
    for z, row in enumerate(supports):
        for y in row:
            K[z][y] = F(1, len(row))
    p = [F(1, 2)] * 5
    boundary = [(3,), (0, 3), (1, 3), (2, 3), (3, 4)]
    starts = [tuple(int(z in A) for z in range(5)) for A in boundary]
    laws, _ = cc.rn.stopped(starts, lambda s: cc.rn.arw_step(K, p, s),
                           lambda s: all(a <= 0 for a in s))
    Q = [[laws[st].get(tuple(-int(z in B) for z in range(5)), F(0))
          for B in boundary] for st in starts]
    inverse = cc.rn.inverse(Q)
    value = inverse[boundary.index((1, 3))][boundary.index((3, 4))]
    assert value == F(1, 80636)
    # An independent root with p=1/2 multiplies every boundary entry by 1/2.
    embedded_inverse = cc.rn.inverse([[q/2 for q in row] for row in Q])
    assert embedded_inverse[2][4] == F(1, 40318)
    return dict(positive_off_diagonal_inverse=str(value),
                independent_second_root_value=str(embedded_inverse[2][4]))


def reward_endpoint(K,p):
 n=len(p);G=cc.rn.inverse([[F(i==j)-(1-p[i])*K[i][j] for j in range(n)] for i in range(n)])
 H={}
 for a in range(n):
  inv=cc.rn.inverse([[F(i==j)-(F(1) if i==a else 1-p[i])*K[i][j] for j in range(n)] for i in range(n)])
  H[a]=[[inv[x][y]*p[y] if y!=a else F(0) for y in range(n)] for x in range(n)]
 hit=[[G[a][y]/G[y][y] for y in range(n)] for a in range(n)]
 T=[[H[a][a][y]*hit[a][y] for y in range(n)] for a in range(n)]
 invT=cc.rn.inverse([[F(a==b)-T[a][b] for b in range(n)] for a in range(n)])
 credit={}
 for v in range(n):
  other=[a for a in range(n) if a!=v];kill=cc.rn.inverse([[F(a==b)-(1-p[a])*K[a][b] for b in other] for a in other])
  for z in other:
   N=[F(0)]*n
   for a in other:
    N[a]=kill[other.index(a)][other.index(z)]*(p[z]+(1-p[z])*K[z][n])
   reward=[H[a][a][v]*N[a] for a in range(n)]
   d=[sum(w*b for w,b in zip(row,reward)) for row in invT]
   for x in range(n):
    credit[x,v,z]=sum(H[v][x][y]*hit[v][y]*d[y] for y in range(n))/(G[v][v]*p[v])
 return credit

def transport_audit(K,p,bound=2):
 n=len(p);bound=min(bound,n)
 _,data=cc.probe(K,p,bound)
 S,L,labels,_,W,starts,laws=data
 g={(B,x,v):sum(q for C,q in laws[starts[B,x]].items() if C[v]==-1)-F(v in B)
    for B in S for x in range(n) for v in range(n)}
 credit=reward_endpoint(K,p)
 controlled=0;uniform_bad=[];credit_bad=[];certificate_entries=0
 for v in range(n):
  occ=[B for B in S if v in B]
  chosen=[(z,C) for z,C in labels if z!=v and v in C]
  M=[]
  for B in occ:
   row=[]
   for z,C in chosen:
    val=F(B==C)-(1-p[z])/p[z]*W.get((B,z),{}).get(C,F(0))
    val-=(1-p[v])/p[v]*(1-p[z])/p[z]*sum(
        w*W.get((D,z),{}).get(C,F(0)) for D,w in W[B,v].items())
    row.append(val)
   M.append(row)
  for x in range(n):
   r=[g[B,x,v]+(1-p[v])/p[v]*sum(w*g[C,x,v] for C,w in W[B,v].items()) for B in occ]
   A=[row+[F(i==j) for j in range(len(occ))] for i,row in enumerate(M)]
   value,solution,_=feasible(A,r)
   assert value==0,(x,v,"controlled target-cycle certificate failed",value)
   a=dict(zip(chosen,solution[:len(chosen)]))
   h={}
   for z,C in labels:
    if z==v:
     h[z,C]=g[C,x,v]+sum((1-p[y])/p[y]*sum(
         W.get((C,y),{}).get(D,F(0))*weight for (mark,D),weight in a.items() if mark==y)
         for y in C)
    else:h[z,C]=a.get((z,C),F(0))
   assert all(weight>=0 for weight in h.values())
   for i,B in enumerate(S):
    remainder=g[B,x,v]-sum(L[j][i]*h[lab] for j,lab in enumerate(labels))
    assert remainder>=0,(x,v,B,remainder)
    certificate_entries+=1
   controlled+=1
   if n>bound:
    aggregate=[[sum(M[i][j] for j,(z,D) in enumerate(chosen) if D==C) for C in occ]
               for i,B in enumerate(occ)]
    inverse=cc.rn.inverse(aggregate)
    assert all(a>=0 for row in inverse for a in row)
    uniform=[sum(a*b for a,b in zip(row,r)) for row in inverse]
    for C,val in zip(occ,uniform):
     if val<0:uniform_bad.append((x,v,C,str(val)))
   i0=occ.index((v,))
   assert sum(credit[x,v,z] for z in range(n) if z!=v)==r[i0]
   fixed={j:credit[x,v,z]/M[i0][j] for j,(z,C) in enumerate(chosen) if C==(v,)}
   rem=[r[i]-sum(M[i][j]*w for j,w in fixed.items()) for i in range(len(occ))]
   free=[j for j in range(len(chosen)) if j not in fixed]
   matrix=[[M[i][j] for j in free]+[F(i==j) for j in range(len(occ))]
           for i in range(len(occ))]
   val,_,_=feasible(matrix,rem)
   if val<0:credit_bad.append((x,v,str(val)))
 return {"n":n,"controlled_certificates":controlled,"checked_entries":certificate_entries,
         "uniform_failures":len(uniform_bad),"uniform_first":uniform_bad[:1],
         "last_vertex_credit_failures":credit_bad}

def strict_cone_audit(case):
 K=[[F(q) for q in row] for row in case["K"]];p=list(map(F,case["p"]))
 _,data=cc.probe(K,p,2);S,L,_,_,_,starts,laws=data
 mu={():F(5,21),(0,):F(8,21),(1,):F(8,21),(0,1):F(0)}
 assert all(sum(row[i]*mu[B] for i,B in enumerate(S))>=0 for row in L)
 oldS=[tuple(i for i,b in enumerate(bits) if b) for bits in case["subsets"]]
 E=[[F(q) for q in row] for row in case["inverse_matrix"]]
 coefficient=sum(mu[B]*E[oldS.index(B)][oldS.index(())] for B in S)
 assert coefficient==-F(1,7)
 P=[[F(q) for q in row] for row in case["stabilization_matrix"]]
 for v in range(2):
  rho=dict(zip(oldS,P[oldS.index((v,))]))
  scale=mu[(v,)]/rho[(v,)]
  assert all(mu[B]-scale*rho[B]>=0 for B in S)
  assert mu[(v,)]-scale*rho[(v,)]==0
 return {"strict_cone_example":True,"signed_empty_input_coefficient":str(coefficient)}

def forced_background_obstruction():
 K=[[F(0),F(1,2),F(0),F(0),F(1,2)],
    [F(0),F(0),F(1,3),F(1,3),F(1,3)],
    [F(1,3),F(0),F(0),F(1,3),F(1,3)],
    [F(1,2),F(0),F(0),F(0),F(1,2)]]
 p=[F(2,3),F(1,2),F(9,10),F(9,10)];v=1
 lower=(0,1,0,1);upper=(0,1,1,1)
 laws,_=cc.rn.stopped([lower,upper],lambda s:cc.rn.arw_step(K,p,s),lambda s:all(a<=0 for a in s))
 _,_,_,W=cc.rn.positive_operator(K,p,3)
 result={}
 for A in [lower,upper]:
  rho={};occupation=F(0)
  for final,w in laws[A].items():
   B=tuple(z for z,a in enumerate(final) if a==-1)
   if v in B:
    occupation+=w;rho[B]=w*(1+(1-p[v])/p[v]*sum(W[B,v].values()))
  assert sum(rho.values())==1
  result[A]=(occupation,sum(w for B,w in rho.items() if 3 in B))
 diff=result[upper][1]-result[lower][1]
 assert diff==-F(5116827833076768435777,211074938042774429635582400)
 assert result[upper][0]>result[lower][0]
 return {"last_arrival_weighted_order_difference":str(diff),
         "actual_target_increment_positive":True}

def cases():
 result=[("archived_"+str(i),[[F(q) for q in row] for row in c["K"]],list(map(F,c["p"])))
         for i,c in enumerate(cc.rn.CASES) if len(c["p"])>=2]
 result.extend([
 ("forced_return",[[0,F(1,2),0,F(1,2)],[1,0,0,0],[F(1,2),0,0,F(1,2)]],[F(1,2)]*3),
 ("conditional_order",[[0,F(1,2),0,0,F(1,2)],[0,0,F(1,3),F(1,3),F(1,3)],[F(1,3),0,0,F(1,3),F(1,3)],[F(1,2),0,0,0,F(1,2)]],[F(2,3),F(1,2),F(9,10),F(9,10)]),
 ("uniform_routing",[[0,0,1,0,0],[F(1,2),0,0,0,F(1,2)],[0,0,0,F(1,2),F(1,2)],[F(1,2),F(1,2),0,0,0]],[F(1,16),F(1,2),F(1,2),F(1,16)]),
 ("creator_capacity",[[0,1,0,0],[F(1,2),0,0,F(1,2)],[0,1,0,0]],[F(1,2),F(99,100),F(1,2)])
 ])
 return [(label,[[F(a) for a in row] for row in K],p) for label,K,p in result]

if __name__=="__main__":
 for case in cc.rn.CASES[1:]:
  print(root_hitting_audit(case),flush=True)
  print(root_exit_audit(case),flush=True)
 print(strict_cone_audit(cc.rn.CASES[1]),flush=True)
 for label,K,p in cases():
  print(label,transport_audit(K,p),flush=True)
 print(forced_background_obstruction(),flush=True)
 print(entrance_bias_obstruction(),flush=True)
 print(source_green_obstruction(),flush=True)
