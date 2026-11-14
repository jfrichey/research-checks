"""Exact audits of the creation-residual cone in response-renewal.md.

Standard library only. Finite LP certificates do not prove the proposed
arbitrary-graph positivity statement. AI-written research code.
"""
import importlib.util
from fractions import Fraction as F
from itertools import product
from pathlib import Path
spec=importlib.util.spec_from_file_location('rn',str(Path(__file__).with_name('verify_response_renewal.py')))
rn=importlib.util.module_from_spec(spec); spec.loader.exec_module(rn)
def simplex(A,b,c):
 m=len(b); n=len(c)
 tab=[[F(a) for a in A[i]]+[F(i==j) for j in range(m)]+[F(b[i])] for i in range(m)]
 tab.append([-F(a) for a in c]+[F(0)]*(m+1)); basis=list(range(n,n+m))
 for it in range(10000):
  j=next((j for j,a in enumerate(tab[-1][:-1]) if a<0),None)
  if j is None:
   sol=[F(0)]*n
   for i,k in enumerate(basis):
    if k<n: sol[k]=tab[i][-1]
   return tab[-1][-1],sol,it,[tab[-1][n+i] for i in range(m)]
  rows=[i for i in range(m) if tab[i][j]>0]
  assert rows,'unbounded'
  i=min(rows,key=lambda i:(tab[i][-1]/tab[i][j],basis[i])); pivot=tab[i][j]
  tab[i]=[a/pivot for a in tab[i]]
  for r in range(m+1):
   if r!=i and tab[r][j]:
    scale=tab[r][j]; tab[r]=[a-scale*b for a,b in zip(tab[r],tab[i])]
  basis[i]=j
 raise RuntimeError('iteration guard')
def probe(K,p,bound=2):
 n=len(p)
 S=[tuple(i for i,b in enumerate(bits) if b) for bits in product((0,1),repeat=n) if sum(bits)<=bound]
 _,T,Z,W=rn.positive_operator(K,p,bound)
 constraints=[]; labels=[]
 for x in range(n):
  for C in S:
   if x in C: continue
   labels.append((x,C))
   constraints.append([F(B==C)-(1-p[x])/p[x]*W.get((B,x),{}).get(C,F(0)) for B in S])
 A=[[-a for a in row] for row in constraints]+[[F(1)]*len(S)]; b=[F(0)]*len(constraints)+[F(1)]
 starts={}
 for B in S:
  for x in range(n):
   state=[-int(z in B) for z in range(n)]; state[x]=2 if state[x]==-1 else 1
   starts[B,x]=tuple(state)
 laws,_=rn.stopped(starts.values(),lambda s:rn.arw_step(K,p,s),lambda s:all(a<=0 for a in s))
 worst=(F(0),None,None,None); maxit=0; certificates=[]
 for x in range(n):
  for v in range(n):
   g=[sum(q for C,q in laws[starts[B,x]].items() if C[v]==-1)-F(v in B) for B in S]
   value,sol,it,dual=simplex(A,b,[-a for a in g]); maxit=max(maxit,it)
   assert all(sum(a*w for a,w in zip(row,sol))>=0 for row in constraints)
   assert -sum(a*w for a,w in zip(g,sol))==value
   assert all(w>=0 for w in dual)
   assert sum(w*rhs for w,rhs in zip(dual,b))==value
   assert all(sum(row[j]*w for row,w in zip(A,dual))>=-g[j] for j in range(len(S)))
   certificates.append((x,v,dual))
   if value>worst[0]: worst=(value,x,v,[(B,str(w)) for B,w in zip(S,sol) if w])
 return {'n':n,'bound':bound,'max_negative_response':str(worst[0]),'source':worst[1],'target':worst[2],'mu':worst[3],'max_pivots':maxit},(S,constraints,labels,certificates,W,starts,laws)

def intertwining_audit(K,p,bound=2):
 n=len(p); all_bound=min(n,bound+1)
 subsets=[tuple(i for i,b in enumerate(bits) if b) for bits in product((0,1),repeat=n) if sum(bits)<=all_bound]
 initial=[B for B in subsets if len(B)<=bound]
 _,_,_,W=rn.positive_operator(K,p,all_bound)
 def deletion(B,y,C):
  return F(B==C) if y not in B else -(1-p[y])/p[y]*W[B,y].get(C,F(0))
 starts={}
 for B in initial:
  for x in range(n):
   state=[-int(z in B) for z in range(n)]; state[x]=2 if state[x]==-1 else 1
   starts[B,x]=tuple(state)
 laws,_=rn.stopped(starts.values(),lambda s:rn.arw_step(K,p,s),lambda s:all(a<=0 for a in s))
 plus={(B,x):{tuple(z for z,a in enumerate(C) if a==-1):w for C,w in laws[s].items()} for (B,x),s in starts.items()}
 checks=0
 for y in range(n):
  columns=[C for C in subsets if y not in C]
  fresh=[C for C in initial if y not in C]
  kills,_=rn.stopped([starts[C,x] for C in fresh for x in range(n)],
                    lambda s:rn.arw_step(K,p,s),lambda s:s[y]>0 or all(a<=0 for a in s))
  killed={(C,x):{tuple(z for z,a in enumerate(D) if a==-1):w for D,w in kills[starts[C,x]].items() if D[y]<=0}
          for C in fresh for x in range(n)}
  for B in initial:
   for x in range(n):
    for D in columns:
     lhs=sum(w*deletion(E,y,D) for E,w in plus[B,x].items())
     rhs=sum(deletion(B,y,C)*killed[C,x].get(D,F(0)) for C in fresh)
     assert lhs==rhs,(B,x,y,D,lhs,rhs)
     checks+=1
 return {"intertwining_entries":checks,"verified":True}

def main():
 cases=[("archived_"+str(i),[[F(q) for q in row] for row in c["K"]],list(map(F,c["p"]))) for i,c in enumerate(rn.CASES)]
 cases.extend([
 ("forced_return",[[0,F(1,2),0,F(1,2)],[1,0,0,0],[F(1,2),0,0,F(1,2)]],[F(1,2)]*3),
 ("conditional_order",[[0,F(1,2),0,0,F(1,2)],[0,0,F(1,3),F(1,3),F(1,3)],[F(1,3),0,0,F(1,3),F(1,3)],[F(1,2),0,0,0,F(1,2)]],[F(2,3),F(1,2),F(9,10),F(9,10)]),
 ("uniform_routing",[[0,0,1,0,0],[F(1,2),0,0,0,F(1,2)],[0,0,0,F(1,2),F(1,2)],[F(1,2),F(1,2),0,0,0]],[F(1,16),F(1,2),F(1,2),F(1,16)]),
 ("creator_capacity",[[0,1,0,0],[F(1,2),0,0,F(1,2)],[0,1,0,0]],[F(1,2),F(99,100),F(1,2)])
 ])
 comparisons=0
 for label,K,p in cases:
  K=[[F(a) for a in row] for row in K]
  result,_=probe(K,p,min(2,len(p)))
  print(label,result,flush=True)
  comparisons+=len(p)**2
 K,p=cases[3][1:]
 print(intertwining_audit(K,p),flush=True)
 print({"named_cases":len(cases),"source_target_comparisons":comparisons,
        "scope":"finite cone certificates only; general sign remains open"},flush=True)
 if "--screen" in __import__("sys").argv:
  import random
  rng=random.Random(20261015)
  for t in range(10):
   n=3 if t<4 else 4; K=[]
   for i in range(n):
    row=[0 if j==i or rng.random()<.35 else rng.randint(1,5) for j in range(n)]+[rng.randint(1,3)]
    total=sum(row); K.append([F(a,total) for a in row])
   p=[rng.choice([F(1,20),F(1,3),F(1,2),F(4,5),F(19,20)]) for i in range(n)]
   result,_=probe(K,p)
   print("seeded_"+str(t),result,flush=True)
   if F(result["max_negative_response"])>0:
    print({"K":[list(map(str,row)) for row in K],"p":list(map(str,p))},flush=True)
    break

if __name__=="__main__":
 main()
