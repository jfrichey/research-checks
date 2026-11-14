"""Generate and verify the exact three-site occupation numerator certificate.

Requires SymPy (audited with 1.14.0). No floating point, parameter sampling, or output files. The saved certificate is checked against
the freshly generated coefficients. See response-renewal.md Section 25.
The fast physical identity audit is verify_three_site_renewal.py.
AI-written research code; no independent user-review or novelty claim.
"""
from collections import defaultdict
from pathlib import Path
import gzip
import hashlib
import json
from math import comb, factorial
import time
import sympy as s

started=time.monotonic()
def say(*args):
    print(round(time.monotonic()-started,2),*args,flush=True)

def build_numerator():
    a,b,c,d,e,f = edges = s.symbols('a b c d e f')
    r0,r1,r2 = r = s.symbols('r0 r1 r2')
    gens = edges+r
    K = s.Matrix([[0,a,b],[c,0,d],[e,f,0]])
    k = [1-a-b,1-c-d,1-e-f]
    ds = [1-d*f,1-b*e,1-a*c]
    detK = s.expand((s.eye(3)-K).det())
    Z = [1+ds[i]*(r[j]+r[l]+r[j]*r[l])
         for i in range(3) for j,l in [[z for z in range(3) if z != i]]]
    N = s.zeros(3)
    for i in range(3):
        N[i,i] = Z[i]
        for j in range(3):
            if i != j:
                l = 3-i-j
                N[i,j] = -r[i]*(K[i,j]*(1+ds[i]*r[l])+K[i,l]*K[l,j])
    detN = s.Poly(N.det(), *gens)
    adjN = N.adjugate()
    D1 = s.expand((s.eye(3)+s.diag(*r)*(s.eye(3)-K)).det())
    say('detN monomials',len(detN.terms()),'sleep degree',s.Poly(detN.as_expr(),*r).total_degree())
    
    # Weighted two-site no-hit law on {1,2}, with 0 killed.
    lam = 1+r1+r2+ds[0]*r1*r2
    U1 = (1+r2)*(1+r1*k[1])+r1*d*(1+r2*k[2])
    U2 = (1+r1)*(1+r2*k[2])+r2*f*(1+r1*k[1])
    B1 = r2*((1+r1)*k[2]+f*k[1])
    B2 = r1*((1+r2)*k[1]+d*k[2])
    Alpha = s.expand(lam*(1-r1*r2*k[1]*k[2])+B1*U1+B2*U2)
    L1,L2 = k[1]+d*k[2],k[2]+f*k[1]
    T1 = s.expand(r1*(ds[0]*Alpha-L1*U2*Z[0]))
    T2 = s.expand(r2*(ds[0]*Alpha-L2*U1*Z[0]))
    F2 = s.expand(Z[0]*(a*U1+b*U2)-b*(Alpha+T1)-a*(k[1]*U1*Z[0]+d*Alpha))
    F1 = s.expand(Z[0]*(a*U1+b*U2)-a*(Alpha+T2)-b*(k[2]*U2*Z[0]+f*Alpha))
    Fstar = s.expand(ds[0]*Z[0]*(a*U1+b*U2)-Alpha*(a*L1+b*L2)
                     -ds[0]*(a*T2*k[1]+b*T1*k[2]))
    say('local numerator monomials', {name:len(s.Poly(poly,*gens).terms())
                                      for name,poly in [('F1',F1),('F2',F2),('Fstar',Fstar)]})
    
    P = s.prod(ds)
    Zbar = s.expand(P+detK*(P*(r0*r1+r0*r2+r1*r2+r0*r1*r2)
                             +sum(r[i]*s.prod(ds[j] for j in range(3) if j != i) for i in range(3))))
    bbar=[]
    for i in range(3):
        j,l=[z for z in range(3) if z != i]
        bbar.append(s.expand(r[i]*(detK*s.prod(ds[z] for z in range(3) if z != i)
                      +P*ds[i]*(k[i]*(r[j]+r[l]+r[j]*r[l])
                         +r[l]*K[i,j]*k[j]+r[j]*K[i,l]*k[l]))))
    source=s.Matrix(bbar)-Zbar*s.eye(3)[:,1]
    Frow=s.Matrix([[0,F1,F2]])
    say('building reduced numerator by sparse polynomial arithmetic')
    H=s.Poly(ds[1]*ds[2],*gens)*s.Poly(Fstar,*gens)*detN
    say('initial H monomials',len(H.terms()))
    for i,Fi in [(1,F1),(2,F2)]:
        for j in range(3):
            H+=s.Poly(Fi,*gens)*s.Poly(adjN[i,j],*gens)*s.Poly(source[j],*gens)
            say('H partial',i,j,'monomials',len(H.terms()))
    say('H monomials',len(H.terms()),'sleep monomials',len({m[-3:] for m in H.monoms()}))
    cancelled=[]
    for name,fac in [('d0',ds[0]),('d1',ds[1]),('d2',ds[2]),('Z0',Z[0]),('D1',D1)]:
        divisor=s.Poly(fac,*gens)
        multiplicity=0
        while True:
            quo,rem=H.div(divisor)
            if not rem.is_zero:break
            H=quo;multiplicity+=1
        if multiplicity:cancelled.append((name,multiplicity))
    say('cancelled factors',cancelled,'remaining monomials',len(H.terms()))
    data=dict(gens=gens,K=K,k=k,r=r,ds=ds,detK=detK,Z=Z,N=N,detN=detN,D1=D1,
              lam=lam,U1=U1,U2=U2,Alpha=Alpha,T1=T1,T2=T2,F1=F1,F2=F2,Fstar=Fstar,
              P=P,Zbar=Zbar,bbar=bbar,H=H,cancelled=cancelled)
    return data

def positive_rate_certificate(data):
    H=data['H']
    blocks=defaultdict(dict)
    for mon,coeff in H.terms():blocks[mon[-3:]][mon[:6]]=int(coeff)
    assert set(blocks)=={(i,j,k) for i in range(3) for j in range(4) for k in range(5)}
    say('sleep blocks',len(blocks),'sleep degrees',[max(m[i] for m in blocks) for i in range(3)])
    power_cache={}
    def trinomial(n):
        if n not in power_cache:
            power_cache[n]=[(i,j,factorial(n)//(factorial(i)*factorial(j)*factorial(n-i-j)))
                            for i in range(n+1) for j in range(n-i+1)]
        return power_cache[n]
    def transform(poly,row,degree):
        result=defaultdict(int)
        a,b=2*row,2*row+1
        for mon,c in poly.items():
            for i,j,weight in trinomial(degree-mon[a]-mon[b]):
                key=list(mon);key[a]+=i;key[b]+=j
                result[tuple(key)]+=c*weight
        return {m:c for m,c in result.items() if c}
    cert={}
    negative=[]
    for index,(sleep,original) in enumerate(sorted(blocks.items())):
        degrees=tuple(max(m[2*i]+m[2*i+1] for m in original) for i in range(3))
        result=original
        for row in range(3):result=transform(result,row,degrees[row])
        neg={m:c for m,c in result.items() if c<0}
        cert[sleep]=(degrees,result)
        if neg:negative.append((sleep,len(neg),next(iter(neg.items()))))
        say('block',sleep,'row degrees',degrees,'monomials',len(result),'negative',len(neg))
    say("positive rate certificate",sum(len(p) for _,p in cert.values()),"coefficients; negative blocks",len(negative))
    return dict(certificate=cert,negative=negative,cancelled=data["cancelled"])

def independent_audit(data,cert):
    gens=data['gens']; r=data['r']; K=data['K']; k=data['k']; ds=data['ds']
    poly=lambda q:s.Poly(q,*gens)
    N=data['N']; adj=N.adjugate()
    adjN=[[poly(adj[i,j]) for j in range(3)] for i in range(3)]
    one=s.eye(3)+s.diag(*r)*(s.eye(3)-K)
    D1=poly(one.det()); onecol=[poly(one.adjugate()[i,0]) for i in range(3)]
    E21=s.zeros(3)
    for i in range(3):
        j,l=[z for z in range(3) if z!=i]
        E21[i,i]=r[j]*r[l]*(K[j,i]*k[l]+K[l,i]*k[j]+K[j,i]*K[l,i]*k[i])
        E21[i,j]=-r[l]*((1+r[j])*k[l]+K[l,j]*k[j])
        E21[i,l]=-r[j]*((1+r[l])*k[j]+K[j,l]*k[l])
    rhs=[(D1 if i!=0 else poly(0))-sum((poly(E21[i,j])*onecol[j] for j in range(3)),poly(0))
         for i in range(3)]
    E31=[]
    for i in range(3):
        j,l=[z for z in range(3) if z!=i]
        E31.append(poly(r[j]*r[l]*(r[i]*k[j]*k[l]
                     +(k[j]+K[j,i]*k[i])*(k[l]+K[l,i]*k[i]))))
    numerator=poly(data['P'])*data['detN']*(D1-sum((E31[i]*onecol[i] for i in range(3)),poly(0)))
    source=[poly(data['bbar'][i]-(data['Zbar'] if i==1 else 0)) for i in range(3)]
    say('complete-block independent numerator started')
    for i in range(3):
        for j in range(3):numerator+=rhs[i]*adjN[i][j]*source[j]
        say('complete-block row',i,'terms',len(numerator.terms()))
    expected=poly(r[0]*ds[1])*data['H']
    assert numerator==expected
    say('EXACT polynomial identity verified; numerator terms',len(numerator.terms()))
    
    assert not cert['negative']
    assert len(cert['certificate'])==60
    total=0
    for sleep,(degrees,coefficients) in cert['certificate'].items():
        assert all(isinstance(c,int) and c>0 for c in coefficients.values())
        total+=len(coefficients)
    assert total==330334
    say('positive certificate coefficients',total)
    
    # Independently invert each row homogenization using t=1+a+b:
    # C(a,b)=t^m P(a/t,b/t). Hence P(A,B)=(1-A-B)^m
    # C(A/(1-A-B),B/(1-A-B)). This checks exact integer coefficients.
    from collections import defaultdict
    from math import factorial
    blocks=defaultdict(dict)
    for mon,coeff in data['H'].terms():blocks[mon[-3:]][mon[:6]]=int(coeff)
    cache={}
    def trinomial(n):
        if n not in cache:
            cache[n]=[(i,j,(-1)**(i+j)*factorial(n)//(factorial(i)*factorial(j)*factorial(n-i-j)))
                      for i in range(n+1) for j in range(n-i+1)]
        return cache[n]
    def undo_row(original,row,degree):
        result=defaultdict(int)
        a,b=2*row,2*row+1
        for mon,c in original.items():
            assert mon[a]+mon[b]<=degree
            for i,j,weight in trinomial(degree-mon[a]-mon[b]):
                m=list(mon);m[a]+=i;m[b]+=j
                result[tuple(m)]+=c*weight
        return {m:c for m,c in result.items() if c}
    for sleep,(degrees,coefficients) in cert['certificate'].items():
        restored=coefficients
        for row in (2,1,0):restored=undo_row(restored,row,degrees[row])
        assert restored==blocks[sleep],sleep
    say('all 60 positive rate blocks inverted exactly to the original numerator')
    return dict(numerator_terms=len(data["H"].terms()),sleep_blocks=len(cert["certificate"]),positive_rate_coefficients=total)

def verify_saved_certificate(cert):
    path=Path(__file__).with_name("three-site-positive-certificate.json.gz")
    encoded=path.read_bytes()
    artifact=json.loads(gzip.decompress(encoded))
    assert artifact["target"]==0 and artifact["added"]==1
    assert artifact["rate_variables"]==["A01","A02","A10","A12","A20","A21"]
    assert artifact["reduced_numerator_monomials"]==23312
    assert artifact["positive_rate_monomials"]==330334
    blocks={tuple(block["sleep_exponents"]):
            (tuple(block["row_degrees"]),{tuple(m):c for m,c in block["terms"]})
            for block in artifact["blocks"]}
    assert len(blocks)==len(artifact["blocks"])==60
    assert blocks==cert["certificate"]
    digest=hashlib.sha256(encoded).hexdigest()
    say("saved coefficient certificate matches",digest)
    return digest


def structural_audit(data, cert):
    """Endpoint factors and the precise mutual-reachability support."""
    H=data['H']; gens=data['gens']; K=data['K']; k=data['k']; r=data['r']
    d0,d1,d2=data['ds']
    G01=K[0,1]+K[0,2]*K[2,1]
    G10=K[1,0]+K[1,2]*K[2,0]
    L1=k[1]+K[1,2]*k[2]
    L2=k[2]+K[2,0]*k[0]
    lam2=1+r[0]+r[1]+d2*r[0]*r[1]
    common=G01*G10*(1+r[1]*L1)*lam2*data['Z'][2]
    layers=[s.Poly.from_dict({m[:-1]:c for m,c in H.terms() if m[-1]==j},
                            gens[:-1]) for j in (0,4)]
    assert layers[0]==s.Poly(d2*common,*gens[:-1])
    assert layers[1]==s.Poly(d0*data['detK']*L2*common,*gens[:-1])
    motifs=((1,0,1,0,0,0),(1,0,0,1,1,0),
            (0,1,1,0,0,1),(0,1,0,1,1,1))
    def has_motif(mon):
        return any(all(a>=b for a,b in zip(mon,need)) for need in motifs)
    assert all(has_motif(m[:6]) for m in H.monoms())
    for odds,(_,terms) in cert['certificate'].items():
        weight=comb(2,odds[0])*comb(3,odds[1])*comb(4,odds[2])
        assert all(terms[motif]==weight for motif in motifs)
        assert all(has_motif(m) for m in terms)
    for path in (G01,G10):
        assert not H.div(s.Poly(path,*gens))[1].is_zero
    say('endpoint factors and all mutual-reachability motifs verified')
    return dict(endpoint_layers=2,mutual_reachability_motifs=len(motifs))


if __name__ == "__main__":
    data=build_numerator()
    certificate=positive_rate_certificate(data)
    verify_saved_certificate(certificate)
    print(independent_audit(data,certificate),flush=True)
    print(structural_audit(data,certificate),flush=True)
