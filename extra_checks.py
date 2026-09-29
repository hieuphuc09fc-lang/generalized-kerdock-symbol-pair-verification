import json, math
from pathlib import Path
import numpy as np
from verify import Field,binary_distance,binary_words

out={}
for m,coeff in [(3,[1,1,0]),(5,[1,0,1,0,0])]:
    f=Field(2,m,coeff);seen=set();rows=[]
    # Frobenius conjugate choices give the same ordered code: relabel trace parameters.
    for t in range(2,f.n):
        if t in seen:continue
        orbit=[];v=t
        while v not in orbit:orbit.append(v);v=int(f.sq[v])
        seen.update(orbit)
        r=binary_distance(f,t);words,_=binary_words(f,t);mask=(1<<(2*f.n))-1
        weights={ (x|((x<<1)&mask)|(x>>(2*f.n-1))).bit_count() for x in words }
        r['min_weight']=min(weights-{0});r['frobenius_orbit']=orbit;rows.append(r)
        print('primitive orbit',m,orbit,'weight',r['min_weight'],'distance',r['min_distance'],flush=True)
    out[str(m)]=rows

# Direct Z4 polynomial-ring construction, independent of the quadratic trace formula.
def ring_binary_check(m,coeff):
    N=2**m
    def add(a,b):return tuple((x+y)%4 for x,y in zip(a,b))
    def mul(a,b):
        v=[0]*(2*m-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):v[i+j]=(v[i+j]+x*y)%4
        for k in range(2*m-2,m-1,-1):
            for j,c in enumerate(coeff):v[k-m+j]=(v[k-m+j]-v[k]*c)%4
        return tuple(v[:m])
    def power(a,n):
        v=(1,)+(0,)*(m-1)
        while n:
            if n&1:v=mul(v,a)
            a=mul(a,a);n>>=1
        return v
    lifts=[power(tuple((i>>j)&1 for j in range(m)),N) for i in range(N)]
    theta=lifts[2];one=lifts[1];zero=lifts[0];p=[one]
    for _ in range(N-2):p.append(mul(p[-1],theta))
    p.append(zero)
    traces={}
    for a in lifts:
        s=zero
        for j in range(m):s=add(s,power(a,2**j))
        assert all(v==0 for v in s[1:]);traces[a]=s[0]
    words=set()
    for a in lifts:
        ta=[traces[mul(a,x)] for x in p]
        for b in lifts:
            tb=[traces[mul(b,x)] for x in p]
            for c in range(4):
                w=0
                for j,(x,y) in enumerate(zip(ta,tb)):
                    v=(x+2*y+c)%4
                    w|=(v//2)<<(2*j);w|=((v//2)^(v%2))<<(2*j+1)
                words.add(w)
    expected=set(binary_words(Field(2,m,coeff))[0]);assert words==expected
    return {'m':m,'ring':'Z4[X]/f(X), Teichmuller lifts via a^(2^m)','distinct_codewords':len(words),'matches_quadratic_construction':True}
out['direct_ring']=[ring_binary_check(3,[1,1,0]),ring_binary_check(5,[1,0,1,0,0])]
Path('verification/extra_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print('direct ring cross-checks passed',flush=True)
