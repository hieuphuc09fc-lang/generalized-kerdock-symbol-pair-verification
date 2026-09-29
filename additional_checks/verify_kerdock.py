import sys, json, hashlib
from collections import Counter
from pathlib import Path
import numpy as np

def field(q,m,coeff):
    ell=q.bit_length()-1; n=q**m
    def base(a,b):
        v=0
        while b:
            if b&1:v^=a
            b>>=1;a<<=1
            if a&q:a^={2:3,4:7,8:11}[q]
        return v
    def mul(a,b):
        c=[0]*(2*m-1)
        for i in range(m):
            for j in range(m):c[i+j]^=base((a>>(ell*i))&(q-1),(b>>(ell*j))&(q-1))
        for i in range(2*m-2,m-1,-1):
            for j in range(m):c[i-m+j]^=base(c[i],coeff[j])
        return sum(c[i]<<(ell*i) for i in range(m))
    powers=[];a=1
    for i in range(n-1):powers.append(a);a=mul(a,q)
    assert a==1 and len(set(powers))==n-1
    ex=np.array(powers+[1],dtype=np.int32)
    log=np.full(n,-1,dtype=np.int32)
    log[powers]=np.arange(n-1)
    def mt(a,b):
        a=np.asarray(a);b=np.asarray(b)
        return np.where((a==0)|(b==0),0,ex[(log[a]+log[b])%(n-1)])
    x=np.arange(n);tr=np.zeros(n,dtype=np.int32);Q=tr.copy();y=x.copy()
    for i in range(m):
        Q^=mt(tr,y);tr^=y
        for j in range(ell):y=mt(y,y)
    assert np.max(tr)<q and np.max(Q)<q
    return powers,mt,tr,Q

def spectrum(q,m,coeff):
    powers,mt,tr,Q=field(q,m,coeff);n=q**m;p=np.array(powers)
    sqrt=np.array([next(b for b in range(q) if int(mt(b,b))==a) for a in range(q)])
    omega=q-1;T=tr[p];quad=Q[p];weights=Counter();Ns=set()
    # Nonunit xi: gamma denotes xi_1, and xi_0=0.
    for gamma in range(n):
        L=tr[mt(gamma,np.array(powers+[0]))]
        for c0 in range(q):
            for c1 in range(q):
                u1=L^c1;F=u1^int(mt(omega,c0));G=u1
                N=int(np.count_nonzero((F==0)&(np.roll(G,-1)==0)))
                nu=int(np.count_nonzero(u1==0)) if c0==0 else 0
                weights[q*n-(q-1)*nu-N]+=1
    # xi_0 != 0. Normalize y=xi_0*x and gamma=xi_1/xi_0.
    for gamma in range(n):
        L=tr[mt(gamma,p)]
        for c0 in range(q):
            u0=T^c0
            core=sqrt[quad]^L^sqrt[mt(c0,T)]
            for c1 in range(q):
                u1=core^c1;F=u1^mt(omega,u0);G=u1
                cycle=int(np.count_nonzero((F==0)&(np.roll(G,-1)==0)))
                # Each entry corresponds to xi_0=p[j], last y=p[j-1].
                tail=(np.roll(F,1)==0);first=(G==0)
                N=cycle-(tail&first).astype(int)
                N+=(tail&(c1==0)).astype(int)
                N+=((c1^int(mt(omega,c0))==0)&first).astype(int)
                nu=int(np.count_nonzero((u0==0)&(u1==0)))+int(c0==c1==0)
                vals,cnts=np.unique(q*n-(q-1)*nu-N,return_counts=True)
                for v,cnt in zip(vals,cnts):weights[int(v)]+=int(cnt)
                Ns.update(map(int,N))
    assert sum(weights.values())==q**(2*m+2)
    result={'q':q,'m':m,'size':sum(weights.values()),'weights':sorted(weights),'N_unit':sorted(Ns),'frequencies':dict(sorted(weights.items()))}
    if q==2 and m in (3,5):
        # Independent direct construction and exhaustive pair distances.
        x=np.array(powers+[0]);words=[];params=[]
        for a in range(n):
            z=mt(a,x);L0=tr[z]
            for b in range(n):
                L1=tr[mt(b,x)]
                for c0 in range(2):
                    for c1 in range(2):
                        v0=L0^c0;v1=Q[z]^L1^(c0*L0)^c1
                        bits=np.column_stack((v1,v1^v0)).reshape(-1)
                        word=sum(int(v)<<i for i,v in enumerate(bits))
                        words.append(word);params.append([a,b,c0,c1])
        assert len(set(words))==len(words)
        length=2*n;mask=(1<<length)-1
        def wt(v):return (v|((v>>1)|((v&1)<<(length-1)))).bit_count()
        assert Counter(map(wt,words))==weights
        best=length;witness=None
        for i,a in enumerate(words):
            for j in range(i):
                d=wt(a^words[j])
                if d<best:best=d;witness=[params[j],params[i]]
        result['distance']=best;result['witness']=witness
        result['word_checksum']=hashlib.sha256(b''.join(v.to_bytes((length+7)//8,'little') for v in words)).hexdigest()
    return result

if __name__=='__main__':
    cases=[(2,3,[1,1,0]),(2,5,[1,0,1,0,0]),(2,7,[1,1,0,0,0,0,0]),(2,9,[1,0,0,0,1,0,0,0,0]),(2,11,[1,0,1,0,0,0,0,0,0,0,0]),(4,3,[2,1,1]),(4,5,[2,1,0,1,0]),(8,3,[2,7,2])]
    out=[]
    for q,m,c in cases:
        r=spectrum(q,m,c);out.append(r)
        print(json.dumps(r),flush=True)
        Path(__file__).with_name('verification-results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
