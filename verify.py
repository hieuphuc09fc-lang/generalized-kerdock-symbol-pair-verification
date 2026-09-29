"""Independent finite-field audit. Python 3 + numpy; no Galois-ring package.
Field elements are base-q coefficient integers. Constants are 0,...,q-1.
The exhaustive unit-parameter sweep uses z=xi0*x; it covers every xi0 != 0,
gamma=xi1/xi0,c0,c1, and corrects the broken multiplicative cycle exactly.
"""
import json, time
from pathlib import Path
import numpy as np

class Field:
    def __init__(self,q,m,coeff):
        self.q,self.m,self.n=q,m,q**m
        self.bits=(q.bit_length()-1)*m
        small=np.zeros((q,q),dtype=np.int32)
        mod={2:3,4:7,8:11}[q]
        for a in range(q):
            for b in range(q):
                x,y,v=a,b,0
                while y:
                    if y&1:v^=x
                    y>>=1;x<<=1
                    if x&q:x^=mod
                small[a,b]=v
        def mul(a,b):
            v=0
            for _ in range(m):
                digit=b%q;b//=q
                for j in range(m):v^=int(small[(a//q**j)%q,digit])*q**j
                top=a//q**(m-1);a=(a*q)%self.n
                for j,c in enumerate(coeff):a^=int(small[top,c])*q**j
            return v
        vals=np.arange(self.n,dtype=np.int32)
        self.mul=np.zeros((self.n,self.n),dtype=np.int32)
        for k in range(self.bits):
            row=np.array([mul(1<<k,int(b)) for b in vals])
            self.mul^=np.where((vals[:,None]&(1<<k))!=0,row[None,:],0)
        self.sq=self.mul[vals,vals]
        f=vals.copy(); conj=[]
        for _ in range(m):
            conj.append(f)
            for _ in range(q.bit_length()-1):f=self.sq[f]
        assert np.array_equal(f,vals),'not a field with the stated degree'
        self.tr=np.bitwise_xor.reduce(conj,axis=0)
        assert max(self.tr)<q
        self.Q=np.zeros(self.n,dtype=np.int32)
        for i in range(m):
            for j in range(i):self.Q^=self.mul[conj[i],conj[j]]
        assert max(self.Q)<q
        self.sqrt=np.arange(q,dtype=np.int32)
        for _ in range(q.bit_length()-2):self.sqrt=self.sq[self.sqrt]
        # sqrt(t)=t^(q/2); q=2 has exponent 1.
        powers=[1]
        while True:
            z=int(self.mul[powers[-1],q])
            if z==1:break
            if z in powers:raise AssertionError('cycle does not return to 1')
            powers.append(z)
        assert len(powers)==self.n-1,('nonprimitive polynomial',q,m,len(powers))
        self.order=np.array(powers+[0],dtype=np.int32)
        self.prev=self.mul[powers[-1]]
        self.next=self.mul[q]

def audit_weights(f,omega):
    q,m,N=f.q,f.m,f.n; weights={0,q*N}; nset=set(); hist={}; example={}
    for gamma in range(N):
        L=f.tr[f.mul[gamma]]
        for c0 in range(q):
            b0=f.tr^c0
            b1base=f.sqrt[f.Q]^L^f.sqrt[f.mul[c0,f.tr]]
            for c1 in range(q):
                b1=b1base^c1
                A=b1^f.mul[omega,b0];G=b1
                az=A==0;gz=G==0
                nu0=int(np.count_nonzero((b0==0)&gz))
                full=int(np.count_nonzero(az&gz[f.next]))
                p=az[f.prev[1:]];s=gz[1:]
                correction=-(p&s).astype(int)+(p&(c1==0)).astype(int)+((c1^int(f.mul[omega,c0])==0)&s).astype(int)
                counts=full-int(c0==0 and c1==0)+correction
                unique,cnts=np.unique(counts,return_counts=True)
                for v,cnt in zip(unique,cnts):
                    v=int(v);w=q*N-(q-1)*nu0-v
                    nset.add(v);weights.add(w);hist[w]=hist.get(w,0)+int(cnt)
                    example.setdefault(w,[int(np.nonzero(counts==v)[0][0])+1,gamma,c0,c1])
    # xi0=0, xi1 arbitrary; direct coordinate-level block count.
    for xi1 in range(N):
        L=f.tr[f.mul[xi1,f.order]]
        for c0 in range(q):
            for c1 in range(q):
                b1=L^c1
                nu0=int(np.count_nonzero(b1==0)) if c0==0 else 0
                end=b1^int(f.mul[omega,c0])
                v=int(np.count_nonzero((end==0)&(np.roll(b1,-1)==0)))
                w=q*N-(q-1)*nu0-v
                weights.add(w);hist[w]=hist.get(w,0)+1
    assert sum(hist.values())==(q*N)**2
    return dict(q=q,m=m,omega=omega,N_delta=sorted(nset),weights=sorted(weights),
                min_weight=min(weights-{0}),weight_frequencies=hist,unit_witnesses=example,
                covered_codewords=sum(hist.values()))

def binary_words(f,theta=None):
    N=f.n
    order=f.order
    if theta is not None:
        p=[1]
        for _ in range(N-2):p.append(int(f.mul[p[-1],theta]))
        assert len(set(p))==N-1
        order=np.array(p+[0])
    words=[];params=[]
    for a in range(N):
        z=f.mul[a,order];L=f.tr[z];Q=f.Q[z]
        for b in range(N):
            v=Q^f.tr[f.mul[b,order]]
            for c0 in range(2):
                for c1 in range(2):
                    first=v^c1^(c0*L);second=first^L^c0
                    bits=np.empty(2*N,dtype=np.uint8);bits[::2]=first;bits[1::2]=second
                    words.append(int.from_bytes(np.packbits(bits,bitorder='little').tobytes(),'little'))
                    params.append([a,b,c0,c1])
    assert len(set(words))==4*N*N
    return words,params

def binary_distance(f,theta=None):
    words,params=binary_words(f,theta);n=2*f.n;mask=(1<<n)-1;minimum=n;wit=None
    for i,a in enumerate(words):
        for j in range(i):
            v=a^words[j];rot=((v<<1)&mask)|(v>>(n-1))
            d=(v|rot).bit_count()
            if d<minimum:minimum=d;wit=[params[i],params[j],format(a,f'0{n}b')[::-1],format(words[j],f'0{n}b')[::-1]]
    return dict(q=2,m=f.m,theta=int(theta or 2),min_distance=minimum,witness=wit,pairs=len(words)*(len(words)-1)//2)

if __name__=='__main__':
    start=time.time();results={}
    for q,m,coeff in [(2,3,[1,1,0]),(2,5,[1,0,1,0,0]),(2,7,[1,1,0,0,0,0,0]),(2,9,[1,0,0,0,1,0,0,0,0]),(2,11,[1,0,1,0,0,0,0,0,0,0,0]),(4,3,[2,1,1]),(4,5,[2,1,0,1,0])]:
        f=Field(q,m,coeff)
        for om in ([1] if q==2 else [2,3]):
            r=audit_weights(f,om);results[f'weights_{q}_{m}_{om}']=r
            print(q,m,om,'min',r['min_weight'],'#weights',len(r['weights']),'N',r['N_delta'],flush=True)
        if q==2 and m in [3,5]:
            r=binary_distance(f);results[f'distance_2_{m}']=r;print(r,flush=True)
        Path('verification/results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print('elapsed',time.time()-start,flush=True)
