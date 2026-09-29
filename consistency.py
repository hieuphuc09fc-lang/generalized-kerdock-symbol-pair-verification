import json
from pathlib import Path

data=json.loads(Path('verification/results.json').read_text())
data['q8']=json.loads(Path('verification/q8_results.json').read_text())
out={}
for key,d in data.items():
    if not key.startswith('weights') and key!='q8':continue
    q,m=d['q'],d['m'];la=(m-1)//2;n=q**(m+1);A=q**(m-2);S=q**(la-1)
    if q==2:omega=set(range(-2,3))
    else:
        omega={(q-1)*r+s for r in range(-q,q+1) for s in range(-q,q+1) if abs(r)+abs(s)<=q}
        omega|={(q-1)*r+s+t for r in range(1-q,q) for s in range(1-q,q) for t in range(-2,3) if abs(r)+abs(s)<=q-1}
    C={A+S*k+b for k in omega for b in [-1,0,1] if A+S*k+b>=0}
    candidates={0,n,n-(q-1)*q**(m-1)-A,n-(q-1)*q**(m-1)-A+1,n-A,n-A-1}
    if q>2:candidates.add(n-A+1)
    for delta in [-1,1]:
        for v in C:
            candidates.add(n-(q-1)*(A-delta*S)-v)
            candidates.add(n-(q-1)*(A-delta*(q-1)*S)-v)
    assert set(d['weights'])<=candidates
    expected={0:1,n:q-1,(q-1)*n//q:q*(n-1)}
    # Check integer frequency sums for the cited Hamming distribution, merging duplicates.
    for delta in [-1,1]:
        w=(q-1)*(n+delta*q**(la+1))//q
        freq=q*q*(q**m-1)*(q**(m-1)-delta*q**la)//2
        expected[w]=expected.get(w,0)+freq
        w=(q-1)*n//q+delta*q**la
        freq=(q-1)*q*q*(q**m-1)*(q**(m-1)+delta*q**la)//2
        expected[w]=expected.get(w,0)+freq
    assert sum(expected.values())==n*n
    out[key]={'candidate_count':len(candidates),'candidate_set':sorted(candidates),'unattained':sorted(candidates-set(d['weights'])),'N_delta_in_candidate':set(d['N_delta'])<=C,'candidate_N_count':len(C),'frequency_total':sum(expected.values())}
Path('verification/consistency_results.json').write_text(json.dumps(out,indent=2))
print({k:(v['candidate_count'],v['N_delta_in_candidate']) for k,v in out.items()})
