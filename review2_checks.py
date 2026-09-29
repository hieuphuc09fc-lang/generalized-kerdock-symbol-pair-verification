from pathlib import Path
import json,math
from field_audit import field,spectrum
out=Path(__file__).resolve().parent
scan=[]
p,mt,tr,Q=field(4,3,[2,1,1])
for k in range(1,63):
    if math.gcd(k,63)!=1:continue
    omega=p[(21*k)%63];omega2=int(mt(omega,omega))
    for name,end in [('A',omega2),('B',omega)]:
        r=spectrum(4,3,[2,1,1],k,end)
        scan.append(dict(k=k,order=name,omega=omega,endpoint=end,weights=r['weights'],N_unit=r['N_unit'],count=len(r['weights']),minimum=min(w for w in r['weights'] if w)))
(out/'q4_m3_all_orders.json').write_text(json.dumps(scan,indent=2),encoding='utf-8')
print('base theta:',[(r['order'],r['count'],r['N_unit']) for r in scan if r['k']==1],flush=True)
from collections import Counter
print('counts',Counter(r['count'] for r in scan),'Nsets',Counter(tuple(r['N_unit']) for r in scan),flush=True)
print('missing both',[k for k in {r['k'] for r in scan} if all(248 not in r['weights'] for r in scan if r['k']==k)],flush=True)
for q,m,co in [(4,5,[2,1,0,1,0]),(8,3,[2,7,2])]:
    r=spectrum(q,m,co)
    O0={(q-1)*a+b for a in range(-q,q+1) for b in range(-q,q+1) if abs(a)+abs(b)<=q}
    O1={(q-1)*a+b+t for a in range(-q,q+1) for b in range(-q,q+1) if abs(a)+abs(b)<=q-1 for t in range(-2,3)}
    lam=(m-1)//2
    N={q**(m-2)+q**(lam-1)*k+b for k in O0|O1 for b in [-1,0,1] if q**(m-2)+q**(lam-1)*k+b>=0}
    length=q**(m+1)
    C={0,length}|{length-(q-1)*q**(m-1)-x for x in [q**(m-2)-1,q**(m-2)]}|{length-x for x in [0,q**(m-2)-1,q**(m-2),q**(m-2)+1]}
    for nu in [q**(m-2)+d*t*q**(lam-1) for d in [-1,1] for t in [1,q-1]]:
        C|={length-(q-1)*nu-x for x in N}
    r['candidate_weights']=sorted(C);r['candidate_count']=len(C);r['actual_count']=len(r['weights']);r['unattained_candidates']=sorted(C-set(r['weights']))
    assert set(r['weights'])<=C
    (out/f'q{q}_m{m}_rechecked.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
    print(q,m,'actual',r['actual_count'],'candidate',r['candidate_count'],'minimum',min(w for w in r['weights'] if w),flush=True)
