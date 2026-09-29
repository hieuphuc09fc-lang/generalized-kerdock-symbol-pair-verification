"""Directly recheck EVERY primitive element (no Frobenius reduction)."""
import json
from pathlib import Path
from verify import Field,binary_distance

result={}
for m,coeff in [(3,[1,1,0]),(5,[1,0,1,0,0])]:
    f=Field(2,m,coeff);rows=[]
    # Group orders 7 and 31 are prime, hence precisely 2,...,2^m-1 are primitive.
    for theta in range(2,f.n):
        r=binary_distance(f,theta)
        rows.append(r)
        print(m,theta,r['min_distance'],flush=True)
    result[str(m)]=rows
    Path('verification/all_primitive_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('ALL PRIMITIVE ELEMENTS CHECKED',flush=True)
