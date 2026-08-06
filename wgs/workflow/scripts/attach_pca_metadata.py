#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--eigenvec',required=True); p.add_argument('--metadata',required=True); p.add_argument('--out',required=True); a=p.parse_args()
meta={}
with open(a.metadata,encoding='utf-8-sig') as h:
    for r in csv.DictReader(h,delimiter='\t'): meta[r['sample_id']]=r
lines=[x.split() for x in open(a.eigenvec) if x.strip()]
header=[x.lstrip('#') for x in lines[0]]
if not any(x.startswith('PC') for x in header):
    n=len(lines[0])-2; header=['FID','IID']+[f'PC{i}' for i in range(1,n+1)]; data=lines
else: data=lines[1:]
fields=['sample_id','source_cohort']+[x for x in header if x.startswith('PC')]
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
with open(a.out,'w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader()
    for vals in data:
        r=dict(zip(header,vals)); sid=r.get('IID',r.get('FID')); row={'sample_id':sid,'source_cohort':meta[sid]['source_cohort']}; row.update({k:r[k] for k in fields if k.startswith('PC')}); w.writerow(row)
