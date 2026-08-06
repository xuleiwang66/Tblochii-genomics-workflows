#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--inputs',nargs='+',required=True); p.add_argument('--out',required=True); a=p.parse_args()
rows=[]; fields=None
for f in a.inputs:
    with open(f,encoding='utf-8-sig') as h:
        r=csv.DictReader(h,delimiter='\t')
        if fields is None: fields=r.fieldnames
        elif fields!=r.fieldnames: raise SystemExit(f'Header mismatch: {f}')
        rows.extend(r)
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
with open(a.out,'w',encoding='utf-8',newline='') as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
