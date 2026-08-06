#!/usr/bin/env python3
import argparse,re,csv
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--logs',nargs='+',required=True); p.add_argument('--out',required=True); a=p.parse_args()
rows=[]
for f in a.logs:
    text=Path(f).read_text(errors='replace')
    m=re.findall(r'CV error \(K=(\d+)\):\s*([0-9.eE+-]+)',text)
    if not m: raise SystemExit(f'CV error not found: {f}')
    k,v=m[-1]; rows.append({'K':int(k),'CV_error':float(v)})
rows.sort(key=lambda r:r['K'])
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
with open(a.out,'w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=['K','CV_error'],delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
