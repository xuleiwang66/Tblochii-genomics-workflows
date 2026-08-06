#!/usr/bin/env python3
import argparse
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--raw',required=True); p.add_argument('--missingness',required=True); p.add_argument('--common',required=True); p.add_argument('--pruned',required=True); p.add_argument('--out',required=True); a=p.parse_args()
def counts(prefix):
    psam=Path(prefix+'.psam'); pvar=Path(prefix+'.pvar')
    n=sum(1 for x in psam.open() if x.strip() and not x.startswith('#'))
    m=sum(1 for x in pvar.open() if x.strip() and not x.startswith('#'))
    return n,m
rows=[]
for name,prefix in [('raw',a.raw),('after_mind_geno',a.missingness),('common_qc',a.common),('ld_pruned',a.pruned)]:
    n,m=counts(prefix); rows.append((name,n,m))
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
with open(a.out,'w') as h:
    h.write('dataset\tsamples\tvariants\n'); [h.write(f'{x}\t{n}\t{m}\n') for x,n,m in rows]
