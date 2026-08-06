#!/usr/bin/env python3
import argparse
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--fai',required=True); p.add_argument('--chromosomes',required=True); p.add_argument('--bed',required=True); p.add_argument('--validated',required=True); a=p.parse_args()
chroms=[x.strip() for x in open(a.chromosomes) if x.strip() and not x.startswith('#')]
lengths={}
for line in open(a.fai):
    f=line.rstrip().split('\t'); lengths[f[0]]=int(f[1])
missing=[c for c in chroms if c not in lengths]
if missing: raise SystemExit(f'Configured chromosomes absent from FASTA index: {missing}')
Path(a.bed).parent.mkdir(parents=True,exist_ok=True)
Path(a.bed).write_text(''.join(f'{c}\t0\t{lengths[c]}\n' for c in chroms),encoding='utf-8')
Path(a.validated).write_text(f'validated_chromosomes\t{len(chroms)}\n',encoding='utf-8')
