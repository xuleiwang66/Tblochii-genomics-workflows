#!/usr/bin/env python3
import argparse
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--chromosomes',required=True); p.add_argument('--prefix-pattern',required=True); p.add_argument('--out',required=True); a=p.parse_args()
chroms=[x.strip() for x in open(a.chromosomes) if x.strip() and not x.startswith('#')]
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
Path(a.out).write_text(''.join(a.prefix_pattern.format(chrom=c)+'\n' for c in chroms),encoding='utf-8')
