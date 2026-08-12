#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--samples',required=True); p.add_argument('--gvcf-dir',required=True); p.add_argument('--out',required=True); a=p.parse_args()
with open(a.samples,encoding='utf-8-sig') as h:
    rows=list(csv.DictReader((x for x in h if not x.startswith('#')),delimiter='\t'))
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
with open(a.out,'w',encoding='utf-8',newline='') as h:
    for r in rows: h.write(f"{r['sample_id']}\t{Path(a.gvcf_dir).resolve()/(r['sample_id']+'.g.vcf.gz')}\n")
