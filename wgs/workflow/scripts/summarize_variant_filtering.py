#!/usr/bin/env python3
import argparse,subprocess,csv
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--chromosomes',required=True); p.add_argument('--raw','--raw-pattern',dest='raw',required=True); p.add_argument('--biallelic',required=True); p.add_argument('--pass-pattern',required=True); p.add_argument('--final',required=True); p.add_argument('--out-by-chr',required=True); p.add_argument('--out-summary',required=True); a=p.parse_args()
chroms=[x.strip() for x in open(a.chromosomes) if x.strip() and not x.startswith('#')]
def count(path):
    cmd=f"bcftools view -H {path!s} | wc -l"
    return int(subprocess.check_output(cmd,shell=True,text=True).strip())
rows=[]
for c in chroms:
    vals={
      'chromosome':c,
      'raw_variants':count(a.raw.format(chrom=c)),
      'biallelic_snps':count(a.biallelic.format(chrom=c)),
      'pass_snps':count(a.pass_pattern.format(chrom=c)),
      'callrate_filtered_snps':count(a.final.format(chrom=c)),
    }
    vals['hard_filter_removed']=vals['biallelic_snps']-vals['pass_snps']
    vals['callrate_removed']=vals['pass_snps']-vals['callrate_filtered_snps']
    rows.append(vals)
Path(a.out_by_chr).parent.mkdir(parents=True,exist_ok=True)
with open(a.out_by_chr,'w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
tot={k:sum(r[k] for r in rows) for k in rows[0] if k!='chromosome'}
with open(a.out_summary,'w',encoding='utf-8') as h:
    h.write('metric\tvalue\n'); [h.write(f'{k}\t{v}\n') for k,v in tot.items()]
