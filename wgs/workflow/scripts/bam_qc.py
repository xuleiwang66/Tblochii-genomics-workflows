#!/usr/bin/env python3
import argparse,csv,re
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--sample',required=True); p.add_argument('--flagstat',required=True); p.add_argument('--depth',required=True); p.add_argument('--dup-metrics',required=True); p.add_argument('--out',required=True); a=p.parse_args()
text=Path(a.flagstat).read_text(errors='replace')
def pct(pattern):
    m=re.search(pattern,text,re.M); return float(m.group(1)) if m else float('nan')
total_m=re.search(r'^(\d+) \+ \d+ in total',text,re.M)
metrics={'sample_id':a.sample,'total_reads':int(total_m.group(1)) if total_m else 'NA','mapped_rate_pct':pct(r'mapped \(([0-9.]+)%'),'properly_paired_rate_pct':pct(r'properly paired \(([0-9.]+)%')}
with open(a.depth) as h:
    d=dict(csv.reader(h,delimiter='\t'))
for k in ['mean_depth_chr01_to_chr24','breadth_1x_chr01_to_chr24_pct','breadth_3x_chr01_to_chr24_pct','breadth_5x_chr01_to_chr24_pct','breadth_10x_chr01_to_chr24_pct']:
    metrics[k]=d.get(k,'NA')
dup='NA'
lines=Path(a.dup_metrics).read_text(errors='replace').splitlines()
for i,line in enumerate(lines):
    if line.startswith('LIBRARY\t') and 'PERCENT_DUPLICATION' in line and i+1<len(lines):
        hdr=line.split('\t'); vals=lines[i+1].split('\t'); dup=float(vals[hdr.index('PERCENT_DUPLICATION')])*100; break
metrics['duplicate_rate_pct']=dup
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
with open(a.out,'w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=list(metrics),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerow(metrics)
