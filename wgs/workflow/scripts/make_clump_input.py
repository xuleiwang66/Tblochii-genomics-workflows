#!/usr/bin/env python3
import argparse,csv
p=argparse.ArgumentParser(); p.add_argument('--candidates',required=True); p.add_argument('--out',required=True); a=p.parse_args()
with open(a.candidates,encoding='utf-8-sig') as h: rows=list(csv.DictReader(h,delimiter='\t'))
with open(a.out,'w',encoding='utf-8') as h:
    h.write('SNP\tP\n'); [h.write(f"{r['snp_id']}\t{r['p_value']}\n") for r in rows]
