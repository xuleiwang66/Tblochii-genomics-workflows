#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--bim',required=True); p.add_argument('--candidate-files',nargs='+',required=True); p.add_argument('--effect-files',nargs='*',default=[]); p.add_argument('--models',nargs='+',required=True); p.add_argument('--out',required=True); a=p.parse_args()
ref={}
with open(a.bim,encoding='utf-8',errors='replace') as h:
    for line in h:
        if not line.strip(): continue
        f=line.split()
        if len(f)>=6: ref[f[1]]={'chromosome':f[0],'position':f[3],'A1':f[4],'A2':f[5]}
effects={}
for model,path in zip(a.models,a.effect_files):
    if not Path(path).exists(): continue
    with open(path,encoding='utf-8-sig') as h:
        for r in csv.DictReader(h,delimiter='\t'): effects[(model,r['snp_id'])]=r
rows=[]
for model,path in zip(a.models,a.candidate_files):
    with open(path,encoding='utf-8-sig') as h:
        for r in csv.DictReader(h,delimiter='\t'):
            snp=r['snp_id']; b=ref.get(snp); fx=effects.get((model,snp),{})
            if model=='continuous': native_formal=r.get('formal_effect','NA'); formal_se=r.get('formal_effect_se','NA'); formal_p=r.get('p_value','NA')
            elif model=='binary': native_formal=fx.get('beta_log_or','NA'); formal_se=fx.get('se','NA'); formal_p=fx.get('p_wald','NA')
            else: native_formal=fx.get('beta_log_hr','NA'); formal_se=fx.get('se','NA'); formal_p=fx.get('p_wald','NA')
            if not b:
                sign=None; status='UNMAPPED'; common_ea=common_oa='NA'
            else:
                common_ea,common_oa=b['A1'],b['A2']; ea=r['effect_allele']; oa=r['other_allele']
                if ea==common_ea and oa==common_oa: sign=1; status='DIRECT'
                elif ea==common_oa and oa==common_ea: sign=-1; status='REVERSED'
                else: sign=None; status='ALLELE_MISMATCH'
            def align(x, cox_flip=False):
                try:
                    v=float(x)
                    if cox_flip: v=-v
                    return f'{sign*v:.12g}' if sign is not None else 'NA'
                except:return 'NA'
            aligned_formal=align(native_formal,cox_flip=(model=='cox'))
            rows.append({**r,'model':model,'common_effect_allele':common_ea,'common_other_allele':common_oa,'allele_alignment':status,'aligned_tolerance_effect':align(r.get('tolerance_effect')),'aligned_tolerance_test_z':align(r.get('tolerance_test_z')),'native_formal_effect':native_formal,'aligned_tolerance_formal_effect':aligned_formal,'formal_effect_se':formal_se,'formal_effect_p':formal_p})
fields=[]
for r in rows:
    for k in r:
        if k not in fields: fields.append(k)
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
with open(a.out,'w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
