#!/usr/bin/env python3
import argparse,csv,gzip,math
from pathlib import Path
from statistics import median
from scipy.stats import chi2
p=argparse.ArgumentParser(); p.add_argument('--model',choices=['continuous','binary','cox'],required=True); p.add_argument('--input',required=True); p.add_argument('--out',required=True); p.add_argument('--candidates',required=True); p.add_argument('--summary',required=True); p.add_argument('--threshold',type=float,default=1e-5); a=p.parse_args()

def first(d,names):
    for n in names:
        if n in d and d[n] not in {'','NA','.','nan','NaN'}: return d[n]
    return None

def opener(path,mode='rt'):
    return gzip.open(path,mode) if str(path).endswith('.gz') else open(path,mode,encoding='utf-8')

def iter_dict_rows(path):
    with opener(path) as h:
        header_line = h.readline().rstrip("\n\r")
        if not header_line:
            return
        if "\t" in header_line:
            header = header_line.split("\t")
            for line in h:
                if line.strip():
                    values = line.rstrip("\n\r").split("\t")
                    yield dict(zip(header, values))
        else:
            header = header_line.split()
            for line in h:
                if line.strip():
                    values = line.split()
                    yield dict(zip(header, values))

rows=[]; stats=[]
for r in iter_dict_rows(a.input):
    if a.model=='continuous':
        chrom=first(r,['chr','CHR','#CHROM']); pos=first(r,['ps','POS','position']); snp=first(r,['rs','SNP','ID']); pval=first(r,['p_wald','P','PVAL']); eff=first(r,['beta']); se=first(r,['se']); ea=first(r,['allele1']); oa=first(r,['allele0']); af=first(r,['af'])
        if pval is None: continue
        pv=float(pval); stat=chi2.isf(max(pv,1e-300),1); tol_eff=float(eff) if eff is not None else math.nan; score='NA'; test_z=(float(eff)/float(se)) if eff is not None and se is not None and float(se)>0 else math.nan
    elif a.model=='binary':
        chrom=first(r,['CHR','chr']); pos=first(r,['POS','position']); snp=first(r,['SNP','ID']); pval=first(r,['PVAL','P','p']); score=first(r,['SCORE']); var=first(r,['VAR']); ea=first(r,['A2','effect_allele']); oa=first(r,['A1','other_allele']); af=first(r,['AF']); se='NA'; eff='NA'
        if pval is None or score is None or var is None: continue
        pv=float(pval); stat=float(score)**2/float(var); tol_eff=float(score); test_z=float(score)/math.sqrt(float(var))
    else:
        chrom=first(r,['original_CHR','CHR','chromosome']); pos=first(r,['original_POS','POS','position']); snp=first(r,['mapped_SNP','SNP','snp.id','index']); pval=first(r,['p','P','PVAL']); score=first(r,['score','cox_sensitivity_score']); score_test=first(r,['score_test']); ea=first(r,['counted_allele','coxmeg_allele1']); oa=first(r,['coxmeg_allele2','other_allele']); af=first(r,['counted_allele_frequency','afreq_inc']); se='NA'; eff='NA'
        if pval is None or score is None: continue
        pv=float(pval); stat=float(score_test) if score_test is not None else chi2.isf(max(pv,1e-300),1); tol_eff=-float(score); test_z=-math.copysign(math.sqrt(max(stat,0.0)),float(score))
    if not (0 <= pv <= 1): raise SystemExit(f'Invalid P value for {snp}: {pv}')
    stats.append(stat)
    rows.append({'model':a.model,'chromosome':chrom,'position':pos,'snp_id':snp,'effect_allele':ea or 'NA','other_allele':oa or 'NA','effect_allele_frequency':af or 'NA','p_value':f'{pv:.12g}','score_statistic':score if a.model!='continuous' else 'NA','formal_effect':eff if a.model=='continuous' else 'NA','formal_effect_se':se if a.model=='continuous' else 'NA','tolerance_effect':f'{tol_eff:.12g}' if math.isfinite(tol_eff) else 'NA','tolerance_test_z':f'{test_z:.12g}' if math.isfinite(test_z) else 'NA','formal_effect_available':'YES' if a.model=='continuous' else 'NO'})
if not rows: raise SystemExit('No association rows parsed')
lam=median(stats)/chi2.ppf(0.5,1)
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
with gzip.open(a.out,'wt',encoding='utf-8',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
cand=[r for r in rows if float(r['p_value']) <= a.threshold]
with open(a.candidates,'w',encoding='utf-8',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(cand)
with open(a.summary,'w',encoding='utf-8') as h:
    h.write('metric\tvalue\n'); h.write(f'model\t{a.model}\n'); h.write(f'tested_variants\t{len(rows)}\n'); h.write(f'candidate_variants\t{len(cand)}\n'); h.write(f'lambda_GC\t{lam:.12g}\n'); h.write(f'minimum_P\t{min(float(r["p_value"]) for r in rows):.12g}\n')
