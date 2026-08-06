#!/usr/bin/env python3
import argparse, csv, math
from pathlib import Path

MISSING = {"", ".", "NA", "NaN", "nan", "NULL"}

def read_table(path):
    with open(path, encoding="utf-8-sig", newline="") as h:
        rows=[r for r in csv.DictReader((x for x in h if not x.startswith('#')), delimiter='\t')]
    return rows

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--phenotypes', required=True)
    p.add_argument('--samples', required=True)
    p.add_argument('--clean', required=True)
    p.add_argument('--analysis-samples', required=True)
    p.add_argument('--pre-target-deaths', required=True)
    p.add_argument('--censor-time', type=float, default=24.0)
    p.add_argument('--exclude-status', nargs='*', type=int, default=[2])
    a=p.parse_args()

    sample_rows=read_table(a.samples)
    fastq_ids=[r['sample_id'].strip() for r in sample_rows]
    if len(fastq_ids) != len(set(fastq_ids)):
        raise SystemExit('Duplicate sample_id values in samples.tsv')
    pheno=read_table(a.phenotypes)
    required=['sample_id','fish_id','source_cohort','LOE_time_h','LOE_status','body_weight_g']
    missing=[x for x in required if not pheno or x not in pheno[0]]
    if missing: raise SystemExit(f'Missing phenotype columns: {missing}')
    if len(pheno) != len({r['sample_id'].strip() for r in pheno}):
        raise SystemExit('Duplicate sample_id values in phenotypes.tsv')
    pmap={r['sample_id'].strip(): r for r in pheno}
    if set(fastq_ids) != set(pmap):
        raise SystemExit(f'Sample-set mismatch: FASTQ-only={sorted(set(fastq_ids)-set(pmap))[:10]}, phenotype-only={sorted(set(pmap)-set(fastq_ids))[:10]}')

    out=[]; analysis=[]; pre=[]
    for sid in fastq_ids:
        r=pmap[sid]
        status=int(float(r['LOE_status']))
        if status not in {0,1,2}: raise SystemExit(f'Invalid LOE_status for {sid}: {status}')
        time=float(r['LOE_time_h']) if r['LOE_time_h'] not in MISSING else math.nan
        weight=float(r['body_weight_g'])
        if weight <= 0: raise SystemExit(f'body_weight_g must be positive for {sid}')
        if status == 0: time=a.censor_time
        if status == 2: time=math.nan
        event=1 if status==1 else 0
        valid=status not in set(a.exclude_status)
        cohort=r['source_cohort'].strip()
        if cohort not in {'C','L','Q'}: raise SystemExit(f'Unexpected source_cohort for {sid}: {cohort}')
        if status == 1 and (math.isnan(time) or time <= 0 or time > a.censor_time): raise SystemExit(f'Invalid event LOE_time_h for {sid}: {time}')
        out.append({
            'sample_id':sid,'fish_id':r['fish_id'].strip(),'source_cohort':cohort,
            'LOE_time_h':'NA' if math.isnan(time) else f'{time:.12g}',
            'LOE_status':status,'time_to_event_h':'NA' if math.isnan(time) else f'{time:.12g}',
            'event':event,'body_weight_g':f'{weight:.12g}',
            'source_L':int(cohort=='L'),'source_Q':int(cohort=='Q'),
            'valid_downstream':int(valid)
        })
        if valid: analysis.append(sid)
        else: pre.append(sid)

    Path(a.clean).parent.mkdir(parents=True, exist_ok=True)
    fields=list(out[0])
    with open(a.clean,'w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(out)
    Path(a.analysis_samples).write_text('\n'.join(analysis)+'\n',encoding='utf-8')
    Path(a.pre_target_deaths).write_text('\n'.join(pre)+('\n' if pre else ''),encoding='utf-8')
if __name__=='__main__': main()
