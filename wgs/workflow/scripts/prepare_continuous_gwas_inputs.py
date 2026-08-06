#!/usr/bin/env python3
import argparse,csv,numpy as np
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--psam',required=True); p.add_argument('--metadata',required=True); p.add_argument('--pca',required=True); p.add_argument('--kinship-order',required=True); p.add_argument('--phenotype',required=True); p.add_argument('--covariates',required=True); p.add_argument('--annotated',required=True); p.add_argument('--sample-order',required=True); a=p.parse_args()
def rows(path):
    with open(path,encoding='utf-8-sig') as h: return list(csv.DictReader(h,delimiter='\t'))
ps=rows(a.psam); ids=[r.get('IID') or r.get('#IID') or r.get('sample_id') for r in ps]
meta={r['sample_id']:r for r in rows(a.metadata)}; pca={r['sample_id']:r for r in rows(a.pca)}
ko=rows(a.kinship_order); kids=[r.get('IID') or r.get('sample_id') for r in ko]
if ids!=kids: raise SystemExit('Kinship order does not match PSAM order')
missing=[x for x in ids if x not in meta or x not in pca]
if missing: raise SystemExit(f'Missing metadata/PCA rows: {missing[:10]}')
X=[]; y=[]; ann=[]
for i,sid in enumerate(ids,1):
    m=meta[sid]; q=pca[sid]
    yy=float(m['LOE_time_h']); row=[1.0,float(m['body_weight_g']),float(m['source_L']),float(m['source_Q']),float(q['PC1']),float(q['PC2']),float(q['PC3'])]
    y.append(yy); X.append(row); ann.append([i,sid,yy]+row[1:])
if np.linalg.matrix_rank(np.asarray(X,dtype=float)) != 7: raise SystemExit('Continuous-GWAS covariate matrix is not full rank')
for path in [a.phenotype,a.covariates,a.annotated,a.sample_order]: Path(path).parent.mkdir(parents=True,exist_ok=True)
Path(a.phenotype).write_text(''.join(f'{x:.12g}\n' for x in y),encoding='utf-8')
Path(a.covariates).write_text(''.join('\t'.join(f'{v:.12g}' for v in r)+'\n' for r in X),encoding='utf-8')
with open(a.sample_order,'w',newline='',encoding='utf-8') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['order_index','IID']); [w.writerow([i,s]) for i,s in enumerate(ids,1)]
with open(a.annotated,'w',newline='',encoding='utf-8') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['order_index','sample_id','LOE_time_h','body_weight_g','source_L','source_Q','PC1','PC2','PC3']); w.writerows(ann)
