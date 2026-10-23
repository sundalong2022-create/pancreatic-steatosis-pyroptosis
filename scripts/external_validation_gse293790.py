import pandas as pd, numpy as np, os
from scipy.stats import mannwhitneyu
os.makedirs('results/external_validation',exist_ok=True)
x=pd.read_csv('data/external/GSE293790_counts.txt.gz',sep='\t').drop_duplicates('SYMBOL').set_index('SYMBOL')
core=['Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18']; cols=[c for c in x.columns if c!='ENTREZID']; m=np.log2(x.loc[core,cols].astype(float)+1)
z=m.sub(m.mean(axis=1),axis=0).div(m.std(axis=1,ddof=1).replace(0,np.nan),axis=0).fillna(0)
out=pd.DataFrame({'sample':z.columns,'PSPS_core':z.mean(axis=0).values}); out['group']=out['sample'].str.extract(r'^(CDDIO|DIO|CD)')[0]
for g in core: out[g]=m.loc[g].values
rows=[]
for a,b in [('DIO','CD'),('CDDIO','CD'),('CDDIO','DIO')]:
 xa=out.loc[out.group==a,'PSPS_core']; xb=out.loc[out.group==b,'PSPS_core']; _,p=mannwhitneyu(xa,xb,alternative='two-sided'); rows.append({'comparison':a+'_vs_'+b,'n_a':len(xa),'n_b':len(xb),'mean_a':xa.mean(),'mean_b':xb.mean(),'delta':xa.mean()-xb.mean(),'mw_p':p})
out.to_csv('results/external_validation/GSE293790_PSPS_sample_scores.csv',index=False); pd.DataFrame(rows).to_csv('results/external_validation/GSE293790_PSPS_group_tests.csv',index=False)
long=[]
for g in core:
 for grp,sub in out.groupby('group'):
  long.append({'gene':g,'group':grp,'n':len(sub),'mean_log2_norm':sub[g].mean()})
pd.DataFrame(long).to_csv('results/external_validation/GSE293790_core_gene_group_means.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
