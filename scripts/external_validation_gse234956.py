import pandas as pd,numpy as np,os
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
os.makedirs('results/external_validation',exist_ok=True)
x=pd.read_csv('data/external/GSE234956/counts_matrix.csv.gz')
meta=pd.read_csv('data/external/GSE234956/sample_metadata.csv')
# Matrix columns are ORSAM IDs; map to GEO sample metadata through title token.
meta['matrix_id']=meta['title'].str.extract(r'\((ORSAM\d+)\)')[0]
cols=[c for c in x.columns if c.startswith('ORSAM')]
meta=meta[meta.matrix_id.isin(cols)].copy(); meta=meta.set_index('matrix_id').loc[cols].reset_index()
mat=x.drop_duplicates('Gene_symbol').set_index('Gene_symbol')[cols].astype(float)
core=['Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18']
print('matrix',mat.shape,'meta',meta.shape,'missing', [g for g in core if g not in mat.index])
lib=mat.sum(axis=0); logcpm=np.log2(mat.div(lib,axis=1)*1e6+1)
z=logcpm.loc[core].sub(logcpm.loc[core].mean(axis=1),axis=0).div(logcpm.loc[core].std(axis=1,ddof=1).replace(0,np.nan),axis=0).fillna(0)
out=pd.DataFrame({'matrix_id':cols,'PSPS_core':z.mean(axis=0).values}).merge(meta[['matrix_id','gsm','diet','strain','sex']],on='matrix_id')
out['HFHS']=(out.diet.str.startswith('44')).astype(int)
for g in core: out[g]=logcpm.loc[g,cols].values
# project immune markers, human? mouse islet set
sets={'leukocyte':['Ptprc','H2-Aa','H2-Ab1','H2-Eb1','Cd74','Lsp1'],'t_cell':['Cd3d','Cd3e','Cd3g','Cd4','Cd8a','Il7r','Ltb'],'cytotoxic':['Nkg7','Gzmb','Gzmk','Prf1','Ctsw','Gzma'],'b_cell':['Ms4a1','Cd79a','Cd79b','Cd37','Cd74','H2-Aa','H2-Ab1','Cd19','Cd83'],'myeloid':['Lst1','Tyrobp','Aif1','Ctss','Ctsd','Lgals3','S100a10'],'macrophage_TAM':['Adgre1','Cd68','Csf1r','C1qa','C1qb','C1qc','Ms4a7','Mertk','Mafb','Cd163','Mrc1','Marco','Lgals3'],'inflammation_noncore':['Tnf','Il6','Cxcl1','Cxcl2','Cxcl5','Ccl2','Ccl3','Ccl4','Ccl5','Nfkb1','Rela','Tlr4','Tgfb1','Ptgs2','S100a8','S100a9','Nfkbia','Tnfaip3']}
for nm,gs in sets.items():
 p=[g for g in gs if g in logcpm.index and g not in core]
 if not p: out[nm]=np.nan; continue
 m=logcpm.loc[p]; zz=m.sub(m.mean(axis=1),axis=0).div(m.std(axis=1,ddof=1).replace(0,np.nan),axis=0).fillna(0); out[nm]=zz.mean(axis=0).reindex(cols).values
# overall model and per-strain contrasts
rows=[]
for feat in ['PSPS_core']+core+list(sets):
 try:
  fit=smf.ols(f'{feat} ~ HFHS + C(strain) + C(sex)',data=out).fit()
  ci=fit.conf_int().loc['HFHS']; rows.append({'feature':feat,'n':len(out),'effect_HFHS':fit.params['HFHS'],'ci_low':ci[0],'ci_high':ci[1],'p_value':fit.pvalues['HFHS']})
 except Exception as e: rows.append({'feature':feat,'n':len(out),'effect_HFHS':np.nan,'ci_low':np.nan,'ci_high':np.nan,'p_value':np.nan})
res=pd.DataFrame(rows); res['FDR_BH']=multipletests(res.p_value.fillna(1),method='fdr_bh')[1]
# per strain PSPS
ps=[]
for s,g in out.groupby('strain'):
 a=g[g.HFHS==1].PSPS_core; b=g[g.HFHS==0].PSPS_core; _,p=mannwhitneyu(a,b,alternative='two-sided'); ps.append({'strain':s,'n_HFHS':len(a),'n_control':len(b),'mean_HFHS':a.mean(),'mean_control':b.mean(),'delta':a.mean()-b.mean(),'p_value':p})
pd.DataFrame(ps).to_csv('results/external_validation/GSE234956_PSPS_by_strain.csv',index=False)
out.to_csv('results/external_validation/GSE234956_PSPS_sample_scores.csv',index=False); res.to_csv('results/external_validation/GSE234956_HFHS_effects_adjusted_strain_sex.csv',index=False)
print(res[res.feature.isin(['PSPS_core']+core)].to_string(index=False)); print(pd.DataFrame(ps).to_string(index=False)); print('coverage', {nm:sum(g in logcpm.index for g in gs) for nm,gs in sets.items()})
