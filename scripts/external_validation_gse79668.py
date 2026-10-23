import pandas as pd, numpy as np, os
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
root='.'; outdir='results/external_validation'; os.makedirs(outdir,exist_ok=True)
raw=pd.read_csv('data/external/GSE79668_counts.txt.gz',sep='\t')
raw['SYMBOL']=raw['Gene_EnsembleID'].astype(str).str.split('_ENSG',n=1).str[0].str.upper(); raw=raw[raw.SYMBOL.notna()]
sample_cols=[c for c in raw.columns if c not in ['Gene_EnsembleID','SYMBOL']]; counts=raw.groupby('SYMBOL')[sample_cols].sum(); lib=counts.sum(axis=0); logcpm=np.log2(counts.div(lib,axis=1)*1e6+1)
core=['TXNIP','NLRP3','CASP1','GSDMD','IL1B','IL18']
sets={'leukocyte':['PTPRC','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','LST1'],'t_cell':['CD3D','CD3E','CD3G','CD4','CD8A','IL7R','LTB'],'cytotoxic':['NKG7','GZMB','GZMK','PRF1','CTSW','GZMA'],'b_cell':['MS4A1','CD79A','CD79B','CD37','CD74','HLA-DRA','HLA-DRB1','CD19','CD83'],'myeloid':['LST1','TYROBP','AIF1','CTSS','CTSD','LGALS3','S100A10'],'macrophage_TAM':['ADGRE1','CD68','CSF1R','C1QA','C1QB','C1QC','MS4A7','MERTK','MAFB','CD163','MRC1','MARCO','LGALS3'],'neutrophil':['S100A8','S100A9','LCN2','CXCR2','CSF3R','NGP','CAMP'],'inflammation_noncore':['TNF','IL6','CXCL1','CXCL2','CXCL5','CCL2','CCL3','CCL4','CCL5','NFKB1','RELA','TLR4','TGFB1','PTGS2','S100A8','S100A9','NFKBIA','TNFAIP3']}
def score(gs):
 p=[g for g in gs if g in logcpm.index and g not in core]; m=logcpm.loc[p]; z=m.sub(m.mean(axis=1),axis=0).div(m.std(axis=1,ddof=1).replace(0,np.nan),axis=0).fillna(0); return z.mean(axis=0),p
m=logcpm.loc[core]; z=m.sub(m.mean(axis=1),axis=0).div(m.std(axis=1,ddof=1).replace(0,np.nan),axis=0).fillna(0); scores=pd.DataFrame({'sample':logcpm.columns,'PSPS_core':z.mean(axis=0).values})
rows=[]; cov=[]
for nm,gs in sets.items():
 s,p=score(gs); scores[nm]=s.to_numpy(); cov.append({'module':nm,'n_requested':len(gs),'n_present':len(p),'coverage':len(p)/len(gs),'genes_present':';'.join(p),'genes_missing':';'.join([g for g in gs if g not in p])}); ok=np.isfinite(scores['PSPS_core'].to_numpy())&np.isfinite(s.to_numpy()); rho,pv=spearmanr(scores.loc[ok,'PSPS_core'].to_numpy(),s.to_numpy()[ok]); rows.append({'module':nm,'n':int(ok.sum()),'rho':rho,'p_value':pv})
res=pd.DataFrame(rows); res['FDR_BH']=multipletests(res.p_value,method='fdr_bh')[1]
scores.to_csv(outdir+'/GSE79668_PSPS_immune_scores.csv',index=False); res.to_csv(outdir+'/GSE79668_PSPS_immune_correlations.csv',index=False); pd.DataFrame(cov).to_csv(outdir+'/GSE79668_immune_module_coverage.csv',index=False); logcpm.loc[core].T.reset_index(names='sample').to_csv(outdir+'/GSE79668_core_log2cpm.csv',index=False)
print(res.to_string(index=False)); print(pd.DataFrame(cov).to_string(index=False))
