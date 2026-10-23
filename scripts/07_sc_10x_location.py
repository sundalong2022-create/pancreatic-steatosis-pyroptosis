#!/usr/bin/env python3
"""Run the same marker/module localization on a standard 10x matrix."""
import argparse, gzip
from pathlib import Path
import numpy as np, pandas as pd
from scipy.io import mmread

CORE=["TXNIP","NLRP3","CASP1","GSDMD","IL1B","IL18"]
EXT=CORE+["CASP4","CASP5","PYCARD","AIM2","NLRC4","NLRP1","GSDME","GBP1","GBP5"]
MARKERS={"macrophage":["LST1","TYROBP","FCER1G","CTSS","LILRB1","AIF1","C1QC","CD68","MS4A7"],"acinar":["PRSS1","PRSS2","REG1A","REG1B","AMY2A","AMY2B","CPA1","CTRB1","CTRB2"],"ductal":["KRT19","KRT8","KRT18","KRT7","MUC1","KRT17","KRT23","CFTR","SPP1","KRT16"],"stellate":["COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","SPARC","DCN","LUM","DES"],"beta":["INS","IAPP","PCSK1","ISL1","MAFA","DLK1","ADCYAP1"],"alpha":["GCG","TTR","LOXL4"],"endothelial":["KDR","ESAM","EMCN","PECAM1","VWF","RAMP2"],"t_cell":["CD3D","CD3E","TRBC1","TRBC2","IL7R","LTB"],"epithelial_tumor":["EPCAM","KRT8","KRT18","KRT19","MSLN","CEACAM6","MUC1","KRT17"]}
def op(p): return gzip.open(p,"rt",errors="replace") if str(p).endswith(".gz") else open(p)
def canon(x): return str(x).strip().split(".")[0].upper()
def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--matrix",required=True); ap.add_argument("--features",required=True); ap.add_argument("--barcodes",required=True); ap.add_argument("--outdir",default="results/sc_location"); a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
 genes=[]
 with op(a.features) as f:
  for line in f: genes.append(canon(line.rstrip("\n").split("\t")[1] if "\t" in line else line.rstrip("\n")))
 cells=[x.rstrip("\n") for x in op(a.barcodes)]
 wanted=set(CORE+EXT+[g for x in MARKERS.values() for g in x]); idx={g:i for i,g in enumerate(genes) if g in wanted}; mat=mmread(a.matrix).tocsr(); rows=sorted(idx.items(),key=lambda z:z[1]); raw=pd.DataFrame(mat[[i for _,i in rows],:].toarray(),index=[g for g,_ in rows],columns=cells)
 lib=raw.sum(axis=0).replace(0,np.nan); logcpm=np.log1p(raw.div(lib,axis=1)*1e4).fillna(0); s=pd.DataFrame(index=cells)
 for n,gs in MARKERS.items():
  p=[g for g in gs if g in logcpm.index]; s[n+"_score"]=logcpm.loc[p].mean(axis=0) if p else 0.0
 for n,gs in [("PSPS_core",CORE),("pyroptosis_extended",EXT)]:
  p=[g for g in gs if g in logcpm.index]; s[n]=logcpm.loc[p].mean(axis=0) if p else 0.0
 mc=[n+"_score" for n in MARKERS]; s["cell_type"]=s[mc].idxmax(axis=1).str.replace("_score","",regex=False); s["top_score"]=s[mc].max(axis=1); s["second_score"]=s[mc].apply(lambda x:np.sort(x.to_numpy())[-2],axis=1); s["assignment_margin"]=s.top_score-s.second_score; s["low_information"]=s.top_score<=0
 s.to_csv(out/"per_cell_scores.csv"); raw.to_csv(out/"selected_gene_counts.csv"); g=s.groupby("cell_type",observed=True); g[["PSPS_core","pyroptosis_extended","top_score","assignment_margin"]].agg(["count","mean","median","std"]).to_csv(out/"cell_type_score_summary.csv"); (raw>0).T.assign(cell_type=s.cell_type.to_numpy()).groupby("cell_type",observed=True).mean().to_csv(out/"cell_type_gene_detection.csv"); g.agg(n_cells=("PSPS_core","size"),mean_PSPS=("PSPS_core","mean"),median_PSPS=("PSPS_core","median"),mean_extended=("pyroptosis_extended","mean"),frac_low_information=("low_information","mean")).sort_values("mean_PSPS",ascending=False).to_csv(out/"cell_type_pyroptosis_ranking.csv"); pd.DataFrame({"gene":EXT,"present_in_matrix":[x in raw.index for x in EXT]}).to_csv(out/"gene_presence.csv",index=False); print("cells",len(cells),"selected_genes",len(raw))
if __name__=="__main__": main()
