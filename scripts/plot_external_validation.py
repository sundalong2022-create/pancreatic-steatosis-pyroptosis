import pandas as pd, matplotlib.pyplot as plt, seaborn as sns, os
os.makedirs('results/external_validation/figures',exist_ok=True)
sns.set_theme(style='whitegrid',context='paper')
# Human scatter panel
h=pd.read_csv('results/external_validation/GSE79668_PSPS_immune_scores.csv')
mods=['leukocyte','myeloid','macrophage_TAM','inflammation_noncore']
fig,axs=plt.subplots(2,2,figsize=(8,7))
for ax,m in zip(axs.flat,mods):
 sns.regplot(data=h,x='PSPS_core',y=m,ax=ax,scatter_kws={'s':20,'alpha':.75},line_kws={'color':'#D95F0E'})
 ax.set_title(m.replace('_',' ')); ax.set_xlabel('PSPS_core'); ax.set_ylabel('module score')
fig.tight_layout(); fig.savefig('results/external_validation/figures/GSE79668_PSPS_immune_scatter.pdf'); fig.savefig('results/external_validation/figures/GSE79668_PSPS_immune_scatter.png',dpi=300); plt.close(fig)
# Mouse boxplot
m=pd.read_csv('results/external_validation/GSE293790_PSPS_sample_scores.csv')
fig,ax=plt.subplots(figsize=(5,4)); sns.boxplot(data=m,x='group',y='PSPS_core',order=['CD','DIO','CDDIO'],ax=ax,color='#9ecae1',fliersize=0); sns.stripplot(data=m,x='group',y='PSPS_core',order=['CD','DIO','CDDIO'],ax=ax,color='black',size=5,jitter=.08); ax.set_xlabel('diet group'); ax.set_ylabel('PSPS_core'); ax.set_title('Independent mouse pancreatic cohort'); fig.tight_layout(); fig.savefig('results/external_validation/figures/GSE293790_PSPS_boxplot.pdf'); fig.savefig('results/external_validation/figures/GSE293790_PSPS_boxplot.png',dpi=300); plt.close(fig)
