import sys
print(sys.path)
print(sys.executable)
import scanpy as sc
import anndata as ad
import os
import seaborn as sns
import matplotlib as plt
import pyreadr
import leidenalg
import rpy2
import rpy2.robjects as ro
import scipy
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob
import warnings
warnings.filterwarnings('ignore')
import logging
from scipy.stats import median_abs_deviation
import rpy2
from rpy2.robjects import pandas2ri
from rpy2.robjects import numpy2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import anndata2ri
anndata2ri.activate()
%load_ext rpy2.ipython

import pathlib
from pathlib import Path

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],      
    "savefig.dpi": 300,               
    "savefig.format": "pdf",           
    "savefig.facecolor": "white",      
    "savefig.edgecolor": "black",      
})

sc.set_figure_params(
    dpi=300,              
    dpi_save=300,         
    format="pdf",         
    frameon=True,         
    transparent=True,     
    facecolor="white",    
)

ADMicroglia = sc.read_h5ad("~/Merged/Finalised/ADMicroglia.h5ad")

# Cluster vs cluster
sc.tl.dendrogram(ADMicroglia, groupby="leiden0_7")
sc.tl.rank_genes_groups(
    ADMicroglia, 
    groupby="leiden0_7", 
    method="wilcoxon", 
    corr_method="benjamini-hochberg",
    key_added="deg_leiden07")
sc.tl.filter_rank_genes_groups(
    ADMicroglia,
    min_in_group_fraction=0.2,
    key="deg_leiden07",
    key_added="deg_f_leiden07"
)
sc.pl.rank_genes_groups_dotplot(
    ADMicroglia,
    groupby="leiden0_7",
    standard_scale="var",
    n_genes=8,
    key="deg_f_leiden07",
    save="1_AnnoLeiden07_Dotplot"
)
sc.pl.dotplot(
    ADMicroglia,
    groupby="leiden0_7",
    var_names = marker_genes_in_data,
    standard_scale = 'var',
    save = '2_AnnoLeiden07_CellIdentities'
)
sc.pl.umap(ADMicroglia,color=["leiden0_7"],size=5,save="3_UMAP_scVI_L07")
deg_res_leiden07 = ADMicroglia.uns['deg_leiden07']
ADMicroglia.obs['leiden0_7'].value_counts()

conditions = [
    ADMicroglia.obs['leiden0_7'] == '0',
    ADMicroglia.obs['leiden0_7'] == '1',
    ADMicroglia.obs['leiden0_7'] == '2',
    ADMicroglia.obs['leiden0_7'] == '3',
    ADMicroglia.obs['leiden0_7'] == '4',
    ADMicroglia.obs['leiden0_7'] == '5',
    ADMicroglia.obs['leiden0_7'] == '6',
    ADMicroglia.obs['leiden0_7'] == '7',
    ADMicroglia.obs['leiden0_7'] == '8',
    ADMicroglia.obs['leiden0_7'] == '9',
    ADMicroglia.obs['leiden0_7'] == '10',
]

choices = [
    'MG0','MG1','MG2','MG3','MG4',
    'MG5','MG6','MG7','MG8','MG9','MG10'
]

ADMicroglia.obs['final_resolution'] = np.select(
    conditions,
    choices,
    default=ADMicroglia.obs['leiden0_7']
)

sc.pl.umap(ADMicroglia, color="final_resolution", legend_fontsize=10, legend_fontoutline=2, frameon=True, size=5, save="13_ADMicroglia_leiden07")
ADMicroglia.obs['final_resolution'].value_counts()

clusters = deg_res_leiden07['names'].dtype.names
deg_tables = {}
for cluster in clusters:
    gene_names = deg_res_leiden07['names'][cluster]
    scores = deg_res_leiden07['scores'][cluster]
    pvals = deg_res_leiden07['pvals'][cluster]
    padj = deg_res_leiden07['pvals_adj'][cluster]
    logFC = deg_res_leiden07['logfoldchanges'][cluster]
    deg_table = pd.DataFrame({
        'gene': gene_names,
        'score': scores,
        'pval': pvals,
        'padj': padj,
        'logFC': logFC,
    })
    deg_table = deg_table.sort_values(by='padj', ascending=True)
    deg_tables[cluster] = deg_table
    deg_table.to_csv(f"/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/FinalisedData/FinalMerge/7_leiden07/4_DEG_ClusvOtherClus/DEG_{cluster}.csv", index=False)

# Iron Enrichment scoring
iron_genes = ["HFE","SLC25A39","B2M","BMP6","TF","TFAP2A","TFR2","CCNB1","ABCB6","GLRX3","SLC46A1","CP","TMEM199","HJV","ERFE","CYB561","TMPRSS6","ALAS2","ABCB7","CYB561A3","ATP13A2","ISCU","SLC39A14","ATP6V0A2","FXN","FTH1","FTL","STEAP2","FBXL5","GDF2","SCARA5","FLVCR1","SLC40A1","HIF1A","HMOX1","HPX","HEPHL1","IFNG","IREB2","BOLA3","FTH1P19","FRRS1","LCN2","MIR210","SMAD1","SMAD4","SMAD5","SMAD9","MYC","NUBP1","NEO1","ACO1","SLC11A2","BOLA1","MECR","GLRX5","SLC22A17","ATP6V1A","ATP6AP1","FTHL17","EGLN1","BOLA2","CISD1","TTC7A","HAMP","SLC39A8","BOLA2B","SLC11A1","SOD1","TFRC","STEAP4","CYBRD1","NCOA4","NDFIP1","PICALM","CCDC115","ATP6V0D1","ACVR2B","FTMT","ATP6V1G1","ABCC5","PGRMC2","FLVCR2","SLC48A1","ASIC3","MCOLN2","SLC25A37","MCOLN1","TTYH1","SLC25A28","MMGT1","HEPH","CLTC","DNM2","LMTK2","REP15","ARHGAP1","LTF","MELTF","STEAP3","NECTIN1","SLC41A2","RAB11B","SFXN1","BTBD9","EPAS1","EPB42","FECH","EIF2AK1","NAGLU","RHAG","SOD2","HYAL2","PDX1","SNCA","MAP1LC3A","BECN1","ATG5","CPOX","CYP1A1","ABAT","DRD2","ALAD","G6PD","MDM2","CCND1","BCL2","SLC6A3","TFF1","ACP5","MIOX","RRM2","CDO1","EGLN2","EGLN3","DNAJC24","ALKBH2","ALKBH3","HAAO","DPH3","PHYH","TET2","HIF1AN","FTO","ALKBH1","P3H3","FAXDC2","CYP46A1","CYP4F8","CYP2U1","CYGB","CYP2R1","TPH2","ISCA2","CYP4F22","TYW5","CYP1A2","CYP1B1","CYP2A6","CYP2A7","CYP3A7","CYP2A13","CYP2B6","CYP2C19","CYP2C8","CYP2C9","CYP2C18","CYP2D6","CYP2E1","CYP2F1","CYP2J2","CYP3A4","CYP3A5","CYP4A11","CYP4B1","CYP7A1","CYP8B1","CYP11A1","CYP11B1","CYP11B2","CYP17A1","CYP19A1","CYP21A2","CYP24A1","CYP26A1","CYP27A1","CYP27B1","CYP51A1","CYP4Z2P","CYP4Z1","FDX1","CYP2G1P","PHF8","JMJD6","ETHE1","ALOX12","ALOX5","ALOX12B","ALOX15","ALOX15B","CYP4X1","NFU1","P4HA3","CYP4A22","CYP4V2","CYP2S1","HBA1","HBA2","HBQ1","HBZ","AOX1","CYP27C1","CYP26C1","AGMO","CYP4F3","ACO2","P4HA1","PAH","CYP39A1","PHF2","PLOD1","PLOD2","P4HTM","PPEF2","PPEF1","CYP2W1","PPP1CA","P3H2","TMLHE","OGFOD1","ADI1","KDM3A","CYP26B1","CIAPIN1","PTGIS","CYP20A1","CYP4F11","ALOXE3","ABCE1","MSMO1","SC5D","SCD","P3H1","CYP3A43","CYP4F12","TBXAS1","TH","TPH1","XDH","FA2H","OGFOD2","RIOX1","OGFOD3","SCD5","TET1","KDM7A","CALR","DOHH","BBOX1","ADO","CYP4F2","P4HA2","PLOD3","CH25H","ALKBH8","CYP7B1","COL7A1","STAB1","FOXP1","RACGAP1","TRNT1","PIGA","PKLR","BCS1L","CARD9","SKIC2","BMP2","KIF23","SKIC3","HBB","SLC25A38","SLC30A10","ABCD3","WDR45","CRAT","DLAT","VPS13A","VPS41","NUDT2","TBCE","PANK2","COASY","C19orf12","PLA2G6","REPS1","ATP6V0A1","ATP6V1H","CUL1","ATP6V1D","HMOX2","ATP6V0A4","TCIRG1","CAND1","SKP1","ATP6V0E1","ATP6V1B1","ATP6V0B","ABCG2","ATP6V1F","NEDD8","ATP6V1E1","ATP6V1C2","RPS27A","ATP6V1B2","ATP6V0D2","UBC","ATP6V1G3","ATP6V1C1","UBB","ATP6V0E2","ATP6V0C","ATP6V1G2","UBA52","ATP6V1E2","UQCRFS1","ATG7","SPTLC1","ATG4A","ATG14","ATG2A","MTOR","WIPI2","PIK3R4","ACACA","MECP2","PIK3C3","GTPBP2","ATG16L1","WIPI1","PRKAA1","RPTOR","RHEB","ATG101","SCP2","MLST8","ATG3","DEPTOR","STK11","TSC1","TSC2","DCAF17","ATG10","ULK1","AKT1S1","ATG12","ATG13","RB1CC1"]

MicIronGenes = [gene for gene in iron_genes if gene in ADMicroglia.var_names]
sc.tl.score_genes(ADMicroglia, MicIronGenes, score_name='iron_score', ctrl_size=min(50, len(MicIronGenes)*10))
with plt.rc_context({"figure.figsize": (8,4)}):
    sc.pl.violin(ADMicroglia, ['iron_score'], groupby='leiden0_7',save = '5_ADMicroglia_IronScores')

import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, permutation_test
from statsmodels.stats.multitest import multipletests

def analyze_iron_scores(ADMicroglia, score_col='iron_score', group_col='leiden0_7', n_perm=1000):
    """
    Comprehensive stats: Mann-Whitney (parametric p), Cohen's d. and FDR-corrected for multiple clusters.
    """
    clusters = ADMicroglia.obs[group_col].cat.categories
    results = []
    
    for cluster in clusters:
        mask = ADMicroglia.obs[group_col] == cluster
        scores_in = ADMicroglia.obs.loc[mask, score_col].values
        scores_out = ADMicroglia.obs.loc[~mask, score_col].values
        
        # Mann-Whitney (Wilcoxon rank-sum)
        mw_stat, mw_p = mannwhitneyu(scores_in, scores_out, alternative='two-sided')
        
        # Effect size (Cohen's d)
        mean_in, mean_out = np.mean(scores_in), np.mean(scores_out)
        pooled_std = np.sqrt((np.var(scores_in, ddof=1) + np.var(scores_out, ddof=1)) / 2)
        cohens_d = (mean_in - mean_out) / pooled_std if pooled_std > 0 else 0
        
        results.append({
            'cluster': cluster,
            'n_in': len(scores_in),
            'n_out': len(scores_out),
            'mean_in': mean_in,
            'mean_out': mean_out,
            'cohens_d': cohens_d,
            'mw_p': mw_p,
        })
    
    results_df = pd.DataFrame(results)
    
    # FDR correction
    _, results_df['mw_p_fdr'], _, _ = multipletests(results_df['mw_p'], method='fdr_bh')
    _, results_df['perm_p_fdr'], _, _ = multipletests(results_df['perm_p'], method='fdr_bh')
    
    # Significant: FDR<0.05 & |d|>0.5 (medium effect)
    results_df['significant'] = (results_df['mw_p_fdr'] < 0.05) & (np.abs(results_df['cohens_d']) > 0.5)
    
    return results_df.sort_values('mw_p_fdr')

# Run
stats_df = analyze_iron_scores(ADMicroglia)
print(stats_df.round(4))
print('\nIron-enriched clusters (|d|>0.5, FDR<0.05):\n', stats_df[stats_df.significant & (stats_df.cohens_d > 0)]['cluster'].tolist())

stats_df.to_csv("Leiden07_IronScores.csv", sep='\t', index=False, quotechar='"')

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

stats_df = stats_df.assign(
    cluster_num=stats_df['cluster'].astype(int)
).sort_values('cluster_num').drop(columns='cluster_num')

# Plot the iron scores in a heatmap
heatmap_data = stats_df.set_index('cluster')['cohens_d'].values.reshape(1, -1)

fig, ax = plt.subplots(figsize=(12, 3.5))
sns.heatmap(heatmap_data, 
            xticklabels=stats_df['cluster'], 
            yticklabels=['Iron Score'], 
            annot=True, fmt='.2f', cmap='RdBu_r', center=0, ax=ax,
            cbar_kws={'label': "Cohen's d"},
            linewidths=1, linecolor='white')
ax.grid(False)
plt.title('Iron Score Enrichment (Cohen\'s d): Positive = Enriched (3 & 4 highlighted)')
plt.xlabel('Leiden0_7 Clusters')
plt.ylabel('Signature')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

sc.pl.umap(ADMicroglia, color='iron_score', size=20, save='10_UMAPIronScores')
ADMicroglia.X = ADMicroglia.layers['raw_counts'].copy()
ADMicroglia.layers['scaled'] = sc.pp.scale(ADMicroglia, copy=True).X

fig, ax = plt.subplots(figsize=(8,4))
sns.boxplot(data=ADMicroglia.obs, x='leiden0_7', y='iron_score', palette=sns.color_palette())
plt.savefig('5a_ADMicroglia_IronScores.pdf')

# Senescence enrichment scoring
senmayo = ["ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2","BMP6","C3","CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26","CCL3","CCL4","CCL5","CCL7","CCL8","CD55","CD9","CSF1","CSF2","CSF2RB","CST4","CTNNB1","CTSB","CXCL1","CXCL10","CXCL12","CXCL16","CXCL2","CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR","EREG","ESM1","ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG","HGF","HMGB1","ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15","IL18","IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7","INHA","IQGAP2","ITGA2","ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13","MMP14","MMP2","MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT","PLAU","PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG","SEMA3F","SERPINB4","SERPINE1","SERPINE2","SPP1","SPX","TIMP2","TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B","TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2","CDC27","E2F2","SCMH1","MRE11","MAP2K3","MAPK9","ANAPC4","MAP2K4","MAP4K4","RPS6KA2","UBE2D1","EED","MAP2K7","TNRC6C","MAPKAPK5","ANAPC5","TNRC6A","TINF2","AGO1","CDC23","CABIN1","MAPK1","HIRA","TNRC6B","E2F1","RBBP7","MAPK3","ACD","NBN","CCNE1","FZR1","ERF","CDK6","H2AZ2","EZH2","MAPK8","UBE2S","MAP2K6","NFKB1","MAPK10","ANAPC15","CDKN1B","PHC1","ASF1A","MAPK14","E2F3","LMNB1","RAD50","TFDP2","MAPKAPK3","RPS6KA1","UBN1","RNF2","CDKN2C","CDK2","H1-3","H1-1","H2BC11","CDKN1A","ID1","AGO3","POT1","CDKN2D","CDC16","H3-3B","KDM6B","TERF2","CCNA1","PHC2","AGO4","ETS1","CDK4","MDM2","TXN","HMGA1","RB1","MINK1","TP53","ANAPC11","CBX8","CBX4","RPS27A","CCNA2","H2BC1","TERF1","CDKN2B","CDKN2A","ATM","HMGA2","UBC","VENTX","ANAPC1","TNIK","MOV10","H2BC5","H4C8","RBBP4","MAPKAPK2","H3-3A","ANAPC10","ANAPC16","MAPK7","TERF2IP","H3-4","BMI1","H1-4","STAT3","UBE2E1","UBB","FOS","IFNB1","CEBPB","KAT5","RELA","PHC3","CBX2","UBE2C","CCNE2","ANAPC2","CDC26","RPS6KA3","SUZ12","H2AC6","H2BC4","EHMT1","EP400","H3C13","CBX6","H2AC20","H1-5","H2BC21","H2BC13","MAPK11","SP1","H1-2","H2AX","H1-0","ANAPC7","H2AC7","H2BC26","H4C3","H3C12","H4C11","H3C4","MAP3K5","H4C16","H2BC12","TFDP1","MDM4","H3C14","H3C15","RING1","EHMT2","UBA52","H2BC15","H2BC12L","H2AJ","H4C15","H4C14","H4C12","H2BC14","H2BC8","H3C8","H2AB1","H2BC6","H4C6","H2BC17","H3C6","H4C13","H3C11","H2BC9","H3C1","H4C9","H2AC14","H2BC3","H4C5","H2AC8","H4C4","H2BC7","H3C7","H2AC4","H2BC10","H4C1","H4C2","H3C10","MIR24-2","MIR24-1","H3C2","H3C3","H2AC18","H2AC19"]

MicSenMayoGenes = [gene for gene in senmayo if gene in ADMicroglia.var_names]
sc.tl.score_genes(ADMicroglia, MicSenMayoGenes, score_name='SenMayo_score', ctrl_size=min(50, len(MicSenMayoGenes)*10))
with plt.rc_context({"figure.figsize": (10,8)}):
    sc.pl.violin(ADMicroglia, ['SenMayo_score'], groupby='leiden0_7',save = '15_ADMicroglia_SenMayoScores')

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, permutation_test
from statsmodels.stats.multitest import multipletests

def analyze_sen_scores(ADMicroglia, score_col='SenMayo_score', group_col='leiden0_7', n_perm=1000):
    """
    Comprehensive stats: Mann-Whitney (parametric p), Cohen's d. and FDR-corrected for multiple clusters.
    """
    clusters = ADMicroglia.obs[group_col].cat.categories
    results = []
    
    for cluster in clusters:
        mask = ADMicroglia.obs[group_col] == cluster
        scores_in = ADMicroglia.obs.loc[mask, score_col].values
        scores_out = ADMicroglia.obs.loc[~mask, score_col].values
        
        # Mann-Whitney (Wilcoxon rank-sum)
        mw_stat, mw_p = mannwhitneyu(scores_in, scores_out, alternative='two-sided')
        
        # Effect size (Cohen's d)
        mean_in, mean_out = np.mean(scores_in), np.mean(scores_out)
        pooled_std = np.sqrt((np.var(scores_in, ddof=1) + np.var(scores_out, ddof=1)) / 2)
        cohens_d = (mean_in - mean_out) / pooled_std if pooled_std > 0 else 0
        
        results.append({
            'cluster': cluster,
            'n_in': len(scores_in),
            'n_out': len(scores_out),
            'mean_in': mean_in,
            'mean_out': mean_out,
            'cohens_d': cohens_d,
            'mw_p': mw_p,
            'perm_p': perm_res.pvalue
        })
    
    results_df = pd.DataFrame(results)
    
    # FDR correction
    _, results_df['mw_p_fdr'], _, _ = multipletests(results_df['mw_p'], method='fdr_bh')
    _, results_df['perm_p_fdr'], _, _ = multipletests(results_df['perm_p'], method='fdr_bh')
    
    # Significant: FDR<0.05 & |d|>0.5 (medium effect)
    results_df['significant'] = (results_df['mw_p_fdr'] < 0.05) & (np.abs(results_df['cohens_d']) > 0.5)
    
    return results_df.sort_values('mw_p_fdr')

# Run
stats_df = analyze_sen_scores(ADMicroglia)
print(stats_df.round(4))
print('\nsen-enriched clusters (|d|>0.5, FDR<0.05):\n', stats_df[stats_df.significant & (stats_df.cohens_d > 0)]['cluster'].tolist())

stats_df.to_csv("leiden07_senMayoScores.csv", sep='\t', index=False, quotechar='"')

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

stats_df = stats_df.assign(
    cluster_num=stats_df['cluster'].astype(int)
).sort_values('cluster_num').drop(columns='cluster_num')

# plotting the enrichment heatmap for senmayo enrichment scores
heatmap_data = stats_df.set_index('cluster')['cohens_d'].values.reshape(1, -1)

fig, ax = plt.subplots(figsize=(12, 3.5))
sns.heatmap(heatmap_data, 
            xticklabels=stats_df['cluster'], 
            yticklabels=['SenMayo Score'], 
            annot=True, fmt='.2f', cmap='RdBu_r', center=0, ax=ax,
            cbar_kws={'label': "Cohen's d"},
            linewidths=1, linecolor='white')
ax.grid(False)
plt.title('SenMayo Score Enrichment (Cohen\'s d): Positive = Enriched (3 highlighted)')
plt.xlabel('Leiden0_7 Clusters')
plt.ylabel('Signature')
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()
