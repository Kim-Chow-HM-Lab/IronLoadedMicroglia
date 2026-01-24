#### 1-Processing of Individual Datasets analysed ####
# Analysis was performed on Python 3.10.19

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

#### Diagrams parameters ####
import pathLib
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

#### Example Preprocessing Code using the snRNA-seq dataset available from Lau et al, 2020 ####
# All individual datasets were processed in the same manner, hence an exemplery code is shown to demonstrate. The data is from Lau et al, 2020 (https://pubmed.ncbi.nlm.nih.gov/32989152/; PMID: 32989152) available on GEO Omnibus at GSE175827. 
main_directory = "~/Lau2020/rawData"
data_dirs = {}
for sample_dir in os.listdir(main_directory):
    if sample_dir.startswith("GSM"):
        sample_path = os.path.join(main_directory, sample_dir)
        required_files = {"barcodes.tsv", "genes.tsv", "matrix.mtx"}
        files_in_dir = set(os.listdir(sample_path))
        if required_files.issubset(files_in_dir):
            data_dirs[sample_dir] = sample_path
        else:
            missing_files = required_files - files_in_dir
            print(f"Warning: Missing files {missing_files} in {sample_dir}")
print(f"Found {len(data_dirs)} valid samples out of {len(os.listdir(main_directory))} directories")
adatas = {}
for sample_name, sample_path in data_dirs.items():
    adata = sc.read_10x_mtx(
        sample_path,
        var_names = 'gene_symbols',
        cache=True
    )
    adata.obs['sample'] = sample_name
    adatas[sample_name] = adata
Lau = sc.concat(adatas.values(), label = 'sample', keys = list(adatas.keys()))
metadata = pd.read_table("~/Lau2020/sampleInfo.txt", sep=",")
metadata = metadata.set_index("Accession")
Lau.obs = Lau.obs.merge(metadata, left_on = 'sample', right_index = True, how = 'left')
Lau
Lau.write_h5ad("~/Lau2020/H5AD/Lau2020_Raw.h5ad")

#### 1. Filtering ####
Lau.var['mt'] = Lau.var_names.str.startswith('MT-')
Lau.var['ribo'] = Lau.var_names.str.startswith(('RPS','RPL'))
Lau.var['hb'] = Lau.var_names.str.contains('^HB[^(P)]')
sc.pp.calculate_qc_metrics(Lau, qc_vars = ['mt','ribo','hb'], inplace = True, percent_top=[20], log1p=True)

# QC Plots
p1 = sns.displot(Lau.obs['total_counts'], bins=100,kde=False)
save_plot(p1.figure, '1_Lau2020_Totalcounts')
sc.pl.violin(Lau, 'total_counts', save = '2_Lau2020_VlnTotalCounts')
sc.pl.violin(Lau, "pct_counts_mt", save = '3_Lau2020_VlnPCTMT')
sc.pl.scatter(Lau, 'n_genes_by_counts', 'total_counts', color='pct_counts_mt', save='4_Lau2020_PCTScatter')
sc.pl.violin(
    Lau,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
    jitter=0.4,
    multi_panel=True,
    rotation=0,
    save = '5_Lau2020_QC'
)

# Here we write a function to define the outliers, taking a metric (column in the .obs) and the number of median absolute deviations (MADs) that is still permissive within the filtering strategy. We filter out the counts within a threshold of 5 MADs to the 'log1p_total_counts', 'log1p_n_genes_by_counts' and 'pct_counts_in_top_20_genes'. 
def is_outlier(Lau,  metric: str, nmads: int):
    M = Lau.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) +  nmads * median_abs_deviation(M) <  M
    )
    return outlier

Lau.obs['outlier'] = (
    is_outlier(Lau, "log1p_total_counts",  5)
    | is_outlier(Lau, "log1p_n_genes_by_counts", 5)
    | is_outlier(Lau, "pct_counts_in_top_20_genes", 5)
)
Lau.obs.outlier.value_counts()
#outlier
#False    169792
#True       9600
#Name: count, dtype: int64

# 'pct_counts_Mt' is  filtered with 3 MADs. Additionally, cells with a percetange of mitochondrial counts exceeding 20% are filtered out.
Lau.obs['mt_outlier'] = is_outlier(Lau, 'pct_counts_mt', 3) | (Lau.obs['pct_counts_mt'] > 20)
Lau.obs.mt_outlier.value_counts()
print(f"Total number of cells: {Lau.n_obs}")
Lau = Lau[(~Lau.obs.outlier) & (~Lau.obs.mt_outlier)].copy()
print(f"Number of cells after filtering: {Lau.n_obs}")
#Total number of cells:  179392
#Number of cells after filtering: 141401
#78% of cells remain.

# QC Plots after filtering
p5 = sns.displot(Lau.obs['total_counts'], bins=100,kde=False)
save_plot(p1.figure, '6_Lau2020_FilteredTotalcounts')
sc.pl.violin(Lau, 'total_counts', save = '7_Lau2020_FilteredVlnTotalCounts')
sc.pl.violin(Lau, "pct_counts_mt", save = '8_Lau2020_FilteredVlnPCTMT')
sc.pl.scatter(Lau, 'n_genes_by_counts', 'total_counts', color='pct_counts_mt', save='9_Lau2020_FilteredPCTScatter')
sc.pl.violin(
    Lau,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
    jitter=0.4,
    multi_panel=True,
    rotation=0,
    save = '10_Sun2023_QCFiltered'
)

#### 2. Doublet Detection ####
# This was performed using the rpy2 package to convert anndata objects to use R interchangeably within the python environment to perform doublet detection using scDblFinder. 
```{python}
%%R
library(Seurat)
library(scater)
library(scDblFinder)
library(BiocParallel)
```

subject_counts = Lau.obs['sample'].value_counts()
valid_subjects = subject_counts[subject_counts >= 10].index.tolist()
Lau = Lau[Lau.obs['sample'].isin(valid_subjects)].copy()
Lau
Lau.layers['raw_counts'] = Lau.X.copy()
data_mat = Lau.X.T
data_mat = data_mat.toarray()
meta = Lau.obs

```{python}
%%R -i data_mat -i meta -o doublet_score -o doublet_class
set.seed(123)
sce = SingleCellExperiment(
    list(counts=data_mat),
    colData = meta
)
sce = scDblFinder(sce, samples = 'sample', verbose=TRUE, BPPARAM=MulticoreParam(62))
doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class 
```

Lau.obs['scDblFinder_score'] = doublet_score
Lau.obs['scDblFinder_class'] = doublet_class
Lau.obs.scDblFinder_class.value_counts()
#scDblFinder_class
#singlet    126579
#doublet     14819
#Name: count, dtype: int64

#### 3. Normalization ####
scale_counts = sc.pp.normalize_total(Lau, target_sum = None, inplace=False)
Lau.layers['shifted_log'] = sc.pp.log1p(scale_counts['X'], copy=True)
fig, axes = plt.subplots(1,2,figsize=(10,5))
p1 = sns.histplot(Lau.obs['total_counts'], bins = 100, kde=False,ax=axes[0])
axes[0].set_title('Total Counts')
p2 = sns.histplot(Lau.layers['shifted_log'].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title('Shifted Logarithm Normalization')
plt.show()
save_plot(fig, '11_Lau2020_ShiftedLogNormalized')

#### 4. Feature Selection using scry ####
Lau.X = Lau.layers['raw_counts'].copy()

```{python}
%%R
library(scry)
```

ro.globalenv['Lau'] = Lau

```{python}
%%R
sce = devianceFeatureSelection(Lau, assay = 'X')
```

binomial_deviance = ro.r('rowData(sce)$binomial_deviance').T
idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(Lau.var_names.shape, dtype = bool)
mask[idx] = True
Lau.var['highly_deviant'] = mask
Lau.var['binomial_deviance'] = binomial_deviance
sc.pp.highly_variable_genes(Lau, layer = "shifted_log")

ax = sns.scatterplot(
    data = Lau.var, x="means", y="dispersions",hue="highly_deviant", s=5
)
ax.set_xlim(None, 5)
ax.set_ylim(None, 5)
plt.show()

#### 5. Dimensionality Reduction and Clustering ####
# PCA
Lau.X = Lau.layers['shifted_log']
Lau.var['highly_variable'] = Lau.var['highly_deviant']
sc.pp.pca(Lau, svd_solver = 'arpack', use_highly_variable = True)
sc.pl.pca_scatter(Lau, color = 'total_counts', frameon = True, save = '13_Lau2020_PCA')

# UMAP
sc.pp.neighbors(Lau)
sc.tl.umap(Lau)
sc.pl.umap(Lau, color = 'total_counts', frameon=True, save = '14_Lau2023_UMAP')

# Clustering
sc.pl.pca_variance_ratio(Lau, n_pcs=50, log=True, save = "15_Lau2020_Variance")

# clustree analysis using pyclustree
from pyclustree import clustree

resolutions = np.arange(0.1, 1.5, 0.1).round(1)
for resolution in resolutions:
    sc.tl.leiden(
        Lau,
        resolution = resolution,
        flavor = 'leidenalg',
        n_iterations=2,
        key_added=f"leiden{str(resolution).replace('.','_')}"
    )

fig = clustree(
    Lau,
    [f"leiden{str(resolution).replace('.','_')}" for resolution in resolutions],
    title = 'Clustree of Lau2020',
    edge_weight_threshold = 0.00,
    show_fraction=True,
    x_spacing = 25
)
fig.set_size_inches(15,10)
fig.set_dpi(100)
save_plot(fig, "16_Lau2020_ClusTree")

#### 6. Automated Cell Annotation using CellTypist ####
import celltypist
from celltypist import models
models.download_models(force_update=True, model="Adult_Human_PrefrontalCortex.pkl")
model_PFC = models.Model.load(model="Adult_Human_PrefrontalCortex.pkl")
model_PFC.cell_types
Lau.obs_names_make_unique()
Lau.var_names_make_unique()
Lau.X = Lau.layers['raw_counts'].copy()
sc.pp.normalize_total(Lau, target_sum=1e4)
sc.pp.log1p(Lau)
predictions = celltypist.annotate(Lau, model=model_PFC, majority_voting=True)
predictions_adata = predictions.to_adata()
Lau.obs['CT_Predictions'] = predictions_adata.obs['predicted_labels']
Lau.obs['CT_overclustering'] = predictions_adata.obs['over_clustering']
Lau.obs['CT_majority_vote']  = predictions_adata.obs['majority_voting']
sc.pl.umap(Lau, color = ['CT_majority_vote'], frameon= True, sort_order = False, wspace=0.4, save =  '17_Lau2020_CellTypist')
Lau.obs['CT_majority_vote'].value_counts()

# Validation using identified cell markers for major cell type groups
marker_genes = {
    'Excitatory': ['SLC17A7','SATB2','NPTX1','CDH22','SLC17A6','CDH13','CHRM2','TRPC6','HPCAL1','CBLN4','NPY','RORB','HTR2C','FEZF2','SYT6','CTGF','CPLX3','OX2R'],
    'Inhibitory':['GAD1','GAD2','ERBB4','LHX6','ADARB2'],
    'Oligodendrocytes': ['SOX10','OLIG2','MBP','PLP1','MAG'],
    'Endothelial': ['PECAM1','FLT1','VWF','EPAS1','PDGFRB'],
    'Astrocytes': ['GFAP','ADH1L1','GJA1','AQP4'],
    'Microglia': ['CSF1R','P2RY12','TMEM119','SALL1'],
    'OPC': ['NG2','PDGFRA','VCAN','PCDH15']
}
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in Lau.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

Lau.X = Lau.layers['shifted_log'].copy()

sc.pl.dotplot(
    Lau,
    groupby="CT_majority_vote",
    var_names = marker_genes_in_data,
    standard_scale = 'var',
    save = '18_Lau2020_CellTypeValidate'
)

conditions = [
    Lau.obs['CT_majority_vote'].str.contains('Astro', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('COP', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('Endo', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('InN', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('L2-3', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('L2-4', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('L3-5', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('L5', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('L5-6', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('L6', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('Macro', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('Micro', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('OPC', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('Oligo', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('PC', case=False, na=False),
    Lau.obs['CT_majority_vote'].str.contains('VLMC', case=False, na=False),
]


choices = [
    'Astrocytes',
    'OPC',
    'Endothelial',
    'Inhibitory',
    'Excitatory',
    'Excitatory',
    'Excitatory',
    'Excitatory',
    'Excitatory',
    'Excitatory',
    'Microglia',
    'Microglia',
    'OPC',
    'Oligo',
    'Endothelial',
    'Vascular'
]

Lau.obs['celltypes'] = np.select(conditions,choices,default = Lau.obs['CT_majority_vote'])


Lau.obs['celltypes'].value_counts()
#celltypes
#Excitatory     45440
#Oligo          41831
#Inhibitory     17912
#Astrocytes     15567
#OPC            11984
#Microglia       6009
#Endothelial     2470
#Vascular         185
#Name: count, dtype: int64

sc.pl.umap(Lau, color = 'celltypes', frameon=False, sort_order = False, legend_fontsize=6, save='19_Lau2020_CellTypes')
ax = sc.pl.umap(Lau, color=['celltypes'], groups = ['Microglia'], show=False, save="20_Lau2020_MicrogliaUMAP")
legend_texts = ax.get_legend().get_texts()
for legend_text in legend_texts:
    if legend_text.get_text() == "NA":
        legend_text.set_text("Other Cell Types")

Lau_Microglia = Lau[Lau.obs['celltypes'] == 'Microglia']
Lau_Microglia
Lau_Microglia.write_h5ad("~/Lau2020/H5AD/Lau2020_Microglia.h5ad")







