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

# Loading the processed anndata for microglial nuclei 
Lau2020 = sc.read_h5ad("~/Merged/Finalised/Lau2020_PFCMicrogliaFinal.h5ad")
Lau2020
Dharshini2025 = sc.read_h5ad("~/Merged/Finalised/Dharshini2025_PFCMicrogliaFinal.h5ad")
Dharshini2025
Liu2025 = sc.read_h5ad("~/Merged/Finalised/Liu2025_PFCMicrogliaFinal.h5ad")
Liu2025
Mathys2019 = sc.read_h5ad("~/Merged/Finalised/Mathys2019_PFCMicrogliaFinal.h5ad")
Mathys2019
Mathys2023 = sc.read_h5ad("~/Merged/Finalised/Mathys2023_PFCMicrogliaFinal.h5ad")
Mathys2023
Morabito2021 = sc.read_h5ad("~/Merged/Finalised/Morabito2021_PFCMicrogliaFinal.h5ad")
Morabito2021
Sun2023 = sc.read_h5ad("~/Merged/Finalised/Sun2023_PFCMicrogliaFinal.h5ad")
Sun2023

# Concatenation of the individual datasets
adatas = [Mathys2019, Lau2020, Morabito2021, Sun2023, Mathys2023, Liu2025, Dharshini2025]
for i, adata in enumerate(adatas):
    dups = adata.var_names[adata.var_names.duplicated()].unique()
    if len(dups) > 0:
        print(f"Dataset {i+1} has {len(dups)} duplicate genes: {dups[:5]}...")

ADMicroglia = ad.concat(
    adatas,
    axis=0,
    join="outer",
    index_unique=None,
)
ADMicroglia

# Unintegrated data
ADMicroglia.obs['APOEStatus'].value_counts()
ADMicroglia.obs['Dataset'].value_counts()
ADMicroglia.obs['Sex'].value_counts()
ADMicroglia.obs['Disease'].value_counts()
ADMicroglia.X = ADMicroglia.layers['raw_counts'].copy()
sc.pp.filter_genes(ADMicroglia, min_cells=1)
ADMicroglia

ADMicroglia.obs['Dataset'] = pd.Categorical(ADMicroglia.obs['Dataset'], categories = ['Mathys2019', 'Lau2020', 'Morabito2021', 'Sun2023', 'Mathys2023', 'Liu2025', 'Dharshini2025'], ordered=True)

# Quality Control Plots
sc.pl.violin(ADMicroglia,'total_counts',groupby = 'Dataset', rotation = 45)
plt.gca().invert_yaxis()
plt.show()
ADMicroglia.X = ADMicroglia.layers['raw_counts'].copy()
sc.pp.normalize_total(ADMicroglia)
sc.pp.log1p(ADMicroglia)
ADMicroglia.layers['log_counts'] = ADMicroglia.X.copy()
sc.pp.highly_variable_genes(ADMicroglia)
sc.tl.pca(ADMicroglia)
sc.pp.neighbors(ADMicroglia)
sc.tl.umap(ADMicroglia)
ADMicroglia.obsm['X_umap_unintegrated'] = ADMicroglia.obsm['X_umap'].copy()
ADMicroglia
sc.pl.umap(ADMicroglia, color='Dataset', size=5,frameon=True, title="ADMicroglia: UMAP Unintegrated",save="2_ADMic_UMAPUnin_Dataset")

# Batch aware feature selection
# We perform batch-aware highly variable gene selection by setting the batch_key argument in the scanpy function, which will then calculate the HVGs for each batch separately and combine the results by selecting those genes that are highly variable in the highest number of batches.
ADMicroglia.X = ADMicroglia.layers['raw_counts'].copy()
sc.pp.highly_variable_genes(
    ADMicroglia,
    n_top_genes=4000,
    flavor="seurat_v3",
    layer="raw_counts",
    subset=False,
    batch_key="Dataset",
)
ADMicroglia.var
n_batches=ADMicroglia.var['highly_variable_nbatches'].value_counts()
ax=n_batches.plot(kind="bar")
n_batches

# Integration
# We assume that all the data is completely unlabelled and we wish to find common axes of variation between the two datasets. A few integration methods will be used and the integration metrics will be examined using scib-metrics after completion. 

# scVI
import os
import scvi
import torch
from rich import print
from scib_metrics.benchmark import Benchmarker
scvi.settings.seed=0
print("Last run with scvi-tools version:", scvi.__version__)
for col in ["Sex", "Sample", "PMI"]:
    if col in ADMicroglia.obs.columns and ADMicroglia.obs[col].dtype.name == "category":
        values = ADMicroglia.obs[col].astype(str)
        ADMicroglia.obs[col] = pd.Categorical(values, categories=sorted(values.unique()))
      
scvi.model.SCVI.setup_anndata(ADMicroglia, layer="raw_counts", batch_key="Dataset", categorical_covariate_keys=["Sex","Sample","PMI"])
model_scvi = scvi.model.SCVI(ADMicroglia, n_layers=2, n_latent=30, gene_likelihood="nb", dispersion="gene-batch")
torch.cuda.is_available()
torch.set_float32_matmul_precision('highest')
model_scvi.train(accelerator="gpu")
ADMicroglia.obsm['X_scVI'] = model_scvi.get_latent_representation()
sc.pp.neighbors(ADMicroglia, use_rep="X_scVI")

# pyclustree
from pyclustree import clustree

resolutions = np.arange(0.1, 1.5,0.1).round(1)
for resolution in resolutions:
    sc.tl.leiden(
        ADMicroglia,
        resolution = resolution,
        flavor = 'leidenalg',
        n_iterations=2,
        key_added=f"leiden{str(resolution).replace('.','_')}",
    )

fig = clustree(
    ADMicroglia,
    [f"leiden{str(resolution).replace('.','_')}" for resolution in resolutions],
    title='Clustree of All AD Microglia',
    edge_weight_threshold=0.00,
    show_fraction=False,
    x_spacing=50,      
    y_spacing=6,     
)

fig.set_size_inches(14, 18)  
fig.set_dpi(300)
save_plot(fig, "4_ADMicroglia_pyClusTree")

sc.tl.umap(ADMicroglia, min_dist=0.3)
ADMicroglia.obsm['X_scVI_umap'] = ADMicroglia.obsm['X_umap'].copy()
sc.pl.umap(ADMicroglia, color='Dataset', size=5,frameon=True, save="5a_UMAP_scVI_ADMicroglia", title="UMAP of ADMicroglia: scVI")

# scANVI
model_scanvi = scvi.model.SCANVI.from_scvi_model(model_scvi, labels_key='leiden0_7', unlabeled_category='unlabelled')
print(model_scanvi)
model_scanvi.view_anndata_setup()
model_scanvi.train(max_epochs=20, accelerator='gpu', n_samples_per_label=100)
ADMicroglia.obsm['X_scANVI'] = model_scanvi.get_latent_representation()

sc.pp.neighbors(ADMicroglia, use_rep="X_scANVI")
sc.tl.umap(ADMicroglia, min_dist = 0.3)
ADMicroglia.obsm['X_scANVI_umap'] = ADMicroglia.obsm['X_umap'].copy()

sc.pl.umap(
    ADMicroglia,
    color = ['leiden0_7'],
    size=5,
    title = 'UMAP: ADMicroglia, scANVI',
    save='16_UMAP_scANVI'
)

sc.pl.umap(
    ADMicroglia,
    color = ['Dataset'],
    size=5,
    title = 'UMAP: Datasets, scANVI',
    save='17_UMAP_scANVIDatasets'
)

# BBKNN
# Batch balanced kNN alters the kNN procedure to identify each cell's top neighbours in each batch separately instead of the entire cell pool with no accounting for batch. The nearest neighbours for each batch are then merged to create a final list of neighbours for the cells, aligning batches in a quick and lightweight manner. 
sc.tl.pca(ADMicroglia)
import bbknn
sc.external.pp.bbknn(
    ADMicroglia,
    batch_key="Dataset",
    use_rep = "X_pca"
)
sc.tl.umap(ADMicroglia)
ADMicroglia.obsm['X_bbknn'] = ADMicroglia.obsm['X_umap'].copy()
ADMicroglia

sc.pl.umap(
    ADMicroglia,
    color = ['Dataset'],
    size=5,
    title = 'Datasets, BBKNN',
    save='18_UMAP_BBKNN_Datasets'
)

# Harmony
import scanpy.external as sce
import harmonypy as hm
sc.pp.pca(ADMicroglia)
data_mat = ADMicroglia.obsm['X_pca'].T
meta_data = pd.DataFrame(index=ADMicroglia.obs.index, data=ADMicroglia.obs[['Dataset']])
ho = hm.run_harmony(
    data_mat, meta_data,
    vars_use = ['Dataset'],
    max_iter_harmony = 10,
    max_iter_kmeans=20,
    epsilon_cluster=1e-05,
    epsilon_harmony=0.0001,
    nclust=100,
    block_size=0.05,
    lamb=1.0,
    theta=2,
    sigma=0.1,
    verbose=True,
    random_state=0
)
ADMicroglia.obsm['X_pca_harmony'] = ho.Z_corr
ADMicroglia

sc.pp.neighbors(ADMicroglia, use_rep = 'X_pca_harmony', n_pcs=50)
sc.tl.umap(ADMicroglia)
ADMicroglia.obsm['X_umap_harmony'] = ADMicroglia.obsm['X_umap'].copy()
sc.pl.umap(
    ADMicroglia,
    color = ['Dataset'],
    size=5,
    title = 'Datasets: Harmony',
    save='19_UMAP_Harmony_Datasets'
)

# Compute integration metrics
bm = Benchmarker(
    ADMicroglia,
    batch_key="Dataset", 
    label_key = "leiden0_7",
    embedding_obsm_keys=["X_pca", 'X_umap_unintegrated',"X_scVI", "X_scANVI", 'X_bbknn', 'X_umap_harmony'],
    n_jobs=-1,
)
bm.benchmark()
bm.plot_results_table(min_max_scale=False)

df = bm.get_results(min_max_scale=False)
print(df)
df.to_csv("SCIBMetrics.csv", sep='\t', index=False, quotechar='"')

