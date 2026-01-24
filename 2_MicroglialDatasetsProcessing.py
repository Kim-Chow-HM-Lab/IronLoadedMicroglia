#### 2-Processing of microglial fractions extracted from individual datasets ####
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

#### Mathys et al, 2019 ####
Mathys2019 = sc.read_h5ad("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/Merged/H5AD/Mathys2019_Microglia.h5ad")
Mathys2019.obs['Dataset'] = "Mathys2019"
Mathys2019
Mathys2019_meta = Mathys2019.obs
Mathys2019_meta.to_csv("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/Merged/H5AD/Mathys2019_meta.csv", index=True)
Mathys2019.obs.rename(
    columns={
        "projid": "Sample",
        "msex": "Sex",
        "apoe_genotype": "APOEStatus",
        "age_death": "Age",
        "pmi": "PMI",
        "braaksc": "BraakStage",
        "ceradsc": "CERADScore",
    },
    inplace=True
)
Mathys2019.obs['Sex'].value_counts()
mapping_gender = {
    0.0: "Female",
    1.0: "Male",
}
Mathys2019.obs['Sex'] = Mathys2019.obs['Sex'].map(mapping_gender)
Mathys2019.obs['Sex'].value_counts()
ROSMAP_Clinical = pd.read_table("~/MetaData/clinical.status.txt")
Mathys2019_Meta = Mathys2019.obs.copy()
Mathys2019_Meta = Mathys2019_Meta.reset_index()
ROSMAP_Clinical = ROSMAP_Clinical.reset_index()
ROSMAP_Status = ROSMAP_Clinical[["index","Status","individualID"]]
Mathys2019_Merged = pd.merge(Mathys2019_Meta, ROSMAP_Status, left_on="individualID", right_on="individualID", how='left')
Mathys2019_Merged = Mathys2019_Merged.set_index("TAG")
if list(Mathys2019_Merged.index) == list(Mathys2019.obs_names):
    Mathys2019.obs = Mathys2019_Merged
Mathys2019.obs['Status'].value_counts()
Mathys2019.obs["Disease"] = Mathys2019.obs["Status"].map({
    'Dementia_Two Pathology': 'AD',
    'Normal': 'NC',
    'Non Dementia_Neurofibrillary Tangle': 'NC',
    'Non Dementia_Two Pathology': 'NC',
    'Dementia_No Pathology': 'AD',
    'Dementia_Neurofibrillary Tangle': 'AD',
    'Non Dementia_Neuritic Plaques': 'NC'
})
Mathys2019.obs['Disease'].value_counts()
Mathys2019.obs["APOEStatus"].value_counts()
Mathys2019.obs['APOEStatus'] = Mathys2019.obs['APOEStatus'].map(mapping)
Mathys2019.obs['APOEStatus'].value_counts()
Mathys2019
columns_to_remove = ['leiden0_1', 'leiden0_2', 'leiden0_3', 'leiden0_4', 'leiden0_5', 'leiden0_6', 'leiden0_7', 'leiden0_8', 'leiden0_9', 'leiden1_0', 'leiden1_1', 'leiden1_2', 'leiden1_3', 'leiden1_4', 'predicted_labels', 'over_clustering', 'majority_voting', 'conf_score', 'CT_Predictions', 'CT_overclustering', 'CT_majority_vote', 'celltypes']
Mathys2019.obs = Mathys2019.obs.drop(columns = columns_to_remove)
del Mathys2019.uns
del Mathys2019.obsm
del Mathys2019.varm
del Mathys2019.layers['shifted_log']
del Mathys2019.obsp

Mathys2019.write_h5ad("~/Merged/Finalised/Mathys2019_PFCMicrogliaFinal.h5ad")

#### Lau et al, 2020 ####
Lau2020 = sc.read_h5ad("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/Merged/H5AD/Lau2020_Microglia.h5ad")
Lau2020.obs['Dataset'] = "Lau2020"
Lau2020
Lau2020_meta = Lau2020.obs
Lau2020_meta.to_csv("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/Merged/H5AD/Lau2020_meta.csv", index=True)
Lau2020
Lau2020.obs.rename(
    columns={
        "sample": "Sample",
        "pmi": "PMI",
        "APOE": "APOEStatus"
    },
    inplace=True
)
Lau2020.obs['Sex'].value_counts()
Lau2020.obs['Diagnosis'].value_counts()
Lau2020.obs['Disease'] = Lau2020.obs['Diagnosis'].map({"Control":"NC","Alzheimer":"AD"})
Lau2020.obs['Disease'].value_counts()
Lau2020.obs['APOEStatus'].value_counts()
mapping = {
    33.0: "E3/E3",
    34.0: "E3/E4",
    23.0: "E2/E3",
    24.0: "E2/E4",
    44.0: "E4/E4",
}

Lau2020.obs["APOEStatus"] = Lau2020.obs["APOEStatus"].map(mapping)
Lau2020.obs["APOEStatus"].value_counts()
Lau2020
columns_to_remove = ['leiden0_1', 'leiden0_2', 'leiden0_3', 'leiden0_4', 'leiden0_5', 'leiden0_6', 'leiden0_7', 'leiden0_8', 'leiden0_9', 'leiden1_0', 'leiden1_1', 'leiden1_2', 'leiden1_3', 'leiden1_4', 'predicted_labels', 'over_clustering', 'majority_voting', 'conf_score', 'CT_Predictions', 'CT_overclustering', 'CT_majority_vote', 'celltypes']
Lau2020.obs = Lau2020.obs.drop(columns=columns_to_remove)
del Lau2020.uns
del Lau2020.obsm
del Lau2020.varm
del Lau2020.obsp
del Lau2020.layers['shifted_log']
Lau2020.write_h5ad("~/Merged/Finalised/Lau2020_PFCMicrogliaFinal.h5ad")

#### Morabito et al, 2021 ####
Morabito2021 = sc.read_h5ad("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/Merged/H5AD/Morabito2021_Microglia.h5ad")
Morabito2021.obs['Dataset'] = "Morabito2021"
Morabito2021
Morabito2021_meta = Morabito2021.obs
Morabito2021_meta.to_csv("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/Merged/H5AD/Morabito2021_meta.csv", index=True)
Morabito2021.obs['Diagnosis'].value_counts()
Morabito2021.obs['Disease'] = Morabito2021.obs['Diagnosis'].map({
    'AD': 'AD',
    'Control': 'NC'
})
Morabito2021.obs['Disease'].value_counts()
Morabito2021.obs.rename(
    columns={
        'SampleID': 'Sample',
    }, inplace=True
)
Morabito2021
columns_to_remove = ['leiden0_1', 'leiden0_2', 'leiden0_3', 'leiden0_4', 'leiden0_5', 'leiden0_6', 'leiden0_7', 'leiden0_8', 'leiden0_9', 'leiden1_0', 'leiden1_1', 'leiden1_2', 'leiden1_3', 'leiden1_4', 'predicted_labels', 'over_clustering', 'majority_voting', 'conf_score', 'CT_Predictions', 'CT_overclustering', 'CT_majority_vote', 'celltypes']
Morabito2021.obs = Morabito2021.obs.drop(columns = columns_to_remove)
del Morabito2021.uns
del Morabito2021.obsm
del Morabito2021.varm
del Morabito2021.layers['shifted_log']
del Morabito2021.obsp
Morabito2021.obs['Sex'].value_counts()
Morabito2021.obs['Sex'] = Morabito2021.obs['Sex'].map({
    'M': 'Male',
    'F': 'Female'
})
Morabito2021.write_h5ad("~/Merged/Finalised/Morabito2021_PFCMicrogliaFinal.h5ad")


#### Mathys et al, 2023 ####
Mathys2023 = sc.read_h5ad("~/Merged/H5AD/Mathys2023_Microglia.h5ad")
Mathys2023.obs['Dataset'] = "Mathys2023"
Mathys2023
Mathys2023_meta = Mathys2023.obs
Mathys2023_meta.to_csv("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/Merged/H5AD/Mathys2023_meta.csv", index=True)
Mathys2023.obs.rename(
    columns={
        "projid": "Sample",
        "apoe_genotype": "APOEStatus",
        "Gender": "Sex",
        "pmi": "PMI",
        "braaksc": "BraakStage",
        "ceradsc": "CERADScore",
        "age_death": "Age",
    },
    inplace=True
)
Mathys2023.obs['Sex'].value_counts()
Mathys2023.obs['Status'].value_counts()
Mathys2023_Meta = Mathys2023.obs.copy()
Mathys2023_Meta = Mathys2023_Meta.reset_index()
ROSMAP_Clinical = ROSMAP_Clinical.reset_index()
ROSMAP_Status = ROSMAP_Clinical[["index","Status","individualID"]]
Mathys2023_Merged = pd.merge(Mathys2023_Meta, ROSMAP_Status, left_on="individualID", right_on="individualID", how='left')
Mathys2023_Merged = Mathys2023_Merged.set_index("index_x")
if list(Mathys2023_Merged.index) == list(Mathys2023.obs_names):
    Mathys2023.obs = Mathys2023_Merged
Mathys2023.obs.rename(
    columns={
        'Status_y': 'Diagnosis'
    },
    inplace=True
)
Mathys2023.obs['Diagnosis'].value_counts()
Mathys2023.obs["Disease"] = Mathys2023.obs["Diagnosis"].map({
    'Dementia_Two Pathology': 'AD',
    'Normal': 'NC',
    'Non Dementia_Neurofibrillary Tangle': 'NC',
    'Non Dementia_Two Pathology': 'NC',
    'Dementia_No Pathology': 'AD',
    'Dementia_Neurofibrillary Tangle': 'AD',
    'Non Dementia_Neuritic Plaques': 'NC'
})
Mathys2023.obs['Disease'].value_counts()
Mathys2023.obs['APOEStatus'].value_counts()
Mathys2023.obs["APOEStatus"] = Mathys2023.obs["APOEStatus"].map(mapping)
Mathys2023.obs['APOEStatus'].value_counts()
Mathys2023
del Mathys2023.layers['shifted_log']
Mathys2023.write_h5ad("~/Merged/Finalised/Mathys2023_PFCMicrogliaFinal.h5ad")

#### Sun et al, 2023 ####
Sun2023 = sc.read_h5ad("~/Merged/H5AD/Sun2023_PFCMicroglia.h5ad")
Sun2023
ROSMAP_Full = pd.read_csv("~/MetaData/ROSMAP_Metadata_All.csv")
ROSMAP_Clinical = pd.read_table("~/MetaData/clinical.status.txt")
ROSMAP_Full_dedup = ROSMAP_Full[~ROSMAP_Full.index.duplicated(keep='first')]
Sun2023_Meta = Sun2023.obs.copy()
Sun2023_Merged = pd.merge(Sun2023_Meta, ROSMAP_Full_dedup, on="subject", how='left')

Sun2023_FullMeta = Sun2023_Merged.merge(ROSMAP_Clinical, left_on="individualID", right_on="individualID", how="left")
Sun2023_Meta = Sun2023.obs.copy()
Sun2023_FullMeta = Sun2023_FullMeta.set_index("barcode_x")
if list(Sun2023_FullMeta.index) == list(Sun2023.obs_names):
    Sun2023.obs = Sun2023_FullMeta
Sun2023.write_h5ad("~/H5AD/Sun2023_FullFinal.h5ad")

Sun2023 = Sun2023[Sun2023.obs['brainRegion'] == "PFC"]
Sun2023.write_h5ad("/home/jacquelyne/Documents/Aging/202512_IronLoadedMicroglia/Merged/H5AD/Sun2023_PFCMicroglia.h5ad")
Sun2023

columns_to_drop = ['orig.ident','age_death_x', 'msex_x','pmi_x','msex_y']
Sun2023.obs.drop(columns=columns_to_drop, inplace=True)
Sun2023
Sun2023.obs.rename(
    columns={
        "subject": "Sample",
        "Gender": "Sex",
        "pmi_y": "PMI",
        "braaksc": "BraakStage",
        "ceradsc": "CERADScore",
    },
    inplace=True
)
Sun2023.obs["Sex"].value_counts()
Sun2023.obs["ADdiag3types"].value_counts()
Sun2023.obs['Disease'] = Sun2023.obs['ADdiag3types'].map({
    'nonAD': 'NC',
    'earlyAD':"AD",
    'lateAD':'AD'
})
Sun2023.obs['Disease'].value_counts()
Sun2023.obs['apoe_genotype'].value_counts()
mapping = {
    33.0: "E3/E3",
    34.0: "E3/E4",
    23.0: "E2/E3",
    24.0: "E2/E4",
    44.0: "E4/E4",
    22.0: "E2/E2",
}

Sun2023.obs["APOEStatus"] = Sun2023.obs["apoe_genotype"].map(mapping)
Sun2023.obs["APOEStatus"].value_counts()
Sun2023
columns_to_remove = ['leiden0_1', 'leiden0_2', 'leiden0_3', 'leiden0_4', 'leiden0_5', 'leiden0_6', 'leiden0_7', 'leiden0_8', 'leiden0_9', 'leiden1_0', 'leiden1_1', 'leiden1_2', 'leiden1_3', 'leiden1_4']
Sun2023.obs = Sun2023.obs.drop(columns = columns_to_remove)
del Sun2023.obsm
del Sun2023.obsp
del Sun2023.layers['shifted_log']
del Sun2023.varm
columns_to_remove = ['n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_deviant', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm']
Sun2023.var = Sun2023.var.drop(columns=columns_to_remove)
del Sun2023.uns
Sun2023.var_names 
Sun2023.obs['Sex'].value_counts()
Sun2023.write_h5ad("~/Merged/Finalised/Sun2023_PFCMicrogliaFinal.h5ad")

#### Liu et al, 2025 ####
Liu2025 = sc.read_h5ad("~/Merged/H5AD/Liu2025_PFCMicroglia.h5ad")
Liu2025.obs['Dataset'] = "Liu2025"
Liu2025
ROSMAP_Clinical = pd.read_table("~/MetaData/clinical.status.txt")
Liu2025_Meta = Liu2025.obs.copy()
Liu2025_Meta = Liu2025_Meta.reset_index()
Liu2025_Merged = pd.merge(Liu2025_Meta, ROSMAP_Clinical, left_on="ROSMAP_IndividualID", right_on="individualID", how='left')
Liu2025_Merged = Liu2025_Merged.set_index("index")
if list(Liu2025_Merged.index) == list(Liu2025.obs_names):
    Liu2025.obs = Liu2025_Merged
columns_to_drop = ['age_death_x','pmi_x','Study_y']
Liu2025.obs.drop(columns=columns_to_drop, inplace=True)
Liu2025
Liu2025.obs.rename(
    columns={
        "pmi_y": "PMI",
        "braaksc": "BraakStage",
        "ceradsc": "CERADScore",
        "apoe_genotype": "APOEStatus"
    },
    inplace=True
)
Liu2025.obs["Sex"].value_counts()
Liu2025.obs["Pathology"].value_counts()
Liu2025.obs['Disease'] = Liu2025.obs['Pathology'].map({
    'nonAD': 'NC',
    'earlyAD':"AD",
    'lateAD':'AD'
})
Liu2025.obs['Disease'].value_counts()
Liu2025.obs['APOEStatus'].value_counts()
mapping = {
    33.0: "E3/E3",
    34.0: "E3/E4",
    23.0: "E2/E3",
    24.0: "E2/E4",
    44.0: "E4/E4",
    22.0: "E2/E2",
}

Liu2025.obs["APOEStatus"] = Liu2025.obs["APOEStatus"].map(mapping)
Liu2025.obs["APOEStatus"].value_counts()
Liu2025
del Liu2025.obsm
del Liu2025.uns
del Liu2025.varm
del Liu2025.obsp
del Liu2025.layers['shifted_log']
Liu2025.write_h5ad(~/Merged/Finalised/Liu2025_PFCMicrogliaFinal.h5ad")

#### Dharshini et al, 2025 ####
Dharshini2025 = sc.read_h5ad("~/Merged/H5AD/Dharshini2025_PFCMicroglia.h5ad")
Dharshini2025.obs['Dataset'] = "Dharshini2025"
Dharshini2025
Dharshini2025_meta = Dharshini2025.obs
Dharshini2025_meta.to_csv("~/Merged/H5AD/Dharshini2025_meta.csv", index=True)
Dharshini2025
columns_to_drop = ['Brain.Region', 'SORT','tissue_ontology_term_id', 'assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'disease_ontology_term_id', 'sex_ontology_term_id', 'is_primary_data', 'organism_ontology_term_id', 'suspension_type','tissue_type', 'cell_type', 'assay','organism','tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid',]
Dharshini2025.obs.drop(columns=columns_to_drop, inplace=True)
Dharshini2025
Dharshini2025.obs.rename(
    columns={
        "APOE_class": "APOE_Status",
        "Braak.stage": "BraakStage",
        "PMI.hr.": "PMI",
        "disease": "Disease",
        "sex": "Sex",
        "sample": "Sample"},
        inplace=True
)
Dharshini2025.obs['Sex'].value_counts()
Dharshini2025.obs['Sex'] = Dharshini2025.obs['Sex'].map({"female":"Female", "male":"Male"})
Dharshini2025.obs['Disease'].value_counts()
Dharshini2025.obs['Disease'] = Dharshini2025.obs['Disease'].map({"Alzheimer disease":"AD", "normal":"NC"})
Dharshini2025.obs['APOEStatus'].value_counts()
Dharshini2025
columns_to_remove = ['n_genes_by_counts', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'outlier', 'mt_outlier']
Dharshini2025.obs = Dharshini2025.obs.drop(columns=columns_to_remove)
del Dharshini2025.uns
del Dharshini2025.obsm
del Dharshini2025.layers['shifted_log']
del Dharshini2025.layers['spliced']
del Dharshini2025.layers['unspliced']
Dharshini2025.write_h5ad("~/Merged/Finalised/Dharshini2025_PFCMicrogliaFinal.h5ad")
