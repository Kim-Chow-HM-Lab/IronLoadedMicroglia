import os
import seaborn as sns
import palantir
import scFates as scf

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

ADMicroglia = sc.read_h5ad("~/H5AD/ADMicroglia_ProteinCodingGenes.h5ad")

# load, normalize and log-transform count data
ADMicroglia.X = ADMicroglia.layers['raw_counts'].copy()
counts = ADMicroglia.to_df()
norm_df = sc.pp.normalize_per_cell(counts, copy=True)
norm_df = palantir.preprocess.log_transform(norm_df)

# HVG selection
ADMic_scf = sc.AnnData(norm_df)
sc.pp.highly_variable_genes(ADMic_scf, n_top_genes = 1500, flavor = 'cell_ranger')
ADMic_scf.obs = ADMicroglia.obs.copy()
sc.pp.pca(ADMic_scf)
projections = pd.DataFrame(ADMic_scf.obsm['X_pca'], index=ADMic_scf.obs_names)

# running Palantir to obtain multiscale diffusion maps
dm_res = palantir.utils.run_diffusion_maps(projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs=8)

# generate embedding from the multiscale diffusion space
ADMic_scf.obsm['X_palantir'] = ms_data.values
sc.pp.neighbors(ADMic_scf, n_neighbors=30, use_rep = "X_palantir")
ADMic_scf.obsm['X_palantirpca2d'] = ADMic_scf.obsm['X_pca'][:,:2]
sc.tl.draw_graph(ADMic_scf, init_pos='X_palantirpca2d')
sc.pl.draw_graph(ADMic_scf, color='P2RY12', color_map='RdBu_r', save='2_palantir_p2ry12expression', size=20)
sc.pl.draw_graph(ADMic_scf, color="leiden0_7", save="3_palantir_leiden07")

# tree learning with elastic principal graph
scf.tl.tree(ADMic_scf, method = "epg", Nodes =50, use_rep="palantir", device="cpu", seed=123)
scf.pl.graph(ADMic_scf, save = "4_palantir_trajectorynodes")
scf.tl.root(ADMic_scf, 6)
scf.tl.pseudotime(ADMic_scf, n_jobs=20, n_map=100, seed=123)
scf.pl.trajectory(ADMic_scf, save = "4_palantir_trajectorypseudotime", show_info = True)
sc.pl.draw_graph(ADMic_scf, color="leiden0_7", size=5, legend_fontoutline=True, show_info = True, save="5_palantir_clustertrajectories")
sc.pl.draw_graph(ADMic_scf, color = "t", size=5, save = "10_pseudotime", cmap = "viridis")
scf.tl.dendrogram(ADMic_scf)
with plt.rc_context({"figure.figsize": (4,8)}):
    scf.pl.dendrogram(ADMic_scf, linewidth_seg = 1, tree_behind=True, color="leiden0_7", size=10, save="6_palantir_dendrogram")

sc.set_figure_params(figsize=(3,6),frameon=False,dpi_save=300)
scf.pl.dendrogram(ADMic_scf, color="t", show_info=False, linewidth_seg=1, tree_behind=True, save="7_time.pdf", cmap="viridis")
scf.pl.dendrogram(ADMic_scf, linewidth_seg=1, tree_behind=True, color="leiden0_7", size=20, save="8_palantir_dendrogram")
scf.pl.trends(ADMic_scf, root_milestone = "6")
scf.pl.trajectory(ADMic_scf, basis="draw_graph_fa", color_cells="leiden0_7", size=20, show_info = True, scale_path = 0.5, rasterized = False, save="9_palantir_trajectory_umap")

