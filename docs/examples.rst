Examples
========

This section provides examples of how to use the Nico SC-SP library for spatial cell type annotation, niche interaction inference, and covariation analysis.


Tutorials
=====================
Please prepare the input files with scRNA-seq count data and cell type annotation (cluster partition), spatial count data, and spatial
cell coordinates to run the complete NiCo tutorials.

`Nico tutorials for Xenium, MERSCOPE, SEQFISH and SlideSeqV2 spatial technologies are available here <https://github.com/ankitbioinfo/nico_tutorial>`_


Basic Usage
-----------

Below is a complete example demonstrating the usage of the main steps (modules) in NiCo.

### Importing Modules and Setting Up

```console
import Annotations as sann
import Interactions as sint
import Covariations as scov
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

# Configure matplotlib for publication-quality plots
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Tahoma', 'DejaVu Sans', 'Lucida Grande', 'Verdana']
plt.rcParams['pdf.fonttype'] = 42  # Embed fonts in PDF files

import warnings
warnings.filterwarnings("ignore")

# Parameters for saving plots
saveas = 'png'
transparent_mode = False

ref_datapath = './inputRef/'
query_datapath = './inputQuery/'
output_nico_dir = './nico_analysis/'
output_annotation_dir = None
annotation_save_fname = 'nico_celltype_annotation.h5ad'
inputRadius = 0
ref_cluster_tag = 'cluster'  # scRNAseq cell type slot
annotation_slot = 'nico_ct'  # spatial cell type slot
```

Module A: Perform Cell Type Annotation of Spatial Data


# Find anchor cells between reference and query datasets

anchors_and_neighbors_info = sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir
)

# Perform NiCo-based annotation
output_info = sann.nico_based_annotation(
anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.4',
across_spatial_clusters_dispersion_cutoff=0.15,
resolved_tie_issue_with_weighted_nearest_neighbor='No'
)

# Clean up temporary files
sann.delete_files(output_info)

# Save annotations to the spatial object
sann.save_annotations_in_spatial_object(output_info, anndata_object_name=annotation_save_fname)



Module A: Cell Type Annotation Visualization

print('\n\nModule A visualization')
# Visualize UMAP and cell coordinates with all cell types
sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag=annotation_slot,
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,
transparent_mode=transparent_mode
)

# Visualize UMAP and cell coordinates with selected cell types
sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag=annotation_slot,
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
choose_celltypes=[],
saveas=saveas,
transparent_mode=transparent_mode
)



Module B: Infer Significant Niche Cell Type Interactions
print('\n\nModule B')
do_not_use_following_CT_in_niche = ['Basophils', 'Cycling/GC B cell', 'pDC']

niche_pred_output = sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag=annotation_slot,
removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff = 0.1

sint.plot_niche_interactions_with_edge_weight(
niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode
)

sint.plot_niche_interactions_without_edge_weight(
niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode
)

sint.find_interacting_cell_types(
niche_pred_output,
choose_celltypes=[],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
coeff_cutoff=30,
saveas=saveas,
transparent_mode=transparent_mode,
figsize=(4.0, 2.0)
)

sint.plot_confusion_matrix(
niche_pred_output,
saveas=saveas,
transparent_mode=transparent_mode
)

sint.plot_coefficient_matrix(
niche_pred_output,
saveas=saveas,
transparent_mode=transparent_mode
)

sint.plot_evaluation_scores(
niche_pred_output,
saveas=saveas,
transparent_mode=transparent_mode,
figsize=(4, 3)
)



Module C: Perform Niche Cell State Covariation Analysis Using Latent Factors
print('\n\nModule C')
cov_out = scov.gene_covariation_analysis(
iNMFmode=True,
Radius=inputRadius,
no_of_factors=3,
spatial_integration_modality='double',
refpath=ref_datapath,
quepath=query_datapath,
output_niche_prediction_dir=output_nico_dir,
ref_cluster_tag=ref_cluster_tag
)

# Visualize the correlation of genes from NMF
scov.plot_cosine_and_spearman_correlation_to_factors(
cov_out,
choose_celltypes=[],
NOG_Fa=30,
saveas=saveas,
transparent_mode=transparent_mode,
figsize=(15, 10)
)

scov.make_excel_sheet_for_gene_correlation(cov_out)




Module D: Cell Type Covariation Visualization
print('\n\nModule D')
scov.plot_significant_regression_covariations_as_circleplot(
cov_out,
choose_celltypes=[],
pvalue_cutoff=0.05,
mention_pvalue=True,
saveas=saveas,
transparent_mode=transparent_mode,
figsize=(6, 1.25)
)




Module E: Analysis of Ligand-Receptor Interactions Within the Cell Type Covariation State
print('\n\nModule E')
scov.save_LR_interactions_in_excelsheet_and_regression_summary_in_textfile_for_interacting_cell_types(
cov_out,
pvalueCutoff=0.05,
correlation_with_spearman=True,
LR_plot_NMF_Fa_thres=0.1,
LR_plot_Exp_thres=0.1,
number_of_top_genes_to_print=5
)

scov.find_LR_interactions_in_interacting_cell_types(
cov_out,
choose_interacting_celltype_pair=[],
choose_factors_id=[],
pvalueCutoff=0.05,
LR_plot_NMF_Fa_thres=0.2,
LR_plot_Exp_thres=0.2,
saveas=saveas,
transparent_mode=transparent_mode,
figsize=(12, 10)
)




Module F: Perform Functional Enrichment Analysis for Genes Associated with Latent Factors
print('\n\nModule F')
scov.pathway_analysis(
cov_out,
choose_celltypes=[],
NOG_pathway=50,
choose_factors_id=[],
savefigure=True,
positively_correlated=True,
saveas='pdf',
rps_rpl_mt_genes_included=False
)




Module G: Visualization of Top Genes Across Cell Type and Factors as Dotplot
print('\n\nModule G')
scov.plot_top_genes_for_a_given_celltype_from_all_three_factors(
cov_out,
choose_celltypes=[],
top_NOG=20,
saveas=saveas,
transparent_mode=transparent_mode
)

scov.plot_top_genes_for_pair_of_celltypes_from_two_chosen_factors(
cov_out,
choose_interacting_celltype_pair=['Stem/TA', 'Paneth'],
visualize_factors_id=[1, 1],
top_NOG=20,
saveas=saveas,
transparent_mode=transparent_mode
)



Module H: Visualize Factor Values in the UMAP
print('\n\nModule H')

scov.visualize_factors_in_spatial_umap(
cov_out,
visualize_factors_id=[1, 1],
choose_interacting_celltype_pair=['Stem/TA', 'Paneth'],
saveas=saveas,
transparent_mode=transparent_mode,
figsize=(8, 3.5)
)

scov.visualize_factors_in_scRNAseq_umap(
cov_out,
choose_interacting_celltype_pair=['Stem/TA', 'Paneth'],
visualize_factors_id=[1, 1],
saveas=saveas,
transparent_mode=transparent_mode,
figsize=(8, 3.5)
)
