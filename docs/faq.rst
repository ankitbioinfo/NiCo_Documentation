Frequently Asked Questions (FAQ)
=================================

This section provides answers to some common questions regarding the NiCo library.

Q. What is NiCo?

A. NiCo is a Python library designed for spatial transcriptomics analysis. It provides modules for cell type annotation, niche interaction inference, and covariation analysis within the niche cell types in spatial transcriptomics data.

Q. How do I install NiCo?

A. You can install NiCo using pip or just keep the python scripts available at github in your analysis path. For detailed installation instructions, refer to the :doc:`installation` section.

Q. What types of data does NiCo support?

A. NiCo supports imaging- or sequencing-based spatial transcriptomics data. It requires an expression matrix in scTransform-like normalization, typically stored in an AnnData object.

Q. How can I annotate cell types in my spatial data?

A. You can use the `Annotations` module provided by NiCo. The `find_anchor_cells_between_ref_and_query` and `nico_based_annotation` functions help perform cell type annotation. Refer to the :doc:`examples` section for a detailed example.

Q. How do I visualize cell type annotations?

A. NiCo offers several visualization functions in the `Annotations` module, such as `visualize_umap_and_cell_coordinates_with_all_celltypes` and `visualize_umap_and_cell_coordinates_with_selected_celltypes`. These functions help visualize expression UMAP and cell coordinates for the annotated cell types.

Q. Can NiCo infer niche interactions?

A. Yes, the `Interactions` module of NiCo provides functions for inferring significant niche cell type interactions, such as `spatial_neighborhood_analysis` and `plot_niche_interactions_with_edge_weight`.

Q. What is the purpose of covariation analysis?

A. Covariation analysis identifies significance cell state covariations between a central cell type and its niche cell types using latent factors. The `Covariations` module in NiCo offers tools like `gene_covariation_analysis` and `plot_cosine_and_spearman_correlation_to_factors` for this purpose.

Q. How can I perform ligand-receptor interaction analysis?

A. The `Covariations` module provides functions like `save_LR_interactions_in_excelsheet_and_regression_summary_in_textfile_for_interacting_cell_types` and `find_LR_interactions_in_interacting_cell_types` for analyzing ligand-receptor interactions between covarying cell types.

Q. Where can I find detailed usage examples and tutorials?

A. Refer to the :doc:`examples` section for detailed usage examples covering the main functionalities of NiCo and github links for the tutorials.

Q. Who can I contact for support?

A. For support, you can open an issue on the project's GitHub repository or contact the maintainers directly through the provided contact information on the repository page.

Q. Is there a way to contribute to NiCo?

A. Yes, contributions are welcome! You can contribute by reporting issues, suggesting new features, or submitting pull requests on the project's GitHub repository.
