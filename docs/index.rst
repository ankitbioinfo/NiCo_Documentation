.. NiCo documentation master file, created by
   sphinx-quickstart on Mon Dec 11 23:14:58 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NiCo: Niche covariation analysis of spatial transcriptomics data
================================================================



.. figure:: /_static/Figure1.png
    :align: center
    :alt: NiCo spatial analysis

    Infer cellular crosstalk from spatial transcriptomics and scRNAseq data



The Niche Covariation (NiCo) package is developed for the integration of single-cell resolution
spatial transcriptomics and scRNA-seq data (or from sequencing-based spatial transcriptomics data alone)
to (1) perform cell type annotations in the spatial modality by label transfer, (2) predict niche cell type 
interactions within local neighborhoods, and (3) infer cell state covariation and the underlying molecular 
crosstalk in the niche. NiCo infers factors capturing cell state variability in both modalities and 
identifies genes correlated to these latent factors for the prediction of ligand-receptor interactions 
and factor-associated pathways.

Highlights of NiCo
==================
1. Annotations of cell types in spatial data by label transfer
2. Prediction of niche interactions using neighborhood analysis
3. Covariation analysis of latent factors across niche cell types
4. Prediction of ligand-receptor interactions mediating niche crosstalk
5. Inference of pathways aassociated with covarying cell states


Installation
=======================


.. note:: Please install using following commands:

      | conda create -n nicoUser python=3.11
      | conda activate nicoUser
      | pip install nico-sc-sp
      | pip install jupyterlab


For more details, follow the python package index guidelines from `nico-sc-sp pypi <https://pypi.org/project/nico-sc-sp/>`_

Tutorials
=====================
Please prepare the input files with scRNA-seq count data and cell type annotation (cluster partition), spatial count data, and spatial
cell coordinates to run the complete NiCo tutorials.

`Nico tutorials for imaging-based spatial transcriptomics (Xenium, MERSCOPE, seqFISH) or sequencing-based methods, e.g., Slide-seqV2, are available here <https://github.com/ankitbioinfo/nico_tutorial>`_



.. toctree::
   :maxdepth: 3
   :caption: Contents:


   introduction
   installation
   examples
   tutorial0
   tutorial1
   tutorial2
   api_reference
   faq


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
