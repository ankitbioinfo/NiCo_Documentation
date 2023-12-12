.. NiCo documentation master file, created by
   sphinx-quickstart on Mon Dec 11 23:14:58 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NiCo: Niche analysis of single-cell resolution of spatial cells
================================================================



.. figure:: /_static/Figure1.png
    :align: center
    :alt: NiCo spatial analysis

    Infer cellular crosstalk from imaging-based spatial transcriptomics and scRNAseq data



The Niche Covariation (NiCo) package is developed for the integration of single-cell resolution of
spatial transcriptomics and scRNA-seq data to perform the annotations in the spatial cells using the
label transfer, predict the niche cell type interactions using the neighborhood, and find the covariation
relationship in the niche using the latent factors of both modalities. The genes in the latent factors
are used to find the ligand-receptor and pathway analysis.

Highlights of NiCo
==================
1. Annotations of spatial cells using label transfer
2. Prediction of niche using neighborhood analysis
3. Covariation analysis of the interacted niche in the latent factors
4. Ligand-receptor analysis in the latent factors
5. Pathway analysis in the latent factors


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
Please prepare the scRNAseq count data, cell type annotation cluster, spatial count data, and spatial
cell coordinates files to run the complete NiCo tutorials.
`Nico tutorials are available here <https://github.com/ankitbioinfo/nico_tutorial>`_



.. toctree::
   :maxdepth: 3
   :caption: Contents:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
