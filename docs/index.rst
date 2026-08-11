```rst
.. NiCo documentation master file, created by
   sphinx-quickstart on Mon Dec 11 23:14:58 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NiCo: Niche covariation analysis of spatial transcriptomics data
================================================================

.. important::

   **🚀 NEW USERS — START HERE!**

   For **new NiCo analyses**, please use **`nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_**
   as the **preferred way to apply NiCo**.

   `nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_ provides streamlined
   **command-line and Python interfaces** for the main NiCo workflow while continuing
   to use **NiCo for the underlying scientific computations**.

   **👉 `Go to nico-wrapper documentation and GitHub repository <https://github.com/gruenlab/nico-wrapper>`_**

   The direct NiCo installation and tutorials documented on this website are primarily
   retained for **existing users, legacy analyses, and reproducibility**.

.. figure:: /_static/Figure1.png
    :align: center
    :alt: NiCo spatial analysis

    Infer cellular crosstalk from spatial transcriptomics and scRNAseq data


.. important::

   **🆕 Recommended workflow for new analyses**

   If you are starting a **new project with NiCo**, please do **not** start with the
   legacy installation instructions below.

   **Use `nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_ instead.**

   `nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_ provides a streamlined
   interface to the NiCo workflow and supports both **CLI and Python-based analysis**,
   while NiCo remains the underlying scientific computation framework.

   **👉 `Start with nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_**

   Existing users who need to reproduce previous analyses can continue using the
   NiCo package and the tutorials provided in this documentation.

The Niche Covariation (NiCo) package is developed for the integration of single-cell resolution
spatial transcriptomics and scRNA-seq data (or from sequencing-based spatial transcriptomics data alone)
to (1) perform cell type annotations in the spatial modality by label transfer, (2) predict niche cell type
interactions within local neighborhoods, and (3) infer cell state covariation and the underlying molecular
crosstalk in the niche. NiCo infers factors capturing cell state variability in both modalities and
identifies genes correlated to these latent factors for the prediction of ligand-receptor interactions
and factor-associated pathways.

For **new analyses**, we recommend using
`nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_ as the preferred interface
to the NiCo workflow. The wrapper simplifies the execution of the main NiCo workflow
while retaining NiCo's underlying scientific methods and computational functionality.


Highlights of NiCo
==================

1. Annotation of cell types in spatial data by label transfer
2. Prediction of niche interactions using neighborhood analysis
3. Covariation analysis of latent factors across niche cell types
4. Prediction of ligand-receptor interactions mediating niche crosstalk
5. Inference of pathways associated with covarying cell states


Recommended workflow for new users
==================================

.. note::

   **Starting a new NiCo analysis?**

   We strongly recommend starting with
   `nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_.

   It provides:

   * A streamlined command-line interface
   * A Python interface for programmatic workflows
   * Simplified execution of the main NiCo workflow
   * Access to the established NiCo scientific computations

   **`Visit nico-wrapper on GitHub → <https://github.com/gruenlab/nico-wrapper>`_**


Installation
============

.. warning::

   The direct NiCo installation below is provided primarily for
   **existing users, legacy workflows, and reproducibility**.

   For **new analyses**, please use
   `nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_.

.. note:: Please install using the following commands:

      | conda create -n nicoUser python=3.11
      | conda activate nicoUser
      | pip install nico-sc-sp
      | pip install jupyterlab


For more details, follow the Python Package Index guidelines from
`nico-sc-sp PyPI <https://pypi.org/project/nico-sc-sp/>`_.


Tutorials
=========

.. note::

   **New users:** For new analyses, please start with
   `nico-wrapper <https://github.com/gruenlab/nico-wrapper>`_.

   The tutorials below describe the direct NiCo workflow and are retained for
   existing users, previous analyses, and reproducibility.

Please prepare the input files with scRNA-seq count data and cell type annotation
(cluster partition), spatial count data, and spatial cell coordinates to run the
complete NiCo tutorials.

`NiCo tutorials for imaging-based spatial transcriptomics
(Xenium, MERSCOPE, seqFISH) or sequencing-based methods, e.g., Slide-seqV2,
are available here
<https://github.com/ankitbioinfo/nico_tutorial>`_


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
```
