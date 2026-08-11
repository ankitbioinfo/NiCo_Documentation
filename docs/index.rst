.. NiCo documentation master file, created by
   sphinx-quickstart on Mon Dec 11 23:14:58 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NiCo: Niche covariation analysis of spatial transcriptomics data
=====================================================================

.. important::

   **🚀 NEW USERS — START HERE!**

   For **new NiCo analyses**, please use `nico-wrapper`_ as the
   **preferred way to apply NiCo**.

   Since August 1, 2026, `nico-wrapper`_ provides streamlined **command-line and Python interfaces**
   for the main NiCo workflow while continuing to use **NiCo for the
   underlying scientific computations**.

   **👉 Go to `nico-wrapper`_ on Gruenlab GitHub**

   The direct NiCo installation and tutorials on this website are primarily
   retained for **existing users, legacy analyses, and reproducibility**.


.. figure:: /_static/Figure1.png
    :align: center
    :alt: NiCo spatial analysis

    Infer cellular crosstalk from spatial transcriptomics and scRNA-seq data


.. important::

   **🆕 Recommended workflow for new analyses**

   If you are starting a **new project with NiCo**, please start with
   `nico-wrapper`_ rather than the legacy direct NiCo installation.

   `nico-wrapper`_ provides a streamlined interface to the main NiCo workflow
   through both **CLI and Python interfaces**, while NiCo remains the
   underlying scientific computation framework.

   **👉 Start with `nico-wrapper`_**

   Existing users who need to reproduce previous analyses can continue using
   the direct NiCo package and the tutorials provided in this documentation.


The Niche Covariation (NiCo) package is developed for the integration of single-cell resolution
spatial transcriptomics and scRNA-seq data (or from sequencing-based spatial transcriptomics data alone)
to (1) perform cell type annotations in the spatial modality by label transfer, (2) predict niche cell type
interactions within local neighborhoods, and (3) infer cell state covariation and the underlying molecular
crosstalk in the niche. NiCo infers factors capturing cell state variability in both modalities and
identifies genes correlated to these latent factors for the prediction of ligand-receptor interactions
and factor-associated pathways.

For **new analyses**, we recommend using `nico-wrapper`_ as the preferred interface
to the NiCo workflow.

Highlights of NiCo
==================
1. Annotations of cell types in spatial data by label transfer
2. Prediction of niche interactions using neighborhood analysis
3. Covariation analysis of latent factors across niche cell types
4. Prediction of ligand-receptor interactions mediating niche crosstalk
5. Inference of pathways aassociated with covarying cell states


Installation
=======================

.. warning::

   **For new analyses, please use `nico-wrapper`_ instead of the direct NiCo
   installation below.**

   Since August 1, 2026, `nico-wrapper`_ provides a streamlined command-line
   and Python interface for the main NiCo workflow while continuing to use
   NiCo as the underlying scientific computation framework.

.. note:: Please install using following commands:

      | conda create -n nicoUser python=3.11
      | conda activate nicoUser
      | pip install nico-sc-sp
      | pip install jupyterlab


For more details, follow the python package index guidelines from `nico-sc-sp pypi <https://pypi.org/project/nico-sc-sp/>`_

Tutorials
=====================

.. note::

   **New users:** For new analyses, please start with `nico-wrapper`_.

   The tutorials below describe the direct NiCo workflow and are retained for
   existing users, previous analyses, and reproducibility.

Please prepare the input files with scRNA-seq count data and cell type annotation
(cluster partition), spatial count data, and spatial cell coordinates to run the
complete NiCo tutorials.

`NiCo tutorials for imaging-based spatial transcriptomics (Xenium, MERSCOPE, seqFISH)
or sequencing-based methods, e.g., Slide-seqV2, are available here <nico-tutorial-github>`_



.. _nico-wrapper: https://github.com/gruenlab/nico-wrapper
.. _nico-wrapper-github: https://github.com/gruenlab/nico-wrapper
.. _nico-tutorial-github: https://github.com/ankitbioinfo/nico_tutorial


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
