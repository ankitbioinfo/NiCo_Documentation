Installation
============

This section provides instructions on how to install NiCo on your system.

Prerequisites
-------------

Before installing NiCo, ensure you have the following dependencies installed:

- Python 3.11 or later
- Anaconda or Miniconda

Installing via Pip
------------------

The recommended way to install NiCo is via Pip. Follow these steps:

1. Create a new Conda environment:

.. code-block:: console

   conda create -n nicoUser python=3.11


2. Activate the environment:

.. code-block:: console

   conda activate nicoUser


3. Install the required dependencies and NiCo:

.. code-block:: console

   conda install -c conda-forge pygraphviz
   pip install nico-sc-sp
   pip install jupyterlab



Verifying the Installation
--------------------------

To verify that NiCo is installed correctly, run the following command:

.. code-block:: console

   >>>import nico
   >>>


If the installation is successful, you should see a blank ```>>>``` prompt.


If you do not wish to install NiCo, you can use the scripts directly,
`available here, <https://github.com/ankitbioinfo/nico_tutorial/tree/main/NiCo>`_ and call them at the beginning of the tutorials.

.. code:: ipython3

    import Annotations as sann
    import Interactions as sint
    import Covariations as scov


Updating NiCo
-------------------

To update NiCo to the latest version, navigate to the cloned directory and pull the latest changes.



Troubleshooting
---------------

If you encounter any issues during installation, please refer to the FAQ section or contact support at ankitplusplus at gmail.com.
