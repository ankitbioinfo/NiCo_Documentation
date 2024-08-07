Installation
============

This section provides instructions on how to install Nico SC-SP on your system.

Prerequisites
-------------

Before installing NiCo, ensure you have the following dependencies installed:

- Python 3.11 or later
- Anaconda or Miniconda

Installing via Pip
--------------------

The recommended way to install NiCo is via Pip. Follow these steps:

1. Create a new Conda environment:
`conda create -n nicoUser python=3.11`

2. Activate the environment:
`conda activate nicoUser`

3. Install the required dependencies and NiCo:
`conda install -c conda-forge pygraphviz`
`pip install nico-sc-sp`
`pip install jupyterlab`


Verifying the Installation
--------------------------

To verify that NiCo is installed correctly, run the following command:
>>>import nico
>>>


If the installation is successful, you should see a blank >>> prompt.

Updating Nico SC-SP
-------------------

To update NiCo to the latest version, navigate to the cloned directory and pull the latest changes:



Troubleshooting
---------------

If you encounter any issues during installation, please refer to the [FAQ](faq.html) section or contact support at ankitplusplus at gmail.com.
