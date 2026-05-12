Installation
============

VGTK uses a **conda environment** to install the dependency packages.  
The ``environment.yml`` file contains all the packages and tools necessary to run VGTK.

Download Conda
--------------

You can install Miniconda from the official website:

`Miniconda Installation Guide <https://www.anaconda.com/docs/getting-started/miniconda/install>`_

Add Required Channels
---------------------

After installation, add the following channels and set strict priority:

.. code-block:: bash

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict


Installing VGTK Packages
________________________

To set up VGTK, create the conda environment using the provided ``environment.yml`` file:

.. code-block:: bash

   conda env create --file environment.yml

This command will create a conda environment named **VGTK** with all required dependencies.  
Before running the VGTK tool, activate the environment:

.. code-block:: bash

   conda activate VGTK


