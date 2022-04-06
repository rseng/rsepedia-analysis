CanDI - A global cancer data integrator
=======================================

|Documentation Status|

Package Installation
--------------------

First, you need to clone this repository to use CanDI.

.. code:: bash

   git clone https://github.com/GilbertLabUCSF/CanDI.git

We suggest to use `Conda <https://docs.conda.io/en/latest/>`__ as a
package manager and environment management system. You can create a
fresh conda environment with all ``CanDI``\ ’s requirements using bellow
command:

.. code:: bash

   conda env create -f CanDI/candi.yml -n candi

Prepare Datasets
~~~~~~~~~~~~~~~~

The python command from CanDI will automatically download and modify
datasets.

.. code:: bash

   python CanDI/CanDI/setup/install.py

Downloaded and formatted datasets would organize this way:

.. code::

   .
   ├── config.ini # modified after Installation 
   ├── depmap
   │   ├── CCLE_expression.csv
   │   ├── CCLE_fusions.csv
   │   ├── CCLE_gene_cn.csv
   │   ├── CCLE_mutations.csv
   │   ├── CCLE_RNAseq_reads.csv
   │   ├── CRISPR_gene_dependency.csv
   │   ├── CRISPR_gene_effect.csv
   │   └── sample_info.csv
   ├── genes
   │   └── gene_info.csv
   └── locations
       └── merged_locations.csv

Package Usage
-------------

Import CanDI into python
~~~~~~~~~~~~~~~~~~~~~~~~

To import ``CanDI``, your active directory in python must be same as the
cloned folder.

.. code:: python

   from CanDI import candi

**OR**, you can add path to the `CanDI` directory if you want to use it from other directories.

.. code:: python

   import sys
   sys.path.append("path-to-candi-directory")

   from CanDI import candi

CanDI Objects
~~~~~~~~~~~~~

-  ``data`` : Container for all candi datasets. All access to datasets
   go through data object.
-  ``Gene`` : Provides cross dataset indexing from the gene perspective.
-  ``CellLine`` : Provides cross dataset indexing from the cell line
   perspective.
-  ``Cancer`` : Provides cross dataset indexing by a group of cell lines
   that are all the same tissue.
-  ``Organelle``: Provides cross dataset indexing for a group of genes
   whose proteins localize to the same organelle.
-  ``CellLineCluster`` : Provides cross dataset indexing for a group of
   user defined cell lines.
-  ``GeneCluster`` : Provides cross dataset indexing for a group of user
   defined genes.

.. |Documentation Status| image:: https://readthedocs.org/projects/candi/badge/?version=latest
   :target: https://candi.readthedocs.io/en/latest/?badge=latest
Module Structure:
=======================
The CanDI data integrator is a python library built on top of the Pandas that is specialized in integrating the publicly available data from:

- The Cancer Dependency Map
- The Cancer Cell Line Encyclopedia
- The Comprehensive Resource of Mammalian Protein Complexes (CORUM)
- Protein Localization Data from: The Cell Atlas, Map of the Cell, and The In Silico Surfaceome.

Access to all datasets is controlled via a python class called Data. Upon import the data class reads the config file established during installation and defines unique paths to each dataset and automatically loads the cell line index table and the gene index table.
Installation of CanDI, configuration, and data retrieval is handled by a manager class that is accessed indirectly through installation scripts and the Data class.
Interactions with this data are controlled through a parent Entity class and several handlers.
The biologically relevant abstraction classes (Gene, CellLine, Cancer, Organelle, GeneCluster, CellLineCluster) inherit their methods from Entity. Entity methods are wrappers for hidden data handler classes who perform specific transformations, such as data indexing and high throughput filtering.

.. image:: _static/imgs/structure.png
   :width: 750

CanDI.candi module
------------------
User facing CanDI classes (Gene, CellLine, Organelle, Cancer, GeneCluster, CellLineCluster) all inherit from a parent Entity Class.

.. automodule:: CanDI.structures.entity
   :members:
   :undoc-members:
   :show-inheritance:

The following are the main CanDI classes. These are what users will use to access and cross reference data.

.. automodule:: CanDI.candi.candi
   :members: Gene, CellLine, Organelle, Cancer, CellLineCluster, GeneCluster
   :undoc-members: SubsetHandler
   :show-inheritance: Gene, CellLine, Organelle, Cancer, CellLineCluster, GeneCluster

CanDI.data module
------------------
The data class is instantiated at import. This class contains paths to all data downloaded with CanDI.
It has internal methods for loading datasets into memory as pandas dataframes.
There are 3 index tables that candi relies on for fetch all data:

- cell_lines
- genes
- locations

These tables are automatically loaded as pandas dataframes upon import of CanDI
It is highly recommended the user familiarize themself with the columns and indexes of these tables.
All candi classes operate through these index tables.

.. automodule:: CanDI.candi.data
   :members:
   :undoc-members:
   :show-inheritance:
.. CanDI documentation master file, created by
   sphinx-quickstart on Wed May 26 13:53:57 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CanDI's documentation!
=================================
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   README

.. image:: _static/imgs/Poster.png
   :width: 400

.. toctree::
   :maxdepth: 4
   :caption: Package:

   CanDI

.. toctree::
   :maxdepth: 2
   :caption: Notebooks:

   get-started.ipynb
   brca_heatmap.ipynb
   kras_egfr_scatter.ipynb
   deseq_setup.ipynb

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
CanDI - A global cancer data integrator
=======================================

|Documentation Status|

Package Installation
--------------------

First, you need to clone this repository to use CanDI.

.. code:: bash

   git clone https://github.com/GilbertLabUCSF/CanDI.git

We suggest to use `Conda <https://docs.conda.io/en/latest/>`__ as a
package manager and environment management system. You can create a
fresh conda environment with all ``CanDI``\ ’s requirements using bellow
command:

.. code:: bash

   conda env create -f CanDI/candi.yml -n candi

Prepare Datasets
~~~~~~~~~~~~~~~~

The python command from CanDI will automatically download and modify
datasets.

.. code:: bash

   python CanDI/CanDI/setup/install.py

Downloaded and formatted datasets would organize this way:

.. code::

   .
   ├── config.ini # modified after Installation 
   ├── depmap
   │   ├── CCLE_expression.csv
   │   ├── CCLE_fusions.csv
   │   ├── CCLE_gene_cn.csv
   │   ├── CCLE_mutations.csv
   │   ├── CCLE_RNAseq_reads.csv
   │   ├── CRISPR_gene_dependency.csv
   │   ├── CRISPR_gene_effect.csv
   │   └── sample_info.csv
   ├── genes
   │   └── gene_info.csv
   └── locations
       └── merged_locations.csv

Package Usage
-------------

Import CanDI into python
~~~~~~~~~~~~~~~~~~~~~~~~

To import ``CanDI``, your active directory in python must be same as the
cloned folder.

.. code:: python

   from CanDI import candi

**OR**, you can add path to the `CanDI` directory if you want to use it from other directories.

.. code:: python

   import sys
   sys.path.append("path-to-candi-directory")

   from CanDI import candi

CanDI Objects
~~~~~~~~~~~~~

-  ``data`` : Container for all candi datasets. All access to datasets
   go through data object.
-  ``Gene`` : Provides cross dataset indexing from the gene perspective.
-  ``CellLine`` : Provides cross dataset indexing from the cell line
   perspective.
-  ``Cancer`` : Provides cross dataset indexing by a group of cell lines
   that are all the same tissue.
-  ``Organelle``: Provides cross dataset indexing for a group of genes
   whose proteins localize to the same organelle.
-  ``CellLineCluster`` : Provides cross dataset indexing for a group of
   user defined cell lines.
-  ``GeneCluster`` : Provides cross dataset indexing for a group of user
   defined genes.

.. |Documentation Status| image:: https://readthedocs.org/projects/candi/badge/?version=latest
   :target: https://candi.readthedocs.io/en/latest/?badge=latest
