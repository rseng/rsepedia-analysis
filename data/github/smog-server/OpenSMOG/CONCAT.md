========
OpenSMOG
========

|Citing OpenSMOG|
|PyPI|
|conda-forge|
|ReadTheDocs|
|SMOG server|
|Update|
|GitHub-Stars|

.. |Citing OpenSMOG| image:: https://img.shields.io/badge/cite-OpenSMOG-informational
   :target: https://opensmog.readthedocs.io/en/latest/Reference/citing.html
.. |PyPI| image:: https://img.shields.io/pypi/v/OpenSMOG.svg
   :target: https://pypi.org/project/OpenSMOG/
.. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/OpenSMOG.svg
   :target: https://anaconda.org/conda-forge/OpenSMOG
.. |ReadTheDocs| image:: https://readthedocs.org/projects/opensmog/badge/?version=latest
   :target: https://opensmog.readthedocs.io/en/latest/
.. |SMOG server| image:: https://img.shields.io/badge/SMOG-Server-informational
   :target: https://smog-server.org/
.. |Update| image:: https://anaconda.org/conda-forge/opensmog/badges/latest_release_date.svg
   :target: https://anaconda.org/conda-forge/opensmog
.. |GitHub-Stars| image:: https://img.shields.io/github/stars/junioreif/OpenSMOG.svg?style=social
   :target: https://github.com/junioreif/OpenSMOG


`Documentation <https://opensmog.readthedocs.io/>`__
| `Install <https://opensmog.readthedocs.io/en/latest/GettingStarted/install.html>`__
| `Tutorials <https://opensmog.readthedocs.io/en/latest/Tutorials/SBM_CA.html>`__

Overview
========

OpenSMOG is a Python library for performing molecular dynamics simulations using Structure-Based Models. OpenSMOG uses the  `OpenMM <http://openmm.org/>`_. Python API, which supports a wide variety of potential energy functions, including those that are commonly employed in C-alpha and all-atom models.
While it is possible to use this library in a standalone fashion, it is expected that users will generate input files using the SMOG2 software (version 2.4, or later, with the flag :code:`-OpenSMOG`). Details on how to generate OpenSMOG-compatible force fields files can be found in the `SMOG2 User Manual <https://smog-server.org/smog2/>`__.

.. raw:: html

    <p align="center">
    <img align="center" src="./docs/source/images/OpenSMOG_pipeline.jpg">
    </p>



Citation
========

When using **OpenSMOG** and **SMOG2**, please `use the following references
<https://opensmog.readthedocs.io/en/latest/Reference/citing.html>`__.



Installation
============

The **OpenSMOG** library can be installed via `conda <https://conda.io/projects/conda/>`_ or `pip <https://pypi.org/>`_, or compiled from `source (GitHub) <https://github.com/junioreif/OpenSMOG>`_.

Install via conda
-----------------

The code below will install **OpenSMOG** from `conda-forge <https://anaconda.org/conda-forge/OpenSMOG>`_.

.. code-block:: bash

    conda install -c conda-forge OpenSMOG

Install via pip
-----------------

The code below will install **OpenSMOG** from `PyPI <https://pypi.org/project/OpenSMOG/>`_.

.. code-block:: bash

    pip install OpenSMOG

Install OpenMM
--------------

The **OpenSMOG** library uses `OpenMM <http://openmm.org/>`_ API to run the molecular dynamics simulations.
These requirements can be met by installing the following packages from the `conda-forge channel <https://conda-forge.org/>`__:

.. code-block:: bash

    conda install -c conda-forge openmm
    
The following are libraries **required** for installing **OpenSMOG**:

    - `Python <https://www.python.org/>`__ (>=3.6)
    - `NumPy <https://www.numpy.org/>`__ (>=1.14)
    - `lxml <https://lxml.de/>`__ (>=4.6.2)

Installing SMOG2
================

The inputs **OpenSMOG** simulations are generated using `SMOG2 <https://smog-server.org/smog2>`_ (version 2.4 and later). Here, there is a quick installation guide based on `conda <https://conda.io/projects/conda/>`_ (Linux and Windows-WSL only).

First, download SMOG 2 (v2.4, or later) at `smog-server.org <https://smog-server.org/smog2/>`__

Second, create a new conda environment and activate it (called smog2.4, but name as appropriate):

.. code-block:: bash

    conda create --name smog2.4 perl
    
.. code-block:: bash

    conda activate smog2.4

Next, it is necessary to install a few **Perl** modules:

.. code-block:: bash

    conda install -c bioconda perl-xml-simple perl-xml-libxml java-jdk

.. code-block:: bash

    conda install -c eumetsat perl-pdl

.. code-block:: bash

    perl -MCPAN -e 'install XML::Validator::Schema'

Add the **Perl** and **smog2** path into the configure.smog2 file (described in the README that comes with SMOG 2).

.. hint:: Use the following command line to find out which installed **Perl** is being used.

.. code-block:: bash

    which perl

Then load and test the **smog2** installation:

.. code-block:: bash

    source configure.smog2
    
.. code-block:: bash

    ./test-config
    
As described in the SMOG 2 manual, it is highly recommended that you also download smog-check and run all checks before using the SMOG 2 software.


Resources
=========

- `Reference Documentation <https://opensmog.readthedocs.io/>`__: Examples, tutorials, and class details.
- `Installing OpenSMOG <https://opensmog.readthedocs.io/en/latest/GettingStarted/install.html#installing-opensmog>`__: Instructions for installing **OpenSMOG**.
- `Installing SMOG2 <https://opensmog.readthedocs.io/en/latest/GettingStarted/install.html#installing-smog2>`__: Instructions for installing **SMOG2**.
- `Issue tracker <https://github.com/smog-server/OpenSMOG/issues>`__: Report issues/bugs or request features.
Table of Contents
=================

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   GettingStarted/install
   GettingStarted/intro

.. toctree::
   :maxdepth: 2
   :caption: OpenSMOG Modules

   OpenSMOG

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   Tutorials/SMOG2_usage
   Tutorials/SBM_CA
   Tutorials/SBM_AA
   
.. toctree::
   :maxdepth: 2
   :caption: Reference

   Reference/citing
   Reference/references
   Reference/license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
OpenSMOG
========

.. automodule:: OpenSMOG
   :members:
   :undoc-members:
   :show-inheritance:
   .. _SMOG2_usage:

=====================================================================
Using SMOG2 to generate  C-alpha and All-Atom Structure-Based models
=====================================================================


This tutorial should take between 5 to 10 minutes to complete. Here, we will use the **SMOG2** software pacakge to generate the SBM (Structure-Based Model) input files that will be used to perform a simulation with **OpenSMOG**. To install SMOG2, please check the installation notes in the SMOG 2 user manual, or use the guide `here <https://opensmog.readthedocs.io/en/latest/GettingStarted/install.html#installing-smog2>`_. Details of SMOG2 usage and options are described in the `manual <https://smog-server.org/smog2/>`_. It is assumed that the executable **smog2** is in your path.


Preparing your PDB file
===============================

The following instructions will use a PDB file of CI2 protein (`2ci2.pdb <https://www.rcsb.org/structure/2CI2>`_).

First, download the PDB file:

.. code-block:: bash

    wget http://www.rcsb.org/pdb/files/2ci2.pdb


Then, it is necessary to clean up the file and only keep information needed to define a structure-based model. In this case, let us keep only the ATOM lines:

.. code-block:: bash

    grep "^ATOM" 2ci2.pdb > 2ci2.atoms.pdb

.. note:: Sometimes, you also want HETATMs. This is up to the user. HETATMs can be things that we don't want to include (e.g. HOH), or things that we may want to included (e.g. posttranslational modifications). In this case, we only want ATOM lines.


Next, add an END line to the file 2ci2.atoms.pdb:

.. code-block:: bash

    sed -i -e "\$aEND" 2ci2.atoms.pdb

Adjust the file, so that the naming convention conforms with the default SMOG models: 

.. code-block:: bash

    smog_adjustPDB -i 2ci2.atoms.pdb -o 2ci2.adj.pdb

Generate OpenSMOG input files for a C-alpha model
==================================    
Use the adjusted file to generate your input CA model:

.. code-block:: bash

    smog2 -i 2ci2.adj.pdb -CA -dname 2ci2.CA -OpenSMOG

Generate OpenSMOG input files for an all-atom model
==================================

To generate input files for the all-atom model, you only need to change the flag -CA to -AA:

.. code-block:: bash

    smog2 -i 2ci2.adj.pdb -AA -dname 2ci2.AA -OpenSMOG

.. note:: When running the simulation in OpenSMOG, there are differences in the simulation protocols and settings. For example, in the case of AA, the cutoff is typically much shorter than the values used with the CA model. However, larger timesteps can typically be used with the AA model. Please, check the `C-alpha <https://opensmog.readthedocs.io/en/latest/Tutorials/SBM_CA.html>`_  and `All-Atom <https://opensmog.readthedocs.io/en/latest/Tutorials/SBM_AA.html>`_ simulation tutorial pages.
========================
How to cite **OpenSMOG**
========================

When using the SMOG2-OpenSMOG framework for publications, please use the following citations: 


.. code-block:: bibtex

    @article{SMOG2, 
    title = {SMOG 2: A Versatile Software Package for Generating Structure-Based Models},
    author = {Noel, Jeffrey K. and 
    Levi, Mariana and 
    Raghunathan, Mohit and 
    Lammert, Heiko and 
    Hayes, Ryan L. and 
    Onuchic, Jos{\'{e}} N. and 
    Whitford, Paul C.},
    journal = {PloS Comp Biol},
    pages = {e1004794},
    volume = {12},
    year = {2016}
    }

    @article{OpenSMOG,
    title={SMOG 2 and openSMOG: Extending the limits of structure-based models},
    author={Oliveira Jr, Antonio B  and Contessoto, Vinicius G and 
    Hassan, Asem and
    Byju, Sandra and
    Wang, Ailun and
    Wang, Yang and
    Dorero-Rojas, Esteban and
    Mohanty, Udayan and
    Noel, Jeffrey K and
    Onuchic, Jose N and
    Whitford, Paul C},
    journal={bioRxiv},
    doi={10.1101/2021.08.15.456423}
    }   
    
==========
References
==========

.. bibliography:: :all: 
    :style: unsrt
    :list: bullet.. _install:

============
Installation
============

Installing OpenSMOG
===================

The **OpenSMOG** library can be installed via `conda <https://conda.io/projects/conda/>`_ or `pip <https://pypi.org/>`_, or compiled from `source (GitHub) <https://github.com/smog-server/OpenSMOG>`_.

Install via conda
-----------------

The code below will install **OpenSMOG** from `conda-forge <https://anaconda.org/conda-forge/OpenSMOG>`_.

.. code-block:: bash

    conda install -c conda-forge OpenSMOG

Install via pip
-----------------

The code below will install **OpenSMOG** from `PyPI <https://pypi.org/project/OpenSMOG/>`_.

.. code-block:: bash

    pip install OpenSMOG

Install OpenMM
--------------

The **OpenSMOG** library uses `OpenMM <http://openmm.org/>`_ API to run the molecular dynamics simulations.
**OpenMM**  may be installed from the `conda-forge channel <https://conda-forge.org/>`__:

.. code-block:: bash

    conda install -c conda-forge openmm
    
The following libraries are **required** when installing **OpenSMOG**:

    - `Python <https://www.python.org/>`__ (>=3.6)
    - `NumPy <https://www.numpy.org/>`__ (>=1.14)
    - `ElementTree XML <https://docs.python.org/3/library/xml.etree.elementtree.html>`__ (>=2.2.0)

Installing SMOG2
================

The input files for **OpenSMOG** simulations are generated using `SMOG2 <https://smog-server.org/smog2>`_ (version 2.4, or newer). Here, there is a quick installation guide based on `conda <https://conda.io/projects/conda/>`_ (Linux and Windows-WSL only). Alternate installation options are described in the SMOG2 manual.

First, create a new environment and activate it:

.. code-block:: bash

    conda create --name smog2.4 perl
    
.. code-block:: bash

    conda activate smog2.4

Next, it is necessary to instal a few **Perl** modules (a complete list of modules is in the SMOG2 README file):

.. code-block:: bash

    conda install -c bioconda perl-xml-simple perl-xml-libxml java-jdk

.. code-block:: bash

    conda install -c eumetsat perl-pdl

.. code-block:: bash

    perl -MCPAN -e 'install XML::Validator::Schema'

Enter the **Perl** and **smog2** paths in the configure.smog2 file.

Then load and test your **smog2** configuration:

.. code-block:: bash

    source configure.smog2
    
.. code-block:: bash

    ./test-config
    
It is also **STRONGLY** recommended that you download **smog-check** (available at smog-server.org) and run all tests before using **smog2** for production calculations.
.. _intro:

============
Introduction
============

OpenSMOG is a Python library for performing molecular dynamics simulations using Structure-Based Models :cite:p:`noel2016smog`. OpenSMOG uses the `OpenMM <http://openmm.org/>`_ Python API that supports a wide variety of potential forms, which includes the commonly employed C-alpha :cite:p:`clementi2000topological` and All-Atom :cite:p:`whitford2009all` models.
The input files are generated using the SMOG2 software package with the flag :code:`-OpenSMOG`. Details on SMOG2 usage can be found in the `SMOG2 User Manual <https://smog-server.org/smog2/>`__.

.. image:: ../images/OpenSMOG_pipeline.jpg

.. bibliography::
    :style: unsrt
