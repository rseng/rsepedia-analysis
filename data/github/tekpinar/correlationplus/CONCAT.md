[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/correlationplus)](https://pypi.org/project/correlationplus/)
[![PyPI](https://img.shields.io/pypi/v/correlationplus)](https://pypi.org/project/correlationplus/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/correlationplus/README.html)
[![Open Source License: GPL v3](https://img.shields.io/badge/License-LGPLv3-blue.svg)](https://opensource.org/licenses/LGPL-3.0)
[![Doc](https://readthedocs.org/projects/correlationplus/badge/?version=latest)](http://correlationplus.readthedocs.org/en/latest/#)
[![Docker Image Version (tag latest semver)](https://img.shields.io/docker/v/structuraldynamicslab/correlationplus/latest)](https://hub.docker.com/repository/docker/structuraldynamicslab/correlationplus)
![Conda](https://img.shields.io/conda/pn/bioconda/correlationplus)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/tekpinar/correlationplus/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/tekpinar/correlationplus)

# correlationplus

A Python package to calculate, visualize and analyze correlation maps of proteins.

correlationplus contains four main scripts that you can use to calculate, visualize
and analyze correlation maps of proteins. 
These correlations can be dynamical cross-correlations, linear mutual
information, sequence based covariation/coevolution or any other pairwise coupling metric. 
You can use elastic network models or your molecular dynamics trajectories for calculation 
of dynamical correlations.  

## Video Tutorials
* Installation: https://www.youtube.com/watch?v=Fc_xpnrbbWU
* Calculate Module: https://www.youtube.com/watch?v=04b7mdulHW8
* Visualize Module: https://www.youtube.com/watch?v=HgUVAV1unXs
* Analyze Module: https://www.youtube.com/watch?v=04BwJDauOn4
* Sequence Conservation Analysis with CoeViz and correlationplus (Part 1): https://www.youtube.com/watch?v=ieu91glEd0s
* Sequence Conservation Analysis with CoeViz and correlationplus (Part 2): https://www.youtube.com/watch?v=BuxQNFoid1A

## Installation

We recommend one of the installation methods for regular users:


### with pip

The pip version required by some dependencies is >= 21.0.1, which is not the pip version bundle with python 3.(6,7,8)
So, you have to update pip before installing *correlationplus*. Otherwise, you will have trouble during *MDAnalysis* dependency installation.
For this reason, we **strongly** encourage you to install correlationplus in a [virtualenv](https://virtualenv.pypa.io/en/latest/).

```bash
python3 -m venv correlationplus
cd correlationplus
source bin/activate
python3 -m pip install -U pip
python3 -m pip install correlationplus
```

if you want to install it without using a virtualenv
and encounter an error related to llvmlite, you can
solve it as follows:
```bash
python3 -m pip install llvmlite --ignore-installed
python3 -m pip install correlationplus
```

or if you do not have administration rights
```bash
python3 -m pip install --user llvmlite --ignore-installed
python3 -m pip install --user correlationplus
```

### with conda
```bash
conda install -c bioconda correlationplus

```

Most of the time, at least one these methods will be sufficient for the installation.
However, if these two methods didn't work for any reason, you can take a look 
to 'Advanced Installation' instructions at
https://correlationplus.readthedocs.io/en/latest/installation.html#advanced-installation.


## Quickstart
There are four main scrips: 
* calculate
* visualize
* analyze
* paths

You can find details of each script with the following commands:

```bash
correlationplus calculate -h
correlationplus visualize -h
correlationplus analyze -h
correlationplus paths -h
```

Go to our [readthedocs](https://correlationplus.readthedocs.io/en/latest/quickstart.html) for 
detailed usage examples for each script.

## Cite
If you use correlationplus, please cite us:

**Extracting Dynamical Correlations and Identifying Key Residues for Allosteric Communication in Proteins by correlationplus**
Mustafa Tekpinar, Bertrand Neron, and Marc Delarue
Journal of Chemical Information and Modeling Article ASAP
[DOI: 10.1021/acs.jcim.1c00742](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00742)


## Licensing

*correlationplus* is developed and released under [GNU Lesser GPL Licence](https://www.gnu.org/licenses/lgpl-3.0.en.html). 
Please read to the **COPYING** and **COPYING.LESSER** files to know more. 
correlationplus.scripts package
===============================

Submodules
----------

correlationplus.scripts.analyze module
--------------------------------------

.. automodule:: correlationplus.scripts.analyze
   :members:
   :undoc-members:
   :show-inheritance:

correlationplus.scripts.calculate module
----------------------------------------

.. automodule:: correlationplus.scripts.calculate
   :members:
   :undoc-members:
   :show-inheritance:

correlationplus.scripts.correlationplus module
----------------------------------------------

.. automodule:: correlationplus.scripts.correlationplus
   :members:
   :undoc-members:
   :show-inheritance:

correlationplus.scripts.diffMap module
--------------------------------------

.. automodule:: correlationplus.scripts.diffMap
   :members:
   :undoc-members:
   :show-inheritance:

correlationplus.scripts.visualize module
----------------------------------------

.. automodule:: correlationplus.scripts.visualize
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: correlationplus.scripts
   :members:
   :undoc-members:
   :show-inheritance:
Quickstart
==========

There are four main scrips of correlationplus package:

* calculate
* visualize
* analyze
* paths

You can find more information about each script as follows::

    correlationplus calculate -h

    correlationplus visualize -h

    correlationplus analyze -h
    
    correlationplus paths -h

Here are some examples of the correlationplus commandline interface.
You can find all required files in the examples folder at our `github page <https://github.com/tekpinar/correlationplus>`_

Correlation Data Types
----------------------
Correlationplus can handle the following correlation/coupling data types:

* dcc: Dynamical cross-correlations in full matrix format.
* ndcc: Normalized dynamical cross-correlations in full matrix format.
* absndcc: Absolute values normalized dynamical cross-correlations in full matrix format
* omegacc: Normalized Pearson correlations of backbone dihedral angle omega in full matrix format.
* phicc: Normalized Pearson correlations of backbone dihedral angle phi in full matrix format.
* psicc: Normalized Pearson correlations of backbone dihedral angle psi in full matrix format.
* lmi: Linear mutual information in full matrix format or output of g_correlation program.
* nlmi: Normalized linear mutual information in full matrix format or output of g_correlation program. 
* coeviz: After removing the header lines, the data is in matrix format. 
* evcouplings: The sequence coupling csv files obtained from https://evcouplings.org/ can be parsed directly. 
* generic: If you have some coupling data (from dynamics, sequences or any other data) in full matrix format, select this option. 

You can control data types in correlationplus script with '-t' option. For example, '-t generic' will tell the script that it is 
a generic data, etc.


**calculate** script
--------------------
With this module, you can calculate dynamical cross-correlation and linear mutual information from
elastic network models and molecular dynamics trajectories. 

The only input file needed is a PDB file. The file can contain single or multiple chains. In all of 
the computations via script interfaces, only Calpha atoms are selected and used.    

To calculate **normalized dynamical cross-correlations** with **Gaussian** network model::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb -m GNM -o ndcc-6lu7-gnm.dat

To calculate **normalized dynamical cross-correlations** with **Anisotropic** network model::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb -m ANM -o ndcc-6lu7-anm.dat

To calculate **normalized dynamical cross-correlations** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr\
			                      -o ndcc-6lu7-md.dat

To calculate **normalized linear mutual informations** with **Anisotropic** network model::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb -t nlmi -o nlmi-6lu7-anm.dat

To calculate **normalized linear mutual informations** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr -t nlmi\
			                      -o nlmi-6lu7-md.dat
To calculate **normalized Pearson cross-correlations of backbone omega dihedral angles** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr\
			    -t omegacc -o omegacc-6lu7-md.dat

To calculate **normalized Pearson cross-correlations of backbone phi dihedral angles** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr\
			    -t phicc -o phicc-6lu7-md.dat

To calculate **normalized Pearson cross-correlations of backbone psi dihedral angles** from a molecular dynamics trajectory (in dcd, xtc or trr format)::

  correlationplus calculate -p 6lu7_dimer_with_N3_protein_sim1.pdb \
                            -f 6lu7_dimer_with_N3_protein_sim1_short.trr\
			    -t psicc -o psicc-6lu7-md.dat

Sometimes, there are not dihedral angles for some residues at the beginning/end of the chains (See https://userguide.mdanalysis.org/1.1.1/examples/analysis/structure/dihedrals.html for details). If there are some missing atoms, you may not also be able to calculate some dihedral angles. To avoid these problems, we fill these missing values with 1.0 to maintain a one-to-one correspondence with the number of residues. Due to this reason, pay attention to the highly correlated values that may be due to this artificial filling. 

**visualize** script
--------------------
Visualize module plots all 2D correlation maps. It prepares tcl and pml files so that the correlated residue pairs can
be visualized with the help of **VMD** and **PyMOL** programs. This interface needs only a pdb file with N residues and
a square matrix of NxN. The correlation data has to be in matrix format, where only A(i,j) values are 
listed in a square matrix format. LMI matrices produced by g_correlation program of Lange and Grubmuller
can also be parsed. 


To run a simple example of visualization, you can use the data and pdb files in the examples folder::

  correlationplus visualize -i ndcc-6lu7-anm.dat -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t absndcc -v 0.75

In addition, the command above will produce plots of absolute values of dynamical cross correlations vs interresidue distances.
This information can be quite useful if you are particularly looking for long-distance interactions. 

The visualize app will produce an output for overall structure and all individual intra-chain correlations, if exist. 
Moreover, the program will give you inter-chain correlations, if you have more than one chain. 

You can analyze the correlations with `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ just by loading the tcl files produced by 
visualize app.  To reduce the clutter, the command above will only dump the correlations greater than 0.75 to your tcl or pml file.
If you would like to visualize an interval, you can specify the maximal value as well with '-x ' parameter.

You can call VMD and go to *Extensions->Tk Console menu*. 
Write the following command to see the correlations::

  source correlation-interchain-chainsA-B.tcl

If you prefer to do the tcl loading in a single command::

  vmd -e correlation-interchain-chainsA-B.tcl

Please, beware that the loading can take some time depending on protein size,
number of correlations and the min-max correlation limits that you imposed. 

Additionally, vmd command has to be in your path if you want to do this 
with the command above.

If you would like to use PyMOL, the following command will be sufficient::
  
  pymol correlation-interchain-chainsA-B.pml



Sometimes, we may need to plot difference map of two correlation maps. 
You may want to see the differences of linear mutual information 
maps of produced with two different methods, conditions etc. The correlations
of ligated vs unligated simulations are some common examples.  
The difference maps can be produced with diffMap app as follows::

  correlationplus diffMap -i 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
                          -j 6lu7_dimer_no_ligand_protein_sim1-lmi.dat\
			  -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

**analyze** script
------------------
This module can be used to perform centrality analysis of the correlation maps.
Centrality analysis is used to deduce active sites, binding sites, key mutation
sites and allosteric residues. 

The script can compute degree, closeness, betweenness, current flow closeness, 
current flow betweenness, eigenvector centralities and major communities. The following
command will do all of the above analysis::

  correlationplus analyze -i 6lu7_dimer_with_N3_protein_sim1-lmi.dat\
                          -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

If you would like to calculate only a certain centrality like betweenness::

  correlationplus analyze -i 6lu7_dimer_with_N3_protein_sim1-lmi.dat\
                          -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb
			  -c betweenness -t lmi

After the calculation, the centrality values will be inserted into **Bfactor** 
column of a pdb file. You can load the pdb files with your favorite visualization 
software and color according to **Bfactors**. If you prefer **VMD** - as we do-, 
the app will produce tcl files so that you can visualize the key residues with **VMD**.
The tcl script highlights the residues with the highest 10% of the selected centrality
in VDW representation.::

  vmd -e correlation_degree.tcl

With PyMol::
  
  pymol correlation_degree.pml

**paths** script
------------------
To calculate suboptimal paths between two active site residues in chain A and chain B of 
SARS-CoV2 main protease::

    correlationplus paths -i ndcc-6lu7-anm.dat\
              		  -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb\
              		  -b A41 -e B41
   
This command will only produce the optimal path and print out the path length. If you would like
to calculate suboptimal paths as well, you can append -n argument. Here is the example command to 
calculate 10 paths between residue 41 of chain A and residue 41 of chain B::

    correlationplus paths -i ndcc-6lu7-anm.dat\
              		  -p 6lu7_dimer_with_N3_protein_sim1_ca.pdb\
              		  -b A41 -e B41 -n 10



Ipython Interface
-----------------
For a detailed analysis, script interfaces provided by calculate, visualize, analyze, paths and 
diffMap scripts may not be sufficient. Therefore, you can use IPython 
to load the modules and do a detailed analysis as follows. 

``from correlationplus.calculate import *``

``from correlationplus.visualize import *``
 
You can get help for each function with

``help(intraChainCorrelationMaps)``

You can check different valueFilters, distanceFilters for your analysis. 
Also, you can scan a range of values by calling the functions in a 
loop. 

There is a minor but important difference between the scripts and the modules for centrality and path 
analyses. If you want to use the module for centrality analysis:

``from correlationplus.centralityAnalysis import *``

Please notice that the name of the script was 'analyze' but the name of the module is 'centralityAnalysis'. 

Similarly, the name of the path analysis script is 'paths' while the name of the module is 'pathAnalysis'. 
Therefore, you have to call path analysis module interactively as follows:

``from correlationplus.pathAnalysis import *``





correlationplus package
=======================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   correlationplus.scripts

Submodules
----------

correlationplus.calculate module
--------------------------------

.. automodule:: correlationplus.calculate
   :members:
   :undoc-members:
   :show-inheritance:

correlationplus.centralityAnalysis module
-----------------------------------------

.. automodule:: correlationplus.centralityAnalysis
   :members:
   :undoc-members:
   :show-inheritance:

correlationplus.visualize module
--------------------------------

.. automodule:: correlationplus.visualize
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: correlationplus
   :members:
   :undoc-members:
   :show-inheritance:
Installation
============

Basic Installation
------------------
We recommend installing correlationplus with pip or conda for regular users:

with pip
~~~~~~~~
We **strongly** encourage you to install correlationplus in a virtualenv::

	python3 -m venv correlationplus
	cd correlationplus
	source bin/activate

Then, you can upgrade pip and install it as follows::

	python3 -m pip install -U pip
	python3 -m pip install correlationplus


or if you do not have administration rights::

	python3 -m pip install --user -U pip
	python3 -m pip install --user correlationplus

If you want to install it without using a virtualenv and encounter an error related to llvmlite, 
you can solve it as follows::

	python3 -m pip install llvmlite --ignore-installed
	python3 -m pip install correlationplus

with conda
~~~~~~~~~~

You can also install correlationplus with conda as follows::

    conda install -c bioconda correlationplus
    
Most of the time, at least one these methods will be sufficient for the installation.
However, if these two methods didn't work for any reason, you can take a look 
at the 'Advanced Installation' instructions to install it from the source.

Advanced Installation
---------------------
If standard installation procedure didn't work for you for any reason, you can 
try one of the methods detailed below:

for developers
~~~~~~~~~~~~~~
The pip version required by some dependencies is >= 21.0.1, which is not the pip version bundled with python 3.(6,7,8)
So, you have to update pip before installing *correlationplus*. Otherwise, you will have trouble during *MDAnalysis* dependency installation.
For this reason, we **strongly** encourage you to install correlationplus in a `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ ::

	python3 -m venv correlationplus
	cd correlationplus
	source bin/activate
	python3 -m pip install -U pip
	python3 -m pip install correlationplus

If you want to install the latest version from the source ::

	python3 -m venv correlationplus
	cd correlationplus
	source bin/activate
	python3 -m pip install -U pip
	mkdir src
	cd src
	git clone https://github.com/tekpinar/correlationplus.git # or git@github.com:tekpinar/correlationplus.git``
	cd correlationplus
	pip install -e .

from Docker image
~~~~~~~~~~~~~~~~~

If you don't want to install or can't install correlatioplus for any reason, you can use docker images. These images do not require any installation. 
You can find docker images for correlationplus  at `<https://hub.docker.com/r/structuraldynamicslab/correlationplus>`_

The computation inside the container will be performed under **correlationplus** id in **/home/correlationplus** directory.
So before running a **correlationplus** container,
do not forget to create and mount a shared directory in the container. 

This directory must be writable to this user. So, you have two possibilities:

1. Map your id on the host to the *correlationplus* user in the container::

       ``-u $(id -u ${USER}):$(id -g ${USER})``
2. Make the shared directory writable to anyone::

       ``chmod 777 shared_dir``

option 1

.. code-block:: shell

    mkdir shared_dir
    cp 6lu7_dimer_no_ligand_protein_sim1-lmi.dat shared_dir
    cp 6lu7_dimer_with_N3_protein_sim1-lmi.dat shared_dir
    cp 6lu7_dimer_with_N3_protein_sim1_ca.pdb shared_dir
    cd shared_dir
    docker run -v $PWD:/home/correlationplus -u $(id -u ${USER}):$(id -g ${USER}) structuraldynamicslab/correlation_plus diffMap\
    				-i 6lu7_dimer_no_ligand_protein_sim1-lmi.dat \
				-j 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
				-p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

option 2

.. code-block:: shell

    mkdir shared_dir
    cp 6lu7_dimer_no_ligand_protein_sim1-lmi.dat shared_dir
    cp 6lu7_dimer_with_N3_protein_sim1-lmi.dat shared_dir
    cp 6lu7_dimer_with_N3_protein_sim1_ca.pdb shared_dir
    chmod 777 shared_dir
    cd shared_dir
    docker run -v $PWD:/home/correlationplus structuraldynamicslab/correlation_plus diffMap \
    						-i 6lu7_dimer_no_ligand_protein_sim1-lmi.dat \
						-j 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
						-p 6lu7_dimer_with_N3_protein_sim1-lmi.dat -t lmi


It is also possible to run an ipython interactive session::

    docker run -v $PWD:/home/correlationplus --entrypoint /bin/bash -it structuraldynamicslab/correlationplus:0.1.4rc2

then once in the container

``ipython``

from Singularity image
~~~~~~~~~~~~~~~~~~~~~~

As the docker image is registered in dockerhub you can also use it directly with `Singularity <https://sylabs.io/docs/>`_ ::

    singularity run docker://structuraldynamicslab/correlationplus diffMap \
    					-i 6lu7_dimer_no_ligand_protein_sim1-lmi.dat \
					-j 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
					-p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

or in 2 steps ::

    singularity pull correlationplus.simg docker://structuraldynamicslab/correlation_plus
    ./correlationplus.simg diffMap \
    			-i 6lu7_dimer_no_ligand_protein_sim1-lmi.dat \
			-j 6lu7_dimer_with_N3_protein_sim1-lmi.dat \
			-p 6lu7_dimer_with_N3_protein_sim1_ca.pdb -t lmi

Unlike docker, you do not have to worry about shared directory, your *home* and */tmp* are automatically shared.
You can also run an *ipython* interactive session ::

    singularity shell correlationplus.simg
correlationplus
===============

.. toctree::
   :maxdepth: 4

   correlationplus
.. correlationplus documentation master file, created by
   sphinx-quickstart on Fri Jan 29 14:29:41 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to correlationplus's documentation!
===========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   license

	   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
