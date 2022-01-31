# iScore

**Support Vector Machine on Graph Kernel for Ranking Protein-Protein Docking Models**

[![PyPI](https://img.shields.io/pypi/v/iscore)](https://pypi.org/project/iscore/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2630566.svg)](https://doi.org/10.5281/zenodo.2630566)
[![RSD](https://img.shields.io/badge/RSD-iScore-red)](https://research-software.nl/software/iscore)
![Build_Test](https://github.com/DeepRank/iScore/workflows/Build_Test/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/DeepRank/iScore/badge.svg?branch=master)](https://coveralls.io/github/DeepRank/iScore?branch=master)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/3e277331d8fe4de0a22b630bdce6a121)](https://www.codacy.com/gh/DeepRank/iScore/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DeepRank/iScore&amp;utm_campaign=Badge_Grade)
[![Documentation Status](https://readthedocs.org/projects/iscoredoc/badge/?version=latest)](http://iscoredoc.readthedocs.io/?badge=latest)

![alt text](./image/workflow.png)

## 1. Installation

Minimal information to install the module
- Check if command `mpiexec` is available or not in your console. If not, download and install [openmpi](https://www.open-mpi.org/) or [mpich](https://www.mpich.org/).
- Install iScore using `pip install iScore`

Possible problems:
- If `pip install iScore` gives problems on installing `mpi4py`, try to first install `mpi4py` using `conda install mpi4py` and then `pip install iScore`.

## 2. Documentation

The documentation of the package can be found at:
- https://iscoredoc.readthedocs.io


## 3. Quick Examples

iScore offers simple solutions to classify protein-protein interfaces using a support vector machine approach on graph kernels. The simplest way to use iScore is through dedicated binaries that hide the complexity of the approach and allows access to the code with simple command line interfaces. The two binaries are `iscore.train` and `iscore.predict` (`iscore.train.mpi` and `iscore.predict.mpi` for parallel running) that respectively train a model using a training set and use this model to rank the docking models of a protein-protein complex.

### Requirements for preparing data:

 - Use the following file structure
      ```
      root/
      |__train/
      |    |__ pdb/
      |    |__ pssm/
      |    |__ class.lst
      |__test/
            |__pdb/
            |__pssm/
            |__ class.lst (optional)
      ```
      The `pdb` folder contains the PDB files of docking models, and `pssm` contains the PSSM files. The `class.lst` is a list of class ID and PDB file name for each docking model, like `0 7CEI_10w`.

            Check the package subfolders `example/train` and `example/test` to see how to prepare your files.

- PDB files and PSSM files must have consistent sequences.
[PSSMGen](https://github.com/DeepRank/PSSMGen) can be used to get consistent PSSM and PDB files. It is already installed along with iScore. Check [README](https://github.com/DeepRank/PSSMGen) to see how to use it.

### Example 1. Use our trained model

You can directly use our trained model to score your docking conformations. The model we provide is trained on docking benchmark version 4 ([BM4](https://zlab.umassmed.edu/benchmark/)) data, in total 234 different structures were used (117 positive and 117 negative). More details see [this paper](https://doi.org/10.1093/bioinformatics/btz496).
You can find the model in the package subfolder `model/training_set.tar.gz`.

To use this model go into your `test` subfolder and type:

```bash
# Without MPI
iScore.predict

# With MPI
mpiexec -n ${NPROC} iScore.predict.mpi
```

The code will automatically detect the path of the model.

This binary will output the binary class and decision value of the conformations in the test set in a text file `iScorePredict.txt`.

      For the predicted iScore values, the lower value, the better quality of the conformation.


### Example 2. Train your own model

To train the model simply go to your `train` subfolder and type:

```bash
# Without MPI
iScore.train

# With MPI
mpiexec -n ${NPROC} iScore.train.mpi
```

This binary will generate a archive file called by default `training_set.tar.gz` that contains all the information needed to predict binary classes of a test set using the trained model.

To use this model go into your `test` subfolder and type:

```bash
# Without MPI
iScore.predict --archive ../train/training_set.tar.gz

# With MPI
mpiexec -n ${NPROC} iScore.predict.mpi --archive ../train/training_set.tar.gz
```

## 4. Citation

If you use iScore software, please cite the following articles:

1. *Cunliang Geng, Yong Jung, Nicolas Renaud, Vasant Honavar, Alexandre M J J Bonvin, and Li C Xue.* “**iScore: A Novel Graph Kernel-Based Function for Scoring Protein-Protein Docking Models.**” Bioinformatics, 2019, https://doi.org/10.1093/bioinformatics/btz496.
2. *Nicolas Renaud, Yong Jung, Vasant Honavar, Cunliang Geng, Alexandre M. J. J. Bonvin, and Li C. Xue.* “**iScore: An MPI Supported Software for Ranking Protein–Protein Docking Models Based on a Random Walk Graph Kernel and Support Vector Machines.**” SoftwareX, 2020, https://doi.org/10.1016/j.softx.2020.100462.---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**Environment:**
- OS system:
- Version:
- Branch commit ID: 
- Inputs:

**To Reproduce**
Steps/commands to reproduce the behaviour:
 1.  
 2.
 3. 

**Expected Results**
A clear and concise description of what you expected to happen.

**Actual Results or Error Info**
If applicable, add screenshots to help explain your problem.

**Additional Context**
Add any other context about the problem here.
iScore Workflow
========================

One of the mainfeature of the software are the serial and MPI binaries that fully automatize the workflow and that can be used directly from the command line. To illustrate the use of these binaries go to the folder ``iScore/example/training_set/``. This folder contains the subfolders ``pdb/`` and ``pssm/`` that contain the PDB and PSSM files of our training set. The binary class corresponding to these PDBs are specified in the file 'caseID.lst'.

Training a model using iScore can be done in a single line using MPI binaries with the command :


``$ mpiexec -n 2 iScore.train.mpi``


This command will first generate the graphs of the conformations stored in ``pdb/`` using the PSSM contained in ``pssm/`` as features. These graphs will be stored as pickle file  in ``graph/``. The command  will then compute the pairwise kernels of these graphs and store the kernel files in ``kernel/``. Finally it will train a SVM model using the kernel files and the ``caseID.lst`` file that contains the binary class of the model.

The calculated graphs and the svm model are stored in a single tar file called here ``training_set.tar.gz``. This file contains all the information needed to predict binary classes of a test set using the trained model.

To predict binary classes (and decision values) of new conformations go to the subfoler ``test/``. Here 5 conformations are specified by the PDB and PSSM files stored in ``pdb/`` and ``pssm/`` that we want to use as a test set. Ranking these conformations can be done in a single command using :

``$ mpiexec -n 2 iScore.predict.mpi --archive ../training_set.tar.gz``

This command will use first compute the graph of the comformation in the test set and store them in `graph/`. The binary will then compute the pair wise kernels of each graph in the test set with all the graph contained in the training set that are stored in the tar file. These kernels will be stored in ``kernel/``. Finally the binary will use the trained SVM model contained in the tar file to predict the binary class and decision value of the conformations in the test set. The results are then stored in a text file and a pickle file ``iScorePredict.pkl`` and ``iScorePredict.txt``. Opening the text file you will see :

+--------+--------+---------+-------------------+
|Name    |   label|     pred|     decision_value|
+--------+--------+---------+-------------------+
|1ACB_2w |   None |       0 |           -0.994  |
+--------+--------+---------+-------------------+
|1ACB_3w |   None |       0 |           -0.994  |
+--------+--------+---------+-------------------+
|1ACB_1w |   None |       0 |           -0.994  |
+--------+--------+---------+-------------------+
|1ACB_4w |   None |       0 |           -0.994  |
+--------+--------+---------+-------------------+
|1ACB_5w |   None |       0 |           -0.994  |
+--------+--------+---------+-------------------+


The ground truth label are here all None because they were not provided in the test set. This can simply be done by adding a ``caseID.lst`` in the ``test/`` subfolder.


Serial Binaries
------------------------

Serial binaries are also provided and can be used in a similar way than  the MPI binaries : ``iscore.train`` and ``iscore.predict``



Visualizing the connection graphs
======================================

iSore allows to easily visualize the connection graphs using the HDF5 browser provided with the software and pymol. First the connections graphs must be stored in a HDF5 file. To do that simply generate the graphs as following:


>>> from iScore.graphrank.graph import iscore_graph
>>> iscore_graph(pdb_path=<pdb_path>,
>>>              pssm_path=<pssm_path>,
>>>              export_hdf5=True)

where you have to specify the folder containing the PDB files abd PSSM files in pdb_path and pssm_path. By default this are simply ``./pdb/`` and ``./pssm/``. The script above will create a HDF5 file containing the graph.

This HDF5 cile can be explored using the the dedicated HDF5 browser. Go to the ``./h5x/`` folder and type:

``./h5x.py``

This will open the hdf5 browser. You can open a hdf5 file by clicking on the file icon in the bottom  left of the browser. Once opened, you will see the content of the file in the browser. Right-click on the name of a conformation and choose ``3D Plot``. This will open PyMol and allow you to visualize the connecton graph

.. image :: h5x_iscore.pngGPU Kernels
============================

TO DOGraphs and Kernels
===============================


Generating the Graphs :
----------------------------


The first step in iSCore is to generate the connections graph of the itnerface. In this graph each node is represented by the PSSM of a residue. The nodes are connected if they form a contact pair between the two proteins.

To create the graph one needs the PDB file of the interface and the two PSSM files (one for each chain) created by the PSSMGen tool. To generate the graph simply use :

>>> from iScore.graph import GenGraph, Graph
>>> 
>>> pdb = name.pdb
>>> pssm = {'A':'name.A.pdb.pssm','B':'name.B.pdb.pssm'}
>>> 
>>> g = GenGraph(pdb,pssm)
>>> g.construct_graph()
>>> g.export_graph('name.pkl')

This simple example will construct the connection graph and export it in a pickle file. A working example can be found in ``example/graph/create_graph.py``

The function ``iscore_graph()`` facilitate the generation of a large number of conformations. By default this function will create the graphs of all the conformations stored in the subfolder ``./pdb/`` using the pssm files stored in the subfolder ``./pssm/``. The resulting graphs will be stored in the subfolder ``./graph/``.

Generating the Graph Kernels :
-------------------------------------

Once we have calculated the graphs of multiple conformation we can simply compute the kernel of the different pairs using iScore. An example can be found at ``example/kernel/create_kernel.py``

>>> from iScore.graph import Graph, iscore_graph
>>> from iScore.kernel import Kernel
>>> 
>>> # create all the graphs of the pdb in ./pdb/
>>> iscore_graph()
>>> 
>>> #init and load the data
>>> ker = Kernel()
>>> ker.import_from_mat()
>>> 
>>> # run the calculations
>>> ker.run(lamb=1.0,walk=4,check=checkfile)

The kernel between the two graphs computed above is calculated with the class `Kernel()`. By default the method `Kernel.import_from_mat()` will read all the graphs stored in the subfolder `graph/`. To compute all the pairwise kernels of the graphs loaded above we can simply use the method `Kernel.run()`. We can here specify the value of lambda and the length of the walk.Quick example
==================================

iScore is most easily used via a set of binaries that are located in the directory `iscore/bin/`. Thes binaries allows to create graphs, compute kernel traina SVM model and use it in just a few command lines. The only restriction is that they have to be called in a well defined directory structure.

We therefore recomend to have the following file arborescence:

::
    root/
     |__train/
     |    |__pdb/
     |    |__pssm/
     |    |__caseID.lst
     |__predict/
          |__pdb/
          |__pssm/
          |__caseID.lst (optional)


The subfolders are ``train/pdb/`` and ``train/pssm/`` contains the pdb and pssm files we want to use to create a training set. The PSSM files can be obtained via a dedicated Python module `PSSMGen` (see section 'Computing PSSM Files'). To train a SVM model on these conformation simply go in the ``train`` subfolder and use the `iScore.train` binary:


  * ``cd root/train/``
  * ``mpiexec -n 4 iScore.train.mpi``


This binary will automatically generates the graphs of the interfaces, compute their pairwise kernels and train a SVM model using those kernels. The binary will produce an archive file called ``training_set.tar.gz`` that can be used to predict the near native character of the new conformations.

To do that go in the ``root/predict/`` subfolder and use the ``iSore.predcit`` binary:


  * ``cd root/train/``
  * ``mpiexec -n 4 iScore.predict.mpi --archive ../train/training_set.tar.gz``

This will use the training set to evaluate the near-native character of all the conformations present in the test set. The result will be outputed in a text file ``iSorePredict.dat``.Support Vector Machine
==================================

TO DO.Computing PSSM files
=============================

As a prepocessign step one must compute the PSSM files corespondng to the PDB files in the training/testing dataset. Thiscan be acheived with the PisBLast library (https://ncbiinsights.ncbi.nlm.nih.gov/2017/10/27/blast-2-7-1-now-available/). The library BioPython allows ane asy use of these libraries.


iScore contains wrapper that allows to compute the PSSM data, map them to the PDB files and format them for further processing. The only input needed is the PDB file of the decoy. To compute the PSSM file one can simply use :


>>> from iscore.pssm.pssm import PSSM
>>>
>>> gen = PSSM('1AK4')
>>>
>>> # generates the FASTA query
>>> gen.get_fasta()
>>>
>>> # configure the generator
>>> gen.configure(blast=<path to blast binary>, database=<path to the blast db>)
>>>
>>> # generates the PSSM
>>> gen.get_pssm()
>>>
>>> # map the pssm to the pdb
>>> gen.map_pssm().. iScore documentation master file, created by
   sphinx-quickstart on Tue Nov  6 11:07:16 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to iScore's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   intro
   install
   workflow
   pssm
   graph
   viz
   doc

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Documentation
**************

All the prototypes of the class/methods are here specified

.. toctree::
   :maxdepth: 4
   :caption: Documentation:

.. automodule:: iScore


Graph
----------------------------------------

.. automodule:: iScore.graph
    :members:
    :undoc-members:

Kernel
----------------------------------------

.. automodule:: iScore.kernel
    :members:
    :undoc-members:


SVM
----------------------------------------

.. automodule:: iScore.rank
    :members:
    :undoc-members:Installation
==============================

Install via pip
-------------------

In most cases you can use pip to install iScore. Simply use:

`pip install iScore`

Install from source
-------------------

If the pip install fails or if you want to modify the code you can install iScore manually. The code is hosted on Github_ (https://github.com/DeepRank/iScore)

.. _Github: https://github.com/DeepRank/iScore

To install the code

 * clone the repository ``git clone https://github.com/DeepRank/iScore.git``
 * go there ``cd iScore``
 * install the module ``pip install -e ./``

To test the module go to the test folder ``cd ./test`` and execute the following test : ``pytest``




.. highlight:: rst

Introduction
=============================

**Support Vector Machine on Graph Kernels for Protein-Protein Docking Scoring**

The software supports the publication of the following articles:

C. Geng *et al.*, *iScore: A novel graph kernel-based function for scoring protein-protein docking models*, bioRxiv 2018,  https://doi.org/10.1101/498584


iScore uses a support vector machine (SVM) approach to rank protein-protein interfaces. Each interface is represented by a connection graph in which each node represents a contact residue and each edge the connection between two contact residues of different proterin chain. As feature, the node contains the Position Specific Similarity Matrix (PSSM) of the corresponding residue.

To measure the similarity between two graphs, iScore use a random walk graph kernel (RWGK) approach. These RWGKs are then used as input of the SVM model to either train the model on a training set or use a pretrained model to rank new protein-protein interface.

.. image :: comp.png







