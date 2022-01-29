Developer Guideline
===================

Here we provide some details about guideline for internal developers. Most of the
choices are explained in the [guide](https://guide.esciencecenter.nl).

For a quick reference on software development, we refer to [the software
guide checklist](https://guide.esciencecenter.nl/best_practices/checklist.html).

Version control
---------------

We use [git](http://git-scm.com/) and [github](https://github.com/) for developing
[DeepRank](https://github.com/DeepRank/deeprank). And [Github Desktop](https://desktop.github.com/) is recommennded to simplify development workflow.

To get onboard, first ask project manager to get access to the DeepRank repo (short for repository), then clone the repo to your local machine.

**The development workflow is as the following**:
1. Create an issue in [the repo Issues](https://github.com/DeepRank/deeprank/issues) for any idea, bug, feature, improvement, etc.
2. Comment and mention (i.e. @) someone in the issue for discussion if necessary.
3. Assign the issue to yourself or someone after asking
4. Create a new branch from the `development` branch for this issue, and name it in a format of `issue{ID}_*`, e.g. `issue7_format_docstring`.
5. Work in this issue branch to solve the issue. Check this [guide](https://google.github.io/eng-practices/review/developer/) to see how to make good commits.
6. Create a Pull Request to merge this issue branch to the `development` branch when last step is completed.
7. Assign someone but not yourself to review the pull request. For reviewers, check this [guide](https://google.github.io/eng-practices/review/reviewer/) to see how to do a code review.
8. Follow reviewer's comments to fix the code until reviewer approve you to merge. Check this [guide](https://google.github.io/eng-practices/review/developer/handling-comments.html) to see how to handle reviewer comments.
9. Merge the issue branch to `development` branch and delete the issue branch.
10. Close the issue after leavning a comment of the related pull request ID.

Repeat the 1-10 steps for next issue.

Note that try to keep the issue and pull request small for efficient code review. In principle, it should take ≤30 mins to review a pull request.


Package management and dependencies
-----------------------------------

You can use either `pip` or `conda` for
installing dependencies and package management. This repository does not
force you to use one or the other, as project requirements differ. For
advice on what to use, please check [the relevant section of the
guide](https://guide.esciencecenter.nl/best_practices/language_guides/python.html#dependencies-and-package-management).

-   Dependencies should be added to `setup.py` in the
    `install_requires` list.

Packaging/One command install
-----------------------------

You can distribute your code using pipy or conda. Again, the project
template does not enforce the use of either one. [The
guide](https://guide.esciencecenter.nl/best_practices/language_guides/python.html#building-and-packaging-code)
can help you decide which tool to use for packaging.

If you decide to use pypi for distributing you code, you can configure
travis to upload to pypi when you make a release. If you specified your
pypi user name during generation of this package, the `.travis.yml` file
contains a section that looks like:

``` {.sourceCode .yaml}
deploy:
  provider: pypi
  user: no
  password:
    secure: FIXME; see README for more info
 on:
    tags: true
    branch: master
```

Before this actually works, you need to add an encrypted password for
your pypi account. The [travis
documentation](https://docs.travis-ci.com/user/deployment/pypi/)
specifies how to do this.

Testing and code coverage
-------------------------

-   Tests should be put in the `test` folder.
-   The `test` folder contains:
    -   Example tests that you should replace with your own meaningful
        tests (file: `test_learn.py`)
-   The testing framework used is [PyTest](https://pytest.org)
    -   [PyTest
        introduction](http://pythontesting.net/framework/pytest/pytest-introduction/)
-   Tests can be run with `python setup.py test`
    -   This is configured in `setup.py` and `setup.cfg`
-   Use [Travis CI](https://travis-ci.com/) to automatically run tests
    and to test using multiple Python versions
    -   Configuration can be found in `.travis.yml`
    -   [Getting started with Travis
        CI](https://docs.travis-ci.com/user/getting-started/)
-   TODO: add something about code quality/coverage tool?
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/best_practices/language_guides/python.html#testing)

Documentation
-------------

-   Documentation should be put in the `docs` folder. The contents have
    been generated using `sphinx-quickstart` (Sphinx version 1.6.5).
-   We recommend writing the documentation using Restructured Text
    (reST) and Google style docstrings.
    -   [Restructured Text (reST) and Sphinx
        CheatSheet](http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html)
    -   [Google style docstring
        examples](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).
-   The documentation is set up with the Read the Docs Sphinx Theme.
    -   Check out the [configuration
        options](https://sphinx-rtd-theme.readthedocs.io/en/latest/).
-   To generate html documentation run `python setup.py build_sphinx`
    -   This is configured in `setup.cfg`
    -   Alternatively, run `make html` in the `docs` folder.
-   The `docs/_templates` directory contains an (empty) `.gitignore`
    file, to be able to add it to the repository. This file can be
    safely removed (or you can just leave it there).
-   To put the documentation on [Read the
    Docs](https://readthedocs.org), log in to your Read the Docs
    account, and import the repository (under 'My Projects').
    -   Include the link to the documentation in this [README](https://guide.esciencecenter.nl/best_practices/language_guides/python.html#writingdocumentation).
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/best_practices/language_guides/python.html#writingdocumentation)

Coding style conventions and code quality
-----------------------------------------

-   Check your code style with `prospector`
-   You may need run `pip install .[dev]` first, to install the required
    dependencies
-   You can use `autopep8` to fix the readability of your code style and
    `isort` to format and group your imports
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/best_practices/language_guides/python.html#coding-style-conventions)

Package version number
----------------------

-   We recommend using [semantic
    versioning](https://guide.esciencecenter.nl/best_practices/releases.html#semantic-versioning).
-   For convenience, the package version is stored in a single place:
    `deeprank/__version__.py`. For updating the
    version number, you only have to change this file.
-   Don't forget to update the version number before [making a
    release](https://guide.esciencecenter.nl/best_practices/releases.html)!

Logging
-------

-   We recommend using the `logging` module for getting
    useful information from your module (instead of using
    `print`).
-   The project is set up with a logging example.
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/best_practices/language_guides/python.html#logging)

CHANGELOG.rst
-------------

-   Document changes to your software package
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/software/releases.html#changelogmd)

CITATION.cff
------------

-   To allow others to cite your software, add a `CITATION.cff` file
-   It only makes sense to do this once there is something to cite
    (e.g., a software release with a DOI).
-   Follow the [making software
    citable](https://guide.esciencecenter.nl/citable_software/making_software_citable.html)
    section in the guide.

CODE_OF_CONDUCT.rst
---------------------

-   Information about how to behave professionally
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/software/documentation.html#code-of-conduct)

CONTRIBUTING.rst
----------------

-   Information about how to contribute to this software package
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/software/documentation.html#contribution-guidelines)

MANIFEST.in
-----------

-   List non-Python files that should be included in a source
    distribution
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/best_practices/language_guides/python.html#building-and-packaging-code)

NOTICE
------

-   List of attributions of this project and Apache-license dependencies
-   [Relevant section in the
    guide](https://guide.esciencecenter.nl/best_practices/licensing.html#notice)
# Change Log

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

## v0.2.0
This is the first stable release.# DeepRank
[![PyPI](https://img.shields.io/pypi/v/deeprank)](https://pypi.org/project/deeprank/)
[![Documentation Status](https://readthedocs.org/projects/deeprank/badge/?version=latest)](http://deeprank.readthedocs.io/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3735042.svg)](https://doi.org/10.5281/zenodo.3735042)
![Build](https://github.com/DeepRank/deeprank/workflows/Build/badge.svg)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/4254dd4798bf4cfa9f8f6fe0079de144)](https://www.codacy.com/gh/DeepRank/deeprank/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DeepRank/deeprank&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/DeepRank/deeprank/badge.svg?branch=master)](https://coveralls.io/github/DeepRank/deeprank?branch=master)


### Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Tutorial](#Tutorial)
- [Documentation](https://deeprank.readthedocs.io/)
- [License](./LICENSE)
- [Issues & Contributing](#Issues-and-Contributing)

## Overview
![alt-text](./pics/deeprank.png)

DeepRank is a general, configurable deep learning framework for data mining protein-protein interactions (PPIs) using 3D convolutional neural networks (CNNs).

DeepRank contains useful APIs for pre-processing PPIs data, computing features and targets, as well as training and testing CNN models.

#### Features:

- Predefined atom-level and residue-level PPI feature types
   - *e.g. atomic density, vdw energy, residue contacts, PSSM, etc.*
- Predefined target types
   - *e.g. binary class, CAPRI categories, DockQ, RMSD, FNAT, etc.*
- Flexible definition of both new features and targets
- 3D grid feature mapping
- Efficient data storage in HDF5 format
- Support both classification and regression (based on PyTorch)

## Installation

DeepRank requires a Python version 3.7 or 3.8 on Linux and MacOS.

#### Stable Release

DeepRank is available in stable releases on [PyPI](https://pypi.org/project/deeprank/):
-  Install the module `pip install deeprank`

#### Development Version

You can also install the under development source code from the branch `development`

- Clone the repository `git clone --branch development https://github.com/DeepRank/deeprank.git`
- Go there             `cd deeprank`
- Install the package  `pip install -e ./`

To check if installation is successful, you can run a test
- Go into the test directory `cd test`
- Run the test suite         `pytest`


## Tutorial

We give here the tutorial like introduction to the DeepRank machinery. More informatoin can be found in the documentation <http://deeprank.readthedocs.io/en/latest/>.  We quickly illsutrate here the two main steps of Deeprank:

-   the generation of the data
-   running deep leaning experiments.

### A . Generate the data set (using MPI)

The generation of the data require only require PDBs files of decoys and their native and the PSSM if needed. All the features/targets and mapped features onto grid points will be auomatically calculated and store in a HDF5 file.

```python
from deeprank.generate import *
from mpi4py import MPI

comm = MPI.COMM_WORLD

# let's put this sample script in the test folder, so the working path will be deeprank/test/
# name of the hdf5 to generate
h5file = './hdf5/1ak4.hdf5'

# for each hdf5 file where to find the pdbs
pdb_source = ['./1AK4/decoys/']


# where to find the native conformations
# pdb_native is only used to calculate i-RMSD, dockQ and so on.
# The native pdb files will not be saved in the hdf5 file
pdb_native = ['./1AK4/native/']


# where to find the pssm
pssm_source = './1AK4/pssm_new/'


# initialize the database
database = DataGenerator(
    chain1='C', chain2='D',
    pdb_source=pdb_source,
    pdb_native=pdb_native,
    pssm_source=pssm_source,
    data_augmentation=0,
    compute_targets=[
        'deeprank.targets.dockQ',
        'deeprank.targets.binary_class'],
    compute_features=[
        'deeprank.features.AtomicFeature',
        'deeprank.features.FullPSSM',
        'deeprank.features.PSSM_IC',
        'deeprank.features.BSA',
        'deeprank.features.ResidueDensity'],
    hdf5=h5file,
    mpi_comm=comm)


# create the database
# compute features/targets for all complexes
database.create_database(prog_bar=True)


# define the 3D grid
 grid_info = {
   'number_of_points': [30,30,30],
   'resolution': [1.,1.,1.],
   'atomic_densities': {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8},
 }

# Map the features
database.map_features(grid_info,try_sparse=True, time=False, prog_bar=True)

```

This script can be exectuted using for example 4 MPI processes with the command:

```
    NP=4
    mpiexec -n $NP python generate.py
```

In  the first part of the script we define the path where to find the PDBs of the decoys and natives that we want to have in the dataset. All the .pdb files present in _pdb_source_ will be used in the dataset. We need to specify where to find the native conformations to be able to compute RMSD and the dockQ score. For each pdb file detected in _pdb_source_, the code will try to find a native conformation in _pdb_native_.

We then initialize the `DataGenerator` object. This object (defined in `deeprank/generate/DataGenerator.py`) needs a few input parameters:

-   pdb_source: where to find the pdb to include in the dataset
-   pdb_native: where to find the corresponding native conformations
-   compute_targets: list of modules used to compute the targets
-   compute_features: list of modules used to compute the features
-   hdf5: Name of the HDF5 file to store the data set

We then create the data base with the command `database.create_database()`. This function autmatically create an HDF5 files where each pdb has its own group. In each group we can find the pdb of the complex and its native form, the calculated features and the calculated targets. We can now mapped the features to a grid. This is done via the command `database.map_features()`. As you can see this method requires a dictionary as input. The dictionary contains the instruction to map the data.

-   number_of_points: the number of points in each direction
-   resolution: the resolution in Angs
-   atomic_densities: {'atom_name': vvdw_radius} the atomic densities required

The atomic densities are mapped following the [protein-ligand paper](https://arxiv.org/abs/1612.02751). The other features are mapped to the grid points using a Gaussian function (other modes are possible but somehow hard coded)

#### Visualization of the mapped features

To explore the HDf5 file and vizualize the features you can use the dedicated browser <https://github.com/DeepRank/DeepXplorer>. This tool saloows to dig through the hdf5 file and to directly generate the files required to vizualie the features in VMD or PyMol. An iPython comsole is also embedded to analyze the feature values, plot them etc ....

### B . Deep Learning

The HDF5 files generated above can be used as input for deep learning experiments. You can take a look at the file `test/test_learn.py` for some examples. We give here a quick overview of the process.

```python
from deeprank.learn import *
from deeprank.learn.model3d import cnn_reg
import torch.optim as optim
import numpy as np

# input database
database = '1ak4.hdf5'

# output directory
out = './my_DL_test/'

# declare the dataset instance
data_set = DataSet(database,
            chain1='C',
            chain2='D',
            grid_info={
                'number_of_points': (10, 10, 10),
                'resolution': (3, 3, 3)},
            select_feature={
                'AtomicDensities': {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8},
                'Features': ['coulomb', 'vdwaals', 'charge', 'PSSM_*']},
            select_target='DOCKQ',
            normalize_features = True, normalize_targets=True,
            pair_chain_feature=np.add,
            dict_filter={'DOCKQ':'<1'})


# create the network
model = NeuralNet(data_set,cnn_reg,model_type='3d',task='reg',
                  cuda=False,plot=True,outdir=out)

# change the optimizer (optional)
model.optimizer = optim.SGD(model.net.parameters(),
                            lr=0.001,
                            momentum=0.9,
                            weight_decay=0.005)

# start the training
model.train(nepoch = 50,divide_trainset=0.8, train_batch_size = 5,num_workers=0)
```

In the first part of the script we create a Torch database from the HDF5 file. We can specify one or several HDF5 files and even select some conformations using the `dict_filter` argument. Other options of `DataSet` can be used to specify the features/targets the normalization, etc ...

We then create a `NeuralNet` instance that takes the dataset as input argument. Several options are available to specify the task to do, the GPU use, etc ... We then have simply to train the model. Simple !

## Issues and Contributing

If you have questions or find a bug, please report the issue in the [Github issue channel](https://github.com/DeepRank/deeprank/issues).

If you want to change or further develop DeepRank code, please check the [Developer Guideline](./developer_guideline.md) to see how to conduct further development.
To train on Cartesius with 10-fold CV you need to use the command 
```sbatch --array=1-10 kfold_train-balance_nofilt.slurm```
# Run DeepRank pre-trained models on your own data

We herein provide the best models discussed in the DeepRank paper (doi: 10.1038/s41467-021-27396-0) for the following applications:
- 3DeepFace
- Re-scoring of docking models

Example of codes to run DeepRank on new data with the pre-trained models are provided (*test.py* in each subdirectory).

## Usage: 
`python test.py `

The data used to train the model are available on [SBGrid](https://data.sbgrid.org/dataset/843/)

Further information can be found in the DeepRank [online documentation ](https://deeprank.readthedocs.io/en/latest/)

## Please cite

Renaud, N., Geng, C., Georgievska, S., F. Ambrosetti, L. Ridder, D. F. Marzella, M. F. Réau, A. M. J. J. Bonvin & L. C. Xue . DeepRank: a deep learning framework for data mining 3D protein-protein interfaces.
Nat Commun 12, 7068 (2021). https://doi-org.proxy.library.uu.nl/10.1038/s41467-021-27396-0


---
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
