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
############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

#. you have a question;
#. you think you may have found a bug (including unexpected behavior);
#. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

#. use the search functionality `here <https://github.com/DeepRank/deeprank/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/DeepRank/deeprank/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
#. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

#. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
#. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
#. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and `here <https://help.github.com/articles/syncing-a-fork/>`__);
#. make sure the existing tests still work by running ``python setup.py test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the DeepRank repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###############################################################################
Contributor Covenant Code of Conduct
###############################################################################

Our Pledge
**********

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
*************

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

Our Responsibilities
********************

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
*****

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
***********

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at n.renaud@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
***********

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
Learning
========

This section describes how to prepare training/validation/test datasets with specific features and targets, how to configure neural network architecture, and how to train and test a neural network with DeepRank.

This tutorial uses the example from the test file ``test/test_learn.py``.

Let's start from importing the necessary DeepRank modules,

>>> from deeprank.learn import Dataset, NeuralNet

The ``Dataset`` module is used for preparing training/validation/test datasets and the ``NeuralNet`` module for training and testing neural networks.

Creating training/validation/test datasets
------------------------------------------

First we need to provide the HDF5 file we generated in the `Data Generation`_ step,

.. _Data Generation: tutorial2_dataGeneration.html

>>> database = '1ak4.hdf5'

You can also provide multiple HDF5 files to a database, e.g.

>>> database = ['1ak4.hdf5','native.hdf5']

Then we will use the data from this ``database`` to create the training/validation/test  datasets,

>>> data_set = DataSet(train_database=database,
>>>                   valid_database=None,
>>>                   test_database=None,
>>>                   chain1='C',
>>>                   chain2='D',
>>>                   grid_shape=(30, 30, 30),
>>>                   select_feature={
>>>                     'AtomicDensities_ind': 'all',
>>>                     'Feature_ind': ['coulomb', 'vdwaals', 'charge','PSSM_*']},
>>>                    select_target='BIN_CLASS',
>>>                    normalize_features=True,
>>>                    normalize_targets=False,
>>>                    pair_chain_feature=np.add,
>>>                    dict_filter={'IRMSD':'<4. or >10.'})

Here we have only one database, so we must provide it to the ``train_database`` parameter, and later you will set how to split this database to training, validation and test data. When independent validation and test databases exist, you can then set the corresponding ``valid_database`` and ``test_database`` parameters.

You must also specify which features and targets to use for the neural network by setting the ``select_feature`` and ``select_target`` parameters. By default, all features exist in the HDF5 files will be selected.

When you specify some of the features, you have to use a *dictionary*. As you can see in the HDF5 file, the ``mapped_features`` subgroup of each complex contains two groups ``AtomicDensities_ind`` and ``Feature_ind``. We'll forget about the ``_ind`` that is for legacy reasons, but these two groups contains the atomic densities and other features respectively. Therefore the value used in the example above:

>>> select_feature={'AtomicDensities': 'all',
>>>                 'Features': ['coulomb', 'vdwaals', 'charge', 'PSSM_*']}

It means that we want to use all the atomic densities features, the ``coulomb``, ``vdwaals`` and ``charge`` features as well as the PSSM features with a name starting with ``PSSM_``.

With the ``normalize_features`` and ``normalize_targets`` parameters, you can use the normalized values of features or targets to train a neural network.

As you can see in the HDF5 file for each feature, the data of ``chain1`` and ``chain2`` are stored separately. But it is possible to combine the feature of both chains into a single channel via the parameter ``pair_chain_feature``, e.g.,

>>> pair_chain_feature = np.add

which means that feature value of ``chain1`` and ``chain2`` will be summed up to create a single channel.

You can further filter the data with specific conditions,

>>> dict_filter={'IRMSD':'<4. or >10.'}

it will only selects the complexes whose IRMSD are lower than 4Å or larger than 10Å.


Creating neural network architecture
------------------------------------

The architecture of the neural network has to be well defined before training and DeepRank provides a very useful tool ``modelGenerator`` to facilitate this process.

Here is a simple example showing how to generate NN architecture:

>>> from deeprank.learn.modelGenerator import *
>>>
>>> conv_layers = []
>>> conv_layers.append(conv(output_size=4,kernel_size=2,post='relu'))
>>> conv_layers.append(pool(kernel_size=2))
>>> conv_layers.append(conv(input_size=4,output_size=5,kernel_size=2,post='relu'))
>>> conv_layers.append(pool(kernel_size=2))
>>>
>>> fc_layers = []
>>> fc_layers.append(fc(output_size=84,post='relu'))
>>> fc_layers.append(fc(input_size=84,output_size=2))
>>>
>>> gen = NetworkGenerator(name='cnn_demo',
>>>                       fname='cnn_arch_demo.py',
>>>                       conv_layers=conv_layers,
>>>                       fc_layers=fc_layers)
>>> gen.print()
>>> gen.write()

It will print out the human readable summary of the architecture,

.. code-block:: python

    #----------------------------------------------------------------------
    # Network Structure
    #----------------------------------------------------------------------
    #conv layer   0: conv | input -1  output  4  kernel  2  post relu
    #conv layer   1: pool | kernel  2  post None
    #conv layer   2: conv | input  4  output  5  kernel  2  post relu
    #conv layer   3: pool | kernel  2  post None
    #fc   layer   0: fc   | input -1  output  84  post relu
    #fc   layer   1: fc   | input  84  output  2  post None
    #----------------------------------------------------------------------

and also generate a ``cnn_arch_demo.py`` file with a class ``cnn_demo`` that defines the NN architecture.

We also provide the predefined NN architectures in ``learn/model3d.py`` and ``learn/model2d.py``.

Training a neural network
-------------------------

We are now all set to start the deep learning experiments. DeepRank supports both classification and regression tasks.


Classification with 3D CNN
^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> from deeprank.learn.model3d import cnn_class
>>>
>>> model = NeuralNet(data_set=data_set,
>>>                  model=cnn_class,
>>>                  model_type='3d',
>>>                  task='class')

``data_set`` is the dataset created above and ``cnn_class`` is the predefined NN architecture.
We also need to specify the ``model_type`` and the learning ``task``.

Then we can start the training process,

>>> model.train(nepoch=50,
>>>             divide_trainset=[0.7, 0.2, 0.1]
>>>             train_batch_size=5,
>>>             num_workers=1,
>>>             hdf5='epoch_data_class.hdf5')

We specify here the number of epoch, the fraction of data for training/validation/test sets, the batch size and the number of workers (CPU threads) in charge of batch preparation, and the output HDF5 file for training results. The model will be saved to ``.pth.tar`` files, e.g. ``model_epoch_0001.pth.tar``.

Regression with 3D CNN
^^^^^^^^^^^^^^^^^^^^^^

To train a regression model, the steps are same as the classification above. But you need to provide the regression NN architecture and set the correct task type, e.g.

>>> from deeprank.learn.model3d import cnn_reg
>>>
>>> model = NeuralNet(data_set=data_set,
>>>                  model=cnn_reg,
>>>                  model_type='3d',
>>>                  task='reg')
>>>
>>> model.train(nepoch=50,
>>>             divide_trainset=[0.7, 0.2, 0.1]
>>>             train_batch_size=5,
>>>             num_workers=1,
>>>             hdf5='epoch_data_reg.hdf5')

2D CNN
^^^^^^

DeepRank allows to transform the 3D volumetric data to 2D data by slicing planes of the data and using each plane as given channel. Very little modification of the code are necessary to do so. The creation of the dataset is identical to the 3D case, and you must specify ``model_type=2d`` for ``NeuralNet``,

>>> from deeprank.learn.model2d import cnn
>>>
>>> model = NeuralNet(data_set=data_set_2d,
>>>                  model=cnn,
>>>                  model_type='2d',
>>>                  task='reg',
>>>                  proj2d=0)

The ``proj2d`` parameter defines how to slice the 3D volumetric data. Value of: 0, 1, 2 are possible to slice along the YZ, XZ or XY plane respectively.

Testing a neural network
------------------------

In many cases after you've trained the NN model, you would like to use the model to do prediction or to test the model's performance on new data. DeepRank provide a very easy way to do that, let's say we have got the trained classification model ``model.pth.tar``,

>>> from deeprank.learn.model3d import cnn_class
>>>
>>> database = '1AK4_test.hdf5'
>>> model = NeuralNet(database, cnn_class, pretrained_model='model.pth.tar')
>>> model.test()

Note that here the database is simply the name of the hdf5 file we want to test the model on. All the processing of the dataset will be automatically done in the exact same way as it was done during the training of the model. Hence you do not have to copy the ``select_features`` and ``select_target`` parameters, all that is done for you automatically.
Learning
========

.. automodule:: deeprank.learn

This module contains all the tools for deep learning in DeepRank. The two main modules are ``deeprank.learn.DataSet`` and ``deeprank.learn.NeuralNet``. The ``DataSet`` class allows to process several hdf5 files created by the ``deeprank.generate`` toolset for use by pyTorch. This is done by creating several ``torch.data_utils.DataLoader`` for the training, valiation and test of the model. Several options are possible to specify and filter which conformations should be used in the dataset. The ``NeuralNet`` class is in charge of the deep learning part.There as well several options are possible to specify the task to be performed, the architecture of the neural network etc ....

Example:

>>> from deeprank.learn import *
>>> from model3d import cnn_class
>>>
>>> database = '1ak4.hdf5'
>>>
>>> # declare the dataset instance
>>> data_set = DataSet(database,
>>>                    chain1='C',
>>>                    chain2='D',
>>>                    select_feature='all',
>>>                    select_target='IRMSD',
>>>                    dict_filter={'IRMSD':'<4. or >10.'})
>>>
>>>
>>> # create the network
>>> model = NeuralNet(data_set, cnn_class, model_type='3d', task='class')
>>>
>>> # start the training
>>> model.train(nepoch = 250,divide_trainset=0.8, train_batch_size = 50, num_workers=8)
>>>
>>> # save the model
>>> model.save_model()

The details of the submodules are presented here. The two main ones are ``deeprank.learn.DataSet`` and ``deeprank.learn.NeuralNet``.

:note: The module ``deeprank.learn.modelGenerator`` can automatically create the file defining the neural network architecture.

DataSet: create a torch dataset
-------------------------------

.. automodule:: deeprank.learn.DataSet
    :members:
    :undoc-members:
    :private-members:


NeuralNet: perform deep learning
--------------------------------

.. automodule:: deeprank.learn.NeuralNet
    :members:
    :undoc-members:
    :private-members:

modelGenerator: generate NN architecture
----------------------------------------

.. automodule:: deeprank.learn.modelGenerator
    :members:
    :undoc-members:
    :private-members:
Targets
=======


.. automodule:: deeprank.targets


This module contains all the functions to compute target values for molecular structures. The implemented targets at the moment are:
    - ``binary_class``: Binary class ID
    - ``capri_class``: CAPRI class
    - ``dockQ``: DockQ
    - ``rmsd_fnat``: CAPRI metric IRMSD, LRMSD or FNAT

As you can see in the source each python file contained a ``__compute_target__`` function. This is the function called in ``deeprank.generate``.

Here are detailed the function in charge of target calculations.

Binary class
------------

.. autofunction:: deeprank.targets.binary_class.__compute_target__

CAPRI class
-----------

.. autofunction:: deeprank.targets.capri_class.__compute_target__

DockQ
-----

.. autofunction:: deeprank.targets.dockQ.__compute_target__


RMSDs & FNAT
------------

.. autofunction:: deeprank.targets.rmsd_fnat.__compute_target__

Data Generation
===============

This section describes how to generate feature and target data using PDB files
and/or other relevant raw data, which will be directly fed to the `learning`_ step.

.. _learning: tutorial3_learning.html

The generated data is stored in a single HDF5 file. The use of HDF5 format for
data storage not only allows to save disk space through compression but also reduce I/O time during the deep learning step.

The following tutorial uses the example from the test file ``test/test_generate.py``.
All the required data, e.g. PDB files and PSSM files, can be found in the ``test/1AK4/`` directory. It's assumed that the working directory is ``test/``.


Let's start:

First import the necessary modules,

>>> from mpi4py import MPI
>>> from deeprank.generate import *

The second line imports all the submodules required for data generation.

Then, we need to set the MPI communicator, that will be used later

>>> comm = MPI.COMM_WORLD

Calculating features and targets
--------------------------------

The data generation requires the PDB files of docking decoys for which we want to compute features and targets. We use ``pdb_source`` to specify where these decoy files are located,

>>> pdb_source = ['./1AK4/decoys']

which contains 5 docking decoys of the 1AK4 complex. The structure information of
these 5 decoys will be copied to the output HDF5 file.

We also need to specify the PDB files of native structures, which are required to
calculate targets like RMSD, FNAT, etc.

>>> pdb_native = ['./1AK4/native']

DeepRank will automatically look for native structure for each docking decoy by
name matching. For example, the native structure for the decoy ``./1AK4/decoys/1AK4_100w.pdb`` will be ``./1AK4/native/1AK4.pdb``.

Then, if you want to compute PSSM-related features like ``PSSM_IC``, you must specify the path to the PSSM files. The PSSM file must be named as ``<PDB_ID>.<Chain_ID>.pssm`` .
To produce PSSM files that are coherent with your pdb files and have already the correct file names, you can check out our module `PSSMGen <https://github.com/DeepRank/PSSMGen>`_.

>>> pssm_source = ['./1AK4/pssm_new/']

Finally we must specify the name of the output HDF5 file,

>>> h5out = '1ak4.hdf5'

We are now ready to initialize the ``DataGenerator`` class,

>>> database = DataGenerator(pdb_source=pdb_source,
                            pdb_native=pdb_native,
                            pssm_source=pssm_source,
                            chain1='C',
                            chain2='D',
                            compute_features=['deeprank.features.AtomicFeature',
                                              'deeprank.features.FullPSSM',
                                              'deeprank.features.PSSM_IC',
                                              'deeprank.features.BSA'],
                            compute_targets=['deeprank.targets.dockQ'],
                            hdf5=h5out)

The ``compute_features`` and ``compute_targets`` are used to set which features
and targets we want to calculate by providing a list of python class names. The
available predefined features and targets can be found in the `API Reference`_.

.. _API Reference: Documentation.html

The above statement initializes a class instance but does not compute anything yet.
To actually execute the calculation, we must use the ``create_database()`` method,

>>> database.create_database(prog_bar=True)

Once you punch that DeepRank will fo through all the protein complexes specified
as input and compute all the features and targets required.

Mapping features to 3D grid
---------------------------
The next step consists in mapping the features calculated above to a grid of points centered around the molecule interface. Before mapping the features, we must define the grid first,

>>> grid_info = {'number_of_points': [30, 30, 30],
                 'resolution': [1., 1., 1.],
                 'atomic_densities': {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8}
                }

Here we define a grid with 30 points in x/y/z direction and a distance interval of 1Å between two neighbor points. We also provide the atomic radius for the element `C`, `N`, `O` and `S` to calculate atomic densities. Atomic density is also a type of features which have to be calculated during the mapping.

Then we start mapping the features to 3D grid,

>>> database.map_features(grid_info, try_sparse=True, prog_bar=True)

By setting ``try_sparse`` to ``True``, DeepRank will try to store the 3D grids with a built-in sparse format, which can seriously save the storage space.

Finally, it should generate the HDF5 file ``1ak4.hdf5`` that contains all the raw features and targets data and the grid-mapped data which will be used for the learning step. To easily explore these data, you could try the `DeepXplorer`_ tool.

.. _DeepXplorer: https://github.com/DeepRank/DeepXplorer

Appending features/targets
--------------------------

Suppose you've finished generating a huge HDF5 database and just realize you forgot to compute some specific features or targets. Do you have to recompute everything?

No! You can append more features and targets to the existing HDF5 database in a very simple way:

>>> h5file = '1ak4.hdf5'
>>> database = DataGenerator(compute_targets=['deeprank.targets.binary_class'],
>>>                          compute_features=['deeprank.features.ResidueDensity'],
>>>                          hdf5=h5file)
>>>
>>> # add targets
>>> database.add_target()
>>>
>>> # add features
>>> database.add_feature()
>>>
>>> # map features
>>> database.map_features()

Voilà! Here we simply specify the name of the existing HDF5 file we generated above, and set the new features/targets to add to this database. The methods ``add_target`` and ``add_feature`` are then called to calculate the corresponding targets and features. Don't forget to map the new features afterwards. Note that you don't have to provide any grid information for the mapping, because DeepRank will automatically detect and use the grid info that exist in the HDF5 file.
Advanced Tutorial
=================

DeepRank allows users to define and create new features and targets, and this section will show how to that.

Create your own features
------------------------

To create your own feature you simply have to create a feature class that must subclass the ``FeatureClass`` contained in ``features/FeatureClass.py``. As an example we will create here a new feature that maps the carbon alpha of the contact residue. The first thing we need to do is to import ``pdb2sql`` and the ``FeatureClass`` superclass,

>>> import pdb2sql
>>> from deeprank.features import FeatureClass

We then have to define the class and its initialization. Here we will simply initialize the class with the pdb information of the molecule we want to process. This is therefore given by

>>> # a new class based on the FeatureClass
>>> class CarbonAlphaFeature(FeatureClass):
>>>
>>>     # init the class
>>>     def __init__(self,pdbfile):
>>>         super().__init__('Atomic')
>>>         self.pdb = pdb
>>>
>>>     # the feature extractor
>>>     def get_feature(self):
>>>
>>>         # create a sql database
>>>         db = pdb2sql(self.pdb)
>>>
>>>         # get the contact atoms
>>>         indA,indB = db.get_contact_atoms()
>>>         contact = indA + indB
>>>
>>>         # extract the atom keys and xyz of the contact CA atoms
>>>         ca_keys = db.get('chainID,resName,resSeq,name',name='CA',rowID=contact)
>>>         ca_xyz = db.get('x,y,z',name='CA',rowID=contact)
>>>
>>>         # create the dictionary of human readable and xyz-val data
>>>         hread, xyzval = {},{}
>>>         for key,xyz in zip(ca_keys,ca_xyz):
>>>
>>>             # human readable
>>>             # { (chainID,resName,resSeq,name): [val] }
>>>             hread[tuple(key)] = [1.0]
>>>
>>>             # xyz-val
>>>             # { (0|1,x,y,z): [val] }
>>>             chain = [{'A':0,'B':1}[key[0]]]
>>>             k = tuple( chain + xyz)
>>>             xyzval[k] = [1.0]
>>>
>>>         self.feature_data['CA'] = hread
>>>         self.feature_data_xyz['CA'] = xyzval

As you can see we must initialize the superclass. Since we use here the argument ``Atomic``, the feature is an atomic based feature. If we would create a residue based feature (e.g. PSSM, residue contact, ...) we would have used here the argument ``Residue``. This argument specifies the printing format of the feature in a human readable format and doesn't affect the mapping. From the super class the new class inherit two methods

>>> export_data_hdf5(self,featgrp)
>>> export_dataxyz_hdf5(self,featgrp)

which are used to store feature values in the HDF5 file.

The class also inherits two variables

>>> self.feature_data = {}
>>> self.feature_data_xyz = {}

where we must store the feature in human readable format and xyz-val format.

To extract the feature value we are going to write a method in the class in charge of the feature extraction. This method is simply going to locate the carbon alpha and gives a value of 1 at the corresponding xyz positions.

In this function we exploit ``pdb2sql`` to locate the carbon alpha that make a contact. We then create two dictionaries where we store the feature value. The format of the human readable and xyz-val are given in comment. These two dictionaries are then added to the superclass variable ``feature_data`` and ``feature_data_xyz``.

We now must use this new class that we've just created in DeepRank. To do that we must create a function called:

>>> def __compute_feature__(pdb_data,featgrp,featgrp_raw)

Several example can be found in the feature already included in DeepRank. The location of the this function doesn't matter as we will provide the python file as value of the ``compute_features`` argument in the ``DataGenerator``. For the feature we have just created we can define that function as

>>> def __compute_feature__(pdb_data,featgrp,featgrp_raw):
>>>
>>>        cafeat = CarbonAlphaFeature(pdb_data)
>>>        cafeat.get_features()
>>>
>>>        # export in the hdf5 file
>>>        cafeat.export_dataxyz_hdf5(featgrp)
>>>        cafeat.export_data_hdf5(featgrp_raw)
>>>
>>>        # close
>>>        cafeat.db._close()

Finally to compute this feature we must call it during the data generation process. Let's assume that the file containing the ``__compute_feature__`` function is in the local folder and is called ``CAfeature.py``. To use this new feature in the generation we can simply pass the name of this file in the DataGenerator as

>>> database = DataGenerator(pdb_source=pdb_source, pdb_native=pdb_native,
>>>                          compute_features = ['CAFeature', ....], ...)


Create your own targets
-----------------------

The creation of new targets is similar to the process of creating new features but simpler. The targets don't need to be mapped on a grid and therefore don't need any fancy formatting. We simply need to create a new dataset in the target group of the molecule concerned. For example let's say we want to associate a random number to each conformation. To do that we can use the following code:

>>> import numpy as np
>>>
>>> def get_random_number():
>>>     return np.random.rand()
>>>
>>> def __compute_target__(pdb_data, targrp):
>>>
>>>     target = get_random_number()
>>>     targrp.create_dataset('FOO', data=np.array(target))

As for the features, the new target must be called in a function with a very precise name convention:

>>> def __compute_target__(pdb_data, targrp)

If as before we assume that the file containing this function is in the local folder and is called ``random.py`` we can compute the target by calling the ``DataGenerator`` with:

>>> database = DataGenerator(pdb_source=pdb_source, pdb_native=pdb_native,
>>>                          compute_targets = ['random',....], ...)
Features
========

This module contains all the tools to compute feature values for molecular structure. Each submodule must be subclass ``deeprank.features.FeatureClass`` to inherit the export function. At the moment a few features have already been implemented. These are:
    - ``AtomicFeatures``:Coulomb, van der Waals interactions and atomic charges
    - ``BSA`` : Burried Surface area
    - ``FullPSSM`` : Complete PSSM data
    - ``PSSM_IC`` : Information content of the PSSM
    - ``ResidueDensity`` : The residue density for polar/apolar/charged pairs

As you can see in the source each python file contained a ``__compute_feature__`` function. This is the function called in ``deeprank.generate``.



Here are detailed the class in charge of feature calculations.


Atomic Feature
--------------

.. automodule:: deeprank.features.AtomicFeature
    :members:
    :undoc-members:
    :private-members:

.. autofunction:: deeprank.features.AtomicFeature.__compute_feature__

Buried Surface Area
-------------------

.. automodule:: deeprank.features.BSA
    :members:
    :undoc-members:
    :private-members:

.. autofunction:: deeprank.features.BSA.__compute_feature__

FullPSSM
--------

.. automodule:: deeprank.features.FullPSSM
    :members:
    :undoc-members:
    :private-members:

.. autofunction:: deeprank.features.FullPSSM.__compute_feature__


PSSM Information Content
------------------------

.. automodule:: deeprank.features.PSSM_IC
    :members:
    :undoc-members:
    :private-members:

.. autofunction:: deeprank.features.PSSM_IC.__compute_feature__


Contact Residue Density
-----------------------

.. automodule:: deeprank.features.ResidueDensity
    :members:
    :undoc-members:
    :private-members:

.. autofunction:: deeprank.features.ResidueDensity.__compute_feature__

Generic Feature Class
---------------------

.. automodule:: deeprank.features.FeatureClass
    :members:
    :undoc-members:
    :private-members:
Tools
=======================

.. automodule:: deeprank.tools

This module contains a series of core tools used for the feature calculations in DeepRank. These tools are at the moment:

    - sasa: a simple solvent surface area calculator
    - sparse: a 3D sparse matrix engine



Here are the details of the submodule given in tools.

Solvent Accessible Surface Area
---------------------------------

.. automodule:: deeprank.tools.sasa
    :members:
    :undoc-members:
    :private-members:


Sparse 3D Matrix
------------------------------

.. automodule:: deeprank.tools.sparse
    :members:
    :undoc-members:
    :private-members:
.. highlight:: rst

Installation
============

Installing via pypi
-------------------
To install DeepRank package via pypi_:

.. code-block:: bash

  pip install deeprank

.. _pypi: https://pypi.org/project/deeprank

Installing from source code
----------------------------

To install DeepRank using GitHub_ source code:

.. _GitHub: https://github.com/DeepRank/deeprank

.. code-block:: bash

  # clone the repository
  git clone https://github.com/DeepRank/deeprank.git

  # enter the deeprank directory
  cd deeprank

  # install the module
  pip install -e ./Tutorial
========

.. toctree::
   :maxdepth: 3

   tutorial1_installing
   tutorial2_dataGeneration
   tutorial3_learning
   tutorial4_advanced.. DeepRank documentation master file, created by
   sphinx-quickstart on Mon Feb 26 14:44:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DeepRank
========

`DeepRank`_ is a general, configurable deep learning framework for data mining
protein-protein interactions (PPIs) using 3D convolutional neural networks (CNNs).

DeepRank contains useful APIs for pre-processing PPIs data, computing features and
targets, as well as training and testing CNN models.

.. _`DeepRank`: https://github.com/DeepRank/deeprank

**DeepRank highlights**:

- Predefined atom-level and residue-level PPI feature types
   - *e.g. atomic density, vdw energy, residue contacts, PSSM, etc.*
- Predefined target types
   - *e.g. binary class, CAPRI categories, DockQ, RMSD, FNAT, etc.*
- Flexible definition of both new features and targets
- 3D grid feature mapping
- Efficient data storage in HDF5 format
- Support both classification and regression (based on PyTorch)

Tutorial
--------

.. toctree::
   :maxdepth: 3

   tutorial

API Reference
-------------

.. toctree::
   :maxdepth: 3

   Documentation

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
Data Generation
===============

.. automodule:: deeprank.generate

This module contains all the tools to compute the features and targets and to map the features onto a grid of points. The main class used for the data generation is ``deeprank.generate.DataGenerator``. Through this class you can specify the molecules you want to consider, the features and the targets that need to be computed and the way to map the features on the grid. The data is stored in a single HDF5 file. In this file, each conformation has its own group that contains all the information related to the conformation. This includes the pdb data, the value of the feature (in human readable format and xyz-val format), the value of the targe values, the grid points and the mapped features on the grid.


At the moment a number of features are already implemented. This include:

    - Atomic densities
    - Coulomb & vd Waals interactions
    - Atomic charges
    - PSSM data
    - Information content
    - Buried surface area
    - Contact Residue Densities

More features can be easily implemented and integrated in the data generation workflow. You can see example here. The calculation of a number of target values have also been implemented:

    - i-RMSD
    - l-RMSD
    - FNAT
    - DockQ
    - binary class

There as well new targets can be implemented and integrated to the workflow.

Normalization of the data can be time consuming as the dataset becomes large. As an attempt to alleviate this problem, the class ``deeprank.generate.NormalizeData`` has been created. This class directly compute and store the standard deviation and mean value of each feature within a given hdf5 file.

Example:

>>> from deeprank.generate import *
>>> from time import time
>>>
>>> pdb_source     = ['./1AK4/decoys/']
>>> pdb_native     = ['./1AK4/native/']
>>> pssm_source    = './1AK4/pssm_new/'
>>> h5file = '1ak4.hdf5'
>>>
>>> #init the data assembler
>>> database = DataGenerator(chain1='C',chain2='D',
>>>                          pdb_source=pdb_source,pdb_native=pdb_native,pssm_source=pssm_source,
>>>                          data_augmentation=None,
>>>                          compute_targets  = ['deeprank.targets.dockQ'],
>>>                          compute_features = ['deeprank.features.AtomicFeature',
>>>                                              'deeprank.features.FullPSSM',
>>>                                              'deeprank.features.PSSM_IC',
>>>                                              'deeprank.features.BSA'],
>>>                          hdf5=h5file)
>>>
>>> t0 = time()
>>> #create new files
>>> database.create_database(prog_bar=True)
>>>
>>> # map the features
>>> grid_info = {
>>>     'number_of_points': [30,30,30],
>>>     'resolution': [1.,1.,1.],
>>>     'atomic_densities': {'CA':3.5,'N':3.5,'O':3.5,'C':3.5},
>>> }
>>> database.map_features(grid_info,try_sparse=True,time=False,prog_bar=True)
>>>
>>> # add a new target
>>> database.add_target(prog_bar=True)
>>> print(' '*25 + '--> Done in %f s.' %(time()-t0))
>>>
>>> # get the normalization
>>> norm = NormalizeData(h5file)
>>> norm.get()

The details of the different submodule are listed here. The only module that really needs to be used is ``DataGenerator`` and ``NormalizeData``. The ``GridTools`` class should not be directly used by inexperienced users.


Structure Alignement
----------------------------------------

All the complexes contained in the dataset can be aligned similarly to facilitate and improve the training of the model. This can easily be done using the `align` option of the `DataGenerator` for example to align all the complexes along the 'z' direction one can use:

>>> database = DataGenerator(chain1='C',chain2='D',
>>>                          pdb_source=pdb_source, pdb_native=pdb_native, pssm_source=pssm_source,
>>>                          align={"axis":'z'}, data_augmentation=2,
>>>                          compute_targets=[ ... ], compute_features=[ ... ], ... )


Other options are possbile, for example if you would like to have the alignement done only using a subpart of the complex, say the chains A and B you can use:

>>> database = DataGenerator(chain1='C',chain2='D',
>>>                          pdb_source=pdb_source, pdb_native=pdb_native, pssm_source=pssm_source,
>>>                          align={"axis":'z', "selection": {"chainID":["A","B"]} }, data_augmentation=2,
>>>                          compute_targets=[ ... ], compute_features=[ ... ], ... )

All the selection offered by `pdb2sql` can be used in the `align` dictionnary e.g.: "resId":[1,2,3], "resName":['VAL','LEU'], ... Only the atoms selected will be aligned in the give direction.

You can also try to align the interface between two chains in a given plane. This can be done using:

>>> database = DataGenerator(chain1='C',chain2='D',
>>>                          pdb_source=pdb_source, pdb_native=pdb_native, pssm_source=pssm_source,
>>>                          align={"plane":'xy', "selection":"interface"}, data_augmentation=2,
>>>                          compute_targets=[ ... ], compute_features=[ ... ], ... )

which by default will use the interface between the first two chains. If you have more than two chains in the complex and want to specify which chains are forming the interface to be aligned you can use:

>>> database = DataGenerator(chain1='C',chain2='D',
>>>                          pdb_source=pdb_source, pdb_native=pdb_native, pssm_source=pssm_source,
>>>                          align={"plane":'xy', "selection":"interface", "chain1":'A', "chain2":'C'}, data_augmentation=2,
>>>                          compute_targets=[ ... ], compute_features=[ ... ], ... )

DataGenerator
----------------------------------------

.. automodule:: deeprank.generate.DataGenerator
    :members:
    :undoc-members:
    :show-inheritance:
    :private-members:

NormalizeData
----------------------------------------

.. automodule:: deeprank.generate.NormalizeData
    :members:
    :undoc-members:
    :private-members:

GridTools
------------------------------------

.. automodule:: deeprank.generate.GridTools
    :members:
    :undoc-members:
    :private-members:API Reference
=============

.. toctree::
   :maxdepth: 3

   deeprank.generate
   deeprank.learn
   deeprank.features
   deeprank.targets
   deeprank.tools
