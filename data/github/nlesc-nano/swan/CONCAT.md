# Change Log

# 0.7.0 [Unreleased]


# 0.6.0 [21/06/2021]
## New
* Add interface to scikit regressors (#85)
* Add interface to HDF5 to store the training results (#88)
* Add interface to GPyTorch (#90)

## Changed
* Fix prediction functionality (#81)
* Return predicted and expected values scaled back (#90)

# 0.5.0 [04/05/2021]
## Changed
* Fix graph neural network implementation (#59)
* Rename the graph neural network to MPNN (message passing NN)

## New
* Interface to [se3-transformer](https://www.dgl.ai/pages/start.html) (#57)
* Calculate a guess to the molecular coordinates using a force field optimization
* Add bond distance as feature for the GNN using the optimized geometries (#59)
* Add early stopping functionality

# 0.4.0 [02/10/2020]

## Changed
* Removed duplicate CAT functionality (#36)
* Moved the properties computation and filtering to its own repo(#44, #45)

# 0.3.0 [21/08/2020]

## New
* Introduce Pipeline to filter ligands (#13, #26)
* Use [SCScore](https://pubs.acs.org/doi/10.1021/acs.jcim.7b00622)
* Use [Horovod](https://github.com/horovod/horovod) to distribute the training 
* Add [mypy](http://mypy-lang.org/) test

# 0.2.0 [25/02/2020]

## New
* Allow to  train in **GPU**.
* Use [Pytorch-geometric](https://github.com/rusty1s/pytorch_geometric) to create molecular graph convolutional networks.

## Changed
* Replace [deepchem](https://deepchem.io/) with [pytorch](https://pytorch.org)
* Use [rdkit](https://www.rdkit.org) to compute the fingerprints.

# 0.1.0

## New

* Interface to [deepchem](https://deepchem.io/) to generate regression models.
* Interface to [CAT](https://github.com/nlesc-nano/CAT) to compute solvation energies and activity coefficients.
* Train or predict values for a set of smiles using a statistical model.
---
name: Bug report
about: Something doesn't work like I expected.
title: ''
labels: bug
assignees: ''

---

Use your best judgment to provide a useful level of information. Depending on the nature of the issue, consider including, e.g.

- Which version of the software you're running
- Any logging or error output you see
---
name: Help wanted
about: I need some help with the code
title: ''
labels: help wanted
assignees: ''

---
---
name: Feature request
about: I would like a new feature to be included in the library
title: ''
labels: enhancement
assignees: ''

---

Tell us what would be a nice feature to make the **swan** library even better!
### All Submissions:

* [ ] Have you followed the guidelines in our Contributing document?
* [ ] Have you checked to ensure there aren't other open [Pull Requests](https://github.com/SCM-NV/qmflows-namd/pulls) for the same update/change?


### New Feature Submissions:

1. [ ] Does your submission pass tests?
2. [ ] Have you check your code quality (e.g. using [flake8](http://flake8.pycqa.org/en/latest/)) prior to submission?

### Changes to Core Features:

* [ ] Have you added an explanation of what your changes do and why you'd like us to include them?
* [ ] Have you written new tests for your core changes, as applicable?
* [ ] Have you successfully ran tests with your changes locally?
# Machine learning model to predict CDFT properties

The following table shows the **mean rvalue** resulting from the linear regression
between the predicted properties as function of the ground true computed properties. Therefore, the closer
these values are to 1 the better the model is to predict the properties *with
the available amount of training data*.

The mean values where computing with the results of **5** training/validation datasets picked randomly
in a 0.8/0.2 proportion, respectively.

The hyperparameters use to perform the training can be found at [model's hyperparameters](Training_hyperparameters.md).


## Amines
The [amines CDFT dataset](https://github.com/nlesc-nano/swan/blob/main/data/Amines/CDFT/all_amines.csv) has the following properties:
* Number of molecules in dataset: **3109**
* Number of labels: **15**

The [geometries file](https://github.com/nlesc-nano/swan/blob/main/data/Amines/CDFT/all_amines.csv) contains all the optimized geometries for the previous dataset.


|          Property name             | FingerprintFullyConnected | MPNN | SE3Transformer|
|:----------------------------------:|:-------------------------:|:----:|:-------------:|
| Dissocation energy (nucleofuge)    | 0.73   | 0.70   | 0.78 |
| Dissociation energy (electrofuge)  | 0.39	  | 0.32   | 0.39 |
| Electroaccepting power(w+)         | 0.42	  | 0.36   | 0.55 |
| Electrodonating power (w-)         | 0.52	  | 0.54   | 0.66 |
| Electronegativity (chi=-mu)        | 0.54	  | 0.50   | 0.61 |
| Electronic chemical potential (mu) | 0.54	  | 0.53   | 0.59 |
| Electronic chemical potential (mu+)| 0.38	  | 0.25   | 0.34 |
| Electronic chemical potential (mu-)| 0.78	  | 0.81   | 0.85 |
| Electrophilicity index (w=omega)   | 0.51	  | 0.53   | 0.62 |
| Global Dual Descriptor Deltaf+     | 0.33	  | 0.26   | 0.35 |
| Global Dual Descriptor Deltaf-     | 0.35	  | 0.26   | 0.36 |
| Hardness (eta)                     | 0.53	  | 0.53   | 0.56 |
| Hyperhardness (gamma)              | 0.39	  | 0.29   | 0.41 |
| Net Electrophilicity               | 0.47	  | 0.55   | 0.62 |
| Softness (S)                       | 0.42	  | 0.07   | 0.20 |



## Carboxylic acids
The [carboxylic acids CDFT dataset](https://github.com/nlesc-nano/swan/blob/main/data/Carboxylic_acids/CDFT/all_carboxylics.csv) has the following properties:
* Number of molecules in dataset: **11044**
* Number of labels: **15**

The [geometries file](https://github.com/nlesc-nano/swan/blob/main/data/Carboxylic_acids/CDFT/all_geometries_carboxylics.json) contains all the optimized geometries for the previous dataset.


|          Property name             | FingerprintFullyConnected | MPNN | SE3Transformer|
|:----------------------------------:|:-------------------------:|:----:|:-------------:|
| Dissocation energy (nucleofuge)    | 0.85   | 0.85 |  0.87 |
| Dissociation energy (electrofuge)  | 0.75   | 0.83 | 	0.83 |
| Electroaccepting power(w+)         | 0.63	  | 0.64 | 	0.66 |
| Electrodonating power (w-)         | 0.79	  | 0.79 | 	0.84 |
| Electronegativity (chi=-mu)        | 0.81	  | 0.85 | 	0.88 |
| Electronic chemical potential (mu) | 0.81	  | 0.80 | 	0.88 |
| Electronic chemical potential (mu+)| 0.72	  | 0.78 | 	0.77 |
| Electronic chemical potential (mu-)| 0.86	  | 0.86 | 	0.91 |
| Electrophilicity index (w=omega)   | 0.72	  | 0.75 | 	0.79 |
| Global Dual Descriptor Deltaf+     | 0.59	  | 0.59 | 	0.61 |
| Global Dual Descriptor Deltaf-     | 0.59	  | 0.59 | 	0.60 |
| Hardness (eta)                     | 0.80	  | 0.74 | 	0.83 |
| Hyperhardness (gamma)              | 0.67	  | 0.68 | 	0.70 |
| Net Electrophilicity               | 0.72	  | 0.78 | 	0.80 |
| Softness (S)                       | 0.41	  | 0.57 | 	0.44 |

# Hyperparameters used to train the NN

The hyperparameters are taken from the following references:

* [Simple Feedforward NN](https://arxiv.org/abs/1509.09292) (No convolution is performed)
* [Neural Message Passing for Quantum Chemistry](https://arxiv.org/abs/1704.01212)
* [SE(3)-Transformer](https://arxiv.org/abs/2006.10503)

## FullyConnectedFingerPrint
| Parameter       | value |
| :-------------: | :---: |
| hidden_cells    | 100   |
| batch_size      | 100   |
| learning rate   | 5e-4  |


## MPNN

|     Parameter   | value |
| :-------------: | :---: |
| num_labels      | 1     |
| num_iterations  | 3     |
| output_channels | 10    |
| batch_size      | 20    |
| learning rate   | 5e-4  |


## SE3Transformer

| Parameter       | value  |
| :-------------: | :----: |
| num_layers      | 4      |
| num_channels    | 16     |
| num_nlayers     | 0      |
| num_degrees     | 4      |
| div             | 4      |
| pooling         | 'avg'  |
| n_heads         | 1      |
| batch_size      | 32     |
| learning rate   | 1e-3   |
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

#. use the search functionality `here <https://github.com/nlesc-nano/swan/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/nlesc-nano/swan/issues>`__ to see if someone already filed the same issue;
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
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the Screening Workflows And Nanomaterials repository on GitHub;
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
reported by contacting the project team at f.zapata@esciencecenter.nl. All
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

.. image:: https://github.com/nlesc-nano/swan/workflows/build%20with%20conda/badge.svg
   :target: https://github.com/nlesc-nano/swan/actions
.. image:: https://codecov.io/gh/nlesc-nano/swan/branch/main/graph/badge.svg?token=1527ficjjx
   :target: https://codecov.io/gh/nlesc-nano/swan
.. image:: https://zenodo.org/badge/191957101.svg
   :target: https://zenodo.org/badge/latestdoi/191957101
.. image:: https://readthedocs.org/projects/swan/badge/?version=latest
   :target: https://swan.readthedocs.io/en/latest/?badge=latest
	    
#####################################
Screening Workflows And Nanomaterials
#####################################

🦢 **Swan** is a Python pacakge to create statistical models using machine learning to predict molecular properties. See Documentation_.


🛠 Installation
===============

- Download miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_.

- Create a new virtual environment using the following commands:

  - ``conda create -n swan``

- Activate the new virtual environment

  - ``conda activate swan``

To exit the virtual environment type  ``conda deactivate``.


.. _dependecies:

Dependencies installation
-------------------------

- Type in your terminal:

  ``conda activate swan``

Using the conda environment the following packages should be installed:


- install RDKit_ and H5PY_:

  - `conda install -y -q -c conda-forge h5py rdkit`

- install Pytorch_ according to this_ recipe

- install `Pytorch_Geometric dependencies <https://github.com/rusty1s/pytorch_geometric#installation>`_.

- install `DGL using conda <https://www.dgl.ai/pages/start.html>`_


.. _installation:

Package installation
--------------------
Finally install the package:

- Install **swan** using pip:
  - ``pip install git+https://github.com/nlesc-nano/swan.git``

Now you are ready to use *swan*.


  **Notes:**

  - Once the libraries and the virtual environment are installed, you only need to type
    ``conda activate swan`` each time that you want to use the software.

.. _Documentation: https://swan.readthedocs.io/en/latest/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _anaconda: https://www.anaconda.com/distribution/#download-section
.. _installConda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _Pytorch: https://pytorch.org
.. _RDKit: https://www.rdkit.org
.. _H5PY: https://www.h5py.org/
.. _this: https://pytorch.org/get-started/locally/
Tutorial to Generate Statistical Models
=======================================
In this tutorial we explore how to create and train statistical models to predict
molecular properties using the Pytorch_ library. We will use smiles_ to represent the molecules
and use the csv_ file format to manipulate the molecules and their properties.

As an example, we will predict the `activity coefficient_` for a subset of carboxylic acids taken
from the `GDB-13 database_`. Firstly, We randomly takes a 1000 smiles_ from the database and
compute the `activity coefficient_` using the `COSMO approach_`. We store the values in the `thousand.csv_`
file.

A peek into the file will show you something like: ::

  smiles,E_solv,gammas
  OC(=O)C1OC(C#C)C2NC1C=C2,-11.05439751550119,8.816417146193844
  OC(=O)C1C2NC3C(=O)C2CC13O,-8.98188869016993,52.806217658944995
  OC(=O)C=C(C#C)C1NC1C1CN1,-11.386853547889574,6.413128231164093
  OC(=O)C1=CCCCC2CC2C#C1,-10.578966144649726,1.426566948888662

Where the first column contains the index of the row, the second the solvation energy and finally the
`activity coefficients_` denoted as *gammas*. Once we have the data we can start exploring different statistical methods.

`swan` offers a thin interface to Pytorch_. It takes yaml_ file as input and either train an statistical model or
generates a prediction using a previously trained model. Let's briefly explore the `swan` input.

Simulation input
****************
A typical `swan` input file looks like: ::

  dataset_file:
    tests/test_files/thousand.csv
  properties:
    - gammas

  use_cuda: True

  featurizer:
    fingerprint: atompair

  model:
    name: FingerprintFullyConnected
    parameters:
       input_features: 2048  # Fingerprint size
       hidden_cells: 200
       output_features: 1  # We are predicting a single property 

  torch_config:
    epochs: 100
    batch_size: 100
    optimizer:
      name: sgd
      lr: 0.002

   
**dataset_file**: A csv_ file with the smiles_ and other molecular properties.

**properties**: the columns names of hte csv_ file representing the molecular properties to fit.

**featurizer**: The type of transformation to apply to the smiles_ to generates the features_. Could be either **fingerprint** or **graph**.

Have a look at the  :ref:`available models`.

Training a model
****************
In order to run the training, run the following command: ::

  modeller --mode train -i input.yml

`swan` will generate a log file called  `output.log` with a timestamp for the different steps during the training.
Finally, you can see in your `cwd` a folder called *swan_models* containing the parameters of your statistical model.

It is possible to restart the training procedure by providing the ``--restart`` option like: ::

    modeller --mode train -i input.yml --restart

Predicting new data
*******************
To predict new data you need to provide some smiles for which you want to compute the properties of interest, in this
case the `activity coefficient_`. For doing so, you need to provide in the `dataset_file` entry of the *input.yml*
file the path to a csv_ file containing the smiles, like the `smiles.csv_`: ::

  ,smiles
  0,OC(=O)C1CNC2C3C4CC2C1N34
  1,OC(=O)C1CNC2COC1(C2)C#C
  2,OC(=O)CN1CC(=C)C(C=C)C1=N

Then run the command: ::

  modeler --mode predict -i input.yml

`swan` will look for a *swan_model.pt* file with the previously trained model and will load it.

Finally, you will find a file called "predicted.csv" with the predicted values for the activity coefficients.

..  _deepchem: https://deepchem.io/
.. _smiles: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
.. _activity coefficient: https://en.wikipedia.org/wiki/Activity_coefficient
.. _GDB-13 database_`: https://pubs.acs.org/doi/abs/10.1021/ja902302h
.. _COSMO approach: https://www.scm.com/doc/ADF/Input/COSMO.html
.. _thousand.csv: https://github.com/nlesc-nano/swan/blob/master/tests/test_files/thousand.csv
.. _features: https://en.wikipedia.org/wiki/Feature_(machine_learning)
.. _smiles.csv: https://github.com/nlesc-nano/swan/blob/master/tests/test_files/smiles.csv
.. _yaml: https://yaml.org
.. _csv: https://en.wikipedia.org/wiki/Comma-separated_values
.. _Pytorch: https://pytorch.org
.. _available models:

Available models
================
Currently **Swan** Implements the following models:

Fully Connected Neural Network
******************************
A standard fully connected neural network that takes *fingerprints* as
input features. To use the model you need to specify in the ``model`` section
of the input YAML file the following: ::

  model:
    name: FingerprintFullyConnected
    parameters:
      input_features: 2048
      hidden_cells: 100
      output_features: 1

The model takes 3 additional optional parameters:
* ``input_features``: fingerprint size. Default 2048.
* ``hidden_cells``: Hiden number of cell(or nodes). Default 100.
* ``num_labels``: the amount of labels to predict. Default 1.

Also, the model requires as a ``featurizer`` a fingerprint calculator that can be provided like: ::

  featurizer:
    fingerprint: atompair

Available fingerprints algorithms are: ``atompair`` (default), ``morgan`` or ``torsion``. These
algorithms are provided by `RDKIT descriptor package <https://rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html>`_.


Message Passing Neural Network
******************************
Implementation of the message passing neural network (MPNN) reported at `<https://arxiv.org/abs/1704.01212>`_.
If you don't have an idea what a MPNN is have a look at
`this introduction to Graph Neural Networks <https://www.youtube.com/watch?v=zCEYiCxrL_0&list=PLVqPBNulzDDg8ieQZ2G643UFbHm-qWW7Z&index=1&t=2239s>`_.

To train your model using the MPNN you need to provide the following section in the YAML input file: ::

  model:
    name: MPNN
    parameters:
      output_channels: 10
      num_labels: 1
      batch_size: 128
      num_iterations: 3

The optional parameters for the model are: ::
* ``output_channels`` Channels in the Convolution. default 10.
* ``num_labels``: the amount of labels to predict. Default 1.
*  ``batch_size``: the size of the batch used to train the model. Default 128.
* ``num_iterations``: number of steps to interchange messages for each epoch. Default 3.

Additionally the model requires the use of the following featurizer: ::

  featurizer:
   graph: molecular
   file_geometries: geometries.json

Where ``file_geometries`` is a JSON file containing an array of molecules on PDB format. Check
`the example file <https://github.com/nlesc-nano/swan/blob/main/tests/files/cdft_geometries.json>`_
If the ``file_geometries`` is not set in the input the model will try to use the RDKit geometries.
API Data Representation
=======================

.. currentmodule:: swan.dataset

Data Base Class
###############
.. autoclass:: swan.dataset.swan_data_base.SwanDataBase
    :members:

Graph Data Base Class
#####################
.. autoclass:: swan.dataset.data_graph_base.SwanGraphData
    :members:

Fingerprints Data
#################
.. autoclass:: FingerprintsData
    :members:

Torch Geometric Data
####################
.. autoclass:: TorchGeometricGraphData
    :members:

DGL Data
########
.. autoclass:: DGLGraphData
    :members:
.. Screening Workflows And Nanomaterials documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SWAN!
================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   includereadme
   tutorial_models
   available_models
   training_and_validation
   api_data
   api_models


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. include:: ../README.rst
API Statistical Models
======================
Available models

.. currentmodule:: swan.modeller.models

Deep Feedforward Network
########################


.. autoclass:: FingerprintFullyConnected
   :members:


Message Passing Graph Neural Network
####################################

.. autoclass:: MPNN
   :members:


Equivariant Neural Networks
###########################

.. autoclass:: TFN

.. autoclass:: SE3TransformerTraining and validation
=======================
The training and validation functionality is implemented by the *Modeller* class.

.. currentmodule:: swan.modeller.modeller

.. autoclass:: swan.modeller.Modeller
   :members:

