# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

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

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at j.m.hoch@uu.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# How to contribute

This python-package is a first outcome of an interdisciplinary project aimed at understanding the complex interplay between conflict and climate and environment.
As such, the presented code and functionalities can only be seen as a first step towards a fully-fledged model.
We therefore strongly encourage other users to contribute to this project!

## General notes

When contributing to this repository, please first discuss the change you wish to make via issue, email, or any other method with the owners of this repository before making a change.

Please note we have a code of conduct, please follow it in all your interactions with the project, the project owners, and users of the project.

## Getting Started

* Make sure you have a GitHub account.
* Fork the repository on GitHub.

## Making Changes

* Create a topic branch from where you want to base your work.
  * This is usually the dev branch.
  * Only target release branches if you are certain your fix must be on that
    branch.
  * To quickly create a topic branch based on master, run `git checkout -b
    fix/dev/my_contribution master`. Please avoid working directly on the
    `dev` (or `master`) branch.
* Make commits of logical and atomic units. Write a [good commit message][commit]!
* Make sure you have added the necessary tests for your changes.

## Submitting Changes

* Push your changes to a topic branch in your fork of the repository.
* Submit a pull request to the repository.
* The core team looks at pull requests as soon as possible, but no maximum waiting time can be given here.
* After feedback has been given we expect responses within two weeks. After two
  weeks we may close the pull request if it isn't showing any activity.

[commit]: http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html---
title: 'CoPro: a data-driven modelling framework for conflict risk projections'
tags:
  - Python
  - climate change
  - projections
  - conflict
  - climate security
  - water
  - risk
authors:
  - name: Jannis M. Hoch^[corresponding author]
    orcid: 0000-0003-3570-6436
    affiliation: 1
  - name: Sophie de Bruin
    orcid: 0000-0003-3429-349X
    affiliation: "1, 2"
  - name: Niko Wanders
    orcid: 0000-0002-7102-5454
    affiliation: 1
affiliations:
 - name: Department of Physical Geography, Utrecht University, Utrecht, the Netherlands
   index: 1
 - name: PBL Netherlands Environmental Assessment Agency, the Hague, the Netherlands
   index: 2
date: 18 February 2021
bibliography: bibliography.bib
---

# Summary

Climate change and environmental degradation are increasingly recognized as factors that can contribute to conflict risk under specific conditions.
In light of predicted shifts in climate patterns and the potentially resulting battle for increasingly scarce resources, it is widely acknowledged that there is an actual risk of increased armed conflict. To efficiently plan and implement adaptation and mitigation measures, it is key to first obtain an understanding of conflict drivers and spatial conflict risk distribution. And second, conflict risk needs to be projected to a given point in the future to be able to prepare accordingly. With CoPro, building and running models investigating the interplay between conflict and climate is made easier. By means of a clear workflow, maps of conflict risk for today as well as the future can be produced. Despite the structured workflow, CoPro caters for a variety of settings and input data, thereby capturing the multitude of facets of the climate-environment-conflict nexus.

# Statement of need 

There is increasing consensus that climate change can exacerbate the risk of (armed) conflict [@koubi2019climate; @mach2019climate]. Nevertheless, making (operational) projections of conflict risk is still challenging due to several reasons [@cederman2017predicting]. Building upon recent, similar approaches to use data-driven models [@colaresi2017robot] and statistical approaches [@witmer2017subnational; @hegre2016forecasting], CoPro is a novel, fully open, and extensible Python-model facilitating the set-up, execution, and evaluation of machine-learning models predicting conflict risk. CoPro provides a structured workflow including pre- and post-processing tools, making it accessible to all levels of experience. Such a user-friendly tool is needed not only to integrate the different disciplines, but also to extend the modeling approach with new insights and data - after all, the established links between climate and societal factors with conflict are still weak [@koubi2019climate; @mach2019climate]. In addition to scholarly explorations of the inter-dependencies and the importance of various conflict drivers, model output such as maps of spatially-disaggregated projected conflict risk can be an invaluable input to inform the decision-making process in affected regions.

Since conflicts are of all times and not limited to specific regions or countries, CoPro is designed with user-flexibility in mind. Therefore, the number and variables provided to the model is not specified, allowing for bespoke model designs. Depending on the modeling exercise and data used, several machine-learning models and pre-processing algorithms are available in CoPro. In its current form, the supervised learning techniques support vector classifier, k-neighbors classifier, and random-forest classifier are implemented. Catering for different model designs is of added value because of the non-linear and sometimes irrational - 'law-breaking' [@cederman2017predicting] - nature of conflicts. On top of that, the analyses can be run at any spatial scale, allowing for better identification of sub-national drivers of conflict risk. After all, conflict onset and conflicts are often limited to specific areas where driving factors coincide. 

Since the replicability of scientific results is important when developing forecast and projection models [@hegre2017introduction], CoPro produces reproducible output using transparent models. Hence, by making model code openly available and by including dedicated features in the model, we hope to advance the existing body of tools developed to project conflict risk.

These functionalities altogether make CoPro suited for both 'quick-and-dirty' and in-depth analyses of the relative importances of climate, environmental, and societal drivers as well as for assessments how conflict risk can change both in time and space.

# Acknowledgements
This research was supported by a Pathways to Sustainability Acceleration Grant from the Utrecht University.
We kindly acknowledge the valuable contributions from all partners at PBL, PRIO (Peace Research Institute Oslo), Uppsala University, and Utrecht University.

# References
===============
CoPro
===============

Welcome to CoPro, a machine-learning tool for conflict risk projections based on climate, environmental, and societal drivers.

.. image:: https://travis-ci.com/JannisHoch/copro.svg?branch=dev
    :target: https://travis-ci.com/JannisHoch/copro

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://github.com/JannisHoch/copro/blob/dev/LICENSE

.. image:: https://readthedocs.org/projects/copro/badge/?version=latest
    :target: https://copro.readthedocs.io/en/latest/?badge=latest

.. image:: https://img.shields.io/github/v/release/JannisHoch/copro
    :target: https://github.com/JannisHoch/copro/releases/tag/v0.0.8

.. image:: https://zenodo.org/badge/254407279.svg
    :target: https://zenodo.org/badge/latestdoi/254407279

.. image:: https://badges.frapsoft.com/os/v2/open-source.svg?v=103
    :target: https://github.com/ellerbrock/open-source-badges/

.. image:: https://joss.theoj.org/papers/1f03334e56413ff71f65092ecc609aa4/status.svg
    :target: https://joss.theoj.org/papers/1f03334e56413ff71f65092ecc609aa4

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/JannisHoch/copro/dev?filepath=%2Fexample%2Fnb_binder.ipynb

Model purpose
--------------

As primary model output, CoPro provides maps of conflict risk.

To that end, it employs observed conflicts as target data together with (user-provided) socio-economic and environmental sample data to train different classifiers (RFClassifier, kNearestClassifier, and Support Vector Classifier).
While the samples have the units of the data, the target value is converted to Boolean, where a 0 indicates no conflict occurrence and 1 indicates occurrence.
To capture the geographical variability of conflict and socio-environmental drivers, the model is spatially explicit and calculates conflict risk at a (user-specified) aggregation level.
This way, the model can also capture the relevant sub-national variability of conflict and conflict drivers.
Model robustness is determined using a split-sample test where a part of the data is used to train the model, while the other part is used to evaluate the outcome. 
Throughout this process, the geographical unit is tracked to be able to map the resulting conflict risk to the correct areas.

In addition to the calculation of conflict risk, can the model, for instance, be used to make scenario projections, evaluate the relative feature importances, or benchmark different datasets.

All in all, CoPro supports the mapping of current and future areas at risk of conflict, while also facilitating obtaining a better understanding of the underlying processes.

Installation
----------------

To install copro, first clone the code from GitHub. It is advised to create an individual python environment first. 
You can then install the model package into this environment.

To do so, you need to have Anaconda or Miniconda installed. For installation guidelines, see `here <https://docs.anaconda.com/anaconda/install/>`_.

.. code-block:: console

    $ git clone https://github.com/JannisHoch/copro.git
    $ cd path/to/copro
    $ conda env create -f environment.yml
    $ conda activate copro

To install CoPro in editable mode in this environment, run this command next in the CoPro-folder:

.. code-block:: console

    $ pip install -e .

When using Jupyter Notebook, it can be handy to have the copro environment available. It can be installed into Jupyter Notebook with the following command:

.. code-block:: console

    $ python -m ipykernel install --name=copro

Command-line script
--------------------

To be able to run the model, the conda environment has to be activated first.

.. code-block:: console

    $ conda activate copro

To run the model from command line, a command line script is provided. The usage of the script is as follows:

.. code-block:: console

    Usage: copro_runner [OPTIONS] CFG

    Main command line script to execute the model. 
    All settings are read from cfg-file.
    One cfg-file is required argument to train, test, and evaluate the model.
    Multiple classifiers are trained based on different train-test data combinations.
    Additional cfg-files for multiple projections can be provided as optional arguments, whereby each file corresponds to one projection to be made.
    Per projection, each classifiers is used to create separate projection outcomes per time step (year).
    All outcomes are combined after each time step to obtain the common projection outcome.

    Args:     CFG (str): (relative) path to cfg-file

    Options:
    -plt, --make_plots        add additional output plots
    -v, --verbose             command line switch to turn on verbose mode

This help information can be also accessed with

.. code-block:: console

    $ copro_runner --help

All data and settings are retrieved from the settings-file (cfg-file) which needs to be provided as inline argument.

In case issues occur, updating ``setuptools`` may be required.

.. code-block:: console

    $ pip3 install --upgrade pip setuptools

Example data
----------------

Example data for demonstration purposes can be downloaded from `Zenodo <https://zenodo.org/record/4297295>`_.
To facilitate this process, the bash-script ``download_example_data.sh`` can be called in the example folder under `/_scripts`.

With this (or other) data, the provided configuration-files (cfg-files) can be used to perform a reference run or a projection run. 
All output is stored in the output directory specified in the cfg-files. 
In the output directory, two folders are created: one name `_REF` for output from the reference run, and `_PROJ` for output for projections.

Jupyter notebooks
^^^^^^^^^^^^^^^^^^

There are multiple jupyter notebooks available to guide you through the model application process step-by-step.

It is possible to execute the notebooks cell-by-cell and explore the full range of possibilities.
Note that in this case the notebooks need to be run in the right order as some temporary files will be saved to file in one notebook and loaded in another!
This is due to the re-initalization of the model at the beginning of each notebook and resulting deletion of all files in existing output folders.

The notebooks are also used to exemplify the `Workflow <https://copro.readthedocs.io/en/latest/examples/index.html>`_ of CoPro.

Command-line
^^^^^^^^^^^^^^^^^^

While the notebooks are great for exploring, the command line script is the envisaged way to use CoPro.

To only test the model for the reference situation and one projection, the cfg-file for the reference run is the required argument.
This cfg-file needs to point to the cfg-file of the projection in turn.

.. code-block:: console

    $ cd path/to/copro/example
    $ copro_runner example_settings.cfg

Alternatively, the same commands can be executed using a bash-file.

.. code-block:: console

    $ cd path/to/copro/example/_scripts
    $ sh run_command_line_script.sh

Validation
^^^^^^^^^^^^^^^^^^

The reference model makes use of the `UCDP Georeferenced Event Dataset <https://ucdp.uu.se/downloads/index.html#ged_global>`_ for observed conflict. 
The selected classifier is trained and validated against this data.

Main validation metrics are the ROC-AUC score as well as accuracy, precision, and recall. 
All metrics are reported and written to file per model evaluation.

With the example data downloadable from `Zenodo <https://zenodo.org/record/4297295>`_, a ROC-AUC score of above 0.8 can be obtained. 
Note that with additional and more explanatory sample data, the score will most likely increase.

.. figure:: docs/_static/roc_curve.png

Additional ways to validate the model are showcased in the `Workflow <https://copro.readthedocs.io/en/latest/examples/index.html>`_.

Documentation
---------------

Extensive model documentation including full model API description can be found at http://copro.rtfd.io/

Code of conduct and Contributing
---------------------------------

The project welcomes contributions from everyone! 
To make collaborations as pleasant as possible, we expect contributors to the project to abide by the Code of Conduct.

License
--------

CoPro is released under the MIT license.

Authors
----------------

* Jannis M. Hoch (Utrecht University)
* Sophie de Bruin (Utrecht University, PBL)
* Niko Wanders (Utrecht University)

Corresponding author: Jannis M. Hoch (j.m.hoch@uu.nl)
Model execution
=========================

To be able to run the model, the conda environment has to be activated first.

.. code-block:: console

    $ conda activate copro

.. _script:

Runner script
----------------

To run the model, a command line script is provided. The usage of the script is as follows:

.. code-block:: console

    Usage: copro_runner [OPTIONS] CFG

    Main command line script to execute the model. 
    All settings are read from cfg-file.
    One cfg-file is required argument to train, test, and evaluate the model.
    Multiple classifiers are trained based on different train-test data combinations.
    Additional cfg-files for multiple projections can be provided as optional arguments, whereby each file corresponds to one projection to be made.
    Per projection, each classifiers is used to create separate projection outcomes per time step (year).
    All outcomes are combined after each time step to obtain the common projection outcome.

    Args:     CFG (str): (relative) path to cfg-file

    Options:
    -plt, --make_plots        add additional output plots
    -v, --verbose             command line switch to turn on verbose mode

Help information can be accessed with

.. code-block:: console

    $ copro_runner --help

All data and settings are retrieved from the configuration-file (``cfg-file``, see :ref:`Settings` ) which needs to be provided as command line argument.
In the cfg-file, the various settings of the simulation are defined.

A typical command would thus look like this:

.. code-block:: console

    $ copro_runner settings.cfg

In case issues occur, updating ``setuptools`` may be required.

.. code-block:: console

    $ pip3 install --upgrade pip setuptools

Binder
--------

There is also a notebook running on `Binder <https://mybinder.org/v2/gh/JannisHoch/copro/dev?filepath=%2Fexample%2Fnb_binder.ipynb>`_. 

Please check it out to go through the model execution step-by-step and interactively explore the functionalities of CoPro.
.. _settings:

Settings
=========================

The cfg-file
----------------

The main model settings need to be specified in a configuration file (``cfg-file``). 
This file looks like this.

.. code-block:: console

    [general]
    input_dir=./path/to/input_data
    output_dir=./path/to/store/output
    # 1: all data. 2: leave-one-out model. 3: single variable model. 4: dubbelsteenmodel
    # Note that only 1 supports sensitivity_analysis
    model=1
    verbose=True

    [settings]
    # start year
    y_start=2000
    # end year
    y_end=2012

    [PROJ_files]
    # cfg-files
    proj_nr_1=./path/to/projection/settings_proj.cfg

    [pre_calc]
    # if nothing is specified, the XY array will be stored in output_dir
    # if XY already pre-calculated, then provide path to npy-file
    XY=

    [extent]
    shp=folder/with/polygons.shp

    [conflict]
    # either specify path to file or state 'download' to download latest PRIO/UCDP dataset
    conflict_file=folder/with/conflict_data.csv
    min_nr_casualties=1
    # 1=state-based armed conflict. 2=non-state conflict. 3=one-sided violence
    type_of_violence=1,2,3

    [climate]
    shp=folder/with/climate_zones.shp
    # define either one or more classes (use abbreviations!) or specify nothing for not filtering
    zones=
    code2class=folder/with/classification_codes.txt

    [data]
    # specify the path to the nc-file, whether the variable shall be log-transformed (True, False), and which statistical function should be applied
    # these three settings need to be separated by a comma
    # NOTE: variable name here needs to be identical with variable name in nc-file
    # NOTE: only statistical functions supported by rasterstats are valid
    precipitation=folder/with/precipitation_data.nc,False,mean
    temperature=folder/with/temperature_data.nc,False,mean
    population=folder/with/population_data.nc,True,sum

    [machine_learning]
    # choose from: MinMaxScaler, StandardScaler, RobustScaler, QuantileTransformer
    scaler=QuantileTransformer
    # choose from: NuSVC, KNeighborsClassifier, RFClassifier
    model=RFClassifier
    train_fraction=0.7
    # number of repetitions
    n_runs=10

.. note::

    All paths for ``input_dir``, ``output_dir``, and in ``[PROJ_files]`` are relative to the location of the cfg-file. 

.. important::

    Empty spaces should be avoided in the cfg-file, besides for those lines commented out with '#'.

The sections
----------------

Here, the different sections are explained briefly. 

[general]
^^^^^^^^^^^^^^^^

``input_dir``: (relative) path to the directory where the input data is stored. This requires all input data to be stored in one main folder, sub-folders are possible.

``output_dir``: (relative) path to the directory where output will be stored. 
If the folder does not exist yet, it will be created. 
CoPro will automatically create the sub-folders ``_REF`` for output for the reference run, and ``_PROJ`` for output from the (various) projection runs.

``model``: the type of simulation to be run can be specified here. Currently, for different models are available:

    1. 'all data': all variable values are used to fit the model and predict results.
    2. 'leave one out': values of each variable are left out once, resulting in n-1 runs with n being the number of variables. This model can be used to identify the relative influence of one variable within the variable set/
    3. 'single variables': each variable is used as sole predictor once. With this model, the explanatory power of each variable on its own can be assessed.
    4. 'dubbelsteen': the relation between variables and conflict are abolished by shuffling the binary conflict data randomly. By doing so, the lower boundary of the model can be estimated.

.. note::

    All model types except 'all_data' will be deprecated in a future release.

``verbose``: if True, additional messages will be printed.

[settings]
^^^^^^^^^^^^^^^^

``y_start``: the start year of the reference run.

``y_end``: the end year of the reference run. 
The period between ``y_start`` and ``y_end`` will be used to train and test the model.

``y_proj``: the end year of the projection run.
The period between ``y_end`` and ``y_proj`` will be used to make annual projections.

[PROJ_files]
^^^^^^^^^^^^^^^^

A key section. Here, one (slightly different) cfg-file per projection needs to be provided. 
This way, multiple projection runs can be defined from within the "main" cfg-file.

The conversion is that the projection name is defined as value here.
For example, the projections "SSP1" and "SSP2" would be defined as

.. code-block:: console

    SSP1=/path/to/ssp1.cfg
    SSP2=/path/to/ssp2.cfg

A cfg-file for a projection is shorter than the main cfg-file used as command line argument and looks like this:

.. code-block:: console

    [general]
    input_dir=./path/to/input_data
    verbose=True

    [settings]
    # year for which projection is to be made
    y_proj=2050

    [data]
    # specify the path to the nc-file, whether the variable shall be log-transformed (True, False), and which statistical function should be applied
    # these three settings need to be separated by a comma
    # NOTE: variable name here needs to be identical with variable name in nc-file
    # NOTE: only statistical functions supported by rasterstats are valid
    precipitation=folder/with/precipitation_data.nc,False,mean
    temperature=folder/with/temperature_data.nc,False,mean
    population=folder/with/population_data.nc,True,sum

[pre_calc]
^^^^^^^^^^^^^^^^

``XY``: if the XY-data was already pre-computed in a previous run and stored as npy-file, it can be specified here and will be loaded from file to save time. 
If nothing is specified, the model will save the XY-data by default to the output directory as ``XY.npy``.

[extent]
^^^^^^^^^^^^^^^^

``shp``: the provided shape-file defines the boundaries for which the model is applied. 
At the same time, it also defines at which aggregation level the output is determined.

.. note:: 

    The shp-file can contain multiple polygons covering the study area. Their size defines the output aggregation level. It is also possible to provide only one polygon, but model behaviour is not well tested for this case.

[conflict]
^^^^^^^^^^^^^^^^

``conflict_file``: path to the csv-file containing the conflict dataset. 
It is also possible to define ``download``, then the latest conflict dataset (currently version 20.1) is downloaded and used as input.

``min_nr_casualties``: minimum number of reported casualties required for a conflict to be considered in the model.

``type_of_violence``: the types of violence to be considered can be specified here. 
Multiple values can be specified. Types of violence are:

    1. 'state-based armed conflict': a contested incompatibility that concerns government and/or territory where the use of armed force between two parties, of which at least one is the government of a state, results in at least 25 battle-related deaths in one calendar year.
    2. 'non-state conflict': the use of armed force between two organized armed groups, neither of which is the government of a state, which results in at least 25 battle-related deaths in a year.
    3. 'one-sided violence': the deliberate use of armed force by the government of a state or by a formally organized group against civilians which results in at least 25 deaths in a year.

.. important::

    CoPro currently only works with UCDP data.

[climate]
^^^^^^^^^^^^^^^^

``shp``: the provided shape-file defines the areas of the different Köppen-Geiger climate zones.

``zones``: abbreviations of the climate zones to be considered in the model.
Can either be 'None' or one or multiple abbreviations.

``code2class``: converting the abbreviations to class-numbers used in the shp-file.

.. warning:: 

    The code2class-file should not be altered!

[data]
^^^^^^^^^^^^^^^^

In this section, all variables to be used in the model need to be provided. 
The paths are relative to ``input_dir``.
Only netCDF-files with annual data are supported.

The main convention is that the name of the file agrees with the variable name in the file.
For example, if the variable ``precipitation`` is provided in a nc-file, this should be noted as follows

.. code-block:: console

    [data]
    precipitation=folder/with/precipitation_data.nc

CoPro furthermore requires information whether the values sampled from a file are ought to be log-transformed.

Besides, it is possible to define a statistical function that is applied when sampling from file per polygon of the ``shp-file``.
CoPro makes use of the ``zonal_stats`` function available within `rasterstats <https://pythonhosted.org/rasterstats/rasterstats.html>`_.

To determine the log-scaled mean value of precipitation per polygon, the following notation is required:

.. code-block:: console

    [data]
    precipitation=folder/with/precipitation_data.nc,False,mean

[machine_learning]
^^^^^^^^^^^^^^^^^^^^

``scaler``: the scaling algorithm used to scale the variable values to comparable scales. 
Currently supported are 

    - `MinMaxScaler <https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html>`_;
    - `StandardScaler <https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html>`_;
    - `RobustScaler <https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.RobustScaler.html>`_;
    - `QuantileTransformer <https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.QuantileTransformer.html>`_.

``model``: the machine learning algorithm to be applied. 
Currently supported are 

    - `NuSVC <https://scikit-learn.org/stable/modules/generated/sklearn.svm.NuSVC.html>`_; 
    - `KNeighborsClassifier <https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html>`_;
    - `RFClassifier <https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html>`_.

``train_fraction``: the fraction of the XY-data to be used to train the model. 
The remaining data (1-train_fraction) will be used to predict and evaluate the model.

``n_runs``: the number of classifiers to use.Installation
=========================

From GitHub
------------

To install CoPro from GitHub, first clone the code. It is advised to create a separate environment first. 

.. note::

    We recommend to use Anaconda or Miniconda to install CoPro as this was used to develop and test the model.
    For installation instructions, see `here <https://docs.anaconda.com/anaconda/install/>`_.

.. code-block:: console

    $ git clone https://github.com/JannisHoch/copro.git
    $ cd path/to/copro
    $ conda env create -f environment.yml

It is now possible to activate this environment with

.. code-block:: console

    $ conda activate copro

To install CoPro in editable mode in this environment, run this command next in the CoPro-folder:

.. code-block:: console

    $ pip install -e .

From PyPI
------------

To install CoPro directly from PyPI, use the following command.

.. code-block:: console

    pip install copro

From conda
------------

.. todo::

    This is not yet supported. Feel invited to provide a pull request enabling installation via conda.Output
=========================

Output folder structure
---------------------------

All output is stored in the output folder as specified in the configurations-file (cfg-file) under [general].

.. code-block:: console

    [general]
    output_dir=./path/to/store/output

By default, CoPro creates two sub-folders: ``_REF`` and ``_PROJ``. In the latter, another sub-folder will be created per projection defined in the cfg-file.
In the example below, this would be the folders ``/_PROJ/SSP1`` and ``/_PROJ/SSP2``.

.. code-block:: console

    [PROJ_files]    
    SSP1=/path/to/ssp1.cfg
    SSP2=/path/to/ssp2.cfg

List of output files
---------------------------

.. important:: 

    Not all model types provide the output mentioned below. If the 'leave-one-out' or 'single variable' model are selected, only the metrics are stored to a csv-file.

_REF
^^^^^^

In addition to the output files listed below, the cfg-file is automatically copied to the _REF folder.

``selected_polygons.shp``: Shapefile containing all remaining polygons after selection procedure.

``selected_conflicts.shp``: Shapefile containing all remaining conflict points after selection procedure,

``XY.npy``: NumPy-array containing geometry, ID, and scaled data of sample (X) and target data (Y). 
Can be provided in cfg-file to safe time in next run; file can be loaded with numpy.load().

``raw_output_data.npy``: NumPy-array containing each single prediction made in the reference run.
Will contain multiple predictions per polygon. File can be loaded with numpy.load().

``evaluation_metrics.csv``: Various evaluation metrics determined per repetition of the split-sample tests.
File can e.g. be loaded with pandas.read_csv().

``feature_importance.csv``: Importance of each model variable in making projections.
This is a property of RF Classifiers and thus only obtainable if RF Classifier is used.

``permutation_importance.csv``: Mean permutation importance per model variable.
Computed with sklearn.inspection.permutation_importance_.

``ROC_data_tprs.csv`` and ``ROC_data_aucs.csv``: False-positive rates respectively Area-under-curve values per repetition of the split-sample test.
Files can e.g. be loaded with pandas.read_csv() and can be used to later plot ROC-curve.

``output_for_REF.geojson``: GeoJSON-file containing resulting conflict risk estimates per polygon based on out-of-sample projections of _REF run.

.. _sklearn.inspection.permutation_importance: https://scikit-learn.org/stable/modules/generated/sklearn.inspection.permutation_importance.html

Conflict risk per polygon
""""""""""""""""""""""""""

At the end of all model repetitions, the resulting ``raw_output_data.npy`` file contains multiple out-of-sample predictions per polygon.
By aggregating results per polygon, it is possible to assess model output spatially as stored in ``output_for_REF.geojson``. 

The main output metrics are calculated per polygon and saved to ``output_per_polygon.shp``:

1. nr_predictions: the number of predictions made;
2. nr_correct_predictions: the number of correct predictions made;
3. nr_observed_conflicts: the number of observed conflict events;
4. nr_predicted_conflicts: the number of predicted conflicts;
5. min_prob_1: minimum probability of conflict in all repetitions;
6. probability_of_conflict (POC): probability of conflict averaged over all repetitions;
7. max_prob_1: maximum probability of conflict in all repetitions;
8. fraction_correct_predictions (FOP): ratio of the number of correct predictions over the total number of predictions made;
9. chance_of_conflict: ratio of the number of conflict predictions over the total number of predictions made.

_PROJ
^^^^^^

Per projection, CoPro creates one output file per projection year.

``output_in_<YEAR>``: GeoJSON-file containing model output per polygon averaged over all classifier instances per YEAR of the projection.
The number of instances is set with ``n_runs`` in ``[machine_learning]`` section.

Conflict risk per polygon
""""""""""""""""""""""""""

During the projection run, each classifier instances produces its own output per YEAR.
CoPro merges these outputs into one ``output_in_<YEAR>.geojson`` file. 

As there are no observations available for the projection period, the output metrics differ from the reference run:

1. nr_predictions: the number of predictions made, ie. number of classifier instances;
2. nr_predicted_conflicts: the number of predicted conflicts.
3. min_prob_1: minimum probability of conflict in all outputs of classifier instances.
4. probability_of_conflict (POC): probability of conflict averaged over all outputs of classifier instances.
5. max_prob_1: maximum probability of conflict in all outputs of classifier instances;
6. chance_of_conflict: ratio of the number of conflict predictions over the total number of predictions made.
Postprocessing
=========================

There are several command line scripts available for post-processing. 
In addition to quick plots to evaluate model output, they also produce files for use in bespoke plotting and analysis scripts.

The scripts are located under ``/copro/scripts/postprocessing``.

The here shown help print-outs can always be accessed with 

.. code-block:: console

    python <SCRIPT_FILE_NAME> --help

plot_value_over_time.py
------------------------

.. code-block:: console

    Usage: python plot_value_over_time.py [OPTIONS] INPUT_DIR OUTPUT_DIR

        Quick and dirty function to plot the develoment of a column in the
        outputted geoJSON-files over time. The script uses all geoJSON-files
        located in input-dir and retrieves values from them. Possible to plot
        obtain development for multiple polygons (indicated via their ID) or
        entire study area. If the latter, then different statistics can be chosen
        (mean, max, min, std).

        Args:     
            input-dir (str): path to input directory with geoJSON-files located per projection year. 
            output-dir (str): path to directory where output will be stored.

        Output:     
            a csv-file containing values per time step.     
            a png-file showing development over time.

        Options:
            -id, --polygon-id TEXT
            -s, --statistics TEXT     which statistical method to use (mean, max, min,
                                        std). note: has only effect if with "-id all"!

            -c, --column TEXT         column name
            -t, --title TEXT          title for plot and file_object name
            --verbose / --no-verbose  verbose on/off

avg_over_time.py
-----------------

.. code-block:: console

    Usage: python avg_over_time.py [OPTIONS] INPUT_DIR OUTPUT_DIR SELECTED_POLYGONS

        Post-processing script to calculate average model output over a user-
        specifeid period or all output geoJSON-files stored in input-dir.
        Computed average values can be outputted as geoJSON-file or png-file or both.

        Args:     
            input_dir: path to input directory with geoJSON-files located per projection year.     
            output_dir (str): path to directory where output will be stored.     
            selected_polygons (str): path to a shp-file with all polygons used in a CoPro run.

        Output:     
            geoJSON-file with average column value per polygon (if geojson is set).     
            png-file with plot of average column value per polygon (if png is set)

        Options:
            -t0, --start-year INTEGER
            -t1, --end-year INTEGER
            -c, --column TEXT          column name
            --geojson / --no-geojson   save output to geojson or not
            --png / --no-png           save output to png or not
            --verbose / --no-verbose   verbose on/off

plot_polygon_vals.py
-----------------------

.. code-block:: console

    Usage: python plot_polygon_vals.py [OPTIONS] FILE_OBJECT OUTPUT_DIR

        Quick and dirty function to plot the column values of a geojson file with
        minimum user input, and save plot. Mainly used for quick inspection of
        model output in specific years.

        Args:     
            file-object (str): path to geoJSON-file whose values are to be plotted.     
            output-dir (str): path to directory where plot will be saved.

        Output:     
            a png-file of values per polygon.

        Options:
            -c, --column TEXT           column name
            -t, --title TEXT            title for plot and file_object name
            -v0, --minimum-value FLOAT
            -v1, --maximum-value FLOAT
            -cmap, --color-map TEXT

geojson2gif.py
---------------

.. code-block:: console

    Usage: python geojson2gif.py [OPTIONS] INPUT_DIR OUTPUT_DIR

        Function to convert column values of all geoJSON-files in a directory into
        one GIF-file. The function provides several options to modify the design
        of the GIF-file. The GIF-file is based on png-files of column value per
        geoJSON-file.  It is possible to keep these png-file as simple plots of
        values per time step.

    Args:     
        input-dir (str): path to directory where geoJSON-files are stored.     
        output_dir (str): path to directory where GIF-file will be stored.

    Output:     
        GIF-file with animated column values per input geoJSON-file.

    Options:
        -c, --column TEXT           column name
        -cmap, --color-map TEXT
        -v0, --minimum-value FLOAT
        -v1, --maximum-value FLOAT
        --delete / --no-delete      whether or not to delete png-files

CoPro
=========================

This is the documentation of CoPro, a machine-learning tool for conflict risk projections.

A software description paper was published in `JOSS <https://doi.org/10.21105/joss.02855>`_.

.. image:: https://travis-ci.com/JannisHoch/copro.svg?branch=dev
    :target: https://travis-ci.com/JannisHoch/copro

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://github.com/JannisHoch/copro/blob/dev/LICENSE

.. image:: https://readthedocs.org/projects/copro/badge/?version=latest
    :target: https://copro.readthedocs.io/en/latest/?badge=latest

.. image:: https://img.shields.io/github/v/release/JannisHoch/copro
    :target: https://github.com/JannisHoch/copro/releases/tag/v0.0.8

.. image:: https://zenodo.org/badge/254407279.svg
    :target: https://zenodo.org/badge/latestdoi/254407279

.. image:: https://badges.frapsoft.com/os/v2/open-source.svg?v=103
    :target: https://github.com/ellerbrock/open-source-badges/

.. image:: https://joss.theoj.org/papers/10.21105/joss.02855/status.svg
    :target: https://doi.org/10.21105/joss.02855

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/JannisHoch/copro/dev?filepath=%2Fexample%2Fnb_binder.ipynb

Main goal
---------------
With CoPro it is possible to apply machine-learning techniques to make projections of future areas at risk. CoPro was developed with a rather clear
application in mind, unravelling the interplay of socio-economic development, climate change, and conflict occurrence. 
Nevertheless, we put a lot of emphasis on making it flexible. 
We hope that other, related questions about climate and conflict can be tackled as well, and that process understanding is deepened further.

Contents
---------------
.. toctree::
   :numbered:
   :maxdepth: 2

   Installation <Installation>
   Execution <Execution>
   Settings <Settings>
   Workflow <examples/index>
   Output <Output>
   Post-processing <Postprocessing>
   API <api/index>

Authors
----------------

* Jannis M. Hoch (Utrecht University)
* Sophie de Bruin (Utrecht University, PBL)
* Niko Wanders (Utrecht University)

Corresponding author: Jannis M. Hoch (j.m.hoch@uu.nl)

Indices and tables
-------------------

* :ref:`search`
* :ref:`genindex`
.. * :ref:`modindex`

.. _workflow:

Workflow
=========

This page provides a short example workflow in Jupyter Notebooks. 
It is designed such that the main steps, features, assumptions, and outcomes of CoPro become clear.

As model input data, the data set downloadable from `Zenodo <https://zenodo.org/record/4297295>`_ was used. 

Even though the model can be perfectly executed using notebooks, the main (and more convenient) way of model execution is the command line script (see :ref:`script`).

An interactive version of the content shown here can be accessed via `Binder <https://mybinder.org/v2/gh/JannisHoch/copro/dev?filepath=%2Fexample%2Fnb_binder.ipynb>`_.

.. toctree::
    :maxdepth: 1
    :numbered:
    :glob:

    ./*Model evaluation
=================================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   evaluation.init_out_dict
   evaluation.fill_out_dict
   evaluation.init_out_df
   evaluation.fill_out_df
   evaluation.evaluate_prediction
   evaluation.polygon_model_accuracy
   evaluation.init_out_ROC_curve
   evaluation.save_out_ROC_curve
   evaluation.calc_correlation_matrix
   evaluation.get_feature_importance
   evaluation.get_permutation_importanceXY-Data
=================================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   data.initiate_XY_data
   data.initiate_X_data
   data.fill_XY
   data.fill_X_sample
   data.fill_X_conflict
   data.split_XY_data
   data.neighboring_polys
   data.find_neighborsSelecting polygons and conflicts
=================================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   selection.select
   selection.filter_conflict_properties
   selection.select_period
   selection.clip_to_extent
   selection.climate_zoningMachine learning
=================================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   machine_learning.define_scaling
   machine_learning.define_model
   machine_learning.split_scale_train_test_split
   machine_learning.fit_predict
   machine_learning.pickle_clf
   machine_learning.load_clfsAuxiliary functions
=================================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   utils.print_model_info
   utils.get_geodataframe
   utils.show_versions
   utils.parse_settings
   utils.parse_projection_settings
   utils.determine_projection_period
   utils.make_output_dir
   utils.download_UCDP
   utils.initiate_setup
   utils.create_artificial_Y
   utils.global_ID_geom_info
   utils.get_conflict_datapoints_only
   utils.save_to_csv
   utils.save_to_npyWork with conflict data
=================================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   conflict.conflict_in_year_bool
   conflict.conflict_in_previous_year
   conflict.read_projected_conflict
   conflict.calc_conflicts_nb
   conflict.get_poly_ID
   conflict.get_poly_geometry
   conflict.split_conflict_geom_data
   conflict.get_pred_conflict_geometryThe model pipeline
====================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   pipeline.create_XY
   pipeline.prepare_ML
   pipeline.run_reference
   pipeline.run_predictionThe various models
====================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   models.all_data
   models.leave_one_out
   models.single_variables
   models.dubbelsteen
   models.predictive

.. note::

    The 'leave_one_out', 'single_variables', and 'dubbelsteen' models are only tested in beta-state. 
    They will most likely be deprecated in near future.API docs
========

This section contains the Documentation of the Application Programming
Interface (API) of 'copro'.

.. toctree::
   :numbered:
   :maxdepth: 1

   pipeline
   models
   selection
   machine_learning
   variables
   XYdata
   conflict
   evaluation
   plotting
   utilsPlotting
=================================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   plots.selected_polygons
   plots.selected_conflicts
   plots.metrics_distribution
   plots.correlation_matrix
   plots.plot_ROC_curve_n_times
   plots.plot_ROC_curve_n_meanVariable values
=================================

.. currentmodule:: copro

.. autosummary::
   :toctree: generated/
   :nosignatures:

   variables.nc_with_float_timestamp
   variables.nc_with_continous_datetime_timestamp

.. warning::

   Reading files with a float timestamp will most likely be deprecated in near future.