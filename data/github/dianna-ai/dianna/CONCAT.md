Planned DIANNA developments
| fair-software.eu recommendations | |
| :-- | :--  |
| license                      | [![github license badge](https://img.shields.io/github/license/dianna-ai/dianna)](https://github.com/dianna-ai/dianna) |
| community registry           | [![RSD](https://img.shields.io/badge/rsd-dianna-00a3e3.svg)](https://www.research-software.nl/software/dianna) [![workflow pypi badge](https://img.shields.io/pypi/v/dianna.svg?colorB=blue)](https://pypi.python.org/project/dianna/) |
| citation                     | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5592606.svg)](https://doi.org/10.5281/zenodo.5592606) |
| howfairis                          | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) |
| **Other best practices**           | &nbsp; |
| Static analysis                    | [![workflow scq badge](https://sonarcloud.io/api/project_badges/measure?project=dianna-ai_dianna&metric=alert_status)](https://sonarcloud.io/dashboard?id=dianna-ai_dianna) |
| Coverage                           | [![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=dianna-ai_dianna&metric=coverage)](https://sonarcloud.io/dashboard?id=dianna-ai_dianna) |
| **GitHub Actions**                 | &nbsp; |
| Build                              | [![build](https://github.com/dianna-ai/dianna/actions/workflows/build.yml/badge.svg)](https://github.com/dianna-ai/dianna/actions/workflows/build.yml) |
| Citation data consistency               | [![cffconvert](https://github.com/dianna-ai/dianna/actions/workflows/cffconvert.yml/badge.svg)](https://github.com/dianna-ai/dianna/actions/workflows/cffconvert.yml) |
| SonarCloud                         | [![sonarcloud](https://github.com/dianna-ai/dianna/actions/workflows/sonarcloud.yml/badge.svg)](https://github.com/dianna-ai/dianna/actions/workflows/sonarcloud.yml) |
| MarkDown link checker              | [![markdown-link-check](https://github.com/dianna-ai/dianna/actions/workflows/markdown-link-check.yml/badge.svg)](https://github.com/dianna-ai/dianna/actions/workflows/markdown-link-check.yml) |
<!--
title: 'DIANNA: Deep Insight And Neural Network Analysis'
tags:
  - Python
  - explainable AI
  - deep neural networks
  - ONNX
  - benchmark sets
authors:
  - name: Elena Ranguelova^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-9834-1756
    affiliation: 1
  - name: Patrick Bos^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-6033-960X
    affiliation: 1
  - name: Yang Liu^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-1966-8460
    affiliation: 1
  - name: Christiaan Meijer^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-5529-5761
    affiliation: 1
  - name: Leon Oostrum^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0001-8724-8372
    affiliation: 1
affiliations:
 - name: Netherlands eScience Center, Amsterdam, the Netherlands
   index: 1
-->

[![build](https://github.com/dianna-ai/dianna/actions/workflows/build.yml/badge.svg)](https://github.com/dianna-ai/dianna/actions/workflows/build.yml)
[![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=dianna-ai_dianna&metric=coverage)](https://sonarcloud.io/dashboard?id=dianna-ai_dianna)
[![workflow pypi badge](https://img.shields.io/pypi/v/dianna.svg?colorB=blue)](https://pypi.python.org/project/dianna/)
[![supported python versions](https://img.shields.io/pypi/pyversions/dianna)](https://pypi.python.org/project/dianna/)


[![Documentation Status](https://readthedocs.org/projects/dianna/badge/?version=latest)](https://dianna.readthedocs.io/en/latest/?badge=latest)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/5542/badge)](https://bestpractices.coreinfrastructure.org/projects/5542)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5592607.svg)](https://zenodo.org/record/5592607)
[![more badges badge](https://img.shields.io/badge/more-badges-lightgrey)](badges.md)

# DIANNA: Deep Insight And Neural Network Analysis
<!-- TODO: add main points and then expand, see issue https://github.com/dianna-ai/dianna/issues/137 -->

## Why DIANNA? 
<!-- TO DO: edit the proposal text into something clear and simpler -->

Issues:
1.	The properties of the heatmaps are not studied and the human interpretation is intertwined with the XAI’s. Suitable datasets are lacking: the popular MNIST benchmark is too complex for the task (10 classes and no structural content variation). The XAI literature does not consider simple scientific “benchmarks”.
2.	Which is the “best” explainability method? There is no agreement in the XAI community. The libraries offer different subsets of XAI methods not chosen systematically.
3.	The available OSS (not for all methods) implementations support a single DNN format/framework, e.g. iNNvestigate supports only Keras, while Captum supports PyTorch. 
4.	Not many demonstrators of XAI exist, except from LRP and RISE.

Solutions:
1.	To demonstrate the usefulness and properties of the heatmaps on an intuitive level we propose: simple geometrical and simple scientific– subset of LeafSnap datasets. Tree species classification on the LeafSnap data is a good example of problem tackled with both classical Computer Vision and the superior DL method.
2.	Recently, several systematically defined criteria for evaluation of the XAI approaches have been proposed with LIME analyzed as example. Analysis of the state-of-the-art XAI methods will highlight the best.
3.	DIANNA is a library conforming with the ONNX standard. There are many ONNX tools available as OSS including the ONNX model zoo and ONNX converters from Keras and TensorFlow. PyTorch also offers built-in PyTorch to ONNX export.
4.	A web demonstrator will be created in a next phase of the project. 

## Installation 

To install DIANNA directly from the GitHub repository, do:

```console
python3 -m pip install git+https://github.com/dianna-ai/dianna.git
```

For development purposes, when you first clone the repository locally, it may be more convenient to install in editable mode using pip's `-e` flag:

```console
git clone https://github.com/dianna-ai/dianna.git
cd dianna
python3 -m pip install -e .
```


## Datasets
DIANNA comes with simple datasets. Their main goal is to provide intuitive insight into the working of the XAI methods. They can be used as benchmarks for evaluation and comparison of existing and new XAI methods.

### Images
|Dataset|Description|Examples|Generation|
|:-----|:----|:---|:----|
|Binary MNIST | Greyscale images of the digits "1" and "0" - a 2-class subset from the famous [MNIST dataset](http://yann.lecun.com/exdb/mnist/) for handwritten digit classification. |<img width="120" alt="BinaryMNIST" src="https://user-images.githubusercontent.com/3244249/150808267-3d27eae0-78f2-45f8-8569-cb2561f2c2e9.png">| [Binary MNIST dataset generation](https://github.com/dianna-ai/dianna-exploration/tree/main/example_data/dataset_preparation/MNIST)|
|[Simple Geometric (circles and triangles)](https://doi.org/10.5281/zenodo.5012824) <img width="20" alt="Simple Geometric Logo" src="https://user-images.githubusercontent.com/3244249/150808842-d35d741e-294a-4ede-bbe9-58e859483589.png"> | Images of circles and triangles for 2-class geometric shape classificaiton. The shapes of varying size and orientation and the background have varying uniform gray levels.  | <img width="130" alt="SimpleGeometric" src="https://user-images.githubusercontent.com/3244249/150808125-e1576237-47fa-4e51-b01e-180904b7c7f6.png">| [Simple geometric shapes dataset generation](https://github.com/dianna-ai/dianna-exploration/tree/main/example_data/dataset_preparation/geometric_shapes) | 
|[Simple Scientific (LeafSnap30)](https://zenodo.org/record/5061353/)<img width="20" alt="LeafSnap30 Logo" src="https://user-images.githubusercontent.com/3244249/150815639-2da560d4-8b26-4eeb-9ab4-dabf221a264a.png"> | Color images of tree leaves - a 30-class post-processed subset from the LeafSnap dataset for automatic identification of North American tree species.|<img width="600" alt="LeafSnap" src="https://user-images.githubusercontent.com/3244249/150804246-f714e517-641d-48b2-af26-2f04166870d6.png">| [LeafSnap30 dataset generation](https://github.com/dianna-ai/dianna-exploration/blob/main/example_data/dataset_preparation/LeafSnap/)|

### Text
|Dataset|Description|Examples|Generation|
|:-----|:----|:---|:----|
| [Stanford sentiment treebank](https://nlp.stanford.edu/sentiment/index.html)|Dataset for predicting the sentiment, positive or negative, of movie reviews. | _This movie was actually neither that funny, nor super witty._|[Sentiment treebank](https://nlp.stanford.edu/sentiment/treebank.html)|

## ONNX models
<!-- TODO: Add all links, see issue https://github.com/dianna-ai/dianna/issues/135 -->

**We work with ONNX!** ONNX is a great unified neural network standard which can be used to boost reproducible science. Using ONXX for your model also gives you a boost in performance! In case your models are still in another popular DNN (deep neural network) format, here are some simple recipes to convert them:
* pytorch
* tensorflow
* keras
* scikit-learn

And here are links to notebooks showing how we created our models on the benchmark datasets:

### Images
* Binary MNIST model
* Simple Geometric model
* Simple Scientific model

### Text
* Movie reviews model

**_We envision the birth of the ONNX Scientific models zoo soon..._**

## Tutorials
DIANNA supports different data modalities and XAI methods. The table contains links to the relevant XAI method's papers. There are DIANNA [tutorials](./tutorials) covering each supported method and data modality on a least one dataset. Our future plans to expand DIANNA with more data modalities and XAI methods are given at the [ROADMAP.md](./ROADMAP.md).

<!-- see issue: https://github.com/dianna-ai/dianna/issues/142, also related issue: https://github.com/dianna-ai/dianna/issues/148 -->

|Data \ XAI|[RISE](http://bmvc2018.org/contents/papers/1064.pdf)|[LIME](https://www.kdd.org/kdd2016/papers/files/rfp0573-ribeiroA.pdf)|[KernelSHAP](https://proceedings.neurips.cc/paper/2017/file/8a20a8621978632d76c43dfd28b67767-Paper.pdf)|
|:-----|:---|:---|:---|
|Images|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|Text|:white_check_mark:|:white_check_mark:|planned|
|Embedding|coming soon|coming soon|coming soon|
|Timeseries|planned|planned|planned|
|Tabular|planned|planned|planned|

[LRP](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0130140&type=printable) and [PatternAttribution](https://arxiv.org/pdf/1705.05598.pdf) also feature in the top 5 of our thoroughly evaluated XAI methods using objective critera (details in coming blog-post). **Contributing by adding these and more (new) post-hoc explainability methods on ONNX models is very welcome!**

## Reference documentation 

For detailed information on using specific DIANNA functions, please visit the [Sphinx documentation page hosted at Readthedocs](https://dianna.readthedocs.io/en/latest).

## Contributing

If you want to contribute to the development of DIANNA,
have a look at the [contribution guidelines](https://github.com/dianna-ai/dianna/blob/main/CONTRIBUTING.md).

## How to cite us 

If you use this package for your scientific work, please consider citing it as:

    Ranguelova, Elena, Bos, Patrick, Liu, Yang, Meijer, Christiaan, & Oostrum, Leon. (2021). dianna (*[VERSION YOU USED]*). Zenodo. https://zenodo.org/record/5592607

See also the [Zenodo page](https://zenodo.org/record/5592607) for exporting the citation to BibTteX and other formats.

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [NLeSC/python-template](https://github.com/NLeSC/python-template).
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

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
reported by contacting the project team at dianna-ai@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/dianna-ai/dianna/issues) to see if someone already filed the same issue;
2. if your issue search did not yield any relevant results, make a new issue;
3. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/dianna-ai/dianna/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``pytest``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. push your feature branch to (your fork of) the dianna repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
We also have a [developer documentation chapter](https://dianna.readthedocs.io/en/latest/developer_info.html) that you can use for reference.
## Tutorials
This folder contains tutorial notebooks for DIANNA.
A tutorial is avaiable for several combinations of data modality and explanability method.
Click one of the links in the table to directly go to a notebook:

|Data \ XAI|RISE|LIME|KernelSHAP|
|:-----|:---|:---|:---|
|Images|[]()|[link](lime_images.ipynb)|[]()|
|Text|[]()|[link](lime_text.ipynb)|[]()|
************************************************
Developer documentation for the DIANNA project
************************************************

This chapter lists the tools and practices we typically use during development.
Most of our choices are motivated in the `NL eScience Center Guide <https://guide.esciencecenter.nl>`__, specifically the `Python chapter <https://guide.esciencecenter.nl/#/best_practices/language_guides/python>`__ and our `software guide checklist <https://guide.esciencecenter.nl/#/best_practices/checklist>`__.
If you're considering contributing to DIANNA, please have a look and be sure to let us know (e.g. via an issue) if you have any questions.


Development install
-------------------

.. code:: shell

   # Create a virtual environment, e.g. with
   python3 -m venv env

   # activate virtual environment
   source env/bin/activate

   # make sure to have a recent version of pip and setuptools
   python3 -m pip install --upgrade pip setuptools

   # (from the project root directory)
   # install dianna as an editable package and install development dependencies
   python3 -m pip install --no-cache-dir --editable .[dev]

It's also possible to use ``conda`` for maintaining virtual environments; in that case you can still install DIANNA and dependencies with ``pip``.

Dependencies and Package management
-----------------------------------

DIANNA aims to support all Python 3 minor versions that are still
actively maintained, currently:

-  3.7
-  3.8
-  3.9
-  3.10 is still missing pending the availibility of `onnxruntime` in 3.10.

Add or remove Python versions based on availability of dependencies in
all versions. See `the
guide <https://guide.esciencecenter.nl/#/best_practices/language_guides/python>`__
for more information about Python versions.

When adding new dependencies, make sure to do so as follows:

-  Runtime dependencies should be added to ``setup.cfg`` in the
   ``install_requires`` list under ``[options]``.
-  Development dependencies should be added to ``setup.cfg`` in one of
   the lists under ``[options.extras_require]``.

Testing and code coverage
-------------------------

-  Tests should be put in the ``tests`` folder.
-  Take a look at the existing tests and add your own meaningful tests
   (file: ``test_my_module.py``) when you add a feature.
-  The testing framework used is `PyTest <https://pytest.org>`__
-  The project uses `GitHub action
   workflows <https://docs.github.com/en/actions>`__ to automatically
   run tests on GitHub infrastructure against multiple Python versions

   -  Workflows can be found in
      `.github/workflows <https:://github.com/dianna-ai/dianna/.github/workflows/>`__

-  `Relevant section in the
   guide <https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=testing>`__

Running the tests
~~~~~~~~~~~~~~~~~

There are two ways to run tests.

The first way requires an activated virtual environment with the
development tools installed:

.. code:: shell

   pytest -v

The second is to use ``tox``, which must be installed separately (e.g. with ``pip install tox``), but then builds the necessary virtual environments itself by simply running:

.. code:: shell

   tox

Testing with ``tox`` allows for keeping the testing environment separate from your development environment.
The development environment will typically accumulate (old) packages during development that interfere with testing; this problem is avoided by testing with ``tox``.

Running linters locally
-----------------------

For linting we use
`prospector <https://pypi.org/project/prospector/>`__ and to sort
imports we use `isort <https://pycqa.github.io/isort/>`__. Running
the linters requires an activated virtual environment with the
development tools installed.

.. code:: shell

   # linter
   prospector

   # recursively check import style for the dianna module only
   isort --recursive --check-only dianna

   # recursively check import style for the dianna module only and show
   # any proposed changes as a diff
   isort --recursive --check-only --diff dianna

   # recursively fix import style for the dianna module only
   isort --recursive dianna

You can enable automatic linting with ``prospector`` and ``isort`` on
commit by enabling the git hook from ``.githooks/pre-commit``, like so:

.. code:: shell

   git config --local core.hooksPath .githooks

We also check linting errors in a GitHub Actions CI workflow.

Documentation
-------------

-  Documentation should be put in the ``docs/`` directory in the repository.
-  We use Restructured Text (reST) and Google style docstrings.

   -  `Restructured Text (reST)
      primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__
   -  `Google style docstring
      examples <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`__.

-  The documentation is set up with the ReadTheDocs Sphinx theme.

   -  Check out its `configuration
      options <https://sphinx-rtd-theme.readthedocs.io/en/latest/>`__.

-  `AutoAPI <https://sphinx-autoapi.readthedocs.io/>`__ is used to
   generate documentation for the package Python objects.
-  ``.readthedocs.yaml`` is the ReadTheDocs configuration file. When
   ReadTheDocs is building the documentation this package and its
   development dependencies are installed so the API reference can be
   rendered.
-  `Relevant section in the
   guide <https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=writingdocumentation>`__

Generating documentation
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   cd docs
   make html

The documentation will be in ``docs/_build/html``

If you do not have ``make`` use

.. code:: shell

   sphinx-build -b html docs docs/_build/html

To find undocumented Python objects you can run

.. code:: shell

   cd docs
   make coverage
   cat _build/coverage/python.txt

We also check for undocumented functionality in a GitHub Actions CI workflow.

To `test
snippets <https://www.sphinx-doc.org/en/master/usage/extensions/doctest.html>`__
in documentation run

.. code:: shell

   cd docs
   make doctest

Versioning
----------

Bumping the version across all files is done with
`bumpversion <https://github.com/c4urself/bump2version>`__, e.g.

.. code:: shell

   bumpversion major
   bumpversion minor
   bumpversion patch

Making a release
----------------

This section describes how to make a release in 4 steps:

1. Verify that the information in ``CITATION.cff`` is correct.
2. Make sure the `version has been updated <#versioning>`__.
3. Run the unit tests with ``pytest -v`` or ``tox``.
4. *If applicable:* list non-Python files that should be included in the distribution in ``MANIFEST.in``.
5. Make a `release on GitHub <https://github.com/dianna-ai/dianna/releases/new>`__.
   This will trigger the release workflow, which will build and upload DIANNA as a package to PyPI.
   It will also trigger Zenodo into making a snapshot of the repository and sticking a DOI on it.

Note that the build is uploaded to both pypi and test-pypi.
If you trigger the workflow manually, it's only uploaded to test-pypi, which can be useful for testing... dianna documentation master file, created by
   sphinx-quickstart on Wed May  5 22:45:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dianna's documentation!
==========================================================

.. toctree::
  :maxdepth: 2
  :caption: Contents:

  developer_info.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
