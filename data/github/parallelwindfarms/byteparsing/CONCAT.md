# Byteparsing

![Python package](https://github.com/parallelwindfarms/byteparsing/workflows/Python%20package/badge.svg)
[![PyPI version](https://img.shields.io/pypi/v/byteparsing.svg?colorB=blue)](https://pypi.python.org/pypi/byteparsing/)
[![codecov](https://codecov.io/gh/parallelwindfarms/byteparsing/graph/badge.svg)](https://codecov.io/gh/parallelwindfarms/byteparsing)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

Parser for mixed ASCII/binary data. Features:

- Works extremely well with memory-mapped Numpy arrays
- Included parsers:
    - OpenFOAM

The project setup is documented in [a separate
document](project_setup.rst).
See also the [extended tutorial](https://parallelwindfarms.github.io/byteparsing/functional.html).

## Installation

### With pip

To install the latest release of byteparsing, do:

```{.console}
pip install byteparsing
```

### With GitHub

To install the latest version of byteparsing, do:

```{.console}
git clone https://github.com/parallelwindfarms/byteparsing.git
cd byteparsing
pip install .
```

Run tests (including coverage) with:

``` {.console}
python setup.py test
```

### Contributing

If you want to contribute to the development of byteparsing, have a look
at the [contribution guidelines](CONTRIBUTING.rst).

### License

Copyright (c) 2019, Netherlands eScience Center, University of Groningen

Licensed under the Apache License, Version 2.0 (the \"License\"); you
may not use this file except in compliance with the License. You may
obtain a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an \"AS IS\" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

### Credits

This package was created with
[Cookiecutter](https://github.com/audreyr/cookiecutter) and the
[NLeSC/python-template](https://github.com/NLeSC/python-template).
# Architecture

## Functional parsers
### Trampoline

## Cursors

## Memory mapping


<!-- vim: set ft=markdown: -->
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

#. use the search functionality `here <https://github.com/parallelwindfarms/byteparsing/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/parallelwindfarms/byteparsing/issues>`__ to see if someone already filed the same issue;
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
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the byteparsing repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###########
Change Log
###########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

[Unreleased]
************

Added
-----

* Empty Python project directory structure
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
reported by contacting the project team at j.hidding@esciencecenter.nl. All
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
Project Setup
*************

Here we provide some details about the project setup. Most of the choices are explained in the `guide <https://guide.esciencecenter.nl>`_. Links to the relevant sections are included below.
Feel free to remove this text when the development of the software package takes off.

For a quick reference on software development, we refer to `the software guide checklist <https://guide.esciencecenter.nl/best_practices/checklist.html>`_.

Version control
---------------

Once your Python package is created, put it under
`version control <https://guide.esciencecenter.nl/best_practices/version_control.html>`_!
We recommend using `git <http://git-scm.com/>`_ and `github <https://github.com/>`_.

.. code-block:: console

  cd byteparsing
  git init
  git add -A
  git commit

To put your code on github, follow `this tutorial <https://help.github.com/articles/adding-an-existing-project-to-github-using-the-command-line/>`_.

Python versions
---------------

This repository is set up with Python versions:

* 3.6
* 3.7

Add or remove Python versions based on project requirements. See `the guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html>`_ for more information about Python versions.

Package management and dependencies
-----------------------------------

You can use either `pip` or `conda` for installing dependencies and package management. This repository does not force you to use one or the other, as project requirements differ. For advice on what to use, please check `the relevant section of the guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html#dependencies-and-package-management>`_.

* Dependencies should be added to `setup.py` in the `install_requires` list.

Packaging/One command install
-----------------------------

You can distribute your code using pipy or conda. Again, the project template does not enforce the use of either one. `The guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html#building-and-packaging-code>`_ can help you decide which tool to use for packaging.

Testing and code coverage
-------------------------

* Tests should be put in the ``tests`` folder.
* The ``tests`` folder contains:

  - Example tests that you should replace with your own meaningful tests (file: ``test_byteparsing``)
  - A test that checks whether your code conforms to the Python style guide (PEP 8) (file: ``test_lint.py``)

* The testing framework used is `PyTest <https://pytest.org>`_

  - `PyTest introduction <http://pythontesting.net/framework/pytest/pytest-introduction/>`_

* Tests can be run with ``python setup.py test``

  - This is configured in ``setup.py`` and ``setup.cfg``

* Use `GitHub actions <https://github.com/parallelwindfarms/byteparsing/actions>`_ to automatically run tests and to test using multiple Python versions

  - Configuration can be found in ``.github/workflows/python-package.yml``

* TODO: add something about code quality/coverage tool?
* `Relevant section in the guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html#testing>`_

Documentation
-------------

* Documentation should be put in the ``docs`` folder. The contents have been generated using ``sphinx-quickstart`` (Sphinx version 1.6.5).
* We recommend writing the documentation using Restructured Text (reST) and Google style docstrings.

  - `Restructured Text (reST) and Sphinx CheatSheet <http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html>`_
  - `Google style docstring examples <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_.

* The documentation is set up with the Read the Docs Sphinx Theme.

  - Check out the `configuration options <https://sphinx-rtd-theme.readthedocs.io/en/latest/>`_.

* To generate html documentation run ``python setup.py build_sphinx``

  - This is configured in ``setup.cfg``
  - Alternatively, run ``make html`` in the ``docs`` folder.

* The ``docs/_templates`` directory contains an (empty) ``.gitignore`` file, to be able to add it to the repository. This file can be safely removed (or you can just leave it there).
* To put the documentation on `Read the Docs <https://readthedocs.org>`_, log in to your Read the Docs account, and import the repository (under 'My Projects').

  - Include the link to the documentation in this README.

* `Relevant section in the guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html#writingdocumentation>`_

Coding style conventions and code quality
-----------------------------------------

* Check your code style with ``prospector``
* You may need run ``pip install .[dev]`` first, to install the required dependencies
* You can use ``yapf`` to fix the readability of your code style and ``isort`` to format and group your imports
* `Relevant section in the guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html#coding-style-conventions>`_

Package version number
----------------------

* We recommend using `semantic versioning <https://guide.esciencecenter.nl/best_practices/releases.html#semantic-versioning>`_.
* For convenience, the package version is stored in a single place: ``byteparsing/__version__.py``. For updating the version number, you only have to change this file.
* Don't forget to update the version number before `making a release <https://guide.esciencecenter.nl/best_practices/releases.html>`_!


Logging
-------

* We recommend using the `logging` module for getting useful information from your module (instead of using `print`).
* The project is set up with a logging example.
* `Relevant section in the guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html#logging>`_

CHANGELOG.rst
-------------

* Document changes to your software package
* `Relevant section in the guide <https://guide.esciencecenter.nl/software/releases.html#changelogmd>`_

CITATION.cff
------------

* To allow others to cite your software, add a ``CITATION.cff`` file
* It only makes sense to do this once there is something to cite (e.g., a software release with a DOI).
* Follow the `making software citable <https://guide.esciencecenter.nl/citable_software/making_software_citable.html>`_ section in the guide.

CODE_OF_CONDUCT.rst
-------------------

* Information about how to behave professionally
* `Relevant section in the guide <https://guide.esciencecenter.nl/software/documentation.html#code-of-conduct>`_

CONTRIBUTING.rst
----------------

* Information about how to contribute to this software package
* `Relevant section in the guide <https://guide.esciencecenter.nl/software/documentation.html#contribution-guidelines>`_

MANIFEST.in
-----------

* List non-Python files that should be included in a source distribution
* `Relevant section in the guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html#building-and-packaging-code>`_

NOTICE
------

* List of attributions of this project and Apache-license dependencies
* `Relevant section in the guide <https://guide.esciencecenter.nl/best_practices/licensing.html#notice>`_
.. byteparsing documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to byteparsing's documentation!
==========================================================

Byteparsing is a package for **mixed text and binary** parsing in Python. The main driver for developing this package was to write a parser for binary OpenFOAM files. The binary file format in OpenFOAM is "special". It is the same as the ASCII based text format, except where large blocks of floating point data are concerned.

When **not to use** `byteparsing`:

* You just need to parse some text: use `pyparsing`, it is the industry's standard.

**Do use** `byteparsing` if:

* You need to tinker with large binary OpenFOAM files directly from Python.
* There is a different package that does not adhere to data standards and hacked together its own mixed ASCII/binary file format. You will have to roll out your own parser. Byteparsing can make this easier.

**Coolest feature**:

* Works with :py:mod:`mmap` and :py:mod:`numpy`! This means you can open the file without reading it entirely into memory, change the NumPy array data and the changes are automatically saved to disk.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Introduction <self>
   install
   cursors
   parsers
   functional
   examples
   ply
   architecture
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
API Documentation
=================

.. automodule:: byteparsing.cursor
   :members:

.. automodule:: byteparsing.failure
   :members:

.. automodule:: byteparsing.parsers
   :members:

.. automodule:: byteparsing.trampoline
   :members:

