<img src="https://raw.githubusercontent.com/yaml2sbml-dev/yaml2sbml/main/doc/logo/logo_yaml2sbml_long.svg" alt="yaml2sbml logo"/>

[![CI](https://github.com/yaml2sbml-dev/yaml2sbml/workflows/CI/badge.svg)](https://github.com/yaml2sbml-dev/yaml2sbml/actions)
[![codecov](https://codecov.io/gh/yaml2sbml-dev/yaml2sbml/branch/master/graph/badge.svg)](https://codecov.io/gh/yaml2sbml-dev/yaml2sbml)
[![Documentation Status](https://readthedocs.org/projects/yaml2sbml/badge/?version=latest)](https://yaml2sbml.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/632acdc8d4ef4f50bf69892b8862fd24)](https://www.codacy.com/gh/yaml2sbml-dev/yaml2sbml/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=yaml2sbml-dev/yaml2sbml&amp;utm_campaign=Badge_Grade)
[![PyPI](https://badge.fury.io/py/yaml2sbml.svg)](https://badge.fury.io/py/yaml2sbml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03215/status.svg)](https://doi.org/10.21105/joss.03215)


## Table of contents

* [About](#about)

* [Installation](#installation)

* [Getting Started](#getting-started)

* [Basic Syntax](#basic-syntax)

* [How to cite](#how-to-cite)

* [Contact](#contact)


## About

`yaml2sbml` is a small package in Python to convert an ODE model specified in a YAML file into 
[**SBML**](http://www.sbml.org/) for ODE simulation and into 
[**PEtab**](https://github.com/PEtab-dev/PEtab/) for parameter fitting. `yaml2sbml` offers:


* a translator of ODE models specified in YAML into SBML/PEtab via a Python and a command-line interface;
* a format validator for the input YAML; and
* a model editor, which provides a simplified interface to generate, import and export YAML models.

## Installation

`yaml2sbml` can be installed from PyPI with

```shell
pip install yaml2sbml
```
For more info see the [docs](https://yaml2sbml.readthedocs.io/en/latest/).

## Getting Started

* The documentation of the tool is hosted on [Read the Docs](https://yaml2sbml.readthedocs.io/en/latest/).
* The [format documentation](https://yaml2sbml.readthedocs.io/en/latest/format_specification.html) describes the input YAML. 

* Jupyter notebooks containing examples can be found under [`doc/examples`](https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/examples).  Most notably:
    * [Lotka_Volterra.ipynb](https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/examples/Lotka_Volterra/Lotka_Volterra_python/Lotka_Volterra.ipynb) showing the Python package,
    * [Lotka_Volterra_CLI.ipynb](https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/examples/Lotka_Volterra/Lotka_Volterra_CLI/Lotka_Volterra_CLI.ipynb) showing the command-line interface, and
    * [Lotka_Volterra_Model_Editor.ipynb](https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/examples/Lotka_Volterra/Lotka_Volterra_Model_Editor/Lotka_Volterra_Model_Editor.ipynb) demonstrates the Model Editor.

## Basic Syntax

### Python

A YAML model can be translated to SBML/PEtab in Python via
```python
import yaml2sbml

# SBML conversion
yaml2sbml.yaml2sbml(yaml_dir, sbml_dir)

#PEtab conversion
yaml2sbml.yaml2petab(yaml_dir, 
                     output_dir,
                     sbml_name)
```
### Command-Line Interface

and in the command-line via 
```shell
# SBML conversion
yaml2sbml <yaml_dir> <sbml_dir>

#PEtab conversion
yaml2petab <yaml_dir> <output_dir> <sbml_name>
```

### Format Validation

Format validation is possible in Python via `yaml2sbml.validate_yaml` and in the command-line via `yaml2sbml_validate`.

## How to cite

`yaml2sbml` is published in the [Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.03215). 

When using `yaml2sbml` in your project, please cite

* Vanhoefer J., Matos, M. R. A., Pathirana, D., Schälte, Y. and Hasenauer, J. (2021). yaml2sbml: Human-readable and -writable specification of ODE models and their conversion to SBML. Journal of Open Source Software, 6(61), 3215, https://doi.org/10.21105/joss.03215


```
@article{Vanhoefer2021,
  doi = {10.21105/joss.03215},
  url = {https://doi.org/10.21105/joss.03215},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {61},
  pages = {3215},
  author = {Jakob Vanhoefer and Marta R. A. Matos and Dilan Pathirana and Yannik Schälte and Jan Hasenauer},
  title = {yaml2sbml: Human-readable and -writable specification of ODE models and their conversion to SBML},
  journal = {Journal of Open Source Software}
}
```
## Contact
If you have a question regarding the tool: Please drop us an [issue](https://github.com/yaml2sbml-dev/yaml2sbml/issues/new) or a [mail](mailto:jakob.vanhoefer@uni-bonn.de), we are happy to help.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**What are the steps to reproduce the bug?**
A clear and concise description of how to reproduce the bug.

**Which output did you expect?**
A clear and concise description of what you expected to happen.

**What happened instead?**
A clear and concise description of what you observed instead.

**What could be a solution?**
If possible, a clear and concise description of a possible solution to the problem.

**Environment**
 - Operating system, e.g. Ubuntu 20.04, macOS 11.2.1  
 - `yaml2sbml` version (`pip show yaml2sbml`), e.g. 0.2.1
 - Python version (`python3 -V`)
 - If relevant: versions of other tools involved (e.g. `libsbml)

**Additional context**
Add any other context about the problem here.
# yaml2sbml examples

The [`yaml2sbml`](https://github.com/yaml2sbml-dev/yaml2sbml) package translates ODEs specified in YAML into [SBML](http://sbml.org/) for model simulation and [PEtab](https://github.com/PEtab-dev/PEtab) for parameter fitting. This folder contains several example notebooks for the usage of `yaml2sbml`.

## Scope of the Notebooks

*   _Lotka Volterra Notebooks_
    *   [Lotka Volterra (Python)](./Lotka_Volterra/Lotka_Volterra_python/Lotka_Volterra.ipynb)
        *   Introduces the input format, syntax & capabilities of `yaml2sbml`
        *   Showcases ODE simulation & parameter fitting using [AMICI](https://github.com/AMICI-dev/AMICI) & [pyPESTO](https://github.com/ICB-DCM/pyPESTO)
    *   [Lotka Volterra CLI](./Lotka_Volterra/Lotka_Volterra_CLI/Lotka_Volterra_CLI.ipynb)
        *   Shows the command line interface of `yaml2sbml`
    *   [Lotka Volterra Model Editor](./Lotka_Volterra/Lotka_Volterra_Model_Editor/Lotka_Volterra_Model_Editor.ipynb)
        *   Shows model construction using the model editor of `yaml2sbml`
*   [Sorensen](./Sorensen/yaml2sbml_Sorensen.ipynb)
    *   Application example, model of glucose and insulin metabolism
    *   Shows the extension of an existing YAML model by the model editor
*   [Finite State Projection](./Finite_State_Projection/Finite_State_Projection.ipynb)
    *   Application example of a stochastic gene transcription model having hundreds of states.
    *   Generates a complex and realistic model within a few lines of Python.

*   [Format Features](./Format_Features/Format_Features.ipynb)
    *   Several didactic examples, that show individual features of `yaml2sbml`:
        *   Time-dependent right hand sides
        *   Step functions in the right hand side
        *   Function definitions
                
## Running the examples

The notebooks come with additional dependencies, mainly for running the ODE-simulation and parameter fitting. Information on the installation of AMICI is given in its [installation guide](http://sbml.org/Special/Software/libSBML/docs/formatted/python-api/libsbml-math.html). Further dependencies ( e.g. [pyPESTO](https://github.com/ICB-DCM/pyPESTO)) can be installed via `pip install yaml2sbml[examples]`.
Contribute
==========

Documentation
-------------

To make yaml2sbml easily usable, we try to provide good documentation,
including code annotation and usage examples.
The documentation is hosted on
`yaml2sbml.readthedocs.io <https://yaml2sbml.readthedocs.io>`_
and updated automatically every time the ``main`` branch is updated.
To create the documentation locally, first install the requirements via::

    pip install .[doc]

and then compile the documentation via::

    cd doc
    make html

Test environment
----------------

We use the virtual testing tool `tox <https://tox.readthedocs.io/en/latest/>`_
for all unit tests, format and quality checks and notebooks.
Its configuration is specified in ``tox.ini``. To run it locally, first
install::

    pip install tox

and then simply execute::

    tox

To run only selected tests (see ``tox.ini`` for what is tested), use e.g.::

    tox -e pyroma,flake8

For continuous integration testing we use GitHub Actions. All tests are run
there on pull requests and required to pass. The configuration is specified
in `.github/workflows/ci.yml`.

Unit tests
----------

Unit tests are located in the ``tests`` folder. All files starting with
``test_`` contain tests and are automatically run on GitHub Actions.
Run them locally via::

    tox -e unittests

which boils mostly down to::

    python3 -m pytest tests

You can also run only specific tests.

Unit tests can be written with `pytest <https://docs.pytest.org/en/latest/>`_
or `unittest <https://docs.python.org/3/library/unittest.html>`_.

PEP8
----

We try to respect the `PEP8 <https://www.python.org/dev/peps/pep-0008>`_
coding standards. We run `flake8 <https://flake8.pycqa.org>`_ as part of the
tests. The flake8 plugins used are specified in ``tox.ini`` and the flake8
configuration given in ``.flake8``.
You can run it locally via::

    tox -e flake8

Workflow
--------

If you start working on a new feature or a fix, please create an issue on
GitHub shortly describing the issue and assign yourself.

To get your code merged, please:

1. create a pull request to develop
2. if not already done in a commit message, use the pull request
   description to reference and automatically close the respective issue
   (see https://help.github.com/articles/closing-issues-using-keywords/)
3. check that all tests pass
4. check that the documentation is up-to-date
5. request a code review
Installation
============

This package requires Python 3.6 or later. Miniconda_ provides a small
installation.

Install from PyPI
-----------------

Install yaml2sbml from PyPI_ via::

    pip install yaml2sbml

Install from GitHub
-------------------

To work with the latest development version, install yaml2sbml from
GitHub_ via::

    pip install https://github.com/yaml2sbml-dev/yaml2sbml/archive/develop.zip

or clone the repository and install from local via::

    git clone https://github.com/yaml2sbml-dev/yaml2sbml
    cd yaml2sbml
    git checkout develop
    pip install -e .

where ``-e`` is short for ``--editable`` and links the installed package to
the current location, such that changes there take immediate effect.

Additional dependencies for running the examples
------------------------------------------------

The notebooks come with additional dependencies. Information on the
installation of the ODE simulator `AMICI <https://github.com/AMICI-dev/AMICI>`_ is given in its
`installation guide <https://github.com/AMICI-dev/AMICI/blob/master/INSTALL.md>`_.
Further dependencies can be installed via::

    pip install yaml2sbml[examples]

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PyPI: https://pypi.org/project/yaml2sbml
.. _GitHub: https://github.com/yaml2sbml-dev/yaml2sbml
Release Notes
=============


0.2 series
..........

0.2.5 (2021-07-30)
------------------
* remove unit assignment `dimensionless` from SBML creation. (#136)
* remove checks for consistency of units. (#135)

0.2.4 (2021-05-28)
------------------

* JOSS badge in the README (#133)
* Citation section in the README (#133)

0.2.3 (2021-05-25)
------------------

* Add additional names to the license (#128).
* Add information about the installation of dependencies for the example notebooks (#128).

0.2.2 (2021-03-19)
------------------

* Add checks in SBML conversion, e.g. to catch invalid identifiers and equations (# 118, #123).
* Add tests for MacOS and Windows (#115).
* Fix pydocstyle and swig (#120).
* Add a logo license. (#116)
* Smaller fixes in notebooks (#113, #117).

0.2.1 (2021-02-22)
------------------

* New option to write observables as assignments in `yaml2sbml` (#105).
* Set SBML file name as model id in the SBML (#105).
* Clarify docs and warnings, e.g. if formula starts with a minus (#109).
* Add issue template (#96).
* Restructure README (#55, #103, #104).
* Testing via tox (#88, #94, #95).
* Replace `requirements.txt` and `setup.py` by `setup.cfg` (#87).

0.2.0 (2021-02-03)
------------------

* Initial release on `pypi`.
Deploy
======

New production releases should be created every time the ``main`` branch is
updated.

Versions
--------

On every merge to ``main``, the version number in ``yaml2sbml/version.py`` should
be incremented. We use version numbers ``A.B.C``, where roughly

* ``C`` is increased for bug fixes,
* ``B`` is increased for new features and minor API breaking changes,
* ``A`` is increased for major API breaking changes.

as suggested by the
`Python packaging guide <https://packaging.python.org>`_.

Create a new release
--------------------

After new commits have been added via pull requests to the ``develop`` branch,
changes can be merged to ``main`` and a new version of yaml2sbml can be
released.

Merge into main
~~~~~~~~~~~~~~~

1. create a pull request from ``develop`` to ``main``,
2. check that all tests pass,
3. check that the documentation is up-to-date,
4. update the version number in ``yaml2sbml/version.py`` (see above),
5. update the release notes in ``doc/release_notes.rst``,
6. request a code review.

To be able to actually perform the merge, sufficient rights may be required.
Also, at least one review is required.

Create a release on GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~

After merging into ``main``, create a new release on GitHub. This can be done
either directly on the project GitHub website, or via the CLI as described
in
`Git Basics - Tagging <https://git-scm.com/book/en/v2/Git-Basics-Tagging>`_.
In the release form,

* specify a tag with the new version as specified in ``yaml2sbml/version.py``,
* include the latest additions to ``doc/release_notes.rst`` in the release
  description.

Upload to PyPI
--------------

The upload to the Python package index PyPI has been automated via GitHub
Actions, specified in ``.github/workflows/deploy.yml``,
and is triggered whenever a new release tag is created.

To manually upload a new version to PyPI, proceed as follows:
First, create a so-called "wheel" via::

    python setup.py sdist bdist_wheel

A wheel is essentially a ZIP archive which contains the source code
and the binaries (if any).

Then upload the archive::

    twine upload dist/yaml2sbml-x.y.z-py3-non-any.wheel

replacing x.y.z by the respective version number.

See also the
`section on distributing packages 
<https://packaging.python.org/tutorials/distributing-packages>`_
of the Python packaging guide.
Contact
=======

Feel free to submit an
`issue <https://github.com/yaml2sbml-dev/yaml2sbml/issues>`_
or send us an `e-mail <mailto:marta.ra.matos@gmail.com;jakob.vanhoefer@uni-bonn.de>`_.
Contribute
==========

Documentation
-------------

To make yaml2sbml easily usable, we try to provide good documentation,
including code annotation and usage examples.
The documentation is hosted on
`yaml2sbml.readthedocs.io <https://yaml2sbml.readthedocs.io>`_
and updated automatically every time the ``main`` branch is updated.
To create the documentation locally, first install the requirements via::

    pip install .[doc]

and then compile the documentation via::

    cd doc
    make html

Test environment
----------------

We use the virtual testing tool `tox <https://tox.readthedocs.io/en/latest/>`_
for all unit tests, format and quality checks and notebooks.
Its configuration is specified in ``tox.ini``. To run it locally, first
install::

    pip install tox

and then simply execute::

    tox

To run only selected tests (see ``tox.ini`` for what is tested), use e.g.::

    tox -e pyroma,flake8

For continuous integration testing we use GitHub Actions. All tests are run
there on pull requests and required to pass. The configuration is specified
in `.github/workflows/ci.yml`.

Unit tests
----------

Unit tests are located in the ``tests`` folder. All files starting with
``test_`` contain tests and are automatically run on GitHub Actions.
Run them locally via::

    tox -e unittests

which boils mostly down to::

    python3 -m pytest tests

You can also run only specific tests.

Unit tests can be written with `pytest <https://docs.pytest.org/en/latest/>`_
or `unittest <https://docs.python.org/3/library/unittest.html>`_.

PEP8
----

We try to respect the `PEP8 <https://www.python.org/dev/peps/pep-0008>`_
coding standards. We run `flake8 <https://flake8.pycqa.org>`_ as part of the
tests. The flake8 plugins used are specified in ``tox.ini`` and the flake8
configuration given in ``.flake8``.
You can run it locally via::

    tox -e flake8

Workflow
--------

If you start working on a new feature or a fix, please create an issue on
GitHub shortly describing the issue and assign yourself.

To get your code merged, please:

1. create a pull request to develop
2. if not already done in a commit message, use the pull request
   description to reference and automatically close the respective issue
   (see https://help.github.com/articles/closing-issues-using-keywords/)
3. check that all tests pass
4. check that the documentation is up-to-date
5. request a code review
API documentation
===================


Convert YAML to SBML
---------------------

.. automodule:: yaml2sbml.yaml2sbml
    :members: 


Validate YAML file
---------------------

.. autofunction:: yaml2sbml.validate_yaml


Convert YAML to PEtab
----------------------------------

.. autofunction:: yaml2sbml.yaml2petab


Validate PEtab tables
----------------------------------

.. autofunction:: yaml2sbml.validate_petab_tables


Model editor
----------------------------------
.. autoclass:: yaml2sbml.YamlModel
    :members:
Input Format for yaml2sbml
==========================


General scope
-------------

*  `yaml2sbml`: translate ODEs (initial value problems) of the form `x' = f(t, x, p)` with time `t`, states `x` and (potentially) unknown parameters `p` into an SBML file for simulation purpose.

*  `yaml2PEtab`: define a fitting problem of the form `y(t_i) = h(x(t_i), p) + eps_i` with independent normal- or Laplace-distributed error terms `eps`. `h` denotes the mapping from system states to observables. PEtab allows one to formulate MLE and MAP based fitting problems.

General remarks
---------------

* All identifiers of states, parameters etc. need to be valid SBML identifiers. Therefore, identifiers must consist of only upper and lower case letters, digits and underscores, and must not start with a digit.
* Mathematical equations are parsed by `libsbml`s `parseL3Formula`. Hence for correct syntax see its `documentation <http://sbml.org/Special/Software/libSBML/docs/formatted/python-api/namespacelibsbml.html#ae79acc3be958963c55f1d03944add36b>`_ and the corresponding section of the format specification.
* Equations starting with a minus must be surrounded by brackets or quotation marks, since a leading minus also has a syntactic meaning in YAML and the YAML file will not be valid otherwise.

time \[optional\]
-----------------


.. code-block:: yaml

  time:
    variable: t


Define the **time variable**, in case the right hand side of the ODE is time-dependent.
  

  
parameters \[optional\]
-----------------------

.. code-block:: yaml

  parameters: 
    - parameterId: p_1
      nominalValue: 1
    
    - parameterId: p_2
      ...     


Define **parameters**. `nominalValue` is optional for SBML/PEtab generation, but will be needed for model simulation. Further optional entries are `parameterName, parameterScale, lowerBound, upperBound, estimate` and entries regarding priors. These entries will be written in the corresponding column of the _parameter table_ by `yaml2PEtab.`

For a detailed description see the documentation of the `PEtab parameter table <https://github.com/PEtab-dev/PEtab/blob/master/doc/documentation_data_format.rst#parameter-table>`_ "PEtab parameter table documentation"). 

Further entries are possible and will be written to the _parameter table_ as well but are currently not part of the PEtab standard. 



odes
----

.. code-block:: yaml

  odes:
      - stateId: x_1
        rightHandSide: p_1 * x_1
        initialValue: 1

      - stateId: x_2
        ...      


Define **ODEs** (and states). An ODE consists of a `stateId` (string), a `rightHandSide` (string, encoding a mathematical expression), and an `initial value`. Initial values can be either numerical values or parameter Ids.

For a more detailed description of the parsing of mathematical expressions ( for  `rightHandSide`) we refer to the :ref:`corresponding section<Parsing of mathematical equations>` of this documentation.



assignments \[optional\]
------------------------

.. code-block:: yaml

  assignments:
      - assignmentId: sum_of_states
        formula: x_1 + x_2

      - assignmentId: ...
        ...


**Assign** the mathematical expression `formula` to the term `assignmentId`. The value is dynamically updated and can depend on parameters, states and time. In SBML, assignments are represented via parameter assignment rules.

For a more detailed description of the parsing of mathematical expressions (e.g. for `formula`) we refer to the [corresponding section](#parsing-of-mathematical-equations) of this documentation.



functions \[optional\]
----------------------

.. code-block:: yaml

  functions:
      - functionId: g_1
        arguments: x_1, s
        formula: s * x_1 + 1

      - functionId: g_2
        ...

Define **functions**, which can be called in other parts of the ODE definitions, e.g. in the example above via  `g_1(x_1, s)`.

**Please note** that all unknowns appearing in the formula (e.g. also parameters or the time variable) also have to be arguments of the function.

For a more detailed description of the parsing of mathematical expressions (e.g. for  `formula`) we refer to the [corresponding section](#parsing-of-mathematical-equations) of this documentation.



observables \[optional\]
------------------------

.. code-block:: yaml

  observables:

      - observableId: Obs_1
        observableFormula: x_1 + x_2

        noiseFormula: noiseParameter1
        noiseDistribution: normal

      - observableId: Obs_2
        ...

Define **observables**. Observables are not part of the SBML standard. If the SBML is generated via the `yaml2sbml.yaml2sbml` command and the `observables_as_assignments` flag is set to `True`, observables are represented as assignments to parameters of the form observable_<observable_id>.
If the SBML is created via `yaml2sbml.yaml2petab`, observables are represented in the PEtab observables table. The entries are written to the corresponding columns of the PEtab observable table. According to the PEtab standard, an observable table can take the following entries:  `observableId, observableName, observableFormula, observableTransformation, noiseFormula, noiseDistribution`.

For a detailed discussion see the corresponding part of the PEtab documentation <https://github.com/PEtab-dev/PEtab/blob/master/doc/documentation_data_format.rst#observables-table>`_.



conditions \[optional\]
----------------------

.. code-block:: yaml

  conditions:

      - conditionId: condition1
        p_1: 1
        x_1: 2
        ...


Conditions allow one to set parameters or initial conditions of states to a numeric value/unknown parameter. This allows for the specification of different experimental setups in the data generation (e.g. different initial conditions for different runs of an experiment).

The "trivial condition table" (if only one setup exists) is generated by:

.. code-block:: yaml

  conditions:
        - conditionId: condition1


For a detailed discussion see the corresponding part of the `PEtab documentation <https://github.com/PEtab-dev/PEtab/blob/master/doc/documentation_data_format.rst#condition-table>`_.



Parsing of mathematical equations
---------------------------------

 Throughout `yaml2sbml` formulas are parsed by `libsbml's` `parseL3Formula` function. Further information on the syntax are given by:

*  the `working with math <http://sbml.org/Special/Software/libSBML/docs/formatted/python-api/libsbml-math.html>`_ - section of the `libsbml` documentation.
*  the `documentation <http://sbml.org/Special/Software/libSBML/docs/formatted/python-api/namespacelibsbml.html#ae79acc3be958963c55f1d03944add36b>`_ of `libsbml.parseL3Formula`.


This gives access to e.g.:

*  `+`, `-`, `*`, `/`, and `power;
*  trigonometric/hyperbolic functions;
*  exponential/logarithmic functions;
*  piecewise defined functions (via `piecewise`); and
*  boolean expressions like "<".
Logo
====

.. image:: logo/logo_yaml2sbml_long.png
   :alt: yaml2sbml logo
   :align: center

yaml2sbml's logo can be found in multiple variants in the
`doc/logo <https://github.com/yaml2sbml-dev/yaml2sbml/tree/main/doc/logo>`_
folder, in svg and png format.
It is made available under a Creative Commons CC0 1.0 Universal (CC0 1.0)
license, with the terms given in `doc/logo/LICENSE.md <https://github.com/yaml2sbml-dev/yaml2sbml/blob/master/doc/logo/LICENSE.md>`_.
We encourage to use it e.g. in presentations and posters.

We thank Elba Raimúndez for her contributions to the logo design.
yaml2sbml
=========

.. image:: https://github.com/yaml2sbml-dev/yaml2sbml/workflows/CI/badge.svg
   :target: https://github.com/yaml2sbml-dev/yaml2sbml/actions
   :alt: Build status
.. image:: https://codecov.io/gh/yaml2sbml-dev/yaml2sbml/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/yaml2sbml-dev/yaml2sbml
   :alt: Code coverage
.. image:: https://readthedocs.org/projects/yaml2sbml/badge/?version=latest
   :target: https://yaml2sbml.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://app.codacy.com/project/badge/Grade/632acdc8d4ef4f50bf69892b8862fd24
   :target: https://www.codacy.com/gh/yaml2sbml-dev/yaml2sbml/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=yaml2sbml-dev/yaml2sbml&amp;utm_campaign=Badge_Grade
   :alt: Code quality
.. image:: https://badge.fury.io/py/yaml2sbml.svg
   :target: https://badge.fury.io/py/yaml2sbml
   :alt: PyPI

:Release: |version|
:Source code: https://github.com/yaml2sbml-dev/yaml2sbml

.. image:: logo/logo_yaml2sbml_long.png
   :alt: yaml2sbml logo
   :align: center

**yaml2sbml** allows the user to convert a system of ODEs specified in a YAML_
file to SBML_.
In addition, if experimental data is provided in the YAML file, it can also be
converted to PEtab_.

.. toctree::
   :caption: Main
   :maxdepth: 1
   :hidden:

   install
   format_specification
   examples/examples
   api_doc

.. toctree::
   :caption: About
   :maxdepth: 1
   :hidden:

   release_notes
   license
   log
   logo
   contact

.. toctree::
   :caption: Developers
   :maxdepth: 1
   :hidden:

   contribute
   deploy

.. _YAML: https://yaml.org
.. _SBML: http://sbml.org
.. _PEtab: https://petab.readthedocs.io/en/stable/
Installation
============

This package requires Python 3.6 or later. Miniconda_ provides a small
installation.

Install from PyPI
-----------------

Install yaml2sbml from PyPI_ via::

    pip install yaml2sbml

Install from GitHub
-------------------

To work with the latest development version, install yaml2sbml from
GitHub_ via::

    pip install https://github.com/yaml2sbml-dev/yaml2sbml/archive/develop.zip

or clone the repository and install from local via::

    git clone https://github.com/yaml2sbml-dev/yaml2sbml
    cd yaml2sbml
    git checkout develop
    pip install -e .

where ``-e`` is short for ``--editable`` and links the installed package to
the current location, such that changes there take immediate effect.

Additional dependencies for running the examples
------------------------------------------------

The notebooks come with additional dependencies. Information on the
installation of the ODE simulator `AMICI <https://github.com/AMICI-dev/AMICI>`_ is given in its
`installation guide <https://github.com/AMICI-dev/AMICI/blob/master/INSTALL.md>`_.
Further dependencies can be installed via::

    pip install yaml2sbml[examples]

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PyPI: https://pypi.org/project/yaml2sbml
.. _GitHub: https://github.com/yaml2sbml-dev/yaml2sbml
Examples
===============================================

We provide the following Jupyter notebooks with examples:

* Three notebooks use the Lotka Volterra equations as a running example for the Python tool, CLI and the Model Editor.
* A notebook demonstrating format features with minimal examples.
* A notebook showing the Sorensen model as a real world application.
* A notebook implementing the Finite State Projection.

These examples can also be found on Github `here <https://github.com/yaml2sbml-dev/yaml2sbml_examples>`_.

 .. toctree::
   :maxdepth: 1

   Lotka_Volterra/Lotka_Volterra_python/Lotka_Volterra
   Lotka_Volterra/Lotka_Volterra_CLI/Lotka_Volterra_CLI
   Lotka_Volterra/Lotka_Volterra_Model_Editor/Lotka_Volterra_Model_Editor
   Format_Features/Format_Features
   Sorensen/yaml2sbml_Sorensen
   Finite_State_Projection/Finite_State_Projection







