[![Build Status](https://img.shields.io/github/workflow/status/fitbenchmarking/fitbenchmarking/Build%20and%20Publish?style=flat-square)](https://github.com/fitbenchmarking/fitbenchmarking/actions/workflows/release.yml)
[![Tests Status](https://img.shields.io/github/workflow/status/fitbenchmarking/fitbenchmarking/Tests?label=tests&style=flat-square)](https://github.com/fitbenchmarking/fitbenchmarking/actions/workflows/main.yml)
[![Documentation Status](https://img.shields.io/readthedocs/fitbenchmarking?style=flat-square)](https://fitbenchmarking.readthedocs.io/en/latest)
[![Coverage Status](https://img.shields.io/coveralls/github/fitbenchmarking/fitbenchmarking.svg?style=flat-square)](https://coveralls.io/github/fitbenchmarking/fitbenchmarking)
![Windows Supported](https://img.shields.io/badge/win10-support-blue.svg?style=flat-square&logo=windows)
![Ubuntu Supported](https://img.shields.io/badge/18.04-support-orange.svg?style=flat-square&logo=ubuntu)
[![Chat](https://img.shields.io/badge/chat-CompareFitMinimizers-lightgrey.svg?style=flat-square&logo=slack)](https://slack.com/)
# FitBenchmarking

FitBenchmarking is an open source tool for comparing different minimizers/fitting frameworks. FitBenchmarking is cross platform and we support Windows, Linux and Mac OS. For questions, feature requests or any other inquiries, please open an issue on GitHub, or send us an e-mail at support@fitbenchmarking.com.

- **Installation Instructions:** https://fitbenchmarking.readthedocs.io/en/latest/users/install_instructions/index.html
- **User Documentation & Example Usage:** https://fitbenchmarking.readthedocs.io/en/latest/users/index.html
- **Community Guidelines:** https://fitbenchmarking.readthedocs.io/en/latest/contributors/guidelines.html
- **Automated Tests:** Run via GitHub Actions, https://github.com/fitbenchmarking/fitbenchmarking/actions, and tests are documented at https://fitbenchmarking.readthedocs.io/en/latest/users/tests.html

The package is the result of a collaboration between STFC’s Scientific Computing Department and ISIS Neutron and Muon Facility and the Diamond Light Source. We also would like to acknowledge support from:

* EU SINE2020 WP-10, which received funding from the European Union’s Horizon2020 research and innovation programme under grant agreement No 654000.
* EPSRC Grant EP/M025179/1  Least Squares: Fit for the Future.
* The Ada Lovelace Centre (ALC). ALC is an integrated, cross-disciplinary data intensive science centre, for better exploitation of research carried out at our large scale National Facilities including the Diamond Light Source (DLS), the ISIS Neutron and Muon Facility, the Central Laser Facility (CLF) and the Culham Centre for Fusion Energy (CCFE).
---
title: '`FitBenchmarking`: an open source `Python` package comparing data fitting software'
tags:
  - Python
  - fitting
  - non-linear least squares
authors:
  - name: Anders Markvardsen
    affiliation: 1
  - name: Tyrone Rees
    affiliation: 1
  - name: Michael Wathen
    affiliation: 1
  - name: Andrew Lister
    affiliation: 1
  - name: Patrick Odagiu
    affiliation: 1
  - name: Atijit Anuchitanukul
    affiliation: 1
  - name: Tom Farmer
    affiliation: 1
  - name: Anthony Lim
    affiliation: 1
  - name: Federico Montesino
    affiliation: 1
  - name: Tim Snow
    affiliation: 2
  - name:  Andrew McCluskey
    affiliation: 2
affiliations:
 - name: Science and Technology Facilities Council, Rutherford Appleton Laboratory, Harwell Campus, Didcot, Oxfordshire, OX11 0QX
   index: 1
 - name: Diamond Light Source Ltd, Diamond House, Harwell Campus, Didcot, Oxfordshire, OX11 0DE
   index: 2
date: October 2020
bibliography: paper.bib
---
# Summary

Fitting a mathematical model to data is a fundamental task across all scientific disciplines. [`FitBenchmarking`](https://fitbenchmarking.com/) has been designed to help:

* Scientists, who want to know the best algorithm for fitting their data to a given model using specific hardware.
* Scientific software developers, who want to identify the best fitting algorithms and implementations. This allows them to recommend a default solver, to see if it is worth adding a new minimizer, and to test their implementation.
* Mathematicians and numerical software developers, who want to understand the types of problems on which current algorithms do not perform well, and to have a route to expose newly developed methods to users.

Representatives of each of these communities have got together to build `FitBenchmarking`. We hope this tool will help foster fruitful interactions and collaborations across the disciplines.

![Benchmarking paradigm: associating fitting problems represented in individual scientific software packages (top cycle) to optimization software packages (bottom cycle), and bringing these closer together. \label{fig:concept}](figures/FitBenchmarkingConcept.png){width=60%}

`FitBenchmarking` is easy to install via `pip` and our [documentation](https://fitbenchmarking.com/) guides users through the installation of some external packages we support. We provide several data sets from a range of applications and adding new data in these formats is as easy as dropping the data into a new folder. The data and fitting packages currently supported are shown in Figure \ref{fig:concept}. A key part of `FitBenchmarking` is the ease of which a user, with a basic knowledge of `Python`, can add new fitting software, data formats and different fitting comparison output metrics.

# State of the field

Fitting data to models is a form of optimization, and 
`CUTEst` [@cutest], and its predecessors, has been the standard tool to
benchmark optimization packages for some time. `CUTEst` can benchmark any problem
written in a custom `SIF` format.  However, only the hooks to run the same problem are
provided, the user must provide their own data analysis.  Tools such as
`Paver` [@paver], part of the COIN-OR initiative, can be used
alongside `CUTEst` (or other tools) for this purpose.
The packages `Olympus` [@olympus] and `Benchopt` [@benchopt] have been recently
developed as benchmarking and analysis frameworks for optimization problems.
`Olympus` is designed for experiment planning and provides analytic benchmark problems,
experimental datasets, and emulated datasets, but could be adapted to be applied to
any optimization (or data-fitting) problem.
`Benchopt`, on the other hand, is currently primarily used to benchmark data fitting
using a range of cost functions.  `Benchopt` ships with a limited number of example data
sets, but it is well documented how to write new benchmarks with custom data and
objective functions.


# Statement of need

While there is some overlap between `FitBenchmarking` and the rest of the field,
what makes our software unique is:

* It is designed to interface directly to the source of data, be that a scientific
  software package or an academic data set.  Our `parser` class can be extended
  to make it clear what a developer needs to do to get data into `FitBenchmarking`.
* While being easy to extend using new software, or new data from currently supported
  packages, `FitBechmarking` ships with open datasets that all can use for testing.
* FitBenchmarking tests implementations of algorithms, not just algorithms.
  A growing number of optimization packages that can be used for data fitting are
  supported, and it is straightforward to extend our `controller` class to add new
  software.
* `FitBenchmarking` performs its own data processing and analysis and, if needed,
  the output generated can be customized for new data sets and/or minimizers.  

As far as we are aware, `FitBenchmarking` is the only package that is designed
specifically to interface directly with optimization  packages and individual
scientific software packages to test different implementations of fitting algorithms.
`FitBenchmarking` originally started as a tool to benchmark fitting algorithms in the data reduction package `Mantid` [@mantid], which is used to process neutron scattering and muon spectroscopy data. `FitBenchmarking` has since been significantly extended to take data and models from other real world applications and data analysis / modelling / treatment packages, such as `SasView` [@sasview] and `CUTEst` [@cutest]. It fits models to the data by using a range of data fitting and nonlinear optimization software packages, and present comparisons through a variety of different metrics. These include comparison tables and performance profile plots.

`FitBenchmarking` compares how different fitting algorithms perform for the same data, model and initial guess. The best parameters for the model are found by solving a nonlinear least-squares problem, which can either be solved using a dedicated optimisation software package or using a fitting algorithm implementation within a scientific software package. Figure \ref{fig:sample} displays a data set from `FitBenchmarking` where the crosses are the data points and the two curves are the fits found by two optimization algorithms implemented in `GSL` [@gsl]. From Figure \ref{fig:sample}, it is clear that the solution given by lmsder is better. As the volume of data increases, and we do more and more scientific analysis algorithmically, it is increasingly important that we apply the best available algorithm for a given category of fitting problems. `FitBenchmarking` generates HTML output that makes it easy to compare minimizers on a given problem set.

![A sample fit: this problem is shipped with `FitBenchmarking`. The data was collected from an instrument named VESUVIO at the ISIS Neutron and Muon Source and has a difficult initial guess. \label{fig:sample}](figures/nmsimplex2_fit_for_EVS14188-90_processed_Gaussian_peaks_1_1.png){width=70%}

`FitBenchmarking` will help the scientist make an informed choice by comparing runtime and accuracy of all available minimizers, on their specific hardware, on problems from their science area.

`FitBenchmarking` will help the scientific software developer ensure that the most robust and quickest algorithms for the type of data analysis they support are available in their software.

`FitBenchmarking` will help mathematicians see what the state of the art is, and what kinds of data are problematic. It will give them access to real data, and will give a route for novel methods to quickly make it into production.

# Acknowledgements

We would like to acknowledge funding support from:

* European Union’s Horizon2020 research and innovation programme, EU SINE2020 WP-10,
* EPSRC Grant EP/M025179/1 -- Least Squares: Fit for the Future.
* The Ada Lovelace Centre (ALC).

We would also like to thank Nick Draper, Roman Tolchenov, Nick Gould and Jaroslav Fowkes for their helpful comments and advice.

# References
#### Description of Work

Fixes


#### Testing Instructions

1.
2.
3.

Function: Does the change do what it's supposed to?

Tests: Does it pass? Is there adequate coverage for new code?

Style: Is the coding style consistent? Is anything overly confusing?

Documentation: Is there a suitable change to documentation for this change?
---
name: Test
about: Specify a test that needs to be implemented
title: ''
labels: Testing

---

**Which module and class/method/function does this relate to?**
Please provide the module and class, method or function name(s) that the tests will apply to.

**What aspect requires additional tests?**
Please provide a description of the specifics of the test - what is the behaviour that the tests will inspect.  Ideally this will include a description of the setup (e.g. Pytest fixtures) and the expected results.

**Is this a unit, system or functional test?**
Simply state what type of test you are expecting is required.

**Additional context**
Add any other context about the tests here.
---
name: Documentation
about: Suggest a improvement to user, developer or other documentation
title: ''
labels: Documentation

---

**Description of the documentation**
A description of what area can benefit from better documentation.

**Additional details**
Add any other context or screenshots about this request.
---
name: Bug report
about: Report an error which requires fixing
title: ''
labels: Bug

---

**Description of the error**
A clear and concise description of what the problem is, including steps to reproduce it and the environment you are running FitBenchmarking in.

**Describe the expected result**
What is the result you expect from running the steps described above?

**Describe the actual result**
What was the actual result.

**Suggested fix**

**Additional details**
Add any other context or screenshots about the feature request here.
---
name: Feature request
about: Suggest a new feature
title: ''
labels: Enhancement

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
This folder contains problem definition files each defining a
fit benchmarking problem that minimizers can be benchmarked against.

More specifically a folder or sub-folder of this directory contain a problem set.

See https://fitbenchmarking.readthedocs.io and section on Problem Definition Files
for the formats supported.

A data file name specified within a definition file is recommended to be stored in a sub-folder named `data_files` relative to the location of the definition file.

Examples of Problem Definition folders include (please note this list is stadily changes and hence below may get out of sync from time to time):

* SAS_modelling : fitting problems specific relevant to fitting SAS (Small Angle Scattering) data
* Muon : generic fitting problems relevant to fitting Muon data collected at a Muon facility such as ISIS Neutron and Muon Facility
* Neutron : generic fitting problems relevant to fitting Neutron data collected at a Neutron facility such as ISIS Neutron and Muon Facility
* NIST : set of made up fitting problem (not against measured data with error bars) as described [here](https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml)
* CUTEes : fitting problems relevant to this tool included in [CUTEet](http://epubs.stfc.ac.uk/bitstream/9327/RAL-TR-2013-005.pdf)

This folder contains scripts that have been used at various times to create more generic type data formatted files from more software specific type ones,
thereby allowing problem definition files referencing such data to be made easier available for benchmarking across minimizer libraries.
Of course, this is not aways feasible, but where it is it is recommended.
A particular simple format is column ascii format, for 2D data this is X and Y columns with an optional E column,
where the E column contains the errors of the Y values.This folder includes images etc. used in the main FitBenchmarking [README](../README.md). Other documentation can also be found on the FitBenchmarking [Wiki](https://github.com/fitbenchmarking/fitbenchmarking/wiki).

Furthermore, sphinx documentation is being setup within the source folder and which can be build with Makefile and make.bat. Also this docs is linked up with [Read the Docs](https://readthedocs.org/) and the automatically build sphinx is viewable from https://fitbenchmarking.readthedocs.io
acc: Start

The accuracy results are calculated from the final chi squared value:

.. math:: \min_p \sum_{i=1}^n \left( \frac{y_i - f(x_i, p)}{e_i} \right)^2

where :math:`n` data points :math:`(x_i,y_i)`, associated errors :math:`e_i`, and a model function :math:`f(x,p)`.

acc: End

compare: Start

The combined results show the accuracy in the first line of the cell, and the runtime on the second line of the cell measured in seconds (:math:`s`).

compare: End

local_min: Start

The local min results shows a ``True`` or ``False`` value together with :math:`\frac{|| J^T r||}{||r||}`. The ``True`` or ``False`` indicates whether the software finds a minimum with respect to the following criteria:


- :math:`||r|| \leq \mbox{RES\_TOL}`,
- :math:`|| J^T r|| \leq \mbox{GRAD\_TOL}`,
- :math:`\frac{|| J^T r||}{||r||} \leq \mbox{GRAD\_TOL}`,

where :math:`J` and :math:`r` are the Jacobian and residual of :math:`f(x, p)`, respectively. The tolerances can be found in the results object.

local_min: End

runtime: Start

The timing results are calculated from an average using the `timeit <https://docs.python.org/2/library/timeit.html>`_  module in python. The number of runs can be set in :ref:`options`.

runtime: End

abs: Start

Absolute values are displayed in the table.

abs: End

rel: Start

Relative values are displayed in the table.

rel: End

both: Start

Absolute and relative values are displayed in the table in the format: ``abs (rel)``, for example ``5.1 (1)`` where 5.1 is abs and 1 is rel values respectively

both: End
.. _mainindex:

.. FitBenchmarking documentation master file, created by
   sphinx-quickstart on Wed Sep 11 09:17:28 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

###########################################
Welcome to FitBenchmarking's documentation!
###########################################

.. mdinclude:: ../../README.md

Table Of Contents
-----------------
  
.. toctree::
   :maxdepth: 2
   :titlesonly:

   Concept <concept/index>
   Users <users/index>
   Extending <extending/index>
   Contributors <contributors/index>



.. _workflow:

############
Git Workflow
############


======
Issues
======

All new work should start with a
`new GitHub issue <https://github.com/fitbenchmarking/fitbenchmarking/issues/new/choose>`_
being filed.
This should clearly explain what the change to the code will do.
There are templates for *Bug report*, *Documentation*,
*Feature request* and *Test* issues on GitHub, and you can also
open a blank issue if none of these work.

If issues help meet a piece of work agreed with our funders, it
is linked to the appropriate `Milestone <https://github.com/fitbenchmarking/fitbenchmarking/milestones>`_ in GitHub.

===============
Adding new code
===============

The first step in adding new code is to create a branch, where the work
will be done. Branches should be named according to the convention
`<nnn>-description_of_work`, where `<nnn>` is the issue number.

Please ensure our :ref:`guidelines` are adhered to throughout
the branch.

When you think your new code is ready to be merged into the codebase,
you should open a pull request (PR) to master.
The description should contain the
words `Fixes #<nnn>`, where `<nnn>` is the issue number; this will ensure
the issue is closed when the code is merged into master.  At this point
the automated tests will trigger, and you can see if the code passes on
an independent system.

Sometimes it is desirable to open a PR when the code is not
quite ready to be merged.  This is a good idea, for example, if you want
to get an early opinion on a coding decision.  If this is the case, you
should mark the PR as a *draft* on GitHub.

Once the work is ready to be reviewed, you may want to assign a reviewer,
if you think someone would be well suited to review this change.  It is worth
messaging them on, for example, Slack, as well as requesting their review on
GitHub.

================
Release branches
================

Branches named `release-*` are protected branches; code must be approved by
a reviewer before being added to them, and automated tests will be run on
PRs to these branches.  If code is to be included in the release, it
must be pulled into this branch from master.

Release branches should have the format `release-major.minor.x`, starting from
`release-0.1.x`.  When the code is released, we will tag that commit with
a version number `v0.1.0`.  Any hotfixes will increment `x` by one, and a new tag will
be created accordingly.  If at some point we don't want to provide hot-fixes
to a given minor release, then the corresponding release branch may be deleted.

All changes must be initially merged into master.
There is a `backport-candidate` label, which must be put on PRs
that in addition must be merged into the release branch.

The recommended mechanism for merging PRs labeled with `backport-candidate` into
master is to use the `Squash and merge` commit option:

.. figure:: ../../images/squash-and-merge.png
   :alt: Squashing and merging


After such a PR (with label `backport-candidate`) has been merged into master, it
must then subsequently be merged into the release branch as soon as possible.
It is the responsibility of the person merging such PRs to also perform this
merge into the release branch.

This can be done using git cherry
pick:

.. code-block:: bash

   git checkout release-x.x.x
   git cherry-pick -x <commit-id>
   git push

If you didn't do a squash merge, you will have to cherry pick each commit in
the PR that this being backported separately.

If you encounter problems with cherry picking into release branch please
don't hesitate to speak to an experienced member of the FitBenchmarking team.


================
Creating a release
================
In order to create a new release for FitBenchmarking, there are a few manual steps.
These have been streamlined as much as possible.

First checkout the branch to create the release from.  Releases should only be made from a `release-x-x` branch, not a development branch or master.

From the root of the repo run the "ci/prep_and_tag_release.sh" script with the new version number.
The version number will be rejected if it is not of the expected form.
We expect a "v" followed by the major, minor, and patch numbers,
and an optional numbered label to mark the type of release.

Possible labels are:
 - -beta (release for testing)
 - -rc (release candidate)

This script will create a new commit with the docs and testing links updated, tag it,
and revert the change in a second commit so that the links point back to the latest versions.

These commits will need to be pushed to github.

Finally, you will need to create a release on github.
This can be done by navigating to the releases page, selecting new release
and typing in the tag that was given to the release
(it should tell you the tag exists at this point!).

For example, For a first beta version of release 0.1.0, one would run:

.. code-block:: bash

   git checkout release-0.1.x
   ci/prep_and_tag_release.sh v0.1.0-beta1
   git push origin release-0.1.x

   <And make the release on GitHub>

===================
Adding New Datasets
===================

Users or developers are encouraged to add new data sets following
the instructions :ref:`here<adding_data>`.  Someone in SCD's 
Computational Mathematics Group must make this publically available
by:

- adding the `zip` and `tar.gz` archives to `powell:/var/www/html/numerical-www/fitbenchmarking/`

- adding the datasets to the master `examples.zip` and `examples.tar.gz` folders, and updating the versions on `powell`

Please note that the maximum file size allowed by GitHub is 100MB, and the
total repository size is recommended to be kept below 1GB.  Please bear
this in mind when advising users whether or not they should also add their
data to the `examples/benchmark_problems` directory of FitBenchmarking.
.. _install_instructions:

#####################################
Install Instructions for Contributors
#####################################

Please see below for the recommended install instructions
for new contributors:

1. Before installing, it is recommended that you create an empty
   virtual environment. For more information about installing
   packages using virtual environments please see 
   `here <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__.

2. To ensure that you are working with the latest version of the
   code, you should install Fitbenchmarking from source. For
   details on how to do this, please see :ref:`installing_from_source`.

   .. note::
        The `editable` flag should be used when installing from source, like 
        ``pip install -e ...``, so that changes made to the cloned code 
        are immediately reflected in the imported package.

3. So that you are able to run all tests locally, all extra dependencies and
   external software should be installed. Please see :ref:`extra_dependencies`
   and :ref:`external-instructions` for instructions on how to do this.

   .. note::
        Please note that if you are using WSL, Matlab is likely
        to be difficult to install so it is recommended that you use an alternative
        operating system if you are working on Matlab related issues. 

        Additionally, the tests run through Github do not currently include Matlab
        fitting software, so please ensure tests are run locally before submitting
        code.

4. You should check that the versions of packages such as pylint are up to date with
   those used listed in requirements.txt. This can be done by running the command
   ``pip install -r requirements.txt`` from  within the ``fitbenchmarking`` directory... _structure:

********************
Repository Structure
********************

At the root of the repository there are six directories:

 - build
 - ci
 - Docker
 - docs
 - examples
 - fitbenchmarking


#################
Build (``build``)
#################

This directory contains scripts to allow for installing packages such as Mantid
through setuptools.

#########################
CI (``ci``)
#########################

We use `GitHub Actions <https://github.com/fitbenchmarking/fitbenchmarking/actions>`__
to run our Continuous Integration tests.
The specific tests run are defined in a series of Bash scripts,
which are stored in this folder.

###################
Docker (``Docker``)
###################

The continuous integration process on Github Actions currently run on a Docker container,
and this directory holds the Dockerfiles.  The Docker containers are hosted on
Dockerhub.

``BasicInstall`` holds the Dockerfile that is pushed to the repository ``fitbenchmarking/fitbenchmarking-deps``, the latest of which should have the tag ``latest``.  This contains a basic Ubuntu install, with just the minimal infrastructure needed to run the tests.

``FullInstall`` holds the Dockerfile that is pushed to the repository ``fitbenchmarking/fitbenchmarking-extras``, the latest of which should have the tag ``latest``.  This is built on top of the basic container, and includes optional third party software that FitBenchmarking can work with.

The versions on Docker Hub can be updated from a connected account by issuing the commands:

.. code-block:: bash
		
		docker build --tag fitbenchmarking-<type>:<tag>
		docker tag fitbenchmarking-<type>:<tag> fitbenchmarking/fitbenchmarking-<type>:<tag>
		docker push fitbenchmarking/fitbenchmarking-<type>:<tag>

where ``<type>`` is, e.g., ``deps`` or ``extras``, and ``<tag>`` is, e.g., ``latest``.

########################
Documentation (``docs``)
########################

The documentation for FitBenchmarking is stored in this folder under
``source``.
A local copy of the documentation can be build using ``make html`` in the
``build`` directory.


#######################
Examples (``examples``)
#######################

Examples is used to store sample problem files and options files.

A collection of problem files can be found organised into datasets within the
``examples/benchmark_problems/`` directory.

An options template file and a prewritten options file to run a full set of
minimizers is also made available in the ``examples/`` directory.


###################################################
FitBenchmarking Package (:py:mod:`fitbenchmarking`)
###################################################

The main FitBenchmarking package is split across several directories
with the intention that it is easily extensible.
The majority of these directories are source code, with exceptions being
Templates, Mock Problems, and System Tests.

Each file that contains source code will have a directory inside it called
``tests``, which contains all of the tests for that section of the code.


Benchmark Problems (``benchmark_problems``)
===========================================

This is a copy of the NIST benchmark problems from `examples/benchmark_problems`.
These are the default problems that are run if no problem is passed to the
``fitbenchmarking`` command, and is copied here so that it is distributed
with the package when installed using, say, `pip`.



CLI (:py:mod:`~fitbenchmarking.cli`)
====================================

The CLI directory is used to define all of the entry points into the software.
Currently this is just the main `fitbenchmarking` command.


Controllers (:py:mod:`~fitbenchmarking.controllers`)
====================================================

In FitBenchmarking, controllers are used to interface with third party
minimizers.

The controllers directory holds a base controller class
(:class:`~fitbenchmarking.controllers.base_controller.Controller`) and all its subclasses,
each of which of which interfaces with a different fitting package.  The controllers
currently implemented are described in :ref:`fitting_option` and :ref:`minimizer_option`.

New controllers can be added by following the instructions in :ref:`controllers`.


Core (:py:mod:`~fitbenchmarking.core`)
======================================

This directory holds all code central to FitBenchmarking.
For example, this manages calling the correct parser and controller, as well as
compiling the results into a data object.

Hessian (:py:mod:`~fitbenchmarking.Hessian`)
==============================================

This directory holds the :class:`~fitbenchmarking.hessian.base_hessian.Hessian` class,
and subclasses, which are used by the controllers to approximate second derivatives.
Currently available options are described in :ref:`fitting_option`, and new
Hessians can be added by following the instructions in :ref:`hessian_extend`.

Jacobian (:py:mod:`~fitbenchmarking.jacobian`)
==============================================

This directory holds the :class:`~fitbenchmarking.jacobian.base_jacobian.Jacobian` class,
and subclasses, which are used by the controllers to approximate derivatives.
Currently available options are described in :ref:`jacobian_option`, and new
numerical Jacobians can be added by following the instructions in
:ref:`jacobian_extend`.


Mock Problems (``mock_problems``)
=================================

The mock problems are used in some tests where full problem files are required.
These are here so that the examples can be moved without breaking the tests.


Parsing (:py:mod:`~fitbenchmarking.parsing`)
============================================

The parsers read raw data into a format that FitBenchmarking can use.
This directory holds a base parser,
:class:`~fitbenchmarking.parsing.base_parser.Parser` and all its subclasses.
Each subclass implements a parser for a specific file format.
Information about existing parsers can be found in :ref:`problem_def`, and
see :ref:`parsers` for instructions on extending these.


Results Processing (:py:mod:`~fitbenchmarking.results_processing`)
==================================================================

All files that are used to generate output are stored here.
This includes index pages, text/html tables, plots, and support pages.
Information about the tables we provide can be found in
:ref:`output`, and instructions on how to add further tables and change
the formatting of the displayed information can be found in :ref:`extending_outputs`.

System Tests (``systests``)
===========================

FitBenchmarking runs regression tests to check that the
accuracy results do not change with updates to the code.
These tests run fitbenchmarking against a subset of problems
(in subdirectories of `/fitbenchmarking/mock_problems/`),
and compares the text output with that stored in
`/fitbenchmarking/systests/expected_results/`.

Templates (``templates``)
=========================

Files in Templates are used to create the resulting html pages, and are a
combination of css, html, and python files.
The python files in this directory are scripts to update the css and html
assets.
Instructions on updating these can be found in :ref:`templates`.

Utils (:py:mod:`~fitbenchmarking.utils`)
========================================

This directory contains utility functions that do not fit into the
above sections.
This includes the :class:`~fitbenchmarking.utils.options.Options`
class (see :ref:`options_extend` to extend) 
and :class:`~fitbenchmarking.utils.fitbm_result.FittingResult` class,
as well as functions for logging and directory creation.
.. _contributors_index:

#########################################
FitBenchmarking Contributor Documentation
#########################################

Thank you for being a contributor to the FitBenchmarking project.
Here you will find all you need in order to get started.

.. toctree::
    :titlesonly:
    :maxdepth: 2
    :caption: Contents:

    install_instructions
    guidelines
    workflow
    structure
    Module Index <module_index/fitbenchmarking>

.. _guidelines:

################
Coding Standards
################

All code submitted must meet certain standards, outlined below, before
it can be merged into the master branch.  It is the contributor's
job to ensure that the following is satisfied, and the reviewer's
role to check that these guidelines have been followed.

The workflow to be used for submitting new code/issues is described in
:ref:`workflow`.

=======
Linting
=======

All pull requests should be `PEP 8 compliant <https://www.python.org/dev/peps/pep-0008/>`_.
We suggest running code through `flake8 <https://flake8.pycqa.org/en/latest/>`_ and
`pylint <https://www.pylint.org/>`_ before submitting to check for this.


=============
Documentation
=============

Any new code will be accepted only if the documentation, written in
`sphinx <https://www.sphinx-doc.org/en/master/>`_ and found in `docs/`,
has been updated accordingly, and the docstrings in the code
have been updated where necessary.

=======
Testing
=======

All tests should pass before submitting code.
Tests are written using `pytest <https://docs.pytest.org/en/stable/>`_.

The following should be checked before any code is merged:

 - Function: Does the change do what it's supposed to?
 - Tests: Does it pass? Is there adequate coverage for new code?
 - Style: Is the coding style consistent? Is anything overly confusing?
 - Documentation: Is there a suitable change to documentation for this change?

=======
Logging
=======

Code should use the logging in ``utils.log``. This uses Python's built in
`logging module <https://docs.python.org/3.8/library/logging.html>`__,
and should be used in place of any print statements to ensure that persistent
logs are kept after runs.
.. _jacobian_extend:

####################
Adding new Jacobians
####################

*FitBenchmarking allows the use of custom methods for evaluating the
Jacobian of the model.  This section describes how to add further
methods within FitBenchmarking*

In order to add a new Jacobian evaluation method, ``<jac_method>``,
you will need to:

1. Create ``fitbenchmarking/jacobian/<jac_method>_jacobian.py``,
   which contains a new subclass of
   :class:`~fitbenchmarking.jacobian.base_jacobian.Jacobian`.
   Then implement the method:

    -  .. automethod:: fitbenchmarking.jacobian.base_jacobian.Jacobian.eval()
              :noindex:
		 
   The numerical method is set sequentially within
   :meth:`~fitbenchmarking.core.fitting_benchmarking.loop_over_jacobians()` by
   using the ``method`` attribute of the class.

2. Enable the new method as an option in :ref:`fitting_option`,
   following the instructions in :ref:`options_extend`.  Specifically:
   
   * Amend the ``VALID_FITTING`` dictionary so that the element associated
     with the ``jac_method`` key contains the new ``<jac_method>``.
     
   * Extend the ``VALID_JACOBIAN`` dictionary to have a new
     key ``<jac_method>``, with the associated element being a list of
     valid options for this Jacobian.
     
   * Extend the ``DEFAULT_JACOBIAN`` dictionary to have a new key
     ``<jac_method>``, with the associated element being a subset of the
     valid options added in ``VALID_JACOBIAN`` in the previous step.

   * Amend the file ``fitbenchmarking/utils/tests/test_options_jacobian.py`` to
     include tests for the new options.

3. Document the available Jacobians by:

  * adding a list of available ``method`` options to the docs for :ref:`jacobian_option`,
    and including licencing information, if appropriate.
  * updating any example files in the ``examples`` directory

4. Create tests for the Jacobian evaluation in
   ``fitbenchmarking/jacobian/tests/test_jacobians.py``.


The :class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When adding new Jacobian, you will find it helpful to make use of the
following members of the :class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem`
class:

.. currentmodule:: fitbenchmarking.parsing.fitting_problem
.. autoclass:: fitbenchmarking.parsing.fitting_problem.FittingProblem
          :members: eval_model, data_x, data_y, data_e
          :noindex:
.. _adding_data:

################
Adding More Data
################

We encourage users to contribute more data sets to FitBenchmarking;
the more data we have available to test against, the more useful
FitBenchmarking is. First, please ensure that there is a
:ref:`parser<problem_def>` available that can read your dataset, and if
necessary follow the :ref:`instructions to add this
functionality<parsers>` to FitBenchmarking.  Once this is done,
follow the steps below and add to a pull request to make the
data available to others.

1. Create a directory that contains:
   
   - the data sets to be included
     
   - a file `META.txt` containing metadata about the dataset.
     
   - a subfolder `data_files` which contains any supplemental data
   needed by the data parser.  We particularly encourage analytic
   derivative information, if available.
     
2. Update the :ref:`BenchmarkProblems` page to include a description of
   the dataset.  As well as information about the source of the data, this
   should include:

   - information about how many parameters and how many data points
   are to be fitted in the dataset
   
   - details of any external software that needs to be installed to load these
   datasets.
   
3. Create `zip` and `tar.gz` archives of these directories, and pass along
   to one of the core developers to put on the webspace.  They will pass you a
   location of the dataset to include in the description page, and update the
   folder containing all examples to contain your data set.

4. If the data is to be distributed with the GitHub source, add the directory to the
   `examples/benchmark_problems` folder and commit to the repository.  Please note
   that the maximum size of a file supported by GitHub is 100MB, and so datasets
   with files larger than this cannot be added to the repository.  Also, GitHub
   recommends that the size of a Git repository is not more than 1GB.



.. _hessian_extend:

###################
Adding new Hessians
###################

*This section describes how to add further methods to approximate the
Hessian within FitBenchmarking*

In order to add a new Hessian evaluation method, ``<hes_method>``,
you will need to:

1. Create ``fitbenchmarking/hessian/<hes_method>_hessian.py``,
   which contains a new subclass of
   :class:`~fitbenchmarking.hessian.base_hessian.hessian`.
   Then implement the method:

    -  .. automethod:: fitbenchmarking.hessian.base_hessian.Hessian.eval()
              :noindex:

   The method is set sequentially within
   :meth:`~fitbenchmarking.core.fitting_benchmarking.loop_over_hessians()` by
   using the ``method`` attribute of the class.

2. Enable the new method as an option in :ref:`fitting_option`,
   following the instructions in :ref:`options_extend`.  Specifically:
   
   * Amend the ``VALID_FITTING`` dictionary so that the element associated
     with the ``hes_method`` key contains the new ``<hes_method>``.

3. Document the available Hessians by:

  * adding to the list of available ``method`` options under ``hes_method`` in :ref:`fitting_option`.
  * updating any example files in the ``examples`` directory

4. Create tests for the Hessian evaluation in
   ``fitbenchmarking/hessian/tests/test_hessians.py``.


The :class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem` and :class:`~fitbenchmarking.jacobian.base_jacobian.Jacobian`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When adding new Hessian, you will find it helpful to make use of the
following members of the :class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem`
and subclasses of :class:`~fitbenchmarking.jacobian.base_jacobian.Jacobian`:

.. currentmodule:: fitbenchmarking.parsing.fitting_problem
.. autoclass:: fitbenchmarking.parsing.fitting_problem.FittingProblem
          :members: eval_model, data_x, data_y, data_e
          :noindex:

.. currentmodule:: fitbenchmarking.jacobian.base_jacobian
.. autoclass:: fitbenchmarking.jacobian.base_jacobian.Jacobian
          :members: eval
          :noindex:
.. _parsers:

#######################################
Adding Fitting Problem Definition Types
#######################################

The problem definition types we currently support are listed in the page :ref:`problem_def`.

To add a new fitting problem type, the parser name
must be derived from the file to be parsed.
For current file formats by including it as the first line
in the file. e.g ``# Fitbenchmark Problem`` or ``NIST/ITL StRD``, or by checking
the file extension.

To add a new fitting problem definition type, complete the following steps:

1. Give the format a name (``<format_name>``).
   This should be a single word or string of alphanumeric characters,
   and must be unique ignoring case.
2. Create a parser in the ``fitbenchmarking/parsing`` directory.
   This parser must satisfy the following:

   - The filename should be of the form ``"<format_name>_parser.py"``
   - The parser must be a subclass of the base parser, :class:`~fitbenchmarking.parsing.base_parser.Parser`
   - The parser must implement ``parse(self)`` method which takes only ``self``
     and returns a populated :class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem`

   Note: File opening and closing is handled automatically.

3. If the format is unable to accommodate the current convention of
   starting with the ``<format_name>``, you will need to edit
   :class:`~fitbenchmarking.parsing.parser_factory.ParserFactory`.
   This should be done in such a way that the type is inferred from the file.

4. Document the parser (see :ref:`problem_def`), being sure to include any licencing
   information.

5. Create the files to test the new parser.
   Automated tests are run against the parsers in FitBenchmarking,
   which work by using test files in 
   ``fitbenchmarking/parsing/tests/<format_name>``.
   In the :meth:`test_parsers.generate_test_cases` function,
   one needs to add the new parser's
   name to the variable ``formats``,
   based on whether or not the parser is ``pip`` installable.
   There are 2 types of test files needed:

   - **Generic tests**: ``fitbenchmarking/parsing/tests/expected/`` contains
     two files, ``basic.json`` and ``start_end_x.json``.
     You must write two input files in the new file format,
     which will be parsed using the new parser to check that the entries
     in the generated fitting problem match the values expected.
     These must be called ``basic.<ext>``, ``start_end_x.<ext>``, where ``<ext>``
     is the extension of the new file format, and they must be placed in
     ``fitbenchmarking/parsing/tests/<format_name>/``.

   - **Function tests**: A file named ``function_evaluations.json``
     must also be provided in
     ``fitbenchmarking/parsing/tests/<format_name>/``, which tests that the 
     function evaluation behaves as expected. This file must be in json format and
     contain a string of the form::

       {"file_name1": [[[x11,x12,...,x1n], [param11, param12,...,param1m], [result11,result12,...,result1n]],
                       [[x21,x22,...,x2n], [param21, param22,...,param2m], [result21,result22,...,result2n]],
                       ...],
       {"file_name2": [...],
        ...}

     The test will then parse the files ``file_name<x>`` in turn evaluate the function
     at the given ``xx`` values and ``params``. If the result is not suitably close to
     the specified value the test will fail.

   - **Jacobian tests**: *If the parser you add has analytic Jacobian
     information*, then in ``test_parsers.py`` add 
     ``<format_name>`` to the ``JACOBIAN_ENABLED_PARSERS`` global variable.
     Then add a file ``jacobian_evaluations.json`` to
     ``fitbenchmarking/parsing/tests/<format_name>/``, which tests that the Jacobian evaluation behaves as expected.
     This file should have the same file structure as `function_evaluations.json`,
     and works in a similar way. 

   - **Hessian tests**: *If the parser you add has analytic Hessian
     information*, then in ``test_parsers.py`` add 
     ``<format_name>`` to the ``HESSIAN_ENABLED_PARSERS`` global variable.
     Then add a file ``hessian_evaluations.json`` to
     ``fitbenchmarking/parsing/tests/<format_name>/``, which tests that the Hessian evaluation behaves as expected.
     This file should have the same file structure as `function_evaluations.json`,
     and works in a similar way. 

   - **Integration tests**: Add an example to the directory
     ``fitbenchmarking/mock_problems/all_parser_set/``.
     This will be used to verify that the problem can be run by scipy, and that
     accuracy results do not change unexpectedly in future updates.
     If the software used for the new parser is pip-installable, and the
     installation is done via FitBenchmarking's ``setup.py``, then add the
     same example to ``fitbenchmarking/mock_problems/default_parsers/``.

     As part of this, the ``systests/expected_results/all_parsers.txt`` file,
      and if necessary the ``systests/expected_results/default_parsers.txt`` file,
      will need to be updated. This is done by running the systests::

       pytest fitbenchmarking/systests

     and then checking that the only difference between the results table and the
     expected value is the new problem, and updating the expected file with the result.

6. Verify that your tests have been found and are successful by running
   `pytest -vv fitbenchmarking/parsing/tests/test_parsers.py`

Once the new parser is added, please add some examples that use this
problem definition type following the instructions at :ref:`adding_data`.
.. _options_extend:

##################
Adding new Options
##################

Default options are set by the class
:class:`~fitbenchmarking.utils.options.Options`, which is defined in
the file `fitbenchmarking/utils/options.py`.

The options used can be changed using an `.ini` formatted file
(`see here <https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`__). :ref:`options` gives examples of how this is currently implemented in
FitBenchmarking.

To add a new option to one of the five sections ``FITTING``, ``MINIMIZERS``,
``JACOBIAN``, ``PLOTTING`` and ``LOGGING``, follow the steps below.
We'll illustrate the steps using ``<SECTION>``, which could be any of the
sections above.

1. Amend the dictionary ``DEFAULT_<SECTION>`` in :class:`~fitbenchmarking.utils.options.Options` to include any new default options.

2. If the option amended is to be checked for validity, add accepted option values to the ``VALID_<SECTION>`` dictionary in :class:`~fitbenchmarking.utils.options.Options`.
   
3. Using the :meth:`~fitbenchmarking.utils.options.Options.read_value` function,
   add your new option to the class, following the examples already in
   :class:`~fitbenchmarking.utils.options.Options`.  The syntax of this function is:

   .. automethod:: fitbenchmarking.utils.options.Options.read_value
		    :noindex:

4. Add tests in the following way:

    - Each of the sections has it's own test file, for example,
      ``test_option_fitting`` has tests for the ``FITTING`` section.

    - Add default tests to the class called ``<SECTION>OptionTests``.

    - Add user defined tests to the class called ``User<SECTION>OptionTests``. These
      should check that the user added option is valid and raise an ``OptionsError``
      if not.
      
5. Add relevant documentation for the new option in :ref:`options`.

Adding new Sections is also possible.  To do this you'll need to extend
``VALID_SECTIONS`` with the new section, and follow the same structure as the
other SECTIONS.
.. _cost_function:

#########################
Adding new cost functions
#########################

*This section describes how to add cost functions to benchmarking in FitBenchmarking*

In order to add a new cost function, ``<cost_func>``,
you will need to:

1. Create ``fitbenchmarking/cost_func/<cost_func>_cost_func.py``,
   which contains a new subclass of
   :class:`~fitbenchmarking.cost_func.base_cost_func.CostFunc`.
   Then implement the methods:

    -  .. automethod:: fitbenchmarking.cost_func.base_cost_func.CostFunc.eval_cost()
              :noindex:

    -  .. automethod:: fitbenchmarking.cost_func.base_cost_func.CostFunc.jac_res()
              :noindex:

    -  .. automethod:: fitbenchmarking.cost_func.base_cost_func.CostFunc.jac_cost()
              :noindex:

    -  .. automethod:: fitbenchmarking.cost_func.base_cost_func.CostFunc.hes_res()
              :noindex:

    -  .. automethod:: fitbenchmarking.cost_func.base_cost_func.CostFunc.hes_cost()
              :noindex:

2. Document the available cost functions by:

  * adding ``<cost_func>`` to the ``cost_func_type`` option in :ref:`fitting_option`.
  * updating any example files in the ``examples`` directory
  * adding the new cost function to the :ref:`cost_func` user docs.

3. Create tests for the cost function in
   ``fitbenchmarking/cost_func/tests/test_cost_func.py``.

The :class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem` and :class:`~fitbenchmarking.cost_func.base_cost_func.CostFunc` classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When adding new cost functions, you will find it helpful to make use of the
following members of the :class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem`
class.

.. currentmodule:: fitbenchmarking.parsing.fitting_problem
.. autoclass:: fitbenchmarking.parsing.fitting_problem.FittingProblem
          :members: eval_model, data_x, data_y, data_e
          :noindex:

You will also find it useful to implement the subclass members of
:class:`~fitbenchmarking.cost_func.base_cost_func.CostFunc`,
:class:`~fitbenchmarking.jacobian.base_jacobian.Jacobian` and
:class:`~fitbenchmarking.hessian.base_hessian.Hessian`.

.. currentmodule:: fitbenchmarking.cost_func.base_cost_func
.. autoclass:: fitbenchmarking.cost_func.base_cost_func.CostFunc
          :members: eval_cost, jac_res, jac_cost, hes_res, hes_cost
          :noindex:

.. currentmodule:: fitbenchmarking.jacobian.base_jacobian
.. autoclass:: fitbenchmarking.jacobian.base_jacobian.Jacobian
          :members: eval
          :noindex:

.. currentmodule:: fitbenchmarking.hessian.base_hessian
.. autoclass:: fitbenchmarking.hessian.base_hessian.Hessian
          :members: eval
          :noindex:
#######################################
Extending FitBenchmarking Documentation
#######################################

FitBenchmarking is designed to be easily extendable to add new features
to the software.  Below we outline instructions for doing this.

.. toctree::
    :titlesonly:
    :caption: Contents:

    parsers
    adding_data
    cost_function
    hessian_extend
    jacobian_extend
    controllers
    options_extend
    outputs/index


.. _controllers:

#######################
Adding Fitting Software
#######################

Controllers are used to interface FitBenchmarking with the various fitting
packages. Controllers are responsible for converting the problem into a format
that the fitting software can use, and converting the result back to a
standardised format (numpy arrays). As well as this, the controller must be
written so that the fitting is separated from the preparation wherever possible
in order to give accurate timings for the fitting. Supported controllers are
found in ``fitbenchmarking/controllers/``.

In order to add a new controller, you will need to:

1. Give the software a name ``<software_name>``.  This will be used by users when
   selecting this software.

2. Create ``fitbenchmarking/controllers/<software_name>_controller.py``
   which contains a new subclass of
   :class:`~fitbenchmarking.controllers.base_controller.Controller`.

   .. note::
      Please note that if the fitting package being added uses Matlab, then the
      new controller should also inherit from the mixin class
      :class:`~fitbenchmarking.controllers.matlab_mixin.MatlabMixin`.

      .. autoclass:: fitbenchmarking.controllers.matlab_mixin.MatlabMixin
          :members: py_to_mat
          :noindex:

   The new controller should implement four functions, as well as initializing the dictionary ``algorithm_check``:

  - .. autoattribute:: fitbenchmarking.controllers.base_controller.Controller.algorithm_check
               :noindex:  

  -  .. automethod:: fitbenchmarking.controllers.base_controller.Controller.__init__()
                     :noindex:
  -  .. automethod:: fitbenchmarking.controllers.base_controller.Controller.setup()
              :noindex:
  -  .. automethod:: fitbenchmarking.controllers.base_controller.Controller.fit()
              :noindex:
  -  .. automethod:: fitbenchmarking.controllers.base_controller.Controller.cleanup()
              :noindex:

   By default, a controller does not accept Jacobian or Hessian information. If the controller
   being added can use hessians and/or jacobians, then the following controller attributes should
   be set:

   - .. autoattribute:: fitbenchmarking.controllers.base_controller.Controller.jacobian_enabled_solvers()
               :noindex:  
   - .. autoattribute:: fitbenchmarking.controllers.base_controller.Controller.hessian_enabled_solvers()
               :noindex: 

3. Add the new software to the default options, following the instructions in
   :ref:`options_extend`.

Your new software is now fully hooked in with FitBenchmarking, and you can compare
it with the current software.  You are encouraged to contribute this to the
repository so that other can use this package.  To do this need to follow our
:ref:`guidelines` and our :ref:`workflow`, and you'll also need to

4. Document the available minimizers (see :ref:`fitting_option`, :ref:`minimizer_option`),
   including licencing information, if appropriate.
   Note: make sure that you use ``<software_name>`` in these places so that the
   software links in the HTML tables link correctly to the documentation.
   Add the software to ``examples/all_software.ini``.

   You should also ensure that the available minimizers are catagorised correctly in ``self.algorithm_check``
   using the :ref:`algorithm type <algorithm_type>` options. Please refer to the :ref:`algorithms`
   page for more information about each algorithm type.

5. Create tests for the software in
   ``fitbenchmarking/controllers/tests/test_controllers.py``. If the package
   is ``pip`` installable then add the tests to the ``DefaultControllerTests`` class
   and if not add to the ``ExternalControllerTests`` class.
   Unless the new controller is more complicated than the currently available
   controllers, this can be done by following the example of the others.

6. If the software is deterministic, add the software to the regression tests in
   ``fitbenchmarking/systests/test_regression.py``. 
   
7. If `pip` installable add to ``install_requires`` in ``setup.py`` and
   add to the installation step in ``.github/workflows/release.yml``.
   If not, document the installation procedure in :ref:`external-instructions`
   and update the ``FullInstall`` Docker Container -- the main developers will
   help you with this.

.. note::
   For ease of maintenance, please add new controllers to a list of
   software in alphabetical order.


The :class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem`,  :class:`~fitbenchmarking.cost_func.base_cost_func.CostFunc` and :class:`~fitbenchmarking.jacobian.base_jacobian.Jacobian` classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When adding new minimizers, you will find it helpful to make use of the
following members of the
:class:`~fitbenchmarking.parsing.fitting_problem.FittingProblem`, subclasses of
:class:`~fitbenchmarking.cost_func.base_cost_func.CostFunc` and subclasses of
:class:`~fitbenchmarking.jacobian.base_jacobian.Jacobian` classes:

.. currentmodule:: fitbenchmarking.parsing.fitting_problem
.. autoclass:: fitbenchmarking.parsing.fitting_problem.FittingProblem
          :members: eval_model, data_x, data_y, data_e, set_value_ranges
          :noindex:

.. currentmodule:: fitbenchmarking.cost_func.base_cost_func
.. autoclass:: fitbenchmarking.cost_func.base_cost_func.CostFunc
          :members: eval_cost
          :noindex:

.. currentmodule:: fitbenchmarking.jacobian.base_jacobian
.. autoclass:: fitbenchmarking.jacobian.base_jacobian.Jacobian
          :members: eval, eval_cost
          :noindex:
.. _tables:

#####################
Adding further Tables
#####################


The tables that are currently supported are listed in :ref:`output`.
In order to add a new table, you will need to:

1. Give the table a name ``<table_name>``. This will be used by users when
   selecting this output from FitBenchmarking.
2. Create ``fitbenchmarking/results_processing/<table_name>_table.py``
   which contains a new subclass of
   :class:`~fitbenchmarking.results_processing.base_table.Table`.
   The main functions to change are:

   - .. automethod:: fitbenchmarking.results_processing.base_table.Table.get_value
        :noindex:

   - .. automethod:: fitbenchmarking.results_processing.base_table.Table.display_str
        :noindex:

   Additional functions that may need to be overridden are:
   
   - .. automethod:: fitbenchmarking.results_processing.base_table.Table.get_error_str
        :noindex:

   - .. automethod:: fitbenchmarking.results_processing.base_table.Table.get_link_str
        :noindex:

   - .. automethod:: fitbenchmarking.results_processing.base_table.Table.vals_to_colour
        :noindex:

3. Extend the ``table_type`` option in ``PLOTTING`` following the instructions in
   :ref:`options_extend`.
	   
4. Document the new table class is by setting the docstring to be
   the description of the table, and add to :ref:`output`.
   
5. Create tests for the table in
   ``fitbenchmarking/results_processing/tests/test_tables.py``. This is done
   by generating, ahead of time using the results problems constructed in
   ``fitbenchmarking/results_processing/tests/test_tables.generate_mock_results``, both a HTML and text table output as the expected
   result and adding the new table name to the global variable
   ``SORTED_TABLE_NAMES``. This will automatically run the comparison tests for the tables.

.. _templates:

##################
HTML/CSS Templates
##################

In FitBenchmarking, templates are used to generate all html output files,
and can be found in the `fitbenchmarking/templates` directory.

==============
HTML Templates
==============
HTML templates allow for easily adaptable outputs.
In the simple case of rearranging the page or making static changes
(changes that don't vary with results), this can be done by editing the
template, and it will appear in every subsequent HTML page.

For more complicated changes such as adding information that is page dependent,
we use `jinja <https://jinja.palletsprojects.com/en/2.11.x/>`__.
Jinja allows code to be added to templates which can contain conditional
blocks, replicated blocks, or substitutions for values from python.

Changes to what information is displayed in the page will usually involve
editing the python source code to ensure you pass the appropriate values to
jinja. In practice the amount of effort required to do this can range from
one line of code to many depending on the object to be added.

=============
CSS Templates
=============
The CSS files contained in the `templates` directory are used to format the
HTML pages.

If you wish to change the style of the output results
(e.g. to match a website), this is where they can be changed.
.. _extending_outputs:

################################
Amending FitBenchmarking Outputs
################################

Here we describe how to add ways of displaying results obtained by FitBenchmarking.

.. toctree::
    :titlesonly:
    :maxdepth: 2
    :caption: Contents:

    tables
    templates
.. _how:

##############################
How does FitBenchmarking work?
##############################

FitBenchmarking takes data and models from real world applications
and data analysis packages.
It fits the data to the models by casting them as a
nonlinear least-squares problem.
We fit the data using a range of data fitting and
nonlinear optimization software, and present comparisons on the
accuracy and timings.

.. figure:: ../../images/FitBenchmarkingConcept.png
   :alt: FitBenchmarking Concept
   :width: 100.0%

*************************
The Benchmarking Paradigm
*************************

FitBenchmarking can compare against any of the supported minmizers listed in
:ref:`minimizer_option`.  We've also made it straightforward to add new software by
following the instructions in :ref:`controllers` -- the software just needs
to be callable from  Python.


Once you have chosen which minimizers you want to compare for a given problem,
running FitBenchmarking will give you a comparison to indicate the
minimizer that performs best.

There are a number of options that you can pick to customize what your tests
are comparing, or how they are run.  A full list of these options, and how to
select them, is given in the section :ref:`options`.

FitBenchmarking creates tables, as given in the section :ref:`output`,
which show a comparison between the different minimizers available.
An example of a table is:

.. figure:: ../../images/example_table.png
   :alt: Example Table

This is the result of FitBenchmarking for a selection of software/minimizers
and different problem definition types supported in FitBenchmarking.
Both the raw chi squared values, and the values normalised with respect
to the best minimizer per problem, are given.
The problem names link to html pages that display plots of the
data and the fit that was performed, together with initial and final
values of the parameters. Here is an example of the final plot fit:

.. figure:: ../../images/example_plot.png
   :alt: Example Plot

.. _performance_profile:

Performance Profile
-------------------

With each test FitBenchmarking also produces a Dolan-Moré performance profile:

.. figure:: ../../images/example_pp.png
	    :alt: Example Performance Profile

Fits are taken from all benchmarks, so if FitBenchmarking is run with
``n`` problems and ``m`` cost functions, the resulting profile plots will have
``n*m`` steps on the y-axis.

The solvers appearing in the top left corner may be considered the best
performing on this test set.
See `Dolan and Moré (2001) <https://link.springer.com/article/10.1007/s101070100263>`_
for more information. 





.. _why:

#################################
Why is FitBenchmarking important?
#################################

Fitting a mathematical model to data is a fundamental task across all
scientific disciplines.  (At least) three groups of people have an interest
in fitting software:

-  **Scientists**, who want to know what is the best algorithm for fitting
   their model to data they might encounter, on their specific hardware;

-  **Scientific software developers**, who want to know what is the
   state-of-the-art in fitting algorithms and implementations,
   what they should recommend as their default solver, and if they should
   implement a new method in their software; and

-  **Mathematicians** and **numerical software developers**, who want to understand the
   types of problems on which current algorithms do not perform well,
   and to have a route to expose newly developed methods to users.

Representatives of each of these communities have got together to build FitBenchmarking.
We hope this tool will help foster fruitful interactions and collaborations across the disciplines.


Example workflow
----------------

The black crosses on the plot below are data obtained from an experiment
at the VESUVIO beamline at ISIS Neutron and Muon source:

.. figure:: ../../images/start_for_EVS14188-90_processed_Gaussian_peaks_1_1.png  
   :alt: VESUVIO experiment data
   :width: 100.0%

   VESUVIO experiment data

The scientist needs to interpret this data, and will typically
use a data analysis package to help with this. Such packages are
written by specialist scientific software developers, who are experts in
analysing the kind of data produced by a given experiment;
examples include `Mantid <https://mantidproject.org/>`_,
`SasView <https://www.sasview.org>`_, and `Horace <https://horace.isis.rl.ac.uk>`_.

These packages include mathematical models, which depend on parameters,
that can describe the data.
We need to find values for the parameters in these models which
best fit the data -- for more background, see this `Wikipedia article <https://en.wikipedia.org/wiki/Goodness_of_fit>`_.
The usual way this is done is by finding parameters that minimize the
(weighted) squares of the error in the data, or :math:`\chi^2` value.
This is equivalent to formulating a nonlinear
least-squares problem; specifically, given :math:`n` data points
:math:`(x_i, y_i)` (the crosses in the figure above), together
with estimates of the errors on the values of :math:`y_i`, :math:`\sigma_i`, we solve

.. math:: {\boldsymbol{\beta}}^* = \arg \min_{{\boldsymbol{\beta}}} \underbrace{\sum_i \left( \frac{y_i - f({\boldsymbol{\beta}};x_i)}{\sigma_i} \right)^2}_{\chi^2({\boldsymbol{\beta}})},\label{eq:chi2}

where :math:`f({\boldsymbol{\beta}};x)` is the model we’re trying to
fit, and :math:`\boldsymbol{\beta}` are the parameters we're trying to
find.

Usually the scientist will supply a starting guess,
:math:`{\boldsymbol{\beta}}_0` (the pink curve in the graph above),
which describes where they think the solution might be.
She then has to *choose which algorithm to use to fit the curve*
from the selection available in the analysis software.
Different algorithms may be more or
less suited to a problem, depending on factors such as the architecture
of the machine, the availability of first and second derivatives, the
amount of data, the type of model used, etc.

Below we show the data overlayed by a blue curve, which is a model fitted using the
implementation of the Levenberg-Marquardt algorithm from the GNU Scientific Library (:code:`lmsder`).
The algorithm claims to have found a local minimum with a Chi-squared error of 
0.4771 in 1.9 seconds.

.. figure:: ../../images/lmsder_fit_for_EVS14188-90_processed_Gaussian_peaks_1_1.png
   :alt: VESUVIO experiment data: :code:`lmsder`
   :width: 100.0%

   GSL's :code:`lmsder` (Levenberg-Marquardt) algorithm on the data

We also solved the nonlinear least squares problem using GSL's implementation of
a Nelder-Mead simplex algorithm (:code:`nmsimplex2`), which again claimed to solve
the problem, this time in a faster 1.5 seconds.  However, this time the Chi-squared error was
0.8505, and we plot the curve obtained in green below.  The previous curve
is in dotted-blue, for comparison.
   
.. figure:: ../../images/nmsimplex2_fit_for_EVS14188-90_processed_Gaussian_peaks_1_1.png
   :alt: VESUVIO experiment data: :code:`nmsimplex2`
   :width: 100.0%

   GSL's :code:`nmsimplex2` (Nelder-Mead Simplex) algorithm on the data

By eye it is clear that the solution given by :code:`lmsder` is better.
As the volume of data increases, and we do more and more data analysis
algorithmically, it is increasingly important that we have the best algorithm
without needing to check it by eye.  

FitBenchmarking will help the **scientist** make an informed choice by
comparing runtime and accuracy of all available minimizers, on their
specific hardware, on problems from their science area, which will
ensure they are using the most appropriate minimizer. 

FitBenchmarking will help the **scientific software developer** ensure
that the most robust and quickest algorithms for the type of data
analysis they support are available in their software.

FitBenchmarking will help **mathematicians** see what the state of the
art is, and what kinds of data are problematic.  It will give
them access to real data, and will give a route for novel methods to
quickly make it into production.

A workflow as described above plays a crucial role in the processing and analysis of
data at large research facilities in tasks as diverse as instrument
calibration, refinement of structures, and data analysis methods specific
to different scientific techniques. FitBenchmarking will ensure that, across
all areas that utilise least-squares fitting, scientists can be confident they are
using the best tool for the job.
   
We discuss the specific
FitBenchmarking paradigm in the Section :ref:`how`
#####################################
FitBenchmarking Concept Documentation
#####################################

Here we outline why we built the fitbenchmarking software, and how the software
benchmarks minimizers with the goal of highlighting the best tool for different
types of data.

.. toctree::
    :titlesonly:
    :caption: Contents:

    why
    how
.. _tests:

#####################
FitBenchmarking Tests
#####################

The tests for FitBenchmarking require ``pytest>=3.6``. We have split the tests
into three categories:

* ``default``: denotes tests involving ``pip`` installable
  :ref:`software packages<getting-started>`,
* ``all``: in addition to ``default``, also runs tests on
  :ref:`external packages <external-instructions>`, with the
  exception of matlab.
* ``matlab``: Runs tests for matlab fitting software. Please
  note that these tests can currently only be run locally
  through pytest.

Unit tests
----------

Each module directory in FitBenchmarking (e.g. ``controllers``) contains a
test folder which has the ``unit`` tests for that module.
One can run the tests for a module by:

.. code-block:: bash

   pytest fitbenchmarking/<MODULE_DIR> --test-type <TEST_TYPE>

where <TEST_TYPE> is either ``default`` or ``all``.
If ``--test-type`` argument is not given the default is ``all``

System tests
------------

System tests can be found in the ``systests`` directory in FitBenchmarking.
As with the unit tests, these can be run via:

.. code-block:: bash

   pytest fitbenchmarking/systests --test-type <TEST_TYPE>


.. warning::
   The files in the expected results subdirectory of the ``systests``
   directory are generated to check consistency in our automated tests via
   `GitHub Actions <https://github.com/fitbenchmarking/fitbenchmarking/actions>`__.
   They might not pass on your local operating system due to, for example,
   different software package versions being installed.

GitHub Actions tests
---------------------

The scripts that are used for our automated tests via
`GitHub Actions <https://github.com/fitbenchmarking/fitbenchmarking/actions>`__
are located in the ``ci`` folder.
These give an example of how to run both the unit and system tests within
FitBenchmarking.
.. _notes:

############
Known Issues
############

This page is used to detail any known issues or unexpected behaviour
within the software.


************************************
Problem-Format/Software Combinations
************************************

When comparing minimizer options from one software package
(e.g., comparing all `scipy_ls` minimizers), we are not aware of any issues.
However, the general problem of comparing minimizers from multiple software
packages, and with different problem-formats, on truly equal terms is harder to
achieve.

The following list details all cases where we are aware of a possible bias:

- **Using native FitBenchmarking problems with the Mantid software and fitting using Mantid.**

  With Mantid data, the function evaluation is slightly faster for Mantid minimizers
  than for all other minimizers. You should account for this when interpreting the
  results obtained in this case.

- **Using non-scalar ties in native FitBenchmarking problems with the Mantid software.**

  Mantid allows parameters to be tied to expressions - e.g. X0=5.0 or X0=X1*2.
  While scalar ties are now supported for all minimizers the more complicated
  expressions are not supported. If you need this feature please get in touch
  with the development team with your use case.

- **Running Mantid problems with Matlab fitting software.**

  To run problems with Matlab fitting software through FitBenchmarking, within
  the Matlab Controller the dynamically created `cost_func.eval_model` function
  is serialized and then loaded in the Matlab Engine workspace. However for
  Mantid problems, this function is not picklable resulting in the problem
  being skipped over.

In all cases, the stopping criterion of each minimizer is set to the default
value.
An experienced user can change this.


***************************************
Specific Problem/Minimizer Combinations
***************************************

- **CrystalField Example with Mantid - DampedGaussNewton Minimizer.**

  With this combination, GSL is known to crash during Mantid's fitting.
  This causes python to exit without completing any remaining runs or
  generating output files.
  More information may be available via
  `the issue on Mantid's github page <https://github.com/mantidproject/mantid/issues/31176>`__.
.. _cost_func:

==============
Cost functions
==============

Fitbenchmarking supports multiple cost functions. These can be set via the ``cost_func_type`` option in :ref:`fitting_option`.

Fitbenchmarking is designed to work with problems that have the form

.. math::

   \min_p F(r(x,y,p)).
   
The function :math:`F(\cdot)` is known as the cost function,
while the function :math:`r(x,u,p)` is known as the *residual* of the cost function.
The residual will generally be zero if the fit was perfect.
Both of these quantities together define a cost function in FitBenchmarking.

The cost functions that are currently supported are:

- Non-linear least squares cost function

    .. currentmodule:: fitbenchmarking.cost_func.nlls_cost_func
    .. autoclass:: fitbenchmarking.cost_func.nlls_cost_func.NLLSCostFunc
               :noindex:
- Weighted non-linear least squares cost function

    .. currentmodule:: fitbenchmarking.cost_func.weighted_nlls_cost_func
    .. autoclass:: fitbenchmarking.cost_func.weighted_nlls_cost_func.WeightedNLLSCostFunc
               :noindex:
- Hellinger non-linear least squares cost function

    .. currentmodule:: fitbenchmarking.cost_func.hellinger_nlls_cost_func
    .. autoclass:: fitbenchmarking.cost_func.hellinger_nlls_cost_func.HellingerNLLSCostFunc
               :noindex:

- Poisson deviance cost function

    .. currentmodule:: fitbenchmarking.cost_func.poisson_cost_func
    .. autoclass:: fitbenchmarking.cost_func.poisson_cost_func.PoissonCostFunc
               :noindex:
.. _BenchmarkProblems:

====================
 Benchmark problems
====================

To help choose between the different minimizers, we have made some curated
problems available to use with FitBenchmarking.  It is also straightforward to
add custom data sets to the benchmark, if that is more appropriate; see
:ref:`problem_def` for specifics of how to add additional problems in a
supported file format.

.. topic:: Downloads

    **You can download a folder containing all examples here:**
    :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/examples.zip>`
    or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/examples.tar.gz>`

    Individual problem sets are also available to download below.

We supply some standard nonlinear least-squares test problems in the
form of the `NIST nonlinear regression set <https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml>`_
and the relevant problems from the `CUTEst problem set <https://github.com/ralna/CUTEst/wiki>`_,
together with some real-world 
data sets that have been extracted from `Mantid <https://www.mantidproject.org>`__ and
`SASView <https://www.sasview.org>`__ usage examples and system tests.
We've made it possible to extend this list by following the steps in 
:ref:`parsers`.

Each of the test problems contain:

* a data set consisting of points :math:`(x_i, y_i)` (with optional errors on :math:`y_i`, :math:`\sigma_i`);
* a definition of the fitting function, :math:`f({\boldsymbol{\beta}};x)`; and
* (at least) one set of initial values for the function parameters :math:`{\boldsymbol{\beta}}_0`.
  
If a problem doesn't have observational
errors (e.g., the NIST problem set), then FitBenchmarking can
approximate errors by taking :math:`\sigma_i = \sqrt{y_i}`.
Alternatively, there is an option to disregard errors and solve the
unweighted nonlinear least-squares problem, setting
:math:`\sigma_i = 1.0` irrespective of what has been passed in with the
problem data.

As we work with scientists in other areas, we will extend the problem
suite to encompass new categories. The FitBenchmarking framework has
been designed to make it easy to integrate new problem sets, and any
additional data added to the framework can be tested with any and all of
the available fitting methods.

Currently FitBenchmarking ships with data from the following sources:


CrystalField Data (Mantid)
==========================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/CrystalField.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/CrystalField.tar.gz>`

This folder (also found in `examples/benchmark_problems/CrystalField`) contains
a test set for inelastic neutron scattering measurements of transitions between
crystal field energy levels.

This problem has 8 parameters, and fits around 200 data points. 

.. warning::
    |MantidWarning|

CUTEst (NIST files)
===================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/CUTEst.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/CUTEst.tar.gz>`

This folder (also found in `examples/benchmark_problems/CUTEst`) contains
several problems from the `CUTEst <https://github.com/ralna/CUTEst>`_
continuous optimization testing environment which have been converted to the NIST
format.

These problems all have 8 unknown parameters, and fit around 15 data points
with the exception of ``VESUVIOLS`` which fits around 1000.

Data Assimilation
=================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/Data_Assimilation.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/Data_Assimilation.tar.gz>`

This folder (also found in `examples/benchmark_problems/Data_Assimilation`) contains
two examples using the data assimilation problem definition in fitbenchmarking.
These examples follow the method set out in 
`this paper <https://www.researchgate.net/publication/324956488_Data_assimilation_approach_to_analysing_systems_of_ordinary_differential_equations>`_.

These data files are synthetic and have been generated as an initial test of
the minimizers. We plan to extend this with time series data which is more
representative of the expectations for data assimilation in future updates.

These problems have either 2 or 3 unknown parameters, and fit either 100 or
1000 data points for ``Simplified ANAC`` and ``Lorentz`` problems respectively.


Powder Diffraction Data (SIF files)
===================================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/DIAMOND_SIF.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/DIAMOND_SIF.tar.gz>`

These problems (also found in the folder `examples/benchmark_problems/DIAMOND_SIF`)
contain data from powder diffraction experiments.  The data supplied comes
from the `I14 Hard X-Ray Nanoprobe <https://www.diamond.ac.uk/Instruments/Imaging-and-Microscopy/I14.html>`_ beamline at
the Diamond Light source, and has been supplied in the SIF
format used by `CUTEst <https://github.com/ralna/CUTEst>`_.

These problems have either 66 or 99 unknown parameters, and fit around 5,000 data points.


.. warning::
    |CUTEstWarning|

   
MultiFit Data (Mantid)
======================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/MultiFit.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/MultiFit.tar.gz>`

These problems (also found in the folder `examples/benchmark_problems/MultiFit`)
contain data
for testing the MultiFit functionality of Mantid.  This contains
a simple data set, on which two fits are done, and a calibration
dataset from the `MuSR <https://www.isis.stfc.ac.uk/Pages/musr.aspx>`_
spectrometer at ISIS, on which there are four fits available.
See :ref:`The MultiFit documentation<multifit>` for more details.

Basic Multifit has 3 unknown parameters, and fits 40 data points.
MUSR62260 has 18 unknown parameters, and fits around 8000 data points.

.. warning::
    |MantidWarning|
   
    This will also only work using the :ref:`Mantid Minimizers<MantidMinimizers>`.

Muon Data (Mantid)
==================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/Muon.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/Muon.tar.gz>`


These problems (also found in the folder `examples/benchmark_problems/Muon`)
contain data from Muon spectrometers.  The data supplied comes
from the `HiFi <https://www.isis.stfc.ac.uk/Pages/hifi.aspx>`_ and 
`EMU <https://www.isis.stfc.ac.uk/Pages/EMU.aspx>`_ instruments at
STFC's ISIS Neutron and Muon source, and has been supplied in the
format that `Mantid <https://mantidproject.org/>`__ uses to process
the data.

These problems have between 5 and 13 unknown parameters, and fit around 1,000 data points.

.. warning::
    |MantidWarning|


Neutron Data (Mantid)
=====================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/Neutron.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/Neutron.tar.gz>`

These problems (also found in the folder `examples/benchmark_problems/Neutron`)
contain
data from Neutron scattering experiments.  The data supplied comes
from the `Engin-X <https://www.isis.stfc.ac.uk/Pages/Engin-X.aspx>`_,
`GEM <https://www.isis.stfc.ac.uk/Pages/gem.aspx>`_,
`eVS <https://www.isis.stfc.ac.uk/Pages/Vesuvio.aspx>`_, and
`WISH <https://www.isis.stfc.ac.uk/Pages/wish.aspx>`_ instruments at
STFC's ISIS Neutron and Muon source, and has been supplied in the
format that `Mantid <https://mantidproject.org/>`__ uses to process
the data.

The size of these problems differ massively.
The Engin-X calibration problems find 7 unknown parameters, and fit to
56-67 data points.
The Engin-X vanadium problems find 4 unknown parameters, and fit to around 14,168
data points.
The eVS problems find 8 unknown parameters, and fit to 1,025 data points.
The GEM problem finds 105 unknown parameters, and fits to 1,314 data points.
The WISH problems find 5 unknown parameters, and fit to 512 data points.

.. warning::
    |MantidWarning|


NIST
====

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/NIST.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/NIST.tar.gz>`

These problems (also found in the folder `examples/benchmark_problems/NIST`) contain
data from the `NIST Nonlinear Regression <https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml>`_ test set.

These problems are split into low, average and high difficulty.
They have between 2 and 9 unknown parameters, and
fit between 6 and 250 data points.


Poisson Data
============

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/Poisson.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/Poisson.tar.gz>`

These problems (also found in the folder `examples/benchmark_problems/Poisson`) contain
both simulated and real data measuring particle counts. The real data is ISIS
muon data, and the simulated datasets have been made to represent counts using
models provided by both Mantid and Bumps.

These problems have between 4 and 6 unknown parameters, and around 350, 800,
and 2000 data points for simulated bumps, HIFI_160973, and simulated mantid
respectively.

.. warning::
    |MantidWarning|

Small Angle Scattering (SASView)
================================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/SAS_modelling.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/SAS_modelling.tar.gz>`


These problems (also found in the folder `examples/benchmark_problems/SAS_modelling/1D`) are
two data sets from small angle scattering experiments.
These are from fitting data to a
`cylinder <https://www.sasview.org/docs/user/models/cylinder.html>`_,
and have been supplied in the format that `SASView <https://www.sasview.org>`__
uses to process the data.

These have 6 unknown parameters, and fit to either 20 or 54 data points.

.. warning::
    The external package ``sasmodels`` must be installed to run this data
    set.  See :ref:`external-instructions` for details.


CUTEst (SIF files)
==================

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/SIF.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/SIF.tar.gz>`

This directory (also found in the folder `examples/benchmark_problems/SIF`) contain
`SIF files <https://github.com/ralna/SIFDecode>`_
encoding least squares problems 
from the `CUTEst <https://github.com/ralna/CUTEst>`_
continuous optimization testing environment.

These are from a wide range of applications.  They have between
2 and 9 unknown parameters, and for the most part fit between
6 and 250 data points, although the `VESUVIO` examples (from
the `VESUVIO <https://www.isis.stfc.ac.uk/Pages/Vesuvio.aspx>`_
instrument at ISIS) have 1,025 data points (with 8 unknown parameters).

.. warning::
    |CUTEstWarning|


SIF_GO
======

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/SIF_GO.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/SIF_GO.tar.gz>`

This directory (also found in the folder `examples/benchmark_problems/SIF_GO`) contains
`SIF files <https://github.com/ralna/SIFDecode>`_
encoding least squares problems 
from the `CUTEst <https://github.com/ralna/CUTEst>`_
continuous optimization testing environment.

All of these problems have been modified, with finite bounds added for all parameters,
making the problems appropriate for testing global optimization solvers. The bounds that
have been added to each problem are the same as those used in SciPy's
`global optimization benchmark functions <https://github.com/scipy/scipy/tree/master/benchmarks/benchmarks/go_benchmark_functions>`_.

These problems have between 3 and 7 unknown parameters, and fit between 9 and 37 data points.

.. warning::
    |CUTEstWarning|


Simple tests
============

**Download** :download:`.zip <https://numerical.rl.ac.uk/fitbenchmarking/simple_tests.zip>`
or :download:`.tar.gz <https://numerical.rl.ac.uk/fitbenchmarking/simple_tests.tar.gz>`

This folder (also found in `examples/benchmark_problems/simple_tests`) contains
a number of simple tests with known, and easy to obtain,
answers.  We recommend that this is used to test any new minimizers
that are added, and also that any new parsers reimplement these
data sets and models (if possible).

These problems have 3 or 4 unknown parameters, and around 100 data points.

.. |CUTEstWarning| replace::
    The external packages CUTEst and pycutest must be installed to run
    this data set.   See :ref:`external-instructions` for details.

.. |MantidWarning| replace::
    The external package Mantid must be installed to run
    this data set.  See :ref:`external-instructions` for details.
.. _running:

#######################
Running FitBenchmarking
#######################

Once installed, issuing the command

.. code-block:: bash

   fitbenchmarking

will run the NIST average difficulty problem set on SciPy minmizers.

Running alternative problems
----------------------------

Other problems written in a :ref:`supported file format <problem_def>`
can be analyzed with FitBenchmarking by
passing the path using the ``-p`` or ``--problem-sets`` argument.
Example problems can be downloaded from
:ref:`BenchmarkProblems`, and they can also be found in the
``fitbenchmarking/examples`` directory of the code.

For example, to run the NIST low difficulty set from the base directory
of the source, type into the terminal:

.. code-block:: bash
		
   fitbenchmarking -p examples/benchmark_problems/NIST/low_difficulty

Changing the options
--------------------
   
An options file can also be passed with the ``-o`` or ``--options-file`` argument. 
For example, the template file can be run by issuing the command

.. code-block:: bash

   fitbenchmarking -o examples/options_template.ini \
   -p examples/benchmark_problems/NIST/low_difficulty

Details about how the options file must be formatted are given in :ref:`options`.

.. _change_results_directory:

Changing the results directory
------------------------------

The default directory where the results are saved can be changed using the ``-r``
or ``--results-dir`` argument. The :ref:`results directory option <results_directory_option>`
can also be changed in the options file.

.. code-block:: bash

   fitbenchmarking -r new_results/

The default results directory is ``fitbenchmarking_results``.
##################################
FitBenchmarking User Documentation
##################################

In these pages we describe how to install, run, and set options to use
the FitBenchmarking software.

.. toctree::
    :titlesonly:
    :maxdepth: 2
    :caption: Contents:

    install_instructions/index
    running
    algorithms/index
    cost_func
    output/index
    BenchmarkProblems
    options/index
    problem_definition_files/index
    tests
    notes

.. _minimizer_types:

***************************************
Algorithm Types of Available Minimizers
***************************************

Least Squares (``ls``):

.. algorithmcheckdocs::
    :key: ls

Deriv-Free (``deriv_free``):

.. algorithmcheckdocs::
    :key: deriv_free

General (``general``):

.. algorithmcheckdocs::
    :key: general

Simplex (``simplex``):

.. algorithmcheckdocs::
    :key: simplex

Trust Region (``trust_region``):

.. algorithmcheckdocs::
    :key: trust_region

Levenberg-Marquardt (``levenberg-marquardt``):

.. algorithmcheckdocs::
    :key: levenberg-marquardt

Gauss Newton (``gauss_newton``):

.. algorithmcheckdocs::
    :key: gauss_newton

BFGS (``bfgs``):

.. algorithmcheckdocs::
    :key: bfgs

Conjugate Gradient (``conjugate_gradient``):

.. algorithmcheckdocs::
    :key: conjugate_gradient

Steepest Descent (``steepest_descent``):

.. algorithmcheckdocs::
    :key: steepest_descent

Global Optimization (``global_optimization``): 

.. algorithmcheckdocs::
    :key: global_optimization.. _line_search:

*******************
Line Search Methods
*******************
In line search methods, each iteration is given by :math:`x_{k+1} = x_k + \alpha_k p_k`, where :math:`p_k` is the search direction and :math:`\alpha_k` is the step length.

The search direction is often of the form :math:`p_k = -B_k^{-1} \nabla f_k` where :math:`B_k` is a symmetric and non-singular matrix. The form of :math:`p_k` is dependent on algorithm choice.

The ideal step length would be :math:`min_{\alpha>0} f(x_k + \alpha p_k)` but this is generally too expensive to calculate. Instead an inexact line search condition such as the Wolfe Conditions can be used:

.. math::
    f(x_k + \alpha p_k) \leq f(x_k) + c_1 \alpha \nabla f_k^T p_k \\
    f(x_k + \alpha_k p_k)^T p_k \geq c_2 \nabla f_k^T p_k

With :math:`0<c_1<c_2<1`. Here, the first condition ensures that :math:`\alpha_k` gives a sufficient decrease in :math:`f`, whilst the second condition rules out unacceptably short steps. [Nocedal]_

.. _steepest_descent:

Steepest Descent (``steepest_descent``)
***************************************
Simple method where search direction :math:`p_k` is set to be :math:`-\nabla f_k`, i.e. the direction along which :math:`f` decreases most rapidly.

**Advantages**:
    - Low storage requirements
    - Easy to compute

**Disadvantages**:
    - Slow convergence for nonlinear problems
[Nocedal]_

.. _conjugate_gradient:

Conjugate Gradient (``conjugate_gradient``)
*******************************************
Conjugate Gradient methods have a faster convergence rate than Steepest Descent but avoid the high computational cost of methods where the inverse Hessian is calculated.

Given an iterate :math:`x_0`, evaluate :math:`f_0 = f(x_0), \nabla f_0 = \nabla f(x_0)`.

Set :math:`p_0 \leftarrow - \nabla f_0, k \leftarrow 0`

Then while :math:`\nabla f_k \neq 0`:

Carry out a line search to compute the next iterate, then evaluate :math:`\nabla f_{k+1}` and use this to determine the subsequent conjugate direction :math:`p_{k+1} = - \nabla f(x_{k+1}) + \beta_k p_k`

Different variations of the Conjugate Gradient algorithm use different formulas for :math:`\beta_k`, for example:

Fletcher-Reeves: :math:`\beta_{k+1} = \frac{f_{k+1}^T \nabla f_{k+1}}{\nabla f_k^T \nabla f_k}`
Polak-Ribiere:  :math:`\beta_{k+1} = \frac{ \nabla f_{k+1}^T ( \nabla f_{k+1} - \nabla f_k)}{\|\nabla f_k\|^2}`

[Nocedal]_ [Poczos]_

**Advantages**:
    - Considered to be one of the best general purpose methods.
    - Faster convergence rate compared to Steepest Descent and only requires evaluation of objective function and it's gradient - no matrix operations.

**Disadvantages**:
    - For Fletcher-Reeves method it can be shown that if the method generates a bad direction and step, then the next direction and step are also likely to be bad. However, this is not the case with the Polak Ribiere method.
    - Generally, the Polak Ribiere method is more efficient that the Fletcher-Reeves method but it has the disadvantage is requiring one more vector of storage.
[Nocedal]_

.. _bfgs:

BFGS (``bfgs``)
***************
Most popular quasi-Newton method, which uses an approximate Hessian rather than the true Hessian which is used in a Newton line search method.

Starting with an initial Hessian approximation :math:`H_0` and starting point :math:`x_0`:

While :math:`\| \nabla f_k \| > \epsilon`:

Compute the search direction :math:`p_k = -H_k \nabla f_k`

Then find next iterate :math:`x_{k+1}` by performing a line search.

Next, define :math:`s_k = x_{k+1}-x_k` and :math:`y_k = \nabla f_{k+1} - \nabla f_k`, then compute

.. math::
    H_{k+1} = (I - \rho_k s_k y_k^T)H_k(I - \rho_k y_k s_K^T) + \rho_k s_k s_k^T

with :math:`\rho_k = \frac{1}{y_k^T s_k}`

**Advantages**:
    - Superlinear rate of convergence
    - Has self-correcting properties - if there is a bad estimate for :math:`H_k`, then it will tend to correct itself within a few iterations.
    - No need to compute the Jacobian or Hessian.

**Disadvantages**:
    - Newton's method has quadratic convergence but this is lost with BFGS.
[Nocedal]_

.. _gauss_newton:

Gauss Newton (``gauss_newton``)
*******************************
Modified Newton's method with line search. Instead of solving standard Newton equations

.. math::
    \nabla^2 f(x_k)p = -\nabla f(x_k),
solve the system

.. math::
    J_k^T J_k p_k^{GN} = - J_k^T r_k

(where :math:`J_k` is the Jacobian) to obtain the search direction :math:`p_k^{GN}`. The next iterate is then set as :math:`x_{k+1} = x_k + p_k^{GN}`.

Here, the approximation of the Hessian :math:`\nabla^2 f_k \approx J_k^T J_k` has been made, which helps to save on computation time as second derivatives are not calculated.

**Advantages**:
    - Calculation of second derivatives is not required.
    - If residuals or their second order partial derivatives are small, then :math:`J_k^T J_k` is a close approximation to :math:`\nabla^2 f_k` and convergence of Gauss-Newton is fast.
    - The search direction :math:`p_J^{GN}` is always a descent direction as long as :math:`J_k` has full rank and the gradient :math:`\nabla f_k` is nonzero.

**Disadvantages**:
    - Without a good initial guess, or if the matrix :math:`J_k^T J_k` is ill-conditioned, the Gauss Newton Algorithm is very slow to converge to a solution.
    - If relative residuals are large, then large amounts of information will be lost.
    - :math:`J_k` must be full rank.
[Nocedal]_ [Floater]_

.. [Nocedal] Jorge Nocedal, Stephen J. Wright (2006), Numerical Optimization

.. [Poczos] Barnabas Poczos, Ryan Tibshirani (2012), Lecture 10: Optimization, School of Computer Science, Carnegie Mellon University

.. [Floater] Michael S. Floater (2018), Lecture 13: Non-linear least squares and the Gauss-Newton method, University of Oslo
.. _trust_region:

*************
Trust Region
*************

Trust region approach involves constructing a model function :math:`m_k` that approximates the function :math:`f` in the region, :math:`\Delta`, near the current point :math:`x_k`. 
The model :math:`m_k` is often a quadratic obtained by a Taylor Series expansion of the function around :math:`x_k`.

.. math::
    m_k(p) = f_k + \nabla f(x_k)^T p + \frac{1}{2} p^T B_k p

where :math:`B_k` is an approximation of the Hessian.

The subproblem to be solved at each iteration in order to find the step length is :math:`\min_p m_k(p)`, subject to :math:`\|p\| \leq \Delta_k`. [Nocedal]_

To select all minimizers in fitbenchmarking that use a trust region approach, use the algorithm type ``trust_region``.

.. _levenberg_marquardt:

Levenberg-Marquardt (``levenberg-Marquardt``)
*********************************************
Most widely used optimization algorithm, which uses the same Hessian approximation as Gauss-Newton but uses a trust region strategy instead of line search. As the Hessian approximation is the same as Gauss-Newton, convergence rate is similar.

For Levenberg-Marquardt, the model function :math:`m_k`, is chosen to be

.. math::
    m_k(p) = \frac{1}{2} \|r_k\|^2 + p^T J_k^T r_k + \frac{1}{2} p^T J_k^T J_k p

So, for a spherical trust region, the subproblem to be solved at each iteration is :math:`\min_p \frac{1}{2} \|J_k p + r_k\|^2`, subject to :math:`\|p\| \leq \Delta_k`.

Levenberg-Marquardt uses a combination of gradient descent and Gauss-Newton method. When the solution :math:`p^{GN}` lies inside of the trust region :math:`\Delta`, then :math:`p^{GN}` also solves the sub-problem. Otherwise, the current iteration is far from the optimal value and so the search direction is determined using steepest descent, which performs better than Gauss-Newton when far from the minimum.

**Advantages**:
    - Robust (more so than Gauss-Newton).
    - Avoids the weakness with Gauss-Newton that Jacobian must be full rank.
    - Fast to converge.
    - Good initial guess not required.

**Disadvantages**:
    - Similarly to Gauss-Newton, not good for large residual problems.
    - Can be slow to converge if a problem has many parameters
[Nocedal]_ [Ranganathan]_

.. [Nocedal] Jorge Nocedal, Stephen J. Wright (2006), Numerical Optimization

.. [Ranganathan] Ananth Ranganathan (2004), The Levenberg-Marquardt Algorithm, University of California, Santa Barbara
.. _algorithms:

#######################
Optimization Algorithms
#######################

The different minimizers used in Fitbenchmarking implement various numerical optimization algorithms. These are catagorised
by the algorithm type list, with options detailed below, along with advantages/disadvantages of each algorithm.

.. toctree::
    :titlesonly:
    :maxdepth: 2
    :caption: Algorithm Types:

    deriv_free
    line_search
    trust_region

    minimizer_types

.. _deriv_free:

****************
Derivative Free
****************

Derivative Free methods do not compute the gradient of a function and so are often used to minimize problems with
nondifferentiable functions. Some derivative free methods will attempt to approximate the gradient using a finite difference
approach, another class of methods constructs a linear or quadratic model of the objective functions and uses a trust
region approach to find the next iterate. Another widely used derivative free method is the Nelder-Mead simplex method. [Nocedal]_

To select all minimizers in fitbenchmarking that use derivative free methods, use the algorithm type ``deriv_free``.

.. _simplex:

Simplex (``simplex``)
*********************
Nelder-Mead is a simplex based algorithm, with a simplex :math:`S` in :math:`{\rm I\!R}` being defined as the convex hull of :math:`n+1` vertices :math:`\{x_1, ..., x_{n+1}\}`.

In an iteration of the algorithm, the idea is to remove the vertex with the worst function value. It is then replaced with another point with a better value. An iteration consists of the following steps:

1. **Ordering** the vertices of :math:`S` so that :math:`f(x_1) \leq f(x_2) \leq ... \leq f(x_{n+1})`

2. Calculating the **centroid**, :math:`\bar{x}` of the best :math:`n` points :math:`\bar{x} = \frac{1}{n} \sum_{i=1}^n x_i`

3. Carry out a **transformation** to compute the new simplex. Try to replace only the worst vertex :math:`x_{n+1}` with a better point, which then becomes the new vertex. If this fails, the simplex is shrunk towards the best vertex :math:`x_1` and :math:`n` new vertices are computed.
   The algorithm terminates when the simplex :math:`S` is sufficiently small. [Singer]_

**Advantages**: 
    - Method only requires 1 or 2 functions evaluations per iteration.
    - Gives significant improvements in first few iterations - quick to produce satisfactory results.

**Disadvantages**:
    - Stagnation can occur at non-optimal points, with large numbers of iterations leading to negligible improvement even when nowhere near a minimum.
    - If numerical computation of function derivative can be trusted, then other algorithms are more robust.
[Singer]_

.. [Nocedal] Jorge Nocedal, Stephen J. Wright (2006), Numerical Optimization

.. [Singer] Saša Singer, John Nelder (2009) Nelder-Mead algorithm. Scholarpedia, 4(7):2928.
.. _compare:

################
Comparison Table
################

.. currentmodule:: fitbenchmarking.results_processing.compare_table
.. autoclass:: fitbenchmarking.results_processing.compare_table.CompareTable
	       :noindex:
.. _problem_summary_page:

====================
Problem Summary Page
====================

The problem summary page can be used to give an overview of the problem and
solutions obtained.

Problem Outline
***************

First is the initial problem. Here you will see information about the function
being fit and the set of initial parameters used for the fitting.
If plots are enabled (see :ref:`MakePlots`), you will also see a scatter plot
of the data to fit with a line of the initial fit given to the minimizer.


Comparison
**********

The main plot on the page shows a comparison of all fits at once.
This can be used to compare how cost functions perform for a problem accross
all minimizers.

This uses colours to identify the cost function for each fit and shows all fits
on a single graph. The best minimizer for each cost function is more pronounced
on the plot.

This should not be used to identify the best individual fit, but can be a good
indication of whether cost functions are biased to certain datapoints in the
input.

Best Plots
**********

The page ends with an expandable section for each cost function tested, which
gives the parameter values and plot of the best fit obtained for that cost
function.
.. _acc:

##############
Accuracy Table
##############

.. currentmodule:: fitbenchmarking.results_processing.acc_table
.. autoclass:: fitbenchmarking.results_processing.acc_table.AccTable
	       :noindex:
.. fitting_report:

==============
Fitting Report
==============

The fitting report pages can be used to see more information about the problem
and a given fit.

Each page represents a single fitting combination and is split into 2 sections.

Problem Outline
***************

First is the initial problem. Here you will see information about the function
being fit and the set of initial parameters used for the fitting.
If plots are enabled (see :ref:`MakePlots`), you will also see a scatter plot
of the data to fit with a line of the initial fit given to the minimizer.

Fitting Results
***************

The second section focusses on the results of the fitting. Here you will find
the minimizer name and the final parameters for the fit found by the minimizer.
A plot of the fit is also shown with an overlaid best fit from whichever
minimizer was found to produce the smallest error.
.. _runtime:

#############
Runtime Table
#############

.. currentmodule:: fitbenchmarking.results_processing.runtime_table
.. autoclass:: fitbenchmarking.results_processing.runtime_table.RuntimeTable
	       :noindex:
.. _local_min:

#####################
Local Minimizer Table
#####################

.. currentmodule:: fitbenchmarking.results_processing.local_min_table
.. autoclass:: fitbenchmarking.results_processing.local_min_table.LocalMinTable
	       :noindex:
.. _output:

######################
FitBenchmarking Output
######################

FitBenchmarking produces tables and reports called support pages as outputs.
The links below give descriptions of these outputs.

Tables
******

.. toctree::
    :titlesonly:
    :maxdepth: 2

    compare
    acc
    runtime
    local_min

Display modes
-------------

The tables for ``accuracy``, ``runtime`` and ``compare`` have three display
modes:

.. prettyprintmodulevalue::
   :module: fitbenchmarking.results_processing.base_table
   :var: FORMAT_DESCRIPTION

This can be set in the option file using the
:ref:`Comparison Mode <ComparisonOption>` option.

The :ref:`local_min` table is formatted differently, and doesn't use this
convention.

Performance profile
-------------------

Below the table there is a :ref:`performance_profile`.

Support Pages
*************

In each of the tables, a :ref:`fitting_report` for an individual result can be accessed
by clicking on the associated table cell. Clicking the problem name at the
start of a row will open a :ref:`problem_summary_page` for the problem as a whole.

.. toctree::
    :titlesonly:
    :hidden:
    :maxdepth: 2

    fitting_report
    problem_summary_pages
.. _hessian_option:

###############
Hessian Options
###############

The Hessian section allows you to control which methods for computing Hessians the software uses.

Analytic (:code:`analytic`)
---------------------------

Analytic Hessians can only be used for specific :ref:`problem_def`. Currently
the supported formats are cutest and NIST. The only option is:

* ``default`` - use the analytic derivative provided by a supported format.

Default is ``default``

.. code-block:: rst

    [HESSIAN]
    analytic: default

.. _scipy-hes:

SciPy (:code:`scipy`)
---------------------

Calculates the Hessian from the Jacobian using the finite differencing in
SciPy, this uses ``scipy.optimize._numdiff.approx_derivative``. The supported
options are:

* ``2-point`` - use the first order accuracy forward or backward difference.
* ``3-point`` - use central difference in interior points and the second order accuracy forward or backward difference near the boundary.
* ``cs`` - use a complex-step finite difference scheme. This assumes that the user function is real-valued and can be analytically continued to the complex plane. Otherwise, produces bogus results.

Default is ``2-point``

**Licence** SciPy is available under a `3-clause BSD Licence <https://github.com/scipy/scipy/blob/master/LICENSE.txt>`__.  Individual packages may have their own (compatible) licences, as listed `here <https://github.com/scipy/scipy/blob/master/LICENSES_bundled.txt>`__.

.. code-block:: rst

    [HESSIAN]
    scipy: 2-point

.. _defaulthessian:

Default Hessian (:code:`default`)
---------------------------------

Hessian information is not passed to minimizers. The only option is:

* ``default`` - don't pass Hessian information to minimizers.

Default is ``default``

.. code-block:: rst

    [HESSIAN]
    default: default

.. _numdifftools-hes:

Numdifftools (:code:`numdifftools`)
-----------------------------------

Calculates the Hessian from the Jacobian using the python package :code:`numdifftools`.
We allow the user to change the method used, but other options
(e.g, the step size generator and the order of the approximation) are set to the defaults.
The supported options are:

* ``central`` - central differencing.  Almost as accurate as complex, but with no restriction on the type of function.
* ``forward`` - forward differencing.
* ``backward`` - backward differencing.
* ``complex`` - based on the complex-step derivative method of `Lyness and Moler <http://epubs.siam.org/doi/abs/10.1137/0704019>`__.  Usually the most accurate, provided the function is analytic.
* ``multicomplex`` - extends complex method using multicomplex numbers. (see, e.g., `Lantoine, Russell, Dargent (2012) <https://dl.acm.org/doi/10.1145/2168773.2168774>`__).

Default is ``central``.

**Licence** :code:`numdifftools` is available under a `3-clause BSD Licence <https://github.com/pbrod/numdifftools/blob/master/LICENSE.txt>`__.

.. code-block:: rst

    [HESSIAN]
    numdifftools: central
.. _output_option:

##############
Output Options
##############

The output section contains options to control how results are outputted.

.. _results_directory_option:

Results directory (:code:`results_dir`)
---------------------------------------

This is used to select where the output should be saved. If the 
:ref:`results directory command line argument <change_results_directory>` 
is provided, this option is overridden.

Default is ``fitbenchmarking_results``

.. code-block:: rst

    [OUTPUT]
    results_dir: fitbenchmarking_results.. _plotting_option:

################
Plotting Options
################

The plotting section contains options to control how results are presented.

.. _MakePlots:

Make plots (:code:`make_plots`)
-------------------------------

This allows the user to decide whether or not to create plots during runtime.
Toggling this to False will be much faster on large data sets.

Default is ``True`` (``yes``/``no`` can also be used)

.. code-block:: rst

    [PLOTTING]
    make_plots: yes

Colourmap (:code:`colour_map`)
------------------------------
Specifies the name of the colourmap the user wishes to use, e.g. ``magma``, ``viridis``, ``OrRd``. Options are:

* Any colourmap from the library in ``matplotlib``, see the complete library `here <https://matplotlib.org/stable/gallery/color/colormap_reference.html>`_.
* Appending ``_r`` to the end of the name will reverse the colourmap.
* The following sequential colourmaps are recommended:

.. figure:: ../../../images/recommended_perceptual_cmaps.png

.. figure:: ../../../images/recommended_sequential_cmaps.png

Default colourmap is ``magma_r``

.. code-block:: rst

    [PLOTTING]
    colour_map: magma_r



Colourmap Range (:code:`cmap_range`)
------------------------------------
A two-element list used to specify the lower and upper limit of the chosen colourmap. Options are:

* ``[lower_limit, upper_limit]`` where limits consider the full colourscale limits to be 0 and 1, so any pair of values must fall within this range.

* Limits should be introduced to make the **white text** readable, see the following example.

.. figure:: ../../../images/example_cmaps.png

Default for ``magma`` is ``[0.2, 0.8]`` (suitability depends on colourmap)

.. code-block:: rst

    [PLOTTING]
    colour_map: magma_r
    cmap_range: [0.2, 0.8] 

Colour Upper Limit (:code:`colour_ulim`)
----------------------------------------

Controls how relatively poorly a minimizer has to perform in order to receive the `worst` colour. For example,
a value of 100 would mean that any performance greater than or equal to 100 times worse than the best
minimizer would receive the worst colour. This ensures that colour scale is not compromised by especially 
poor relative results. Options are:

* Any float between ``1`` and ``np.inf``
* Recommended value ``100``

Default is ``100``

.. code-block:: rst

    [PLOTTING]
    colour_map: magma_r
    cmap_range: [0.2, 0.8] 
    colour_ulim: 100


.. _ComparisonOption:

Comparison mode (:code:`comparison_mode`)
-----------------------------------------

This selects the mode for displaying values in the resulting table
options are ``abs``, ``rel``, ``both``:

* ``abs`` indicates that the absolute values should be displayed
* ``rel`` indicates that the values should all be relative to the best result
* ``both`` will show data in the form "abs (rel)"

Default is ``both``

.. code-block:: rst

    [PLOTTING]
    comparison_mode: both


Table type (:code:`table_type`)
-------------------------------

This selects the types of tables to be produced in FitBenchmarking.
Options are:

* ``acc`` indicates that the resulting table should contain the chi squared values for each of the minimizers.
* ``runtime`` indicates that the resulting table should contain the runtime values for each of the minimizers.
* ``compare`` indicates that the resulting table should contain both the chi squared value and runtime value for each of the minimizers. The tables produced have the chi squared values on the top line of the cell and the runtime on the bottom line of the cell.
* ``local_min`` indicates that the resulting table should return true if a local minimum was found, or false otherwise.
  The value of :math:`\frac{|| J^T r||}{||r||}` for those parameters is also returned.
  The output looks like ``{bool} (norm_value)``, and the colouring is red for false and cream for true.
  This option is only meaningful for least-squares cost functions.

Default is ``acc``, ``runtime``, ``compare``, and ``local_min``.

.. code-block:: rst

    [PLOTTING]
    table_type: acc
                runtime
                compare
                local_min
.. _jacobian_option:

################
Jacobian Options
################

The Jacobian section allows you to control which methods for computing Jacobians the software uses.

Analytic (:code:`analytic`)
---------------------------

Analytic Jacobians can only be used for specific :ref:`problem_def`. Currently
the supported formats are cutest and NIST. The only option is:

* ``default`` - use the analytic derivative provided by a supported format.

Default is ``default``

.. code-block:: rst

    [JACOBIAN]
    analytic: default

.. _scipy-jac:

SciPy (:code:`scipy`)
---------------------

Calculates the Jacobian using the numerical Jacobian in
SciPy, this uses ``scipy.optimize._numdiff.approx_derivative``. The supported
options are:

* ``2-point`` - use the first order accuracy forward or backward difference.
* ``3-point`` - use central difference in interior points and the second order accuracy forward or backward difference near the boundary.
* ``cs`` - use a complex-step finite difference scheme. This assumes that the user function is real-valued and can be analytically continued to the complex plane. Otherwise, produces bogus results.

Default is ``2-point``

**Licence** SciPy is available under a `3-clause BSD Licence <https://github.com/scipy/scipy/blob/master/LICENSE.txt>`__.  Individual packages may have their own (compatible) licences, as listed `here <https://github.com/scipy/scipy/blob/master/LICENSES_bundled.txt>`__.

.. code-block:: rst

    [JACOBIAN]
    scipy: 2-point

.. _defaultjacobian:

Solver Default Jacobian (:code:`default`)
--------------------------------------------

This uses the approximation of the Jacobian that is used by default in the minimizer,
and will vary between solvers.  If the minimizer requires the user to pass a Jacobian,
a warning will be printed to the screen and the :ref:`scipy-jac` 2-point
approximation will be used.  The only option is:

* ``default`` - use the default derivative approximation provided by the software.

Default is ``default``

.. code-block:: rst

    [JACOBIAN]
    default: default

.. _numdifftools-jac:

Numdifftools (:code:`numdifftools`)
-----------------------------------

Calculates the Jacobian using the python package :code:`numdifftools`.
We allow the user to change the method used, but other options
(e.g, the step size generator and the order of the approximation) are set the defaults.
The supported options are:

* ``central`` - central differencing.  Almost as accurate as complex, but with no restriction on the type of function.
* ``forward`` - forward differencing.
* ``backward`` - backward differencing.
* ``complex`` - based on the complex-step derivative method of `Lyness and Moler <http://epubs.siam.org/doi/abs/10.1137/0704019>`__.  Usually the most accurate, provided the function is analytic.
* ``multicomplex`` - extends complex method using multicomplex numbers. (see, e.g., `Lantoine, Russell, Dargent (2012) <https://dl.acm.org/doi/10.1145/2168773.2168774>`__).

Default is ``central``.

**Licence** :code:`numdifftools` is available under a `3-clause BSD Licence <https://github.com/pbrod/numdifftools/blob/master/LICENSE.txt>`__.

.. code-block:: rst

    [JACOBIAN]
    numdifftools: central
.. _minimizer_option:

===================
 Minimizer Options
===================

This section is used to declare the minimizers to use for each fitting
software. If a fitting software has been selected in :ref:`fitting_option`
then a default set of minimizers for that solver will be run unless alternative
minimizer options have been set. All minimizers for a software are included on
the default list of minimizers unless otherwise stated.

.. warning::

   Options set in this section will only have an effect if the related
   software is also set in :ref:`fitting_option` (either explicitly, or
   as a default option).

.. _bumps:

Bumps (:code:`bumps`)
=====================

`Bumps <https://bumps.readthedocs.io>`__ is a set of data fitting (and Bayesian uncertainty analysis) routines.
It came out of the University of Maryland and NIST as part of the DANSE
(*Distributed Data Analysis of Neutron Scattering Experiments*) project.

FitBenchmarking currently supports the Bumps minimizers:

* `Nelder-Mead Simplex <https://bumps.readthedocs.io/en/latest/guide/optimizer.html#nelder-mead-simplex>`__ (:code:`amoeba`)

* `Levenberg-Marquardt <https://bumps.readthedocs.io/en/latest/guide/optimizer.html#fit-lm>`__  (:code:`lm`)

* `Quasi-Newton BFGS <https://bumps.readthedocs.io/en/latest/guide/optimizer.html#quasi-newton-bfgs>`__ (:code:`newton`)

* `Differential Evolution <https://bumps.readthedocs.io/en/latest/guide/optimizer.html#differential-evolution>`__ (:code:`de`)

* `MINPACK <https://github.com/bumps/bumps/blob/96b5100fc3d5b9485bd4a444c83a33617b74aa9d/bumps/mpfit.py>`__ (:code:`mp`)  This is a translation of `MINPACK` to Python.

**Licence** The main licence file for Bumps is `here <https://github.com/bumps/bumps/blob/master/LICENSE.txt>`__.  Individual files have their own copyright and licence
-- if you plan to incorporate this in your own software you should first check that the
licences used are compatible.
  
**Links** `GitHub - bumps <https://github.com/bumps/bumps>`__

The Bumps minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    bumps: amoeba
           lm-bumps
           newton
           de
           mp

.. warning::
   The additional dependency Bumps must be installed for this to be available;
   See :ref:`extra_dependencies`.	 

.. note::
   `de` is not included in the default list of minimizers for bumps. To run this solver, you must
   explicitly set the minimizer as seen above.
	   
.. _dfo:

DFO (``dfo``)
=============

There are two Derivative-Free Optimization packages, `DFO-LS <http://people.maths.ox.ac.uk/robertsl/dfols/userguide.html>`__ and
`DFO-GN <http://people.maths.ox.ac.uk/robertsl/dfogn/userguide.html>`__.
They are derivative free optimization solvers that were developed by Lindon Roberts at the University
of Oxford, in conjunction with NAG.  They are particularly well suited for solving noisy problems.

FitBenchmarking currently supports the DFO minimizers:

* `Derivative-Free Optimizer for Least Squares <http://people.maths.ox.ac.uk/robertsl/dfols/userguide.html>`__ (:code:`dfols`)

* `Derivative-Free Gauss-Newton Solver <http://people.maths.ox.ac.uk/robertsl/dfogn/userguide.html>`__ (:code:`dfogn`)

**Licence** Both `DFO-GN <https://github.com/numericalalgorithmsgroup/dfogn/blob/master/LICENSE.txt>`__ and `DFO-LS <https://github.com/numericalalgorithmsgroup/dfols/blob/master/LICENSE.txt>`__ are available under the GPL-3 licence.  A proprietary licence is also available from `NAG <https://www.nag.com/content/worldwide-contact-information>`__ . 
  
**Links** `GitHub - DFO-GN <https://github.com/numericalalgorithmsgroup/dfogn>`__ `GitHub - DFO-LS <https://github.com/numericalalgorithmsgroup/dfols>`__

The DFO minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    dfo: dfols
         dfogn

.. warning::
   Additional dependencies `DFO-GN` and `DFO-LS` must be installed for
   these to be available;
   See :ref:`extra_dependencies`.

.. _gradient-free:

Gradient-Free-Optimizers (``gradient_free``)
============================================

`Gradient-Free-Optimizers <https://github.com/SimonBlanke/Gradient-Free-Optimizers>`__ are a collection of
gradient-free methods capable of solving various optimization problems. Please note that Gradient-Free-Optimizers
must be run with problems that have finite bounds on all parameters.

*  Hill Climbing (:code:`HillClimbingOptimizer`)

*  Repulsing Hill Climbing (:code:`RepulsingHillClimbingOptimizer`)
                   
*  Simulated Annealing (:code:`SimulatedAnnealingOptimizer`)

*  Random Search (:code:`RandomSearchOptimizer`)
                   
*  Random Restart Hill Climbing (:code:`RandomRestartHillClimbingOptimizer`)
               
*  Random Annealing (:code:`RandomAnnealingOptimizer`)
   
*  Parallel Tempering (:code:`ParallelTemperingOptimizer`)
   
*  Particle Swarm (:code:`ParticleSwarmOptimizer`)
                   
*  Evolution Strategy (:code:`EvolutionStrategyOptimizer`)
                   
*  Bayesian (:code:`BayesianOptimizer`)

*  Tree Structured Parzen Estimators (:code:`TreeStructuredParzenEstimators`)
                   
*  Decision Tree (:code:`DecisionTreeOptimizer`)

**Licence** The Gradient-Free-Optimizers package is available under an `MIT Licence <https://github.com/SimonBlanke/Gradient-Free-Optimizers/blob/master/LICENSE>`__ .

   
The `gradient_free` minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    gradient_free: HillClimbingOptimizer
                   RepulsingHillClimbingOptimizer
                   SimulatedAnnealingOptimizer
                   RandomSearchOptimizer
                   RandomRestartHillClimbingOptimizer
                   RandomAnnealingOptimizer
                   ParallelTemperingOptimizer
                   ParticleSwarmOptimizer
                   EvolutionStrategyOptimizer
                   BayesianOptimizer
                   TreeStructuredParzenEstimators
                   DecisionTreeOptimizer

.. warning::
   The additional dependency Gradient-Free-Optimizers must be installed for this to be available;
   See :ref:`extra_dependencies`.

.. note::
   BayesianOptimizer, TreeStructuredParzenEstimators and DecisionTreeOptimizer may be slow running and
   so are not run by default when `gradient_free` software is selected. To run these minimizers you must
   explicity set them as seen above.

.. _gsl:
	 
GSL (``gsl``)
=============

The `GNU Scientific Library <https://www.gnu.org/software/gsl/>`__ is a numerical library that
provides a wide range of mathematical routines.  We call GSL using  the `pyGSL Python interface
<https://sourceforge.net/projects/pygsl/>`__.

The GSL routines have a number of parameters that need to be chosen, often without default suggestions.
We have taken the values as used by Mantid.

We provide implementations for the following
packages in the `multiminimize <https://www.gnu.org/software/gsl/doc/html/multimin.html>`__ and `multifit <https://www.gnu.org/software/gsl/doc/html/nls.html>`__ sections of the library:


* `Levenberg-Marquardt (unscaled) <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multifit__nlin.lmder>`__ (:code:`lmder`)

* `Levenberg-Marquardt (scaled) <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multifit_nlin.lmsder>`__ (:code:`lmsder`)

* `Nelder-Mead Simplex Algorithm <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multiminimize.nmsimplex>`__ (:code:`nmsimplex`)

* `Nelder-Mead Simplex Algorithm (version 2) <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multiminimize.nmsimplex2>`__ (:code:`nmsimplex2`)

* `Polak-Ribiere Conjugate Gradient Algorithm <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multiminimize.conjugate_pr>`__ (:code:`conjugate_pr`)

* `Fletcher-Reeves Conjugate-Gradient <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multiminimize.conjugate_fr>`__ (:code:`conjugate_fr`)

* `The vector quasi-Newton BFGS method <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multiminimize.vector_bfgs>`__ (:code:`vector_bfgs`)

* `The vector quasi-Newton BFGS method (version 2) <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multiminimize.vector_bfgs2>`__ (:code:`vector_bfgs2`)

* `Steepest Descent <http://pygsl.sourceforge.net/api/pygsl.html#pygsl.multiminimize.steepest_descent>`__ (:code:`steepest_descent`)

**Links** `SourceForge PyGSL <http://pygsl.sourceforge.net/>`__

**Licence** The GNU Scientific Library is available under the `GPL-3 licence <https://www.gnu.org/licenses/gpl-3.0.html>`__ .

The GSL minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    gsl: lmsder
         lmder
         nmsimplex
         nmsimplex2
         conjugate_pr
         conjugate_fr
         vector_bfgs
         vector_bfgs2
         steepest_descent
	 
.. warning::
   The external packages GSL and pygsl must be installed to use these minimizers.

.. _mantid:

Mantid (``mantid``)
===================

`Mantid <https://www.mantidproject.org>`__ is a framework created to
manipulate and analyze neutron scattering and muon spectroscopy data.
It has support for a number of minimizers, most of which are from GSL.

* `BFGS <https://docs.mantidproject.org/nightly/fitting/fitminimizers/BFGS.html>`__ (:code:`BFGS`)

* `Conjugate gradient (Fletcher-Reeves) <https://docs.mantidproject.org/nightly/fitting/fitminimizers/FletcherReeves.html>`__ (:code:`Conjugate gradient (Fletcher-Reeves imp.)`)

* `Conjugate gradient (Polak-Ribiere) <https://docs.mantidproject.org/nightly/fitting/fitminimizers/PolakRibiere.html>`__ (:code:`Conjugate gradient (Polak-Ribiere imp.)`)

* `Damped GaussNewton <https://docs.mantidproject.org/nightly/fitting/fitminimizers/DampedGaussNewton.html>`__ (:code:`Damped GaussNewton`)

* `FABADA <https://docs.mantidproject.org/nightly/concepts/FABADA.html>`__ (:code:`FABADA`)
 
* `Levenberg-Marquardt algorithm <https://docs.mantidproject.org/nightly/fitting/fitminimizers/LevenbergMarquardt.html>`__ (:code:`Levenberg-Marquardt`)

* `Levenberg-Marquardt MD <https://docs.mantidproject.org/nightly/fitting/fitminimizers/LevenbergMarquardtMD.html>`__ (:code:`Levenberg-MarquardtMD`) - An implementation of Levenberg-Marquardt intended for MD workspaces, where work is divided into chunks to achieve a greater efficiency for a large number of data points.

* `Simplex <https://docs.mantidproject.org/nightly/fitting/fitminimizers/Simplex.html>`__ (:code:`Simplex`)

* `SteepestDescent <https://docs.mantidproject.org/nightly/fitting/fitminimizers/GradientDescent.html>`__ (:code:`SteepestDescent`)

* `Trust Region <https://docs.mantidproject.org/nightly/fitting/fitminimizers/TrustRegion.html>`__ (:code:`Trust Region`) - An implementation of one of the algorithms available in RALFit.

 **Links** `GitHub - Mantid <https://github.com/mantidproject/mantid>`__ `Mantid's Fitting Docs <https://docs.mantidproject.org/nightly/algorithms/Fit-v1.html>`__

**Licence** Mantid is available under the `GPL-3 licence <https://github.com/mantidproject/mantid/blob/master/LICENSE.txt>`__ .

 
The Mantid minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    mantid: BFGS
            Conjugate gradient (Fletcher-Reeves imp.)
            Conjugate gradient (Polak-Ribiere imp.)
            Damped GaussNewton
	    FABADA
            Levenberg-Marquardt
            Levenberg-MarquardtMD
            Simplex
            SteepestDescent
            Trust Region

.. warning::
   The external package Mantid must be installed to use these minimizers.

.. _levmar:

Levmar (``levmar``)
===================

The `levmar <http://users.ics.forth.gr/~lourakis/levmar/>`__ package
which implements the Levenberg-Marquardt method for nonlinear least-squares.
We interface via the python interface `available on PyPI <https://pypi.org/project/levmar/>`__.

* Levenberg-Marquardt with supplied Jacobian (:code:`levmar`)  - the Levenberg-Marquardt method

**Licence** Levmar is available under the `GPL-3 licence <http://www.gnu.org/copyleft/gpl.html>`__ .  A paid licence for proprietary commerical use is `available from the author <http://users.ics.forth.gr/~lourakis/levmar/faq.html#Q37>`__ .
  
The `levmar` minimizer is set as follows:

.. code-block:: rst

   [MINIMIZERS]
   levmar: levmar


.. warning::
   The additional dependency levmar must be installed for this to be available;
   See :ref:`extra_dependencies`. This package also requires the BLAS and LAPACK
   libraries to be present on the system.

.. _matlab:

Matlab (``matlab``)
===================

We call the `fminsearch <https://uk.mathworks.com/help/matlab/ref/fminsearch.html>`__
function from `MATLAB <https://uk.mathworks.com/products/matlab.html>`__, using the
MATLAB Engine API for Python.

* Nelder-Mead Simplex (:code:`Nelder-Mead Simplex`)

**Licence** Matlab is a `proprietary product <https://www.mathworks.com/pricing-licensing.html>`__ .
  
The `matlab` minimizer is set as follows:

.. code-block:: rst

   [MINIMIZERS]
   matlab: Nelder-Mead Simplex

.. warning::
   MATLAB must be installed for this to be available; See :ref:`external-instructions`.

.. _matlab-curve:

Matlab Curve Fitting Toolbox (``matlab_curve``)
===============================================

We call the `fit <https://uk.mathworks.com/help/curvefit/fit.html>`_
function from the `MATLAB Curve Fitting Toolbox <https://uk.mathworks.com/help/curvefit/index.html>`_,
using the MATLAB Engine API for Python.

* Levenberg-Marquardt (:code:`Levenberg-Marquardt`)
* Trust-Region (:code:`Trust-Region`)

**Licence** Matlab and the Curve Fitting Toolbox are both `proprietary products <https://www.mathworks.com/pricing-licensing.html>`__ .
  
The `matlab_curve` minimizers are set as follows:

.. code-block:: rst

   [MINIMIZERS]
   matlab_curve: Levenberg-Marquardt
                 Trust-Region

.. warning::
   MATLAB Curve Fitting Toolbox must be installed for this to be available; See :ref:`external-instructions`.

.. _matlab-opt:

Matlab Optimization Toolbox (``matlab_opt``)
============================================

We call the `lsqcurvefit <https://uk.mathworks.com/help/optim/ug/lsqcurvefit.html>`__
function from the `MATLAB Optimization Toolbox <https://uk.mathworks.com/products/optimization.html>`__,
using the MATLAB Engine API for Python.

* Levenberg-Marquardt (:code:`levenberg-marquardt`)
* Trust-Region-Reflective (:code:`trust-region-reflective`)

**Licence** Matlab and the Optimization Toolbox are both `proprietary products <https://www.mathworks.com/pricing-licensing.html>`__ .
  
The `matlab_opt` minimizers are set as follows:

.. code-block:: rst

   [MINIMIZERS]
   matlab_opt: levenberg-marquardt
               trust-region-reflective

.. warning::
   MATLAB Optimization Toolbox must be installed for this to be available; See :ref:`external-instructions`.

.. _matlab-stats:

Matlab Statistics Toolbox (``matlab_stats``)
============================================


We call the `nlinfit <https://uk.mathworks.com/help/stats/nlinfit.html>`__
function from the `MATLAB Statistics Toolbox <https://uk.mathworks.com/products/statistics.html>`__,
using the MATLAB Engine API for Python.

* Levenberg-Marquardt (:code:`Levenberg-Marquardt`)

**Licence** Matlab and the Statistics Toolbox are both `proprietary products <https://www.mathworks.com/pricing-licensing.html>`__ .
  
The `matlab_stats` minimizer is set as follows:

.. code-block:: rst
  
  [MINIMIZERS]
  matlab_stats: Levenberg-Marquardt

.. warning::
   MATLAB Statistics Toolbox must be installed for this to be available; See :ref:`external-instructions`.

.. _minuit:
	   
Minuit (``minuit``)
===================

CERN developed the `Minuit 2 <https://root.cern.ch/doc/master/Minuit2Page.html>`__ package
to find the minimum value of a multi-parameter function, and also to compute the
uncertainties.
We interface via the python interface `iminuit <https://iminuit.readthedocs.io>`__ with
support for the 2.x series. 

* `Minuit's MIGRAD <https://root.cern.ch/root/htmldoc/guides/minuit2/Minuit2.pdf>`__ (:code:`minuit`)

**Links** `Github - iminuit <https://github.com/scikit-hep/iminuit>`__

**Licence** iminuit is released under the `MIT licence <https://github.com/scikit-hep/iminuit/blob/develop/LICENSE>`__, while Minuit 2 is `LGPL v2 <https://github.com/root-project/root/blob/master/LICENSE>`__ .

The Minuit minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    minuit: minuit

.. warning::
   The additional dependency Minuit must be installed for this to be available;
   See :ref:`extra_dependencies`.	 

.. _ralfit:
    
RALFit (``ralfit``)
===================

`RALFit <https://ralfit.readthedocs.io/projects/Fortran/en/latest/>`__
is a nonlinear least-squares solver, the development of which was funded
by the EPSRC grant `Least-Squares: Fit for the Future`.  RALFit is designed to be able
to take advantage of higher order derivatives, although only first
order derivatives are currently utilized in FitBenchmarking.

* Gauss-Newton, trust region method (:code:`gn`)
* Hybrid Newton/Gauss-Newton, trust region method (:code:`hybrid`)
* Gauss-Newton, regularization (:code:`gn_reg`)
* Hybrid Newton/Gauss-Newton, regularization (:code:`hybrid_reg`)

**Links** `Github - RALFit <https://github.com/ralna/ralfit/>`__. RALFit's Documentation on: `Gauss-Newton/Hybrid models <https://ralfit.readthedocs.io/projects/Fortran/en/latest/method.html#the-models>`__,  `the trust region method <https://ralfit.readthedocs.io/projects/Fortran/en/latest/method.html#the-trust-region-method>`__ and  `The regularization method <https://ralfit.readthedocs.io/projects/C/en/latest/method.html#regularization>`__

**Licence** RALFit is available under a `3-clause BSD Licence <https://github.com/ralna/RALFit/blob/master/LICENCE>`__

The RALFit minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    ralfit: gn
            gn_reg
            hybrid
            hybrid_reg

.. warning::
   The external package RALFit must be installed to use these minimizers.

.. _scipy:

SciPy (``scipy``)
=================

`SciPy <https://www.scipy.org>`__ is the standard python package for mathematical
software.  In particular, we use the `minimize <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`__
solver for general minimization problems from the optimization chapter of
SciPy's library. Currently we only use the algorithms that do not require
Hessian information as inputs.

* `Nelder-Mead algorithm <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html>`__ (:code:`Nelder-Mead`)
* `Powell algorithm <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-powell.html>`__ (:code:`Powell`)
* `Conjugate gradient algorithm <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-cg.html>`__ (:code:`CG`)
* `BFGS algorithm <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-bfgs.html>`__ (:code:`BFGS`)
* `Newton-CG algorithm <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-newtoncg.html>`__  (:code:`Newton-CG`)
* `L-BFGS-B algorithm <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html>`__ (:code:`L-BFGS-B`)
* `Truncated Newton (TNC) algorithm <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-tnc.html>`__ (:code:`TNC`)
* `Sequential Least SQuares Programming <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html>`__ (:code:`SLSQP`)

**Links** `Github - SciPy minimize <https://github.com/scipy/scipy/blob/master/scipy/optimize/_minimize.py>`__

**Licence** Scipy is available under a `3-clause BSD Licence <https://github.com/scipy/scipy/blob/master/LICENSE.txt>`__.  Individual packages may have their own (compatible) licences, as listed `here <https://github.com/scipy/scipy/blob/master/LICENSES_bundled.txt>`__.

The SciPy minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    scipy: Nelder-Mead
           Powell
           CG
           BFGS
           Newton-CG
           L-BFGS-B
           TNC
           SLSQP

.. _scipy-ls:

SciPy LS (``scipy_ls``)
=======================

`SciPy <https://www.scipy.org>`__ is the standard python package for mathematical
software.  In particular, we use the `least_squares <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html#scipy.optimize.least_squares>`__
solver for Least-Squares minimization problems from the optimization chapter
of SciPy's library.

* Levenberg-Marquardt with supplied Jacobian (:code:`lm-scipy`)  - a wrapper around MINPACK
* The Trust Region Reflective algorithm (:code:`trf`)
* A dogleg algorithm with rectangular trust regions (:code:`dogbox`)

**Links** `Github - SciPy least_squares <https://github.com/scipy/scipy/blob/master/scipy/optimize/_lsq/least_squares.py>`__

**Licence** Scipy is available under a `3-clause BSD Licence <https://github.com/scipy/scipy/blob/master/LICENSE.txt>`__.  Individual packages many have their own (compatible) licences, as listed `here <https://github.com/scipy/scipy/blob/master/LICENSES_bundled.txt>`__.

The SciPy least squares minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    scipy_ls: lm-scipy
              trf
              dogbox

.. _scipy-go:

SciPy GO (``scipy_go``)
=======================

`SciPy <https://www.scipy.org>`__ is the standard python package for mathematical
software.  In particular, we use the `Global Optimization <https://docs.scipy.org/doc/scipy/reference/optimize.html#global-optimization>`__
solvers for global optimization problems from the optimization chapter
of SciPy's library.

* `Differential Evolution (derivative-free) <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html#scipy.optimize.differential_evolution>`__ (:code:`differential_evolution`)
* `Simplicial Homology Global Optimization (SHGO) <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html#scipy.optimize.shgo>`__ (:code:`shgo`)
* `Dual Annealing <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.dual_annealing.html#scipy.optimize.dual_annealing>`__ (:code:`dual_annealing`)

**Links** `Github - SciPy optimization <https://github.com/scipy/scipy/blob/master/scipy/optimize/>`__

**Licence** Scipy is available under a `3-clause BSD Licence <https://github.com/scipy/scipy/blob/master/LICENSE.txt>`__.  Individual packages may have their own (compatible) licences, as listed `here <https://github.com/scipy/scipy/blob/master/LICENSES_bundled.txt>`__.

The SciPy global optimization minimizers are set as follows:

.. code-block:: rst

    [MINIMIZERS]
    scipy_go: differential_evolution
              shgo
              dual_annealing

.. note::
   The shgo solver is particularly slow running and should generally be avoided. As a result, this solver is
   not run by default when `scipy_go` software is selected. In order to run this minimizer, you must explicitly
   set it as above.
.. _fitting_option:

###############
Fitting Options
###############

Options that control the benchmarking process are set here.


Software (:code:`software`)
---------------------------

Software is used to select the fitting software to benchmark, this should be
a newline-separated list. Available options are:

* ``bumps`` (default software)
* ``dfo`` (default software)
* ``gradient_free`` (default software)
* ``gsl`` (external software -- see :ref:`external-instructions`)
* ``levmar`` (external software -- see :ref:`extra_dependencies`)
* ``mantid`` (external software -- see :ref:`external-instructions`)
* ``matlab`` (external software -- see :ref:`external-instructions`)
* ``matlab_curve`` (external software -- see :ref:`external-instructions`)
* ``matlab_opt`` (external software -- see :ref:`external-instructions`)
* ``matlab_stats`` (external software -- see :ref:`external-instructions`)
* ``minuit`` (default software)
* ``ralfit`` (external software -- see :ref:`external-instructions`)
* ``scipy`` (default software)
* ``scipy_ls`` (default software)
* ``scipy_go`` (default software)


Default are ``bumps``, ``dfo``, ``gradient_free``, ``minuit``, ``scipy``, ``scipy_ls`` and ``scipy_go``

.. code-block:: rst

    [FITTING]
    software: bumps
              dfo
              minuit
              scipy
              scipy_ls
              scipy_go

.. warning::

   Software must be listed to be here to be run.
   Any minimizers set in :ref:`minimizer_option` will not be run if the software is not also
   present in this list.


Number of minimizer runs (:code:`num_runs`)
-------------------------------------------

Sets the number of runs to average each fit over.

Default is ``5``

.. code-block:: rst

    [FITTING]
    num_runs: 5

.. _algorithm_type:

Algorithm type (:code:`algorithm_type`)
---------------------------------------

This is used to select what type of algorithm is used within a specific software.
For a full list of available minimizers for each algorithm type, see :ref:`minimizer_types`.
The options are:

* ``all`` - all minimizers
* ``ls`` - least-squares fitting algorithms
* ``deriv_free`` - derivative free algorithms (these are algorithms that cannot use
  information about derivatives -- e.g., the ``Simplex`` method in ``Mantid``),
  see :ref:`deriv_free`.
* ``general`` - minimizers which solve a generic `min f(x)`
* ``simplex`` - derivative free simplex based algorithms e.g. Nelder-Mead, see :ref:`Simplex <simplex>`
* ``trust_region`` - algorithms which employ a trust region approach,  see :ref:`trust_region`
* ``levenberg-marquardt`` - minimizers that use the Levenberg Marquardt algorithm, see :ref:`Levenberg-Marquardt <levenberg_marquardt>`.
* ``gauss_newton`` - minimizers that use the Gauss Newton algorithm, see :ref:`Gauss-Newton <gauss_newton>`
* ``bfgs`` - minimizers that use the BFGS algorithm, see :ref:`BFGS <bfgs>`
* ``conjugate_gradient`` - Conjugate Gradient algorithms, see :ref:`Conjugate Gradient <conjugate_gradient>`
* ``steepest_descent`` - Steepest Descent algorithms, see :ref:`Steepest Descent <steepest_descent>`
* ``global_optimization`` - Global Optimization algorithms

Default is ``all``

.. code-block:: rst

    [FITTING]
    algorithm_type: all

.. warning::

   Choosing an option other than ``all`` may deselect certain
   minimizers set in the options file


Jacobian method (:code:`jac_method`)
------------------------------------

This sets the Jacobian used.
Choosing multiple options via a new line seperated list will result in all
combinations being benchmarked.
Current Jacobian methods are:

* ``analytic`` - uses the analytic Jacobian extracted from the fitting problem.
* ``scipy`` -  uses :ref:`SciPy's finite difference Jacobian approximations <scipy-jac>`.
* ``default`` - uses the default derivative approximation implemented in the minimizer.
* ``numdifftools`` - uses the python package :ref:`numdifftools <numdifftools-jac>`.

Default is ``default``

.. code-block:: rst

    [FITTING]
    jac_method: scipy

.. warning::

   Currently analytic Jacobians are only available for
   problems that use the cutest and NIST parsers.


Hessian method (:code:`hes_method`)
------------------------------------

This sets the Hessian used.
Choosing multiple options via a new line seperated list will result in all
combinations being benchmarked.
Current Hessian methods are:

* ``default`` - Hessian information is not passed to minimizers
* ``analytic`` - uses the analytic Hessian extracted from the fitting problem.
* ``scipy`` -  uses :ref:`SciPy's finite difference approximations <scipy-hes>`.
* ``numdifftools`` - uses the python package :ref:`numdifftools <numdifftools-hes>`.

Default is ``default``

.. code-block:: rst

    [FITTING]
    hes_method: default

.. warning::

   Currently analytic Hessians are only available for
   problems that use the cutest and NIST parsers.

Cost function (:code:`cost_func_type`)
--------------------------------------

This sets the cost functions to be used for the given data.
Choosing multiple options via a new line seperated list will result in all
combinations being benchmarked.
Currently supported cost functions are:

* ``nlls`` - This sets the cost function to be non-weighted non-linear least squares, :class:`~fitbenchmarking.cost_func.nlls_cost_func.NLLSCostFunc`.

* ``weighted_nlls`` - This sets the cost function to be weighted non-linear least squares, :class:`~fitbenchmarking.cost_func.weighted_nlls_cost_func.WeightedNLLSCostFunc`.

* ``hellinger_nlls`` - This sets the cost function to be the Hellinger cost function, :class:`~fitbenchmarking.cost_func.hellinger_nlls_cost_func.HellingerNLLSCostFunc`.

* ``poisson`` - This sets the cost function to be the Poisson Deviation cost function, :class:`~fitbenchmarking.cost_func.poisson_cost_func.PoissonCostFunc`.


Default is ``weighted_nlls``

.. code-block:: rst

    [FITTING]
    cost_func_type: weighted_nlls

Maximum Runtime (:code:`max_runtime`)
--------------------------------------

This sets the maximum runtime a minimizer has to solve one benchmark
problem `num_runs` number of times, where `num_runs` is another option a
user can set. If the minimizer is still running after the maximum time
has elapsed, then this result will be skipped and FitBenchmarking will move
on to the next minimizer / benchmark dataset combination. The main purpose
of this option is to get to result tables quicker by limit the runtime.

`max_runtime` is set by specifying a number in unit of seconds. Please note
that depending on platform the time specified with `max_runtime` may not
match entirely with the absolute run-times specified in tables. Hence you
may have to experiment a bit with this option to get the cutoff you want.

Default is 600 seconds

.. code-block:: rst

    [FITTING]
    max_runtime: 600
.. _options:

#######################
FitBenchmarking Options
#######################

The default behaviour of FitBenchmarking can be changed by
supplying an options file.  The default values of these options,
and how to override them, are given in the pages below.

.. toctree::
    :titlesonly:
    :maxdepth: 2
    :caption: Contents:

    fitting_option
    minimizer_option
    jacobian_option
    hessian_option
    plotting_option
    output_option
    logging_option


The options file must be a ``.ini`` formatted file
(`see here <https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`__).
Some example files can be found in the ``examples`` folder of the source, which is also
available to download at :ref:`BenchmarkProblems`.


.. _logging_option:

###############
Logging Options
###############

The logging section contains options to control how fitbenchmarking logs
information.

Logging file name (:code:`file_name`)
-------------------------------------

This specifies the file path to write the logs to.

Default is ``fitbenchmarking.log``

.. code-block:: rst

    [LOGGING]
    file_name: fitbenchmarking.log

Logging append (:code:`append`)
-------------------------------

This specifies whether to log in append mode or not.
If append mode is active, the log file will be extended with each subsequent
run, otherwise the log will be cleared after each run.

Default is ``False`` (``yes``/``no`` can also be used)

.. code-block:: rst

    [LOGGING]
    append: no


Logging level (:code:`level`)
-----------------------------------------

This specifies the minimum level of logging to display on console during
runtime.
Options are (from most logging to least):

* ``NOTSET``
* ``DEBUG``
* ``INFO``
* ``WARNING``
* ``ERROR``
* ``CRITICAL``

Default is ``INFO``

.. code-block:: rst

    [LOGGING]
    level: INFO

Logging external output (:code:`external_output`)
-------------------------------------------------

This selects the amount of information displayed from third-parties.
There are 3 options:

* ``display``: Print information from third-parties to the stdout stream during a run.
* ``log_only``: Print information to the log file but not the stdout stream.
* ``debug``: Do not intercept third-party use of output streams.

Default is ``log_only``

.. code-block:: rst

    [LOGGING]
    append: log_only
.. _multifit:

************************************
Native File Format (Mantid MultiFit)
************************************

As part of the Mantid parsing we also offer limited support for Mantid's
`MultiFit <https://docs.mantidproject.org/nightly/algorithms/Fit-v1.html?highlight=fit#multiple-fit>`__
functionality.

Here we outline how to use Mantid's MultiFit with FitBenchmarking,
in which some options differ from the standard  :ref:`native`.

.. warning::
   Due to the way Mantid uses ties (a central feature of MultiFit),
   MultiFit problems can only be used with Mantid minimizers.

In this format, data is separated from the function. This allows running the
same dataset against multiple different models to assess which is the most
appropriate.

An example of a multifit problem is:

.. literalinclude:: ../../../../examples/benchmark_problems/MultiFit/MUSR62260.txt


Below we outline the differences between this and the :ref:`native`.

software
  Must be `Mantid`.

name
  As in :ref:`native`.

description
  As in :ref:`native`.

input_file
  As in :ref:`native`, but you must pass in a list of
  data files (see above example).

function
  As in :ref:`native`.
  
  When fitting, this function will be used for each of the ``input_files`` given
  simultaneously.

ties
  This entry is used to define global variables by tieing a variable across
  input files.

  Each string in the list should reference a parameter in the function using
  Mantid's convention of ``f<i>.<name>`` where ``i`` is the position of the
  function in the function string, and ``name`` is the global parameter.

  For example to run a fit which has a shared background and peak height,
  the function and ties fields might look like::

     function='name=LinearBackground, A0=0, A1=0; name=Gaussian, Height=0.01, PeakCentre=0.00037, Sigma=1e-05'
     ties=['f0.A0', 'f0.A1', 'f1.Height']

fit_ranges
  As in :ref:`native`.
====================
 CUTEst File Format
====================

The CUTEst file format in FitBenchmarking is a slight modification of the
`SIF format <http://www.numerical.rl.ac.uk/lancelot/sif/sif.html>`_.
Specifically, the data points, errors, and the number of variables
must be defined in such a way to allow FitBenchmarking to access this data; see below.
In FitBenchmarking, all SIF files are assumed to be CUTEst problems.

These problems are a subset of the problems in the
`CUTEr/st Test Problem Set <http://www.cuter.rl.ac.uk/Problems/mastsif.shtml>`_,
which may have been adapted to work with FitBenchmarking.

The SIF file format is very powerful, and CUTEst will work with arbitrary
variable names, however for FitBenchmarking, these must match a set of expected
variable names.

**Licence** This file format needs PyCUTEst and the packages ARCHDefs, CUTEst and
SIFDECODE to be installed.
PyCUTEst is available under the
`GPL-3 <https://github.com/jfowkes/pycutest/blob/master/LICENSE>`__ licence.
`ARCHDEFS <https://github.com/ralna/ARCHDefs/blob/master/LICENSE>`__,
`CUTEst <https://github.com/ralna/CUTEst/blob/master/LICENSE>`__ and
`SIFDECODE <https://github.com/ralna/SIFDecode/blob/master/LICENSE>`__
are available under an LGPL (v2.1 or later) licence.

Modifications to the SIF format for FitBenchmarking problems
============================================================

In order for FitBenchmarking to access the data, the SIF files must
be written using the following conventions.

Defining Data
-------------

Data should be defined using the format::

     RE X<idx>        <val_x>
     RE Y<idx>        <val_y>
     RE E<idx>        <val_error>

where ``<idx>`` is the index of the data point, and ``<val_x>``, ``<val_y>``,
and ``<val_error>`` are the values attributed to it.

Usually, ``<idx>`` will range from 1 to ``<num_x>``, with that defined as::

     IE M             <num_x>

If ``<idx>`` does not start at 1, the following lines can be used to specify
the range::

     IE MLOWER        <min_idx>
     IE MUPPER        <max_idx>

Defining Variables
------------------

For the free variables in functions, we use the convention::

     IE N             <num_vars>

This is used to tell FitBenchmarking how many degrees of freedom we need to
fit.
In some cases variables will be vectors, and the number of degrees of freedom
will be greater, most problems use ``NVEC`` as a convention to input the number
of vectors.

Support for Bounds
==================

Parameter ranges can be added to SIF files using the `BOUNDS <https://www.numerical.rl.ac.uk/lancelot/sif/node26.html>`_
indicator card.

Currently in Fitbenchmarking, problems with parameter ranges can be handled by SciPy, Bumps, Minuit, Mantid,
Matlab Optimization Toolbox, DFO, Levmar and RALFit fitting software. Please note that the following Mantid
minimizers currently throw an exception when parameter ranges are used: BFGS, Conjugate gradient
(Fletcher-Reeves imp.), Conjugate gradient (Polak-Ribiere imp.) and SteepestDescent.
.. _native:

******************
Native File Format
******************

In FitBenchmarking, the native file format is used to read IVP, Mantid, and
SASView problems.

In this format, data is separated from the function. This allows running the
same dataset against multiple different models to assess which is the most
appropriate.

Examples of native problems are:

.. literalinclude:: ../../../../examples/benchmark_problems/Data_Assimilation/lorentz.txt

.. literalinclude:: ../../../../examples/benchmark_problems/Muon/Muon_HIFI_113856.txt

.. literalinclude:: ../../../../examples/benchmark_problems/SAS_modelling/SASView_Simple_Shapes_1D/1D_cylinder_neutron_def0.txt

These examples show the basic structure in which the file starts with a comment
indicating it is a FitBenchmark problem followed by key-value pairs. Available
keys are described below:

software
  Either 'IVP', 'Mantid', or 'SasView' (case insensitive).
  
  This defines whether to use an IVP format, Mantid, or SasView to generate the model.
  The 'Mantid' software also supports Mantid's MultiFit functionality, which
  requires the parameters listed here to be defined slightly differently.
  More information can be found in :ref:`multifit`.

  **Licence** Mantid is available under a
  `GPL-3 Licence <https://github.com/mantidproject/mantid/blob/master/LICENSE.txt>`_.
  The component of SasView we use is SasModels, which is available under a
  `BSD 3-clause <https://github.com/SasView/sasmodels/blob/master/LICENSE.txt>`_ licence.

name
  The name of the problem.

  This will be used as a unique reference so should not match other names in the
  dataset. A sanitised version of this name will also be used in filenames with
  commas stripped out and spaces replaced by underscores.

description
  A description of the dataset.

  This will be displayed in the Problem Summary Pages and Fitting Reports produced by a
  benchmark.

input_file
  The name of a file containing the data to fit.

  The file must be in a subdirectory named ``data_files``, and should have the form::

     header

     x11 [x12 [x13 ...]] y11 [y12 [y13 ...]] [e11 [e12 ...]]
     x21 [x22 [x23 ...]] y21 [y22 [y23 ...]] [e21 [e22 ...]]
     ...

  Mantid uses the convention of ``# X Y E`` as the header and SASView uses
  the convention ``<X>   <Y>   <E>``, although neither of these are enforced.
  The error column is optional in this format.

  If the data contains multiple inputs or outputs, the header must be written
  in one of the above conventions with the labels as "x", "y", or "e" followed by
  a number. An example of this can be seen in
  ``examples/benchmark_problems/Data_Assimilation/data_files/lorentz.txt``

function
  This defines the function that will be used as a model for the fitting.

  Inside FitBenchmarking, this is passed on to the specified software and, as
  such, the format is specific to the package we wish to use, as described below.

  **IVP**

  The IVP parser allows a user to define ``f`` in the following equation:

  .. math:: x' = f(t, x, *args)

  To do this we use a python module to define the function. As in the above
  formula, the function can take the following arguments:

  - *t* (float): The time to evaluate at
  - *x* (np.array): A value for x to evaluate at
  - *\*args* (floats): The parameters to fit

  To link to this function we use a function string with the following
  parameters:

  - *module*: The path to the module
  - *func*: The name of the function within the module
  - *step*: The time step that the input data uses
    (currently only fixed steps are supported - if you need
    varying time steps please raise an issue on our GitHub)
  - *\*args*: Starting values for the parameters

  **Mantid**

  A Mantid function consists of one or more base functions separated by a semicolon.
  This allows for a powerful way of describing problems, which may have multiple
  components such as more than one Gaussian and a linear background.

  To use one of the base functions in Mantid, please see the list available
  `here <https://docs.mantidproject.org/nightly/fitting/fitfunctions/categories/FitFunctions.html>`__.

  *Note: Any non-standard arguments (e.g. ties, constraints, fixes, ...) will
  only work with Mantid fitting software. Using other minimizers to fit these
  problems will result in the non-standard arguments being ignored.*

  **SASView**

  SASView functions can be any of
  `these <http://www.sasview.org/docs/user/qtgui/Perspectives/Fitting/models/index.html>`__.

fit_ranges
  This specifies the region to be fit.

  It takes the form shown in the example, where the first number
  is the minimum in the range and the second is the maximum.

parameter_ranges
  An optional setting which specifies upper and lower bounds for 
  parameters in the problem.

  Similarly to ``fit_ranges``, it takes the form where the first number
  is the minimum in the range and the second is the maximum.

  Currently in Fitbenchmarking, problems with `parameter_ranges` can
  be handled by SciPy, Bumps, Minuit, Mantid, Matlab Optimization Toolbox,
  DFO, Levmar and RALFit fitting software. Please note that the following
  Mantid minimizers currently throw an exception when `parameter_ranges`
  are used: BFGS, Conjugate gradient (Fletcher-Reeves imp.),
  Conjugate gradient (Polak-Ribiere imp.) and SteepestDescent.
.. _problem_def:

************************
Problem Definition Files
************************

In FitBenchmarking, problems can be defined using several file formats.
The ``examples/benchmark_problems`` directory holds a collection of these
that can be used for reference.

More information on the supported formats can be found on the following pages.

.. toctree::
    :titlesonly:
    :maxdepth: 2
    :caption: File Formats:

    cutest
    native
    multifit
    nist

Detecting problem file type
===========================

FitBenchmarking detects which parser to use in two ways:

    - For the CUTEst file format we check that the extension of the data file is `sif`
    - For native and NIST file formats we check the first line of the file

        - `# FitBenchmark Problem` corresponds to the native format
        - `# NIST/ITL StRD` corresponds to the NIST format
***********
NIST Format
***********

The NIST file format is based on the `nonlinear regression <https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml>`__ problems found at the `NIST Standard Reference Database <https://www.itl.nist.gov/div898/strd/>`__. Documentation and background of these problems can be found `here <https://www.itl.nist.gov/div898/strd/general/bkground.html>`__.

We note that FitBenchmarking recognizes the NIST file type by checking the first line of the file starts with `# NIST/ITL StRD`.
.. _getting-started:

#######################################################
Installing FitBenchmarking and default fitting packages
#######################################################

We recommend using |Python 3.7.1+| for running/installing Fitbenchmarking.
The easiest way to install FitBenchmarking is by using the Python package manager,
`pip <https://pip.pypa.io/en/stable/>`__.


Installing via ``pip``
----------------------

FitBenchmarking can be installed via the command line by entering:

.. code-block:: bash

      python -m pip install fitbenchmarking[bumps,DFO,gradient_free,minuit,SAS,numdifftools]


This will install the latest stable version of FitBenchmarking.
For all available versions please visit the FitBenchmarking
`PyPI project <https://pypi.org/project/fitbenchmarking/>`__.
FitBenchmarking can also use additional software that cannot be installed
using pip; please see :ref:`external-instructions` for details.

.. note::

    This install will include additional optional packages --
    see :ref:`extra_dependencies`.
    Any of the dependencies in the square brackets can be omitted, if required,
    and that package will not be available for Benchmarking, or will use the
    version of the package already on your system, if appropriate.

.. _installing_from_source:

Installing from source
----------------------

You may instead wish to install from source, e.g., to get the very latest version
of the code that is still in development.

1. Download this repository or clone it using
   `git <https://git-scm.com/>`__:
   ``git clone https://github.com/fitbenchmarking/fitbenchmarking.git``
2. Open up a terminal (command prompt) and go into the
   ``fitbenchmarking`` directory.
3. Once you are in the right directory, we recommend that you type

   .. code-block:: bash

      python -m pip install .[bumps,DFO,gradient_free,minuit,SAS,numdifftools]

4. Additional software that cannot be installed via pip can also be used
   with FitBenchmarking.  Follow the instructions at
   :ref:`external-instructions`.

.. _extra_dependencies:

Extra dependencies
------------------

In addition to the external packages described at :ref:`external-instructions`,
some optional dependencies can be installed directly by FitBenchmarking.
These are installed by issuing the commands

.. code-block:: bash

   python -m pip install fitbenchmarking['option-1','option-2',...]

or

.. code-block:: bash

   python -m pip install .['option-1','option-2',...]

where valid strings ``option-x`` are:

* ``bumps``-- installs the `Bumps <https://bumps.readthedocs.io>`_ fitting package.
* ``DFO`` -- installs the `DFO-LS <http://people.maths.ox.ac.uk/robertsl/dfols/userguide.html>`_ and `DFO-GN <http://people.maths.ox.ac.uk/robertsl/dfogn/userguide.html>`_ fitting packages.
* ``gradient_free`` -- installs the `Gradient-Free-Optimizers <https://github.com/SimonBlanke/Gradient-Free-Optimizers>`_ fitting package 
* ``levmar`` -- installs the `levmar <http://users.ics.forth.gr/~lourakis/levmar/>`_ fitting package.  Note that the interface we use also requires BLAS and LAPLACK to be installed on the system, and calls to this minimizer will fail if these libraries are not present.
* ``mantid`` -- installs the `h5py <https://pypi.org/project/h5py/>`_ and `pyyaml <https://pypi.org/project/PyYAML/>`_ modules.
* ``matlab`` -- installs the `dill <https://pypi.org/project/dill/>`_ module required to run matlab controllers in fitbenchmarking
* ``minuit`` -- installs the `Minuit <http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/>`_ fitting package.
* ``SAS`` -- installs the `Sasmodels <https://github.com/SasView/sasmodels>`_ fitting package and the `tinycc <https://pypi.org/project/tinycc/>`_ module.
* ``numdifftools`` -- installs the `numdifftools <https://numdifftools.readthedocs.io/en/latest/index.html>`_ numerical differentiation package.


.. |Python 3.7.1+| image:: https://img.shields.io/badge/python-3.7.1+-blue.svg
   :alt: Python 3.7.1+
   :target: https://www.python.org/downloads/

.. _external-instructions:

############################
Installing External Software
############################

Fitbenchmarking will install all packages that are available through pip.

To enable Fitbenchmarking with the other supported software,
they need to be installed and available on your machine.  We give
pointers outlining how to do this below, and you can find install scripts
for Ubuntu 18.04 in the directory `/build/<software>/`

CUTEst
------

CUTEst is used to parse SIF files in FitBenchmarking, and is called via the
PyCUTEst interface.

Currently this is only supported for Mac and Linux, and can be installed by
following the instructions outlined on the `pycutest documentation <https://jfowkes.github.io/pycutest/_build/html/install.html>`_

Please note that the ``PYCUTEST_CACHE`` environment variable must be set, and it must be
in the ``PYTHONPATH``.

GSL
---

GSL is used as a fitting software in FitBenchmarking, and is called via the
pyGSL interface.

Install instructions can be found at the `pyGSL docs <http://pygsl.sourceforge.net/>`__.
This package is also installable via pip, provided GSL is available on your system;
see our example build script in `build/gsl`.

Note: pyGSL may not be installable with the latest versions of pip. We have found that 20.0.2 works for our tests.

Mantid
------

Mantid is used both as fitting software, and to parse data files.

Instructions on how to install Mantid for a range of systems are available
at `<https://download.mantidproject.org/>`_.

MATLAB
------

MATLAB is available to use as fitting software in FitBenchmarking, and is
called via the MATLAB Engine API for Python.

To use this fitting software, both MATLAB and the MATLAB engine must be
installed. Installation instructions for MATLAB are available at
`<https://uk.mathworks.com/help/install/ug/install-products-with-internet-connection.html>`_,
and instructions for installing and setting up the MATLAB engine are
here: `<https://uk.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html>`_

RALFit
------

RALFit is available to use as fitting software.

Instructions on how to build the python interface are at `<https://ralfit.readthedocs.io/projects/Python/en/latest/install.html>`_

.. _install_instructions:

############
Installation
############

Fitbenchmarking will install all packages that can be installed through pip.
This includes the minimizers from SciPy, bumps, DFO-GN/LS, Minuit, and
also the SASModels package.

To enable Fitbenchmarking with the other supported software,
you must install them with the external software instructions.


.. toctree::
    :titlesonly:
    :maxdepth: 2
    :caption: Software:

    fitbenchmarking
    externals
