[![Build status](https://github.com/kutaslab/fitgrid/actions/workflows/fitgrid-cid.yml/badge.svg)](https://github.com/kutaslab/fitgrid/actions)
[![Coverage](https://codecov.io/gh/kutaslab/fitgrid/branch/main/graph/badge.svg)](https://codecov.io/gh/kutaslab/fitgrid)
[![DOI](https://zenodo.org/badge/147436563.svg)](https://zenodo.org/badge/latestdoi/147436563)

# fitgrid

A Python library for regression modeling time-varying patterns of activity in sensor-array data streams on a 2-D grid.

We gratefully acknowledge the support of grant NICHD 5R01HD022614 for the development of these routines.

## Documentation

User guide, installation instructions, workflow and usage examples are available [here](https://kutaslab.github.io/fitgrid).

## Demo

Click this button to launch a demo notebook:

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/kutaslab/fitgrid/main?filepath=notebooks/Demo.ipynb)
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
reported by contacting the project team at turbach@ucsd.edu. All
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
Gold standard test data files have been moved to
`fitgrid/fitgrid/data`.

This directory is retained as landing site for temp fake data files
generated during testing.

---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

Thank you for helping to improve fitgrid.

## Before you submit an issue please check your fitgrid installation.

Problems are unavoidable when fitgrid is installed with incompatible or missing
Python or R packages.

Please try to replicate the issue after first using mamba or conda to install the latest fitgrid [stable version](https://kutaslab.github.io/fitgrid/installation.html#fitgrid-stable-release) into a newly created conda virtual environment.

If the problem persists, please try to replicate the issue after installing the latest [pre-release version](https://kutaslab.github.io/fitgrid/installation.html#fitgrid-development-version) into a newly created conda virtual environment.

**Warning: `pip install fitgrid` is officially not supported.**


## Please provide the following information

### 1. Description
A clear and concise description of what the problem is, specifically:
- What you expected to happen and what actually happened.
- Anything you tried to solve the issue.

### 2. Minimal reproducible example
These are the shortest steps that reconstruct the problem.
- The **exact** character-for-character command(s) you ran to to install fitgrid, copy-paste is best.
- A Python code snippet or shell commands, the shorter the better, that runs and exhibits the issue.

### 3. Conda environment
Activate the conda environment that has fitgrid installed, run the following command in a terminal window, and upload the `fitgrid_issue.txt` file as an attachment with your issue.
```
conda list --explicit > fitgrid_issue.txt
```

### 4. System information
Please provide the specifics about your computer hardware architecture and
operating system version. For example:

- Linux, in a terminal window

```
	$ uname -mprsv
	Linux 3.10.0-957.21.3.el7.x86_64 #1 SMP Tue Jun 18 16:35:19 UTC 2019 x86_64 x86_64
	
	$ cat /etc/*-release
	CentOS Linux release 7.3.1611 (Core) 
	NAME="CentOS Linux"
	VERSION="7 (Core)"
	ID="centos"
	ID_LIKE="rhel fedora"
	VERSION_ID="7"
	PRETTY_NAME="CentOS Linux 7 (Core)"
	ANSI_COLOR="0;31"
	CPE_NAME="cpe:/o:centos:centos:7"
	HOME_URL="https://www.centos.org/"
	BUG_REPORT_URL="https://bugs.centos.org/"
	
	CENTOS_MANTISBT_PROJECT="CentOS-7"
	CENTOS_MANTISBT_PROJECT_VERSION="7"
	REDHAT_SUPPORT_PRODUCT="centos"
	REDHAT_SUPPORT_PRODUCT_VERSION="7"
	
	CentOS Linux release 7.3.1611 (Core) 
	CentOS Linux release 7.3.1611 (Core) 
```


- Mac OSX, in a terminal window

```
	$ uname -mprsv
	Darwin 19.6.0 Darwin Kernel Version 19.6.0: Thu Jun 18 20:49:00 PDT 2020; root:xnu-6153.141.1~1/RELEASE_X86_64 x86_64 i386

	$ sw_vers
	ProductName:	Mac OS X
	ProductVersion:	10.15.6
	BuildVersion:	19G2021
```

  
- [TODO: not officially supported] Windows, in a Windows Command Window

```
	C:\Users\some_user> systeminfo

	Host Name:                 DESKTOP-G57OVSM
	OS Name:                   Microsoft Windows 10 Home
	OS Version:                10.0.18362 N/A Build 18362
	OS Manufacturer:           Microsoft Corporation
	OS Configuration:          Standalone Workstation
	OS Build Type:             Multiprocessor Free
	Registered Owner:          some_user
	Registered Organization:
	Product ID:                00326-00840-79774-AAOEM
	Original Install Date:     7/29/2019, 6:13:59 AM
	System Boot Time:          9/9/2020, 9:07:46 AM
	System Manufacturer:       System manufacturer
	System Model:              System Product Name
	System Type:               x64-based PC
	Processor(s):              1 Processor(s) Installed.
    [01]: Intel64 Family 6 Model 158 Stepping 12 GenuineIntel ~3600 Mhz
	BIOS Version:              American Megatrends Inc. 0606, 8/31/2018
	Windows Directory:         C:\WINDOWS
	System Directory:          C:\WINDOWS\system32
	Boot Device:               \Device\HarddiskVolume2
	System Locale:             en-us;English (United States)
	Input Locale:              en-us;English (United States)
	Time Zone:                 (UTC-08:00) Pacific Time (US & Canada)
	Total Physical Memory:     16,305 MB
	Available Physical Memory: 13,598 MB
	Virtual Memory: Max Size:  18,737 MB
	Virtual Memory: Available: 14,542 MB
	Virtual Memory: In Use:    4,195 MB
```
.. _how_to_contribute:

############
Contributing
############


Field reports and ideas large and small for how to improve ``fitgrid``
are welcome. Please post what you have in mind on GitHub in
``fitgrid`` `Issues <https://github.com/kutaslab/fitgrid/issues>`_ in
accord with the `Code of Conduct
<https://github.com/kutaslab/fitgrid/blob/main/CODE_OF_CONDUCT.md>`_
to start a discussion of next steps and plan the approach. If you
think you have encountered a bug, please follow the bug report
guidlines when


========
Overview
========

The ``fitgrid`` code is written in Python, and requires many
open-source scientific computing packages including :std:doc:`numpy
<numpy:index>`, :std:doc:`pandas <pandas:index>`, :std:doc:`matplotlib
<matplotlib:index>`, and :std:doc:`statsmodels
<statsmodels:index>`. The linear mixed-effects modeling further
requires the :std:doc:`pymer4 <pymer4:index>` and :std:doc:`rpy2
<rpy2:index>` Python packages as well as the `R
<https://www.r-project.org/other-docs.html>`_ language, and many R
packages including `lme4
<https://cran.r-project.org/web/packages/lme4/index.html>`_, and
`lmerTest
<https://cran.r-project.org/web/packages/lmerTest/index.html>`_. Version
control is managed with `git <https://git-scm.com/doc>`_ and the
primary source repository is hosted on GitHub in
https://github.com/kutaslab.

The stable and pre-release :std:doc:`fitgrid <fitgrid:index>` packages
are available on https://anaconda.org/kutaslab/fitgrid for easy
installation into conda virtual environments with :std:doc:`conda
<conda:index>` or :std:doc:`mamba <mamba:index>`. The Python package
installer :std:doc:`pip <pip:index>` is not supported because of the R
dependencies.

The documentation is generated with :std:doc:`sphinx
<sphinx:contents>` and :std:doc:`sphinx-gallery
<sphinx-gallery:index>`. 


Versions
  ``fitgrid`` is semantically versioned following a simplified 
  Python `PEP 440 <https://www.python.org/dev/peps/pep-0440>`_ scheme.

  * **Syntax.** Legal ``fitgrid`` version strings have three required
    numeric segments, ``M.N.P`` for major, minor, patch versions and
    an optional ``.devX`` suffix with numeric ``X``. As
    usual, the numeric values increment monotonically except that subordinate
    segments reset to 0 when a superordinate segment
    increments. Version strings match the regular expression,
    ``\d+\.\d+\.\d+(\.dev\d+){0,1}``.

  * **Semantics.** ``M.N.P`` designates a stable version and ``M.N.P.devX``
    designates a development version.
    
  .. note::
     In this scheme, the stable package version ``M.N.P`` sorts higher
     than the development version ``M.N.P.devX``, for example,
     ``0.5.1`` > ``0.5.1.dev1`` > ``0.5.1.dev0`` > ``0.5.0``

GitHub reserved branches and version strings: https://github.com/kutsalab/fitgrid
  * The `main` and `dev` branches are reserved for deploying package
    releases, only maintainers make pull requests to these branches.
  * The branch `main` = the latest stable release, version ``M.N.P`` with tag ``vM.N.P``
  * The branch `dev` = the latest development version ``M.N.P.devX``.
  * Names and version strings for working branches other than `main`
    and `dev` are not strictly defined in this scheme. However, the
    natural version sequence increments the ``devX`` segment to
    ``M.N.P.devX+1`` following the pre-release package upload of
    ``M.N.P.devX`` and increments the patch and resets the ``devX``
    segment to ``M.N.P+1.dev0`` following the stable release of
    ``vM.N.P``.

Stable releases: ``vM.N.P``
  Stable versions are tagged ``vM.N.P`` and released manually on GitHub
  (``fitgrid`` `Releases
  <https://github.com/kutaslab/fitgrid/releases>`_). The stable
  version ``vM.N.P`` source code is frozen and subsequent
  modifications require the version to increment at least the patch
  segment. Development versions ``M.N.P.devX`` are not released on
  GitHub but they are deployed as conda packages (see next).

Conda packages and channels: stable ``M.N.P`` and pre-release ``M.N.P.devX``
  Conda packages are deployed for stable releases (``vM.N.P`` on branch
  `main`) and for development versions (``M.N.P.devX`` on branch
  `dev`). The stable release deploys to conda channel `kutaslab/label/main
  <https://anaconda.org/kutaslab/fitgrid/files>`_ and is the
  default for ``fitgrid`` conda installation. The development
  package deploys to channel `kutaslab/label/pre-release
  <https://anaconda.org/kutaslab/fitgrid/files>`_
  so the latest features and bug-fixes can be installed in conda
  environments with conda package dependency resolution.

Sphinx and sphinx-gallery documentation
  Documentation for the latest stable conda package ``vM.N.P`` is
  deployed to `gh-pages
  <https://github.com/kutaslab/fitgrid/tree/gh-pages>`_ and available
  online at https://kutaslab.github.io/fitgrid. Documentation for the
  latest development version ``M.N.P.devX`` is deployed to
  `gh-pages-dev
  <https://github.com/kutaslab/fitgrid-dev/tree/gh-pages-dev>`_ and
  available online at https://kutaslab.github.io/fitgrid-dev-docs.
  

Continuous Integration and Deployment (CID)
  The ``fitgrid`` CID is implemented in a single-pass GitHub Action
  workflow, `figrid-cid.yml
  <https://github.com/kutaslab/fitgrid/blob/main/.github/workflows/fitgrid-cid.yml>`_.
  The continuous integration workflow is triggered by push, pull
  request and manual release events on GitHub. The deploy phase
  selectively uploads the conda packages and documentation for
  development version pre-releases and stable releases. This scheme
  allows conda or mamba installation of both stable and development
  versions and automatically synchronizes the stable release version
  string and source code across the GitHub repository at
  `github.com/kutaslab/fitgrid
  <https://github.com/kutaslab/fitgrid>`_, the conda packages at
  `anaconda.org/kutaslab/fitgrid <https://anaconda.org>`_ , the online
  `sphinx documentation <https:kutaslab.github.io/fitgrid>`_, and the
  Zenodo source code archive at `DOI 10.5281/zenodo.3581496
  <https://doi.org/10.5281/zenodo.3581496>`_.


  .. _cid-figure:

  .. figure:: _static/fitgrid_cid_scheme.png

     Continuous Integration and Deployment Scheme

	    
  **Continuous Integration.** The conda package is built from the source
  on the triggering branch and installed into a newly created conda
  test environment.  The pytests in `fitgrid/tests/test_*.py` are run
  and the Sphinx html documentation is generated, including the
  sphinx-gallery `*.py` examples, in the test environment with the
  just-built package as installed.

  **Deployment**. If the CI passes, workflows triggered on branch `dev`
  with version string of the form ``M.N.P.devX`` or triggered by a
  GitHub manual releases tagged ``vM.N.P`` on branch main auto-upload
  the just-built conda package and Sphinx documentation to the
  appropriate destination repositories.

  * Pre-release: ``M.N.P.devX``

    * Conda packages: `--channel kutaslab/label/pre-release <https://anaconda.org/kutaslab/fitgrid/files>`_
    * Sphinx documentation: `kutaslab.github.io/fitgrid-dev-docs <https://kutaslab.github.io/fitgrid-dev-docs>`_
      
  * Stable release: ``vM.N.P``

    * Conda packages: `--channel kutaslab <https://anaconda.org/kutaslab/fitgrid/files>`_
    * Sphinx documentation: `kutaslab.github.io/fitgrid <https://kutaslab.github.io/fitgrid>`_
    * Zenodo archive DOI: `10.5281/zenodo.3581496 <https://doi.org/10.5281/zenodo.3581496>`_


Developing new features, bug fixes, and docs
  Updates to ``fitgrid`` source and docs are committed to working
  branches typically derived from the `kutaslab/fitgrid/dev` branch and not
  directly to the `main` or `dev` branches which are reserved for
  deploying conda packages and documentation. As development on the
  working branches progesses (magenta in the :ref:`cid-figure`),
  maintainers periodically pull the changes to the `dev` branch in
  order to deploy a pre-release package for installation into conda
  environments. When development is ready for a stable release,
  maintainers pull `dev` to the `main` branch and manually issue a
  stable release on GitHub tagged ``vM.N.P``. The tagged release
  uploads the ``M.N.P`` conda packages and sphinx documentation and
  archives the ``M.N.P`` source on Zenodo.


====================
Development workflow
====================

It is generally advisable to develop, test, and document new work
on a local computer in an active conda environment populated with the
latest compatible ``fitgrid`` dependencies along with :std:doc:`pytest
<pytest:index>`, the :std:doc:`black <black:index>` code formatter,
and sphinx documentation generation packages because that's what the
continuous integration workflow does.

The following illustrates the steps for a hypothetical working branch
called `new-feature` in the `github.com/kutaslab/fitgrid
<https://github.com/kutaslab/fitgrid>`_ GitHub repo. It assumes the
``git``, ``conda``, and ``mamba`` executables are already installed on
the local computer and the commands are executed in a bash(-like)
shell.


---------
Git setup
---------

#. Sign in to GitHub and create a fork of `github.com/kutaslab/fitgrid
   <https://github.com/kutaslab/fitgrid>`_ in your GitHub account.

#. On the local computer where you plan to work, ``git clone`` the
   fork.
   
   .. code-block:: bash

      $ git clone https://github.com/<your_github_username>/fitgrid

   By default, the local repo created this way will include the `main`
   branch only. Alternatively, the repo can be cloned with a specific
   working branch such as `new-feature` like so:

   .. code-block:: bash

      $ git clone https://github.com/<your_github_username>/fitgrid \
        --single-branch --branch new-feature


---------------------------------------
Development environment setup
---------------------------------------

#. Create a new named conda development environment for working on the
   feature, fix, or docs by installing the latest ``fitgrid``
   pre-release conda package, document generation, and development
   tools:

   .. code-block:: bash

      $ mamba create --name fg-new-feature \
           -c conda-forge -c ejolly -c kutaslab/label/pre-release \
           fitgrid
      $ mamba install --name fg-new-feature \
           black pytest sphinx sphinx-gallery sphinx_rtd_theme


#. Navigate to the top-level directory of your local fitgrid git
   repository, activate the new development environment, and install
   ``fitgrid`` from the local source in editable (a.k.a "develop")
   mode:

   .. code-block:: bash

      $ cd ~/path/to/fitgrid
      $ conda activate fg-new-feature
      (fg-new-feature) $ pip install --no-deps -e .

Why? Because installing the pre-release ``fitgrid`` conda package
automatically populates the just-created environment with the latest
compatible versions of the hundreds of Python, R, and matrix math
dependencies that the latest version of ``fitgrid`` needs to run. Then
``pip`` replaces the just-downloaded-and-installed ``fitgrid`` conda
package located in your
`~/path/to/conda/envs/fg-new-feature/path/to/site-packages/fitgrid`
with a link to your `~/path/to/fitgrid` local git repo. This way, the
files you modify are loaded when ``fitgrid`` modules are imported by
the pytests and sphinx document generators and your changes are
version-controlled by git.

.. note::

   Experience indicates this is the **only** time ``pip install``
   should be used while developing ``fitgrid`` on pain of corrupting
   the conda environment. If you want to add other packages to the
   development environment use ``mamba install`` or ``conda install``.


.. _dev_doc_test:

-----------------------
Develop, test, document
-----------------------

#. Activate the `fg-new-feature` development environment.

   .. code-block:: bash

      $ conda activate fg-new-feature
      (fg-new-feature) $

#. Checkout the git working branch. If it doesn't exist locally,
   ``git`` should automagically set it to track the remote working
   branch in your GitHub fork, make sure it does.

   .. code-block:: bash

      $ git checkout new-feature

#. Ensure the commit history of the `new-feature` branch in your
   GitHub fork and local repo are both up to date with the branch in
   the upstream GitHub repo `github.com/kutaslab/fitgrid
   <https://github.com/kutaslab/fitgrid>`_ where you will make the
   pull request (PR), i.e., `new-feature` in this example.  This helps
   reduce risk of merge conflicts later when changes are pulled back
   into the upstream repository.

#. Make the changes to the source code .py or docs .rst.

#. Document the .py source files with `numpy-style docstrings
   <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

#. For new functionality add or update pytests in
   `fitgrid/tests/test_*.py` to cover the changes.

#. If it is useful, add or update a working `.py` example in the
   relevant `fitgrid/docs/gallery` subdirectories for display in the
   :ref:`gallery`.

#. Run a code checker such as `flake8` or `pylint` on the .py files.

#. Navigate to the top level of the ``fitgrid`` repository, run the
   code formatter and pytests the same way the GitHub Action CI does:

   .. code-block:: bash

      (fg-new-feature) $ black -S --line-length 79 .
      (fg-new-feature) $ pytest

#. When pytests pass, navigate to the top level of the ``fitgrid``
   repository and build the sphinx documentation the same way the
   GitHub Action CI does:

   .. code-block:: bash

      (fg-new-feature) $ make -C docs clean; make -C docs html

   Monitor the docs building for errors and warnings, then open the
   local file `~/path/to/fitgrid/docs/build/html/index.html` in your
   web browser and visually verify that the .rst docs and
   sphinx-gallery ``*.py`` Python examples in the subdirectories of
   `~/path/to/fitgrid/docs/source/gallery` produce the expected
   results and figures.

-------------------
Pull requests (PRs)
-------------------

#. When pytests pass and documentation builds locally, commit the
   changes on branch `new-feature` and push the working branch to your
   forked ``fitgrid`` repository on GitHub.

#. Sign in to GitHub, navigate to your fork's Action tab and verify
   that the push on branch `new-feature` triggered an Action
   workflow that runs without error.

#. If the workflow fails, inspect the Action log, diagnose the
   problem, go back to :ref:`dev_doc_test`, fix the problem in the
   local repo, commit the changes, and push them to the forked
   repository.

#. When the CI workflow for branch `new-feature` passes in the forked
   repository, make a pull request back to the upstream working branch.


====   
Tips
====

* Run ``conda list fitgrid`` to confirm it is installed in editable
  mode in the active development environment. It should look something
  like this:
  
  .. code-block:: bash

     (fg-new-feature) userid@machine$ conda list fitgrid
     # packages in environment at /home/userid/miniconda3/envs/fg-new-feature:
     #
     # Name                    Version                   Build  Channel
     fitgrid                   0.5.1.dev5                dev_0    <develop>


  Check that Version matches the version string in
  `fitgrid/__init__.py` in your local source git repo and the conda
  Channel is `<develop>`.

* If you plan to use :std:doc:`Jupyter <jupyter:index>` or
  :std:doc:`JupyterLab <jupyterlab:index>` to develop code or
  documentation examples things may go more smoothly if you ``mamba
  install`` or ``conda install`` the package into the development
  enviroment where you are working on ``fitgrid``.

* If working in a Jupyter notebook, you can use
  :py:func:`importlib.reload` to load modified source code between
  kernel restarts.

* You can rebuild the .rst documentation quickly without running the sphinx-gallery 
  Python examples by running this command in the top-level repository directory: ::

    make -C docs html-noexec 
.. _about_fitgrid:

#################
About ``fitgrid``
#################


============================
EEG and signal-averaged ERPs
============================

In the late 1920's Berger demonstrated that some of the human brain
activity related to external stimulation and internal mental events
could be measured at the surface of the scalp as tiny time-varying
electrical potential waveforms on the order of tens of microvolts
peak-to-peak, the human electroencephalogram (EEG) [Berger1930]_. In
the early 1950s Dawson presented a demonstration that even tinier
brain responses to external stimulation that were too small to be seen
by the naked eye in the EEG could, however, be observed by repeating
the stimulation multiple times, aligning fixed-length segments
("epochs") of the EEG recordings to the onset of the stimulus and
summing the recordings together at each time-point ([Dawson1951]_,
[Dawson1954]_). The idea of aggregating several noisy measurements so
the pluses and minuses of random variation cancel to yield a better
estimate the "true" value of the sum or mean was already
well-known. The conceptual breakthrough for EEG data analysis that
Dawson credited to Hunt was to sweep the familiar noise-reduction
trick along the time-aligned EEG epochs to supress larger variation in
the EEG (noise) and reveal the time course of the much smaller brain
response (signal) to the stimulation, a waveform on the order of
microvolts peak-to-peak. Laboratory experiments in subsequent decades
found that the Hunt-Dawson aggregation procedure could reveal a
variety systematic brain responses before, during, and after sensory
stimulation, motor responses, and internal mental events. With the
advance of computer hardware and software, oscilloscopic sweep
averagers were replaced by analog-to-digital conversion and
sum-and-divide averaging in software on general purpose
computers. Since the 1970s, this discrete time series average
event-related brain potential (ERP) has been a cornerstone of
experimental EEG research on human sensation, perception, and
cognition. For a compendium see the Oxford Handbook of Event-related
Potentials, [LucKap2011]_.


===============
regression ERPs
===============

In 2015 Smith and Kutas published a paper noting
that the average of a set of values, :math:`y`, is identical to
the estimated constant, :math:`\hat{\beta}_{0}` for the linear model

.. math::

  y = \beta_{0} + e

fit by minimizing squared error (ordinary least squares)
[SmiKut2015]_. They pointed out that this makes the average ERP a
special case of sweeping a linear regression model along the EEG at
time, *t*, and generalized this to more complex multiple regression
models,

.. math::

   y(t) = \beta_{0}(t) + \beta_{1}X_{1}(t) \ldots \beta_{n}X_{i}(t) + e(t)

Sweeping any such model along the EEG time point by time point
likewise produces time series of estimates for the intercept and
regressor coefficients, the :math:`\hat{\beta}_{i}` they dubbed the
"regression ERP" (rERP) waveforms. More generally it produces a time
series for each of the the model estimates and derived quantities,
such as coefficient standard errors, residuals, :math:`R^2`, likelihood,
Akiake's information criterion, and so forth.

This insight extends sum-and-divide Hunt-Dawson aggregation and embeds
event-related EEG data analysis in a general framework for
discovering, evaluating, and comparing a wide range of models to
account for systematic variation in the time course of EEG responses
using well-established methods of applied regression. With
this shift, however, comes a new problem.

==========================================================
Modeling: fit, diagnose, compare, evaluate, revise, repeat
==========================================================

These days specifying and fitting a linear regression model is a
matter of organizing the data into a table of rows (observations) and
columns (variables), typing a model specification formula like
:math:`\mathsf{1 + a + b + a:b}` and pressing Enter. While **fitting** a model is
relatively easy and mechanical, **modeling**, by contrast, is a laborious
process that iterates cycles of data quality control, fitting,
data diagnosis, model evaluation, comparison, and selection with numerous
decision points that require thought and judgment along the way.

Modeling EEG data as regression ERPs at each time point and data
channel multiplies the iterative cycles in a combinatorial explosion
of time points :math:`\times` channels :math:`\times` models
:math:`\times` comparisons. For example, at a digital sampling rate of
250 samples per second, there are 750 time points in 3 seconds of EEG
data. For 32 EEG channels, this makes 750 timepoints x 32 channels =
24,000 data sets. To fit three candidate models requires 72,000
separate model fits where the size of the data set might range
anywhere from a few dozens of observations for a single subject to
tens of thousands of observations for a large scale experiment.

Nothing can prevent the combinatorial explosion; ``fitgrid``
is designed to contain it.


=========================================
``fitgrid``: Modeling :math:`\times` 10e4
=========================================

The ``fitgrid`` package allows researchers generally familiar with
regression modeling and model specification formulas in Python
(:std:doc:`statsmodels.formula.api <statsmodels:index>` via
:std:doc:`patsy <patsy:index>`) or R (``lm``, ``lme4``, ``lmerTest``)
to use these tools to readily and reproducibly fit ordinary least
squares and linear mixed-effects regression models of multichannel
event-related time series recordings, at scale, with a few lines of
scripted Python.

With one function call, ``fitgrid`` sweeps a model formula across the
data observations at each time and channel (in parallel on multiple
CPU cores if supported by hardware) and collects the resulting fit
objects returned by ``statsmodels.ols`` or ``lme4::lmer`` via
``pymer4`` in a single :py:class:`FitGrid[times, channels]
<fitgrid.fitgrid.FitGrid>` Python object.

The ``fitgrid`` can be sliced by time and channel like a dataframe, and
the results for a fit attribute are queried for the entire grid with
the same syntax as single fit: ``results = FitGrid.<attr>``. The
results include the time-series of coefficient estimates comprising
the regression ERPs, including, but not restricted to, the special
case average ERP.  Equally important for modeling, the results also include
everything else in the bundle of information comprising the fit object
such as coefficient standard errors, model log likelihood, Akiake's
information criterion, model and error mean squares, and so forth. The
results are returned as tidy Time x Channel dataframes for handy
visualization and analysis in Python and data interchange across
scientific computing platforms as illustrated in
:ref:`workflow` and the :ref:`gallery`. For examles of ``fitgrid``
at work see [TroUrbKut2020]_ and [UrbachEtAl2020]_.

================================
``fitgrid`` design: How it works
================================

Ordinary least squares models are fit in Python using the
:std:doc:`statsmodels <statsmodels:index>` [SeaSkiPer2010]_
statistical modeling package via the :std:doc:`patsy <patsy:index>` formula
language interface [Smith2020]_.  Linear mixed effects models are
shipped out of Python and into R via Eshin Jolly's
:py:class:`pymer4.models.Lmer` interface [Jolly2018]_ and fit with
`lme4::lmer
<https://cran.r-project.org/web/packages/lme4/index.html>`_ (see
[BatesEtAl2015]_).

For illustration with ``patsy`` and ``statsmodels``, suppose you have a
:py:class:`pandas.DataFrame` ``data`` with independent variables ``x``
and ``a``, where ``x`` is continuous and ``a`` is categorical. Suppose
also ``channel`` is your continuous dependent variable.  Here's how
you would run an ordinary least squares linear regression of
``channel`` on ``x + a`` using ``statsmodels``::

    from statsmodels.formula.api import ols

    fit = ols('channel ~ x + a', data).fit()

Now this ``fit`` object contains all the fit and diagnostic information,
mirroring what is provided by ``lm`` in R. This information can be retrieved by
accessing various attributes of ``fit``. For example, the betas::

    betas = fit.params

or the t-values::

    tvalues = fit.tvalues

or :math:`Pr(>|t|)`::

    pvalues = fit.pvalues

Compare to R, where this is usually done by calling functions like ``summary``
or ``coef``.

Now the issue with using that interface for single trial rERP analyses
is of course the dimensionality: instead of fitting a single model, we
need to fit :math:`m \times n` models, where :math:`m` is the number
of discrete time points and :math:`n` is the number of channels.

This can be handled using ``for`` loops of the form::

    for channel in channels:
        for timepoint in timepoints:
            # run regression 'channel ~ x + a', save fit object somewhere

And to access some particular kind of fit information, the exact same two
nested ``for`` loops are required::

    for channel in channels:
        for timepoint in timepoints:
            # extract diagnostic or fit measure, save it somewhere


``fitgrid`` abstracts this complexity away and handles the iteration and
storage of the data behind the scenes. The first loop above is now replaced
with::

    lm_grid = fitgrid.lm(epochs, RHS='x + a')

and the second loop with::

    betas = lm_grid.params

or::

    tvalues = lm_grid.tvalues

or::

    pvalues = lm_grid.pvalues

The crux of the approach conceived and implemented by Andrey Portnoy
is that ``lm_grid``, a :py:class:`FitGrid[times, channels]
<fitgrid.fitgrid.FitGrid>` object, can be queried for the exact same
attributes as a regular ``statsmodels`` fit object as above.

The result is most often a :py:class:`pandas.DataFrame`, sometimes
another :py:class:`FitGrid[times, channels]
<fitgrid.fitgrid.FitGrid>`. In other words, if you are running linear
regression, the attributes of a fit object documented in the
``statsmodels`` :py:class:`linear_model.RegressionResults
<statsmodels.regression.linear_model.RegressionResults>` API, can be
used to query a :py:class:`FitGrid[times, channels]
<fitgrid.fitgrid.FitGrid>`.

``statsmodels``::

    fit.rsquared

``fitgrid``::

    lm_grid.rsquared

Some of the attributes are methods. For example, influence diagnostics
in ``statsmodels`` are stored in a separate object that is created by
calling the ``get_influence`` method. So Cook's distance measures can
be retrieved as follows::

    influence = fit.get_influence()
    cooks_d = influence.cooks_distance

The exact same approach works in ``fitgrid``::

    influence = lm_grid.get_influence()
    cooks_d = influence.cooks_distance


============================
``fitgrid`` in other domains
============================

Although the origins of ``fitgrid`` are in EEG data analysis,
``fitgrid`` can also be used with sensor array time-series data from
other domains where event-related signal averaging and and regression
modeling is appropriate. The :ref:`gallery` includes hourly NOAA tide
and atmospheric data to illustrate event-related time-domain
aggregation to detect lunar atmospheric tides, an approach first
attempted by Laplace in the early 19th century [LinCha1969]_.

.. _user_guide:


.. module:: fitgrid
   :noindex:


##########
User Guide
##########

TL;DR These are notes and highlights. For usage see
the :ref:`workflow`, :ref:`gallery`, and :ref:`api`


======================================================
``fitgrid`` :py:class:`Epochs <fitgrid.epochs.Epochs>`
======================================================

Fitting linear regression models in Python and R with formulas like
``y ~ a + b + a:b`` in Python (:std:doc:`patsy <patsy:index>`,
:std:doc:`statsmodels.formula.api <statsmodels:api>`) or R (``lm``,
``lme4::lmer``) assumes the data are represented as a 2D array with the
variables in named columns, ``y``, ``a``, ``b`` and values
("observations", "scores") in rows.

``fitgrid`` follows this format with the additional assumption that
the data are vertically stacked fixed-length time-series "epochs", so
the user must specify two additional row indexes that together
uniquely identify the epoch and timestamp of each data row.


.. _epochs_data_format:

------------------
Epochs data format
------------------

Data for ``fitgrid`` modeling should be prepared as a
single :py:class:`pandas.DataFrame` with a
:py:class:`pandas.MultiIndex` as follows:

MultiIndex names and data types
  Each epoch must have a unique integer identifier in the ``epoch_id``
  index. The index values need not be ordered and gaps are allowed but
  duplicate epoch indices are not. The ``time`` index integer
  timestamps must be the same for all epochs. Gaps and irregular
  intervals between timestamps are allowed, duplicate timestamps within an
  epoch are not.
  
  * ``epoch_id`` (default): integer
  * ``time`` (default): integer

  The MultiIndex names default to ``epoch_id`` and ``time`` but any
  index column may be designated as the epoch or time index, provided
  the conditions are met. Event-relative epochs are often time-stamped
  so the event is at time=0 but this is not a ``fitgrid`` requirement
  and the temporal resolution of the timeseries is not specified.  For
  data epochs with time index values -10, 0, 10, 20, the
  :py:class:`FitGrid[-10:20, channels] <fitgrid.fitgrid.FitGrid>`
  object will be the same whether the timestamps denote milliseconds
  or months.

Data columns
  * channel data columns: numeric
  * predictor variable columns: numeric, string, boolean

Example:
  A canonical source of data epochs for ``fitgrid`` are multichannel
  "strip charts" as in EEG and MEG recordings. In this case, the
  epochs are regularly sampled fixed-length segments extracted from
  the strip chart and time-stamped relative to an experimental event.

.. image:: _static/eeg_epochs.png


--------------------------------------
:ref:`Data Ingestion <data_ingestion>`
--------------------------------------

Rows and columns epochs data can be loaded into a `fitgrid`
:py:class:`Epochs <fitgrid:epochs.Epochs>` object from a
:py:class:`pandas.DataFrame` in memory or read from files in `feather
<https://arrow.apache.org/docs/python/feather.html>`_ or `HDF5
<https://portal.hdfgroup.org/display/HDF5/HDF5>`_ format. For details
on using these data formats see :py:func:`pandas.read_feather` and and
:py:func:`pandas.read_hdf`.


----------------------------------------
:ref:`Data Simulation <data_simulation>`
----------------------------------------

``fitgrid`` has a built-in method :py:meth:`fitgrid.fake_data.generate` that
returns a :py:class:`Epochs <fitgrid:epochs.Epochs>` data objects or
:py:class:`pandas.DataFrame` for testing.


---------------
EEG Sample Data
---------------

Optional sample files containing previously prepared EEG data epochs
are availble for download from the Zenodo archive at
`eeg-workshops/mkpy_data_examples/data
<https://zenodo.org/record/3968485>`_.

The files can be installed into the fitgrid package with
:py:meth:`fitgrid.sample_data.get_file` and accessed thereafter as
``fitgrid.DATA_DIR(<filename>)`` or downloaded manually to another
location.

The sample data with `msNNN.epochs.*` filenames contain segmented EEG
epochs and the `msNNNN` infix gives the length of the epoch in
millesconds. Files with the shortest epochs (100 ms) are suitable for
testing, those with longer epochs (1500, 3000 ms) are more
representative of actual experimental EEG data. The feather format
versions `*.epochs.feather` are recommended for use with pandas
:py:func:`pandas.read_feather`.

Additional information about the experimental designs for these sample
files is available online at https://eeg-workshops.github.io/mkpy_data_examples.


===============
Fitting a model
===============


The following methods populate the :py:class:`FitGrid[times,channel]
<fitgrid.fitgrid.FitGrid>` object with `statsmodels` results for OLS
model fits and `lme4::lmer` for linear mixed-effects fits.

* Ordinary least squares: :py:meth:`fitgrid.lm <fitgrid.models.lm>`

  .. code-block:: python

     lm_grid = fitgrid.lm(
         epochs_fg,
         RHS='1 + categorical + continuous'
     )



* Linear mixed-effects: :py:meth:`fitgrid.lmer <fitgrid.models.lmer>`

  .. code-block:: python

     lmer_grid = fitgrid.lmer(
         epochs_fg,
         RHS='1 + continuous + (continuous | categorical)'
     )



* User-defined (experimental): :py:meth:`fitgrid.run_model <fitgrid.models.run_model>`


==================================================================
The :py:class:`FitGrid[times, channels] <fitgrid.fitgrid.FitGrid>`
==================================================================


--------------------------
Slice by `time`, `channel`
--------------------------


Slice the :py:class:`FitGrid[times, channels] <fitgrid.fitgrid.FitGrid>`
with :py:class:`pandas.DataFrame` range ``:`` and label slicers.  The
range includes the upper bound.

.. code-block:: python

   lm_grid[:, ["MiCe", "MiPa"]]
   lm_grid[-100:300, :]
   lm_grid[0, "MiPa"]


--------------
Access results
--------------


Query the :py:class:`FitGrid[times, channels] <fitgrid.fitgrid.FitGrid>`
results like a single fit object. Result grids are returned as a
`pandas.DataFrame` or another :py:class:`FitGrid[times, channels]
<fitgrid.fitgrid.FitGrid>` which can be queried the same way.

.. code-block:: python

   lmg_grid.params
   lmg_grid.llf


----------------
Slice and access
----------------

.. code-block:: python

   lm_grid[-100:300, ["MiCe", "MiPa"].params


------------------------------------
:py:meth:`fitgrid.LMFitGrid` methods
------------------------------------

The fitted OLS grid provides time-series plots of selected model
results: estimated coefficients :py:meth:`fitgrid.lm.plot_betas` and
adjusted :math:`R^2` :py:meth:`fitgrid.lm.plot_adj_rsquared` (see also
:py:meth:`fitgrid.utils` for additional model summary wrappers).


========================
Saving and loading grids
========================

Running models on large datasets can take a long time. `fitgrid` lets
you save your grid to disk so you can restore them later without
having to refit the models. However, saving and loading large grids
may still be slow and generate very large files.

Suppose you run `lmer` like so::

    grid = fitgrid.lmer(epochs, RHS='x + (x|a)')

Save the ``grid``::

    grid.save('lmer_results')

Later you can reload the ``grid``::

    grid = fitgrid.load_grid('lmer_results')


.. warning::

   Fitted grids are saved and loaded with Python `pickle` which is not
   guaranteed to be portable across different versions of Python.
   Unpickling unknown files **is not secure** (for details see the
   Python `docs
   <https://docs.python.org/3/library/pickle.html>`_). Only load grids
   you trust such as those you saved yourself. For reproducibility and
   portability fit the grid, collect the results you need, and export
   the dataframe to a standard data interchange format.



.. _guide_summaries:

===============================
Model comparisons and summaries
===============================

To reduce memory demands when comparing sets of models, `fitgrid`
provides a convenience wrapper, `fitgrid.utils.summarize`, that
iteratively fits a list of models and collects a lightweight summary
dataframe with key results for model interpretation and
comparison. Unlike :py:class:`FitGrid[times, channels]
<fitgrid.fitgrid.FitGrid>` objects, the summary dataframe format is the
same for `fitgrid.lm` and `fitgrid.lmer`. Some helper functions are
available for visualizing selected summary results.


.. _diagnostics:

================================
Model and data diagnostics (WIP)
================================

Model and data diagnostics in the `fitgrid` framework is work in
progress. For ordinary least squares fitting, there is some support
for the native `statsmodels` OLS diagnostic measures.  Diagostics that
can be computed analytically from a single model fit, e.g., via the
hat matrix diagonal, may be useable but many are not for realistically
large data sets. The per-observation diagnostic measures, e.g., the
influence of observations on estimated parameters, are the same size
as the original data multiplied by the number of model parameters
which may overload memory and measures that require on
leave-one-observation-out model refitting take intractably long for
large data sets. A minimal effort is made to guard the user from known
trouble but the general policy is `fitgrid` stays out of the way
so you can try what you want. If it works great, if it chokes, that's
the nature of the beast you are modeling.

Support for linear-mixed effects diagnostics in `fitgrid` is limited
to a variance inflation factor computation implemented in Python as a
proof-of-concept. `fitgrid` does not interface with mixed-effect model
diagnostics libraries in R and plans are to improve
support for mixed-effects modeling in Python rather than expand further
into the R ecosystem.



========================
`fitgrid` under the hood
========================


--------------------------------
How mixed effects models are run
--------------------------------

Mixed effects models do not have a complete implementation in Python, so we
interface with R from Python and use `lme4` in R. The results that you get when
fitting mixed effects models in `fitgrid` are the same as if you used `lme4`
directly, because we use `lme4` (indirectly).


-----------------------
Multicore model fitting
-----------------------

On a multicore machine, it may be possible to significantly speed
fitting by computing the models in parallel especially for
:py:meth:`fitgrid.lmer <fitgrid.models.lmer>`. For least squares
models, :py:meth:`fitgrid.lm <fitgrid.models.lm>` uses :py:mod:`statsmodels`
under the hood, which in turn employs :py:mod:`numpy` for calculations.
:py:mod:`numpy` itself depends on linear algebra libraries that might be
configured to use multiple threads by default. This means that on a 48
core machine, common linear algebra calculations might use 24 cores
automatically, without any explicit parallelization. So when you
explicitly parallelize your calculations using Python processes (say 4
of them), each process might start 24 threads. In this situation, 96
CPU bound threads are wrestling each other for time on the 48 core
CPU. This is called oversubscription and results in *slower*
computations.

To deal with this when running :py:meth:`fitgrid.lm
<fitgrid.models.lm>`, we try to instruct the linear algebra libraries
your :py:mod:`numpy` distribution depends on to only use a single
thread in every computation. This then lets you control the number of
CPU cores being used by setting the ``n_cores`` parameter in
:py:meth:`fitgrid.lm <fitgrid.models.lm>` and :py:meth:`fitgrid.lmer
<fitgrid.models.lmer>`.

If you are using your own 8-core laptop, you might want to use all but
one core, so set something like ``n_cores=7``. On a shared machine,
it's a good idea to run on half or 3/4 of the cores if no one else is
running heavy computations.

Note that fitgrid parallel processing counts the "logical" cores
available to the operating system and this may differ from the number
of physical cores, depending on the system hardware and setting, e.g.,
Intel CPUs with hyperthreading enabled. The Python package
`psutil <https://psutil.readthedocs.io/en/latest/>`_ and
``psutil.cpu_count(logical=True)`` and
``psutil.cpu_count(logical=False)`` may be useful for interrogating
the system about the available resources.
.. _bibliography:

============
Bibliography
============

.. _fitgrid_reports:

Research reports using `fitgrid`
--------------------------------

.. [TroUrbKut2020] Troyer, M., Urbach, T.P., &
      Kutas, M. (2020). Toward dissociating general reading experience
      and domain-specific knowledge sources during RSVP reading: An
      exploratory rERP analysis. Virtual poster presentation, Society
      for the Neurobiology of Language Annual Meeting,
      October 2020. https://doi.org/10.17605/OSF.IO/YNV9K

.. [UrbachEtAl2020] Urbach, T. P., DeLong, K. A., Chan, W-H., &
      Kutas, M. (2020). An exploratory data analysis of word form
      prediction during word-by-word reading. Proceedings of the
      National Academy of
      Sciences, 201922028. https://doi.org/10.1073/pnas.1922028117


.. _references:

References
----------

.. [BatesEtAl2015] Douglas Bates, Martin Maechler, Ben Bolker,
       Steve Walker (2015). Fitting Linear Mixed-Effects Models Using
       lme4. Journal of Statistical Software, 67(1),
       1-48. https://doi:10.18637/jss.v067.i01.

.. [BenHoc1995] Benjamini, Y., & Hochberg, Y. (1995). Controlling the
      false discovery rate: A practical and powerful approach to
      multiple testing. Journal of the Royal Statistical
      Society. Series B (Methodological), 57, 289-300.

.. [BenYek2001] Benjamini, Y., & Yekutieli, D. (2001). The control of
       the false discovery rate in multiple testing under
       dependency.The Annals of Statistics, 29, 1165-1188.

.. [Berger1930] Berger, H. (1930). Electroencephalogram of man
      (P. Gloor, Trans.). In P. Gloor (Ed.), Hans Berger on the
      Electroencephalogram of Man. The Fourteen Original Reports on
      the Human Electroencephalogram. 1969. Electroencephalography and
      Clinical Neurophysiology, supplement 28. New York: Elsevier.

.. [BurAnd2004] Burnham, K. P., & Anderson, D. R. (2004). Multimodel
       inference - understanding AIC and BIC in model
       selection. Sociological Methods & Research, 33(2),
       261-304. https://doi:10.1177/0049124104268644

.. [Dawson1951] Dawson, G. D. (1951). A summation technique for
       detecting small signals in a large irregular
       background. Journal of Physiology, 115(1), P2-P3.

.. [Dawson1954] Dawson, G. D. (1954). A summation technique for the
       detection of small evoked potentials. Electroencephalography
       and Clinical Neurophysiology, 6(1),
       65-84. https://doi:10.1016/0013-4694(54)90007-3

.. [Gautier2021] Laurent Gautier. rpy2 - R in
       Python. https://rpy2.github.io

.. [Jolly2018] Jolly, (2018). Pymer4: Connecting R and Python for
       Linear Mixed Modeling. Journal of Open Source Software, 3(31),
       862, https://doi.org/10.21105/joss.00862

.. [KuzBroChr2017] Kuznetsova A, Brockhoff PB, Christensen RHB
      (2017). lmerTest Package: Tests in Linear Mixed Effects
      Models. Journal of Statistical Software, 82(13), 1?26. doi:
      10.18637/jss.v082.i13.

.. [LinCha1969] Lindzen, R. S., & Chapman, S. (1969). Atmospheric
       tides. Space science reviews, 10(1), 3-188.

.. [LucKap2011] Luck, S. J., & Kappenman, E. S. (Eds.). (2011). The
       Oxford handbook of event-related potential components. Oxford
       University Press. https://doi:10.1093/oxfordhb/9780195374148.001.0001

.. [NOS-CO-OPS2] Computational Technniques for tidal datums handbook. N0AA Special Publication
       NOS CO-OPS2, National Oceanic and Atmospheric Administration,
       National Ocean Service, Center for Operational Oceanographic
       Products and Services. URL:
       https://tidesandcurrents.noaa.gov/publications/Computational_Techniques_for_Tidal_Datums_handbook.pdf.

.. [NieGroPel2012] Nieuwenhuis, R., Grotenhuis, M., &
       Pelzer, B. (2012).  influence.ME: tools for detecting
       influential data in mixed effects models.  R journal, 4(2),
       38-47.

.. [R2020] R Core Team (2020). R: A language and environment for
      statistical computing. R Foundation for Statistical Computing,
      Vienna, Austria.  URL: <https://www.R-project.org>.

.. [SeaSkiPer2010] Seabold, Skipper, and Josef Perktold. statsmodels:
       Econometric and statistical modeling with python.  Proceedings
       of the 9th Python in Science Conference. 2010. Documentation:
       <https://www.statsmodels.org>.

.. [SmiKut2015] Smith, N. J., & Kutas, M. (2015). Regression-based
      estimation of ERP waveforms: I. The rERP
      framework. Psychophysiology. doi:10.1111/psyp.12317. `[PubMed
      open access]
      <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5308234/>`_

.. [Smith2020] Smith, N. J., patsy - Describing statistical models in
       Python. Documentation: <https://patsy.readthedocs.io>.

       
.. _gallery:

Examples Gallery
================

.. toctree::
   :maxdepth: 2
   :glob:

   auto_gallery/workflow
   auto_gallery/1_epochs_data/index
   auto_gallery/2_model_fitting/index
   auto_gallery/3_model_evaluation/index

This is a list of tools available in ``fitgrid``.

.. module:: fitgrid
   :noindex:



.. _data_ingestion:

==============
Data Ingestion
==============

Functions that read epochs tables and create ``Epochs`` and load ``FitGrid``
objects.

.. autofunction:: epochs_from_dataframe
   :noindex:

.. autofunction:: epochs_from_hdf
   :noindex:

.. autofunction:: load_grid
   :noindex:



.. _data_simulation:

===============
Data Simulation
===============

``fitgrid`` has a built-in function that generates data and creates ``Epochs``:

.. autofunction:: generate
   :noindex:


==================
``Epochs`` methods
==================

Models and plotting.

.. autofunction:: fitgrid.epochs.Epochs.plot_averages
   :noindex:

=============
Model running
=============

.. autofunction:: fitgrid.lm
   :noindex:

.. autofunction:: fitgrid.lmer
   :noindex:

.. autofunction:: fitgrid.run_model
   :noindex:

===================
``FitGrid`` methods
===================


.. autofunction:: fitgrid.fitgrid.FitGrid.save
   :noindex:


=====================
``LMFitGrid`` methods
=====================

Plotting and statistics.

.. autofunction:: fitgrid.fitgrid.LMFitGrid.influential_epochs
   :noindex:

.. autofunction:: fitgrid.fitgrid.LMFitGrid.plot_betas
   :noindex:

.. autofunction:: fitgrid.fitgrid.LMFitGrid.plot_adj_rsquared
   :noindex:


.. _model-diagnostic-utilities:

=========
Utilities
=========

-------------------
model fit summaries
-------------------

.. autofunction:: fitgrid.utils.summary.summarize
   :noindex:

.. autofunction:: fitgrid.utils.summary.plot_betas
   :noindex:

.. autofunction:: fitgrid.utils.summary.plot_AICmin_deltas
   :noindex:


--------------
lm diagnostics
--------------

.. autofunction:: fitgrid.utils.lm.get_vifs
   :noindex:

.. autofunction:: fitgrid.utils.lm.list_diagnostics
   :noindex:

.. autofunction:: fitgrid.utils.lm.get_diagnostic
   :noindex:

.. autofunction:: fitgrid.utils.lm.filter_diagnostic
   :noindex:


----------------
lmer diagnostics
----------------

.. autofunction:: fitgrid.utils.lmer.get_lmer_dfbetas
   :noindex:

.. _installation:

############
Installation
############

**TL;DR** Use :std:doc:`mamba <mamba:index>` or :std:doc:`conda
<conda:index>` to install the ``fitgrid`` conda package along with your
other packages into a fresh conda environment on a fast multicore
x64_86 computer with gobs of RAM.


================================
About conda virtual environments
================================

``fitgrid`` is packaged on `anaconda.org/kutaslab/fitgrid
<https://anaconda.org/kutaslab/fitgrid>`_ for installation into conda
"virtual" environments using the `conda <https://conda.io>`_ or
:std:doc:`mamba <mamba:index>` package manager. A virtual environment
isolates the ``fitgrid`` installation to prevent clashes with what is
already installed elsewhere in your system and other virtual
environments. When the package manager installs ``fitgrid`` in a
virtual environment it also automatically installs compatible versions
of the hundreds of Python and R packages ``fitgrid`` requires to run
including :std:doc:`numpy <numpy:index>`, :std:doc:`pandas
<pandas:index>`, :std:doc:`matplotlib <matplotlib:index>`,
:std:doc:`statsmodels <statsmodels:index>`, and :std:doc:`pymer4
<pymer4:index>`, :std:doc:`rpy2 <rpy2:index>`, `R
<https://www.r-project.org/other-docs.html>`_, `lme4
<https://cran.r-project.org/web/packages/lme4/index.html>`_, and
`lmerTest
<https://cran.r-project.org/web/packages/lmerTest/index.html>`_ to
name a few. You can also install other conda packages in addition to
``fitgrid`` as needed for the task at hand.

The steps for creating conda environments and installing ``fitgrid``
are straightforward but it is prudent to have a general understanding
of conda virtual environments and at least these commands: ``conda
create ...``, ``conda install -c ...``, ``conda activate ...``, and
``conda deactivate``. See the `Conda Cheat Sheet
<https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html>`_
for a summary. For fine-tuning conda environments and working around
incompatible package versions refer to the core `conda tasks
<https://conda.io/projects/conda/en/latest/user-guide/tasks/index.html>`_
especially `managing conda channels and channel priority
<https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html>`_
and `installing packages
<https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-pkgs.html#installing-packages>`_.
The :std:doc:`mamba <mamba:index>` package installer is an alternative
to ``conda``. At present, the ``mamba create ...`` and ``mamba install
...`` commands tend to resolve the complex ``fitgrid`` package
dependencies substantially faster than ``conda``.

For working with ``fitgrid`` and mixed Python and R conda environments
generally, it is important to attend to the difference between the
``conda`` default channels `anaconda.org/main
<https://anaconda.org.main>`_ and `anaconda.org/r
<https://anaconda.org/r>`_ where packages are maintained by the
Anaconda, Inc. team and the not-always-compatible parallel universe of
the `conda-forge <https://conda-forge.org/>`_ channel where many of
the same-named conda packages are maintained by the open-source
community. Choosing suitable channels for installing conda packages
can be tricky because the specific versions of packages required for
compatibility and best performance depend on the computing hardware,
operating systems, compilers, and the version requirements of packages
that are already in the environment or will be. The conda-forge
maintainers recommend setting the .condarc `configuration file
<https://docs.conda.io/projects/conda/en/master/user-guide/configuration/use-condarc.html#using-the-condarc-conda-configuration-file>`_
to strict conda-forge channel priority.  However, revising the default
channel priority may not be appropriate for all users, in which case
the command-line options ``--channel conda-forge
--strict-channel-priority`` may be used. The examples below illustrate
command line options for a few common installation scenarios
encountered in practice.

Our current recommended best practice for working with conda
environments is to install the lightweight `miniconda3
<https://docs.conda.io/en/latest/miniconda.html>`_ and then avoid
polluting the "base" conda environment with data analysis and
application packages like ``fitgrid``.  Instead, create separate new
working environments, each populated with the packages needed for a
given project. The :std:doc:`mamba <mamba:index>` package is an
exception to this rule. If you elect to use ``mamba`` follow the
`installation instructions
<https://mamba.readthedocs.io/en/latest/installation.html>`_
carefully.


.. _conda_install_fitgrid:

==========================
How to install ``fitgrid``
==========================

These examples show how to install ``fitgrid`` into a new conda
working environment from the conda base environment with a shell
command in a linux or Mac terminal window.  They assume the ``conda``
and ``mamba`` executables are already installed in the base
environment and the users's channel configuration is the minconda3
default shown here:

.. code-block:: bash

   (base) $ which conda mamba
   /home/your_userid/miniconda3/bin/conda
   /home/your_userid/miniconda3/bin/mamba
   (base) $ conda config --show channels default_channels channel_priority
   channels:
     - defaults
   default_channels:
     - https://repo.anaconda.com/pkgs/main
     - https://repo.anaconda.com/pkgs/r
   channel_priority: flexible

.. note::

   The example installation commands are broken into separate lines for
   readability. If you do this, make sure the \\ is the last character on each line.
   Alternatively you can enter the command as a single line without any \\.

~~~~~~~~~~~~~~
with ``mamba``
~~~~~~~~~~~~~~

``fitgrid`` stable release
--------------------------

This is a typical installation of the latest stable release of
``fitgrid`` into a fresh conda environment named ``fg_012021``. This
pattern is likely to be compatible with recent versions of other conda
packages for x86_64 linux platforms and recent Intel Mac OSX.

.. code-block:: bash

   (base) $ mamba create --name fg_012021 \
       -c conda-forge -c ejolly -c kutaslab \
       fitgrid

.. note::

   This installation currently defaults to OpenBLAS builds of matrix
   math and linear algebra libraries so execution time on some Intel
   CPUs may be substantially longer than for the Intel Math
   Kernel (MKL) builds of the libraries. For a workaround see
   :ref:`mkl_v_openblas` below.


``fitgrid`` development version
-------------------------------

At times, the development version of ``fitgrid`` runs ahead of the latest
stable release and includes bug fixes and new features. The
latest development version may be installed by overriding the default
`kutaslab` conda channel with `kutaslab/label/pre-release` like so:

.. code-block:: bash

   (base) $ mamba create --name fg_012021 \
       -c conda-forge -c ejolly -c kutaslab/label/pre-release \
       fitgrid



Selecting a Python version
--------------------------

Specific versions of Python and other packages can be selected for
installation with the conda package specification syntax. This example
installs ``fitgrid`` with the most recent version of Python 3.8.

.. code-block:: bash

   (base) $ mamba create --name fg_012021 \
       -c conda-forge -c ejolly -c kutaslab \
       fitgrid python=3.8



.. _mkl_v_openblas:


       
Selecting MKL or OpenBLAS
-------------------------

On Intel CPUs, the `Intel Math Kernel Library (MKL)
<https://en.wikipedia.org/wiki/Math_Kernel_Library>`_ builds of
optimized math libraries like the Basic Linear Algebra Subprograms
(BLAS) may offer a substantial performance advantage over `OpenBLAS
<https://en.wikipedia.org/wiki/OpenBLAS>`_. For AMD CPUs OpenBLAS may
outperform MKL. This example shows how to enforce installation of the
MKL build and use ``conda list`` to inspect the installed packages.  To
select OpenBLAS builds, replace ``mkl`` with ``openblas`` in the first
command.

.. code-block:: bash

   (base) $ mamba create --name fg_012021 \
       -c conda-forge -c ejolly -c kutaslab \
       fitgrid "blas=*=mkl*"
   (base) $ activate fg_012021
   (fg_012021) $ conda list | egrep "(mkl|blas|liblapack)"
   # packages in environment at /home/userid/miniconda3/envs/fg_012021:
   blas                      2.109                       mkl    conda-forge
   blas-devel                3.9.0                     9_mkl    conda-forge
   libblas                   3.9.0                     9_mkl    conda-forge
   libcblas                  3.9.0                     9_mkl    conda-forge
   liblapack                 3.9.0                     9_mkl    conda-forge
   liblapacke                3.9.0                     9_mkl    conda-forge
   mkl                       2021.2.0           h06a4308_296  
   mkl-devel                 2021.2.0           h66538d2_296  
   mkl-include               2021.2.0           h06a4308_296  



Install fitgrid and run Examples Gallery notebooks
--------------------------------------------------
   
To run the notebooks in the :ref:`gallery` install `JupyterLab or
Jupyter <https://jupyter.org/>`_ in the same conda environment as
``fitgrid`` and it launch like so:

.. code-block:: bash

   (base) $ mamba create --name fg_012021 \
       -c conda-forge -c ejolly -c kutaslab \
       fitgrid jupyterlab
   (base) $ conda activate fg_012021
   (fg_012021) $ jupyter lab


Prioritize anaconda.org default channels over conda-forge
---------------------------------------------------------

This example shows how to install fitgrid into an environment
populated primarily with the stale-but-stable packages from the
Anaconda default channels. The explicit ``-c conda-forge`` channel is
necessary here because not all dependencies are available on the
default conda channels.

.. code-block:: bash

   (base) $ mamba create --name fg_012021 \
       -c defaults -c conda-forge -c ejolly -c kutaslab \
       fitgrid


~~~~~~~~~~~~~~
with ``conda``
~~~~~~~~~~~~~~

The ``conda`` installer may be used in place of ``mamba`` as shown in
the next example, although dependency resolution may be substantially
slower.


.. code-block:: bash

   (base) $ conda create --name fg_012021 \
       -c conda-forge -c ejolly -c kutaslab \
       fitgrid

.. note::

   The ``conda`` and ``mamba`` dependency resolution algorithms are not
   identical and may arrive at different solutions.


~~~~~~~~~~~~~~~~~~~~~~~~
``pip`` is not supported
~~~~~~~~~~~~~~~~~~~~~~~~

Since ``fitgrid`` requires numerous R packages, installing with the
Python package installer, :std:doc:`pip <pip:index>` is no longer
supported and is not recommended for general use.


===================
System requirements
===================

The platform of choice is linux. Minimum system requirements are not
known but obviously large scale regression modeling with millions of
data points is computationally demanding. Current versions of fitgrid
are developed and used in Ubuntu 20.04 running on a high-performance
multicore server with Intel CPUs (72 cores/144 threads, 1TB RAM);
continuous integrations tests run on ubuntu-latest and macos-10.15 on
GitHub Actions `hosted runners
<https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources>`_.
Previous versions of ``fitgrid`` were developed and used in CentOS 7
with Intel CPUs (24 cores/48 threads, 256-512 GB RAM). We are unable
to test the Windows 64-bit conda package, field reports are welcome,
see :ref:`Contributing <how_to_contribute>` for more information.

====
Tips
====

* Use ``conda list`` to inspect package versions and the channels they come
  from when constructing conda enviroments.

* To help avoid package version conflicts and speed up the dependency
  solver it can be useful to specify the Python version and install
  ``fitgrid`` along with the other conda packages you want into a
  fresh environment in one fell swoop. The package installers cannot
  see into the future. If packages are installed one by one, the next
  package version you want may not be compatible with what is already
  in the environment.

* ``mamba create`` and ``mamba install`` are not exact drop in
  replacements for ``conda create`` and ``conda install`` because
  ``conda`` has an affinity for packages on default conda channels and
  ``mamba`` has an affinity for packages on conda-forge and they may
  resolve dependencies differently.

* What works and what doesn't when creating conda environments and
  installing packages depends greatly on the *combinations* of
  packages you wish to install. Not all combinations of platforms,
  Python versions, installers, channel priority, and packages are
  compatible.

* Depending on your computer hardware, you may see a significant
  performance difference between the Intel MKL and OpenBLAS builds of
  the Basic Linear Algebra Support (BLAS) and Linear Algebra Package
  (LAPACK) libraries, particularly for fitting mixed-effects models.
#########
`fitgrid`
#########

``fitgrid`` is a Python package originally designed to streamline the
computation and interpretation of regression ERPs (rERPs) as described
by Smith and Kutas ([SmiKut2015]_). That report articulates the
conceptual foundation of rERPs as time-series of estimated regression
variable coefficients, :math:`\hat{\beta}_i`, for linear models of the
form :math:`\beta_0 + \beta_1 X_1 + \ldots + \beta_n X_n + e` and
their role in modeling EEG and other types of time-series data.

The ``fitgrid`` package is intended for researchers with a basic working
knowledge of scientific computing in Python. It provides access to the
multichannel time-series regression modeling computations with one
line of code and the familar ordinary least squares (OLS) and linear
mixed-effects regression (LMER) modeling formulas shared by Python and
R (``patsy`` [Smith2020]_; ``lm`` [R2020]_, ``lme4::lmer``
[BatesEtAl2015]_). The fit results across time and channels are
available on demand and returned as a tidy indexed
:py:class:`pandas.DataFrame` with one line of code and the same syntax
used to access results in a single fit object. These interfaces allow
researchers to conduct this type of modeling flexibly, efficiently,
informatively, and reproducibly with familiar tools and minimal
programming in Python data analysis workflows.

For a summary of the problem ``fitgrid`` solves, why it is worth
solving, and how it is solved, see the :ref:`about_fitgrid`. The
:ref:`user_guide` has information on specific topics including how the
OLS models are fit in Python ``statsmodels`` [SeaSkiPer2010]_ and the
LMER models are fit in R (``lme4::lmer``, ``lmerTest`` [KuzBroChr2017]_)
via ``pymer4`` [Jolly2018]_ and ``rpy2`` [Gautier2021]_ under the
hood. The :ref:`gallery` contains executable ``fitgrid`` vignettes
with simulated data, experimental EEG recordings, and NOAA tide and
atmospheric observations. The :ref:`workflow` example illustrates the
``fitgrid`` modeling pipeline. The examples can be downloaded as
Python scripts or Jupyter notebooks thanks to :std:doc:`sphinx-gallery
<sphinx-gallery:index>`. The :ref:`api` is a complete listing of
``fitgrid`` classes, methods, attributes, and functions.  The
:ref:`bibliography` includes :ref:`references` and
:ref:`fitgrid_reports`.

.. image:: ../source/_static/fitgrid_overview.png

.. toctree::
    :hidden:
    :glob:

    about_fitgrid
    installation
    auto_gallery/quickstart
    user_guide
    gallery
    api
    bibliography
    contributing



.. _api:

###
API
###


.. _quick_reference:

===============
Quick reference
===============

.. toctree::
   :maxdepth: 2

   quick_reference

=======
Indexes
=======

* :ref:`Alphabetical Index <genindex>`
* :ref:`Module Index <modindex>`


fitgrid.errors module
=====================

.. automodule:: fitgrid.errors
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.utils.lm module
=======================

.. automodule:: fitgrid.utils.lm
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.models module
=====================

.. automodule:: fitgrid.models
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.sample\_data module
===========================

.. automodule:: fitgrid.sample_data
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.tools module
====================

.. automodule:: fitgrid.tools
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.fitgrid module
======================

.. automodule:: fitgrid.fitgrid
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.plots module
====================

.. automodule:: fitgrid.plots
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.fake\_data module
=========================

.. automodule:: fitgrid.fake_data
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.utils.summary module
============================

.. automodule:: fitgrid.utils.summary
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.utils.lmer module
=========================

.. automodule:: fitgrid.utils.lmer
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.utils package
=====================

Submodules
----------

.. toctree::
   :maxdepth: 4

   fitgrid.utils.lm
   fitgrid.utils.lmer
   fitgrid.utils.summary

Module contents
---------------

.. automodule:: fitgrid.utils
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.epochs module
=====================

.. automodule:: fitgrid.epochs
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid package
===============

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   fitgrid.utils

Submodules
----------

.. toctree::
   :maxdepth: 4

   fitgrid.defaults
   fitgrid.epochs
   fitgrid.errors
   fitgrid.fake_data
   fitgrid.fitgrid
   fitgrid.io
   fitgrid.models
   fitgrid.plots
   fitgrid.sample_data
   fitgrid.tools

Module contents
---------------

.. automodule:: fitgrid
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.defaults module
=======================

.. automodule:: fitgrid.defaults
   :members:
   :undoc-members:
   :show-inheritance:
fitgrid.io module
=================

.. automodule:: fitgrid.io
   :members:
   :undoc-members:
   :show-inheritance:
:orphan:

Examples
========
:orphan:

Epochs data
===========

:orphan:

Model fitting
=============
:orphan:

Model evaluation
================



