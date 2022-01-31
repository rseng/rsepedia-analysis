[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Documentation Status](https://readthedocs.org/projects/esmvaltool/badge/?version=latest)](https://esmvaltool.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3401363.svg)](https://doi.org/10.5281/zenodo.3401363)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ESMValGroup?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![CircleCI](https://circleci.com/gh/ESMValGroup/ESMValTool/tree/main.svg?style=svg)](https://circleci.com/gh/ESMValGroup/ESMValTool/tree/main)
[![example branch parameter](https://github.com/github/docs/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/ESMValGroup/ESMValTool/actions)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/79bf6932c2e844eea15d0fb1ed7e415c)](https://www.codacy.com/gh/ESMValGroup/ESMValTool?utm_source=github.com&utm_medium=referral&utm_content=ESMValGroup/ESMValTool&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/79bf6932c2e844eea15d0fb1ed7e415c)](https://www.codacy.com/gh/ESMValGroup/ESMValTool?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ESMValGroup/ESMValTool&amp;utm_campaign=Badge_Grade)
[![Docker Build Status](https://img.shields.io/docker/cloud/build/esmvalgroup/esmvaltool.svg)](https://hub.docker.com/r/esmvalgroup/esmvaltool/)
[![Anaconda-Server Badge](https://anaconda.org/esmvalgroup/esmvaltool/badges/installer/conda.svg)](https://conda.anaconda.org/esmvalgroup)

![esmvaltoollogo](https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/figures/ESMValTool-logo-2.png)

- [**Documentation**](https://docs.esmvaltool.org/en/latest/)
- [**ESMValTool Website**](https://www.esmvaltool.org/)
- [**ESMValTool Tutorial**](https://esmvalgroup.github.io/ESMValTool_Tutorial/index.html)
- [**ESMValGroup Project on GitHub**](https://github.com/ESMValGroup)
- [**Gallery**](https://docs.esmvaltool.org/en/latest/gallery.html)

# Introduction

ESMValTool is a community-developed climate model diagnostics and evaluation software package, driven
both by computational performance and scientific accuracy and reproducibility. ESMValTool is open to both
users and developers, encouraging open exchange of diagnostic source code and evaluation results from the
Coupled Model Intercomparison Project [CMIP](https://www.wcrp-climate.org/wgcm-cmip) ensemble. For a
comprehensive introduction to ESMValTool please visit our
[documentation](https://docs.esmvaltool.org/en/latest/introduction.html) page.

# Running esmvaltool

Diagnostics from ESMValTool are run using [recipe](https://docs.esmvaltool.org/en/latest/recipes/index.html)
files that contain pointers to the requested data types, directives for the preprocessing steps that data
will be subject to, and directives for the actual diagnostics that will be run with the now preprocessed data.
Data preprocessing is done via the [ESMValCore](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/index.html) package, a pure Python, highly-optimized scientific library, developed by the ESMValTool core developers,
and that performs a number of common analysis tasks
such as regridding, masking, levels extraction etc. [Diagnostics](https://docs.esmvaltool.org/en/latest/develop/diagnostic.html) are written in a variety of programming languages (Python, NCL, R, Julia) and are developed by the wider
scientific community, and included after a scientific and technical review process.

# Input data

ESMValTool can run with the following types of data as input:

- CMIP5
- CMIP6
- OBS, OBS6
- obs4MIPs
- ana4mips
- CORDEX

# Getting started

Please see [getting started](https://docs.esmvaltool.org/en/latest/quickstart/index.html) on readthedocs as well as [ESMValTool tutorial](https://esmvalgroup.github.io/ESMValTool_Tutorial/index.html). The tutorial is a set of lessons that together teach skills needed to work with ESMValTool in climate-related domains.

## Getting help

The easiest way to get help if you cannot find the answer in the documentation on [readthedocs](https://docs.esmvaltool.org), is to open an [issue on GitHub](https://github.com/ESMValGroup/ESMValTool/issues).

## Contributing

If you would like to contribute a new diagnostic or feature, please have a look at our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/index.html).
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
reported by contacting the project team at birgit.hassler@dlr.de or
alistair.sellar@metoffice.gov.uk. All
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
# Contributions are very welcome

Please read our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/index.html).
<!--
    Thank you for contributing to our project!

    Please do not delete this text completely, but read the text below and keep
    items that seem relevant. If in doubt, just keep everything and add your
    own text at the top, a reviewer will update the checklist for you.

    While the checklist is intended to be filled in by the technical and scientific
    reviewers, it is the responsibility of the author of the pull request to make
    sure all items on it are properly implemented.
-->

## Description

<!--
    Please describe your changes here, especially focusing on why this pull request makes
    ESMValTool better and what problem it solves.

    Before you start, please read our contribution guidelines: https://docs.esmvaltool.org/en/latest/community/

    Please fill in the GitHub issue that is closed by this pull request, e.g. Closes #1903
-->
- Closes #issue_number
- Link to documentation:

* * *

## Before you get started

<!--
    Please discuss your idea with the development team before getting started,
    to avoid disappointment or unnecessary work later. The way to do this is
    to open a new issue on GitHub.
-->

- [ ] [☝ Create an issue](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#contributing-code-and-documentation) to discuss what you are going to do

## Checklist

It is the responsibility of the author to make sure the pull request is ready to review. The icons indicate whether the item will be subject to the [🛠 Technical][1] or [🧪 Scientific][2] review.

<!-- The next two lines turn the 🛠 and 🧪 below into hyperlinks -->
[1]: https://docs.esmvaltool.org/en/latest/community/review.html#technical-review
[2]: https://docs.esmvaltool.org/en/latest/community/review.html#scientific-review

- [ ] [🛠][1] This pull request has a [descriptive title](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#pull-request-title)
- [ ] [🛠][1] Code is written according to the [code quality guidelines](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#code-quality)
- [ ] [🛠][1] [Documentation](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#documentation) is available
- [ ] [🛠][1] [Tests](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#tests) run successfully
- [ ] [🛠][1] The [list of authors](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#list-of-authors) is up to date
- [ ] [🛠][1] Any changed dependencies have been [added or removed correctly](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#dependencies)
- [ ] [🛠][1] All [checks below this pull request](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#pull-request-checks) were successful

### [New or updated recipe/diagnostic](https://docs.esmvaltool.org/en/latest/community/diagnostic.html)

- [ ] [🧪][2] [Recipe runs successfully](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#testing-recipes)
- [ ] [🧪][2] [Recipe is well documented](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recipe-and-diagnostic-documentation)
- [ ] [🧪][2] [Figure(s) and data](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#diagnostic-output) look as expected from literature
- [ ] [🛠][1] [Provenance information](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recording-provenance) has been added

### [New or updated data reformatting script](https://docs.esmvaltool.org/en/latest/develop/dataset.html)

- [ ] [🛠][1] [Documentation](https://docs.esmvaltool.org/en/latest/community/dataset.html#dataset-documentation) is available
- [ ] [🛠][1] The dataset has been [added to the CMOR check recipe](https://docs.esmvaltool.org/en/latest/community/dataset.html#testing)
- [ ] [🧪][2] Numbers and units of the data look [physically meaningful](https://docs.esmvaltool.org/en/latest/community/dataset.html#scientific-sanity-check)

***

To help with the number of pull requests:

-  🙏 We kindly ask you to [review](https://docs.esmvaltool.org/en/latest/community/review.html#review-of-pull-requests) two other [open pull requests](https://github.com/ESMValGroup/ESMValTool/pulls) in this repository

<!--
If you need help with any of the items on the checklists above, please do not hesitate to ask by commenting in the issue or pull request.
-->
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is. If you are developing a new diagnostic script, please provide a link to the code/branch on GitHub that you are working in.

**Please attach**
  - The recipe that you are trying to run, you can find a copy in the `run` directory in the output directory
  - The `main_log_debug.txt` file, this can also be found in the `run` directory in the output directory
---
name: New diagnostic
about: Develop a new diagnostic.
title: ''
labels: diagnostic
assignees: ''

---

**Short description of the diagnostic**
Add a short description of the diagnostic that you would like to add.

**Branch and pull request**
Once you've started working, add the branch (and pull request)
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is.

**Would you be able to help out?**
Would you have the time and skills to implement the solution yourself?
ESMVal python library
=====================

General structure
-----------------

This directory contains usefull scripts and tools for ESMVal which are to support the writing of diagnostics in a clean and short way.

In general you have two options to include code into this library.

A) you have some small code snippet that is needed somewhere else (e.g. a special function).
In this case you should put your code into either a new file in the 'python' directory or an already existing file.
Let's say that you have developed some statistics test functionality, then it would be good to put it in a file statistics.py

However, more importantly it is important to recognize that python is in general organized in modules (http://docs.python.org/2/tutorial/modules.html#more-on-modules).
Modules shall ever been used if you are working with bigger projects or code which encapsulated some functionality in different code parts (e.g. a bigger statistical library, a plotting library ...)

B) Modules are simply organized in directories which and contain an additional __init__.py file which tells the python interpreter which files are part of the module.
Using modules allows you to load entire packages as a whole or parts of it into your code using e.g.

from scipy import stats

Thus the ESMVal python library is supposed to look like e.g.

/lib+
    |_python+
            |_my_cool_collection_of_snippets.py
            |_plots.py
            |_statlib+
                     |___init__.py
                     |_linear_regressions.py
                     |_pdf.py
                     |_stattests.py


Further reading
---------------
For further reading reagarding modules in python look at the following links:

* http://docs.python.org/2/tutorial/modules.html
Contributing to the climdex.pcic.ncdf R package
==========================================

Getting Started
---------------

- Create a `Github account`_.
- Fork the repository on Github at https://github.com/pacificclimate/climdex.pcic.ncdf.
- Work on the code (see the `next section`_)
- Send us a `pull request`_.

.. _Github account: https://github.com/signup/free
.. _pull request: https://help.github.com/articles/using-pull-requests/
.. _next section: #how-to-set-up-a-development-environment

How to set up a development environment
---------------------------------------

You don't need much to get started for development. You'll need to have installed:

- R (ensure that all of the "Depends", "Imports", and "Suggests" pacakges are also installed)
- Any C++ build environment supported by the `CRAN package checking`_
- git
- your text editor of choice

That's it!

Once you have the required software installed, create a local clone of the repository.
::
    $ git clone https://github.com/[your_user]/climdex.pcic.ncdf.git

Build the docs (which builds the auto-generated NAMESPACE file needed to build). See `below <#how-to-build-the-docs>`_.

Then make sure that everything builds out of the box
::
    $ R CMD build climdex.pcic.ncdf/

.. _CRAN package checking: http://cran.r-project.org/web/checks/check_flavors.html

How to run the tests
--------------------

Running the tests can be done with one command:
::
    james@basalt ~/code/git $ R CMD check climdex.pcic.ncdf/

You'll see a bunch of package building spew that has nothing to do with the tests. But towards the end, you see something like this:
::
   * checking for unstated dependencies in tests ... OK
   * checking tests ...
       Running ‘bootstrap.R’
       Running ‘test_basic_file_funcs.R’
       Running ‘test_var_meta.R’
    OK

Bug reports
-----------

If there are problems with our package or bugs in the code, please let us know! We welcome bug reports. To submit one:

- `Create a new issue`_ on our GitHub page.
- Tag/label the issue as a bug
- Leave it unassigned

Then please follow these guidelines for writing your report:

- Please describe in as much detail as possible
- Include a complete description of:

  - Exactly what you did (i.e. "steps to reproduce")
  - What you expected to happen?
  - What did happen?

- Include *all* output from the terminal.
- Run R's ``sessionInfo()`` function and include the full output.

I cannot stress enough how important it is to contrast what you expected to happen, with what actually happened. When executing the code does not produce the *advertised* result, there is a bug in the package. When the code does not produce the result that you *wished* it had, this is *not* a bug. We receive far too many reports in the latter category.

.. _Create a new issue: https://github.com/pacificclimate/climdex.pcic.ncdf/issues/new

.. _build-the-docs:

How to build the docs
---------------------

The package documentation is inline in the code. All of the manual pages are built by using ``roxygen2``. Make sure that you have ``roxygen2`` installed and loaded:
::
   james@basalt ~/code/git/climdex.pcic.ncdf $ R

   R version 3.0.3 (2014-03-06) -- "Warm Puppy"
   Copyright (C) 2014 The R Foundation for Statistical Computing
   Platform: x86_64-pc-linux-gnu (64-bit)

   R is free software and comes with ABSOLUTELY NO WARRANTY.
   You are welcome to redistribute it under certain conditions.
   Type 'license()' or 'licence()' for distribution details.

   Natural language support but running in an English locale

   R is a collaborative project with many contributors.
   Type 'contributors()' for more information and
   'citation()' on how to cite R or R packages in publications.

   Type 'demo()' for some demos, 'help()' for on-line help, or
   'help.start()' for an HTML browser interface to help.
   Type 'q()' to quit R.

   > library(roxygen2)

Then call ``roxygenize()`` to build the docs.
::
   > roxygenize()
   First time using roxygen2 4.0. Upgrading automatically...
   Loading required package: PCICt
   Loading required package: ncdf4
   Loading required package: climdex.pcic
   Loading required package: ncdf4.helpers
   Loading required package: snow
   Loading required package: udunits2
   Loading required package: functional
   Loading required package: proj4
   Writing NAMESPACE
   Writing climdex.pcic.ncdf.Rd
   Writing create.climdex.cmip5.filenames.Rd
   Writing get.climdex.variable.list.Rd
   Writing get.climdex.functions.Rd
   Writing get.climdex.variable.metadata.Rd
   Writing create.ncdf.output.files.Rd
   Writing compute.climdex.indices.Rd
   Writing flatten.dims.Rd
   Writing get.data.Rd
   Writing get.northern.hemisphere.booleans.Rd
   Writing get.quantiles.object.Rd
   Writing compute.indices.for.stripe.Rd
   Writing get.thresholds.chunk.Rd
   Writing write.climdex.results.Rd
   Writing get.quantiles.for.stripe.Rd
   Writing create.thresholds.file.Rd
   Writing get.var.file.idx.Rd
   Writing create.file.metadata.Rd
   Writing get.thresholds.metadata.Rd
   Writing create.thresholds.from.file.Rd
   Writing thresholds.open.Rd
   Writing thresholds.close.Rd
   Writing create.indices.from.files.Rd


Submitting pull requests
------------------------

We would love help from the greater climate community in developing the package and we welcome contributions to climdex.pcic.ncdf package.

- Please write tests for any functionality that you may add.
- Please modify tests for any functionality that you change.
- In short, please make sure that all of the tests pass.

After you are *positive* that everything is completely tested with passing test suite, we would love to see your pull request. If you are not familiar with the process, please follow the GitHub's help page for submitting `pull request`_.

Don't code? No problem!
-----------------------

Even if you don't program for a living there are plenty of ways to help. Not only is the code open and collaborative, but so is the documentation and issue tracking. Anyone can help with these. If you can't program, consider helping with the following:

- If the documentation doesn't answer your questions, it probably doesn't answer many people's questions. Help us all out and write something that does.
- Take a look through the outstanding `"help wanted" issues`_, and see if you know any of the answers.
- If there are `open bug reports`_, see if you can reproduce the problem and verify that it exists. Having bug reports validated and/or clarified by multiple parties is extremely valuable.
- Tell us your story. If ``climdex.pcic.ncdf`` has helped your project to better understand climate extremes, we would love to hear about it. Write a blog post and/or send an e-mail to the `package maintainer`_.

.. _"help wanted" issues: https://github.com/pacificclimate/climdex.pcic.ncdf/labels/help%20wanted
.. _open bug reports: https://github.com/pacificclimate/climdex.pcic.ncdf/labels/bug
.. _package maintainer: mailto:hiebert@uvic.ca
What is climdex.pcic.ncdf?
=====================

* `climdex.pcic.ncdf` is a companion library for `climdex.pcic` which helps in using NetCDF input grids and writing to NetCDF output files when computing the `27 core indices of extreme climate`_. The code allows for parallel computation of indices using either a SOCK or MPI cluster. It was written for the `R statistical programming language`_ by the `Pacific Climate Impacts Consortium`_.

.. _27 core indices of extreme climate: http://etccdi.pacificclimate.org/list_27_indices.shtml
.. _R statistical programming language: http://www.r-project.org/
.. _Pacific Climate Impacts Consortium: http://pacificclimate.org/

Getting Help
============

New to programming or to R?
---------------------------

* Read the the `Software Carpentry`_  `Programming in R`_ lessons
* Read one of the man `R Manuals`_.
* Attend an `R Users Group`_ meeting.

.. _Software Carpentry: http://software-carpentry.org/index.html
.. _Programming in R: http://software-carpentry.org/v5/novice/r/index.html
.. _R Manuals: http://cran.r-project.org/manuals.html
.. _R Users Group: http://r-users-group.meetup.com/

Looking for code?
-----------------

* Get the latest `climdex.pcic.ncdf release from our website`_.
* Explore the `development repository`_.
* Install it with devtools ::

    > library(devtools)
    > install_github('pacificclimate/climdex.pcic.ncdf', ref='release')

.. _climdex.pcic.ncdf release from our website: http://www.pacificclimate.org/sites/default/files/climdex.pcic_.ncdf_0.5-4.tar_.gz
.. _development repository: https://github.com/pacificclimate/climdex.pcic.ncdf/

Need help using the package?
----------------------------

* Read the manual ::

    > library(climdex.pcic.ncdf)
    Loading required package: PCICt
    > ?climdex.pcic.ncdf

* Create a `new issue`_ on the `package issue tracker`_ and label it "help wanted"[1]_.

.. _new issue: https://github.com/pacificclimate/climdex.pcic.ncdf/issues/new

Want to contribute?
-------------------

* To report a bug in pcic.climdex use the `package issue tracker`_ (after you've read the `bug reporting guide`_).
* To help with development read through the `contributor's guide`_

.. _bug reporting guide: https://github.com/pacificclimate/climdex.pcic.ncdf/blob/master/CONTRIBUTING.rst#bug-reports
.. _package issue tracker: https://github.com/pacificclimate/climdex.pcic.ncdf/issues
.. _contributor's guide: https://github.com/pacificclimate/climdex.pcic.ncdf/blob/master/CONTRIBUTING.rst

Still need help?
----------------

* Contact climate@uvic.ca and let us know what we can do.

.. [1] Please know that the pool of people who can provide support for the package is extremely small and time is limited.  We don't necessarily have the capacity for long, open-ended user support. If you keep your questions short, specific and direct, there's a greater probability that someone will take on the ticket.
.. _utils:

Utilities
*********

This section provides information on tools that are useful when developing
ESMValTool.
Tools that are specific to ESMValTool live in the
`esmvaltool/utils <https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/utils>`_
directory, while others can be installed using the usual package managers.

.. _pre-commit:

Pre-commit
==========

`pre-commit <https://pre-commit.com/>`__ is a handy tool that can run many
tools for checking code quality with a single command.
Usually it is used just before committing, to avoid accidentally committing
mistakes.
It knows knows which tool to run for each filetype, and therefore provides
a convenient way to check your code!


To run ``pre-commit`` on your code, go to the ESMValTool directory
(``cd ESMValTool``) and run

::

   pre-commit run

By default, pre-commit will only run on the files that have been changed,
meaning those that have been staged in git (i.e. after
``git add your_script.py``).

To make it only check some specific files, use

::

   pre-commit run --files your_script.py

or

::

   pre-commit run --files your_script.R

Alternatively, you can configure ``pre-commit`` to run on the staged files before
every commit (i.e. ``git commit``), by installing it as a `git hook <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>`__ using

::

   pre-commit install

Pre-commit hooks are used to inspect the code that is about to be committed. The
commit will be aborted if files are changed or if any issues are found that
cannot be fixed automatically. Some issues cannot be fixed (easily), so to
bypass the check, run

::

   git commit --no-verify

or

::

   git commit -n

or uninstall the pre-commit hook

::

   pre-commit uninstall


Note that the configuration of pre-commit lives in
`.pre-commit-config.yaml <https://github.com/ESMValGroup/ESMValTool/blob/main/.pre-commit-config.yaml>`_.

.. _nclcodestyle:

nclcodestyle
============

A tool for checking the style of NCL code, based on pycodestyle.
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash

    nclcodestyle /path/to/file.ncl

.. _recipe_test_tool:

Colormap samples
================
Tool to generate colormap samples for ESMValTool's default Python and NCL colormaps.

Run

.. code-block:: bash

    esmvaltool colortables python

or

.. code-block:: bash

    esmvaltool colortables ncl

to generate the samples.

Testing recipes
===============

Tools for testing recipes.

Run all recipes
===============

A cylc suite for running all recipes is available in
`esmvaltool/utils/testing/regression <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/utils/testing/regression>`__

To prepare for using this tool:

#. Log in to Mistral or another system that uses slurm
#. Make sure the required CMIP and observational datasets are available and configured in config-user.yml
#. Make sure the required auxiliary data is available (see recipe documentation)
#. Install ESMValTool
#. Update config-user.yml so it points to the right data locations

Next, get started with `cylc <https://cylc.github.io/cylc-doc/stable/html/tutorial.html>`_:

#. Run ``module load cylc``
#. Register the suite with cylc ``cylc register run-esmvaltool-recipes ~/ESMValTool/esmvaltool/utils/testing/regression``
#. Edit the suite if needed, this allows e.g. choosing which recipes will be run
#. Validate the suite ``cylc validate run-esmvaltool-recipes --verbose``, this will e.g. list the recipes in the suite
#. Run all recipes ``cylc run run-esmvaltool-recipes``
#. View progress ``cylc log run-esmvaltool-recipes``, use e.g. ``cylc log run-all-esmvaltool-recipes examples-recipe_python_yml.1 --stdout`` to see the log of an individual esmvaltool run. Once the suite has finished running, you will see the message "WARNING - suite stalled" in the log.
#. Stop the cylc run once everything is done ``cylc stop run-esmvaltool-recipes``.
#. Create the index.html overview page by running ``python esmvaltool/utils/testing/regression/summarize.py ~/esmvaltool_output/``

Test recipe settings
--------------------

A tool for generating recipes with various diagnostic settings, to test of those work.
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash

    test_recipe --help


.. _draft_release_notes.py:

draft_release_notes.py
======================

`draft_release_notes.py <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/utils/draft_release_notes.py>`__
is a script for drafting release notes based on the titles and labels of
the GitHub pull requests that have been merged since the previous release.

To use it, install the package pygithub_:

.. code-block:: bash

   pip install pygithub

Create a `GitHub access token`_ (leave all boxes for additional
permissions unchecked) and store it in the file ``~/.github_api_key``.

Edit the script and update the date and time of the previous release and run
the script:

.. code-block:: bash

   python esmvaltool/utils/draft_release_notes.py ${REPOSITORY}

``REPOSITORY`` can be either ``esmvalcore`` or ``esmvaltool`` depending on the
release notes you want to create.

Review the resulting output (in ``.rst`` format) and if anything needs changing,
change it on GitHub and re-run the script until the changelog looks acceptable.
In particular, make sure that pull requests have the correct label, so they are
listed in the correct category.
Finally, copy and paste the generated content at the top of the changelog.

Converting Version 1 Namelists to Version 2 Recipes
===================================================

The
`xml2yml <https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/utils/xml2yml>`_
converter can turn the old xml namelists into new-style yml
recipes. It is implemented as a xslt stylesheet that needs a processor
that is xslt 2.0 capable. With this, you simply process your old
namelist with the stylesheet xml2yml.xsl to produce a new yml recipe.

After the conversion you need to manually check the mip information in
the variables! Also, check the caveats below!

Howto
-----

One freely available processor is the Java based
`saxon <http://saxon.sourceforge.net/>`__. You can download the free he
edition
`here <https://sourceforge.net/projects/saxon/files/latest/download?source=files>`__.
Unpack the zip file into a new directory. Then, provided you have Java
installed, you can convert your namelist simply with:

::

   java -jar $SAXONDIR/saxon9he.jar -xsl:xml2yml.xsl -s:namelist.xml -o:recipe.yml

Caveats/Known Limitations
-------------------------

-  At the moment, not all model schemes (OBS, CMIP5, CMIP5_ETHZ…) are
   supported. They are, however, relatively easy to add, so if you need
   help adding a new one, please let me know!
-  The documentation section (namelist_summary in the old file) is not
   automatically converted.
-  In version 1, one could name an exclude, similar to the reference
   model. This is no longer possible and the way to do it is to include
   the models with another ``additional_models`` tag in the variable
   section. That conversion is not performed by this tool.

Authored by **Klaus Zimmermann**, direct questions and comments to
klaus.zimmermann@smhi.se

.. _GitHub access token: https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line
.. _pygithub: https://pygithub.readthedocs.io/en/latest/introduction.html


Recipe filler
=============

If you need to fill in a blank recipe with additional datasets, you can do that with
the command `recipe_filler`. This runs a tool to obtain a set of additional datasets when
given a blank recipe, and you can give an arbitrary number of data parameters. The blank recipe
should contain, to the very least, a list of diagnostics, each with their variable(s).
Example of running the tool:

.. code-block:: bash

    recipe_filler recipe.yml

where `recipe.yml` is the recipe that needs to be filled with additional datasets; a minimal
example of this recipe could be:

.. code-block:: yaml

    diagnostics:
      diagnostic:
        variables:
          ta:
            mip: Amon  # required
            start_year: 1850  # required
            end_year: 1900  # required


Key features
------------

- you can add as many variable parameters as are needed; if not added, the
  tool will use the ``"*"`` wildcard and find all available combinations;
- you can restrict the number of datasets to be looked for with the ``dataset:``
  key for each variable, pass a list of datasets as value, e.g.
  ``dataset: [MPI-ESM1-2-LR, MPI-ESM-LR]``;
- you can specify a pair of experiments, e.g. ``exp: [historical, rcp85]``
  for each variable; this will look for each available dataset per experiment
  and assemble an aggregated data stretch from each experiment to complete
  for the total data length specified by ``start_year`` and ``end_year``; equivalent to
  ESMValTool's syntax on multiple experiments; this option needs an ensemble
  to be declared explicitly; it will return no entry if there are gaps in data;
- ``start_year`` and ``end_year`` are required and are used to filter out the
  datasets that don't have data in the interval; as noted above, the tool will not
  return datasets with partial coverage from ``start_year`` to ``end_year``;
  if you want all possible years hence no filtering on years just use ``"*"``
  for start and end years;
- ``config-user: rootpath: CMIPX`` may be a list, rootpath lists are supported;
- all major DRS paths (including ``default``, ``BADC``, ``ETHZ`` etc) are supported;
- speedup is achieved through CMIP mip tables lookup, so ``mip`` is required in recipe;

Caveats
-------

- the tool doesn't yet work with derived variables; it will not return any available datasets;
- operation restricted to CMIP data only, OBS lookup is not available yet.


Extracting a list of input files from the provenance
====================================================

There is a small tool available to extract just the list of input files used to generate
a figure from the ``*_provenance.xml`` files (see :ref:`recording-provenance` for more
information).

To use it, install ESMValTool from source and run

.. code-block:: bash

    python esmvaltool/utils/prov2files.py /path/to/result_provenance.xml

The tool is based on the `prov <https://prov.readthedocs.io/en/latest/readme.html>`_
library, a useful library for working with provenance files.
With minor adaptations, this script could also print out global attributes
of the input NetCDF files, e.g. the tracking_id.
.. _faq:

Frequently Asked Questions
**************************

Is there a mailing list?
========================

Yes, you can subscribe to the ESMValTool user mailing list and join the discussion on general topics (installation, configuration, etc). See :ref:`mailing-list`.

What is YAML?
=============

While ``.yaml`` or ``.yml`` is a relatively common format, users may not have
encountered this language before. The key information about this format is:

- yaml is a human friendly markup language;
- yaml is commonly used for configuration files (gradually replacing the
  venerable ``.ini``);
- the syntax is relatively straightforward;
- indentation matters a lot (like ``Python``)!
- yaml is case sensitive;

More information can be found in the `yaml tutorial
<https://learnxinyminutes.com/docs/yaml/>`_ and `yaml quick reference card
<https://yaml.org/refcard.html>`_. ESMValTool uses the `yamllint
<http://www.yamllint.com>`_ linter tool to check recipe syntax.


.. _rerunning:

Re-running diagnostics
======================

If a diagnostic fails, you will get the message

.. code:: bash

   INFO    To re-run this diagnostic script, run:

If you run the command in the stdout you will be able to re-run the
diagnostic without having to re-run the whole preprocessor. If you add the ``-f``
argument (available only for Python diagnostics, check your options with ``--help``)
that will force an overwrite, and it will delete not just the failed diagnostic,
but the contents of its ``work_dir`` and ``plot_dir`` directories - this is useful when needing to
redo the whole work. Adding ``-i`` or ``--ignore-existing`` will not delete any existing files,
and it can be used to skip work that was already done successfully, provided
that the diagnostic script supports this.


Enter interactive mode with iPython
===================================

Sometimes it is useful to enter an interactive session to have a look what's going on.
Insert a single line in the code where you want to enter IPython:
``import IPython; IPython.embed()``

This is a useful functionality because it allows the user to `fix` things on-the-fly and after
quitting the Ipython console, code execution continues as per normal.


Use multiple config-user.yml files
==================================

The user selects the configuration yaml file at run time. It's possible to
have several configurations files. For instance, it may be practical to have one
config file for debugging runs and another for production runs.

Create a symbolic link to the latest output directory
=====================================================

When running multiple times the same recipe, the tool creates separate output directories
sorted by the time tag that they were created at; sometimes, when running quite a few times,
it is not straightforward to detect which one is the `latest` output directory, so a symbolic
link attached to it would make things more clear e.g.:

.. code:: bash

   recipe_example_20190905_163431
   recipe_example_20190905_163519
   recipe_example_latest -> recipe_example_20190905_163519


You can do that by running the tool using the latest output as basis and creating
a symbolic link to it so it gets picked up at every re-run iteration:

.. code:: bash

   esmvaltool run recipe_example.yml; \
   ln -sfT $(ls -1d ~/esmvaltool_output/recipe_example_* | tail -1) ~/esmvaltool_output/recipe_example_latest


.. uncomment when feature plopped in main
.. # Running a dry run
.. =================

.. You can run in dry-run mode with

.. .. code:: bash

..   esmvaltool run recipe_xxx.yml --dry-run


.. This mode activated will run through the data finding and CMOR checks and fixes
.. and will highlight on screen and in `run/main_log.txt` every time certain data is
.. missing or there are issues with the CMOR checks; note that no data is written
.. to disk and no diagnostics are run; you don't have to modify your recipe in any
.. way to have this mode run. The information provided will help you obtain any data
.. that is missing and/or create fixes for the datasets and variables that failed the
.. CMOR checks and could not be fixed on the fly.
.. ESMValTool documentation master file, created by
   sphinx-quickstart on Tue Jun  2 11:34:13 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ESMValTool's documentation!
======================================

.. include:: _sidebar.rst.inc

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

Introduction
************

About
=====

The Earth System Model Evaluation Tool (ESMValTool) is a
community-development that aims at improving diagnosing and
understanding of the causes and effects of model biases and inter-model
spread. The ESMValTool is open to both users and developers encouraging
open exchange of diagnostic source code and evaluation results from the
Coupled Model Intercomparison Project (CMIP) ensemble. This will
facilitate and improve ESM evaluation beyond the state-of-the-art and
aims at supporting the activities within CMIP and at individual
modelling centers. We envisage running the ESMValTool routinely on the
CMIP model output utilizing observations available through the Earth
System Grid Federation (ESGF) in standard formats (obs4MIPs) or made
available at ESGF nodes.

The goal is to develop a benchmarking and evaluation tool that produces
well-established analyses as soon as model output from CMIP simulations
becomes available, e.g., at one of the central repositories of the ESGF.
This is realized through standard recipes that reproduce a certain set
of diagnostics and performance metrics that have demonstrated its
importance in benchmarking Earth System Models (ESMs) in a paper or
assessment report, such as Chapter 9 of the Intergovernmental Panel on
Climate Change (IPCC) Fifth Assessment Report (AR5) (Flato et al.,
2013). The expectation is that in this way a routine and systematic
evaluation of model results can be made more efficient, thereby enabling
scientists to focus on developing more innovative methods of analysis
rather than constantly having to "reinvent the wheel".

In parallel to standardization of model output, the ESGF also hosts
observations for Model Intercomparison Projects (obs4MIPs) and
reanalyses data (ana4MIPs). obs4MIPs provides open access data sets of
satellite data that are comparable in terms of variables, temporal and
spatial frequency, and periods to CMIP model output (Taylor et al.,
2012). The ESMValTool utilizes these observations and reanalyses from
ana4MIPs plus additionally available observations in order to evaluate
the models performance. In many diagnostics and metrics, more than one
observational data set or meteorological reanalysis is used to assess
uncertainties in observations.

The main idea of the ESMValTool is to provide a broad suite of
diagnostics which can be performed easily when new model simulations are
run. The suite of diagnostics needs to be broad enough to reflect the
diversity and complexity of Earth System Models, but must also be robust
enough to be run routinely or semi-operationally. In order the address
these challenging objectives the ESMValTool is conceived as a framework
which allows community contributions to be bound into a coherent
framework.

.. _contact:

Contact
=======

See `www.esmvaltool.org <https://www.esmvaltool.org>`_ for general contact information.

.. _mailing-list:

User mailing list
-----------------

Subscribe to the ESMValTool announcements mailing list
`esmvaltool@listserv.dfn.de <mailto:esmvaltool@listserv.dfn.de>`__
to stay up to date about new releases, monthly online meetings, upcoming workshops, and trainings.

To subscribe, send an email to
`sympa@listserv.dfn.de <mailto:sympa@listserv.dfn.de?subject=subscribe%20esmvaltool>`_
with the following subject line:

-  *subscribe esmvaltool*

or

-  *subscribe esmvaltool YOUR_FIRSTNAME YOUR_LASTNAME*

The mailing list also has a `public archive <https://www.listserv.dfn.de/sympa/arc/esmvaltool>`_ online.

.. _discussions_page:

Discussions page
----------------

The `ESMValTool Discussions page <https://github.com/ESMValGroup/ESMValTool/discussions>`__
is open for all general and technical questions on the ESMValTool: installation, application, development, or any other question or comment you may have.

.. _core-team:

Monthly meetings
----------------

We have monthly online meetings using `zoom <https://zoom.us/>`__, anyone with
an interest in the ESMValTool is welcome to join these meetings to connect with
the community.
These meetings are always announced in an issue
on the `ESMValTool <https://github.com/ESMValGroup/ESMValTool/issues>`_
repository and on the mailing-list_.

Core development team
---------------------

-  Deutsches Zentrum für Luft- und Raumfahrt (DLR), Institut für Physik
   der Atmosphäre, Germany (Co-PI)

   - ESMValTool Core Co-PI and Developer: contact for requests to use the ESMValTool and for collaboration with the development team, access to the PRIVATE GitHub repository.

-  Met Office, United Kingdom (Co-PI)
-  Alfred Wegener institute (AWI) Bremerhaven, Germany
-  Barcelona Supercomputing Center (BSC), Spain
-  Netherlands eScience Center (NLeSC), The Netherlands
-  Ludwig Maximilian University of Munich, Germany
-  Plymouth Marine Laboratory (PML), United Kingdom
-  Swedish Meteorological and Hydrological Institute (SMHI), Sweden
-  University of Bremen, Germany
-  University of Reading, United Kingdom

Recipes and diagnostics
-----------------------

Contacts for specific diagnostic sets are the respective authors, as
listed in the corresponding :ref:`recipe and diagnostic documentation<recipes>`
and in the source code.


License
=======

The ESMValTool is released under the Apache License, version 2.0.
Citation of the ESMValTool paper ("Software Documentation Paper") is
kindly requested upon use, alongside with the software DOI for
ESMValTool
(`doi:10.5281/zenodo.3401363 <https://doi.org/10.5281/zenodo.3401363>`__)
and ESMValCore
(`doi:10.5281/zenodo.3387139 <https://doi.org/10.5281/zenodo.3387139>`__)
and version number:

-  Righi, M., Andela, B., Eyring, V., Lauer, A., Predoi, V., Schlund,
   M., Vegas-Regidor, J., Bock, L., Brötz, B., de Mora, L., Diblen, F.,
   Dreyer, L., Drost, N., Earnshaw, P., Hassler, B., Koldunov, N.,
   Little, B., Loosveldt Tomas, S., and Zimmermann, K.: Earth System
   Model Evaluation Tool (ESMValTool) v2.0 – technical overview, Geosci.
   Model Dev., 13, 1179–1199, https://doi.org/10.5194/gmd-13-1179-2020,
   2020.

Besides the above citation, users are kindly asked to register any
journal articles (or other scientific documents) that use the software
at the ESMValTool webpage (http://www.esmvaltool.org/). Citing the
Software Documentation Paper and registering your paper(s) will serve to
document the scientific impact of the Software, which is of vital
importance for securing future funding. You should consider this an
obligation if you have taken advantage of the ESMValTool, which
represents the end product of considerable effort by the development
team.

What ESMValTool can do for you
==============================

The ESMValTool applies a great variety of standard diagnostics and
metrics, and produces a collection of netCDF and graphical files
(plots). Thus, the tool needs a certain amount of input from the user so
that it can:

-  establish the correct input and output parameters and the structured
   workflow;
-  acquire the correct data;
-  execute the workflow; and
-  output the desired collective data and media.

To facilitate these four steps, the user has control over the tool via
two main input files: the :ref:`user configuration file <config-user>`
and the :ref:`recipe <esmvalcore:recipe>`. The configuration file sets
user and site-specific parameters (like input and output paths, desired
output graphical formats, logging level, etc.), whereas the recipe file
sets data, preprocessing and diagnostic-specific parameters (data
parameters grouped in the datasets sections, preprocessing steps for
various preprocessors sections, variables' parameters and
diagnostic-specific instructions grouped in the diagnostics sections).
The configuration file may be used for a very large number of runs with
very minimal changes since most of the parameters it sets are
recyclable; the recipe file can be used for a large number of
applications, since it may include as many datasets, preprocessors and
diagnostics sections as the user deems useful.

Once the user configuration files and the recipe are at hand, the user
can start the tool. A schematic overview of the ESMValTool workflow is
depited in the figure below.

.. container::
   :name: figarch

   .. figure:: figures/schematic.png
      :alt: Schematic of the system architecture.
      :figclass: align-center

      Schematic of the system architecture.

For a generalized run scenario, the tool will perform the following
ordered procedures.

Data finding
------------

-  read the data requirements from the :ref:`datasets section
   <esmvalcore:Datasets>` of the recipe and assemble the data request to
   locate the data;
-  find the data using the specified root paths and DRS types in the
   configuration file (note the flexibility allowed by the
   :ref:`data finder
   <esmvalcore:findingdata>`);

Data selection
--------------

-  data selection is performed using the parameters specified in the
   :ref:`datasets section <esmvalcore:Datasets>` (including e.g. type of
   experiment, type of ensemble, time boundaries etc); data will be
   retrieved and selected for each variable that is specified in the
   :ref:`diagnostics <esmvalcore:Diagnostics>` section of the recipe;

Data fixing
-----------

-  the ESMValTool requires data to be in CMOR format; since errors in
   the data formatting are not uncommon, the ESMValTool performs
   :ref:`checks against the
   CMOR library and fixes small irregularities <esmvalcore:CMOR check and
   dataset-specific fixes>` (note that the degree of leniency is not
   very high).

Variable derivation
-------------------

-  :ref:`variable derivation <esmvalcore:Variable derivation>` (in the
   case of non CMOR-standard variables, most likely associated with
   observational datasets) is performed automatically before running the
   preprocessor;
-  if the variable definitions are already in the database then the user
   will just have to specify the variableto be derived in the
   :ref:`diagnostics
   <esmvalcore:Diagnostics>` section (as any other standard variable,
   just setting ``derive: true``).

Run the preprocessor
--------------------

-  if any :ref:`preprocessor section <esmvalcore:preprocessor>` is
   specified in the recipe file, then data will be loaded in memory as
   iris cubes and passed through the preprocessing steps required by the
   user and specified in the preprocessor section, using the specific
   preprocessing step parameters provided by the user as keys (for the
   parameter name) and values (for the parameter value); the
   preprocessing order is very imprtant since a number of steps depend
   on prior execution of other steps (e.g. :ref:`multimodel
   statistics <esmvalcore:Multi-model statistics>` can not be computed
   unless all models are on a common grid, hence a prior
   :ref:`regridding
   <esmvalcore:Horizontal regridding>` on a common grid is necessary);
   the preprocessor steps order can be set by the user as custom or the
   default order can be used;
-  once preprocessing has finished, the tool writes the data output to
   disk as netCDF files so that the diagnostics can pick it up and use
   it; the user will also be provided with a metadata file containing a
   summary of the preprocessing and pointers to its output. Note that
   writing data to disk between the preprocessing and the diagnostic
   phase is required to ensure multi-language support for the latter.

Run the diagnostics
-------------------

-  the last and most important phase can now be run: using output files
   from the preprocessor, the diagnostic scripts are executed using the
   provided diagnostics parameters.
.. _inputdata:

********************
Obtaining input data
********************

ESMValTool supports input data from climate models participating in
`CMIP6 <https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip6>`__,
`CMIP5 <https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip5>`__,
`CMIP3 <https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip3>`__, and
`CORDEX <https://cordex.org/>`__
as well as observations, reanalysis, and any other data, provided that it
adheres to the
`CF conventions <https://cfconventions.org/>`__
and the data is described in a
`CMOR table <http://pcmdi.github.io/software/cmorTable/index.html>`__
as used in the various
`Climate Model Intercomparison Projects <http://pcmdi.github.io/mips/>`__.
This section provides some guidelines for unfamiliar users.

Because the amount of data required by ESMValTool is typically large, it is
recommended that you use the tool on a compute cluster where the data is
already available, for example because it is connected to an
`ESGF node <https://esgf.llnl.gov/index.html>`__.
Examples of such compute clusters are
`Mistral <https://www.dkrz.de/up/systems/mistral>`__
and
`Jasmin <https://www.jasmin.ac.uk/>`__,
but many more exist around the world.

If you do not have access to such a facility through your institute or the
project you are working on, you can request access by applying for the
`ENES Climate Analytics Service <https://portal.enes.org/data/data-metadata-service/climate-analytics-service>`__
or, if you need longer term access or more computational resources, the
`IS-ENES3 Trans-national Access call <https://portal.enes.org/data/data-metadata-service/analysis-platforms>`__.

If the options above are not available to you, ESMValTool also offers a feature
to make it easy to download CMIP6, CMIP5, CMIP3, CORDEX, and obs4MIPs from ESGF.

The chapter in the ESMValCore documentation on
:ref:`finding data <esmvalcore:findingdata>` explains how to
configure the ESMValTool so it can find locally available data and/or
download it from ESGF if it isn't available locally yet.

Models
======

If you do not have access to a compute cluster with the data already mounted,
the ESMValTool can automatically download any required data that is available on ESGF.
This is the recommended approach for first-time users to obtain some data for
running ESMValTool.
For example, run

.. code-block:: bash

    esmvaltool run --offline=False examples/recipe_python.yml

to run the default example recipe and automatically download the required data
to the directory ``~/climate_data``.
The data only needs to be downloaded once, every following run will re-use
previously downloaded data stored in this directory.
See :ref:`esmvalcore:config-esgf` for a more in depth explanation and the
available configuration options.

Alternatively, you can use an external tool called
`Synda <http://prodiguer.github.io/synda/index.html>`__
to maintain your own collection of ESGF data.

Observations
============

Observational and reanalysis products in the standard CF/CMOR format used in CMIP and required by the ESMValTool are available via the obs4MIPs and ana4mips projects at the ESGF (e.g., https://esgf-data.dkrz.de/projects/esgf-dkrz/). Their use is strongly recommended, when possible.

Other datasets not available in these archives can be obtained by the user from the respective sources and reformatted to the CF/CMOR standard. ESMValTool currently support two ways to perform this reformatting (aka 'CMORization'). The first is to use a CMORizer script to generate a local pool of reformatted data that can readily be used by the ESMValTool. The second way is to implement specific 'fixes' for your dataset. In that case, the reformatting is performed 'on the fly' during the execution of an ESMValTool recipe (note that one of the first preprocessor tasks is 'CMOR checks and fixes'). Below, both methods are explained in more detail.

Using a CMORizer script
-----------------------

ESMValTool comes with a set of CMORizers readily available.
The CMORizers are dataset-specific scripts that can be run once to generate
a local pool of CMOR-compliant data. The necessary information to download
and process the data is provided in the header of each CMORizing script.
These scripts also serve as template to create new CMORizers for datasets not
yet included.
Note that datasets CMORized for ESMValTool v1 may not be working with v2, due
to the much stronger constraints on metadata set by the iris library.

To CMORize one or more datasets, run:

.. code-block:: bash

    cmorize_obs -c [CONFIG_FILE] -o [DATASET_LIST]

The path to the raw data to be CMORized must be specified in the
:ref:`user configuration file<config-user>` as RAWOBS.
Within this path, the data are expected to be organized in subdirectories
corresponding to the data tier: Tier2 for freely-available datasets (other
than obs4MIPs and ana4mips) and Tier3 for restricted datasets (i.e., dataset
which requires a registration to be retrieved or provided upon request to
the respective contact or PI).
The CMORization follows the
`CMIP5 CMOR tables <https://github.com/PCMDI/cmip5-cmor-tables>`_ or
`CMIP6 CMOR tables <https://github.com/PCMDI/cmip6-cmor-tables>`_ for the
OBS and OBS6 projects respectively.
The resulting output is saved in the output_dir, again following the Tier
structure.
The output file names follow the definition given in
:ref:`config-developer file <esmvalcore:config-developer>` for the ``OBS``
project:

.. code-block::

    [project]_[dataset]_[type]_[version]_[mip]_[short_name]_YYYYMM_YYYYMM.nc

where ``project`` may be OBS (CMIP5 format) or OBS6 (CMIP6 format), ``type``
may be ``sat`` (satellite data), ``reanaly`` (reanalysis data),
``ground`` (ground observations), ``clim`` (derived climatologies),
``campaign`` (aircraft campaign).

At the moment, cmorize_obs supports Python and NCL scripts.

.. _cmorization_as_fix:

CMORization as a fix
--------------------
ESMValCore also provides support for some datasets in their native format.
In this case, the steps needed to reformat the data are executed as datasets
fixes during the execution of an ESMValTool recipe, as one of the first
preprocessor steps, see :ref:`fixing data <esmvalcore:fixing_data>`.
Compared to the workflow described above, this has the advantage that the user
does not need to store a duplicate (CMORized) copy of the data.
Instead, the CMORization is performed 'on the fly' when running a recipe.
The native6 project supports files named according to the format defined in
the :ref:`config-developer file <esmvalcore:config-developer>`.
Some of ERA5, ERA5-Land and MSWEP data are currently supported, see
:ref:`supported datasets <supported_datasets>`.

To use this functionality, users need to provide a path for the ``native6``
project data in the :ref:`user configuration file<config-user>`.
Then, in the recipe, they can refer to the native6 project.
For example:

.. code-block:: yaml

    datasets:
    - {dataset: ERA5, project: native6, type: reanaly, version: '1', tier: 3, start_year: 1990, end_year: 1990}

More examples can be found in the diagnostics ``ERA5_native6`` in the recipe
`examples/recipe_check_obs.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_check_obs.yml>`_.

.. _supported_datasets:

Supported datasets
------------------
A list of the datasets for which a CMORizers is available is provided in the following table.

.. tabularcolumns:: |p{3cm}|p{6cm}|p{3cm}|p{3cm}|

+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Dataset                      | Variables (MIP)                                                                                      | Tier | Script language |
+==============================+======================================================================================================+======+=================+
| APHRO-MA                     | pr, tas (day), pr, tas (Amon)                                                                        |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| AURA-TES                     | tro3 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| BerkelyEarth                 | tas, tasa (Amon), sftlf (fx)                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CALIPSO-GOCCP                | clcalipso (cfMon)                                                                                    |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-SATELLITE-ALBEDO         | bdalb (Lmon), bhalb (Lmon)                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-SATELLITE-LAI-FAPAR      | fapar (Lmon), lai (Lmon)                                                                             |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-SATELLITE-SOIL-MOISTURE  | sm (day), sm (Lmon)                                                                                  |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-UERRA                    | sm (E6hr)                                                                                            |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-XCH4                     | xch4 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-XCO2                     | xco2 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CERES-EBAF                   | rlut, rlutcs, rsut, rsutcs (Amon)                                                                    |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CERES-SYN1deg                | rlds, rldscs, rlus, rluscs, rlut, rlutcs, rsds, rsdscs, rsus, rsuscs, rsut, rsutcs (3hr)             |   3  | NCL             |
|                              | rlds, rldscs, rlus, rlut, rlutcs, rsds, rsdt, rsus, rsut, rsutcs (Amon)                              |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CLARA-AVHRR                  | clt, clivi, lwp (Amon)                                                                               |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CowtanWay                    | tasa (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CRU                          | tas, pr (Amon)                                                                                       |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CT2019                       | co2s (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Duveiller2018                | albDiffiTr13                                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| E-OBS                        | tas, tasmin, tasmax, pr, psl (day, Amon)                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Eppley-VGPM-MODIS            | intpp (Omon)                                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA5 [#note1]_               | clt, evspsbl, evspsblpot, mrro, pr, prsn, ps, psl, ptype, rls, rlds, rlns, rlus [#note2]_, rsds,     |   3  | n/a             |
|                              | rsns, rsus [#note2]_, rsdt, rss, uas, vas, tas, tasmax, tasmin, tdps, ts, tsn (E1hr/Amon), orog (fx) |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA5-Land [#note1]_          | pr                                                                                                   |   3  | n/a             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA-Interim                  | clivi, clt, clwvi, evspsbl, hur, hus, pr, prsn, prw, ps, psl, rlds, rsds, rsdt, ta, tas, tauu, tauv, |   3  | Python          |
|                              | ts, ua, uas, va, vas, wap, zg (Amon), ps, rsdt (CFday), clt, pr, prsn, psl, rsds, rss, ta, tas,      |      |                 |
|                              | tasmax, tasmin, uas, va, vas, zg (day), evspsbl, tdps, ts, tsn, rss, tdps (Eday), tsn (LImon), hfds, |      |                 |
|                              | tos (Omon), orog, sftlf (fx)                                                                         |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA-Interim-Land             | sm (Lmon)                                                                                            |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-AEROSOL               | abs550aer, od550aer, od550aerStderr, od550lt1aer, od870aer, od870aerStderr (aero)                    |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-CLOUD                 | clivi, clt, cltStderr, lwp, rlut, rlutcs, rsut, rsutcs, rsdt, rlus, rsus, rsuscs (Amon)              |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-FIRE                  | burntArea (Lmon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-LANDCOVER             | baresoilFrac, cropFrac, grassFrac, shrubFrac, treeFrac (Lmon)                                        |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-LST                   | ts (Amon)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OC                    | chl (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OZONE                 | toz, tozStderr, tro3prof, tro3profStderr (Amon)                                                      |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SEA-SURFACE-SALINITY  | sos (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SOILMOISTURE          | dos, dosStderr, sm, smStderr (Lmon)                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SST                   | ts, tsStderr (Amon)                                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-WATERVAPOUR           | prw (Amon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESRL                         | co2s (Amon)                                                                                          |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| FLUXCOM                      | gpp (Lmon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GCP2018                      | fgco2 (Omon), nbp (Lmon)                                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GCP2020                      | fgco2 (Omon), nbp (Lmon)                                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GHCN                         | pr (Amon)                                                                                            |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GHCN-CAMS                    | tas (Amon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GISTEMP                      | tasa (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GLODAP                       | dissic, ph, talk (Oyr)                                                                               |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GPCC                         | pr (Amon)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GRACE                        | lweGrace (Lmon)                                                                                      |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT3                     | tas, tasa (Amon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT4                     | tas, tasa (Amon), tasConf5, tasConf95                                                                |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT5                     | tas (Amon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadISST                      | sic (OImon), tos (Omon), ts (Amon)                                                                   |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HALOE                        | tro3, hus (Amon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HWSD                         | cSoil (Lmon), areacella (fx), sftlf (fx)                                                             |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ISCCP-FH                     | alb, prw, ps, rlds, rlus, rlut, rlutcs, rsds, rsdt, rsus, rsut, rsutcs, tas, ts (Amon)               |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| JMA-TRANSCOM                 | nbp (Lmon), fgco2 (Omon)                                                                             |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| LAI3g                        | lai (Lmon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| LandFlux-EVAL                | et, etStderr (Lmon)                                                                                  |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Landschuetzer2016            | dpco2, fgco2, spco2 (Omon)                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MAC-LWP                      | lwp, lwpStderr (Amon)                                                                                |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MERRA2                       | sm (Lmon)                                                                                            |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MLS-AURA                     | hur, hurStderr (day)                                                                                 |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MODIS                        | cliwi, clt, clwvi, iwpStderr, lwpStderr (Amon), od550aer (aero)                                      |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MSWEP [#note1]_              | pr                                                                                                   |   3  | n/a             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MTE                          | gpp, gppStderr (Lmon)                                                                                |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NCEP                         | hur, hus, pr, ta, tas, ua, va, wap, zg (Amon)                                                        |   2  | NCL             |
|                              | pr, rlut, ua, va (day)                                                                               |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NDP                          | cVeg (Lmon)                                                                                          |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NIWA-BS                      | toz, tozStderr (Amon)                                                                                |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NSIDC-0116-[nh|sh]           | usi, vsi (day)                                                                                       |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| OSI-450-[nh|sh]              | sic (OImon), sic (day)                                                                               |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PATMOS-x                     | clt (Amon)                                                                                           |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PERSIANN-CDR                 | pr (Amon), pr (day)                                                                                  |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PHC                          | thetao, so                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PIOMAS                       | sit (day)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| REGEN                        | pr (day, Amon)                                                                                       |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Scripps-CO2-KUM              | co2s (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| UWisc                        | clwvi, lwpStderr (Amon)                                                                              |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| WFDE5                        | tas, pr (Amon, day)                                                                                  |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| WOA                          | thetao, so, tos, sos (Omon)                                                                          |   2  | Python          |
|                              | no3, o2, po4, si (Oyr)                                                                               |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+

.. [#note1] CMORization is built into ESMValTool through the native6 project, so there is no separate CMORizer script.

.. [#note2] Derived on the fly from down & net radiation.
.. _changelog:

Changelog
=========

.. _changelog-v2-4-0:

v2.4.0
------

Highlights
~~~~~~~~~~

- ESMValTool is moving from Conda to Mamba as the preferred installation method. This will speed up the
  installation and comes with some improvements behind the scenes.
  Read more about it at :ref:`Move to Mamba<move-to-mamba>` and in :ref:`the installation guide<install>`.

Please also note the highlights from the corresponding ESMValCore release :ref:`here<esmvalcore:changelog-v2-4-0>`.
Thanks to that ESMValTool has gained the following features:

- Download any missing data that is available on the ESGF automatically.
- Resume previous runs, reusing expensive pre-processing results.


This release includes

Bug fixes
~~~~~~~~~

-  Fixed `recipe_meehl20sciadv.yml` for ESMValCore 2.3 (`#2253 <https://github.com/ESMValGroup/ESMValTool/pull/2253>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix provenance of NCL figures created using the log_provenance function (`#2279 <https://github.com/ESMValGroup/ESMValTool/pull/2279>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix bug in ClimWIP brunner19 recipe when plotting (`#2226 <https://github.com/ESMValGroup/ESMValTool/pull/2226>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Pin docutils <0.17 to fix sphinx build with rtd theme (`#2312 <https://github.com/ESMValGroup/ESMValTool/pull/2312>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix example recipes (`#2338 <https://github.com/ESMValGroup/ESMValTool/pull/2338>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Do not add bounds to plev (plev19) in era interim cmorizer (`#2328 <https://github.com/ESMValGroup/ESMValTool/pull/2328>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix problem with pip 21.3 that prevents installation from source (`#2344 <https://github.com/ESMValGroup/ESMValTool/pull/2344>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add title to recipe embedded in test_diagnostic_run.py (`#2353 <https://github.com/ESMValGroup/ESMValTool/pull/2353>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix capitalization of obs4MIPs (`#2368 <https://github.com/ESMValGroup/ESMValTool/pull/2368>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Specify that areacella is needed for area statistics in the Python example recipe (`#2371 <https://github.com/ESMValGroup/ESMValTool/pull/2371>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Enabling variable `obs550lt1aer` in recipes (`#2388 <https://github.com/ESMValGroup/ESMValTool/pull/2388>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Update a diagnostic to new Iris version (`#2390 <https://github.com/ESMValGroup/ESMValTool/pull/2390>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fixed bug in provenance tracking of ecs_scatter.ncl (`#2391 <https://github.com/ESMValGroup/ESMValTool/pull/2391>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix provenance issue in pv_capacity_factor.R (`#2392 <https://github.com/ESMValGroup/ESMValTool/pull/2392>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Remove obsolete write_plots option from R diagnostics (`#2395 <https://github.com/ESMValGroup/ESMValTool/pull/2395>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix arctic ocean diagnostic (`#2397 <https://github.com/ESMValGroup/ESMValTool/pull/2397>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix sea ice drift recipe and script (`#2404 <https://github.com/ESMValGroup/ESMValTool/pull/2404>`__) `sloosvel <https://github.com/sloosvel>`__
-  Adapt diagnostic script to new version of iris (`#2403 <https://github.com/ESMValGroup/ESMValTool/pull/2403>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix ocean multimap (`#2406 <https://github.com/ESMValGroup/ESMValTool/pull/2406>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix diagnostic that uses `xarray`: `dtype` correctly set and harmonize `xarray` and `matplotlib` (`#2409 <https://github.com/ESMValGroup/ESMValTool/pull/2409>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Deactivate provenance logging for plots in thermodyn toolbox (`#2414 <https://github.com/ESMValGroup/ESMValTool/pull/2414>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Deprecations
~~~~~~~~~~~~

-  Removed write_plots and write_netcdf from some NCL diagnostics (`#2293 <https://github.com/ESMValGroup/ESMValTool/pull/2293>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fixed provenance logging of all python diagnostics by removing 'plot_file' entry (`#2296 <https://github.com/ESMValGroup/ESMValTool/pull/2296>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Do not deprecate classes Variable, Variables and Datasets on a specific version (`#2286 <https://github.com/ESMValGroup/ESMValTool/pull/2286>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Remove obsolete write_netcdf option from ncl diagnostic scripts (`#2387 <https://github.com/ESMValGroup/ESMValTool/pull/2387>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Remove write plots from ocean diagnostics (`#2393 <https://github.com/ESMValGroup/ESMValTool/pull/2393>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  More removals of instances of `write_plots` from Python diagnostics (appears to be the final removal from Py diags) (`#2394 <https://github.com/ESMValGroup/ESMValTool/pull/2394>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Documentation
~~~~~~~~~~~~~

-  List Manuel Schlund as release manager for v2.5 (`#2268 <https://github.com/ESMValGroup/ESMValTool/pull/2268>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  GlobWat fix download links and gdal command (`#2334 <https://github.com/ESMValGroup/ESMValTool/pull/2334>`__) `Banafsheh Abdollahi <https://github.com/babdollahi>`__
-  Add titles to recipes authored by `predoi_valeriu` (`#2333 <https://github.com/ESMValGroup/ESMValTool/pull/2333>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Added titles to recipes maintained by lauer_axel (`#2332 <https://github.com/ESMValGroup/ESMValTool/pull/2332>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Update the documentation of the GRACE CMORizer (`#2349 <https://github.com/ESMValGroup/ESMValTool/pull/2349>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Add titles in BSC recipes (`#2351 <https://github.com/ESMValGroup/ESMValTool/pull/2351>`__) `sloosvel <https://github.com/sloosvel>`__
-  Update esmvalcore dependency to 2.4.0rc1 (`#2348 <https://github.com/ESMValGroup/ESMValTool/pull/2348>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add titles to recipes maintained by Peter Kalverla (`#2356 <https://github.com/ESMValGroup/ESMValTool/pull/2356>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Adding titles to the recipes with maintainer hb326 (`#2358 <https://github.com/ESMValGroup/ESMValTool/pull/2358>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Add title for zmnam as for #2354 (`#2363 <https://github.com/ESMValGroup/ESMValTool/pull/2363>`__) `fserva <https://github.com/fserva>`__
-  Added recipe titles the the ocean recipes.  (`#2364 <https://github.com/ESMValGroup/ESMValTool/pull/2364>`__) `Lee de Mora <https://github.com/ledm>`__
-  Update recipe_thermodyn_diagtool.yml - add title (`#2365 <https://github.com/ESMValGroup/ESMValTool/pull/2365>`__) `ValerioLembo <https://github.com/ValerioLembo>`__
-  Fix provenance of figures of several R diagnostics (`#2300 <https://github.com/ESMValGroup/ESMValTool/pull/2300>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Adding titles to Mattia's recipes (`#2367 <https://github.com/ESMValGroup/ESMValTool/pull/2367>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Adding titles to wenzel recipes (`#2366 <https://github.com/ESMValGroup/ESMValTool/pull/2366>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Fix formatting of some recipe titles merged from PR 2364 (`#2372 <https://github.com/ESMValGroup/ESMValTool/pull/2372>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Adding titles to Bjoern's recipes (`#2369 <https://github.com/ESMValGroup/ESMValTool/pull/2369>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Add titles to ocean recipes (maintainer Lovato) (`#2375 <https://github.com/ESMValGroup/ESMValTool/pull/2375>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Add titles for three c3s-magic recipes (`#2378 <https://github.com/ESMValGroup/ESMValTool/pull/2378>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add title for recipe maintained by Ruth Lorenz (`#2379 <https://github.com/ESMValGroup/ESMValTool/pull/2379>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix toymodel recipe (`#2381 <https://github.com/ESMValGroup/ESMValTool/pull/2381>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Added titles for recipes of maintainer `schlund_manuel` (`#2377 <https://github.com/ESMValGroup/ESMValTool/pull/2377>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Write_plots and titles for deangelis15nat, li17natcc, martin18grl, pv_capacity_factor (`#2382 <https://github.com/ESMValGroup/ESMValTool/pull/2382>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Add titles for some recipes (`#2383 <https://github.com/ESMValGroup/ESMValTool/pull/2383>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Adding titles for recipes by von Hardenberg and Arnone (`#2384 <https://github.com/ESMValGroup/ESMValTool/pull/2384>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Last two missing titles (`#2386 <https://github.com/ESMValGroup/ESMValTool/pull/2386>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update documentation on downloading data (`#2370 <https://github.com/ESMValGroup/ESMValTool/pull/2370>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix installation instructions for Julia (`#2335 <https://github.com/ESMValGroup/ESMValTool/pull/2335>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix provenance of Julia example diagnostic (`#2289 <https://github.com/ESMValGroup/ESMValTool/pull/2289>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Added notes on use of mamba in the installation documentation chapter (`#2236 <https://github.com/ESMValGroup/ESMValTool/pull/2236>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update version number for 2.4.0 release (`#2410 <https://github.com/ESMValGroup/ESMValTool/pull/2410>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update release schedule for 2.4.0 (`#2412 <https://github.com/ESMValGroup/ESMValTool/pull/2412>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update changelog for 2.4.0 release (`#2411 <https://github.com/ESMValGroup/ESMValTool/pull/2411>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Diagnostics
~~~~~~~~~~~

-  Add all available CMIP5 and CMIP6 models to recipe_impact.yml (`#2251 <https://github.com/ESMValGroup/ESMValTool/pull/2251>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add Fig. 6, 7 and 9 of Bock20jgr (`#2252 <https://github.com/ESMValGroup/ESMValTool/pull/2252>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Generalize `recipe_validation*` diagnostic to work with identical control and experiment dataset names (`#2284 <https://github.com/ESMValGroup/ESMValTool/pull/2284>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add missing preprocessor to recipe_gier2020bg and adapt to available data (`#2399 <https://github.com/ESMValGroup/ESMValTool/pull/2399>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Removed custom version of `AtmosphereSigmaFactory` in diagnostics (`#2405 <https://github.com/ESMValGroup/ESMValTool/pull/2405>`__) `Manuel Schlund <https://github.com/schlunma>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Replace recipe_era5.yml with recipe_daily_era5.yml (`#2182 <https://github.com/ESMValGroup/ESMValTool/pull/2182>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update WOA cmorizer for WOA18 and WOA13v2 (`#1812 <https://github.com/ESMValGroup/ESMValTool/pull/1812>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  GLODAP v2.2016 ocean data cmorizer (`#2185 <https://github.com/ESMValGroup/ESMValTool/pull/2185>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Updated GCP CMORizer (`#2295 <https://github.com/ESMValGroup/ESMValTool/pull/2295>`__) `Manuel Schlund <https://github.com/schlunma>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Add a cylc suite to run all recipes (`#2219 <https://github.com/ESMValGroup/ESMValTool/pull/2219>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Retire test with Python 3.6 from full development Github Actions test (`#2229 <https://github.com/ESMValGroup/ESMValTool/pull/2229>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Remove Python 3.6 tests from GitHub Actions (`#2264 <https://github.com/ESMValGroup/ESMValTool/pull/2264>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Unpin upper bound for iris (previously was at <3.0.4) (`#2266 <https://github.com/ESMValGroup/ESMValTool/pull/2266>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin latest esmvalcore to allow use of the bugfix release 2.3.1 always (`#2269 <https://github.com/ESMValGroup/ESMValTool/pull/2269>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add apt update so Julia gets found and installed by Docker (`#2290 <https://github.com/ESMValGroup/ESMValTool/pull/2290>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Use mamba for environment update and creation in the Docker container build on DockerHub (`#2297 <https://github.com/ESMValGroup/ESMValTool/pull/2297>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Docker container experimental - run a full env solve with mamba instead of a conda update (`#2306 <https://github.com/ESMValGroup/ESMValTool/pull/2306>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Full use of mamba in Github Actions source install test and use generic Python 3.7 (removing the very specific 3.7.10) (`#2287 <https://github.com/ESMValGroup/ESMValTool/pull/2287>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Replace use of conda with mamba for conda_install test on Circle CI (`#2237 <https://github.com/ESMValGroup/ESMValTool/pull/2237>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update circleci configuration (`#2357 <https://github.com/ESMValGroup/ESMValTool/pull/2357>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Installation
~~~~~~~~~~~~

-  Remove `mpich` from conda dependencies list (`#2343 <https://github.com/ESMValGroup/ESMValTool/pull/2343>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Add script for extracting a list of input files from the provenance (`#2278 <https://github.com/ESMValGroup/ESMValTool/pull/2278>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update github actions (`#2360 <https://github.com/ESMValGroup/ESMValTool/pull/2360>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Removed 'write_plots' from all NCL diagnostics (`#2331 <https://github.com/ESMValGroup/ESMValTool/pull/2331>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Update and modernize `config-user-example.yml` (`#2374 <https://github.com/ESMValGroup/ESMValTool/pull/2374>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__


.. _changelog-v2-3-0:

v2.3.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Indent block to pick up and raise exception if cmorizer data not found (TierX dir is not there) (`#1877 <https://github.com/ESMValGroup/ESMValTool/pull/1877>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Skip recipe filler tests until we have a new release since GA tests are failing (`#2089 <https://github.com/ESMValGroup/ESMValTool/pull/2089>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fixed broken link to contributions in README (`#2102 <https://github.com/ESMValGroup/ESMValTool/pull/2102>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix recipe filler for the case the variable doesn't contain short_name (`#2104 <https://github.com/ESMValGroup/ESMValTool/pull/2104>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add fix for iris longitude bug to ClimWIP (`#2107 <https://github.com/ESMValGroup/ESMValTool/pull/2107>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Update for outdated link to reference Déandreis et al. (2014). (`#2076 <https://github.com/ESMValGroup/ESMValTool/pull/2076>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fixed recipes for ESMValCore 2.3.0 (`#2203 <https://github.com/ESMValGroup/ESMValTool/pull/2203>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix the WFDE5 cmorizer (`#2211 <https://github.com/ESMValGroup/ESMValTool/pull/2211>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fix broken CMORizer log message if no Tier directory exists (`#2207 <https://github.com/ESMValGroup/ESMValTool/pull/2207>`__) `jmrgonza <https://github.com/jmrgonza>`__
-  Fix bug in ClimWIP basic test recipe when plotting (`#2225 <https://github.com/ESMValGroup/ESMValTool/pull/2225>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Fix bug in ClimWIP advanced test recipe when plotting (`#2227 <https://github.com/ESMValGroup/ESMValTool/pull/2227>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Adjust time range for the `WDFE5` dataset in the `recipe_check_obs.yml` (`#2232 <https://github.com/ESMValGroup/ESMValTool/pull/2232>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fix plot and provenance of recipe_consecdrydays (`#2244 <https://github.com/ESMValGroup/ESMValTool/pull/2244>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

-  Improving the README.md file with a more appealing look and bit more info (`#2065 <https://github.com/ESMValGroup/ESMValTool/pull/2065>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update plot title martin18grl (`#2080 <https://github.com/ESMValGroup/ESMValTool/pull/2080>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Update contribution guidelines (`#2031 <https://github.com/ESMValGroup/ESMValTool/pull/2031>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update links in pull request template to point to latest documentation (`#2083 <https://github.com/ESMValGroup/ESMValTool/pull/2083>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update release schedule (`#2081 <https://github.com/ESMValGroup/ESMValTool/pull/2081>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Updates to contribution guidelines (`#2092 <https://github.com/ESMValGroup/ESMValTool/pull/2092>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update documentation for ERA5 with new variables (`#2111 <https://github.com/ESMValGroup/ESMValTool/pull/2111>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Add OSX installation instructions to docs (`#2115 <https://github.com/ESMValGroup/ESMValTool/pull/2115>`__) `Barbara Vreede <https://github.com/bvreede>`__
-  Instructions to use pre-installed versions on HPC clusters (`#2197 <https://github.com/ESMValGroup/ESMValTool/pull/2197>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Add functional Autoassess diagnostics: land surface metrics: permafrost, soil moisture, surface radiation (`#2170 <https://github.com/ESMValGroup/ESMValTool/pull/2170>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add citation info in `recipe_eady_growth_rate.yml` (`#2188 <https://github.com/ESMValGroup/ESMValTool/pull/2188>`__) `sloosvel <https://github.com/sloosvel>`__
-  Update version number to 2.3.0 (`#2213 <https://github.com/ESMValGroup/ESMValTool/pull/2213>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update release schedule for 2.3.0 (`#2247 <https://github.com/ESMValGroup/ESMValTool/pull/2247>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Changelog update to v2.3.0 (`#2214 <https://github.com/ESMValGroup/ESMValTool/pull/2214>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Diagnostics
~~~~~~~~~~~

-  Added figures 8 and 10 to recipe_bock20jgr.yml (`#2074 <https://github.com/ESMValGroup/ESMValTool/pull/2074>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add hydrological forcing comparison recipe (`#2013 <https://github.com/ESMValGroup/ESMValTool/pull/2013>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Added recipe for Meehl et al., Sci. Adv. (2020) (`#2094 <https://github.com/ESMValGroup/ESMValTool/pull/2094>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add GlobWat recipe and diagnostic  (`#1808 <https://github.com/ESMValGroup/ESMValTool/pull/1808>`__) `Banafsheh Abdollahi <https://github.com/babdollahi>`__
-  Add ClimWIP recipe to reproduce Brunner et al. 2019 (`#2109 <https://github.com/ESMValGroup/ESMValTool/pull/2109>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Update Climwip recipe to reproduce brunner2020esd (`#1859 <https://github.com/ESMValGroup/ESMValTool/pull/1859>`__) `Ruth Lorenz <https://github.com/ruthlorenz>`__
-  Update recipe_thermodyn_diagtool.yml: code improvements and more user options (`#1391 <https://github.com/ESMValGroup/ESMValTool/pull/1391>`__) `ValerioLembo <https://github.com/ValerioLembo>`__
-  Remove model AWI-CM-1-1-MR from recipe_impact.yml (`#2238 <https://github.com/ESMValGroup/ESMValTool/pull/2238>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  PV capacity factor for ESMValTool GMD paper  (`#2153 <https://github.com/ESMValGroup/ESMValTool/pull/2153>`__) `katjaweigel <https://github.com/katjaweigel>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorize wfde5 (`#1991 <https://github.com/ESMValGroup/ESMValTool/pull/1991>`__) `mwjury <https://github.com/mwjury>`__
-  Make cmorizer utils funcs public in utilities.py and add some numpy style docstrings (`#2206 <https://github.com/ESMValGroup/ESMValTool/pull/2206>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  CMORizer for CLARA-AVHRR cloud data (`#2101 <https://github.com/ESMValGroup/ESMValTool/pull/2101>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Update of ESACCI-CLOUD CMORizer (`#2144 <https://github.com/ESMValGroup/ESMValTool/pull/2144>`__) `Axel Lauer <https://github.com/axel-lauer>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Force latest Python in empty environment in conda install CI test (`#2069 <https://github.com/ESMValGroup/ESMValTool/pull/2069>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Removed imports from private sklearn modules and improved test coverage of custom_sklearn.py (`#2078 <https://github.com/ESMValGroup/ESMValTool/pull/2078>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Move private _(global)_stock_cube from esmvacore.preprocessor._regrid to cmorizer (`#2087 <https://github.com/ESMValGroup/ESMValTool/pull/2087>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Try mamba install esmvaltool (`#2125 <https://github.com/ESMValGroup/ESMValTool/pull/2125>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Reinstate OSX Github Action tests (`#2110 <https://github.com/ESMValGroup/ESMValTool/pull/2110>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin mpich to avoid default install of 3.4.1 and 3.4.2 with external_0 builds (`#2220 <https://github.com/ESMValGroup/ESMValTool/pull/2220>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Include test sources in distribution (`#2234 <https://github.com/ESMValGroup/ESMValTool/pull/2234>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Pin `iris<3.0.4` to ensure we still (sort of) support Python 3.6 (`#2246 <https://github.com/ESMValGroup/ESMValTool/pull/2246>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Installation
~~~~~~~~~~~~

-  Fix conda build by skipping documentation test (`#2058 <https://github.com/ESMValGroup/ESMValTool/pull/2058>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Update pin on esmvalcore pick up esmvalcore=2.3.0 (`#2200 <https://github.com/ESMValGroup/ESMValTool/pull/2200>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin Python to 3.9 for development installation (`#2208 <https://github.com/ESMValGroup/ESMValTool/pull/2208>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Improvements
~~~~~~~~~~~~

-  Add EUCP and IS-ENES3 projects to config-references (`#2066 <https://github.com/ESMValGroup/ESMValTool/pull/2066>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Fix flake8 tests on CircleCI (`#2070 <https://github.com/ESMValGroup/ESMValTool/pull/2070>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Added recipe filler. (`#1707 <https://github.com/ESMValGroup/ESMValTool/pull/1707>`__) `Lee de Mora <https://github.com/ledm>`__
-  Update use of fx vars to new syntax  (`#2145 <https://github.com/ESMValGroup/ESMValTool/pull/2145>`__) `sloosvel <https://github.com/sloosvel>`__
-  Add recipe for climate impact research (`#2072 <https://github.com/ESMValGroup/ESMValTool/pull/2072>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Update references "master" to "main" (`#2172 <https://github.com/ESMValGroup/ESMValTool/pull/2172>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Force git to ignore VSCode workspace files (`#2186 <https://github.com/ESMValGroup/ESMValTool/pull/2186>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Update to new ESMValTool logo (`#2168 <https://github.com/ESMValGroup/ESMValTool/pull/2168>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Python cmorizers for CDR1 and CDR2 ESACCI H2O (TCWV=prw) data. (`#2152 <https://github.com/ESMValGroup/ESMValTool/pull/2152>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Remove obsolete conda package (closes #2100) (`#2103 <https://github.com/ESMValGroup/ESMValTool/pull/2103>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

.. _changelog-v2-2-0:

v2.2.0
------

Highlights
~~~~~~~~~~

ESMValTool is now using the recently released `Iris 3 <https://scitools-iris.readthedocs.io/en/latest/whatsnew/3.0.html>`__.
We acknowledge that this change may impact your work, as Iris 3 introduces
several changes that are not backward-compatible, but we think that moving forward is the best
decision for the tool in the long term.


This release includes

Bug fixes
~~~~~~~~~

-  Bugfix: time weights in time_operations (`#1956 <https://github.com/ESMValGroup/ESMValTool/pull/1956>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Fix issues with bibtex references (`#1955 <https://github.com/ESMValGroup/ESMValTool/pull/1955>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Fix ImportError for `configure_logging` (`#1976 <https://github.com/ESMValGroup/ESMValTool/pull/1976>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add required functional parameters for extract time in recipe_er5.yml (`#1978 <https://github.com/ESMValGroup/ESMValTool/pull/1978>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Revert "Fix ImportError for `configure_logging`" (`#1992 <https://github.com/ESMValGroup/ESMValTool/pull/1992>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix import of esmvalcore _logging module in cmorize_obs.py (`#2020 <https://github.com/ESMValGroup/ESMValTool/pull/2020>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix logging import in cmorize_obs again since last merge was nulled by pre-commit hooks (`#2022 <https://github.com/ESMValGroup/ESMValTool/pull/2022>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Refactor the functions in derive_evspsblpot due to new iris (`#2023 <https://github.com/ESMValGroup/ESMValTool/pull/2023>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Avoid importing private ESMValCore functions in CMORizer (`#2027 <https://github.com/ESMValGroup/ESMValTool/pull/2027>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix extract_seasons in validation recipe  (`#2054 <https://github.com/ESMValGroup/ESMValTool/pull/2054>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Deprecations
~~~~~~~~~~~~

-  Deprecate classes Variable, Variables and Datasets (`#1944 <https://github.com/ESMValGroup/ESMValTool/pull/1944>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Python 3.9: remove pynio as dependency and replace with rasterio and pin Matplotlib>3.3.1 and pin cartopy>=0.18 (`#1997 <https://github.com/ESMValGroup/ESMValTool/pull/1997>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Removed write_plots and write_netcdf in some python diagnostics (`#2036 <https://github.com/ESMValGroup/ESMValTool/pull/2036>`__) `Manuel Schlund <https://github.com/schlunma>`__

Documentation
~~~~~~~~~~~~~

-  Update instructions on making a release (`#1867 <https://github.com/ESMValGroup/ESMValTool/pull/1867>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update review.rst (`#1917 <https://github.com/ESMValGroup/ESMValTool/pull/1917>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add guidance on how to review a pull request (`#1872 <https://github.com/ESMValGroup/ESMValTool/pull/1872>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Adding tutorial links to documentation (`#1927 <https://github.com/ESMValGroup/ESMValTool/pull/1927>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Added bibtex file for schlund20jgr (`#1928 <https://github.com/ESMValGroup/ESMValTool/pull/1928>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Documentation contact added the actual email for the mailing list (`#1938 <https://github.com/ESMValGroup/ESMValTool/pull/1938>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Make CircleCI badge specific to main branch (`#1831 <https://github.com/ESMValGroup/ESMValTool/pull/1831>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Documentation on how to move code from a private repository to a public repository (`#1920 <https://github.com/ESMValGroup/ESMValTool/pull/1920>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Refine pull request review guidelines (`#1924 <https://github.com/ESMValGroup/ESMValTool/pull/1924>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Update release schedule (`#1948 <https://github.com/ESMValGroup/ESMValTool/pull/1948>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Improve contact info and move to more prominent location (`#1950 <https://github.com/ESMValGroup/ESMValTool/pull/1950>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add some maintainers to some recipes that are missing them (`#1970 <https://github.com/ESMValGroup/ESMValTool/pull/1970>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update core team info (`#1973 <https://github.com/ESMValGroup/ESMValTool/pull/1973>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Combine installation from source instructions and add common issues (`#1971 <https://github.com/ESMValGroup/ESMValTool/pull/1971>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update iris documentation URL for sphinx (`#2003 <https://github.com/ESMValGroup/ESMValTool/pull/2003>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix iris documentation link(s) with new iris3 location on readthedocs (`#2012 <https://github.com/ESMValGroup/ESMValTool/pull/2012>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Document how to run tests for installation verification  (`#1847 <https://github.com/ESMValGroup/ESMValTool/pull/1847>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  List Remi Kazeroni as a code owner and sole merger of CMORizers (`#2017 <https://github.com/ESMValGroup/ESMValTool/pull/2017>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Install documentation: mention that we build conda package with python>=3.7 (`#2030 <https://github.com/ESMValGroup/ESMValTool/pull/2030>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Recipe and documentation update for ERA5-Land. (`#1906 <https://github.com/ESMValGroup/ESMValTool/pull/1906>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Update changelog and changelog tool for v2.2.0 (`#2043 <https://github.com/ESMValGroup/ESMValTool/pull/2043>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Final update to the changelog for v2.2.0 (`#2056 <https://github.com/ESMValGroup/ESMValTool/pull/2056>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Diagnostics
~~~~~~~~~~~

-  Add mapplot diagnostic to ClimWIP (`#1864 <https://github.com/ESMValGroup/ESMValTool/pull/1864>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Add the option to weight variable groups in ClimWIP (`#1856 <https://github.com/ESMValGroup/ESMValTool/pull/1856>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Implementation of ensemble member recognition to the ClimWIP diagnostic (`#1852 <https://github.com/ESMValGroup/ESMValTool/pull/1852>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Restructure ClimWIP (`#1919 <https://github.com/ESMValGroup/ESMValTool/pull/1919>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Diagnostic for recipe_eyring13jgr.yml Fig. 12 (`#1922 <https://github.com/ESMValGroup/ESMValTool/pull/1922>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Added changes in shared functions necessary for schlund20esd (`#1967 <https://github.com/ESMValGroup/ESMValTool/pull/1967>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Adding recipe and diagnostics for Gier et al 2020 (`#1914 <https://github.com/ESMValGroup/ESMValTool/pull/1914>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Added recipe, diagnostics and documentation for Schlund et al., ESD (2020) (`#2015 <https://github.com/ESMValGroup/ESMValTool/pull/2015>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add PRIMAVERA Eady Growth Rate diagnostic (`#1285 <https://github.com/ESMValGroup/ESMValTool/pull/1285>`__) `sloosvel <https://github.com/sloosvel>`__
-  Implement shape parameter calibration for ClimWIP (`#1905 <https://github.com/ESMValGroup/ESMValTool/pull/1905>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Extended ESRL cmorizer (`#1937 <https://github.com/ESMValGroup/ESMValTool/pull/1937>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Cmorizer for GRACE data (`#1694 <https://github.com/ESMValGroup/ESMValTool/pull/1694>`__) `bascrezee <https://github.com/bascrezee>`__
-  Cmorizer for latest ESACCI-SST data (`#1895 <https://github.com/ESMValGroup/ESMValTool/pull/1895>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix longitude in ESRL cmorizer (`#1988 <https://github.com/ESMValGroup/ESMValTool/pull/1988>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Selectively turn off fixing bounds for coordinates during cmorization with utilities.py (`#2014 <https://github.com/ESMValGroup/ESMValTool/pull/2014>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Cmorize hadcrut5 (`#1977 <https://github.com/ESMValGroup/ESMValTool/pull/1977>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize gpcc masking (`#1995 <https://github.com/ESMValGroup/ESMValTool/pull/1995>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize_utils_save_1mon_Amon (`#1990 <https://github.com/ESMValGroup/ESMValTool/pull/1990>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize gpcc fix (`#1982 <https://github.com/ESMValGroup/ESMValTool/pull/1982>`__) `mwjury <https://github.com/mwjury>`__
-  Fix flake8 raised by develop test in cmorize_obs_gpcc.py (`#2038 <https://github.com/ESMValGroup/ESMValTool/pull/2038>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Switched miniconda conda setup hooks for Github Actions workflows (`#1913 <https://github.com/ESMValGroup/ESMValTool/pull/1913>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix style issue (`#1929 <https://github.com/ESMValGroup/ESMValTool/pull/1929>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix mlr test with solution that works for CentOS too (`#1936 <https://github.com/ESMValGroup/ESMValTool/pull/1936>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Temporary deactivation Github Actions on OSX (`#1939 <https://github.com/ESMValGroup/ESMValTool/pull/1939>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix conda installation test on CircleCI (`#1952 <https://github.com/ESMValGroup/ESMValTool/pull/1952>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Github Actions: change time for cron job that installs from conda (`#1969 <https://github.com/ESMValGroup/ESMValTool/pull/1969>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  CI upload relevant artifacts for test job (`#1999 <https://github.com/ESMValGroup/ESMValTool/pull/1999>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Github Actions test that runs with the latest ESMValCore main (`#1989 <https://github.com/ESMValGroup/ESMValTool/pull/1989>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Introduce python 39 in Github Actions tests (`#2029 <https://github.com/ESMValGroup/ESMValTool/pull/2029>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Remove test for conda package installation on Python 3.6 (`#2033 <https://github.com/ESMValGroup/ESMValTool/pull/2033>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update codacy coverage reporter to fix coverage (`#2039 <https://github.com/ESMValGroup/ESMValTool/pull/2039>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Installation
~~~~~~~~~~~~

-  Simplify installation of R development dependencies (`#1930 <https://github.com/ESMValGroup/ESMValTool/pull/1930>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix docker build (`#1934 <https://github.com/ESMValGroup/ESMValTool/pull/1934>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Use new conda environment for installing ESMValTool in Docker containers (`#1993 <https://github.com/ESMValGroup/ESMValTool/pull/1993>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix conda build (`#2026 <https://github.com/ESMValGroup/ESMValTool/pull/2026>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Improvements
~~~~~~~~~~~~

-  Allow multiple references for a cmorizer script (`#1953 <https://github.com/ESMValGroup/ESMValTool/pull/1953>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Add GRACE to the recipe check_obs (`#1963 <https://github.com/ESMValGroup/ESMValTool/pull/1963>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Align ESMValTool to ESMValCore=2.2.0 (adopt iris3, fix environment for new Core release) (`#1874 <https://github.com/ESMValGroup/ESMValTool/pull/1874>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Make it possible to use write_plots and write_netcdf from recipe instead of config-user.yml (`#2018 <https://github.com/ESMValGroup/ESMValTool/pull/2018>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Revise lisflood and hype recipes (`#2035 <https://github.com/ESMValGroup/ESMValTool/pull/2035>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Set version to 2.2.0 (`#2042 <https://github.com/ESMValGroup/ESMValTool/pull/2042>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

.. _changelog-v2-1-1:

v2.1.1
------

This release includes

Improvements
~~~~~~~~~~~~

- Fix the conda build on CircleCI (`#1883 <https://github.com/ESMValGroup/ESMValTool/pull/1883>`__) `Bouwe Andela <https://github.com/bouweandela>`__
- Pin matplotlib to <3.3 and add compilers (`#1898 <https://github.com/ESMValGroup/ESMValTool/pull/1898>`__) `Bouwe Andela <https://github.com/bouweandela>`__
- Pin esmvaltool subpackages to the same version and build as the esmvaltool conda package (`#1899 <https://github.com/ESMValGroup/ESMValTool/pull/1899>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

- Release notes v2.1.1 (`#1932 <https://github.com/ESMValGroup/ESMValTool/pull/1932>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

.. _changelog-v2-1-0:

v2.1.0
------

This release includes

Diagnostics
~~~~~~~~~~~

-  Add extra steps to diagnostic to make output of hydrology/recipe_lisflood.yml usable by the LISFLOOD model (`#1737 <https://github.com/ESMValGroup/ESMValTool/pull/1737>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__
-  Recipe to reproduce the 2014 KNMI Climate Scenarios (kcs). (`#1667 <https://github.com/ESMValGroup/ESMValTool/pull/1667>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Implement the climwip weighting scheme in a recipe and diagnostic (`#1648 <https://github.com/ESMValGroup/ESMValTool/pull/1648>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__
-  Remove unreviewed autoassess recipes (`#1840 <https://github.com/ESMValGroup/ESMValTool/pull/1840>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Changes in shared scripts for Schlund et al., JGR: Biogeosciences, 2020 (`#1845 <https://github.com/ESMValGroup/ESMValTool/pull/1845>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Updated derivation test recipe (`#1790 <https://github.com/ESMValGroup/ESMValTool/pull/1790>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Support for multiple model occurrence in perf main (`#1649 <https://github.com/ESMValGroup/ESMValTool/pull/1649>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Add recipe and diagnostics for Schlund et al., JGR: Biogeosciences, 2020 (`#1860 <https://github.com/ESMValGroup/ESMValTool/pull/1860>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Adjust recipe_extract_shape.yml to recent changes in the example diagnostic.py (`#1880 <https://github.com/ESMValGroup/ESMValTool/pull/1880>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

-  Add pip installation instructions (`#1783 <https://github.com/ESMValGroup/ESMValTool/pull/1783>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add installation instruction for R and Julia dependencies tot pip install (`#1787 <https://github.com/ESMValGroup/ESMValTool/pull/1787>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Avoid autodocsumm 0.2.0 and update documentation build dependencies (`#1794 <https://github.com/ESMValGroup/ESMValTool/pull/1794>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more information on working on cluster attached to ESGF node (`#1821 <https://github.com/ESMValGroup/ESMValTool/pull/1821>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add release strategy to community documentation (`#1809 <https://github.com/ESMValGroup/ESMValTool/pull/1809>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update esmvaltool run command everywhere in documentation (`#1820 <https://github.com/ESMValGroup/ESMValTool/pull/1820>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more info on documenting a recipe (`#1795 <https://github.com/ESMValGroup/ESMValTool/pull/1795>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve the Python example diagnostic and documentation (`#1827 <https://github.com/ESMValGroup/ESMValTool/pull/1827>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve description of how to use draft_release_notes.py (`#1848 <https://github.com/ESMValGroup/ESMValTool/pull/1848>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update changelog for release 2.1 (`#1886 <https://github.com/ESMValGroup/ESMValTool/pull/1886>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Fix R installation in WSL (`#1789 <https://github.com/ESMValGroup/ESMValTool/pull/1789>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add pre-commit for linting/formatting (`#1796 <https://github.com/ESMValGroup/ESMValTool/pull/1796>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Speed up tests on CircleCI and use pytest to run them (`#1804 <https://github.com/ESMValGroup/ESMValTool/pull/1804>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Move pre-commit excludes to top-level and correct order of lintr and styler (`#1805 <https://github.com/ESMValGroup/ESMValTool/pull/1805>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Remove isort setup to fix formatting conflict with yapf (`#1815 <https://github.com/ESMValGroup/ESMValTool/pull/1815>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  GitHub Actions (`#1806 <https://github.com/ESMValGroup/ESMValTool/pull/1806>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix yapf-isort import formatting conflict (`#1822 <https://github.com/ESMValGroup/ESMValTool/pull/1822>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Replace vmprof with vprof as the default profiler (`#1829 <https://github.com/ESMValGroup/ESMValTool/pull/1829>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update ESMValCore v2.1.0 requirement (`#1839 <https://github.com/ESMValGroup/ESMValTool/pull/1839>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Pin iris to version 2 (`#1881 <https://github.com/ESMValGroup/ESMValTool/pull/1881>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Pin eccodes to not use eccodes=2.19.0 for cdo to work fine (`#1869 <https://github.com/ESMValGroup/ESMValTool/pull/1869>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Increase version to 2.1.0 and add release notes (`#1868 <https://github.com/ESMValGroup/ESMValTool/pull/1868>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Github Actions Build Packages and Deploy tests (conda and PyPi) (`#1858 <https://github.com/ESMValGroup/ESMValTool/pull/1858>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Added CMORizer for Scripps-CO2-KUM (`#1857 <https://github.com/ESMValGroup/ESMValTool/pull/1857>`__) `Manuel Schlund <https://github.com/schlunma>`__

.. _changelog-v2-0-0:

v2.0.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Fix pep8-naming errors and fix zmnam diagnostic (`#1702 <https://github.com/ESMValGroup/ESMValTool/pull/1702>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix keyword argument in cmorize_obs (`#1721 <https://github.com/ESMValGroup/ESMValTool/pull/1721>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fixed JMA-TRANSCOM CMORizer (`#1735 <https://github.com/ESMValGroup/ESMValTool/pull/1735>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix bug in extract_doi_value (`#1734 <https://github.com/ESMValGroup/ESMValTool/pull/1734>`__) `bascrezee <https://github.com/bascrezee>`__
-  Fix small errors in the arctic_ocean diagnostic (`#1722 <https://github.com/ESMValGroup/ESMValTool/pull/1722>`__) `Nikolay Koldunov <https://github.com/koldunovn>`__
-  Flatten ancestor lists for diag_spei.R and diag_spi.R. (`#1745 <https://github.com/ESMValGroup/ESMValTool/pull/1745>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fix for recipe_ocean_ice_extent.yml (`#1744 <https://github.com/ESMValGroup/ESMValTool/pull/1744>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix recipe_combined_indices.yml provenance (`#1746 <https://github.com/ESMValGroup/ESMValTool/pull/1746>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix provenance in recipe_multimodel_products (`#1747 <https://github.com/ESMValGroup/ESMValTool/pull/1747>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Exclude FGOALS-g2 due to ESMValCore issue #728 (`#1749 <https://github.com/ESMValGroup/ESMValTool/pull/1749>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix recipe_modes_of_variability (`#1753 <https://github.com/ESMValGroup/ESMValTool/pull/1753>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Flatten lists for ancestors for hyint to prevent nested lists. (`#1752 <https://github.com/ESMValGroup/ESMValTool/pull/1752>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fix bug in cmorize_obs_eppley_vgpm_modis.py (#1729) (`#1759 <https://github.com/ESMValGroup/ESMValTool/pull/1759>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Correct mip for clltkisccp in example derive preprocessor recipe (`#1768 <https://github.com/ESMValGroup/ESMValTool/pull/1768>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update date conversion in recipe_hype.yml (`#1769 <https://github.com/ESMValGroup/ESMValTool/pull/1769>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix recipe_correlation.yml (`#1767 <https://github.com/ESMValGroup/ESMValTool/pull/1767>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add attribute positive: down to plev coordinate in ERA-Interim CMORizer (`#1771 <https://github.com/ESMValGroup/ESMValTool/pull/1771>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix sispeed in recipe_preprocessor_derive_test (`#1772 <https://github.com/ESMValGroup/ESMValTool/pull/1772>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix extreme events and extreme index ancestors (`#1774 <https://github.com/ESMValGroup/ESMValTool/pull/1774>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Correct date in output filenames of ERA5 CMORizer recipe (`#1773 <https://github.com/ESMValGroup/ESMValTool/pull/1773>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Exclude WOA from multi-model stats in recipe_ocean_bgc (`#1778 <https://github.com/ESMValGroup/ESMValTool/pull/1778>`__) `Mattia Righi <https://github.com/mattiarighi>`__

Diagnostics
~~~~~~~~~~~

-  Enhancement of the hyint recipe to include etccdi indices (`#1133 <https://github.com/ESMValGroup/ESMValTool/pull/1133>`__) `Enrico Arnone <https://github.com/earnone>`__
-  Add lazy regridding for wflow diagnostic (`#1630 <https://github.com/ESMValGroup/ESMValTool/pull/1630>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Miles default domains to include lat=0 (`#1626 <https://github.com/ESMValGroup/ESMValTool/pull/1626>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  Miles: selection of reference dataset based on experiment (`#1632 <https://github.com/ESMValGroup/ESMValTool/pull/1632>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  New recipe/diagnostic:  recipe_li17natcc.yml for Axels GMD Paper (`#1567 <https://github.com/ESMValGroup/ESMValTool/pull/1567>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  New recipe/diagnostics: recipe_deangelis_for_gmdpart4.yml for Axels GMD Paper (`#1576 <https://github.com/ESMValGroup/ESMValTool/pull/1576>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  EWaterCycle: Add recipe to prepare input for LISFLOOD (`#1298 <https://github.com/ESMValGroup/ESMValTool/pull/1298>`__) `Stefan Verhoeven <https://github.com/sverhoeven>`__
-  Use area weighted regridding in wflow diagnostic (`#1643 <https://github.com/ESMValGroup/ESMValTool/pull/1643>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Workaround for permetrics recipe until Iris3 (`#1674 <https://github.com/ESMValGroup/ESMValTool/pull/1674>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  C3S_511_MPQB_bas-features (`#1465 <https://github.com/ESMValGroup/ESMValTool/pull/1465>`__) `bascrezee <https://github.com/bascrezee>`__
-  Additional Land perfmetrics (`#1641 <https://github.com/ESMValGroup/ESMValTool/pull/1641>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Necessary diagnostic from eyring06jgr for the release of version2 (`#1686 <https://github.com/ESMValGroup/ESMValTool/pull/1686>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Drought characteristics based on Martin2018 and SPI for gmd paper (`#1689 <https://github.com/ESMValGroup/ESMValTool/pull/1689>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Additional features and bugfixes for recipe anav13clim (`#1723 <https://github.com/ESMValGroup/ESMValTool/pull/1723>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Gmd laueretal2020 revisions (`#1725 <https://github.com/ESMValGroup/ESMValTool/pull/1725>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Wenzel16nature (`#1692 <https://github.com/ESMValGroup/ESMValTool/pull/1692>`__) `zechlau <https://github.com/zechlau>`__
-  Add mask albedolandcover (`#1673 <https://github.com/ESMValGroup/ESMValTool/pull/1673>`__) `bascrezee <https://github.com/bascrezee>`__
-  IPCC AR5 fig. 9.3 (seasonality) (`#1726 <https://github.com/ESMValGroup/ESMValTool/pull/1726>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Added additional emergent constraints on ECS (`#1585 <https://github.com/ESMValGroup/ESMValTool/pull/1585>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  A diagnostic to evaluate the turnover times of land ecosystem carbon (`#1395 <https://github.com/ESMValGroup/ESMValTool/pull/1395>`__) `koir-su <https://github.com/koir-su>`__
-  Removed multi_model_statistics step in recipe_oceans_example.yml as a workaround (`#1779 <https://github.com/ESMValGroup/ESMValTool/pull/1779>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Documentation
~~~~~~~~~~~~~

-  Extend getting started instructions to obtain config-user.yml (`#1642 <https://github.com/ESMValGroup/ESMValTool/pull/1642>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Extend information about native6 support on RTD (`#1652 <https://github.com/ESMValGroup/ESMValTool/pull/1652>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Update citation of ESMValTool paper in the doc (`#1664 <https://github.com/ESMValGroup/ESMValTool/pull/1664>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Updated references to documentation (now docs.esmvaltool.org) (`#1679 <https://github.com/ESMValGroup/ESMValTool/pull/1679>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Replace dead link with ESGF link. (`#1681 <https://github.com/ESMValGroup/ESMValTool/pull/1681>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Add all European grants to Zenodo (`#1682 <https://github.com/ESMValGroup/ESMValTool/pull/1682>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update Sphinx to v3 or later (`#1685 <https://github.com/ESMValGroup/ESMValTool/pull/1685>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Small fix to number of models in ensclus documentation (`#1691 <https://github.com/ESMValGroup/ESMValTool/pull/1691>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  Move draft_release_notes.py from ESMValCore to here and update (`#1701 <https://github.com/ESMValGroup/ESMValTool/pull/1701>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve the installation instructions (`#1634 <https://github.com/ESMValGroup/ESMValTool/pull/1634>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Improve description of how to implement provenance in diagnostic (`#1750 <https://github.com/ESMValGroup/ESMValTool/pull/1750>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update command line interface documentation and add links to ESMValCore configuration documentation (`#1776 <https://github.com/ESMValGroup/ESMValTool/pull/1776>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Documentation on how to find shapefiles for hydrology recipes (`#1777 <https://github.com/ESMValGroup/ESMValTool/pull/1777>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__

Improvements
~~~~~~~~~~~~

-  Pin flake8<3.8.0 (`#1635 <https://github.com/ESMValGroup/ESMValTool/pull/1635>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update conda package path in more places (`#1636 <https://github.com/ESMValGroup/ESMValTool/pull/1636>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Remove curly brackets around issue number in pull request template (`#1637 <https://github.com/ESMValGroup/ESMValTool/pull/1637>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix style issue in test (`#1639 <https://github.com/ESMValGroup/ESMValTool/pull/1639>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update Codacy badges (`#1662 <https://github.com/ESMValGroup/ESMValTool/pull/1662>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Support extra installation methods in R (`#1360 <https://github.com/ESMValGroup/ESMValTool/pull/1360>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add ncdf4.helpers package as a dependency again (`#1678 <https://github.com/ESMValGroup/ESMValTool/pull/1678>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Speed up conda installation (`#1677 <https://github.com/ESMValGroup/ESMValTool/pull/1677>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update CMORizers and recipes for ESMValCore v2.0.0 (`#1699 <https://github.com/ESMValGroup/ESMValTool/pull/1699>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update setup.py for PyPI package (`#1700 <https://github.com/ESMValGroup/ESMValTool/pull/1700>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Cleanup recipe headers before the release (`#1740 <https://github.com/ESMValGroup/ESMValTool/pull/1740>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-    Add colortables as esmvaltool subcommand (`#1666 <https://github.com/ESMValGroup/ESMValTool/pull/1666>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Increase version to v2.0.0 (`#1756 <https://github.com/ESMValGroup/ESMValTool/pull/1756>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update job script (`#1757 <https://github.com/ESMValGroup/ESMValTool/pull/1757>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Read authors and description from .zenodo.json (`#1758 <https://github.com/ESMValGroup/ESMValTool/pull/1758>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update docker recipe to install from source (`#1651 <https://github.com/ESMValGroup/ESMValTool/pull/1651>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorize aphro ma (`#1555 <https://github.com/ESMValGroup/ESMValTool/pull/1555>`__) `mwjury <https://github.com/mwjury>`__
-  Respectable testing for cmorizers/obs/utilities.py and cmorizers/obs/cmorize_obs.py (`#1517 <https://github.com/ESMValGroup/ESMValTool/pull/1517>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix start year in recipe_check_obs (`#1638 <https://github.com/ESMValGroup/ESMValTool/pull/1638>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Cmorizer for the PERSIANN-CDR precipitation data (`#1633 <https://github.com/ESMValGroup/ESMValTool/pull/1633>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Cmorize eobs (`#1554 <https://github.com/ESMValGroup/ESMValTool/pull/1554>`__) `mwjury <https://github.com/mwjury>`__
-  Update download cds satellite lai fapar (`#1654 <https://github.com/ESMValGroup/ESMValTool/pull/1654>`__) `bascrezee <https://github.com/bascrezee>`__
-  Added monthly mean vars (ta, va, zg) to era5 cmorizer via recipe (`#1644 <https://github.com/ESMValGroup/ESMValTool/pull/1644>`__) `Evgenia Galytska <https://github.com/egalytska>`__
-  Make format time check more flexible (`#1661 <https://github.com/ESMValGroup/ESMValTool/pull/1661>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Exclude od550lt1aer from recipe_check_obs.yml (`#1720 <https://github.com/ESMValGroup/ESMValTool/pull/1720>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  PERSIANN-CDR cmorizer update: adding the capability to save monthly mean files (`#1728 <https://github.com/ESMValGroup/ESMValTool/pull/1728>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Add standard_name attribute to lon and lat in cmorize_obs_esacci_oc.py (`#1760 <https://github.com/ESMValGroup/ESMValTool/pull/1760>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Allow for incomplete months on daily frequency in cmorizer ncl utilities (`#1754 <https://github.com/ESMValGroup/ESMValTool/pull/1754>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix AURA-TES cmorizer (`#1766 <https://github.com/ESMValGroup/ESMValTool/pull/1766>`__) `Mattia Righi <https://github.com/mattiarighi>`__

.. _changelog-v2-0-0b4:

v2.0.0b4
--------

This release includes

Bug fixes
~~~~~~~~~

-  Fix HALOE plev coordinate (`#1590 <https://github.com/ESMValGroup/ESMValTool/pull/1590>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix tro3 units in HALOE (`#1591 <https://github.com/ESMValGroup/ESMValTool/pull/1591>`__) `Mattia Righi <https://github.com/mattiarighi>`__

Diagnostics
~~~~~~~~~~~

-  Applicate sea ice negative feedback (`#1299 <https://github.com/ESMValGroup/ESMValTool/pull/1299>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add Russell18jgr ocean diagnostics (`#1592 <https://github.com/ESMValGroup/ESMValTool/pull/1592>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Refactor marrmot recipe and diagnostic to use ERA5 daily data made by new cmorizer (`#1600 <https://github.com/ESMValGroup/ESMValTool/pull/1600>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  In recipe_wflow, use daily ERA5 data from the new cmorizer. (`#1599 <https://github.com/ESMValGroup/ESMValTool/pull/1599>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  In wflow diagnostic, calculate PET after(!) interpolation and lapse rate correction (`#1618 <https://github.com/ESMValGroup/ESMValTool/pull/1618>`__) `Jerom Aerts <https://github.com/jeromaerts>`__
-  Fixed wenz14jgr (`#1562 <https://github.com/ESMValGroup/ESMValTool/pull/1562>`__) `zechlau <https://github.com/zechlau>`__
-  Update portrait_plot.ncl (`#1625 <https://github.com/ESMValGroup/ESMValTool/pull/1625>`__) `Bettina Gier <https://github.com/bettina-gier>`__

Documentation
~~~~~~~~~~~~~

-  Restructure documentation (`#1587 <https://github.com/ESMValGroup/ESMValTool/pull/1587>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more links to documentation (`#1595 <https://github.com/ESMValGroup/ESMValTool/pull/1595>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update links in readme (`#1598 <https://github.com/ESMValGroup/ESMValTool/pull/1598>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Minor improvements to installation documentation (`#1608 <https://github.com/ESMValGroup/ESMValTool/pull/1608>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add info for new mailing list to documentation. (`#1607 <https://github.com/ESMValGroup/ESMValTool/pull/1607>`__) `Björn Brötz <https://github.com/bjoernbroetz>`__
-  Update making a release documentation (`#1627 <https://github.com/ESMValGroup/ESMValTool/pull/1627>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Improvements
~~~~~~~~~~~~

-  Avoid broken pytest-html plugin (`#1583 <https://github.com/ESMValGroup/ESMValTool/pull/1583>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Remove reference section in config-references.yml (`#1545 <https://github.com/ESMValGroup/ESMValTool/pull/1545>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Various improvements to development infrastructure (`#1570 <https://github.com/ESMValGroup/ESMValTool/pull/1570>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Install scikit-learn from conda, remove libunwind as a direct dependency (`#1611 <https://github.com/ESMValGroup/ESMValTool/pull/1611>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Create conda subpackages and enable tests (`#1624 <https://github.com/ESMValGroup/ESMValTool/pull/1624>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorizer for HALOE (`#1581 <https://github.com/ESMValGroup/ESMValTool/pull/1581>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Add CMORizer for CT2019 (`#1604 <https://github.com/ESMValGroup/ESMValTool/pull/1604>`__) `Manuel Schlund <https://github.com/schlunma>`__

For older releases, see the release notes on https://github.com/ESMValGroup/ESMValTool/releases.
Recipe
******

Writing a basic recipe
======================
The user will need to write a basic recipe to be able to run their own personal diagnostic.
An example of such a recipe is found in `esmvaltool/recipes/recipe_my_personal_diagnostic.yml`.
For general guidelines with regards to ESMValTool recipes please consult the User Guide;
the specific parameters needed by a recipe that runs a personal diagnostic are:

.. code-block:: yaml

  scripts:
    my_diagnostic:
    script: /path/to/your/my_little_diagnostic.py

i.e. the full path to the personal diagnostic that the user needs to run.

There is also a lesson available in the 
`ESMValTool tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`_
that describes in a step-by-step procedure how to write your own recipe. It can be found
`here <https://esmvalgroup.github.io/ESMValTool_Tutorial/06-preprocessor/index.html>`_.
.. _new-cmorizer:

Writing a CMORizer script for an additional dataset
***************************************************

ESMValTool is designed to work with `CF compliant <http://cfconventions.org/>`_
data and follows the CMOR tables from the CMIP data request, therefore
the observational datasets need to be CMORized for usage in ESMValTool.
The following steps are necessary to prepare an observational
data set for the use in ESMValTool.

| `1. Check if your variable is CMOR standard`_
| `2. Edit your configuration file`_
| `3. Store your dataset in the right place`_
| `4. Create a cmorizer for the dataset`_
| `4.1 Cmorizer script written in python`_
| `4.2 Cmorizer script written in NCL`_
| `5. Run the cmorizing script`_
| `6. Naming convention of the observational data files`_
| `7. Test the cmorized dataset`_

.. note::
  **CMORization as a fix.** As of early 2020, we've started implementing cmorization as
  *fixes*. As compared to the workflow described below, this has the advantage that
  the user does not need to store a duplicate (CMORized) copy of the data. Instead, the
  CMORization is performed 'on the fly' when running a recipe. **ERA5** is the first dataset
  for which this 'CMORization on the fly' is supported. For more information, see:
  :ref:`cmorization_as_fix`.


1. Check if your variable is CMOR standard
==========================================

Most variables are defined in the CMIP data request and can be found in the
CMOR tables in the folder `/esmvalcore/cmor/tables/cmip6/Tables/
<https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/cmor/tables/cmip6/Tables>`_,
differentiated according to the ``MIP`` they belong to. The tables are a
copy of the `PCMDI <https://github.com/PCMDI>`_ guidelines. If you find the
variable in one of these tables, you can proceed to the next section.

If your variable is not available in the standard CMOR tables,
you need to write a custom CMOR table entry for the variable
as outlined below and add it to `/esmvalcore/cmor/tables/custom/
<https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/cmor/tables/custom>`_.

To create a new custom CMOR table you need to follow these
guidelines:

- Provide the ``variable_entry``;
- Provide the ``modeling_realm``;
- Provide the variable attributes, but leave ``standard_name`` blank. Necessary
  variable attributes are: ``units``, ``cell_methods``, ``cell_measures``,
  ``long_name``, ``comment``.
- Provide some additional variable attributes. Necessary additional variable
  attributes are: ``dimensions``, ``out_name``, ``type``. There are also
  additional variable attributes that can be defined here (see the already
  available cmorizers).

It is recommended to use an existing custom table as a template, to edit the
content and save it as ``CMOR_<short_name>.dat``.

2. Edit your configuration file
===============================

Make sure that beside the paths to the model simulations and observations, also
the path to raw observational data to be cmorized (``RAWOBS``) is present in
your configuration file.

3. Store your dataset in the right place
========================================

The folder ``RAWOBS`` needs the subdirectories ``Tier1``, ``Tier2`` and
``Tier3``. The different tiers describe the different levels of restrictions
for downloading (e.g. providing contact information, licence agreements)
and using the observations. The unformatted (raw) observations
should then be stored then in the appropriate of these three folders.

4. Create a cmorizer for the dataset
====================================

There are many cmorizing scripts available in `/esmvaltool/cmorizers/obs/
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/>`_
where solutions to many kinds of format issues with observational data are
addressed. Most of these scripts are written in NCL at the moment, but more
and more examples for Python-based cmorizing scripts become available.

.. note::
  NCL support will terminate soon, so new cmorizer scripts should preferably be
  written in Python.

How much cmorizing an observational data set needs is strongly dependent on
the original NetCDF file and how close the original formatting already is to
the strict CMOR standard.

In the following two subsections two cmorizing scripts, one written in Python
and one written in NCL, are explained in more detail.

4.1 Cmorizer script written in python
-------------------------------------

Find here an example of a cmorizing script, written for the ``MTE`` dataset
that is available at the MPI for Biogeochemistry in Jena: `cmorize_obs_mte.py
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/cmorize_obs_mte.py>`_.

All the necessary information about the dataset to write the filename
correctly, and which variable is of interest, is stored in a separate
configuration file: `MTE.yml
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/cmor_config/MTE.yml>`_
in the directory ``ESMValTool/esmvaltool/cmorizers/obs/cmor_config/``. Note
that the name of this configuration file has to be identical to the name of
your data set. It is recommended that you set ``project`` to ``OBS6`` in the
configuration file. That way, the variables defined in the CMIP6 CMOR table,
augmented with the custom variables described above, are available to your script.

The first part of this configuration file defines the filename of the raw
observations file. The second part defines the common global attributes for
the cmorizer output, e.g. information that is needed to piece together the
final observations file name in the correct structure (see Section `6. Naming convention of the observational data files`_).
Another global attribute is ``reference`` which includes a ``doi`` related to the dataset.
Please see the section `adding references
<https://docs.esmvaltool.org/en/latest/community/diagnostic.html#adding-references>`_
on how to add reference tags to the ``reference`` section in the configuration file.
If a single dataset has more than one reference,
it is possible to add tags as a list e.g. ``reference: ['tag1', 'tag2']``.
The third part in the configuration file defines the variables that are supposed to be cmorized.

The actual cmorizing script ``cmorize_obs_mte.py`` consists of a header with
information on where and how to download the data, and noting the last access
of the data webpage.

The main body of the CMORizer script must contain a function called

.. code-block:: python

   def cmorization(in_dir, out_dir, cfg, config_user):

with this exact call signature. Here, ``in_dir`` corresponds to the input
directory of the raw files, ``out_dir`` to the output directory of final
reformatted data set and ``cfg`` to the configuration dictionary given by
the  ``.yml`` configuration file. The return value of this function is ignored. All
the work, i.e. loading of the raw files, processing them and saving the final
output, has to be performed inside its body. To simplify this process, ESMValTool
provides a set of predefined utilities.py_, which can be imported into your CMORizer
by

.. code-block:: python

   from . import utilities as utils

Apart from a function to easily save data, this module contains different kinds
of small fixes to the data attributes, coordinates, and metadata which are
necessary for the data field to be CMOR-compliant.

Note that this specific CMORizer script contains several subroutines in order
to make the code clearer and more readable (we strongly recommend to follow
that code style). For example, the function ``_get_filepath`` converts the raw
filepath to the correct one and the function ``_extract_variable`` extracts and
saves a single variable from the raw data.

.. _utilities.py: https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/utilities.py


4.2 Cmorizer script written in NCL
----------------------------------

Find here an example of a cmorizing script, written for the ``ESACCI XCH4``
dataset that is available on the Copernicus Climate Data Store:
`cmorize_obs_cds_xch4.ncl
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/cmorize_obs_cds_xch4.ncl>`_.

The first part of the script collects all the information about the dataset
that are necessary to write the filename correctly and to understand which
variable is of interest here. Please make sure to provide the correct
information for following key words: DIAG_SCRIPT, VAR, NAME, MIP, FREQ,
CMOR_TABLE.

- **Note:** the fields ``VAR``, ``NAME``, ``MIP`` and ``FREQ`` all ask for one
  or more entries. If more than one entry is provided, make sure that the order
  of the entries is the same for all four fields! (for example, that the first
  entry in all four fields describe the variable ``xch4`` that you would like
  to extract);
- **Note:** some functions in the script are NCL-specific and are available
  through the loading of the script interface.ncl_. There are similar
  functions available for python scripts.

.. _interface.ncl: https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/interface.ncl

.. _utilities.ncl: https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/utilities.ncl

In the second part of the script each variable defined in ``VAR`` is separately
extracted from the original data file and processed. Most parts of the code are
commented, and therefore it should be easy to follow. ESMValTool provides a set
of predefined utilities.ncl_, which can be imported into your CMORizer
by

.. code-block:: NCL

   loadscript(getenv("esmvaltool_root") + "/esmvaltool/cmorizers/obs/utilities.ncl")

This module contains different kinds of small fixes to the data attributes,
coordinates, and metadata which are necessary for the data field to be
CMOR-compliant.

5. Run the cmorizing script
===========================

The cmorizing script for the given dataset can be run with:

.. code-block:: console

 cmorize_obs -c <config-user.yml> -o <dataset-name>


.. note::

   The output path given in the configuration file is the path where
   your cmorized dataset will be stored. The ESMValTool will create a folder
   with the correct tier information (see Section `2. Edit your configuration file`_) if that tier folder is not
   already available, and then a folder named after the data set. In this
   folder the cmorized data set will be stored as a netCDF file.

If your run was successful, one or more NetCDF files are produced in your
output directory.


6. Naming convention of the observational data files
====================================================

For the ESMValTool to be able to read the observations from the NetCDF file,
the file name needs a very specific structure and order of information parts
(very similar to the naming convention for observations in ESMValTool
v1.0). The file name will be automatically correctly created if a cmorizing
script has been used to create the netCDF file.

The correct structure of an observational data set is defined in
`config-developer.yml
<https://github.com/ESMValGroup/ESMValCore/blob/main/esmvalcore/config-developer.yml>`_,
and looks like the following:

.. code-block:: console

  OBS_[dataset]_[type]_[version]_[mip]_[short_name]_YYYYMM-YYYYMM.nc

For the example of the ``CDS-XCH4`` data set, the correct structure of the
file name looks then like this:

.. code-block:: console

  OBS_CDS-XCH4_sat_L3_Amon_xch4_200301-201612.nc

The different parts of the name are explained in more detail here:

- OBS: describes what kind of data can be expected in the file, in this case
  ``observations``;
- CDS-XCH4: that is the name of the dataset. It has been named this way for
  illustration purposes (so that everybody understands it is the xch4 dataset
  downloaded from the CDS), but a better name would indeed be ``ESACCI-XCH4``
  since it is a ESA-CCI dataset;
- sat: describes the source of the data, here we are looking at satellite data
  (therefore ``sat``), could also be ``reanaly`` for reanalyses;
- L3: describes the version of the dataset:
- Amon: is the information in which ``mip`` the variable is to be expected, and
  what kind of temporal resolution it has; here we expect ``xch4`` to be part
  of the atmosphere (``A``) and we have the dataset in a monthly resolution
  (``mon``);
- xch4: Is the name of the variable. Each observational data file is supposed
  to only include one variable per file;
- 200301-201612: Is the period the dataset spans with ``200301`` being the
  start year and month, and ``201612`` being the end year and month;

.. note::
   There is a different naming convention for ``obs4MIPs`` data (see the exact
   specifications for the obs4MIPs data file naming convention in the
   ``config-developer.yml`` file).

7. Test the cmorized dataset
======================================

To verify that the cmorized data file is indeed correctly formatted, you can
run a dedicated test recipe, that does not include any diagnostic, but only
reads in the data file and has it processed in the preprocessor. Such a recipe
is called ``recipes/examples/recipe_check_obs.yml``. You just need to add a
diagnostic for your dataset following the existing entries.
Only the diagnostic of interest needs to be run, the others should be commented
out for testing.
Making a recipe or diagnostic
=============================

.. toctree::
   :maxdepth: 1

		Introduction <introduction>
		Recipe <recipe>
		Diagnostic <diagnostic>
		Dataset <dataset>
Introduction
============

This chapter contains instructions for developing your own recipes and/or diagnostics.
It also contains a section describing how to use additional datasets with ESMValTool.
While it is possible to use just the ESMValCore package and run any recipes/diagnostics you develop with just this package, it is highly recommended that you consider contributing the work you do back to the ESMValTool community.
Among the advantages of contributing to the community are improved visibility of your work and support by the community with making and maintaining your diagnostic.
See the :ref:`Community <community>` chapter for a guide on how to contribute to the community.
Diagnostic
**********

Instructions for personal diagnostic
====================================

Anyone can run a personal diagnostic, no matter where the location of it;
there is no need to install esmvaltool in developer mode nor is it to
git push or for that matter, do any git operations; the example recipe

.. code-block:: console

  esmvaltool/recipes/recipe_my_personal_diagnostic.yml

shows the use of running a personal diagnostic; the example

.. code-block:: console

  esmvaltool/diag_scripts/examples/my_little_diagnostic.py

and any of its alterations may be used as training wheels for the future ESMValTool
diagnostic developer. The purpose of this example is to familiarize the user with
the framework of ESMValTool without the constraints of installing and running the
tool as developer.

Functionality
=============

`my_little_diagnostic` (or whatever the user will call their diagnostic) makes full use
of ESMValTool's preprocessor output (both phyisical files and run variables); this output
comes in form of a nested dictionary, or config dictionary, see an example below;
it also makes full use of the ability to call any of the preprocessor's functions,
note that relative imports of modules from the esmvaltool package are allowed and
work without altering the $PYTHONPATH.

The user may parse this dictionary so that they execute a number of operations on the
preprocessed data; for example the `my_little_diagnostic.plot_time_series` grabs the
preprocessed data output, computes global area averages for each model, then plots
a time-series for each model. Different manipulation functionalities for grouping,
sorting etc of the data in the config dictionary are available,
please consult ESMValTool User Manual.


Example of config dictionary
============================
To be added (use python-style code-block).
.. _contributing_code_docs:

Contributing code and documentation
===================================

If you would like to contribute a new diagnostic and recipe or a new feature,
please discuss your idea with the development team before getting started, to
avoid double work and/or disappointment later.
A good way to do this is to open an
`issue on GitHub <https://github.com/ESMValGroup/ESMValTool/issues>`__.
This is also a good way to get help with the implementation.

We value the time you invest in contributing and strive to make the process as
easy as possible.
If you have suggestions for improving the process of contributing, please do
not hesitate to propose them, for example by starting a discussion on our
`discussions page <https://github.com/ESMValGroup/ESMValTool/discussions>`__.

Getting started
---------------

See :ref:`install_from_source` for instructions on how to set up a development
installation.

New development should preferably be done in the
`ESMValTool <https://github.com/ESMValGroup/ESMValTool>`__
GitHub repository.
However, for scientists requiring confidentiality, private repositories are
available, see :ref:`private_repository` for more information.
The default git branch is ``main``. Use
this branch to create a new feature branch from and make a pull request
against.
This
`page <https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow>`__
offers a good introduction to git branches, but it was written for
BitBucket while we use GitHub, so replace the word BitBucket by GitHub
whenever you read it.

It is recommended that you open a `draft pull
request <https://github.blog/2019-02-14-introducing-draft-pull-requests/>`__
early, as this will cause :ref:`CircleCI to run the unit tests <tests>`,
:ref:`Codacy to analyse your code <code_quality>`, and
:ref:`readthedocs to build the documentation <documentation>`.
It’s also easier to get help from other developers if
your code is visible in a pull request.

Please review the results of the automatic checks below your pull request.
If one of the tests shows a red cross instead of a green checkmark, please click
the ``Details`` link behind the failing check and try to solve the issue.
Ask `@ESMValGroup/tech-reviewers`_ for help if you do not know how to fix the
failing check.
Note that this kind of automated checks make it easier to
:ref:`review code <reviewing>`, but they are not flawless.
Preferably Codacy code quality checks pass, however a few remaining hard to
solve Codacy issues are still acceptable.
If you suspect Codacy may be wrong, please ask by commenting on your pull
request.

.. _pull_request_checklist:

Checklist for pull requests
---------------------------

To clearly communicate up front what is expected from a pull request, we have
the following checklist.
Please try to do everything on the list before requesting a review.
If you are unsure about something on the list, please ask the
`@ESMValGroup/tech-reviewers`_ or `@ESMValGroup/science-reviewers`_ for help
by commenting on your (draft) pull request or by starting a new
`discussion <https://github.com/ESMValGroup/ESMValTool/discussions>`__.

In the ESMValTool community we use
:ref:`pull request reviews <reviewing>` to ensure all code and
documentation contributions are of good quality.
The icons indicate whether the item will be checked during the
:ref:`🛠 Technical review <technical_review>` or
:ref:`🧪 Scientific review <scientific_review>`.

All pull requests
~~~~~~~~~~~~~~~~~

- 🛠 :ref:`The pull request has a descriptive title <descriptive_pr_title>`
- 🛠 Code is written according to the :ref:`code quality guidelines <code_quality>`
- 🛠 Documentation_ is available
- 🛠 Tests_ run successfully
- 🛠 The :ref:`list of authors <authors>` is up to date
- 🛠 Changed dependencies are :ref:`added or removed correctly <dependencies>`
- 🛠 The :ref:`checks shown below the pull request <pull_request_checks>` are successful

New or updated recipe and/or diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See :ref:`new-diagnostic` for detailed instructions.

- 🧪 :ref:`Recipe runs successfully <testing_recipes>`
- 🧪 :ref:`recipe_documentation` is available
- 🧪 :ref:`Figure(s) and data <diagnostic_output>` look as expected from literature
- 🛠 :ref:`Provenance information <recording-provenance>` has been added

New or updated data reformatting script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See :ref:`new dataset <new-dataset>` for detailed instructions.

- 🛠 :ref:`dataset-documentation` is available
- 🛠 The dataset has been :ref:`added to the CMOR check recipe <dataset-test>`
- 🧪 Numbers and units of the data look :ref:`physically meaningful <dataset-sanity-check>`

.. _descriptive_pr_title:

Pull request title
------------------

The title of a pull request should clearly describe what the pull request changes.
If you need more text to describe what the pull request does, please add it in
the description.
The titles of pull requests are used to compile the :ref:`changelog`, therefore
it is important that they are easy to understand for people who are not
familiar with the code or people in the project.
Descriptive pull request titles also makes it easier to find back what was
changed when, which is useful in case a bug was introduced.

.. _code_quality:

Code quality
------------

To increase the readability and maintainability or the ESMValTool source
code, we aim to adhere to best practices and coding standards.
For code in all languages, it is highly recommended that you split your code up
in functions that are short enough to view without scrolling, e.g. no more than
50 lines long.

We include checks for Python, R, NCL, and yaml files, most of which are
described in more detail in the sections below.
This includes checks for invalid syntax and formatting errors.
:ref:`pre-commit` is a handy tool that can run all of these checks automatically
just before you commit your code.
It knows knows which tool to run for each filetype, and therefore provides
a convenient way to check your code.

Python
~~~~~~

The standard document on best practices for Python code is
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ and there is
`PEP257 <https://www.python.org/dev/peps/pep-0257/>`__ for code documentation.
We make use of
`numpy style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`__
to document Python functions that are visible on
`readthedocs <https://docs.esmvaltool.org>`__.

To check if your code adheres to the standard, go to the directory where
the repository is cloned, e.g. ``cd ESMValTool``, and run `prospector <http://prospector.landscape.io/>`_

::

   prospector esmvaltool/diag_scripts/your_diagnostic/your_script.py

In addition to prospector, we also use `flake8 <https://flake8.pycqa.org/en/latest/>`_
to automatically check for obvious bugs and formatting mistakes.

When you make a pull request, adherence to the Python development best practices
is checked in two ways:

#. As part of the unit tests, flake8_ is run by
   `CircleCI <https://app.circleci.com/pipelines/github/ESMValGroup/ESMValTool>`_,
   see the section on Tests_ for more information.
#. `Codacy <https://app.codacy.com/gh/ESMValGroup/ESMValTool/pullRequests>`_
   is a service that runs prospector (and other code quality tools) on changed
   files and reports the results.
   Click the 'Details' link behind the Codacy check entry and then click
   'View more details on Codacy Production' to see the results of the static
   code analysis done by Codacy_.
   If you need to log in, you can do so using your GitHub account.

A pull request should preferably not introduce any new prospector issues.
However, we understand that there is a limit to how much time can be spent on
polishing code, so up to 10 new (non-trivial) issues is still an acceptable
amount.
Formatting issues are considered trivial and need to be addressed.
Note that the automatic code quality checks by prospector are really helpful to
improve the quality of your code, but they are not flawless.
If you suspect prospector or Codacy may be wrong, please ask the
`@ESMValGroup/tech-reviewers`_ by commenting on your pull request.

Note that running prospector locally will give you quicker and sometimes more
accurate results than waiting for Codacy.

Most formatting issues in Python code can be fixed automatically by
running the commands

::

   isort some_file.py

to sort the imports in `the standard way <https://www.python.org/dev/peps/pep-0008/#imports>`__
using `isort <https://pycqa.github.io/isort/>`__ and

::

   yapf -i some_file.py

to add/remove whitespace as required by the standard using `yapf <https://github.com/google/yapf>`__,

::

   docformatter -i some_file.py

to run `docformatter <https://github.com/myint/docformatter>`__ which helps
formatting the docstrings (such as line length, spaces).

NCL
~~~

Because there is no standard best practices document for NCL, we use
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ for NCL code as
well, with some minor adjustments to accommodate for differences in the
languages. The most important difference is that for NCL code the
indentation should be 2 spaces instead of 4.
Use the command ``nclcodestyle /path/to/file.ncl`` to check if your code
follows the style guide.
More information on the ``nclcodestyle`` command can be found
:ref:`here <nclcodestyle>`.

R
~

Best practices for R code are described in `The tidyverse style
guide <https://style.tidyverse.org/>`__. We check adherence to this
style guide by using
`lintr <https://cran.r-project.org/web/packages/lintr/index.html>`__ on
CircleCI. Please use `styler <https://styler.r-lib.org/>`__ to
automatically format your code according to this style guide. In the
future we would also like to make use of
`goodpractice <https://cran.r-project.org/web/packages/goodpractice/index.html>`__
to assess the quality of R code.

YAML
~~~~

Please use `yamllint <https://yamllint.readthedocs.io>`_ to check that your
YAML files do not contain mistakes.
``yamllint`` checks for valid syntax, common mistakes like key repetition and
cosmetic problems such as line length, trailing spaces, wrong indentation, etc.
When the tool complains about the maximum line length or too many spaces, please
use your own best judgement about whether solving the issue will make your
recipe more readable.

Any text file
~~~~~~~~~~~~~

A generic tool to check for common spelling mistakes is
`codespell <https://pypi.org/project/codespell/>`__.

.. _documentation:

Documentation
-------------

The documentation lives on `docs.esmvaltool.org <https://docs.esmvaltool.org>`_
and is built using `Sphinx <https://www.sphinx-doc.org>`_.
There are two main ways of adding documentation:

#. As written text in the directory
   `doc/sphinx/source <https://github.com/ESMValGroup/ESMValTool/tree/main/doc/sphinx/source>`__.
   When writing
   `reStructuredText <https://www.sphinx-doc.org/en/main/usage/restructuredtext/basics.html>`_
   (``.rst``) files, please try to limit the line length to 80 characters and
   always start a sentence on a new line.
   This makes it easier to review changes to documentation on GitHub.

#. As docstrings or comments in code.
   For Python code, the
   `docstrings <https://www.python.org/dev/peps/pep-0257/>`__
   of Python modules, classes, and functions
   that are mentioned in
   `doc/sphinx/source/api <https://github.com/ESMValGroup/ESMValTool/tree/main/doc/sphinx/source/api>`__
   are used to generate documentation.
   This results in the :ref:`api`.

.. _doc_howto:

What should be documented
~~~~~~~~~~~~~~~~~~~~~~~~~

See also :ref:`recipe_documentation` and :ref:`dataset-documentation`.

Any code documentation that is visible on `docs.esmvaltool.org`_
should be well written and adhere to the standards for documentation for the
respective language.
Note that there is no need to write extensive documentation for functions that
are not visible in the online documentation.
However, a short description in the docstring helps other contributors to
understand what a function is intended to do and and what its capabilities are.
For short functions, a one-line docstring is usually sufficient, but more
complex functions might require slightly more extensive documentation.

How to build and view the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whenever you make a pull request or push new commits to an existing pull
request, readthedocs will automatically build the documentation.
The link to the documentation will be shown in the list of checks below your
pull request, click 'Details' behind the check
``docs/readthedocs.org:esmvaltool`` to preview the documentation.
If all checks were successful, you may need to click 'Show all checks' to see
the individual checks.

To build the documentation on your own computer, go to the directory where the
repository was cloned and run

::

   python setup.py build_sphinx

or

::

   python setup.py build_sphinx -Ea

to build it from scratch.
Make sure that your newly added documentation builds without warnings or
errors and looks correctly formatted.
CircleCI_ will build the documentation with the command

.. code-block:: bash

   python setup.py build_sphinx --warning-is-error

to catch mistakes that can be detected automatically.

The configuration file for Sphinx_ is
`doc/shinx/source/conf.py <https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/conf.py>`_.

When reviewing a pull request, always check that the documentation checks
shown below the pull request were successful.
Successful checks have a green ✓ in front, a ❌ means the test job failed.

.. _esmvalcore-documentation-integration:

Integration with the ESMValCore documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The `ESMValCore documentation <https://docs.esmvaltool.org/projects/esmvalcore>`_
is hosted as a
`subproject <https://docs.readthedocs.io/en/stable/subprojects.html>`_
of the ESMValTool documentation on readthedocs.
To link to a section from the ESMValCore documentation from the reStructuredText
(``.rst``) files, use the usual ``:ref:`` but prefix the reference with
``esmvalcore:``.
For example, ``:ref:`esmvalcore:recipe``` to link to
:ref:`esmvalcore:recipe`.

There is a script that generates the navigation menu shown on the left when
you view the documentation.
This script is called
`doc/sphinx/source/gensidebar.py <https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/gensidebar.py>`_
in the ESMValTool repository and it should be identical to
`doc/gensidebar.py <https://github.com/ESMValGroup/ESMValCore/blob/main/doc/gensidebar.py>`_
in the ESMValCore repository, or the sidebar will change when navigating from
the ESMValTool documentation to the ESMValCore documentation and vice-versa.

.. _tests:

Tests
-----

To check various aspects of the recipes and code, there tests available in the
`tests <https://github.com/ESMValGroup/ESMValTool/tree/main/tests>`__
directory.

Whenever you make a pull request or push new commits to an existing pull
request, these tests will be run automatically on CircleCI_.
The results appear at the bottom of the pull request.
Click on 'Details' for more information on a specific test job.
To see some of the results on CircleCI, you may need to log in.
You can do so using your GitHub account.

To run the tests on your own computer, go to the directory where the repository
is cloned and run the command ``pytest``.

Have a look at :ref:`testing_recipes` for information on testing recipes.

Every night, more extensive tests are run to make sure that problems with the
installation of the tool are discovered by the development team before users
encounter them.
These nightly tests have been designed to mimic the installation procedures
described in the documentation, e.g. in the :ref:`install` chapter.
The nightly tests are run using both CircleCI and GitHub Actions, the
result of the tests ran by CircleCI can be seen on the
`CircleCI project page <https://app.circleci.com/pipelines/github/ESMValGroup/ESMValTool?branch=main>`__
and the result of the tests ran by GitHub Actions can be viewed on the
`Actions tab <https://github.com/ESMValGroup/ESMValTool/actions>`__
of the repository.

The configuration of the tests run by CircleCI can be found in the directory
`.circleci <https://github.com/ESMValGroup/ESMValTool/blob/main/.circleci>`__,
while the configuration of the tests run by GitHub Actions can be found in the
directory
`.github/workflows <https://github.com/ESMValGroup/ESMValTool/blob/main/.github/workflows>`__.

When reviewing a pull request, always check that all test jobs on CircleCI_ were
successful.
Successful test jobs have a green ✓ in front, a ❌ means the test job failed.

.. _authors:

List of authors
---------------

If you make a contribution to ESMValTool and you would like to be listed as an
author (e.g. on `Zenodo <https://zenodo.org/record/4562215>`__), please add your
name to the list of authors in ``CITATION.cff`` and generate the entry for the
``.zenodo.json`` file by running the commands

::

   pip install cffconvert
   cffconvert --ignore-suspect-keys --outputformat zenodo --outfile .zenodo.json

Note that authors of recipes and/or diagnostics also need to be added to the file
`esmvaltool/config-references.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/config-references.yml>`__,
see :ref:`recording-provenance` for more information.

.. _dependencies:

Dependencies
------------

Before considering adding a new dependency, carefully check that the
`license <https://the-turing-way.netlify.app/reproducible-research/licensing/licensing-software.html>`__
of the dependency you want to add and any of its dependencies are
`compatible <https://the-turing-way.netlify.app/reproducible-research/licensing/licensing-compatibility.html>`__
with the
`Apache 2.0 <https://github.com/ESMValGroup/ESMValTool/blob/main/LICENSE/>`_
license that applies to the ESMValTool.
Note that GPL version 2 license is considered incompatible with the Apache 2.0
license, while the compatibility of GPL version 3 license with the Apache 2.0
license is questionable.
See this `statement <https://www.apache.org/licenses/GPL-compatibility.html>`__
by the authors of the Apache 2.0 license for more information.

When adding or removing dependencies, please consider applying the changes in
the following files:

- ``environment.yml``
  contains development dependencies that cannot be installed from
  `PyPI <https://pypi.org/>`__/`CRAN <https://cran.r-project.org/>`__/`Julia package registry <https://github.com/JuliaRegistries/General>`__
- ``environment_osx.yml``
  contains development dependencies for MacOSX. Should be the same as ``environment.yml``,
  but currently without multi language support.
- ``docs/sphinx/source/requirements.txt``
  contains Python dependencies needed to build the documentation that can be
  installed from PyPI
- ``docs/sphinx/source/conf.py``
  contains a list of Python dependencies needed to build the documentation that
  cannot be installed from PyPI and need to be mocked when building the
  documentation.
  (We do not use conda to build the documentation because this is too time
  consuming.)
- ``esmvaltool/install/R/r_requirements.txt``
  contains R dependencies that can be installed from CRAN
- ``esmvaltool/install/Julia/Project.toml``
  contains Julia dependencies that can be installed from the default Julia
  package registry
- ``setup.py``
  contains all Python dependencies, regardless of their installation source
- ``package/meta.yaml``
  contains dependencies for the conda package; all Python and compiled
  dependencies that can be installed from conda should be listed here, but no Julia
  dependencies because doing that would make it impossible to solve the conda environment

Note that packages may have a different name on
`conda-forge <https://conda-forge.org/>`__ than on PyPI or CRAN.

Several test jobs on CircleCI_ related to the installation of the tool will only
run if you change the dependencies.
These will be skipped for most pull requests.

When reviewing a pull request where dependencies are added or removed, always
check that the changes have been applied in all relevant files.

.. _pull_request_checks:

Pull request checks
-------------------

To check that a pull request is up to standard, several automatic checks are
run when you make a pull request.
Read more about it in the Tests_ and Documentation_ sections.
Successful checks have a green ✓ in front, a ❌ means the check failed.

If you need help with the checks, please ask the technical reviewer of your pull
request for help.
Ask `@ESMValGroup/tech-reviewers`_ if you do not have a technical reviewer yet.

If the checks are broken because of something unrelated to the current
pull request, please check if there is an open issue that reports the problem
and create one if there is no issue yet.
You can attract the attention of the `@ESMValGroup/esmvaltool-coreteam`_ by
mentioning them in the issue if it looks like no-one is working on solving the
problem yet.
The issue needs to be fixed in a separate pull request first.
After that has been merged into the ``main`` branch and all checks are green
again on the ``main`` branch, merge it into your own branch to get the tests
to pass.

When reviewing a pull request, always make sure that all checks were successful.
If the Codacy check keeps failing, please run prospector locally.
If necessary, ask the pull request author to do the same and to address the
reported issues.
See the section on code_quality_ for more information.
Never merge a pull request with failing CircleCI or readthedocs checks.

.. _`@ESMValGroup/esmvaltool-coreteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-coreteam
.. _`@ESMValGroup/esmvaltool-developmentteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam
.. _`@ESMValGroup/tech-reviewers`: https://github.com/orgs/ESMValGroup/teams/tech-reviewers
.. _`@ESMValGroup/science-reviewers`: https://github.com/orgs/ESMValGroup/teams/science-reviewers
.. _git-repository:

***************
GitHub Workflow
***************

Basics
======

The source code of the ESMValTool is hosted on GitHub. The following description gives an overview of the typical workflow and usage for implementing new diagnostics or technical changes into the ESMValTool. For general information on Git, see e.g. the online documentation at https://www.git-scm.com/doc.

There are *two* ESMValTool GitHub repositories available:

#. The **PUBLIC** GitHub repository is open to the public. The ESMValTool is released as open-source software under the Apache License 2.0. Use of the software constitutes acceptance of this license and terms. The PUBLIC ESMValTool repository is located at https://github.com/ESMValGroup/ESMValTool

#. The **PRIVATE** GitHub repository is restricted to the ESMValTool Development Team. This repository is only accessible to ESMValTool developers that have accepted the terms of use for the ESMValTool development environment. The use of the ESMValTool software and access to the private ESMValTool GitHub repository constitutes acceptance of these terms. *When you fork or copy this repository, you must ensure that you do not copy the PRIVATE repository into an open domain!* The PRIVATE ESMValTool repository for the ESMValTool development team is located at https://github.com/ESMValGroup/ESMValTool-private

All developments can be made in either of the two repositories. The creation of *FEATURE BRANCHES* (see below), however, is restricted to registered ESMValTool developers in both repositories. We encourage all developers to join the ESMValTool development team. Please contact the :ref:`ESMValTool Core Development Team <core-team>` if you want to join the ESMValTool development team.
The PRIVATE GitHub repository offers a central protected environment for ESMValTool developers who would like to keep their contributions undisclosed (e.g., unpublished scientific work, work in progress by PhD students) while at the same time benefiting from the possibilities of collaborating with other ESMValTool developers and having a backup of their work. *FEATURE BRANCHES* created in the PRIVATE repository are only visible to the ESMValTool development team but not to the public. The concept of a PRIVATE repository has proven to be very useful to efficiently share code during the development across institutions and projects in a common repository without having the contributions immediately accessible to the public.

Both, the PUBLIC and the PRIVATE repository, contain the following kinds of branches:

* *MAIN BRANCH* (official releases),
* *DEVELOPMENT BRANCH* (includes approved new contributions but version is not yet fully tested),
* *FEATURE BRANCH* (development branches for new features and diagnostics created by developers, the naming convention for *FEATURE BRANCHES* is <Project>_<myfeature>).

Access rights
=============

* Write access to the *MAIN* and *DEVELOPMENT BRANCH* in both, the PUBLIC and the PRIVATE GitHub repositories, is restricted to the :ref:`ESMValTool Core Development Team <core-team>`.
* *FEATURE BRANCHES* in both the PUBLIC and the PRIVATE repository can be created by all members of the ESMValTool development team (i.e. members in the GitHub organization "ESMValGroup"). If needed, branches can be individually write-protected within each repository so that other developers cannot accidently push changes to these branches.

The *MAIN BRANCH* of the PRIVATE repository will be regularly synchronized with the *MAIN BRANCH* of the PUBLIC repository (currently by hand). This ensures that they are identical at all times (see schematic in Figure :numref:`fig-git`). The recommended workflow for members of the ESMValTool development team is to create additional *FEATURE BRANCHES* in either the PUBLIC or the PRIVATE repository, see further instructions below.

.. _fig-git:

.. figure:: /figures/git_diagram.png
   :width: 10cm
   :align: center

   Schematic diagram of the ESMValTool GitHub repositories.

Workflow
========

The following description gives an overview of the typical workflow and usage for implementing new diagnostics or technical changes into the ESMValTool. The description assumes that your local development machine is running a Unix-like operating system. For a general introduction to Git tutorials such as, for instance, https://www.git-scm.com/docs/gittutorial are recommended.

Getting started
---------------

First make sure that you have Git installed on your development machine. On shared machines, software is usually installed using the environment modules. Try e.g.

.. code:: bash

   module avail git

if this is the case. You can ask your system administrator for assistance. You can test this with the command:

.. code:: bash

   git --version

In order to properly identify your contributions to the ESMValTool you need to configure your local Git with some personal data. This can be done with the following commands:

.. code:: bash

   git config --global user.name "YOUR NAME"
   git config --global user.email "YOUR EMAIL"

.. note:: For working on GitHub you need to create an account and login to https://github.com/.

Working with the ESMValTool GitHub Repositories
-----------------------------------------------

As a member of the ESMValTool development team you can create *FEATURE BRANCHES* in the PUBLIC as well as in the PRIVATE repository. We encourage all ESMValTool developers to use the following workflow for long-lived developments (>2 weeks).

* Login to GitHub.com
* On GitHub, go to the website of the ESMValTool repository (https://github.com/ESMValGroup/ESMValTool-private or https://github.com/ESMValGroup/ESMValTool)
* Click on the button create *FEATURE BRANCH*
* Select the *"DEVELOPMENT" BRANCH* and create a new *FEATURE BRANCH* for the diagnostic/feature you want to implement. Please follow the following naming convention for your new *FEATURE BRANCH*: <Project>_<myfeature>.

.. figure::  /figures/git_branch.png
   :align:   center
   :width:   6cm

* Click the button "Clone or Download" and copy the URL shown there
* Open a terminal window and go to the folder where you would like to store your local copy of the ESMValTool source
* Type git clone, and paste the URL:

.. code:: bash

   git clone <URL_FROM_CLIPBOARD>

This will clone the ESMValTool repository at GitHub to a local folder. You can now query the status of your local working copy with:

.. code:: bash

   git status

You will see that you are on a branch called main and your local working copy is up to date with the remote repository. With

.. code:: bash

   git branch --all

you can list all available remote and local branches. Now switch to your feature branch by:

.. code:: bash

   git checkout <NAME_OF_YOUR_FEATURE_BRANCH>

You can now start coding. To check your current developments you can use the command

.. code:: bash

   git status

You can add new files and folders that you want to have tracked by Git using:

.. code:: bash

   git add <NEW_FILE|FOLDER>

Commit your tracked changes to your local working copy via:

.. code:: bash

   git commit -m "YOUR COMMIT MESSAGE"

You can inspect your changes with (use man git-log for all options):

.. code:: bash

   git log

To share your work and to have an online backup, push your local development to your *FEATURE BRANCH* on GitHub:

.. code:: bash

   git push origin <YOUR_FEATURE_BRANCH>

.. note:: An overview on Git commands and best practices can be found e.g. here: https://zeroturnaround.com/rebellabs/git-commands-and-best-practices-cheat-sheet/

Pull requests
-------------

Once your development is completely finished, go to the GitHub website of the ESMValTool repository and switch to your *FEATURE BRANCH*. You can then initiate a pull request by clicking on the button "New pull request". Select the *DEVELOPMENT BRANCH* as "base branch" and click on "Create pull request". Your pull request will then be tested, discussed and implemented into the *DEVELPOMENT BRANCH* by the :ref:`ESMValTool Core Development Team <core-team>`.

.. attention:: When creating a pull request, please carefully review the requirements and recommendations in CONTRIBUTING.md and try to implement those (see also checklist in the pull request template). It is recommended that you create a draft pull request early in the development process, when it is still possible to implement feedback. Do not wait until shortly before the deadline of the project you are working on. If you are unsure how to implement any of the requirements, please do not hesitate to ask for help in the pull request.

GitHub issues
-------------

In case you encounter a bug of if you have a feature request or something similar you can open an issue on the PUBLIC ESMValTool GitHub repository.

General do-s and don't-s
========================

Do-s
----

* Create a *FEATURE BRANCH* and use exclusively this branch for developing the ESMValTool. The naming convention for *FEATURE BRANCHES* is <Project>_<myfeature>.
* Comment your code as much as possible and in English.
* Use short but self-explanatory variable names (e.g., model_input and reference_input instead of xm and xr).
* Consider a modular/functional programming style. This often makes code easier to read and deletes intermediate variables immediately. If possible, separate diagnostic calculations from plotting routines.
* Consider reusing or extending existing code. General-purpose code can be found in esmvaltool/diag_scripts/shared/.
* Comment all switches and parameters including a list of all possible settings/options in the header section of your code (see also ...).
* Use templates for recipes (see ...) and diagnostics (see ...) to help with proper documentation.
* Keep your *FEATURE BRANCH* regularly synchronized with the *DEVELOPMENT BRANCH* (git merge).
* Keep developments / modifications of the ESMValTool framework / backend / basic structure separate from developments of diagnostics by creating different *FEATURE BRANCHES* for these two kinds of developments. Create *FEATURE BRANCHES* for changes / modifications of the ESMValTool framework only in the *PUBLIC* repository.

Don't-s
-------

* Do not use other programming languages than the ones currently supported (Python, R, NCL, Julia). If you are unsure what language to use, Python is probably the best choice, because it has very good libraries available and is supported by a large community. Contact the :ref:`ESMValTool Core Development Team <core-team>` if you wish to use another language, but remember that only open-source languages are supported by the ESMValTool.
* Do not develop without proper version control (see do-s above).
* Avoid large (memory, disk space) intermediate results. Delete intermediate files/variables or see modular/functional programming style.
* Do not use hard-coded pathnames or filenames.
* Do not mix developments / modifications of the ESMValTool framework and developments / modifications of diagnostics in the same *FEATURE BRANCH*.

Release schedule and procedure for ESMValCore and ESMValTool
============================================================

This document describes the process for the release of ESMValCore
and ESMValTool.
By following a defined process, we streamline the work, reduce
uncertainty about required actions, and clarify the state of the code for the
user.

ESMValTool follows a strategy of timed releases.
That means that we do releases with a regular frequency and all features
that are implemented up to a certain cut-off-point can go
into the upcoming release; those that are not are deferred to the next
release.
This means that generally no release will be delayed due to a pending feature.
Instead, the regular nature of the release guarantees that every feature can be
released in a timely manner even if a specific target release is missed.

Because of limited resources, only the latest released versions of ESMValTool and ESMValCore is maintained.
If your project requires longer maintenance or you have other concerns about
the release strategy, please contact the ESMValTool core development team, see
:ref:`contact`.


Overall Procedure
-----------------

Timeline
~~~~~~~~~

.. figure::  /figures/release-timeline.png
   :align:   center

   Example of a Release Timeline (in this case for 2.1.0)

1. Contributors assign issues (and pull requests) that they intend to finish before the due date, there is a separate milestone for ESMValCore and ESMValTool
2. The ESMValCore feature freeze takes place on the ESMValCore due date
3. Some additional testing of ESMValCore takes place
4. ESMValCore release
5. The ESMValTool feature freeze takes place
6. Some additional testing of ESMValTool takes place
7. ESMValTool release
8. Soon after the release, the core development team meets to coordinate the content of the milestone for the next release

.. _release_schedule:

Release schedule
~~~~~~~~~~~~~~~~

With the following release schedule, we strive to have three releases per year and to avoid releases too close to holidays, as well as avoiding weekends.

Upcoming releases
^^^^^^^^^^^^^^^^^

- 2.5.0 (Coordinating Release Manager: `Axel Lauer`_, team members: `Manuel Schlund`_, `Rémi Kazeroni`_)

+------------+--------------------------+
| 2022-02-07 |ESMValCore feature freeze |
+------------+--------------------------+
| 2022-02-14 |ESMValCore release        |
+------------+--------------------------+
| 2022-02-21 |ESMValTool feature freeze |
+------------+--------------------------+
| 2022-02-28 |ESMValTool release        |
+------------+--------------------------+

- 2.6.0 (Release Manager: TBD)

+------------+--------------------------+
| 2022-06-06 |ESMValCore feature freeze |
+------------+--------------------------+
| 2022-06-13 |ESMValCore release        |
+------------+--------------------------+
| 2022-06-20 |ESMValTool feature freeze |
+------------+--------------------------+
| 2022-06-27 |ESMValTool release        |
+------------+--------------------------+

Past releases
^^^^^^^^^^^^^

- 2.4.0 (Release Manager: `Klaus Zimmermann`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2021-10-04 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-10-11 | 2021-11-08 | `ESMValCore Release 2.4.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.4.0>`_ | :ref:`esmvalcore:changelog-v2-4-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-10-18 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-10-25 | 2021-11-09 | `ESMValTool Release 2.4.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.4.0>`_ |      :ref:`changelog-v2-4-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.3.1 (Bugfix, Release Manager: `Klaus Zimmermann`_)

+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|    Done    |                                            Event                                            |             Changelog              |
+============+=============================================================================================+====================================+
| 2021-07-23 | `ESMValCore Release 2.3.1 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.3.1>`_ | :ref:`esmvalcore:changelog-v2-3-1` |
+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.3.0 (Release Manager: `Klaus Zimmermann`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2021-06-07 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-06-14 | 2021-06-14 | `ESMValCore Release 2.3.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.3.0>`_ | :ref:`esmvalcore:changelog-v2-3-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-06-21 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-06-28 | 2021-07-27 | `ESMValTool Release 2.3.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.3.0>`_ |      :ref:`changelog-v2-3-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.2.0 (Release Manager: `Javier Vegas-Regidor`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2021-02-01 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-02-07 | 2021-02-09 | `ESMValCore Release 2.2.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.2.0>`_ | :ref:`esmvalcore:changelog-v2-2-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-02-14 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2021-02-21 | 2021-02-25 | `ESMValTool Release 2.2.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.2.0>`_ |      :ref:`changelog-v2-2-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.1.1 (Bugfix, Release Manager: `Valeriu Predoi`_)

+------------+---------------------------------------------------------------------------------------------+-------------------------+
|    Done    |                                            Event                                            |        Changelog        |
+============+=============================================================================================+=========================+
| 2020-12-01 | `ESMValTool Release 2.1.1 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.1.1>`_ | :ref:`changelog-v2-1-1` |
+------------+---------------------------------------------------------------------------------------------+-------------------------+

- 2.1.0 (Release Manager: `Valeriu Predoi`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2020-10-05 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-10-12 | 2020-10-12 | `ESMValCore Release 2.1.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.1.0>`_ | :ref:`esmvalcore:changelog-v2-1-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-10-19 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-10-26 | 2020-10-26 | `ESMValTool Release 2.1.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.1.0>`_ |      :ref:`changelog-v2-1-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+

- 2.0.0 (Release Manager: `Bouwe Andela`_)

+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
|  Planned   |    Done    |                                            Event                                            |             Changelog              |
+============+============+=============================================================================================+====================================+
| 2020-07-01 |            |                                  ESMValCore Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-07-20 | 2020-07-20 | `ESMValCore Release 2.0.0 <https://github.com/ESMValGroup/ESMValCore/releases/tag/v2.0.0>`_ | :ref:`esmvalcore:changelog-v2-0-0` |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-07-22 |            |                                  ESMValTool Feature Freeze                                  |                                    |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+
| 2020-08-03 | 2020-08-03 | `ESMValTool Release 2.0.0 <https://github.com/ESMValGroup/ESMValTool/releases/tag/v2.0.0>`_ |      :ref:`changelog-v2-0-0`       |
+------------+------------+---------------------------------------------------------------------------------------------+------------------------------------+



Detailed timeline steps
~~~~~~~~~~~~~~~~~~~~~~~

These are the detailed steps to take to make a release.

1. Populate the milestone

   - The core development team will make sure it adds issues that it intends to work on as early as possible.
   - Any contributor is welcome to add issues or pull requests that they intend to work on themselves to a milestone.


2. ESMValCore feature freeze

   - A release branch is created and branch protection rules are set up so only the release manager (i.e. the person in charge of the release branch) can push commits to that branch.
   - The creation of the release branch is announced to the ESMValTool development team along with the procedures to use the branch for testing and making last-minute changes (see next step)


3. Some additional testing of ESMValCore

   - Run all the recipes (optionally with a reduced amount of data) to check that they still work
   - If a bug is discovered that needs to be fixed before the release, a pull request can be made to the main branch to fix the bug. The person making the pull request can then ask the release manager to cherry-pick that commit into the release branch.


4. ESMValCore release

   - Make the release by following the :ref:`ESMValCore release instructions <esmvalcore:how-to-make-a-release>`.
   - Ask the user engagement team to announce the release to the user mailing list, the development team mailing list, on twitter


5. ESMValTool feature freeze

   - A release branch is created and branch protection rules are set up so only the release manager (i.e. the person in charge of the release branch) can push commits to that branch.
   - The creation of the release branch is announced to the ESMValTool development team along with the procedures to use the branch for testing and making last-minute changes (see next step)


6. Some additional testing of ESMValTool

   - Run all the recipes to check that they still work and ask authors to review the plots
   - If a bug is discovered that needs to be fixed before the release, a pull request can be made to the main branch to fix the bug. The person making the pull request can then ask the release manager to cherry-pick that commit into the release branch.


7. ESMValTool release

   - Make the release by following :ref:`How to make a release`
   - Ask the user engagement team to announce the release to the user mailing list, the development team mailing list, and on twitter


8. Core development team meets to coordinate the content of next milestone

   - Create a doodle for the meeting or even better, have the meeting during an ESMValTool workshop
   - Prepare the meeting by filling the milestone
   - At the meeting, discuss

     - If the proposed issues cover everything we would like to accomplish
     - Are there things we need to change about the release process
     - Who will be the release manager(s) for the next release

Bugfix releases
---------------

Next to the feature releases described above, it is also possible to have bugfix releases (2.0.1, 2.0.2, etc). In general bugfix releases will only be done on the latest release, and may include ESMValCore, ESMValTool, or both.


Procedure
~~~~~~~~~

1. One or more issues are resolved that are deemed (by the core development team) to warrant a bugfix release.
2. A release branch is created from the last release tag and the commit that fixes the bug/commits that fix the bugs are cherry-picked into it from the main branch.
3. Some additional testing of the release branch takes place.
4. The release takes place.

Compatibility between ESMValTool and ESMValCore is ensured by the appropriate version pinning of ESMValCore by ESMValTool.

Glossary
--------

Feature freeze
~~~~~~~~~~~~~~
The date on which no new features may be submitted for the upcoming release. After this date, only critical bug fixes can still be included.

Milestone
~~~~~~~~~
A milestone is a list of issues and pull-request on GitHub. It has a due date, this date is the date of the feature freeze. Adding an issue or pull request indicates the intent to finish the work on this issue before the due date of the milestone. If the due date is missed, the issue can be included in the next milestone.

Release manager
~~~~~~~~~~~~~~~
The person in charge of making the release, both technically and organizationally. Appointed for a single release.

Release branch
~~~~~~~~~~~~~~
The release branch can be used to do some additional testing before the release, while normal development work continues in the main branch. It will be branched off from the main branch after the feature freeze and will be used to make the release on the release date. The only way to still get something included in the release after the feature freeze is to ask the release manager to cherry-pick a commit from the main branch into this branch.


.. _How to make a release:

How to make an ESMValTool release
---------------------------------

The release manager makes the release, assisted by the release manager of the
previous release, or if that person is not available, another previous release
manager. Perform the steps listed below with two persons, to reduce the risk of
error.

To make a new release of the package, follow these steps:

1. Check the tests on GitHub Actions and CircleCI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check the ``nightly``
`build on CircleCI <https://circleci.com/gh/ESMValGroup/ESMValTool/tree/main>`__
and the
`GitHub Actions run <https://github.com/ESMValGroup/ESMValTool/actions>`__.
All tests should pass before making a release (branch).

2. Increase the version number
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The version number is stored in ``esmvaltool/__init__.py``,
``package/meta.yaml``, ``CITATION.cff``. Make sure to update all files.
Also update the release date in ``CITATION.cff``.
See https://semver.org for more information on choosing a version number.
Make a pull request and get it merged into ``main``.

3. Add release notes
~~~~~~~~~~~~~~~~~~~~
Use the script :ref:`draft_release_notes.py` to create create a draft of the
release notes.
This script uses the titles and labels of merged pull requests since the
previous release.
Review the results, and if anything needs changing, change it on GitHub and
re-run the script until the changelog looks acceptable.
Copy the result to the file ``doc/sphinx/source/changelog.rst``.
Make a pull request and get it merged into ``main``.

4. Create a release branch
~~~~~~~~~~~~~~~~~~~~~~~~~~
Create a branch off the ``main`` branch and push it to GitHub.
Ask someone with administrative permissions to set up branch protection rules
for it so only you and the person helping you with the release can push to it.
Announce the name of the branch in an issue and ask the members of the
`ESMValTool development team <https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam>`__
to run their favourite recipe using this branch.

5. Cherry pick bugfixes into the release branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If a bug is found and fixed (i.e. pull request merged into the
``main`` branch) during the period of testing, use the command
``git cherry-pick COMMIT_HASH``, where ``COMMIT_HASH`` is the commit hash of the
commit that needs to be cherry-picked, to include the commit for this bugfix
into the release branch.
Cherry-pick any new contributions in the order they were merged, to avoid
conflicts.
When the testing period is over, make a pull request to update
the release notes with the latest changes (do not forget to include the pull
request itself into the changelog), get it merged into ``main`` and
cherry-pick it into the release branch.

6. Make the release on GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Do a final check that all tests on CircleCI and GitHub Actions completed
successfully.
Then click the
`releases tab <https://github.com/ESMValGroup/ESMValTool/releases>`__
and create the new release from the release branch (i.e. not from ``main``).
The release tag always starts with the letter ``v`` followed by the version
number, e.g. ``v2.1.0``.

7. Create and upload the PyPI package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The package is automatically uploaded to the
`PyPI <https://pypi.org/project/ESMValTool/>`__
by a GitHub action.
If has failed for some reason, build and upload the package manually by
following the instructions below.

Follow these steps to create a new Python package:

-  Check out the tag corresponding to the release,
   e.g. ``git checkout tags/v2.1.0``
-  Make sure your current working directory is clean by checking the output
   of ``git status`` and by running ``git clean -xdf`` to remove any files
   ignored by git.
-  Install the required packages:
   ``python3 -m pip install --upgrade pep517 twine``
-  Build the package:
   ``python3 -m pep517.build --source --binary --out-dir dist/ .``
   This command should generate two files in the ``dist`` directory, e.g.
   ``ESMValTool-2.1.0-py3-none-any.whl`` and ``ESMValTool-2.1.0.tar.gz``.
-  Upload the package:
   ``python3 -m twine upload dist/*``
   You will be prompted for an API token if you have not set this up
   before, see
   `here <https://pypi.org/help/#apitoken>`__ for more information.

You can read more about this in
`Packaging Python Projects <https://packaging.python.org/tutorials/packaging-projects/>`__.

8. Update the conda-forge packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The upload to PyPI will automatically trigger an update PR on the
esmvaltool-suite-feedstock_. Check that it builds correctly and merge
the PR to update the conda-forge packages.

.. _esmvaltool-suite-feedstock: https://github.com/conda-forge/esmvaltool-suite-feedstock

Changelog
---------

- 2020-09-09 Converted to rst and added to repository (future changes tracked by git)
- 2020-09-03 Update during video conference (present: Bouwe Andela, Niels Drost, Javier Vegas, Valeriu Predoi, Klaus Zimmermann)
- 2020-07-27 Update including tidying up and Glossary by Klaus Zimmermann and Bouwe Andela
- 2020-07-23 Update to timeline format by Bouwe Andela and Klaus Zimmermann
- 2020-06-08 First draft by Klaus Zimmermann and Bouwe Andela

.. _Bouwe Andela: https://github.com/bouweandela
.. _Rémi Kazeroni: https://github.com/remi-kazeroni
.. _Axel Lauer: https://github.com/axel-lauer
.. _Valeriu Predoi: https://github.com/valeriupredoi
.. _Manuel Schlund: https://github.com/schlunma
.. _Javier Vegas-Regidor: https://github.com/jvegasbsc
.. _Klaus Zimmermann: https://github.com/zklaus
.. _upgrading:

************************************************************
Upgrading a namelist (recipe) or diagnostic to ESMValTool v2
************************************************************

This guide summarizes the main steps to be taken in order to port an ESMValTool namelist (now called **recipe**) and the corresponding diagnostic(s) from v1.0 to v2.0, hereafter also referred as the *"old"* and the *"new version"*, respectively. The new ESMValTool version is being developed in the public git branch ``main``. An identical version of this branch is maintained in the private repository as well and kept synchronized on an hourly basis.

In the following, it is assumed that the user has successfully installed ESMValTool v2 and has a rough overview of its structure (see `Technical Overview <http://www.esmvaltool.org/download/Righi_ESMValTool2-TechnicalOverview.pdf>`_).

Create a github issue
=====================

Create an issue in the public repository to keep track of your work and inform other developers. See an example `here <https://github.com/ESMValGroup/ESMValTool/issues/293>`_. Use the following title for the issue: "PORTING <recipe> into v2.0".
Do not forget to assign it to yourself.

Create your own branch
======================

Create your own branch from ``main`` for each namelist (recipe) to be ported:

.. code-block:: bash

    git checkout main
    git pull
    git checkout -b <recipe>

``main`` contains only v2.0 under the ``./esmvaltool/`` directory.

Convert xml to yml
==================

In ESMValTool v2.0, the namelist (now recipe) is written in yaml format (`Yet Another Markup Language format <http://www.yaml.org/>`_). It may be useful to activate the yaml syntax highlighting for the editor in use. This improves the readability of the recipe file and facilitates the editing, especially concerning the indentations which are essential in this format (like in python). Instructions can be easily found online, for example for `emacs <https://www.emacswiki.org/emacs/YamlMode>`_ and `vim <http://www.vim.org/scripts/script.php?script_id=739>`_.

A xml2yml converter is available in ``esmvaltool/utils/xml2yml/``, please refer to the corresponding README file for detailed instructions on how to use it.

Once the recipe is converted, a first attempt to run it can be done, possibly starting with a few datasets and one diagnostics and proceed gradually. The recipe file ``./esmvaltool/recipes/recipe_perfmetrics_CMIP5.yml`` can be used as an example, as it covers most of the common cases.

Do not forget to also rewrite the recipe header in a ``documentation`` section using the yaml syntax and, if possible, to add  themes and realms item to each diagnostic section. All keys and tags used for this part must be defined in ``./esmvaltool/config-references.yml``. See ``./esmvaltool/recipes/recipe_perfmetrics_CMIP5.yml`` for an example.

Create a copy of the diag script in v2.0
========================================

The diagnostic script to be ported goes into the directory ./esmvaltool/diag_script/. It is recommended to get a copy of the very last version of the script to be ported from the ``version1`` branch (either in the public or in the private repository). Just create a local (offline) copy of this file from the repository and add it to ../esmvaltool/diag_script/ as a new file.

Note that (in general) this is not necessary for plot scripts and for the libraries in ``./esmvaltool/diag_script/ncl/lib/``, which have already been ported. Changes may however still be necessary, especially in the plot scripts which have not yet been fully tested with all diagnostics.

Check and apply renamings
=========================

The new ESMValTool version includes a completely revised interface, handling the communication between the python workflow and the (NCL) scripts. This required several variables and functions to be renamed or removed. These changes are listed in the following table and have to be applied to the diagnostic code before starting with testing.

.. tabularcolumns:: |p{6cm}|p{6cm}|p{3cm}|

+-------------------------------------------------+-----------------------------------------------------------+------------------+
| Name in v1.0                                    | Name in v2.0                                              | Affected code    |
+=================================================+===========================================================+==================+
| ``getenv("ESMValTool_wrk_dir")``                | ``config_user_info@work_dir``                             | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``getenv(ESMValTool_att)``                      | ``diag_script_info@att`` or                               | all .ncl scripts |
|                                                 | ``config_user_info@att``                                  |                  |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``xml``                                         | ``yml``                                                   | all scripts      |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``var_attr_ref(0)``                             | ``variable_info@reference_dataset``                       | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``var_attr_ref(1)``                             | ``variable_info@alternative_dataset``                     | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``models``                                      | ``input_file_info``                                       | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``models@name``                                 | ``input_file_info@dataset``                               | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``verbosity``                                   | ``config_user_info@log_level``                            | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``isfilepresent_esmval``                        | ``fileexists``                                            | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``messaging.ncl``                               | ``logging.ncl``                                           | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``info_output(arg1, arg2, arg3)``               | ``log_info(arg1)`` if ``arg3=1``                          | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``info_output(arg1, arg2, arg3)``               | ``log_debug(arg1)`` if ``arg3>1``                         | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``verbosity = config_user_info@verbosity``      | remove this statement                                     | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``enter_msg(arg1, arg2, arg3)``                 | ``enter_msg(arg1, arg2)``                                 | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``leave_msg(arg1, arg2, arg3)``                 | ``leave_msg(arg1, arg2)``                                 | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``noop()``                                      | appropriate ``if-else`` statement                         | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``nooperation()``                               | appropriate ``if-else`` stsatement                        | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``fullpaths``                                   | ``input_file_info@filename``                              | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``get_output_dir(arg1, arg2)``                  | ``config_user_info@plot_dir``                             | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``get_work_dir``                                | ``config_user_info@work_dir``                             | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``inlist(arg1, arg2)``                          | ``any(arg1.eq.arg2)``                                     | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``load interface_scripts/*.ncl``                | ``load $diag_scripts/../interface_scripts/interface.ncl`` | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``<varname>_info.tmp``                          | ``<varname>_info.ncl`` in ``preproc`` dir                 | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``ncl.interface``                               | ``settings.ncl`` in ``run_dir`` and                       | all .ncl scripts |
|                                                 | ``interface_scripts/interface.ncl``                       |                  |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``load diag_scripts/lib/ncl/``                  | ``load $diag_scripts/shared/``                            | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``load plot_scripts/ncl/``                      | ``load $diag_scripts/shared/plot/``                       | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``load diag_scripts/lib/ncl/rgb/``              | ``load $diag_scripts/shared/plot/rgb/``                   | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``load diag_scripts/lib/ncl/styles/``           | ``load $diag_scripts/shared/plot/styles``                 | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``load diag_scripts/lib/ncl/misc_function.ncl`` | ``load $diag_scripts/shared/plot/misc_function.ncl``      | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``LW_CRE``, ``SW_CRE``                          | ``lwcre``, ``swcre``                                      | some yml recipes |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``check_min_max_models``                        | ``check_min_max_datasets``                                | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``get_ref_model_idx``                           | ``get_ref_dataset_idx``                                   | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+
| ``get_model_minus_ref``                         | ``get_dataset_minus_ref``                                 | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------------+------------------+

The following changes may also have to be considered:

- namelists are now called recipes and collected in ``esmvaltool/recipes``;
- models are now called datasets and all files have been updated accordingly, including NCL functions (see table above);
- ``run_dir`` (previous ``interface_data``), ``plot_dir``, ``work_dir`` are now unique to each diagnostic script, so it is no longer necessary to define specific paths in the diagnostic scripts to prevent file collision;
- ``input_file_info`` is now a list of a list of logicals, where each element describes one dataset and one variable. Convenience functions to extract the required elements (e.g., all datasets of a given variable) are provided in ``esmvaltool/interface_scripts/interface.ncl``;
- the interface functions ``interface_get_*`` and ``get_figure_filename`` are no longer available: their functionalities can be easily reproduced using the ``input_file_info`` and the convenience functions in ``esmvaltool/interface_scripts/interface.ncl`` to access the required attributes;
- there are now only 4 log levels (``debug``, ``info``, ``warning``, and ``error``) instead of (infinite) numerical values in ``verbosity``
- diagnostic scripts are now organized in subdirectories in ``esmvaltool/diag_scripts/``: all scripts belonging to the same diagnostics are to be collected in a single subdirectory (see ``esmvaltool/diag_scripts/perfmetrics/`` for example). This applies also to the ``aux_`` scripts, unless they are shared among multiple diagnostics (in this case they go in ``shared/``);
- the relevant input_file_info items required by a plot routine should be passed as argument to the routine itself;
- upper case characters have to be avoided in script names, if possible.

As for the recipe, the diagnostic script ``./esmvaltool/diag_scripts/perfmetrics/main.ncl`` can be followed as working example.

Move preprocessing from the diagnostic script to the backend
============================================================

Many operations previously performed by the diagnostic scripts, are now included in the backend, including level extraction, regridding, masking, and multi-model statistics. If the diagnostics to be ported contains code performing any of such operations, the corresponding code has to be removed from the diagnostic script and the respective backend functionality can be used instead.

The backend operations are fully controlled by the ``preprocessors`` section in the recipe. Here, a number of preprocessor sets can be defined, with different options for each of the operations. The sets defined in this section are applied in the ``diagnostics`` section to preprocess a given variable.

It is recommended to proceed step by step, porting and testing each operation separately before proceeding with the next one. A useful setting in the user configuration file (``config-private.yml``) called ``write_intermediary_cube`` allows writing out the variable field after each preprocessing step, thus facilitating the comparison with the old version (e.g., after CMORization, level selection, after regridding, etc.). The CMORization step of the new backend exactly corresponds to the operation performed by the old backend (and stored in the ``climo`` directory, now called ``preprec``): this is the very first step to be checked, by simply comparing the intermediary file produced by the new backend after CMORization with the output of the old backend in the ``climo`` directorsy (see "Testing" below for instructions).

The new backend also performs variable derivation, replacing the ``calculate`` function in the ``variable_defs`` scripts. If the recipe which is being ported makes use of derived variables, the corresponding calculation must be ported from the ``./variable_defs/<variable>.ncl`` file to ``./esmvaltool/preprocessor/_derive.py``.

Note that the Python library ``esmval_lib``, containing the ``ESMValProject`` class is no longer available in version 2. Most functionalities have been moved to the new preprocessor. If you miss a feature, please open an issue on github [https://github.com/ESMValGroup/ESMValTool/issues].

Move diagnostic- and variable-specific settings to the recipe
===============================================================

In the new version, all settings are centralized in the recipe, completely replacing the diagnostic-specific settings in ``./nml/cfg_files/`` (passed as ``diag_script_info`` to the diagnostic scripts) and the variable-specific settings in ``variable_defs/<variable>.ncl`` (passed as ``variable_info``). There is also no distinction anymore between diagnostic- and variable-specific settings: they are collectively defined in the ``scripts`` dictionary of each diagnostic in the recipe and passed as ``diag_script_info`` attributes by the new ESMValTool interface. Note that the ``variable_info`` logical still exists, but it is used to pass variable information as given in the corresponding dictionary of the recipe.

Make sure the diagnostic script writes NetCDF output
======================================================

Each diagnostic script is required to write the output of the anaylsis in one or more NetCDF files. This is to give the user the possibility to further look into the results, besides the plots, but (most importantly) for tagging purposes when publishing the data in a report and/or on a website.

For each of the plot produced by the diagnostic script a single NetCDF file has to be generated. The variable saved in this file should also contain all the necessary metadata that documents the plot (dataset names, units, statistical methods, etc.).
The files have to be saved in the work directory (defined in `cfg['work_dir']` and `config_user_info@work_dir`, for the python and NCL diagnostics, respectively).

Test the recipe/diagnostic in the new version
===============================================

Once complete, the porting of the diagnostic script can be tested. Most of the diagnostic script allows writing the output in a NetCDF file before calling the plotting routine. This output can be used to check whether the results of v1.0 are correctly reproduced. As a reference for v1.0, it is recommended to use the development branch.

There are two methods for comparing NetCDF files: ``cdo`` and ``ncdiff``. The first method is applied with the command:

.. code-block:: bash

    cdo diffv old_output.nc new_output.nc

which will print a log on the stdout, reporting how many records of the file differ and the absolute/relative differences.

The second method produces a NetCDF file (e.g., ``diff.nc``) with the difference between two given files:

.. code-block:: bash

    ncdiff old_output.nc new_output.nc diff.nc

This file can be opened with ``ncview`` to visually inspect the differences.

In general, binary identical results cannot be expected, due to the use of different languages and algorithms in the two versions, especially for complex operations such as regridding. However, difference within machine precision are desirable. At this stage, it is essential to test all datasets in the recipe and not just a subset of them.

It is also recommended to compare the graphical output (this may be necessary if the ported diagnostic does not produce a NetCDF output). For this comparison, the PostScript format is preferable, since it is easy to directly compare two PostScript files with the standard ``diff`` command in Linux:

.. code-block:: bash

   diff old_graphic.ps new_graphic.ps

but it is very unlikely to produce no differences, therefore visual inspection of the output may also be required.

Clean the code
==============

Before submitting a pull request, the code should be cleaned to adhere to the coding standard, which are somehow stricter in v2.0. This check is performed automatically on GitHub (CircleCI and Codacy) when opening a pull request on the public repository. A code-style checker (``nclcodestyle``) is available in the tool to check NCL scripts and installed alongside the tool itself. When checking NCL code style, the following should be considered in addition to the warning issued by the style checker:

- two-space instead of four-space indentation is now adopted for NCL as per NCL standard;
- ``load`` statements for NCL standard libraries should be removed: these are automatically loaded since NCL v6.4.0 (see `NCL documentation <http://www.ncl.ucar.edu/current_release.shtml#PreloadedScripts6.4.0>`_);
- the description of diagnostic- and variable-specific settings can be moved from the header of the diagnostic script to the recipe, since the settings are now defined there (see above);
- NCL ``print`` and ``printVarSummary`` statements must be avoided and replaced by the ``log_info`` and ``log_debug`` functions;
- for error and warning statements, the ``error_msg`` function can be used, which automatically include an exit statement.

Update the documentation
========================

If necessary, add or update the documentation for your recipes in the corrsponding rst file, which is now in ``doc\sphinx\source\recipes``. Do not forget to also add the documentation file to the list in ``doc\sphinx\source\annex_c`` to make sure it actually appears in the documentation.

Open a pull request
===================

Create a pull request on github to merge your branch back to ``main``, provide a short description of what has been done and nominate one or more reviewers.
.. _new-dataset:

Making a new dataset
********************

If you are contributing a new dataset, please have a look at :ref:`new-cmorizer` for how to do so.
If you need the new dataset for a new recipe, please make a separate pull
for the CMORizer script.

.. _dataset-documentation:

Dataset documentation
=====================

The documentation required for a CMORizer script is the following:

- Make sure that the new dataset is added to the list of
  :ref:`supported_datasets`
- The in code documentation should contain clear instructions on how to obtain
  the data
- A BibTeX file named ``<dataset>.bibtex`` defining the reference for the new
  dataset should be placed in the directory ``esmvaltool/references/``, see
  :ref:`adding_references` for detailed instructions.

For more general information on writing documentation, see :ref:`documentation`.

.. _dataset-test:

Testing
=======

When contributing a new script, add an entry for the CMORized data to
`recipes/examples/recipe_check_obs.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_check_obs.yml>`__
and run the recipe, to make sure the CMOR checks pass without warnings or errors.

To test a pull request for a new CMORizer script:

#. Download the data following the instructions included in the script and place
   it in the ``RAWOBS`` path specified in your ``config-user.yml``
#. Run the CMORizer script by running ``cmorize_obs -c <config-file> -o <dataset>``
#. Copy the resulting data to the ``OBS`` (for CMIP5 compliant data) or ``OBS6``
   (for CMIP6 compliant data) path specified in your
   ``config-user.yml``
#. Run ``recipes/examples/recipe_check_obs.yml`` with the new dataset to check that
   the data can be used

.. _dataset-sanity-check:

Scientific sanity check
=======================

When contributing a new dataset, we expect that the numbers and units of the dataset look physically meaningful.
The scientific reviewer needs to check this.

Data availability
=================

Once your pull request has been approved by the reviewers, ask
`@remi-kazeroni <https://github.com/remi-kazeroni>`_
to add the new dataset to the data pool at DKRZ and CEDA-Jasmin.
He is also the person in charge of merging CMORizer pull requests.

.. _dataset_checklist:

Detailed checklist for reviews
==============================

This (non-exhaustive) checklist provides ideas for things to check when reviewing
pull requests for new or updated CMORizer scripts.

Dataset description
-------------------

Check that new dataset has been added to the table of observations defined in
the ESMValTool guide user’s guide in section :ref:`inputdata`
(generated from ``doc/sphinx/source/input.rst``).

BibTeX info file
----------------

Check that a BibTeX file, i.e. ``<dataset>.bibtex`` defining the reference for
the new dataset has been created in ``esmvaltool/references/``.

recipe_check_obs.yml
--------------------

Check that new dataset has been added to the testing recipe
``esmvaltool/recipes/examples/recipe_check_obs.yml``

CMORizer script
---------------

Check that the new CMORizer script
``esmvaltool/cmorizers/obs/cmorize_obs_<dataset>.{py,ncl}``
meets standards.
This includes the following items:

* In-code documentation (header) contains

  1. Download instructions
  2. Reference(s)

* Code quality checks

  1. Code quality (e.g. no hardcoded pathnames)
  2. No Codacy errors reported


Config file
-----------

If present, check config file ``<dataset>.yml`` in
``esmvaltool/cmorizers/obs/cmor_config/`` for correctness.
Use ``yamllint`` to check for syntax errors and common mistakes.

Run CMORizer
------------

Make sure CMORizer is working by running ``cmorize_obs -c <config-file> -o <dataset>``

Check output of CMORizer
------------------------

After successfully running the new CMORizer, check that:

* Output contains (some) valid values (e.g. not only nan or zeros)
* Metadata is defined properly

Run ``esmvaltool/recipes/examples/recipe_check_obs.yml`` for new dataset.


RAW data
--------

Contact person in charge of ESMValTool data pool (`@remi-kazeroni`_) and
request to copy RAW data to RAWOBS/Tier2 (Tier3).


CMORized data
-------------

Contact person in charge of ESMValTool data pool (`@remi-kazeroni`_) and
request to

* Merge the pull request
* Copy CMORized dataset to OBS/Tier2 (Tier3)
* Set file access rights for new dataset
.. _private_repository:

Moving work from the private to the public repository
*****************************************************

In case you develop a new diagnostic with the ESMValTool, and you plan on publishing the results of the diagnostic in a peer-reviewed paper, you might want to develop the diagnostic in a slightly less open setting than the ESMValTool-public repository. That is what the ESMValTool-private repository is for. It would be great, though, if you would make the diagnostic available for the whole community after your paper was accepted. The steps that you need to take to develop a diagnostic in the private repository and then open a pull request for it in the public repository are described in the following:

1. Clone the private repository
===============================
For example, to clone a repository called esmvaltool-private, you would run:

``git clone git@github.com:esmvalgroup/esmvaltool-private``

or

``git clone https://github.com/esmvalgroup/esmvaltool-private``


2. Make a branch to develop your recipe and diagnostic
======================================================
``git checkout main``

``git pull``

``git checkout -b my-awesome-diagnostic``


3. Develop your diagnostic in that branch and push it to the private repository
===============================================================================
``git push -u origin my-awesome-diagnostic``

the first time and

``git push``

any other time


4. Write and submit your paper
==============================

5. Push your branch to the public repository
============================================
Add the public repository as a remote

``git remote add public git@github.com:esmvalgroup/esmvaltool``

or

``git remote add public https://github.com/esmvalgroup/esmvaltool``

and push your branch to the public repository

``git push -u public my-awesome-diagnostic``


6. Make a pull request in the public repository
===============================================
Go to https://github.com/esmalgroup/esmvaltool/pulls and click the 'New pull request button'.
Process reviewer comments and get it merged as described in :ref:`reviewing`.

7. Obtain a DOI for your code and add it to your paper
======================================================
Wait for a new release of ESMValTool. Releases are scheduled normally every four months. Find the release schedule here: :ref:`release_schedule`.
With the next release, your diagnostic recipe and source code will automatically be included in the archive on Zenodo and you can add the DOI from Zenodo to your paper: https://zenodo.org/record/3698045
.. _community:

Contributing to the community
=============================

**Contributions are very welcome!**

This chapter explains how to contribute to ESMValTool.
We greatly value contributions of any kind.
Contributions could include, but are not limited to documentation improvements,
bug reports, new or improved diagnostic code, scientific and technical code
reviews, infrastructure improvements, mailing list and chat participation,
community help/building, education and outreach.

If you have a bug or other issue to report, please open an issue on the
`issues tab on the ESMValTool github
repository <https://github.com/ESMValGroup/ESMValTool/issues>`__.

In case anything is unclear feel free to contact us for more information and
help, e.g. on our
`GitHub Discussions page <https://github.com/ESMValGroup/ESMValTool/discussions>`__.

.. toctree::
   :maxdepth: 1

		Contributing code and documentation <code_documentation>
		Contributing a diagnostic or recipe <diagnostic>
		Contributing a dataset <dataset>
		Contributing a review <review>
		Upgrading a namelist to a recipe <upgrading>
		GitHub workflow <repository>
		Moving work from the private to the public repository <private_repository>
		Release schedule and procedure <release_strategy>
.. _new-diagnostic:

Making a new diagnostic or recipe
*********************************

Getting started
===============

Please discuss your idea for a new diagnostic or recipe with the development team before getting started,
to avoid disappointment later. A good way to do this is to open an
`issue on GitHub <https://github.com/ESMValGroup/ESMValTool/issues>`_.
This is also a good way to get help.

.. _diagnostic_from_example:

Creating a recipe and diagnostic script(s)
==========================================
First create a recipe in esmvaltool/recipes to define the input data your analysis script needs
and optionally preprocessing and other settings.
Also create a script in the
`esmvaltool/diag_scripts <https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts>`_
directory and make sure it is referenced from your recipe.
The easiest way to do this is probably to copy the example recipe and diagnostic
script and adjust those to your needs.

If you have no preferred programming language yet, Python 3 is highly recommended, because it is most well supported.
However, NCL, R, and Julia scripts are also supported.

Good example recipes for the different languages are:

-  python: `esmvaltool/recipes/examples/recipe_python.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_python.yml>`_
-  R: `esmvaltool/recipes/examples/recipe_r.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_r.yml>`_
-  julia: `esmvaltool/recipes/examples/recipe_julia.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_julia.yml>`_
-  ncl: `esmvaltool/recipes/examples/recipe_ncl.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_ncl.yml>`_

Good example diagnostics are:

-  python: `esmvaltool/diag_scripts/examples/diagnostic.py <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/diag_scripts/examples/diagnostic.py>`_
-  R: `esmvaltool/diag_scripts/examples/diagnostic.R <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/diag_scripts/examples/diagnostic.R>`_
-  julia: `esmvaltool/diag_scripts/examples/diagnostic.jl <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/diag_scripts/examples/diagnostic.jl>`_
-  ncl: `esmvaltool/diag_scripts/examples/diagnostic.ncl <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/diag_scripts/examples/diagnostic.ncl>`_

For an explanation of the recipe format, you might want to read about the
:ref:`ESMValTool recipe <esmvalcore:recipe_overview>` and have a look at the
available :ref:`preprocessor functions <esmvalcore:preprocessor>`.
For further inspiration, check out the already
:ref:`available recipes and diagnostics <recipes>`.

There is a directory
`esmvaltool/diag_scripts/shared <https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/shared>`_
for code that is shared by many diagnostics.
This directory contains code for creating common plot types, generating output
file names, selecting input data, and other commonly needed functions.
See :ref:`api_shared` for the documentation of the shared Python code.

Re-using existing code
======================
Always make sure your code is or can be released under a license that is compatible with the Apache 2.0 license.

If you have existing code in a supported scripting language, you have two options for re-using it. If it is fairly
mature and a large amount of code, the preferred way is to package and publish it on the
official package repository for that language and add it as a dependency of ESMValTool.
If it is just a few simple scripts or packaging is not possible (i.e. for NCL) you can simply copy
and paste the source code into the ``esmvaltool/diag_scripts`` directory.

If you have existing code in a compiled language like
C, C++, or Fortran that you want to re-use, the recommended way to proceed is to add Python bindings and publish
the package on PyPI so it can be installed as a Python dependency. You can then call the functions it provides
using a Python diagnostic.

.. _recipe_documentation:

Recipe and diagnostic documentation
===================================

This section describes how to document a recipe.
For more general information on writing documentation, see :ref:`documentation`.

On readthedocs
--------------

Recipes should have a page in the :ref:`recipes` chapter which describes what
the recipe/diagnostic calculates.

When adding a completely new recipe, please start by copying
`doc/sphinx/source/recipes/recipe_template.rst.template <https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/recipes/recipe_template.rst.template>`_
to a new file ``doc/sphinx/source/recipes/recipe_<name of diagnostic>.rst``
and do not forget to add your recipe to the
`index <https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/recipes/index.rst>`_.

Fill all sections from the template:

- Add a brief description of the method
- Add references
- Document recipe options for the diagnostic scripts
- Fill in the list of variables required to run the recipe
- Add example images

An example image for each type of plot produced by the recipe should be added
to the documentation page to show the kind of output the recipe produces.
The '.png' files can be stored in a subdirectory specific for the recipe under
`doc/sphinx/source/recipes/figures <https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/recipes/figures>`_
and linked from the recipe documentation page.
A resolution of 150 `dpi <https://en.wikipedia.org/wiki/Dots_per_inch>`_ is
recommended for these image files, as this is high enough for the images to look
good on the documentation webpage, but not so high that the files become large.

In the recipe
-------------
Fill in the ``documentation`` section of the recipe as described in
:ref:`esmvalcore:recipe_documentation` and add a ``description`` to each
diagnostic entry.
When reviewing a recipe, check that these entries have been filled with
descriptive content.

In the diagnostic scripts
-------------------------
Functions implementing scientific formula should contain comments with
references to the source paper(s) and formula number(s).

When reviewing diagnostic code, check that formulas are implemented according
to the referenced paper(s) and/or other resources and that the computed numbers
look as expected from literature.

.. _diagnostic_output:

Diagnostic output
=================

Typically, diagnostic scripts create plots, but any other output such as e.g.
text files or tables is also possible.
Figures should be saved in the ``plot_dir``, either in both ``.pdf`` and
``.png`` format (preferred), or
respect the ``output_file_type`` specified in the
:ref:`esmvalcore:user configuration file`.
Data should be saved in the ``work_dir``, preferably as a ``.nc``
(`NetCDF <https://www.unidata.ucar.edu/software/netcdf/>`__) file, following the
`CF-Conventions <https://cfconventions.org/>`__ as much as possible.

Have a look at the :ref:`example scripts <diagnostic_from_example>` for how to
access the value of ``work_dir``, ``plot_dir``, and ``output_file_type`` from
the diagnostic script code.
More information on the interface between ESMValCore and the diagnostic script
is available :ref:`here <esmvalcore:interface_esmvalcore_diagnostic>` and
the description of the :ref:`outputdata` may also help to understand this.

If a diagnostic script creates plots, it should save the data used to create
those plots also to a NetCDF file.
If at all possible, there will be one NetCDF file for each plot the diagnostic
script creates.
There are several reasons why it is useful to have the plotted data available
in a NetCDF file:

- for interactive visualization of the recipe on a website
- for automated regression tests, e.g. checking that the numbers are still the
  same with newer versions of libraries

If the output data is prohibitively large, diagnostics authors can choose to
implement a ``write_netcdf: false`` diagnostic script option, so writing the
NetCDF files can be disabled from the recipe.

When doing a scientific review, please check that the figures and data look as
expected from the literature and that appropriate references have been added.

.. _recording-provenance:

Recording provenance
====================

When ESMValCore (the ``esmvaltool`` command) runs a recipe,
it will first find all data and run the default preprocessor steps plus any
additional preprocessing steps defined in the recipe. Next it will run the diagnostic script defined in the recipe
and finally it will store provenance information. Provenance information is stored in the
`W3C PROV XML format <https://www.w3.org/TR/prov-xml/>`_
and provided that the provenance tree is small, also plotted in an SVG file for
human inspection.
In addition to provenance information, a caption is also added to the plots.
When contributing a diagnostic, please make sure it records the provenance,
and that no warnings related to provenance are generated when running the recipe.
To allow the ESMValCore to keep track of provenance (e.g. which input files
were used to create what plots by the diagnostic script), it needs the
:ref:`esmvalcore:interface_diagnostic_esmvalcore`.

.. note::

    Provenance is recorded by the ``esmvaltool`` command provided by the
    ESMValCore package.
    No ``*_provenance.xml`` files will be generated when re-running just
    the diagnostic script with the command that is displayed on the screen
    during a recipe run, because that will only run the diagnostic script.

Provenance items provided by the recipe
---------------------------------------
Provenance tags can be added in several places in the recipe.
The :ref:`esmvalcore:recipe_documentation` section provides information about
the entire recipe.

For each diagnostic in the recipe, ESMValCore supports the following additional information:

- :code:`realms` a list of high-level modeling components
- :code:`themes` a list of themes

Please see the (installed version of the) file
`esmvaltool/config-references.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/config-references.yml>`_
for all available information on each item.

Provenance items provided by the diagnostic script
--------------------------------------------------
For each output file produced by the diagnostic script, ESMValCore supports the following additional information:

- :code:`ancestors` a list of input files used to create the plot.
- :code:`caption` a caption text for the plot

Note that the level of detail is limited, the only valid choices for ``ancestors`` are files produced by
:ref:`ancestor tasks<esmvalcore:ancestor-tasks>`.

It is also possible to add more information for the implemented diagnostics using the following items:

- :code:`authors` a list of authors
- :code:`references` a list of references, see :ref:`adding_references` below
- :code:`projects` a list of projects
- :code:`domains` a list of spatial coverage of the dataset
- :code:`plot_types` a list of plot types if the diagnostic created a plot, e.g. error bar
- :code:`statistics` a list of types of the statistic, e.g. anomaly

Arbitrarily named other items are also supported.

Please see the (installed version of the) file
`esmvaltool/config-references.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/config-references.yml>`_
for all available information on each item, see :ref:`esmvalcore:config-ref` for
an introduction.
In this file, the information is written in the form of ``key: value``.
Note that we add the keys to the diagnostics.
The keys will automatically be replaced by their values in the final provenance records.
For example, in the ``config-references.yml`` there is a category for types of the plots:

.. code-block:: console

  plot_types:
    errorbar: error bar plot

In the diagnostics, we add the key as:
:code:`plot_types: [errorbar]`
It is also possible to add custom provenance information by adding items to each category in this file.

In order to communicate with the diagnostic script, two interfaces have been defined,
which are described in the `ESMValCore documentation <https://docs.esmvaltool.org/projects/esmvalcore/en/latest/interfaces.html>`_.
Note that for Python and NCL diagnostics much more convenient methods are available than
directly reading and writing the interface files. For other languages these are not implemented (yet).

Depending on your preferred programming language for developing a diagnostic,
see the instructions and examples below on how to add provenance information:

Recording provenance in a Python diagnostic script
--------------------------------------------------
Always use :meth:`esmvaltool.diag_scripts.shared.run_diagnostic` at the end of your script:

.. code-block:: python

  if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

And make use of a :class:`esmvaltool.diag_scripts.shared.ProvenanceLogger` to log provenance:

.. code-block:: python

  with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)

The ``diagnostic_file`` can be obtained using :class:`esmvaltool.diag_scripts.shared.get_diagnostic_filename`.

The ``provenance_record`` is a dictionary of provenance items, for example:

.. code-block:: python

  provenance_record = {
        'ancestors': ancestor_files,
        'authors': [
            'andela_bouwe',
            'righi_mattia',
        ],
        'caption': caption,
        'domains': ['global'],
        'plot_types': ['zonal'],
        'references': [
            'acknow_project',
        ],
        'statistics': ['mean'],
      }

Have a look at the example Python diagnostic in
`esmvaltool/diag_scripts/examples/diagnostic.py <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/diag_scripts/examples/diagnostic.py>`_
for a complete example.

Recording provenance in an NCL diagnostic script
------------------------------------------------
Always call the ``log_provenance`` procedure after plotting from your NCL diag_script:

.. code-block:: console

  log_provenance(nc-file,plot_file,caption,statistics,domain,plottype,authors,references,input-files)

For example:

.. code-block:: console

  log_provenance(ncdf_outfile, \
                 map@outfile, \
                 "Mean of variable: " + var0, \
                 "mean", \
                 "global", \
                 "geo", \
                 (/"righi_mattia", "gottschaldt_klaus-dirk"/), \
                 (/"acknow_author"/), \
                 metadata_att_as_array(info0, "filename"))

Have a look at the example NCL diagnostic in
`esmvaltool/diag_scripts/examples/diagnostic.ncl <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/diag_scripts/examples/diagnostic.ncl>`_
for a complete example.

Recording provenance in a Julia diagnostic script
-------------------------------------------------
The provenance information is written in a ``diagnostic_provenance.yml`` that will be located in ``run_dir``.
For example a ``provenance_record`` can be stored in a yaml file as:

.. code-block:: julia

  provenance_file = string(run_dir, "/diagnostic_provenance.yml")

  open(provenance_file, "w") do io
      JSON.print(io, provenance_records, 4)
  end

The ``provenance_records`` can be defined as a dictionary of provenance items.
For example:

.. code-block:: julia

  provenance_records = Dict()

  provenance_record = Dict(
      "ancestors" => [input_file],
      "authors" => ["vonhardenberg_jost", "arnone_enrico"],
      "caption" => "Example diagnostic in Julia",
      "domains" => ["global"],
      "projects" => ["crescendo", "c3s-magic"],
      "references" => ["zhang11wcc"],
      "statistics" => ["other"],
  )

  provenance_records[output_file] = provenance_record

Have a look at the example Julia diagnostic in
`esmvaltool/diag_scripts/examples/diagnostic.jl <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/diag_scripts/examples/diagnostic.jl>`_
for a complete example.

Recording provenance in an R diagnostic script
----------------------------------------------
The provenance information is written in a ``diagnostic_provenance.yml`` that will be located in ``run_dir``.
For example a ``provenance_record`` can be stored in a yaml file as:

.. code-block:: R

  provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
  write_yaml(provenance_records, provenance_file)

The ``provenance_records`` can be defined as a list of provenance items.
For example:

.. code-block:: R

  provenance_records <- list()

  provenance_record <- list(
    ancestors = input_filenames,
    authors = list("hunter_alasdair", "perez-zanon_nuria"),
    caption = title,
    projects = list("c3s-magic"),
    statistics = list("other"),
  )

  provenance_records[[output_file]] <- provenance_record

.. _adding_references:

Adding references
=================
Recipes and diagnostic scripts can include references.
When a recipe is run, citation information is stored in `BibTeX <https://en.wikipedia.org/wiki/BibTeX>`__ format.
Follow the steps below to add a reference to a recipe (or a diagnostic):

-  make a ``tag`` that is representative of the reference entry.
   For example, ``righi15gmd`` shows the last name of the first author, year and journal abbreviation.
-  add the ``tag`` to the ``references`` section in the recipe (or the diagnostic script provenance, see recording-provenance_).
-  make a BibTeX file for the reference entry. There are some online tools to convert a doi to BibTeX format like https://doi2bib.org/
-  rename the file to the ``tag``, keep the ``.bibtex`` extension.
-  add the file to the folder ``esmvaltool/references``.

Note: the ``references`` section in ``config-references.yaml`` has been replaced by the folder ``esmvaltool/references``.

.. _testing_recipes:

Testing recipes
===============

To test a recipe, you can run it yourself on your local infrastructure or you
can ask the `@esmvalbot <https://github.com/apps/esmvalbot>`_ to run it for you.
To request a run of ``recipe_xyz.yml``, write the following comment below a pull
request:

::

   @esmvalbot Please run recipe_xyz.yml

Note that only members of the `@ESMValGroup/esmvaltool-developmentteam`_
can request runs. The memory of the `@esmvalbot`_ is limited to 16 GB and it only
has access to data available at DKRZ.

When reviewing a pull request, at the very least check that a recipes runs
without any modifications.
For a more thorough check, you might want to try out different datasets or
changing some settings if the diagnostic scripts support those.
A simple :ref:`tool <recipe_test_tool>` is available for testing recipes
with various settings.

.. _diagnostic_checklist:

Detailed checklist for reviews
==============================

This (non-exhaustive) checklist provides ideas for things to check when reviewing
pull requests for new or updated recipes and/or diagnostic scripts.

Technical reviews
-----------------

Documentation
~~~~~~~~~~~~~

Check that the scientific documentation of the new diagnostic has been added to
the user’s guide:

* A file ``doc/sphinx/source/recipes/recipe_<diagnostic>.rst`` exists
* New documentation is included in ``doc/sphinx/source/recipes/index.rst``
* Documentation follows template `doc/sphinx/source/recipes/recipe_template.rst.template`_
* Description of configuration options
* Description of variables
* Valid image files
* Resolution of image files (~150 dpi is usually enough; file size should be
  kept small)

Recipe
~~~~~~

Check yaml syntax (with ``yamllint``) and that new recipe contains:

* Documentation: description, authors, maintainer, references, projects
* Provenance tags: themes, realms

Diagnostic script
~~~~~~~~~~~~~~~~~

Check that the new diagnostic script(s) meet(s) standards.
This includes the following items:

* In-code documentation (comments, docstrings)
* Code quality (e.g. no hardcoded pathnames)
* No Codacy errors reported
* Re-use of existing functions whenever possible
* Provenance implemented

Run recipe
~~~~~~~~~~

Make sure new diagnostic(s) is working by running the ESMValTool with the recipe.

Check output of diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~

After successfully running the new recipe, check that:

* NetCDF output has been written
* Output contains (some) valid values (e.g. not only nan or zeros)
* Provenance information has been written

Check automated tests
~~~~~~~~~~~~~~~~~~~~~

Check for errors reported by automated tests

* Codacy
* CircleCI
* Documentation build

Scientific reviews
------------------

Documentation added to user’s guide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check that the scientific documentation of the new diagnostic
in ``doc/sphinx/source/recipes/recipe_<diagnostic>.rst``:

* Meets scientific documentation standard and
* Contains brief description of method
* Contains references
* Check for typos / broken text
* Documentation is complete and written in an understandable language
* References are complete

Recipe
~~~~~~

Check that new recipe contains valid:

* Documentation: description, references
* Provenance tags: themes, realms

Diagnostic script
~~~~~~~~~~~~~~~~~

Check that the new diagnostic script(s) meet(s) scientific standards.
This can include the following items:

* Clear and understandable in-code documentation including brief description of
  diagnostic
* References
* Method / equations match reference(s) given

Run recipe
~~~~~~~~~~

Make sure new diagnostic(s) is working by running the ESMValTool.

Check output of diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~

After successfully running the new recipe, check that:

* Output contains (some) valid values (e.g. not only nan or zeros)
* If applicable, check plots and compare with corresponding plots in the
  paper(s) cited


.. _`@ESMValGroup/esmvaltool-developmentteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam
.. _reviewing:

Review of pull requests
=======================

In the ESMValTool community we use pull request reviews to ensure all code and
documentation contributions are of good quality.
An introduction to code reviews can be found in `The Turing Way`_, including
`why code reviews are important`_ and advice on
`how to have constructive reviews`_.

Most pull requests will need two reviews before they can be merged.
First a technical review takes place and then a scientific review.
Once both reviewers have approved a pull request, it can be merged.
These three steps are described in more detail below.
If a pull request contains only technical changes, e.g. a pull request that
corrects some spelling errors in the documentation or a pull request that
fixes some installation problem, a scientific review is not needed.

If you are a regular contributor, please try to review a bit more than two
other pull requests for every pull request you create yourself, to make sure
that each pull request gets the attention it deserves.

.. _technical_review:

1. Technical review
-------------------

Technical reviews are done by the technical review team.
This team consists of regular contributors that have a strong interest and
experience in software engineering.

Technical reviewers use the technical checklist from the
`pull request template`_ to make sure the pull request follows the standards we
would like to uphold as a community.
The technical reviewer also keeps an eye on the design and checks that no major
design changes are made without the approval from the technical lead development
team.
If needed, the technical reviewer can help with programming questions, design
questions, and other technical issues.

The technical review team can be contacted by writing
`@ESMValGroup/tech-reviewers`_ in a comment on an issue or pull request on
GitHub.

.. _scientific_review:

2. Scientific review
--------------------

Scientific reviews are done by the scientific review team.
This team consists of contributors that have a strong interest and
experience in climate science or related domains.

Scientific reviewers use the scientific checklist from the
`pull request template`_ to make sure the pull request follows the standards we
would like to uphold as a community.

The scientific review team can be contacted by writing
`@ESMValGroup/science-reviewers`_ in a comment on an issue or pull request on
GitHub.

3. Merge
--------

Pull requests are merged by the `@ESMValGroup/esmvaltool-coreteam`_.
Specifically, pull requests containing a :ref:`CMORizer script<new-dataset>` can only be merged by
`@remi-kazeroni`_, who will then add the CMORized data to the OBS data pool at
DKRZ and CEDA-Jasmin.
The team member who does the merge first checks that both the technical and
scientific reviewer approved the pull request and that the reviews were
conducted thoroughly.
He or she looks at the list of files that were changed in the pull request and
checks that all relevant checkboxes from the checklist in the pull request
template have been added and ticked.
Finally, he or she checks that the :ref:`pull_request_checks` passed and
merges the pull request.
The person doing the merge commit edits the merge commit message so it
contains a concise and meaningful text.

Any issues that were solved by the pull request can be closed after merging.
It is always a good idea to check with the author of an issue and ask if it is
completely solved by the related pull request before closing the issue.

The core development team can be contacted by writing `@ESMValGroup/esmvaltool-coreteam`_
in a comment on an issue or pull request on GitHub.

Frequently asked questions
--------------------------

How do I request a review of my pull request?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you know a suitable reviewer, e.g. because your pull request fixes an issue
that they opened or they are otherwise interested in the work you are
contributing, you can ask them for a review by clicking the cogwheel next to
'Reviewers' on the pull request 'Conversation' tab and clicking on that person.
When changing code, it is a good idea to ask the original authors of that code
for a review.
An easy way to find out who previously worked on a particular piece of code is
to use `git blame`_.
GitHub will also suggest reviewers based on who previously worked on the files
changed in a pull request.
Every recipe has a maintainer and authors listed in the recipe, it is a good
idea to ask these people for a review.

If there is no obvious reviewer, you can attract the attention of the relevant
team of reviewers by writing to `@ESMValGroup/tech-reviewers`_ or
`@ESMValGroup/science-reviewers`_ in a comment on your pull request.
You can also label your pull request with one of the labels
`looking for technical reviewer <https://github.com/ESMValGroup/ESMValTool/labels/looking%20for%20technical%20reviewer>`_
or
`looking for scientific reviewer <https://github.com/ESMValGroup/ESMValTool/labels/looking%20for%20scientific%20reviewer>`_,
though asking people for a review directly is probably more effective.

.. _easy_review:

How do I optimize for a fast review?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When authoring a pull request, please keep in mind that it is easier and
faster to review a pull request that does not contain many changes.
Try to add one new feature per pull request and change only a few files.
For the ESMValTool repository, try to limit changes to a few hundred lines of
code and new diagnostics to not much more than a thousand lines of code.
For the ESMValCore repository, a pull request should ideally change no more
than about a hundred lines of existing code, though adding more lines for unit
tests and documentation is fine.

If you are a regular contributor, make sure you regularly review other people's
pull requests, that way they will be more inclined to return the favor by
reviewing your pull request.

How do I find a pull request to review?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please pick pull requests to review yourself based on your interest or
expertise.
We try to be self organizing, so there is no central authority that will assign
you to review anything.
People may advertise that they are looking for a reviewer by applying the label
`looking for technical reviewer <https://github.com/ESMValGroup/ESMValTool/labels/looking%20for%20technical%20reviewer>`_
or `looking for scientific reviewer <https://github.com/ESMValGroup/ESMValTool/labels/looking%20for%20scientific%20reviewer>`_.
If someone knows you have expertise on a certain topic, they might request your
review on a pull request though.
If your review is requested, please try to respond within a few days if at all
possible.
If you do not have the time to review the pull request, notify the author and
try to find a replacement reviewer.

How do I actually do a review?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To do a review, go to the pull request on GitHub, the list of all pull requests
is available here https://github.com/ESMValGroup/ESMValCore/pulls for the ESMValCore
and here https://github.com/ESMValGroup/ESMValTool/pulls for the ESMValTool, click the
pull request you would like to review.

The top comment should contain (a selection of) the checklist available in the
`pull request template`_.
If it is not there, copy the relevant items from the `pull request template`_.
Which items from the checklist are relevant, depends on which files are changed
in the pull request.
To see which files have changed, click the tab 'Files changed'.
Please make sure you are familiar with all items from the checklist by reading
the content linked from :ref:`pull_request_checklist` and check all items
that are relevant.
Checklists with some of the items to check are available:
:ref:`recipe and diagnostic checklist <diagnostic_checklist>` and
:ref:`dataset checklist <dataset_checklist>`.

In addition to the items from the checklist, good questions to start a review
with are 'Do I understand why these changes improve the tool?' (if not, ask the
author to improve the documentation contained in the pull request and/or the
description of the pull request on GitHub) and 'What could possibly go wrong if
I run this code?'.

To comment on specific lines of code or documentation, click the 'plus' icon
next to a line of code and write your comment.
When you are done reviewing, use the 'Review changes' button in the top right
corner to comment on, request changes to, or approve the pull request.

What if the author and reviewer disagree?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the author and the reviewer of a pull request have difficulty agreeing
on what needs to be done before the pull request can be approved, it is usually
both more pleasant and more efficient to schedule a meeting or co-working
session, for example using `Google meet`_ or `Jitsi meet`_.

When reviewing a pull request, try to refrain from making changes to the pull
request yourself, unless the author specifically agrees to those changes, as
this could potentially be perceived as offensive.

If talking about the pull requests in a meeting still does not resolve the
disagreement, ask a member of the `@ESMValGroup/esmvaltool-coreteam`_ for
their opinion and try to find a solution.


.. _`The Turing Way`: https://the-turing-way.netlify.app/reproducible-research/reviewing.html
.. _`why code reviews are important`: https://the-turing-way.netlify.app/reproducible-research/reviewing/reviewing-motivation.html
.. _`how to have constructive reviews`: https://the-turing-way.netlify.app/reproducible-research/reviewing/reviewing-recommend.html
.. _`@ESMValGroup/tech-reviewers`: https://github.com/orgs/ESMValGroup/teams/tech-reviewers
.. _`@ESMValGroup/science-reviewers`: https://github.com/orgs/ESMValGroup/teams/science-reviewers
.. _`@ESMValGroup/esmvaltool-coreteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-coreteam
.. _`@remi-kazeroni`: https://github.com/remi-kazeroni
.. _`pull request template`: https://raw.githubusercontent.com/ESMValGroup/ESMValTool/main/.github/pull_request_template.md
.. _`Google meet`: https://meet.google.com
.. _`Jitsi meet`: https://meet.jit.si
.. _`git blame`: https://www.freecodecamp.org/news/git-blame-explained-with-examples/
.. _recipes_smpi:

Single Model Performance Index (SMPI)
=====================================

Overview
--------

This diagnostic calculates the Single Model Performance Index (SMPI) following Reichler and Kim (2008). The SMPI (called "I\ :sup:`2`") is based on the comparison of several different climate variables (atmospheric, surface and oceanic) between climate model simulations and observations or reanalyses, and it focuses on the validation of the time-mean state of climate. For I\ :sup:`2` to be determined, the differences between the climatological mean of each model variable and observations at each of the available data grid points are calculated, and scaled to the interannual variance from the validating observations. This interannual variability is determined by performing a bootstrapping method (random selection with replacement) for the creation of a large synthetic ensemble of observational climatologies. The results are then scaled to the average error from a reference ensemble of models, and in a final step the mean over all climate variables and one model is calculated. The plot shows the I\ :sup:`2` values for each model (orange circles) and the multi-model mean (black circle), with the diameter of each circle representing the range of I\ :sup:`2` values encompassed by the 5th and 95th percentiles of the bootstrap ensemble. The I\ :sup:`2` values vary around one, with values greater than one for underperforming models, and values less than one for more accurate models.

Note: The SMPI diagnostic needs all indicated variables from all added models for exactly the same time period to be calculated correctly. If one model does not provide a specific variable, either that model cannot be added to the SMPI calculations, or the missing variable has to be removed from the diagnostics all together.

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_smpi.yml
* recipe_smpi_4cds.yml

Diagnostics are stored in diag_scripts/perfmetrics/

* main.ncl: calculates and (optionally) plots annual/seasonal cycles, zonal means, lat-lon fields and time-lat-lon fields. The calculated fields can also be plotted as difference w.r.t. a given reference dataset. main.ncl also calculates RMSD, bias and taylor metrics. Input data have to be regridded to a common grid in the preprocessor. Each plot type is created by a separated routine, as detailed below.
* cycle_zonal.ncl: calculates single model perfomance index (Reichler and Kim, 2008). It requires fields precalculated by main.ncl.
* collect.ncl: collects the metrics previously calculated by cycle_latlon.ncl and passes them to the plotting functions.

User settings
-------------

#. perfmetrics/main.ncl

   *Required settings for script*

   * plot_type: only "cycle_latlon (time, lat, lon)" and "cycle_zonal (time, plev, lat)" available for SMPI; usage is defined in the recipe and is dependent on the used variable (2D variable: cycle_latlon, 3D variable: cycle_zonal)
   * time_avg: type of time average (only "yearly" allowed for SMPI, any other settings are not supported for this diagnostic)
   * region: selected region (only "global" allowed for SMPI, any other settings are not supported for this diagnostic)
   * normalization: metric normalization ("CMIP5" for analysis of CMIP5 simulations; to be adjusted accordingly for a different CMIP phase)
   * calc_grading: calculates grading metrics (has to be set to "true" in the recipe)
   * metric: chosen grading metric(s) (if calc_grading is True; has to be set to "SMPI")
   * smpi_n_bootstrap: number of bootstrapping members used to determine uncertainties on model-reference differences (typical number of bootstrapping members: 100)

   *Required settings for variables*

   * reference_dataset: reference dataset to compare with (usually the observations).

These settings are passed to the other scripts by main.ncl, depending on the selected plot_type.

#. collect.ncl

   *Required settings for script*

   * metric: selected metric (has to be "SMPI")


Variables
---------

* hfds (ocean, monthly mean, longitude latitude time)
* hus (atmos, monthly mean, longitude latitude lev time)
* pr (atmos, monthly mean, longitude latitude time)
* psl (atmos, monthly mean, longitude latitude time)
* sic (ocean-ice, monthly mean, longitude latitude time)
* ta (atmos, monthly mean, longitude latitude lev time)
* tas (atmos, monthly mean, longitude latitude time)
* tauu (atmos, monthly mean, longitude latitude time)
* tauv (atmos, monthly mean, longitude latitude time)
* tos (ocean, monthly mean, longitude latitude time)
* ua (atmos, monthly mean, longitude latitude lev time)
* va (atmos, monthly mean, longitude latitude lev time)


Observations and reformat scripts
---------------------------------

The following list shows the currently used observational data sets for this recipe with their variable names and the reference to their respective reformat scripts in parentheses. Please note that obs4MIPs data can be used directly without any reformating. For non-obs4MIPs data see headers of cmorization scripts (in `/esmvaltool/cmorizers/obs/
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/>`_) for downloading and processing instructions.

* ERA-Interim (hfds, hus, psl, ta, tas, tauu, tauv, ua, va - esmvaltool/utils/cmorizers/obs/cmorize_obs_ERA-Interim.ncl)
* HadISST (sic, tos - reformat_scripts/obs/reformat_obs_HadISST.ncl)
* GPCP-SG (pr - obs4MIPs)

References
----------

* Reichler, T. and J. Kim, How well do coupled models simulate today's climate? Bull. Amer. Meteor. Soc., 89, 303-311, doi: 10.1175/BAMS-89-3-303, 2008.

Example plots
-------------

.. figure:: /recipes/figures/smpi/reichlerkim08bams_smpi.png
   :width: 70 %

   Performance index I\ :sup:`2` for individual models (circles). Circle sizes indicate the length of the 95% confidence intervals. The black circle indicates the I\ :sup:`2` of the multi-model mean (similar to Reichler and Kim (2008), Figure 1).
.. _recipe_climwip:

Climate model Weighting by Independence and Performance (ClimWIP)
=================================================================

Overview
--------

Projections of future climate change are often based on multi-model
ensembles of global climate models such as CMIP6. To condense the
information from these models they are often combined into
probabilistic estimates such as mean and a related uncertainty range
(such as the standard deviation). However, not all models in a given
multi-model ensemble are always equally ‘fit for purpose’ and it can
make sense to weight models based on their ability to simulate
observed quantities related to the target. In addition, multi-model
ensembles, such as CMIP can contain several models based on a very
similar code-base (sharing of components, only differences in
resolution etc.) leading to complex inter-dependencies between the
models. Adjusting for this by weighting models according to their
independence helps to adjust for this.


This recipe implements the **Climate model Weighting by Independence and Performance
(ClimWIP)** method. It is based on work by `Knutti et al. (2017) <https://doi.org/10.1002/2016GL072012>`_,
`Lorenz et al. (2018) <https://doi.org/10.1029/2017JD027992>`_,
`Brunner et al. (2019) <https://doi.org/10.1088/1748-9326/ab492f>`_,
`Merrifield et al. (2020) <https://doi.org/10.5194/esd-11-807-2020>`_,
`Brunner et al. (2020) <https://doi.org/10.5194/esd-11-995-2020>`_. Weights are
calculated based on historical model performance in several metrics (which can be
defined by the ``performance_contributions`` parameter) as well as by their independence
to all the other models in the ensemble based on their output fields in several metrics
(which can be defined by the ``independence_contributions`` parameter). These weights
can be used in subsequent evaluation scripts (some of which are implemented as part of
this diagnostic).

**Note**: this recipe is still being developed! A more comprehensive (yet older)
implementation can be found on GitHub:  https://github.com/lukasbrunner/ClimWIP


Using shapefiles for cutting scientific regions
-----------------------------------------------

To use shapefiles for selecting SREX or AR6 regions by name it is necessary to download them, e.g.,
from the sources below and reference the file using the `shapefile` parameter. This can either be a
absolute or a relative path. In the example recipes they are stored in  a subfolder `shapefiles`
in the `auxiliary_data_dir` (with is specified in the
`config-user.yml <https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#user-configuration-file>`_).

SREX regions (AR5 reference regions): http://www.ipcc-data.org/guidelines/pages/ar5_regions.html

AR6 reference regions: https://github.com/SantanderMetGroup/ATLAS/tree/v1.6/reference-regions

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * ``recipe_climwip_test_basic.yml``: Basic sample recipe using only a few models
    * ``recipe_climwip_test_performance_sigma.yml``: Advanced sample recipe for testing the perfect model test in particular
    * ``recipe_climwip_brunner2019_med.yml``: Slightly modified results for one region from `Brunner et al. (2019) <https://doi.org/10.1088/1748-9326/ab492f>`_ (to change regions see below)
    * ``recipe_climwip_brunner2020esd.yml``: Slightly modified results for `Brunner et al. (2020) <https://doi.org/10.5194/esd-11-995-2020>`_

Diagnostics are stored in esmvaltool/diag_scripts/weighting/climwip/

    * ``main.py``: Compute weights for each input dataset
    * ``calibrate_sigmas.py``: Compute the sigma values on the fly
    * ``core_functions.py``: A collection of core functions used by the scripts
    * ``io_functions.py``: A collection of input/output functions used by the scripts

Plot scripts are stored in esmvaltool/diag_scripts/weighting/

    * ``weighted_temperature_graph.py``: Show the difference between weighted and non-weighted temperature anomalies as time series.
    * ``weighted_temperature_map.py``: Show the difference between weighted and non-weighted temperature anomalies on a map.
    * ``plot_utilities.py``: A collection of functions used by the plot scripts.


User settings in recipe
-----------------------

1. Script ``main.py``

  *Required settings for script*
    * ``performance_sigma`` xor ``calibrate_performance_sigma``: If ``performance_contributions`` is given exactly one of the two
      has to be given. Otherwise they can be skipped or not set.

        * ``performance_sigma``: float setting the shape parameter for the performance weights calculation (determined offline).
        * ``calibrate_performance_sigma``: dictionary setting the performance sigma calibration. Has to contain at least the
          key-value pair specifying ``target``: ``variable_group``. Optional parameters for adjusting the calibration are not
          yet implemented. **Warning:** It is highly recommended to visually inspect the graphical output of the calibration to
          check if everything worked as intended. In case the calibration fails, the best performance sigma will still be
          indicated in the figure (see example :numref:`fig_climwip_5` below) but not automatically picked - the user can decide
          to use it anyway by setting it in the recipe (not recommenced).
    * ``independence_sigma``: float setting the shape parameter for the independence weights calculation (determined offline).
      Can be skipped or not set if ``independence_contributions`` is skipped or not set. A on-the-fly calculation of the
      independence sigma is not yet implemented
    * ``performance_contributions``: dictionary where the keys represent the variable groups to be included in the performance
      calculation. The values give the relative contribution of each group, with 0 being equivalent to not including the group.
      Can be skipped or not set then weights will be based purely on model independence (this is mutually exclusive with
      ``independence_contributions`` being skipped or not set).
    * ``independence_contributions``: dictionary where the keys represent the variable groups to be included in the independence
      calculation. The values give the relative contribution of each group, with 0 being equivalent to not including the group.
      If skipped or not set weights will be based purely on model performance (this is mutually exclusive with
      ``performance_contributions`` being skipped or not set).
    * ``combine_ensemble_members``: set to true if ensemble members of the same model should be combined during the processing
      (leads to identical weights for all ensemble members of the same model). Recommended if running with many (>10) ensemble
      members per model. If set to false, the model independence weighting will still (partly) account for the (very high)
      dependence between members of the same model. The success of this will depend on the case and the selected parameters.
      See `Merrifield et al. (2020) <https://doi.org/10.5194/esd-11-807-2020>`_ for an in-depth discussion.
    * ``obs_data``: list of project names to specify which are the observational data. The rest is assumed to be model data.

  *Required settings for variables*
  * This script takes multiple variables as input as long as they're available for all models
  * ``start_year``: provide the period for which to compute performance and independence.
  * ``end_year``: provide the period for which to compute performance and independence.
  * ``mip``: typically Amon
  * ``preprocessor``: e.g., climatological_mean
  * ``additional_datasets``: this should be ``*obs_data`` and is only needed for variables used in ``performance_contributions``.

  *Required settings for preprocessor*
    * Different combinations of preprocessor functions can be used, but the end result should always be aggregated over the time
      dimension, i.e. the input for the diagnostic script should be 2d (lat/lon).

  *Optional settings for preprocessor*
    * ``extract_region`` or ``extract_shape`` can be used to crop the input data.
    * ``extract_season`` can be used to focus on a single season.
    * different climate statistics can be used to calculate mean, (detrended) std_dev, or trend.

2. Script ``weighted_temperature_graph.py``

  *Required settings for script*
    * ``ancestors``: must include weights from previous diagnostic
    * ``weights``: the filename of the weights: 'weights.nc'
    * ``settings``: a list of plot settings: ``start_year`` (integer), ``end_year`` (integer), ``central_estimate`` ('mean' or integer between 0 and 100 giving the percentile), ``lower_bound`` (integer between 0 and 100), ``upper_bound`` (integer between 0 and 100)

  *Required settings for variables*
   * This script only takes temperature (tas) as input
   * ``start_year``: provide the period for which to plot a temperature change graph.
   * ``end_year``: provide the period for which to plot a temperature change graph.
   * ``mip``: typically Amon
   * ``preprocessor``: temperature_anomalies

  *Required settings for preprocessor*
    * Different combinations of preprocessor functions can be used, but the end result should always be aggregated over the
      latitude and longitude dimensions, i.e. the input for the diagnostic script should be 1d (time).

  *Optional settings for preprocessor*
    * Can be a global mean or focus on a point, region or shape
    * Anomalies can be calculated with respect to a custom reference period
    * Monthly, annual or seasonal average/extraction can be used

3. Script ``weighted_temperature_map.py``

   *Required settings for script*
     * ``ancestors``: must include weights from previous diagnostic
     * ``weights``: the filename of the weights: 'weights_combined.nc'

   *Optional settings for script*
     * ``model_aggregation``: how to aggregate the models: mean (default), median, integer between 0 and 100 representing a percentile
     * ``xticks``: positions to draw xticks at
     * ``yticks``: positions to draw yticks at

   *Required settings for variables*
     * This script takes temperature (tas) as input
     * ``start_year``: provide the period for which to plot a temperature change graph.
     * ``end_year``: provide the period for which to plot a temperature change graph.
     * ``mip``: typically Amon
     * ``preprocessor``: temperature_anomalies

   *Optional settings for variables*
     * A second variable is optional: temperature reference (tas_reference). If given, maps of temperature change to
       the reference are drawn, otherwise absolute temperatures are drawn.
     * tas_reference takes the same fields as tas


Updating the Brunner et al. (2019) recipe for new regions
---------------------------------------------------------

``recipe_climwip_brunner2019_med.yml`` demonstrates a very similar setup to `Brunner et al. (2019) <https://doi.org/10.1088/1748-9326/ab492f>`_
but only for one region (the Mediterranean). To calculated weights for other regions the recipe needs to be updated in two places:

.. code-block:: yaml

    extract_shape:
       shapefile: shapefiles/srex.shp
       decomposed: True
       method: contains
       crop: true
       ids:
         - 'South Europe/Mediterranean [MED:13]'

The ``ids`` field takes any valid `SREX <http://www.ipcc-data.org/guidelines/pages/ar5_regions.html>`_ region
key or any valid `AR6 <https://github.com/SantanderMetGroup/ATLAS/tree/v1.6/reference-regions>`_ region key
(depending on the shapefile). Note that this needs to be the full string here (not the abbreviation).

The sigma parameters need to be set according to the selected region. The sigma values for the regions
used in `Brunner et al. (2019) <https://doi.org/10.1088/1748-9326/ab492f>`_ can be found in table 1 of the paper.

.. code-block:: yaml

    performance_sigma: 0.546
    independence_sigma: 0.643

**Warning:** if a new region is used the sigma values should be recalculated! This can be done by commenting
out the sigma values (lines above) and commenting in the blocks defining the target of the weighting:

.. code-block:: yaml

    CLIM_future:
       short_name: tas
       start_year: 2081
       end_year: 2100
       mip: Amon
       preprocessor: region_mean

as well as

.. code-block:: yaml

    calibrate_performance_sigma:
       target: CLIM_future

In this case ClimWIP will attempt to perform an on-the-fly perfect model test to estimate the lowest
performance sigma (strongest weighting) which does not lead to overconfident weighting. **Important:**
the user should always check the test output for unusual behaviour. For most cases the performance sigma
should lie around 0.5. In cases where the perfect model test fails (no appropriate performance sigma
can be found) the test will still produce graphical output before raising a ValueError. The user can then decide
to manually set the performance sigma to the most appropriate value (based on the output) - **this is
not recommended** and should only be done with care! The perfect model test failing can be a hint for
one of the following: (1) not enough models in the ensemble for a robust distribution (normally >20
models should be used) or (2) the performance metrics used are not relevant for the target.

An on-the-fly calibration for the independence sigma is not yet implemented. For most cases we recommend to
use the same setup as in `Brunner et al. (2020) <https://doi.org/10.5194/esd-11-995-2020>`_ or
`Merrifield et al. (2020) <https://doi.org/10.5194/esd-11-807-2020>`_ (global or hemispherical
temperature and sea level pressure climatologies as metrics and independence sigma values between 0.2
and 0.5).

**Warning:** if a new region or target is used the provided metrics to establish the weights
might no longer be appropriate. Using unrelated metrics with no correlation and/or physical
relation to the target will reduce the skill of the weighting and ultimately render it useless! In
such cases the perfect model test might fail. This means the performance metrics should be updated.


Brunner et al. (2020) recipe and example independence weighting
---------------------------------------------------------------

``recipe_climwip_brunner2020esd.yml`` implements the weighting used in `Brunner et al. (2020) <https://doi.org/10.5194/esd-11-995-2020>`_. Compared to the paper there are minor differences due to two models which had to be excluded due to errors in the ESMValTool pre-processor: CAMS-CSM1-0 and MPI-ESM1-2-HR (r2) as well as the use of only one observational dataset (ERA5).

The recipe uses an additional step between pre-processor and weight calculation to calculate anomalies relative to the global mean (e.g., tas_ANOM = tas_CLIM - global_mean(tas_CLIM)). This means we do not use the absolute temperatures of a model as performance criterion but rather the horizontal temperature distribution (see `Brunner et al. 2020 <https://doi.org/10.5194/esd-11-995-2020>`_ for a discussion).

This recipe also implements a somewhat general independence weighting for CMIP6. In contrast to model performance (which should be case specific) model independence can largely be seen as only dependet on the multi-model ensemble in use but not the target variable or region. This means that the configuration used should be valid for similar subsets of CMIP6 as used in this recipe:


.. code-block:: yaml

   combine_ensemble_members: true
   independence_sigma: 0.54
   independence_contributions:
       tas_CLIM_i: 1
       psl_CLIM_i: 1

Note that this approach weights ensemble members of the same model with a 1/N independence scaling (combine_ensemble_members: true) as well as different models with an output-based independence weighting. Different approaches to handle ensemble members are discussed in `Merrifield et al. (2020) <https://doi.org/10.5194/esd-11-807-2020>`_. Note that, unlike for performance, climatologies are used for independence (i.e., the global mean is **not** removed for independence). **Warning:** Using only the independence weighting without any performance weighting might not always lead to meaningful results! The independence weighting is based on model output, which means that if a model is very different from all other models as well as the observations it will get a very high independence weight (and also total weight in absence of any performance weighting). This might not reflect the actual independence. It is therefore recommended to use weights based on both independence and performance for most cases.


Variables
---------

* pr (atmos, monthly mean, longitude latitude time)
* tas (atmos, monthly mean, longitude latitude time)
* psl (atmos, monthly mean, longitude latitude time)
* rsus, rsds, rlus, rlds, rsns, rlns (atmos, monthly mean, longitude latitude time)
* more variables can be added if available for all datasets.


Observations and reformat scripts
---------------------------------

Observation data is defined in a separate section in the recipe and may include
multiple datasets.

References
----------

* `Brunner et al. (2020) <https://doi.org/10.5194/esd-11-995-2020>`_, Earth Syst. Dynam., 11, 995-1012
* `Merrifield et al. (2020) <https://doi.org/10.5194/esd-11-807-2020>`_, Earth Syst. Dynam., 11, 807-834
* `Brunner et al. (2019) <https://doi.org/10.1088/1748-9326/ab492f>`_, Environ. Res. Lett., 14, 124010
* `Lorenz et al. (2018) <https://doi.org/10.1029/2017JD027992>`_, J. Geophys. Res.: Atmos., 9, 4509-4526
* `Knutti et al. (2017) <https://doi.org/10.1002/2016GL072012>`_, Geophys. Res. Lett., 44, 1909-1918

Example plots
-------------

.. _fig_climwip_1:
.. figure::  /recipes/figures/climwip/independence_tas.png
   :align:   center

   Distance matrix for temperature, providing the independence metric.

.. _fig_climwip_2:
.. figure::  /recipes/figures/climwip/performance_pr.png
   :align:   center

   Distance of preciptation relative to observations, providing the performance metric.

.. _fig_climwip_3:
.. figure::  /recipes/figures/climwip/weights_tas.png
   :align:   center

   Weights determined by combining independence and performance metrics for tas.

   .. _fig_climwip_4:
.. figure::  /recipes/figures/climwip/temperature_anomaly_graph.png
   :align:   center

   Interquartile range of temperature anomalies relative to 1981-2010, weighted versus non-weighted.

   .. _fig_climwip_5:
.. figure::  /recipes/figures/climwip/performance_sigma_calibration.png
   :align:   center

   Performance sigma calibration: The thick black line gives the reliability (c.f., weather forecast verification) which should
   reach at least 80%. The thick grey line gives the mean change in spread between the unweighted and weighted 80% ranges as an
   indication of the weighting strength (if it reaches 1, the weighting has no effect on uncertainty). The smallest sigma (i.e.,
   strongest weighting) which is not overconfident (reliability >= 80%) is selected. If the test fails (like in this example) the
   smallest sigma which comes closest to 80% will be indicated in the legend (but NOT automatically selected).

   .. _fig_climwip_6:
.. figure::  /recipes/figures/climwip/temperature_change_weighted_map.png
   :align:   center

   Map of weighted mean temperature change 2081-2100 relative to 1995-2014
.. _recipes_radiation_budget:

Radiation Budget
================

Overview
--------

The aim of monitoring the energy budget is to understand the (im)balance
of energy flux between the atmosphere and the surface of a model due to its
link with the hydrological cycle and climate change.

This diagnostic analyses the radiation budget by separating top-of-atmosphere
fluxes into clear-sky and cloud forcing components, and surface fluxes into
downwelling and upwelling components. Model predictions are compared against
three observational estimates, one of which (Stephens et al. 2012) includes
uncertainty estimates. When the black error bars overlap the zero line, the
model is consistent with observations according to Stephens et al. (2012).

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_radiation_budget.yml

Diagnostics are stored in esmvaltool/diag_scripts/radiation_budget/

    * radiation_budget.py: Plot the global radiation budget.
    * seasonal_radiation_budget.py: Write the global climatological seasonal radiation budget to a text file.



User settings in recipe
-----------------------

None


Variables
---------

* rss (atmos, monthly mean, longitude latitude time)
* rsdt (atmos, monthly mean, longitude latitude time)
* rsut (atmos, monthly mean, longitude latitude time)
* rsutcs (atmos, monthly mean, longitude latitude time)
* rsds (atmos, monthly mean, longitude latitude time)
* rls (atmos, monthly mean, longitude latitude time)
* rlut (atmos, monthly mean, longitude latitude time)
* rlutcs (atmos, monthly mean, longitude latitude time)
* rlds (atmos, monthly mean, longitude latitude time)
* hfss (atmos, monthly mean, longitude latitude time)
* hfls (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
instructions.*

* CERES-EBAF (rlut, rlutcs, rsut, rsutcs - obs4MIPs)
* Demory observations can be found in esmvaltool/diag_scripts/radiation_budget/Demory_et_al_2014_obs_Energy_Budget.yml and are from Figure 2 in Demory et al. (2014).
* Stephens observations can be found in esmvaltool/diag_scripts/radiation_budget/Stephens_et_al_2012_obs_Energy_Budget.yml from figure 1b in Stephens et al. (2012).


References
----------

* Demory, ME., Vidale, P.L., Roberts, M.J. et al. The role of horizontal resolution in simulating drivers of the global hydrological cycle. Clim Dyn 42, 2201–2225 (2014). https://doi.org/10.1007/s00382-013-1924-4
* Stephens, G., Li, J., Wild, M. et al. An update on Earth's energy balance in light of the latest global observations. Nature Geosci 5, 691–696 (2012). https://doi.org/10.1038/ngeo1580


Example plots
-------------

.. _fig_radiation_budget_1:
.. figure::  /recipes/figures/radiation_budget/UKESM1-0-LL.png
   :align:   center

   Radiation budget for UKESM1-0-LL
.. _recipes_seaice_drift:

Seaice drift
============

Overview
--------
This recipe allows to quantify the relationships between Arctic sea-ice drift
speed, concentration and thickness (Docquier et al., 2017). A decrease in
concentration or thickness, as observed in recent decades in the Arctic Ocean
(Kwok, 2018; Stroeve and Notz, 2018), leads to reduced sea-ice strength and
internal stress, and thus larger sea-ice drift speed (Rampal et al., 2011).
This in turn could provide higher export of sea ice out of the Arctic Basin,
resulting in lower sea-ice concentration and further thinning. Olason and
Notz (2014) investigate the relationships between Arctic sea-ice drift speed,
concentration and thickness using satellite and buoy observations.
They show that both seasonal and recent long-term changes in sea ice drift are
primarily correlated to changes in sea ice concentration and thickness.
This recipe allows to quantify these relationships in climate models.

In this recipe, four process-based metrics are computed based on the multi-year
monthly mean sea-ice drift speed, concentration and thickness, averaged over
the Central Arctic.

The first metric is the ratio between the modelled drift-concentration slope
and the observed drift-concentration slope. The second metric is similar to the
first one, except that sea-ice thickness is involved instead of sea-ice
concentration. The third metric is the normalised distance between the model
and observations in the drift-concentration space. The fourth metric is similar
to the third one, except that sea-ice thickness is involved instead of sea-ice
concentration.

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_seaice_drift.yml


Diagnostics are stored in diag_scripts/seaice_drift/

    * seaice_drift.py: Compute metrics and plot results


User settings in recipe
-----------------------

#. Script diag_shapeselect.py

   *Required settings (scripts)*

    One of the following two combinations is required:

    1. Latitude threshold:

        * latitude_threshold: metric will be computed north of this latitude value

    2. Polygon:

        * polygon: metric will be computed inside the give polygon. Polygon is defined as a list of (lon, lat) tuple

        * polygon_name: name of the region defined by the polygon


Variables
---------

* sispeed, sithick, siconc (daily)

Example plots
-------------

.. _fig_seaice_drift:
.. figure::  /recipes/figures/seaice_drift/drift-strength.png
   :align:   center

Scatter plots of modelled (red) and observed (blue) monthly mean
sea-ice drift speed against sea-ice concentration (left panel) and sea-ice
thickness (right panel) temporally averaged over the period 1979–2005 and
spatially averaged over the SCICEX box.
.. _recipes_impact:

Quick insights for climate impact researchers
=============================================

Overview
--------

Many impact researchers do not have the time and finances to use a large
ensemble of climate model runs for their impact analysis. To get an idea of the
range of impacts of climate change it also suffices to use a small number of
climate model runs. In case a system is only sensitive to annual temperature,
one can select a run with a high change and one with a low change of annual
temperature, preferably both with a low bias.

This recipe calculates the bias with respect to observations, and the change
with respect to a reference period, for a wide range of (CMIP) models. These
metrics are tabulated and also visualized in a diagram.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_impact.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * impact/bias_and_change.py: tabulate and visualize bias and change.


User settings in recipe
-----------------------

#. Script ``impact.py``

   *Required settings for variables*

   * tag: ``'model'`` or ``'observations'``, so the diagnostic script knows which datasets to use for the bias calculation. This must be specified for each dataset.

   *Optional settings for preprocessor*

   * Region and time settings (both for the future and reference period) can be changed at will.


Variables
---------

* tas (atmos, mon, longitude latitude time)
* pr (atmos, mon, longitude latitude time)
* any other variables of interest


Observations and reformat scripts
---------------------------------

* ERA5 data can be used via the native6 project.

References
----------

* None

Example plots
-------------

.. _fig_impact_1:
.. figure::  /recipes/figures/impact/bias_vs_change.png
   :align:   center

   "Bias and change for each variable"

.. raw:: html
    <embed>

        <style>#T_caeac_ .index_name{text-align:right}#T_caeac_
        .row_heading{text-align:right}#T_caeac_ td{padding:3px
        25px}#T_caeac_row0_col0,#T_caeac_row16_col1,#T_caeac_row16_col3{background-color:#bd1726;color:#f1f1f1}#T_caeac_row0_col1{background-color:#fffdbc;color:#000}#T_caeac_row0_col2,#T_caeac_row9_col3{background-color:#f7844e;color:#000}#T_caeac_row0_col3,#T_caeac_row10_col1{background-color:#e14430;color:#f1f1f1}#T_caeac_row1_col0{background-color:#e95538;color:#000}#T_caeac_row1_col1{background-color:#fdc171;color:#000}#T_caeac_row1_col2{background-color:#f99153;color:#000}#T_caeac_row15_col0,#T_caeac_row1_col3{background-color:#e24731;color:#f1f1f1}#T_caeac_row2_col0{background-color:#dc3b2c;color:#f1f1f1}#T_caeac_row10_col3,#T_caeac_row2_col1{background-color:#ea5739;color:#000}#T_caeac_row2_col2{background-color:#fed07e;color:#000}#T_caeac_row2_col3{background-color:#f88c51;color:#000}#T_caeac_row10_col0,#T_caeac_row3_col0{background-color:#b10b26;color:#f1f1f1}#T_caeac_row3_col1,#T_caeac_row8_col2,#T_caeac_row8_col3,#T_caeac_row9_col0{background-color:#feffbe;color:#000}#T_caeac_row3_col2{background-color:#fdb163;color:#000}#T_caeac_row11_col1,#T_caeac_row3_col3{background-color:#d93429;color:#f1f1f1}#T_caeac_row4_col0{background-color:#ab0626;color:#f1f1f1}#T_caeac_row4_col1,#T_caeac_row7_col2{background-color:#f67c4a;color:#000}#T_caeac_row4_col2{background-color:#fca55d;color:#000}#T_caeac_row15_col1,#T_caeac_row4_col3{background-color:#fa9857;color:#000}#T_caeac_row13_col3,#T_caeac_row5_col0{background-color:#ed5f3c;color:#000}#T_caeac_row5_col1{background-color:#dd3d2d;color:#f1f1f1}#T_caeac_row5_col2{background-color:#fdb365;color:#000}#T_caeac_row12_col1,#T_caeac_row5_col3,#T_caeac_row6_col3{background-color:#f67a49;color:#000}#T_caeac_row11_col2,#T_caeac_row6_col0{background-color:#f47044;color:#000}#T_caeac_row6_col1{background-color:#fba35c;color:#000}#T_caeac_row14_col0,#T_caeac_row17_col1,#T_caeac_row17_col3,#T_caeac_row6_col2,#T_caeac_row8_col0{background-color:#a50026;color:#f1f1f1}#T_caeac_row13_col0,#T_caeac_row7_col0,#T_caeac_row8_col1{background-color:#ad0826;color:#f1f1f1}#T_caeac_row16_col0,#T_caeac_row7_col1{background-color:#b50f26;color:#f1f1f1}#T_caeac_row7_col3{background-color:#d62f27;color:#f1f1f1}#T_caeac_row9_col1{background-color:#f57748;color:#000}#T_caeac_row9_col2{background-color:#ec5c3b;color:#000}#T_caeac_row10_col2{background-color:#e75337;color:#000}#T_caeac_row11_col0{background-color:#e54e35;color:#000}#T_caeac_row11_col3,#T_caeac_row13_col2{background-color:#d22b27;color:#f1f1f1}#T_caeac_row12_col0{background-color:#af0926;color:#f1f1f1}#T_caeac_row12_col2{background-color:#d42d27;color:#f1f1f1}#T_caeac_row12_col3{background-color:#e65036;color:#000}#T_caeac_row13_col1{background-color:#f16640;color:#000}#T_caeac_row14_col1,#T_caeac_row16_col2{background-color:#c62027;color:#f1f1f1}#T_caeac_row14_col2,#T_caeac_row14_col3{background-color:#f7814c;color:#000}#T_caeac_row15_col2{background-color:#fff3ac;color:#000}#T_caeac_row15_col3{background-color:#f36b42;color:#000}#T_caeac_row17_col0{background-color:#a70226;color:#f1f1f1}#T_caeac_row17_col2{background-color:#ca2427;color:#f1f1f1}</style><table
        id=T_caeac_><thead><tr><th class="level0 index_name">metric<th
        class="level0 col_heading col0"colspan=2>Bias (RMSD of all
        gridpoints)<th class="level0 col_heading col2"colspan=2>Mean change
        (Future - Reference)<tr><th class="level1 index_name">variable<th
        class="col0 col_heading level1">Temperature (K)<th class="col1
        col_heading level1">Precipitation (kg/m2/s)<th class="col2 col_heading
        level1">Temperature (K)<th class="col3 col_heading level1">Precipitation
        (kg/m2/s)<tr><th class="level0 index_name">dataset<th class=blank><th
        class=blank><th class=blank><th class=blank><tbody><tr><th class="level0
        row_heading row0"id=T_caeac_level0_row0>CMIP5_ACCESS1-0_r1i1p1<td
        class="data col0 row0"id=T_caeac_row0_col0>3.19e+00<td class="data col1
        row0"id=T_caeac_row0_col1>1.96e-05<td class="data col2
        row0"id=T_caeac_row0_col2>2.36e+00<td class="data col3
        row0"id=T_caeac_row0_col3>8.00e-09<tr><th class="level0 row_heading
        row1"id=T_caeac_level0_row1>CMIP5_BNU-ESM_r1i1p1<td class="data col0
        row1"id=T_caeac_row1_col0>4.08e+00<td class="data col1
        row1"id=T_caeac_row1_col1>1.87e-05<td class="data col2
        row1"id=T_caeac_row1_col2>2.44e+00<td class="data col3
        row1"id=T_caeac_row1_col3>2.96e-08<tr><th class="level0 row_heading
        row2"id=T_caeac_level0_row2>CMIP6_ACCESS-CM2_r1i1p1f1<td class="data
        col0 row2"id=T_caeac_row2_col0>3.75e+00<td class="data col1
        row2"id=T_caeac_row2_col1>1.77e-05<td class="data col2
        row2"id=T_caeac_row2_col2>2.87e+00<td class="data col3
        row2"id=T_caeac_row2_col3>6.63e-07<tr><th class="level0 row_heading
        row3"id=T_caeac_level0_row3>CMIP6_ACCESS-ESM1-5_r1i1p1f1<td class="data
        col0 row3"id=T_caeac_row3_col0>3.01e+00<td class="data col1
        row3"id=T_caeac_row3_col1>1.96e-05<td class="data col2
        row3"id=T_caeac_row3_col2>2.63e+00<td class="data col3
        row3"id=T_caeac_row3_col3>-1.39e-07<tr><th class="level0 row_heading
        row4"id=T_caeac_level0_row4>CMIP6_AWI-CM-1-1-MR_r1i1p1f1<td class="data
        col0 row4"id=T_caeac_row4_col0>2.91e+00<td class="data col1
        row4"id=T_caeac_row4_col1>1.80e-05<td class="data col2
        row4"id=T_caeac_row4_col2>2.56e+00<td class="data col3
        row4"id=T_caeac_row4_col3>7.67e-07<tr><th class="level0 row_heading
        row5"id=T_caeac_level0_row5>CMIP6_BCC-CSM2-MR_r1i1p1f1<td class="data
        col0 row5"id=T_caeac_row5_col0>4.22e+00<td class="data col1
        row5"id=T_caeac_row5_col1>1.74e-05<td class="data col2
        row5"id=T_caeac_row5_col2>2.64e+00<td class="data col3
        row5"id=T_caeac_row5_col3>5.02e-07<tr><th class="level0 row_heading
        row6"id=T_caeac_level0_row6>CMIP6_CAMS-CSM1-0_r1i1p1f1<td class="data
        col0 row6"id=T_caeac_row6_col0>4.43e+00<td class="data col1
        row6"id=T_caeac_row6_col1>1.84e-05<td class="data col2
        row6"id=T_caeac_row6_col2>1.48e+00<td class="data col3
        row6"id=T_caeac_row6_col3>4.89e-07<tr><th class="level0 row_heading
        row7"id=T_caeac_level0_row7>CMIP6_CESM2-WACCM_r1i1p1f1<td class="data
        col0 row7"id=T_caeac_row7_col0>2.95e+00<td class="data col1
        row7"id=T_caeac_row7_col1>1.69e-05<td class="data col2
        row7"id=T_caeac_row7_col2>2.33e+00<td class="data col3
        row7"id=T_caeac_row7_col3>-1.91e-07<tr><th class="level0 row_heading
        row8"id=T_caeac_level0_row8>CMIP6_CanESM5_r1i1p1f1<td class="data col0
        row8"id=T_caeac_row8_col0>2.81e+00<td class="data col1
        row8"id=T_caeac_row8_col1>1.69e-05<td class="data col2
        row8"id=T_caeac_row8_col2>3.36e+00<td class="data col3
        row8"id=T_caeac_row8_col3>2.10e-06<tr><th class="level0 row_heading
        row9"id=T_caeac_level0_row9>CMIP6_FGOALS-g3_r1i1p1f1<td class="data col0
        row9"id=T_caeac_row9_col0>6.74e+00<td class="data col1
        row9"id=T_caeac_row9_col1>1.80e-05<td class="data col2
        row9"id=T_caeac_row9_col2>2.13e+00<td class="data col3
        row9"id=T_caeac_row9_col3>5.95e-07<tr><th class="level0 row_heading
        row10"id=T_caeac_level0_row10>CMIP6_FIO-ESM-2-0_r1i1p1f1<td class="data
        col0 row10"id=T_caeac_row10_col0>3.02e+00<td class="data col1
        row10"id=T_caeac_row10_col1>1.75e-05<td class="data col2
        row10"id=T_caeac_row10_col2>2.07e+00<td class="data col3
        row10"id=T_caeac_row10_col3>1.89e-07<tr><th class="level0 row_heading
        row11"id=T_caeac_level0_row11>CMIP6_MIROC6_r1i1p1f1<td class="data col0
        row11"id=T_caeac_row11_col0>4.00e+00<td class="data col1
        row11"id=T_caeac_row11_col1>1.74e-05<td class="data col2
        row11"id=T_caeac_row11_col2>2.25e+00<td class="data col3
        row11"id=T_caeac_row11_col3>-2.45e-07<tr><th class="level0 row_heading
        row12"id=T_caeac_level0_row12>CMIP6_MPI-ESM1-2-HR_r1i1p1f1<td
        class="data col0 row12"id=T_caeac_row12_col0>2.98e+00<td class="data
        col1 row12"id=T_caeac_row12_col1>1.80e-05<td class="data col2
        row12"id=T_caeac_row12_col2>1.84e+00<td class="data col3
        row12"id=T_caeac_row12_col3>1.18e-07<tr><th class="level0 row_heading
        row13"id=T_caeac_level0_row13>CMIP6_MPI-ESM1-2-LR_r1i1p1f1<td
        class="data col0 row13"id=T_caeac_row13_col0>2.95e+00<td class="data
        col1 row13"id=T_caeac_row13_col1>1.78e-05<td class="data col2
        row13"id=T_caeac_row13_col2>1.82e+00<td class="data col3
        row13"id=T_caeac_row13_col3>2.52e-07<tr><th class="level0 row_heading
        row14"id=T_caeac_level0_row14>CMIP6_MRI-ESM2-0_r1i1p1f1<td class="data
        col0 row14"id=T_caeac_row14_col0>2.81e+00<td class="data col1
        row14"id=T_caeac_row14_col1>1.71e-05<td class="data col2
        row14"id=T_caeac_row14_col2>2.36e+00<td class="data col3
        row14"id=T_caeac_row14_col3>5.75e-07<tr><th class="level0 row_heading
        row15"id=T_caeac_level0_row15>CMIP6_NESM3_r1i1p1f1<td class="data col0
        row15"id=T_caeac_row15_col0>3.90e+00<td class="data col1
        row15"id=T_caeac_row15_col1>1.83e-05<td class="data col2
        row15"id=T_caeac_row15_col2>3.22e+00<td class="data col3
        row15"id=T_caeac_row15_col3>3.60e-07<tr><th class="level0 row_heading
        row16"id=T_caeac_level0_row16>CMIP6_NorESM2-LM_r1i1p1f1<td class="data
        col0 row16"id=T_caeac_row16_col0>3.08e+00<td class="data col1
        row16"id=T_caeac_row16_col1>1.70e-05<td class="data col2
        row16"id=T_caeac_row16_col2>1.74e+00<td class="data col3
        row16"id=T_caeac_row16_col3>-4.97e-07<tr><th class="level0 row_heading
        row17"id=T_caeac_level0_row17>CMIP6_NorESM2-MM_r1i1p1f1<td class="data
        col0 row17"id=T_caeac_row17_col0>2.86e+00<td class="data col1
        row17"id=T_caeac_row17_col1>1.67e-05<td class="data col2
        row17"id=T_caeac_row17_col2>1.76e+00<td class="data col3
        row17"id=T_caeac_row17_col3>-7.65e-07</table>

    </embed>
.. _recipes_hydro_forcing:

Hydro forcing comparison
========================

Overview
--------

This recipe can be used to assess the agreement between forcing datasets
(i.e. MSWEP, ERA5, ERA-Interim) for a defined catchment. The recipe can be used
to:

1. Plot a timeseries of the raw daily data
2. Plot monthly aggregrated data over a defined period
3. Plot the monthly / daily climatology statistics over a defined period


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/hydrology

    * ``recipe_hydro_forcing.yml``

Diagnostics are stored in esmvaltool/diag_scripts/hydrology/

    * ``hydro_forcing.py``: Compares and plots precipitation for MSWEP / ERA5 / ERA-5 interim datasets


User settings in recipe
-----------------------

All hydrological recipes require a shapefile as an input to select forcing data. This shapefile determines the shape of the basin for which the data will be cut out and processed. All recipes are tested with `the shapefiles <https://github.com/eWaterCycle/recipes_auxiliary_datasets/tree/main/>`_ from HydroSHEDS that are used for the eWaterCycle project. In principle any shapefile can be used, for example, the freely available basin shapefiles from the `HydroSHEDS project <https://www.hydrosheds.org/>`_.

#. recipe ``hydrology/hydro_forcing.yml``

  *Optional preprocessor settings:*

    * ``extract_shape``: The region specified here should match the catchment

  *Required settings for script:*

    * ``plot_type``: Define which plot function to run. Choices:

      * ``timeseries``: Plot a timeseries for the variable data over the defined period
      * ``climatology``: Plot the climate statistics over the defined period

  *Required settings for ``timeseries`` plots:*

    * ``time_period``: Defines the period of the output for the correct captions/labels. This value should match the period used for the preprocessor. Choices: ``day``, ``month``.



Variables
---------

* pr (atmos, daily or monthly, longitude, latitude, time)


Observations
------------

All data can be used directly without any preprocessing.

*  ERA-Interim
*  ERA5
*  MSWEP

.. References
.. ----------

.. * xxx

Example plots
-------------

.. _fig_hydro_forcing_1:
.. figure::  /recipes/figures/hydrology/Precipitation_day_plot.png
  :align:   center

  Precipitation per day for 2015-01-01:2016-12-31.

.. _fig_hydro_forcing_2:
.. figure::  /recipes/figures/hydrology/Precipitation_month_plot.png
  :align:   center

  Precipitation per month for 2015-01:2016-12.

.. _fig_hydro_forcing_3:
.. figure::  /recipes/figures/hydrology/Precipitation_climatology_month_number_plot.png
  :align:   center

  Precipitation climatology statistics per month number.

.. _fig_hydro_forcing_4:
.. figure::  /recipes/figures/hydrology/Precipitation_climatology_day_of_year_plot.png
  :align:   center

  Precipitation climatology statistics per day of year.
.. _recipes_quantilebias:

Precipitation quantile bias
===========================


Overview
--------

Precipitation is a dominant component of the hydrological cycle, and as such a main driver of the climate system and human development. The reliability of climate projections and water resources strategies therefore depends on how well precipitation can be reproduced by the models used for simulations. While global circulation models from the CMIP5 project observations can reproduce the main patterns of mean precipitation, they often show shortages and biases in the ability to reproduce the strong precipitation tails of the distribution. Most models underestimate precipitation over arid regions and overestimate it over regions of complex topography, and these shortages are amplified at high quantile precipitation. The quantilebias recipe implements calculation of the quantile bias to allow evaluation of the precipitation bias based on a user defined quantile in models as compared to a reference dataset following Mehran et al. (2014). The quantile bias (QB) is defined as the ratio of monthly precipitation amounts in each simulation to that of the reference dataset (GPCP observations in the example) above a speciﬁed threshold t (e.g., the 75th percentile of all the local monthly values). A quantile bias equal to 1 indicates no bias in the simulations, whereas a value above (below) 1 corresponds to a climate model's overestimation (underestimation) of the precipitation amount above the specified threshold t, with respect to that of the reference dataset.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_quantilebias.yml

Diagnostics are stored in diag_scripts/quantilebias/

* quantilebias.R


User settings
-------------

*Required settings for script*

* perc_lev: quantile (in %), e.g. 50


Variables
---------

* pr (atmos, monthly, longitude latitude time)


Observations and reformat scripts
---------------------------------

* GPCP-SG observations (accessible via the obs4MIPs project)


References
----------

* Mehran, A. et al.: Journal of Geophysical Research: Atmospheres, Volume 119, Issue 4, pp. 1695-1707, 2014.

Example plots
-------------

.. figure:: /recipes/figures/quantilebias/quantilebias.png
   :width: 10cm

   Quantile bias, as defined in Mehran et al. 2014, with threshold t=75th percentile, evaluated for the CanESM2 model over the 1979-2005 period, adopting GPCP-SG v 2.3 gridded precipitation as a reference dataset. The optimal reference value is 1. Both datasets have been regridded onto a 2° regular grid.
.. _recipes_martin18grl:

Drought characteristics following Martin (2018)
===============================================

Overview
--------


Following `Martin (2018)`_ drought characteristics are calculated based on the standard precipitation index (SPI), see `Mckee et al. (1993)`_. These characteristics are frequency, average duration, SPI index and severity index of drought events.

.. _`Martin (2018)`: https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018GL079807
.. _`Mckee et al. (1993)`: https://www.nature.com/articles/nclimate3387


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_martin18grl.yml


Diagnostics are stored in diag_scripts/

   * droughtindex/diag_save_spi.R
   * droughtindex/collect_drought_obs_multi.py
   * droughtindex/collect_drought_model.py
   * droughtindex/collect_drought_func.py


User settings in recipe
-----------------------

The recipe can be run with different CMIP5 and CMIP6 models and one observational or reanalysis data set.

The droughtindex/diag_save_spi.R script calculates the SPI index for any given time series. It is based on droughtindex/diag_spi.R but saves the SPI index and does not plot the histogram. The distribution and the representative time scale (smooth_month) can be set by the user, the values used in Martin (2018) are smooth_month: 6 and distribution: 'Gamma' for SPI.

There are two python diagnostics, which can use the SPI data to calculate the drought characteristics (frequency, average duration, SPI index and severity index of drought events) based on Martin (2018):

* To compare these characteristics between model data and observations or renanalysis data use droughtindex/collect_drought_obs_multi.py
  Here, the user can set:
  * indexname: Necessary to identify data produced by droughtindex/diag_save_spi.R as well as write captions and filenames. At the moment only indexname: 'SPI' is supported.
  * threshold: Threshold for this index below which an event is considered to be a drought, the setting for SPI should be usually threshold: -2.0 but any other value will be accepted. Values should not be < - 3.0 or > 3.0 for SPI (else it will identify none/always drought conditions).

* To compare these ccharacteristics between different time periods in model data use droughtindex/collect_drought_model.py
  Here, the user can set:
  * indexname: Necessary to identify data produced by droughtindex/diag_save_spi.R as well as write captions and filenames. At the moment only indexname: 'SPI' is supported.
  * threshold: Threshold for this index below which an event is considered to be a drought, the setting for SPI should be usually threshold: -2.0 but any other value will be accepted. Values should not be < - 3.0 or > 3.0 for SPI (else it will identify none/always drought conditions).
  * start_year: Needs to be equal or larger than the start_year for droughtindex/diag_save_spi.R.
  * end_year: Needs to be equal or smaller than the end_year for droughtindex/diag_save_spi.R.
  * comparison_period: should be < (end_year - start_year)/2 to have non overlapping time series in the comparison.

The third diagnostic droughtindex/collect_drought_func.py contains functions both ones above use.

Variables
---------

* *pr* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Martin, E.R. (2018). Future Projections of Global Pluvial and Drought Event Characteristics. Geophysical Research Letters, 45, 11913-11920.

* McKee, T. B., Doesken, N. J., & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. In Proceedings of the 8th Conference on Applied Climatology (Vol. 17, No. 22, pp. 179-183). Boston, MA: American Meteorological Society.

Example plots
-------------

.. _martin18grl_fig1:
.. figure:: /recipes/figures/droughtindex/martin18grl_fig1.png
   :align: center
   :width: 50%

   Global map of the percentage difference between multi-model mean of 15 CMIP models and the CRU data for the number of drought events [%] based on SPI.

.. _martin18grl_fig2:
.. figure:: /recipes/figures/droughtindex/martin18grl_fig2.png
   :align: center
   :width: 50%

   Global map of the percentage difference between multi-model mean for RCP8.5 scenarios (2050-2100) runs and historical data (1950-2000) for 15 CMIP models for the number of drought events [%] based on SPI.


.. _recipes_miles:

Blocking metrics and indices, teleconnections and weather regimes (MiLES)
=========================================================================


Overview
--------

Atmospheric blocking is a recurrent mid-latitude weather pattern identified by a large-amplitude, quasi-stationary, long-lasting, high-pressure anomaly that ‘‘blocks’’ the westerly flow forcing the jet stream to split or meander
`(Rex, 1950) <https://onlinelibrary.wiley.com/action/showCitFormats?doi=10.1111%2Fj.2153-3490.1950.tb00331.x>`_.

It is typically initiated by the breaking of a Rossby wave in a diffluence region at the exit of the storm track, where it amplifies the underlying stationary ridge `(Tibaldi and Molteni, 1990) <https://doi.org/10.1034/j.1600-0870.1990.t01-2-00003.x>`_.
Blocking occurs more frequently in the Northern Hemisphere cold season, with larger frequencies observed over the Euro-Atlantic and North Pacific sectors. Its lifetime oscillates from a few days up to several weeks `(Davini et al., 2012)  <https://doi.org/10.1175/JCLI-D-12-00032.1)>`_ sometimes leading to winter cold spells or summer heat waves.

To this end, the MId-Latitude Evaluation System (MiLES) was developed as stand-alone package (https://github.com/oloapinivad/MiLES) to support analysis of mid-latitude weather patterns in terms of atmospheric blocking, teleconnections and weather regimes. The package was then implemented as recipe for ESMValTool.

The tool works on daily 500hPa geopotential height data (with data interpolated on a common 2.5x2.5 grid) and calculates the following diagnostics:

1D Atmospheric Blocking
***********************
`Tibaldi and Molteni (1990) <https://doi.org/10.1034/j.1600-0870.1990.t01-2-00003.x>`_ index for Northern Hemisphere. Computed at fixed latitude of 60N, with delta of -5,-2.5,0,2.5,5 deg, fiN=80N and fiS=40N. Full timeseries and climatologies are provided in NetCDF4 Zip format.

2D Atmospheric blocking
***********************
Following the index by `Davini et al. (2012) <https://doi.org/10.1175/JCLI-D-12-00032.1>`_. It is a 2D version of `Tibaldi and Molteni (1990) <https://doi.org/10.1034/j.1600-0870.1990.t01-2-00003.x>`_ for Northern Hemisphere atmospheric blocking evaluating meridional gradient reversal at 500hPa. It computes both Instantaneous Blocking and Blocking Events frequency, where the latter allows the estimation of the each blocking duration. It includes also two blocking intensity indices, i.e. the Meridional Gradient Index and the Blocking Intensity index. In addition the orientation (i.e. cyclonic or anticyclonic) of the Rossby wave breaking is computed. A supplementary Instantaneous Blocking index with the GHGS2 condition (see `Davini et al., 2012 <https://doi.org/10.1175/JCLI-D-12-00032.1>`_) is also evaluated.
Full timeseries and climatologies are provided in NetCDF4 Zip format.

Z500 Empirical Orthogonal Functions
***********************************
Based on SVD. The first 4 EOFs for North Atlantic (over the 90W-40E 20N-85N box) and Northern Hemisphere (20N-85N) or a custom region are computed. North Atlantic Oscillation, East Atlantic Pattern, and Arctic Oscillation can be evaluated.
Figures showing linear regression of PCs on monthly Z500 are provided. PCs and eigenvectors, as well as the variances explained are provided in NetCDF4 Zip format.

North Atlantic Weather Regimes
******************************
Following k-means clustering of 500hPa geopotential height. 4 weather regimes over North Atlantic (80W-40E 30N-87.5N) are evaluated using anomalies from daily seasonal cycle. This is done retaining the first North Atlantic EOFs which explains the 80% of the variance to reduce the phase-space dimensions and then applying k-means clustering using Hartigan-Wong algorithm with k=4. Figures report patterns and frequencies of occurrence. NetCDF4 Zip data are saved. Only 4 regimes and DJF supported so far.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_miles_block.yml
* recipe_miles_eof.yml
* recipe_miles_regimes.yml

Diagnostics are stored in diag_scripts/miles/

* miles_block.R
* miles_eof.R
* miles_regimes.R

and subroutines

* basis_functions.R
* block_figures.R
* eof_figures.R
* regimes_figures.R
* block_fast.R
* eof_fast.R
* miles_parameters.R
* regimes_fast.R

`miles_parameters.R` contains additional internal parameters which affect plot sizes, colortables etc.


User settings
-------------

#. miles_block.R

   *Required settings for variables*

   * reference_dataset: reference dataset for comparison
   * reference_exp: optional reference experiment for comparison (to use when comparing two experiments of the same dataset)

   *Required settings for script*

   * seasons: Selected season('DJF','MAM','JJA','SON','ALL') or your period as e.g. 'Jan_Feb_Mar'

#. miles_eof.R

   *Required settings for variables*

   * reference_dataset: reference dataset for comparison
   * reference_exp: optional reference experiment for comparison (to use when comparing two experiments of the same dataset)

   *Required settings for script*

   * seasons: Selected season('DJF','MAM','JJA','SON','ALL') or your period as e.g. 'Jan_Feb_Mar'
   * teles: Select EOFs ('NAO','AO','PNA') or specify custom area as "lon1_lon2_lat1_lat2"

#. miles_regimes.R

    *Required settings for variables*

    * reference_dataset: reference dataset
    * reference_exp: optional reference experiment for comparison (to use when comparing two experiments of the same dataset)

    *Required or optional settings for script*

    * None (the two parameters seasons and nclusters in the recipe should not be changed)


Variables
---------

* zg (atmos, daily mean, longitude latitude time)


Observations and reformat scripts
---------------------------------
* ERA-INTERIM


References
----------
* REX, D. F. (1950), Blocking Action in the Middle Troposphere and its Effect upon Regional Climate. Tellus, 2: 196-211. doi: http://doi.org/10.1111/j.2153-3490.1950.tb00331.x
* Davini, P., C. Cagnazzo, S. Gualdi, and A. Navarra (2012): Bidimensional Diagnostics, Variability, and Trends of Northern Hemisphere Blocking. J. Climate, 25, 6496–6509, doi: http://doi.org/10.1175/JCLI-D-12-00032.1.
* Tibaldi S, Molteni F.: On the operational predictability of blocking. Tellus A 42(3): 343–365, doi: 10.1034/j.1600- 0870.1990.t01- 2- 00003.x, 1990. https://doi.org/10.1034/j.1600-0870.1990.t01-2-00003.x
* Paolo Davini. (2018, April 30). MiLES - Mid Latitude Evaluation System (Version v0.51). Zenodo. http://doi.org/10.5281/zenodo.1237838


Example plots
-------------

.. figure:: /recipes/figures/miles/miles_block.png
   :width: 14cm

   Blocking Events frequency for a CMIP5 EC-Earth historical run (DJF 1980-1989), compared to ERA-Interim. Units are percentage of blocked days per season.

.. figure:: /recipes/figures/miles/miles_eof1.png
   :width: 14cm
   
   North Atlantic Oscillation for a CMIP5 EC-Earth historical run (DJF 1980-1989) compared to ERA-Interim, shown as the linear regression of the monthly Z500 against the first Principal Component (PC1) of the North Atlantic region.
.. _recipes_toymodel:

Toymodel
========

Overview
--------

The goal of this diagnostic is to simulate single-model ensembles from an observational dataset to investigate the effect of observational uncertainty.  For further discussion of this synthetic value generator, its general application to forecasts and its limitations, see Weigel et al. (2008). The output is a netcdf file containing the synthetic observations. Due to the sampling of the perturbations from a Gaussian distribution, running the recipe multiple times, with the same observation dataset and input parameters, will result in different outputs.


Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_toymodel.yml


Diagnostics are stored in diag_scripts/magic_bsc/

* toymodel.R: generates a single model ensemble of synthetic observations




User settings
-------------

User setting files are stored in recipes/

#.	recipe_toymodel.yml

   *Required settings for preprocessor*

	extract_region:

   * start_longitude: minimum longitude
   * end_longitude: maximum longitude
   * start_latitude: minimum longitude
   * end_latitude: maximum latitude

  	extract_levels: (for 3D variables)

   * levels: [50000] # e.g. for 500 hPa level


   *Required settings for script*

   * number_of_members: integer specifying the number of members to be generated
   * beta: the user defined underdispersion (beta >= 0)


Variables
---------

* any variable (atmos/ocean, daily-monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Bellprat, O., Massonnet, F., Siegert, S., Prodhomme, C., Macias-Gómez, D., Guemas, V., & Doblas-Reyes, F. (2017). Uncertainty propagation in observational references to climate model scales. Remote Sensing of Environment, 203, 101-108.

* Massonet, F., Bellprat, O. Guemas, V., & Doblas-Reyes, F. J. (2016). Using climate models to estimate the quality of global observational data sets. Science, aaf6369.

* Weigel, A. P., Liniger, M. A., & Appenzeller, C. (2008). Can multi-model combinations really enhance the prediction skill of probabilistic ensemble forecasts? Quarterly Journal of the Royal Meteorological Society, 134(630), 241-260.


Example plots
-------------

.. _fig_toymodel:
.. figure::  /recipes/figures/toymodel/synthetic_CMIP5_bcc-csm1-1_Amon_rcp45_r1i1p1_psl_2051-2060.jpg

Twenty synthetic single-model ensemble generated by the recipe_toymodel.yml (see Section 3.7.2) for the 2051-2060 monthly data of r1i1p1 RCP 4.5 scenario of BCC_CSM1-1 simulation.
.. _recipes_capacity_factor:

Capacity factor of wind power: Ratio of average estimated power to theoretical maximum power
============================================================================================

Overview
--------

The goal of this diagnostic is to compute the wind capacity factor,  taking as input the daily instantaneous surface wind speed, which is then extrapolated to obtain the  wind speed at a height of 100 m as described in Lledó (2017).

The capacity factor is a normalized indicator of the suitability of wind speed conditions to produce electricity, irrespective of the size and number of installed turbines. This indicator is provided for three different classes of wind turbines (IEC, 2005) that are designed specifically for low, medium and high wind speed conditions.

The user can select the region, temporal range and season of interest.

The output of the recipe is a netcdf file containing the capacity factor for each of the three turbine classes.

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_capacity_factor.yml

Diagnostics are stored in diag_scripts/magic_bsc/

* capacity_factor.R: calculates the capacity factor for the three turbine classes.
* PC.R: calculates the power curves for the three turbine classes.


User settings
-------------

User setting files are stored in recipes/

#. recipe_capacity_factor.yml

   *Required settings for script*

   * power_curves: (should not be changed)

Variables
---------

* sfcWind (atmos, daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

Main features of the selected turbines:

=================  ==================  ================  ==================  =================  ===================
Turbine name       Rotor diameter (m)  Rated power (MW)  Cut-in speed (m/s)  Rated speed (m/s)  Cut-out speed (m/s)

-----------------  ------------------  ----------------  ------------------  -----------------  -------------------
Enercon E70 2.3MW  70                  2.3               2.0                 16.0               25.0
Gamesa G80 2.0MW   80                  2.0               4.0                 17.0               25.0
Gamesa G87 2.0MW   87                  2.0               4.0                 16.0               25.0
Vestas V100 2.0MW  100                 2.0               3.0                 15.0               20.0
Vestas V110 2.0MW  110                 2.0               3.0                 11.5               20.0
=================  ==================  ================  ==================  =================  ===================

References
----------

* IEC. (2005). International Standard IEC 61400-1, third edition, International Electrotechnical Commission. https://webstore.iec.ch/preview/info_iec61400-1%7Bed3.0%7Den.pdf

* Lledó, L. (2017). Computing capacity factor. Technical note BSC-ESS-2017-001, Barcelona Supercomputing Center. Available online at https://earth.bsc.es/wiki/lib/exe/fetch.php?media=library:external:bsc-ess-2017-001-c4e_capacity_factor.pdf [last accessed 11 October 2018]

Example plots
-------------

.. _fig_capfactor1:
.. figure::  /recipes/figures/capacity_factor/capacity_factor_IPSL-CM5A-MR_2021-2050.png
   :align:   center
   :width:   14cm

Wind capacity factor for five turbines: Enercon E70 (top-left), Gamesa G80 (middle-top), Gamesa G87 (top-right), Vestas V100 (bottom-left) and Vestas V110 (middle-bottom) using the IPSL-CM5A-MR simulations for the r1p1i1 ensemble for the rcp8.5 scenario during the period 2021-2050.
Emergent constraints on carbon cycle feedbacks
==============================================

Overview
--------

Figures from Wenzel et al. (2014) are reproduced with recipe_wenzel14jgr.yml. Variables relevant for the carbon cycle - climate feedback such as near surface air temperature (tas), net biosphere productivity (nbp) and carbon flux into the ocean (fgco2) are analyzed for coupled (1pctCO2, here the carbon cycle is fully coupled to the climate response) and uncoupled (esmFixCLim1, here the carbon cycle is uncoupled to the climate response) simulations. The standard namelist includes a comparison of cumulated nbp from coupled and uncoupled simulations and includes a set of routines to diagnose the long-term carbon cycle - climate feedback parameter (GammaLT) from an ensemble of CMIP5 models. Also included in the recipe is a comparison of the interannual variability of nbp and fgco2 for historical simulations used to diagnose the observable sensitivity of CO2 to tropical temperature changes (GammaIAV). As a key figure of this recipe, the diagnosed values from the models GammaLT vs. GammaIAV are compared in a scatter plot constituting an emergent constraint.


Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_wenzel14jgr.yml

Diagnostics are stored in diag_scripts/

* carbon_tsline.ncl: time line plots of annual means for spatial averages
* carbon_gammaHist.ncl: scatter plot of annual mean anomalies of two different variables; diagnosing and saving GammaIAV
* carbon_constraint.ncl: scatter plot of GammaLT vs. GammaIAV + line plot of probability density functions, diagnosing GammaLT


User settings
-------------

User setting files (cfg files) are stored in nml/cfg_carbon/

#. carbon_tsline

   *Required Settings (scripts)*

   * ts_minlat: minimum latitude for area averaging
   * ts_maxlat: maximum latitude for area averaging
   * ts_minlon: minimum longitude for area averaging
   * ts_maxlon: maximum longitude for area averaging
   * ts_maxyear: last year (time range)
   * ts_minyear: first year (time range)
   * plot_units: units to appear on Figure
   * time_avg: currently, only yearly is available
   * area_opper: type of area operation (sum)
   * styleset: Plot style

   *Optional settings (scripts)*

   * multi_model_mean: True for multi-model mean calculation
   * volcanoes: True for marking years with lage volcanic eruptions
   * align: True for aligning models to have the same start year (needed for idealized 2x CO2 simulations)
   * ts_anomaly: calculates anomalies with respect to a defined time range average (anom)
   * ridx_start: if ts_anomaly is True, define start time index for reference period
   * ridx_end: if ts_anomaly is True, define end time index for reference period
   * ref_start: if ts_anomaly is True, define start year for reference period
   * ref_end: if ts_anomaly is True, define end year for reference period

   *Required settings (variables)*

   * reference_dataset: name of reference data set

#. carbon_gammaHist.ncl

   *Required Settings (scripts)*

   * start_year: first year (time range)
   * end_year: last year (time range)
   * plot_units: units to appear on Figure
   * ec_anom: calculates anomalies with respect to the first 10-year average (anom)
   * scatter_log: set logarithmic axes in scatterplot.ncl
   * styleset: Plot style

   *Optional settings (scripts)*

   * ec_volc : exclude 2 years after volcanic erruptions (True/False)

#. carbon_constraint.ncl

   *Required Settings (scripts)*

   * gIAV_diagscript: "gammaHist_Fig3and4"
   * gIAV_start: start year of GammIAV calculation period
   * gIAV_end: end year of GammIAV calculation period
   * ec_anom: True
   * con_units: label string for units, e.g. (GtC/K)
   * nc_infile: specify path to historical gamma values derived by carbon_gammaHist.ncl
   * styleset: Plot style

   *Optional settings (scripts)*

   * reg_models: Explicit naming of individual models to be excluded from the regression


Variables
---------

* tas (atmos, monthly mean, longitude latitude time)
* nbp (land, monthly mean, longitude latitude time)
* fgco2 (ocean, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

* GCP2018: Global Carbon Budget including land (nbp) and ocean (fgco2) carbon fluxes
* NCEP: National Centers for Environmental Prediction reanalysis data for near surface temperature


References
----------

* Cox, P. M., D. B. Pearson, B. B. Booth, P. Friedlingstein, C. C. Huntingford, C. D. B. Jones, and C. M. Luke, 2013, Sensitivity of tropical carbon to climate change constrained by carbon dioxide variability, Nature, 494(7437), 341-344. doi: 10.1038/nature11882
* Wenzel, S., P. M. Cox, V. Eyring, and P. Friedlingstein, 2014, Emergent Constraints on Climate Carbon Cycle Feedbacks in the CMIP5 Earth System Models, JGR Biogeoscience, 119(5), doi: 2013JG002591.


Example plots
-------------

.. figure:: /recipes/figures/wenzel14jgr/tas_Global_CMIP5_1pctCO2_anom__1-1999.png
   :width: 10 cm
   :align: center

   Time series of tropical (30S to 30N) mean near surface temperature (tas) change between year 30 and year 110 for the CMIP5 models simulated with prescribed CO2 (1%/yr CO2 increase) coupled simulation (1pctCO2).


.. figure:: /recipes/figures/wenzel14jgr/corr_tas-nbp_anom_1960-2005.png
   :width: 10 cm
   :align: center

   Correlations between the interannual variability of global co2flux (nbp+fgco2) and tropical temperature for the individual CMIP5 models using esmHistorical simulations, and for observations.


.. figure:: /recipes/figures/wenzel14jgr/constr_tas-nbp_30-1960.000001.png
   :scale: 50 %
   :align: center

   Carbon cycle-climate feedback of tropical land carbon vs. the sensitivity of co2flux to interannual temperature variability in the tropics (30S to 30N). The red line shows the linear best fit of the regression together with the prediction error (orange shading) and the gray shading shows the observed range.


.. figure:: /recipes/figures/wenzel14jgr/constr_tas-nbp_30-1960.000002.png
   :scale: 30 %
   :align: center

   Probability Density Functions for the pure CMIP5 ensemble (black dashed) and after applying the observed constraint to the models (red solid)
.. _recipes_crem:

Cloud Regime Error Metric (CREM)
================================

Overview
--------

The radiative feedback from clouds remains the largest source of uncertainty
in determining the climate sensitivity. Traditionally, cloud has been
evaluated in terms of its impact on the mean top of atmosphere fluxes.
However it is quite possible to achieve good performance on these criteria
through compensating errors, with boundary layer clouds being too reflective
but having insufficient horizontal coverage being a common example (e.g.,
Nam et al., 2012). Williams and Webb (2009) (WW09) propose a Cloud Regime
Error Metric (CREM) which critically tests the ability of a model to
simulate both the relative frequency of occurrence and the radiative
properties correctly for a set of cloud regimes determined by the daily
mean cloud top pressure, cloud albedo and fractional coverage at each
grid-box. WW09 describe in detail how to calculate their metrics and we
have included the CREMpd metric from their paper in ESMValTool, with clear
references in the lodged code to tables in their paper. This has been
applied to those CMIP5 models who have submitted the required diagnostics
for their AMIP simulation (see Figure 8 below). As documented by WW09, a
perfect score with respect to ISCCP would be zero. WW09 also compared
MODIS/ERBE to ISCCP in order to provide an estimate of observational
uncertainty. This was found to be 0.96 and this is marked on Figure 8,
hence a model with a CREM similar to this value could be considered to have
an error comparable with observational uncertainty, although it should be
noted that this does not necessarily mean that the model lies within the
observations for each regime. A limitation of the metric is that it requires
a model to be good enough to simulate each regime. If a model is that poor
that the simulated frequency of occurrence of a particular regime is zero,
then a NaN will be returned from the code and a bar not plotted on the
figure for that model.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_williams09climdyn_CREM.yml

Diagnostics are stored in diag_scripts/crem/

* ww09_esmvaltool.py



User settings
-------------

None.


Variables
---------

* albisccp (atmos, daily mean, longitude latitude time)
* cltisccp (atmos, daily mean, longitude latitude time)
* pctisccp (atmos, daily mean, longitude latitude time)
* rlut (atmos, daily mean, longitude latitude time)
* rlutcs (atmos, daily mean, longitude latitude time)
* rsut (atmos, daily mean, longitude latitude time)
* rsutcs (atmos, daily mean, longitude latitude time)
* sic/siconc (seaice, daily mean, longitude latitude time)
* snc (atmos, daily mean, longitude latitude time)

If snc is not available then snw can be used instead. For AMIP simulations,
sic/siconc is often not submitted as it a boundary condition and effectively
the same for every model. In this case the same daily sic data set can be
used for each model.

**Note: in case of using sic/siconc data from a different model (AMIP), it has to
be checked by the user that the calendar definitions of all data sets are
compatible, in particular whether leap days are included or not.**



Observations and reformat scripts
---------------------------------

All observational data have been pre-processed and included within the
routine. These are ISCCP, ISCCP-FD, MODIS, ERBE. No additional observational
data are required at runtime.



References
----------

* Nam, C., Bony, S., Dufresne, J.-L., and Chepfer, H.: The 'too few, too bright'
  tropical low-cloud problem in CMIP5 models, Geophys. Res. Lett., 39, L21801,
  doi: 10.1029/2012GL053421, 2012.
* Williams, K.D. and Webb, M.J.: A quantitative performance assessment of
  cloud regimes in climate models. Clim. Dyn. 33, 141-157, doi:
  10.1007/s00382-008-0443-1, 2009.


Example plots
-------------

.. figure:: /recipes/figures/crem/crem_error_metric.png
   :width: 10cm
   :alt: xxxxx

   Cloud Regime Error Metrics (CREMpd) from William and Webb (2009) applied
   to those CMIP5 AMIP simulations with the required data in the archive. A
   perfect score with respect to ISCCP is zero; the dashed red line is an
   indication of observational uncertainty.
.. _recipes_clouds:

Clouds
======

Overview
--------

The recipe recipe_lauer13jclim.yml computes the climatology and interannual
variability of climate relevant cloud variables such as cloud radiative forcing
(CRE), liquid water path (lwp), cloud amount (clt), and total precipitation (pr)
reproducing some of the evaluation results of Lauer and Hamilton (2013). The
recipe includes a comparison of the geographical distribution of multi-year
average cloud parameters from individual models and the multi-model mean with
satellite observations. Taylor diagrams are generated that show the multi-year
annual or seasonal average performance of individual models and the multi-model
mean in reproducing satellite observations. The diagnostic also facilitates the
assessment of the bias of the multi-model mean and zonal averages of individual
models compared with satellite observations. Interannual variability is
estimated as the relative temporal standard deviation from multi-year timeseries
of data with the temporal standard deviations calculated from monthly anomalies
after subtracting the climatological mean seasonal cycle.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_lauer13jclim.yml

Diagnostics are stored in diag_scripts/clouds/

    * clouds.ncl: global maps of (multi-year) annual means including multi-model
      mean
    * clouds_bias.ncl: global maps of the multi-model mean and the multi-model
      mean bias
    * clouds_interannual: global maps of the interannual variability
    * clouds_isccp: global maps of multi-model mean minus observations + zonal
      averages of individual models, multi-model mean and observations
    * clouds_taylor.ncl: taylor diagrams


User settings in recipe
-----------------------

#. Script clouds.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * embracesetup: true = 2 plots per line, false = 4 plots per line (default)
   * explicit_cn_levels: explicit contour levels (array)
   * extralegend: plot legend(s) to extra file(s)
   * filename_add: optionally add this string to plot filesnames
   * panel_labels: label individual panels (true, false)
   * PanelTop: manual override for "@gnsPanelTop" used by panel plot(s)
   * projection: map projection for plotting (default =
     "CylindricalEquidistant")
   * showdiff: calculate and plot differences model - reference
     (default = false)
   * rel_diff: if showdiff = true, then plot relative differences (%)
     (default = False)
   * ref_diff_min: lower cutoff value in case of calculating relative
     differences (in units of input variable)
   * region: show only selected geographic region given as latmin, latmax,
     lonmin, lonmax
   * timemean: time averaging - "seasonal" = DJF, MAM, JJA, SON),
     "annual" = annual mean
   * treat_var_as_error: treat variable as error when averaging (true, false);
     true:  avg = sqrt(mean(var*var)), false: avg = mean(var)

   *Required settings (variables)*

   none

   * Optional settings (variables)

   * long_name: variable description
   * reference_dataset: reference dataset; REQUIRED when calculating
     differences (showdiff = True)
   * units: variable units (for labeling plot only)

   *Color tables*

   * variable "lwp": diag_scripts/shared/plot/rgb/qcm3.rgb

#. Script clouds_bias.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * plot_abs_diff: additionally also plot absolute differences (true, false)
   * plot_rel_diff: additionally also plot relative differences (true, false)
   * projection: map projection, e.g., Mollweide, Mercator
   * timemean: time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean)

   * Required settings (variables)*

   * reference_dataset: name of reference datatset

   *Optional settings (variables)*

   * long_name: description of variable

   *Color tables*

   * variable "tas": diag_scripts/shared/plot/rgb/ipcc-tas.rgb,
     diag_scripts/shared/plot/rgb/ipcc-tas-delta.rgb
   * variable "pr-mmday": diag_scripts/shared/plots/rgb/ipcc-precip.rgb,
     diag_scripts/shared/plot/rgb/ipcc-precip-delta.rgb

#. Script clouds_interannual.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * colormap: e.g., WhiteBlueGreenYellowRed, rainbow
   * explicit_cn_levels: use these contour levels for plotting
   * extrafiles: write plots for individual models to separate files
     (true, false)
   * projection: map projection, e.g., Mollweide, Mercator

   *Required settings (variables)*

   none

   *Optional settings (variables)*

   * long_name: description of variable
   * reference_dataset: name of reference datatset

   *Color tables*

   * variable "lwp": diag_scripts/shared/plots/rgb/qcm3.rgb

.. _clouds_ipcc.ncl:

#. Script clouds_ipcc.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * explicit_cn_levels: contour levels
   * mask_ts_sea_ice: true = mask T < 272 K as sea ice (only for variable "ts");
     false = no additional grid cells masked for variable "ts"
   * projection: map projection, e.g., Mollweide, Mercator
   * styleset: style set for zonal mean plot ("CMIP5", "DEFAULT")
   * timemean: time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean)
   * valid_fraction: used for creating sea ice mask (mask_ts_sea_ice = true):
     fraction of valid time steps required to mask grid cell as valid data

   *Required settings (variables)*

   * reference_dataset:  name of reference data set

   *Optional settings (variables)*

   * long_name: description of variable
   * units: variable units

   *Color tables*

   * variables "pr", "pr-mmday": diag_scripts/shared/plot/rgb/ipcc-precip-delta.rgb

#. Script clouds_taylor.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * embracelegend: false (default) = include legend in plot, max. 2 columns
     with dataset names in legend; true = write extra file with legend, max. 7
     dataset names per column in legend, alternative observational dataset(s)
     will be plotted as a red star and labeled "altern. ref. dataset" in legend
     (only if dataset is of class "OBS")
   * estimate_obs_uncertainty: true = estimate observational uncertainties
     from mean values (assuming fractions of obs. RMSE from documentation of
     the obs data); only available for "CERES-EBAF", "MODIS", "MODIS-L3";
     false = do not estimate obs. uncertainties from mean values
   * filename_add: legacy feature: arbitrary string to be added to all
     filenames of plots and netcdf output produced (default = "")
   * mask_ts_sea_ice: true = mask T < 272 K as sea ice (only for variable "ts");
     false = no additional grid cells masked for variable "ts"
   * styleset: "CMIP5", "DEFAULT" (if not set, clouds_taylor.ncl will create a
     color table and symbols for plotting)
   * timemean: time averaging; annualclim (default) = 1 plot annual mean;
     seasonalclim = 4 plots (DJF, MAM, JJA, SON)
   * valid_fraction: used for creating sea ice mask (mask_ts_sea_ice = true):
     fraction of valid time steps required to mask grid cell as valid data

   *Required settings (variables)*

   * reference_dataset: name of reference data set

   *Optional settings (variables)*

   none


Variables
---------

* clwvi (atmos, monthly mean, longitude latitude time)
* clivi (atmos, monthly mean, longitude latitude time)
* clt (atmos, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* rlut, rlutcs (atmos, monthly mean, longitude latitude time)
* rsut, rsutcs (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
instructions.*

* CERES-EBAF (obs4MIPs) - CERES TOA radiation fluxes (used for calculation of
  cloud forcing)
* GPCP-SG (obs4MIPs) - Global Precipitation Climatology Project total
  precipitation
* MODIS (obs4MIPs) - MODIS total cloud fraction
* UWisc - University of Wisconsin-Madison liquid water path climatology, based
  on satellite observbations from TMI, SSM/I, and AMSR-E, reference: O'Dell et
  al. (2008), J. Clim.

  *Reformat script:* reformat_scripts/obs/reformat_obs_UWisc.ncl

References
----------

* Flato, G., J. Marotzke, B. Abiodun, P. Braconnot, S.C. Chou, W. Collins, P.
  Cox, F. Driouech, S. Emori, V. Eyring, C. Forest, P. Gleckler, E. Guilyardi,
  C. Jakob, V. Kattsov, C. Reason and M. Rummukainen, 2013: Evaluation of
  Climate Models. In: Climate Change 2013: The Physical Science Basis.
  Contribution of Working Group I to the Fifth Assessment Report of the
  Intergovernmental Panel on Climate Change [Stocker, T.F., D. Qin, G.-K.
  Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex and
  P.M. Midgley (eds.)]. Cambridge University Press, Cambridge, United Kingdom
  and New York, NY, USA.

* Lauer A., and K. Hamilton (2013), Simulating clouds with global climate
  models: A comparison of CMIP5 results with CMIP3 and satellite data, J. Clim.,
  26, 3823-3845, doi: 10.1175/JCLI-D-12-00451.1.

* O’Dell, C.W., F.J. Wentz, and R. Bennartz (2008), Cloud liquid water path
  from satellite-based passive microwave observations: A new climatology over
  the global oceans, J. Clim., 21, 1721-1739, doi:10.1175/2007JCLI1958.1.

* Pincus, R., S. Platnick, S.A. Ackerman, R.S. Hemler, Robert J. Patrick
  Hofmann (2012), Reconciling simulated and observed views of clouds: MODIS,
  ISCCP, and the limits of instrument simulators. J. Climate, 25, 4699-4720,
  doi: 10.1175/JCLI-D-11-00267.1.


Example plots
-------------

.. _fig_cloud_1:
.. figure::  /recipes/figures/clouds/liq_h2o_path_multi.png
   :align:   center

   The 20-yr average LWP (1986-2005) from the CMIP5 historical model runs and
   the multi-model mean in comparison with the UWisc satellite climatology
   (1988-2007) based on SSM/I, TMI, and AMSR-E (O'Dell et al. 2008).

.. _fig_cloud_2:
.. figure::  /recipes/figures/clouds/liq_h2o_taylor.png
   :align:   center
   :width:   7cm

   Taylor diagram showing the 20-yr annual average performance of CMIP5 models
   for total cloud fraction as compared to MODIS satellite observations.

.. _fig_cloud_3:
.. figure::  /recipes/figures/clouds/cloud_sweffect.png
   :align:   center
   :width:   9cm

.. figure::  /recipes/figures/clouds/cloud_lweffect.png
   :align:   center
   :width:   9cm

.. figure::  /recipes/figures/clouds/cloud_neteffect.png
   :align:   center
   :width:   9cm

   20-year average (1986-2005) annual mean cloud radiative effects of CMIP5
   models against the CERES EBAF (2001–2012). Top row shows the shortwave
   effect; middle row the longwave effect, and bottom row the net effect.
   Multi-model mean biases against CERES EBAF are shown on the left, whereas the
   right panels show zonal averages from CERES EBAF (thick black), the
   individual CMIP5 models (thin gray lines) and the multi-model mean (thick
   red line). Similar to Figure 9.5 of Flato et al. (2013).

.. _fig_cloud_4:
.. figure::  /recipes/figures/clouds/cloud_var_multi.png
   :align:   center

   Interannual variability of modeled and observed (GPCP) precipitation rates
   estimated as relative temporal standard deviation from 20 years (1986-2005)
   of data. The temporal standard devitions are calculated from monthly
   anomalies after subtracting the climatological mean seasonal cycle.

.. _recipe_autoassess_stratosphere.rst:

Stratosphere - Autoassess diagnostics
=====================================

Overview
--------

Polar night jet / easterly jet strengths are defined as the maximum / minimum wind
speed of the climatological zonal mean jet, and measure how realistic the zonal
wind climatology is in the stratosphere.

Extratropical temperature at 50hPa (area averaged poleward of 60 degrees) is important
for polar stratospheric cloud formation (in winter/spring), determining the amount of
heterogeneous ozone depletion simulated by models with interactive chemistry schemes.

The Quasi-Biennial Oscillation (QBO) is a good measure of tropical variability in the
stratosphere.  Zonal mean zonal wind at 30hPa is used to define the period and amplitude
of the QBO.

The tropical tropopause cold point (100hPa, 10S-10N) temperature is an important factor in
determining the stratospheric water vapour concentrations at entry point (70hPa, 10S-10N),
and this in turn is important for the accurate simulation of stratospheric chemistry and
radiative balance.

Prior and current contributors
------------------------------
Met Office:

* Prior to May 2008: Neal Butchart
* May 2008 - May 2016: Steven C Hardiman
* Since May 2016: Alistair Sellar and Paul Earnshaw

ESMValTool:

* Since April 2018: Porting into ESMValTool by Valeriu Predoi


Developers
----------
Met Office:

* Prior to May 2008: Neal Butchart
* May 2008 - May 2016: Steven C Hardiman

ESMValTool:

* Since April 2018: Valeriu Predoi

Review of current port in ESMValTool
------------------------------------
The code and results review of the port from native Autoassess to ESMValTool
was conducted by Alistair Sellar (`<alistair.sellar@matoffice.gov.uk>`_) and
Valeriu Predoi (`<valeriu.predoi@ncas.ac.uk>`_) in July 2019. Review consisted in
comparing results from runs using ESMValTool's port and native Autoassess using
the same models and data stretches.

Metrics and Diagnostics
-----------------------

Performance metrics:

* Polar night jet: northern hem (January) vs. ERA Interim
* Polar night jet: southern hem (July) vs. ERA Interim
* Easterly jet: southern hem (January) vs. ERA Interim
* Easterly jet: northern hem (July) vs. ERA Interim
* 50 hPa temperature: 60N-90N (DJF) vs. ERA Interim
* 50 hPa temperature: 60N-90N (MAM) vs. ERA Interim
* 50 hPa temperature: 90S-60S (JJA) vs. ERA Interim
* 50 hPa temperature: 90S-60S (SON) vs. ERA Interim
* QBO period at 30 hPa vs. ERA Interim
* QBO amplitude at 30 hPa (westward) vs. ERA Interim
* QBO amplitude at 30 hPa (eastward) vs. ERA Interim
* 100 hPa equatorial temp (annual mean) vs. ERA Interim
* 100 hPa equatorial temp (annual cycle strength) vs. ERA Interim
* 70 hPa 10S-10N water vapour (annual mean) vs. ERA-Interim

Diagnostics:

* Age of stratospheric air vs. observations from Andrews et al. (2001) and Engel et al. (2009)


Model Data
----------

===========================   ================== ============== ==============================================
Variable/Field name           realm              frequency      Comment
===========================   ================== ============== ==============================================
Eastward wind (ua)            Atmosphere         monthly mean   original stash: x-wind, no stash
Air temperature (ta)          Atmosphere         monthly mean   original stash: m01s30i204
Specific humidity (hus)       Atmosphere         monthly mean   original stash: m01s30i205
===========================   ================== ============== ==============================================

The recipe takes as input a control model and experimental model, comparisons being made
with these two CMIP models; additionally it can take observational data s input, in the
current implementation ERA-Interim.

Inputs and usage
----------------
The ``stratosphere`` area metric is part of the ``esmvaltool/diag_scripts/autoassess`` diagnostics,
and, as any other ``autoassess`` metric, it uses the ``autoassess_area_base.py`` as general purpose
wrapper. This wrapper accepts a number of input arguments that are read through from the recipe.

This recipe is part of the larger group of Autoassess metrics ported to ESMValTool
from the native Autoassess package from the UK's Met Office. The ``diagnostics`` settings
are almost the same as for the other Atoassess metrics.

.. note::

   **Time gating for autoassess metrics.**

   To preserve the native Autoassess functionalities,
   data loading and selection on time is done somewhat
   differently for ESMValTool's autoassess metrics: the
   time selection is done in the preprocessor as per usual but
   a further time selection is performed as part of the diagnostic.
   For this purpose the user will specify a ``start:`` and ``end:``
   pair of arguments of ``scripts: autoassess_script`` (see below
   for example). These are formatted as ``YYYY/MM/DD``; this is
   necessary since the Autoassess metrics are computed from 1-Dec
   through 1-Dec rather than 1-Jan through 1-Jan. This is a temporary
   implementation to fully replicate the native Autoassess functionality
   and a minor user inconvenience since they need to set an extra set of
   ``start`` and ``end`` arguments in the diagnostic; this will be phased
   when all the native Autoassess metrics have been ported to ESMValTool
   review has completed.

.. note::

   **Polar Night/Easterly Jets Metrics**

   Polar Night Jets (PNJ) metrics require data available at very low air pressures
   ie very high altitudes; both Olar Night Jet and Easterly Jets computations should
   be preformed using ``ta`` and ``ua`` data at ``<< 100 Pa``; the lowest air pressure
   found in atmospheric CMOR mip tables corresponds to ``plev39`` air pressure table,
   and is used in the ``AERmonZ`` mip. If the user requires correct calculations of these
   jets, it is highly advisable to use data from ``AERmonZ``. Note that standard QBO
   calculation is exact for ``plev17`` or ``plev19`` tables.

An example of standard inputs as read by ``autoassess_area_base.py`` and passed
over to the diagnostic/metric is listed below.


.. code-block:: yaml

    scripts:
      autoassess_strato_test_1: &autoassess_strato_test_1_settings
        script: autoassess/autoassess_area_base.py  # the base wrapper
        title: "Autoassess Stratosphere Diagnostic Metric"  # title
        area: stratosphere  # assesment area
        control_model: UKESM1-0-LL-hist  # control dataset name
        exp_model: UKESM1-0-LL-piCont  # experiment dataset name
        obs_models: [ERA-Interim]  # list to hold models that are NOT for metrics but for obs operations
        additional_metrics: [ERA-Interim]  # list to hold additional datasets for metrics
        start: 2004/12/01  # start date in native Autoassess format
        end: 2014/12/01  # end date in native Autoassess format


References
----------
* Andrews, A. E., and Coauthors, 2001: Mean ages of stratospheric air derived from in situ observations of CO2, CH4, and N2O. J. Geophys. Res.,   106 (D23), 32295-32314.
* Dee, D. P., and Coauthors, 2011: The ERA-Interim reanalysis: configuration and performance of the data assimilation system. Q. J. R. Meteorol.  Soc, 137, 553-597, doi:10.1002/qj.828.
* Engel, A., and Coauthors, 2009: Age of stratospheric air unchanged within uncertainties over the past 30 years. Nat. Geosci., 2, 28-31, doi:10  .1038/NGEO388.

Observations Data sets
----------------------

ERA-Interim data (Dee et al., 2011) data can be obtained online from ECMWF and NASA respectively.  Monthly mean zonal mean U and T data are required. CMORized that exists on CEDA-Jasmin or DKRZ (contact Valeriu Predoi (`<valeriu.predoi@ncas.ac.uk>`_) for Jasmin or Mattia Righi (`<mattia.righi@dlr.de>`_ )for DKRZ).

Sample Plots and metrics
------------------------
Below is a set of metrics for  UKESM1-0-LL (historical data); the table
shows a comparison made between running ESMValTool on CMIP6 CMORized
netCDF data freely available on ESGF nodes and the run made using native
Autoassess performed at the Met Office using the pp output of the model.

===============================================     ================     ====================
Metric name                                         UKESM1-0-LL;         UKESM1-0-LL;
                                                    CMIP6: AERmonZ;      pp files;
                                                    historical, ESGF     historical, u-bc179
===============================================     ================     ====================
Polar night jet: northern hem (January)             44.86                44.91
Polar night jet: southern hem (July)                112.09               112.05
Easterly jet: southern hem (January)                76.12                75.85
Easterly jet: northern hem (July)                   55.68                55.74
QBO period at 30 hPa                                41.50                41.00
QBO amplitude at 30 hPa (westward)                  27.39                27.39
QBO amplitude at 30 hPa (eastward)                  17.36                17.36
50 hPa temperature: 60N-90N (DJF)                   27.11                26.85
50 hPa temperature: 60N-90N (MAM)                   40.94                40.92
50 hPa temperature: 90S-60S (JJA)                   11.75                11.30
50 hPa temperature: 90S-60S (SON)                   23.88                23.63
100 hPa equatorial temp (annual mean)               15.29                15.30
100 hPa equatorial temp (annual cycle strength)      1.67                 1.67
100 hPa 10Sto10N temp (annual mean)                 15.48                15.46
100 hPa 10Sto10N temp (annual cycle strength)        1.62                 1.62
70 hPa 10Sto10N wv (annual mean)                     5.75                 5.75
===============================================     ================     ====================

Results from ``u-bc179`` have been obtained by running the native Autoassess/stratosphere
on ``.pp`` data from UKESM1 ``u-bc179`` suite and are listed here to confirm the
compliance between the ported Autoassess metric in ESMValTool and the original native metric.

Another reference run comparing UKESM1-0-LL to the physical model HadGEM3-GC31-LL can be found
`here <https://github.com/NCAS-CMS/NCAS-Useful-Documentation/tree/main/autoassess_review_results/stratosphere_AERmonZ/plots/aa_strato/autoassess_strato_test_1/HadGEM3-GC31-LL_vs_UKESM1-0-LL/stratosphere>`_ .


.. figure:: /recipes/figures/autoassess_stratosphere/metrics.png
   :scale: 50 %
   :alt: metrics.png

   Standard metrics plot comparing standard metrics from UKESM1-0-LL and HadGEM3-GC31.


.. figure:: /recipes/figures/autoassess_stratosphere/UKESM1-0-LL_u_jan.png
   :scale: 50 %
   :alt: UKESM1-0-LL_u_jan.png

   Zonal mean zonal wind in January for UKESM1-0-LL.

.. figure:: /recipes/figures/autoassess_stratosphere/HadGEM3-GC31-LL_u_jan.png
   :scale: 50 %
   :alt: HadGEM3-GC31-LL_u_jan.png

   Zonal mean zonal wind in January for HadGEM3-GC31-LL.

.. figure:: /recipes/figures/autoassess_stratosphere/UKESM1-0-LL_qbo.png
   :scale: 50 %
   :alt: UKESM1-0-LL_qbo.png

   QBO for UKESM1-0-LL.

.. figure:: /recipes/figures/autoassess_stratosphere/HadGEM3-GC31-LL_qbo.png
   :scale: 50 %
   :alt: HadGEM3-GC31-LL_qbo.png

   QBO for HadGEM3-GC31-LL.

.. figure:: /recipes/figures/autoassess_stratosphere/qbo_30hpa.png
   :scale: 50 %
   :alt: qbo_30hpa.png

   QBO at 30hPa comparison between UKESM1-0-LL and HadGEM3-GC31-LL.

.. figure:: /recipes/figures/autoassess_stratosphere/teq_100hpa.png
   :scale: 50 %
   :alt: teq_100hpa.png

   Equatorial temperature at 100hPa, multi annual means.
.. _recipes_seaice_feedback:

Seaice feedback
===============


Overview
--------

In this recipe, one process-based diagnostic named the
Ice Formation Efficiency (IFE) is computed based on monthly mean
sea-ice volume estimated north of 80°N. The choice of this domain
is motivated by the desire to minimize the influence of dynamic
processes but also by the availability of sea-ice thickness measurements.
The diagnostic intends to evaluate the strength of the negative sea-ice
thickness/growth feedback, which causes late-summer negative anomalies
in sea-ice area and volume to be partially recovered during the next
growing season. A chief cause behind the existence of this feedback is
the non-linear inverse dependence between heat conduction fluxes and
sea-ice thickness, which implies that thin sea ice grows faster than thick
sea ice. To estimate the strength of that feedback, anomalies of the annual
minimum of sea-ice volume north of 80°N are first estimated. Then,
the increase in sea-ice volume until the next annual maximum is computed
for each year. The IFE is defined as the regression of this ice volume
production onto the baseline summer volume anomaly.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_seaice_feedback.yml

Diagnostics are stored in diag_scripts/seaice_feedback/

    * negative_seaice_feedback.py: scatterplot showing the feedback between
      seaice volume and seaice growth


User settings
-------------

script negative_seaice_feedback.py

    *Optional settings for script*

    * plot: dictionary containing plot options:

        - point_color: color of the plot points. (Default: black)
        - point_size: size of the plot points. (Default: 10)
        - show_values: show numerical values of feedback in plot. (Default: True)

Variables
---------

* sit (seaice, monthly mean, time latitude longitude)



References
----------

* Massonnet, F., Vancoppenolle, M., Goosse, H., Docquier, D., Fichefet, T. and Blanchard-Wrigglesworth, E., 2018.
  Arctic sea-ice change tied to its mean state through thermodynamic processes. Nature Climate Change, 8: 599-603.

Example plots
-------------

.. _fig_negative_feedback_1:
.. figure::  /recipes/figures/seaice_feedback/negative_feedback.png
   :align:   center
   :width:   14cm

   Seaice negative feedback values (CMIP5 historical experiment 1979-2004).

.. _recipes_cox18nature:

Emergent constraint on equilibrium climate sensitivity from global temperature variability
==========================================================================================

Overview
--------

This recipe reproduces the emergent constraint proposed by `Cox et al. (2018)`_
for the equilibrium climate sensitivity (ECS) using global temperature
variability. The latter is defined by a metric which can be calculated from the
global temperature variance (in time) :math:`\sigma_T` and the one-year-lag
autocorrelation of the global temperature :math:`\alpha_{1T}` by

.. math::

   \psi = \frac{\sigma_T}{\sqrt{-\ln(\alpha_{1T})}}

Using the simple `Hasselmann model`_ they show that this quantity is linearly
correlated with the ECS. Since it only depends on the temporal evolution of the
global surface temperature, there is lots of observational data available which
allows the construction of an emergent relationship. This method predicts an
ECS range of 2.2K to 3.4K (66% confidence limit).

.. _`Cox et al. (2018)`: https://www.nature.com/articles/nature25450
.. _`Hasselmann model`: https://onlinelibrary.wiley.com/doi/10.1111/j.2153-3490.1976.tb00696.x


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_cox18nature.yml


Diagnostics are stored in diag_scripts/

   * emergent_constraints/cox18nature.py
   * climate_metrics/ecs.py
   * climate_metrics/psi.py


User settings in recipe
-----------------------

* Preprocessor

   * ``area_statistics`` (*operation: mean*): Calculate global mean.

* Script emergent_constraints/cox18nature.py

   See
   :ref:`here<api.esmvaltool.diag_scripts.emergent_constraints.cox18nature>`.

* Script climate_metrics/ecs.py

   See :ref:`here<ecs.py>`.

.. _psi.py:

* Script climate_metrics/psi.py

   * ``output_attributes``, *dict*, optional: Write additional attributes to
     all output netcdf files.
   * ``lag``, *int*, optional (default: 1): Lag (in years) for the
     autocorrelation function.
   * ``window_length``, *int*, optional (default: 55): Number of years used for
     the moving window average.


Variables
---------

* *tas* (atmos, monthly, longitude, latitude, time)
* *tasa* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* HadCRUT4_ (*tasa*)

.. _HadCRUT4: https://crudata.uea.ac.uk/cru/data/temperature/


References
----------

* Cox, Peter M., Chris Huntingford, and Mark S. Williamson. "Emergent
  constraint on equilibrium climate sensitivity from global temperature
  variability." Nature 553.7688 (2018): 319.


Example plots
-------------

.. _fig_cox18nature_1:
.. figure:: /recipes/figures/cox18nature/temperature_anomaly_HadCRUT4.png
   :align: center
   :width: 50%

   Simulated change in global temperature from CMIP5 models (coloured lines),
   compared to the global temperature anomaly from the HadCRUT4 dataset (black
   dots). The anomalies are relative to a baseline period of 1961–1990. The model
   lines are colour-coded, with lower-sensitivity models (λ > 1
   Wm\ :sup:`-2`\ K\ :sup:`-1`\ ) shown by green lines and higher-sensitivity
   models (λ < 1 Wm\ :sup:`-2`\ K\ :sup:`-1`\ ) shown by magenta lines.

.. _fig_cox18nature_2:
.. figure:: /recipes/figures/cox18nature/emergent_relationship_HadCRUT4.png
   :align: center
   :width: 50%

   Emergent relationship between ECS and the ψ metric. The black dot-dashed
   line shows the best-fit linear regression across the model ensemble, with
   the prediction error for the fit given by the black dashed lines. The
   vertical blue lines show the observational constraint from the HadCRUT4
   observations: the mean (dot-dashed line) and the mean plus and minus one
   standard deviation (dashed lines).

.. _fig_cox18nature_3:
.. figure:: /recipes/figures/cox18nature/pdf_HadCRUT4.png
   :align: center
   :width: 50%

   The PDF for ECS. The orange histograms (both panels) show the prior
   distributions that arise from equal weighting of the CMIP5 models in 0.5 K
   bins.
.. _recipes_consecdrydays:

Consecutive dry days
====================

Overview
--------
Meteorological drought can in its simplest form be described by a lack of
precipitation.
First, a wet day threshold is set, which can be either a limit related to
measurement accuracy, or more directly a process related to an amount that
would break the drought.
The diagnostic calculates the longest period of consecutive dry days, which
is an indicator of the worst drought in the time series.
Further, the diagnostic calculates the frequency of dry periods longer than a
user defined number of days.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_consecdrydays.yml

Diagnostics are stored in diag_scripts/droughtindex/

    * diag_cdd.py: calculates the longest period of consecutive dry days, and
      the frequency of dry day periods longer than a user defined length


User settings in recipe
-----------------------

#. Script diag_cdd.py

   *Required settings (script)*

   * plim: limit for a day to be considered dry [mm/day]

   * frlim: the shortest number of consecutive dry days for entering statistic on frequency of dry periods.


Variables
---------

* pr      (atmos, daily mean, time latitude longitude)


Example plots
-------------

.. _fig_consecdrydays:
.. figure::  /recipes/figures/consecdrydays/consec_example_freq.png
   :align:   center
   :width:   14cm

   Example of the number of occurrences with consecutive dry days of more than five days in the period 2001 to 2002 for the CMIP5 model bcc-csm1-1-m.
.. _recipe_autoassess_landsurface_surfrad.rst:

Land-surface Surface Radiation - Autoassess diagnostics
=======================================================

Overview
--------

The simulation of surface radiation is central to all aspects of model
performance, and can often reveal compensating errors which are hidden within
top of atmosphere fluxes. This recipe provides metrics that evaluate the skill
of models' spatial and seasonal distribution of surface shortwave and longwave
radiation against the CERES EBAF satellite dataset.

Performance metrics:

* median absolute error (model minus observations) net surface shortwave (SW) radiation
* median absolute error (model minus observations) net surface longwave (LW) radiation

Metrics are calculated using model and observation multi-year climatologies (seasonal means) 
for meteorological seasons:
* December-January-February (djf)
* March-April-May (mam)
* June-July-August (jja)
* September-October-November (son)
* Annual mean (ann)


Plots:

* Normalised assessment metrics plot comparing control and experiment

The recipe takes as input a control model and experimental model, comparisons being made
with these two models.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_autoassess_landsurface_surfrad.yml

Diagnostics are stored in esmvaltool/diag_scripts/autoassess/

    * autoassess_area_base.py: wrapper for autoassess scripts
    * land_surface_surfrad/surfrad.py: script to calculate surface radiation
      metrics
    * plot_autoassess_metrics.py: plot normalised assessment metrics


User settings in recipe
-----------------------

#. Script autoassess_area_base.py

   *Required settings for script*

   * area: must equal land_surface_surfrad to select this diagnostic
   * control_model: name of model to be used as control
   * exp_model: name of model to be used as experiment
   * start: date (YYYY/MM/DD) at which period begins (see note on time gating)
   * end: date (YYYY/MM/DD) at which period ends (see note on time gating)
   * climfiles_root: path to observation climatologies

   *Optional settings for script*

   * title: arbitrary string with name of diagnostic
   * obs_models: unused for this recipe

   *Required settings for variables*

   none

   *Optional settings for variables*

   none


#. Script plot_autoassess_metrics.py

   *Required settings for script*

   * area: must equal land_surface_surfrad to select this diagnostic
   * control_model: name of model to be used as control in metrics plot
   * exp_model: name of model to be used as experiment in metrics plot
   * title: string to use as plot title

   *Optional settings for script*

   none

   *Required settings for variables*

   none

   *Optional settings for variables*

   none


Variables
---------

* rsns (atmos, monthly mean, longitude latitude time)
* rlns (atmos, monthly mean, longitude latitude time)
* sftlf (mask, fixed, longitude latitude)


Observations and reformat scripts
---------------------------------

2001-2012 climatologies (seasonal means) from CERES-EBAF Ed2.7.


References
----------
* Loeb, N. G., D. R. Doelling, H. Wang, W. Su, C. Nguyen, J. G. Corbett, L. Liang, C. Mitrescu, F. G. Rose, and S. Kato, 2018: Clouds and the Earth's Radiant Energy System (CERES) Energy Balanced and Filled (EBAF) Top-of-Atmosphere (TOA) Edition-4.0 Data Product. J. Climate, 31, 895-918, doi: 10.1175/JCLI-D-17-0208.1.

* Kato, S., F. G. Rose, D. A. Rutan, T. E. Thorsen, N. G. Loeb, D. R. Doelling, X. Huang, W. L. Smith, W. Su, and S.-H. Ham, 2018: Surface irradiances of Edition 4.0 Clouds and the Earth's Radiant Energy System (CERES) Energy Balanced and Filled (EBAF) data product, J. Climate, 31, 4501-4527, doi: 10.1175/JCLI-D-17-0523.1



Example plots
-------------

.. figure:: /recipes/figures/autoassess_landsurface/Surfrad_Metrics.png
   :scale: 50 %
   :alt: Surfrad_Metrics.png

   Normalised metrics plot comparing a control and experiment simulation



Inputs and usage
----------------
The ``landsurface_soilmoisture`` area metric is part of the ``esmvaltool/diag_scripts/autoassess`` diagnostics,
and, as any other ``autoassess`` metric, it uses the ``autoassess_area_base.py`` as general purpose
wrapper. This wrapper accepts a number of input arguments that are read through from the recipe.

This recipe is part of the larger group of Autoassess metrics ported to ESMValTool
from the native Autoassess package from the UK's Met Office. The ``diagnostics`` settings
are almost the same as for the other Autoassess metrics.

.. note::

   **Time gating for autoassess metrics.**

   To preserve the native Autoassess functionalities,
   data loading and selection on time is done somewhat
   differently for ESMValTool's autoassess metrics: the
   time selection is done in the preprocessor as per usual but
   a further time selection is performed as part of the diagnostic.
   For this purpose the user will specify a ``start:`` and ``end:``
   pair of arguments of ``scripts: autoassess_script`` (see below
   for example). These are formatted as ``YYYY/MM/DD``; this is
   necessary since the Autoassess metrics are computed from 1-Dec
   through 1-Dec rather than 1-Jan through 1-Jan. This is a temporary
   implementation to fully replicate the native Autoassess functionality
   and a minor user inconvenience since they need to set an extra set of
   ``start`` and ``end`` arguments in the diagnostic; this will be phased
   when all the native Autoassess metrics have been ported to ESMValTool
   review has completed.


An example of standard inputs as read by ``autoassess_area_base.py`` and passed
over to the diagnostic/metric is listed below.


.. code-block:: yaml

    scripts:
      autoassess_landsurf_surfrad: &autoassess_landsurf_surfrad_settings
        script: autoassess/autoassess_area_base.py
        title: "Autoassess Land-Surface Diagnostic Surfrad Metric"
        area: land_surface_surfrad
        control_model: UKESM1-0-LL
        exp_model: UKESM1-0-LL
        obs_models: [CERES-EBAF]
        obs_type: obs4MIPs
        start: 1997/12/01
        end: 2002/12/01
.. _XML_oceans:

Ocean diagnostics
=================

Overview
........

These recipes are used for evaluating the marine component of models of the
earth system. Using these recipes, it should be possible to evaluate both the
physical models and biogeochemistry models. All these recipes use the
ocean diagnostics package.

The ocean diagnostics package contains several diagnostics which produce
figures and statistical information of models of the ocean. The datasets have
been pre-processed by ESMValTool, based on recipes in the recipes directory.
Most of the diagnostics produce two or less types of figure, and several
diagnostics are called by multiple recipes.

Each diagnostic script expects a metadata file, automatically generated by
ESMValTool, and one or more pre-processed dataset. These are passed to the
diagnostic by ESMValTool in the settings.yml and metadata.yml files.

The ocean diagnostics toolkit can not figure out how to plot data by itself.
The current version requires the recipe to produce the correct pre-processed
data for each diagnostic script. ie: to produce a time series plot,
the preprocessor must produce a time-dimensional dataset.

While these tools were built to evaluate the ocean component models, they also
can be used to produce figures for other domains. However, there are some ocean
specific elements, such as the z-direction being positive and reversed, and
some of the map plots have the continents coloured in by default.

As elsewhere, both the model and observational datasets need to be
compliant with the CMOR data.

Available recipes
.................

* recipe_ocean_amoc.yml_
* recipe_ocean_example.yml_
* recipe_ocean_scalar_fields.yml_
* recipe_ocean_bgc.yml_
* recipe_ocean_quadmap.yml_
* recipe_ocean_ice_extent.yml_
* recipe_ocean_multimap.yml_


recipe_ocean_amoc.yml
---------------------

The recipe_ocean_amoc.yml_ is an recipe that produces figures describing the
Atlantic Meridional Overturning Circulation (AMOC) and the drake passage
current.

The recipes produces time series of the AMOC at 26 north and the
drake passage current.

.. centered:: |pic_amoc|

.. |pic_amoc| image:: /recipes/figures/ocean/amoc_fig_1.png


This figure shows the multi model comparison of the AMOC from several CMIP5
historical simulations, with a 6 year moving average (3 years either side of the
central value). A similar figure is produced for each individual model, and
for the Drake Passage current.

This recipe also produces a contour transect and a coloured transect plot
showing the Atlantic stream function for each individual model, and a
multi-model contour is also produced:

.. centered:: |pic_ocean_sf3| |pic_ocean_sf4|

.. |pic_ocean_sf3| image:: /recipes/figures/ocean/stream_function1.png
.. |pic_ocean_sf4| image:: /recipes/figures/ocean/stream_function2.png


recipe_ocean_example.yml
------------------------

The recipe_ocean_example.yml_ is an example recipe which shows several examples
of how to manipulate marine model data using the ocean diagnostics tools.

While several of the diagnostics here have specific uses in evaluating models,
it is meant to be a catch-all recipe demonstrating many different ways to
evaluate models.

All example calculations are performed using the ocean temperature in a three
dimensional field (thetao), or at the surface (tos). This recipe demonstrates
the use of a range of preprocessors in a marine context, and also shows many
of the standard model-only diagnostics (no observational component is included.)

This recipe includes examples of how to manipulate both 2D and 3D fields to
produce:

* Time series:

  * Global surface area weighted mean time series
  * Volume weighted average time series within a specific depth range
  * Area weighted average time series at a specific depth
  * Area weighted average time series at a specific depth in a specific region.
  * Global volume weighted average time series
  * Regional volume weighted average time series

* Maps:

  * Global surface map (from 2D ad 3D initial fields)
  * Global surface map using re-gridding to a regular grid
  * Global map using re-gridding to a regular grid at a specific depth level
  * Regional map using re-gridding to a regular grid at a specific depth level

* Transects:

  * Produce various transect figure showing a re-gridded transect plot, and multi model comparisons

* Profile:

  * Produce a Global area-weighted depth profile figure
  * Produce a regional area-weighted depth profile figure

All the these fields can be expanded using a

recipe_ocean_bgc.yml
--------------------

The recipe_ocean_bgc.yml_ is an example recipe which shows a several simple examples of how to
manipulate marine biogeochemical model data.

This recipe includes the following fields:

* Global total volume-weighted average time series:

  * temperature, salinity, nitrate, oxygen, silicate (vs WOA data) `*`
  * chlorophyll, iron, total alkalinity (no observations)

* Surface area-weighted average time series:

  * temperature, salinity, nitrate, oxygen, silicate (vs WOA data) `*`
  * fgco2 (global total), integrated primary production, chlorophyll,
    iron, total alkalinity (no observations)

* Scalar fields time series:

  * mfo (including stuff like drake passage)

* Profiles:

  * temperature, salinity, nitrate, oxygen, silicate (vs WOA data) `*`
  * chlorophyll, iron, total alkalinity (no observations)

* Maps + contours:

  * temperature, salinity, nitrate, oxygen, silicate (vs WOA data) `*`
  * chlorophyll, iron, total alkalinity (no observations)

* Transects + contours:

  * temperature, salinity, nitrate, oxygen, silicate (vs WOA data) `*`
  * chlorophyll, iron, no observations)

`*` Note that Phosphate is also available as a WOA diagnostic, but I haven't
included it as HadGEM2-ES doesn't include a phosphate field.

This recipe uses the World Ocean Atlas data, which can be downloaded from:
https://www.ncei.noaa.gov/products/world-ocean-atlas
(last access 02/08/2021)

Instructions: Select the "All fields data links (1° grid)" netCDF file,
which contain all fields.


.. recipe_OxygenMinimumZones.yml
.. ------------------------------------------
.. This recipe will appear in a future version.

.. This recipe produces an analysis of Marine oxygen. The diagnostics are based on
.. figure 1 of the following work:
.. Cabré, A., Marinov, I., Bernardello, R., and Bianchi, D.: Oxygen minimum zones
.. in the tropical Pacific across CMIP5 models: mean state differences and climate
.. change trends, Biogeosciences, 12, 5429-5454,
.. https://doi.org/10.5194/bg-12-5429-2015, 2015.


recipe_ocean_quadmap.yml
------------------------

The recipe_ocean_quadmap.yml_ is an example recipe showing the
diagnostic_maps_quad.py_ diagnostic.
This diagnostic produces an image showing four maps. Each of these four maps
show latitude vs longitude and the cube value is used as the colour scale.
The four plots are:

=================   ====================
model1              model 1 minus model2
-----------------   --------------------
model2 minus obs    model1 minus obs
=================   ====================

These figures are also known as Model vs Model vs Obs plots.


The figure produced by this recipe compares two versions of the HadGEM2 model
against ATSR sea surface temperature:

.. centered:: |pic_quad_plot|

.. |pic_quad_plot| image:: /recipes/figures/ocean/ocean_quad_plot1.png

This kind of figure can be very useful when developing a model, as it
allows model developers to quickly see the impact of recent changes
to the model.


recipe_ocean_ice_extent.yml
---------------------------

The recipe_ocean_ice_extent.yml_ recipe produces several metrics describing
the behaviour of sea ice in a model, or in multiple models.

This recipe has four preprocessors, a combinatorial combination of

* Regions: Northern or Southern Hemisphere
* Seasons: December-January-February or June-July-August

Once these seasonal hemispherical fractional ice cover is processed,
the resulting cube is passed 'as is' to the diagnostic_seaice.py_
diagnostic.

This diagnostic produces the plots:

* Polar Stereographic projection Extent plots of individual models years.
* Polar Stereographic projection maps of the ice cover and ice extent for
  individual models.
* A time series of Polar Stereographic projection Extent plots - see below.
* Time series plots of the total ice area and the total ice extent.


The following image shows an example of the sea ice extent plot, showing the
Summer Northern hemisphere ice extent for the HadGEM2-CC model, in the
historical scenario.

.. centered:: |pic_sea_ice1|

.. |pic_sea_ice1| image:: /recipes/figures/ocean/ocean_sea_ice1.png


The sea ice diagnostic is unlike the other diagnostics in the ocean diagnostics
toolkit. The other tools are build to be generic plotting tools which
work with any field (ie ``diagnostic_timeseries.py`` works fine for Temperature,
Chlorophyll, or any other field. On the other hand, the
sea ice diagnostic is the only tool that performs a field specific evaluation.

The diagnostic_seaice.py_ diagnostic is more fully described below.

recipe_ocean_multimap.yml
-------------------------

The recipe_ocean_multimap.yml_ is an example recipe showing the
diagnostic_maps_multimodel.py_ diagnostic.
This diagnostic produces an image showing Model vs Observations maps or
only Model fields when observational data are not provided.
Each map shows latitude vs longitude fields and user defined values are used to set the colour scale.
Plot layout can be modified by modifying the `layout_rowcol` argument.

The figure produced by this recipe compares the ocean surface CO2 fluxes
for 16 different CMIP5 model against Landschuetzer2016 observations.

The diagnostic_maps_multimodel.py_ diagnostic is documented below.


Available diagnostics
........................

Diagnostics are stored in the diag_scripts directory: ocean_.

The following python modules are included in the ocean diagnostics package.
Each module is described in more detail both below and inside the module.

- diagnostic_maps.py
- diagnostic_maps_quad.py
- diagnostic_model_vs_obs.py
- diagnostic_profiles.py
- diagnostic_seaice.py
- diagnostic_timeseries.py
- diagnostic_tools.py
- diagnostic_transects.py
- diagnostic_maps_multimodel.py


diagnostic_maps.py
------------------

The diagnostic_maps.py_ produces a spatial map from a NetCDF. It requires the
input netCDF to have the following dimensions. Either:

- A two dimensional file: latitude, longitude.
- A three dimensional file: depth, latitude, longitude.

In the case of a 3D netCDF file, this diagnostic produces a map for EVERY layer.
For this reason, we recommend extracting a small number of specific layers in
the preprocessor, using the `extract_layer` preprocessor.

This script can not process NetCDFs with multiple time steps. Please use the
`climate_statistics` preprocessor to collapse the time dimension.

This diagnostic also includes the optional arguments, `threshold` and
`thresholds`.

- threshold: a single float.
- thresholds: a list of floats.

Only one of these arguments should be provided at a time. These two arguments
produce a second kind of diagnostic map plot: a contour map showing the spatial
distribution of the threshold value, for each dataset. Alternatively, if the
thresholds argument is used instead of threshold, the single-dataset contour
map shows the contours of all the values in the thresholds list.

If multiple datasets are provided, in addition to the single dataset contour,
a multi-dataset contour map is also produced for each value in the thresholds
list.

Some appropriate preprocessors for this diagnostic would be:

For a  Global 2D field:

  .. code-block:: yaml

      prep_map_1:
	climate_statistics:


For a  regional 2D field:

  .. code-block:: yaml

	prep_map_2:
	    extract_region:
	      start_longitude: -80.
	      end_longitude: 30.
	      start_latitude: -80.
	      end_latitude: 80.
        climate_statistics:
          operator: mean

For a  Global 3D field at the surface and 10m depth:

  .. code-block:: yaml

	prep_map_3:
	  custom_order: true
	  extract_levels:
	    levels: [0., 10.]
	    scheme: linear_horizontal_extrapolate_vertical
      climate_statistics:
        operator: mean


For a multi-model comparison mean of 2D global fields including contour thresholds.

  .. code-block:: yaml

	prep_map_4:
	  custom_order: true
      climate_statistics:
        operator: mean
	  regrid:
	    target_grid: 1x1
	    scheme: linear

And this also requires the threshold key in the diagnostic:

  .. code-block:: yaml

	diagnostic_map:
	  variables:
	    tos: # Temperature ocean surface
	      preprocessor: prep_map_4
	      field: TO2M
	  scripts:
	    Ocean_regrid_map:
	      script: ocean/diagnostic_maps.py
	      thresholds: [5, 10, 15, 20]


diagnostic_maps_quad.py
--------------------------------

The diagnostic_maps_quad.py_ diagnostic produces an image showing four maps.
Each of these four maps show latitude vs longitude and the cube value is used
as the colour scale. The four plots are:

=================   ====================
model1              model 1 minus model2
-----------------   --------------------
model2 minus obs    model1 minus obs
=================   ====================


These figures are also known as Model vs Model vs Obs plots.

This diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cubes received by this diagnostic (via the settings.yml
and metadata.yml files) have no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An appropriate preprocessor for a 2D field would be:

  .. code-block:: yaml

	prep_quad_map:
        climate_statistics:
            operator: mean

and an example of an appropriate diagnostic section of the recipe would be:

  .. code-block:: yaml

	diag_map_1:
	  variables:
	    tos: # Temperature ocean surface
	      preprocessor: prep_quad_map
	      field: TO2Ms
	      mip: Omon
	  additional_datasets:
	#        filename: tos_ATSR_L3_ARC-v1.1.1_199701-201112.nc
	#        download from: https://datashare.is.ed.ac.uk/handle/10283/536
	    - {dataset: ATSR,  project: obs4MIPs,  level: L3,  version: ARC-v1.1.1,  start_year: 2001,  end_year: 2003, tier: 3}
	  scripts:
	    Global_Ocean_map:
	      script: ocean/diagnostic_maps_quad.py
	      control_model: {dataset: HadGEM2-CC, project: CMIP5, mip: Omon, exp: historical, ensemble: r1i1p1}
	      exper_model: {dataset: HadGEM2-ES, project: CMIP5, mip: Omon, exp: historical, ensemble: r1i1p1}
	      observational_dataset: {dataset: ATSR, project: obs4MIPs,}

Note that the details about the control model, the experiment models
and the observational dataset are all provided in the script section of the
recipe.



diagnostic_model_vs_obs.py
--------------------------------

The diagnostic_model_vs_obs.py_ diagnostic makes model vs observations maps
and scatter plots. The map plots shows four latitude vs longitude maps:

========================     =======================
Model                        Observations
------------------------     -----------------------
Model minus Observations     Model over Observations
========================     =======================

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

This diagnostic also includes the optional arguments, `maps_range` and
`diff_range` to manually define plot ranges. Both arguments are a list of two floats
to set plot range minimun and maximum values respectively for Model and Observations
maps (Top panels) and for the Model minus Observations panel (bottom left).
Note that if input data have negative values the Model over Observations map
(bottom right) is not produced.

The scatter plots plot the matched model coordinate on the x axis, and the
observational dataset on the y coordinate, then performs a linear
regression of those data and plots the line of best fit on the plot.
The parameters of the fit are also shown on the figure.

An appropriate preprocessor for a 3D+time field would be:

  .. code-block:: yaml

	preprocessors:
	  prep_map:
	    extract_levels:
	      levels:  [100., ]
	      scheme: linear_extrap
        climate_statistics:
          operator: mean
	    regrid:
	      target_grid: 1x1
	      scheme: linear



diagnostic_maps_multimodel.py
-----------------------------

The diagnostic_maps_multimodel.py_ diagnostic makes model(s) vs observations maps
and if data are not provided it draws only model field.

It is always nessary to define the overall layout trough the argument `layout_rowcol`,
which is a list of two integers indicating respectively the number of rows and columns
to organize the plot. Observations has not be accounted in here as they are automatically
added at the top of the figure.

This diagnostic also includes the optional arguments, `maps_range` and
`diff_range` to manually define plot ranges. Both arguments are a list of two floats
to set plot range minimun and maximum values respectively for variable data and
the Model minus Observations range.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An appropriate preprocessor for a 3D+time field would be:

  .. code-block:: yaml

        preprocessors:
          prep_map:
            extract_levels:
              levels:  [100., ]
              scheme: linear_extrap
        climate_statistics:
          operator: mean
            regrid:
              target_grid: 1x1
              scheme: linear



diagnostic_profiles.py
--------------------------------

The diagnostic_profiles.py_ diagnostic produces images of the profile over time from a cube.
These plots show cube value (ie temperature) on the x-axis, and depth/height
on the y axis. The colour scale is the annual mean of the cube data.
Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has a time component, and depth component, but no
latitude or longitude coordinates.

An appropriate preprocessor for a 3D+time field would be:

  .. code-block:: yaml

	preprocessors:
	  prep_profile:
	    extract_volume:
	      long1: 0.
	      long2:  20.
	      lat1:  -30.
	      lat2:  30.
	      z_min: 0.
	      z_max: 3000.
	    area_statistics:
              operator: mean



diagnostic_timeseries.py
--------------------------------

The diagnostic_timeseries.py_ diagnostic produces images of the time development
of a metric from a cube. These plots show time on the x-axis and cube value
(ie temperature) on the y-axis.

Two types of plots are produced: individual model timeseries plots and
multi model time series plots. The individual plots show the results from a
single cube, even if this cube is a multi-model mean made by the `multimodel`
preprocessor.

The multi model time series plots show several models on the same axes, where
each model is represented by a different line colour. The line colours are
determined by the number of models, their alphabetical order and the `jet`
colour scale. Observational datasets and multimodel means are shown as black
lines.

This diagnostic assumes that the preprocessors do the bulk of the work,
and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) is time-dimensional cube. This means that the pre-processed
netcdf has a time component, no depth component, and no latitude or longitude
coordinates.

Some appropriate preprocessors would be :

For a global area-weighted average 2D field:

  .. code-block:: yaml

	area_statistics:
          operator: mean

For a global volume-weighted average 3D field:

  .. code-block:: yaml

	volume_statistics:
          operator: mean

For a global area-weighted surface of a 3D field:

  .. code-block:: yaml

	extract_levels:
	  levels: [0., ]
	  scheme: linear_horizontal_extrapolate_vertical
	area_statistics:
          operator: mean


An example of the multi-model time series plots can seen here:

.. centered:: |pic_amoc2|

.. |pic_amoc2| image:: /recipes/figures/ocean/amoc_fig_1.png



diagnostic_transects.py
--------------------------------



The diagnostic_transects.py_ diagnostic produces images of a transect,
typically along a constant latitude or longitude.

These plots show 2D plots with either latitude or longitude along the x-axis,
depth along the y-axis and and the cube value is used as the colour scale.


This diagnostic assumes that the preprocessors do the bulk of the hard work,
and that the cube received by this diagnostic (via the settings.yml and
metadata.yml files) has no time component, and one of the latitude or
longitude coordinates has been reduced to a single value.

An appropriate preprocessor for a 3D+time field would be:

  .. code-block:: yaml

    climate_statistics:
      operator: mean
    extract_slice:
      latitude: [-50.,50.]
      longitude: 332.

Here is an example of the transect figure:
.. centered:: |pic_ocean_sf1|

.. |pic_ocean_sf1| image:: /recipes/figures/ocean/stream_function1.png

And here is an example of the multi-model transect contour figure:

.. centered:: |pic_ocean_sf2|

.. |pic_ocean_sf2| image:: /recipes/figures/ocean/stream_function2.png



diagnostic_seaice.py
--------------------------------



The diagnostic_seaice.py_ diagnostic is unique in this module, as it produces
several different kinds of images, including time series, maps, and contours.
It is a good example of a diagnostic where the preprocessor does very little
work, and the diagnostic does a lot of the hard work.

This was done purposely, firstly to demonstrate the flexibility of ESMValTool,
and secondly because Sea Ice is a unique field where several Metrics can be
calculated from the sea ice cover fraction.

The recipe Associated with with diagnostic is the recipe_SeaIceExtent.yml.
This recipe contains 4 preprocessors which all perform approximately the same
calculation. All four preprocessors extract a season:
- December, January and February (DJF)
- June, July and August (JJA)
and they also extract either the North or South hemisphere. The four
preprocessors are combinations of DJF or JJA and North or South hemisphere.

One of the four preprocessors is North Hemisphere Winter ice extent:

.. code-block:: yaml

	timeseries_NHW_ice_extent: # North Hemisphere Winter ice_extent
	  custom_order: true
	  extract_time: &time_anchor # declare time here.
	      start_year: 1960
	      start_month: 12
	      start_day: 1
	      end_year: 2005
	      end_month: 9
	      end_day: 31
	  extract_season:
	    season: DJF
	  extract_region:
	    start_longitude: -180.
	    end_longitude: 180.
	    start_latitude: 0.
	    end_latitude: 90.

Note that the default settings for ESMValTool assume that the year starts on the
first of January. This causes a problem for this preprocessor, as the first
DJF season would not include the first Month, December, and the final would not
include both January and February. For this reason, we also add the
`extract_time` preprocessor.

This preprocessor group produces a 2D field with a time component, allowing
the diagnostic to investigate the time development of the sea ice extend.

The diagnostic section of the recipe should look like this:

.. code-block:: yaml

	diag_ice_NHW:
	  description: North Hemisphere Winter Sea Ice diagnostics
	  variables:
	    sic: # surface ice cover
	      preprocessor: timeseries_NHW_ice_extent
	      field: TO2M
	      mip: OImon
	  scripts:
	    Global_seaice_timeseries:
	      script: ocean/diagnostic_seaice.py
	      threshold: 15.

Note the the threshold here is 15%, which is the standard cut of for the
ice extent.

The sea ice diagnostic script produces three kinds of plots, using the
methods:

- `make_map_extent_plots`: extent maps plots of individual models using a Polar Stereographic project.
- `make_map_plots`: maps plots of individual models using a Polar Stereographic project.
- `make_ts_plots`: time series plots of individual models

There are no multi model comparisons included here (yet).



diagnostic_tools.py
-------------------



The diagnostic_tools.py_ is a module that contains several python tools used
by the ocean diagnostics tools.

These tools are:

- folder: produces a directory at the path provided and returns a string.
- get_input_files: loads a dictionary from the input files in the metadata.yml.
- bgc_units: converts to sensible units where appropriate (ie Celsius, mmol/m3)
- timecoord_to_float: Converts time series to decimal time ie: Midnight on January 1st 1970 is 1970.0
- add_legend_outside_right: a plotting tool, which adds a legend outside the axes.
- get_image_format: loads the image format, as defined in the global user config.yml.
- get_image_path: creates a path for an image output.
- make_cube_layer_dict: makes a dictionary for several layers of a cube.

We just show a simple description here, each individual function is more fully
documented in the diagnostic_tools.py_ module.


A note on the auxiliary data directory
......................................

Some of these diagnostic scripts may not function on machines with no access
to the internet, as cartopy may try to download the shape files. The solution
to this issue is the put the relevant cartopy shapefiles in a directory which
is visible to esmvaltool, then link that path to ESMValTool via
the `auxiliary_data_dir` variable in your config-user.yml file.

The cartopy masking files can be downloaded from:
https://www.naturalearthdata.com/downloads/


In these recipes, cartopy uses the 1:10, physical coastlines and land files::

      110m_coastline.dbf
      110m_coastline.shp
      110m_coastline.shx
      110m_land.dbf
      110m_land.shp
      110m_land.shx


Associated Observational datasets
........................................

The following observations datasets are used by these recipes:

World Ocean ATLAS
-----------------
These data can be downloaded from:
https://www.nodc.noaa.gov/OC5/woa13/woa13data.html
(last access 10/25/2018)
Select the "All fields data links (1° grid)" netCDF file, which contain all
fields.

The following WOA datasets are used by the ocean diagnostics:
 - Temperature
 - Salinity
 - Nitrate
 - Phosphate
 - Silicate
 - Dissolved Oxygen

These files need to be reformatted using the `cmorize_obs_py` script with output name `WOA`.


Landschuetzer 2016
------------------
These data can be downloaded from:
ftp://ftp.nodc.noaa.gov/nodc/archive/arc0105/0160558/1.1/data/0-data/spco2_1998-2011_ETH_SOM-FFN_CDIAC_G05.nc
(last access 02/28/2019)

The following variables are used by the ocean diagnostics:
 - fgco2, Surface Downward Flux of Total CO2
 - spco2, Surface Aqueous Partial Pressure of CO2
 - dpco2, Delta CO2 Partial Pressure

The file needs to be reformatted using the `cmorize_obs_py` script with output name `Landschuetzer2016`.



.. Links:

.. Recipes:
.. _recipe_ocean_amoc.yml: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/recipes/recipe_ocean_amoc.yml
.. _recipe_ocean_example.yml: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/recipes/recipe_ocean_example.yml
.. _recipe_ocean_scalar_fields.yml: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/recipes/recipe_ocean_scalar_fields.yml
.. _recipe_ocean_bgc.yml: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/recipes/recipe_ocean_bgc.yml
.. _recipe_ocean_quadmap.yml: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/recipes/recipe_ocean_quadmap.yml
.. _recipe_ocean_Landschuetzer2016.yml: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/recipes/recipe_ocean_Landschuetzer2016.yml
.. _recipe_ocean_multimap.yml: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/recipes/recipe_ocean_multimap.yml

.. Diagnostics:
.. _ocean: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/:
.. _diagnostic_maps.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_maps.py
.. _diagnostic_maps_quad.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_maps_quad.py
.. _diagnostic_model_vs_obs.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_model_vs_obs.py
.. _diagnostic_maps_multimodel.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_maps_multimodel.py
.. _diagnostic_profiles.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_profiles.py
.. _diagnostic_timeseries.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_timeseries.py
.. _diagnostic_transects.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_transects.py
.. _diagnostic_seaice.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_seaice.py
.. _diagnostic_tools.py: https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/ocean/diagnostic_tools.py
.. _recipe_carvalhais14nat:

Turnover time of carbon over land ecosystems
============================================

Overview
--------

This recipe evaluates the turnover time of carbon over
land ecosystems (tau_ctotal) based on the analysis of
`Carvalhais et al. (2014)`_. In summary, it provides an overview on:

    * Comparisons of global distributions of tau_ctotal from all models against
      observation and other models
    * Variation of tau_ctotal across latitude (zonal distributions)
    * Variation of association of tau_ctotal and climate across latitude
      (zonal correlations)
    * metrics of global tau_ctotal and correlations


.. _tau calculation:

Calculation of turnover time
----------------------------

First, the total carbon content of land ecosystems is calculated as,

.. math::

 ctotal = cSoil + cVeg

where :math:`cSoil` and :math:`cVeg` are the carbon contents in soil and
vegetation. **Note that this is not fully consistent with `Carvalhais et al.
(2014)`_, in which `ctotal` includes all carbon storages that respire to the
atmosphere. Due to inconsistency across models, it resulted in having different
carbon storage components in calculation of ctotal for different models**.

The turnover time of carbon is then calculated as,

.. math::

 \tau_{ctotal} = \frac{ctotal}{gpp}

where `ctotal` and `gpp` are temporal means of total carbon content and
gross primary productivity, respectively. **The equation
is valid for steady state, and is only applicable when both ctotal and gpp
are long-term averages.** Therefore, the recipe should always include the mean
operator of climate_statistics in preprocessor.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_carvalhais14nat.yml


Diagnostics are stored in diag_scripts/

   * land_carbon_cycle/diag_global_turnover.py
   * land_carbon_cycle/diag_zonal_turnover.py
   * land_carbon_cycle/diag_zonal_correlation.py


User settings in recipe
-----------------------

Observation-related details
............................

The settings needed for loading the observational dataset in all diagnostics
are provided in the recipe through `obs_info` within `obs_details` section.

    * ``obs_data_subdir``: subdirectory of auxiliary_data_dir (set in
      config-user file) where observation data are stored {e.g.,
      data_ESMValTool_Carvalhais2014}.
    * ``source_label``: source data label {'Carvalhais2014'}.
    * ``variant_label``: variant of the observation {'BE'} for best estimate.
    * ``grid_label``: label denoting the spatial grid specification {'gn'}.
    * ``frequency``: temporal frequency of the observation data {'fx'}

The observation data file used in the recipe should be changed through the
fields above, as these are used to generate observation file name and
locations. For details, see :ref:`observations` section.

Preprocessor
............

   * ``climate_statistics``: {mean} - calculate the mean over full time period.
   * ``regrid``: {nearest} - nearest neighbor regridding to the selected
     observation resolution.
   * ``mask_landsea``: {sea} - mask out all the data points from sea.
   * ``multi_model_statistics``: {median} - calculate and include the
     multimodel median.


Script land_carbon_cycle/diag_global_turnover.py
................................................

  * Required settings:

    * ``obs_variable``: {``str``} list of the variable(s) to be read from the
      observation files

  * Optional settings:

    * ``ax_fs``: {``float``, 7.1} - fontsize in the figure.
    * ``fill_value``: {``float``, nan} - fill value to be used in analysis and
      plotting.
    * ``x0``: {``float``, 0.02} - X - coordinate of the left edge of the figure.
    * ``y0``: {``float``, 1.0} Y - coordinate of the upper edge of the figure.
    * ``wp``: {``float``, 1 / number of models} - width of each map.
    * ``hp``: {``float``, = wp} - height of each map.
    * ``xsp``: {``float``, 0} - spacing betweeen maps in X - direction.
    * ``ysp``: {``float``, -0.03} - spacing between maps in Y -direction.
      Negative to reduce the spacing below default.
    * ``aspect_map``: {``float``, 0.5} - aspect of the maps.
    * ``xsp_sca``: {``float``, wp / 1.5} - spacing between the scatter plots in
      X - direction.
    * ``ysp_sca``: {``float``, hp / 1.5} - spacing between the scatter plots in
      Y - direction.
    * ``hcolo``: {``float``, 0.0123} - height (thickness for horizontal
      orientation) of the colorbar .
    * ``wcolo``: {``float``, 0.25} - width (length) of the colorbar.
    * ``cb_off_y``: {``float``, 0.06158} - distance of colorbar from top of the
      maps.
    * ``x_colo_d``: {``float``, 0.02} - X - coordinate of the colorbar for maps
      along the diagonal (left).
    * ``x_colo_r``: {``float``, 0.76} - Y - coordinate of the colorbar for
      ratio maps above the diagonal (right).
    * ``y_colo_single``: {``float``, 0.1086} - Y-coordinate of the colorbar in
      the maps per model (separate figures).
    * ``correlation_method``: {``str``, spearman | pearson} - correlation
      method to be used while calculating the correlation displayed in the
      scatter plots.
    * ``tx_y_corr``: {``float``, 1.075} - Y - coordinate of the inset text of
      correlation.
    * ``valrange_sc``: {``tuple``, (2, 256)} - range of turnover times in X -
      and Y - axes of scatter plots.
    * ``obs_global``: {``float``, 23} - global turnover time, provided as
      additional info for map of the observation.  For models, they are
      calculated within the diagnostic.
    * ``gpp_threshold``: {``float``, 0.01} - The threshold of gpp in
      `kg m^{-2} yr^{-1}` below which the grid cells are masked.


Script land_carbon_cycle/diag_zonal_turnover.py
...............................................

  * Required settings:

    * ``obs_variable``: {``str``} list of the variable(s) to be read from the
      observation files

  * Optional settings:

    * ``ax_fs``: {``float``, 7.1} - fontsize in the figure.
    * ``fill_value``: {``float``, nan} - fill value to be used in analysis and
      plotting.
    * ``valrange_x``: {``tuple``, (2, 1000)} - range of turnover values in the
      X - axis.
    * ``valrange_y``: {``tuple``, (-70, 90)} - range of latitudes in the Y -
      axis.
    * ``bandsize``: {``float``, 9.5} - size of the latitudinal rolling window
      in degrees. One latitude row if set to ``None``.
    * ``gpp_threshold``: {``float``, 0.01} - The threshold of gpp in
      `kg m^{-2} yr^{-1}` below which the grid cells are masked.


Script land_carbon_cycle/diag_zonal_correlation.py
..................................................

  * Required settings:

    * ``obs_variable``: {``str``} list of the variable(s) to be read from the
      observation files

  * Optional settings:

    * ``ax_fs``: {``float``, 7.1} - fontsize in the figure.
    * ``fill_value``: {``float``, nan} - fill value to be used in analysis and
      plotting.
    * ``correlation_method``: {``str``, pearson | spearman} - correlation
      method to be used while calculating the zonal correlation.
    * ``min_points_frac: {``float``, 0.125} - minimum fraction of valid points
      within the latitudinal band for calculation of correlation.
    * ``valrange_x``: {``tuple``, (-1, 1)} - range of correlation values in the
      X - axis.
    * ``valrange_y``: {``tuple``, (-70, 90)} - range of latitudes in the Y -
      axis.
    * ``bandsize``: {``float``, 9.5} - size of the latitudinal rolling window
      in degrees. One latitude row if set to ``None``.
    * ``gpp_threshold``: {``float``, 0.01} - The threshold of gpp in
      `kg m^{-2} yr^{-1}` below which the grid cells are masked.


Required Variables
------------------

* *tas* (atmos, monthly, longitude, latitude, time)
* *pr* (atmos, monthly, longitude, latitude, time)
* *gpp* (land, monthly, longitude, latitude, time)
* *cVeg* (land, monthly, longitude, latitude, time)
* *cSoil* (land, monthly, longitude, latitude, time)

.. _observations:

Observations
------------

The observations needed in the diagnostics are publicly available for download
from the `Data Portal of the Max Planck Institute for Biogeochemistry <http://
www.bgc-jena.mpg.de/geodb/BGI/tau4ESMValTool.php>`_ after registration.

Due to inherent dependence of the diagnostic on uncertainty estimates in
observation, the data needed for each diagnostic script are processed at
different spatial resolutions (as in Carvalhais et al., 2014), and provided in
11 different resolutions (see Table 1). Note that the uncertainties were
estimated at the resolution of the selected models, and, thus, only the 
pre-processed observed data can be used with the recipe. 
It is not possible to use regridding functionalities of ESMValTool to regrid 
the observational data to other spatial resolutions, as the uncertainty 
estimates cannot be regridded.

Table 1. A summary of the observation datasets at different resolutions.

+-------------+---------------+-------------+
| Reference   | target_grid   | grid_label* |
+=============+===============+=============+
| Observation |     0.5x0.5   | gn          |
+-------------+---------------+-------------+
| NorESM1-M   |   2.5x1.875   | gr          |
+-------------+---------------+-------------+
| bcc-csm1-1  | 2.812x2.813   | gr1         |
+-------------+---------------+-------------+
| CCSM4       |   1.25x0.937  | gr2         |
+-------------+---------------+-------------+
| CanESM2     | 2.812x2.813   | gr3         |
+-------------+---------------+-------------+
| GFDL-ESM2G  |   2.5x2.0     | gr4         |
+-------------+---------------+-------------+
| HadGEM2-ES  | 1.875x1.241   | gr5         |
+-------------+---------------+-------------+
| inmcm4      |   2.0x1.5     | gr6         |
+-------------+---------------+-------------+
| IPSL-CM5A-MR|   2.5x1.259   | gr7         |
+-------------+---------------+-------------+
| MIROC-ESM   | 2.812x2.813   | gr8         |
+-------------+---------------+-------------+
| MPI-ESM-LR  | 1.875x1.875   | gr9         |
+-------------+---------------+-------------+

\* The grid_label is suffixed with z for data in zonal/latitude coordinates:
the zonal turnover and zonal correlation.

**To change the spatial resolution of the evaluation, change {grid_label} in
obs_details and the corresponding {target_grid} in regrid preprocessor of the
recipe**.


At each spatial resolution, four data files are provided:

  * ``tau_ctotal_fx_Carvalhais2014_BE_gn.nc`` - global data of tau_ctotal
  * ``tau_ctotal_fx_Carvalhais2014_BE_gnz.nc`` - zonal data of tau_ctotal
  * ``r_tau_ctotal_tas_fx_Carvalhais2014_BE_gnz.nc`` - zonal correlation of
    tau_ctotal and tas, controlled for pr
  * ``r_tau_ctotal_pr_fx_Carvalhais2014_BE_gnz.nc`` - zonal correlation of
    tau_ctotal
    and pr, controlled for tas.

The data is produced in obs4MIPs standards, and provided in netCDF4 format.
The filenames use the convention:

``{variable}_{frequency}_{source_label}_{variant_label}_{grid_label}.nc``

  * {variable}: variable name, set in every diagnostic script as obs_variable
  * {frequency}: temporal frequency of data, set from obs_details
  * {source_label}: observational source, set from obs_details
  * {variant_label}: observation variant, set from obs_details
  * {grid_label}: temporal frequency of data, set from obs_details

Refer to the `Obs4MIPs Data Specifications`_  for details of the definitions above.

All data variables have additional variables ({variable}_5 and {variable}_95)
in the same file. These variables are necessary for a successful execution of
the diagnostics.

References
----------

* Carvalhais, N., et al. (2014), Global covariation of carbon turnover times
  with climate in terrestrial ecosystems, Nature, 514(7521), 213-217,
  doi: 10.1038/nature13731.

.. _`Carvalhais et al. (2014)`: https://doi.org/10.1038/nature13731

.. _`Obs4MIPs Data Specifications`:
  https://esgf-node.llnl.gov/site_media/projects/obs4mips/ODSv2p1.pdf


Example plots
-------------

.. _fig_carvalhais14nat_1:
.. figure:: /recipes/figures/carvalhais14nat/r_tau_ctotal_climate_pearson_Carvalhais2014_gnz.png
   :align: center
   :width: 80%

   Comparison of latitudinal (zonal) variations of pearson correlation between
   turnover time and climate: turnover time and precipitation, controlled for
   temperature (left) and vice-versa (right). Reproduces figures 2c and 2d in 
   `Carvalhais et al. (2014)`_.

.. _fig_carvalhais14nat_2:

.. figure:: /recipes/figures/carvalhais14nat/global_matrix_map_ecosystem_carbon_turnover_time_Carvalhais2014_gn.png
   :align: center
   :width: 80%

   Comparison of observation-based and modelled ecosystem carbon turnover time.
   Along the diagnonal, tau_ctotal are plotted, above the bias, and below
   density plots. The inset text in density plots indicate the correlation. 

.. _fig_carvalhais14nat_3:

.. figure:: /recipes/figures/carvalhais14nat/global_multimodelAgreement_ecosystem_carbon_turnover_time_Carvalhais2014_gn.png
   :align: center
   :width: 80%

   Global distributions of multimodel bias and model agreement. Multimodel bias 
   is calculated as the ratio of multimodel median turnover time and that from 
   observation.  Stippling indicates the regions where only less than one 
   quarter of the models fall within the range of observational uncertainties 
   (`5^{th}` and `95^{th}` percentiles). Reproduces figure 3 in `Carvalhais et 
   al. (2014)`_.

.. _fig_carvalhais14nat_4:

.. figure:: /recipes/figures/carvalhais14nat/zonal_mean_ecosystem_carbon_turnover_time_Carvalhais2014_gnz.png
   :align: center
   :width: 80%

   Comparison of latitudinal (zonal) variations of observation-based and 
   modelled ecosystem carbon turnover time. The zonal turnover time is 
   calculated as the ratio of zonal `ctotal` and `gpp`. Reproduces figures 2a 
   and 2b in `Carvalhais et al. (2014)`_.
.. _recipes_bock20jgr:

Quantifying progress across different CMIP phases
=================================================

Overview
--------

The recipe recipe_bock20jgr.yml generates figures to quantify the progress across
different CMIP phases.

.. note::
   The current recipe uses a horizontal 5x5 grid for figure 10, while the
   original plot in the paper shows a 2x2 grid. This is solely done for
   computational reasons (running the recipe with a 2x2 grid for figure 10
   takes considerably more time than running it with a 5x5 grid) and can be
   easily changed in the preprocessor section of the recipe if necessary.



Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/bock20jgr

    * recipe_bock20jgr_fig_1-4.yml
    * recipe_bock20jgr_fig_6-7.yml
    * recipe_bock20jgr_fig_8-10.yml

Diagnostics are stored in diag_scripts/

    Fig. 1:

    * bock20jgr/tsline.ncl: timeseries of global mean surface temperature
      anomalies

    Fig. 2:

    * bock20jgr/tsline_collect.ncl: collect different timeseries from
      tsline.ncl to compare different models ensembles

    Fig. 3 and 4:

    * bock20jgr/model_bias.ncl: global maps of the multi-model mean and the
      multi-model mean bias

    Fig. 6:

    * perfmetrics/main.ncl
    * perfmetrics/collect.ncl

    Fig. 7:

    * bock20jgr/corr_pattern.ncl: calculate pattern correlation
    * bock20jgr/corr_pattern_collect.ncl: create pattern correlation plot

    Fig. 8:

    * climate_metrics/ecs.py
    * climate_metrics/create_barplot.py

    Fig. 9:

    * clouds/clouds_ipcc.ncl

    Fig. 10:

    * climate_metrics/feedback_parameters.py


User settings in recipe
-----------------------

#. Script tsline.ncl

   *Required settings (scripts)*

   * styleset: as in diag_scripts/shared/plot/style.ncl functions

   *Optional settings (scripts)*

   * time_avg: type of time average (currently only "yearly" and "monthly" are
     available).
   * ts_anomaly: calculates anomalies with respect to the defined reference
     period; for each gird point by removing the mean for the given
     calendar month (requiring at least 50% of the data to be
     non-missing)
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * ref_value: if true, right panel with mean values is attached
   * ref_mask: if true, model fields will be masked by reference fields
   * region: name of domain
   * plot_units: variable unit for plotting
   * y_min: set min of y-axis
   * y_max: set max of y-axis
   * mean_nh_sh: if true, calculate first NH and SH mean
   * volcanoes: if true, lines of main volcanic eruptions will be added
   * header: if true, use region name as header
   * write_stat: if true, write multi-model statistics to nc-file

   *Required settings (variables)*

   none

   * Optional settings (variables)

   none

#. Script tsline_collect.ncl

   *Required settings (scripts)*

   * styleset: as in diag_scripts/shared/plot/style.ncl functions

   *Optional settings (scripts)*

   * time_avg: type of time average (currently only "yearly" and "monthly" are
     available).
   * ts_anomaly: calculates anomalies with respect to the defined period
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * region: name of domain
   * plot_units: variable unit for plotting
   * y_min: set min of y-axis
   * y_max: set max of y-axis
   * order: order in which experiments should be plotted
   * header: if true, region name as header
   * stat_shading: if true: shading of statistic range
   * ref_shading: if true: shading of reference period


   *Required settings (variables)*

   none

   * Optional settings (variables)

   none

#. Script model_bias.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * projection: map projection, e.g., Mollweide, Mercator
   * timemean: time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean)

   * Required settings (variables)*

   * reference_dataset: name of reference datatset

   *Optional settings (variables)*

   * long_name: description of variable

   *Color tables*

   * variable "tas": diag_scripts/shared/plot/rgb/ipcc-ar6_temperature_div.rgb,
   * variable "pr-mmday": diag_scripts/shared/plots/rgb/ipcc-ar6_precipitation_seq.rgb
     diag_scripts/shared/plot/rgb/ipcc-ar6_precipitation_div.rgb

#. Script perfmetrics_main.ncl

   See :ref:`here<perf-main.ncl>`.

#. Script perfmetrics_collect.ncl

   See :ref:`here<perf-collect.ncl>`.

#. Script corr_pattern.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * plot_median

   *Required settings (variables)*

   * reference_dataset

   *Optional settings (variables)*

   * alternative_dataset

#. Script corr_pattern_collect.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * diag_order

   *Color tables*

   * diag_scripts/shared/plot/rgb/ipcc-ar6_line_03.rgb

#. Script ecs.py

   See :ref:`here<ecs.py>`.

#. Script create_barplot.py

   See :ref:`here<create_barplot.py>`.

#. Script clouds_ipcc.ncl

   See :ref:`here<clouds_ipcc.ncl>`.

#. Script feedback_parameters.py

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * calculate_mmm: *bool* (default: ``True``). Calculate multi-model means.
   * only_consider_mmm: *bool* (default: ``False``). Only consider multi-model
     mean dataset. This automatically sets ``calculate_mmm`` to ``True``. For
     large multi-dimensional datasets, this might significantly reduce the
     computation time if only the multi-model mean dataset is relevant.
   * output_attributes: *dict*. Write additional attributes to netcdf files.
   * seaborn_settings: *dict*. Options for :func:`seaborn.set` (affects all
     plots).


Variables
---------

* clt (atmos, monthly, longitude latitude time)
* hus (atmos, monthly, longitude latitude lev time)
* pr (atmos, monthly, longitude latitude time)
* psl (atmos, monthly, longitude latitude time)
* rlut (atmos, monthly, longitude latitude time)
* rsdt (atmos, monthly, longitude latitude time)
* rsut (atmos, monthly, longitude latitude time)
* rtmt (atmos, monthly, longitude latitude time)
* rlutcs (atmos, monthly, longitude latitude time)
* rsutcs (atmos, monthly, longitude latitude time)
* ta (atmos, monthly, longitude latitude lev time)
* tas (atmos, monthly, longitude latitude time)
* ts (atmos, monthly, longitude latitude time)
* ua (atmos, monthly, longitude latitude lev time)
* va (atmos, monthly, longitude latitude lev time)
* zg (atmos, monthly, longitude latitude time)


Observations and reformat scripts
---------------------------------

* AIRS (obs4MIPs) - specific humidity

* CERES-EBAF (obs4MIPs) - CERES TOA radiation fluxes (used for calculation of
  cloud forcing)

* ERA-Interim - reanalysis of surface temperature, sea surface pressure

  *Reformat script:* recipes/cmorizers/recipe_era5.yml

* ERA5 - reanalysis of surface temperature

  *Reformat script:* recipes/cmorizers/recipe_era5.yml

* ESACCI-CLOUD - total cloud cover

  *Reformat script:* cmorizers/obs/cmorize_obs_esacci_cloud.ncl

* ESACCI-SST - sea surface temperature

  *Reformat script:* cmorizers/obs/cmorize_obs_esacci_sst.ncl

* GHCN - Global Historical Climatology Network-Monthly gridded land precipitation

  *Reformat script:* cmorizers/obs/cmorize_obs_ghcn.ncl

* GPCP-SG (obs4MIPs) - Global Precipitation Climatology Project total
  precipitation

* HadCRUT4 - surface temperature anomalies

  *Reformat script:* cmorizers/obs/cmorize_obs_hadcrut4.ncl

* HadISST - surface temperature

  *Reformat script:* cmorizers/obs/cmorize_obs_hadisst.ncl

* JRA-55 (ana4mips) - reanalysis of sea surface pressure

* NCEP - reanalysis of surface temperature

  *Reformat script:* cmorizers/obs/cmorize_obs_NCEP.ncl

* PATMOS-x - total cloud cover

  *Reformat script:* cmorizers/obs/cmorize_obs_patmos_x.ncl


References
----------

* Bock, L., Lauer, A., Schlund, M., Barreiro, M., Bellouin, N., Jones, C.,
  Predoi, V., Meehl, G., Roberts, M., and Eyring, V.: Quantifying progress
  across different CMIP phases with the ESMValTool, Journal of Geophysical
  Research: Atmospheres, 125, e2019JD032321. https://doi.org/10.1029/2019JD032321

* Copernicus Climate Change Service (C3S), 2017: ERA5: Fifth generation of
  ECMWF atmospheric reanalyses of the global climate, edited, Copernicus
  Climate Change Service Climate Data Store (CDS).
  https://cds.climate.copernicus.eu/cdsapp#!/home

* Flato, G., J. Marotzke, B. Abiodun, P. Braconnot, S.C. Chou, W. Collins, P.
  Cox, F. Driouech, S. Emori, V. Eyring, C. Forest, P. Gleckler, E. Guilyardi,
  C. Jakob, V. Kattsov, C. Reason and M. Rummukainen, 2013: Evaluation of
  Climate Models. In: Climate Change 2013: The Physical Science Basis.
  Contribution of Working Group I to the Fifth Assessment Report of the
  Intergovernmental Panel on Climate Change [Stocker, T.F., D. Qin, G.-K.
  Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex and
  P.M. Midgley (eds.)]. Cambridge University Press, Cambridge, United Kingdom
  and New York, NY, USA.

* Morice, C. P., Kennedy, J. J., Rayner, N. A., & Jones, P., 2012: Quantifying
  uncertainties in global and regional temperature change using an ensemble of
  observational estimates: The HadCRUT4 data set, Journal of Geophysical
  Research, 117, D08101. https://doi.org/10.1029/2011JD017187


Example plots
-------------

.. _fig_bock20jgr_1:
.. figure::  /recipes/figures/bock20jgr/tas_Global_CMIP6_historical_anom_1850-2014.png
   :align:   center

   Observed and simulated time series of the anomalies in annual and global mean
   surface temperature. All anomalies are differences from the 1850-1900 time
   mean of each individual time series (Fig. 1).

.. _fig_bock20jgr_2:
.. figure::  /recipes/figures/bock20jgr/tas_Global_multimodel_anom_1850-2017.png
   :align:   center
   :width:   7cm

   Observed and simulated time series of the anomalies in annual
   and global mean surface temperature as in Figure 1; all anomalies are
   calculated by subtracting the 1850-1900 time mean from the time series.
   Displayed are the multimodel means of all three CMIP ensembles with
   shaded range of the respective standard deviation. In black the HadCRUT4
   data set (HadCRUT4; Morice et al., 2012). Gray shading shows the 5% to
   95% confidence interval of the combined effects of all the uncertainties
   described in the HadCRUT4 error model (measurement and sampling, bias,
   and coverage uncertainties) (Morice et al., 2012) (Fig. 2).

.. _fig_bock20jgr_3:
.. figure::  /recipes/figures/bock20jgr/model_bias_tas_annual_CMIP6.png
   :align:   center
   :width:   9cm

   Annual mean near‐surface (2 m) air temperature (°C). (a) Multimodel (ensemble)
   mean constructed with one realization of CMIP6 historical experiments for the
   period 1995-2014. Multimodel‐mean bias of (b) CMIP6 (1995-2014) compared to
   the corresponding time period of the climatology from ERA5
   (Copernicus Climate Change Service (C3S), 2017). (Fig. 3)

.. _fig_bock20jgr_4:
.. figure::  /recipes/figures/bock20jgr/ta850-global_to_swcre-global_RMSD.png
   :align:   center
   :width:   9cm

   Relative space-time root-mean-square deviation (RMSD) calculated from the 
   climatological seasonal cycle of the CMIP3, CMIP5, and CMIP6 simulations 
   (1980-1999) compared to observational data sets (Table 5). A relative 
   performance is displayed, with blue shading being better and red shading 
   worse than the median RMSD of all model results of all ensembles. A diagonal 
   split of a grid square shows the relative error with respect to the reference 
   data set (lower right triangle) and the alternative data set (upper left 
   triangle) which are marked in Table 5. White boxes are used when data are not 
   available for a given model and variable (Fig. 6).

.. _fig_bock20jgr_5:
.. figure::  /recipes/figures/bock20jgr/patterncor.png
   :align:   center
   :width:   9cm

   Centered pattern correlations between models and observations for the annual 
   mean climatology over the period 1980–1999 (Fig. 7). 
.. _nml_perfmetrics:

Performance metrics for essential climate parameters
====================================================

Overview
--------

The goal is to create a standard recipe for the calculation of performance metrics to quantify the ability of the models to reproduce the climatological mean annual cycle for selected "Essential Climate Variables" (ECVs) plus some additional corresponding diagnostics and plots to better understand and interpret the results.

The recipe can be used to calculate performance metrics at different vertical levels (e.g., 5, 30, 200, 850 hPa as in `Gleckler et al. (2008) <http://dx.doi.org/10.1029/2007JD008972>`_ and in different regions. As an additional reference, we consider `Righi et al. (2015) <https://doi.org/10.5194/gmd-8-733-2015>`_.

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_perfmetrics_CMIP5.yml
* recipe_perfmetrics_CMIP5_cds.yml
* recipe_perfmetrics_land_CMIP5.yml

Diagnostics are stored in diag_scripts/perfmetrics/

* main.ncl: calculates and (optionally) plots annual/seasonal cycles, zonal means, lat-lon fields and time-lat-lon fields. The calculated fields can also be plotted as difference w.r.t. a given reference dataset. main.ncl also calculates RMSD, bias and taylor metrics. Input data have to be regridded to a common grid in the preprocessor. Each plot type is created by a separated routine, as detailed below.
* cycle.ncl: creates an annual/seasonal cycle plot.
* zonal.ncl: creates a zonal (lat-pressure) plot.
* latlon.ncl: creates a lat-lon plot.
* cycle_latlon.ncl: precalculates the metrics for a time-lat-lon field, with different options for normalization.
* collect.ncl: collects and plots the metrics previously calculated by cycle_latlon.ncl.

User settings in recipe
-----------------------

.. _perf-main.ncl:

#. Script main.ncl

   *Required settings (scripts)*

   * plot_type: cycle (time), zonal (plev, lat), latlon (lat, lon), cycle_latlon (time, lat, lon), cycle_zonal (time, plev, lat)
   * time_avg: type of time average (monthlyclim, seasonalclim, annualclim)
   * region: selected region (global, trop, nhext, shext, nhtrop, shtrop, nh, sh, nhmidlat, shmidlat, nhpolar, shpolar, eq)

   *Optional settings (scripts)*

   * styleset: for plot_type cycle only (cmip5, righi15gmd, cmip6, default)
   * plot_stddev: for plot_type cycle only, plots standard deviation as shading
   * legend_outside: for plot_type cycle only, plots the legend in a separate file
   * t_test: for plot_type zonal or latlon, calculates t-test in difference plots (default: False)
   * conf_level: for plot_type zonal or latlon, adds the confidence level for the t-test to the plot (default: False)
   * projection: map projection for plot_type latlon (default: CylindricalEquidistant)
   * plot_diff: draws difference plots (default: False)
   * calc_grading: calculates grading metrics (default: False)
   * stippling: uses stippling to mark statistically significant differences (default: False = mask out non-significant differences in gray)
   * show_global_avg: diplays the global avaerage of the input field as string at the top-right of lat-lon plots (default: False)
   * metric: chosen grading metric(s) (if calc_grading is True)
   * normalization: metric normalization (for RMSD and BIAS metrics only)
   * abs_levs: list of contour levels for absolute plot
   * diff_levs: list of contour levels for difference plot
   * zonal_cmap: for plot_type zonal only, chosen color table (default: "amwg_blueyellowred")
   * zonal_ymin: for plot_type zonal only, minimum pressure level on the y-axis (default: 5. hPa)
   * latlon_cmap: for plot_type latlon only, chosen color table (default: "amwg_blueyellowred")
   * plot_units: plotting units (if different from standard CMOR units)

   *Required settings (variables)*

   * reference_dataset: reference dataset to compare with (usually the observations).

   *Optional settings (variables)*

   * alternative_dataset: a second dataset to compare with.

   These settings are passed to the other scripts by main.ncl, depending on the selected plot_type.

.. _perf-collect.ncl:

#. Script collect.ncl

   *Required settings (scripts)*

   * metric: selected metric (RMSD, BIAS or taylor)
   * label_bounds: for RMSD and BIAS metrics, min and max of the labelbar
   * label_scale: for RMSD and BIAS metrics, bin width of the labelbar
   * colormap: for RMSD and BIAS metrics, color table of the labelbar

   *Optional settings (scripts)*

   * label_lo: adds lower triange for values outside range
   * label_hi: adds upper triange for values outside range
   * cm_interval: min and max color of the color table
   * cm_reverse: reverses the color table
   * sort: sorts datasets in alphabetic order (excluding MMM)
   * diag_order: sort diagnostics in a specific order (name = 'diagnostic'-'region')
   * title: plots title
   * scale_font: scaling factor applied to the default font size
   * disp_values: switches on/off the grading values on the plot
   * disp_rankings: switches on/off the rankings on the plot
   * rank_order: displays rankings in increasing (1) or decreasing (-1) order

Variables
---------
#.  recipe_perfmetrics_CMIP5.yml

    * clt (atmos, monthly mean, longitude latitude time)
    * hus (atmos, monthly mean, longitude latitude lev time)
    * od550aer, od870aer, od550abs, od550lt1aer (aero, monthly mean, longitude latitude time)
    * pr (atmos, monthly mean, longitude latitude time)
    * rlut, rlutcs, rsut, rsutcs (atmos, monthly mean, longitude latitude time)
    * sm (land, monthly mean, longitude latitude time)
    * ta (atmos, monthly mean, longitude latitude lev time)
    * tas (atmos, monthly mean, longitude latitude time)
    * toz (atmos, monthly mean, longitude latitude time)
    * ts (atmos, monthly mean, longitude latitude time)
    * ua (atmos, monthly mean, longitude latitude lev time)
    * va (atmos, monthly mean, longitude latitude lev time)
    * zg (atmos, monthly mean, longitude latitude lev time)

#. recipe_perfmetrics_land_CMIP5.yml

    * sm (land, monthly mean, longitude latitude time)
    * nbp (land, monthly mean, longitude latitude time)
    * gpp (land, monthly mean, longitude latitude time)
    * lai (land, monthly mean, longitude latitude time)
    * fgco2 (ocean, monthly mean, longitude latitude time)
    * et (land, monthly mean, longitude latitude time)
    * rlus, rlds, rsus, rdsd (atmos, monthly mean, longitude latitude time)

Observations and reformat scripts
---------------------------------

The following list shows the currently used observational data sets for this recipe with their variable names and the reference to their respective reformat scripts in parentheses. Please note that obs4MIPs data can be used directly without any reformating. For non-obs4MIPs data see headers of cmorization scripts (in `/esmvaltool/cmorizers/obs/
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/>`_) for downloading and processing instructions.
#.  recipe_perfmetrics_CMIP5.yml

    * AIRS (hus - obs4MIPs)
    * CERES-EBAF (rlut, rlutcs, rsut, rsutcs - obs4MIPs)
    * ERA-Interim (tas, ta, ua, va, zg, hus - esmvaltool/cmorizers/obs/cmorize_obs_ERA-Interim.ncl)
    * ESACCI-AEROSOL (od550aer, od870aer, od550abs, od550lt1aer - esmvaltool/cmorizers/obs/cmorize_obs_ESACCI-AEROSOL.ncl)
    * ESACCI-CLOUD (clt - esmvaltool/cmorizers/obs/cmorize_obs_ESACCI-CLOUD.ncl)
    * ESACCI-OZONE (toz - esmvaltool/cmorizers/obs/cmorize_obs_ESACCI-OZONE.ncl)
    * ESACCI-SOILMOISTURE (sm - esmvaltool/cmorizers/obs/cmorize_obs_ESACCI-SOILMOISTURE.ncl)
    * ESACCI-SST (ts - esmvaltool/ucmorizers/obs/cmorize_obs_ESACCI-SST.ncl)
    * GPCP-SG (pr - obs4MIPs)
    * HadISST (ts - esmvaltool/cmorizers/obs/cmorize_obs_HadISST.ncl)
    * MODIS (od550aer - esmvaltool/cmorizers/obs/cmorize_obs_MODIS.ncl)
    * NCEP (tas, ta, ua, va, zg - esmvaltool/cmorizers/obs/cmorize_obs_NCEP.ncl)
    * NIWA-BS (toz - esmvaltool/cmorizers/obs/cmorize_obs_NIWA-BS.ncl)
    * PATMOS-x (clt - esmvaltool/cmorizers/obs/cmorize_obs_PATMOS-x.ncl)

#. recipe_perfmetrics_land_CMIP5.yml

    * CERES-EBAF (rlus, rlds, rsus, rsds - obs4MIPs)
    * ESACCI-SOILMOISTURE (sm - esmvaltool/cmorizers/obs/cmorize_obs_ESACCI-SOILMOISTURE.ncl)
    * FLUXCOM (gpp - esmvaltool/cmorizers/obs/cmorize_obs_fluxcom.py)
    * JMA-TRANSCOM (nbp, fgco2 - esmvaltool/cmorizers/obs/cmorize_obs_jma_transcom.py)
    * LAI3d (lai - esmvaltool/cmorizers/obs/cmorize_obs_lai3g.py)
    * LandFlux-EVAL (et - esmvaltool/cmorizers/obs/cmorize_obs_landflux_eval.py)
    * Landschuetzer2016 (fgco2 - esmvaltool/cmorizers/obs/cmorize_obs_landschuetzer2016.py)
    * MTE (gpp - esmvaltool/cmorizers/obs/cmorize_obs_mte.py)

References
----------

* Gleckler, P. J., K. E. Taylor, and C. Doutriaux, Performance metrics for climate models, J. Geophys. Res., 113, D06104, doi: 10.1029/2007JD008972 (2008).

* Righi, M., Eyring, V., Klinger, C., Frank, F., Gottschaldt, K.-D., Jöckel, P., and Cionni, I.: Quantitative evaluation of oone and selected climate parameters in a set of EMAC simulations, Geosci. Model Dev., 8, 733, doi: 10.5194/gmd-8-733-2015 (2015).

Example plots
-------------

.. figure:: /recipes/figures/perfmetrics/perfmetrics_fig_1.png
   :width: 90%

   Annual cycle of globally averaged temperature at 850 hPa (time period 1980-2005) for different CMIP5 models (historical simulation) (thin colored lines) in comparison to ERA-Interim (thick yellow line) and NCEP (thick black dashed line) reanalysis data.

.. figure:: /recipes/figures/perfmetrics/perfmetrics_fig_2.png
   :width: 90%

   Taylor diagram of globally averaged temperature at 850 hPa (ta) and longwave cloud radiative effect (lwcre) for different CMIP5 models (historical simulation, 1980-2005). Reference data (REF) are ERA-Interim for temperature (1980-2005) and CERES-EBAF (2001-2012) for longwave cloud radiative effect.

.. figure:: /recipes/figures/perfmetrics/perfmetrics_fig_3.png
   :width: 90%

   Difference in annual mean of zonally averaged temperature (time period 1980-2005) between the CMIP5 model MPI-ESM-MR (historical simulation) and ERA-Interim. Stippled areas indicdate differences that are statistically significant at a 95% confidence level.

.. figure:: /recipes/figures/perfmetrics/perfmetrics_fig_4.png
   :width: 90%

   Annual mean (2001-2012) of the shortwave cloud radiative effect from CERES-EBAF.

.. figure:: /recipes/figures/perfmetrics/perfmetrics_fig_5.png
   :width: 90%
   :align: center

   Relative space-time root-mean-square deviation (RMSD) calculated from the climatological seasonal cycle of CMIP5 simulations. A relative performance is displayed, with blue shading indicating better and red shading indicating worse performance than the median of all model results. A diagonal split of a grid square shows the relative error with respect to the reference data set (lower right triangle) and the alternative data set (upper left triangle). White boxes are used when data are not available for a given model and variable.
.. _recipes_eady_growth_rate:

Eady growth rate
================

Overview
--------

This recipe computes the maximum Eady Growth Rate and performs the annual and seasonal means, storing 
the results for each dataset. 
For the seasonal means, the results are plotted over the North-Atlantic region for the selected
pressure levels.


Available recipes and diagnostics
---------------------------------

Recipes are stored in ``esmvaltool/recipes/``

    * ``recipe_eady_growth_rate.yml``

Diagnostics are stored in ``esmvaltool/diag_scripts/eady_growth_rate/``

    * ``eady_growth_rate.py``: Computes and stores the eady growth rate. 
      Plots can be produced for the seasonal mean over the North Atlantic region.


User settings in recipe
-----------------------

#. Script ``eady_growth_rate.py``

   *Required settings for script*

   * ``time_statistic``: Set to `'annual'` to compute the annual mean. Set to `'seasonal'` to compute the seasonal mean.

   *Optional settings for script*

   * ``plot_levels``: list of pressure levels to be plotted for the seasonal mean. If not specified, all levels will be plotted.


Variables
---------

* ta (atmos, monthly mean, longitude latitude level time)
* zg (atmos, monthly mean, longitude latitude level time)
* ua (atmos, monthly mean, longitude latitude level time) 

References
----------
* Moreno-Chamarro, E., Caron, L-P., Ortega, P., Loosveldt Tomas, S., and Roberts, M. J., Can we trust CMIP5/6 future projections of European winter precipitation?. Environ. Res. Lett. 16 054063
* Brian J Hoskins and Paul J Valdes. On the existence of storm-tracks. Journal of the atmospheric sciences, 47(15):1854–1864, 1990.

Example plots
-------------

.. _fig_eady_growth_rate:
.. figure::  /recipes/figures/eady_growth_rate/HadGEM3-GC31-LM_winter_eady_growth_rate_70000.png 
   :align:   center

   Eady Growth Rate values over the North-Atlantic region at 70000 Pa.
.. _recipes_deangelis15nat:

Evaluate water vapor short wave radiance absorption schemes of ESMs with the observations.
==========================================================================================================================

Overview
--------


The recipe reproduces figures from `DeAngelis et al. (2015)`_:
Figure 1b to 4 from the main part as well as extended data figure 1 and 2.
This paper compares models with different schemes for water vapor short wave radiance absorption with the observations.
Schemes using pseudo-k-distributions with more than 20 exponential terms show the best results.

.. _`DeAngelis et al. (2015)`: https://www.nature.com/articles/nature15770


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_deangelis15nat.yml

Diagnostics are stored in diag_scripts/

   * deangelis15nat/deangelisf1b.py
   * deangelis15nat/deangelisf2ext.py
   * deangelis15nat/deangelisf3f4.py


User settings in recipe
-----------------------

The recipe can be run with different CMIP5 and CMIP6 models.
deangelisf1b.py:
Several flux variables (W m\ :sup:`-2`\) and up to 6 different model exeriements can be handeled.
Each variable needs to be given for each model experiment. The same experiments must
be given for all models.
In `DeAngelis et al. (2015)`_
150 year means are used but the recipe can handle any duration.

deangelisf2ext.py:

deangelisf3f4.py:
For each model, two experiments must be given:
a pre industrial control run, and a scenario with 4 times CO\ :sub:`2`\.
Possibly, 150 years should be given, but shorter time series work as well.


Variables
---------

deangelisf1b.py:
Tested for:

* *rsnst* (atmos, monthly, longitude, latitude, time)
* *rlnst* (atmos, monthly, longitude, latitude, time)
* *lvp* (atmos, monthly, longitude, latitude, time)
* *hfss* (atmos, monthly, longitude, latitude, time)

any flux variable (W m\ :sup:`-2`\) should be possible.

deangelisf2ext.py:

* *rsnst* (atmos, monthly, longitude, latitude, time)
* *rlnst* (atmos, monthly, longitude, latitude, time)
* *rsnstcs* (atmos, monthly, longitude, latitude, time)
* *rlnstcs* (atmos, monthly, longitude, latitude, time)
* *lvp* (atmos, monthly, longitude, latitude, time)
* *hfss* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)

deangelisf3f4.py:
* *rsnstcs* (atmos, monthly, longitude, latitude, time)
* *rsnstcsnorm* (atmos, monthly, longitude, latitude, time)
* *prw* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

deangelisf1b.py:
* None

deangelisf2ext.py:
* None

deangelisf3f4.py:

* *rsnstcs*:
   CERES-EBAF

* *prw*
   ERA-Interim, SSMI


References
----------

* DeAngelis, A. M., Qu, X., Zelinka, M. D., and Hall, A.: An observational radiative constraint on hydrologic cycle intensification, Nature, 528, 249, 2015.


Example plots
-------------


.. _bar_all:
.. figure:: /recipes/figures/deangelis15nat/bar_all.png
   :align: center
   :width: 50%

   Global average multi-model mean comparing different model experiments for the sum of upward long wave flux at TOA and net downward long wave flux at the surface (rlnst),  heating from short wave absorption (rsnst), latent heat release from precipitation (lvp), and sensible heat flux (hfss). The panel shows three model experiments, namely the pre-industrial control simulation averaged over 150 years (blue), the RCP8.5 scenario averaged over 2091-2100 (orange) and the abrupt quadrupled CO\ :sub:`2`\  scenario averaged over the years 141-150 after CO\ :sub:`2`\  quadrupling in all models except CNRM-CM5-2 and IPSL-CM5A-MR, where the average is calculated over the years 131-140 (gray). The figure shows that energy sources and sinks readjust in reply to an increase in greenhouse gases, leading to a decrease in the sensible heat flux and an increase in the other fluxes.

.. _exfig2a:
.. figure:: /recipes/figures/deangelis15nat/exfig2a.png
   :align: center
   :width: 50%

   The temperature-mediated response of each atmospheric energy budget term for each model as blue circles and the model mean as a red cross. The numbers above the abscissa are the cross-model correlations between dlvp/dtas and each other temperature-mediated response.'

.. _fig3b:
.. figure:: /recipes/figures/deangelis15nat/fig3b.png
   :align: center
   :width: 50%

   Scatter plot and regression line the between the ratio of the change of net short wave radiation (rsnst) and the change of the Water Vapor Path (prw) against the ratio of the change of netshort wave radiation for clear skye (rsnstcs) and the the change of surface temperature (tas). The width of horizontal shading for models and the vertical dashed lines for observations (Obs.) represent statistical uncertainties of the ratio, as the 95% confidence interval (CI) of the regression slope to the rsnst versus prw curve. For the observations the minimum of the lower bounds of all CIs to the maximum of the upper bounds of all CIs is shown.
.. _recipes_wenzel16jclim:

Multiple ensemble diagnostic regression (MDER) for constraining future austral jet position
===========================================================================================

Overview
--------

`Wenzel et al. (2016)`_ use multiple ensemble diagnostic regression (MDER) to
constrain the CMIP5 future projection of the summer austral jet position with
several historical process-oriented diagnostics and respective observations.

The following plots are reproduced:

* Absolute correlation between the target variable and the diagnostics.
* Scatterplot between the target variable and the MDER-calculated linear
  combination of diagnostics.
* Boxplot of RMSE for the unweighted multi-model mean and the (MDER) weighted
  multi-model mean of the target variable in a pseudo-reality setup.
* Time series of the target variable for all models, observations and MDER
  predictions.
* Errorbar plots for all diagnostics.
* Scatterplots between the target variable and all diagnostics.

.. _`Wenzel et al. (2016)`: https://journals.ametsoc.org/doi/full/10.1175/JCLI-D-15-0412.1


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_wenzel16jclim.yml


Diagnostics are stored in diag_scripts/

   * austral_jet/asr.ncl
   * austral_jet/main.ncl
   * mder/absolute_correlation.ncl
   * mder/regression_stepwise.ncl
   * mder/select_for_mder.ncl


User settings in recipe
-----------------------

#. Preprocessor

   * ``extract_region``: Region extraction.
   * ``extract_levels``: Pressure level extraction.
   * ``area_statistics``: Spatial average calculations.

#. Script austral_jet/asr.ncl

   * ``season``, *str*: Season.
   * ``average_ens``, *bool*, optional (default: ``False``): Average over all
     given ensemble members of a climate model.
   * ``wdiag``, *array of str*, optional: Names of the diagnostic for MDER
     output.  Necessary when MDER output is desired.
   * ``wdiag_title``, *array of str*, optional: Names of the diagnostic in
     plots.

#. Script austral_jet/main.ncl

   * ``styleset``, *str*: Style set used for plotting the multi-model plots.
   * ``season``, *str*: Season.
   * ``average_ens``, *bool*, optional (default: ``False``): Average over all
     given ensemble members of a climate model.
   * ``rsondes``, *array of str*, optional: Additional observations used in the
     plot but not for MDER output.
   * ``rsondes_file``, *array of str*, optional: Paths to the additional
     observations Necessary when ``rsondes`` is given.
   * ``rsondes_yr_min``, *int*, optional: Minimum year for additional
     observations. Necessary when ``rsondes`` is given.
   * ``rsondes_yr_max``, *int*, optional: Maximum year for additional
     observations. Necessary when ``rsondes`` is given.
   * ``wdiag``, *array of str*, optional: Names of the diagnostic for MDER
     output.  Necessary when MDER output is desired.
   * ``wdiag_title``, *array of str*, optional: Names of the diagnostic in
     plots.
   * ``derive_var``, *str*, optional: Derive variables using NCL functions.
     Must be one of ``"tpp"``, ``"mmstf"``.
   * ``derive_latrange``, *array of float*, optional: Latitude range for
     variable derivation.  Necessary if ``derive_var`` is given.
   * ``derive_lev``, *float*, optional: Pressure level (given in *Pa*) for
     variable derivation.  Necessary if ``derive_var`` is given.

#. Script mder/absolute_correlation.ncl

   * ``p_time``, *array of int*: Start years for future projections.
   * ``p_step``, *int*: Time range for future projections (in years).
   * ``scal_time``, *array of int*: Time range for base period (in years) for
     anomaly calculations used when ``calc_type = "trend"``.
   * ``time_oper``, *str*: Operation used in NCL ``time_operation`` function.
   * ``time_opt``, *str*: Option used in NCL ``time_operation`` function.
   * ``calc_type``, *str*: Calculation type for the target variable. Must be
     one of ``"trend"``, ``"pos"``, ``"int"``.
   * ``domain``, *str*: Domain tag for provenance tracking.
   * ``average_ens``, *bool*, optional (default: ``False``): Average over all
     given ensemble members of a climate model.
   * ``region``, *str*, optional: Region used for area aggregation. Necessary
     if input of target variable is multidimensional.
   * ``area_oper``, *str*, optional: Operation used in NCL ``area_operation``
     function. Necessary if multidimensional is given.
   * ``plot_units``, *str*, optional (attribute for ``variable_info``): Units
     for the target variable used in the plots.

#. Script mder/regression_stepwise.ncl

   * ``p_time``, *array of int*: Start years for future projections.
   * ``p_step``, *int*: Time range for future projections (in years).
   * ``scal_time``, *array of int*: Time range for base period (in years) for
     anomaly calculations used when ``calc_type = "trend"``.
   * ``time_oper``, *str*: Operation used in NCL ``time_operation`` function.
   * ``time_opt``, *str*: Option used in NCL ``time_operation`` function.
   * ``calc_type``, *str*: Calculation type for the target variable. Must be
     one of ``"trend"``, ``"pos"``, ``"int"``.
   * ``domain``, *str*: Domain tag for provenance tracking.
   * ``average_ens``, *bool*, optional (default: ``False``): Average over all
     given ensemble members of a climate model.
   * ``smooth``, *bool*, optional (default: ``False``): Smooth time period with
     1-2-1 filter.
   * ``iter``, *int*, optional: Number of iterations for smoothing. Necessary
     when ``smooth`` is given.
   * ``cross_validation_mode``, *bool*, optional (default: ``False``): Perform
     cross-validation.
   * ``region``, *str*, optional: Region used for area aggregation. Necessary
     if input of target variable is multidimensional.
   * ``area_oper``, *str*, optional: Operation used in NCL ``area_operation``
     function. Necessary if multidimensional is given.
   * ``plot_units``, *str*, optional (attribute for ``variable_info``): Units
     for the target variable used in the plots.

#. Script mder/select_for_mder.ncl

   * ``wdiag``, *array of str*: Names of the diagnostic for MDER output.
     Necessary when MDER output is desired.
   * ``domain``, *str*: Domain tag for provenance tracking.
   * ``ref_dataset``, *str*: Style set used for plotting the multi-model plots.
   * ``average_ens``, *bool*, optional (default: ``False``): Average over all
     given ensemble members of a climate model.
   * ``derive_var``, *str*, optional: Derive variables using NCL functions.
     Must be one of ``"tpp"``, ``"mmstf"``.


Variables
---------

* *ta* (atmos, monthly, longitude, latitude, pressure level, time)
* *uajet* (atmos, monthly, time)
* *va* (atmos, monthly, longitude, latitude, pressure level, time)
* *ps* (atmos, monthly, longitude, latitude, time)
* *asr* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* ERA-Intermin (*ta*, *uajet*, *va*, *ps*)
* CERES-EBAF (*asr*)


References
----------

* Wenzel, S., V. Eyring, E.P. Gerber, and A.Y. Karpechko: Constraining Future
  Summer Austral Jet Stream Positions in the CMIP5 Ensemble by Process-Oriented
  Multiple Diagnostic Regression. J. Climate, 29, 673–687,
  doi:10.1175/JCLI-D-15-0412.1, 2016.


Example plots
-------------

.. _fig_wenzel16jclim_1:
.. figure:: /recipes/figures/wenzel16jclim/CMPI5_uajet-pos_rcp45_20ystep_FIG1.png
   :align: center
   :width: 80%

   Time series of the the target variable (future austral jet position in the RCP
   4.5 scenario) for the CMIP5 ensemble, observations, unweighted multi-model mean
   projections and (MDER) weighted multi-model mean projections.

.. _fig_wenzel16jclim_2:
.. figure:: /recipes/figures/wenzel16jclim/CMPI5_uajet-pos_rcp45_20ystep_FIG2b.png
   :align: center
   :width: 80%

   Scatterplot of the target variable (future austral jet position in the RCP
   4.5 scenario) vs. the MDER-determined linear combination of diagnostics for the
   CMIP5 ensemble.

.. _fig_wenzel16jclim_3:
.. figure:: /recipes/figures/wenzel16jclim/CMPI5_uajet-pos_rcp45_20ystep_FIG3.png
   :align: center
   :width: 80%

   Boxplot for the RMSE of the target variable for the unweighted and (MDER)
   weighted multi-model mean projections in a pseudo-reality setup.

.. _fig_wenzel16jclim_4:
.. figure:: /recipes/figures/wenzel16jclim/ta_trop250_ta_DJF_trend.png
   :align: center
   :width: 80%

   Trends in tropical DJF temperature at 250hPa for different CMIP5 models and
   observations.

.. _fig_wenzel16jclim_5:
.. figure:: /recipes/figures/wenzel16jclim/uajet_H-SH_c.png
   :align: center
   :width: 80%

   Scatterplot of the target variable (future austral jet position in the RCP
   4.5 scenario) vs. a single diagnostic, the historical location of the
   Southern hemisphere Hadley cell boundary for the CMIP5 ensemble.
.. _recipes_cvdp:

Climate Variability Diagnostics Package (CVDP)
==============================================

Overview
--------
The Climate Variability Diagnostics Package (CVDP) developed by NCAR's Climate Analysis Section is an analysis tool that documents the major modes of climate variability in models and observations, including ENSO, Pacific Decadal Oscillation, Atlantic Multi-decadal Oscillation, Northern and Southern Annular Modes, North Atlantic Oscillation, Pacific North and South American teleconnection patterns. For details please refer to the [1] and [2].

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_cvdp.yml

Diagnostics are stored in diag_scripts/cvdp/

    * cvdp_wrapper.py

User settings in recipe
-----------------------

The recipe can be run with several data sets including different model ensembles, multi-model mean statistics are currently not supported.

Variables
---------

* ts (atmos, monthly mean, longitude latitude time)
* tas (atmos, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* psl (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

None. 

References
----------
[1] http://www.cesm.ucar.edu/working_groups/CVC/cvdp/

[2] https://github.com/NCAR/CVDP-ncl

Example plots
-------------

.. figure::  /recipes/figures/cvdp/nam.prreg.ann.png
   :align:   center

   Regression of the precipitation anomalies (PR) onto the Northern Annular 
   Mode (NAM) index for the time period 1900-2005 for 30 CMIP5 models and observations (GPCP (pr) / IFS-Cy31r2 (psl); time period 1984-2005).
.. _recipes_li17natcc:

Constraining future Indian Summer Monsoon projections with the present-day precipitation over the tropical western Pacific
==========================================================================================================================

Overview
--------


Following `Li et al. (2017)`_ the change between present-day and future Indian Summer Monsoon (ISM) precipitation is constrained
using the precipitation over the tropical western Pacific compared to
a fixed, observed amount of 6 mm d\ :sup:`-1` from Global Precipitation Climatology Project (GPCP) `(Adler et al., 2003)`_ for 1980-2009.
For CMIP6, historical data for 1980-2009 should be used. For CMIP5 historical data from 1980-2005 should be used, due to the length of the data sets.
At the moment it is not possible to use a combined ``['historical', 'rcp']`` data set, because the diagnostic requires that a historical data set is given.

.. _`(Adler et al., 2003)`: https://journals.ametsoc.org/doi/abs/10.1175/1525-7541%282003%29004%3C1147%3ATVGPCP%3E2.0.CO%3B2
.. _`Li et al. (2017)`: https://www.nature.com/articles/nclimate3387


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_li17natcc.yml


Diagnostics are stored in diag_scripts/

   * emergent_constraints/lif1f2.py


User settings in recipe
-----------------------

The recipe can be run with different CMIP5 and CMIP6 models. For each model, two experiments must be given: 
one historical run, possibly between 1980-2009 and one other model experiment. The user can choose the other model experiment, 
but it needs to be the same for all given models. 
The start and end year for the second data set can be choosen by the user, but should be consistent for all models 
(the same for future scenarios, the same length for other experiments). Different ensemble members are not possible, yet.


Variables
---------

* *pr* (atmos, monthly, longitude, latitude, time)
* *ua* (atmos, monthly, longitude, latitude, plev, time)
* *va* (atmos, monthly, longitude, latitude, plev, time)
* *ts* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Li, G., Xie, S. P., He, C., and Chen, Z. S.: Western Pacific emergent constraint lowers projected increase in Indian summer monsoon rainfall, Nat Clim Change, 7, 708-+, 2017


Example plots
-------------

.. _li17natcc_fig2a:
.. figure:: /recipes/figures/emergent_constraints/li17natcc_fig2a.png
   :align: center
   :width: 50%

   Scatter plot of the simulated tropical western Pacific precipitation (mm d\ :sup:`-1`\ ) versus projected average ISM (Indian Summer Monsoon) rainfall changes under the ssp585 scenario. The red line denotes the observed present-day western Pacific precipitation and the inter-model correlation (r) is shown. (CMIP6).

.. _li17natcc_fig2b:
.. figure:: /recipes/figures/emergent_constraints/li17natcc_fig2b.png
   :align: center
   :width: 50%

   Scatter plot of the uncorrected versus corrected average ISM (Indian Summer Monsoon) rainfall change ratios (% per degree Celsius of global SST warming). The error bars for the Multi-model mean indicate the standard deviation spread among models and the 2:1 line (y = 0.5x) is used to illustrate the Multi-model mean reduction in projected rainfall increase. (CMIP6).

.. _li17natcc_fig2c:
.. figure:: /recipes/figures/emergent_constraints/li17natcc_fig2c.png
   :align: center
   :width: 50%

   Multi-model mean rainfall change due to model error. Box displays the area used to define the average ISM (Indian Summer Monsoon) rainfall. Precipitation changes are normalized by the corresponding global mean SST increase for each model. (CMIP6).

.. _li17natcc_fig2d:
.. figure:: /recipes/figures/emergent_constraints/li17natcc_fig2d.png
   :align: center
   :width: 50%

   Corrected multi-model mean rainfall change. Box displays the area used to define the average ISM (Indian Summer Monsoon) rainfall. Precipitation changes are normalized by the corresponding global mean SST increase for each model. (CMIP6).
.. _recipe_kcs:

KNMI Climate Scenarios 2014
===========================

Overview
--------

This recipe implements the method described in `Lenderink et al., 2014 <https://doi.org/10.1088/1748-9326/9/11/115008>`_, to prepare the 2014 KNMI Climate Scenarios (KCS) for the Netherlands. A set of 8 global climate projections from EC-Earth were downscaled with the RACMO regional climate model. Since the EC-Earth ensemble is not readily representative for the spread in the full CMIP ensemble, this method recombines 5-year segments from the EC-Earth ensemble to obtain a large suite of "resamples". Subsequently, 8 new resamples are selected that cover the spread in CMIP much better than the original set.

The original method created 8 resampled datasets:

* 2 main scenarios: Moderate (M) and Warm (W) (Lenderink 2014 uses "G" instead of "M").
* 2 'sub'scenarios: Relatively high (H) or low (L) changes in seasonal temperature and precipitation
* 2 time horizons: Mid-century (MOC; 2050) and end-of-century (EOC; 2085)
* Each scenario consists of changes calculated between 2 periods: Control (1981-2010) and future (variable).

The configuration settings for these resamples can be found in table 1 of Lenderink 2014's `supplementary data <https://iopscience.iop.org/1748-9326/9/11/115008/media/erl503687suppdata.pdf>`_.

Implementation
--------------

The implementation is such that application to other datasets, regions, etc. is relatively straightforward. The description below focuses on the reference use case of Lenderink et al., 2014, where the target model was EC-Earth. An external set of EC-Earth data (all RCP85) was used, for which 3D fields for downscaling were available as well. In the recipe shipped with ESMValTool, however, the target model is CCSM4, so that it works out of the box with ESGF data only.

In the first diagnostic, the spread of the full CMIP ensemble is used to obtain 4 values of a *global* :math:`{\Delta}T_{CMIP}`, corresponding to the 10th and 90th percentiles for the M and W scenarios, respectively, for both MOC and EOC. Subsequently, for each of these 4 *steering parameters*, 30-year periods are selected from the target model ensemble, where :math:`{\Delta}T_{target}{\approx}{\Delta}T_{CMIP}`.

In the second diagnostic, for both the control and future periods, the N target model ensemble members are split into 6 segments of 5 years each. Out of all :math:`N^6` possible re-combinations of these 5-year segments, eventually M new 'resamples' are selected based on *local* changes in seasonal temperature and precipitation. This is done in the following steps:

1. Select 1000 samples for the control period, and 2 x 1000 samples for the future period (one for each subscenario). Step 1 poses a constraint on winter precipitation. For the control period, winter precipitation must still closely represent the average of the original ensemble. For the two future periods, the change in winter precipitation with respect to the control period must approximately equal 4% per degree :math:`{\Delta}T` (subscenario L) or 8% per degree :math:`{\Delta}T` (subscenario H).
2. Further constrain the selection by picking samples that represent either high or low changes in summer precipitation and summer and winter temperature, by limiting the remaining samples to certain percentile ranges: relatively wet/cold in the control and dry/warm in the future, or vice versa. The percentile ranges are listed in table 1 of Lenderink 2014's supplement. This should result is approximately 50 remaining samples for each scenario, for both control and future.
3. Use a Monte-Carlo method to make a final selection of 8 resamples with minimal reuse of the same ensemble member/segment.

Datasets have been split in two parts: the CMIP datasets and the target model datasets. An example use case for this recipe is to compare between CMIP5 and CMIP6, for example. The recipe can work with a target model that is not part of CMIP, provided that the data are CMOR compatible, and using the same data referece syntax as the CMIP data. Note that you can specify :ref:`multiple data paths<config-user-rootpath>` in the user configuration file.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

- recipe_kcs.yml

Diagnostics are stored in diag_scripts/kcs/

- global_matching.py
- local_resampling.py

.. note::
    We highly recommend using the options described in :ref:`rerunning`. The speed bottleneck for the first diagnostic is the preprocessor. In the second diagnostic, step 1 is most time consuming, whereas steps 2 and 3 are likely to be repeated several times. Therefore, intermediate files are saved after step 1, and the diagnostic will automatically detect and use them if the ``-i`` flag is used.

User settings
-------------

1. Script <global_matching.py>

  *Required settings for script*

  * ``scenario_years``: a list of time horizons. Default: ``[2050, 2085]``
  * ``scenario_percentiles``: a list of percentiles for the steering table. Default: ``[p10, p90]``

  *Required settings for preprocessor*
  This diagnostic needs global mean temperature anomalies for each dataset, both CMIP and the target model. Additionally, the multimodel statistics preprocessor must be used to produce the percentiles specified in the setting for the script above.

2. Script <local_resampling.py>

  *Required settings for script*

  * ``control_period``: the control period shared between all scenarios. Default: ``[1981, 2010]``
  * ``n_samples``: the final number of recombinations to be selected. Default: ``8``
  * ``scenarios``: a scenario name and list of options. The default setting is a single scenario:

    .. code-block:: yaml

        scenarios:
          ML_MOC:  # scenario name; can be chosen by the user
            description: "Moderate / low changes in seasonal temperature & precipitation"
            global_dT: 1.0
            scenario_year: 2050
            resampling_period: [2021, 2050]
            dpr_winter: 4
            pr_summer_control: [25, 55]
            pr_summer_future: [45, 75]
            tas_winter_control: [50, 80]
            tas_winter_future: [20, 50]
            tas_summer_control: [0, 100]
            tas_summer_future: [0, 50]

    These values are taken from table 1 in the Lenderink 2014's supplementary material. Multiple scenarios can be processed at once by appending more configurations below the default one. For new applications, ``global_dT``, ``resampling_period`` and ``dpr_winter`` are informed by the output of the first diagnostic. The percentile bounds in the scenario settings (e.g. ``tas_winter_control`` and ``tas_winter_future``) are to be tuned until a satisfactory scenario spread over the full CMIP ensemble is achieved.

  *Required settings for preprocessor*

  This diagnostic requires data on a single point. However, the ``extract_point`` preprocessor can be changed to ``extract_shape`` or ``extract_region``, in conjunction with an area mean. And of course, the coordinates can be changed to analyze a different region.

Variables
---------

Variables are precipitation and temperature, specified separately for the target model and the CMIP ensemble:

* pr_target (atmos, monthly mean, longitude latitude time)
* tas_target (atmos, monthly mean, longitude latitude time)
* pr_cmip (atmos, monthly mean, longitude latitude time)
* tas_cmip (atmos, monthly mean, longitude latitude time)

References
----------

* `Lenderink et al. 2014, Environ. Res. Lett., 9, 115008 <https://doi.org/10.1088/1748-9326/9/11/115008>`_.

Example output
--------------

The diagnostic ``global_matching`` produces a scenarios table like the one below

.. code-block:: python

       year percentile  cmip_dt period_bounds  target_dt  pattern_scaling_factor
    0  2050        P10     0.98  [2019, 2048]       0.99                    1.00
    1  2050        P90     2.01  [2045, 2074]       2.02                    0.99
    2  2085        P10     1.38  [2030, 2059]       1.38                    1.00
    3  2085        P90     3.89  [2071, 2100]       3.28                    1.18


which is printed to the log file and also saved as a csv-file ``scenarios.csv``.
Additionally, a figure is created showing the CMIP spread in global temperature change,
AND highlighting the selected steering parameters and resampling periods:

.. _fig_kcs_global_matching:
.. figure::  /recipes/figures/kcs/global_matching.png
   :align:   center

The diagnostic ``local_resampling`` procudes a number of output files:

* ``season_means_<scenario>.nc``: intermediate results, containing the season means for each segment of the original target model ensemble.
* ``top1000_<scenario>.csv``: intermediate results, containing the 1000 combinations that have been selected based on winter mean precipitation.
* ``indices_<scenario>.csv``: showing the final set of resamples as a table:

  .. code-block:: python

                      control                                                      future
                    Segment 0 Segment 1 Segment 2 Segment 3 Segment 4 Segment 5 Segment 0 Segment 1 Segment 2 Segment 3 Segment 4 Segment 5
     Combination 0          5         7         6         3         1         3         2         4         2         4         7         7
     Combination 1          0         3         0         4         3         2         4         1         6         1         3         0
     Combination 2          2         4         3         7         4         2         5         4         6         6         4         2
     Combination 3          1         4         7         2         3         6         5         3         1         7         4         1
     Combination 4          5         7         6         3         1         3         2         3         0         6         1         7
     Combination 5          7         2         1         4         5         1         6         0         4         2         3         3
     Combination 6          7         2         2         0         6         6         5         2         1         5         4         2
     Combination 7          6         3         2         1         6         1         2         1         0         2         1         3


* Provenance information: bibtex, xml, and/or text files containing citation information are stored alongside the final result and the final figure.
  The final combinations only derive from the target model data, whereas the figure also uses CMIP data.
* A figure used to validate the final result, reproducing figures 5 and 6 from Lenderink et al.:

.. _fig_kcs_local_validation:
.. figure::  /recipes/figures/kcs/local_validation_2085.png
   :align:   center
.. _recipes_modes_of_variability:

Modes of variability
====================

Overview
--------

The goal of this recipe is to compute modes of variability from a reference or observational dataset and from a set of climate projections and calculate the root-mean-square error between the mean anomalies obtained for the clusters from the reference and projection data sets.
This is done through K-means or hierarchical clustering applied either directly to the spatial data or after computing the EOFs.

The user can specify the number of clusters to be computed.

The recipe's output consist of three netcdf files for both the observed and projected weather regimes and the RMSE between them.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_modes_of_variability.yml


Diagnostics are stored in diag_scripts/magic_bsc/

* WeatherRegime.R - function for computing the EOFs and k-means and hierarchical clusters.

* weather_regime.R - applies the above weather regimes function to the datasets



User settings
-------------

User setting files are stored in recipes/

#. recipe_modes_of_variability.yml

   *Required settings for script*

   * plot type: rectangular or polar
   * ncenters: number of centers to be computed by the clustering algorithm (maximum 4)
   * cluster_method: kmeans (only psl variable) or hierarchical clustering (for psl or sic variables) 
   * detrend_order: the order of the polynomial detrending to be applied (0, 1 or 2)
   * EOFs: logical indicating wether the k-means clustering algorithm is applied directly to the spatial data ('false') or to the EOFs ('true')
   * frequency: select the month (format: JAN, FEB, ...) or season (format: JJA, SON, MAM, DJF) for the diagnostic to be computed for (does not work yet for MAM with daily data).


Variables
---------

* psl (atmos, monthly/daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Dawson, A., T. N. Palmer, and S. Corti, 2012: Simulating regime structures in weather and climate prediction models. Geophysical Research Letters, 39 (21), https://doi.org/10.1029/2012GL053284.

* Ferranti, L., S. Corti, and M. Janousek, 2015: Flow-dependent verification of the ECMWF ensemble over the Euro-Atlantic sector. Quarterly Journal of the Royal Meteorological Society, 141 (688), 916-924, https://doi.org/10.1002/qj.2411.

* Grams, C. M., Beerli, R., Pfenninger, S., Staffell, I., & Wernli, H. (2017). Balancing Europe's wind-power output through spatial deployment informed by weather regimes. Nature climate change, 7(8), 557, https://doi.org/10.1038/nclimate3338.

* Hannachi, A., D. M. Straus, C. L. E. Franzke, S. Corti, and T. Woollings, 2017: Low Frequency Nonlinearity and Regime Behavior in the Northern Hemisphere Extra-Tropical Atmosphere. Reviews of Geophysics, https://doi.org/10.1002/2015RG000509.

* Michelangeli, P.-A., R. Vautard, and B. Legras, 1995: Weather regimes: Recurrence and quasi stationarity. Journal of the atmospheric sciences, 52 (8), 1237-1256, doi: `10.1175/1520-0469(1995)052<1237:WRRAQS>2.0.CO <https://journals.ametsoc.org/doi/10.1175/1520-0469%281995%29052%3C1237%3AWRRAQS%3E2.0.CO%3B2>`_. 

* Vautard, R., 1990: Multiple weather regimes over the North Atlantic: Analysis of precursors and successors. Monthly weather review, 118 (10), 2056-2081, doi: `10.1175/1520-0493(1990)118<2056:MWROTN>2.0.CO;2 <https://journals.ametsoc.org/doi/10.1175/1520-0493%281990%29118%3C2056%3AMWROTN%3E2.0.CO%3B2>`_.

* Yiou, P., K. Goubanova, Z. X. Li, and M. Nogaj, 2008: Weather regime dependence of extreme value statistics for summer temperature and precipitation. Nonlinear Processes in Geophysics, 15 (3), 365-378, https://doi.org/10.5194/npg-15-365-2008.




Example plots
-------------

.. _fig_modesofvar:
.. figure::  /recipes/figures/modes_of_variability/SON-psl_predicted_regimes.png
   :align:   center
   :width:   14cm

Four modes of variability for autumn (September-October-November) in the North Atlantic European Sector for the RCP 8.5 scenario using BCC-CSM1-1 future projection during the period 2020-2075. The frequency of occurrence of each variability mode is indicated in the title of each map.


.. _recipes_albedolandcover:

Landcover - Albedo
==================


Overview
--------

The diagnostic determines the coefficients of multiple linear regressions fitted between the albedo values and the tree, shrub, short vegetation (crops and grasses) fractions of each grid cell within spatially moving windows encompassing 5x5 model grid cells. Solving these regressions provides the albedo values for trees, shrubs and short vegetation (crops and grasses) from which the albedo changes associated with transitions between these three landcover types are derived. The diagnostic distinguishes between snow-free and snow-covered grid cells.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_albedolandcover.yml

Diagnostics are stored in diag_scripts/landcover/

    * albedolandcover.py


User settings
-------------

Several parameters can be set in the recipe


Variables
---------

* rsus (atmos, monthly mean, time latitude longitude)
* rsds (atmos, monthly mean, time latitude longitude)
* snc (landice, monthly mean, time latitude longitude)
* grassFrac (land, monthly mean, time latitude longitude)
* treeFrac (land, monthly mean, time latitude longitude)
* shrubFrac (land, monthly mean, time latitude longitude)
* cropFrac (land, monthly mean, time latitude longitude)
* pastureFrac (land, monthly mean, time latitude longitude)


Observations and reformat scripts
---------------------------------

A reformatting script for observational data is available here:
    * cmorize_obs_duveiller2018.py


References
----------

* Duveiller, G., Hooker, J. and Cescatti, A., 2018a. A dataset mapping the potential biophysical effects of vegetation cover change. Scientific Data, 5: 180014.

* Duveiller, G., Hooker, J. and Cescatti, A., 2018b. The mark of vegetation change on Earth’s surface energy balance. Nature communications, 9(1): 679.

Example plots
-------------

.. _fig_landcoveralbedo_CMIP5_MPI-ESM-LR:
.. figure::  /recipes/figures/albedolandcover/MPI-ESM-LR_albedo_change_from_tree_to_crop-grass.png
   :align:   center
   :width:   14cm

   Example of albedo change from tree to crop and grass for the CMIP5 model MPI-ESM-LR derived for the month of July and averaged over the years 2000 to 2004.
.. _recipes_anav13jclim:

Land and ocean components of the global carbon cycle
====================================================

Overview
--------

This recipe reproduces most of the figures of `Anav et al. (2013)`_:

* Timeseries plot for different regions
* Seasonal cycle plot for different regions
* Errorbar plot for different regions showing mean and standard deviation
* Scatterplot for different regions showing mean vs. interannual variability
* 3D-scatterplot for different regions showing mean vs. linear trend and the
  model variability index (MVI) as a third dimension (color coded)
* Scatterplot for different regions comparing two variable against each other
  (*cSoil* vs. *cVeg*)

In addition, performance metrics are calculated for all variables using the
performance metric diagnostics (see details in :ref:`nml_perfmetrics`).


.. _mvi calculation:

MVI calculation
---------------

The Model variability index (MVI) on a single grid point (calculated in
``carbon_cycle/mvi.ncl`` is defined as

.. math::

   MVI = \left( \frac{s^M}{s^O} - \frac{s^O}{s^M} \right)^2

where :math:`s^M` and :math:`s^O` are the standard deviations of the annual
time series on a single grid point of a climate model :math:`M` and the
reference observation :math:`O`. In order to get a global or regional result,
this index is simple averaged over the respective domain.

In its given form, this equation is prone to small standard deviations close to
zero. For example, values of :math:`s^M = 10^{-5} \mu` and :math:`s^O = 10^{-7}
\mu` (where :math:`\mu` is the mean of :math:`s^O` over all grid cells) results
in a MVI of the order of :math:`10^4` for this single grid cell even
though the two standard deviations are close to zero and negligible compared to
other grid cells. Due to the use of the arithmetic mean, a single high value is
able to distort the overall MVI.

In the original publication, the maximum MVI is in the order of 10 (for the
variable `gpp`). However, a naive application of the MVI definition yields
values over :math:`10^9` for some models. Unfortunately, `Anav et al. (2013)`_
do not provide an explanation on how to deal with this problem. Nevertheless,
this script provides two configuration options to avoid high MVI values, but
they are not related to the original paper or any other peer-revied study and
should be used with great caution (see :ref:`user settings`).

.. _`Anav et al. (2013)`: https://journals.ametsoc.org/doi/full/10.1175/JCLI-D-12-00417.1


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_anav13jclim.yml


Diagnostics are stored in diag_scripts/

   * carbon_cycle/main.ncl
   * carbon_cycle/mvi.ncl
   * carbon_cycle/two_variables.ncl
   * perfmetrics/main.ncl
   * perfmetrics/collect.ncl


.. _user settings:

User settings in recipe
-----------------------

#. Preprocessor

   * ``mask_fillvalues``: Mask common missing values on different datasets.
   * ``mask_landsea``: Mask land/ocean.
   * ``regrid``: Regridding.
   * ``weighting_landsea_fraction``: Land/ocean fraction weighting.

#. Script carbon_cycle/main.ncl

   * ``region``, *str*: Region to be averaged.
   * ``legend_outside``, *bool*: Plot legend in a separate file (does not
     affect errorbar plot and evolution plot)
   * ``seasonal_cycle_plot``, *bool*: Draw seasonal cycle plot.
   * ``errorbar_plot``, *bool*: Draw errorbar plot.
   * ``mean_IAV_plot``, *bool*: Draw Mean (x-axis), IAV (y-axis) plot.
   * ``evolution_plot``, *bool*: Draw time evolution of a variable comparing
     a reference dataset to multi-dataset mean; requires ref_dataset in recipe.
   * ``sort``, *bool*, optional (default: ``False``): Sort dataset in
     alphabetical order.
   * ``anav_month``, *bool*, optional (default: ``False``): Conversion of
     y-axis to PgC/month instead of /year.
   * ``evolution_plot_ref_dataset``, *str*, optional: Reference dataset for
     evolution_plot. Required when ``evolution_plot`` is ``True``.
   * ``evolution_plot_anomaly``, *str*, optional (default: ``False``): Plot
     anomalies in evolution plot.
   * ``evolution_plot_ignore``, *list*, optional: Datasets to ignore in
     evolution plot.
   * ``evolution_plot_volcanoes``, *bool*, optional (default: ``False``): Turns
     on/off lines of volcano eruptions in evolution plot.
   * ``evolution_plot_color``, *int*, optional (default: ``0``): Hue of the
     contours in the evolution plot.
   * ``ensemble_name``, *string*, optional: Name of ensemble for use in evolution plot legend

#. Script carbon_cycle/mvi.ncl

   * ``region``, *str*: Region to be averaged.
   * ``reference_dataset``, *str*: Reference dataset for the MVI calculation
     specified for each variable seperately.
   * ``mean_time_range``, *list*, optional: Time period over which the mean is
     calculated (if not given, use whole time span).
   * ``trend_time_range``, *list*, optional: Time period over which the trend
     is calculated (if not given, use whole time span).
   * ``mvi_time_range``, *list*, optional: Time period over which the MVI is
     calculated (if not given, use whole time span).
   * ``stddev_threshold``, *float*, optional (default: ``1e-2``): Threshold to
     ignore low standard deviations (relative to the mean) in the MVI
     calculations. See also :ref:`mvi calculation`.
   * ``mask_below``, *float*, optional: Threshold to mask low absolute values
     (relative to the mean) in the input data (not used by default). See also
     :ref:`mvi calculation`.

#. Script carbon_cycle/two_variables.ncl

   * ``region``, *str*: Region to be averaged.

#. Script perfmetrics/main.ncl

   See :ref:`nml_perfmetrics`.

#. Script perfmetrics/collect.ncl

   See :ref:`nml_perfmetrics`.


Variables
---------

* *tas* (atmos, monthly, longitude, latitude, time)
* *pr* (atmos, monthly, longitude, latitude, time)
* *nbp* (land, monthly, longitude, latitude, time)
* *gpp* (land, monthly, longitude, latitude, time)
* *lai* (land, monthly, longitude, latitude, time)
* *cveg* (land, monthly, longitude, latitude, time)
* *csoil* (land, monthly, longitude, latitude, time)
* *tos* (ocean, monthly, longitude, latitude, time)
* *fgco2* (ocean, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* CRU (*tas*, *pr*)
* JMA-TRANSCOM (*nbp*, *fgco2*)
* MTE (*gpp*)
* LAI3g (*lai*)
* NDP (*cveg*)
* HWSD (*csoil*)
* HadISST (*tos*)


References
----------

* Anav, A. et al.: Evaluating the land and ocean components of the global
  carbon cycle in the CMIP5 Earth System Models, J. Climate, 26, 6901-6843,
  doi: 10.1175/JCLI-D-12-00417.1, 2013.


Example plots
-------------

.. _fig_anav13jclim_1:
.. figure:: /recipes/figures/anav13jclim/nbp_evolution_global.png
   :align: center
   :width: 80%

   Time series of global net biome productivity (NBP) over the period
   1901-2005. Similar to Anav et al.  (2013), Figure 5.

.. _fig_anav13jclim_2:
.. figure:: /recipes/figures/anav13jclim/gpp_cycle_nh.png
   :align: center
   :width: 80%

   Seasonal cycle plot for nothern hemisphere gross primary production (GPP)
   over the period 1986-2005. Similar to Anav et al. (2013), Figure 9.

.. _fig_anav13jclim_3:
.. figure:: /recipes/figures/anav13jclim/gpp_errorbar_trop.png
   :align: center
   :width: 80%

   Errorbar plot for tropical gross primary production (GPP) over the period
   1986-2005.

.. _fig_anav13jclim_4:
.. figure:: /recipes/figures/anav13jclim/tos_scatter_global.png
   :align: center
   :width: 80%

   Scatterplot for interannual variability and mean of global sea surface
   temperature (TOS) over the period 1986-2005.

.. _fig_anav13jclim_5:
.. figure:: /recipes/figures/anav13jclim/tas_global.png
   :align: center
   :width: 80%

   Scatterplot for multiyear average of 2m surface temperature (TAS) in x axis,
   its linear trend in y axis, and MVI. Similar to Anav et al. (2013) Figure 1
   (bottom).

.. _fig_anav13jclim_6:
.. figure:: /recipes/figures/anav13jclim/cSoil-cVeg_scatter_global.png
   :align: center
   :width: 80%

   Scatterplot for vegetation carbon content (cVeg) and soil carbon content
   (cSoil) over the period 1986-2005. Similar to Anav et al. (2013), Figure 12.

.. _fig_anav13jclim_7:
.. figure:: /recipes/figures/anav13jclim/diag_grading_pr-global_to_diag_grading_gpp-global_RMSD.png
   :align: center
   :width: 80%

   Performance metrics plot for carbon-cycle-relevant diagnostics.
.. _recipes_gier20bg:

Spatially resolved evaluation of ESMs with satellite column-averaged CO\ :sub:`2`
=================================================================================

Overview
--------

This recipe reproduces the figures of Gier et al. (2020). It uses satellite
column-averaged CO\ :sub:`2` data to evaluate ESMs by plotting several
quantities such as timeseries, seasonal cycle and growth rate in different
areas.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_gier20bg.yml

Diagnostics are stored in diag_scripts/

Diagnostics are stored in esmvaltool/diag_scripts/xco2_analysis/

    * carbon_plots.ncl: plot script for panel plots
    * delta_T.ncl: IAV of growth rate against growing season temperature - Figure C1
    * global_maps.ncl: global maps for seasonal cycle amplitude - Figures 5, 6
    * main.ncl: Timeseries and histogram - Figures 3, 4
    * panel_plots.ncl: scatter plot of SCA/GR vs variable - Figures 7, 9, B1, B2
    * sat_masks.ncl: data coverage of input data - Figures 1, 8
    * stat.ncl: auxiliary functions for GR, SCA computation
    * station_comparison.ncl: - comparison of surface and column data - Figure 2


User settings in recipe
-----------------------

#. Preprocessor

    * ``conv_units``: converts units to plot-units
    * ``mmm_ref``: calculates multi-model mean and regrids to ref dataset
    * ``mmm_2x2``: computes multi-model mean on 2x2 grid
    * ``mmm``: computes multi-model mean for 3D variable, 5x5 grid with specific
      pressure levels

#. Script xco2_analysis/delta_T.ncl

    * Required diag_script_info attributes:
        * ``region``: region to average over
        * ``masking``: the kind of masking to apply prior to region average
          (possible options: obs, land, sciamachy, gosat, none)
        * ``var_order``: First main variable, then temperature variable to compare

    * Optional diag_script_info attributes:
        * ``styleset``: styleset for color coding panels
        * ``output_file_type``: output file type for plots, default: config_user -> png
        * ``var_plotname``: NCL string formatting how variable should be named in plots
          defaults to short_name if not assigned.

#. Script xco2_analysis/global_maps.ncl:

    * Required diag_script_info attributes:
        * ``contour_max_level``: maximum value displayed for seasonal cycle
          amplitude contour plot

    * Optional diag_script_info attributes:
        * ``output_file_type``: output file type for plots, default: config_user -> png

#. Script xco2_analysis/main.ncl:

    * Required diag_script_info attributes:
        * ``styleset``: styleset to use for plotting colors, linestyles...
        * ``region``: latitude range for averaging
        * ``masking``: different masking options are available to use on dataset:
          (possible options: none, obs)
        * ``ensemble_mean``: if true calculates multi-model mean only
          accounting for the ensemble member named in "ensemble_refs"

    * Optional diag_script_info attributes:
        * ``output_file_type``: output file type for plots, default: config_user -> png
        * ``ensemble_refs``: list of model-ensemble pairs to denote which ensemble
          member to use for calculating multi-model mean. required if
          ensemble_mean = true
        * ``var_plotname``: String formatting how variable should be named in plots
          defaults to short_name if not assigned

#. Script xco2_analysis/panel_plots.ncl:

    * Required diag_script_info attributes:
        * ``styleset``: styleset to use for plotting colors, linestyles...
        * ``region``: latitude range for averaging
        * ``masking``: different masking options are available to use on dataset:
          (possible options: obs, land, sciamachy, gosat, none)
        * ``obs_in_panel``: True if observations should be included in plot
        * ``area_avg``: Type of area averaging: "full-area" normal area-average
          "lat-first" calculate zonal means first, then average these
        * ``plot_var2_mean``: If True adds mean of seasonal cycle to panel as string.

    * Optional diag_script_info attributes:
        * ``output_file_type``: output file type for plots, default: config_user -> png
        * ``var_plotname``: String formatting how variable should be named in plots
          defaults to short_name if not assigned

#. Script xco2_analysis/sat_masks.ncl:

    * Optional diag_script_info attributes:
        * ``output_file_type``: output file type for plots, default: config_user -> png
        * ``var_plotname``: String formatting how variable should be named in plots
          defaults to short_name if not assigned
        * ``c3s_plots``: Missing value plots seperated by timeseries of c3s satellites

#. Script xco2_analysis/station_comparison.ncl:

    * Required diag_script_info attributes:
        * ``var_order``: in this case xco2, co2, co2s - column averaged with obs dataset
          first, then 2D variable, followed by surface stations

    * Optional diag_script_info attributes:
        * ``output_file_type``: output file type for plots, default: config_user -> png
        * ``var_plotnames``: String formatting how variables should be named in plots
          defaults to short_name if not assigned
        * ``overwrite_altitudes``: Give other altitude values than the ones attached in
          the station data. Valid if altitude changes and
          timeseries spans range with different sample
          altitude. Caveat: If used, need to give altitude
          values for all stations.
        * ``output_map``: boolean if stations to be displayed on map. As this requires
          finetuning, currently only implemented for station set of
          (ASK, CGO, HUN, LEF, WIS) following the paper. Change for different
          plot inset locations, if others are desired.

Variables
---------

* *xco2* (atmos, monthly, longitude, latitude, time)
* *co2s* (atmos, monthly, longitude, latitude, time)
* *co2* (atmos, monthly, pressure, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)
* *tasa* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* CDS-XCO2_ (*xco2*)
* ESRL_ (*co2s*)
* GISTEMP_ (*tasa*)
* MODIS_ (land cover map, auxiliary data folder)

.. _ESRL: https://www.esrl.noaa.gov/gmd/dv/data/
.. _GISTEMP: https://data.giss.nasa.gov/gistemp/
.. _CDS-XCO2: https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-carbon-dioxide?tab=form
.. _MODIS: https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=968

References
----------

* Gier, B. K., Buchwitz, M., Reuter, M., Cox, P. M., Friedlingstein, P.,
  and Eyring, V.: Spatially resolved evaluation of Earth system models with
  satellite column-averaged CO2, Biogeosciences, 17, 6115–6144,
  https://doi.org/10.5194/bg-17-6115-2020, 2020.

Example plots
-------------

.. _fig_gier20bg_1:
.. figure::  /recipes/figures/gier20bg/fig01.png
   :align:   center
   :width: 80%

   Mean fractional coverage of monthly satellite data.

.. _fig_gier20bg_2:
.. figure::  /recipes/figures/gier20bg/fig02.png
   :align:   center
   :width: 80%

   Comparison of time series from satellite, in situ, and models sampled
   accordingly. Caveat: inset plot positions are hardcoded.

.. _fig_gier20bg_3:
.. figure::  /recipes/figures/gier20bg/fig03.png
   :align:   center
   :width: 70%

   Timeseries with panels depicting growth rate and seasonal cycle.

.. _fig_gier20bg_4:
.. figure::  /recipes/figures/gier20bg/fig04.png
   :align:   center
   :width: 50%

   Barplot of the growth rate, averaged over all years, with standard deviation
   of interannual variability.

.. _fig_gier20bg_5:
.. figure::  /recipes/figures/gier20bg/fig05.png
   :align:   center
   :width: 80%

   Panel plot of spatially resolved seasonal cycle amplitude for all models,
   including a zonal average sidepanel.

.. _fig_gier20bg_6:
.. figure::  /recipes/figures/gier20bg/fig06.png
   :align:   center
   :width: 60%

   Seasonal cycle amplitude map comparing influence of sampling, and difference
   to observations.

.. _fig_gier20bg_7:
.. figure::  /recipes/figures/gier20bg/fig07.png
   :align:   center
   :width: 50%

   Panel plots showing seasonal cycle amplitude against XCO\ :sub:`2`, includes
   regression line and p-value.

.. _fig_gier20bg_8:
.. figure::  /recipes/figures/gier20bg/fig08.png
   :align:   center
   :width: 50%

   Mean spatial data coverage for different satellites.
.. _nml_seaice:

Sea Ice
=======

Overview
--------
The sea ice diagnostics include:

(1) time series of Arctic and Antarctic sea ice area and extent
    (calculated as the total area (km\ :sup:`2`\) of grid cells with sea ice concentrations
    (sic) of at least 15%).
(2) ice extent trend distributions for the Arctic in September and the Antarctic in February.
(3) calculation of year of near disappearance of Arctic sea ice
(4) scatter plots of (a) historical trend in September Arctic sea ice extent (SSIE) vs
    historical long-term mean SSIE; (b) historical SSIE mean vs 1st year of disappearance
    (YOD) RCP8.5; (c) historical SSIE trend vs YOD RCP8.5.

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_seaice.yml

Diagnostics are stored in diag_scripts/seaice/

* seaice_aux.ncl: contains functions for calculating sea ice area or extent from sea ice
  concentration and first year of disappearance
* seaice_ecs.ncl: scatter plots of mean/trend of historical September Arctic sea ice extent
  vs 1st year of disappearance (RCP8.5) (similar to IPCC AR5 Chapter 12, Fig. 12.31a)
* seaice_trends.ncl: calculates ice extent trend distributions
  (similar to IPCC AR5 Chapter 9, Fig. 9.24c/d)
* seaice_tsline.ncl: creates a time series line plots of total sea ice area and extent (accumulated)
  for northern and southern hemispheres with optional multi-model mean and standard deviation. One
  value is used per model per year, either annual mean or the mean value of a selected month
  (similar to IPCC AR5 Chapter 9, Fig. 9.24a/b)
* seaice_yod.ncl: calculation of year of near disappearance of Arctic sea ice

User settings in recipe
-----------------------
#. Script seaice_ecs.ncl

   *Required settings (scripts)*

   * hist_exp: name of historical experiment (string)
   * month: selected month (1, 2, ..., 12) or annual mean ("A")
   * rcp_exp: name of RCP experiment (string)
   * region: region to be analyzed ( "Arctic" or "Antarctic")

   *Optional settings (scripts)*

   * fill_pole_hole: fill observational hole at North pole (default: False)
   * styleset: color style (e.g. "CMIP5")

   *Optional settings (variables)*

   * reference_dataset: reference dataset

#. Script seaice_trends.ncl

   *Required settings (scripts)*

   * month: selected month (1, 2, ..., 12) or annual mean ("A")
   * region: region to be analyzed ( "Arctic" or "Antarctic")

   *Optional settings (scripts)*

   * fill_pole_hole: fill observational hole at North pole, Default: False

   *Optional settings (variables)*

   * ref_model: array of references plotted as vertical lines

#. Script seaice_tsline.ncl

   *Required settings (scripts)*

   * region: Arctic, Antarctic
   * month: annual mean (A), or month number (3 = March, for Antarctic; 9 = September for Arctic)

   *Optional settings (scripts)*

   * styleset: for plot_type cycle only (cmip5, cmip6, default)
   * multi_model_mean: plot multi-model mean and standard deviation (default: False)
   * EMs_in_lg: create a legend label for individual ensemble members (default: False)
   * fill_pole_hole: fill polar hole (typically in satellite data) with sic = 1 (default: False)

#. Script seaice_yod.ncl

   *Required settings (scripts)*

   * month: selected month (1, 2, ..., 12) or annual mean ("A")
   * region: region to be analyzed ( "Arctic" or "Antarctic")

   *Optional settings (scripts)*

   * fill_pole_hole: fill observational hole at North pole, Default: False
   * wgt_file: netCDF containing pre-determined model weights

   *Optional settings (variables)*

   * ref_model: array of references plotted as vertical lines

Variables
---------

* sic (ocean-ice, monthly mean, longitude latitude time)
* areacello (fx, longitude latitude)

Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing; (2) see headers of cmorization scripts (in esmvaltool/utils/cmorizers/obs) for non-obs4MIPs data for download instructions.*

* HadISST (sic - esmvaltool/utils/cmorizers/obs/cmorize_obs_HadISST.ncl)

References
----------

* Massonnet, F. et al., The Cryosphere, 6, 1383-1394, doi: 10.5194/tc-6-1383-2012, 2012.
* Stroeve, J. et al., Geophys. Res. Lett., 34, L09501, doi:10.1029/2007GL029703, 2007.

Example plots
-------------

.. figure::  /recipes/figures/seaice/trend_sic_extend_Arctic_September_histogram.png
   :align:   center
   :width:   9cm

   Sea ice extent trend distribution for the Arctic in September
   (similar to IPCC AR5 Chapter 9, Fig. 9.24c). [seaice_trends.ncl]

.. figure::  /recipes/figures/seaice/extent_sic_Arctic_September_1960-2005.png
   :align:   center
   :width:   12cm

   Time series of total sea ice area and extent (accumulated) for the Arctic in September
   including multi-model mean and standard deviation (similar to IPCC AR5 Chapter 9, Fig. 9.24a).
   [seaice_tsline.ncl]

.. figure::  /recipes/figures/seaice/timeseries_rcp85.png
   :align:   center
   :width:   12cm

   Time series of September Arctic sea ice extent for individual CMIP5 models,
   multi-model mean and multi-model standard deviation, year of disappearance
   (similar to IPCC AR5 Chapter 12, Fig. 12.31e). [seaice_yod.ncl]

.. figure::  /recipes/figures/seaice/SSIE-MEAN_vs_YOD_sic_extend_Arctic_September_1960-2100.png
   :align:   center
   :width:   9cm

   Scatter plot of mean historical September Arctic sea ice extent vs 1st year of disappearance
   (RCP8.5) (similar to IPCC AR5 Chapter 12, Fig. 12.31a). [seaice_ecs.ncl]
.. _recipes_ensclus:

Ensemble Clustering - a cluster analysis tool for climate model simulations (EnsClus)
=====================================================================================


Overview
--------
EnsClus is a cluster analysis tool in Python, based on the k-means algorithm, for ensembles of climate model simulations.

Multi-model studies allow to investigate climate processes beyond the limitations of individual models by means of inter-comparison or averages of several members of an ensemble. With large ensembles, it is often an advantage to be able to group members according to similar characteristics and to select the most representative member for each cluster.

The user chooses which feature of the data is used to group the ensemble members by clustering: time mean, maximum, a certain percentile (e.g., 75% as in the examples below), standard deviation and trend over the time period. For each ensemble member this value is computed at each grid point, obtaining N lat-lon maps, where N is the number of ensemble members. The anomaly is computed subtracting the ensemble mean of these maps to each of the single maps. The anomaly is therefore computed with respect to the ensemble members (and not with respect to the time) and the Empirical Orthogonal Function (EOF) analysis is applied to these anomaly maps.

Regarding the EOF analysis, the user can choose either how many Principal Components (PCs) to retain or the percentage of explained variance to keep. After reducing dimensionality via EOF analysis, k-means analysis is applied using the desired subset of PCs.

The major final outputs are the classification in clusters, i.e. which member belongs to which cluster (in k-means analysis the number k of clusters needs to be defined prior to the analysis) and the most representative member for each cluster, which is the closest member to the cluster centroid.

Other outputs refer to the statistics of clustering: in the PC space, the minimum and the maximum distance between a member in a cluster and the cluster centroid (i.e. the closest and the furthest member), the intra-cluster standard deviation for each cluster (i.e. how much the cluster is compact).


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_ensclus.yml

Diagnostics are stored in diag_scripts/ensclus/

* ensclus.py

and subroutines

* ens_anom.py
* ens_eof_kmeans.py
* ens_plots.py
* eof_tool.py
* read_netcdf.py
* sel_season_area.py


User settings
-------------

*Required settings for script*

* season: season over which to perform seasonal averaging (DJF, DJFM, NDJFM, JJA)
* area: region of interest (EAT=Euro-Atlantic, PNA=Pacific North American, NH=Northern Hemisphere, EU=Europe)
* extreme: extreme to consider: XXth_percentile (XX can be set arbitrarily, e.g. 75th_percentile), mean (mean value over the period), maximum (maximum value over the period), std (standard deviation), trend (linear trend over the period)
* numclus: number of clusters to be computed
* perc: percentage of variance to be explained by PCs (select either this or numpcs, default=80)
* numpcs: number of PCs to retain (has priority over perc unless it is set to 0 (default))

*Optional settings for script*

* max_plot_panels: maximum number of panels (datasets) in a plot. When exceeded multiple plots are created. Default: 72


Variables
---------

* chosen by user (e.g., precipitation as in the example)


Observations and reformat scripts
---------------------------------

None.


References
----------

* Straus, D. M., S. Corti, and F. Molteni: Circulation regimes: Chaotic variability vs. SST forced predictability. J. Climate, 20, 2251–2272, 2007. https://doi.org/10.1175/JCLI4070.1


Example plots
-------------

.. figure:: /recipes/figures/ensclus/ensclus.png
   :width: 10cm

   Clustering based on the 75th percentile of historical summer (JJA) precipitation rate for CMIP5 models over 1900-2005. 3 clusters are computed, based on the principal components explaining 80% of the variance. The 32 models are grouped in three different clusters. The green cluster is the most populated with 16 ensemble members mostly characterized by a positive anomaly over central-north Europe. The red cluster counts 12 elements that exhibit a negative anomaly centered over southern Europe. The third cluster – labelled in blue- includes only 4 models showing a north-south dipolar precipitation anomaly, with a wetter than average Mediterranean counteracting dryer North-Europe. Ensemble members No.9, No.26 and No.19 are the “specimen” of each cluster, i.e. the model simulations that better represent the main features of that cluster. These ensemble members can eventually be used as representative of the whole possible outcomes of the multi-model ensemble distribution associated to the 32 CMIP5 historical integrations for the summer precipitation rate 75 th percentile over Europe when these outcomes are reduced from 32 to 3. The number of ensemble members of each cluster might provide a measure of the probability of occurrence of each cluster. 
.. _recipes_validation:

Zonal and Meridional Means
==========================

Overview
--------

This functional diagnostic takes two models designated by CONTROL and EXPERIMENT and compares them via a number of
analyses. Optionally a number of observational datasets can be added for processing. There are three types of standard analysis:
lat_lon, meridional_mean and zonal_mean. Each of these diagnostics can be run on a separate basis (each an entry to diagnostics/scripts).
The lat_lon analysis produces the following plots: a simple global plot for each variable for each dataset, a global plot for the
difference between CONTROL and EXPERIMENT, a global plot for the difference between CONTROL and each of the observational datasets.
The meridional_mean and zonal_mean produce variable vs coordinate (``latitude`` or ``longitude``) with both ``CONTROL`` and ``EXPERIMENT`` curves
in each plot, for the entire duration of time specified and also, if the user wishes, for each season (seasonal means): winter DJF, spring MAM, summer JJA, autumn SON (by setting ``seasonal_analysis: true`` in the recipe).

At least regridding on a common grid for all model and observational datasets should be performed in preprocessing (if datasets
are on different grids). Also note that it is allowed to use the same dataset (with varying parameters like experiment
or ensemble or mip) for both CONTROL and EXPERIMENT (as long as at least one data parameter is different).

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_validation.yml (CMIP5)
* recipe_validation_CMIP6.yml (CMIP6)

Diagnostics are stored in diag_scripts/

* validation.py
* shared/_validation.py

User settings
-------------

#. validation.py

   *Required settings for script*

   * title: title of the analysis, user defined;
   * control_model: control dataset name e.g. UKESM1-0-LL;
   * exper_model: experiment dataset name e.g. IPSL-CM6A-LR;
   * observational_datasets: list of at least one element; if no OBS wanted comment out; e.g. ['ERA-Interim'];
   * analysis_type: use any of: lat_lon, meridional_mean, zonal_mean;
   * seasonal_analysis: boolean, if seasonal means are needed e.g. ``true``;
   * save_cubes: boolean, save each of the plotted cubes in ``/work``; 

Variables
---------

* any variable

Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs or OBS or ana4mips can be used.*

* any observations
* it is important to note that all observational data should go through the same preprocessing as model data

References
----------

* none, basic technical analysis

Example plots
-------------

.. figure:: /recipes/figures/validation/Merid_Mean_DJF_longitude_tas_UKESM1-0-LL_vs_IPSL-CM6A-LR.png
   :width: 70 %

   Meridional seasonal mean for winter (DJF) comparison beween CMIP6 UKESM1 and IPSL models.

.. figure:: /recipes/figures/validation/Zonal_Mean_DJF_latitude_tas_UKESM1-0-LL_vs_IPSL-CM6A-LR.png
   :width: 70 %

   Zonal seasonal mean for winter (DJF) comparison beween CMIP6 UKESM1 and IPSL models.
.. _recipes_heatwaves_coldwaves:

Heat wave and cold wave duration
================================

Overview
--------

The goal of this diagnostic is to estimate the relative change in heat/cold wave characteristics  in future climates compared to a reference period using daily maximum or minimum temperatures.

The user can select whether to compute the frequency of exceedances or non-exceedances, which corresponds to extreme high or extreme low temperature events, respectively. The user can also select the minimum duration for an event to be classified as a heat/cold wave and the season of interest.

The diagnostic calculates the number of days in which the temperature exceeds or does not exceeds the necessary threshold for a consecutive number of days in future climate projections. The result is an annual time series of the total number of heat/cold wave days for the selected season at each grid point. The final output is the average number of heat/cold wave days for the selected season in the future climate projections.

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_heatwaves_coldwaves.yml

Diagnostics are stored in diag_scripts/magic_bsc/

* extreme_spells.R: calculates the heatwave or coldwave duration.


User settings
-------------

User setting files are stored in recipes/

#. recipe_heatwaves_coldwaves.yml

   *Required settings for script*

   * quantile: quantile defining the exceedance/non-exceedance threshold
   * min_duration: Min duration in days of a heatwave/coldwave event
   * Operator: either '>' for exceedances or '<' for non-exceedances
   * season: 'summer' or 'winter

Variables
---------

* tasmax or tasmin (atmos, daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Cardoso, S., Marta-Almeida, M., Carvalho, A.C., & Rocha, A. (2017). Heat wave and cold spell changes in Iberia for a future climate scenario. International Journal of Climatology, 37(15), 5192-5205. https://doi.org/10.1002/joc.5158

* Ouzeau, G., Soubeyroux, J.-M., Schneider, M., Vautard, R., & Planton, S. (2016). Heat waves analysis over France in present and future climate: Application of a new method on the EURO-CORDEX ensemble. Climate Services, 4, 1-12. https://doi.org/10.1016/J.CLISER.2016.09.002

* Wang, Y., Shi, L., Zanobetti, A., & Schwartz, J. D. (2016). Estimating and projecting the effect of cold waves on mortality in 209 US cities. Environment International, 94, 141-149. https://doi.org/10.1016/j.envint.2016.05.008

* Zhang, X., Hegerl, G., Zwiers, F. W., & Kenyon, J. (2005). Avoiding inhomogeneity in percentile-based indices of temperature extremes. Journal of Climate, 18(11), 1641-1651. https://doi.org/10.1175/JCLI3366.1


Example plots
-------------

.. _fig_heatwaves:
.. figure::  /recipes/figures/heatwaves/tasmax_extreme_spell_durationsummer_IPSL-CM5A-LR_rcp85_2020_2040.png
   :align:   center
   :width:   14cm

Mean number of summer days during the period 2060-2080 when the daily maximum near-surface air temperature exceeds the 80th quantile of the 1971-2000 reference period. The results are based on one RCP 8.5 scenario simulated by BCC-CSM1-1.
.. _recipes_multimodel_products:

Multi-model products
====================

Overview
--------

The goal of this diagnostic is to compute the multi-model ensemble mean for a set of models selected by the user for individual variables and different temporal resolutions (annual, seasonal, monthly).

After selecting the region (defined by the lowermost and uppermost longitudes and latitudes), the mean for the selected reference period is subtracted from the projections in order to obtain the anomalies for the desired period. In addition, the recipe computes the percentage of models agreeing on the sign of this anomaly, thus providing some indication on the robustness of the climate signal.

The output of the recipe consists of a colored map showing the time average of the multi-model mean anomaly and stippling to indicate locations where the percentage of models agreeing on the sign of the anomaly exceeds a threshold selected by the user. Furthermore, a time series of the area-weighted mean anomaly for the projections is plotted. For the plots, the user can select the length of the running window for temporal smoothing and choose to display the ensemble mean with a light shading to represent the spread of the ensemble or choose to display each individual models.



Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_multimodel_products.yml


Diagnostics are stored in diag_scripts/magic_bsc/

* multimodel_products.R - script for computing multimodel anomalies and their agreement.




User settings
-------------

User setting files are stored in recipes/

#. recipe_multimodel_products.yml

   *Required settings for script*

   * colorbar_lim: positive number specifying the range (-colorbar_lim ... +colorbar_lim) of the colorbar
     (0 = automatic colorbar scaling)
   * moninf: integer specifying the first month of the seasonal mean period to be computed
   * monsup: integer specifying the last month of the seasonal mean period to be computed, if it's null the anomaly of month indicated in moninf will be computed
   * agreement_threshold: integer between 0 and 100 indicating the threshold in percent for the minimum agreement between models on the sign of the multi-model mean anomaly for the stipling to be plotted
   * running_mean: integer indictating the length of the window for the running mean to be computed
   * time_series_plot: Either single or maxmin (plot the individual or the mean with shading between the max and min).


Variables
---------

* any Amon variable (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Hagedorn, R., Doblas-Reyes, F. J., Palmer, T. N., Nat E H Ag E D O R N, R. E., & Pa, T. N. (2005). The rationale behind the success of multi-model ensembles in seasonal forecasting-I. Basic concept, 57, 219–233. https://doi.org/10.3402/tellusa.v57i3.14657

* Weigel, A. P., Liniger, M. A., & Appenzeller, C. (2008). Can multi-model combination really enhance the prediction skill of probabilistic ensemble forecasts? Quarterly Journal of the Royal Meteorological Society, 134(630), 241–260. https://doi.org/10.1002/qj.210






Example plots
-------------

.. _fig_multimodprod:
.. figure::  /recipes/figures/multimodel_products/tas_JUN_multimodel-anomaly_2006_2099_1961_1990.png

Multi-model mean anomaly of 2-m air temperature during the future projection 2006-2099 in June considering the reference period 1961-1990 (colours). Crosses indicate that the 80% of models agree in the sign of the multi-model mean anomaly. The models selected are BCC-CSM1-1, MPI-ESM-MR and MIROC5 in the r1i1p1 ensembles for the RCP 2.6 scenario.
.. _recipes_schlund20esd:

Emergent constraints on equilibrium climate sensitivity in CMIP5: do they hold for CMIP6?
=========================================================================================

Overview
--------

This recipe reproduces the analysis of `Schlund et al., Earth Sys. Dyn.
(2020)`_. In this paper, emergent constraints on the equilibrium climate
sensitivity are evaluated on CMIP5 and CMIP6 models. Since none of the
considered emergent constraints have been developed on the CMIP6 ensemble, this
allows an out-of-sample testing of the emergent constraints. Most emergent
constraints show a reduced skill in CMIP6 when compared to CMIP5.

.. _`Schlund et al., Earth Sys. Dyn. (2020)`: https://doi.org/10.5194/esd-11-1233-2020


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_schlund20esd.yml

Diagnostics are stored in diag_scripts/

   * :ref:`climate_metrics/ecs.py<ecs.py>`
   * :ref:`climate_metrics/psi.py<psi.py>`
   * :ref:`emergent_constraints/ecs_scatter.ncl<ecs_scatter.ncl>`
   * :ref:`emergent_constraints/ecs_scatter.py<api.esmvaltool.diag_scripts.emergent_constraints.ecs_scatter>`
   * :ref:`emergent_constraints/multiple_constraints.py<api.esmvaltool.diag_scripts.emergent_constraints.multiple_constraints>`

More details on the emergent constraint module are given in the API
documentation which is available
:ref:`here<api.esmvaltool.diag_scripts.emergent_constraints>`.


Variables
---------

* *cl* (atmos, monthly, longitude, latitude, level, time)
* *clt* (atmos, monthly, longitude, latitude, time)
* *hur* (atmos, monthly, longitude, latitude, level, time)
* *hus* (atmos, monthly, longitude, latitude, level, time)
* *pr* (atmos, monthly, longitude, latitude, time)
* *rsdt* (atmos, monthly, longitude, latitude, time)
* *rsut* (atmos, monthly, longitude, latitude, time)
* *rsutcs* (atmos, monthly, longitude, latitude, time)
* *rtnt* or *rtmt* (atmos, monthly, longitude, latitude, time)
* *ta* (atmos, monthly, longitude, latitude, level, time)
* *tas* (atmos, monthly, longitude, latitude, time)
* *tasa* (atmos, monthly, longitude, latitude, time)
* *tos* (atmos, monthly, longitude, latitude, time)
* *ts* (atmos, monthly, longitude, latitude, time)
* *va* (atmos, monthly, longitude, latitude, level, time)
* *wap* (atmos, monthly, longitude, latitude, level, time)


Observations and reformat scripts
---------------------------------

* AIRS_ (*hur*, *hus*)
* CERES-EBAF_ (*rsut*, *rsutcs*, *rsdt*)
* ERA-Interim_ (*hur*, *ta*, *va*, *wap*)
* GPCP-SG_ (*pr*)
* HadCRUT4_ (*tasa*)
* HadISST_ (*ts*)
* MLS-AURA_ (*hur*)

.. _AIRS: https://opendata.dwd.de/climate_environment/GPCC/html/fulldata-monthly_v2018_doi_download.html
.. _CERES-EBAF: https://opendata.dwd.de/climate_environment/GPCC/html/fulldata-monthly_v2018_doi_download.html
.. _ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda/
.. _GPCP-SG: https://opendata.dwd.de/climate_environment/GPCC/html/fulldata-monthly_v2018_doi_download.html
.. _HadCRUT4: https://crudata.uea.ac.uk/cru/data/temperature/
.. _HadISST: http://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
.. _MLS-AURA: https://disc.gsfc.nasa.gov/datasets/ML2RHI_004/summary


References
----------

* Schlund, M., Lauer, A., Gentine, P., Sherwood, S. C., and Eyring, V.:
  Emergent constraints on equilibrium climate sensitivity in CMIP5: do they
  hold for CMIP6?, Earth Syst. Dynam., 11, 1233–1258,
  `<https://doi.org/10.5194/esd-11-1233-2020>`_, 2020.


Example plots
-------------

.. _fig_schlund20esd_1:
.. figure:: /recipes/figures/schlund20esd/SHL_scatter.png
   :align: center
   :width: 50%

   Emergent relationship (solid blue and orange lines) of the `Sherwood et al.
   (2014) <https://doi.org/10.1038/nature12829>`_ emergent constraint, which is
   based on the lower tropospheric mixing index (LTMI). The numbers correspond
   to individual CMIP models. The shaded area around the regression line
   corresponds to the standard prediction error, which defines the error in the
   regression model itself. The vertical dashed black line corresponds to the
   observational reference with its uncertainty range given as standard error
   (gray shaded area). The horizontal dashed lines show the best estimates of
   the constrained ECS for CMIP5 (blue) and CMIP6 (orange). The colored dots
   mark the CMIP5 (blue) and CMIP6 (orange) multi-model means.

.. _fig_schlund20esd_2:
.. figure:: /recipes/figures/schlund20esd/SHL_pdf.png
   :align: center
   :width: 50%

   Probability densities for the constrained ECS (solid lines) and the
   unconstrained model ensembles (histograms) of the emergent relationship
   shown in the figure above.

.. _fig_schlund20esd_3:
.. figure:: /recipes/figures/schlund20esd/ZHA_scatter.png
   :align: center
   :width: 50%

   Emergent relationship of the `Zhai et al. (2015)
   <https://doi.org/10.1002/2015GL065911>`_ emergent constraint for different
   subsets of CMIP5 models. Blue circles show the 15 CMIP5 models used in the
   original publication (except for CESM1-CAM5); the solid blue line and blue
   shaded area show the emergent relationships evaluated on these models
   including the uncertainty range. In this study, 11 more CMIP5 models have
   been added (red circles). The corresponding emergent relationship that
   considers all available CMIP5 models is shown in red colors. This
   relationship shows a considerably lower coefficient of determination
   (:math:`R^2`) and higher *p*-value than the relationship using the original
   subset of CMIP5 models. The vertical dashed line and shaded area correspond
   to the observational reference, and the horizontal dashed lines show the
   corresponding ECS constraints using this observation.
.. _recipes_extreme_events:

Extreme Events Indices (ETCCDI)
===============================


Overview
--------

This diagnostic uses the standard climdex.pcic.ncdf R library to
compute the 27 climate change indices specified by
the joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices http://etccdi.pacificclimate.org/.
The needed input fields are daily average precipitation flux and minimum, maximum and average daily surface temperatures.
The recipe reproduces panels of figure 9.37 of the IPCC AR5 report, producing both a Gleckler plot,
with relative error metrics for the CMIP5 temperature and precipitation extreme indices,
and timeseries plots comparing the ensemble spread with observations.
For plotting 1 to 4 observational reference datasets are supported. If no observational reference datasets are given, the plotting routines do not work, however, index generation without plotting is still possible.
All datasets are regridded to a common grid and considered only over land.

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_extreme_events.yml

Diagnostics are stored in diag_scripts/extreme_events/

* ExtremeEvents.r

and subroutines

* common_climdex_preprocessing_for_plots.r
* make_Glecker_plot2.r
* make_timeseries_plot.r
* cfg_climdex.r
* cfg_extreme.r

User settings
-------------

*Required settings for script*

* reference_datasets: list containing the reference datasets to compare with
* timeseries_idx: list of indices to compute for timeseries plot.
  The syntax is "XXXETCCDI_TT", where "TT" can be either "yr" or "mon"
  (yearly or monthly indices are computed) and "XXX" can be one of the following:
  "altcdd", "altcsdi", "altcwd", "altwsdi", "cdd", "csdi", "cwd",
  "dtr", "fd", "gsl", "id", "prcptot", "r10mm", "r1mm", "r20mm",
  "r95p", "r99p", "rx1day", "rx5day", "sdii", "su", "tn10p",
  "tn90p", "tnn", "tnx", "tr", "tx10p", "tx90p", "txn", "txx", "wsdi".
  The option "mon" for "TT" can be only used in combination with one of:
  "txx", "tnx", "txn", "tnn", tn10p", "tx10p", "tn90p", "tx90p", "dtr", "rx1day", "rx5day".
* gleckler_idx: list of indices to compute for Gleckler plot. Same syntax as above.
  The diagnostic computes all unique indices specified in either ``gleckler_idx`` or ``timeseries_idx``.
  If at least one "mon" index is selected, the indices are computed but no plots are produced.
* base_range: a list of two years to specify the range to be used as "base range" for climdex
  (the period in which for example reference percentiles are computed)

*Optional settings for script*

* regrid_dataset: name of dataset to be used as common target for regridding. If missing the first reference dataset is used
* mip_name: string containing the name of the model ensemble, used for titles and labels in the plots (default: "CMIP")
* analysis_range: a list of two years to specify the range to be used for the analysis in the plots.
  The input data will need to cover both ``analysis_range`` and ``base_range``. If missing the full period covered by the
  input datasets will be used.
* ts_plt: (logical) if to produce the timeseries or not (default: true)
* glc_plt: (logical) if to produce the Gleckler or not (default: true)
* climdex_parallel: number of parallel threads to be used for climdex calculation (default: 4). Also the logical ``false`` can be passed to switch off parallel computation.
* normalize: (logical) if to detrend and normalize with the standard deviation for the datasets for use in the timeseries plot. When this option is used the data for the following indices  are detrended and normalized in the timeseries plots: "altcdd", "altcsdi", "altcwd", "altwsdi", "cdd",  "cwd","dtr", "fd", "gsl", "id", "prcptot", "r10mm", "r1mm", "r20mm", "r95p", "r99p", "rx1day", "rx5day", "sdii", "su", "tnn", "tnx", "tr", "txn","txn","txx" (default: false)

Additional optional setting controlling the plots:

* Timeseries plots:

  * ts_png_width: width for png figures (dafult: 640)
  * ts_png_height: height for png figures (default: 480)
  * ts_png_units: units for figure size (default: "px")
  * ts_png_pointsize: fontsize (default: 12)
  * ts_png_bg: background color (default: "white")
  * ts_col_list: list of colors for lines (default: ["dodgerblue2", "darkgreen", "firebrick2", "darkorchid", "aquamarine3"]``)
  * ts_lty_list: list of linetypes (default: [1, 4, 2, 3, 5])
  * ts_lwd_list: list of linewidths (default: [2, 2, 2, 2, 2])

* Gleckler plot:

  * gl_png_res: height for png figures (default: 480).
    The width of the figure is computed automatically.
  * gl_png_units: units for figure size (default: "px")
  * gl_png_pointsize: fontsize (default: 12)
  * gl_png_bg: background color (default: "white")
  * gl_mar_par: page margins vector (default: [10, 4, 3, 14])
  * gl_rmsespacer: spacing of RMSE column (default: 0.01)
  * gl_scaling_factor: scaling factor for colorscale height (default: 0.9)
  * gl_text_scaling_factor: scaling factor for text size (default: 1.0)
  * gl_xscale_spacer_rmse: horizontal posizion of coloured colorbar (default: 0.05)
  * gl_xscale_spacer_rmsestd: horizontal posizion of gray colorbar (default: 0.05)
  * gl_symb_scaling_factor: scaling factor for white "symbol" square explaining the partition (default: 1.0)
  * gl_symb_xshift: horizontal position of the symbol box (default: 0.2)
  * gl_symb_yshift: vertical position of the symbol box (default: 0.275)
  * gl_text_symb_scaling_factor: scaling factor for text to be used for symbol box (default: 0.5)

Variables
---------

* tas (atmos, daily mean, longitude latitude time)
* tasmin (atmos, daily minimum, longitude latitude time)
* tasmax (atmos, daily maximum, longitude latitude time)
* pr (atmos, daily mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

None.


References
----------

* Zhang, X., Alexander, L., Hegerl, G. C., Jones, P., Klein Tank, A., Peterson, T. C., Trewin, B., Zwiers, F. W., Indices for monitoring changes in extremes based on daily temperature and precipitation data, WIREs Clim. Change, doi:10.1002/wcc.147, 2011

* Sillmann, J., V. V. Kharin, X. Zhang, and F. W. Zwiers, Climate extreme indices in the CMIP5 multi-model ensemble. Part 1: Model evaluation in the present climate. J. Geophys. Res., doi:10.1029/2012JD018390, 2013


Example plots
-------------

.. figure:: /recipes/figures/extreme_events/gleckler.png
   :width: 12cm

   Portrait plot of relative error metrics for the CMIP5 temperature and precipitation extreme indices evaluated over 1981-2000. Reproduces Fig. 9.37 of the IPCC AR5 report, Chapter 9.

.. figure:: /recipes/figures/extreme_events/cdd_timeseries.png
   :width: 10cm

   Timeseries of the Consecutive Dry Days index over 1981-2000 for a selection of CMIP5 models, the CMIP5 multi-model mean (CMIP) and ERA-Interim. Shading is used to reproduce the multi-model spread.
.. _recipes_cmug_h2o:

Evaluate water vapor short wave radiance absorption schemes of ESMs with the observations, including ESACCI data.
==========================================================================================================================

Overview
--------

The recipe contains several diagnostics to use ESACCI water vapour data to evaluate CMIP models.

The diagnostic deangelisf3f4.py reproduces figures 3 and 4 from `DeAngelis et al. (2015)`_:
See also doc/sphinx/source/recipes/recipe_deangelis15nat.rst
This paper compares models with different schemes for water vapor short wave radiance absorption with the observations.
Schemes using pseudo-k-distributions with more than 20 exponential terms show the best results.

The diagnostic diag_tropopause.py plots given variable at cold point tropopause height,
here Specific Humidity (hus) is used. This will be calculated from the ESACCI water vapour data CDR-4, which are planed to consist of
three-dimensional vertically resolved monthly mean water vapour data (in ppmv) with
spatial resolution of 100 km, covering the troposphere and lower stratosphere.
The envisaged coverage is 2010-2014. The calculation of hus from water vapour in ppmv will be part of the cmorizer.
Here, ERA-Interim data are used.

The diagnostic diag_tropopause_zonalmean.py plots zonal mean for given variable for
all pressure levels between 250 and 1hPa and at cold point tropopause height.
Here Specific Humidity (hus) is used. This will be calculated from the
ESACCI water vapour data CDR-3, which are planed to contain
the vertically resolved water vapour ECV in units of ppmv (volume mixing ratio) and will be provided as
zonal monthly means on the SPARC Data Initiative latitude/pressure level grid
(SPARC, 2017; Hegglin et al., 2013). It covers the vertical range between 250 hPa and 1 hPa,
and the time period 1985 to the end of 2019. The calculation of hus from water vapour in ppmv will be
part of the cmorizer. Here, ERA-Interim  data are used.


.. _`DeAngelis et al. (2015)`: https://www.nature.com/articles/nature15770


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_cmug_h2o.yml

Diagnostics are stored in diag_scripts/

   * deangelis15nat/deangelisf3f4.py

   * cmug_h2o/diag_tropopause.py

   * cmug_h2o/diag_tropopause_zonalmean.py


User settings in recipe
-----------------------

The recipe can be run with different CMIP5 and CMIP6 models.

deangelisf3f4.py:
For each model, two experiments must be given:
a pre industrial control run, and a scenario with 4 times CO\ :sub:`2`\.
Possibly, 150 years should be given, but shorter time series work as well.
Currently, HOAPS data are included as place holder for expected ESACCI-WV data, type CDR-2:
Gridded monthly time series of TCWV in units of kg/m2 (corresponds to prw)
that cover the global land and ocean areas with a spatial resolution of 0.05° / 0.5°
for the period July 2002 to December 2017.


Variables
---------

deangelisf3f4.py:
* *rsnstcs* (atmos, monthly, longitude, latitude, time)
* *rsnstcsnorm* (atmos, monthly, longitude, latitude, time)
* *prw* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)


diag_tropopause.py:
* *hus* (atmos, monthly, longitude, latitude, time, plev)
* *ta* (atmos, monthly, longitude, latitude, time, plev)


diag_tropopause_zonalmean.py:
* *hus* (atmos, monthly, longitude, latitude, time, plev)
* *ta* (atmos, monthly, longitude, latitude, time, plev)


Observations and reformat scripts
---------------------------------

deangelisf3f4.py:

* *rsnstcs*:
   CERES-EBAF

* *prw*
   HOAPS, planed for ESACCI-WV data, type CDR-2

diag_tropopause.py:

* *hus*
   ERA-Interim, ESACCI water vapour paned

diag_tropopause_zonalmean.py:

* *hus*
   ERA-Interim, ESACCI water vapour paned


References
----------

* DeAngelis, A. M., Qu, X., Zelinka, M. D., and Hall, A.: An observational radiative constraint on hydrologic cycle intensification, Nature, 528, 249, 2015.


Example plots
-------------



.. _fig_deangelis_cmug_cdr2:
.. figure:: /recipes/figures/deangelis15nat/fig_deangelis_cmug_cdr2.png
   :align: center
   :width: 50%

   Scatter plot and regression line computed between the ratio of the change of net short wave radiation (rsnst) and the change of the Water Vapor Path (prw) against the ratio of the change of netshort wave radiation for clear skye (rsnstcs) and the the change of surface temperature (tas). The width of horizontal shading for models and the vertical dashed lines for observations (Obs.) represent statistical uncertainties of the ratio, as the 95% confidence interval (CI) of the regression slope to the rsnst versus prw curve. For the prw observations ESACCI CDR-2 data from 2003 to 2014 are used.

.. _fig_ERA-Interim_Cold_point_tropopause_Specific_Humidity_map:
.. figure:: /recipes/figures/cmug_h2o/fig_ERA-Interim_Cold_point_tropopause_Specific_Humidity_map.png
   :align: center
   :width: 50%

   Map of the average Specific Humidity (hus) at the cold point tropopause from ERA-Interim data. The diagnostic averages the complete time series, here 2010-2014.

.. _fig_ERA-Interim_Cold_point_tropopause_Specific_Humidity:
.. figure:: /recipes/figures/cmug_h2o/fig_ERA-Interim_Cold_point_tropopause_Specific_Humidity.png
   :align: center
   :width: 50%

   Latitude versus time plot of the Specific Humidity (hus) at the cold point tropopause from ERA-Interim data.

.. _fig_ERA-Interim_Zonal_mean_Specific_Humidity:
.. figure:: /recipes/figures/cmug_h2o/fig_ERA-Interim_Zonal_mean_Specific_Humidity.png
   :align: center
   :width: 50%

   Zonal average Specific Humidity (hus) between 250 and 1 hPa from ERA-Interim data. The diagnostic averages the complete time series, here 1985-2014.

.. _fig_profile_Specific_Humidity:
.. figure:: /recipes/figures/cmug_h2o/fig_profile_Specific_Humidity.png
   :align: center
   :width: 50%

   Average Specific Humidity (hus) profile between 250 and 1 hPa from ERA-Interim and CMIP6 model data. The diagnostic averages the complete time series, here 1985-2014.


.. _recipes_landcover:

Landcover diagnostics
=====================


Overview
--------

The diagnostic computes the accumulated and fractional extent of major land cover classes,
namely bare soil, crops, grasses, shrubs and trees. The numbers are compiled for the whole
land surface as well as separated into Tropics, northern Extratropics and southern Extratropics.
The cover fractions are compared to ESA-CCI land cover data.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_landcover.yml

Diagnostics are stored in diag_scripts/landcover/

    * landcover.py: bar plots showing the accumulated area and mean fractional coverage for five land
      cover classes for all experiments as well as their bias compared to observations.


User settings
-------------

script landcover.py

   *Required settings for script*

   * reference_dataset: land cover extent dataset for comparison. The script was developed using
     ESACCI-LANDCOVER observations.

   *Optional settings for script*

   * comparison: [variable, model] Choose whether one plot per land cover class is generated comparing
     the different experiments (default) or one plot per model comparing the different
     land cover classes.
   * colorscheme: Plotstyle used for the bar plots. A list of available style is found at
     https://matplotlib.org/gallery/style_sheets/style_sheets_reference.html. Seaborn is used as default.


Variables
---------

* baresoilFrac (land, monthly mean, time latitude longitude)
* grassFrac    (land, monthly mean, time latitude longitude)
* treeFrac     (land, monthly mean, time latitude longitude)
* shrubFrac    (land, monthly mean, time latitude longitude)
* cropFrac     (land, monthly mean, time latitude longitude)


Observations and reformat scripts
---------------------------------

ESA-CCI land cover data (Defourny et al., 2015) needs to be downloaded manually by the user and converted to netCDF files
containing the grid cell fractions for the five major land cover types. The data and a conversion tool
are available at https://maps.elie.ucl.ac.be/CCI/viewer/ upon registration. After obtaining the data and the user
tool, the remapping to 0.5 degree can be done with::

  ./bin/aggregate-map.sh
  -PgridName=GEOGRAPHIC_LAT_LON
  -PnumRows=360
  -PoutputLCCSClasses=true
  -PnumMajorityClasses=0
  ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7b.nc

Next, the data needs to be aggregated into the five major classes (PFT) similar to the study of Georgievski & Hagemann (2018)
and converted from grid cell fraction into percentage.

+--------------+-------------------------------------------------------------------------------------------------------------+
| PFT          | ESA-CCI Landcover Classes                                                                                   |
+==============+=============================================================================================================+
| baresoilFrac | Bare_Soil                                                                                                   |
+--------------+-------------------------------------------------------------------------------------------------------------+
| cropFrac     | Managed_Grass                                                                                               |
+--------------+-------------------------------------------------------------------------------------------------------------+
| grassFrac    | Natural_Grass                                                                                               |
+--------------+-------------------------------------------------------------------------------------------------------------+
| shrubFrac    | Shrub_Broadleaf_Deciduous + Shrub_Broadleaf_Evergreen + Shrub_Needleleaf_Evergreen                          |
+--------------+-------------------------------------------------------------------------------------------------------------+
| treeFrac     | Tree_Broadleaf_Deciduous + Tree_Broadleaf_Evergreen + Tree_Needleleaf_Deciduous + Tree_Needleleaf_Evergreen |
+--------------+-------------------------------------------------------------------------------------------------------------+

Finally, it might be necessary to adapt the grid structure to the experiments files, e.g converting the -180 --> 180 degree grid
to 0 --> 360 degree and inverting the order of latitudes. Note, that all experiments will be regridded onto the grid of the land
cover observations, thus it is recommended to convert to the coarses resolution which is sufficient for the planned study.
For the script development, ESA-CCI data on 0.5 degree resolution was used with land cover data averaged over the
2008-2012 period.


References
----------

* Defourny et al. (2015): ESA Land Cover Climate Change Initiative (ESA LC_cci) data:
  ESACCI-LC-L4-LCCS-Map-300m-P5Y-[2000,2005,2010]-v1.6.1 via Centre for Environmental Data Analysis
* Georgievski, G. & Hagemann, S. Characterizing uncertainties in the ESA-CCI land cover map of the epoch 2010 and their impacts on MPI-ESM climate simulations,
  Theor Appl Climatol (2018). https://doi.org/10.1007/s00704-018-2675-2


Example plots
-------------

.. _fig_landcover_1:
.. figure::  /recipes/figures/landcover/area_treeFrac.png
   :align:   center
   :width:   14cm

   Accumulated tree covered area for different regions and experiments.

.. _fig_landcover_2:
.. figure::  /recipes/figures/landcover/frac_grassFrac.png
   :align:   center
   :width:   14cm

   Average grass cover fraction for different regions and experiments

.. _fig_landcover_3:
.. figure::  /recipes/figures/landcover/bias_CMIP5_MPI-ESM-LR_rcp85_r1i1p1.png
   :align:   center
   :width:   14cm

   Biases in five major land cover fractions for different regions and one experiment.
.. _recipes_spei:

Standardized Precipitation-Evapotranspiration Index (SPEI)
==========================================================

Overview
--------
Droughts can be separated into three main types: meteorological, hydrological, and agricultural drought.

Common for all types is that a drought needs to be put in context of local and seasonal characteristics, i.e. a drought should not be defined with an absolute threshold, but as an anomalous condition.

Meteorological droughts are often described using the standardized precipitation index (SPI; McKee et al, 1993), which in a standardized way describes local precipitation anomalies. It is calculated on monthly mean precipitation, and is therefore not accounting for the intensity of precipitation and the runoff process. Because SPI does not account for evaporation from the ground, it lacks one component of the water fluxes at the surface and is therefore not compatible with the concept of hydrological drought.

A hydrological drought occurs when low water supply becomes evident, especially in streams, reservoirs, and groundwater levels, usually after extended periods of meteorological drought. GCMs normally do not simulate hydrological processes in sufficient detail to give deeper insights into hydrological drought processes. Neither do they properly describe agricultural droughts, when crops become affected by the hydrological drought. However, hydrological drought can be estimated by accounting for evapotranspiration, and thereby estimate the surface retention of water. The standardized precipitation-evapotranspiration index (SPEI; Vicente-Serrano et al., 2010) has been developed to also account for temperature effects on the surface water fluxes. Evapotranspiration is not normally calculated in GCMs, so SPEI often takes other inputs to estimate the evapotranspiration. Here, the Thornthwaite (Thornthwaite, 1948) method based on temperature is applied.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_spei.yml


Diagnostics are stored in diag_scripts/droughtindex/

    * diag_spi.R: calculate the SPI index

    * diag_spei.R: calculate the SPEI index


User settings
-------------

#. Script diag_spi.py

   *Required settings (script)*

   * reference_dataset: dataset_name
     The reference data set acts as a baseline for calculating model bias.

#. Script diag_spei.py

   *Required settings (script)*

   * reference_dataset: dataset_name
     The reference data set acts as a baseline for calculating model bias.


Variables
---------

* pr      (atmos, monthly mean, time latitude longitude)
* tas     (atmos, monthly mean, time latitude longitude)


References
----------
* McKee, T. B., Doesken, N. J., & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. In Proceedings of the 8th Conference on Applied Climatology (Vol. 17, No. 22, pp. 179-183). Boston, MA: American Meteorological Society.

* Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010). A multiscalar drought index sensitive to global warming: the standardized precipitation evapotranspiration index. Journal of climate, 23(7), 1696-1718.


Example plots
-------------

.. _fig_spei:
.. figure::  /recipes/figures/spei/histogram_spei.png
   :align:   center
   :width:   14cm

   (top) Probability distribution of the standardized precipitation-evapotranspiration index of a sub-set of the CMIP5 models, and (bottom) bias relative to the CRU reference data set.

.. _fig_spi:
.. figure::  /recipes/figures/spei/histogram_spi.png
   :align:   center
   :width:   14cm

   (top) Probability distribution of the standardized precipitation index of a sub-set of the CMIP5 models, and (bottom) bias relative to the CRU reference data set.
.. _nml_oceanmetrics:

Ocean metrics
=============

Overview
--------

The Southern Ocean is central to the global climate and the global carbon cycle, and to the climate’s response to increasing levels of atmospheric greenhouse gases. Global coupled climate models and earth system models, however, vary widely in their simulations of the Southern Ocean and its role in, and response to, the ongoing anthropogenic trend. Observationally-based metrics are critical for discerning processes and mechanisms, and for validating and comparing climate and earth system models. New observations and understanding have allowed for progress in the creation of observationally-based data/model metrics for the Southern Ocean.

The metrics presented in this recipe provide a means to assess multiple simulations relative to the best available observations and observational products. Climate models that perform better according to these metrics also better simulate the uptake of heat and carbon by the Southern Ocean. Russell et al. 2018 assessed only a few of the available CMIP5 simulations, but most of the available CMIP5 and CMIP6 climate models can be analyzed with these recipes.

The goal is to create a recipe for recreation of metrics in Russell, J.L., et al., 2018, J. Geophys. Res. – Oceans, 123, 3120-3143, doi: 10.1002/2017JC013461.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_russell18jgr.yml

Diagnostics are stored in diag_scripts/russell18jgr/

* russell18jgr-polar.ncl (figures 1, 7, 8): calculates and plots annual-mean variables (tauu, sic, fgco2, pH) as polar contour map.
* russell18jgr-fig2.ncl:  calculates and plots The zonal and annual means of the zonal wind stress (N/m\ :sup:`2`\).
* russell18jgr-fig3b.ncl: calculates and plots the latitudinal position of Subantarctic Front. Using definitions from Orsi et al (1995).
* russell18jgr-fig3b-2.ncl: calculates and plots the latitudinal position of Polar Front. Using definitions from Orsi et al (1995).
* russell18jgr-fig4.ncl:  calculates and plots the zonal velocity through Drake Passage (at 69W) and total transport through the passage if the volcello file is available.
* russell18jgr-fig5.ncl:  calculates and plots the mean extent of sea ice for September(max) in blue and mean extent of sea ice for February(min) in red. 
* russell18jgr-fig5g.ncl: calculates and plots the annual cycle of sea ice area in southern ocean.
* russell18jgr-fig6a.ncl: calculates and plots the density layer based volume transport(in Sv) across 30S based on the layer definitions in Talley (2008).
* russell18jgr-fig6b.ncl: calculates and plots the Density layer based heat transport(in PW) across 30S based on the layer definitions in Talley (2008).
* russell18jgr-fig7h.ncl: calculates and plots the zonal mean flux of fgco2 in gC/(yr * m\ :sup:`2`\). 
* russell18jgr-fig7i.ncl: calculates and plots the cumulative integral of the net CO2 flux from 90S to 30S (in PgC/yr).
* russell18jgr-fig9a.ncl: calculates and plots the scatter plot of the width of the Southern Hemisphere westerly wind band against the annual-mean integrated heat uptake south of 30S (in PW), along with the line of best fit.
* russell18jgr-fig9b.ncl: calculates and plots the scatter plot of the width of the Southern Hemisphere westerly wind band against the annual-mean integrated carbon uptake south of 30S (in Pg C/yr), along with the line of best fit.
* russell18jgr-fig9c.ncl: calculates and plots the scatter plot of the net heat uptake south of 30S (in PW) against the annual-mean integrated carbon uptake south of 30S (in Pg C/yr), along with the line of best fit.

User settings in recipe
-----------------------

#. Script russell18jgr-polar.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.
   * max_lat   : -30.0

   *Optional settings (scripts)*

   * grid_max  :  0.4 (figure 1),  30 (figure 7), 8.2 (figure 8)
   * grid_min  : -0.4 (figure 1), -30 (figure 7), 8.0 (figure 8)
   * grid_step :  0.1 (figure 1), 2.5 (figure 7), 0.1 (figure 8)
   * colormap  : BlWhRe (figure 7)
   * colors    : [[237.6, 237.6, 0.], [ 255, 255, 66.4], [255, 255, 119.6], [255, 255, 191.8], [223.8, 191.8, 223.8], [192.8, 127.5, 190.8], [161.6, 65.3, 158.6], [129.5, 1.0, 126.5] ] (figure 1)
     [[132,12,127], [147,5,153], [172,12,173], [195,33,196], [203,63,209], [215,89,225], [229,117,230], [243,129,238], [253,155,247], [255,178,254], [255,255,255],
     [255,255,255], [126,240,138], [134,234,138], [95,219,89], [57,201,54], [39,182,57], [33,161,36], [16,139,22], [0,123,10], [6,96,6], [12,77,9.0] ]      (figure 8)
   * max_vert  :  1 - 4 (user preference)
   * max_hori  :  1 - 4 (user preference)
   * grid_color:  blue4 (figure 8)
   * labelBar_end_type:  ExcludeOuterBoxes (figure 1), both_triangle (figure 7, 8)
   * unitCorrectionalFactor: -3.154e+10 (figure 7)
   * new_units : "gC/ (m~S~2~N~ * yr)" (figure 7)

   *Required settings (variables)*

   * additional_dataset: datasets to plot.

   *Optional settings (variables)*

   * none


#. Script russell18jgr-fig2.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig3b.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig3b-2.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig4.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * max_vert  :  1 - 4 (user preference)
   * max_hori  :  1 - 4 (user preference)
   * unitCorrectionalFactor: 100 (m/s to cm/s)
   * new_units : "cm/s"


#. Script russell18jgr-fig5.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.
   * max_lat  : -45.0

   *Optional settings (scripts)*

   * max_vert  :  1 - 4 (user preference)
   * max_hori  :  1 - 4 (user preference)


#. Script russell18jgr-fig5g.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig6a.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig6b.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig7h.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig7i.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none

#. Script russell18jgr-fig9a.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig9b.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none


#. Script russell18jgr-fig9c.ncl

   *Required settings (scripts)*

   * styleset : CMIP5(recommended), default, etc.
   * ncdf     : default(recommended), CMIP5, etc.

   *Optional settings (scripts)*

   * none



Variables
---------

* tauu (atmos, monthly mean, longitude latitude time)
* tauuo, hfds, fgco2 (ocean, monthly mean, longitude latitude time)
* thetao, so, vo (ocean, monthly mean, longitude latitude lev time)
* pH (ocnBgchem, monthly mean, longitude latitude time)
* uo (ocean, monthly mean, longitude latitude lev time)
* sic (seaIce, monthly mean, longitude latitude time)

Observations and reformat scripts
---------------------------------

Note: WOA data has not been tested with reciepe_russell18jgr.yml and
      corresponding diagnostic scripts.

* WOA (thetao, so - esmvaltool/utils/cmorizers/obs/cmorize_obs_woa.py)

References
----------

* Russell, J.L., et al., 2018, J. Geophys. Res. – Oceans, 123, 3120-3143. https://doi.org/10.1002/2017JC013461

* Talley, L.D., 2003. Shallow,intermediate and deep overturning components of the global heat budget. Journal of Physical Oceanography 33, 530–560


Example plots
-------------

.. _fig_russell_1:
.. figure::  /recipes/figures/russell18jgr/Fig1_polar-contour_tauu_1986-2005.png
   :align:   center
   :width: 50%

   Figure 1: Annual-mean zonal wind stress (tauu - N/m\ :sup:`2`\) with eastward wind stress as positive plotted as a polar contour map. 

.. _fig_russell_2:
.. figure::  /recipes/figures/russell18jgr/Fig2_1986-2005.png
   :align:   center
   :width: 50%

   Figure 2: The zonal and annual means of the zonal wind stress (N/m\ :sup:`2`\) plotted in a line plot.

.. _fig_russell_3a:
.. figure::  /recipes/figures/russell18jgr/Fig3_Polar-Front.png
   :align:   center
   :width: 50%

   Figure 3a: The latitudinal position of Subantarctic Front using definitions from Orsi et al (1995).

.. _fig_russell_3b:
.. figure::  /recipes/figures/russell18jgr/Fig3_Subantarctic-Fronts.png
   :align:   center
   :width: 50%

   Figure 3b: The latitudinal position of Polar Front using definitions from Orsi et al (1995).

.. _fig_russell_4:
.. figure::  /recipes/figures/russell18jgr/Fig4_Drake_passage.png
   :align:   center
   :width: 50%

   Figure 4: Time averaged zonal velocity through Drake Passage (at 69W, in cm/s, eastward is positive). The total transport by the ACC is calculated if volcello file is available.

.. _fig_russell_5:
.. figure::  /recipes/figures/russell18jgr/Fig5_sic-max-min.png
   :align:   center
   :width: 50%

   Figure 5: Mean extent of sea ice for September(max) in blue and February(min) in red plotted as polar contour map.


.. _fig_russell_5g:
.. figure::  /recipes/figures/russell18jgr/Fig5g_sic-line.png
   :align:   center
   :width: 50%

   Figure 5g: Annual cycle of sea ice area in southern ocean as a line plot (monthly climatology).

.. _fig_russell_6a:
.. figure::  /recipes/figures/russell18jgr/Fig6a.png
   :align:   center
   :width: 50%

   Figure 6a: Density layer based volume transport (in Sv) across 30S based on the layer definitions in Talley (2008).

.. _fig_russell_6b:
.. figure::  /recipes/figures/russell18jgr/Fig6b.png
   :align:   center
   :width: 50%

   Figure 6b: Density layer based heat transport(in PW) across 30S based on the layer definitions in Talley (2008).


.. _fig_russell_7:
.. figure::  /recipes/figures/russell18jgr/Fig7_fgco2_polar.png
   :align:   center
   :width: 50%

   Figure 7: Annual mean CO\ :sub:`2`\  flux (sea to air, gC/(yr * m\ :sup:`2`\), positive (red) is out of the ocean) as a polar contour map.

.. _fig_russell_7h:
.. figure:: /recipes/figures/russell18jgr/Fig7h_fgco2_zonal-flux.png
   :align:   center
   :width: 50%

   Figure 7h: the time and zonal mean flux of CO\ :sub:`2`\  in gC/(yr * m\ :sup:`2`\) plotted as a line plot.


.. _fig_russell_7i:
.. figure::  /recipes/figures/russell18jgr/Fig7i_fgco2_integrated-flux.png
   :align:   center
   :width: 50%

   Figure 7i is the cumulative integral of the net CO\ :sub:`2`\  flux from 90S to 30S (in PgC/yr) plotted as a line plot. 

.. _fig_russell_8:
.. figure::  /recipes/figures/russell18jgr/Fig8_polar-ph.png
   :align:   center
   :width: 50%

   Figure 8: Annual-mean surface pH plotted as a polar contour map.

.. _fig_russell_9a:
.. figure::  /recipes/figures/russell18jgr/Fig9a.png
   :align:   center
   :width: 50%

   Figure 9a: Scatter plot of the width of the Southern Hemisphere westerly wind band (in degrees of latitude) against the annual-mean integrated heat uptake south of 30S (in PW—negative uptake is heat lost from the ocean) along with the best fit line.

.. _fig_russell_9b:
.. figure::  /recipes/figures/russell18jgr/Fig9b.png
   :align:   center
   :width: 50%

   Figure 9b: Scatter plot of the width of the Southern Hemisphere westerly wind band (in degrees of latitude) against the annual-mean integrated carbon uptake south of 30S (in Pg C/yr), along with the best fit line.

.. _fig_russell_9c:
.. figure:: /recipes/figures/russell18jgr/Fig9c.png
   :align:   center
   :width: 50%

   Figure 9c: Scatter plot of the net heat uptake south of 30S (in PW) against the annual-mean integrated carbon uptake south of 30S (in Pg C/yr), along with the best fit line.
.. _recipes_snowalbedo:

Emergent constraint on snow-albedo effect
=========================================

Overview
--------

The recipe recipe_snowalbedo.yml computes the springtime snow-albedo
feedback values in climate change versus springtime values in the seasonal
cycle in transient climate change experiments following Hall and Qu (2006).
The strength of the snow-albedo effect is quantified by the variation in net
incoming shortwave radiation (Q) with surface air temperature (T\ :sub:`s`\) due
to changes in surface albedo :math:`\alpha_s`:

.. math::

   \left( \frac{\partial Q}{\partial T_s} \right) = -I_t \cdot \frac{\partial \alpha_p}{\partial \alpha_s} \cdot \frac{\Delta \alpha_s}{\Delta T_s}

The diagnostic produces scatterplots of simulated springtime
:math:`\Delta \alpha_s`/:math:`\Delta T_s` values in climate change (ordinate)
vs. simulated springtime :math:`\Delta \alpha_s`/:math:`\Delta T_s` values in the
seasonal cycle (abscissa).

Ordinate values: the change in April :math:`\alpha_s` (future projection - historical)
averaged over NH land masses poleward of 30°N is divided by the change in
April T\ :sub:`s` (future projection - historical) averaged over the same region.
The change in :math:`\alpha_s` (or T\ :sub:`s`) is defined as the difference between
22nd-century-mean :math:`\alpha_s`: (T\ :sub:`s`) and 20th-century-mean :math:`\alpha_s`. Values of
:math:`\alpha_s` are weighted by April incoming insolation (I\ :sub:`t`) prior to averaging.

Abscissa values: the seasonal cycle :math:`\Delta \alpha_s`/:math:`\Delta T_s`
values, based on 20th century climatological means, are calculated by
dividing the difference between April and May :math:`\alpha_s`: averaged over NH continents
poleward of 30°N by the difference between April and May T\ :sub:`s` averaged over the
same area. Values of :math:`\alpha_s`: are weighted by April incoming insolation prior to
averaging.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_snowalbedo.yml

Diagnostics are stored in diag_scripts/emergent_constraints/

    * snowalbedo.ncl: springtime snow-albedo feedback values vs. seasonal cycle


User settings in recipe
-----------------------

#. Script snowalbedo.ncl

   *Required settings for script*

   * exp_presentday: name of present-day experiment (e.g. "historical")
   * exp_future: name of climate change experiment (e.g. "rcp45")

   *Optional settings for script*

   * diagminmax: observational uncertainty (min and max)
   * legend_outside: create extra file with legend (true, false)
   * styleset: e.g. "CMIP5" (if not set, this diagnostic will create its own
     color table and symbols for plotting)
   * suffix: string to be added to output filenames
   * xmax: upper limit of x-axis (default = automatic)
   * xmin: lower limit of x-axis (default = automatic)
   * ymax: upper limit of y-axis (default = automatic)
   * ymin: lower limit of y-axis (default = automatic)

   *Required settings for variables*

   * ref_model: name of reference data set

   *Optional settings for variables*

   none

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*


Variables
---------

* tas (atmos, monthly mean, longitude latitude time)
* rsdt (atmos, monthly mean, longitude latitude time)
* rsuscs, rsdscs (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

* ERA-Interim (tas - esmvaltool/utils/cmorizers/obs/cmorize_obs_ERA-Interim.ncl)
* ISCCP-FH (rsuscs, rsdscs, rsdt - esmvaltool/utils/cmorizers/obs/cmorize_obs_isccp_fh.ncl)


References
----------

* Flato, G., J. Marotzke, B. Abiodun, P. Braconnot, S.C. Chou, W. Collins, P.
  Cox, F. Driouech, S. Emori, V. Eyring, C. Forest, P. Gleckler, E. Guilyardi,
  C. Jakob, V. Kattsov, C. Reason and M. Rummukainen, 2013: Evaluation of
  Climate Models. In: Climate Change 2013: The Physical Science Basis.
  Contribution of Working Group I to the Fifth Assessment Report of the
  Intergovernmental Panel on Climate Change [Stocker, T.F., D. Qin, G.-K.
  Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex and
  P.M. Midgley (eds.)]. Cambridge University Press, Cambridge, United Kingdom
  and New York, NY, USA.

* Hall, A., and X. Qu, 2006: Using the current seasonal cycle to constrain
  snow albedo feedback in future climate change, Geophys. Res. Lett., 33,
  L03502, doi:10.1029/2005GL025127.


Example plots
-------------

.. figure:: /recipes/figures/flato13ipcc/fig-9-45a.png
   :align: center

   Scatterplot of springtime snow-albedo effect values in climate
   change vs. springtime :math:`\Delta \alpha_s`/:math:`\Delta T_s` values in
   the seasonal cycle in transient climate change experiments (CMIP5 historical
   experiments: 1901-2000, RCP4.5 experiments: 2101-2200). Similar to IPCC AR5
   Chapter 9 (Flato et al., 2013), Figure 9.45a.
.. _recipe_autoassess_landsurface_permafrost.rst:

Land-surface Permafrost - Autoassess diagnostics
================================================

Overview
--------

Permafrost thaw is an important impact of climate change, and is the source of
a potentially strong Earth system feedback through the release of soil carbon
into the atmosphere. This recipe provides metrics that evaluate the
climatological performance of models in simulating soil temperatures that
control permafrost. Performance metrics (with observation-based estimates in brackets):

* permafrost area (17.46 million square km)
* fractional area of permafrost northwards of zero degree isotherm (0.47)
* soil temperature at 1m minus soil temperature at surface (-0.53 degrees C)
* soil temperature at surface minus air temperature (6.15 degrees C)
* annual amplitude at 1m / annual amplitude at the surface (0.40 unitless)
* annual amplitude at the surface / annual air temperature (0.57 unitless)


Plots:

* Maps of permafrost extent and zero degC isotherm
* Normalised assessment metrics plot comparing control and experiment

The recipe takes as input a control model and experimental model, comparisons being made
with these two models.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_autoassess_landsurface_permafrost.yml

Diagnostics are stored in esmvaltool/diag_scripts/autoassess/

    * autoassess_area_base.py: wrapper for autoassess scripts
    * land_surface_permafrost/permafrost.py: script to calculate permafrost
      metrics
    * plot_autoassess_metrics.py: plot normalised assessment metrics


User settings in recipe
-----------------------

#. Script autoassess_area_base.py

   *Required settings for script*

   * area: must equal land_surface_permafrost to select this diagnostic
   * control_model: name of model to be used as control
   * exp_model: name of model to be used as experiment
   * start: date (YYYY/MM/DD) at which period begins (see note on time gating)
   * end: date (YYYY/MM/DD) at which period ends (see note on time gating)
   * climfiles_root: path to observation climatologies

   *Optional settings for script*

   * title: arbitrary string with name of diagnostic
   * obs_models: unused for this recipe

   *Required settings for variables*

   none

   *Optional settings for variables*

   none


#. Script plot_autoassess_metrics.py

   *Required settings for script*

   * area: must equal land_surface_permafrost to select this diagnostic
   * control_model: name of model to be used as control in metrics plot
   * exp_model: name of model to be used as experiment in metrics plot
   * title: string to use as plot title

   *Optional settings for script*

   none

   *Required settings for variables*

   none

   *Optional settings for variables*

   none



Variables
---------

* tas (atmos, monthly mean, longitude latitude time)
* tsl (land, monthly mean, longitude latitude time)
* mrsos (land, monthly mean, longitude latitude time)
* sftlf (mask, fixed, longitude latitude)


Observations and reformat scripts
---------------------------------

None


References
----------

* Observed permafrost extent is from http://nsidc.org/data/ggd318.html: Brown, J.,
  O. Ferrians, J. A. Heginbottom, and E. Melnikov. 2002. Circum-Arctic Map of
  Permafrost and Ground-Ice Conditions, Version 2. Boulder, Colorado USA. NSIDC:
  National Snow and Ice Data Center.  When calculating the global area of
  permafrost the grid cells are weighted by the proportion of permafrost within
  them.

* Annual mean air temperature is from: Legates, D. R., and C. J. Willmott, 1990:
  Mean seasonal and spatial variability in global surface air temperature. Theor.
  Appl. Climatol., 41, 11-21.  The annual mean is calculated from the seasonal
  mean data available at the Met Office.

* The soil temperature metrics are calcuated following: Charles D. Koven, William
  J. Riley, and Alex Stern, 2013: Analysis of Permafrost Thermal Dynamics and
  Response to Climate Change in the CMIP5 Earth System Models. J. Climate, 26. 
  (Table 3) http://dx.doi.org/10.1175/JCLI-D-12-00228.1 The
  locations used for Table 3 were extracted from the model and the modelled
  metrics calculated.


Example plots
-------------

.. figure:: /recipes/figures/autoassess_landsurface/pf_extent_north_america_ACCESS-CM2.png
   :scale: 50 %
   :alt: pf_extent_north_america_ACCESS-CM2.png

   Permafrost extent and zero degC isotherm, showing North America

.. figure:: /recipes/figures/autoassess_landsurface/pf_extent_asia_ACCESS-CM2.png
   :scale: 50 %
   :alt: pf_extent_asia_ACCESS-CM2.png

   Permafrost extent and zero degC isotherm, showing Asia and Europe

.. figure:: /recipes/figures/autoassess_landsurface/Permafrost_Metrics.png
   :scale: 50 %
   :alt: Permafrost_Metrics.png

   Normalised metrics plot comparing a control and experiment simulation


Additional notes on usage
-------------------------
The ``landsurface_permafrost`` area metric is part of the ``esmvaltool/diag_scripts/autoassess`` diagnostics,
and, as any other ``autoassess`` metric, it uses the ``autoassess_area_base.py`` as general purpose
wrapper. This wrapper accepts a number of input arguments that are read through from the recipe.

This recipe is part of the larger group of Autoassess metrics ported to ESMValTool
from the native Autoassess package from the UK's Met Office. The ``diagnostics`` settings
are almost the same as for the other Autoassess metrics.

.. note::

   **Time gating for autoassess metrics.**

   To preserve the native Autoassess functionalities,
   data loading and selection on time is done somewhat
   differently for ESMValTool's autoassess metrics: the
   time selection is done in the preprocessor as per usual but
   a further time selection is performed as part of the diagnostic.
   For this purpose the user will specify a ``start:`` and ``end:``
   pair of arguments of ``scripts: autoassess_script`` (see below
   for example). These are formatted as ``YYYY/MM/DD``; this is
   necessary since the Autoassess metrics are computed from 1-Dec
   through 1-Dec rather than 1-Jan through 1-Jan. This is a temporary
   implementation to fully replicate the native Autoassess functionality
   and a minor user inconvenience since they need to set an extra set of
   ``start`` and ``end`` arguments in the diagnostic; this will be phased
   when all the native Autoassess metrics have been ported to ESMValTool
   review has completed.


An example of standard inputs as read by ``autoassess_area_base.py`` and passed
over to the diagnostic/metric is listed below.

.. code-block:: yaml

    scripts:
      plot_landsurf_permafrost: &plot_landsurf_permafrost_settings
        <<: *autoassess_landsurf_permafrost_settings
        control_model: MPI-ESM-LR
        exp_model: MPI-ESM-MR
        script: autoassess/plot_autoassess_metrics.py
        ancestors: ['*/autoassess_landsurf_permafrost']
        title: "Plot Land-Surface Permafrost Metrics"
        plot_name: "Permafrost_Metrics"
        diag_tag: aa_landsurf_permafrost
        diag_name: autoassess_landsurf_permafrost
.. _recipes_combined_indices:

Nino indices, North Atlantic Oscillation (NAO), Souther Oscillation Index (SOI)
===============================================================================

Overview
--------

The goal of this diagnostic is to compute indices based on area averages.

In recipe_combined_indices.yml, after defining the period (historical or
future projection), the variable is selected. The predefined areas are:

* Nino 3
* Nino 3.4
* Nino 4
* North Atlantic Oscillation (NAO)
* Southern Oscillation Index (SOI)

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_combined_indices.yml

Diagnostics are stored in diag_scripts/magic_bsc/

* combined_indices.R : calculates the area-weighted means and multi-model means, with or without weights



User settings
-------------

User setting files are stored in recipes/

#.	recipe_combined_indices.yml

   *Required settings for script*

   * region: one of the following strings Nino3, Nino3.4, Nino4, NAO, SOI
   * running_mean: an integer specifying the length of the window (in months) to be used for computing the running mean.
   * moninf: an integer can be given to determine the first month of the seasonal mean to be computed (from 1 to 12, corresponding to January to December respectively).
   * monsup: an integer specifying the last month to be computed (from 1 to 12, corresponding to January to December respectively).
   * standardized: ‘true’ or ‘false’ to specify whether to compute the standarization of the variable.


     *Required settings for preprocessor (only for 3D variables)*
     
	  extract_levels:
   *   levels: [50000] # e.g. for 500 hPa level
   *   scheme: nearest
   
Variables
---------

* all variables (atmos/ocean, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Trenberth, Kevin & National Center for Atmospheric Research Staff (Eds). Last modified 11 Jan 2019. "The Climate Data Guide: Nino SST Indices (Nino 1+2, 3, 3.4, 4; ONI and TNI)." Retrieved from https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni.


Example plots
-------------

.. _fig_combinedindices1:
.. figure::  /recipes/figures/Combined_Indices_Area_Average/Nino3.4_tos_Dec-Feb_running-mean__1950-2005.png
   :align:   center
   :width:   14cm

Time series of the standardized sea surface temperature (tos) area averaged over the Nino 3.4 region during the boreal winter (December-January-February). The time series correspond to the MPI-ESM-MR (red) and BCC-CSM1-1 (blue) models and their mean (black) during the period 1950-2005 for the ensemble r1p1i1 of the historical simulations... _recipes_shapeselect:

Shapeselect
===========

Overview
--------
Impact modelers are often interested in data for irregular regions best defined by a shapefile. With the shapefile selector tool, the user can extract time series or CII data for a user defined region. The region is defined by a user provided shapefile that includes one or several polygons. For each polygon, a new timeseries, or CII, is produced with only one time series per polygon. The spatial information is reduced to a representative point for the polygon ('representative') or as an average of all grid points within the polygon boundaries ('mean_inside'). If there are no grid points strictly inside the polygon, the 'mean_inside' method defaults to 'representative' for that polygon. An option for displaying the grid points together with the shapefile polygon allows the user to assess which method is most optimal. In case interpolation to a high input grid is necessary, this can be provided in a pre-processing stage. Outputs are in the form of a NetCDF file, or as ascii code in csv format.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_shapeselect.yml


Diagnostics are stored in diag_scripts/shapeselect/

    * diag_shapeselect.py: calculate the average of grid points inside the
      user provided shapefile and returns the result as a NetCDF or Excel sheet.


User settings in recipe
-----------------------

#. Script diag_shapeselect.py

   *Required settings (scripts)*

   * shapefile: path to the user provided shapefile. A relative path is relative to the auxiliary_data_dir as configured in config-user.yml.

   * weighting_method: the preferred weighting method 'mean_inside' - mean of all grid points inside polygon; 'representative' - one point inside or close to the polygon is used to represent the complete area.

   * write_xlsx: true or false to write output as Excel sheet or not.

   * write_netcdf: true or false to write output as NetCDF or not.

Variables
---------

* pr,tas      (daily)

Example plots
-------------

.. _fig_shapeselect:
.. figure::  /recipes/figures/shapeselect/shapeselect.png
   :align:   center
   :width:   14cm

   Example of the selection of model grid points falling within (blue pluses) and without (red dots) a provided shapefile (blue contour).
.. _recipes_hyint:

Hydroclimatic intensity and extremes (HyInt)
============================================


Overview
--------
The HyInt tool calculates a suite of hydroclimatic and climate extremes indices to perform a multi-index evaluation of climate models. The tool firstly computes a set of 6 indices that allow to evaluate the response of the hydrological cycle to global warming with a joint view of both wet and dry extremes. The indices were selected following Giorgi et al. (2014) and include the simple precipitation intensity index (SDII) and extreme precipitation index (R95), the maximum dry spell length (DSL) and wet spell length (WSL), the hydroclimatic intensity index (HY-INT), which is a measure of the overall behaviour of the hydroclimatic cycle (Giorgi et al., 2011), and the precipitation area (PA), i.e. the area over which at any given day precipitation occurs, (Giorgi et al., 2014). Secondly, a selection of the 27 temperature and precipitation -based indices of extremes from the Expert Team on Climate Change Detection and Indices (ETCCDI) produced by the climdex (https://www.climdex.org) library can be ingested to produce a multi-index analysis. The tool allows then to perform a subsequent analysis of the selected indices calculating timeseries and trends over predefined continental areas, normalized to a reference period. Trends are calculated using the R `lm` function and significance testing performed with a Student T test on non-null coefficients hypothesis. Trend coefficients are stored together with their statistics which include standard error, t value and Pr(>|t|). The tool can then produce a variety of types of plots including global and regional maps, maps of comparison between models and a reference dataset, timeseries with their spread, trend lines and summary plots of trend coefficients.

The hydroclimatic indices calculated by the recipe_hyint.yml and included in the output are defined as follows:

* PRY = mean annual precipitation
* INT = mean annual precipitation intensity (intensity during wet days, or simple precipitation intensity index SDII)
* WSL = mean annual wet spell length (number of consecutive days during each wet spell)
* DSL = mean annual dry spell lenght (number of consecutive days during each dry spell)
* PA  = precipitation area (area over which of any given day precipitation occurs)
* R95 = heavy precipitation index (percent of total precipitation above the 95% percentile of the reference distribution)
* HY-INT = hydroclimatic intensity. HY-INT = normalized(INT) x normalized(DSL).

The recipe_hyint_extreme_events.yml includes an additional call to the :ref:`recipes_extreme_events` diagnostics, which allows to calculate the ETCCDI indices and include them in the subsequent analysis together with the hydroclimatic indices. All of the selected indices are then stored in output files and figures.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_hyint.yml (evaluating the 6 hydroclimatic indices, performing trend analysis and plotting)
* recipe_hyint_extreme_events.yml (similar to the recipe_hyint.yml but with an additional call to the :ref:`recipes_extreme_events` diagnostic for calculation of ETCCDI indices and inclusion of them in the trend analysis and plotting)

Diagnostics are stored in diag_scripts/hyint/

* hyint.R

and subroutines

* hyint_diagnostic.R
* hyint_functions.R
* hyint_parameters.R
* hyint_plot_trends.R
* hyint_etccdi_preproc.R
* hyint_metadata.R
* hyint_plot_maps.R
* hyint_preproc.R
* hyint_trends.R

See details of the extreme_events diagnostics under recipe_extreme_events.yml.

Known issues
------------

*recipe_hyint_extreme_events.yml*

Call to the :ref:`recipes_extreme_events` diagnostic requires the ncdf4.helpers library, which is currently unavailable on CRAN. Users need therefore to install the library manually, e.g. through the following commands to download the package tarball from CRAN archive, install it and remove the package tarball:

  * url <- "https://cran.r-project.org/src/contrib/Archive/ncdf4.helpers/ncdf4.helpers_0.3-3.tar.gz"
  * pkgFile <- "ncdf4.helpers_0.3-3.tar.gz"
  * download.file(url = url, destfile = pkgFile)
  * install.packages(pkgs=pkgFile, type="source", repos=NULL)
  * unlink(pkgFile)

User settings
-------------

*Required settings for script*

* norm_years: first and last year of reference normalization period to be used for normalized indices

* select_indices: indices to be analysed and plotted. Select one or more fields from the following list (order-sensitive): "pa_norm", "hyint",  "int_norm", "r95_norm", "wsl_norm", "dsl_norm", "int", "dsl", "wsl"

* select_regions: Select regions for timeseries and maps from the following list: GL=Globe, GL60=Global 60S/60N, TR=Tropics (30S/30N), SA=South America, AF=Africa, NA=North America, IN=India, EU=Europe, EA=East-Asia, AU=Australia

* plot_type: type of figures to be plotted. Select one or more from: 1=lon/lat maps per individual field/exp/multi-year mean, 2=lon/lat maps per individual field exp-ref-diff/multi-year mean, 3=lon/lat maps multi-field/exp-ref-diff/multi-year mean, 11=timeseries over required individual region/exp, 12=timeseries over multiple regions/exp, 13=timeseries with multiple models, 14=summary trend coefficients multiple regions, 15=summary trend coefficients multiple models


*Additional settings for recipe_hyint_extreme_events.yml*

* call to the extreme_events diagnostics: see details in recipe_extreme_events.yml. Make sure that the base_range for extreme_events coincides with the norm_range of hyint and that all ETCCDI indices that are required to be imported in hyint are calculated by the extreme_events diagnostics.

* etccdi_preproc: set to true to pre-process and include ETCCDI indices in hyint

* etccdi_list_import: specify the list of ETCCDI indices to be imported, e.g.: "tn10pETCCDI", "tn90pETCCDI", "tx10pETCCDI", "tx90pETCCDI"

* select_indices: this required settings should here be revised to include the imported indices, e.g.: "pa_norm", "hyint", "tn10pETCCDI", "tn90pETCCDI", "tx10pETCCDI", "tx90pETCCDI"


*Optional settings for script (with default setting)*

#. Data

   * rgrid (false): Define whether model data should be regridded. (a) false to keep original resolution; (b) set desired regridding resolution in cdo format e.g., "r320x160"; (c) "REF" to use resolution of reference model

#. Plotting

   * npancol (2): number of columns in timeseries/trends multipanel figures
   * npanrow (3): number of rows in timeseries/trends multipanel figures
   * autolevels (true): select automated (true) or pre-set (false) range of values in plots
   * autolevels_scale (1): factor multiplying automated range for maps and timeseries
   * autolevels_scale_t (1.5): factor multiplying automated range for trend coefficients

#. Maps

   * oplot_grid (false): plot grid points over maps
   * boxregion (false): !=0 plot region boxes over global maps with thickness = abs(boxregion); white (>0) or grey (<0).
   * removedesert (false) remove (flag as NA) grid points with mean annual pr < 0.5 mm/day (deserts, Giorgi2014). This affects timeseries and trends calculations too.

#. Timeseries and trends

   * weight_tseries (true): adopt area weights in timeseries
   * trend_years (false): (a) false = apply trend to all years in dataset; (b) [year1, year2] to apply trend calculation and plotting only to a limited time interval
   * add_trend (true): add linear trend to plot
   * add_trend_sd (false): add dashed lines of stdev range to timeseries
   * add_trend_sd_shade (false): add shade of stdev range to timeseries
   * add_tseries_lines (true): plot lines connecting timeseries points
   * add_zeroline (true): plot a dashed line at y=0
   * trend_years_only (false): limit timeseries plotting to the time interval adopted for trend calculation (excluding the normalization period)
   * scale100years (true): plot trends scaled as 1/100 years
   * scalepercent (false): plot trends as percent change


Variables
---------

* pr (atmos, daily mean, longitude latitude time)

*Additional variables for recipe_hyint_extreme_events.yml*

* tas (atmos, daily mean, longitude latitude time)
* tasmin (atmos, daily mean, longitude latitude time)
* tasmax (atmos, daily mean, longitude latitude time)

Observations and reformat scripts
---------------------------------

None.


References
----------

* Giorgi et al., 2014, J. Geophys. Res. Atmos., 119, 11,695–11,708, doi:10.1002/ 2014JD022238
* Giorgi et al., 2011, J. Climate 24, 5309-5324, doi:10.1175/2011JCLI3979.1


Example plots
-------------

.. figure:: /recipes/figures/hyint/hyint_maps.png
   :width: 10cm

   Mean hydroclimatic intensity for the EC-EARTH model, for the historical + RCP8.5 projection in the period 1976-2099

.. figure:: /recipes/figures/hyint/hyint_timeseries.png
   :width: 12cm

   Timeseries for multiple indices and regions for the ACCESS1-0 model, for the historical + RCP8.5 projection in the period 1976-2099, normalized to the 1976-2005 historical period.

.. figure:: /recipes/figures/hyint/hyint_trends.png
   :width: 12cm

   Multi-model trend coefficients over selected indices for CMIP5 models in the RCP8.5 2006-2099 projection, normalized to the 1976-2005 historical period.
.. _recipes_meehl20sciadv:

Context for interpreting equilibrium climate sensitivity and transient climate response from the CMIP6 Earth system models
==========================================================================================================================

Overview
--------

This recipe reproduces the analysis of `Meehl et al., Sci. Adv. (2020)`_. In
this paper, the equilibrium climate sensitivity (ECS) and transient climate
response (TCR) are evaluated for the CMIP6 models and put into historical
context.

.. _`Meehl et al., Sci. Adv. (2020)`: https://advances.sciencemag.org/content/6/26/eaba1981


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_meehl20sciadv.yml

Diagnostics are stored in diag_scripts/

   * climate_metrics/ecs.py
   * climate_metrics/tcr.py
   * climate_metrics/create_table.py
   * ipcc_ar5/ch09_fig09_42b.py


User settings in recipe
-----------------------

* Script climate_metrics/ecs.py

   See :ref:`here<ecs.py>`.


* Script climate_metrics/tcr.py

   See :ref:`here<tcr.py>`.


* Script climate_metrics/create_table.py

   * ``calculate_mean``, *bool*, optional (default: ``True``): Calculate
     mean over all datasets and add it to table.
   * ``calculate_std``, *bool*, optional (default: ``True``): Calculate
     standard deviation over all datasets and add it to table.
   * ``exclude_datasets``, *list of str*, optional (default:
     ``['MultiModelMean']``): Exclude certain datasets when calculating
     statistics over all datasets and for assigning an index.
   * ``patterns``, *list of str*, optional: Patterns to filter list of input
     data.
   * ``round_output``, *int*, optional: If given, round output to given
     number of decimals.


* Script ipcc_ar5/ch09_fig09_42b.py

   See :ref:`here<ch09_fig09_42b.py>`.


Variables
---------

* *rlut* (atmos, monthly, longitude, latitude, time)
* *rsdt* (atmos, monthly, longitude, latitude, time)
* *rsut* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)



References
----------

* Meehl, G. A., Senior, C. A., Eyring, V., Flato, G., Lamarque, J. F.,
  Stouffer, R. J., Taylor, K. E. and Schlund, M., *Context for interpreting
  equilibrium climate sensitivity and transient climate response from the CMIP6
  Earth system models*, Science Advances, 6(26), eaba1981,
  `<https://doi.org/10.1126/sciadv.aba1981>`_, 2020.


Example plots
-------------

.. _fig_meehl20sciadv_1:
.. figure:: /recipes/figures/meehl20sciadv/cmip6_gregory_regression.png
   :align: center
   :width: 50%

   ECS calculated for the CMIP6 models using the Gregory method over different
   time scales. Using the entire 150-year 4xCO2 experiment (black line), there
   is an ECS value of 3.8 K; using only the first 20 years (blue dots and blue
   line), there is an ECS of 3.4 K; and using the last 130 years, there is an
   ECS of 4.1 K (orange dots and orange line).

.. _fig_meehl20sciadv_2:
.. figure:: /recipes/figures/meehl20sciadv/cmip6_tcr_vs_ecs.png
   :align: center
   :width: 50%

   TCR as a function of ECS for the CMIP6 models (black line is a linear fit).
   The :math:`R^2` values are given in the upper left parts of each panel. The
   numbers denote individual CMIP6 models.
.. _recipe_mpqb_xch4:

Diagnostics of integrated atmospheric methane (XCH4)
====================================================

Overview
--------

This recipe ``recipe_mpqb_xch4.yml`` allows the comparison of integrated atmospheric methane
between CMIP6 model simulations and observations, and produces lineplots of monthly mean
methane values, annual cycles and annual growth rates:

* Monthly mean time series of XCH4 for pre-defined regions (global, Northern Hemisphere, Southern Hemisphere)
* Annual cycles of XCH4 for pre-defined regions (global, Northern Hemisphere, Southern Hemisphere)
* Annual growth rates of XCH4 for pre-defined regions (global, Northern Hemisphere, Southern Hemisphere)

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/mpqb/

* recipe_mpqb_xch4.yml

Diagnostics are stored in esmvaltool/diag_scripts/mpqb/

* mpqb_lineplot.py
* mpqb_lineplot_anncyc.py
* mpqb_lineplot_growthrate.py

Observations and reformat scripts
---------------------------------

Observations used in this recipe are:

* CDS-XCH4 (ESA CCI dataset served on the Copernicus Climate data store)

A cmorizing script for this dataset is available (``cmorize_obs_cds_xch4.ncl``).

XCH4 is a derived variable that needs to be calculated from four different variables (ch4, hus, zg, ps).
A derivation script is included in the ESMValCore.


User settings in recipe
-----------------------
#. Preprocessor

   * ``pp_lineplots_xx_mon``: Regridding, masking all missing values from all used datasets, area-mean ('xx' can ge replaced by 'gl'=global, 'sh'=southern hemisphere, 'nh'=northern hemisphere), units converted to [ppbv] to obtain one time series of monthly mean values for the selected region (global, southern hemisphere, northern hemisphere)
   * ``pp_lineplots_xx_ann``: Regridding, masking all missing values from all used datasets, area-mean ('xx' can ge replaced by 'gl'=global, 'sh'=southern hemisphere, 'nh'=northern hemisphere), units converted to [ppbv] to obtain one time series of annual mean values for the selected region (global, southern hemisphere, northern hemisphere)   
   * ``pp_lineplots_anncyc_xx:`` : Regridding, masking all missing values from all used datasets, area-mean ('xx' can ge replaced by 'gl'=global, 'sh'=southern hemisphere, 'nh'=northern hemisphere), units converted to [ppbv], monthly climate statistics applied to one annual cycle for the whole chosen time period and for the selected region (global, southern hemisphere, northern hemisphere)
   * ``xch4_def_xx``: defining the time period over which the analysis should be calculated; options are "cmip6" which overlapping period of the observations and the CMIP6 historical simulations, and "future" which covers the time period of CMIP6 scenarios

#. Additional needed files
   
   * ``mpqb_cfg_xch4.yml``: In this file additional information for the used datasets are defined and stored, e.g. alias of the dataset name and the color that is used to display the dataset in the figures
   * ``mpqb_utils.yml``: In this file the preparations for the dataset displays are made.

#. Script <mpqb_lineplot.py>

   *Required settings for script*

   * no additional settings required

   *Optional settings for script*
   
   * no optional settings available

   *Required settings for variables*
   
   * no settings for the variable required

#. Script <mpqb_lineplot_anncyc.py>

   *Required settings for script*

   * no additional settings required

   *Optional settings for script*
   
   * no optional settings available

   *Required settings for variables*
   
   * no settings for the variable required

#. Script <mpqb_lineplot_growthrate.py>

   *Required settings for script*

   * no additional settings required

   *Optional settings for script*
   
   * no optional settings available

   *Required settings for variables*
   
   * no settings for the variable required


Variables
---------

* ch4 (atmos, monthly mean, longitude latitude level time)
* hus (atmos, monthly mean, longitude latitude level time)
* zg (atmos, monthly mean, longitude latitude level time)
* ps (atmos, monthly mean, longitude latitude time)

All variables are necessary to calculate the derived variable xch4.


Example plots
-------------

.. _lineplot_xch4_2003-2014_monmean: 
.. figure::  /recipes/figures/mpqb/lineplot_xch4_2003-2014_monmean.png
   :align:   center

   Monthly mean time series of XCH4, calculated over the whole globe, for individual CMIP6 model simulations.
   
   
   
.. _recipes_schlund20jgr:

Constraining uncertainty in projected gross primary production (GPP) with machine learning
==========================================================================================

Overview
--------

These recipes reproduce the analysis of `Schlund et al., JGR: Biogeosciences
(2020)`_. In this paper, a machine learning regression (MLR) approach (using
the MLR algorithm `Gradient Boosted Regression Trees, GBRT`_) is proposed to
constrain uncertainties in projected gross primary production (GPP) in the RCP
8.5 scenario using observations of process-based diagnostics.

.. _`Gradient Boosted Regression Trees, GBRT`: https://scikit-learn.org/stable/modules/ensemble.html#gradient-tree-boosting
.. _`Schlund et al., JGR: Biogeosciences (2020)`: https://doi.org/10.1029/2019JG005619


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * schlund20jgr/recipe_schlund20jgr_gpp_abs_rcp85.yml
   * schlund20jgr/recipe_schlund20jgr_gpp_change_1pct.yml
   * schlund20jgr/recipe_schlund20jgr_gpp_change_rcp85.yml

Diagnostics are stored in diag_scripts/

   * :ref:`mlr/evaluate_residuals.py<api.esmvaltool.diag_scripts.mlr.evaluate_residuals>`
   * :ref:`mlr/main.py<api.esmvaltool.diag_scripts.mlr.main>`
   * :ref:`mlr/mmm.py<api.esmvaltool.diag_scripts.mlr.mmm>`
   * :ref:`mlr/plot.py<api.esmvaltool.diag_scripts.mlr.plot>`
   * :ref:`mlr/postprocess.py<api.esmvaltool.diag_scripts.mlr.postprocess>`
   * :ref:`mlr/preprocess.py<api.esmvaltool.diag_scripts.mlr.preprocess>`
   * :ref:`mlr/rescale_with_emergent_constraint.py<api.esmvaltool.diag_scripts.mlr.rescale_with_emergent_constraint>`

General information (including an example and more details) on machine learning
regression (MLR) diagnostics is given
:ref:`here<api.esmvaltool.diag_scripts.mlr.models>`. The API documentation is
available :ref:`here<api.esmvaltool.diag_scripts.mlr>`.


Variables
---------

* *co2s* (atmos, monthly, longitude, latitude, time)
* *gpp* (land, monthly, longitude, latitude, time)
* *gppStderr* (land, monthly, longitude, latitude, time)
* *lai* (land, monthly, longitude, latitude, time)
* *pr* (atmos, monthly, longitude, latitude, time)
* *rsds* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* CRU_ (*pr*, *tas*)
* ERA-Interim_ (*rsds*)
* LAI3g_ (*lai*)
* MTE_ (*gpp*, *gppStderr*)
* Scripps-CO2-KUM_ (*co2s*)

.. _CRU: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/cruts.1811131722.v4.02/
.. _ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda/
.. _LAI3g: http://cliveg.bu.edu/modismisr/lai3g-fpar3g.html
.. _MTE: http://www.bgc-jena.mpg.de/geodb/BGI/Home
.. _Scripps-CO2-KUM: https://scrippsco2.ucsd.edu/data/atmospheric_co2/kum.html


References
----------

* Schlund, M., Eyring, V., Camps‐Valls, G., Friedlingstein, P., Gentine, P., &
  Reichstein, M. (2020). Constraining uncertainty in projected gross primary
  production with machine learning. Journal of Geophysical Research:
  Biogeosciences, 125, e2019JG005619,
  `<https://doi.org/10.1029/2019JG005619>`_.


Example plots
-------------

.. _fig_schlund20jgr_1:
.. figure:: /recipes/figures/schlund20jgr/map_prediction_output___GBRT_change.png
   :align: center
   :width: 50%

   GBRT-based prediction of the fractional GPP change over the 21st century (=
   GPP(2091-2100) / GPP(1991-2000)).

.. _fig_schlund20jgr_2:
.. figure:: /recipes/figures/schlund20jgr/map_prediction_output_error___GBRT_change.png
   :align: center
   :width: 50%

   Corresponding error of the GBRT-based prediction of the fractional GPP
   change over the 21st century (considering errors in the MLR model and errors
   in the predictors).

.. _fig_schlund20jgr_3:
.. figure:: /recipes/figures/schlund20jgr/map_prediction_output___GBRT_abs.png
   :align: center
   :width: 50%

   GBRT-based prediction of the absolute GPP at the end of the 21st century
   (2091-2100).

.. _fig_schlund20jgr_4:
.. figure:: /recipes/figures/schlund20jgr/map_prediction_output_error___GBRT_abs.png
   :align: center
   :width: 50%

   Corresponding error of the GBRT-based prediction of the absolute GPP at the
   end of the 21st century (considering errors in the MLR model and errors in
   the predictors).

.. _fig_schlund20jgr_5:
.. figure:: /recipes/figures/schlund20jgr/rmse_plot.png
   :align: center
   :width: 50%

   Boxplot of the root mean square error of prediction (RMSEP) distributions
   for six different statistical models used to predict future absolute GPP
   (2091-2100) using a leave-one-model-out cross-validation approach. The
   distribution for each statistical model contains seven points (black dots,
   one for each climate model used as truth) and is represented in the
   following way: the lower and upper limit of the blue boxes correspond to the
   25% and 75% quantiles, respectively. The central line in the box shows the
   median, the black "x" the mean of the distribution. The whiskers outside the
   box represent the range of the distribution

.. _fig_schlund20jgr_6:
.. figure:: /recipes/figures/schlund20jgr/feature_importance.png
   :align: center
   :width: 50%

   Global feature importance of the GBRT model for prediction of the absolute
   GPP at the end of the 21st century (2091-2100).

.. _fig_schlund20jgr_7:
.. figure:: /recipes/figures/schlund20jgr/residuals_distribution.png
   :align: center
   :width: 50%

   Distribution of the residuals of the GBRT model for the prediction of
   absolute GPP at the end of the 21st century (2091-2100) for the training
   data (blue) and test data excluded from training (green).

.. _fig_schlund20jgr_8:
.. figure:: /recipes/figures/schlund20jgr/training_progress.png
   :align: center
   :width: 50%

   Training progress of the GBRT model for the prediction of absolute GPP at
   the end of the 21st century (2091-2100) evaluated as normalized root mean
   square error on the training data (blue) and test data excluded from
   training (green).
.. _recipe_diurnal_temperature_index:

Diurnal temperature range
=========================

Overview
--------

The goal of this diagnostic is to compute a vulnerability indicator for the diurnal temperature range (DTR); the maximum variation in temperature within a period of 24 hours at a given location.  This indicator was first proposed by the energy sector, to identify locations which may experience increased diurnal temperature variation in the future, which would put additional stress on the operational management of district heating systems. This indicator was defined as the DTR exceeding 5 degrees celsius at a given location and day of the year (Deandreis et al., N.D.). Projections of this indicator currently present high uncertainties, uncertainties associated to both Tmax and Tmin in future climate projections.

As well as being of use to the energy sector, the global‐average DTR has been evaluated using both observations and climate model simulations (Braganza et. al., 2004) and changes in the mean and variability of the DTR have been shown to have a wide range of impacts on society, such as on the transmission of diseases (Lambrechts et al., 2011;  Paaijmans et al., 2010).

The recipe recipe_diurnal_temperature_index.yml computes first a mean DTR for a reference period using historical simulations and then, the number of days when the DTR from the future climate projections exceeds that of the reference period by 5 degrees or more. The user can define both the reference and projection periods, and the region to be considered.  The output produced by this recipe consists of a four panel plot showing the maps of the projected mean DTR indicator for each season and a netcdf file containing the corresponding data.



Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_diurnal_temperature_index.yml

Diagnostics are stored in diag_scripts/magic_bsc/

* diurnal_temp_index.R : calculates the diaurnal temperature vulnerability index.


User settings
-------------

User setting files are stored in recipes/

#. recipe_diurnal_temperature_index.yml

   *Required settings for script*

   * None

Variables
---------

* tasmin and tasmax (atmos, daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Amiri, S. (2013). Economic and Environmental Benefits of CHP-based District Heating Systems in Sweden. Retrieved from http://www.sgc.se/ckfinder/userfiles/files/sokmotor/LiU67.pdf

* Braganza, K., Karoly, D. J., & Arblaster, J. M. (2004). Diurnal temperature range as an index of global climate change during the twentieth century. Geophysical Research Letters, 31(13), n/a – n/a. https://doi.org/10.1029/2004GL019998

* Déandreis, C. (IPSL), Braconnot, P. (IPSL), and Planton, S.(CNRMGAME)(2014). Impact du changement climatique sur la gestion des réseaux de chaleur. DALKIA, Étude réalisée pour l'entreprise DALKIA. Last access 24.02.2021. https://docplayer.fr/9496504-Impact-du-changement-climatique-sur-la-gestion-des-reseaux-de-chaleur.html

* Lambrechts, L., Paaijmans, K. P., Fansiri, T., Carrington, L. B., Kramer, L. D., Thomas, M. B., & Scott, T. W. (2011). Impact of daily temperature fluctuations on dengue virus transmission by Aedes aegypti. Proceedings of the National Academy of Sciences of the United States of America, 108(18), 7460–7465. https://doi.org/10.1073/pnas.1101377108

* Paaijmans, K. P., Blanford, S., Bell, A. S., Blanford, J. I., Read, A. F., & Thomas, M. B. (2010). Influence of climate on malaria transmission depends on daily temperature variation. Proceedings of the National Academy of Sciences of the United States of America, 107(34), 15135–15139. https://doi.org/10.1073/pnas.1006422107

* Kalnay, E., & Cai, M. (2003). Impact of urbanization and land-use change on climate. Nature, 423(6939), 528–531. https://doi.org/10.1038/nature01675

* Thyholt, M., & Hestnes, A. G. (2008). Heat supply to low-energy buildings in district heating areas: Analyses of CO2 emissions and electricity supply security. Energy and Buildings, 40(2), 131–139. https://doi.org/10.1016/J.ENBUILD.2007.01.016

Example plots
-------------

.. _fig_diurnal:
.. figure::  /recipes/figures/diurnal_temp_index/Seasonal_DTRindicator_MPI-ESM-MR_2030_2080_1961_1990.png
   :align:   center
   :width:   14cm

Mean number of days exceeding the Diurnal Temperature Range (DTR) simulated during the historical period (1961-1990) by 5 degrees during the period 2030-2080. The result is derived from one RCP 8.5 scenario simulated by MPI-ESM-MR.
.. _recipes_flato13ipcc:

IPCC AR5 Chapter 9 (selected figures)
=====================================

Overview
--------

The goal of this recipe is to collect diagnostics to reproduce Chapter 9 of AR5,
so that the plots can be readily reproduced and compared to previous CMIP
versions. In this way we can next time start with what was available in the
previous round and can focus on developing more innovative methods of analysis
rather than constantly having to "re-invent the wheel".

The plots are produced collecting the diagnostics from individual recipes. The
following figures from Flato et al. (2013) can currently be reproduced:

    * Figure 9.2 a,b,c: Annual-mean surface air temperature for the period
      1980-2005. a) multi-model mean, b) bias as the difference between the
      CMIP5 multi-model mean and the climatology from ERA-Interim
      (Dee et al., 2011), c) mean absolute model error with respect to the
      climatology from ERA-Interim.

    * Figure 9.3: Seasonality (December-January-February minus June-July-August)
      of surface (2 m) air temperature (°C) for the period 1980-2005.
      (a) Multi-model mean for the historical experiment. (b) Multi-model mean
      of absolute seasonality. (c) Difference between the multi-model mean
      and the ERA-Interim reanalysis seasonality. (d) Difference between the
      multi-model mean and the ERA-Interim absolute seasonality.

    * Figure 9.4: Annual-mean precipitation rate (mm day-1) for the period
      1980-2005. a) multi-model mean, b) bias as the difference between the
      CMIP5 multi-model mean and the climatology from the Global Precipitation
      Climatology Project (Adler et al., 2003), c) multi-model mean absolute
      error with respect to observations, and d) multi-model mean error
      relative to the multi-model mean precipitation ifself.

    * Figure 9.5: Climatological (1985-2005) annual-mean cloud radiative
      effects in Wm-2 for the CMIP5 models against CERES EBAF (2001-2011) in
      Wm-2. Top row shows the shortwave effect; middle row the longwave effect,
      and bottom row the net effect. Multi-model-mean biases against CERES
      EBAF 2.6 are shown on the left, whereas the right panels show zonal
      averages from CERES EBAF 2.6 (black), the individual CMIP5 models (thin
      gray lines), and the multi-model mean (thick red line).

    * Figure 9.6: Centred pattern correlations between models and observations
      for the annual mean climatology over the period 1980–1999. Results are
      shown for individual CMIP3 (black) and CMIP5 (blue) models as thin
      dashes, along with the corresponding ensemble average (thick dash) and
      median (open circle). The four variables shown are surface air
      temperature (TAS), top of the atmosphere (TOA) outgoing longwave
      radiation (RLUT), precipitation (PR) and TOA shortwave cloud radiative
      effect (SW CRE). The correlations between the reference and alternate
      observations are also shown (solid green circles).

    * Figure 9.8: Observed and simulated time series of the anomalies in annual
      and global mean surface temperature. All anomalies are differences from
      the 1961-1990 time-mean of each individual time series. The reference
      period 1961-1990 is indicated by yellow shading; vertical dashed grey
      lines represent times of major volcanic eruptions. Single simulations
      for CMIP5 models (thin lines); multi-model mean (thick red line);
      different observations (thick black lines). Dataset pre-processing like
      described in Jones et al., 2013.

    * Figure 9.14: Sea surface temperature plots for zonal mean error, equatorial
      (5 deg north to 5 deg south) mean error, and multi model mean for zonal error
      and equatorial mean.

    * Figure 9.24: Time series of (a) Arctic and (b) Antarctic sea ice extent;
      trend distributions of (c) September Arctic and (d) February Antarctic
      sea ice extent.

    * Figure 9.26: Ensemble-mean global ocean carbon uptake (a) and global land
      carbon uptake (b) in the CMIP5 ESMs for the historical period 1900–2005.
      For comparison, the observation-based estimates provided by the Global
      Carbon Project (GCP) are also shown (thick black line). The confidence
      limits on the ensemble mean are derived by assuming that the CMIP5 models
      are drawn from a t-distribution. The grey areas show the range of annual mean
      fluxes simulated across the model ensemble. This figure includes results
      from all CMIP5 models that reported land CO2 fluxes, ocean CO2 fluxes, or
      both (Anav et al., 2013).

    * Figure 9.27: Simulation of global mean (a) atmosphere–ocean CO2 fluxes
      ("fgCO2") and (b) net atmosphere–land CO2 fluxes ("NBP"), by ESMs for the
      period 1986–2005. For comparison, the observation-based estimates
      provided by Global Carbon Project (GCP) and the Japanese Meteorological
      Agency (JMA) atmospheric inversion are also shown. The error bars for the
      ESMs and observations represent interannual variability in the fluxes,
      calculated as the standard deviation of the annual means over the period
      1986–2005.

    * Figure 9.42a: Equilibrium climate sensitivity (ECS) against the global
      mean surface air temperature, both for the period 1961-1990 and for the
      pre-industrial control runs.

    * Figure 9.42b: Transient climate response (TCR) against equilibrium climate
      sensitivity (ECS).

    * Figure 9.45a: Scatterplot of springtime snow-albedo effect values in climate
      change vs. springtime d(alpha\ :sub:`s`\)/d(T\ :sub:`s`\) values in the seasonal
      cycle in transient climate change experiments (Hall and Qu, 2006).

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_flato13ipcc.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * carbon_cycle/main.ncl: See :ref:`here<recipes_anav13jclim>`.
    * climate_metrics/ecs.py: See :ref:`here<ecs.py>`.
    * clouds/clouds_bias.ncl: global maps of the multi-model mean and the multi-model
      mean bias (Fig. 9.2, 9.4)
    * clouds/clouds_isccp: global maps of multi-model mean minus observations + zonal
      averages of individual models, multi-model mean and observations (Fig. 9.5)
    * ipcc_ar5/ch09_fig09_3.ncl: multi-model mean seasonality of near-surface
      temperature (Fig. 9.3)
    * ipcc_ar5/ch09_fig09_6.ncl: calculating pattern correlations of annual mean
      climatologies for one variable (Fig 9.6 preprocessing)
    * ipcc_ar5/ch09_fig09_6_collect.ncl: collecting pattern correlation for each
      variable and plotting correlation plot (Fig 9.6)
    * ipcc_ar5/tsline.ncl: time series of the global mean (anomaly) (Fig. 9.8)
    * ipcc_ar5/ch09_fig09_14.py: Zonally averaged and equatorial SST (Fig. 9.14)
    * seaice/seaice_tsline.ncl: Time series of sea ice extent (Fig. 9.24a/b)
    * seaice/seaice_trends.ncl: Trend distributions of sea ice extent (Fig 9.24c/d)
    * ipcc_ar5/ch09_fig09_42a.py: ECS vs. surface air temperature (Fig. 9.42a)
    * ipcc_ar5/ch09_fig09_42b.py: TCR vs. ECS (Fig. 9.42b)
    * emergent_constraints/snowalbedo.ncl: snow-albedo effect (Fig. 9.45a)

User settings in recipe
-----------------------

#. Script carbon_cycle/main.ncl

   See :ref:`here<recipes_anav13jclim>`.

#. Script climate_metrics/ecs.py

   See :ref:`here<ecs.py>`.

#. Script clouds/clouds_bias.ncl

#. Script clouds_bias.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * plot_abs_diff: additionally also plot absolute differences (true, false)
   * plot_rel_diff: additionally also plot relative differences (true, false)
   * projection: map projection, e.g., Mollweide, Mercator
   * timemean: time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean)

   * Required settings (variables)*

   * reference_dataset: name of reference datatset

   *Optional settings (variables)*

   * long_name: description of variable

   *Color tables*

   * variable "tas": diag_scripts/shared/plot/rgb/ipcc-tas.rgb,
     diag_scripts/shared/plot/rgb/ipcc-tas-delta.rgb
   * variable "pr-mmday": diag_scripts/shared/plots/rgb/ipcc-precip.rgb,
     diag_scripts/shared/plot/rgb/ipcc-precip-delta.rgb

#. Script clouds/clouds_ipcc.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * explicit_cn_levels: contour levels
   * mask_ts_sea_ice: true = mask T < 272 K as sea ice (only for variable "ts");
     false = no additional grid cells masked for variable "ts"
   * projection: map projection, e.g., Mollweide, Mercator
   * styleset: style set for zonal mean plot ("CMIP5", "DEFAULT")
   * timemean: time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean)
   * valid_fraction: used for creating sea ice mask (mask_ts_sea_ice = true):
     fraction of valid time steps required to mask grid cell as valid data

   *Required settings (variables)*

   * reference_dataset:  name of reference data set

   *Optional settings (variables)*

   * long_name: description of variable
   * units: variable units

   *Color tables*

   * variables "pr", "pr-mmday": diag_scripts/shared/plot/rgb/ipcc-precip-delta.rgb

#. Script ipcc_ar5/tsline.ncl

   *Required settings for script*

   * styleset: as in diag_scripts/shared/plot/style.ncl functions

   *Optional settings for script*

   * time_avg: type of time average (currently only "yearly" and "monthly" are
     available).
   * ts_anomaly: calculates anomalies with respect to the defined period; for
     each gird point by removing the mean for the given calendar month
     (requiring at least 50% of the data to be non-missing)
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * ref_value: if true, right panel with mean values is attached
   * ref_mask: if true, model fields will be masked by reference fields
   * region: name of domain
   * plot_units: variable unit for plotting
   * y-min: set min of y-axis
   * y-max: set max of y-axis
   * mean_nh_sh: if true, calculate first NH and SH mean
   * volcanoes: if true, lines of main volcanic eruptions will be added
   * run_ave: if not equal 0 than calculate running mean over this number of
     years
   * header: if true, region name as header

   *Required settings for variables*

   none

   *Optional settings for variables*

   * reference_dataset: reference dataset; REQUIRED when calculating
     anomalies

   *Color tables*

   * e.g. diag_scripts/shared/plot/styles/cmip5.style

#. Script ipcc_ar5/ch09_fig09_3.ncl

   *Required settings for script*

   none

   *Optional settings for script*

   * projection: map projection, e.g., Mollweide, Mercator (default = Robinson)

   *Required settings for variables*

   * reference_dataset: name of reference observation

   *Optional settings for variables*

   * map_diff_levels: explicit contour levels for plotting

#. Script ipcc_ar5/ch09_fig09_6.ncl

   *Required settings for variables*

   * reference_dataset: name of reference observation

   *Optional settings for variables*

   * alternative_dataset: name of alternative observations

#. Script ipcc_ar5/ch09_fig09_6_collect.ncl

   *Required settings for script*

   none

   *Optional settings for script*

   * diag_order: List of diagnostic names in the order variables
     should appear on x-axis

#. Script seaice/seaice_trends.ncl

   *Required settings (scripts)*

   * month: selected month (1, 2, ..., 12) or annual mean ("A")
   * region: region to be analyzed ( "Arctic" or "Antarctic")

   *Optional settings (scripts)*

   * fill_pole_hole: fill observational hole at North pole, Default: False

   *Optional settings (variables)*

   * ref_model: array of references plotted as vertical lines

#. Script seaice/seaice_tsline.ncl

   *Required settings (scripts)*

   * region: Arctic, Antarctic
   * month: annual mean (A), or month number (3 = March, for Antarctic; 9 = September for Arctic)

   *Optional settings (scripts)*

   * styleset: for plot_type cycle only (cmip5, cmip6, default)
   * multi_model_mean: plot multi-model mean and standard deviation (default: False)
   * EMs_in_lg: create a legend label for individual ensemble members (default: False)
   * fill_pole_hole: fill polar hole (typically in satellite data) with sic = 1 (default: False)

#. Script ipcc_ar5/ch09_fig09_42a.py

   *Required settings for script*

   none

   *Optional settings for script*

   * axes_functions: :obj:`dict` containing methods executed for the plot's
     :class:`matplotlib.axes.Axes` object.
   * dataset_style: name of the style file (located in
     :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
   * matplotlib_style: name of the matplotlib style file (located in
     :mod:`esmvaltool.diag_scripts.shared.plot.styles_python.matplotlib`).
   * save: :obj:`dict` containing keyword arguments for the function
     :func:`matplotlib.pyplot.savefig`.
   * seaborn_settings: Options for :func:`seaborn.set` (affects all plots).

.. _ch09_fig09_42b.py:

#. Script ipcc_ar5/ch09_fig09_42b.py

   *Required settings for script*

   none

   *Optional settings for script*

   * dataset_style: Dataset style file (located in
     :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`). The entry
     ``marker`` is ignored when ``marker_file`` is given.
   * log_x: Apply logarithm to X axis (ECS).
   * log_y: Apply logarithm to Y axis (TCR).
   * marker_column: Name of the column to look up markers in ``marker_file``.
   * marker_file: CSV file with markers (can also be integers). Must have the
     columns ``dataset`` and ``marker`` (or the column specified by
     ``marker_column``).  If a relative path is given, assumes that this is a
     pattern to search for ancestor files.
   * savefig_kwargs: Keyword arguments for :func:`matplotlib.pyplot.savefig`.
   * seaborn_settings: Options for :func:`seaborn.set` (affects all plots).
   * x_lim: Plot limits for X axis (ECS).
   * y_lim: Plot limits for Y axis (TCR).

#. Script emergent_constraints/snowalbedo.ncl

   *Required settings for script*

   * exp_presentday: name of present-day experiment (e.g. "historical")
   * exp_future: name of climate change experiment (e.g. "rcp45")

   *Optional settings for script*

   * diagminmax: observational uncertainty (min and max)
   * legend_outside: create extra file with legend (true, false)
   * styleset: e.g. "CMIP5" (if not set, this diagnostic will create its own
     color table and symbols for plotting)
   * suffix: string to be added to output filenames
   * xmax: upper limit of x-axis (default = automatic)
   * xmin: lower limit of x-axis (default = automatic)
   * ymax: upper limit of y-axis (default = automatic)
   * ymin: lower limit of y-axis (default = automatic)

   *Required settings for variables*

   * ref_model: name of reference data set

   *Optional settings for variables*

   none

Variables
---------

* areacello (fx, longitude latitude)
* fgco2 (ocean, monthly mean, longitude latitude time)
* nbp (ocean, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* rlut, rlutcs (atmos, monthly mean, longitude latitude time)
* rsdt (atmos, monthly mean, longitude latitude time)
* rsuscs, rsdscs (atmos, monthly mean, longitude latitude time)
* rsut, rsutcs (atmos, monthly mean, longitude latitude time)
* sic (ocean-ice, monthly mean, longitude latitude time)
* tas (atmos, monthly mean, longitude latitude time)
* tos (ocean, monthly mean, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
instructions.*

* CERES-EBAF (rlut, rlutcs, rsut, rsutcs - obs4MIPs)
* ERA-Interim (tas, ta, ua, va, zg, hus - esmvaltool/cmorizers/obs/cmorize_obs_ERA-Interim.ncl)
* GCP2018 (fgco2, nbp - esmvaltool/cmorizers/obs/cmorize_obs_gcp2018.py)
* GPCP-SG (pr - obs4MIPs)
* JMA-TRANSCOM (fgco2, nbp - esmvaltool/cmorizers/obs/cmorize_obs_jma_transcom.py)
* HadCRUT4 (tas - esmvaltool/cmorizers/obs/cmorize_obs_hadcrut4.ncl)
* HadISST (sic, tos - esmvaltool/cmorizers/obs/cmorize_obs_hadisst.ncl)
* ISCCP-FH (rsuscs, rsdscs, rsdt - esmvaltool/cmorizers/obs/cmorize_obs_isccp_fh.ncl)


References
----------

* Flato, G., J. Marotzke, B. Abiodun, P. Braconnot, S.C. Chou, W. Collins, P.
  Cox, F. Driouech, S. Emori, V. Eyring, C. Forest, P. Gleckler, E. Guilyardi,
  C. Jakob, V. Kattsov, C. Reason and M. Rummukainen, 2013: Evaluation of
  Climate Models. In: Climate Change 2013: The Physical Science Basis.
  Contribution of Working Group I to the Fifth Assessment Report of the
  Intergovernmental Panel on Climate Change [Stocker, T.F., D. Qin, G.-K.
  Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex and
  P.M. Midgley (eds.)]. Cambridge University Press, Cambridge, United Kingdom
  and New York, NY, USA.

* Hall, A., and X. Qu, 2006: Using the current seasonal cycle to constrain
  snow albedo feedback in future climate change, Geophys. Res. Lett., 33,
  L03502, doi:10.1029/2005GL025127.

* Jones et al., 2013: Attribution of observed historical near-surface temperature
  variations to anthropogenic and natural causes using CMIP5 simulations. Journal
  of Geophysical Research: Atmosphere, 118, 4001-4024, doi:10.1002/jgrd.50239.


Example plots
-------------

.. figure::  /recipes/figures/flato13ipcc/fig-9-2.png
   :align:   center

   Figure 9.2 a,b,c: Annual-mean surface air temperature for the period
   1980-2005. a) multi-model mean, b) bias as the difference between the
   CMIP5 multi-model mean and the climatology from ERA-Interim
   (Dee et al., 2011), c) mean absolute model error with respect to the
   climatology from ERA-Interim.

.. figure::  /recipes/figures/flato13ipcc/fig-9-3.png
   :align:   center

   Figure 9.3: Multi model values for seasonality of near-surface temperature,
   from top left to bottom right: mean, mean of absolute seasonality, mean bias
   in seasonality, mean bias in absolute seasonality. Reference dataset:
   ERA-Interim.

.. figure::  /recipes/figures/flato13ipcc/fig-9-4.png
   :align:   center

   Figure 9.4: Annual-mean precipitation rate (mm day-1) for the period
   1980-2005. a) multi-model mean, b) bias as the difference between the
   CMIP5 multi-model mean and the climatology from the Global Precipitation
   Climatology Project (Adler et al., 2003), c) multi-model mean absolute
   error with respect to observations, and d) multi-model mean error
   relative to the multi-model mean precipitation ifself.

.. figure::  /recipes/figures/flato13ipcc/fig-9-5.png
   :align:   center

   Figure 9.5: Climatological (1985-2005) annual-mean cloud radiative
   effects in Wm-2 for the CMIP5 models against CERES EBAF (2001-2011) in
   Wm-2. Top row shows the shortwave effect; middle row the longwave effect,
   and bottom row the net effect. Multi-model-mean biases against CERES
   EBAF 2.6 are shown on the left, whereas the right panels show zonal
   averages from CERES EBAF 2.6 (black), the individual CMIP5 models (thin
   gray lines), and the multi-model mean (thick red line).

.. figure::  /recipes/figures/flato13ipcc/fig-9-6.png
   :align:   center

   Figure 9.6: Centred pattern correlations between models and observations
   for the annual mean climatology over the period 1980–1999. Results are
   shown for individual CMIP3 (black) and CMIP5 (blue) models as thin
   dashes, along with the corresponding ensemble average (thick dash) and
   median (open circle). The four variables shown are surface air
   temperature (TAS), top of the atmosphere (TOA) outgoing longwave
   radiation (RLUT), precipitation (PR) and TOA shortwave cloud radiative
   effect (SW CRE). The correlations between the reference and alternate
   observations are also shown (solid green circles).

.. figure::  /recipes/figures/flato13ipcc/fig-9-8.png
   :align:   center

   Figure 9.8: Observed and simulated time series of the anomalies in annual
   and global mean surface temperature. All anomalies are differences from
   the 1961-1990 time-mean of each individual time series. The reference
   period 1961-1990 is indicated by yellow shading; vertical dashed grey
   lines represent times of major volcanic eruptions. Single simulations
   for CMIP5 models (thin lines); multi-model mean (thick red line);
   different observations (thick black lines). Dataset pre-processing like
   described in Jones et al., 2013.

.. figure:: /recipes/figures/flato13ipcc/fig-9-14.png
   :align: center

   Figure 9.14: (a) Zonally averaged sea surface temperature (SST) error
   in CMIP5 models. (b) Equatorial SST error in CMIP5 models. (c) Zonally
   averaged multi-model mean SST error for CMIP5 together with
   inter-model standard deviation (shading). (d) Equatorial multi-model
   mean SST in CMIP5 together with inter-model standard deviation
   (shading) and observations (black).  Model climatologies are derived
   from the 1979-1999 mean of the historical simulations. The Hadley
   Centre Sea Ice and Sea Surface Temperature (HadISST) (Rayner et
   al., 2003) observational climatology for 1979-1999 is used as a
   reference for the error calculation (a), (b), and (c); and for
   observations in (d).

.. figure::  /recipes/figures/seaice/trend_sic_extend_Arctic_September_histogram.png
   :align:   center
   :width:   9cm

   Figure 9.24c: Sea ice extent trend distribution for the Arctic in September.

.. figure::  /recipes/figures/seaice/extent_sic_Arctic_September_1960-2005.png
   :align:   center
   :width:   12cm

   Figure 9.24a: Time series of total sea ice area and extent (accumulated) for the Arctic
   in September including multi-model mean and standard deviation.

.. figure:: /recipes/figures/flato13ipcc/fig-9-26.png
   :align: center

   Figure 9.26 (bottom): Ensemble-mean global land carbon uptake in the CMIP5
   ESMs for the historical period 1900–2005.  For comparison, the
   observation-based estimates provided by the Global Carbon Project (GCP) are
   also shown (black line). The confidence limits on the ensemble mean are
   derived by assuming that the CMIP5 models come from a t-distribution. The
   grey areas show the range of annual mean fluxes simulated across the model
   ensemble.

.. figure:: /recipes/figures/flato13ipcc/fig-9-27.png
   :align: center

   Figure 9.27 (top): Simulation of global mean atmosphere–ocean CO2 fluxes
   ("fgCO2") by ESMs for the period 1986–2005. For comparison, the
   observation-based estimates provided by Global Carbon Project (GCP) are also
   shown. The error bars for the ESMs and observations represent interannual
   variability in the fluxes, calculated as the standard deviation of the
   annual means over the period 1986–2005.

.. figure:: /recipes/figures/flato13ipcc/fig-9-42a.png
   :align: center

   Figure 9.42a: Equilibrium climate sensitivity (ECS) against the global mean
   surface air temperature of CMIP5 models, both for the period 1961-1990
   (larger symbols) and for the pre-industrial control runs (smaller symbols).

.. figure:: /recipes/figures/flato13ipcc/fig-9-42b.png
   :align: center

   Figure 9.42b: Transient climate response (TCR) against equilibrium climate
   sensitivity (ECS) for CMIP5 models.

.. figure:: /recipes/figures/flato13ipcc/fig-9-45a.png
   :align: center

   Figure 9.45a: Scatterplot of springtime snow-albedo effect values in climate
   change vs. springtime :math:`\Delta \alpha_s`/:math:`\Delta T_s` values in
   the seasonal cycle in transient climate change experiments (CMIP5 historical
   experiments: 1901-2000, RCP4.5 experiments: 2101-2200).
.. _recipes_wenzel16nat:

Projected land photosynthesis constrained by changes in the seasonal cycle of atmospheric CO\ :sub:`2`
======================================================================================================

Overview
--------

Selected figures from `Wenzel et al. (2016)`_ are reproduced with recipe_wenzel16nat.yml. Gross primary productivity (gpp) and atmospheric CO\ :sub:`2` concentrations at the surface  (co2s) are analyzed for the carbon cycle - concentration feedback in the historical (esmHistorical) and uncoupled (esmFixCLim1, here the carbon cycle is uncoupled to the climate response) simulations. The recipe includes a set of routines to diagnose the long-term carbon cycle - concentration feedback parameter (beta) from an ensemble of CMIP5 models and the observable change in the CO\ :sub:`2` seasonal cycle amplitude due to rising atmospheric CO\ :sub:`2` levels. As a key figure of this recipe, the diagnosed values from the models beta vs. the change in CO\ :sub:`2` amplitude are compared in a scatter plot constituting an emergent constraint.

.. _`Wenzel et al. (2016)`: https://www.nature.com/articles/nature19772

Available recipe and diagnostics
-----------------------------------

Recipes are stored in recipes/

    * recipe_wenzel16nat.yml

Diagnostics are stored in diag_scripts/carbon_ec/

    * carbon_beta: (1) scatter plot of annual gpp vs. annual CO\ :sub:`2` and
      (2) barchart of gpp(2xCO\ :sub:`2`)/gpp(1xCO\ :sub:`2`); calculates beta
      for emergent constraint (carbon_co2_cycle.ncl)
    * carbon_co2_cycle.ncl: (1) scatter plot of CO\ :sub:`2` amplitude vs.
      annual CO\ :sub:`2`, (2) barchart of sensitivity of CO\ :sub:`2` amplitude
      to CO\ :sub:`2`, (3) emergent constraint:
      gpp(2xCO\ :sub:`2`)/gpp(1xCO\ :sub:`2`) vs. sensitivity of CO\ :sub:`2`
      amplitude to CO\ :sub:`2`, (4) probability density function of constrained
      and unconstrained sensitivity of CO\ :sub:`2` amplitude to CO\ :sub:`2`


User settings
-------------

#. Script carbon_beta.ncl

   *Required Settings (scripts)*

   * styleset: project style for lines, colors and symbols

   *Optional Settings (scripts)*

   * bc_xmax_year: end year to calculate beta (default: use last available year of all models)
   * bc_xmin_year: start year to calculate beta (default: use first available year of all models)

   *Required settings (variables)*

   none

   *Optional settings (variables)*

   none

#. Script carbon_co2_cycle.ncl 

   *Required Settings (scripts)*

   * nc_infile: path of netCDF file containing beta (output from carbon_beta.ncl)
   * styleset: project style for lines, colors and symbols

   *Optional Settings (scripts)*

   * bc_xmax_year: end year (default = last year of all model datasets available)
   * bc_xmin_year: start year (default = first year of all model datasets available)

   *Required settings (variables)*

   * reference_dataset: name of reference datatset (observations)

   *Optional settings (variables)*

   none


Variables
---------

* co2s (atmos, monthly mean, plev longitude latitude time)
* gpp (land, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

* ESRL: Earth System Research Laboratory, ground-based CO\ :sub:`2` measurements


References
----------

* Wenzel, S., Cox, P., Eyring, V. et al., 2016, Projected land photosynthesis constrained by changes in the seasonal cycle of atmospheric CO\ :sub:`2`. Nature 538, 499501, doi: doi.org/10.1038/nature19772


Example plots
-------------

.. figure:: /recipes/figures/wenzel16nat/fig_1.png
   :width: 12 cm 
   :align: center
   
   Comparison of CO\ :sub:`2` seasonal amplitudes for CMIP5 historical simulations and observations showing annual mean atmospheric CO\ :sub:`2` versus the amplitudes of the CO\ :sub:`2` seasonal cycle at Pt. Barrow, Alaska (produced with carbon_co2_cycle.ncl, similar to Fig. 1a from Wenzel et al. (2016)).
      
.. figure:: /recipes/figures/wenzel16nat/fig_2.png
   :width: 12 cm 
   :align: center
   
   Barchart showing the gradient of the linear correlations for the comparison of CO\ :sub:`2` seasonal amplitudes for CMIP5 historical for at Pt. Barrow, Alaska (produced with carbon_co2_cycle.ncl, similar to Fig. 1b from Wenzel et al. (2016)).

.. figure:: /recipes/figures/wenzel16nat/fig_3.png
   :width: 12 cm
   :align: center

   Emergent constraint on the relative increase of large-scale GPP for a doubling of CO\ :sub:`2`, showing the correlations between the sensitivity of the CO\ :sub:`2` amplitude to annual mean CO\ :sub:`2` increases at Pt. Barrow (x-axis) and the high-latitude (60N - 90N) CO\ :sub:`2` fertilization on GPP at 2xCO\ :sub:`2`. The red line shows the linear best fit of the regression together with the prediction error (orange shading), the gray shading shows the observed range (produced with carbon_co2_cycle.ncl, similar to Fig. 3a from Wenzel et al. (2016)).
.. _recipes:

Recipes
-------

.. toctree::
   :maxdepth: 1

Atmosphere
^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   recipe_miles
   recipe_clouds
   recipe_cmug_h2o
   recipe_crem
   recipe_combined_climate_extreme_index
   recipe_consecdrydays
   recipe_deangelis15nat
   recipe_diurnal_temperature_index
   recipe_eady_growth_rate
   recipe_extreme_events
   recipe_eyring06jgr
   recipe_eyring13jgr
   recipe_gier20bg
   recipe_heatwaves_coldwaves
   recipe_hyint
   recipe_impact
   recipe_modes_of_variability
   recipe_mpqb_xch4
   recipe_quantilebias
   recipe_bock20jgr
   recipe_spei
   recipe_martin18grl
   recipe_autoassess_stratosphere
   recipe_autoassess_landsurface_permafrost
   recipe_autoassess_landsurface_surfrad
   recipe_autoassess_landsurface_soilmoisture
   recipe_zmnam
   recipe_thermodyn_diagtool
   recipe_validation
   recipe_radiation_budget

Climate metrics
^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   recipe_perfmetrics
   recipe_smpi

Future projections
^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   recipe_climwip
   recipe_li17natcc
   recipe_schlund20jgr
   recipe_meehl20sciadv
   recipe_emergent_constraints
   recipe_wenzel14jgr
   recipe_schlund20esd
   recipe_cox18nature
   recipe_snowalbedo
   recipe_ecs
   recipe_kcs
   recipe_wenzel16jclim
   recipe_wenzel16nat
   recipe_tcr

IPCC
^^^^
.. toctree::
   :maxdepth: 1

   recipe_flato13ipcc
   recipe_collins13ipcc

Land
^^^^
.. toctree::
   :maxdepth: 1

   recipe_albedolandcover
   recipe_carvalhais14nat
   recipe_hydrology
   recipe_hydro_forcing
   recipe_landcover
   recipe_anav13jclim
   recipe_runoff_et

Ocean
^^^^^
.. toctree::
   :maxdepth: 1

   recipe_arctic_ocean
   recipe_cvdp
   recipe_combined_indices
   recipe_esacci_oc
   recipe_oceans
   recipe_sea_surface_salinity
   recipe_russell18jgr

Other
^^^^^
.. toctree::
   :maxdepth: 1

   recipe_examples
   recipe_capacity_factor
   recipe_cmorizers
   recipe_ensclus
   recipe_esacci_lst
   recipe_multimodel_products
   recipe_rainfarm
   recipe_pv_capacity_factor
   recipe_seaice_feedback
   recipe_seaice
   recipe_seaice_drift
   recipe_shapeselect
   recipe_toymodel
.. _recipe_cmorizers:

CMORizer recipes
=================

Overview
--------

These are CMORizer recipes calling CMORizer diagnostic scripts.

ESMValCore supports ERA5 hourly and monthly datasets in their native
format, see :ref:`CMORization as a fix <esmvaltool:cmorization_as_fix>`
and `ERA5 data documentation <https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation>`_.
It may be useful in some cases to create ERA5 daily CMORized data. This can be
achieved by using a CMORizer *recipe*,
see `recipe_daily_era5.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/cmorizers/recipe_daily_era5.yml>`_.
This recipe reads native, hourly ERA5 data, performs a daily aggregation
preprocessor, and then calls a diagnostic that operates on the data. In this
example, the diagnostic renames the files to the standard OBS6 file names. The output
are thus daily, CMORized ERA5 data, that can be used through the OBS6 project.
As such, this example recipe creates a local pool of CMORized data. The advantage, in this
case, is that the daily aggregation is performed only once, which can save a lot
of time and compute if it is used often.

The example CMORizer recipe can be run like any other ESMValTool recipe:

.. code-block:: bash

    esmvaltool run cmorizers/recipe_daily_era5.yml

Note that the ``recipe_daily_era5.yml`` adds the next day of the new year to
the input data. This is because one of the fixes needed for the ERA5 data is to
shift the time axis of non-instantaneous variables half an hour back in time, resulting in a missing
record on the last day of the year. ERA5 data can be downloaded using `era5cli <https://era5cli.readthedocs.io>`_.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * cmorizers/recipe_daily_era5.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * cmorizers/era5.py: generates output filename


User settings in recipe
-----------------------

#. cmorizers/recipe_daily_era5.yml

   *Required add_one_day preprocessor settings:*

   * start_year: 1990
   * start_month: 1
   * start_day: 1
   * end_year: 1991
   * end_month: 1
   * end_day: 1

These settings should not be changed
   * daily_mean:
         operator: mean
   * daily_min:
         operator: min
   * daily_max:
         operator: max

Variables
---------

#. cmorizers/recipe_daily_era5.yml

   * clt
   * evspsbl
   * evspsblpot
   * mrro
   * pr
   * prsn
   * ps
   * psl
   * rlds
   * rls
   * rsds
   * rsdt
   * rss
   * tas
   * tasmax
   * tasmin
   * tdps
   * ts
   * tsn
   * uas
   * vas

References
----------

* Hersbach, H., et al., Quarterly Journal of the Royal Meteorological Society, 730, 1999-2049, doi:10.1002/qj.3803, 2020.
.. _recipe_autoassess_landsurface_soilmoisture.rst:

Land-surface Soil Moisture - Autoassess diagnostics
===================================================

Overview
--------

Soil moisture is a critical component of the land system, controling surface
energy fluxes in many areas of the world. This recipe provides metrics that
evaluate the skill of models' spatial and seasonal distribution of soil
moisture against the ESA CCI soil moisture ECV.

Performance metrics:

* median absolute error (model minus observations)

Metrics are calculated using model and observation multi-year climatologies (seasonal means) 
for meteorological seasons:
* December-January-February (djf)
* March-April-May (mam)
* June-July-August (jja)
* September-October-November (son)

Plots:

* Normalised assessment metrics plot comparing control and experiment

The recipe takes as input a control model and experimental model, comparisons being made
with these two models.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_autoassess_landsurface_soilmoisture.yml

Diagnostics are stored in esmvaltool/diag_scripts/autoassess/

    * autoassess_area_base.py: wrapper for autoassess scripts
    * land_surface_soilmoisture/soilmoisture.py: script to calculate soil moisture
      metrics
    * plot_autoassess_metrics.py: plot normalised assessment metrics


User settings in recipe
-----------------------

#. Script autoassess_area_base.py

   *Required settings for script*

   * area: must equal land_surface_soilmoisture to select this diagnostic
   * control_model: name of model to be used as control
   * exp_model: name of model to be used as experiment
   * start: date (YYYY/MM/DD) at which period begins (see note on time gating)
   * end: date (YYYY/MM/DD) at which period ends (see note on time gating)
   * climfiles_root: path to observation climatologies

   *Optional settings for script*

   * title: arbitrary string with name of diagnostic
   * obs_models: unused for this recipe

   *Required settings for variables*

   none

   *Optional settings for variables*

   none


#. Script plot_autoassess_metrics.py

   *Required settings for script*

   * area: must equal land_surface_soilmoisture to select this diagnostic
   * control_model: name of model to be used as control in metrics plot
   * exp_model: name of model to be used as experiment in metrics plot
   * title: string to use as plot title

   *Optional settings for script*

   none

   *Required settings for variables*

   none

   *Optional settings for variables*

   none


Variables
---------

* mrsos (land, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

1999-2008 climatologies (seasonal means) from ESA ECV Soil Moisture Dataset v1.
Produced by the ESA CCI soil moisture project: https://www.esa-soilmoisture-cci.org/node/93


References
----------
* Dorigo, W.A., Wagner, W., Albergel, C., Albrecht, F.,  Balsamo, G., Brocca, L., Chung, D., Ertl, M., Forkel, M., Gruber, A., Haas, E., Hamer, D. P. Hirschi, M., Ikonen, J., De Jeu, R. Kidd, R.  Lahoz, W., Liu, Y.Y., Miralles, D., Lecomte, P. (2017).  ESA CCI Soil Moisture for improved Earth system understanding: State-of-the art and future directions. In Remote Sensing of Environment, 2017,  ISSN 0034-4257, https://doi.org/10.1016/j.rse.2017.07.001.

* Gruber, A., Scanlon, T., van der Schalie, R., Wagner, W., Dorigo, W. (2019). Evolution of the ESA CCI Soil Moisture Climate Data Records and their underlying merging methodology. Earth System Science Data 11, 717-739, https://doi.org/10.5194/essd-11-717-2019


Example plots
-------------

.. figure:: /recipes/figures/autoassess_landsurface/Soilmoisture_Metrics.png
   :scale: 50 %
   :alt: Soilmoisture_Metrics.png

   Normalised metrics plot comparing a control and experiment simulation


Additional notes on usage
-------------------------
The ``landsurface_soilmoisture`` area metric is part of the ``esmvaltool/diag_scripts/autoassess`` diagnostics,
and, as any other ``autoassess`` metric, it uses the ``autoassess_area_base.py`` as general purpose
wrapper. This wrapper accepts a number of input arguments that are read through from the recipe.

This recipe is part of the larger group of Autoassess metrics ported to ESMValTool
from the native Autoassess package from the UK's Met Office. The ``diagnostics`` settings
are almost the same as for the other Autoassess metrics.

.. note::

   **Time gating for autoassess metrics.**

   To preserve the native Autoassess functionalities,
   data loading and selection on time is done somewhat
   differently for ESMValTool's autoassess metrics: the
   time selection is done in the preprocessor as per usual but
   a further time selection is performed as part of the diagnostic.
   For this purpose the user will specify a ``start:`` and ``end:``
   pair of arguments of ``scripts: autoassess_script`` (see below
   for example). These are formatted as ``YYYY/MM/DD``; this is
   necessary since the Autoassess metrics are computed from 1-Dec
   through 1-Dec rather than 1-Jan through 1-Jan. This is a temporary
   implementation to fully replicate the native Autoassess functionality
   and a minor user inconvenience since they need to set an extra set of
   ``start`` and ``end`` arguments in the diagnostic; this will be phased
   when all the native Autoassess metrics have been ported to ESMValTool
   review has completed.


An example of standard inputs as read by ``autoassess_area_base.py`` and passed
over to the diagnostic/metric is listed below.


.. code-block:: yaml

    scripts:
      autoassess_landsurf_soilmoisture: &autoassess_landsurf_soilmoisture_settings
        script: autoassess/autoassess_area_base.py
        title: "Autoassess Land-Surface Soilmoisture Diagnostic"
        area: land_surface_soilmoisture
        control_model: IPSL-CM5A-LR
        exp_model: inmcm4
        obs_models: []
        start: 1997/12/01
        end: 2002/12/01
        climfiles_root: '/gws/nopw/j04/esmeval/autoassess_specific_files/files'  # on JASMIN
.. _recipe_examples:

Example recipes
===============

Overview
--------

These are example recipes calling example diagnostic scripts.

The recipe examples/recipe_python.yml produces time series plots of global mean
temperature and for the temperature in Amsterdam.
It also produces a map of global temperature in January 2020.

The recipe examples/recipe_extract_shape.yml produces a map of the mean
temperature in the Elbe catchment over the years 2000 to 2002.
Some example shapefiles for use with this recipe are available
`here <https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/diag_scripts/shapeselect/testdata>`__,
make sure to download all files with the same name but different extensions.

The recipe examples/recipe_julia.yml produces a map plot with the mean temperature
over the year 1997 plus a number that is configurable from the recipe.

For detailed instructions on obtaining input data, please refer to
:ref:`inputdata`. However, in case you just quickly want to run through the
example, you can use the following links to obtain the data from ESGF:

  * `BCC-ESM1 <http://esgf3.dkrz.de/thredds/fileServer/cmip6/CMIP/BCC/BCC-ESM1/historical/r1i1p1f1/Amon/tas/gn/v20181214/tas_Amon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc>`_
  * `CanESM2 <http://esgf2.dkrz.de/thredds/fileServer/lta_dataroot/cmip5/output1/CCCma/CanESM2/historical/mon/atmos/Amon/r1i1p1/v20120718/tas/tas_Amon_CanESM2_historical_r1i1p1_185001-200512.nc>`_

Please refer to the terms of use for `CMIP5
<https://pcmdi.llnl.gov/mips/cmip5/terms-of-use.html>`_ and `CMIP6
<https://pcmdi.llnl.gov/CMIP6/TermsOfUse/TermsOfUse6-1.html>`_ data.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * examples/recipe_python.yml
    * examples/recipe_extract_shape.yml
    * examples/recipe_jula.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * examples/diagnostic.py: visualize results and store provenance information
    * examples/diagnostic.jl: visualize results and store provenance information

User settings in recipe
-----------------------

#. Script ``examples/diagnostic.py``

   *Required settings for script*

   * ``quickplot: plot_type``: which of the :py:mod:`iris.quickplot` functions to use.
     Arguments that are accepted by these functions can also be specified here, e.g. ``cmap``.
     Preprocessors need to be configured such that the resulting data matches the plot type, e.g. a timeseries or a map.

   *Optional settings for script*

   * ``write_netcdf``: ``true`` (default) or ``false``.
     This can be used to disable writing the results to netcdf files.

#. Script ``examples/diagnostic.jl``

   *Required settings for script*

   * ``parameter1``: example parameter, this number will be added to the mean (over time) value of the input data.

Variables
---------

* tas (atmos, monthly, longitude, latitude, time)

Example plots
-------------

.. _global_map:
.. figure::  /recipes/figures/examples/map.png
   :align:   center

   Air temperature in January 2000 (BCC-ESM1 CMIP6).

.. _timeseries:
.. figure::  /recipes/figures/examples/timeseries.png
   :align:   center

   Amsterdam air temperature (multimodel mean of CMIP5 CanESM2 and CMIP6 BCC-ESM1).

.. _elbe:
.. figure::  /recipes/figures/examples/elbe.png
   :align:   center

   Mean air temperature over the Elbe catchment during 2000-2002 according to CMIP5 CanESM2.
.. _recipes_pv_capacity_factor:

Capacity factor for solar photovoltaic (PV) systems
===================================================

Overview
--------

This diagnostic computes the photovoltaic (PV) capacity factor, 
a measure of the fraction of the 
maximum possible energy produced per PV grid cell. It uses the daily incoming 
surface solar radiation and the surface temperature with a method described
in `Bett and Thornton (2016)`_. The user can select temporal
range, season, and region of interest.


.. _`Bett and Thornton (2016)`: https://doi.org/10.1016/j.renene.2015.10.006


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_pv_capacity_factor.yml

Diagnostics are stored in diag_scripts/pv_capacityfactor/

    * pv_capacity_factor.R: prepares data and plots results.
    * PV_CF.R: calculates the daily capacity factor.


User settings
-------------

User setting files are stored in recipes/

#. recipe_capacity_factor.yml

   *Required settings for script*

   * season: String to include shortcut for season in plot title and name (e.g. "djf").
     It will be converted to upper case. This season should be the one set in the preprocessor,
     since it is only used as a string and does not affect the data in the diagnostic.
     In the default recipe this is solved through a node anchor.
   
   *Optional settings for script*
   
   * maxval_colorbar: Optional upper limit for the colorbar.

Variables
---------

* tas (atmos, daily, longitude, latitude, time)
* rsds (atmos, daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* ERA-Interim

References
----------

* Bett, P. E. and Thornton, H. E.: The climatological relationships between wind and solar energy supply in Britain, Renew. Energ., 87, 96–110, https://doi.org/10.1016/j.renene.2015.10.006, 2016.


Example plots
-------------

.. _fig_pv_capfactor1:
.. figure::  /recipes/figures/pv_capacity_factor/capacity_factor_IPSL-CM5A-MR_1980-2005_DJF.png
   :align:   center
   :width:   14cm

PV capacity factor calculated from IPSL-CM5-MR during the DJF season for 1980–2005... _recipe_eyring06jgr:

Diagnostics of stratospheric dynamics and chemistry
===================================================

Overview
--------

This recipe reproduces the figures of `Eyring et al. (2006)`_
The following plots are reproduced:

* Vertical profile climatological mean bias of climatological mean for selected seasons and latitudinal region.
* Vertical and latitudinal profile of climatological mean for selected seasons this figure and setting is valid for figure 5 (CH4) figure 6 (H2O) figure 11 (HCL) figure 13 (tro3).
* Total ozone anomalies at different latitudinal band and seasons.

.. _`Eyring et al. (2006)`: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2006JD007327

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_eyring06jgr.yml

Diagnostics are stored in esmvaltool/diag_scripts/eyring06jgr/

* eyring06jgr_fig01.ncl
* eyring06jgr_fig05a.ncl
* eyring06jgr_fig05b.ncl
* eyring06jgr_fig15.ncl

User settings in recipe
-----------------------
#. Preprocessor

   * ``regrid_interp_lev_zonal``: Regridding and interpolation reference_dataset levels used by eyring06jgr_fig01 and eyring06jgr_fig05
   * ``zonal`` : Regridding and zonal mean used by eyring06jgr_fig15


#. Script <eyring06jgr_fig01.ncl>

   *Required settings for script*

   * ``latmin``: array of float, min lat where variable is averaged, i.e. [60., 60., -90., -90. ]
   * ``latmax``: array of float,and max lat where variable is averaged, i.e. [90., 90., -60., -60. ]
   * ``season``: array of string., season when variable is averaged, i.e. ["DJF", "MAM", "JJA", "SON"]
   * ``XMin``: array of float, min limit X axis [-30., -30., -30., -30.]
   * ``XMax``: array of float, max limit X axis [20., 20., 20., 20.]
   * ``levmin``: array of float, min limit Y axis [1., 1., 1., 1.]
   * ``levmax``: array of float, max limit Y axis [350., 350., 350., 350.]


   *Optional settings for script*
   
   * ``start_year``: int,  year when start the climatology calculation [1980] (default max among the models start year).
   * ``end_year``:int, year when end  the climatology calculation [1999] (default min among the models end year).
   * ``multimean``: bool, calculate multi-model mean, (i.e. False/True) (default False).

   *Required settings for variables*
   
   * ``preprocessor``: regrid_interp_lev.
   * ``reference_dataset``: name of the reference model or observation for regridding and bias calculation (e.g. ERA-Interim").
   *  ``mip``:  Amon.



Variables
---------

*  ta (atmos, monthly mean, longitude latitude level time)



Example plots
-------------

.. _fig_eyring06jgr_01:
.. figure::  /recipes/figures/eyring06jgr/fig_diagn01.png
   :align:   center

   Climatological mean temperature biases for (top) 60–90N and (bottom) 60–90S for the (left) winter and (right) spring seasons. The climatological means for the CCMs and ERA-Interim data from 1980 to 1999 are included. Biases are calculated relative to ERA-Interim reanalyses. The grey area shows ERA-Interim plus and minus 1 standard deviation (s) about the climatological mean. The turquoise area shows plus and minus 1 standard deviation about the multi-model mean.
.. _recipes_thermodyn_diagtool:

Thermodynamics of the Climate System - The Diagnostic Tool TheDiaTo v1.0
========================================================================

Overview
--------

The tool allows to compute TOA, atmospheric and surface energy budgets, latent energy and water mass budgets,
meridional heat transports, the Lorenz Energy Cycle (LEC), the material entropy production with the direct
and indirect method.

The energy budgets are computed from monthly mean radiative and heat fluxes at the TOA and at the surface
(cfr. Wild et al., 2013). The meridional heat transports are obtained from the latitudinal integration
of the zonal mean energy budgets. When a land-sea mask is provided, results are also available for
land and oceans, separately.

The water mass budget is obtained from monthly mean latent heat fluxes (for evaporation), total and snowfall
precipitation (cfr. Liepert et al., 2012). The latent energy budget is obtained multiplying each component of
the water mass budget by the respective latent heat constant. When a land-sea mask is provided, results are
also available for land and oceans, separately.

The LEC is computed from 3D fields of daily mean velocity and temperature fields in the troposphere over
pressure levels. The analysis is carried on in spectral fields, converting lonlat grids in Fourier coefficients.
The components of the LEC are computed as in Ulbrich and Speth, 1991. In order to account for possible gaps
in pressure levels, the daily fields of 2D near-surface temperature and horizontal velocities are needed. These are
required to perform a vertical interpolation, substituting data in pressure levels where surface pressure is
lower than the respective level and fields are not stored as an output of the analysed model.

The material entropy production is computed by using the indirect or the direct method (or both). The former 
method relies on the convergence of radiative heat in the atmosphere (cfr. Lucarini et al., 2011; Pascale et al., 2011),
the latter on all viscous and non-viscous dissipative processes occurring in the atmosphere
(namely the sensible heat fluxes, the hydrological cycle with its components and the kinetic energy dissipation).

For a comprehensive report on the methods used and some descriptive results, please refer to Lembo et al., 2019.


In order to account for possible gaps in pressure levels, the daily
fields of 2D near-surface temperature and horizontal velocities.'

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_thermodyn_diagtool.yml

Diagnostics are stored in diag_scripts/thermodyn_diagtool/

    * thermodyn_diagnostics.py: the main script, handling input files, calling computation and plotting scricpts;

    * computations.py: a module containing all the main computations that are carried out by the program;

    * fluxogram.py: a module for the retrieval of the block diagrams displaying the reservoirs and conversion terms of the LEC

    * fourier_coefficients.py: a module for the computation of the Fourier coefficients from the lonlat input grid

    * lorenz_cycle.py: a module for the computation of the LEC components in Fourier coefficients

    * mkthe.py: a module for the computation of indirect variables obtained from the input fields, such as LCL height, boundary layer top height and temperature, potential temperature

    * plot_script.py: a module for the computation of maps, scatter plots, time series and meridional sections of some derived quantities for each model in the ensemble. The meridional heat and water mass transports are also computed here, as well as the peak magnitudes and locations;

    * provenance_meta.py: a module for collecting metadata and writing them to produced outputs;

User settings
-------------

Besides the datasets, to be set according to usual ESMValTool convention, the user can set the following optional variables in the recipe_Thermodyn_diagtool.yml:

   * wat: if set to 'true', computations are performed of the water mass and latent energy budgets and transports
   * lsm: if set to true, the computations of the energy budgets, meridional energy transports, water mass and latent energy budgets and transports are performed separately over land and oceans
   * lec: if set to 'true', computation of the LEC are performed
   * entr: if set to 'true', computations of the material entropy production are performed
   * met (1, 2 or 3): the computation of the material entropy production must be performed with the indirect method (1), the direct method (2), or both methods. If 2 or 3 options are chosen, the intensity of the LEC is needed for the entropy production related to the kinetic energy dissipation. If lec is set to 'false', a default value is provided.

   These options apply to all models provided for the multi-model ensemble computations


Variables
---------

Default variables needed for computation of energy budgets and transports:

* hfls    (atmos,  monthly mean, time latitude longitude)
* hfss    (atmos,  monthly mean, time latitude longitude)
* rlds    (atmos,  monthly mean, time latitude longitude)
* rlus    (atmos,  monthly mean, time latitude longitude)
* rlut    (atmos,  monthly mean, time latitude longitude)
* rsds    (atmos,  monthly mean, time latitude longitude)
* rsdt    (atmos,  monthly mean, time latitude longitude)
* rsus    (atmos,  monthly mean, time latitude longitude)
* rsut    (atmos,  monthly mean, time latitude longitude)

Additional variables needed for water mass and latent energy computation (optional, with 'wat' set to 'true'):

* pr      (atmos,  monthly mean, time latitude longitude)
* prsn    (atmos,  monthly mean, time latitude longitude)

Additional variable needed for LEC computations (optional, with 'lec' set to 'true'):

* ta      (atmos,  daily   mean, time plev latitude longitude)
* tas     (atmos,  daily   mean, time latitude longitude)
* ua      (atmos,  daily   mean, time plev latitude longitude)
* uas     (atmos,  daily   mean, time latitude longitude)
* va      (atmos,  daily   mean, time plev latitude longitude)
* vas     (atmos,  daily   mean, time latitude longitude)
* wap     (atmos,  daily   mean, time plev latitude longitude)

Additional variables needed for material entropy production computations with direct method (optional, with 'entr' set to 'true' and 'mep' to '2' or '3'):

* hus     (atmos,  monthly mean, time plev latitude longitude)
* pr      (atmos,  monthly mean, time latitude longitude)
* prsn    (atmos,  monthly mean, time latitude longitude)
* ps      (atmos,  monthly mean, time latitude longitude)
* ts      (atmos,  monthly mean, time latitude longitude)

Additional variables needed for material entropy production computations with indirect method (optional, with 'entr' set to 'true' and 'mep' to '1' or '3'):

* tas     (atmos,  daily   mean, time latitude longitude)
* uas     (atmos,  daily   mean, time latitude longitude)
* vas     (atmos,  daily   mean, time latitude longitude)

Depending on the user's options, variables listed above must be provided. All other variables shall be commented in the recipe file.


References
----------
* Lembo V, Lunkeit F, Lucarini V (2019) A new diagnostic tool for diagnosing water, energy and entropy budgets in climate models. Geophys Mod Dev Disc. doi:10.5194/gmd-12-3805-2019
* Liepert BG, Previdi M (2012) Inter-model variability and biases of the global water cycle in CMIP3 coupled climate models. Environ Res Lett 7:014006. doi: 10.1088/1748-9326/7/1/014006
* Lorenz EN (1955) Available Potential Energy and the Maintenance of the General Circulation. Tellus 7:157–167. doi: 10.1111/j.2153-3490.1955.tb01148.x
* Lucarini V, Fraedrich K, Ragone F (2010) New Results on the Thermodynamical Properties of the Climate System. J Atmo 68:. doi: 10.1175/2011JAS3713.1
* Lucarini V, Blender R, Herbert C, et al (2014) Reviews of Geophysics Mathematical and physical ideas for climate science. doi: 10.1002/2013RG000446
* Pascale S, Gregory JM, Ambaum M, Tailleux R (2011) Climate entropy budget of the HadCM3 atmosphere–ocean general circulation model and of FAMOUS, its low-resolution version. Clim Dyn 36:1189–1206. doi: 10.1007/s00382-009-0718-1
* Ulbrich U, Speth P (1991) The global energy cycle of stationary and transient atmospheric waves: Results from ECMWF analyses. Meteorol Atmos Phys 45:125–138. doi: 10.1007/BF01029650
* Wild M, Folini D, Schär C, et al (2013) The global energy balance from a surface perspective. Clim Dyn 40:3107–3134. doi: 10.1007/s00382-012-1569-8


Example plots
-------------

.. _fig_1:
.. figure:: /recipes/figures/thermodyn_diagtool/meridional_transp.png
   :align:   left
   :width:   14cm

.. _fig_2:
.. figure:: /recipes/figures/thermodyn_diagtool/CanESM2_wmb_transp.png
   :align:   right
   :width:   14cm
.. _recipes_ecs:

Equilibrium climate sensitivity
===============================

Overview
--------


Equilibrium climate sensitivity is defined as the change in global mean
temperature as a result of a doubling of the atmospheric CO\ :sub:`2`
concentration compared to pre-industrial times after the climate system has
reached a new equilibrium. This recipe uses a regression method based on
`Gregory et al. (2004)`_ to calculate it for several CMIP models.

.. _`Gregory et al. (2004)`: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003GL018747


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_ecs.yml


Diagnostics are stored in diag_scripts/

   * climate_metrics/ecs.py
   * climate_metrics/create_barplot.py
   * climate_metrics/create_scatterplot.py


User settings in recipe
-----------------------

* Preprocessor

   * ``area_statistics`` (*operation: mean*): Calculate global mean.

.. _ecs.py:

* Script climate_metrics/ecs.py

   * ``calculate_mmm``, *bool*, optional (default: ``True``): Calculate
     multi-model mean ECS.
   * ``complex_gregory_plot``, *bool*, optional (default: ``False``): Plot
     complex Gregory plot (also add response for first ``sep_year`` years and
     last 150 - ``sep_year`` years, default: ``sep_year=20``) if ``True``.
   * ``output_attributes``, *dict*, optional: Write additional attributes to
     netcdf files.
   * ``read_external_file``, *str*, optional: Read ECS and feedback parameters
     from external file. The path can be given relative to this diagnostic
     script or as absolute path.
   * ``savefig_kwargs``, *dict*, optional: Keyword arguments for
     :func:`matplotlib.pyplot.savefig`.
   * ``seaborn_settings``, *dict*, optional: Options for :func:`seaborn.set`
     (affects all plots).
   * ``sep_year``, *int*, optional (default: ``20``): Year to separate
     regressions of complex Gregory plot. Only effective if
     ``complex_gregory_plot`` is ``True``.
   * ``x_lim``, *list of float*, optional (default: ``[1.5, 6.0]``): Plot
     limits for X axis of Gregory regression plot (T).
   * ``y_lim``, *list of float*, optional (default: ``[0.5, 3.5]``): Plot
     limits for Y axis of Gregory regression plot (N).

.. _create_barplot.py:

* Script climate_metrics/create_barplot.py

   * ``add_mean``, *str*, optional: Add a bar representing the mean for each
     class.
   * ``label_attribute``, *str*, optional: Cube attribute which is used as
     label for different input files.
   * ``order``, *list of str*, optional: Specify the order of the different
     classes in the barplot by giving the ``label``, makes most sense when
     combined with ``label_attribute``.
   * ``patterns``, *list of str*, optional: Patterns to filter list of input
     data.
   * ``savefig_kwargs``, *dict*, optional: Keyword arguments for
     :func:`matplotlib.pyplot.savefig`.
   * ``seaborn_settings``, *dict*, optional: Options for :func:`seaborn.set`
     (affects all plots).
   * ``sort_ascending``, *bool*, optional (default: ``False``): Sort bars in
     ascending order.
   * ``sort_descending``, *bool*, optional (default: ``False``): Sort bars in
     descending order.
   * ``subplots_kwargs``, *dict*, optional: Keyword arguments for
     :func:`matplotlib.pyplot.subplots`.
   * ``value_labels``, *bool*, optional (default: ``False``): Label bars with
     value of that bar.
   * ``y_range``, *list of float*, optional: Range for the Y axis of the plot.

.. _create_scatterplot.py:

* Script climate_metrics/create_scatterplot.py

   * ``dataset_style``, *str*, optional: Name of the style file (located in
     :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
   * ``pattern``, *str*, optional: Pattern to filter list of input files.
   * ``seaborn_settings``, *dict*, optional: Options for :func:`seaborn.set`
     (affects all plots).
   * ``y_range``, *list of float*, optional: Range for the Y axis of the plot.


Variables
---------

* *rlut* (atmos, monthly, longitude, latitude, time)
* *rsdt* (atmos, monthly, longitude, latitude, time)
* *rsut* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Gregory, Jonathan M., et al. "A new method for diagnosing radiative forcing
  and climate sensitivity." Geophysical research letters 31.3 (2004).


Example plots
-------------

.. _fig_ecs_1:
.. figure:: /recipes/figures/ecs/CanESM2.png
   :align: center
   :width: 50%

   Scatterplot between TOA radiance and global mean surface temperature anomaly
   for 150 years of the abrupt 4x CO2 experiment including linear regression to
   calculate ECS for CanESM2 (CMIP5).
.. _recipes_arctic_ocean:

Recipe for evaluating Arctic Ocean
==================================

Overview
........

The Arctic ocean is one of the areas of the Earth where the effects of climate change are especially visible today. Two most prominent processes are Arctic amplification [e.g. Serreze and Barry, 2011] and decrease of the sea ice area and thickness. Both receive good coverage in the literature and are already well-studied. Much less attention is paid to the interior of the Arctic Ocean itself. In order to increase our confidence in projections of the Arctic climate future proper representation of the Arctic Ocean hydrography is necessary.

The main focus of this diagnostics is evaluation of ocean components of climate models in the Arctic Ocean, however most of the diagnostics are implemented in a way that can be easily expanded to other parts of the World Ocean. Most of the diagnostics aim at model comparison to climatological data (PHC3), so we target historical CMIP simulations. However scenario runs also can be analysed to have an impression of how Arcti Ocean hydrography will chnage in the future.

At present only the subset of CMIP models can be used in particular because our analysis is limited to z coordinate models.

Available recipes
.................
Recipe is stored in recipes/

* recipe_arctic_ocean.yml : contains all setting nessesary to run diagnostics and metrics.

Currenly the workflow do not allow to easily separate diagnostics from each other, since some of the diagnostics rely on the results of other diagnostics. The recipe currently do not use preprocessor options, so input files are CMORised monthly mean 3D ocean varibales on original grid.

The following plots will be produced by the recipe:

Hovmoeller diagrams
-------------------

The characteristics of vertical TS distribution can change with time, and consequently the vertical TS distribution is an important indicator of the behaviour of the coupled ocean-sea ice-atmosphere system in the North Atlantic and Arctic Oceans. One way to evaluate these changes is by using Hovmoller diagrams. Hovmoller diagrams for two main Arctic Ocean basins – Eurasian and Amerasian with T and S spatially averaged on a monthly basis for every vertical level are available. This diagnostic allows the temporal evolution of vertical ocean potential temperature distribution to be assessed.

Related settings in the recipe:

  .. code-block:: yaml

	# Define regions, as a list.
	# 'EB' - Eurasian Basin of the Arctic Ocean
	# 'AB' - Amerasian Basin of the Arctic Ocean
	# 'Barents_sea' - Barrents Sea
	# 'North_sea'   - North Sea
	hofm_regions: ["AB" ,  'EB']
	# Define variables to use, should also be in "variables"
	# entry of your diagnostic
	hofm_vars: ['thetao', 'so']
	# Maximum depth of Hovmoeller and vertical profiles
	hofm_depth: 1500
	# Define if Hovmoeller diagrams will be ploted.
	hofm_plot: True
	# Define colormap (as a list, same size as list with variables)
	# Only cmaps from matplotlib and cmocean are supported.
	# Additional cmap - 'custom_salinity1'.
	hofm_cmap: ['Spectral_r', 'custom_salinity1']
	# Data limits for plots,
	# List of the same size as the list of the variables
	# each entry is [vmin, vmax, number of levels, rounding limit]
	hofm_limits: [[-2, 2.3, 41, 1], [30.5, 35.1, 47, 2]]
	# Number of columns in the plot
	hofm_ncol: 3

.. _fig_hofm:
.. figure::  /recipes/figures/arctic_ocean/hofm.png
   :align:   center
   :width:   9cm

   Hovmoller diagram of monthly spatially averaged potential temperature in the Eurasian Basin of the Arctic Ocean for selected CMIP5 climate models (1970-2005).

Vertical profiles
-----------------

The vertical structure of temperature and salinity (T and S) in the ocean model is a key diagnostic that is used for ocean model evaluation. Realistic T and S distributions means that model properly represent dynamic and thermodynamic processes in the ocean. Different ocean basins have different hydrological regimes so it is important to perform analysis of vertical TS distribution for different basins separately. The basic diagnostic in this sense is mean vertical profiles of temperature and salinity over some basin averaged for relatively long period of time. In addition to individual vertical profiles for every model, we also show the mean over all participating models and similar profile from climatological data (PHC3).

Several settings for vertical profiles (region, variables, maximum depths) will be determined by the Hovmoeller diagram settings. The reason is that vertical profiles are calculated from Hovmoeller diagram data. Mean vertical profile is calculated by lineraly interpolating data on standard WOA/PHC depths.

Related settings in the recipe:

  .. code-block:: yaml

	# Define regions, as a list.
	# 'EB' - Eurasian Basin of the Arctic Ocean
	# 'AB' - Amerasian Basin of the Arctic Ocean
	# 'Barents_sea' - Barrents Sea
	# 'North_sea'   - North Sea
	hofm_regions: ["AB" ,  'EB']
	# Define variables to use, should also be in "variables" entry of your diagnostic
	hofm_vars: ['thetao', 'so']
	# Maximum depth of Hovmoeller and vertical profiles
	hofm_depth: 1500

.. _fig_vertical:
.. figure::  /recipes/figures/arctic_ocean/vertical.png
   :align:   center
   :width:   9cm

   Mean (1970-2005) vertical potential temperature distribution in the Eurasian basin for participating CMIP5 coupled ocean models, PHC3 climatology (dotted red line) and multi-model mean (dotted black line).

Spatial distribution maps of variables
--------------------------------------

The spatial distribution of basic oceanographic variables characterises the properties and spreading of ocean water masses. For the coupled models, capturing the spatial distribution of oceanographic variables is especially important in order to correctly represent the ocean-ice-atmosphere interface. We have implemented plots with spatial maps of temperature and salinity at original model levels.

Plots spatial distribution of variables at selected depths in North Polar projection on original model grid.
For plotting the model depths that are closest to provided `plot2d_depths` will be selected. Settings allow to define color maps and limits for each variable individually. Color maps should be ehter part of standard matplotlib set or one of the cmocean color maps. Additional colormap `custom_salinity1` is provided.

Related settings in the recipe:

  .. code-block:: yaml

	# Depths for spatial distribution maps
	plot2d_depths: [10, 100]
	# Variables to plot spatial distribution maps
	plot2d_vars: ['thetao', 'so']
	# Define colormap (as a list, same size as list with variables)
	# Only cmaps from matplotlib and cmocean are supported.
	# Additional cmap - 'custom_salinity1'.
	plot2d_cmap: ['Spectral_r', 'custom_salinity1']
	# Data limits for plots,
	# List of the same size as the list of the variables
	# each entry is [vmin, vmax, number of levels, rounding limit]
	plot2d_limits: [[-2, 4, 20, 1], [30.5, 35.1, 47, 2]]
	# number of columns for plots
	plot2d_ncol: 3

.. _fig_spatial:
.. figure::  /recipes/figures/arctic_ocean/spatial.png
   :align:   center
   :width:   9cm

   Mean (1970-2005) salinity distribution at 100 meters.

Spatial distribution maps of biases
-----------------------------------

For temperature and salinity, we have implemented spatial maps of model biases from the observed climatology. For the model biases, values from the original model levels are linearly interpolated to the climatology and then spatially interpolated from the model grid to the regular PHC (climatology) grid. Resulting fields show model performance in simulating spatial distribution of temperature and salinity.

Related settings in the recipe:

  .. code-block:: yaml

	plot2d_bias_depths: [10, 100]
	# Variables to plot spatial distribution of the bias for.
	plot2d_bias_vars: ['thetao', 'so']
	# Color map names for every variable
	plot2d_bias_cmap: ['balance', 'balance']
	# Data limits for plots,
	# List of the same size as the list of the variables
	# each entry is [vmin, vmax, number of levels, rounding limit]
	plot2d_bias_limits: [[-3, 3, 20, 1], [-2, 2, 47, 2]]
	# number of columns in the bias plots
	plot2d_bias_ncol: 3

.. _fig_bias:
.. figure::  /recipes/figures/arctic_ocean/bias.png
   :align:   center
   :width:   9cm

   Mean (1970-2005) salinity bias at 100m relative to PHC3 climatology

Transects
---------
Vertical transects through arbitrary sections are important for analysis of vertical distribution of ocean water properties and especially useful when exchange between different ocean basins is evaluated. We have implemented diagnostics that allow for the definition of an arbitrary ocean section by providing set of points on the ocean surface. For each point, a vertical profile on the original model levels is interpolated. All profiles are then connected to form a transect. The great-circle distance between the points is calculated and used as along-track distance.

One of the main use cases is to create vertical sections across ocean passages, for example Fram Strait.

Plots transect maps for pre-defined set of transects (defined in `regions.py`, see below). The `transect_depth` defines maximum depth of the transect. Transects are calculated from data averaged over the whole time period.

Related settings in the recipe:

  .. code-block:: yaml

	# Select regions (transects) to plot
	# Available options are:
	# AWpath - transect along the path of the Atlantic Water
	# Fram - Fram strait
	transects_regions: ["AWpath", "Fram"]
	# Variables to plot on transects
	transects_vars: ['thetao', 'so']
	# Color maps for every variable
	transects_cmap: ['Spectral_r', 'custom_salinity1']
	# Data limits for plots,
	# List of the same size as the list of the variables
	# each entry is [vmin, vmax, number of levels, rounding limit]
	transects_limits: [[-2, 4, 20, 1], [30.5, 35.1, 47, 2]]
	# Maximum depth to plot the data
	transects_depth: 1500
	# number of columns
	transects_ncol: 3

.. _fig_transect:
.. figure::  /recipes/figures/arctic_ocean/transect.png
   :align:   center
   :width:   9cm

   Mean (1970-2005) potential temperature across the Fram strait.

Atlantic Water core depth and temperature
-----------------------------------------

Atlantic water is a key water mass of the Arctic Ocean and its proper representation is one of the main challenges in Arctic Ocean modelling. We have created two metrics by which models can be easily compared in terms of Atlantic water simulation. The temperature of the Atlantic Water core is calculated for every model as the maximum potential temperature between 200 and 1000 meters depth in the Eurasian Basin. The depth of the Atlantic Water core is calculated as the model level depth where the maximum temperature is found in Eurasian Basin (Atlantic water core temperature).

The AW core depth and temperature will be calculated from data generated for Hovmoeller diagrams for `EB` region, so it should be selected in the Hovmoeller diagrams settings as one of the `hofm_regions`.

In order to evaluate the spatial distribution of Atlantic water in different climate models we also provide diagnostics with maps of the spatial temperature distribution at model’s Atlantic Water depth.

.. _fig_aw_temp:
.. figure::  /recipes/figures/arctic_ocean/aw_temp.png
   :align:   center
   :width:   9cm

   Mean (1970-2005) Atlantic Water core temperature. PHC33 is an observed climatology.

TS-diagrams
-----------

T-S diagrams combine temperature and salinity, which allows the analysis of water masses and their potential for mixing. The lines of constant density for specific ranges of temperature and salinity are shown on the background of the T-S diagram. The dots on the diagram are individual grid points from specified region at all model levels within user specified depth range.

Related settings in the recipe:

  .. code-block:: yaml

	tsdiag_regions: ["AB" ,  'EB']
	# Maximum depth to consider data for TS diagrams
	tsdiag_depth: 1500
	# Number of columns
	tsdiag_ncol: 3

.. _fig_ts:
.. figure::  /recipes/figures/arctic_ocean/ts.png
   :align:   center
   :width:   9cm

   Mean (1970-2005) T-S diagrams for Eurasian Basin of the Arctic Ocean.

Available diagnostics
.....................

The following python modules are included in the diagnostics package:

* arctic_ocean.py : Reads settings from the recipe and call functions to do analysis and plots.
* getdata.py : Deals with data preparation.
* interpolation.py	: Include horizontal and vertical interpolation functions specific for ocean models.
* plotting.py : Ocean specific plotting functions
* regions.py : Contains code to select specific regions, and definition of the regions themselves.
* utils.py : Helpful utilites.

Diagnostics are stored in diag_scripts/arctic_ocean/


Variables
---------

* thetao (ocean, monthly, longitude, latitude, time)
* so (ocean, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* PHC3 climatology

References
----------

* Ilıcak, M. et al., An assessment of the Arctic Ocean in a suite of interannual CORE-II simulations. Part III: Hydrography and fluxes, Ocean Modelling, Volume 100, April 2016, Pages 141-161, ISSN 1463-5003, doi.org/10.1016/j.ocemod.2016.02.004

* Steele, M., Morley, R., & Ermold, W. (2001). PHC: A global ocean hydrography with a high-quality Arctic Ocean. Journal of Climate, 14(9), 2079-2087.

* Wang, Q., et al., An assessment of the Arctic Ocean in a suite of interannual CORE-II simulations. Part I: Sea ice and solid freshwater, Ocean Modelling, Volume 99, March 2016, Pages 110-132, ISSN 1463-5003, doi.org/10.1016/j.ocemod.2015.12.008

* Wang, Q., Ilicak, M., Gerdes, R., Drange, H., Aksenov, Y., Bailey, D. A., ... & Cassou, C. (2016). An assessment of the Arctic Ocean in a suite of interannual CORE-II simulations. Part II: Liquid freshwater. Ocean Modelling, 99, 86-109, doi.org/10.1016/j.ocemod.2015.12.009
.. _recipes_tcr:

Transient Climate Response
==========================

Overview
--------


The transient climate response (TCR) is defined as the global and annual mean
surface air temperature anomaly in the *1pctCO2* scenario (1% CO\ :sub:`2`
increase per year) for a 20 year period centered at the time of CO\ :sub:`2`
doubling, i.e. using the years 61 to 80 after the start of the simulation. We
calculate the temperature anomaly by subtracting a linear fit of the
*piControl* run for all 140 years of the *1pctCO2* experiment prior to the TCR
calculation (see `Gregory and Forster, 2008`_).

.. _`Gregory and Forster, 2008`: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JD010405


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_tcr.yml


Diagnostics are stored in diag_scripts/

   * climate_metrics/tcr.py
   * climate_metrics/create_barplot.py
   * climate_metrics/create_scatterplot.py


User settings in recipe
-----------------------

* Preprocessor

   * ``area_statistics`` (*operation: mean*): Calculate global mean.

.. _tcr.py:

* Script climate_metrics/tcr.py

   * ``calculate_mmm``, *bool*, optional (default: ``True``): Calculate
     multi-model mean TCR.
   * ``plot``, *bool*, optional (default: ``True``): Plot temperature vs. time.
   * ``read_external_file``, *str*, optional: Read TCR from external file. The
     path can be given relative to this diagnostic script or as absolute path.
   * ``savefig_kwargs``, *dict*, optional: Keyword arguments for
     :func:`matplotlib.pyplot.savefig`.
   * ``seaborn_settings``, *dict*, optional: Options for :func:`seaborn.set`
     (affects all plots).

* Script climate_metrics/create_barplot.py

   See :ref:`here<create_barplot.py>`.

* Script climate_metrics/create_scatterplot.py

   See :ref:`here<create_scatterplot.py>`.


Variables
---------

* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Gregory, J. M., and P. M. Forster. "Transient climate response estimated from
  radiative forcing and observed temperature change." Journal of Geophysical
  Research: Atmospheres 113.D23 (2008).


Example plots
-------------

.. _fig_tcr_1:
.. figure:: /recipes/figures/tcr/CanESM2.png
   :align: center
   :width: 50%

   Time series of the global mean surface air temperature anomaly (relative to
   the linear fit of the pre-industrial control run) of CanESM2 (CMIP5) for the
   1% CO\ :sub:`2` increase per year experiment. The horizontal dashed line
   indicates the transient climate response (TCR) defined as the 20 year
   average temperature anomaly centered at the time of CO\ :sub:`2` doubling
   (vertical dashed lines).
.. _recipes_sea_surface_salinity:

Sea Surface Salinity Evaluation
===============================

Overview
--------

This recipe compares the regional means of sea surface salinity with a
reference dataset (ESACCI-SEA-SURFACE-SALINITY v1 or v2 by default).
To do this, the recipe generates plots for the timeseries of each region and
a radar plot showing the correlation between dataset and reference timeseries for
each region during the time they both exist.

Preprocessor requirements:
--------------------------

The recipe is created in a way that should make possible (although is not
tested) to use it for other variables and datasets, even for more than one at
a time. The diagnostic only expects variables with dimensions `time` and `depth_id`
and it does not assume any other constraint.

It is therefore mandatory to keep the `extract_shape` preprocessor for more than
one region and any form of region operation (`mean`, `max`, `min` ...) to collapse
the `latitude` and `longitude` coordinates. In case you want to try with variables
that have extra dimensions (i.e. `depth`) you must add an extra preprocessor
call to collapse them (i.e. `depth_integration`)

The recipe can be used with any shapefile. As it is, it uses the IHO Sea Areas
(version 3) downloaded from https://marineregions.org/downloads.php, but any
shapefile containing marine regions can be used.

Any number of regions can be choosed also, even though plots may look odd if
too few or too many are selected.

Regions available on IHO Sea Areas file:
----------------------------------------

- Adriatic Sea
- Aegean Sea
- Alboran Sea
- Andaman or Burma Sea
- Arabian Sea
- Arafura Sea
- Arctic Ocean
- Baffin Bay
- Balearic (Iberian Sea)
- Bali Sea
- Baltic Sea
- Banda Sea
- Barentsz Sea
- Bass Strait
- Bay of Bengal
- Bay of Biscay
- Bay of Fundy
- Beaufort Sea
- Bering Sea
- Bismarck Sea
- Black Sea
- Bristol Channel
- Caribbean Sea
- Celebes Sea
- Celtic Sea
- Ceram Sea
- Chukchi Sea
- Coral Sea
- Davis Strait
- East Siberian Sea
- Eastern China Sea
- English Channel
- Flores Sea
- Great Australian Bight
- Greenland Sea
- Gulf of Aden
- Gulf of Alaska
- Gulf of Aqaba
- Gulf of Boni
- Gulf of Bothnia
- Gulf of California
- Gulf of Finland
- Gulf of Guinea
- Gulf of Mexico
- Gulf of Oman
- Gulf of Riga
- Gulf of St. Lawrence
- Gulf of Suez
- Gulf of Thailand
- Gulf of Tomini
- Halmahera Sea
- Hudson Bay
- Hudson Strait
- Indian Ocean
- Inner Seas off the West Coast of Scotland
- Ionian Sea
- Irish Sea and St. George's Channel
- Japan Sea
- Java Sea
- Kara Sea
- Kattegat
- Labrador Sea
- Laccadive Sea
- Laptev Sea
- Ligurian Sea
- Lincoln Sea
- Makassar Strait
- Malacca Strait
- Mediterranean Sea - Eastern Basin
- Mediterranean Sea - Western Basin
- Molukka Sea
- Mozambique Channel
- North Atlantic Ocean
- North Pacific Ocean
- North Sea
- Norwegian Sea
- Persian Gulf
- Philippine Sea
- Red Sea
- Rio de La Plata
- Savu Sea
- Sea of Azov
- Sea of Marmara
- Sea of Okhotsk
- Seto Naikai or Inland Sea
- Singapore Strait
- Skagerrak
- Solomon Sea
- South Atlantic Ocean
- South China Sea
- South Pacific Ocean
- Southern Ocean
- Strait of Gibraltar
- Sulu Sea
- Tasman Sea
- The Coastal Waters of Southeast Alaska and British Columbia
- The Northwestern Passages
- Timor Sea
- Tyrrhenian Sea
- White Sea
- Yellow Sea


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_sea_surface_salinity.yml

Diagnostics are stored in diag_scripts/sea_surface_salinity/

    * compare_salinity.py: plot timeseries for each region and generate radar
      plot.


User settings in recipe
-----------------------

#. compare_salinity.py

   *Required settings for script*

   none

   *Optional settings for script*

   none

   *Required settings for variables*

   * ref_model: name of reference data set

   *Optional settings for variables*

   none


Variables
---------

* sos (ocean, monthly, time latitude longitude)


Observations and reformat scripts
---------------------------------

* ESACCI-SEA-SURFACE-SALINITY (sos)


References
----------

* Diagnostic: please contact authors

* ESACCI-SEA-SURFACE-SALINITY dataset: Boutin, J., J.-L. Vergely, J. Koehler,
  F. Rouffi, N. Reul: ESA Sea Surface Salinity Climate Change Initiative
  (Sea_Surface_Salinity_cci): Version 1.8 data collection. Centre for
  Environmental Data Analysis, 25 November 2019. doi:
  10.5285/9ef0ebf847564c2eabe62cac4899ec41.
  http://dx.doi.org/10.5285/9ef0ebf847564c2eabe62cac4899ec41


Example plots
-------------

.. figure:: /recipes/figures/sea_surface_salinity/radar.png
   :align: center

   Radar plot showing correlation of average sea surface salinity for multiple
   regions with the observations
.. _recipes_hydrology:

Hydrological models - data pre-processing
=========================================

Overview
--------

We provide a collection of scripts that pre-processes environmental data for use in several hydrological models:

PCR-GLOBWB
**********
PCR-GLOBWB (PCRaster Global Water Balance) is a large-scale hydrological model intended for global to regional studies and developed at the Department of Physical Geography, Utrecht University (Netherlands). The recipe pre-processes ERA-Interim reanalyses data for use in the PCR-GLOBWB.

MARRMoT
**********
MARRMoT (Modular Assessment of Rainfall-Runoff Models Toolbox) is a rainfall-runoff model comparison framework that allows objective comparison between different conceptual hydrological model structures https://github.com/wknoben/MARRMoT. The recipe pre-processes ERA-Interim and ERA5 reanalyses data for use in the MARRMoT.

MARRMoT requires potential evapotranspiration (evspsblpot). The variable evspsblpot is not available in ERA-Interim. Thus, we use the debruin function (De Bruin et al. 2016) to obtain evspsblpot using both ERA-Interim and ERA5. This function needs the variables tas, psl, rsds, and rsdt as input.

wflow_sbm and wflow_topoflex
****************************
Forcing data for the `wflow_sbm <https://wflow.readthedocs.io/en/latest/wflow_sbm.html>`_
and `wflow_topoflex <https://wflow.readthedocs.io/en/latest/wflow_topoflex.html>`_
hydrological models can be prepared using recipe_wflow.yml.
If PET is not available from the source data (e.g. ERA-Interim), then it can be derived from psl, rsds and rsdt using De Bruin's 2016 formula (De Bruin et al. 2016). For daily ERA5 data, the time points of these variables are shifted 30 minutes with respect to one another. This is because in ERA5, accumulated variables are recorded over the past hour, and in the process of cmorization, we shift the time coordinates to the middle of the interval over which is accumulated. However, computing daily statistics then averages the times, which results in 12:00 UTC for accumulated variables and 11:30 UTC for instantaneous variables. Therefore, in this diagnostic, the time coordinates of the daily instantaneous variables are shifted 30 minutes forward in time.

LISFLOOD
********
`LISFLOOD <https://ec-jrc.github.io/lisflood-model/>`_ is a spatially distributed water resources model, developed by the Joint Research Centre (JRC) of the European Commission since 1997. We provide a recipe to produce meteorological forcing data for the Python 3 version of LISFLOOD.

LISFLOOD has a separate preprocessor LISVAP that derives some additional variables. We don't replace LISVAP. Rather, we provide input files that can readily be passed to LISVAP and then to LISFLOOD.


HYPE
****

The hydrological catchment model HYPE simulates water flow and substances on their way from precipitation through soil, river and lakes to the river outlet.
HYPE is developed at the Swedish Meteorological and Hydrological Institute. The recipe pre-processes ERA-Interim and ERA5 data for use in HYPE.

GlobWat
*******
GlobWat is a soil water balance model that has been provided by the Food and Agriculture Organization (FAO) to assess water use in irrigated agriculture (http://www.fao.org/nr/water/aquamaps). The recipe pre-processes ERA-Interim and ERA5 reanalyses data for use in the GlobWat model. GlobWat requires potential evapotranspiration (evspsblpot) as input. The variable evspsblpot is not available in ERA-Interim. Thus, we use debruin function (De Bruin et al. 2016) or the langbein method (Langbein et al. 1949) to obtain evspsblpot using both ERA-Interim and ERA5. The Langbein function needs a variable tas and the debruin function besides that needs the variables psl, rsds, and rsdt as input. In order to calculate monthly/daily pet with Langbein method we assumed that tas is constant over time and the average value is equal to the annual average.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/hydrology

    * recipe_pcrglobwb.yml
    * recipe_marrmot.yml
    * recipe_wflow.yml
    * recipe_lisflood.yml
    * recipe_hype.yml
    * recipe_globwat.yml

Diagnostics are stored in esmvaltool/diag_scripts/hydrology

    * pcrglobwb.py
    * marrmot.py
    * wflow.py
    * lisflood.py
    * hype.py
    * globwat.py 


User settings in recipe
-----------------------

All hydrological recipes require a shapefile as an input to produce forcing data. This shapefile determines the shape of the basin for which the data will be cut out and processed. All recipes are tested with `the shapefiles <https://github.com/eWaterCycle/recipes_auxiliary_datasets/tree/main/>`_  that are used for the eWaterCycle project. In principle any shapefile can be used, for example, the freely available basin shapefiles from the `HydroSHEDS project <https://www.hydrosheds.org/>`_. 

#. recipe_pcrglobwb.yml

   *Required preprocessor settings:*

   * start_year: 1979
   * end_year: 1979

#. recipe_marrmot.yml

   There is one diagnostic ``diagnostic_daily`` for using daily data.

   *Required preprocessor settings:*

      The settings below should not be changed.

      *extract_shape:*

         * shapefile: Meuse.shp (MARRMoT is a hydrological Lumped model that needs catchment-aggregated forcing data. The catchment is provided as a shapefile, the path can be relative to ``auxiliary_data_dir`` as defined in config-user.yml.).
         * method: contains
         * crop: true

   *Required diagnostic script settings:*

      * basin: Name of the catchment

#. recipe_wflow.yml

   *Optional preprocessor settings:*

      * extract_region: the region specified here should match the catchment

   *Required diagnostic script settings:*

	    * basin: name of the catchment
	    * dem_file: netcdf file containing a digital elevation model with
	      elevation in meters and coordinates latitude and longitude.
	    * regrid: the regridding scheme for regridding to the digital elevation model. Choose ``area_weighted`` (slow) or ``linear``.

#. recipe_lisflood.yml

   *Required preprocessor settings:*

      * extract_region: A region bounding box slightly larger than the shapefile. This is run prior to regridding, to save memory.
      * extract_shape:*

         * shapefile: A shapefile that specifies the extents of the catchment.

         These settings should not be changed

         * method: contains
         * crop: true

      * regrid:*

         * target_grid: Grid of LISFLOOD input files

         These settings should not be changed

         * lon_offset: true
         * lat_offset: true
         * scheme: linear

   There is one diagnostic ``diagnostic_daily`` for using daily data.

   *Required diagnostic script settings:*

      * catchment: Name of the catchment, used in output filenames

#. recipe_hype.yml

   *Required preprocessor settings:*

   * start_year: 1979
   * end_year: 1979
   * shapefile: Meuse_HYPE.shp (expects shapefile with subcatchments)

   These settings should not be changed

   * method: contains
   * decomposed: true

#. recipe_globwat.yml

   *Required preprocessor settings:*

   * start_year: 2004
   * end_year: 2004
   * target_grid_file: grid of globwat input files. A target file has been generated from one of the GlobWat models sample files (prc01wb.asc) for regridding ERA5 and ERA-Interim datasets. The ASCII file can be found at: https://storage.googleapis.com/fao-maps-catalog-data/geonetwork/aquamaps/GlobWat-InputP1_prec.zip. You can use the GDAL translator to convert the file from ASCII format to NetCDF format by entering the following command into the terminal: gdal_translate -of netCDF prc01wb.asc globwat_target_grid.nc

   *Optional preprocessor settings:*

   * area_selection: A region bounding box to extract the data for a specific region. The area selection preprocessor can be used by users to process the data for their desired region. The data will be processed at the global scale if the preprocessor in the recipe is commented.
   * regrid_scheme: The area-weighted regridding scheme is used as a default regridding scheme to ensure that the total volume of water is consistent before and after regridding.
   * langbein_pet: Can be set to True to use langbein function for calculating evspsblpot (default is de bruin method)


Variables
---------

#. recipe_pcrglobwb.yml

   * tas (atmos, daily, longitude, latitude, time)
   * pr (atmos, daily, longitude, latitude, time)

#. recipe_marrmot.yml

   * pr (atmos, daily or hourly mean, longitude, latitude, time)
   * psl (atmos, daily or hourly mean, longitude, latitude, time)
   * rsds (atmos, daily or hourly mean, longitude, latitude, time)
   * rsdt (atmos, daily or hourly mean, longitude, latitude, time)
   * tas (atmos, daily or hourly mean, longitude, latitude, time)

#. recipe_wflow.yml

   * orog (fx, longitude, latitude)
   * pr (atmos, daily or hourly mean, longitude, latitude, time)
   * tas (atmos, daily or hourly mean, longitude, latitude, time)

   Either potential evapotranspiration can be provided:

   * evspsblpot(atmos, daily or hourly mean, longitude, latitude, time)

   or it can be derived from tas, psl, rsds, and rsdt using the De Bruin formula, in that case the following variables need to be provided:

   * psl (atmos, daily or hourly mean, longitude, latitude, time)
   * rsds (atmos, daily or hourly mean, longitude, latitude, time)
   * rsdt (atmos, daily or hourly mean, longitude, latitude, time)

#. recipe_lisflood.yml

   * pr (atmos, daily, longitude, latitude, time)
   * tas (atmos, daily, longitude, latitude, time)
   * tasmax (atmos, daily, longitude, latitude, time)
   * tasmin (atmos, daily, longitude, latitude, time)
   * tdps (atmos, daily, longitude, latitude, time)
   * uas (atmos, daily, longitude, latitude, time)
   * vas (atmos, daily, longitude, latitude, time)
   * rsds (atmos, daily, longitude, latitude, time)

#. recipe_hype.yml

   * tas (atmos, daily or hourly, longitude, latitude, time)
   * tasmin (atmos, daily or hourly, longitude, latitude, time)
   * tasmax (atmos, daily or hourly, longitude, latitude, time)
   * pr (atmos, daily or hourly, longitude, latitude, time)

#. recipe_globwat.yml

   * pr (atmos, daily or monthly, longitude, latitude, time)
   * tas (atmos, daily or monthly, longitude, latitude, time)
   * psl (atmos, daily or monthly, longitude, latitude, time)
   * rsds (atmos, daily or monthly, longitude, latitude, time)
   * rsdt (atmos, daily or monthly , longitude, latitude, time)

Observations and reformat scripts
---------------------------------
*Note: see headers of cmorization scripts (in esmvaltool/cmorizers/obs) for download instructions.*

*  ERA-Interim (esmvaltool/cmorizers/obs/cmorize_obs_era_interim.py)
*  ERA5 (esmvaltool/cmorizers/obs/cmorize_obs_era5.py)

Output
---------

#. recipe_pcrglobwb.yml

#. recipe_marrmot.yml

    The forcing data, the start and end times of the forcing data, the latitude and longitude of the catchment are saved in a .mat file as a data structure readable by MATLAB or Octave.

#. recipe_wflow.yml

	The forcing data, stored in a single NetCDF file.

#. recipe_lisflood.yml

   The forcing data, stored in separate files per variable.

#. recipe_globwat.yml

   The forcing data, stored in separate files per timestep and variable.

References
----------

* Sutanudjaja, E. H., van Beek, R., Wanders, N., Wada, Y., Bosmans, J. H. C., Drost, N., van der Ent, R. J., de Graaf, I. E. M., Hoch, J. M., de Jong, K., Karssenberg, D., López López, P., Peßenteiner, S., Schmitz, O., Straatsma, M. W., Vannametee, E., Wisser, D., and Bierkens, M. F. P.: PCR-GLOBWB 2: a 5 arcmin global hydrological and water resources model, Geosci. Model Dev., 11, 2429-2453, https://doi.org/10.5194/gmd-11-2429-2018, 2018.
* De Bruin, H. A. R., Trigo, I. F., Bosveld, F. C., Meirink, J. F.: A Thermodynamically Based Model for Actual Evapotranspiration of an Extensive Grass Field Close to FAO Reference, Suitable for Remote Sensing Application, American Meteorological Society, 17, 1373-1382, DOI: 10.1175/JHM-D-15-0006.1, 2016.
* Arheimer, B., Lindström, G., Pers, C., Rosberg, J. och J. Strömqvist, 2008. Development and test of a new Swedish water quality model for small-scale and large-scale applications. XXV Nordic Hydrological Conference, Reykjavik, August 11-13, 2008. NHP Report No. 50, pp. 483-492.
* Lindström, G., Pers, C.P., Rosberg, R., Strömqvist, J., Arheimer, B. 2010. Development and test of the HYPE (Hydrological Predictions for the Environment) model – A water quality model for different spatial scales. Hydrology Research 41.3-4:295-319.
* van der Knijff, J. M., Younis, J. and de Roo, A. P. J.: LISFLOOD: A GIS-based distributed model for river basin scale water balance and flood simulation, Int. J. Geogr. Inf. Sci., 24(2), 189–212, 2010.
* Hoogeveen, J., Faurès, J. M., Peiser, L., Burke, J., de Giesen, N. V.: GlobWat--a global water balance model to assess water use in irrigated agriculture, Hydrology & Earth System Sciences Discussions, 2015 Jan 1;12(1), Doi:10.5194/hess-19-3829-2015.
* Langbein, W.B., 1949. Annual runoff in the United States. US Geol. Surv.(https://pubs.usgs.gov/circ/1949/0052/report.pdf)
.. _recipe_ecs_scatter:

Emergent constraints for equilibrium climate sensitivity
========================================================

Overview
--------

Calculates equilibrium climate sensitivity (ECS) versus

1) S index, D index and lower tropospheric mixing index (LTMI); similar to fig. 5 from Sherwood et al. (2014)
2) southern ITCZ index and tropical mid-tropospheric humidity asymmetry index; similar to fig. 2 and 4 from Tian (2015)
3) covariance of shortwave cloud reflection (Brient and Schneider, 2016)
4) climatological Hadley cell extent (Lipat et al., 2017)
5) temperature variability metric; similar to fig. 2 from Cox et al. (2018)
6) total cloud fraction difference between tropics and mid-latitudes; similar to fig. 3 from Volodin (2008)
7) response of marine boundary layer cloud (MBLC) fraction changes to sea surface temperature (SST); similar to fig. 3 of Zhai et al. (2015)
8) Cloud shallowness index (Brient et al., 2016)
9) Error in vertically-resolved tropospheric zonal average relative humidity (Su et al., 2014)

The results are displayed as scatterplots.

.. note:: The recipe ``recipe_ecs_scatter.yml`` requires pre-calulation of the
   equilibrium climate sensitivites (ECS) for all models. The ECS values are
   calculated with recipe_ecs.yml. The netcdf file containing the ECS values
   (path and filename) is specified by diag_script_info@ecs_file.
   Alternatively, the netcdf file containing the ECS values can be generated
   with the cdl-script
   $diag_scripts/emergent_constraints/ecs_cmip.cdl (recommended method):

   1) save script given at the end of this recipe as ecs_cmip.cdl
   2) run command: ncgen -o ecs_cmip.nc ecs_cmip.cdl
   3) copy ecs_cmip.nc to directory given by diag_script_info@ecs_file
      (e.g. $diag_scripts/emergent_constraints/ecs_cmip.nc)


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_ecs_scatter.yml
    * recipe_ecs_constraints.yml

Diagnostics are stored in diag_scripts

    * emergent_constraints/ecs_scatter.ncl: calculate emergent constraints for ECS
    * emergent_constraints/ecs_scatter.py: calculate further emergent constraints for ECS
    * emergent_constraints/single_constraint.py: create scatterplots for emergent constraints
    * climate_metrics/psi.py: calculate temperature variabililty metric (Cox et al., 2018)


User settings in recipe
-----------------------

.. _ecs_scatter.ncl:

* Script emergent_constraints/ecs_scatter.ncl

   *Required settings (scripts)*

   * diag: emergent constraint to calculate ("itczidx", "humidx", "ltmi",
     "covrefl", "shhc", "sherwood_d", "sherwood_s")
   * ecs_file: path and filename of netCDF containing precalculated
     ECS values (see note above)

   *Optional settings (scripts)*

   * calcmm: calculate multi-model mean (True, False)
   * legend_outside: plot legend outside of scatterplots (True, False)
   * output_diag_only: Only write netcdf files for X axis (True) or write all
     plots (False)
   * output_models_only: Only write models (no reference datasets) to netcdf
     files (True, False)
   * output_attributes: Additonal attributes for all output netcdf files
   * predef_minmax: use predefined internal min/max values for axes
     (True, False)
   * styleset: "CMIP5" (if not set, diagnostic will create a color table
     and symbols for plotting)
   * suffix: string to add to output filenames (e.g."cmip3")

   *Required settings (variables)*

   * reference_dataset: name of reference data set

   *Optional settings (variables)*

   none

   *Color tables*

   none


* Script emergent_constraints/ecs_scatter.py

   See
   :ref:`here<api.esmvaltool.diag_scripts.emergent_constraints.ecs_scatter>`.


* Script emergent_constraints/single_constraint.py

   See
   :ref:`here<api.esmvaltool.diag_scripts.emergent_constraints.single_constraint>`.


* Script climate_metrics/psi.py

   See :ref:`here<psi.py>`.


Variables
---------

* cl (atmos, monthly mean, longitude latitude level time)
* clt (atmos, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* hur (atmos, monthly mean, longitude latitude level time)
* hus (atmos, monthly mean, longitude latitude level time)
* rsdt (atmos, monthly mean, longitude latitude time)
* rsut (atmos, monthly mean, longitude latitude time)
* rsutcs (atmos, monthly mean, longitude latitude time)
* rtnt or rtmt (atmos, monthly mean, longitude latitude time)
* ta (atmos, monthly mean, longitude latitude level time)
* tas (atmos, monthly mean, longitude latitude time)
* tasa (atmos, monthly mean, longitude latitude time)
* tos (atmos, monthly mean, longitude latitude time)
* ts (atmos, monthly mean, longitude latitude time)
* va (atmos, monthly mean, longitude latitude level time)
* wap (atmos, monthly mean, longitude latitude level time)
* zg (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

.. note:: (1) Obs4mips data can be used directly without any preprocessing.
          (2) See headers of reformat scripts for non-obs4MIPs data for download instructions.

* AIRS (obs4MIPs): hus, husStderr
* AIRS-2-0 (obs4MIPs): hur
* CERES-EBAF (obs4MIPs): rsdt, rsut, rsutcs
* ERA-Interim (OBS6): hur, ta, va, wap
* GPCP-SG (obs4MIPs): pr
* HadCRUT4 (OBS): tasa
* HadISST (OBS): ts
* MLS-AURA (OBS6): hur
* TRMM-L3 (obs4MIPs): pr, prStderr


References
----------

* Brient, F., and T. Schneider, J. Climate, 29, 5821-5835, doi:10.1175/JCLIM-D-15-0897.1, 2016.
* Brient et al., Clim. Dyn., 47, doi:10.1007/s00382-015-2846-0, 2016.
* Cox et al., Nature, 553, doi:10.1038/nature25450, 2018.
* Gregory et al., Geophys. Res. Lett., 31,  doi:10.1029/2003GL018747, 2004.
* Lipat et al., Geophys. Res. Lett., 44, 5739-5748, doi:10.1002/2017GL73151, 2017.
* Sherwood et al., nature, 505, 37-42, doi:10.1038/nature12829, 2014.
* Su, et al., J. Geophys. Res. Atmos., 119, doi:10.1002/2014JD021642, 2014.
* Tian, Geophys. Res. Lett., 42, 4133-4141, doi:10.1002/2015GL064119, 2015.
* Volodin, Izvestiya, Atmospheric and Oceanic Physics, 44, 288-299, doi:10.1134/S0001433808030043, 2008.
* Zhai, et al., Geophys. Res. Lett., 42,  doi:10.1002/2015GL065911, 2015.

Example plots
-------------

.. _fig_ec_ecs_1:
.. figure::  /recipes/figures/emergent_constraints/ltmi.png
   :align:   center

   Lower tropospheric mixing index (LTMI; Sherwood et al., 2014) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_2:
.. figure::  /recipes/figures/emergent_constraints/shhc.png
   :align:   center

   Climatological Hadley cell extent (Lipat et al., 2017) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_3:
.. figure::  /recipes/figures/emergent_constraints/humidx.png
   :align:   center

   Tropical mid-tropospheric humidity asymmetry index (Tian, 2015) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_4:
.. figure::  /recipes/figures/emergent_constraints/itczidx.png
   :align:   center

   Southern ITCZ index (Tian, 2015) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_5:
.. figure::  /recipes/figures/emergent_constraints/covrefl.png
   :align:   center

   Covariance of shortwave cloud reflection (Brient and Schneider, 2016) vs.
   equilibrium climate sensitivity from CMIP5 models.

.. _fig_ec_ecs_6:
.. figure::  /recipes/figures/emergent_constraints/volodin.png
   :align:   center

   Difference in total cloud fraction between tropics (28°S - 28°N) and
   Southern midlatitudes (56°S - 36°S) (Volodin, 2008) vs. equilibrium climate
   sensitivity from CMIP5 models.
.. _recipes_runoff_et:

Runoff, Precipitation, Evapotranspiration
=========================================

Overview
--------
This diagnostic calculates biases of long-term climatological annual means of total runoff R,
precipitation P and evapotranspiration E for 12 large-scale catchments on different continents
and climates. For total runoff, catchment averaged model values are compared to climatological
GRDC station observations of river runoff (Duemenil Gates et al., 2000). Due to the incompleteness
of these station data, a year-to-year correspondence of data cannot be achieved in a generalized way,
so that only climatological data are considered, such it has been done in Hagemann, et al. (2013).
For precipitation, catchment-averaged WFDEI precipitation data (Weedon et al., 2014) from 1979-2010
is used as reference. For evapotranspiration, observations are estimated using the difference of the
above mentioned precipitation reference minus the climatological GRDC river runoff.

The catchments are Amazon, Congo, Danube, Ganges-Brahmaputra, Lena, Mackenzie, Mississippi, Murray,
Niger, Nile, Parana and Yangtze-Kiang. Variable names are expected to follow CMOR standard, e.g.
precipitation as pr, total runoff as mrro and evapotranspiration as evspsbl with all fluxes given in
kg m-2 s-1 . Evapotranspiration furthermore has to be defined positive upwards.

The diagnostic produces text files with absolute and relative bias to the observations, as well as the
respective absolute values. Furthermore it creates a bar plot for relative and absolute bias,
calculates and plots biases in runoff coefficient (R/P) and evapotranspiration coefficient (E/P) and
saves everything as one pdf file per model or one png file per model and analysis.

The bias of the runoff coefficient is calculated via:
:math:`C_R = \frac{R_{model}}{P_{model}} - \frac{R_{GRDC}}{P_{WFDEI}}` and similar for the
evapotranspiration coefficient. In a very first approximation, evapotranspiration
and runoff are determined only by precipitation. In other words :math:`R = P - E`. Hence, the runoff coefficient
(and similar the evapotranspiration coefficient) tells you how important runoff (or evapotranspiration)
is in this region. By plotting the bias of the runoff coefficient against the evapotranspiration coefficient
we can immediately see whether there is a shift from runoff to evapotranspiration. On the other hand, by
plotting the bias of the runoff coefficient against the relative bias of precipitation we can see whether
an error in runoff is due to an error in precipitation.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_runoff_et.yml

Diagnostics are stored in diag_scripts/runoff_et/

    * catchment_analysis.py: bar and scatter plots for catchment averages of
      runoff, evapotranspiration and precipitation


User settings in recipe
-----------------------

#. Script catchment_analysis.py

   *Required settings (scripts)*

   * catchmentmask: netCDF file indicating the grid cell for a specific catchment. Modus of
     distribution not yet clearified. ESGF?

   *Optional settings (variables)*

   * reference_dataset: dataset_name
     Datasets can be used as reference instead of defaults provided with the diagnostics.
     Must be identical for all variables.


Variables
---------

* evspsbl (atmos, monthly mean, time latitude longitude)
* pr      (atmos, monthly mean, time latitude longitude)
* mrro    (land,  monthly mean, time latitude longitude)


Observations and reformat scripts
---------------------------------

Default reference data based on GRDC and WFDEI are included in the diagnostic script
as catchment averages. They can be replaced with any gridded dataset by defining a
reference_dataset. The necessary catchment mask is available at

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2025776.svg
   :target: https://doi.org/10.5281/zenodo.2025776

All other datasets are remapped onto the catchment mask grid as part
of the diagnostics.


References
----------
* Duemenil Gates, L., S. Hagemann and C. Golz,
  Observed historical discharge data from major rivers for climate model validation.
  Max Planck Institute for Meteorology Report 307, Hamburg, Germany, 2000.

* Hagemann, S., A. Loew, A. Andersson,
  Combined evaluation of MPI-ESM land surface water and energy fluxes
  J. Adv. Model. Earth Syst., 5, doi:10.1029/2012MS000173, 2013.

* Weedon, G. P., G. Balsamo, N. Bellouin, S. Gomes, M. J. Best, and P. Viterbo,
  The WFDEI meteorological forcing data set: WATCH Forcing Data methodology applied
  to ERA‐Interim reanalysis data,
  Water Resour. Res., 50, 7505–7514, doi: 10.1002/2014WR015638, 2014


Example plots
-------------

.. _fig_runoff_et_1:
.. figure::  /recipes/figures/runoff_et/catchments.png
   :align:   center
   :width:   14cm

   Catchment definitions used in the diagnostics.

.. _fig_runoff_et_2:
.. figure::  /recipes/figures/runoff_et/MPI-ESM-LR_historical_r1i1p1_bias-plot_mrro.png
   :align:   center
   :width:   14cm

   Barplot indicating the absolute and relative bias in annual runoff between MPI-ESM-LR (1970-2000)
   and long term GRDC data for specific catchments.

.. _fig_runoff_et_3:
.. figure::  /recipes/figures/runoff_et/MPI-ESM-LR_historical_r1i1p1_rocoef-vs-relprbias.png
   :align:   center
   :width:   14cm

   Biases in runoff coefficient (runoff/precipitation) and precipitation for major catchments of
   the globe. The MPI-ESM-LR historical simulation (1970-2000) is used as an example.
.. _recipes_esacci_lst:

ESA CCI LST comparison to Historical Models
===========================================

Overview
--------

This diagnostic compares ESA CCI LST to multiple historical emsemble members of CMIP models.
It does this over a defined region for monthly values of the land surface temperature.
The result is a plot showing the mean differnce of CCI LST to model average LST, with a region of +/- one standard deviation of the model mean LST given as a measure of model variability.

The recipe and diagnostic need the all time average monthly LST from the CCI data.
We use the L3C single sensor monthy data.
A CMORizing script calculates the mean of the day time, and night time overpasses to give the all time average LST.
This is so that the Amon output from CMIP models can be used.
We created such a dataset from the Aqua MODIS data from CCI.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * ``recipe_esacci_lst.yml``

Diagnostics are stored in esmvaltool/diag_scripts/lst/

    * ``lst.py``


User settings in recipe
-----------------------

#. Script ``recipe_esacci_lst.yml``

   *No required settings for script*
  
   *No user defined inputs to the diagnostic*

   *Required settings for variables*
    
    * The diagnostic works with all data sources on having the same start_year and end_year, and hence that data is also available.

   *Required settings for preprocessor*
     
    * start_longitude, end_longitude The western and eastern bounds of the region to work with.
    * start_latitude, end_latitude The southern and northern bounds of the region to work with.
    * target_grid This should be one of the model grids.
   

Variables
---------

* ts (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

This recipe and diagnostic is written to work with data created from the CMORizer esmvaltool/cmorizers/obs/cmorize_obs_esacci_lst.py.
This takes the orginal ESA CCI LST files for the L3C data from Aqua MODIS DAY and NIGHT files and creates a the all time mean data this diagnostic uses.
Advice from the CCI LST team is to use the monthly not daily files to create the all time average to avoid th epossibility of biasing towards night time LST values being more prevalent because of how the cloud screening algorithms work.

References
----------

* ESA CCI LST project https://climate.esa.int/en/projects/land-surface-temperature/

Example plots
-------------

.. _fig_lst_example:
.. figure::  /recipes/figures/lst/lst_example.png
   :align:   center

   Timeseries of the ESA CCI LST minus mean of CMIP6 ensembles. The selected region is 35E-175E, 55N-70N.
   The black line is the mean difference, and the blue shaded area denotes one standard deviation either way of the individual ensemble member's differecen in LST.
   Models used for this are UKESM1 members r1i1p1f2 and r2i1p1f2, and CESM members r2i1p1f1 and r3i1p1f1.
   We have used the entire timeseries of available CCI data 2004-2014 inclusive, noting we have not written the CMORizer to process the incomplete year of 2003 for the Aqua MODIS data.
.. _recipes_esaccioc:
   
Ocean chlorophyll in ESMs compared to ESA-CCI observations.
===========================================================

Overview
--------

This recipe compares monthly surface chlorophyll from CMIP models to ESA CCI ocean colour chlorophyll (ESACCI-OC). The observations are the merged sensor geographic monthly L3S chlor_a data Sathyendranath et al. (2019). Multiple models and different observational versions can be used by the script.

The recipe_esacci_oc.yml produces an image showing four maps. Each of these four maps shows latitude vs longitude and the chlorophyll value. The four plots are: ESACCI-OC v5.0 chlorophyll, the CMIP6 model, the bias model-observation and the ratio model/observations. The script also produces a scatter plot for all coordinates with the model on the x-axis and the observations on the y axis and a line of best fit with the parameter values given in the panel.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/ocean/

    * recipe_esacci_oc.yml

Diagnostics are stored in esmvaltool/diag_scripts/ocean/

    * diagnostic_model_vs_obs.py


User settings in recipe
-----------------------

#. Script diagnostic_model_vs_obs.py

   *Required settings for script*

   * observational_dataset: name of reference dataset (e.g. {dataset: ESACCI-OC,})


Variables
---------

* chl (ocean, monthly mean, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* ESACCI-OC (chl)

  *Reformat script:* reformat_scripts/obs/reformat_obs_esacci_oc.py


References
----------

* Sathyendranath, S., et al. (2019), An ocean-colour time series for use in climate studies: the experience of the Ocean-Colour Climate Change Initiative (OC-CCI). Sensors: 19, 4285. doi:10.3390/s19194285.
* ESACCI-OC dataset: http://dx.doi.org/10.5285/00b5fc99f9384782976a4453b0148f49

Example plots
-------------

.. _fig_ocdiag_maps:
.. figure::  /recipes/figures/ocean/model_vs_obs_MassConcentrationofTotalPhytoplanktonExpressedasChlorophyllinSeaWater_NorESM2-LM_ESACCI-OC__maps.png
   :align:   center
   :width:   12cm

   Surface chlorophyll from ESACCI-OC ocean colour data version 5.0 and the
   CMIP6 model NorESM2-LM. This model overestimates chlorophyll compared to
   the observations.

.. _fig_ocdiag_scatter:
.. figure::  /recipes/figures/ocean/model_vs_obs_MassConcentrationofTotalPhytoplanktonExpressedasChlorophyllinSeaWater_NorESM2-LM_ESACCI-OC__scatter.png
   :align:   center
   :width:   8cm

   Scatter plot of surface chlorophyll from ESACCI-OC ocean colour data
   version 5.0 and the CMIP6 model NorESM2-LM.
.. _recipe_eyring13jgr:

Ozone and associated climate impacts
====================================

Overview
--------

This recipe is implemented into the ESMValTool to evaluate atmospheric chemistry and the climate impact of stratospheric ozone changes. It reproduces selected plots from Eyring et al. (2013).

The following plots are reproduced:

* Zonal mean of long-term zonal wind with linear trend

.. _`Eyring et al. (2013)`: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/jgrd.50316

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_eyring13jgr_12.yml

Diagnostics are stored in esmvaltool/diag_scripts/eyring13jgr/

* eyring13jgr_fig12.ncl

User settings in recipe
-----------------------
#. Preprocessor

   * ``zonal`` : Regridding and zonal mean used by eyring13jgr_fig12

#. Script <eyring13jgr_fig12.ncl>

   *Required settings for script*

   * ``e13fig12_exp_MMM``: name of the experiments for the MMM

   *Optional settings for script*

   * ``e13fig12_start_year``: year when to start the climatology calculation
   * ``e13fig12_end_year``: year when to end the climatology calculation
   * ``e13fig12_multimean``: calculate multimodel mean (default: False)
   * ``e13fig12_season``: season (default: ANN (annual))

   *Required settings for variables*
   
   * ``preprocessor``: zonal
   * ``reference_dataset``: name of the reference model or observation for regridding and bias calculation (e.g. ERA5).
   *  ``mip``:  Amon.

Variables
---------

*  ua (atmos, monthly mean, longitude latitude level time)

Observations and reformat scripts
---------------------------------

* ERA5
  *Reformatting with:* recipes/cmorizers/recipe_era5.yml


Example plots
-------------

.. _fig_eyring13jgr_12:
.. figure::  /recipes/figures/eyring13jgr/fig_eyr13jgr_12.png
   :align:   center

   Long-term mean (thin black contour) and linear trend (colour) of zonal mean DJF zonal winds for the multi-model mean CMIP6 over 1995-2014
.. _recipes_extreme_index:

Combined Climate Extreme Index
==============================

Overview
--------

The goal of this diagnostic is to compute time series of a number of extreme events: heatwave, coldwave, heavy precipitation, drought and high wind. Then, the user can combine these different components (with or without weights). The result is an index similar to the Climate Extremes Index (CEI; Karl et al., 1996), the modified CEI (mCEI; Gleason et al., 2008) or the Actuaries Climate Index (ACI; American Academy of Actuaries, 2018). The output consists of a netcdf file containing the area-weighted and multi-model multi-metric index. This recipe can be applied to data with any temporal resolution, and the running average is computed based on the user-defined window length (e.g. a window length of 5 would compute the 5-day running mean when applied to data, or 5-month running mean when applied to monthly data).

In recipe_extreme_index.yml, after defining the area and reference and projection period, the weigths for each metric selected. The options are

* weight_t90p the weight of the number of days when the maximum temperature exceeds the 90th percentile,

* weight_t10p the weight of the number of days when the minimum temperature falls below the 10th percentile,

* weight_Wx the weight of the number of days when wind power (third power of wind speed) exceeds the 90th percentile,

* weight_cdd the weight of the maximum length of a dry spell, defined as the maximum number of consecutive days when the daily precipitation is lower than 1 mm, and

* weight_rx5day the weight of the maximum precipitation accumulated during 5 consecutive days.

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_extreme_index.yml

Diagnostics are stored in diag_scripts/magic_bsc/

* extreme_index.R


User settings
-------------

User setting files are stored in recipes/

#. recipe_extreme_index.yml

   *Required settings for script*

   *   weight_t90p: 0.2 (from 0 to 1, the total sum of the weight should be 1)
   *   weight_t10p: 0.2 (from 0 to 1, the total sum of the weight should be 1)
   *   weight_Wx: 0.2 (from 0 to 1, the total sum of the weight should be 1)
   *   weight_rx5day: 0.2 (from 0 to 1, the total sum of the weight should be 1)
   *   weight_cdd: 0.2 (from 0 to 1, the total sum of the weight should be 1)
   *   running_mean: 5 (depends on the length of the future projection period selected, but recommended not greater than 11)

Variables
---------

* tasmax (atmos, daily, longitude, latitude, time)
* tasmin (atmos, daily, longitude, latitude, time)
* sfcWind (atmos, daily, longitude, latitude, time)
* pr (atmos, daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Alexander L.V.  and Coauthors (2006). Global observed changes in daily climate extremes of temperature and precipitation. J. Geophys. Res., 111, D05109. https://doi.org/10.1029/2005JD006290

* American Academy of Actuaries, Canadian Institute of Actuaries, Casualty Actuarial Society and Society of Actuaries. Actuaries Climate Index. http://actuariesclimateindex.org (2018-10-06).

* Donat, M., and Coauthors (2013). Updated analyses of temperature and precipitation extreme indices since the beginning of the twentieth century: The HadEX2 dataset. J.  Geophys. Res., 118, 2098–2118, https://doi.org/10.1002/jgrd.50150.

* Fouillet, A., Rey, G., Laurent, F., Pavillon, G. Bellec, S., Guihenneuc-Jouyaux, C., Clavel J., Jougla, E. and Hémon, D. (2006) Excess mortality related to the August 2003 heat wave in France. Int. Arch. Occup. Environ. Health, 80, 16–24. https://doi.org/10.1007/s00420-006-0089-4

* Gleason, K.L., J.H. Lawrimore, D.H. Levinson, T.R. Karl, and D.J. Karoly (2008). A Revised U.S. Climate Extremes Index. J. Climate, 21, 2124-2137 https://doi.org/10.1175/2007JCLI1883.1

* Meehl, G. A., and Coauthors (2000). An introduction to trends inextreme weather and climate events: Observations, socio-economic impacts, terrestrial ecological impacts, and model projections. Bull. Amer. Meteor. Soc., 81, 413–416. `doi: 10.1175/1520-0477(2000)081<0413:AITTIE>2.3.CO;2 <https://journals.ametsoc.org/doi/abs/10.1175/1520-0477%282000%29081%3C0413%3AAITTIE%3E2.3.CO%3B2>`_

* Whitman, S., G. Good, E. R. Donoghue, N. Benbow, W. Y. Shou and S. X. Mou (1997). Mortality in Chicago attributed to the July 1995 heat wave. Amer. J. Public Health, 87, 1515–1518. https://doi.org/10.2105/AJPH.87.9.1515

* Zhang, Y., M. Nitschke, and P. Bi (2013). Risk factors for direct heat-related hospitalization during the 2009 Adelaide heat-wave: A case crossover study. Sci. Total Environ., 442, 1–5. https://doi.org/10.1016/j.scitotenv.2012.10.042

* Zhang, X. , Alexander, L. , Hegerl, G. C., Jones, P. , Tank, A. K.,  Peterson, T. C., Trewin, B.  and Zwiers, F. W. (2011). Indices for  monitoring changes in extremes based on daily temperature and  precipitation data. WIREs Clim Change, 2: 851-870. doi:10.1002/wcc.147. https://doi.org/10.1002/wcc.147



Example plots
-------------

.. _fig_combinedindices1:
.. figure::  /recipes/figures/combined_climate_extreme_index/t90p_IPSL-CM5A-LR_rcp85_2020_2040.png
   :align:   center
   :width:   14cm

Average change in the heat component (t90p metric) of the Combined Climate Extreme Index for the 2020-2040 compared to the 1971-2000 reference period for the RCP 8.5 scenario simulated by MPI-ESM-MR.
.. _recipes_rainfarm:

RainFARM stochastic downscaling
===============================


Overview
--------

Precipitation extremes and small-scale variability are essential drivers in many climate change impact studies. However, the spatial resolution currently achieved by global and regional climate models is still insufficient to correctly identify the fine structure of precipitation intensity fields. In the absence of a proper physically based representation, this scale gap can be at least temporarily bridged by adopting a stochastic rainfall downscaling technique (Rebora et al, 2006). With this aim, the Rainfall Filtered Autoregressive Model (RainFARM) was developed to apply the stochastic precipitation downscaling method to climate models. The RainFARM Julia library and command-line tool version (https://github.com/jhardenberg/RainFARM.jl) was implemented as recipe. The stochastic method allows to predict climate variables at local scale from information simulated by climate models at regional scale: It first evaluates the statistical distribution of precipitation fields at regional scale and then applies the relationship to the boundary conditions of the climate model to produce synthetic fields at the requested higher resolution. RainFARM exploits the nonlinear transformation of a Gaussian random precipitation field, conserving the information present in the fields at larger scale (Rebora et al., 2006; D’Onofrio et al., 2014).


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_rainfarm.yml

Diagnostics are stored in diag_scripts/rainfarm/

* rainfarm.jl


User settings
-------------

*Required settings for script*

* slope: spatial spectral slope (set to 0 to compute automatically from large scales)
* nens: number of ensemble members to be calculated
* nf: number of subdivisions for downscaling (e.g. 8 will produce output fields with linear resolution increased by a factor 8)
* conserv_glob: logical, if to conserve precipitation over full domain
* conserv_smooth: logical, if to conserve precipitation using convolution (if neither conserv_glob or conserv_smooth is chosen, box conservation is used)
* weights_climo: set to false or omit if no orographic weights are to be used, else set it to the path to a fine-scale precipitation climatology file. If a relative file path is used, `auxiliary_data_dir` will be searched for this file. The file is expected to be in NetCDF format and should contain at least one precipitation field. If several fields at different times are provided, a climatology is derived by time averaging. Suitable climatology files could be for example a fine-scale precipitation climatology from a high-resolution regional climate model (see e.g. Terzago et al. 2018), a local high-resolution gridded climatology from observations, or a reconstruction such as those which can be downloaded from the WORLDCLIM (http://www.worldclim.org) or CHELSA (http://chelsa-climate.org) websites. The latter data will need to be converted to NetCDF format before being used (see for example the GDAL tools (https://www.gdal.org).


Variables
---------

* pr (atmos, daily mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

None.


References
----------

* Terzago et al. 2018, Nat. Hazards Earth Syst. Sci., 18, 2825-2840
* D'Onofrio et al. 2014, J of Hydrometeorology 15, 830-843
* Rebora et. al 2006, JHM 7, 724

Example plots
-------------

.. figure:: /recipes/figures/rainfarm/rainfarm.png
   :width: 14cm

   Example of daily cumulated precipitation from the CMIP5 EC-EARTH model on a specific day, downscaled using RainFARM from its original resolution (1.125°) (left panel), increasing spatial resolution by a factor of 8 to 0.14°; Two stochastic realizations are shown (central and right panel). A fixed spectral slope of s=1.7 was used. Notice how the downscaled fields introduce fine scale precipitation structures, while still maintaining on average the original coarse-resolution precipitation. Different stochastic realizations are shown to demonstrate how an ensemble of realizations can be used to reproduce unresolved subgrid variability. (N.B.: this plot was not produced by ESMValTool - the recipe output is netcdf only). 
.. _recipes_zmnam:

Stratosphere-troposphere coupling and annular modes indices (ZMNAM)
===================================================================


Overview
--------

The current generation of climate models include the representation of stratospheric processes, as the vertical coupling with the troposphere is important for the weather and climate at the surface (e.g., `Baldwin and Dunkerton, 2001 <https://doi.org/10.1126/science.1063315>`_).

The recipe recipe_zmnam.yml can be used to evaluate the representation of the Northern Annular Mode (NAM, e.g., `Wallace, 2000 <https://doi.org/10.1002/qj.49712656402>`_) in climate simulations, using reanalysis datasets as reference.

The calculation is based on the “zonal mean algorithm” of `Baldwin and Thompson (2009) <https://doi.org/10.1002/qj.479>`_, and is alternative to pressure based or height-dependent methods.

This approach provides a robust description of the stratosphere-troposphere coupling on daily timescales, requiring less subjective choices and a reduced amount of input data.
Starting from daily mean geopotential height on pressure levels, the leading empirical orthogonal function/principal component are computed from zonal mean daily anomalies, with the leading principal component representing the zonal mean NAM index. The regression of the monthly mean geopotential height onto this monthly averaged index represents the NAM pattern for each selected pressure level.

The outputs of the procedure are the monthly time series and the histogram of the daily zonal-mean NAM index, and the monthly regression maps for selected pressure levels. The users can select the specific datasets (climate model simulation and/or reanalysis) to be evaluated, and a subset of pressure levels of interest.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_zmnam.yml

Diagnostics are stored in diag_scripts/zmnam/

* zmnam.py

and subroutines

* zmnam_calc.py
* zmnam_plot.py
* zmnam_preproc.py


User settings
-------------

None.


Variables
---------

* zg (atmos, daily mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

None.


References
----------

* Baldwin, M. P. and Thompson, D. W. (2009), A critical comparison of stratosphere–troposphere coupling indices. Q.J.R. Meteorol. Soc., 135: 1661-1672. `doi:10.1002/qj.479 <https://doi.org/10.1002/qj.479>`_.
* Baldwin, M. P and Dunkerton, T. J. (2001), Stratospheric Harbingers of Anomalous Weather Regimes. Science  294 (5542): 581-584. `doi:10.1126/science.1063315 <https://doi.org/10.1126/science.1063315>`_.
* Wallace, J. M. (2000), North Atlantic Oscillation/annular mode: Two paradigms-one phenomenon. Q.J.R. Meteorol. Soc., 126 (564): 791-805. `doi:10.1002/qj.49712656402 <https://doi.org/10.1002/qj.49712656402>`_.



Example plots
-------------

.. figure:: /recipes/figures/zmnam/zmnam_reg.png
   :width: 10cm

   Regression map of the zonal-mean NAM index onto geopotential height, for a selected pressure level (250 hPa) for the MPI-ESM-MR model (CMIP5 AMIP experiment, period 1979-2008). Negative values are shaded in grey.

.. figure:: /recipes/figures/zmnam/zmnam_ts.png
   :width: 10cm

   Time series of the zonal-mean NAM index for a selected pressure level (250 hPa) for the MPI-ESM-MR model (CMIP5 AMIP experiment, period 1979-2008).
.. _nml_collins:

IPCC AR5 Chapter 12 (selected figures)
======================================

Overview
--------

The goal is to create a standard recipe for creating selected Figures from
IPCC AR5 Chapter 12 on "Long-term Climate Change: Projections, Commitments
and Irreversibility". These include figures showing the change in a variable
between historical and future periods, e.g. maps (2D variables), zonal means
(3D variables), timeseries showing the change in certain variables from
historical to future periods for multiple scenarios, and maps visualizing
change in variables normalized by global mean temperature change (pattern
scaling) as in Collins et al., 2013.


Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_collins13ipcc.yml

Diagnostics are stored in diag_scripts/

* ipcc_ar5/ch12_map_diff_each_model_fig12-9.ncl: calculates the difference between
  future and historical runs for one scenario for each given model
  individually on their native grid and plots all of them in one Figure.
  As in Figure 12.9 in AR5.
* ipcc_ar5/ch12_ts_line_mean_spread.ncl: calculates time series for one variable,
  change in future relative to base period in historical, multi-model mean as
  well as spread around it (as standard deviation).
* ipcc_ar5/ch12_plot_ts_line_mean_spread.ncl: plots the timeseries multi-model mean
  and spread calculated above. As in Figure 12.5 in AR5.
* ipcc_ar5/ch12_calc_IAV_for_stippandhatch.ncl: calculates the interannual variability
  over piControl runs, either over the whole time period or in chunks over
  some years.
* ipcc_ar5/ch12_calc_map_diff_mmm_stippandhatch.ncl: calculates the difference between
  future and historical periods for each given model and then calculates
  multi-model mean as well as significance. Significant is where the
  multi-model mean change is greater than two standard deviations of the
  internal variability and where at least 90% of the models agree on the
  sign of change. Not significant is where the multi-model mean change is
  less than one standard deviation of internal variability.
* ipcc_ar5/ch12_plot_map_diff_mmm_stipp.ncl: plots multi-model mean maps calculated
  above including stippling where significant and hatching where not
  significant. As in Figure 12.11 in AR5.
* ipcc_ar5/ch12_calc_zonal_cont_diff_mmm_stippandhatch.ncl: calculates zonal means
  and the difference between future and historical periods for each given
  model and then calculates multi-model mean as well as significance as above.
* ipcc_ar5/ch12_plot_zonal_diff_mmm_stipp.ncl: plots the multi-model mean zonal plots
  calculated above including stippling where significant and hatching where
  not significant. As in Figure 12.12 in AR5.
* ipcc_ar5/ch12_calc_map_diff_scaleT_mmm_stipp.ncl: calculates the change in variable
  between future and historical period normalized by gloabl mean temperature
  change of each given model and scenario. Then averages over all realizations
  and calculates significance. Significant is where the mean change averaged
  over all realizations is larger than the 95% percentile of the distribution
  of models (assumed to be gaussian). Can be plotted using
  ipcc_ar5/ch12_plot_map_diff_mmm_stipp.ncl.
* seaice/seaice_ecs.ncl: scatter plot of historical trend in September
  Arctic sea ice extent (SSIE) vs historical long-term mean SSIE (similar to
  Fig. 12.31a in AR5) and historical SSIE trend vs YOD RCP8.5 (similar to Fig. 12.31d
  in AR5).
* seaice/seaice_yod.ncl: calculation of year of near disappearance of Arctic sea ice
  (similar to Fig 12.31e in AR5)
* ipcc_ar5/ch12_snw_area_change_fig12-32.ncl: calculate snow area extent in a region
  (e.g Northern Hemisphere) and season (e.g. Northern Hemisphere spring March
  & April) relative to a reference period (e.g 1986-2005) and spread over
  models as in Fig. 12.32 of IPCC AR5. Can be plotted using
  ipcc_ar5/ch12_plot_ts_line_mean_spread.ncl.

User settings
-------------

#. Script ipcc_ar5/ch12_map_diff_each_model_fig12-9.ncl

   *Required settings (script)*

   * time_avg: time averaging ("annualclim", "seasonalclim")
   * experiment: IPCC Scenario, used to pair historical and rcp runs from
     same model

   *Optional settings (script)*

   * projection: map projection, any valid ncl projection, default = Robinson
   * max_vert: maximum number of plots in vertical
   * max_hori: maximum number of plots in horizontal
   * title: plot title
   * colormap: alternative colormap, path to rgb file or ncl name
   * diff_levs: list with contour levels for plots
   * span: span whole colormap? (True, False, default = False)

   *Required settings (variables)*

   * project: CMIP5 (or CMIP6)
   * mip: variable mip, generally Amon or Omon

#. Script ipcc_ar5/ch12_ts_line_mean_spread.ncl

   *Required settings (script)*

   * scenarios: list with scenarios included in figure
   * syears: list with start years in time periods (e.g. start of historical
     period and rcps)
   * eyears: list with end years in time periods (end year of historical runs
     and rcps)
   * begin_ref_year: start year of reference period (e.g. 1986)
   * end_ref_year: end year of reference period (e.g 2005)
   * label: list with labels to use in legend depending on scenarios

   *Optional settings (script)*

   * spread: how many standard deviations to calculate the spread with?
     default is 1., ipcc tas used 1.64
   * model_nr: save number of model runs per period and scenario in netcdf
     to print in plot? (True, False, default = False)
   * ts_minlat: minimum latitude if not global
   * ts_maxlat: maximum latitude if not global
   * ts_minlon: minimum longitude if not global
   * ts_maxlon: maximum longitude if not global

   *Required settings (variables)*

   * project: CMIP5 (or CMIP6)
   * mip: variable mip, generally Amon or Omon

#. Script ipcc_ar5/ch12_plot_ts_line_mean_spread.ncl:

   *Required settings (script)*

   * ancestors: variable and diagnostics that calculated data to be plotted

   *Optional settings (script)*

   * title: specify plot title
   * yaxis: specify y-axis title
   * ymin: minimim value on y-axis, default calculated from data
   * ymax: maximum value on y-axis
   * colormap: alternative colormap, path to rgb file or ncl name

#. Script ipcc_ar5/ch12_calc_IAV_for_stippandhatch.ncl:

   *Required settings (script)*

   * time_avg: time averaging ("annualclim", "seasonalclim"), needs to be
     consistent with calculation in ch12_calc_map_diff_mmm_stippandhatch.ncl

   *Optional settings (script)*

   * periodlength: length of period in years to calculate variability over,
     default is total time period
   * iavmode: calculate IAV from multi-model mean or save individual models
     ("each": save individual models, "mmm": multi-model mean, default),
     needs to be consistent with ch12_calc_map_diff_mmm_stippandhatch.ncl

   *Required settings (variables)*

   * project: CMIP5 (or CMIP6)
   * mip: variable mip, generally Amon or Omon
   * exp: piControl
   * preprocessor: which preprocessor to use, depends on dimension of variable,
     for 2D preprocessor only needs to regrid, for 3D we need to extract levels
     either based on reference_dataset or specify levels.

   *Optional settings (variables)*

   * reference_dataset: the reference dataset for level extraction in case of
     3D variables.

#. Script ipcc_ar5/ch12_calc_map_diff_mmm_stippandhatch.ncl:

   *Required settings (script)*

   * ancestors: variable and diagnostics that calculated interannual
     variability for stippling and hatching
   * time_avg: time averaging ("annualclim", "seasonalclim")
   * scenarios: list with scenarios to be included
   * periods: list with start years of periods to be included
   * label: list with labels to use in legend depending on scenarios

   *Optional settings (script)*

   * seasons: list with seasons index if time_avg "seasonalclim" (then
     required),  DJF:0, MAM:1, JJA:2, SON:3
   * iavmode: calculate IAV from multi-model mean or save individual models
     ("each": save individual models, "mmm": multi-model mean, default),
     needs to be consistent with ch12_calc_IAV_for_stippandhatch.ncl
   * percent: determines if difference expressed in percent (0, 1, default = 0)

   *Required settings (variables)*

   * project: CMIP5 (or CMIP6)
   * mip: variable mip, generally Amon or Omon
   * preprocessor: which preprocessor to use, preprocessor only needs to regrid

#. Script ipcc_ar5/ch12_plot_map_diff_mmm_stipp.ncl:

   *Required settings (script)*

   * ancestors: variable and diagnostics that calculated field to be plotted

   *Optional settings (script)*

   * projection: map projection, any valid ncl projection, default = Robinson
   * diff_levs: list with explicit levels for all contour plots
   * max_vert: maximum number of plots in vertical
   * max_hori: maximum number of plots in horizontal
   * model_nr: save number of model runs per period and scenario in netcdf to
     print in plot? (True, False, default = False)
   * colormap: alternative colormap, path to rgb file or ncl name
   * span: span whole colormap? (True, False, default = True)
   * sig: plot stippling for significance? (True, False)
   * not_sig: plot hatching for uncertainty? (True, False)
   * pltname: alternative name for output plot, default is diagnostic +
     varname + time_avg
   * units: units written next to colorbar, e.g (~F35~J~F~C)

#. Script ipcc_ar5/ch12_calc_zonal_cont_diff_mmm_stippandhatch.ncl:

   *Required settings (script)*

   * ancestors: variable and diagnostics that calculated interannual
     variability for stippling and hatching
   * time_avg: time averaging ("annualclim", "seasonalclim")
   * scenarios: list with scenarios to be included
   * periods: list with start years of periods to be included
   * label: list with labels to use in legend depending on scenarios

   *Optional settings (script)*

   * base_cn: if want contours of base period as contour lines, need to save
     base period field (True, False)
   * seasons: list with seasons index if time_avg "seasonalclim" (then
     required),  DJF:0, MAM:1, JJA:2, SON:3
   * iavmode: calculate IAV from multi-model mean or save individual models
     ("each": save individual models, "mmm": multi-model mean, default),
     needs to be consistent with ch12_calc_IAV_for_stippandhatch.ncl
   * percent: determines if difference expressed in percent (0, 1, default = 0)

   *Required settings (variables)*

   * project: CMIP5 (or CMIP6)
   * mip: variable mip, generally Amon or Omon
   * preprocessor: which preprocessor to use, preprocessor needs to regrid,
     extract leves and calculate the zonal mean.

   *Optional settings (variables)*

   * reference_dataset: the reference dataset for level extraction

#. Script ipcc_ar5/ch12_plot_zonal_diff_mmm_stipp.ncl:

   *Required settings (script)*

   * ancestors: variable and diagnostics that calculated field to be plotted

   *Optional settings (script)*

   * diff_levs: list with explicit levels for all contour plots
   * max_vert: maximum number of plots in vertical
   * max_hori: maximum number of plots in horizontal
   * model_nr: save number of model runs per period and scenario in netcdf to
     print in plot? (True, False, default = False)
   * colormap: alternative colormap, path to rgb file or ncl name
   * span: span whole colormap? (True, False, default = True)
   * sig: plot stippling for significance? (True, False)
   * not_sig: plot hatching for uncertainty? (True, False)
   * pltname: alternative name for output plot, default is diagnostic +
     varname + time_avg
   * units: units written next to colorbar in ncl strings, e.g (m s~S~-1~N~)
   * if base_cn: True in ch12_calc_zonal_cont_diff_mmm_stippandhatch.ncl
     further settings to control contour lines:

     * base_cnLevelSpacing: spacing between contour levels
     * base_cnMinLevel: minimum contour line
     * base_cnMaxLevel: maximum contour line

#. Script ipcc_ar5/ch12_calc_map_diff_scaleT_mmm_stipp.ncl:

   *Required settings (script)*

   * time_avg: time averaging ("annualclim", "seasonalclim")
   * scenarios: list with scenarios to be included
   * periods: list with start years of periods to be included
   * label: list with labels to use in legend depending on scenarios

   *Optional settings (script)*

   * seasons: list with seasons index if time_avg "seasonalclim"
     (then required),  DJF:0, MAM:1, JJA:2, SON:3
   * percent: determines if difference expressed in percent (0, 1, default = 0)

   *Required settings (variables)*

   * project: CMIP5 (or CMIP6)
   * mip: variable mip, generally Amon or Omon
   * preprocessor: which preprocessor to use, preprocessor only needs to regrid

#. Script ipcc_ar5/ch12_snw_area_change_fig12-32.ncl:

   *Required settings (script)*

   * scenarios: list with scenarios included in figure
   * syears: list with start years in time periods (e.g. start of historical
     period and rcps)
   * eyears: list with end years in time periods (end year of historical runs
     and rcps)
   * begin_ref_year: start year of reference period (e.g. 1986)
   * end_ref_year: end year of reference period (e.g 2005)
   * months: first letters of  months included in analysis? e.g. for MA
     (March + April) for Northern Hemisphere
   * label: list with labels to use in legend depending on scenarios

   *Optional settings (script)*

   * spread: how many standard deviations to calculate the spread with?
     default is 1., ipcc tas used 1.64
   * model_nr: save number of model runs per period and scenario in netcdf
     to print in plot? (True, False, default = False)
   * colormap: alternative colormap, path to rgb file or ncl name
   * ts_minlat: minimum latitude if not global
   * ts_maxlat: maximum latitude if not global
   * ts_minlon: minimum longitude if not global
   * ts_maxlon: maximum longitude if not global

   *Required settings (variables)*

   * project: CMIP5 (or CMIP6)
   * mip: variable mip, LImon
   * fx_files: [sftlf, sftgif]

#. Script seaice/seaice_ecs.ncl

   *Required settings (scripts)*

   * hist_exp: name of historical experiment (string)
   * month: selected month (1, 2, ..., 12) or annual mean ("A")
   * rcp_exp: name of RCP experiment (string)
   * region: region to be analyzed ( "Arctic" or "Antarctic")

   *Optional settings (scripts)*

   * fill_pole_hole: fill observational hole at North pole (default: False)
   * styleset: color style (e.g. "CMIP5")

   *Optional settings (variables)*

   * reference_dataset: reference dataset

#. Script seaice/seaice_yod.ncl

   *Required settings (scripts)*

   * month: selected month (1, 2, ..., 12) or annual mean ("A")
   * region: region to be analyzed ( "Arctic" or "Antarctic")

   *Optional settings (scripts)*

   * fill_pole_hole: fill observational hole at North pole, Default: False
   * wgt_file: netCDF containing pre-determined model weights

   *Optional settings (variables)*

   * ref_model: array of references plotted as vertical lines


Variables
---------

*Note: These are the variables tested and used in IPCC AR5. However, the code is flexible and in theory other variables of the same kind can be used.*

* areacello (fx, longitude latitude)
* clt (atmos, monthly mean, longitude latitude time)
* evspsbl (atmos, monthly mean, longitude latitude time)
* hurs (atmos, monthly mean, longitude latitude time)
* mrro (land, monthly mean, longitude latitude time)
* mrsos (land, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* psl (atmos, monthly mean, longitude latitude time)
* rlut, rsut, rtmt (atmos, monthly mean, longitude latitude time)
* sic (ocean-ice, monthly mean, longitude latitude time)
* snw (land, monthly mean, longitude latitude time)
* sos (ocean, monthly mean, longitude latitude time)
* ta (atmos, monthly mean, longitude latitude lev time)
* tas (atmos, monthly mean, longitude latitude time)
* thetao (ocean, monthly mean, longitude latitude lev time)
* ua (atmos, monthly mean, longitude latitude lev time)

Observations and reformat scripts
---------------------------------

* HadISST (sic - esmvaltool/utils/cmorizers/obs/cmorize_obs_HadISST.ncl)

Reference
---------

* Collins, M., R. Knutti, J. Arblaster, J.-L. Dufresne, T. Fichefet, P.
  Friedlingstein, X. Gao, W.J. Gutowski, T. Johns, G. Krinner, M. Shongwe, C.
  Tebaldi, A.J. Weaver and M. Wehner, 2013: Long-term Climate Change:
  Projections, Commitments and Irreversibility. In: Climate Change 2013: The
  Physical Science Basis. Contribution of Working Group I to the Fifth
  Assessment Report of the Intergovernmental Panel on Climate Change [Stocker,
  T.F., D. Qin, G.-K. Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels,
  \Y. Xia, V. Bex and P.M. Midgley (eds.)]. Cambridge University Press,
  Cambridge, United Kingdom and New York, NY, USA.


Example plots
-------------

.. figure:: /recipes/figures/collins13ipcc/collins_fig_1.png
   :width: 85%
   :align: center

   Surface air temperature change in 2081–2100 displayed as anomalies with
   respect to 1986–2005 for RCP4.5 from individual CMIP5 models.


.. figure:: /recipes/figures/collins13ipcc/collins_fig_2.png
   :width: 50%
   :align: center

   Time series of global annual mean surface air temperature anomalie
   (relative to 1986–2005) from CMIP5 concentration-driven experiments.

.. figure:: /recipes/figures/collins13ipcc/collins_fig_4.png
   :width: 70%
   :align: center

   Multi-model CMIP5 average percentage change in seasonal mean precipitation
   relative to the reference period 1986–2005 averaged over the periods
   2081–2100 and 2181–2200 under the RCP8.5 forcing scenario. Hatching
   indicates regions where the multi-model mean change is less than one
   standard deviation of internal variability. Stippling indicates regions
   where the multi-model mean change is greater than two standard deviations
   of internal variability and where at least 90% of models agree on the sign
   of change

.. figure:: /recipes/figures/collins13ipcc/collins_fig_3.png
   :width: 70%
   :align: center

   Temperature change patterns scaled to 1°C of global mean surface
   temperature change.

.. figure::  /recipes/figures/seaice/SSIE-MEAN_vs_YOD_sic_extend_Arctic_September_1960-2100.png
   :align:   center
   :width:   9cm

   Scatter plot of mean historical September Arctic sea ice extent vs 1st year of disappearance
   (RCP8.5) (similar to IPCC AR5 Chapter 12, Fig. 12.31a).

.. figure::  /recipes/figures/seaice/timeseries_rcp85.png
   :align:   center
   :width:   12cm

   Time series of September Arctic sea ice extent for individual CMIP5 models,
   multi-model mean and multi-model standard deviation, year of disappearance
   (similar to IPCC AR5 Chapter 12, Fig. 12.31e).
Ocean diagnostics toolkit
=============================

Welcome to the API documentation for the ocean diagnostics tool kit.
This toolkit is built to assist in the evaluation of models of the ocean.

This toolkit is part of ESMValTool v2.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk

.. toctree::

.. automodule:: esmvaltool.diag_scripts.ocean.diagnostic_maps
.. automodule:: esmvaltool.diag_scripts.ocean.diagnostic_maps_quad
.. automodule:: esmvaltool.diag_scripts.ocean.diagnostic_model_vs_obs
.. automodule:: esmvaltool.diag_scripts.ocean.diagnostic_profiles
.. automodule:: esmvaltool.diag_scripts.ocean.diagnostic_timeseries
.. automodule:: esmvaltool.diag_scripts.ocean.diagnostic_transects
.. automodule:: esmvaltool.diag_scripts.ocean.diagnostic_seaice
.. automodule:: esmvaltool.diag_scripts.ocean.diagnostic_tools
.. _api.esmvaltool.diag_scripts.mlr:

Machine Learning Regression (MLR) diagnostics
=============================================

This module provides various tools to create and evaluate MLR models for
arbitrary input variables.


Examples
--------

* :ref:`recipes_schlund20jgr`: Use Gradient Boosted Regression Tree (GBRT)
  algorithm to constrain projected Gross Primary Production (GPP) in RCP 8.5
  scenario using observations of process-based predictors.


Diagnostic scripts
------------------
.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.mlr/evaluate_residuals
   esmvaltool.diag_scripts.mlr/main
   esmvaltool.diag_scripts.mlr/mmm
   esmvaltool.diag_scripts.mlr/plot
   esmvaltool.diag_scripts.mlr/postprocess
   esmvaltool.diag_scripts.mlr/preprocess
   esmvaltool.diag_scripts.mlr/rescale_with_emergent_constraint


Auxiliary scripts
-----------------
.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.mlr/init
   esmvaltool.diag_scripts.mlr/custom_sklearn
   esmvaltool.diag_scripts.mlr/models
   esmvaltool.diag_scripts.mlr/models.gbr_base
   esmvaltool.diag_scripts.mlr/models.linear_base


.. _availableMLRModels:

Available MLR models
--------------------
.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.mlr/models.gbr_sklearn
   esmvaltool.diag_scripts.mlr/models.gbr_xgboost
   esmvaltool.diag_scripts.mlr/models.gpr_sklearn
   esmvaltool.diag_scripts.mlr/models.huber
   esmvaltool.diag_scripts.mlr/models.krr
   esmvaltool.diag_scripts.mlr/models.lasso
   esmvaltool.diag_scripts.mlr/models.lasso_cv
   esmvaltool.diag_scripts.mlr/models.lasso_lars_cv
   esmvaltool.diag_scripts.mlr/models.linear
   esmvaltool.diag_scripts.mlr/models.rfr
   esmvaltool.diag_scripts.mlr/models.ridge
   esmvaltool.diag_scripts.mlr/models.ridge_cv
   esmvaltool.diag_scripts.mlr/models.svr
.. _api:

ESMValTool Code API Documentation
=================================

ESMValTool is mostly used as a command line tool. However, it is also possible
to use (parts of) ESMValTool as a library. This section documents the public
API of ESMValTool.

.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.shared
   esmvaltool.diag_scripts
   esmvaltool.diag_scripts.ocean
   esmvaltool.diag_scripts.emergent_constraints
   esmvaltool.diag_scripts.mlr
Diagnostic scripts
==================

Various diagnostic packages exist as part of ESMValTool.

.. automodule:: esmvaltool.diag_scripts

.. automodule:: esmvaltool.diag_scripts.ocean
.. _api_shared:

Shared diagnostic script code
=============================

.. automodule:: esmvaltool.diag_scripts.shared

Plotting
--------

.. automodule:: esmvaltool.diag_scripts.shared.plot
.. _api.esmvaltool.diag_scripts.emergent_constraints:

Emergent constraints diagnostics
================================

This module provides various tools to evaluate emergent constraints for
arbitrary input variables.


Examples
--------

* :ref:`recipe_ecs_scatter`
* :ref:`recipes_cox18nature`
* :ref:`recipes_schlund20esd`


Diagnostic scripts
------------------
.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.emergent_constraints/cox18nature
   esmvaltool.diag_scripts.emergent_constraints/ecs_scatter
   esmvaltool.diag_scripts.emergent_constraints/multiple_constraints
   esmvaltool.diag_scripts.emergent_constraints/single_constraint


Auxiliary scripts
-----------------
.. toctree::
   :maxdepth: 1

   esmvaltool.diag_scripts.emergent_constraints/init
.. _api.esmvaltool.diag_scripts.mlr.plot:

Plotting functionalities
========================

.. automodule:: esmvaltool.diag_scripts.mlr.plot
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.mlr.evaluate_residuals:

Evaluate residuals
==================

.. automodule:: esmvaltool.diag_scripts.mlr.evaluate_residuals
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.mlr.models.gbr_xgboost:


Gradient Boosted Regression Trees (xgboost implementation)
==========================================================

.. automodule:: esmvaltool.diag_scripts.mlr.models.gbr_xgboost
.. _api.esmvaltool.diag_scripts.mlr.models.rfr:

Random Forest Regression
========================

.. automodule:: esmvaltool.diag_scripts.mlr.models.rfr
.. _api.esmvaltool.diag_scripts.mlr.custom_sklearn:

Custom extensions of sklearn functionalities
============================================

.. automodule:: esmvaltool.diag_scripts.mlr.custom_sklearn
.. _api.esmvaltool.diag_scripts.mlr.models.lasso:

LASSO Regression
================

.. automodule:: esmvaltool.diag_scripts.mlr.models.lasso
.. _api.esmvaltool.diag_scripts.mlr.models.krr:

Kernel Ridge Regression
=======================

.. automodule:: esmvaltool.diag_scripts.mlr.models.krr
.. _api.esmvaltool.diag_scripts.mlr.models.lasso_lars_cv:

LASSO Regression (using Least-angle Regression algorithm) with built-in CV
==========================================================================

.. automodule:: esmvaltool.diag_scripts.mlr.models.lasso_lars_cv
.. _api.esmvaltool.diag_scripts.mlr.models.ridge_cv:

Ridge Regression with built-in CV
=================================

.. automodule:: esmvaltool.diag_scripts.mlr.models.ridge_cv
.. _api.esmvaltool.diag_scripts.mlr.models.lasso_cv:

LASSO Regression with built-in CV
=================================

.. automodule:: esmvaltool.diag_scripts.mlr.models.lasso_cv
.. _api.esmvaltool.diag_scripts.mlr.models.gpr_sklearn:

Gaussian Process Regression (sklearn implementation)
====================================================

.. automodule:: esmvaltool.diag_scripts.mlr.models.gpr_sklearn
.. _api.esmvaltool.diag_scripts.mlr.models.linear_base:

Base class for Linear models
============================

.. automodule:: esmvaltool.diag_scripts.mlr.models.linear_base
.. _api.esmvaltool.diag_scripts.mlr.models.linear:

Linear Regression
=================

.. automodule:: esmvaltool.diag_scripts.mlr.models.linear
.. _api.esmvaltool.diag_scripts.mlr.preprocess:

Preprocessing functionalities
=============================

.. automodule:: esmvaltool.diag_scripts.mlr.preprocess
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.mlr.models.huber:

Huber Regression
================

.. automodule:: esmvaltool.diag_scripts.mlr.models.huber
.. _api.esmvaltool.diag_scripts.mlr.rescale_with_emergent_constraint:

Rescale data with emergent constraints
=======================================

.. automodule:: esmvaltool.diag_scripts.mlr.rescale_with_emergent_constraint
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.mlr.models.gbr_sklearn:

Gradient Boosted Regression Trees (sklearn implementation)
==========================================================

.. automodule:: esmvaltool.diag_scripts.mlr.models.gbr_sklearn
.. _api.esmvaltool.diag_scripts.mlr.mmm:

Multi-model means (MMM)
=======================

.. automodule:: esmvaltool.diag_scripts.mlr.mmm
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.mlr.models.ridge:

Ridge Regression
================

.. automodule:: esmvaltool.diag_scripts.mlr.models.ridge
.. _api.esmvaltool.diag_scripts.mlr.models:

MLRModel base class
===================

.. automodule:: esmvaltool.diag_scripts.mlr.models
.. _api.esmvaltool.diag_scripts.mlr.init:

Auxiliary functions for MLR scripts
===================================

.. automodule:: esmvaltool.diag_scripts.mlr
.. _api.esmvaltool.diag_scripts.mlr.models.gbr_base:

Base class for Gradient Boosted Regression models
=================================================

.. automodule:: esmvaltool.diag_scripts.mlr.models.gbr_base
.. _api.esmvaltool.diag_scripts.mlr.postprocess:

Postprocessing functionalities
==============================

.. automodule:: esmvaltool.diag_scripts.mlr.postprocess
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.mlr.main:

MLR main diagnostic
===================

.. automodule:: esmvaltool.diag_scripts.mlr.main
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.mlr.models.svr:

Support Vector Regression
=========================

.. automodule:: esmvaltool.diag_scripts.mlr.models.svr
.. _api.esmvaltool.diag_scripts.emergent_constraints.cox18nature:

Emergent constraint on ECS from global temperature variability
==============================================================

.. automodule:: esmvaltool.diag_scripts.emergent_constraints.cox18nature
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.emergent_constraints.ecs_scatter:

Calculation of emergent constraints on ECS
==========================================

.. automodule:: esmvaltool.diag_scripts.emergent_constraints.ecs_scatter
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.emergent_constraints.single_constraint:

Evaluate single emergent constraint
===================================

.. automodule:: esmvaltool.diag_scripts.emergent_constraints.single_constraint
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.emergent_constraints.multiple_constraints:

Evaluate multiple emergent constraints simultaneously
=====================================================

.. automodule:: esmvaltool.diag_scripts.emergent_constraints.multiple_constraints
   :no-members:
   :no-inherited-members:
   :no-show-inheritance:
.. _api.esmvaltool.diag_scripts.emergent_constraints.init:

Auxiliary functions for emergent constraints scripts
====================================================

.. automodule:: esmvaltool.diag_scripts.emergent_constraints
.. _outputdata:

******
Output
******

ESMValTool automatically generates a new output directory with every run. The
location is determined by the output_dir option  in the config-user.yml file,
the recipe name, and the date and time, using the the format: YYYYMMDD_HHMMSS.

For instance, a typical output location would be:
output_directory/recipe_ocean_amoc_20190118_1027/

This is effectively produced by the combination:
output_dir/recipe_name_YYYYMMDD_HHMMSS/


This directory will contain 4 further subdirectories:

1. `Diagnostic output`_ (work): A place for any diagnostic script results that are not plots, e.g. files in NetCDF format (depends on the diagnostics).

2. `Plots`_: The location for all the plots, split by individual diagnostics and fields.

3. `Run`_: This directory includes all log files, a copy of the recipe, a summary of the resource usage, and the `settings.yml`_ interface files and temporary files created by the diagnostic scripts.

4. `Preprocessed datasets`_ (preproc): This directory contains all the preprocessed netcdfs data and the `metadata.yml`_ interface files. Note that by default this directory will be deleted after each run, because most users will only need the results from the diagnostic scripts.


Preprocessed datasets
=====================

The preprocessed datasets will be stored to the preproc/ directory.
Each variable in each diagnostic will have its own the `metadata.yml`_
interface files saved in the preproc directory.

If the option ``save_intermediary_cubes`` is set to ``true`` in the
config-user.yml file, then the intermediary cubes will also be saved here.
This option is set to false in the default ``config-user.yml`` file.

If the option ``remove_preproc_dir`` is set to ``true`` in the config-user.yml
file, then the preproc directory will be deleted after the run completes. This
option is set to true in the default  ``config-user.yml`` file.


Run
===

The log files in the run directory are automatically generated by ESMValTool
and create a record of the output messages produced by ESMValTool and they are
saved in the run directory. They can be helpful for debugging or monitoring the
job, but also allow a record of the job output to screen after the job has been
completed.

The run directory will also contain a copy of the recipe and the
`settings.yml`_ file, described below.
The run directory is also where the diagnostics are executed, and may also
contain several temporary files while diagnostics are running.

Diagnostic output
=================

The work/ directory will contain all files that are output at the diagnostic
stage. Ie, the model data is preprocessed by ESMValTool and stored in the
preproc/ directory. These files are opened by the diagnostic script, then some
processing is applied. Once the diagnostic level processing has been applied,
the results should be saved to the work directory.


Plots
=====

The plots directory is where diagnostics save their output figures. These
plots are saved in the format requested by the option `output_file_type` in the
config-user.yml file.


Settings.yml
============

The settings.yml file is automatically generated by ESMValCore. For each diagnostic,
a unique settings.yml file will be produced.

The settings.yml file passes several global level keys to diagnostic scripts.
This includes several flags from the config-user.yml file (such as
'write_netcdf', 'write_plots', etc...), several paths which are specific to the
diagnostic being run (such as 'plot_dir' and 'run_dir') and the location on
disk of the metadata.yml file (described below).

.. code-block:: yaml

    input_files:[[...]recipe_ocean_bgc_20190118_134855/preproc/diag_timeseries_scalars/mfo/metadata.yml]
    log_level: debug
    output_file_type: png
    plot_dir: [...]recipe_ocean_bgc_20190118_134855/plots/diag_timeseries_scalars/Scalar_timeseries
    profile_diagnostic: false
    recipe: recipe_ocean_bgc.yml
    run_dir: [...]recipe_ocean_bgc_20190118_134855/run/diag_timeseries_scalars/Scalar_timeseries
    script: Scalar_timeseries
    version: 2.0a1
    work_dir: [...]recipe_ocean_bgc_20190118_134855/work/diag_timeseries_scalars/Scalar_timeseries
    write_netcdf: true
    write_plots: true

The first item in the settings file will be a list of `Metadata.yml`_ files.
There is a metadata.yml file generated for each field in each diagnostic.


Metadata.yml
============

The metadata.yml files is automatically generated by ESMValTool. Along with the
settings.yml file, it passes all the paths, boolean flags, and additional
arguments that your diagnostic needs to know in order to run.

The metadata is loaded from cfg as a dictionary object in python diagnostics.

Here is an example metadata.yml file:

.. code-block:: yaml

  ?
    [...]/recipe_ocean_bgc_20190118_134855/preproc/diag_timeseries_scalars/mfo/CMIP5_HadGEM2-ES_Omon_historical_r1i1p1_TO0M_mfo_2002-2004.nc
    : cmor_table: CMIP5
    dataset: HadGEM2-ES
    diagnostic: diag_timeseries_scalars
    end_year: 2004
    ensemble: r1i1p1
    exp: historical
    field: TO0M
    filename: [...]recipe_ocean_bgc_20190118_134855/preproc/diag_timeseries_scalars/mfo/CMIP5_HadGEM2-ES_Omon_historical_r1i1p1_TO0M_mfo_2002-2004.nc
    frequency: mon
    institute: [INPE, MOHC]
    long_name: Sea Water Transport
    mip: Omon
    modeling_realm: [ocean]
    preprocessor: prep_timeseries_scalar
    project: CMIP5
    recipe_dataset_index: 0
    short_name: mfo
    standard_name: sea_water_transport_across_line
    start_year: 2002
    units: kg s-1
    variable_group: mfo


As you can see, this is effectively a dictionary with several items including
data paths, metadata and other information.

There are  several tools available in python which are built to read and parse
these files. The tools are avaialbe in the shared directory in the diagnostics
directory.
.. _config-user:

*************
Configuration
*************

The ``esmvaltool`` command is provided by the ESMValCore package, the
documentation on configuring ESMValCore can be found
:ref:`here <esmvalcore:config>`.
In particular, it is recommended to read the section on the
:ref:`User configuration file <esmvalcore:user configuration file>`
and the section on
:ref:`Finding data <esmvalcore:findingdata>`.

To install the default configuration file in the default location, run

 .. code:: bash

	  esmvaltool config get_config_user

Note that this file needs to be customized using the instructions above, so
the ``esmvaltool`` command can find the data on your system, before it can run
a recipe.

There is a lesson available in the
`ESMValTool tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`_
that describes how to personalize the configuration file. It can be found
`at this site <https://esmvalgroup.github.io/ESMValTool_Tutorial/03-configuration/index.html>`_.
.. _install:

************
Installation
************

.. note::
   ESMValTool now uses `mamba` instead of `conda` for the recommended installation.
   For more information about the change, have a look at :ref:`Move to Mamba<move-to-mamba>`.

ESMValTool 2.0 requires a Unix(-like) operating system and Python 3.7+.

The ESMValTool can be installed in multiple ways.

Recommended installation methods:

* On Linux, please install via the :ref:`mamba package manager<install_with_mamba>` (see https://anaconda.com);

* For MacOSX, please follow separate instructions for :ref:`installation on MacOSX<install_on_macosx>`.

Further options for installation are:

* :ref:`Installation with pip and mamba<install_with_pip>` (see https://pypi.org);

* :ref:`Deployment through a Docker container<install_with_docker>` (see https://www.docker.com);

* :ref:`Deployment through a Singularity container<install_with_singularity>` (see https://sylabs.io/guides/latest/user-guide/);

* :ref:`From the source code<install_from_source>` available at https://github.com/ESMValGroup/ESMValTool;

* :ref:`From pre-installed versions on HPC clusters<install_on_hpc>`.

The next sections will detail the procedure to install ESMValTool through each
of these methods.

Note that there is also a
`Tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`__
available with more explanations and guidance if the installation instructions
here are too concise.
See `common installation issues`_ if you run into trouble.

.. _install_with_mamba:

Mamba installation
==================

In order to install the `conda <https://docs.conda.io>`_ package, you will need
mamba pre-installed.

For a minimal mamba installation (recommended) go to
https://mamba.readthedocs.io/en/latest/installation.html.
It is recommended that you always use the latest version of mamba, as problems
have been reported when trying to use older versions.

First download the installation file for
`Linux <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh>`_
or
`MacOSX <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh>`_.
After downloading the installation file from one of the links above, execute it
by running (Linux example):

.. code-block:: bash

    bash Mambaforge-Linux-x86_64.sh

and follow the instructions on your screen.
Immediately update mamba after installing it:

.. code-block:: bash

    mamba update --name base mamba

You can check that mamba installed correctly by running

.. code-block:: bash

    which mamba

this should show the path to your mamba executable, e.g.
``~/mambaforge/bin/mamba``.

ESMValTool installation
-----------------------

Once you have installed the above prerequisites, you can install the entire
ESMValTool package by running:

.. code-block:: bash

    mamba create --name esmvaltool esmvaltool 'python=3.10'

Here ``mamba`` is the executable calling the mamba package manager to install
``esmvaltool``. The reason why we are also specifying ``'python=3.10'`` is that
it will make it easier for mamba to find a working combination of all required
packages, see `Mamba fails to solve the environment`_ in `common installation
issues`_ for an in-depth explanation. Python 3.7, 3.8 and 3.9 are also
supported, in case you prefer to work with an older version of Python.

This will create a new
`conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_
and install ESMValTool into it with a single command.


.. code-block:: bash

    conda activate esmvaltool

Of course it is also possible to choose a different name than ``esmvaltool`` for the environment.

The next step is to check that the installation works properly.
To do this, run the tool with:

.. code-block:: bash

    esmvaltool --help

If everything was installed properly, ESMValTool should have printed a help
message to the console.

.. note::

    Creating a new conda environment is often much faster and more reliable than
    trying to update an existing conda environment.

Julia installation
------------------

If you want to use the ESMValTool Julia functionality, you will also need to
install Julia. If you are just getting started, we suggest that you
come back to this step later when and if you need it.
To perform the Julia installation, make sure that your conda
environment is activated and then execute

.. code-block:: bash

    mamba install julia

.. _conda subpackages:

Installation of subpackages
---------------------------

The diagnostics bundled in ESMValTool are scripts in four different programming
languages: Python, NCL, R, and Julia.

There are three language specific packages available:

* ``esmvaltool-ncl``
* ``esmvaltool-python``
* ``esmvaltool-r``

The main ``esmvaltool`` package contains all four subpackages listed above. If
you only need to run a recipe with diagnostics in some of these languages, it is
possible to install only the dependencies needed to do just that. The diagnostic
script(s) used in each recipe, are documented in :ref:`recipes`. The extension
of the diagnostic script can be used to see in which language a diagnostic
script is written.

To install support for diagnostics written in Python and NCL into an existing
environment, run

.. code-block:: bash

    mamba install esmvaltool-python esmvaltool-ncl

Some of the CMORization scripts are written in Python, while others are written
in NCL. Therefore, both ``esmvaltool-python`` and ``esmvaltool-ncl`` need to be
installed in order to be able to run all CMORization scripts.

Note that the ESMValTool source code is contained in the ``esmvaltool-python``
package, so this package will always be installed as a dependency if you install
one or more of the packages for other languages.

There is also a lesson available in the
`ESMValTool tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`_
that describes the installation of the ESMValTool in more detail. It can be found
`here <https://esmvalgroup.github.io/ESMValTool_Tutorial/02-installation/index.html>`_.

.. _install_on_macosx:

Installation on MacOSX
======================

The Python diagnostics of the ESMValTool are supported on MacOSX, but Julia, NCL,
and R are not. If any of these are needed, deployment through a :ref:`Docker<install_with_docker>`
container is advised.

The ``esmvaltool-python`` diagnostics can be installed as follows:

First, ensure mamba is pre-installed (see `Mamba installation`_ for more details).

Create a new environment with the ``esmvaltool-python`` package:

.. code-block:: bash

    mamba create --name esmvaltool esmvaltool-python 'python=3.10'

Activate the new environment:

.. code-block:: bash

    conda activate esmvaltool

Confirm that the ESMValTool is working with:

.. code-block:: bash

    esmvaltool --help

Note that some recipes may depend on the OpenMP library, which does not
install via mamba on MacOSX. To install this library, run:

.. code-block:: bash

    brew install libomp

to install the library with Homebrew. In case you do not have Homebrew, follow
installation instructions `here <https://brew.sh/>`__.

.. _install_with_pip:

Pip installation
================

It is also possible to install ESMValTool from `PyPI <https://pypi.org/project/ESMValTool/>`_.
However, this requires first installing dependencies that are not available on PyPI in some other way.
By far the easiest way to install these dependencies is to use `mamba`.
For a minimal mamba installation (recommended) go to https://mamba.readthedocs.io/en/latest/installation.html.

After installing mamba, download
`the file with the list of dependencies <https://raw.githubusercontent.com/ESMValGroup/ESMValTool/main/environment.yml>`_:

.. code-block:: bash

    wget https://raw.githubusercontent.com/ESMValGroup/ESMValTool/main/environment.yml

and install these dependencies into a new conda environment with the command

.. code-block:: bash

    mamba env create --name esmvaltool --file environment.yml


Finally, activate the newly created environment

.. code-block:: bash

    conda activate esmvaltool

and install ESMValTool as well as any remaining Python dependencies with the command:

.. code-block:: bash

    pip install esmvaltool

If you would like to run Julia diagnostic scripts, you will also need to
install the Julia dependencies:

.. code-block:: bash

    esmvaltool install Julia

If you would like to run R diagnostic scripts, you will also need to install the R
dependencies:

.. code-block:: bash

    esmvaltool install R

.. _install_with_docker:

Docker installation
===================

ESMValTool is also provided through `DockerHub <https://hub.docker.com/u/esmvalgroup/>`_
in the form of docker containers.
See https://docs.docker.com for more information about docker containers and how to
run them.

You can get the latest release with

.. code-block:: bash

   docker pull esmvalgroup/esmvaltool:stable

If you want to use the current main branch, use

.. code-block:: bash

   docker pull esmvalgroup/esmvaltool:latest

To run a container using those images, use:

.. code-block:: bash

   docker run esmvalgroup/esmvaltool:stable --help

Note that the container does not see the data or environmental variables
available in the host by default. You can make data available with
``-v /path:/path/in/container`` and environmental variables with ``-e VARNAME``.

For example, the following command would run a recipe

.. code-block:: bash

   docker run -e HOME -v "$HOME":"$HOME" -v /data:/data esmvalgroup/esmvaltool:stable run examples/recipe_python.yml

with the environmental variable ``$HOME`` available inside the container and
the data in the directories ``$HOME`` and ``/data``, so these can be used to
find the configuration file, recipe, and data.

It might be useful to define a `bash alias
<https://opensource.com/article/19/7/bash-aliases>`_
or script to abbreviate the above command, for example

.. code-block:: bash

   alias esmvaltool="docker run -e HOME -v $HOME:$HOME -v /data:/data esmvalgroup/esmvaltool:stable"

would allow using the ``esmvaltool`` command without even noticing that the
tool is running inside a Docker container.

.. _install_with_singularity:

Singularity installation
========================

Docker is usually forbidden in clusters due to security reasons. However,
there is a more secure alternative to run containers that is usually available
on them: `Singularity <https://sylabs.io/guides/3.0/user-guide/quick_start.html>`_.

Singularity can use docker containers directly from DockerHub with the
following command

.. code-block:: bash

   singularity run docker://esmvalgroup/esmvaltool:stable run examples/recipe_python.yml

Note that the container does not see the data available in the host by default.
You can make host data available with ``-B /path:/path/in/container``.

It might be useful to define a `bash alias
<https://opensource.com/article/19/7/bash-aliases>`_
or script to abbreviate the above command, for example

.. code-block:: bash

   alias esmvaltool="singularity run -B $HOME:$HOME -B /data:/data docker://esmvalgroup/esmvaltool:stable"

would allow using the ``esmvaltool`` command without even noticing that the
tool is running inside a Singularity container.

Some clusters may not allow to connect to external services, in those cases
you can first create a singularity image locally:

.. code-block:: bash

   singularity build esmvaltool.sif docker://esmvalgroup/esmvaltool:stable

and then upload the image file ``esmvaltool.sif`` to the cluster.
To run the container using the image file ``esmvaltool.sif`` use:

.. code-block:: bash

   singularity run esmvaltool.sif run examples/recipe_python.yml

.. _install_from_source:

Install from source
===================

Installing the tool from source is recommended if you need the very latest
features or if you would like to contribute to its development.

Obtaining the source code
-------------------------

The ESMValTool source code is available on a public GitHub repository:
https://github.com/ESMValGroup/ESMValTool

The easiest way to obtain it is to clone the repository using git
(see https://git-scm.com/). To clone the public repository:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool

or

.. code-block:: bash

    git clone git@github.com:ESMValGroup/ESMValTool

if you prefer to connect to the repository over SSH.

The command above will create a folder called ``ESMValTool``
containing the source code of the tool in the current working directory.

.. note::
    Using SSH is much more convenient if you push to the repository regularly
    (recommended to back up your work), because then you do not need to type
    your password over and over again.
    See
    `this guide <https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account>`__
    for information on how to set it up if you have not done so yet.
    If you are developing ESMValTool on a shared compute cluster, you can set up
    `SSH agent forwarding <https://docs.github.com/en/free-pro-team@latest/developers/overview/using-ssh-agent-forwarding>`__
    to use your local SSH keys also from the remote machine.

It is also possible to work in one of the ESMValTool private repositories, e.g.:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValTool-private

GitHub also allows one to download the source code in as a ``tar.gz`` or ``zip``
file.
If you choose to use this option, download the compressed file and extract its
contents at the desired location.


Prerequisites
-------------

It is recommended to use mamba to manage ESMValTool dependencies.
For a minimal mamba installation go to https://mamba.readthedocs.io/en/latest/installation.html.
To simplify the installation process, an environment definition file is provided
in the repository (``environment.yml`` in the root folder).

.. attention::
    Some systems provide a preinstalled version of conda (e.g., via the module environment).
    However, several users reported problems when installing NCL with such versions. It is
    therefore preferable to use a local, fully user-controlled mamba installation.
    Using an older version of mamba can also be a source of problems, so if you have mamba
    installed already, make sure it is up to date by running ``mamba update -n base mamba``.

To enable the ``mamba`` command, please source the appropriate configuration file
from your ``~/.bashrc``  file:

.. code-block:: bash

    source <prefix>/etc/profile.d/conda.sh

or ``~/.cshrc``/``~/.tcshrc`` file:

.. code-block:: bash

    source <prefix>/etc/profile.d/conda.csh

where ``<prefix>`` is the install location of your anaconda or miniconda
(e.g. ``/home/$USER/mambaforge``, ``/home/$USER/anaconda3`` or ``/home/$USER/miniconda3``).


.. note::
    Note that during the installation, mamba will ask you
    if you want the installation to be automatically sourced from your
    ``.bashrc`` or ``.bash-profile`` files; if you answered yes, then mamba
    will write bash directives to those files and every time you get to your
    shell, you will automatically be inside conda's ``(base)`` environment.
    To deactivate this feature, look for the ``# >>> conda initialize >>>``
    code block in your ``.bashrc`` or ``.bash-profile`` and comment the whole block out.


The ESMValTool conda environment file can also be used as a requirements list
for those cases in which a mamba installation is not possible or advisable.
From now on, we will assume that the installation is going to be done through
mamba.

Ideally, you should create a separate conda environment for ESMValTool, so it is
independent from any other Python tools present in the system.

Note that it is advisable to update mamba to the latest version before
installing ESMValTool, using the command (as mentioned above)

.. code-block:: bash

    mamba update --name base mamba

To create an environment, go to the directory containing the ESMValTool source
code (called ``ESMValTool`` if you did not choose a different name)

.. code-block:: bash

    cd ESMValTool

and (when on Linux) create a new environment called ``esmvaltool``
containing just Python with the command

.. code-block:: bash

    mamba create --name esmvaltool 'python=3.10'

if needed, older versions of Python can also be selected.
Next, install many of the required dependencies, including the ESMValCore package
and Python, R, and NCL interpreters, into this environment by running

.. code-block:: bash

    mamba env update --name esmvaltool --file environment.yml

**MacOSX note:** ESMValTool functionalities in Julia, NCL, and R are not
supported on MacOSX, due to conflicts in the conda environment. To install a
conda environment on MacOSX, use the dedicated environment file:

.. code-block:: bash

    mamba env create --name esmvaltool --file environment_osx.yml

The environment is called ``esmvaltool`` by default, but it is possible to use
the option ``--name SOME_ENVIRONMENT_NAME`` to define a custom name. You should
then activate the environment using the command:

.. code-block:: bash

    conda activate esmvaltool

It is also possible to update an existing environment from the environment
file. This may be useful when updating an older installation of ESMValTool:

.. code-block:: bash

    mamba env update --name esmvaltool --file environment.yml

but if you run into trouble, please try creating a new environment.

.. attention::
    From now on, we assume that the conda environment for ESMValTool is
    activated.

Software installation
---------------------

Once all prerequisites are fulfilled, ESMValTool can be installed by running
the following commands in the directory containing the ESMValTool source code
(called ``ESMValTool`` if you did not choose a different name):

.. code-block:: bash

    pip install --editable '.[develop]'

Using the ``--editable`` flag will cause the installer to create a symbolic link
from the installation location to your source code, so any changes you make to
the source code will immediately be available in the installed version of the
tool.
This command will also install extra development dependencies needed for
building the documentation, running the unit tests, etc.

If you would like to run Julia diagnostic scripts, you will need to
install the ESMValTool Julia dependencies:

.. code-block:: bash

    esmvaltool install Julia

If you would like to run R diagnostic scripts, you will also need to install the R
dependencies. Install the R dependency packages:

.. code-block:: bash

    esmvaltool install R

The next step is to check that the installation works properly.
To do this, run the tool with:

.. code-block:: bash

    esmvaltool --help

If everything was installed properly, ESMValTool should have printed a
help message to the console.

**MacOSX note:** some recipes may depend on the OpenMP library, which does not
install via mamba on MacOSX. Instead run

.. code-block:: bash

    brew install libomp

to install the library with Homebrew. In case you do not have Homebrew, follow
installation instructions `here <https://brew.sh/>`__.

For a more complete installation verification, run the automated tests and
confirm that no errors are reported:

.. code-block:: bash

    pytest -m "not installation"

or if you want to run the full test suite remove the ``-m "not installation"`` flag;
also if you want to run the tests on multiple threads, making the run faster, use
the `-n N` flag where N is the number of available threads e.g:

.. code-block:: bash

    pytest -n 4


.. _esmvalcore-development-installation:

Using the development version of the ESMValCore package
-------------------------------------------------------

If you need the latest developments of the ESMValCore package, you
can install it from source into the same conda environment.

.. attention::
    The recipes and diagnostics in the ESMValTool repository are compatible
    with the latest released version of the ESMValCore.
    Using the development version of the ESMValCore package is only recommended
    if you are planning to develop new features for the ESMValCore, e.g.
    you want to implement a new preprocessor function.

First follow all steps above.
Next, go to the place where you would like to keep the source code and clone the
ESMValCore github repository:

.. code-block:: bash

    git clone https://github.com/ESMValGroup/ESMValCore

or

.. code-block:: bash

    git clone git@github.com:ESMValGroup/ESMValCore

The command above will create a folder called ``ESMValCore``
containing the source code of the tool in the current working directory.

Go into the folder you just downloaded

.. code-block:: bash

    cd ESMValCore

and then install ESMValCore in development mode

.. code-block:: bash

    pip install --editable '.[develop]'

To check that the installation was successful, run

.. code-block:: bash

    python -c 'import esmvalcore; print(esmvalcore.__path__[0])'

this should show the directory of the source code that you just downloaded.

If the command above shows a directory inside your conda environment instead,
e.g. ``~/mamba/envs/esmvaltool/lib/python3.8/site-packages/esmvalcore``, you
may need to manually remove that directory and run
```pip install -e '.[develop]'``
again.

.. _install_on_hpc:

Pre-installed versions on HPC clusters
======================================

The ESMValTool is also available on the HPC clusters CEDA-JASMIN and DKRZ-MISTRAL and there will be no need
to install it yourself if you are just running diagnostics:

 - CEDA-JASMIN: `esmvaltool` is available on the scientific compute nodes (`sciX.jasmin.ac.uk` where
   `X = 1, 2, 3, 4, 5`) after login and module loading via `module load esmvaltool`; see the helper page at
   `CEDA <https://help.jasmin.ac.uk/article/4955-community-software-esmvaltool>`__ ;
 - DKRZ-Mistral: `esmvaltool` is available on login nodes (`mistral.dkrz.de`) and pre- and post-processing
   nodes (`mistralpp.dkrz.de`) after login and module loading via `module load esmvaltool`; the command
   `module help esmvaltool` provides some information about the module.

.. _common installation issues:

Common installation problems and their solutions
================================================

Mamba fails to solve the environment
------------------------------------
If you see the text ``Solving environment:`` with the characters ``-\|/`` rotating
behind it for more than 10 minutes, mamba may be having problems finding a
working combination of versions of the packages that the ESMValTool depends on.
Because the ESMValTool is a community tool, there is no strict selection of
which tools can be used and installing the ESMValTool requires installing almost
any package that is available for processing climate data.
To help mamba solve the environment, you can try the following.

Always use the latest version of mamba, as problems have been reported by people
using older versions, to update, run:

.. code-block:: bash

    mamba update --name base mamba

Usually mamba is much better at solving new environments than updating older
environments, so it is often a good idea to create a new environment if updating
does not work.

It can help mamba if you let it know what version of certain packages you want,
for example by running

.. code-block:: bash

    mamba create -n esmvaltool esmvaltool 'python=3.10'

you ask for Python 3.10 specifically and that makes it much easier for mamba to
solve the environment, because now it can ignore any packages that were built
for other Python versions. Note that, since the esmvaltool package is built
with Python>=3.7, asking for an older Python version, e.g. `python=3.6`, in
this way, it will result in installation failure.

Problems with proxies
---------------------
If you are installing ESMValTool from source from behind a proxy that does not
trust the usual PyPI URLs you can declare them with the option
``--trusted-host``, e.g.

.. code-block:: bash

    pip install --trusted-host=pypi.python.org --trusted-host=pypi.org --trusted-host=files.pythonhosted.org -e .[develop]

If R packages fail to download, you might be able to solve this by
setting the environment variable ``http_proxy`` to the correct value, e.g.
in bash:

.. code-block:: bash

    export http_proxy=http://user:pass@proxy_server:port

the username and password can be omitted if they are not required. See e.g.
`here <https://support.rstudio.com/hc/en-us/articles/200488488-Configuring-R-to-Use-an-HTTP-or-HTTPS-Proxy>`__
for more information.

Anaconda servers connection issues
----------------------------------
HTTP connection errors (of e.g. type 404) to the Anaconda servers are rather common, and usually a retry
will solve the problem.

Installation of R packages fails
--------------------------------
Problems have been reported if the ``R`` interpreter was made available
through the ``module load`` command in addition to installation from mamba.
If your ESMValTool conda environment is called ``esmvaltool`` and you want to
use the R interpreter installed from mamba, the path to the R interpreter should
end with ``mamba/envs/esmvaltool/bin/R`` or ``conda/envs/esmvaltool/bin/R``.
When the conda environment for ESMValTool is activated, you can check which R
interpreter is used by running

.. code-block:: bash

    which R

The Modules package is often used by system administrators to make software
available to users of scientific compute clusters.
To list any currently loaded modules run ``module list``, run ``module help``
or ``man module`` for more information about the Modules package.

Problems when using ssh
-----------------------
If you log in to a cluster or other device via SSH and your origin
machine sends the ``locale`` environment via the SSH connection,
make sure the environment is set correctly, specifically ``LANG`` and
``LC_ALL`` are set correctly (for GB English UTF-8 encoding these
variables must be set to ``en_GB.UTF-8``; you can set them by adding
``export LANG=en_GB.UTF-8`` and ``export LC_ALL=en_GB.UTF-8``) in your
origin or login machines’ ``.profile``.

Problems when updating the conda environment
--------------------------------------------
Usually mamba is much better at solving new environments than updating older
environments, so it is often a good idea to create a new environment if updating
does not work. See also `Mamba fails to solve the environment`_.

Do not run ``mamba update --update-all`` in the ``esmvaltool``
environment since that will update some packages that are pinned to
specific versions for the correct functionality of the tool.


.. _move-to-mamba:

Move to Mamba
=============

Mamba is a much faster alternative to `conda`, and environment creation and updating
benefits from the use of a much faster (C++ backend) dependency solver; tests have been performed
to verify the integrity of the `esmvaltool` environment built with `mamba`, and we are
now confident that the change will not affect the way ESMValTool is installed and run, whether it be on a Linux or OSX platform.
From the user's perspective, it is a straightforward use change: the CLI (command line
interface) of `mamba` is identical to `conda`: any command that was run with `conda` before
will now be run with `mamba` instead, keeping all the other command line arguments and
flags as they were before. The only place where `conda` should not be replaced with `mamba`
at command line level is at the environment activation point: `conda activate` will still
have to be used.
.. _running:

*******
Running
*******

ESMValTool is mostly used as a command line tool. Whenever your
conda environment for ESMValTool is active, you can just run the command
``esmvaltool``. See
:ref:`running esmvaltool <esmvalcore:running>`
in the ESMValCore documentation for a short introduction.

Running a recipe
================

An
`example recipe <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_python.yml>`_
is available in the ESMValTool installation folder as
``examples/recipe_python.yml``.

This recipe finds data from BCC-ESM1 and CanESM2 and creates two plot types:

- a global map plot that shows the monthly mean 2m surface air temperature in
  January 2000.
- a time series plot that shows the globally averaged annual mean 2m surface
  air temperature and compares it to the one in Amsterdam.

You can download the recipe from :download:`here 
<../../../../esmvaltool/recipes/examples/recipe_python.yml>` and save it in
your project directory as (e.g.) ``recipe_python.yml`` and then run ESMValTool
with

.. code:: bash

   esmvaltool run --offline=False recipe_python.yml

The ``--offline=False`` option tells ESMValTool to automatically download
missing data. This might require additional :ref:`configuration
<esmvalcore:config-esgf>` to work. The data only needs to be downloaded once,
every following run will re-use previously downloaded data.

ESMValTool will also find recipes that are stored in its installation directory.
A copy of the example recipe is shipped with ESMValTool as:
``/path/to/installation/esmvaltool/recipes/examples/recipe_python.yml``.
Thus, the following also works:

.. code:: bash

	esmvaltool run examples/recipe_python.yml

Note that this command does not automatically download missing data. The
required data should thus be located in the directories specified in your user
configuration file.  Recall that the chapter :ref:`Configuring ESMValTool
<config-user>` provides an explanation of how to create your own
config-user.yml file.

To get help on additional commands, please use

.. code:: bash

	esmvaltool --help

It is also possible to get help on specific commands, e.g.

.. code:: bash

	esmvaltool run --help

will display the help message with all options for the ``run`` command.

There is a step-by-step description available in the
`ESMValTool tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`_
on how to run your first recipe. It can be found
`here <https://esmvalgroup.github.io/ESMValTool_Tutorial/04-recipe/index.html>`_.

Available diagnostics and metrics
=================================

See Section :doc:`Recipes <../recipes/index>` for a description of all
available recipes.

To see a list of installed recipes run

.. code:: bash

	esmvaltool recipes list

Running multiple recipes
========================

It is possible to run more than one recipe in one go: currently this relies on the user
having access to a HPC that has ``rose`` and ``cylc`` installed since the procedure involves
installing and submitting a Rose suite. the utility that allows you to do this is
``esmvaltool/utils/rose-cylc/esmvt_rose_wrapper.py``.

Base suite:
-----------
The base suite to run esmvaltool via rose-cylc is `u-bd684`; you can find
this suite in the Met Office Rose repository at:

https://code.metoffice.gov.uk/svn/roses-u/b/d/6/8/4/trunk/

When ``rose`` will be working with python3.x, this location will become
default and the pipeline will aceess it independently of user, unless, of
course the user will specify ``-s $SUITE_LOCATION``; until then the user needs
to grab a copy of it in ``$HOME`` or specify the default location via ``-s`` option.

Environment:
------------
We will move to a unified and centrally-installed esmvaltool environment;
until then, the user will have to alter the env_setup script:

``u-bd684/app/esmvaltool/env_setup``

with the correct pointers to esmvaltool installation, if desired.

To be able to submit to cylc, you need to have the `/metomi/` suite in path
AND use a `python2.7` environment. Use the Jasmin-example below for guidance.

Jasmin-example:
---------------
This shows how to interact with rose-cylc and run esmvaltool under cylc
using this script:

.. code:: bash

   export PATH=/apps/contrib/metomi/bin:$PATH
   export PATH=/home/users/valeriu/miniconda2/bin:$PATH
   mkdir esmvaltool_rose
   cd esmvaltool_rose
   cp ESMValTool/esmvaltool/utils/rose-cylc/esmvt_rose_wrapper.py .
   svn checkout https://code.metoffice.gov.uk/svn/roses-u/b/d/6/8/4/trunk/ ~/u-bd684
   [enter Met Office password]
   [configure ~/u-bd684/rose_suite.conf]
   [configure ~/u-bd684/app/esmvaltool/env_setup]
   python esmvt_rose_wrapper.py -c config-user.yml \
   -r recipe_autoassess_stratosphere.yml recipe_OceanPhysics.yml \
   -d $HOME/esmvaltool_rose
   rose suite-run u-bd684

Note that you need to pass FULL PATHS to cylc, no `.` or `..` because all
operations are done remotely on different nodes.

A practical actual example of running the tool can be found on JASMIN:
``/home/users/valeriu/esmvaltool_rose``.
There you will find the run shell: ``run_example``, as well as an example
how to set the configuration file. If you don't have Met Office credentials,
a copy of `u-bd684` is always located in ``/home/users/valeriu/roses/u-bd684`` on Jasmin.
Getting started
***************

.. toctree::
   :maxdepth: 1

    Installation <installation>
    Configuration <configuration>
    Running <running>
    Output <output>
