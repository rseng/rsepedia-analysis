# ESMValCore package

[![Documentation Status](https://readthedocs.org/projects/esmvaltool/badge/?version=latest)](https://esmvaltool.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3387139.svg)](https://doi.org/10.5281/zenodo.3387139)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ESMValGroup?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![CircleCI](https://circleci.com/gh/ESMValGroup/ESMValCore/tree/main.svg?style=svg)](https://circleci.com/gh/ESMValGroup/ESMValCore/tree/main)
[![codecov](https://codecov.io/gh/ESMValGroup/ESMValCore/branch/main/graph/badge.svg?token=wQnDzguwq6)](https://codecov.io/gh/ESMValGroup/ESMValCore)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/5d496dea9ef64ec68e448a6df5a65783)](https://www.codacy.com/gh/ESMValGroup/ESMValCore?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ESMValGroup/ESMValCore&amp;utm_campaign=Badge_Grade)
[![Docker Build Status](https://img.shields.io/docker/cloud/build/esmvalgroup/esmvalcore)](https://hub.docker.com/r/esmvalgroup/esmvalcore/)
[![Anaconda-Server Badge](https://anaconda.org/esmvalgroup/esmvalcore/badges/installer/conda.svg)](https://conda.anaconda.org/esmvalgroup)

![esmvaltoollogo](https://github.com/ESMValGroup/ESMValCore/blob/main/doc/figures/ESMValTool-logo-2.png)

ESMValCore: core functionalities for the ESMValTool, a community diagnostic
and performance metrics tool for routine evaluation of Earth System Models
in the Climate Model Intercomparison Project (CMIP).

# Getting started

Please have a look at the
[documentation](https://docs.esmvaltool.org/projects/esmvalcore/en/latest/quickstart/install.html)
to get started.

## Using the ESMValCore package to run recipes

The ESMValCore package provides the `esmvaltool` command, which can be used to run
[recipes](https://docs.esmvaltool.org/projects/esmvalcore/en/latest/recipe/overview.html)
for working with CMIP-like data.
A large collection of ready to use
[recipes and diagnostics](https://docs.esmvaltool.org/en/latest/recipes/index.html)
is provided by the
[ESMValTool](https://github.com/ESMValGroup/ESMValTool)
package.

## Using ESMValCore as a Python library

The ESMValCore package provides various functions for:

-   Finding data in a directory structure typically used for CMIP data.

-   Reading CMIP/CMOR tables and using those to check model and observational data.

-   ESMValTool preprocessor functions based on
    [iris](https://scitools-iris.readthedocs.io) for e.g. regridding,
    vertical interpolation, statistics, correcting (meta)data errors, extracting
    a time range, etcetera.

read all about it in the
[API documentation](https://docs.esmvaltool.org/projects/esmvalcore/en/latest/api/esmvalcore.html).

## Getting help

The easiest way to get help if you cannot find the answer in the documentation
on [readthedocs](https://docs.esmvaltool.org), is to open an
[issue on GitHub](https://github.com/ESMValGroup/ESMValCore/issues).

## Contributing

Contributions are very welcome, please read our
[contribution guidelines](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html)
to get started.
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

Please read our [contribution guidelines](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html).
<!--
    Thank you for contributing to our project!

    Please do not delete this text completely, but read the text below and keep
    items that seem relevant. If in doubt, just keep everything and add your
    own text at the top, a reviewer will update the checklist for you.

-->

## Description

<!--
    Please describe your changes here, especially focusing on why this pull
    request makes ESMValCore better and what problem it solves.

    Before you start, please read our contribution guidelines: https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html

    Please fill in the GitHub issue that is closed by this pull request,
    e.g. Closes #1903
-->

Closes #issue_number

Link to documentation:

***

## [Before you get started](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#getting-started)

- [ ] [☝ Create an issue](https://github.com/ESMValGroup/ESMValCore/issues) to discuss what you are going to do

## [Checklist](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#checklist-for-pull-requests)

It is the responsibility of the author to make sure the pull request is ready to review. The icons indicate whether the item will be subject to the [🛠 Technical][1] or [🧪 Scientific][2] review.

<!-- The next two lines turn the 🛠 and 🧪 below into hyperlinks -->
[1]: https://docs.esmvaltool.org/en/latest/community/review.html#technical-review
[2]: https://docs.esmvaltool.org/en/latest/community/review.html#scientific-review

- [ ] [🧪][2] The new functionality is [relevant and scientifically sound](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#scientific-relevance)
- [ ] [🛠][1] This pull request has a [descriptive title and labels](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#pull-request-title-and-label)
- [ ] [🛠][1] Code is written according to the [code quality guidelines](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#code-quality)
- [ ] [🧪][2] and [🛠][1] [Documentation](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#documentation) is available
- [ ] [🛠][1] [Unit tests](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#tests) have been added
- [ ] [🛠][1] Changes are [backward compatible](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#backward-compatibility)
- [ ] [🛠][1] Any changed [dependencies have been added or removed](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#dependencies) correctly
- [ ] [🛠][1] The [list of authors](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#list-of-authors) is up to date
- [ ] [🛠][1] All [checks below this pull request](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#pull-request-checks) were successful

***

To help with the number pull requests:

- 🙏 We kindly ask you to [review](https://docs.esmvaltool.org/en/latest/community/review.html#review-of-pull-requests) two other [open pull requests](https://github.com/ESMValGroup/ESMValCore/pulls) in this repository
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
name: Dataset issue report
about: Create a report about a problematic dataset
title: 'Dataset problem: '
labels: ''
assignees: valeriupredoi, zklaus

---

**Describe the dataset issue**
A clear and concise description of what the problem with the data is is. Here are some guidelines that you can use
right out the box by filling in the blanks (note that, if needed we will assist you with writing a fix; also note that
if the problem is with a CMIP6 dataset, we will alert the ESGF and/or the model developers for them to fix the problem
in the long run; fixes for CMIP5 by the model developers are no longer possible):

- Data file has been changed by you in any way (if you answer yes the issue will be void and closed, we are not
supporting custom-made datasets that present problems, it is your resposability to fix those problems):
- Project (CMIP5/CMIP6/OBS/obs4MIPs/etc):
- Full dataset description (fill out as many as you know, please):
  - dataset:
  - experiment:
  - mip:
  - grid:
  - type:
  - version:
  - tier:
  - variable used:
  - data version:
- Problems encountered (please describe what the actual problems are: e.g. wrong standard name, issue with dimensions):
- Pointer to existing copy of data on ESGF node (it would be very useful if you could provide a physical
fullpath to the file(s) that are causing the problem, e.g. on CEDA Jasmin or DKRZ):
- Other notes and mentions:

**Assign the right people**
If you are already a member of the ESMValTool GitHub project, please assign Valeriu Predoi (valeriupredoi) and
Klaus Zimmermann (zklaus) to the issue. They will then check the issue raised and propagate it further to the
data model developers.

**Please attach**
  - The recipe that you are trying to run, you can find a copy in the `run` directory in the output directory
  - The `main_log_debug.txt` file, this can also be found in the `run` directory in the output directory
  - Attach a text file with the output of `ncdump -h /path/to/file.nc`

Cheers :beer:
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
# cmip3-cmor-tables
Text Tables for CMOR1 to create a CMIP3 dataset
# cmip6-cmor-tables

Note:
----
Only `CMOR` and `PrePARE` should use the `JSON` file of concatenated `CMIP6` Control Vocabulary [CMIP6_CV.json](https://github.com/PCMDI/cmip6-cmor-tables/blob/master/Tables/CMIP6_CV.json) found in this repository.  

All other software should get the CV's directly from the [WCRP](https://github.com/WCRP-CMIP/CMIP6_CVs) repository. 

CMIP6 Data Request can be found in the following [mip tables](http://clipc-services.ceda.ac.uk/dreq/index/miptable.html).
.. _contributing:

Contributions are very welcome
==============================

We greatly value contributions of any kind.
Contributions could include, but are not limited to documentation improvements, bug reports, new or improved code, scientific and technical code reviews, infrastructure improvements, mailing list and chat participation, community help/building, education and outreach.
We value the time you invest in contributing and strive to make the process as easy as possible.
If you have suggestions for improving the process of contributing, please do not hesitate to propose them.

If you have a bug or other issue to report or just need help, please open an issue on the
`issues tab on the ESMValCore github repository <https://github.com/ESMValGroup/ESMValCore/issues>`__.

If you would like to contribute a new
:ref:`preprocessor function <preprocessor_function>`,
:ref:`derived variable <derivation>`, :ref:`fix for a dataset <fixing_data>`, or
another new feature, please discuss your idea with the development team before
getting started, to avoid double work and/or disappointment later.
A good way to do this is to open an
`issue <https://github.com/ESMValGroup/ESMValCore/issues>`_ on GitHub.

Getting started
---------------

See :ref:`installation-from-source` for instructions on how to set up a development
installation.

New development should preferably be done in the
`ESMValCore <https://github.com/ESMValGroup/ESMValCore>`__
GitHub repository.
The default git branch is ``main``.
Use this branch to create a new feature branch from and make a pull request
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
It’s also easier to get help from other developers if your code is visible in a
pull request.

:ref:`Make small pull requests <easy_review>`, the ideal pull requests changes
just a few files and adds/changes no more than 100 lines of production code.
The amount of test code added can be more extensive, but changes to existing
test code should be made sparingly.

Design considerations
~~~~~~~~~~~~~~~~~~~~~

When making changes, try to respect the current structure of the program.
If you need to make major changes to the structure of program to add a feature,
chances are that you have either not come up with the
most optimal design or the feature is not a very good fit for the tool.
Discuss your feature with the `@ESMValGroup/esmvaltool-coreteam`_ in an issue_
to find a solution.

Please keep the following considerations in mind when programming:

- Changes should preferably be :ref:`backward compatible <backward_compatibility>`.
- Apply changes gradually and change no more than a few files in a single pull
  request, but do make sure every pull request in itself brings a meaningful
  improvement.
  This reduces the risk of breaking existing functionality and making
  :ref:`backward incompatible <backward_compatibility>` changes, because it
  helps you as well as the reviewers of your pull request to better understand
  what exactly is being changed.
- :ref:`preprocessor_functions` are Python functions (and not classes) so they
  are easy to understand and implement for scientific contributors.
- No additional CMOR checks should be implemented inside preprocessor functions.
  The input cube is fixed and confirmed to follow the specification in
  `esmvalcore/cmor/tables <https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/cmor/tables>`__
  before applying any other preprocessor functions.
  This design helps to keep the preprocessor functions and diagnostics scripts
  that use the preprocessed data from the tool simple and reliable.
  See :ref:`cmor_table_configuration` for the mapping from ``project`` in the
  recipe to the relevant CMOR table.
- The ESMValCore package is based on :ref:`iris <iris_docs>`.
  Preprocessor functions should preferably be small and just call the relevant
  iris code.
  Code that is more involved and more broadly applicable than just in the
  ESMValCore, should be implemented in iris instead.
- Any settings in the recipe that can be checked before loading the data should
  be checked at the :ref:`task creation stage <Diagnostics>`.
  This avoids that users run a recipe for several hours before finding out they
  made a mistake in the recipe.
  No data should be processed or files written while creating the tasks.
- CMOR checks should provide a good balance between reliability of the tool
  and ease of use.
  Several :ref:`levels of strictness of the checks <cmor_check_strictness>`
  are available to facilitate this.
- Keep your code short and simple: we would like to make contributing as easy as
  possible.
  For example, avoid implementing complicated class inheritance structures and
  `boilerplate <https://stackoverflow.com/questions/3992199/what-is-boilerplate-code>`__
  code.
- If you find yourself copy-pasting a piece of code and making minor changes
  to every copy, instead put the repeated bit of code in a function that you can
  re-use, and provide the changed bits as function arguments.
- Be careful when changing existing unit tests to make your new feature work.
  You might be breaking existing features if you have to change existing tests.

Finally, if you would like to improve the design of the tool, discuss your plans
with the `@ESMValGroup/esmvaltool-coreteam`_ to make sure you understand the
current functionality and you all agree on the new design.

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
:ref:`pull request reviews <esmvaltool:reviewing>` to ensure all code and
documentation contributions are of good quality.
The icons indicate whether the item will be checked during the
:ref:`🛠 Technical review <technical_review>` or
:ref:`🧪 Scientific review <scientific_review>`.

- 🧪 The new functionality is :ref:`relevant and scientifically sound<scientific_relevance>`
- 🛠 :ref:`The pull request has a descriptive title and labels <descriptive_pr_title>`
- 🛠 Code is written according to the :ref:`code quality guidelines <code_quality>`
- 🧪 and 🛠 Documentation_ is available
- 🛠 Unit tests_ have been added
- 🛠 Changes are :ref:`backward compatible <backward_compatibility>`
- 🛠 Changed :ref:`dependencies have been added or removed correctly <dependencies>`
- 🛠 The :ref:`list of authors <authors>` is up to date
- 🛠 The :ref:`checks shown below the pull request <pull_request_checks>` are successful

.. _scientific_relevance:

Scientific relevance
--------------------

The proposed changes should be relevant for the larger scientific community.
The implementation of new features should be scientifically sound; e.g.
the formulas used in new preprocesssor functions should be accompanied by the
relevant references and checked for correctness by the scientific reviewer.
The `CF Conventions <https://cfconventions.org/>`_ as well as additional
standards imposed by `CMIP <https://www.wcrp-climate.org/wgcm-cmip>`_ should be
followed whenever possible.

.. _descriptive_pr_title:

Pull request title and label
----------------------------

The title of a pull request should clearly describe what the pull request changes.
If you need more text to describe what the pull request does, please add it in
the description.
`Add one or more labels <https://docs.github.com/en/github/managing-your-work-on-github/managing-labels#applying-labels-to-issues-and-pull-requests>`__
to your pull request to indicate the type of change.
At least one of the following
`labels <https://github.com/ESMValGroup/ESMValCore/labels>`__ should be used:
`bug`, `deprecated feature`, `fix for dataset`, `preprocessor`, `cmor`, `api`,
`testing`, `documentation` or `enhancement`.

The titles and labels of pull requests are used to compile the :ref:`changelog`,
therefore it is important that they are easy to understand for people who are
not familiar with the code or people in the project.
Descriptive pull request titles also makes it easier to find back what was
changed when, which is useful in case a bug was introduced.

.. _code_quality:

Code quality
------------

To increase the readability and maintainability or the ESMValCore source
code, we aim to adhere to best practices and coding standards.

We include checks for Python and yaml files, most of which are described in more
detail in the sections below.
This includes checks for invalid syntax and formatting errors.
:ref:`esmvaltool:pre-commit` is a handy tool that can run all of these checks
automatically just before you commit your code.
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
`readthedocs <https://docs.esmvaltool.org>`_.

To check if your code adheres to the standard, go to the directory where
the repository is cloned, e.g. ``cd ESMValCore``, and run `prospector <http://prospector.landscape.io/>`_

::

   prospector esmvalcore/preprocessor/_regrid.py

In addition to prospector, we use `flake8 <https://flake8.pycqa.org/en/latest/>`_
to automatically check for bugs and formatting mistakes and
`mypy <https://mypy.readthedocs.io>`_ for checking that
`type hints <https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html>`_ are
correct.
Note that `type hints`_ are completely optional, but if you do choose to add
them, they should be correct.

When you make a pull request, adherence to the Python development best practices
is checked in two ways:

#. As part of the unit tests, flake8_ and mypy_ are run by
   `CircleCI <https://app.circleci.com/pipelines/github/ESMValGroup/ESMValCore>`_,
   see the section on Tests_ for more information.
#. `Codacy <https://app.codacy.com/gh/ESMValGroup/ESMValCore/pullRequests>`_
   is a service that runs prospector (and other code quality tools) on changed
   files and reports the results.
   Click the 'Details' link behind the Codacy check entry and then click
   'View more details on Codacy Production' to see the results of the static
   code analysis done by Codacy_.
   If you need to log in, you can do so using your GitHub account.

The automatic code quality checks by prospector are really helpful to improve
the quality of your code, but they are not flawless.
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

YAML
~~~~

Please use `yamllint <https://yamllint.readthedocs.io>`_ to check that your
YAML files do not contain mistakes.
``yamllint`` checks for valid syntax, common mistakes like key repetition and
cosmetic problems such as line length, trailing spaces, wrong indentation, etc.

Any text file
~~~~~~~~~~~~~

A generic tool to check for common spelling mistakes is
`codespell <https://pypi.org/project/codespell/>`__.

.. _documentation:

Documentation
-------------

The documentation lives on `docs.esmvaltool.org <https://docs.esmvaltool.org>`_.

Adding documentation
~~~~~~~~~~~~~~~~~~~~

The documentation is built by readthedocs_ using `Sphinx <https://www.sphinx-doc.org>`_.
There are two main ways of adding documentation:

#. As written text in the directory
   `doc <https://github.com/ESMValGroup/ESMValCore/tree/main/doc/>`__.
   When writing
   `reStructuredText <https://www.sphinx-doc.org/en/main/usage/restructuredtext/basics.html>`_
   (``.rst``) files, please try to limit the line length to 80 characters and
   always start a sentence on a new line.
   This makes it easier to review changes to documentation on GitHub.

#. As docstrings or comments in code.
   For Python code, only the
   `docstrings <https://www.python.org/dev/peps/pep-0257/>`__
   of Python modules, classes, and functions
   that are mentioned in
   `doc/api <https://github.com/ESMValGroup/ESMValCore/tree/main/doc/api>`__
   are used to generate the online documentation.
   This results in the :ref:`api`.
   The standard document with best practices on writing docstrings is
   `PEP257 <https://www.python.org/dev/peps/pep-0257/>`__.
   For the API documentation, we make use of
   `numpy style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`__.

What should be documented
~~~~~~~~~~~~~~~~~~~~~~~~~

Functionality that is visible to users should be documented.
Any documentation that is visible on readthedocs_ should be well
written and adhere to the standards for documentation.
Examples of this include:

- The :ref:`recipe <recipe_overview>`
- Preprocessor :ref:`functions <preprocessor_functions>` and their
  :ref:`use from the recipe <preprocessor>`
- :ref:`Configuration options <config>`
- :ref:`Installation <install>`
- :ref:`Output files <outputdata>`
- :ref:`Command line interface <running>`
- :ref:`Diagnostic script interfaces <interfaces>`
- :ref:`The experimental Python interface <experimental_api>`

Note that:

- For functions that compute scientific results, comments with references to
  papers and/or other resources as well as formula numbers should be included.
- When making changes to/introducing a new preprocessor function, also update the
  :ref:`preprocessor documentation <preprocessor>`.
- There is no need to write complete numpy style documentation for functions that
  are not visible in the :ref:`api` chapter on readthedocs.
  However, adding a docstring describing what a function does is always a good
  idea.
  For short functions, a one-line docstring is usually sufficient, but more
  complex functions might require slightly more extensive documentation.

When reviewing a pull request, always check that documentation is easy to
understand and available in all expected places.

How to build and view the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whenever you make a pull request or push new commits to an existing pull
request, readthedocs will automatically build the documentation.
The link to the documentation will be shown in the list of checks below your
pull request.
Click 'Details' behind the check ``docs/readthedocs.org:esmvaltool`` to preview
the documentation.
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
CircleCI_ will build the documentation with the command:

.. code-block:: bash

   python setup.py build_sphinx --warning-is-error

This will catch mistakes that can be detected automatically.

The configuration file for Sphinx_ is
`doc/shinx/source/conf.py <https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/conf.py>`_.

See :ref:`esmvaltool:esmvalcore-documentation-integration` for information on
how the ESMValCore documentation is integrated into the complete ESMValTool
project documentation on readthedocs.

When reviewing a pull request, always check that the documentation checks
shown below the pull request were successful.

.. _tests:

Tests
-----

To check that the code works correctly, there tests available in the
`tests <https://github.com/ESMValGroup/ESMValCore/tree/main/tests>`__ directory.
We use `pytest <https://docs.pytest.org>`_ to write and run our tests.

Contributions to ESMValCore should be
`covered by unit tests <https://the-turing-way.netlify.app/reproducible-research/testing/testing-guidance.html#aim-to-have-a-good-code-coverage>`_.
Have a look at the existing tests in the ``tests`` directory for inspiration on
how to write your own tests.
If you do not know how to start with writing unit tests, ask the
`@ESMValGroup/tech-reviewers`_ for help by commenting on the pull request and
they will try to help you.
It is also recommended that you have a look at the pytest_ documentation at some
point when you start writing your own tests.

Running tests
~~~~~~~~~~~~~

To run the tests on your own computer, go to the directory where the repository
is cloned and run the command

.. code-block:: bash

   pytest

Optionally you can skip tests which require additional dependencies for
supported diagnostic script languages by adding ``-m 'not installation'`` to the
previous command. To only run tests from a single file, run the command

.. code-block:: bash

   pytest tests/unit/test_some_file.py

If you would like to avoid loading the default pytest configuration from
`setup.cfg <https://github.com/ESMValGroup/ESMValCore/blob/main/setup.cfg>`_
because this can be a bit slow for running just a few tests, use

.. code-block:: bash

   pytest -c /dev/null tests/unit/test_some_file.py

Use

.. code-block:: bash

    pytest --help

for more information on the available commands.

Whenever you make a pull request or push new commits to an existing pull
request, the tests in the ``tests`` directory of the branch associated with the
pull request will be run automatically on CircleCI_.
The results appear at the bottom of the pull request.
Click on 'Details' for more information on a specific test job.

When reviewing a pull request, always check that all test jobs on CircleCI_ were
successful.

Test coverage
~~~~~~~~~~~~~

To check which parts of your code are `covered by unit tests`_, open the file
``test-reports/coverage_html/index.html`` (available after running a ``pytest``
command) and browse to the relevant file.

CircleCI will upload the coverage results from running the tests to codecov and
Codacy.
`codecov <https://app.codecov.io/gh/ESMValGroup/ESMValCore/pulls>`_ is a service
that will comment on pull requests with a summary of the test coverage.
If codecov_ reports that the coverage has decreased, check the report and add
additional tests.
Alternatively, it is also possible to view code coverage on Codacy_ (click the
Files tab) and CircleCI_ (open the ``tests`` job and click the ARTIFACTS tab).
To see some of the results on CircleCI, Codacy, or codecov, you may need to log
in; you can do so using your GitHub account.

When reviewing a pull request, always check that new code is covered by unit
tests and codecov_ reports an increased coverage.

.. _sample_data_tests:

Sample data
~~~~~~~~~~~

New or modified preprocessor functions should preferably also be tested using
the sample data.
These tests are located in
`tests/sample_data <https://github.com/ESMValGroup/ESMValCore/tree/main/tests/sample_data>`__.
Please mark new tests that use the sample data with the
`decorator <https://docs.python.org/3/glossary.html#term-decorator>`__
``@pytest.mark.use_sample_data``.

The `ESMValTool_sample_data <https://github.com/ESMValGroup/ESMValTool_sample_data>`_
repository contains samples of CMIP6 data for testing ESMValCore.
The `ESMValTool-sample-data <https://pypi.org/project/ESMValTool-sample-data/>`_
package is installed as part of the developer dependencies.
The size of the package is relatively small (~ 100 MB), so it can be easily
downloaded and distributed.

Preprocessing the sample data can be time-consuming, so some
intermediate results are cached by pytest to make the tests run faster.
If you suspect the tests are failing because the cache is invalid, clear it by
running

.. code-block:: bash

   pytest --cache-clear

To avoid running the time consuming tests that use sample data altogether, run

.. code-block:: bash

   pytest -m "not use_sample_data"


Automated testing
~~~~~~~~~~~~~~~~~

Whenever you make a pull request or push new commits to an existing pull
request, the tests in the ``tests`` of the branch associated with the
pull request will be run automatically on CircleCI_.

Every night, more extensive tests are run to make sure that problems with the
installation of the tool are discovered by the development team before users
encounter them.
These nightly tests have been designed to follow the installation procedures
described in the documentation, e.g. in the :ref:`install` chapter.
The nightly tests are run using both CircleCI and GitHub Actions.
The result of the tests ran by CircleCI can be seen on the
`CircleCI project page <https://app.circleci.com/pipelines/github/ESMValGroup/ESMValCore?branch=main>`__
and the result of the tests ran by GitHub Actions can be viewed on the
`Actions tab <https://github.com/ESMValGroup/ESMValCore/actions>`__
of the repository.

The configuration of the tests run by CircleCI can be found in the directory
`.circleci <https://github.com/ESMValGroup/ESMValCore/blob/main/.circleci>`__,
while the configuration of the tests run by GitHub Actions can be found in the
directory
`.github/workflows <https://github.com/ESMValGroup/ESMValCore/blob/main/.github/workflows>`__.

.. _backward_compatibility:

Backward compatibility
----------------------

The ESMValCore package is used by many people to run their recipes.
Many of these recipes are maintained in the public
`ESMValTool <https://github.com/ESMValGroup/ESMValTool>`_ repository, but
there are also users who choose not to share their work there.
While our commitment is first and foremost to users who do share their recipes
in the ESMValTool repository, we still try to be nice to all of the ESMValCore
users.
When making changes, e.g. to the :ref:`recipe format <recipe_overview>`, the
:ref:`diagnostic script interface <interfaces>`, the public
:ref:`Python API <api>`, or the :ref:`configuration file format <config>`,
keep in mind that this may affect many users.
To keep the tool user friendly, try to avoid making changes that are not
backward compatible, i.e. changes that require users to change their existing
recipes, diagnostics, configuration files, or scripts.

If you really must change the public interfaces of the tool, always discuss this
with the `@ESMValGroup/esmvaltool-coreteam`_.
Try to deprecate the feature first by issuing a :py:class:`DeprecationWarning`
using the :py:mod:`warnings` module and schedule it for removal three
`minor versions <https://semver.org/>`__ from the latest released version.
For example, when you deprecate a feature in a pull request that will be
included in version 2.3, that feature could be removed in version 2.5.
Mention the version in which the feature will be removed in the deprecation
message.
Label the pull request with the
`deprecated feature <https://github.com/ESMValGroup/ESMValCore/labels/deprecated%20feature>`__
label.
When deprecating a feature, please follow up by actually removing the feature
in due course.

If you must make backward incompatible changes, you need to update the available
recipes in ESMValTool and link the ESMValTool pull request(s) in the ESMValCore
pull request description.
You can ask the `@ESMValGroup/esmvaltool-recipe-maintainers`_ for help with
updating existing recipes, but please be considerate of their time.

When reviewing a pull request, always check for backward incompatible changes
and make sure they are needed and have been discussed with the
`@ESMValGroup/esmvaltool-coreteam`_.
Also, make sure the author of the pull request has created the accompanying pull
request(s) to update the ESMValTool, before merging the ESMValCore pull request.

.. _dependencies:

Dependencies
------------

Before considering adding a new dependency, carefully check that the
`license <https://the-turing-way.netlify.app/reproducible-research/licensing/licensing-software.html>`__
of the dependency you want to add and any of its dependencies are
`compatible <https://the-turing-way.netlify.app/reproducible-research/licensing/licensing-compatibility.html>`__
with the
`Apache 2.0 <https://github.com/ESMValGroup/ESMValCore/blob/main/LICENSE/>`_
license that applies to the ESMValCore.
Note that GPL version 2 license is considered incompatible with the Apache 2.0
license, while the compatibility of GPL version 3 license with the Apache 2.0
license is questionable.
See this `statement <https://www.apache.org/licenses/GPL-compatibility.html>`__
by the authors of the Apache 2.0 license for more information.

When adding or removing dependencies, please consider applying the changes in
the following files:

- ``environment.yml``
  contains development dependencies that cannot be installed from
  `PyPI <https://pypi.org/>`_
- ``docs/requirements.txt``
  contains Python dependencies needed to build the documentation that can be
  installed from PyPI
- ``docs/conf.py``
  contains a list of Python dependencies needed to build the documentation that
  cannot be installed from PyPI and need to be mocked when building the
  documentation.
  (We do not use conda to build the documentation because this is too time
  consuming.)
- ``setup.py``
  contains all Python dependencies, regardless of their installation source
- ``package/meta.yaml``
  contains dependencies for the conda package; all Python and compiled
  dependencies that can be installed from conda should be listed here

Note that packages may have a different name on
`conda-forge <https://conda-forge.org/>`__ than on PyPI_.

Several test jobs on CircleCI_ related to the installation of the tool will only
run if you change the dependencies.
These will be skipped for most pull requests.

When reviewing a pull request where dependencies are added or removed, always
check that the changes have been applied in all relevant files.

.. _authors:

List of authors
---------------

If you make a contribution to ESMValCore and you would like to be listed as an
author (e.g. on `Zenodo <https://zenodo.org/record/4525749>`__), please add your
name to the list of authors in ``CITATION.cff`` and generate the entry for the
``.zenodo.json`` file by running the commands

::

   pip install cffconvert
   cffconvert --ignore-suspect-keys --outputformat zenodo --outfile .zenodo.json

Presently, this method unfortunately discards entries `communities`
and `grants` from that file; please restore them manually, or
alternately proceed with the addition manually

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
pull request, please check if there is an open issue that reports the problem.
Create one if there is no issue yet.
You can attract the attention of the `@ESMValGroup/esmvaltool-coreteam`_ by
mentioning them in the issue if it looks like no-one is working on solving the
problem yet.
The issue needs to be fixed in a separate pull request first.
After that has been merged into the ``main`` branch and all checks on this
branch are green again, merge it into your own branch to get the tests to pass.

When reviewing a pull request, always make sure that all checks were successful.
If the Codacy check keeps failing, please run prospector locally.
If necessary, ask the pull request author to do the same and to address the
reported issues.
See the section on code_quality_ for more information.
Never merge a pull request with failing CircleCI or readthedocs checks.


.. _how-to-make-a-release:

Making a release
----------------

The release manager makes the release, assisted by the release manager of the
previous release, or if that person is not available, another previous release
manager. Perform the steps listed below with two persons, to reduce the risk of
error.

To make a new release of the package, follow these steps:

1. Check the tests on GitHub Actions and CircleCI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check the ``nightly``
`build on CircleCI <https://circleci.com/gh/ESMValGroup/ESMValCore/tree/main>`__
and the
`GitHub Actions run <https://github.com/ESMValGroup/ESMValCore/actions>`__.
All tests should pass before making a release (branch).

2. Create a release branch
~~~~~~~~~~~~~~~~~~~~~~~~~~
Create a branch off the ``main`` branch and push it to GitHub.
Ask someone with administrative permissions to set up branch protection rules
for it so only you and the person helping you with the release can push to it.
Announce the name of the branch in an issue and ask the members of the
`ESMValTool development team <https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam>`__
to run their favourite recipe using this branch.

3. Increase the version number
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The version number is stored in ``esmvalcore/_version.py``,
``package/meta.yaml``, ``CITATION.cff``. Make sure to update all files.
Also update the release date in ``CITATION.cff``.
See https://semver.org for more information on choosing a version number.
Make a pull request and get it merged into ``main`` and cherry pick it into
the release branch.

4. Add release notes
~~~~~~~~~~~~~~~~~~~~
Use the script
:ref:`esmvaltool/utils/draft_release_notes.py <esmvaltool:draft_release_notes.py>`
to create create a draft of the release notes.
This script uses the titles and labels of merged pull requests since the
previous release.
Review the results, and if anything needs changing, change it on GitHub and
re-run the script until the changelog looks acceptable.
Copy the result to the file ``doc/changelog.rst``.
Make a pull request and get it merged into ``main`` and cherry pick it into
the release branch..

5. Cherry pick bugfixes into the release branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If a bug is found and fixed (i.e. pull request merged into the
``main`` branch) during the period of testing, use the command
``git cherry-pick`` to include the commit for this bugfix into
the release branch.
When the testing period is over, make a pull request to update
the release notes with the latest changes, get it merged into
``main`` and cherry-pick it into the release branch.

6. Make the release on GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Do a final check that all tests on CircleCI and GitHub Actions completed
successfully.
Then click the
`releases tab <https://github.com/ESMValGroup/ESMValCore/releases>`__
and create the new release from the release branch (i.e. not from ``main``).

7. Create and upload the Conda package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The package is automatically uploaded to the
`ESMValGroup conda channel <https://anaconda.org/esmvalgroup/esmvalcore>`__
by a GitHub action (note that this is an obsolete procedure for the main package upload,
since the main package is now uploaded to
`conda-forge conda channel <https://anaconda.org/conda-forge>`__ via
the upload to PyPi, but we still upload to the esmvalgroup channel as a backup option;
also the upload to esmvalcore gives us a chance to verify it immediately after upload).
If this has failed for some reason, build and upload the package manually by
following the instructions below.

Follow these steps to create a new conda package:

-  Check out the tag corresponding to the release,
   e.g. ``git checkout tags/v2.1.0``
-  Make sure your current working directory is clean by checking the output
   of ``git status`` and by running ``git clean -xdf`` to remove any files
   ignored by git.
-  Edit ``package/meta.yaml`` and uncomment the lines starting with ``git_rev`` and
   ``git_url``, remove the line starting with ``path`` in the ``source``
   section.
-  Activate the base environment ``conda activate base``
-  Install the required packages:
   ``conda install -y conda-build conda-verify ripgrep anaconda-client``
-  Run ``conda build package -c conda-forge`` to build the
   conda package
-  If the build was successful, upload the package to the esmvalgroup
   conda channel, e.g.
   ``anaconda upload --user esmvalgroup /path/to/conda/conda-bld/noarch/esmvalcore-2.3.1-py_0.tar.bz2``.

8. Create and upload the PyPI package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The package is automatically uploaded to the
`PyPI <https://pypi.org/project/ESMValCore/>`__
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
   ``ESMValCore-2.3.1-py3-none-any.whl`` and ``ESMValCore-2.3.1.tar.gz``.
-  Upload the package:
   ``python3 -m twine upload dist/*``
   You will be prompted for an API token if you have not set this up
   before, see
   `here <https://pypi.org/help/#apitoken>`__ for more information.

You can read more about this in
`Packaging Python Projects <https://packaging.python.org/tutorials/packaging-projects/>`__.


.. _`@ESMValGroup/esmvaltool-coreteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-coreteam
.. _`@ESMValGroup/esmvaltool-developmentteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam
.. _`@ESMValGroup/tech-reviewers`: https://github.com/orgs/ESMValGroup/teams/tech-reviewers
.. _`@ESMValGroup/science-reviewers`: https://github.com/orgs/ESMValGroup/teams/science-reviewers
.. _`@ESMValGroup/esmvaltool-recipe-maintainers`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-recipe-maintainers
.. _interfaces:

Diagnostic script interfaces
============================

In order to communicate with diagnostic scripts, ESMValCore uses YAML files.
The YAML files provided by ESMValCore to the diagnostic script tell the diagnostic script the settings that were provided in the recipe and where to find the pre-processed input data.
On the other hand, the YAML file provided by the diagnostic script to ESMValCore tells ESMValCore which pre-processed data was used to create what plots.
The latter is optional, but needed for recording provenance.

Provenance
----------
When ESMValCore (the ``esmvaltool`` command) runs a recipe, it will first find all data and run the default preprocessor steps plus any
additional preprocessing steps defined in the recipe. Next it will run the diagnostic script defined in the recipe
and finally it will store provenance information. Provenance information is stored in the
`W3C PROV XML format <https://www.w3.org/TR/prov-xml/>`_.
To read in and extract information, or to plot these files, the
`prov <https://prov.readthedocs.io>`_ Python package can be used.
In addition to provenance information, a caption is also added to the plots.

.. _interface_esmvalcore_diagnostic:

Information provided by ESMValCore to the diagnostic script
-----------------------------------------------------------
To provide the diagnostic script with the information it needs to run (e.g. location of input data, various settings),
the ESMValCore creates a YAML file called settings.yml and provides the path to this file as the first command line
argument to the diagnostic script.

The most interesting settings provided in this file are

.. code-block:: yaml

  run_dir:  /path/to/recipe_output/run/diagnostic_name/script_name
  work_dir: /path/to/recipe_output/work/diagnostic_name/script_name
  plot_dir: /path/to/recipe_output/plots/diagnostic_name/script_name
  input_files:
    - /path/to/recipe_output/preproc/diagnostic_name/ta/metadata.yml
    - /path/to/recipe_output/preproc/diagnostic_name/pr/metadata.yml

Custom settings in the script section of the recipe will also be made available in this file.

There are three directories defined:

- :code:`run_dir` use this for storing temporary files
- :code:`work_dir` use this for storing NetCDF files containing the data used to make a plot
- :code:`plot_dir` use this for storing plots

Finally :code:`input_files` is a list of YAML files, containing a description of the preprocessed data. Each entry in these
YAML files is a path to a preprocessed file in NetCDF format, with a list of various attributes.
An example preprocessor metadata.yml file could look like this:

.. code-block:: yaml

  ? /path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc
  : alias: GFDL-ESM2G
    cmor_table: CMIP5
    dataset: GFDL-ESM2G
    diagnostic: diagnostic_name
    end_year: 2002
    ensemble: r1i1p1
    exp: historical
    filename: /path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc
    frequency: mon
    institute: [NOAA-GFDL]
    long_name: Precipitation
    mip: Amon
    modeling_realm: [atmos]
    preprocessor: preprocessor_name
    project: CMIP5
    recipe_dataset_index: 1
    reference_dataset: MPI-ESM-LR
    short_name: pr
    standard_name: precipitation_flux
    start_year: 2000
    units: kg m-2 s-1
    variable_group: pr
  ? /path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc
  : alias: MPI-ESM-LR
    cmor_table: CMIP5
    dataset: MPI-ESM-LR
    diagnostic: diagnostic_name
    end_year: 2002
    ensemble: r1i1p1
    exp: historical
    filename: /path/to/recipe_output/preproc/diagnostic1/pr/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc
    frequency: mon
    institute: [MPI-M]
    long_name: Precipitation
    mip: Amon
    modeling_realm: [atmos]
    preprocessor: preprocessor_name
    project: CMIP5
    recipe_dataset_index: 2
    reference_dataset: MPI-ESM-LR
    short_name: pr
    standard_name: precipitation_flux
    start_year: 2000
    units: kg m-2 s-1
    variable_group: pr


.. _interface_diagnostic_esmvalcore:

Information provided by the diagnostic script to ESMValCore
-----------------------------------------------------------

After the diagnostic script has finished running, ESMValCore will try to store provenance information. In order to
link the produced files to input data, the diagnostic script needs to store a YAML file called :code:`diagnostic_provenance.yml`
in its :code:`run_dir`.

For every output file (netCDF files, plot files, etc.) produced by the diagnostic script, there should be an entry in the :code:`diagnostic_provenance.yml` file.
The name of each entry should be the path to the file.
Each output file entry should at least contain the following items:

- :code:`ancestors` a list of input files used to create the plot.
- :code:`caption` a caption text for the plot.

Each file entry can also contain items from the categories defined in the file :code:`esmvaltool/config_references.yml`.
The short entries will automatically be replaced by their longer equivalent in the final provenance records.
It is possible to add custom provenance information by adding custom items to entries.

An example :code:`diagnostic_provenance.yml` file could look like this

.. code-block:: yaml

  ? /path/to/recipe_output/work/diagnostic_name/script_name/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_pr_2000-2002_mean.nc
  : ancestors:[/path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_pr_2000-2002.nc]
    authors: [andela_bouwe, righi_mattia]
    caption: Average Precipitation between 2000 and 2002 according to GFDL-ESM2G.
    domains: [global]
    plot_types: [zonal]
    references: [acknow_project]
    statistics: [mean]

  ? /path/to/recipe_output/plots/diagnostic_name/script_name/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_pr_2000-2002_mean.png
  : ancestors:[/path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_pr_2000-2002.nc]
    authors: [andela_bouwe, righi_mattia]
    caption: Average Precipitation between 2000 and 2002 according to GFDL-ESM2G.
    domains: [global]
    plot_types: ['zonal']
    references: [acknow_project]
    statistics: [mean]

You can check whether your diagnostic script successfully provided the provenance information to the ESMValCore by
checking the following points:

  - for each output file in the ``work_dir`` and ``plot_dir``, a file with the same
    name, but ending with ``_provenance.xml`` is created
  - the output file is shown on the ``index.html`` page
  - there were no warning messages in the log related to provenance

See :ref:`esmvaltool:recording-provenance` for more extensive usage notes.
Welcome to ESMValTool's documentation!
======================================

.. include:: _sidebar.rst.inc

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

Changelog
=========

.. _changelog-v2-4-0:

v2.4.0
------

Highlights
~~~~~~~~~~

- ESMValCore now has the ability to automatically download missing data from ESGF. For details, see :ref:`Data Retrieval<data-retrieval>`.
- ESMValCore now also can resume an earlier run. This is useful to re-use expensive preprocessor results. For details, see :ref:`Running<running>`.

This release includes

Bug fixes
~~~~~~~~~

-  Crop on the ID-selected region(s) and not on the whole shapefile (`#1151 <https://github.com/ESMValGroup/ESMValCore/pull/1151>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add 'comment' to list of removed attributes (`#1244 <https://github.com/ESMValGroup/ESMValCore/pull/1244>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Speed up multimodel statistics and fix bug in peak computation (`#1301 <https://github.com/ESMValGroup/ESMValCore/pull/1301>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  No longer make plots of provenance (`#1307 <https://github.com/ESMValGroup/ESMValCore/pull/1307>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  No longer embed provenance in output files (`#1306 <https://github.com/ESMValGroup/ESMValCore/pull/1306>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Removed automatic addition of areacello to obs4mips datasets (`#1316 <https://github.com/ESMValGroup/ESMValCore/pull/1316>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Pin docutils <0.17 to fix bullet lists on readthedocs (`#1320 <https://github.com/ESMValGroup/ESMValCore/pull/1320>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix obs4MIPs capitalization (`#1328 <https://github.com/ESMValGroup/ESMValCore/pull/1328>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix Python 3.7 tests (`#1330 <https://github.com/ESMValGroup/ESMValCore/pull/1330>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Handle fx variables in `extract_levels` and some time operations (`#1269 <https://github.com/ESMValGroup/ESMValCore/pull/1269>`__) `sloosvel <https://github.com/sloosvel>`__
-  Refactored mask regridding for irregular grids (fixes #772) (`#865 <https://github.com/ESMValGroup/ESMValCore/pull/865>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix `da.broadcast_to` call when the fx cube has different shape than target data cube (`#1350 <https://github.com/ESMValGroup/ESMValCore/pull/1350>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add tests for _aggregate_time_fx (`#1354 <https://github.com/ESMValGroup/ESMValCore/pull/1354>`__) `sloosvel <https://github.com/sloosvel>`__
-  Fix extra facets (`#1360 <https://github.com/ESMValGroup/ESMValCore/pull/1360>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Pin pip!=21.3 to avoid pypa/pip#10573 with editable installs (`#1359 <https://github.com/ESMValGroup/ESMValCore/pull/1359>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add a custom `date2num` function to deal with changes in cftime (`#1373 <https://github.com/ESMValGroup/ESMValCore/pull/1373>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Removed custom version of `AtmosphereSigmaFactory` (`#1382 <https://github.com/ESMValGroup/ESMValCore/pull/1382>`__) `Manuel Schlund <https://github.com/schlunma>`__

Deprecations
~~~~~~~~~~~~

-  Remove write_netcdf and write_plots from config-user.yml (`#1300 <https://github.com/ESMValGroup/ESMValCore/pull/1300>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

-  Add link to plot directory in index.html (`#1256 <https://github.com/ESMValGroup/ESMValCore/pull/1256>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Work around issue with yapf not following PEP8 (`#1277 <https://github.com/ESMValGroup/ESMValCore/pull/1277>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update the core development team (`#1278 <https://github.com/ESMValGroup/ESMValCore/pull/1278>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update the documentation of the provenance interface (`#1305 <https://github.com/ESMValGroup/ESMValCore/pull/1305>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update version number to first release candidate 2.4.0rc1 (`#1363 <https://github.com/ESMValGroup/ESMValCore/pull/1363>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update to new ESMValTool logo (`#1374 <https://github.com/ESMValGroup/ESMValCore/pull/1374>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update version number for third release candidate 2.4.0rc3 (`#1384 <https://github.com/ESMValGroup/ESMValCore/pull/1384>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update changelog for 2.4.0rc3 (`#1385 <https://github.com/ESMValGroup/ESMValCore/pull/1385>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update version number to final 2.4.0 release (`#1389 <https://github.com/ESMValGroup/ESMValCore/pull/1389>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update changelog for 2.4.0 (`#1366 <https://github.com/ESMValGroup/ESMValCore/pull/1366>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Fixes for datasets
~~~~~~~~~~~~~~~~~~

-  Add fix for differing latitude coordinate between historical and ssp585 in MPI-ESM1-2-HR r2i1p1f1 (`#1292 <https://github.com/ESMValGroup/ESMValCore/pull/1292>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add fixes for time and latitude coordinate of EC-Earth3 r3i1p1f1 (`#1290 <https://github.com/ESMValGroup/ESMValCore/pull/1290>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Apply latitude fix to all CCSM4 variables (`#1295 <https://github.com/ESMValGroup/ESMValCore/pull/1295>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix lat and lon bounds for FGOALS-g3 mrsos (`#1289 <https://github.com/ESMValGroup/ESMValCore/pull/1289>`__) `Thomas Crocker <https://github.com/thomascrocker>`__
-  Add grid fix for tos in fgoals-f3-l (`#1326 <https://github.com/ESMValGroup/ESMValCore/pull/1326>`__) `sloosvel <https://github.com/sloosvel>`__
-  Add fix for CIESM pr (`#1344 <https://github.com/ESMValGroup/ESMValCore/pull/1344>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix DRS for IPSLCM : split attribute 'freq' into : 'out' and 'freq' (`#1304 <https://github.com/ESMValGroup/ESMValCore/pull/1304>`__) `Stéphane Sénési - work <https://github.com/senesis>`__

CMOR standard
~~~~~~~~~~~~~

-  Remove history attribute from coords (`#1276 <https://github.com/ESMValGroup/ESMValCore/pull/1276>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Increased flexibility of CMOR checks for datasets with generic alevel coordinates (`#1032 <https://github.com/ESMValGroup/ESMValCore/pull/1032>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Automatically fix small deviations in vertical levels (`#1177 <https://github.com/ESMValGroup/ESMValCore/pull/1177>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Adding standard names to the custom tables of the `rlns` and `rsns` variables (`#1386 <https://github.com/ESMValGroup/ESMValCore/pull/1386>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__

Preprocessor
~~~~~~~~~~~~

-  Implemented fully lazy climate_statistics (`#1194 <https://github.com/ESMValGroup/ESMValCore/pull/1194>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Run the multimodel statistics preprocessor last (`#1299 <https://github.com/ESMValGroup/ESMValCore/pull/1299>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Pin Python to 3.9 in environment.yml on CircleCI and skip mypy tests in conda build (`#1176 <https://github.com/ESMValGroup/ESMValCore/pull/1176>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improving test coverage for _task.py (`#514 <https://github.com/ESMValGroup/ESMValCore/pull/514>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Upload coverage to codecov (`#1190 <https://github.com/ESMValGroup/ESMValCore/pull/1190>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve codecov status checks (`#1195 <https://github.com/ESMValGroup/ESMValCore/pull/1195>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix curl install in CircleCI (`#1228 <https://github.com/ESMValGroup/ESMValCore/pull/1228>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Drop support for Python 3.6 (`#1200 <https://github.com/ESMValGroup/ESMValCore/pull/1200>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Allow more recent version of `scipy` (`#1182 <https://github.com/ESMValGroup/ESMValCore/pull/1182>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Speed up conda build `conda_build` Circle test by using `mamba` solver via `boa` (and use it for Github Actions test too) (`#1243 <https://github.com/ESMValGroup/ESMValCore/pull/1243>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix numpy deprecation warnings (`#1274 <https://github.com/ESMValGroup/ESMValCore/pull/1274>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Unpin upper bound for iris (previously was at <3.0.4)  (`#1275 <https://github.com/ESMValGroup/ESMValCore/pull/1275>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Modernize `conda_install` test on Circle CI by installing from conda-forge with Python 3.9 and change install instructions in documentation (`#1280 <https://github.com/ESMValGroup/ESMValCore/pull/1280>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Run a nightly Github Actions workflow to monitor tests memory per test (configurable for other metrics too) (`#1284 <https://github.com/ESMValGroup/ESMValCore/pull/1284>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Speed up tests of tasks (`#1302 <https://github.com/ESMValGroup/ESMValCore/pull/1302>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix upper case to lower case variables and functions for flake compliance in `tests/unit/preprocessor/_regrid/test_extract_levels.py` (`#1347 <https://github.com/ESMValGroup/ESMValCore/pull/1347>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Cleaned up a bit Github Actions workflows (`#1345 <https://github.com/ESMValGroup/ESMValCore/pull/1345>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update circleci jobs: renaming tests to more descriptive names and removing conda build test (`#1351 <https://github.com/ESMValGroup/ESMValCore/pull/1351>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Pin iris to latest `>=3.1.0` (`#1341 <https://github.com/ESMValGroup/ESMValCore/pull/1341>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Installation
~~~~~~~~~~~~

-  Pin esmpy to anything but 8.1.0 since that particular one changes the CPU affinity (`#1310 <https://github.com/ESMValGroup/ESMValCore/pull/1310>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Add a more friendly and useful message when using default config file (`#1233 <https://github.com/ESMValGroup/ESMValCore/pull/1233>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Replace os.walk by glob.glob in data finder (only look for data in the specified locations) (`#1261 <https://github.com/ESMValGroup/ESMValCore/pull/1261>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Machine-specific directories for auxiliary data in the `config-user.yml` file (`#1268 <https://github.com/ESMValGroup/ESMValCore/pull/1268>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Add an option to download missing data from ESGF (`#1217 <https://github.com/ESMValGroup/ESMValCore/pull/1217>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Speed up provenance recording (`#1327 <https://github.com/ESMValGroup/ESMValCore/pull/1327>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve results web page (`#1332 <https://github.com/ESMValGroup/ESMValCore/pull/1332>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Move institutes from config-developer.yml to default extra facets config and add wildcard support for extra facets (`#1259 <https://github.com/ESMValGroup/ESMValCore/pull/1259>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add support for re-using preprocessor output from previous runs (`#1321 <https://github.com/ESMValGroup/ESMValCore/pull/1321>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Log fewer messages to screen and hide stack trace for known recipe errors (`#1296 <https://github.com/ESMValGroup/ESMValCore/pull/1296>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Log ESMValCore and ESMValTool versions when running (`#1263 <https://github.com/ESMValGroup/ESMValCore/pull/1263>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add "grid" as a tag to the output file template for CMIP6 (`#1356 <https://github.com/ESMValGroup/ESMValCore/pull/1356>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Implemented ICON project to read native ICON model output (`#1079 <https://github.com/ESMValGroup/ESMValCore/pull/1079>`__) `Brei Soliño <https://github.com/bsolino>`__


.. _changelog-v2-3-1:

v2.3.1
------

This release includes

Bug fixes
~~~~~~~~~

-  Update config-user.yml template with correct drs entries for CEDA-JASMIN (`#1184 <https://github.com/ESMValGroup/ESMValCore/pull/1184>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Enhancing MIROC5 fix for hfls and evspsbl (`#1192 <https://github.com/ESMValGroup/ESMValCore/pull/1192>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fix alignment of daily data with inconsistent calendars in multimodel statistics (`#1212 <https://github.com/ESMValGroup/ESMValCore/pull/1212>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Pin cf-units, remove github actions test for Python 3.6 and fix test_access1_0 and test_access1_3 to use cf-units for comparisons (`#1197 <https://github.com/ESMValGroup/ESMValCore/pull/1197>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fixed search for fx files when no ``mip`` is given (`#1216 <https://github.com/ESMValGroup/ESMValCore/pull/1216>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Make sure climate statistics always returns original dtype (`#1237 <https://github.com/ESMValGroup/ESMValCore/pull/1237>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Bugfix for regional regridding when non-integer range is passed (`#1231 <https://github.com/ESMValGroup/ESMValCore/pull/1231>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Make sure area_statistics preprocessor always returns original dtype (`#1239 <https://github.com/ESMValGroup/ESMValCore/pull/1239>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add "." (dot) as allowed separation character for the time range group (`#1248 <https://github.com/ESMValGroup/ESMValCore/pull/1248>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Documentation
~~~~~~~~~~~~~

-  Add a link to the instructions to use pre-installed versions on HPC clusters (`#1186 <https://github.com/ESMValGroup/ESMValCore/pull/1186>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Bugfix release: set version to 2.3.1 (`#1253 <https://github.com/ESMValGroup/ESMValCore/pull/1253>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Fixes for datasets
~~~~~~~~~~~~~~~~~~

-  Set circular attribute in MCM-UA-1-0 fix (`#1178 <https://github.com/ESMValGroup/ESMValCore/pull/1178>`__) `sloosvel <https://github.com/sloosvel>`__
-  Fixed time coordinate of MIROC-ESM (`#1188 <https://github.com/ESMValGroup/ESMValCore/pull/1188>`__) `Manuel Schlund <https://github.com/schlunma>`__

Preprocessor
~~~~~~~~~~~~

-  Filter warnings about collapsing multi-model dimension in multimodel statistics preprocessor function (`#1215 <https://github.com/ESMValGroup/ESMValCore/pull/1215>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Remove fx variables before computing multimodel statistics (`#1220 <https://github.com/ESMValGroup/ESMValCore/pull/1220>`__) `sloosvel <https://github.com/sloosvel>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Pin Python to 3.9 in environment.yml on CircleCI and skip mypy tests in conda build (`#1176 <https://github.com/ESMValGroup/ESMValCore/pull/1176>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Installation
~~~~~~~~~~~~

-  Pin lower bound for iris to 3.0.2 (`#1206 <https://github.com/ESMValGroup/ESMValCore/pull/1206>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin `iris<3.0.4` to ensure we still (sort of) support Python 3.6 (`#1252 <https://github.com/ESMValGroup/ESMValCore/pull/1252>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Add test to verify behaviour for scalar height coord for tas in multi-model (`#1209 <https://github.com/ESMValGroup/ESMValCore/pull/1209>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Sort missing years in "No input data available for years" message (`#1225 <https://github.com/ESMValGroup/ESMValCore/pull/1225>`__) `Lee de Mora <https://github.com/ledm>`__


.. _changelog-v2-3-0:

v2.3.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Extend preprocessor multi_model_statistics to handle data with "altitude" coordinate (`#1010 <https://github.com/ESMValGroup/ESMValCore/pull/1010>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Remove scripts included with CMOR tables (`#1011 <https://github.com/ESMValGroup/ESMValCore/pull/1011>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Avoid side effects in extract_season (`#1019 <https://github.com/ESMValGroup/ESMValCore/pull/1019>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Use nearest scheme to avoid interpolation errors with masked data in regression test (`#1021 <https://github.com/ESMValGroup/ESMValCore/pull/1021>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Move _get_time_bounds from preprocessor._time to cmor.check to avoid circular import with cmor module (`#1037 <https://github.com/ESMValGroup/ESMValCore/pull/1037>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix test that makes conda build fail (`#1046 <https://github.com/ESMValGroup/ESMValCore/pull/1046>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix 'positive' attribute for rsns/rlns variables (`#1051 <https://github.com/ESMValGroup/ESMValCore/pull/1051>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Added preprocessor mask_multimodel (`#767 <https://github.com/ESMValGroup/ESMValCore/pull/767>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix bug when fixing bounds after fixing longitude values (`#1057 <https://github.com/ESMValGroup/ESMValCore/pull/1057>`__) `sloosvel <https://github.com/sloosvel>`__
-  Run conda build parallel AND sequential tests (`#1065 <https://github.com/ESMValGroup/ESMValCore/pull/1065>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add key to id_prop (`#1071 <https://github.com/ESMValGroup/ESMValCore/pull/1071>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Fix bounds after reversing coordinate values (`#1061 <https://github.com/ESMValGroup/ESMValCore/pull/1061>`__) `sloosvel <https://github.com/sloosvel>`__
-  Fixed --skip-nonexistent option (`#1093 <https://github.com/ESMValGroup/ESMValCore/pull/1093>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Do not consider CMIP5 variable sit to be the same as sithick from CMIP6 (`#1033 <https://github.com/ESMValGroup/ESMValCore/pull/1033>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve finding date range in filenames (enforces separators) (`#1145 <https://github.com/ESMValGroup/ESMValCore/pull/1145>`__) `Stéphane Sénési - work <https://github.com/senesis>`__
-  Review fx handling (`#1147 <https://github.com/ESMValGroup/ESMValCore/pull/1147>`__) `sloosvel <https://github.com/sloosvel>`__
-  Fix lru cache decorator with explicit call to method (`#1172 <https://github.com/ESMValGroup/ESMValCore/pull/1172>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update _volume.py (`#1174 <https://github.com/ESMValGroup/ESMValCore/pull/1174>`__) `Lee de Mora <https://github.com/ledm>`__

Deprecations
~~~~~~~~~~~~



Documentation
~~~~~~~~~~~~~

-  Final changelog for 2.3.0 (`#1163 <https://github.com/ESMValGroup/ESMValCore/pull/1163>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Set version to 2.3.0 (`#1162 <https://github.com/ESMValGroup/ESMValCore/pull/1162>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix documentation build (`#1006 <https://github.com/ESMValGroup/ESMValCore/pull/1006>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add labels required for linking from ESMValTool docs (`#1038 <https://github.com/ESMValGroup/ESMValCore/pull/1038>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update contribution guidelines (`#1047 <https://github.com/ESMValGroup/ESMValCore/pull/1047>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix basestring references in documentation (`#1106 <https://github.com/ESMValGroup/ESMValCore/pull/1106>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Updated references master to main (`#1132 <https://github.com/ESMValGroup/ESMValCore/pull/1132>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add instructions how to use the central installation at DKRZ-Mistral (`#1155 <https://github.com/ESMValGroup/ESMValCore/pull/1155>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__

Fixes for datasets
~~~~~~~~~~~~~~~~~~

-  Added fixes for various CMIP5 datasets, variable cl (3-dim cloud fraction) (`#1017 <https://github.com/ESMValGroup/ESMValCore/pull/1017>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Added fixes for hybrid level coordinates of CESM2 models (`#882 <https://github.com/ESMValGroup/ESMValCore/pull/882>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Extending LWP fix for CMIP6 models (`#1049 <https://github.com/ESMValGroup/ESMValCore/pull/1049>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add fixes for the net & up radiation variables from ERA5 (`#1052 <https://github.com/ESMValGroup/ESMValCore/pull/1052>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Add derived variable rsus (`#1053 <https://github.com/ESMValGroup/ESMValCore/pull/1053>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Supported `mip`-level fixes (`#1095 <https://github.com/ESMValGroup/ESMValCore/pull/1095>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix erroneous use of `grid_latitude` and `grid_longitude` and cleaned ocean grid fixes (`#1092 <https://github.com/ESMValGroup/ESMValCore/pull/1092>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix for pr of miroc5 (`#1110 <https://github.com/ESMValGroup/ESMValCore/pull/1110>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Ocean depth fix for cnrm_esm2_1, gfdl_esm4, ipsl_cm6a_lr datasets +  mcm_ua_1_0 (`#1098 <https://github.com/ESMValGroup/ESMValCore/pull/1098>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Fix for uas variable of the MCM_UA_1_0 dataset (`#1102 <https://github.com/ESMValGroup/ESMValCore/pull/1102>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fixes for sos and siconc of BCC models (`#1090 <https://github.com/ESMValGroup/ESMValCore/pull/1090>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Run fgco2 fix for all CESM2 models (`#1108 <https://github.com/ESMValGroup/ESMValCore/pull/1108>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Fixes for the siconc variable of CMIP6 models (`#1105 <https://github.com/ESMValGroup/ESMValCore/pull/1105>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fix wrong sign for land surface flux (`#1113 <https://github.com/ESMValGroup/ESMValCore/pull/1113>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Fix for pr of EC_EARTH (`#1116 <https://github.com/ESMValGroup/ESMValCore/pull/1116>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__

CMOR standard
~~~~~~~~~~~~~

-  Format cmor related files (`#976 <https://github.com/ESMValGroup/ESMValCore/pull/976>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Check presence of time bounds and guess them if needed (`#849 <https://github.com/ESMValGroup/ESMValCore/pull/849>`__) `sloosvel <https://github.com/sloosvel>`__
-  Add custom variable "tasaga" (`#1118 <https://github.com/ESMValGroup/ESMValCore/pull/1118>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Find files for CMIP6 DCPP startdates (`#771 <https://github.com/ESMValGroup/ESMValCore/pull/771>`__) `sloosvel <https://github.com/sloosvel>`__

Preprocessor
~~~~~~~~~~~~

-  Update tests for multimodel statistics preprocessor (`#1023 <https://github.com/ESMValGroup/ESMValCore/pull/1023>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Raise in extract_season and extract_month if result is None (`#1041 <https://github.com/ESMValGroup/ESMValCore/pull/1041>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Allow selection of shapes in extract_shape (`#764 <https://github.com/ESMValGroup/ESMValCore/pull/764>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add option for regional regridding to regrid preprocessor (`#1034 <https://github.com/ESMValGroup/ESMValCore/pull/1034>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Load fx variables as cube cell measures / ancillary variables (`#999 <https://github.com/ESMValGroup/ESMValCore/pull/999>`__) `sloosvel <https://github.com/sloosvel>`__
-  Check horizontal grid before regridding (`#507 <https://github.com/ESMValGroup/ESMValCore/pull/507>`__) `Benjamin Müller <https://github.com/BenMGeo>`__
-  Clip irregular grids (`#245 <https://github.com/ESMValGroup/ESMValCore/pull/245>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Use native iris functions in multi-model statistics (`#1150 <https://github.com/ESMValGroup/ESMValCore/pull/1150>`__) `Peter Kalverla <https://github.com/Peter9192>`__

Notebook API (experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~



Automatic testing
~~~~~~~~~~~~~~~~~

-  Report coverage for tests that run on any pull request (`#994 <https://github.com/ESMValGroup/ESMValCore/pull/994>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Install ESMValTool sample data from PyPI (`#998 <https://github.com/ESMValGroup/ESMValCore/pull/998>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix tests for multi-processing with spawn method (i.e. macOSX with Python>3.8) (`#1003 <https://github.com/ESMValGroup/ESMValCore/pull/1003>`__) `Barbara Vreede <https://github.com/bvreede>`__
-  Switch to running the Github Action test workflow every 3 hours in single thread mode to observe if Sementation Faults occur (`#1022 <https://github.com/ESMValGroup/ESMValCore/pull/1022>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Revert to original Github Actions test workflow removing the 3-hourly test run with -n 1 (`#1025 <https://github.com/ESMValGroup/ESMValCore/pull/1025>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Avoid stale cache for multimodel statistics regression tests (`#1030 <https://github.com/ESMValGroup/ESMValCore/pull/1030>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add newer Python versions in OSX to Github Actions (`#1035 <https://github.com/ESMValGroup/ESMValCore/pull/1035>`__) `Barbara Vreede <https://github.com/bvreede>`__
-  Add tests for type annotations with mypy (`#1042 <https://github.com/ESMValGroup/ESMValCore/pull/1042>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Run problematic cmor tests sequentially to avoid segmentation faults on CircleCI (`#1064 <https://github.com/ESMValGroup/ESMValCore/pull/1064>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Test installation of esmvalcore from conda-forge (`#1075 <https://github.com/ESMValGroup/ESMValCore/pull/1075>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Added additional test cases for integration tests of data_finder.py (`#1087 <https://github.com/ESMValGroup/ESMValCore/pull/1087>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Pin cf-units and fix tests (cf-units>=2.1.5) (`#1140 <https://github.com/ESMValGroup/ESMValCore/pull/1140>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix failing CircleCI tests (`#1167 <https://github.com/ESMValGroup/ESMValCore/pull/1167>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix test failing due to fx files chosen differently on different OS's (`#1169 <https://github.com/ESMValGroup/ESMValCore/pull/1169>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Compare datetimes instead of strings in _fixes/cmip5/test_access1_X.py (`#1173 <https://github.com/ESMValGroup/ESMValCore/pull/1173>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin Python to 3.9 in environment.yml on CircleCI and skip mypy tests in conda build (`#1176 <https://github.com/ESMValGroup/ESMValCore/pull/1176>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Installation
~~~~~~~~~~~~

-  Update yamale to version 3 (`#1059 <https://github.com/ESMValGroup/ESMValCore/pull/1059>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Improvements
~~~~~~~~~~~~

-  Refactor diagnostics / tags management (`#939 <https://github.com/ESMValGroup/ESMValCore/pull/939>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Support multiple paths in input_dir (`#1000 <https://github.com/ESMValGroup/ESMValCore/pull/1000>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Generate HTML report with recipe output (`#991 <https://github.com/ESMValGroup/ESMValCore/pull/991>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add timeout to requests.get in _citation.py (`#1091 <https://github.com/ESMValGroup/ESMValCore/pull/1091>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Add SYNDA drs for CMIP5 and CMIP6 (closes #582) (`#583 <https://github.com/ESMValGroup/ESMValCore/pull/583>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add basic support for variable mappings (`#1124 <https://github.com/ESMValGroup/ESMValCore/pull/1124>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Handle IPSL-CM6  (`#1153 <https://github.com/ESMValGroup/ESMValCore/pull/1153>`__) `Stéphane Sénési - work <https://github.com/senesis>`__


.. _changelog-v2-2-0:

v2.2.0
------

Highlights
~~~~~~~~~~

ESMValCore is now using the recently released `Iris 3 <https://scitools-iris.readthedocs.io/en/latest/whatsnew/3.0.html>`__.
We acknowledge that this change may impact your work, as Iris 3 introduces
several changes that are not backward-compatible, but we think that moving forward is the best
decision for the tool in the long term.

This release is also the first one including support for downloading CMIP6 data
using Synda and we have also started supporting Python 3.9. Give it a try!


This release includes

Bug fixes
~~~~~~~~~

-  Fix path settings for DKRZ/Mistral (`#852 <https://github.com/ESMValGroup/ESMValCore/pull/852>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Change logic for calling the diagnostic script to avoid problems with scripts where the executable bit is accidentally set (`#877 <https://github.com/ESMValGroup/ESMValCore/pull/877>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix overwriting in generic level check (`#886 <https://github.com/ESMValGroup/ESMValCore/pull/886>`__) `sloosvel <https://github.com/sloosvel>`__
-  Add double quotes to script args in rerun screen message when using vprof profiling (`#897 <https://github.com/ESMValGroup/ESMValCore/pull/897>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Simplify time handling in multi-model statistics preprocessor (`#685 <https://github.com/ESMValGroup/ESMValCore/pull/685>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Fix links to Iris documentation (`#966 <https://github.com/ESMValGroup/ESMValCore/pull/966>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Bugfix: Fix units for MSWEP data (`#986 <https://github.com/ESMValGroup/ESMValCore/pull/986>`__) `Stef Smeets <https://github.com/stefsmeets>`__

Deprecations
~~~~~~~~~~~~

-  Deprecate defining write_plots and write_netcdf in config-user file (`#808 <https://github.com/ESMValGroup/ESMValCore/pull/808>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

-  Fix numbering of steps in release instructions (`#838 <https://github.com/ESMValGroup/ESMValCore/pull/838>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add labels to changelogs of individual versions for easy reference (`#899 <https://github.com/ESMValGroup/ESMValCore/pull/899>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Make CircleCI badge specific to main branch (`#902 <https://github.com/ESMValGroup/ESMValCore/pull/902>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix docker build badge url (`#906 <https://github.com/ESMValGroup/ESMValCore/pull/906>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Update github PR template (`#909 <https://github.com/ESMValGroup/ESMValCore/pull/909>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Refer to ESMValTool GitHub discussions page in the error message (`#900 <https://github.com/ESMValGroup/ESMValCore/pull/900>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Support automatically closing issues (`#922 <https://github.com/ESMValGroup/ESMValCore/pull/922>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix checkboxes in PR template (`#931 <https://github.com/ESMValGroup/ESMValCore/pull/931>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Change in config-user defaults and documentation with new location for esmeval OBS data on JASMIN (`#958 <https://github.com/ESMValGroup/ESMValCore/pull/958>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update Core Team info (`#942 <https://github.com/ESMValGroup/ESMValCore/pull/942>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Update iris documentation URL for sphinx (`#964 <https://github.com/ESMValGroup/ESMValCore/pull/964>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Set version to 2.2.0 (`#977 <https://github.com/ESMValGroup/ESMValCore/pull/977>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add first draft of v2.2.0 changelog (`#983 <https://github.com/ESMValGroup/ESMValCore/pull/983>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add checkbox in PR template to assign labels (`#985 <https://github.com/ESMValGroup/ESMValCore/pull/985>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Update install.rst (`#848 <https://github.com/ESMValGroup/ESMValCore/pull/848>`__) `bascrezee <https://github.com/bascrezee>`__
-  Change the order of the publication steps (`#984 <https://github.com/ESMValGroup/ESMValCore/pull/984>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add instructions how to use esmvaltool from HPC central installations (`#841 <https://github.com/ESMValGroup/ESMValCore/pull/841>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Fixes for datasets
~~~~~~~~~~~~~~~~~~

-  Fixing unit for derived variable rsnstcsnorm to prevent overcorrection2 (`#846 <https://github.com/ESMValGroup/ESMValCore/pull/846>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Cmip6 fix awi cm 1 1 mr (`#822 <https://github.com/ESMValGroup/ESMValCore/pull/822>`__) `mwjury <https://github.com/mwjury>`__
-  Cmip6 fix ec earth3 veg (`#836 <https://github.com/ESMValGroup/ESMValCore/pull/836>`__) `mwjury <https://github.com/mwjury>`__
-  Changed latitude longitude fix from Tas to AllVars. (`#916 <https://github.com/ESMValGroup/ESMValCore/pull/916>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fix for precipitation (pr) to use ERA5-Land cmorizer (`#879 <https://github.com/ESMValGroup/ESMValCore/pull/879>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Cmip6 fix ec earth3 (`#837 <https://github.com/ESMValGroup/ESMValCore/pull/837>`__) `mwjury <https://github.com/mwjury>`__
-  Cmip6_fix_fgoals_f3_l_Amon_time_bnds (`#831 <https://github.com/ESMValGroup/ESMValCore/pull/831>`__) `mwjury <https://github.com/mwjury>`__
-  Fix for FGOALS-f3-L sftlf (`#667 <https://github.com/ESMValGroup/ESMValCore/pull/667>`__) `mwjury <https://github.com/mwjury>`__
-  Improve ACCESS-CM2 and ACCESS-ESM1-5 fixes and add CIESM and CESM2-WACCM-FV2 fixes for cl, clw and cli (`#635 <https://github.com/ESMValGroup/ESMValCore/pull/635>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add  fixes for cl, cli, clw and tas for several CMIP6 models (`#955 <https://github.com/ESMValGroup/ESMValCore/pull/955>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Dataset fixes for MSWEP (`#969 <https://github.com/ESMValGroup/ESMValCore/pull/969>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Dataset fixes for: ACCESS-ESM1-5, CanESM5, CanESM5 for carbon cycle (`#947 <https://github.com/ESMValGroup/ESMValCore/pull/947>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Fixes for KIOST-ESM (CMIP6) (`#904 <https://github.com/ESMValGroup/ESMValCore/pull/904>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fixes for AWI-ESM-1-1-LR (CMIP6, piControl) (`#911 <https://github.com/ESMValGroup/ESMValCore/pull/911>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__

CMOR standard
~~~~~~~~~~~~~

-  CMOR check generic level coordinates in CMIP6 (`#598 <https://github.com/ESMValGroup/ESMValCore/pull/598>`__) `sloosvel <https://github.com/sloosvel>`__
-  Update CMIP6 tables to 6.9.33 (`#919 <https://github.com/ESMValGroup/ESMValCore/pull/919>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Adding custom variables for tas uncertainty (`#924 <https://github.com/ESMValGroup/ESMValCore/pull/924>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Remove monotonicity coordinate check for unstructured grids (`#965 <https://github.com/ESMValGroup/ESMValCore/pull/965>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Preprocessor
~~~~~~~~~~~~

-  Added clip_start_end_year preprocessor (`#796 <https://github.com/ESMValGroup/ESMValCore/pull/796>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add support for downloading CMIP6 data with Synda (`#699 <https://github.com/ESMValGroup/ESMValCore/pull/699>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add multimodel tests using real data (`#856 <https://github.com/ESMValGroup/ESMValCore/pull/856>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add plev/altitude conversion to extract_levels (`#892 <https://github.com/ESMValGroup/ESMValCore/pull/892>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add possibility of custom season extraction. (`#247 <https://github.com/ESMValGroup/ESMValCore/pull/247>`__) `mwjury <https://github.com/mwjury>`__
-  Adding the ability to derive xch4  (`#783 <https://github.com/ESMValGroup/ESMValCore/pull/783>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Add preprocessor function to resample time and compute x-hourly statistics (`#696 <https://github.com/ESMValGroup/ESMValCore/pull/696>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix duplication in preprocessors DEFAULT_ORDER introduced in `#696 <https://github.com/ESMValGroup/ESMValCore/pull/696>`__  (`#973 <https://github.com/ESMValGroup/ESMValCore/pull/973>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Use consistent precision in multi-model statistics calculation and update reference data for tests (`#941 <https://github.com/ESMValGroup/ESMValCore/pull/941>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Refactor multi-model statistics code to facilitate ensemble stats and lazy evaluation (`#949 <https://github.com/ESMValGroup/ESMValCore/pull/949>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Add option to exclude input cubes in output of multimodel statistics to solve an issue introduced by `#949 <https://github.com/ESMValGroup/ESMValCore/pull/949>`__ (`#978 <https://github.com/ESMValGroup/ESMValCore/pull/978>`__) `Peter Kalverla <https://github.com/Peter9192>`__


Automatic testing
~~~~~~~~~~~~~~~~~

-  Pin cftime>=1.3.0 to have newer string formatting in and fix two tests (`#878 <https://github.com/ESMValGroup/ESMValCore/pull/878>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Switched miniconda conda setup hooks for Github Actions workflows (`#873 <https://github.com/ESMValGroup/ESMValCore/pull/873>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add test for latest version resolver (`#874 <https://github.com/ESMValGroup/ESMValCore/pull/874>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Update codacy coverage reporter to fix coverage (`#905 <https://github.com/ESMValGroup/ESMValCore/pull/905>`__) `Niels Drost <https://github.com/nielsdrost>`__
-  Avoid hardcoded year in tests and add improvement to plev test case (`#921 <https://github.com/ESMValGroup/ESMValCore/pull/921>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Pin scipy to less than 1.6.0 until ESMValGroup/ESMValCore/issues/927 gets resolved (`#928 <https://github.com/ESMValGroup/ESMValCore/pull/928>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Github Actions: change time when conda install test runs (`#930 <https://github.com/ESMValGroup/ESMValCore/pull/930>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Remove redundant test line from test_utils.py (`#935 <https://github.com/ESMValGroup/ESMValCore/pull/935>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Removed netCDF4 package from integration tests of fixes (`#938 <https://github.com/ESMValGroup/ESMValCore/pull/938>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Use new conda environment for installing ESMValCore in Docker containers (`#951 <https://github.com/ESMValGroup/ESMValCore/pull/951>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Notebook API (experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Implement importable config object in experimental API submodule (`#868 <https://github.com/ESMValGroup/ESMValCore/pull/868>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add loading and running recipes to the notebook API (`#907 <https://github.com/ESMValGroup/ESMValCore/pull/907>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add displaying and loading of recipe output to the notebook API (`#957 <https://github.com/ESMValGroup/ESMValCore/pull/957>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add functionality to run single diagnostic task to notebook API (`#962 <https://github.com/ESMValGroup/ESMValCore/pull/962>`__) `Stef Smeets <https://github.com/stefsmeets>`__

Improvements
~~~~~~~~~~~~

-  Create CODEOWNERS file (`#809 <https://github.com/ESMValGroup/ESMValCore/pull/809>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Remove code needed for Python <3.6 (`#844 <https://github.com/ESMValGroup/ESMValCore/pull/844>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add requests as a dependency (`#850 <https://github.com/ESMValGroup/ESMValCore/pull/850>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Pin Python to less than 3.9 (`#870 <https://github.com/ESMValGroup/ESMValCore/pull/870>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Remove numba dependency (`#880 <https://github.com/ESMValGroup/ESMValCore/pull/880>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add Listing and finding recipes to the experimental notebook API (`#901 <https://github.com/ESMValGroup/ESMValCore/pull/901>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Skip variables that don't have dataset or additional_dataset keys (`#860 <https://github.com/ESMValGroup/ESMValCore/pull/860>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Refactor logging configuration (`#933 <https://github.com/ESMValGroup/ESMValCore/pull/933>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Xco2 derivation (`#913 <https://github.com/ESMValGroup/ESMValCore/pull/913>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Working environment for Python 3.9 (pin to !=3.9.0) (`#885 <https://github.com/ESMValGroup/ESMValCore/pull/885>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Print source file when using config get_config_user command (`#960 <https://github.com/ESMValGroup/ESMValCore/pull/960>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Switch to Iris 3 (`#819 <https://github.com/ESMValGroup/ESMValCore/pull/819>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Refactor tasks (`#959 <https://github.com/ESMValGroup/ESMValCore/pull/959>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Restore task summary in debug log after #959 (`#981 <https://github.com/ESMValGroup/ESMValCore/pull/981>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Pin pre-commit hooks (`#974 <https://github.com/ESMValGroup/ESMValCore/pull/974>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Improve error messages when data is missing (`#917 <https://github.com/ESMValGroup/ESMValCore/pull/917>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Set remove_preproc_dir to false in default config-user (`#979 <https://github.com/ESMValGroup/ESMValCore/pull/979>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Move fiona to be installed from conda forge (`#987 <https://github.com/ESMValGroup/ESMValCore/pull/987>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Re-added fiona in setup.py (`#990 <https://github.com/ESMValGroup/ESMValCore/pull/990>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

.. _changelog-v2-1-0:

v2.1.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Set unit=1 if anomalies are standardized (`#727 <https://github.com/ESMValGroup/ESMValCore/pull/727>`__) `bascrezee <https://github.com/bascrezee>`__
-  Fix crash for FGOALS-g2 variables without longitude coordinate (`#729 <https://github.com/ESMValGroup/ESMValCore/pull/729>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve variable alias management (`#595 <https://github.com/ESMValGroup/ESMValCore/pull/595>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix area_statistics fx files loading (`#798 <https://github.com/ESMValGroup/ESMValCore/pull/798>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix units after derivation (`#754 <https://github.com/ESMValGroup/ESMValCore/pull/754>`__) `Manuel Schlund <https://github.com/schlunma>`__

Documentation
~~~~~~~~~~~~~

-  Update v2.0.0 release notes with final additions (`#722 <https://github.com/ESMValGroup/ESMValCore/pull/722>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update package description in setup.py (`#725 <https://github.com/ESMValGroup/ESMValCore/pull/725>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Add installation instructions for pip installation (`#735 <https://github.com/ESMValGroup/ESMValCore/pull/735>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve config-user documentation (`#740 <https://github.com/ESMValGroup/ESMValCore/pull/740>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update the zenodo file with contributors (`#807 <https://github.com/ESMValGroup/ESMValCore/pull/807>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Improve command line run documentation (`#721 <https://github.com/ESMValGroup/ESMValCore/pull/721>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Update the zenodo file with contributors (continued) (`#810 <https://github.com/ESMValGroup/ESMValCore/pull/810>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Reduce size of docker image (`#723 <https://github.com/ESMValGroup/ESMValCore/pull/723>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add 'test' extra to installation, used by docker development tag (`#733 <https://github.com/ESMValGroup/ESMValCore/pull/733>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Correct dockerhub link (`#736 <https://github.com/ESMValGroup/ESMValCore/pull/736>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Create action-install-from-pypi.yml (`#734 <https://github.com/ESMValGroup/ESMValCore/pull/734>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add pre-commit for linting/formatting (`#766 <https://github.com/ESMValGroup/ESMValCore/pull/766>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Run tests in parallel and when building conda package (`#745 <https://github.com/ESMValGroup/ESMValCore/pull/745>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Readable exclude pattern for pre-commit (`#770 <https://github.com/ESMValGroup/ESMValCore/pull/770>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Github Actions Tests (`#732 <https://github.com/ESMValGroup/ESMValCore/pull/732>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Remove isort setup to fix formatting conflict with yapf (`#778 <https://github.com/ESMValGroup/ESMValCore/pull/778>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Fix yapf-isort import formatting conflict (Fixes #777) (`#784 <https://github.com/ESMValGroup/ESMValCore/pull/784>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Sorted output for `esmvaltool recipes list` (`#790 <https://github.com/ESMValGroup/ESMValCore/pull/790>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Replace vmprof with vprof (`#780 <https://github.com/ESMValGroup/ESMValCore/pull/780>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update CMIP6 tables to 6.9.32 (`#706 <https://github.com/ESMValGroup/ESMValCore/pull/706>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Default config-user path now set in config-user read function (`#791 <https://github.com/ESMValGroup/ESMValCore/pull/791>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add custom variable lweGrace (`#692 <https://github.com/ESMValGroup/ESMValCore/pull/692>`__) `bascrezee <https://github.com/bascrezee>`__
- Create Github Actions workflow to build and deploy on Test PyPi and PyPi (`#820 <https://github.com/ESMValGroup/ESMValCore/pull/820>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
- Build and publish the esmvalcore package to conda via Github Actions workflow (`#825 <https://github.com/ESMValGroup/ESMValCore/pull/825>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Fixes for datasets
~~~~~~~~~~~~~~~~~~

-  Fix cmip6 models (`#629 <https://github.com/ESMValGroup/ESMValCore/pull/629>`__) `npgillett <https://github.com/npgillett>`__
-  Fix siconca variable in EC-Earth3 and EC-Earth3-Veg models in amip simulation (`#702 <https://github.com/ESMValGroup/ESMValCore/pull/702>`__) `Evgenia Galytska <https://github.com/egalytska>`__

Preprocessor
~~~~~~~~~~~~

-  Move cmor_check_data to early in preprocessing chain (`#743 <https://github.com/ESMValGroup/ESMValCore/pull/743>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add RMS iris analysis operator to statistics preprocessor functions (`#747 <https://github.com/ESMValGroup/ESMValCore/pull/747>`__) `Pep Cos <https://github.com/pcosbsc>`__
-  Add surface chlorophyll concentration as a derived variable (`#720 <https://github.com/ESMValGroup/ESMValCore/pull/720>`__) `sloosvel <https://github.com/sloosvel>`__
-  Use dask to reduce memory consumption of extract_levels for masked data (`#776 <https://github.com/ESMValGroup/ESMValCore/pull/776>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

.. _changelog-v2-0-0:

v2.0.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Fixed derivation of co2s (`#594 <https://github.com/ESMValGroup/ESMValCore/pull/594>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Padding while cropping needs to stay within sane bounds for shapefiles that span the whole Earth (`#626 <https://github.com/ESMValGroup/ESMValCore/pull/626>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix concatenation of a single cube (`#655 <https://github.com/ESMValGroup/ESMValCore/pull/655>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix mask fx dict handling not to fail if empty list in values (`#661 <https://github.com/ESMValGroup/ESMValCore/pull/661>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Preserve metadata during anomalies computation when using iris cubes difference (`#652 <https://github.com/ESMValGroup/ESMValCore/pull/652>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Avoid crashing when there is directory 'esmvaltool' in the current working directory (`#672 <https://github.com/ESMValGroup/ESMValCore/pull/672>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Solve bug in ACCESS1 dataset fix for calendar.  (`#671 <https://github.com/ESMValGroup/ESMValCore/pull/671>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Fix the syntax for adding multiple ensemble members from the same dataset (`#678 <https://github.com/ESMValGroup/ESMValCore/pull/678>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Fix bug that made preprocessor with fx files fail in rare cases (`#670 <https://github.com/ESMValGroup/ESMValCore/pull/670>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add support for string coordinates (`#657 <https://github.com/ESMValGroup/ESMValCore/pull/657>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fixed the shape extraction to account for wraparound shapefile coords (`#319 <https://github.com/ESMValGroup/ESMValCore/pull/319>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fixed bug in time weights calculation (`#695 <https://github.com/ESMValGroup/ESMValCore/pull/695>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix diagnostic filter (`#713 <https://github.com/ESMValGroup/ESMValCore/pull/713>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Documentation
~~~~~~~~~~~~~

-  Add pandas as a requirement for building the documentation (`#607 <https://github.com/ESMValGroup/ESMValCore/pull/607>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Document default order in which preprocessor functions are applied (`#633 <https://github.com/ESMValGroup/ESMValCore/pull/633>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add pointers about data loading and CF standards to documentation (`#571 <https://github.com/ESMValGroup/ESMValCore/pull/571>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Config file populated with site-specific data paths examples (`#619 <https://github.com/ESMValGroup/ESMValCore/pull/619>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update Codacy badges (`#643 <https://github.com/ESMValGroup/ESMValCore/pull/643>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update copyright info on readthedocs (`#668 <https://github.com/ESMValGroup/ESMValCore/pull/668>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Updated references to documentation (now docs.esmvaltool.org) (`#675 <https://github.com/ESMValGroup/ESMValCore/pull/675>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add all European grants to Zenodo (`#680 <https://github.com/ESMValGroup/ESMValCore/pull/680>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update Sphinx to v3 or later (`#683 <https://github.com/ESMValGroup/ESMValCore/pull/683>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Increase version to 2.0.0 and add release notes (`#691 <https://github.com/ESMValGroup/ESMValCore/pull/691>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update setup.py and README.md for use on PyPI (`#693 <https://github.com/ESMValGroup/ESMValCore/pull/693>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Suggested Documentation changes (`#690 <https://github.com/ESMValGroup/ESMValCore/pull/690>`__) `Steve Smith <https://github.com/ssmithClimate>`__

Improvements
~~~~~~~~~~~~

-  Reduce the size of conda package (`#606 <https://github.com/ESMValGroup/ESMValCore/pull/606>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add a few unit tests for DiagnosticTask (`#613 <https://github.com/ESMValGroup/ESMValCore/pull/613>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Make ncl or R tests not fail if package not installed (`#610 <https://github.com/ESMValGroup/ESMValCore/pull/610>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin flake8<3.8.0 (`#623 <https://github.com/ESMValGroup/ESMValCore/pull/623>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Log warnings for likely errors in provenance record (`#592 <https://github.com/ESMValGroup/ESMValCore/pull/592>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Unpin flake8 (`#646 <https://github.com/ESMValGroup/ESMValCore/pull/646>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  More flexible native6 default DRS (`#645 <https://github.com/ESMValGroup/ESMValCore/pull/645>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Try to use the same python for running diagnostics as for esmvaltool (`#656 <https://github.com/ESMValGroup/ESMValCore/pull/656>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix test for lower python version and add note on lxml (`#659 <https://github.com/ESMValGroup/ESMValCore/pull/659>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Added 1m deep average soil moisture variable (`#664 <https://github.com/ESMValGroup/ESMValCore/pull/664>`__) `bascrezee <https://github.com/bascrezee>`__
-  Update docker recipe (`#603 <https://github.com/ESMValGroup/ESMValCore/pull/603>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Improve command line interface (`#605 <https://github.com/ESMValGroup/ESMValCore/pull/605>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Remove utils directory (`#697 <https://github.com/ESMValGroup/ESMValCore/pull/697>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Avoid pytest version that crashes (`#707 <https://github.com/ESMValGroup/ESMValCore/pull/707>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Options arg in read_config_user_file now optional (`#716 <https://github.com/ESMValGroup/ESMValCore/pull/716>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Produce a readable warning if ancestors are a string instead of a list. (`#711 <https://github.com/ESMValGroup/ESMValCore/pull/711>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Pin Yamale to v2 (`#718 <https://github.com/ESMValGroup/ESMValCore/pull/718>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Expanded cmor public API (`#714 <https://github.com/ESMValGroup/ESMValCore/pull/714>`__) `Manuel Schlund <https://github.com/schlunma>`__

Fixes for datasets
~~~~~~~~~~~~~~~~~~

-  Added various fixes for hybrid height coordinates (`#562 <https://github.com/ESMValGroup/ESMValCore/pull/562>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Extended fix for cl-like variables of CESM2 models (`#604 <https://github.com/ESMValGroup/ESMValCore/pull/604>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Added fix to convert "geopotential" to "geopotential height" for ERA5 (`#640 <https://github.com/ESMValGroup/ESMValCore/pull/640>`__) `Evgenia Galytska <https://github.com/egalytska>`__
-  Do not fix longitude values if they are too far from valid range (`#636 <https://github.com/ESMValGroup/ESMValCore/pull/636>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Preprocessor
~~~~~~~~~~~~

-  Implemented concatenation of cubes with derived coordinates (`#546 <https://github.com/ESMValGroup/ESMValCore/pull/546>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix derived variable ctotal calculation depending on project and standard name (`#620 <https://github.com/ESMValGroup/ESMValCore/pull/620>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  State of the art FX variables handling without preprocessing (`#557 <https://github.com/ESMValGroup/ESMValCore/pull/557>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add max, min and std operators to multimodel (`#602 <https://github.com/ESMValGroup/ESMValCore/pull/602>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Added preprocessor to extract amplitude of cycles (`#597 <https://github.com/ESMValGroup/ESMValCore/pull/597>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Overhaul concatenation and allow for correct concatenation of multiple overlapping datasets (`#615 <https://github.com/ESMValGroup/ESMValCore/pull/615>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Change volume stats to handle and output masked array result (`#618 <https://github.com/ESMValGroup/ESMValCore/pull/618>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Area_weights for cordex in area_statistics (`#631 <https://github.com/ESMValGroup/ESMValCore/pull/631>`__) `mwjury <https://github.com/mwjury>`__
-  Accept cubes as input in multimodel (`#637 <https://github.com/ESMValGroup/ESMValCore/pull/637>`__) `sloosvel <https://github.com/sloosvel>`__
-  Make multimodel work correctly with yearly data (`#677 <https://github.com/ESMValGroup/ESMValCore/pull/677>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Optimize time weights in time preprocessor for climate statistics (`#684 <https://github.com/ESMValGroup/ESMValCore/pull/684>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add percentiles to multi-model stats (`#679 <https://github.com/ESMValGroup/ESMValCore/pull/679>`__) `Peter Kalverla <https://github.com/Peter9192>`__

.. _changelog-v2-0-0b9:

v2.0.0b9
--------

This release includes

Bug fixes
~~~~~~~~~

-  Cast dtype float32 to output from zonal and meridional area preprocessors (`#581 <https://github.com/ESMValGroup/ESMValCore/pull/581>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Unpin on Python<3.8 for conda package (run) (`#570 <https://github.com/ESMValGroup/ESMValCore/pull/570>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update pytest installation marker (`#572 <https://github.com/ESMValGroup/ESMValCore/pull/572>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Remove vmrh2o (`#573 <https://github.com/ESMValGroup/ESMValCore/pull/573>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Restructure documentation (`#575 <https://github.com/ESMValGroup/ESMValCore/pull/575>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix mask in land variables for CCSM4 (`#579 <https://github.com/ESMValGroup/ESMValCore/pull/579>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix derive scripts wrt required method (`#585 <https://github.com/ESMValGroup/ESMValCore/pull/585>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Check coordinates do not have repeated standard names (`#558 <https://github.com/ESMValGroup/ESMValCore/pull/558>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Added derivation script for co2s (`#587 <https://github.com/ESMValGroup/ESMValCore/pull/587>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Adapted custom co2s table to match CMIP6 version (`#588 <https://github.com/ESMValGroup/ESMValCore/pull/588>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Increase version to v2.0.0b9 (`#593 <https://github.com/ESMValGroup/ESMValCore/pull/593>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add a method to save citation information (`#402 <https://github.com/ESMValGroup/ESMValCore/pull/402>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__

For older releases, see the release notes on https://github.com/ESMValGroup/ESMValCore/releases.
.. _preprocessor:

************
Preprocessor
************

In this section, each of the preprocessor modules is described,
roughly following the default order in which preprocessor functions are applied:

* :ref:`Variable derivation`
* :ref:`CMOR check and dataset-specific fixes`
* :ref:`Fx variables as cell measures or ancillary variables`
* :ref:`Vertical interpolation`
* :ref:`Weighting`
* :ref:`Land/Sea/Ice masking`
* :ref:`Horizontal regridding`
* :ref:`Masking of missing values`
* :ref:`Multi-model statistics`
* :ref:`Time operations`
* :ref:`Area operations`
* :ref:`Volume operations`
* :ref:`Cycles`
* :ref:`Trend`
* :ref:`Detrend`
* :ref:`Unit conversion`
* :ref:`Bias`
* :ref:`Other`

See :ref:`preprocessor_functions` for implementation details and the exact default order.

Overview
========

..
   ESMValTool is a modular ``Python 3.7+`` software package possessing capabilities
   of executing a large number of diagnostic routines that can be written in a
   number of programming languages (Python, NCL, R, Julia). The modular nature
   benefits the users and developers in different key areas: a new feature
   developed specifically for version 2.0 is the preprocessing core  or the
   preprocessor (esmvalcore) that executes the bulk of standardized data
   operations and is highly optimized for maximum performance in data-intensive
   tasks. The main objective of the preprocessor is to integrate as many
   standardizable data analysis functions as possible so that the diagnostics can
   focus on the specific scientific tasks they carry. The preprocessor is linked
   to the diagnostics library and the diagnostic execution is seamlessly performed
   after the preprocessor has completed its steps. The benefit of having a
   preprocessing unit separate from the diagnostics library include:

   * ease of integration of new preprocessing routines;
   * ease of maintenance (including unit and integration testing) of existing
     routines;
   * a straightforward manner of importing and using the preprocessing routines as
     part  of the overall usage of the software and, as a special case, the use
     during diagnostic execution;
   * shifting the effort for the scientific diagnostic developer from implementing
     both standard and diagnostic-specific functionalities to allowing them to
     dedicate most of the effort to developing scientifically-relevant diagnostics
     and metrics;
   * a more strict code review process, given the smaller code base than for
     diagnostics.

The ESMValTool preprocessor can be used to perform a broad range of operations
on the input data before diagnostics or metrics are applied. The preprocessor
performs these operations in a centralized, documented and efficient way, thus
reducing the data processing load on the diagnostics side.  For an overview of
the preprocessor structure see the :ref:`Preprocessors`.

Each of the preprocessor operations is written in a dedicated python module and
all of them receive and return an instance of
:obj:`iris.cube.Cube`, working
sequentially on the data with no interactions between them. The order in which
the preprocessor operations is applied is set by default to minimize
the loss of information due to, for example, temporal and spatial subsetting or
multi-model averaging. Nevertheless, the user is free to change such order to
address specific scientific requirements, but keeping in mind that some
operations must be necessarily performed in a specific order. This is the case,
for instance, for multi-model statistics, which required the model to be on a
common grid and therefore has to be called after the regridding module.


.. _Variable derivation:

Variable derivation
===================
The variable derivation module allows to derive variables which are not in the
CMIP standard data request using standard variables as input. The typical use
case of this operation is the evaluation of a variable which is only available
in an observational dataset but not in the models. In this case a derivation
function is provided by the ESMValTool in order to calculate the variable and
perform the comparison. For example, several observational datasets deliver
total column ozone as observed variable (`toz`), but CMIP models only provide
the ozone 3D field. In this case, a derivation function is provided to
vertically integrate the ozone and obtain total column ozone for direct
comparison with the observations.

To contribute a new derived variable, it is also necessary to define a name for
it and to provide the corresponding CMOR table. This is to guarantee the proper
metadata definition is attached to the derived data. Such custom CMOR tables
are collected as part of the `ESMValCore package
<https://github.com/ESMValGroup/ESMValCore>`_. By default, the variable
derivation will be applied only if the variable is not already available in the
input data, but the derivation can be forced by setting the appropriate flag.

.. code-block:: yaml

  variables:
    toz:
      derive: true
      force_derivation: false

The required arguments for this module are two boolean switches:

* ``derive``: activate variable derivation
* ``force_derivation``: force variable derivation even if the variable is
  directly available in the input data.

See also :func:`esmvalcore.preprocessor.derive`. To get an overview on
derivation scripts and how to implement new ones, please go to
:ref:`derivation`.


.. _CMOR check and dataset-specific fixes:

CMORization and dataset-specific fixes
======================================

Data checking
-------------

Data preprocessed by ESMValTool is automatically checked against its
cmor definition. To reduce the impact of this check while maintaining
it as reliable as possible, it is split in two parts: one will check
the metadata and will be done just after loading and concatenating the
data and the other one will check the data itself and will be applied
after all extracting operations are applied to reduce the amount of
data to process.

Checks include, but are not limited to:

   - Requested coordinates are present and comply with their definition.
   - Correctness of variable names, units and other metadata.
   - Compliance with the valid minimum and maximum values allowed if defined.

The most relevant (i.e. a missing coordinate) will raise an error while others
(i.e an incorrect long name) will be reported as a warning.

Some of those issues will be fixed automatically by the tool, including the
following:

    - Incorrect standard or long names.
    - Incorrect units, if they can be converted to the correct ones.
    - Direction of coordinates.
    - Automatic clipping of longitude to 0 - 360 interval.
    - Minute differences between the required and actual vertical coordinate values


Dataset specific fixes
----------------------

Sometimes, the checker will detect errors that it can not fix by itself.
ESMValTool deals with those issues by applying specific fixes for those
datasets that require them. Fixes are applied at three different preprocessor
steps:

    - fix_file: apply fixes directly to a copy of the file. Copying the files
      is costly, so only errors that prevent Iris to load the file are fixed
      here. See :func:`esmvalcore.preprocessor.fix_file`

    - fix_metadata: metadata fixes are done just before concatenating the cubes
      loaded from different files in the final one. Automatic metadata fixes
      are also applied at this step. See
      :func:`esmvalcore.preprocessor.fix_metadata`

    - fix_data: data fixes are applied before starting any operation that will
      alter the data itself. Automatic data fixes are also applied at this step.
      See :func:`esmvalcore.preprocessor.fix_data`

To get an overview on data fixes and how to implement new ones, please go to
:ref:`fixing_data`.

.. _Fx variables as cell measures or ancillary variables:

Fx variables as cell measures or ancillary variables
====================================================
The following preprocessors may require the use of ``fx_variables`` to be able
to perform the computations:

============================================================== =====================
Preprocessor                                                   Default fx variables
============================================================== =====================
:ref:`area_statistics<area_statistics>`                        ``areacella``, ``areacello``
:ref:`mask_landsea<land/sea/ice masking>`                      ``sftlf``, ``sftof``
:ref:`mask_landseaice<ice masking>`                            ``sftgif``
:ref:`volume_statistics<volume_statistics>`                    ``volcello``
:ref:`weighting_landsea_fraction<land/sea fraction weighting>` ``sftlf``, ``sftof``
============================================================== =====================

If the option ``fx_variables`` is not explicitly specified for these
preprocessors, the default fx variables in the second column are automatically
used. If given, the ``fx_variables`` argument specifies the fx variables that
the user wishes to input to the corresponding preprocessor function. The user
may specify these by simply adding the names of the variables, e.g.,

.. code-block:: yaml

    fx_variables:
      areacello:
      volcello:

or by additionally specifying further keys that are used to define the fx
datasets, e.g.,

.. code-block:: yaml

    fx_variables:
      areacello:
        mip: Ofx
        exp: piControl
      volcello:
        mip: Omon

This might be useful to select fx files from a specific ``mip`` table or from a
specific ``exp`` in case not all experiments provide the fx variable.

Alternatively, the ``fx_variables`` argument can also be specified as a list:

.. code-block:: yaml

    fx_variables: ['areacello', 'volcello']

or as a list of dictionaries:

.. code-block:: yaml

    fx_variables: [{'short_name': 'areacello', 'mip': 'Ofx', 'exp': 'piControl'}, {'short_name': 'volcello', 'mip': 'Omon'}]

The recipe parser will automatically find the data files that are associated
with these variables and pass them to the function for loading and processing.

If ``mip`` is not given, ESMValTool will search for the fx variable in all
available tables of the specified project.

.. warning::
   Some fx variables exist in more than one table (e.g., ``volcello`` exists in
   CMIP6's ``Odec``, ``Ofx``, ``Omon``, and ``Oyr`` tables; ``sftgif`` exists
   in CMIP6's ``fx``, ``IyrAnt`` and ``IyrGre``, and ``LImon`` tables). If (for
   a given dataset) fx files are found in more than one table, ``mip`` needs to
   be specified, otherwise an error is raised.

.. note::
   To explicitly **not** use any fx variables in a preprocessor, use
   ``fx_variables: null``.  While some of the preprocessors mentioned above do
   work without fx variables (e.g., ``area_statistics`` or ``mask_landsea``
   with datasets that have regular latitude/longitude grids), using this option
   is **not** recommended.

Internally, the required ``fx_variables`` are automatically loaded by the
preprocessor step ``add_fx_variables`` which also checks them against CMOR
standards and adds them either as ``cell_measure`` (see `CF conventions on cell
measures
<https://cfconventions.org/cf-conventions/cf-conventions.html#cell-measures>`_
and :class:`iris.coords.CellMeasure`) or ``ancillary_variable`` (see `CF
conventions on ancillary variables
<https://cfconventions.org/cf-conventions/cf-conventions.html#ancillary-data>`_
and :class:`iris.coords.AncillaryVariable`) inside the cube data. This ensures
that the defined preprocessor chain is applied to both ``variables`` and
``fx_variables``.

Note that when calling steps that require ``fx_variables`` inside diagnostic
scripts, the variables are expected to contain the required ``cell_measures`` or
``ancillary_variables``. If missing, they can be added using the following functions:

.. code-block::

    from esmvalcore.preprocessor import (add_cell_measure, add_ancillary_variable)

    cube_with_area_measure = add_cell_measure(cube, area_cube, 'area')

    cube_with_volume_measure = add_cell_measure(cube, volume_cube, 'volume)

    cube_with_ancillary_sftlf = add_ancillary_variable(cube, sftlf_cube)

    cube_with_ancillary_sftgif = add_ancillary_variable(cube, sftgif_cube)

  Details on the arguments needed for each step can be found in the following sections.

.. _Vertical interpolation:

Vertical interpolation
======================
Vertical level selection is an important aspect of data preprocessing since it
allows the scientist to perform a number of metrics specific to certain levels
(whether it be air pressure or depth, e.g. the Quasi-Biennial-Oscillation (QBO)
u30 is computed at 30 hPa). Dataset native vertical grids may not come with the
desired set of levels, so an interpolation operation will be needed to regrid
the data vertically. ESMValTool can perform this vertical interpolation via the
``extract_levels`` preprocessor. Level extraction may be done in a number of
ways.

Level extraction can be done at specific values passed to ``extract_levels`` as
``levels:`` with its value a list of levels (note that the units are
CMOR-standard, Pascals (Pa)):

.. code-block:: yaml

    preprocessors:
      preproc_select_levels_from_list:
        extract_levels:
          levels: [100000., 50000., 3000., 1000.]
          scheme: linear

It is also possible to extract the CMIP-specific, CMOR levels as they appear in
the CMOR table, e.g. ``plev10`` or ``plev17`` or ``plev19`` etc:

.. code-block:: yaml

    preprocessors:
      preproc_select_levels_from_cmip_table:
        extract_levels:
          levels: {cmor_table: CMIP6, coordinate: plev10}
          scheme: nearest

Of good use is also the level extraction with values specific to a certain
dataset, without the user actually polling the dataset of interest to find out
the specific levels: e.g. in the example below we offer two alternatives to
extract the levels and vertically regrid onto the vertical levels of
``ERA-Interim``:

.. code-block:: yaml

    preprocessors:
      preproc_select_levels_from_dataset:
        extract_levels:
          levels: ERA-Interim
          # This also works, but allows specifying the pressure coordinate name
          # levels: {dataset: ERA-Interim, coordinate: air_pressure}
          scheme: linear_extrapolate

By default, vertical interpolation is performed in the dimension coordinate of
the z axis. If you want to explicitly declare the z axis coordinate to use
(for example, ``air_pressure``' in variables that are provided in model levels
and not pressure levels) you can override that automatic choice by providing
the name of the desired coordinate:

.. code-block:: yaml

    preprocessors:
      preproc_select_levels_from_dataset:
        extract_levels:
          levels: ERA-Interim
          scheme: linear_extrapolate
          coordinate: air_pressure

If ``coordinate`` is specified, pressure levels (if present) can be converted
to height levels and vice versa using the US standard atmosphere. E.g.
``coordinate = altitude`` will convert existing pressure levels
(air_pressure) to height levels (altitude);
``coordinate = air_pressure`` will convert existing height levels
(altitude) to pressure levels (air_pressure).

If the requested levels are very close to the values in the input data,
the function will just select the available levels instead of interpolating.
The meaning of 'very close' can be changed by providing the parameters:

* ``rtol``
    Relative tolerance for comparing the levels in the input data to the requested
    levels. If the levels are sufficiently close, the requested levels
    will be assigned to the vertical coordinate and no interpolation will take place.
    The default value is 10^-7.
* ``atol``
    Absolute tolerance for comparing the levels in the input data to the requested
    levels. If the levels are sufficiently close, the requested levels
    will be assigned to the vertical coordinate and no interpolation will take place.
    By default, `atol` will be set to 10^-7 times the mean value of
    of the available levels.

Schemes for vertical interpolation and extrapolation
----------------------------------------------------

The vertical interpolation currently supports the following schemes:

* ``linear``: Linear interpolation without extrapolation, i.e., extrapolation
  points will be masked even if the source data is not a masked array.
* ``linear_extrapolate``: Linear interpolation with **nearest-neighbour**
  extrapolation, i.e., extrapolation points will take their value from the
  nearest source point.
* ``nearest``: Nearest-neighbour interpolation without extrapolation, i.e.,
  extrapolation points will be masked even if the source data is not a masked
  array.
* ``nearest_extrapolate``: Nearest-neighbour interpolation with nearest-neighbour
  extrapolation, i.e., extrapolation points will take their value from the
  nearest source point.

.. note::
   Previous versions of ESMValCore (<2.5.0) supported the schemes
   ``linear_horizontal_extrapolate_vertical`` and
   ``nearest_horizontal_extrapolate_vertical``. These schemes have been renamed
   to ``linear_extrapolate`` and ``nearest_extrapolate``, respectively, in
   version 2.5.0 and are identical to the new schemes. Support for the old
   names will be removed in version 2.7.0.

* See also :func:`esmvalcore.preprocessor.extract_levels`.
* See also :func:`esmvalcore.preprocessor.get_cmor_levels`.

.. note::

   Controlling the extrapolation mode allows us to avoid situations where
   extrapolating values makes little physical sense (e.g. extrapolating beyond
   the last data point).


.. _weighting:

Weighting
=========

.. _land/sea fraction weighting:

Land/sea fraction weighting
---------------------------

This preprocessor allows weighting of data by land or sea fractions. In other
words, this function multiplies the given input field by a fraction in the range 0-1 to
account for the fact that not all grid points are completely land- or sea-covered.

The application of this preprocessor is very important for most carbon cycle variables (and
other land surface outputs), which are e.g. reported in units of
:math:`kgC~m^{-2}`. Here, the surface unit actually refers to 'square meter of land/sea' and
NOT 'square meter of gridbox'. In order to integrate these globally or
regionally one has to weight by both the surface quantity and the
land/sea fraction.

For example, to weight an input field with the land fraction, the following
preprocessor can be used:

.. code-block:: yaml

    preprocessors:
      preproc_weighting:
        weighting_landsea_fraction:
          area_type: land
          exclude: ['CanESM2', 'reference_dataset']

Allowed arguments for the keyword ``area_type`` are ``land`` (fraction is 1
for grid cells with only land surface, 0 for grid cells with only sea surface
and values in between 0 and 1 for coastal regions) and ``sea`` (1 for
sea, 0 for land, in between for coastal regions). The optional argument
``exclude`` allows to exclude specific datasets from this preprocessor, which
is for example useful for climate models which do not offer land/sea fraction
files. This arguments also accepts the special dataset specifiers
``reference_dataset`` and ``alternative_dataset``.

Optionally you can specify your own custom fx variable to be used in cases when
e.g. a certain experiment is preferred for fx data retrieval:

.. code-block:: yaml

    preprocessors:
      preproc_weighting:
        weighting_landsea_fraction:
          area_type: land
          exclude: ['CanESM2', 'reference_dataset']
          fx_variables:
            sftlf:
              exp: piControl
            sftof:
              exp: piControl

or alternatively:

.. code-block:: yaml

    preprocessors:
      preproc_weighting:
        weighting_landsea_fraction:
          area_type: land
          exclude: ['CanESM2', 'reference_dataset']
          fx_variables: [
            {'short_name': 'sftlf', 'exp': 'piControl'},
            {'short_name': 'sftof', 'exp': 'piControl'}
            ]

More details on the argument ``fx_variables`` and its default values are given
in :ref:`Fx variables as cell measures or ancillary variables`.

See also :func:`esmvalcore.preprocessor.weighting_landsea_fraction`.


.. _masking:

Masking
=======

Introduction to masking
-----------------------

Certain metrics and diagnostics need to be computed and performed on specific
domains on the globe. The ESMValTool preprocessor supports filtering
the input data on continents, oceans/seas and ice. This is achieved by masking
the model data and keeping only the values associated with grid points that
correspond to, e.g., land, ocean or ice surfaces, as specified by the
user. Where possible, the masking is realized using the standard mask files
provided together with the model data as part of the CMIP data request (the
so-called fx variable). In the absence of these files, the Natural Earth masks
are used: although these are not model-specific, they represent a good
approximation since they have a much higher resolution than most of the models
and they are regularly updated with changing geographical features.

.. _land/sea/ice masking:

Land-sea masking
----------------

In ESMValTool, land-sea-ice masking can be done in two places: in the
preprocessor, to apply a mask on the data before any subsequent preprocessing
step and before running the diagnostic, or in the diagnostic scripts
themselves. We present both these implementations below.

To mask out a certain domain (e.g., sea) in the preprocessor,
``mask_landsea`` can be used:

.. code-block:: yaml

    preprocessors:
      preproc_mask:
        mask_landsea:
          mask_out: sea

and requires only one argument: ``mask_out``: either ``land`` or ``sea``.

Optionally you can specify your own custom fx variable to be used in cases when e.g. a certain
experiment is preferred for fx data retrieval. Note that it is possible to specify as many tags
for the fx variable as required:


.. code-block:: yaml

    preprocessors:
      landmask:
        mask_landsea:
          mask_out: sea
          fx_variables:
            sftlf:
              exp: piControl
            sftof:
              exp: piControl
              ensemble: r2i1p1f1

or alternatively:

.. code-block:: yaml

    preprocessors:
      landmask:
        mask_landsea:
          mask_out: sea
          fx_variables: [
            {'short_name': 'sftlf', 'exp': 'piControl'},
            {'short_name': 'sftof', 'exp': 'piControl', 'ensemble': 'r2i1p1f1'}
            ]

More details on the argument ``fx_variables`` and its default values are given
in :ref:`Fx variables as cell measures or ancillary variables`.

If the corresponding fx file is not found (which is
the case for some models and almost all observational datasets), the
preprocessor attempts to mask the data using Natural Earth mask files (that are
vectorized rasters). As mentioned above, the spatial resolution of the the
Natural Earth masks are much higher than any typical global model (10m for
land and glaciated areas and 50m for ocean masks).

See also :func:`esmvalcore.preprocessor.mask_landsea`.

.. _ice masking:

Ice masking
-----------

Note that for masking out ice sheets, the preprocessor uses a different
function, to ensure that both land and sea or ice can be masked out without
losing generality. To mask ice out, ``mask_landseaice`` can be used:

.. code-block:: yaml

  preprocessors:
    preproc_mask:
      mask_landseaice:
        mask_out: ice

and requires only one argument: ``mask_out``: either ``landsea`` or ``ice``.

Optionally you can specify your own custom fx variable to be used in cases when
e.g. a certain experiment is preferred for fx data retrieval:


.. code-block:: yaml

    preprocessors:
      landseaicemask:
        mask_landseaice:
          mask_out: sea
          fx_variables:
            sftgif:
              exp: piControl

or alternatively:

.. code-block:: yaml

    preprocessors:
      landseaicemask:
        mask_landseaice:
          mask_out: sea
          fx_variables: [{'short_name': 'sftgif', 'exp': 'piControl'}]

More details on the argument ``fx_variables`` and its default values are given
in :ref:`Fx variables as cell measures or ancillary variables`.

See also :func:`esmvalcore.preprocessor.mask_landseaice`.

Glaciated masking
-----------------

For masking out glaciated areas a Natural Earth shapefile is used. To mask
glaciated areas out, ``mask_glaciated`` can be used:

.. code-block:: yaml

  preprocessors:
    preproc_mask:
      mask_glaciated:
        mask_out: glaciated

and it requires only one argument: ``mask_out``: only ``glaciated``.

See also :func:`esmvalcore.preprocessor.mask_landseaice`.

.. _masking of missing values:

Missing values masks
--------------------

Missing (masked) values can be a nuisance especially when dealing with
multi-model ensembles and having to compute multi-model statistics; different
numbers of missing data from dataset to dataset may introduce biases and
artificially assign more weight to the datasets that have less missing data.
This is handled in ESMValTool via the missing values masks: two types of such
masks are available, one for the multi-model case and another for the single
model case.

The multi-model missing values mask (``mask_fillvalues``) is a preprocessor step
that usually comes after all the single-model steps (regridding, area selection
etc) have been performed; in a nutshell, it combines missing values masks from
individual models into a multi-model missing values mask; the individual model
masks are built according to common criteria: the user chooses a time window in
which missing data points are counted, and if the number of missing data points
relative to the number of total data points in a window is less than a chosen
fractional threshold, the window is discarded i.e. all the points in the window
are masked (set to missing).

.. code-block:: yaml

    preprocessors:
      missing_values_preprocessor:
        mask_fillvalues:
          threshold_fraction: 0.95
          min_value: 19.0
          time_window: 10.0

In the example above, the fractional threshold for missing data vs. total data
is set to 95% and the time window is set to 10.0 (units of the time coordinate
units). Optionally, a minimum value threshold can be applied, in this case it
is set to 19.0 (in units of the variable units).

See also :func:`esmvalcore.preprocessor.mask_fillvalues`.

Common mask for multiple models
-------------------------------

To create a combined multi-model mask (all the masks from all the analyzed
datasets combined into a single mask using a logical OR), the preprocessor
``mask_multimodel`` can be used. In contrast to ``mask_fillvalues``,
``mask_multimodel`` does not expect that the datasets have a ``time``
coordinate, but works on datasets with arbitrary (but identical) coordinates.
After ``mask_multimodel``, all involved datasets have an identical mask.

See also :func:`esmvalcore.preprocessor.mask_multimodel`.

Minimum, maximum and interval masking
-------------------------------------

Thresholding on minimum and maximum accepted data values can also be performed:
masks are constructed based on the results of thresholding; inside and outside
interval thresholding and masking can also be performed. These functions are
``mask_above_threshold``, ``mask_below_threshold``, ``mask_inside_range``, and
``mask_outside_range``.

These functions always take a cube as first argument and either ``threshold``
for threshold masking or the pair ``minimum``, ``maximum`` for interval masking.

See also :func:`esmvalcore.preprocessor.mask_above_threshold` and related
functions.


.. _Horizontal regridding:

Horizontal regridding
=====================

Regridding is necessary when various datasets are available on a variety of
`lat-lon` grids and they need to be brought together on a common grid (for
various statistical operations e.g. multi-model statistics or for e.g. direct
inter-comparison or comparison with observational datasets). Regridding is
conceptually a very similar process to interpolation (in fact, the regridder
engine uses interpolation and extrapolation, with various schemes). The primary
difference is that interpolation is based on sample data points, while
regridding is based on the horizontal grid of another cube (the reference
grid). If the horizontal grids of a cube and its reference grid are sufficiently
the same, regridding is automatically and silently skipped for performance reasons.

The underlying regridding mechanism in ESMValTool uses
:obj:`iris.cube.Cube.regrid`
from Iris.

The use of the horizontal regridding functionality is flexible depending on
what type of reference grid and what interpolation scheme is preferred. Below
we show a few examples.

Regridding on a reference dataset grid
--------------------------------------

The example below shows how to regrid on the reference dataset
``ERA-Interim`` (observational data, but just as well CMIP, obs4MIPs,
or ana4mips datasets can be used); in this case the `scheme` is
`linear`.

.. code-block:: yaml

    preprocessors:
      regrid_preprocessor:
        regrid:
          target_grid: ERA-Interim
          scheme: linear

Regridding on an ``MxN`` grid specification
-------------------------------------------

The example below shows how to regrid on a reference grid with a cell
specification of ``2.5x2.5`` degrees. This is similar to regridding on
reference datasets, but in the previous case the reference dataset grid cell
specifications are not necessarily known a priori. Regridding on an ``MxN``
cell specification is oftentimes used when operating on localized data.

.. code-block:: yaml

    preprocessors:
      regrid_preprocessor:
        regrid:
          target_grid: 2.5x2.5
          scheme: nearest

In this case the ``NearestNeighbour`` interpolation scheme is used (see below
for scheme definitions).

When using a ``MxN`` type of grid it is possible to offset the grid cell
centrepoints using the `lat_offset` and ``lon_offset`` arguments:

* ``lat_offset``: offsets the grid centers of the latitude coordinate w.r.t. the
  pole by half a grid step;
* ``lon_offset``: offsets the grid centers of the longitude coordinate
  w.r.t. Greenwich meridian by half a grid step.

.. code-block:: yaml

    preprocessors:
      regrid_preprocessor:
        regrid:
          target_grid: 2.5x2.5
          lon_offset: True
          lat_offset: True
          scheme: nearest

Regridding to a regional target grid specification
--------------------------------------------------

This example shows how to regrid to a regional target grid specification.
This is useful if both a ``regrid`` and ``extract_region`` step are necessary.

.. code-block:: yaml

    preprocessors:
      regrid_preprocessor:
        regrid:
          target_grid:
            start_longitude: 40
            end_longitude: 60
            step_longitude: 2
            start_latitude: -10
            end_latitude: 30
            step_latitude: 2
          scheme: nearest

This defines a grid ranging from 40° to 60° longitude with 2° steps,
and -10° to 30° latitude with 2° steps. If ``end_longitude`` or ``end_latitude`` do
not fall on the grid (e.g., ``end_longitude: 61``), it cuts off at the nearest
previous value (e.g. ``60``).

The longitude coordinates will wrap around the globe if necessary, i.e.
``start_longitude: 350``, ``end_longitude: 370`` is valid input.

The arguments are defined below:

* ``start_latitude``: Latitude value of the first grid cell center (start point).
  The grid includes this value.
* ``end_latitude``: Latitude value of the last grid cell center (end point).
  The grid includes this value only if it falls on a grid point.
  Otherwise, it cuts off at the previous value.
* ``step_latitude``: Latitude distance between the centers of two neighbouring cells.
* ``start_longitude``: Latitude value of the first grid cell center (start point).
  The grid includes this value.
* ``end_longitude``: Longitude value of the last grid cell center (end point).
  The grid includes this value only if it falls on a grid point.
  Otherwise, it cuts off at the previous value.
* ``step_longitude``: Longitude distance between the centers of two neighbouring cells.

Regridding (interpolation, extrapolation) schemes
-------------------------------------------------

The schemes used for the interpolation and extrapolation operations needed by
the horizontal regridding functionality directly map to their corresponding
implementations in :mod:`iris`:

* ``linear``: Linear interpolation without extrapolation, i.e., extrapolation
  points will be masked even if the source data is not a masked array (uses
  ``Linear(extrapolation_mode='mask')``, see :obj:`iris.analysis.Linear`).
* ``linear_extrapolate``: Linear interpolation with extrapolation, i.e.,
  extrapolation points will be calculated by extending the gradient of the
  closest two points (uses ``Linear(extrapolation_mode='extrapolate')``, see
  :obj:`iris.analysis.Linear`).
* ``nearest``: Nearest-neighbour interpolation without extrapolation, i.e.,
  extrapolation points will be masked even if the source data is not a masked
  array (uses ``Nearest(extrapolation_mode='mask')``, see
  :obj:`iris.analysis.Nearest`).
* ``area_weighted``: Area-weighted regridding (uses ``AreaWeighted()``, see
  :obj:`iris.analysis.AreaWeighted`).
* ``unstructured_nearest``: Nearest-neighbour interpolation for unstructured
  grids (uses ``UnstructuredNearest()``, see
  :obj:`iris.analysis.UnstructuredNearest`).

See also :func:`esmvalcore.preprocessor.regrid`

.. note::

   Controlling the extrapolation mode allows us to avoid situations where
   extrapolating values makes little physical sense (e.g. extrapolating beyond
   the last data point).

.. note::

   The regridding mechanism is (at the moment) done with fully realized data in
   memory, so depending on how fine the target grid is, it may use a rather
   large amount of memory. Empirically target grids of up to ``0.5x0.5``
   degrees should not produce any memory-related issues, but be advised that
   for resolutions of ``< 0.5`` degrees the regridding becomes very slow and
   will use a lot of memory.


.. _multi-model statistics:

Multi-model statistics
======================
Computing multi-model statistics is an integral part of model analysis and
evaluation: individual models display a variety of biases depending on model
set-up, initial conditions, forcings and implementation; comparing model data to
observational data, these biases have a significantly lower statistical impact
when using a multi-model ensemble. ESMValTool has the capability of computing a
number of multi-model statistical measures: using the preprocessor module
``multi_model_statistics`` will enable the user to ask for either a multi-model
``mean``, ``median``, ``max``, ``min``, ``std``, and / or ``pXX.YY`` with a set
of argument parameters passed to ``multi_model_statistics``. Percentiles can be
specified like ``p1.5`` or ``p95``. The decimal point will be replaced by a dash
in the output file.

Restrictive computation is also available by excluding  any set of models that
the user will not want to include in the statistics (by setting ``exclude:
[excluded models list]`` argument). The implementation has a few restrictions
that apply to the input data: model datasets must have consistent shapes, apart
from the time dimension; and cubes with more than four dimensions (time,
vertical axis, two horizontal axes) are not supported.

Input datasets may have different time coordinates. Statistics can be computed
across overlapping times only (``span: overlap``) or across the full time span
of the combined models (``span: full``). The preprocessor sets a common time
coordinate on all datasets. As the number of days in a year may vary between
calendars, (sub-)daily data with different calendars are not supported.

Input datasets may have different time coordinates. The multi-model statistics
preprocessor sets a common time coordinate on all datasets. As the number of
days in a year may vary between calendars, (sub-)daily data are not supported.

.. code-block:: yaml

    preprocessors:
      multi_model_preprocessor:
        multi_model_statistics:
          span: overlap
          statistics: [mean, median]
          exclude: [NCEP]

see also :func:`esmvalcore.preprocessor.multi_model_statistics`.

When calling the module inside diagnostic scripts, the input must be given
as a list of cubes. The output will be saved in a dictionary where each
entry contains the resulting cube with the requested statistic operations.

.. code-block::

    from esmvalcore.preprocessor import multi_model_statistics
    statistics = multi_model_statistics([cube1,...,cubeN], 'overlap', ['mean', 'median'])
    mean_cube = statistics['mean']
    median_cube = statistics['median']

.. note::

   The multi-model array operations can be rather memory-intensive (since they
   are not performed lazily as yet). The Section on :ref:`Memory use` details
   the memory intake for different run scenarios, but as a thumb rule, for the
   multi-model preprocessor, the expected maximum memory intake could be
   approximated as the number of datasets multiplied by the average size in
   memory for one dataset.

.. _time operations:

Time manipulation
=================
The ``_time.py`` module contains the following preprocessor functions:

* extract_time_: Extract a time range from a cube.
* extract_season_: Extract only the times that occur within a specific season.
* extract_month_: Extract only the times that occur within a specific month.
* hourly_statistics_: Compute intra-day statistics
* daily_statistics_: Compute statistics for each day
* monthly_statistics_: Compute statistics for each month
* seasonal_statistics_: Compute statistics for each season
* annual_statistics_: Compute statistics for each year
* decadal_statistics_: Compute statistics for each decade
* climate_statistics_: Compute statistics for the full period
* resample_time_: Resample data
* resample_hours_: Convert between N-hourly frequencies by resampling
* anomalies_: Compute (standardized) anomalies
* regrid_time_: Aligns the time axis of each dataset to have common time
  points and calendars.
* timeseries_filter_: Allows application of a filter to the time-series data.

Statistics functions are applied by default in the order they appear in the
list. For example, the following example applied to hourly data will retrieve
the minimum values for the full period (by season) of the monthly mean of the
daily maximum of any given variable.

.. code-block:: yaml

    daily_statistics:
      operator: max

    monthly_statistics:
      operator: mean

    climate_statistics:
      operator: min
      period: season


.. _extract_time:

``extract_time``
----------------

This function subsets a dataset between two points in times. It removes all
times in the dataset before the first time and after the last time point.
The required arguments are relatively self explanatory:

* ``start_year``
* ``start_month``
* ``start_day``
* ``end_year``
* ``end_month``
* ``end_day``

These start and end points are set using the datasets native calendar.
All six arguments should be given as integers - the named month string
will not be accepted.

See also :func:`esmvalcore.preprocessor.extract_time`.

.. _extract_season:

``extract_season``
------------------

Extract only the times that occur within a specific season.

This function only has one argument: ``season``. This is the named season to
extract, i.e. DJF, MAM, JJA, SON, but also all other sequentially correct
combinations, e.g. JJAS.

Note that this function does not change the time resolution. If your original
data is in monthly time resolution, then this function will return three
monthly datapoints per year.

If you want the seasonal average, then this function needs to be combined with
the seasonal_mean function, below.

See also :func:`esmvalcore.preprocessor.extract_season`.

.. _extract_month:

``extract_month``
-----------------

The function extracts the times that occur within a specific month.
This function only has one argument: ``month``. This value should be an integer
between 1 and 12 as the named month string will not be accepted.

See also :func:`esmvalcore.preprocessor.extract_month`.

.. _hourly_statistics:

``hourly_statistics``
---------------------

This function produces statistics at a x-hourly frequency.

Parameters:
    * every_n_hours: frequency to use to compute the statistics. Must be a divisor of
      24.

    * operator: operation to apply. Accepted values are 'mean',
      'median', 'std_dev', 'min', 'max' and 'sum'. Default is 'mean'

See also :func:`esmvalcore.preprocessor.daily_statistics`.

.. _daily_statistics:

``daily_statistics``
--------------------

This function produces statistics for each day in the dataset.

Parameters:
    * operator: operation to apply. Accepted values are 'mean',
      'median', 'std_dev', 'min', 'max', 'sum' and 'rms'. Default is 'mean'

See also :func:`esmvalcore.preprocessor.daily_statistics`.

.. _monthly_statistics:

``monthly_statistics``
----------------------

This function produces statistics for each month in the dataset.

Parameters:
    * operator: operation to apply. Accepted values are 'mean',
      'median', 'std_dev', 'min', 'max', 'sum' and 'rms'. Default is 'mean'

See also :func:`esmvalcore.preprocessor.monthly_statistics`.

.. _seasonal_statistics:

``seasonal_statistics``
-----------------------

This function produces statistics for each season (default: ``[DJF, MAM, JJA,
SON]`` or custom seasons e.g. ``[JJAS, ONDJFMAM]``) in the dataset. Note that
this function will not check for missing time points. For instance, if you are
looking at the DJF field, but your datasets starts on January 1st, the first
DJF field will only contain data from January and February.

We recommend using the extract_time to start the dataset from the following
December and remove such biased initial datapoints.

Parameters:
    * operator: operation to apply. Accepted values are 'mean',
      'median', 'std_dev', 'min', 'max', 'sum' and 'rms'. Default is 'mean'

    * seasons: seasons to build statistics.
      Default is '[DJF, MAM, JJA, SON]'

See also :func:`esmvalcore.preprocessor.seasonal_statistics`.

.. _annual_statistics:

``annual_statistics``
---------------------

This function produces statistics for each year.

Parameters:
    * operator: operation to apply. Accepted values are 'mean',
      'median', 'std_dev', 'min', 'max', 'sum' and 'rms'. Default is 'mean'

See also :func:`esmvalcore.preprocessor.annual_statistics`.

.. _decadal_statistics:

``decadal_statistics``
----------------------

This function produces statistics for each decade.

Parameters:
    * operator: operation to apply. Accepted values are 'mean',
      'median', 'std_dev', 'min', 'max', 'sum' and 'rms'. Default is 'mean'

See also :func:`esmvalcore.preprocessor.decadal_statistics`.

.. _climate_statistics:

``climate_statistics``
----------------------

This function produces statistics for the whole dataset. It can produce scalars
(if the full period is chosen) or daily, monthly or seasonal statistics.

Parameters:
    * operator: operation to apply. Accepted values are 'mean', 'median',
      'std_dev', 'min', 'max', 'sum' and 'rms'. Default is 'mean'

    * period: define the granularity of the statistics: get values for the
      full period, for each month or day of year.
      Available periods: 'full', 'season', 'seasonal', 'monthly', 'month',
      'mon', 'daily', 'day'. Default is 'full'

    * seasons: if period 'seasonal' or 'season' allows to set custom seasons.
      Default is '[DJF, MAM, JJA, SON]'

Examples:
    * Monthly climatology:

        .. code-block:: yaml

            climate_statistics:
                operator: mean
                period: month

    * Daily maximum for the full period:

        .. code-block:: yaml

            climate_statistics:
              operator: max
              period: day

    * Minimum value in the period:

        .. code-block:: yaml

            climate_statistics:
              operator: min
              period: full

See also :func:`esmvalcore.preprocessor.climate_statistics`.

.. _resample_time:

``resample_time``
-----------------

This function changes the frequency of the data in the cube by extracting the
timesteps that meet the criteria. It is important to note that it is mainly
meant to be used with instantaneous data.

Parameters:
    * month: Extract only timesteps from the given month or do nothing if None.
      Default is `None`
    * day: Extract only timesteps from the given day of month or do nothing if
      None. Default is `None`
    * hour: Extract only timesteps from the given hour or do nothing if None.
      Default is `None`

Examples:
    * Hourly data to daily:

        .. code-block:: yaml

            resample_time:
              hour: 12

    * Hourly data to monthly:

        .. code-block:: yaml

            resample_time:
              hour: 12
              day: 15

    * Daily data to monthly:

        .. code-block:: yaml

            resample_time:
              day: 15

See also :func:`esmvalcore.preprocessor.resample_time`.


resample_hours:

``resample_hours``
------------------

This function changes the frequency of the data in the cube by extracting the
timesteps that belongs to the desired frequency. It is important to note that
it is mainly mean to be used with instantaneous data

Parameters:
    * interval: New frequency of the data. Must be a divisor of 24
    * offset: First desired hour. Default 0. Must be lower than the interval

Examples:
    * Convert to 12-hourly, by getting timesteps at 0:00 and 12:00:

        .. code-block:: yaml

            resample_hours:
              hours: 12

    * Convert to 12-hourly, by getting timesteps at 6:00 and 18:00:

        .. code-block:: yaml

            resample_hours:
              hours: 12
	      offset: 6

See also :func:`esmvalcore.preprocessor.resample_hours`.

.. _anomalies:

``anomalies``
----------------------

This function computes the anomalies for the whole dataset. It can compute
anomalies from the full, seasonal, monthly and daily climatologies. Optionally
standardized anomalies can be calculated.

Parameters:
    * period: define the granularity of the climatology to use:
      full period, seasonal, monthly or daily.
      Available periods: 'full', 'season', 'seasonal', 'monthly', 'month',
      'mon', 'daily', 'day'. Default is 'full'
    * reference: Time slice to use as the reference to compute the climatology
      on. Can be 'null' to use the full cube or a dictionary with the
      parameters from extract_time_. Default is null
    * standardize: if true calculate standardized anomalies (default: false)
    * seasons: if period 'seasonal' or 'season' allows to set custom seasons.
      Default is '[DJF, MAM, JJA, SON]'
Examples:
    * Anomalies from the full period climatology:

        .. code-block:: yaml

            anomalies:

    * Anomalies from the full period monthly climatology:

        .. code-block:: yaml

            anomalies:
              period: month

    * Standardized anomalies from the full period climatology:

        .. code-block:: yaml

            anomalies:
              standardized: true


     * Standardized Anomalies from the 1979-2000 monthly climatology:

        .. code-block:: yaml

            anomalies:
              period: month
              reference:
                start_year: 1979
                start_month: 1
                start_day: 1
                end_year: 2000
                end_month: 12
                end_day: 31
              standardize: true

See also :func:`esmvalcore.preprocessor.anomalies`.


.. _regrid_time:

``regrid_time``
---------------

This function aligns the time points of each component dataset so that the Iris
cubes from different datasets can be subtracted. The operation makes the
datasets time points common; it also resets the time
bounds and auxiliary coordinates to reflect the artificially shifted time
points. Current implementation for monthly and daily data; the ``frequency`` is
set automatically from the variable CMOR table unless a custom ``frequency`` is
set manually by the user in recipe.

See also :func:`esmvalcore.preprocessor.regrid_time`.


.. _timeseries_filter:

``timeseries_filter``
---------------------

This function allows the user to apply a filter to the timeseries data. This filter may be
of the user's choice (currently only the ``low-pass`` Lanczos filter is implemented); the
implementation is inspired by this `iris example
<https://scitools-iris.readthedocs.io/en/latest/generated/gallery/general/plot_SOI_filtering.html>`_ and uses aggregation via :obj:`iris.cube.Cube.rolling_window`.

Parameters:
    * window: the length of the filter window (in units of cube time coordinate).
    * span: period (number of months/days, depending on data frequency) on which
      weights should be computed e.g. for 2-yearly: span = 24 (2 x 12 months).
      Make sure span has the same units as the data cube time coordinate.
    * filter_type: the type of filter to be applied; default 'lowpass'.
      Available types: 'lowpass'.
    * filter_stats: the type of statistic to aggregate on the rolling window;
      default 'sum'. Available operators: 'mean', 'median', 'std_dev', 'sum', 'min', 'max', 'rms'.

Examples:
    * Lowpass filter with a monthly mean as operator:

        .. code-block:: yaml

            timeseries_filter:
                window: 3  # 3-monthly filter window
                span: 12   # weights computed on the first year
                filter_type: lowpass  # low-pass filter
                filter_stats: mean    # 3-monthly mean lowpass filter

See also :func:`esmvalcore.preprocessor.timeseries_filter`.

.. _area operations:

Area manipulation
=================
The area manipulation module contains the following preprocessor functions:

* extract_region_: Extract a region from a cube based on ``lat/lon``
  corners.
* extract_named_regions_: Extract a specific region from in the region
  coordinate.
* extract_shape_: Extract a region defined by a shapefile.
* extract_point_: Extract a single point (with interpolation)
* extract_location_: Extract a single point by its location (with interpolation)
* zonal_statistics_: Compute zonal statistics.
* meridional_statistics_: Compute meridional statistics.
* area_statistics_: Compute area statistics.


``extract_region``
------------------

This function returns a subset of the data on the rectangular region requested.
The boundaries of the region are provided as latitude and longitude coordinates
in the arguments:

* ``start_longitude``
* ``end_longitude``
* ``start_latitude``
* ``end_latitude``

Note that this function can only be used to extract a rectangular region. Use
``extract_shape`` to extract any other shaped region from a shapefile.

If the grid is irregular, the returned region retains the original coordinates,
but is cropped to a rectangular bounding box defined by the start/end
coordinates. The deselected area inside the region is masked.

See also :func:`esmvalcore.preprocessor.extract_region`.


``extract_named_regions``
-------------------------

This function extracts a specific named region from the data. This function
takes the following argument: ``regions`` which is either a string or a list
of strings of named regions. Note that the dataset must have a ``region``
coordinate which includes a list of strings as values. This function then
matches the named regions against the requested string.

See also :func:`esmvalcore.preprocessor.extract_named_regions`.


``extract_shape``
-------------------------

Extract a shape or a representative point for this shape from
the data.

Parameters:
  * ``shapefile``: path to the shapefile containing the geometry of the
    region to be extracted. If the file contains multiple shapes behaviour
    depends on the decomposed parameter. This path can be relative to
    ``auxiliary_data_dir`` defined in the :ref:`user configuration file`.
  * ``method``: the method to select the region, selecting either all points
	  contained by the shape or a single representative point. Choose either
	  'contains' or 'representative'. If not a single grid point is contained
	  in the shape, a representative point will be selected.
  * ``crop``: by default extract_region_ will be used to crop the data to a
	  minimal rectangular region containing the shape. Set to ``false`` to only
	  mask data outside the shape. Data on irregular grids will not be cropped.
  * ``decomposed``: by default ``false``, in this case the union of all the
    regions in the shape file is masked out. If ``true``, the regions in the
    shapefiles are masked out separately, generating an auxiliary dimension
    for the cube for this.
  * ``ids``: by default, ``[]``, in this case all the shapes in the file will
    be used. If a list of IDs is provided, only the shapes matching them will
    be used. The IDs are assigned from the ``name`` or ``id`` attributes (in
    that order of priority) if present in the file or from the reading order
    if otherwise not present. So, for example, if a file has both ```name``
    and ``id`` attributes, the ids will be assigned from ``name``. If the file
    only has the ``id`` attribute, it will be taken from it and if no ``name``
    nor ``id`` attributes are present, an integer id starting from 1 will be
    assigned automatically when reading the shapes. We discourage to rely on
    this last behaviour as we can not assure that the reading order will be the
    same in different platforms, so we encourage you to modify the file to add
    a proper id attribute. If the file has an id attribute with a name that is
    not supported, please open an issue so we can add support for it.

Examples:
    * Extract the shape of the river Elbe from a shapefile:

        .. code-block:: yaml

            extract_shape:
              shapefile: Elbe.shp
              method: contains

    * Extract the shape of several countries:

        .. code-block:: yaml

            extract_shape:
            shapefile: NaturalEarth/Countries/ne_110m_admin_0_countries.shp
            decomposed: True
            method: contains
            ids:
              - Spain
              - France
              - Italy
              - United Kingdom
              - Taiwan

See also :func:`esmvalcore.preprocessor.extract_shape`.


``extract_point``
-----------------

Extract a single point from the data. This is done using either
nearest or linear interpolation.

Returns a cube with the extracted point(s), and with adjusted latitude
and longitude coordinates (see below).

Multiple points can also be extracted, by supplying an array of
latitude and/or longitude coordinates. The resulting point cube will
match the respective latitude and longitude coordinate to those of the
input coordinates. If the input coordinate is a scalar, the dimension
will be missing in the output cube (that is, it will be a scalar).

Parameters:
  * ``cube``: the input dataset cube.
  * ``latitude``, ``longitude``: coordinates (as floating point
    values) of the point to be extracted. Either (or both) can also
    be an array of floating point values.
  * ``scheme``: interpolation scheme: either ``'linear'`` or
    ``'nearest'``. There is no default.
    
See also :func:`esmvalcore.preprocessor.extract_point`.


``extract_location``
--------------------

Extract a single point using a location name, with interpolation
(either linear or nearest). This preprocessor extracts a single
location point from a cube, according to the given interpolation
scheme ``scheme``. The function retrieves the coordinates of the
location and then calls the :func:`esmvalcore.preprocessor.extract_point`
preprocessor. It can be used to locate cities and villages,
but also mountains or other geographical locations.

.. note::
   Note that this function's geolocator application needs a
   working internet connection.

Parameters
  * ``cube``: the input dataset cube to extract a point from.
  * ``location``: the reference location. Examples: 'mount everest',
    'romania', 'new york, usa'. Raises ValueError if none supplied.
  * ``scheme`` : interpolation scheme. ``'linear'`` or ``'nearest'``.
    There is no default, raises ValueError if none supplied.

See also :func:`esmvalcore.preprocessor.extract_location`.


``zonal_statistics``
--------------------

The function calculates the zonal statistics by applying an operator
along the longitude coordinate. This function takes one argument:

* ``operator``: Which operation to apply: mean, std_dev, median, min, max, sum or rms.

See also :func:`esmvalcore.preprocessor.zonal_means`.


``meridional_statistics``
-------------------------

The function calculates the meridional statistics by applying an
operator along the latitude coordinate. This function takes one
argument:

* ``operator``: Which operation to apply: mean, std_dev, median, min, max, sum or rms.

See also :func:`esmvalcore.preprocessor.meridional_means`.


.. _area_statistics:

``area_statistics``
-------------------

This function calculates the average value over a region - weighted by the cell
areas of the region. This function takes the argument, ``operator``: the name
of the operation to apply.

This function can be used to apply several different operations in the
horizontal plane: mean, standard deviation, median, variance, minimum, maximum and root mean square.

Note that this function is applied over the entire dataset. If only a specific
region, depth layer or time period is required, then those regions need to be
removed using other preprocessor operations in advance.

The optional ``fx_variables`` argument specifies the fx variables that the user
wishes to input to the function. More details on this are given in :ref:`Fx
variables as cell measures or ancillary variables`.

See also :func:`esmvalcore.preprocessor.area_statistics`.


.. _volume operations:

Volume manipulation
===================
The ``_volume.py`` module contains the following preprocessor functions:

* ``extract_volume``: Extract a specific depth range from a cube.
* ``volume_statistics``: Calculate the volume-weighted average.
* ``depth_integration``: Integrate over the depth dimension.
* ``extract_transect``: Extract data along a line of constant latitude or
  longitude.
* ``extract_trajectory``: Extract data along a specified trajectory.


``extract_volume``
------------------

Extract a specific range in the `z`-direction from a cube.  This function
takes two arguments, a minimum and a maximum (``z_min`` and ``z_max``,
respectively) in the `z`-direction.

Note that this requires the requested `z`-coordinate range to be the same sign
as the Iris cube. That is, if the cube has `z`-coordinate as negative, then
``z_min`` and ``z_max`` need to be negative numbers.

See also :func:`esmvalcore.preprocessor.extract_volume`.


.. _volume_statistics:

``volume_statistics``
---------------------

This function calculates the volume-weighted average across three dimensions,
but maintains the time dimension.

This function takes the argument: ``operator``, which defines the operation to
apply over the volume.

No depth coordinate is required as this is determined by Iris. This function
works best when the ``fx_variables`` provide the cell volume. The optional
``fx_variables`` argument specifies the fx variables that the user wishes to
input to the function. More details on this are given in :ref:`Fx variables as
cell measures or ancillary variables`.

See also :func:`esmvalcore.preprocessor.volume_statistics`.


``depth_integration``
---------------------

This function integrates over the depth dimension. This function does a
weighted sum along the `z`-coordinate, and removes the `z` direction of the
output cube. This preprocessor takes no arguments.

See also :func:`esmvalcore.preprocessor.depth_integration`.


``extract_transect``
--------------------

This function extracts data along a line of constant latitude or longitude.
This function takes two arguments, although only one is strictly required.
The two arguments are ``latitude`` and ``longitude``. One of these arguments
needs to be set to a float, and the other can then be either ignored or set to
a minimum or maximum value.

For example, if we set latitude to 0 N and leave longitude blank, it would
produce a cube along the Equator. On the other hand, if we set latitude to 0
and then set longitude to ``[40., 100.]`` this will produce a transect of the
Equator in the Indian Ocean.

See also :func:`esmvalcore.preprocessor.extract_transect`.


``extract_trajectory``
----------------------

This function extract data along a specified trajectory.
The three arguments are: ``latitudes``, ``longitudes`` and number of point
needed for extrapolation ``number_points``.

If two points are provided, the ``number_points`` argument is used to set a
the number of places to extract between the two end points.

If more than two points are provided, then ``extract_trajectory`` will produce
a cube which has extrapolated the data of the cube to those points, and
``number_points`` is not needed.

Note that this function uses the expensive ``interpolate`` method from
``Iris.analysis.trajectory``, but it may be necessary for irregular grids.

See also :func:`esmvalcore.preprocessor.extract_trajectory`.


.. _cycles:

Cycles
======

The ``_cycles.py`` module contains the following preprocessor functions:

* ``amplitude``: Extract the peak-to-peak amplitude of a cycle aggregated over
  specified coordinates.

``amplitude``
-------------

This function extracts the peak-to-peak amplitude (maximum value minus minimum
value) of a field aggregated over specified coordinates. Its only argument is
``coords``, which can either be a single coordinate (given as :obj:`str`) or
multiple coordinates (given as :obj:`list` of :obj:`str`). Usually, these
coordinates refer to temporal categorised coordinates
:obj:`iris.coord_categorisation`
like `year`, `month`, `day of year`, etc. For example, to extract the amplitude
of the annual cycle for every single year in the data, use ``coords: year``; to
extract the amplitude of the diurnal cycle for every single day in the data,
use ``coords: [year, day_of_year]``.

See also :func:`esmvalcore.preprocessor.amplitude`.


.. _trend:

Trend
=====

The trend module contains the following preprocessor functions:

* ``linear_trend``: Calculate linear trend along a specified coordinate.
* ``linear_trend_stderr``: Calculate standard error of linear trend along a
  specified coordinate.

``linear_trend``
----------------

This function calculates the linear trend of a dataset (defined as slope of an
ordinary linear regression) along a specified coordinate. The only argument of
this preprocessor is ``coordinate`` (given as :obj:`str`; default value is
``'time'``).

See also :func:`esmvalcore.preprocessor.linear_trend`.

``linear_trend_stderr``
-----------------------

This function calculates the standard error of the linear trend of a dataset
(defined as the standard error of the slope in an ordinary linear regression)
along a specified coordinate. The only argument of this preprocessor is
``coordinate`` (given as :obj:`str`; default value is ``'time'``). Note that
the standard error is **not** identical to a confidence interval.

See also :func:`esmvalcore.preprocessor.linear_trend_stderr`.


.. _detrend:

Detrend
=======

ESMValTool also supports detrending along any dimension using
the preprocessor function 'detrend'.
This function has two parameters:

* ``dimension``: dimension to apply detrend on. Default: "time"
* ``method``: It can be ``linear`` or ``constant``. Default: ``linear``

If method is ``linear``, detrend will calculate the linear trend along the
selected axis and subtract it to the data. For example, this can be used to
remove the linear trend caused by climate change on some variables is selected
dimension is time.

If method is ``constant``, detrend will compute the mean along that dimension
and subtract it from the data

See also :func:`esmvalcore.preprocessor.detrend`.


.. _unit conversion:

Unit conversion
===============

Converting units is also supported. This is particularly useful in
cases where different datasets might have different units, for example
when comparing CMIP5 and CMIP6 variables where the units have changed
or in case of observational datasets that are delivered in different
units.

In these cases, having a unit conversion at the end of the processing
will guarantee homogeneous input for the diagnostics.

.. note::
   Conversion is only supported between compatible units! In other
   words, converting temperature units from ``degC`` to ``Kelvin`` works
   fine, changing precipitation units from a rate based unit to an
   amount based unit is not supported at the moment.

See also :func:`esmvalcore.preprocessor.convert_units`.


.. _bias:

Bias
====

The bias module contains the following preprocessor functions:

* ``bias``: Calculate absolute or relative biases with respect to a reference
  dataset

``bias``
--------

This function calculates biases with respect to a given reference dataset. For
this, exactly one input dataset needs to be declared as ``reference_for_bias:
true`` in the recipe, e.g.,

.. code-block:: yaml

  datasets:
    - {dataset: CanESM5, project: CMIP6, ensemble: r1i1p1f1, grid: gn}
    - {dataset: CESM2,   project: CMIP6, ensemble: r1i1p1f1, grid: gn}
    - {dataset: MIROC6,  project: CMIP6, ensemble: r1i1p1f1, grid: gn}
    - {dataset: ERA-Interim, project: OBS6, tier: 3, type: reanaly, version: 1,
       reference_for_bias: true}

In the example above, ERA-Interim is used as reference dataset for the bias
calculation. For this preprocessor, all input datasets need to have identical
dimensional coordinates. This can for example be ensured with the preprocessors
:func:`esmvalcore.preprocessor.regrid` and/or
:func:`esmvalcore.preprocessor.regrid_time`.

The ``bias`` preprocessor supports 4 optional arguments:

   * ``bias_type`` (:obj:`str`, default: ``'absolute'``): Bias type that is
     calculated. Can be ``'absolute'`` (i.e., calculate bias for dataset
     :math:`X` and reference :math:`R` as :math:`X - R`) or ``relative`` (i.e,
     calculate bias as :math:`\frac{X - R}{R}`).
   * ``denominator_mask_threshold`` (:obj:`float`, default: ``1e-3``):
     Threshold to mask values close to zero in the denominator (i.e., the
     reference dataset) during the calculation of relative biases. All values
     in the reference dataset with absolute value less than the given threshold
     are masked out. This setting is ignored when ``bias_type`` is set to
     ``'absolute'``. Please note that for some variables with very small
     absolute values (e.g., carbon cycle fluxes, which are usually :math:`<
     10^{-6}` kg m :math:`^{-2}` s :math:`^{-1}`) it is absolutely essential to
     change the default value in order to get reasonable results.
   * ``keep_reference_dataset`` (:obj:`bool`, default: ``False``): If
     ``True``, keep the reference dataset in the output. If ``False``, drop the
     reference dataset.
   * ``exclude`` (:obj:`list` of :obj:`str`): Exclude specific datasets from
     this preprocessor. Note that this option is only available in the recipe,
     not when using :func:`esmvalcore.preprocessor.bias` directly (e.g., in
     another python script). If the reference dataset has been excluded, an
     error is raised.

Example:

.. code-block:: yaml

    preprocessors:
      preproc_bias:
        bias:
          bias_type: relative
          denominator_mask_threshold: 1e-8
          keep_reference_dataset: true
          exclude: [CanESM2]

See also :func:`esmvalcore.preprocessor.bias`.


.. _Memory use:

Information on maximum memory required
======================================
In the most general case, we can set upper limits on the maximum memory the
analysis will require:


``Ms = (R + N) x F_eff - F_eff`` - when no multi-model analysis is performed;

``Mm = (2R + N) x F_eff - 2F_eff`` - when multi-model analysis is performed;

where

* ``Ms``: maximum memory for non-multimodel module
* ``Mm``: maximum memory for multi-model module
* ``R``: computational efficiency of module; `R` is typically 2-3
* ``N``: number of datasets
* ``F_eff``: average size of data per dataset where ``F_eff = e x f x F``
  where ``e`` is the factor that describes how lazy the data is (``e = 1`` for
  fully realized data) and ``f`` describes how much the data was shrunk by the
  immediately previous module, e.g. time extraction, area selection or level
  extraction; note that for fix_data ``f`` relates only to the time extraction,
  if data is exact in time (no time selection) ``f = 1`` for fix_data so for
  cases when we deal with a lot of datasets ``R + N \approx N``, data is fully
  realized, assuming an average size of 1.5GB for 10 years of `3D` netCDF data,
  ``N`` datasets will require:


``Ms = 1.5 x (N - 1)`` GB

``Mm = 1.5 x (N - 2)`` GB

As a rule of thumb, the maximum required memory at a certain time for
multi-model analysis could be estimated by multiplying the number of datasets by
the average file size of all the datasets; this memory intake is high but also
assumes that all data is fully realized in memory; this aspect will gradually
change and the amount of realized data will decrease with the increase of
``dask`` use.

.. _Other:

Other
=====

Miscellaneous functions that do not belong to any of the other categories.

Clip
----

This function clips data values to a certain minimum, maximum or range. The function takes two
arguments:

* ``minimum``: Lower bound of range. Default: ``None``
* ``maximum``: Upper bound of range. Default: ``None``

The example below shows how to set all values below zero to zero.


.. code-block:: yaml

    preprocessors:
      clip:
        minimum: 0
        maximum: null
.. _recipe_overview:

Overview
********

After ``config-user.yml``, the ``recipe.yml`` is the second file the user needs
to pass to ``esmvaltool`` as command line option, at each run time point.
Recipes contain the data and data analysis information and instructions needed
to run the diagnostic(s), as well as specific diagnostic-related instructions.

Broadly, recipes contain a general section summarizing the provenance and
functionality of the diagnostics, the datasets which need to be run, the
preprocessors that need to be applied, and the diagnostics which need to be run
over the preprocessed data. This information is provided to ESMValTool in four
main recipe sections: :ref:`Documentation <recipe_documentation>`, Datasets_,
Preprocessors_, and Diagnostics_, respectively.

.. _recipe_documentation:

Recipe section: ``documentation``
=================================

The documentation section includes:

- The recipe's author's user name (``authors``, matching the definitions in the
  :ref:`config-ref`)
- The recipe's maintainer's user name (``maintainer``, matching the definitions in the
  :ref:`config-ref`)
- The title of the recipe (``title``)
- A description of the recipe (``description``, written in MarkDown format)
- A list of scientific references (``references``, matching the definitions in
  the :ref:`config-ref`)
- the project or projects associated with the recipe (``projects``, matching
  the definitions in the :ref:`config-ref`)

For example, the documentation section of ``recipes/recipe_ocean_amoc.yml`` is
the following:

.. code-block:: yaml

    documentation:
      title: Atlantic Meridional Overturning Circulation (AMOC) and the drake passage current
      description: |
        Recipe to produce time series figures of the derived variable, the
        Atlantic meridional overturning circulation (AMOC).
        This recipe also produces transect figures of the stream functions for
        the years 2001-2004.

      authors:
        - demo_le

      maintainer:
        - demo_le

      references:
        - demora2018gmd

      projects:
        - ukesm

.. note::

   Note that all authors, projects, and references mentioned in the description
   section of the recipe need to be included in the (locally installed copy of the) file
   `esmvaltool/config-references.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/config-references.yml>`_,
   see :ref:`config-ref`.
   The author name uses the format: ``surname_name``. For instance, John
   Doe would be: ``doe_john``. This information can be omitted by new users
   whose name is not yet included in ``config-references.yml``.

.. _Datasets:

Recipe section: ``datasets``
============================

The ``datasets`` section includes dictionaries that, via key-value pairs, define standardized
data specifications:

- dataset name (key ``dataset``, value e.g. ``MPI-ESM-LR`` or ``UKESM1-0-LL``)
- project (key ``project``, value ``CMIP5`` or ``CMIP6`` for CMIP data,
  ``OBS`` for observational data, ``ana4mips`` for ana4mips data,
  ``obs4MIPs`` for obs4MIPs data, ``EMAC`` for EMAC data)
- experiment (key ``exp``, value e.g. ``historical``, ``amip``, ``piControl``,
  ``RCP8.5``)
- mip (for CMIP data, key ``mip``, value e.g. ``Amon``, ``Omon``, ``LImon``)
- ensemble member (key ``ensemble``, value e.g. ``r1i1p1``, ``r1i1p1f1``)
- sub-experiment id (key `sub_experiment`, value e.g. `s2000`, `s(2000:2002)`,
  for DCPP data only)
- time range (e.g. key-value ``start_year: 1982``, ``end_year: 1990``.
  Please note that `yaml`_ interprets numbers with a leading ``0`` as octal numbers,
  so we recommend to avoid them. For example, use ``128`` to specify the year
  128 instead of ``0128``.)
  Alternatively, the time range can be specified in `ISO 8601 format <https://en.wikipedia.org/wiki/ISO_8601>`_, for both dates
  and periods. Or as a wildcard to work with all available data, set the starting
  point at the first available year, and set the ending point at the last available
  year. The starting point and end point must be separated with '/'.
  (e.g key-value ``timerange: '1982/1990'``)
- model grid (native grid ``grid: gn`` or regridded grid ``grid: gr``, for
  CMIP6 data only).

For example, a datasets section could be:

.. code-block:: yaml

    datasets:
      - {dataset: CanESM2, project: CMIP5, exp: historical, ensemble: r1i1p1, start_year: 2001, end_year: 2004}
      - {dataset: UKESM1-0-LL, project: CMIP6, exp: historical, ensemble: r1i1p1f2, start_year: 2001, end_year: 2004, grid: gn}
      - {dataset: ACCESS-CM2, project: CMIP6, exp: historical, ensemble: r1i1p1f2, timerange: 'P5Y/*', grid: gn}
      - {dataset: EC-EARTH3, alias: custom_alias, project: CMIP6, exp: historical, ensemble: r1i1p1f1, start_year: 2001, end_year: 2004, grid: gn}
      - {dataset: CMCC-CM2-SR5, project: CMIP6, exp: historical, ensemble: r1i1p1f1, timerange: '2001/P10Y', grid: gn}
      - {dataset: HadGEM3-GC31-MM, project: CMIP6, exp: dcppA-hindcast, ensemble: r1i1p1f1, sub_experiment: s2000, grid: gn, start_year: 2000, end_year, 2002}
      - {dataset: BCC-CSM2-MR, project: CMIP6, exp: dcppA-hindcast, ensemble: r1i1p1f1, sub_experiment: s2000, grid: gn, timerange: '*'}

It is possible to define the experiment as a list to concatenate two experiments.
Here it is an example concatenating the `historical` experiment with `rcp85`

.. code-block:: yaml

    datasets:
      - {dataset: CanESM2, project: CMIP5, exp: [historical, rcp85], ensemble: r1i1p1, start_year: 2001, end_year: 2004}

It is also possible to define the ensemble as a list when the two experiments have different ensemble names.
In this case, the specified datasets are concatenated into a single cube:

.. code-block:: yaml

    datasets:
      - {dataset: CanESM2, project: CMIP5, exp: [historical, rcp85], ensemble: [r1i1p1, r1i2p1], start_year: 2001, end_year: 2004}

ESMValTool also supports a simplified syntax to add multiple ensemble members from the same dataset.
In the ensemble key, any element in the form `(x:y)` will be replaced with all numbers from x to y (both inclusive),
adding a dataset entry for each replacement. For example, to add ensemble members r1i1p1 to r10i1p1
you can use the following abbreviated syntax:

.. code-block:: yaml

    datasets:
      - {dataset: CanESM2, project: CMIP5, exp: historical, ensemble: "r(1:10)i1p1", start_year: 2001, end_year: 2004}

It can be included multiple times in one definition. For example, to generate the datasets definitions
for the ensemble members r1i1p1 to r5i1p1 and from r1i2p1 to r5i1p1 you can use:

.. code-block:: yaml

    datasets:
      - {dataset: CanESM2, project: CMIP5, exp: historical, ensemble: "r(1:5)i(1:2)p1", start_year: 2001, end_year: 2004}

Please, bear in mind that this syntax can only be used in the ensemble tag.
Also, note that the combination of multiple experiments and ensembles, like
exp: [historical, rcp85], ensemble: [r1i1p1, "r(2:3)i1p1"] is not supported and will raise an error.

The same simplified syntax can be used to add multiple sub-experiment ids:

.. code-block:: yaml

    datasets:
      - {dataset: MIROC6, project: CMIP6, exp: dcppA-hindcast, ensemble: r1i1p1f1, sub_experiment: s(2000:2002), grid: gn, start_year: 2003, end_year: 2004}

When using the ``timerange`` tag to specify the start and end points, possible values can be as follows:

  - A start and end point specified with a resolution up to seconds (YYYYMMDDThhmmss)
    * ``timerange: '1980/1982'``. Spans from 01/01/1980 to 31/12/1980.
    * ``timerange: '198002/198205'``. Spans from 01/02/1980 to 31/05/1982.
    * ``timerange: '19800302/19820403'``. Spans from 02/03/1980 to 03/04/1982.
    * ``timerange: '19800504T100000/19800504T110000'``. Spans from 04/05/1980 at 10h to 11h.

  - A start point or end point, and a relative period with a resolution up to second (P[n]Y[n]M[n]DT[n]H[n]M[n]S).
    * ``timerange: '1980/P5Y'``. Starting from 01/01/1980, spans 5 years.
    * ``timerange: 'P2Y5M/198202``. Ending at 28/02/1982, spans 2 years and 5 months.
  - A wildcard to load all available years, the first available start point or the last available end point.
    * ``timerange: '*'``. Finds all available years.
    * ``timerange: '*/1982``. Finds first available point, spans to 31/12/1982.
    * ``timerange: '*/P6Y``. Finds first available point, spans 6 years from it.
    * ``timerange: '198003/*``. Starting from 01/03/1980, spans until the last available point.
    * ``timerange: 'P5M/*``. Finds last available point, spans 5 months backwards from it.


Note that this section is not required, as datasets can also be provided in the
Diagnostics_ section.

.. _`yaml`: https://yaml.org/refcard.html

.. _Preprocessors:

Recipe section: ``preprocessors``
=================================

The preprocessor section of the recipe includes one or more preprocessors, each
of which may call the execution of one or several preprocessor functions.

Each preprocessor section includes:

- A preprocessor name (any name, under ``preprocessors``);
- A list of preprocessor steps to be executed (choose from the API);
- Any or none arguments given to the preprocessor steps;
- The order that the preprocessor steps are applied can also be specified using
  the ``custom_order`` preprocessor function.

The following snippet is an example of a preprocessor named ``prep_map`` that
contains multiple preprocessing steps (:ref:`Horizontal regridding` with two
arguments, :ref:`Time operations` with no arguments (i.e., calculating the
average over the time dimension) and :ref:`Multi-model statistics` with two
arguments):

.. code-block:: yaml

    preprocessors:
      prep_map:
        regrid:
          target_grid: 1x1
          scheme: linear
        climate_statistics:
          operator: mean
        multi_model_statistics:
          span: overlap
          statistics: [mean]

.. note::

   In this case no ``preprocessors`` section is needed the workflow will apply
   a ``default`` preprocessor consisting of only basic operations like: loading
   data, applying CMOR checks and fixes (:ref:`CMOR check and dataset-specific
   fixes`) and saving the data to disk.

Preprocessor operations will be applied using the default order
as listed in :ref:`preprocessor_functions`.
Preprocessor tasks can be set to run in the order they are listed in the recipe
by adding ``custom_order: true`` to the preprocessor definition.

.. _Diagnostics:

Recipe section: ``diagnostics``
===============================

The diagnostics section includes one or more diagnostics. Each diagnostic
section will include:

- the variable(s) to preprocess, including the preprocessor to be applied to each variable;
- the diagnostic script(s) to be run;
- a description of the diagnostic and lists of themes and realms that it applies to;
- an optional ``additional_datasets`` section.
- an optional ``title`` and ``description``, used to generate the title and description
  of the ``index.html`` output file.

.. _tasks:

The diagnostics section defines tasks
-------------------------------------
The diagnostic section(s) define the tasks that will be executed when running the recipe.
For each variable a preprocessing task will be defined and for each diagnostic script a
diagnostic task will be defined. If variables need to be derived
from other variables, a preprocessing task for each of the variables
needed to derive that variable will be defined as well. These tasks can be viewed
in the main_log_debug.txt file that is produced every run. Each task has a unique
name that defines the subdirectory where the results of that task are stored. Task
names start with the name of the diagnostic section followed by a '/' and then
the name of the variable section for a preprocessing task or the name of the diagnostic
script section for a diagnostic task.

A (simplified) example diagnostics section could look like

.. code-block:: yaml

  diagnostics:
    diagnostic_name:
      title: Air temperature tutorial diagnostic
      description: A longer description can be added here.
      themes:
        - phys
      realms:
        - atmos
      variables:
        variable_name:
          short_name: ta
          preprocessor: preprocessor_name
          mip: Amon
      scripts:
        script_name:
          script: examples/diagnostic.py


Note that the example recipe above contains a single diagnostic section
called ``diagnostic_name`` and will result in two tasks:

- a preprocessing task called ``diagnostic_name/variable_name`` that will preprocess
  air temperature data for each dataset in the Datasets_ section of the recipe (not shown).
- a diagnostic task called ``diagnostic_name/script_name``

The path to the script provided in the ``script`` option should be
either the absolute path to the script, or the path relative to the
``esmvaltool/diag_scripts`` directory.

Depending on the installation configuration, you may get an error of
"file does not exist" when the system tries to run the diagnostic script
using relative paths. If this happens, use an absolute path instead.

Note that the script should either have the extension for a supported language,
i.e. ``.py``, ``.R``, ``.ncl``, or ``.jl`` for Python, R, NCL, and Julia diagnostics
respectively, or be executable if it is written in any other language.

.. _ancestor-tasks:

Ancestor tasks
--------------
Some tasks require the result of other tasks to be ready before they can start,
e.g. a diagnostic script needs the preprocessed variable data to start. Thus
each tasks has zero or more ancestor tasks. By default, each diagnostic task
in a diagnostic section has all variable preprocessing tasks in that same section
as ancestors. However, this can be changed using the ``ancestors`` keyword. Note
that wildcard expansion can be used to define ancestors.

.. code-block:: yaml

  diagnostics:
    diagnostic_1:
      variables:
        airtemp:
          short_name: ta
          preprocessor: preprocessor_name
          mip: Amon
      scripts:
        script_a:
          script: diagnostic_a.py
    diagnostic_2:
      variables:
        precip:
          short_name: pr
          preprocessor: preprocessor_name
          mip: Amon
      scripts:
        script_b:
          script: diagnostic_b.py
          ancestors: [diagnostic_1/script_a, precip]


The example recipe above will result in four tasks:

- a preprocessing task called ``diagnostic_1/airtemp``
- a diagnostic task called ``diagnostic_1/script_a``
- a preprocessing task called ``diagnostic_2/precip``
- a diagnostic task called ``diagnostic_2/script_b``

the preprocessing tasks do not have any ancestors, while the diagnostic_a.py
script will receive the preprocessed air temperature data
(has ancestor ``diagnostic_1/airtemp``) and the diagnostic_b.py
script will receive the results of diagnostic_a.py and the preprocessed precipitation
data (has ancestors ``diagnostic_1/script_a`` and ``diagnostic_2/precip``).

Task priority
-------------
Tasks are assigned a priority, with tasks appearing earlier on in the recipe
getting higher priority. The tasks will be executed sequentially or in parellel,
depending on the setting of ``max_parallel_tasks`` in the :ref:`user configuration file`.
When there are fewer than ``max_parallel_tasks`` running, tasks will be started
according to their priority. For obvious reasons, only tasks that are not waiting for
ancestor tasks can be started. This feature makes it possible to
reduce the processing time of recipes with many tasks, by placing tasks that
take relatively long near the top of the recipe. Of course this only works when
settings ``max_parallel_tasks`` to a value larger than 1. The current priority
and run time of individual tasks can be seen in the log messages shown when
running the tool (a lower number means higher priority).

Variable and dataset definitions
--------------------------------
To define a variable/dataset combination that corresponds to an actual
variable from a dataset, the keys in each variable section
are combined with the keys of each dataset definition. If two versions of the same
key are provided, then the key in the datasets section will take precedence
over the keys in variables section. For many recipes it makes more sense to
define the ``start_year`` and ``end_year`` items in the variable section,
because the diagnostic script assumes that all the data has the same time
range.

Variable short names usually do not change between datasets supported by
ESMValCore, as they are usually changed to match CMIP. Nevertheless, there are
small changes in variable names in CMIP6 with respect to CMIP5 (i.e. sea ice
concentration changed from ``sic`` to ``siconc``). ESMValCore is aware of some
of them and can do the automatic translation when needed. It will even do the
translation in the preprocessed file so the diagnostic does not have to deal
with this complexity, setting the short name in all files to match the one used
by the recipe. For example, if ``sic`` is requested, ESMValCore will
find ``sic`` or ``siconc`` depending on the project, but all preprocessed files
while use ``sic`` as their short_name. If the recipe requested ``siconc``, the
preprocessed files will be identical except that they will use the short_name
``siconc`` instead.

Diagnostic and variable specific datasets
-----------------------------------------
The ``additional_datasets`` option can be used to add datasets beyond those
listed in the Datasets_ section. This is useful if specific datasets need to
be used only by a specific diagnostic or variable, i.e. it can be added both
at diagnostic level, where it will apply to all variables in that diagnostic
section or at individual variable level. For example, this can be a good way
to add observational datasets, which are usually variable-specific.

Running a simple diagnostic
---------------------------
The following example, taken from ``recipe_ocean_example.yml``, shows a
diagnostic named `diag_map`, which loads the temperature at the ocean surface
between the years 2001 and 2003 and then passes it to the ``prep_map``
preprocessor. The result of this process is then passed to the ocean diagnostic
map script, ``ocean/diagnostic_maps.py``.

.. code-block:: yaml

  diagnostics:

    diag_map:
      title: Global Ocean Surface regridded temperature map
      description: Add a longer description here.
      variables:
        tos: # Temperature at the ocean surface
          preprocessor: prep_map
          start_year: 2001
          end_year: 2003
      scripts:
        Global_Ocean_Surface_regrid_map:
          script: ocean/diagnostic_maps.py

Passing arguments to a diagnostic script
----------------------------------------
The diagnostic script section(s) may include custom arguments that can be used by
the diagnostic script; these arguments are stored at runtime in a dictionary
that is then made available to the diagnostic script via the interface link,
independent of the language the diagnostic script is written in. Here is an
example of such groups of arguments:

.. code-block:: yaml

    scripts:
      autoassess_strato_test_1: &autoassess_strato_test_1_settings
        script: autoassess/autoassess_area_base.py
        title: "Autoassess Stratosphere Diagnostic Metric MPI-MPI"
        area: stratosphere
        control_model: MPI-ESM-LR
        exp_model: MPI-ESM-MR
        obs_models: [ERA-Interim]  # list to hold models that are NOT for metrics but for obs operations
        additional_metrics: [ERA-Interim, inmcm4]  # list to hold additional datasets for metrics

In this example, apart from specifying the diagnostic script ``script:
autoassess/autoassess_area_base.py``, we pass a suite of parameters to be used
by the script (``area``, ``control_model`` etc). These parameters are stored in
key-value pairs in the diagnostic configuration file, an interface file that
can be used by importing the ``run_diagnostic`` utility:

.. code-block:: python

   from esmvaltool.diag_scripts.shared import run_diagnostic

   # write the diagnostic code here e.g.
   def run_some_diagnostic(my_area, my_control_model, my_exp_model):
       """Diagnostic to be run."""
       if my_area == 'stratosphere':
           diag = my_control_model / my_exp_model
           return diag

   def main(cfg):
       """Main diagnostic run function."""
       my_area = cfg['area']
       my_control_model = cfg['control_model']
       my_exp_model = cfg['exp_model']
       run_some_diagnostic(my_area, my_control_model, my_exp_model)

   if __name__ == '__main__':

       with run_diagnostic() as config:
           main(config)

This way a lot of the optional arguments necessary to a diagnostic are at the
user's control via the recipe.

Running your own diagnostic
---------------------------
If the user wants to test a newly-developed ``my_first_diagnostic.py`` which
is not yet part of the ESMValTool diagnostics library, he/she do it by passing
the absolute path to the diagnostic:

.. code-block:: yaml

  diagnostics:

    myFirstDiag:
      title: Let's do some science!
      description: John Doe wrote a funny diagnostic
      variables:
        tos: # Temperature at the ocean surface
          preprocessor: prep_map
          start_year: 2001
          end_year: 2003
      scripts:
        JoeDiagFunny:
          script: /home/users/john_doe/esmvaltool_testing/my_first_diagnostic.py

This way the user may test a new diagnostic thoroughly before committing to the
GitHub repository and including it in the ESMValTool diagnostics library.

Re-using parameters from one ``script`` to another
--------------------------------------------------
Due to ``yaml`` features it is possible to recycle entire diagnostics sections
for use with other diagnostics. Here is an example:

.. code-block:: yaml

    scripts:
      cycle: &cycle_settings
        script: perfmetrics/main.ncl
        plot_type: cycle
        time_avg: monthlyclim
      grading: &grading_settings
        <<: *cycle_settings
        plot_type: cycle_latlon
        calc_grading: true
        normalization: [centered_median, none]

In this example the hook ``&cycle_settings`` can be used to pass the ``cycle:``
parameters to ``grading:`` via the shortcut ``<<: *cycle_settings``.
.. _recipe:

The recipe format
*****************

.. toctree::
   :maxdepth: 1

    Overview <overview>
    Preprocessor <preprocessor>
 .. _preprocessor_function:

Preprocessor function
*********************

Preprocessor functions are located in :py:mod:`esmvalcore.preprocessor`.
To add a new preprocessor function, start by finding a likely looking file to
add your function to in
`esmvalcore/preprocessor <https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/preprocessor>`_.
Create a new file in that directory if you cannot find a suitable place.

The function should look like this:


.. code-block:: python

    def example_preprocessor_function(
        cube,
        example_argument,
        example_optional_argument=5,
    ):
        """Compute an example quantity.

        A more extensive explanation of the computation can be added here. Add
        references to scientific literature if available.

        Parameters
        ----------
        cube: iris.cube.Cube
           Input cube.

        example_argument: str
           Example argument, the value of this argument can be provided in the
           recipe. Describe what valid values are here. In this case, a valid
           argument is the name of a dimension of the input cube.

        example_optional_argument: int, optional
           Another example argument, the value of this argument can optionally
           be provided in the recipe. Describe what valid values are here.

        Returns
        -------
        iris.cube.Cube
          The result of the example computation.
        """

        # Replace this with your own computation
        cube = cube.collapsed(example_argument, iris.analysis.MEAN)

        return cube


The above function needs to be imported in the file
`esmvalcore/preprocessor/__init__.py <https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/preprocessor/__init__.py>`__:

.. code-block:: python

    from ._example_module import example_preprocessor_function

    __all__ = [
    ...
    'example_preprocessor_function',
    ...
    ]

The location in the ``__all__`` list above determines the default order in which
preprocessor functions are applied, so carefully consider where you put it
and ask for advice if needed.

The preprocessor function above can then be used from the :ref:`preprocessors`
like this:

.. code-block:: yaml

   preprocessors:
     example_preprocessor:
       example_preprocessor_function:
         example_argument: median
         example_optional_argument: 6

The optional argument (in this example: ``example_optional_argument``) can be
omitted in the recipe.

Lazy and real data
==================

Preprocessor functions should support both
:ref:`real and lazy data <iris:real_and_lazy_data>`.
This is vital for supporting the large datasets that are typically used with
the ESMValCore.
If the data of the incoming cube has been realized (i.e. ``cube.has_lazy_data()``
returns ``False`` so ``cube.core_data()`` is a `NumPy <https://numpy.org/>`__
array), the returned cube should also have realized data.
Conversely, if the incoming cube has lazy data (i.e. ``cube.has_lazy_data()``
returns ``True`` so ``cube.core_data()`` is a
`Dask array <https://docs.dask.org/en/latest/array.html>`__), the returned
cube should also have lazy data.
Note that NumPy functions will often call their Dask equivalent if it exists
and if their input array is a Dask array, and vice versa.

Note that preprocessor functions should preferably be small and just call the
relevant :ref:`iris <iris_docs>` code.
Code that is more involved, e.g. lots of work with Numpy and Dask arrays,
and more broadly applicable, should be implemented in iris instead.

Documentation
=============

The documentation in the function docstring will be shown in
the :ref:`preprocessor_functions` chapter.
In addition, you should add documentation on how to use the new preprocessor
function from the recipe in
`doc/recipe/preprocessor.rst <https://github.com/ESMValGroup/ESMValCore/tree/main/doc/recipe/preprocessor.rst>`__
so it is shown in the :ref:`preprocessor` chapter.
See the introduction to :ref:`documentation` for more information on how to
best write documentation.

Tests
=====

Tests are should be implemented for new or modified preprocessor functions.
For an introduction to the topic, see :ref:`tests`.

Unit tests
----------

To add a unit test for the preprocessor function from the example above, create
a file called
``tests/unit/preprocessor/_example_module/test_example_preprocessor_function.py``
and add the following content:

.. code-block:: python

    """Test function `esmvalcore.preprocessor.example_preprocessor_function`."""
    import cf_units
    import dask.array as da
    import iris
    import numpy as np
    import pytest

    from esmvalcore.preprocessor import example_preprocessor_function


    @pytest.mark.parametrize('lazy', [True, False])
    def test_example_preprocessor_function(lazy):
        """Test that the computed result is as expected."""

        # Construct the input cube
        data = np.array([1, 2], dtype=np.float32)
        if lazy:
            data = da.asarray(data, chunks=(1, ))
        cube = iris.cube.Cube(
            data,
            var_name='tas',
            units='K',
        )
        cube.add_dim_coord(
            iris.coords.DimCoord(
                np.array([0.5, 1.5], dtype=np.float64),
                bounds=np.array([[0, 1], [1, 2]], dtype=np.float64),
                standard_name='time',
                units=cf_units.Unit('days since 1950-01-01 00:00:00',
                                    calendar='gregorian'),
            ),
            0,
        )

        # Compute the result
        result = example_preprocessor_function(cube, example_argument='time')

        # Check that lazy data is returned if and only if the input is lazy
        assert result.has_lazy_data() is lazy

        # Construct the expected result cube
        expected = iris.cube.Cube(
            np.array(1.5, dtype=np.float32),
            var_name='tas',
            units='K',
        )
        expected.add_aux_coord(
            iris.coords.AuxCoord(
                np.array([1], dtype=np.float64),
                bounds=np.array([[0, 2]], dtype=np.float64),
                standard_name='time',
                units=cf_units.Unit('days since 1950-01-01 00:00:00',
                                    calendar='gregorian'),
            ))
        expected.add_cell_method(
            iris.coords.CellMethod(method='mean', coords=('time', )))

        # Compare the result of the computation with the expected result
        print('result:', result)
        print('expected result:', expected)
        assert result == expected


In this test we used the decorator
`pytest.mark.parametrize <https://docs.pytest.org/en/stable/parametrize.html>`_
to test two scenarios, with both lazy and realized data, with a single test.


Sample data tests
-----------------

The idea of adding :ref:`sample data tests <sample_data_tests>` is to check that
preprocessor functions work with realistic data.
This also provides an easy way to add regression tests, though these should
preferably be implemented as unit tests instead, because using the sample data
for this purpose is slow.
To add a test using the sample data, create a file
``tests/sample_data/preprocessor/example_preprocessor_function/test_example_preprocessor_function.py``
and add the following content:

.. code-block:: python

    """Test function `esmvalcore.preprocessor.example_preprocessor_function`."""
    from pathlib import Path

    import esmvaltool_sample_data
    import iris
    import pytest

    from esmvalcore.preprocessor import example_preprocessor_function


    @pytest.mark.use_sample_data
    def test_example_preprocessor_function():
        """Regression test to check that the computed result is as expected."""
        # Load an example input cube
        cube = esmvaltool_sample_data.load_timeseries_cubes(mip_table='Amon')[0]

        # Compute the result
        result = example_preprocessor_function(cube, example_argument='time')

        filename = Path(__file__).with_name('example_preprocessor_function.nc')
        if not filename.exists():
            # Create the file the expected result if it doesn't exist
            iris.save(result, target=str(filename))
            raise FileNotFoundError(
                f'Reference data was missing, wrote new copy to {filename}')

        # Load the expected result cube
        expected = iris.load_cube(str(filename))

        # Compare the result of the computation with the expected result
        print('result:', result)
        print('expected result:', expected)
        assert result == expected


This will use a file from the sample data repository as input.
The first time you run the test, the computed result will be stored in the file
``tests/sample_data/preprocessor/example_preprocessor_function/example_preprocessor_function.nc``
Any subsequent runs will re-load the data from file and check that it did not
change.
Make sure the stored results are small, i.e. smaller than 100 kilobytes, to
keep the size of the ESMValCore repository small.

Using multiple datasets as input
================================

The name of the first argument of the preprocessor function should in almost all
cases be ``cube``.
Only when implementing a preprocessor function that uses all datasets as input,
the name of the first argument should be ``products``.
If you would like to implement this type of preprocessor function, start by
having a look at the existing functions, e.g.
:py:func:`esmvalcore.preprocessor.multi_model_statistics` or
:py:func:`esmvalcore.preprocessor.mask_fillvalues`.
.. _fixing_data:

***********
Fixing data
***********

The baseline case for ESMValCore input data is CMOR fully compliant
data that is read using Iris' :func:`iris:iris.load_raw`.
ESMValCore also allows for some departures from compliance (see
:ref:`cmor_check_strictness`). Beyond that situation, some datasets
(either model or observations) contain (known) errors that would
normally prevent them from being processed. The issues can be in
the metadata describing the dataset and/or in the actual data.
Typical examples of such errors are missing or wrong attributes (e.g.
attribute ''units'' says 1e-9 but data are actually in 1e-6), missing or
mislabeled coordinates (e.g. ''lev'' instead of ''plev'' or missing
coordinate bounds like ''lat_bnds'') or problems with the actual data
(e.g. cloud liquid water only instead of sum of liquid + ice as specified by the CMIP data request).

As an extreme case, some data sources simply are not NetCDF
files and must go through some other data load function.

The ESMValCore can apply on the fly fixes to such datasets when
issues can be fixed automatically. This is implemented for a set
of `Natively supported non-CMIP datasets`_. The following provides
details on how to design such fixes.

.. note::

  **CMORizer scripts**. Support for many observational and reanalysis
  datasets is also possible through a priori reformatting by
  :ref:`CMORizer scripts in the ESMValTool <esmvaltool:new-dataset>`,
  which are rather relevant for datasets of small volume

.. _fix_structure:

Fix structure
=============

Fixes are Python classes stored in
``esmvalcore/cmor/_fixes/[PROJECT]/[DATASET].py`` that derive from
:class:`esmvalcore.cmor._fixes.fix.Fix` and are named after the short name of
the variable they fix. You can also use the names of ``mip`` tables (e.g.,
``Amon``, ``Lmon``, ``Omon``, etc.) if you want the fix to be applied to all
variables of that table in the dataset or ``AllVars`` if you want the fix to be
applied to the whole dataset.

.. warning::
    Be careful to replace any ``-`` with ``_`` in your dataset name.
    We need this replacement to have proper python module names.

The fixes are automatically loaded and applied when the dataset is preprocessed.
They are a special type of :ref:`preprocessor function <preprocessor_function>`,
called by the preprocessor functions
:py:func:`esmvalcore.preprocessor.fix_file`,
:py:func:`esmvalcore.preprocessor.fix_metadata`, and
:py:func:`esmvalcore.preprocessor.fix_data`.

Fixing a dataset
================

To illustrate the process of creating a fix we are going to construct a new
one from scratch for a fictional dataset. We need to fix a CMIPX model
called PERFECT-MODEL that is reporting a missing latitude coordinate for
variable tas.

Check the output
----------------

Next to the error message, you should see some info about the iris cube: size,
coordinates. In our example it looks like this:

.. code-block:: python

    air_temperature/ (K) (time: 312; altitude: 90; longitude: 180)
        Dimension coordinates:
            time                                     x              -              -
            altitude                                 -              x              -
            longitude                                -              -              x
        Auxiliary coordinates:
            day_of_month                             x              -              -
            day_of_year                              x              -              -
            month_number                             x              -              -
            year                                     x              -              -
        Attributes:
            {'cmor_table': 'CMIPX', 'mip': 'Amon', 'short_name': 'tas', 'frequency': 'mon'})


So now the mistake is clear: the latitude coordinate is badly named and the
fix should just rename it.

Create the fix
--------------

We start by creating the module file. In our example the path will be
``esmvalcore/cmor/_fixes/CMIPX/PERFECT_MODEL.py``. If it already exists
just add the class to the file, there is no limit in the number of fixes
we can have in any given file.

Then we have to create the class for the fix deriving from
:class:`esmvalcore.cmor._fixes.Fix`

.. code-block:: python

    """Fixes for PERFECT-MODEL."""
    from esmvalcore.cmor.fix import Fix

    class tas(Fix):
         """Fixes for tas variable.""""

Next we must choose the method to use between the ones offered by the
Fix class:

- ``fix_file`` : should be used only to fix errors that prevent data loading.
  As a rule of thumb, you should only use it if the execution halts before
  reaching the checks.

- ``fix_metadata`` : you want to change something in the cube that is not
  the data (e.g variable or coordinate names, data units).

- ``fix_data``: you need to fix the data. Beware: coordinates data values are
  part of the metadata.

In our case we need to rename the coordinate ``altitude`` to ``latitude``,
so we will implement the ``fix_metadata`` method:

.. code-block:: python

    """Fixes for PERFECT-MODEL."""
    from esmvalcore.cmor.fix import Fix

    class tas(Fix):
        """Fixes for tas variable.""""

        def fix_metadata(self, cubes):
            """
            Fix metadata for tas.

            Fix the name of the latitude coordinate, which is called altitude
            in the original file.
            """"
            # Sometimes Iris will interpret the data as multiple cubes.
            # Good CMOR datasets will only show one but we support the
            # multiple cubes case to be able to fix the errors that are
            # leading to that extra cubes.
            # In our case this means that we can safely assume that the
            # tas cube is the first one
            tas_cube = cubes[0]
            latitude = tas_cube.coord('altitude')

            # Fix the names. Latitude values, units and
            latitude.short_name = 'lat'
            latitude.standard_name = 'latitude'
            latitude.long_name = 'latitude'
            return cubes

This will fix the error. The next time you run ESMValTool you will find that the error
is fixed on the fly and, hopefully, your recipe will run free of errors.
The ``cubes`` argument to the ``fix_metadata`` method will contain all cubes
loaded from a single input file.
Some care may need to be taken that the right cube is selected and fixed in case
multiple cubes are created.
Usually this happens when a coordinate is mistakenly loaded as a cube, because
the input data does not follow the
`CF Conventions <https://cfconventions.org/>`__.

Sometimes other errors can appear after you fix the first one because they were
hidden by it. In our case, the latitude coordinate could have bad units or
values outside the valid range for example. Just extend your fix to address those
errors.

Finishing
---------

Chances are that you are not the only one that wants to use that dataset and
variable. Other users could take advantage of your fixes as
soon as possible. Please, create a separated pull request for the fix and
submit it.

It will also be very helpful if you just scan a couple of other variables from
the same dataset and check if they share this error. In case that you find that
it is a general one, you can change the fix name to the corresponding ``mip``
table name (e.g., ``Amon``, ``Lmon``, ``Omon``, etc.) so it gets executed for
all variables in that table in the dataset or to ``AllVars`` so it gets
executed for all variables in the dataset. If you find that this is shared only
by a handful of similar vars you can just make the fix for those new vars
derive from the one you just created:

.. code-block:: python

    """Fixes for PERFECT-MODEL."""
    from esmvalcore.cmor.fix import Fix

    class tas(Fix):
        """Fixes for tas variable.""""

        def fix_metadata(self, cubes):
            """
            Fix metadata for tas.

            Fix the name of the latitude coordinate, which is called altitude
            in the original file.
            """"
            # Sometimes Iris will interpret the data as multiple cubes.
            # Good CMOR datasets will only show one but we support the
            # multiple cubes case to be able to fix the errors that are
            # leading to that extra cubes.
            # In our case this means that we can safely assume that the
            # tas cube is the first one
            tas_cube = cubes[0]
            latitude = tas_cube.coord('altitude')

            # Fix the names. Latitude values, units and
            latitude.short_name = 'lat'
            latitude.standard_name = 'latitude'
            latitude.long_name = 'latitude'
            return cubes


    class ps(tas):
        """Fixes for ps variable."""


Common errors
=============

The above example covers one of the most common cases: variables / coordinates that
have names that do not match the expected. But there are some others that use
to appear frequently. This section describes the most common cases.

Bad units declared
------------------

It is quite common that a variable declares to be using some units but the data
is stored in another. This can be solved by overwriting the units attribute
with the actual data units.

.. code-block:: python

    def fix_metadata(self, cubes):
        cube.units = 'real_units'


Detecting this error can be tricky if the units are similar enough. It also
has a good chance of going undetected until you notice strange results in
your diagnostic.

For the above example, it can be useful to access the variable definition
and associated coordinate definitions as provided by the CMOR table.
For example:

.. code-block:: python

    def fix_metadata(self, cubes):
        cube.units = self.vardef.units

To learn more about what is available in these definitions, see:
:class:`esmvalcore.cmor.table.VariableInfo` and
:class:`esmvalcore.cmor.table.CoordinateInfo`.



Coordinates missing
-------------------

Another common error is to have missing coordinates. Usually it just means
that the file does not follow the CF-conventions and Iris can therefore not interpret it.

If this is the case, you should see a warning from the ESMValTool about
discarding some cubes in the fix metadata step. Just before that warning you
should see the full list of cubes as read by Iris. If that list contains your
missing coordinate you can create a fix for this model:

.. code-block:: bash

    def fix_metadata(self, cubes):
        coord_cube = cubes.extract_strict('COORDINATE_NAME')
        # Usually this will correspond to an auxiliary coordinate
        # because the most common error is to forget adding it to the
        # coordinates attribute
        coord = iris.coords.AuxCoord(
            coord_cube.data,
            var_name=coord_cube.var_name,
            standard_name=coord_cube.standard_name,
            long_name=coord_cube.long_name,
            units=coord_cube.units,
        }

        # It may also have bounds as another cube
        coord.bounds = cubes.extract_strict('BOUNDS_NAME').data

        data_cube = cubes.extract_strict('VAR_NAME')
        data_cube.add_aux_coord(coord, DIMENSIONS_INDEX_TUPLE)
        return [data_cube]


.. _cmor_check_strictness:

Customizing checker strictness
==============================

The data checker classifies its issues using four different levels of
severity. From highest to lowest:

 - ``CRITICAL``: issues that most of the time will have severe consequences.
 - ``ERROR``: issues that usually lead to unexpected errors, but can be safely
   ignored sometimes.
 - ``WARNING``: something is not up to the standard but is unlikely to have
   consequences later.
 - ``DEBUG``: any info that the checker wants to communicate. Regardless of
   checker strictness, those will always be reported as debug messages.

Users can have control about which levels of issues are interpreted as errors,
and therefore make the checker fail or warnings or debug messages.
For this purpose there is an optional command line option `--check-level`
that can take a number of values, listed below from the lowest level of
strictness to the highest:

- ``ignore``: all issues, regardless of severity, will be reported as
  warnings. Checker will never fail. Use this at your own risk.
- ``relaxed``: only CRITICAL issues are treated as errors. We recommend not to
  rely on this mode, although it can be useful if there are errors preventing
  the run that you are sure you can manage on the diagnostics or that will
  not affect you.
- ``default``: fail if there are any CRITICAL or ERROR issues (DEFAULT); Provides
  a good measure of safety.
- ``strict``: fail if there are any warnings, this is the highest level of
  strictness. Mostly useful for checking datasets that you have produced, to
  be sure that future users will not be distracted by inoffensive warnings.


Natively supported non-CMIP datasets
====================================

Some fixed datasets and native models formats are supported through
the ``native6`` project or through a dedicated project.

Observational Datasets
----------------------
Put the files containing the data in the directory that you have configured
for the ``native6`` project in your :ref:`user configuration file`, in a
subdirectory called ``Tier{tier}/{dataset}/{version}/{frequency}/{short_name}``.
Replace the items in curly braces by the values used in the variable/dataset
definition in the :ref:`recipe <recipe_overview>`.
Below is a list of datasets currently supported.

ERA5
~~~~

- Supported variables: ``clt``, ``evspsbl``, ``evspsblpot``, ``mrro``, ``pr``, ``prsn``, ``ps``, ``psl``, ``ptype``, ``rls``, ``rlds``, ``rsds``, ``rsdt``, ``rss``, ``uas``, ``vas``, ``tas``, ``tasmax``, ``tasmin``, ``tdps``, ``ts``, ``tsn`` (``E1hr``/``Amon``), ``orog`` (``fx``)
- Tier: 3

MSWEP
~~~~~

- Supported variables: ``pr``
- Supported frequencies: ``mon``, ``day``, ``3hr``.
- Tier: 3

For example for monthly data, place the files in the ``/Tier3/MSWEP/latestversion/mon/pr`` subdirectory of your ``native6`` project location.

.. note::
  For monthly data (V220), the data must be postfixed with the date, i.e. rename ``global_monthly_050deg.nc`` to ``global_monthly_050deg_197901-201710.nc``

For more info: http://www.gloh2o.org/

.. _fixing_native_models:

Native models
-------------

The following models are natively supported through the procedure described
above (:ref:`fix_structure`) and at :ref:`configure_native_models`:

ICON
~~~~

The ESMValTool is able to read native `ICON
<https://code.mpimet.mpg.de/projects/iconpublic>`_ model output. Example
dataset entries could look like this:

.. code-block:: yaml

  datasets:
    - {project: ICON, dataset: ICON, component: atm, version: 2.6.1,
       exp: amip, grid: R2B5, ensemble: r1v1i1p1l1f1, mip: Amon,
       short_name: tas, var_type: atm_2d_ml, start_year: 2000, end_year: 2014}
    - {project: ICON, dataset: ICON, component: atm, version: 2.6.1,
       exp: amip, grid: R2B5, ensemble: r1v1i1p1l1f1, mip: Amon,
       short_name: ta, var_type: atm_3d_ml, start_year: 2000, end_year: 2014}

Please note the duplication of the name ``ICON`` in ``project`` and
``dataset``, which is necessary to comply with ESMValTool's data finding and
CMORizing functionalities.

Similar to any other fix, the ICON fix allows the use of :ref:`extra
facets<extra_facets>`. By default, the file :download:`icon-mapping.yml
</../esmvalcore/_config/extra_facets/icon-mapping.yml>` is used for that
purpose. For some variables, extra facets are necessary; otherwise ESMValTool
cannot read them properly. Supported keys for extra facets are:

============= ===============================================================
Key           Description
============= ===============================================================
``latitude``  Standard name of the latitude coordinate in the raw input file
``longitude`` Standard name of the longitude coordinate in the raw input file
``raw_name``  Variable name of the variables in the raw input file
============= ===============================================================

IPSL-CM6
~~~~~~~~

Both output formats (i.e. the ``Output`` and the ``Analyse / Time series``
formats) are supported, and should be configured in recipes as e.g.:

.. code-block:: yaml

  datasets:
    - {simulation: CM61-LR-hist-03.1950, exp: piControl, out: Analyse, freq: TS_MO,
       account: p86caub,  status: PROD, dataset: IPSL-CM6, project: IPSLCM,
       root: /thredds/tgcc/store}
    - {simulation: CM61-LR-hist-03.1950, exp: historical, out: Output, freq: MO,
       account: p86caub,  status: PROD, dataset: IPSL-CM6, project: IPSLCM,
       root: /thredds/tgcc/store}

.. _ipslcm_extra_facets_example:

The ``Output`` format is an example of a case where variables are grouped in
multi-variable files, which name cannot be computed directly from datasets
attributes alone but requires to use an extra_facets file, which principles are
explained in :ref:`extra_facets`, and which content is :download:`available here
</../esmvalcore/_config/extra_facets/ipslcm-mappings.yml>`. These multi-variable
files must also undergo some data selection.

.. _extra-facets-fixes:

Use of extra facets in fixes
============================
Extra facets are a mechanism to provide additional information for certain kinds
of data. The general approach is described in :ref:`extra_facets`. Here, we
describe how they can be used in fixes to mold data into the form required by
the applicable standard. For example, if the input data is part of an
observational product that delivers surface temperature with a variable name of
`t2m` inside a file named `2m_temperature_1950_monthly.nc`, but the same
variable is called `tas` in the applicable standard, a fix can be created that
reads the original variable from the correct file, and provides a renamed
variable to the rest of the processing chain.

Normally, the applicable standard for variables is CMIP6.

For more details, refer to existing uses of this feature as examples,
as e.g. :ref:`for IPSL-CM6<ipslcm_extra_facets_example>`.
Development
***********

To get started developing, have a look at our
:ref:`contribution guidelines <contributing>`.
This chapter describes how to implement the most commonly contributed new
features.

.. toctree::
   :maxdepth: 1

    Preprocessor function <preprocessor_function>
    Fixing data <fixing_data>
    Deriving a variable <derivation>
.. _derivation:

*******************
Deriving a variable
*******************

The variable derivation preprocessor module allows to derive variables which are
not in the CMIP standard data request using standard variables as input.
This is a special type of :ref:`preprocessor function <preprocessor_function>`.
All derivation scripts are located in
`esmvalcore/preprocessor/_derive/ <https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/preprocessor/_derive>`_.
A typical example looks like this:

.. code-block:: py

   """Derivation of variable `dummy`."""
   from ._baseclass import DerivedVariableBase


   class DerivedVariable(DerivedVariableBase):
       """Derivation of variable `dummy`."""

       @staticmethod
       def required(project):
           """Declare the variables needed for derivation."""
           mip = 'fx'
           if project == 'CMIP6':
               mip = 'Ofx'
           required = [
               {'short_name': 'var_a'},
               {'short_name': 'var_b', 'mip': mip, 'optional': True},
           ]
           return required

       @staticmethod
       def calculate(cubes):
           """Compute `dummy`."""

           # `cubes` is a CubeList containing all required variables.
           cube = do_something_with(cubes)

           # Return single cube at the end
           return cube

The static function ``required(project)`` returns a ``list`` of ``dict``
containing all required variables for deriving the derived variable. Its only
argument is the ``project`` of the specific dataset. In this particular
example script, the derived variable ``dummy`` is derived from ``var_a`` and
``var_b``. It is possible to specify arbitrary attributes for each required
variable, e.g. ``var_b`` uses the mip ``fx`` (or ``Ofx`` in the case of
CMIP6) instead of the original one of ``dummy``. Note that you can also declare
a required variable as ``optional=True``, which allows the skipping of this
particular variable during data extraction. For example, this is useful for
fx variables which are often not available for observational datasets.
Otherwise, the tool will fail if not all required variables are available for
all datasets.

The actual derivation takes place in the static function ``calculate(cubes)``
which returns a single ``cube`` containing the derived variable. Its only
argument ``cubes`` is a ``CubeList`` containing all required variables.
.. _api_recipe:

Recipes
=======

This section describes the :py:mod:`~esmvalcore.experimental.recipe` submodule of the API (:py:mod:`esmvalcore.experimental`).

Recipe metadata
***************

:py:class:`~esmvalcore.experimental.recipe.Recipe` is a class that holds metadata from a recipe.

.. code-block:: python

    >>> Recipe('path/to/recipe_python.yml')
    recipe = Recipe('Recipe Python')

Printing the recipe will give a nice overview of the recipe:

.. code-block:: python

    >>> print(recipe)
    ## Recipe python

    Example recipe that plots a map and timeseries of temperature.

    ### Authors
     - Bouwe Andela (NLeSC, Netherlands; https://orcid.org/0000-0001-9005-8940)
     - Mattia Righi (DLR, Germany; https://orcid.org/0000-0003-3827-5950)

    ### Maintainers
     - Manuel Schlund (DLR, Germany; https://orcid.org/0000-0001-5251-0158)

    ### Projects
     - DLR project ESMVal
     - Copernicus Climate Change Service 34a Lot 2 (MAGIC) project

    ### References
     - Please acknowledge the project(s).

Running a recipe
****************

To run the recipe, call the :py:meth:`~esmvalcore.experimental.recipe.Recipe.run` method.

.. code-block:: python

    >>> output = recipe.run()
    <log messages>

By default, a new :py:class:`~esmvalcore.experimental.config.Session` is automatically created, so that data are never overwritten.
Data are stored in the ``esmvaltool_output`` directory specified in the config.
Sessions can also be explicitly specified.

.. code-block:: python

    >>> from esmvalcore.experimental import CFG
    >>> session = CFG.start_session('my_session')
    >>> output = recipe.run(session)
    <log messages>

:py:meth:`~esmvalcore.experimental.recipe.Recipe.run` returns an dictionary of objects that can be used to inspect
the output of the recipe. The output is an instance of :py:class:`~esmvalcore.experimental.recipe_output.ImageFile` or
:py:class:`~esmvalcore.experimental.recipe_output.DataFile` depending on its type.

For working with recipe output, see: :ref:`api_recipe_output`.

Running a single diagnostic or preprocessor task
************************************************

The python example recipe contains 5 tasks:

Preprocessors:

- ``timeseries/tas_amsterdam``
- ``timeseries/script1``
- ``map/tas``

Diagnostics:

- ``timeseries/tas_global``
- ``map/script1``

To run a single diagnostic or preprocessor, the name of the task can be passed
as an argument to :py:meth:`~esmvalcore.experimental.recipe.Recipe.run`. If a diagnostic
is passed, all ancestors will automatically be run too.

.. code-block:: python

    >>> output = recipe.run('map/script1')
    >>> output
    map/script1:
      DataFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.nc')
      DataFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.nc')
      ImageFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.png')
      ImageFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.png')

It is also possible to run a single preprocessor task:

.. code-block:: python

    >>> output = recipe.run('map/tas')
    >>> output
    map/tas:
      DataFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.nc')
      DataFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.nc')


API reference
*************

.. automodule:: esmvalcore.experimental.recipe
Find and download files from ESGF
=================================

This module provides the function :py:func:`esmvalcore.esgf.find_files`
for searching for files on ESGF using the ESMValTool vocabulary.
It returns :py:class:`esmvalcore.esgf.ESGFFile` objects, which have a convenient
:py:meth:`esmvalcore.esgf.ESGFFile.download` method for downloading the files.

See :ref:`config-esgf` for instructions on configuring this module.

esmvalcore.esgf
---------------
.. autofunction:: esmvalcore.esgf.find_files
.. autofunction:: esmvalcore.esgf.download
.. autoclass:: esmvalcore.esgf.ESGFFile

esmvalcore.esgf.facets
----------------------
.. automodule:: esmvalcore.esgf.facets
.. _experimental_api:

Experimental API
================

This page describes the new ESMValCore API.
The API module is available in the submodule ``esmvalcore.experimental``.
The API is under development, so use at your own risk!

.. toctree::

   esmvalcore.api.config
   esmvalcore.api.recipe
   esmvalcore.api.recipe_output
   esmvalcore.api.recipe_metadata
   esmvalcore.api.utils
.. _api_utils:



Utils
=====

This section describes the :py:class:`~esmvalcore.experimental.utils` submodule of the API (:py:mod:`esmvalcore.experimental`).


Finding recipes
***************

One of the first thing we may want to do, is to simply get one of the recipes available in ``ESMValTool``

If you already know which recipe you want to load, call :py:func:`~esmvalcore.experimental.utils.get_recipe`.

.. code-block:: python

    from esmvalcore.experimental import get_recipe
    >>> get_recipe('examples/recipe_python')
    Recipe('Recipe python')

Call the :py:func:`~esmvalcore.experimental.utils.get_all_recipes` function to get a list of all available recipes.

.. code-block:: python

    >>> from esmvalcore.experimental import get_all_recipes
    >>> recipes = get_all_recipes()
    >>> recipes
    [Recipe('Recipe perfmetrics cmip5 4cds'),
     Recipe('Recipe martin18grl'),
     ...
     Recipe('Recipe wflow'),
     Recipe('Recipe pcrglobwb')]

To search for a specific recipe, you can use the :py:meth:`~esmvalcore.experimental.utils.RecipeList.find` method. This takes a search query that looks through the recipe metadata and returns any matches. The query can be a regex pattern, so you can make it as complex as you like.

.. code-block:: python

    >>> results = recipes.find('climwip')
    [Recipe('Recipe climwip')]

The recipes are loaded in a :py:class:`~esmvalcore.experimental.recipe.Recipe` object, which knows about the documentation, authors, project, and related references of the recipe. It resolves all the tags, so that it knows which institute an author belongs to and which references are associated with the recipe.

This means you can search for something like this:

.. code-block:: python

    >>> recipes.find('Geophysical Research Letters')
    [Recipe('Recipe martin18grl'),
     Recipe('Recipe climwip'),
     Recipe('Recipe ecs constraints'),
     Recipe('Recipe ecs scatter'),
     Recipe('Recipe ecs'),
     Recipe('Recipe seaice')]


API reference
*************

.. automodule:: esmvalcore.experimental.utils
    :no-inherited-members:
    :no-show-inheritance:
.. _api_config:

Configuration
=============

This section describes the :py:class:`~esmvalcore.experimental.config` submodule of the API (:py:mod:`esmvalcore.experimental`).

Config
******

Configuration of ESMValCore/Tool is done via the :py:class:`~esmvalcore.experimental.config.Config` object.
The global configuration can be imported from the :py:mod:`esmvalcore.experimental` module as :py:data:`~esmvalcore.experimental.CFG`:

.. code-block:: python

    >>> from esmvalcore.experimental import CFG
    >>> CFG
    Config({'auxiliary_data_dir': PosixPath('/home/user/auxiliary_data'),
            'compress_netcdf': False,
            'config_developer_file': None,
            'config_file': PosixPath('/home/user/.esmvaltool/config-user.yml'),
            'drs': {'CMIP5': 'default', 'CMIP6': 'default'},
            'exit_on_warning': False,
            'log_level': 'info',
            'max_parallel_tasks': None,
            'output_dir': PosixPath('/home/user/esmvaltool_output'),
            'output_file_type': 'png',
            'profile_diagnostic': False,
            'remove_preproc_dir': True,
            'rootpath': {'CMIP5': '~/default_inputpath',
                         'CMIP6': '~/default_inputpath',
                         'default': '~/default_inputpath'},
            'save_intermediary_cubes': False)

The parameters for the user configuration file are listed :ref:`here <user configuration file>`.

:py:data:`~esmvalcore.experimental.CFG` is essentially a python dictionary with a few extra functions, similar to :py:mod:`matplotlib.rcParams`.
This means that values can be updated like this:

.. code-block:: python

    >>> CFG['output_dir'] = '~/esmvaltool_output'
    >>> CFG['output_dir']
    PosixPath('/home/user/esmvaltool_output')

Notice that :py:data:`~esmvalcore.experimental.CFG` automatically converts the path to an instance of ``pathlib.Path`` and expands the home directory.
All values entered into the config are validated to prevent mistakes, for example, it will warn you if you make a typo in the key:

.. code-block:: python

    >>> CFG['output_directory'] = '~/esmvaltool_output'
    InvalidConfigParameter: `output_directory` is not a valid config parameter.

Or, if the value entered cannot be converted to the expected type:

.. code-block:: python

    >>> CFG['max_parallel_tasks'] = '🐜'
    InvalidConfigParameter: Key `max_parallel_tasks`: Could not convert '🐜' to int

:py:class:`~esmvalcore.experimental.config.Config` is also flexible, so it tries to correct the type of your input if possible:

.. code-block:: python

    >>> CFG['max_parallel_tasks'] = '8'  # str
    >>> type(CFG['max_parallel_tasks'])
    int

By default, the config is loaded from the default location (``/home/user/.esmvaltool/config-user.yml``).
If it does not exist, it falls back to the default values.
to load a different file:

.. code-block:: python

    >>> CFG.load_from_file('~/my-config.yml')

Or to reload the current config:

.. code-block:: python

    >>> CFG.reload()


Session
*******

Recipes and diagnostics will be run in their own directories.
This behaviour can be controlled via the :py:data:`~esmvalcore.experimental.config.Session` object.
A :py:data:`~esmvalcore.experimental.config.Session` can be initiated from the global :py:class:`~esmvalcore.experimental.config.Config`.

.. code-block:: python

    >>> session = CFG.start_session(name='my_session')

A :py:data:`~esmvalcore.experimental.config.Session` is very similar to the config.
It is also a dictionary, and copies all the keys from the :py:class:`~esmvalcore.experimental.config.Config`.
At this moment, ``session`` is essentially a copy of :py:data:`~esmvalcore.experimental.CFG`:

.. code-block:: python

    >>> print(session == CFG)
    True
    >>> session['output_dir'] = '~/my_output_dir'
    >>> print(session == CFG)  # False
    False

A :py:data:`~esmvalcore.experimental.config.Session` also knows about the directories where the data will stored.
The session name is used to prefix the directories.

.. code-block:: python

    >>> session.session_dir
    /home/user/my_output_dir/my_session_20201203_155821
    >>> session.run_dir
    /home/user/my_output_dir/my_session_20201203_155821/run
    >>> session.work_dir
    /home/user/my_output_dir/my_session_20201203_155821/work
    >>> session.preproc_dir
    /home/user/my_output_dir/my_session_20201203_155821/preproc
    >>> session.plot_dir
    /home/user/my_output_dir/my_session_20201203_155821/plots

Unlike the global configuration, of which only one can exist, multiple sessions can be initiated from :py:class:`~esmvalcore.experimental.config.Config`.


API reference
*************

.. automodule:: esmvalcore.experimental.config
    :no-inherited-members:
    :no-show-inheritance:
.. _api_recipe_metadata:

Recipe Metadata
===============

This section describes the :py:mod:`~esmvalcore.experimental.recipe_metadata` submodule of the API (:py:mod:`esmvalcore.experimental`).

API reference
*************

.. automodule:: esmvalcore.experimental.recipe_metadata
.. _preprocessor_functions:

Preprocessor functions
======================

.. autodata:: esmvalcore.preprocessor.DEFAULT_ORDER
.. automodule:: esmvalcore.preprocessor
.. _api_recipe_output:

Recipe output
=============

This section describes the :py:mod:`~esmvalcore.experimental.recipe_output` submodule of the API (:py:mod:`esmvalcore.experimental`).

After running a recipe, output is returned by the :py:meth:`~esmvalcore.experimental.recipe.Recipe.run` method. Alternatively, it can be retrieved using the :py:meth:`~esmvalcore.experimental.recipe.Recipe.get_output` method.

.. code:: python

    >>> recipe_output = recipe.get_output()

``recipe_output`` is a mapping of the individual tasks and their output
filenames (data and image files) with a set of attributes describing the
data.

.. code:: python

    >>> recipe_output
    timeseries/script1:
      DataFile('tas_amsterdam_CMIP5_CanESM2_Amon_historical_r1i1p1_tas_1850-2000.nc')
      DataFile('tas_amsterdam_CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_1850-2000.nc')
      DataFile('tas_amsterdam_MultiModelMean_Amon_tas_1850-2000.nc')
      DataFile('tas_global_CMIP5_CanESM2_Amon_historical_r1i1p1_tas_1850-2000.nc')
      DataFile('tas_global_CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_1850-2000.nc')
      ImageFile('tas_amsterdam_CMIP5_CanESM2_Amon_historical_r1i1p1_tas_1850-2000.png')
      ImageFile('tas_amsterdam_CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_1850-2000.png')
      ImageFile('tas_amsterdam_MultiModelMean_Amon_tas_1850-2000.png')
      ImageFile('tas_global_CMIP5_CanESM2_Amon_historical_r1i1p1_tas_1850-2000.png')
      ImageFile('tas_global_CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_1850-2000.png')

    map/script1:
      DataFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.nc')
      DataFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.nc')
      ImageFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.png')
      ImageFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.png')


Output is grouped by the task that produced them. They can be accessed like
a dictionary.

.. code:: python

    >>> task_output = recipe_output['map/script1']
    >>> task_output
    map/script1:
      DataFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.nc')
      DataFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.nc')
      ImageFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.png')
      ImageFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.png')


The task output has a list of files associated with them, usually image
(``.png``) or data files (``.nc``). To get a list of all files, use
:py:meth:`~esmvalcore.experimental.recipe_output.TaskOutput.files`.

.. code:: python

    >>> print(task_output.files)
    (DataFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.nc'),
    ..., ImageFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.png'))


It is also possible to select the image (:py:meth:`~esmvalcore.experimental.recipe_output.TaskOutput.image_files`) files or data files (:py:meth:`~esmvalcore.experimental.recipe_output.TaskOutput.data_files`) only.

.. code:: python

    >>> for image_file in task_output.image_files:
    >>>     print(image_file)
    ImageFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.png')
    ImageFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.png')

    >>> for data_file in task_output.data_files:
    >>>     print(data_file)
    DataFile('CMIP5_CanESM2_Amon_historical_r1i1p1_tas_2000-2000.nc')
    DataFile('CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_2000-2000.nc')


Working with output files
*************************

Output comes in two kinds, :py:class:`~esmvalcore.experimental.recipe_output.DataFile` corresponds to data
files in ``.nc`` format and :py:class:`~esmvalcore.experimental.recipe_output.ImageFile` corresponds to plots
in ``.png`` format (see below). Both object are derived from the same base class
(:py:class:`~esmvalcore.experimental.recipe_output.OutputFile`) and therefore share most of the functionality.

For example, author information can be accessed as instances of :py:class:`~esmvalcore.experimental.recipe_metadata.Contributor`  via

.. code:: python

    >>> output_file = task_output[0]
    >>> output_file.authors
    (Contributor('Andela, Bouwe', institute='NLeSC, Netherlands', orcid='https://orcid.org/0000-0001-9005-8940'),
     Contributor('Righi, Mattia', institute='DLR, Germany', orcid='https://orcid.org/0000-0003-3827-5950'))

And associated references as instances of :py:class:`~esmvalcore.experimental.recipe_metadata.Reference` via

.. code:: python

    >>> output_file.references
    (Reference('acknow_project'),)

:py:class:`~esmvalcore.experimental.recipe_output.OutputFile` also knows about associated files

.. code:: python

    >>> data_file.citation_file
    Path('.../tas_global_CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_1850-2000_citation.bibtex')
    >>> data_file.data_citation_file
    Path('.../tas_global_CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_1850-2000_data_citation_info.txt')
    >>> data_file.provenance_svg_file
    Path('.../tas_global_CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_1850-2000_provenance.svg')
    >>> data_file.provenance_xml_file
    Path('.../tas_global_CMIP6_BCC-ESM1_Amon_historical_r1i1p1f1_tas_1850-2000_provenance.xml')



Working with image files
************************

Image output uses IPython magic to plot themselves in a notebook
environment.

.. code:: python

    >>> image_file = recipe_output['map/script1'].image_files[0]
    >>> image_file

For example:

.. image:: /figures/api_recipe_output.png
   :width: 600

Using :py:mod:`IPython.display`, it is possible to show all image files.

.. code:: python

    >>> from IPython.display import display
    >>>
    >>> task = recipe_output['map/script1']
    >>> for image_file in task.image_files:
    >>>      display(image_file)


Working with data files
***********************

Data files can be easily loaded using ``xarray``:

.. code:: python

    >>> data_file = recipe_output['timeseries/script1'].data_files[0]
    >>> data = data_file.load_xarray()
    >>> type(data)
    xarray.core.dataset.Dataset


Or ``iris``:

.. code:: python

    >>> cube = data_file.load_iris()
    >>> type(cube)
    iris.cube.CubeList


API reference
*************

.. automodule:: esmvalcore.experimental.recipe_output
CMOR functions
==============

.. automodule:: esmvalcore.cmor

Checking compliance
-------------------

.. automodule:: esmvalcore.cmor.check

Automatically fixing issues
---------------------------

.. automodule:: esmvalcore.cmor.fix

Functions for fixing issues
---------------------------

.. automodule:: esmvalcore.cmor.fixes

Using CMOR tables
-----------------

.. automodule:: esmvalcore.cmor.table
.. _api:

ESMValCore API Reference
========================

ESMValCore is mostly used as a commandline tool. However, it is also possibly to use (parts of) ESMValTool as a
library. This section documents the public API of ESMValCore.

.. toctree::

   esmvalcore.cmor
   esmvalcore.esgf
   esmvalcore.exceptions
   esmvalcore.preprocessor
   esmvalcore.api
Exceptions
==========

.. automodule:: esmvalcore.exceptions
    :no-inherited-members:
.. _recipes:

Working with the installed recipes
**********************************

Although ESMValTool can be used just to simplify the management of data
and the creation of your own analysis code, one of its main strengths is the
continuously growing set of diagnostics and metrics that it directly provides to
the user. These metrics and diagnostics are provided as a set of preconfigured
recipes that users can run or customize for their own analysis.
The latest list of available recipes can be found :ref:`here <esmvaltool:recipes>`.

In order to make the management of these installed recipes easier, ESMValTool
provides the ``recipes`` command group with utilities that help the users in
discovering and customizing the provided recipes.

The first command in this group allows users to get the complete list of installed
recipes printed to the console:

.. code:: bash

	esmvaltool recipes list

If the user then wants to explore any one of this recipes, they can be printed
using the following command

.. code:: bash

	esmvaltool recipes show recipe_name.yml

And finally, to get a local copy that can then be customized and run, users can
use the following command

.. code:: bash

	esmvaltool recipes get recipe_name.yml
.. _outputdata:

Output
******

ESMValTool automatically generates a new output directory with every run. The
location is determined by the output_dir option  in the config-user.yml file,
the recipe name, and the date and time, using the the format: ``YYYYMMDD_HHMMSS``.

For instance, a typical output location would be:
``output_directory/recipe_ocean_amoc_20190118_1027/``

This is effectively produced by the combination:
``output_dir/recipe_name_YYYYMMDD_HHMMSS/``

This directory will contain 4 further subdirectories:

1. `Diagnostic output`_ (``work``): A place for any diagnostic script results that are not plots, e.g. files in NetCDF format (depends on the diagnostics).

2. `Plots`_ (``plots``): The location for all the plots, split by individual diagnostics and fields.

3. `Run`_ (``run``): This directory includes all log files, a copy of the recipe, a summary of the resource usage, and the `settings.yml`_ interface files and temporary files created by the diagnostic scripts.

4. `Preprocessed datasets`_ (``preproc``): This directory contains all the preprocessed netcdfs data and the `metadata.yml`_ interface files. Note that by default this directory will be deleted after each run, because most users will only need the results from the diagnostic scripts.

A summary of the output is produced in the file:
``index.html``


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

The settings.yml file is automatically generated by ESMValTool. Each diagnostic
will produce a unique settings.yml file.

The settings.yml file passes several global level keys to diagnostic scripts.
This includes several flags from the config-user.yml file (such as
'log_level'), several paths which are specific to the
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

The first item in the settings file will be a list of `Metadata.yml`_ files.
There is a metadata.yml file generated for each field in each diagnostic.


Metadata.yml
============

The metadata.yml files is automatically generated by ESMValTool. Along with the
settings.yml file, it passes all the paths, boolean flags, and additional
arguments that your diagnostic needs to know in order to run.

The metadata is loaded from cfg as a dictionairy object in python diagnostics.

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
these files. The tools are available in the shared directory in the diagnostics
directory.
.. _findingdata:

************
Input data
************

Overview
========
Data discovery and retrieval is the first step in any evaluation process;
ESMValTool uses a `semi-automated` data finding mechanism with inputs from both
the user configuration file and the recipe file: this means that the user will
have to provide the tool with a set of parameters related to the data needed
and once these parameters have been provided, the tool will automatically find
the right data. We will detail below the data finding and retrieval process and
the input the user needs to specify, giving examples on how to use the data
finding routine under different scenarios.

Data types
==========

.. _CMOR-DRS:

CMIP data
---------
CMIP data is widely available via the Earth System Grid Federation
(`ESGF <https://esgf.llnl.gov/>`_) and is accessible to users either
via automatic download by ``esmvaltool`` or through the ESGF data nodes hosted
by large computing facilities (like CEDA-Jasmin, DKRZ, etc). This data
adheres to, among other standards, the DRS and Controlled Vocabulary
standard for naming files and structured paths; the `DRS
<https://www.ecmwf.int/sites/default/files/elibrary/2014/13713-data-reference-syntax-governing-standards-within-climate-research-data-archived-esgf.pdf>`_
ensures that files and paths to them are named according to a
standardized convention. Examples of this convention, also used by
ESMValTool for file discovery and data retrieval, include:

* CMIP6 file: ``[variable_short_name]_[mip]_[dataset_name]_[experiment]_[ensemble]_[grid]_[start-date]-[end-date].nc``
* CMIP5 file: ``[variable_short_name]_[mip]_[dataset_name]_[experiment]_[ensemble]_[start-date]-[end-date].nc``
* OBS file: ``[project]_[dataset_name]_[type]_[version]_[mip]_[short_name]_[start-date]-[end-date].nc``

Similar standards exist for the standard paths (input directories); for the
ESGF data nodes, these paths differ slightly, for example:

* CMIP6 path for BADC: ``ROOT-BADC/[institute]/[dataset_name]/[experiment]/[ensemble]/[mip]/
  [variable_short_name]/[grid]``;
* CMIP6 path for ETHZ: ``ROOT-ETHZ/[experiment]/[mip]/[variable_short_name]/[dataset_name]/[ensemble]/[grid]``

From the ESMValTool user perspective the number of data input parameters is
optimized to allow for ease of use. We detail this procedure in the next
section.

Native model data
-----------------
Support for native model data that is not formatted according to a CMIP
data request is quite easy using basic
:ref:`ESMValCore fix procedure <fixing_data>` and has been implemented
for some models :ref:`as described here <fixing_native_models>`

Observational data
------------------
Part of observational data is retrieved in the same manner as CMIP data, for example
using the ``OBS`` root path set to:

  .. code-block:: yaml

    OBS: /gws/nopw/j04/esmeval/obsdata-v2

and the dataset:

  .. code-block:: yaml

    - {dataset: ERA-Interim, project: OBS, type: reanaly, version: 1, start_year: 2014, end_year: 2015, tier: 3}

in ``recipe.yml`` in ``datasets`` or ``additional_datasets``, the rules set in
CMOR-DRS_ are used again and the file will be automatically found:

.. code-block::

  /gws/nopw/j04/esmeval/obsdata-v2/Tier3/ERA-Interim/OBS_ERA-Interim_reanaly_1_Amon_ta_201401-201412.nc

Since observational data are organized in Tiers depending on their level of
public availability, the ``default`` directory must be structured accordingly
with sub-directories ``TierX`` (``Tier1``, ``Tier2`` or ``Tier3``), even when
``drs: default``.

.. _data-retrieval:

Data retrieval
==============
Data retrieval in ESMValTool has two main aspects from the user's point of
view:

* data can be found by the tool, subject to availability on disk or `ESGF <https://esgf.llnl.gov/>`_;
* it is the user's responsibility to set the correct data retrieval parameters;

The first point is self-explanatory: if the user runs the tool on a machine
that has access to a data repository or multiple data repositories, then
ESMValTool will look for and find the available data requested by the user.
If the files are not found locally, the tool can search the ESGF_ and download
the missing files, provided that they are available.

The second point underlines the fact that the user has full control over what
type and the amount of data is needed for the analyses. Setting the data
retrieval parameters is explained below.

Enabling automatic downloads from the ESGF
------------------------------------------
To enable automatic downloads from ESGF, set ``offline: false`` in
the :ref:`user configuration file` or provide the command line argument
``--offline=False`` when running the recipe.
The files will be stored in the ``download_dir`` set in
the :ref:`user configuration file`.

Setting the correct root paths
------------------------------
The first step towards providing ESMValTool the correct set of parameters for
data retrieval is setting the root paths to the data. This is done in the user
configuration file ``config-user.yml``. The two sections where the user will
set the paths are ``rootpath`` and ``drs``. ``rootpath`` contains pointers to
``CMIP``, ``OBS``, ``default`` and ``RAWOBS`` root paths; ``drs`` sets the type
of directory structure the root paths are structured by. It is important to
first discuss the ``drs`` parameter: as we've seen in the previous section, the
DRS as a standard is used for both file naming conventions and for directory
structures.

Synda
-----

If the `synda install <https://prodiguer.github.io/synda/sdt/user_guide.html#synda-install>`_ command is used to download data,
it maintains the directory structure as on ESGF. To find data downloaded by
synda, use the ``SYNDA`` ``drs`` parameter.

.. code-block:: yaml

 drs:
   CMIP6: SYNDA
   CMIP5: SYNDA

.. _config-user-drs:

Explaining ``config-user/drs: CMIP5:`` or ``config-user/drs: CMIP6:``
---------------------------------------------------------------------
Whereas ESMValTool will **always** use the CMOR standard for file naming (please
refer above), by setting the ``drs`` parameter the user tells the tool what
type of root paths they need the data from, e.g.:

  .. code-block:: yaml

   drs:
     CMIP6: BADC

will tell the tool that the user needs data from a repository structured
according to the BADC DRS structure, i.e.:

``ROOT/[institute]/[dataset_name]/[experiment]/[ensemble]/[mip]/[variable_short_name]/[grid]``;

setting the ``ROOT`` parameter is explained below. This is a
strictly-structured repository tree and if there are any sort of irregularities
(e.g. there is no ``[mip]`` directory) the data will not be found! ``BADC`` can
be replaced with ``DKRZ`` or ``ETHZ`` depending on the existing ``ROOT``
directory structure.
The snippet

  .. code-block:: yaml

   drs:
     CMIP6: default

is another way to retrieve data from a ``ROOT`` directory that has no DRS-like
structure; ``default`` indicates that the data lies in a directory that
contains all the files without any structure.

.. note::
   When using ``CMIP6: default`` or ``CMIP5: default`` it is important to
   remember that all the needed files must be in the same top-level directory
   set by ``default`` (see below how to set ``default``).

.. _config-user-rootpath:

Explaining ``config-user/rootpath:``
------------------------------------

``rootpath`` identifies the root directory for different data types (``ROOT`` as we used it above):

* ``CMIP`` e.g. ``CMIP5`` or ``CMIP6``: this is the `root` path(s) to where the
  CMIP files are stored; it can be a single path or a list of paths; it can
  point to an ESGF node or it can point to a user private repository. Example
  for a CMIP5 root path pointing to the ESGF node on CEDA-Jasmin (formerly
  known as BADC):

  .. code-block:: yaml

    CMIP5: /badc/cmip5/data/cmip5/output1

  Example for a CMIP6 root path pointing to the ESGF node on CEDA-Jasmin:

  .. code-block:: yaml

    CMIP6: /badc/cmip6/data/CMIP6/CMIP

  Example for a mix of CMIP6 root path pointing to the ESGF node on CEDA-Jasmin
  and a user-specific data repository for extra data:

  .. code-block:: yaml

    CMIP6: [/badc/cmip6/data/CMIP6/CMIP, /home/users/johndoe/cmip_data]

* ``OBS``: this is the `root` path(s) to where the observational datasets are
  stored; again, this could be a single path or a list of paths, just like for
  CMIP data. Example for the OBS path for a large cache of observation datasets
  on CEDA-Jasmin:

  .. code-block:: yaml

    OBS: /gws/nopw/j04/esmeval/obsdata-v2

* ``default``: this is the `root` path(s) where the tool will look for data
  from projects that do not have their own rootpath set.

* ``RAWOBS``: this is the `root` path(s) to where the raw observational data
  files are stored; this is used by ``cmorize_obs``.

Dataset definitions in ``recipe``
---------------------------------
Once the correct paths have been established, ESMValTool collects the
information on the specific datasets that are needed for the analysis. This
information, together with the CMOR convention for naming files (see CMOR-DRS_)
will allow the tool to search and find the right files. The specific
datasets are listed in any recipe, under either the ``datasets`` and/or
``additional_datasets`` sections, e.g.

.. code-block:: yaml

  datasets:
    - {dataset: HadGEM2-CC, project: CMIP5, exp: historical, ensemble: r1i1p1, start_year: 2001, end_year: 2004}
    - {dataset: UKESM1-0-LL, project: CMIP6, exp: historical, ensemble: r1i1p1f2, grid: gn, start_year: 2004, end_year: 2014}

``_data_finder`` will use this information to find data for **all** the variables specified in ``diagnostics/variables``.

Recap and example
=================
Let us look at a practical example for a recap of the information above:
suppose you are using a ``config-user.yml`` that has the following entries for
data finding:

.. code-block:: yaml

  rootpath:  # running on CEDA-Jasmin
    CMIP6: /badc/cmip6/data/CMIP6/CMIP
  drs:
    CMIP6: BADC  # since you are on CEDA-Jasmin

and the dataset you need is specified in your ``recipe.yml`` as:

.. code-block:: yaml

  - {dataset: UKESM1-0-LL, project: CMIP6, mip: Amon, exp: historical, grid: gn, ensemble: r1i1p1f2, start_year: 2004, end_year: 2014}

for a variable, e.g.:

.. code-block:: yaml

  diagnostics:
    some_diagnostic:
      description: some_description
      variables:
        ta:
          preprocessor: some_preprocessor

The tool will then use the root path ``/badc/cmip6/data/CMIP6/CMIP`` and the
dataset information and will assemble the full DRS path using information from
CMOR-DRS_ and establish the path to the files as:

.. code-block:: bash

  /badc/cmip6/data/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Amon

then look for variable ``ta`` and specifically the latest version of the data
file:

.. code-block:: bash

  /badc/cmip6/data/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Amon/ta/gn/latest/

and finally, using the file naming definition from CMOR-DRS_ find the file:

.. code-block:: bash

  /badc/cmip6/data/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Amon/ta/gn/latest/ta_Amon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc

.. _observations:


Data loading
============

Data loading is done using the data load functionality of `iris`; we will not go into too much detail
about this since we can point the user to the specific functionality
`here <https://scitools-iris.readthedocs.io/en/latest/userguide/loading_iris_cubes.html>`_ but we will underline
that the initial loading is done by adhering to the CF Conventions that `iris` operates by as well (see
`CF Conventions Document <http://cfconventions.org/cf-conventions/cf-conventions.html>`_ and the search
page for CF `standard names <http://cfconventions.org/standard-names.html>`_).

Data concatenation from multiple sources
========================================

Oftentimes data retrieving results in assembling a continuous data stream from
multiple files or even, multiple experiments. The internal mechanism through which
the assembly is done is via cube concatenation. One peculiarity of iris concatenation
(see `iris cube concatenation <https://scitools-iris.readthedocs.io/en/latest/userguide/merge_and_concat.html>`_)
is that it doesn't allow for concatenating time-overlapping cubes; this case is rather
frequent with data from models overlapping in time, and is accounted for by a function that performs a
flexible concatenation between two cubes, depending on the particular setup:

* cubes overlap in time: resulting cube is made up of the overlapping data plus left and
  right hand sides on each side of the overlapping data; note that in the case of the cubes
  coming from different experiments the resulting concatenated cube will have composite data
  made up from multiple experiments: assume [cube1: exp1, cube2: exp2] and cube1 starts before cube2,
  and cube2 finishes after cube1, then the concatenated cube will be made up of cube2: exp2 plus the
  section of cube1: exp1 that contains data not provided in cube2: exp2;
* cubes don't overlap in time: data from the two cubes is bolted together;

Note that two cube concatenation is the base operation of an iterative process of reducing multiple cubes
from multiple data segments via cube concatenation ie if there is no time-overlapping data, the
cubes concatenation is performed in one step.

.. _extra-facets-data-finder:

Use of extra facets in the datafinder
=====================================
Extra facets are a mechanism to provide additional information for certain kinds
of data. The general approach is described in :ref:`extra_facets`. Here, we
describe how they can be used to locate data files within the datafinder
framework.
This is useful to build paths for directory structures and file names
that require more information than what is provided in the recipe.
A common application is the location of variables in multi-variable files as
often found in climate models' native output formats.

Another use case is files that use different names for variables in their
file name than for the netCDF4 variable name.

To apply the extra facets for this purpose, simply use the corresponding tag in
the applicable DRS inside the `config-developer.yml` file. For example, given
the extra facets in :ref:`extra-facets-example-1`, one might write the
following.

.. _extra-facets-example-2:

.. code-block:: yaml
   :caption: Example drs use in `config-developer.yml`

   native6:
     input_file:
       default: '{name_in_filename}*.nc'

The same replacement mechanism can be employed everywhere where tags can be
used, particularly in `input_dir` and `input_file`.
.. _running:

Running
*******

The ESMValCore package provides the ``esmvaltool`` command line tool, which can
be used to run a :doc:`recipe <../recipe/index>`.

To run a recipe, call ``esmvaltool run`` with the desired recipe:

.. code:: bash

	esmvaltool run recipe_python.yml

If the configuration file is not in the default location
``~/.esmvaltool/config-user.yml``, you can pass its path explicitly:

.. code:: bash

	esmvaltool run --config_file /path/to/config-user.yml recipe_python.yml

It is also possible to explicitly change values from the config file using flags:

.. code:: bash

	esmvaltool run --argument_name argument_value recipe_python.yml

To automatically download the files required to run a recipe from ESGF, set
``offline`` to ``false`` in the :ref:`user configuration file`
or run the tool with the command

.. code:: bash

    esmvaltool run --offline=False recipe_python.yml

This feature is available for projects that are hosted on the ESGF, i.e.
CMIP3, CMIP5, CMIP6, CORDEX, and obs4MIPs.

To control the strictness of the CMOR checker, use the flag ``--check_level``:

.. code:: bash

	esmvaltool run --check_level=relaxed recipe_python.yml

Possible values are:

  - `ignore`: all errors will be reported as warnings
  - `relaxed`: only fail if there are critical errors
  - `default`: fail if there are any errors
  - `strict`: fail if there are any warnings

To re-use pre-processed files from a previous run of the same recipe, you can
use

.. code:: bash

    esmvaltool run recipe_python.yml --resume_from ~/esmvaltool_output/recipe_python_20210930_123907

Multiple directories can be specified for re-use, make sure to quote them:

.. code:: bash

    esmvaltool run recipe_python.yml --resume_from "~/esmvaltool_output/recipe_python_20210930_101007 ~/esmvaltool_output/recipe_python_20210930_123907"

The first preprocessor directory containing the required data will be used.

This feature can be useful when developing new diagnostics, because it avoids
the need to re-run the preprocessor.
Another potential use case is running the preprocessing part of a recipe on
one or more machines that have access to a lot of data and then running the
diagnostics on a machine without access to data.

To run only the preprocessor tasks from a recipe, use

.. code:: bash

    esmvaltool run recipe_python.yml --remove_preproc_dir=False --run_diagnostic=False

.. note::

    Only preprocessing :ref:`tasks <tasks>` that completed successfully
    can be re-used with the ``--resume_from`` option.
    Preprocessing tasks that completed successfully, contain a file called
    :ref:`metadata.yml <interface_esmvalcore_diagnostic>` in their output
    directory.

To run a reduced version of the recipe, usually for testing purpose you can use

.. code:: bash

	esmvaltool run --max_datasets=NDATASETS --max_years=NYEARS recipe_python.yml

In this case, the recipe will limit the number of datasets per variable to
NDATASETS and the total amount of years loaded to NYEARS. They can also be used
separately.
Note that diagnostics may require specific combinations of available data, so
use the above two flags at your own risk and for testing purposes only.

To run a recipe, even if some datasets are not available, use

.. code:: bash

    esmvaltool run --skip_nonexistent=True recipe_python.yml

It is also possible to select only specific diagnostics to be run. To tun only
one, just specify its name. To provide more than one diagnostic to filter use
the syntax 'diag1 diag2/script1' or '("diag1", "diag2/script1")' and pay
attention to the quotes.

.. code:: bash

    esmvaltool run --diagnostics=diagnostic1 recipe_python.yml



To get help on additional commands, please use

.. code:: bash

	esmvaltool --help



.. note::

	ESMValTool command line interface is created using the Fire python package.
	This package supports the creation of completion scripts for the Bash and
	Fish shells. Go to https://google.github.io/python-fire/using-cli/#python-fires-flags
	to learn how to set up them.
Getting started
***************

.. toctree::
   :maxdepth: 1

		Installation <install>
    Configuration <configure>
    Input data <find_data>
    Installed recipes <recipes>
		Running <run>
		Output <output>
.. _config:

*******************
Configuration files
*******************

Overview
========

There are several configuration files in ESMValCore:

* ``config-user.yml``: sets a number of user-specific options like desired
  graphical output format, root paths to data, etc.;
* ``config-developer.yml``: sets a number of standardized file-naming and paths
  to data formatting;

and one configuration file which is distributed with ESMValTool:

* ``config-references.yml``: stores information on diagnostic and recipe authors and
  scientific journals references;

.. _user configuration file:

User configuration file
=======================


The ``config-user.yml`` configuration file contains all the global level
information needed by ESMValTool. It can be reused as many times the user needs
to before changing any of the options stored in it. This file is essentially
the gateway between the user and the machine-specific instructions to
``esmvaltool``. By default, esmvaltool looks for it in the home directory,
inside the ``.esmvaltool`` folder.

Users can get a copy of this file with default values by running

.. code-block:: bash

  esmvaltool config get-config-user --path=${TARGET_FOLDER}

If the option ``--path`` is omitted, the file will be created in
``${HOME}/.esmvaltool``

The following shows the default settings from the ``config-user.yml`` file
with explanations in a commented line above each option:

.. code-block:: yaml

  # Destination directory where all output will be written
  # including log files and performance stats
  output_dir: ./esmvaltool_output

  # Directory for storing downloaded climate data
  download_dir: ~/climate_data

  # Disable the automatic download of missing CMIP3, CMIP5, CMIP6, CORDEX,
  # and obs4MIPs data from ESGF [true]/false.
  offline: true

  # Auxiliary data directory, used by some recipes to look for additional datasets
  auxiliary_data_dir: ~/auxiliary_data

  # Rootpaths to the data from different projects (lists are also possible)
  # these are generic entries to better allow you to enter your own
  # For site-specific entries, see the default config-user.yml file
  # that can be installed with the command `esmvaltool config get_config_user`.
  # For each project, this can be either a single path or a list of paths.
  # Comment out these when using a site-specific path
  rootpath:
    default: ~/climate_data

  # Directory structure for input data: [default]/ESGF/BADC/DKRZ/ETHZ/etc
  # See config-developer.yml for definitions.
  # comment out/replace as per needed
  drs:
    CMIP3: ESGF
    CMIP5: ESGF
    CMIP6: ESGF
    CORDEX: ESGF
    obs4MIPs: ESGF

  # Run at most this many tasks in parallel [null]/1/2/3/4/..
  # Set to null to use the number of available CPUs.
  # If you run out of memory, try setting max_parallel_tasks to 1 and check the
  # amount of memory you need for that by inspecting the file
  # run/resource_usage.txt in the output directory. Using the number there you
  # can increase the number of parallel tasks again to a reasonable number for
  # the amount of memory available in your system.
  max_parallel_tasks: null

  # Set the console log level debug, [info], warning, error
  # for much more information printed to screen set log_level: debug
  log_level: info

  # Exit on warning (only for NCL diagnostic scripts)? true/[false]
  exit_on_warning: false

  # Plot file format? [png]/pdf/ps/eps/epsi
  output_file_type: png

  # Remove the ``preproc`` dir if the run was successful
  # By default this option is set to "true", so all preprocessor output files
  # will be removed after a successful run. Set to "false" if you need those files.
  remove_preproc_dir: true

  # Use netCDF compression true/[false]
  compress_netcdf: false

  # Save intermediary cubes in the preprocessor true/[false]
  # set to true will save the output cube from each preprocessing step
  # these files are numbered according to the preprocessing order
  save_intermediary_cubes: false

  # Use a profiling tool for the diagnostic run [false]/true
  # A profiler tells you which functions in your code take most time to run.
  # For this purpose we use vprof, see below for notes
  # Only available for Python diagnostics
  profile_diagnostic: false

  # Path to custom config-developer file, to customise project configurations.
  # See config-developer.yml for an example. Set to "null" to use the default
  config_developer_file: null

.. code-block:: yaml

  # Auxiliary data directory (used for some additional datasets)
  auxiliary_data_dir: ~/auxiliary_data

The ``offline`` setting can be used to disable or enable automatic downloads from ESGF.
If ``offline`` is set to ``false``, the tool will automatically download
any CMIP3, CMIP5, CMIP6, CORDEX, and obs4MIPs data that is required to run a recipe
but not available locally and store it in ``download_dir`` using the ``ESGF``
directory structure defined in the :ref:`config-developer`.

The ``auxiliary_data_dir`` setting is the path to place any required
additional auxiliary data files. This is necessary because certain
Python toolkits, such as cartopy, will attempt to download data files at run
time, typically geographic data files such as coastlines or land surface maps.
This can fail if the machine does not have access to the wider internet. This
location allows the user to specify where to find such files if they can not be
downloaded at runtime. The example user configuration file already contains two valid
locations for ``auxiliary_data_dir`` directories on CEDA-JASMIN and DKRZ, and a number
of such maps and shapefiles (used by current diagnostics) are already there. You will
need ``esmeval`` group workspace membership to access the JASMIN one (see
`instructions <https://help.jasmin.ac.uk/article/199-introduction-to-group-workspaces>`_
how to gain access to the group workspace.

.. warning::

   This setting is not for model or observational datasets, rather it is for
   extra data files such as shapefiles or other data sources needed by the diagnostics.

The ``profile_diagnostic`` setting triggers profiling of Python diagnostics,
this will tell you which functions in the diagnostic took most time to run.
For this purpose we use `vprof <https://github.com/nvdv/vprof>`_.
For each diagnostic script in the recipe, the profiler writes a ``.json`` file
that can be used to plot a
`flame graph <https://queue.acm.org/detail.cfm?id=2927301>`__
of the profiling information by running

.. code-block:: bash

  vprof --input-file esmvaltool_output/recipe_output/run/diagnostic/script/profile.json

Note that it is also possible to use vprof to understand other resources used
while running the diagnostic, including execution time of different code blocks
and memory usage.

A detailed explanation of the data finding-related sections of the
``config-user.yml`` (``rootpath`` and ``drs``) is presented in the
:ref:`data-retrieval` section. This section relates directly to the data
finding capabilities  of ESMValTool and are very important to be understood by
the user.

.. note::

   You can choose your ``config-user.yml`` file at run time, so you could have several of
   them available with different purposes. One for a formalised run, another for
   debugging, etc. You can even provide any config user value as a run flag
   ``--argument_name argument_value``


.. _config-esgf:

ESGF configuration
==================

The ``esmvaltool run`` command can automatically download the files required
to run a recipe from ESGF for the projects CMIP3, CMIP5, CMIP6, CORDEX, and obs4MIPs.
The downloaded files will be stored in the ``download_dir`` specified in the
:ref:`user configuration file`.
To enable automatic downloads from ESGF, set ``offline: false`` in
the :ref:`user configuration file` or provide the command line argument
``--offline=False`` when running the recipe.

.. note::

   When running a recipe that uses many or large datasets on a machine that
   does not have any data available locally, the amount of data that will be
   downloaded can be in the range of a few hundred gigabyte to a few terrabyte.
   See :ref:`esmvaltool:inputdata` for advice on getting access to machines
   with large datasets already available.

   A log message will be displayed with the total amount of data that will
   be downloaded before starting the download.
   If you see that this is more than you would like to download, stop the
   tool by pressing the ``Ctrl`` and ``C`` keys on your keyboard simultaneously
   several times, edit the recipe so it contains fewer datasets and try again.

For downloading some files (e.g. those produced by the CORDEX project),
you need to log in to be able to download the data.

See the
`ESGF user guide <https://esgf.github.io/esgf-user-support/user_guide.html>`_
for instructions on how to create an ESGF OpenID account if you do not have
one yet.
Note that the OpenID account consists of 3 components instead of the usual
two, in addition a username and password you also need the hostname of the
provider of the ID; for example
`esgf-data.dkrz.de <https://esgf-data.dkrz.de/user/add/?next=http://esgf-data.dkrz.de/projects/esgf-dkrz/>`_.
Even though the account is issued by a particular host, the same OpenID
account can be used to download data from all hosts in the ESGF.

Next, configure your system so the ``esmvaltool`` can use your credentials.
This can be done using the keyring_ package or they can be stored in a
:ref:`configuration file <config_esgf_pyclient>`.

.. _keyring:

Storing credentials in keyring
------------------------------
First install the keyring package. Note that this requires a supported
backend that may not be available on compute clusters, see the
`keyring documentation <https://pypi.org/project/keyring>`__ for more
information.

.. code-block:: bash

    pip install keyring

Next, set your username and password by running the commands:

.. code-block:: bash

    keyring set ESGF hostname
    keyring set ESGF username
    keyring set ESGF password

for example, if you created an account on the host `esgf-data.dkrz.de`_ with username
'cookiemonster' and password 'Welcome01', run the command

.. code-block:: bash

    keyring set ESGF hostname

this will display the text

.. code-block:: bash

    Password for 'hostname' in 'ESGF':

type ``esgf-data.dkrz.de`` (the characters will not be shown) and press ``Enter``.
Repeat the same procedure with ``keyring set ESGF username``, type ``cookiemonster``
and press ``Enter`` and ``keyring set ESGF password``, type ``Welcome01`` and
press ``Enter``.

To check that you entered your credentials correctly, run:

.. code-block:: bash

    keyring get ESGF hostname
    keyring get ESGF username
    keyring get ESGF password

.. _config_esgf_pyclient:

Configuration file
------------------
An optional configuration file can be created for configuring how the tool uses
`esgf-pyclient <https://esgf-pyclient.readthedocs.io>`_
to find and download data.
The name of this file is ``~/.esmvaltool/esgf-pyclient.yml``.

Logon
`````
In the ``logon`` section you can provide arguments that will be passed on to
:py:meth:`pyesgf.logon.LogonManager.logon`.
For example, you can store the hostname, username, and password or your OpenID
account in the file like this:

.. code-block:: yaml

    logon:
      hostname: "your-hostname"
      username: "your-username"
      password: "your-password"

for example

.. code-block:: yaml

    logon:
      hostname: "esgf-data.dkrz.de"
      username: "cookiemonster"
      password: "Welcome01"

if you created an account on the host `esgf-data.dkrz.de`_ with username
'cookiemonster' and password 'Welcome01'.
Alternatively, you can configure an interactive log in:

.. code-block:: yaml

    logon:
      interactive: true

Note that storing your password in plain text in the configuration
file is less secure.
On shared systems, make sure the permissions of the file are set so
only you and administrators can read it, i.e.

.. code-block:: bash

    ls -l ~/.esmvaltool/esgf-pyclient.yml

shows permissions ``-rw-------``.

Search
``````
Any arguments to :py:obj:`pyesgf.search.connection.SearchConnection` can
be provided in the section ``search_connection``, for example:

.. code-block:: yaml

    search_connection:
      url: "http://esgf-index1.ceda.ac.uk/esg-search"

to choose the CEDA index node or

.. code-block:: yaml

    search_connection:
      expire_after: 2592000  # the number of seconds in a month

to keep cached search results for a month.

The default settings are:

.. code-block:: yaml

    url: 'http://esgf-node.llnl.gov/esg-search'
    distrib: true
    timeout: 120  # seconds
    cache: '~/.esmvaltool/cache/pyesgf-search-results'
    expire_after: 86400  # cache expires after 1 day

If you experience errors while searching, it sometimes helps to delete the
cached results.

Download statistics
-------------------
The tool will maintain statistics of how fast data can be downloaded
from what host in the file ~/.esmvaltool/cache/esgf-hosts.yml and
automatically select hosts that are faster.
There is no need to manually edit this file, though it can be useful
to delete it if you move your computer to a location that is very
different from the place where you previously downloaded data.
An entry in the file might look like this:

.. code-block:: yaml

    esgf2.dkrz.de:
      duration (s): 8
      error: false
      size (bytes): 69067460
      speed (MB/s): 7.9

The tool only uses the duration and size to determine the download speed,
the speed shown in the file is not used.
If ``error`` is set to ``true``, the most recent download request to that
host failed and the tool will automatically try this host only as a last
resort.

.. _config-developer:

Developer configuration file
============================

Most users and diagnostic developers will not need to change this file,
but it may be useful to understand its content.
It will be installed along with ESMValCore and can also be viewed on GitHub:
`esmvalcore/config-developer.yml
<https://github.com/ESMValGroup/ESMValCore/blob/main/esmvalcore/config-developer.yml>`_.
This configuration file describes the file system structure and CMOR tables for several
key projects (CMIP6, CMIP5, obs4MIPs, OBS6, OBS) on several key machines (e.g. BADC, CP4CDS, DKRZ,
ETHZ, SMHI, BSC), and for native output data for some
models (ICON, IPSL, ... see :ref:`configure_native_models`).
CMIP data is stored as part of the Earth System Grid
Federation (ESGF) and the standards for file naming and paths to files are set
out by CMOR and DRS. For a detailed description of these standards and their
adoption in ESMValCore, we refer the user to :ref:`CMOR-DRS` section where we
relate these standards to the data retrieval mechanism of the ESMValCore.

By default, esmvaltool looks for it in the home directory,
inside the '.esmvaltool' folder.

Users can get a copy of this file with default values by running

.. code-block:: bash

  esmvaltool config get-config-developer --path=${TARGET_FOLDER}

If the option ``--path`` is omitted, the file will be created in
```${HOME}/.esmvaltool``.

.. note::

  Remember to change your config-user file if you want to use a custom
  config-developer.

Example of the CMIP6 project configuration:

.. code-block:: yaml

   CMIP6:
     input_dir:
       default: '/'
       BADC: '{activity}/{institute}/{dataset}/{exp}/{ensemble}/{mip}/{short_name}/{grid}/{latestversion}'
       DKRZ: '{activity}/{institute}/{dataset}/{exp}/{ensemble}/{mip}/{short_name}/{grid}/{latestversion}'
       ETHZ: '{exp}/{mip}/{short_name}/{dataset}/{ensemble}/{grid}/'
     input_file: '{short_name}_{mip}_{dataset}_{exp}_{ensemble}_{grid}*.nc'
     output_file: '{project}_{dataset}_{mip}_{exp}_{ensemble}_{short_name}'
     cmor_type: 'CMIP6'
     cmor_strict: true

Input file paths
----------------

When looking for input files, the ``esmvaltool`` command provided by
ESMValCore replaces the placeholders ``{item}`` in
``input_dir`` and ``input_file`` with the values supplied in the recipe.
ESMValCore will try to automatically fill in the values for institute, frequency,
and modeling_realm based on the information provided in the CMOR tables
and/or extra_facets_ when reading the recipe.
If this fails for some reason, these values can be provided in the recipe too.

The data directory structure of the CMIP projects is set up differently
at each site. As an example, the CMIP6 directory path on BADC would be:

.. code-block:: yaml

   '{activity}/{institute}/{dataset}/{exp}/{ensemble}/{mip}/{short_name}/{grid}/{latestversion}'

The resulting directory path would look something like this:

.. code-block:: bash

    CMIP/MOHC/HadGEM3-GC31-LL/historical/r1i1p1f3/Omon/tos/gn/latest

Please, bear in mind that ``input_dirs`` can also be a list for those  cases in
which may be needed:

.. code-block:: yaml

  - '{exp}/{ensemble}/original/{mip}/{short_name}/{grid}/{latestversion}'
  - '{exp}/{ensemble}/computed/{mip}/{short_name}/{grid}/{latestversion}'

In that case, the resultant directories will be:

.. code-block:: bash

  historical/r1i1p1f3/original/Omon/tos/gn/latest
  historical/r1i1p1f3/computed/Omon/tos/gn/latest

For a more in-depth description of how to configure ESMValCore so it can find
your data please see :ref:`CMOR-DRS`.

Preprocessor output files
-------------------------

The filename to use for preprocessed data is configured in a similar manner
using ``output_file``. Note that the extension ``.nc`` (and if applicable,
a start and end time) will automatically be appended to the filename.

.. _cmor_table_configuration:

Project CMOR table configuration
--------------------------------

ESMValCore comes bundled with several CMOR tables, which are stored in the directory
`esmvalcore/cmor/tables <https://github.com/ESMValGroup/ESMValCore/tree/main/esmvalcore/cmor/tables>`_.
These are copies of the tables available from `PCMDI <https://github.com/PCMDI>`_.

For every ``project`` that can be used in the recipe, there are four settings
related to CMOR table settings available:

* ``cmor_type``: can be ``CMIP5`` if the CMOR table is in the same format as the
  CMIP5 table or ``CMIP6`` if the table is in the same format as the CMIP6 table.
* ``cmor_strict``: if this is set to ``false``, the CMOR table will be
  extended with variables from the ``esmvalcore/cmor/tables/custom`` directory
  and it is possible to use variables with a ``mip`` which is different from
  the MIP table in which they are defined.
* ``cmor_path``: path to the CMOR table.
  Relative paths are with respect to `esmvalcore/cmor/tables`_.
  Defaults to the value provided in ``cmor_type`` written in lower case.
* ``cmor_default_table_prefix``: Prefix that needs to be added to the ``mip``
  to get the name of the file containing the ``mip`` table.
  Defaults to the value provided in ``cmor_type``.

.. _configure_native_models:

Configuring native models and observation data sets
----------------------------------------------------

ESMValCore can be configured for handling native model output formats
and specific
observation data sets without preliminary reformatting. You can choose
to host this new data source either under a dedicated project or under
project ``native6``; when choosing the latter, such a configuration
involves the following steps:

  - allowing for ESMValTool to locate the data files:

    - entry ``native6`` of ``config-developer.yml`` should be
      complemented with sub-entries for ``input_dir`` and ``input_file``
      that goes under a new key representing the
      data organization (such as ``MY_DATA_ORG``), and these sub-entries can
      use an arbitrary list of ``{placeholders}``. Example :

      .. code-block:: yaml

        native6:
          ...
          input_dir:
             default: 'Tier{tier}/{dataset}/{latestversion}/{frequency}/{short_name}'
             MY_DATA_ORG: '{model}/{exp}/{simulation}/{version}/{type}'
          input_file:
            default: '*.nc'
            MY_DATA_ORG: '{simulation}_*.nc'
          ...

    - if necessary, provide a so-called ``extra facets file`` which
      allows to cope e.g. with variable naming issues for finding
      files. See :ref:`extra_facets` and :download:`this example of
      such a file for IPSL-CM6
      <../../esmvalcore/_config/extra_facets/ipslcm-mappings.yml>`.

  - ensuring that ESMValCore get the right metadata and data out of
    your data files: this is described in :ref:`fixing_data`


.. _config-ref:

References configuration file
=============================

The `esmvaltool/config-references.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/config-references.yml>`__ file contains the list of ESMValTool diagnostic and recipe authors,
references and projects. Each author, project and reference referred to in the
documentation section of a recipe needs to be in this file in the relevant
section.

For instance, the recipe ``recipe_ocean_example.yml`` file contains the
following documentation section:

.. code-block:: yaml

  documentation:
    authors:
      - demo_le

    maintainer:
      - demo_le

    references:
      - demora2018gmd

    projects:
      - ukesm


These four items here are named people, references and projects listed in the
``config-references.yml`` file.

.. _extra_facets:

Extra Facets
============

It can be useful to automatically add extra key-value pairs to variables
or datasets in the recipe.
These key-value pairs can be used for :ref:`finding data <findingdata>`
or for providing extra information to the functions that
:ref:`fix data <extra-facets-fixes>` before passing it on to the preprocessor.

To support this, we provide the extra facets facilities. Facets are the
key-value pairs described in :ref:`Datasets`. Extra facets allows for the
addition of more details per project, dataset, mip table, and variable name.

More precisely, one can provide this information in an extra yaml file, named
`{project}-something.yml`, where `{project}` corresponds to the project as used
by ESMValTool in :ref:`Datasets` and "something" is arbitrary.

Format of the extra facets files
--------------------------------
The extra facets are given in a yaml file, whose file name identifies the
project. Inside the file there is a hierarchy of nested dictionaries with the
following levels. At the top there is the `dataset` facet, followed by the `mip`
table, and finally the `short_name`. The leaf dictionary placed here gives the
extra facets that will be made available to data finder and the fix
infrastructure. The following example illustrates the concept.

.. _extra-facets-example-1:

.. code-block:: yaml
   :caption: Extra facet example file `native6-era5.yml`

   ERA5:
     Amon:
       tas: {source_var_name: "t2m", cds_var_name: "2m_temperature"}

The three levels of keys in this mapping can contain
`Unix shell-style wildcards <https://en.wikipedia.org/wiki/Glob_(programming)#Syntax>`_.
The special characters used in shell-style wildcards are:

+------------+----------------------------------------+
|Pattern     | Meaning                                |
+============+========================================+
| ``*``      |   matches everything                   |
+------------+----------------------------------------+
| ``?``      |   matches any single character         |
+------------+----------------------------------------+
| ``[seq]``  |   matches any character in ``seq``     |
+------------+----------------------------------------+
| ``[!seq]`` |   matches any character not in ``seq`` |
+------------+----------------------------------------+

where ``seq`` can either be a sequence of characters or just a bunch of characters,
for example ``[A-C]`` matches the characters ``A``, ``B``, and ``C``,
while ``[AC]`` matches the characters ``A`` and ``C``.

For example, this is used to automatically add ``product: output1`` to any
variable of any CMIP5 dataset that does not have a ``product`` key yet:

.. code-block:: yaml
   :caption: Extra facet example file `cmip5-product.yml <https://github.com/ESMValGroup/ESMValCore/blob/main/esmvalcore/_config/extra_facets/cmip5-product.yml>`_

   '*':
     '*':
       '*': {product: output1}

Location of the extra facets files
----------------------------------
Extra facets files can be placed in several different places. When we use them
to support a particular use-case within the ESMValTool project, they will be
provided in the sub-folder `extra_facets` inside the package
`esmvalcore._config`. If they are used from the user side, they can be either
placed in `~/.esmvaltool/extra_facets` or in any other directory of the users
choosing. In that case this directory must be added to the `config-user.yml`
file under the `extra_facets_dir` setting, which can take a single directory or
a list of directories.

The order in which the directories are searched is

1. The internal directory `esmvalcore._config/extra_facets`
2. The default user directory `~/.esmvaltool/extra_facets`
3. The custom user directories in the order in which they are given in
   `config-user.yml`.

The extra facets files within each of these directories are processed in
lexicographical order according to their file name.

In all cases it is allowed to supersede information from earlier files in later
files. This makes it possible for the user to effectively override even internal
default facets, for example to deal with local particularities in the data
handling.

Use of extra facets
-------------------
For extra facets to be useful, the information that they provide must be
applied. There are fundamentally two places where this comes into play. One is
:ref:`the datafinder<extra-facets-data-finder>`, the other are
:ref:`fixes<extra-facets-fixes>`.
.. _install:

Installation
============

Conda installation
------------------

In order to install the Conda package, you will need to install `Conda <https://docs.conda.io>`_ first.
For a minimal conda installation (recommended) go to https://conda.io/miniconda.html.
It is recommended that you always use the latest version of conda, as problems have been reported when trying to use older versions.

Once you have installed conda, you can install ESMValCore by running:

.. code-block:: bash

    conda install -c conda-forge esmvalcore

It is also possible to create a new
`Conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_
and install ESMValCore into it with a single command:

.. code-block:: bash

    conda create --name esmvalcore -c conda-forge esmvalcore 'python=3.10'

Don't forget to activate the newly created environment after the installation:

.. code-block:: bash

    conda activate esmvalcore

Of course it is also possible to choose a different name than ``esmvalcore`` for the environment.

.. note::

	  Creating a new Conda environment is often much faster and more reliable than trying to update an existing Conda environment.

Pip installation
-----------------

It is also possible to install ESMValCore from `PyPI <https://pypi.org/project/ESMValCore/>`_.
However, this requires first installing dependencies that are not available on PyPI in some other way.
By far the easiest way to install these dependencies is to use conda_.
For a minimal conda installation (recommended) go to https://conda.io/miniconda.html.

After installing Conda, download
`the file with the list of dependencies <https://raw.githubusercontent.com/ESMValGroup/ESMValCore/main/environment.yml>`_:

.. code-block:: bash

    wget https://raw.githubusercontent.com/ESMValGroup/ESMValCore/main/environment.yml

and install these dependencies into a new conda environment with the command

.. code-block:: bash

    conda env create --name esmvalcore --file environment.yml

Finally, activate the newly created environment

.. code-block:: bash

    conda activate esmvalcore

and install ESMValCore as well as any remaining dependencies with the command:

.. code-block:: bash

    pip install esmvalcore


Docker installation
-----------------------

ESMValCore is also provided through `DockerHub <https://hub.docker.com/u/esmvalgroup/>`_
in the form of docker containers.
See https://docs.docker.com for more information about docker containers and how to
run them.

You can get the latest release with

.. code-block:: bash

   docker pull esmvalgroup/esmvalcore:stable

If you want to use the current main branch, use

.. code-block:: bash

   docker pull esmvalgroup/esmvalcore:latest

To run a container using those images, use:

.. code-block:: bash

   docker run esmvalgroup/esmvalcore:stable --help

Note that the container does not see the data or environmental variables available in the host by default.
You can make data available with ``-v /path:/path/in/container`` and environmental variables with ``-e VARNAME``.

For example, the following command would run a recipe

.. code-block:: bash

   docker run -e HOME -v "$HOME":"$HOME" -v /data:/data esmvalgroup/esmvalcore:stable -c ~/config-user.yml ~/recipes/recipe_example.yml

with the environmental variable ``$HOME`` available inside the container and the data
in the directories ``$HOME`` and ``/data``, so these can be used to find the configuration file, recipe, and data.

It might be useful to define a `bash alias
<https://opensource.com/article/19/7/bash-aliases>`_
or script to abbreviate the above command, for example

.. code-block:: bash

	 alias esmvaltool="docker run -e HOME -v $HOME:$HOME -v /data:/data esmvalgroup/esmvalcore:stable"

would allow using the ``esmvaltool`` command without even noticing that the tool is running inside a Docker container.


Singularity installation
----------------------------

Docker is usually forbidden in clusters due to security reasons. However,
there is a more secure alternative to run containers that is usually available
on them: `Singularity <https://sylabs.io/guides/3.0/user-guide/quick_start.html>`_.

Singularity can use docker containers directly from DockerHub with the
following command

.. code-block:: bash

   singularity run docker://esmvalgroup/esmvalcore:stable -c ~/config-user.yml ~/recipes/recipe_example.yml

Note that the container does not see the data available in the host by default.
You can make host data available with ``-B /path:/path/in/container``.

It might be useful to define a `bash alias
<https://opensource.com/article/19/7/bash-aliases>`_
or script to abbreviate the above command, for example

.. code-block:: bash

	 alias esmvaltool="singularity run -B $HOME:$HOME -B /data:/data docker://esmvalgroup/esmvalcore:stable"

would allow using the ``esmvaltool`` command without even noticing that the tool is running inside a Singularity container.

Some clusters may not allow to connect to external services, in those cases
you can first create a singularity image locally:

.. code-block:: bash

   singularity build esmvalcore.sif docker://esmvalgroup/esmvalcore:stable

and then upload the image file ``esmvalcore.sif`` to the cluster.
To run the container using the image file ``esmvalcore.sif`` use:

.. code-block:: bash

   singularity run esmvalcore.sif -c ~/config-user.yml ~/recipes/recipe_example.yml

.. _installation-from-source:

Installation from source
------------------------

.. note::
    If you would like to install the development version of ESMValCore alongside
    ESMValTool, please have a look at
    :ref:`these instructions <esmvaltool:esmvalcore-development-installation>`.

To install from source for development, follow these instructions.

-  `Download and install
   conda <https://conda.io/projects/conda/en/latest/user-guide/install/linux.html>`__
   (this should be done even if the system in use already has a
   preinstalled version of conda, as problems have been reported with
   using older versions of conda)
-  To make the ``conda`` command available, add
   ``source <prefix>/etc/profile.d/conda.sh`` to your ``.bashrc`` file
   and restart your shell. If using (t)csh shell, add
   ``source <prefix>/etc/profile.d/conda.csh`` to your
   ``.cshrc``/``.tcshrc`` file instead.
-  Update conda: ``conda update -y conda``
-  Clone the ESMValCore Git repository:
   ``git clone https://github.com/ESMValGroup/ESMValCore.git``
-  Go to the source code directory: ``cd ESMValCore``
-  Create the esmvalcore conda environment
   ``conda env create --name esmvalcore --file environment.yml``
-  Activate the esmvalcore environment: ``conda activate esmvalcore``
-  Install in development mode: ``pip install -e '.[develop]'``. If you
   are installing behind a proxy that does not trust the usual pip-urls
   you can declare them with the option ``--trusted-host``,
   e.g. ``pip install --trusted-host=pypi.python.org --trusted-host=pypi.org --trusted-host=files.pythonhosted.org -e .[develop]``
-  Test that your installation was successful by running
   ``esmvaltool -h``.

Pre-installed versions on HPC clusters
--------------------------------------

You will find the tool available on HPC clusters and there will be no need to install it
yourself if you are just running diagnostics:

 - CEDA-JASMIN: `esmvaltool` is available on the scientific compute nodes (`sciX.jasmin.ac.uk` where
   `X = 1, 2, 3, 4, 5`) after login and module loading via `module load esmvaltool`; see the helper page at
   `CEDA <https://help.jasmin.ac.uk/article/4955-community-software-esmvaltool>`__ ;
 - DKRZ-Mistral: `esmvaltool` is available on login nodes (`mistral.dkrz.de`) and pre- and post-processing
   nodes (`mistralpp.dkrz.de`) after login and module loading via `module load esmvaltool`; the command
   `module help esmvaltool` provides some information about the module.

.. note::
    If you would like to use pre-installed versions on HPC clusters (currently CEDA-JASMIN and DKRZ-MISTRAL),
    please have a look at
    :ref:`these instructions <esmvaltool:install_on_hpc>`.


Installation from the conda lock file
-------------------------------------

A fast conda environment creation is possible using the provided conda lock files. This is a secure alternative
to the installation from source, whenever the conda environment can not be created for some reason. A conda lock file
is an explicit environment file that contains pointers to dependency packages as they are hosted on the Anaconda cloud;
these have frozen version numbers, build hashes, and channel names, parameters established at the time
of the conda lock file creation, so may be obsolete after a while,
but they allow for a robust environment creation while they're still up-to-date.
We regenerate these lock files every 10 days, so to minimize the risk of dependencies becoming obsolete.
Conda environment creation from a lock file is done just like with any other environment file:

.. code-block:: bash

   conda create --name esmvaltool --file conda-linux-64.lock


.. note::
   `pip` and `conda` are NOT installed, so you will have to install them in the new environment: use conda-forge as channel): ``conda install -c conda-forge pip`` at the very minimum so we can install `esmvalcore` afterwards.


Creating a conda lock file
--------------------------

We provide a conda lock file for Linux-based operating systems, but if you prefer to
build a conda lock file yourself, install the `conda-lock` package first:

.. code-block:: bash

   conda install -c conda-forge conda-lock

then run

.. code-block:: bash

   conda-lock lock --platform linux-64 -f environment.yml --mamba

(mamba activated for speed) to create a conda lock file for Linux platforms,
or run

.. code-block:: bash

   conda-lock lock --platform osx-64 -f environment.yml --mamba

to create a lock file for OSX platforms. Note however, that using conda lock files on OSX is still problematic!
