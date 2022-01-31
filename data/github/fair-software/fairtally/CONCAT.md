# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 09-March-2021

First release of fairtally.

[Unreleased]: https://github.com/fair-software/fairtally/compare/0.1.0...HEAD
[0.1.0]: https://github.com/fair-software/fairtally/releases/tag/0.1.0
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

#. use the search functionality `here <https://github.com/fair-software/fairtally/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/fair-software/fairtally/issues>`__ to see if someone already filed the same issue;
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
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the fairtally repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
``fairtally`` developer documentation
=====================================

If you're looking for user documentation, go `here <README.rst>`_.


The project setup is documented in `a separate document <project_setup.rst>`_. Feel free to remove this document (and/or the link to this document) if you don't need it.



Development install
-------------------

.. code:: shell

    # Create a virtualenv, e.g. with
    python3 -m venv venv3

    # activate virtualenv
    source venv3/bin/activate

    # make sure to have a recent version of pip
    pip install --upgrade pip

    # (from the project root directory)
    # install fairtally as an editable package
    pip install --no-cache-dir --editable .
    # install development dependencies
    pip install --no-cache-dir --editable .[dev]

Afterwards check that the install directory was added to the ``PATH``
environment variable. You should then be able to call the executable,
like so:

.. code:: shell

    fairtally --help


Running the tests
-----------------

Running the tests requires an activated virtualenv with the development tools installed.

.. code:: shell

    # unit tests with mocked representations of repository behavior
    pytest

The test coverage report is available in `htmlcov/index.html`.
Running linters locally
-----------------------

Running the linters requires an activated virtualenv with the development tools installed.

.. code:: shell

    # linter
    prospector

    # recursively check import style for the fairtally module only
    isort --recursive --check-only fairtally

    # recursively check import style for the fairtally module only and show
    # any proposed changes as a diff
    isort --recursive --check-only --diff fairtally

    # recursively fix import style for the fairtally module only
    isort --recursive fairtally


You can enable automatic linting with ``prospector`` and ``isort`` on commit like so:

.. code:: shell

    git config --local core.hooksPath .githooks

Versioning
----------

Bumping the version across all files is done with bump2version, e.g.

.. code:: shell

    bump2version minor


Example report
--------------

The `example HTML report <https://fair-software.github.io/fairtally/_static/fairtally_example.html>`_ and its `screenshot <https://github.com/fair-software/fairtally/raw/main/docs/_static/fairtally_example.png>`_ need to be updated when the User Interface or example command changes.

Update example HTML report by running

.. code:: shell

  fairtally -o docs/_static/fairtally_example.html https://github.com/fair-software/fairtally https://github.com/fair-software/howfairis

Update screenshot (``docs/_static/fairtally_example.png``), for example by running Google Chrome screenshot command with:

.. code:: shell

  google-chrome --headless --disable-gpu --window-size=1150,280 --screenshot=docs/_static/fairtally_example.png docs/_static/fairtally_example.html

If size of report changed then adjust the ``--window-size`` argument value accordingly.

Making a release
----------------

Preparation
^^^^^^^^^^^

1. Update the ``CHANGELOG.md``
2. Verify that the information in ``CITATION.cff`` is correct, and that ``.zenodo.json`` contains equivalent data
3. Make sure the version has been updated.
4. Run the unit tests with ``pytest``
5. Make sure `example report <#example-report>`_ is up to date

PyPI
^^^^

In a new terminal, without an activated virtual environment or a venv3 directory:

.. code:: shell

    # prepare a new directory
    cd $(mktemp -d --tmpdir fairtally.XXXXXX)

    # fresh git clone ensures the release has the state of origin/main branch
    git clone https://github.com/fair-software/fairtally.git .

    # prepare a clean virtual environment and activate it
    python3 -m venv venv3
    source venv3/bin/activate

    # make sure to have a recent version of pip
    pip install --upgrade pip

    # install runtime dependencies and publishing dependencies
    pip install --no-cache-dir .
    pip install --no-cache-dir .[publishing]

    # clean up any previously generated artefacts
    rm -rf fairtally.egg-info
    rm -rf dist

    # create the source distribution and the wheel
    python setup.py sdist bdist_wheel

    # upload to test pypi instance (requires credentials)
    twine upload --repository-url https://test.pypi.org/legacy/ dist/*

In a new terminal, without an activated virtual environment or a venv3 directory:

.. code:: shell

    cd $(mktemp -d --tmpdir fairtally-test.XXXXXX)

    # check you don't have an existing fairtally
    which fairtally
    python3 -m pip uninstall fairtally

    # install in user space from test pypi instance:
    python3 -m pip -v install --user --no-cache-dir \
    --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple fairtally

Check that the package works as it should when installed from pypitest.

Then upload to pypi.org with:

.. code:: shell

    # Back to the first terminal,
    # FINAL STEP: upload to PyPI (requires credentials)
    twine upload dist/*

GitHub
^^^^^^

Don't forget to also `make a release on GitHub <https://github.com/fair-software/fairtally/releases/new>`_.

DockerHub
^^^^^^^^^

To build the image, run:

.. code:: shell

    docker build -t fairsoftware/fairtally:latest .

.. code:: shell

    VERSION=<your-version>
    docker tag fairsoftware/fairtally:latest fairsoftware/fairtally:${VERSION}

Check that you have the tags you want with:

.. code:: shell

    docker images

To push the image to DockerHub, run:

.. code:: shell

    # (requires credentials)
    docker login
    docker push fairsoftware/fairtally:${VERSION}
    docker push fairsoftware/fairtally:latest

The new image and its tags should now be listed here https://hub.docker.com/r/fairsoftware/fairtally/tags?page=1&ordering=last_updated.
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
reported by contacting the project team at j.spaaks@esciencecenter.nl. All
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

  cd fairtally
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

If you decide to use pypi for distributing you code, you can configure travis to upload to pypi when you make a release. If you specified your pypi user name during generation of this package, the ``.travis.yml`` file contains a section that looks like:

.. code-block:: yaml

  deploy:
    provider: pypi
    user: howfairis
    password:
      secure: FIXME; see README for more info
   on:
      tags: true
      branch: master

Before this actually works, you need to add an encrypted password for your pypi account. The `travis documentation <https://docs.travis-ci.com/user/deployment/pypi/>`_ specifies how to do this.

Testing and code coverage
-------------------------

* Tests should be put in the ``tests`` folder.
* The ``tests`` folder contains:

  - Example tests that you should replace with your own meaningful tests (file: ``test_fairtally``)

* The testing framework used is `PyTest <https://pytest.org>`_

  - `PyTest introduction <http://pythontesting.net/framework/pytest/pytest-introduction/>`_

* Tests can be run with ``pytest``

  - This is configured in ``setup.py`` and ``setup.cfg``

* Use `Travis CI <https://travis-ci.com/>`_ to automatically run tests and to test using multiple Python versions

  - Configuration can be found in ``.travis.yml``
  - `Getting started with Travis CI <https://docs.travis-ci.com/user/getting-started/>`_

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

  - Include the link to the documentation in this README_.

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
* For convenience, the package version is stored in a single place: ``fairtally/__version__.py``. For updating the version number, you only have to change this file.
* Don't forget to update the version number before `making a release <https://guide.esciencecenter.nl/best_practices/releases.html>`_!


Logging
-------

* We recommend using the `logging` module for getting useful information from your module (instead of using `print`).
* The project is set up with a logging example.
* `Relevant section in the guide <https://guide.esciencecenter.nl/best_practices/language_guides/python.html#logging>`_

CHANGELOG.md
------------

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
################################################################################
fairtally
################################################################################

Python application to analyze multiple GitHub and GitLab repositories compliance with the `fair-software.eu <fair-software.eu>`_ recommendations.

.. list-table::
   :widths: 25 25
   :header-rows: 1

   * - fair-software.nl recommendations
     - Badges
   * - \1. Code repository
     - |GitHub Badge|
   * - \2. License
     - |License Badge|
   * - \3. Community Registry
     - |PyPI Badge| |Research Software Directory Badge|
   * - \4. Enable Citation
     - |Zenodo Badge|
   * - \5. Checklist
     - |CII Best Practices Badge|
   * - **Other best practices**
     -
   * - Continuous integration
     - |Python Build| |Linter|
   * - DockerHub
     - |dockerhub badge|

.. |GitHub Badge| image:: https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue
   :target: https://github.com/fair-software/fairtally
   :alt: GitHub Badge

.. |License Badge| image:: https://img.shields.io/github/license/fair-software/fairtally
   :target: https://github.com/fair-software/fairtally
   :alt: License Badge

.. |PyPI Badge| image:: https://img.shields.io/pypi/v/fairtally.svg?colorB=blue
   :target: https://pypi.python.org/project/fairtally/
   :alt: PyPI Badge
.. |Research Software Directory Badge| image:: https://img.shields.io/badge/rsd-fairtally-00a3e3.svg
   :target: https://www.research-software.nl/software/fairtally
   :alt: Research Software Directory Badge

..
    Goto https://zenodo.org/account/settings/github/ to enable Zenodo/GitHub integration.
    After creation of a GitHub release at https://github.com/fair-software/fairtally/releases
    there will be a Zenodo upload created at https://zenodo.org/deposit with a DOI, this DOI can be put in the Zenodo badge urls.
    In the README, we prefer to use the concept DOI over versioned DOI, see https://help.zenodo.org/#versioning.
.. |Zenodo Badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4590882.svg
   :target: https://doi.org/10.5281/zenodo.4590882
   :alt: Zenodo Badge

..
    A CII Best Practices project can be created at https://bestpractices.coreinfrastructure.org/en/projects/new
.. |CII Best Practices Badge| image:: https://bestpractices.coreinfrastructure.org/projects/4690/badge
   :target: https://bestpractices.coreinfrastructure.org/en/projects/4690
   :alt: CII Best Practices Badge

.. |Python Build| image:: https://github.com/fair-software/fairtally/actions/workflows/build.yml/badge.svg
   :target: https://github.com/fair-software/fairtally/actions?query=workflow%3A%22build%22
   :alt: Python Build

.. |Linter| image:: https://github.com/fair-software/fairtally/actions/workflows/linting.yml/badge.svg
   :target: https://github.com/fair-software/fairtally/actions?query=workflow%3A%22Linting%22
   :alt: Linter

.. |dockerhub badge| image:: https://img.shields.io/docker/pulls/fairsoftware/fairtally
   :target: https://hub.docker.com/r/fairsoftware/fairtally
   :alt: Docker Pulls

Installation
------------

To install fairtally, do:

.. code-block:: console

  pip3 install --user fairtally

Usage
-----

To use fairtally to check the compliance of multiple repositories, one can use the command below.

.. code-block:: console

  fairtally https://github.com/fair-software/fairtally https://github.com/fair-software/howfairis

This command will generate a html report called `tally.html` which will contain the results of the checks for each repository.

Then open the analysis in a web-browser, for example Firefox:

.. code-block:: console

  firefox tally.html

The report will look similar to the example below:

.. image:: docs/_static/fairtally_example.png
  :target: https://fair-software.github.io/fairtally/_static/fairtally_example.html

You can sort the table by clicking on the table headers. The purple plus signs provide access to log messages of each repository.

  Checking many repositories will quickly exceed the rate limit of the APIs of GitLab and GitHub and resulting in all remaining repositories to be fully non-compliantly.
  See `howfairis docs <https://github.com/fair-software/howfairis/#rate-limit>`_ how setup environment variables to increase the rate limit.

Using Docker image
------------------

You can run fairtally Docker image using the command below.

.. code:: console

    docker pull fairsoftware/fairtally

You can run fairtally Docker image using the command below.

.. code:: console

    docker run --rm fairsoftware/fairtally --help

`--rm` argument will remove Docker container after execution.

To tally 2 URLs and save the report as `tally.html` in the current working directory you can run the command below.

.. code:: console

    docker run --rm fairsoftware/fairtally -o - https://github.com/fair-software/fairtally https://github.com/fair-software/howfairis > tally.html

See developer documentation to learn how to modify the Docker image.

Research Software Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To FAIR tally the software listed on the `Research Software Directory of the Netherlands eScience Center <https://research-software.nl/>`_.

First download a list of software by calling `RSD API <https://github.com/research-software-directory/research-software-directory/blob/master/docs/documentation-for-developers.md#api>`_

.. code-block:: console

  curl https://research-software.nl/api/software > software.json

Next, extract the repository URLs with `jq <https://stedolan.github.io/jq/>`_.

.. code-block:: console

  cat software.json | jq -r '[.[].repositoryURLs.github] | flatten | .[]' > urls.txt

Finally run fairtally to generate a report.

.. code-block:: console

  fairtally --output-file report.html --input-file urls.txt

Documentation
*************

Command line interface help can be retrieved with

.. code-block:: console

  fairtally --help

The output of the command will be something like:

.. code-block:: console

  Usage: fairtally [OPTIONS] [URLS]...

  Options:
    -o, --output-file TEXT     Filename of where to write the results. Use `-`
                               to write to standard out.  [default: tally.html]

    -i, --input-file FILENAME  Check URLs in file. One URL per line. Use `-` to
                               read from standard input.

    --format [html|json]       Format of output.  [default: html]
    --version                  Show the version and exit.
    --help                     Show this message and exit.

Contributing
************

If you want to contribute to the development of fairtally,
have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

License
*******

Copyright (c) 2021, Netherlands eScience Center

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Credits
*******

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.

Instructions for developers
***************************

The developer documentation can be found in `README.dev.rst <README.dev.rst>`_.
.. fairtally documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to fairtally's documentation!
==========================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
