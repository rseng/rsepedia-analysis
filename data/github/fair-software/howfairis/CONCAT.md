# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.14.1] - 2021-Mar-09

### Added

- Describe how to get api keys [#319](https://github.com/fair-software/howfairis/issues/319)
- Dont show cant-remove-comment warning for rst file without comments [#272](https://github.com/fair-software/howfairis/issues/272)

## [0.14.0] - 2021-Mar-02

### Added

- Can now ask instances of `Compliance` for their color [#301](https://github.com/fair-software/howfairis/issues/301)
- Can now ask instances of `Compliance` for their badge image url [#304](https://github.com/fair-software/howfairis/issues/304)
- Now optionally uses authenticated requests when making requests to github.com and gitlab.com 
- Rate limits are now configurable and use exponential backoff and retry (adds [ratelimit](https://pypi.org/project/ratelimit/) and [backoff](https://pypi.org/project/backoff/) dependencies) [PR#286](https://github.com/fair-software/howfairis/pull/286)
- Now warns about GitHub's caching when using READMEs that were recently changed [PR#153](https://github.com/fair-software/howfairis/pull/153)
- Directory structure of tests was updated for conceptually more meaningful scenarios, improved consistency between platforms, and directory-level mocked API calls using `pytest`'s standard `conftest.py` pattern. [PR#285](https://github.com/fair-software/howfairis/pull/285)
- More tests, e.g. [PR#305](https://github.com/fair-software/howfairis/pull/305), [PR#293](https://github.com/fair-software/howfairis/pull/293)

## [0.13.0] - 2021-Feb-18

### Added
- docstrings for public API
- documentation hosted on readthedocs [#51](https://github.com/fair-software/howfairis/issues/51)
- code of conduct [#87](https://github.com/fair-software/howfairis/issues/87)
- contributing guide [#74](https://github.com/fair-software/howfairis/issues/74)
- Docker image [#62](https://github.com/fair-software/howfairis/issues/62)
- developer documentation [#83](https://github.com/fair-software/howfairis/issues/83)
- adhere to fair-software recommendations [#50](https://github.com/fair-software/howfairis/issues/50) [#53](https://github.com/fair-software/howfairis/issues/53) [#137](https://github.com/fair-software/howfairis/issues/137) [#151](https://github.com/fair-software/howfairis/pull/151)
- support for more anaconda badges [#124](https://github.com/fair-software/howfairis/issues/124)
- warning for commented badges in README.rst [#72](https://github.com/fair-software/howfairis/issues/72)
- quiet mode to the cli [#182](https://github.com/fair-software/howfairis/issues/182)
- Readme.get_compliance() [#94](https://github.com/fair-software/howfairis/issues/94)
- retrieve the default branch [#48](https://github.com/fair-software/howfairis/issues/48)

### Changed
- automated tests for Python 3.6, 3.7, 3.8, 3.9 [#80](https://github.com/fair-software/howfairis/issues/80)
- rename configuration keys [#164](https://github.com/fair-software/howfairis/issues/164) [#179](https://github.com/fair-software/howfairis/issues/179)
- users can now add a reason if they want to skip a check [#179](https://github.com/fair-software/howfairis/issues/179)
- Config class is merged into the Checker class [#172](https://github.com/fair-software/howfairis/issues/172) [#194](https://github.com/fair-software/howfairis/issues/194)
- Checker.check_five_recommendations() now returns Compliance object [#145](https://github.com/fair-software/howfairis/issues/145)
- moved badge generation to Compliance class [#94](https://github.com/fair-software/howfairis/issues/94)
- renamed config related argument names of cli [#172](https://github.com/fair-software/howfairis/issues/172) [#194](https://github.com/fair-software/howfairis/issues/194)

### Removed
- option to set compliant symbol [#178](https://github.com/fair-software/howfairis/issues/178)
- config argument from the Repo constructor [#194](https://github.com/fair-software/howfairis/issues/194)

## [0.12.0] - 2020-December-09

We started to keep a changelog after this release.


[Unreleased]: https://github.com/fair-software/howfairis/compare/0.14.1...HEAD
[0.14.1]: https://github.com/fair-software/howfairis/compare/0.14.0...0.14.1
[0.14.0]: https://github.com/fair-software/howfairis/compare/0.13.0...0.14.0
[0.13.0]: https://github.com/fair-software/howfairis/compare/0.12.0...0.13.0
[0.12.0]: https://github.com/fair-software/howfairis/releases/tag/0.12.0
**List of related issues or pull requests**

Refs: #ISSUE_NUMBER

**Describe the changes made in this pull request**

**Instructions to review the pull request**
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

#. use the search functionality `here <https://github.com/fair-software/howfairis/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/fair-software/howfairis/issues>`__ to see if someone already filed the same issue;
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
#. if cloned locally then perform a `development installation <README.dev.rst#development-install>`_;
#. make sure the existing tests still work by running ``pytest``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the howfairis repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
``howfairis`` developer documentation
=====================================

If you're looking for user documentation, go `here <README.rst>`_.

|
|

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
    # install howfairis as an editable package
    pip install --no-cache-dir --editable .
    # install development dependencies
    pip install --no-cache-dir --editable .[dev]

Afterwards check that the install directory was added to the ``PATH``
environment variable. You should then be able to call the executable,
like so:

.. code:: shell

    howfairis https://github.com/<owner>/<repo>

Running the tests
-----------------

Running the tests requires an activated virtualenv with the development tools installed.

.. code:: shell

    # unit tests with mocked representations of repository behavior
    pytest
    pytest tests/
    
    # live tests with actual repository behavior (slow, prone to HttpError too many requests)
    pytest livetests/
    
    # command line interface tests
    bash clitests/script.sh

Running linters locally
-----------------------

Running the linters requires an activated virtualenv with the development tools installed.

.. code:: shell

    # linter
    prospector

    # recursively check import style for the howfairis module only
    isort --recursive --check-only howfairis

    # recursively check import style for the howfairis module only and show
    # any proposed changes as a diff
    isort --recursive --check-only --diff howfairis

    # recursively fix  import style for the howfairis module only
    isort --recursive howfairis

.. code:: shell

    # requires activated virtualenv with development tools
    prospector && isort --recursive --check-only howfairis

You can enable automatic linting with ``prospector`` and ``isort`` on commit like so:

.. code:: shell

    git config --local core.hooksPath .githooks

Versioning
----------

Bumping the version across all files is done with bump2version, e.g.

.. code:: shell

    bump2version minor


Making a release
----------------

Preparation
^^^^^^^^^^^

1. Update the ``CHANGELOG.md``
2. Verify that the information in ``CITATION.cff`` is correct, and that ``.zenodo.json`` contains equivalent data
3. Make sure the version has been updated.
4. Run the unit tests with ``pytest tests/``
5. Run the live tests with ``pytest livetests/``
6. Run the clitests with ``bash clitests/script.sh``

PyPI
^^^^

In a new terminal, without an activated virtual environment or a venv3 directory:

.. code:: shell

    # prepare a new directory
    cd $(mktemp -d --tmpdir howfairis.XXXXXX)
    
    # fresh git clone ensures the release has the state of origin/main branch
    git clone https://github.com/fair-software/howfairis.git .
    
    # prepare a clean virtual environment and activate it
    python3 -m venv venv3
    source venv3/bin/activate
    
    # make sure to have a recent version of pip
    pip install --upgrade pip 

    # install runtime dependencies and publishing dependencies
    pip install --no-cache-dir .
    pip install --no-cache-dir .[publishing]
    
    # clean up any previously generated artefacts 
    rm -rf howfairis.egg-info
    rm -rf dist
    
    # create the source distribution and the wheel
    python setup.py sdist bdist_wheel

    # upload to test pypi instance (requires credentials)
    twine upload --repository-url https://test.pypi.org/legacy/ dist/*

In a new terminal, without an activated virtual environment or a venv3 directory:

.. code:: shell
    
    cd $(mktemp -d --tmpdir howfairis-test.XXXXXX)

    # check you don't have an existing howfairis
    which howfairis
    python3 -m pip uninstall howfairis

    # install in user space from test pypi instance:
    python3 -m pip -v install --user --no-cache-dir \
    --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple howfairis

Check that the package works as it should when installed from pypitest.

Then upload to pypi.org with:

.. code:: shell

    # Back to the first terminal,
    # FINAL STEP: upload to PyPI (requires credentials)
    twine upload dist/*

GitHub
^^^^^^

Don't forget to also make a release on GitHub.

DockerHub
^^^^^^^^^

To build the image, run:

.. code:: shell

    docker build -t fairsoftware/howfairis:latest .
    
.. code:: shell

    VERSION=<your-version>
    docker tag fairsoftware/howfairis:latest fairsoftware/howfairis:${VERSION}

Check that you have the tags you want with:

.. code:: shell

    docker images

To push the image to DockerHub, run:

.. code:: shell

    # (requires credentials)  
    docker login
    docker push fairsoftware/howfairis:${VERSION}
    docker push fairsoftware/howfairis:latest    
    
The new image and its tags should now be listed here https://hub.docker.com/r/fairsoftware/howfairis/tags?page=1&ordering=last_updated.
    
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
reported by contacting the project team at gfairsoftware@esciencecenter.nl. All
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
howfairis
=========

|

Python package to analyze a GitHub or GitLab repository's compliance with the
fair-software.eu_ recommendations.

Badges
------


====================================================== ============================
fair-software.nl recommendations
====================================================== ============================
(1/5) code repository                                  |github repo badge|
(2/5) license                                          |github license badge|
(3/5) community registry                               |pypi badge|
(4/5) citation                                         |zenodo badge|
(5/5) checklist                                        |core infrastructures badge|
overall                                                |fair-software badge|
**Other best practices**
Documentation                                          |readthedocs badge|
Supported Python versions                              |python versions badge| 
Code quality                                           |sonarcloud quality badge|
Code coverage of unit tests                            |sonarcloud coverage badge|
DockerHub                                              |dockerhub badge|
**GitHub Actions**
Citation metadata consistency                          |workflow cffconvert badge|
Unit tests                                             |workflow tests badge|
Live tests (triggered manually)                        |workflow livetests badge|
====================================================== ============================

.. |github repo badge| image:: https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue
   :target: https://github.com/fair-software/howfairis

.. |github license badge| image:: https://img.shields.io/github/license/fair-software/howfairis
   :target: https://github.com/fair-software/howfairis

.. |pypi badge| image:: https://img.shields.io/pypi/v/howfairis.svg?colorB=blue
   :target: https://pypi.python.org/pypi/howfairis/

.. |zenodo badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4017908.svg
   :target: https://doi.org/10.5281/zenodo.4017908
   
.. |core infrastructures badge| image:: https://bestpractices.coreinfrastructure.org/projects/4630/badge
   :target: https://bestpractices.coreinfrastructure.org/en/projects/4630

.. |fair-software badge| image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green
   :target: https://fair-software.eu
   
.. |readthedocs badge| image:: https://readthedocs.org/projects/howfairis/badge/?version=latest
   :target: https://howfairis.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
   
.. |python versions badge| image:: https://img.shields.io/pypi/pyversions/howfairis.svg
   :target: https://pypi.python.org/pypi/howfairis   

.. |sonarcloud quality badge| image:: https://sonarcloud.io/api/project_badges/measure?project=fair-software_howfairis&metric=alert_status
   :target: https://sonarcloud.io/dashboard?id=fair-software_howfairis
   :alt: Quality Gate Status

.. |sonarcloud coverage badge| image:: https://sonarcloud.io/api/project_badges/measure?project=fair-software_howfairis&metric=coverage
   :target: https://sonarcloud.io/dashboard?id=fair-software_howfairis
   :alt: Coverage

.. |dockerhub badge| image:: https://img.shields.io/docker/pulls/fairsoftware/howfairis
   :target: https://hub.docker.com/r/fairsoftware/howfairis
   :alt: Docker Pulls

.. |workflow tests badge| image:: https://github.com/fair-software/howfairis/workflows/tests/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3Atests

.. |workflow livetests badge| image:: https://github.com/fair-software/howfairis/workflows/livetests/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3Alivetests

.. |workflow cffconvert badge| image:: https://github.com/fair-software/howfairis/workflows/metadata%20consistency/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3A%22metadata+consistency%22

Install
-------

.. code:: console

    pip3 install --user howfairis

Verify that the install directory is on the ``PATH`` environment variable. If so,
you should be able to call the executable, like so:

.. code:: console

    howfairis https://github.com/<owner>/<repo>


``howfairis`` supports URLs from the following code repository platforms:

1. ``https://github.com``
2. ``https://gitlab.com`` (not including self-hosted instances)

Docker
---------------

You can run howfairis Docker image using the command below.

.. code:: console

    docker pull fairsoftware/howfairis

You can run howfairis Docker image using the command below.

.. code:: console

    docker run --rm fairsoftware/howfairis --help

`--rm` argument will remove Docker container after execution.

See developer documentation to learn how to modify the Docker image.

Expected output
---------------

Depending on which repository you are doing the analysis for, the output
looks something like this:

.. code:: console

    Checking compliance with fair-software.eu...
    url: https://github.com/fair-software/badge-test
    (1/5) repository
          ✓ has_open_repository
    (2/5) license
          ✓ has_license
    (3/5) registry
          × has_ascl_badge
          × has_bintray_badge
          × has_conda_badge
          × has_cran_badge
          × has_crates_badge
          × has_maven_badge
          × has_npm_badge
          ✓ has_pypi_badge
          × has_rsd_badge
          × is_on_github_marketplace
    (4/5) citation
          × has_citation_file
          × has_citationcff_file
          × has_codemeta_file
          ✓ has_zenodo_badge
          × has_zenodo_metadata_file
    (5/5) checklist
          ✓ has_core_infrastructures_badge

If your README already has the fair-software badge, you'll see some output like this:

.. code:: console

    Calculated compliance: ● ● ○ ● ●

    Expected badge is equal to the actual badge. It's all good.

If your README doesn't have the fair-software badge yet, or its compliance is different from what's been calculated,
you'll see output like this:

.. code:: console

    Calculated compliance: ● ● ○ ○ ○

    It seems you have not yet added the fair-software.eu badge to
    your README.md. You can do so by pasting the following snippet:

    [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8B%20%20%E2%97%8B-orange)](https://fair-software.eu)

When you get this message, just copy-and-paste the suggested badge into your README.

Some examples of badges
-----------------------

The color of the badge depends on the level of compliance; the pattern of filled and empty circles will vary depending
on which recommendations the repository complies with.

Each circle represents one of the recommendations, meaning the first symbol represents the first recommendation, *Use a
publicly accessible repository with version control*, the second symbol represents the second recommendation, and so on.
You can find more information about the recommendations on fair-software.eu_.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8B%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8B-red

The state of the third circle indicates the software has been registered in a community registry. Since the repository
only complies with one of the recommendations, this badge gets a red color.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-orange

The repository with this badge complies with 3 out of 5 recommendations, hence its color is orange. From the open/closed
state of the circles, it is a publicly accessible repository with version control. It has been registered in a community
registry, and it contains citation information. There is no license in this repository, and the project does not use a
checklist.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow

Almost complete compliance yields a yellow badge. The corresponding repository meets all the recommendations except
the one that calls for adding a checklist.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green

Perfect compliance!

More options
------------

There are some command line options to the executable. You can see them using:

.. code:: console

    howfairis --help

Which then shows something like:

.. code:: console

    Usage: howfairis [OPTIONS] [URL]

      Determine compliance with recommendations from fair-software.eu for the
      repository at URL. The following code repository platforms are supported:

      * https://github.com

      * https://gitlab.com (not including any self-hosted instances)

    Options:
      -b, --branch TEXT               Which git branch to use. Also accepts other
                                      git references like SHA or tag.

      -u, --user-config-filename PATH
                                      Name of the configuration file to control
                                      howfairis'es behavior. The configuration
                                      file needs to be present on the local system
                                      and can include a relative path.

      -d, --show-default-config       Show default configuration and exit.
      -i, --ignore-repo-config        Ignore any configuration files on the
                                      remote.

      -p, --path TEXT                 Relative path (on the remote). Use this if
                                      you want howfairis to look for a README and
                                      a configuration file in a subdirectory.

      -q, --quiet                     Use this flag to disable all printing except
                                      errors.

      -r, --repo-config-filename TEXT
                                      Name of the configuration file to control
                                      howfairis'es behavior. The configuration
                                      file needs to be on the remote, and takes
                                      into account the value of --branch and
                                      --path. Default: .howfairis.yml

      -t, --show-trace                Show full traceback on errors.
      -v, --version                   Show version and exit.
      -h, --help                      Show this message and exit.

Configuration file
^^^^^^^^^^^^^^^^^^

Each category of checks can be skipped using a configuration file. This file needs to be present at ``URL``, taking into
account the values passed with ``--path`` and with ``--repo-config-filename``.

The configuration file should follow the voluptuous_ schema laid out in schema.py_:

.. code:: python

    schema = {
        Optional("skip_repository_checks_reason"): Any(str, None),
        Optional("skip_license_checks_reason"): Any(str, None),
        Optional("skip_registry_checks_reason"): Any(str, None),
        Optional("skip_citation_checks_reason"): Any(str, None),
        Optional("skip_checklist_checks_reason"): Any(str, None),
        Optional("ignore_commented_badges"): Any(bool, None)
    }

For example, the following is a valid configuration file document:

.. code:: yaml

    ## Uncomment a line if you want to skip a given category of checks

    #skip_repository_checks_reason: <reason for skipping goes here>
    #skip_license_checks_reason: <reason for skipping goes here>
    #skip_registry_checks_reason: <reason for skipping goes here>
    #skip_citation_checks_reason: <reason for skipping goes here>
    skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

    ignore_commented_badges: false


The manual override will be reflected in the output, as follows:

.. code:: console

    (1/5) repository
          ✓ has_open_repository
    (2/5) license
          ✓ has_license
    (3/5) registry
          × has_ascl_badge
          × has_bintray_badge
          × has_conda_badge
          × has_cran_badge
          × has_crates_badge
          × has_maven_badge
          × has_npm_badge
          ✓ has_pypi_badge
          × has_rsd_badge
          × is_on_github_marketplace
    (4/5) citation
          × has_citation_file
          ✓ has_citationcff_file
          × has_codemeta_file
          ✓ has_zenodo_badge
          ✓ has_zenodo_metadata_file
    (5/5) checklist
          ✓ skipped (reason: I'm using the Codacy dashboard to guide my development)

Rate limit
^^^^^^^^^^

By default ``howfairis`` uses anonymous requests to the API of the source code platforms.
However when a lot of repositories are checked you will exceed the rate limit of those APIs and checks will fail.
To increase the rate limit you need to use authenticated requests.
Your username and token can be passed to ``howfairis`` using environment variables called ``APIKEY_GITHUB`` and ``APIKEY_GITLAB``.
The format of the environment variable values are:

.. code-block:: shell

  export APIKEY_GITHUB=<user who made the token>:<personal access token>
  export APIKEY_GITLAB=<user who made the token>:<personal access token>

Generation of personal access tokens are explained on `GitHub documentation <https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token>`_ and `GitLab documentation <https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html#creating-a-personal-access-token>`_.
No scopes have to be selected, being authenticated is enough to get higher rate limit.

Contributing
------------

If you want to contribute to the development of howfairis, have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

If you're looking for developer documentation, go `here <README.dev.rst>`_.

.. _fair-software.eu: https://fair-software.eu
.. _voluptuous: https://pypi.org/project/voluptuous/
.. _schema.py: https://github.com/fair-software/howfairis/blob/master/howfairis/schema.py

Credits
-------

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.
..
   This is a comment in rst

   The comment continues

This is normal text

This is more normal text
This is normal text

This is more normal text..
   This is a comment in rst

These badges are nested deeper in the DOM than regular text:

1. .. image:: https://img.shields.io/badge/ascl-1410.001-red
      :target: https://ascl.net/1410.001
2. .. image:: https://bestpractices.coreinfrastructure.org/projects/4630/badge
      :target: https://bestpractices.coreinfrastructure.org/en/projects/4630
howfairis

Python package to analyze a GitHub or GitLab repository's compliance with the
fair-software.eu_ recommendations.

fair-software.eu_

Python package to analyze a GitHub or GitLab repository's compliance with the
fair-software.eu_ recommendations.

Badges

fair-software.nl recommendations

(1/5) code repository

|github repo badge|

(2/5) license

|github license badge|

(3/5) community registry

|pypi badge|

(4/5) citation

|zenodo badge|

(5/5) checklist

|core infrastructures badge|

overall

|fair-software badge|

**Other best practices**

Documentation

|readthedocs badge|

Supported Python versions

|python versions badge|

Code quality

|sonarcloud quality badge|

Code coverage of unit tests

|sonarcloud coverage badge|

DockerHub

|dockerhub badge|

**GitHub Actions**

Citation metadata consistency

|workflow cffconvert badge|

Unit tests

|workflow tests badge|

Live tests (triggered manually)

|workflow livetests badge|

image:: https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue
   :target: https://github.com/fair-software/howfairis

image:: https://img.shields.io/github/license/fair-software/howfairis
   :target: https://github.com/fair-software/howfairis

image:: https://img.shields.io/pypi/v/howfairis.svg?colorB=blue
   :target: https://pypi.python.org/pypi/howfairis/

image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4017908.svg
   :target: https://doi.org/10.5281/zenodo.4017908

image:: https://bestpractices.coreinfrastructure.org/projects/4630/badge
   :target: https://bestpractices.coreinfrastructure.org/en/projects/4630

image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green
   :target: https://fair-software.eu

image:: https://readthedocs.org/projects/howfairis/badge/?version=latest
   :target: https://howfairis.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

image:: https://img.shields.io/pypi/pyversions/howfairis.svg
   :target: https://pypi.python.org/pypi/howfairis

image:: https://sonarcloud.io/api/project_badges/measure?project=fair-software_howfairis&metric=alert_status
   :target: https://sonarcloud.io/dashboard?id=fair-software_howfairis
   :alt: Quality Gate Status

image:: https://sonarcloud.io/api/project_badges/measure?project=fair-software_howfairis&metric=coverage
   :target: https://sonarcloud.io/dashboard?id=fair-software_howfairis
   :alt: Coverage

image:: https://img.shields.io/docker/pulls/fairsoftware/howfairis
   :target: https://hub.docker.com/r/fairsoftware/howfairis
   :alt: Docker Pulls

image:: https://github.com/fair-software/howfairis/workflows/tests/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3Atests

image:: https://github.com/fair-software/howfairis/workflows/livetests/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3Alivetests

image:: https://github.com/fair-software/howfairis/workflows/metadata%20consistency/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3A%22metadata+consistency%22

Install

pip3 install --user howfairis

Verify that the install directory is on the ``PATH`` environment variable. If so,
you should be able to call the executable, like so:

``PATH``

Verify that the install directory is on the ``PATH`` environment variable. If so,
you should be able to call the executable, like so:

howfairis https://github.com/<owner>/<repo>

``howfairis``

``howfairis`` supports URLs from the following code repository platforms:

``https://github.com``

``https://gitlab.com``

``https://gitlab.com`` (not including self-hosted instances)

Docker

You can run howfairis Docker image using the command below.

docker pull fairsoftware/howfairis

You can run howfairis Docker image using the command below.

docker run --rm fairsoftware/howfairis --help

`--rm`

`--rm` argument will remove Docker container after execution.

See developer documentation to learn how to modify the Docker image.

Expected output

Depending on which repository you are doing the analysis for, the output
looks something like this:

Checking compliance with fair-software.eu...
url: https://github.com/fair-software/badge-test


(1/5)

Checking compliance with fair-software.eu...
url: https://github.com/fair-software/badge-test
(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      × has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      × has_zenodo_metadata_file
(5/5) checklist
      ✓ has_core_infrastructures_badge

repository
      ✓ has_open_repository


(2/5)

Checking compliance with fair-software.eu...
url: https://github.com/fair-software/badge-test
(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      × has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      × has_zenodo_metadata_file
(5/5) checklist
      ✓ has_core_infrastructures_badge

license
      ✓ has_license


(3/5)

Checking compliance with fair-software.eu...
url: https://github.com/fair-software/badge-test
(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      × has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      × has_zenodo_metadata_file
(5/5) checklist
      ✓ has_core_infrastructures_badge

registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace


(4/5)

Checking compliance with fair-software.eu...
url: https://github.com/fair-software/badge-test
(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      × has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      × has_zenodo_metadata_file
(5/5) checklist
      ✓ has_core_infrastructures_badge

citation
      × has_citation_file
      × has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      × has_zenodo_metadata_file


(5/5)

Checking compliance with fair-software.eu...
url: https://github.com/fair-software/badge-test
(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      × has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      × has_zenodo_metadata_file
(5/5) checklist
      ✓ has_core_infrastructures_badge

checklist
      ✓ has_core_infrastructures_badge

If your README already has the fair-software badge, you'll see some output like this:

Calculated compliance: ● ● ○ ● ●

Expected badge is equal to the actual badge. It's all good.

If your README doesn't have the fair-software badge yet, or its compliance is different from what's been calculated,
you'll see output like this:

Calculated compliance: ● ● ○ ○ ○

It seems you have not yet added the fair-software.eu badge to
your README.md. You can do so by pasting the following snippet:

[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8B%20%20%E2%97%8B-orange)](https://fair-software.eu)

When you get this message, just copy-and-paste the suggested badge into your README.

Some examples of badges

The color of the badge depends on the level of compliance; the pattern of filled and empty circles will vary depending
on which recommendations the repository complies with.

Each circle represents one of the recommendations, meaning the first symbol represents the first recommendation, *Use a
publicly accessible repository with version control*, the second symbol represents the second recommendation, and so on.
You can find more information about the recommendations on fair-software.eu_.

*Use a
publicly accessible repository with version control*

Each circle represents one of the recommendations, meaning the first symbol represents the first recommendation, *Use a
publicly accessible repository with version control*, the second symbol represents the second recommendation, and so on.
You can find more information about the recommendations on fair-software.eu_.

fair-software.eu_

Each circle represents one of the recommendations, meaning the first symbol represents the first recommendation, *Use a
publicly accessible repository with version control*, the second symbol represents the second recommendation, and so on.
You can find more information about the recommendations on fair-software.eu_.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8B%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8B-red


The state of the third circle indicates the software has been registered in a community registry. Since the repository
only complies with one of the recommendations, this badge gets a red color.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-orange


The repository with this badge complies with 3 out of 5 recommendations, hence its color is orange. From the open/closed
state of the circles, it is a publicly accessible repository with version control. It has been registered in a community
registry, and it contains citation information. There is no license in this repository, and the project does not use a
checklist.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow


Almost complete compliance yields a yellow badge. The corresponding repository meets all the recommendations except
the one that calls for adding a checklist.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green


Perfect compliance!

More options

There are some command line options to the executable. You can see them using:

howfairis --help

Which then shows something like:

Usage: howfairis [OPTIONS] [URL]

  Determine compliance with recommendations from fair-software.eu for the
  repository at URL. The following code repository platforms are supported:

  * https://github.com

  * https://gitlab.com (not including any self-hosted instances)

Options:
  -b, --branch TEXT               Which git branch to use. Also accepts other
                                  git references like SHA or tag.

  -u, --user-config-filename PATH
                                  Name of the configuration file to control
                                  howfairis'es behavior. The configuration
                                  file needs to be present on the local system
                                  and can include a relative path.

  -d, --show-default-config       Show default configuration and exit.
  -i, --ignore-repo-config        Ignore any configuration files on the
                                  remote.

  -p, --path TEXT                 Relative path (on the remote). Use this if
                                  you want howfairis to look for a README and
                                  a configuration file in a subdirectory.

  -q, --quiet                     Use this flag to disable all printing except
                                  errors.

  -r, --repo-config-filename TEXT
                                  Name of the configuration file to control
                                  howfairis'es behavior. The configuration
                                  file needs to be on the remote, and takes
                                  into account the value of --branch and
                                  --path. Default: .howfairis.yml

  -t, --show-trace                Show full traceback on errors.
  -v, --version                   Show version and exit.
  -h, --help                      Show this message and exit.

Configuration file

Each category of checks can be skipped using a configuration file. This file needs to be present at ``URL``, taking into
account the values passed with ``--path`` and with ``--repo-config-filename``.

``URL``

Each category of checks can be skipped using a configuration file. This file needs to be present at ``URL``, taking into
account the values passed with ``--path`` and with ``--repo-config-filename``.

``--path``

Each category of checks can be skipped using a configuration file. This file needs to be present at ``URL``, taking into
account the values passed with ``--path`` and with ``--repo-config-filename``.

``--repo-config-filename``

Each category of checks can be skipped using a configuration file. This file needs to be present at ``URL``, taking into
account the values passed with ``--path`` and with ``--repo-config-filename``.

The configuration file should follow the voluptuous_ schema laid out in schema.py_:

voluptuous_

The configuration file should follow the voluptuous_ schema laid out in schema.py_:

schema.py_

The configuration file should follow the voluptuous_ schema laid out in schema.py_:

schema

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

=

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

{

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Optional

(

"skip_repository_checks_reason"

):

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Any

(

str

,

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

None

),

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Optional

(

"skip_license_checks_reason"

):

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Any

(

str

,

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

None

),

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Optional

(

"skip_registry_checks_reason"

):

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Any

(

str

,

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

None

),

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Optional

(

"skip_citation_checks_reason"

):

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Any

(

str

,

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

None

),

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Optional

(

"skip_checklist_checks_reason"

):

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Any

(

str

,

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

None

),

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Optional

(

"ignore_commented_badges"

):

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

Any

(

bool

,

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

None

)

schema = {
    Optional("skip_repository_checks_reason"): Any(str, None),
    Optional("skip_license_checks_reason"): Any(str, None),
    Optional("skip_registry_checks_reason"): Any(str, None),
    Optional("skip_citation_checks_reason"): Any(str, None),
    Optional("skip_checklist_checks_reason"): Any(str, None),
    Optional("ignore_commented_badges"): Any(bool, None)
}

}

For example, the following is a valid configuration file document:

## Uncomment a line if you want to skip a given category of checks

## Uncomment a line if you want to skip a given category of checks

#skip_repository_checks_reason: <reason for skipping goes here>
#skip_license_checks_reason: <reason for skipping goes here>
#skip_registry_checks_reason: <reason for skipping goes here>
#skip_citation_checks_reason: <reason for skipping goes here>
skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

ignore_commented_badges: false

#skip_repository_checks_reason: <reason for skipping goes here>

## Uncomment a line if you want to skip a given category of checks

#skip_repository_checks_reason: <reason for skipping goes here>
#skip_license_checks_reason: <reason for skipping goes here>
#skip_registry_checks_reason: <reason for skipping goes here>
#skip_citation_checks_reason: <reason for skipping goes here>
skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

ignore_commented_badges: false

#skip_license_checks_reason: <reason for skipping goes here>

## Uncomment a line if you want to skip a given category of checks

#skip_repository_checks_reason: <reason for skipping goes here>
#skip_license_checks_reason: <reason for skipping goes here>
#skip_registry_checks_reason: <reason for skipping goes here>
#skip_citation_checks_reason: <reason for skipping goes here>
skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

ignore_commented_badges: false

#skip_registry_checks_reason: <reason for skipping goes here>

## Uncomment a line if you want to skip a given category of checks

#skip_repository_checks_reason: <reason for skipping goes here>
#skip_license_checks_reason: <reason for skipping goes here>
#skip_registry_checks_reason: <reason for skipping goes here>
#skip_citation_checks_reason: <reason for skipping goes here>
skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

ignore_commented_badges: false

#skip_citation_checks_reason: <reason for skipping goes here>

## Uncomment a line if you want to skip a given category of checks

#skip_repository_checks_reason: <reason for skipping goes here>
#skip_license_checks_reason: <reason for skipping goes here>
#skip_registry_checks_reason: <reason for skipping goes here>
#skip_citation_checks_reason: <reason for skipping goes here>
skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

ignore_commented_badges: false

skip_checklist_checks_reason

:

## Uncomment a line if you want to skip a given category of checks

#skip_repository_checks_reason: <reason for skipping goes here>
#skip_license_checks_reason: <reason for skipping goes here>
#skip_registry_checks_reason: <reason for skipping goes here>
#skip_citation_checks_reason: <reason for skipping goes here>
skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

ignore_commented_badges: false

"I'm

 

using

 

the

 

Codacy

 

dashboard

 

to

 

guide

 

my

 

development"

## Uncomment a line if you want to skip a given category of checks

#skip_repository_checks_reason: <reason for skipping goes here>
#skip_license_checks_reason: <reason for skipping goes here>
#skip_registry_checks_reason: <reason for skipping goes here>
#skip_citation_checks_reason: <reason for skipping goes here>
skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

ignore_commented_badges: false

ignore_commented_badges

:

## Uncomment a line if you want to skip a given category of checks

#skip_repository_checks_reason: <reason for skipping goes here>
#skip_license_checks_reason: <reason for skipping goes here>
#skip_registry_checks_reason: <reason for skipping goes here>
#skip_citation_checks_reason: <reason for skipping goes here>
skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

ignore_commented_badges: false

false

The manual override will be reflected in the output, as follows:

(1/5)

(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      ✓ has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      ✓ has_zenodo_metadata_file
(5/5) checklist
      ✓ skipped (reason: I'm using the Codacy dashboard to guide my development)

repository
      ✓ has_open_repository


(2/5)

(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      ✓ has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      ✓ has_zenodo_metadata_file
(5/5) checklist
      ✓ skipped (reason: I'm using the Codacy dashboard to guide my development)

license
      ✓ has_license


(3/5)

(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      ✓ has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      ✓ has_zenodo_metadata_file
(5/5) checklist
      ✓ skipped (reason: I'm using the Codacy dashboard to guide my development)

registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace


(4/5)

(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      ✓ has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      ✓ has_zenodo_metadata_file
(5/5) checklist
      ✓ skipped (reason: I'm using the Codacy dashboard to guide my development)

citation
      × has_citation_file
      ✓ has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      ✓ has_zenodo_metadata_file


(5/5)

(1/5) repository
      ✓ has_open_repository
(2/5) license
      ✓ has_license
(3/5) registry
      × has_ascl_badge
      × has_bintray_badge
      × has_conda_badge
      × has_cran_badge
      × has_crates_badge
      × has_maven_badge
      × has_npm_badge
      ✓ has_pypi_badge
      × has_rsd_badge
      × is_on_github_marketplace
(4/5) citation
      × has_citation_file
      ✓ has_citationcff_file
      × has_codemeta_file
      ✓ has_zenodo_badge
      ✓ has_zenodo_metadata_file
(5/5) checklist
      ✓ skipped (reason: I'm using the Codacy dashboard to guide my development)

checklist
      ✓ skipped (reason: I'm using the Codacy dashboard to guide my development)

Contributing

If you want to contribute to the development of howfairis, have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

`contribution guidelines <CONTRIBUTING.rst>`_

 <CONTRIBUTING.rst>

If you want to contribute to the development of howfairis, have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

If you're looking for developer documentation, go `here <README.dev.rst>`_.

`here <README.dev.rst>`_

 <README.dev.rst>

If you're looking for developer documentation, go `here <README.dev.rst>`_.

.. _fair-software.eu: https://fair-software.eu

.. _voluptuous: https://pypi.org/project/voluptuous/

.. _schema.py: https://github.com/fair-software/howfairis/blob/master/howfairis/schema.py

Credits

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.

`Cookiecutter <https://github.com/audreyr/cookiecutter>`_

 <https://github.com/audreyr/cookiecutter>

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.

`NLeSC/python-template <https://github.com/NLeSC/python-template>`_

 <https://github.com/NLeSC/python-template>

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.This is normal text

This is more normal text..
   This is a comment in rst

This is normal text

This is more normal text
howfairis
=========

|

Python package to analyze a GitHub or GitLab repository's compliance with the
fair-software.eu_ recommendations.

Badges
------

..
    This is a comment in rst

====================================================== ============================
fair-software.nl recommendations
====================================================== ============================
(1/5) code repository                                  |github repo badge|
(2/5) license                                          |github license badge|
(3/5) community registry                               |pypi badge|
(4/5) citation                                         |zenodo badge|
(5/5) checklist                                        |core infrastructures badge|
overall                                                |fair-software badge|
**Other best practices**
Documentation                                          |readthedocs badge|
Supported Python versions                              |python versions badge|
Code quality                                           |sonarcloud quality badge|
Code coverage of unit tests                            |sonarcloud coverage badge|
DockerHub                                              |dockerhub badge|
**GitHub Actions**
Citation metadata consistency                          |workflow cffconvert badge|
Unit tests                                             |workflow tests badge|
Live tests (triggered manually)                        |workflow livetests badge|
====================================================== ============================

.. |github repo badge| image:: https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue
   :target: https://github.com/fair-software/howfairis

.. |github license badge| image:: https://img.shields.io/github/license/fair-software/howfairis
   :target: https://github.com/fair-software/howfairis

.. |pypi badge| image:: https://img.shields.io/pypi/v/howfairis.svg?colorB=blue
   :target: https://pypi.python.org/pypi/howfairis/

.. |zenodo badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4017908.svg
   :target: https://doi.org/10.5281/zenodo.4017908

.. |core infrastructures badge| image:: https://bestpractices.coreinfrastructure.org/projects/4630/badge
   :target: https://bestpractices.coreinfrastructure.org/en/projects/4630

.. |fair-software badge| image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green
   :target: https://fair-software.eu

.. |readthedocs badge| image:: https://readthedocs.org/projects/howfairis/badge/?version=latest
   :target: https://howfairis.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |python versions badge| image:: https://img.shields.io/pypi/pyversions/howfairis.svg
   :target: https://pypi.python.org/pypi/howfairis

.. |sonarcloud quality badge| image:: https://sonarcloud.io/api/project_badges/measure?project=fair-software_howfairis&metric=alert_status
   :target: https://sonarcloud.io/dashboard?id=fair-software_howfairis
   :alt: Quality Gate Status

.. |sonarcloud coverage badge| image:: https://sonarcloud.io/api/project_badges/measure?project=fair-software_howfairis&metric=coverage
   :target: https://sonarcloud.io/dashboard?id=fair-software_howfairis
   :alt: Coverage

.. |dockerhub badge| image:: https://img.shields.io/docker/pulls/fairsoftware/howfairis
   :target: https://hub.docker.com/r/fairsoftware/howfairis
   :alt: Docker Pulls

.. |workflow tests badge| image:: https://github.com/fair-software/howfairis/workflows/tests/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3Atests

.. |workflow livetests badge| image:: https://github.com/fair-software/howfairis/workflows/livetests/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3Alivetests

.. |workflow cffconvert badge| image:: https://github.com/fair-software/howfairis/workflows/metadata%20consistency/badge.svg
   :target: https://github.com/fair-software/howfairis/actions?query=workflow%3A%22metadata+consistency%22

Install
-------

.. code:: console

    pip3 install --user howfairis

Verify that the install directory is on the ``PATH`` environment variable. If so,
you should be able to call the executable, like so:

.. code:: console

    howfairis https://github.com/<owner>/<repo>


``howfairis`` supports URLs from the following code repository platforms:

1. ``https://github.com``
2. ``https://gitlab.com`` (not including self-hosted instances)

Docker
---------------

You can run howfairis Docker image using the command below.

.. code:: console

    docker pull fairsoftware/howfairis

You can run howfairis Docker image using the command below.

.. code:: console

    docker run --rm fairsoftware/howfairis --help

`--rm` argument will remove Docker container after execution.

See developer documentation to learn how to modify the Docker image.

Expected output
---------------

Depending on which repository you are doing the analysis for, the output
looks something like this:

.. code:: console

    Checking compliance with fair-software.eu...
    url: https://github.com/fair-software/badge-test
    (1/5) repository
          ✓ has_open_repository
    (2/5) license
          ✓ has_license
    (3/5) registry
          × has_ascl_badge
          × has_bintray_badge
          × has_conda_badge
          × has_cran_badge
          × has_crates_badge
          × has_maven_badge
          × has_npm_badge
          ✓ has_pypi_badge
          × has_rsd_badge
          × is_on_github_marketplace
    (4/5) citation
          × has_citation_file
          × has_citationcff_file
          × has_codemeta_file
          ✓ has_zenodo_badge
          × has_zenodo_metadata_file
    (5/5) checklist
          ✓ has_core_infrastructures_badge

If your README already has the fair-software badge, you'll see some output like this:

.. code:: console

    Calculated compliance: ● ● ○ ● ●

    Expected badge is equal to the actual badge. It's all good.

If your README doesn't have the fair-software badge yet, or its compliance is different from what's been calculated,
you'll see output like this:

.. code:: console

    Calculated compliance: ● ● ○ ○ ○

    It seems you have not yet added the fair-software.eu badge to
    your README.md. You can do so by pasting the following snippet:

    [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8B%20%20%E2%97%8B-orange)](https://fair-software.eu)

When you get this message, just copy-and-paste the suggested badge into your README.

Some examples of badges
-----------------------

The color of the badge depends on the level of compliance; the pattern of filled and empty circles will vary depending
on which recommendations the repository complies with.

Each circle represents one of the recommendations, meaning the first symbol represents the first recommendation, *Use a
publicly accessible repository with version control*, the second symbol represents the second recommendation, and so on.
You can find more information about the recommendations on fair-software.eu_.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8B%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8B-red

The state of the third circle indicates the software has been registered in a community registry. Since the repository
only complies with one of the recommendations, this badge gets a red color.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-orange

The repository with this badge complies with 3 out of 5 recommendations, hence its color is orange. From the open/closed
state of the circles, it is a publicly accessible repository with version control. It has been registered in a community
registry, and it contains citation information. There is no license in this repository, and the project does not use a
checklist.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow

Almost complete compliance yields a yellow badge. The corresponding repository meets all the recommendations except
the one that calls for adding a checklist.

.. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green

Perfect compliance!

More options
------------

There are some command line options to the executable. You can see them using:

.. code:: console

    howfairis --help

Which then shows something like:

.. code:: console

    Usage: howfairis [OPTIONS] [URL]

      Determine compliance with recommendations from fair-software.eu for the
      repository at URL. The following code repository platforms are supported:

      * https://github.com

      * https://gitlab.com (not including any self-hosted instances)

    Options:
      -b, --branch TEXT               Which git branch to use. Also accepts other
                                      git references like SHA or tag.

      -u, --user-config-filename PATH
                                      Name of the configuration file to control
                                      howfairis'es behavior. The configuration
                                      file needs to be present on the local system
                                      and can include a relative path.

      -d, --show-default-config       Show default configuration and exit.
      -i, --ignore-repo-config        Ignore any configuration files on the
                                      remote.

      -p, --path TEXT                 Relative path (on the remote). Use this if
                                      you want howfairis to look for a README and
                                      a configuration file in a subdirectory.

      -q, --quiet                     Use this flag to disable all printing except
                                      errors.

      -r, --repo-config-filename TEXT
                                      Name of the configuration file to control
                                      howfairis'es behavior. The configuration
                                      file needs to be on the remote, and takes
                                      into account the value of --branch and
                                      --path. Default: .howfairis.yml

      -t, --show-trace                Show full traceback on errors.
      -v, --version                   Show version and exit.
      -h, --help                      Show this message and exit.

Configuration file
^^^^^^^^^^^^^^^^^^

Each category of checks can be skipped using a configuration file. This file needs to be present at ``URL``, taking into
account the values passed with ``--path`` and with ``--repo-config-filename``.

The configuration file should follow the voluptuous_ schema laid out in schema.py_:

.. code:: python

    schema = {
        Optional("skip_repository_checks_reason"): Any(str, None),
        Optional("skip_license_checks_reason"): Any(str, None),
        Optional("skip_registry_checks_reason"): Any(str, None),
        Optional("skip_citation_checks_reason"): Any(str, None),
        Optional("skip_checklist_checks_reason"): Any(str, None),
        Optional("ignore_commented_badges"): Any(bool, None)
    }

For example, the following is a valid configuration file document:

.. code:: yaml

    ## Uncomment a line if you want to skip a given category of checks

    #skip_repository_checks_reason: <reason for skipping goes here>
    #skip_license_checks_reason: <reason for skipping goes here>
    #skip_registry_checks_reason: <reason for skipping goes here>
    #skip_citation_checks_reason: <reason for skipping goes here>
    skip_checklist_checks_reason: "I'm using the Codacy dashboard to guide my development"

    ignore_commented_badges: false


The manual override will be reflected in the output, as follows:

.. code:: console

    (1/5) repository
          ✓ has_open_repository
    (2/5) license
          ✓ has_license
    (3/5) registry
          × has_ascl_badge
          × has_bintray_badge
          × has_conda_badge
          × has_cran_badge
          × has_crates_badge
          × has_maven_badge
          × has_npm_badge
          ✓ has_pypi_badge
          × has_rsd_badge
          × is_on_github_marketplace
    (4/5) citation
          × has_citation_file
          ✓ has_citationcff_file
          × has_codemeta_file
          ✓ has_zenodo_badge
          ✓ has_zenodo_metadata_file
    (5/5) checklist
          ✓ skipped (reason: I'm using the Codacy dashboard to guide my development)

Contributing
------------

If you want to contribute to the development of howfairis, have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

If you're looking for developer documentation, go `here <README.dev.rst>`_.

.. _fair-software.eu: https://fair-software.eu
.. _voluptuous: https://pypi.org/project/voluptuous/
.. _schema.py: https://github.com/fair-software/howfairis/blob/master/howfairis/schema.py

Credits
-------

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.
This is normal text

This is more normal textThese badges are nested deeper in the DOM than regular text:

.. image:: https://img.shields.io/badge/ascl-1410.001-red
   :target: https://ascl.net/1410.001

.. image:: https://bestpractices.coreinfrastructure.org/projects/4630/badge
   :target: https://bestpractices.coreinfrastructure.org/en/projects/4630..
   _This: is a comment!

..
   [and] this!

..
   this:: too!

..
   |even| this:: !

.. .. |pypi badge| image:: https://img.shields.io/pypi/v/howfairis.svg?colorB=blue
..    :target: https://pypi.python.org/pypi/howfairis/

This is normal text

This is more normal text

..
   Examples from https://stackoverflow.com/a/4783854.. This is a comment in rst

This is normal text

This is more normal text
This is normal text

This is more normal text# README.rst mocked content with all badges

1. .. image:: https://img.shields.io/badge/ascl-1410.001-red
      :target: https://ascl.net/1410.001
2. .. image:: https://bestpractices.coreinfrastructure.org/projects/4630/badge
      :target: https://bestpractices.coreinfrastructure.org/en/projects/4630
3. .. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4017908.svg
      :target: https://doi.org/10.5281/zenodo.4017908
4. .. image:: https://anaconda.org/nlesc/matchms/badges/installer/conda.svg
      :target: https://conda.anaconda.org/nlesc
5. .. image:: https://img.shields.io/badge/rsd-matchms-00a3e3.svg
      :target: https://www.research-software.nl/software/matchms
6. .. image:: https://img.shields.io/pypi/v/matchms?color=blue
      :target: https://pypi.org/project/matchms/
7. .. image:: https://api.bintray.com/packages/nlesc/xenon/xenon/images/download.svg
      :target: https://bintray.com/nlesc/xenon/xenon/_latestVersion
8. .. image:: https://cranlogs.r-pkg.org/badges/grand-total/GGIR
      :target: https://cran.r-project.org/package=GGIR
9. .. image:: https://img.shields.io/crates/v/tensorflow.svg
      :target: https://crates.io/crates/tensorflow
10. .. image:: https://img.shields.io/maven-central/v/com.spotify/dockerfile-maven.svg
      :target: https://search.maven.org/#search%7Cga%7C1%7Cg%3A%22com.spotify%22%20dockerfile-maven
11. .. image:: https://img.shields.io/npm/v/cnpm.svg?style=flat
      :target: https://npmjs.org/package/cnpm
12. .. image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green
      :target: https://fair-software.eu
CLI Reference
=============

.. click:: howfairis.cli.cli:cli
   :prog: howfairis
   :nested: full
Welcome to howfairis's documentation!
=====================================

.. toctree::
   :maxdepth: 1

   cli
   API Reference <apidocs/howfairis.rst>

.. include:: ../README.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
