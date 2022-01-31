# Project Setup

Here we provide some details about the project setup. Most of the choices are explained in the
[guide](https://guide.esciencecenter.nl). Links to the relevant sections are included below. Feel free to remove this
text when the development of the software package takes off.

For a quick reference on software development, we refer to [the software guide
checklist](https://guide.esciencecenter.nl/#/best_practices/checklist).

## Version control

Once your Python package is created, put it under [version
control](https://guide.esciencecenter.nl/#/best_practices/version_control)! We recommend using
[git](http://git-scm.com/) and [github](https://github.com/).

```shell
cd machine-learning
git init
git add --all
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/online-behaviour/machine-learning
```

Go to
[https://github.com/online-behaviour?tab=repositories](https://github.com/online-behaviour?tab=repositories)
and create a new repository named machine-learning as an empty repository, then:

```shell
git push --set-upstream origin main
```

## Python versions

This repository is set up with Python versions:

- 3.6
- 3.7
- 3.8
- 3.9

Add or remove Python versions based on project requirements. See [the
guide](https://guide.esciencecenter.nl/#/best_practices/language_guides/python) for more information about Python
versions.

## Package management and dependencies

You can use either pip or conda for installing dependencies and package management. This repository does not force you
to use one or the other, as project requirements differ. For advice on what to use, please check [the relevant section
of the
guide](https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=dependencies-and-package-management).

- Runtime dependencies should be added to `setup.cfg` in the `install_requires` list under `[options]`.
- Development dependencies should be added to `setup.cfg` in one of the lists under `[options.extras_require]`.

## Packaging/One command install

You can distribute your code using PyPI.
[The guide](https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=building-and-packaging-code) can
help you decide which tool to use for packaging.

## Testing and code coverage

- Tests should be put in the `tests` folder.
- The `tests` folder contains:
  - Example tests that you should replace with your own meaningful tests (file: `test_my_module.py`)
- The testing framework used is [PyTest](https://pytest.org)
  - [PyTest introduction](http://pythontesting.net/framework/pytest/pytest-introduction/)
  - PyTest is listed as a development dependency, and can thus be installed with `pip3 install --editable .[dev]`
- Tests can be run with `pytest`
  - This is configured in `setup.cfg`
- The project uses [GitHub action workflows](https://docs.github.com/en/actions) to automatically run tests on GitHub infrastructure against multiple Python versions
  - Workflows can be found in [`.github/workflows`](.github/workflows/)
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=testing)

## Documentation

- Documentation should be put in the [`docs/`](docs/) directory. The contents have been generated using `sphinx-quickstart` (Sphinx version 1.6.5).
- We recommend writing the documentation using Restructured Text (reST) and Google style docstrings.
  - [Restructured Text (reST) and Sphinx CheatSheet](http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html)
  - [Google style docstring examples](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).
- The documentation is set up with the ReadTheDocs Sphinx theme.
  - Check out its [configuration options](https://sphinx-rtd-theme.readthedocs.io/en/latest/).
- To generate HTML documentation, run `make html` in the `docs/` folder.
- To put the documentation on [ReadTheDocs](https://readthedocs.org), log in to your ReadTheDocs account, and import
  the repository (under 'My Projects').
  - Include the link to the documentation in your project's [README.md](README.md).
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=writingdocumentation)

## Coding style conventions and code quality

- Check your code style with `prospector`
- You may need run `pip install --editable .[dev]` first, to install the required dependencies
- You can use `yapf` to fix the readability of your code style and `isort` to format and group your imports
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=coding-style-conventions)

## Continuous code quality

- [Sonarcloud](https://sonarcloud.io/) is used to perform quality analysis and code coverage report on each push
- Sonarcloud must be configured for the analysis to work
  1. go to [Sonarcloud](https://sonarcloud.io/projects/create)
  2. login with your GitHub account
  3. add organization or reuse existing one
  4. set up repository
  5. go to [new code definition administration page](https://sonarcloud.io/project/new_code?id=online-behaviour_machine-learning) and select `Number of days` option
- The analysis will be run by [GitHub Action workflow](.github/workflows/sonarcloud.yml)
- To be able to run the analysis, a token must be created at [Sonarcloud account](https://sonarcloud.io/account/security/) and this token must be added as `SONAR_TOKEN` to [secrets on GitHub](https://github.com/online-behaviour/machine-learning/settings/secrets/actions)

## Package version number

- We recommend using [semantic versioning](https://guide.esciencecenter.nl/#/best_practices/releases?id=semantic-versioning).
- For convenience, the package version is stored in a single place: `machine-learning/.bumpversion.cfg`.
  For updating the version number, make sure the dev dependencies are installed and run `bumpversion patch`,
  `bumpversion minor`, or `bumpversion major` as appropriate.
- Don't forget to update the version number before [making a release](https://guide.esciencecenter.nl/#/best_practices/releases)!

## Publish on Python Package Index (PyPI)

To publish your package on PyPI, you need to create a [PyPI API token](https://pypi.org/help#apitoken) and
save it as a secret called `PYPI_TOKEN` on [Settings page](https://github.com/online-behaviour/machine-learning/settings/secrets/actions)

[Creating a release](https://github.com/online-behaviour/machine-learning/releases/new) on GitHub will trigger a [GitHub action workflow](.github/workflows/publish.yml) to publish the release on PyPI for you.

## Logging

- We recommend using the logging module for getting useful information from your module (instead of using print).
- The project is set up with a logging example.
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=logging)

## CHANGELOG.md

- Document changes to your software package
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/releases?id=changelogmd)

## CITATION.cff

- To allow others to cite your software, add a `CITATION.cff` file
- It only makes sense to do this once there is something to cite (e.g., a software release with a DOI).
- Follow the [making software citable](https://guide.esciencecenter.nl/#/citable_software/making_software_citable) section in the guide.

## CODE_OF_CONDUCT.md

- Information about how to behave professionally
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/documentation?id=code-of-conduct)

## CONTRIBUTING.md

- Information about how to contribute to this software package
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/documentation?id=contribution-guidelines)

## MANIFEST.in

- List non-Python files that should be included in a source distribution
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/language_guides/python?id=building-and-packaging-code)

## NOTICE

- List of attributions of this project and Apache-license dependencies
- [Relevant section in the guide](https://guide.esciencecenter.nl/#/best_practices/licensing?id=notice)
# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.1] - 1900-12-31

### Added

### Removed

### Changed

[Unreleased]: https://github.com/olivierlacan/keep-a-changelog/compare/v1.0.0...HEAD
[0.0.1]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/v0.0.1
# `machine_learning` developer documentation

If you're looking for user documentation, go [here](README.md).

## Development install

```shell
# Create a virtual environment, e.g. with
python3 -m venv env

# activate virtual environment
source env/bin/activate

# make sure to have a recent version of pip and setuptools
python3 -m pip install --upgrade pip setuptools

# (from the project root directory)
# install machine_learning as an editable package
python3 -m pip install --no-cache-dir --editable .
# install development dependencies
python3 -m pip install --no-cache-dir --editable .[dev]
```

Afterwards check that the install directory is present in the `PATH` environment variable.

## Running the tests

Running the tests requires an activated virtual environment with the development tools installed.

```shell
# unit tests
pytest
pytest tests/
```

## Running linters locally

For linting we will use [prospector](https://pypi.org/project/prospector/) and to sort imports we will use
[isort](https://pycqa.github.io/isort/). Running the linters requires an activated virtual environment with the
development tools installed.

```shell
# linter
prospector

# recursively check import style for the machine_learning module only
isort --recursive --check-only machine_learning

# recursively check import style for the machine_learning module only and show
# any proposed changes as a diff
isort --recursive --check-only --diff machine_learning

# recursively fix import style for the machine_learning module only
isort --recursive machine_learning
```

You can enable automatic linting with `prospector` and `isort` on commit by enabling the git hook from `.githooks/pre-commit`, like so:

```shell
git config --local core.hooksPath .githooks
```

## Generating the API docs

```shell
cd docs
make html
```

The documentation will be in `docs/_build/`

## Versioning

Bumping the version across all files is done with bumpversion, e.g.

```shell
bumpversion major
bumpversion minor
bumpversion patch
```

## Making a release

This section describes how to make a release in 3 parts:

1. preparation
1. making a release on PyPI
1. making a release on GitHub

### (1/3) Preparation

1.  Update the `CHANGELOG.md`
2.  Verify that the information in `CITATION.cff` is correct, and that `.zenodo.json` contains equivalent data
3.  Make sure the version has been updated.
4.  Run the unit tests with `pytest tests/`

### (2/3) PyPI

In a new terminal, without an activated virtual environment or an env directory:

```shell
# prepare a new directory
cd $(mktemp -d --tmpdir machine_learning.XXXXXX)

# fresh git clone ensures the release has the state of origin/main branch
git clone https://github.com/online-behaviour/machine-learning .

# prepare a clean virtual environment and activate it
python3 -m venv env
source env/bin/activate

# make sure to have a recent version of pip and setuptools
python3 -m pip install --upgrade pip setuptools

# install runtime dependencies and publishing dependencies
python3 -m pip install --no-cache-dir .
python3 -m pip install --no-cache-dir .[publishing]

# clean up any previously generated artefacts 
rm -rf machine_learning.egg-info
rm -rf dist

# create the source distribution and the wheel
python3 setup.py sdist bdist_wheel

# upload to test pypi instance (requires credentials)
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

Visit
[https://test.pypi.org/project/machine_learning](https://test.pypi.org/project/machine_learning)
and verify that your package was uploaded successfully. Keep the terminal open, we'll need it later.

In a new terminal, without an activated virtual environment or an env directory:

```shell
cd $(mktemp -d --tmpdir machine_learning-test.XXXXXX)

# prepare a clean virtual environment and activate it
python3 -m venv env
source env/bin/activate

# make sure to have a recent version of pip and setuptools
pip install --upgrade pip setuptools

# install from test pypi instance:
python3 -m pip -v install --no-cache-dir \
--index-url https://test.pypi.org/simple/ \
--extra-index-url https://pypi.org/simple machine_learning
```

Check that the package works as it should when installed from pypitest.

Then upload to pypi.org with:

```shell
# Back to the first terminal,
# FINAL STEP: upload to PyPI (requires credentials)
twine upload dist/*
```

### (3/3) GitHub

Don't forget to also make a release on GitHub. If your repository uses the GitHub-Zenodo integration this will also
trigger Zenodo into making a snapshot of your repository and sticking a DOI on it.
# Hercules: Tweet Retrieval and Machine Learning

This directory contains the software developed in the main
Machine Learning part of the project [Automated Analysis of
Online Behaviour on Social
Media](https://www.esciencecenter.nl/project/automated-analysis-of-online-behaviour-on-social-media),
a cooperation of the University of Groningen and the
Netherlands eScience Center. This project also has a
software repository regarding [finding
journalists](https://github.com/online-behaviour/find-journalists).

The software consist of a collection of Python scripts.
These can be divided in three groups:

1. tweet-fetching scripts
1. scripts related to the IEEE paper (Auckland)
1. scripts related to the Casablanca paper

There are also several unrelated scripts, which have been
left undocumented.

## Tweet-fetching scripts

### getTweetsUser.py

The Python script getTweetsUser.py can be used for obtaining
tweets from certain users. Run like:

```
./getTweetsUser.py barackobama realdonaldtrump > file
```

It will retrieve all available tweets from the specified
users and store these in the specified file. The tweets are
stored in the data format
[JSON](https://en.wikipedia.org/wiki/JSON). The command may
require several minutes to complete. 

The script needs two things to run:

First: you will need to install the Twitter package from:
[https://github.com/sixohsix/twitter](https://github.com/sixohsix/twitter)
The commands for this on Linux and MacOSX are:

```
git clone https://github.com/sixohsix/twitter
cd twitter
python setup.py build
sudo python setup.py install
```

Second: you need to store your Twitter account data in a file named
"definitions.py" in the same directory as getTweetsUser.py. The file
should contain the following lines:

```
# twitter.com authentication keys
token = "???"
token_secret = "???"
consumer_key = "???"
consumer_secret = "???"
```

Replace the strings "???" with the key information from 
https://apps.twitter.com , see
https://www.slickremix.com/docs/how-to-get-api-keys-and-tokens-for-twitter/
for instructions

### getTweetText.py

The Python script getTweetText.py can be used for extracting
the tweets from the JSON output of getTweetsUser.py:

```
./getTweetText.py < getTweetsUser.py.out > file
```

## Scripts related to the IEEE paper (Auckland)

> Erik Tjong Kim Sang, Herbert Kruitbosch, Marcel Broersma and
> Marc Esteve del Valle, Determining the function of political
> tweets. In: Proceedings of the 13th IEEE International
> Conference on eScience (eScience 2017), IEEE, Auckland, New
> Zealand, 2017, pages 438-439, ISBN 978-1-5386-2686-3,
> doi:10.1109/eScience.2017.60. 
> ([PDF](https://ifarm.nl/erikt/papers/2017-escience.pdf),
> [bibtex](https://ifarm.nl/erikt/papers/2017-escience.txt)]

First, the data needs to be converted to the format required
by the machine learner [fasttext](https://github.com/facebookresearch/fastText).
We use tokenized text preceded by the class label, for 
example *__label__1 this is a tweet !*:

```
for FILE in test train
do
   ./expandReplies.py -t dutch-2012.$FILE.csv -r EMPTY |\
      cut -d' ' -f1,4- | sed 's/ RAWTEXT /*$/' > dutch-2012.$FILE.txt
done
```

Note that the data files with annotated tweets 
(dutch-2012.*) are unavailable.

Next, [fasttext](https://github.com/facebookresearch/fastText)
can be applied to the data:

```
fasttext supervised -input dutch-2012.train.txt -output MODEL \
   -dim 5 -minCount 300
fasttext predict MODEL.bin dutch-2012.test.txt |\
   paste -d ' ' - dutch-2012.test.txt | cut -d' ' -f1,2 |
      ./eval.py | head -1 | rev | sed 's/^ *//' | cut -d' ' -f1 | rev
```

For most of the experiments mentioned in Table II of the
paper, these two commands can be reused with a different
training file. Only the language modeling experiments
require an extra step, for creating the language models:

```
fasttext skipgram -input EXTRADATA -output VECTORS -dim 5 \
   -minCount 300
fasttext supervised -input dutch-2012.train.txt -output MODEL \
   -dim 5 -minCount 300 -pretrainedVectors VECTORS.vec
fasttext predict MODEL.bin dutch-2012.test.txt |\
   paste -d ' ' - dutch-2012.test.txt | cut -d' ' -f1,2 |\
   ./eval.py | head -1 | rev | sed 's/^ *//' | cut -d' ' -f1 | rev
```

We always remove the labels from the EXTRADATA files.

## Scripts related to the Casablanca paper

> Erik Tjong Kim Sang, Marc Esteve del Valle, Herbert Kruitbosch,
> and Marcel Broersma, Active Learning for Classifying Political 
> Tweets. In: Proceedings of the International Conference on
> Natural Language, Signal and Speech Processing (ICNLSSP),
> Casablanca, Morocco, 2017.

The experiments related to Figure 1 and Table 1 of the
paper, were performed with the bash script `run.sh`.

After annotating a file for active learning, the next data
file was generated with the bash script `run-make-batch`.

## Contact

Erik Tjong Kim Sang, e.tjongkimsang(at)esciencecenter.nl

# Information added by the python template

## Badges

| fair-software.eu recommendations | |
| :-- | :--  |
| (1/5) code repository              | [![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/online-behaviour/machine-learning) |
| (2/5) license                      | [![github license badge](https://img.shields.io/github/license/online-behaviour/machine-learning)](https://github.com/online-behaviour/machine-learning) |
| (3/5) community registry           | [![Research Software Directory](https://img.shields.io/badge/rsd-Research%20Software%20Directory-00a3e3.svg)](https://www.research-software.nl/software/online-behaviour-machine-learning) |
| (4/5) citation                     | [![DOI](https://zenodo.org/badge/87834727.svg)](https://zenodo.org/badge/latestdoi/87834727) |
| (5/5) checklist                    | [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/4837/badge)](https://bestpractices.coreinfrastructure.org/projects/4837) |
| howfairis                            | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu) |
| **Other best practices**           | &nbsp; |
| Static analysis              | [![workflow scq badge](https://sonarcloud.io/api/project_badges/measure?project=online-behaviour_machine-learning&metric=alert_status)](https://sonarcloud.io/dashboard?id=online-behaviour_machine-learning) |
| Coverage              | [![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=online-behaviour_machine-learning&metric=coverage)](https://sonarcloud.io/dashboard?id=online-behaviour_machine-learning) |
| **GitHub Actions**                 | &nbsp; |
| Build                              | [![build](https://github.com/online-behaviour/machine-learning/actions/workflows/build.yml/badge.svg)](https://github.com/online-behaviour/machine-learning/actions/workflows/build.yml) |
|  Metadata consistency              | [![cffconvert](https://github.com/online-behaviour/machine-learning/actions/workflows/cffconvert.yml/badge.svg)](https://github.com/online-behaviour/machine-learning/actions/workflows/cffconvert.yml) |
| Lint                               | [![lint](https://github.com/online-behaviour/machine-learning/actions/workflows/lint.yml/badge.svg)](https://github.com/online-behaviour/machine-learning/actions/workflows/lint.yml) |
| Publish                            | [![publish](https://github.com/online-behaviour/machine-learning/actions/workflows/publish.yml/badge.svg)](https://github.com/online-behaviour/machine-learning/actions/workflows/publish.yml) |
| SonarCloud                         | [![sonarcloud](https://github.com/online-behaviour/machine-learning/actions/workflows/sonarcloud.yml/badge.svg)](https://github.com/online-behaviour/machine-learning/actions/workflows/sonarcloud.yml) |
| MarkDown link checker              | [![markdown-link-check](https://github.com/online-behaviour/machine-learning/actions/workflows/markdown-link-check.yml/badge.svg)](https://github.com/online-behaviour/machine-learning/actions/workflows/markdown-link-check.yml) |

## How to use machine_learning



The project setup is documented in [project_setup.md](project_setup.md). Feel free to remove this document (and/or the link to this document) if you don't need it.

## Installation

To install machine_learning from GitHub repository, do:

```console
git clone https://github.com/online-behaviour/machine-learning.git
cd machine-learning
python3 -m pip install .
```

## Documentation

Include a link to your project's full documentation here.

## Contributing

If you want to contribute to the development of machine-learning,
have a look at the [contribution guidelines](CONTRIBUTING.md).

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
reported by contacting the project team at e.tjongkimsang@esciencecenter.nl. All
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
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation);
1. you want to make a new release of the code base.

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/online-behaviour/machine-learning/issues) to see if someone already filed the same issue;
2. if your issue search did not yield any relevant results, make a new issue;
3. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/online-behaviour/machine-learning/issues) to see if someone already filed the same issue;
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
1. update the `CHANGELOG.md` file with change;
1. push your feature branch to (your fork of) the machine_learning repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
.. machine_learning documentation main file, created by
  sphinx-quickstart on Thu Jun 21 11:07:11 2018.
  You can adapt this file completely to your liking, but it should at least
  contain the root `toctree` directive.

Welcome to machine_learning's documentation!
==========================================================

.. toctree::
  :maxdepth: 2
  :caption: Contents:

API Reference
=============

.. toctree::
  :maxdepth: 2

  machine_learning <apidocs/machine_learning.rst>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
