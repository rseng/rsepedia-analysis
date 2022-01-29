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
cd notebooks
git init
git add --all
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/puregome/notebooks
```

Go to
[https://github.com/puregome?tab=repositories](https://github.com/puregome?tab=repositories)
and create a new repository named notebooks as an empty repository, then:

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
  5. go to [new code definition administration page](https://sonarcloud.io/project/new_code?id=puregome_notebooks) and select `Number of days` option
- The analysis will be run by [GitHub Action workflow](.github/workflows/sonarcloud.yml)
- To be able to run the analysis, a token must be created at [Sonarcloud account](https://sonarcloud.io/account/security/) and this token must be added as `SONAR_TOKEN` to [secrets on GitHub](https://github.com/puregome/notebooks/settings/secrets/actions)

## Package version number

- We recommend using [semantic versioning](https://guide.esciencecenter.nl/#/best_practices/releases?id=semantic-versioning).
- For convenience, the package version is stored in a single place: `notebooks/.bumpversion.cfg`.
  For updating the version number, make sure the dev dependencies are installed and run `bumpversion patch`,
  `bumpversion minor`, or `bumpversion major` as appropriate.
- Don't forget to update the version number before [making a release](https://guide.esciencecenter.nl/#/best_practices/releases)!

## Publish on Python Package Index (PyPI)

To publish your package on PyPI, you need to create a [PyPI API token](https://pypi.org/help#apitoken) and
save it as a secret called `PYPI_TOKEN` on [Settings page](https://github.com/puregome/notebooks/settings/secrets/actions)

[Creating a release](https://github.com/puregome/notebooks/releases/new) on GitHub will trigger a [GitHub action workflow](.github/workflows/publish.yml) to publish the release on PyPI for you.

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
# `notebooks` developer documentation

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
# install notebooks as an editable package
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

# recursively check import style for the notebooks module only
isort --recursive --check-only notebooks

# recursively check import style for the notebooks module only and show
# any proposed changes as a diff
isort --recursive --check-only --diff notebooks

# recursively fix import style for the notebooks module only
isort --recursive notebooks
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
cd $(mktemp -d --tmpdir notebooks.XXXXXX)

# fresh git clone ensures the release has the state of origin/main branch
git clone https://github.com/puregome/notebooks .

# prepare a clean virtual environment and activate it
python3 -m venv env
source env/bin/activate

# make sure to have a recent version of pip and setuptools
python3 -m pip install --upgrade pip setuptools

# install runtime dependencies and publishing dependencies
python3 -m pip install --no-cache-dir .
python3 -m pip install --no-cache-dir .[publishing]

# clean up any previously generated artefacts 
rm -rf notebooks.egg-info
rm -rf dist

# create the source distribution and the wheel
python3 setup.py sdist bdist_wheel

# upload to test pypi instance (requires credentials)
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

Visit
[https://test.pypi.org/project/notebooks](https://test.pypi.org/project/notebooks)
and verify that your package was uploaded successfully. Keep the terminal open, we'll need it later.

In a new terminal, without an activated virtual environment or an env directory:

```shell
cd $(mktemp -d --tmpdir notebooks-test.XXXXXX)

# prepare a clean virtual environment and activate it
python3 -m venv env
source env/bin/activate

# make sure to have a recent version of pip and setuptools
pip install --upgrade pip setuptools

# install from test pypi instance:
python3 -m pip -v install --no-cache-dir \
--index-url https://test.pypi.org/simple/ \
--extra-index-url https://pypi.org/simple notebooks
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
# PuReGoMe

[PuReGoMe](https://research-software.nl/projects/4728) is a research project of [Utrecht University](https://www.uu.nl/en/research/intelligent-software-systems/intelligent-systems) and the [Netherlands eScience Center](https://www.esciencecenter.nl/). We analyze Dutch social media messages to assess the opinion of the public towards the COVID-19 pandemic measures taken by the Dutch government.


## Data

PuReGoMe uses three data sources for social media analysis. Our main data source are Dutch tweets from [Twitter](https://twitter.com/). We use Dutch posts from [Reddit](https://www.reddit.com/) and comments from [Nu.nl](https://www.nu.nl) as sources of verification of the results obtained by the tweet analysis.

After a month has ended, the data from the month are collected. Here are the steps taken for each data source:


### Twitter

1. the tweets are taken from the crawler of [twiqs.nl](http://twiqs.nl). They are stored in json format in hourly files. (this data set is not publicly available)
2. the tweets are extracted from the json files and and stored in csv files (six columns: id\_str, in\_reply\_to\_status\_id\_str, user, verified, text, location). The conversion is performed by the script [query-text.py](https://github.com/puregome/queries/blob/master/query-text.py)
3. duplicate tweets are removed from from the csv files produced in step 2 by the script [text-unique.py](https://github.com/puregome/scripts/blob/master/text-unique.py)

(it would be useful to combine steps 2 and 3 in the future)


### Reddit

1. In a new month directory, run the script [get\_subreddit\_ids.py](https://github.com/puregome/scripts/blob/master/get_subreddit_ids.py) to automatically retrieve the subreddits of the Dutch corona reddits
2. Copy the list of ids of the subreddit [Megathread Coronavirus COVID-19 in Nederland](https://www.reddit.com/r/thenetherlands/search?q=Megathread+Coronavirus+COVID-19+in+Nederland&restrict_sr=on&sort=new&t=all) (submissions_ids_thenetherlands.txt) from a previous month directory and manually add the ids of recent subreddits
3. Run the script [coronamessagesnl.py](https://github.com/puregome/scripts/blob/master/coronamessagesnl.py) on the files `submissions_ids_*` to automatically retrieve the posts in the found subreddits
4. Run the notebook [reddit.ipynb](reddit.ipynb) to get all the posts from the monthly `downloads` directory and store them in the directory `text`


### Nu.nl

1. Run code blocks 1, 3 and 4 of the notebook selenium-test.ipynb, after updating the name of the file in URLFILE in code block 1
2. Restart the notebook and run code blocks 1 and 6, after changing the name of the new downloads directory in DATADIROUT in code block 6. This process takes many hours (even days) to complete
3. The notebook can be copied and several copies can be run in parallel
4. When the notebooks have finished: delete all pairs of files of sizes 1 and 3 in this month's directory (keep the ones with only size 1) and rerun the notebooks
5. Repeat step 4 until no articles with comments are found
6. go to (cd) the directory data/nunl
7. run the script [../../scripts/combineNunlComments.sh](https://github.com/puregome/scripts/blob/master/combineNunlComments.sh) to update the files in the main directory `downloads`
8. run the notebook [nunl-convert-data.ipynb](nunl-convert-data.ipynb) to regenerate the data files in the directory `text`
9. for fetching the article texts: run code blocks 1 and 2 of the notebook [selenium-test.ipynb](selenium-test.ipynb), after updating the variables `URLFILE` and `OUTFILEMETADONELIST`


## Analysis

PuReGoMe performs analysis on three different levels: by counting messages, by determining their polarity (sentiment) and by determining their stance with respect to anti-pandemic government measures.


### Frequency analysis

Frequency analysis of tweets is performed in the notebook [tweet-counts.ipynb](tweet-counts.ipynb). The notebook defines several pandemic queries, for example face mask, lockdown, social distancing and pandemic, where pandemic is a combination of 60+ relevant terms. The notebook produces a graph with the absolute daily frequencies of the tweets matching each of these pandemic queries:

![tweet frequencies](tweet-frequencies.png)

Frequency analysis of the Nu.nl and Reddit data is included in the respective data generation notebooks [nunl-convert-data.ipynb](nunl-convert-data.ipynb) and [reddit.ipynb](reddit.ipynb)


### Polarity analysis

Polarity analysis is the same as sentiment analysis. This analysis is performed by the notebook [sentiment-pattern.ipynb](sentiment-pattern.ipynb) which uses the Python package [Pattern](https://github.com/clips/pattern) for sentiment analysis of Dutch text (De Smedt &amp; Daelemans, 2011). The notebook requires two types of input files: the csv files in the text directories of each of the data sources and sentiment score files which should be generated from these csv files with the script [../../scripts/sentiment-pattern-text.py](https://github.com/puregome/scripts/blob/master/sentiment-pattern-text.py) The polarity analysis of the different topics takes a lot of time and can be run in parallel. It produces time series graphs for all tweets, all pandemic tweets and several individual pandemic topics.

![polarity for face masks, social distancing and lockdown over time](sentiment-all.png)


### Stance analysis

Stance analysis is performed by the notebook [fasttext.ipynb](fasttext.ipynb). The analysis originates from a model trained by [fastText](https://github.com/facebookresearch/fastText) on manually labeled tweets. The notebook contains a section for searching for the best parameters of fastText using grid search but when the training data is unchanged this section can be skipped. The notebook has two main modes related to topics: analysis related to the social distancing policy and analysis related to the former (April 2020) face mask policy. These are the only two topics for which we have manually labeled training data. The graphs combine analysis for all three data sources used in the project: Twitter, Nu.nl and Reddit.

![stance for social distancing](social-distancing-all.png)


### Other analyses

The notebook [topic-analysis.ipynb](topic-analysis.ipynb) is used to find new topics in a day of tweets. It compares the vocabulary of a day of tweets with the vocabulary of the preceding day.

[echo-chambers.ipynb](echo-chambers.ipynb) is used for finding groups of users which collectively retweet similar content. The notebook found a group of a few hundred users retweeting right-wing propaganda. Further study needs to be done to check if this content has any effect on the findings of this project.

[geo-analysis.ipynb](geo-analysis.ipynb) and [geo-classification.ipynb](geo-classification.ipynb) can be used to divide the tweets in groups depending on on the location of the tweeter. This only works for about half of the data set. Next, maps representing tweet data can be created with the notebook [maps.ipynb](maps.ipynb).

[corona-nl-totals.ipynb](corona-nl-totals.ipynb) creates graphs of the number infections, hospitalizations and deaths in The Netherlands based on data provided by the health organization [RIVM](https://www.rivm.nl).


## Publications, talks and media coverage

Erik Tjong Kim Sang, Shihan Wang, Marijn Schraagen and Mehdi Dastani, [**Extracting Stances on Pandemic Measures from Social Media Data**](https://ifarm.nl/erikt/papers/ieee2021.pdf). 17th IEEE eScience Conference (poster), 2021.

Erik Tjong Kim Sang, Marijn Schraagen, Shihan Wang and  Mehdi Dastani, [**Transfer Learning for Stance Analysis in COVID-19 Tweets**](https://ifarm.nl/erikt/papers/clin-20210419.pdf). CLIN 2021. ([data annotations](http://145.100.59.103/downloads/clin2021-annotations.zip))

Erik Tjong Kim Sang, Marijn Schraagen, Mehdi Dastani and Shihan Wang, [**Discovering Pandemic Topics on Twitter**](https://ifarm.nl/erikt/papers/dhbenelux-20210509.pdf). DHBenelux 2021.

Erik Tjong Kim Sang [**PuReGoMe: Social Media Analysis of the Pandemic**](https://ifarm.nl/erikt/talks/20210211-escience.pdf). Lunch talk, Netherlands eScience Center, Amsterdam, The Netherlands, 11 February 2021.

Shihan Wang, Marijn Schraagen, Erik Tjong Kim Sang and Mehdi Dastani, [**Dutch General Public Reaction on Governmental COVID-19 Measures and Announcements in Twitter Data**](https://arxiv.org/abs/2006.07283). Preprint report on arXiv.org, 21 December 2020.

Shihan Wang, Marijn Schraagen, Erik Tjong Kim Sang and Mehdi Dastani, [**Public Sentiment on Governmental COVID-19 Measures in Dutch Social Media**](https://www.aclweb.org/anthology/2020.nlpcovid19-2.17/). In: Workshop on NLP for COVID-19 (Part 2) at EMNLP 2020
(NLP-COVID19-EMNLP), 20 November 2020.

Erik Tjong Kim Sang [**PuReGoMe: Dutch Public Reaction on Governmental COVID-19 Measures and Announcements**](https://ifarm.nl/erikt/talks/20200626-escience.pdf). Lunch talk, Netherlands eScience Center, Amsterdam, The Netherlands, 26 June 2020.

Shihan Wang, **Public Sentiment during COVID-19 -Data Mining on Twitter data**. Talk at the [CLARIN Café](https://www.clarin.eu/content/clarin-cafe), Utrecht, The Netherlands, 27 May 2020.

Redactie Emerce, [**Onderzoekers leiden publieke opinie coronamaatregelen af uit social media-data**](https://www.emerce.nl/nieuws/onderzoekers-leiden-publieke-opinie-coronamaatregelen-af-uit-social-mediadata). Emerce, 12 May 2020 (in Dutch).

Nienke Vergunst, [**Researchers use social media data to analyse public sentiment about Coronavirus measures**](https://www.uu.nl/en/news/researchers-use-social-media-data-to-analyse-public-sentiment-about-coronavirus-measures). University of Utrecht news message, 11 May 2020.

# Information added by the Python template

## Badges

| fair-software.eu recommendations | |
| :-- | :--  |
| (1/5) code repository              | [![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/puregome/notebooks) |
| (2/5) license                      | [![github license badge](https://img.shields.io/github/license/puregome/notebooks)](https://github.com/puregome/notebooks) |
| (3/5) community registry           | [![Research Software Directory](https://img.shields.io/badge/rsd-Research%20Software%20Directory-00a3e3.svg)](https://www.research-software.nl/software/puregome) |
| (4/5) citation                     | [![DOI](https://zenodo.org/badge/256566641.svg)](https://zenodo.org/badge/latestdoi/256566641) |
| (5/5) checklist                    | [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/4357/badge)](https://bestpractices.coreinfrastructure.org/projects/4357) |
| howfairis                            | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu) |
| **Other best practices**           | &nbsp; |
| Static analysis              | [![workflow scq badge](https://sonarcloud.io/api/project_badges/measure?project=puregome_notebooks&metric=alert_status)](https://sonarcloud.io/dashboard?id=puregome_notebooks) |
| Coverage              | [![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=puregome_notebooks&metric=coverage)](https://sonarcloud.io/dashboard?id=puregome_notebooks) |
| **GitHub Actions**                 | &nbsp; |
| Build                              | [![build](https://github.com/puregome/notebooks/actions/workflows/build.yml/badge.svg)](https://github.com/puregome/notebooks/actions/workflows/build.yml) |
|  Metadata consistency              | [![cffconvert](https://github.com/puregome/notebooks/actions/workflows/cffconvert.yml/badge.svg)](https://github.com/puregome/notebooks/actions/workflows/cffconvert.yml) |
| Lint                               | [![lint](https://github.com/puregome/notebooks/actions/workflows/lint.yml/badge.svg)](https://github.com/puregome/notebooks/actions/workflows/lint.yml) |
| Publish                            | [![publish](https://github.com/puregome/notebooks/actions/workflows/publish.yml/badge.svg)](https://github.com/puregome/notebooks/actions/workflows/publish.yml) |
| SonarCloud                         | [![sonarcloud](https://github.com/puregome/notebooks/actions/workflows/sonarcloud.yml/badge.svg)](https://github.com/puregome/notebooks/actions/workflows/sonarcloud.yml) |
| MarkDown link checker              | [![markdown-link-check](https://github.com/puregome/notebooks/actions/workflows/markdown-link-check.yml/badge.svg)](https://github.com/puregome/notebooks/actions/workflows/markdown-link-check.yml) |

## How to use notebooks

The project setup is documented in [project_setup.md](project_setup.md). Feel free to remove this document (and/or the link to this document) if you don't need it.

## Installation

To install notebooks from GitHub repository, do:

```console
git clone https://github.com/puregome/notebooks.git
cd notebooks
python3 -m pip install .
```

## Documentation

Include a link to your project's full documentation here.

## Contributing

If you want to contribute to the development of notebooks,
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

1. use the search functionality [here](https://github.com/puregome/notebooks/issues) to see if someone already filed the same issue;
2. if your issue search did not yield any relevant results, make a new issue;
3. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/puregome/notebooks/issues) to see if someone already filed the same issue;
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
1. push your feature branch to (your fork of) the notebooks repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
