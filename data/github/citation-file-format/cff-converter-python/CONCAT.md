# 2.0.0

## CLI

- added APA output (PR [#149](https://github.com/citation-file-format/cff-converter-python/pull/149); thanks [@wleoncio](https://github.com/wleoncio))
- added support for validation and conversion of `CITATION.cff` files with `cff-version: 1.2.0`
- argument `--outputformat` was renamed to `--format`
- argument `-ig`, `--ignore-suspect-keys` was removed
- argument `--verbose` was removed
- argument `--show-trace` was added
 
## Library

- added APA output (PR [#149](https://github.com/citation-file-format/cff-converter-python/pull/149); thanks [@wleoncio](https://github.com/wleoncio))
- added support for validation and conversion of `CITATION.cff` files with `cff-version: 1.2.0`
- simplified the `Citation` class and its interface
- `cli` is no longer part of the public interface of the library
- URLs are now constructed from `identifiers`, `repository`, `repository-artifact`, `repository-code`, or `url`, with a transparent mechanism to choose what to use given the data that is available from a given `CITATION.cff` file
- Authors are now constructed from `given-names`, `family-names` (including `name-particle` and `name-suffix`), `alias`, `name`, `affiliation` and `orcid`, with a transparent mechanism to choose what to use given the data that is available from a given `CITATION.cff` file

## Other

- switched to static configuration (setup.cfg over setup.py)
- dependencies are now in `setup.cfg` as opposed to `requirements[-dev].txt`
- updated version ranges for dependencies
- tests are no longer `unittest.TestCase` based, but pytest with fixtures
- added jsonschema based validation for CITATION.cff files with `cff-version: 1.2.0`
- implemented _State pattern_ for `Citation` to help it deal with multiple behaviors under past and future versions of the Citation File Format.
- switched from TravisCI to GitHub Actions workflows, added linting and publishing workflows
- CI is now testing against Python 3.6, 3.7, 3.8, and 3.9 on Mac, Linux and Windows
- copies of the relevant schemas are now bundled with the package 
- organized the tests to be more orthogonal to each other / less overlap between tests

# 1.3.3

-   With recent changes to the release process, the schema will be in a
    different place than before. This release fixes
    <https://github.com/citation-file-format/cff-converter-python/issues/119>).

# 1.3.2

-   the ruamel.yaml dependency was not specified tightly enough,
    `requirements.txt` has been updated as have the notes for
    maintainers.

# 1.3.1

-   'cff-version: 1.0.3' is now interpreted as 1.0.3-1 (the latest
    schema version that implements the spec 1.0.3). This will fix some
    problems with the list of SPDX license abbreviations. These
    additional licenses should now work:
    -   `AGPL-3.0-only`
    -   `AGPL-3.0-or-later`
    -   `BSD-1-Clause`
    -   `BSD-2-Clause-Patent`
    -   `CDLA-Permissive-1.0`
    -   `CDLA-Sharing-1.0`
    -   `EPL-2.0`
    -   `EUPL-1.2`
    -   `GFDL-1.1-only`
    -   `GFDL-1.1-or-later`
    -   `GFDL-1.2-only`
    -   `GFDL-1.2-or-later`
    -   `GFDL-1.3-only`
    -   `GFDL-1.3-or-later`
    -   `GPL-1.0-only`
    -   `GPL-1.0-or-later`
    -   `GPL-2.0-only`
    -   `GPL-2.0-or-later`
    -   `GPL-3.0-only`
    -   `GPL-3.0-or-later`
    -   `LGPL-2.0-only`
    -   `LGPL-2.0-or-later`
    -   `LGPL-2.1-only`
    -   `LGPL-2.1-or-later`
    -   `LGPL-3.0-only`
    -   `LGPL-3.0-or-later`

# 1.3.0

-   added schema.org converter method

# 1.2.2

-   added documentation for the Google Cloud Function interface

# 1.2.1

-   setup.py no longer includes test dependencies as install
    dependencies

# 1.2.0

-   corrected an error where cffconvert could not raise an error during
    validation
    (<https://github.com/citation-file-format/cff-converter-python/issues/94>).

# 1.1.0

-   replaced pykwalifire with its parent pykwalify
-   now works for python 3.7 (refs \#80)
-   not using PyYAML anymore (but it still comes along with pykwalify
    for some reason)
-   added a function that can be used as Google Cloud function
-   hopefully fixed parsing of strings that should have been entered as
    dates (the new validator does not find that offensive, hence I had
    to fix it myself)

# 1.0.4

-   replaced PyYAML dependency with ruamel.yaml

# 1.0.3

-   security bugfix by updating requests from 2.18.4 to 2.20.0

# 1.0.2

-   fixed bug
    <https://github.com/citation-file-format/cff-converter-python/issues/82>
    (warnings on stdout)

# 1.0.1

-   fixed bug
    <https://github.com/citation-file-format/cff-converter-python/issues/73>
    (orcid format in zenodo export)

# 1.0.0

-   first stable release
-   solved bug
    <https://github.com/citation-file-format/cff-converter-python/issues/59>
    (cffconvert creates local file `data.yaml` and `schema.yaml` on
    validate)

# 0.0.5

-   Minor changes

# 0.0.4

-   added optional validation of CITATION.cff files using pykwalifire
    (`--validate`)
-   added printing the CITATION.cff contents from the command line
-   added unit tests for command line interface
-   added integration with sonarcloud code quality monitoring
-   removed shorthand command line argument `-v` (represented both
    `--validate` and --verbose)
-   added showing its own version (`--version`)
-   command line argument `--ignore-suspect-keys` no longer needs to be
    assigned a value, it's simply a flag

# Documentation for developers

## Install

```shell
# get a copy of the cff-converter-python software
git clone https://github.com/citation-file-format/cff-converter-python.git
# change directory into cff-converter-python
cd cff-converter-python
# make a virtual environment named env
python3 -m venv env
# activate the virtual environment
source env/bin/activate
# upgrade pip, wheel, setuptools
python3 -m pip install --upgrade pip wheel setuptools
# install cffconvert and the 'dev' set of additional dependencies
python3 -m pip install --editable .[dev]
```

## Running tests

Running the tests requires an activated virtual environment with the development tools installed.

```shell
# (from the project root)

# run all tests
python3 -m pytest test/

# tests for consistent file naming
bash test/test_consistent_file_naming.sh dir=test/
bash test/test_consistent_file_naming.sh dir=livetest/

# tests for consistent versioning
python3 -m pytest test/test_consistent_versioning.py
```

## Running linters locally

For linting we use [prospector](https://pypi.org/project/prospector/) and to sort imports we will use
[isort](https://pycqa.github.io/isort/). Running the linters requires an activated virtual environment with the
development tools installed.

```shell
# linter
prospector

# recursively check import style for the cffconvert module only
isort --check-only cffconvert

# recursively check import style for the cffconvert module only and show
# any proposed changes as a diff
isort --check-only --diff cffconvert

# recursively fix import style for the cffconvert module only
isort cffconvert
```

Developers should consider enabling automatic linting with `prospector` and `isort` on commit by enabling the git hook from `.githooks/pre-commit`, like so:

```shell
git config --local core.hooksPath .githooks
```


## For maintainers

### Making a release


1. make sure the release notes are up to date
1. preparation

    ```shell
    # remove old cffconvert from your system if you have it
    python3 -m pip uninstall cffconvert

    # this next command should now return empty
    which cffconvert

    # install the package to user space, using no caching (can bring to light dependency problems)
    python3 -m pip install --user --no-cache-dir .
    # check if cffconvert works, e.g.
    cffconvert --version
    
    # run the tests, make sure they pass
    python3 -m pip pytest test

    # git push everything, merge into main as appropriate
    ```
    
1. publishing on test instance of PyPI

    ```shell
    # verify that everything has been pushed and merged by testing as follows
    cd $(mktemp -d --tmpdir cffconvert-release.XXXXXX)
    git clone https://github.com/citation-file-format/cff-converter-python.git .
    python3 -m venv env
    source env/bin/activate
    python3 -m pip install --upgrade pip wheel setuptools
    python3 -m pip install --no-cache-dir .

    # register with PyPI test instance https://test.pypi.org

    # remove these directories if you have them
    rm -rf dist
    rm -rf cffconvert-egg.info
    # make a source distribution:
    python setup.py sdist
    # make a wheel
    python setup.py bdist_wheel 
    # install the 'upload to pypi/testpypi tool' aka twine
    pip install .[publishing]
    # upload the contents of the source distribution we just made (requires credentials for test.pypi.org)
    twine upload --repository-url https://test.pypi.org/legacy/ dist/*
    ```
    
1. Checking the package

    Open another shell but keep the other one. We'll return to the first shell momentarily.
    
    Verify that there is a new version of the package on Test PyPI https://test.pypi.org/project/cffconvert/

    ```shell
    python3 -m pip -v install --user --no-cache-dir \
    --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple cffconvert

    # check that the package works as it should when installed from pypitest
    ```
1. FINAL STEP: upload to PyPI

    Go back to the first shell, then (requires credentials for pypi.org)

    ```shell
    twine upload dist/*
    ```
1. Make the release on GitHub
1. Go to Zenodo, log in to inspect the draft. Then click `Publish` to finalize it.

### Building the docker image

```shell
# (requires 2.0.0 to be downloadable from PyPI)
docker build --tag cffconvert:2.0.0 .
docker build --tag cffconvert:latest .
```

See if the Docker image works as expected:
```shell
docker run --rm -v $PWD:/app cffconvert --validate
docker run --rm -v $PWD:/app cffconvert --version
docker run --rm -v $PWD:/app cffconvert
# etc
```

### Publishing on DockerHub

See <https://docs.docker.com/docker-hub/repos/#pushing-a-docker-container-image-to-docker-hub> for more information on publishing.

```shell
# log out of any dockerhub credentials
docker logout

# log back in with username 'citationcff' credentials
docker login

# re-tag existing images
docker tag cffconvert:2.0.0 citationcff/cffconvert:2.0.0
docker tag cffconvert:latest citationcff/cffconvert:latest

# publish
docker push citationcff/cffconvert:2.0.0
docker push citationcff/cffconvert:latest
```
# `cffconvert`

[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1162057.svg)](https://doi.org/10.5281/zenodo.1162057)
[![testing](https://github.com/citation-file-format/cff-converter-python/actions/workflows/testing.yml/badge.svg)](https://github.com/citation-file-format/cff-converter-python/actions/workflows/testing.yml)
[![linting](https://github.com/citation-file-format/cff-converter-python/actions/workflows/linting.yml/badge.svg)](https://github.com/citation-file-format/cff-converter-python/actions/workflows/linting.yml)
[![Code Smells](https://sonarcloud.io/api/project_badges/measure?project=citation-file-format_cff-converter-python&metric=code_smells)](https://sonarcloud.io/dashboard?id=citation-file-format_cff-converter-python)
[![PyPI Badge](https://img.shields.io/pypi/v/cffconvert.svg?colorB=blue)](https://pypi.python.org/pypi/cffconvert/)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/1811/badge)](https://bestpractices.coreinfrastructure.org/projects/1811)
[![Research Software Directory](https://img.shields.io/badge/rsd-cffconvert-00a3e3.svg)](https://www.research-software.nl/software/cff-converter-python)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![Docker Pulls](https://img.shields.io/docker/pulls/citationcff/cffconvert)](https://hub.docker.com/r/citationcff/cffconvert)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/cffconvert)](https://pypistats.org/packages/cffconvert)

Command line program to validate and convert [`CITATION.cff`](https://github.com/citation-file-format/citation-file-format) files.

## Supported input versions of the Citation File Format

| Citation File Format schema version | Link to Zenodo release |
| --- | --- |
| `1.2.0` | [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5171937.svg)](https://doi.org/10.5281/zenodo.5171937) |
| `1.1.0` | [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4813122.svg)](https://doi.org/10.5281/zenodo.4813122) |
| `1.0.3` | [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1222163.svg)](https://doi.org/10.5281/zenodo.1222163) |
| `1.0.2` | [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1120256.svg)](https://doi.org/10.5281/zenodo.1120256) |
| `1.0.1` | [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1117789.svg)](https://doi.org/10.5281/zenodo.1117789) |

## Supported output formats

1.  APA-like plaintext
1.  BibTeX
3.  CodeMeta
4.  EndNote
1.  RIS
1.  schema.org JSON
1.  Zenodo JSON

`cffconvert` does not support converting items from `references` or `preferred-citation` keys at the moment.

## Installing

To install in user space, 

```shell
python3 -m pip install --user cffconvert
```
Ensure that the user space directory `~/.local/bin/` is on the `PATH`.

```shell
which cffconvert
```
should now return the location of the program.

See [docs/alternative-install-options.md](docs/alternative-install-options.md) for alternative install options.

## Docker

`cffconvert` is available from DockerHub: https://hub.docker.com/r/citationcff/cffconvert

Example usage:

```shell
docker run --rm -v $PWD:/app citationcff/cffconvert --validate
docker run --rm -v $PWD:/app citationcff/cffconvert --version
docker run --rm -v $PWD:/app citationcff/cffconvert --help
# etc
```

## Command line interface

See `cffconvert`'s options:

```shell
cffconvert --help
```

Shows:

```shell
Usage: cffconvert [OPTIONS]

Options:
  -i, --infile PATH               Path to the CITATION.cff input file. If this
                                  option is omitted, './CITATION.cff' is used.
  -o, --outfile PATH              Path to the output file.
  -f, --format [apalike|bibtex|cff|codemeta|endnote|ris|schema.org|zenodo]
                                  Output format.
  -u, --url TEXT                  URL to the CITATION.cff input file.
  -h, --help                      Show help and exit.
  --show-trace                    Show error trace.
  --validate                      Validate the CITATION.cff file and exit.
  --version                       Print version and exit.
```

## Example usage

### Validating a local CITATION.cff file

```shell
cffconvert --validate
cffconvert --validate -i CITATION.cff
cffconvert --validate -i ${PWD}/CITATION.cff
cffconvert --validate -i ../some-other-dir/CITATION.cff
```

### Validating a remote CITATION.cff file

```shell
cffconvert --validate --url https://github.com/<org>/<repo>
cffconvert --validate --url https://github.com/<org>/<repo>/commit/<sha>
cffconvert --validate --url https://github.com/<org>/<repo>/tree/<sha>
cffconvert --validate --url https://github.com/<org>/<repo>/tree/<tag>
cffconvert --validate --url https://github.com/<org>/<repo>/tree/<branch>
```


### Converting metadata to other formats

If there is a valid `CITATION.cff` file in the current directory, you can convert to various other formats and 
print the result on standard out with:

```shell
cffconvert -f bibtex
cffconvert -f codemeta
cffconvert -f endnote
cffconvert -f ris
cffconvert -f schema.org
cffconvert -f zenodo
cffconvert -f apalike
```

### Writing to a file

```shell
# with i/o redirection:
cffconvert -f bibtex > bibtex.bib
cffconvert -f zenodo > zenodo.json
cffconvert -f endnote > ${PWD}/endnote.enw
# etc

# without i/o redirection
cffconvert -f bibtex -o bibtex.bib
cffconvert -f zenodo -o zenodo.json
cffconvert -f endnote -o ${PWD}/endnote.enw
# etc
```
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or
question to a full fledged [pull
request](https://help.github.com/articles/about-pull-requests/). However, we ask
that you read and follow this organization's [Code of Conduct](https://github.com/citation-file-format/citation-file-format/blob/master/CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add
   a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality
   [here](https://github.com/citation-file-format/cff-converter-python/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality
   [here](https://github.com/citation-file-format/cff-converter-python/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue,
   making sure to provide enough information to the rest of the community to
   understand the cause and context of the problem. Depending on the issue, you may
   want to include:
    - the [SHA
      hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas)
      of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're
      using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you
   start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea
   being a good idea;
1. if needed, fork the repository to your own Github profile and create your own
   feature branch off of the latest master commit. While working on your feature
   branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``pytest test/`` and ``pytest livetest/``;
1. add your own tests (if applicable);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your
   fork of) the ``cff-converter-python`` repository on GitHub;
1. create the pull request, e.g. following the instructions
   [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you have a valuable contribution to make, but you don't know
how to write or run tests for it, or how to create the documentation: don't
let this discourage you from making the pull request; we can help you! Just go
ahead and submit the pull request, but keep in mind that you might be asked to
append additional commits to your pull request.

The subdirectories here contain various combinations of author properties, mirroring the setup in https://github.com/citation-file-format/cff-converter-python/blob/3265d14b2db2c35b5670247c653556f9e6286093/cffconvert/behavior_shared/schemaorg_author.py 

- `authors[i].given-names`: omitted, valid
- `authors[i].name-particle`: omitted, valid
- `authors[i].family-names`: omitted, valid
- `authors[i].name-suffix`: omitted, valid
- `authors[i].name`: omitted, valid
- `authors[i].alias`: omitted, valid
- `authors[i].affiliation`: omitted, valid
- `authors[i].orcid`: omitted, valid

There are tests for a singular author (`./one`) and tests for when there are two authors (`./two`).

The table below lists the naming convention based on the state of an author in `CITATION.cff`:

| subdirectory | has_given_name | has_family_name | has_alias | has_name | has_affiliation | has_orcid | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `GFANAO` | True | True | True | True | True | True | |
| `GFANA_` | True | True | True | True | True | False | |
| `GFAN_O` | True | True | True | True | False | True | |
| `GFAN__` | True | True | True | True | False | False | |
| `GFA_AO` | True | True | True | False | True | True | |
| `GFA_A_` | True | True | True | False | True | False | |
| `GFA__O` | True | True | True | False | False | True | |
| `GFA___` | True | True | True | False | False | False | |
| `GF_NAO` | True | True | False | True | True | True | |
| `GF_NA_` | True | True | False | True | True | False | |
| `GF_N_O` | True | True | False | True | False | True | |
| `GF_N__` | True | True | False | True | False | False | |
| `GF__AO` | True | True | False | False | True | True | |
| `GF__A_` | True | True | False | False | True | False | |
| `GF___O` | True | True | False | False | False | True | |
| `GF____` | True | True | False | False | False | False | |
| `G_ANAO` | True | False | True | True | True | True | |
| `G_ANA_` | True | False | True | True | True | False | |
| `G_AN_O` | True | False | True | True | False | True | |
| `G_AN__` | True | False | True | True | False | False | |
| `G_A_AO` | True | False | True | False | True | True | |
| `G_A_A_` | True | False | True | False | True | False | |
| `G_A__O` | True | False | True | False | False | True | |
| `G_A___` | True | False | True | False | False | False | |
| `G__NAO` | True | False | False | True | True | True | |
| `G__NA_` | True | False | False | True | True | False | |
| `G__N_O` | True | False | False | True | False | True | |
| `G__N__` | True | False | False | True | False | False | |
| `G___AO` | True | False | False | False | True | True | |
| `G___A_` | True | False | False | False | True | False | |
| `G____O` | True | False | False | False | False | True | |
| `G_____` | True | False | False | False | False | False | |
| `_FANAO` | False | True | True | True | True | True | |
| `_FANA_` | False | True | True | True | True | False | |
| `_FAN_O` | False | True | True | True | False | True | |
| `_FAN__` | False | True | True | True | False | False | |
| `_FA_AO` | False | True | True | False | True | True | |
| `_FA_A_` | False | True | True | False | True | False | |
| `_FA__O` | False | True | True | False | False | True | |
| `_FA___` | False | True | True | False | False | False | |
| `_F_NAO` | False | True | False | True | True | True | |
| `_F_NA_` | False | True | False | True | True | False | |
| `_F_N_O` | False | True | False | True | False | True | |
| `_F_N__` | False | True | False | True | False | False | |
| `_F__AO` | False | True | False | False | True | True | |
| `_F__A_` | False | True | False | False | True | False | |
| `_F___O` | False | True | False | False | False | True | |
| `_F____` | False | True | False | False | False | False | |
| `__ANAO` | False | False | True | True | True | True | |
| `__ANA_` | False | False | True | True | True | False | |
| `__AN_O` | False | False | True | True | False | True | |
| `__AN__` | False | False | True | True | False | False | |
| `__A_AO` | False | False | True | False | True | True | |
| `__A_A_` | False | False | True | False | True | False | |
| `__A__O` | False | False | True | False | False | True | |
| `__A___` | False | False | True | False | False | False | |
| `___NAO` | False | False | False | True | True | True | |
| `___NA_` | False | False | False | True | True | False | |
| `___N_O` | False | False | False | True | False | True | |
| `___N__` | False | False | False | True | False | False | |
| `____AO` | False | False | False | False | True | True | |
| `____A_` | False | False | False | False | True | False | |
| `_____O` | False | False | False | False | False | True | |
| `______` | False | False | False | False | False | False | |
The subdirectories here contain various combinations of sources for a doi 

- `doi`: omitted, valid
- `identifiers[i].type==doi`: omitted, valid

| subdirectory | has_doi   | has_identifiers_doi | notes |
| --- | --- | --- | --- |
| `__` | False | False | |
| `_I` | False | True | |
| `D_` | True | False | |
| `DI` | True | True | |
The subdirectories here contain various combinations of sources for a URL

- `identifiers[i].type==url`: omitted, valid
- `repository`: omitted, valid
- `repository-artifact`: omitted, valid
- `repository-code`: omitted, valid
- `url`: omitted, valid


| subdirectory | has_indentifiers_url | has_repository | has_repository_artifact | has_repository_code |  has_url | notes |
| --- | --- | --- | --- | --- | --- | --- |
| `IRACU` | True | True | True |  True |  True | |
| `IRAC_` | True | True | True |  True |  False | |
| `IRA_U` | True | True | True |  False |  True | |
| `IRA__` | True | True | True |  False |  False | |
| `IR_CU` | True | True | False |  True |  True | |
| `IR_C_` | True | True | False |  True |  False | |
| `IR__U` | True | True | False |  False |  True | |
| `IR___` | True | True | False |  False |  False | |
| `I_ACU` | True | False | True |  True |  True | |
| `I_AC_` | True | False | True |  True |  False | |
| `I_A_U` | True | False | True |  False |  True | |
| `I_A__` | True | False | True |  False |  False | |
| `I__CU` | True | False | False |  True |  True | |
| `I__C_` | True | False | False |  True |  False | |
| `I___U` | True | False | False |  False |  True | |
| `I____` | True | False | False |  False |  False | |
| `_RACU` | False | True | True |  True |  True | |
| `_RAC_` | False | True | True |  True |  False | |
| `_RA_U` | False | True | True |  False |  True | |
| `_RA__` | False | True | True |  False |  False | |
| `_R_CU` | False | True | False |  True |  True | |
| `_R_C_` | False | True | False |  True |  False | |
| `_R__U` | False | True | False |  False |  True | |
| `_R___` | False | True | False |  False |  False | |
| `__ACU` | False | False | True |  True |  True | |
| `__AC_` | False | False | True |  True |  False | |
| `__A_U` | False | False | True |  False |  True | |
| `__A__` | False | False | True |  False |  False | |
| `___CU` | False | False | False |  True |  True | |
| `___C_` | False | False | False |  True |  False | |
| `____U` | False | False | False |  False |  True | |
| `_____` | False | False | False |  False |  False | |
# Alternative install options for using `cffconvert`

## Install in virtual environment

```shell
# use venv to make a virtual environment named env
python3 -m venv env

# activate the environment
source env/bin/activate

# install cffconvert in it
pip install cffconvert
```

## Install globally

Note: this option needs sudo rights.

```shell
sudo -H python3 -m pip install cffconvert
```

## Install with conda

Make an environment definition file `environment.yml` with the following contents:

```yaml
name: env
channels:
  - conda-forge
  - defaults
dependencies:
  - pip
  - pip:
    - cffconvert
```

Then run:

```shell
conda env create --file environment.yml
conda activate env
```

## No-install options

### Using `cffconvert` as a Google Cloud Function

`cffconvert` comes with [an interface](/cffconvert/gcloud.py) for
running as a Google Cloud Function. We set it up here
<https://bit.ly/cffconvert> for the time being / as long as we have
enough credits on the Google Cloud Function platform.

Really, all the Google Cloud interface does is get any supplied URL
parameters, and use them as if they had been entered as command line
arguments. For more detailed explanation and examples, see
<https://bit.ly/cffconvert>.

On Google Cloud Function, set `requirements.txt` to:

```text
cffconvert[gcloud]
```

and use the following as `main.py`:

```python
from cffconvert.gcloud.gcloud import cffconvert

def main(request):
   return cffconvert(request)
```

### Docker

Build the Docker container 

```shell
cd <project root>
docker build --tag cffconvert:2.0.0 .
docker build --tag cffconvert:latest .
```

Build the Docker container 

```shell
cd <where your CITATION.cff is>
docker run --rm -ti -v ${PWD}:/app cffconvert
```

# BibTeX example result

```bibtex
@misc{YourReferenceHere,
author = {Jurriaan H. Spaaks},
doi = {10.5281/zenodo.1162057},
month = {11},
title = {cffconvert},
url = {https://github.com/citation-file-format/cff-converter-python},
year = {2019}
}
```
# CodeMeta example result

```jsonld
{
   "@context": "https://doi.org/10.5063/schema/codemeta-2.0", 
   "@type": "SoftwareSourceCode", 
   "author": [
      {
         "@type": "Person", 
         "familyName": "Spaaks", 
         "givenName": "Jurriaan H."
      }
   ], 
   "codeRepository": "https://github.com/citation-file-format/cff-converter-python", 
   "datePublished": "2019-11-12", 
   "identifier": "https://doi.org/10.5281/zenodo.1162057", 
   "name": "cffconvert", 
   "version": "1.3.3"
}
```
# Zenodo JSON example result

```json
{
   "creators": [
      {
         "name": "Spaaks, Jurriaan H."
      }
   ], 
   "doi": "10.5281/zenodo.1162057", 
   "publication_date": "2019-11-12", 
   "title": "cffconvert", 
   "version": "1.3.3"
}
```
# schema.org example result

```jsonld
{
   "@context": "https://schema.org", 
   "@type": "SoftwareSourceCode", 
   "author": [
      {
         "@type": "Person", 
         "familyName": "Spaaks", 
         "givenName": "Jurriaan H."
      }
   ], 
   "codeRepository": "https://github.com/citation-file-format/cff-converter-python", 
   "datePublished": "2019-11-12", 
   "identifier": "https://doi.org/10.5281/zenodo.1162057", 
   "name": "cffconvert", 
   "version": "1.3.3"
}
```
# APA-like example result

Spaaks J.H. (2019). cffconvert (version 1.3.3). DOI: http://doi.org/10.5281/zenodo.1162057 URL: https://github.com/citation-file-format/cff-converter-python
# RIS example result

```text
TY  - GEN
AU  - Spaaks, Jurriaan H.
DA  - 2019-11-12
DO  - 10.5281/zenodo.1162057
PY  - 2019
TI  - cffconvert
UR  - https://github.com/citation-file-format/cff-converter-python
ER
```
# EndNote example result

```text
%0 Generic
%A Spaaks, Jurriaan H.
%D 2019
%R 10.5281/zenodo.1162057
%T cffconvert
%U https://github.com/citation-file-format/cff-converter-python
%9 source code
```
