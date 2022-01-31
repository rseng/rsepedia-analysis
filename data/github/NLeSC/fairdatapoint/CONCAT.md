# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.7.2] - 2021-02-01
### Changed
- Changed the url of Swagger UI from `/ui` to the base url `/`

## [0.7.1] - 2020-09-15
### Added
- PUT method used to update metadata

### Changed
- Replaced `flask-restplus` with `connexion` package for OpenAPI generation
- Refactored code based on OpenAPI specification
- Removed ending '/' in endpoints, e.g. '/catalog/' was changed to '/catalog'.
- Removed Bottle and Paste dependecies
- Removed docker compose file for production run


## [0.7.0] - 2020-06-12
### Added
- POST and DELETE methods
- SHACL validator to validate different layers of metadata
- Support for SPARQL database, such as Virtuoso RDF Triple Store
- Docker compose file
- Running in production mode

### Changed
- Improved GET and POST methods
- Updated support for common RDF serializations
- Removed loading metadata for the command tool `fdp-run`
- Removed support for INI-based configuration files
- Updated docker file
- Improved unit testing
- Updated README
- Use Bottle and Paste for production run

## [0.6.0] - 2020-02-10
### Added
- Loading meta-data from RDF/Turtle file (TTL).
- Automatic Swagger API (OpenAPI) generation (via `flask-restplus`).

### Changed
- Documentation in README
- Use Flask instead of Bottle

### Fixed
- Test syntax
# FAIR Data Point (FDP)

[![PyPI](https://img.shields.io/pypi/v/fairdatapoint)](https://pypi.org/project/fairdatapoint/)
[![Docker Image Version (latest by date)](https://img.shields.io/docker/v/nlesc/fairdatapoint?label=Docker)](https://hub.docker.com/r/nlesc/fairdatapoint)
[![DOI](https://zenodo.org/badge/37470907.svg)](https://zenodo.org/badge/latestdoi/37470907)
[![Research Software Directory](https://img.shields.io/badge/RSD-FAIRDataPoint-red)](https://research-software.nl/software/fairdatapoint)
[![Build_Test](https://github.com/fair-data/fairdatapoint/actions/workflows/build_test.yml/badge.svg)](https://github.com/fair-data/fairdatapoint/actions/workflows/build_test.yml)
[![Coverage Status](https://coveralls.io/repos/github/fair-data/fairdatapoint/badge.svg?branch=master)](https://coveralls.io/github/fair-data/fairdatapoint?branch=master)

## Overview
Python implementation of FAIR Data Point.

FDP is a RESTful web service that enables data owners to describe and to expose their datasets (metadata) as well as data users to discover more information about available datasets according to the [FAIR Data Guiding Principles](http://www.force11.org/group/fairgroup/fairprinciples). In particular, FDP addresses the findability or discoverability of data by providing machine-readable descriptions (metadata) at four hierarchical levels:

*FDP -> catalogs -> datasets -> distributions*

FDP software specification can be found [here](https://github.com/FAIRDataTeam/FAIRDataPoint-Spec/blob/master/spec.md).
Other implementations are also available, e.g. [Java implementation](https://github.com/DTL-FAIRData/FAIRDataPoint)

### Demo server
A demo server of this Python implementation is http://fdp.fairdatapoint.nl/

## Installation

To install FDP, do

From pypi
```bash
pip install fairdatapoint
```

Or from this repo, but note that the in-development version might be unstable,
```bash
git clone https://github.com/fair-data/fairdatapoint.git
cd fairdatapoint
pip install .
```

## Running
```bash
fdp-run localhost 80
```

The [Swagger UI](https://swagger.io/tools/swagger-ui/) is enabled for FDP service, and you can have a try by visiting http://localhost.

## Unit testing
Run tests (including coverage) with:

```bash
pip install .[tests]
pytest
```

## Deploy with Docker

Check [fairdatapoint-service](https://github.com/CunliangGeng/fairdatapoint-service).

## Deploy without Docker

Before deploying FDP, it's necessary to first have a running SPARQL database which can be used to store metadata.

```
pip install fairdatapoint

# fdp-run <host> <port> --db=<sparql-endpoint>
# Let's assume your <host> is 'example.com' and <sparql-endpoint> is 'http://example.com/sparql', then
fdp-run example.com 80 --db='http://example.com/sparql'
```

## Web API documentation

FAIR Data Point (FDP) exposes the following endpoints (URL paths):

| Endpoint |  GET  | POST |  PUT | DELETE     |
|--------------|:--------------:|:-----------------:|:--------------:|:--------------:
| fdp | Output fdp metadata | Create new fdp metadata | Update fdp metadata | Not Allowed |
| catalog     | Output all catalog IDs   | Create new catalog metadata| Not Allowed | Not Allowed |
| dataset     | Output all dataset IDs   | Create new dataset metadata| Not Allowed | Not Allowed |
| distribution  | Output all distribution IDs  | Create new distribution metadata| Not Allowed | Not Allowed |
| catalog/\<catalogID\> | Output \<catalogID\> metadata | Not Allowed | Update \<catalogID\> metadata | Remove \<catalogID\> metadata |
| dataset/\<datasetID\> | Output \<datasetID\> metadata | Not Allowed | Update \<datasetID\> metadata | Remove \<datasetID\> metadata |
| distribution/\<distributionID\> | Output \<distributionID\> metadata | Not Allowed | Update \<distributionID\> metadata | Remove \<distributionID\> metadata |


### Access endpoints to request metadata programmatically

FDP: `curl -iH 'Accept: text/turtle' [BASE URL]/fdp`

Catalog: `curl -iH 'Accept: text/turtle' [BASE URL]/catalog/catalog01`

Dataset: `curl -iH 'Accept: text/turtle' [BASE URL]/dataset/dataset01`

Distribution: `curl -iH 'Accept: text/turtle' [BASE URL]/distribution/dist01`

### FDP supports the following RDF serializations (MIME-types):
* Turtle: `text/turtle`
* N-Triples: `application/n-triples`
* N3: `text/n3`
* RDF/XML: `application/rdf+xml`
* JSON-LD: `application/ld+json`


## Issues and Contributing
If you have questions or find a bug, please report the issue in the
[Github issue channel](https://github.com/fair-data/fairdatapoint/issues).

If you want to contribute to the development of FDP, have a look at the
[contribution guidelines](CONTRIBUTING.rst).############################
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

#. use the search functionality `here <https://github.com/fair-data/fairdatapoint/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/fair-data/fairdatapoint/issues>`__ to see if someone already filed the same issue;
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
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the FAIR Data Point repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
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
reported by contacting the project team at c.martinez@esciencecenter.nl. All
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
.. FAIR Data Point documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FAIR Data Point's documentation!
==========================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
