[![PyPI](https://img.shields.io/pypi/v/fairdatapoint-client)](https://pypi.org/project/fairdatapoint-client/)
[![Documentation Status](https://readthedocs.org/projects/fairdatapoint-client/badge/?version=latest)](https://fairdatapoint-client.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/334084536.svg)](https://zenodo.org/badge/latestdoi/334084536)
[![Build_Test](https://github.com/fair-data/fairdatapoint-client/actions/workflows/build_test.yml/badge.svg)](https://github.com/fair-data/fairdatapoint-client/actions/workflows/build_test.yml)
[![Coverage Status](https://coveralls.io/repos/github/fair-data/fairdatapoint-client/badge.svg?branch=master)](https://coveralls.io/github/fair-data/fairdatapoint-client?branch=master)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=fair-data_fairdatapoint-client&metric=alert_status)](https://sonarcloud.io/dashboard?id=fair-data_fairdatapoint-client)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)


# fairdatapoint-client

### Contents

-   [Overview](#overview)
-   [Installation](#installation)
-   [Quick Tutorial](#Tutorial)
-   [Issues & Contributing](#Issues-and-Contributing)
-   [License](./LICENSE)

## Overview

fairdatapoint-client is a simple and elegant library to interact with
[FAIR Data Point](https://github.com/fair-data/fairdatapoint) resources from
Python, e.g. read and write catalogs, datasets and distributions in an FDP server.

The supported APIs are listed below:

| FDP Layers   | Path Endpoint               | Specific Resource Endpoint              |
|--------------|-----------------------------|-----------------------------------------|
| fdp          | [baseURL] or [baseURL]/fdp  |                                         |
| catalog      | [baseURL]/catalog           | [baseURL]/catalog/[catalogID]           |
| dataset      | [baseURL]/dataset           | [baseURL]/dataset/[datasetID]           |
| distribution | [baseURL]/distribution      | [baseURL]/distribution/[distributionID] |

## Installation

It requires a Python version of 3.7, 3.8 or 3.9.

#### Stable Release

The fairdatapoint-client is available on [PyPI](https://pypi.org/project/fairdatapoint-client/),
you can install it using:

`pip install fairdatapoint-client`

#### Development Version

You can also install from the latest source code, but note that the
in-development version might be unstable:

```{.sourceCode .console}
git clone https://github.com/fair-data/fairdatapoint-client.git
cd fairdatapoint-client
pip install .
```

To run tests (including coverage):

```{.sourceCode .console}
pip install '.[tests]'
pytest
```


## Tutorial

### Using Client
```python
from fdpclient.client import Client

# create a client with base URL
client = Client('http://example.org')

# create metadata
with open('catalog01.ttl') as f:
    data = f.read()
client.create_catalog(data)

# let's assume the catalogID was assigned as 'catalog01'
# read metadata, return a RDF graph
r = client.read_catalog('catalog01')
print(r.serialize(format="turtle").decode("utf-8"))

# update metadata
with open('catalog01_update.ttl') as f:
    data_update = f.read()
client.update_catalog('catalog01', data_update)

# delete metadata
client.delete_catalog('catalog01')
```

### Using operation functions
```python
from fdpclient import operations

# create metadata
with open('catalog01.ttl') as f:
    data = f.read()
operations.create('http://example.org/catalog', data)

# read metadata, return a RDF graph
r = operations.read('http://example.org/catalog/catalog01')
print(r.serialize(format="turtle").decode("utf-8"))

# update metadata
with open('catalog01_update.ttl') as f:
    data_update = f.read()
operations.update('http://example.org/catalog/catalog01', data_update)

# delete metadata
operations.delete('http://example.org/catalog/catalog01')
```

## Issues and Contributing
If you have questions or find a bug, please report the issue in the
[Github issue channel](https://github.com/fair-data/fairdatapoint-client/issues).

If you want to contribute to the development of fairdatapoint-client, have a
look at the [contribution guidelines](CONTRIBUTING.rst).#######################
Contributing guidelines
#######################

We welcome any kind of contribution to our software, from simple comment or
question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_.
Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

#. you have a question;
#. you think you may have found a bug (including unexpected behavior);
#. you want to make some kind of change to the code base (e.g. to fix a bug, to
   add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

#. use the search functionality `here <https://github.com/fair-data/fairdatapoint-client/issues>`__
   to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/fair-data/fairdatapoint-client/issues>`__
   to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue,
   making sure to provide enough information to the rest of the community to
   understand the cause and context of the problem. Depending on the issue, you
   may want to include:

    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_
      of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies
      you're using;
    - information about the operating system;

#. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

#. (**important**) announce your plan to the rest of the community *before you
   start working*. This announcement should be in the form of a (new) issue;
#. (**important**) wait until some kind of consensus is reached about your idea
   being a good idea;
#. if needed, fork the repository to your own Github profile and create your
   own feature branch off of the latest master commit. While working on your
   feature branch, make sure to stay up to date with the master branch by
   pulling in changes, possibly from the 'upstream' repository (follow the
   instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__
   and `here <https://help.github.com/articles/syncing-a-fork/>`__);
#. make sure the existing tests still work by running ``python setup.py test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to
   (your fork of) the fairdatapoint-client repository on GitHub;
#. create the pull request, e.g. following the instructions
   `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know
how to write or run tests for it, or how to generate the documentation: don't
let this discourage you from making the pull request; we can help you! Just go
ahead and submit the pull request, but keep in mind that you might be asked to
append additional commits to your pull request.###############################################################################
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
reported by contacting the project team at c.geng@esciencecenter.nl. All
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
###########
Change Log
###########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

[0.1.0]
*******

* First release==========================
FAIR Data Point API Client
==========================
.. contents::
    :depth: 3

.. meta::
   :description: A collection of tutorials and guides for fairdatapoint-client package.
   :keywords: fairdatapoint, metadata, guide, tutorial, doc

fairdatapoint-client is a simple and elegant library to interact with
`FAIR Data Point <https://github.com/fair-data/fairdatapoint>`_ resources from
Python, e.g. read and write catalogs, datasets and distributions in an FDP server.

The supported APIs are listed below:

.. list-table::
   :widths: 20 30 40
   :header-rows: 1

   * - FDP Layers
     - Path Endpoint
     - Specific Resource Endpoint
   * - fdp
     - [baseURL] or [baseURL]/fdp
     -
   * - catalog
     - [baseURL]/catalog
     - [baseURL]/catalog/[catalogID]
   * - dataset
     - [baseURL]/dataset
     - [baseURL]/dataset/[datasetID]
   * - distribution
     - [baseURL]/distribution
     - [baseURL]/distribution/[distributionID]

Quick Start
-----------

Using Client

.. code-block:: python

    from fdpclient.client import Client

    # create a client with base URL
    client = Client('http://example.org')

    # create metadata
    with open('catalog01.ttl') as f:
        data = f.read()
    client.create_catalog(data)

    # let's assume the catalogID was assigned as 'catalog01'
    # read metadata, return a RDF graph
    r = client.read_catalog('catalog01')
    print(r.serialize(format="turtle").decode("utf-8"))

    # update metadata
    with open('catalog01_update.ttl') as f:
        data_update = f.read()
    client.update_catalog('catalog01', data_update)

    # delete metadata
    client.delete_catalog('catalog01')

Using operation functions

.. code-block:: python


    from fdpclient import operations

    # create metadata
    with open('catalog01.ttl') as f:
    data = f.read()
    operations.create('http://example.org/catalog', data)

    # read metadata, return a RDF graph
    r = operations.read('http://example.org/catalog/catalog01')
    print(r.serialize(format="turtle").decode("utf-8"))

    # update metadata
    with open('catalog01_update.ttl') as f:
        data_update = f.read()
    operations.update('http://example.org/catalog/catalog01', data_update)

    # delete metadata
    operations.delete('http://example.org/catalog/catalog01')


List of Methods/Functions
-------------------------

Client methods
::::::::::::::
.. currentmodule:: fdpclient.client.Client

You can use :class:`fdpclient.client.Client` class and its methods to process FDP metadata.

.. autosummary::

    create_fdp
    create_catalog
    create_dataset
    create_distribution

.. autosummary::

    read_fdp
    read_catalog
    read_dataset
    read_distribution

.. autosummary::

    update_fdp
    update_catalog
    update_dataset
    update_distribution

.. autosummary::

    delete_catalog
    delete_dataset
    delete_distribution

Operation functions
:::::::::::::::::::

You can also use :mod:`fdpclient.operations` functions to process FDP metadata.

.. currentmodule:: fdpclient.operations

.. autosummary::

    create
    read
    update
    delete


Indices and tables
------------------
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`APIs
====

Client
------
.. automodule:: fdpclient.client
    :members:
    :undoc-members:
    :show-inheritance:
    :private-members:

Operations
----------
.. automodule:: fdpclient.operations
    :members:
    :undoc-members:
    :show-inheritance:
    :private-members:

Global Variables
----------------
.. autodata:: fdpclient.config.DATA_FORMATS
