## Badges

| fair-software.eu recommendations | |
| :-- | :--  |
| code repository              | [![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY) |
| license                      | [![github license badge](https://img.shields.io/github/license/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY)](https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY) |
| community registry           | [![RSD](https://img.shields.io/badge/rsd-FHIR_TO_CAPACITY-00a3e3.svg)](https://www.research-software.nl/software/fhir-to-capacity)  |
| citation                     | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4818370.svg)](https://doi.org/10.5281/zenodo.4818370) |
| howfairis                          | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) |
| **Other best practices**           | &nbsp; |
| Coverage                           | [![Coverage Status](https://coveralls.io/repos/github/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY/badge.svg?branch=master)](https://coveralls.io/github/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY?branch=master)|
| **GitHub Actions**                 | &nbsp; |
| Build                              | [![Build](https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY/actions/workflows/build.yml/badge.svg)](https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY/actions/workflows/build.yml)|

## Installation

To install FHIR-to-CAPACITY from GitHub repository, do:

```console
git clone https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY.git
cd FHIR-to-CAPACITY
python3 -m pip install .
```



## Usage

To upload records from a FHIR server to a CAPACITY registry, execute the following command:
```console
fhir-to-capacity FHIR_BASE_URL CAPACITY_API_URL CAPACITY_TOKEN
```

Using the `--help` option will give you more information:
```console
fhir-to-capacity --help

Usage: scripts/fhir_to_capacity fhir-base-url capacity-url capacity-token

Queries a FHIR server for patients, converts this to records matching the CAPACITY codebook, and uploads the result to a REDCAP server.

Arguments:
  fhir-base-url
  capacity-url
  capacity-token

Other actions:
  -h, --help       Show the help

```

## Development
There is an additional command available for filling a FHIR server with synthetic data. It can 
be run as follows:
```console
fill-server FHIR_BASE_URL
```

The `--help` option will show you more information:
```console
Usage: scripts/fill_server fhir-base [n]

Fills the server with n ramdom patients.

Arguments:
  fhir-base    the url to the FHIR api
  n            number of patients to create (type: INT, default: 10)

Other actions:
  -h, --help   Show the help

```
## Contributing

If you want to contribute to the development of FHIR-to-CAPACITY,
have a look at the [contribution guidelines](CONTRIBUTING.rst).

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [NLeSC/python-template](https://github.com/NLeSC/python-template).
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

#. use the search functionality `here <https://github.com/FAIR-data-for-CAPACITY /FHIR-to-CAPACITY/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/FAIR-data-for-CAPACITY /FHIR-to-CAPACITY/issues>`__ to see if someone already filed the same issue;
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
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the FHIR-to-CAPACITY repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###########
Change Log
###########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

[Unreleased]
************

Added
-----

* Empty Python project directory structure
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
reported by contacting the project team at d.smits@esciencecenter.nl. All
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
