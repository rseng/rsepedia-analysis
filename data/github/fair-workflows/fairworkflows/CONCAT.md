# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2021-06-25

### Added
* Added CHANGELOG.md, CODE_OF_CONDUCT.md and CONTRIBUTING.md

### Changed
* Upgrade to nanopub v1.2.7 (among other things, to fix click bug)
![Build Status](https://github.com/fair-workflows/fairworkflows/workflows/Python%20application/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/fairworkflows/badge/?version=latest)](https://fairworkflows.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/fair-workflows/fairworkflows/badge.svg?branch=main)](https://coveralls.io/github/fair-workflows/fairworkflows?branch=main)
[![PyPI version](https://badge.fury.io/py/fairworkflows.svg)](https://badge.fury.io/py/fairworkflows)
[![DOI](https://zenodo.org/badge/244369407.svg)](https://zenodo.org/badge/latestdoi/244369407)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fair-workflows/fairworkflows/main?filepath=examples%2Ffairworkflows-quick-start.ipynb)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/4804/badge)](https://bestpractices.coreinfrastructure.org/projects/4804)

# ```fairworkflows``` python library
`fairworkflows` is a high-level, user-friendly python library that supports the construction,
manipulation and publishing of FAIR scientific workflows using semantic technologies. 

## Background
`fairworkflows` is developed as a component of the FAIR Workbench, as part of the FAIR is FAIR project. 

The focus is on description of workflows consisting of manual and computational steps using semantic technology, 
such as the ontology described in the publication:

_Celebi, R., Moreira, J. R., Hassan, A. A., Ayyar, S., Ridder, L., Kuhn, T., & Dumontier, M. (2019). Towards FAIR protocols and workflows: The OpenPREDICT case study._ [_arXiv:1911.09531._](https://arxiv.org/abs/1911.09531)

The goals of the project are:
1. To facilitate the construction of RDF descriptions of a variety of scientific 'workflows', in the most general sense. This includes experimental procedures, ipython notebooks, computational analysis of results, etc.
2. To allow validation and publication of the resultant RDF (for example, by means of nanopublications).
3. Re-use of previously published steps, in new workflows.
4. FAIR data flow from end-to-end.

We seek to provide an easy-to-use python interface for achieving the above.

## Documentation
Checkout the [user documentation](https://fairworkflows.readthedocs.io/).

## Installation

The most recent release can be installed from the python package index using ```pip```:

```
pip install fairworkflows
```

To publish workflows to the nanopub server you need to setup your nanopub profile. This
allows the nanopub server to identify you. Run the following in the terminal after installation:
```
setup_nanopub_profile
```
This will add and store RSA keys to sign your nanopublications, publish a
nanopublication with your name and ORCID iD to declare that you are
using using these RSA keys, and store your ORCID iD to automatically add
as author to the provenance of any nanopublication you will publish
using this library.

## Quick demo
Try out the library in this online executable notebook: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fair-workflows/fairworkflows/main?filepath=examples%2Ffairworkflows-quick-start.ipynb)

## Quick Start
### Import from `fairworkflows` library
```python
from fairworkflows import is_fairworkflow, is_fairstep, FairWorkflow
```

### Define a step for your workflow
Mark a function as a FAIR step using the `is_fairstep` decorator.
Use keyword arguments to semantically annotate the step. 
In this example to provide a label and describe that this is a script task.
```python
@is_fairstep(label='Addition', is_script_task=True)
def add(x: float, y: float) -> float:
    """Adding up numbers."""
    return x + y
```
### Define your workflow
Define your workflow by calling previously defined step functions. 
Mark the function as a workflow using the `is_fairworkflow` decorator.
```python
@is_fairworkflow(label='My Workflow')
def my_workflow(in1, in2):
    """
    A simple workflow
    """
    t1 = add(in1, in2)
    return t1
```
### Construct and publish a workflow
Construct a FairWorkflow object from the function defining the workflow and publish as nanopublication.
```python
workflow = FairWorkflow.from_function(my_workflow)
workflow.publish_as_nanopub(use_test_server=True, publish_steps=True)
```

### Execute the workflow
Execute the workflow and inspect the prospective provenance
```python
result, prov = workflow.execute(1, 4)
print(prov)
```

### Example notebook
* See [examples/fairworkflows-quick-start.ipynb](examples/fairworkflows-quick-start.ipynb) for a current example of using the fairworkflows library to build a workflow using plex rdf

## How is the ```fairworkflows``` library expected to be used?
While this library could be used as a standalone tool to build/publish RDF workflows,
it is intended more as a component to be used in a variety of other tools that seek to add FAIR elements to workflows. At present the library is used in the following tools:

* [FAIRWorkflowsExtension](https://github.com/fair-workflows/FAIRWorkflowsExtension): A Jupyter Lab extension that adds a widget for searching for previously published FairSteps or FairWorkflows. These can then be loaded into the notebook for modification or combination into new workflows.

It is expected that the library will soon interact with FAIR Data Points as well e.g. [fairdatapoint](https://github.com/NLeSC/fairdatapoint).

## Relation to existing workflow formats/engines (e.g. CWL, WDL, Snakemake etc)
This library is not intended to replace or compete with the hundreds of existing computational workflow formats, but rather to aid in RDF description and comparison of workflows in the most general sense of the term (including manual experiemental steps, notebooks, and so on). Steps in a FAIRWorkflow may very well be 'run this CWL workflow' or 'run this script', so such workflows are expected to sit more on a meta-level, describing the before-and-after of running one of these fully automated computational workflows as well.
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
reported by contacting the project team at r.richardson@esciencecenter.nl. All
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
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/fair-workflows/fairworkflows/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/fair-workflows/fairworkflows/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. install the dependencies required to run `fairworkflows` in development:

    ```bash
    pip install -r requirements.txt
    ```

1. make sure the existing tests still work by running ``pytest``. Note that any pull requests to the fairworkflows repository on github will automatically trigger running of the test suite;
1. check that the code is in accordance with the PEP8 style guide, by running ``flake8 . --count --show-source --statistics``,
configuration is in `tox.ini`.
1. add your own tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the fairworkflows repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
