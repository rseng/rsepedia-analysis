

## Removed
 
 * `noodles.file` unused module
 
 
## Changed

 * Added `SerPath` class to serial namespace


## Fixed
 * `Fail` class serialization
---
title: Noodles - parallel programming in Python
---

[![Travis](https://travis-ci.org/NLeSC/noodles.svg?branch=master)](https://travis-ci.org/NLeSC/noodles)
[![Zenodo DOI](https://zenodo.org/badge/45391130.svg)](https://zenodo.org/badge/latestdoi/45391130)
[![Code coverage](https://codecov.io/gh/NLeSC/noodles/branch/master/graph/badge.svg)](https://codecov.io/gh/NLeSC/noodles)
[![Documentation](https://readthedocs.org/projects/noodles/badge/?version=latest)](https://noodles.readthedocs.io/en/latest/?badge=latest)

::: {.splash}
* Write readable code
* Parallelise with a dash of Noodle sauce!
* Scale your applications from laptop to HPC using Xenon
    + [Learn more about Xenon](https://xenon-middleware.github.io/xenon)
* Read our [documentation](https://noodles.rtfd.io/), including tutorials on:
    + [Creating parallel programs](https://noodles.readthedocs.io/en/latest/poetry_tutorial.html)
    + [Circumventing the global interpreter lock](https://noodles.readthedocs.io/en/latest/prime_numbers.html)
    + [Handling errors in a meaningful way](https://noodles.readthedocs.io/en/latest/errors.html)
    + [Serialising your data](https://noodles.readthedocs.io/en/latest/serialisation.html)
    + [Functional programming and flow control](https://noodles.readthedocs.io/en/latest/control_your_flow.html)
:::

# What is Noodles?

Noodles is a task-based parallel programming model in Python that offers the same intuitive interface when running complex workflows on your laptop or on large computer clusters.

# Installation
To install the latest version from PyPI:

```
pip install noodles
```

To enable the Xenon backend for remote job execution,

```
pip install noodles[xenon]
```

This requires a Java Runtime to be installed, you may check this by running

```
java --version
```

which should print the version of the currently installed JRE.


# Documentation
All the latest documentation is available on [Read the Docs](https://noodles.rtfd.io/).

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
reported by contacting the project team at j.hidding@esciencecenter.nl. All
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
Thank you for being interested to contribute to Noodles!

## Issues / Pull requests
If you find a bug or unexpected behaviour in Noodles you are welcome to open an issue.
We'll do our best to address the issue if it is within our capacity to do so.

Pull requests are certainly welcome. Please first open an issue outlining the bug or feature request that is being addressed.

All contributions will be integrated per the Apache License agreement (see LICENSE); contributing authors will so be attributed.
# Howto run

From this folder say:
    > PYTHONPATH=".." python3 example<n>.py

