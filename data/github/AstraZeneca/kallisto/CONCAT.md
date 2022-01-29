---
title: 'kallisto: A command-line interface to simplify computational modelling and the generation of atomic features'
tags:
  - Python
  - Poetry
  - Atomic Features
  - Machine Learning
  - Computational Chemistry
  - Pharmaceutical Science

authors:
  - name: Eike Caldeweyher^[Corresponding author]
    orcid: 0000-0002-3985-595X
    affiliation: "1" 

affiliations:
 - name: Data Science and Modelling, Pharmaceutical Science, R&D, AstraZeneca, Gothenburg, Sweden
   index: 1

date: 8 March 2021

bibliography: paper.bib

---

# Statement of Need

Machine learning (ML) has recently become very popular within pharmaceutical industry [@roy2015; @sprous2010].
Tasks as, e.g., building predictive models, performing virtual screening, or predicting compound activities are potential use cases for such ML applications [@li2021; @simm2018].
Traditionally, ML models often rely on the quantitative structure-activity relationship (QSAR) that has been popularized by medicinal chemists and statisticians to relate bioactivities to specific functional group manipulations [@dudek2006; @verma2010].
This QSAR approach decreases the dimensionality of the underlying problem and projects the molecular structure into a space spanned by the physicochemical features.
While early approaches relied more on linear regression, modern approaches combine such features with non-linear ML algorithms.

Chemoinformatic packages like RDKit [@landrum2006] enable the fast calculation of atomic/molecular features based on structural information like the molecular graph, while recently an extended Hueckel package has been added as well [@landrum2019].
However, frequently we want to go beyond a structure-only approach thus incorporating electronic structure effects as obtained, e.g., by a (higher-level) quantum mechanical (QM) treatment.
The calculation of QM-based features relies often on well-established quantum chemistry methods like Kohn-Sham density functional theory (DFT) that is currently the workhorse of computational chemistry [@parr1980; @kohn1999].
However, generating the feature space by DFT is computationally demanding and can become the computational bottleneck especially when aiming for high-throughput experiments with several hundred to thousands of molecules.

Since there exists a critical need for an efficient yet accurate featurizer, we developed the ``kallisto`` command-line interface that is able to calculate QM-based atomic features for atoms and molecules efficiently (whole periodic table up to Radon).
The features are either interpolating high-level references (e.g., static/dynamic polarizabilities with time-dependent DFT data) or are parametrized [@caldeweyher2019] to reproduce QM references (e.g., DFT Hirshfeld [@hirshfeld1977] atomic partial charges).
Molecular geometries need to have an [``xmol``](https://en.wikipedia.org/wiki/XYZ_file_format) or a [``Turbomole``](https://www.turbomole.org/wp-content/uploads/2019/11/Turbomole_Manual_7-4-1.pdf) like format to be processed by ``kallisto``.
Besides, we implemented several computational modelling helpers to simplify the development of high-throughput procedures.
Some of those modelling helpers depend on the open-source [xtb](https://github.com/grimme-lab/xtb) tight-binding scheme that has been developed by Stefan Grimme and co-worker [@bannwarth2020].
The ``kallisto`` software depends on the scientific libraries Numpy [@harris2020] and SciPy [@virtanen2020].
The [online documentation](https://ehjc.gitbook.io/kallisto/) covers all high-level functionalizations of this software mostly in terms of copy-paste recipes.
Furthermore, we cover bits of the underlying theory and compare to experimental data as well as to other modern deep learning models.

# Atomic and Molecular Features

The following atomic and molecular features are available for all atoms up to Radon

- Coordination numbers [@grimme2010; @caldeweyher2019]
- Proximity shells 
- Environment-dependent electronegativity equilibration partial charges [@caldeweyher2019]
- Environment- and charge-dependent dynamic polarizabilities [@grimme2010; @caldeweyher2019]
- Environment- and charge-dependent van-der-Waals radii [@fedorov2018; @rahm2017; @mantina2009]
- Sterimol descriptors (L, Bmin, Bmax) [@brethome2019]

# Modelling Helpers

The following modelling helper are implemented

- Breadth first sorting
- Root mean squared deviation (quaternions) [@coutsias2004]
- Substructure identifier
- Substructure exchanger 

# Acknowledgements

EC acknowledges contributions from Philipp Pracht (`@pprcht`) and thanks Kjell Jorner (`@kjelljorner`) for sharing his Sterimol algorithm.

# References
<div align="center">
<img src="./assets/logo.svg" alt="Kallisto" width="300">
</div>

##

![PyPI - Python Version](https://img.shields.io/pypi/pyversions/kallisto)
[![Documentation](https://img.shields.io/badge/GitBook-Docu-lightgrey)](https://app.gitbook.com/@ehjc/s/kallisto/)
[![Maturity Level](https://img.shields.io/badge/Maturity%20Level-Under%20Development-orange)](https://img.shields.io/badge/Maturity%20Level-Under%20Development-orange)
[![Tests](https://github.com/AstraZeneca/kallisto/workflows/Tests/badge.svg)](https://github.com/AstraZeneca/kallisto/actions?workflow=Tests)
[![codecov](https://codecov.io/gh/AstraZeneca/kallisto/branch/master/graph/badge.svg?token=HI0U0R96X8)](https://codecov.io/gh/AstraZeneca/kallisto)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/AstraZeneca/kallisto.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/AstraZeneca/kallisto/context:python)
[![status](https://joss.theoj.org/papers/16126cbcfb826bf4810d243a009a6b02/status.svg)](https://joss.theoj.org/papers/16126cbcfb826bf4810d243a009a6b02)

# Table of Contents

- Full Author List
- Introduction
- Installation
- Testing suite
- Reference

# Full Author List

- Developer [Eike Caldeweyher](https://scholar.google.com/citations?user=25n8C3wAAAAJ&hl)
- Developer [Rocco Meli](https://scholar.google.com/citations?hl=de&user=s8cVcvYAAAAJ)
- Developer [Philipp Pracht](https://scholar.google.com/citations?user=PJiGPk0AAAAJ&hl)

# Introduction

We developed the `kallisto` program for the efficient and robust calculation of atomic features using molecular geometries either in a ``xmol`` or a ``Turbomole`` format.
Furthermore, several modelling tools are implemented, e.g., to calculate root-mean squared deviations via quaternions (including rotation matrices), sorting of molecular geometries and many more. All features of ``kallisto`` are described in detail within our [documentation](https://app.gitbook.com/@ehjc/s/kallisto/) ([GitBook repository](https://github.com/f3rmion/gitbook-kallisto)).

Main dependencies
-----------------

```bash
click 7.1.2 Composable command line interface toolkit
numpy 1.20.1 NumPy is the fundamental package for array computing with Python.
scipy 1.6.0 SciPy: Scientific Library for Python
└── numpy >=1.16.5
```

For a list of all dependencies have a look at the pyproject.toml file.

Installation from PyPI
----------------------

To install ``kallisto`` via `pip` use our published PyPI package
```bash
pip install kallisto
```

Installation from Source
------------------------

Requirements to install ``kallisto``from sources:
- [poetry](https://python-poetry.org/docs/#installation)
- [pyenv](https://github.com/pyenv/pyenv#installation) or [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- python >=3.7

First check that ``poetry`` is running correctly (v1.0.10 at the time of writing)

```bash
> poetry --version
Poetry version 1.0.10
```

Create a virtual environment (via ``pyenv`` or ``conda``) and activate it. Afterwards, clone the ``kallisto`` project from GitHub and install it using ``poetry``

```bash
> git clone git@github.com:AstraZeneca/kallisto.git
> cd kallisto
> poetry install
```

Testing suite
-------------

The ``kallisto`` project uses [nox](https://nox.thea.codes/en/stable/tutorial.html#installation) as an automated unit test suite, which is therefore an additional dependency.

### Default nox session

The default session includes: linting (lint), type checks (mypy, pytype), and unit tests (tests). 

```bash
> nox
```

When everything runs smoothly through, you are ready to go! After one successful nox run, we can reuse the created virtual environment via the ``-r`` flag.

```bash
> nox -r
```

Different unit test sessions are implemented (check the noxfile.py). They can be called separately via the run session ``-rs`` flag.

### Tests

Run all unit tests that are defined in the /tests directory.

```bash 
> nox -rs tests
```

### Lint

``kallisto`` uses the [flake8](https://flake8.pycqa.org/en/latest/) linter (check the .flake8 config file).

```bash
> nox -rs lint
```

### Black

``kallisto`` uses the [black](https://github.com/psf/black) code formatter.

```bash 
> nox -rs black
```

### Safety

``kallisto`` checks the security of dependencies via [safety](https://pyup.io/safety/).

```bash
> nox -rs safety
```

### Mypy

``kallisto`` checks for static types via [mypy](https://github.com/python/mypy) (check the mypy.ini config file).

```bash
> nox -rs mypy
```

### Pytype

``kallisto`` furthermore uses [pytype](https://github.com/google/pytype) for type checks.

```bash
> nox -rs pytype
```

### Coverage

Unit test [coverage](https://coverage.readthedocs.io/en/coverage-5.4/) can be checked as well.


```bash
> nox -rs coverage
```

Reference
---------

Always cite:

Eike Caldeweyher, J. Open Source Softw., *2021*, 6, 3050. DOI: [10.21105/joss.03050](https://doi.org/10.21105/joss.03050)

```
@article{Caldeweyher2021,
  doi = {10.21105/joss.03050},
  url = {https://doi.org/10.21105/joss.03050},
  year = {2021},
  volume = {6},
  number = {60},
  pages = {3050},
  author = {Eike Caldeweyher},
  title = {kallisto: A command-line interface to simplify computational modelling and the generation of atomic features},
  journal = {J. Open Source Softw.}
}
```
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
reported by contacting the project team at hello@eikecaldeweyher.de. All
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
# Contributing to kallisto

Please take a moment to review this guidelines to make the contribution process
simple and effective for all involved.

Respecting these guidelines helps communicate that you respect the time of
the developers who manage and develop this open source project.
In return, they should return this respect by addressing your problem,
evaluating changes, and helping you handle your pull requests.

## Reporting a Bug

A bug is a *demonstratable problem* caused by the code in this repository.

Before opening a bug report:

1. Check if the issue has already been reported.
2. Check if it still is an issue or has already been fixed?
   Try to reproduce it with the latest version from the `master` branch.
3. Isolate the problem and create a reduced test case.

A good bug report should not leave others needing to chase you up for more
information. So please try to be as detailed as possible in your report,
answer at least these questions:


1. Which version of `kallisto` are you using? The current version is always
   a subject to change, so be more specific.
2. What is your environment (your laptop, the cluster of the university)?
3. What steps will reproduce the issue?
   We have to reproduce the issue, so we need all the input files.
4. What would be the expected outcome?
5. What did you see instead?

All these details will help people to fix any potential bugs.


## Suggesting a New Feature

Feature requests are welcome. But take a moment to find out if your idea fits
the scope and goals of the project. It is up to you to provide a strong
argument to convince the project's developers of the benefits of this feature.
Please provide as much detail and context as possible.


## Implementing a New Feature


Contributions are welcome via Github pull requests.


- Each pull request should implement *one* feature or fix *one* bug.
  If you want to add or fix more than one thing, submit more than one
  pull request.
- Do not commit changes to files that are irrelevant to your feature or
  bugfix (*e.g.* `.gitignore`).
- Be willing to accept criticism and work on improving your code.
- Do not add third-party dependencies, the `kallisto` program should be usable as
  standalone.
- Make sure the tests run successful on more than
  your local machine (*e.g.* on cluster of your university).
  
  Please sign-off your commits
  
  
### For New Contributors

If you never created a pull request before you can learn how to do this 
from [this great tutorial](https://app.egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github)

Don't know where to start?
You can start by looking through these [help-wanted issues](https://github.com/f3rmion/kallisto/labels/help%20wanted).

## Sign Your Work

The sign-off is a simple line at the end of the explanation for a commit. All 
commits needs to be signed. Your signature certifies that you wrote the patch or
otherwise have the right to contribute the material. The rules are pretty simple,
if you can certify the below (from [developercertificate.org](https://developercertificate.org/)):

```
Developer Certificate of Origin
Version 1.1

Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
1 Letterman Drive
Suite D4700
San Francisco, CA, 94129

Everyone is permitted to copy and distribute verbatim copies of this
license document, but changing it is not allowed.

Developer's Certificate of Origin 1.1

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I
    have the right to submit it under the open source license
    indicated in the file; or

(b) The contribution is based upon previous work that, to the best
    of my knowledge, is covered under an appropriate open source
    license and I have the right under that license to submit that
    work with modifications, whether created in whole or in part
    by me, under the same open source license (unless I am
    permitted to submit under a different license), as indicated
    in the file; or

(c) The contribution was provided directly to me by some other
    person who certified (a), (b) or (c) and I have not modified
    it.

(d) I understand and agree that this project and the contribution
    are public and that a record of the contribution (including all
    personal information I submit with it, including my sign-off) is
    maintained indefinitely and may be redistributed consistent with
    this project or the open source license(s) involved.
```

Then you just add a line to every git commit message:

    Signed-off-by: Joe Smith <joe.smith@example.com>

Use your real name (no pseudonyms or anonymous contributions possible.)

If you set your `user.name` and `user.email` git configs, you can sign your
commit automatically with `git commit -s`.

## Contributors

The `kallisto` program has been developed by:

- E. Caldeweyher (@f3rmion)
- R. Meli (@RMeli)
- P. Pracht (@pprcht)
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
