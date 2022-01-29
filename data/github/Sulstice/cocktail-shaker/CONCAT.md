---
title: 'Cocktail Shaker: An open source drug expansion and enumeration library for peptides'
tags:
  - Python
  - Cheminformatics
  - RDKit
  - Peptides
authors:
  - name: Suliman Sharif
    orcid: 0000-0002-1342-9258
    affiliation: 1
affiliations:
 - name: None
   index: 1
date: 23 November 2019
bibliography: paper.bib
---

# Introduction

Without expensive software, the rapid creation and design of peptide ligand libraries has been a
challenge for many drug discovery scientists. Currently, protein- and peptide-based therapeutics constitute 10% of the 
pharmaceutical market and will make up a larger proportion of the market in the future [@Craig:2013-1; @Bruno:2013-2]. With the demand for designing new peptide therapeutics on the rise, new high throughput peptide-specific informatic tools are needed. Currently, two platforms for this purpose exist: Molecular Operating 
Environment (MOE) [@Reynolds:2010-3] and Rapid Peptides Generator (RPG) [@Maillet:2020-4] but both have drawbacks. MOE works efficiently for creating peptide molecule 3D chemical 
files in one particular format (mol2), but at the high cost for a licence. RPG, although free of charge, does not account
for non-natural amino acids and production of multiple chemical files. In this study, I present the first open source python package,
 ```Cocktail Shaker```, developed for exploring, expanding, and synthesizing chemical peptide data.

# Methodology and Implementation

```Cocktail Shaker``` operates within the RDKit platform [@Landrum:2019-5] and is designed for the chemically-oriented
computational research community. RDKit utilizes C++-based functions for speed and rapid creation of molecule objects.
The toolkit offers a variety of utilities that includes: parsing and producing ready-to-use scientific files
designed for any chemical software, employing click chemistry methods for ease of exchange compounds, chemical data writing, 
and chemical representation enumeration employed for machine learning.

```Cocktail Shaker``` consists of three major class objects available to the user: PeptideMolecule, CocktailShaker, and FileWriter.

Using string manipulation PeptideMolecule can build SMILES strings with allocated slots defined by the user. The user can
then enter the produced SMILES into the CocktailShaker object with a library of ligands represented by smiles and optional arguments
of whether to include generation of stereoisomers and/or natural amino acids. ```Cocktail Shaker``` will generate all combinations of the library and allocate them to a slot within the peptide. This process of 
string manipulation is presented in ```Figure 1```.

![Full string manipulation diagram of how ```Cocktail Shaker``` works with a ligand library of just bromine and iodine. 1D representations are labeled above with their 2D depictions displayed below.](https://raw.githubusercontent.com/Sulstice/cocktail-shaker/master/images/figure_1.png)
  

```Cocktail Shaker``` also allows for File Writing of the molecules into a wide array of chemical formats (found in the documentation).
```Cocktail Shaker ``` uses RDKit to convert from 1D to 2D and the CIR Resolver built from webchem to convert 1D SMILES to 3D. At the request 
of the user, the data is saved into one large data file or separated using the keyword fragmentation. This additional API
allows the user the flexibility to write a variety of files to implement in their respective chemical software.  

# Conclusion

Using  ```Cocktail Shaker```, individual research groups and companies can quickly construct private compound collections and progressively improve public
libraries with increased variations of chemical compound data.

```Cocktail Shaker``` with its first version release provides a basis for drug expansion and enumeration. For future
releases ```Cocktail Shaker``` will be expanding into specifying shapes of compounds and, recently partnered with MolPort,
vendor information on any compound generated. It was presented at the RDKit UGM conference at the University of Hamburg
to the cheminformatics community with positive feedback with its second version 1.0.1. With incorporated feedback it will now be
released with version 1.1.8. With more contributions ```Cocktail Shaker``` will be an exciting tool for drug library creation 
and drug discovery for scientists and engineers alike. 

# Acknowledgements

We acknowledge contributions from Ryland Forsythe as an academic consultation, Marvin Corro for quality assurance, Rose Gierth
for technical documentation, and Elena Chow for her work on the graphics. 

# References
Cocktail Shaker: Drug Expansion and Enumeration for Peptides!
=============================================================

[![Build](https://travis-ci.org/Sulstice/Cocktail-Shaker.svg?branch=master)](https://travis-ci.org/Sulstice/Cocktail-Shaker)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Coverage](https://coveralls.io/repos/github/Sulstice/Cocktail-Shaker/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/Cocktail-Shaker?branch=master)
![Python](https://img.shields.io/badge/python-3.6-blue.svg)
[![Gitter](https://badges.gitter.im/Cocktail-Shaker/community.svg)](https://gitter.im/Cocktail-Shaker/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Zenodo](https://zenodo.org/badge/170644606.svg)](https://zenodo.org/badge/latestdoi/170644606)
[![Documentation Status](https://readthedocs.org/projects/cocktail-shaker/badge/?version=latest)](https://cocktail-shaker.readthedocs.io/en/latest/?badge=latest)
[![status](https://joss.theoj.org/papers/c2e1d3c408a5729d832b34ac680d6305/status.svg)](https://joss.theoj.org/papers/c2e1d3c408a5729d832b34ac680d6305)
[![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)


<p align="center">
  <img width="200" height="300" src="images/logoshaker.png">
</p>

cocktail-shaker is a **high-performance drug enumeration and expansion
library**. cocktail-shaker leverages the computational power of
**RDKit** to create and enumerate large volumes of drug compounds. 

Announcements
=============

-   **Release!** Version 1.0.0-beta, August 26, 2019
-   **RDkit UGM 2019**: talk at September 25th at the University of
    Hamburg, Germany.

Using Cocktail Shaker
=====================

[![asciicast](https://asciinema.org/a/7h2s2BQKBWeVuASnutnNjagDW.svg)](https://asciinema.org/a/7h2s2BQKBWeVuASnutnNjagDW)


cocktail-shaker is a young library under heavy development at this time.
It targets two categories of users:

1.  **Users familiar with RDKit**, or those willing to learn RDKit, who want to
    create fast sets of data for high throughput screening or machine
    learning.
2.  **Open-Science Scientists without any knowledge of RDKit** who are
    seeking a a high-level wrapper to create chemical files for their
    software.

If you're in the first category, then you can already start using RDKit.
cocktail-shaker offers a Pythonic, easy-to-use library and you can start
channeling molecules in the expansion library. Instead of validating the
sanity of the data, cocktail-shaker takes care of that for you. With
each molecule being generated it will head into a 1D and/or 2D
validation check (3D not supported yet).

If you're in the second category, we're starting to build experimental
high-level python code to take care a lot of the underpinnings of RDKit.

Installation 
==================

Cocktail-shaker runs on Python 3.5+ and within a conda environment due to the RDKit dependency. A conda installation is coming soon (as soon as I figure it out). 

For the time being, here is the following:

- [Anaconda](https://docs.anaconda.com/anaconda/install/)
- [RDKit](https://www.rdkit.org/docs/Install.html) Version: 2019.09.1+

Cocktail shaker is distributed through PyPi and can be installed via:

`your_env/bin/python -m pip install cocktail-shaker`



Development Installation
========================

As cocktail-shaker is under heavy development at this time, we highly
recommend developers to use the development version on Github (master
branch). You need to clone the repository and install cocktail-shaker with

`python setup.py install`.

As a one-liner, assuming git is installed:

    git clone https://github.com/Sulstice/Cocktail-Shaker.git

This will automatically install the latest version of cocktail-shaker.

Structure of cocktail-shaker
============================

Currently, the main subpackages are:

-   **cocktail_shaker**: Contains a lot of the high level functionality; request
    handling, file parsing/writing, enumeration, and expansion.
-   **docs**: An access point for the readthedocs implementation.
-   **cocktail_shaker/datasources**: This is where the system stores its data on
    predfined functional groups and/or shapes (coming soon).
-   **tests**: Tests that are for the file handling, requests, and
    testing molecule pattern recognition.

The API of all public interfaces are subject to change in the future,
although **datasources** are *relatively* stable at this point.

Genesis
=======

cocktail-shaker began when one developer/scientist wanted an open source
drug library.

- Lead Developer [Suliman sharif](http://sulstice.github.io/)
- Artwork [Elena Chow](http://www.chowelena.com/)
- Technical Documentation [Rose Gierth](https://www.linkedin.com/in/rose-gierth-69a4a083/)
- QA Tester [Marvin Corro](https://www.linkedin.com/in/marvincorro/)

Now cocktail-shaker looks to build on the expertise of these
developers/scientists and the broader open-science community to build an
effective drug library.

* * * * *

External links
==============

-   [Documentation](http://cocktail-shaker.readthedocs.org)

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
reported by contacting the project team at shairsuliman1@gmail.com. All
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
