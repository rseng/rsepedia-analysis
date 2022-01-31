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
.. Cocktail Shaker documentation master file, created by
   sphinx-quickstart on Wed Aug 28 23:41:51 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Cocktail Shaker's documentation!
===========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. Cocktail Shaker documentation master file, created by sphinx-quickstart on Tue Mar 24 16:12:38 2015.

Cocktail Shaker
===============

.. sectionauthor:: Suliman Sharif <sharifsuliman1@gmail.com>

Cocktail Shaker is a **high-performance drug enumeration and expansion library**.
Cocktail Shaker leverages the computational power of  **RDKit** to create and enumerate large volumes of drug compounds.

>>> from cocktail_shaker import Cocktail, FileWriter
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> FileWriter('example', combinations, 'mol2')

Cocktail Shaker makes your drug enumeration and expansion life easy. It also generates your files for you in as many
formats needed for any cheminformatics software.

Features
--------


- File parsing of TXT, SDF, and Chemical Smiles.
- File writing in a variety of formats some of which include: cif, sdf, pdb, mol, mol2 and many others
- Ability to recognize and expand libraries of compounds some of which include halogens, acyl halides, aldehydes.
- Ability to enumerate in 1D, and 2D structures and produce those compounds.
- Supports Python versions 3.3+.
- Released under the `MPL license`_.

User guide
----------

A step-by-step guide to getting started with Cocktail Shaker.

.. toctree::
    :maxdepth: 2

    guide/installation
    guide/quickstart

API documentation
-----------------

Comprehensive API documentation with information on every function, class and method.

.. toctree::
    :maxdepth: 2

    guide/file_formats
    guide/file_handling
    guide/functional_groups
    guide/amino_acids
    guide/cocktail
    guide/peptide_builder
    guide/contributing

Cocktail Shaker's license
-------------------------

Cocktail Shaker is released under the Mozilla Public License 2.0. This is a short, permissive software license that allows commercial use,
modifications, distribution, sublicensing and private use. Basically, you can do whatever you want with Cocktail Shaker as long as
you include the original copyright and license in any copies or derivative projects.

.. _`MPL license`: https://github.com/Sulstice/Cocktail-Shaker/blob/master/LICENSE
.. _contributing:

Contributing
============

.. sectionauthor:: Suliman Sharif <sharifsuliman1@gmail.com>

Contributing
------------

Contributions of any kind to cocktail-shaker are greatly appreciated! All contributions are welcome, no matter how big or small.

If you are able to contribute changes yourself, just fork the `source code`_ on GitHub, make changes, and file a pull
request.

Quick Guide to Contributing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. `Fork the Cocktail-Shaker repository on GitHub`_, then clone the fork to your local machine::

    git clone https://github.com/<username>/Cocktail-Shaker.git

2. Install the development requirements::

    cd Cocktail-Shaker
    pip install -r requirements/development.txt

3. Create a new branch for your changes::

    git checkout -b <name-for-changes>

4. Make your changes or additions. Ideally create some tests and ensure they pass.

5. Commit your changes and push to your fork on GitHub::

    git add .
    git commit -m "<description-of-changes>"
    git push origin <name-for-changes>

4. `Submit a pull request`_.

Tips
~~~~

- Follow the `PEP8`_ style guide.
- Include docstrings as described in `PEP257`_.
- Try and include tests that cover your changes.
- Try to write `good commit messages`_.
- Consider `squashing your commits`_ with rebase.
- Read the GitHub help page on `Using pull requests`_.

.. _`source code`: https://github.com/Sulstice/Cocktail-Shaker.git
.. _`Fork the Cocktail-Shaker repository on GitHub`: https://github.com/Sulstice/Cocktail-Shaker/fork
.. _`Submit a pull request`: https://github.com/Sulstice/Cocktail-Shaker/compare/
.. _`squashing your commits`: http://gitready.com/advanced/2009/02/10/squashing-commits-with-rebase.html
.. _`PEP8`: https://www.python.org/dev/peps/pep-0008
.. _`PEP257`: https://www.python.org/dev/peps/pep-0257
.. _`good commit messages`: http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html
.. _`Using pull requests`: https://help.github.com/articles/using-pull-requests
.. _fileformats:

File Formats
============

Input can be parsed in a variety of file formats that are autodetected::

    >>> FileParser('compound.sdf')

    ['RDKit MolObject']

The full list of file formats that cocktail-shaker supports parsing::

    sdf       # 2D Structure file format

Output can be returned in a variety of file formats that are specified your enumerated compounds and your
desired file extension::

    >>> FileWriter('test'', compounds, '.mol2')

    Writes all the compounds to 1 mol2 file.


The full list of file formats::

    alc         # Alchemy format
    cdxml       # CambridgeSoft ChemDraw XML format
    cerius      # MSI Cerius II format
    charmm      # Chemistry at Harvard Macromolecular Mechanics file format
    cif         # Crystallographic Information File
    cml         # Chemical Markup Language
    gjf         # Gaussian input data file
    gromacs     # GROMACS file format
    hyperchem   # HyperChem file format
    jme         # Java Molecule Editor format
    maestro     # Schroedinger MacroModel structure file format
    mol         # Symyx molecule file
    mol2        # Tripos Sybyl MOL2 format
    pdb         # Protein Data Bank
    sdf         # 2D formatted structure data files
    sdf3000     # Symyx Structure Data Format 3000
    sln         # SYBYL Line Notation
    xyz         # xyz file format
.. _peptidebuilder:

Peptide Builder API
===================

This page introduces the functionality of the peptide molecule object and provides a deeper look into what it can accomplish in the future.


The PeptideBuilder Class
-------------------------

The peptide molecule class is the preliminary peptide builder work before you pass it into Cocktail Shaker. It will build
the peptides smiles for you (woo!)

You instantiate a ``PeptideBuilder`` object by parsing in the length of peptide marked by how many amino acid sides you
would like. Peptide Molecule will handle any SMARTS and validation under the hood so you will not have too. A simple example
is detailed below:

>>> from cocktail_shaker import PeptideBuilder
>>> peptide_molecule = PeptideBuilder(2)
>>> print (peptide_molecule)
>>> NCC(NC([*:1])C(NC([*:2])C(O)=O)=O)=O

Here we have two slots ready for our combinations in the cocktail. Alternatively, you can include the ``include_proline``
function on the n-terminus that will provide a SMILES string like so.

>>> from cocktail_shaker import PeptideBuilder
>>> peptide_molecule = PeptideBuilder(2, include_proline=True)
>>> print (peptide_molecule)
>>> N2CCCC2C(NC([*:1])C(NC([*:2])C(O)=O)=O)=O

Generation of a circular peptide

>>> from cocktail_shaker import PeptideBuilder
>>> peptide_molecule = PeptideBuilder(3, circular=True)
>>> print (peptide_molecule)
>>> O=C1C([*:1])NC(=O)C([*:2])NC(=O)C([*:3])N1

.. attribute:: length_of_peptide

  Length of the peptide the user would like generated

.. attribute:: include_proline (optional)

  Whether the user would like to include proline on the N-Terminus. Defaults to False.

.. attribute:: circular (optional)

  If you would like the peptide to be circular or not... _gettingstarted:

Getting Started
===============

This page instructs you on how to get started with cocktail-shaker. To begin, make sure you have performed
:ref:`installing cocktail_shaker <install>`.

Basic Usage
-----------

The simplest way to use cocktail-shaker is to create a Peptide Builder object to generate your "base" peptide string,
then to create cocktail object passing the peptide backbone and with the ``shake`` function you can generate your combinations:

    >>> from cocktail_shaker import Cocktail
    >>> from cocktail_shaker import PeptideBuilder
    >>> peptide_backbone = PeptideBuilder(2)
    >>> cocktail = Cocktail(
    >>>                     peptide_backbone,
    >>>                     ligand_library = ['Br', 'I']
    >>>                     )
    >>> combinations = cocktail.shake()
    >>> print (combinations)
    >>> ['NC(Br)C(=O)NC(I)C(=O)NCC(=O)O', 'NC(I)C(=O)NC(Br)C(=O)NCC(=O)O']

Write the new compounds into an SDF file:

    >>> FileWriter('new_compounds', new_compounds, 'sdf')

In this example, we have taken one SMILES string and expanded the compounds into a variety of variations into one SDF file.

Cocktail Shaker uses the PeptideBuilder class to generate base peptide backbone strings that can be then passed as into the
cocktail shaker object. Alternatively, a user can generate their own peptide string as long as it conforms with the cocktail shaker
requirements. See the peptide builder documentation for more information.

More Examples Peptide Builder
-----------------------------

Generation of a circular peptide

>>> from cocktail_shaker import PeptideBuilder
>>> peptide_molecule = PeptideBuilder(3, circular=True)
>>> print (peptide_molecule)
>>> O=C1C([*:1])NC(=O)C([*:2])NC(=O)C([*:3])N1

Using the stereoisomer function with Cocktail Shaker

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(1)
>>> cocktail = Cocktail(
>>>     peptide_backbone,
>>>     ligand_library = ['Br'],
>>>     enable_isomers = True
>>> )
>>> combinations = cocktail.shake()
>>> print (combinations)
>>> ['N[C@H](Br)C(=O)NCC(=O)O', 'N[C@@H](Br)C(=O)NCC(=O)O']


More Examples of Cocktail Shaker
--------------------------------

Using the include amino acid function

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(1)
>>> cocktail = Cocktail(
>>>     peptide_backbone,
>>>     include_amino_acids = True
>>> )
>>> combinations = cocktail.shake()
>>> print (combinations)
>>> ['NCCCCC(N)C(=O)NCC(=O)O', 'CC(O)C(N)C(=O)NCC(=O)O', 'NC(Cc1ccc(O)cc1)C(=O)NCC(=O)O', 'NC(=O)CCC(N)C(=O)NCC(=O)O',
>>>  'CCC(C)C(N)C(=O)NCC(=O)O', 'NC(CS)C(=O)NCC(=O)O', 'NC(CC(=O)O)C(=O)NCC(=O)O', 'N=C(N)NCCCC(N)C(=O)NCC(=O)O',
>>>  'NC(Cc1c[nH]cn1)C(=O)NCC(=O)O', 'NC(Cc1ccccc1)C(=O)NCC(=O)O', 'CC(N)C(=O)NCC(=O)O', 'NC(CCC(=O)O)C(=O)NCC(=O)O',
>>>  '[H]C(N)C(=O)NCC(=O)O', 'CSCCC(N)C(=O)NCC(=O)O', 'NC(CO)C(=O)NCC(=O)O', 'CC(C)C(N)C(=O)NCC(=O)O',
>>>  'NC(CCc1c[nH]c2ccccc12)C(=O)NCC(=O)O', 'CC(C)CC(N)C(=O)NCC(=O)O']

Using the Cocktail Shaker to generate a library of halogens & single atoms with then converting into a Pandas DataFrame
Note to use this example you will need to install an extra dependency of 'tables' for handling h4 data.

>>> from peptide_builder import PeptideBuilder
>>> from functional_group_enumerator import Cocktail
>>> import pandas as pd
>>>
>>>
>>> peptide_molecule = PeptideBuilder(length_of_peptide=7)
>>> cocktail = Cocktail(peptide_backbone=peptide_molecule,
>>>                     ligand_library=["Br", "Cl", "I", "F", "O", "N", "C"],
>>>                     include_amino_acids=False,
>>>                     enable_isomers=False)
>>> molecules = cocktail.shake()
>>>
>>> dataframe = pd.DataFrame(molecules, columns=["Smiles"])
>>> dataframe.to_hdf('data.h5', key='s', mode='w').. _install:

Installation
============

cocktail-shaker supports Python versions 3.3+. Pyhton's required dependencies, most notably RDKit, must be installed.

Dependencies can be located in either `setup.py`, `requirements.txt`, or a pdf generated graph of the cocktail shaker dependency.

.. _`Dependency Graph`: https://github.com/Sulstice/Cocktail-Shaker/dependencies_cocktail_shaker.pdf

Option 1 (Recommended): Use Pip 
-------------------------------

The easiest way to install cocktail-shaker is using pip::

    pip install cocktail_shaker

This will download the latest version of cocktail-shaker and place it in your `site-packages` folder so it is automatically
available to all your Python scripts.

If you do not have pip installed yet, you can `install it using get-pip.py`_::

       curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
       python get-pip.py

Option 2: Download the Latest Release
-------------------------------------

Alternatively, you can get cocktail-shaker by manually `download the latest release`_ and installing it yourself::

    tar -xzvf cocktail_shaker-1.1.6tar.gz
    cd cocktail_shaker-1.0.0
    python setup.py install

The setup.py command will install cocktail-shaker in your `site-packages` folder so it is automatically available to all your
Python scripts.

Option 3: Clone the Repository
------------------------------

Lastly, the latest development version of cocktail-shaker is always `available on GitHub`_. The version on GitHub is not guaranteed to be stable, but may include new features that have not yet been released. 

Simply clone the repository and install as usual::

    git clone https://github.com/Sulstice/Cocktail-Shaker.git
    cd Cocktail-Shaker
    python setup.py install

.. _`install it using get-pip.py`: http://www.pip-installer.org/en/latest/installing.html
.. _`download the latest release`: https://github.com/mcs07/cocktail_shaker/releases/
.. _`available on GitHub`: https://github.com/Sulstice/Cocktail-Shaker
.. _functionalgroups:

Functional Groups
=================

Below is a list of the functional groups that cocktail-shaker currently supports.

+------------------+--------------+---------------+--------------------------------------------------------------+
| Functional Group | Class        | SMILES        | Smarts                                                       |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Bromine          | Halogens     | Br            | [Br]                                                         |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Chlorine         | Halogens     | Cl            | [Cl]                                                         |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Fluorine         | Halogens     | F             | [F]                                                          |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Acyl Bromide     | Acyl Halides | C(=O)Br       | [CX3](=[OX1])[Br]                                            |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Acyl Chloride    | Acyl Halides | C(=O)Cl       | [CX3](=[OX1])[Cl]                                            |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Acyl Fluoride    | Acyl Halides | C(=O)F        | [CX3](=[OX1])[F]                                             |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Acyl Iodide      | Acyl Halides | C(=O)I        | [CX3](=[OX1])[I]                                             |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Primary Alcohol  | Alcohols     | O             | [OX2H]                                                       |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Ketone           | Ketones      | C(=O)OC       | [CX3]=[OX1]                                                  |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Carboxylic Acid  | Acids        | C(=O)O        | [CX3](=O)[OX1H0-,OX2H1]                                      |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Acetic Anhydride | Anhydrides   | CC(=O)OC(=O)C | [CX3](=[OX1])[OX2][CX3](=[OX1])                              |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Primary Amine    | Amines       | N             | [NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]                       |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Secondary Amine  | Amines       | NC            | [NX3;H2,H1;!$(NC=O)]                                         |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Enamine          | Amines       | N             | [NX3][CX3]=[CX3]                                             |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Amide            | Amides       | C(=O)N        | [NX3][CX3](=[OX1])[#6]                                       |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Nitro            | Nitros       | [N+](=O)[O-]  | [$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]                      |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Sulfoxide        | Sulfoxides   | S(=O)(=O)     | [$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])] |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Ether            | Ethers       | COC           | [OD2]([#6])[#6]                                              |
+------------------+--------------+---------------+--------------------------------------------------------------+
| Azide            | Azides       | C([N-][N+]#N) | [$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]          |
+------------------+--------------+---------------+--------------------------------------------------------------+

If you would like to add to this list, the functional group information is stored in the **datasource** directory under
R_Groups.yml.

Each smart pattern recognition is tested thoroughly before moving into the cocktail-shaker package, so please PR your
test as well. You can follow the schema outlined in the function **test_primary_finding_r_groups**.
.. _filehandling:

File API Documentation
======================

This page provides an overview of how cocktail-shaker operates reading and writing files. To see
the files that cocktail-shaker supports, refer to :ref:`file formats <fileformats>`.

The FileParser Module
---------------------

You instantiate a ``FileParser``
by providing the exact path to a file with the extension included.
cocktail-shaker is smart enough to detect the file extension and allocate its specific parsing.
If the file being parsed is *not* supported, then ``FileNotSupportedError`` will be raised instead.

>>> from cocktail_shaker import FileParser
>>> molecules = FileParser('compounds.sdf')
>>> print (molecules)
>>> ['c1cc(CCCO)ccc1']

FileParser will then return a list of SMILES.

.. attribute:: file_path

   The path to a file


The FileWriter Module
---------------------

The FileWriter module uses the NIH web cactus resolver under the hood to generate the 2D or 3D data. By issuing requests
of the compound in SMILES and returning a XML formatted data struct of the 2D/3D coordinates in whatever chemical format
indicated by the user. Unfortunately, the cactus resolver can only process one api call at a time and not bulk requests.
This means that sending combinations worth of 1000+ -> 1000+ api requests might not be in your best interest and
eventually the server will deny any requests.

This feature is still in development so use at your own discretion.


You instantiate a ``FileWriter``
by providing the path to the file, compounds to be written, and the extension you would like the files in.
If the file being written is *not* supported, then ``FileNotSupportedError`` will be raised instead.

>>> from cocktail_shaker import Cocktail, FileWriter
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> FileWriter('new_compounds', combinations, 'sdf')
Generates an SDF File....

However, if you would like to generate the files into separate files, then you can pass in the fragmentation parameter.

>>> from cocktail_shaker import Cocktail, FileWriter
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> FileWriter('new_compounds', combinations, 'sdf', fragmentation=2)
Generates 2 SDF Files....

.. attribute:: name

   The path to a file

.. attribute:: molecules

   List of RDKit molecules you would like to write.

.. attribute:: option

   The extension of the file you would like to write

.. attribute:: fragmentation (optional)

   How many files you would like to produce.
.. _aminoacids:

Amino Acids
===========

Below is a list of the amino acids that cocktail-shaker currently supports.

+------------------+---------------------------+
| Amino Acid       | SMILES                    |
+------------------+---------------------------+
| Alanine          | C                         |
+------------------+---------------------------+
| Arginine         | CCCCNC(N)=N               |
+------------------+---------------------------+
| Asparagine       | CCC(N)=O                  |
+------------------+---------------------------+
| Aspartic Acid    | CC(O)=O                   |
+------------------+---------------------------+
| Cysteine         | CS                        |
+------------------+---------------------------+
| Glutamic Acid    | CCC(O)=O                  |
+------------------+---------------------------+
| Glutamine        | CCC(N)=O                  |
+------------------+---------------------------+
| Glycine          | [H]                       |
+------------------+---------------------------+
| Histidine        | CC1=CNC=N1                |
+------------------+---------------------------+
| Isoleucine       | C(CC)([H])C               |
+------------------+---------------------------+
| Leucine          | CC(C)C                    |
+------------------+---------------------------+
| Lysine           | CCCCN                     |
+------------------+---------------------------+
| Methionine       | CCSC                      |
+------------------+---------------------------+
| PhenylAlanine    | CC1=CC=CC=C1              |
+------------------+---------------------------+
| Proline          | C2CCCN2                   |
+------------------+---------------------------+
| Serine           | CO                        |
+------------------+---------------------------+
| Threonine        | C(C)([H])O                |
+------------------+---------------------------+
| Tryptophan       | CCC1=CNC2=C1C=CC=C2       |
+------------------+---------------------------+
| Tyrosine         | CC1=CC=C(O)C=C1           |
+------------------+---------------------------+
| Valine            | C(C)C                    |
+------------------+---------------------------+.. _cocktail:

Cocktail API
============

This page introduces the functionality of the cocktail object and provides a deeper look into what it can accomplish in the future.


The GlobalChem Class
--------------------

.. attribute:: amino_acid_side_chains

  The peptide backbone string generated from PeptideBuilder

.. attribute:: ligand_library

  The list of the ligands you would like installed on the peptide. It can be of any order.

.. attribute:: enable_isomers (optional)

  Include stereochemistry and stereoisomers in the the results. Defaults to False.

.. attribute:: include_amino_acids (optional)

  Include the natural amino acids except for Proline (TBD), list of smiles found in :ref:`amino acids <aminoacids>`.
  Defaults to False.
