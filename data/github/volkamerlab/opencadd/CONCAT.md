# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age,
body size, disability, ethnicity, gender identity and expression, level of
experience, nationality, personal appearance, race, religion, or sexual
identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Moreover, project maintainers will strive to offer feedback and advice to
ensure quality and consistency of contributions to the code.  Contributions
from outside the group of project maintainers are strongly welcomed but the
final decision as to whether commits are merged into the codebase rests with
the team of project maintainers.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an
appointed representative at an online or offline event. Representation of a
project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at 'jaime.rodriguez@charite.de'. The project team will
review and investigate all complaints, and will respond in a way that it deems
appropriate to the circumstances. The project team is obligated to maintain
confidentiality with regard to the reporter of an incident. Further details of
specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.4, available at
[http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# OpenCADD

[//]: # "Badges"

[![GH Actions Status](https://github.com/volkamerlab/opencadd/workflows/CI/badge.svg)](https://github.com/volkamerlab/opencadd/actions?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/volkamerlab/opencadd/branch/master/graph/badge.svg)](https://codecov.io/gh/volkamerlab/opencadd/branch/master)
[![Documentation Status](https://readthedocs.org/projects/opencadd/badge/?version=latest)](https://opencadd.readthedocs.io)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/opencadd.svg)](https://anaconda.org/conda-forge/opencadd)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03951/status.svg)](https://doi.org/10.21105/joss.03951) [![DOI](https://zenodo.org/badge/237947037.svg)](https://zenodo.org/badge/latestdoi/237947037)

![OpenCADD](/docs/_static/opencadd.png)

A Python library for structural cheminformatics.

## Overview

> Some modules of this library are still in early stages of development as indicated below.

- `databases.klifs`: utilities to query the KLIFS database, offline or online.
- :construction: `io`: read and write molecules from/to files.
- :construction: `structure.pocket`: identification and analysis of protein (sub)pockets.
- :construction: `structure.superposition` (formerly `superposer`): superimpose macromolecules using sequence and structural information.

## Documentation

The documentation is available [here](https://opencadd.readthedocs.io/en/latest/).

## Citation

If you are using the OpenCADD-KLIFS module, please cite our [JOSS publication](https://joss.theoj.org/papers/10.21105/joss.03951):

```
@article{Sydow2022,
  doi = {10.21105/joss.03951},
  url = {https://doi.org/10.21105/joss.03951},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {70},
  pages = {3951},
  author = {Dominique Sydow and Jaime Rodríguez-Guerra and Andrea Volkamer},
  title = {OpenCADD-KLIFS: A Python package to fetch kinase data from the KLIFS database},
  journal = {Journal of Open Source Software}
}
```

If you are using other modules of the OpenCADD package, please cite our [Zenodo entry](https://doi.org/10.5281/zenodo.5653542).

## License

`opencadd` is free software and is licensed under the MIT license. Copyright (c) 2020, Volkamer Lab

## Authors

`opencadd` is the cumulative work of several members of the [Volkamer Lab](https://volkamerlab.org), as well as contributions from students that have participated in our lab. In no particular order:

- Jaime Rodríguez-Guerra, PhD
- Dominique Sydow
- Dennis Köser, Annie Pham, Enes Kurnaz, Julian Pipart (structural superposition, 2020)

# Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
## Description
Provide a brief description of the PR's purpose here.

## Todos
Notable points that this PR has either accomplished or will accomplish.
  - [ ] TODO 1

## Questions
- [ ] Question1

## Status
- [ ] Ready to go# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into this opencadd.

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](https://help.github.com/articles/fork-a-repo/) this repository on GitHub.
* On your local machine,
  [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
  the repository.

## Making Changes

* Add some really awesome code to your local fork.  It's usually a [good
  idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
  to make changes on a
  [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
  with the branch name relating to the feature you are going to add.
* When you are ready for others to examine and comment on your new feature,
  navigate to your fork of opencadd on GitHub and open a [pull
  request](https://help.github.com/articles/using-pull-requests/) (PR). Note that
  after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.  Each commit added to the PR will be validated for
  mergability, compilation and test suite compliance; the results of these tests
  will be visible on the PR page.
* If you're providing a new feature, you must add test cases and documentation.
* When the code is ready to go, make sure you run the test suite using pytest.
* When you're ready to be considered for merging, check the "Ready to go"
  box on the PR page to let the opencadd devs know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration returns checkmarks,
  and multiple core developers give "Approved" reviews.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
# Compiling opencadd's Documentation

The docs for this project are built with [Sphinx](http://www.sphinx-doc.org/en/master/).
To compile the docs, first ensure that Sphinx and the ReadTheDocs theme are installed.


```bash
conda install sphinx sphinx_rtd_theme
```


Once installed, you can use the `Makefile` in this directory to compile static HTML pages by
```bash
make html
```

The compiled docs will be in the `_build` directory and can be viewed by opening `index.html` (which may itself
be inside a directory called `html/` depending on what version of Sphinx is installed).# Static Doc Directory

Add any paths that contain custom static files (such as style sheets) here,
relative to the `conf.py` file's directory. 
They are copied after the builtin static files,
so a file named "default.css" will overwrite the builtin "default.css".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
templates_path = ['_static']
```

## Examples of file to add to this directory
* Custom Cascading Style Sheets
* Custom JavaScript code
* Static logo images
# Templates Doc Directory

Add any paths that contain templates here, relative to  
the `conf.py` file's directory.
They are copied after the builtin template files,
so a file named "page.html" will overwrite the builtin "page.html".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
html_static_path = ['_templates']
```

## Examples of file to add to this directory
* HTML extensions of stock pages like `page.html` or `layout.html`
# Exploratory notebooks

This folder contains exploratory notebooks that may use the opencadd package and/or explore datasets like the KLIFS database.
# Data

- `2itz.pdb`
  - Example _pdb_ file (complex)
  - Taken from the PDB: https://www.rcsb.org/structure/2ITZ
- `2itz_chainA_protein.mol2`
  - Example _mol2_ file with 10 columns (protein)
  - Taken from KLIFS: https://klifs.vu-compmedchem.nl/details.php?structure_id=789
- `2itz_chainA_ligand.mol2`
  - Example _mol2_ file with 10 columns (ligand)
  - Taken from KLIFS: https://klifs.vu-compmedchem.nl/details.php?structure_id=789
- `2itz_protein.mol2`
  - Example _mol2_ file with 9 columns (protein)
  - Taken from the scPDB: http://bioinfo-pharma.u-strasbg.fr/scPDB/SITE=2itz_1
# Sample Package Data

This directory contains sample additional data you may want to include with your package.
This is a place where non-code related additional information (such as data files, molecular structures,  etc.) can 
go that you want to ship alongside your code.

Please note that it is not recommended to place large files in your git directory. If your project requires files larger
than a few megabytes in size it is recommended to host these files elsewhere. This is especially true for binary files
as the `git` structure is unable to correctly take updates to these files and will store a complete copy of every version
in your `git` history which can quickly add up. As a note most `git` hosting services like GitHub have a 1 GB per repository
cap.

## Including package data

Modify your package's `setup.py` file and the `setup()` command. Include the 
[`package_data`](http://setuptools.readthedocs.io/en/latest/setuptools.html#basic-use) keyword and point it at the 
correct files.

## Manifest

* `look_and_say.dat`: first entries of the "Look and Say" integer series, sequence [A005150](https://oeis.org/A005150)
* `klifs_ids.YYYYMMDD.csv.zip`: KLIFS structure and kinase IDs that are available at wwww.klifs.net by the time of YYYYMMDD.
---
title: 'OpenCADD-KLIFS: A Python package to fetch kinase data from the KLIFS database'
tags:
  - Python
  - KLIFS
  - kinase
authors:
  - name: Dominique Sydow^[corresponding author]
    orcid: 0000-0003-4205-8705
    affiliation: 1
  - name: Jaime Rodríguez-Guerra
    orcid: 0000-0001-8974-1566
    affiliation: 1
  - name: Andrea Volkamer^[corresponding author]
    affiliation: 1
    orcid: 0000-0002-3760-580X
affiliations:
 - name: _In Silico_ Toxicology and Structural Bioinformatics, Institute of Physiology, Charité – Universitätsmedizin Berlin, corporate member of Freie Universität Berlin and Humboldt-Universität zu Berlin, Augustenburger Platz 1, 13353 Berlin, Germany
   index: 1
date: 27 October 2021
bibliography: paper.bib
---

# Summary

Protein kinases are involved in most aspects of cell life due to their role in signal transduction. Dysregulated kinases can cause severe diseases such as cancer, inflammation, and neurodegeneration, which has made them a frequent target in drug discovery for the last decades [@Cohen:2021].
The immense research on kinases has led to an increasing amount of kinase resources [@Kooistra:2017].
Among them is the KLIFS database, which focuses on storing and analyzing structural data on kinases and interacting ligands [@Kanev:2021].
The OpenCADD-KLIFS Python module offers a convenient integration of the KLIFS data into workflows to facilitate computational kinase research.

[OpenCADD-KLIFS](https://opencadd.readthedocs.io/en/latest/databases_klifs.html) (``opencadd.databases.klifs``) is part of the [OpenCADD](https://opencadd.readthedocs.io/) package, a collection of Python modules for structural cheminformatics.

# Statement of need

The KLIFS resource [@Kanev:2021] contains information about kinases, structures, ligands, interaction fingerprints, and bioactivities. 
KLIFS thereby focuses especially on the ATP binding site, defined as a set of 85 residues and aligned across all structures using a multiple sequence alignment [@vanLinden:2014].
Fetching, filtering, and integrating the KLIFS content on a larger scale into Python-based pipelines is currently not straight-forward, especially for users without a background in online queries. Furthermore, switching between data queries from a _local_ KLIFS download and the _remote_ KLIFS database is not readily possible.

OpenCADD-KLIFS is aimed at current and future users of the KLIFS database who seek to 
integrate kinase resources into Python-based research projects.
With OpenCADD-KLIFS, KLIFS data can be queried either locally from a KLIFS download or remotely from the KLIFS webserver. 
The presented module provides identical APIs for the remote and local queries and streamlines all output into 
standardized Pandas DataFrames [@pandas] to allow for easy and quick downstream data analyses (\autoref{fig:opencadd_klifs_toc}). This Pandas-focused setup is ideal if you work with Jupyter notebooks [@Kluyver:2016]. 

![OpenCADD-KLIFS fetches KLIFS data [@Kanev:2021] offline from a local KLIFS download or online from the KLIFS database and formats the output as user-friendly Pandas DataFrames [@pandas].\label{fig:opencadd_klifs_toc}](opencadd_klifs_toc.png)

# State of the field

The KLIFS database is unique in the structure-based kinase field in terms of integrating and annotating different data resources in a kinase- and pocket-focused manner. Kinases, structures, and ligands have unique identifiers in KLIFS, which makes it possible to fetch and filter cross-referenced information for a query kinase, structure, or ligand.

- Kinase structures are fetched from the PDB, split by chains and alternate models, annotated with the KLIFS pocket of 85 residues, and aligned across the fully structurally covered kinome.
- Kinase-ligand interactions seen in experimental structures are annotated for the 85 pocket residues in the form of the KLIFS interaction fingerprint (KLIFS IFP).
- Bioactivity data measured against kinases are fetched from ChEMBL [@Mendez:2018] and linked to kinases, structures, and ligands available in KLIFS.
- Kinase inhibitor metadata are fetched from the PKIDB [@Carles:2018] and linked to co-crystallized ligands available in KLIFS.

The KLIFS data integrations and annotations can be accessed in different ways, which are all open source:

- Manually via the [KLIFS website](https://klifs.net/) interface: This mode is preferable when searching for information on a specific structure or smaller set of structures.
- Automated via the [KLIFS KNIME](https://github.com/3D-e-Chem/knime-klifs) nodes [@McGuire:2017; @Kooistra:2018]: This mode is extremely useful if the users' projects are embedded in KNIME workflows; programming is not needed.
- Programmatically using the REST API and KLIFS OpenAPI specifications: This mode is needed for users who seek to perform larger scale queries or to integrate different queries into programmatic workflows. In the following, we will discuss this mode in context of Python-based projects and explain how OpenCADD-KLIFS improves the user experience.

The KLIFS database offers standardized URL schemes (REST API), which allows users to query data by defined URLs, using e.g., the Python package requests [@requests]. Instead of writing customized scripts to generate such KLIFS URLs, the KLIFS OpenAPI specifications, a document that defines the KLIFS REST API scheme, can be used to generate a Python client, using e.g., the Python package bravado [@bravado]. This client offers a Python API to send requests and receive responses.
This setup is already extremely useful, however, it has a few drawbacks: the setup is technical; the output is not easily readable for humans and not ready for immediate downstream integrations, requiring similar but not identical reformatting functions for different query results; and switching from remote requests to local KLIFS download queries is not possible. Facilitating and streamlining these tasks is the purpose of OpenCADD-KLIFS as discussed in more detail in the next section.

# Key Features

The KLIFS database offers a REST API compliant with the OpenAPI specification [@klifsswagger]. Our module OpenCADD-KLIFS uses bravado to dynamically generate a Python client based on the OpenAPI definitions and adds wrappers to enable the following functionalities:

- A session is set up automatically, which allows access to various KLIFS *data sources* by different *identifiers* with the API ``session.data_source.by_identifier``. *Data sources* currently include kinases, structures and annotated conformations, modified residues, pockets, ligands, drugs, and bioactivities; *identifiers* refer to kinase names, PDB IDs, KLIFS IDs, and more. For example, ``session.structures.by_kinase_name`` fetches information on all structures for a query kinase.
- The same API is used for local and remote sessions, i.e., interacting with data from a KLIFS download folder and from the KLIFS website, respectively.
- The returned data follows the same schema regardless of the session type (local/remote); all results obtained with bravado are formatted as Pandas DataFrames with standardized column names, data types, and handling of missing data.
- Files with the structural 3D coordinates deposited on KLIFS include full complexes or selections such as proteins, pockets, ligands, and more. These files can be downloaded to disc or loaded via biopandas [@Raschka:2017] or RDKit [@rdkit].

OpenCADD-KLIFS is especially convenient whenever users are interested in multiple or more complex queries such as "fetching all structures for the kinase EGFR in the DFG-in conformation" or "fetching the measured bioactivity profiles for all ligands that are structurally resolved in complex with EGFR". Formatting the output as DataFrames facilitates subsequent filtering steps and DataFrame merges in case multiple KLIFS datasets need to be combined.

OpenCADD-KLIFS is currently used in several projects from the Volkamer Lab [@volkamerlab] including TeachOpenCADD [@teachopencadd], OpenCADD-pocket [@opencadd_pocket], KiSSim [@kissim], KinoML [@kinoml], and PLIPify [@plipify].
For example, OpenCADD-KLIFS is applied in a [TeachOpenCADD tutorial](https://projects.volkamerlab.org/teachopencadd/talktorials/T012_query_klifs.html) to demonstrate how to fetch all kinase-ligand interaction profiles for all available EGFR kinase structures to visualize the per-residue interaction types and frequencies with only a few lines of code.

# Acknowledgements

We thank the whole KLIFS team for providing such a great kinase resource with an easy-to-use API and especially Albert Kooistra for his help with questions and wishes regarding the KLIFS database. 
We thank David Schaller for his feedback on the OpenCADD-KLIFS module.
We acknowledge the contributors involved in software programs and packages used by OpenCADD-KLIFS, such as bravado, RDKit, Pandas, Jupyter, and Pytest, and Sphinx. 

# References
---
title: 'Superposer: A Python library to use different methods for structural alignment.'
tags:
- python
- bioinformatics
authors:
- name: Julian Pipart
  orcid: 0000-0002-1089-5931
- name: Jaime Rodríguez-Guerra
  orcid: 0000-0001-8974-1566
- name: Andrea Volkamer
- name: Dennis Köser
- name: Annie Pham
- name: Enes Kurnaz
date: 29 April 2020
---

# Superposer

A Python library to use different methods for structural alignment.

Contributors:

- Jaime Rodríguez-Guerra
- Dennis Köser
- Enes Kurnaz
- Annie Pham
- Julian Pipart
- Andrea Volkamer

# Abstract

**Motivation**: There are a lot of different methods to align structures with each other. But most of the time you do not know which method is the best for a specific situation, or you want to build a python script around these methods. Some methods are not even implemented in python. That makes it difficult to import these methods in custom scripts or software projects. In order to use these different methods you need to install them separately. They also return output which is not consistent over the different methods.
**Results**: superposer has been developed to bring multiple methods for structural alignment together and provide consistency in the output as well as input for any Python 3.7 or Python 3.8 interpreter. It provides the possibility of interactive programming and makes it easier to use third-party software for structural alignment.
**Availability and implementation**: superposer is MIT-licensed and available at https://github.com/volkamerlab/superposer
**Contact**: julian.pipart@fu-berlin.de or jaime.rodriguez@charite.de


# 1 Introduction

Python is a very powerful programming language. No matter if it is used for projects in big companies or for small scripts. It is especially useful for solving scientific or computational problems because of its variety of different packages to use. This includes problems of structural bioinformatics and chemistry. There is already a lot of software to solve specific problems and most of them use the same idea but a different approach. They often accept not the same input and provide different outputs. This can be really hindering for the research or the process of software development.


# 2 Materials and methods

superposer is a Python package which brings multiple methods for structural alignment together. The user can specify the method to use and does not have to worry about the dependencies and different input parameters for different methods. To ensure consistency across all methods superposer uses Atomium ([Ireland, S. M., & Martin, A. C., 2020](#references)). Atomium is able to catch the structures the user gives as input from the RCSB. With the help of Atomium superposer first transforms the structure to an atomium.model and passes this model on to the method of choice. All methods calculate a transformation matrix with their own approaches. This matrix is applied to the original model. By that there will be no metadata lost. The output of superposer is also an atomium.model for all methods implemented. The model can be saved as a Protein Data Bank file (PDB) or similar filetypes. Because most methods do not know what an atomium.model is, they save the model as a temporary PDB file and parse it. This makes the package very easy to use in development of software. It also assures the ability to be further extended. The correctness is ensured by the Github CI workflow. If a method is written in a different programming language and is open source, it is wrapped in Python. Other methods, like Chimera matchmaker ([Pettersen, E. F., Goddard, T. D., Huang, C. C., Couch, G. S., Greenblatt, D. M., Meng, E. C., & Ferrin, T. E. 2004](#references)) has been reimplemented in this library.

# 3 Results

superposer is a command line application with support for the latest Linux distributions (it is tested for Python 3.7 and Python 3.8 on Ubuntu 18.04) and Mac OS X (it is tested for Python 3.7 and Python 3.8 on Mac OS X 10.15.4). By installing superposer, all dependencies needed are installed with it. Only Python 3.7 or later is required before installation. The installation is made possible through conda with the command: `conda install superposer`. It is the same command for both platforms. A great advantage is, that superposer uses Python 3. That guarantees good compatibility with most modern software written in Python.

**3.1 Usage**

Running superposer in the terminal can be performed with the command: `superposer TARGET_STRUCTURE STRUCTURE_TO_ALIGN --method METHOD_OF_CHOICE --method-options METHOD_OPTIONS`. superposer features the methods MMLigner ([Collier, J. H., Allison, L., Lesk, A. M., Stuckey, P. J., Garcia de la Banda, M., & Konagurthu, A. S., 2017](#references)), Theseus ([Theobald, D. L., & Wuttke, D. S., 2006](#references)) and Chimera's MatchMaker. It is also possible to perform multiple pairwise alignments by adding more structures to align to the command above. All structures will therefore be aligned against the target structure. The output in the terminal is the RMSD.

**3.2 Demonstrating the Benefits of 'superposer'**

superposer brings different methods for structural alignment together, so that it is not necessary to have all the stand alone programs installed. The methods also have different approaches.
For the more traditional approach, one can use the method matchmaker from chimera

```
INPUT: superposer 2GZ9 5R8T --method=matchmaker
OUTPUT:
  _____ _                   _                   _
  / ____| |                 | |                 | |
 | (___ | |_ _ __ _   _  ___| |_ _   _ _ __ __ _| |
  \___ \| __| '__| | | |/ __| __| | | | '__/ _` | |
  ____) | |_| |  | |_| | (__| |_| |_| | | | (_| | |
 |_____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_|
     /\   | (_)                                | |
    /  \  | |_  __ _ _ __  _ __ ___   ___ _ __ | |_
   / /\ \ | | |/ _` | '_ \| '_ ` _ \ / _ \ '_ \| __|
  / ____ \| | | (_| | | | | | | | | |  __/ | | | |_
 /_/    \_\_|_|\__, |_| |_|_| |_| |_|\___|_| |_|\__|
                __/ |
               |___/

 v0+untagged.91.g31533ad · Brought to you by @volkamerlab

RMSD for alignment between `2GZ9` and `5R8T` is 0.4Å
```

If one prefer a statistical approach, one can use for example Theseus.

```
INPUT: superposer 2GZ9 5R8T --method=theseus
OUTPUT:
  _____ _                   _                   _
  / ____| |                 | |                 | |
 | (___ | |_ _ __ _   _  ___| |_ _   _ _ __ __ _| |
  \___ \| __| '__| | | |/ __| __| | | | '__/ _` | |
  ____) | |_| |  | |_| | (__| |_| |_| | | | (_| | |
 |_____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_|
     /\   | (_)                                | |
    /  \  | |_  __ _ _ __  _ __ ___   ___ _ __ | |_
   / /\ \ | | |/ _` | '_ \| '_ ` _ \ / _ \ '_ \| __|
  / ____ \| | | (_| | | | | | | | | |  __/ | | | |_
 /_/    \_\_|_|\__, |_| |_|_| |_| |_|\___|_| |_|\__|
                __/ |
               |___/

 v0+untagged.91.g31533ad · Brought to you by @volkamerlab

RMSD for alignment between `2GZ9` and `5R8T` is 1.7Å
```

**3.3 Current limitations**

For now, the only method, that supports a benchmark is MMLigner in form of the ivalue. The other methods do not provide a benchmark yet. It is not clear if the ivalue is a good and consistent benchmark over different methods. This means it is not obvious whether an alignment is good or not, compared to the other alignments. It is also not possible to perform a multiple alignment at the moment. One can perform multiple pairwise alignments successively.

**3.4 Future work**

In the future more alignment methods will be added. Because superposer is open source, every developer can participate in making superposer better and more versatile. The MMLigner method already provides a benchmark with the ivalue. The developer team is working on providing a benchmark for every alignment method. This benchmark should be consistent over the different methods. In addition to the benchmarking the team is working on a solution to visualize the aligned structures in this library, without the need of any third-party program.

# Acknowledgments

The authors want to thank the team of Volkamerlab, and more particularly Jaime Rodríguez-Guerra and Andrea Volkamer for making this project possible and for leading the development process.

# References

Collier, J. H., Allison, L., Lesk, A. M., Stuckey, P. J., Garcia de la Banda, M., & Konagurthu, A. S. (2017). Statistical inference of protein structural alignments using information and compression. Bioinformatics, 33(7), 1005-1013.

Ireland, S. M., & Martin, A. C. (2020). atomium—a Python structure parser. Bioinformatics.

Theobald, D. L., & Wuttke, D. S. (2006). THESEUS: maximum likelihood superpositioning and analysis of macromolecular structures. Bioinformatics, 22(17), 2171-2172.

Pettersen, E. F., Goddard, T. D., Huang, C. C., Couch, G. S., Greenblatt, D. M., Meng, E. C., & Ferrin, T. E. (2004). UCSF Chimera—a visualization system for exploratory research and analysis. Journal of computational chemistry, 25(13), 1605-1612.

# Development, testing, and deployment tools

This directory contains a collection of tools for running Continuous Integration (CI) tests,
conda installation, and other development tools not directly related to the coding process.


## Manifest

### Continuous Integration

You should test your code, but do not feel compelled to use these specific programs. You also may not need Unix and
Windows testing if you only plan to deploy on specific platforms. These are just to help you get started

* `travis-ci`: Linux and OSX based testing through [Travis-CI](https://about.travis-ci.com/)
  * `before_install.sh`: Pip/Miniconda pre-package installation script for Travis
* `appveyor`: Windows based testing through [AppVeyor](https://www.appveyor.com/) (there are no files directly related to this)

### Conda Environment:

This directory contains the files to setup the Conda environment for testing purposes

* `conda-envs`: directory containing the YAML file(s) which fully describe Conda Environments, their dependencies, and those dependency provenance's
  * `test_env.yaml`: Simple test environment file with base dependencies. Channels are not specified here and therefore respect global Conda configuration

### Additional Scripts:

This directory contains OS agnostic helper scripts which don't fall in any of the previous categories
* `scripts`
  * `create_conda_env.py`: Helper program for spinning up new conda environments based on a starter file with Python Version and Env. Name command-line options


## How to contribute changes
- Clone the repository if you have write access to the main repo, fork the repository if you are a collaborator.
- Make a new branch with `git checkout -b {your branch name}`
- Make changes and test your code
- Ensure that the test environment dependencies (`conda-envs`) line up with the build and deploy dependencies (`conda-recipe/meta.yaml`)
- Push the branch to the repo (either the main or your fork) with `git push -u origin {your branch name}`
  * Note that `origin` is the default name assigned to the remote, yours may be different
- Make a PR on GitHub with your changes
- We'll review the changes and get your code into the repo after lively discussion!


## Checklist for updates
- [ ] Make sure there is an/are issue(s) opened for your specific update
- [ ] Create the PR, referencing the issue
- [ ] Debug the PR as needed until tests pass
- [ ] Tag the final, debugged version
   *  `git tag -a X.Y.Z [latest pushed commit] && git push --follow-tags`
- [ ] Get the PR merged in

## Versioneer Auto-version
[Versioneer](https://github.com/warner/python-versioneer) will automatically infer what version
is installed by looking at the `git` tags and how many commits ahead this version is. The format follows
[PEP 440](https://www.python.org/dev/peps/pep-0440/) and has the regular expression of:
```regexp
\d+.\d+.\d+(?\+\d+-[a-z0-9]+)
```
If the version of this commit is the same as a `git` tag, the installed version is the same as the tag,
e.g. `opencadd-0.1.2`, otherwise it will be appended with `+X` where `X` is the number of commits
ahead from the last tag, and then `-YYYYYY` where the `Y`'s are replaced with the `git` commit hash.
Installing
==========

.. note::

    We are assuming you have a working ``mamba`` installation in your computer. 
    If this is not the case, please refer to their `official documentation <https://mamba.readthedocs.io/en/latest/installation.html#mamba>`_. 

    If you installed ``mamba`` into an existing ``conda`` installation, also make sure that the ``conda-forge`` channel is configured by running ``conda config --add channels conda-forge``.


Install from the conda package
------------------------------

1. Create a new conda environment called ``opencadd`` with the ``opencadd`` package and all its dependencies installed::

    mamba create -n opencadd opencadd

2. Activate the new conda environment::

    conda activate opencadd

.. 3. Test that your installation works::

    superposer -h


Install from the latest development snapshot
--------------------------------------------

Install the latest development snapshot from the `GitHub repository's master branch <https://github.com/volkamerlab/opencadd>`_.


1. Create a new conda environment called ``opencadd``::

    mamba env create -f https://raw.githubusercontent.com/volkamerlab/opencadd/master/devtools/conda-envs/user_env.yaml

2. Activate the new conda environment::

    conda activate opencadd

3. Install ``opencadd`` package via pip::

    pip install https://github.com/volkamerlab/opencadd/archive/master.tar.gz

.. 4. Test that your installation works::

    superposer -h


Development version
-------------------

To install the development version of OpenCADD, you can run::

    # Install development environment (incl. packages for testing and documentation)
    mamba env create -f https://raw.githubusercontent.com/volkamerlab/opencadd/master/devtools/conda-envs/test_env.yaml -n opencadd-dev
    conda activate opencadd-dev
    
    # Download repository and install `opencadd` package in editable mode
    git clone git@github.com:volkamerlab/opencadd.git
    pip install -e opencadd

    # Let's change into the repo's folder
    cd opencadd
    
    # Run tests like this
    pytest -v opencadd/tests/

    # Before you push changes to GitHub, lint and style-check your code
    pylint opencadd
    black --check -l 99 opencadd

    # Check how the documentation renders
    cd docs
    make html

Note: If you add new dependencies to ``opencadd``, please update the 
`test <https://github.com/volkamerlab/opencadd/blob/master/devtools/conda-envs/test_env.yaml>`_ and 
`user <https://github.com/volkamerlab/opencadd/blob/master/devtools/conda-envs/user_env.yaml>`_ environment, 
and leave a note in the 
`dependencies <https://github.com/volkamerlab/opencadd/blob/master/docs/installing.rst#dependencies>`_ section.


Dependencies
------------

``opencadd`` is supported for at least Python 3.7 and needs the following modules: 

+---------------------+---------------------------+--------------------+--+--+
| Package             | Minimal supported version | Used by submodules |  |  |
+=====================+===========================+====================+==+==+
| ``pandas``          | 1.3.4                     | 1, 2, 4            |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``biopandas``       | 0.2.9                     | 1, 2               |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``bravado``         | 11.0.3                    | 1                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``rdkit``           | 2021.09.2                 | 1                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``tqdm``            | 4.62.3                    | 1                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``biopython``       | = 1.77                    | 2                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``mdanalysis``      | 2.0.0                     | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``biotite``         | 0.31.0                    | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``mmligner``        | 1.0.2                     | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``muscle``          | 3.8.1551                  | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``theseus``         | 3.3.0                     | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``matplotlib-base`` | 3.5.0                     | 4                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``nglview``         | 3.0.3                     | 4                  |  |  |
+---------------------+---------------------------+--------------------+--+--+


Some packages are only needed for a subset of the following modules: [1] ``opencadd.databases.klifs``, 
[2] ``opencadd.io``, 
[3] ``opencadd.structure.superposition``, 
[4] ``opencadd.structure.pocket``

This list of minimal supported versions is based on `this CI run <https://github.com/volkamerlab/opencadd/runs/4462667598?check_suite_focus=true#step:6:42>`_.API Documentation
=================

.. autosummary::
   :toctree: generated

   opencadd
   opencadd.structure.core
   opencadd.structure.superposition
   opencadd.structure.superposition.api
   opencadd.structure.superposition.cli
   opencadd.structure.superposition.engines.base
   opencadd.structure.superposition.engines.mda
   opencadd.structure.superposition.engines.mmligner
   opencadd.structure.superposition.engines.theseus
   opencadd.structure.superposition.sequences
   opencadd.utils

Structural superposition
========================

.. todo::

    Consider using ``opencadd superpose`` as a subcommand.

Once you have installed the package, you will have access to a ``superposer`` executable (among other utilities!). Some quick examples:

Two structures superposed with mmligner without any additional options::

    superposer 1lol 5XME
    superposer 1lol 5XME --method theseus

Input structures can be PDB identifiers, or paths. to local files. If more than two structures are provided,
the first will be considered the target structure (it will remain static). The other structure will be superposed
on top of this one, creating multiple structure alignments.

``superposer`` will create several PDB files on disk for you to check with your favourite visualization tool.
Structure: Pocket
=================

Once you have installed the package, you will have access (among others) 
to the ``opencadd.structure.pocket`` module.

This module offers a simple API to define and visualize 
(``nglview``) subpockets and regions within a pocket.

``nglview``: http://nglviewer.org/nglview/latest/

Pocket
------

The ``Pocket`` class is the main API to work with. A ``Pocket`` object is initialized from a 
protein structure file, a name, pocket residue PDB IDs and optionally pocket residue labels.
For instance, pocket residue labels can be used to label residues with alignment IDs, 
such as the KLIFS residue ID (pocket residues 1-85 are aligned across all kinases).


.. code-block:: python

    from opencadd.structure.pocket import Pocket

    pocket = Pocket.from_file(
        filepath="protein.mol2",  # Dummy file path
        pocket_residue_pdb_ids=[50, 51, ..., 198], 
        name="example kinase", 
        pocket_residue_labels=[1, 2, ..., 85]
    )


Subpockets
----------

Currently, the main focus of this module is the definition of subpockets within a pocket.
Once a ``Pocket`` object is initialized, subpockets can be added one-by-one.
The user can define each subpocket by so-called anchor residues. The centroid of all anchor
residues' CA atoms is the subpocket center and will be visualized as a subpocket sphere 
(together with the anchor residues' CA atoms as small spheres). 
Anchor residues need to be selected manually, in a way that their centroid will cover the subpocket 
center as intended.

As an example, we add here a kinase subpocket called ``"hinge"``, which will be 
visualized in magenta and whose center is calculated based on the CA atoms of residues 73, 128, and
193 (residue PDB IDs). These residue PDB IDs correspond to the KLIFS alignment IDs 16, 47, and 80.


.. code-block:: python

    pocket.add_subpocket(
        name="hinge", 
        anchor_residue_pdb_ids=[73, 128, 193], 
        color="magenta", 
        anchor_residue_labels=[16, 47, 80]  # Optionally
    )


Regions
------- 

Usually, it is helpful to highlight the pocket region, which inspired the subpocket choice, 
hence the module allows definitions of important pocket regions. 

In our example, the subpocket ``"hinge"`` is intended to target one of the most important
regions in kinases, the hinge region, where kinase inhibitors form hydrogen bonds with the kinase.
We use the hinge region residues to add a region to our ``Pocket`` object.

.. code-block:: python

    pocket.add_region(
        name="hinge region", 
        residue_pdb_ids=[127, 128, 129], 
        color="magenta", 
        residue_labels=[46, 47, 48]  # Optionally
    )


Visualize the pocket
--------------------

Now we can visualize the pocket using:

.. code-block:: python

    pocket.visualize()


Check out our tutorial to find out more!
Statement of need
================= 

The KLIFS resource [Kanev_2021]_ contains information about kinases, structures, ligands, 
interaction fingerprints, and bioactivities. 
KLIFS thereby focuses especially on the ATP binding site, defined as a set of 85 residues and 
aligned across all structures using a multiple sequence alignment [vanLinden_2014]_.
Fetching, filtering, and integrating the KLIFS content on a larger scale into Python-based 
pipelines is currently not straight-forward, especially for users without a background in 
online queries. 
Furthermore, switching between data queries from a *local* KLIFS download and 
the *remote* KLIFS database is not readily possible.

OpenCADD-KLIFS is aimed at current and future users of the KLIFS database who seek to 
integrate kinase resources into Python-based research projects.
With OpenCADD-KLIFS, KLIFS data can be queried either locally from a KLIFS download or 
remotely from the KLIFS webserver. 
The presented module provides identical APIs for the remote and local queries and 
streamlines all output into standardized Pandas DataFrames 
`Pandas <https://doi.org/10.5281/zenodo.5574486>`_  to allow for easy and quick 
downstream data analyses (Figure 1). 
This Pandas-focused setup is ideal if you work with Jupyter notebooks [Kluyver_2016]_.

.. raw:: html

   <p align="center">
   <img src="_static/opencadd_klifs_toc.png" alt="OpenCADD-KLIFS" width="600"/>
   </p>

*Figure 1*: OpenCADD-KLIFS fetches KLIFS data offline from a KLIFS download or 
online from the KLIFS database and formats the output as user-friendly Pandas DataFrames.

.. [Kanev_2021] Kanev et al., (2021),
   KLIFS: an overhaul after the first 5 years of supporting kinase research,
   Nucleic Acids Research, 
   49(D1), D562–D569, doi:10.1093/nar/gkaa895.
.. [vanLinden_2014] van Linden et al., (2014)
   KLIFS: A Knowledge-Based Structural Database To Navigate Kinase–Ligand 
   Interaction Space, 
   Journal of Medicinal Chemistry, 
   57(2), 249-277, doi:10.1021/jm400378w.
.. [Kluyver_2016] Kluyver et al., (2016),
   Jupyter Notebooks – a publishing format for reproducible computational workflows,
   In Positioning and Power in Academic Publishing: Players, Agents and Agendas. IOS Press. pp. 87-90,
   doi:10.3233/978-1-61499-649-1-87... opencadd documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OpenCADD's documentation!
=========================================================

OpenCADD is a Python package for structural cheminformatics!

.. image::
   https://github.com/volkamerlab/opencadd/workflows/CI/badge.svg
   :target: https://github.com/volkamerlab/opencadd/actions
.. image::
   https://codecov.io/gh/volkamerlab/opencadd/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/volkamerlab/opencadd/branch/master
.. image::
   https://readthedocs.org/projects/opencadd/badge/?version=latest
   :target: https://opencadd.readthedocs.io/en/latest/
.. image::
   https://img.shields.io/conda/vn/conda-forge/opencadd.svg
   :target: https://anaconda.org/conda-forge/opencadd
.. image::
   https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT

.. raw:: html

   <p align="center">
   <img src="_static/opencadd.png" alt="Subpocket-based structural fingerprint for kinase pockets" width="600"/>
   </p>

.. toctree::
   :maxdepth: 1
   :caption: User guide

   installing
   installing_opencadd_klifs

.. toctree::
   :maxdepth: 1
   :caption: IO formats

   io
   tutorials/io

.. toctree::
   :maxdepth: 1
   :caption: OpenCADD-superposition

   superposition
   tutorials/mda
   tutorials/mmligner
   tutorials/theseus

.. toctree::
   :maxdepth: 1
   :caption: OpenCADD-pocket

   structure_pocket
   tutorials/structure_pocket

.. toctree::
   :maxdepth: 1
   :caption: OpenCADD-KLIFS

   databases_klifs
   databases_klifs_statement_of_need
   tutorials/databases_klifs

.. toctree::
   :maxdepth: 1
   :caption: Developers

   developers
   api




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
IO
==

Once you have installed the package, you will have access (among others) 
to the ``opencadd.io`` module.

This module offers a simple API to load a variety of input formats in the form of a variety of
output formats.

This module is very much work-in-progress and offers so far only a first round of input/output 
options.

Check out the following tutorial to explore the API.
OpenCADD-KLIFS
================

Once you have installed the ``opencadd`` package, you will have access (among others) 
to the ``opencadd.databases.klifs`` module (OpenCADD-KLIFS).
In case you wish to install only the dependencies relevant to OpenCADD-KLIFS, please follow the installation instructions `here <https://opencadd.readthedocs.io/en/latest/installing_opencadd_klifs.html>`_.

OpenCADD-KLIFS offers a simple API to interact with data from KLIFS remotely and locally.

Find a detailed tutorial at the `TeachOpenCADD platform <https://projects.volkamerlab.org/teachopencadd/talktorials/T012_query_klifs.html>`_ on the KLIFS database and on how to apply the module OpenCADD-KLIFS to an example research question.

What is KLIFS and who created it?
---------------------------------

  KLIFS is a kinase database that dissects experimental structures of catalytic kinase domains and the way kinase inhibitors interact with them. The KLIFS structural alignment enables the comparison of all structures and ligands to each other. Moreover, the KLIFS residue numbering scheme capturing the catalytic cleft with 85 residues enables the comparison of the interaction patterns of kinase-inhibitors, for example, to identify crucial interactions determining kinase-inhibitor selectivity.

- KLIFS database: https://klifs.net (official), https://dev.klifs.net/ (developmental)
- KLIFS online service: https://klifs.net/swagger (official), https://dev.klifs.net/swagger_v2 (developmental, used here) 
- KLIFS citation: `Nucleic Acids Res. (2021), 49, D1, D562–D569 <https://academic.oup.com/nar/article/49/D1/D562/5934416>`_

What does ``opencadd.databases.klifs`` offer?
---------------------------------------------

This module allows you to access KLIFS data such as information about 
kinases, structures, structural conformations, modified residues, ligands, drugs, interaction fingerprints, and bioactivities.

On the one hand, you can query the KLIFS webserver directly. 
On the other hand, you can query your local KLIFS download.
We provide identical APIs for the remote and local queries and streamline all output into standardized ``pandas`` DataFrames for easy and quick downstream data analyses.

Work with KLIFS data from KLIFS server (remotely)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``opencadd.databases.klifs.remote`` submodule offers you to access KLIFS data from the KLIFS server.

Our API relies on the REST API and OpenAPI (formerly Swagger API) specification at https://dev.klifs.net/swagger_v2/ to dynamically generate a Python client with ``bravado``.

Example for ``opencadd``'s API to access remote data:

.. code-block:: python

    from opencadd.databases.klifs import setup_remote

    # Set up remote session
    session = setup_remote()

    # Get all kinases that are available remotely
    session.kinases.all_kinases()

    # Get kinases by kinase name
    session.kinases.by_kinase_name(["EGFR", "BRAF"])

Work with KLIFS data from disc (locally)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``opencadd.databases.klifs.local`` submodule offers you to access KLIFS data from your KLIFS download. 
In order to make use of the module's functionality, you need a KLIFS download folder ``KLIFS_download`` with the following structure (files downloaded from `KLIFS <from https://klifs.net>`_):

.. code-block:: console 

    └── KLIFS_download 
        ├── KLIFS_export.csv           # Metadata file 
        ├── overview.csv               # Metadata file 
        └── HUMAN     	               # Species name 
            ├── AAK1                   # Kinase name 
            │   ├── 4wsq_altA_chainA   # PDB ID, alternate model ID, chain ID 
            │   │   ├── complex.mol2 
            │   │   ├── ligand.mol2 
            │   │   ├── protein.mol2 
            │   │   ├── pocket.mol2 
            │   │   └── water.mol2 
            │   └── ... 
            └── ... 

Example for ``opencadd``'s API to access local data:

.. code-block:: python

    from opencadd.databases.klifs import setup_local

    # Set up local session
    session = setup_local("../../opencadd/tests/databases/data/KLIFS_download")

    # Get all kinases that are available locally
    session.kinases.all_kinases()

    # Get kinases by kinase name
    session.kinases.by_kinase_name(["EGFR", "BRAF"])


How is ``opencadd.databases.klifs`` structured?
----------------------------------------------------------

The module's structure looks like this, using the same API for both modules ``local`` and ``remote`` whenever possible:

.. code-block:: console 

    opencadd/ 
        └── databases/
            └── klifs/
                ├── api.py         # Defines the main API for local and remote sessions.
                ├── session.py     # Defines a KLIFS session.
                ├── core.py        # Defines the parent classes used in the local and remote modules.
                ├── local.py       # Defines the API for local queries.
                ├── remote.py      # Defines the API for remote queries.
                ├── schema.py      # Defines the schema for class method return values.
                ├── fields.py      # Defines the different KLIFS data fields and their names/dtypes in ``opencadd``.
                ├── utils.py       # Defines utility functions.
                └── exceptions.py  # Defines exceptions.

This structure mirrors the KLIFS OpenAPI structure in the following way to access different kinds of information both remotely and locally:

- ``kinases``  

  - Get information about kinases (groups, families, names).  
  - In KLIFS OpenAPI called ``Information``: https://dev.klifs.net/swagger_v2/#/Information

- ``ligands``  

  - Get ligand information.  
  - In KLIFS OpenAPI called ``Ligands``: https://dev.klifs.net/swagger_v2/#/Ligands

- ``structures``

  - Get structure information.  
  - In KLIFS OpenAPI called ``Structures``: https://dev.klifs.net/swagger_v2/#/Structures  

- ``bioactivities``  

  - Get bioactivity information.  
  - In KLIFS OpenAPI part of ``Ligands``: https://dev.klifs.net/swagger_v2/#/Ligands  

- ``interactions``  

  - Get interaction information.  
  - In KLIFS OpenAPI called ``Interactions``: https://dev.klifs.net/swagger_v2/#/Interactions  

- ``pocket``  

  - Get interaction information.  
  - In KLIFS OpenAPI part of ``Interactions``: https://dev.klifs.net/swagger_v2/#/Interactions/get_interactions_match_residues 

- ``coordinates``  

  - Get structural data (structure coordinates).
  - In KLIFS OpenAPI part of ``Structures``: https://dev.klifs.net/swagger_v2/#/Structures 

- ``conformations``

  - Get information on structure conformations.
  - In KLIFS OpenAPI part of ``Structures``: https://dev.klifs.net/swagger_v2/#/Structures/get_structure_conformation

- ``modified_residues``

  - Get information on residue modifications in structures.
  - In KLIFS OpenAPI part of ``Structures``: https://dev.klifs.net/swagger_v2/#/Structures/get_structure_modified_residues


Installing OpenCADD-KLIFS only
==============================

In case you would like to install the dependencies for the OpenCADD-KLIFS module only, please follow these instructions.

.. note::

    We are assuming you have a working ``mamba`` installation in your computer. 
    If this is not the case, please refer to their `official documentation <https://mamba.readthedocs.io/en/latest/installation.html#mamba>`_. 

    If you installed ``mamba`` into an existing ``conda`` installation, also make sure that the ``conda-forge`` channel is configured by running ``conda config --add channels conda-forge``.


1. Create a new conda environment called ``opencadd-klifs`` only with the dependencies needed for the OpenCADD-KLIFS module::

    mamba create -n opencadd-klifs bravado pandas tqdm rdkit biopandas

2. Activate the new conda environment::

    conda activate opencadd-klifs

3. Install the ``opencadd`` package without any dependencies (all OpenCADD-KLIFS-relevant dependencies have been installed in step 1 already)::

    mamba install opencadd --no-deps

   If you are planning on working with Jupyter notebooks, install JupyterLab and IPyWidgets::

    mamba install jupyterlab ipywidgets
For developers
==============


How the package will be structured
----------------------------------

Contributing to opencadd
-----------------------------------

First of all, thank you for considering making this project better.

You can help to make opencadd better by raising an issue to report a bug or suggest a
feature which is not implemented yet.
You can also create a new feature and add it via a Pull Request.

How Pull Requests work
----------------------

1. Fork the ``opencadd`` repository into your profile. Now you have your own copy of the repository.
2. Clone the repository. Run ``clone https://github.com/<your_username>/opencadd`` in your terminal.
3. Create a new branch to work on: ``git checkout -b meaningful-branch-name``.
4. Make your changes.
5. Add, commit and push your change to your forked repository.
6. Click on ``Create pull request`` Button.
7. Follow the template given there. Be aware of the requirements_.

The ``opencadd`` team will review the pull request and if there are no flaws it will be merged
to the ``master`` branch. If there are some flaws, the team will request you to fix them.

.. _requirements:

**************************************
What should each Pull Request contain?
**************************************

* Documentation (should be in ``rst`` format)
* Tests_ with Pytest
* Examples for the use of the new functions
* Little benchmark
* If 3rd party dependencies have been added, rationale behind that decision
* It should be formatted by black with ``-l 99``
* Short summary of the changes that were made
* If you implemented a new method or bigger feature, please porvide a tutorial for this method. For this you can follow this template.

Some people have trouble with NGLview. If that is the case for you, follow this `troubleshooting guide
<https://github.com/SBRG/ssbio/wiki/Troubleshooting#tips-for-nglview>`_.

LINK IS MISSING, WILL BE ADDED AS SOON AS TEMPLATE IS IN MASTER.

Formatting with black
---------------------

**1. Option:**

* Apply "black -l 99" on your code before committing::

        $> black -l 99 opencadd/


**2. Option:**

* Configuring your IDE to that effect
* Example in VS Code:

    * go to the settings
    * search for "python formatting"
    * set "Python › Formatting: Provider" to black
    * add "-l 99" to "Python › Formatting: Black Args"
    * activate "Editor: Format On Save"


.. _Tests:

How to add a new test
---------------------

- write a unit test for the new (or changed) function with `pytest
  <https://docs.pytest.org/en/latest/>`_.
- add new dependencies to ``test_env.yaml``


Steps made by the Github Actions Workflow
-----------------------------------------

The actions are executed automatically for every Pull Request submitted,
and for every commit pushed to ``master``. Steps run in Ubuntu and MacOS are:

* Report additional information about the test-build.
* Fixing conda in MacOS (to get the project).
* Creating the environment and getting all necessary dependencies.
* Installing the package in this environment.
* Running the tests.

The formating check is done in ubuntu.

* Checkout the code.
* Installing the linter (pylint) and the formatter (black).
* Running pylint (using  configuration at ``.pylintrc``).
* Running ``black -l 99 --check`` (check mode).
{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block attributes %}
   {% if attributes %}
   .. HACK -- the point here is that we don't want this to appear in the output, but the autosummary should still generate the pages.
      .. autosummary::
         :toctree:
      {% for item in all_attributes %}
         {%- if not item.startswith('_') %}
         {{ name }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% if methods %}
   .. HACK -- the point here is that we don't want this to appear in the output, but the autosummary should still generate the pages.
      .. autosummary::
         :toctree:
      {% for item in all_methods %}
         {%- if not item.startswith('_') or item in ['__call__'] %}
         {{ name }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
   {% endif %}
   {% endblock %}
:orphan:

{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}{{ fullname | escape | underline }}

.. rubric:: Description

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}

{% if classes %}
.. rubric:: Classes

.. autosummary::
    :toctree: .
    {% for class in classes %}
    {{ class }}
    {% endfor %}

{% endif %}

{% if functions %}
.. rubric:: Functions

.. autosummary::
    :toctree: .
    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% endif %}