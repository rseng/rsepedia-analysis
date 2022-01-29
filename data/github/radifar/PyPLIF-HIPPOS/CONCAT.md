# PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results of AutoDock Vina and PLANTS

[![GitHub Actions Build Status](https://github.com/radifar/PyPLIF-HIPPOS/workflows/CI/badge.svg)](https://github.com/radifar/PyPLIF-HIPPOS/actions?query=workflow%3ACI)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/radifar/PyPLIF-HIPPOS.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/radifar/PyPLIF-HIPPOS/context:python)
[![codecov](https://codecov.io/gh/radifar/pyplif-hippos/branch/main/graph/badge.svg)](https://codecov.io/gh/radifar/pyplif-hippos/branch/main)   
[![Anaconda-Server Badge](https://img.shields.io/badge/Install%20with-conda-green.svg?style=flat)](https://anaconda.org/conda-forge/pyplif-hippos)
[![Documentation Status](https://readthedocs.org/projects/pyplif-hippos/badge/?version=latest&style=flat)](https://pyplif-hippos.readthedocs.io/en/latest/)
[![DOI:10.1021/acs.jcim.0c00305](https://zenodo.org/badge/DOI/10.1021/acs.jcim.0c00305.svg)](https://doi.org/10.1021/acs.jcim.0c00305)   
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/conda-forge/pyplif-hippos?color=green)](https://anaconda.org/conda-forge/pyplif-hippos)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fradifar%2Fpyplif-hippos&title=visitor%20today%2Ftotal)](https://hits.seeyoufarm.com)

<p align="center">
  <img alt="Icons made by Freepik from Flaticon is licensed by CC 3.0 BY" src="docs/source/hippopotamus_small.png">
</p>

<p align="center">Icons made by <a href="https://www.freepik.com/">Freepik</a> from <a href="http://www.flaticon.com">Flaticon</a> is licensed by CC 3.0 BY</p>

Welcome to PyPLIF-HIPPOS's project page. PyPLIF-HIPPOS is an upgraded version of [PyPLIF](https://github.com/radifar/pyplif/) (**Python-based Protein-Ligand Interaction Fingerprinting**), a tool for molecular docking post-analysis. It will translate the 3D coordinates of both ligand(s) (generated from docking simulation) and protein into a series of *interaction bitstring* (also known as *Interaction Fingerprint*) (see image below). **HIPPOS** (/ˌhipoʊz/) is a recursive acronym of **HIPPOS Is PyPLIF On Steroids**. From this point forward, PyPLIF-HIPPOS is simplified to HIPPOS.

Compared to PyPLIF, HIPPOS is not only faster and able to generate more customized interaction bitstring, but also supports both [PLANTS](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/pharmazie-und-biochemie/pharmazie/pharmazeutische-chemie/pd-dr-t-exner/research/plants/) & [Vina](http://vina.scripps.edu/)! More over, unlike its predecessor it is (far) more well-documented.

<p align="center">
  <img alt="Table of Content Abstract Graphic JCIM" src="docs/source/toc-abstract-graphics_small.png">
</p>

<p align="center">Reprinted with permission from https://doi.org/10.1021/acs.jcim.0c00305. Copyright 2020 American Chemical Society.</p>

<p align="center">
  <img alt="PyPLIF output from PyPLIF publication" src="docs/source/pyplif.png">
</p>

<p align="center">Illustration by Radifar et al (2013) from <a href="http://www.bioinformation.net/009/97320630009325.htm">Bioinformation.net</a> is licensed by <a href="http://creativecommons.org/licenses/by/4.0">CC 4.0 BY</a>

## Quick Installation

The easiest way to install HIPPOS is using [Anaconda or Miniconda](https://docs.anaconda.com/anaconda/install/).
If you have Anaconda or Miniconda ready in your machine, you can start with
creating new environment (recommended):

`conda create -n hippos python=3.6`

Then activate the environment and install HIPPOS:

`conda activate hippos`  
`conda install -c conda-forge pyplif-hippos`

next you can try run HIPPOS and HIPPOS-GENREF with the following command:

`hippos`  
`hippos-genref` 

## How to Use HIPPOS

So I already installed HIPPOS, now what? Well you could start with how to generate
the [reference bitstring](https://pyplif-hippos.readthedocs.io/en/latest/getting-started-genref.html)
and Getting Started tutorial for [AutoDock Vina](https://pyplif-hippos.readthedocs.io/en/latest/getting-started-vina.html)
or [PLANTS](https://pyplif-hippos.readthedocs.io/en/latest/getting-started-plants.html).

## Ideas for Improvement? Found Bug(s)?

If you have any idea for improvement or found bug to report feel free to write them [here](https://github.com/radifar/PyPLIF-HIPPOS/issues).

## Citing HIPPOS

If you are using HIPPOS please cite this paper:

Istyastono, E., Radifar, M., Yuniarti, N., Prasasty, V. and Mungkasi, S., 2020. 
PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results 
of AutoDock Vina and PLANTS. Journal of Chemical Information and Modeling, 60(8), pp.3697-3702.
https://doi.org/10.1021/acs.jcim.0c00305

## Acknowledgment

This project has received funding from the [National Agency for Research and Innovation](https://international.ristekdikti.go.id/) (Indonesia)
under grant agreement No. 807.7/LL5/PG/2020. This project has been restructured based on the
[MOLSSI Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.3, and benefited greatly from [MOLSSI Python Package Development Best Practices](https://molssi.org/2020/04/20/may-webinar-series-python-package-development/)
workshop.


-----

&copy; Copyright 2021, Muhammad Radifar & Enade Perdana Istyastono
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
reported by contacting the project team at 'muhammad.radifar@picomps.org'. The project team will
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
e.g. `pyplif_hippos-0.1.2`, otherwise it will be appended with `+X` where `X` is the number of commits 
ahead from the last tag, and then `-YYYYYY` where the `Y`'s are replaced with the `git` commit hash.
# Sample Package Data

This directory contain the data required for testing. As integrating the content
into test function became impractical (e.g. inserting protein structure in Python
script). Also separating the test data with the script make the test code cleaner.
