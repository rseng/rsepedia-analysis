# <img align="center" src="https://github.com/carlosborca/CrystaLattE/blob/master/media/logo/Logo.png" height=260>

Automated calculation of crystal lattice energies with the many-body expansion.

| Category | Badges |
|-------------|-------------|
| **Status** | [![Travis Build Status](https://travis-ci.com/carlosborca/CrystaLattE.svg?branch=master)](https://travis-ci.org/carlosborca/CrystaLattE) [![codecov](https://codecov.io/gh/carlosborca/CrystaLattE/branch/master/graph/badge.svg)](https://codecov.io/gh/carlosborca/CrystaLattE/branch/master) [![Total alerts](https://img.shields.io/lgtm/alerts/g/carlosborca/CrystaLattE.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/carlosborca/CrystaLattE/alerts/) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/carlosborca/CrystaLattE.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/carlosborca/CrystaLattE/context:python) |
| **Foundation** | [![License](https://img.shields.io/github/license/carlosborca/CrystaLattE.svg)](https://opensource.org/licenses/LGPL-3.0) [![GitHub Top Languages](https://img.shields.io/github/languages/top/carlosborca/CrystaLattE)](https://github.com/carlosborca/CrystaLattE/) |
| **GitHub Info** | [![GitHub Code Size](https://img.shields.io/github/languages/code-size/carlosborca/CrystaLattE)](https://github.com/carlosborca/CrystaLattE/) [![GitHub Commits per Month](https://img.shields.io/github/commit-activity/m/carlosborca/CrystaLattE)](https://github.com/carlosborca/CrystaLattE/) [![GitHub Last Commit](https://img.shields.io/github/last-commit/carlosborca/CrystaLattE)](https://github.com/carlosborca/CrystaLattE/) |
| **Citation** | [![doi](https://img.shields.io/badge/DOI-10.1063%2F1.5120520-blue)](http://dx.doi.org/10.1063/1.5120520) |

## Overview

CrystaLattE is a software that automates the computation of crystal lattice energies using the many-body cluster expansion. The required computations on dimers, trimers, etc., within the crystal are independent of each other, leading to a naturally parallel approach. The algorithm exploits the long-range three-dimensional periodic order of crystals to automatically detect and avoid redundant or unnecessary computations.

## General Information

CrystaLattE has an interface with the quantum chemistry package PSI4. To run, the code requires a crystallographic information file containing structural information of the crystal and an input file specifying execution details. Work continues in the creation of a CrystaLattE `pip` package. So, for the moment, the instructions to download and install CrystaLattE and to create a _conda environment_ that includes PSI4 are presented below. 

### Installation

Minimal set of commands to install CrystaLattE on Linux, MacOS, or Windows (with the Windows Subsystem for Linux). Last tested on 4 October 2019:

#### 1. Install Miniconda:

If you have an installation of _Conda_ in your system, please skip to step 2. Otherwise, _Miniconda_ is required and the installer is available from the the Anaconda website. To download the installer from the terminal (in Linux or the Windows Subsystem for Linux): 

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

(_Note_) If you use MacOS, replace `Linux` by `MacOSX` on the previous command.

Run the installer following the on-screen instructions:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

After the installation is complete, close the terminal and start a new shell.

(_Optional_) Disable automatic activation of the _base_ conda environment:

```
conda config --set auto_activate_base false
```

#### 2. Create a _Conda Environment_ for CrystaLattE

CrystaLattE requires PSI4 and PyCIFRW. Conda offers the possibility of creating an _environment_ that contains all the dependencies required by CrystaLattE. To download and install PSI4 and other related software tools in a new _cle_ environment execute the command below and follow the on-screen instructions:

```
conda create -n cle python=3.7 psi4 pycifrw -c psi4/label/dev -c psi4 -c conda-forge
```

#### 3. Activate the _cle_ environment

To use the recently created _cle_ environment, activate it:

```
conda activate cle
```

#### 4. Clone CrystaLattE from its GitHub repository:

In you file system navigate to the location where you would like to place the root directory of CrystaLattE and clone it from its corresponding GitHub repository:

```
git clone https://github.com/carlosborca/CrystaLattE.git
```

#### 5. Test CrystaLattE

Go to the root directory of the CrystaLatte repository in the file system and run the test suite:

```
pytest
```

(_Note_) It is known that the `v2rdm_casscf` plugin of PSI4 may cause execution errors when trying to import PSI4 from CrystaLattE. If such error is encountered during testing, remove it.

```
conda remove v2rdm_casscf
```

### How to run CrystaLattE

CrystaLattE requires a crystallographic information file (.cif) and an options input file (.cle). CIF files can be obtained from multiple sources. For example, from the Cambridge Crystallographic Data Centre (CCDC) website. The options file must be generated by the user. A template is presented below:

#### Preparing the options input

Example of an options input file for CrystaLattE.

```
# This is a typical CrystaLattE input template: input.cle

# Blank lines and lines starting with the hash character are omitted.
# Padding blank spaces are also ignored.

# The format is:
# Keywords  equal  Value

cif_input       =  ../MyCrystals/OneCrystal.cif
cif_output      =  ../MyCrystals/SuperCell.xyz
cif_a           =  5
cif_b           =  5
cif_c           =  5
bfs_thresh      =  1.2
uniq_filter     =  ChSEV
nmers_up_to     =  3
r_cut_com       =  12.0
r_cut_monomer   =  15.0
r_cut_dimer     =  9.0
r_cut_trimer    =  12.0
r_cut_tetramer  =  6.0
r_cut_pentamer  =  4.0
cle_run_type    =  psi4api + quiet
psi4_method     =  MP2/aug-cc-pV[TQ]Z + D:FNO-CCSD(T)/aug-cc-pVDZ
psi4_bsse       =  cp
psi4_memory     =  8 GB
verbose         =  2
```

#### Keywords

Description of keywords and their acceptable values.

| Keyword          | Acceptable Values                              | Default      | Example      | Notes                |
|------------------|------------------------------------------------|--------------|--------------|----------------------|
| `cif_input`      | \<RelativePath\>\/\<FileName\>.cif             | *No default* | Benzene.cif  | Must be a .cif       |
| `cif_output`     | \<RelativePath\>\/\<FileName\>.xyz             | sc.xyz       | Benzene.xyz  | Must be a .xyz       |
| `cif_a`          | Odd positive integer                           | 5            | 9            |                      |
| `cif_b`          | Odd positive integer                           | 5            | 7            |                      |
| `cif_c`          | Odd positive integer                           | 5            | 11           |                      |
| `bfs_thresh`     | Positive float                                 | 1.2          | 1.3          | vdW radii multiplier |
| `uniq_filter`    | ChSEV, Dreamaligner                            | ChSEV        | Dreamaligner |                      |
| `nmers_up_to`    | 2, 3, 4, 5                                     | 2            | 3            | Dimers, Trimers...   |
| `r_cut_com`      | Positive float                                 | 10.0         | 12.0         | Angstroms            |
| `r_cut_monomer`  | Positive float                                 | 12.0         | 15.0         | Angstroms            |
| `r_cut_dimer`    | Positive float                                 | 10.0         | 8.0          | Angstroms            |
| `r_cut_trimer`   | Positive float                                 | 8.0          | 10.0         | Angstroms            |
| `r_cut_tetramer` | Positive float                                 | 6.0          | 5.0          | Angstroms            |
| `r_cut_pentamer` | Positive float                                 | 4.0          | 3.0          | Angstroms            |
| `cle_run_type`   | psi4api, psithon, makefp, test, quiet, timings | psi4api      | test + quiet | Separate with +      |
| `psi4_method`    | String                                         | HF/STO-3G    | HF-3c        | See PSI4 manual      |
| `psi4_bsse`      | vmfc, cp, nocp                                 | cp           | nocp         |                      |
| `psi4_memory`    | String                                         | 500 MB       | 2 GB         |                      |
| `verbose`        | 0, 1, 2                                        | 1            | 2            |                      |

#### Running CrystaLattE

Finally, to execute the code:

```
./crystalatte.py input.cle
```

#### Known Issues

Precision issues may arise when computing a large number of structures at default energy and density convergence criteria if using Psi4 1.3.2 and newer versions. Although the default energy and density convergence criteria for Psi4 calculations has been increased in CrystaLattE, the user should be careful.

#### Copyright

Copyright (c) 2020, Carlos H. Borca


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
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
reported by contacting the project team at 'carlosborca@gmail.com'. The project team will
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
e.g. `crystalatte-0.1.2`, otherwise it will be appended with `+X` where `X` is the number of commits 
ahead from the last tag, and then `-YYYYYY` where the `Y`'s are replaced with the `git` commit hash.
# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into this crystalatte. 

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
  navigate to your fork of crystalatte on GitHub and open a [pull
  request](https://help.github.com/articles/using-pull-requests/) (PR). Note that
  after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.  Each commit added to the PR will be validated for
  mergability, compilation and test suite compliance; the results of these tests
  will be visible on the PR page.
* If you're providing a new feature, you must add test cases and documentation.
* When the code is ready to go, make sure you run the test suite using pytest.
* When you're ready to be considered for merging, check the "Ready to go"
  box on the PR page to let the crystalatte devs know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration returns checkmarks,
  and multiple core developers give "Approved" reviews.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
## Description
Provide a brief description of the PR's purpose here.

## Todos
Notable points that this PR has either accomplished or will accomplish.
  - [ ] TODO 1

## Questions
- [ ] Question1

## Status
- [ ] Ready to go# Sample Package Data

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
# Compiling CrystaLattE's Documentation

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
be inside a directory called `html/` depending on what version of Sphinx is installed).# Templates Doc Directory

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
# Static Doc Directory

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
