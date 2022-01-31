**Dev**

**1.6.1**

- Add option --ignore-CH3s (-igch3) which ignore CH pairs when they belongs to a CH3 group.
- Update doc and notebook04 with option -igch3
- Fix beta H reconstruction for Berger and GROMOSCPK POPS
- Update .readthedocs.yml 

**1.6.0**

- Avoid output trajectory rewind when writing box dimensions
- Switch to MDAnalysis 2.0
- Add support of Python 3.9
- Improve docstrings

**1.5.0**

- Write box dimensions in the requested trajectory output
- Fix write duplicate 1st frame when a trajectory output is requested
- Avoid using universe.trajectory.time on a single pdb
- Limit Python version >= 3.6 <=3.8 (for MDAnalysis compatibility)
- Add support for: Berger DOPC/DPPC/POPS, GROMOS-CKP POPC/POPS, GROMOS-53A6L DPPC, CHARMM36UA
- Force Python 3.8 for doc building

**1.4.0**

- Add -v / --version option
- Reorganize doc
- Add Notebook04 (launch buildH as a module)
- Support Berger cholesterol
- Add Notebook05 (mixture POPC / cholesterol)
- Create buildH logo and add it to doc
- Add paper for JOSS
- Add community guidelines

**1.3.1**

- Fix setup.cfg to include json files in python package archive

**1.3.0**

- Complete documentation
- Accelerate functions within geometry.py with Numba
- Implement the use of buildH as a module
- Simplify calculation of CH on an sp3 carbon
- Use MyST parser for documentation (handles latex equations)
- Clarify some error messages
- Fix residue number exceeding 9999
- Add POPE def and json files
- Add Notebook01 (basic buildH analysis on a Berger traj)
- Add Notebook02 (+trajectory output)
- Add Notebook03 (analysis on a mixture POPC/POPE)
- Move CHARMM36 POPC validation to Zenodo

**1.2.0**

- Build docs
- Rename '-x/--xtc' flag to -t/--traj' one to be more generic
- Replace mandatory topology argument to '-c/--coord' flag
- Improve performance of control functions.
- Move misc functions to a module utils.py
- Improve Exception handling & add proper exits
- Improve PEP8 & PEP257 compliance
- Improve test coverage
- Fix bug when a trajectory was written when only a pdb was provided.
- Add sanity checks for the various input files
- Use json files instead of python module to read lipid topologies.
- Optimize package for better performance

**1.1.0**

- Create Python package structure
- Create conda environment
- Fix tests
- Separate entry point
- Update README for dev version installation
- Handle version with bump2version
# buildH

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03521/status.svg)](https://doi.org/10.21105/joss.03521)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4676217.svg)](https://doi.org/10.5281/zenodo.4676217)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/patrickfuchs/buildH/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/patrickfuchs/buildH/)
[![License: BSD](https://img.shields.io/badge/License-BSD-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/patrickfuchs/buildH/master?urlpath=lab)
[![Code CI Status](https://github.com/patrickfuchs/buildH/workflows/GitHub%20CI%20code/badge.svg)](https://github.com/patrickfuchs/buildH/actions?query=workflow%3A%22GitHub+CI+code%22)
[![Doc CI Status](https://github.com/patrickfuchs/buildH/workflows/GitHub%20CI%20doc/badge.svg)](https://github.com/patrickfuchs/buildH/actions?query=workflow%3A%22GitHub+CI+doc%22)
[![Documentation Status](https://readthedocs.org/projects/buildh/badge/?version=latest)](https://buildh.readthedocs.io/en/latest/?badge=latest)
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![buildH version on PyPI](https://badge.fury.io/py/buildh.svg)](https://pypi.python.org/pypi/buildh)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/buildh/badges/version.svg)](https://anaconda.org/bioconda/buildh)

![buildH_logo](docs/img/buildH_logo.png)

> Build hydrogen atoms from united-atom molecular dynamics of lipids and calculate the order parameters

## Features

**buildH** can :
  - Reconstruct hydrogens from a **united-atom** structure file (pdb, gro) or trajectory (e.g. xtc).
  - Calculate the order parameters based on the reconstructed hydrogens.
  - Write new structure trajectory files with the reconstructed hydrogens.

**buildH** works in two modes :
  1. A slow mode when an output trajectory (in xtc format) is requested by
     the user. In this case, the whole trajectory including newly built
     hydrogens are written to this trajectory.
  2. A fast mode without any output trajectory.

In both modes, the order parameters are calculated. All calculations are accelerated with [Numba](https://numba.pydata.org/). As a CPU cost indication, running **buildH** on a trajectory of 2500 frames with 128 POPC (without trajectory output) takes ~ 7' on a single core Xeon @ 3.60GHz.

## Requirements

Python >= 3.6 is mandatory for running buildH.

**buildH** is written in Python 3 and needs the modules numpy, pandas, MDAnalysis and Numba.

## Installation

### Pip

A simple installation with pip will do the trick:

```
python3 -m pip install buildh
```

All dependencies (modules) will be installed automatically by pip.


### Conda

**buildH** is also available through the [Bioconda](https://anaconda.org/bioconda/buildh) channel:

```
conda install buildh -c bioconda -c conda-forge
```

More details on installation [here](https://buildh.readthedocs.io/en/latest/installation.html).

For installing a development version, see [here](devtools/install_dev.md).

## Running buildH

Once installed, a simple invocation of the `buildH` command will run the program (`$` represents the Unix prompt):

```
$ buildH
usage: buildH [-h] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]] -d DEFOP
              [-opx OPDBXTC] [-o OUT] [-b BEGIN] [-e END] [-pi PICKLE] [-igch3]
buildH: error: the following arguments are required: -c/--coord, -l/--lipid, -d/--defop
```

The minimal command for running **buildH** can resemble this:

```
$ buildH -c start_128popc.pdb -t popc0-25ns_dt1000.xtc -l Berger_POPC -d Berger_POPC.def
```

The different arguments mean the following: `-c start_128popc.pdb` is a pdb file with 128 POPC molecules, `-t popc0-25ns_dt1000.xtc` is a trajectory with 25 frames, `-l Berger_POPC` indicates the united-atom force field and the type of lipid to be analyzed, `-d Berger_POPC.def` indicates what C-H are considered for H building and order parameter calculation (the structure and trajectory files can be found [here](https://github.com/patrickfuchs/buildH/tree/master/docs/Berger_POPC_test_case)). The def file can be found [here](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def). The final order parameters averaged over the trajectory will be written to the default output name `OP_buildH.out`

Other detailed examples and Jupyter Notebooks can be found in the documentation at [Read the Docs](https://buildh.readthedocs.io/en/latest/index.html).

**Important**: sometimes, when performing MD, some molecules are split over periodic boundary conditions (PBC). **buildH** takes as input whole structures (pdb, gro, xtc, etc.). If broken molecules are supplied, it will most likely generate nonsense results. So it is up to the user to take care of making molecules whole before running **buildH** (e.g. by using a tool like [trjconv](https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html) in GROMACS with flag `-pbc mol`).

Invoking **buildH** with the `-h` flag will display some help to the screen and tell which lipids are supported.

```
$ buildH -h
usage: buildH [-h] [-v] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]] -d
              DEFOP [-opx OPDBXTC] [-o OUT] [-b BEGIN] [-e END] [-igch3]
[...]
The list of supported lipids (-l option) are: Berger_CHOL, Berger_DOPC, Berger_DPPC, Berger_POP, Berger_POPC, Berger_PLA, Berger_POPE, Berger_POPS, CHARMM36UA_DPPC, CHARMM36UA_DPUC, CHARMM36_POPC, GROMOS53A6L_DPPC, GROMOSCKP_POPC, GROMOSCKP_POPS. More documentation can be found at https://buildh.readthedocs.io.
```

## Documentation

The full documentation is available at [Read the Docs](https://buildh.readthedocs.io/en/latest/index.html).

## Contributors

- Hubert Santuz
- Amélie Bâcle
- Pierre Poulain
- Patrick Fuchs

## License

**buildH** is licensed under the [BSD License](LICENSE.txt).


## Contributing

If you want to report a bug, request a feature, or propose an improvement use the [GitHub issue system](https://github.com/patrickfuchs/buildH/issues/).

Please, see also the [CONTRIBUTING](CONTRIBUTING.md) file.

Note that this project is released with a [Contributor Code of
Conduct](http://contributor-covenant.org/). By participating in this project you
agree to abide by its terms. See the [CODE_OF_CONDUCT](CODE_OF_CONDUCT.md) file.

## Citing **buildH**

If you use **buildH** for your research, please cite :

```
Santuz et al., (2021). buildH: Build hydrogen atoms from united-atom molecular dynamics of lipids and calculate the order parameters. Journal of Open Source Software, 6(65), 3521, https://doi.org/10.21105/joss.03521
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
reported by contacting the project team at patrick.fuchs@u-paris.fr . All
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
# Contributing Guidelines

First, **thank you** for contributing!

Before starting, you should read, agree to, and follow these three points:

  - [How to contribute?](#how-to-contribute)
    - [Report bugs](#report-bugs)
    - [Submit feedback](#submit-feedback)
  - [Pull Request Guidelines](#pull-request-guidelines)
  - [Code of Conduct](CODE_OF_CONDUCT.md)

This page is largely inspired from a [blog post of William Durand](https://williamdurand.fr/2013/11/20/on-creating-pull-requests/).

---

## How to contribute?

### Report bugs

Report bugs at: https://github.com/patrickfuchs/buildH/issues/new.

When reporting a bug, please include:

* Any details about your local setup which might be helpful in troubleshooting
* Detailed steps to reproduce the bug. Where possible, please write a test case

If you are not able to do that, that's fine! Open an issue anyway and let us
know as much information as you can. We will get back to you to determine the
problem, and (hopefully) fix it.


### Submit feedback

The best way to send feedback is to [create a new
issue](https://github.com/patrickfuchs/buildH/issues/new) on GitHub.

If you are proposing a feature:

* Explain how you see it working. Try to be as detailed as you can.
* Try to keep the scope as narrow as possible. This will help make it easier to
  implement.
* Feel free to include any code you might already have, even if it is
  just a rough idea. This is a volunteer-driven project, and contributions are
  welcome :)


## Pull Request Guidelines

Here are a few rules to follow in order to make code reviews and discussions go
more smoothly before maintainers accept and merge your work:

* You MUST run the test suite.
* You MUST write (or update) unit tests.
* You SHOULD write documentation.

Please, write [commit messages that make
sense](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html),
and [rebase your branch](http://git-scm.com/book/en/Git-Branching-Rebasing)
before submitting your Pull Request.

You may be asked to [squash your
commits](http://gitready.com/advanced/2009/02/10/squashing-commits-with-rebase.html)
too. This is to "clean" your Pull Request before merging it (we try to avoid
commits such as `fix tests`, `fix 2`, `fix 3`, etc.).

Also, while creating your Pull Request on GitHub, you MUST write a description
which gives the context and/or explains why you are creating it.

For further information about creating a Pull Request, please read [this blog
post](http://williamdurand.fr/2013/11/20/on-creating-pull-requests/).

Thank you!
This directory contains def files which tell **buildH** which C-H to analyze and give them a generic name for a more readable output. Initially, the use of def files comes from the [NMRlipids](https://nmrlipids.blogspot.com/) project, some def files for all-atom force fields can be found [here](https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs).

The global format is the following:

```
name            residue   carbon atom    hydrogen atom 
 
gamma1_1         POPC        C13             H13A 
gamma1_2         POPC        C13             H13B 
gamma1_3         POPC        C13             H13C 
. 
. 
.
```

Each column has to be separated by any combination of whitespaces (at least one).

Last, in the [NMRlipids](https://nmrlipids.blogspot.com/) project, the def files are generally named like this: `order_parameter_definitions_MODEL_CHARMM36_POPC.def`. In **buildH** we use merely `CHARMM36_POPC.def`, or more generally `Forcefield_Lipid.def` where `Lipid` is a residue name as found in the pdb or gro file.


More on def files can be found in the [documentation](https://buildh.readthedocs.io/en/latest/def_format.html).
---
title: 'buildH: Build hydrogen atoms from united-atom molecular dynamics of lipids and calculate the order parameters'
tags:
- python
- molecular-dynamics-simulation
- order-parameters
- lipids
- united-atom
authors:
- name: Hubert Santuz
  orcid: 0000-0001-6149-9480
  affiliation: "1, 2"
- name: Amélie Bacle
  orcid: 0000-0002-3317-9110
  affiliation: 3
- name: Pierre Poulain
  orcid: 0000-0003-4177-3619
  affiliation: 4
- name: Patrick F.J. Fuchs^[corresponding author]
  orcid: 0000-0001-7117-994X
  affiliation: "5, 6"
affiliations:
- name: CNRS, Université de Paris, UPR 9080, Laboratoire de Biochimie Théorique, 13 Rue Pierre et Marie Curie, F-75005 Paris, France
  index: 1
- name: Institut de Biologie Physico-Chimique–Fondation Edmond de Rothschild, PSL Research University, Paris, France
  index: 2
- name: Laboratoire Coopératif "Lipotoxicity and Channelopathies - ConicMeds", Université de Poitiers, F-86000 Poitiers, France
  index: 3
- name: Université de Paris, CNRS, Institut Jacques Monod, F-75006, Paris, France
  index: 4
- name: Sorbonne Université, Ecole Normale Supérieure, PSL Research University, CNRS, Laboratoire des Biomolécules (LBM), F-75005 Paris, France
  index: 5
- name: Université de Paris, UFR Sciences du Vivant, F-75013 Paris, France
  index: 6
date: 27 May 2021
bibliography: paper.bib
---

# Background

Molecular dynamics (MD) simulations of lipids are widely used to understand the complex structure and dynamics of biological or model membranes [@Tieleman1997; @Feller2000; @Lyubartsev2011]. They are very complementary to biophysical experiments and have thus become important to get insights at the microscopic scale. Many of them are performed using all-atom (AA) or united-atom (UA) representations. In AA force fields (such as CHARMM36 [@Klauda2010]), all the atoms are considered whereas in UA force fields (such as Berger [@Berger1997]) the aliphatic hydrogen atoms (Hs) are merged with their parent carbon into a larger particle representing a CH, CH2 or CH3 (e.g. a methyl group is represented by a single CH3 particle). In simulations of phospholipids, the use of UA representations allows one to reduce the number of particles that are simulated to almost a third of the true number of atoms. This is because lipid molecules contain many aliphatic Hs. This simplification thus reduces the computational cost without losing important chemical details.

MD simulations of lipids are usually validated against experimental data [@Klauda2010] or used to help interpret experiments [@Feller2007]. One type of experiment which is often used for that is $^2H$ NMR. In this type of experiment, aliphatic Hs are replaced by deuterons. $^2H$ NMR allows one to measure the order parameter of a given C-H bond (where the H is replaced by a deuteron):

$$S_{CH} = \frac{1}{2} \left \langle 3cos^2(\theta) -1 \right \rangle$$

where $\theta$ is the angle between the C-H bond and the magnetic field. The symbol $\langle ... \rangle$ means averaging over time and molecules.

This order parameter is useful because it is directly related to the flexibility of the given C-H bond. It describes the amount of possible orientations visited by the C-H bond, ranges between -0.5 and 1 (inclusive), and is unitless. Values close to 0 indicate high mobility, while values with higher magnitudes (either positive or negative) indicate that the C-H bond is less mobile. Although it varies theoretically between -0.5 and 1, its absolute value $\lvert S_{CH} \rvert$ is often reported because its sign is usually difficult to measure experimentally [@Ollila2014]. In MD simulations, $\theta$ is the angle between the C-H bond and the normal to the membrane (usually the $z$ axis). For AA simulations, $S_{CH}$ is trivial to calculate. However, it is more difficult for UA simulations since the aliphatic Hs are not present. There are two main strategies to compute $S_{CH}$ from UA simulations [@Piggot2017]: i) expressing $S_{CH}$ as a function of the coordinates of other atoms [@Douliez1995], ii) reconstructing hydrogen coordinates and calculating $S_{CH}$ as in AA simulations. The trend in the recent years has been towards strategy ii), such as in the NMRlipids project [@Botan2015; @Catte2016; @Antila2019; @Bacle2021]. NMRlipids is an open science project which uses MD simulations and experimental $S_{CH}$ with the goal of improving lipid force fields or conducting fundamental research on the structure and dynamics of lipid membranes.

# Statement of Need

Reconstructing Hs from the heavy atom coordinates can be done using standard geometric rules respecting stereochemistry. In the first NMRlipids project [@Botan2015], H reconstruction was performed using a program called `g_protonate` from the software GROMACS version 3.* [@Berendsen1995]. However, `g_protonate` has been removed from GROMACS version 4.* or higher. So, there was a need to find other solutions. Currently, there are many programs in the field of chemoinformatics that are able to build Hs such as OpenBabel [@OBoyle2011] or other proprietary softwares. It is also possible to use `pdb2gmx` from the GROMACS software [@Abraham2015] to build Hs, but it is not initially intended for that since its main purpose is to build a topology. Many of these solutions remain workarounds where one uses a sophisticated software for just doing one basic task. Using these open or proprietary software usually has a slow learning curve. Moreover, it is not always easy to use them in the Unix command line when they have complicated graphical interfaces, thus complicating their use in automatic pipelines of trajectory analyses.
Here, we propose the `buildH` software to meet this need. `buildH` is very light and usable in the Unix command line or as a Python module, making it a tool of choice to integrate in analysis pipelines or Jupyter notebooks. It can be easily extended or customized if needed. `buildH` is currently widely used in the NMRlipids project IVb dealing with the conformational plasticity of lipid headgroups in cellular membranes and protein-lipid complexes [@Bacle2021]. In addition, it is planned to use `buildH` in the next NMRlipids project dealing with a databank containing MD trajectories of lipids [@Kiirikki2021].

# Overview

`buildH` is a Python software that enables automatic analyses of order parameter calculations from UA trajectories of lipids. The software has the following features:

- It reads a single structure or a trajectory of lipids in a UA representation.
- It reconstructs the aliphatic Hs using standard geometric rules.
- From the reconstructed Hs, it calculates and outputs the order parameters on each requested C-H bond.
- Optionally, it outputs a structure (in pdb format) and a trajectory (in xtc format) with all reconstructed Hs.

Beyond order parameter calculations, the trajectory with Hs can be used for any further analyses (e.g. precise molecular volume calculation).

`buildH` has been natively developed for a use in the Unix command line. It possesses a minimum number of options, making it easy to use. It is also possible to utilize it as a Python module, which may be convenient in some cases.

To reconstruct H atoms, `buildH` uses standard geometric rules. These rules require so-called *helper* atoms. For example, the reconstruction of the two Hs of a CH2 on carbon $C_i$, requires two helpers which are $C_{i-1}$ and $C_{i+1}$, that is, the two neighbors of $C_i$ in the chain (note that helpers can also be other heavy atoms such as oxygen or nitrogen). The list of helpers used for the reconstruction of each H is written in a json file. Many json files are already present on the `buildH` repository representing the major lipids: phosphatidylcholine (PC), phosphatidylethanolamine (PE), phosphatidylglycerol (PG), phosphatidylserine (PS) for the polar heads and palmitoyl, myristoyl, oleoyl for the aliphatic chains, as well as cholesterol. Major UA force fields are also represented (Berger [@Berger1997], GROMOS-CKP [@Piggot2012], CHARMM-UA [@Lee2014]). In case a user wants to analyze a lipid which is not present in `buildH`, a step-by-step documentation guides the user in the process of creating and supplying his/her lipid description as a json file.

All structure and trajectory read / write operations are handled by the MDAnalysis module [@Michaud-Agrawal2011; @Gowers2016]. Mathematical vector operations are performed by Numpy [@Harris2020] and accelerated with Numba [@Lam2015], leading to very decent performances. For example, the reconstruction of all Hs and order parameter calculation on a trajectory of 2500 frames with 128 POPC molecules can be handled in approximately 7 minutes using a single core Xeon processor at 3.60 GHz, whereas it was almost 30 minutes long without Numba.

`buildH` has been implemented with good practices of software development in mind [@Jimenez2017; @Taschuk2017]:

- version controlled repository on GitHub [https://github.com/patrickfuchs/buildH](https://github.com/patrickfuchs/buildH),
- open-source license (BSD-3-Clause),
- continuous integration through tests,
- and documentation [https://buildh.readthedocs.io/](https://buildh.readthedocs.io/).

Some notebooks are provided in the GitHub repository to explain how `buildH` works and how to analyze the data produced. In case of trouble, any user can post an issue on GitHub.

`buildH` is available in the Python Package Index (PyPI) as well as in the Bioconda repository. The current version 1.6.0 of `buildH` is archived in the Zenodo repository ([https://zenodo.org/record/5356246](https://zenodo.org/record/5356246)) and in the Software Heritage archive ([swh:1:dir:4c63d5ca3497726a1e54ac152ce1667d7c004d2b](https://archive.softwareheritage.org/swh:1:dir:4c63d5ca3497726a1e54ac152ce1667d7c004d2b;origin=https://github.com/patrickfuchs/buildH/;visit=swh:1:snp:a63a8d07dbebeb442a06707be476817cec44ac72;anchor=swh:1:rev:9f05672515e1cdb0064eeb34f63844296193bc0d)).

# Acknowledgements

The authors thank the community of [NMRlipids](http://nmrlipids.blogspot.com/) for useful discussions, especially Samuli Ollila.

# References

In this file is explained how to install buildH for developers.

## Installation (development)

1. Install conda (either with Miniconda or Anaconda, we recommend Miniconda)

2. Clone this GitHub repository:
```
git clone https://github.com/patrickfuchs/buildH.git
cd buildH
```

3. Create conda environment:
```
conda env create -f binder/environment.yml
conda activate buildh
```

If needed, update your conda env with
```
conda env update -f binder/environment.yml
```

4. Install the dev version of buildH:
```
pip install -e .
```
# How to release


## Setup

Install required packages:
```
$ conda env create -f binder/environment.yml
```

Or if needed, update your conda environment:
```
$ conda env update -f binder/environment.yml
```

For Zenodo integration, see [Making Your Code Citable](https://guides.github.com/activities/citable-code/).

To publish a package in [PyPI](https://pypi.org/):

- Create an [account](https://pypi.org/account/register/).
- Create a new API token in the [account settings](https://pypi.org/manage/account/#api-tokens). Copy this token because you won't be able to see it again.
- Paste this token in the [GitHub secrets for your repo](https://docs.github.com/en/actions/reference/encrypted-secrets#creating-encrypted-secrets-for-a-repository) with the name `PYPI_API_TOKEN`.

## Tests

Before publishing any release, double-check all tests had run successfully:
```
$ make tests
```


## Update version number

We use `bump2version` to update and synchronize the version number across different files.

For patch update (x.y.z → x.y.**z+1**):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg patch
```

For minor update (x.y.z → x.**y+1**.0):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg minor
```

For major update (x.y.z → **x+1**.0.0):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg major
```

Remark:

1. For a dry run with `bump2version`, use option `-n`.
2. `bump2version` will fail if the git working directory is not clean, i.e. all changes are not commited.

Once version number is updated, push everything to GitHub:
```
$ git push origin
$ git push origin --tags
```


## Add new release on GitHub

On [GitHub release page](https://github.com/patrickfuchs/buildH/releases) :

- Click the *Draft a release* button.
- Select the latest version as *tag version*.
- Add release version as *Release title* (e.g.: v1.3.7).
- Copy recent changes from `CHANGELOG.md` and paste them in the *Describe this release* field.
- Hit the *Publish Release* button :rocket:


## Zenodo integration

After the creation of the new release on GitHub, check a new archive has been created on [Zenodo](https://doi.org/10.5281/zenodo.4676217).


## PyPI package

After the creation of the new release on GitHub, check a new package has been published on [PyPI](https://pypi.org/project/buildh/).


If you need to manually build and upload your package on PyPI, run the following commands:

```bash
$ make build
$ make upload-to-pypi
```

Enter your username and password upon request.


## Bioconda package

The `meta.yaml` file used to create the conda package is under the folder `bioconda`.

Here a link to the merged PR which add buildH to the bioconda channel : https://github.com/bioconda/bioconda-recipes/pull/28673.

Once a new release is made on GitHub, it should be automatically updated on the bioconda channel (https://bioconda.github.io/contributor/updating.html).
# Usage


## Simple examples

The examples below are based on a simple test case using Berger POPC. The files can be found on github in the directories [docs/Berger_POPC_test_case](https://github.com/patrickfuchs/buildH/tree/master/docs/Berger_POPC_test_case) and [def_files](https://github.com/patrickfuchs/buildH/tree/master/def_files). You will need 3 files :

- [`start_128popc.pdb`](https://github.com/patrickfuchs/buildH/blob/master/docs/Berger_POPC_test_case/start_128popc.pdb): contains 128 POPC.
- [`popc0-25ns_dt1000.xtc`](https://github.com/patrickfuchs/buildH/blob/master/docs/Berger_POPC_test_case/popc0-25ns_dt1000.xtc): contains a small trajectory of 25 frames.
- [`Berger_POPC.def`](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def): contains a list of C-H which tells **buildH** what hydrogens to reconstruct, what C-H to calculate the order parameters on.


Here are some examples on how to run **buildH** with these 3 files:

### Basic run on a single structure

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def
```

**buildH** can be used on a single structure (OK not very common for research, but useful for debugging ;-)). The pdb structure is passed with option `-c` (it also works with gro files), the def file with `-d`. The flag `-l` is mandatory, it tells **buildH** what force field and lipid to use: here it is `Berger_POPC`. The order parameters will be written to `OP_buildH.out` which is the default name.

### Same but with a chosen output name

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def \
-o my_OP_buildH.out
```

Here we add a `-o` flag which tells **buildH** to output the results in a file named `my_OP_buildH.out`.

### Run on a trajectory

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def \
-t popc0-25ns_dt1000.xtc
```

Here the flag `-t` indicates a trajectory. The final order parameters will be averaged over all lipids and all frames for each C-H present in the def file. More can be found on how the averaging is done [here](order_parameter.md). The default name `OP_buildH.out` will be used.

### Same with an output trajectory with reconstructed hydrogens

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def \
-t popc0-25ns_dt1000.xtc -opx popc0-25ns_dt1000_with_H
```

Here we added the flag `-opx` to request a pdb and an xtc file of the system with all the reconstructed hydrogens. Note that the flag takes a base name without extension since it will create a pdb and an xtc, here `popc0-25ns_dt1000_with_H.pdb` and `popc0-25ns_dt1000_with_H.xtc`. The use of this flag `-opx` requires the `.def` file to contain **all possible pairs of C-H** to reconstruct (since the trajectory with all Hs will be reconstructed). Importantly, the newly built hydrogens in the output pdb will be named according to the names written in the def file. See more about this [here](def_format.md). The order parameters will be written in `OP_buildH.out` (default name).

### Get a single pdb file with reconstructed hydrogens

If you do not provide a trajectory with the `-t` flag and you use the `opx` flag, **buildH** will only output a pdb file with hydrogens (no xtc will be produced):

```bash
buildH -c start_128popc.pdb -l Berger_POPC -d Berger_POPC.def \
-opx start_128popc_wH
```

In this case, the file `start_128popc_wH.pdb` with reconstructed hydrogens will be created as well as `OP_buildH.out` with the order parameters.

## Additional details

### Why do I need a def file?

A def file looks like this:

```
gamma1_1 POPC C1  H11
gamma1_2 POPC C1  H12
gamma1_3 POPC C1  H13
[...]
```

Each line corresponds to a given C-H. The 4 columns correspond to the generic name, residue name, carbon name and hydrogen name, respectively, for that C-H.

In **buildH**, the def file has three main purposes:

- Tell what are the C-H we want to consider for H reconstruction and order parameter calculation.
- Give a generic name to each C-H (which will appear in the output) and make the correspondance with the PDB names (e.g. `gamma1_1` stands for the C-H which have `C1` and `H11` atom names in the pdb file.
- If an output file with the newly built hydrogens is requested, their names will follow the 4th column of the def file. For example, the 3 hydrogens reconstructed on atom `C1` will be named `H11`, `H12` and `H13`.


In the following example dealing with a Berger POPC, the order parameters will be calculated on the polar head only (excluding the CH3s of choline):

```
beta1  POPC C5  H51
beta2  POPC C5  H52
alpha1 POPC C6  H61
alpha2 POPC C6  H62
g3_1   POPC C12 H121
g3_2   POPC C12 H122
g2_1   POPC C13 H131
g1_1   POPC C32 H321
g1_2   POPC C32 H322
```

Using the [Berger POPC trajectory](https://github.com/patrickfuchs/buildH/tree/master/docs/Berger_POPC_test_case) of 25 frames, the output `OP_buildH.out` will contain the order parameters of the C-H specified in the def file:

```
# OP_name            resname atom1 atom2  OP_mean OP_stddev OP_stem
#--------------------------------------------------------------------
beta1                POPC    C5    H51    0.04934  0.11999  0.01061
beta2                POPC    C5    H52    0.07162  0.12108  0.01070
alpha1               POPC    C6    H61    0.11839  0.15261  0.01349
alpha2               POPC    C6    H62    0.13903  0.19003  0.01680
g3_1                 POPC    C12   H121  -0.28674  0.09135  0.00807
g3_2                 POPC    C12   H122  -0.16195  0.14832  0.01311
g2_1                 POPC    C13   H131  -0.15159  0.14511  0.01283
g1_1                 POPC    C32   H321   0.21133  0.22491  0.01988
g1_2                 POPC    C32   H322   0.09638  0.16189  0.01431
```

The def files of the lipids supported by **buildH** can be found [here](https://github.com/patrickfuchs/buildH/tree/master/def_files).

More on def files and creating your own ones can be found [here](def_format.md).

### Supported lipids

The list of supported lipids can be requested with `buildH -h`. This command will throw a detailed help to the screen, the list will be indicated at the last line. If you want to analyze a lipid that is not present in **buildH**, you will have to create your own def file as well as a json file which explains to **buildH** how the hydrogens will be reconstructed. This user json file is passed with option `-lt`. Here is more documentation on how to [create your own def file](def_format.md) and how to [create your own json file](json_format.md).

### What about polar hydrogens?

When a lipid contains polar hydrogens, such as the 3 Hs of ethanolamine in PE or the H of the hydroxyl group in cholesterol, these Hs are handled explicitely by the force field. Thus they already exist in the input pdb (and possibly xtc) given as input to **buildH**. In this case, those Hs will be ignored by **buildH** and no order parameter will be calculated on these ones. Usually, these Hs are exchangeable and we do not have experimental order parameters for them. If an output trajectory is requested, **buildH** will just copy the coordinates of these Hs as it does for the heavy atoms.

There is an exception for the force field CHARMM36UA. In this force field, only the apolar Hs of the sn-1 and sn-2 aliphatic tails are in a united-atom representation (starting from the 3rd carbon up to the end of the chain). The other apolar Hs (choline, glycerol, second carbon of sn-1 and sn-2) are explicit. Thus for these latter, **buildH** will ignore them as it does for polar Hs as explained above. Again, if an output pdb (or xtc) is requested, those Hs will be copied to the output pdb and xtc files.

### Mixtures of lipids

If you have a mixture of lipids, you will have to run **buildH** for each lipid separately. If you request an output trajectory, this will have to be done iteratively as well. A guided example on a POPC/POPE mixture can be found in [Notebook03](notebooks/Notebook_03_mixture_POPC_POPE.ipynb). Another one on a POPC/cholesterol mixture can be found in [Notebook05](notebooks/Notebook_05_mixture_POPC_cholesterol.ipynb).

### Periodic boundary conditions

Sometimes, when performing MD, some molecules are split over periodic boundary conditions (PBC). **buildH** takes as input whole structures (pdb, gro, xtc, etc.). If broken molecules are supplied, it will most likely generate nonsense results. So it is up to the user to take care of making molecules whole before running **buildH** (e.g. by using a tool like [trjconv](https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html) in GROMACS with flag `-pbc mol`).


### BuildH as a module

buildH is intended to be used mainly in the Unix command line. It is also possible to use it as a module but to a lesser extent.
The features available are minimal: you can just call the main function (`buildh.launch()`) and result files are still written.

It's not a proper API but more a way to call buildH inside larger analysis python scripts.

A guided example can be found on [Notebook04](notebooks/Notebook_04_library.ipynb).

# Lipid json file format

To build new hydrogens on a united-atom lipid, we need different informations that are read by **buildH** in a json file. By default, some standard lipids are present in **buildH** (in the directory [`buildh/lipids/`](https://github.com/patrickfuchs/buildH/tree/master/buildh/lipids)). However, it is possible for the user to supply his/her own json file. Here we explain the format of these json files. All images in this page were generated with [Pymol](https://www.pymol.org/).

## Generality on the lipid json format

The convention for naming the json file is `Forcefield_Lipid.json`. `Forcefield` obviously specifies the force field and `Lipid` is the residue name of the lipid in the pdb or gro file. One such name can be for example `Berger_POPC.json`.

The format of the json file ressembles a Python dictionnary. For example, if we look at the `Berger_POPC.json` file (already present in **buildH**), we have the following:

```
{
  "resname": ["POPC", "PLA", "POP"],
  "C1": ["CH3", "N4", "C5"],
  "C2": ["CH3", "N4", "C5"],
  "C3": ["CH3", "N4", "C5"],
  "C5": ["CH2", "N4", "C6"],
  "C6": ["CH2", "C5", "O7"],
  "C12": ["CH2", "O11", "C13"],
  "C13": ["CH", "C12", "C32", "O14"],
[...]
  "C24": ["CHdoublebond", "C23", "C25"],
  "C25": ["CHdoublebond", "C24", "C26"],
  "C26": ["CH2", "C25", "C27"],
[...]
  "C48": ["CH2", "C47", "C49"],
  "C49": ["CH2", "C48", "C50"],
  "C50": ["CH3", "C49", "C48"]
  }
```

Each lines has a `key: value` pattern as in a Python dictionnary. The `value` ressembles a Python list `[value1, value2, ...]`. Each couple `key: value` is separated by a comma. At the beginning and end of the file we have curly braces `{}`. 

**Important**: note that in the last line (atom `"C50"`), **the comma is not present** at the end of the line.

The first line with a key `"resname"` indicates some possible residue names for the lipid described in the file. Here for example, it can be called `"POPC"`, `"PLA"` or `"POP"` in the pdb or gro file. Do not forget the quotes for each element. Thanks to this line, one can then use `Berger_POPC`, `Berger_PLA`, or `Berger_POP` with the `-l` argument when launching **buildH** at the Unix command line.

In the the next lines, each `key` is basically a carbon atom name (between quotes) on which one wants to reconstruct hydrogens. This is the same atom name as found in the pdb or gro file. The corresponding `value` is a list containing 3 or 4 strings separated by a comma:

- The first string can be either `"CH"`, `"CH2"`, `"CH3"` or `"CHdoublebond"`. It indicates to buildH if we want to reconstruct one H, 2 Hs, 3 Hs (sp3 carbon) or one H of a double bond (sp2 carbon) respectively. In fact, it represents the type of carbon on which we want to build new hydrogens.
- The next strings are called helpers (see below) and are atom names between quotes. We have 2 helpers for `"CH2"`, `"CH3"` or `"CHdoublebond"`, and 3 helpers for `"CH"`.

So the general syntax is `[type of C on which we want to build Hs, name of helper1, name of helper2, ...]`. The choice of helpers and their order in the json file depends on the type of carbon. Everything is described below.

## CH3 

In the figure below is shown a resconstruction of 3 hydrogens (methyl group) on atom `C1`. In the json file, it corresponds to the line `"C1": ["CH3", "N4", "C5"],`. The first helper (helper1) needs to be the one connected to `C1` (thus `N4`), and the second helper (helper2) is connected to `N4` and 2 atoms away from `C1` (thus `C5`) along the main chain. 

![Reconstruction of a CH3](img/build_CH3.png)

The names of the 3 reconstructed H, here `H11`, `H12` and `H13`, are infered from the def file supplied with option `-d`. In this example, [Berger_POPC.def](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def) was used, the names for the 3 H of `C1` were read from the 4th column of these lines:

```
gamma1_1 POPC C1  H11
gamma1_2 POPC C1  H12
gamma1_3 POPC C1  H13
```

In this example, the use of `C2` or `C3` as helper2 would have worked too. However, we decided to use `C5` because it stands along the main chain of the lipid.

## CH2

In the figure below is shown the resconstruction of 2 hydrogens (methylene group) on atom `C26`. On the left is shown a CH2 reconstruction coming from the line `"C26": ["CH2", "C25", "C27"],` in the json file. `"CH2"` means we want to reconstruct 2 hydrogens, `"C25"` is helper1, `"C27"` is helper2. With `C25` being up and `C27` being down, the new hydrogens reconstructed are arranged in space so that `H261` comes towards us and `H262` goes backwards.

On the right, we show the other case where we swapped the order of helper1 and helper2. One can see that the two reconstructed hydrogens are also swapped. **This shows that the order of helpers matters** for H reconstruction on CH2!

![Reconstruction of a CH2](img/build_CH2.png)

The names `H261` or `H262` come again from the 4th column of the def file ([Berger_POPC.def](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def)). The relevant lines are:

```
oleoyl_C11a POPC C26  H261
oleoyl_C11b POPC C26  H262
```

**TODO**: tell which H is pro-R and pro-S.

## CH

For a CH, we want to reconstruct a single hydrogen on a carbon connected to 3 other heavy atoms. In this case, the carbon can be asymetric. This is the case, for example, in phospholipids for the second carbon of the glycerol as shown in the figure below. There is shown the resconstruction of a unique hydrogen on atom `C13`. In this case, we have 3 helpers which are merely the 3 heavy atoms (`C12`, `C32`, and `O14`) connected to that carbon. Note that the order of helpers in the json file `"C13": ["CH", "C12", "C32", "O14"],` does not matter in this case. `"C12"`, `"C32"` and `"O14"` can be put in any order in this list, the H reconstruction will be strictly identical.

![Reconstruction of a CH](img/build_CH.png)

The name `H131` comes again from the 4th column of the def file ([Berger_POPC.def](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def)). The relevant line is:

```
g2_1 POPC C13 H131
```

## CH of a double bond

When a carbon is involved in a double bond, we want to reconstruct a single H which respects the sp2 geometry. Below is shown an example on which the double bond stands between `C24` and `C25` and we want to build the single H on `C25`. The line in the json file for such a case is `"C25": ["CHdoublebond", "C24", "C26"],`, where the first string in the list is now `"CHdoublebond"`. The two helpers are `C24` and `C26` which are the atoms directly bonded to `C25`. Note that the order of helpers in the list does not matter in this case, `"C24", "C26"` or `"C26", "C24"` will work the same.

![Reconstruction of a CH in a double bond](img/build_CHdoublebond.png)

The name `H251` comes again from the 4th column of the def file ([Berger_POPC.def](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def)). The relevant line is:

```
oleoyl_C10a POPC C25  H251
```

## Summary

We have explained here the format of the json file which tells buildH what are the carbons on which we want to reconstruct hydrogens and the corresponding helpers. Finally, we draw again your attention about the order of helpers within the json file:

- the order of helpers in each list does not matter in the case of a CH or CHdoublebond;
- the order of helpers in each list **does matter** in the case of a CH3 or CH2 reconstruction:
  - for CH3, helper1 is bonded to the carbon on which we reconstruct hydrogens and helper2 is two atoms away;
  - for CH2, both helper1 and helper2 are bonded to the carbon on which we reconstruct hydrogens, but their order will determine which reconstructed H is pro-R or pro-S.

## A guided example for writing a lipid json file

We show here how to build your own json file on the simple molecule of butane. We start with a pdb file of the molecule `butane.pdb`:

```
ATOM      1  C1  BUTA    1      -1.890   0.170   0.100  1.00  0.00
ATOM      2  C2  BUTA    1      -0.560  -0.550  -0.100  1.00  0.00
ATOM      3  C3  BUTA    1       0.540   0.520  -0.110  1.00  0.00
ATOM      4  C4  BUTA    1       1.910  -0.140   0.100  1.00  0.00
```

![Butane without hydrogen](img/butane.png)

Now we need to build the json file. According to the rules above, we have the following:

- `C1` is of type CH3. Helper1 is connected to `C1` (thus `C2`), helper2 is two atoms away (thus `C3`).
- `C2` is of type CH2. Helper1 is the carbon before in the chain (thus `C1`), helper2 is the atom after in the chain (thus `C3`).
- `C3` is also of type CH2, so following the same rule, helper1 is `C2` and helper2 is `C4`. 
- `C4` is also of type CH3, so following the same rule, helper1 is `C3` and helper2 is `C2`. 

In summary, this would give the following file:

```
{
    "resname": ["BUTA", "BUT"],
    "C1": ["CH3", "C2", "C3"],
    "C2": ["CH2", "C1", "C3"],
    "C3": ["CH2", "C2", "C4"],
    "C4": ["CH3", "C3", "C2"]
}
```

We have to name it with the convention `Forcefield_Residue.json`. So let us imagine we use Berger force field, we can call it `Berger_BUTA.json`.

We also have to create the def file (see [here](def_format.md#a-guided-example-for-writing-a-lipid-def-file) on how to do that). We can use the following `Berger_BUTA.def`:

```
butane_C1a BUTA C1 H11
butane_C1b BUTA C1 H12
butane_C1c BUTA C1 H13
butane_C2a BUTA C2 H21
butane_C2b BUTA C2 H22
butane_C3a BUTA C3 H31
butane_C3b BUTA C3 H32
butane_C4a BUTA C4 H41
butane_C4b BUTA C4 H42
butane_C4c BUTA C4 H43
```

With those 3 files, we can launch **buildH**:

```bash
buildH -c butane.pdb -l Berger_BUTA -lt Berger_BUTA.json -d Berger_BUTA.def -opx butane_wH
```

So we used here the option `-lt` to supply our own `Berger_BUTA.json` file. We also requested an ouput with option `-opx` which will generate the pdb with hydrogens `butane_wH.pdb`. Below is shown the generated pdb and structure.

```
ATOM      1  C1  BUTA    1      -1.890   0.170   0.100  1.00  0.00             C
ATOM      2  H11 BUTA    1      -2.700  -0.560   0.113  1.00  0.00             H
ATOM      3  H12 BUTA    1      -2.048   0.874  -0.717  1.00  0.00             H
ATOM      4  H13 BUTA    1      -1.872   0.710   1.047  1.00  0.00             H
ATOM      5  C2  BUTA    1      -0.560  -0.550  -0.100  1.00  0.00             C
ATOM      6  H21 BUTA    1      -0.566  -1.088  -1.048  1.00  0.00             H
ATOM      7  H22 BUTA    1      -0.390  -1.253   0.716  1.00  0.00             H
ATOM      8  C3  BUTA    1       0.540   0.520  -0.110  1.00  0.00             C
ATOM      9  H31 BUTA    1       0.356   1.235   0.692  1.00  0.00             H
ATOM     10  H32 BUTA    1       0.531   1.039  -1.069  1.00  0.00             H
ATOM     11  C4  BUTA    1       1.910  -0.140   0.100  1.00  0.00             C
ATOM     12  H41 BUTA    1       2.687   0.625   0.092  1.00  0.00             H
ATOM     13  H42 BUTA    1       2.096  -0.855  -0.702  1.00  0.00             H
ATOM     14  H43 BUTA    1       1.920  -0.658   1.059  1.00  0.00             H
```

The pdb name of the newly built hydrogens (`H11`, `H12`, etc.) were infered from the 4th column of the def file.

![Butane with hydrogens](img/butane_wH.png)

**Last advices**

- We showed you a simple example on butane. Although this molecule is very simple, you can see that it is easy to make a mistake. So we recommend to triple check the json file before using it for production. **buildH** makes for you a lot of checks and will throw an error if something is wrong, but it cannot detect all types of mistakes. Any spelling error on atom names, inversion, etc., may lead to aberrant results. So before going to production, do test on a single molecule and check thoroughly the molecule has all the hydrogens in good place. 

- The main lipids are already included in **buildH** (in the directory `buildh/lipids`) so you might not need to build your own json file. You can have a list of the supported lipids by invoking **buildH** with option `-h`:

```
$ buildH -h
usage: buildH [-h] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]] -d DEFOP
[...]
The list of supported lipids (-l option) are: Berger_CHOL, Berger_DOPC, Berger_DPPC, Berger_PLA, Berger_POPC, Berger_POP, Berger_POPE, Berger_POPS, CHARMM36UA_DPPC, CHARMM36UA_DPUC, CHARMM36_POPC, GROMOS53A6L_DPPC, GROMOSCKP_POPC, GROMOSCKP_POPS. More documentation can be found at https://buildh.readthedocs.io.
```

- Last, one other project developed by us, called [autoLipMap](https://github.com/patrickfuchs/autoLipMap), can build automatically def and json files for the main known lipids.

In case of problem, you can post an issue on github.
# Installation

## Requirements and compatibility

buildH requires at least Python 3.6 and needs the following modules:
  - numpy
  - pandas
  - numba
  - MDAnalysis (with support of 2.0)

All the instructions below have been tested on Unix like platforms (e.g. Ubuntu), which we recommend for running **buildH**. We do not provide support for other platforms, but since **buildH** has been written in pure Python, it should work there provided its dependencies are supported. 

## Simple installation

A simple installation with pip will do the trick:

```
python3 -m pip install buildh
```

All dependencies (modules) will be installed automatically by pip.

Note that this way of proceeding will install **buildH** and its dependencies within the python of your Unix system, which may lead to conflicts of version if you have other scientific packages installed. To avoid this you may want to create a specific conda or virtual environment for **buildH** (see below).

## Installation within a conda environment

**buildH** is also available through the [Bioconda](https://anaconda.org/bioconda/buildh) channel.

We recommend to install **buildH** within a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). First create a new conda env:

```
conda create -n buildH "python>=3.6"
```

Then activate your environment:

```
conda activate buildH
```

Last, install **buildH** within that environment:

```
conda config --add channels conda-forge
conda config --add channels bioconda
conda install buildh
```

## Building from source

We recommend to use a specific environment, either by using [venv](https://docs.python.org/3/library/venv.html) or [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands). All packages will be installed within that environment which avoids conflicts of version.

If you still do not want to create a specific environment for **buildH**, you can skip the first section `Create an environment` below.

In any case, the python version should be >= 3.6.

### Create an environment

First, create a new environment (we call it `env4buildH`):

- If you chose conda: `conda create -n env4buildH python=3.8`
- If you chose venv: `python3 -m venv /path/to/env4buildH`

Activate your environment:

- If you chose conda: `conda activate env4buildH`
- If you chose venv: `source /path/to/env4buildH/bin/activate`

### Install **buildH** from source

Clone the **buildH** repository:

```
git clone https://github.com/patrickfuchs/buildH.git
cd buildH
```

Install with `pip` the packages required by **buildH**, namely numpy, pandas, MDAnalysis and Numba, which are all specified in the file [requirements.txt](https://github.com/patrickfuchs/buildH/blob/master/requirements.txt):

```
pip install -r requirements.txt
```

Install **buildH** from source with `pip`:

```
pip install -e .
```

## For developers

For installing a development version from source with the full environment (allowing building the doc, launching tests, etc.), see [here](https://github.com/patrickfuchs/buildH/tree/master/devtools/install_dev.md).

## Testing

The tests rely on the `pytest` package. You need first to install a development version (see above). Once done, you can run the tests with just:

```
cd buildh # if it's not already the case
pytest
```

All tests should pass. If anything fails, please [open an issue](https://github.com/patrickfuchs/buildH/issues).
# Algorithms for building hydrogens

**buildH** builds hydrogens using general geometric rules which are explained in this document. All the Python functions implementing these reconstructions are written in [`hydrogens.py`](https://github.com/patrickfuchs/buildH/blob/master/buildh/hydrogens.py). These functions are largely inspired from [a code of Jon Kapla](https://github.com/kaplajon/trajman/blob/master/module_trajop.f90#L242) originally written in fortran. All mathematical functions (vector operations, rotations, etc.) are written in [`geometry.py`](https://github.com/patrickfuchs/buildH/blob/master/buildh/hydrogens.py) and accelarated using [Numba](https://numba.pydata.org/). In this page, we use the following conventions:

- the word `helper` describes heavy atoms connected (or not far away) to the carbon on which we reconstruct Hs, their position is used in the process of H reconstruction;
- all represented vectors are unit vectors;
- $\theta$ is the [tetrahedral bond angle](https://en.wikipedia.org/wiki/Tetrahedron) which equals arccos(-1/3) rad ~ 1.9106 rad ~ 109.47°; 
- all the angles will be described in rad;
- $l_{CH}$ is the [carbon-hydrogen bond length](https://en.wikipedia.org/wiki/Carbon%E2%80%93hydrogen_bond) which equals ~ 1.09 Å.

All images in this page were generated using [VMD](http://www.ks.uiuc.edu/Research/vmd/).

## Building CH3

Building a methyl on a primary carbon (`C1`) requires two helpers: i) helper1 (`N4`) is connected to that carbon, ii) helper2 (`C5`) is connected to helper1. As described [here](json_format.md#ch3), helper1 has to be a heavy atom connected to the primary carbon, helper2 is a heavy atoms 2 atoms away from the primary carbon. In some cases, there may be multiple choices for helper2 (such as in the CH3s of choline in PC lipids). Any choice is fine, as long as the rule explained here is followed (helper1 connected to the carbon, helper2 is 2 atoms away).

We start with the first hydrogen reconstruction as explained in the figure below.

![CH3_building](img/how_CH3_building1.png)

`vect1` (red) is first computed as the vector product between vectors "carbon -> helper2" and "carbon -> helper1", it will be our rotation axis in the next step. `vect2` (blue) is then computed by rotating vector "carbon -> helper1" about `vect1` by $\theta$ rad. The first H will be obtained by translating a point located at the "carbon" along `vect2` of $l_{CH}$ Å. Note that the newly built H is in a *trans* configuration with respect to helper2.

We then go on with the reconstruction of the two other Hs as explained in the figure below.

![CH_building](img/how_CH3_building2.png)

`vect3` (green) is obtained by rotating `vect2` of $+\frac{2\pi}{3}$ rad about vector "carbon -> helper1". 
`vect4` (magenta) is obtained by rotating `vect2` of $-\frac{2\pi}{3}$ rad about vector "carbon -> helper1". 

The second and third H will be obtained by translating of $l_{CH}$ Å a point located at the "carbon" along `vect3` and `vect4` respectively.

## Building CH2

The building of 2 hydrogens on a secondary carbon involves a few geometrical procedures that are explained in the figure below.

![CH2_building](img/how_CH2_building.png)

We start with the 3 atoms, the central carbon on which we want to reconstruct hydrogens (`C26`), helper1 (`C25`) and helper2 (`C27`) which are heavy atoms connected to the central carbon. The two helpers will help us build the new hydrogens following standard [tetrahedral geometry](https://en.wikipedia.org/wiki/Tetrahedral_molecular_geometry). 

On the left panel, we first show how to construct 3 vectors:

- `vect1` (red) is normal to the plane of the 3 atoms. It is calculated as the cross product between vectors "central carbon -> helper2" and "central carbon -> helper1".
- `vect2` (blue) will be our **rotation axis** used later. It is calculated as vector "central carbon -> helper1" minus vector "central carbon -> helper2".
- `vect3` (green) is a vector that will be rotated in the next step. It is the cross product between `vect1` and `vect2`.

On the right panel, we go on to construct 2 other vectors:

- `vect4` (magenta) is obtained by rotating `vect3` of $\frac{\theta}{2}$ rad about `vect2`. The first H will be obtained by translating a point located at the central carbon along `vect4` of $l_{CH}$ Å.
- `vect5` (orange) is obtained by rotating `vect3` of $-\frac{\theta}{2}$ rad about `vect2`. The second H will be obtained by translating a point located at the central carbon along vect5 of $l_{CH}$ Å.

**Important**: last, we want to stress that the **order of helpers** matters as described [here](json_format.md#ch2). If the two helpers are switched, the two reconstructed Hs will be switched too, so this will have a consequence on which one is pro-R and pro-S.

## Building CH

The building of 1 hydrogen on a tertiary carbon is quite simple. This carbon can be asymetric. In a POPC lipid there is a single case of tertiary carbon on the second carbon of glycerol. In natural phospholipids, this carbon is in *R* configuration. 

The figure below  shows an example using a Berger POPC. The central carbon on which we want to reconstruct one H is the second carbon of the glycerol (`C13`). Helper1 (`C13`) is the third glycerol carbon (towards sn-3, that is, the polar head),  helper2 (`C32`) is the fist glycerol carbon (towards sn-1) and helper3  (`O14`) is the glycerol oxygen towards sn-2. Importantly, the order of helpers (what is helper1, helper2 and helper3) does not matter in this case, any combination will lead to the same H resconstruction.

![CH_building](img/how_CH_building.png)

We first compute the red vector `vect1` as the sum of the 3 vectors "central carbon -> helper1" + "central carbon -> helper2" + "central carbon -> helper3". We see that this `vect1` defines a [median](https://en.wikipedia.org/wiki/Median_(geometry)#Tetrahedron) of the tetrahedron. The blue vector `vect2` is merely the opposite of `vect1` and gives the direction of the C-H bond. The new H will be obtained by translating a point located at the central carbon along `vect2` of $l_{CH}$ Å.

## Building CH on a double bond

For the H to reconstruct on a carbon involved in a double bond, **buildH** uses the following strategy as explained in the figure below.

![CH_building](img/how_CHdoublebond_building.png)

First we compute the angle $\gamma$ between atoms "helper1-central carbon-helper2". `vect1` is next calculated as the cross product between vectors "central carbon -> helper1" and "central carbon -> helper2". `vect1` will be used in the next step as a rotation axis. `vect2` is finally obtained by rotating vector "central carbon -> helper2" of $\pi-\frac{\gamma}{2}$ rad about `vect1`. Why this angle value? In fact, **buildH** uses here the bisection strategy. `vect2` is along the same axis as the vector bisecting the angle "helper1-central carbon-helper2" but on the opposite direction. To obtain this bisecting vector, we would need to rotate vector "central carbon -> helper2" of $-\frac{\gamma}{2}$. Since `vect2` is on the opposite direction, we simply add $\pi$ to that value obtaining thus $\pi-\frac{\gamma}{2}$ rad. The new H will be obtained by translating a point located at the central carbon along `vect2` of $l_{CH}$ Å.

Importantly, the order of helpers (what atom is helper1 or helper2) does not matter, any combination will lead to the same H reconstruction.

The double bond case is somewhat more complicated than the other ones (CH3/CH2/CH) as largely discussed in an article published in JCTC by [Piggot at al.](https://doi.org/10.1021/acs.jctc.7b00643). One of the problem is that the ideal CCC angle (helper1-carbon-helper2) of 120° can vary according to the functional groups and some details in the force field. One solution is to adapt the value for each force field and each double bond at play which requires to do some testing for each case. In **buildH** we use the bisection strategy (as described above) so that the H reconstruction will adapt itself to any value of ${\gamma}$. Comparing the results with H reconstructed like this vs the real ones from the oleoyl double-bond of an all-atom CHARMM36 snapshot of 256 POPC yielded a reasonable (averaged) difference in $S_{CH}$ of 0.008 and 0.016 (see [here](https://zenodo.org/record/4715962)). This difference was even reduced upon averaging over a trajectory. In summary, the bisection appears as an acceptable tradeoff.
# Command line options

We explain in this document the details of the different options on the command line. 

## General usage

When **buildH** is invoked with the flag `-h`, it displays some quite detailed help to the screen:

```
usage: buildH [-h] [-v] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]] -d DEFOP [-opx OPDBXTC] [-o OUT] [-b BEGIN] [-e END]
              [-igch3]

This program builds hydrogens and calculates the order parameters (OP) from a united-atom trajectory of lipids. If -opx is requested, pdb and xtc
output files with hydrogens are created but OP calculation will be slow. If no trajectory output is requested (no use of flag -opx), it uses a fast
procedure to build hydrogens and calculate the OP.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -c COORD, --coord COORD
                        Coordinate file (pdb or gro format).
  -t TRAJ, --traj TRAJ  Input trajectory file. Could be in XTC, TRR or DCD format.
  -l LIPID, --lipid LIPID
                        Combinaison of ForceField name and residue name for the lipid to calculate the OP on (e.g. Berger_POPC).It must match with
                        the internal topology files or the one(s) supplied.A list of supported terms is printed when calling the help.
  -lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...], --lipid_topology LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]
                        User topology lipid json file(s).
  -d DEFOP, --defop DEFOP
                        Order parameter definition file. Can be found on https://github.com/patrickfuchs/buildH/tree/master/def_files.
  -opx OPDBXTC, --opdbxtc OPDBXTC
                        Base name for trajectory output with hydrogens. File extension will be automatically added. For example -opx trajH will
                        generate trajH.pdb and trajH.xtc. So far only xtc is supported.
  -o OUT, --out OUT     Output file name for storing order parameters. Default name is OP_buildH.out.
  -b BEGIN, --begin BEGIN
                        The first frame (ps) to read from the trajectory.
  -e END, --end END     The last frame (ps) to read from the trajectory.
  -igch3, --ignore-CH3s
                        Ignore CH3s groups for the construction of hydrogens and the calculation of the OP.

The list of supported lipids (-l option) are: Berger_CHOL, Berger_DOPC, Berger_DPPC, Berger_POPC, Berger_PLA, Berger_POP, Berger_POPE, Berger_POPS,
CHARMM36UA_DPPC, CHARMM36UA_DPUC, CHARMM36_POPC, GROMOS53A6L_DPPC, GROMOSCKP_POPC, GROMOSCKP_POPS. More documentation can be found at
https://buildh.readthedocs.io.
```

Importantly, the `-h` option also displays the list of supported lipids at the end. Note that they are always written using the naming convention `ForceField_Lipid`.

## Description of options

Each option is explained below. Some of them are mandatory.

### Help

`-h` or `--help`: display the help message shown above and exit.

(**optional flag**)

### Coordinates

`-c COORD` or `--coord COORD`: `COORD` is the main coordinate file in pdb or gro format describing the system. This is strictly required.

(**mandatory flag**)

### Trajectory

`-t TRAJ` or `--traj TRAJ`: `TRAJ` is an input trajectory file in XTC, TRR or DCD format. If not provided, the H reconstruction and order parameter calculation will be done solely on the `COORD` file. If a trajecotry is provided with `-t` flag, the resulting order parameters will be averaged over that trajectory. More on how this averaging is performed can be found [here](https://buildh.readthedocs.io/en/latest/buildh.html#order-parameters-and-statistics).

(**optional flag**)

### Lipid requested

`-l LIPID` or `--lipid LIPID`: `LIPID` is the name of the lipid to calculate the OP on. It must follow the naming convention `ForceField_Lipid`, for example `Berger_POPC`. The list of supported lipids can be queried with the `-h` flag.

(**mandatory flag**)

### User lipid (json file)

`-lt LIPID_TOPOLOGY` or `--lipid_topology LIPID_TOPOLOGY`: `LIPID_TOPOLOGY` is a user supplied topology lipid json file(s). When you want to analyze a lipid not present in **buildH** you can [build your own json file](json_format.md) and supply it with this option. Again, it has to follow the naming convention `ForceField_Lipid.json`. For example, if you build your own json file for butane with the Berger force field, you can use `-lt Berger_BUTA.json`; in this case, you will have to use also `-l Berger_BUTA` flag. 

(**optional flag**)

### Definition file

`-d DEFOP` or `--defop DEFOP`: `DEFOP` is the order parameter definition file. It tells **buildH** what C-H will be considered for H reconstruction and order parameter calculation. You can find some [def files](https://github.com/patrickfuchs/buildH/tree/master/def_files) on the **buildH** repository for the supported lipids. You can also create your [own def file](def_format.md).

(**mandatory flag**)

### Output trajectory

`-opx OPDBXTC` or `--opdbxtc OPDBXTC`: if you want a trajectory with all hydrogens reconstructed, `OPDBXTC` is a base name. If you supply `-opx traj_wH`, two files will be created: `traj_wH.pdb` and `traj_wH.xtc` both containing all hydrogens. File extension will be automatically added. If no trajectory is supplied with the option `-t`, **buildH** will only create a pdb but not an xtc. So far only pdb and xtc are supported. Note that this option is slow.

(**optional flag**)

### Output name for the order parameter

`-o OUT` or `--out OUT`: `OUT` will be the output name for storing the calculated order parameters. If this option is not supplied, the default output name is `OP_buildH.out`.

(**optional flag**)

### Precising beginning and end of trajectory

`-b BEGIN` or `--begin BEGIN`: The first frame (in ps) to read from the trajectory.

`-e END` or `--end END`: The last frame (in ps) to read from the trajectory.

**buildH** checks whether the `BEGIN` and `END` make sense with the supplied trajectory.

(**optional flag**)

### Ignoring CH3 in the output

`-igch3` or `--ignore-CH3s`: when this flag is on, the OP of each C-H belonging to any CH3 will not be computed and not be written in the ouput file even if they are present in the def file. Individual reconstruted C-H of methyl groups are not very interesting for OP calculation since we cannot precisely know where they are because of the methyl rotation. With this option, one can avoid their evaluation without having to remove these C-H in the def file (recall, the [def files from the **buildH** website](https://github.com/patrickfuchs/buildH/tree/master/def_files) contain all possible C-Hs). However, note that this `-igch3` option is not usable with the `-opx` option (which outputs the structure pdb and xtc files with hydrogens) since this latter needs all hydrogen atoms to reconstruct.
# Lipid def file format

## Role of def files

The def file has three main functions:

- Tell **buildH** which C-H bonds are considered for H reconstruction and order parameter calculation.
- Give a generic name to each C-H for the output of the order parameter (the default name for that file is `OP_buildH.out`); a generic name can be for example `beta1` or `beta2`.
- Give a specific pdb name to each newly built hydrogen when a trajectory is requested.

Initially, this type of file has been created in the [NMRlipids](https://nmrlipids.blogspot.com/) project. Starting from all-atom trajectories, the script [calcOrderParameters.py](https://github.com/NMRLipids/MATCH/blob/master/scripts/calcOrderParameters.py) takes such def files as input to infer what C-H are considered for the calculation. But it is also useful for assigning a unique generic name for each C-H along the lipid regardless of the force field considered. For example, the two C-H of of the beta carbon of the choline are called `beta1` and `beta2` in all PC lipids whatever the force field. Last, it is also useful to specify a specific pdb name to each reconstructed hydrogen.

Many of these `.def` files can be found on the [MATCH repository](https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs) for all-atom force fields. As in **buildH** we work with united-atom trajectories, there is a directory [def_files](https://github.com/patrickfuchs/buildH/tree/master/def_files) where you will find all the def files for the lipids supported.

## Format of def files

The general naming convention of def files on NMRLipids is `order_parameters_definitions_MODEL_X_Y.def`, where `X` is the force field and `Y` the lipid considered. You can name it whatever you like, but we recommend to indicate at least the force field and the lipid. In **buildH**, the [def files](https://github.com/patrickfuchs/buildH/tree/master/def_files) present in the repo have mere names folllowing the convention `Forcefield_Lipid.def`, for example `Berger_POPC.def`.

Each line represents a given C-H with 4 columns:

- Column 1 is the generic name of the order parameter considered. We recommend to use the same convention as those already present in the [def_files](https://github.com/patrickfuchs/buildH/tree/master/def_files) provided in **buildH**.
- Column 2 is the residue name in the pdb or gro file.
- Column 3 is the carbon name in the pdb or gro file for that C-H.
- Column 4 is the H name for that C-H, it is used if an output pdb file is requested.

Here is an example for the polar head of [Berger POPC](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def):

```
beta1  POPC C5  H51
beta2  POPC C5  H52
alpha1 POPC C6  H61
alpha2 POPC C6  H62
g3_1   POPC C12 H121
g3_2   POPC C12 H122
g2_1   POPC C13 H131
g1_1   POPC C32 H321
g1_2   POPC C32 H322
```

Each column has to be separated by any combination of white space(s) (at least one).

## Trajectory output and def files

If an output trajectory (option `-opx`) is requested, the `.def` file provided **must** contain all possible C-H in that lipid (since the whole trajectory with Hs will be reconstructed). The newly built hydrogens will follow the names written in the 4th column of the def file. This option is slow, we do not recommend it if an output xtc file is not wanted.

If no option `-opx` is used, **buildH** uses fast indexing. In this case the `.def` file can contain any subset of all possible C-H pairs. For example, if one wants to get the order parameters of the polar head only (Berger POPC), the `.def` will be the following:

```
beta1  POPC C5  H51
beta2  POPC C5  H52
alpha1 POPC C6  H61
alpha2 POPC C6  H62
g3_1   POPC C12 H121
g3_2   POPC C12 H122
g2_1   POPC C13 H131
g1_1   POPC C32 H321
g1_2   POPC C32 H322
```

Of course, a lower number of Hs to reconstruct will make **buildH** run faster.

## A guided example for writing a lipid def file

We show here how to build your own def file on the simple molecule of butane. We start with a pdb file of the molecule `butane.pdb`:

```
ATOM      1  C1  BUTA    1      -1.890   0.170   0.100  1.00  0.00
ATOM      2  C2  BUTA    1      -0.560  -0.550  -0.100  1.00  0.00
ATOM      3  C3  BUTA    1       0.540   0.520  -0.110  1.00  0.00
ATOM      4  C4  BUTA    1       1.910  -0.140   0.100  1.00  0.00
```

![Butane without hydrogen](img/butane.png)

Let us consider first `C1` which has 3 hydrogens to reconstruct. We can call the new hydrogens `H11`, `H12` and `H13`. The residue name is `BUTA`. Since the butane molecule does not exist in known def files, we can create any generic name we want for each C-H. Following the same model as [Berger_POPC.def](https://github.com/patrickfuchs/buildH/blob/master/def_files/Berger_POPC.def), we could name them `butane_C1a`, `butane_C1a` and  `butane_C1c`. Applying this to all possible C-H, we obtain:

```
butane_C1a BUTA C1 H11
butane_C1b BUTA C1 H12
butane_C1c BUTA C1 H13
butane_C2a BUTA C2 H21
butane_C2b BUTA C2 H22
butane_C3a BUTA C3 H31
butane_C3b BUTA C3 H32
butane_C4a BUTA C4 H41
butane_C4b BUTA C4 H42
butane_C4c BUTA C4 H43
```

Following the file naming convention `Forcefield_Lipid.def`, this file should be called `Berger_BUTA.def`.

We also have to create the json file (see [here](json_format.md#a-guided-example-for-writing-a-lipid-json-file) on how to do that). We can use the following `Berger_BUTA.json`:

```
{
    "resname": ["BUTA", "BUT"],
    "C1": ["CH3", "C2", "C3"],
    "C2": ["CH2", "C1", "C3"],
    "C3": ["CH2", "C2", "C4"],
    "C4": ["CH3", "C3", "C2"]
}
```

With those 3 files, we can launch buildH:

```bash
buildH -c butane.pdb -l Berger_BUTA -lt Berger_BUTA.json -d Berger_BUTA.def -opx butane_wH
```

We supplied our newly created def file to **buildH** with the option `-d`. We also used our own `Berger_BUTA.json` file with the option `-lt`. We requested an ouput with the option `-opx` which will generate the pdb with hydrogens `butane_wH.pdb`. Below is shown the generated pdb and structure.

```
ATOM      1  C1  BUTA    1      -1.890   0.170   0.100  1.00  0.00             C
ATOM      2  H11 BUTA    1      -2.700  -0.560   0.113  1.00  0.00             H
ATOM      3  H12 BUTA    1      -2.048   0.874  -0.717  1.00  0.00             H
ATOM      4  H13 BUTA    1      -1.872   0.710   1.047  1.00  0.00             H
ATOM      5  C2  BUTA    1      -0.560  -0.550  -0.100  1.00  0.00             C
ATOM      6  H21 BUTA    1      -0.566  -1.088  -1.048  1.00  0.00             H
ATOM      7  H22 BUTA    1      -0.390  -1.253   0.716  1.00  0.00             H
ATOM      8  C3  BUTA    1       0.540   0.520  -0.110  1.00  0.00             C
ATOM      9  H31 BUTA    1       0.356   1.235   0.692  1.00  0.00             H
ATOM     10  H32 BUTA    1       0.531   1.039  -1.069  1.00  0.00             H
ATOM     11  C4  BUTA    1       1.910  -0.140   0.100  1.00  0.00             C
ATOM     12  H41 BUTA    1       2.687   0.625   0.092  1.00  0.00             H
ATOM     13  H42 BUTA    1       2.096  -0.855  -0.702  1.00  0.00             H
ATOM     14  H43 BUTA    1       1.920  -0.658   1.059  1.00  0.00             H
```

![Butane with hydrogens](img/butane_wH.png)

**Last advices**

- We showed you a simple example on butane. Although this molecule is very simple, you can see that it is easy to make a mistake. Although this def file is less critical than the [json lipid file](json_format.md), we recommend to double check it before using it for production. **buildH** makes for you a lot of checks and will throw an error if something is wrong, but it cannot detect all types of mistakes. Any spelling error on atom names, inversion, etc., may lead to nonsense results. So before going to production, do test on a single molecule and check thoroughly on a couple of examples whether the output looks right.

- The main lipids are already included in **buildH** (in the directory [`buildh/lipids`](https://github.com/patrickfuchs/buildH/tree/master/def_files)) so you might not need to build your own def file. You can have a list of the supported lipids by invoking **buildH** with option `-h`:

```
$ buildH -h
usage: buildH [-h] -c COORD [-t TRAJ] -l LIPID [-lt LIPID_TOPOLOGY [LIPID_TOPOLOGY ...]] -d DEFOP
[...]
The list of supported lipids (-l option) are: Berger_POP, Berger_POPC, Berger_PLA, Berger_POPE, CHARMM36_POPC.
[...]
```

- Last, one other project developed by us, called [autoLipMap](https://github.com/patrickfuchs/autoLipMap), can build automatically def and json files for the main known lipids.

In case of problem, you can post an issue on github.
# Order parameters and statistics

The mean order parameter of bond $CH_j$ (i.e. the $j^{th}$ C-H bond) is calculated using the standard formula:

$$\overline{S_{CH_j}} = \frac{1}{2} \left \langle 3cos^2(\theta) -1 \right \rangle$$

where $\theta$ is the angle between the $CH_j$ bond and the normal to the membrane (usually the *z* axis), $\langle ... \rangle$ means averaging over molecules and frames. $S_{CH}$ can be measured by NMR which is useful to validate simulation results, as largely described in the [NMRlipids project](http://nmrlipids.blogspot.com).

The order parameter output of buildH (default name `OP_buildH.out`) looks like this:

```
# OP_name            resname atom1 atom2  OP_mean OP_stddev OP_stem
#--------------------------------------------------------------------
gamma1_1             POPC    C1    H11    0.01304  0.12090  0.01069
gamma1_2             POPC    C1    H12    0.00666  0.09279  0.00820
gamma1_3             POPC    C1    H13   -0.01531  0.09141  0.00808
[...]
```

Each line corresponds to a given CH. The 4 first columns contain the generic name, residue name, carbon and hydrogen names respectively. The other columns contains different statistics on order parameters (OP):

- `OP_mean`, also written $\overline{S_{CH_j}}$ as described above, is the mean OP of bond $CH_j$ averaged over all lipids and all frames of the trajectory.
- `OP_stddev` is the standard deviation of the OP over residues, we shall write it $\sigma(S_{CH_j})$; first we average each OP of bond $CH_j$ (e.g. the C-H of beta1) of residue $i$ (i.e. lipid $i$) over the whole trajectory:

$$ \overline{S_{CH_j}(i)} = \frac{1}{nframes} \sum_{t=0}^{t=nframes} S_{CH_j}(i)(t) $$

where $nframes$ is the total number of frames, then we calculate the standard deviation of those means over all residues:

$$ \sigma(S_{CH_j}) =
\sqrt{
\frac{1}{nres} \sum_{i=1}^{i=nres} (\overline{S_{CH_j}(i)} - \overline{S_{CH_j}})^2
}$$

where $nres$ is the total number of residues (i.e. lipids).
- `OP_stem` is the standard error of the mean averaged in the same spirit, let's call it $err(S_{CH_j})$:

$$err(S_{CH_j}) = \frac{\sigma(S_{CH_j})}{\sqrt{nres}}$$
# Test case

Here are some examples on how to use buildH with a simple system made of POPC Berger lipids. The file `Berger_POPC.def` comes from the [MATCH repository](https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs).

All output files (`OUT*`) were obtained by lauching buildH in the following way:

```bash

# On a file with a single POPC (1POPC.pdb)
buildH -c 1POPC.pdb -l Berger_POPC \
	-d Berger_POPC.def \
	-o OUT.buildH.1POPC.pdb.out

# On a file with 128 POPC (start_128popc.pdb).
buildH -c start_128popc.pdb -l Berger_POPC \
	-d Berger_POPC.def \
	-o OUT.buildH.start_128popc.pdb.out

# On a small trajectory with 25 frames (popc0-25ns_dt1000.xtc).
buildH -c start_128popc.pdb -l Berger_POPC \
	-d Berger_POPC.def \
	-t popc0-25ns_dt1000.xtc -o OUT.buildH.popc0-25ns_dt1000.xtc.out

```
Here are some files that are used in Notebook 03. They describe an example of the use of **buildH** on a POPC/POPE mixture.
# Validation of buildH 

**buildH** reconstructs hydrogens from a united-atom trajectory and calculates the order parameter on each reconstructed C-H bond. To validate **buildH**, we took an all-atom POPC trajectory generated with the CHARMM36 force field. First, we removed the hydrogens and reconstructed them with **buildH**. Then we compared the H reconstruction and the order parameter values calculated with **buildH** to the real ones from the all-atom trajectory. The output of **buildH** was also compared to two  scripts made by Josef Melcr and Angel Pineiro.

Here is a [report](https://github.com/patrickfuchs/buildH/blob/master/docs/CHARMM36_POPC_validation/report_buildH.pdf) made in August 2019 describing this validation.

All the files used for making this validation have been deposited on [Zenodo](https://zenodo.org/record/4715962) with the following DOI: 10.5281/zenodo.4715962.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4715962.svg)](https://doi.org/10.5281/zenodo.4715962)

Below is shown an animated gif highlighting the tiny difference between the hydrogens in a CHARMM36 POPC vs those reconstructed by **buildH**.

![CHARMM36_vs_buildH.gif](CHARMM36_vs_buildH.gif)

Note that this validation was made in August 2019. You can retrieve the corresponding old version of **buildH** [here](https://github.com/patrickfuchs/buildH/tree/7cf8a331b1758abffd03ebb9737704dee3f12a88).
API Reference
=============

This page provides the documentation generated from the source files.


.. toctree::
   :maxdepth: 1

   ./api/UI
   ./api/core
   ./api/geometry
   ./api/hydrogens
   ./api/init_dics
   ./api/lipids
   ./api/utils
   ./api/writers
Tutorials
=========

.. toctree::
    :maxdepth: 1

    notebooks/Notebook_01_buildH_calc_OP
    notebooks/Notebook_02_buildH_calc_OP_outputwH
    notebooks/Notebook_03_mixture_POPC_POPE
    notebooks/Notebook_04_library
    notebooks/Notebook_05_mixture_POPC_cholesterol
==================================
Welcome to buildH's documentation!
==================================

**Version** |release|

    *Build hydrogen atoms from a united-atom MD of lipids and calculate the order parameters.*

**buildH** is a software that reads a united-atom (UA) trajectory of lipids, builds the hydrogens on it and calculates the order parameter on each C-H bond. **buildH** also allows to output the trajectory with the new reconstructed hydrogens.

**buildH** works in two modes:

1. A slow mode when an output trajectory is requested by the user. In this case, the whole trajectory including newly built hydrogens is written to this trajectory file. If other molecules are present (e.g. water, ions, etc.), they will just be copied to the output trajectory with the same coordinates. So far, only the xtc format is supported.
2. A fast mode without any output trajectory.

In both modes, the order parameters are calculated.

It is possible to select only a part of the lipid on which **buildH** will do his job (e.g. the polar head, the sn-1 aliphatic chain, etc) thanks to the `def file <def_format.md>`_.


**buildH** has been carefully validated as explained in :doc:`Validation of buildH <CHARMM36_POPC_validation/validation>`.
The algorithms used to reconstruct hydrogens are detailed in :doc:`Algorithms for building hydrogens <algorithms_Hbuilding>`
and the formulas for computing the order parameters in :doc:`Order parameters and statistics <order_parameter>`.

All basic geometrical operations in **buildH** are accelerated using `Numba <https://numba.pydata.org>`_. **buildH** is hosted on `Github <https://github.com/patrickfuchs/buildH>`_.

Motivation
==========

The initial motivation comes from the `NMRlipids <https://nmrlipids.blogspot.com/>`_ project.
As stated in this `post <https://nmrlipids.blogspot.com/2019/04/nmrlipids-ivb-assembling-pe-pg-results.html>`_,
there was a lack of suitable program for reconstructing hydrogens.
In the past, we used to use `g_protonate` in GROMACS 3 but this program has been removed in recent versions.

Our idea was to build our own implementation in Python using libraries such as ``MDAnalysis``, ``Numpy`` and ``Pandas``.

**buildH** is used actively in the recent projects of NMRlipids such as `NMRlipidsIVPEandPG <https://github.com/NMRLipids/NMRlipidsIVPEandPG>`_ or `Databank <https://github.com/NMRLipids/Databank>`_.
**buildH** can also be used by anyone willing to analyze the order parameters from a UA trajectory, or if one needs to have explicit hydrogens for some further analyzes.

Citations
=========

If you use buildH, please cite:

    Santuz et al., (2021). buildH: Build hydrogen atoms from united-atom molecular dynamics of lipids and calculate the order parameters. Journal of Open Source Software, 6(65), 3521, https://doi.org/10.21105/joss.03521


License
=======

buildH is licensed under `BSD 3-Clause <https://github.com/patrickfuchs/buildH/blob/master/LICENSE.txt>`_.


Content
=======

User manual
-----------
.. toctree::
    :maxdepth: 2

    installation
    usage
    command_line_options
    def_format
    json_format

Tutorials
---------
.. toctree::
   :maxdepth: 2

   tutorials

Algorithms, OP calculations & validation
----------------------------------------

.. toctree::
   :maxdepth: 2

   algorithms


API documentation
-----------------
.. toctree::
   :maxdepth: 1

   api_reference

Changelog
---------

.. toctree::
   :maxdepth: 1

   changelog
Algorithms & validation
=======================

.. toctree::
   :maxdepth: 1

   algorithms_Hbuilding
   order_parameter
   CHARMM36_POPC_validation/validation
Changelog
=========

.. include:: ../CHANGELOG.md
API reference for `writers` submodule
=====================================

.. automodule:: buildh.writers
    :members:
API reference for `hydrogens` submodule
=======================================

.. automodule:: buildh.hydrogens
    :members:
API reference for `utils` submodule
===================================

.. automodule:: buildh.utils
    :members:
API reference for `geometry` submodule
======================================

.. automodule:: buildh.geometry
    :members:
API reference for `init_dics` submodule
=======================================

.. automodule:: buildh.init_dics
    :members:
API reference for `core` submodule
==================================

.. automodule:: buildh.core
    :members:
API reference for `lipids` submodule
====================================

.. automodule:: buildh.lipids
    :members:
API reference for `UI` submodule
==================================

.. automodule:: buildh.UI
    :members:
