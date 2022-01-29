[![Documentation Status](https://readthedocs.org/projects/turtlefsi2/badge/?version=latest)](https://turtlefsi2.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/KVSlab/turtleFSI.svg?branch=master)](https://travis-ci.org/KVSlab/turtleFSI)
[![status](https://joss.theoj.org/papers/b7febdaa2709205d40b51227091c3b0b/status.svg)](https://joss.theoj.org/papers/b7febdaa2709205d40b51227091c3b0b)

# turtleFSI - a Fluid-Structure Interaction Solver

<p align="center">
    <img src="figs/turtleFSI_swim.gif" width="288" height="200" alt="turtleFSI_swim"/>
    <img src="figs/turek_benchmark.gif" width="404" height="200" alt="turtleFSI_swim"/>
</p>
<p align="center">
  To the left we show a turtle swimming (in turtleFSI), and to the right, the classical Turek benchmark (FSI2).
</p>


Description
-----------
turtleFSI is a monolithic fluid-structure interaction solver written in FEniCS, and has out-of-the-box high performance capabilities. The goal of turtleFSI is to provide research groups, and other individuals, with a simple, but robust solver to investigate fluid structure interaction problems.


Authors
-------
turtleFSI is developed by:

  * Andreas Slyngstad
  * Sebastian Gjertsen
  * Aslak W. Bergersen
  * Alban Souche
  * Kristian Valen-Sendstad


Licence
-------
turtleFSI is licensed under the GNU GPL, version 3 or (at your option) any
later version. turtleFSI is Copyright (2016-2019) by the authors.


Documentation
-------------
For an introduction to turtleFSI, and tutorials, please refer to the [documentation](https://turtlefsi2.readthedocs.io/en/latest/).

If you wish to use turtleFSI for journal publications, please refer to the [JOSS publication](https://joss.theoj.org/papers/10.21105/joss.02089#):

Bergersen et al., (2020). turtleFSI: A Robust and Monolithic FEniCS-based Fluid-Structure Interaction Solver. Journal of Open Source Software, 5(50), 2089, https://doi.org/10.21105/joss.02089

Installation
------------
turtleFSI is build upon the open source Finite Elements FEniCS project (version 2018.1.0 or 2019.1.0).
Please refer to the respective FEniCS documentation for installing the dependencies on your system.

However, if you are using Linux or MaxOSX you can install turtleFSI through anaconda::

        conda create -n your_environment -c conda-forge turtleFSI

You can then activate your environment by runing ``source activate your_environment``.
You are now all set, and can start running fluid-structure interaction simulations. If 
you would like to use ``save_deg > 1``, then please install fenicstools by runningthe following::

	pip install git+https://github.com/mikaem/fenicstools

Use
---
Run turtleFSI with all the default parameters::
   ``turtleFSI``

See all the command line parameters run the following command::
  ``turtleFSI -h``

Run a specific problem file::
  ``turtleFSI --problem [path_to_problem]``

When calling a specific problem file, turtleFSI will first look for the file name locally, then check if the file name is present in the directory "/turtleFSI/problems/".
Please refere to the [documentation](https://turtlefsi2.readthedocs.io/en/latest/) to learn how to define a new problem file and for a more complete description of usage.


Contact
-------
The latest version of this software can be obtained from

  https://github.com/KVSlab/turtleFSI

Please report bugs and other issues through the issue tracker at:

  https://github.com/KVSlab/turtleFSI/issues
# Contributing to turtleFSI

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change.

### Testing

Please provide unit tests for the new code you create, testing the main functionality or feature to be submitted. We can always use more test coverage!

### Submitting changes

In order to submit you changes, please send a [GitHub Pull Request to turtleFSI](https://github.com/KVSlab/turtleFSI/pull/new/master) with a clear list of what you've done (read more about [pull requests](http://help.github.com/pull-requests/)). Please follow our coding conventions (below) and make sure all of your commits are atomic (one feature per commit).

Always write a clear commit message when submitting your changes. One-line messages are fine for small changes, but bigger changes should look like this:

    > $ git commit -m "A brief summary of the commit
    >
    > A paragraph describing what changed and its impact."

### Coding conventions

#### Formatting
* Avoid inline comments.
* Break long lines after 120 characters.
* Delete trailing whitespace.
* Don't include spaces after `(`, `[` or before `]`, `)`.
* Don't misspell.
* Use 4 space indentation.
* Use an empty line between methods.
* Use spaces around operators, except for unary operators, such as `!`.
* Use spaces after commas, after colons and semicolons, around `{` and before
  `}`.

#### Naming

* Avoid abbreviations.
* Use snake case for variables and methods.
* Use camel case for classes.
* Name variables, methods, and classes to reveal intent.

#### Organization

* Order methods so that caller methods are earlier in the file than the methods
  they call.
* Place methods receiving command line arguments at the bottom of the file, but above the top-level script environment check.
* Separate local and global imports of modules.

### Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

### Our Standards

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

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at alban@simula.no. All complaints will be reviewed and
investigated and will result in a response that is deemed necessary and
appropriate to the circumstances. The project team is obligated to maintain
confidentiality with regard to the reporter of an incident. Further details of
specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
---
title: 'turtleFSI: A Robust and Monolithic FEniCS-based Fluid-Structure Interaction Solver'
tags:
  - fluid-structure interaction
  - FEniCS
  - finite-elements
  - numerical methods
authors:
  - name: Aslak W. Bergersen
    orcid: 0000-0001-5063-3680
    affiliation: 1
  - name: Andreas Slyngstad
    affiliation: 1
  - name: Sebastian Gjertsen
    affiliation: 1
  - name: Alban Souche
    orcid: 0000-0001-7547-7979
    affiliation: 1
  - name: Kristian Valen-Sendstad
    orcid: 0000-0002-2907-0171
    affiliation: 1
affiliations:
  - name: Department of Computational Physiology, Simula Research Laboratory, Fornebu, Norway
    index: 1
date: 15 June 2020
bibliography: paper.bib
---

# Summary

It is often sufficient to study fluids [@Moin:1998] and solids [@Holzapfel:2002] in isolation to gain fundamental insights into a physical problem, as other factors may play a secondary role and can be neglected.  On the other hand, there are certain phenomena or situations where the stresses on or by a fluid or a solid can lead to large deformations, and the interaction between fluids and solids are essential [@LeTallec:2001]. Computational fluid-structure interaction (FSI) is an active field of research with much focus on numerical accuracy, stability,  and convergence rates. At the same time, there is also a sweet spot in between these areas of research where there is a need to experiment with FSI without having an in-depth, bottom-up mathematical understanding of the problem, but where a physical insight might suffice. Therefore, the aim was to develop a fully monolithic and robust entry-level research code with ease-of-use targeted towards students, educators, and researchers.

FEniCS [@Logg:2012] has emerged as one of the leading platforms for development of scientific software due to the close connection between mathematical notation and compact computer implementation, where highly efficient C++ code is compiled during execution of a program. Combined with the out-of-the-box entry-level high-performance computing capabilities, FEniCS was a natural choice of computing environment. Compared to other open-source FSI solvers [@Malinen:2013; @Heil:2006; @Jasak:2007], turtleFSI is written in only a couple of hundred lines of high-level Python code, in contrast to tens of thousands of lines of low-level C++ code. This provides full transparency and a unique opportunity for researchers and educators to modify and experiment with the code, while still providing out of the box entry-level high-performance computing capabilities. Furthermore, because of the close resemblance between mathematics and code in FEniCS, users can make additions or modifications with ease.

The turtleFSI solver relies on a fully monolithic approach in the classical arbitrary Lagrangian-Eulerian formulation, and we used the generalized theta scheme for temporal discretization and P2P1P2 elements for velocity, pressure, and displacement, respectively. We implemented and evaluated four different mesh lifting operators, ranging from a simple and efficient second-order Laplace equation, most suitable for small deformations, to more sophisticated and computationally expensive 4th order bi-harmonic equations that can handle larger mesh deformations. We used The Method of Manufactured Solutions to verify the implementation. The obtained results are formally second-order accurate (L2) in space and time [@Wick:2011], respectively, and we demonstrate that all building blocks of code exhibit desired properties. The solver's validity was confirmed using the classical Turek Flag benchmark case [@Turek:2006] with a good agreement – including a diverged numerical solution for long term evolution under certain conditions, as expected. For a complete justification of computational approaches and further details, we refer to [@Slyngstad:2017; @Gjertsen:2017]. We demonstrate adequate strong scaling up to 64 cores (from one cluster node), although the latter is problem size-dependent. In the online documentation, we provide benchmarks, tutorials, and simple demos. The naive FEniCS implementation provides full transparency with compact code, which can easily be adapted to other 2D or 3D FSI problems.

In conclusion, turtleFSI is not a superior FSI solver in terms of speed, but it is a robust entry-level FSI solver and performs exactly as designed and intended; ‘slow and steady wins the race’.


# turtleFSI in Action

turtleFSI comes with several problem files, found under /turtleFSI/problems/, to illustrate the usage and document the Turek flag benchmarks used to validate the implementation of the solver. Here are some illustrations of the execution and outputs expected from the solver.

![Fluid_Turek*='#center'](./cfd_illu.png){ width=100% }\
**Figure 1:**
  Fluid dynamics benchmark snapshot. Simulation executed with the command:
  ```
  turtleFSI --problem TF_cfd
  ```

![Solid_Turek*='#center'](./csm_illu.png){ width=100% }\
**Figure 2:**
  Solid mechanics benchmark snapshots. Simulation executed with the command:
  ```
  turtleFSI --problem TF_csm
  ```

![FSI_Turek*='#center'](./fsi_illu.png){ width=100% }\
**Figure 3:**
  Full fluid-structure interaction benchmark snapshot. Simulation executed with the command:
  ```
  turtleFSI --problem TF_fsi
  ```

# Acknowledgements
The study was supported by The Research Council of Norway through the Center for Biomedical Computing (grant 179578), the Centre for Cardiological Innovation (grant number 203489), and the SIMMIS project (grant number 262827). Simulations were performed on the Abel Cluster (University of Oslo and the Norwegian metacenter for High Performance Computing (NOTUR), project nn9316k), and the Experimental Infrastructure for Exploration of Exascale Computing (eX3) cluster (Norwegian Research Council grant 270053).

# References
