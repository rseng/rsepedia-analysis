---
title: 'ESVM: an open-source finite volume Electrostatic Vlasov-Maxwell code'
tags:
  - Fortran
  - OpenMP
  - Python
  - Electrostatic 
  - Collisionless
  - Plasma
  - Poisson 
  - Maxwell-Gauss
  - Maxwell-Ampere
  - Vlasov
  - Advection
  - Finite volume
  - Donor-cell
  - Lax-Wendroff
  - Beam-Warming
  - Fromm 
  - minmod
  - superbee
  - Van Leer
  - MUSCL
  - Landau damping
  - Two-stream instability
  - Electrostatic wakefield
authors:
  - name: Michaël J. Touati
    orcid: 0000-0001-7590-0941
    affiliation: "1, 2, 3"
affiliations:
 - name: Department of Electrical Engineering, University of California Los Angeles, Los Angeles, CA 90095, USA
   index: 1
 - name: Group of Lasers and Plasmas, IPFN, IST, Universidade de Lisboa, Lisbon, Portugal
   index: 2
 - name: Centro de Láseres Pulsados de Salamanca (CLPU), Edificio M5, Parque Cientfico, C/ Adaja 8, 37185 Villamayor, Salamanca, Spain (current affiliation)
   index: 3
date: 5 August 2021
bibliography: paper.bib
---

# Summary

A plasma is a set of charged particles consisting of electrons and ionized atoms whose quantity is sufficiently large to behave collectively through the long-distance electromagnetic fields they produce. It is thought that more than 99.9% of visible matter in the universe is in the plasma state. In a collisionless plasma  consisting in an ionized gas composed of electrons moving inbetween much heavier ions, any electrostatic field is rapidly screened by the plasma electrons over the Debye screening distance [@DebyeHuckel:1923]. When the number of electrons in these Debye spheres can be assumed to be infinite, the plasma electron population is correctly described by the Vlasov equation [@Vlasov:1938] that neglects all correlations between particles such as the binary Coulomb collisions between them. In addition to being simple, the resulting Vlasov-Maxwell set of equations is extremely rich in physics and has many applications ranging from astrophysics and theoretical plasma physics to intense laser-matter interaction experiments. [ESVM](https://github.com/michaeltouati/ESVM) (ElectroStatic Vlasov-Maxwell) is a Vlasov-Maxwell Fortran 95 standard-compliant code, parallelized with OpenMP and using Python 3 for post-processing, that allows for the study of these collisionless plasmas. Many finite volume advection schemes [@Godunov:1959] are implemented in order to discretize the Vlasov equation, namely:

- the donor-cell scheme, i.e., the downwind / upwind scheme [@Courant:1952] depending on the advection direction in each phase-space cell, 
- the Lax-Wendroff scheme [@LaxWendroff:1960], 
- the Fromm scheme [@Fromm:1968],
- the Beam-Warming scheme [@BeamWarming:1976],
- the Van Leer scheme [@VanLeerIII:1977],
- the minmod scheme [@Roe:1986], 
- the superbee scheme [@Roe:1986], and 
- two Monotonic Upwind-centered Scheme for Conservation Laws (MUSCL) [@VanLeerV:1977] schemes MUSCL2 [@Crouseilles:2004] and MUSCL1 [@Duclous:2009]. 

Unlike the linear second order Lax-Wendroff, Fromm, and Beam-Warming schemes, the non-linear second order minmod, superbee, Van Leer, and MUSCL schemes make use of a Total Variation Diminishing (TVD) non-linear flux limiter with the price of becoming a first order scheme in some phase-space cells to limit the numerical oscillations. The donor-cell scheme is a first order method and has the pros of limiting such eventual oscillations but the cons of being numerically less consistent and more diffusive. In ESVM, the discretized Vlasov equation is coupled with the self-consistent Maxwell-Gauss equation or equivalently with the Maxwell-Ampere equation with Maxwell-Gauss equation computed at the first time step, only. While the second order Maxwell-Gauss solver needs a computationally expensive inversion of a tridiagonal matrix for the computation of the Poisson equation, the Maxwell-Ampere equation solver makes use of a faster first order numerical scheme (in time). Both absorbing and periodic boundary conditions for both the particles and the fields are implemented. Python scripts, using the Matplotlib and Numpy packages, are provided to automatically extract and plot the stored simulation results. The simulation parameters are described in the [input deck](https://github.com/michaeltouati/ESVM/blob/master/input-deck) and they can be modified without having to recompile the code. Compilation rules can be modified in the [makefile](https://github.com/michaeltouati/ESVM/blob/master/makefile) depending on the user's compiler preferences. Classical plasma physics academic case simulations that need less than one CPU$\times$hour each, tools for testing the compilation of the code, and tools for checking the simulation parameters are provided. 

# Statement of need

[ESVM](https://github.com/michaeltouati/ESVM) has been developed in order to adapt simulations to specific plasma physics problems by chosing the more adequate finite volume numerical advection scheme in order to compute the Vlasov equation phase-space advection derivatives and to chose between computing the Maxwell-Gauss equation or the Maxwell-Ampere equation with Maxwell-Gauss equation computed at the first time step, only. The code aims at beeing used by the open-source highly parallel computing (HPC) plasma physics community, ranging from under or post-graduate students to teachers and researchers who usually use Particle-In-Cell (PIC) codes [@Dawson:1962] to study collisionless plasmas. Indeed, the PIC method may prohibit the study of plasma physical processes on large time scales and/or for very dense collisionless plasmas due to the statistical and numerical fluctuations of the computed quantities imposed by the use of a finite number of particles. Also, plasma instabilities naturally develop in PIC codes, seeded by the available fluctuations spatial spectrum k-vector for which the instability growth rate is maximum and some small amplitude plasma physical processes may be hidden under the fluctuactions level. Compared to the many open source PIC codes such as [Smilei](https://github.com/SmileiPIC/Smilei) [@Derouillat:2018] and semi-Lagrangian codes such as [vmf90](https://github.com/pdebuyl/vmf90) [@Debuyl:2014], there are no open source finite volume Vlasov-Maxwell codes in the literature that are not based on an expansion method such as [OSHUN](https://github.com/UCLA-Plasma-Simulation-Group/OSHUN) [@Tzoufras:2011],  [AMoRE](https://github.com/michaeltouati/AMoRE) [@Touati:2014], or [Vlapy](https://github.com/joglekara/VlaPy) [@Joglekar:2020]. Finally, since the Vlasov equation is a conservation equation of the probable number of particles in the phase-space, using a finite volume method in order to compute the Vlasov equation has the advantage of allowing for the use of numerical schemes that are numerically flux conserving and/or that ensure the distribution function positivity compared to other numerical methods. ESVM has already been used during courses for under and post-graduate students about the "numerical tools for laser-plasma interaction Physics" and it is currently used for theoretical Plasma Physics investigations.

# Future work

The author plans, in the near future, to:

1. provide another plasma physics academic simulation about one non-linear BGK electron plasma wave, from the name of its finder I. B. Bernstein, J. M. Greene, and M. D. Kruskal [@BernsteinGreenKruskal:1957]
2. provide another plasma physics academic simulation about the echo of two plasma electron waves [@Gould:1967]
3. implement non-equally spaced phase-space cells
4. implement high order Weighted Essentially Non-Oscillatory (WENO) advection schemes [@Liu:1994]
5. compute the plasma ion Vlasov equation to allow for the ions to be mobile 
6. store the simulation results in hdf5 files instead of text files
7. implement MPI parallelization
8. implement vectorization
9. extend the code to the relativistic regime: ESVM $\Rightarrow$ RESVM for open source Relativistic ElectroStatic Vlasov-Maxwell code
10. implement a BGK collision operator, from the name of its finder P. L. Bhatnagar, E. P. Gross,  and M. Krook [@BhatnagarGrossKrook:1954]
11. extend the code to 1D-2V and 1D-3V phase-space electrostatic plasma simulations
12. implement the Landau [@Landau:1936] and Belaiev-Budker [@BelaievBudker:1957] relativistic collision operators using the Rosenbluth potentials [@Rosenbluth:1957] and their relativistic Braams-Karney extension [@BraamsKarney:1987]: (R)ESVM $\Rightarrow$ (R)ESVFPM for open source (Relativistic) ElectroStatic Vlasov-Fokker-Planck-Maxwell code
13. extend the code to electrostatic 2D-1V, 2D-2V, and 2D-3V phase-space plasma simulations: (R)ESV(FP)M $\Rightarrow$ (R)ESV(FP)M2 for open source (Relativistic) ElectroStatic Vlasov-(Fokker-Planck-)Maxwell in 2D
14. extend the code with the second order finite difference Yee scheme [@Yee:1966] to electromagnetic 2D-1V, 2D-2V, and 2D-3V phase-space plasma simulations: (R)ESV(FP)M(2) $\Rightarrow$ (R)EMV(FP)M(2) for open source (Relativistic) ElectroMagnetic Vlasov-(Fokker-Planck)-Maxwell (in 2D)
15. implement the Perfectly Matched Layer (PML) technique [@Berenger:1994] to absorb the electromagnetic fields at the spatial simulation box boundaries
16. deploy the code to GPU architectures.

# References
[<img src='https://img.shields.io/badge/Fortran-%23734F96.svg?style=for-the-badge&logo=fortran&logoColor=white' height="20">](https://fortran-lang.org/)
[<img src='https://img.shields.io/badge/GNU%20Make-black?style=for-the-badge&logo=gnu&logoColor=#7D929E' height="20">](https://www.gnu.org/software/make/)
[<img src='https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white' height="20">](https://www.gnu.org/software/bash/)
[<img src='https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54' height="20">](https://www.python.org/)
[<img src='https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white' height="20">](https://numpy.org/)
[<img src='https://matplotlib.org/_static/logo2_compressed.svg' height="20">](https://matplotlib.org/stable/index.html#)
[<img src='https://img.shields.io/badge/latex-%23008080.svg?style=for-the-badge&logo=latex&logoColor=white' height="20">](https://www.latex-project.org//)
[![JOSS : status](https://joss.theoj.org/papers/d0f6a6298710d69b63a8a62d843f8f88/status.svg)](https://joss.theoj.org/papers/d0f6a6298710d69b63a8a62d843f8f88)
[![License : GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![Clones : 14 days](https://img.shields.io/badge/dynamic/json?color=success&label=Clones%20(<15%20days)&query=count&url=https://github.com/michaeltouati/ESVM/blob/master/.github/clone.json?raw=True&logo=github)](https://github.com/michaeltouati/ESVM/actions/workflows/clones.yml)
[![Compilation check](https://github.com/michaeltouati/ESVM/actions/workflows/compilation.yml/badge.svg?branch=master)](https://github.com/michaeltouati/ESVM/actions/workflows/compilation.yml)
[![Tests check](https://github.com/michaeltouati/ESVM/actions/workflows/tests.yml/badge.svg?branch=master)](https://github.com/michaeltouati/ESVM/actions/workflows/tests.yml)
[![Plotting tools check](https://github.com/michaeltouati/ESVM/actions/workflows/plots.yml/badge.svg?branch=master)](https://github.com/michaeltouati/ESVM/actions/workflows/plots.yml)
[![Open issues](https://img.shields.io/github/issues/michaeltouati/ESVM)](https://github.com/michaeltouati/ESVM/issues)
[![Closed issues](https://img.shields.io/github/issues-closed/michaeltouati/ESVM)](https://github.com/michaeltouati/ESVM/issues)
[![Open pull requests](https://img.shields.io/github/issues-pr/michaeltouati/ESVM)](https://github.com/michaeltouati/ESVM/pulls)
[![Closed pull requests](https://img.shields.io/github/issues-pr-closed/michaeltouati/ESVM)](https://github.com/michaeltouati/ESVM/pulls)
<!-- [![Views : 14 days](https://img.shields.io/badge/dynamic/json?color=success&label=Views%20(<15%20days)&query=count&url=https://github.com/michaeltouati/ESVM/blob/master/.github/view.json?raw=True&logo=github)](https://github.com/michaeltouati/ESVM/actions/workflows/views.yml) -->
<!-- [![Downloads](https://img.shields.io/github/downloads/michaeltouati/ESVM/total)](https://github.com/michaeltouati/ESVM/releases) -->

## Written by Michaël J TOUATI

[ESVM](https://github.com/michaeltouati/ESVM) (ElectroStatic Vlasov-Maxwell) is a Vlasov-Maxwell [Fortran 95](https://fortran-lang.org/) standard-compliant code, parallelized with OpenMP and using Python 3 for post-processing, that allows for the study of collisionless plasmas. Vlasov equation is coupled with the self-consistent Maxwell-Gauss equation, or equivalently with the Maxwell-Ampere equation with Maxwell-Gauss equation computed at the first time step, only. Both absorbing and periodic boundary conditions for both the particles and the fields are implemented. More pieces of information can be found in the [esvm.pdf](https://github.com/michaeltouati/ESVM/blob/master/esvm.pdf) peer-reviewed article draft. [Python](https://www.python.org/) scripts, using the [Matplotlib](https://matplotlib.org/stable/index.html#) and [Numpy](https://numpy.org/) packages, are provided to automatically extract and plot the stored simulation results. The simulation parameters are described in the [input-deck](https://github.com/michaeltouati/ESVM/blob/master/input-deck) and they can be modified without having to recompile the code. Compilation rules can be modified in the [makefile](https://github.com/michaeltouati/ESVM/blob/master/makefile) depending on the user compiler preferences. Classical Plasma Physics academic case simulations that need less than one CPUxhour each, tools for testing the compilation of the code and tools for checking the simulation parameters are provided.

# Simulation plot examples

<p align="center">
  <img width="400" height="325" src="test-cases/Wakefield-Emission/Ex.png">
  <img width="400" height="325" src="test-cases/Linear-Landau-Damping/energy_log.png">
</p>

<p align="center">
  <img width="400" height="325" src="test-cases/Two-Stream-Instability/fe_81.png">
  <img width="400" height="325" src="test-cases/Non-Linear-Landau-Damping/logfe_124.png">
</p>

# Code units

The code units consist in the commonly used electrostatic units : the electron mass for masses, the elementary charge for electrical charges, the inverse of the Langmuir plasma electron angular frequency for times, the Debye electron screening length for spaces and the average plasma electron density for spatial densities. The initial plasma electron velocity distribution standard deviation is therefore an important unit parameter of normalization since it fixes indirectly the spatial unit.

# Compiling the code

Modify the [makefile](https://github.com/michaeltouati/ESVM/blob/master/makefile) as a function of the wished compilation options and the Fortran compiler installed on your computer and then type

```sh
make
```
The compilation can be tested by typing
```sh
make test
```
The tests consist in comparing file1 and file2 where :
* file1 is one test simulation terminal output performed with an input deck located in the directory 'test-cases/Tests/' and
* file2 is the terminal output of the corresponding simulation already performed by the developper also located in 'test-cases/Tests/'.

# Running a simulation

Fill the wished [input-deck](https://github.com/michaeltouati/ESVM/blob/master/input-deck) (all parameters are described inside), eventually check them by typing
```sh
./check-input-deck
```
or
```sh
make check
```
and then type
```sh
./esvm
```
or
```sh
make run
```
to run the simulation.

# Plotting the simulation results

All simulation results are stored in files located in the directory 'results'. 
Python scripts allowing to extract and plot the simulation results are located in the directory 'sources/plot'.
They can be used by simply typing :
```sh
make plot
```
to plot all the results, even when the simulation is still running. The resulting plots will be stored in the directory 'figures'. It can also be plotted separately :
- the energies scalar plots by typing :
```sh
make plot_energies
```
or
```sh
python3 sources/plot/plot_logfe.py
```
- the 1D plasma electron hydrodynamic moments space-time density maps by typing :
```sh
make plot_hydro2D
```
or
```sh
python3 sources/plot/plot_hydro2d.py
```
- the 1D plasma electron hydrodynamic moments scalar plots by typing : 
```sh
make plot_hydro1D
```
or
```sh
python3 sources/plot/plot_hydro1d.py
```
- or the 1D1V plasma electron distribution function phase-space density maps by typing :
```sh
make plot_fe
```
or
```sh
python3 sources/plot/plot_fe.py
```
If you need to plot the 1D1V plasma electron distribution function phase-space density maps in logarithmic scale instead, type :
```sh
make plot_logfe
```
or
```sh
python3 sources/plot/plot_logfe.py
```

# Cleaning the directory

If you want to remove from the [ESVM](https://github.com/michaeltouati/ESVM) directory :
- the compilation files and executables, type :
```sh
make distclean
```
- the directory 'figures' containing all simulations results plots, type :
```sh
make figclean
```
- the directory 'results' containing all simulations results data files, type :
```sh
make resclean
```
- the three previous ones, type :
```sh
make clean
```
Be careful, the three latters will remove all simulations results and or figures. Store them elsewhere if you don't want to lose them.

# License
[ESVM](https://github.com/michaeltouati/ESVM) is distributed under the terms of the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) license. 
## Written by Michaël J TOUATI

“Physics is like sex: sure, it may give some practical results, but that's not why we do it.”-Richard P. Feynman

“We need science education to produce scientists, but we need it equally to create literacy in the public. Man has a fundamental urge to comprehend the world about him, and science gives today the only world picture which we can consider as valid. It gives an understanding of the inside of the atom and of the whole universe, or the peculiar properties of the chemical substances and of the manner in which genes duplicate in biology. An educated layman can, of course, not contribute to science, but can enjoy and participate in many scientific discoveries which as constantly made. Such participation was quite common in the 19th century, but has unhappily declined. Literacy in science will enrich a person's life.”-Hans Bethe

# Contributing to ESVM

Thank you very much for your interest in contributing to [ESVM](https://github.com/michaeltouati/ESVM). Sharing Richard P. Feynmann enthusiasm to study Mother Nature Physics laws (and to play drum too 😄) while also deeply convinced by this Hans Bethe statement, it seemed to me natural to develop [ESVM](https://github.com/michaeltouati/ESVM) with Github in order to try to understand the complexity of plasmas, numerical analysis and High Parallel Computing (HPC) and to make accessible this small piece of Modern theoretical and computational Physics knowledge to the largest Public as possible at the same time.
In the following can be found how to :
- report a bug you faced when compiling or using [ESVM](https://github.com/michaeltouati/ESVM)
- propose new features that would make [ESVM](https://github.com/michaeltouati/ESVM) better
- fix a bug or submit a new feature

# Report a bug

For problems related with plotting, please make sure first LaTeX, dvipng and ghostscript softwares are each working and on your PATH for the Matplotlib Python package to be able to render tex fonts. If the bug persists or if it is related to another problem, follow these steps :
1) Go to ['Issues'](https://github.com/michaeltouati/ESVM/issues) on the [ESVM master branch](https://github.com/michaeltouati/ESVM) 
2) Click on 'New issue' and then Get started (Bug report)
4) Describe the bug the more "clear and concise" as possible in the title starting with "Bug :"
5) Describe with details the bug providing the more pieces of information as possible such as input-deck, terminal output and/or screenshot, etc... but more importantly, your environment :
- OS: [e.g., Ubuntu 20.04]
- Fortran compiler [e.g., gfortran 11.2.0]
- Python version [e.g., Python 3.7.11]
- Matplotlib Python package version [e.g., Matplotlib 3.4.3]
- Numpy Python package version [e.g., Numpy 1.21.2]
6) Click on 'Submit new issue'

# Propose a new feature

If you want to propose a new feature in [ESVM](https://github.com/michaeltouati/ESVM) that is not already in the perspectives of the code mentioned in the last section of the documentation [esvm.pdf](https://github.com/michaeltouati/ESVM/blob/master/esvm.pdf), follow these steps :
1) Go to ['Issues'](https://github.com/michaeltouati/ESVM/issues) on the [ESVM master branch](https://github.com/michaeltouati/ESVM)
2) Click on 'New issue' and then Get started (Feature request)
4) Describe the new feature request the more "clear and concise" as possible in the title starting with "Feature request :"
5) Describe with details the requested feature and
6) Click on 'Submit new issue'

It will be a pleasure to discuss about it.

# Fix a bug or submit a new feature

[ESVM](https://github.com/michaeltouati/ESVM) uses the [Fork and pull model](https://docs.github.com/en/github/collaborating-with-pull-requests/getting-started/about-collaborative-development-models). In order to fix a bug or to submit a new feature that you've added in [ESVM](https://github.com/michaeltouati/ESVM), follow these steps:

1) Fork the repo and create your own branch from the [ESVM 'master' branch](https://github.com/michaeltouati/ESVM).
2) Add your code. Please, keep the code Fortran 95 standard compliant with same coding style by : 
- respecting 2 spaces for indentation rather than tabs
- not using object oriented Fortran 2003 features
- not forgetting to deallocate allocated arrays
- not using allocatable arrays with assumed-shape dummy argument in subroutines or functions for bounds-preservation
- etc...
3) If you've added a new feature that needs new simulation parameters in the [input-deck](https://github.com/michaeltouati/ESVM/blob/master/input-deck) and updated correspondingly the source file [input.f90](https://github.com/michaeltouati/ESVM/blob/master/sources/input.f90), add their descriptions in the [input-deck](https://github.com/michaeltouati/ESVM/blob/master/input-deck) following the same style
```sh
##                                                                   ##
## T  = electron temperature in eV                                   ##
##                                                                   ##
#######################################################################
#
#T  1000.
#
```
4) If you've added a new feature, let's call it 'my-new-feature', create a directory 'my-new-feature' in 'test-cases/Tests/New-features' for testing and add a typical input deck inside that uses your new feature with the simulation name `#simu New-features` and parameters such as :
- `#x_min 0.`, `#x_max 5.`, `#d_x 0.25`, `#vx_min -5.`, `#vx_max 5.`, `#d_vx 0.1` and
- `#cfl 9.e-1`, `#L_t 5.`, `#dt_diag 0.25`
for the test to be fast and a file entitled 'output' containing the terminal output generated by the corresponding simulation. For this, you can just type from the directory that you've created 'test-cases/Tests/New-features/my-new-feature'
```sh
../../../.././esvm > output
```
and remove the generated directory 'results'. Finally, add in the makefile at the section ```## TESTING ##``` (after the line ```TEST_DIR=' Maxwell/Ampere'```) the following line
```sh
TEST_DIR += New-features/my-new-feature
```
5) Please, mention your contribution following the same style at the beginning of :
- the modified Fortran source file
```
!! Initial commit written by Michaël J TOUATI - Dec. 2015
!! my-new-feature or my-bug-fix commit by my-name - my-add-date
```
- the modified Python script
```
# Initial commit written by Michaël J Touati - Dec. 2015
# my-new-feature or my-bug-fix commit by my-name - my-add-date
```
- the input-deck if 3.
```sh
##   Initial commit written by Michaël J TOUATI - Dec. 2015          ##
##   my-new-feature commit by my-name - my-add-date                  ##
```
- the makefile if 4.
```sh
## Initial commit written by Michaël J TOUATI - Dec. 2015
## my-new-feature commit by my-name - my-add-date
```
6) Make sure your code passes all the tests by typing on your terminal:
```sh
make test
```
7) Ensure the code compiles with the following compiler debugging options:
- `-g -traceback -fopenmp -r8 -std95 -fpe0 -debug all -debug-parameters all -C` for the INTEL compiler ifort
- `-fdefault-real-8 -O -g -fopenmp -Wall -fcheck=all -fbacktrace -std=f95 -fall-intrinsics -ffpe-trap=invalid,zero,overflow` for the GNU compiler gfortran
8) Save the simulation results and/or figures that you want to keep somewhere else, clean the repo by typing
```sh
make clean
```
and issue that pull request! :relieved:

# License
When you submit code changes, your submissions are understood to be under the same [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) that covers [ESVM](https://github.com/michaeltouati/ESVM). 
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
michaelj.touati@gmail.com.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
---
name: Bug report
about: Create a report to help us improve
title: 'Bug :'
labels: ''
assignees: michaeltouati

---

**Describe the bug**
A clear and concise description of what the bug is.

**input-deck**
to reproduce the behavior

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g., Ubuntu 20.04]
- Fortran compiler [e.g., gfortran 11.2.0]
- Python version [e.g., Python 3.7.11]
- Matplotlib Python package version [e.g., Matplotlib 3.4.3]
- Numpy Python package version [e.g., Numpy 1.21.2]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: 'Feature request :'
labels: ''
assignees: michaeltouati

---

**Is your feature request related to a particular Plasma Physics problem, additional Physics, code improvement or a new numerical scheme? Please describe.**
A clear and concise description of what the problem is. Ex. 

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Additional context**
Add any other context or screenshots about the feature request here.
