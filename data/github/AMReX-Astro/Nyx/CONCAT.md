### Citation

If you use Nyx, we appreciate you citing the main code paper, [Almgren et al. (2013)](https://iopscience.iop.org/article/10.1088/0004-637X/765/1/39) [[PDF]](https://iopscience.iop.org/article/10.1088/0004-637X/765/1/39/pdf)


Citation in BibTeX:
```
@ARTICLE{Nyx,
         author = {{Almgren}, A.~S. and {Bell}, J.~B. and {Lijewski}, M.~J. and {Luki{\'c}}, Z. and {Van Andel}, E.},
         title = "{Nyx: A Massively Parallel AMR Code for Computational Cosmology}",
         journal = {The Astrophysical Journal},
         archivePrefix = "arXiv",
         eprint = {1301.4498},
         keywords = {gravitation, hydrodynamics, methods: numerical},
         year = 2013,
         month = mar,
         volume = 765,
         eid = {39},
         pages = {39},
         doi = {10.1088/0004-637X/765/1/39}
}
```
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5059767.svg)](https://doi.org/10.5281/zenodo.5059767)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03068/status.svg)](https://doi.org/10.21105/joss.03068)
[![AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io)


![Nyx](https://github.com/AMReX-Astro/Nyx/blob/development/Util/banner.jpeg)

*An adaptive mesh, massively-parallel, cosmological simulation code*

******

## About

Nyx code solves equations of compressible hydrodynamics on an adaptive grid
hierarchy coupled with an N-body treatment of dark matter. The gas dynamics in
Nyx uses a finite volume methodology on a set of 3-D Eulerian grids;
dark matter is represented as discrete particles moving under the influence of
gravity. Particles are evolved via a particle-mesh method, using Cloud-in-Cell
deposition/interpolation scheme. Both baryonic and dark matter contribute to
the gravitational field. In addition, Nyx includes physics needed to
accurately model the intergalactic medium: in optically thin limit and assuming
ionization equilibrium, the code calculates heating and cooling processes of the
primordial-composition gas in an ionizing ultraviolet background radiation field.
Additional physics capabilities are under development.

While Nyx can run on any Linux system in general, we particularly focus on supercomputer systems.
Nyx is parallelized with MPI + X, where X can be OpenMP on multicore architectures and
CUDA/HIP/DPC++ on hybrid CPU/GPU architectures.
In the OpenMP regime, Nyx has been successfully run at parallel concurrency
of up to 2,097,152 on NERSC's Cori-KNL. With Cuda implementation, it was run on up to
13,824 GPUs on OLCF's Summit.

More information on Nyx can be found at the [main web page](http://amrex-astro.github.io/Nyx/) and
the [online documentation](https://amrex-astro.github.io/Nyx/docs_html/).

## Standards and dependencies

To compile the code we require C++11 compliant compilers that support MPI-2 or
higher implementation.  If threads or accelerators are used, we require 
OpenMP 4.5 or higher, Cuda 9 or higher, or HIP-Clang.
To use Nyx, you also need [AMReX](https://github.com/AMReX-codes/amrex).

For example, to compile the Lyman alpha (LyA) executable on Summit:
```sh
$ module load gcc/6.4.0 cuda/11.0.3

$ git clone https://github.com/AMReX-Codes/amrex.git
$ git clone https://github.com/AMReX-astro/Nyx.git

$ cd Nyx/Exec/LyA
$ make -j 12 USE_CUDA=TRUE
```

See the [Getting Started section](https://amrex-astro.github.io/Nyx/docs_html/NyxGettingStarted.html) for more information.

## Development model

Please see CONTRIBUTING.md for details on how to contribute to AMReX.

## Outputs

Nyx outputs certain global diagnostics at each timestep and plot files at regular
intervals, or at user-specified redshifts. Visualization packages
[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit),
[Paraview](https://www.paraview.org/),
[yt](http://yt-project.org/),
and [Amrvis](https://github.com/AMReX-Codes/amrvis)
have built-in support for the AMReX file format used by Nyx.

In addition, Nyx interfaces with two post-processing suites, Reeber and Gimlet. Reeber
uses topological methods to construct merge trees of scalar fields, which is in
turn used to find halos.  Gimlet computes a variety of quantities
related to the Lyman-alpha forest science.  These suites are fully MPI-parallel and can
be run either *in situ* or in-transit, or with a combination of both
(see [Friesen et al. 2016](https://comp-astrophys-cosmol.springeropen.com/articles/10.1186/s40668-016-0017-2)).

## License

Nyx Copyright (c) 2017, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy).  All rights
reserved.

Details of the license can be found in [license.txt](license.txt) file.

If you have questions about your rights to use or distribute this software, 
please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

## Contact

For questions, comments, suggestions, contact Jean Sexton (JMSexton@lbl.gov)
or Zarija Lukic (zarija@lbl.gov).
## Development Model

Development generally follows the following ideas:

  * New features are merged into to the `development` branch using
    Pull Requests (PRs).

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

  * Bug fixes, questions and contributions of new features are welcome!

       * Bugs should be reported through GitHub issues
       * We suggest asking questions through GitHub issues as well
       * *Any contributions of new features that have the potential
         to change answers should be done via pull requests.*
         A pull request should be generated from your fork of
         Nyx and target the `development` branch. See below for
         details on how this process works.

         In general we squash commits upon merge to have a clean history.
         *Please ensure that your PR title and first post are descriptive,
         since these will be used for a squashed commit message.*

         Please note the following:
            If you choose to make contributions to the code 
            then you hereby grant a non-exclusive, royalty-free perpetual license 
            to install, use, modify, prepare derivative works, 
            incorporate into other computer software,
            distribute, and sublicense such enhancements or derivative works
            thereof, in binary and source code form.

  * On the first workday of each month (or soon thereafer), we make a tagged release.  

## Git workflow

Nyx uses [git](https://git-scm.com) for version control. If you
are new to git, you can follow one of these tutorials:
- [Learn git with bitbucket](https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud)
- [git - the simple guide](http://rogerdudler.github.io/git-guide/)

### Make your own fork and create a branch on it

The basic workflow is:
- Fork the main repo (or update it if you already created it).
- Implement your changes and push them on a new branch `<branch_name>` on
your fork.
- Create a Pull Request from branch `<branch_name>` on your fork to branch
`development` on the main Nyx repository.

First, let us setup your local git repo. Make your own fork of the main
(`upstream`) repository:
on the [Nyx Github page](https://github.com/AMReX-Astro/Nyx), press the
fork button. 

Then, you can execute:
```
# Clone your fork on your local computer. You can get this address on your fork's Github page.
git clone --branch development https://github.com/<myGithubUsername>/Nyx.git

# Navigate into your repo, add a new remote for the main Nyx repo, and fetch it
cd Nyx
git remote add upstream https://github.com/AMReX-Astro/Nyx
git remote set-url --push upstream https://github.com/<myGithubUsername>/Nyx.git
git fetch upstream

# We recommend setting your development branch to track the upstream one instead of your fork:
git branch -u upstream/development
```
Now you are free to play with your fork (for additional information, you can visit the
[Github fork help page](https://help.github.com/en/articles/fork-a-repo)).

> Note: you do not have to re-do the setup above every time.
> Instead, in the future, you need to update the `development` branch
> on your fork with
> ```
> git checkout development
> git pull
> ```

Make sure you are on the `development` branch with
```
git checkout development
```
in the Nyx directory.

Create a branch `<branch_name>` (the branch name should reflect the piece
of code you want to add, like `new_kind_of_particle`) with
```
git checkout -b <branch_name>
```
and do the coding you want.
Add the files you work on to the git staging area with
```
git add <file_I_created> <and_file_I_modified>
```
### Commit & push your changes

Periodically commit your changes with
```
git commit -m "This is a 50-char description to explain my work"
```

The commit message (between quotation marks) is super important in order to
follow the developments and identify bugs.

For the moment, commits are on your local repo only. You can push them to
your fork with
```
git push -u origin <branch_name>
```

If you want to synchronize your branch with the `development` branch (this is useful
when `development` is being modified while you are working on
`<branch_name>`), you can use
```
git pull upstream development
```
and fix any conflicts that may occur.

Do not merge your branch for PR into your local `development` branch,
because it will make your local `development` branch diverge from the
matching branch in the main repository after your PR is merged.

### Submit a Pull Request

A Pull Request is the way to efficiently visualize the changes you made
and to propose your new feature/improvement/fix to the Nyx project.
Right after you push changes, a banner should appear on the Github page of
your fork, with your `<branch_name>`.
- Click on the `compare & pull request` button to prepare your PR.
- It is time to communicate your changes: write a title and a description for
your PR. People who review your PR are happy to know
  * what feature/fix you propose, and why
  * how you made it (created a new class than inherits from...)
  * and anything relevant to your PR (performance tests, images, *etc.*)
- Press `Create pull request`. Now you can navigate through your PR, which
highlights the changes you made.

Please DO NOT write large Pull Requests, as they are very difficult and
time-consuming to review. As much as possible, split them into small,
targeted PRs.
For example, if find typos in the documentation open a pull request that only fixes typos.
If you want to fix a bug, make a small pull request that only fixes a bug.
If you want to implement a large feature, write helper functionality first, test it and submit those as a first pull request.
If you want to implement a feature and are not too sure how to split it, just open an issue about your plans and ping other AMReX developers on it to chime in.

Even before your work is ready to merge, it can be convenient to create a PR
(so you can use Github tools to visualize your changes). In this case, please
make a "draft" PR using the drop-down menu next to the "Create pull request" button.

Once your pull request is made, we will review and potentially merge it.
We recommend always creating a new branch for each pull request, as per the above instructions.
Once your pull request is merged, you can delete your local PR branch with
```
git branch -D <branch_name>
```

and you can delete the remote one on your fork with
```
git push origin --delete <branch_name>
```

Generally speaking, you want to follow the following rules.

  * Do not merge your branch for PR into your local `development` branch that tracks Nyx
    `development` branch.  Otherwise your local `development` branch will diverge from Nyx
    `development` branch.

  * Do not commit in your `development` branch that tracks Nyx `development` branch.

  * Always create a new branch based off `development` branch for each pull request, unless you are
    going to use git to fix it later.

If you have accidentally committed in `development` branch, you can fix it as follows,
```
git checkout -b new_branch
git checkout development
git reset HEAD~2  # Here 2 is the number of commits you have accidentally committed in development
git checkout .
```
After this, the local `development` should be in sync with Nyx `development` and your recent
commits have been saved in `new_branch` branch.

If for some reason your PR branch has diverged from Nyx, you can try to fix it as follows.  Before
you try it, you should back up your code in case things might go wrong.
```
git fetch upstream   # assuming upstream is the remote name for the official amrex repo
git checkout -b xxx upstream/development  # replace xxx with whatever name you like
git branch -D development
git checkout -b development upstream/development
git checkout xxx
git merge yyy  # here yyy is your PR branch with unclean history
git rebase -i upstream/development
```
You will see something like below in your editor,
```
pick 7451d9d commit message a
pick c4c2459 commit message b
pick 6fj3g90 commit message c
```
This now requires a bit of knowledge on what those commits are, which commits have been merged,
which commits are actually new.  However, you should only see your only commits.  So it should be
easy to figure out which commits have already been merged.  Assuming the first two commits have been
merged, you can drop them by replace `pick` with `drop`,
```
drop 7451d9d commit message a
drop c4c2459 commit message b
pick 6fj3g90 commit message c
```
After saving and then exiting the editor, `git log` should show a clean history based on top of
`development` branch.  You can also do `git diff yyy..xxx` to make sure nothing new was dropped.  If
all goes well, you can submit a PR using `xxx` branch.
Don't worry, if something goes wrong during the rebase, you an always `git rebase --abort` and start over.
---
title: 'Nyx: A Massively Parallel AMR Code for Computational Cosmology'
tags:
  - C++
  - cosmology
  - hydrodynamics
  - dark matter
  - N-body
authors:
  - name: Jean Sexton
    orcid: 0000-0003-2551-1678
    affiliation: 1
  - name: Zarija Lukic
    affiliation: 2
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    affiliation: 1
  - name: Chris Daley
    affiliation: 3
  - name: Brian Friesen
    orcid: 0000-0002-1572-1631
    affiliation: 3
  - name: Andrew Myers
    orcid: 0000-0001-8427-8330
    affiliation: 1
  - name: Weiqun Zhang
    orcid: 0000-0001-8092-1974
    affiliation: 1
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory
    index: 1
  - name: Computational Cosmology Center, Lawrence Berkeley National Laboratory
    index: 2
  - name: National Energy Research Scientific Computing Center (NERSC), Berkeley, CA, USA
    index: 3
date: February 2021
bibliography: paper.bib
---

# Summary
Nyx is a highly parallel, adaptive mesh, finite-volume
N-body compressible hydrodynamics solver for cosmological simulations.
It has been used to simulate different cosmological scenarios with
a recent focus on the intergalactic medium and Lyman alpha forest.
Together, Nyx, the compressible astrophysical simulation code, Castro [@castro],
and the low Mach number code MAESTROeX [@maestroex], make up the
AMReX-Astrophysics Suite of open-source, adaptive mesh, performance-portable
astrophysical simulation codes.
Other examples of cosmological simulation research codes include Enzo [@Bryan2014], Enzo-P/Cello [@Bordner2018; @Norman2018], RAMSES [@Teyssier2002], ART [@Kravtsov1997],
FLASH [@Fryxell2000], Cholla [@Villasenor2021], as well as Gadget [@Springel2020], Gasoline [@Wadsley2017], Arepo [@Weinberger2020], Gizmo [@Hopkins2014], and SWIFT [@Schaller2016].

The core hydrodynamics solver in Nyx [@nyx-paper1] is based on the
directionally unsplit corner transport upwind method of @ctu with
piecewise parabolic reconstruction [@ppm].  In Nyx, we have
several modes of coupling the stiff heating-cooling source terms to the hydro.
The simplest method is the traditional operator splitting approach,
using Strang splitting [@strang1968] to achieve second-order accuracy in time. However,
this coupling can break down, and we have an alternative to Strang splitting
based on spectral deferred corrections (SDC), a method
that aims to prevent the hydro and stiff source terms from becoming decoupled.
The simplified SDC method uses the CTU PPM hydro together with an
iterative scheme to fully couple the source terms and hydro, still to
second-order accuracy in time [@simple_sdc].

Nyx has a set of additional physics necessary to model the intergalactic medium
using heating-cooling source terms.
The code follows the abundance of six species: neutral and ionized hydrogen,
neutral, once and twice ionized helium, and free electrons. For these species,
all relevant atomic processes - ionization, recombination, and free-free transitions are
modeled in the code. Heating and cooling source terms are calculated using
a sub-cycled approach in order to avoid running the whole code on a short, cooling timescale.
Cosmological reionization is accounted for via a spatially uniform,
but time-varying ultraviolet background (UVB) radiation field,
inputted to the code as a list of photoionization and photoheating
rates that vary with redshift. Nyx also has the capability to model flash reionization
as well as inhomogeneous reionization [@onorbe2019].

The evolution of baryonic gas is coupled with an N-body treatment of the dark
matter in an expanding universe. The mesh-based hydrodynamical baryonic gas
evolution is coupled through gravity to the particle-based representation of
dark matter. The dark matter particles are moved with a move-kick-drift algorithm
[@movekickdrift]. The Poisson equation for self-gravity of the baryonic gas and dark
matter is solved using the geometric multigrid method. Nyx simulations can optionally
model neutrino particle effects and active galactic nuclei feedback.

Nyx is built on the AMReX [@AMReX] adaptive mesh refinement (AMR)
library and is written in C++.
AMR levels are advanced at their own timestep (sub-cycling)
and jumps by factors of 2 and 4 are supported between levels.  We use
MPI to distribute AMR grids across nodes and use logical tiling with
OpenMP to divide a grid across threads for multi-core CPU machines
(exposing coarse-grained parallelism) or CUDA/HIP/DPC++ to spread the work across
GPU threads on GPU-based machines (fine-grained parallelism).  All of
the core physics can run on GPUs and have been shown to scale well.
For performance portability, we use the same source code
for both CPUs and GPUs; additionally, we implement our parallel loops
in an abstraction
layer provided by AMReX. An abstract parallel for loop accepts as arguments
a range of indices and the body of the loop to execute for a given index,
and the AMReX backend dispatches the work appropriately (e.g., one cell per
GPU thread). This strategy is similar to the way the Kokkos [@Kokkos] and
RAJA [@RAJA] abstraction models provide performance portability in C++.

# Statement of Need

While there are a number of cosmological simulation codes, Nyx
offers a few unique features.  The original motivation for developing
Nyx was to build a simulation code based on a modern,
well-supported AMR library (BoxLib which evolved to AMReX), using
unsplit integration techniques.
The large developer community contributing to AMReX
helps in Nyx continually gaining optimizations for new architectures
and various operating systems.

At present, Nyx is mostly used to model the cosmological evolution of the
intergalactic medium (IGM) - a reservoir of low density gas that fills the space
between galaxies. Different physical effects, ranging from the nature of
dark matter to different reionization scenarios related to the radiation
from star-forming galaxies, set the observable
properties of IGM, making it a powerful probe of cosmology and astrophysics.
But in order to extract scientific insights, confronting observations of
the IGM (usually through the Lyman alpha forest) against simulated models is a necessity, and that is where Nyx steps in.
Incoming observations, for example, Dark Energy Spectroscopic Instrument (DESI) or
high-resolution spectrographs like the one mounted on the Keck telescope, are noticeably
reducing the statistical errors of measurements; the community needs tools
capable of producing model universes that are of the same or better accuracy as observations.
Nyx's scalability allows modeling of
representative cosmological volumes, while maintaining the resolution needed
to resolve small fluctuations in the intergalactic gas. Nyx includes physics to allow simulations
of different cosmological and reionization scenarios, enabling users to produce
mock universes for a variety of physically relevant models.

Our main targets are high-performance computer architectures and massively parallel simulations
needed for cosmological and astrophysical research.
Given these targets, Nyx's computational optimizations focus on large simulations on HPC systems,
although smaller simulations can be run on Linux distributions and macOS using AMReX's
build system support.


# Acknowledgements

The work at LBNL was supported by the U.S. Department of Energy
under contract No. DE-AC02-05CH11231.
Nyx development was further supported by
the Exascale Computing Project (17-SC-20-SC), a collaborative effort
of the U.S. Department of Energy Office of Science and the National
Nuclear Security Administration.

# References

# Example diagnostic/analysis tool

*****

## About

This is a simple example of opening an AMReX `plt` file, reading and manipulating data in it.
It is meant to be a starting place for building analysis tools for user's specific needs.
In this particular example, given a N-dimensional data set, the program computes a one-dimensional 
profile in direction *dir* where each value in the profile is the average of the data
in the non-*dir* directions.   This can work with multiple components at one time.

## Usage

This is a stand-alone application built on top of AMReX.  Set the compiler value in the `GNUmakefile` (variable `COMP`)
and make an executable with:
```sh
$ make -j 8
```

To run, you can execute with runtime parameters in the `inputs` file
```sh
$ AmrDerive3d.Linux.gcc.gfortran.MPI.ex inputs
```
where the file `inputs` contains, for example:
```sh
infile  = plt00000
dir     = 0
nComp   = 1
sComp   = 0
```

or you can specify runtime parameters via command-line:
```sh
AmrDerive3d.Linux.gcc.gfortran.MPI.ex infile=plt00000 dir=0 nComp=2 sComp=0
```

You can also set verbose = 1 (via command line or inputs file) in order to get a quite verbose output.
Nyx: Zeldovich Test Case
========================

This Exec sets up and runs the Zeldovich problem.

Create Initial Conditions
-------------------------

Generate the initial conditions by compiling and running
``generate_zeldovich_ics.cpp``. You may change some of the parameters in
that file, such as:

 - redshifts: initial and caustic
 - number of particles
 - box length
 - number of sheets (determines perturbation wavelength)
 - perturbation offset
 - sheet normal vector (determines direction of the perturbations)
 - Hubble rate (make sure it matches your cosmology in ``probin``!)

Then compile the ICs generator with:

    $ cd ics
    $ g++ generate_zeldovich_ics.cpp -o generate_zeldovich_ics.ex
    $ ./generate_zeldovich_ics.ex

This creates files ``zeldovich_particles.ascii`` and ``zeldovich_params.csv``.
The csv file contains the parameters used to generate the particle ascii file
and is used during analysis so you don't have to worry about adjusting params in
the analysis scripts manually.

Checking Initial Conditions
---------------------------

You might want to check the file with a quick ``head zeldovich_particles.ascii``
to make sure the number of particles in the first line is fine, and the particle
format ``x y z mass xdot ydot zdot`` is good.

You can also run

    $ cd analysis
    $ python check_ics.py

which will test ``ics/zeldovich_particles.ascii`` to make sure the sheet normal
is indeed a unit vector, that the velocities perpendicular to that are 0, and
that the density in the box matches the critical (matter) density given by the
cosmology. It also plots the particle distribution in comoving space (3d
projection) and comoving phase space for you to check by eye.

We might make this automatic in the future, but this is sufficient for now.

Run Nyx
-------

Compile the Nyx executable by checking the make options in ``build/GNUMakefile``
and running ``make`` in the build directory.

Check the ``inputs`` and ``probin`` files to make sure you match the conditions
set in your ICs. Check that ``prob_hi`` matches the size of the box and ``h`` is
the same.

If you want to run on hopper, run the normal ``qsub pbs_hopper``, adjusting for
however many cores you think your run requires. If you want to run locally, use
something like ``mpirun -np 4 ./Nyx3d.Darwin.gcc.gfortran.MPI.ex inputs >&
out``.

Analyze
-------

    $ cd analysis
    $ python process_zeldovich.py

will make a ``figs`` directory, open the ``particles.ascii`` file in each
directory beginning with ``plt``, and analyze the results. For each plt ascii
particle file, it will plot the phase space with zeldovich predictions, the 3d
projection of particle positions, and a histogram of velocities perpendicular to
the normal (should be close to 0).
*****************************
Radiative Heating and Cooling
*****************************

Nyx provides the capability to compute local heating and cooling effects due to radiation.
The motivation and algorithm for the heating and cooling components is documented in :raw-latex:`\cite{lukic15}`, and the relevant code is located in the ``Source/HeatCool`` subdirectory.
The code is activated through the ``USE_HEATCOOL=TRUE`` option in the ``GNUmakefile``.
Mathematically, the heating and cooling can be described by a single ODE in each cell, to be integrated per time step :math:`\Delta t`.
This ODE exhibits a sensitive relationship to quantities such as temperature and free electron density, and consequently it often requires sophisticated integration techniques to compute correctly.

Nyx provides a few different parallelization strategies techniques for solving this ODE.

For legacy reasons, the integration scheme is selected via the ``nyx.heat_cool_type`` input parameter.
One method is to use ``USE_HEATCOOL=FALSE`` with no ODE solve (selected with ``nyx.heat_cool_type=0``).
The other method is to use ``USE_HEATCOOL=TRUE`` with a vectorized CVODE solve (selected with ``nyx.heat_cool_type=11``).

Users should note that, when compiling with GNUmake , CVODE must be compiled as a separate library; instructions for compiling CVODE are provided in Getting Started.
To link the external CVODE solver into Nyx, one must set ``USE_HEATCOOL=TRUE`` as well as ``USE_SUNDIALS=TRUE`` in the ``GNUmakefile``.

This vectorized option uses CVODE while treating groups of ODEs in different cells as a single system of coupled ODEs.
The purpose of this approach is to enable the evaluation of multiple RHSs simultaneously, for improved parallel performance.
This approach can lead to a significant performance gain in the ODE integration (which is the among the most expensive computational kernels in Nyx).
The group of cells integrated together is controlled by the tilesize used in the CVODE integration loop. If tiling is not used, the entire box of cells is integrated together.
CVODE integrates this resulting system of ODEs with using adaptive time-step control which selects
timesteps for the whole system. We set the absolute tolerance required for the ODE integration in CVODE to ``1e-4``, scaled by the initial value of each cells independent variable in the ODE and
the relative tolerance required for the ODE integration in CVODE to ``1e-4``.
These tolerances, in particular the relative tolerance, have different effects depending on whether one is integrating a single ODE at a time, or a system of ODEs simultaneously.
One should be mindful of the numerical differences which arise from these, which can be observed with the ``fcompare`` tool in AMReX.

Input flags with affect the CVODE integration include:

- ``nyx.use_sundials_constraint`` which when non-zero requires that the internal energy (while the problem is evolved) is positive
- ``nyx.use_sundials_fused`` which when non-zero uses Sundials's GPU fused operations (which are mathematically equivalent, but reduces GPU kernel launch time overhead)
- ``nyx.sundials_alloc_type`` which has up to 5 different vector memory allocation strategies and only affects executables built for GPUs
- ``nyx.use_typical_steps`` which when non-zero sets CVODE's adaptive step-size selection (which substeps the total CVODE integration time) to be ``dt / old_max_steps``. This maximum was over the entire problem domain. (In the strang case, the old maximum is taken from the same phase in the previous time-step)
- ``nyx.sundials_use_tiling`` which controls whether the MFIter loop that iterates over the CVODE integration uses tiling
- ``nyx.sundials_tile_size`` which controls the tile size of the MFIter loop that iterates over the CVODE integration
- ``nyx.hydro_tile_size`` which controls the tile size of the MFIter loop that iterates over the hydro advance that makes hydro_src for SDC integration
- ``nyx.sundials_reltol`` which controls the relative tolerance given for the CVODE integration
- ``nyx.sundials_abstol`` which controls the absolute tolerance factor which multiplies the local internal energy to set the absolute tolerance vector for the CVODE integration

A ``HeatCoolTests`` directory allows for the testing of the CVODE integration alone, given data from a full step of the ``LyA`` executable.
At present, it assumes the same tilesizes, boxarray, and distribution mapping between both runs. Note, if you supply your own filename paths
you must ensure the directory structure exists. Some additional CVODE statistics will be displayed if ``nyx.v>1``

Input flags to control these tests include:

- ``hctest_example_write`` which when non-zero writes all relevant data for 1 step, including an inputs file (should be used with ``LyA`` or similar)
- ``hctest_example_write_proc`` which when set will only write data for boxes owned by that processor number (which doesn't currently have a read mode that will take this smaller data set)
- ``hctest_example_index`` which represents the part of the filename that should typically denote step number, if unset defaults to nStep()
- ``hctest_example_read`` which when non-zero reads all relevant data for 1 step (should be used with ``HeatCoolTests``)

- ``hctest_filename_inputs`` which sets the name of the inputs file to write, if unset, hctest/inputs.${hctest_example_index}
- ``hctest_filename_badmap`` which sets the name of the BoxArray and DistributionMapping file to write, if unset, hctest/BADMAP.${hctest_example_index}
- ``hctest_filename_chunk`` which sets the name of the Chunk data file prefix to use, if unset, hctest/Chunk.${hctest_example_index} (the mfi.index() is used to index the boxes for the full filepath)
*********************
Dark matter particles
*********************

Dark matter particles are included in the simulation by setting ::

    nyx.do_dm_particles = true

in the inputs file.

When dark matter particles are present, one has the option to run with or without the baryons; 
to build with no hydro capability set::

    NO_HYDRO=TRUE

It is also possible to build with ::

    NO_HYDRO=FALSE

but turn off hydro at run-time by setting ::

  nyx.do_hydro = 0

in the inputs file.

Our default is to use double precision for the mesh data and single precision for the particles; this is set in the 
GNUMakefile by::

  PRECISION = DOUBLE
  USE_SINGLE_PRECISION_PARTICLES = TRUE

We do not recommend using single precision for the mesh data, as this is not tested and will potentially degrade the quality of the hydrodynamical solve. To use single precision for the mesh use::

  PRECISION = FLOAT

Equations
=========

If we define :math:`{\mathbf x}_i` and :math:`{\bf u}_i` as the location and velocity of particle :math:`i`, respectively, then we wish
to solve

.. math::

   \begin{aligned}
   \frac{d {\mathbf x}_i}{d t} &=& \frac{1}{a} {\mathbf u}_i \\
   \frac{d (a {\mathbf u}_i) }{d t} &=& {\mathbf g}_i\end{aligned}

where :math:`{\mathbf g}_i` is the gravitational force evaluated at the location of particle :math:`i`, i.e.,
:math:`{\mathbf g}_i = {\mathbf g}({\mathbf x}_i,t).`


Particle Time Stepping: Move-Kick-Drift Algorithm
=================================================

In each time step:

-  Solve for :math:`{\mathbf g}^n` (only if multilevel, otherwise use :math:`{\mathbf g}^{n+1}` from previous step)

-  :math:`{\mathbf u}_i^{{n+\frac{1}{2}}} = \frac{1}{a^{{n+\frac{1}{2}}}} ( (a^n {\mathbf u}^n_i) + \frac{{\Delta t}}{2} \; {\mathbf g}^n_i )`

-  :math:`{\mathbf x}_i^{n+1 } = {\mathbf x}^n_i +  \frac{{\Delta t}}{a^{{n+\frac{1}{2}}}}  {\mathbf u}_i^{{n+\frac{1}{2}}}`

-  Solve for :math:`{\mathbf g}^{n+1}` using :math:`{\mathbf x}_i^{n+1}`

-  :math:`{\mathbf u}_i^{n+1} = \frac{1}{a^{n+1}} ( (a^{{n+\frac{1}{2}}} {\mathbf u}^{{n+\frac{1}{2}}}_i) + \frac{{\Delta t}}{2} \; {\mathbf g}^{n+1}_i )`

Note that at the end of the timestep :math:`{\bf x}_i^{n+1}` is consistent with :math:`{\bf g}^{n+1}` becasue
we have not advanced the positions after computing the new-time gravity. This has the benefit that
we perform only one gravity solve per timestep (in a single-level calculation with no hydro) because
the particles are only moved once.

Computing **g**
~~~~~~~~~~~~~~~

We solve for the gravitational vector as follows:

-  Assign the mass of the particles onto the grid in the form of density, :math:`\rho_{DM}`.
   The mass of each particle is assumed to be uniformly distributed over a cube of side :math:`\Delta x`,
   centered at what we call the position of the particle. We distribute the mass of each
   particle to the cells on the grid in proportion to the volume of the intersection of each cell
   with the particle’s cube. We then divide these cell values by :math:`\Delta x^3` so that the
   right hand side of the Poisson solve will be in units of density rather than mass.
   Note that this is the *comoving* density.

-  Solve :math:`\nabla^2 \phi = \frac{4 \pi G}{a} \rho_{DM}`.
   We discretize with the standard 7-point Laplacian (5-point in 2D)
   and use multigrid with Gauss-Seidel red-black relaxation to solve the equation for :math:`\phi` at cell centers.

-  Compute the normal component of :math:`{\bf g} = -\nabla \phi` at cell faces by differencing the adjacent values of :math:`\phi,`
   e.g. if :math:`{\bf g} = (g_x, g_y, g_z),` then we define :math:`g_x` on cell faces with a normal in the x-direction by computing
   :math:`g_{x,i-{\frac{1}{2}},j,k} = -(\phi_{i,j,k} - \phi_{i-1,j,k}) / \Delta x.`

-  Interpolate each component of :math:`{\bf g}` from normal cell faces onto each particle position using
   linear interpolation in the normal direction.

Output Format
=============

Checkpoint Files
----------------

  The particle positions and velocities are stored in a binary file in each checkpoint directory.
  This format is designed for being read by the code at restart rather than for diagnostics.
  We note that the value of :math:`a` is also written in each checkpoint directory,
  in a separate ASCII file called *comoving_a*, containing only the single value.

Particle Data in Plot Files
----------------------------

The particle positions and velocities will be written in binary files in each plotfile directory.
Dark matter particles will be in DM, active galactic nuclei particles will be in AGN,
neutrino particles will be in NPC.

  In addition, we can also
  visualize the particle locations as represented on the grid. There are multiple “derived quantities”
  which represent the particles. Including particle variables in the derived variables will make them
  be written as plotfile fields on the grid, i.e.::
  
    amr.derive_plot_vars = particle_count particle_mass_density 

  in the inputs file will generate plotfiles with only two variables.
  **particle_count** represents the number of particles in a grid cell;
  **particle_mass_density** is the density on the grid resulting from the particles.

  The same naming convention follows for particle velocities on the grid: **particle_x_velocity**,
  **particle_y_velocity**, **particle_z_velocity**

  Derived variables with **particle_** represent quantities from the Dark Matter Particle Container.
  Similar variables from the AGN particle container, and the Neutrino Particle Container
  are named **agn_** and **neutrino_**. Note these are particle fields written to the grid,
  which are distinct from the **density** field in the plotfile, which is baryonic density on the grid.

  We note that the value of :math:`a` is also written in each plotfile directory,
  in a separate ASCII file called *comoving_a*, containing only the single value.

.. role:: cpp(code)
   :language: c++

.. _sec:grid_creation:

Grid Creation
-------------

To run Nyx you must specify :cpp:`n_cell` in the inputs file -- 
this is the number of cells spanning the domain in each coordinate direction at level 0.

Users often specify :cpp:`max_grid_size` as well. The default load balancing algorithm then divides the 
domain in every direction so that each grid is no longer than :cpp:`max_grid_size` in that direction.
If not specified by the user, :cpp:`max_grid_size` defaults to 128 in 2D and 32 in 3D (in each coordinate direction).

Another popular input is :cpp:`blocking_factor`.  The value of :cpp:`blocking_factor` 
constrains grid creation in that in that each grid must be divisible by :cpp:`blocking_factor`.  
Note that both the domain (at each level) and :cpp:`max_grid_size` must be divisible by :cpp:`blocking_factor`
If not specified by the user, :cpp:`blocking_factor` defaults to 8 in each coordinate direction.
The typical purpose of :cpp:`blocking_factor` is to ensure that the grids will be 
sufficiently coarsenable for good multigrid performance.

There is one more default behavior to be aware of.  There is a boolean :cpp:`refine_grid_layout` 
that defaults to true but can be over-ridden at run-time.
If :cpp:`refine_grid_layout` is true and the number of grids created is less than the number of processors 
(Ngrids < Nprocs), then grids will be further subdivided until Ngrids >= Nprocs.

Caveat: if subdividing the grids to achieve Ngrids >= Nprocs would violate the 
:cpp:`blocking_factor` criterion then additional grids are not created and the 
number of grids will remain less than the number of processors

Note that :cpp:`n_cell` must be given as three separate integers, one for each coordinate direction.

However, :cpp:`max_grid_size` and :cpp:`blocking_factor` can be specified as a single value 
applying to all coordinate directions, or as separate values for each direction.  

 - if :cpp:`max_grid_size` (or :cpp:`blocking_factor`) is specified as multiple integers then the first 
   integer applies to level 0, the second to level 1, etc.  If you don't specify as many
   integers as there are levels, the final value will be used for the remaining levels.

 - if different values of :cpp:`max_grid_size` (or :cpp:`blocking_factor`) are wanted for each coordinate direction, 
   then :cpp:`max_grid_size_x`, :cpp:`max_grid_size_y` and :cpp:`max_grid_size_z` 
   (or :cpp:`blocking_factor_x`, :cpp:`blocking_factor_y` and :cpp:`blocking_factor_z`) must be used.  
   If you don't specify as many integers as there are levels, the final value will be used for the remaining levels.

 - note that :cpp:`max_grid_size` is just an upper bound; with :cpp:`n_cell = 48` 
   and :cpp:`max_grid_size = 32`, we will typically have one grid of length 32 and one of length 16.

The grid creation process at level 0 proceeds as follows (if not using the KD-tree approach):

#. The domain is initially defined by a single grid of size :cpp:`n_cell`.

#. If :cpp:`n_cell` is greater than :cpp:`max_grid_size` then the grids are subdivided until
   each grid is no longer than  :cpp:`max_grid_size` cells on each side.  The :cpp:`blocking_factor` criterion
   (ie that the length of each side of each grid is divisible by :cpp:`blocking_factor` in that direction)
   is satisfied during this process.

#. Next, if :cpp:`refine_grid_layout = true` and there are more processors than grids
   at this level, then the grids at this level are further divided until Ngrids >= Nprocs
   (unless doing so would violate the :cpp:`blocking_factor` criterion).

===================
Refinement Criteria
===================

The dynamic creation and destruction of grid levels is a fundamental part of `Nyx`'s capabilities. The
process for this is described in some detail in the `AMReX` documentation, but we summarize the key points
here.

At regular intervals (set by the user), each Amr level that is not the finest allowed for the run
will invoke a "regrid" operation.  When invoked, a list of error tagging functions is traversed. For each,
a field specific to that function is derived from the state over the level, and passed through a kernel
that "set"'s or "clear"'s a flag on each cell.  The field and function for each error tagging quantity is
identified in the setup phase of the code where the state descriptors are defined (i.e., in `Nyx_setup.cpp`).
Each function in the list adds or removes to the list of cells tagged for refinement. This final list of tagged
cells is sent to a grid generation routine, which uses the Berger-Rigoutsos algorithm to create rectangular grids
which will define a new finer level (or set of levels).  State data is filled over these new grids, copying where
possible, and interpolating from coarser level when no fine data is available.  Once this process is complete,
the existing Amr level(s) is removed, the new one is inserted into the hierarchy, and the time integration
continues.

The traditional `AMReX` approach to setting up and controlling the regrid process involves explicitly
creating ("hard coding") a number of functions directly into `Nyx`'s setup code. (Consult the source code
and `AMReX` documentation for precisely how this is done).  `Nyx` provides a limited capability to augment
the standard set of error functions that is based entirely on runtime data specified in the inputs (ParmParse)
data.  The following example portion of a ParmParse'd input file demonstrates the usage of this feature:

::

   amr.refinement_indicators = denerr dengrad velgrad

   amr.denerr.max_level = 3
   amr.denerr.value_greater = 3
   amr.denerr.field_name = density
   amr.dengrad.max_level = 1
   amr.dengrad.adjacent_difference_greater = 0.01
   amr.dengrad.field_name = density

   amr.velgrad.max_level = 2
   amr.velgrad.adjacent_difference_greater = 0.01
   amr.velgrad.field_name = x_velocity

Here, we have added five new custom-named criteria -- ``denerr``: cells with the density defined on the mesh greater than 3; ``dengrad``: cells having a density difference of 0.01 from that of their
immediate neighbor and ``velgrad``: cells having a x_velocity difference of 0.01 from that of their
immediate neighbor. 
The first will trigger up to Amr level 3, the second only to level 1, and the third to level 2.

An example of a more specific derived type is overden, a rescaling of the baryonic gas density divided by the average total density:

.. math::

   \begin{aligned}
   \mathtt{overden} &=&\frac{\rho_b}{ \overline{\rho} * \mathtt{tagging\_base}^{\mathtt{level+1}}} .\cr\cr \end{aligned}

where :math:`\overline{\rho}` is the average of :math:`\rho=\rho_b+\rho_{dm}` over the entire domain.

::

   amr.refinement_indicators = density
   amr.density.value_greater = 1
   amr.density.field_name = overden
   nyx.tagging_base = 1.1

Note that these additional user-created criteria operate in place of those defined as defaults.  Also note that
these can be modified between restarts of the code.  By default, the new criteria will take effect at the next
scheduled regrid operation.  Alternatively, the user may restart with ``amr.regrid_on_restart = 1`` in order to
do a full (all-levels) regrid after reading the checkpoint data and before advancing any cells.


==============================================
Hydrodynamic Equations in Comoving Coordinates
==============================================

Conservative Form
-----------------

We solve the equations of gas dynamics in a coordinate system that is comoving
with the expanding universe, with expansion factor, :math:`a,` related to the redshift, :math:`z`, by :math:`a = 1 / (1 + z).`
The continuity equation is written,

.. math::

   \label{eq:dens}
   \frac{\partial \rho_b}{\partial t} = - \frac{1}{a} \nabla \cdot (\rho_b {\bf U})  , \\

where :math:`\rho_b` is the comoving baryonic density, related to the proper density by :math:`\rho_b = a^3 \rho_{proper},`
and :math:`{\bf U}` is the proper peculiar baryonic velocity.

The momentum evolution equation can be expressed as

.. math::

   \begin{aligned}
   \frac{\partial (\rho_b {\bf U})}{\partial t} &=&  \frac{1}{a} \left(
   - \nabla \cdot (\rho_b {\bf U} {\bf U}) 
   - \nabla p 
   + \rho_b {\bf g} 
   + {\bf S}_{\rho {\bf U}}
   - \dot{a} \rho_b {\bf U} \right)  , \end{aligned}

or equivalently,

.. math::

   \begin{aligned}
   \label{eq:momt}
   \frac{\partial (a \rho_b {\bf U})}{\partial t} &=& 
   -             \nabla \cdot (\rho_b {\bf U} {\bf U}) 
   -             \nabla p 
   +             \rho_b {\bf g} 
   +             {\bf S}_{\rho {\bf U}}  , \end{aligned}

where the pressure, :math:`p`, that appears in the
evolution equations is related to the proper pressure, :math:`p_{proper},` by :math:`p = a^3 p_{proper}.`
Here :math:`{\bf g} = - \nabla \phi` is the gravitational acceleration vector, and
:math:`{\bf S}_{\rho {\bf U}}` represents any external forcing terms.

The energy equation can be written,

.. math::

   \begin{aligned}
   \frac{\partial (\rho_b E)}{\partial t} &=& \frac{1}{a} \left[
   - \nabla \cdot (\rho_b {\bf U} E + p {\bf U})
   + ( \rho_b {\bf U} \cdot {\bf g} +  S_{\rho E} ) 
   - \dot{a} ( 3 (\gamma - 1) \rho_b e + \rho_b ( {\bf U} \cdot {\bf U}) ) \right]  . \end{aligned}

or equivalently,

.. math::

   \begin{aligned}
   \label{eq:energy}
   \frac{\partial (a^2 \rho_b E)}{\partial t} &=& a \left[
   - \nabla \cdot (\rho_b {\bf U} E + p {\bf U})
   +  \rho_b {\bf U} \cdot {\bf g} 
   +  S_{\rho E}  
   +  \dot{a} ( \; ( 2 - 3 (\gamma - 1) ) \; \rho_b e ) \right]  . \end{aligned}

Here :math:`E = e + {\bf U} \cdot {\bf U} / 2` is the total energy per unit mass,
where :math:`e` is the specific internal energy.
:math:`S_{\rho E} = S_{\rho e} + {\bf U} \cdot {\bf S}_{\rho {\bf U}}`
where :math:`S_{\rho e} = \Lambda^H - \Lambda^C` represents the heating and cooling terms, respectively.
Additionally, :math:`{\bf S}_{\rho {E}}` may include any external forcing terms on the total energy, for example as in the stochastic forcing application.
We can write the evolution equation for internal energy as

.. math::

   \begin{aligned}
   \frac{\partial (\rho_b e)}{\partial t} &=& \frac{1}{a} \left[
   - \nabla \cdot (\rho_b {\bf U} e)
   - p \nabla \cdot {\bf U}
   - \dot{a} ( 3 (\gamma - 1) \rho_b e )
   + S_{\rho e}  \right]  . \end{aligned}

or equivalently,

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho_b e)}{\partial t} &=&  a \left[
   - \nabla \cdot (\rho_b {\bf U} e)
   - p \nabla \cdot {\bf U}
   + S_{\rho e} 
   + \dot{a} ( \; ( 2 - 3 (\gamma - 1) ) \; \rho_b e ) \right]  . \end{aligned}

Note that for a gamma-law gas with :math:`\gamma = 5/3,` we can write

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho_b E)}{\partial t} &=&  a \left[
    -\nabla \cdot (\rho_b {\bf U} E + p {\bf U})
   +  \rho_b {\bf U} \cdot {\bf g} 
   +  S_{\rho e}  \right]   . \end{aligned}

and

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho_b e)}{\partial t} &=& a \left[ 
   - \nabla \cdot (\rho_b {\bf U} e)
   -  p \nabla \cdot {\bf U}
   +  S_{\rho e}  \right]   . \end{aligned}

Work in Progress
================

Here we list work in progress; code pieces related to these efforts are already in the GitHub repo
but the functionality has not been finished yet and is thus not fomally relesed to the
community.  We provide this list for completness and clarity of the source code.


Active Galactic Nuclei
----------------------

We are implmenting Active Galactic Nuclei (AGN) subgrid model based on Bondi accretion
of gas onto a blackhole, and feedback of energy into the surrounding medium.  Ongoing
development is (mostly) contained inside Source/AGN/ folder.


Massive neutrinos
-----------------

Simulating massive neutrinos as particles is implemented in Nyx via separate
particle container.  Tests are ongoing.


Units and Conventions
=====================

Nyx uses a cosmological system of units based on Mpc, M\ :math:`_\odot`, and km/s for length, mass, and velocity,
respectively.  Temperature is given in degrees Kelvin.  The unit of time is derived from velocity and length units.
All inputs and problem initialization should be specified consistently with these units,
and the outputs should be interpreted with these units in mind.
In the equation of state calls, there is an internal
conversion between cosmological and CGS units.

:numref:`table:units` shows 
some of the common symbols / names used throughout the code documentation and papers.

.. _table:units:

.. table:: Common quantities and units.

   +-----------------------+--------------------------------------------------------+-----------------------+
   | name                  | units                                                  | description           |
   +=======================+========================================================+=======================+
   | :math:`t`             | (Mpc/km) s                                             | time                  |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`\rho`          | M\ :math:`_\odot` / Mpc\ :math:`^3`                    | mass density          |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`\ub`           | km/s                                                   | velocity vector       |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`p`             | M\ :math:`_\odot` (km/s)\ :math:`^2` / Mpc\ :math:`^3` | pressure              |
   |                       |                                                        |                       |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`\gb`           | (km/s)2 / Mpc                                          | gravitational         |
   |                       |                                                        | acceleration          |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`E`             | (km/s)\ :math:`^2`                                     | specific total energy |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`e`             | (km/s)\ :math:`^2`                                     | specific internal     |
   |                       |                                                        | energy                |
   +-----------------------+--------------------------------------------------------+-----------------------+
   | :math:`T`             | :math:`K`                                              | temperature           |
   +-----------------------+--------------------------------------------------------+-----------------------+


In :numref:`table:constants` we list the values used for physical constants in cosmological units.
Note that :math:`\Omega_m`, :math:`\Omega_b`, :math:`\Omega_r`  and :math:`h` are set in the inputs file :ref:`Comoving: List of Parameters<list-of-parameters-13>`.
Full list of constants and conversion factors is set in Source/Driver/constants_cosmo.H.

.. _table:constants:
.. table::
	   Physical constant values

   +-------------------------------------+-----------------------------------------------------------+
   | Constant                            | Cosmological units                                        |       
   +=====================================+===========================================================+
   | Gravitational constant (:math:`G`)  | 4.3019425e-9 Mpc (km/s)\ :math:`^2` M\ :math:`_\odot^{-1}`|
   +-------------------------------------+-----------------------------------------------------------+
   | Avogadro’s number (:math:`n_A`)     | 1.1977558e57 M\ :math:`_\odot^{-1}`                       |
   +-------------------------------------+-----------------------------------------------------------+
   | Boltzmann’s constant (:math:`k_B`)  | 0.6941701e-59 M\ :math:`_\odot` (km/s)\ :math:`^2` / K    |
   +-------------------------------------+-----------------------------------------------------------+
   | Hubble constant (:math:`H`)         | 100 (km/s) / Mpc                                          |
   +-------------------------------------+-----------------------------------------------------------+

The only other place that dimensional numbers are used in the code is in the tracing and Riemann solve.
We set three *small* numbers which need to be consistent with the data specified
We set one *large* number which needs to be consistent with the data specified.
Each of these can be specified in the inputs file.

-  small_dens – small value for density

-  small_press – small value for pressure

-  small_temp – small value for temperature

-  large_temp – large value for temperature

These are the places that each is used in the code:

-  **small_dens**

   -  | **enforce_minimum_density** (called after subroutine consup) – there are two choices for this. In the flooring routine, 
      | **enforce_minimum_density_floor** – density is set to small_dens, (rho e) and (rho E) are computed from small_temp,
      | and momenta are set to zero.  In the conservative routine, **enforce_minimum_density_cons**, an iterative procedure 
      | is used to create diffusive fluxes that adjusts all the variables conservatively until density is greater than small_dens.

   -  | **tracexy / tracez / tracexy_ppm / tracez_ppm**:
      | qxp = max(qxp,small_dens)
      | qxm = max(qxm,small_dens)
      | and analogously for qyp/qym and qzp/qzm. This modifies the primitive density and pressure inside the tracing, not the underlying state variables

   -  **riemannus** – we set

      wsmall = small_dens \* csmall

      and then

      | wl = max(wsmall, sqrt(gaml \* pl \* rl))
      | wr = max(wsmall, sqrt(gamr \* pr \* rr))

      Also, we set

      ro = max(small_dens,ro)

      where ro = 0.5 \* (rl + rr) – this state is only chosen when ustar = 0, and

      rstar = max(small_dens,rstar)

      where rstar = ro + (pstar-po)/co:math:`^2`

-  **small_temp**:

   -  | **compute_new_temp**: if :math:`\rho e < 0`, then
      | call nyx_eos_given_RT (e,...,small_temp,...) in order to compute a new energy, :math:`e`.
      | This energy is then used to define a new :math:`E = e + ke`
      | Coming out of this the temperature is equal to small_temp and the energy :math:`e` has been reset.

   -  | **reset_internal_energy**: if :math:`e < 0` and :math:`E - ke < 0` then
      | call nyx_eos_given_RT (e,...,small_temp,...) in order to compute a new energy, :math:`e`. This energy is also used to define a new :math:`E = e + ke`

-  **large_temp**:

   -  | **compute_new_temp**: if :math:`T > \mathrm{large\_temp}`, and the input flag ``nyx.local_max_temp_dt=1`` then
      | set :math:`T = \mathrm{large\_temp}` and call nyx_eos_given_RT (e,...,large_temp,...) in order to compute a new energy, :math:`e`.
      | This energy is then used to define a new :math:`E = e + ke`
      | Coming out of this the temperature is equal to large_temp and the energy :math:`e` has been reset.

-  **small_pres**:

   -  | **tracexy / tracez / tracexy_ppm / tracez_ppm**:
      | qpres = max(qpres,small_pres)
      | for qxp/qyp, qyp/qym and qzp/qzm. This modifies the primitive density and pressure inside the tracing, not the underlying state variables

   -  **riemannus** – we set

      | pstar = max(small_pres,pstar)
      | pgdnv = max(small_pres,pgdnv). Note that pgdnv is the pressure explicitly used in the fluxes.

   -  **uflatten** – small_pres is used to keep the denominator away from zero

   -  Everywhere we define values of pressure on a face, we set that value to be at least small_pres.
.. role:: cpp(code)
   :language: c++

.. _sec:load_balancing:

Load Balancing
--------------

The process of load balancing is typically independent of the process of grid creation; 
the inputs to load balancing are a given set of grids with a set of weights 
assigned to each grid.

Single-level load balancing algorithms are sequentially applied to each AMR level independently, 
and the resulting distributions are mapped onto the ranks taking into account the weights 
already assigned to them (assign heaviest set of grids to the least loaded rank)

Options supported by AMReX include the following; the default is SFC:

- Knapsack: the default weight of a grid in the knapsack algorithm is the number of grid cells, 
  but AMReX supports the option to pass an array of weights – one per grid – or alternatively 
  to pass in a MultiFab of weights per cell which is used to compute the weight per grid

- SFC: enumerate grids with a space-filling Z-morton curve, then partition the 
  resulting ordering across ranks in a way that balances the load

- Round-robin: sort grids and assign them to ranks in round-robin fashion -- specifically
  FAB *i* is owned by CPU *i* %N where N is the total number of MPI ranks.

Load Balancing the Hydrodynamical Mesh
--------------------------------------

For Nyx, the DistributionMapping defaults to using a SFC strategy.
Setting the DistributionMapping to use a different strategy is controlled
by flags beginning with DistributionMapping::

  DistributionMapping.strategy = {KNAPSACK, SFC, ROUNDROBIN}
  DistributionMapping.verbose = {0, 1}

These flags take effect whenever the Regrid operation is called on the mesh,
typically either when ``amr.regrid_int`` is reached in a multilevel simulation or if
``amr.regrid_on_restart=1``.

Load Balancing the Dark Matter Particles
----------------------------------------

For Nyx, the flags which affect how the active particles get their DistributionMapping change are::

  nyx.load_balance_int = -1
  nyx.load_balance_start_z = 15
  nyx.load_balance_wgt_strategy = {0, 1, 2}
  nyx.load_balance_wgt_nmax = -1
  nyx.load_balance_strategy = {KNAPSACK, SFC, ROUNDROBIN}

The ``wgt_strategy`` uses either the 0=number of cells, 1=number of dark matter particles, or 2=total count of dark matter particles to determine the cost of each box. 1 and 2 should be equivalent, and 1 should be cheaper.

This load balancing is not on by default, and has a flag to determine based on redshift z when to begin affecting the calculation because this introduces additional communication penalties when the particles interact with mesh data. This additional communication is a trade-off to balance the amount of particle memory and particle work between ranks.

This strategy is not currently implemented to work for multilevel simulations.
The AGN Model in Nyx
====================

In the AGN model, super-massive black hole (SMBH) particles are formed
at *haloes*, where each halo is defined by a connected mass enclosed by
a user-defined density isocontour. In order to find haloes, we use the
Reeber package described in Section \ `2 <#sec:Reeber>`__. Each AGN
particle has the standard dark matter particle attributes of position,
velocity, and mass, as well as two additional attributes, its stored
accretion energy and its mass accretion rate.

.. table:: Parameters of the AGN model

   ================ ============================ ======================== ===========================================
   In “probin” file Parameter                    Fiducial value           Explanation
   ================ ============================ ======================== ===========================================
   \*               :math:`M_{\rm h, min}`       :math:`10^{10}~ M_\odot` Minimum halo mass for SMBH placement
   \*               :math:`M_{\rm seed}`         :math:`10^5~ M_\odot`    Seed mass of SMBH
   T_min            :math:`T_{\rm min}`          :math:`10^7` K           Minimum heating of the surrounding gas
   bondi_boost      :math:`\alpha`               100                      Bondi accretion boost factor
   max_frac_removed :math:`f_{\rm max, removed}` 0.5                      Maximum fraction of mass removed from gas
   eps_rad          :math:`\epsilon_{\rm r}`     0.1                      Radiation efficiency
   eps_coupling     :math:`\epsilon{\rm c}`      0.15                     Coupling efficiency
   eps_kinetic      :math:`\epsilon_{\rm kin}`   0.1                      Kinetic feedback efficiency
   frac_kinetic     :math:`f_{\rm kin}`          0                        Fraction of feedback energy that is kinetic
   ================ ============================ ======================== ===========================================

[tab:agn_params1]

\* :math:`M_{\rm h, min}` and :math:`M_{\rm seed}` are not set in the
“probin” file, but in the inputs file, by respectively Nyx.mass_halo_min
and Nyx.mass_seed.

Creating AGN Particles from Haloes
----------------------------------

Each halo with threshold mass of :math:`M_h \geqslant M_{\rm h, min}`
that does not already host a black hole particle is seeded with a black
hole of mass :math:`M_{\rm seed}`. The initial position of this AGN
particle is the center of the cell where the density is highest in the
halo.

When an AGN particle is created, the density in its cell is reduced by
the amount required for mass to be conserved, and the velocity of the
AGN particle is initialized so that momentum is conserved. The accretion
energy and mass accretion rate are initialized to zero.

Merging AGN Particles
---------------------

Two AGN particles merge when both of these conditions obtain:

#. The distance between them, :math:`l`, is less than the mesh spacing,
   :math:`h`.

#. [velocity-item] The difference of their velocities,
   :math:`v_{\rm rel}`, is less than the circular velocity at distance
   :math:`l`:

   .. math:: v_{\rm rel} < \sqrt{GM_{\rm BH}/l}

   where :math:`M_{\rm BH}` is the mass of the more massive SMBH in the
   pair, and :math:`G` is the gravitational constant.

Criterion \ `[velocity-item] <#velocity-item>`__ above is necessary in
order to prevent AGN particles from merging during a fly-through
encounter of two haloes, as this could lead to AGN particles being
quickly removed from the host halo due to momentum conservation.

The merger of two AGN particles is implemented as the less massive one
being removed, and its mass and momentum being transferred to the more
massive one.

Accretion
---------

For an AGN particle of mass :math:`M_{\rm BH}`, the Bondi–Hoyle
accretion rate is

.. math::

   \dot{M}_{\rm B} = \alpha
   \frac{4 \pi G^2 M_{\rm BH}^2 \overline{\rho}}{(\overline{c_s^2} + \overline{u^2})^{3/2}} ,

where :math:`\overline{\rho}`, :math:`\overline{c_s^2}`, and
:math:`\overline{u^2}` are volume averages with a cloud-in-cell stencil
of the gas’s density, squared sound speed, and squared velocity,
respectively, in the neighborhood of the particle.

The maximum black hole accretion rate is the Eddington limit,

.. math::

   \dot{M}_{\rm Edd} = 
   \frac{4 \pi G M_{\rm BH} m_{\rm p}}{\epsilon_{\rm r} \sigma_{\rm T} c} \, ,

with proton mass :math:`m_{\rm p}`, Thomson cross section
:math:`\sigma_{\rm T}`, and speed of light :math:`c`.

The mass accretion rate of the SMBH is the smaller of the two rates
above:
:math:`\dot{M}_{\rm acc} = {\rm min} \{ \dot{M}_{\rm B}, \dot{M}_{\rm Edd} \}`.
Then the gas will lose mass :math:`\dot{M}_{\rm acc} \Delta t`, where
:math:`\Delta t` is the length of the time step. However,
:math:`\dot{M}_{\rm acc}` is adjusted downward if necessary so that when
cloud-in-cell stencil weights are applied in the neighborhood of the
particle, the fraction of gas density removed from any cell of the
stencil is at most :math:`f_{\rm max, removed}`.

The mass of the AGN particle increases by
:math:`(1-\epsilon_{\rm r}) \dot{M}_{\rm acc} \Delta t`, while
:math:`\dot{M}_{\rm acc} \Delta t` amount of gas mass is removed from
the grid according to cloud-in-cell stencil weights in the neighborhood
of the particle. The momentum transfer can be computed by assuming the
velocity of the gas is unchanged; thus the gas in each cell loses
momentum in proportion to its mass loss, and the particle gains the sum
of the gas momentum loss multiplied by :math:`(1-\epsilon_{\rm r})`.

Feedback Energy
---------------

Feedback energy is stored in an AGN particle variable
:math:`E_{\rm AGN}`, and is accumulated over time until released. The
fraction :math:`f_{\rm kin}` goes to kinetic energy, and the rest to
thermal energy.

Thermal Feedback
~~~~~~~~~~~~~~~~

We increment :math:`E_{\rm AGN}` by thermal feedback energy, calculated
from the mass accretion rate as

.. math::

   % high energy
   \Delta E_{\rm thermal} = (1 - f_{\rm kin})
   \epsilon_{\rm c} \epsilon_{\rm r} \dot{M}_{\rm acc} c^2 \Delta t .

Kinetic/Momentum Feedback
~~~~~~~~~~~~~~~~~~~~~~~~~

We increment :math:`E_{\rm AGN}` by the kinetic feedback energy

.. math::

   % low energy
   \Delta E_{\rm kinetic} =
   f_{\rm kin} \epsilon_{\rm kin} \dot{M}_{\rm acc} c^2 \Delta t .

We also need to adjust the energy density and momentum density of the
gas. We do this by computing a jet velocity

.. math:: \vec{v}_{\rm jet} = \sqrt{\frac{2 \Delta E_{\rm kinetic}}{m_{\rm g}}} \vec{n}

where :math:`m_{\rm g}` is the total gas mass inside the cloud-in-cell
local environment, and :math:`\vec{n}` is a *randomly* chosen unit
vector. We add :math:`\rho \vec{v}` to the momentum density
:math:`\vec{p}` of the gas, and :math:`\vec{v}_{\rm jet} \cdot \vec{p}`
to its energy density, both of these weighted by the cloud-in-cell
stencil of the particle.

Releasing Feedback Energy
~~~~~~~~~~~~~~~~~~~~~~~~~

The accumulated energy is released when

.. math::

   E_{\rm AGN} > m_{\rm g} \overline{e}
   \label{eq:E_agn}

where :math:`\overline{e}` is the average specific internal energy of
the gas over the cloud-in-cell stencil, obtained from the equation of
state using temperature :math:`T_{\rm min}` and average density of the
gas over the same stencil, and :math:`m_{\rm g}` is the total gas mass
inside the cloud-in-cell local environment.

.. _sec:Reeber:

The Reeber Package
==================

Reeber is a separate package with a halo finder. Here are the Reeber
parameters that are assigned in the input file.

=================================== =================================== ===================================== =========
Parameter                           Definition                          Acceptable Values                     Default
=================================== =================================== ===================================== =========
**reeber.halo_int**                 timesteps between halo finder calls Integer                               -1 (none)
**reeber.negate**                   allow negative values for analysis  0 if false, 1 if true                 1
**reeber.halo_density_vars**        density variable list               density, particle_mass_density        “density”
**reeber.halo_extrema_threshold**   extrema threshold for haloes        Real                                  200.
**reeber.halo_component_threshold** component threshold for haloes      Real                                  82.
**reeber.absolute_halo_thresholds** are halo thresholds absolute        0 if multiples of mean, 1 if absolute 0
=================================== =================================== ===================================== =========

[Table:Reeber-inputs]
*******
Gravity
*******
Introduction
============

In order to use gravity, one must set
  
  nyx.do_grav= 1
  
in the inputs file. See :ref:`Gravity: List of Parameters<list-of-parameters-11>` for relevant input flags.

Poisson Approximation
=====================

In Nyx we always compute gravity by solving a Poisson equation on the mesh hierarchy.
To define the gravitational vector we set

  .. math:: \mathbf{g}(\mathbf{x},t) = -\nabla \phi

  where

  .. math:: \mathbf{\Delta} \phi = \frac{4 \pi G}{a} (\rho - \overline{\rho}) \label{eq:Self Gravity}

where :math:`\overline{\rho}` is the average of :math:`\rho` over the entire domain if we assume triply periodic boundary conditions,
and :math:`a(t)` is the scale of the universe as a function of time.

Time Integration Strategy
=========================

Nyx uses subcycling to integrate levels at different timesteps.
The gravity algorithm needs to respect this. Self-gravity is computed
via multigrid. At coarse-fine interfaces, the stencil used in the
Laplacian understands the coarse-fine interface and is different than
the stencil used in the interior.

There are two types of solves that we discuss with AMR:

-  *composite solve* : This is a multilevel solve, starting at
   a coarse level (usually level 0) and solving for the potential on
   all levels up to the finest level.

-  *level solve* : This solves for the potential only on
   a particular level. Finer levels are ignored. At coarse-fine
   interfaces, the data from the coarse levels acts as Dirichlet
   boundary conditions for the current-level-solve.

Briefly:

-  At the beginning of a simulation, we do a multilevel composite
   solve (if ``gravity.no_composite`` = 0).

   We also do a multilevel composite solve after each regrid.

-  The old-time gravity on the coarse level is defined based on
   this composite solve, but we also do a level solve on the coarse
   level, and use it to determine the difference between the composite
   solve and the level solve, and store that in a MultiFab.

-  After the hydro advance on the coarse level, we do another level
   solve, and use the (level solve - compositive solve) as a lagged
   predictor of how much we need to add back to that level solve to get
   an effective new-time value for phi on the coarse level, and that’s
   what defines the phi used for the new-time gravity

-  Then we do the fine grid timestep(s), each using the same
   strategy

-  At an AMR synchronization step across levels, if we’re
   choosing to synchronize the gravitational field across levels
   (``gravity.no_sync`` = 0) we then do a solve starting from the coarse
   grid that adjusts for the mismatch between the fine-grid phi and
   the coarse-grid phi, as well as the mismatch between the fine-grid
   density fluxes and the coarse-grid density fluxes, and add the
   resulting sync solve phi to both the coarse and the fine level

   Thus, to within the gravity error tolerance, you get the same final
   result as if you had done a full composite solve at the end of the
   timestep (assuming ``gravity.no_sync`` = 0).

If you do ``gravity.no_composite`` = 1, then you never do a full
multilevel solve, and the gravity on any level is defined only by the
solve on that level. The only time this would be appropriate is if
the fine level(s) cover essentially all of the mass on the grid for
all time.

Hydrodynamics Source Terms
==========================

  We use a standard predictor-corrector formalism for updating the momentum and
  energy. Specifically, our first update is equal to :math:`\Delta t
  \times \mathbf{S}^n` , where :math:`\mathbf{S}^n` is the value of
  the source terms at the old-time (which is usually called time-level
  :math:`n`). At the end of the timestep, we do a corrector step where
  we subtract off :math:`\Delta t / 2 \times \mathbf{S}^n` and add on
  :math:`\Delta t / 2 \times \mathbf{S}^{n+1}`, so that at the end of
  the timestep the source term is properly time centered.
.. role:: cpp(code)
   :language: c++

Gridding and Load Balancing
===========================

Nyx has a great deal of flexibility when it comes to how to decompose the 
computational domain into individual rectangular grids, and how to distribute
those grids to MPI ranks.  There can be grids of different sizes, 
more than one grid per MPI rank, and different strategies for distributing the grids to MPI ranks.

We use the phrase "load balancing" here to refer to the combined process
of grid creation (and re-creation when regridding) and distribution of grids to MPI ranks.

See :ref:`sec:grid_creation` for grids are created, i.e. how the :cpp:`BoxArray` on which 
:cpp:`MultiFabs` will be built is defined at each level.

See :ref:`sec:load_balancing` for the strategies AMReX supports for distributing
grids to MPI ranks, i.e. defining the :cpp:`DistributionMapping` with which 
:cpp:`MultiFabs` at that level will be built.  

When running on multicore machines with OpenMP, we can also control the distribution of 
work by setting the size of grid tiles (by defining :cpp:`fabarray_mfiter.tile_size`).
We can also specify the strategy for assigning tiles to OpenMP threads.  
See the AMReX documentation for more about tiling.

.. toctree::
   :maxdepth: 1

   GridCreation
   LoadBalancing

 .. role:: cpp(code)
    :language: c++

.. _InSitu:

In situ Analysis
================

Nyx supports in situ visualization using Ascent, and can leverage AMReX's Sensei visualization pipeline. Additionally, Nyx supports in situ halo finding using Reeber2. This is useful both as
in situ analysis tool, and for subgrid models, like AGN (under development).

These insitu calls are controlled by input flags, in this example analysis happens starting after step 100, for step 199, 299, and so on and so forth::

  insitu.int = 100
  insitu.start = 100

Additionally, we can request analysis at specific redshifts, and Nyx will adjust the time-stepping to reach those red-shifts exactly::

  nyx.analysis_z_values = 7.0 6.0 5.0 4.0 3.0 2.0

Ascent visualization
--------------------

The primary example of this functionality is in `Nyx/Exec/LyA`. To compile with Ascent, add in GNUmakefile.summit::

  USE_ASCENT_INSITU = TRUE

The ``ascent_actions.yaml`` file will determine while the code is running what the Ascent publish action does. The ``ascent_actions_slicefile.yaml`` file gives an example of saving slice data, rather than an image file, while the simulation is running.

To build Ascent with Conduit for other machines and configurations, see `Building Ascent <https://ascent.readthedocs.io/en/latest/BuildingAscent.html>`_. For further information about using Ascent with AMReX-based applications, see `AMReX Blueprint Tutorial <https://amrex-codes.github.io/amrex/tutorials_html/Blueprint_Tutorial.html>`_ and `WarpX Ascent InSitu Documentation <https://warpx.readthedocs.io/en/latest/visualization/ascent.html>`_. 

Sensei interface
----------------

See AMReX documentation: `SENSEI <https://amrex-codes.github.io/amrex/docs_html/Visualization.html#sensei>`_

Halo finding
------------

To find halos in situ while Nyx is running we use the Reeber package.
To compile with Reeber, add in GNUmakefile::

  REEBER = TRUE

and set the location of Boost library::

  BOOST_INLUDE_DIR := ${BOOST_ROOT}/include/boost

GNUmake will default to::

  DIY_INCLUDE_DIR ?= ../../../diy/include
  REEBER_EXT_HOME ?= ../../../Reeber2

Note that these codes are in separate repositories and are not included in Nyx distribution.
If you intend to use in situ halo finding, you should clone Reeber from its
`GitHub page <https://github.com/mrzv/reeber>`_ and follow the installation instructions provided there.

In the Nyx inputs file, one should specify the time step interval for halo finder, fields which will be
used (to use the total density, one should specify both (gas) ``density`` and ``particle_mass_density``),
and thresholds of the boundary and extrema values::

  # Halo Finder
  reeber.halo_int = 1
  reeber.negate = 1
  reeber.halo_density_vars = density particle_mass_density
  reeber.halo_extrema_threshold = 20
  reeber.halo_component_threshold = 10
  
  # Call reeber insitu analysis
  insitu.reeber_int = 100


.. _note:
  These instructions are based on Reeber hash 8a274d35a415f7b15d8308a30763f52c4eeb7c7b, diy hash 88eca5107935b2d50eb352d99a6b0ed109b9c31c, and Nyx hash 33006ce18b1f945053c05a7cade0f4aba63378b5
============================================
Hydrodynamical and Heating-Cooling Splitting
============================================

We solve the equations of gas dynamics in a coordinate system that is comoving
with the expanding universe, with expansion factor, :math:`a,` related to the redshift, :math:`z`, by :math:`a = 1 / (1 + z).`

We describe the state of the gas
as :math:`\overline{{\bf U}} = (\rho, a \rho {\bf U}, a^2 \rho E, a^2 \rho e),`
then write the evolution of the gas as

.. math:: \frac{\partial\overline{{\bf U}}}{\partial t} = -\nabla\cdot{\bf F}+ S_e + S_g + S_{HC},

where :math:`{\bf F}= (1/a \; \rho {\bf U}, \rho {\bf U}{\bf U}, a (\rho {\bf U}E + p {\bf U}), a \rho {\bf U}e)`
is the flux vector,
:math:`S_e = (0, 0, 0, -a p \nabla \cdot {\bf U})` represents the additional term in the evolution
equation for internal energy, :math:`S_g = (0, \rho_b {\bf g}, a \rho {\bf U}\cdot {\bf g}, 0)`
represents the gravitational source terms,
and :math:`S_{HC} = (0, 0, a \rho \Lambda_{HC}, a \rho \Lambda_{HC})`
represents the combined heating and cooling source terms. The state, :math:`\overline{{\bf U}},` and
all source terms are defined at cell centers; the fluxes are defined on cell faces.

We compute :math:`{\bf F}`
using an unsplit Godunov method with characteristic tracing and full
corner coupling.

We track derived variables such as temperature and electron fraction ne based on the internal energy. The concentrations
of different isotopes of hydrogen and helium then depend on a Newton solve of the heating/cooling equations of state, which
are closely tied to the internal energy.

Strang Splitting
----------------

The original splitting used in Nyx is Strang splitting, where a half-step of the heating-cooling
is evolved, then a full step of the hydrodynamical terms, followed by a half-step of the heating-cooling.
This algorithm is the classical Strang splitting, adapted to include gravity and other forcing terms.

.. math::

   \label{eq:dens}
   \frac{\partial \rho}{\partial t} = - \frac{1}{a} \nabla \cdot (\rho {\bf U}) , \\

.. math::

   \begin{aligned}
   \label{eq:momt}
   \frac{\partial (a \rho {\bf U})}{\partial t} &=& 
   -             \nabla \cdot (\rho {\bf U}{\bf U}) 
   -             \nabla p 
   +             \rho {\bf g}, \end{aligned}

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho E)}{\partial t} &=&  a \left[
    -\nabla \cdot (\rho {\bf U}E + p {\bf U})
   +  \rho {\bf U}\cdot {\bf g}
   +  \rho \Lambda_{HC}  \right]  . \end{aligned}

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho e)}{\partial t} &=& a \left[ 
   - \nabla \cdot (\rho {\bf U}e)
   -  p \nabla \cdot {\bf U}
   +  \rho \Lambda_{HC}  \right]  . \end{aligned}

The algorithm at a single level of refinement begins by computing the time step, :math:`\Delta t,`
and advancing :math:`a` from :math:`t^n` to :math:`t^{n+1} = t^n + \Delta t`. The rest of the time step is
then composed of the following steps:

Step 1:
   Compute :math:`{\phi}^n` and :math:`{\bf g}^n` using :math:`\rho_b^n` and :math:`\rho_{dm}^n`,
   where :math:`\rho_{dm}^{n}` is computed from the particles at :math:`{\bf x}_i^{n}`.

   We note that in the single-level algorithm we can instead use :math:`{\bf g}` as computed at the
   end of the previous step because there have been no changes to :math:`{\bf x}_i`
   or :math:`\rho` since then.

Step 2:
   Interpolate :math:`{\bf g}^n` from the grid to the particle locations, then
   advance the particle velocities by :math:`\Delta t/ 2` and particle positions by :math:`\Delta t`.

   .. math::

      \begin{aligned}
           {\bf u}_i^{{n+\frac{1}{2}}} &=& \frac{1}{a^{{n+\frac{1}{2}}}} ((a^n {\bf u}^n_i) + \frac{\Delta t}{2} \; {\bf g}^n_i) \\
           {\bf x}_i^{n+1}  &=& {\bf x}^n_i + \frac{\Delta t}{a^{{n+\frac{1}{2}}}} {\bf u}_i^{{n+\frac{1}{2}}}\end{aligned}

Step 3:
   Advance :math:`{\bf U}` by :math:`\frac{1}{2}\Delta t` for first strang

   We advance :math:`e` by integrating the source terms in time for :math:`\frac{1}{2}\Delta t`

   .. math::

      \begin{aligned}
           ( e)^{n,\ast} &=& ( e)^n +  \int \Lambda_{HC} \; dt^\prime  .\end{aligned}

   We update :math:`(\rho e)^{n,\ast}=(\rho e)^{n,\ast}+\rho^{n}\left((e)^{n,\ast}-(e)^{n}\right)`

   We update :math:`(\rho E)^{n,\ast}=(\rho E)^{n,\ast}+\rho^{n}\left((e)^{n,\ast}-(e)^{n}\right)`

Step 4:
   Advance :math:`{\bf U}` by :math:`\Delta t` for advective terms

   Advance the solution using time-centered fluxes and :math:`S_e`
   and an explicit representation of :math:`S_g` at time :math:`t^n`:

   .. math:: {\bf U}^{n+1,\ast} = {\bf U}^{n,\ast} + \Delta tA^{n+\frac{1}{2}}+ \Delta tS_g^n

   where :math:`A^{n+\frac{1}{2}}` is computed by predicting from the :math:`{\bf U}^{n,\ast}` states.

   After adding the advective update terms :math:`A^{n+\frac{1}{2}}` and before calculating :math:`S_{g}`, apply a correction to :math:`{\bf U}^{n+1,\ast}` to enforce :math:`\rho > 1.1 \times small_dense`
		 
Step 5: 
   Second Strang step

   We advance :math:`e` by integrating the source terms in time for :math:`\frac{1}{2}\Delta t`

   .. math::

        \begin{aligned}
        ( e)^{n+1} &=& ( e)^{n+1,\ast } +  \int \Lambda_{HC} \; dt^\prime .\end{aligned}

   We update :math:`(\rho e)^{n+1}=(\rho e)^{n+1,\ast}+\rho^{n+1}\left((e)^{n+1}-(e)^{n+1,\ast}\right)`

   We update :math:`(\rho E)^{n+1}=(\rho E)^{n+1,\ast}+\rho^{n+1}\left((e)^{n+1}-(e)^{n+1,\ast}\right)`

   We store Ne and Temp based on eos\_ hc updates from :math:`(e)^{n+1}`

Step 6:
   Compute :math:`{\phi}^{n+1}` and :math:`{\bf g}^{n+1}` using
   :math:`\rho^{n+1,*}` and :math:`\rho_{dm}^{n+1}`, where :math:`\rho_{dm}^{n+1}`
   is computed from the particles at :math:`{\bf x}_i^{n+1}`.

   Here we can use :math:`{\phi}^n` as an initial guess for :math:`{\phi}^{n+1}` in order to reduce the time
   spent in multigrid to reach the specified tolerance.

Step 7:
   Correct :math:`{\bf U}` with time-centered source terms, and replace :math:`e` by
   :math:`E - \frac{1}{2}U^2` as appropriate.

   We time-center the
   gravitational source terms only,

   .. math:: {\bf U}^{n+1} = {\bf U}^{n+1} + \frac{\Delta t}{2} (S_g^{n+1} - S_g^n)

Step 8:
   Interpolate :math:`{\bf g}^{n+1}` from the grid to the particle locations, then
   update the particle velocities, :math:`{\bf u}_i^{n+1}`

   .. math::

      \begin{aligned}
          {\bf u}_i^{n+1} &=& \frac{1}{a^{n+1}}
                          \left( \left( a^{{n+\frac{1}{2}}} {\bf u}^{{n+\frac{1}{2}}}_i \right)
                               + \frac{\Delta t}{2} \; {\bf g}^{n+1}_i \right)  \end{aligned}

Step \**:
   in post\_ timestep, do a reset and compute\_ new\_ temp after syncing the gravity sources

Deferred-Correction Splitting Algorithm
---------------------------------------

This algorithm is based on the version of SDC introduced by Nonaka et. al. ``\cite{Nonaka2012}``

.. math::

   \frac{\partial \rho}{\partial t} = - \frac{1}{a} \nabla \cdot (\rho {\bf U}) , \\

.. math::

   \begin{aligned}
   \frac{\partial (a \rho {\bf U})}{\partial t} &=& 
   -             \nabla \cdot (\rho {\bf U}{\bf U}) 
   -             \nabla p 
   +             \rho {\bf g} , \end{aligned}

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho E)}{\partial t} &=&  a \left[
    -\nabla \cdot (\rho {\bf U}E + p {\bf U})
   +  \rho {\bf U}\cdot {\bf g}
   +  \rho \Lambda_{HC}  \right] . \end{aligned}

.. math::

   \begin{aligned}
   \frac{\partial (a^2 \rho e)}{\partial t} &=& a \left[ 
   - \nabla \cdot (\rho {\bf U}e)
   -  p \nabla \cdot {\bf U}
   +  \rho \Lambda_{HC}  \right] . \end{aligned}

The algorithm at a single level of refinement begins by computing the time step, :math:`\Delta t,`
and advancing :math:`a` from :math:`t^n` to :math:`t^{n+1} = t^n + \Delta t`. The rest of the time step is
then composed of the following steps:

Step 1:
   Compute :math:`{\phi}^n` and :math:`{\bf g}^n` using :math:`\rho^n` and :math:`\rho_{dm}^n`,
   where :math:`\rho_{dm}^{n}` is computed from the particles at :math:`{\bf x}_i^{n}`.

   We note that in the single-level algorithm we can instead use :math:`{\bf g}` as computed at the
   end of the previous step because there have been no changes to :math:`{\bf x}_i`
   or :math:`\rho` since then.

Step 2:
   Interpolate :math:`{\bf g}^n` from the grid to the particle locations, then
   advance the particle velocities by :math:`\Delta t/ 2` and particle positions by :math:`\Delta t`.

   .. math::

      \begin{aligned}
           {\bf u}_i^{{n+\frac{1}{2}}} &=& \frac{1}{a^{{n+\frac{1}{2}}}} ((a^n {\bf u}^n_i) + \frac{\Delta t}{2} \; {\bf g}^n_i) \\
           {\bf x}_i^{n+1}  &=& {\bf x}^n_i + \frac{\Delta t}{a^{{n+\frac{1}{2}}}} {\bf u}_i^{{n+\frac{1}{2}}}\end{aligned}

Step 3:
   Construct advective update terms using :math:`I_R` from last timestep as source

   .. math::

      \begin{aligned}
      A_{\rho} & = & -\frac{1}{a}\nabla\cdot(\rho{\bf U})\\
      A_{\rho u} & = & -\nabla\cdot\left(\rho uu\right)-\nabla p\\
      A_{\rho E} & = & a\left[-\nabla\cdot(\rho{\bf U}E+p{\bf U})\right]\\
      A_{\rho e} & = & \frac{1}{a} \left[
      - \nabla \cdot (\rho_b {\bf U}e)
      - p \nabla \cdot {\bf U}) \right]\end{aligned}

Step 4:
   Update momentum and :math:`\rho E`

   .. math::

      \begin{aligned}
      S_{g} & = & \rho g\\
            &  & \rho{\bf U}\cdot{\bf g}\end{aligned}

   .. math:: u^{n+1,\ast} = u^{n} + \Delta tA^{n+\frac{1}{2}}+ \Delta tS_g^n

   .. math:: \left(\rho E\right)^{n+1,\ast }=\left(\rho E\right)^{n}+ \Delta tA_{\rho E}^{n+1/2} + \Delta tS_g

   After adding the advective update terms :math:`A_\ast` to the appropriate components of :math:`{\bf U}^{n+1,\ast}` and before calculating :math:`S_{g}`, apply a correction to :math:`{\bf U}^{n+1,\ast}` and `A_{\rho}` to enforce :math:`\rho > 1.1 \times small_dens`

Step 5:
   Simultaneously solve heating-cooling:

   .. math::

      \begin{aligned}
      \rho^{n+1,\ast} & = & \rho^{n}+\int_{t^{n}}^{t^{n+1}}A_{\rho}dt^{\prime}\\
      e^{n+1,\ast} & = & e^{n}+\int_{t^{n}}^{t^{n+1}} \left(A_{e}+\Lambda_{HC}\right) dt^{\prime}\end{aligned}

   where :math:`A_{e}=\frac{1}{\Delta t}\left(\left(\left[\frac{1}{a^{n+1}}\right]^{2}\left(\left[a^{n}\right]^{2}\left(\rho e\right)^{n}+\Delta t*A_{\rho e}\right)+A_{reset}\right)/\left(\rho^{n}+\Delta tA_{\rho}\right)-e^{n}\right)`

Step 6:
   We define

   .. math::

      \begin{aligned}
      I_{R_{\left(\rho e\right)}} & = & \left( \left[a^{n+1}\right]^{2}\rho^{n+1,\ast}e^{n+1,\ast}-\left(\left[a^{n}\right]^{2}\rho^{n}e^{n}+\Delta tA_{\rho e}\right)\right)/\left[\Delta t\left(\frac{a^{n}+a^{n+1}}{2}\right)\right]\\
      & & -\left[a^{n+1}\right]^{2}A_{reset}/\left[\Delta t\left(\frac{a^{n}+a^{n+1}}{2}\right)\right]\end{aligned}

Step 7:
   We update internal and total energy using this forcing:
   :math:`\left(\rho e\right)^{n+1,\ast}=\left(\rho e\right)^{n+1,\ast} + \left(\frac{a^{n}+a^{n+1}}{2}\right) \left(\frac{1}{a^{n+1}}\right)^2 \Delta tI_{R_{\rho e}}`
	  
   :math:`\left(\rho E\right)^{n+1,\ast}=\left(\rho E\right)^{n+1,\ast} + \left(\frac{a^{n}+a^{n+1}}{2}\right) \left(\frac{1}{a^{n+1}}\right)^2 \Delta tI_{R_{\rho e}}`

   We store Ne and Temp based on eos\_ hc updates from :math:`(e)^{n+1}`

Step 8:
   Repeat step 3-7

Step 9:
   Compute :math:`{\phi}^{n+1}` and :math:`{\bf g}^{n+1}` using
   :math:`\rho^{n+1,*}` and :math:`\rho_{dm}^{n+1}`, where :math:`\rho_{dm}^{n+1}`
   is computed from the particles at :math:`{\bf x}_i^{n+1}`.

   Here we can use :math:`{\phi}^n` as an initial guess for :math:`{\phi}^{n+1}` in order to reduce the time
   spent in multigrid to reach the specified tolerance.

Step 10:
   Correct :math:`{\bf U}` with time-centered source terms, and replace :math:`e` by
   :math:`E - \frac{1}{2}U^2` as appropriate.

   We time-center the
   gravitational source terms only,

   .. math:: {\bf U}^{n+1} = {\bf U}^{n+1} + \frac{\Delta t}{2} (S_g^{n+1} - S_g^n)

Step 11:
   Interpolate :math:`{\bf g}^{n+1}` from the grid to the particle locations, then
   update the particle velocities, :math:`{\bf u}_i^{n+1}`

   .. math::

      \begin{aligned}
          {\bf u}_i^{n+1} &=& \frac{1}{a^{n+1}}
                          \left( \left( a^{{n+\frac{1}{2}}} {\bf u}^{{n+\frac{1}{2}}}_i \right)
                               + \frac{\Delta t}{2} \; {\bf g}^{n+1}_i \right)  \end{aligned}

Step \**:
   in post\_ timestep, do a reset and compute\_ new\_ temp after syncing the gravity sources
.. _Chap:NightlyTesting :

Nightly Regression Tests
========================

The following regression tests are run nightly with Nyx.   The plotfiles generated in each night's test 
are compared with the benchmark plotfiles using the AMReX :cpp:`fcompare` utility to compare the mesh data
and :cpp:`particle_compare` to compare the particle data.

The results of these tests can be found at https://ccse.lbl.gov/pub/RegressionTesting/Nyx.
These tests are also run on an NVIDIA GPU and those results can be found at https://ccse.lbl.gov/pub/GpuRegressionTesting/Nyx.

We have a number of compile-time options -- this chart summarizes which regression tests
test which compile-time options.  (Note that the USE_GRAV and AMREX_USE_PARTICLES macros
which are used in the code source are both set to TRUE in the build system.)

+---------------------------+----------+--------------+
| Test                      | NO_HYDRO | USE_HEATCOOL |
+===========================+==========+==============+
| AMR-density               | FALSE    |   TRUE       |
+---------------------------+----------+--------------+
| DR_restart                | FALSE    |   FALSE      |
+---------------------------+----------+--------------+
| DoubleRarefaction         | FALSE    |   FALSE      |
+---------------------------+----------+--------------+
| DrivenTurbulence          | FALSE    |   FALSE      |
+---------------------------+----------+--------------+
| DrivenTurbulence-OMP      | FALSE    |   FALSE      |
|  (CPU only)               |          |              |
+---------------------------+----------+--------------+
| LyA                       | FALSE    |   TRUE       |
+---------------------------+----------+--------------+
| LyA-adiabatic             | FALSE    |   FALSE      |
+---------------------------+----------+--------------+
| LyA-OMP                   | FALSE    |   TRUE       |
|  (CPU only)               |          |              |
+---------------------------+----------+--------------+
| LyA_Neutrinos             | FALSE    |   TRUE       |
|  (CPU only)               |          |              |
+---------------------------+----------+--------------+
| MiniSB                    | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| MiniSB-ppm                | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| MiniSB-ref                | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| MiniSB-ref-ppm            | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line                | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line-nbody          | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line-nomesh         | TRUE     |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line_restart        | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line_restart-nbody  | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-2line_restart-nomesh | TRUE     |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.nosub           | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.nosub-nbody     | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.nosub-nomesh    | TRUE     |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.sub             | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.sub-nbody       | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.sub-nbody       | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Part-mass.sub-nomesh      | TRUE     |  FALSE       |
+---------------------------+----------+--------------+
| Sedov                     | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Sedov-ppm                 | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Sod                       | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| Sod-ppm                   | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| StrongShockTube           | FALSE    |  FALSE       |
+---------------------------+----------+--------------+
| StrongShockTube-ppm       | FALSE    |  FALSE       |
+---------------------------+----------+--------------+

Below nx,ny,nz are the number of cells in each coordinate direction at the coarsest level;
Npa = number of particles, Np = number of MPI ranks / number of OMP threads per rank.
(If just a single number then pure MPI.)

+---------------------------+----------+----+--------+-----+----------------------+
| Test                      | nx ny nz | Nl | Npa    | Np  | What does this test? |
+===========================+==========+====+========+=====+======================+
| AMR-density               | 64 64 64 | 3  | 262144 | 4   | Start from chk00300  |
+---------------------------+----------+----+--------+-----+----------------------+
| DR_restart                | 32  4  4 | 3  | 0      | 2/2 | Hydro only - restart |
+---------------------------+----------+----+--------+-----+----------------------+
| DoubleRarefaction         | 32  4  4 | 3  | 0      | 2/2 | Hydro only           |
+---------------------------+----------+----+--------+-----+----------------------+
| DrivenTurbulence          | 32 32 32 | 1  | 0      | 4   | Turbulent forcing    |
+---------------------------+----------+----+--------+-----+----------------------+
| DrivenTurbulence-OMP      | 32 32 32 | 1  | 0      | 1/4 |  Turbulent forcing   |
|  (CPU only)               |          |    |        |     |  with OMP            |
+---------------------------+----------+----+--------+-----+----------------------+
| LyA                       | 32 32 32 | 1  | 32768  | 4   | Includes h/c         |
+---------------------------+----------+----+--------+-----+----------------------+
| Lya-OMP                   | 32 32 32 | 1  | 32768  | 1/4 | LyA with OMP         |
|  (CPU only)               |          |    |        |     |                      |
+---------------------------+----------+----+--------+-----+----------------------+
| LyA-adiabatic             | 32 32 32 | 1  | 32768  | 4   | No h/c               |
+---------------------------+----------+----+--------+-----+----------------------+
| LyA_Neutrinos             | 32 32 32 | 1  | 32768  | 1   | LyA with OMP         |
|  (CPU only)               |          |    |        |     |                      |  
+---------------------------+----------+----+--------+-----+----------------------+
| MiniSB                    | 32 32 32 | 1  | 2500   | 2   | Small version of SB  |
+---------------------------+----------+----+--------+-----+----------------------+
| MiniSB-ppm                | 32 32 32 | 1  | 2500   | 2   | Small version of SB  |
+---------------------------+----------+----+--------+-----+----------------------+
| MiniSB-ref                | 32 32 32 | 2  | 2500   | 2   | Small version of SB  |
|                           |          |    |        |     | with refinement      |
+---------------------------+----------+----+--------+-----+----------------------+
| MiniSB-ref-ppm            | 32 32 32 | 2  | 2500   | 2   | Small version of SB  |
|                           |          |    |        |     | with refinement      |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line                | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line-nbody          | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line-nomesh         | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line_restart        | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line_restart-nbody  | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-2line_restart-nomesh | 16 16 16 | 3  | 16     | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.nosub           | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.nosub-nbody     | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.nosub-nomesh    | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.sub             | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.sub-nbody       | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Part-mass.sub-nomesh      | 16 16 16 | 3  | 8      | 2/2 | Particle-only        |
+---------------------------+----------+----+--------+-----+----------------------+
| Sedov                     | 32 32 32 | 1  | 0      | 1   | Hydro only (PLM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| Sedov-ppm                 | 32 32 32 | 1  | 0      | 1   | Hydro only (PPM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| Sod                       | 4  32 4  | 3  | 0      | 1   | Hydro only (PLM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| Sod-ppm                   | 4  32 4  | 3  | 0      | 1/2 | Hydro only (PPM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| StrongShockTube           | 32 4  4  | 3  | 0      | 1/2 | Hydro only (PLM)     |
+---------------------------+----------+----+--------+-----+----------------------+
| StrongShockTube-ppm       | 32 4  4  | 3  | 0      | 1   | Hydro only (PPM)     |
+---------------------------+----------+----+--------+-----+----------------------+

*******
Preface
*******

Nyx is a adaptive mesh, N-body/hydro code for cosmological simulation on massively parallel
computers.  It couples the compressible hydrodynamic equations on a grid with a particle represenation
of dark matter.

The major capabilities:

  * 3-dimensional unsplit, 2nd-order hydrodynamics

  * adaptive mesh refinement with subcycling; jumps of 2x and 4x between levels

  * full Poisson gravity (with triply periodic boundary conditions)

  * hybrid parallelization strategy with MPI + X, where X = OpenMP on multicore architectures
and CUDA/HIP/DPC++ on hybrid CPU/GPU architectures.

Nyx uses an Eulerian grid for the hydrodynamics solver and incorporates adaptive mesh refinement (AMR).
Our approach to AMR uses a nested hierarchy of logically-rectangular grids with simultaneous 
refinement in both space and time, utilizing the AMReX library.

.. role:: cpp(code)
  :language: c++

******
Inputs
******
.. toctree::
   :maxdepth: 1

The Nyx executable reads run-time information from an “inputs” file which you put on the command line. 
This section describes the inputs which can be specified either in the inputs file or on the command line.
If a value is specified on the command line, that value will override a value specified in the inputs file.

Problem Geometry
================

List of Parameters
------------------

+--------------------------+-----------------+-----------------+-------------+
| Parameter                | Definition      | Acceptable      | Default     |
|                          |                 | Values          |             |
+==========================+=================+=================+=============+
| **geometry.prob_lo**     | physical        | Real            | must be set |
|                          | location of low |                 |             |
|                          | corner of the   |                 |             |
|                          | domain          |                 |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.prob_hi**     | physical        | Real            | must be set |
|                          | location of     |                 |             |
|                          | high corner of  |                 |             |
|                          | the domain      |                 |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.coord_sys**   | coordinate      | 0 = Cartesian,  | must be set |
|                          | system          | 1 = r-z, 2 =    |             |
|                          |                 | spherical       |             |
+--------------------------+-----------------+-----------------+-------------+
| **geometry.is_periodic** | is the domain   | 0 if false, 1   | 0 0 0       |
|                          | periodic in     | if true         |             |
|                          | this direction  |                 |             |
+--------------------------+-----------------+-----------------+-------------+

Examples of Usage
-----------------

-  **geometry.prob_lo** = 0 0 0
   defines the low corner of the domain at (0,0,0) in physical space.

-  **geometry.prob_hi** = 1.e8 2.e8 2.e8
   defines the high corner of the domain at (1.e8,2.e8,2.e8) in
   physical space.

-  **geometry.coord_sys** = 0
   defines the coordinate system as Cartesian

-  **geometry.is_periodic** = 0 1 0
   says the domain is periodic in the y-direction only.

Domain Boundary Conditions
==========================

.. _list-of-parameters-1:

List of Parameters
------------------

+---------------+---------------------------------+-------------------+-------------+
| Parameter     | Definition                      | Acceptable Values | Default     |
+===============+=================================+===================+=============+
| **nyx.lo_bc** | boundary type of each low face  | 0,1,2,3,4,5       | must be set |
+---------------+---------------------------------+-------------------+-------------+
| **nyx.hi_bc** | boundary type of each high face | 0,1,2,3,4,5       | must be set |
+---------------+---------------------------------+-------------------+-------------+

[Table:BC]

Notes
-----

Boundary types are:

======================= ================
0 – Interior / Periodic 3 – Symmetry      
1 – Inflow              4 – Slip Wall     
2 – Outflow             5 – No Slip Wall  
======================= ================

Note – **nyx.lo_bc** and **nyx.hi_bc** must be consistent with
**geometry.is_periodic** – if the domain is periodic in a particular
direction then the low and high bc’s must be set to 0 for that
direction.

.. _examples-of-usage-1:

Examples of Usage
-----------------

-  **nyx.lo_bc** = 1 4 0

-  **nyx.hi_bc** = 2 4 0

-  **geometry.is_periodic** = 0 0 1

would define a problem with inflow (1) in the low-x direction,
outflow(2) in the high-x direction, slip wall (4) on the low and high
y-faces, and periodic in the z-direction.

Resolution
==========

.. _list-of-parameters-2:

List of Parameters
------------------

+---------------------------+-----------------+-----------------+-------------+
| Parameter                 | Definition      | Acceptable      | Default     |
|                           |                 | Values          |             |
+===========================+=================+=================+=============+
| **amr.n_cell**            | number of cells | Integer > 0     | must be set |
|                           | in each         |                 |             |
|                           | direction at    |                 |             |
|                           | the coarsest    |                 |             |
|                           | level           |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.max_level**         | number of       | Integer >= 0    | must be set |
|                           | levels of       |                 |             |
|                           | refinement      |                 |             |
|                           | above the       |                 |             |
|                           | coarsest level  |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.ref_ratio**         | ratio of coarse | 2 or 4          | must be set |
|                           | to fine grid    |                 |             |
|                           | spacing between |                 |             |
|                           | subsequent      |                 |             |
|                           | levels          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_int**        | how often to    | Integer > 0     | must be set |
|                           | regrid          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_on_restart** | should we       | 0 or 1          | 0           |
|                           | regrid          |                 |             |
|                           | immediately     |                 |             |
|                           | after           |                 |             |
|                           | restarting      |                 |             |
+---------------------------+-----------------+-----------------+-------------+

[Table:ResInputs]

Note: if **amr.max_level** = 0 then you do not need to set
**amr.ref_ratio** or **amr.regrid_int**.

.. _examples-of-usage-2:

Examples of Usage
-----------------

-  **amr.n_cell** = 32 64 64

   would define the domain to have 32 cells in the x-direction, 64 cells
   in the y-direction, and 64 cells in the z-direction *at the coarsest
   level*. (If this line appears in a 2D inputs file then the final
   number will be ignored.)

-  | **amr.max_level** = 2
   | would allow a maximum of 2 refined levels in addition to the coarse
     level. Note that these additional levels will only be created only
     if the tagging criteria are such that cells are flagged as needing
     refinement. The number of refined levels in a calculation must be
     :math:`\leq` **amr.max_level**, but can change in time and need not
     always be equal to **amr.max_level**.

-  | **amr.ref_ratio** = 2 4
   | would set factor 2 refinement between levels 0 and 1, and factor 4
     refinement between levels 1 and 2. Note that you must have at least
     **amr.>ax_level** values of **amr.ref_ratio** (Additional values
     may appear in that line and they will be ignored).

-  | **amr.regrid_int** = 2 2
   | tells the code to regrid every 2 steps. Thus in this example, new
     level-1 grids will be created every 2 level-0 time steps, and new
     level-2 grids will be created every 2 level-1 time steps.

Regridding
==========

Overview
--------

The details of the regridding strategy are described in Section
`5.5 <#subsec:grid-generation>`__, but first we cover how the input
parameters can control the gridding.

As described later, the user defines Fortran subroutines which tag
individual cells at a given level if they need refinement. This list of
tagged cells is sent to a grid generation routine, which uses the
Berger–Rigoutsos algorithm to create rectangular grids that contain the
tagged cells.

.. _list-of-parameters-4:

List of Parameters
------------------

+----------------------------+----------------+----------------+----------------+
| Parameter                  | Definition     | Acceptable     | Default        |
|                            |                | Values         |                |
+============================+================+================+================+
| **amr.regrid_file**        | name of file   | text           | no file        |
|                            | from which to  |                |                |
|                            | read the grids |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.grid_eff**           | grid           | Real > 0, < 1  | 0.7            |
|                            | efficiency at  |                |                |
|                            | coarse level   |                |                |
|                            | at which grids |                |                |
|                            | are created    |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.n_error_buf**        | radius of      | Integer >= 0   | 1              |
|                            | additional     |                |                |
|                            | tagging around |                |                |
|                            | already tagged |                |                |
|                            | cells          |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.max_grid_size**      | maximum size   | Integer > 0    | 128 in 2D, 32  |
|                            | of a grid in   |                | in 3D          |
|                            | any direction  |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.max_grid_size**      | maximum size   | Integer        | 128 in 2D, 32  |
+----------------------------+----------------+----------------+----------------+
| **amr.blocking_factor**    | grid size must | Integer > 0    | 2              |
|                            | be a multiple  |                |                |
|                            | of this        |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.refine_grid_layout** | refine grids   | 0 if false, 1  | 1              |
|                            | more if # of   | if true        |                |
|                            | processors     |                |                |
|                            | :math:`>` # of |                |                |
|                            | grids          |                |                |
+----------------------------+----------------+----------------+----------------+

[Table:GriddingInputs]

.. _notes-2:

Notes
-----

-  **amr.n_error_buf**, **amr.max_grid_size** and
   **amr.blocking_factor** can be read in as a single value which is
   assigned to every level, or as multiple values, one for each level

-  **amr.max_grid_size** at every level must be even

-  **amr.blocking_factor** at every level must be a power of 2

-  the domain size **amr.n_cell** must be a multiple of
   **amr.blocking_factor** at level 0

-  **amr.max_grid_size** must be a multiple of **amr.blocking_factor**
   at every level

.. _examples-of-usage-3:

Examples of Usage
-----------------

-  | **amr.regrid_file** = *fixed_grids*
   | In this case the list of grids at each fine level are contained in
     the file *fixed_grids*, which will be read during the gridding
     procedure. These grids must not violate the **amr.max_grid_size**
     criterion. The rest of the gridding procedure described below will
     not occur if **amr.regrid_file** is set.

-  | **amr.grid_eff** = 0.9
   | During the grid creation process, at least 90% of the cells in each
     grid at the level at which the grid creation occurs must be tagged
     cells. Note that this is applied at the coarsened level at which
     the grids are actually made, and before **amr.max_grid_size** is
     imposed.

-  | **amr.max_grid_size** = 64
   | The final grids will be no longer than 64 cells on a side at every
     level.

-  | **amr.max_grid_size** = 64 32 16
   | The final grids will be no longer than 64 cells on a side at level
     0, 32 cells on a side at level 1, and 16 cells on a side at level
     2.

-  | **amr.blocking_factor** = 32
   | The dimensions of all the final grids will be multiples of 32 at
     all levels.

-  | **amr.blocking_factor** = 32 16 8
   | The dimensions of all the final grids will be multiples of 32 at
     level 0, multiples of 16 at level 1, and multiples of 8 at level 2.

   Having grids that are large enough to coarsen multiple levels in a
   V-cycle is essential for good multigrid performance in simulations
   that use self-gravity.

.. _subsec:grid-generation:

How Grids are Created
---------------------

The gridding algorithm proceeds in this order:

#. If at level 0, the domain is initially defined by **n_cell** as
   specified in the inputs file. If at level greater than 0, grids are
   created using the Berger–Rigoutsis clustering algorithm applied to
   the tagged cells, modified to ensure that the lengths of all new fine
   grids are divisible by **blocking_factor**.

#. Next, the grid list is chopped up if any grids have length longer
   than **max_grid_size**. Note that because **max_grid_size** is a
   multiple of **blocking_factor** (as long as **max_grid_size** is
   greater than **blocking_factor**), the **blocking_factor** criterion
   is still satisfied.

#. Next, if **refine_grid_layout** = 1 and there are more processors
   than grids at this level, then the grids at this level are further
   divided in order to ensure that no processor has fewer than one grid
   (at each level).

   -  if **max_grid_size** / 2 in the **BL_SPACEDIM** direction is a
      multiple of **blocking_factor**, then chop the grids in the
      **BL_SPACEDIM** direction so that none of the grids are longer in
      that direction than **max_grid_size / 2**

   -  If there are still fewer grids than processes, repeat the
      procedure in the **BL_SPACEDIM-1** direction, and again in the
      **BL_SPACEDIM-2** direction if necessary

   -  If after completing a sweep in all coordinate directions with
      **max_grid_size / 2**, there are still fewer grids than processes,
      repeat the steps above with **max_grid_size / 4**.

Simulation Time
===============

.. _list-of-parameters-5:

List of Parameters
------------------

+-----------------+---------------------------+--------------+---------+
| Parameter       | Definition                | Acceptable   | Default |
|                 |                           | Values       |         |
+=================+===========================+==============+=========+
| **max_step**    | maximum number of level 0 | Integer >= 0 | -1      |
|                 | time steps                |              |         |
+-----------------+---------------------------+--------------+---------+
| **stop_time**   | final simulation          | Real >= 0    | -1.0    |
|                 | time                      |              |         |
+-----------------+---------------------------+--------------+---------+
| **nyx.final_a** | final value of a          | Real > 0     | -1.0    |
+-----------------+---------------------------+--------------+---------+
| **nyx.final_z** | final value of z          | Real > 0     | -1.0    |
+-----------------+---------------------------+--------------+---------+

[Table:TimeInputs]

.. _notes-3:

Notes
-----

To control the number of time steps, you can limit by the maximum number
of level-0 time steps (**max_step**), or the final simulation time
(**stop_time**), or both. The code will stop at whichever criterion
comes first. Note that if the code reaches **stop_time** then the final
time step will be shortened so as to end exactly at **stop_time**, not
pass it.

If running in comoving coordinates you can also set a final value of
:math:`a` by setting **nyx.final_a**, or a final value of :math:`z` by
setting **nyx.final_z**. You may only specify one or the other of these.
Once this value of :math:`a` or :math:`z` is reached in a time step, the
code will stop at the end of this full coarse time step. (Note it does
not stop at :math:`a` (or :math:`z`) exactly equal to the final value,
rather it stops once :math:`a` is greater than (or :math:`z` is less
than) this value.)

.. _examples-of-usage-4:

Examples of Usage
-----------------

-  **max_step** = 1000

-  **stop_time** = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level-0 steps taken equals 1000, whichever comes first.

Time Step
=========

-  If **nyx.do_hydro**\ :math:`= 1`, then typically the code chooses a
   time step based on the CFL number (dt = cfl \* dx / max(u+c) ).

-  If **nyx.do_hydro**\ :math:`= 0` and we are running with dark matter
   particles, then we use a time step based on the velocity of the
   particles, i.e., we set :math:`\Delta t` so that the particle goes no
   further than :math:`f \Delta t` in a coordinate direction where
   :math:`0 \leq f \leq 1.` The value for :math:`f` is currently
   hard-wired in Particles.H, but it will become an inputs parameter.

.. _list-of-parameters-6:

List of Parameters
------------------

+---------------------+----------------+----------------+----------------+
| Parameter           | Definition     | Acceptable     | Default        |
|                     |                | Values         |                |
+=====================+================+================+================+
| **nyx.cfl**         | CFL number for | Real > 0 and   | 0.8            |
|                     | hydro          | <= 1           |                |
|                     |                |                |                |
|                     |                |                |                |
+---------------------+----------------+----------------+----------------+
| **particles.cfl**   | CFL number for | Real > 0 and   | 0.5            |
|                     | particles      | <= 1           |                |
|                     |                |                |                |
|                     |                |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.init_shrink** | factor by      | Real > 0 and   | 1.0            |
|                     | which to       | <= 1           |                |
|                     | shrink the     |                |                |
|                     | initial time   |                |                |
|                     | step           |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.change_max**  | factor by      | Real >= 1      | 1.1            |
|                     | which the time |                |                |
|                     | step can grow  |                |                |
|                     | in subsequent  |                |                |
|                     | steps          |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.fixed_dt**    | level-0 time   | Real > 0       | unused if not  |
|                     | step           |                | set            |
|                     | regardless of  |                |                |
|                     | cfl or other   |                |                |
|                     | settings       |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.initial_dt**  | initial        | Real > 0       | unused if not  |
|                     | level-0 time   |                | set            |
|                     | step           |                |                |
|                     | regardless of  |                |                |
|                     | other settings |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.dt_cutoff**   | time step      | Real > 0       | 0.0            |
|                     | below which    |                |                |
|                     | calculation    |                |                |
|                     | will abort     |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.dt_binpow**   | time step      | Real >=  0     | -1.0           |
|                     | chosen to be   |                |                |
|                     | a power of a   |                |                |
|                     | half times the |                |                |
|                     | comoving time  |                |                |
|                     | step           |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.relative_max**| max da/dt      | Real > 0       | 0.01           |
| **_change_a**       |                |                |                |
+---------------------+----------------+----------------+----------------+
| **nyx.absolute_max**| a_new-a_old    | Real > 0       | -1.0           |
| **_change_a**       |                |                |                |
+---------------------+----------------+----------------+----------------+

[Table:TimeStepInputs]

.. _examples-of-usage-5:

Examples of Usage
-----------------

-  | **nyx.cfl** = 0.9
   | defines the timestep as dt = cfl \* dx / umax_hydro.

-  | **particles.cfl** = 0.9
   | defines the timestep as dt = cfl \* dx / umax_particles where
     umax_particles is the maximum velocity of any particle in the
     domain.

-  | **nyx.init_shrink** = 0.01
   | sets the initial time step to 1% of what it would be otherwise.

-  | **nyx.change_max** = 1.1
   | allows the time step to increase by no more than 10% in this case.
     Note that the time step can shrink by any factor; this only
     controls the extent to which it can grow.

-  | **nyx.fixed_dt** = 1.e-4
   | sets the level-0 time step to be 1.e-4 for the entire simulation,
     ignoring the other timestep controls. Note that if
     **nyx.init_shrink** :math:`\neq 1` then the first time step will in
     fact be **nyx.init_shrink** \* **nyx.fixed_dt**.

-  | **nyx.initial_dt** = 1.e-4
   | sets the *initial* level-0 time step to be 1.e-4 regardless of
     **nyx.cfl** or **nyx.fixed_dt**. The time step can grow in
     subsequent steps by a factor of **nyx.change_max** each step.

-  | **nyx.dt_cutoff** = 1.e-20
   | tells the code to abort if the time step ever gets below 1.e-20.
     This is a safety mechanism so that if things go nuts you don’t burn
     through your entire computer allocation because you don’t realize
     the code is misbehaving.

-  | **nyx.dt_binpow** = 1.0
   | sets :math:`\mathrm{dt}=\left(\frac{1}{2}\right)^{n}\mathrm{dt}_{\mathrm{a}}|n:\mathrm{dt}_{\mathrm{cfl}}>\left(\frac{1}{2}\right)^{n}\mathrm{dt_{a}}`
     where :math:`\mathrm{dt}_{\mathrm{cfl}}` is determined by the more
     restrictive timestep of **nyx.cfl** and **particles.cfl**, and
     where :math:`\mathrm{dt}_{\mathrm{a}}` is determined by the
     **relative_max_change_a**, **absolute_max_change_a**, and the
     evolution of :math:`\frac{da}{dt}`

Subcycling
==========

 supports a number of different modes for subcycling in time.

-  If **amr.subcycling_mode**\ :math:`=`\ Auto (default), then the code
   will run with equal refinement in space and time. In other words, if
   level :math:`n+1` is a factor of 2 refinement above level :math:`n`,
   then :math:`n+1` will take 2 steps of half the duration for every
   level :math:`n` step.

-  If **amr.subcycling_mode**\ :math:`=`\ None, then the code will not
   refine in time. All levels will advance together with a timestep
   dictated by the level with the strictest :math:`dt`. Note that this
   is identical to the deprecated command **amr.nosub = 1**.

-  If **amr.subcycling_mode**\ :math:`=`\ Manual, then the code will
   subcycle according to the values supplied by
   **subcycling_iterations**.

.. _list-of-parameters-7:

List of Parameters
------------------

+----------------+----------------+----------------+----------------+
| Parameter      | Definition     | Acceptable     | Default        |
|                |                | Values         |                |
+================+================+================+================+
| **amr.sub      | How shall we   | Auto, None or  | Auto           |
| cycling_mode** | subcycle       | Manual         |                |
+----------------+----------------+----------------+----------------+
| **amr.subcycli | Number of      | 1 or           | must be set in |
| g_iterations** | cycles at each | ``ref_ratio``  | Manual mode    |
|                | level          |                |                |
+----------------+----------------+----------------+----------------+

.. _examples-of-usage-6:

Examples of Usage
-----------------

-  | **amr.subcycling_mode**\ :math:`=`\ Manual
   | Subcycle in manual mode with largest allowable timestep.

-  | **amr.subcycling_iterations** = 1 2 1 2
   | Take 1 level-0 timestep at a time (required). Take 2 level-1
     timesteps for each level-0 step, 1 timestep at level 2 for each
     level-1 step, and take 2 timesteps at level 3 for each level 2
     step.

-  | **amr.subcycling_iterations** = 2
   | Alternative form. Subcycle twice at every level (except level 0).

Restart Capability
==================

|  has a standard sort of checkpointing and restarting capability. In
  the inputs file, the following options control the generation of
  checkpoint files (which are really directories):

.. _list-of-parameters-8:

List of Parameters
------------------

+---------------------------------+----------------+----------------+----------------+
| Parameter                       | Definition     | Acceptable     | Default        |
|                                 |                | Values         |                |
+=================================+================+================+================+
| **amr.check_file**              | prefix for     | String         | “*chk*”        |
|                                 | restart files  |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_int**               | how often (by  | Integer        | -1             |
|                                 | level-0 time   | :math:`> 0`    |                |
|                                 | steps) to      |                |                |
|                                 | write restart  |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_per**               | how often (by  | Real           | -1.0           |
|                                 | simulation     | :math:`> 0`    |                |
|                                 | time) to write |                |                |
|                                 | restart files  |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.restart**                 | name of the    | String         | not used if    |
|                                 | file           |                | not set        |
|                                 | (directory)    |                |                |
|                                 | from which to  |                |                |
|                                 | restart        |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.checkpoint_files_output** | should we      | 0 or 1         | 1              |
|                                 | write          |                |                |
|                                 | checkpoint     |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.check_nfiles**            | how parallel   | Integer        | 64             |
|                                 | is the writing | :math:`\geq 1` |                |
|                                 | of the         |                |                |
|                                 | checkpoint     |                |                |
|                                 | files          |                |                |
+---------------------------------+----------------+----------------+----------------+
| **amr.checkpoint_on_restart**   | should we      | 0 or 1         | 0              |
|                                 | write a        |                |                |
|                                 | checkpoint     |                |                |
|                                 | immediately    |                |                |
|                                 | after          |                |                |
|                                 | restarting     |                |                |
+---------------------------------+----------------+----------------+----------------+

.. _notes-4:

Notes
-----

-  You should specify either **amr.check_int** or **amr.check_per**. Do
   not try to specify both.

-  Note that if **amr.check_per** is used then in order to hit that
   exact time the code may modify the time step slightly, which will
   change your results ever so slightly than if you didn’t set this
   flag.

-  Note that **amr.plotfile_on_restart** and
   **amr.checkpoint_on_restart** only take effect if
   **amr.regrid_on_restart** is in effect.

-  See the Software Section for more details on parallel I/O and the
   **amr.check_nfiles** parameter.

-  If you are doing a scaling study then set
   **amr.checkpoint_files_output** = 0 so you can test scaling of the
   algorithm without I/O.

.. _examples-of-usage-7:

Examples of Usage
-----------------

-  **amr.check_file** = *chk_run*

-  **amr.check_int** = 10

   means that restart files (really directories) starting with the
   prefix “*chk_run*” will be generated every 10 level-0 time steps. The
   directory names will be *chk_run00000*, *chk_run00010*,
   *chk_run00020*, etc.

If instead you specify

-  **amr.check_file** = *chk_run*

-  **amr.check_per** = 0.5

   then restart files (really directories) starting with the prefix
   “*chk_run*” will be generated every 0.1 units of simulation time. The
   directory names will be *chk_run00000*, *chk_run00043*,
   *chk_run00061*, etc, where :math:`t = 0.1` after 43 level-0 steps,
   :math:`t = 0.2` after 61 level-0 steps, etc.

To restart from *chk_run00061*,for example, then set

-  **amr.restart** = *chk_run00061*

.. _sec:PlotFiles:

Controlling PlotFile Generation
===============================

The main output from  is in the form of plotfiles (which are really
directories). The following options in the inputs file control the
generation of plotfiles

.. _list-of-parameters-9:

List of Parameters
------------------

+-----------------------------+------------------+------------------+---------+
| Parameter                   | Definition       | Acceptable       | Default |
|                             |                  | Values           |         |
+=============================+==================+==================+=========+
| **amr.plot_file**           | prefix for       | String           | “*plt*” |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_int**            | how often (by    | Integer          | -1      |
|                             | level-0 time     | :math:`> 0`      |         |
|                             | steps) to write  |                  |         |
|                             | plot files       |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_per**            | how often (by    | Real :math:`> 0` | -1.0    |
|                             | simulation time) |                  |         |
|                             | to write plot    |                  |         |
|                             | files            |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_vars**           | name of state    | ALL, NONE or     | ALL     |
|                             | variables to     | list             |         |
|                             | include in       |                  |         |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.derive_plot_vars**    | name of derived  | ALL, NONE or     | NONE    |
|                             | variables to     | list             |         |
|                             | include in       |                  |         |
|                             | plotfiles        |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_files_output**   | should we write  | 0 or 1           | 1       |
|                             | plot files       |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plotfile_on_restart** | should we write  | 0 or 1           | 0       |
|                             | a plotfile       |                  |         |
|                             | immediately      |                  |         |
|                             | after restarting |                  |         |
+-----------------------------+------------------+------------------+---------+
| **amr.plot_nfiles**         | how parallel is  | Integer          | 64      |
|                             | the writing of   | :math:`\geq 1`   |         |
|                             | the plotfiles    |                  |         |
+-----------------------------+------------------+------------------+---------+
| **nyx.plot_rank**           | should we plot   | True / False     | False   |
|                             | the processor ID |                  |         |
|                             | in the plotfiles |                  |         |
+-----------------------------+------------------+------------------+---------+
| **fab.format**              | Should we write  | NATIVE or IEEE32 | NATIVE  |
|                             | the plotfile in  |                  |         |
|                             | double or single |                  |         |
|                             | precision        |                  |         |
+-----------------------------+------------------+------------------+---------+

All the options for **amr.derive_plot_vars** are kept in ``derive_lst``
in ``Nyx_setup.cpp``. Feel free to look at it and see what’s there.

.. _notes-5:

Notes
-----

-  You should specify either **amr.plot_int** or **amr.plot_per**. Do
   not try to specify both.

-  Note that if **amr.plot_per** is used then in order to hit that exact
   time the code may modify the time step slightly, which will change
   your results ever so slightly than if you didn’t set this flag.

-  See the Software Section for more details on parallel I/O and the
   **amr.plot_nfiles** parameter.

-  If you are doing a scaling study then set **amr.plot_files_output** =
   0 so you can test scaling of the algorithm without I/O.

-  By default, plotfiles are written in double precision (NATIVE
   format). If you want to save space by writing them in single
   precision, set the fab.format flag to IEEE32.

.. _examples-of-usage-8:

Examples of Usage
-----------------

-  **amr.plot_file** = *plt_run*

-  **amr.plot_int** = 10

   means that plot files (really directories) starting with the prefix
   “*plt_run*” will be generated every 10 level-0 time steps. The
   directory names will be *plt_run00000*, *plt_run00010*,
   *plt_run00020*, etc.

If instead you specify

-  **amr.plot_file** = *plt_run*

-  **amr.plot_per** = 0.5

   then restart files (really directories) starting with the prefix
   “plt_run” will be generated every 0.1 units of simulation time. The
   directory names will be *plt_run00000*, *plt_run00043*,
   *plt_run00061*, etc, where :math:`t = 0.1` after 43 level-0 steps,
   :math:`t = 0.2` after 61 level-0 steps, etc.

Plotfile Variables
------------------

Native variables
^^^^^^^^^^^^^^^^

These variables come directly from the ``StateData``, either the
``State_Type`` (for the hydrodynamic variables), ``DiagEOS_Type``
(for the nuclear energy generation quantities). ``PhiGrav_Type`` and
``Gravity_Type`` (for the gravity quantities)

+-----------------------------------+---------------------------------------------------+--------------------------------------+
| variable name                     | description                                       | units                                |
+===================================+===================================================+======================================+
| ``density``                       | Baryonic mass density, :math:`\rho`               | M\ :math:`_\odot` / Mpc\ :math:`^3`  |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``xmom``                          | x-momentum, :math:`(\rho u)`                      | :math:`{\rm g~km^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``ymom``                          | y-momentum, :math:`(\rho v)`                      | :math:`{\rm g~km^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``zmom``                          | z-momentum, :math:`(\rho w)`                      | :math:`{\rm g~km^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_E``                         | Total energy density                              | :math:`{\rm erg~km^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_e``                         | Internal energy density                           | :math:`{\rm erg~km^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Temp``                          | Temperature                                       | :math:`{\rm K}`                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Ne``                            | Number density of electrons                       | dimensionless                        |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_X``                         | Mass density of species X (only valid for non-    | dimensionless                        |
| (where X is H or He, the species  | constant species)                                 |                                      |
| defined in the network)           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiGrav``                       | Gravitational potential                           | :math:`{\rm erg~g^{-1}}`             |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``grav_x``, ``grav_y``,           | Gravitational acceleration                        | :math:`{\rm km~s^{-2}}`              |
| ``grav_z``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+

Derived variables
^^^^^^^^^^^^^^^^^

+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| variable name                     | description                                       | derive routine              | units                                   |
+===================================+===================================================+=============================+=========================================+
| ``divu``                          | :math:`\nabla \cdot \ub`                          | ``derdivu``                 | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_e``                        | Specific internal energy computed from the        | ``dereint2``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | conserved :math:`(\rho e)` state variable as      |                             |                                         |
|                                   | :math:`e = (\rho e)/\rho`                         |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_E``                        | Specific internal energy computed from the        | ``dereint1``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | total energy and momentum conserved state as      |                             |                                         |
|                                   | :math:`e=[(\rho E)-\frac{1}{2}(\rho \ub^2)]/\rho` |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``kineng``                        | Kinetic energy density,                           | ``derkineng``               | :math:`{\rm erg~km^{-3}}`               |
|                                   | :math:`K = \frac{1}{2} |(\rho \ub)|^2`            |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``logden``                        | :math:`\log_{10} \rho`                            | ``derlogden``               | M\ :math:`_\odot` / Mpc\ :math:`^3`     |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``MachNumber``                    | Fluid Mach number, :math:`|\ub|/c_s`              | ``dermachnumber``           | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``maggrav``                       | Gravitational acceleration magnitude              | ``dermaggrav``              | :math:`{\rm km~s^{-2}}`                 |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magmom``                        | Momentum density magnitude,                       | ``dermagmom``               | :math:`{\rm g~km^{-2}~s^{-1}}`          |
|                                   | :math:`|\rho \ub|`                                |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvel``                        | Velocity magnitude, :math:`|\ub|`                 | ``dermagvel``               | :math:`\mathrm{km~s^{-1}}`              |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvort``                       | Vorticity magnitude, :math:`|\nabla\times\ub|`    | ``dermagvort``              | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``pressure``                      | Total pressure, including ions and electrons      | ``derpres``                 | :math:`{\rm dyn~km^{-2}}`               |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``soundspeed``                    | Sound speed                                       | ``dersoundspeed``           | :math:`\mathrm{km~s^{-1}}`              |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``H`` or ``He``                   | Mass fraction of species H or He                  | ``derspec``                 | --                                      |
|                                   |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``x_velocity``,                   | Fluid velocity,                                   | ``dervel``                  | :math:`\mathrm{km~s^{-1}}`              |
| ``y_velocity``,                   | :math:`\ub = (\rho \ub)/\rho`                     |                             |                                         |
| ``z_velocity``                    |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
	 
Screen Output
=============

.. _list-of-parameters-10:

List of Parameters
------------------

+----------------------------+------------------+----------------+----------------+
| Parameter                  | Definition       | Acceptable     | Default        |
|                            |                  | Values         |                |
+============================+==================+================+================+
| **amr.v**                  | verbosity of     | 0 or 1         | 0              |
|                            | Amr.cpp          |                |                |
+----------------------------+------------------+----------------+----------------+
| **nyx.v**                  | verbosity of     | 0 or 1         | 0              |
|                            | Nyx.cpp          |                |                |
+----------------------------+------------------+----------------+----------------+
| **gravity.v**              | verbosity of     | 0 or 1         | 0              |
|                            | Gravity.cpp      |                |                |
+----------------------------+------------------+----------------+----------------+
| **mg.v**                   | verbosity of     | 0,1,2,3,4      | 0              |
|                            | multigrid        |                |                |
|                            | solver (for      |                |                |
|                            | gravity)         |                |                |
+----------------------------+------------------+----------------+----------------+
| **particles.v**            | verbosity of     | 0,1,2,3,4      | 0              |
|                            | particle-related |                |                |
|                            | processes        |                |                |
+----------------------------+------------------+----------------+----------------+
| **amr.grid_log**           | name of the      | String         | not used if    |
|                            | file to which    |                | not set        |
|                            | the grids are    |                |                |
|                            | written          |                |                |
+----------------------------+------------------+----------------+----------------+
| **amr.run_log**            | name of the      | String         | not used if    |
|                            | file to which    |                | not set        |
|                            | certain output   |                |                |
|                            | is written       |                |                |
+----------------------------+------------------+----------------+----------------+
| **amr.run_log_terse**      | name of the      | String         | not used if    |
|                            | file to which    |                | not set        |
|                            | certain          |                |                |
|                            | (terser)         |                |                |
|                            | output is        |                |                |
|                            | written          |                |                |
+----------------------------+------------------+----------------+----------------+
| **amr.sum_interval**       | if               |                |                |
|                            | :math:`> 0,`     |                |                |
|                            | how often (in    |                |                |
|                            | level-0 time     |                |                |
|                            | steps)           |                |                |
+----------------------------+------------------+----------------+----------------+
|                            | to compute and   | Integer        | -1             |
|                            | print integral   |                |                |
|                            | quantities       |                |                |
+----------------------------+------------------+----------------+----------------+

.. _examples-of-usage-9:

Examples of Usage
-----------------

-  | **amr.grid_log** = *grdlog*
   | Every time the code regrids it prints a list of grids at all
     relevant levels. Here the code will write these grids lists into
     the file *grdlog*.

-  | **amr.run_log** = *runlog*
   | Every time step the code prints certain statements to the screen
     (if **amr.v** = 1), such as
   | STEP = 1 TIME = 1.91717746 DT = 1.91717746
   | PLOTFILE: file = plt00001
   | Here these statements will be written into *runlog* as well.

-  | **amr.run_log_terse** = *runlogterse*
   | This file, *runlogterse*, differs from *runlog* in that it only
     contains lines of the form
   | 10 0.2 0.005
   | in which “10” is the number of steps taken, “0.2” is the simulation
     time, and “0.005” is the level-0 time step. This file can be
     plotted very easily to monitor the time step.

-  | **nyx.sum_interval** = 2
   | if **nyx.sum_interval** :math:`> 0` then the code computes and
     prints certain integral quantities, such as total mass, momentum
     and energy in the domain every **nyx.sum_interval** level-0 steps.
     In this example the code will print these quantities every two
     coarse time steps. The print statements have the form
   | TIME= 1.91717746 MASS= 1.792410279e+34
   | for example. If this line is commented out then it will not compute
     and print these quanitities.

Gravity
=======

.. _list-of-parameters-11:

List of Parameters
------------------

+--------------------------+------------------+----------------+------------------+
| Parameter                | Definition       | Acceptable     | Default          |
|                          |                  | Values         |                  |
+==========================+==================+================+==================+
| **nyx.do_grav**          | Include          | 0 if false     | must be set      |
|                          | gravity as a     | 1 if true      |                  |
|                          | forcing term     |                |                  |
+--------------------------+------------------+----------------+------------------+
| **gravity.no_sync**      | whether to       | 0 if false     |  0               |
|                          | perform the      | 1 if true      |                  |
|                          | “sync solve”     |                |                  |
+--------------------------+------------------+----------------+------------------+
| **gravity.no_composite** | whether to       | 0 if false     |  0               |
|                          | perform a        | 1 if true      |                  |
|                          | composite        |                |                  |
|                          | solve            |                |                  |
+--------------------------+------------------+----------------+------------------+

.. _notes-6:

Notes
-----

-  To include gravity you must set **nyx.do_grav** = 1 in the inputs file

Physics
=======

.. _list-of-parameters-12:

List of Parameters
------------------

+----------------------------------+------------------+-----------------+-------------+
| Parameter                        | Definition       | Acceptable      | Default     |
|                                  |                  | Values          |             |
+==================================+==================+=================+=============+
| **nyx.do_hydro**                 | Time-advance     | 0 if false, 1   | must be set |
|                                  | the fluid        | if true         |             |
|                                  | dynamical        |                 |             |
|                                  | equations        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.ppm_type**                 | Use PPM or       | 0 for PLM       | 1 (PPM)     |
|                                  | PLM for hydro    | 1 for PPM       |             |
|                                  | advance          |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.enforce_min_density_type** | how to enforce   | "floor"         | "floor"     |
|                                  | rho greater than | "cons"          |             |
|                                  | small_dens       |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.strang_split**             | Use strang       | 0 if false, 1   | 1           |
|                                  | splitting        | if true         |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.sdc_split**                | Use sdc          | 0 if false, 1   | 0           |
|                                  | splitting        | if true         |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.strang_grown_box**         | Use growntilebox | 0 if false, 1   | 1           |
|                                  | to avoid comms   | if true         |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.add_ext_src**              | Include          | 0 if false, 1   | 0           |
|                                  | additional       | if true         |             |
|                                  | user-specified   |                 |             |
|                                  | source term      |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.nghost_state**             | Set number of    | {1,2,3,4}       | 1           |
|                                  | ghost cells for  |                 |             |
|                                  | state variables  |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.use_const_species**        | If 1 then read   | 0 or 1          | 0           |
|                                  | h_species and    |                 |             |
|                                  | he_species       |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.h_species**                | Concentration    | 0 :math:`<` X   | 0           |
|                                  | of H             | :math:`<` 1     |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.he_species**               | Concentration    | 0 :math:`<` X   | 0           |
|                                  | of He            | :math:`<` 1     |             |
+----------------------------------+------------------+-----------------+-------------+


Cosmology
=========

.. _list-of-parameters-13:

List of Parameters
------------------

+----------------------------------+--------------------+-----------------+-------------+
| Parameter                        | Definition         | Acceptable      | Default     |
|                                  |                    | Values          |             |
+==================================+====================+=================+=============+
| **nyx.comoving_OmM**             | Relative (total)   |  0 :math:`<` X  | must be set |
|                                  | mass density       |  :math:`<` 1    |             |
+----------------------------------+--------------------+-----------------+-------------+
| **nyx.comoving_OmB**             | Relative baryon    |  0 :math:`<` X  | must be set |
|                                  | density            |  :math:`<` 1    |             |
+----------------------------------+--------------------+-----------------+-------------+
| **nyx.comoving_OmR**             | Relative           |  0 :math:`<` X  | must be set |
|                                  | radiation density  |  :math:`<` 1    |             |
+----------------------------------+--------------------+-----------------+-------------+
| **nyx.comoving_h**               | Dimensionless      |  0 :math:`<` X  | must be set |
|                                  | Hubble parameter   |  :math:`<` 1    |             |
+----------------------------------+--------------------+-----------------+-------------+
| **nyx.gamma**                    | Dimensionless      |  0 :math:`<` X  | :math:`5/3` |
|                                  | factor relating    |  :math:`<` 2    |             |
|                                  | :math:`p, \rho, e` |                 |             |
+----------------------------------+--------------------+-----------------+-------------+

Examples of Usage
-----------------

-  | **nyx.gamma** This changes :math:`\gamma` in the :math:`\gamma` law gas: :math:`p = (\gamma - 1) \rho e.`

Reionization models
===================

.. _list-of-parameters-14:

List of Parameters
------------------

+----------------------------------+------------------+-----------------+-------------+
| Parameter                        | Definition       | Acceptable      | Default     |
|                                  |                  | Values          |             |
+==================================+==================+=================+=============+
| **uvb_rates_file**               | Name of the UVB  |  string         | must be set |
|                                  | (TREECOOL) file  |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **uvb_density_A**                | Density-dependent|  real           | 1.0         |
|                                  | heating          |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **uvb_density_B**                | Density dependent|  real           | 0.0         |
|                                  | heating          |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **reionization_zHI_flash**       | Redshift of      |  0 :math:`<` X  | -1.0        |
|                                  | "flash" H reion. |  or -1 if off   |             |
+----------------------------------+------------------+-----------------+-------------+
| **reionization_zHeII_flash**     | Redshift of      |  0 :math:`<` X  | -1.0        |
|                                  | "flash" He reion.|  of -1 if off   |             |
+----------------------------------+------------------+-----------------+-------------+
| **inhomo_reion**                 | Inhomogeneous    |  0 or 1         | 0           |
|                                  | reionization     |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **inhomo_zhi_file**              | File with        |  string         | must be set |
|                                  | reionization map |                 | (if used)   |
+----------------------------------+------------------+-----------------+-------------+
| **inhomo_grid**                  | Size of the      |  integer        | must be set |
|                                  | reionization grid|                 | (if used)   |
+----------------------------------+------------------+-----------------+-------------+
| **reionization_T_zHI**           | H reionization   |  real           | 2.0e4       |
|                                  | heat input       |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **reionization_T_zHeII**         | He reionization  |  real           | 1.5e4       |
|                                  | heat input       |                 |             |
+----------------------------------+------------------+-----------------+-------------+



Multigrid Inputs
================

The following inputs can be set directly in the AMReX solver classes but we set them via the Nyx gravity routines.

These must be preceded by "gravity" in the inputs file:

+----------------------+---------------------------------------------------+-----------+--------------+
|                      | Description                                       | Type      | Default      |
+======================+===================================================+===========+==============+
| v                    |  Verbosity of Gravity class                       |  Int      |   0          |
+----------------------+---------------------------------------------------+-----------+--------------+
| ml_tol               |  Relative tolerance for multilevel solves         |  Real     |   1.e-12     |
+----------------------+---------------------------------------------------+-----------+--------------+
| sl_tol               |  Relative tolerance for single-level solves       |  Real     |   1.e-12     |
+----------------------+---------------------------------------------------+-----------+--------------+
| delta_tol            |  Relative tolerance for synchronization solves    |  Real     |   1.e-12     |
+----------------------+---------------------------------------------------+-----------+--------------+
| mlmg_agglomeration   |  Should we agglomerate deep in the V-cycle        |  Int      |   1          |
+----------------------+---------------------------------------------------+-----------+--------------+
| mlmg_consolidation   |  Should we consolidate deep in the V-cycle        |  Int      |   1          |
+----------------------+---------------------------------------------------+-----------+--------------+
| dirichlet_bcs        |  Should we use homogeneous Dirichlet bcs in the   |  Int      |   0          |
|                      |  gravity solves (used for testing only)           |           |              |
+----------------------+---------------------------------------------------+-----------+--------------+

These must be preceded by "mg" in the inputs file:

+----------------------+-----------------------------------------------------+-------------+--------------+
|                      | Description                                         |  Type       | Default      |
+======================+=====================================================+=============+==============+
| v                    |  Verbosity of multigrid solver                      |  Int        |   0          |
+----------------------+-----------------------------------------------------+-------------+--------------+
| bottom_solver        |  What is the bottom solver?                         |  String     |   "bicg"     |
|                      |  Options include "bicg", "smoother", "hypre", etc   |             |              |
+----------------------+-----------------------------------------------------+-------------+--------------+
| max_fmg_iter         |  Maximum number of F-cycles to do before            |  Int        |   0          |
|                      |  continuing with V-cycles in a multigrid solve      |             |              |   
+----------------------+-----------------------------------------------------+-------------+--------------+

There are a number of additional inputs that can be used to control the multigrid solver.  

See the `AMReX Multigrid documentation`_ for more details.

.. _AMReX Multigrid documentation: https://amrex-codes.github.io/amrex/docs_html/LinearSolvers_Chapter.html

Memory Optimization
===================

.. _list-of-parameters-15:

List of Parameters
------------------

+----------------------------------+------------------+-----------------+-------------+
| Parameter                        | Definition       | Acceptable      | Default     |
|                                  |                  | Values          |             |
+==================================+==================+=================+=============+
| **nyx.shrink_to_fit**            | Shrink Particle  | 0 if false, 1   |             |
|                                  | vector to save   | if true         |             |
|                                  | memory           |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.minimize_memory**          | Use less         | 0 if false, 1   |             |
|                                  | temporary scratch| if true         |             |
|                                  | memory in hydro  |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_int**         | How often to     | Int < 0 if never| -1          |
|                                  | load-balance     | Int > 0         |             |
|                                  | particles        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_start_z**     | Redshift to start| Real > 0        | 7.0         |
|                                  | load-balancing   |                 |             |
|                                  | particles        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_wgt_stategy** | Weight strategy  | {0, 1, 2}       | 0           |
|                                  | to load-balance  |                 |             |
|                                  | particles        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_wgt_nmax**    | Max ranks to     | 0 < Int < Ranks | -1          |
|                                  | load-balance     |                 |             |
|                                  | particles        |                 |             |
+----------------------------------+------------------+-----------------+-------------+
| **nyx.load_balance_stategy**     | Dmap strategy    | {KNAPSACK,      | SFC         |
|                                  | type for particle|  SFC,           |             |
|                                  | load-balancing   |  ROUNDROBIN     |             |
+----------------------------------+------------------+-----------------+-------------+
Initial conditions
===================

There are a couple main ways in which initial conditions can be set in Nyx:
using an ASCII particle file, using binary particle file(s), using
a uniform particle setup, using a random particle setup, using a binary mesh file setup, or
using an analytic setup.
As said in the Units section of the documentation, the units are: Mpc, M\ :math:`_\odot`, and km/s,
and particle velocities should be peculiar proper velocities.


Start from an ASCII file
------------------------

To enable this option, set::
  
  nyx.particle_init_type = AsciiFile
  nyx.ascii_particle_file = *particle_file*

Here *particle_file* is the user-specified name of the file. The first line in this file is
(long) assumed to contain the number of particles. Each line after that contains

x y z mass vx vy vz



Start from a binary file
------------------------

To enable this option, set::

  nyx.particle_init_type = BinaryFile
  nyx.binary_particle_file = *particle_file*
  
With binary file, the header should have 3 numbers:
(long) NP, which is the total number of particles in the file
followed by the (int) DM=3 (number of dimensions), and (int) NX=4 (number of "extra" fields).
Following this small header, 7 float numbers should be listed for each particle as before:
x y z mass vx vy vz.

The main difference between the ASCII and binary format thus amounts to a different header.
Here is an example C++ code which writes a Nyx-readable binary file::

      fwrite(&npart, sizeof(long), 1, outFile);
      fwrite(&DM, sizeof(int), 1, outFile);
      fwrite(&NX, sizeof(int), 1, outFile);
      for(i=0; i<Npart; i++) {
         fwrite(&x[i], sizeof(float), 1, outFile);
         fwrite(&y[i], sizeof(float), 1, outFile);
         fwrite(&z[i], sizeof(float), 1, outFile);
         fwrite(&mass[i], sizeof(float), 1, outFile);
         fwrite(&vx[i], sizeof(float), 1, outFile);
         fwrite(&vy[i], sizeof(float), 1, outFile);
         fwrite(&vz[i], sizeof(float), 1, outFile);
      }


Start from a binary "meta" file
-------------------------------

This option allows you to read particles from a series of files rather than
just a single file. This is very convenient for large simulations.
To enable this option, set::

  nyx.particle_init_type = BinaryMetaFile
  nyx.binary_particle_file =*particle file*

In this case the *particle_file* you specify is an ASCII file specifying a
list of file names with full paths. Each of the files in this list is assumed
to be binary and is read sequentially (individual files are read in parallel) in
the order listed.

Since individual files are read sequentially, more particles should be read before
redistributing across MPI ranks. This is set by optimizing the maximum number of
readers and increasing the number of particles per read::

  amr.nreaders
  amr.nparts_per_read


Start from a plotfile or checkpoint
-----------------------------------

To enable this option, set::

  nyx.particle_init_type = Restart
  nyx.restart_particle_file = *plot_file*

In this case the *plot_file* should contain particles in directory *DM*. Testing of this
functionality is mainly for the current default *Version_Two_Dot_Zero_single*.


Reading SPH particles
---------------------

The above initialization from a single particle specie assumes that baryons trace dark matter.
Baryon density and velocity is set by CIC-mapping particles onto the Eulerian grid.
Alternatively, one can initialize the baryonic gas from the SPH
particles. To enable this option, you must set::

    nyx.do_santa_barbara = 1
    nyx.init_with_sph_particles =1

The SPH-type particles can then be read in by setting
where *sph_particle_file* is the user-specified name of the file
containing the SPH particles. The type of *sph_particle_file*
must be the same (Ascii, Binary or BinaryMeta) as the dark matter particle
file as specified by
The SPH particles will be discarded by the code once the grid data has been initialized.


Initial conditions for testing purposes
---------------------------------------

The following are used for code testing purposes and will not result in a meaningful cosmological simulation.


Random placement
^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.particle_init_type = Random
  
There are then a number of parameters to set, for example::
  
  nyx.particle_initrandom_count = 100000
  nyx.particle_initrandom_mass_total = 100000
  nyx.particle_initrandom_iseed = 15
  nyx.fix_random_seed = 0


Random placement (1 particle per grid cell)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.particle_init_type = RandomPerCell
  
Then only set the mass per particle::

  nyx.particle_initrandom_mass = 1

Note to increase the number of cells and keep the problem domain size 
and total mass fixed, the mass per particle must decrease proportionally.
An alternative is to set the total mass of all particles in the simulation.
This will cause Nyx to scale the mass per particle to fit the number of cells.
The following will all have the same total mass::

  nyx.particle_initrandom_mass = 1
  amr.n_cell = 64 64 64

  nyx.particle_initrandom_mass = 1
  nyx.particle_initrandom_mass_total = 262144
  amr.n_cell = 64 64 64

  nyx.particle_initrandom_mass = -1
  nyx.particle_initrandom_mass_total = 262144
  amr.n_cell = 64 64 64

  nyx.particle_initrandom_mass = 0.125
  amr.n_cell = 128 128 128

  nyx.particle_initrandom_mass = 0.125
  nyx.particle_initrandom_mass_total = 262144
  amr.n_cell = 128 128 128

  nyx.particle_initrandom_mass = -1
  nyx.particle_initrandom_mass_total = 262144
  amr.n_cell = 128 128 128


Uniform placement
^^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.particle_init_type = OnePerCell
  
There are then a number of parameters to set, for example::
  
  nyx.particle_inituniform_mass = 1
  nyx.particle_inituniform_vx = -1
  nyx.particle_inituniform_vy = 1
  nyx.particle_inituniform_vz = 1

Initial Multifab-based setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.particle_init_type = Cosmological
  nyx.do_readinics = 1

Then set the directory name of the MultiFab to restart the state variables from::

  nyx.readin_ics_fname = "mf"
  
Initial Analytic Problem Setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To enable this option, set::

  nyx.do_santa_barbara = 0
  nyx.init_sb_vels = 0

For different executable directories, the ``Prob.cpp`` setup can be further customised
with ``prob.`` input flags. For the HydroTests directory, ``prob.prob_type=0`` corresponds
to Sod, StrongShockTube and DoubleRarefaction type tests, and ``prob.prob_type!=0``
corresponds to the Sedov type tests.
.. Nyx documentation master file, created by
   sphinx-quickstart on Fri Dec  7 12:26:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Nyx's documentation!
===============================

.. toctree::
   :maxdepth: 1
   :caption: Nyx basics

   NyxPreface
   NyxCitation
   NyxGettingStarted
   NyxUnits
   NyxInputs
   ICs
   Refinement
   ManagingGridHierarchy_Chapter
   NyxGravity
   NyxEquations
   NyxSplitting
   NyxHeatCool
   NyxForcing
   Particles
   InSitu
   PostProcessing
   NightlyTests
   Debugging
   Development


.. toctree::      
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
******************
Stochastic Forcing
******************

There is an option to apply a stochastic force field. 

See Nyx/Exec/DrivenTurbulence for an example; note that ::

    nyx.do_forcing = 1

must be set in the inputs file.

The external forcing term in the momentum equation 
(`[eq:momt] <#eq:momt>`__) is then given by

  .. math:: {\bf S}_{\rho \Ub} = \rho_b \mathbf{f}

  where the acceleration field :math:`\mathbf{f}(\mathbf{x},t)` is
  computed as inverse Fourier transform of the forcing spectrum
  :math:`\widehat{\mathbf{f}}(\mathbf{k},t`). The time evolution of each
  wave mode is given by an Ornstein-Uhlenbeck process (see
  :raw-latex:`\cite{SchmHille06,Schmidt14}` for details). Since the real
  space forcing acts on large scales :math:`L`, non-zero modes are
  confined to a narrow window of small wave numbers with a prescribed
  shape (the forcing profile). The resulting flow reaches a
  statistically stationary and isotropic state with a root-mean-square
  velocity of the order :math:`V=L/T`, where the integral time scale
  :math:`T` (also known as large-eddy turn-over time) is usually set
  equal to the autocorrelation time of the forcing. It is possible to
  vary the force field from solenoidal (divergence-free) if the weight
  parameter :math:`\zeta=1` to dilational (rotation-free) if
  :math:`\zeta=0`.

To maintain a nearly constant root-mean-square Mach number, a simple
model for radiative heating and cooling around a given equilibrium
temperature :math:`T_0` is applied in the energy
equation (`[eq:energy] <#eq:energy>`__):

.. math:: S_{\rho E} = S_{\rho e} + \Ub \cdot {\bf S}_{\rho \Ub} = -\frac{\alpha k_{\rm B}(T-T_0)}{\mu m_{\rm H}(\gamma-1)} + \rho_b\Ub\cdot\mathbf{f}

The parameters :math:`T_0` and :math:`\alpha` correspond to temp0 and
alpha, respectively, in the probin file (along with rho0 for the mean
density, which is unity by default). While the gas is adiabatic for
:math:`\alpha=0`, it becomes nearly isothermal if the cooling time scale
given by :math:`1/\alpha` is chosen sufficiently short compared to
:math:`T`. For performance reasons, a constant composition
(corresponding to constant molecular weight :math:`\mu`) is assumed.

List of Parameters
==================

+-------------------------+--------------------+-----------------+-------------+
| Parameter               | Definition         | Acceptable      | Default     |
|                         |                    | Values          |             |
+=========================+====================+=================+=============+
| **forcing.seed**        | seed of the        | Integer         | 27011974    |
|                         | random number      | :math:`>0`      |             |
|                         | generator          |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.profile**     | shape of           | 1 (plane), 2    | 3           |
|                         | forcing            | (band), 3       |             |
|                         | spectrum           | (parabolic)     |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.alpha**       | ratio of domain    | Integer         | 2 2 2       |
|                         | size :math:`X`     | :math:`>0`      |             |
|                         | to integral        |                 |             |
|                         | length             |                 |             |
|                         | :math:`L=X/\alpha` |                 |             |
|                         |                    |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.band_width**  | band width of      | Real            | 1.0 1.0 1.0 |
|                         | the forcing        | :math:`\ge 0`   |             |
|                         | spectrum           | and             |             |
|                         | relative to        | :math:`\le 1`   |             |
|                         | alpha              |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.intgr_vel**   | characteristic     | Real            | must be set |
|                         | velocity           | :math:`> 0`     |             |
|                         | :math:`V`          |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.auto_corrl**  | autocorrelation    | Real            | 1.0 1.0 1.0 |
|                         | time in units      | :math:`> 0`     |             |
|                         | of                 |                 |             |
|                         | :math:`T=L/V`      |                 |             |
+-------------------------+--------------------+-----------------+-------------+
| **forcing.soln_weight** | weight             | Real            | 1.0         |
|                         | :math:`\zeta`      | :math:`\ge 0`   |             |
|                         | of solenoidal      | and             |             |
|                         | relative to        | :math:`\le 1`   |             |
|                         | dilatational       |                 |             |
|                         | modes              |                 |             |
+-------------------------+--------------------+-----------------+-------------+

Triples for forcing.alpha, forcing.band_width, forcing.intgr_vel, and
forcing.auto_corrl correspond to the three spatial dimensions.
.. _Chap:GettingStartedNew:

Getting Started
===============

The Nyx source code currently lives in a
`github repository <https://github.com/AMReX-Astro/Nyx.git>`_ that
can be cloned by using git:

.. code-block:: shell

   > git clone https://github.com/AMReX-Astro/Nyx.git --recursive

.. note::

   We recommend Git version 1.7.x or higher.

Once you have obtained the source code, the following sections describe the
source code contents, compiling, running a simple simulation, and visualizing
the simulations results.

.. toctree::
   :maxdepth: 1

   Source directory overview <getting_started/Structure>
   Building Nyx with GNU Make <getting_started/BuildingGMake>
   Building Nyx with CMake <getting_started/BuildingCMake>
   Running the code <getting_started/RunningTheCode>
   Compiling Nyx with SUNDIALS 5 <getting_started/NyxSundials>

Additional details about compiling Nyx in different environments can be found
in the `AMReX build documentation <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX_Chapter.html>`_ .
Containerized Ubuntu builds with CMake are run as 
`GitHub Actions <https://github.com/AMReX-Astro/Nyx/actions/workflows/linux.yml>`_
whose workflow files are stored in Nyx/.github/workflows .
**********
Citation
**********

If you use Nyx, we appreciate you citing the main code paper:
  * A. S. Almgren, J. B. Bell, M.J. Lijewski, Z. Lukic, E. Van Andel, **"Nyx: A Massively Parallel AMR Code for Computational Cosmology"** Astrophysical Journal, **765**, 39, 2013. [[Journal]](https://iopscience.iop.org/article/10.1088/0004-637X/765/1/39) [[PDF]](https://iopscience.iop.org/article/10.1088/0004-637X/765/1/39/pdf)

  In BibTeX::

    @ARTICLE{Nyx,
               author = {{Almgren}, A.~S. and {Bell}, J.~B. and {Lijewski}, M.~J. and {Luki{\'c}}, Z. and {Van Andel}, E.},
               title = "{Nyx: A Massively Parallel AMR Code for Computational Cosmology}",
               journal = {The Astrophysical Journal},
               archivePrefix = "arXiv",
               eprint = {1301.4498},
               keywords = {gravitation, hydrodynamics, methods: numerical },
               year = 2013,
               month = mar,
               volume = 765,
               eid = {39},
               pages = {39},
               doi = {10.1088/0004-637X/765/1/39}
    }

The development of AMReX library is led by the
Center for Computational Sciences and Engineering / Lawrence Berkeley
National Laboratory. Nyx development is done collaboratively, including the Computational Cosmology Center and CCSE. 

All of Nyx's development is done in the public github repository—anyone can see the updates as they are done.  
We are always happy to have new developers as part of the Nyx team. 

Fork the Nyx git repository on github, make your changes, and issue a pull request against the development branch 
of the github repo. Any level of changes are welcome, including documentation, bug fixes, new test problems, and new solvers.

 .. todo::
    Describe developer/contributor/author list


 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _PostProcessing:

Analyzing Outputs
=================

For analysis and visualizations purposes, Nyx outputs plotfiles with user-defined grid quantities and/or
particle data.  This file format is native to the AMReX framework, but there is also an option to use HDF5
outputs (this is implemented and is in a testing/optimization phase now).


Visualization
-------------

There are several visualization tools that can be used for AMReX plotfiles,
specifically VisIt, ParaView, and yt.

See the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/Visualization_Chapter.html>`_ for
guidance on how to use each of these tools.

The following inputs control when plotfiles are written

+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
|                         | Description                                                           |   Type      | Default   |
+=========================+=======================================================================+=============+===========+
| amr.plot_int            | Frequency of plotfile output;                                         |    Int      | -1        |
|                         | if -1 then no plotfiles will be written                               |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| amr.plotfile_on_restart | Should we write a plotfile when we restart (only used if plot_int>0)  |   Bool      | False     |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| amr.plot_file           | Prefix to use for plotfile output                                     |  String     | plt       |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+

and the following control which variables appear in the plotfile

+----------------------------+---------------------------------------------------+------------+-----------+
|                            | Description                                       |   Type     | Default   |
+============================+===================================================+============+===========+
| amr.plot_vars              | Name of state variables to be in plotfiles        |   Strings  | All       |
+----------------------------+---------------------------------------------------+------------+-----------+
| amr.plot_dervive_plot_vars | Name of derived variables to be in plotfiles      |   Strings  | All       |
+----------------------------+---------------------------------------------------+------------+-----------+


Nyx also easily interfaces with two post-processing suites, Reeber used for halo finding
and Gimlet, used for calculating different summary statistics.


Reeber
------

Reeber uses topological methods to construct merge trees of scalar fields.
These trees are effectively parameter-independent and contain a complete
description of the field topology. In the context of Nyx, the field of interest
is usually the total matter density. Nyx then queries the merge tree with user-defined
runtime parameters in order to locate the volume-averaged center of dark matter
halos. The same tree can be queried with any number of such parameters to find
halos with different mass/density thresholds.  Reeber is publicly available on
`GitHub <https://github.com/mrzv/reeber>`_.


Gimlet
------

Gimlet computes a variety of quantities about the simulation, including optical
depths, Lyman-alpha fluxes, power spectra (both 1-D "line-of-sight" as well as
fully 3-D), and probability distribution functions. These suites are fully
MPI-parallel and can be run either "in situ" or "in-transit", or with a
combination of both. Gimlet code is not yet publicly available, but we are working
on releasing it.


Do It Yourself
--------------

Nyx and AMReX provide the capability for the user to execute an arbitrary
post-processing workflow.  In ``Util/Diagnostics/`` we provide a simple example
of opening an AMReX plotfile, reading and manipulating data in it.  That can be a
starting place for building analysis tools for any specific need.
.. role:: cpp(code)
   :language: c++

.. _Chap:Debugging:

Debugging
=========

Debugging is an art.  Everyone has their own favorite method.  Here we
offer a few tips we have found to be useful.

Compiling in debug mode (e.g., :cpp:`make DEBUG=TRUE` for gmake users;
:cpp:`cmake -DDEBUG` for cmake users) and running with
:cpp:`amrex.fpe_trap_invalid=1` in the inputs file can be helpful.
In debug mode, many compiler debugging flags are turned on and all
:cpp:`MultiFab` s are initialized to signaling NaNs.  The
:cpp:`amrex.fpe_trap_invalid` parameter will result in backtrace files
when a floating point exception occurs.  One can then examine those
files to track down the origin of the issue.

Several other ways to look at the data include:

1) Writing a :cpp:`MultiFab` to disk with

.. highlight:: c++

::

    VisMF::Write(const FabArray<FArrayBox>& mf, const std::string& name);

and examining it with ``Amrvis`` (section :ref:`sec:amrvis` in the AMReX documentation).

2) You can also use the :cpp:`print_state` routine: 

.. highlight:: c++

::

    void print_state(const MultiFab& mf, const IntVect& cell, const int n=-1);

which outputs the data for a single cell.

3) If you want to compare old and new plotfiles, 

.. highlight:: c++

::

    fcompare --infile1 plt00000_run1 --infile2 plt00000_run2 --diffvar u_g

will print out the maximum absolute and relative differences between the two plotfiles
for each variable and will also create a new plotfile "diffs" that contains the difference
in u_g (in this case) between the two plotfiles.

The :cpp:`fcompare` executable can be built in AMReX (go to amrex/Tools/Plotfile and type "make").

4) Valgrind is another useful debugging tool.  Note that for runs using
more than one MPI process, one can tell valgrind to output to different 
files for different processes.  For example,

.. highlight:: console

::

    mpiexec -n 4 valgrind --leak-check=yes --track-origins=yes --log-file=vallog.%p ./Nyx3d.exe ...
Building Nyx with CMake
=======================

CMake build is a two-step process. First ``cmake`` is invoked to create
configuration files and makefiles in a chosen directory (``builddir``).
Next, the actual build is performed by invoking ``make`` from within ``builddir``.
If you are new to CMake, `this short tutorial <https://hsf-training.github.io/hsf-training-cmake-webpage/>`_
from the HEP Software foundation is the perfect place to get started with it.

The CMake build process for Nyx is summarized as follows:

#. Create and enter the build directory:

   .. highlight:: console

   ::

       mkdir /path/to/builddir
       cd    /path/to/builddir

#. Perform the configuration step:

   .. highlight:: console

   ::

      cmake [nyx options] [dependencies options] -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] /path/to/Nyx/repo


   In the above snippet, ``[nyx options]`` indicates one or more of the options
   for the customization of the build listed in the :ref:`table <tab:nyxcmakeoptions>` below,
   whereas ``[dependecies options]`` are one or more of the CMake configuration options for the Nyx dependencies,
   i.e. `AMReX CMake options <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake>`_,
   and SUNDIALS. If the option ``CMAKE_BUILD_TYPE`` is omitted, ``CMAKE_BUILD_TYPE=Release`` is assumed.
   For example, to enable AMReX profiling capabilities in Nyx, configure as follows:

   .. code:: shell

             > cmake [nyx options] -DAMReX_TINY_PROFILE=yes -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] ..

#. Build the executable you are interested in (a small
   version of the Santa Barbara problem in this example):

   .. highlight:: console

   ::

      cd /path/to/builddir/Nyx/Exec/MiniSB
      make


   The resulting executable is ``nyx_MiniSB``.

.. note::
   **Nyx requires CMake 3.14 or higher.**


.. raw:: latex

   \begin{center}

.. _tab:nyxcmakeoptions:

.. table:: Nyx configuration options

   +-----------------+------------------------------+------------------+-------------+
   | Option name     | Description                  | Possible values  | Default     |
   |                 |                              |                  | value       |
   +=================+==============================+==================+=============+
   | CMAKE\_CXX\     | User-defined C++ flags       | valid C++        | None        |
   | _FLAGS          |                              | compiler flags   |             |
   +-----------------+------------------------------+------------------+-------------+
   | CMAKE\_CUDA\    | User-defined CUDA flags      | valid CUDA       | None        |
   | _FLAGS          |                              | compiler flags   |             |
   +-----------------+------------------------------+------------------+-------------+
   | BUILD\_SHARED\  | Build dependencies as shared | no/yes           | no          |
   | _LIBS           | libraries                    |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_MPI        | Enable build with MPI        | no/yes           | yes         |
   |                 |                              |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_MPI\_THREAD| Concurrent MPI calls from    | no/yes           | yes         |
   | \_MULTIPLE      | multiple threads             |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_OMP        | Enable build with OpenMP     | no/yes           | no          |
   |                 |                              |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_GPU\_      | On-node, accelerated GPU \   | NONE             | NONE,SYCL,\ |
   | BACKEND         | backend                      |                  | CUDA,HIP    |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_HYDRO      | Run with baryon hydro solve  | no/yes           | yes         |
   |                 | and fields                   |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_HEATCOOL   | Run with H and He heating-   | no/yes           | no          |
   |                 | cooling effects              |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_CONST\_    | Don't evolve H and He, treat | no/yes           | no          |
   | SPECIES         | them as constant             |                  |             |
   +-----------------+------------------------------+------------------+-------------+
   | Nyx\_CGS        | Evolve quantities in CGS     | no/yes           | no          |
   |                 | units instead of code units  |                  |             |
   +-----------------+------------------------------+------------------+-------------+
.. raw:: latex

   \end{center}


Working with Git submodules
~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default Nyx CMake searches the system for existing installations of the required dependencies
(AMReX is always required, SUNDIALS may be required depending on the configuration options).
If the required dependencies are not found on the system, Nyx CMake will automatically checkout
and build those dependencies as part of its build process. To this end, Nyx CMake relies on git
submodules to checkout the AMReX and/or SUNDIALS git repositories. In what follows, we will
focus on the AMReX git submodule only, but the same concepts apply unchanged to the
SUNDIAL submodule as well.


If the AMReX submodule is not initialized, Nyx CMake will initialize it and checkout
the commit the submodule is pointing to. If instead the AMReX  submodule has already
been manually initialized and a custom commit has been checked out, that commit will
be used. For Nyx development or testing, you may need to build with a different
branch or release of AMReX.

The ``subprojects/amrex`` directory is a Git repo, so use all standard Git
commands. Either ``cd subprojects/amrex`` to run Git commands in the ``amrex``
directory, or use ``git -C subprojects/amrex`` in the Nyx repo. For
instance, to build with the ``my-amrex-branch`` branch of the AMReX repo:

.. code:: shell

    > git -C subprojects/amrex checkout my-amrex-branch
    > git status
    ...
    modified:   subprojects/amrex (new commits)

The branch ``my-amrex-branch`` will then be used when building Nyx.

To revert to the default version of the AMReX submodule, run ``git submodule
update``:

.. code:: shell

    > git submodule update
    > git status
    ...
    nothing to commit, working tree clean

You can edit, commit, pull, and push AMReX changes from ``subprojects/amrex``.
AMReX development is outside the scope of this document. Run ``git status`` in
the top-level Nyx repo to see whether the AMReX submodule has new commits,
modified files, or untracked files.

To update the AMReX submodule referenced by Nyx:

.. code:: shell

    > git -C subprojects/amrex checkout UPDATED_AMREX_COMMIT_SHA1
    > git add subprojects/amrex
    > git commit -m 'Updating AMReX version'

This will only update the AMReX SHA-1 referenced by Nyx. Uncommitted AMReX
changes and untracked AMReX files under ``subprojects/amrex`` are not added by
``git add subprojects/amrex``. (To commit to the AMReX repo, change directories
to ``subprojects/amrex`` and run Git commands there, before ``git add
subprojects/amrex``.)

.. note::

    Only update the AMReX submodule reference in coordination with the other
    Nyx developers!


.. _sec:build:external:

Using existing installations of required dependencies
-----------------------------------------------------

You may prefer to build Nyx against an AMReX and/or SUNDIALS installation already
present on your system. Unless these installations are located in standard system
paths, you need to tell Nyx CMake where to look for them.

.. code:: shell

    > cmake -DCMAKE_BUILD_TYPE=[Debug|Release|RelWithDebInfo|MinSizeRel] [nyx options] -DAMReX_ROOT=/absolute/path/to/amrex/installdir /path/to/Nyx/repo

In the example above, ``-DAMReX_ROOT=/absolute/path/to/amrex/installdir`` instructs CMake to search
``/absolute/path/to/amrex/installdir`` before searching system paths for an available AMReX installation.
``AMReX_ROOT`` can also be set as an environmental variable instead of passing it as a command line option.
Similarly, you can define a ``SUNDIALS_ROOT`` variable, either via command line or the environment, to
teach CMake where to look for SUNDIALS. Choose one of the ``CMAKE_BUILD_TYPE`` to control the level of
optimization, the option ``-DCMAKE_BUILD_TYPE=Release`` will give the same defaults as GMake.


Few more notes on building Nyx
-----------------------------------

The system defaults compilers can be overwritten as follows:

.. code:: shell

    > cmake -DCMAKE_CXX_COMPILER=<c++-compiler> -DCMAKE_Fortran_COMPILER=<f90-compiler> [options]  ..

When building on a platform that uses the ``module`` utility, use either
the above command (with full path to the compilers) or the following:

.. code:: shell

    > cmake -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn [options] ..

Nyx uses the same compiler flags used to build AMReX, unless
``CMAKE_Fortran_FLAGS``/``CMAKE_CXX_FLAGS`` is explicitly provided, or
the environmental variables ``FFLAGS``/``CXXFLAGS`` are set.


For GPU builds, Nyx relies on the `AMReX GPU build infrastructure <https://amrex-codes.github.io/amrex/docs_html/GPU.html#building-with-cmake>`_
. The target architecture to build for can be specified via the AMReX configuration option ``-DAMReX_CUDA_ARCH=<target-architecture>``,
or by defining the *environmental variable* ``AMREX_CUDA_ARCH`` (all caps). If no GPU architecture is specified,
CMake will try to determine which GPU is supported by the system.


Building Nyx for Cori (NERSC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard build
--------------

For the Cori cluster at NERSC, you first need to load/unload modules required to build Nyx.

.. code:: shell

    > module unload altd
    > module unload darshan
    > module load cmake/3.14.0

The default options for Cori are the **Haswell** architecture and **Intel** compiler, if you want to compile with the **Knight's Landing (KNL)** architecture:

.. code:: shell

    > module swap craype-haswell craype-mic-knl

Or use the **GNU** compiler:

.. code:: shell

    > module swap PrgEnv-intel PrgEnv-gnu

Now Nyx can be built.

.. note::

    The load/unload modules options could be saved in the `~/.bash_profile.ext`


GPU build
---------

To compile on the GPU nodes in Cori, you first need to purge your modules, most of which won't work on the GPU nodes

.. code:: shell

    > module purge

Then, you need to load the following modules:

.. code:: shell

    > module load cgpu gcc/7.3.0 cuda/11.1.1 openmpi/4.0.3 cmake/3.14.4

Currently, you need to use OpenMPI; mvapich2 seems not to work.

Then, you need to use slurm to request access to a GPU node:

.. code:: shell

    > salloc -N 1 -t 02:00:00 -c 10 -C gpu -A m1759 --gres=gpu:8 --exclusive

This reservers an entire GPU node for your job. Note that you can’t cross-compile for the GPU nodes - you have to log on to one and then build your software.

Finally, navigate to the base of the Nyx repository and compile in GPU mode:

.. code:: shell

    > cd Nyx
    > mkdir build
    > cd build
    > cmake -DNyx_GPU_BACKEND=CUDA -DAMReX_CUDA_ARCH=Volta -DCMAKE_CXX_COMPILER=g++ ..
    > make -j

For more information about GPU nodes in Cori -- `<https://docs-dev.nersc.gov/cgpu/>`_

Building Nyx for Summit (OLCF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the Summit cluster at OLCF, you first need to load/unload modules required to build Nyx.

.. code:: shell

    > module load gcc
    > module load cmake/3.14.0

To build Nyx for GPUs, you need to load cuda module:

.. code:: shell

    > module load cuda/11.0.3

To compile:

.. code:: shell

    > cd Nyx
    > mdkir build
    > cd build
    > cmake -DNyx_GPU_BACKEND=CUDA -DAMReX_CUDA_ARCH=Volta -DCMAKE_C_COMPILER=$(which gcc)  -DCMAKE_CXX_COMPILER=$(which g++)   -DCMAKE_CUDA_HOST_COMPILER=$(which g++)  -DCMAKE_CUDA_ARCHITECTURES=70 ..
    > make -j

For more information about the Summit cluster: `<https://www.olcf.ornl.gov/for-users/system-user-guides/summit/>`_

.. note::

    The load/unload modules options could be saved in the `~/.profile`
Running the Code
================

The Nyx executable reads run-time information from an “inputs” file which you designate on the command line.
Values can be specified either in the inputs file or on the command line.
If a value is specified on the command line, that value will override a value specified in the inputs file.

See the Inputs section for more detail about what parameters can be specified at run-time.

#. From within the build directory type:

   ::

      <executable-name> <inputs-name>

   ``<executable-name>`` is  ``Nyx3d.Linux.gnu.ex1`` if you built Nyx with GNU Make, or
   ``nyx_<name-of-build-directory>`` if you built Nyx with CMake.

   ``<inputs-name>`` for a small test problem is  ``inputs.32`` for the MiniSB example,
   and ``inputs.rt`` for the LyA example. Most executable directories have an ``inputs``
   for a larger problem, and an ``inputs.rt`` or ``inputs.regtest`` for regression-test
   sized problems.

   .. note::
      For certain HPC systems, you may want to have a run directory separate from your compile / build directory.

      In that case, copy the executable and inputs file from the build directory to your run directory on scratch.
      Runs starting from a ``binary_particle_file`` need an absolute path in <inputs-name> or a symlink in the run directory.
      Runs with heating-cooling must have access to the TREECOOL_middle file in the run directory, and ascent in-situ runs
      need access to ascent_actions.yaml.
      

#. You will notice that running the code generates directories that look like
   plt00000, plt00020, etc.,
   and chk00000, chk00020, etc. These are “plotfiles” and
   “checkpoint” files. The plotfiles are used for visualization,
   and the checkpoint files for restarting the code.

See the Visualization chapter for how to visualize these plotfiles.
.. _Chap:Structure:

Directory overview
==================

+---------------+--------------------------------------------------+
| DIrectory     | Description                                      |
+===============+==================================================+
| cmake         | CMake module files                               |
+---------------+--------------------------------------------------+
| Docs          | Source code for building the documentation       |
+---------------+--------------------------------------------------+
| Exec          | Directory for building with GNU make             |
+---------------+--------------------------------------------------+
| Source        | Source files                                     |
+---------------+--------------------------------------------------+
| subprojects   | Directory for Git submodules (CMake only)        |
+---------------+--------------------------------------------------+
.. role:: cpp(code)
   :language: c++

.. _SUNDIALS:

Compiling Nyx with SUNDIALS 6
===============================

The following steps describe how to compile Nyx with
SUNDIALS 6 support. This library is necessary for non-adiabatic
heating-cooling Nyx runs, where ``USE_HEATCOOL=TRUE``, such as
the ``Nyx/Exec/LyA`` directory.

In order to use SUNDIALS:

#. We suggest using the Github mirror:
   https://github.com/LLNL/sundials and picking the type of
   parallelism that is appropriate for your architecture.

   To install with cuda and openmp support:
   
   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      INSTALL_PREFIX=$(pwd)/instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_C_COMPILER=$(which gcc)  \
      -DCMAKE_CXX_COMPILER=$(which g++)   \
      -DCMAKE_CUDA_HOST_COMPILER=$(which g++)    \
      -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CUDA_FLAGS="-DSUNDIALS_DEBUG_CUDA_LASTERROR" \
      -DSUNDIALS_BUILD_PACKAGE_FUSED_KERNELS=ON \
      -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG" \
      -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"  \
      -DCUDA_ENABLE=ON  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=ON   \
      -DF2003_INTERFACE_ENABLE=OFF   \
      -DSUNDIALS_INDEX_SIZE:INT=32   \
      -DCUDA_ARCH=sm_70 ../
      make -j8
      make install -j8

   To install with openmp and no cuda support:
         
   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      INSTALL_PREFIX=$(pwd)/instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_C_COMPILER=$(which gcc)  \
      -DCMAKE_CXX_COMPILER=$(which g++)   \
      -DCMAKE_CUDA_HOST_COMPILER=$(which g++)    \
      -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG" \
      -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"  \
      -DCUDA_ENABLE=OFF  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=ON   \
      -DF2003_INTERFACE_ENABLE=OFF   \
      -DSUNDIALS_INDEX_SIZE:INT=32 ../
      make -j8
      make install -j8

   To install with HIP support (with ROCm 4.5):

   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=$(pwd)/../instdir  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_C_COMPILER=$(which clang) \
      -DCMAKE_CXX_COMPILER=$(which hipcc) \
      -DEXAMPLES_INSTALL_PATH=$(pwd)/../instdir/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_FLAGS_RELEASE=-O3 \
      -DCMAKE_CXX_FLAGS_RELEASE=-O3 \
      -DCUDA_ENABLE=OFF  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=OFF   \
      -DF2003_INTERFACE_ENABLE=OFF \
      -DENABLE_HIP=ON \
      -DEXAMPLES_INSTALL=OFF ../
      make -j8
      make install -j8


   To install with SYCL support:

   ::

      #!/bin/bash
      set -e
      git clone https://github.com/LLNL/sundials
      cd sundials
      mkdir builddir instdir
      cd builddir
      cmake \
      -DCMAKE_INSTALL_PREFIX=$(pwd)/../instdir  \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -DCMAKE_CXX_COMPILER=$(which dpcpp)  \
      -DCMAKE_CXX_STANDARD=17 \
      -DEXAMPLES_INSTALL_PATH=$(pwd)/../instdir/examples \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS_RELEASE=-O3  \
      -DCUDA_ENABLE=OFF  \
      -DMPI_ENABLE=OFF  \
      -DOPENMP_ENABLE=OFF   \
      -DF2003_INTERFACE_ENABLE=OFF \
      -DENABLE_SYCL=ON ../
      make -j8
      make install -j8

#. Note that we give these examples for the gnu compiler or the appropriate parallel compiler.
   The compiler chosen needs to be consistent with Nyx's GNUMakefile
   variable COMP to ensure matching OMP runtime libraries for use with the OpenMP NVector. 

#. ``CUDA_ARCH`` must be set to the appropriate value for the GPU being targeted

#. For more detailed instructions for installing SUNDIALS with different flags and versions see
   the `SUNDIALS documentation <https://computing.llnl.gov/projects/sundials/sundials-software>`_.

#. In the ``GNUmakefile`` for the application which uses the interface to SUNDIALS, add
   ``USE_SUNDIALS = TRUE`` and ``SUNDIALS_ROOT=${INSTALL_PREFIX}``. Note that one must define the
   ``SUNDIALS_LIB_DIR`` make variable to point to the location where the libraries are installed
   if they are not installed in ``${INSTALL_PREFIX}/lib``. Note the default location
   for 64 is ``${INSTALL_PREFIX}/lib64``, which we override with ``-DCMAKE_INSTALL_LIBDIR=lib``.

#. If the application uses the SUNDIALS CVODE time integrator package, then the variable
   ``USE_CVODE_LIBS = TRUE`` should also be added in the ``GNUmakefile`` for the application.
   If the application used the SUNDIALS ARKode time integrator package, then the variable
   ``USE_ARKODE_LIBS = TRUE`` should be added.

Note that SUNDIALS can also be installed via Spack:

   ::
      
      spack install sundials+cuda+openmp
  
Building Nyx with GNU Make
============================

Nyx is built on top of the AMReX framework so you must
download AMReX in order to build Nyx with GNU Make.

.. raw:: latex

   \vspace{.1in}

#. Clone/fork the AMReX repository:

   ::

       git clone https://github.com/AMReX-Codes/amrex

   You will want to periodically update AMReX by typing

   ::

       git pull

   in the ``amrex/`` directory.

.. note::
   When you check out AMReX (and Nyx), you will get the development
   branch.  Active development is done on this branch; monthly
   tags are made of this version that are compatible with the same
   monthly tag of AMReX itself.

#. Set the environment variable ``AMREX_HOME`` to point to
   your local copy of the AMReX repository.
   You can add this to your ``.bashrc`` as:

   ::

       export AMREX_HOME=/path/to/amrex/

   or to your ``.cshrc`` as:

   ::

       setenv AMREX_HOME /path/to/amrex/

   where ``/path/to/amrex/`` is the full path to the
   amrex directory.

#. Clone/fork the Nyx repository:

   ::

       git clone https://github.com/AMReX-Astro/Nyx

   As with AMReX development on Nyx is done in the
   ``development`` branch, so you should work there if you want
   the latest source.


#. Choose which executable to compile, for example MiniSB or LyA.

   MiniSB is a small version of the Santa Barbara problem, and LyA is a Lyman-:math:`\alpha` 
   forest simulation for investigating the large-scale structure formation of the universe.

   For some executables, namely LyA, AMR-density, AMR-zoom, and LyA_Neutrinos, Nyx uses a more complicated model for the heating-cooling.
   This requires you to also install a matching sundials installation to support the ODE solve. To install a stand-alone copy of Sundials, see :ref:`Compiling Nyx with SUNDIALS 5<sundials>`

#. From the directory in which you checked out Nyx, change directory to your build directory by typing for the small Santa-Barbara problem:

   ::

       cd Nyx/Exec/MiniSB

   or for the Lyman-:math:`\alpha` problem:

   ::

       cd Nyx/Exec/LyA

#. In your build directory, edit the GNUmakefile, and set

   ``COMP = your favorite compiler (e.g, gnu, Intel)``

   ``DEBUG = FALSE``

   We like ``COMP = gnu``.

   More information on the AMReX GNU Make setup can be found
   `here <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_.

   All executables should work with MPI+CUDA by setting ``USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=TRUE``.

   HIP and DPC++ builds are under development and can be tested by compiling with ``USE_HIP=TRUE``  and ``USE_DPCPP=TRUE``  respectively.

   .. note::
      For executables with ``USE_HEATCOOL=TRUE`` in their GNUmakefile, a matching Sundials implementation is required. If Sundials is built with ``-DSUNDIALS_BUILD_PACKAGE_FUSED_KERNELS=ON``, Nyx should be built with ``USE_FUSED=TRUE``.
      The flag ``USE_FUSED`` tells the Nyx compile whether you compiled Sundials with fused cuda kernels. The default assumption is that non-cuda Nyx compiles set ``USE_FUSED=FALSE`` to match Sundials being built without fused cuda kernels.
      Starting with Sundials version 5.7.0, set ``USE_SUNDIALS_SUNMEMORY=TRUE`` to compile the optional Sundials SunMemory to AMReX Arena interface for GPU memory reuse.

#. Now type “make”. The resulting executable will look something like
   “Nyx3d.Linux.gnu.ex”, which means this is a 3-d version of the code,
   made on a Linux machine, with ``COMP = gnu``.

   Note that if you build with ``USE_MPI = TRUE`` in the GNUMakefile, then the
   name of the code will be something like “Nyx3d.Linux.gnu.MPI.ex”
