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
