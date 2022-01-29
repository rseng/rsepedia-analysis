METALWALLS
==========

[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://gitlab.com/ampere2/metalwalls/-/wikis/home)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02373/status.svg)](https://doi.org/10.21105/joss.02373)
[![Repository](https://img.shields.io/badge/Zenodo-10.5281/zenodo.4912611-blue)](https://doi.org/10.5281/zenodo.4912611)

MetalWalls (MW) is a molecular dynamics code dedicated to the modelling of electrochemical systems. Its main originality is the inclusion of a series of methods allowing to apply a constant potential within the electrode materials.

**Extended documentation is provided in the [WIKI](https://gitlab.com/ampere2/metalwalls/-/wikis/home) section of the gitlab project.**
Details of the implemented force fields, thermodynamic ensembles and models, a description of the installation process and of input and output files, and guidelines for developers are given in detail.

In the following we reproduce installation instructions. To report bugs or contact the developpers, please raise an issue on the [Gitlab page](https://gitlab.com/ampere2/metalwalls/-/wikis/home).

# Reference

[A. Marin-Laflèche, M. Haefele, L. Scalfi, A. Coretti, T. Dufils, G. Jeanmairet, S. Reed, A. Serva, R. Berthin, C. Bacon, S. Bonella, B. Rotenberg, P.A. Madden, and M. Salanne. MetalWalls: A Classical Molecular Dynamics Software Dedicated to the Simulation of Electrochemical Systems. Journal of Open Source Software, 5, 2373, DOI:10.21105/joss.02373 (2020)](https://dx.doi.org/10.21105/joss.02373)

# Compiling

MW requires a Fortran compiler and the LAPACK library. The installation is based on a Makefile. A few machine dependent variables must be defined in the file *./config.mk* prior to invoking the `make` utility. The structure of the *config.mk* file is as follow:

```make
# Compilation options
F90 := fortran-compiler
F90FLAGS := compilation-flags
FPPFLAGS := preprocessor-flags
LDFLAGS := linker-flags
F2PY := path-to-f2py
F90WRAP := path-to-f90wrap
FCOMPILER := f2py-option (intel, intelem, gnu95...)
J := flag-to-specify-modfiles-output-dir (gnu: -J, intel: -module )

# Path to pFUnit (Unit testing Framework) -- optional
PFUNIT := path-to-pfunit
```

Some examples are provided for common use cases in the *./computers/* directory.
For example to compile MW on a Linux machine with the GNU compiler one can use the following parameters.
Here we assume that the MPI compiler wrapper is in the PATH of the user.

```make
# Compilation options
F90 := mpif90
F90FLAGS := -O2 -g
FPPFLAGS := -cpp
LDFLAGS := -llapack
F2PY := f2py
F90WRAP := f90wrap
FCOMPILER := gnu95
J := -J
# Path to pFUnit (Unit testing Framework)
PFUNIT := /opt/pfunit/pfunit-parallel
```

On a typical cluster with Intel Skylake processors and with Intel compiler
```make
# Compilation options
F90 := mpiifort
F90STDFLAGS := -g
F90OPTFLAGS := -O2 -xCORE-AVX512 -align array64byte
F90REPORTFLAGS :=
F90FLAGS := $(F90STDFLAGS) $(F90OPTFLAGS) $(F90REPORTFLAGS)
FPPFLAGS := -fpp
LDFLAGS := -mkl=cluster
F2PY := f2py
F90WRAP := f90wrap
FCOMPILER := intelem
J := -module
# Path to pFUnit (Unit testing Framework)
PFUNIT := $(ALL_CCCHOME)/opt/pfunit/pfunit-parallel
```

On a typical cluster with Nvidia V100 GPUs and pgfortran compiler
```make
# Compilation options
F90 := pgfortran
F90FLAGS := -fast -Mvect -m64 -Minfo=ccff -Mpreprocess -g
F90FLAGS := -tp=px -Minfo=ccff -Mpreprocess
F90FLAGS += -acc -ta=tesla:managed -Minline
F90FLAGS += -Minfo=accel,inline
FPPFLAGS := -DMW_SERIAL
LDFLAGS := -llapack -lblas
J := -module
```
**Warning** Continuous integration tests are not run on GPUs so you should always check the calculation on a CPU architecture before starting production.

Some internal flags can be defined in `F90FLAGS` to activate certain features:

- `-DMW_USE_PLUMED` to compile with the Plumed library
- `-DMW_SERIAL` to compile in serial mode
- `-DMW_CI` to allow unit tests

Note that on GPUs `-DMW_SERIAL` is enforced since the parallelization is made with OpenACC.

The Makefile is located in the root directory. The command to compile the code is simply `make`. It will create the object files in a dedicated build directory and produce the *mw* executable in the root directory.


# PLUMED

MW can be run with the PLUMED biased-MD library. To achieve this use the option `-DMW_USE_PLUMED` to be specified in the `FPPFLAGS` of the *./config.mk* file. PLUMED should be installed and compiled externally to MW, using the procedure on the [PLUMED website](https://www.plumed.org/). To then link PLUMED to MW, from the root directory of the code `./` type:
```bash
plumed patch --new mw2
plumed patch --patch --shared --engine mw2
```
See the `./example/plumed/` directory for an example of a MW run coupled to PLUMED.


# Python interface

It is possible (but not necessary) to compile MW as a python library using **f2py** and **f90wrap**. For this, you should define the variables F2PY, F90WRAP and FCOMPILER in the *./config.mk*, add the flag `-fPIC` to the `F90FLAGS` and use the command `make python`. Please be aware that this compiler flag may cause a decrease in the performance. For more information, see the [python interface page](https://gitlab.com/ampere2/metalwalls/-/wikis/python-interface).


# Testing
Two test suites are available in the *./tests/* folder, one that includes unit tests to check individual subroutines and another regression tests that run the code as a whole. Both suites are independent and can be run separately.
For user purposes, we recommend to only use regression tests.

## Regression tests

To run the regression test suite, you will need a working python interpreter with the **numpy** package installed.

Regression tests are reference test cases against which the code is compared.
To run the regression tests, type the command in the *./tests/* directory:

```bash
python regression_tests.py
```

To run a reduced version of the tests one can use:

```bash
python regression_tests.py -r
```

To only run a subset of the tests one can use:

```bash
python regression_tests.py -s <subset>
```

To run tests that use the python interface, one can specify the path to the python executable using:

```bash
python regression_tests.py -s python_interface -py path-to-python
```

The various tests subsets are:
* *nist*: energy comparison for the NIST validation case
* *benchmark*: forces, energies and charges comparison with LAMMPS for several systems
* *tosi_fumi*: comparison of forces, energies and stress tensor with PIM results for a NaCl system
* *pim*: comparison with PIM results
* *aim*: comparison with PIMAIM results
* *dihedrals*: comparison with reference data
* *matrix_inversion*: comparison of matrix, forces, energies, dipoles and charges with reference data
* *maze*: comparison of dipoles and electrodes calculated with MaZe method and with the conjugate gradient and matrix inversion method  
* *charge_neutrality*: comparison of charge calculation with conjugate gradient with symmetric and asymmetric potential difference
* *non_neutral*: comparison of matrix, forces, energies and charges with reference data for non neutral electrolyte
* *plumed*: plumed test
* *external_field*: comparison with reference data
* *thomas_fermi*: comparison with reference data
* *python_interface*: *nist* test case run with the python interface
* *steele*: comparison with reference data
* *piston*: comparison with reference data
* *four_site_model*: comparison with reference data
* *efg*: comparison with reference data
* *unwrap*: comparison with reference data
* *dip_plus_elec*: comparison of simultaneous calculation of dipoles and electrodes using conjugate gradient with reference data
* *dump_per_species*: comparison with reference data

## Unit tests

Unit tests aim at validating individual components of the code, we stress that they are used to monitor the code by developpers but are recommended for normal users.

In order to properly run the unit test suite you will need to install pFunit (the unit test framework).
To install pFUnit, download the version 3.3 from the [project page on Github](https://github.com/Goddard-Fortran-Ecosystem/pFUnit/releases/tag/3.3.3)
and follow the installation instructions in the *README.md* which we briefly reproduce here for a *bash* example. It is necessary to compile the MPI enabled version of pFUnit.

```bash
wget https://github.com/Goddard-Fortran-Ecosystem/pFUnit/archive/3.3.3.tar.gz
tar -zxvf 3.3.3.tar.gz
cd pFUnit-3.3.3/
export F90=gfortran
export F90_VENDOR=GNU
export MPIF90=mpif90
make tests MPI=YES
make install INSTALL_DIR=/path-to-pfunit/pfunit-parallel
```

To run the tests, MW has to be compiled with the compiling option `-DMW_CI`, to be specified in the *config.mk* file.
To launch the tests, do

```bash
make
make check
```
which will compile and run the unit tests. The `make check` command is a shortcut for

```bash
make mw_tests
cd tests/pFUnit
../../mw_tests
```

On batch processing system it is required to build the *mw_tests*
executable and create the appropriate job submission script.

# Running MW

## Running

Running a MW simulation requires two input files: a configuration file, *runtime.inpt*, and a data file, *data.inpt*.
Their format is described in the [system configuration page](https://gitlab.com/ampere2/metalwalls/-/wikis/system-configuration) and the [data input file page](https://gitlab.com/ampere2/metalwalls/-/wikis/data-input-file-format).
In a folder with the *data.inpt* and *runtime.inpt* files, run the executable with an MPI wrapper, for example:

```bash
mpirun -np 4 ./mw
```
Different flags can tune the executable behavior:

```bash
  -h, --help                 show a help message
  -v, --version              show program version number and exit
      --output-rank=VALUE    enables output on some of the ranks
                             possible VALUE are:
                               root     - only rank 0 performs output (default)
                               all      - all mpi processes perform output
                               r1[,r2]* - comma separated list of rank ids which perform output
```

## Restarting

The MW-generated restart files have the same format as the data files so one only has to rename them *data.inpt* and run the simulation again.
Be careful to the velocity creation keyword in the system configuration file *runtime.inpt*.

If you are running simulations using the Mass-Zero method (*maze*) with matrix inversion, the computed matrix can be given as input *maze_matrix.inpt*. Similarly, if you are using the *matrix_inversion* algorithm, the computed matrix can be given as input *hessian_matrix.inpt*. When restarting the simulation, these matrices will be read instead of computed from scratch.

# Funding

The development of MW has received support from:
* [EoCoE](http://www.eocoe.eu), a project funded by the European Union Contracts No. H2020-EINFRA-2015-1-676629 and H2020-INFRAEDI-2018-824158.
* European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant Agreement No. 771294).
* French National Research Agency (Labex STORE-EX, Grant No. ANR-10-LABX-0076).
# Git essential command lines #

## Getting some help ##

`git help <verb>` get the manpage help for the `git <verb>` command.

## Customize your Git environment ##

Use `git config` to customize your Git environment.

### Set up your identity ###

git config --global user.name = "John Doe"
git config --global user.email = "johndoe@example.com"

### Set up your editor ###

git config --global core.editor vi
git config --global core.pager less

## Checking status ##

`git status`    gives output to determine which files are in which states
`git status -s` gives a simplified git status output

`git diff`          shows diff log between local and staged area
`git diff --staged` shows diff log between staged area and last commit

`git log`                                 lists the commits in reverse chronological order (most
     					  recent first)
`git log -p -2`                           lists the last 2 commits with diff following each entry
`git log --stat`                          lists the commits with statistics (list of modified files,
                                          number of lines added/removed in those files) following
					  each entry.
`git log --pretty=format:"%h %s" --graph` lists commits history with abbreviated commit hash (%h)
                                          and commit subject (%s) on one line with a nice little
					  ASCII graph showing branch and merge history

## Committing changes ##

`git commit`                       launches editor to fill commit messages (empty message will abort
                                   commit)
`git commit -m "commit message"`   commits with in-line commit message
`git commit -a`                    add all tracked and modified files to the commit, skipping
                                   staging step
`git commit --amend`               amend previous commit (add newly staged file and modify commit message)

## Removing files ##

`git rm FILE`          stops tracking FILE and removes it from the working tree
`git rm --cached FILE` stops tracking FILE but keeps it in the working tree

# Working with remotes #

`git remote -v` shows remotes shortnames and addresses
`git remote show <name>` shows more information about remote <name>
`git remote add <name> <url>` adds a remote named <name> for the repository at <url>
`git remote rename <oldname> <newname>` rename a remote named <oldname> to <newname>
`git remote remove <name>` removes a remote named <name>

`git fetch <name>` Fetch branches and/or tags from the remote repository named <name>
`git push <repository> <refspec>` Updates remote refs using local refs




