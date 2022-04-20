# How to install ELSI

## Prerequisites

The installation of ELSI makes use of the CMake software. Minimum requirements:

* CMake (3.0 or newer)
* Fortran compiler (Fortran 2003 compliant)
* C compiler (C99 compliant)
* MPI (MPI-3)
* BLAS, LAPACK, ScaLAPACK (with PBLAS and BLACS)

Enabling the PEXSI solver (highly recommended) requires:

* C++ compiler (C++11 compliant)

Enabling the EigenExa solver requires:

* EigenExa (2.11 or newer)

Enabling the SLEPc-SIPs solver requires:

* SLEPc (3.9 or newer)

Enabling the MAGMA solver requires:

* MAGMA (2.5 or newer)

Enabling the CUDA-based GPU acceleration in the ELPA solver requires:

* CMake (3.8 or newer)
* CUDA (10.0 or newer recommended)

## Quick start

We recommend preparing and editing configuration settings in a toolchain file
that can be read by CMake. Edit one of the templates provided in the
[`./toolchains`](./toolchains) directory. Then follow the steps below:

    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_TOOLCHAIN_FILE=YOUR_TOOLCHAIN_FILE ..
      ...
      ...
      -- Generating done
      -- Build files have been written to: /current/dir
    $ make [-j np]
    $ [make install]

`YOUR_TOOLCHAIN_FILE` should be the user's toolchain file. Commands in square
brackets are optional.

## Configuration

### 1) Compilers

CMake automatically detects and sets default compilers. The choices made by
CMake often work, but not necessarily guarantee the optimal performance. In some
cases, the compilers selected by CMake may not be the ones desired by the user.
Therefore, it is mandatory that the user explicitly sets the identification of
compilers in the following keywords:

* `CMAKE_Fortran_COMPILER`
* `CMAKE_C_COMPILER`
* `CMAKE_CXX_COMPILER` (only needed if building ELSI with PEXSI support)

It is highly recommended to specify the compiler flags, in particular the
optimization flags:

* `CMAKE_Fortran_FLAGS`
* `CMAKE_C_FLAGS`
* `CMAKE_CXX_FLAGS`

### 2) Solvers

The ELPA, libOMM, NTPoly, BSEPACK, and PEXSI solver libraries, as well as the
SuperLU\_DIST and PT-SCOTCH libraries (both required by PEXSI) are redistributed
through this ELSI package. Experienced users are encouraged to link the ELSI
interface against externally installed, better optimized solver libraries.
Relevant options are:

* `USE_EXTERNAL_ELPA`
* `USE_EXTERNAL_OMM`
* `USE_EXTERNAL_PEXSI`
* `USE_EXTERNAL_NTPOLY`
* `USE_EXTERNAL_BSEPACK`

All external libraries and include paths should be set via:

* `INC_PATHS`
* `LIB_PATHS`
* `LIBS`

Each of the above is a list separated by ` ` (space) or `;` (semicolon).
`INC_PATHS` and `LIB_PATHS` should be absolute paths. `LIBS` accepts three
formats:

* `-lelpa;-lpexsi;-lblas` (link line style)
* `elpa;pexsi;blas` (name of library)
* `libelpa.a;libpexsi.a;libblas.so` (full name of library)

On computers with NVIDIA GPUs, GPU acceleration based on CUDA and cuBLAS may be
enabled for the ELPA solver by `USE_GPU_CUDA`. The CUDA runtime and cuBLAS
libraries must be available in `LIB_PATHS` and `LIBS`. The nvcc compiler and
flags to compile the CUDA sources may be specified by setting these keywords:

* `CMAKE_CUDA_COMPILER`
* `CMAKE_CUDA_FLAGS`

Please note that in the current version of ELSI, the redistributed PEXSI and
BSEPACK solvers are not enabled by default. They may be switched on by
`ENABLE_PEXSI` and `ENABLE_BSEPACK`, respectively. In addition, `ENABLE_SIPS`,
`ELABLE_EIGENEXA`, and `ENABLE_MAGMA` may be used to enable support for the
SLEPc, EigenExa, and MAGMA solvers, respectively. These libraries are not
redistributed with ELSI, thus must be installed separately by the user.

### 3) Tests

Building ELSI test programs may be enabled by `ENABLE_TESTS`. Then, the
compilation may be verified by

    $ make test

or

    $ ctest

The tests may fail to start if launching MPI jobs is prohibited on the user's
working platform.

## Beyond this guide

A complete description of the ELSI build system is available in
[`./doc/elsi_manual.pdf`](./doc/elsi_manual.pdf).

## Troubleshooting

For comments, feedback, and suggestions, please
[contact the ELSI team](mailto:elsi-team@duke.edu).

Copyright (c) 2015-2021, the ELSI team. All rights reserved.
# ELSI changelog

## v2.9.0 (April 2022)

The indices have been updated to integer8 instead of integer4. 
This allows us to use a relatively small number of MPI ranks for relatively large matrices. 
This update is essential for accelerators (GPU,TPU,...). 
The interface is updated so that both integer8 and integer4 should work.

Add Cholesky extrapolation for density matrix when used with elsi restart files. 
ovlp_old and dm are read from disk and not stored in the elsi handle, thus, a separate interface was needed.

The single precision is used for forward and backward transformation to the standard eigenvalue problem. 
This further save the computational cost in the mixed precision calculations.

We updated the ELPA version to 2021.11.001, the old ELPA version is still the default. 
While the new version can be switched on with the cmake option USE_ELPA_2021.

We provide a description and examples for elsipy.

We update EigenExa to 2.11.

## v2.8.3 (July 2021)
* suggest BLACS distribution

## v2.8.2 (June 2021)
* Detailed error message for wrong occupation number: bug fix

## v2.8.1 (March 2021)
* Cython interface for ELSI -- elsipy
* Detailed error message for wrong occupation number

## v2.8.0 (March 2021)
* Interface for the non-Aufbau occupation.

## v2.7.1 (March 2021)

### ELSI interface
* Fixed bugs in the frozen core approximation code.

### Known issues
* The ELPA code cannot be compiled with the NAG Fortran compiler, due to the
  use of GNU extensions in ELPA.
* Depending on the choice of k-points, the complex PEXSI solver may randomly
  fail at the inertia counting stage.

## v2.7.0 (March 2021)

### ELSI interface
* Added support for frozen core approximation when using the dense eigensolver
  interfaces with ELPA and LAPACK.

## v2.6.4 (November 2020)

### ELSI interface
* Computation of density matrix from eigenvectors was made more robust.

## v2.6.3 (November 2020)

### ELSI interface
* Computation of chemical potential and occupation numbers was made more robust.

### PEXSI
* Updated redistributed (PT-)SCOTCH source code to version 6.1.0.

### NTPoly
* Updated redistributed NTPoly source code to version 2.5.1.

### EigenExa
* Interface compatible with EigenExa 2.6.

## v2.6.2 (July 2020)

### ELPA
* Fixed a performance regression of the ELPA2 generic kernel.

## v2.6.1 (June 2020)

### PEXSI
* Removed an improper abort from the error handling code of PEXSI.

## v2.6.0 (June 2020)

### ELSI interface
* C compiler and MPI-3 have become mandatory to build ELSI.
* Added an option to choose which sparsity pattern to use when converting input
  dense matrices to the sparse format used by the solver.

### ELPA
* Updated redistributed ELPA source code to version 2020.05.001, which supports
  single-precision calculations, autotuning of runtime parameters, and (NVIDIA)
  GPU acceleration.

### PEXSI
* AAA method has become the default pole expansion method in PEXSI.
* Increased default number of poles from 20 to 30.
* Improved accuracy of pole expansion based on minimax rational approximation.
* Updated redistributed (PT-)SCOTCH source code to version 6.0.9.

### NTPoly
* Updated redistributed NTPoly source code to version 2.5.0.

### SLEPc-SIPs
* Interface compatible with PETSc 3.13 and SLEPc 3.13.

## v2.5.0 (February 2020)

### ELSI interface
* Added utility subroutines to retrieve the internally computed eigenvalues,
  eigenvectors, and occupation numbers when using the density matrix solver
  interfaces with an eigensolver.
* Fixed Marzari-Vanderbilt broadening.

### Solvers
* Added support for the Bethe-Salpeter eigensolvers in the BSEPACK library.

### ELPA
* Interface for externally linked ELPA compatible with ELPA 2019.11.

### NTPoly
* Updated redistributed NTPoly source code to version 2.4.0.

### PEXSI
* Updated redistributed SuperLU\_DIST source code to version 6.2.0.
* Added support for computing the electronic entropy via the free energy density
  matrix.

### BSEPACK
* Redistributed source code of BSEPACK 0.1.
* Added parallel BSE eigensolvers PDBSEIG and PZBSEIG.

## v2.4.1 (November 2019)

### ELSI interface
* Fixed energy-weighted density matrix computation when using ELPA and Fermi
  broadening.

### NTPoly
* Updated redistributed NTPoly source code to version 2.3.2.

## v2.4.0 (November 2019)

### ELSI interface
* Fixed density matrix computation when using ELPA and Methfessel-Paxton
  broadening.

### Solvers
* Added support for the tridiagonalization and pentadiagonalization eigensolvers
  implemented in the EigenExa library.
* Added support for the one-stage and two-stage tridiagonalization eigensolvers
  implemented in the MAGMA library.

### EigenExa
* Interface compatible with EigeExa 2.4.
* Added tridiagonalization eigensolver eigen\_s and pentadiagonalization
  eigensolver eigen\_sx.

### SLEPc-SIPs
* Interface compatible with PETSc 3.12 and SLEPc 3.12.

### MAGMA
* Interface compatible with MAGMA 2.5.
* Added one-stage and two-stage eigensolvers.

## v2.3.1 (July 2019)

### SLEPc-SIPs
* Fixed memory leaks in the redistributed SIPs code.
* Interface compatible with PETSc 3.11 and SLEPc 3.11.

## v2.3.0 (June 2019)

### ELSI interface
* Added density matrix extrapolation subroutines for sparse matrices.
* Extended the test suite to increase code coverage.

### ELPA
* Interface for externally linked ELPA compatible with ELPA 2019.05.

### NTPoly
* Updated redistributed NTPoly source code to version 2.3.1.
* Fixed complex matrix conversion from BLACS\_DENSE and GENERIC\_COO to NTPoly.

### PEXSI
* Fixed complex energy-weighted density matrix.
* Added support for linking ELSI against an externally compiled PEXSI.

## v2.2.1 (March 2019)

### ELSI interface
* Fixed LAPACK eigensolver for ill-conditioned overlap matrices.

### PEXSI
* Updated redistributed SuperLU\_DIST source code to version 6.1.1.

## v2.2.0 (February 2019)

### ELSI interface
* Added utility subroutines for geometry optimization and molecular dynamics
  calculations, including reinitialization of ELSI between geometry steps,
  density matrix extrapolation, and Gram-Schmidt orthogonalization of
  eigenvectors.
* Extended the test suite to increase code coverage.

### Matrix formats
* Added arbitrarily distributed coordinate (GENERIC\_COO) format.

### ELPA
* Interface for externally linked ELPA compatible with ELPA 2018.11.
* Fixed single-precision calculations with externally linked ELPA.
* Fixed internal ELPA two-stage real solver with AVX512 kernel.

### NTPoly
* Updated redistributed NTPoly source code to version 2.2.

### OMM
* Fixed libOMM Cholesky flavor with random initial guess.

### PEXSI
* Updated redistributed PEXSI source code to version 1.2.0.
* Updated redistributed SuperLU\_DIST source code to version 6.1.0.

## v2.1.0 (October 2018)

### ELSI interface
* Adopted literature definition of the electronic entropy.
* Added subroutines to query the version number and date stamp of ELSI.

### Solvers
* Added support for the density matrix purification methods implemented in the
  NTPoly library. Implementation of the same methods with dense linear algebra
  has been removed.

### ELPA
* For externally linked ELPA, added options to perform single-precision
  calculations and to automatically tune the internal runtime parameters of the
  solver.
* Interface for externally linked ELPA compatible with ELPA 2018.05.

### NTPoly
* Redistributed source code of NTPoly 2.0.
* Added canonical purification, trace correcting purification, 4th order trace
  resetting purification, and generalized hole-particle canonical purification
  methods.

### PEXSI
* Updated redistributed PEXSI source code to version 1.0.3, which returns the
  complex density matrix and energy-weighted density matrix instead of their
  transpose.

### SLEPc-SIPs
* Updated interface to support PETSc 3.9 and SLEPc 3.9.

## v2.0.2 (June 2018)

### PEXSI
* Updated redistributed PEXSI source code to version 1.0.1, which fixes the
  complex Fermi operator expansion routine.
* Downgraded redistributed (PT-)SCOTCH source code to version 6.0.0, as newer
  versions seem to be incompatible with PEXSI.

## v2.0.1 (June 2018)

### ELSI interface
* Switched to the [semantic versioning scheme](http://semver.org).
* Fixed building ELSI as a shared library with tests enabled.
* Improved stability when calling PBLAS routines pdtran and pztranc.

## v2.0.0 (May 2018)

### ELSI interface
* CMake build system has replaced the Makefile-based system.
* Added support for building ELSI as a shared library.
* Added support for spin channels and k-points.
* Added support for energy-weighted density matrix.
* Added support for electronic entropy calculations.
* Added support for complex sparse matrix formats for eigensolver and density
  matrix solver interfaces.
* Removed optional variables from mutator subroutines.
* Added matrix I/O subroutines using the MPI I/O standard.
* Removed TOMATO dependency for the test suite.
* Added a unified JSON output framework via the FortJSON library.

### Solvers
* Added support for the SLEPc-SIPs solver (PETSc 3.8 and SLEPc 3.8 required).
* Implemented density matrix purification with dense linear algebra operations.

### Matrix formats
* Added 1D block-cyclic compressed sparse column (SIESTA\_CSC) format.

### ELPA
* Updated redistributed ELPA source code to version 2016.11.001.
* Added AVX512 kernel.
* Made the two-stage solver default for all matrix sizes.
* Updated the interface for externally linked ELPA to the AEO version (ELPA
  release 2017.05 or later). GPU acceleration and GPU kernels may be enabled
  through the ELSI interface for externally linked ELPA.

### PEXSI
* Updated redistributed PEXSI source code to version 1.0.0.
* Reduced the default number of poles to 20 without sacrificing accuracy.
* Switched to the PT-SCOTCH library as the default sparse matrix reordering
  software.
* Redistributed SuperLU\_DIST 5.3.0 and (PT-)SCOTCH 6.0.5a libraries. Users may
  still provide their own SuperLU\_DIST library linked against any compatible
  sparse matrix reordering library.
* Removed ParMETIS as a mandatory external dependency for PEXSI.

## v1.0.0 (May 2017)

### Solvers
* ELPA (version 2016.11.001.pre)
* libOMM
* PEXSI (version 0.10.2)

### Matrix formats
* 2D block-cyclic dense (BLACS format)
* 1D block compressed sparse column (PEXSI format)
# Contributing to ELSI

Contributions to the ELSI project are very welcome.

## The first step

[Open an issue](https://gitlab.com/elsi_project/elsi_interface/-/issues) to
briefly describe the proposed changes to the code. This is to avoid the
situation where two or more developers work on something similar without knowing
each other.

## Suggested workflow

1. Fork the project on `gitlab.com`
2. Clone the forked repository
3. Create a new branch (with a meaningful name) in it
4. Commit any changes to the newly created branch
5. Push the branch to the forked repository
6. Submit a merge request to the ELSI repository

## Merge requests

* Put `WIP:` in the merge request title to indicate work in progress.
* An opened merge request can be updated by pushing new commits to the branch or
  amending existing commits. The branch to be merged must
  1. Be up to date with the upstream master branch
  2. Not contain any revert commits
  3. Not contain any merge commits
  4. Pass the continuous integration tests
* Always rebase the branch onto the upstream master, instead of merging the
  upstream master into the branch. This helps maintain a linear, readable, and
  therefore more useful git history.
* It is highly recommended to clean up the commits in a merge request by an
  interactive rebase, i.e., `git rebase -i upstream/master`.
* Please make sure that each merge request introduces only one logical change to
  the code.

## Commit messages

* Write meaningful commit messages in this format:
  1. A single "title" line of 50 characters or less summarizing the commit. The
     title should be written in the imperative mood, e.g., "Fix typo in doc"
     instead of "This commit fixed typo in doc" or "typo fixed". This aligns
     with the choice made by git itself.
  2. An empty line.
  3. An optional "body" containing multiple lines, each line containing 72
     characters or less.
* When applicable, refer to existing issues and/or merge requests.

## Coding style

* In general, avoid tabs, trailing whitespaces, and consecutive blank lines.
* Maximal line length is 80 characters. Use continuation for longer lines.
* Nested blocks are indented by 3 white spaces.
* Variable names follow the `lower_case_with_underscore` convention. Do not use
  upper case except for constants, e.g., `MPI_COMM_WORLD`.
* Subroutines that ar part of the public API should be documented with the
  [Doxygen](https://www.doxygen.nl/index.html) tool.
# ELSI - ELectronic Structure Infrastructure (v2.9.0)

## About

ELSI is a unified software interface designed for electronic structure codes to
connect with various high-performance eigensolvers and density matrix solvers.
For more information, visit the [ELSI interchange](https://elsi-interchange.org)
website.

## Installation

The standard installation of ELSI requires:

* CMake (3.0 or newer)
* Fortran compiler (Fortran 2003)
* C compiler (C99)
* C++ compiler (C++11, optional)
* MPI (MPI-3)
* BLAS, LAPACK, ScaLAPACK
* CUDA (optional)

Installation with recent versions of Cray, GNU, IBM, Intel, and NVIDIA (formerly
PGI) compilers has been tested. For a complete description of the installation
process, please refer to [`./INSTALL.md`](./INSTALL.md).

## More

A User's Guide is available at [`./doc/elsi_manual.pdf`](./doc/elsi_manual.pdf).
For comments, feedback, and suggestions, please
[contact the ELSI team](mailto:elsi-team@duke.edu).

Copyright (c) 2015-2022, the ELSI team. All rights reserved.
BSEPACK
=======

Contributers
------------
### Meiyue Shao and Chao Yang

Scalable Solvers Group
Computational Research Division
Lawrence Berkeley National Laboratory

Last update: version 0.1, August 2016


1. Introduction
---------------
BSEPACK is a parallel ScaLAPACK-style software library for computing all
eigenpairs of definite Bethe--Salpeter Hamiltonian matrices of the form

    H = [        A,        B
          -conj(B), -conj(A) ],

where

    [       A,       B
      conj(B), conj(A) ]

is Hermitian and positive definite.

The source code can be downloaded from the software homepage
<https://sites.google.com/a/lbl.gov/bsepack/>


2. How to Install
-----------------
The library is written in fixed format Fortran 90.  In addition to a Fortran
90/95 compiler, the following libraries are required.

  1. MPI, e.g., OpenMPI or MPICH.
  2. An optimized BLAS library, e.g., ATLAS or OpenBLAS,
     see <http://www.netlib.org/blas/> for a reference implementation.
  3. LAPACK, see <http://www.netlib.org/lapack/>
  4. ScaLAPACK (including BLACS and PBLAS),
     see <http://www.netlib.org/scalapack/>

Follow the instruction below to build the library:

  1. Download and unpack the BSEPACK archive.
  2. Modify the file `make.inc` to setup the compilers, linkers, and
     libraries.  Templates of this file are provided in the directory
     `MAKE_INC/`.
  3. Type `make all` to build the library and test programs.
  4. Run the test programs to check whether the installation is successful.

More detailed instructions regarding installation and testing can be found in
*BSEPACK User's Guide*.


3. Functionalities
------------------
The BSEPACK library provides two main functionalities:

  1. Diagonalizing the definite Bethe--Salpeter Hamiltonian matrix.
  2. Compute or estimate the absorption spectrum.

Solvers based on Tamm--Dancoff approximation (TDA) for both cases are also
provided.  Current release supports both real and complex data types, but only
for double precision.  Support to single precision real and complex data types
are *not* planned in future releases unless strong requests are received.

Examples for computing the eigenvalues and absorption spectrum are provided in
the directory `EXAMPLES/`.  The outputs are valid Matlab/Octave scripts that
can be used directly in Matlab/Octave for postprocessing.  To use your own
matrices, replace the input files (`EXAMPLES/input_real.txt` and/or
`EXAMPLES/input_complex.txt`) by yours.  Detailed usage is described in
*BSEPACK User's Guide*.


4. Comments/Questions/Bug Reports
---------------------------------
Please send your request to <myshao@lbl.gov>.


5. Selected References
----------------------
The dense eigensolver in BSEPACK is implemented based on [1].
The estimation of absorption spectrum uses algorithms in [2,3].

  [1] M. Shao, F. H. da Jornada, C. Yang, J. Deslippe, and S. G. Louie.
      Structure preserving parallel algorithms for solving the
      Bethe--Salpeter eigenvalue problem.
      Linear Algebra Appl., 488:148--167, 2016.
      DOI: 10.1016/j.laa.2015.09.036

  [2] J. Brabec, L. Lin, M. Shao, N. Govind, Y. Saad, C. Yang, and E. G. Ng.
      Efficient algorithms for estimating the absorption spectrum within
      linear response TDDFT.
      J. Chem. Theory Comput., 11(11):5197--5208, 2015.
      DOI: 10.1021/acs.jctc.5b00887

  [3] M. Shao, F. H. da Jornada, L. Lin, C. Yang, J. Deslippe, and S. G. Louie.
      A structure preserving Lanczos algorithm for computing the optical
      absorption spectrum.
      Avaliable as arXiv:1611.02348, 2016.
Project Overview
================================================================================

[![Build Status](https://travis-ci.org/william-dawson/NTPoly.svg?branch=travis-ci)](https://travis-ci.org/william-dawson/NTPoly)

NTPoly is a massively parallel library for computing the functions of sparse,
symmetric matrices based on polynomial expansions. For sufficiently sparse
matrices, most of the matrix functions in NTPoly can be computed in linear
time.

Set Up Guide
--------------------------------------------------------------------------------
NTPoly is freely available and open source under the MIT license. It can be
downloaded from the [Github](https://github.com/william-dawson/NTPoly)
repository. We of course recommend that you download a
[release version](https://github.com/william-dawson/NTPoly/releases)
to get started.

Installing NTPoly requires the following software:

* A Fortran Compiler.
* An MPI Installation (MPI-3 Standard+).
* CMake (Version 3.2+).

The following optional software can greatly enhance the NTPoly experience:

* BLAS: for multiplying dense matrices, if they emerge in the calculation.
* A C++ Compiler for building C++ bindings.
* Doxygen: for building documentation.
* Python (Version 2.7+): for testing.
* MPI4PY: for testing.
* SciPy: for testing.
* NumPy: for testing.
* SWIG (Version 3.0+): for building the Python bindings.

NTPoly uses CMake as a build system. First, take a look in the Targets
directory. You'll find a list of `.cmake` files which have example configurations
on popular systems. You should copy one of these files, and create your own
mymachine.cmake file. Then, cd into the Build directory, and type:
> cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/mymachine.cmake ..

There are a few options you can pass to CMake to modify the build. You can set
`-DCMAKE_BUILD_TYPE=Debug` for debugging purposes. You can set the install
directory using the standard `-DCMAKE_INSTALL_PREFIX=/path/to/dir`. You can
also set `-DFORTRAN_ONLY=YES` if you want to only build the Fortran interface.
Note that with just the Fortran interface, it is not possible to perform local
tests.

After that you can build using:
> make

And for the documentation:
> make doc

[Online documentation](https://william-dawson.github.io/NTPoly/documentation/) is also
available. Further details about the library can be found on the
[Wiki](https://github.com/william-dawson/NTPoly/wiki).
If you aren't cross compiling, you can perform local tests using:
> make test

Basic Theory
--------------------------------------------------------------------------------
The theory of matrix functions is a long studied branch of matrix algebra.
Matrix functions have a wide range of applications, including graph problems,
differential equations, and materials science. Common examples of matrix
functions include the matrix exponential:

> f(A) = e^A.

from the study of networks, or the inverse square root:

> f(A) = A^(-1/2)

from quantum chemistry. NTPoly is a massively parallel library that can be used
to compute a variety of matrix using polynomial expansions. Consider for example
the Taylor series expansion of a function *f(x)* .

> f(x) = f(0) + f'(0)x + f''(0)x^2/2! + ...

We can imagine expanding this from the function of a single variable, to a
function of a matrix:

> f(A) = f(0) + f'(0)A + f''(0)A^2/2! + ...

where matrices can be summed using matrix addition, and raised to a power
using matrix multiplication. At the heart of NTPoly are polynomial expansions
like this. We implement not only Taylor expansions, but also Chebyshev
polynomial expansions, and other specialized expansions based on the function
of interest.

When the input matrix *A* and the output matrix *f(A)* are sparse, we can
replace the dense matrix addition and multiplication routines with sparse
matrix routines. This allows us to use NTPoly to efficiently compute many
functions of sparse matrices.

Getting Start With Examples
--------------------------------------------------------------------------------
In the examples directory, there are a number of different example programs that
use NTPoly. You can check the ReadMe.md file in each example directory to
learn how to build and run each example. The simplest example is PremadeMatrix,
which includes sample output you can compare to.

Feature Outline
--------------------------------------------------------------------------------
The following features and methods have been implemented in NTPoly:

* General Polynomials
    * Standard Polynomials
    * Chebyshev Polynomials
    * Hermite Polynomials
* Transcendental Functions
    * Trigonometric Functions
    * Exponential and Logarithm
* Matrix Roots
    * Square Root and Inverse Square Root
    * Matrix *p* th Root
* Quantum Chemistry
    * Density Matrix Purification
    * Chemical Potential Calculation
    * Geometry Optimization
* Other
    * Matrix Inverse/Moore-Penrose Pseudo Inverse
    * Sign Function/Polar Decomposition
    * Load Balancing Matrices
    * File I/O

Citation
--------------------------------------------------------------------------------
A description of the techniques used in NTPoly can be found in the following
Computer Physics Communications paper:

> Dawson, William, and Takahito Nakajima. "Massively parallel sparse matrix
> function calculations with NTPoly." Computer Physics Communications (2017).

Please cite this paper in accordance to the practices in your field.

How To Contribute
--------------------------------------------------------------------------------
To begin contributing to NTPoly, take a look at the
[Wiki](https://github.com/william-dawson/NTPoly/wiki) pages. The
[Contributing Guide](https://github.com/william-dawson/NTPoly/blob/master/CONTRIBUTING.md)
provides an overview of best development practices. Additionally, there is a
[Adding New Functionality](https://github.com/william-dawson/NTPoly/wiki/Adding-New-Functionality-(Example))
page which documents how one would go about adding a matrix function to NTPoly.
# FortJSON #

## Overview ##

FortJSON is a JSON library written in Fortran 2003. It is designed with portability across HPC architectures in mind.

Requirements:

- Fortran 2003 compiler
- cmake 3.0+

Recommended:

- doxygen 1.8.11+

## Installation ##

FortJSON uses a standard cmake (out-of-source) build system for installation:

1. Create a build subdirectory in the root FortJSON directory and change directory into it
2. Generate the build system using cmake.  There are multiple ways to do this, but two common ways are:
   - Specify a toolchain file containing compilation settings, then run `cmake -DCMAKE_TOOLCHAIN_FILE=/path/to/toolchain/file ..` to generate the build system.  An example toolchain file named "the-gibson.cmake" is provided in the root directory.
   - Run `cmake ..` to generate the build system using settings automatically determined by cmake.  This approach may select the wrong Fortran compiler and will likely attempt to install FortJSON to a folder that you do not have write access to (for example, /usr/local/).
3. Run `make` to compile FortJSON in the build subdirectory
4. Run `make install` to install FortJSON into the install directory

FortJSON uses doxygen to generate documentation:

1. Enter into the doc subdirectory in the root FortJSON directory.
2. Run `doxygen Doxyfile` to generate HTML and LaTeX versions of the documentation.

## How to Use ##

Example program:

```fortran
program hello_world

   use FortJSON

   implicit none

   type(fjson_handle) :: fj_h

   integer(kind=i4) :: my_age
   real(kind=r8)    :: my_weight

   ! Open the file
   call fjson_open_file(fj_h, 66, "hello_world.json")
   call fjson_start_array(fj_h)

   ! Rex the Dog
   my_age = 1_i4
   my_weight = 5.0_r8

   call fjson_start_object(fj_h)
   call fjson_write_name_value(fj_h, "My Age", my_age)
   call fjson_write_name_value(fj_h, "My Weight", my_weight)
   call fjson_write_name_value(fj_h, "My Name", "Rex Jr.")
   call fjson_start_name_array(fj_h, "My Favorite Things")
   call fjson_write_value(fj_h, "Meat")
   call fjson_finish_array(fj_h)
   call fjson_start_name_object(fj_h, "My Favorite People")
   call fjson_write_name_value(fj_h, "Mom", "Precious")
   call fjson_finish_object(fj_h)
   call fjson_finish_object(fj_h)

   ! Tiger the Cat
   my_age = 15_i4
   my_weight = 15.0_r8

   call fjson_start_object(fj_h)
   call fjson_write_name_value(fj_h, "My Age", my_age)
   call fjson_write_name_value(fj_h, "My Weight", my_weight)
   call fjson_write_name_value(fj_h, "My Name", "Tiger")
   call fjson_start_name_array(fj_h, "My Favorite Things")
   call fjson_finish_array(fj_h)
   call fjson_start_name_object(fj_h, "My Favorite People")
   call fjson_finish_object(fj_h)
   call fjson_finish_object(fj_h)

   call fjson_write_value(fj_h, &
                          fjson_error_message(fjson_get_last_error(fj_h)))

   ! Close the file
   call fjson_finish_array(fj_h)
   call fjson_close_file(fj_h)

end program hello_world
```

Example output:

```json
[
  {
    "My Age": 1,
    "My Weight": 0.50000000E+01,
    "My Name": "Rex Jr.",
    "My Favorite Things": [
      "Meat",
    ],
    "My Favorite People": {
      "Mom": "Precious",
    }
  },
  {
    "My Age": 15,
    "My Weight": 0.15000000E+02,
    "My Name": "Tiger",
    "My Favorite Things": [
    ],
    "My Favorite People": {
    }
  },
  "FortJSON Error:  No error"
]
```

## Additional Information ##

- The JSON standard used for FortJSON is "ECMA-404 The JSON Data Interchange Standard".  This standard is 16 pages long and contains mostly whitespace and pictures.
- A simple description of the JSON grammer, as well as other JSON libraries, may be found at <http://json.org/>.
- Google maintains a JSON style guide at <https://google.github.io/styleguide/jsoncstyleguide.xml>.
# [Eigenvalue SoLvers for Petaflop-Applications (ELPA)](http://elpa.mpcdf.mpg.de)

## Current Release ##

The current release is ELPA 2020.05.001.rc1 The current supported API version
is 20190501. This release supports the earliest API version 20170403.

The old, obsolete legacy API will be deprecated in the future !
Already now, all new features of ELPA are only available with the new API. Thus, there
is no reason to keep the legacy API arround for too long.

The release ELPA 2018.11.001 was the last release, where the legacy API has been
enabled by default (and can be disabled at build time).
With release ELPA 2019.05.001 the legacy API is disabled by default, however,
can be still switched on at build time.
With the release ELPA 2019.11.001 the legacy API has been deprecated and support has been droped.

[![Build
status](https://gitlab.mpcdf.mpg.de/elpa/elpa/badges/master/build.svg)](https://gitlab.mpcdf.mpg.de/elpa/elpa/commits/master)

[![Code
coverage](https://gitlab.mpcdf.mpg.de/elpa/badges/master/coverage.svg)](http://elpa.pages.mpcdf.de/elpa/coverage_summary)

![License LGPL v3][license-badge]

[license-badge]: https://img.shields.io/badge/License-LGPL%20v3-blue.svg


## About *ELPA* ##

The computation of selected or all eigenvalues and eigenvectors of a symmetric
(Hermitian) matrix has high relevance for various scientific disciplines.
For the calculation of a significant part of the eigensystem typically direct
eigensolvers are used. For large problems, the eigensystem calculations with
existing solvers can become the computational bottleneck.

As a consequence, the *ELPA* project was initiated with the aim to develop and
implement an efficient eigenvalue solver for petaflop applications, supported
by the German Federal Government, through BMBF Grant 01IH08007, from
Dec 2008 to Nov 2011.

The challenging task has been addressed through a multi-disciplinary consortium
of partners with complementary skills in different areas.

The *ELPA* library was originally created by the *ELPA* consortium,
consisting of the following organizations:

- Max Planck Computing and Data Facility (MPCDF), fomerly known as
  Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
- Bergische Universität Wuppertal, Lehrstuhl für angewandte
  Informatik,
- Technische Universität München, Lehrstuhl für Informatik mit
  Schwerpunkt Wissenschaftliches Rechnen ,
- Fritz-Haber-Institut, Berlin, Abt. Theorie,
- Max-Plack-Institut für Mathematik in den Naturwissenschaften,
  Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
  and
- IBM Deutschland GmbH

*ELPA* is distributed under the terms of version 3 of the license of the
GNU Lesser General Public License as published by the Free Software Foundation.

## Obtaining *ELPA*

There exist several ways to obtain the *ELPA* library either as sources or pre-compiled packages:

- official release tar-gz sources from the *ELPA* [webpage](https://elpa.mpcdf.mpg.de/elpa-tar-archive)
- from the *ELPA* [git repository](https://gitlab.mpcdf.mpg.de/elpa/elpa)
- as packaged software for several Linux distributions (e.g. Debian, Fedora, OpenSuse)

## Terms of usage

Your are free to obtain and use the *ELPA* library, as long as you respect the terms
of version 3 of the license of the GNU Lesser General Public License.

No other conditions have to be met.

Nonetheless, we are grateful if you cite the following publications:

  If you use ELPA in general:

  T. Auckenthaler, V. Blum, H.-J. Bungartz, T. Huckle, R. Johanni,
  L. Krämer, B. Lang, H. Lederer, and P. R. Willems,
  "Parallel solution of partial symmetric eigenvalue problems from
  electronic structure calculations",
  Parallel Computing 37, 783-794 (2011).
  doi:10.1016/j.parco.2011.05.002.

  Marek, A.; Blum, V.; Johanni, R.; Havu, V.; Lang, B.; Auckenthaler,
  T.; Heinecke, A.; Bungartz, H.-J.; Lederer, H.
  "The ELPA library: scalable parallel eigenvalue solutions for electronic
  structure theory and computational science",
  Journal of Physics Condensed Matter, 26 (2014)
  doi:10.1088/0953-8984/26/21/213201

  If you use the GPU version of ELPA:

  Kus, P; Marek, A.; Lederer, H.
  "GPU Optimization of Large-Scale Eigenvalue Solver",
  In: Radu F., Kumar K., Berre I., Nordbotten J., Pop I. (eds)
  Numerical Mathematics and Advanced Applications ENUMATH 2017. ENUMATH 2017.
  Lecture Notes in Computational Science and Engineering, vol 126. Springer, Cham

  Yu, V.; Moussa, J.; Kus, P.; Marek, A.; Messmer, P.; Yoon, M.; Lederer, H.; Blum, V.
  "GPU-Acceleration of the ELPA2 Distributed Eigensolver for Dense Symmetric and Hermitian Eigenproblems",
  https://arxiv.org/abs/2002.10991

  If you use the new API and/or autotuning:

  Kus, P.; Marek, A.; Koecher, S. S.; Kowalski H.-H.; Carbogno, Ch.; Scheurer, Ch.; Reuter, K.; Scheffler, M.; Lederer, H.
  "Optimizations of the Eigenvaluesolvers in the ELPA Library",
  Parallel Computing 85, 167-177 (2019)

  If you use the new support for skew-symmetric matrices:
  Benner, P.; Draxl, C.; Marek, A.; Penke C.; Vorwerk, C.;
  "High Performance Solution of Skew-symmetric Eigenvalue Problems with Applications in Solving the Bethe-Salpeter Eigenvalue Problem",
  https://arxiv.org/abs/1912.04062, submitted to Parallel Computing


## Installation of the *ELPA* library

*ELPA* is shipped with a standard autotools automake installation infrastructure.
Some other libraries are needed to install *ELPA* (the details depend on how you
configure *ELPA*):

  - Basic Linear Algebra Subroutines (BLAS)
  - Lapack routines
  - Basic Linear Algebra Communication Subroutines (BLACS)
  - Scalapack routines
  - a working MPI library

Please refer to the [INSTALL document](INSTALL.md) on details of the installation process and
the possible configure options.

## Using *ELPA*

Please have a look at the [USERS_GUIDE](USERS_GUIDE.md) file, to get a documentation or at the [online](http://elpa.mpcdf.mpg.de/html/Documentation/ELPA-2020.05.001.rc1/html/index.html) doxygen documentation, where you find the definition of the interfaces.

## Contributing to *ELPA*

It has been, and is, a tremendous effort to develop and maintain the
*ELPA* library. A lot of things can still be done, but our man-power is limited.

Thus every effort and help to improve the *ELPA* library is highly appreciated.
For details please see the [CONTRIBUTING](CONTRIBUTING.md) document.


PEXSI: Pole EXpansion and Selected Inversion
============================================

The Pole EXpansion and Selected Inversion (PEXSI) method is a fast method for
electronic structure calculation based on Kohn-Sham density functional theory.
It efficiently evaluates certain selected elements of matrix functions, e.g.,
the Fermi-Dirac function of the KS Hamiltonian, which yields a density matrix.
It can be used as an alternative to diagonalization methods for obtaining the
density, energy and forces in electronic structure calculations. The PEXSI
library is written in C++, and uses message passing interface (MPI) to
parallelize the computation on distributed memory computing systems and achieve
scalability on more than 10,000 processors.

From numerical linear algebra perspective, the PEXSI library can be used as a
general tool for evaluating certain selected elements of a matrix function, and
therefore has application beyond electronic structure calculation as well.

The documentation of PEXSI is compiled by Sphinx hosted on

https://pexsi.readthedocs.io

For installation instructions please (be patient for dependecies) see

https://pexsi.readthedocs.io/en/latest/install.html
# [Eigenvalue SoLvers for Petaflop-Applications (ELPA)](http://elpa.mpcdf.mpg.de)

## Current Release ##

The current release is ELPA 2020.05.001.rc1 The current supported API version
is 20190501. This release supports the earliest API version 20170403.

The old, obsolete legacy API will be deprecated in the future !
Already now, all new features of ELPA are only available with the new API. Thus, there
is no reason to keep the legacy API arround for too long.

The release ELPA 2018.11.001 was the last release, where the legacy API has been
enabled by default (and can be disabled at build time).
With release ELPA 2019.05.001 the legacy API is disabled by default, however,
can be still switched on at build time.
With the release ELPA 2019.11.001 the legacy API has been deprecated and support has been droped.

[![Build
status](https://gitlab.mpcdf.mpg.de/elpa/elpa/badges/master/build.svg)](https://gitlab.mpcdf.mpg.de/elpa/elpa/commits/master)

[![Code
coverage](https://gitlab.mpcdf.mpg.de/elpa/badges/master/coverage.svg)](http://elpa.pages.mpcdf.de/elpa/coverage_summary)

![License LGPL v3][license-badge]

[license-badge]: https://img.shields.io/badge/License-LGPL%20v3-blue.svg


## About *ELPA* ##

The computation of selected or all eigenvalues and eigenvectors of a symmetric
(Hermitian) matrix has high relevance for various scientific disciplines.
For the calculation of a significant part of the eigensystem typically direct
eigensolvers are used. For large problems, the eigensystem calculations with
existing solvers can become the computational bottleneck.

As a consequence, the *ELPA* project was initiated with the aim to develop and
implement an efficient eigenvalue solver for petaflop applications, supported
by the German Federal Government, through BMBF Grant 01IH08007, from
Dec 2008 to Nov 2011.

The challenging task has been addressed through a multi-disciplinary consortium
of partners with complementary skills in different areas.

The *ELPA* library was originally created by the *ELPA* consortium,
consisting of the following organizations:

- Max Planck Computing and Data Facility (MPCDF), fomerly known as
  Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
- Bergische Universität Wuppertal, Lehrstuhl für angewandte
  Informatik,
- Technische Universität München, Lehrstuhl für Informatik mit
  Schwerpunkt Wissenschaftliches Rechnen ,
- Fritz-Haber-Institut, Berlin, Abt. Theorie,
- Max-Plack-Institut für Mathematik in den Naturwissenschaften,
  Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
  and
- IBM Deutschland GmbH

*ELPA* is distributed under the terms of version 3 of the license of the
GNU Lesser General Public License as published by the Free Software Foundation.

## Obtaining *ELPA*

There exist several ways to obtain the *ELPA* library either as sources or pre-compiled packages:

- official release tar-gz sources from the *ELPA* [webpage](https://elpa.mpcdf.mpg.de/elpa-tar-archive)
- from the *ELPA* [git repository](https://gitlab.mpcdf.mpg.de/elpa/elpa)
- as packaged software for several Linux distributions (e.g. Debian, Fedora, OpenSuse)

## Terms of usage

Your are free to obtain and use the *ELPA* library, as long as you respect the terms
of version 3 of the license of the GNU Lesser General Public License.

No other conditions have to be met.

Nonetheless, we are grateful if you cite the following publications:

  If you use ELPA in general:

  T. Auckenthaler, V. Blum, H.-J. Bungartz, T. Huckle, R. Johanni,
  L. Krämer, B. Lang, H. Lederer, and P. R. Willems,
  "Parallel solution of partial symmetric eigenvalue problems from
  electronic structure calculations",
  Parallel Computing 37, 783-794 (2011).
  doi:10.1016/j.parco.2011.05.002.

  Marek, A.; Blum, V.; Johanni, R.; Havu, V.; Lang, B.; Auckenthaler,
  T.; Heinecke, A.; Bungartz, H.-J.; Lederer, H.
  "The ELPA library: scalable parallel eigenvalue solutions for electronic
  structure theory and computational science",
  Journal of Physics Condensed Matter, 26 (2014)
  doi:10.1088/0953-8984/26/21/213201

  If you use the GPU version of ELPA:

  Kus, P; Marek, A.; Lederer, H.
  "GPU Optimization of Large-Scale Eigenvalue Solver",
  In: Radu F., Kumar K., Berre I., Nordbotten J., Pop I. (eds)
  Numerical Mathematics and Advanced Applications ENUMATH 2017. ENUMATH 2017.
  Lecture Notes in Computational Science and Engineering, vol 126. Springer, Cham

  Yu, V.; Moussa, J.; Kus, P.; Marek, A.; Messmer, P.; Yoon, M.; Lederer, H.; Blum, V.
  "GPU-Acceleration of the ELPA2 Distributed Eigensolver for Dense Symmetric and Hermitian Eigenproblems",
  https://arxiv.org/abs/2002.10991

  If you use the new API and/or autotuning:

  Kus, P.; Marek, A.; Koecher, S. S.; Kowalski H.-H.; Carbogno, Ch.; Scheurer, Ch.; Reuter, K.; Scheffler, M.; Lederer, H.
  "Optimizations of the Eigenvaluesolvers in the ELPA Library",
  Parallel Computing 85, 167-177 (2019)

  If you use the new support for skew-symmetric matrices:
  Benner, P.; Draxl, C.; Marek, A.; Penke C.; Vorwerk, C.;
  "High Performance Solution of Skew-symmetric Eigenvalue Problems with Applications in Solving the Bethe-Salpeter Eigenvalue Problem",
  https://arxiv.org/abs/1912.04062, submitted to Parallel Computing


## Installation of the *ELPA* library

*ELPA* is shipped with a standard autotools automake installation infrastructure.
Some other libraries are needed to install *ELPA* (the details depend on how you
configure *ELPA*):

  - Basic Linear Algebra Subroutines (BLAS)
  - Lapack routines
  - Basic Linear Algebra Communication Subroutines (BLACS)
  - Scalapack routines
  - a working MPI library

Please refer to the [INSTALL document](INSTALL.md) on details of the installation process and
the possible configure options.

## Using *ELPA*

Please have a look at the [USERS_GUIDE](USERS_GUIDE.md) file, to get a documentation or at the [online](http://elpa.mpcdf.mpg.de/html/Documentation/ELPA-2020.05.001.rc1/html/index.html) doxygen documentation, where you find the definition of the interfaces.

## Contributing to *ELPA*

It has been, and is, a tremendous effort to develop and maintain the
*ELPA* library. A lot of things can still be done, but our man-power is limited.

Thus every effort and help to improve the *ELPA* library is highly appreciated.
For details please see the [CONTRIBUTING](CONTRIBUTING.md) document.


# SuperLU_DIST (version 6.2)

[![Build Status](https://travis-ci.org/xiaoyeli/superlu_dist.svg?branch=master)](https://travis-ci.org/xiaoyeli/superlu_dist)
[Nightly tests](http://my.cdash.org/index.php?project=superlu_dist)

SuperLU_DIST contains a set of subroutines to solve a sparse linear system
A*X=B. It uses Gaussian elimination with static pivoting (GESP).
Static pivoting is a technique that combines the numerical stability of
partial pivoting with the scalability of Cholesky (no pivoting),
to run accurately and efficiently on large numbers of processors.

SuperLU_DIST is a parallel extension to the serial SuperLU library.
It is targeted for the distributed memory parallel machines.
SuperLU_DIST is implemented in ANSI C, and MPI for communications.
Currently, the LU factorization and triangular solution routines,
which are the most time-consuming part of the solution process,
are parallelized. The other routines, such as static pivoting and
column preordering for sparsity are performed sequentially.
This "alpha" release contains double-precision real and double-precision
complex data types.

### The distribution contains the following directory structure:

```
SuperLU_DIST/README    instructions on installation
SuperLU_DIST/CBLAS/    needed BLAS routines in C, not necessarily fast
	 	       (NOTE: this version is single threaded. If you use the
		       library with multiple OpenMP threads, performance
		       relies on a good multithreaded BLAS implementation.)
SuperLU_DIST/DOC/      the Users' Guide
SuperLU_DIST/EXAMPLE/  example programs
SuperLU_DIST/INSTALL/  test machine dependent parameters
SuperLU_DIST/SRC/      C source code, to be compiled into libsuperlu_dist.a
SuperLU_DIST/TEST/     testing code
SuperLU_DIST/lib/      contains library archive libsuperlu_dist.a
SuperLU_DIST/Makefile  top-level Makefile that does installation and testing
SuperLU_DIST/make.inc  compiler, compiler flags, library definitions and C
	               preprocessor definitions, included in all Makefiles.
	               (You may need to edit it to suit your system
	               before compiling the whole package.)
SuperLU_DIST/MAKE_INC/ sample machine-specific make.inc files
```

## INSTALLATION

There are two ways to install the package. One requires users to
edit makefile manually, the other uses CMake automatic build system.
The procedures are described below.

### Installation option 1: Using CMake build system.
You will need to create a build tree from which to invoke CMake.

First, in order to use parallel symbolic factorization function, you
need to install ParMETIS parallel ordering package and define the
two environment variables: PARMETIS_ROOT and PARMETIS_BUILD_DIR

```
export PARMETIS_ROOT=<Prefix directory of the ParMETIS installation>
export PARMETIS_BUILD_DIR=${PARMETIS_ROOT}/build/Linux-x86_64
```

Second, in order to use parallel weighted matching AWPM for numerical
pre-pivoting, you need to install CombBLAS and define the environment
variable:

```
export COMBBLAS_ROOT=<Prefix directory of the CombBLAS installation>
export COMBBLAS_BUILD_DIR=${COMBBLAS_ROOT}/_build
```

Once these needed third-party libraries are in place, SuperLU installation
can be done as follows from the top level directory:

For a simple installation with default setting, do:
(ParMETIS is needed, i.e., TPL_ENABLE_PARMETISLIB=ON)
```
mkdir build ; cd build;
cmake .. \
    -DTPL_PARMETIS_INCLUDE_DIRS="${PARMETIS_ROOT}/include;${PARMETIS_ROOT}/metis/include" \
    -DTPL_PARMETIS_LIBRARIES="${PARMETIS_BUILD_DIR}/libparmetis/libparmetis.a;${PARMETIS_BUILD_DIR}/libmetis/libmetis.a" \
```

For a more sophisticated installation including third-part libraries, do:
```
cmake .. \
    -DTPL_PARMETIS_INCLUDE_DIRS="${PARMETIS_ROOT}/include;${PARMETIS_ROOT}/metis/include" \
    -DTPL_PARMETIS_LIBRARIES="${PARMETIS_BUILD_DIR}/libparmetis/libparmetis.a;${PARMETIS_BUILD_DIR}/libmetis/libmetis.a" \
    -DTPL_ENABLE_COMBBLASLIB=ON \
    -DTPL_COMBBLAS_INCLUDE_DIRS="${COMBBLAS_ROOT}/_install/include;${COMBBLAS_R\
OOT}/Applications/BipartiteMatchings" \
    -DTPL_COMBBLAS_LIBRARIES="${COMBBLAS_BUILD_DIR}/libCombBLAS.a" \
    -DCMAKE_C_FLAGS="-std=c99 -g -DPRNTlevel=0 -DDEBUGlevel=0" \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_CXX_FLAGS="-std=c++11" \
    -DTPL_ENABLE_BLASLIB=OFF \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_INSTALL_PREFIX=.

( see example cmake script: run_cmake_build.sh )
```
You can disable LAPACK, ParMetis or CombBLAS with the following cmake option:
`-DTPL_ENABLE_LAPACKLIB=FALSE`
`-DTPL_ENABLE_PARMETISLIB=FALSE`
`-DTPL_ENABLE_COMBBLASLIB=FALSE`

To actually build (compile), type:
`make`

To install the libraries, type:
`make install`

To run the installation test, type:
`ctest`
(The outputs are in file: `build/Testing/Temporary/LastTest.log`)
or,
`ctest -D Experimental`
or,
`ctest -D Nightly`

**NOTE:**
The parallel execution in ctest is invoked by "mpiexec" command which is
from MPICH environment. If your MPI is not MPICH/mpiexec based, the test
execution may fail. You can always go to TEST/ directory to perform
testing manually.

**Note on the C-Fortran name mangling handled by C preprocessor definition:**
In the default setting, we assume that Fortran expects a C routine
to have an underscore postfixed to the name. Depending on the
compiler, you may need to define one of the following flags in
during the cmake build to overwrite default setting:
```
cmake .. -DCMAKE_C_FLAGS="-DNoChange"
cmake .. -DCMAKE_C_FLAGS="-DUpCase"
```


### Installation option 2: Manual installation with makefile.
Before installing the package, please examine the three things dependent
on your system setup:

#### 1.1 Edit the make.inc include file.

This make include file is referenced inside each of the Makefiles
in the various subdirectories. As a result, there is no need to
edit the Makefiles in the subdirectories. All information that is
machine specific has been defined in this include file.

Sample machine-specific make.inc are provided in the MAKE_INC/
directory for several platforms, such as Cray XT5, Linux, Mac-OS, and CUDA.
When you have selected the machine to which you wish to install
SuperLU_DIST, copy the appropriate sample include file
(if one is present) into make.inc.

For example, if you wish to run SuperLU_DIST on a Cray XT5,  you can do
`cp MAKE_INC/make.xt5  make.inc`

For the systems other than listed above, some porting effort is needed
for parallel factorization routines. Please refer to the Users' Guide
for detailed instructions on porting.

The following CPP definitions can be set in CFLAGS.
```
-DXSDK_INDEX_SIZE=64
use 64-bit integers for indexing sparse matrices. (default 32 bit)

-DPRNTlevel=[0,1,2,...]
printing level to show solver's execution details. (default 0)

-DDEBUGlevel=[0,1,2,...]
diagnostic printing level for debugging purpose. (default 0)
```

#### 1.2. The BLAS library.

The parallel routines in SuperLU_DIST use some BLAS routines on each MPI
process. Moreover, if you enable OpenMP with multiple threads, you need to
link with a multithreaded BLAS library. Otherwise performance will be poor.
A good public domain BLAS library is OpenBLAS (http://www.openblas.net),
which has OpenMP support.

If you have a BLAS library your machine, you may define the following in
the file make.inc:
```
BLASDEF = -DUSE_VENDOR_BLAS
BLASLIB = <BLAS library you wish to link with>
```
The CBLAS/ subdirectory contains the part of the C BLAS (single threaded)
needed by SuperLU_DIST package. However, these codes are intended for use
only if there is no faster implementation of the BLAS already
available on your machine. In this case, you should go to the
top-level SuperLU_DIST/ directory and do the following:

1) In make.inc, undefine (comment out) BLASDEF, and define:
` BLASLIB = ../lib/libblas$(PLAT).a`

2) Type: `make blaslib`
to make the BLAS library from the routines in the
` CBLAS/ subdirectory.`

#### 1.3. External libraries.

  ##### 1.3.1 LAPACK.
  Starting Version 6.0, the triangular solve routine can perform explicit
  inversion on the diagonal blocks, using LAPACK's xTRTRI inversion routine.
  To use this feature, you should define the following in make.inc:
```
SLU_HAVE_LAPACK = TRUE
LAPACKLIB = <lapack library you wish to link with>
```
You can disable LAPACK with the following line in SRC/superlu_dist_config.h:
```
#undef SLU_HAVE_LAPACK
```

  ##### 1.3.2 Metis and ParMetis.

If you will use Metis or ParMetis for sparsity ordering, you will
need to install them yourself. Since ParMetis package already
contains the source code for the Metis library, you can just
download and compile ParMetis from:
[http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)

After you have installed it, you should define the following in make.inc:
```
HAVE_PARMETIS = TRUE
METISLIB = -L<metis directory> -lmetis
PARMETISLIB = -L<parmetis directory> -lparmetis
I_PARMETIS = -I<parmetis directory>/include -I<parmetis directory>/metis/include
```
You can disable ParMetis with the following line in SRC/superlu_dist_config.h:
```
#undef HAVE_PARMETIS
```

 ##### 1.3.3 CombBLAS.

You can use parallel approximate weight perfect matching (AWPM) algorithm
to perform numerical pre-pivoting for stability. The default pre-pivoting
is to use MC64 provided internally, which is an exact algorithm, but serial.
In order to use AWPM, you will need to install CombBLAS yourself, at the
download site:
[https://people.eecs.berkeley.edu/~aydin/CombBLAS/html/index.html](https://people.eecs.berkeley.edu/~aydin/CombBLAS/html/index.html)

After you have installed it, you should define the following in make.inc:
```
HAVE_COMBBLAS = TRUE
COMBBLASLIB = <combblas root>/_build/libCombBLAS.a
I_COMBBLAS=-I<combblas root>/_install/include -I<combblas root>/Applications/BipartiteMatchings
```
You can disable CombBLAS with the following line in SRC/superlu_dist_config.h:
```
#undef HAVE_COMBBLAS
```


#### 1.4. C preprocessor definition CDEFS. (Replaced by cmake module FortranCInterface.)

In the header file SRC/Cnames.h, we use macros to determine how
C routines should be named so that they are callable by Fortran.
(Some vendor-supplied BLAS libraries do not have C interfaces. So the
re-naming is needed in order for the SuperLU BLAS calls (in C) to
interface with the Fortran-style BLAS.)
The possible options for CDEFS are:
```
-DAdd_: Fortran expects a C routine to have an underscore
  postfixed to the name;
  (This is set as the default)
-DNoChange: Fortran expects a C routine name to be identical to
      that compiled by C;
-DUpCase: Fortran expects a C routine name to be all uppercase.
```

#### 1.5. Multicore and GPU (optional).

To use OpenMP parallelism, need to link with an OpenMP library, and
set the number of threads you wish to use as follows (bash):

`export OMP_NUM_THREADS=<##>`

To enable NVIDIA GPU access, need to take the following 2 step:
1) Set the following Linux environment variable:
`export ACC=GPU`

2) Add the CUDA library location in make.inc:
```
ifeq "${ACC}" "GPU"
CFLAGS += -DGPU_ACC
INCS += -I<CUDA directory>/include
LIBS += -L<CUDA directory>/lib64 -lcublas -lcudart
endif
```
A Makefile is provided in each subdirectory. The installation can be done
completely automatically by simply typing "make" at the top level.



## Windows Usage
Prerequisites: CMake, Visual Studio, Microsoft HPC Pack
This has been tested with Visual Studio 2017, without Parmetis,
without Fortran, and with OpenMP disabled.

The cmake configuration line used was
```
'/winsame/contrib-vs2017/cmake-3.9.4-ser/bin/cmake' \
  -DCMAKE_INSTALL_PREFIX:PATH=C:/winsame/volatile-vs2017/superlu_dist-master.r147-parcomm \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_COLOR_MAKEFILE:BOOL=FALSE \
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
  -Denable_openmp:BOOL=FALSE \
  -DCMAKE_C_COMPILER:FILEPATH='C:/Program Files (x86)/Microsoft Visual Studio/2017/Professional/VC/Tools/MSVC/14.11.25503/bin/HostX64/x64/cl.exe' \
  -DCMAKE_C_FLAGS:STRING='/DWIN32 /D_WINDOWS /W3' \
  -DTPL_ENABLE_PARMETISLIB:BOOL=FALSE \
  -DXSDK_ENABLE_Fortran=OFF \
  -G 'NMake Makefiles JOM' \
  C:/path/to/superlu_dist
```

After configuring, simply do
```
  jom # or nmake
  jom install  # or nmake install
```

Libraries will be installed under
C:/winsame/volatile-vs2017/superlu_dist-master.r147-parcomm/lib
for the above configuration.

If you wish to test:
  `ctest`

## READING SPARSE MATRIX FILES

The SRC/ directory contains the following routines to read different file
formats, they all have the similar calling sequence.
```
$ ls -l dread*.c
dreadMM.c              : Matrix Market, files with suffix .mtx
dreadhb.c              : Harrell-Boeing, files with suffix .rua
dreadrb.c              : Rutherford-Boeing, files with suffix .rb
dreadtriple.c          : triplet, with header
dreadtriple_noheader.c : triplet, no header, which is also readable in Matlab
```

## REFERENCES

**[1]** X.S. Li and J.W. Demmel, "SuperLU_DIST: A Scalable Distributed-Memory
 Sparse Direct Solver for Unsymmetric Linear Systems", ACM Trans. on Math.
 Software, Vol. 29, No. 2, June 2003, pp. 110-140.
**[2]** L. Grigori, J. Demmel and X.S. Li, "Parallel Symbolic Factorization
 for Sparse LU with Static Pivoting", SIAM J. Sci. Comp., Vol. 29, Issue 3,
 1289-1314, 2007.
**[3]** P. Sao, R. Vuduc and X.S. Li, "A distributed CPU-GPU sparse direct
 solver", Proc. of EuroPar-2014 Parallel Processing, August 25-29, 2014.
 Porto, Portugal.
**[4]** P. Sao, X.S. Li, R. Vuduc, “A Communication-Avoiding 3D Factorization
 for Sparse Matrices”, Proc. of IPDPS, May 21–25, 2018, Vancouver.
**[5]** Y. Liu, M. Jacquelin, P. Ghysels and X.S. Li, “Highly scalable
 distributed-memory sparse triangular solution algorithms”, Proc. of
 SIAM workshop on Combinatorial Scientific Computing, June 6-8, 2018,
 Bergen, Norway.

**Xiaoye S. Li**, Lawrence Berkeley National Lab, [xsli@lbl.gov](xsli@lbl.gov)
**Gustavo Chavez**, Lawrence Berkeley National Lab, [gichavez@lbl.gov](gichavez@lbl.gov)
**Laura Grigori**, INRIA, France, [laura.grigori@inria.fr](laura.grigori@inria.fr)
**Yang Liu**, Lawrence Berkeley National Lab, [liuyangzhuan@lbl.gov](liuyangzhuan@lbl.gov)
**Meiyue Shao**, Lawrence Berkeley National Lab, [myshao@lbl.gov](myshao@lbl.gov)
**Piyush Sao**, Georgia Institute of Technology, [piyush.feynman@gmail.com](piyush.feynman@gmail.com)
**Ichitaro Yamazaki**, Univ. of Tennessee, [ic.yamazaki@gmail.com](ic.yamazaki@gmail.com)
**Jim Demmel**, UC Berkeley, [demmel@cs.berkeley.edu](demmel@cs.berkeley.edu)
**John Gilbert**, UC Santa Barbara, [gilbert@cs.ucsb.edu](gilbert@cs.ucsb.edu)



## RELEASE VERSIONS
```
October 15, 2003    Version 2.0
October 1,  2007    Version 2.1
Feburary 20, 2008   Version 2.2
October 15, 2008    Version 2.3
June 9, 2010        Version 2.4
November 23, 2010   Version 2.5
March 31, 2013      Version 3.3
October 1, 2014     Version 4.0
July 15, 2014       Version 4.1
September 25, 2015  Version 4.2
December 31, 2015   Version 4.3
April 8, 2016       Version 5.0.0
May 15, 2016        Version 5.1.0
October 4, 2016     Version 5.1.1
December 31, 2016   Version 5.1.3
September 30, 2017  Version 5.2.0
January 28, 2018    Version 5.3.0
June 1, 2018        Version 5.4.0
```
Welcome to the API documentation of ELSI {#mainpage}
========================================

ELSI provides and enhances scalable, open-source software library solutions for electronic structure calculations in materials science, condensed matter physics, chemistry, and many other fields. ELSI focuses on methods that solve or circumvent eigenproblems in electronic structure theory. For more information, please visit the [ELSI Interchange](https://elsi-interchange.org) website, or contact us via [email](mailto:elsi-team@duke.edu).

The eigensolvers and density matrix solvers supported in ELSI include:

* Shared-memory eigensolvers

    * [LAPACK](https://www.netlib.org/lapack): One-stage and two-stage tridiagonalization-based dense eigensolvers.

    * [MAGMA](https://icl.utk.edu/magma): GPU-accelerated one-stage and two-stage tridiagonalization-based dense eigensolvers.

* Distributed-memory eigensolvers

    * [ELPA](https://elpa.mpcdf.mpg.de): One-stage and two-stage tridiagonalization-based dense eigensolvers.

    * [EigenEXA](https://www.r-ccs.riken.jp/labs/lpnctrt/en/projects/eigenexa): One-stage tridiagonalization-based and pentadiagonalization-based dense eigensolvers.

    * [SLEPc](https://slepc.upv.es): Sparse eigensolver based on parallel spectrum slicing.

* Distributed-memory density matrix solvers

    * [libOMM](https://esl.cecam.org/LibOMM): Orbital minimization method based on dense linear algebra.

    * [PEXSI](https://pexsi.org): Pole expansion and selected inversion based on sparse linear algebra.

    * [NTPoly](https://william-dawson.github.io/NTPoly): Linear scaling density matrix purification based on sparse linear algebra.

The design of ELSI focues on portability from laptop-type computers all the way up to the most efficient massively parallel supercomputers and new architectures. Work is in progress to support additional solver libraries, providing electronic structure code developers and users with a flexible, customizable choice of solution for the central algebraic problems in large-scale electronic structure simulations.

The ELSI software is adopted by electronic structure code projects such as [DFTB+](https://www.dftbplus.org), DGDFT, [FHI-aims](https://aimsclub.fhi-berlin.mpg.de), and [SIESTA](https://departments.icmab.es/leem/siesta).
elsipy provide a general parallel python interface to ELSI via mpi4py.

Steps to install elsipy, I use miniconda as package manager

1) conda create -n elsi_env python=3.9
2) conda activate elsi_env
3) pip install numpy scipy cython nose
4) conda install mpich scalapack mpi4py
5) After modify the paths in setup.py; python setup.py install

An example is given in test.py
