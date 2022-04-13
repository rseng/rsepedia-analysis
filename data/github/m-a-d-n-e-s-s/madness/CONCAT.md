madness
=======

Multiresolution Adaptive Numerical Environment for Scientific Simulation

# Summary

MADNESS provides a high-level environment for the solution of integral and differential equations in many dimensions using adaptive, fast methods with guaranteed precision based on multi-resolution analysis and novel separated representations. There are three main components to MADNESS. At the lowest level is a new petascale parallel programming environment that increases programmer productivity and code performance/scalability while maintaining backward compatibility with current programming tools such as MPI and Global Arrays. The numerical capabilities built upon the parallel tools provide a high-level environment for composing and solving numerical problems in many (1-6+) dimensions. Finally, built upon the numerical tools are new applications with initial focus upon chemistry, atomic and molecular physics, material science, and nuclear structure.

Please look in the [wiki](https://github.com/m-a-d-n-e-s-s/madness/wiki) for more information and project activity.

Here's a [video](http://www.youtube.com/watch?v=dBwWjmf5Tic) about MADNESS.

# Funding
The developers gratefully acknowledge the support of the Department of Energy, Office of Science, Office of Basic Energy Sciences and Office of Advanced Scientific Computing Research, under contract DE-AC05-00OR22725 with Oak Ridge National Laboratory.

The developers gratefully acknowledge the support of the National Science Foundation under grant 0509410 to the University of Tennessee in collaboration with The Ohio State University (P. Sadayappan). The MADNESS parallel runtime and parallel tree-algorithms include concepts and software developed under this project.

The developers gratefully acknowledge the support of the National Science Foundation under grant NSF OCI-0904972 to the University of Tennessee. The solid state physics and multiconfiguration SCF capabilities are being developed by this project.

The developers gratefully acknowledge the support of the National Science Foundation under grant NSF CHE-0625598 to the University of Tennessee, in collaboration with UIUC/NCSA. Some of the multi-threading and preliminary GPGPU ports were developed by this project.

The developers gratefully acknowledge the support of the Defense Advanced Research Projects Agency (DARPA) under subcontract from Argonne National Laboratory as part of the High-Productivity Computer Systems (HPCS) language evaluation project.

<p align="left" style="padding: 20px">
<img src="http://libelemental.org/_static/elemental.png">
</p>

**Elemental** is a modern C++ library for distributed-memory dense linear 
algebra.
The library was initially released in
[Elemental: A new framework for distributed memory dense linear algebra](https://dl.acm.org/citation.cfm?doid=2427023.2427030) 
and is the key building block for the distributed-memory sparse-direct solver 
[Clique](http://www.github.com/poulson/Clique.git).

Please visit [the download page](http://libelemental.org/download/) for 
download instructions.

### Documentation

The [documentation for the development version of Elemental](http://libelemental.org/documentation) is built using [Sphinx](http://sphinx.pocoo.org).

### Related open-source packages

Implementations:

1. [DPLASMA](http://icl.eecs.utk.edu/dplasma/)
2. [PLAPACK](http://www.cs.utexas.edu/~plapack)
3. [ScaLAPACK](http://www.netlib.org/scalapack) (and the add-on, [ELPA](http://elpa.rzg.mpg.de/))

Wrappers:

1. [PETSc](https://www.mcs.anl.gov/petsc/)
2. [Trilinos](http://trilinos.sandia.gov)

Note that [PETSc](https://www.mcs.anl.gov/petsc/) contains interfaces for both 
[Elemental](http://github.com/elemental/Elemental.git) and
[Clique](http://github.com/poulson/Clique.git).

### Elemental's root directory

This is the root directory of the entire project, and it contains:

-  `AUTHORS`: the list of source code contributors
-  `cmake/`: auxiliary files for CMake configuration
-  `CMakeLists.txt`: the CMake configuration file
-  `doc/`: Sphinx documentation 
-  `examples/`: various concise examples of Elemental's functionality
-  `experimental/`: experimental code which is not yet library quality
-  `external/`: non-standard external code which Elemental builds on top of
-  `include/`: Elemental's header files; most of the library resides here
-  `LICENSE`: the New BSD License file
-  `octave/`: pedagogical versions of algorithms used in Elemental
-  `PUBLICATIONS`: publications directly related to this source code
-  `README.md`: this file
-  `REFERENCES`: some publications referenced in the creation of this software
-  `src/`: Elemental's source files; a small portion of the library is here
-  `tests/`: programs meant to test the accuracy of Elemental
-  `TODO`: a list of near/long-term goals of the project
-  `vagrant/`: lightweight preconfigured virtual-machines for Elemental
### `cmake/`

This folder contains Elemental's auxiliary files for its CMake configuration
system:

-  `config.h.cmake`: The preprocessor definitions passed to Elemental
-  `ElemSub.cmake`: A CMake include file for simplifying the usage of
   Elemental as a CMake subproject
-  `ElemVars.cmake`: A Makefile include meant to simplify the usage of
   Elemental as a dependency in a project directly using 'make'
-  `tests/`: Various configuration tests
-  `toolchains/`: Toolchain files for various specialized architectures 
   (e.g., Blue Gene/Q)
### `cmake/toolchains/`

This folder contains a variety of toolchain files (several user-contributed)
which are in different states of repair. They are provided in the hopes that
they might be useful but it is likely that at least one is non-functional.
This folder will slowly accumulate various [Vagrantfiles](http://vagrantup.com) 
for developing Elemental. For instance, to load a 32-bit development environment
for Elemental within Ubuntu 13.10 (Saucy), simply run: 
``cd saucy32; vagrant up``.
### `include/`

This is the root of the include folder, which is where most of the library is 
contained. In addition to this README, this folder holds:

-  `elemental/`: the supporting header files
-  `elemental.hpp`: the catch-all header file (users should start by just 
   including this)
-  `elemental-lite.hpp`: the minimal header file, which can be used in 
   combination with manually including header files from `elemental/` in order
   to decrease build times
### `include/elemental/`

-  `blas-like/`: BLAS-like functionality (e.g., GEMM)
-  `control/`: solvers for control theory (e.g., Sylvester equations)
-  `convex/`: convex optimization routines (e.g., SVT)
-  `core/`: core data structures (e.g., `DistMatrix`)
-  `io/`: input/output functionality (e.g., printing and visualization)
-  `lapack-like/`: LAPACK-like functionality (e.g., SVD)
-  `matrices/`: special matrices (e.g., uniform random, Wilkinson, etc.)
### `include/elemental/convex/`

- `decl.hpp`: include all declarations of convex opt. related routines
- `impl.hpp`: include all implementations of convex opt. related routines

This folder contains a few utilities for convex optimization; the only 
nontrivial one is for Singular Value soft-Thresholding (SVT).

-  `LogBarrier.hpp`: negative log of the determinant of an HPD matrix
-  `LogDetDiv.hpp`: divergence between two HPD matrices
-  `SoftThreshold.hpp`: soft-threshold each entry of a matrix
-  `SVT.hpp`: soft-threshold the singular values of a matrix
### `include/elemental/core/`

#### Axpy interface

-  `AxpyInterface/`:
-  `AxpyInterface.hpp`:

#### Complex arithmetic 

-  `Complex/`:
-  `Complex.hpp`:

#### `BlockDistMatrix`

-  `BlockDistMatrix/`:
-  `BlockDistMatrix.hpp`:

#### `DistMatrix`

-  `DistMatrix/`:
-  `DistMatrix.hpp`:

#### Basic environment

-  `environment/`:

#### Process grid

-  `Grid/`:

#### Imported functionality

-  `imports/`:

#### Indexing utilities

-  `indexing/`:

#### `Matrix`

-  `Matrix/`:
-  `Matrix.hpp`:

#### Memory

-  `Memory/`:

#### Random number generation

-  `random/`:

#### Matrix view manipulation

-  `views/`:

#### Timings

-  `Timer/`:

#### Datatypes

-  `types/`:
### `include/elemental/core/dist_matrix/`

This folder contains the header files for the various partial specializations of
the `DistMatrix` class; please see `src/core/dist_matrix` for the 
corresponding source files. Each specialization involves choosing a 
sensical pairing of distributions for the rows and columns of the matrix:

-  `CIRC`/"o": Only give the data to a single process
-  `STAR`/"\*": Give the data to every process
-  `MC`: Distribute round-robin within each column of the 2D process grid (*M*atrix *C*olumn)
-  `MR`: Distribute round-robin within each row of the 2D process grid (*M*atrix *R*ow)
-  `VC`: Distribute round-robin within a column-major ordering of the entire 
   2D process grid (*V*ector *C*olumn)
-  `VR`: Distribute round-robin within a row-major ordering of the entire
   2D process grid (*V*ector *R*ow)
-  `MD`: Distribute round-robin over a diagonal of the tiling of the 2D process
   grid (*M*atrix *D*iagonal)

The valid pairings are:

| Distribution | ColComm | RowComm | DistComm  | RedundantComm | CrossComm |
|:------------:|:-------:|:-------:|:---------:|:-------------:|:---------:|
| `(o ,o )`    | self    | self    | self      | self          | `VC`      |
| `(* ,* )`    | self    | self    | self      | `VC`          | self      |
| `(MD,* )`    | `MD`    | self    | `MD`      | self          | `MDPerp`  |
| `(* ,MD)`    | self    | `MD`    | `MD`      | self          | `MDPerp`  |
| `(MC,MR)`    | `MC`    | `MR`    | `VC`      | self          | self      |
| `(MR,MC)`    | `MR`    | `MC`    | `VR`      | self          | self      |
| `(MC,* )`    | `MC`    | self    | `MC`      | `MR`          | self      |
| `(* ,MC)`    | self    | `MC`    | `MC`      | `MR`          | self      |
| `(MR,* )`    | `MR`    | self    | `MR`      | `MC`          | self      |
| `(* ,MR)`    | self    | `MR`    | `MR`      | `MC`          | self      |
| `(VC,* )`    | `VC`    | self    | `VC`      | self          | self      |
| `(* ,VC)`    | self    | `VC`    | `VC`      | self          | self      |
| `(VR,* )`    | `VR`    | self    | `VR`      | self          | self      |
| `(* ,VR)`    | self    | `VR`    | `VR`      | self          | self      |

where `DistComm` refers to the communicator that the entire matrix (rather than
just the rows or columns) is distributed over. When the matrix is distributed
over a communicator which only involves only a subset of the processes, it is
possible to either assign the data to just that subset or redundantly store 
the entire matrix on each such subset of processes (e.g., within each row of a 
2D arrangement of the set of processes). The `RedundantComm` refers to the 
communicator where each member process stores the same information, and the 
`CrossComm` is the communicator where only a single process (the *root*) is 
assigned any data.

To make this discussion more precise, each valid matrix distribution for 
`DistMatrix` logically arranges the set of `p` processes of the `r` by `c` 
process grid into a 4D mesh: `ColComm` x `RowComm` x `RedundantComm` x `CrossComm`, where `DistComm` is equal to `ColComm` x `RowComm`.

We are now ready to describe the contents of this folder (in addition to this
file):

-  `Abstract.hpp`: The underlying distribution-agnostic base class
-  `CIRC_CIRC.hpp`: The `<T,CIRC,CIRC>` specialization, which provides a
   distributed matrix where only one process owns data. It provides a simple
   mechanism for forming a matrix on a single process and then redistributing
   into another distribution, e.g., `(MC,MR)`.
-  `MC_MR.hpp`: The standard matrix distribution
-  `MC_STAR.hpp`: Only distribute each column like a standard matrix
   distribution
-  `MD_STAR.hpp`: Distribute each column like the diagonal of the standard
   matrix distribution
-  `MR_MC.hpp`: The transpose of the standard matrix distribution
-  `MR_STAR.hpp`: Distribute each column like the row of a standard matrix
   distribution
-  `STAR_MC.hpp`: Distribute each row like a column of the standard matrix
   distribution
-  `STAR_MD.hpp`: Distribute each row like the diagonal of a standard matrix
   distribution
-  `STAR_MR.hpp`: Distribute each row like a standard matrix distribution
-  `STAR_STAR.hpp`: Give each process a full copy of the matrix
-  `STAR_VC.hpp`: Distribute each row using a round-robin wrapping over a
   column-major ordering of the process grid
-  `STAR_VR.hpp`: Distribute each row using a round-robin wrapping over a
   row-major ordering of the process grid
-  `VC_STAR.hpp`: Distribute each column using a round-robin wrapping over a
   column-major ordering of the process grid
-  `VR_STAR.hpp`: Distribute each column using a round-robin wrapping over a
   row-major ordering of the process grid
### `include/elemental/core/dist_matrix/`

This folder contains the header files for the various partial specializations of
the `DistMatrix` class; please see `src/core/dist_matrix` for the 
corresponding source files. Each specialization involves choosing a 
sensical pairing of distributions for the rows and columns of the matrix:

-  `CIRC`/"o": Only give the data to a single process
-  `STAR`/"\*": Give the data to every process
-  `MC`: Distribute round-robin within each column of the 2D process grid (*M*atrix *C*olumn)
-  `MR`: Distribute round-robin within each row of the 2D process grid (*M*atrix *R*ow)
-  `VC`: Distribute round-robin within a column-major ordering of the entire 
   2D process grid (*V*ector *C*olumn)
-  `VR`: Distribute round-robin within a row-major ordering of the entire
   2D process grid (*V*ector *R*ow)
-  `MD`: Distribute round-robin over a diagonal of the tiling of the 2D process
   grid (*M*atrix *D*iagonal)

The valid pairings are:

| Distribution | ColComm | RowComm | DistComm  | RedundantComm | CrossComm |
|:------------:|:-------:|:-------:|:---------:|:-------------:|:---------:|
| `(o ,o )`    | self    | self    | self      | self          | `VC`      |
| `(* ,* )`    | self    | self    | self      | `VC`          | self      |
| `(MD,* )`    | `MD`    | self    | `MD`      | self          | `MDPerp`  |
| `(* ,MD)`    | self    | `MD`    | `MD`      | self          | `MDPerp`  |
| `(MC,MR)`    | `MC`    | `MR`    | `VC`      | self          | self      |
| `(MR,MC)`    | `MR`    | `MC`    | `VR`      | self          | self      |
| `(MC,* )`    | `MC`    | self    | `MC`      | `MR`          | self      |
| `(* ,MC)`    | self    | `MC`    | `MC`      | `MR`          | self      |
| `(MR,* )`    | `MR`    | self    | `MR`      | `MC`          | self      |
| `(* ,MR)`    | self    | `MR`    | `MR`      | `MC`          | self      |
| `(VC,* )`    | `VC`    | self    | `VC`      | self          | self      |
| `(* ,VC)`    | self    | `VC`    | `VC`      | self          | self      |
| `(VR,* )`    | `VR`    | self    | `VR`      | self          | self      |
| `(* ,VR)`    | self    | `VR`    | `VR`      | self          | self      |

where `DistComm` refers to the communicator that the entire matrix (rather than
just the rows or columns) is distributed over. When the matrix is distributed
over a communicator which only involves only a subset of the processes, it is
possible to either assign the data to just that subset or redundantly store 
the entire matrix on each such subset of processes (e.g., within each row of a 
2D arrangement of the set of processes). The `RedundantComm` refers to the 
communicator where each member process stores the same information, and the 
`CrossComm` is the communicator where only a single process (the *root*) is 
assigned any data.

To make this discussion more precise, each valid matrix distribution for 
`DistMatrix` logically arranges the set of `p` processes of the `r` by `c` 
process grid into a 4D mesh: `ColComm` x `RowComm` x `RedundantComm` x `CrossComm`, where `DistComm` is equal to `ColComm` x `RowComm`.

We are now ready to describe the contents of this folder (in addition to this
file):

-  `Abstract.hpp`: The underlying distribution-agnostic base class
-  `CIRC_CIRC.hpp`: The `<T,CIRC,CIRC>` specialization, which provides a
   distributed matrix where only one process owns data. It provides a simple
   mechanism for forming a matrix on a single process and then redistributing
   into another distribution, e.g., `(MC,MR)`.
-  `MC_MR.hpp`: The standard matrix distribution
-  `MC_STAR.hpp`: Only distribute each column like a standard matrix 
   distribution
-  `MD_STAR.hpp`: Distribute each column like the diagonal of the standard
   matrix distribution
-  `MR_MC.hpp`: The transpose of the standard matrix distribution
-  `MR_STAR.hpp`: Distribute each column like the row of a standard matrix 
   distribution
-  `STAR_MC.hpp`: Distribute each row like a column of the standard matrix 
   distribution
-  `STAR_MD.hpp`: Distribute each row like the diagonal of a standard matrix
   distribution
-  `STAR_MR.hpp`: Distribute each row like a standard matrix distribution
-  `STAR_STAR.hpp`: Give each process a full copy of the matrix
-  `STAR_VC.hpp`: Distribute each row using a round-robin wrapping over a 
   column-major ordering of the process grid
-  `STAR_VR.hpp`: Distribute each row using a round-robin wrapping over a 
   row-major ordering of the process grid
-  `VC_STAR.hpp`: Distribute each column using a round-robin wrapping over a
   column-major ordering of the process grid
-  `VR_STAR.hpp`: Distribute each column using a round-robin wrapping over a 
   row-major ordering of the process grid
### `include/elemental/core/imports/`

This folder contains the header files for Elemental's wrappers
for external libraries, such as BLAS, LAPACK, and PMRRR. Please see 
`src/core/imports` for the corresponding source files.

In addition to this file, the folder contains:

-  `blas.hpp`: wrappers for Basic Linear Algebra Subprograms (BLAS)
-  `choice.hpp`: command-line processing for sequential programs
-  `flame.hpp`: wrappers for FLAME's QR-based bidiagonal SVD
-  `lapack.hpp`: wrappers for Linear Algebra PACKage (LAPACK)
-  `mpi.hpp`: wrappers for the Message Passing Interface (MPI)
-  `mpi_choice.hpp`: command-line processing for MPI-based programs
-  `pmrrr.hpp`: wrappers for Parallel Multiple Relatively Robust Representations
   (PMRRR)
### `include/elemental/blas-like/`

- `decl.hpp`: include all declarations of BLAS-like routines
- `impl.hpp`: include all implementations of BLAS-like routines

Elemental's BLAS-like functionality is categorized into the typical 'levels':

-  level 1: operations with essentially no data reuse (e.g., DOT)
-  level 2: matrix/vector-like operations (e.g., GEMV)
-  level 3: matrix/matrix-like operations (e.g., GEMM)
### `include/elemental/blas-like/level1/`

The level-1 BLAS-like routines implemented in this folder are:

-  `Adjoint.hpp`: form the adjoint (conjugate-transpose) of a matrix
-  `Axpy.hpp`: Y becomes alpha X + Y (AXPY)
-  `AxpyTriangle.hpp`: Axpy a triangular matrix onto a triangular matrix
-  `Conjugate.hpp`: form the conjugate of a matrix
-  `Copy.hpp`: form a copy of a matrix
-  `DiagonalScale.hpp`: apply a diagonal matrix
-  `DiagonalSolve.hpp`: apply the inverse of a diagonal matrix
-  `Dot.hpp`: form the dot product of two vectors
-  `Dotu.hpp`: form the unconjugated dot product of two vectors
-  `MakeHermitian.hpp`: force a matrix to be Hermitian
-  `MakeReal.hpp`: force a matrix to be real
-  `MakeSymmetric.hpp`: force a matrix to be symmetric
-  `MakeTrapezoidal.hpp`: force a matrix to be trapezoidal
-  `MakeTriangular.hpp`: force a matrix to be triangular
-  `Max.hpp`: return the maximum entry in a vector or matrix
-  `Nrm2.hpp`: compute the two-norm of a vector
-  `QuasiDiagonalScale.hpp`: apply a quasi-diagonal matrix
-  `QuasiDiagonalSolve.hpp`: apply the inverse of a quasi-diagonal matrix
-  `Scale.hpp`: scale a matrix or vector 
-  `ScaleTrapezoid.hpp`: scale a trapezoid of a matrix or vector
-  `SetDiagonal.hpp`: set the entire diagonal of a matrix to a particular value
-  `Swap.hpp`: swap the contents of two matrices, or swap two rows or columns of
   a particular matrix
-  `Symmetric2x2Inv.hpp`: form the inverse of a symmetric 2x2 matrix
-  `Symmetric2x2Scale.hpp`: apply a symmetric 2x2 matrix
-  `Symmetric2x2Solve.hpp`: apply the inverse of a symmetric 2x2 matrix
-  `Transpose.hpp`: form the transpose of a matrix
-  `UpdateDiagonal.hpp`: add a constant value to every diagonal entry of a 
   matrix
-  `Zero.hpp`: set every entry in a matrix to zero
### `include/elemental/blas-like/level2/`

The level-2 BLAS-like routines implemented in this folder are:

-  `ApplyColumnPivots.hpp`: apply a permutation matrix from the right
-  `ApplyRowPivots.hpp`: apply a permutation matrix from the left
-  `ApplySymmetricPivots.hpp`: apply a symmetric permutation to a symmetric
   matrix
-  `ComposePivots.hpp`: form the explicit image and preimage of a permutation
-  `Gemv/`: supporting routines for GEneral Matrix/Vector multiplication (GEMV)
-  `Gemv.hpp`: interface for GEMV
-  `Ger.hpp`: GEneral Rank-one update (GER)
-  `Geru.hpp`: unconjugated version of GER
-  `Hemv.hpp`: HErmitian Matrix/Vector multiplication (HEMV)
-  `Her2.hpp`: HErmitian Rank-2 update (HER2)
-  `Her.hpp`: HErmitian Rank-1 update (HER)
-  `Symv/`: supporting routines for SYmmetric Matrix/Vector multiplication 
   (SYMV)
-  `Symv.hpp`: interface for SYMV
-  `Syr2.hpp`: SYmmetric Rank-2 update (SYR2)
-  `Syr.hpp`: SYmmetric Rank-1 update (SYR)
-  `Trmv.hpp`: TRiangular Matrix/Vector multiplication (TRMV)
-  `Trr2.hpp`: TRiangular Rank-2 update (TRR2)
-  `Trr.hpp`: TRiangular Rank-1 update (TRR)
-  `Trsv/`: supporting routines for TRiangular Solve against Vector (TRSV)
-  `Trsv.hpp`: interface for TRSV

NOTE: The symmetric pivot applications are not yet optimized.
### `include/elemental/blas-like/level3/`

The level-3 BLAS-like routines implemented in this folder are:

-  `Gemm/`: supporting routines for GEneral Matrix/Matrix multiplication (GEMM)
-  `Gemm.hpp`: interface for GEMM
-  `Hemm.hpp`: HErmitian Matrix/Matrix multiplication (HEMM)
-  `Her2k.hpp`: HErmitian Rank-2k update (HER2K)
-  `Herk.hpp`: HErmitian Rank-k update (HERK)
-  `Symm/`: supporting routines for SYmmetric Matrix/Matrix multiplication 
   (SYMM)
-  `Symm.hpp`: interface for SYMM
-  `Syr2k/`: supporting routines for SYmmetric Rank-2k update (SYR2K)
-  `Syr2k.hpp`: interface for SYR2K
-  `Syrk/`: supporting routines for SYmmetric Rank-k update (SYRK)
-  `Syrk.hpp`: interface for SYRK
-  `Trdtrmm/`: supporting routiens for TRiangular/Diagonal/TRiangular Matrix 
   Multiplication (TRDTRMM)
-  `Trdtrmm.hpp`: interface for TRDTRMM
-  `Trmm/`: supporting routines for TRiangular Matrix Multiplication (TRMM)
-  `Trmm.hpp`: interface for TRMM
-  `Trsm/`: supporting routines for TRiangular Solve of Matrix (TRSM)
-  `Trsm.hpp`: interface for TRSM
-  `Trtrmm/`: supporting routines for TRiangular TRiangular Matrix 
   Multiplication (TRTRMM)
-  `Trtrmm.hpp`: interface for TRTRMM
-  `Trstrm/`: supporting routines for TRiangular Solve against TRiangular Matrix
   (TRSTRM)
-  `Trstrm.hpp`: interface for TRSTRM
-  `TwoSidedTrmm/`: supporting routines for overwriting a Hermitian matrix with
   its congruence relative to a given triangular matrix
-  `TwoSidedTrmm.hpp`: interface for `TwoSidedTrmm`
-  `TwoSidedTrsm/`: supporting routines for overwriting a Hermitian matrix with
   its congruence relative to the inverse of a given triangular matrix
-  `TwoSidedTrsm.hpp`: interface for `TwoSidedTrsm`
### `include/elemental/control/`

- `decl.hpp`: include all declarations of control-related routines
- `impl.hpp`: include all implementations of control-related routines

A few matrix sign function based solvers for control theory:

-  `Lyapunov.hpp`: Solves A X + X A' = C for X when A has its eigenvalues
   in the open right-half plane
-  `Ricatti.hpp`: Solves X K X - A' X - X A = L for X when K and L are 
   Hermitian.
-  `Sylvester.hpp`: Solves A X + X B = C for X when A and B both have all of 
   their eigenvalues in the open right-half plane

#### TODO

Implement algorithms from Benner, Quintana-Orti, and Quintana-Orti's 
"Solving Stable Sylvester Equations via Rational Iterative Schemes".
### `tests/`

This folder contains accuracy tests for some of Elemental's routines.
It is divided into the following subfolders:

-  `blas-like/`: BLAS-like functionality
-  `core/`: core data structures
-  `lapack-like/`: LAPACK-like functionality
### `tests/convex`

This folder stores the correctness tests for Elemental's functionality meant
to support convex optimization. It currently only contains the following test:

-  `TSSVT.cpp`: A test for Tall-Skinny Singular Value soft-Thresholding
### `tests/lapack-like`

This folder contains correctness tests of a few of Elemental's LAPACK-like 
routines. More details will hopefully follow soon.

-  `ApplyPackedReflectors.cpp`
-  `Cholesky.cpp`
-  `CholeskyQR.cpp`
-  `HermitianEig.cpp`
-  `HermitianGenDefiniteEig.cpp`
-  `HermitianTridiag.cpp`
-  `LDL.cpp`
-  `LQ.cpp`
-  `LU.cpp`
-  `QR.cpp`
-  `RQ.cpp`
-  `SequentialLU.cpp`
-  `TriangularInverse.cpp`
-  `TSQR.cpp`
### `tests/core`

This folder stores the correctness tests for Elemental's core functionality:

-  `AxpyInterface.cpp`: Tests the local-to-global and global-to-local Axpy 
   (y := alpha x plus y)  interface
-  `DifferentGrids.cpp`: Tests a redistribution between different process grids
-  `DistMatrix.cpp`: Tests various redistributions for the DistMatrix class
-  `Matrix.cpp`: Tests buffer attachment for the Matrix class
-  `Version.cpp`: Prints the version information of this Elemental build
### `tests/blas-like`

This folder contains correctness tests of a few of Elemental's BLAS-like 
routines. More details will hopefully follow soon.

-  `Gemm.cpp`
-  `Hemm.cpp`
-  `Her2k.cpp`
-  `Herk.cpp`
-  `Symm.cpp`
-  `Symv.cpp`
-  `Syr2k.cpp`
-  `Syrk.cpp`
-  `Trmm.cpp`
-  `Trsm.cpp`
-  `Trsv.cpp`
-  `TwoSidedTrmm.cpp`
-  `TwoSidedTrsm.cpp`
### `src/`

This folder contains Elemental's routines which are to be directly instantiated 
rather than keeping them datatype-agnostic. Several commonly-used large routines
are kept here in order to keep the build times reasonable. This folder contains 
the following subfolders:

-  `blas-like/`: BLAS-like routines
-  `core/`: Elemental's core data structures and environment
-  `io/`: input/output, such as Qt5 graphics
-  `lapack-like/`: LAPACK-like routines
### `src/io/`

This folder contains the source-level implementations of Elemental's 
input/output functionality. Please see `include/elemental/io` for the 
header-level implementations. In addition to this file, this folder contains:

-  `ComplexDisplayWindow.cpp`: a Qt5-based graphical display of a complex matrix
-  `DisplayWindow.cpp`: a Qt5-based graphical display of a real matrix
-  `SpyWindow.cpp`: a Qt5-based graphical display of the nonzeros of a matrix
### `src/lapack-like/`

This folder contains Elemental's source-level implementations of LAPACK-like
routines. Most such routines are implemented in `include/elemental/lapack-like`,
but the exceptions are `HermitianEig` and `HermitianTridiag`, which both involve
a substantial amount of code.

In addition to this file, this folder contains:

-  `HermitianEig.cpp`: Implementation of Hermitian eigensolvers
-  `HermitianTridiag/`: Underlying implementations of the reduction of a 
   Hermitian matrix to real symmetric tridiagonal form
-  `HermitianTridiag.cpp`: High-level interface for the reduction of a Hermitian
   matrix to real symmetric tridiagonal form
### `src/lapack-like/HermitianTridiag/`

This folder contains Elemental's code for reducing a Hermitian matrix to 
real symmetric tridiagonal form. The various pieces are organized as follows:

-  `L.hpp`: Lower-triangular storage
-  `LSquare.hpp`: Lower-triangular storage specialized to
   square process grids
-  `LPan.hpp`: Panel portion of a blocked algorithm for lower-triangular 
   storage
-  `LPanSquare.hpp`: Panel portion of a blocked algorithm for lower-triangular
   storage specialized to square process grids
-  `U.hpp`: Upper-triangular storage
-  `USquare.hpp`: Upper-triangular storage specialized to square process grids
-  `UPan.hpp`: Panel portion of a blocked algorithm for upper-triangular 
   storage
-  `UPanSquare.hpp`: Panel portion of a blocked algorithm for upper-triangular
   storage specialized to square process grids
### `src/core/`

This folder contains the directly-instantiated portions of Elemental's core
functionality, such as the `Matrix` and `DistMatrix` classes, and wrappers of
external libraries, such as BLAS and LAPACK. Please see 
`include/elemental/core` for the corresponding header-level implementations and
prototypes.

In addition to this file, the folder currently contains:

-  `dist_matrix/`: distributed matrix class (`DistMatrix`)
-  `global.cpp`: initialization/finalization, call-stack manipulation, etc.
-  `imports/`: wrappers for external software
-  `matrix.cpp`: sequential matrix class (`Matrix`)
-  `mpi_register.cpp`: custom MPI datatypes and reduction operations
### `src/core/dist_matrix/`

This folder contains the source code for the various partial specializations of
the `DistMatrix` class; please see `include/elemental/core/dist_matrix/` for the
corresponding header-level prototypes and a detailed README explaining the 
various data distributions.
### `src/core/dist_matrix/`

This folder contains the source code for the various partial specializations of
the `DistMatrix` class; please see `include/elemental/core/dist_matrix/` for the
corresponding header-level prototypes and a detailed README explaining the 
various data distributions.
### `src/core/imports/`

This folder contains the directly-instantiated portions of Elemental's wrappers
for external libraries, such as BLAS, LAPACK, and PMRRR. Please see 
`include/elemental/core/imports` for the corresponding header-level prototypes
and implementations.

In addition to this file, the folder contains:

-  `blas.cpp`: wrappers for Basic Linear Algebra Subprograms (BLAS)
-  `flame.cpp`: wrappers for FLAME's QR-based bidiagonal SVD
-  `lapack.cpp`: wrappers for Linear Algebra PACKage (LAPACK)
-  `mpi.cpp`: wrappers for the Message Passing Interface (MPI)
-  `pmrrr.cpp`: wrappers for Parallel Multiple Relatively Robust Representations
   (PMRRR)
### `src/blas-like/`

This folder contains Elemental's source files for BLAS-like routines. 
The vast majority are in the `include/elemental/blas-like` folder, but 
the following are directly instantiated:

-  `Trr2k/`: underlying implementations of rank-2k triangular updates
-  `Trr2k.cpp`: the high-level interface to rank-2k triangular updates
-  `Trrk/`: underlying implementations of rank-k triangular updates
-  `Trrk.cpp`: the high-level interface to rank-k triangular updates
### `src/blas-like/Trrk/`

This folder contains the underlying implementations of rank-k triangular 
updates. In particular, in addition to this file, it holds:

-  `Local.hpp`: all sequential implementations
-  `NN.hpp`: parallel normal/normal implementations
-  `NT.hpp`: parallel normal/transposed implementations
-  `TN.hpp`: parallel transposed/normal implementations
-  `TT.hpp`: parallel transposed/transposed implementations

#### Notes

There is currently one TODO item related to this folder:

1. Making the orientation options of LocalTrrk more consistent with Trrk
### `src/blas-like/Trr2k/`

This folder contains the underlying implementations of rank-2k triangular 
updates. In particular, in addition to this file, it holds:

-  `Local.hpp`: all sequential implementations
-  `NNNN.hpp`: parallel normal/normal/normal/normal implementations
-  `NNNT.hpp`: parallel normal/normal/normal/transpose implementations
-  `NNTN.hpp`: parallel normal/normal/transpose/normal implementations
-  `NNTT.hpp`: etc.
-  `NTNN.hpp`
-  `NTNT.hpp`
-  `NTTN.hpp`
-  `NTTT.hpp`
-  `TNNN.hpp`
-  `TNNT.hpp`
-  `TNTN.hpp`
-  `TNTT.hpp`
-  `TTNN.hpp`
-  `TTNT.hpp`
-  `TTTN.hpp`
-  `TTTT.hpp`

#### Notes

There are currently two TODO items related to this folder:

1. Making the orientation options of LocalTrr2k more consistent with Trr2k
2. Implementing sequential versions of Trr2k
### `examples/`

This folder contains (hopefully) concise examples of Elemental's functionality.
It is divided into the following subfolders:

-  `blas-like/`: BLAS-like functionality
-  `convex/`: convex optimization
-  `core/`: core data structures
-  `lapack-like/`: LAPACK-like functionality
-  `matrices/`: special matrices
### `examples/convex`

This folder contains a few examples of using Elemental for tasks related to
convex optimization:

-  `LogDetDivergence.cpp`: Compute a divergence based upon the log of the 
   determinant of an HPD matrix
-  `RPCA.cpp`: A simple Robust Principal Component Analysis example which allows
   for choosing between several different SVT algorithms
### `examples/matrices`

This folder contains several examples of Elemental's special matrices.
Details behind them will hopefully be added soon.

-  `Cauchy.cpp`:
-  `CauchLike.cpp`:
-  `Circulant.cpp`:
-  `Diagonal.cpp`:
-  `Egorov.cpp`:
-  `Fourier.cpp`:
-  `Hankel.cpp`:
-  `Helmholtz1D.cpp`:
-  `Helmholtz2D.cpp`:
-  `Helmholtz3D.cpp`:
-  `HermitianUniformSpectrum.cpp`:
-  `Hilbert.cpp`:
-  `Identity.cpp`:
-  `Kahan.cpp`:
-  `Legendre.cpp`:
-  `LehmerParterRis.cpp`:
-  `NormalUniformSpectrum.cpp`:
-  `Ones.cpp`:
-  `OneTwoOne.cpp`:
-  `PSFW.cpp`:
-  `RiemannRedhefferGCD.cpp`:
-  `Toeplitz.cpp`:
-  `Uniform.cpp`:
-  `Walsh.cpp`:
-  `Wilkinson.cpp`:
-  `Zeros.cpp`:
### `examples/lapack-like`

This folder contains several examples of Elemental's LAPACK-like functionality:

-  `BunchKaufman.cpp`: Accurate symmetric/Hermitian-indefinite factorization
-  `BusingerGolub.cpp`: Column-pivoted QR decomposition
-  `ComplexHermitianFunction.cpp`: Applies a complex function to the eigenvalues
   of a Hermitian matrix
-  `GaussianElimination.cpp`: Solves systems of equations via Gaussian elim.
-  `HermitianEig.cpp`: Computes the eigen{values/pairs} of a Hermitian matrix
-  `HermitianEigFromSequential.cpp`: Distributes a sequential Hermitian matrix, computes its EVD, and then gathers the result back to the original process
-  `HermitianPseudoinverse.cpp`: Forms the pseudoinverse of a Hermitian matrix
-  `HermitianQDWH.cpp`: A variant of the QDWH algorithm for the polar 
   decomposition which is specialized for Hermitian matrices
-  `HermitianSDC.cpp`: Spectral Divide and Conquer eigensolver for Hermitian 
   matrices
-  `HermitianSVD.cpp`: Singular Value Decomposition of a Hermitian matrix
-  `HPDInverse.cpp`: Inverts a Hermitian Positive-Definite matrix
-  `HPSDCholesky.cpp`: Computes the (non-unique) Cholesky decomposition of a 
   Hermitian Positive-SemiDefinite matrix via its eigenvalue decomposition
-  `HPSDSquareRoot.cpp`: Computes the square-root of a Hermitian 
   Positive-SemiDefinite matrix via its eigenvalue decomposition
-  `ID.cpp`: Computes an Interpolate Decomposition 
   (closely related to pivoted QR)
-  `KyFanAndSchatten.cpp`: Compute Ky Fan and Schatten norms
-  `LDL.cpp`: Unpivoted LDL^T/LDL^H factorization
-  `LDLInverse.cpp`: Invert a symmetric/Hermitian matrix via a pivoted 
   symmetric factorization (e.g., Bunch-Kaufman)
-  `LeastSquares.cpp`: Solve a least-squares problem via a QR decomposition
-  `Polar.cpp`: Compute a polar decomposition (unitary times HPD)
-  `Pseudoinverse.cpp`: Compute the pseudoinverse of an arbitrary matrix
-  `QDWH.cpp`: Compute the polar factor of an arbitrary matrix via the QDWH 
   algorithm
-  `QR.cpp`: Compute a QR decomposition
-  `RealHermitianFunction.cpp`: Apply a real function to the eigenvalues of a
   Hermitian matrix
-  `RealSchur.cpp`: Compute the Schur decomposition of a real matrix
-  `RealSymmetricFunction.cpp`: Apply a real function to the eigenvalues of a 
   (real) symmetric matrix
-  `Schur.cpp`: Compute the Schur decomposition of a matrix
-  `SequentialBunchKaufman.cpp`: Test the sequential algorithm for Bunch-Kaufman
-  `SequentialQR.cpp`: Test the sequential algorithm for QR decomposition
-  `SequentialSVD.cpp`: Test the sequential algorithm for SVD
-  `Sign.cpp`: Test the matrix sign function (maps eigenvalues to {-1,+1})
-  `SimpleSVD.cpp`: An extremely simple SVD driver
-  `Skeleton.cpp`: Compute a matrix skeleton
-  `SkewHermitianEig.cpp`: Compute the EVD of a skew-Hermitian matrix
-  `SVD.cpp`: Compute the SVD of an arbitrary matrix
### `examples/core`

This (currently sparse) folder contains:

-  `Constructors.cpp`: A few examples of DistMatrix constructors
### `examples/blas-like`

This folder contains a few examples of Elemental's BLAS-like functionality:

-  `Cannon.cpp`: An unoptimized implementation of Cannon's algorithm for 
   matrix-matrix multiplication
-  `Gemv.cpp`: Matrix-vector multiplication
-  `Gemm.cpp`: Matrix-matrix multiplication
# QCSchema Json

Integration of Madness into QCArchive requires to ability to translate data
generated by Madness into a standard quantum chemistry variable defined in 
`QCSchema`.  Within QCEngine, the madness `harvester.py` code is responible 
for reading in madness output files and harvesting relevant output data.  
This is done using the `regex` python api which pattern matches.  

`json.hpp` is included to allow applications to write json files with
relevant Quantum Chemistry data.  The ultimate goal is to be able
to generate data requested by QCArchive according to the `QCSchema` as
explained (QCSchema)['https://molssi-qc-schema.readthedocs.io/en/latest/auto_props.html'].

This ways we generate `json` output files that can easily be read within 
QCA.  


## Calculation Information 

`calc_info.json`

A list of fields that involve basic information of the requested computation.

- `calcinfo_nbasis`	The number of basis functions for the computation.	number
    - Not relevant in MADNESS
- `calcinfo_nmo`	The number of molecular orbitals for the computation.	number
    - MOLDFT nmo_alpha + nmo_beta ... derived
- `calcinfo_nalpha`	The number of alpha electrons in the computation.	number
    - nalpha
- `calcinfo_nbeta`	The number of beta electrons in the computation.	number
    - nbeta
- `calcinfo_natom`	The number of atoms in the computation.	number
    - molecule .size
- `return_energy`	The energy of the requested method, identical to return_value for energy computations.	number

## Self-Consistent Field Information

A list of fields added at the SCF level.  (HF and DFT)

`scf_info.json`

- scf_one_electron_energy	The one-electron (core Hamiltonian) energy contribution to the total SCF energy.	number
- scf_two_electron_energy	The two-electron energy contribution to the total SCF energy.	number
- nuclear_repulsion_energy	The nuclear repulsion energy contribution to the total SCF energy.	number
- scf_vv10_energy	The VV10 functional energy contribution to the total SCF energy.	number
- scf_xc_energy	The functional energy contribution to the total SCF energy.	number
- scf_dispersion_correction_energy	The dispersion correction appended to an underlying functional when a DFT-D method is requested.	number
- scf_dipole_moment	The X, Y, and Z dipole components.	array[number]
- scf_total_energy	The total electronic energy of the SCF stage of the calculation. This is represented as the sum of the … quantities.	number
- scf_iterations

## Wavefunction Schema

A list of valid quantum chemistry wavefunction properties tracked byu the schema.  
Matrices are in column-major order. AO basis functions are ordered according to the CCA standard as implemented in libint.




