![build](https://github.com/openmopac/mopac/actions/workflows/CI-CD.yaml/badge.svg)
[![codecov](https://codecov.io/gh/openmopac/mopac/branch/main/graph/badge.svg?token=qM2KeRvw06)](https://codecov.io/gh/openmopac/mopac)

The modern open-source version of the Molecular Orbital PACkage (MOPAC).

This version contains a CMake build script and compiles a library in addition to an executable for integration with other software.

This project is presently in a pre-release state, and we are planning a more official release in late 2021.
Fortran 90 source files for the MOPAC executable

Subdirectories correspond to different features and subsystems within MOPAC

Modules for storing global data are labeled "*_C.F90"

Interfaces for subroutines with arguments are labeled "*_I.F90"

Larger subroutines are in source files named after the subroutine
MOZYME linear-scaling electronic structure solver
molecular geometry calculations & transformations
INDO model & calculations
semiempirical models & their tabulated parameters
Explicit shared library interface for MOPAC. For use with Windows DLLs, the interface must be explicitly
defined and marked by 'dllexport' preprocessor statements. Since most of the input and output of MOPAC is
through the disk, its programming API is very spartan at the present time.
methods for extracting an effective electrostatic potential (ESP & PMEP)
output file writing
geometry optimization
properties calculations using the output of converged SCF calculations
COSMO solvation model
self-consistent field cycle
dense matrix operations

NOTE: remove diag_for_GPU & re-introduce diag, either as a historical method w/ a special keyword, or a new form of pseudo-diagonalization

NOTE: MKL-specific routines have been removed from the rest of the code, but some remain here.
The code had been forcing the maximum number of threads using:

    call mkl_set_num_threads(mkl_get_max_threads())

but this is not a good idea and has been removed. It prevents single-thread use of MOPAC, and MKL
has a sensible default (1 thread per core) if the number of threads isn't specified by MKL_NUM_THREADS.

NOTE: some files in this folder should be migrated back to other components once their matrix operations have been isolated
and implemented in a separate subroutine.

NOTE: there are subroutines inside larger files that need to be extracted & re-implemented w/ BLAS & LAPACK calls.
They will be listed here as they are identified and before they are fixed:

densf @ polar.F90
bmakuf @ polar.F90
bdenup @ polar.F90
bdenin @ polar.F90
aval @ polar.F90
trugud @ polar.F90
trudgu @ polar.F90
trugdu @ polar.F90
trsub @ polar.F90
transf @ polar.F90
tf @ polar.F90
coscl1 @ cosmo.F90
coscl2 @ cosmo.F90
determinant @ charst.F90 (recursive, exponential-scaling matrix determinant, YIKES)

NOTE: there are large temporary matrix workspaces throughout MOPAC that perhaps should be consolidated somehow.
They will be listed here for consideration:

polar @ polar.F90
symmetry detection & operations
Configuration Interaction
input file reading
simple utilities that are not associated with a specific feature of MOPAC
parameter fitting for new semiempirical models

These files are mostly self-contained to PARAM, but update is called by datin within MOPAC.

psort.F90 contains some hard-coded paths that need to be updated in some way.
atomic force & vibrational calculations & derivative routines that support them
Chemical detection & modification
matrix elements for semiempirical models
transition state searching
interatomic potential corrections

The D3 model is contained in the files:

copyc6.F90
dftd3.F90
dftd3_bits.F90
gdisp.F90
==================
MOPAC Contributors
==================

The development of MOPAC began in 1981 within the Dewar group at the University of Texas at Austin.
Prof. Michael Dewar and his group had been developing semiempirical models since the late 1960's,
and this activity had produced a substantial amount of internal software as well as two programs,
MINDO/3 and MNDO, released through the Quantum Chemistry Program Exchange (QCPE). The main developer
of MOPAC, James (Jimmy) J. P. Stewart, was visting the Dewar group on sabbatical leave from the
University of Strathclyde and took on the task of refactoring and consolidating this collection of
software into a more cohesive and user-friendly program with a unified system of input and output
files that was released as QCPE Program #455 in 1983.

The overwhelming majority of development and maintenance of MOPAC from 1981-2021 has been performed
by Jimmy Stewart. The longstanding policy of MOPAC has been to make source code freely available to
interested developers and to accept the donation of new features, provided either in a modified
version of MOPAC or a self-contained reference program to be integrated into MOPAC. One of the
stipulations of these donations was the prior academic publication of the new feature, and the
historic contributions listed below include a DOI link to the relevant publications in these cases.

Pre-MOPAC Contributors
======================

MINDO/3 Developers
------------------

The 1975 methodology paper [`DOI:10.1021/ja00839a001 <https://doi.org/10.1021/ja00839a001>`_]
was authored by Richard C. Bingham, Michael J. S. Dewar, and Donald H. Lo.
The MINDO/3 program [QCPE Program #279 (1975)] developers according to the QCPE listing were
M. J. S. Dewar, H. Metiu, P. J. Student, A. Brown, R. C. Bingham, D. H. Lo, C. A. Ramsden,
H. Kollmar, P. Werner, and P. K. Bischof.

MNDO Developers
---------------

The 1977 methodology paper [`DOI:10.1021/ja00457a004 <https://doi.org/10.1021/ja00457a004>`_]
was authored by Michael J. S. Dewar and Walter Thiel.
The MNDO program [QCPE Program #428 (1981)] developers according to the QCPE listing were
W. Thiel, P. Weiner, J. Stewart, and M. J. S. Dewar.

Historic MOPAC Contributors
===========================

Peter Pulay
   design & optimization of pseudodiagonalization
   [`DOI:10.1002/jcc.540030214 <https://doi.org/10.1002/jcc.540030214>`_]

Harry King & R. Nicholas Camp
   original implementation of the Camp-King SCF converger
   [`DOI:10.1063/1.441834 <https://doi.org/10.1063/1.441834>`_]

John McKelvey
   adaptation of the Camp-King converger for MOPAC & improved output formatting

James McIver, Jr. & Andrew Komornicki
   POWSQ geometry optimizer
   [`DOI:10.1021/ja00763a011 <https://doi.org/10.1021/ja00763a011>`_]

Roger Sargent, Dimitris Agrafiotis, & Henry Rzepa
   implementation of the Broyden-Fletcher-Goldfarb-Shanno (BFGS) optimizer.

Larry Davis & Larry Burggraf
   design of the dynamic reaction coordinate (DRC) & intrinsic reaction coordinate (IRC) features
   [`DOI:10.1002/jcc.540080808 <https://doi.org/10.1002/jcc.540080808>`_]

Frank Jensen
   efficiency improvements to the Eigenvector-Following (EF) method

Juan Carlos Paniagua
   improvements to the orbital localization procedure
   [`DOI:10.1002/qua.560260307 <https://doi.org/10.1002/qua.560260307>`_]

Jorge Medrano
   expanded bonding analysis
   [`DOI:10.1002/jcc.540060205 <https://doi.org/10.1002/jcc.540060205>`_]

Santiago Olivella
   semiempirical energy partitioning (ENPART subroutine)
   [`DOI:10.1002/jhet.5570180625 <https://doi.org/10.1002/jhet.5570180625>`_]

Tsuneo Hirano
   revision of energy partitioning & thermodynamic corrections

James Friedheim
   testing & bug hunting

Eamonn Healy
   testing & feature validation

James Ritchie
   bug fixes (SCF restarting)

Masamoto Togashi, Jerzy Rudzinski, Zdenek Slanina, & Eiji Osawa
   bug fixes (vibrational analysis)

Michael Frisch
   bug fixes (density matrix)

Patrick Redington
   bug fixes (heavy atom matrix elements)

Ernest Davidson
   improvements to the 2-electron matrix elements

Daniel Liotard
   partial analytical derivatives of the density matrix & 2-electron matrix elements

Yukio Yamaguchi
   partial support for analytical derivatives
   [`DOI:10.1016/0097-8485(78)80005-9 <https://doi.org/10.1016/0097-8485(78)80005-9>`_]

George Purvis III
   expanded STO-6G orbital implementation up to principal quantum number 6
   for use in analytical derivatives

Henry Kurtz
   implementation of polarizability and hyperpolarizability
   [`DOI:10.1002/jcc.540110110 <https://doi.org/10.1002/jcc.540110110>`_]

Prakashan Korambath
   frequency dependence of hyperpolarizability
   [`DOI:10.1021/bk-1996-0628.ch007 <https://doi.org/10.1021/bk-1996-0628.ch007>`_]

David Danovich
   implementation of point-group symmetry & Green's function corrections to ionization potentials
   [`DOI:10.1039/P29930000321 <https://doi.org/10.1039/P29930000321>`_]

Michael Coolidge
   use of symmetry to accelerate vibrational analysis
   [`DOI:10.1002/jcc.540120807 <https://doi.org/10.1002/jcc.540120807>`_]

Andreas Klamt
   implementation of the COSMO solvation model
   [`DOI:10.1039/P29930000799 <https://doi.org/10.1039/P29930000799>`_]

Anna Stewart
   copyediting of MOPAC documentation

Victor Danilov
   edited the MOPAC7 manual & identified bugs in the MECI feature

John Simmie
   conversion of the MOPAC7 manual to LaTeX

Walter Thiel & Alexander Voityuk
   reference implementation of semiempirical models with d orbitals
   [`DOI:10.1007/BF01134863 <https://doi.org/10.1007/BF01134863>`_]

Brent Besler & Kenneth Merz, Jr.
   implementation of atomic charge model for electrostatic potentials (ESP)
   [`DOI:10.1002/jcc.540110404 <https://doi.org/10.1002/jcc.540110404>`_]

Bingze Wang
   implementation of parametric electrostatic potentials (PMEP)
   [`DOI:10.1002/jcc.540150210 <https://doi.org/10.1002/jcc.540150210>`_]

Stephan Grimme
   reference implementation of the D3 dispersion model
   [`DOI:10.1063/1.3382344 <https://doi.org/10.1063/1.3382344>`_]

Jan Rezac
   expanded implementation of classical energy corrections (hydrogen bonding, halogen bonding, dispersion)
   [`DOI:10.1021/ct200751e <https://doi.org/10.1021/ct200751e>`_]

Gerd Rocha
   expanded BLAS/LAPACK support, Intel MKL for multi-threading, & cuBLAS/MAGMA for GPU acceleration
   [`DOI:10.1021/ct3004645 <https://doi.org/10.1021/ct3004645>`_]

Rebecca Gieseking
   implementation of the INDO/S spectroscopy model
   [`DOI:10.1002/jcc.26455 <https://doi.org/10.1002/jcc.26455>`_]

Open-Source MOPAC Contributors
==============================

Jonathan Moussa
   reorganization & clean-up of the codebase, portability testing & debugging, minor performance tuning
=====================
Contributing to MOPAC
=====================

MOPAC is an open source project and, as such, contributions are welcome and appreciated.
The standard, minimal legal stipulations apply: all contributions must be provided under
the LGPL 3.0 license and all contributors must agree to the DCO by using their real name
and a valid email address in their git commit statements.

How to contribute
=================

The preferred method is to fork the project on `GitHub <https://github.com/openmopac/MOPAC/>`_,
make your changes starting from the master branch, and then create a pull request.
MOPAC does not presently have strong stylistic or structural requirements, but
contributions are encouraged to conform to the existing codebase as much as possible.

While not a strict requirement, contributions should pass all available tests,
either as they are or with their reference outputs adjusted with an explanation
for why they needed to be adjusted. Contributions may be delayed or rejected if they
fail a test in a meaningful way that is not straightforward to correct.

What to contribute
==================

Many types of contributions, including non-software, are welcome:

- bug reports (preferably with a representative example)
- new tests
- new features
- performance improvements
- code quality improvements
- feature requests

Attribution
===========

Contributors are welcome to add their names to the `AUTHORS.rst` file as part
of their pull request if they are not already listed as a contributor.
Affiliations at the time of the contribution can be included if desired.
For major contributions, a brief summary of the contribution can also be included.

Developer Certificate of Origin
===============================

Contributors to MOPAC are assumed to agree with the Linux Foundation's 
`Developer Certificate of Origin <https://developercertificate.org/>`::

    Developer Certificate of Origin
    Version 1.1
    
    Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
    1 Letterman Drive
    Suite D4700
    San Francisco, CA, 94129
    
    Everyone is permitted to copy and distribute verbatim copies of this
    license document, but changing it is not allowed.
    
    
    Developer's Certificate of Origin 1.1
    
    By making a contribution to this project, I certify that:
    
    (a) The contribution was created in whole or in part by me and I
        have the right to submit it under the open source license
        indicated in the file; or
    
    (b) The contribution is based upon previous work that, to the best
        of my knowledge, is covered under an appropriate open source
        license and I have the right under that license to submit that
        work with modifications, whether created in whole or in part
        by me, under the same open source license (unless I am
        permitted to submit under a different license), as indicated
        in the file; or
    
    (c) The contribution was provided directly to me by some other
        person who certified (a), (b) or (c) and I have not modified
        it.
    
    (d) I understand and agree that this project and the contribution
        are public and that a record of the contribution (including all
        personal information I submit with it, including my sign-off) is
        maintained indefinitely and may be redistributed consistent with
        this project or the open source license(s) involved.
