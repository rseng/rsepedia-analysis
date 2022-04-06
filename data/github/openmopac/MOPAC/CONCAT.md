![build](https://github.com/openmopac/mopac/actions/workflows/CI-CD.yaml/badge.svg)
[![codecov](https://codecov.io/gh/openmopac/mopac/branch/main/graph/badge.svg?token=qM2KeRvw06)](https://codecov.io/gh/openmopac/mopac)

The modern open-source version of the Molecular Orbital PACkage (MOPAC).

This version contains a CMake build script and compiles a library in addition to an executable for integration with other software.

This project is presently in a pre-release state, and we are planning a more official release in late 2021.
# MOPAC source code

This directory contains all Fortran source code associated with MOPAC and accompanying programs
(PARAM, MAKPOL, & BZ). Historically, MOPAC was developed with all source code contained in a single,
unorganized source directory. The organization of source code into subdirectories here is an attempt
to identify and reverse-engineer the substructure of distinct features and contributions to MOPAC,
to make it easier for future open-source developers to understand the codebase.

Most source files are named after the first subroutine in the file, and all data modules are contained
in separate source files with the suffix `*_C.F90`.
# Self-consistent field cycle

This directory contains MOPAC's self-consistent field (SCF) cycle and supporting functionality. The
core logic of the SCF cycle is contained in `iter.F90`.
# MOPAC API

This directory contains the application programming interface (API) to MOPAC, which presently has
very limited functionality. MOPAC's primary interface is disk-based rather than API-based, with
input passed in through input files and output read from output files (the AUX file being
machine-readable). The only functionality available right now is access to version information
and an ability to enable/disable expanded disk-based output for GUI support.

To allow MOPAC to function with Windows DLL's, the API subroutines must be marked by `dllexport` preprocessor statements that are used by the Intel Fortran compiler to build a DLL symbol table.
# Matrix operations

This directory contains all of the files in MOPAC that have been identified as mostly containing
dense linear algebra operations. They have been centrally collected in preparation for their
consolidation into a more unified subsystem for dense linear algebra within MOPAC. In particular,
the prior support for GPUs in MOPAC, which is not presently functioning, needs to be re-introduced
in a more uniform and maintainable way. Some of the files in this directory may be relocated when
their matrix operations have been appropriately encapsulated.

`ddiag.F90` is from the EISPACK library of eigensolvers. `linpack.F90` is from LINPACK library of
linear solvers. `minv.F90` is from the IBM-SSP library of mathematical and statistical subroutines.
`esp_utilities.F90` contains subroutines from LINPACK and BLAS. `interp.F90` was derived and adapted
from software written by R. Nicholas Camp and Harry F. King.

There are also cases where dense matrix storage and operations are directly embedded in other files
that have not been relocated to this directory. The present list of such subroutines and files is:

- `densf @ polar.F90`
- `bmakuf @ polar.F90`
- `bdenup @ polar.F90`
- `bdenin @ polar.F90`
- `aval @ polar.F90`
- `trugud @ polar.F90`
- `trudgu @ polar.F90`
- `trugdu @ polar.F90`
- `trsub @ polar.F90`
- `transf @ polar.F90`
- `tf @ polar.F90`
- `coscl1 @ cosmo.F90`
- `coscl2 @ cosmo.F90`
- `determinant @ charst.F90`
- `polar @ polar.F90`
# Geometry optimization

This directory contains various optimizers that are primarily used for geometry optimization. These
are all established optimization algorithms, and many of the implementations in MOPAC were derived
from standard, publicly-available implementations.

`powsq.F90` was originally integrated into the pre-MOPAC MINDO/3 codebase within the Dewar group by
one of its original authors, Andrew Komornicki. `flepo.F90` and `nllsq.F90` seem to have been from
the original MNDO program, but this provenance is unclear. `lbfgs.F90` was derived and adapted from
Algorithm 778 of the ACM Transactions on Mathematical Software [https://doi.org/10.1145/279232.279236].
# Chemical detection & modification

This directory contains subroutines for determining and manipulating Lewis dot structures of
organic molecules (with some support for organometallic molecules). Its functionality includes
an ability to add hydrogen atoms, which are often missing from experimental structures determined
by x-ray crystallography. It also identifies standard functional groups and can label them according
to the PDB file format.

The optimizer of water molecule orientation in `orient_water.F90` is derived and adapted from the Constrained Optimization BY Linear Approximation (COBYLA) algorithm and software originally written
by M. J. D. Powell [https://www.zhangzk.net/software.html].
# MAKPOL

This directory contains the MAKPOL program for generating supercells for use by MOPAC.
It is a useful tool for converging finite-size effects for  periodic calculations in MOPAC,
since there is no support in MOPAC for Brillouin zone sampling (all periodic calculations
are performed at the Gamma point). There is some overlap between the source code of MOPAC
and MAKPOL, but they are being kept separate for now (merge efforts are welcome).
# Property calculations

This directory contains various property calculations based on the output of MNDO-form calculations.

`static_polarizability.F90` was written by Henry Kurtz, and `polar.F90` was written by Henry Kurtz and
Prakashan Korambath. `CPE_Energy.F90` was incorporated into MOPAC as part of its interface with the
Sybyl GUI and was presumably written by the developers of Sybyl.
# Analytical & numerical derivatives for forces & vibrations

This directory contains most of the functionality for calculating 1st and 2nd derivatives
of the Fock matrix and interatomic potential corrections with respect to atomic coordinates
for use in force and vibrational calculations.

A significant portion of the analytical derivatives were derived and adapted from
Daniel Liotard's contributions to AMPAC version 2.1 (QCPE No. 506).
# MNDO Hamiltonian matrix elements

This directory contains an implementation of the Hamiltonian matrix elements for the
Modified Neglect of Differential Overlap (MNDO) model. The original MNDO model was developed
specifically for s and p orbitals, and was expanded to d orbitals in the 1990's.

The s/p-orbital components were derived and adapted from the original MNDO program, and the
d-orbital components were derived and adapted from the MNDO96 program written by Walter Thiel.

Several subroutines in this directory are attributed to E. R. Davidson, but their origin is unclear.
# Transition states

This directory contains several different methods for calculating transition states and reaction
coordinates.

`ef.F90` was derived and adapted from software written by Frank Jensen.
# Classical molecular-mechanics corrections

This directory contains all of the classical interatomic potential terms that are used to
correct the quantum mechanical calculations performed by MOPAC. The oldest corrections were
introduced to fix very specific problems in certain covalent bonds (e.g. triple bonds), and
the more recent corrections have focused on weak intermolecular interactions, mainly
dispersion and hydrogen bonding.

The D3 dispersion model implementation was derived and adapted from the DFTD3 library
[https://github.com/dftbplus/dftd3-lib] with permission from Stephan Grimme. The files
derived from this source are `copyc6.F90`, `dftd3.F90`, `dftd3_bits.F90`, and `gdisp.F90`.

The H4 model implementation was derived and adapted from the h_bonds4 library written by
Jan Rezac, available either as a standalone C implementation [https://www.rezacovi.cz]
or in Ruby as part of the Cuby4 framework [http://cuby4.molecular.cz]. The file derived from
this source is `H_bonds4.F90`.
# MOZYME

This directory contains the MOZYME feature of MOPAC, which performs reduced-scaling electronic
structure calculations based directly on localized molecular orbitals (LMOs) instead of canonical
orbitals, which allows for the use of sparse linear algebra rather than dense linear algebra. Many
of the core computational steps are then reduced from the cubic-scaling cost of dense linear algebra
to a linear-scaling cost in the number of atoms, but there are still steps and prefactors that cause
overall super-linear scaling, notably electrostatics are still performed with brute-force pair
summation rather than any fast hierarchical method and the SCF cycle and geometric relaxation cycle
inevitably retains some scaling with system size, depending on the degree of electronic and atomic
polarizability in a system.

This feature of MOPAC was formerly covered by a US patent owned by Fujitsu [5,604,686],
which has now expired.
# COSMO solvation model

This directory contains MOPAC's implementation of the COSMO solvation model,
which was developed and written by Andreas Klamt.
# Point-group symmetry

This directory contains functionality for detecting and manipulating point-group symmetry
of molecules, which was written by David Danovich.
# INDO spectroscopy model

This directory contains an implementation of the INDO/S model and excited state calculations
using the INDO/S model Hamiltonian with electron correlation treated at various levels of theory.

The software in this directory was written by Rebecca Gieseking. It was derived and adapted from
the CNDO/INDO program written by Jeff Reimers [DOI:10.4231/D3R49G96G], which was in turn derived
and adapted from the CNDO/S program [QCPE 174] written by J. Del Bene, H. H. Jaffe, R. L. Ellis,
and G. Kuehnlenz.
# Electrostatic potentials

This directory contains two methods (ESP & PMEP) for mapping the density matrix produced by AM1
calculations to an electrostatic potential that attempts to match *ab initio* Hartree-Fock 
calculations. This methods are considered to be outdated and may be removed from future versions
of MOPAC.

`esp.F90` was written by Kenneth Merz, Jr. and `pmep.F90` was written by Bingze Wang.
# Output file processing

This directory contains functionality for the processing of output files. Much of the standard
`.out` file is printed by `writmo.F90`, and the machine-readable `.aux` file is printed by
`to_screen.F90`.

The energy decomposition in `enpart.F90` was a contribution to the original MNDO program written by
Santiago Olivella.
# BZ

This directory contains the Windows BZ program for visualizing band structures and Fermi
surfaces from periodic calculations performed by MOPAC. It generates a Brillouin zone by
downfolding the supercells created by MAKPOL to recover the original unit cell. This program
is presently available for Windows only as it uses the QuickWin library for visualization.
There is some overlap between the source code of MOPAC and BZ, but they are being kept
separate for now (merge efforts are welcome).

`blas.F90` is from the BLAS library of basic linear algebra subroutines. `minv.F90` is from
the IBM-SSP library of mathematical and statistical subroutines. `rsp.F90` is derived and
adapted from the EISPACK library of eigenvalue problem solvers.
# Configuration Interaction

This directory contains the inter-determinant matrix elements and core computational routines for
performing truncated configuration interaction (CI) calculations. While not physically meaningful
according to the strictest interpretation of a thermochemitry model (that is only intended to model
heats of formation with Hartree-Fock calculations), these types of calculations are useful for
understanding the sufficiency of these models and the effects of strong correlation in certain
situations such as bond breaking and multi-radicals.

As annotated, several subroutines were derived and adapted from Daniel Liotard's contributions to
AMPAC version 2.1 (QCPE No. 506).
# Miscellaneous utilities

This directory contains simple, self-contained utilities not directly associated with any of the
main features of MOPAC.
# PARAM

This directory contains subroutines and modules specific to the PARAM program for optimizing
parameters of MNDO-form models. The rest of its source code is shared with MOPAC.

*NOTE: psort.F90 contains some hard-coded paths that need to be updated in some way.*
# Semiempirical model parameters

This directory contains hard-coded tables for the core semiempirical parameters that define the
MNDO-form model Hamiltonian for the various model fits supported by MOPAC. It also contains
functionality for rerouting data to switch between models.
# Input file processing

This directory contains functionality for the processing and parsing of input files. Note that
the input file is copied to a scratch file by `getdat.F90` before any other processing is performed.
The primary parsing of keywords occurs in `wrtkey.F90`, and any new keywords should be introduced
there before adding any other keyword-based control logic associated with that keyword.
# Molecular geometry calculations & transformations

This directory contains some basic subroutines for transforming molecular geometries
between coordinate systems and calculating some derived geometric quantities.

As noted in the source, `dihed.F90` was derived and adapted from software written by
Walter Thiel, and `gmetry.F90` was derived and adapted from software written by Michael Dewar.
The provenance of this software is not completely clear, but it is likely to have come from
either the MINDO/3 or MNDO programs.
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
by Jimmy Stewart. The longstanding policy of MOPAC has been to make its source code freely available
to interested developers and to accept the donation of new features, provided either in a modified
version of MOPAC or a self-contained reference program to be integrated into MOPAC. One of the
stipulations of these donations was the prior academic publication of the new feature, and the
historic contributions listed below include a DOI link to the relevant publications in these cases.

As MOPAC has an old codebase, attribution will inevitably be imperfect and have error bars. This
contributor list is an attempt at being maximally inclusive, and requests for corrections are welcome.
From a minimally inclusive perspective, Jimmy Stewart was the pre-open-source copyright owner of
MOPAC (donations to the code were made through a copyright transfer agreement), wrote the majority
of the code, and has borne all responsibility for maintenance, support, and promoting the interests
of MOPAC. Between this minimum and maximum, an earnest attempt at a formal author list is contained
in `CITATION.cff`, which signifies major contributions that persistent in the present MOPAC codebase
that are not part of some other original work of which MOPAC is a derivative work.

Pre-MOPAC Contributors
======================

MINDO/3 Authors
---------------

The 1975 methodology paper [`DOI:10.1021/ja00839a001 <https://doi.org/10.1021/ja00839a001>`_]
was authored by Richard C. Bingham, Michael J. S. Dewar, and Donald H. Lo.
The MINDO/3 program [QCPE Program #279 (1975)] authors according to the QCPE listing were
M. J. S. Dewar, H. Metiu, P. J. Student, A. Brown, R. C. Bingham, D. H. Lo, C. A. Ramsden,
H. Kollmar, P. Werner, and P. K. Bischof.

MNDO Authors
------------

The 1977 methodology paper [`DOI:10.1021/ja00457a004 <https://doi.org/10.1021/ja00457a004>`_]
was authored by Michael J. S. Dewar and Walter Thiel.
The MNDO program [QCPE Program #428 (1981)] authors according to the QCPE listing were
W. Thiel, P. Weiner, J. Stewart, and M. J. S. Dewar.

*NOTE: As with much of the software previously distributed through the QCPE, the original MINDO/3
and MNDO programs are no longer distributed and are likely to be permanently lost. If anyone has
a copy of either of these programs, then please consider sharing them on GitHub for the purpose
of historical preservation and better understanding of the origins of MOPAC.*

Historic MOPAC Contributors
===========================

Andrew Komornicki
   adapted the POWSQ geometry optimizer
   [`DOI:10.1021/ja00763a011 <https://doi.org/10.1021/ja00763a011>`_]
   to the MINDO/3 program after its QPCE release
   [`DOI:10.1021/ja00444a012 <https://doi.org/10.1021/ja00444a012>`_]
   but before the development of MOPAC itself

Santiago Olivella
   semiempirical energy partitioning (ENPART subroutine)
   [`DOI:10.1002/jhet.5570180625 <https://doi.org/10.1002/jhet.5570180625>`_]

Peter Pulay
   design & optimization of pseudodiagonalization
   [`DOI:10.1002/jcc.540030214 <https://doi.org/10.1002/jcc.540030214>`_]

Harry King & R. Nicholas Camp
   original implementation of the Camp-King SCF converger
   [`DOI:10.1063/1.441834 <https://doi.org/10.1063/1.441834>`_]

John McKelvey
   adaptation of the Camp-King converger for MOPAC & improved output formatting

Roger Sargent, Dimitris Agrafiotis, & Henry Rzepa
   implementation of the Broyden-Fletcher-Goldfarb-Shanno (BFGS) optimizer.

Larry Davis & Larry Burggraf
   design of the dynamic reaction coordinate (DRC) & intrinsic reaction coordinate (IRC) features
   [`DOI:10.1002/jcc.540080808 <https://doi.org/10.1002/jcc.540080808>`_]

Frank Jensen
   efficiency improvements to the Eigenvector-Following (EF) method
   [`DOI:10.1063/1.469144 <https://doi.org/10.1063/1.469144>`_]

Juan Carlos Paniagua
   improvements to the orbital localization procedure
   [`DOI:10.1002/qua.560260307 <https://doi.org/10.1002/qua.560260307>`_]

Jorge Medrano
   expanded bonding analysis
   [`DOI:10.1002/jcc.540060205 <https://doi.org/10.1002/jcc.540060205>`_]

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
   [`DOI:10.1016/0166-1280(90)85012-C <https://doi.org/10.1016/0166-1280(90)85012-C>`]

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

Kenneth Merz, Jr.
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

Jonathan Moussa
   open-source transition: reorganization & clean-up of the codebase, portability testing & debugging,
   minor performance tuning, transition to CMake-based build system, automation of continuous integration & deployment

Open-Source MOPAC Contributors
==============================

One of the major benefits of modern open-source software development is that contributions are passively recorded by
git version control. As such, refer to the commit records for a complete list of contributions. Major new feature
contributions will continue to be added to this list over time.
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

Historically, major feature contributions to MOPAC were required to be preceded by a journal
publication to maintain clear attribution. While this is no longer necessary for open-source
software developed under automated version control, it is still encouraged.

Major feature contributions should include new tests that cover standard use cases of the
new feature, and corresponding additions to the manual [https://openmopac.github.io],
particularly if new keywords are introduced.

What to contribute
==================

Many types of contributions are welcome:

- bug reports (preferably with a representative example)
- new tests
- new features
- performance improvements
- code quality improvements
- feature requests

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
