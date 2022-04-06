Nightly runs on full testset: https://testboard.org/cdash/index.php?project=Dalton

# Copyright and License
Dalton: A Molecular Electronic Structure Program
Copyright (C) by the authors of Dalton.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License version 2.1 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

If a copy of the GNU LGPL v2.1 was not distributed with this
code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.


# Dalton Links

- [Home page](https://daltonprogram.org/)
- [User Support](https://gitlab.com/dalton/user-support)
- [Article](https://doi.org/10.1002/wcms.1172)


# Quick Start

Note that it is currently not practical to download the source using the
download button on GitLab, because it will not include the submodules that are
required to build Dalton. Instead you should clone the repository as described
below.

Clone the repository:
```
$ git clone --recursive https://gitlab.com/dalton/dalton.git
```

This will fetch the entire repository in a directory called *dalton*. By default
it checks out the master branch which is the main development branch. To
checkout a specific release version, run the following commands from inside the
*dalton* directory:
```
$ git checkout Dalton2020.1
$ git submodule update
```
where you replace *Dalton2020.1* by the release version that you are
interested in. The list of past releases available in this repository can be
found here: https://gitlab.com/dalton/dalton/-/releases.

You can also clone the release version directly as:
```
$ git clone --recursive -b Dalton2020.1 https://gitlab.com/dalton/dalton.git
```

In case you did not include the `--recursive` argument when you cloned the
repository, it is necessary to run the following two commands:
```
$ git submodule update --init --recursive
```

To build the code, perform the following steps:
```
$ ./setup
$ cd build
$ make [-j4]
```

There are several setup options available, e.g., for setting up an MPI build.
To see the available options run:
```
$ ./setup --help
```

Once the build is complete, you can run the test set as:
```
$ ctest [-j4] -L dalton
```

To switch branch (or release tag), run the following two commands from the *dalton* directory:
```
$ git checkout feature-branch
$ git submodule update
```
This can also be achieved in one step when you clone the repository:
```
$ git clone --recursive -b feature-branch https://gitlab.com/dalton/dalton.git
```
# DALTON Change Log -- All notable changes to the DALTON program will be documented in this file.

## [2022.0-dev] (unreleased)

### New features added
- EKT, i.e. ionization potentials using the extended Koopmans' theorem for MCSCF and MC-srDFT. Is always activated for these types of wave functions. (Martin J. R. Jensen and H. J. Aa. Jensen)
- SOPPA (Luna Zamok)
  - Block Lanczos RPA eigenvalue solver for computing mean excitation energy and the dipole oscillator strength sums in Lanczos basis at RPA level. The solver is a part of AO-SOPPA and its driver is called from ABACUS driver.
  - Added PE-AO-SOPPA/RPA (Peter Reinholdt)
- Atomic Integrals (Juanjo Aucar)
  - Divergence of the zz Electric Field Gradient component
  - Laplacian of xx,yy and zz Electric Field Gradient components
- LRESC (Juanjo Aucar)
  - Implementation for the Electric Field Gradient at first order in 1/c2
- The MP3 model has been added to the CC module for the calculation of ground-state energies. (Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal, Nicolai Machholdt Høyer, and Kurt V. Mikkelsen)
- ENUE - print effective number of unpaired electrons ENUE = sum(i) (2\*n\_i - n\_i\*\*2) in final wave function output (H. J. Aa. Jensen)


## [2020.1] (2022-01-20)

### New features added
- New features available through the Polarizable Embedding library (PElib)
  - Fast multipole method (FMM) for linear-scaling evaluation of static/induced multipole fields (P. Reinholdt and J. M. H. Olsen)
    - M. Scheurer, P. Reinholdt, J. M. H. Olsen, A. Dreuw, J. Kongsted, J. Chem. Theory Comput. 17, 3445-3454 (2021)
  - Additional solvents available for the FixSol continuum solvation model (see 26.1.5 in the manual for a list of available solvents)
  - Fast approximate environment coupling for PE (P. Reinholdt)
    - P. Reinholdt, J. Kongsted, and F. Lipparini, J. Chem. Theory Comput. 18, 344-356 (2022)
- Allow combination of DFT with CI/MCSCF, i.e. use DFT orbitals in CI or as initial guess in MCSCF (H. J. Aa. Jensen)

### Fixed
- Performance improvement: disable DFT\_SPINDNS if true when singlet (H. J. Aa. Jensen)
- Fixed bug for molden.inp: empty orbitals got arbitrary occupation numbers. (H. J. Aa. Jensen)
- Fixed bug for \*CUBE (crashed always) (H. J. Aa. Jensen)
- Fixed error in 2-el. integrals for FCKTRA with MPI and MC-srDFT (H. J. Aa. Jensen)
- Fixed error leading to test energy\_selected\_localize failed (H. J. Aa. Jensen)
- Fixed PElib compile error when using explicit math libs
- Fixed error for SDKE, relativistic kinetic energy corrections to spin-dipole. Has never worked since it was introduced in 2004. (H. J. Aa. Jensen)
- Quit nicely if orbital relaxation is on for excited-state first-order property CC calculations


## [2020.0] (2020-10-20)

### Major new features added
- MC-srDFT code (H. J. Aa. Jensen, E. Fromager, S. Knecht, E. R. Kjellgren and others)
  - HF-srDFT, MP2-srDFT, CAS-srDFT, RAS-srDFT, NEVPT2-srDFT ground state energies and wave functions
  - State specific CAS-srDFT and RAS-srDFT for excited states
  - Singlet and triplet excitation energies and transition moments for HF-srDFT and MC-srDFT (linear response module)
  - Singlet and triplet excitation energies and transition moments from SOPPA-srDFT, using MP2-srDFT
  - srLDA and srPBEGWS short-range functionals
- Added ".TDA TRIPLET" keyword for invoking Tamm-Dancoff approximation for triplet response properties under \*\*PROPERTIES.
  Useful for avoiding (near-)triplet-instability problems, for example in DFT or MC-srDFT calculations of spin-spin coupling constants. (H. J. Aa. Jensen)
- Added ".TDA SINGLET" keyword for invoking Tamm-Dancoff approximation for singlet response properties under \*\*PROPERTIES. (H. J. Aa. Jensen)
- Added the ability in QFITLIB to fit up to and including quadrupoles (C. Steinmann)
- Added the core-valence separation (CVS) approximation for CC calculations of core-excited states (S. Coriani et al.)
- Added the possibility to calculate triplet-triplet excited state moments using the EOM-CC approximation (R. Faber)
- New features available through the Polarizable Embedding library (PElib)
  - Polarizable density embedding (PDE) model (use -DENABLE\_PDE=ON during setup to enable it [requires HDF5])
    - J. M. H. Olsen, C. Steinmann, K. Ruud, and J. Kongsted, J. Phys. Chem. A 119, 5344 (2015)
    - P. Reinholdt, J. Kongsted, and J. M. H. Olsen, J. Phys. Chem. Lett. 8, 5949 (2017)
  - PDE-CC2, PDE-CCSD, and PDE-CCSDR(3) including linear and quadratic response (also enables PE-CC through PElib)
    - D. Hrsak, J. M. H. Olsen, and J. Kongsted, J. Chem. Theory Comput. 14, 1351 (2018)
  - FixSol continuum solvation with FIXPVA2 cavity tesselation
    - M. S. Nørby, C. Steinmann, J. M. H. Olsen, H. Li, and J. Kongsted, J. Chem. Theory Comput. 12, 5050 (2016)
    - N. M. Thellamurege and H. Li, J. Chem. Phys. 137, 246101 (2012)
  - Enabled cubic response for PE-HF/DFT and PDE-HF/DFT
    - J. M. H. Olsen and J. Kongsted, Adv. Quantum Chem. 61, 107 (2011)
  - Effective external field (EEF) can now be enabled for all dipole properties
    - N. H. List, H. J. Aa. Jensen, and J. Kongsted. Phys. Chem. Chem. Phys. 18, 10070 (2016)
  - Added support for AMOEBA potential (P. Reinholdt & J. M. H. Olsen)
  - Added pseudopotentials for avoiding electron spill-out
    - A. M. Khah, P. Reinholdt, J. M. H. Olsen, J. Kongsted and C. Hattig, J. Chem. Theory. Comput. 16, 1373-1381 (2020)
- Added easy ghost atom input in the molecule input (P. Reinholdt & H. J. Aa. Jensen)

### Other new features added
- Added the possiblitly to create the Dalton pdf manual with "make pdfmanual" in the build directory. (H. J. Aa. Jensen)
- Added the possiblitly to create the Dalton html manual with "make htmlmanual" in the build directory. (H. J. Aa. Jensen)
- Added the ability to run SOPPA linear response calculations via the AOSOPPA module (R. Faber et al.)
- Added possibility to optimize MCSCF singlet wave functions with CSFs when used for triplet properties,
  both in \*\*RESPONS and \*\*PROPERTIES. Previously .DETERMINANTS in wave function optimization was required. (H. J. Aa. Jensen)
- Extended "super matrix integrals" to work with >255 basis functions. (H. J. Aa. Jensen)
- Added information about .MS2 input option to manual, quit if invalid value specified. (H. J. Aa. Jensen)

### Fixed
- .FCKTRA error: would sometimes skip 2-el. integral transformation when MOs had changed. (H. J. Aa. Jensen)
- Compilation with 64-bit integer and linking with a 32-bit integer MPI. (H. J. Aa. Jensen)
- Errors when running DFT with 64-bit integers and MPI. (P. Reinholdt and H. J. Aa. Jensen)
- Errors for Fermi-contact (FC) labels on APROPER and therefore FC properties in \*\*RESPONS when more than 99 atoms (H. J. Aa. Jensen)
- Error for \*ESR spin-dipole properties when more than 33 atoms (H. J. Aa. Jensen)
- Singlet totally-symmetric excitation energies for MCSCF in \*\*RESPONS with super-symmetry activated (.SUPSYM keyword). (H. J. Aa. Jensen)
- Make sure we include all (near-)degenerate diagonal elements for linear response excitation energies
  via \*\*PROPERTIES .EXCITA or via \*\*RESPONS \*LINEAR .SINGLE (increase .NROOTS if needed).
  Otherwise the calculation will probably exhibit spin and/or space symmetry contamination proportional
  to the convergence threshold. (H. J. Aa. Jensen)
- Never use plus combinations of determinants as start guess for singlet linear response excitation energies
  when reference wave function is not singlet (we do not want singlet states then). (H. J. Aa. Jensen)
- Dalton script: fix for using input files located in subfolders
- Fixed error from March 2015 which meant that double-hybrid DFT was not working correctly (MP2 part was ignored).
- Fixed error for MC-TDA excitation energies for RASSCF (CASSCF was OK).
- Fixed implementation of RPBEx functional with functional derivatives generated with the Python SymPy library

### Additions and fixes in enclosed basis set files
- Added pcH-n and aug-pcH-n basis sets levels 1-4 to Dalton basis set library.
- Error in diffuse d-orbital exponents for Aluminum and Silicon (factor 10 too big) in aug-cc-pV(D+d)Z basis sets (H. J. Aa. Jensen)
- Error in diffuse f-orbital exponents for Aluminum and Silicon (factor 10 too big) in aug-cc-pV(Q+d)Z basis sets (H. J. Aa. Jensen)
- Error in diffuse f-orbital exponent for Aluminum (factor 10 too big) in aug-cc-pV(T+d)Z basis sets (H. J. Aa. Jensen)

### Changed
- Allow basis set(s) after BASIS in line 1 of .mol file (instead of on second line). (H. J. Aa. Jensen)
- Moved .SUPSYM and .THRSSY options to \*OPTIMIZATION from \*ORBITAL INPUT; updated Sirius part of manual (H. J. Aa. Jensen)
- Removed .NOSUPSYM option (H. J. Aa. Jensen)
- Renamed .PEQM and \*PEQM to .PELIB and \*PELIB, respectively. The former is still allowed but is deprecated.


## [2018.2] (2019-03-17)

### Fixed
- Fixed error in AO-direct CC3 response calculations causing segmentation faults
- Fixed calculation of DSO contribution to spin-spin coupling for MCSCF and HSROHF when no symmetry
- Dalton script: do not set OMP\_NUM\_THREADS=1 if not MPI parallel (better performance
  for sequential calculations if threaded blas is used, e.g. MKL or openBLAS)
- Dalton script: stop if user asks for MPI run with a sequential dalton.x
- More robust .STEX input specification (old failed in some situations with gfortran 8); changed documentation accordingly

### Added
- Dalton script: -gb and -ngb options for specifying work memory in gigabytes
- Dalton script: -np as an alternative to -N for specifying number of MPI nodes

## [2018.1] (2019-01-14)

### Fixed
- Error in code for 2-el integral transformation level 4 (used in some cases for MCSCF). Error was not in Dalton2016.
- Error in export for FDE fixed.
- Compilation with PGI compilers now possible (but code with pelib or qfit using gen1int is not working).

## [2018.0] (2018-11-19)

### New features added
- New parallel 2-electron integral transformation .FCKTRA  (H. J. Aa. Jensen).
- Triplet excitation energies and polarizabilities with the AO-SOPPA code (P. A. B. Haase).
- Analytical PE-HF/DFT molecular gradients (see e.g. test/pehf\_geoopt for example of input).
  - Reference: N. H. List, M. T. B. Beerepoot, J. M. H. Olsen, B. Gao, K. Ruud, H. J. Aa. Jensen, and J. Kongsted. J. Chem. Phys. 142, 034119 (2015).
- Freezing atoms in geometry optimization (N. H. List and H. J. Aa. Jensen).
  (See e.g. test/geoopt\_freeze for example of input.)
- Effective external field (EEF) for one- and two-photon absorption in PE-HF/DFT calculations.
  - Reference: N. H. List, H. J. Aa. Jensen, and J. Kongsted. Phys. Chem. Chem. Phys. 18, 10070 (2016).
- Remove the most diffuse virtual orbitals after SCF or MCSCF (new .VIRTRUNC option).
- Add purely classical multipole-multipole interaction energy in PE-QM calculations (which can be skipped using the .SKIPMUL keyword under the \*PEQM section).
- Add basic frozen density embedding (FDE) functionality (A. Gomes, C. Jacob, L. Visscher).
- Dipole velocity complex linear polarizability with test rsp\_cpp\_veloci (N. H. List).
- Resonant-convergent (damped) cubic response at HF/DFT levels (T. Fahleson and P. Norman)
  - Reference: T. Fahleson and P. Norman. J. Chem. Phys. 147, 144109 (2017).

### Fixed
- Nuclear model keyword .NUCMOD was ignored, now the Gaussian nuclear model used in the Dirac program can be used in Dalton.
- Open-shell DFT is not implemented for many derivative properties in \*\*PROPERTIES, dalton now quits.
- Bugfix for .MNF\_SO (mean-field spin-orbit, AMFI) when basis set has big exponents (>10^9).
- Open-shell doublet ROKS DFT geometry optimization.
- Fix of .GSPOL for parallel PE-QM quadratic response.
- Corrected text about elimination of some two-photon transitions between excited states
  because they were duplicates (text had "Third order" instead of "Second order").
- Bugfixes for CC2 in environments (incl. external fields) by Ove Christiansen.
- Bugfix for .TDA for MCSCF (i.e. zero B matrix in linear response in \*\*RESPONSE).
- Bugfix for .GASCI after .HSROHF
- Fix of several library basis sets that were not read correctly for some atoms, which caused Dalton to abort.
- Now a .G-TENSOR calculation in \*\*RESPONSE:\*ESR module with a CI wave function does not abort.

### Changed
- OK to run ECD or OECD with SOPPA.
- More documentation of .STEX in manual.
- Default induced-dipole solver in polarizable embedding (through .PEQM keyword) is changed to JI/DIIS method, which improves parallel scaling performance, and default convergence threshold for induced dipoles is changed $`1.0\cdot10^{-8}>|\mu^{[k]}-\mu^{[k-1]}|`$ where $`\mu`$ is a vector containing all induced dipoles and $`k`$ is the iteration index.
- Minimum CMake version is now v3.1.

### Deprecated
- Environment variable DALTON\_NUM\_MPI\_PROCS is deprecated and will be removed in future releases, use DALTON\_LAUNCHER instead


## [2016.2] (2016-07-12)

### Added
- Added and documented Basis=INTGRL option for ATOMBASIS in .mol file.
- Included Be in cc-pV5Z basis set

### Fixed
- More robust code for reading exponents and contraction coefficients in Dalton-type basis set files, incl. such files from EMSL
- Work-around for Intel 15 compiler I/O problem in some response calculations
- Fix for spin-orbit coupling (SOC) between S/T excited states of same symmetry (problem reported on daltonforum.org)
- Further fixes of MCSCF in \*\*PROPERTIES for more than 255 basis functions - hopefully it is OK now for all requests.
- Fixed an error in the manual for spin-dipole (problem reported on daltonforum.org)
- Fix of open-shell Hartree-Fock occupation output (only output, not the calculation, was wrong if ROHF was followed by MCSCF)
- Fix of Douglas-Kroll post-SCF with less than 256 contracted basis functions, but more than 255 uncontracted basis functions
- Fix of an insufficient memory error for construction of 2-el. integrals in Dirac format with more than 255 basis functions
- Removed OpenACC CMake variable (currently no OpenACC directives in Dalton).


## [2016.1] (2016-04-07)

### Added
- Possibility to read basis set files as made by the EMSL web site
  (this makes it possible to also read basis set files in basis/ based on
  emsl output as e.g. aug-pcseg-1; only LSDALTON has been able to read them so far)

### Fixed
- MCSCF in \*\*PROPERTIES for more than 255 basis functions (fixes problem with MCSCF shielding reported on daltonforum.org)
- Make sure molecule is not moved in ADDSYM during numerical differentiation
- Fixed error in the printing of the cpu/wall time used in Sirius
- Fixed error in PBEc functional: gave NaN when rho was zero.
- Polished some format statements to reduce number of compiler warnings
- Fixed error in memory addressing for MCSCF g-tensor calculations
- Fixed 2 errors in author list for WIRE dalton publication in dalton output
- Removed unsupported configure options.


## [2016.0] (2015-12-22)

### Changed
- Separated Dalton and LSDalton.

### Added
- New faster CC3 module (Main authors: Rolf H. Myhre and Henrik Koch)
- Parallel AO-SOPPA (Main authors: Frederik Beyer Kjær Hansen and Rasmus Faber)
- QFIT library - fitting charges and, if desired, dipoles to simulate the electrostatic potential from a QM wave function: (Main author: Casper Steinmann)
- Vibrational averaging of NMR coupling constants at SOPPA/MP2 level (Main author: Rasmus Faber)

### Fixed (since version 2015.1)
- Ahlrichs-TZV basis: Fixed error in an exponent for Boron.
- ANO-RCC basis: Fixed Carbon basis set (wrong contraction coefficients, see [MOLCAS ANO-RCC](http://www.molcas.org/ANO/).
- ANO-RCC basis: Modified the 3 Th h-functions by replacing them with the 3 Ac h-functions to Th.
                 (A mistake was made in the original work when the 3 Th h-functions were made,
                  and this has never been redone. They are too diffuse, exponents
                  (0.3140887600, 0.1256355100, 0.0502542000) for Th, Z=90, compared to
                  (0.7947153600, 0.3149038200, 0.1259615200) for Ac, Z=89, and
                  (0.8411791300, 0.3310795400, 0.1324318200) for Pa, Z=91.
                  We have selected to just replace the 3 Th h-functions with those from the Ac basis set,
                  because the Ac g-functions are quite close to the Th g-functions, closer than Ac g-functions,
                  and therefore differences in results compared to optimized Th h-functions should be minimal.)
                  Thanks to Kirk Peterson for pointing out the Th problem on http://daltonforum.org.
- Fixed reading of ANO-RCC basis set library file.
- Bug fix for when more than 30 excitation energies requested (EIGENVALUES NOT PAIRED problem reported by Frank Jensen).
- Fixed some bugs for two byte packing of derivative and spin-orbit two-electron integrals.
- Fixed .NEWTRA integral transformation for 32 bit integers and exactly n\*256 orbitals and no integer overflow test
  (the first 32 bits of (n\*256)\*\*4 are zero !!!).
- Improved performance of .NEWTRA integral transformation for response calculations.
- Do not include floating orbitals in calculation of smallest atom-atom distance.
- Enable Tamm-Dancoff approximation (.TDA) for embedding models, e.g. PE, PCM etc.
- Provide date and time stamp also for Darwin (i.e. MacOSX).
- Assume nobody uses gfortran version 4.0.2 any more (removed special test for that).


## [2015.1] (2015-07-20)

### Common
- Added ANO-RCC basis set.

### DALTON
- Fixed a bug in an LRESC correction.
- Improved calculation of one LRESC correction.
- Update PElib (v.1.2.3): Workaround for faulty system detection using Macports CMake
- Fixed a bug with Intel Compiler 15 during initialization of Cauchy-Schwarz parameters
- Fixed a bug for parallel build on some systems
- Fixed a segmentation fault for approx. 3000 basis functions
  (an array was allocated in stack memory, and became too big for default size of stack memory).
- Fixed a bug in CC sumrules.
- Fixed a bug in the preoptimization, i.e. when using smaller basis sets first to geometry optimize molecule.
- Fixed some far from optimal defaults for preoptimization.
- Fixed geometry optimization for HS-ROHF and HS-RODFT with symmetry - .SINGLY input option under '\*SCF INPUT'
  (use numerical gradients as analytical gradients are only implemented without symmetry).
- Some minor corrections to the Dalton manual.
- More reasonable output for TPCD.

### LSDALTON
- Fixed .UNCONT for EMSLs Dalton basis set format.


## [2015.0] (2015-02-18)

### Common
- Read EMSL Dalton format
- New included basis sets: pc-seg by Frank Jensen
- Many small improvements here and there

### DALTON
- CPP (Complex Polarization Propagator): MChD, NSCD
- Extensions of the polarizable embedding model (QM/MM via PE library):
  - PE-MCSCF wave function and linear response
  - PE-CPP damped linear response
  - PE for magnetic linear response using London atomic orbitals (LAOs)
- PCM-SOPPA excitation energies
- QFIT (electrostatic potential fitted charges)
- QM/CMM approach
- DFT-D3 and DFT-D3(BJ) dispersion energy corrections

### LSDALTON
- Geometry optimizations:
  - Quasi-Newton transition state optimization
  - HOPE algorithm
- Automated Counterpoise corrected DFT, HF, DEC-MP2 and CCSD interaction energies
- Dynamics: Nose-Hoover thermostat
- DFT-D3 and DFT-D3(BJ) dispersion energy corrections
- Performance improvements:
  - Charge-constrained ADMM exchange (energy+gradients)
  - Optimized Complex Polarization Propagator (CPP) solver for LSresponse
- Quadratic Response calculation to compute full dipole moment matrices in LSDalton


## [2013.4] (2014-07-10)

### DALTON
- Memory bugfix for serial PCM calculations (segmentation fault for large PCM cavities).

### LSDALTON
- Fixed a bug in the basis set reading. This bugfix affects almost no basis sets,
  and none of the standard basis sets, but a very few general contracted basis sets
  where the first contracted function had much smaller number of
  primitives compared to the last: Basis sets such as the pcS-1 basis set.
- Reduced the memory requirements for internal MPI buffer handling.


## [2013.3] (2014-06-11)

### Common
- aug-cc-pVTZ-lresc basis set added to $BASDIR.

### DALTON
- Default DIIS space increased from 5 to 8, often resulting in 1-2 fewer SCF iterations.
- Removed the maximum of 20 excitations in summary output for second and third order transition moments.
- Warning is issued when orbitals are deleted due to linear dependencies (before SCF),
  AngPso (a 0th order LRESC diamagnetic corr) is not calculated in this case.
- Bugfix for parallel calculations and some type of geometry optimizations with ANO basis sets
  (this bug resulted in aborted calculations, not in wrong results).
- Print irrep names together with symmetry numbers for easier interpretation of output.
- More important output with '@' in column 1 (can be obtained with 'grep @' on the output).
- Environment variable DALTON\_USE\_GLOBAL\_SCRATCH disables copying of binaries to worker nodes.
- Environment variable DALTON\_LAUNCHER introduced.
- Fixed output information about number of MPI processes and number of OpenMP threads.
- Added information in the error messages when values in maxorb.h are exceeded (which values to increase).
- Increased some of the values in the common blocks:
  MXSHEL 1000 -> 1500; MXCORB 2400 -> 5000; MXPRIM 8000 -> 15000;
  MAXOCC 800 -> 1500; MXCENT 200 -> 500; MXCENT\_QM 200 -> 500
  (the static size of dalton.x went from 100 MB to 165 MB).
- Do not print garbage non-zero transition moments and oscillator strengths for triplet excitations (\*EXCITA module).
- Corrected input description for transition moments between excited states (\*QUADRA with .DOUBLE RESIDUE).
- Fix for \*\*RESPONSE .EXMOM .ISPABC=1,0,1 (only half the excited state spin-orbit transition moments were calculated).
- Fix for Molden file when exponent greater than 1.0D8.
- Fix for MNF-SO (amfi) if more than 40 nuclei.
- Bugfix in quadratic response function using CPP in the tensor contraction routine of the A[2] terms.
- Added interface to ChemShell.
- Bugfix for small non-default WORK array sizes. For specific small custom values of the WORK array size
  KBLOCK was larger than MXBLCK leading to unpredictable results due to array length mismatch in DALTON/abacus/herrdn.F.

### LSDALTON
- Environment variable LSDALTON\_LAUNCHER introduced.


## [2013.2] (2014-03-05)

### Common
- Recognize CYGWIN as a LINUX and UNIX system, for proper definition of compilation flags.
- Define M\_PI in C-code if not already defined (problem seen with Cygwin).
- Added setup option --blacs to be used in combination with --scalapack; defaults to --blacs=intelmpi.

### DALTON
- Fixed a bug in printing results in CPP-QRF.
- New CPP solver works also for non-direct calculation.
- More efficient evaluation of numerical Hessian when C1 symmetry
  (in each geometry step start wave function optimization from a
  converged wave function from a neighboring geometry rather than from scratch each time).
- Fix of error which sometimes caused a geometry optimization to stop with "\*\*\* ERROR, Wrong interval in WLKBIS".
- Fix of a bug which occasionally caused DALTON to abort a .STEX calculation.
- Print final geometry in xyz format (angstrom). File called "final\_geometry.xyz" is put into the restart tarball.
- Append PID to scratch directory to avoid multiple tests running in the same directory.
- Improved manual for two-photon and non-adiabatic coupling.
- Updated/corrected g-factors for Ag, Nd, and Tl (thanks to M. Jaszunski).

### LSDALTON
- Print sensible error message when running out of memory.
- Added functionality to search through several basis-set libraries.
- Increased max length of WRKDIR from 60 to 200.
- Fixed a bug related to improper shutdown of MPI calculation. In the case
  of wrong LSDALTON.INP for instance the calculation will issue a error
  statement and afterward hang forever in a MPI call.
- Fixed an OpenMP bug in the calculation of how much memory there should be used during
  an exchange-correlation calculation - resulting in huge memory usage for large molecular system.


## [2013.1] (2013-12-19)

### DALTON
- Correct the printout of relativistic corrections to the shielding (thanks to M. Jaszunski).
- Compilation fix for DALTON/abacus/rma\_windows.F90 (Intel 10.0.011).
- Fix of error where basis set names were changed to upper case and could not be found (reported by Yurij Rusakov).
- Each MPI slave sleeps 10 millisecond between tests for new task
  (only Intel; should enable turbomode in sequential parts of DALTON, and more efficient use of threaded MKL when combined with MPI).
- added metric scaled output of orbital response vectors in \*\*RESPONS
  (for easier interpretation of excitation operators).

### LSDALTON
- Fixed a bug in Jengine, related to screening for nonsymmetric density matrices.
  This may affect CCSD and some response calculation.
- Modified the input section of the manual concerning
  Casida-Salahub asymptotic correction CS00 (thanks to Raul Crespo).
- Changed defaults for Casida-Salahub asymptotic correction CS00 (thanks to Raul Crespo).
- Fixed errors in the MCD B terms output files (.dat files) now one file is generated
  for each B term and each A term (thanks to Raul Crespo).
- Modified the input section of the manual concerning MCD B terms. Added description of MCDEXSTATES.
- Fixed a bug for LSDALTON geometry optimization and dynamics related to
  screening. The initial Cauchy-Schwartz screening matrices were incorrectly
  used in each subsequent geometry step


## [2013.0] (2013-11-11)

### DALTON
- Subsystems CC using Cholesky decomposition
- Multiscale modeling using the Polarizable Embedding (PE) library.
- Static exchange (STEX) for X-ray spectroscopy
- Damped response via Complex Polarization Propagator (CPP)
- Quadratic response for open-shell DFT
- Relativistic corrections to nuclear shielding constants
- Empirical dispersion corrections DFT-D2, DFT-D3 and DFT-D3BJ
- Various performance improvements and a few bug fixes

### LSDALTON
- Dynamics using HF and DFT
- More response properties, including some magnetic properties like MCD
- DEC-MP2 energy, density and gradient
- Local orbitals
- Improved SCF optimization routines
- Massively parallel CCSD
- MPI parallelization of HF and DFT
- Improved integral code including MPI parallelism
- Matrix operation parallelization using PBLAS/SCALAPACK
- ADMM exchange (energy, gradients)

# Folder description

Scripts to generate srDFT functionals for Dalton.

* functionals.py; normal DFT functionals.
* functionals_mu.py; range-separated functionals.
* functionals_special.py; special case of some functionals. Might be changed later, to just be directly embedded into the functionals that include these special cases.
* print_functional_to_DALTON.py; Script to write functional and it derivatives to Dalton format.
* print_functional_to_DALTON_TPSSc_special_case.py; Script to write TPSSc to Dalton format, special case because of MAX() function. This should be changed later for shorter code and faster code generation.
* constants.txt; Constants for the functionals.
* constants_libxc.txt; Constants for the functionals, that are identical to those used in libxc.
* write_dalton_file.py; This is where all the ugly stuff is written in!


# Example of implementing code-generation of a functional (spin-PBE)

```python
E = mu_func.PBEc_mu(parameters).subs({rho_s: 0})*rho_c
d1E_rhoc = E.diff(rho_c)
d1E_gammacc = E.diff(gamma_cc)
d2E_rhoc2 = d1E_rhoc.diff(rho_c)
d2E_rhocgammacc = d1E_rhoc.diff(gamma_cc)
d2E_gammacc2 = d1E_gammacc.diff(gamma_cc)

Kernel = [E]
description = ["Implemented by E.R. Kjellgren.\n"]
diff_order = [1]
diff_idx = [[0]]
dalprint.dalton_functional_printer(Kernel, "ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[daltonfile, f2pyfile])
# # #
Kernel = [E,d1E_rhoc,d1E_gammacc]
description = ["Implemented by E.R. Kjellgren.\n"]
diff_order = [1,2]
diff_idx = [[0],[1,3]]
dalprint.dalton_functional_printer(Kernel, "D1ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[daltonfile, f2pyfile])
# # #
Kernel = [E,d1E_rhoc,d1E_gammacc,d2E_rhoc2,d2E_rhocgammacc,d2E_gammacc2]
description = ["Implemented by E.R. Kjellgren.\n"]
diff_order = [1,2,3]
diff_idx = [[0],[1,3],[1,4,6]]
dalprint.dalton_functional_printer(Kernel, "D2ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[daltonfile, f2pyfile])
```

## Specify the functional and create the derivatives. (having put the functionals into the right python functional files (just open the files to see how that works))

```python
E = mu_func.PBEc_mu(parameters).subs({rho_s: 0})*rho_c
d1E_rhoc = E.diff(rho_c)
d1E_gammacc = E.diff(gamma_cc)
d2E_rhoc2 = d1E_rhoc.diff(rho_c)
d2E_rhocgammacc = d1E_rhoc.diff(gamma_cc)
d2E_gammacc2 = d1E_gammacc.diff(gamma_cc)
```

## Specify all the values for the dalton_printer_function.

```python
Kernel = [E,d1E_rhoc,d1E_gammacc,d2E_rhoc2,d2E_rhocgammacc,d2E_gammacc2]
description = ["Implemented by E.R. Kjellgren.\n"]
diff_order = [1,2,3]
diff_idx = [[0],[1,3],[1,4,6]]
dalprint.dalton_functional_printer(Kernel, "D2ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[daltonfile, f2pyfile])
```

* Kernel; a list of all the terms that is generated. Have to be order by differentiation order.
* description; some notes about the functional if wanted.
* diff_order; specify how many terms there is at every order of differentiation. [1,2,3] means one term for zeroth order. Two terms for first order. Three terms for second order.
* diff_idx; a list of lists of indexes to tell where to put the derivatives inside the Dalton vectors. [[0],[1,3],[1,4,6]] mean that for the zeroth order derivative it is put into place 0, zero meaning it is a scalar. [1,3] means that the first derivatives are put into index 1 and index 3 in the d1E() inside Dalton. 
* shortrange=True; specify that the functionals should have "mu" in their input. 
* output_files=[daltonfile, f2pyfile]; specify what files to write the output to.

The "diff_order" and "diff_idx" specification is a little tedious, but is needed because it is hard to know Python side what is supposed to be put where inside Dalton.




Where does CMake search math libraries if you specify --blas/lapack=auto?
-------------------------------------------------------------------------

CMake will look in the environment variable MATH_ROOT.

For instance my .bashrc contains::

  source /opt/intel/bin/compilervars.sh intel64
  export MATH_ROOT=/opt/intel/mkl


Order of math libraries
-----------------------

Order is set by MATH_LIB_SEARCH_ORDER in MathLibs.cmake.
You can override this order by setting BLAS_TYPE and/or LAPACK_TYPE
for example to ATLAS or some other library that you prefer.


What to edit if your math library is not found although you have set MATH_ROOT?
-------------------------------------------------------------------------------

Normally you only need to edit MathLibs.cmake to add new libraries
or edit existing ones.

Since a vendor can provide libraries with different "fingerprints"
(example MKL), you can define different combinations (up to 9), for instance::

  set(MKL_BLAS_LIBS  ...)
  set(MKL_BLAS_LIBS2 ...)
  set(MKL_BLAS_LIBS3 ...)
  set(MKL_BLAS_LIBS4 ...)
  set(MKL_BLAS_LIBS5 ...)

Then CMake will first try MKL_BLAS_LIBS, then MKL_BLAS_LIBS2, etc.
The first pattern that will match will be linked against.
