Quantum ESPRESSO GPU
====================

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

This repository also contains the GPU-accelerated version of Quantum ESPRESSO.

Installation
============

This version requires the nvfortran (previously PGI) compiler from the
freely available NVidia HPC SDK. You are advised to use the most recent
version of NVidia software you can find. While any version later than 17.4
should work, many glitches are known to exist in older versions. 
The `configure` script checks for the presence of the nvfortran compiler and 
of a few cuda libraries. For this reason the path pointing to the cuda toolkit
must be present in `LD_LIBRARY_PATH`.

A template for the configure command is:

```
./configure --with-cuda=XX --with-cuda-runtime=YY --with-cuda-cc=ZZ --enable-openmp [--enable-openacc] [ --with-scalapack=no ]
```

where `XX` is the location of the CUDA Toolkit (in HPC environments is 
generally `$CUDA_HOME`), `YY` is the version of the cuda toolkit and `ZZ`
is the compute capability of the card. You can get those numbers from
command `nvaccelinfo`, if you have a properly configured HPC SDK:
```
$ nvaccelinfo | grep -e 'Target' -e 'Driver'
CUDA Driver Version:           11000
Default Target:                cc70
...
```
The version is returned as (1000 major + 10 minor). For example, CUDA 9.2 
would be represented by 9020. For the above case, configure QE with:
```
./configure --with-cuda=$CUDA_HOME --with-cuda-cc=70 --with-cuda-runtime=11.0
```
Alternatively, you may use the (deprecated) tool `get_device_props.py` in
directory `dev-tools/`.

It is generally a good idea to disable Scalapack when running small test
cases since the serial GPU eigensolver outperforms the parallel CPU
eigensolver in many circumstances.

From time to time PGI links to the wrong CUDA libraries and fails reporting a 
problem in `cusolver` missing `GOmp` (GNU Openmp). This problem can be solved
by removing the cuda toolkit from the `LD_LIBRARY_PATH` before compiling.

Serial compilation is also supported.

Execution
=========

By default, GPU support is active. The following message will appear at
the beginning of the output

```
     GPU acceleration is ACTIVE.
```

GPU acceleration can be switched off by setting the following environment
variable:

```
$ export USEGPU=no
```


Testing
=======

The current GPU version passes all tests with both parallel and serial 
compilation.
# How to contribute
You can contribute to this project, even as an ordinary user, by:

1. subscribing and partecipating to the [users' mailing list](https://lists.quantum-espresso.org/mailman/listinfo/users), answering other people's questions
and doubts;

2. suggesting improvements and reporting bugs, either to the [developers mailing list](https://lists.quantum-espresso.org/mailman/listinfo/developers),
or using the [Issue section](https://gitlab.com/QEF/q-e/issues) of the gitlab repository.

Information on the mailing lists and on how to report bugs can be found in the
[Quantum ESPRESSO web site](https://www.quantum-espresso.org), "Contacts" section.

You can contribute even more by:

3.  preparing new tests for the test suite: see the
[Wiki page](https://gitlab.com/QEF/q-e/wikis/Developers/Test-suite-and-test-farm);

4.  improving the documentation;

5.  porting to new/unsupported architectures or configurations;

6.  implementing new features, or improving existing ones.

For more information, consult the documentation for developers in the
[Wiki](https://gitlab.com/QEF/q-e/wikis/home)
![q-e-logo](logo.jpg)

This is the distribution of the Quantum ESPRESSO suite of codes (ESPRESSO:
opEn-Source Package for Research in Electronic Structure, Simulation, and
Optimization)

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## USAGE
Quick installation instructions for the impatient. Go to the directory 
where this file is. Using "make"
(`[]` means "optional"):
```
./configure [options]
make all
```
"make" alone prints a list of acceptable targets. Optionally,
`make -jN` runs parallel compilation on `N` processors.
Link to binaries are found in bin/.

Using "CMake" (v.3.14 or later):

```
mkdir ./build
cd ./build
cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc [-DCMAKE_INSTALL_PREFIX=/path/to/install] ..
make [-jN]
[make install]
```
Although CMake has the capability to guess compilers, it is strongly recommended to specify
the intended compilers or MPI compiler wrappers as `CMAKE_Fortran_COMPILER` and `CMAKE_C_COMPILER`.
"make" builds all targets. Link to binaries are found in build/bin.
If `make install` is invoked, directory `CMAKE_INSTALL_PREFIX`
is prepended onto all install directories.

For more information, see the general documentation in directory Doc/, 
package-specific documentation in \*/Doc/, and the web site 
http://www.quantum-espresso.org/. Documentation for developers 
can be found on [Wiki page on gitlab](https://gitlab.com/QEF/q-e/-/wikis/home).

## PACKAGES

- PWscf: structural optimisation and molecular dynamics on the electronic ground state, with self-consistent solution of DFT equations;
- CP: Car-Parrinello molecular dynamics;
- PHonon: vibrational and dielectric properties from DFPT (Density-Functional Perturbation Theory);
- TD-DFPT: spectra from Time-dependent DFPT;
- HP: calculation of Hubbard parameters from DFPT;
- EPW: calculation of electron-phonon coefficients, carrier transport, phonon-limited superconductivity and phonon-assisted optical processes;
- PWCOND: ballistic transport;
- XSpectra: calculation of X-ray absorption spectra;
- PWneb: reaction pathways and transition states with the Nudged Elastic Band method;
- GWL: many-body perturbation theory in the GW approach using ultra-localised Wannier functions and Lanczos chains;
- QEHeat: energy current in insulators for thermal transport calculations in DFT.

## Modular libraries
The following libraries have been isolated and partially encapsulated in view of their release for usage in other codes as well:

- UtilXlib: performing basic MPI handling, error handling, timing handling.
- FFTXlib: parallel (MPI and OpenMP) distributed three-dimensional FFTs, performing also load-balanced distribution of data (plane waves, G-vectors and real-space grids) across processors.
- LAXlib: parallel distributed dense-matrix diagonalization, using ELPA, SCALapack, or a custom algorithm.
- KS Solvers: parallel iterative diagonalization for the Kohn-Sham Hamiltonian (represented as an operator),using block Davidson and band-by-band or block Conjugate-Gradient algorithms.
- LRlib: performs a variety of tasks connected with (time-dependent) DFPT, to be used also in connection with Many-Body Perturbation Theory.
- upflib: pseudopotential-related code.

## GPU-enabled version
Since Feb.2021 this repository also works for GPU's (currently only NVIDIA). See file [README_GPU.md](README_GPU.md).

## Contributing
Quantum ESPRESSO is an open project: contributions are welcome.
Read the [Contribution Guidelines](CONTRIBUTING.md) to see how you
can contribute.

## LICENSE

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
##NOTES##

* Routines or utilities existing in two copies, one for QE and one for upflib:
  - randy
    in module uspp, file uspp.f90
  - invmat
    simplified version, in module upf_invmat, file upf_invmat.f90
  - capital, lowercase, isnumeric, matches, version_compare
    in module upf_utils, file upf_utils.f90
  - errore and infomsg
    as upf_error, in file upf_error.f90

* Module variables that have been (partially) duplicated:
   - kinds      => upf_kinds  (only dp)
   - constants  => upf_const
   - nsp in ions_base points to nsp in uspp_param

* TO BE DONE: 
  - set the correct value of nsp in uspp_param when allocate_uspp is called,
    use it ONLY inside upflib, remove link of nsp in ions_base to uspp_param
  - nh(:) is allocated in init_uspp_dims, but maybe it should allocated
    together with upf(:), when upf is read? Or maybe nh should be part of upf?
    It is used in many many places, though!
  - Merge pseudopotential_indexes from CPV/src/pseudopot_sub.f90 with the
    uspp initialization in upflib (init_us_1 etc); merge qvan2b and qvan2
    (requires merge of interpolation tables qrad and qradb)
  - upf_ions now contains just a function n_atom_wfc: move somewhere else?
  - upf_spinorb contains just two variables: merge into uspp? add to it
    the two functions spinor and sph_ind used only for spin-orbit?
  - lmaxq should be "the maximum value of L in Q functions", not "... + 1"
    and should be used to dimension arrays where l=0,...,L. The dimension
    of spherical harmonics (2*lmaxkb+1)^2 is something different and should
    be stored in a different variable (something like ylmdim, or maxlm)
  - Names of interpolation tables and related routines are random:
      CP     PW      new name? 	     contains 	         computed in
    betagx   tab     tab_beta      beta(G) functions	compute_betagx,
             tab_d2y tab_beta_d2y  splines for beta(G)	init_tab_beta
    dbetagx             	   dbeta(G)/dG 		compute_betagx
    qradx    qrad    tab_q         Q(G) for  USPP/PAW	init_tab_qrad
    dqradx              	   dQ(G)/dG  		compute_qradx
             tab_at  tab_atwfc     atomic R_nl(G)	init_tab_atwfc
                    (tab_atrho     atomic rho(G)	to be done)
                    (tab_vloc      local potential	to be done)

* upflib restructuring:
  - shall we keep just one src folder ? or structure it a bit more, such as

  upflib/baselib     all basic data structures and io  
  upflib/advlib      advanced initializations (init_us_0,1,2)
  upflib/tools       tools from upftools
  

# Library of pseudopotential code

This directory contains a library of pseudopotential-related code,
extracted from the Quantum ESPRESSO distribution. This library depends only
upon module mp.f90 of UtilXlib and upon a few modules and routines of devXlib;
upon a few LAPACK routines; requires a suitable `../make.inc` file in Makefile. 
Other than this, it can be independently compiled.

Currently, it includes
- basic definitions of the UPF (Unified Pseudopotential File) format
- basic I/O operations on UPF files
- setup of the interpolation tables and of other basic variables 
- interpolation of pseudopotentials
- generation of various pseudopotentials matrix elements
- utilities: spherical harmonics and Bessel functions, integration routines
Old UPF specifications can be found here:
http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format
The xml schema for the newer UPF definition can be found here:
http://www.quantum-espresso.org/ns/qes/qe_pp-1.0.xsd

In addition to the `libupf.a` library, executable utilities are produced:

- `upfconv.x`, converting pseudopotentials in other formats into UPF:
   see `upfconv.x -h` for more

- `virtual_v2.x`, courtesy Jingyang Wang (jw598@cornell.edu), generates
   an averaged pseudopotential suitable for Virtual Crystal Approximation

- `casino2upf.x`, courtesy Mike Towler (see below)

A python script `fixfile.py` is also present, to remove undesired `&`
characters from UPF files that hinder their parsing by xml tools.

## CASINO and QE pseudopotentials

The following notes are kept for reference (they might be obsolete).
Code `upfconv.x -c` should replace code `upf2casino2.x` mentioned below.
Code `casino2upf.x` was moved to upflib/ and works (?) again, at least
for the example provided by Jake Muff, since v.6.8. Old notes start here:

Two utilities are provided with the Quantum Espresso distribution to 
enable the PWscf code to be used in conjunction with the CASINO quantum 
Monte Carlo code.

Of course all pseudopotentials generated via these automatic tools should 
be tested before being used for production runs.

It should be noted that ultrasoft and PAW pseudopotentials cannot be used
with the CASINO code. Currently only UPF files containing norm-conserving 
pseudopotentials can be converted using these utilities.

### casino2upf.x

The first of these is casino2upf.x . This utility takes a given CASINO 
tabulated pseudopotential file and one or more awfn.data files specifying 
the pseudoatomic wavefunctions to be used in creating the 
Kleinman-Bylander projectors. A UPF file containing the projectors and the 
local potential is then written to the file name specified in inputpp. Any
errors are communicated to the user via stderr.

Usage:
	
        ./casino2upf.x < inputpp

A sample inputpp file for converting a Trail and Needs pseudopotential 
would be:

```
inputpp:
	&inputpp
		pp_data='pp.data'
		upf_file='my_pseudo_potential.UPF'
	/
	3
	awfn.data_s1_2S
	awfn.data_p1_2P
	awfn.data_d1_2D
```

Here pp_data specifies the name and location of the file containing the 
CASINO pseudopotential. The utility then expects an input card after 
&inputpp consisting of the number of awfn.data files supplied (in this 
case 3) and then their names. The files are searched sequentially so the 
first s wavefunction found will be used for the s projector, first p for 
the p projector and so on.


*A note on the radial grid*

The utility currently performs no interpolation and attempts to use the 
same radial grid as the original pseudopotential. It therefore assumes 
that the grid will be of the standard form used by Trail and Needs.

If this is not the case the flag tn_grid=.false. can be set in the input 
file. The standard logarithmic form, r(i)=exp(xmin + i*dx) / Z is then 
assumed. Values for xmin and dx can also be specified in the input file in 
the usual way.

If interpolation from a different non-standard grid is required then the 
current recommended route is to use the casino2gon utility supplied with 
the CASINO distribution. This produces the older GON format that is 
(currently) still read by PWscf.


*Ghost states*

The Kleinman-Bylander form can unfortunately introduce ghost states into 
some calculations. If this does occur we recommend that the 
pseudopotential is re-converted using a different local channel. The local 
channel can be specified in the original CASINO pp.data file and is read 
in automatically by casino2upf.x .

### up2casino.x

This utility takes a standard UPF pseudopotential from standard input and 
writes a CASINO tabulated pseudopotential file to standard output. Any 
errors are communicated via stderr.

Usage:
	
	./up2casino.x < pseudo.UPF > pp.data

Care must be taken that the resulting pseudopotential file spec fies the 
required local channel. Also this utility should only be used with 
norm-conserving pseudopotentials.


project: FFTXlib
project_dir: .
output_dir: ./Doc
predocmark: >
docmark_alt: #
predocmark_alt: <
display: public
         private
exclude: fft_scalar.DFTI.f90
         fft_scalar.ESSL.f90
         fft_scalar.FFTW3.f90
         fft_scalar.FFTW.f90
         fft_scalar.SX6.f90
         mpif.h
include: /cineca/prod/compilers/intel/cs-xe-2015/binary/impi_5.0.2/include64/
graph: true

A self-contained library for handling FFT in QE
# FFTXlib

Implements real space grid parallelization of FFT and task groups. 

## Testing and Benchmarking

This library also provides a testing and timing code to asses the performance of your FFT, estimate the
scalability and the optimal parameters for your simulation.

To compile the test program, once you have properly configure QE within a parallel environment,
go inside the directory FFTXlib and type:

    make TEST

Then you can run your FFT tests using command like:

    mpirun -np 4 ./fft_test.x -ecutwfc 80 -alat 20  -nbnd 128 -ntg 4

Command line arguments:

    -ecutwfc  Plane wave energy cut off
    -alat     Lattice parameter (for hard coded lattice structure)
    -nbnd     Number of bands (fft cycles)
    -ntg      Number of task groups
    -av1  x y z    First lattice vector, in atomic units. N.B.: when using -av1, -alat is ignored!
    -av2  x y z    Second lattice vector, in atomic units. N.B.: when using -av2, -alat is ignored!
    -av3  x y z    Third lattice vector, in atomic units. N.B.: when using -av3, -alat is ignored!
    -kmax kx ky kz    Reciprocal lattice vector inside the BZ with maximum norm. Used to calculate max(|G+K|). (2pi/a)^2 units.

A python script to extract the parameters from an output file of pw.x is also available. Example usage:

    $ python gen_test_params.py a_pw_output
    To analize performances run with:
    mpirun -np X ./fft_test.x -ntg Y -ecutwfc 36.7500 -ecutrho 147.0000 -av1 36.6048 0.0 0.0 -av2 -18.3024 31.70067192 0.0 -av3 0.0 0.0 18.3024 -nbnd 400 -gamma .true.

Replace `X` and `Y` with appropriate values for your simualtion.
    
## Files
Compile time parameters:

    fft_param.f90

Descriptor types:

    stick_base.f90
    fft_types.f90
    fft_smallbox_type.f90

Parallel execution routines:

    fft_interfaces.f90 fft_fwinv.f90
      fft_parallel.f90
        scatter_mod.f90
          tg_gather.f90
      fft_interpolate.f90
      fft_smallbox.f90

Low level library wrappers:

    fft_scalar.f90
    fft_scalar.DFTI.f90
    fft_scalar.ESSL.f90
    fft_scalar.FFTW.f90 fftw_interfaces.f90
      fft_stick.c fftw.c fftw_dp.c fftw_sp.c fftw_dp.h fftw.h fftw_sp.h konst.h
    fft_scalar.FFTW3.f90
    fft_scalar.SX6.f90

Misc. helper routines:

    fft_ggen.f90
    fft_error.f90
    fft_helper_subroutines.f90
    fft_support.f90

Tests:

    test0.f90
    test.f90

## Release checklist 

1. run all examples; ideally also check for discrepancies, or at least for crashes
2. verify that all README, README.md, etc. files contain updated information on the content of the relative package and directory
3. verify that all documentation files, in particular Doc/developer-man.tex and all user_guide.tex files, contain updated information. In particular, verify that there are no references to removed or obsolete software and no missing references to new or changed software.
4. update the release number in developer_man.tex and in all user_guide.tex and other documentation that contains references to version number
5. verify that input documentation (files INPUT_*.def) is updated
6. update Doc/release-notes with the release number and with updated information on what is new, changed, removed, etc.
7. Re-generate new documentation with "make doc"
8. verify that install/configure is updated and aligned with install/configure.ac
9. update version number in Modules/version.f90
10. set a git tag "qe-x.y[.z]" for version x.y[.z]
11. align master to develop, github to gitlab
12. make packages on gitlab and github
13. if there are changes to the schema, copy the new schema to 
quantumespresso@qe.safevps.it:/storage/vhosts/quantum-espresso.org/ns/qes
14. update the web site: add a piece of news, update pages Downloads, Roadmap, and any other page that needs to be updated, copy the new documentation to directory quantumespresso@qe.safevps.it:/storage/vhosts/quantum-espresso.org/htdocs/Doc
15. send a message to the mailing list, post to twitter, facebook, and whatnot
# Dev Tools

This directory contains several tools that may be useful for developers

- `mem_counter`. A script that tracks all calls to `allocate` and `deallocate`,
   appending a call to subroutine `UtilXlib/mem_counter.f90`.
   Calls python script `mem_counter.py`, written by Pietro Bonfà (CINECA)
   and improved by Samuel Poncé. 
   `mem_counter -h` gives information on how to use it.
   BEWARE: you may still need to manually edit some files that do not compile.
-  `mem_analyse.py` is a python script, by Samuel Poncé, locating memory leaks.
   See the script header for directions on how to use it.
- `src-normal`. A script that "normalizes" the fortran syntax to QE style
   (see below). Calls python script `src-normal.py`, written by Norbert Nemec.

   Usage: `src-normal file1.f90 [file2.f90 ...]` or `src-normal`
- Utilities for PWgui:
  * `check_gui` (called via `Makefile`)
  * `diff_gui_help`
  * `guihelp.xsl`
  * `update_gui_help`
- Utilities for helpdoc (see `README.helpdoc`):
  * `helpdoc`
  * `helpdoc.d`
  * `helpdoc.schema`
  * `input_xx.xsl`
- Utilities for emacs_mode:
  * `gen-emacs-mode`
  * `gen-emacs-mode.tcl`

## Obsolescent utilities
- GPU utilities by Pietro Bonfà:
  * `get_device_props.py`
  * `device_props.c`
- Other utilities:
  * `calltree.pl`
   A perl script, to be run from the root QE directory, producing in the
   standard output the tree of called routines
  * `callhtml.pl`
   As above, producing a html page with the tree of called routines

## Coding style
These are some basic rules for Fortran codes enforced by `src_normal`:
* Use spaces for indentation instead of tabs (tab width 8 characters).
* Trailing whitespaces at the end the line should be removed.
* Normalize multiword keywords (e.g. END DO).
* Use capitalize version of the intrisic keywords (IF, DO, SUBROUTINE, etc.).
* Use the newest version of the comparison operators (==, >, etc.) instead of the old one (.eq., .gt., etc.)
title: EPW
src_dir: ./src
output_dir: ./doc
project_website: http://epw.org.uk/
summary: EPW is the short name for "Electron-phonon Wannier". EPW is an open-source F90/MPI code which calculates properties related to the electron-phonon interaction using Density-Functional Perturbation Theory and Maximally Localized Wannier Functions.
authors: Samuel Poncé
         Roxana Margine
         Carla Verdi
         Feliciano Giustino
author_description: The EPW project is mainly developed at the university of Oxford.
github: https://github.com/sponce24
email: samuel.pon@gmail.com
predocmark: >
media_dir: ./media
page_dir: ./Ford
docmark_alt: #
predocmark_alt: <
display: public
         private
source: false
graph: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: GNU
extra_filetypes: sh #

EPW is the short name for "Electron-phonon Wannier". EPW is an open-source F90/MPI
 code which calculates properties related to the electron-phonon interaction
 using [Density-Functional Perturbation Theory](http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.73.515) 
and [Maximally Localized Wannier Functions](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.12847). 
EPW is developed and maintained by [Samuel Poncé](http://giustino.materials.ox.ac.uk/index.php/Site/SamuelPonc%e9), [Roxana Margine](http://www.binghamton.edu/physics/people/margine.html), [Carla Verdi](http://giustino.materials.ox.ac.uk/index.php/Site/CarlaVerdi), and [Feliciano Giustino](http://giustino.materials.ox.ac.uk/).

The reference technical manuscript for the latest pubic release is:
[EPW: Electron-phonon coupling, transport and superconducting properties using maximally localized Wannier functions](http://arxiv.org/abs/1604.03525)
by S. Poncé, E. R. Margine, C. Verdi, and F. Giustino. 


@Note
Since 26 April 2016 EPW is distributed as part of the [Quantum ESPRESSO](http://www.quantum-espresso.org/) suite. 

The code was written by Feliciano Giustino (EPW v1) and Jesse Noffsinger (EPW v2) while
 at the University of California, Berkeley. Brad Malone (Harvard) and Cheol-Hwan Park
 (Seoul National University) contributed with tests and benchmarks. 
Roxana Margine implemented the anisotropic Eliashberg theory while at the University of Oxford (EPW v3). 
Samuel Poncé (Oxford) made the code compatible with the latest version of Quantum Espresso v5 
in the latest release EPW v4. Carla Verdi (Oxford) developed the electron-phonon interpolation 
for polar materials including Froehlich correction (released within EPW v4).

EPW is based on the method introduced in F. Giustino et al, [Phys. Rev. B 76, 165108 (2007)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.165108). 
An extended description of the first public release has been published in 
J. Noffsinger et al, [Comput. Phys. Comm. 181, 2140 (2010)](http://www.sciencedirect.com/science/article/pii/S0010465510003218). The extension of EPW to include the
 anisotropic Midgal-Eliashberg theory is based on the method described in
 E. R. Margine et al, [Phys. Rev. B 87, 024505 (2013)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.024505). 
The latest release of the code is described in S. Poncé et al, [arXiv:1604.03525](http://arxiv.org/abs/1604.03525). 
title: EPW overview
author: Samuel Poncé
date: 01-06-2016


<div style="text-align:center"><img src ="http://epw.org.uk/figures/logo_v7.png" width="600"></div>


EPW is an open-source F90/MPI code which calculates properties related to the electron-phonon interaction
 using [Density-Functional Perturbation Theory](http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.73.515)
and [Maximally Localized Wannier Functions](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.12847).

EPW is licensed under a [[GNU General Public License]](http://www.gnu.org/licenses/gpl-3.0.en.html)

EPW is part of the [[Quantum Espresso]](http://www.quantum-espresso.org/) software package.

## Installation 

The EPW software is only tested and intened to run on Linux (Mac OS might work but not tested).

* Download the latest version of [[Quantum-ESPRESSO]](http://www.qe-forge.org/gf/project/q-e/frs/?action=FrsReleaseBrowse&frs_package_id=18).

* Unpack and configure Quantum-ESPRESSO
```bash
tar -xvf espresso-5.4.0.tar.gz && cd espresso-5.4.0 && ./configure
``` 

* Compile EPW (this will also compile pwscf, phonon, and wannier90)
```bash
make -j 4 pwall
make -j 4 ph
make -j 4 epw
```

* The executable will be available in espresso-5.4.0/bin/epw.x or espresso-5.4.0/EPW/bin/epw.x
title: Quantum Espresso PWSCF
src_dir: ./src/
output_dir: ./doc
include: ./../LAXlib
         ./../KS_Solvers
exclude: laxlib.fh
         ks_solver_interfaces.fh
project_website: http://www.quantum-espresso.org/
summary: Quantum Espresso is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials.
authors: Paolo Giannozzi et al. 
author_description: The Quantum Espresso project is developed by the QE community.
github: https://github.com/QEF/q-e
gitlab: https://gitlab.com/QEF/q-e
email: users@lists.quantum-espresso.org
predocmark: >
media_dir: ./media
page_dir: ./Ford
docmark_alt: #
predocmark_alt: <
display: public
         private
source: false
graph: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: GNU
extra_filetypes: sh #

[Quantum ESPRESSO](http://www.quantum-espresso.org/) is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale.
It is based on density-functional theory, plane waves, and pseudopotentials.

Quantum ESPRESSO has evolved into a distribution of inter-operable codes in the spirit of an open-source project. The Quantum ESPRESSO distribution consists of a core set of components, a set of additional packages performing more advanced tasks, and a number of third-party packages designed to be inter-operable with the core components. Researchers active in the field of electronic-structure calculations are encouraged to participate in the project by contributing their own codes or by implementing their own ideas into existing codes.

Quantum ESPRESSO is an open initiative, in collaboration with many groups world-wide, coordinated by the Quantum ESPRESSO Foundation. Present members of the latter include Scuola Internazionale Superiore di Studi Avanzati, the Abdus Salam International Centre for Theoretical Physics (Trieste), the CINECA National Supercomputing Center (Bologna), the Ecole Polytechnique Fédérale de Lausanne, the University of North Texas (Dallas), the Oxford University. Courses on modern electronic-structure theory with hands-on tutorials on the Quantum ESPRESSO codes are offered on a regular basis in collaboration with the Abdus Salam International Centre for Theoretical Physics in Trieste.

The reference technical manuscript for the latest public release is:
[P. Giannozzi et al., J.Phys.:Condens.Matter 29, 465901 (2017)](http://iopscience.iop.org/article/10.1088/1361-648X/aa8f79)


title: Quantum ESPRESSO overview
author: Samuel Poncé
date: 10-01-2017


<div style="text-align:center"><img src ="https://www.quantum-espresso.org/project/logos/Quantum_espresso_logo.jpg" width="600"></div>


[[Quantum Espresso]](http://www.quantum-espresso.org/) is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials.

QE is licensed under a [[GNU General Public License]](http://www.gnu.org/licenses/gpl-3.0.en.html)


## Installation 

Quantum ESPRESSO is currently distributed as source packages, but selected binary packages for Linux and Windows are also available. Since v.6.0 the distribution is split into sources, tests, examples, plus a few specialized packages and plug-ins. These packages can be automatically downloaded and installed on demand after the installation of the base distribution: see the general documentation for more information.


    Copyright (c) Targacept, Inc.

    Targacept, Inc., 
    200 East First Street, 
    Suite 300, 
    Winston-Salem, NC, USA 27101 
    atp@targacept.com

This file describes the Autopilot Feature Suite as introduced and used by 
Targacept, Inc.   This documentation accompanies free software; The software
is subject to the terms of the GNU General Public License as published by the 
Free Software Foundation; either version 2 of the License, or (at your option) 
any later version. See the GNU General Public License at 
www.gnu.or/copyleft/gpl.txt for more details.

This documentation, like the software it accompanies, is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY FOR A PARTICULAR PURPOSE.  

[[_TOC_]]

--------------------------------------------------------------------------------
AUTOPILOT DOCUMENTATION
--------------------------------------------------------------------------------

The Autopilot Feature Suite is a user level enhancement for directing 
Car-Parrinello simulations based on CP.X packaged in ESPRESSO. 

The following features are incorporated: 

 - Auto Restart Mode 
 - Autopilot Course Configuration (Dynamic Rules) 
 - Autopilot Course Correction (Steering) 


--------------------------------------------------------------------------------
Auto Restart Mode 
--------------------------------------------------------------------------------

Auto Restart Mode is an extension of restart_mode declared in the CONTROL section 
of the input file. When restart mode is set to "auto", control determines if the 
current run is "from_scratch" or a valid "restart" based on the presence of a 
restart file associated with unit NDR. When NDR, the unit number for input, and 
NDW, the unit number for output, are the same, a simulation that is system 
terminated can be restarted without significant loss, providing that ISAVE, the 
parameter that indicates the frequency at which intermediate data are saved, is 
not large. 

Auto Restart Mode implements an effective "upto" mode and is also designed for 
use on remote machines where simulations may frequently be terminated and 
restarted.  Auto Restart Mode is especially useful in connection with Autopilot's 
Dynamic Rules capability. When they are used together, only one segment of a 
simulation is necessary, thereby reducing run_script volume and errors, and 
placing more control with the user. 

    restart_mode   CHARACTER ( default = 'restart' )
           from_scratch = from scratch.  NEB only: the starting path is 
                             obtained with a linear interpolation between 
                             the images specified in the ATOMIC_POSITIONS 
                             card.  Note that in the linear interpolation,
                             periodic boundary conditions ARE NOT USED.
           restart  = continue a previous simulation and perform  
                             "nstep" new steps.
           reset_counters  = continue a previous simulation, perform  
                             "nstep" new steps, resetting the counter 
                             and averages.
           auto = automatically detect "from_scratch" or "restart"; 
                             continue any previous simulation, and stop 
                             when the counter value is equal to "nstep".




--------------------------------------------------------------------------------
Autopilot Course Configuration (Dynamic Rules) 
--------------------------------------------------------------------------------

Autopilot Course Configuration (Dynamic Rules) is a method that allows select 
input parameters (Autopilot variables) to change during the course of a 
simulation. This method allows the user to create a more concise set of 
instructions that are easier to read and maintain and enables a more continuous 
execution on remote resources. 

Typically and historically, a user issues a run_script that creates a sequence of
input files, each with fixed parameter values. This run_script then calls cp.x 
against each input file in the sequence, such that, after the first, each 
execution continues with the next input file as well as restart information from
the previous execution. 

The Autopilot Course Configuration effectively consolidates multiple input files 
into one, allowing the user to specify at what time step a parameter should change
along with its new value. Thus a run_script becomes much shorter, and the user can 
easily see the projected path of the simulation. 

The Autopilot Course Configuration feature is implemented by adding a new card type 
to the "CARDS" section of the input file. The Autopilot card must be placed after 
the "NAMELIST" section but otherwise may appear before or after any other card. 
A favorable place is as the first card. 

Sytnax is as follows: 

    CARDS ...

    AUTOPILOT

      optional card :  read dynamic rules to set parameters on an absolute
                       timestep (iteration) from either standard input or mailbox 
                       (pilot.mb)
      Syntax:

    AUTOPILOT
      ON_STEP = ith_event_STEP : varname = value  
      ON_STEP = jth_event_STEP : varname = value  
      ON_STEP = jth_event_STEP : varname = value  

    ...
      ON_STEP = nth_event_STEP : varname = value  
      ENDRULES

Description:

 - `ON_STEP` 	   LABEL, must be in numerical timestep order, otherwise rule 
    is ignored

 - `ith_event_STEP`   INTEGER, iteration (NFI) when rule is to be employed

 - `varname`	   Autopilot variable, currently limited to one of the 
    following: isave,iprint,dt,emass, electron_dynamics, 
    electron_damping, ion_dynamics, ion_damping, 
    ion_temperature, tempw.

 - `value`            Must be valid value of variable type 
    (for example: isave, iprint must have a value of type 
    INTEGER, while dt must have a value of type REAL)

 - `ENDRULES`         Required only for input (STDIN) if other cards follow.

The event specification (`ON_STEP`) should precede the variable assignment. The 
colon separator between the event assignment and the variable assignment is 
required, as are the equal signs. No semi-colon or comma should appear after the
variable assignment. There can be multiple rules per event but only one variable 
assignment per rule (and only one rule per line). Within one event, there should 
be only one assignment per variable.  If multiple assignments are made for the 
same variable for the same event, only the last assignment will be accepted. 
Rules for which event specifications are not in numerical order will be ignored. 
If syntax errors are present in the AUTOPILOT card during start-up, the 
execution will stop. 

Example Syntax: 

      AUTOPILOT
        ON_STEP = 200 : tempw = 500.0
        ON_STEP = 200 : dt = 3.0
        ON_STEP = 250 : ISAVE = 50
      ENDRULES

Currently there is a maximum of 32 supported events and 10 supported Autopilot 
variables. Events that are out of timestep order are ignored. A user may establish 
up to 10 rules (one for each Autopilot variable) per event. Currently implemented 
Autopilot variables are: isave, iprint, dt, emass, electron_dynamics, 
electron_damping, ion_dynamics, ion_damping, ion_temperature, and tempw. 
If desired, users may implement other Autopilot variables. See Appendix below for 
an explanation of "Adding an Autopilot Variable". 

IMPORTANT: Variables should have values in accordance with their TYPE, or a 
runtime error may occur.





--------------------------------------------------------------------------------
Autopilot Course Correction (Steering) 
--------------------------------------------------------------------------------

Autopilot Course Correction (Steering) provides a run-time method of changing 
Autopilot variables on the fly, after the simulation is underway. Autopilot 
Course Correction (Steering) can be applied through any of the following 
sub-features: New Course (power steering), Manual Steering, and Pause. 

Steering utilizes a new mailbox file:  pilot.mb. This file can be created via the
user's favorite text editor and can be "mailed" by placing the file in the 
"results" directory. The user can also very quickly implement a single course 
correction command with UNIX redirect to the pilot.mb file. 

When a pilot.mb mailbox file is detected, the current event table is cleared to 
prepare for the new course.  The mailbox file is then parsed, and Autopilot 
processes the command(s) before deleting the mailbox file. If Autopilot cannot 
parse a command, it issues a warning and goes into PAUSE mode (see below). 

The Steering subfeatures, including pilot.mb syntax are described here: 

 - New Course or 'power steering' is implemented with the same syntax as the 
   INPUT file card for Autopilot. Remember that ON_STEP represents an absolute 
   iteration (NFI) step. 

   For example: 

         AUTOPILOT                               -required
         ON_STEP=400 : ISAVE = 50                -events must be ordered by step
         ON_STEP=400 : DT = 5.0                  -use valid variable types (or die)
         ON_STEP = 600:IONS_TEMPERATURE='damped' -indention optional     
         ON_STEP = 600: TEMPW=350.0              -white spaces are ignored
         ENDRULES                                -optional

   In this example, when NFI reaches 400, the value of ISAVE will be reset to 50 
   and the value of DT to 5.0.  Then, when NFI reaches 600, IONS_TEMPERATURE and 
   TEMPW will be reset to the indicated values.

 - Manual Steering is implemented with a similar syntax except that the card type 
   is PILOT instead of AUTOPILOT and the user specifies a timestep relative to 
   the time the mailbox is read, rather than an absolute timestep. The relative 
   timestep allows the user to set a rule for a near future event without having 
   to judge the current absolute NFI value. The user may also pre-write multiple 
   mailboxes using relative event steps without regard to absolute iteration (NFI) 
   values.

   For example, assume mailbox contents are:

         NOW:ISAVE=50
         NOW+100:TEMPW=600.0

   Assume further that the mailbox is saved to the "results" directory and then 
   read when the NFI is 380.  Manual Steering will reset the value of ISAVE on 
   the next event that is modulo 50, and an ISAVE event will occur twice 
   (at 400 and again at 450) before TEMPW is reset to 600.0 on step 480. Compare 
   this with the syntax that specifies an absolute timestep:

         ON_STEP=400:ISAVE=50
         ON_STEP=500;TEMPW=600.0


   In this example, if the NFI is less than 400 when the mailbox is read, ISAVE 
   becomes 50 on step 400 and TEMPW becomes 600.0 on step 500, and ISAVE is 
   performed twice before TEMPW is reset, just as in the previous example that 
   uses relative indexing.  

   However, if the user misjudges the momentary NFI, and it is 530 when the 
   mailbox is read, then both rules are implemented immediately and simultaneously.
   Furthermore, the ISAVE rule takes effect after the NFI specified. Neither of 
   these effects may have been intended by the user.

   Following is an example of a Manual Steering mailbox to change temperature from 
   a relative iteration (NFI) step: 

   Example syntax for a Manual Steering mailbox is as follows: 

         PILOT                                -optional for single line
           NOW : ISAVE = 50                   -events must be ordered
           NOW : DT = 5.0                     -use valid variable types (or die)
           NOW+50 :IONS_TEMPERATURE='damped'  -offsets from NOW are supported  
           NOW + 150: TEMPW=350.0             -white spaces are ignored
         ENDRULES                             -optional

  Example format for a quick mailbox change using a single rule is as follows:

                                      	  -defaults to PILOT
     NOW + 250: TEMPW=450.0               -single line with NOW 

c) Pause is a steering sub-feature that allows the user to suspend the simulation
until the user can decide on a future course. Pause is very helpful when the user 
knows that a change should be imposed but needs time to establish rules and 
create an appropriate mailbox. Steering then resumes as AUTOPILOT or PILOT upon 
receiving another pilot.mb mailbox. The syntax is a single line with one of the 
following: 

      PAUSE
      SLEEP
      HOLD
      HOVER
      WAIT

All of the above perform the same PAUSE mechanism.  The user can issue the command
quickly through UNIX redirect:

  >echo "PAUSE" > results/pilot.mb

Any mailbox not correctly identified with a `AUTOPILOT`, `PILOT`, `NOW`, or a `PAUSE` 
command, will result in a warning to standard output (STDOUT), and the simulation 
will pause. 





--------------------------------------------------------------------------------
TESTING 
--------------------------------------------------------------------------------

The entire Autopilot Feature Suite issues directives to slave nodes under 
MPI, with the fewest broadcasted parameters. All features have been tested 
under Intel 8.1 with MKL 7.0.1 libraries on a Linux-32 single processor and under 
PGI 5.2 with MPI on Linux-64 with 1, 2 and 4 processors. 





--------------------------------------------------------------------------------
ADDING AN AUTOPILOT VARIABLE
--------------------------------------------------------------------------------

See `autopilot.f90` for examples.
- Select the input parameter from the list in file INPUT_CP 
- Identify parameter dependencies, initializations, assignments, etc 
- Edit autopilot.f90 to add the following, 
    where VARNAME is the name of the new Autopilot variable:
        VARTYPE :: rule_VARNAME(max_event_step) at module scope 
        LOGICAL :: event_VARNAME(max_event_step) at module scope
* Remember to add to the PUBLIC block as well
        event_VARNAME(:) = .false. to init_autopilot subroutine 
        rule_VARNAME(:) = VARDEFAULT to init_autopilot subroutine 
* Import VARNAME with USE to employ_rules subroutine
* In employ_rules, add conditional clause on event_VARNAME to assign VARNAME: 
         ! VARNAME
         if (event_VARNAME(event_index)) then
           VARNAME  = rule_VARNAME(event_index)
           CALL init_other_VARNAME_dependent_variables( VARNAME)
           write(*,*) 'RULE EVENT: VARNAME', VARNAME
         endif
* Import VARNAME with USE to assign_rule subroutine
* In assign_rule, add condition clause matching the VARNAME create rule as so: 
         ELSEIF ( matches( "VARNAME", var ) ) THEN
                     read(value, *) VARTYPE_value
                     rule_VARNAME(event)  = VARTYPE_value
                     event_VARNAME(event) = .true.
* TEST  

WARNING: Some Autopilot variables may create "side-effects".  For example, the 
inclusion of a rule for TEMPW rules invokes a side-effect call to ions_nose_init.  
The user is cautioned to be aware of possible side-effects when adding other 
Autopilot variables. 

[[_TOC_]]

Introduction
============

This guide covers the usage of the `CP` package, version 7.0, a core
component of the Quantum ESPRESSO distribution. Further documentation,
beyond what is provided in this guide, can be found in the directory
`CPV/Doc/`, containing a copy of this guide.

This guide assumes that you know the physics that `CP` describes and the
methods it implements. It also assumes that you have already installed,
or know how to install, Quantum ESPRESSO. If not, please read the
general User's Guide for Quantum ESPRESSO, found in directory `Doc/` two
levels above the one containing this guide; or consult the web site:
`http://www.quantum-espresso.org`.

People who want to modify or contribute to `CP` should read the
Developer Manual: `https://gitlab.com/QEF/q-e/-/wikis/home`.

`CP` can perform Car-Parrinello molecular dynamics, including
variable-cell dynamics. The `CP` package is based on the original code
written by Roberto Car
and Michele Parrinello. `CP` was developed by Alfredo Pasquarello (EPF
Lausanne), Kari Laasonen (Oulu), Andrea Trave, Roberto Car (Princeton),
Nicola Marzari (EPF Lausanne), Paolo Giannozzi, and others. FPMD, later
merged with `CP`, was developed by Carlo Cavazzoni (Leonardo), Gerardo
Ballabio (CINECA), Sandro Scandolo (ICTP), Guido Chiarotti, Paolo Focher,
and others. We quote in particular:

-   Sergio Orlandini (CINECA) for completing the CUDA Fortran acceleration
    started by Carlo Cavazzoni

-   Fabio Affinito and Maruella Ippolito (CINECA) for testing and benchmarking

-   Ivan Carnimeo and Pietro Delugas (SISSA) for further openACC acceleration

-   Riccardo Bertossa (SISSA) for extensive refactoring of ensemble dynamics /
    conjugate gradient part

-   Federico Grasselli and Riccardo Bertossa (SISSA) for bug fixes,
    extensions to Autopilot;

-   Biswajit Santra, Hsin-Yu Ko, Marcus Calegari Andrade (Princeton) for
    various contribution, notably the SCAN functional;

-   Robert DiStasio (Cornell)), Biswajit Santra, and Hsin-Yu Ko for
    hybrid functionals with MLWF; (maximally localized Wannier
    functions);

-   Manu Sharma (Princeton) and Yudong Wu (Princeton) for dynamics with
    MLWF;

-   Paolo Umari (Univ. Padua) for finite electric fields and conjugate
    gradients;

-   Paolo Umari and Ismaila Dabo (Penn State) for ensemble-DFT;

-   Xiaofei Wang (Princeton) for META-GGA;

-   The Autopilot feature was implemented by Targacept, Inc.

The original version of this guide was mostly written by Gerardo Ballabio
and Carlo Cavazzoni.

`CP` is free software, released under the GNU General Public License.\
See `http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt`, or the file
`License` in the distribution.

We shall greatly appreciate if scientific work done using the Quantum
ESPRESSO distribution will contain an acknowledgment to the following
references:

> P. Giannozzi, S. Baroni, N. Bonini, M. Calandra, R. Car, C. Cavazzoni,
> D. Ceresoli, G. L. Chiarotti, M. Cococcioni, I. Dabo, A. Dal Corso, S.
> Fabris, G. Fratesi, S. de Gironcoli, R. Gebauer, U. Gerstmann, C.
> Gougoussis, A. Kokalj, M. Lazzeri, L. Martin-Samos, N. Marzari, F.
> Mauri, R. Mazzarello, S. Paolini, A. Pasquarello, L. Paulatto, C.
> Sbraccia, S. Scandolo, G. Sclauzero, A. P. Seitsonen, A. Smogunov, P.
> Umari, R. M. Wentzcovitch, J.Phys.: Condens.Matter 21, 395502 (2009)

and

> P. Giannozzi, O. Andreussi, T. Brumme, O. Bunau, M. Buongiorno
> Nardelli, M. Calandra, R. Car, C. Cavazzoni, D. Ceresoli, M.
> Cococcioni, N. Colonna, I. Carnimeo, A. Dal Corso, S. de Gironcoli, P.
> Delugas, R. A. DiStasio Jr, A. Ferretti, A. Floris, G. Fratesi, G.
> Fugallo, R. Gebauer, U. Gerstmann, F. Giustino, T. Gorni, J Jia, M.
> Kawamura, H.-Y. Ko, A. Kokalj, E. Küçükbenli, M .Lazzeri, M. Marsili,
> N. Marzari, F. Mauri, N. L. Nguyen, H.-V. Nguyen, A. Otero-de-la-Roza,
> L. Paulatto, S. Poncé, D. Rocca, R. Sabatini, B. Santra, M. Schlipf,
> A. P. Seitsonen, A. Smogunov, I. Timrov, T. Thonhauser, P. Umari, N.
> Vast, X. Wu, S. Baroni, J.Phys.: Condens.Matter 29, 465901 (2017)

Users of the GPU-enabled version should also cite the following paper:

> P. Giannozzi, O. Baseggio, P. Bonfà, D. Brunato, R. Car, I. Carnimeo,
> C. Cavazzoni, S. de Gironcoli, P. Delugas, F. Ferrari Ruffino, A.
> Ferretti, N. Marzari, I. Timrov, A. Urru, S. Baroni, J. Chem. Phys.
> 152, 154105 (2020)

Note the form `Quantum ESPRESSO` (in small caps) for textual citations
of the code. Please also see other package-specific documentation for
further recommended citations. Pseudopotentials should be cited as
(for instance)

> \[ \] We used the pseudopotentials C.pbe-rrjkus.UPF and O.pbe-vbc.UPF
> from `http://www.quantum-espresso.org`.

Compilation
===========

`CP` is included in the core Quantum ESPRESSO distribution. Instruction
on how to install it can be found in the general documentation (User's
Guide) for Quantum ESPRESSO.

Typing `make cp` from the main Quantum ESPRESSO directory or `make` from
the `CPV/` subdirectory produces the following codes in `CPV/src`:

-   `cp.x`: Car-Parrinello Molecular Dynamics code

-   `cppp.x`: postprocessing code for `cp.x`. See `Doc/INPUT_CPPP.*` for
    input variables.

-   `wfdd.x`: utility code for finding maximally localized Wannier
    functions using damped dynamics.

Symlinks to executable programs will be placed in the `bin/`
subdirectory.

As a final check that compilation was successful, you may want to run
some or all of the tests and examples. Automated tests for `cp.x` are in
directory `test-suite/` and can be run via the `Makefile` found there.
Please see the general User's Guide for their setup.

You may take the tests and examples distributed with `CP` as templates
for writing your own input files. Input files for tests are contained in
subdirectories `test-suite/cp_*` with file type `*.in1`, `*.in2`, \... .
Input files for examples are produced, if you run the examples, in the
`results/` subdirectories, with names ending with `.in`.

For general information on parallelism and how to run in parallel
execution, please see the general User's Guide. `CP` currently can take
advantage of both MPI and OpenMP parallelization and on GPU acceleration.
The "plane-wave", "linear-algebra" and "task-group" parallelization levels
are implemented.

Input data
==========

Input data for `cp.x` is organized into several namelists, followed by
other fields ("cards") introduced by keywords. The namelists are:

>  &CONTROL:           general variables controlling the run\
>  &SYSTEM:            structural information on the system under investigation\
>  &ELECTRONS:         electronic variables, electron dynamics\
>  &IONS :             ionic variables, ionic dynamics\
>  &CELL (optional):   variable-cell dynamics\

The `&CELL` namelist may be omitted for fixed-cell calculations. This
depends on the value of variable `calculation` in namelist &CONTROL.
Most variables in namelists have default values. Only he following
variables in &SYSTEM must always be specified:

>  `ibrav`     (integer)             Bravais-lattice index\
>  `celldm`    (real, dimension 6)   crystallographic constants\
>  `nat`       (integer)             number of atoms in the unit cell\
>  `ntyp`      (integer)             number of types of atoms in the unit cell\
>  `ecutwfc`   (real)                kinetic energy cutoff (Ry) for wavefunctions

Explanations for the meaning of variables `ibrav` and `celldm`, as well
as on alternative ways to input structural data, are contained in files
`Doc/INPUT_CP.*`. These files are the reference for input data and
describe a large number of other variables as well. Almost all variables
have default values, which may or may not fit your needs.

After the namelists, you have several fields ("cards") introduced by
keywords with self-explanatory names:

> ATOMIC\_SPECIES\
> ATOMIC\_POSITIONS\
> CELL\_PARAMETERS (optional)\
> OCCUPATIONS (optional)

The keywords may be followed on the same line by an option. Unknown
fields are ignored. See the files mentioned above for details on the
available "cards".

Comment lines in namelists can be introduced by a \"!\", exactly as in
fortran code. Comments lines in "cards" can be introduced by either a "!"
or a "\#" character in the first position of a line.

Data files
----------

The output data files are written in the directory specified by variable
`outdir`, with names specified by variable `prefix` (a string that is
prepended to all file names, whose default value is `prefix=’cp_$ndw’`,
where `ndw` is an integer specified in input).
In order to use the data on a different machine, you may need to
compile `CP` with HDF5 enabled.

The execution stops if you create a file `prefix.EXIT` either in the
working directory (i.e. where the program is executed), or in the
`outdir` directory. Note that with some versions of MPI, the working
directory is the directory where the executable is! The advantage of
this procedure is that all files are properly closed, whereas just
killing the process may leave data and output files in an unusable
state.

The format of arrays containing charge density, potential, etc.
is described in the developer manual.

Output files
==========

The `cp.x` code produces many output files, that together build up the trajectory.

You have a file for the positions, called `prefix.pos`, where `prefix` is defined in
the input file, that is formatted like:

           10    0.00157227
           0.48652245874924E+01     0.38015905345591E+01     0.37361508020082E+01
           0.40077990926697E+01     0.59541011690914E+01     0.34691399577808E+01
           0.43874410242643E+01     0.38553718662714E+01     0.59039702898524E+01
           20    0.00641004
           0.49677092782926E+01     0.38629427979469E+01     0.37777995137803E+01
           0.42395189282719E+01     0.55766875434652E+01     0.31291744042209E+01
           0.45445534106843E+01     0.36049553522533E+01     0.55864387532281E+01

where the first line contains the step number and elapsed time, in ps, at this
step; the following lines contain the positions, in Bohr radii,  of all the
atoms (3 in this examples), in the same order as in the input file (since v6.6
-- previously, atoms were sorted by type; the type must be deduced from the
input file). The same structure is repeated for the second step and so on. 
The printout is made every `iprint` steps (10 in this case, so at step 10, 20,
etc.). Note that the atomic coordinates are not wrapped into the simulation
cell, so it is possible that they lie outside it.

The velocities are written in a similar file named `prefix.vel`, where `prefix`
is defined in the input file, that is formatted like the `.pos` file. The units
are the usual Hartree atomic units (note that the velocities in the `pw.x` code
are in _Rydberg_ a.u. and differ by a factor 2).

The `prefix.for` file, formatted like the previous two, contains the computed
forces, in Hartree atomic units as well. It is written only if a molecular
dynamics calculation is performed, or if `tprnfor = .true.` is set in input.

The simulation cell is written in a file named `prefix.cel` with the same header as the previous
described files, and the cell matrix is then listed. NB: **THE CELL MATRIX IN THE
OUTPUT IS TRANSPOSED** that means that if you want to reuse it again for a new input file,
you have to pick the one that you find in `prefix.cel` and write in the input file
after inverting rows and columns.

The file `prefix.evp` has one line per printed step and contains some
thermodynamical data.
The first line of the file names the columns:
```
#   nfi  time(ps)  ekinc  Tcell(K)  Tion(K)  etot  enthal  econs  econt  Volume  Pressure(GPa)
```
where:
   - `ekinc` is the electrons fictitious kinetic energy, $`K_{ELECTRONS}`$
   - `enthal` is the enthalpy, $`E_{DFT}+PV`$
   - `etot` is the DFT (potential) energy of the system, $`E_{DFT}`$
   - `econs` is a physically meaningful constant of motion, $`E_{DFT} + K_{NUCLEI}`$,
   in the limit of zero electronic fictitious mass
   - `econt` is the constant of motion of the lagrangian$`E_{DFT} + K_{IONS} + K_{ELECTRONS}`$ t.
   If the time step `dt` is small enough this will be up to a very good precision a constant.
   It is not a physical quantity, since $`K_{ELECTRONS}`$ has _nothing_ to do with the quantum
   kinetic energy of the electrons.


Using `CP`
==========

It is important to understand that a CP simulation is a sequence of
different runs, some of them used to \"prepare\" the initial state of
the system, and other performed to collect statistics, or to modify the
state of the system itself, i.e. to modify the temperature or the pressure.

To prepare and run a CP simulation you should first of all define the
system:

> atomic positions\
> system cell\
> pseudopotentials\
> cut-offs\
> number of electrons and bands (optional)\
> FFT grids (optional)

An example of input file (Benzene Molecule):

             &control
                title = 'Benzene Molecule',
                calculation = 'cp',
                restart_mode = 'from_scratch',
                ndr = 51,
                ndw = 51,
                nstep = 100,
                iprint = 10,
                isave = 100,
                tstress = .TRUE.,
                tprnfor = .TRUE.,
                dt    = 5.0d0,
                etot_conv_thr = 1.d-9,
                ekin_conv_thr = 1.d-4,
                prefix = 'c6h6',
                pseudo_dir='/scratch/benzene/',
                outdir='/scratch/benzene/Out/'
             /
             &system
                ibrav = 14,
                celldm(1) = 16.0,
                celldm(2) = 1.0,
                celldm(3) = 0.5,
                celldm(4) = 0.0,
                celldm(5) = 0.0,
                celldm(6) = 0.0,
                nat = 12,
                ntyp = 2,
                nbnd = 15,
                ecutwfc = 40.0,
                nr1b= 10, nr2b = 10, nr3b = 10,
                input_dft = 'BLYP'
             /
             &electrons
                emass = 400.d0,
                emass_cutoff = 2.5d0,
                electron_dynamics = 'sd'
             /
             &ions
                ion_dynamics = 'none'
             /
             &cell
                cell_dynamics = 'none',
                press = 0.0d0,
              /
              ATOMIC_SPECIES
              C 12.0d0 c_blyp_gia.pp
              H 1.00d0 h.ps
              ATOMIC_POSITIONS (bohr)
              C     2.6 0.0 0.0
              C     1.3 -1.3 0.0
              C    -1.3 -1.3 0.0
              C    -2.6 0.0 0.0
              C    -1.3 1.3 0.0
              C     1.3 1.3 0.0
              H     4.4 0.0 0.0
              H     2.2 -2.2 0.0
              H    -2.2 -2.2 0.0
              H    -4.4 0.0 0.0
              H    -2.2 2.2 0.0
              H     2.2 2.2 0.0

You can find the description of the input variables in file `Doc/INPUT_CP.*`.

Reaching the electronic ground state
------------------------------------

The first run, when starting from scratch, is always an electronic
minimization, with fixed ions and cell, to bring the electronic system
on the ground state (GS) relative to the starting atomic configuration.
This step is conceptually very similar to self-consistency in a
`pw.x` run.

Sometimes a single run is not enough to reach the GS. In this case, you
need to re-run the electronic minimization stage. Use the input of the
first run, changing `restart_mode = ’from_scratch’` to
`restart_mode = ’restart’`.

NOTA BENE: Unless you are already experienced with the system you are
studying or with the internals of the code, you will usually need to
tune some input parameters, like `emass`, `dt`, and cut-offs. For this
purpose, a few trial runs could be useful: you can perform short
minimizations (say, 10 steps) changing and adjusting these parameters to
fit your needs. You can specify the degree of convergence with these two
thresholds:

> `etot_conv_thr`: total energy difference between two consecutive
> steps\
> `ekin_conv_thr`: value of the fictitious kinetic energy of the
> electrons.

Usually we consider the system on the GS when `ekin_conv_thr`
$`< 10^{-5}`$. You could check the value of the fictitious kinetic energy
on the standard output (column EKINC).

Different strategies are available to minimize electrons, but the most
frequently used is _damped dynamics_: `electron_dynamics = ’damp’` and
`electron_damping` = a number typically ranging from 0.1 and 0.5.
See the input description to compute the optimal damping factor.
Steepest descent: `electron_dynamics = ’sd’`, is also available but it
is typicallyslower than damped dynamics and should be used only to
start the minimization.

Relax the system
----------------

Once your system is in the GS, depending on how you have prepared the
starting atomic configuration:

1.  if you have set the atomic positions \"by hand\" and/or from a
    classical code, check the forces on atoms, and if they are large
    ($`\sim 0.1 \div 1.0`$ atomic units), you should perform an ionic
    minimization, otherwise the system could break up during the
    dynamics.

2.  if you have taken the positions from a previous run or a previous
    ab-initio simulation, check the forces, and if they are too small
    ($`\sim 10^{-4}`$ atomic units), this means that atoms are already in
    equilibrium positions and, even if left free, they will not move.
    Then you need to randomize positions a little bit (see below).

Let us consider case 1). There are different strategies to relax the
system, but the most used are again steepest-descent or damped-dynamics
for ions and electrons. You could also mix electronic and ionic
minimization scheme freely, i.e. ions in steepest-descent and electron
in with damped-dynamics or vice versa.

-   suppose we want to perform steepest-descent for ions. Then we should
    specify the following section for ions:

                 &ions
                   ion_dynamics = 'sd'
                 /

    Change also the ionic masses to accelerate the minimization:

                 ATOMIC_SPECIES
                  C 2.0d0 c_blyp_gia.pp
                  H 2.00d0 h.ps

    while leaving other input parameters unchanged. *Note* that if the
    forces are really high ($`> 1.0`$ atomic units), you should always use
    steepest descent for the first ($`\sim 100`$ relaxation steps.

-   As the system approaches the equilibrium positions, the steepest
    descent scheme slows down, so is better to switch to damped
    dynamics:

                 &ions
                   ion_dynamics = 'damp',
                   ion_damping = 0.2,
                   ion_velocities = 'zero'
                 /

    A value of `ion_damping` around 0.05 is good for many systems. It is
    also better to specify to restart with zero ionic and electronic
    velocities, since we have changed the masses.

    Change further the ionic masses to accelerate the minimization:

                   ATOMIC_SPECIES
                   C 0.1d0 c_blyp_gia.pp
                   H 0.1d0 h.ps

-   when the system is really close to the equilibrium, the damped
    dynamics slow down too, especially because, since we are moving
    electron and ions together, the ionic forces are not properly
    correct, then it is often better to perform a ionic step every N
    electronic steps, or to move ions only when electron are in their GS
    (within the chosen threshold).

    This can be specified by adding, in the ionic section, the
    `ion_nstepe` parameter, then the &IONS namelist become as follows:

                 &ions
                   ion_dynamics = 'damp',
                   ion_damping = 0.2,
                   ion_velocities = 'zero',
                   ion_nstepe = 10
                 /

    Then we specify in the &CONTROL namelist:

                   etot_conv_thr = 1.d-6,
                   ekin_conv_thr = 1.d-5,
                   forc_conv_thr = 1.d-3

    As a result, the code checks every 10 electronic steps whether the
    electronic system satisfies the two thresholds `etot_conv_thr`,
    `ekin_conv_thr`: if it does, the ions are advanced by one step. The
    process thus continues until the forces become smaller than
    `forc_conv_thr`.

    *Note* that to fully relax the system you need many runs, and
    different strategies, that you should mix and change in order to
    speed-up the convergence. The process is not automatic, but is
    strongly based on experience, and trial and error.

    Remember also that the convergence to the equilibrium positions
    depends on the energy threshold for the electronic GS, in fact
    correct forces (required to move ions toward the minimum) are
    obtained only when electrons are in their GS. Then a small threshold
    on forces could not be satisfied, if you do not require an even
    smaller threshold on total energy.

Let us now move to case 2: randomization of positions.

If you have relaxed the system or if the starting system is already in
the equilibrium positions, then you need to displace ions from the
equilibrium positions, otherwise they will not move in a dynamics
simulation. After the randomization you should bring electrons on the GS
again, in order to start a dynamic with the correct forces and with
electrons in the GS. Then you should switch off the ionic dynamics and
activate the randomization for each species, specifying the amplitude of
the randomization itself. This could be done with the following &IONS
namelist:

              &ions
                ion_dynamics = 'none',
                tranp(1) = .TRUE.,
                tranp(2) = .TRUE.,
                amprp(1) = 0.01
                amprp(2) = 0.01
              /

In this way a random displacement (of max 0.01 a.u.) is added to atoms
of species 1 and 2. All other input parameters could remain the same.
Note that the difference in the total energy (etot) between relaxed and
randomized positions can be used to estimate the temperature that will
be reached by the system. In fact, starting with zero ionic velocities,
all the difference is potential energy, but in a dynamics simulation,
the energy will be equipartitioned between kinetic and potential, then
to estimate the temperature take the difference in energy (de), convert
it in Kelvin, divide for the number of atoms and multiply by 2/3.
Randomization could be useful also while we are relaxing the system,
especially when we suspect that the ions are in a local minimum or in an
energy plateau.

CP dynamics
-----------

At this point after having minimized the electrons, and with ions
displaced from their equilibrium positions, we are ready to start a CP
dynamics. We need to specify `’verlet’` both in ionic and electronic
dynamics. The threshold in control input section will be ignored, like
any parameter related to minimization strategy. The first time we
perform a CP run after a minimization, it is always better to put
velocities equal to zero, unless we have velocities, from a previous
simulation, to specify in the input file. Restore the proper masses for
the ions. In this way we will sample the microcanonical ensemble. The
input section changes as follow:

               &electrons
                  emass = 400.d0,
                  emass_cutoff = 2.5d0,
                  electron_dynamics = 'verlet',
                  electron_velocities = 'zero'
               /
               &ions
                  ion_dynamics = 'verlet',
                  ion_velocities = 'zero'
               /
               ATOMIC_SPECIES
               C 12.0d0 c_blyp_gia.pp
               H 1.00d0 h.ps

If you want to specify the initial velocities for ions, you have to set
`ion_velocities =’from_input’`, and add the ATOMIC\_VELOCITIES card,
after the ATOMIC\_POSITION card, with the list of velocities in atomic
units.

NOTA BENE: in restarting the dynamics after the first CP run, remember
to remove or comment the velocities parameters:

               &electrons
                  emass = 400.d0,
                  emass_cutoff = 2.5d0,
                  electron_dynamics = 'verlet'
                  ! electron_velocities = 'zero'
               /
               &ions
                  ion_dynamics = 'verlet'
                  ! ion_velocities = 'zero'
               /

otherwise you will quench the system interrupting the sampling of the
microcanonical ensemble.

####  Varying the temperature 

It is possible to change the temperature of the system or to sample the
canonical ensemble fixing the average temperature, this is done using
the Nosé thermostat. To activate this thermostat for ions you have to
specify in namelist &IONS:

               &ions
                  ion_dynamics = 'verlet',
                  ion_temperature = 'nose',
                  fnosep = 60.0,
                  tempw = 300.0
               /  

where `fnosep` is the frequency of the thermostat in THz, that should be
chosen to be comparable with the center of the vibrational spectrum of
the system, in order to excite as many vibrational modes as possible.
`tempw` is the desired average temperature in Kelvin.

*Note:* to avoid a strong coupling between the Nosé thermostat and the
system, proceed step by step. Don't switch on the thermostat from a
completely relaxed configuration: adding a random displacement is
strongly recommended. Check which is the average temperature via a few
steps of a microcanonical simulation. Don't increase the temperature too
much. Finally switch on the thermostat. In the case of molecular system,
different modes have to be thermalized: it is better to use a chain of
thermostat or equivalently running different simulations with different
frequencies.

####  Nośe thermostat for electrons 

It is possible to specify also the thermostat for the electrons. This is
usually activated in metals or in systems where we have a transfer of
energy between ionic and electronic degrees of freedom. Beware: the
usage of electronic thermostats is quite delicate. The following
information comes from K. Kudin:

"The main issue is that there is usually some \"natural\" fictitious
kinetic energy that electrons gain from the ionic motion (\"drag\"). One
could easily quantify how much of the fictitious energy comes from this
drag by doing a CP run, then a couple of CG (same as BO) steps, and then
going back to CP. The fictitious electronic energy at the last CP
restart will be purely due to the drag effect."

"The thermostat on electrons will either try to overexcite the otherwise
\"cold\" electrons, or it will try to take them down to an unnaturally
cold state where their fictitious kinetic energy is even below what
would be just due pure drag. Neither of this is good."

"I think the only workable regime with an electronic thermostat is a
mild overexcitation of the electrons, however, to do this one will need
to know rather precisely what is the fictitious kinetic energy due to
the drag."

Advanced usage
--------------

### Autopilot features

For changing variables while the simulation is running see
[the autopilot guide](autopilot_guide.md)

###  Self-interaction Correction 

The self-interaction correction (SIC) included in the `CP` package is
based on the Constrained Local-Spin-Density approach proposed my F.
Mauri and coworkers (M. D'Avezac et al. PRB 71, 205210 (2005)). It was
used for the first time in Quantum ESPRESSO by F. Baletto, C. Cavazzoni
and S.Scandolo (PRL 95, 176801 (2005)).

This approach is a simple and nice way to treat ONE, and only one,
excess charge. It is moreover necessary to check a priori that the
spin-up and spin-down eigenvalues are not too different, for the
corresponding neutral system, working in the Local-Spin-Density
Approximation (setting `nspin = 2`). If these two conditions are
satisfied and you are interest in charged systems, you can apply the
SIC. This approach is a on-the-fly method to correct the
self-interaction with the excess charge with itself.

Briefly, both the Hartree and the XC part have been corrected to avoid
the interaction of the excess charge with itself.

For example, for the Boron atoms, where we have an even number of
electrons (valence electrons = 3), the parameters for working with the
SIC are:

               &system
               nbnd= 2,
               tot_magnetization=1,
               sic_alpha = 1.d0,
               sic_epsilon = 1.0d0,
               sic = 'sic_mac',
               force_pairing = .true.,

The two main parameters are:

> `force_pairing = .true.`, which forces the paired electrons to be the
> same;\
> `sic=’sic_mac’`, which instructs the code to use Mauri's correction.

**Warning**: This approach has known problems for dissociation mechanism
driven by excess electrons.

Comment 1: Two parameters, `sic_alpha` and `sic_epsilon’`, have been
introduced following the suggestion of M. Sprik (ICR(05)) to treat the
radical (OH)-H$`_2`$O. In any case, a complete ab-initio approach is
followed using `sic_alpha=1`, `sic_epsilon=1`.

Comment 2: When you apply this SIC scheme to a molecule or to an atom,
which are neutral, remember to add the correction to the energy level as
proposed by Landau: in a neutral system, subtracting the
self-interaction, the unpaired electron feels a charged system, even if
using a compensating positive background. For a cubic box, the
correction term due to the Madelung energy is approx. given by
$`1.4186/L_{box} - 1.047/(L_{box})^3`$, where $`L_{box}`$ is the linear
dimension of your box (=celldm(1)). The Madelung coefficient is taken
from I. Dabo et al. PRB 77, 115139 (2007). (info by F. Baletto,
francesca.baletto\@kcl.ac.uk)

###  ensemble-DFT 

The ensemble-DFT (eDFT) is a robust method to simulate the metals in the
framework of "ab-initio" molecular dynamics. It was introduced in 1997
by Marzari et al.

The specific subroutines for the eDFT are in `CPV/src/ensemble_dft.f90`
where you define all the quantities of interest. The subroutine
`CPV/src/inner_loop_cold.f90` called by `cg_sub.f90`, control the inner
loop, and so the minimization of the free energy $`A`$ with respect to the
occupation matrix.

To select a eDFT calculations, the user has to set:

                calculation = 'cp'
                occupations= 'ensemble' 
                tcg = .true.
                passop= 0.3
                maxiter = 250

to use the CG procedure. In the eDFT it is also the outer loop, where
the energy is minimized with respect to the wavefunction keeping fixed
the occupation matrix. While the specific parameters for the inner loop.
Since eDFT was born to treat metals, keep in mind that we want to
describe the broadening of the occupations around the Fermi energy.
Below the new parameters in the electrons list, are listed.

-   `smearing`: used to select the occupation distribution; there are
    two options: Fermi-Dirac smearing='fd', cold-smearing smearing='cs'
    (recommended)

-   `degauss`: is the electronic temperature; it controls the broadening
    of the occupation numbers around the Fermi energy.

-   `ninner`: is the number of iterative cycles in the inner loop, done
    to minimize the free energy $`A`$ with respect the occupation numbers.
    The typical range is 2-8.

-   `conv_thr`: is the threshold value to stop the search of the
    'minimum' free energy.

-   `niter_cold_restart`: controls the frequency at which a full
    iterative inner cycle is done. It is in the range $`1\div`$ `ninner`.
    It is a trick to speed up the calculation.

-   `lambda_cold`: is the length step along the search line for the best
    value for $`A`$, when the iterative cycle is not performed. The value
    is close to 0.03, smaller for large and complicated metallic
    systems.

*NOTE:* `degauss` is in Hartree, while in `PWscf`is in Ry (!!!). The
typical range is 0.01-0.02 Ha.

The input for an Al surface is:

                &CONTROL
                 calculation = 'cp',
                 restart_mode = 'from_scratch',
                 nstep  = 10,
                 iprint = 5,
                 isave  = 5,
                 dt    = 125.0d0,
                 prefix = 'Aluminum_surface',
                 pseudo_dir = '~/UPF/',
                 outdir = '/scratch/'
                 ndr=50
                 ndw=51
                /
                &SYSTEM
                 ibrav=  14,
                 celldm(1)= 21.694d0, celldm(2)= 1.00D0, celldm(3)= 2.121D0,
                 celldm(4)= 0.0d0,   celldm(5)= 0.0d0, celldm(6)= 0.0d0,
                 nat= 96,
                 ntyp= 1,
                 nspin=1,
                 ecutwfc= 15,
                 nbnd=160,
                 input_dft = 'pbe'
                 occupations= 'ensemble',
                 smearing='cs',
                 degauss=0.018,
                /
                &ELECTRONS
                 orthogonalization = 'Gram-Schmidt',
                 startingwfc = 'random',
                 ampre = 0.02,
                 tcg = .true.,
                 passop= 0.3,
                 maxiter = 250,
                 emass_cutoff = 3.00,
                 conv_thr=1.d-6
                 n_inner = 2,
                 lambda_cold = 0.03,
                 niter_cold_restart = 2,
                /
                &IONS
                 ion_dynamics  = 'verlet',
                 ion_temperature = 'nose'
                 fnosep = 4.0d0,
                 tempw = 500.d0
                /
                ATOMIC_SPECIES
                 Al 26.89 Al.pbe.UPF

*NOTA1* remember that the time step is to integrate the ionic dynamics,
so you can choose something in the range of 1-5 fs.\
*NOTA2* with eDFT you are simulating metals or systems for which the
occupation number is also fractional, so the number of band, `nbnd`, has
to be chosen such as to have some empty states. As a rule of thumb,
start with an initial occupation number of about 1.6-1.8 (the more bands
you consider, the more the calculation is accurate, but it also takes
longer. The CPU time scales almost linearly with the number of bands.)\
*NOTA3* the parameter `emass_cutoff` is used in the preconditioning and
it has a completely different meaning with respect to plain CP. It
ranges between 4 and 7.

All the other parameters have the same meaning in the usual `CP` input,
and they are discussed above.

### Treatment of USPPs

The cutoff `ecutrho` defines the resolution on the real space FFT mesh
(as expressed by `nr1`, `nr2` and `nr3`, that the code left on its own
sets automatically). In the USPP case we refer to this mesh as the
\"hard\" mesh, since it is denser than the smooth mesh that is needed to
represent the square of the non-norm-conserving wavefunctions.

On this \"hard\", fine-spaced mesh, you need to determine the size of
the cube that will encompass the largest of the augmentation charges -
this is what `nr1b`, `nr2b`, `nr3b` are. hey are independent of the
system size, but dependent on the size of the augmentation charge (an
atomic property that doesn't vary that much for different systems) and
on the real-space resolution needed by augmentation charges (rule of
thumb: `ecutrho` is between 6 and 12 times `ecutwfc`).

The small boxes should be set as small as possible, but large enough to
contain the core of the largest element in your system. The formula for
estimating the box size is quite simple:

> `nr1b` = $`2 R_c / L_x \times`$ `nr1`

and the like, where $`R_{cut}`$ is largest cut-off radius among the
various atom types present in the system, $`L_x`$ is the physical length
of your box along the $`x`$ axis. You have to round your result to the
nearest larger integer. In practice, `nr1b` etc. are often in the region
of 20-24-28; testing seems again a necessity.

The core charge is in principle finite only at the core region (as
defined by some $`R_{rcut}`$ ) and vanishes out side the core. Numerically
the charge is represented in a Fourier series which may give rise to
small charge oscillations outside the core and even to negative charge
density, but only if the cut-off is too low. Having these small boxes
removes the charge oscillations problem (at least outside the box) and
also offers some numerical advantages in going to higher cut-offs.\"
(info by Nicola Marzari)

### Hybrid functional calculations using maximally localized Wannier functions

In this section, we illustrate some guidelines to perform exact exchange
(EXX) calculations using Wannier functions efficiently.

The references for this algorithm are:

-   Theory: X. Wu , A. Selloni, and R. Car, Phys. Rev. B 79, 085102
    (2009).

-   Implementation: H.-Y. Ko, B. Santra, R. A. DiStasio, L. Kong, Z.
    Li, X. Wu, and R. Car, arxiv.

The parallelization scheme in this algorithm is based upon the number of
electronic states. In the current implementation, there are certain
restrictions on the choice of the number of MPI tasks. Also slightly
different algorithms are employed depending on whether the number of MPI
tasks used in the calculation are greater or less than the number of
electronic states. We highly recommend users to follow the notes below.
This algorithm can be used most efficiently if the numbers of electronic
states are uniformly distributed over the number of MPI tasks. For a
system having N electronic states the optimum numbers of MPI tasks
(nproc) are the following:

-   In case of nproc $`\leq`$ N, the optimum choices are N/m, where m is
    any positive integer.

    -   Robustness: Can be used for odd and even number of electronic
        states.

    -   OpenMP threads: Can be used.

    -   Taskgroup: Only the default value of the task group (-ntg 1) is
        allowed.

-   In case of nproc $`>`$ N, the optimum choices are N\*m, where m is any
    positive integer.

    -   Robustness: Can be used for even number of electronic states.

    -   Largest value of m: As long as nj\_max (see output) is greater
        than 1, however beyond m=8 the scaling may become poor. The
        scaling should be tested by users.

    -   OpenMP threads: Can be used and highly recommended. We have
        tested number of threads starting from 2 up to 64. More threads
        are also allowed. For very large calculations (nproc $`>`$ 1000 )
        efficiency can largely depend on the computer architecture and
        the balance between the MPI tasks and the OpenMP threads. User
        should test for an optimal balance. Reasonably good scaling can
        be achieved by using m=6-8 and OpenMP threads=2-16.

    -   Taskgroup: Can be greater than 1 and users should choose the
        largest possible value for ntg. To estimate ntg, find the value
        of nr3x in the output and compute nproc/nr3x and take the
        integer value. We have tested the value of ntg as $`2^m`$, where m
        is any positive integer. Other values of ntg should be used with
        caution.

    -   Ndiag: Use -ndiag X option in the execution of cp.x. Without
        this option jobs may crash on certain architectures. Set X to
        any perfect square number which is equal to or less than N.

-   DEBUG: The EXX calculations always work when number of MPI tasks =
    number of electronic states. In case of any uncertainty, the EXX
    energy computed using different numbers of MPI tasks can be checked
    by performing test calculations using number of MPI tasks = number
    of electronic states.

An example input is listed as following:

    &CONTROL
      calculation       = 'cp-wf',
      title             = "(H2O)32 Molecule: electron minimization PBE0",
      restart_mode      = "from_scratch",
      pseudo_dir        = './',
      outdir            = './',
      prefix            = "water",
      nstep             = 220,
      iprint            = 100,
      isave             = 100,
      dt                = 4.D0,
      ekin_conv_thr     = 1.D-5,
      etot_conv_thr     = 1.D-5,
    /
    &SYSTEM
      ibrav             = 1,
      celldm(1)         = 18.6655, 
      nat               = 96,
      ntyp              = 2,
      ecutwfc           = 85.D0,
      input_dft         = 'pbe0',
    /
    &ELECTRONS
      emass             = 400.D0,
      emass_cutoff      = 3.D0,
      ortho_eps         = 1.D-8,
      ortho_max         = 300,
      electron_dynamics = "damp",
      electron_damping  = 0.1D0,
    /
    &IONS
      ion_dynamics      = "none", 
    /
    &WANNIER
      nit               = 60,
      calwf             = 3,
      tolw              = 1.D-6,
      nsteps            = 20,
      adapt             = .FALSE.
      wfdt              = 4.D0,
      wf_q              = 500,
      wf_friction       = 0.3D0,
      exx_neigh         = 60,     ! exx related optional
      exx_dis_cutoff    = 8.0D0,  ! exx related optional
      exx_ps_rcut_self  = 6.0D0,  ! exx related optional
      exx_ps_rcut_pair  = 5.0D0,  ! exx related optional
      exx_me_rcut_self  = 9.3D0,  ! exx related optional
      exx_me_rcut_pair  = 7.0D0,  ! exx related optional
      exx_poisson_eps   = 1.D-6,  ! exx related optional
    /
    ATOMIC_SPECIES
    O 16.0D0 O_HSCV_PBE-1.0.UPF
    H  2.0D0 H_HSCV_PBE-1.0.UPF 

Parallel Performances
=====================

`cp.x` can run in principle on any number of processors. The
effectiveness of parallelization is ultimately judged by the "scaling",
i.e. how the time needed to perform a job scales with the number of
processors. Ideally one would like to have linear scaling, i.e.
$`T \sim T_0/N_p`$ for $`N_p`$ processors, where $`T_0`$ is the estimated
time for serial execution. In addition, one would like to have linear
scaling of the RAM per processor: $`O_N \sim O_0/N_p`$, so that large-memory
systems fit into the RAM of each processor.

We refer to the "Parallelization" section of the general User's Guide for
a description of MPI and OpenMP parallelization paradigms, of the various
MPI parallelization levels, and on how to activate them. 

A judicious choice of the various levels of parallelization, together
with the availability of suitable hardware (e.g. fast communications)
is fundamental to reach good performances._VERY IMPORTANT_: For each
system there is an optimal range of number of processors on which to
run the job. A too large number of processors or a bad parallelization
style will yield performance degradation. 

For `CP` with hybrid functionals, see the related section above this one. 
For all other cases, the relevant MPI parallelization levels are:

- "plane waves" (PW);
- "tasks" (activated by command-line option `-nt N`);
- "linear algebra" (`-nd N`);
- "bands" parallelization (`-nb N`), to be used only in
special cases;
- "images" parallelization (`-ni N`), used only in code `manycp.x`
(see the header of `CPV/src/manycp.f90` for documentation).

As a rule of thumb:
- start with PW parallelization only (e.g. `mpirun -np N cp.x ...` with
no other parallelization options); the code will scale well unless `N` 
exceeds the third FFT dimensions `nr3` and/or `nr3s`.
- To further increase the number of processors, use "task groups",
typically 4 to 8 (e.g. `mpirun -np N cp.x -nt 8 ...`).
- Alternatively, or in addition, you may compile with OpenMP:
`./configure --enable-openmp ...`, then `export OMP_NUM_THREADS=n`
and run on `n` threads (4 to 8 typically).
_Beware conflicts between MPI and OpenMP threads_!
don't do this unless you know what you are doing.
- Finally, the optimal number of processors for \"linear-algebra\"
parallelization can be found by observing the performances of `ortho` 
in the final time report for different numbers of processors in the 
linear-algebra group (must be a square integer, not larger than the 
number of processoris for plane-wave parallelization). Linear-algebra
parallelization distributes `M\times M`$ matrices, with `M` number of
bands, so it may be useful if memory-constrained.

Note: optimal serial performances are achieved when the data are as much
as possible kept into the cache. As a side effect, PW parallelization may
yield superlinear (better than linear) scaling, thanks to the increase in
serial speed coming from the reduction of data size (making it easier for
the machine to keep data in the cache).

UtilXlib
========

This library implements various basic tasks such as timing, tracing,
optimized memory accesses and an abstraction layer for the MPI subroutines.

The following pre-processor directives can be used to enable/disable some
features:

* `__MPI` : activates MPI support.
* `__TRACE` : activates verbose output for debugging purposes
* `__CUDA` : activates CUDA Fortran based interfaces.
* `__GPU_MPI` : use CUDA aware MPI calls instead of standard sync-send-update method (experimental).


Usage of wrapper interfaces for MPI
===================================

This library offers a number of interfaces to abstract the MPI APIs and 
to optionally relax the dependency on a MPI library.

`mp_*` interfaces present in the library can only be called after the 
initialization performed by the subroutine `mp_start` and before the
finalization done by `mp_end`.
All rules have exceptions and indeed subroutines `mp_count_nodes`,
`mp_type_create_column_section` and `mp_type_free` can also be called 
outside the aforementioned window.

If CUDA Fortran support is enabled, almost all interfaces accept input
data declared with the `device` attribute. Note however that CUDA Fortran
support should be considered experimental.


CUDA specific notes
===================

All calls to message passing interfaces are synchronous with respect to
both MPI and CUDA streams. The code will synchronize the device before
starting the communication, also in those cases where communication
may be avoided (for example in serial version).
A different behaviour may be observed when the default stream 
synchronization behaviour is overridden by the user (see `cudaStreamCreateWithFlags`).

Be careful when using CUDA-aware MPI. Some implementations are not
complete. The library will not check for the CUDA-aware MPI APIs during
the initialization, but may report failure codes during the execution.
If you encounter problems when adding the flag `__GPU_MPI` it might
be that the MPI library does not support some CUDA-aware APIs.


Testing
=======

Partial unit testing is available in the `tests` sub-directory. See the 
README.md file in that directory for further information.
# UtilXlib testing suite

A set of tests for UtilXlib also showing the functionalities of the library.

In order to run the tests first run `./configure` in QE topdir and generate a
valid make.inc
After that edit the variables at the top of `compile_and_run_tests.sh` and run with:

    bash compile_and_run_tests.sh -sm[cn]

Options meaning:

 * -s: serial compilation and execution
 * -m: MPI compilation and execution
 * -c: CUDA-Fortran interface with data transfer done with memory on the *host*
 * -n: CUDA-Fortran interface with data transfer done with memory on the *device*

The tester module is part of the fortran_tester project avaiulable at:
https://github.com/pdebuyl/fortran_tester/archive/master.zip

TODO:

 * Get rid of seeds printed on screen. It's pointless.
title: Quantum Espresso PHonon
src_dir: ./PH/
output_dir: ./doc
project_website: http://www.quantum-espresso.org/
summary: Quantum Espresso is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials.
authors: Paolo Giannozzi et al. 
author_description: The Quantum Espresso project is developed by the QE community.
github: https://github.com/QEF/q-e
gitlab: https://gitlab.com/QEF/q-e
email: users@lists.quantum-espresso.org
predocmark: >
media_dir: ./media
page_dir: ./Ford
docmark_alt: #
predocmark_alt: <
display: public
         private
source: false
graph: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: GNU
extra_filetypes: sh #

[Quantum ESPRESSO](http://www.quantum-espresso.org/) is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale.
It is based on density-functional theory, plane waves, and pseudopotentials.

Quantum ESPRESSO has evolved into a distribution of inter-operable codes in the spirit of an open-source project. The Quantum ESPRESSO distribution consists of a core set of components, a set of additional packages performing more advanced tasks, and a number of third-party packages designed to be inter-operable with the core components. Researchers active in the field of electronic-structure calculations are encouraged to participate in the project by contributing their own codes or by implementing their own ideas into existing codes.

Quantum ESPRESSO is an open initiative, in collaboration with many groups world-wide, coordinated by the Quantum ESPRESSO Foundation. Present members of the latter include Scuola Internazionale Superiore di Studi Avanzati, the Abdus Salam International Centre for Theoretical Physics (Trieste), the CINECA National Supercomputing Center (Bologna), the Ecole Polytechnique Fédérale de Lausanne, the University of North Texas (Dallas), the Oxford University. Courses on modern electronic-structure theory with hands-on tutorials on the Quantum ESPRESSO codes are offered on a regular basis in collaboration with the Abdus Salam International Centre for Theoretical Physics in Trieste.

The reference technical manuscript for the latest public release is:
[P. Giannozzi et al., J.Phys.:Condens.Matter 29, 465901 (2017)](http://iopscience.iop.org/article/10.1088/1361-648X/aa8f79)


title: Quantum ESPRESSO overview
author: Samuel Poncé
date: 10-01-2017


<div style="text-align:center"><img src ="https://www.quantum-espresso.org/project/logos/Quantum_espresso_logo.jpg" width="600"></div>


[[Quantum Espresso]](http://www.quantum-espresso.org/) is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials.

QE is licensed under a [[GNU General Public License]](http://www.gnu.org/licenses/gpl-3.0.en.html)


## Installation 

Quantum ESPRESSO is currently distributed as source packages, but selected binary packages for Linux and Windows are also available. Since v.6.0 the distribution is split into sources, tests, examples, plus a few specialized packages and plug-ins. These packages can be automatically downloaded and installed on demand after the installation of the base distribution: see the general documentation for more information.


title: Quantum Espresso PWSCF
src_dir: ./../
output_dir: ./doc
project_website: http://www.quantum-espresso.org/
summary: Quantum Espresso is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials.
authors: Paolo Giannozzi et al. 
author_description: The Quantum Espresso project is developed by the QE community.
github: https://github.com/QEF/q-e
gitlab: https://gitlab.com/QEF/q-e
email: users@lists.quantum-espresso.org
predocmark: >
media_dir: ./media
page_dir: ./pagedir
docmark_alt: #
predocmark_alt: <
display: public
         private
source: false
graph: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: GNU
extra_filetypes: sh #

[Quantum ESPRESSO](http://www.quantum-espresso.org/) is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale.
It is based on density-functional theory, plane waves, and pseudopotentials.

Quantum ESPRESSO has evolved into a distribution of inter-operable codes in the spirit of an open-source project. The Quantum ESPRESSO distribution consists of a core set of components, a set of additional packages performing more advanced tasks, and a number of third-party packages designed to be inter-operable with the core components. Researchers active in the field of electronic-structure calculations are encouraged to participate in the project by contributing their own codes or by implementing their own ideas into existing codes.

Quantum ESPRESSO is an open initiative, in collaboration with many groups world-wide, coordinated by the Quantum ESPRESSO Foundation. Present members of the latter include Scuola Internazionale Superiore di Studi Avanzati, the Abdus Salam International Centre for Theoretical Physics (Trieste), the CINECA National Supercomputing Center (Bologna), the Ecole Polytechnique Fédérale de Lausanne, the University of North Texas (Dallas), the Oxford University. Courses on modern electronic-structure theory with hands-on tutorials on the Quantum ESPRESSO codes are offered on a regular basis in collaboration with the Abdus Salam International Centre for Theoretical Physics in Trieste.

The reference technical manuscript for the latest public release is:
[P. Giannozzi et al., J.Phys.:Condens.Matter 29, 465901 (2017)](http://iopscience.iop.org/article/10.1088/1361-648X/aa8f79)


title: Quantum ESPRESSO overview
author: Samuel Poncé
date: 10-01-2017


<div style="text-align:center"><img src ="https://www.quantum-espresso.org/project/logos/Quantum_espresso_logo.jpg" width="600"></div>


[[Quantum Espresso]](http://www.quantum-espresso.org/) is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials.

QE is licensed under a [[GNU General Public License]](http://www.gnu.org/licenses/gpl-3.0.en.html)


## Installation 

Quantum ESPRESSO is currently distributed as source packages, but selected binary packages for Linux and Windows are also available. Since v.6.0 the distribution is split into sources, tests, examples, plus a few specialized packages and plug-ins. These packages can be automatically downloaded and installed on demand after the installation of the base distribution: see the general documentation for more information.


# LAXlib testing suite

In order to run the tests first run `./configure` in QE topdir and generate a
valid `make.inc`.
You may also download large eigenvalue problems (in this directory) with:

    wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1EAB3xkoD-i9p4nW6NJDED3WaEK8ZCcf4' -O SiGeK1.bin
    wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=13lFkDbv99V8fqiXER1N2IzoJ_EhuGtt9' -O SiGeK2.bin

Finally `make` and run all `.x` executable files.
# Notes in the implementation of the KLI aproximation

# INSTALLATION OVERVIEW

 Before running the examples, you need to be sure that QEHeat is installed.
 To compile, we suggest the following procedure. As usual in a QuantumEspresso installation, enter the distribution folder and run autoconf :

```
 > ./configure
```

 Than from the distribution folder:

```
 >  make all_currents
```

 This will produce the executable `all_currents.x`, the executable for QEHeat, in the respective src and bin folders.


# EXAMPLES

 See also `../Doc/INPUT_ALL_CURRENTS.html` for a description of the inputs.
 To run an example, enter in the respective trajectory and execute the script `run_example.sh`. Modify it if necessary.
 In general, running a Qeheat calculation just needs the execution of the command: 
  `all_currents.x -in input_energycurrent` , 
 after the input file has been prepared.

 Each example comes with a reference folder where the output files can be compared with the ones produced by a new installation/run.
 Pseudopotentials can be downloaded from [http://www.quantum-simulation.org/potentials/sg15_oncv/]()
 For the examples, the following pseudos were used:
 H_HSCV_PBE-1.0.upf,  O_HSCV_PBE-1.0.upf,  O_ONCV_PBE-1.0.upf,  Si_ONCV_PBE-1.1.upf 
 , which should be present in the folder `pseudo`

 Example 1 and 2 need a parallel installation to finish in a reasonable time. Example 1 was run in the reference calculation on 4 cores and Example 2 on 12. 

 Example 3 can be easily run on a single core (serial) installation. Note that Example 3 requires the program `cp.x` to be installed, for example with `make cp`
 from the main distribution folder.



# 1. example_SiO2_single  



Here we compute the energy current and its indivial components for a single snapshot of Silica. 
For this purpose, one needs the additional namelist in the input file :

```
&energy_current
    delta_t=   0.500,
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
 /

& IONS
    ion_velocities = 'from_input',
```

and the CARD ATOMIC_VELOCITIES in the input file must be filled as well with the istantenous atomic velocities.
These are the ingredients needed for a basic single snapshot calculation.
 
The files produced, apart from the standard output, are :

- `file_output.dat` : this reports the total energy current, the eletronic current and the center of mass velocity for each species.
Only file_output.dat needs to be used to evaluate the thermal conductivity coefficient.

- `file_output` : this reports the total energy current divided in individual components, if a more specific analysis is needed. This file is mainly thought for development purposes.

`file_output.dat` comes with a header specifying the output units. The same units are used in the more detailed current decomposion given in `file_output`.
 See also the description ../Doc/INPUT_ALL_CURRENTS.html 



# 2. example_H2O_trajectory



Here we evaluate the energy current from a previously computed Car-Parrinello (CP) trajectory, which is provided together with the input file. The trajectory provided 
comes from a 125 water molecule simulation.

We calculate the energy current for every timestep of the trajectory located in  `${trajdir}.pos` and `${trajdir}.vel` (velocities are in CP units in this example). 
For this purpose we need to insert some additional keywords in the energy_current namelists :

```
 &energy_current
    delta_t=   0.500,
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
    trajdir='traj'
    first_step=1,
    vel_input_units='CP'
 /
```

note that the only different keywords with respect to a single snapshot calculations are `trajdir` and `first_step` in the `energy_current` namelist. Still, in the IONS namelist 
the keyword ion_velocities='from_input' must be set and the ATOMIC_VELOCITIES card must be filled.

The output with the istantenous energy currents is written in the files with names `${file_output}` and `${file_output}.dat`, as in example 1. 
The same format of the single snapshot calculation is kept, data from all the snapshots of the trajectory being appended sequentially.

If computational time is an issue, for this example we suggest to insert the flag last_step=xxx in the energy_current namelist,
replacing xxx with the desired index step, and than execute the run script.
Note that the keywords  `last_step` and  `first_step` refer to the indexes reported in the files  `${trajdir}.pos` and `${trajdir}.vel` and are not sequential indexes. The snapshot
in the input file is assigned an index 0. As a concrete example, the combination of `first_step=1`  and `last_step=953008` will skip the snapshot of the input file and
evaluate only the first snapshot of the trajectory because 953008 is the first index that appears in the trajectory file.



# 3. example_small_H2O_trajectory



This example is very similar to the previous one, but a Car-Parrinello trajectory is computed on-the-fly via the cp.x program of the just installed QE distribution. It produces the  
trajectory of a single water molecule and therefore the calculation is suited for a serial environment. 

This example requires the program cp.x to be installed. If this is not the case, you can enter the distribution folder and run "make cp".

Note that the trajectory produced by cp.x will be probably different due to the stochasticity inherent in the Car-Parrinello molecular dynamics simulation. For exact comparison 
with a novel installation one can substitute `trajdir='reference/traj/cp'` and comment in the run_example_water script the call to cp.x. This way the files produced by `file_output`
should be comparable with the reference, up to numerical noise that is always present in the finite difference derivative with non perfectly converged wavefunctions.




testcode
========

testcode is a python module for testing for regression errors in numerical
(principally scientific) software.  Essentially testcode runs a set of
calculations, and compares the output data to that generated by a previous
calculation (which is regarded to be "correct").  It is designed to be
lightweight and highly portable.

testcode can run a set of tests and check the calculated data is within a the
desired tolerance of results contained in previous output (using an internal
data extraction engine, a user-supplied data extraction program or
a user-supplied verification program).  The programs to be tested can be run in
serial and in parallel and tests can be run in either locally or submitted to
a compute cluster running a queueing system such as PBS.  Previous tests can be
compared and diffed against other tests or benchmarks.

Documentation
-------------

Full documentation can be found in the ``docs/`` subdirectory and in the
appropriate docstrings.  Documentation can be compiled using `sphinx
<http://sphinx.pocoo.org/>`_.

Documentation can also be viewed at `readthedocs
<http://testcode.readthedocs.org>`_.

Author
------

James Spencer, Imperial College London.

Contributions and suggestions from:

Keith Refson, Science and Technology Facilities Council.

Shawn Chin, Science and Technology Facilities Council.

LICENSE
-------

Modified BSD license; see LICENSE for more details.

See also
--------

`testcode_buildbot.py <https://gist.github.com/shawnchin/5678957>`_: a custom buildbot BuildStep for running testcode by Shawn Chin.
.. _verification:

Test verification
=================

testcode compares selected data from an output with previously obtained output
(the 'benchmark'); a test passes if all data is within a desired tolerance.
The data can be compared using an absolute tolerance and/or a relative
tolerance.  testcode needs some way of knowing what data from the output files
should be validated.  There are three options.

* label output with a 'data tag'

  If a data tag is supplied, then testcode will search each output file for
  lines starting with that tag.  The first numerical entry on those lines will
  then be checked against the benchmark.  For example, if the data tag is set
  to be '[QA]', and the line

      [QA] Energy = 1.23456 eV

  appears in the test output, then testcode will ensure the value 1.23456 is
  identical (within the specified tolerance) to the equivalent line in the
  benchmark output.  The text preceding the value is used to label that data
  item; lines with identical text but different values are handled but it is
  assumed that such lines always come in the same (relative) order.

* user-supplied data extraction program

  An external program can be used to extract data from the test and benchmark
  output.  The program must print the data to be compared in an output file in
  either a tabular format (default) or in a YAML format to standard output.
  Using YAML format requires the `PyYAML <http://pyyaml.org>`_ module to be
  installed.

  tabular format
      A row of text is assumed to start a table.  Multiple tables are permitted,
      but each table must be square (i.e.  no gaps and the same number of elements
      on each row) and hence each column heading must contain no spaces.  For
      example, a single table is of the format::

        val_1   val_2   val3
         1.2     2      3.32
         8.7     4      17.2

      and a table containing multiple subtables::

        val_1   val_2   val3
         1.2     2      3.32
         8.7     4      17.2
        val_4   val_5
        11.22   221.0

      Tables need not be beautifully presented: the amount of whitespace
      between each table cell is not important, so long as there's at least one
      space separating adjacent cells.

      Column headings are used to label the data in the subsequent rows.  These
      labels can be used to specify different tolerances for different types of
      data.

  YAML format
      The format accepted is a very restricted subset of YAML.  Specifically,
      only one YAML document is accepted and that document must contain
      a single block mapping.  Each key in the block mapping can contain
      a single data element to be compared or block sequence containing
      a series of data elements to be compared.  However, block sequences may
      not be nested.  The equivalent YAML formats for the two examples given
      above are::

          val_1:
            - 1.2
            - 8.7
          val_2:
            - 2
            - 4
          val_3:
            - 3.32
            - 17.2

      and::

          val_1:
            - 1.2
            - 8.7
          val_2:
            - 2
            - 4
          val_3:
            - 3.32
            - 17.2
          val_4: 11.22
          val_5: 221.0

      See the `PyYAML documentation
      <http://pyyaml.org/wiki/PyYAMLDocumentation>`_ for more details.

  Non-numerical values apart from the column headings in tabular ouput are
  required to be equal (within python's definition of equality for a given
  object).

* user-supplied verification program

  An external program can be used to validate the test output; the program must
  set an exit status of 0 to indicate the test passed and a non-zero value to
  indicate failure.
.. _config:

Configuration files
===================

For convenience, tests can be specified via configuration files rather than
using the testcode API directly.  These configuration files are required for
work with the command-line interface.

The two configuration files are, by default, :ref:`jobconfig` and
:ref:`userconfig` in the working directory.  Different names and/or paths can
be specified if required.

Both configuration files take options in the ini format (as understood by
Python's `configparser <http://docs.python.org/library/configparser.html>`_ module).  For example::

    [section_1]
    a = 2
    b = test_option

    [section_2]
    v = 4.5
    hello = world

defines an ini file with two sections (named 'section_1' and 'section_2'), each
with two variables set.

.. note::

    Any paths can either be absolute or relative to the directory containing
    the configuration file.  The full path need not be given for any program
    which exists on the user's PATH.  Environment variables in **program** names will be expanded.
.. _jobconfig:

jobconfig
=========

The jobconfig file defines the tests to run.  If a section named 'categories'
exists, then it gives labels to sets of tests.  All other sections are assumed
to individually define a test.

Tests
-----

A test is assumed to reside in the directory given by the name of the test
section.  For example::

    [carbon_dioxide_ccsd]
    inputs_args = ('co2.inp','')

would define a test in the ``carbon_dioxide_ccsd`` subdirectory relative to the
``jobconfig`` configuration file, with the input file as ``co2.inp`` (in the
``carbon_dioxide_ccsd`` subdirectory) with no additional arguments to be passed
to the test program.  All input and output files related to the test are
assumed to be contained within the test subdirectory.

The following options are permitted:

inputs_args [inputs and arguments format (see :ref:`below <inputs>`)]
    Input filename and associated arguments to be passed to the test program.
    No default.
min_nprocs [integer]
    Minimum number of processors to run test on.  Cannot be overridden by the
    '--processors' command-line option.  Default: 0.
max_nprocs [integer]
    Maximum number of processors to run test on.  Cannot be overridden by the
    '--processors' command-line option.  Default: 2^31-1 or 2^63-1.
nprocs [integer]
    Number of processors to run the test on.  Zero indicates to run the test
    purely in serial, without using an external program such as mpirun to
    launch the test program.  Default: 0.
output [string]
    Filename to which the output is written if the output is not written to
    standard output.  The output file is moved to the specific testcode test
    filename at the end of the calculation before the test output is validated
    against the benchmark output.  Wildcards are allowed so long as the pattern
    only matches a single file at the end of the calculation.  Default:
    inherits from setting in :ref:`userconfig`.
path [string]
    Set path (relative to the directory containing the ``jobconfig``
    configuration file) of the test.  The test is run in this directory and so
    input filenames need to be relative to it.  If the given path contains
    wildcards, then this is expanded and an individual test is created for each
    path that maches the pattern.  Note that Python's configparser restricts
    the use of special characters in section names and hence some patterns can
    only be accomplished by explicitly using the path option.  Default: test
    name (i.e.  the name of the section defining the test).
run_concurrent [boolean]
    If true then subtests defined by the inputs_args option are allowed to run
    concurrently rather than consecutively, assuming enough processors are
    available.  Default: false.
submit_template [string]
    Path to a template of a submit script used to submit jobs to a queueing
    system.  testcode will replace the string given in submit_pattern with the
    command(s) to run the test.  The submit script must do all other actions (e.g.
    setting environment variables, loading modules, copying files from the test
    directory to a local disk and copying files back afterwards).  No default.
program [string]
    Program name (appropriate section heading in :ref:`userconfig`) to use to
    run the test.  Default: specified in the [user] section of
    :ref:`userconfig`.
tolerance [tolerance format (see :ref:`tolerance`)]
    Tolerances for comparing test output to the benchmark output.  Default:
    inherits from the settings in :ref:`userconfig`.

If a test is defined via a category/path containing wildcards and explicitly,
then the explicit category will inherit any settings from the wildcard
definition.  For example, given the subdirectories ``t1`` and ``t2``, each
containing tests, the definition::

    [t*]
    inputs_args = ('test.in', '')
    [t1]
    nprocs = 2

is entirely equivalent to::

    [t1]
    nprocs = 2
    inputs_args = ('test.in', '')
    [t2]
    inputs_args = ('test.in', '')

.. note::

    Explicitly defining a test multiple times, e.g.::

        [t1]
        inputs_args = ('inp1', '')
        [t1]
        inputs_args = ('inp2', '')

    is not permitted and the resultant settings are not uniquely defined.

Test categories
---------------

For the purposes of selecting a subset of the tests in :ref:`testcode.py`, each
test is automatically placed in two separate categories, one labelled by the
test's name and the other by the test's path.  A test can hence be referred to
by either its path or by its name (which are identical by default).  

Additional categories can be specified in the [categories] section.  This makes
it very easy to select subsets of the tests to run.  For example::

    [categories]
    cat1 = t1 t2
    cat2 = t3 t4
    cat3 = cat1 t3

defines three categories (`cat`, `cat2` and `cat3`), each containing a subset
of the overall tests.  A category may contain another category so long as
circular dependencies are avoided.  There are two special categories, `_all_`
and `_default_`.  The `_all_` category contains, by default, all tests and
should not be changed under any circumstances.  The `_default_` category can
be set; if it is not specified then it is set to be the `_all_` category.

.. _inputs:

Program inputs and arguments
----------------------------

The inputs and arguments must be given in a specific format.  As with the
:ref:`tolerance format <tolerance>`,  the inputs and arguments are specified
using a comma-separated list of python tuples.  Each tuple (basically
a comma-separated list enclosed in parantheses) contains two elements: the name
of an input file and the associated arguments, in that order, represents
a subtest belonging to the given test.  Both elements must be quoted.  If the
input filename contains wildcard, then those wildcards are expanded to find all
files in the test subdirectory which match that pattern; the expanded list is
sorted in alphanumerical order.  A separate subtest (with the same arguments
string) is then created for each file matching the pattern.  used to construct
the command to run.  A null string (``''``) should be used to represent the
absence of an input file or arguments.  By default subtests run in the order
they are specified.  For example::

    inputs_args = ('test.inp', '')

defines a single subtest, with input filename ``test.inp`` and no arguments,

::

    inputs_args = ('test.inp', ''), ('test2.inp', '--verbose')

defines two subtests, with an additional argument for the second subtest, and

::

    inputs_args = ('test*.inp', '')

defines a subtest for each file matching the pattern ``test*inp`` in the
subdirectory of the test.
.. _userconfig:

userconfig
==========

The userconfig file must contain at least two sections.  One section must be
entitled 'user' and contains various user settings.  Any other section is
assumed to define a program to be tested, where the program is referred to
internally by its section name.  This makes it possible for a set of tests to
cover multiple, heavily intertwined, programs.  It is, however, far better to
have a distinct set of tests for each program where possible.

[user] section
--------------

The following options are allowed in the [user] section:

benchmark [string]
    Specify the ID of the benchmark to compare to.  This should be set running

    .. code-block bash

        $ testcode.py make-benchmarks

    The format of the benchmark files is'benchmark.out.ID.inp=INPUT_FILE.arg=ARGS'.  
    The 'inp' and/or 'arg' section is not included if it is empty.

    Multiple benchmarks can be used by providing a space-separated list of IDs.  The first
    ID in the list which corresponds to an existing benchmark filename is used to
    validate the test.
date_fmt [string]
    Format of the date string used to uniquely label test outputs.  This must
    be a valid date format string (see `Python documenation
    <http://docs.python.org/library/time.html>`_).  Default: %d%m%Y.
default_program [string]
    Default program used to run each test.  Only needs to be set if
    multiple program sections are specified.  No default.
diff [string]
    Program used to diff test and benchmark outputs.  Default: diff.
tolerance [tolerance format (see :ref:`below <tolerance>`.)]
    Default tolerance(s) used to compare all tests to their respective
    benchmarks.  Default: absolute tolerance 10^-10; no relative tolerance set.

[program_name] section(s)
-------------------------

The following options are allowed to specify a program (called 'program_name')
to be tested:

data_tag [string]
    Data tag to be used to extract data from test and benchmark output.  See
    :ref:`verification` for more details.  No default.
ignore_fields [space-separated list of strings]
    Specify the fields (e.g. column headings in the output from the extraction
    program) to ignore.  This can be used to include, say, timing information
    in the test output for performance comparison without causing failure of
    tests.  No default.
exe [string]
    Path to the program executable.  No default.
extract_args [string]
    Arguments to supply to the extraction program.  Default: null string. 
extract_cmd_template [string]
    Template of command used to extract data from output(s) with the following
    substitutions made:

        tc.extract
            replaced with the extraction program.
        tc.args
            replaced with extract_args.
        tc.file
            replaced with (as required) the filename of the test output or the
            filename of the benchmark output.
        tc.bench
            replaced with the filename of the benchmark output.
        tc.test
            replaced with the filename of the test output.

    Default: tc.extract tc.args tc.file if verify is False and
    tc.extract tc.args tc.test tc.bench if verify is True.
extract_program [string]
    Path to program to use to extract data from test and benchmark output.
    See :ref:`verification` for more details.  No default.
extract_fmt [string]
    Format of the data returned by extraction program. See :ref:`verification`
    for more details.  Can only take values table or yaml.  Default: table.
launch_parallel [string]
    Command template used to run the test program in parallel.  tc.nprocs is
    replaced with the number of processors a test uses (see run_cmd_template).
    If tc.nprocs does not appear, then testcode has no control over the number
    of processors a test is run on.  Default: mpirun -np tc.nprocs.
run_cmd_template [string]
    Template of command used to run the program on the test with the following
    substitutions made:

        tc.program
            replaced with the program to be tested.
        tc.args
            replaced with the arguments of the test.
        tc.input
            replaced with the input filename of the test.
        tc.output
            replaced with the filename for the standard output.  The filename
            is selected at runtime.
        tc.error
            replaced with the filename for the error output.  The filename is
            selected at runtime.
        tc.nprocs
            replaced with the number of processors the test is run on.

    Default: 'tc.program tc.args tc.input > tc.output 2> tc.error' in serial
    and 'launch_command tc.program tc.args tc.input > tc.output 2> tc.error' in
    parallel, where launch_command is specified above.  The parallel version is
    only used if the number of processors to run a test on is greater than
    zero.
skip_args [string]
    Arguments to supply to the program to test whether to skip the comparison
    of the test and benchmark.  Default: null string.
skip_cmd_template [string]
    Template of command used to test whether test was successfully run or
    whether the comparison of the benchmark and test output should be skipped.
    See :ref:`below <skip>` for more details.  The following strings in the
    template are replaced:

        tc.skip
            replaced with skip_program.
        tc.args
            replaced with skip_args.
        tc.test
            replaced with the filename of the test output.

    Default: tc.skip tc.args tc.test.
skip_program [string]
    Path to the program to test whether to skip the comparison of the test and
    benchmark.  If null, then this test is not performed.  Default: null string.
submit_pattern [string]
    String in the submit template to be replaced by the run command.  Default:
    testcode.run_cmd.
tolerance [tolerance format (see :ref:`below <tolerance>`.)]
    Default tolerance for tests of this type.  Default: inherits from
    [user].
verify [boolean]
    True if the extraction program compares the benchmark and test
    outputs directly.  See :ref:`verification` for more details.  Default:
    False.
vcs [string]
    Version control system used for the source code.  This is used to
    label the benchmarks.  The program binary is assumed to be in the same
    directory tree as the source code.  Supported values are: hg, git and svn
    and None.  If vcs is set to None, then the version id of the program is
    requested interactively when benchmarks are produced.  Default: None.

Most settings are optional and need only be set if certain functionality is
required or the default is not appropriate.  Note that either data_tag or
extract_program must be supplied.

In addition, the following variables are used, if present, as default settings
for all tests of this type:

* inputs_args (no default)
* nprocs (default: 0)
* min_nprocs (default: 0)
* max_nprocs (default: 2^31-1 or 2^63-1)
* output (no default)
* run_concurrent (defailt: false)
* submit_template
 
See :ref:`jobconfig` for more details.

All other settings are assumed to be paths to other versions of the program
(e.g. a stable version).  Using one of these versions instead of the one listed
under the 'exe' variable can be selected by an option to :ref:`testcode.py`.

.. _tolerance:

Tolerance format
----------------

The format for the tolerance for the data is very specific.  Individual
tolerance elements are specified in a comma-separated list.  Each individual
tolerance element is a python tuple (essentially a comma-separated list
enclosed in parentheses) consisting of, in order, the absolute tolerance, the
relative tolerance, the label of the field to which the tolerances apply and
a boolean value specifying the strictness of the tolerance (see below).  The
labels must be quoted.  If no label is supplied (or is set to None) then the
setting is taken to be the default tolerance to be applied to all data.  If the
strictness value is not given, the tolerance is assumed to be strict.  For
example, the setting::

    (1e-8, 1.e-6), (1.e-4, 1.e-4, 'Force')

uses an absolute tolerance of 10^-8 and a relative tolerance of 10^-6 by
default and an absolte tolerance and a relative tolerance of 10^-4 for data
items labelled with 'Force' (i.e. in columns headed by 'Force' using an
external data extraction program or labelled 'Force' by the internal data
extraction program using data tags).  If a tolerance is set to None, then it is
ignored.  At least one of the tolerances must be set.

A strict tolerance requires both the test value to be within the absolute and
relative tolerance of the benchmark value in order to be considered to pass.
This is the default behaviour.  A non-strict tolerance only requires the test
value to be within the absolute or relative tolerance of the benchmark value.
For example::

    (1e-8, 1e-6, None, False), (1e-10, 1e-10, 'Energy')

sets the default absolute and relative tolerances to be 10^-8 and 10^-6
respectively and sets the default tolerance to be non-strict except for the
'Energy' values, which have a strict absolute and relative tolerances of
10^-10.  If only one of the tolerances is set, then the strict and non-strict
settings are equivalent.

Alternatively, the tolerance can be labelled by a regular expression, in which case any
data labels which match the regular expression will use that tolerance unless there is
a tolerance with that specific label (i.e. exact matches override a regular
expression match).  Note that this is the case even if the tolerance using the exact
tolerance is defined in :ref:`userconfig` and the regular expression match is
defined in :ref:`jobconfig`.

.. _skip:

Skipping tests
--------------

Sometimes a test should not be compared to the benchmark---for example, if the
version of the program does not support a given feature or can only be run in
parallel.  testcode supports this by running a command to detect whether a test
should be skipped.

If the skipped program is set, then the skipped command is ran before
extracting data from output files.  For example, if

skip_program = grep
skip_args = "is not implemented."

are set, then testcode will run:

.. code-block:: bash

    grep "is not implemented." test_file

where test_file is the test output file.  If grep returns 0 (i.e.
test_file contains the string "is not implemented") then the test is
marked as skipped and the test file is not compared to the benchmark.
.. _testcode.py:

testcode.py
===========

.. only:: html

    testcode.py - a command-line interface to testcode.

Synopsis
--------

testcode.py [options] [action1 [action2...]]

Description
-----------

Run a set of actions on a set of tests.

Requires two configuration files, :ref:`jobconfig` and :ref:`userconfig`.  See
testcode documentation for further details.

testcode.py provides a command-line interface to testcode, a simple framework
for comparing output from (principally numeric) programs to previous output to
reveal regression errors or miscompilation.

Actions
-------

''run'' is th default action.

compare
    compare set of test outputs from a previous testcode run against the
    benchmark outputs.
diff
    diff set of test outputs from a previous testcode run against the benchmark
    outputs.
make-benchmarks
    create a new set of benchmarks and update the :ref:`userconfig` file with
    the new benchmark id.  Also runs the 'run' action unless the 'compare'
    action is also given.
recheck
    compare set of test outputs from a previous testcode run against
    benchmark outputs and rerun any failed tests.
run
    run a set of tests and compare against the benchmark outputs.
tidy
    Remove files from previous testcode runs from the test directories.

Options
-------

-h, --help
    show this help message and exit
-b BENCHMARK, --benchmark=BENCHMARK
    Set the file ID of the benchmark files.  If BENCHMARK is in the format
    t:ID, then the test files with the corresponding ID are used.  This
    allows two sets of tests to be compared.  Default: specified in the [user]
    section of the :ref:`userconfig` file.
-c CATEGORY, --category=CATEGORY
    Select the category/group of tests.  Can be specified multiple times.
    Wildcards or parent directories can be used to select multiple directories
    by their path.  Default: use the `_default_` category if run is an action
    unless make-benchmarks is an action.  All other cases use the `_all_`
    category by default.  The `_default_` category contains all  tests unless
    otherwise set in the :ref:`jobconfig` file.
-e EXECUTABLE, --executable=EXECUTABLE
    Set the executable(s) to be used to run the tests.  Can be  a path or name
    of an option in the :ref:`userconfig` file, in which case all test programs are
    set to use that value, or in the format program_name=value, which affects
    only the specified program.  Only relevant to the run action.  Default: exe
    variable set for each program listed in the :ref:`userconfig` file.
-i, --insert
    Insert the new benchmark into the existing list of benchmarks in userconfig
    rather than overwriting it.  Only relevant to the make-benchmarks action.
    Default: False.
--jobconfig=JOBCONFIG
    Set path to the job configuration file.  Default: jobconfig.
--job-option=JOB_OPTION
    Override/add setting to :ref:`jobconfig`.  Takes three arguments.  Format:
    section_name option_name value.  Default: none.
--older-than=OLDER_THAN
    Set the age (in days) of files to remove.  Only relevant to the tidy
    action.  Default: 14 days.
-p NPROCS, --processors=NPROCS
    Set the number of processors to run each test on.  Only relevant to the run
    action.  Default: run tests as serial jobs.
-q, --quiet
    Print only minimal output.  Default: False.
-s QUEUE_SYSTEM, --submit=QUEUE_SYSTEM
    Submit tests to a queueing system of the specified type.  Only PBS system
    is currently implemented.  Only relevant to the run action.  Default: none.
-t TEST_ID, --test-id=TEST_ID
    Set the file ID of the test outputs.  If TEST_ID is in the format b:ID, then
    the benchmark files with the corresponding ID are used.  This allows two
    sets of benchmarks to be compared.  Default: unique filename based upon
    date if running tests and most recent test_id if comparing tests.
--total-processors=TOT_NPROCS
    Set the total number of processors to use to run as many tests as possible
    at the same time.  Relevant only to the run option.  Default: run all tests
    concurrently run if --submit is used; run tests sequentially otherwise.
--userconfig=USERCONFIG
    Set path to the user configuration file.  Default: userconfig.
--user-option=USER_OPTION
    Override/add setting to :ref:`userconfig`.  Takes three arguments.  Format:
    section_name option_name value.  Default: none.
-v, --verbose
    Increase verbosity of output.  Can be specified up to two times.
    The default behaviour is to print out the test and its status.  (See the
    --quiet option to suppress even this.)  Specify -v or --verbose once to
    show (if relevant) which data values caused warnings or failures.
    Specify -v or --verbose twice to see all (external) commands run and all
    data extracted from running the tests.  Using the maximum verbosity level
    is highly recommended for debugging.

Exit status
-----------

1 if one or more tests fail (run and compare actions only) and 0 otherwise.

License
-------

Modified BSD License.  See LICENSE in the source code for more details.

Bugs
----

Contact James Spencer (j.spencer@imperial.ac.uk) regarding bug reports,
suggestions for improvements or code contributions.
testcode
========

testcode is a python module for testing for regression errors in numerical
(principally scientific) software.  Essentially testcode runs a set of
calculations, and compares the output data to that generated by a previous
calculation (which is regarded to be "correct").  It is designed to be
lightweight and highly portable: it can be used both as part of the development
process and to verify the correctness of a binary on a new architecture.
testcode requires python 2.4-3.4.  If these are not available, then `pypy
<http://www.pypy.org>`_ is recommended---for this purpose pypy serves as
a portable, self-contained python implementation but this is a tiny aspect of
the pypy project.

testcode can run a set of tests and check the calculated data is within a the
desired tolerance of results contained in previous output (using an internal
data extraction engine, a user-supplied data extraction program or
a user-supplied verification program).  The programs to be tested can be run in
serial and in parallel and tests can be run in either locally or submitted to
a compute cluster running a queueing system such as PBS.  Previous tests can be
compared and diffed against other tests or benchmarks.

testcode provides access to these features via an API.  The supplied
command-line interface, :ref:`testcode.py`, should be sufficient for most
purposes.  The command-line interface utilises simple :ref:`configuration files
<config>`, wich makes it easy to customise to the local environment and to add
new tests.

.. toctree::
   :maxdepth: 1

   installation
   configuration_files
   jobconfig
   userconfig
   verification
   testcode.py

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Installation
============

testcode2 is designed to be very lightweight and portable, so it can easily and
quickly be used on a variety of machines.  Typically only downloading the
testcode2 package is required.

If the :ref:`testcode.py` script is used, then no additional installation steps
are required assuming the directory structure is preserved.  If the
``testcode2`` module is used or the files are split up and installed elsewhere,
then the ``testcode2`` module must be able to be found by python (i.e. exists
on $PYTHONPATH).
