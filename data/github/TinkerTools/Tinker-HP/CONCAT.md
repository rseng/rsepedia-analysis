# Tinker-HP: High-Performance Massively Parallel Evolution of Tinker on CPUs & GPUs
-----------------------------------------------------------------------------------------------------------------------------------------------
A phase advance update (GPUs) has been pushed to GitHub (28/10/2020) to support COVID-19 research: the Tinker-HP (multi)-GPUs plateform is now available: https://github.com/TinkerTools/tinker-hp/tree/master/GPU 

!!! Check out the JCTC paper (Open Access) about the native GPUs version of Tinker-HP (03/23/2021) : 
Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems. O. Adjoua, L. Lagardère, L.-H. Jolly, Arnaud Durocher, Z. Wang, T. Very, I. Dupays, T. Jaffrelot Inizan, F. Célerse, P. Ren, J. Ponder, J-P. Piquemal, J. Chem. Theory. Comput., 2021, 17 (4), 2034–2053 (Open Access) https://doi.org/10.1021/acs.jctc.0c01164 
-----------------------------------------------------------------------------------------------------------------------------------------------
<H2><B>Versions</B></H2>
<B>Update 04/12/2021 : In addition to GitHub, a GPUs container (quick install!) is available thanks to NVIDIA on the NVIDIA NGC's website: https://ngc.nvidia.com/catalog/containers/hpc:tinkerhp </B>


<B>Update 03/02/2021 : PLUMED Support for GPUs !</B>

<B>Update 02/11/2021 : PLUMED Support for version 1.2 (CPUs) !</B>

<B>Update 11/24/2020 : all versions have been pushed to GitHub</B>

Current Github version: 1.1v (enhanced AVX512 vectorized CPUs version), 1.2 (CPUs) + 1.2/phase advance (multi)-GPUs

Current Development version: 1.3 (CPUs + multi-GPUs)


All releases of the Tinker-HP code are now being performed on Github. For news, benchmarks and additional tutorials, please visit the Tinker-HP website
http://tinker-hp.ip2ct.upmc.fr/   (new website design to appear)

Tinker-HP is a CPUs and GPUs based, multi-precision, MPI massively parallel package dedicated to long polarizable molecular dynamics simulations and to polarizable QM/MM. Tinker-HP is an evolution of the popular Tinker package that conserves it simplicity of use but brings new 
capabilities allowing performing very long molecular dynamics simulations on modern supercomputers that use thousands of cores. 
The Tinker-HP approach offers various strategies using domain decomposition techniques for periodic boundary conditions in the 
framework of the (n)log(n) Smooth Particle Mesh Ewald. Tinker-HP proposes a high performance scalable computing environment for 
polarizable (AMOEBA, Amberpol...) and classical (Amber, Charmm, OPLS...) force fields giving access to large systems up to millions of atoms. It can be used on supercomputers as well as on lab clusters. Tinker-HP supports Intel (AVX5212 enhanced version) and AMD CPUs platforms as well as NVIDIA GPUs (1080, 2080, 3090, P100, V100, A100). 

Tinker-HP is available free of charge for ALL Academic Institutions, National Laboratories and supercomputer centers through the global Tinker license (https://dasher.wustl.edu/tinker/downloads/license.pdf).
Non-academic entities (e.g., companies, for profit organizations) should contact the managing universities (see license).

If you want support, support is only provided to registered users:

i) <B>Please fill in the form at:</B>
https://dasher.wustl.edu/tinker/downloads/license.pdf

and sent it back to TinkerHP_Download@ip2ct.upmc.fr to be added to the user support.

ii) <B>Please cite:</B>

- Tinker-HP: a Massively Parallel Molecular Dynamics Package for Multiscale Simulations of Large Complex Systems 
with Advanced Polarizable Force Fields.
L. Lagardère, L.-H. Jolly, F. Lipparini, F. Aviat, B. Stamm, Z. F. Jing, M. Harger, H. Torabifard, G. A. Cisneros, 
M. J. Schnieders, N. Gresh, Y. Maday, P. Ren, J. W. Ponder, J.-P. Piquemal, Chem. Sci., 2018, 9, 956-972 (Open Access)
https://doi.org/10.1039/C7SC04531J

- if you are using the AVX512 vectorized version (1.1v) dedicated to Intel's CPUs (Skylake, CascadeLake etc...), please also cite:
Raising the Performance of the Tinker-HP Molecular Modeling Package [Article v1.0].
L. H. Jolly, A. Duran, L. Lagardère, J. W. Ponder, P. Y. Ren, J.-P. Piquemal, LiveCoMS, 2019, 1 (2), 10409  (Open Access)
 https://doi.org/10.33011/livecoms.1.2.10409
 
- if you are using the GPUs version, please also cite:
Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems
Olivier Adjoua, Louis Lagardère, Luc-Henri Jolly, Arnaud Durocher, Thibaut Very, Isabelle Dupays, Zhi Wang, Théo Jaffrelot Inizan, Frédéric Célerse, Pengyu Ren, Jay W. Ponder, Jean-Philip Piquemal, J. Chem. Theory. Comput., 2021, 17 (4), 2034–2053(Open Access) https://doi.org/10.1021/acs.jctc.0c01164

iii) <b>Tinkertools</b>
Tinker-HP is part of the Tinker distribution and uses the same tools as Tinker. These tools can be found here : https://github.com/TinkerTools/tinker
if you use the Tinkertools please cite :

Tinker 8: Software Tools for Molecular Design.
J. A. Rackers, Z. Wang, C. Lu, M. L. Maury, L. Lagardère, M. J. Schnieders, J.-P. Piquemal, P. Ren, J. W. Ponder,  J. Chem. Theory. Comput., 2018, 14 (10), 5273–5289 DOI: http://dx.doi.org/10.1021/acs.jctc.8b00529 ; PMC free text : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6335969/

iii) <B>Support:</B>

We provide support to registered users only (see i) ).

Email: TinkerHP_Support@ip2ct.upmc.fr


# Prerequisites

## Hardware
A relatively recent Nvidia GPU is mandatory for the GPU code. One with at least a compute capability above 5.0.

Nothing special is needed for the CPU code.


## Compilers, Operating Systems and Environment
   - A GNU or Intel FORTRAN compiler for the CPU code (Not recommended at this point)
   - The most recent [PGI compiler](https://www.pgroup.com/products/community.htm).
   - The nvidia CUDA C/C++ compiler for CUDA Files in GPU code
   - MPI Fortran wrapper on the selected compiler and compiled to be CUDA-Aware for device communications during multi-GPUs execution.
   - Linux or Windows 10 (Windows Subsystem for Linux) are preferred  
     __Note :__ PGI compiler does not support the GPU code on macOS.
   - [Nvidia HPC Package](https://developer.nvidia.com/nvidia-hpc-sdk-releases) Any version under __21.2 !__  
     An issue involving "atomic operations" has been discovered with higher versions of the package (21.[3-7]). Luckily it is possible to bypass it by adding a compiler flag during the build step. Look for `add_options` variable configuration in `ci/install.sh` to be revert to the old behavior and fix the problem. However be aware that this is just a temporary fix.  
     This package contains everything you need to build Tinker-HP (GPU part). Previously listed items are already available within it. However we need to make sure the installation has been correctly done and the environment variables (PATH; CPATH & NVHPC) are correctly set and updated. Follow instructions described in `$PREFIX/modulefiles` with `$PREFIX` the installation directory if not.


## Mandatory Libraries
   - FFt libraries
      - `libfftw3` or `libfftw3f` in single or mixed precision mode.
        You will need to provide in that case FFTW install directory to both Tinker's and 2decompfft's Makefile.
        Exporting `FFTW=/path_to_fftw3` as an environment variable should be enough to make.
      - It is also possible to use the generic fft provided with your compiler (Default behavior).
   - Intel MKL Library (Only for host/CPU build!)
      - Export `MKLROOT=/path_to_mkl_library` inside your shell environment
      - `lmkl_intel_lp64 lmkl_sequential lmkl_core` are required
   - `2decomp_fft` library already provided inside the repository.
   - `libstdc++` Standard C++ library
      - Export `GNUROOT=/path_to_gnu`.
      This is not mandatory for a first and unique installation with `install.sh` but it is required for developers or when you need to recompile in an other configuration.
      e.g. `export GNUROOT=/usr` for conventional Linux systems
   - `CUDA` (Toolkit and libraries)
      - Cuda(C/C++) compiler `nvcc`
      - Recent CUDA libraries are required by the GPU code. Every CUDA version since 9.1 has been tested.
        They are delivered with NVIDIA HPC-SDK packages

### Notes
It is important to make sure that your environments (both build and run) are version consistent. In most cases, we have to pay attention our native CUDA version does not differ from PGI provided CUDA and/or NVIDIA driver installed. This matter will not happen with Nvidia HPC Package installed and loaded in your environment.
# Build Tinker-HP (GPU)

### Easy Build
A relatively easy way to build GPU version of Tinker-HP is to use the installation bash script. After setting your environment according to the Prerequisites, you can proceed to installation by typing in your shell.
```
$> pwd
#> /home/user/.../tinker-hp/GPU
$> ci/install.sh
```

As long as you stand for one configuration at the time, you can familiarize with the __source/Makefile.pgi's__ targets to customize your own build along.
In the event of an installation scipt failure, we recommend to follow the regular build procedure given below.

### Configuration options
   All those options are disabled by default.
   - `NVTX_SUPPORT` enables NVTX markers during profiling if set to 1. This option is not useful in release construct.
   - `NVSHMEM_SUPPORT` enables build with nvshmem library if set to 1. This option require `NVSHMEM_HOME` variable to be set on nvshmem install directory. It has not been tested with recent Nvidia HPC package.
   - `NO_MUTATION` disables soft-core computation if set to 1
   - `FPA_SUPPORT` enables fixed point arithmetic in non-Bonded Forces and energy reduction if set to 1. It requires mixed precision to be enabled through the Makefile precision variable.

You can force options at compile by editing in `source/Makefile.pgi:l62`. A proper way will be introduced in next developments.

### Using Makefile
You can almost do any construct in the code with `source/Makefile` linked to `source/Makefile.gpi`. Here are some useful tips to configure your build's type

##### configuration variables

  - `arch=(host|[device]) .or. (cpu|[gpu])` allows you to select target architecture
  - `prec=(single|mixed|[double]) .or. (s|m|[d])` changes the precision build
  - `prefix=(path/|[../bin])` controls the installation directory
  - `prog_suffix=(any_string|[])` append suffix to install binary. This variable may allow you to have multiple binaries in the installation directory.  
  For instance `analyze` or `analyze.gpu`
  - `compute_capability=(list|[60,70])` selects the device's compute capability for the construct. this variable shall accept a list of compute capability separated by a comma.   
  e.g. `compute_capability=35,60,75`
  - `cuda_version=([10.1])` contains the cuda version to be used for construct.  
    _PGI 19.10_ compilers only supports _cuda 9.2 10.0 and 10.1_ for OpenACC. To be consistent with the construct, we strongly recommend to use the same version of CUDA compiler as your OpenACC compiler.
  - `opt=(debug|debug+|[release])` decides the optimization level of the compilation

##### Some targets

  - `create_build`  
    Used In association with BUILD_DIR will create a new build directory
    with links to source files and _source/Makefile.pgi_. When BUILD_DIR is not specified, the default directory is `build/`. This target become useful when you want to keep object's files from a different build.  
    N.B. Depending on the precision you might have to rebuild the 22decomp_fft library.  
    e.g.  `make create_build BUIL_DIR=../build-mix`

  - `2decomp_fft_rebuild` ;  
    This phony target allow you to rebuild 2decomp_fft library in double precision
    Since the compilation uses the same directory for 2decomp_fft, you might want to
    be sure that your library is being compiled in the correct precision

  - `2decomp_fft_rebuild_single` ;  
    Same as previous target but excpet from the precision's construct which happen here to in single. There is only two precision modes for 2decomp_fft Library. Tinker mixed precision building requires 2decomp_fft single precision library.

  - `thrust_lib_rebuild` ;  
    [Re]Build the wrapper on CUDA thrust library

  - `all` ;  
    compile link and install `bin/dynamic` `bin/analyze` and `bin/minimize`. If you just want one program to be build just enter `make <program>`.  
    e.g. `make analyze dynamic`

  - `[dynamic|analyze|minimize].mixed` ;  
    compiles link and install analyze|dynamic|minimize.mixed. Those specific targets calls behind `make` on their main target. Only used one at the time to avoid over building.  
    e.g. `make analyze.mixed -j4` is fine

  - `[dynamic|analyze|minimize].single` ;  
    compiles link and install analyze|dynamic|minimize.single.

  - `[dynamic|analyze|minimize].cpu` ;  
    compiles link and install analyze|dynamic|minimize.cpu. CPU binaries

  - `libtinker` builds `libtinker.a`


### Custom Build
Let suggest that we want CPU binaries in double precision combined with an out-of-source build by hand. After linking `source/Makefile.pgi` to `source/Makefile`, we should follow the next script for an out-of-source build :
```
$> pwd
#> /home/user/PME
$> mkdir -p bin
$> cd source
$> make arch=cpu 2decomp_fft_rebuild          # build lib2decompfft.a
$> make thrust_lib_rebuild                    # build libwapper.a
$> make create_build BUILD_DIR=../build_CPU   # create build directory
$> cd ../build_CPU
$> make arch=cpu all -j6                      # build Tinker's program for CPU
$> ls ../bin                                  # list binaries
#> analyze dynamic minimize bar
```
or the following one for an in-source build
```
$> cd source
$> make arch=cpu prog_suffix=.cpu all -j4
$> ls ../bin
#> analyze.cpu dynamic.cpu minimize.cpu bar.cpu
```

You can also create a configuration bash file that store your configuration and run it with your desire targets.  
```
$> cat << EOF >Tconf.sh
   make FPA_SUPPORT=1 prec=mixed arch=device compute_capability=60 prog_suffix=.gmix $@
   EOF
$> chmod 740 Tconf.sh
$> cd source
$> ln -s ../Tconf.sh T
$> ./T 2decomp_fft_rebuild_single  # Rebuild 2decomp_fft library in single precision
$> ./T thrust_lib_rebuild
$> ./T all -j4    # Build binaries (*.gfix) for 60 compute_capability device using fixed precision
```
_______________________________
Tinker-HP: High Performance Multi-GPUs Massively Parallel Evolution of Tinker
==================================================================


<b>This phase-advance GPU version (1.2 ++) is not (yet) an official release of Tinker-HP but is made freely available in link with the COVID-19 HPC community effort.</b>

This work will be part of a larger 2021 Tinker-HP 1.3 official release.
In addition to GitHub, a GPUs container (quick install!) is available thanks to NVIDIA on the NVIDIA NGC's website: https://ngc.nvidia.com/catalog/containers/hpc:tinkerhp

# Getting started with Tinker-HP
   - Installation Guide

## Installation Guide
   -  [Prerequisites](Prerequisites.md)
   -  [Build Tinker-HP (GPU version)](build.md)

## Run Tinker-HP (CPU/GPU)
There is no difference between the use of Tinker-HP and Tinker-HP (GPU version) as long as the feature you are looking for is available on the GPU version. The present version is optimized to accelerate simulations using the AMOEBA polarizable force field. Some minimal non-polarizable capabilities are present (enhanced support will be available in 2021). The code has been extensively tested on 1080, 2080, 3090, P100, V100 and A100 NVIDIA GPU cards and support multi-GPUs computations. It will be part of the major Tinker-HP 1.3 2021 release but this present version will continue to evolve. 

### GPU available features
   - dynamic analyze minimize and bar programs
   - Integrators (RESPA, RESPA1, BAOAB, BAOAB-RESPA1, VERLET)
   - Amoeba polarizable force field, classical force fields (AMBER/CHARMM/OPLS)
   - New implementation of PCG and DC-DIIS solver for polarization (DC-DIIS is not adapted to the device, use PCG instead)
   - Bussi Thermostat for NVT simulations  (it is default)
   - Montecarlo and Berendsen barostat for NPT simulations (default is Berendsen)
   - Accelerate Molecular Dynamics : aMD and GaMD Simulations
   - Steered Molecular Dynamics (SMD)
   - Orthogonal and Octahedron PBC Box Shapes (the latest to be used only on a single MPI process for now)
   - Plumed support available (updated : 03/2021)
   - More to come

   For more detailed informations on how to use the application, see section V to VIII of the readme of the CPU 1.2 version: https://github.com/TinkerTools/tinker-hp/blob/master/v1.2/Readme_v1.2.pdf
   Beware that not all the features available on CPU are available on GPU (see above, for example the TCG solver is only available on CPU as is the Langevin Piston barostat).
   
   <B>If you use the code please cite :</B>
   
   Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems.
O. Adjoua,  L. Lagardère, L.-H. Jolly, Arnaud Durocher, Z. Wang, T. Very, I. Dupays, T. Jaffrelot Inizan, F. Célerse, P. Ren, J. Ponder, J-P. Piquemal, J. Chem. Theory. Comput., 2021, 17 (4), 2034–2053 (Open Access) https://doi.org/10.1021/acs.jctc.0c01164
   
   and 
   
Tinker-HP: a Massively Parallel Molecular Dynamics Package for Multiscale Simulations of Large Complex Systems with Advanced Polarizable Force Fields.
L. Lagardère, L.-H. Jolly, F. Lipparini, F. Aviat, B. Stamm, Z. F. Jing, M. Harger, H. Torabifard, G. A. Cisneros, M. J. Schnieders, N. Gresh, Y. Maday, P. Ren, J. W. Ponder, J.-P. Piquemal, Chem. Sci., 2018, 9, 956-972 (Open Access) https://doi.org/10.1039/C7SC04531J

<B>License :</B> 

Tinker-HP is available free of charge for ALL Academic Institutions, National Laboratories and supercomputer centers through the global Tinker license (https://dasher.wustl.edu/tinker/downloads/license.pdf). Non-academic entities (e.g., companies, for profit organizations) should contact the managing universities (see license).

<B>Tinkertools :</B> Tinker-HP is part of the Tinker distribution and uses the same tools as Tinker. These tools can be found here : https://github.com/TinkerTools/tinker if you use the Tinkertools please cite :

Tinker 8: Software Tools for Molecular Design. J. A. Rackers, Z. Wang, C. Lu, M. L. Maury, L. Lagardère, M. J. Schnieders, J.-P. Piquemal, P. Ren, J. W. Ponder, J. Chem. Theory. Comput., 2018, 14 (10), 5273–5289 DOI: http://dx.doi.org/10.1021/acs.jctc.8b00529 PMC free text : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6335969/

<B>Support :</B>

We provide support to registered users only (http://tinker-hp.ip2ct.upmc.fr/?Download-instructions).

Email: TinkerHP_Support@ip2ct.upmc.fr

<B>Funding :</B> 
- this work has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 810367), project EMC2 (see preprint for full acknowledgments)

- we thank GENCI, NVIDIA and HPE as well as the engineering team of the IDRIS Supercomputer center (CNRS/GENCI, France). 
[![Homepage](https://img.shields.io/badge/Home-plumed.org-green.svg)](http://www.plumed.org)
[![Homepage](https://img.shields.io/badge/Google_group-plumed--users-green.svg)](http://groups.google.com/forum/#!forum/plumed-users)
[![codecov](https://codecov.io/gh/plumed/plumed2/branch/master/graph/badge.svg)](https://codecov.io/gh/plumed/plumed2)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/plumed/plumed2.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/plumed/plumed2/context:python)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/plumed/plumed2.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/plumed/plumed2/context:cpp)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](http://www.gnu.org/licenses/lgpl-3.0)
[![Github Releases](https://img.shields.io/github/release/plumed/plumed2.svg)](https://github.com/plumed/plumed2/releases)
[![MacPorts package](https://repology.org/badge/version-for-repo/macports/plumed.svg)](https://repology.org/project/plumed/versions)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/plumed/badges/version.svg)](https://anaconda.org/conda-forge/plumed)
[![AUR package](https://repology.org/badge/version-for-repo/aur/plumed.svg)](https://repology.org/project/plumed/versions)
[![Twitter Follow](https://img.shields.io/twitter/follow/plumed_org.svg?style=social&label=Follow)](https://twitter.com/plumed_org)

Branches and releases
---------------------

Several branches and tags are stored on the git repository.

Branches named `v2.X` correspond to release branches.

Master branch may contain non tested features and is not expected to be used by non-developers.
It typically contains features that will be available on the next release.

Tags named `v2.XbY` correspond to beta releases, use it with care.
Tags named `v2.X.Y` correspond to official releases, use the latest available.

In addition, the repository contains a number of other branches related to specific features.
Please contact the developers that are committing on those branches before basing your work
there, since they might contain temporary work and might be rebased later.
For instance, branch `testdoc` is setup so as to push a test copy of the manual
and is often force pushed.

To report problems found on beta or official releases, use the normal
[plumed-users@googlegroups.com](mailto:plumed-users@googlegroups.com)
mailing list. Please state exactly which version you are using.
To report problems found on `master` branch, use the
[plumed2-git@googlegroups.com](plumed2-git@googlegroups.com) mailing list.
This is also the correct place for discussions about new features etc.
When reporting please provide the git hash (you can obtain it with `git rev-parse HEAD`).

Status
------

Below you find the status on [Travis-CI](http://travis-ci.org/plumed/plumed2) for the release branches.

| Branch   |      Status   | First stable release (year) | Still supported |
|:--------:|:-------------:|:--------:|:------:|
| master   | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=master)](https://travis-ci.org/plumed/plumed2) | 2020 (expected) | / |
| v2.6     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.6)](https://travis-ci.org/plumed/plumed2)   | 2019 | yes |
| v2.5     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.5)](https://travis-ci.org/plumed/plumed2)   | 2018 | yes |
| v2.4     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.4)](https://travis-ci.org/plumed/plumed2)   | 2017 | no |
| v2.3     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.3)](https://travis-ci.org/plumed/plumed2)   | 2016 | no |
| v2.2     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.2)](https://travis-ci.org/plumed/plumed2)   | 2015 | no |
| v2.1     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.1)](https://travis-ci.org/plumed/plumed2)   | 2014 | no |
| v2.0     | Not available | 2013 | no |

Content
-------

Here's a description of the content of each file and directory in the root PLUMED directory.

    CHANGES          : change log
    COPYING.LESSER   : license
    Makefile         : makefile
    Makefile.conf.in : template configuration makefile
    PEOPLE           : list of authors
    README.md        : this file
    VERSION          : version file
    astyle           : a local version of astyle, used to format code
    configure        : configuration script
    configure.ac     : configuration script (autoconf)
    developer-doc    : developer documentation
    docker           : directory where Docker is generated
    macports         : directory where Portfiles are generated
    patches          : patch scripts
    python           : python stuff
    regtest          : regression tests, including reference results
    release.sh       : developer utility to publish releases
    scripts          : shell tools
    sourceme.sh.in   : template configuration script
    src              : source code
    test             : examples
    user-doc         : user documentation
    vim              : directory where vim syntax is generated

Required software
-----------------

Required software:

* GNU make.
* C/c++ compiler (c++11 support is required as of version 2.4).
* A modern version of the `patch` command line tool.
* Support for POSIX library `dirent.h`.
* `xxd` (present in most UNIX distributions).

Suggested software (libraries are checked by `./configure` and enabled if available):

* MPI library to run parallel simulations. It should be the same library used by your MD code.
* Optimized blas and lapack libraries. They are automatically replaced by an internal version if not available.
* [VMD molfile plugins](http://www.ks.uiuc.edu/Research/vmd/plugins) to read arbitrary file formats. They are automatically replaced by an internal version supporting a few formats if not available.
* [Zlib library](http://zlib.net/) to use compressed data files.
* [Xdrfile library](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library) to have read/write access to gromacs
  trajectory files.
* [Doxygen](http:://www.doxygen.org) to build user manual. Doxygen might need the following packages:
  * Latex to build the pdf user manual.
  * [Graphviz](http://www.graphviz.org) to show class hierarchy in
    developer manual.

Quick compilation instructions
------------------------------

Extensive installation instructions are in the [user documentation](http://www.plumed.org/documentation).
Quick instructions:

    ./configure --prefix=$HOME/opt
    make
    make doc # optional
    make test # optional

User documentation can be found at `user-doc/html/index.html`.
Developer documentation can be found at `developer-doc/html/index.html`.
[Pre-compiled documentation](http://www.plumed.org/documentation) is available online, so this is only required
if you are working with a modified version of the code!

In order to run PLUMED without installing it you should type `source sourceme.sh`. However,
we recomment installing PLUMED. 
To install it in `$HOME/opt` (directory should be set during `./configure`):

    umask 022
    make install
    
Now you will be able to run plumed using e.g.

    plumed help

If you compiled your own documentation, paths to the installed documentation can be found with command `plumed info --user-doc`.

A sample modulefile with environment variable will be placed in
`$HOME/opt/lib/plumed/src/lib/modulefile`. This can be useful if you want to
install multiple PLUMED versions side by side and select them with env modules.
Instructions for using Artistic Style are included in the *doc* directory.

The file **install.html** contains instructions for compiling and
installing Artistic Style.

The file **astyle.html**' contains information on using Artistic Style.

The files **news.html** and **notes.html** contain information on changes
made to the various releases.
ANN (Artificial Neural Network) function for plumed
====================

This is plumed ANN function (annfunc) module.  It implements `ANN` class, which is a subclass of `Function` class.  `ANN` class takes multi-dimensional arrays as inputs for a fully-connected feedforward neural network with specified neural network weights and generates corresponding outputs.  The `ANN` outputs can be used as collective variables, inputs for other collective variables, or inputs for data analysis tools.  

## Installation

Enable compilation by adding the `--enable-modules=annfunc` to the configure command.

## Usage

It is used in a similar way to [other plumed functions](https://www.plumed.org/doc-v2.5/user-doc/html/_function.html).  To define an `ANN` function object, we need to define following keywords:

- `ARG` (string array): input variable names for the fully-connected feedforward neural network

- `NUM_LAYERS` (int): number of layers for the neural network

- `NUM_NODES` (int array): number of nodes in all layers of the neural network

- `ACTIVATIONS` (string array): types of activation functions of layers, currently we have implemented "Linear", "Tanh", "Circular" layers, it should be straightforward to add other types as well

- `WEIGHTS` (numbered keyword, double array): this is a numbered keyword, `WEIGHTS0` represents flattened weight array connecting layer 0 and layer 1, `WEIGHTS1` represents flattened weight array connecting layer 1 and layer 2, ...  An example is given in the next section.

- `BIASES` (numbered keyword, double array): this is a numbered keyword, BIASES0 represents bias array for layer 1, BIASES1 represents bias array for layer 2, ...

Assuming we have an `ANN` function object named `ann`, we use `ann.node-0, ann.node-1, ...` to access component 0, 1, ... of its outputs (used as collective variables, inputs for other collective variables, or data analysis tools).

## Examples

Assume we have an ANN with numbers of nodes being [2, 3, 1], and weights connecting layer 0 and 1 are

```
[[1,2],
[3,4],
[5,6]]
```

weights connecting layer 1 and 2 are

```
[[7,8,9]]
```

Bias for layer 1 and 2 are

```
[10, 11, 12]
```

and 

```
[13]
```

respectively.

All activation functions are `Tanh`.

Then if input variables are `l_0_out_0, l_0_out_1`, the corresponding `ANN` function object can be defined using following plumed script: 

```
ann: ANN ARG=l_0_out_0,l_0_out_1 NUM_LAYERS=3 NUM_NODES=2,3,1 ACTIVATIONS=Tanh,Tanh  WEIGHTS0=1,2,3,4,5,6 WEIGHTS1=7,8,9  BIASES0=10,11,12 BIASES1=13
```

This plumed script can be generated with function `Plumed_helper.get_ANN_expression()` in [this](https://github.com/weiHelloWorld/plumed_helper/blob/master/plumed_helper.py) repository.  Following is the Python code using this function to generate the script above:

```Python
from plumed_helper import Plumed_helper
ANN_weights = [np.array([1,2,3,4,5,6]), np.array([7,8,9])]
ANN_bias = [np.array([10, 11, 12]), np.array([13])]
Plumed_helper.get_ANN_expression('ANN', node_num=[2, 3, 1], 
                                 ANN_weights=ANN_weights, ANN_bias=ANN_bias,
                                 activation_list=['Tanh', 'Tanh'])
```

## Authors

Wei Chen (UIUC, weichen9@illinois.edu) and Andrew Ferguson (University of Chicago, andrewferguson@uchicago.edu)

## Copyright

See ./COPYRIGHT
maze
================================================================================
This version of PLUMED2 has the maze module implemented.

Install
--------------------------------------------------------------------------------
Enable the compilation of maze by adding the `--enable-modules=maze` to the 
configure command.

Page
--------------------------------------------------------------------------------
See [this link](http://maze-code.github.io) for further information.

Documentation
--------------------------------------------------------------------------------
Run `make doc`; the documentation should be in `user-doc/html/_m_a_z_e.html`.

Author
--------------------------------------------------------------------------------
Jakub Rydzewski (Nicolaus Copernicus University) <jr@fizyka.umk.pl>

Copyright
--------------------------------------------------------------------------------
This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU Lesser General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any 
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  

See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along 
with this program.  If not, see <http://www.gnu.org/licenses/>.
Welcome to the plumed2 wiki!

This fork of PLUMED has eABF/DRR implementation.

**Requirements**

Boost::serialization and C++11 compiler

**Compiling instruction:**

After clone this repository, please cd to the plumed2 directory and run:

1. `autoconf`
2. `./configure --enable-boost_serialization --enable-modules=drr`
3. Modify your Makefile.conf and add `-lboost_serialization` in `DYNAMIC_LIBS=`
4. `make` and `sudo make install`

**Usage**

Run `make doc` and the usage is in `user-doc/html/_d_r_r.html`

**Authors**

Chen Haochuan (All files in drr module except colvar_UIestimator.h)

Fu Haohao (colvar_UIestimator.h)

**COPYRIGHT**

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Experiment Directed Simulation (EDS)
====================================


Install
------------------------------------
Enable compilation by adding the `--enable-modules=+eds`
to the configure command.


Documentation
------------------------------------
See the generated documentation for information on
using


Authors
------------------------------------
Glen Hocky (University of Chicago) <hockyg@uchicago.edu>
Andrew White (University of Rochester) <andrew.white@rochester.edu>


Copyright
------------------------------------
See ./COPYRIGHT
@page CHANGES-2-6 Version 2.6
  
## Version 2.6 (Jan 27, 2020)

Changes from version 2.5 which are relevant for users:
- Changes leading to incompatible behavior:
  - PLUMED input file parsing is now case insensitive that is that all directives can be written using uppercase characters (compatible with former versions) as well as lowercase characters (not compatible) internally PLUMED still uses uppercase definitions
  - `plumed partial_tempering` now uses `gawk` instead of `awk`. You might need to install `gawk` for it to work correctly.

- Other changes:
  - Asmjit is now embedded into PLUMED. In order to enable it, it is sufficient to configure with `--enable-asmjit`. See \ref Lepton "this page".
  - Fixed grids so as to decrease memory footprint of derivatives (see \issue{465}).
  - Added option `--idlp4` to \ref driver to read DLPOLY4 HISTORY files (see \issue{478}, thanks to Alin Marin Elena).
  - Added atom selectors using mdtraj/MDAnalysis/VMD syntax, see \ref MOLINFO and \issue{448}.
  - \ref EEFSOLV is now faster in scalar and also mpi/openmp parallel
  - New shortcuts are available for selecting protein atoms: `@sidechain-#`, `@back-#`
  - VIM syntax highlight is now case insensitive. Notice that autocompletion still only works with upper case commands.

- New contributed modules:
  - A new Maze module by Jakub Rydzewski
     - \ref MAZE_LOSS
     - \ref MAZE_MEMETIC_SAMPLING
     - \ref MAZE_RANDOM_ACCELERATION_MD
     - \ref MAZE_RANDOM_WALK
     - \ref MAZE_SIMULATED_ANNEALING
     - \ref MAZE_STEERED_MD
     - \ref MAZE_OPTIMIZER_BIAS
  - A new ANN module by Wei Chen and Andrew Ferguson
     - \ref ANN

- New patches:
  - added support for AMBER PMEMD 18 (contributed by Viktor Drobot, see \issue{486}).

- Changes in the VES module
  - new \ref VES_DELTA_F bias.
  - ves_md_linearexpansion now outputs one-dimensional free energy projections of the potential energy landscape. 

- Changes in the DRR module
  - The MAXFACTOR option now is tunable for each CV in multidimensional cases.
  - Output .zcount file (the same as .czar.count) for compatibility with newer abf_integrate.
  - The citation of DRR module has been updated.

- Changes in the ISDB module
  - in \ref METAINFERENCE we removed the MC_STRIDE keyword
  - in \ref METAINFERENCE the bias value (metainference score) now includes the Jeffrey's prior (values are different, but forces are equal)
  - components were previously named using _ but now they abide to the standard is -
  - removed ADDEXP keywords for \ref JCOUPLING \ref NOE \ref PRE \ref RDC
  - \ref METAINFERENCE performs more check on the input and restart files to ensure a consistent setup
  - \ref SAXS is slightly faster and scales better, removed BESSEL options

- Python module:
  - Removed compatibility with Python 2.
  - Added capability to read and write pandas dataset from PLUMED files (see \issue{496}).

Changes from version 2.5 which are relevant for developers:
  - Components documentation is now enforced
  - `readdir_r` is deprecated and is thus not used by default (can be enabled with `./configure --enable-readdir-r`).

## Version 2.6.1

For users:
- New patches:
  - added gromacs 2019.6 
  - added gromacs 2020.1 (experimental) 

For developers:
- Small fix to avoid unique global symbols (see \issue{549})

@page CHANGES-2-1 Version 2.1

Version 2.1.0 (September 15, 2014)
----------------------------

Version 2.1 contains several improvements with respect to 2.0. Users currently working with 2.0 
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files. In 2.1 we restored more features of 1.3
that were missing in 2.0, so users still working with 1.3 could opt for an upgrade.
A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Below you find a list of all the changes with respect to version 2.0.
Notice that version 2.1 includes already all the fixes in branch 2.0 up to 2.0.4.

Changes from version 2.0 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref COORDINATION now skips pairs of one atom with itself.
  - Labels of quantities calculated by \ref BIASVALUE have changed from <i>label</i>.bias.<i>argname</i> to <i>label</i>.<i>argname</i>_bias, which is more consistent with steered MD
  - Labels of quantities calculated by \ref ABMD have change from <i>label</i>.min_<i>argname</i> to <i>label</i>.<i>argname</i>_min, which is more consistent with steered MD
  - Labels of quantities calculated by \ref PIECEWISE have change from <i>label</i>.<i>argnumber</i> to <i>label</i>.<i>argname</i>_pfunc, which is more consistent with steered MD
  - For multicolvars components calculated with LESS_THAN and MORE_THAN keywords are now labelled lessthan and morethan. This change is necessary as the underscore
character now has a special usage in component names.
  - In \ref CONTACTMAP components are now labelled <i>label</i>.contact-\f$n\f$.
  - The command SPHERE has been replaced by \ref UWALLS.
- New configuration system based on autoconf (use ./configure from root directory).
  Optional packages are detected at compile time and correctly
  enabled or disabled. An internal version of LAPACK and BLAS will be used
  if these libraries are not installed.
- New actions:
  - \ref SPRINT topological collective variables.
  - CH3SHIFTS collective variable.
  - \ref POSITION collective variable.
  - \ref FIT_TO_TEMPLATE.
  - \ref COMMITTOR analysis.
  - \ref LOCAL_AVERAGE.
  - \ref NLINKS.
  - \ref DIHCOR.
  - \ref NOE.
  - \ref RDC. 
  - \ref CLASSICAL_MDS.
  - \ref XDISTANCES.
  - \ref YDISTANCES.
  - \ref ZDISTANCES.
  - \ref DUMPMULTICOLVAR.
  - Crystallization module, including \ref Q3, \ref LOCAL_Q3, \ref Q4, \ref Q6, \ref LOCAL_Q4, \ref LOCAL_Q6, \ref MOLECULES, \ref SIMPLECUBIC, \ref TETRAHEDRAL and \ref FCCUBIC.
  - \ref ENSEMBLE to perform Replica-Averaging on any collective variable.
- New features for existing actions:
  - \ref METAD : WALKERS_MPI flag (multiple walkers in a mpi-based multi-replica framework),
    ACCELERATION flag (calculate on the fly the Metadynamics acceleration factor),
    TAU option (alternative way to set Gaussian height in well-tempered metadynamics),
    GRID_SPACING (alternative to GRID_BIN to set grid spacing).
    Notice that now one can also omit GRID_BIN and GRID_SPACING when using
    fixed size Gaussian, and the grid spacing will be automatically set.
  - \ref DISTANCE : added SCALED_COMPONENTS
  - \ref COORDINATION : if a single group is provided, it avoids permuted atom indexes and runs
    at twice the speed.
  - \ref DUMPATOMS : PRECISION option to set number of digits in output file.
  - \ref GROUP : NDX_FILE and NDX_GROUP options to import atom lists from ndx (gromacs) files.
  - In many multicolvars, MIN and MAX options can be used.
  - \ref HISTOGRAM : GRID_SPACING (alternative to GRID_BIN to set grid spacing),
    FREE-ENERGY flags in addition to standard probability density,
    additional option for KERNEL=DISCRETE to accumulate standard histograms. 
  - \ref sum_hills : added options --spacing (alternative to --bin to set grid spacing)
    and --setmintozero to translate the minimum of the output files to zero.
  - \ref CONTACTMAP : parallelized and added weights.
- New features in MD patches (require re-patch):
  - New patch for Gromacs 5.0
  - Gromacs 4.6.X patch updated to 4.6.7
  - Gromacs 4.6.7 supports \ref COMMITTOR analysis; can be now be used to perform energy minimization;
     now passes temperature to PLUMED (this allows temperature to be omitted in some actions,
     namely \ref METAD and analysis actions).
  .
  Notice that if you use runtime binding it is not compulsory to re-patch,
  and that all combinations should work correctly
  (new/old PLUMED with re-patched/non-re-patched MD code).
- Other new features:
  - \ref driver can now read trajectories in many formats using VMD molfile plugin
    (requires VMD plugins to be compiled and installed). In case VMD plugins are not installed,
    the configuration system falls back to an internal version which implements a minimal list
    of plugins (gromacs and dcd) (kindly provided by T. Giorgino).
  - \ref switchingfunction : added STRETCH flag.
  - Negative strides in atom ranges (e.g. ATOMS=10-1:-3 is expanded to ATOMS=10,7,4,1).
  - \ref COORDINATION and \ref DHENERGY with NLIST now work correctly in replica exchange simulations.
  - Multicolvars with neighbor lists now work correctly in replica exchange simulations.
  - Improved multicolvar neighbor lists.
- Optimization:
  - Root-mean-square deviations with align weights different from displace weights
    are now considerably faster. This will affect \ref RMSD calculations plus
    other variables based on RMSD.
  - \ref WHOLEMOLECULES is slightly faster.
  - \ref COORDINATION is slightly faster when NN and MM are even and D_0=0.
  - Atom scattering with domain decomposition is slightly faster.
  - Link cells are now exploited in some multicolvars.
  - Derivatives are not calculated unless they are specifically required, because for instance you are adding
    a bias.
- Documentation:
  - All tutorial material from the recent plumed meeting in Belfast is now in the manual
  - Improvements to documentation, including lists of quantities that are output by each action that can be referenced 
  - Manual has been re-organized following suggestions received at the plumed meeting.
  - An experimental PDF version of the manual is now provided (a link can be found in the documentation homepage).

Changes from version 2.0 which are relevant for developers:
- Added regtests for plumed as a library (e.g. basic/rt-make-0). plumed command has an additional
  flag (--is-installed) to probe if running from a compilation directory or from a fully installed copy
  (this is needed for regtests to work properly).
- Improved class Communicator. Many operations can now be done directly on Vectors, Tensors, std::vector and PLMD::Matrix.
- Modified class RMSD.
- Patches for GPL codes (Quantum Espresso and Gromacs) now also include
  original code so as to simplify their modification.
- Fixed dependencies among actions such that it is now possible (and reliable)
  to use MPI calls inside Action::prepare()
- colvar/CoordinationBase.cpp has been changed to make it faster. If you devised a class which inherits from here,
  consider that CoordinationBase::pairing now needs _squared_ distance instead of distance
- It is possible to run "make install" from sub-directories (e.g. from src/colvar)
- There is a small script which disables/enables all optional modules (make mod-light/mod-heavy/mod-reset)
- Added "-q" option to plumed patch
- You can now create new metrics to measure distances from a reference configurations. If you do so such
  metrics can then be used in paths straightforwardly
- You can now use multicolvars in tandem with manyrestraints in order to add a large numbers of restraints.
- Can now do multicolvar like things in which each colvar is a vector rather than a scalar.
- Updated script that generated header files so that they properly show years. Notice that the script
  should new be run from within a git repository

This list is likely incomplete, if you are developing in PLUMED you are encouraged to follow changes on github.

Version 2.1.1 (December 15, 2014)
----------------------------------------------

This release includes all the fixes available in branch 2.0 until 2.0.5.

For users:
- New patch for AMBER 14 (sander module only). This patch should be compatible
  with any PLUMED 2 version (including 2.0). It includes most PLUMED features
  with the notable exception of multi-replica framework.
- Changed definition in arbitrary phase of eigenvectors. This will change the result of some
  analysis method where the phase does matter (e.g. \ref CLASSICAL_MDS) and make
  some regression test better reproducible.
- Fixed a portability issue in BG/P where gettimeofday is not implemented.
  Notice that this fix implies that one should execute again ./configure to have
  plumed timing working correctly.
- CS2Backbone: fixed a bug that resulted in only a fraction of the chemical shifts being printed with WRITE_CS and 
  parallel simulations (requires to get the last almost updated from SVN)
- NOE: fixed a bug in the replica-averaging
- Fixed a linking issue with ALMOST, where bz2 was always used to link ALMOST to PLUMED even if it is not compulsory 
  to build ALMOST.
- Fixed a wrong include in the GMX5 patch.
- \ref FUNCPATHMSD can now be used together with \ref CONTACTMAP to define pathways in contact-map space
- Configuration is more verbose, a warning is given if a default option cannot be enabled and an error is given if 
  an option explicitly enabled cannot be enabled.
- Compilation is less verbose (use "make VERBOSE=1" to have old behavior)
- Small fixes in documentation.

For developers:
- Tests are now performed at every single push on travis-ci.org
- Manual is built and pushed to the online server from travis-ci.org (see developer doc)
- Fixes in developer doc.

Version 2.1.2 (Mar 16, 2015)
----------------------------------------------

For users:
- Added two new short tutorials to the manual ( \ref cambridge and \ref munster ).
- Fixed a severe bug on \ref DRMSD - cutoff values were ignored by PLUMED.
  Notice that this bug was introduced in 2.1.0, so that it should not affect the 2.0.x series.
- Fixed a bug affecting LAMMPS patch used with a single processor. Notice that
  the fix is inside PLUMED, thus it does not necessarily requires re-patching.
- Sander patch now works with multiple replica (no replica exchange yet). It also contains
  some fix from J. Swails.
- GMX5 patch was not working for bias-exchange like cases
- Patching system now checks for the availability of shared/static/runtime version of plumed before
  patching
- Configure now check better if compiler flag are accepted by the compiler. This makes
  configure on bluegene more robust.
- Sourceme.sh now sets proper library path in linux also.


Version 2.1.3 (June 30, 2015)
----------------------------------------------

For users:
- Fixed bug in \ref ENSEMBLE derivatives when more than 1 argument was provided
- Fixed bug in \ref GHOST : virial is now computed correctly.
- Fixed a serious bug in virial communicated from plumed to gromacs, for both gromacs versions 4.6 and 5.0.
  See \issue{132}.
  This fix requires gromacs to be re-patched and could be very important if you run biased simulations in the NPT ensemble.
- Fixed a bug in the virial computed with \ref FIT_TO_TEMPLATE when the reference pdb had center non located at the origin.
- Fixed a bug in the the forces computed with \ref FIT_TO_TEMPLATE when used in combination with \ref COM, \ref CENTER, or \ref GHOST
- Fixed a bug that could lead plumed to be stuck with domain decomposition in some extreme case (one domain with all atoms, other domains empty).
- Fixed a bug when \ref COMBINE or \ref MATHEVAL are used with PERIODIC keyword. Now when PERIODIC keyword is used the result
  of the calculation is brought within the periodicity domain. See \issue{139}.
- Fixed a bug related to \ref RANDOM_EXCHANGES followed by \ref INCLUDE
- Fixed bug in derivatives of histogram bead with triangular kernels
- Updated gromacs patch 4.5.5 to 4.5.7
- Updated internal molfile plugins to VMD 1.9.2.
- Included crd and crdbox formats to internal molfile.
- Added --natoms to \ref driver . This is required to read coordinate
  files with VMD plugins when number of atoms is not present (e.g. amber
  crd files)
- Added the checks in the driver to detect cases where molinfo does not provide box information
  (e.g. pdb).
- Added support for readdir_r when available, which makes opening files thread safe.
- CFLAGS now include -fPIC by default
- Added a warning when using \ref METAD without grids with a large number of hills.
- Fixes in user documentation.

For developers:
- Allow external VMD plugins to be detected with --has-external-molfile. This
  is required to enable some regtest with amber files.
- Added --dump-full-virial to \ref driver
- Allow definition of variables where some of the components have derivatives and some haven't (\issue{131}).
- Improved travis tests with more debug options.
- Improved some regtest to check out-of-diagonal virial components
- Improved make cppcheck options.
- Fixes in developer documentation.

Version 2.1.4 (Oct 13, 2015)
-----------------------------

For users:
- Fixed NAMD patch. Masses and charges were not passed correctly, thus resulting in wrong
  \ref COM or \ref CENTER with MASS.
  This fix required re-patching NAMD.
  Notice that this bug was present also in v2.0 but in a different form.
  More information here (\issue{162}), including a workaround that allows masses to be fixed
  without re-patching.
- When installing with PLUMED_LIBSUFFIX an underscore is used as separator instead of a dash.
  E.g. `make install PLUMED_LIBSUFFIX=2.1` will result in an executable named `plumed_v2.1`.
  This fix a potential problem (see \ref Installation).
- Fixed erroneously reported message about MPI at the end of ./configure.
- Changed warning message about undocumented components.
- PLUMED now says in the log file if it was compiled from a dirty git repository.
- Fixed a problem leading to rare random crashes when using \ref METAD with WALKERS_MPI and multiple
  processors per replica.
- Small change in numerical accuracy of lattice reduction. Should be more
  robust when running with highly optimizing compilers.
- Fixed a bug in normalization of kernel functions.  This affects \ref HISTOGRAM
  If these actions were used with previous versions of the code care should be taken when analyzing the 
  results.
- Fixed a bug in derivatives of kernel functions with non-diagonal covariance matrices. This affects the 
  derivatives output by \ref sum_hills

Version 2.1.5 (Jan 18, 2016)
---------------------------------------------

\plumednotmaintained

For users:
- PLUMED now reports an error when using \ref HISTOGRAM with FREE-ENERGY without USE_ALL_DATA. See \issue{175}
- Fixed a bug in configure together with --enable-almost. The check for lbz2 library was not working properly.

@page CHANGES-2-5 Version 2.5

## Version 2.5 (Dec 19, 2018)

This page contains changes that will end up in 2.5

Changes from version 2.4 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref RMSD, \ref MULTI-RMSD, \ref PATHMSD, \ref PROPERTYMAP, \ref PCAVARS, \ref PCARMSD, \ref FIT_TO_TEMPLATE,
    \ref DIPOLE, \ref ALPHARMSD, \ref ANTIBETARMSD, and \ref PARABETARMSD now automatically make molecules whole.
    In case you do not want them to do it, use NOPBC flag,
  - There is some subtle change in the installation layout (see below). There should be no visible effect, however it is now compulsory
    to set correctly the `LD_LIBRARY_PATH` variable for the linux executable to work correctly. The procedure has been tested well on OSX and Linux,
    but could give problems on other platform. Please report possible problems on the mailing list.
  - \ref driver now stops correctly when using \ref COMMITTOR. If you want to continue the analysis, use the `NOSTOP` flag in \ref COMMITTOR.
  - \ref METAD the calculation of the reweighting factor is now activated by CALC_RCT instead of REWEIGHTING_NGRID and REWEIGHTING_NHILLS, the frequency of update can be set 
    by RCT_USTRIDE, the default value is 1 and should be OK for most of the cases
  - Fixed sign in Cartesian components of \ref PUCKERING with 6 membered rings (thanks to Carol Simoes and Javi Iglesias).

- New actions:
  - \ref COLLECT_FRAMES
  - \ref EUCLIDEAN_DISSIMILARITIES
  - \ref HBPAMM_MATRIX
  - \ref HBPAMM_SH
  - \ref LANDMARK_SELECT_FPS
  - \ref LANDMARK_SELECT_RANDOM
  - \ref LANDMARK_SELECT_STAGED
  - \ref LANDMARK_SELECT_STRIDE
  - \ref OUTPUT_ANALYSIS_DATA_TO_COLVAR
  - \ref OUTPUT_ANALYSIS_DATA_TO_PDB
  - \ref OUTPUT_PCA_PROJECTION
  - \ref PAMM
  - \ref PLUMED
  - \ref PRINT_DISSIMILARITY_MATRIX
  - \ref PROJECT_ALL_ANALYSIS_DATA
  - \ref READ_DISSIMILARITY_MATRIX
  - \ref RESELECT_LANDMARKS
  - \ref REWEIGHT_WHAM
  - \ref SKETCHMAP_CONJGRAD
  - \ref SKETCHMAP_POINTWISE
  - \ref SKETCHMAP_READ
  - \ref SKETCHMAP_SMACOF
  - \ref SKETCH_MAP
  - \ref SMACOF_MDS
  - \ref WHAM_HISTOGRAM
  - \ref WHAM_WEIGHTS

- New command line tools:
  - \ref completion (used to generate command line completion scripts).
  - \ref pdbrenumber (see \issue{371}).

- New modules:
  - A new PIV module has been included, contributed by Silvio Pipolo and Fabio Pietrucci.
    This module implements the following collective variable:
    - \ref PIV
  - A new LOGMFD module has been included, contributed by Tetsuya Morishita.
    This module implements the following bias:
    - \ref LOGMFD

- Changes in the ISDB module
  - \ref CS2BACKBONE is now mpi parallelized in particular with DOSCORE and CAMSHIFT
  - \ref SAXS has an additional implementation based on Bessel functions that can be faster for large systems (new keyword BESSEL)
  - \ref SAXS keyword SCEXP has been renamed into SCALEINT
  - \ref SAXS includes the MARTINI bead structure factors for Proteins and Nucleic Acids
  - \ref SAXS includes a GPU implementation based on ArrayFire (need to be linked at compile time) that can be activated with GPU
  - \ref METAINFERENCE and all related methods has a new keyword REGRES_ZERO to scale data using a linear scale fit
  - \ref CALIBER new bias to perform Maximum Caliber replica-averaged restrained simulations 

- Changes in the eABF/DRR module (contributed by Haochuan Chen and Haohao Fu):
  - \ref DRR now supports the extended generalized ABF(egABF) method.
  - \ref DRR accepts different GRID options for CVs and extended variables.
  - The MAXFACTOR option is added in \ref DRR to control the factor of biasing force.
  - \ref drr_tool can calculate the divergence of gradients now. (Maybe useful for future pABF)
  - Fixed conflicts of output files in multiple replicas.

- Changes in the EDS module:
  - \ref EDS implements Levenberg-Marquardt optimization in addition to previous gradient descent. 
  - \ref EDS no longer automatically increases prefactor for bias parameter updates. This results in more stable optimization for the cases tested.
  - \ref EDS now has a larger default RANGE parameter to go with these other changes.

- Other changes:
  - \ref METAD there is a new FLYING_GAUSSIAN keyword to activate the flying gaussian methods by Spiwok (contributed by Spiwok and Hozzova)
  - \ref EXTERNAL can now SCALE the input grid. This allows for more flexibility without modifying the grid file.
  - \ref ALPHABETA can now combine dihedral angles with different coefficients
  - \ref INCLUDE can now be used also before setup actions.
  - \ref CENTER can now be computed using trigonometric functions (PHASES) to simplify its calculation with periodic boundary conditions.
  - Libmatheval is not used anymore. \ref MATHEVAL (and \ref CUSTOM) are still available
    but employ an internal implementation of the lepton library.
    Functions available in libmatheval and absent in the original lepton library have been added so as to have backward compatibility.
    `atan2(y,x)` function has also been added.
    Notice that MATHEVAL (and CUSTOM) \ref switchingfunction "switching functions"
    using the lepton library have been further optimized with respect to PLUMED 2.4.
    Finally, notice that it is possible to use asmjit to optimize performance (see \ref Lepton).
  - Implemented bash autocompletion, see \ref BashAutocompletion.
  - \ref MOLINFO now allows selecting atoms from chains with a numeric ID (see \issue{320}).
  - Removed the patch for GMX 5.1.4
  - LAMMPS patch has been finally removed. Notice that LAMMPS has native support for PLUMED now.
  - AMBER patch has been finally removed. Notice that AMBER (sander module) has native support for PLUMED starting from version 15.
  - \ref RMSD calculation has been optimized. This should positively affect the performances of CVs where
     many RMSD values are computed on small groups of atoms, such as secondary structure variables.
  - In \ref METAD, when using a bias factor equal to one (no bias) the `rct` component is set to zero rather than to one.
  - New shortcuts are available for selecting atoms: `@allatoms` and `@mdatoms` (see \ref atomSpecs).
  - When using \ref MOLINFO, also the following shortcuts are available for selecting atoms: `@nucleic`, `@protein`, `@water`, `@ions`, `@hydrogens`, `@nonhydrogens`.
  - When using \ref MOLINFO, individual atoms can be chosen also from water molecules (e.g. `@OW-100`).
  - Additional switching function COSINUS contributed by Michael King
  - added API to set the number of used openMP threads from the linked code, updated gromacs 2018.3 patch to use it

Changes from version 2.4 which are relevant for developers:
- Code has been cleanup up replacing a number of pointers with `std::unique_ptr`. All `delete` statements
  in the core parts of the code have been eliminated.
- Exceptions cannot be disabled (`--disable-cxx-exceptions` option has been removed from `./configure`).
- Every exception thrown in PLUMED now also writes its message on PLUMED log.
- Runtime loader in `Plumed.c` now works also when linked without `-rdynamic` (that is, 
  its names are not exported). Notice that all the combinations are expected to
  work, that is: `Plumed.c` from <=2.4 or >=2.5 combined with libplumedKernel
  from <=2.4 or >=2.5. In order to achieve this the following changes are implemented:
  - libplumedKernel does not depend anymore on `Plumed.c`. This allows loading it even
    in cases where names in the loader are not visible. The relevant function needed
    to be compatible with `Plumed.c` <=2.4 are found using `dlsym`.
  - `Plumed.c` does not need anymore libplumedKernel to register itself, but rather
    searches the relevant functions using `dlsym`. In addition, if it is not able to
    load `libplumedKernel` since the latter is <=2.4 and needs `Plumed.c` to be visible,
    it just uses as a fallback `libplumed`, which should load properly.
- In addition to the capability mentioned above, the MD-code interface has been significantly
  improved and allows for:
  - Translation of exception (allowing to mix PLUMED and an MD-code linked against a different C++ library).
  - Possibility to choose the path to the PLUMED kernel while instantiating a Plumed object.
  See the developer documentation for more information.
- The installation layout of shared libraries has been modified. In particular,
  both `libplumed.so` and `plumed` links to `libplumedKernel.so`.
  This reduces considerably the size of the installed package. In addition, it allows
  using two-level namespace on OSX. Notice that this implies that on Linux one should
  always set the `LD_LIBRARY_PATH` flag to have a working executable.
- A smaller number of header files is installed. In particular, all the files that were historically generated in subdirectories
  (such as `plumed/core/tools/Vector.h', just including `plumed/tools/Vector.h`) are not installed and the related include
  statements are fixed. This makes the installed package smaller.
- List of preferred compilers (used when `CXX` or `CC` are not set) has been changed. On OSX, `./configure` will try `clang++/clang` as first choices.
- Added `--enable-static-archive` to `./configure` to build a `libplumed.a` static library (yes by default).
- Stop setting `DYLD_LIBRARY_PATH` in `sourceme.sh` and in modulefile. Notice that as of PLUMED v2.3.3
  it should not be needed.
- Coverage scan is not anymore contained in developer manual. It can be found in a separate repository
  `github.com/coverage-branch` (see \issue{348}). In addition, coverage for third-party libraries included in PLUMED
  is reported as well.
- It is not possible anymore to use `make install prefix=/path`. Prefix can only be changed during `./configure` (see \issue{332}).
- Exception class has been rewritten to allow more extensive messages. Now also function name is shown.
- On linux, library is linked with `-Bsymbolic`.
- When launching `plumed`, flags `--no-mpi` and `--mpi` can appear multiple times. The last appearance is the effective one.
- Internal BLAS and LAPACK libraries updated to gromacs 2018.
- Choosing `./configure --prefix=$PWD` does not lead anymore to deletion of all header files.
- A copy of `plumed-runtime` is installed in `prefix/lib/plumed` and can be used for testing.
- Absolute/relative soname/install_name can be configured on linux/OSX. This feature is only
  for testing, the default choice is the typical one used on the respective operating system.
- On OSX, `plumed` and `libplumed.dylib` will find `libplumedKernel.dylib` using `@loader_path`.
- Using CXX compiler to link the main program.
- plumed can be compiled with ArrayFire to enable for gpu code. \ref SAXS collective variable is available as part of the isdb module to provide an example of a gpu implementation for a CV


## Version 2.5.1 (Apr 1, 2019)

For users:
- in \ref SAXS the keyword ADDEXP is removed. Furthemore, SAXS intensities are automatically normalised for I(0)=1, in case experimental data are provided, the intensity is rescaled with the intensity of the lowest q provided. As a consequence SCALEINT is only needed for additional adjustments.
- gromacs patch updated to gromacs 2018.5
- Fixed a bug in gromacs patch that was resulting in incorrect number of threads (0) set when not explicitly using `-ntomp` on the 
  command line or setting `OMP_NUM_THREADS` (see \issue{446}). To apply this fix you need to re-patch gromacs.
  Notice that setting the number of threads to zero might lead to inconsistent results when using secondary structure variables
  or other multicolvars.
- Fixed PLUMED so that when zero threads are selected from gromacs (see previous fix) the number of used threads is set to 1.
  This fix allows to use a GROMACS executable patched with PLUMED 2.5.0 and linked at runtime with PLUMED 2.5.1 without introducing
  errors. However, re-patching is preferred since it selectes the correct number of threads.
- Python wrappers:
  - Fixed building of python interface on MacOS Mojave (see \issue{445}, thanks to Omar Valsson).
  - Numpy is not required anymore at build time (though it is required at runtime for our tests).
  - Raw python arrays can be passed as an alternative to Numpy ndarrays.

## Version 2.5.2 (Jul 19, 2019)

For users:
- New shortcuts are available for selecting protein atoms: `@chi2-#`, `@chi3-#`,`@chi4-#` and `@chi5-#`
- Fixed performance of \ref CUSTOM when having zero derivatives with respect to some arguments.
- New --parse-only option in \ref driver to check the validity of a plumed input file
- New patch for GROMACS 2019.2
- Module VES: Fixed performance of \ref BF_CUSTOM for basis functions with linear terms (e.g. having zero derivatives). 
- Python wrappers:
  - Python module is now always named `plumed` irrespectively of program prefix and suffix. Notice 
    that python module is installed inside the `lib/program_name` directory and thus it is not necessary to
    use `program_name` in order to install multiple modules side by side.
  - Python module can be compiled without compiling PLUMED first.
  - `Plumed` object can be explicitly finalized using `finalize()`. Can be used to make sure all files are closed,
    but it is not necessary if the `Plumed` object gets correctly collected by Python.
  - `Plumed` object can be used in context managers (e.g. `with plumed.Plumed() as p:`).
- Precompiled binaries are available on Anaconda cloud on the [conda-forge channel](https://anaconda.org/conda-forge/plumed).

## Version 2.5.3 (Oct 11, 2019)

For users:
- Fixed a bug with \ref CONVERT_TO_FES and periodic variables, see \issue{441}
- Fixed a bug with \ref FOURIER_TRANSFORM 
- Updated patch for GROMACS 2019.4
- Updated patch for GROMACS 2018.8
- Python module:
  - Fixed building with clang-8.
  - Set `language_level` for cython to the actually used language level.
  - Force using cython when compiling from source. Still using the pre-generated cpp file
    when installing from PyPI, to avoid cython dependency.
  - Using python 2 to create the cpp file uploaded on PyPI (this will change to python 3 in 2.6, see \issue{502}).
- Module VES: Fixed a bug in updating of bias potential in \ref VES_LINEAR_EXPANSION that is present for certain integrators that call the calculation of the bias multiple times (see [here](https://groups.google.com/d/msg/plumed-users/kPZu_tNZtgk/LrkS0EqrCQAJ)) and replica exchange.

## Version 2.5.4 (Jan 27, 2020)

For users:
- Includes all fixes up to 2.4.7

## Version 2.5.5

For developers:
- Small fix to avoid unique global symbols (see \issue{549})

@page CHANGES-2-2 Version 2.2

Version 2.2 (Oct 13, 2015)
----------------------------

Version 2.2 contains several improvements with respect to 2.1. Users currently working with 2.1
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files. In 2.2 we restored more features of 1.3
that were missing in 2.1, so users still working with 1.3 could opt for an upgrade.
A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Below you find a list of all the changes with respect to version 2.1.
Notice that version 2.2 includes already all the fixes in branch 2.1 up to 2.1.4 indicated in \ref CHANGES-2-1 .

Changes from version 2.1 which are relevant for users:
- Changes leading to incompatible behavior:
  - Labels of quantities calculates by \ref SPRINT have changed from <i>label</i>.coord_<i>num</i> to <i>label</i>.coord-<i>num</i>
  - \ref METAD with WALKERS_MPI now writes a single hills file, without suffixes
  - removed the ./configure.sh script of v2.0.x, now plumed can only be configured using autotools (./configure)
  - \ref COM, \ref CENTER, and \ref GYRATION now automatically make molecules whole. In case you do not want them to do it, use NOPBC flag,
    which recovers plumed 2.1 behavior
  - Some MD code could now automatically trigger restart (e.g. gromacs when starting from cpt files). This can be overwritten using
    \ref RESTART NO.
  - Replica suffixes are now added by PLUMED *before* extension (e.g. use plumed.0.dat instead of plumed.dat.0)
  - When using \ref switchingfunction the STRETCH keyword is now implicit. NOSTRETCH is available to enforce the old behavior.
- Module activation can now be controlled during configure with `--enable-modules` option.
- Almost complete refactoring of installation procedure. Now
  DESTDIR and other standard autoconf directories (e.g. bindir) are completely supported.
  Additionally, everything should work properly also when directory names include spaces (\issue{157}).
  Finally, compiler is not invoked on install unless path are explicitly changed (\issue{107}).
- Related to installation refactoring, upon install a previously installed PLUMED is not removed.
  This is to avoid data loss if prefix variable is not properly set
- Several changes have been made in the Makefile.conf that makes it not compatible with those
  packaged with plumed 2.0/2.1. Please use ./configure to generate a new configuration file.
- Added partial OpenMP parallelization, see \ref Openmp
- Added multiple time step integration for bias potentials, see \ref MTS
- Link cells are now used in all multicolvars that involve \ref switchingfunction.  The link cell cutoff is
  set equal to 2.*\f$d_{\textrm{max}}\f$.  Where \f$d_{\textrm{max}}\f$ is the (user-specified) point at which
  the switching function goes to zero. Users should always set this parameter when using a switching function
  in order to achieve optimal performance.
- DHENERGY option is no longer possible within \ref DISTANCES.  You can still calculate the DHENERGY colvar by using \ref DHENERGY
- Reweighting in the manner described in \cite Tiwary_jp504920s is now possible using a combination of the \ref METAD and \ref HISTOGRAM actions.  The relevant keywords in \ref METAD are REWEIGHTING_NGRID and REWEIGHTING_NHILLS.  The \f$c(t)\f$ and the appropriate weight to apply to the configurations are given by the values labeled rct and rbias. 
- News in configure and install:
  - ./configure now allows external BLAS to be used with internal LAPACK. This is done automatically if only BLAS are available,
    and can be enforced with --disable-external-lapack.
  - ./configure supports --program-prefix, --program-suffix, and --program-transform-name.
  - make install supports DESTDIR and prefix.
  - Environment variables PLUMED_LIBSUFFIX and PLUMED_PREFIX are deprecated and will be removed in a later version.
- New actions
  - \ref DUMPMASSCHARGE to dump a file with mass and charges during MD.
  - \ref EFFECTIVE_ENERGY_DRIFT to check that plumed forces are not screwing the MD integrator.
  - \ref EXTENDED_LAGRANGIAN : in combination with  \ref METAD it implements metadynamics with Extended Lagrangian; standalone it implements TAMD/dAFED.
  - \ref DFSCLUSTERING calculate the size of clusters 
  - \ref DUMPMULTICOLVAR print out a multicolvar
  - \ref MFILTER_LESS filter multicolvar by the value of the colvar
  - \ref MFILTER_MORE 
  - \ref MFILTER_BETWEEN
  - \ref PCARMSD PCA collective variables using OPTIMAL rmsd measure
  - \ref PCAVARS PCA collective variables using any one of the measures in reference
  - \ref GRADIENT can be used to calculate the gradient of a quantity.  Used to drive nucleation
  - \ref CAVITY
  - \ref PUCKERING implemented for 5-membered rings (thanks to Alejandro Gil-Ley).
  - \ref WRAPAROUND to fix periodic boundary conditions.
- New features for existing actions:
  - Keywords UPDATE_FROM and UPDATE_UNTIL to limit update step in a defined time window, available only for actions where it would be useful.
  - Keyword UNNORMALIZED for \ref HISTOGRAM.
  - Possibility to use Tiwary-Parrinello reweighting for \ref METAD
  - Keywords for \ref GROUP (REMOVE, SORT, UNIQUE) to allow more flexible editing of groups.
  - \ref DUMPATOMS now supports dumping xtc and trr files (requires xdrfile library).
  - \ref driver can now read xtc and trr files also with xdrfile library.
  - \ref driver accepts a --mc flag to read charges and masses from a file produced during
    molecular dynamics with \ref DUMPMASSCHARGE
  - Possibility to enable or disable \ref RESTART on a per action basis, available only for actions where it would be useful.
  - \ref MOLINFO now supports many more special names for rna and dna (thanks to Alejandro Gil-Ley).
  - VMEAN and VSUM allow one to calculate the sum of a set of vectors calculated by VectorMultiColvar.  Note these
  can also be used in tandem with \ref AROUND or \ref MFILTER_MORE to calculate the average vector within a particular
  part of the cell or the average vector among those that have a magnitude greater than some tolerance
  - New way of calculating the minimum value in multicolvars (ALT_MIN). This is less susceptible to overflow for certain 
    values of \f$\beta\f$.  
  - New keywords for calculating the LOWEST and HIGHEST colvar calculated by a multicolvar
  - Added components to \ref DIPOLE (\issue{160}).
- Other changes:
  - File reader now supports dos newlines as well as files with no endline at the end.

For developers:

- In order to be able to use openMP parallelism within multicolvar, secondarystructure, manyrestraints and crystallisation
we had to make some substantial changes to the code that underlies these routines that is contained within vesselbase. In 
particular we needed to get rid of the derivatives and buffer private variables in the class ActionWithVessel.  As a consequence
the derivatives calculated in the various performTask methods are stored in an object of type MultiValue.  Within multicolvar
this is contained within an object of type AtomValuePack, which stores information on the atom indices.  If you have implemented
a new multicolvar it should be relatively straightforward to translate them so they can exploit this new version of the code.  Look 
at what has been done to the other multicolvars in there for guidance.  Sorry for any inconvenience caused.
- Changed the logic of several PLUMED ifdef macros so as to make them consistent.
  Now every feature based on external libraries is identified by a __PLUMED_HAS_* macro.

Version 2.2.1 (Jan 18, 2016)
---------------------------------------------

For users:
- \ref PBMETAD implement the new Parallel Bias Metadynamics flavor of the Metadynamics sampling method.
- PLUMED now reports an error when using \ref HISTOGRAM with UNNORMALIZED without USE_ALL_DATA. See \issue{175}
- Fixed a bug in configure together with --enable-almost. The check for lbz2 library was not working properly.
- Fixed a bug in install procedure that was introducing an error in linking with CP2K.
- Fixed a bug that sometimes was preventing the printing of a useful error message.

For developers:
- Vector and Tensor now support direct output with `<<`.
- Added some missing matmul operation Vector and Tensor.
- ./configure is automatically relaunched when changing ./configure or Makefile.conf. This makes it more robust
  to switch between branches.

Version 2.2.2 (Apr 13, 2016)
----------------------------------------------

For users:
- \ref MOLINFO for RNA accepts more residue names, see \issue{180}.
- added two mpi barries (one was missing in PBMetaD for multiple walkers) to help synchronized initialisation
- Fixed a bug in internal stopwatches that was making \ref DEBUG logRequestedAtoms not working
- Some multicolvars (including \ref BRIDGE, \ref ANGLES, and \ref INPLANEDISTANCES) now crashes if one
  asks for too many atoms, see \issue{185}.
- Optimisations (activation of the dependencies, secondary structures, DRMSD)
- Fixed a performance regression with RMSD=OPTIMAL-FAST
- Fixed a bug in the normalization of kernel functions (relevant for \ref HISTOGRAM).
- Fixed a regression introduced in v2.2 that was making \ref METAD with non-MPI multiple walkers crash
  if reading frequently. See \issue{190}
- Updated patch for gromacs 5.x. Patches for gromacs 5.0 and 5.1 have been fixed so as to allow
  patching in runtime mode.
- Possibility to control manual generation (including pdf) from ./configure. Pdf manual is now off
  by default. Notice that on travis CI it is still generated.

For developers:
- Fixed a bug in the interpretation of cmd strings. Namely, an erroneous string was not triggering an error.
  This is harmless for MD codes properly patched, but could have introduced problems in MD codes with typoes
  in cmd strings.
- ./configure is not automatically relaunched anymore when doing `make clean`.

Version 2.2.3 (Jun 30, 2016)
----------------------------------------------

For users:
- Updated patches for gromacs 5.1.x and 5.0.x to fix a problem when plumed was trying to write to an already
  closed gromacs log file.
- When looking for a value outside the GRID now the error include the name of the responsible 
  collective variable
- Numerical check in LatticeReduction made less picky. This should solve some of the internal errors reported
  by `LatticeReduction.cpp` when using aggressive compilers.
- Files are now flushed at the correct step. Before this fix, they were flushed at the step before the requested one
  (e.g. with \ref FLUSH STRIDE=100 at step 99, 199, etc).
- In \ref METAD, INTERVAL with periodic variables now report an error.
- \ref LOAD now works also when plumed is installed with a suffix.
- Added `--md-root` option to `plumed patch` which allows it to be run from a directory different from the one
  where the md code is located.
- Wham script in \ref munster tutorial now writes weights in scientific notation.

For developers:
- `./configure` checks if dependencies can be generated. If not, they are disabled.
- Added --disable-dependency-tracking to ./configure
- Added a make target `all_plus_doc` that builds both code and docs.
- Added possibility to set a default location for plumed library in runtime binding.
  If the plumed wrapped is compiled with `-D__PLUMED_DEFAULT_KERNEL=/path/libplumedKernel.so`,
  then if the env var PLUMED_KERNEL is undefined or empty PLUMED will look in the path at compile time.
- Tentative port files are now available at [this link](http://github.com/plumed/ports). 
  They can be used to install PLUMED using MacPorts.

Version 2.2.4 (Dec 12, 2016)
-------------

For users:
- Fix a bug in \ref PBMETAD when biasing periodic and not periodic collective variables at the same time 
- GSL library is now treated by `./configure` in the same way as other libraries, that is `-lgsl -lgslcblas` are only
  added if necessary.
- Fix a bug in \ref METAD when using INTERVAL and ADAPTIVE gaussians at the same time
- Updated gromacs patch for 5.1.x to 5.1.4
- Fix a performance regression in the calculate loop where derivatives and forces were set to zero even if an action
  was not active, this is relevant for postprocessing and for the on-the-fly analysis
- Torsion calculation has been made slightly faster and improved so as to provide correct
  derivatives even for special angles (e.g. +pi/2 and -pi/2).

For developers:
- Macports portile is now tested on travis at every plumed push.

Version 2.2.5 (Mar 31, 2017)
-------------

\plumednotmaintained

For users:
- Fixed a problem with large step numbers in driver (see \issue{209}).
- Fixed a problem leading to crashes when using switching functions without cutoff with some compiler (see \issue{210}).
- Fixed a bug when using \ref FIT_TO_TEMPLATE and domain decomposition (see \issue{214}).
- Added an automatic flush of HILLS files when using \ref METAD with file-based multiple walkers.
- Root dir is logged to allow easier debugging of problems.

@page CHANGES-2-3 Version 2.3

## Version 2.3 (Dec 12, 2016)

Version 2.3 contains several improvements with respect to 2.2. Users currently working with 2.2
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files.

Below you find a list of all the changes with respect to version 2.2.
Notice that version 2.3 includes already all the fixes in branch 2.2 up to 2.2.3 indicated in \ref CHANGES-2-2 .

Changes from version 2.2 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref COMMITTOR can now be used to define multiple basins, but the syntax has been changed
  - Syntax for \ref SPRINT and \ref DFSCLUSTERING has changed.
    We have separated the Actions that calculate the contact matrix from these actions.  These actions thus now take a contact
    matrix as input.  This means that we these actions can be used with contact matrices that measures whether or not a pair of atoms
    are hydrogen bonded.  For more details on this see \ref contactmatrix.  For clustering the output can now be passed to the actions
    \ref CLUSTER_PROPERTIES, \ref CLUSTER_DIAMETER, \ref CLUSTER_NATOMS, \ref OUTPUT_CLUSTER and \ref CLUSTER_DISTRIBUTION.  These
    provide various different kinds of information about the connected components found by clustering 
  - In \ref driver masses and charges are set by default to NaN.
    This makes it less likely to do mistakes trying to compute centers of mass or electrostatic-dependent variables
    when masses or charges were not set. To compute these variables from the driver you are now forced to use
    `--pdb` or `--mc`.
  - In rational switching functions, by default MM is twice NN. This is valid both in \ref switchingfunction with expanded
    syntax and when specifying MM on e.g. \ref COORDINATION
  - Patch script `plumed patch` now patches by default with `--shared`. This should make the procedure more robust (see \issue{186}).
  - Faster \ref GYRATION but new default behavior is not mass weighted
  - When using \ref HISTOGRAM you now output the accumulated grid using \ref DUMPGRID or \ref DUMPCUBE to get the free energy you use
    the method \ref CONVERT_TO_FES.  These changes allow one to use grids calculated within PLUMED in a work flow of tasks similarly to 
    the way that you can currently use Values.
  - The way that reweighting is performed is now different.  There are three separate actions \ref REWEIGHT_BIAS, \ref REWEIGHT_TEMP and
    \ref REWEIGHT_METAD.  These actions calculate the quantities that were calculated using the keywords REWEIGHT_BIAS and REWEIGHT_TEMP that
    used to appear in the old HISTOGRAM method.  Now those these methods can be used in any methods that calculate ensemble averages for
    example \ref HISTOGRAM and \ref AVERAGE
  - Manual is now build with locally compiled plumed
  - Removed CH3SHIFT
  - \ref CS2BACKBONE is now native in PLUMED removing the need to link ALMOST, small syntax differences
  - \ref CS2BACKBONE, \ref NOE, \ref RDC, removed the keyword ENSEMBLE: now ensemble averages can only be calculated using \ref ENSEMBLE
  - \ref RDC, syntax changes
  - It is not possible anymore to select modules using `modulename.on` and `modulename.off` files. Use `./configure --enable-modules` instead.
  - Removed IMD modules. In case someone is interested in restoring it, please contact the PLUMED developers.
- New actions:
  - \ref FIXEDATOM
  - \ref HBOND_MATRIX
  - \ref CLUSTER_PROPERTIES
  - \ref CLUSTER_DIAMETER
  - \ref CLUSTER_NATOMS
  - \ref OUTPUT_CLUSTER
  - \ref CLUSTER_DISTRIBUTION
  - \ref ROWSUMS
  - \ref COLUMNSUMS
  - \ref UPDATE_IF
  - \ref DUMPGRID
  - \ref DUMPCUBE
  - \ref CONVERT_TO_FES
  - \ref INTERPOLATE_GRID
  - \ref FIND_CONTOUR
  - \ref FIND_SPHERICAL_CONTOUR
  - \ref FIND_CONTOUR_SURFACE
  - \ref AVERAGE
  - \ref REWEIGHT_BIAS
  - \ref REWEIGHT_TEMP
  - \ref REWEIGHT_METAD
  - \ref PCA
  - \ref PRE
  - \ref STATS
  - \ref METAINFERENCE
  - \ref LOCALENSEMBLE
  - \ref FRET
  - \ref RESET_CELL
  - \ref JCOUPLING
  - \ref ERMSD
- New features in MD patches (require re-patch):
  - Patch for amber 14 now passes charges with appropriate units (fixes \issue{165}). Notice that
    the patch is still backward compatible with older PLUMED version, but the charges will only be passed
    when using PLUMED 2.3 or later.
  - Patch for GROMACS 5.1 incorporates Hamiltonian replica exchange, see \ref hrex
  - Gromacs 2016, 5.1.x, 5.0.x, flush the plumed output files upon checkpointing
  - Added patch for Gromacs 2016.1
  - gromacs 5.1.x patch updated to 5.1.4
  - Removed the patch for Gromacs 4.6.x 
  - LAMMPS patch updated to support multiple walkers and report plumed bias to LAMMPS (thanks to Pablo Piaggi).
- New features for existing actions:
  - The SPECIES and SPECIESA keyword in MultiColvars can now take a multicolvar as input.  This allows one
    to calculate quantities such as the Q4 parameters for those atoms that have a coordination number greater
    than x.
  - Added MATHEVAL type in \ref switchingfunction
  - Added Q type native contacts in \ref switchingfunction (thanks to Jan Domanski).
  - \ref COMMITTOR can now be used to define multiple basins
  - The number of atoms admitted in \ref BRIDGE has been significantly increased, see \issue{185}.
  - \ref driver now allows --trajectory-stride to be set to zero when reading with --ixtc/--itrr. In this case, step number is read from the trajectory file.
  - \ref METAD and \ref PBMETAD can now be restarted from a GRID 
  - Added keywords TARGET and DAMPFACTOR in \ref METAD
  - When using \ref METAD with file-based multiple walkers and parallel jobs (i.e. mpirun) extra suffix is not added (thanks to Marco De La Pierre).
  - \ref ENSEMBLE added keywords for weighted averages, and calculation of higher momenta
  - \ref MOLINFO now allows single atoms to be picked by name.
  - \ref FIT_TO_TEMPLATE now supports optimal alignment.
  - \ref CONSTANT added the possibility of storing more values as components with or without derivatives
  - \ref PUCKERING now supports 6 membered rings.
  - Extended checkpoint infrastructure, now \ref METAD and \ref PBMETAD will write GRIDS also on checkpoint step (only the GROMACS patch
    is currently using the checkpointing interface)
- Other features:
  - Added a plumed-config command line tool. Can be used to inspect configuration also when cross compiling.
  - Added a `--mpi` option to `plumed`, symmetric to `--no-mpi`. Currently, it has no effect (MPI is initialized by default when available).
  - PLUMED now generate a VIM syntax file, see \ref VimSyntax
  - The backward cycle is now parallelized in MPI/OpenMP in case many collective variables are used.
  - GSL library is now searched by default during `./configure`.
  - Tutorials have been (partially) updated to reflect some of the changes in the syntax
  - Parser now reports errors when passing numbers that cannot be parsed instead of silently replacing their default value. See \issue{104}.
  - More and more documentation
- Bug fixes:
- Fixed a bug in \ref PBMETAD that was preventing the writing of GRIDS if a hill was not added in that same step 

For developers:
- IMPORTANT: BIAS can now be BIASED as well, this changes can lead to some incompatibility: now the "bias" component is always defined automatically
  by the constructor of Bias as a componentWithDerivatives, derivatives are automatically obtained by forces. The main change is that you don't have to define
  the bias component anymore in your constructor and that you can use setBias(value) to set the value of the bias component in calculate. 
- Added new strings for plumed cmd: setMDMassUnits, setMDChargeUnits, readInputLine, performCalcNoUpdate, update and doCheckPoint.
- Easier to add actions with multiple arguments
- New functions to access local quantities in domain decomposition
- Active modules to enable regtests are chosen using `plumed config`.
- A script is available to check if source code complies plumed standard. Notice that this script is run together with cppcheck on travis-ci.
- Cppcheck on travis-ci has been updated to 1.75. Several small issues triggering errors on 1.75 were fixed (e.g. structures passed by value
   are now passed by const ref) and false positives marked as such.
- Added coverage scan.

## Version 2.3.1 (Mar 31, 2017)

- Fix to FIT_TO_TEMPLATE as in 2.2.5. Notice that in 2.3.0 also the case with TYPE=OPTIMAL was affected. This is fixed now.
- small change in \ref CS2BACKBONE to symmetrize the ring current contribution with respect to ring rotations (also faster)
- fixed `plumed-config` that was not working.
- log file points to the `config.txt` files to allow users to check which features were available in that compiled version.
- `make clean` in root dir now also cleans `vim` sub-directory.
- Updated gromacs patch to version 2016.3 

For developers:
- Cppcheck on travis-ci has been updated to 1.77.
- Doxygen on travis-ci has been updated to 1.8.13

## Version 2.3.2 (Jun 12, 2017)

See branch \branch{v2.3} on git repository.

- Resolved problem with nan in \ref SMAC with SPECIESA and SPECIESB involving molecules that are the same
- PDB reader is now able to read files with dos newlines (see \issue{223}).
- Fixed bug in \ref CS2BACKBONE (v2.3.1) related to ring currents of HIS and TRP
- Fixed bug in if condition in \ref PCAVARS so that you can run with only one eigenvector defined in input 
- Fixed bug with timers in \ref sum_hills \issue{194}.
- Fixed bug when using \ref MOVINGRESTRAINT with periodic variables such as \ref TORSION \issue{225}.
- Fixed bug in \ref HBOND_MATRIX that used to appear when you used DONORS and ACCEPTORS with same numbers of atoms 
- Fixed bug in \ref DISTANCES that appears when using BETWEEN and link cells.
- Prevented users from causing segfaults by storing derivatives without LOWMEM flag.  In these cases PLUMED crashes with meaningful errors.
- Fixed bug in \ref HISTOGRAM that causes NaNs when using KERNEL=DISCRETE option
- Fixed a bug in the parser related to braces, see \issue{229}
- Fixed a bug that appeared when using \ref Q3, \ref Q4 and \ref Q6 with LOWEST or HIGHEST flag
- Fixed a bug that appears when you use \ref MFILTER_LESS as input to \ref COORDINATIONNUMBER with SPECIESA and SPECIESB flags
- Fixed a bug that was making flushing when gromacs checkpoints not functional (thanks to Summer Snow).
- Fixed a bug affecting \ref EXTENDED_LAGRANGIAN and \ref METAD with ADAPT=DIFF when using an argument
  with periodicity (min,max) such that min is different from -max.
  This does not affect normal \ref TORSION, but would affect \ref PUCKERING component phi
  with 6-membered rings. In addition, it would affect any variable that is created by the user with a periodicity
  domain not symmetric around zero. See \issue{235} (thanks to Summer Snow for reporting this bug).
- Fixed numerical issue leading to simulations stuck (LatticeReduction problem) with intel compiler and
  large simulation cells.
- Fixed a bug affecting \ref LOCAL_AVERAGE and outputting all multicolvars calculated by \ref Q6 with \ref DUMPMULTICOLVAR
- `plumed info --user-doc` and `plumed info --developer-doc` now fall back to online manual when local doc is not installed,
  see \issue{240}.

For developers:
- IMPORTANT: we started to enforce code formatting using astyle. Check the developer documentation to learn how to
  take care of not-yet-formatted branches.
- plumedcheck validation has been made stricter. All the checks are now described in the developer manual.
- New flag `--disable-libsearch` for `configure`, allowing an easier control of linked libraries when installing PLUMED
  with a package manager such as MacPorts.
- Added `--disable-static-patch` to `./configure` to disable tests related to static patching. It can be used
  when static patching is not needed to make sure a wrong c++ library is not linked by mistake.
- Using `install_name_tool` to fix the name of the installed library on OSX. Allows linking the PLUMED
  shared library without explicitly setting `DYLD_LIBRARY_PATH`.
- Added environment variable `PLUMED_ASYNC_SHARE` to enforce synchronous/asynchronous atom sharing (mostly for debug purpose).
- On travis-ci, using ccache to speedup builds.
- On travis-ci, added a regtest using Docker with gcc6 and MPI.
- On travis-ci, docs for unofficial or unsupported branches are set not to be indexed by search engines (see \issue{239})
- Cppcheck on travis-ci has been updated to 1.79.

## Version 2.3.3 (Oct 3, 2017)

For users:
- Fixed a bug in \ref switchingfunction MATHEVAL, leading to inconsistent results when using OpenMP with multiple threads (see \issue{249}).
- \ref FIT_TO_TEMPLATE now reports when it is used with a reference file with zero weights.
- Fixed logging of \ref UNITS (thanks to Omar Valsson).
- Fixed a possible bug with \ref EFFECTIVE_ENERGY_DRIFT and domain decomposition with a domain containing zero atoms.


For developers:
- Fixed a bug in `./configure --disable-libsearch` when searching for molfile plugins.
- Cppcheck on travis-ci has been updated to 1.80.
- Configure script now has a list of better alternatives to find a working `ld -r -o` tool to merge object files.
  This solves linking issues on some peculiar systems (see \issue{291}, thanks to Massimiliano Culpo). 
- Using `install_name_tool` also on non-installed libraries. This makes it possible to link them and later
  find them without explicitly setting `DYLD_LIBRARY_PATH`. This should also make the `DYLD_LIBRARY_PATH` irrelevant.
  Notice that `DYLD_LIBRARY_PATH` is not well behaved in OSX El Capitan.

## Version 2.3.4 (Dec 15, 2017)

For users:
- GROMACS patch updated to gromacs-2016.4. This patch was also fixed in order to properly work with \ref ENERGY (see \issue{316})
  and to implement `-hrex` option (see \issue{197}).
- Patch for GROMACS 5.1.4 updated to fix an error with \ref ENERGY (see \issue{316}).
- Solved a bug in \ref ERMSD leading to incorrect results when using non-default length units (e.g. with `UNITS LENGTH=A`).

For developers:
- Regtest script also reports when exitcode different from zero is returned.
- Patch script reports errors returning a nonzero exit code.
- cppcheck update to 1.81
- Solved small bug in stored PLUMED_ROOT directory as obtained from statically patched MD codes.
  Namely, the compilation directory was stored rather than the installation one.

## Version 2.3.5 (Mar 2, 2018)

For users:
- Fixed `plumed partial_tempering` to agree with GROMACS conventions for the choice of dihedral angles (see \issue{337}).
  Should be irrelevant for the vast majority of cases.
- Fixed small bug in regexp parser - the part outside the parentheses was just ignored.

For developers:
- Doxygen on travis-ci has been updated to 1.8.14.
- Embedded astyle updated to 3.1.
- `make clean` now correctly removes the `src/lib/plumed` executable.

## Version 2.3.6 (Jul 2, 2018)

For users:
- Fixed a problem leading to NaN derivatives of \ref switchingfunction `Q` when distance between two atoms is large.
- GROMACS patch updated to gromacs-2016.5.
- `./configure` crashes if prefix is set to present working directory (notice that this choice was already leading to issues).
- \ref DUMPATOMS reports an error when trying to write xtc/xdr files without the xdrfile library installed.
- Fixed a bug appearing when using \ref PATH or \ref GPROPERTYMAP with virtual atoms without simultaneously using the same
  atoms in a different action.
- Fixed incorrect format of the pdb file written by \ref PCA (see \issue{363}).
- Fixed behavior of natural units. When an MD code asks for natural units, it is not necessary to also set units within PLUMED using \ref UNITS (see \issue{364}).

For developers:
- Fixed small issue in debug options of \ref driver (see \issue{245}).
- `plumed patch -e` now accepts a name closely matching the patch name (e.g. `plumed patch -e gromacs2016.5` will try to patch
  even if the stored patch is for `gromacs-2016.4`). This simplifies managing Portfiles. Nothing changes when picking the patch
  from the interactive menu.
- Install newer ccache on travis-ci, build faster.
- Small fix in provided env modules (`PLUMED_VIMPATH` is set also when shared libraries are disabled).

## Version 2.3.7 (Oct 5, 2018)

For users:
- Fixed flag DETAILED_TIMERS in \ref DEBUG (flag was ignored and detailed timers always written).
- Small fix in \ref DUMPMASSCHARGE (atoms are now correctly requested only at first step).

## Version 2.3.8 (Dec 19, 2018)

\plumednotmaintained

For users:
- Fixed some openMP regression (some related to the whole codes and some specifics for Coordination and Multicolvar), this were compiler dependent so not all users may have experienced them
- Fixed an issue with \ref CS2BACKBONE when more than 2 chains were used
- Fixed memory leak in \ref RDC.
- Fixed segmentation fault with more than two CVs in reweighting \ref METAD (see \issue{399}, thanks to Fiskissimo).

For developers:
- Small fix in LDFLAGS when enabling coverage.
- Fixed order of flags in tests for static linking done by configure (see \issue{407}).
- Fixed the way paths are hard-coded so as to facilitate conda packaging (see \issue{416}).


*/
@page CHANGES-UNRELEASED Unreleased changes 

This page contains changes that will end up in 2.7

Changes from version 2.6 which are relevant for users:

- New contributed modules:
  - A new Funnel module by Stefano Raniolo and Vittorio Limongelli 
     - \ref FUNNEL_PS 
     - \ref FUNNEL 


For developers:
- small fix in `Plumed.h` too avoid unique global symbols (see \issue{549})


@page CHANGES-2-0 Version 2.0

Version 2.0.0 (September 27, 2013)
----------------------------

Version 2.0 is a complete rewrite, so there is no way to write a complete set of difference
with respect to plumed 1.3. Here is a possibly incomplete summary of the difference:
- The input is simpler, more flexible, and more error proof.
  Many checks are now performed and in this way common errors are avoided. 
- The units are now the same for all MD codes.
  If you want to use a different unit than the default you set it in the input file. 
- The analysis tools are now much more flexible.
  As an example of this it is now possible to write different collective variables with different frequencies.
- Many complex collective variables are considerably faster than they were in plumed1.
  In particular, all variables based on RMSD distances. 
- Centers of mass can be used as if they were atoms.
  Hence, unlike plumed 1.3, you can use center of mass positions in ALL collective variables.
- The virial contribution is now computed and passed to the MD code.
  Plumed can thus now be used to perform biased NPT simulations.
- Variables can be dumped on different files, and are
  computed only when this is necessary.
- PLUMED is now compiled as a separate library. This simplifies the patching
  procedure, but might require some extra work to configure PLUMED properly.
  Since PLUMED can be loaded as a shared library, it is possible to setup
  everything such that PLUMED and MD codes can be updated independently from each
  other.

In addition, it is now much easier to contribute new functionality to the code because: 
- There is a much simpler interface between plumed and the base MD codes.
  This makes it much easier to add plumed to a new MD code. Hopefully, in the future,
  interfaces with MD codes will be maintained by the developers of the MD codes
  independently from PLUMED developers. This will allow more MD codes
  to be compatible with PLUMED.
- There is C++ object oriented programming and full compatibility with the C++ standard library 
- A modular structure.
- New collective variables and methods can be released independently.
- There is an extensive developer documentation.
- User documentation is provided together inside the implementation files.

Caveats:
- PLUMED 2 input file (plumed.dat) has a syntax which is not
  compatible with PLUMED 1.
  Transition should be easy, but cannot
  be done just using the new version with the old input file.
- PLUMED 2 is written in C++, thus requires a C++ compiler
- PLUMED 2 may not include all the features that were available
  in PLUMED 1.

A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Version 2.0.1 (Nov 14, 2013)
----------------------------

For users:
- Fixed a bug in \ref HISTOGRAM with REWEIGHT_BIAS. Reweighting was only done when also temperature-reweighting was enabled.
- Fixed a bug that was sometime crashing code with domain decomposition and
  non-dense simulation boxes (e.g. implicit solvent).
- Performance improvements for \ref GYRATION.
- Flush all files every 10000 steps by default, without need to use \ref FLUSH
- Errors when writing input for \ref switchingfunction are now properly
  recognized.
- Added message when \ref simplemd is used on a non-existing file.
- Fixed `plumed mklib` such that it deletes the target shared library in case
  of compilation error.
- Several small fixes in documentation and log file.

For developers:
- Added possibility to setup replica exchange from MD codes in Fortran (commands "GREX setMPIFIntercomm" and "GREX setMPIFIntracomm").
- cmd("setStopFlag") should now be called after PLUMED initialization.
- Several small fixes in documentation.

Version 2.0.2 (Feb 11, 2014)
----------------------------

For users:
- Fixed bug with \ref METAD with INTERVAL and replica exchange, including bias exchange.
  Now the bias is correctly computed outside the boundaries. Notice that this is different
  from what was done in PLUMED 1.3. Also notice that INTERVAL now works
  correctly with grids and splines.
- Fixed bug with \ref READ and periodic variables.
- Fixed bug with \ref HISTOGRAM (option USE_ALL_DATA was not working properly).
- Gromacs patch updated to 4.6.5.
- Gromacs patch for 4.6 has been modified to allow for better load balancing when
  using GPUs.
- Added option 'plumed info --long-version' and 'plumed info --git-version'.
- Added full reference (page/number) to published paper in doc and log.
- Fixed a bug in file backups (only affecting Windows version - thanks to T. Giorgino).
- Added possibility to search in the documentation.
- Several small fixes in documentation and log file.

For developers:
- Fixed Makefile dependencies in some auxiliary files in src/lib (*cmake and *inc).
- Changed way modules are linked in src/.
  E.g. src/colvar/tools/ is not anymore a symlink to src/colvar but a real directory.
  (Notice that this introduces a regression: when using plumed as an external library
  some include files could not work - this only applies when plumed is installed;
  also notice that this is fixed in 2.0.3)
- Patch for gromacs 4.6 now also include original code so as to simplify its modification.
- Added option 'plumed patch --save-originals'.
- Fixed regtest regtest/secondarystructure/rt32 to avoid problems with NUMERICAL_DERIVATIVES.
- Removed include graphs in the documentation (too large).
- Several small fixes in documentation.

Version 2.0.3 (June 30, 2014)
----------------------------

For users:
- Now compiles on Blue Gene Q with IBM compilers.
- Fixed bug in \ref CENTER where default WEIGHTS were missing. 
- Fixed broken \ref CONTACTMAP with SUM
- Fixed \ref DUMPATOMS with gro file and more than 100k atoms.
- Added CMDIST in \ref CONTACTMAP to emulate plumed1 CMAP.
- Several small fixes in documentation and log file.

For developers:
- Fixed cmd("getBias") to retrieve bias. It was not working with
  single precision codes and it was not converting units properly.
- Fixed a regression in 2.0.2 concerning include files from installed plumed
  (see commit 562d5ea9dfc3).
- Small fix in tools/Random.cpp that allows Random objects to be
  declared as static.
- Small fix in user-doc compilation, so that if plumed is not found
  the sourceme.sh file is sourced
- Fixed non-ANSI syntax in a few points and a non-important memory leakage.
- Split cltools/Driver.cpp to make parallel compilation faster.

Version 2.0.4 (September 15, 2014)
----------------------------------------------

For users:
- Fixed a bug in \ref BIASVALUE that could produce wrong acceptance with replica exchange simulations.
- Fixed a few innocuous memory leaks.
- Fixed reader for xyz files, that now correctly detects missing columns. Also a related regtest has
  been changed.
- Several small fixes in documentation and log file.

For developers:
- Renamed Value.cpp to BiasValue.cpp

Version 2.0.5 (December 15, 2014)
----------------------------------------------

\plumednotmaintained

For users:
- Fixed a bug in replica exchange with different Hamiltonians (either lambda-dynamics
  or plumed XX-hrex branch) possibly occurring when using charge or mass dependent
  variables.
- Fixed a bug in analysis (e.g. \ref HISTOGRAM) leading to wrong accumulation
  of statistics when running a replica exchange simulation.
- Fixed a bug in the calculation of derivatives in histograms. This should
  be harmless since people usually only consider the value in histograms
  and not the derivatives.
- Fixed an issue in Makefile that could results in problems when
  patching an MD code with --shared option (pointed out by Abhi Acharya).
  This fixes a regression introduced in 2.0.2.
- Small fixes in documentation.

For developers:
- Added warning when performing regtests using an instance of plumed from
  a different directory

@page CHANGES-2-4 Version 2.4

## Version 2.4 (Dec 15, 2017)

Version 2.4 contains several improvements with respect to 2.3. Users currently working with 2.3
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files.
Notice that version 2.4 includes already all the fixes in branch 2.3 up to 2.3.3 indicated in \ref CHANGES-2-3 .

Changes from version 2.3 which are relevant for users:
- Changes leading to incompatible behavior:
  - A c++11 compliant compiler is required (see \issue{212}). This should mean:
    - gcc 4.8
    - clang 3.3
    - intel 15
    Since the number of c++11 features that we use is limited, older compilers might work as well.
  - The meaning of `BIASFACTOR=1` in \ref METAD has been modified and can now be used to indicate unbiased
    simulations. Non-well-tempered metadynamics is BIASFACTOR=-1, which is the new default value.
    Notice that this has an implication on the bias factor written in the HILLS file when doing
    non-well-tempered metadynamics.
  - Due to a change in \ref COMMITTOR, the format of its output file has been slightly changed.
  - \ref HISTOGRAM : When using weights default is now to output histogram divided by number of frames from which data was taken.  In addition the 
    UNORMALIZED flag has been replaced with the keyword `NORMALIZATION`, which can be set equal to true, false or ndata.
  - All switching functions are now stretched by default, also when using the "simple syntax" (e.g. `COORDINATION NN=6`).
    Switching functions were already stretched by default when using the advanced syntax (e.g. `COORDINATION SWITCH={}`)
    since version 2.2.  Notice that this will introduce small numerical differences in the computed switching functions.
- New modules:
  - A new PLUMED-ISDB module have been included, this module includes a number of CVs to calculate experimental data with the internal ability
    to also calculate a \ref METAINFERENCE score.
    - New actions include:
      - \ref EMMI
      - \ref SAXS
      - \ref RESCALE, \ref SELECT, \ref SELECTOR
    - Updated actions include:
      - \ref CS2BACKBONE
      - \ref FRET
      - \ref JCOUPLING
      - \ref METAINFERENCE
      - \ref NOE
      - \ref PRE
      - \ref RDC, \ref PCS
      - \ref PBMETAD
  - A new EDS module have been included, contributed by Glen Hocky and Andrew White.
    This module implements the following methods:
    - \ref EDS
  - A new DRR module have been included, contributed by Haochuan Chen and Haohao Fu.
    This module implements the following methods:
    - \ref DRR
    - \ref drr_tool
  - A new VES module have been included, contributed by Omar Valsson.
    This module implements the following methods:
    - \ref BF_CHEBYSHEV
    - \ref BF_COMBINED
    - \ref BF_COSINE
    - \ref BF_CUSTOM
    - \ref BF_FOURIER
    - \ref BF_LEGENDRE
    - \ref BF_POWERS
    - \ref BF_SINE
    - \ref OPT_AVERAGED_SGD
    - \ref OPT_DUMMY
    - \ref TD_CHI
    - \ref TD_CHISQUARED
    - \ref TD_CUSTOM
    - \ref TD_EXPONENTIAL
    - \ref TD_EXPONENTIALLY_MODIFIED_GAUSSIAN
    - \ref TD_GAUSSIAN
    - \ref TD_GENERALIZED_EXTREME_VALUE
    - \ref TD_GENERALIZED_NORMAL
    - \ref TD_GRID
    - \ref TD_LINEAR_COMBINATION
    - \ref TD_PRODUCT_COMBINATION
    - \ref TD_PRODUCT_DISTRIBUTION
    - \ref TD_UNIFORM
    - \ref TD_VONMISES
    - \ref TD_WELLTEMPERED
    - \ref VES_LINEAR_EXPANSION
    - \ref VES_OUTPUT_BASISFUNCTIONS
    - \ref VES_OUTPUT_FES
    - \ref VES_OUTPUT_TARGET_DISTRIBUTION
    - \ref ves_md_linearexpansion
- New collective variables:
  - \ref DIMER (thanks to Marco Nava).
  - \ref EEFSOLV : EEF1 implicit solvent solvation energy
  - \ref ADAPTIVE_PATH : Adaptive path variables using the method from \cite BerndAdaptivePath
- New actions:
  - \ref INENVELOPE
  - \ref TOPOLOGY_MATRIX
  - \ref BOND_DIRECTIONS
  - \ref DUMPGRAPH
  - \ref GRID_TO_XYZ
  - \ref INTEGRATE_GRID
  - \ref LWALLS
  - \ref MAXENT
  - \ref MCOLV_COMBINE
  - \ref MCOLV_PRODUCT
  - \ref POLYMER_ANGLES
  - \ref XANGLES , \ref YANGLES , \ref ZANGLES
  - \ref XYTORSIONS , \ref XZTORSIONS , \ref YXTORSIONS , \ref YZTORSIONS , \ref ZXTORSIONS , and \ref ZYTORSIONS
- New command line tools:
  - \ref pesmd : Tool for performing Langevin dynamics on an energy landscape that is specified using a PLUMED input file
  - \ref pathtools 
- Other changes:
  - Sharing coordinates and applying force is now faster (in some cases these can result in much better scaling of the performances in parallel).
  - \ref COMMITTOR : new flag to use committor to keep track of the visited basins without stopping the simulation
  - \ref PBMETAD : multiple walkers using files (thanks to Marco De La Pierre).
  - \ref PBMETAD : adaptive Gaussian kernels
  - \ref PBMETAD : default names for `GRID` and `FILE` (useful with many collective variables) 
  - \ref METAD : BIASFACTOR=1 is allowed and performs unbiased sampling. HILLS file can be used
    to recover free energy also in this case.
  - \ref METAD : a RECT option is available that allows setting an array of bias factors, one for each replica.
  - \ref METAD : added options to perform Transition Tempered Metadynamics (thanks to James Dama)
  - \ref PATHMSD and \ref PROPERTYMAP now support alignment to a close structure (thanks to Jana Pazurikova)
  - PDB files with more than 100k atoms can now be read using [hybrid 36](http://cci.lbl.gov/hybrid_36/) format,
    see \issue{226}.
  - Added lepton support. Set env var `export PLUMED_USE_LEPTON=yes` to activate lepton as a matheval replacement
    in \ref MATHEVAL, \ref CUSTOM, and \ref switchingfunction "MATHEVAL switching function".
    Notice that in v2.5 matheval support will be dropped and all these keywords will use lepton.
    See \issue{244}.
  - When parsing constants, PLUMED uses lepton library. This allows to pass
    arguments such as `HEIGHT=exp(0.5)` (see \ref parsing-constants).
  - \ref CUSTOM function has been added as an alias to \ref MATHEVAL .
  - Trajectories read in \ref driver also support the usual replica convention, that is if
    trajectory with replica suffix is not found the driver will look for a trajectory without the replica suffix.
  - A new syntax (`@replicas:`) can be used to specify different arguments for different replicas (see \ref special-replica-syntax).
  - Internal molfile implementation has been updated to VMD 1.9.3.
  - Examples in the documentation now have syntax highlighting and links to the documentation of used actions.
  - \ref COORDINATIONNUMBER : Added option to have pairwise distance moments of coordination number in the multicolvar module
  - GROMACS patch updated to gromacs-2016.4
  - Implemented HREX for gromacs-2016.4.
  - Added patch for Quantum ESPRESSO 6.2 (thanks to Ralf Meyer).
  - Fixed a bug in \ref LOCAL_AVERAGE which appears when you use `SPECIESA` and `SPECIESB` keywords instead of just `SPECIES`
  - Added possibility to pass `--kt` from \ref driver.

Changes from version 2.3 which are relevant for developers:
  - A few fixes has been made to improve exception safety. Although we still cannot declare
    PLUMED totally exception safe (there are still many non-safe pointers around),
    this made it possible to add a regtest that actually tests erroneous cmd strings
    and erroneous inputs.
  - Due to the required c++11 support, travis-ci test on Ubuntu Precise has been removed.
  - `gettimeofdate` and `gettime` have been replaced with portable `chrono` classes introduced in c++11.
  - C++ exceptions are enabled by default.
  - A large number of loops have been changed to use the `auto` keyword in order to improve code readability.
  - Stack trace is not written upon error anymore, unless environment variable `PLUMED_STACK_TRACE` is set at runtime.
  - Fixed a potential bug using single precision system BLAS on a mac (notice that currently plumed only uses
    double precision, so it is harmless).
  - Added `--enable-rpath` option for autoconf (off by default).
  - Files related to changelog are now stored as `.md` files. This makes
    it possible to navigate them from github.
  - `configure.ac` has been simplified and improved in order to more easily probe C++ libraries.
  - added `plumed_custom_skip` function to regtests in order to skip specific tests based on specific conditions (e.g. OS).
  - environment variable `LDSO` has been renamed to `LDSHARED`, which is standard in the python community.
  - a `libplumedWrapper.a` library is installed as well, that is used in `--runtime` patching.
  - pkgconfig files are installed.
  - `plumed config makefile_conf` can be used to retrieve `Makefile.conf` file a posteriori.
  - Store `MPIEXEC` variable at configure time and use it later for running regtests. Notice that in case
    `MPIEXEC` is not specified regtests will be run using the command stored in env var `PLUMED_MPIRUN` or, if this is
    also not defined, using `mpirun`.
  - Added canonical Makefile targets `check` and `installcheck`. Notice that `check` runs checks with
    non-installed plumed whereas `installcheck` uses the installed one, including its correct program name if it
    was personalized (e.g. with suffixes). Notice that this modifies the previously available `check` target.
    
    
## Version 2.4.1 (Mar 2, 2018)

For users:
  - Fixed an important bug affecting RMSD calculations with compilers supporting OpenMP 4 (e.g.: intel compiler). Notice that this bug might potentially affect not only
    \ref RMSD variable, but also \ref PATHMSD variables using RMSD, \ref FIT_TO_TEMPLATE, \ref PCAVARS, and possibly other variables based on RMSD calculations and optimal alignments
    (see \issue{343}). Results might depend on the exact architecture and on how aggressive is the compiler. The bug is a consequence of some erroneous SIMD directives introduced in 2.4.0, so it does not affect PLUMED 2.3.x. 
  - Resolved a problem with \ref CS2BACKBONE and glycine atom names.
  - Module VES: Fixed a bug with basis functions that have a constant function different from 1 (e.g. scaled version of the Legendre basis functions, \ref BF_LEGENDRE) that was causing a time-dependent shift in the bias potential.
  - Module VES: In optimizers (\ref OPT_AVERAGED_SGD and \ref OPT_DUMMY) the output of quantities related to the instantaneous gradients are now off by default as these quantities are generally not useful for normal users, their output can instead by re-enabled by using the `MONITOR_INSTANTANEOUS_GRADIENT` keyword. Also added an keyword `MONITOR_AVERAGE_GRADIENT` that allows to monitor the averaged gradient and output quantities related to it. 
  - \ref RMSD variable and other collective variables using reference PDB files now crash when zero weights are passed (see \issue{247}).
  - Using \ref COM with \ref driver without passing masses now triggers an error instead of reporting NaNs (see \issue{251}).

For developers:
  - `plumed patch -p` command can be used twice without triggering an error. This will allow e.g. building again
    on MacPorts in cases where the build was interrupted. Notice that this only works for patches without special
    after/before patch/revert functions.

## Version 2.4.2 (Jul 2, 2018)

For users:
  - All fixes done in version 2.3.6. Notice that \issue{363} in version 2.4 also applies to \ref pathtools.
  - Additional residue names (without the prefix `D`) are now supported by \ref MOLINFO for DNA. See \issue{367}.
  - Solved an important bug appearing in NAMD interface. Notice that the bug was a regression introduced in 2.4.0. As consequence, versions <= 2.3 and versions >=2.4.2
    are expected to work correctly. See \issue{254}.
  - GROMACS patch for gromacs-2018.1.
  - \ref VimSyntax now highlights `__FILL__` strings.
  - \ref METAD and \ref PBMETAD give a warning when one restarts a simulation and the old hills file is not found. See \issue{366}.

For developers:
  - `LDSHARED` is now correctly taken into account when launching `./configure`.
  - Fixed installation with `--disable-shared`.
  - Cppcheck upgraded to 1.84.

## Version 2.4.3 (Oct 5, 2018)

For users:
  - All fixes done in version 2.3.7.
  - Module VES: Fixed a bug in `TD_GRID` for 2D grids where the grid spacing is not the same for both dimensions.
  - GROMACS patch for gromacs-2018.3.

## Version 2.4.4 (Dec 19, 2018)

For users:
  - Fixed some performances regression issue with OpenMP
  - Updated NAMD patches to version 2.12 and 2.13. Old patches have been removed.
  - GROMACS patch for gromacs-2018.4.
  - Fixed a thread safety issue using forces on \ref HISTOGRAM 
  - Fixed error message suggesting wrong actions (see \issue{421}).

For developers:
  - All fixed done in version 2.3.8
  - Cppcheck updated to 1.85

## Version 2.4.5 (Apr 1, 2019)

For users:
  - Fixed an inconsistency in parsing of braces.
    It is now possible to pass individual options
    including spaces (e.g. with `FILE={/path with space/file}`). Notice 
    that this invalidates syntax such as `ATOMS={1}{2}{3}{4}`. See more
    at \issue{434}.
  - Fixed \ref simplemd so as to call "runFinalJobs" at the end of the simulation.
  - GROMACS patch for gromacs-2016.6.
  - GROMACS patch for gromacs-2018.6.
  - Added aliases for some actions/options containing dashes (`-`) in their name. This will improve
    backward compatibility when these actions/options will be removed (see \issue{449}).

## Version 2.4.6 (Jul 19, 2019)

For users:
  - Fixed a bug in \ref COORDINATIONNUMBER where derivatives were wrong when using `R_POWER` > 2, thanks to `@MoleOrbitalHybridAnalyst` for spotting and fixing
  - Fixed a bug in library search, possibly affecting linked blas/lapack on OSX (see \issue{476}).
  - Fixed a bug in \ref METAD with `TARGET` and `GRID_SPARSE` (see \issue{467}).

## Version 2.4.7 (Jan 27, 2020)

\plumednotmaintained

For users:
  - Fixed a bug with \ref CONVERT_TO_FES and periodic variables, see \issue{441} (backported from v2.5.3).
  - More robust backup for output files when running over multiple processes
  - Fixed a regression in the performances of `GEOMETRY` based flexible hills in \ref METAD and \ref PBMETAD
  - Fixed \issue{538}.
  - Fixed potential issue with VMD plugins from 1.9.4 (\issue{545}, thanks to Lixin Sun).
  - Module VES: Fixed an off-by-one bug in the output of target distribution averages. The bug only affects the output and does not affect results. The bug also affected the output of coefficients when using a bias cutoff. 
  - Module VES: Made sure that all relevant output files are written out at the final step when shutting down the simulation. This solves issues reported by @PabloPiaggi with restarting when there is a mismatch been the output of files and the number of MD steps. 


Version 1.0 of Tinker-HP
Version 1.1 of Tinker-HP
[![Homepage](https://img.shields.io/badge/Home-plumed.org-green.svg)](http://www.plumed.org)
[![Homepage](https://img.shields.io/badge/Google_group-plumed--users-green.svg)](http://groups.google.com/forum/#!forum/plumed-users)
[![codecov](https://codecov.io/gh/plumed/plumed2/branch/master/graph/badge.svg)](https://codecov.io/gh/plumed/plumed2)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/plumed/plumed2.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/plumed/plumed2/context:python)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/plumed/plumed2.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/plumed/plumed2/context:cpp)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](http://www.gnu.org/licenses/lgpl-3.0)
[![Github Releases](https://img.shields.io/github/release/plumed/plumed2.svg)](https://github.com/plumed/plumed2/releases)
[![MacPorts package](https://repology.org/badge/version-for-repo/macports/plumed.svg)](https://repology.org/project/plumed/versions)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/plumed/badges/version.svg)](https://anaconda.org/conda-forge/plumed)
[![AUR package](https://repology.org/badge/version-for-repo/aur/plumed.svg)](https://repology.org/project/plumed/versions)
[![Twitter Follow](https://img.shields.io/twitter/follow/plumed_org.svg?style=social&label=Follow)](https://twitter.com/plumed_org)

Branches and releases
---------------------

Several branches and tags are stored on the git repository.

Branches named `v2.X` correspond to release branches.

Master branch may contain non tested features and is not expected to be used by non-developers.
It typically contains features that will be available on the next release.

Tags named `v2.XbY` correspond to beta releases, use it with care.
Tags named `v2.X.Y` correspond to official releases, use the latest available.

In addition, the repository contains a number of other branches related to specific features.
Please contact the developers that are committing on those branches before basing your work
there, since they might contain temporary work and might be rebased later.
For instance, branch `testdoc` is setup so as to push a test copy of the manual
and is often force pushed.

To report problems found on beta or official releases, use the normal
[plumed-users@googlegroups.com](mailto:plumed-users@googlegroups.com)
mailing list. Please state exactly which version you are using.
To report problems found on `master` branch, use the
[plumed2-git@googlegroups.com](plumed2-git@googlegroups.com) mailing list.
This is also the correct place for discussions about new features etc.
When reporting please provide the git hash (you can obtain it with `git rev-parse HEAD`).

Status
------

Below you find the status on [Travis-CI](http://travis-ci.org/plumed/plumed2) for the release branches.

| Branch   |      Status   | First stable release (year) | Still supported |
|:--------:|:-------------:|:--------:|:------:|
| master   | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=master)](https://travis-ci.org/plumed/plumed2) | 2020 (expected) | / |
| v2.6     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.6)](https://travis-ci.org/plumed/plumed2)   | 2019 | yes |
| v2.5     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.5)](https://travis-ci.org/plumed/plumed2)   | 2018 | yes |
| v2.4     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.4)](https://travis-ci.org/plumed/plumed2)   | 2017 | no |
| v2.3     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.3)](https://travis-ci.org/plumed/plumed2)   | 2016 | no |
| v2.2     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.2)](https://travis-ci.org/plumed/plumed2)   | 2015 | no |
| v2.1     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.1)](https://travis-ci.org/plumed/plumed2)   | 2014 | no |
| v2.0     | Not available | 2013 | no |

Content
-------

Here's a description of the content of each file and directory in the root PLUMED directory.

    CHANGES          : change log
    COPYING.LESSER   : license
    Makefile         : makefile
    Makefile.conf.in : template configuration makefile
    PEOPLE           : list of authors
    README.md        : this file
    VERSION          : version file
    astyle           : a local version of astyle, used to format code
    configure        : configuration script
    configure.ac     : configuration script (autoconf)
    developer-doc    : developer documentation
    docker           : directory where Docker is generated
    macports         : directory where Portfiles are generated
    patches          : patch scripts
    python           : python stuff
    regtest          : regression tests, including reference results
    release.sh       : developer utility to publish releases
    scripts          : shell tools
    sourceme.sh.in   : template configuration script
    src              : source code
    test             : examples
    user-doc         : user documentation
    vim              : directory where vim syntax is generated

Required software
-----------------

Required software:

* GNU make.
* C/c++ compiler (c++11 support is required as of version 2.4).
* A modern version of the `patch` command line tool.
* Support for POSIX library `dirent.h`.
* `xxd` (present in most UNIX distributions).

Suggested software (libraries are checked by `./configure` and enabled if available):

* MPI library to run parallel simulations. It should be the same library used by your MD code.
* Optimized blas and lapack libraries. They are automatically replaced by an internal version if not available.
* [VMD molfile plugins](http://www.ks.uiuc.edu/Research/vmd/plugins) to read arbitrary file formats. They are automatically replaced by an internal version supporting a few formats if not available.
* [Zlib library](http://zlib.net/) to use compressed data files.
* [Xdrfile library](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library) to have read/write access to gromacs
  trajectory files.
* [Doxygen](http:://www.doxygen.org) to build user manual. Doxygen might need the following packages:
  * Latex to build the pdf user manual.
  * [Graphviz](http://www.graphviz.org) to show class hierarchy in
    developer manual.

Quick compilation instructions
------------------------------

Extensive installation instructions are in the [user documentation](http://www.plumed.org/documentation).
Quick instructions:

    ./configure --prefix=$HOME/opt
    make
    make doc # optional
    make test # optional

User documentation can be found at `user-doc/html/index.html`.
Developer documentation can be found at `developer-doc/html/index.html`.
[Pre-compiled documentation](http://www.plumed.org/documentation) is available online, so this is only required
if you are working with a modified version of the code!

In order to run PLUMED without installing it you should type `source sourceme.sh`. However,
we recomment installing PLUMED. 
To install it in `$HOME/opt` (directory should be set during `./configure`):

    umask 022
    make install
    
Now you will be able to run plumed using e.g.

    plumed help

If you compiled your own documentation, paths to the installed documentation can be found with command `plumed info --user-doc`.

A sample modulefile with environment variable will be placed in
`$HOME/opt/lib/plumed/src/lib/modulefile`. This can be useful if you want to
install multiple PLUMED versions side by side and select them with env modules.
Instructions for using Artistic Style are included in the *doc* directory.

The file **install.html** contains instructions for compiling and
installing Artistic Style.

The file **astyle.html**' contains information on using Artistic Style.

The files **news.html** and **notes.html** contain information on changes
made to the various releases.
ANN (Artificial Neural Network) function for plumed
====================

This is plumed ANN function (annfunc) module.  It implements `ANN` class, which is a subclass of `Function` class.  `ANN` class takes multi-dimensional arrays as inputs for a fully-connected feedforward neural network with specified neural network weights and generates corresponding outputs.  The `ANN` outputs can be used as collective variables, inputs for other collective variables, or inputs for data analysis tools.  

## Installation

Enable compilation by adding the `--enable-modules=annfunc` to the configure command.

## Usage

It is used in a similar way to [other plumed functions](https://www.plumed.org/doc-v2.5/user-doc/html/_function.html).  To define an `ANN` function object, we need to define following keywords:

- `ARG` (string array): input variable names for the fully-connected feedforward neural network

- `NUM_LAYERS` (int): number of layers for the neural network

- `NUM_NODES` (int array): number of nodes in all layers of the neural network

- `ACTIVATIONS` (string array): types of activation functions of layers, currently we have implemented "Linear", "Tanh", "Circular" layers, it should be straightforward to add other types as well

- `WEIGHTS` (numbered keyword, double array): this is a numbered keyword, `WEIGHTS0` represents flattened weight array connecting layer 0 and layer 1, `WEIGHTS1` represents flattened weight array connecting layer 1 and layer 2, ...  An example is given in the next section.

- `BIASES` (numbered keyword, double array): this is a numbered keyword, BIASES0 represents bias array for layer 1, BIASES1 represents bias array for layer 2, ...

Assuming we have an `ANN` function object named `ann`, we use `ann.node-0, ann.node-1, ...` to access component 0, 1, ... of its outputs (used as collective variables, inputs for other collective variables, or data analysis tools).

## Examples

Assume we have an ANN with numbers of nodes being [2, 3, 1], and weights connecting layer 0 and 1 are

```
[[1,2],
[3,4],
[5,6]]
```

weights connecting layer 1 and 2 are

```
[[7,8,9]]
```

Bias for layer 1 and 2 are

```
[10, 11, 12]
```

and 

```
[13]
```

respectively.

All activation functions are `Tanh`.

Then if input variables are `l_0_out_0, l_0_out_1`, the corresponding `ANN` function object can be defined using following plumed script: 

```
ann: ANN ARG=l_0_out_0,l_0_out_1 NUM_LAYERS=3 NUM_NODES=2,3,1 ACTIVATIONS=Tanh,Tanh  WEIGHTS0=1,2,3,4,5,6 WEIGHTS1=7,8,9  BIASES0=10,11,12 BIASES1=13
```

This plumed script can be generated with function `Plumed_helper.get_ANN_expression()` in [this](https://github.com/weiHelloWorld/plumed_helper/blob/master/plumed_helper.py) repository.  Following is the Python code using this function to generate the script above:

```Python
from plumed_helper import Plumed_helper
ANN_weights = [np.array([1,2,3,4,5,6]), np.array([7,8,9])]
ANN_bias = [np.array([10, 11, 12]), np.array([13])]
Plumed_helper.get_ANN_expression('ANN', node_num=[2, 3, 1], 
                                 ANN_weights=ANN_weights, ANN_bias=ANN_bias,
                                 activation_list=['Tanh', 'Tanh'])
```

## Authors

Wei Chen (UIUC, weichen9@illinois.edu) and Andrew Ferguson (University of Chicago, andrewferguson@uchicago.edu)

## Copyright

See ./COPYRIGHT
maze
================================================================================
This version of PLUMED2 has the maze module implemented.

Install
--------------------------------------------------------------------------------
Enable the compilation of maze by adding the `--enable-modules=maze` to the 
configure command.

Page
--------------------------------------------------------------------------------
See [this link](http://maze-code.github.io) for further information.

Documentation
--------------------------------------------------------------------------------
Run `make doc`; the documentation should be in `user-doc/html/_m_a_z_e.html`.

Author
--------------------------------------------------------------------------------
Jakub Rydzewski (Nicolaus Copernicus University) <jr@fizyka.umk.pl>

Copyright
--------------------------------------------------------------------------------
This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU Lesser General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any 
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  

See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along 
with this program.  If not, see <http://www.gnu.org/licenses/>.
Welcome to the plumed2 wiki!

This fork of PLUMED has eABF/DRR implementation.

**Requirements**

Boost::serialization and C++11 compiler

**Compiling instruction:**

After clone this repository, please cd to the plumed2 directory and run:

1. `autoconf`
2. `./configure --enable-boost_serialization --enable-modules=drr`
3. Modify your Makefile.conf and add `-lboost_serialization` in `DYNAMIC_LIBS=`
4. `make` and `sudo make install`

**Usage**

Run `make doc` and the usage is in `user-doc/html/_d_r_r.html`

**Authors**

Chen Haochuan (All files in drr module except colvar_UIestimator.h)

Fu Haohao (colvar_UIestimator.h)

**COPYRIGHT**

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Experiment Directed Simulation (EDS)
====================================


Install
------------------------------------
Enable compilation by adding the `--enable-modules=+eds`
to the configure command.


Documentation
------------------------------------
See the generated documentation for information on
using


Authors
------------------------------------
Glen Hocky (University of Chicago) <hockyg@uchicago.edu>
Andrew White (University of Rochester) <andrew.white@rochester.edu>


Copyright
------------------------------------
See ./COPYRIGHT
@page CHANGES-2-6 Version 2.6
  
## Version 2.6 (Jan 27, 2020)

Changes from version 2.5 which are relevant for users:
- Changes leading to incompatible behavior:
  - PLUMED input file parsing is now case insensitive that is that all directives can be written using uppercase characters (compatible with former versions) as well as lowercase characters (not compatible) internally PLUMED still uses uppercase definitions
  - `plumed partial_tempering` now uses `gawk` instead of `awk`. You might need to install `gawk` for it to work correctly.

- Other changes:
  - Asmjit is now embedded into PLUMED. In order to enable it, it is sufficient to configure with `--enable-asmjit`. See \ref Lepton "this page".
  - Fixed grids so as to decrease memory footprint of derivatives (see \issue{465}).
  - Added option `--idlp4` to \ref driver to read DLPOLY4 HISTORY files (see \issue{478}, thanks to Alin Marin Elena).
  - Added atom selectors using mdtraj/MDAnalysis/VMD syntax, see \ref MOLINFO and \issue{448}.
  - \ref EEFSOLV is now faster in scalar and also mpi/openmp parallel
  - New shortcuts are available for selecting protein atoms: `@sidechain-#`, `@back-#`
  - VIM syntax highlight is now case insensitive. Notice that autocompletion still only works with upper case commands.

- New contributed modules:
  - A new Maze module by Jakub Rydzewski
     - \ref MAZE_LOSS
     - \ref MAZE_MEMETIC_SAMPLING
     - \ref MAZE_RANDOM_ACCELERATION_MD
     - \ref MAZE_RANDOM_WALK
     - \ref MAZE_SIMULATED_ANNEALING
     - \ref MAZE_STEERED_MD
     - \ref MAZE_OPTIMIZER_BIAS
  - A new ANN module by Wei Chen and Andrew Ferguson
     - \ref ANN

- New patches:
  - added support for AMBER PMEMD 18 (contributed by Viktor Drobot, see \issue{486}).

- Changes in the VES module
  - new \ref VES_DELTA_F bias.
  - ves_md_linearexpansion now outputs one-dimensional free energy projections of the potential energy landscape. 

- Changes in the DRR module
  - The MAXFACTOR option now is tunable for each CV in multidimensional cases.
  - Output .zcount file (the same as .czar.count) for compatibility with newer abf_integrate.
  - The citation of DRR module has been updated.

- Changes in the ISDB module
  - in \ref METAINFERENCE we removed the MC_STRIDE keyword
  - in \ref METAINFERENCE the bias value (metainference score) now includes the Jeffrey's prior (values are different, but forces are equal)
  - components were previously named using _ but now they abide to the standard is -
  - removed ADDEXP keywords for \ref JCOUPLING \ref NOE \ref PRE \ref RDC
  - \ref METAINFERENCE performs more check on the input and restart files to ensure a consistent setup
  - \ref SAXS is slightly faster and scales better, removed BESSEL options

- Python module:
  - Removed compatibility with Python 2.
  - Added capability to read and write pandas dataset from PLUMED files (see \issue{496}).

Changes from version 2.5 which are relevant for developers:
  - Components documentation is now enforced
  - `readdir_r` is deprecated and is thus not used by default (can be enabled with `./configure --enable-readdir-r`).

## Version 2.6.1

For users:
- New patches:
  - added gromacs 2019.6 
  - added gromacs 2020.1 (experimental) 

For developers:
- Small fix to avoid unique global symbols (see \issue{549})

@page CHANGES-2-1 Version 2.1

Version 2.1.0 (September 15, 2014)
----------------------------

Version 2.1 contains several improvements with respect to 2.0. Users currently working with 2.0 
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files. In 2.1 we restored more features of 1.3
that were missing in 2.0, so users still working with 1.3 could opt for an upgrade.
A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Below you find a list of all the changes with respect to version 2.0.
Notice that version 2.1 includes already all the fixes in branch 2.0 up to 2.0.4.

Changes from version 2.0 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref COORDINATION now skips pairs of one atom with itself.
  - Labels of quantities calculated by \ref BIASVALUE have changed from <i>label</i>.bias.<i>argname</i> to <i>label</i>.<i>argname</i>_bias, which is more consistent with steered MD
  - Labels of quantities calculated by \ref ABMD have change from <i>label</i>.min_<i>argname</i> to <i>label</i>.<i>argname</i>_min, which is more consistent with steered MD
  - Labels of quantities calculated by \ref PIECEWISE have change from <i>label</i>.<i>argnumber</i> to <i>label</i>.<i>argname</i>_pfunc, which is more consistent with steered MD
  - For multicolvars components calculated with LESS_THAN and MORE_THAN keywords are now labelled lessthan and morethan. This change is necessary as the underscore
character now has a special usage in component names.
  - In \ref CONTACTMAP components are now labelled <i>label</i>.contact-\f$n\f$.
  - The command SPHERE has been replaced by \ref UWALLS.
- New configuration system based on autoconf (use ./configure from root directory).
  Optional packages are detected at compile time and correctly
  enabled or disabled. An internal version of LAPACK and BLAS will be used
  if these libraries are not installed.
- New actions:
  - \ref SPRINT topological collective variables.
  - CH3SHIFTS collective variable.
  - \ref POSITION collective variable.
  - \ref FIT_TO_TEMPLATE.
  - \ref COMMITTOR analysis.
  - \ref LOCAL_AVERAGE.
  - \ref NLINKS.
  - \ref DIHCOR.
  - \ref NOE.
  - \ref RDC. 
  - \ref CLASSICAL_MDS.
  - \ref XDISTANCES.
  - \ref YDISTANCES.
  - \ref ZDISTANCES.
  - \ref DUMPMULTICOLVAR.
  - Crystallization module, including \ref Q3, \ref LOCAL_Q3, \ref Q4, \ref Q6, \ref LOCAL_Q4, \ref LOCAL_Q6, \ref MOLECULES, \ref SIMPLECUBIC, \ref TETRAHEDRAL and \ref FCCUBIC.
  - \ref ENSEMBLE to perform Replica-Averaging on any collective variable.
- New features for existing actions:
  - \ref METAD : WALKERS_MPI flag (multiple walkers in a mpi-based multi-replica framework),
    ACCELERATION flag (calculate on the fly the Metadynamics acceleration factor),
    TAU option (alternative way to set Gaussian height in well-tempered metadynamics),
    GRID_SPACING (alternative to GRID_BIN to set grid spacing).
    Notice that now one can also omit GRID_BIN and GRID_SPACING when using
    fixed size Gaussian, and the grid spacing will be automatically set.
  - \ref DISTANCE : added SCALED_COMPONENTS
  - \ref COORDINATION : if a single group is provided, it avoids permuted atom indexes and runs
    at twice the speed.
  - \ref DUMPATOMS : PRECISION option to set number of digits in output file.
  - \ref GROUP : NDX_FILE and NDX_GROUP options to import atom lists from ndx (gromacs) files.
  - In many multicolvars, MIN and MAX options can be used.
  - \ref HISTOGRAM : GRID_SPACING (alternative to GRID_BIN to set grid spacing),
    FREE-ENERGY flags in addition to standard probability density,
    additional option for KERNEL=DISCRETE to accumulate standard histograms. 
  - \ref sum_hills : added options --spacing (alternative to --bin to set grid spacing)
    and --setmintozero to translate the minimum of the output files to zero.
  - \ref CONTACTMAP : parallelized and added weights.
- New features in MD patches (require re-patch):
  - New patch for Gromacs 5.0
  - Gromacs 4.6.X patch updated to 4.6.7
  - Gromacs 4.6.7 supports \ref COMMITTOR analysis; can be now be used to perform energy minimization;
     now passes temperature to PLUMED (this allows temperature to be omitted in some actions,
     namely \ref METAD and analysis actions).
  .
  Notice that if you use runtime binding it is not compulsory to re-patch,
  and that all combinations should work correctly
  (new/old PLUMED with re-patched/non-re-patched MD code).
- Other new features:
  - \ref driver can now read trajectories in many formats using VMD molfile plugin
    (requires VMD plugins to be compiled and installed). In case VMD plugins are not installed,
    the configuration system falls back to an internal version which implements a minimal list
    of plugins (gromacs and dcd) (kindly provided by T. Giorgino).
  - \ref switchingfunction : added STRETCH flag.
  - Negative strides in atom ranges (e.g. ATOMS=10-1:-3 is expanded to ATOMS=10,7,4,1).
  - \ref COORDINATION and \ref DHENERGY with NLIST now work correctly in replica exchange simulations.
  - Multicolvars with neighbor lists now work correctly in replica exchange simulations.
  - Improved multicolvar neighbor lists.
- Optimization:
  - Root-mean-square deviations with align weights different from displace weights
    are now considerably faster. This will affect \ref RMSD calculations plus
    other variables based on RMSD.
  - \ref WHOLEMOLECULES is slightly faster.
  - \ref COORDINATION is slightly faster when NN and MM are even and D_0=0.
  - Atom scattering with domain decomposition is slightly faster.
  - Link cells are now exploited in some multicolvars.
  - Derivatives are not calculated unless they are specifically required, because for instance you are adding
    a bias.
- Documentation:
  - All tutorial material from the recent plumed meeting in Belfast is now in the manual
  - Improvements to documentation, including lists of quantities that are output by each action that can be referenced 
  - Manual has been re-organized following suggestions received at the plumed meeting.
  - An experimental PDF version of the manual is now provided (a link can be found in the documentation homepage).

Changes from version 2.0 which are relevant for developers:
- Added regtests for plumed as a library (e.g. basic/rt-make-0). plumed command has an additional
  flag (--is-installed) to probe if running from a compilation directory or from a fully installed copy
  (this is needed for regtests to work properly).
- Improved class Communicator. Many operations can now be done directly on Vectors, Tensors, std::vector and PLMD::Matrix.
- Modified class RMSD.
- Patches for GPL codes (Quantum Espresso and Gromacs) now also include
  original code so as to simplify their modification.
- Fixed dependencies among actions such that it is now possible (and reliable)
  to use MPI calls inside Action::prepare()
- colvar/CoordinationBase.cpp has been changed to make it faster. If you devised a class which inherits from here,
  consider that CoordinationBase::pairing now needs _squared_ distance instead of distance
- It is possible to run "make install" from sub-directories (e.g. from src/colvar)
- There is a small script which disables/enables all optional modules (make mod-light/mod-heavy/mod-reset)
- Added "-q" option to plumed patch
- You can now create new metrics to measure distances from a reference configurations. If you do so such
  metrics can then be used in paths straightforwardly
- You can now use multicolvars in tandem with manyrestraints in order to add a large numbers of restraints.
- Can now do multicolvar like things in which each colvar is a vector rather than a scalar.
- Updated script that generated header files so that they properly show years. Notice that the script
  should new be run from within a git repository

This list is likely incomplete, if you are developing in PLUMED you are encouraged to follow changes on github.

Version 2.1.1 (December 15, 2014)
----------------------------------------------

This release includes all the fixes available in branch 2.0 until 2.0.5.

For users:
- New patch for AMBER 14 (sander module only). This patch should be compatible
  with any PLUMED 2 version (including 2.0). It includes most PLUMED features
  with the notable exception of multi-replica framework.
- Changed definition in arbitrary phase of eigenvectors. This will change the result of some
  analysis method where the phase does matter (e.g. \ref CLASSICAL_MDS) and make
  some regression test better reproducible.
- Fixed a portability issue in BG/P where gettimeofday is not implemented.
  Notice that this fix implies that one should execute again ./configure to have
  plumed timing working correctly.
- CS2Backbone: fixed a bug that resulted in only a fraction of the chemical shifts being printed with WRITE_CS and 
  parallel simulations (requires to get the last almost updated from SVN)
- NOE: fixed a bug in the replica-averaging
- Fixed a linking issue with ALMOST, where bz2 was always used to link ALMOST to PLUMED even if it is not compulsory 
  to build ALMOST.
- Fixed a wrong include in the GMX5 patch.
- \ref FUNCPATHMSD can now be used together with \ref CONTACTMAP to define pathways in contact-map space
- Configuration is more verbose, a warning is given if a default option cannot be enabled and an error is given if 
  an option explicitly enabled cannot be enabled.
- Compilation is less verbose (use "make VERBOSE=1" to have old behavior)
- Small fixes in documentation.

For developers:
- Tests are now performed at every single push on travis-ci.org
- Manual is built and pushed to the online server from travis-ci.org (see developer doc)
- Fixes in developer doc.

Version 2.1.2 (Mar 16, 2015)
----------------------------------------------

For users:
- Added two new short tutorials to the manual ( \ref cambridge and \ref munster ).
- Fixed a severe bug on \ref DRMSD - cutoff values were ignored by PLUMED.
  Notice that this bug was introduced in 2.1.0, so that it should not affect the 2.0.x series.
- Fixed a bug affecting LAMMPS patch used with a single processor. Notice that
  the fix is inside PLUMED, thus it does not necessarily requires re-patching.
- Sander patch now works with multiple replica (no replica exchange yet). It also contains
  some fix from J. Swails.
- GMX5 patch was not working for bias-exchange like cases
- Patching system now checks for the availability of shared/static/runtime version of plumed before
  patching
- Configure now check better if compiler flag are accepted by the compiler. This makes
  configure on bluegene more robust.
- Sourceme.sh now sets proper library path in linux also.


Version 2.1.3 (June 30, 2015)
----------------------------------------------

For users:
- Fixed bug in \ref ENSEMBLE derivatives when more than 1 argument was provided
- Fixed bug in \ref GHOST : virial is now computed correctly.
- Fixed a serious bug in virial communicated from plumed to gromacs, for both gromacs versions 4.6 and 5.0.
  See \issue{132}.
  This fix requires gromacs to be re-patched and could be very important if you run biased simulations in the NPT ensemble.
- Fixed a bug in the virial computed with \ref FIT_TO_TEMPLATE when the reference pdb had center non located at the origin.
- Fixed a bug in the the forces computed with \ref FIT_TO_TEMPLATE when used in combination with \ref COM, \ref CENTER, or \ref GHOST
- Fixed a bug that could lead plumed to be stuck with domain decomposition in some extreme case (one domain with all atoms, other domains empty).
- Fixed a bug when \ref COMBINE or \ref MATHEVAL are used with PERIODIC keyword. Now when PERIODIC keyword is used the result
  of the calculation is brought within the periodicity domain. See \issue{139}.
- Fixed a bug related to \ref RANDOM_EXCHANGES followed by \ref INCLUDE
- Fixed bug in derivatives of histogram bead with triangular kernels
- Updated gromacs patch 4.5.5 to 4.5.7
- Updated internal molfile plugins to VMD 1.9.2.
- Included crd and crdbox formats to internal molfile.
- Added --natoms to \ref driver . This is required to read coordinate
  files with VMD plugins when number of atoms is not present (e.g. amber
  crd files)
- Added the checks in the driver to detect cases where molinfo does not provide box information
  (e.g. pdb).
- Added support for readdir_r when available, which makes opening files thread safe.
- CFLAGS now include -fPIC by default
- Added a warning when using \ref METAD without grids with a large number of hills.
- Fixes in user documentation.

For developers:
- Allow external VMD plugins to be detected with --has-external-molfile. This
  is required to enable some regtest with amber files.
- Added --dump-full-virial to \ref driver
- Allow definition of variables where some of the components have derivatives and some haven't (\issue{131}).
- Improved travis tests with more debug options.
- Improved some regtest to check out-of-diagonal virial components
- Improved make cppcheck options.
- Fixes in developer documentation.

Version 2.1.4 (Oct 13, 2015)
-----------------------------

For users:
- Fixed NAMD patch. Masses and charges were not passed correctly, thus resulting in wrong
  \ref COM or \ref CENTER with MASS.
  This fix required re-patching NAMD.
  Notice that this bug was present also in v2.0 but in a different form.
  More information here (\issue{162}), including a workaround that allows masses to be fixed
  without re-patching.
- When installing with PLUMED_LIBSUFFIX an underscore is used as separator instead of a dash.
  E.g. `make install PLUMED_LIBSUFFIX=2.1` will result in an executable named `plumed_v2.1`.
  This fix a potential problem (see \ref Installation).
- Fixed erroneously reported message about MPI at the end of ./configure.
- Changed warning message about undocumented components.
- PLUMED now says in the log file if it was compiled from a dirty git repository.
- Fixed a problem leading to rare random crashes when using \ref METAD with WALKERS_MPI and multiple
  processors per replica.
- Small change in numerical accuracy of lattice reduction. Should be more
  robust when running with highly optimizing compilers.
- Fixed a bug in normalization of kernel functions.  This affects \ref HISTOGRAM
  If these actions were used with previous versions of the code care should be taken when analyzing the 
  results.
- Fixed a bug in derivatives of kernel functions with non-diagonal covariance matrices. This affects the 
  derivatives output by \ref sum_hills

Version 2.1.5 (Jan 18, 2016)
---------------------------------------------

\plumednotmaintained

For users:
- PLUMED now reports an error when using \ref HISTOGRAM with FREE-ENERGY without USE_ALL_DATA. See \issue{175}
- Fixed a bug in configure together with --enable-almost. The check for lbz2 library was not working properly.

@page CHANGES-2-5 Version 2.5

## Version 2.5 (Dec 19, 2018)

This page contains changes that will end up in 2.5

Changes from version 2.4 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref RMSD, \ref MULTI-RMSD, \ref PATHMSD, \ref PROPERTYMAP, \ref PCAVARS, \ref PCARMSD, \ref FIT_TO_TEMPLATE,
    \ref DIPOLE, \ref ALPHARMSD, \ref ANTIBETARMSD, and \ref PARABETARMSD now automatically make molecules whole.
    In case you do not want them to do it, use NOPBC flag,
  - There is some subtle change in the installation layout (see below). There should be no visible effect, however it is now compulsory
    to set correctly the `LD_LIBRARY_PATH` variable for the linux executable to work correctly. The procedure has been tested well on OSX and Linux,
    but could give problems on other platform. Please report possible problems on the mailing list.
  - \ref driver now stops correctly when using \ref COMMITTOR. If you want to continue the analysis, use the `NOSTOP` flag in \ref COMMITTOR.
  - \ref METAD the calculation of the reweighting factor is now activated by CALC_RCT instead of REWEIGHTING_NGRID and REWEIGHTING_NHILLS, the frequency of update can be set 
    by RCT_USTRIDE, the default value is 1 and should be OK for most of the cases
  - Fixed sign in Cartesian components of \ref PUCKERING with 6 membered rings (thanks to Carol Simoes and Javi Iglesias).

- New actions:
  - \ref COLLECT_FRAMES
  - \ref EUCLIDEAN_DISSIMILARITIES
  - \ref HBPAMM_MATRIX
  - \ref HBPAMM_SH
  - \ref LANDMARK_SELECT_FPS
  - \ref LANDMARK_SELECT_RANDOM
  - \ref LANDMARK_SELECT_STAGED
  - \ref LANDMARK_SELECT_STRIDE
  - \ref OUTPUT_ANALYSIS_DATA_TO_COLVAR
  - \ref OUTPUT_ANALYSIS_DATA_TO_PDB
  - \ref OUTPUT_PCA_PROJECTION
  - \ref PAMM
  - \ref PLUMED
  - \ref PRINT_DISSIMILARITY_MATRIX
  - \ref PROJECT_ALL_ANALYSIS_DATA
  - \ref READ_DISSIMILARITY_MATRIX
  - \ref RESELECT_LANDMARKS
  - \ref REWEIGHT_WHAM
  - \ref SKETCHMAP_CONJGRAD
  - \ref SKETCHMAP_POINTWISE
  - \ref SKETCHMAP_READ
  - \ref SKETCHMAP_SMACOF
  - \ref SKETCH_MAP
  - \ref SMACOF_MDS
  - \ref WHAM_HISTOGRAM
  - \ref WHAM_WEIGHTS

- New command line tools:
  - \ref completion (used to generate command line completion scripts).
  - \ref pdbrenumber (see \issue{371}).

- New modules:
  - A new PIV module has been included, contributed by Silvio Pipolo and Fabio Pietrucci.
    This module implements the following collective variable:
    - \ref PIV
  - A new LOGMFD module has been included, contributed by Tetsuya Morishita.
    This module implements the following bias:
    - \ref LOGMFD

- Changes in the ISDB module
  - \ref CS2BACKBONE is now mpi parallelized in particular with DOSCORE and CAMSHIFT
  - \ref SAXS has an additional implementation based on Bessel functions that can be faster for large systems (new keyword BESSEL)
  - \ref SAXS keyword SCEXP has been renamed into SCALEINT
  - \ref SAXS includes the MARTINI bead structure factors for Proteins and Nucleic Acids
  - \ref SAXS includes a GPU implementation based on ArrayFire (need to be linked at compile time) that can be activated with GPU
  - \ref METAINFERENCE and all related methods has a new keyword REGRES_ZERO to scale data using a linear scale fit
  - \ref CALIBER new bias to perform Maximum Caliber replica-averaged restrained simulations 

- Changes in the eABF/DRR module (contributed by Haochuan Chen and Haohao Fu):
  - \ref DRR now supports the extended generalized ABF(egABF) method.
  - \ref DRR accepts different GRID options for CVs and extended variables.
  - The MAXFACTOR option is added in \ref DRR to control the factor of biasing force.
  - \ref drr_tool can calculate the divergence of gradients now. (Maybe useful for future pABF)
  - Fixed conflicts of output files in multiple replicas.

- Changes in the EDS module:
  - \ref EDS implements Levenberg-Marquardt optimization in addition to previous gradient descent. 
  - \ref EDS no longer automatically increases prefactor for bias parameter updates. This results in more stable optimization for the cases tested.
  - \ref EDS now has a larger default RANGE parameter to go with these other changes.

- Other changes:
  - \ref METAD there is a new FLYING_GAUSSIAN keyword to activate the flying gaussian methods by Spiwok (contributed by Spiwok and Hozzova)
  - \ref EXTERNAL can now SCALE the input grid. This allows for more flexibility without modifying the grid file.
  - \ref ALPHABETA can now combine dihedral angles with different coefficients
  - \ref INCLUDE can now be used also before setup actions.
  - \ref CENTER can now be computed using trigonometric functions (PHASES) to simplify its calculation with periodic boundary conditions.
  - Libmatheval is not used anymore. \ref MATHEVAL (and \ref CUSTOM) are still available
    but employ an internal implementation of the lepton library.
    Functions available in libmatheval and absent in the original lepton library have been added so as to have backward compatibility.
    `atan2(y,x)` function has also been added.
    Notice that MATHEVAL (and CUSTOM) \ref switchingfunction "switching functions"
    using the lepton library have been further optimized with respect to PLUMED 2.4.
    Finally, notice that it is possible to use asmjit to optimize performance (see \ref Lepton).
  - Implemented bash autocompletion, see \ref BashAutocompletion.
  - \ref MOLINFO now allows selecting atoms from chains with a numeric ID (see \issue{320}).
  - Removed the patch for GMX 5.1.4
  - LAMMPS patch has been finally removed. Notice that LAMMPS has native support for PLUMED now.
  - AMBER patch has been finally removed. Notice that AMBER (sander module) has native support for PLUMED starting from version 15.
  - \ref RMSD calculation has been optimized. This should positively affect the performances of CVs where
     many RMSD values are computed on small groups of atoms, such as secondary structure variables.
  - In \ref METAD, when using a bias factor equal to one (no bias) the `rct` component is set to zero rather than to one.
  - New shortcuts are available for selecting atoms: `@allatoms` and `@mdatoms` (see \ref atomSpecs).
  - When using \ref MOLINFO, also the following shortcuts are available for selecting atoms: `@nucleic`, `@protein`, `@water`, `@ions`, `@hydrogens`, `@nonhydrogens`.
  - When using \ref MOLINFO, individual atoms can be chosen also from water molecules (e.g. `@OW-100`).
  - Additional switching function COSINUS contributed by Michael King
  - added API to set the number of used openMP threads from the linked code, updated gromacs 2018.3 patch to use it

Changes from version 2.4 which are relevant for developers:
- Code has been cleanup up replacing a number of pointers with `std::unique_ptr`. All `delete` statements
  in the core parts of the code have been eliminated.
- Exceptions cannot be disabled (`--disable-cxx-exceptions` option has been removed from `./configure`).
- Every exception thrown in PLUMED now also writes its message on PLUMED log.
- Runtime loader in `Plumed.c` now works also when linked without `-rdynamic` (that is, 
  its names are not exported). Notice that all the combinations are expected to
  work, that is: `Plumed.c` from <=2.4 or >=2.5 combined with libplumedKernel
  from <=2.4 or >=2.5. In order to achieve this the following changes are implemented:
  - libplumedKernel does not depend anymore on `Plumed.c`. This allows loading it even
    in cases where names in the loader are not visible. The relevant function needed
    to be compatible with `Plumed.c` <=2.4 are found using `dlsym`.
  - `Plumed.c` does not need anymore libplumedKernel to register itself, but rather
    searches the relevant functions using `dlsym`. In addition, if it is not able to
    load `libplumedKernel` since the latter is <=2.4 and needs `Plumed.c` to be visible,
    it just uses as a fallback `libplumed`, which should load properly.
- In addition to the capability mentioned above, the MD-code interface has been significantly
  improved and allows for:
  - Translation of exception (allowing to mix PLUMED and an MD-code linked against a different C++ library).
  - Possibility to choose the path to the PLUMED kernel while instantiating a Plumed object.
  See the developer documentation for more information.
- The installation layout of shared libraries has been modified. In particular,
  both `libplumed.so` and `plumed` links to `libplumedKernel.so`.
  This reduces considerably the size of the installed package. In addition, it allows
  using two-level namespace on OSX. Notice that this implies that on Linux one should
  always set the `LD_LIBRARY_PATH` flag to have a working executable.
- A smaller number of header files is installed. In particular, all the files that were historically generated in subdirectories
  (such as `plumed/core/tools/Vector.h', just including `plumed/tools/Vector.h`) are not installed and the related include
  statements are fixed. This makes the installed package smaller.
- List of preferred compilers (used when `CXX` or `CC` are not set) has been changed. On OSX, `./configure` will try `clang++/clang` as first choices.
- Added `--enable-static-archive` to `./configure` to build a `libplumed.a` static library (yes by default).
- Stop setting `DYLD_LIBRARY_PATH` in `sourceme.sh` and in modulefile. Notice that as of PLUMED v2.3.3
  it should not be needed.
- Coverage scan is not anymore contained in developer manual. It can be found in a separate repository
  `github.com/coverage-branch` (see \issue{348}). In addition, coverage for third-party libraries included in PLUMED
  is reported as well.
- It is not possible anymore to use `make install prefix=/path`. Prefix can only be changed during `./configure` (see \issue{332}).
- Exception class has been rewritten to allow more extensive messages. Now also function name is shown.
- On linux, library is linked with `-Bsymbolic`.
- When launching `plumed`, flags `--no-mpi` and `--mpi` can appear multiple times. The last appearance is the effective one.
- Internal BLAS and LAPACK libraries updated to gromacs 2018.
- Choosing `./configure --prefix=$PWD` does not lead anymore to deletion of all header files.
- A copy of `plumed-runtime` is installed in `prefix/lib/plumed` and can be used for testing.
- Absolute/relative soname/install_name can be configured on linux/OSX. This feature is only
  for testing, the default choice is the typical one used on the respective operating system.
- On OSX, `plumed` and `libplumed.dylib` will find `libplumedKernel.dylib` using `@loader_path`.
- Using CXX compiler to link the main program.
- plumed can be compiled with ArrayFire to enable for gpu code. \ref SAXS collective variable is available as part of the isdb module to provide an example of a gpu implementation for a CV


## Version 2.5.1 (Apr 1, 2019)

For users:
- in \ref SAXS the keyword ADDEXP is removed. Furthemore, SAXS intensities are automatically normalised for I(0)=1, in case experimental data are provided, the intensity is rescaled with the intensity of the lowest q provided. As a consequence SCALEINT is only needed for additional adjustments.
- gromacs patch updated to gromacs 2018.5
- Fixed a bug in gromacs patch that was resulting in incorrect number of threads (0) set when not explicitly using `-ntomp` on the 
  command line or setting `OMP_NUM_THREADS` (see \issue{446}). To apply this fix you need to re-patch gromacs.
  Notice that setting the number of threads to zero might lead to inconsistent results when using secondary structure variables
  or other multicolvars.
- Fixed PLUMED so that when zero threads are selected from gromacs (see previous fix) the number of used threads is set to 1.
  This fix allows to use a GROMACS executable patched with PLUMED 2.5.0 and linked at runtime with PLUMED 2.5.1 without introducing
  errors. However, re-patching is preferred since it selectes the correct number of threads.
- Python wrappers:
  - Fixed building of python interface on MacOS Mojave (see \issue{445}, thanks to Omar Valsson).
  - Numpy is not required anymore at build time (though it is required at runtime for our tests).
  - Raw python arrays can be passed as an alternative to Numpy ndarrays.

## Version 2.5.2 (Jul 19, 2019)

For users:
- New shortcuts are available for selecting protein atoms: `@chi2-#`, `@chi3-#`,`@chi4-#` and `@chi5-#`
- Fixed performance of \ref CUSTOM when having zero derivatives with respect to some arguments.
- New --parse-only option in \ref driver to check the validity of a plumed input file
- New patch for GROMACS 2019.2
- Module VES: Fixed performance of \ref BF_CUSTOM for basis functions with linear terms (e.g. having zero derivatives). 
- Python wrappers:
  - Python module is now always named `plumed` irrespectively of program prefix and suffix. Notice 
    that python module is installed inside the `lib/program_name` directory and thus it is not necessary to
    use `program_name` in order to install multiple modules side by side.
  - Python module can be compiled without compiling PLUMED first.
  - `Plumed` object can be explicitly finalized using `finalize()`. Can be used to make sure all files are closed,
    but it is not necessary if the `Plumed` object gets correctly collected by Python.
  - `Plumed` object can be used in context managers (e.g. `with plumed.Plumed() as p:`).
- Precompiled binaries are available on Anaconda cloud on the [conda-forge channel](https://anaconda.org/conda-forge/plumed).

## Version 2.5.3 (Oct 11, 2019)

For users:
- Fixed a bug with \ref CONVERT_TO_FES and periodic variables, see \issue{441}
- Fixed a bug with \ref FOURIER_TRANSFORM 
- Updated patch for GROMACS 2019.4
- Updated patch for GROMACS 2018.8
- Python module:
  - Fixed building with clang-8.
  - Set `language_level` for cython to the actually used language level.
  - Force using cython when compiling from source. Still using the pre-generated cpp file
    when installing from PyPI, to avoid cython dependency.
  - Using python 2 to create the cpp file uploaded on PyPI (this will change to python 3 in 2.6, see \issue{502}).
- Module VES: Fixed a bug in updating of bias potential in \ref VES_LINEAR_EXPANSION that is present for certain integrators that call the calculation of the bias multiple times (see [here](https://groups.google.com/d/msg/plumed-users/kPZu_tNZtgk/LrkS0EqrCQAJ)) and replica exchange.

## Version 2.5.4 (Jan 27, 2020)

For users:
- Includes all fixes up to 2.4.7

## Version 2.5.5

For developers:
- Small fix to avoid unique global symbols (see \issue{549})

@page CHANGES-2-2 Version 2.2

Version 2.2 (Oct 13, 2015)
----------------------------

Version 2.2 contains several improvements with respect to 2.1. Users currently working with 2.1
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files. In 2.2 we restored more features of 1.3
that were missing in 2.1, so users still working with 1.3 could opt for an upgrade.
A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Below you find a list of all the changes with respect to version 2.1.
Notice that version 2.2 includes already all the fixes in branch 2.1 up to 2.1.4 indicated in \ref CHANGES-2-1 .

Changes from version 2.1 which are relevant for users:
- Changes leading to incompatible behavior:
  - Labels of quantities calculates by \ref SPRINT have changed from <i>label</i>.coord_<i>num</i> to <i>label</i>.coord-<i>num</i>
  - \ref METAD with WALKERS_MPI now writes a single hills file, without suffixes
  - removed the ./configure.sh script of v2.0.x, now plumed can only be configured using autotools (./configure)
  - \ref COM, \ref CENTER, and \ref GYRATION now automatically make molecules whole. In case you do not want them to do it, use NOPBC flag,
    which recovers plumed 2.1 behavior
  - Some MD code could now automatically trigger restart (e.g. gromacs when starting from cpt files). This can be overwritten using
    \ref RESTART NO.
  - Replica suffixes are now added by PLUMED *before* extension (e.g. use plumed.0.dat instead of plumed.dat.0)
  - When using \ref switchingfunction the STRETCH keyword is now implicit. NOSTRETCH is available to enforce the old behavior.
- Module activation can now be controlled during configure with `--enable-modules` option.
- Almost complete refactoring of installation procedure. Now
  DESTDIR and other standard autoconf directories (e.g. bindir) are completely supported.
  Additionally, everything should work properly also when directory names include spaces (\issue{157}).
  Finally, compiler is not invoked on install unless path are explicitly changed (\issue{107}).
- Related to installation refactoring, upon install a previously installed PLUMED is not removed.
  This is to avoid data loss if prefix variable is not properly set
- Several changes have been made in the Makefile.conf that makes it not compatible with those
  packaged with plumed 2.0/2.1. Please use ./configure to generate a new configuration file.
- Added partial OpenMP parallelization, see \ref Openmp
- Added multiple time step integration for bias potentials, see \ref MTS
- Link cells are now used in all multicolvars that involve \ref switchingfunction.  The link cell cutoff is
  set equal to 2.*\f$d_{\textrm{max}}\f$.  Where \f$d_{\textrm{max}}\f$ is the (user-specified) point at which
  the switching function goes to zero. Users should always set this parameter when using a switching function
  in order to achieve optimal performance.
- DHENERGY option is no longer possible within \ref DISTANCES.  You can still calculate the DHENERGY colvar by using \ref DHENERGY
- Reweighting in the manner described in \cite Tiwary_jp504920s is now possible using a combination of the \ref METAD and \ref HISTOGRAM actions.  The relevant keywords in \ref METAD are REWEIGHTING_NGRID and REWEIGHTING_NHILLS.  The \f$c(t)\f$ and the appropriate weight to apply to the configurations are given by the values labeled rct and rbias. 
- News in configure and install:
  - ./configure now allows external BLAS to be used with internal LAPACK. This is done automatically if only BLAS are available,
    and can be enforced with --disable-external-lapack.
  - ./configure supports --program-prefix, --program-suffix, and --program-transform-name.
  - make install supports DESTDIR and prefix.
  - Environment variables PLUMED_LIBSUFFIX and PLUMED_PREFIX are deprecated and will be removed in a later version.
- New actions
  - \ref DUMPMASSCHARGE to dump a file with mass and charges during MD.
  - \ref EFFECTIVE_ENERGY_DRIFT to check that plumed forces are not screwing the MD integrator.
  - \ref EXTENDED_LAGRANGIAN : in combination with  \ref METAD it implements metadynamics with Extended Lagrangian; standalone it implements TAMD/dAFED.
  - \ref DFSCLUSTERING calculate the size of clusters 
  - \ref DUMPMULTICOLVAR print out a multicolvar
  - \ref MFILTER_LESS filter multicolvar by the value of the colvar
  - \ref MFILTER_MORE 
  - \ref MFILTER_BETWEEN
  - \ref PCARMSD PCA collective variables using OPTIMAL rmsd measure
  - \ref PCAVARS PCA collective variables using any one of the measures in reference
  - \ref GRADIENT can be used to calculate the gradient of a quantity.  Used to drive nucleation
  - \ref CAVITY
  - \ref PUCKERING implemented for 5-membered rings (thanks to Alejandro Gil-Ley).
  - \ref WRAPAROUND to fix periodic boundary conditions.
- New features for existing actions:
  - Keywords UPDATE_FROM and UPDATE_UNTIL to limit update step in a defined time window, available only for actions where it would be useful.
  - Keyword UNNORMALIZED for \ref HISTOGRAM.
  - Possibility to use Tiwary-Parrinello reweighting for \ref METAD
  - Keywords for \ref GROUP (REMOVE, SORT, UNIQUE) to allow more flexible editing of groups.
  - \ref DUMPATOMS now supports dumping xtc and trr files (requires xdrfile library).
  - \ref driver can now read xtc and trr files also with xdrfile library.
  - \ref driver accepts a --mc flag to read charges and masses from a file produced during
    molecular dynamics with \ref DUMPMASSCHARGE
  - Possibility to enable or disable \ref RESTART on a per action basis, available only for actions where it would be useful.
  - \ref MOLINFO now supports many more special names for rna and dna (thanks to Alejandro Gil-Ley).
  - VMEAN and VSUM allow one to calculate the sum of a set of vectors calculated by VectorMultiColvar.  Note these
  can also be used in tandem with \ref AROUND or \ref MFILTER_MORE to calculate the average vector within a particular
  part of the cell or the average vector among those that have a magnitude greater than some tolerance
  - New way of calculating the minimum value in multicolvars (ALT_MIN). This is less susceptible to overflow for certain 
    values of \f$\beta\f$.  
  - New keywords for calculating the LOWEST and HIGHEST colvar calculated by a multicolvar
  - Added components to \ref DIPOLE (\issue{160}).
- Other changes:
  - File reader now supports dos newlines as well as files with no endline at the end.

For developers:

- In order to be able to use openMP parallelism within multicolvar, secondarystructure, manyrestraints and crystallisation
we had to make some substantial changes to the code that underlies these routines that is contained within vesselbase. In 
particular we needed to get rid of the derivatives and buffer private variables in the class ActionWithVessel.  As a consequence
the derivatives calculated in the various performTask methods are stored in an object of type MultiValue.  Within multicolvar
this is contained within an object of type AtomValuePack, which stores information on the atom indices.  If you have implemented
a new multicolvar it should be relatively straightforward to translate them so they can exploit this new version of the code.  Look 
at what has been done to the other multicolvars in there for guidance.  Sorry for any inconvenience caused.
- Changed the logic of several PLUMED ifdef macros so as to make them consistent.
  Now every feature based on external libraries is identified by a __PLUMED_HAS_* macro.

Version 2.2.1 (Jan 18, 2016)
---------------------------------------------

For users:
- \ref PBMETAD implement the new Parallel Bias Metadynamics flavor of the Metadynamics sampling method.
- PLUMED now reports an error when using \ref HISTOGRAM with UNNORMALIZED without USE_ALL_DATA. See \issue{175}
- Fixed a bug in configure together with --enable-almost. The check for lbz2 library was not working properly.
- Fixed a bug in install procedure that was introducing an error in linking with CP2K.
- Fixed a bug that sometimes was preventing the printing of a useful error message.

For developers:
- Vector and Tensor now support direct output with `<<`.
- Added some missing matmul operation Vector and Tensor.
- ./configure is automatically relaunched when changing ./configure or Makefile.conf. This makes it more robust
  to switch between branches.

Version 2.2.2 (Apr 13, 2016)
----------------------------------------------

For users:
- \ref MOLINFO for RNA accepts more residue names, see \issue{180}.
- added two mpi barries (one was missing in PBMetaD for multiple walkers) to help synchronized initialisation
- Fixed a bug in internal stopwatches that was making \ref DEBUG logRequestedAtoms not working
- Some multicolvars (including \ref BRIDGE, \ref ANGLES, and \ref INPLANEDISTANCES) now crashes if one
  asks for too many atoms, see \issue{185}.
- Optimisations (activation of the dependencies, secondary structures, DRMSD)
- Fixed a performance regression with RMSD=OPTIMAL-FAST
- Fixed a bug in the normalization of kernel functions (relevant for \ref HISTOGRAM).
- Fixed a regression introduced in v2.2 that was making \ref METAD with non-MPI multiple walkers crash
  if reading frequently. See \issue{190}
- Updated patch for gromacs 5.x. Patches for gromacs 5.0 and 5.1 have been fixed so as to allow
  patching in runtime mode.
- Possibility to control manual generation (including pdf) from ./configure. Pdf manual is now off
  by default. Notice that on travis CI it is still generated.

For developers:
- Fixed a bug in the interpretation of cmd strings. Namely, an erroneous string was not triggering an error.
  This is harmless for MD codes properly patched, but could have introduced problems in MD codes with typoes
  in cmd strings.
- ./configure is not automatically relaunched anymore when doing `make clean`.

Version 2.2.3 (Jun 30, 2016)
----------------------------------------------

For users:
- Updated patches for gromacs 5.1.x and 5.0.x to fix a problem when plumed was trying to write to an already
  closed gromacs log file.
- When looking for a value outside the GRID now the error include the name of the responsible 
  collective variable
- Numerical check in LatticeReduction made less picky. This should solve some of the internal errors reported
  by `LatticeReduction.cpp` when using aggressive compilers.
- Files are now flushed at the correct step. Before this fix, they were flushed at the step before the requested one
  (e.g. with \ref FLUSH STRIDE=100 at step 99, 199, etc).
- In \ref METAD, INTERVAL with periodic variables now report an error.
- \ref LOAD now works also when plumed is installed with a suffix.
- Added `--md-root` option to `plumed patch` which allows it to be run from a directory different from the one
  where the md code is located.
- Wham script in \ref munster tutorial now writes weights in scientific notation.

For developers:
- `./configure` checks if dependencies can be generated. If not, they are disabled.
- Added --disable-dependency-tracking to ./configure
- Added a make target `all_plus_doc` that builds both code and docs.
- Added possibility to set a default location for plumed library in runtime binding.
  If the plumed wrapped is compiled with `-D__PLUMED_DEFAULT_KERNEL=/path/libplumedKernel.so`,
  then if the env var PLUMED_KERNEL is undefined or empty PLUMED will look in the path at compile time.
- Tentative port files are now available at [this link](http://github.com/plumed/ports). 
  They can be used to install PLUMED using MacPorts.

Version 2.2.4 (Dec 12, 2016)
-------------

For users:
- Fix a bug in \ref PBMETAD when biasing periodic and not periodic collective variables at the same time 
- GSL library is now treated by `./configure` in the same way as other libraries, that is `-lgsl -lgslcblas` are only
  added if necessary.
- Fix a bug in \ref METAD when using INTERVAL and ADAPTIVE gaussians at the same time
- Updated gromacs patch for 5.1.x to 5.1.4
- Fix a performance regression in the calculate loop where derivatives and forces were set to zero even if an action
  was not active, this is relevant for postprocessing and for the on-the-fly analysis
- Torsion calculation has been made slightly faster and improved so as to provide correct
  derivatives even for special angles (e.g. +pi/2 and -pi/2).

For developers:
- Macports portile is now tested on travis at every plumed push.

Version 2.2.5 (Mar 31, 2017)
-------------

\plumednotmaintained

For users:
- Fixed a problem with large step numbers in driver (see \issue{209}).
- Fixed a problem leading to crashes when using switching functions without cutoff with some compiler (see \issue{210}).
- Fixed a bug when using \ref FIT_TO_TEMPLATE and domain decomposition (see \issue{214}).
- Added an automatic flush of HILLS files when using \ref METAD with file-based multiple walkers.
- Root dir is logged to allow easier debugging of problems.

@page CHANGES-2-3 Version 2.3

## Version 2.3 (Dec 12, 2016)

Version 2.3 contains several improvements with respect to 2.2. Users currently working with 2.2
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files.

Below you find a list of all the changes with respect to version 2.2.
Notice that version 2.3 includes already all the fixes in branch 2.2 up to 2.2.3 indicated in \ref CHANGES-2-2 .

Changes from version 2.2 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref COMMITTOR can now be used to define multiple basins, but the syntax has been changed
  - Syntax for \ref SPRINT and \ref DFSCLUSTERING has changed.
    We have separated the Actions that calculate the contact matrix from these actions.  These actions thus now take a contact
    matrix as input.  This means that we these actions can be used with contact matrices that measures whether or not a pair of atoms
    are hydrogen bonded.  For more details on this see \ref contactmatrix.  For clustering the output can now be passed to the actions
    \ref CLUSTER_PROPERTIES, \ref CLUSTER_DIAMETER, \ref CLUSTER_NATOMS, \ref OUTPUT_CLUSTER and \ref CLUSTER_DISTRIBUTION.  These
    provide various different kinds of information about the connected components found by clustering 
  - In \ref driver masses and charges are set by default to NaN.
    This makes it less likely to do mistakes trying to compute centers of mass or electrostatic-dependent variables
    when masses or charges were not set. To compute these variables from the driver you are now forced to use
    `--pdb` or `--mc`.
  - In rational switching functions, by default MM is twice NN. This is valid both in \ref switchingfunction with expanded
    syntax and when specifying MM on e.g. \ref COORDINATION
  - Patch script `plumed patch` now patches by default with `--shared`. This should make the procedure more robust (see \issue{186}).
  - Faster \ref GYRATION but new default behavior is not mass weighted
  - When using \ref HISTOGRAM you now output the accumulated grid using \ref DUMPGRID or \ref DUMPCUBE to get the free energy you use
    the method \ref CONVERT_TO_FES.  These changes allow one to use grids calculated within PLUMED in a work flow of tasks similarly to 
    the way that you can currently use Values.
  - The way that reweighting is performed is now different.  There are three separate actions \ref REWEIGHT_BIAS, \ref REWEIGHT_TEMP and
    \ref REWEIGHT_METAD.  These actions calculate the quantities that were calculated using the keywords REWEIGHT_BIAS and REWEIGHT_TEMP that
    used to appear in the old HISTOGRAM method.  Now those these methods can be used in any methods that calculate ensemble averages for
    example \ref HISTOGRAM and \ref AVERAGE
  - Manual is now build with locally compiled plumed
  - Removed CH3SHIFT
  - \ref CS2BACKBONE is now native in PLUMED removing the need to link ALMOST, small syntax differences
  - \ref CS2BACKBONE, \ref NOE, \ref RDC, removed the keyword ENSEMBLE: now ensemble averages can only be calculated using \ref ENSEMBLE
  - \ref RDC, syntax changes
  - It is not possible anymore to select modules using `modulename.on` and `modulename.off` files. Use `./configure --enable-modules` instead.
  - Removed IMD modules. In case someone is interested in restoring it, please contact the PLUMED developers.
- New actions:
  - \ref FIXEDATOM
  - \ref HBOND_MATRIX
  - \ref CLUSTER_PROPERTIES
  - \ref CLUSTER_DIAMETER
  - \ref CLUSTER_NATOMS
  - \ref OUTPUT_CLUSTER
  - \ref CLUSTER_DISTRIBUTION
  - \ref ROWSUMS
  - \ref COLUMNSUMS
  - \ref UPDATE_IF
  - \ref DUMPGRID
  - \ref DUMPCUBE
  - \ref CONVERT_TO_FES
  - \ref INTERPOLATE_GRID
  - \ref FIND_CONTOUR
  - \ref FIND_SPHERICAL_CONTOUR
  - \ref FIND_CONTOUR_SURFACE
  - \ref AVERAGE
  - \ref REWEIGHT_BIAS
  - \ref REWEIGHT_TEMP
  - \ref REWEIGHT_METAD
  - \ref PCA
  - \ref PRE
  - \ref STATS
  - \ref METAINFERENCE
  - \ref LOCALENSEMBLE
  - \ref FRET
  - \ref RESET_CELL
  - \ref JCOUPLING
  - \ref ERMSD
- New features in MD patches (require re-patch):
  - Patch for amber 14 now passes charges with appropriate units (fixes \issue{165}). Notice that
    the patch is still backward compatible with older PLUMED version, but the charges will only be passed
    when using PLUMED 2.3 or later.
  - Patch for GROMACS 5.1 incorporates Hamiltonian replica exchange, see \ref hrex
  - Gromacs 2016, 5.1.x, 5.0.x, flush the plumed output files upon checkpointing
  - Added patch for Gromacs 2016.1
  - gromacs 5.1.x patch updated to 5.1.4
  - Removed the patch for Gromacs 4.6.x 
  - LAMMPS patch updated to support multiple walkers and report plumed bias to LAMMPS (thanks to Pablo Piaggi).
- New features for existing actions:
  - The SPECIES and SPECIESA keyword in MultiColvars can now take a multicolvar as input.  This allows one
    to calculate quantities such as the Q4 parameters for those atoms that have a coordination number greater
    than x.
  - Added MATHEVAL type in \ref switchingfunction
  - Added Q type native contacts in \ref switchingfunction (thanks to Jan Domanski).
  - \ref COMMITTOR can now be used to define multiple basins
  - The number of atoms admitted in \ref BRIDGE has been significantly increased, see \issue{185}.
  - \ref driver now allows --trajectory-stride to be set to zero when reading with --ixtc/--itrr. In this case, step number is read from the trajectory file.
  - \ref METAD and \ref PBMETAD can now be restarted from a GRID 
  - Added keywords TARGET and DAMPFACTOR in \ref METAD
  - When using \ref METAD with file-based multiple walkers and parallel jobs (i.e. mpirun) extra suffix is not added (thanks to Marco De La Pierre).
  - \ref ENSEMBLE added keywords for weighted averages, and calculation of higher momenta
  - \ref MOLINFO now allows single atoms to be picked by name.
  - \ref FIT_TO_TEMPLATE now supports optimal alignment.
  - \ref CONSTANT added the possibility of storing more values as components with or without derivatives
  - \ref PUCKERING now supports 6 membered rings.
  - Extended checkpoint infrastructure, now \ref METAD and \ref PBMETAD will write GRIDS also on checkpoint step (only the GROMACS patch
    is currently using the checkpointing interface)
- Other features:
  - Added a plumed-config command line tool. Can be used to inspect configuration also when cross compiling.
  - Added a `--mpi` option to `plumed`, symmetric to `--no-mpi`. Currently, it has no effect (MPI is initialized by default when available).
  - PLUMED now generate a VIM syntax file, see \ref VimSyntax
  - The backward cycle is now parallelized in MPI/OpenMP in case many collective variables are used.
  - GSL library is now searched by default during `./configure`.
  - Tutorials have been (partially) updated to reflect some of the changes in the syntax
  - Parser now reports errors when passing numbers that cannot be parsed instead of silently replacing their default value. See \issue{104}.
  - More and more documentation
- Bug fixes:
- Fixed a bug in \ref PBMETAD that was preventing the writing of GRIDS if a hill was not added in that same step 

For developers:
- IMPORTANT: BIAS can now be BIASED as well, this changes can lead to some incompatibility: now the "bias" component is always defined automatically
  by the constructor of Bias as a componentWithDerivatives, derivatives are automatically obtained by forces. The main change is that you don't have to define
  the bias component anymore in your constructor and that you can use setBias(value) to set the value of the bias component in calculate. 
- Added new strings for plumed cmd: setMDMassUnits, setMDChargeUnits, readInputLine, performCalcNoUpdate, update and doCheckPoint.
- Easier to add actions with multiple arguments
- New functions to access local quantities in domain decomposition
- Active modules to enable regtests are chosen using `plumed config`.
- A script is available to check if source code complies plumed standard. Notice that this script is run together with cppcheck on travis-ci.
- Cppcheck on travis-ci has been updated to 1.75. Several small issues triggering errors on 1.75 were fixed (e.g. structures passed by value
   are now passed by const ref) and false positives marked as such.
- Added coverage scan.

## Version 2.3.1 (Mar 31, 2017)

- Fix to FIT_TO_TEMPLATE as in 2.2.5. Notice that in 2.3.0 also the case with TYPE=OPTIMAL was affected. This is fixed now.
- small change in \ref CS2BACKBONE to symmetrize the ring current contribution with respect to ring rotations (also faster)
- fixed `plumed-config` that was not working.
- log file points to the `config.txt` files to allow users to check which features were available in that compiled version.
- `make clean` in root dir now also cleans `vim` sub-directory.
- Updated gromacs patch to version 2016.3 

For developers:
- Cppcheck on travis-ci has been updated to 1.77.
- Doxygen on travis-ci has been updated to 1.8.13

## Version 2.3.2 (Jun 12, 2017)

See branch \branch{v2.3} on git repository.

- Resolved problem with nan in \ref SMAC with SPECIESA and SPECIESB involving molecules that are the same
- PDB reader is now able to read files with dos newlines (see \issue{223}).
- Fixed bug in \ref CS2BACKBONE (v2.3.1) related to ring currents of HIS and TRP
- Fixed bug in if condition in \ref PCAVARS so that you can run with only one eigenvector defined in input 
- Fixed bug with timers in \ref sum_hills \issue{194}.
- Fixed bug when using \ref MOVINGRESTRAINT with periodic variables such as \ref TORSION \issue{225}.
- Fixed bug in \ref HBOND_MATRIX that used to appear when you used DONORS and ACCEPTORS with same numbers of atoms 
- Fixed bug in \ref DISTANCES that appears when using BETWEEN and link cells.
- Prevented users from causing segfaults by storing derivatives without LOWMEM flag.  In these cases PLUMED crashes with meaningful errors.
- Fixed bug in \ref HISTOGRAM that causes NaNs when using KERNEL=DISCRETE option
- Fixed a bug in the parser related to braces, see \issue{229}
- Fixed a bug that appeared when using \ref Q3, \ref Q4 and \ref Q6 with LOWEST or HIGHEST flag
- Fixed a bug that appears when you use \ref MFILTER_LESS as input to \ref COORDINATIONNUMBER with SPECIESA and SPECIESB flags
- Fixed a bug that was making flushing when gromacs checkpoints not functional (thanks to Summer Snow).
- Fixed a bug affecting \ref EXTENDED_LAGRANGIAN and \ref METAD with ADAPT=DIFF when using an argument
  with periodicity (min,max) such that min is different from -max.
  This does not affect normal \ref TORSION, but would affect \ref PUCKERING component phi
  with 6-membered rings. In addition, it would affect any variable that is created by the user with a periodicity
  domain not symmetric around zero. See \issue{235} (thanks to Summer Snow for reporting this bug).
- Fixed numerical issue leading to simulations stuck (LatticeReduction problem) with intel compiler and
  large simulation cells.
- Fixed a bug affecting \ref LOCAL_AVERAGE and outputting all multicolvars calculated by \ref Q6 with \ref DUMPMULTICOLVAR
- `plumed info --user-doc` and `plumed info --developer-doc` now fall back to online manual when local doc is not installed,
  see \issue{240}.

For developers:
- IMPORTANT: we started to enforce code formatting using astyle. Check the developer documentation to learn how to
  take care of not-yet-formatted branches.
- plumedcheck validation has been made stricter. All the checks are now described in the developer manual.
- New flag `--disable-libsearch` for `configure`, allowing an easier control of linked libraries when installing PLUMED
  with a package manager such as MacPorts.
- Added `--disable-static-patch` to `./configure` to disable tests related to static patching. It can be used
  when static patching is not needed to make sure a wrong c++ library is not linked by mistake.
- Using `install_name_tool` to fix the name of the installed library on OSX. Allows linking the PLUMED
  shared library without explicitly setting `DYLD_LIBRARY_PATH`.
- Added environment variable `PLUMED_ASYNC_SHARE` to enforce synchronous/asynchronous atom sharing (mostly for debug purpose).
- On travis-ci, using ccache to speedup builds.
- On travis-ci, added a regtest using Docker with gcc6 and MPI.
- On travis-ci, docs for unofficial or unsupported branches are set not to be indexed by search engines (see \issue{239})
- Cppcheck on travis-ci has been updated to 1.79.

## Version 2.3.3 (Oct 3, 2017)

For users:
- Fixed a bug in \ref switchingfunction MATHEVAL, leading to inconsistent results when using OpenMP with multiple threads (see \issue{249}).
- \ref FIT_TO_TEMPLATE now reports when it is used with a reference file with zero weights.
- Fixed logging of \ref UNITS (thanks to Omar Valsson).
- Fixed a possible bug with \ref EFFECTIVE_ENERGY_DRIFT and domain decomposition with a domain containing zero atoms.


For developers:
- Fixed a bug in `./configure --disable-libsearch` when searching for molfile plugins.
- Cppcheck on travis-ci has been updated to 1.80.
- Configure script now has a list of better alternatives to find a working `ld -r -o` tool to merge object files.
  This solves linking issues on some peculiar systems (see \issue{291}, thanks to Massimiliano Culpo). 
- Using `install_name_tool` also on non-installed libraries. This makes it possible to link them and later
  find them without explicitly setting `DYLD_LIBRARY_PATH`. This should also make the `DYLD_LIBRARY_PATH` irrelevant.
  Notice that `DYLD_LIBRARY_PATH` is not well behaved in OSX El Capitan.

## Version 2.3.4 (Dec 15, 2017)

For users:
- GROMACS patch updated to gromacs-2016.4. This patch was also fixed in order to properly work with \ref ENERGY (see \issue{316})
  and to implement `-hrex` option (see \issue{197}).
- Patch for GROMACS 5.1.4 updated to fix an error with \ref ENERGY (see \issue{316}).
- Solved a bug in \ref ERMSD leading to incorrect results when using non-default length units (e.g. with `UNITS LENGTH=A`).

For developers:
- Regtest script also reports when exitcode different from zero is returned.
- Patch script reports errors returning a nonzero exit code.
- cppcheck update to 1.81
- Solved small bug in stored PLUMED_ROOT directory as obtained from statically patched MD codes.
  Namely, the compilation directory was stored rather than the installation one.

## Version 2.3.5 (Mar 2, 2018)

For users:
- Fixed `plumed partial_tempering` to agree with GROMACS conventions for the choice of dihedral angles (see \issue{337}).
  Should be irrelevant for the vast majority of cases.
- Fixed small bug in regexp parser - the part outside the parentheses was just ignored.

For developers:
- Doxygen on travis-ci has been updated to 1.8.14.
- Embedded astyle updated to 3.1.
- `make clean` now correctly removes the `src/lib/plumed` executable.

## Version 2.3.6 (Jul 2, 2018)

For users:
- Fixed a problem leading to NaN derivatives of \ref switchingfunction `Q` when distance between two atoms is large.
- GROMACS patch updated to gromacs-2016.5.
- `./configure` crashes if prefix is set to present working directory (notice that this choice was already leading to issues).
- \ref DUMPATOMS reports an error when trying to write xtc/xdr files without the xdrfile library installed.
- Fixed a bug appearing when using \ref PATH or \ref GPROPERTYMAP with virtual atoms without simultaneously using the same
  atoms in a different action.
- Fixed incorrect format of the pdb file written by \ref PCA (see \issue{363}).
- Fixed behavior of natural units. When an MD code asks for natural units, it is not necessary to also set units within PLUMED using \ref UNITS (see \issue{364}).

For developers:
- Fixed small issue in debug options of \ref driver (see \issue{245}).
- `plumed patch -e` now accepts a name closely matching the patch name (e.g. `plumed patch -e gromacs2016.5` will try to patch
  even if the stored patch is for `gromacs-2016.4`). This simplifies managing Portfiles. Nothing changes when picking the patch
  from the interactive menu.
- Install newer ccache on travis-ci, build faster.
- Small fix in provided env modules (`PLUMED_VIMPATH` is set also when shared libraries are disabled).

## Version 2.3.7 (Oct 5, 2018)

For users:
- Fixed flag DETAILED_TIMERS in \ref DEBUG (flag was ignored and detailed timers always written).
- Small fix in \ref DUMPMASSCHARGE (atoms are now correctly requested only at first step).

## Version 2.3.8 (Dec 19, 2018)

\plumednotmaintained

For users:
- Fixed some openMP regression (some related to the whole codes and some specifics for Coordination and Multicolvar), this were compiler dependent so not all users may have experienced them
- Fixed an issue with \ref CS2BACKBONE when more than 2 chains were used
- Fixed memory leak in \ref RDC.
- Fixed segmentation fault with more than two CVs in reweighting \ref METAD (see \issue{399}, thanks to Fiskissimo).

For developers:
- Small fix in LDFLAGS when enabling coverage.
- Fixed order of flags in tests for static linking done by configure (see \issue{407}).
- Fixed the way paths are hard-coded so as to facilitate conda packaging (see \issue{416}).


*/
@page CHANGES-UNRELEASED Unreleased changes 

This page contains changes that will end up in 2.7

Changes from version 2.6 which are relevant for users:

- New contributed modules:
  - A new Funnel module by Stefano Raniolo and Vittorio Limongelli 
     - \ref FUNNEL_PS 
     - \ref FUNNEL 


For developers:
- small fix in `Plumed.h` too avoid unique global symbols (see \issue{549})


@page CHANGES-2-0 Version 2.0

Version 2.0.0 (September 27, 2013)
----------------------------

Version 2.0 is a complete rewrite, so there is no way to write a complete set of difference
with respect to plumed 1.3. Here is a possibly incomplete summary of the difference:
- The input is simpler, more flexible, and more error proof.
  Many checks are now performed and in this way common errors are avoided. 
- The units are now the same for all MD codes.
  If you want to use a different unit than the default you set it in the input file. 
- The analysis tools are now much more flexible.
  As an example of this it is now possible to write different collective variables with different frequencies.
- Many complex collective variables are considerably faster than they were in plumed1.
  In particular, all variables based on RMSD distances. 
- Centers of mass can be used as if they were atoms.
  Hence, unlike plumed 1.3, you can use center of mass positions in ALL collective variables.
- The virial contribution is now computed and passed to the MD code.
  Plumed can thus now be used to perform biased NPT simulations.
- Variables can be dumped on different files, and are
  computed only when this is necessary.
- PLUMED is now compiled as a separate library. This simplifies the patching
  procedure, but might require some extra work to configure PLUMED properly.
  Since PLUMED can be loaded as a shared library, it is possible to setup
  everything such that PLUMED and MD codes can be updated independently from each
  other.

In addition, it is now much easier to contribute new functionality to the code because: 
- There is a much simpler interface between plumed and the base MD codes.
  This makes it much easier to add plumed to a new MD code. Hopefully, in the future,
  interfaces with MD codes will be maintained by the developers of the MD codes
  independently from PLUMED developers. This will allow more MD codes
  to be compatible with PLUMED.
- There is C++ object oriented programming and full compatibility with the C++ standard library 
- A modular structure.
- New collective variables and methods can be released independently.
- There is an extensive developer documentation.
- User documentation is provided together inside the implementation files.

Caveats:
- PLUMED 2 input file (plumed.dat) has a syntax which is not
  compatible with PLUMED 1.
  Transition should be easy, but cannot
  be done just using the new version with the old input file.
- PLUMED 2 is written in C++, thus requires a C++ compiler
- PLUMED 2 may not include all the features that were available
  in PLUMED 1.

A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Version 2.0.1 (Nov 14, 2013)
----------------------------

For users:
- Fixed a bug in \ref HISTOGRAM with REWEIGHT_BIAS. Reweighting was only done when also temperature-reweighting was enabled.
- Fixed a bug that was sometime crashing code with domain decomposition and
  non-dense simulation boxes (e.g. implicit solvent).
- Performance improvements for \ref GYRATION.
- Flush all files every 10000 steps by default, without need to use \ref FLUSH
- Errors when writing input for \ref switchingfunction are now properly
  recognized.
- Added message when \ref simplemd is used on a non-existing file.
- Fixed `plumed mklib` such that it deletes the target shared library in case
  of compilation error.
- Several small fixes in documentation and log file.

For developers:
- Added possibility to setup replica exchange from MD codes in Fortran (commands "GREX setMPIFIntercomm" and "GREX setMPIFIntracomm").
- cmd("setStopFlag") should now be called after PLUMED initialization.
- Several small fixes in documentation.

Version 2.0.2 (Feb 11, 2014)
----------------------------

For users:
- Fixed bug with \ref METAD with INTERVAL and replica exchange, including bias exchange.
  Now the bias is correctly computed outside the boundaries. Notice that this is different
  from what was done in PLUMED 1.3. Also notice that INTERVAL now works
  correctly with grids and splines.
- Fixed bug with \ref READ and periodic variables.
- Fixed bug with \ref HISTOGRAM (option USE_ALL_DATA was not working properly).
- Gromacs patch updated to 4.6.5.
- Gromacs patch for 4.6 has been modified to allow for better load balancing when
  using GPUs.
- Added option 'plumed info --long-version' and 'plumed info --git-version'.
- Added full reference (page/number) to published paper in doc and log.
- Fixed a bug in file backups (only affecting Windows version - thanks to T. Giorgino).
- Added possibility to search in the documentation.
- Several small fixes in documentation and log file.

For developers:
- Fixed Makefile dependencies in some auxiliary files in src/lib (*cmake and *inc).
- Changed way modules are linked in src/.
  E.g. src/colvar/tools/ is not anymore a symlink to src/colvar but a real directory.
  (Notice that this introduces a regression: when using plumed as an external library
  some include files could not work - this only applies when plumed is installed;
  also notice that this is fixed in 2.0.3)
- Patch for gromacs 4.6 now also include original code so as to simplify its modification.
- Added option 'plumed patch --save-originals'.
- Fixed regtest regtest/secondarystructure/rt32 to avoid problems with NUMERICAL_DERIVATIVES.
- Removed include graphs in the documentation (too large).
- Several small fixes in documentation.

Version 2.0.3 (June 30, 2014)
----------------------------

For users:
- Now compiles on Blue Gene Q with IBM compilers.
- Fixed bug in \ref CENTER where default WEIGHTS were missing. 
- Fixed broken \ref CONTACTMAP with SUM
- Fixed \ref DUMPATOMS with gro file and more than 100k atoms.
- Added CMDIST in \ref CONTACTMAP to emulate plumed1 CMAP.
- Several small fixes in documentation and log file.

For developers:
- Fixed cmd("getBias") to retrieve bias. It was not working with
  single precision codes and it was not converting units properly.
- Fixed a regression in 2.0.2 concerning include files from installed plumed
  (see commit 562d5ea9dfc3).
- Small fix in tools/Random.cpp that allows Random objects to be
  declared as static.
- Small fix in user-doc compilation, so that if plumed is not found
  the sourceme.sh file is sourced
- Fixed non-ANSI syntax in a few points and a non-important memory leakage.
- Split cltools/Driver.cpp to make parallel compilation faster.

Version 2.0.4 (September 15, 2014)
----------------------------------------------

For users:
- Fixed a bug in \ref BIASVALUE that could produce wrong acceptance with replica exchange simulations.
- Fixed a few innocuous memory leaks.
- Fixed reader for xyz files, that now correctly detects missing columns. Also a related regtest has
  been changed.
- Several small fixes in documentation and log file.

For developers:
- Renamed Value.cpp to BiasValue.cpp

Version 2.0.5 (December 15, 2014)
----------------------------------------------

\plumednotmaintained

For users:
- Fixed a bug in replica exchange with different Hamiltonians (either lambda-dynamics
  or plumed XX-hrex branch) possibly occurring when using charge or mass dependent
  variables.
- Fixed a bug in analysis (e.g. \ref HISTOGRAM) leading to wrong accumulation
  of statistics when running a replica exchange simulation.
- Fixed a bug in the calculation of derivatives in histograms. This should
  be harmless since people usually only consider the value in histograms
  and not the derivatives.
- Fixed an issue in Makefile that could results in problems when
  patching an MD code with --shared option (pointed out by Abhi Acharya).
  This fixes a regression introduced in 2.0.2.
- Small fixes in documentation.

For developers:
- Added warning when performing regtests using an instance of plumed from
  a different directory

@page CHANGES-2-4 Version 2.4

## Version 2.4 (Dec 15, 2017)

Version 2.4 contains several improvements with respect to 2.3. Users currently working with 2.3
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files.
Notice that version 2.4 includes already all the fixes in branch 2.3 up to 2.3.3 indicated in \ref CHANGES-2-3 .

Changes from version 2.3 which are relevant for users:
- Changes leading to incompatible behavior:
  - A c++11 compliant compiler is required (see \issue{212}). This should mean:
    - gcc 4.8
    - clang 3.3
    - intel 15
    Since the number of c++11 features that we use is limited, older compilers might work as well.
  - The meaning of `BIASFACTOR=1` in \ref METAD has been modified and can now be used to indicate unbiased
    simulations. Non-well-tempered metadynamics is BIASFACTOR=-1, which is the new default value.
    Notice that this has an implication on the bias factor written in the HILLS file when doing
    non-well-tempered metadynamics.
  - Due to a change in \ref COMMITTOR, the format of its output file has been slightly changed.
  - \ref HISTOGRAM : When using weights default is now to output histogram divided by number of frames from which data was taken.  In addition the 
    UNORMALIZED flag has been replaced with the keyword `NORMALIZATION`, which can be set equal to true, false or ndata.
  - All switching functions are now stretched by default, also when using the "simple syntax" (e.g. `COORDINATION NN=6`).
    Switching functions were already stretched by default when using the advanced syntax (e.g. `COORDINATION SWITCH={}`)
    since version 2.2.  Notice that this will introduce small numerical differences in the computed switching functions.
- New modules:
  - A new PLUMED-ISDB module have been included, this module includes a number of CVs to calculate experimental data with the internal ability
    to also calculate a \ref METAINFERENCE score.
    - New actions include:
      - \ref EMMI
      - \ref SAXS
      - \ref RESCALE, \ref SELECT, \ref SELECTOR
    - Updated actions include:
      - \ref CS2BACKBONE
      - \ref FRET
      - \ref JCOUPLING
      - \ref METAINFERENCE
      - \ref NOE
      - \ref PRE
      - \ref RDC, \ref PCS
      - \ref PBMETAD
  - A new EDS module have been included, contributed by Glen Hocky and Andrew White.
    This module implements the following methods:
    - \ref EDS
  - A new DRR module have been included, contributed by Haochuan Chen and Haohao Fu.
    This module implements the following methods:
    - \ref DRR
    - \ref drr_tool
  - A new VES module have been included, contributed by Omar Valsson.
    This module implements the following methods:
    - \ref BF_CHEBYSHEV
    - \ref BF_COMBINED
    - \ref BF_COSINE
    - \ref BF_CUSTOM
    - \ref BF_FOURIER
    - \ref BF_LEGENDRE
    - \ref BF_POWERS
    - \ref BF_SINE
    - \ref OPT_AVERAGED_SGD
    - \ref OPT_DUMMY
    - \ref TD_CHI
    - \ref TD_CHISQUARED
    - \ref TD_CUSTOM
    - \ref TD_EXPONENTIAL
    - \ref TD_EXPONENTIALLY_MODIFIED_GAUSSIAN
    - \ref TD_GAUSSIAN
    - \ref TD_GENERALIZED_EXTREME_VALUE
    - \ref TD_GENERALIZED_NORMAL
    - \ref TD_GRID
    - \ref TD_LINEAR_COMBINATION
    - \ref TD_PRODUCT_COMBINATION
    - \ref TD_PRODUCT_DISTRIBUTION
    - \ref TD_UNIFORM
    - \ref TD_VONMISES
    - \ref TD_WELLTEMPERED
    - \ref VES_LINEAR_EXPANSION
    - \ref VES_OUTPUT_BASISFUNCTIONS
    - \ref VES_OUTPUT_FES
    - \ref VES_OUTPUT_TARGET_DISTRIBUTION
    - \ref ves_md_linearexpansion
- New collective variables:
  - \ref DIMER (thanks to Marco Nava).
  - \ref EEFSOLV : EEF1 implicit solvent solvation energy
  - \ref ADAPTIVE_PATH : Adaptive path variables using the method from \cite BerndAdaptivePath
- New actions:
  - \ref INENVELOPE
  - \ref TOPOLOGY_MATRIX
  - \ref BOND_DIRECTIONS
  - \ref DUMPGRAPH
  - \ref GRID_TO_XYZ
  - \ref INTEGRATE_GRID
  - \ref LWALLS
  - \ref MAXENT
  - \ref MCOLV_COMBINE
  - \ref MCOLV_PRODUCT
  - \ref POLYMER_ANGLES
  - \ref XANGLES , \ref YANGLES , \ref ZANGLES
  - \ref XYTORSIONS , \ref XZTORSIONS , \ref YXTORSIONS , \ref YZTORSIONS , \ref ZXTORSIONS , and \ref ZYTORSIONS
- New command line tools:
  - \ref pesmd : Tool for performing Langevin dynamics on an energy landscape that is specified using a PLUMED input file
  - \ref pathtools 
- Other changes:
  - Sharing coordinates and applying force is now faster (in some cases these can result in much better scaling of the performances in parallel).
  - \ref COMMITTOR : new flag to use committor to keep track of the visited basins without stopping the simulation
  - \ref PBMETAD : multiple walkers using files (thanks to Marco De La Pierre).
  - \ref PBMETAD : adaptive Gaussian kernels
  - \ref PBMETAD : default names for `GRID` and `FILE` (useful with many collective variables) 
  - \ref METAD : BIASFACTOR=1 is allowed and performs unbiased sampling. HILLS file can be used
    to recover free energy also in this case.
  - \ref METAD : a RECT option is available that allows setting an array of bias factors, one for each replica.
  - \ref METAD : added options to perform Transition Tempered Metadynamics (thanks to James Dama)
  - \ref PATHMSD and \ref PROPERTYMAP now support alignment to a close structure (thanks to Jana Pazurikova)
  - PDB files with more than 100k atoms can now be read using [hybrid 36](http://cci.lbl.gov/hybrid_36/) format,
    see \issue{226}.
  - Added lepton support. Set env var `export PLUMED_USE_LEPTON=yes` to activate lepton as a matheval replacement
    in \ref MATHEVAL, \ref CUSTOM, and \ref switchingfunction "MATHEVAL switching function".
    Notice that in v2.5 matheval support will be dropped and all these keywords will use lepton.
    See \issue{244}.
  - When parsing constants, PLUMED uses lepton library. This allows to pass
    arguments such as `HEIGHT=exp(0.5)` (see \ref parsing-constants).
  - \ref CUSTOM function has been added as an alias to \ref MATHEVAL .
  - Trajectories read in \ref driver also support the usual replica convention, that is if
    trajectory with replica suffix is not found the driver will look for a trajectory without the replica suffix.
  - A new syntax (`@replicas:`) can be used to specify different arguments for different replicas (see \ref special-replica-syntax).
  - Internal molfile implementation has been updated to VMD 1.9.3.
  - Examples in the documentation now have syntax highlighting and links to the documentation of used actions.
  - \ref COORDINATIONNUMBER : Added option to have pairwise distance moments of coordination number in the multicolvar module
  - GROMACS patch updated to gromacs-2016.4
  - Implemented HREX for gromacs-2016.4.
  - Added patch for Quantum ESPRESSO 6.2 (thanks to Ralf Meyer).
  - Fixed a bug in \ref LOCAL_AVERAGE which appears when you use `SPECIESA` and `SPECIESB` keywords instead of just `SPECIES`
  - Added possibility to pass `--kt` from \ref driver.

Changes from version 2.3 which are relevant for developers:
  - A few fixes has been made to improve exception safety. Although we still cannot declare
    PLUMED totally exception safe (there are still many non-safe pointers around),
    this made it possible to add a regtest that actually tests erroneous cmd strings
    and erroneous inputs.
  - Due to the required c++11 support, travis-ci test on Ubuntu Precise has been removed.
  - `gettimeofdate` and `gettime` have been replaced with portable `chrono` classes introduced in c++11.
  - C++ exceptions are enabled by default.
  - A large number of loops have been changed to use the `auto` keyword in order to improve code readability.
  - Stack trace is not written upon error anymore, unless environment variable `PLUMED_STACK_TRACE` is set at runtime.
  - Fixed a potential bug using single precision system BLAS on a mac (notice that currently plumed only uses
    double precision, so it is harmless).
  - Added `--enable-rpath` option for autoconf (off by default).
  - Files related to changelog are now stored as `.md` files. This makes
    it possible to navigate them from github.
  - `configure.ac` has been simplified and improved in order to more easily probe C++ libraries.
  - added `plumed_custom_skip` function to regtests in order to skip specific tests based on specific conditions (e.g. OS).
  - environment variable `LDSO` has been renamed to `LDSHARED`, which is standard in the python community.
  - a `libplumedWrapper.a` library is installed as well, that is used in `--runtime` patching.
  - pkgconfig files are installed.
  - `plumed config makefile_conf` can be used to retrieve `Makefile.conf` file a posteriori.
  - Store `MPIEXEC` variable at configure time and use it later for running regtests. Notice that in case
    `MPIEXEC` is not specified regtests will be run using the command stored in env var `PLUMED_MPIRUN` or, if this is
    also not defined, using `mpirun`.
  - Added canonical Makefile targets `check` and `installcheck`. Notice that `check` runs checks with
    non-installed plumed whereas `installcheck` uses the installed one, including its correct program name if it
    was personalized (e.g. with suffixes). Notice that this modifies the previously available `check` target.
    
    
## Version 2.4.1 (Mar 2, 2018)

For users:
  - Fixed an important bug affecting RMSD calculations with compilers supporting OpenMP 4 (e.g.: intel compiler). Notice that this bug might potentially affect not only
    \ref RMSD variable, but also \ref PATHMSD variables using RMSD, \ref FIT_TO_TEMPLATE, \ref PCAVARS, and possibly other variables based on RMSD calculations and optimal alignments
    (see \issue{343}). Results might depend on the exact architecture and on how aggressive is the compiler. The bug is a consequence of some erroneous SIMD directives introduced in 2.4.0, so it does not affect PLUMED 2.3.x. 
  - Resolved a problem with \ref CS2BACKBONE and glycine atom names.
  - Module VES: Fixed a bug with basis functions that have a constant function different from 1 (e.g. scaled version of the Legendre basis functions, \ref BF_LEGENDRE) that was causing a time-dependent shift in the bias potential.
  - Module VES: In optimizers (\ref OPT_AVERAGED_SGD and \ref OPT_DUMMY) the output of quantities related to the instantaneous gradients are now off by default as these quantities are generally not useful for normal users, their output can instead by re-enabled by using the `MONITOR_INSTANTANEOUS_GRADIENT` keyword. Also added an keyword `MONITOR_AVERAGE_GRADIENT` that allows to monitor the averaged gradient and output quantities related to it. 
  - \ref RMSD variable and other collective variables using reference PDB files now crash when zero weights are passed (see \issue{247}).
  - Using \ref COM with \ref driver without passing masses now triggers an error instead of reporting NaNs (see \issue{251}).

For developers:
  - `plumed patch -p` command can be used twice without triggering an error. This will allow e.g. building again
    on MacPorts in cases where the build was interrupted. Notice that this only works for patches without special
    after/before patch/revert functions.

## Version 2.4.2 (Jul 2, 2018)

For users:
  - All fixes done in version 2.3.6. Notice that \issue{363} in version 2.4 also applies to \ref pathtools.
  - Additional residue names (without the prefix `D`) are now supported by \ref MOLINFO for DNA. See \issue{367}.
  - Solved an important bug appearing in NAMD interface. Notice that the bug was a regression introduced in 2.4.0. As consequence, versions <= 2.3 and versions >=2.4.2
    are expected to work correctly. See \issue{254}.
  - GROMACS patch for gromacs-2018.1.
  - \ref VimSyntax now highlights `__FILL__` strings.
  - \ref METAD and \ref PBMETAD give a warning when one restarts a simulation and the old hills file is not found. See \issue{366}.

For developers:
  - `LDSHARED` is now correctly taken into account when launching `./configure`.
  - Fixed installation with `--disable-shared`.
  - Cppcheck upgraded to 1.84.

## Version 2.4.3 (Oct 5, 2018)

For users:
  - All fixes done in version 2.3.7.
  - Module VES: Fixed a bug in `TD_GRID` for 2D grids where the grid spacing is not the same for both dimensions.
  - GROMACS patch for gromacs-2018.3.

## Version 2.4.4 (Dec 19, 2018)

For users:
  - Fixed some performances regression issue with OpenMP
  - Updated NAMD patches to version 2.12 and 2.13. Old patches have been removed.
  - GROMACS patch for gromacs-2018.4.
  - Fixed a thread safety issue using forces on \ref HISTOGRAM 
  - Fixed error message suggesting wrong actions (see \issue{421}).

For developers:
  - All fixed done in version 2.3.8
  - Cppcheck updated to 1.85

## Version 2.4.5 (Apr 1, 2019)

For users:
  - Fixed an inconsistency in parsing of braces.
    It is now possible to pass individual options
    including spaces (e.g. with `FILE={/path with space/file}`). Notice 
    that this invalidates syntax such as `ATOMS={1}{2}{3}{4}`. See more
    at \issue{434}.
  - Fixed \ref simplemd so as to call "runFinalJobs" at the end of the simulation.
  - GROMACS patch for gromacs-2016.6.
  - GROMACS patch for gromacs-2018.6.
  - Added aliases for some actions/options containing dashes (`-`) in their name. This will improve
    backward compatibility when these actions/options will be removed (see \issue{449}).

## Version 2.4.6 (Jul 19, 2019)

For users:
  - Fixed a bug in \ref COORDINATIONNUMBER where derivatives were wrong when using `R_POWER` > 2, thanks to `@MoleOrbitalHybridAnalyst` for spotting and fixing
  - Fixed a bug in library search, possibly affecting linked blas/lapack on OSX (see \issue{476}).
  - Fixed a bug in \ref METAD with `TARGET` and `GRID_SPARSE` (see \issue{467}).

## Version 2.4.7 (Jan 27, 2020)

\plumednotmaintained

For users:
  - Fixed a bug with \ref CONVERT_TO_FES and periodic variables, see \issue{441} (backported from v2.5.3).
  - More robust backup for output files when running over multiple processes
  - Fixed a regression in the performances of `GEOMETRY` based flexible hills in \ref METAD and \ref PBMETAD
  - Fixed \issue{538}.
  - Fixed potential issue with VMD plugins from 1.9.4 (\issue{545}, thanks to Lixin Sun).
  - Module VES: Fixed an off-by-one bug in the output of target distribution averages. The bug only affects the output and does not affect results. The bug also affected the output of coefficients when using a bias cutoff. 
  - Module VES: Made sure that all relevant output files are written out at the final step when shutting down the simulation. This solves issues reported by @PabloPiaggi with restarting when there is a mismatch been the output of files and the number of MD steps. 


Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

WARNING: You will need to also build and install the plumed library (see http://www.plumed.org) and make sure the file `libplumedKernel.so` (or `libplumedKernel.dylib`) is available on your system.

You should then make sure the library is found setting the environment variable `PLUMED_KERNEL`::

     export PLUMED_KERNEL=/path/to/libplumedKernel.so
     python
     >>> import plumed
     >>> p=plumed.Plumed()

If you manage multiple plumed versions on your system using tcl environment modules, this should be taken care automatically
by the plumed module.

Alternatively, a pure python solution is::

    >>> import plumed
    >>> os.environ["PLUMED_KERNEL"]="/path/to/libplumedKernel.so"
    >>> p=plumed.Plumed()

Finally, notice that you can set the path to the plumed library directly when declaring a Plumed object::

    >>> import plumed
    >>> p=plumed.Plumed(kernel="/path/to/libplumedKernel.so")

This will allow you to mix different plumed versions in the same python script.

CHANGES: See the PLUMED documentation.
Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

WARNING: You will need to also build and install the plumed library (see http://www.plumed.org) and make sure the file `libplumedKernel.so` (or `libplumedKernel.dylib`) is available on your system.

You should then make sure the library is found setting the environment variable `PLUMED_KERNEL`::

     export PLUMED_KERNEL=/path/to/libplumedKernel.so
     python
     >>> import plumed
     >>> p=plumed.Plumed()

If you manage multiple plumed versions on your system using tcl environment modules, this should be taken care automatically
by the plumed module.

Alternatively, a pure python solution is::

    >>> import plumed
    >>> os.environ["PLUMED_KERNEL"]="/path/to/libplumedKernel.so"
    >>> p=plumed.Plumed()

Finally, notice that you can set the path to the plumed library directly when declaring a Plumed object::

    >>> import plumed
    >>> p=plumed.Plumed(kernel="/path/to/libplumedKernel.so")

This will allow you to mix different plumed versions in the same python script.

CHANGES: See the PLUMED documentation.
