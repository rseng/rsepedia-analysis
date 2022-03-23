## Prerequisites

The following are mandatory top-level prerequisites
- C++ compiler with support for the [C++14 standard](https://www.iso.org/standard/64029.html). This includes the following compilers:
  - [GNU C++](https://gcc.gnu.org/), version 5.0 or higher
  - [Clang](https://clang.llvm.org/), version 3.4 or higher
- [CMake](https://cmake.org/), version 3.1 and higher
- [TiledArray](https://github.com/ValeevGroup/tiledarray), source from the master branch
- [Libint](http://libint.valeyev.net), version 2.5.0-beta.1 or higher
- [Boost libraries](www.boost.org/), any recent version should do (e.g. Travis CI tests pass with version 1.55). The following Boost components are used:
  - [Boost.Algorithm](https://www.boost.org/doc/libs/master/libs/algorithm/doc/html/index.html) -- misc algorithms
  - [Boost.Filesystem](https://www.boost.org/doc/libs/master/libs/filesystem/doc/index.htm) -- only used if C++17 support is not enabled; __N.B.__ this is the only Boost library used by MPQC that must be compiled (i.e. it cannot be used as a header-only library)
    - [Boost.System](https://www.boost.org/doc/libs/master/libs/system/doc/index.html) -- this non-header-only library is a prerequisite of Boost.Filesystem, only needed if C++17 support is not available.
  - [Boost.Locale](https://www.boost.org/doc/libs/master/libs/locale/doc/html/index.html) -- to be replaced by the C++ Standard Library facilities
  - [Boost.Math](https://www.boost.org/doc/libs/master/libs/math/doc/html/index.html) -- misc special functions
  - [Boost.Optional](https://www.boost.org/doc/libs/master/libs/optional/doc/html/index.html) -- to be replaced by C++17 std::optional
  - [Boost.PropertyTree](https://www.boost.org/doc/libs/master/doc/html/property_tree.html) -- used to implement KeyVal class
  - [Boost.Serialization](https://www.boost.org/doc/libs/master/libs/serialization/doc/index.html) -- class GUID registration
- Intel Thread Building Blocks (TBB), available in a [commercial](software.intel.com/tbb‎) or
  an [open-source](https://www.threadingbuildingblocks.org/) form
- (for documentation only) Doxygen

The following are transitive dependencies of the above:
- [MADNESS parallel runtime](https://github.com/m-a-d-n-e-s-s/madness)
- [Eigen](http://eigen.tuxfamily.org), version 3.3 or higher
- BLAS and LAPACK libraries

## Compile TiledArray
Please refer to [TiledArry Wiki](https://github.com/ValeevGroup/tiledarray/wiki)
for information on how to compile TiledArray.

## Compile Libint
- obtain the latest version of Libint library from [here](https://github.com/evaleev/libint/releases) (download the library marked "standard ints only")
- compile and install according to [Libint Wiki](https://github.com/evaleev/libint/wiki#compiling-libint-library).

## Compile MPQC

1. Configure cmake
2. Build
3. (optional) Validate
4. Install

### Configuring MPQC

Run cmake from a build directory in which MPQC will be built. For many reasons you should avoid building in the source directory.
There is a number of variables that you may need to provide to `cmake`, such as the C++ compiler to use, the MPI compiler wrapper to use, etc.
At a minimum, you will need to specify the installed locations of the top-level prerequisites, as shown in the example script:

```
#!/bin/bash

export MPQC_SOURCE=/path/to/mpqc/source/code/directory
export TiledArray_INSTALL = /path/to/tiledarray/install/directory
export LIBINT2_INSTALL = /path/to/libint/install/direcotry

cmake \
    -DTiledArray_INSTALL_DIR= ${TiledArray_INSTALL} \
    -DLIBINT2_INSTALL_DIR=${LIBINT2_INSTALL} \
    -DBOOST_ROOT=/path/to/boost/install/direcotry \
    ${MPQC_SOURCE}
```

In practice we recommend using one of the provided *toolchain* files that come with TiledArray to configure MPQC
using exactly the same compiler and library combination that was used to compile TiledArray.

The most useful MPQC-specific `cmake` variables are listed below:

|Variables            |Description|
|---------------------|-----------|
| `TiledArray_INSTALL_DIR` | path to TiledArray install directory |
| `LIBINT2_INSTALL_DIR` | path to Libint2 install directory |
| `BOOST_ROOT` | root path for Boost |
| `TA_POLICY` |  dense or sparse, default sparse. control which policy to use with TiledArray. Some classes may only support sparse |
| `MPQC_VALIDATION_TEST_PRINT` | default off, control if print output after validation test failed |

### Building MPQC
For simplicity here we assume that `cmake` was used to generate UNIX Makefiles (which is the default). To build, validate, and install
MPQC run the
following commands:
- `make`
- (optional) `make check`
- `make install`

## Platform-Specific Notes

### MacOS

#### System Integrity Protection (SIP)
Intel MKL and TBB libraries on MacOS come with user-configurable [RPATH](https://en.wikipedia.org/wiki/Rpath).
Due to the use of [System Integrity Protection](https://support.apple.com/en-us/HT204899)(SIP)
on recent MacOS platforms it is not sufficient to load the appropriate `mklvars`/`tbbvars` scripts in a shell
to allow MPQC find these libraries. Rather than disabling SIP, it is possible to run MPQC with SIP enabled
using a script that sets the `DYLD_LIBRARY_PATH` environment variable before calling MPQC. Store this in an executable script, named `mpqcrun.sh`,
```
#!/bin/sh

source /opt/intel/mkl/bin/mklvars.sh intel64
source /opt/intel/tbb/bin/tbbvars.sh intel64
export DYLD_LIBRARY_PATH=/path/to/tiledarray/install/directory/lib:$DYLD_LIBRARY_PATH

$*
```
and execute MPQC as `mpqcrun.sh /path/to/mpqc/binary <input_file>`. To use this script when validating MPQC,
set environment variable `MPQC_PRE_CMD` to the full path to this script before executing `make check`.

#### Example: MPQC build script

[A simple shell script](https://github.com/ValeevGroup/mpqc/blob/master/bin/osx-brew-build.sh) that uses HomeBrew to install basic prerequisites.
[![Build Status](https://travis-ci.com/ValeevGroup/mpqc4.svg?token=2pDpbi3swi4zsJxpapq7&branch=master)](https://travis-ci.com/ValeevGroup/mpqc4)
[![codecov](https://codecov.io/gh/ValeevGroup/mpqc4/branch/master/graph/badge.svg?token=ErSgG98B8l)](https://codecov.io/gh/ValeevGroup/mpqc4)
# Synopsis

![MPQC Logo](https://github.com/ValeevGroup/mpqc/wiki/images/mpqc_logo_med.png)

The Massively Parallel Quantum Chemistry (MPQC) platform is a research package for *ab initio* simulation of the electronic structure of molecules and periodic solids.
The current (4th) version of the package, sometimes referred to as __MPQC4__, is a modern reenvisioning of the conceptual design of the original MPQC platform using the massively-parallel tensor framework [TiledArray](https://github.com/ValeevGroup/tiledarray), distributed task-based programming model and runtime [MADWorld](https://github.com/m-a-d-n-e-s-s/madness), and the Gaussian integrals library [Libint](https://github.com/evaleev/libint).

# Developers
MPQC is developed by the [Valeev Group](http://research.valeyev.net) at [Virginia Tech](http://www.vt.edu).

# License

MPQC is freely available under the terms of the GPL v3+ licence. See the the included LICENSE file for details. If you are interested in using MPQC under different licensing terms, please contact us.

# How to Cite

See the CITATION file for the recommended citation.

The development of electronic structure methods in MPQC is partially described in the following publications:
* Cannada A. Lewis , Justus A. Calvin , and Edward F. Valeev, "Clustered Low-Rank Tensor Format: Introduction and Application to Fast Construction of Hartree-Fock Exchange.", *J. Chem. Theor. Comp.*, DOI 10.1021/acs.jctc.6b00884;
* Chong Peng, Justus A. Calvin, Fabijan Pavosevic, Jinmei Zhang, and Edward F. Valeev, "Massively Parallel Implementation of Explicitly Correlated Coupled-Cluster Singles and Doubles Using TiledArray Framework.", *J. Phys. Chem. A*, DOI 10.1021/acs.jpca.6b10150.

# Performance

Excellent strong scaling performance of the electronic structure methods in MPQC is demonstrated below for the coupled-cluster singles and doubles (CCSD) wave function solver. Parallel speed-up of 1 iteration of CCSD solver for uracil trimer in 6-31G* AO basis was measured on ["BlueRidge" cluster](https://secure.hosting.vt.edu/www.arc.vt.edu/computing/blueridge-sandy-bridge/) at Virginia Tech (wall time on 1 16-core node = 1290 sec):

![CCSD:UracilTrimer-speedup](https://github.com/ValeevGroup/tiledarray/wiki/images/uracil-trimer-ccsd-blueridge-speedup.png)

This figure was obtained with the help of an allocation from [Advanced Research Computing](https://secure.hosting.vt.edu/www.arc.vt.edu/) at Virginia Tech.

# Acknowledgements
Development of MPQC and its key components is made possible by past and present contributions from the National Science Foundation (awards CHE-0847295, CHE-0741927, OCI-1047696, CHE-1362655, ACI-1450262, and ACI-1550456), the Alfred P. Sloan Foundation, the Camille and Henry Dreyfus Foundation, the Department of Energy Exascale Computing Project (NWChemEx subproject), and the Department of Energy INCITE Program.
This code example, written in C++, demonstrates how to augment the functionality of MPQC by adding a new class that implements MP2 energy.
This serves as an illustration of how to *plug extra functionality into MPQC*, by
extending an existing abstract class; in this case the MPQC LCAOWavefunction class is extended to
implement new quantum chemical methods. The new class can be used alongside
all other MPQC methods to, for example, compute forces by finite differences, etc.

1. To compile this example you must have compiled MPQC and installed it using `make install` command. This will install the code examples under `$prefix/doc/share/mpqc-${version}/examples directory`, where `$prefix` is the installation prefix provided to configure script during configuration of MPQC (the default is `/usr/local/mpqc/$mpqcversion`). The rest of instructions will assume that you are in the `mp2` subdirectory of that directory.
2. Configure the example by typing `cmake .` (for a clean restart remove file `CMakeCache.txt` and directory `CMakeFiles`); then compile the mp2 example by typing `make mp2`. This will create an executable file called `mp2`. To run it type `./mp2 ./mp2.json`.
3. To compare the MP2 energy against that obtained with MPQC, simply change object type from MP2 to RMP2 in `mp2.json`,
   change keyword "localize" to false (since RMP2 implements canonical MP2 energy only), and run the example again as `./mp2 ./mp2.json` (or you can run the main MPQC executable: `mpqc ./mp2.json`).
This code example, written in C++, demonstrates how to implement a simple standalone program, with its own `main()` function, by using the MPQC infrastructure.
This is an illustration of how to *plug MPQC into an existing codebase*;
specificaly, MPQC AO Integral Factories are used to compute and manipulate AO integrals.

1. To compile this example you must have compiled MPQC and installed it using `make install` command. This will install the code examples under `$prefix/doc/share/mpqc-${version}/examples directory`, where `$prefix` is the installation prefix provided to configure script during configuration of MPQC (the default is `/usr/local/mpqc/$mpqcversion`). The rest of instructions will assume that you are in the `aoints` subdirectory of that directory.
2. Configure the example by typing `cmake .`; then compile the aoints example by typing `make`. This will create an executable file called `aoints`. To run it type `./aoints ./aoints.json`.
