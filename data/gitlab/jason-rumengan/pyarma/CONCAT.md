## Changelog:

### v0.500.0 (11 February 2021)  
* Size-only constructors (i.e. mat(rows, cols), mat(size(X)), cube(rows, cols, slices), cube(size(Q))) are initialised with zeros by default
* Conversion of NumPy arrays require the same data type (i.e. a NumPy array of integers cannot be converted to a PyArmadillo mat) (#10)
* Added pyarma_rng.set_seed(value)
* set_seed_random() is now pyarma_rng.set_seed_random()
* linspace(start, end, N) and logspace(start, end, N) take N as an unsigned integer (#12, #13)
* range(X) has been renamed to spread(X) to prevent conflicts with Python's built-in range() function (#14)
* Fixed bug where U, H = hess(X) returns X as well as U, H (#11)
* Added extra forms for lu(X), qr(X), qr_econ(X), qz(A, B), and svd_econ(X)
* Fixed bug where cube(subcube) threw a TypeError (#16)
* Added pyarma_version for version information (#15)
* Added randu(), randn(), zeros(), ones(), eye() generators
* Added subscripting for size objects
* Removed excess newline when printing matrices and cubes
* Removed excess newline printed by libraries() (#17)
* solve_opts_types and fill_types are now solve_opts.types and fill.types (i.e. solve_opts.fast, fill.randu)
* Internal types are hidden
* Fixed bug where clamp(ucube/icube) took wrong argument types (#18)

### v0.400.0 (1 February 2021)  
Initial public release, adapting dense matrices and cubes from Armadillo

### v0.100.0-0.300.0
Internal development releases### PyArmadillo: linear algebra library for Python
https://pyarma.sourceforge.io

Copyright 2020-2021 Jason Rumengan   (https://www.jasonrumengan.my.id)  
Copyright 2020-2021 Terry Yue Zhuo   (https://terryyz.github.io)  
Copyright 2020-2021 Conrad Sanderson (https://conradsanderson.id.au)  
Copyright 2020-2021 Data61 / CSIRO  

[![PyPI version](https://badge.fury.io/py/pyarma.svg)](https://badge.fury.io/py/pyarma)
[![status](https://joss.theoj.org/papers/9e0ec1df49e009869b5a308f6d319c71/status.svg)](https://joss.theoj.org/papers/9e0ec1df49e009869b5a308f6d319c71)
[![Build status](https://ci.appveyor.com/api/projects/status/8hl857p0vwuilnwh?svg=true)](https://ci.appveyor.com/project/jason-rumengan/pyarma)
---

### Quick Links  

- [download latest release](https://pyarma.sourceforge.io)
- [documentation for functions and classes](https://pyarma.sourceforge.io/docs.html)
- [bug reports & questions](https://pyarma.sourceforge.io/faq.html)  

---

### Contents

1. [Introduction](#1-introduction)
2. [Citation Details](#2-citation-details)
3. [Documentation and Examples](#3-documentation-and-examples)

4. [Installation via pip](#4-installation-via-pip)
5. [Installation from Source](#5-installation-from-source)

6. [Distribution License](#6-distribution-license)
7. [Bug Reports and Frequently Asked Questions](#7-bug-reports-and-frequently-asked-questions)

---

### 1. Introduction

* PyArmadillo is a streamlined linear algebra library (matrix maths) for the Python language, with emphasis on ease of use

* Provides high-level syntax and functionality deliberately similar to Matlab

* Provides classes for matrices and cubes; integer, floating point and complex elements are supported 

* Relies on [Armadillo](http://arma.sourceforge.net) for the underlying C++ implementation of matrix objects

* Authors:
    * Jason Rumengan   - https://www.jasonrumengan.my.id
    * Terry Yue Zhuo   - https://terryyz.github.io
    * Conrad Sanderson - http://conradsanderson.id.au

* Example program:    
    ```python
    from pyarma import *

    A = mat(4, 5, fill.ones)
    B = mat(4, 5, fill.randu)
    
    C = A*B.t()
    
    C.print("C:")
    ```
---

### 2: Citation Details

Please cite the following paper if you use PyArmadillo in your research and/or software.  
Citations are useful for the continued development and maintenance of the library.

  * Jason Rumengan, Terry Yue Zhuo, Conrad Sanderson.  
    PyArmadillo: a streamlined linear algebra library for Python.  
    Journal of Open Source Software, Vol. 6(66), p. 3051, 2021.  
    [https://doi.org/10.21105/joss.03051](https://doi.org/10.21105/joss.03051)

---

### 3: Documentation and Examples

The documentation for PyArmadillo functions and classes is available at:  
https://pyarma.sourceforge.io/docs.html

The documentation is also in the `doc/docs.html` file in the PyArmadillo archive,
which can be viewed with a web browser.

A short example program named ```example.py``` that uses PyArmadillo
is included with the PyArmadillo archive.

---

### 4: Installation via pip

* A precompiled version of PyArmadillo is available via the Python Package Index (PyPI)

* Use the following command for installation:  
  `pip3 install --user pyarma`  
  or  
  `pip3 install pyarma`  

* If `pip3` cannot be found, try using the following alternatives:
  - `python3 -m pip`
  - `py -m pip`

* To upgrade PyArmadillo via pip:  
  `pip3 install --upgrade --user pyarma`  
  or  
  `pip3 install --upgrade pyarma`  
  
  **NOTE:** It's possible that pip may erroneously not find the newest version.
  In that case, try the following command:  
  `pip3 install --no-cache-dir --upgrade --user pyarma`  
  or  
  `pip3 install --no-cache-dir --upgrade pyarma`  
  
  

* More info on the pyarma package at PyPI:  
  https://pypi.org/project/pyarma/


---

### 5: Installation from Source

* **Preliminaries**
  
  - Installing PyArmadillo from source requires:
    - at least Python 3.6; the minimum recommended version is Python 3.8
    - a C++ compiler that supports at least the C++11 standard
    - at least 8 GB of RAM
    - 64-bit CPU, preferably with 4+ cores
    - OpenBLAS and LAPACK (or compatible implementations)
  
  - **Linux** based operating systems (eg. Fedora, Ubuntu, CentOS, Red Hat, Debian, etc)
    - First install OpenBLAS, LAPACK, Python 3, and pip3, along with the corresponding development/header files
    - On CentOS 8 / RHEL 8, the CentOS PowerTools repository may need to be enabled:  
      `dnf config-manager --set-enabled powertools`
    - Recommended packages to install before installing PyArmadillo:
      - Fedora, CentOS, RHEL: gcc-c++, libstdc++-devel, openblas-devel, lapack-devel, python3-devel, python3-pip
      - Ubuntu and Debian: g++, libopenblas-dev, liblapack-dev, python3-dev, python3-pip
    - `pip3` needs to be updated:  
      `pip3 install --user --upgrade pip`
  
  - **macOS**
    - First install Xcode (version 8 or later) and then type the following command in a terminal window:  
      `xcode-select --install`
    - Xcode command-line tools include the Python 3 development files, but `pip3` needs to be updated:  
      `pip3 install --user --upgrade pip`
    - The "Accelerate" framework is used for accessing BLAS and LAPACK functions

  - **Windows (x64)**
    - First install Microsoft Visual Studio (2019 or later)
    - Use the x64 Native Tools Command Prompt
    - PyArmadillo contains pre-compiled [OpenBLAS](https://github.com/xianyi/OpenBLAS/releases/),
      which is used for accessing BLAS and LAPACK functions
    - `pip3` needs to be updated:  
      `py -m pip install --user --upgrade pip`
    - Alternative implementations and/or distributions of BLAS and LAPACK are available at:
      - http://software.intel.com/en-us/intel-mkl/
      - http://icl.cs.utk.edu/lapack-for-windows/lapack/
      - http://ylzhao.blogspot.com.au/2013/10/blas-lapack-precompiled-binaries-for.html
    - **Caveat**: 32-bit Windows (x86) is currently not supported

* **Running the Installer**
  
  - Open a terminal window and change into the directory containing PyArmadillo sources
    - if the source was obtained as a package downloaded from SourceForge:
      
          tar xf pyarmadillo-0.123.4.tar.xz
          cd pyarmadillo-0.123.4
        
      (change `0.123.4` to match the downloaded version)
    
    - if the source was obtained by cloning the GitLab repo:
      
          git clone https://gitlab.com/jason-rumengan/pyarma/
          cd pyarma
        
  - Execute the following command:
  
        pip3 install --user .
  
    **NOTE:** the full stop character at the end is important
  
  - To see the progress of installation, change the above command to `pip3 install --verbose --user .`
  
  - If `pip3` cannot be found, try using the following alternatives:
    - `python3 -m pip`
    - `py -m pip`
  
  - Installation may take **5 to 20 minutes** due to compiling C++ sources that extensively use template metaprogramming;
    the time taken depends on the number of CPU cores and the amount of available memory
  
  - **Caveat:** on systems with low memory (&lt; 8 GB), parallel compilation may fail
    due to template metaprogramming requiring large amounts of memory.
    To avoid parallel compilation, first install `scikit-build` using `pip3 install --user scikit-build`
    and then install PyArmadillo using `python setup.py install -- -- -j1`
    
  - Link-time optimisation (LTO) is off by default.
    LTO reduces the size of PyArmadillo at the expense of considerably longer compilation time.
    To enable LTO, first install `scikit-build` and `ninja`,
    and then enable the `PYARMA_LTO` option during installation:
    ```
    pip3 install --user scikit-build ninja
    python3 setup.py install -DPYARMA_LTO=ON
    ```

* **Support for Intel MKL and Other BLAS/LAPACK Implementations**

  - PyArmadillo can optionally use the Intel Math Kernel Library (MKL) as high-speed replacement for standard BLAS and LAPACK
  
  - Intel MKL should be automatically detected during installation from source
  
  - For other BLAS/LAPACK implementations, minor modifications to the built-in Armadillo sources may be required.
    Specifically `ext/armadillo/include/armadillo_bits/config.hpp` may need to be edited
    to ensure Armadillo uses the same integer sizes and style of function names as used by the replacement libraries.
    The following defines may need to be enabled or disabled:
    
        ARMA_BLAS_CAPITALS  
        ARMA_BLAS_UNDERSCORE  
        ARMA_BLAS_LONG  
        ARMA_BLAS_LONG_LONG  
    
    See the [Armadillo](http://arma.sourceforge.net) site for more information:
      - http://arma.sourceforge.net/faq.html
      - http://arma.sourceforge.net/docs.html#config_hpp
  
  - On Linux-based systems, MKL might be installed in a non-standard location such as `/opt` which can cause problems during linking.
    Before installing PyArmadillo, the system should know where the MKL libraries are located.
    For example, `/opt/intel/mkl/lib/intel64/`.
    This can be achieved by setting the `LD_LIBRARY_PATH` environment variable,
    or for a more permanent solution, adding the directory locations to `/etc/ld.so.conf`.
    It may also be possible to store a text file with the locations in the `/etc/ld.so.conf.d` directory.
    For example, `/etc/ld.so.conf.d/mkl.conf`.
    If `/etc/ld.so.conf` is modified or `/etc/ld.so.conf.d/mkl.conf` is created,
    `/sbin/ldconfig` must be run afterwards.  
    Below is an example of `/etc/ld.so.conf.d/mkl.conf` where Intel MKL is installed in `/opt/intel`

          /opt/intel/lib/intel64  
          /opt/intel/mkl/lib/intel64  

  - If MKL is installed and it is persistently giving problems during linking,
    support for MKL can be disabled by editing `ext/armadillo/CMakeLists.txt`
    and commenting out the line containing `INCLUDE(ARMA_FindMKL)`,
    then deleting `ext/armadillo/CMakeCache.txt`,
    and finally re-running PyArmadillo installation.

---

### 6: Distribution License

PyArmadillo is open source software licensed under the Apache License, Version 2.0 (the "License").
A copy of the License is included in the "LICENSE" file.

Any software that incorporates or distributes PyArmadillo in source or binary form
must include, in the documentation and/or other materials provided with the software,
a readable copy of the attribution notices present in the "NOTICE" file.
See the License for details. The contents of the "NOTICE" file are for
informational purposes only and do not modify the License.

---

### 7: Bug Reports and Frequently Asked Questions

If you find a bug in the library or the documentation, we are interested in hearing about it.
Please make a _small_ and _self-contained_ program which exposes the bug,
and then send the program source and the bug description to the developers.

The contact details are at:  
https://pyarma.sourceforge.io/contact.html

Further information about PyArmadillo is on the frequently asked questions page:  
https://pyarma.sourceforge.io/faq.html
### Armadillo: C++ Library for Linear Algebra & Scientific Computing  
http://arma.sourceforge.net

Copyright 2008-2021 Conrad Sanderson (http://conradsanderson.id.au)  
Copyright 2008-2016 National ICT Australia (NICTA)  
Copyright 2017-2021 Data61 / CSIRO  

---

### Quick Links  

- [download latest stable release](http://arma.sourceforge.net/download.html)
- [documentation for functions and classes](http://arma.sourceforge.net/docs.html)
- [bug reports & questions](http://arma.sourceforge.net/faq.html)  

---
 
### Contents

1.  [Introduction](#1-introduction)
2.  [Citation Details](#2-citation-details)
3.  [Distribution License](#3-distribution-license)

4.  [Compilers and External Dependencies](#4-compilers-and-external-dependencies)

5.  [Linux and macOS: Installation](#5-linux-and-macos-installation)
6.  [Linux and macOS: Compiling and Linking](#6-linux-and-macos-compiling-and-linking)

7.  [Windows: Installation](#7-windows-installation)
8.  [Windows: Compiling and Linking](#8-windows-compiling-and-linking)

9.  [Support for OpenBLAS and Intel MKL](#9-support-for-openblas-and-intel-mkl)
10. [Support for ATLAS](#10-support-for-atlas)
11. [Support for OpenMP](#11-support-for-openmp)

12. [Documentation](#12-documentation)
13. [API Stability and Versioning](#13-api-stability-and-versioning)
14. [Bug Reports and Frequently Asked Questions](#14-bug-reports-and-frequently-asked-questions)

15. [MEX Interface to Octave/Matlab](#15-mex-interface-to-octavematlab)
16. [Related Software Using Armadillo](#16-related-software-using-armadillo)

---

### 1. Introduction

Armadillo is a high quality C++ library for linear algebra and scientific computing,
aiming towards a good balance between speed and ease of use.

It's useful for algorithm development directly in C++,
and/or quick conversion of research code into production environments.
It has high-level syntax and functionality which is deliberately similar to Matlab.

The library provides efficient classes for vectors, matrices and cubes,
as well as 200+ associated functions covering essential and advanced functionality
for data processing and manipulation of matrices.

Various matrix decompositions (eigen, SVD, QR, etc) are provided through
integration with LAPACK, or one of its high performance drop-in replacements
(eg. OpenBLAS, Intel MKL, Apple Accelerate framework, etc).

A sophisticated expression evaluator (via C++ template meta-programming)
automatically combines several operations (at compile time) to increase speed
and efficiency.

The library can be used for machine learning, pattern recognition, computer vision,
signal processing, bioinformatics, statistics, finance, etc.

Authors:
  * Conrad Sanderson - http://conradsanderson.id.au
  * Ryan Curtin      - http://ratml.org

---

### 2: Citation Details

Please cite the following papers if you use Armadillo in your research and/or software.  
Citations are useful for the continued development and maintenance of the library.

  * Conrad Sanderson and Ryan Curtin.  
    Armadillo: a template-based C++ library for linear algebra.  
    Journal of Open Source Software, Vol. 1, pp. 26, 2016.  
  
  * Conrad Sanderson and Ryan Curtin.  
    A User-Friendly Hybrid Sparse Matrix Class in C++.  
    Lecture Notes in Computer Science (LNCS), Vol. 10931, pp. 422-430, 2018.

---

### 3: Distribution License

Armadillo can be used in both open-source and proprietary (closed-source) software.

Armadillo is licensed under the Apache License, Version 2.0 (the "License").
A copy of the License is included in the "LICENSE.txt" file.

Any software that incorporates or distributes Armadillo in source or binary form
must include, in the documentation and/or other materials provided with the software,
a readable copy of the attribution notices present in the "NOTICE.txt" file.
See the License for details. The contents of the "NOTICE.txt" file are for
informational purposes only and do not modify the License.

---

### 4: Compilers and External Dependencies

Armadillo 10.x requires a C++ compiler that supports at least the C++11 standard.
Use Armadillo 9.900 if your compiler only supports the old C++98/C++03 standards.

The functionality of Armadillo is partly dependent on other libraries:
LAPACK, BLAS (preferably OpenBLAS), ARPACK and SuperLU.
LAPACK and BLAS are used for dense matrices,
while ARPACK and SuperLU are used for sparse matrices.

Armadillo can work without the above libraries, but its functionality will be reduced.
Basic functionality will be available (eg. matrix addition and multiplication),
but operations like eigen decomposition or matrix inversion will not be.
Matrix multiplication (mainly for big matrices) may not be as fast.

As Armadillo is a template library, we recommended that optimisation
is enabled during compilation of programs that use Armadillo.
For example, for GCC and Clang compilers use -O2 or -O3

---

### 5: Linux and macOS: Installation

* Step 1:
  Ensure a C++ compiler is installed on your system.

  - On macOS systems install Xcode (version 8 or later)
    and then type the following command in a terminal window:  

    xcode-select --install

* Step 2:
  Ensure the CMake tool is installed on your system.

  - Cmake can be downloaded from http://www.cmake.org
    or (preferably) installed using the package manager on your system.

  - On Linux-based systems, CMake can be installed using dnf, yum, apt, aptitude, ...

  - On macOS systems, CMake can be installed through MacPorts or Homebrew.

* Step 3:
  Ensure that OpenBLAS (or standard BLAS and LAPACK) is installed on your system.
  On macOS, the Accelerate framework can be used for BLAS/LAPACK.

  - On macOS, optionally install OpenBLAS for better performance.

  - If support for sparse matrices is required, also install ARPACK and SuperLU.
    Caveat: only SuperLU version 5.2 can be used!

  - On Linux-based systems, the following libraries are recommended
    to be present: OpenBLAS, LAPACK, SuperLU and ARPACK.
    It is also necessary to install the corresponding development files for each library.
    For example, when installing the "libopenblas" package, also install the "libopenblas-dev" package.
  
* Step 4:
  Run the cmake installer.

  - Open a terminal window and change into the directory that was created
    by unpacking the armadillo archive.
  
  - The simplest case is to run cmake using:

    cmake .

  - NOTE: the full stop separated from "cmake" by a space is important.
  
  - Options to the cmake installer:
  
    - On Linux, to enable the detection of FlexiBLAS, 
      use the additional ALLOW_FLEXIBLAS_LINUX option when running cmake:

      cmake -DALLOW_FLEXIBLAS_LINUX=ON .

    - On macOS, to enable the detection of OpenBLAS, 
      use the additional ALLOW_OPENBLAS_MACOS option when running cmake:

      cmake -DALLOW_OPENBLAS_MACOS=ON .

      Note: depending on your installation, OpenBLAS may masquerade as standard BLAS.
      To detect standard BLAS and LAPACK, use the ALLOW_BLAS_LAPACK_MACOS option:

      cmake -DALLOW_BLAS_LAPACK_MACOS=ON .

    - By default, cmake assumes that the Armadillo library and the
      corresponding header files will be installed in the default 
      system directory (eg. in the /usr hierarchy in Linux-based systems).
      To install the library and headers in an alternative directory,
      use the additional option CMAKE_INSTALL_PREFIX in this form:

      cmake . -DCMAKE_INSTALL_PREFIX:PATH=alternative_directory

  - CMake will detect which relevant libraries are installed on your system
    (eg. OpenBLAS, LAPACK, SuperLU, ARPACK, etc)
    and will modify Armadillo's configuration correspondingly.
    CMake will also generate the Armadillo run-time library,
    which is a wrapper for all the detected libraries.

  - If cmake needs to re-run, it's a good idea to first delete the
    "CMakeCache.txt" file (not "CMakeLists.txt").

  - Caveat: if Armadillo is installed in a non-system directory,
    make sure that the C++ compiler is configured to use the "lib" and "include"
    sub-directories present within this directory.  Note that the "lib"
    directory might be named differently on your system.
    On recent 64 bit Debian & Ubuntu systems it is "lib/x86_64-linux-gnu".
    On recent 64 bit Fedora & RHEL systems it is "lib64".

* Step 5:
  If you and have access to root/administrator/superuser privileges
  (ie. able to use "sudo") and didn't use the CMAKE_INSTALL_PREFIX option,
  type the following command:

    sudo make install

  If you don't have root/administrator/superuser privileges,
  make sure that you use the CMAKE_INSTALL_PREFIX option in Step 4,
  and type the following command:

    make install

---

### 6: Linux and macOS: Compiling and Linking

If you have installed Armadillo via the CMake installer,
use the following command:

    g++ prog.cpp -o prog -std=c++11 -O2 -larmadillo

Otherwise, if you want to use Armadillo without installation
(ie. without the Armadillo runtime library), use the following command:
  
    g++ prog.cpp -o prog -std=c++11 -O2 -I /home/blah/armadillo-7.200.3/include -DARMA_DONT_USE_WRAPPER -lopenblas

The above command assumes that the armadillo archive was unpacked into /home/blah/  
The command needs to be adjusted if the archive was unpacked into a different directory
and/or for each specific version of Armadillo (ie. "7.200.3" needs to be changed).
  
If you don't have OpenBLAS, on Linux change -lopenblas to -lblas -llapack
and on macOS change -lopenblas to -framework Accelerate

See the Questions page for more info on linking:
http://arma.sourceforge.net/faq.html

The "examples" directory contains a short example program that uses the Armadillo library.

---

### 7: Windows: Installation

The installation is comprised of 3 steps:

* Step 1:
  Copy the entire "include" folder to a convenient location
  and tell your compiler to use that location for header files
  (in addition to the locations it uses already).
  Alternatively, the "include" folder can be used directly.

* Step 2:
  Modify "include/armadillo_bits/config.hpp" to indicate which
  libraries are currently available on your system. For example,
  if LAPACK, BLAS (or OpenBLAS), ARPACK and SuperLU present,
  uncomment the following lines:

    #define ARMA_USE_LAPACK  
    #define ARMA_USE_BLAS  
    #define ARMA_USE_ARPACK  
    #define ARMA_USE_SUPERLU  

  If support for sparse matrices is not required,
  don't worry about ARPACK or SuperLU.

* Step 3:
  Configure your compiler to link with LAPACK and BLAS
  (and optionally ARPACK and SuperLU).

---

### 8: Windows: Compiling and Linking

Within the "examples" folder, there is an MSVC project named "example1_win64"
which can be used to compile "example1.cpp". The project needs to be compiled as a
64 bit program: the active solution platform must be set to x64, instead of win32.

The MSVC project was tested on Windows 10 (64 bit) with Visual Studio C++ 2019.
Adaptations may need to be made for 32 bit systems, later versions of Windows
and/or the compiler. For example, options such as ARMA_BLAS_LONG and ARMA_BLAS_UNDERSCORE,
defined in "armadillo_bits/config.hpp", may need to be either enabled or disabled.

The folder "examples/lib_win64" contains a copy of lib and dll files
obtained from a pre-compiled release of OpenBLAS 0.3.10:
https://github.com/xianyi/OpenBLAS/releases/download/v0.3.10/OpenBLAS-0.3.10-x64.zip
The compilation was done by a third party.  USE AT YOUR OWN RISK.

**Caveat:** 
for any high performance scientific/engineering workloads,
we strongly recommend using a Linux based operating system:
  * Fedora  http://fedoraproject.org/
  * Ubuntu  http://www.ubuntu.com/
  * CentOS  http://centos.org/

---

### 9: Support for OpenBLAS and Intel MKL

Armadillo can use OpenBLAS or Intel Math Kernel Library (MKL) as high-speed
replacements for BLAS and LAPACK. In essence this involves linking with the
replacement libraries instead of BLAS and LAPACK.

Minor modifications to include/armadillo_bits/config.hpp may be required
to ensure Armadillo uses the same integer sizes and style of function names
as used by the replacement libraries. Specifically, the following defines
may need to be enabled or disabled:

    ARMA_USE_WRAPPER  
    ARMA_BLAS_CAPITALS  
    ARMA_BLAS_UNDERSCORE  
    ARMA_BLAS_LONG  
    ARMA_BLAS_LONG_LONG  

See the documentation for more information on the above defines.

On Linux-based systems, MKL might be installed in a non-standard location
such as /opt which can cause problems during linking.  Before installing
Armadillo, the system should know where the MKL libraries are located.
For example, /opt/intel/mkl/lib/intel64/.  This can be achieved by setting
the LD_LIBRARY_PATH environment variable, or for a more permanent solution,
adding the directory locations to /etc/ld.so.conf.  It may also be possible
to store a text file with the locations in the /etc/ld.so.conf.d directory.
For example, /etc/ld.so.conf.d/mkl.conf.  If /etc/ld.so.conf is modified
or /etc/ld.so.conf.d/mkl.conf is created, /sbin/ldconfig must be run afterwards.

Below is an example of /etc/ld.so.conf.d/mkl.conf
where Intel MKL is installed in /opt/intel

    /opt/intel/lib/intel64  
    /opt/intel/mkl/lib/intel64  

If MKL is installed and it is persistently giving problems during linking,
Support for MKL can be disabled by editing the CMakeLists.txt file,
deleting CMakeCache.txt and re-running the CMake based installation.
Comment out the line containing:

    INCLUDE(ARMA_FindMKL)

---

### 10: Support for ATLAS

Armadillo can use the ATLAS library for faster versions of a subset
of LAPACK and BLAS functions. LAPACK should still be installed to
obtain full functionality.

Caveat: the minimum recommended version of ATLAS is 3.10;
earlier versions (such as 3.6 and 3.8) can produce incorrect
results and/or corrupt memory, leading to random crashes.

---

### 11: Support for OpenMP

Armadillo can use OpenMP to automatically speed up computationally
expensive element-wise functions such as exp(), log(), cos(), etc.
This requires a C++11/C++14 compiler with OpenMP 3.1+ support.

When using gcc or clang, use the following options to enable both
C++11 and OpenMP:  -std=c++11 -fopenmp

---

### 12: Documentation

The documentation for Armadillo functions and classes is available at:  
http://arma.sourceforge.net/docs.html

The documentation is also in the "docs.html" file in this folder,
which can be viewed with a web browser.

---

### 13: API Stability and Versioning

Each release of Armadillo has its public API (functions, classes, constants)
described in the accompanying API documentation (docs.html) specific
to that release.

Each release of Armadillo has its full version specified as A.B.C,
where A is a major version number, B is a minor version number,
and C is a patch level (indicating bug fixes).

Within a major version (eg. 7), each minor version has a public API that
strongly strives to be backwards compatible (at the source level) with the
public API of preceding minor versions. For example, user code written for
version 7.100 should work with version 7.200, 7.300, 7.400, etc. However,
as later minor versions may have more features (API extensions) than
preceding minor versions, user code specifically written for version 7.400
may not work with 7.300.

An increase in the patch level, while the major and minor versions are retained,
indicates modifications to the code and/or documentation which aim to fix bugs
without altering the public API.

We don't like changes to existing public API and strongly prefer not to break
any user software. However, to allow evolution, we reserve the right to
alter the public API in future major versions of Armadillo while remaining
backwards compatible in as many cases as possible (eg. major version 8 may
have slightly different public API than major version 7).

**CAVEAT:** any function, class, constant or other code _not_ explicitly described
in the public API documentation is considered as part of the underlying internal
implementation details, and may change or be removed without notice.
(In other words, don't use internal functionality).

---

### 14: Bug Reports and Frequently Asked Questions

Armadillo has gone through extensive testing and has been successfully
used in production environments. However, as with almost all software,
it's impossible to guarantee 100% correct functionality.

If you find a bug in the library or the documentation, we are interested
in hearing about it. Please make a _small_ and _self-contained_ program
which exposes the bug, and then send the program source and the bug description
to the developers. The small program must have a main() function and use only
functions/classes from Armadillo and the standard C++ library (no other libraries).

The contact details are at:  
http://arma.sourceforge.net/contact.html

Further information about Armadillo is on the frequently asked questions page:  
http://arma.sourceforge.net/faq.html

---

### 15: MEX Interface to Octave/Matlab

The "mex_interface" folder contains examples of how to interface
Octave/Matlab with C++ code that uses Armadillo matrices.

---

### 16: Related Software Using Armadillo

* MLPACK: extensive library of machine learning algorithms  
  http://mlpack.org

* ensmallen: C++ library of numerical optimisation methods  
  http://ensmallen.org/

* SigPack: C++ signal processing library  
  http://sigpack.sourceforge.net

* RcppArmadillo: integration of Armadillo with the R system and environment  
  http://dirk.eddelbuettel.com/code/rcpp.armadillo.html

* PyArmadillo: linear algebra library for Python  
  https://pyarma.sourceforge.io

## Description

<!-- Include relevant issues or PRs here, describe what changed and why -->


## Suggested changelog entry:

<!-- Fill in the below block with the expected RestructuredText entry. Delete if no entry needed;
     but do not delete header or rst block if an entry is needed! Will be collected via a script. -->

```rst

```

<!-- If the upgrade guide needs updating, note that here too -->
Thank you for your interest in this project! Please refer to the following
sections on how to contribute code and bug reports.

### Reporting bugs

Before submitting a question or bug report, please take a moment of your time
and ensure that your issue isn't already discussed in the project documentation
provided at [pybind11.readthedocs.org][] or in the [issue tracker][]. You can
also check [gitter][] to see if it came up before.

Assuming that you have identified a previously unknown problem or an important
question, it's essential that you submit a self-contained and minimal piece of
code that reproduces the problem. In other words: no external dependencies,
isolate the function(s) that cause breakage, submit matched and complete C++
and Python snippets that can be easily compiled and run in isolation; or
ideally make a small PR with a failing test case that can be used as a starting
point.

## Pull requests

Contributions are submitted, reviewed, and accepted using GitHub pull requests.
Please refer to [this article][using pull requests] for details and adhere to
the following rules to make the process as smooth as possible:

* Make a new branch for every feature you're working on.
* Make small and clean pull requests that are easy to review but make sure they
  do add value by themselves.
* Add tests for any new functionality and run the test suite (`cmake --build
  build --target pytest`) to ensure that no existing features break.
* Please run [`pre-commit`][pre-commit] to check your code matches the
  project style. (Note that `gawk` is required.) Use `pre-commit run
  --all-files` before committing (or use installed-mode, check pre-commit docs)
  to verify your code passes before pushing to save time.
* This project has a strong focus on providing general solutions using a
  minimal amount of code, thus small pull requests are greatly preferred.

### Licensing of contributions

pybind11 is provided under a BSD-style license that can be found in the
``LICENSE`` file. By using, distributing, or contributing to this project, you
agree to the terms and conditions of this license.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to the author of this software, without
imposing a separate written license agreement for such Enhancements, then you
hereby grant the following license: a non-exclusive, royalty-free perpetual
license to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.


## Development of pybind11

To setup an ideal development environment, run the following commands on a
system with CMake 3.14+:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r tests/requirements.txt
cmake -S . -B build -DDOWNLOAD_CATCH=ON -DDOWNLOAD_EIGEN=ON
cmake --build build -j4
```

Tips:

* You can use `virtualenv` (from PyPI) instead of `venv` (which is Python 3
  only).
* You can select any name for your environment folder; if it contains "env" it
  will be ignored by git.
* If you don’t have CMake 3.14+, just add “cmake” to the pip install command.
* You can use `-DPYBIND11_FINDPYTHON=ON` to use FindPython on CMake 3.12+
* In classic mode, you may need to set `-DPYTHON_EXECUTABLE=/path/to/python`.
  FindPython uses `-DPython_ROOT_DIR=/path/to` or
  `-DPython_EXECUTABLE=/path/to/python`.

### Configuration options

In CMake, configuration options are given with “-D”. Options are stored in the
build directory, in the `CMakeCache.txt` file, so they are remembered for each
build directory. Two selections are special - the generator, given with `-G`,
and the compiler, which is selected based on environment variables `CXX` and
similar, or `-DCMAKE_CXX_COMPILER=`. Unlike the others, these cannot be changed
after the initial run.

The valid options are:

* `-DCMAKE_BUILD_TYPE`: Release, Debug, MinSizeRel, RelWithDebInfo
* `-DPYBIND11_FINDPYTHON=ON`: Use CMake 3.12+’s FindPython instead of the
  classic, deprecated, custom FindPythonLibs
* `-DPYBIND11_NOPYTHON=ON`: Disable all Python searching (disables tests)
* `-DBUILD_TESTING=ON`: Enable the tests
* `-DDOWNLOAD_CATCH=ON`: Download catch to build the C++ tests
* `-DOWNLOAD_EIGEN=ON`: Download Eigen for the NumPy tests
* `-DPYBIND11_INSTALL=ON/OFF`: Enable the install target (on by default for the
  master project)
* `-DUSE_PYTHON_INSTALL_DIR=ON`: Try to install into the python dir


<details><summary>A few standard CMake tricks: (click to expand)</summary><p>

* Use `cmake --build build -v` to see the commands used to build the files.
* Use `cmake build -LH` to list the CMake options with help.
* Use `ccmake` if available to see a curses (terminal) gui, or `cmake-gui` for
  a completely graphical interface (not present in the PyPI package).
* Use `cmake --build build -j12` to build with 12 cores (for example).
* Use `-G` and the name of a generator to use something different. `cmake
  --help` lists the generators available.
      - On Unix, setting `CMAKE_GENERATER=Ninja` in your environment will give
        you automatic mulithreading on all your CMake projects!
* Open the `CMakeLists.txt` with QtCreator to generate for that IDE.
* You can use `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON` to generate the `.json` file
  that some tools expect.

</p></details>


To run the tests, you can "build" the check target:

```bash
cmake --build build --target check
```

`--target` can be spelled `-t` in CMake 3.15+. You can also run individual
tests with these targets:

* `pytest`: Python tests only, using the
[pytest](https://docs.pytest.org/en/stable/) framework
* `cpptest`: C++ tests only
* `test_cmake_build`: Install / subdirectory tests

If you want to build just a subset of tests, use
`-DPYBIND11_TEST_OVERRIDE="test_callbacks.cpp;test_pickling.cpp"`. If this is
empty, all tests will be built.

You may also pass flags to the `pytest` target by editing `tests/pytest.ini` or
by using the `PYTEST_ADDOPTS` environment variable
(see [`pytest` docs](https://docs.pytest.org/en/2.7.3/customize.html#adding-default-options)). As an example:

```bash
env PYTEST_ADDOPTS="--capture=no --exitfirst" \
    cmake --build build --target pytest
# Or using abbreviated flags
env PYTEST_ADDOPTS="-s -x" cmake --build build --target pytest
```

### Formatting

All formatting is handled by pre-commit.

Install with brew (macOS) or pip (any OS):

```bash
# Any OS
python3 -m pip install pre-commit

# OR macOS with homebrew:
brew install pre-commit
```

Then, you can run it on the items you've added to your staging area, or all
files:

```bash
pre-commit run
# OR
pre-commit run --all-files
```

And, if you want to always use it, you can install it as a git hook (hence the
name, pre-commit):

```bash
pre-commit install
```

### Clang-Format

As of v2.6.2, pybind11 ships with a [`clang-format`][clang-format]
configuration file at the top level of the repo (the filename is
`.clang-format`). Currently, formatting is NOT applied automatically, but
manually using `clang-format` for newly developed files is highly encouraged.
To check if a file needs formatting:

```bash
clang-format -style=file --dry-run some.cpp
```

The output will show things to be fixed, if any. To actually format the file:

```bash
clang-format -style=file -i some.cpp
```

Note that the `-style-file` option searches the parent directories for the
`.clang-format` file, i.e. the commands above can be run in any subdirectory
of the pybind11 repo.

### Clang-Tidy

[`clang-tidy`][clang-tidy] performs deeper static code analyses and is
more complex to run, compared to `clang-format`, but support for `clang-tidy`
is built into the pybind11 CMake configuration. To run `clang-tidy`, the
following recipe should work. Files will be modified in place, so you can
use git to monitor the changes.

```bash
docker run --rm -v $PWD:/pybind11 -it silkeh/clang:10
apt-get update && apt-get install python3-dev python3-pytest
cmake -S pybind11/ -B build -DCMAKE_CXX_CLANG_TIDY="$(which clang-tidy);-fix"
cmake --build build
```

### Include what you use

To run include what you use, install (`brew install include-what-you-use` on
macOS), then run:

```bash
cmake -S . -B build-iwyu -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE=$(which include-what-you-use)
cmake --build build
```

The report is sent to stderr; you can pipe it into a file if you wish.

### Build recipes

This builds with the Intel compiler (assuming it is in your path, along with a
recent CMake and Python 3):

```bash
python3 -m venv venv
. venv/bin/activate
pip install pytest
cmake -S . -B build-intel -DCMAKE_CXX_COMPILER=$(which icpc) -DDOWNLOAD_CATCH=ON -DDOWNLOAD_EIGEN=ON -DPYBIND11_WERROR=ON
```

This will test the PGI compilers:

```bash
docker run --rm -it -v $PWD:/pybind11 nvcr.io/hpc/pgi-compilers:ce
apt-get update && apt-get install -y python3-dev python3-pip python3-pytest
wget -qO- "https://cmake.org/files/v3.18/cmake-3.18.2-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C /usr/local
cmake -S pybind11/ -B build
cmake --build build
```

### Explanation of the SDist/wheel building design

> These details below are _only_ for packaging the Python sources from git. The
> SDists and wheels created do not have any extra requirements at all and are
> completely normal.

The main objective of the packaging system is to create SDists (Python's source
distribution packages) and wheels (Python's binary distribution packages) that
include everything that is needed to work with pybind11, and which can be
installed without any additional dependencies. This is more complex than it
appears: in order to support CMake as a first class language even when using
the PyPI package, they must include the _generated_ CMake files (so as not to
require CMake when installing the `pybind11` package itself). They should also
provide the option to install to the "standard" location
(`<ENVROOT>/include/pybind11` and `<ENVROOT>/share/cmake/pybind11`) so they are
easy to find with CMake, but this can cause problems if you are not an
environment or using ``pyproject.toml`` requirements. This was solved by having
two packages; the "nice" pybind11 package that stores the includes and CMake
files inside the package, that you get access to via functions in the package,
and a `pybind11-global` package that can be included via `pybind11[global]` if
you want the more invasive but discoverable file locations.

If you want to install or package the GitHub source, it is best to have Pip 10
or newer on Windows, macOS, or Linux (manylinux1 compatible, includes most
distributions).  You can then build the SDists, or run any procedure that makes
SDists internally, like making wheels or installing.


```bash
# Editable development install example
python3 -m pip install -e .
```

Since Pip itself does not have an `sdist` command (it does have `wheel` and
`install`), you may want to use the upcoming `build` package:

```bash
python3 -m pip install build

# Normal package
python3 -m build -s .

# Global extra
PYBIND11_GLOBAL_SDIST=1 python3 -m build -s .
```

If you want to use the classic "direct" usage of `python setup.py`, you will
need CMake 3.15+ and either `make` or `ninja` preinstalled (possibly via `pip
install cmake ninja`), since directly running Python on `setup.py` cannot pick
up and install `pyproject.toml` requirements. As long as you have those two
things, though, everything works the way you would expect:

```bash
# Normal package
python3 setup.py sdist

# Global extra
PYBIND11_GLOBAL_SDIST=1 python3 setup.py sdist
```

A detailed explanation of the build procedure design for developers wanting to
work on or maintain the packaging system is as follows:

#### 1. Building from the source directory

When you invoke any `setup.py` command from the source directory, including
`pip wheel .` and `pip install .`, you will activate a full source build. This
is made of the following steps:

1. If the tool is PEP 518 compliant, like Pip 10+, it will create a temporary
   virtual environment and install the build requirements (mostly CMake) into
   it. (if you are not on Windows, macOS, or a manylinux compliant system, you
   can disable this with `--no-build-isolation` as long as you have CMake 3.15+
   installed)
2. The environment variable `PYBIND11_GLOBAL_SDIST` is checked - if it is set
   and truthy, this will be make the accessory `pybind11-global` package,
   instead of the normal `pybind11` package. This package is used for
   installing the files directly to your environment root directory, using
   `pybind11[global]`.
2. `setup.py` reads the version from `pybind11/_version.py` and verifies it
   matches `includes/pybind11/detail/common.h`.
3. CMake is run with `-DCMAKE_INSTALL_PREIFX=pybind11`. Since the CMake install
   procedure uses only relative paths and is identical on all platforms, these
   files are valid as long as they stay in the correct relative position to the
   includes. `pybind11/share/cmake/pybind11` has the CMake files, and
   `pybind11/include` has the includes. The build directory is discarded.
4. Simpler files are placed in the SDist: `tools/setup_*.py.in`,
   `tools/pyproject.toml` (`main` or `global`)
5. The package is created by running the setup function in the
   `tools/setup_*.py`.  `setup_main.py` fills in Python packages, and
   `setup_global.py` fills in only the data/header slots.
6. A context manager cleans up the temporary CMake install directory (even if
   an error is thrown).

### 2. Building from SDist

Since the SDist has the rendered template files in `tools` along with the
includes and CMake files in the correct locations, the builds are completely
trivial and simple. No extra requirements are required. You can even use Pip 9
if you really want to.


[pre-commit]: https://pre-commit.com
[clang-format]: https://clang.llvm.org/docs/ClangFormat.html
[clang-tidy]: https://clang.llvm.org/extra/clang-tidy/
[pybind11.readthedocs.org]: http://pybind11.readthedocs.org/en/latest
[issue tracker]: https://github.com/pybind/pybind11/issues
[gitter]: https://gitter.im/pybind/Lobby
[using pull requests]: https://help.github.com/articles/using-pull-requests
---
name: Feature Request
about: File an issue about adding a feature
title: "[FEAT] "
---


Make sure you've completed the following steps before submitting your issue -- thank you!

1. Check if your feature has already been mentioned / rejected / planned in other issues.
2. If those resources didn't help, consider asking in the [Gitter chat room][] to see if this is interesting / useful to a larger audience and possible to implement reasonably,
4. If you have a useful feature that passes the previous items (or not suitable for chat), please fill in the details below.

[Gitter chat room]: https://gitter.im/pybind/Lobby

*After reading, remove this checklist.*
---
name: Question
about: File an issue about unexplained behavior
title: "[QUESTION] "
---

If you have a question, please check the following first:

1. Check if your question has already been answered in the [FAQ][] section.
2. Make sure you've read the [documentation][]. Your issue may be addressed there.
3. If those resources didn't help and you only have a short question (not a bug report), consider asking in the [Gitter chat room][]
4. Search the [issue tracker][], including the closed issues, to see if your question has already been asked/answered. +1 or comment if it has been asked but has no answer.
5. If you have a more complex question which is not answered in the previous items (or not suitable for chat), please fill in the details below.
6. Include a self-contained and minimal piece of code that illustrates your question. If that's not possible, try to make the description as clear as possible.

[FAQ]: http://pybind11.readthedocs.io/en/latest/faq.html
[documentation]: https://pybind11.readthedocs.io
[issue tracker]: https://github.com/pybind/pybind11/issues
[Gitter chat room]: https://gitter.im/pybind/Lobby

*After reading, remove this checklist.*
---
name: Bug Report
about: File an issue about a bug
title: "[BUG] "
---


Make sure you've completed the following steps before submitting your issue -- thank you!

1. Make sure you've read the [documentation][]. Your issue may be addressed there.
2. Search the [issue tracker][] to verify that this hasn't already been reported. +1 or comment there if it has.
3. Consider asking first in the [Gitter chat room][].
4. Include a self-contained and minimal piece of code that reproduces the problem. If that's not possible, try to make the description as clear as possible.
    a. If possible, make a PR with a new, failing test to give us a starting point to work on!

[documentation]: https://pybind11.readthedocs.io
[issue tracker]: https://github.com/pybind/pybind11/issues
[Gitter chat room]: https://gitter.im/pybind/Lobby

*After reading, remove this checklist and the template text in parentheses below.*

## Issue description

(Provide a short description, state the expected behavior and what actually happens.)

## Reproducible example code

(The code should be minimal, have no external dependencies, isolate the function(s) that cause breakage. Submit matched and complete C++ and Python snippets that can be easily compiled and run to diagnose the issue.)
---
title: 'PyArmadillo: a streamlined linear algebra library for Python'
tags:
    - linear algebra
    - scientific computing
    - mathematics
    - Python
authors:
    - name: Jason Rumengan
      orcid: 0000-0003-1839-5138
      affiliation: "1, 2"
    - name: Terry Yue Zhuo
      orcid: 0000-0002-5760-5188
      affiliation: 3
    - name: Conrad Sanderson
      orcid: 0000-0002-0049-4501
      affiliation: "1, 4"
affiliations:
    - name: Data61/CSIRO, Australia
      index: 1
    - name: Queensland University of Technology, Australia
      index: 2
    - name: University of New South Wales, Australia
      index: 3
    - name: Griffith University, Australia
      index: 4
date: 12 October 2021
bibliography: paper.bib
---


# Summary

PyArmadillo is a linear algebra library for the Python language,
with the aim of closely mirroring the programming interface of the widely used Armadillo C++ library, 
which in turn is deliberately similar to Matlab.
PyArmadillo hence facilitates algorithm prototyping with Matlab-like syntax directly in Python,
and relatively straightforward conversion of PyArmadillo-based Python code
into performant Armadillo-based C++ code.
The converted code can be used for purposes such as speeding up Python-based programs
in conjunction with pybind11 [@pybind11],
or the integration of algorithms originally prototyped in Python into larger C++ codebases.

PyArmadillo provides objects for matrices and cubes,
as well as over 200 associated functions for manipulating data stored in the objects.
Integer, floating point and complex numbers are supported.
Various matrix factorisations are provided through integration with LAPACK [@anderson1999lapack],
or one of its high performance drop-in replacements such as Intel MKL [@Intel_MKL] or OpenBLAS [@Xianyi_OpenBLAS].


# Statement of Need

Armadillo is a popular linear algebra and scientific computing library for the C++ language [@Sanderson_2016; @Sanderson_2018]
that has three main characteristics:
(i) a high-level programming interface deliberately similar to Matlab,
(ii) an expression evaluator (based on template meta-programming)
that automatically combines several operations to increase speed and efficiency,
and
(iii) an efficient mapper between mathematical expressions and low-level BLAS/LAPACK functions [@Psarras_2021].
Matlab is widely used in both industrial and academic contexts,
providing a programming interface that allows mathematical expressions to be written in a concise and natural manner [@Linge_MatlabOctave_2016], 
especially in comparison to directly using low-level libraries such as LAPACK [@anderson1999lapack].
In industrial settings, algorithms are often first prototyped in Matlab,
before conversion into another language, such as C++, for the purpose of integration into products.
The similarity of the programming interfaces between Armadillo and Matlab
facilitates direct prototyping in C++,
as well as the conversion of research code into production environments.
Armadillo is also often used for implementing performance critical parts of software packages
running under the R environment for statistical computing [@R_manual], 
via the RcppArmadillo bridge [@Eddelbuettel_2014].

Over the past few years, Python has become popular for data science and machine learning.
This partly stems from a rich ecosystem of supporting frameworks and packages,
as well as lack of licensing costs in comparison to Matlab.
Python allows relatively quick prototyping of algorithms,
aided by its dynamically typed nature
and the interpreted execution of user code,
avoiding time-consuming compilation into machine code.
However, for the joint purpose of algorithm prototyping and deployment,
the flexibility of Python comes with two main issues:
(i) slow execution speed due to the interpreted nature of the language,
(ii) difficulty with integration of code written in Python into larger programs and/or frameworks written in another language.
The first issue can be somewhat addressed through conversion of Python-based code into the low-level Cython language [@Cython].
However, since Cython is closely tied with Python, 
conversion of Python code into C++ may be preferred as it also addresses the second issue,
as well as providing a higher-level of abstraction.

PyArmadillo is aimed at:
(i) users that prefer compact Matlab-like syntax rather than the somewhat more verbose syntax provided by NumPy/SciPy [@harris2020array; @2020SciPy-NMeth],
and
(ii) users that would like a straightforward conversion path to performant C++ code.
More specifically,
PyArmadillo aims to closely mirror the programming interface of the Armadillo library,
thereby facilitating the prototyping of algorithms with Matlab-like syntax directly in Python.
Furthermore, PyArmadillo-based Python code can be easily converted
into high-performance Armadillo-based C++ code.
Due to the similarity of the programming interfaces,
the risk of introducing bugs in the conversion process is considerably reduced.
Moreover, conversion into C++ based code allows taking advantage of expression optimisation
performed at compile-time by Armadillo, resulting in further speedups.
The resulting code can be used in larger C++ programs, 
or used as a replacement of performance critical parts within a Python program
with the aid of the pybind11 interface layer [@pybind11].



# Functionality

PyArmadillo provides matrix objects for several distinct element types: integers, single- and double-precision floating point numbers, as well as complex numbers.
In addition to matrices, PyArmadillo also has support for cubes (3 dimensional arrays), where each cube can be treated as an ordered set of matrices.
Multi-dimensional arrays beyond 3 dimensions are explicitly beyond the scope of PyArmadillo.
Over 200 functions are provided for manipulating data stored in the objects, covering the following areas:
fundamental arithmetic operations, 
contiguous and non-contiguous submatrix views,
diagonal views, 
element-wise functions,
scalar/vector/matrix valued functions of matrices,
generation of various vectors/matrices,
statistics,
signal processing,
storage of matrices in files,
matrix decompositions/factorisations,
matrix inverses,
and equation solvers.
See the online documentation at [https://pyarma.sourceforge.io/docs.html](https://pyarma.sourceforge.io/docs.html) for details.
PyArmadillo matrices and cubes are convertible to/from NumPy arrays,
allowing users to tap into the wider Python data science ecosystem,
including plotting tools such as Matplotlib [@Hunter2007].

# Implementation

PyArmadillo relies on pybind11 [@pybind11] for interfacing C++ and Python,
as well as on Armadillo for the underlying C++ implementation of matrix objects and associated functions.
Due to its expressiveness and relatively straightforward use,
pybind11 was selected over other interfacing approaches such as Boost.Python [@osti_815409] and manually writing C++ extensions for Python.
In turn, Armadillo interfaces with low-level routines in BLAS and LAPACK [@anderson1999lapack],
where BLAS is used for matrix multiplication,
and LAPACK is used for various matrix decompositions/factorisations and equation solvers.
As the low-level routines in BLAS and LAPACK are considered as the _de facto_ standard for numerical linear algebra,
it is possible to use high performance drop-in replacements such as Intel MKL [@Intel_MKL] and OpenBLAS [@Xianyi_OpenBLAS].

PyArmadillo is open-source software, distributed under the Apache 2.0 license [@apache],
making it useful in both open-source and proprietary (closed-source) contexts [@Laurent_2008].
It can be obtained at [https://pyarma.sourceforge.io](https://pyarma.sourceforge.io)
or via the Python Package Index in precompiled form.


# Acknowledgements

We would like to thank our colleagues at Data61/CSIRO
(Dan Pagendam, Dan Gladish, Andrew Bolt, Piotr Szul)
for providing feedback and testing.


# References
.. figure:: https://github.com/pybind/pybind11/raw/master/docs/pybind11-logo.png
   :alt: pybind11 logo

**pybind11 — Seamless operability between C++11 and Python**

|Latest Documentation Status| |Stable Documentation Status| |Gitter chat| |CI| |Build status|

|Repology| |PyPI package| |Conda-forge| |Python Versions|

`Setuptools example <https://github.com/pybind/python_example>`_
• `Scikit-build example <https://github.com/pybind/scikit_build_example>`_
• `CMake example <https://github.com/pybind/cmake_example>`_

.. start

.. warning::

   Combining older versions of pybind11 (< 2.6.0) with Python 3.9.0 will
   trigger undefined behavior that typically manifests as crashes during
   interpreter shutdown (but could also destroy your data. **You have been
   warned.**)

   We recommend that you update to the latest patch release of Python (3.9.1),
   which includes a `fix <https://github.com/python/cpython/pull/22670>`_
   that resolves this problem. If you do use Python 3.9.0, please update to
   the latest version of pybind11 (2.6.0 or newer), which includes a temporary
   workaround specifically when Python 3.9.0 is detected at runtime.


**pybind11** is a lightweight header-only library that exposes C++ types
in Python and vice versa, mainly to create Python bindings of existing
C++ code. Its goals and syntax are similar to the excellent
`Boost.Python <http://www.boost.org/doc/libs/1_58_0/libs/python/doc/>`_
library by David Abrahams: to minimize boilerplate code in traditional
extension modules by inferring type information using compile-time
introspection.

The main issue with Boost.Python—and the reason for creating such a
similar project—is Boost. Boost is an enormously large and complex suite
of utility libraries that works with almost every C++ compiler in
existence. This compatibility has its cost: arcane template tricks and
workarounds are necessary to support the oldest and buggiest of compiler
specimens. Now that C++11-compatible compilers are widely available,
this heavy machinery has become an excessively large and unnecessary
dependency.

Think of this library as a tiny self-contained version of Boost.Python
with everything stripped away that isn’t relevant for binding
generation. Without comments, the core header files only require ~4K
lines of code and depend on Python (2.7 or 3.5+, or PyPy) and the C++
standard library. This compact implementation was possible thanks to
some of the new C++11 language features (specifically: tuples, lambda
functions and variadic templates). Since its creation, this library has
grown beyond Boost.Python in many ways, leading to dramatically simpler
binding code in many common situations.

Tutorial and reference documentation is provided at
`pybind11.readthedocs.io <https://pybind11.readthedocs.io/en/latest>`_.
A PDF version of the manual is available
`here <https://pybind11.readthedocs.io/_/downloads/en/latest/pdf/>`_.
And the source code is always available at
`github.com/pybind/pybind11 <https://github.com/pybind/pybind11>`_.


Core features
-------------


pybind11 can map the following core C++ features to Python:

- Functions accepting and returning custom data structures per value,
  reference, or pointer
- Instance methods and static methods
- Overloaded functions
- Instance attributes and static attributes
- Arbitrary exception types
- Enumerations
- Callbacks
- Iterators and ranges
- Custom operators
- Single and multiple inheritance
- STL data structures
- Smart pointers with reference counting like ``std::shared_ptr``
- Internal references with correct reference counting
- C++ classes with virtual (and pure virtual) methods can be extended
  in Python

Goodies
-------

In addition to the core functionality, pybind11 provides some extra
goodies:

- Python 2.7, 3.5+, and PyPy/PyPy3 7.3 are supported with an
  implementation-agnostic interface.

- It is possible to bind C++11 lambda functions with captured
  variables. The lambda capture data is stored inside the resulting
  Python function object.

- pybind11 uses C++11 move constructors and move assignment operators
  whenever possible to efficiently transfer custom data types.

- It’s easy to expose the internal storage of custom data types through
  Pythons’ buffer protocols. This is handy e.g. for fast conversion
  between C++ matrix classes like Eigen and NumPy without expensive
  copy operations.

- pybind11 can automatically vectorize functions so that they are
  transparently applied to all entries of one or more NumPy array
  arguments.

- Python’s slice-based access and assignment operations can be
  supported with just a few lines of code.

- Everything is contained in just a few header files; there is no need
  to link against any additional libraries.

- Binaries are generally smaller by a factor of at least 2 compared to
  equivalent bindings generated by Boost.Python. A recent pybind11
  conversion of PyRosetta, an enormous Boost.Python binding project,
  `reported <http://graylab.jhu.edu/RosettaCon2016/PyRosetta-4.pdf>`_
  a binary size reduction of **5.4x** and compile time reduction by
  **5.8x**.

- Function signatures are precomputed at compile time (using
  ``constexpr``), leading to smaller binaries.

- With little extra effort, C++ types can be pickled and unpickled
  similar to regular Python objects.

Supported compilers
-------------------

1. Clang/LLVM 3.3 or newer (for Apple Xcode’s clang, this is 5.0.0 or
   newer)
2. GCC 4.8 or newer
3. Microsoft Visual Studio 2015 Update 3 or newer
4. Intel classic C++ compiler 18 or newer (ICC 20.2 tested in CI)
5. Cygwin/GCC (previously tested on 2.5.1)
6. NVCC (CUDA 11.0 tested in CI)
7. NVIDIA PGI (20.9 tested in CI)

About
-----

This project was created by `Wenzel
Jakob <http://rgl.epfl.ch/people/wjakob>`_. Significant features and/or
improvements to the code were contributed by Jonas Adler, Lori A. Burns,
Sylvain Corlay, Eric Cousineau, Ralf Grosse-Kunstleve, Trent Houliston, Axel
Huebl, @hulucc, Yannick Jadoul, Sergey Lyskov Johan Mabille, Tomasz Miąsko,
Dean Moldovan, Ben Pritchard, Jason Rhinelander, Boris Schäling,  Pim
Schellart, Henry Schreiner, Ivan Smirnov, Boris Staletic, and Patrick Stewart.

We thank Google for a generous financial contribution to the continuous
integration infrastructure used by this project.


Contributing
~~~~~~~~~~~~

See the `contributing
guide <https://github.com/pybind/pybind11/blob/master/.github/CONTRIBUTING.md>`_
for information on building and contributing to pybind11.

License
~~~~~~~

pybind11 is provided under a BSD-style license that can be found in the
`LICENSE <https://github.com/pybind/pybind11/blob/master/LICENSE>`_
file. By using, distributing, or contributing to this project, you agree
to the terms and conditions of this license.

.. |Latest Documentation Status| image:: https://readthedocs.org/projects/pybind11/badge?version=latest
   :target: http://pybind11.readthedocs.org/en/latest
.. |Stable Documentation Status| image:: https://img.shields.io/badge/docs-stable-blue.svg
   :target: http://pybind11.readthedocs.org/en/stable
.. |Gitter chat| image:: https://img.shields.io/gitter/room/gitterHQ/gitter.svg
   :target: https://gitter.im/pybind/Lobby
.. |CI| image:: https://github.com/pybind/pybind11/workflows/CI/badge.svg
   :target: https://github.com/pybind/pybind11/actions
.. |Build status| image:: https://ci.appveyor.com/api/projects/status/riaj54pn4h08xy40?svg=true
   :target: https://ci.appveyor.com/project/wjakob/pybind11
.. |PyPI package| image:: https://img.shields.io/pypi/v/pybind11.svg
   :target: https://pypi.org/project/pybind11/
.. |Conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/pybind11.svg
   :target: https://github.com/conda-forge/pybind11-feedstock
.. |Repology| image:: https://repology.org/badge/latest-versions/python:pybind11.svg
   :target: https://repology.org/project/python:pybind11/versions
.. |Python Versions| image:: https://img.shields.io/pypi/pyversions/pybind11.svg
   :target: https://pypi.org/project/pybind11/
Upgrade guide
#############

This is a companion guide to the :doc:`changelog`. While the changelog briefly
lists all of the new features, improvements and bug fixes, this upgrade guide
focuses only the subset which directly impacts your experience when upgrading
to a new version. But it goes into more detail. This includes things like
deprecated APIs and their replacements, build system changes, general code
modernization and other useful information.

.. _upgrade-guide-2.6:

v2.6
====

Usage of the ``PYBIND11_OVERLOAD*`` macros and ``get_overload`` function should
be replaced by ``PYBIND11_OVERRIDE*`` and ``get_override``. In the future, the
old macros may be deprecated and removed.

``py::module`` has been renamed ``py::module_``, but a backward compatible
typedef has been included. This change was to avoid a language change in C++20
that requires unqualified ``module`` not be placed at the start of a logical
line. Qualified usage is unaffected and the typedef will remain unless the
C++ language rules change again.

The public constructors of ``py::module_`` have been deprecated. Use
``PYBIND11_MODULE`` or ``module_::create_extension_module`` instead.

An error is now thrown when ``__init__`` is forgotten on subclasses. This was
incorrect before, but was not checked. Add a call to ``__init__`` if it is
missing.

A ``py::type_error`` is now thrown when casting to a subclass (like
``py::bytes`` from ``py::object``) if the conversion is not valid. Make a valid
conversion instead.

The undocumented ``h.get_type()`` method has been deprecated and replaced by
``py::type::of(h)``.

Enums now have a ``__str__`` method pre-defined; if you want to override it,
the simplest fix is to add the new ``py::prepend()`` tag when defining
``"__str__"``.

If ``__eq__`` defined but not ``__hash__``, ``__hash__`` is now set to
``None``, as in normal CPython. You should add ``__hash__`` if you intended the
class to be hashable, possibly using the new ``py::hash`` shortcut.

The constructors for ``py::array`` now always take signed integers for size,
for consistency. This may lead to compiler warnings on some systems. Cast to
``py::ssize_t`` instead of ``std::size_t``.

The ``tools/clang`` submodule and ``tools/mkdoc.py`` have been moved to a
standalone package, `pybind11-mkdoc`_. If you were using those tools, please
use them via a pip install from the new location.

The ``pybind11`` package on PyPI no longer fills the wheel "headers" slot - if
you were using the headers from this slot, they are available by requesting the
``global`` extra, that is, ``pip install "pybind11[global]"``. (Most users will
be unaffected, as the ``pybind11/include`` location is reported by ``python -m
pybind11 --includes`` and ``pybind11.get_include()`` is still correct and has
not changed since 2.5).

.. _pybind11-mkdoc: https://github.com/pybind/pybind11-mkdoc

CMake support:
--------------

The minimum required version of CMake is now 3.4.  Several details of the CMake
support have been deprecated; warnings will be shown if you need to change
something. The changes are:

* ``PYBIND11_CPP_STANDARD=<platform-flag>`` is deprecated, please use
  ``CMAKE_CXX_STANDARD=<number>`` instead, or any other valid CMake CXX or CUDA
  standard selection method, like ``target_compile_features``.

* If you do not request a standard, pybind11 targets will compile with the
  compiler default, but not less than C++11, instead of forcing C++14 always.
  If you depend on the old behavior, please use ``set(CMAKE_CXX_STANDARD 14 CACHE STRING "")``
  instead.

* Direct ``pybind11::module`` usage should always be accompanied by at least
  ``set(CMAKE_CXX_VISIBILITY_PRESET hidden)`` or similar - it used to try to
  manually force this compiler flag (but not correctly on all compilers or with
  CUDA).

* ``pybind11_add_module``'s ``SYSTEM`` argument is deprecated and does nothing;
  linking now behaves like other imported libraries consistently in both
  config and submodule mode, and behaves like a ``SYSTEM`` library by
  default.

* If ``PYTHON_EXECUTABLE`` is not set, virtual environments (``venv``,
  ``virtualenv``, and ``conda``) are prioritized over the standard search
  (similar to the new FindPython mode).

In addition, the following changes may be of interest:

* ``CMAKE_INTERPROCEDURAL_OPTIMIZATION`` will be respected by
  ``pybind11_add_module`` if set instead of linking to ``pybind11::lto`` or
  ``pybind11::thin_lto``.

* Using ``find_package(Python COMPONENTS Interpreter Development)`` before
  pybind11 will cause pybind11 to use the new Python mechanisms instead of its
  own custom search, based on a patched version of classic ``FindPythonInterp``
  / ``FindPythonLibs``. In the future, this may become the default. A recent
  (3.15+ or 3.18.2+) version of CMake is recommended.



v2.5
====

The Python package now includes the headers as data in the package itself, as
well as in the "headers" wheel slot. ``pybind11 --includes`` and
``pybind11.get_include()`` report the new location, which is always correct
regardless of how pybind11 was installed, making the old ``user=`` argument
meaningless. If you are not using the function to get the location already, you
are encouraged to switch to the package location.


v2.2
====

Deprecation of the ``PYBIND11_PLUGIN`` macro
--------------------------------------------

``PYBIND11_MODULE`` is now the preferred way to create module entry points.
The old macro emits a compile-time deprecation warning.

.. code-block:: cpp

    // old
    PYBIND11_PLUGIN(example) {
        py::module m("example", "documentation string");

        m.def("add", [](int a, int b) { return a + b; });

        return m.ptr();
    }

    // new
    PYBIND11_MODULE(example, m) {
        m.doc() = "documentation string"; // optional

        m.def("add", [](int a, int b) { return a + b; });
    }


New API for defining custom constructors and pickling functions
---------------------------------------------------------------

The old placement-new custom constructors have been deprecated. The new approach
uses ``py::init()`` and factory functions to greatly improve type safety.

Placement-new can be called accidentally with an incompatible type (without any
compiler errors or warnings), or it can initialize the same object multiple times
if not careful with the Python-side ``__init__`` calls. The new-style custom
constructors prevent such mistakes. See :ref:`custom_constructors` for details.

.. code-block:: cpp

    // old -- deprecated (runtime warning shown only in debug mode)
    py::class<Foo>(m, "Foo")
        .def("__init__", [](Foo &self, ...) {
            new (&self) Foo(...); // uses placement-new
        });

    // new
    py::class<Foo>(m, "Foo")
        .def(py::init([](...) { // Note: no `self` argument
            return new Foo(...); // return by raw pointer
            // or: return std::make_unique<Foo>(...); // return by holder
            // or: return Foo(...); // return by value (move constructor)
        }));

Mirroring the custom constructor changes, ``py::pickle()`` is now the preferred
way to get and set object state. See :ref:`pickling` for details.

.. code-block:: cpp

    // old -- deprecated (runtime warning shown only in debug mode)
    py::class<Foo>(m, "Foo")
        ...
        .def("__getstate__", [](const Foo &self) {
            return py::make_tuple(self.value1(), self.value2(), ...);
        })
        .def("__setstate__", [](Foo &self, py::tuple t) {
            new (&self) Foo(t[0].cast<std::string>(), ...);
        });

    // new
    py::class<Foo>(m, "Foo")
        ...
        .def(py::pickle(
            [](const Foo &self) { // __getstate__
                return py::make_tuple(self.value1(), self.value2(), ...); // unchanged
            },
            [](py::tuple t) { // __setstate__, note: no `self` argument
                return new Foo(t[0].cast<std::string>(), ...);
                // or: return std::make_unique<Foo>(...); // return by holder
                // or: return Foo(...); // return by value (move constructor)
            }
        ));

For both the constructors and pickling, warnings are shown at module
initialization time (on import, not when the functions are called).
They're only visible when compiled in debug mode. Sample warning:

.. code-block:: none

    pybind11-bound class 'mymodule.Foo' is using an old-style placement-new '__init__'
    which has been deprecated. See the upgrade guide in pybind11's docs.


Stricter enforcement of hidden symbol visibility for pybind11 modules
---------------------------------------------------------------------

pybind11 now tries to actively enforce hidden symbol visibility for modules.
If you're using either one of pybind11's :doc:`CMake or Python build systems
<compiling>` (the two example repositories) and you haven't been exporting any
symbols, there's nothing to be concerned about. All the changes have been done
transparently in the background. If you were building manually or relied on
specific default visibility, read on.

Setting default symbol visibility to *hidden* has always been recommended for
pybind11 (see :ref:`faq:symhidden`). On Linux and macOS, hidden symbol
visibility (in conjunction with the ``strip`` utility) yields much smaller
module binaries. `CPython's extension docs`_ also recommend hiding symbols
by default, with the goal of avoiding symbol name clashes between modules.
Starting with v2.2, pybind11 enforces this more strictly: (1) by declaring
all symbols inside the ``pybind11`` namespace as hidden and (2) by including
the ``-fvisibility=hidden`` flag on Linux and macOS (only for extension
modules, not for embedding the interpreter).

.. _CPython's extension docs: https://docs.python.org/3/extending/extending.html#providing-a-c-api-for-an-extension-module

The namespace-scope hidden visibility is done automatically in pybind11's
headers and it's generally transparent to users. It ensures that:

* Modules compiled with different pybind11 versions don't clash with each other.

* Some new features, like ``py::module_local`` bindings, can work as intended.

The ``-fvisibility=hidden`` flag applies the same visibility to user bindings
outside of the ``pybind11`` namespace. It's now set automatic by pybind11's
CMake and Python build systems, but this needs to be done manually by users
of other build systems. Adding this flag:

* Minimizes the chances of symbol conflicts between modules. E.g. if two
  unrelated modules were statically linked to different (ABI-incompatible)
  versions of the same third-party library, a symbol clash would be likely
  (and would end with unpredictable results).

* Produces smaller binaries on Linux and macOS, as pointed out previously.

Within pybind11's CMake build system, ``pybind11_add_module`` has always been
setting the ``-fvisibility=hidden`` flag in release mode. From now on, it's
being applied unconditionally, even in debug mode and it can no longer be opted
out of with the ``NO_EXTRAS`` option. The ``pybind11::module`` target now also
adds this flag to it's interface. The ``pybind11::embed`` target is unchanged.

The most significant change here is for the ``pybind11::module`` target. If you
were previously relying on default visibility, i.e. if your Python module was
doubling as a shared library with dependents, you'll need to either export
symbols manually (recommended for cross-platform libraries) or factor out the
shared library (and have the Python module link to it like the other
dependents). As a temporary workaround, you can also restore default visibility
using the CMake code below, but this is not recommended in the long run:

.. code-block:: cmake

    target_link_libraries(mymodule PRIVATE pybind11::module)

    add_library(restore_default_visibility INTERFACE)
    target_compile_options(restore_default_visibility INTERFACE -fvisibility=default)
    target_link_libraries(mymodule PRIVATE restore_default_visibility)


Local STL container bindings
----------------------------

Previous pybind11 versions could only bind types globally -- all pybind11
modules, even unrelated ones, would have access to the same exported types.
However, this would also result in a conflict if two modules exported the
same C++ type, which is especially problematic for very common types, e.g.
``std::vector<int>``. :ref:`module_local` were added to resolve this (see
that section for a complete usage guide).

``py::class_`` still defaults to global bindings (because these types are
usually unique across modules), however in order to avoid clashes of opaque
types, ``py::bind_vector`` and ``py::bind_map`` will now bind STL containers
as ``py::module_local`` if their elements are: builtins (``int``, ``float``,
etc.), not bound using ``py::class_``, or bound as ``py::module_local``. For
example, this change allows multiple modules to bind ``std::vector<int>``
without causing conflicts. See :ref:`stl_bind` for more details.

When upgrading to this version, if you have multiple modules which depend on
a single global binding of an STL container, note that all modules can still
accept foreign  ``py::module_local`` types in the direction of Python-to-C++.
The locality only affects the C++-to-Python direction. If this is needed in
multiple modules, you'll need to either:

* Add a copy of the same STL binding to all of the modules which need it.

* Restore the global status of that single binding by marking it
  ``py::module_local(false)``.

The latter is an easy workaround, but in the long run it would be best to
localize all common type bindings in order to avoid conflicts with
third-party modules.


Negative strides for Python buffer objects and numpy arrays
-----------------------------------------------------------

Support for negative strides required changing the integer type from unsigned
to signed in the interfaces of ``py::buffer_info`` and ``py::array``. If you
have compiler warnings enabled, you may notice some new conversion warnings
after upgrading. These can be resolved using ``static_cast``.


Deprecation of some ``py::object`` APIs
---------------------------------------

To compare ``py::object`` instances by pointer, you should now use
``obj1.is(obj2)`` which is equivalent to ``obj1 is obj2`` in Python.
Previously, pybind11 used ``operator==`` for this (``obj1 == obj2``), but
that could be confusing and is now deprecated (so that it can eventually
be replaced with proper rich object comparison in a future release).

For classes which inherit from ``py::object``, ``borrowed`` and ``stolen``
were previously available as protected constructor tags. Now the types
should be used directly instead: ``borrowed_t{}`` and ``stolen_t{}``
(`#771 <https://github.com/pybind/pybind11/pull/771>`_).


Stricter compile-time error checking
------------------------------------

Some error checks have been moved from run time to compile time. Notably,
automatic conversion of ``std::shared_ptr<T>`` is not possible when ``T`` is
not directly registered with ``py::class_<T>`` (e.g. ``std::shared_ptr<int>``
or ``std::shared_ptr<std::vector<T>>`` are not automatically convertible).
Attempting to bind a function with such arguments now results in a compile-time
error instead of waiting to fail at run time.

``py::init<...>()`` constructor definitions are also stricter and now prevent
bindings which could cause unexpected behavior:

.. code-block:: cpp

    struct Example {
        Example(int &);
    };

    py::class_<Example>(m, "Example")
        .def(py::init<int &>()); // OK, exact match
        // .def(py::init<int>()); // compile-time error, mismatch

A non-``const`` lvalue reference is not allowed to bind to an rvalue. However,
note that a constructor taking ``const T &`` can still be registered using
``py::init<T>()`` because a ``const`` lvalue reference can bind to an rvalue.

v2.1
====

Minimum compiler versions are enforced at compile time
------------------------------------------------------

The minimums also apply to v2.0 but the check is now explicit and a compile-time
error is raised if the compiler does not meet the requirements:

* GCC >= 4.8
* clang >= 3.3 (appleclang >= 5.0)
* MSVC >= 2015u3
* Intel C++ >= 15.0


The ``py::metaclass`` attribute is not required for static properties
---------------------------------------------------------------------

Binding classes with static properties is now possible by default. The
zero-parameter version of ``py::metaclass()`` is deprecated. However, a new
one-parameter ``py::metaclass(python_type)`` version was added for rare
cases when a custom metaclass is needed to override pybind11's default.

.. code-block:: cpp

    // old -- emits a deprecation warning
    py::class_<Foo>(m, "Foo", py::metaclass())
        .def_property_readonly_static("foo", ...);

    // new -- static properties work without the attribute
    py::class_<Foo>(m, "Foo")
        .def_property_readonly_static("foo", ...);

    // new -- advanced feature, override pybind11's default metaclass
    py::class_<Bar>(m, "Bar", py::metaclass(custom_python_type))
        ...


v2.0
====

Breaking changes in ``py::class_``
----------------------------------

These changes were necessary to make type definitions in pybind11
future-proof, to support PyPy via its ``cpyext`` mechanism (`#527
<https://github.com/pybind/pybind11/pull/527>`_), and to improve efficiency
(`rev. 86d825 <https://github.com/pybind/pybind11/commit/86d825>`_).

1. Declarations of types that provide access via the buffer protocol must
   now include the ``py::buffer_protocol()`` annotation as an argument to
   the ``py::class_`` constructor.

   .. code-block:: cpp

       py::class_<Matrix>("Matrix", py::buffer_protocol())
           .def(py::init<...>())
           .def_buffer(...);

2. Classes which include static properties (e.g. ``def_readwrite_static()``)
   must now include the ``py::metaclass()`` attribute. Note: this requirement
   has since been removed in v2.1. If you're upgrading from 1.x, it's
   recommended to skip directly to v2.1 or newer.

3. This version of pybind11 uses a redesigned mechanism for instantiating
   trampoline classes that are used to override virtual methods from within
   Python. This led to the following user-visible syntax change:

   .. code-block:: cpp

       // old v1.x syntax
       py::class_<TrampolineClass>("MyClass")
           .alias<MyClass>()
           ...

       // new v2.x syntax
       py::class_<MyClass, TrampolineClass>("MyClass")
           ...

   Importantly, both the original and the trampoline class are now specified
   as arguments to the ``py::class_`` template, and the ``alias<..>()`` call
   is gone. The new scheme has zero overhead in cases when Python doesn't
   override any functions of the underlying C++ class.
   `rev. 86d825 <https://github.com/pybind/pybind11/commit/86d825>`_.

   The class type must be the first template argument given to ``py::class_``
   while the trampoline can be mixed in arbitrary order with other arguments
   (see the following section).


Deprecation of the ``py::base<T>()`` attribute
----------------------------------------------

``py::base<T>()`` was deprecated in favor of specifying ``T`` as a template
argument to ``py::class_``. This new syntax also supports multiple inheritance.
Note that, while the type being exported must be the first argument in the
``py::class_<Class, ...>`` template, the order of the following types (bases,
holder and/or trampoline) is not important.

.. code-block:: cpp

    // old v1.x
    py::class_<Derived>("Derived", py::base<Base>());

    // new v2.x
    py::class_<Derived, Base>("Derived");

    // new -- multiple inheritance
    py::class_<Derived, Base1, Base2>("Derived");

    // new -- apart from `Derived` the argument order can be arbitrary
    py::class_<Derived, Base1, Holder, Base2, Trampoline>("Derived");


Out-of-the-box support for ``std::shared_ptr``
----------------------------------------------

The relevant type caster is now built in, so it's no longer necessary to
include a declaration of the form:

.. code-block:: cpp

    PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

Continuing to do so won’t cause an error or even a deprecation warning,
but it's completely redundant.


Deprecation of a few ``py::object`` APIs
----------------------------------------

All of the old-style calls emit deprecation warnings.

+---------------------------------------+---------------------------------------------+
|  Old syntax                           |  New syntax                                 |
+=======================================+=============================================+
| ``obj.call(args...)``                 | ``obj(args...)``                            |
+---------------------------------------+---------------------------------------------+
| ``obj.str()``                         | ``py::str(obj)``                            |
+---------------------------------------+---------------------------------------------+
| ``auto l = py::list(obj); l.check()`` | ``py::isinstance<py::list>(obj)``           |
+---------------------------------------+---------------------------------------------+
| ``py::object(ptr, true)``             | ``py::reinterpret_borrow<py::object>(ptr)`` |
+---------------------------------------+---------------------------------------------+
| ``py::object(ptr, false)``            | ``py::reinterpret_steal<py::object>(ptr)``  |
+---------------------------------------+---------------------------------------------+
| ``if (obj.attr("foo"))``              | ``if (py::hasattr(obj, "foo"))``            |
+---------------------------------------+---------------------------------------------+
| ``if (obj["bar"])``                   | ``if (obj.contains("bar"))``                |
+---------------------------------------+---------------------------------------------+
.. _installing:

Installing the library
######################

There are several ways to get the pybind11 source, which lives at
`pybind/pybind11 on GitHub <https://github.com/pybind/pybind11>`_. The pybind11
developers recommend one of the first three ways listed here, submodule, PyPI,
or conda-forge, for obtaining pybind11.

.. _include_as_a_submodule:

Include as a submodule
======================

When you are working on a project in Git, you can use the pybind11 repository
as a submodule. From your git repository, use:

.. code-block:: bash

    git submodule add -b stable ../../pybind/pybind11 extern/pybind11
    git submodule update --init

This assumes you are placing your dependencies in ``extern/``, and that you are
using GitHub; if you are not using GitHub, use the full https or ssh URL
instead of the relative URL ``../../pybind/pybind11`` above. Some other servers
also require the ``.git`` extension (GitHub does not).

From here, you can now include ``extern/pybind11/include``, or you can use
the various integration tools (see :ref:`compiling`) pybind11 provides directly
from the local folder.

Include with PyPI
=================

You can download the sources and CMake files as a Python package from PyPI
using Pip. Just use:

.. code-block:: bash

    pip install pybind11

This will provide pybind11 in a standard Python package format. If you want
pybind11 available directly in your environment root, you can use:

.. code-block:: bash

    pip install "pybind11[global]"

This is not recommended if you are installing with your system Python, as it
will add files to ``/usr/local/include/pybind11`` and
``/usr/local/share/cmake/pybind11``, so unless that is what you want, it is
recommended only for use in virtual environments or your ``pyproject.toml``
file (see :ref:`compiling`).

Include with conda-forge
========================

You can use pybind11 with conda packaging via `conda-forge
<https://github.com/conda-forge/pybind11-feedstock>`_:

.. code-block:: bash

    conda install -c conda-forge pybind11


Include with vcpkg
==================
You can download and install pybind11 using the Microsoft `vcpkg
<https://github.com/Microsoft/vcpkg/>`_ dependency manager:

.. code-block:: bash

    git clone https://github.com/Microsoft/vcpkg.git
    cd vcpkg
    ./bootstrap-vcpkg.sh
    ./vcpkg integrate install
    vcpkg install pybind11

The pybind11 port in vcpkg is kept up to date by Microsoft team members and
community contributors. If the version is out of date, please `create an issue
or pull request <https://github.com/Microsoft/vcpkg/>`_ on the vcpkg
repository.

Global install with brew
========================

The brew package manager (Homebrew on macOS, or Linuxbrew on Linux) has a
`pybind11 package
<https://github.com/Homebrew/homebrew-core/blob/master/Formula/pybind11.rb>`_.
To install:

.. code-block:: bash

    brew install pybind11

.. We should list Conan, and possibly a few other C++ package managers (hunter,
.. perhaps). Conan has a very clean CMake integration that would be good to show.

Other options
=============

Other locations you can find pybind11 are `listed here
<https://repology.org/project/python:pybind11/versions>`_; these are maintained
by various packagers and the community.
.. _reference:

.. warning::

    Please be advised that the reference documentation discussing pybind11
    internals is currently incomplete. Please refer to the previous sections
    and the pybind11 header files for the nitty gritty details.

Reference
#########

.. _macros:

Macros
======

.. doxygendefine:: PYBIND11_MODULE

.. _core_types:

Convenience classes for arbitrary Python types
==============================================

Common member functions
-----------------------

.. doxygenclass:: object_api
    :members:

Without reference counting
--------------------------

.. doxygenclass:: handle
    :members:

With reference counting
-----------------------

.. doxygenclass:: object
    :members:

.. doxygenfunction:: reinterpret_borrow

.. doxygenfunction:: reinterpret_steal

Convenience classes for specific Python types
=============================================

.. doxygenclass:: module_
    :members:

.. doxygengroup:: pytypes
    :members:

.. _extras:

Passing extra arguments to ``def`` or ``class_``
================================================

.. doxygengroup:: annotations
    :members:

Embedding the interpreter
=========================

.. doxygendefine:: PYBIND11_EMBEDDED_MODULE

.. doxygenfunction:: initialize_interpreter

.. doxygenfunction:: finalize_interpreter

.. doxygenclass:: scoped_interpreter

Redirecting C++ streams
=======================

.. doxygenclass:: scoped_ostream_redirect

.. doxygenclass:: scoped_estream_redirect

.. doxygenfunction:: add_ostream_redirect

Python built-in functions
=========================

.. doxygengroup:: python_builtins
    :members:

Inheritance
===========

See :doc:`/classes` and :doc:`/advanced/classes` for more detail.

.. doxygendefine:: PYBIND11_OVERRIDE

.. doxygendefine:: PYBIND11_OVERRIDE_PURE

.. doxygendefine:: PYBIND11_OVERRIDE_NAME

.. doxygendefine:: PYBIND11_OVERRIDE_PURE_NAME

.. doxygenfunction:: get_override

Exceptions
==========

.. doxygenclass:: error_already_set
    :members:

.. doxygenclass:: builtin_exception
    :members:


Literals
========

.. doxygennamespace:: literals
.. _classes:

Object-oriented code
####################

Creating bindings for a custom type
===================================

Let's now look at a more complex example where we'll create bindings for a
custom C++ data structure named ``Pet``. Its definition is given below:

.. code-block:: cpp

    struct Pet {
        Pet(const std::string &name) : name(name) { }
        void setName(const std::string &name_) { name = name_; }
        const std::string &getName() const { return name; }

        std::string name;
    };

The binding code for ``Pet`` looks as follows:

.. code-block:: cpp

    #include <pybind11/pybind11.h>

    namespace py = pybind11;

    PYBIND11_MODULE(example, m) {
        py::class_<Pet>(m, "Pet")
            .def(py::init<const std::string &>())
            .def("setName", &Pet::setName)
            .def("getName", &Pet::getName);
    }

:class:`class_` creates bindings for a C++ *class* or *struct*-style data
structure. :func:`init` is a convenience function that takes the types of a
constructor's parameters as template arguments and wraps the corresponding
constructor (see the :ref:`custom_constructors` section for details). An
interactive Python session demonstrating this example is shown below:

.. code-block:: pycon

    % python
    >>> import example
    >>> p = example.Pet('Molly')
    >>> print(p)
    <example.Pet object at 0x10cd98060>
    >>> p.getName()
    u'Molly'
    >>> p.setName('Charly')
    >>> p.getName()
    u'Charly'

.. seealso::

    Static member functions can be bound in the same way using
    :func:`class_::def_static`.

Keyword and default arguments
=============================
It is possible to specify keyword and default arguments using the syntax
discussed in the previous chapter. Refer to the sections :ref:`keyword_args`
and :ref:`default_args` for details.

Binding lambda functions
========================

Note how ``print(p)`` produced a rather useless summary of our data structure in the example above:

.. code-block:: pycon

    >>> print(p)
    <example.Pet object at 0x10cd98060>

To address this, we could bind a utility function that returns a human-readable
summary to the special method slot named ``__repr__``. Unfortunately, there is no
suitable functionality in the ``Pet`` data structure, and it would be nice if
we did not have to change it. This can easily be accomplished by binding a
Lambda function instead:

.. code-block:: cpp

        py::class_<Pet>(m, "Pet")
            .def(py::init<const std::string &>())
            .def("setName", &Pet::setName)
            .def("getName", &Pet::getName)
            .def("__repr__",
                [](const Pet &a) {
                    return "<example.Pet named '" + a.name + "'>";
                }
            );

Both stateless [#f1]_ and stateful lambda closures are supported by pybind11.
With the above change, the same Python code now produces the following output:

.. code-block:: pycon

    >>> print(p)
    <example.Pet named 'Molly'>

.. [#f1] Stateless closures are those with an empty pair of brackets ``[]`` as the capture object.

.. _properties:

Instance and static fields
==========================

We can also directly expose the ``name`` field using the
:func:`class_::def_readwrite` method. A similar :func:`class_::def_readonly`
method also exists for ``const`` fields.

.. code-block:: cpp

        py::class_<Pet>(m, "Pet")
            .def(py::init<const std::string &>())
            .def_readwrite("name", &Pet::name)
            // ... remainder ...

This makes it possible to write

.. code-block:: pycon

    >>> p = example.Pet('Molly')
    >>> p.name
    u'Molly'
    >>> p.name = 'Charly'
    >>> p.name
    u'Charly'

Now suppose that ``Pet::name`` was a private internal variable
that can only be accessed via setters and getters.

.. code-block:: cpp

    class Pet {
    public:
        Pet(const std::string &name) : name(name) { }
        void setName(const std::string &name_) { name = name_; }
        const std::string &getName() const { return name; }
    private:
        std::string name;
    };

In this case, the method :func:`class_::def_property`
(:func:`class_::def_property_readonly` for read-only data) can be used to
provide a field-like interface within Python that will transparently call
the setter and getter functions:

.. code-block:: cpp

        py::class_<Pet>(m, "Pet")
            .def(py::init<const std::string &>())
            .def_property("name", &Pet::getName, &Pet::setName)
            // ... remainder ...

Write only properties can be defined by passing ``nullptr`` as the
input for the read function.

.. seealso::

    Similar functions :func:`class_::def_readwrite_static`,
    :func:`class_::def_readonly_static` :func:`class_::def_property_static`,
    and :func:`class_::def_property_readonly_static` are provided for binding
    static variables and properties. Please also see the section on
    :ref:`static_properties` in the advanced part of the documentation.

Dynamic attributes
==================

Native Python classes can pick up new attributes dynamically:

.. code-block:: pycon

    >>> class Pet:
    ...     name = 'Molly'
    ...
    >>> p = Pet()
    >>> p.name = 'Charly'  # overwrite existing
    >>> p.age = 2  # dynamically add a new attribute

By default, classes exported from C++ do not support this and the only writable
attributes are the ones explicitly defined using :func:`class_::def_readwrite`
or :func:`class_::def_property`.

.. code-block:: cpp

    py::class_<Pet>(m, "Pet")
        .def(py::init<>())
        .def_readwrite("name", &Pet::name);

Trying to set any other attribute results in an error:

.. code-block:: pycon

    >>> p = example.Pet()
    >>> p.name = 'Charly'  # OK, attribute defined in C++
    >>> p.age = 2  # fail
    AttributeError: 'Pet' object has no attribute 'age'

To enable dynamic attributes for C++ classes, the :class:`py::dynamic_attr` tag
must be added to the :class:`py::class_` constructor:

.. code-block:: cpp

    py::class_<Pet>(m, "Pet", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("name", &Pet::name);

Now everything works as expected:

.. code-block:: pycon

    >>> p = example.Pet()
    >>> p.name = 'Charly'  # OK, overwrite value in C++
    >>> p.age = 2  # OK, dynamically add a new attribute
    >>> p.__dict__  # just like a native Python class
    {'age': 2}

Note that there is a small runtime cost for a class with dynamic attributes.
Not only because of the addition of a ``__dict__``, but also because of more
expensive garbage collection tracking which must be activated to resolve
possible circular references. Native Python classes incur this same cost by
default, so this is not anything to worry about. By default, pybind11 classes
are more efficient than native Python classes. Enabling dynamic attributes
just brings them on par.

.. _inheritance:

Inheritance and automatic downcasting
=====================================

Suppose now that the example consists of two data structures with an
inheritance relationship:

.. code-block:: cpp

    struct Pet {
        Pet(const std::string &name) : name(name) { }
        std::string name;
    };

    struct Dog : Pet {
        Dog(const std::string &name) : Pet(name) { }
        std::string bark() const { return "woof!"; }
    };

There are two different ways of indicating a hierarchical relationship to
pybind11: the first specifies the C++ base class as an extra template
parameter of the :class:`class_`:

.. code-block:: cpp

    py::class_<Pet>(m, "Pet")
       .def(py::init<const std::string &>())
       .def_readwrite("name", &Pet::name);

    // Method 1: template parameter:
    py::class_<Dog, Pet /* <- specify C++ parent type */>(m, "Dog")
        .def(py::init<const std::string &>())
        .def("bark", &Dog::bark);

Alternatively, we can also assign a name to the previously bound ``Pet``
:class:`class_` object and reference it when binding the ``Dog`` class:

.. code-block:: cpp

    py::class_<Pet> pet(m, "Pet");
    pet.def(py::init<const std::string &>())
       .def_readwrite("name", &Pet::name);

    // Method 2: pass parent class_ object:
    py::class_<Dog>(m, "Dog", pet /* <- specify Python parent type */)
        .def(py::init<const std::string &>())
        .def("bark", &Dog::bark);

Functionality-wise, both approaches are equivalent. Afterwards, instances will
expose fields and methods of both types:

.. code-block:: pycon

    >>> p = example.Dog('Molly')
    >>> p.name
    u'Molly'
    >>> p.bark()
    u'woof!'

The C++ classes defined above are regular non-polymorphic types with an
inheritance relationship. This is reflected in Python:

.. code-block:: cpp

    // Return a base pointer to a derived instance
    m.def("pet_store", []() { return std::unique_ptr<Pet>(new Dog("Molly")); });

.. code-block:: pycon

    >>> p = example.pet_store()
    >>> type(p)  # `Dog` instance behind `Pet` pointer
    Pet          # no pointer downcasting for regular non-polymorphic types
    >>> p.bark()
    AttributeError: 'Pet' object has no attribute 'bark'

The function returned a ``Dog`` instance, but because it's a non-polymorphic
type behind a base pointer, Python only sees a ``Pet``. In C++, a type is only
considered polymorphic if it has at least one virtual function and pybind11
will automatically recognize this:

.. code-block:: cpp

    struct PolymorphicPet {
        virtual ~PolymorphicPet() = default;
    };

    struct PolymorphicDog : PolymorphicPet {
        std::string bark() const { return "woof!"; }
    };

    // Same binding code
    py::class_<PolymorphicPet>(m, "PolymorphicPet");
    py::class_<PolymorphicDog, PolymorphicPet>(m, "PolymorphicDog")
        .def(py::init<>())
        .def("bark", &PolymorphicDog::bark);

    // Again, return a base pointer to a derived instance
    m.def("pet_store2", []() { return std::unique_ptr<PolymorphicPet>(new PolymorphicDog); });

.. code-block:: pycon

    >>> p = example.pet_store2()
    >>> type(p)
    PolymorphicDog  # automatically downcast
    >>> p.bark()
    u'woof!'

Given a pointer to a polymorphic base, pybind11 performs automatic downcasting
to the actual derived type. Note that this goes beyond the usual situation in
C++: we don't just get access to the virtual functions of the base, we get the
concrete derived type including functions and attributes that the base type may
not even be aware of.

.. seealso::

    For more information about polymorphic behavior see :ref:`overriding_virtuals`.


Overloaded methods
==================

Sometimes there are several overloaded C++ methods with the same name taking
different kinds of input arguments:

.. code-block:: cpp

    struct Pet {
        Pet(const std::string &name, int age) : name(name), age(age) { }

        void set(int age_) { age = age_; }
        void set(const std::string &name_) { name = name_; }

        std::string name;
        int age;
    };

Attempting to bind ``Pet::set`` will cause an error since the compiler does not
know which method the user intended to select. We can disambiguate by casting
them to function pointers. Binding multiple functions to the same Python name
automatically creates a chain of function overloads that will be tried in
sequence.

.. code-block:: cpp

    py::class_<Pet>(m, "Pet")
       .def(py::init<const std::string &, int>())
       .def("set", static_cast<void (Pet::*)(int)>(&Pet::set), "Set the pet's age")
       .def("set", static_cast<void (Pet::*)(const std::string &)>(&Pet::set), "Set the pet's name");

The overload signatures are also visible in the method's docstring:

.. code-block:: pycon

    >>> help(example.Pet)

    class Pet(__builtin__.object)
     |  Methods defined here:
     |
     |  __init__(...)
     |      Signature : (Pet, str, int) -> NoneType
     |
     |  set(...)
     |      1. Signature : (Pet, int) -> NoneType
     |
     |      Set the pet's age
     |
     |      2. Signature : (Pet, str) -> NoneType
     |
     |      Set the pet's name

If you have a C++14 compatible compiler [#cpp14]_, you can use an alternative
syntax to cast the overloaded function:

.. code-block:: cpp

    py::class_<Pet>(m, "Pet")
        .def("set", py::overload_cast<int>(&Pet::set), "Set the pet's age")
        .def("set", py::overload_cast<const std::string &>(&Pet::set), "Set the pet's name");

Here, ``py::overload_cast`` only requires the parameter types to be specified.
The return type and class are deduced. This avoids the additional noise of
``void (Pet::*)()`` as seen in the raw cast. If a function is overloaded based
on constness, the ``py::const_`` tag should be used:

.. code-block:: cpp

    struct Widget {
        int foo(int x, float y);
        int foo(int x, float y) const;
    };

    py::class_<Widget>(m, "Widget")
       .def("foo_mutable", py::overload_cast<int, float>(&Widget::foo))
       .def("foo_const",   py::overload_cast<int, float>(&Widget::foo, py::const_));

If you prefer the ``py::overload_cast`` syntax but have a C++11 compatible compiler only,
you can use ``py::detail::overload_cast_impl`` with an additional set of parentheses:

.. code-block:: cpp

    template <typename... Args>
    using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

    py::class_<Pet>(m, "Pet")
        .def("set", overload_cast_<int>()(&Pet::set), "Set the pet's age")
        .def("set", overload_cast_<const std::string &>()(&Pet::set), "Set the pet's name");

.. [#cpp14] A compiler which supports the ``-std=c++14`` flag
            or Visual Studio 2015 Update 2 and newer.

.. note::

    To define multiple overloaded constructors, simply declare one after the
    other using the ``.def(py::init<...>())`` syntax. The existing machinery
    for specifying keyword and default arguments also works.

Enumerations and internal types
===============================

Let's now suppose that the example class contains an internal enumeration type,
e.g.:

.. code-block:: cpp

    struct Pet {
        enum Kind {
            Dog = 0,
            Cat
        };

        Pet(const std::string &name, Kind type) : name(name), type(type) { }

        std::string name;
        Kind type;
    };

The binding code for this example looks as follows:

.. code-block:: cpp

    py::class_<Pet> pet(m, "Pet");

    pet.def(py::init<const std::string &, Pet::Kind>())
        .def_readwrite("name", &Pet::name)
        .def_readwrite("type", &Pet::type);

    py::enum_<Pet::Kind>(pet, "Kind")
        .value("Dog", Pet::Kind::Dog)
        .value("Cat", Pet::Kind::Cat)
        .export_values();

To ensure that the ``Kind`` type is created within the scope of ``Pet``, the
``pet`` :class:`class_` instance must be supplied to the :class:`enum_`.
constructor. The :func:`enum_::export_values` function exports the enum entries
into the parent scope, which should be skipped for newer C++11-style strongly
typed enums.

.. code-block:: pycon

    >>> p = Pet('Lucy', Pet.Cat)
    >>> p.type
    Kind.Cat
    >>> int(p.type)
    1L

The entries defined by the enumeration type are exposed in the ``__members__`` property:

.. code-block:: pycon

    >>> Pet.Kind.__members__
    {'Dog': Kind.Dog, 'Cat': Kind.Cat}

The ``name`` property returns the name of the enum value as a unicode string.

.. note::

    It is also possible to use ``str(enum)``, however these accomplish different
    goals. The following shows how these two approaches differ.

    .. code-block:: pycon

        >>> p = Pet( "Lucy", Pet.Cat )
        >>> pet_type = p.type
        >>> pet_type
        Pet.Cat
        >>> str(pet_type)
        'Pet.Cat'
        >>> pet_type.name
        'Cat'

.. note::

    When the special tag ``py::arithmetic()`` is specified to the ``enum_``
    constructor, pybind11 creates an enumeration that also supports rudimentary
    arithmetic and bit-level operations like comparisons, and, or, xor, negation,
    etc.

    .. code-block:: cpp

        py::enum_<Pet::Kind>(pet, "Kind", py::arithmetic())
           ...

    By default, these are omitted to conserve space.
Limitations
###########

Design choices
^^^^^^^^^^^^^^

pybind11 strives to be a general solution to binding generation, but it also has
certain limitations:

- pybind11 casts away ``const``-ness in function arguments and return values.
  This is in line with the Python language, which has no concept of ``const``
  values. This means that some additional care is needed to avoid bugs that
  would be caught by the type checker in a traditional C++ program.

- The NumPy interface ``pybind11::array`` greatly simplifies accessing
  numerical data from C++ (and vice versa), but it's not a full-blown array
  class like ``Eigen::Array`` or ``boost.multi_array``. ``Eigen`` objects are
  directly supported, however, with ``pybind11/eigen.h``.

Large but useful features could be implemented in pybind11 but would lead to a
significant increase in complexity. Pybind11 strives to be simple and compact.
Users who require large new features are encouraged to write an extension to
pybind11; see `pybind11_json <https://github.com/pybind/pybind11_json>`_ for an
example.


Known bugs
^^^^^^^^^^

These are issues that hopefully will one day be fixed, but currently are
unsolved. If you know how to help with one of these issues, contributions
are welcome!

- Intel 20.2 is currently having an issue with the test suite.
  `#2573 <https://github.com/pybind/pybind11/pull/2573>`_

- Debug mode Python does not support 1-5 tests in the test suite currently.
  `#2422 <https://github.com/pybind/pybind11/pull/2422>`_

- PyPy3 7.3.1 and 7.3.2 have issues with several tests on 32-bit Windows.

Known limitations
^^^^^^^^^^^^^^^^^

These are issues that are probably solvable, but have not been fixed yet. A
clean, well written patch would likely be accepted to solve them.

- Type casters are not kept alive recursively.
  `#2527 <https://github.com/pybind/pybind11/issues/2527>`_
  One consequence is that containers of ``char *`` are currently not supported.
  `#2245 <https://github.com/pybind/pybind11/issues/2245>`_

- The ``cpptest`` does not run on Windows with Python 3.8 or newer, due to DLL
  loader changes. User code that is correctly installed should not be affected.
  `#2560 <https://github.com/pybind/pybind11/issue/2560>`_

Python 3.9.0 warning
^^^^^^^^^^^^^^^^^^^^

Combining older versions of pybind11 (< 2.6.0) with Python on 3.9.0 will
trigger undefined behavior that typically manifests as crashes during
interpreter shutdown (but could also destroy your data. **You have been
warned**).

This issue has been
`fixed in Python <https://github.com/python/cpython/pull/22670>`_.  As a
mitigation until 3.9.1 is released and commonly used, pybind11 (2.6.0 or newer)
includes a temporary workaround specifically when Python 3.9.0 is detected at
runtime, leaking about 50 bytes of memory when a callback function is garbage
collected. For reference; the pybind11 test suite has about 2,000 such
callbacks, but only 49 are garbage collected before the end-of-process. Wheels
built with Python 3.9.0 will correctly avoid the leak when run in Python 3.9.1.
.. _basics:

First steps
###########

This sections demonstrates the basic features of pybind11. Before getting
started, make sure that development environment is set up to compile the
included set of test cases.


Compiling the test cases
========================

Linux/macOS
-----------

On Linux  you'll need to install the **python-dev** or **python3-dev** packages as
well as **cmake**. On macOS, the included python version works out of the box,
but **cmake** must still be installed.

After installing the prerequisites, run

.. code-block:: bash

   mkdir build
   cd build
   cmake ..
   make check -j 4

The last line will both compile and run the tests.

Windows
-------

On Windows, only **Visual Studio 2015** and newer are supported since pybind11 relies
on various C++11 language features that break older versions of Visual Studio.

.. Note::

    To use the C++17 in Visual Studio 2017 (MSVC 14.1), pybind11 requires the flag
    ``/permissive-`` to be passed to the compiler `to enforce standard conformance`_. When
    building with Visual Studio 2019, this is not strictly necessary, but still advised.

..  _`to enforce standard conformance`: https://docs.microsoft.com/en-us/cpp/build/reference/permissive-standards-conformance?view=vs-2017

To compile and run the tests:

.. code-block:: batch

   mkdir build
   cd build
   cmake ..
   cmake --build . --config Release --target check

This will create a Visual Studio project, compile and run the target, all from the
command line.

.. Note::

    If all tests fail, make sure that the Python binary and the testcases are compiled
    for the same processor type and bitness (i.e. either **i386** or **x86_64**). You
    can specify **x86_64** as the target architecture for the generated Visual Studio
    project using ``cmake -A x64 ..``.

.. seealso::

    Advanced users who are already familiar with Boost.Python may want to skip
    the tutorial and look at the test cases in the :file:`tests` directory,
    which exercise all features of pybind11.

Header and namespace conventions
================================

For brevity, all code examples assume that the following two lines are present:

.. code-block:: cpp

    #include <pybind11/pybind11.h>

    namespace py = pybind11;

Some features may require additional headers, but those will be specified as needed.

.. _simple_example:

Creating bindings for a simple function
=======================================

Let's start by creating Python bindings for an extremely simple function, which
adds two numbers and returns their result:

.. code-block:: cpp

    int add(int i, int j) {
        return i + j;
    }

For simplicity [#f1]_, we'll put both this function and the binding code into
a file named :file:`example.cpp` with the following contents:

.. code-block:: cpp

    #include <pybind11/pybind11.h>

    int add(int i, int j) {
        return i + j;
    }

    PYBIND11_MODULE(example, m) {
        m.doc() = "pybind11 example plugin"; // optional module docstring

        m.def("add", &add, "A function which adds two numbers");
    }

.. [#f1] In practice, implementation and binding code will generally be located
         in separate files.

The :func:`PYBIND11_MODULE` macro creates a function that will be called when an
``import`` statement is issued from within Python. The module name (``example``)
is given as the first macro argument (it should not be in quotes). The second
argument (``m``) defines a variable of type :class:`py::module_ <module>` which
is the main interface for creating bindings. The method :func:`module_::def`
generates binding code that exposes the ``add()`` function to Python.

.. note::

    Notice how little code was needed to expose our function to Python: all
    details regarding the function's parameters and return value were
    automatically inferred using template metaprogramming. This overall
    approach and the used syntax are borrowed from Boost.Python, though the
    underlying implementation is very different.

pybind11 is a header-only library, hence it is not necessary to link against
any special libraries and there are no intermediate (magic) translation steps.
On Linux, the above example can be compiled using the following command:

.. code-block:: bash

    $ c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) example.cpp -o example$(python3-config --extension-suffix)

.. note::

    If you used :ref:`include_as_a_submodule` to get the pybind11 source, then
    use ``$(python3-config --includes) -Iextern/pybind11/include`` instead of
    ``$(python3 -m pybind11 --includes)`` in the above compilation, as
    explained in :ref:`building_manually`.

For more details on the required compiler flags on Linux and macOS, see
:ref:`building_manually`. For complete cross-platform compilation instructions,
refer to the :ref:`compiling` page.

The `python_example`_ and `cmake_example`_ repositories are also a good place
to start. They are both complete project examples with cross-platform build
systems. The only difference between the two is that `python_example`_ uses
Python's ``setuptools`` to build the module, while `cmake_example`_ uses CMake
(which may be preferable for existing C++ projects).

.. _python_example: https://github.com/pybind/python_example
.. _cmake_example: https://github.com/pybind/cmake_example

Building the above C++ code will produce a binary module file that can be
imported to Python. Assuming that the compiled module is located in the
current directory, the following interactive Python session shows how to
load and execute the example:

.. code-block:: pycon

    $ python
    Python 2.7.10 (default, Aug 22 2015, 20:33:39)
    [GCC 4.2.1 Compatible Apple LLVM 7.0.0 (clang-700.0.59.1)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import example
    >>> example.add(1, 2)
    3L
    >>>

.. _keyword_args:

Keyword arguments
=================

With a simple code modification, it is possible to inform Python about the
names of the arguments ("i" and "j" in this case).

.. code-block:: cpp

    m.def("add", &add, "A function which adds two numbers",
          py::arg("i"), py::arg("j"));

:class:`arg` is one of several special tag classes which can be used to pass
metadata into :func:`module_::def`. With this modified binding code, we can now
call the function using keyword arguments, which is a more readable alternative
particularly for functions taking many parameters:

.. code-block:: pycon

    >>> import example
    >>> example.add(i=1, j=2)
    3L

The keyword names also appear in the function signatures within the documentation.

.. code-block:: pycon

    >>> help(example)

    ....

    FUNCTIONS
        add(...)
            Signature : (i: int, j: int) -> int

            A function which adds two numbers

A shorter notation for named arguments is also available:

.. code-block:: cpp

    // regular notation
    m.def("add1", &add, py::arg("i"), py::arg("j"));
    // shorthand
    using namespace pybind11::literals;
    m.def("add2", &add, "i"_a, "j"_a);

The :var:`_a` suffix forms a C++11 literal which is equivalent to :class:`arg`.
Note that the literal operator must first be made visible with the directive
``using namespace pybind11::literals``. This does not bring in anything else
from the ``pybind11`` namespace except for literals.

.. _default_args:

Default arguments
=================

Suppose now that the function to be bound has default arguments, e.g.:

.. code-block:: cpp

    int add(int i = 1, int j = 2) {
        return i + j;
    }

Unfortunately, pybind11 cannot automatically extract these parameters, since they
are not part of the function's type information. However, they are simple to specify
using an extension of :class:`arg`:

.. code-block:: cpp

    m.def("add", &add, "A function which adds two numbers",
          py::arg("i") = 1, py::arg("j") = 2);

The default values also appear within the documentation.

.. code-block:: pycon

    >>> help(example)

    ....

    FUNCTIONS
        add(...)
            Signature : (i: int = 1, j: int = 2) -> int

            A function which adds two numbers

The shorthand notation is also available for default arguments:

.. code-block:: cpp

    // regular notation
    m.def("add1", &add, py::arg("i") = 1, py::arg("j") = 2);
    // shorthand
    m.def("add2", &add, "i"_a=1, "j"_a=2);

Exporting variables
===================

To expose a value from C++, use the ``attr`` function to register it in a
module as shown below. Built-in types and general objects (more on that later)
are automatically converted when assigned as attributes, and can be explicitly
converted using the function ``py::cast``.

.. code-block:: cpp

    PYBIND11_MODULE(example, m) {
        m.attr("the_answer") = 42;
        py::object world = py::cast("World");
        m.attr("what") = world;
    }

These are then accessible from Python:

.. code-block:: pycon

    >>> import example
    >>> example.the_answer
    42
    >>> example.what
    'World'

.. _supported_types:

Supported data types
====================

A large number of data types are supported out of the box and can be used
seamlessly as functions arguments, return values or with ``py::cast`` in general.
For a full overview, see the :doc:`advanced/cast/index` section.
Benchmark
=========

The following is the result of a synthetic benchmark comparing both compilation
time and module size of pybind11 against Boost.Python. A detailed report about a
Boost.Python to pybind11 conversion of a real project is available here: [#f1]_.

.. [#f1] http://graylab.jhu.edu/RosettaCon2016/PyRosetta-4.pdf

Setup
-----

A python script (see the ``docs/benchmark.py`` file) was used to generate a set
of files with dummy classes whose count increases for each successive benchmark
(between 1 and 2048 classes in powers of two). Each class has four methods with
a randomly generated signature with a return value and four arguments. (There
was no particular reason for this setup other than the desire to generate many
unique function signatures whose count could be controlled in a simple way.)

Here is an example of the binding code for one class:

.. code-block:: cpp

    ...
    class cl034 {
    public:
        cl279 *fn_000(cl084 *, cl057 *, cl065 *, cl042 *);
        cl025 *fn_001(cl098 *, cl262 *, cl414 *, cl121 *);
        cl085 *fn_002(cl445 *, cl297 *, cl145 *, cl421 *);
        cl470 *fn_003(cl200 *, cl323 *, cl332 *, cl492 *);
    };
    ...

    PYBIND11_MODULE(example, m) {
        ...
        py::class_<cl034>(m, "cl034")
            .def("fn_000", &cl034::fn_000)
            .def("fn_001", &cl034::fn_001)
            .def("fn_002", &cl034::fn_002)
            .def("fn_003", &cl034::fn_003)
        ...
    }

The Boost.Python version looks almost identical except that a return value
policy had to be specified as an argument to ``def()``. For both libraries,
compilation was done with

.. code-block:: bash

    Apple LLVM version 7.0.2 (clang-700.1.81)

and the following compilation flags

.. code-block:: bash

    g++ -Os -shared -rdynamic -undefined dynamic_lookup -fvisibility=hidden -std=c++14

Compilation time
----------------

The following log-log plot shows how the compilation time grows for an
increasing number of class and function declarations. pybind11 includes many
fewer headers, which initially leads to shorter compilation times, but the
performance is ultimately fairly similar (pybind11 is 19.8 seconds faster for
the largest largest file with 2048 classes and a total of 8192 methods -- a
modest **1.2x** speedup relative to Boost.Python, which required 116.35
seconds).

.. only:: not latex

    .. image:: pybind11_vs_boost_python1.svg

.. only:: latex

    .. image:: pybind11_vs_boost_python1.png

Module size
-----------

Differences between the two libraries become much more pronounced when
considering the file size of the generated Python plugin: for the largest file,
the binary generated by Boost.Python required 16.8 MiB, which was **2.17
times** / **9.1 megabytes** larger than the output generated by pybind11. For
very small inputs, Boost.Python has an edge in the plot below -- however, note
that it stores many definitions in an external library, whose size was not
included here, hence the comparison is slightly shifted in Boost.Python's
favor.

.. only:: not latex

    .. image:: pybind11_vs_boost_python2.svg

.. only:: latex

    .. image:: pybind11_vs_boost_python2.png
Frequently asked questions
##########################

"ImportError: dynamic module does not define init function"
===========================================================

1. Make sure that the name specified in PYBIND11_MODULE is identical to the
filename of the extension library (without suffixes such as .so)

2. If the above did not fix the issue, you are likely using an incompatible
version of Python (for instance, the extension library was compiled against
Python 2, while the interpreter is running on top of some version of Python
3, or vice versa).

"Symbol not found: ``__Py_ZeroStruct`` / ``_PyInstanceMethod_Type``"
========================================================================

See the first answer.

"SystemError: dynamic module not initialized properly"
======================================================

See the first answer.

The Python interpreter immediately crashes when importing my module
===================================================================

See the first answer.

.. _faq_reference_arguments:

Limitations involving reference arguments
=========================================

In C++, it's fairly common to pass arguments using mutable references or
mutable pointers, which allows both read and write access to the value
supplied by the caller. This is sometimes done for efficiency reasons, or to
realize functions that have multiple return values. Here are two very basic
examples:

.. code-block:: cpp

    void increment(int &i) { i++; }
    void increment_ptr(int *i) { (*i)++; }

In Python, all arguments are passed by reference, so there is no general
issue in binding such code from Python.

However, certain basic Python types (like ``str``, ``int``, ``bool``,
``float``, etc.) are **immutable**. This means that the following attempt
to port the function to Python doesn't have the same effect on the value
provided by the caller -- in fact, it does nothing at all.

.. code-block:: python

    def increment(i):
        i += 1 # nope..

pybind11 is also affected by such language-level conventions, which means that
binding ``increment`` or ``increment_ptr`` will also create Python functions
that don't modify their arguments.

Although inconvenient, one workaround is to encapsulate the immutable types in
a custom type that does allow modifications.

An other alternative involves binding a small wrapper lambda function that
returns a tuple with all output arguments (see the remainder of the
documentation for examples on binding lambda functions). An example:

.. code-block:: cpp

    int foo(int &i) { i++; return 123; }

and the binding code

.. code-block:: cpp

   m.def("foo", [](int i) { int rv = foo(i); return std::make_tuple(rv, i); });


How can I reduce the build time?
================================

It's good practice to split binding code over multiple files, as in the
following example:

:file:`example.cpp`:

.. code-block:: cpp

    void init_ex1(py::module_ &);
    void init_ex2(py::module_ &);
    /* ... */

    PYBIND11_MODULE(example, m) {
        init_ex1(m);
        init_ex2(m);
        /* ... */
    }

:file:`ex1.cpp`:

.. code-block:: cpp

    void init_ex1(py::module_ &m) {
        m.def("add", [](int a, int b) { return a + b; });
    }

:file:`ex2.cpp`:

.. code-block:: cpp

    void init_ex2(py::module_ &m) {
        m.def("sub", [](int a, int b) { return a - b; });
    }

:command:`python`:

.. code-block:: pycon

    >>> import example
    >>> example.add(1, 2)
    3
    >>> example.sub(1, 1)
    0

As shown above, the various ``init_ex`` functions should be contained in
separate files that can be compiled independently from one another, and then
linked together into the same final shared object.  Following this approach
will:

1. reduce memory requirements per compilation unit.

2. enable parallel builds (if desired).

3. allow for faster incremental builds. For instance, when a single class
   definition is changed, only a subset of the binding code will generally need
   to be recompiled.

"recursive template instantiation exceeded maximum depth of 256"
================================================================

If you receive an error about excessive recursive template evaluation, try
specifying a larger value, e.g. ``-ftemplate-depth=1024`` on GCC/Clang. The
culprit is generally the generation of function signatures at compile time
using C++14 template metaprogramming.

.. _`faq:hidden_visibility`:

"‘SomeClass’ declared with greater visibility than the type of its field ‘SomeClass::member’ [-Wattributes]"
============================================================================================================

This error typically indicates that you are compiling without the required
``-fvisibility`` flag.  pybind11 code internally forces hidden visibility on
all internal code, but if non-hidden (and thus *exported*) code attempts to
include a pybind type (for example, ``py::object`` or ``py::list``) you can run
into this warning.

To avoid it, make sure you are specifying ``-fvisibility=hidden`` when
compiling pybind code.

As to why ``-fvisibility=hidden`` is necessary, because pybind modules could
have been compiled under different versions of pybind itself, it is also
important that the symbols defined in one module do not clash with the
potentially-incompatible symbols defined in another.  While Python extension
modules are usually loaded with localized symbols (under POSIX systems
typically using ``dlopen`` with the ``RTLD_LOCAL`` flag), this Python default
can be changed, but even if it isn't it is not always enough to guarantee
complete independence of the symbols involved when not using
``-fvisibility=hidden``.

Additionally, ``-fvisiblity=hidden`` can deliver considerably binary size
savings.  (See the following section for more details).


.. _`faq:symhidden`:

How can I create smaller binaries?
==================================

To do its job, pybind11 extensively relies on a programming technique known as
*template metaprogramming*, which is a way of performing computation at compile
time using type information. Template metaprogamming usually instantiates code
involving significant numbers of deeply nested types that are either completely
removed or reduced to just a few instructions during the compiler's optimization
phase. However, due to the nested nature of these types, the resulting symbol
names in the compiled extension library can be extremely long. For instance,
the included test suite contains the following symbol:

.. only:: html

    .. code-block:: none

        _​_​Z​N​8​p​y​b​i​n​d​1​1​1​2​c​p​p​_​f​u​n​c​t​i​o​n​C​1​I​v​8​E​x​a​m​p​l​e​2​J​R​N​S​t​3​_​_​1​6​v​e​c​t​o​r​I​N​S​3​_​1​2​b​a​s​i​c​_​s​t​r​i​n​g​I​w​N​S​3​_​1​1​c​h​a​r​_​t​r​a​i​t​s​I​w​E​E​N​S​3​_​9​a​l​l​o​c​a​t​o​r​I​w​E​E​E​E​N​S​8​_​I​S​A​_​E​E​E​E​E​J​N​S​_​4​n​a​m​e​E​N​S​_​7​s​i​b​l​i​n​g​E​N​S​_​9​i​s​_​m​e​t​h​o​d​E​A​2​8​_​c​E​E​E​M​T​0​_​F​T​_​D​p​T​1​_​E​D​p​R​K​T​2​_

.. only:: not html

    .. code-block:: cpp

        __ZN8pybind1112cpp_functionC1Iv8Example2JRNSt3__16vectorINS3_12basic_stringIwNS3_11char_traitsIwEENS3_9allocatorIwEEEENS8_ISA_EEEEEJNS_4nameENS_7siblingENS_9is_methodEA28_cEEEMT0_FT_DpT1_EDpRKT2_

which is the mangled form of the following function type:

.. code-block:: cpp

    pybind11::cpp_function::cpp_function<void, Example2, std::__1::vector<std::__1::basic_string<wchar_t, std::__1::char_traits<wchar_t>, std::__1::allocator<wchar_t> >, std::__1::allocator<std::__1::basic_string<wchar_t, std::__1::char_traits<wchar_t>, std::__1::allocator<wchar_t> > > >&, pybind11::name, pybind11::sibling, pybind11::is_method, char [28]>(void (Example2::*)(std::__1::vector<std::__1::basic_string<wchar_t, std::__1::char_traits<wchar_t>, std::__1::allocator<wchar_t> >, std::__1::allocator<std::__1::basic_string<wchar_t, std::__1::char_traits<wchar_t>, std::__1::allocator<wchar_t> > > >&), pybind11::name const&, pybind11::sibling const&, pybind11::is_method const&, char const (&) [28])

The memory needed to store just the mangled name of this function (196 bytes)
is larger than the actual piece of code (111 bytes) it represents! On the other
hand, it's silly to even give this function a name -- after all, it's just a
tiny cog in a bigger piece of machinery that is not exposed to the outside
world. So we'll generally only want to export symbols for those functions which
are actually called from the outside.

This can be achieved by specifying the parameter ``-fvisibility=hidden`` to GCC
and Clang, which sets the default symbol visibility to *hidden*, which has a
tremendous impact on the final binary size of the resulting extension library.
(On Visual Studio, symbols are already hidden by default, so nothing needs to
be done there.)

In addition to decreasing binary size, ``-fvisibility=hidden`` also avoids
potential serious issues when loading multiple modules and is required for
proper pybind operation.  See the previous FAQ entry for more details.

Working with ancient Visual Studio 2008 builds on Windows
=========================================================

The official Windows distributions of Python are compiled using truly
ancient versions of Visual Studio that lack good C++11 support. Some users
implicitly assume that it would be impossible to load a plugin built with
Visual Studio 2015 into a Python distribution that was compiled using Visual
Studio 2008. However, no such issue exists: it's perfectly legitimate to
interface DLLs that are built with different compilers and/or C libraries.
Common gotchas to watch out for involve not ``free()``-ing memory region
that that were ``malloc()``-ed in another shared library, using data
structures with incompatible ABIs, and so on. pybind11 is very careful not
to make these types of mistakes.

How can I properly handle Ctrl-C in long-running functions?
===========================================================

Ctrl-C is received by the Python interpreter, and holds it until the GIL
is released, so a long-running function won't be interrupted.

To interrupt from inside your function, you can use the ``PyErr_CheckSignals()``
function, that will tell if a signal has been raised on the Python side.  This
function merely checks a flag, so its impact is negligible. When a signal has
been received, you must either explicitly interrupt execution by throwing
``py::error_already_set`` (which will propagate the existing
``KeyboardInterrupt``), or clear the error (which you usually will not want):

.. code-block:: cpp

    PYBIND11_MODULE(example, m)
    {
        m.def("long running_func", []()
        {
            for (;;) {
                if (PyErr_CheckSignals() != 0)
                    throw py::error_already_set();
                // Long running iteration
            }
        });
    }

CMake doesn't detect the right Python version
=============================================

The CMake-based build system will try to automatically detect the installed
version of Python and link against that. When this fails, or when there are
multiple versions of Python and it finds the wrong one, delete
``CMakeCache.txt`` and then add ``-DPYTHON_EXECUTABLE=$(which python)`` to your
CMake configure line. (Replace ``$(which python)`` with a path to python if
your prefer.)

You can alternatively try ``-DPYBIND11_FINDPYTHON=ON``, which will activate the
new CMake FindPython support instead of pybind11's custom search. Requires
CMake 3.12+, and 3.15+ or 3.18.2+ are even better. You can set this in your
``CMakeLists.txt`` before adding or finding pybind11, as well.

Inconsistent detection of Python version in CMake and pybind11
==============================================================

The functions ``find_package(PythonInterp)`` and ``find_package(PythonLibs)``
provided by CMake for Python version detection are modified by pybind11 due to
unreliability and limitations that make them unsuitable for pybind11's needs.
Instead pybind11 provides its own, more reliable Python detection CMake code.
Conflicts can arise, however, when using pybind11 in a project that *also* uses
the CMake Python detection in a system with several Python versions installed.

This difference may cause inconsistencies and errors if *both* mechanisms are
used in the same project. Consider the following CMake code executed in a
system with Python 2.7 and 3.x installed:

.. code-block:: cmake

    find_package(PythonInterp)
    find_package(PythonLibs)
    find_package(pybind11)

It will detect Python 2.7 and pybind11 will pick it as well.

In contrast this code:

.. code-block:: cmake

    find_package(pybind11)
    find_package(PythonInterp)
    find_package(PythonLibs)

will detect Python 3.x for pybind11 and may crash on
``find_package(PythonLibs)`` afterwards.

There are three possible solutions:

1. Avoid using ``find_package(PythonInterp)`` and ``find_package(PythonLibs)``
   from CMake and rely on pybind11 in detecting Python version. If this is not
   possible, the CMake machinery should be called *before* including pybind11.
2. Set ``PYBIND11_FINDPYTHON`` to ``True`` or use ``find_package(Python
   COMPONENTS Interpreter Development)`` on modern CMake (3.12+, 3.15+ better,
   3.18.2+ best). Pybind11 in these cases uses the new CMake FindPython instead
   of the old, deprecated search tools, and these modules are much better at
   finding the correct Python.
3. Set ``PYBIND11_NOPYTHON`` to ``TRUE``. Pybind11 will not search for Python.
   However, you will have to use the target-based system, and do more setup
   yourself, because it does not know about or include things that depend on
   Python, like ``pybind11_add_module``. This might be ideal for integrating
   into an existing system, like scikit-build's Python helpers.

How to cite this project?
=========================

We suggest the following BibTeX template to cite pybind11 in scientific
discourse:

.. code-block:: bash

    @misc{pybind11,
       author = {Wenzel Jakob and Jason Rhinelander and Dean Moldovan},
       year = {2017},
       note = {https://github.com/pybind/pybind11},
       title = {pybind11 -- Seamless operability between C++11 and Python}
    }
.. only:: latex

   Intro
   =====

.. include:: readme.rst

.. only:: not latex

    Contents:

.. toctree::
   :maxdepth: 1

   changelog
   upgrade

.. toctree::
   :caption: The Basics
   :maxdepth: 2

   installing
   basics
   classes
   compiling

.. toctree::
   :caption: Advanced Topics
   :maxdepth: 2

   advanced/functions
   advanced/classes
   advanced/exceptions
   advanced/smart_ptrs
   advanced/cast/index
   advanced/pycpp/index
   advanced/embedding
   advanced/misc

.. toctree::
   :caption: Extra Information
   :maxdepth: 1

   faq
   benchmark
   limitations
   reference
   cmake/index
.. _compiling:

Build systems
#############

.. _build-setuptools:

Building with setuptools
========================

For projects on PyPI, building with setuptools is the way to go. Sylvain Corlay
has kindly provided an example project which shows how to set up everything,
including automatic generation of documentation using Sphinx. Please refer to
the [python_example]_ repository.

.. [python_example] https://github.com/pybind/python_example

A helper file is provided with pybind11 that can simplify usage with setuptools.

To use pybind11 inside your ``setup.py``, you have to have some system to
ensure that ``pybind11`` is installed when you build your package. There are
four possible ways to do this, and pybind11 supports all four: You can ask all
users to install pybind11 beforehand (bad), you can use
:ref:`setup_helpers-pep518` (good, but very new and requires Pip 10),
:ref:`setup_helpers-setup_requires` (discouraged by Python packagers now that
PEP 518 is available, but it still works everywhere), or you can
:ref:`setup_helpers-copy-manually` (always works but you have to manually sync
your copy to get updates).

An example of a ``setup.py`` using pybind11's helpers:

.. code-block:: python

    from glob import glob
    from setuptools import setup
    from pybind11.setup_helpers import Pybind11Extension

    ext_modules = [
        Pybind11Extension(
            "python_example",
            sorted(glob("src/*.cpp")),  # Sort source files for reproducibility
        ),
    ]

    setup(
        ...,
        ext_modules=ext_modules
    )

If you want to do an automatic search for the highest supported C++ standard,
that is supported via a ``build_ext`` command override; it will only affect
``Pybind11Extensions``:

.. code-block:: python

    from glob import glob
    from setuptools import setup
    from pybind11.setup_helpers import Pybind11Extension, build_ext

    ext_modules = [
        Pybind11Extension(
            "python_example",
            sorted(glob("src/*.cpp")),
        ),
    ]

    setup(
        ...,
        cmdclass={"build_ext": build_ext},
        ext_modules=ext_modules
    )

Since pybind11 does not require NumPy when building, a light-weight replacement
for NumPy's parallel compilation distutils tool is included. Use it like this:

.. code-block:: python

    from pybind11.setup_helpers import ParallelCompile

    # Optional multithreaded build
    ParallelCompile("NPY_NUM_BUILD_JOBS").install()

    setup(...)

The argument is the name of an environment variable to control the number of
threads, such as ``NPY_NUM_BUILD_JOBS`` (as used by NumPy), though you can set
something different if you want; ``CMAKE_BUILD_PARALLEL_LEVEL`` is another choice
a user might expect. You can also pass ``default=N`` to set the default number
of threads (0 will take the number of threads available) and ``max=N``, the
maximum number of threads; if you have a large extension you may want set this
to a memory dependent number.

If you are developing rapidly and have a lot of C++ files, you may want to
avoid rebuilding files that have not changed. For simple cases were you are
using ``pip install -e .`` and do not have local headers, you can skip the
rebuild if a object file is newer than it's source (headers are not checked!)
with the following:

.. code-block:: python

    from pybind11.setup_helpers import ParallelCompile, naive_recompile

    SmartCompile("NPY_NUM_BUILD_JOBS", needs_recompile=naive_recompile).install()


If you have a more complex build, you can implement a smarter function and pass
it to ``needs_recompile``, or you can use [Ccache]_ instead. ``CXX="cache g++"
pip install -e .`` would be the way to use it with GCC, for example. Unlike the
simple solution, this even works even when not compiling in editable mode, but
it does require Ccache to be installed.

Keep in mind that Pip will not even attempt to rebuild if it thinks it has
already built a copy of your code, which it deduces from the version number.
One way to avoid this is to use [setuptools_scm]_, which will generate a
version number that includes the number of commits since your last tag and a
hash for a dirty directory. Another way to force a rebuild is purge your cache
or use Pip's ``--no-cache-dir`` option.

.. [Ccache] https://ccache.dev

.. [setuptools_scm] https://github.com/pypa/setuptools_scm

.. _setup_helpers-pep518:

PEP 518 requirements (Pip 10+ required)
---------------------------------------

If you use `PEP 518's <https://www.python.org/dev/peps/pep-0518/>`_
``pyproject.toml`` file, you can ensure that ``pybind11`` is available during
the compilation of your project.  When this file exists, Pip will make a new
virtual environment, download just the packages listed here in ``requires=``,
and build a wheel (binary Python package). It will then throw away the
environment, and install your wheel.

Your ``pyproject.toml`` file will likely look something like this:

.. code-block:: toml

    [build-system]
    requires = ["setuptools>=42", "wheel", "pybind11~=2.6.1"]
    build-backend = "setuptools.build_meta"

.. note::

    The main drawback to this method is that a `PEP 517`_ compliant build tool,
    such as Pip 10+, is required for this approach to work; older versions of
    Pip completely ignore this file. If you distribute binaries (called wheels
    in Python) using something like `cibuildwheel`_, remember that ``setup.py``
    and ``pyproject.toml`` are not even contained in the wheel, so this high
    Pip requirement is only for source builds, and will not affect users of
    your binary wheels. If you are building SDists and wheels, then
    `pypa-build`_ is the recommended offical tool.

.. _PEP 517: https://www.python.org/dev/peps/pep-0517/
.. _cibuildwheel: https://cibuildwheel.readthedocs.io
.. _pypa-build: https://pypa-build.readthedocs.io/en/latest/

.. _setup_helpers-setup_requires:

Classic ``setup_requires``
--------------------------

If you want to support old versions of Pip with the classic
``setup_requires=["pybind11"]`` keyword argument to setup, which triggers a
two-phase ``setup.py`` run, then you will need to use something like this to
ensure the first pass works (which has not yet installed the ``setup_requires``
packages, since it can't install something it does not know about):

.. code-block:: python

    try:
        from pybind11.setup_helpers import Pybind11Extension
    except ImportError:
        from setuptools import Extension as Pybind11Extension


It doesn't matter that the Extension class is not the enhanced subclass for the
first pass run; and the second pass will have the ``setup_requires``
requirements.

This is obviously more of a hack than the PEP 518 method, but it supports
ancient versions of Pip.

.. _setup_helpers-copy-manually:

Copy manually
-------------

You can also copy ``setup_helpers.py`` directly to your project; it was
designed to be usable standalone, like the old example ``setup.py``. You can
set ``include_pybind11=False`` to skip including the pybind11 package headers,
so you can use it with git submodules and a specific git version. If you use
this, you will need to import from a local file in ``setup.py`` and ensure the
helper file is part of your MANIFEST.


Closely related, if you include pybind11 as a subproject, you can run the
``setup_helpers.py`` inplace. If loaded correctly, this should even pick up
the correct include for pybind11, though you can turn it off as shown above if
you want to input it manually.

Suggested usage if you have pybind11 as a submodule in ``extern/pybind11``:

.. code-block:: python

    DIR = os.path.abspath(os.path.dirname(__file__))

    sys.path.append(os.path.join(DIR, "extern", "pybind11"))
    from pybind11.setup_helpers import Pybind11Extension  # noqa: E402

    del sys.path[-1]


.. versionchanged:: 2.6

    Added ``setup_helpers`` file.

Building with cppimport
========================

[cppimport]_ is a small Python import hook that determines whether there is a C++
source file whose name matches the requested module. If there is, the file is
compiled as a Python extension using pybind11 and placed in the same folder as
the C++ source file. Python is then able to find the module and load it.

.. [cppimport] https://github.com/tbenthompson/cppimport

.. _cmake:

Building with CMake
===================

For C++ codebases that have an existing CMake-based build system, a Python
extension module can be created with just a few lines of code:

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.4...3.18)
    project(example LANGUAGES CXX)

    add_subdirectory(pybind11)
    pybind11_add_module(example example.cpp)

This assumes that the pybind11 repository is located in a subdirectory named
:file:`pybind11` and that the code is located in a file named :file:`example.cpp`.
The CMake command ``add_subdirectory`` will import the pybind11 project which
provides the ``pybind11_add_module`` function. It will take care of all the
details needed to build a Python extension module on any platform.

A working sample project, including a way to invoke CMake from :file:`setup.py` for
PyPI integration, can be found in the [cmake_example]_  repository.

.. [cmake_example] https://github.com/pybind/cmake_example

.. versionchanged:: 2.6
   CMake 3.4+ is required.

Further information can be found at :doc:`cmake/index`.

pybind11_add_module
-------------------

To ease the creation of Python extension modules, pybind11 provides a CMake
function with the following signature:

.. code-block:: cmake

    pybind11_add_module(<name> [MODULE | SHARED] [EXCLUDE_FROM_ALL]
                        [NO_EXTRAS] [THIN_LTO] [OPT_SIZE] source1 [source2 ...])

This function behaves very much like CMake's builtin ``add_library`` (in fact,
it's a wrapper function around that command). It will add a library target
called ``<name>`` to be built from the listed source files. In addition, it
will take care of all the Python-specific compiler and linker flags as well
as the OS- and Python-version-specific file extension. The produced target
``<name>`` can be further manipulated with regular CMake commands.

``MODULE`` or ``SHARED`` may be given to specify the type of library. If no
type is given, ``MODULE`` is used by default which ensures the creation of a
Python-exclusive module. Specifying ``SHARED`` will create a more traditional
dynamic library which can also be linked from elsewhere. ``EXCLUDE_FROM_ALL``
removes this target from the default build (see CMake docs for details).

Since pybind11 is a template library, ``pybind11_add_module`` adds compiler
flags to ensure high quality code generation without bloat arising from long
symbol names and duplication of code in different translation units. It
sets default visibility to *hidden*, which is required for some pybind11
features and functionality when attempting to load multiple pybind11 modules
compiled under different pybind11 versions.  It also adds additional flags
enabling LTO (Link Time Optimization) and strip unneeded symbols. See the
:ref:`FAQ entry <faq:symhidden>` for a more detailed explanation. These
latter optimizations are never applied in ``Debug`` mode.  If ``NO_EXTRAS`` is
given, they will always be disabled, even in ``Release`` mode. However, this
will result in code bloat and is generally not recommended.

As stated above, LTO is enabled by default. Some newer compilers also support
different flavors of LTO such as `ThinLTO`_. Setting ``THIN_LTO`` will cause
the function to prefer this flavor if available. The function falls back to
regular LTO if ``-flto=thin`` is not available. If
``CMAKE_INTERPROCEDURAL_OPTIMIZATION`` is set (either ``ON`` or ``OFF``), then
that will be respected instead of the built-in flag search.

.. note::

   If you want to set the property form on targets or the
   ``CMAKE_INTERPROCEDURAL_OPTIMIZATION_<CONFIG>`` versions of this, you should
   still use ``set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)`` (otherwise a
   no-op) to disable pybind11's ipo flags.

The ``OPT_SIZE`` flag enables size-based optimization equivalent to the
standard ``/Os`` or ``-Os`` compiler flags and the ``MinSizeRel`` build type,
which avoid optimizations that that can substantially increase the size of the
resulting binary. This flag is particularly useful in projects that are split
into performance-critical parts and associated bindings. In this case, we can
compile the project in release mode (and hence, optimize performance globally),
and specify ``OPT_SIZE`` for the binding target, where size might be the main
concern as performance is often less critical here. A ~25% size reduction has
been observed in practice. This flag only changes the optimization behavior at
a per-target level and takes precedence over the global CMake build type
(``Release``, ``RelWithDebInfo``) except for ``Debug`` builds, where
optimizations remain disabled.

.. _ThinLTO: http://clang.llvm.org/docs/ThinLTO.html

Configuration variables
-----------------------

By default, pybind11 will compile modules with the compiler default or the
minimum standard required by pybind11, whichever is higher.  You can set the
standard explicitly with
`CMAKE_CXX_STANDARD <https://cmake.org/cmake/help/latest/variable/CMAKE_CXX_STANDARD.html>`_:

.. code-block:: cmake

    set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ version selection")  # or 11, 14, 17, 20
    set(CMAKE_CXX_STANDARD_REQUIRED ON)  # optional, ensure standard is supported
    set(CMAKE_CXX_EXTENSIONS OFF)  # optional, keep compiler extensionsn off

The variables can also be set when calling CMake from the command line using
the ``-D<variable>=<value>`` flag. You can also manually set ``CXX_STANDARD``
on a target or use ``target_compile_features`` on your targets - anything that
CMake supports.

Classic Python support: The target Python version can be selected by setting
``PYBIND11_PYTHON_VERSION`` or an exact Python installation can be specified
with ``PYTHON_EXECUTABLE``.  For example:

.. code-block:: bash

    cmake -DPYBIND11_PYTHON_VERSION=3.6 ..

    # Another method:
    cmake -DPYTHON_EXECUTABLE=/path/to/python ..

    # This often is a good way to get the current Python, works in environments:
    cmake -DPYTHON_EXECUTABLE=$(python3 -c "import sys; print(sys.executable)") ..


find_package vs. add_subdirectory
---------------------------------

For CMake-based projects that don't include the pybind11 repository internally,
an external installation can be detected through ``find_package(pybind11)``.
See the `Config file`_ docstring for details of relevant CMake variables.

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.4...3.18)
    project(example LANGUAGES CXX)

    find_package(pybind11 REQUIRED)
    pybind11_add_module(example example.cpp)

Note that ``find_package(pybind11)`` will only work correctly if pybind11
has been correctly installed on the system, e. g. after downloading or cloning
the pybind11 repository  :

.. code-block:: bash

    # Classic CMake
    cd pybind11
    mkdir build
    cd build
    cmake ..
    make install

    # CMake 3.15+
    cd pybind11
    cmake -S . -B build
    cmake --build build -j 2  # Build on 2 cores
    cmake --install build

Once detected, the aforementioned ``pybind11_add_module`` can be employed as
before. The function usage and configuration variables are identical no matter
if pybind11 is added as a subdirectory or found as an installed package. You
can refer to the same [cmake_example]_ repository for a full sample project
-- just swap out ``add_subdirectory`` for ``find_package``.

.. _Config file: https://github.com/pybind/pybind11/blob/master/tools/pybind11Config.cmake.in


.. _find-python-mode:

FindPython mode
---------------

CMake 3.12+ (3.15+ recommended, 3.18.2+ ideal) added a new module called
FindPython that had a highly improved search algorithm and modern targets
and tools. If you use FindPython, pybind11 will detect this and use the
existing targets instead:

.. code-block:: cmake

    cmake_minumum_required(VERSION 3.15...3.19)
    project(example LANGUAGES CXX)

    find_package(Python COMPONENTS Interpreter Development REQUIRED)
    find_package(pybind11 CONFIG REQUIRED)
    # or add_subdirectory(pybind11)

    pybind11_add_module(example example.cpp)

You can also use the targets (as listed below) with FindPython. If you define
``PYBIND11_FINDPYTHON``, pybind11 will perform the FindPython step for you
(mostly useful when building pybind11's own tests, or as a way to change search
algorithms from the CMake invocation, with ``-DPYBIND11_FINDPYTHON=ON``.

.. warning::

    If you use FindPython2 and FindPython3 to dual-target Python, use the
    individual targets listed below, and avoid targets that directly include
    Python parts.

There are `many ways to hint or force a discovery of a specific Python
installation <https://cmake.org/cmake/help/latest/module/FindPython.html>`_),
setting ``Python_ROOT_DIR`` may be the most common one (though with
virtualenv/venv support, and Conda support, this tends to find the correct
Python version more often than the old system did).

.. warning::

    When the Python libraries (i.e. ``libpythonXX.a`` and ``libpythonXX.so``
    on Unix) are not available, as is the case on a manylinux image, the
    ``Development`` component will not be resolved by ``FindPython``. When not
    using the embedding functionality, CMake 3.18+ allows you to specify
    ``Development.Module`` instead of ``Development`` to resolve this issue.

.. versionadded:: 2.6

Advanced: interface library targets
-----------------------------------

Pybind11 supports modern CMake usage patterns with a set of interface targets,
available in all modes. The targets provided are:

   ``pybind11::headers``
     Just the pybind11 headers and minimum compile requirements

   ``pybind11::python2_no_register``
     Quiets the warning/error when mixing C++14 or higher and Python 2

   ``pybind11::pybind11``
     Python headers + ``pybind11::headers`` + ``pybind11::python2_no_register`` (Python 2 only)

   ``pybind11::python_link_helper``
     Just the "linking" part of pybind11:module

   ``pybind11::module``
     Everything for extension modules - ``pybind11::pybind11`` + ``Python::Module`` (FindPython CMake 3.15+) or ``pybind11::python_link_helper``

   ``pybind11::embed``
     Everything for embedding the Python interpreter - ``pybind11::pybind11`` + ``Python::Embed`` (FindPython) or Python libs

   ``pybind11::lto`` / ``pybind11::thin_lto``
     An alternative to `INTERPROCEDURAL_OPTIMIZATION` for adding link-time optimization.

   ``pybind11::windows_extras``
     ``/bigobj`` and ``/mp`` for MSVC.

   ``pybind11::opt_size``
     ``/Os`` for MSVC, ``-Os`` for other compilers. Does nothing for debug builds.

Two helper functions are also provided:

    ``pybind11_strip(target)``
      Strips a target (uses ``CMAKE_STRIP`` after the target is built)

    ``pybind11_extension(target)``
      Sets the correct extension (with SOABI) for a target.

You can use these targets to build complex applications. For example, the
``add_python_module`` function is identical to:

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.4)
    project(example LANGUAGES CXX)

    find_package(pybind11 REQUIRED)  # or add_subdirectory(pybind11)

    add_library(example MODULE main.cpp)

    target_link_libraries(example PRIVATE pybind11::module pybind11::lto pybind11::windows_extras)

    pybind11_extension(example)
    pybind11_strip(example)

    set_target_properties(example PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                             CUDA_VISIBILITY_PRESET "hidden")

Instead of setting properties, you can set ``CMAKE_*`` variables to initialize these correctly.

.. warning::

    Since pybind11 is a metatemplate library, it is crucial that certain
    compiler flags are provided to ensure high quality code generation. In
    contrast to the ``pybind11_add_module()`` command, the CMake interface
    provides a *composable* set of targets to ensure that you retain flexibility.
    It can be expecially important to provide or set these properties; the
    :ref:`FAQ <faq:symhidden>` contains an explanation on why these are needed.

.. versionadded:: 2.6

.. _nopython-mode:

Advanced: NOPYTHON mode
-----------------------

If you want complete control, you can set ``PYBIND11_NOPYTHON`` to completely
disable Python integration (this also happens if you run ``FindPython2`` and
``FindPython3`` without running ``FindPython``). This gives you complete
freedom to integrate into an existing system (like `Scikit-Build's
<https://scikit-build.readthedocs.io>`_ ``PythonExtensions``).
``pybind11_add_module`` and ``pybind11_extension`` will be unavailable, and the
targets will be missing any Python specific behavior.

.. versionadded:: 2.6

Embedding the Python interpreter
--------------------------------

In addition to extension modules, pybind11 also supports embedding Python into
a C++ executable or library. In CMake, simply link with the ``pybind11::embed``
target. It provides everything needed to get the interpreter running. The Python
headers and libraries are attached to the target. Unlike ``pybind11::module``,
there is no need to manually set any additional properties here. For more
information about usage in C++, see :doc:`/advanced/embedding`.

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.4...3.18)
    project(example LANGUAGES CXX)

    find_package(pybind11 REQUIRED)  # or add_subdirectory(pybind11)

    add_executable(example main.cpp)
    target_link_libraries(example PRIVATE pybind11::embed)

.. _building_manually:

Building manually
=================

pybind11 is a header-only library, hence it is not necessary to link against
any special libraries and there are no intermediate (magic) translation steps.

On Linux, you can compile an example such as the one given in
:ref:`simple_example` using the following command:

.. code-block:: bash

    $ c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) example.cpp -o example$(python3-config --extension-suffix)

The flags given here assume that you're using Python 3. For Python 2, just
change the executable appropriately (to ``python`` or ``python2``).

The ``python3 -m pybind11 --includes`` command fetches the include paths for
both pybind11 and Python headers. This assumes that pybind11 has been installed
using ``pip`` or ``conda``. If it hasn't, you can also manually specify
``-I <path-to-pybind11>/include`` together with the Python includes path
``python3-config --includes``.

Note that Python 2.7 modules don't use a special suffix, so you should simply
use ``example.so`` instead of ``example$(python3-config --extension-suffix)``.
Besides, the ``--extension-suffix`` option may or may not be available, depending
on the distribution; in the latter case, the module extension can be manually
set to ``.so``.

On macOS: the build command is almost the same but it also requires passing
the ``-undefined dynamic_lookup`` flag so as to ignore missing symbols when
building the module:

.. code-block:: bash

    $ c++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup $(python3 -m pybind11 --includes) example.cpp -o example$(python3-config --extension-suffix)

In general, it is advisable to include several additional build parameters
that can considerably reduce the size of the created binary. Refer to section
:ref:`cmake` for a detailed example of a suitable cross-platform CMake-based
build system that works on all platforms including Windows.

.. note::

    On Linux and macOS, it's better to (intentionally) not link against
    ``libpython``. The symbols will be resolved when the extension library
    is loaded into a Python binary. This is preferable because you might
    have several different installations of a given Python version (e.g. the
    system-provided Python, and one that ships with a piece of commercial
    software). In this way, the plugin will work with both versions, instead
    of possibly importing a second Python library into a process that already
    contains one (which will lead to a segfault).


Building with Bazel
===================

You can build with the Bazel build system using the `pybind11_bazel
<https://github.com/pybind/pybind11_bazel>`_ repository.

Generating binding code automatically
=====================================

The ``Binder`` project is a tool for automatic generation of pybind11 binding
code by introspecting existing C++ codebases using LLVM/Clang. See the
[binder]_ documentation for details.

.. [binder] http://cppbinder.readthedocs.io/en/latest/about.html

[AutoWIG]_ is a Python library that wraps automatically compiled libraries into
high-level languages. It parses C++ code using LLVM/Clang technologies and
generates the wrappers using the Mako templating engine. The approach is automatic,
extensible, and applies to very complex C++ libraries, composed of thousands of
classes or incorporating modern meta-programming constructs.

.. [AutoWIG] https://github.com/StatisKit/AutoWIG

[robotpy-build]_ is a is a pure python, cross platform build tool that aims to
simplify creation of python wheels for pybind11 projects, and provide
cross-project dependency management. Additionally, it is able to autogenerate
customizable pybind11-based wrappers by parsing C++ header files.

.. [robotpy-build] https://robotpy-build.readthedocs.io
On version numbers
^^^^^^^^^^^^^^^^^^

The two version numbers (C++ and Python) must match when combined (checked when
you build the PyPI package), and must be a valid `PEP 440
<https://www.python.org/dev/peps/pep-0440>`_ version when combined.

For example:

.. code-block:: C++

    #define PYBIND11_VERSION_MAJOR X
    #define PYBIND11_VERSION_MINOR Y
    #define PYBIND11_VERSION_PATCH Z.dev1

For beta, ``PYBIND11_VERSION_PATCH`` should be ``Z.b1``. RC's can be ``Z.rc1``.
Always include the dot (even though PEP 440 allows it to be dropped). For a
final release, this must be a simple integer.


To release a new version of pybind11:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Update the version number
  - Update ``PYBIND11_VERSION_MAJOR`` etc. in
    ``include/pybind11/detail/common.h``. PATCH should be a simple integer.
  - Update ``pybind11/_version.py`` (match above)
  - Ensure that all the information in ``setup.cfg`` is up-to-date, like
    supported Python versions.
  - Add release date in ``docs/changelog.rst``.
      - Check to make sure
        `needs-changelog <https://github.com/pybind/pybind11/pulls?q=is%3Apr+is%3Aclosed+label%3A%22needs+changelog%22>`_
        issues are entered in the changelog (clear the label when done).
  - ``git add`` and ``git commit``, ``git push``. **Ensure CI passes**. (If it
    fails due to a known flake issue, either ignore or restart CI.)
- Add a release branch if this is a new minor version, or update the existing release branch if it is a patch version
  - New branch: ``git checkout -b vX.Y``, ``git push -u origin vX.Y``
  - Update branch: ``git checkout vX.Y``, ``git merge <release branch>``, ``git push``
- Update tags (optional; if you skip this, the GitHub release makes a
  non-annotated tag for you)
  - ``git tag -a vX.Y.Z -m 'vX.Y.Z release'``.
  - ``git push --tags``.
- Update stable
    - ``git checkout stable``
    - ``git merge master``
    - ``git push``
- Make a GitHub release (this shows up in the UI, sends new release
  notifications to users watching releases, and also uploads PyPI packages).
  (Note: if you do not use an existing tag, this creates a new lightweight tag
  for you, so you could skip the above step).
  - GUI method: click "Create a new release" on the far right, fill in the tag
    name (if you didn't tag above, it will be made here), fill in a release
    name like "Version X.Y.Z", and optionally copy-and-paste the changelog into
    the description (processed as markdown by Pandoc). Check "pre-release" if
    this is a beta/RC. You can get partway there with
    ``cat docs/changelog.rst | pandsoc -f rst -t markdown``.
  - CLI method: with ``gh`` installed, run ``gh release create vX.Y.Z -t "Version X.Y.Z"``
    If this is a pre-release, add ``-p``.

- Get back to work
  - Make sure you are on master, not somewhere else: ``git checkout master``
  - Update version macros in ``include/pybind11/detail/common.h`` (set PATCH to
    ``0.dev1`` and increment MINOR).
  - Update ``_version.py`` to match
  - Add a spot for in-development updates in ``docs/changelog.rst``.
  - ``git add``, ``git commit``, ``git push``

If a version branch is updated, remember to set PATCH to ``1.dev1``.

If you'd like to bump homebrew, run:

.. code-block::

    brew bump-formula-pr --url https://github.com/pybind/pybind11/archive/vX.Y.Z.tar.gz

Conda-forge should automatically make a PR in a few hours, and automatically
merge it if there are no issues.


Manual packaging
^^^^^^^^^^^^^^^^

If you need to manually upload releases, you can download the releases from the job artifacts and upload them with twine. You can also make the files locally (not recommended in general, as your local directory is more likely to be "dirty" and SDists love picking up random unrelated/hidden files); this is the procedure:

.. code-block:: bash

    python3 -m pip install build
    python3 -m build
    PYBIND11_SDIST_GLOBAL=1 python3 -m build
    twine upload dist/*

This makes SDists and wheels, and the final line uploads them.
.. _changelog:

Changelog
#########

Starting with version 1.8.0, pybind11 releases use a `semantic versioning
<http://semver.org>`_ policy.


v2.6.2 (Jan 26, 2021)
---------------------

Minor missing functionality added:

* enum: add missing Enum.value property.
  `#2739 <https://github.com/pybind/pybind11/pull/2739>`_

* Allow thread termination to be avoided during shutdown for CPython 3.7+ via
  ``.disarm`` for ``gil_scoped_acquire``/``gil_scoped_release``.
  `#2657 <https://github.com/pybind/pybind11/pull/2657>`_

Fixed or improved behavior in a few special cases:

* Fix bug where the constructor of ``object`` subclasses would not throw on
  being passed a Python object of the wrong type.
  `#2701 <https://github.com/pybind/pybind11/pull/2701>`_

* The ``type_caster`` for integers does not convert Python objects with
  ``__int__`` anymore with ``noconvert`` or during the first round of trying
  overloads.
  `#2698 <https://github.com/pybind/pybind11/pull/2698>`_

* When casting to a C++ integer, ``__index__`` is always called and not
  considered as conversion, consistent with Python 3.8+.
  `#2801 <https://github.com/pybind/pybind11/pull/2801>`_

Build improvements:

* Setup helpers: ``extra_compile_args`` and ``extra_link_args`` automatically set by
  Pybind11Extension are now prepended, which allows them to be overridden
  by user-set ``extra_compile_args`` and ``extra_link_args``.
  `#2808 <https://github.com/pybind/pybind11/pull/2808>`_

* Setup helpers: Don't trigger unused parameter warning.
  `#2735 <https://github.com/pybind/pybind11/pull/2735>`_

* CMake: Support running with ``--warn-uninitialized`` active.
  `#2806 <https://github.com/pybind/pybind11/pull/2806>`_

* CMake: Avoid error if included from two submodule directories.
  `#2804 <https://github.com/pybind/pybind11/pull/2804>`_

* CMake: Fix ``STATIC`` / ``SHARED`` being ignored in FindPython mode.
  `#2796 <https://github.com/pybind/pybind11/pull/2796>`_

* CMake: Respect the setting for ``CMAKE_CXX_VISIBILITY_PRESET`` if defined.
  `#2793 <https://github.com/pybind/pybind11/pull/2793>`_

* CMake: Fix issue with FindPython2/FindPython3 not working with ``pybind11::embed``.
  `#2662 <https://github.com/pybind/pybind11/pull/2662>`_

* CMake: mixing local and installed pybind11's would prioritize the installed
  one over the local one (regression in 2.6.0).
  `#2716 <https://github.com/pybind/pybind11/pull/2716>`_


Bug fixes:

* Fixed segfault in multithreaded environments when using
  ``scoped_ostream_redirect``.
  `#2675 <https://github.com/pybind/pybind11/pull/2675>`_

* Leave docstring unset when all docstring-related options are disabled, rather
  than set an empty string.
  `#2745 <https://github.com/pybind/pybind11/pull/2745>`_

* The module key in builtins that pybind11 uses to store its internals changed
  from std::string to a python str type (more natural on Python 2, no change on
  Python 3).
  `#2814 <https://github.com/pybind/pybind11/pull/2814>`_

* Fixed assertion error related to unhandled (later overwritten) exception in
  CPython 3.8 and 3.9 debug builds.
  `#2685 <https://github.com/pybind/pybind11/pull/2685>`_

* Fix ``py::gil_scoped_acquire`` assert with CPython 3.9 debug build.
  `#2683 <https://github.com/pybind/pybind11/pull/2683>`_

* Fix issue with a test failing on PyTest 6.2.
  `#2741 <https://github.com/pybind/pybind11/pull/2741>`_

Warning fixes:

* Fix warning modifying constructor parameter 'flag' that shadows a field of
  'set_flag' ``[-Wshadow-field-in-constructor-modified]``.
  `#2780 <https://github.com/pybind/pybind11/pull/2780>`_

* Suppressed some deprecation warnings about old-style
  ``__init__``/``__setstate__`` in the tests.
  `#2759 <https://github.com/pybind/pybind11/pull/2759>`_

Valgrind work:

* Fix invalid access when calling a pybind11 ``__init__`` on a non-pybind11
  class instance.
  `#2755 <https://github.com/pybind/pybind11/pull/2755>`_

* Fixed various minor memory leaks in pybind11's test suite.
  `#2758 <https://github.com/pybind/pybind11/pull/2758>`_

* Resolved memory leak in cpp_function initialization when exceptions occurred.
  `#2756 <https://github.com/pybind/pybind11/pull/2756>`_

* Added a Valgrind build, checking for leaks and memory-related UB, to CI.
  `#2746 <https://github.com/pybind/pybind11/pull/2746>`_

Compiler support:

* Intel compiler was not activating C++14 support due to a broken define.
  `#2679 <https://github.com/pybind/pybind11/pull/2679>`_

* Support ICC and NVIDIA HPC SDK in C++17 mode.
  `#2729 <https://github.com/pybind/pybind11/pull/2729>`_

* Support Intel OneAPI compiler (ICC 20.2) and add to CI.
  `#2573 <https://github.com/pybind/pybind11/pull/2573>`_



v2.6.1 (Nov 11, 2020)
---------------------

* ``py::exec``, ``py::eval``, and ``py::eval_file`` now add the builtins module
  as ``"__builtins__"`` to their ``globals`` argument, better matching ``exec``
  and ``eval`` in pure Python.
  `#2616 <https://github.com/pybind/pybind11/pull/2616>`_

* ``setup_helpers`` will no longer set a minimum macOS version higher than the
  current version.
  `#2622 <https://github.com/pybind/pybind11/pull/2622>`_

* Allow deleting static properties.
  `#2629 <https://github.com/pybind/pybind11/pull/2629>`_

* Seal a leak in ``def_buffer``, cleaning up the ``capture`` object after the
  ``class_`` object goes out of scope.
  `#2634 <https://github.com/pybind/pybind11/pull/2634>`_

* ``pybind11_INCLUDE_DIRS`` was incorrect, potentially causing a regression if
  it was expected to include ``PYTHON_INCLUDE_DIRS`` (please use targets
  instead).
  `#2636 <https://github.com/pybind/pybind11/pull/2636>`_

* Added parameter names to the ``py::enum_`` constructor and methods, avoiding
  ``arg0`` in the generated docstrings.
  `#2637 <https://github.com/pybind/pybind11/pull/2637>`_

* Added ``needs_recompile`` optional function to the ``ParallelCompiler``
  helper, to allow a recompile to be skipped based on a user-defined function.
  `#2643 <https://github.com/pybind/pybind11/pull/2643>`_


v2.6.0 (Oct 21, 2020)
---------------------

See :ref:`upgrade-guide-2.6` for help upgrading to the new version.

New features:

* Keyword-only arguments supported in Python 2 or 3 with ``py::kw_only()``.
  `#2100 <https://github.com/pybind/pybind11/pull/2100>`_

* Positional-only arguments supported in Python 2 or 3 with ``py::pos_only()``.
  `#2459 <https://github.com/pybind/pybind11/pull/2459>`_

* ``py::is_final()`` class modifier to block subclassing (CPython only).
  `#2151 <https://github.com/pybind/pybind11/pull/2151>`_

* Added ``py::prepend()``, allowing a function to be placed at the beginning of
  the overload chain.
  `#1131 <https://github.com/pybind/pybind11/pull/1131>`_

* Access to the type object now provided with ``py::type::of<T>()`` and
  ``py::type::of(h)``.
  `#2364 <https://github.com/pybind/pybind11/pull/2364>`_

* Perfect forwarding support for methods.
  `#2048 <https://github.com/pybind/pybind11/pull/2048>`_

* Added ``py::error_already_set::discard_as_unraisable()``.
  `#2372 <https://github.com/pybind/pybind11/pull/2372>`_

* ``py::hash`` is now public.
  `#2217 <https://github.com/pybind/pybind11/pull/2217>`_

* ``py::class_<union_type>`` is now supported. Note that writing to one data
  member of the union and reading another (type punning) is UB in C++. Thus
  pybind11-bound enums should never be used for such conversions.
  `#2320 <https://github.com/pybind/pybind11/pull/2320>`_.

* Classes now check local scope when registering members, allowing a subclass
  to have a member with the same name as a parent (such as an enum).
  `#2335 <https://github.com/pybind/pybind11/pull/2335>`_

Code correctness features:

* Error now thrown when ``__init__`` is forgotten on subclasses.
  `#2152 <https://github.com/pybind/pybind11/pull/2152>`_

* Throw error if conversion to a pybind11 type if the Python object isn't a
  valid instance of that type, such as ``py::bytes(o)`` when ``py::object o``
  isn't a bytes instance.
  `#2349 <https://github.com/pybind/pybind11/pull/2349>`_

* Throw if conversion to ``str`` fails.
  `#2477 <https://github.com/pybind/pybind11/pull/2477>`_


API changes:

* ``py::module`` was renamed ``py::module_`` to avoid issues with C++20 when
  used unqualified, but an alias ``py::module`` is provided for backward
  compatibility.
  `#2489 <https://github.com/pybind/pybind11/pull/2489>`_

* Public constructors for ``py::module_`` have been deprecated; please use
  ``pybind11::module_::create_extension_module`` if you were using the public
  constructor (fairly rare after ``PYBIND11_MODULE`` was introduced).
  `#2552 <https://github.com/pybind/pybind11/pull/2552>`_

* ``PYBIND11_OVERLOAD*`` macros and ``get_overload`` function replaced by
  correctly-named ``PYBIND11_OVERRIDE*`` and ``get_override``, fixing
  inconsistencies in the presence of a closing ``;`` in these macros.
  ``get_type_overload`` is deprecated.
  `#2325 <https://github.com/pybind/pybind11/pull/2325>`_

Packaging / building improvements:

* The Python package was reworked to be more powerful and useful.
  `#2433 <https://github.com/pybind/pybind11/pull/2433>`_

  * :ref:`build-setuptools` is easier thanks to a new
    ``pybind11.setup_helpers`` module, which provides utilities to use
    setuptools with pybind11. It can be used via PEP 518, ``setup_requires``,
    or by directly importing or copying ``setup_helpers.py`` into your project.

  * CMake configuration files are now included in the Python package. Use
    ``pybind11.get_cmake_dir()`` or ``python -m pybind11 --cmakedir`` to get
    the directory with the CMake configuration files, or include the
    site-packages location in your ``CMAKE_MODULE_PATH``. Or you can use the
    new ``pybind11[global]`` extra when you install ``pybind11``, which
    installs the CMake files and headers into your base environment in the
    standard location.

  * ``pybind11-config`` is another way to write ``python -m pybind11`` if you
    have your PATH set up.

  * Added external typing support to the helper module, code from
    ``import pybind11`` can now be type checked.
    `#2588 <https://github.com/pybind/pybind11/pull/2588>`_

* Minimum CMake required increased to 3.4.
  `#2338 <https://github.com/pybind/pybind11/pull/2338>`_ and
  `#2370 <https://github.com/pybind/pybind11/pull/2370>`_

  * Full integration with CMake’s C++ standard system and compile features
    replaces ``PYBIND11_CPP_STANDARD``.

  * Generated config file is now portable to different Python/compiler/CMake
    versions.

  * Virtual environments prioritized if ``PYTHON_EXECUTABLE`` is not set
    (``venv``, ``virtualenv``, and ``conda``) (similar to the new FindPython
    mode).

  * Other CMake features now natively supported, like
    ``CMAKE_INTERPROCEDURAL_OPTIMIZATION``, ``set(CMAKE_CXX_VISIBILITY_PRESET
    hidden)``.

  * ``CUDA`` as a language is now supported.

  * Helper functions ``pybind11_strip``, ``pybind11_extension``,
    ``pybind11_find_import`` added, see :doc:`cmake/index`.

  * Optional :ref:`find-python-mode` and :ref:`nopython-mode` with CMake.
    `#2370 <https://github.com/pybind/pybind11/pull/2370>`_

* Uninstall target added.
  `#2265 <https://github.com/pybind/pybind11/pull/2265>`_ and
  `#2346 <https://github.com/pybind/pybind11/pull/2346>`_

* ``pybind11_add_module()`` now accepts an optional ``OPT_SIZE`` flag that
  switches the binding target to size-based optimization if the global build
  type can not always be fixed to ``MinSizeRel`` (except in debug mode, where
  optimizations remain disabled).  ``MinSizeRel`` or this flag reduces binary
  size quite substantially (~25% on some platforms).
  `#2463 <https://github.com/pybind/pybind11/pull/2463>`_

Smaller or developer focused features and fixes:

* Moved ``mkdoc.py`` to a new repo, `pybind11-mkdoc`_. There are no longer
  submodules in the main repo.

* ``py::memoryview`` segfault fix and update, with new
  ``py::memoryview::from_memory`` in Python 3, and documentation.
  `#2223 <https://github.com/pybind/pybind11/pull/2223>`_

* Fix for ``buffer_info`` on Python 2.
  `#2503 <https://github.com/pybind/pybind11/pull/2503>`_

* If ``__eq__`` defined but not ``__hash__``, ``__hash__`` is now set to
  ``None``.
  `#2291 <https://github.com/pybind/pybind11/pull/2291>`_

* ``py::ellipsis`` now also works on Python 2.
  `#2360 <https://github.com/pybind/pybind11/pull/2360>`_

* Pointer to ``std::tuple`` & ``std::pair`` supported in cast.
  `#2334 <https://github.com/pybind/pybind11/pull/2334>`_

* Small fixes in NumPy support. ``py::array`` now uses ``py::ssize_t`` as first
  argument type.
  `#2293 <https://github.com/pybind/pybind11/pull/2293>`_

* Added missing signature for ``py::array``.
  `#2363 <https://github.com/pybind/pybind11/pull/2363>`_

* ``unchecked_mutable_reference`` has access to operator ``()`` and ``[]`` when
  const.
  `#2514 <https://github.com/pybind/pybind11/pull/2514>`_

* ``py::vectorize`` is now supported on functions that return void.
  `#1969 <https://github.com/pybind/pybind11/pull/1969>`_

* ``py::capsule`` supports ``get_pointer`` and ``set_pointer``.
  `#1131 <https://github.com/pybind/pybind11/pull/1131>`_

* Fix crash when different instances share the same pointer of the same type.
  `#2252 <https://github.com/pybind/pybind11/pull/2252>`_

* Fix for ``py::len`` not clearing Python's error state when it fails and throws.
  `#2575 <https://github.com/pybind/pybind11/pull/2575>`_

* Bugfixes related to more extensive testing, new GitHub Actions CI.
  `#2321 <https://github.com/pybind/pybind11/pull/2321>`_

* Bug in timezone issue in Eastern hemisphere midnight fixed.
  `#2438 <https://github.com/pybind/pybind11/pull/2438>`_

* ``std::chrono::time_point`` now works when the resolution is not the same as
  the system.
  `#2481 <https://github.com/pybind/pybind11/pull/2481>`_

* Bug fixed where ``py::array_t`` could accept arrays that did not match the
  requested ordering.
  `#2484 <https://github.com/pybind/pybind11/pull/2484>`_

* Avoid a segfault on some compilers when types are removed in Python.
  `#2564 <https://github.com/pybind/pybind11/pull/2564>`_

* ``py::arg::none()`` is now also respected when passing keyword arguments.
  `#2611 <https://github.com/pybind/pybind11/pull/2611>`_

* PyPy fixes, PyPy 7.3.x now supported, including PyPy3. (Known issue with
  PyPy2 and Windows `#2596 <https://github.com/pybind/pybind11/issues/2596>`_).
  `#2146 <https://github.com/pybind/pybind11/pull/2146>`_

* CPython 3.9.0 workaround for undefined behavior (macOS segfault).
  `#2576 <https://github.com/pybind/pybind11/pull/2576>`_

* CPython 3.9 warning fixes.
  `#2253 <https://github.com/pybind/pybind11/pull/2253>`_

* Improved C++20 support, now tested in CI.
  `#2489 <https://github.com/pybind/pybind11/pull/2489>`_
  `#2599 <https://github.com/pybind/pybind11/pull/2599>`_

* Improved but still incomplete debug Python interpreter support.
  `#2025 <https://github.com/pybind/pybind11/pull/2025>`_

* NVCC (CUDA 11) now supported and tested in CI.
  `#2461 <https://github.com/pybind/pybind11/pull/2461>`_

* NVIDIA PGI compilers now supported and tested in CI.
  `#2475 <https://github.com/pybind/pybind11/pull/2475>`_

* At least Intel 18 now explicitly required when compiling with Intel.
  `#2577 <https://github.com/pybind/pybind11/pull/2577>`_

* Extensive style checking in CI, with `pre-commit`_ support. Code
  modernization, checked by clang-tidy.

* Expanded docs, including new main page, new installing section, and CMake
  helpers page, along with over a dozen new sections on existing pages.

* In GitHub, new docs for contributing and new issue templates.

.. _pre-commit: https://pre-commit.com

.. _pybind11-mkdoc: https://github.com/pybind/pybind11-mkdoc

v2.5.0 (Mar 31, 2020)
-----------------------------------------------------

* Use C++17 fold expressions in type casters, if available. This can
  improve performance during overload resolution when functions have
  multiple arguments.
  `#2043 <https://github.com/pybind/pybind11/pull/2043>`_.

* Changed include directory resolution in ``pybind11/__init__.py``
  and installation in ``setup.py``. This fixes a number of open issues
  where pybind11 headers could not be found in certain environments.
  `#1995 <https://github.com/pybind/pybind11/pull/1995>`_.

* C++20 ``char8_t`` and ``u8string`` support. `#2026
  <https://github.com/pybind/pybind11/pull/2026>`_.

* CMake: search for Python 3.9. `bb9c91
  <https://github.com/pybind/pybind11/commit/bb9c91>`_.

* Fixes for MSYS-based build environments.
  `#2087 <https://github.com/pybind/pybind11/pull/2087>`_,
  `#2053 <https://github.com/pybind/pybind11/pull/2053>`_.

* STL bindings for ``std::vector<...>::clear``. `#2074
  <https://github.com/pybind/pybind11/pull/2074>`_.

* Read-only flag for ``py::buffer``. `#1466
  <https://github.com/pybind/pybind11/pull/1466>`_.

* Exception handling during module initialization.
  `bf2b031 <https://github.com/pybind/pybind11/commit/bf2b031>`_.

* Support linking against a CPython debug build.
  `#2025 <https://github.com/pybind/pybind11/pull/2025>`_.

* Fixed issues involving the availability and use of aligned ``new`` and
  ``delete``. `#1988 <https://github.com/pybind/pybind11/pull/1988>`_,
  `759221 <https://github.com/pybind/pybind11/commit/759221>`_.

* Fixed a resource leak upon interpreter shutdown.
  `#2020 <https://github.com/pybind/pybind11/pull/2020>`_.

* Fixed error handling in the boolean caster.
  `#1976 <https://github.com/pybind/pybind11/pull/1976>`_.

v2.4.3 (Oct 15, 2019)
-----------------------------------------------------

* Adapt pybind11 to a C API convention change in Python 3.8. `#1950
  <https://github.com/pybind/pybind11/pull/1950>`_.

v2.4.2 (Sep 21, 2019)
-----------------------------------------------------

* Replaced usage of a C++14 only construct. `#1929
  <https://github.com/pybind/pybind11/pull/1929>`_.

* Made an ifdef future-proof for Python >= 4. `f3109d
  <https://github.com/pybind/pybind11/commit/f3109d>`_.

v2.4.1 (Sep 20, 2019)
-----------------------------------------------------

* Fixed a problem involving implicit conversion from enumerations to integers
  on Python 3.8. `#1780 <https://github.com/pybind/pybind11/pull/1780>`_.

v2.4.0 (Sep 19, 2019)
-----------------------------------------------------

* Try harder to keep pybind11-internal data structures separate when there
  are potential ABI incompatibilities. Fixes crashes that occurred when loading
  multiple pybind11 extensions that were e.g. compiled by GCC (libstdc++)
  and Clang (libc++).
  `#1588 <https://github.com/pybind/pybind11/pull/1588>`_ and
  `c9f5a <https://github.com/pybind/pybind11/commit/c9f5a>`_.

* Added support for ``__await__``, ``__aiter__``, and ``__anext__`` protocols.
  `#1842 <https://github.com/pybind/pybind11/pull/1842>`_.

* ``pybind11_add_module()``: don't strip symbols when compiling in
  ``RelWithDebInfo`` mode. `#1980
  <https://github.com/pybind/pybind11/pull/1980>`_.

* ``enum_``: Reproduce Python behavior when comparing against invalid values
  (e.g. ``None``, strings, etc.). Add back support for ``__invert__()``.
  `#1912 <https://github.com/pybind/pybind11/pull/1912>`_,
  `#1907 <https://github.com/pybind/pybind11/pull/1907>`_.

* List insertion operation for ``py::list``.
  Added ``.empty()`` to all collection types.
  Added ``py::set::contains()`` and ``py::dict::contains()``.
  `#1887 <https://github.com/pybind/pybind11/pull/1887>`_,
  `#1884 <https://github.com/pybind/pybind11/pull/1884>`_,
  `#1888 <https://github.com/pybind/pybind11/pull/1888>`_.

* ``py::details::overload_cast_impl`` is available in C++11 mode, can be used
  like ``overload_cast`` with an additional set of parantheses.
  `#1581 <https://github.com/pybind/pybind11/pull/1581>`_.

* Fixed ``get_include()`` on Conda.
  `#1877 <https://github.com/pybind/pybind11/pull/1877>`_.

* ``stl_bind.h``: negative indexing support.
  `#1882 <https://github.com/pybind/pybind11/pull/1882>`_.

* Minor CMake fix to add MinGW compatibility.
  `#1851 <https://github.com/pybind/pybind11/pull/1851>`_.

* GIL-related fixes.
  `#1836 <https://github.com/pybind/pybind11/pull/1836>`_,
  `8b90b <https://github.com/pybind/pybind11/commit/8b90b>`_.

* Other very minor/subtle fixes and improvements.
  `#1329 <https://github.com/pybind/pybind11/pull/1329>`_,
  `#1910 <https://github.com/pybind/pybind11/pull/1910>`_,
  `#1863 <https://github.com/pybind/pybind11/pull/1863>`_,
  `#1847 <https://github.com/pybind/pybind11/pull/1847>`_,
  `#1890 <https://github.com/pybind/pybind11/pull/1890>`_,
  `#1860 <https://github.com/pybind/pybind11/pull/1860>`_,
  `#1848 <https://github.com/pybind/pybind11/pull/1848>`_,
  `#1821 <https://github.com/pybind/pybind11/pull/1821>`_,
  `#1837 <https://github.com/pybind/pybind11/pull/1837>`_,
  `#1833 <https://github.com/pybind/pybind11/pull/1833>`_,
  `#1748 <https://github.com/pybind/pybind11/pull/1748>`_,
  `#1852 <https://github.com/pybind/pybind11/pull/1852>`_.

v2.3.0 (June 11, 2019)
-----------------------------------------------------

* Significantly reduced module binary size (10-20%) when compiled in C++11 mode
  with GCC/Clang, or in any mode with MSVC. Function signatures are now always
  precomputed at compile time (this was previously only available in C++14 mode
  for non-MSVC compilers).
  `#934 <https://github.com/pybind/pybind11/pull/934>`_.

* Add basic support for tag-based static polymorphism, where classes
  provide a method to returns the desired type of an instance.
  `#1326 <https://github.com/pybind/pybind11/pull/1326>`_.

* Python type wrappers (``py::handle``, ``py::object``, etc.)
  now support map Python's number protocol onto C++ arithmetic
  operators such as ``operator+``, ``operator/=``, etc.
  `#1511 <https://github.com/pybind/pybind11/pull/1511>`_.

* A number of improvements related to enumerations:

   1. The ``enum_`` implementation was rewritten from scratch to reduce
      code bloat. Rather than instantiating a full implementation for each
      enumeration, most code is now contained in a generic base class.
      `#1511 <https://github.com/pybind/pybind11/pull/1511>`_.

   2. The ``value()``  method of ``py::enum_`` now accepts an optional
      docstring that will be shown in the documentation of the associated
      enumeration. `#1160 <https://github.com/pybind/pybind11/pull/1160>`_.

   3. check for already existing enum value and throw an error if present.
      `#1453 <https://github.com/pybind/pybind11/pull/1453>`_.

* Support for over-aligned type allocation via C++17's aligned ``new``
  statement. `#1582 <https://github.com/pybind/pybind11/pull/1582>`_.

* Added ``py::ellipsis()`` method for slicing of multidimensional NumPy arrays
  `#1502 <https://github.com/pybind/pybind11/pull/1502>`_.

* Numerous Improvements to the ``mkdoc.py`` script for extracting documentation
  from C++ header files.
  `#1788 <https://github.com/pybind/pybind11/pull/1788>`_.

* ``pybind11_add_module()``: allow including Python as a ``SYSTEM`` include path.
  `#1416 <https://github.com/pybind/pybind11/pull/1416>`_.

* ``pybind11/stl.h`` does not convert strings to ``vector<string>`` anymore.
  `#1258 <https://github.com/pybind/pybind11/issues/1258>`_.

* Mark static methods as such to fix auto-generated Sphinx documentation.
  `#1732 <https://github.com/pybind/pybind11/pull/1732>`_.

* Re-throw forced unwind exceptions (e.g. during pthread termination).
  `#1208 <https://github.com/pybind/pybind11/pull/1208>`_.

* Added ``__contains__`` method to the bindings of maps (``std::map``,
  ``std::unordered_map``).
  `#1767 <https://github.com/pybind/pybind11/pull/1767>`_.

* Improvements to ``gil_scoped_acquire``.
  `#1211 <https://github.com/pybind/pybind11/pull/1211>`_.

* Type caster support for ``std::deque<T>``.
  `#1609 <https://github.com/pybind/pybind11/pull/1609>`_.

* Support for ``std::unique_ptr`` holders, whose deleters differ between a base and derived
  class. `#1353 <https://github.com/pybind/pybind11/pull/1353>`_.

* Construction of STL array/vector-like data structures from
  iterators. Added an ``extend()`` operation.
  `#1709 <https://github.com/pybind/pybind11/pull/1709>`_,

* CMake build system improvements for projects that include non-C++
  files (e.g. plain C, CUDA) in ``pybind11_add_module`` et al.
  `#1678 <https://github.com/pybind/pybind11/pull/1678>`_.

* Fixed asynchronous invocation and deallocation of Python functions
  wrapped in ``std::function``.
  `#1595 <https://github.com/pybind/pybind11/pull/1595>`_.

* Fixes regarding return value policy propagation in STL type casters.
  `#1603 <https://github.com/pybind/pybind11/pull/1603>`_.

* Fixed scoped enum comparisons.
  `#1571 <https://github.com/pybind/pybind11/pull/1571>`_.

* Fixed iostream redirection for code that releases the GIL.
  `#1368 <https://github.com/pybind/pybind11/pull/1368>`_,

* A number of CI-related fixes.
  `#1757 <https://github.com/pybind/pybind11/pull/1757>`_,
  `#1744 <https://github.com/pybind/pybind11/pull/1744>`_,
  `#1670 <https://github.com/pybind/pybind11/pull/1670>`_.

v2.2.4 (September 11, 2018)
-----------------------------------------------------

* Use new Python 3.7 Thread Specific Storage (TSS) implementation if available.
  `#1454 <https://github.com/pybind/pybind11/pull/1454>`_,
  `#1517 <https://github.com/pybind/pybind11/pull/1517>`_.

* Fixes for newer MSVC versions and C++17 mode.
  `#1347 <https://github.com/pybind/pybind11/pull/1347>`_,
  `#1462 <https://github.com/pybind/pybind11/pull/1462>`_.

* Propagate return value policies to type-specific casters
  when casting STL containers.
  `#1455 <https://github.com/pybind/pybind11/pull/1455>`_.

* Allow ostream-redirection of more than 1024 characters.
  `#1479 <https://github.com/pybind/pybind11/pull/1479>`_.

* Set ``Py_DEBUG`` define when compiling against a debug Python build.
  `#1438 <https://github.com/pybind/pybind11/pull/1438>`_.

* Untangle integer logic in number type caster to work for custom
  types that may only be castable to a restricted set of builtin types.
  `#1442 <https://github.com/pybind/pybind11/pull/1442>`_.

* CMake build system: Remember Python version in cache file.
  `#1434 <https://github.com/pybind/pybind11/pull/1434>`_.

* Fix for custom smart pointers: use ``std::addressof`` to obtain holder
  address instead of ``operator&``.
  `#1435 <https://github.com/pybind/pybind11/pull/1435>`_.

* Properly report exceptions thrown during module initialization.
  `#1362 <https://github.com/pybind/pybind11/pull/1362>`_.

* Fixed a segmentation fault when creating empty-shaped NumPy array.
  `#1371 <https://github.com/pybind/pybind11/pull/1371>`_.

* The version of Intel C++ compiler must be >= 2017, and this is now checked by
  the header files. `#1363 <https://github.com/pybind/pybind11/pull/1363>`_.

* A few minor typo fixes and improvements to the test suite, and
  patches that silence compiler warnings.

* Vectors now support construction from generators, as well as ``extend()`` from a
  list or generator.
  `#1496 <https://github.com/pybind/pybind11/pull/1496>`_.


v2.2.3 (April 29, 2018)
-----------------------------------------------------

* The pybind11 header location detection was replaced by a new implementation
  that no longer depends on ``pip`` internals (the recently released ``pip``
  10 has restricted access to this API).
  `#1190 <https://github.com/pybind/pybind11/pull/1190>`_.

* Small adjustment to an implementation detail to work around a compiler segmentation fault in Clang 3.3/3.4.
  `#1350 <https://github.com/pybind/pybind11/pull/1350>`_.

* The minimal supported version of the Intel compiler was >= 17.0 since
  pybind11 v2.1. This check is now explicit, and a compile-time error is raised
  if the compiler meet the requirement.
  `#1363 <https://github.com/pybind/pybind11/pull/1363>`_.

* Fixed an endianness-related fault in the test suite.
  `#1287 <https://github.com/pybind/pybind11/pull/1287>`_.

v2.2.2 (February 7, 2018)
-----------------------------------------------------

* Fixed a segfault when combining embedded interpreter
  shutdown/reinitialization with external loaded pybind11 modules.
  `#1092 <https://github.com/pybind/pybind11/pull/1092>`_.

* Eigen support: fixed a bug where Nx1/1xN numpy inputs couldn't be passed as
  arguments to Eigen vectors (which for Eigen are simply compile-time fixed
  Nx1/1xN matrices).
  `#1106 <https://github.com/pybind/pybind11/pull/1106>`_.

* Clarified to license by moving the licensing of contributions from
  ``LICENSE`` into ``CONTRIBUTING.md``: the licensing of contributions is not
  actually part of the software license as distributed.  This isn't meant to be
  a substantial change in the licensing of the project, but addresses concerns
  that the clause made the license non-standard.
  `#1109 <https://github.com/pybind/pybind11/issues/1109>`_.

* Fixed a regression introduced in 2.1 that broke binding functions with lvalue
  character literal arguments.
  `#1128 <https://github.com/pybind/pybind11/pull/1128>`_.

* MSVC: fix for compilation failures under /permissive-, and added the flag to
  the appveyor test suite.
  `#1155 <https://github.com/pybind/pybind11/pull/1155>`_.

* Fixed ``__qualname__`` generation, and in turn, fixes how class names
  (especially nested class names) are shown in generated docstrings.
  `#1171 <https://github.com/pybind/pybind11/pull/1171>`_.

* Updated the FAQ with a suggested project citation reference.
  `#1189 <https://github.com/pybind/pybind11/pull/1189>`_.

* Added fixes for deprecation warnings when compiled under C++17 with
  ``-Wdeprecated`` turned on, and add ``-Wdeprecated`` to the test suite
  compilation flags.
  `#1191 <https://github.com/pybind/pybind11/pull/1191>`_.

* Fixed outdated PyPI URLs in ``setup.py``.
  `#1213 <https://github.com/pybind/pybind11/pull/1213>`_.

* Fixed a refcount leak for arguments that end up in a ``py::args`` argument
  for functions with both fixed positional and ``py::args`` arguments.
  `#1216 <https://github.com/pybind/pybind11/pull/1216>`_.

* Fixed a potential segfault resulting from possible premature destruction of
  ``py::args``/``py::kwargs`` arguments with overloaded functions.
  `#1223 <https://github.com/pybind/pybind11/pull/1223>`_.

* Fixed ``del map[item]`` for a ``stl_bind.h`` bound stl map.
  `#1229 <https://github.com/pybind/pybind11/pull/1229>`_.

* Fixed a regression from v2.1.x where the aggregate initialization could
  unintentionally end up at a constructor taking a templated
  ``std::initializer_list<T>`` argument.
  `#1249 <https://github.com/pybind/pybind11/pull/1249>`_.

* Fixed an issue where calling a function with a keep_alive policy on the same
  nurse/patient pair would cause the internal patient storage to needlessly
  grow (unboundedly, if the nurse is long-lived).
  `#1251 <https://github.com/pybind/pybind11/issues/1251>`_.

* Various other minor fixes.

v2.2.1 (September 14, 2017)
-----------------------------------------------------

* Added ``py::module_::reload()`` member function for reloading a module.
  `#1040 <https://github.com/pybind/pybind11/pull/1040>`_.

* Fixed a reference leak in the number converter.
  `#1078 <https://github.com/pybind/pybind11/pull/1078>`_.

* Fixed compilation with Clang on host GCC < 5 (old libstdc++ which isn't fully
  C++11 compliant). `#1062 <https://github.com/pybind/pybind11/pull/1062>`_.

* Fixed a regression where the automatic ``std::vector<bool>`` caster would
  fail to compile. The same fix also applies to any container which returns
  element proxies instead of references.
  `#1053 <https://github.com/pybind/pybind11/pull/1053>`_.

* Fixed a regression where the ``py::keep_alive`` policy could not be applied
  to constructors. `#1065 <https://github.com/pybind/pybind11/pull/1065>`_.

* Fixed a nullptr dereference when loading a ``py::module_local`` type
  that's only registered in an external module.
  `#1058 <https://github.com/pybind/pybind11/pull/1058>`_.

* Fixed implicit conversion of accessors to types derived from ``py::object``.
  `#1076 <https://github.com/pybind/pybind11/pull/1076>`_.

* The ``name`` in ``PYBIND11_MODULE(name, variable)`` can now be a macro.
  `#1082 <https://github.com/pybind/pybind11/pull/1082>`_.

* Relaxed overly strict ``py::pickle()`` check for matching get and set types.
  `#1064 <https://github.com/pybind/pybind11/pull/1064>`_.

* Conversion errors now try to be more informative when it's likely that
  a missing header is the cause (e.g. forgetting ``<pybind11/stl.h>``).
  `#1077 <https://github.com/pybind/pybind11/pull/1077>`_.

v2.2.0 (August 31, 2017)
-----------------------------------------------------

* Support for embedding the Python interpreter. See the
  :doc:`documentation page </advanced/embedding>` for a
  full overview of the new features.
  `#774 <https://github.com/pybind/pybind11/pull/774>`_,
  `#889 <https://github.com/pybind/pybind11/pull/889>`_,
  `#892 <https://github.com/pybind/pybind11/pull/892>`_,
  `#920 <https://github.com/pybind/pybind11/pull/920>`_.

  .. code-block:: cpp

      #include <pybind11/embed.h>
      namespace py = pybind11;

      int main() {
          py::scoped_interpreter guard{}; // start the interpreter and keep it alive

          py::print("Hello, World!"); // use the Python API
      }

* Support for inheriting from multiple C++ bases in Python.
  `#693 <https://github.com/pybind/pybind11/pull/693>`_.

  .. code-block:: python

      from cpp_module import CppBase1, CppBase2

      class PyDerived(CppBase1, CppBase2):
          def __init__(self):
              CppBase1.__init__(self)  # C++ bases must be initialized explicitly
              CppBase2.__init__(self)

* ``PYBIND11_MODULE`` is now the preferred way to create module entry points.
  ``PYBIND11_PLUGIN`` is deprecated. See :ref:`macros` for details.
  `#879 <https://github.com/pybind/pybind11/pull/879>`_.

  .. code-block:: cpp

      // new
      PYBIND11_MODULE(example, m) {
          m.def("add", [](int a, int b) { return a + b; });
      }

      // old
      PYBIND11_PLUGIN(example) {
          py::module m("example");
          m.def("add", [](int a, int b) { return a + b; });
          return m.ptr();
      }

* pybind11's headers and build system now more strictly enforce hidden symbol
  visibility for extension modules. This should be seamless for most users,
  but see the :doc:`upgrade` if you use a custom build system.
  `#995 <https://github.com/pybind/pybind11/pull/995>`_.

* Support for ``py::module_local`` types which allow multiple modules to
  export the same C++ types without conflicts. This is useful for opaque
  types like ``std::vector<int>``. ``py::bind_vector`` and ``py::bind_map``
  now default to ``py::module_local`` if their elements are builtins or
  local types. See :ref:`module_local` for details.
  `#949 <https://github.com/pybind/pybind11/pull/949>`_,
  `#981 <https://github.com/pybind/pybind11/pull/981>`_,
  `#995 <https://github.com/pybind/pybind11/pull/995>`_,
  `#997 <https://github.com/pybind/pybind11/pull/997>`_.

* Custom constructors can now be added very easily using lambdas or factory
  functions which return a class instance by value, pointer or holder. This
  supersedes the old placement-new ``__init__`` technique.
  See :ref:`custom_constructors` for details.
  `#805 <https://github.com/pybind/pybind11/pull/805>`_,
  `#1014 <https://github.com/pybind/pybind11/pull/1014>`_.

  .. code-block:: cpp

      struct Example {
          Example(std::string);
      };

      py::class_<Example>(m, "Example")
          .def(py::init<std::string>()) // existing constructor
          .def(py::init([](int n) { // custom constructor
              return std::make_unique<Example>(std::to_string(n));
          }));

* Similarly to custom constructors, pickling support functions are now bound
  using the ``py::pickle()`` adaptor which improves type safety. See the
  :doc:`upgrade` and :ref:`pickling` for details.
  `#1038 <https://github.com/pybind/pybind11/pull/1038>`_.

* Builtin support for converting C++17 standard library types and general
  conversion improvements:

  1. C++17 ``std::variant`` is supported right out of the box. C++11/14
     equivalents (e.g. ``boost::variant``) can also be added with a simple
     user-defined specialization. See :ref:`cpp17_container_casters` for details.
     `#811 <https://github.com/pybind/pybind11/pull/811>`_,
     `#845 <https://github.com/pybind/pybind11/pull/845>`_,
     `#989 <https://github.com/pybind/pybind11/pull/989>`_.

  2. Out-of-the-box support for C++17 ``std::string_view``.
     `#906 <https://github.com/pybind/pybind11/pull/906>`_.

  3. Improved compatibility of the builtin ``optional`` converter.
     `#874 <https://github.com/pybind/pybind11/pull/874>`_.

  4. The ``bool`` converter now accepts ``numpy.bool_`` and types which
     define ``__bool__`` (Python 3.x) or ``__nonzero__`` (Python 2.7).
     `#925 <https://github.com/pybind/pybind11/pull/925>`_.

  5. C++-to-Python casters are now more efficient and move elements out
     of rvalue containers whenever possible.
     `#851 <https://github.com/pybind/pybind11/pull/851>`_,
     `#936 <https://github.com/pybind/pybind11/pull/936>`_,
     `#938 <https://github.com/pybind/pybind11/pull/938>`_.

  6. Fixed ``bytes`` to ``std::string/char*`` conversion on Python 3.
     `#817 <https://github.com/pybind/pybind11/pull/817>`_.

  7. Fixed lifetime of temporary C++ objects created in Python-to-C++ conversions.
     `#924 <https://github.com/pybind/pybind11/pull/924>`_.

* Scope guard call policy for RAII types, e.g. ``py::call_guard<py::gil_scoped_release>()``,
  ``py::call_guard<py::scoped_ostream_redirect>()``. See :ref:`call_policies` for details.
  `#740 <https://github.com/pybind/pybind11/pull/740>`_.

* Utility for redirecting C++ streams to Python (e.g. ``std::cout`` ->
  ``sys.stdout``). Scope guard ``py::scoped_ostream_redirect`` in C++ and
  a context manager in Python. See :ref:`ostream_redirect`.
  `#1009 <https://github.com/pybind/pybind11/pull/1009>`_.

* Improved handling of types and exceptions across module boundaries.
  `#915 <https://github.com/pybind/pybind11/pull/915>`_,
  `#951 <https://github.com/pybind/pybind11/pull/951>`_,
  `#995 <https://github.com/pybind/pybind11/pull/995>`_.

* Fixed destruction order of ``py::keep_alive`` nurse/patient objects
  in reference cycles.
  `#856 <https://github.com/pybind/pybind11/pull/856>`_.

* NumPy and buffer protocol related improvements:

  1. Support for negative strides in Python buffer objects/numpy arrays. This
     required changing integers from unsigned to signed for the related C++ APIs.
     Note: If you have compiler warnings enabled, you may notice some new conversion
     warnings after upgrading. These can be resolved with ``static_cast``.
     `#782 <https://github.com/pybind/pybind11/pull/782>`_.

  2. Support ``std::complex`` and arrays inside ``PYBIND11_NUMPY_DTYPE``.
     `#831 <https://github.com/pybind/pybind11/pull/831>`_,
     `#832 <https://github.com/pybind/pybind11/pull/832>`_.

  3. Support for constructing ``py::buffer_info`` and ``py::arrays`` using
     arbitrary containers or iterators instead of requiring a ``std::vector``.
     `#788 <https://github.com/pybind/pybind11/pull/788>`_,
     `#822 <https://github.com/pybind/pybind11/pull/822>`_,
     `#860 <https://github.com/pybind/pybind11/pull/860>`_.

  4. Explicitly check numpy version and require >= 1.7.0.
     `#819 <https://github.com/pybind/pybind11/pull/819>`_.

* Support for allowing/prohibiting ``None`` for specific arguments and improved
  ``None`` overload resolution order. See :ref:`none_arguments` for details.
  `#843 <https://github.com/pybind/pybind11/pull/843>`_.
  `#859 <https://github.com/pybind/pybind11/pull/859>`_.

* Added ``py::exec()`` as a shortcut for ``py::eval<py::eval_statements>()``
  and support for C++11 raw string literals as input. See :ref:`eval`.
  `#766 <https://github.com/pybind/pybind11/pull/766>`_,
  `#827 <https://github.com/pybind/pybind11/pull/827>`_.

* ``py::vectorize()`` ignores non-vectorizable arguments and supports
  member functions.
  `#762 <https://github.com/pybind/pybind11/pull/762>`_.

* Support for bound methods as callbacks (``pybind11/functional.h``).
  `#815 <https://github.com/pybind/pybind11/pull/815>`_.

* Allow aliasing pybind11 methods: ``cls.attr("foo") = cls.attr("bar")``.
  `#802 <https://github.com/pybind/pybind11/pull/802>`_.

* Don't allow mixed static/non-static overloads.
  `#804 <https://github.com/pybind/pybind11/pull/804>`_.

* Fixed overriding static properties in derived classes.
  `#784 <https://github.com/pybind/pybind11/pull/784>`_.

* Added support for write only properties.
  `#1144 <https://github.com/pybind/pybind11/pull/1144>`_.

* Improved deduction of member functions of a derived class when its bases
  aren't registered with pybind11.
  `#855 <https://github.com/pybind/pybind11/pull/855>`_.

  .. code-block:: cpp

      struct Base {
          int foo() { return 42; }
      }

      struct Derived : Base {}

      // Now works, but previously required also binding `Base`
      py::class_<Derived>(m, "Derived")
          .def("foo", &Derived::foo); // function is actually from `Base`

* The implementation of ``py::init<>`` now uses C++11 brace initialization
  syntax to construct instances, which permits binding implicit constructors of
  aggregate types. `#1015 <https://github.com/pybind/pybind11/pull/1015>`_.

    .. code-block:: cpp

        struct Aggregate {
            int a;
            std::string b;
        };

        py::class_<Aggregate>(m, "Aggregate")
            .def(py::init<int, const std::string &>());

* Fixed issues with multiple inheritance with offset base/derived pointers.
  `#812 <https://github.com/pybind/pybind11/pull/812>`_,
  `#866 <https://github.com/pybind/pybind11/pull/866>`_,
  `#960 <https://github.com/pybind/pybind11/pull/960>`_.

* Fixed reference leak of type objects.
  `#1030 <https://github.com/pybind/pybind11/pull/1030>`_.

* Improved support for the ``/std:c++14`` and ``/std:c++latest`` modes
  on MSVC 2017.
  `#841 <https://github.com/pybind/pybind11/pull/841>`_,
  `#999 <https://github.com/pybind/pybind11/pull/999>`_.

* Fixed detection of private operator new on MSVC.
  `#893 <https://github.com/pybind/pybind11/pull/893>`_,
  `#918 <https://github.com/pybind/pybind11/pull/918>`_.

* Intel C++ compiler compatibility fixes.
  `#937 <https://github.com/pybind/pybind11/pull/937>`_.

* Fixed implicit conversion of `py::enum_` to integer types on Python 2.7.
  `#821 <https://github.com/pybind/pybind11/pull/821>`_.

* Added ``py::hash`` to fetch the hash value of Python objects, and
  ``.def(hash(py::self))`` to provide the C++ ``std::hash`` as the Python
  ``__hash__`` method.
  `#1034 <https://github.com/pybind/pybind11/pull/1034>`_.

* Fixed ``__truediv__`` on Python 2 and ``__itruediv__`` on Python 3.
  `#867 <https://github.com/pybind/pybind11/pull/867>`_.

* ``py::capsule`` objects now support the ``name`` attribute. This is useful
  for interfacing with ``scipy.LowLevelCallable``.
  `#902 <https://github.com/pybind/pybind11/pull/902>`_.

* Fixed ``py::make_iterator``'s ``__next__()`` for past-the-end calls.
  `#897 <https://github.com/pybind/pybind11/pull/897>`_.

* Added ``error_already_set::matches()`` for checking Python exceptions.
  `#772 <https://github.com/pybind/pybind11/pull/772>`_.

* Deprecated ``py::error_already_set::clear()``. It's no longer needed
  following a simplification of the ``py::error_already_set`` class.
  `#954 <https://github.com/pybind/pybind11/pull/954>`_.

* Deprecated ``py::handle::operator==()`` in favor of ``py::handle::is()``
  `#825 <https://github.com/pybind/pybind11/pull/825>`_.

* Deprecated ``py::object::borrowed``/``py::object::stolen``.
  Use ``py::object::borrowed_t{}``/``py::object::stolen_t{}`` instead.
  `#771 <https://github.com/pybind/pybind11/pull/771>`_.

* Changed internal data structure versioning to avoid conflicts between
  modules compiled with different revisions of pybind11.
  `#1012 <https://github.com/pybind/pybind11/pull/1012>`_.

* Additional compile-time and run-time error checking and more informative messages.
  `#786 <https://github.com/pybind/pybind11/pull/786>`_,
  `#794 <https://github.com/pybind/pybind11/pull/794>`_,
  `#803 <https://github.com/pybind/pybind11/pull/803>`_.

* Various minor improvements and fixes.
  `#764 <https://github.com/pybind/pybind11/pull/764>`_,
  `#791 <https://github.com/pybind/pybind11/pull/791>`_,
  `#795 <https://github.com/pybind/pybind11/pull/795>`_,
  `#840 <https://github.com/pybind/pybind11/pull/840>`_,
  `#844 <https://github.com/pybind/pybind11/pull/844>`_,
  `#846 <https://github.com/pybind/pybind11/pull/846>`_,
  `#849 <https://github.com/pybind/pybind11/pull/849>`_,
  `#858 <https://github.com/pybind/pybind11/pull/858>`_,
  `#862 <https://github.com/pybind/pybind11/pull/862>`_,
  `#871 <https://github.com/pybind/pybind11/pull/871>`_,
  `#872 <https://github.com/pybind/pybind11/pull/872>`_,
  `#881 <https://github.com/pybind/pybind11/pull/881>`_,
  `#888 <https://github.com/pybind/pybind11/pull/888>`_,
  `#899 <https://github.com/pybind/pybind11/pull/899>`_,
  `#928 <https://github.com/pybind/pybind11/pull/928>`_,
  `#931 <https://github.com/pybind/pybind11/pull/931>`_,
  `#944 <https://github.com/pybind/pybind11/pull/944>`_,
  `#950 <https://github.com/pybind/pybind11/pull/950>`_,
  `#952 <https://github.com/pybind/pybind11/pull/952>`_,
  `#962 <https://github.com/pybind/pybind11/pull/962>`_,
  `#965 <https://github.com/pybind/pybind11/pull/965>`_,
  `#970 <https://github.com/pybind/pybind11/pull/970>`_,
  `#978 <https://github.com/pybind/pybind11/pull/978>`_,
  `#979 <https://github.com/pybind/pybind11/pull/979>`_,
  `#986 <https://github.com/pybind/pybind11/pull/986>`_,
  `#1020 <https://github.com/pybind/pybind11/pull/1020>`_,
  `#1027 <https://github.com/pybind/pybind11/pull/1027>`_,
  `#1037 <https://github.com/pybind/pybind11/pull/1037>`_.

* Testing improvements.
  `#798 <https://github.com/pybind/pybind11/pull/798>`_,
  `#882 <https://github.com/pybind/pybind11/pull/882>`_,
  `#898 <https://github.com/pybind/pybind11/pull/898>`_,
  `#900 <https://github.com/pybind/pybind11/pull/900>`_,
  `#921 <https://github.com/pybind/pybind11/pull/921>`_,
  `#923 <https://github.com/pybind/pybind11/pull/923>`_,
  `#963 <https://github.com/pybind/pybind11/pull/963>`_.

v2.1.1 (April 7, 2017)
-----------------------------------------------------

* Fixed minimum version requirement for MSVC 2015u3
  `#773 <https://github.com/pybind/pybind11/pull/773>`_.

v2.1.0 (March 22, 2017)
-----------------------------------------------------

* pybind11 now performs function overload resolution in two phases. The first
  phase only considers exact type matches, while the second allows for implicit
  conversions to take place. A special ``noconvert()`` syntax can be used to
  completely disable implicit conversions for specific arguments.
  `#643 <https://github.com/pybind/pybind11/pull/643>`_,
  `#634 <https://github.com/pybind/pybind11/pull/634>`_,
  `#650 <https://github.com/pybind/pybind11/pull/650>`_.

* Fixed a regression where static properties no longer worked with classes
  using multiple inheritance. The ``py::metaclass`` attribute is no longer
  necessary (and deprecated as of this release) when binding classes with
  static properties.
  `#679 <https://github.com/pybind/pybind11/pull/679>`_,

* Classes bound using ``pybind11`` can now use custom metaclasses.
  `#679 <https://github.com/pybind/pybind11/pull/679>`_,

* ``py::args`` and ``py::kwargs`` can now be mixed with other positional
  arguments when binding functions using pybind11.
  `#611 <https://github.com/pybind/pybind11/pull/611>`_.

* Improved support for C++11 unicode string and character types; added
  extensive documentation regarding pybind11's string conversion behavior.
  `#624 <https://github.com/pybind/pybind11/pull/624>`_,
  `#636 <https://github.com/pybind/pybind11/pull/636>`_,
  `#715 <https://github.com/pybind/pybind11/pull/715>`_.

* pybind11 can now avoid expensive copies when converting Eigen arrays to NumPy
  arrays (and vice versa). `#610 <https://github.com/pybind/pybind11/pull/610>`_.

* The "fast path" in ``py::vectorize`` now works for any full-size group of C or
  F-contiguous arrays. The non-fast path is also faster since it no longer performs
  copies of the input arguments (except when type conversions are necessary).
  `#610 <https://github.com/pybind/pybind11/pull/610>`_.

* Added fast, unchecked access to NumPy arrays via a proxy object.
  `#746 <https://github.com/pybind/pybind11/pull/746>`_.

* Transparent support for class-specific ``operator new`` and
  ``operator delete`` implementations.
  `#755 <https://github.com/pybind/pybind11/pull/755>`_.

* Slimmer and more efficient STL-compatible iterator interface for sequence types.
  `#662 <https://github.com/pybind/pybind11/pull/662>`_.

* Improved custom holder type support.
  `#607 <https://github.com/pybind/pybind11/pull/607>`_.

* ``nullptr`` to ``None`` conversion fixed in various builtin type casters.
  `#732 <https://github.com/pybind/pybind11/pull/732>`_.

* ``enum_`` now exposes its members via a special ``__members__`` attribute.
  `#666 <https://github.com/pybind/pybind11/pull/666>`_.

* ``std::vector`` bindings created using ``stl_bind.h`` can now optionally
  implement the buffer protocol. `#488 <https://github.com/pybind/pybind11/pull/488>`_.

* Automated C++ reference documentation using doxygen and breathe.
  `#598 <https://github.com/pybind/pybind11/pull/598>`_.

* Added minimum compiler version assertions.
  `#727 <https://github.com/pybind/pybind11/pull/727>`_.

* Improved compatibility with C++1z.
  `#677 <https://github.com/pybind/pybind11/pull/677>`_.

* Improved ``py::capsule`` API. Can be used to implement cleanup
  callbacks that are involved at module destruction time.
  `#752 <https://github.com/pybind/pybind11/pull/752>`_.

* Various minor improvements and fixes.
  `#595 <https://github.com/pybind/pybind11/pull/595>`_,
  `#588 <https://github.com/pybind/pybind11/pull/588>`_,
  `#589 <https://github.com/pybind/pybind11/pull/589>`_,
  `#603 <https://github.com/pybind/pybind11/pull/603>`_,
  `#619 <https://github.com/pybind/pybind11/pull/619>`_,
  `#648 <https://github.com/pybind/pybind11/pull/648>`_,
  `#695 <https://github.com/pybind/pybind11/pull/695>`_,
  `#720 <https://github.com/pybind/pybind11/pull/720>`_,
  `#723 <https://github.com/pybind/pybind11/pull/723>`_,
  `#729 <https://github.com/pybind/pybind11/pull/729>`_,
  `#724 <https://github.com/pybind/pybind11/pull/724>`_,
  `#742 <https://github.com/pybind/pybind11/pull/742>`_,
  `#753 <https://github.com/pybind/pybind11/pull/753>`_.

v2.0.1 (Jan 4, 2017)
-----------------------------------------------------

* Fix pointer to reference error in type_caster on MSVC
  `#583 <https://github.com/pybind/pybind11/pull/583>`_.

* Fixed a segmentation in the test suite due to a typo
  `cd7eac <https://github.com/pybind/pybind11/commit/cd7eac>`_.

v2.0.0 (Jan 1, 2017)
-----------------------------------------------------

* Fixed a reference counting regression affecting types with custom metaclasses
  (introduced in v2.0.0-rc1).
  `#571 <https://github.com/pybind/pybind11/pull/571>`_.

* Quenched a CMake policy warning.
  `#570 <https://github.com/pybind/pybind11/pull/570>`_.

v2.0.0-rc1 (Dec 23, 2016)
-----------------------------------------------------

The pybind11 developers are excited to issue a release candidate of pybind11
with a subsequent v2.0.0 release planned in early January next year.

An incredible amount of effort by went into pybind11 over the last ~5 months,
leading to a release that is jam-packed with exciting new features and numerous
usability improvements. The following list links PRs or individual commits
whenever applicable.

Happy Christmas!

* Support for binding C++ class hierarchies that make use of multiple
  inheritance. `#410 <https://github.com/pybind/pybind11/pull/410>`_.

* PyPy support: pybind11 now supports nightly builds of PyPy and will
  interoperate with the future 5.7 release. No code changes are necessary,
  everything "just" works as usual. Note that we only target the Python 2.7
  branch for now; support for 3.x will be added once its ``cpyext`` extension
  support catches up. A few minor features remain unsupported for the time
  being (notably dynamic attributes in custom types).
  `#527 <https://github.com/pybind/pybind11/pull/527>`_.

* Significant work on the documentation -- in particular, the monolithic
  ``advanced.rst`` file was restructured into a easier to read hierarchical
  organization. `#448 <https://github.com/pybind/pybind11/pull/448>`_.

* Many NumPy-related improvements:

  1. Object-oriented API to access and modify NumPy ``ndarray`` instances,
     replicating much of the corresponding NumPy C API functionality.
     `#402 <https://github.com/pybind/pybind11/pull/402>`_.

  2. NumPy array ``dtype`` array descriptors are now first-class citizens and
     are exposed via a new class ``py::dtype``.

  3. Structured dtypes can be registered using the ``PYBIND11_NUMPY_DTYPE()``
     macro. Special ``array`` constructors accepting dtype objects were also
     added.

     One potential caveat involving this change: format descriptor strings
     should now be accessed via ``format_descriptor::format()`` (however, for
     compatibility purposes, the old syntax ``format_descriptor::value`` will
     still work for non-structured data types). `#308
     <https://github.com/pybind/pybind11/pull/308>`_.

  4. Further improvements to support structured dtypes throughout the system.
     `#472 <https://github.com/pybind/pybind11/pull/472>`_,
     `#474 <https://github.com/pybind/pybind11/pull/474>`_,
     `#459 <https://github.com/pybind/pybind11/pull/459>`_,
     `#453 <https://github.com/pybind/pybind11/pull/453>`_,
     `#452 <https://github.com/pybind/pybind11/pull/452>`_, and
     `#505 <https://github.com/pybind/pybind11/pull/505>`_.

  5. Fast access operators. `#497 <https://github.com/pybind/pybind11/pull/497>`_.

  6. Constructors for arrays whose storage is owned by another object.
     `#440 <https://github.com/pybind/pybind11/pull/440>`_.

  7. Added constructors for ``array`` and ``array_t`` explicitly accepting shape
     and strides; if strides are not provided, they are deduced assuming
     C-contiguity. Also added simplified constructors for 1-dimensional case.

  8. Added buffer/NumPy support for ``char[N]`` and ``std::array<char, N>`` types.

  9. Added ``memoryview`` wrapper type which is constructible from ``buffer_info``.

* Eigen: many additional conversions and support for non-contiguous
  arrays/slices.
  `#427 <https://github.com/pybind/pybind11/pull/427>`_,
  `#315 <https://github.com/pybind/pybind11/pull/315>`_,
  `#316 <https://github.com/pybind/pybind11/pull/316>`_,
  `#312 <https://github.com/pybind/pybind11/pull/312>`_, and
  `#267 <https://github.com/pybind/pybind11/pull/267>`_

* Incompatible changes in ``class_<...>::class_()``:

    1. Declarations of types that provide access via the buffer protocol must
       now include the ``py::buffer_protocol()`` annotation as an argument to
       the ``class_`` constructor.

    2. Declarations of types that require a custom metaclass (i.e. all classes
       which include static properties via commands such as
       ``def_readwrite_static()``) must now include the ``py::metaclass()``
       annotation as an argument to the ``class_`` constructor.

       These two changes were necessary to make type definitions in pybind11
       future-proof, and to support PyPy via its cpyext mechanism. `#527
       <https://github.com/pybind/pybind11/pull/527>`_.


    3. This version of pybind11 uses a redesigned mechanism for instantiating
       trampoline classes that are used to override virtual methods from within
       Python. This led to the following user-visible syntax change: instead of

       .. code-block:: cpp

           py::class_<TrampolineClass>("MyClass")
             .alias<MyClass>()
             ....

       write

       .. code-block:: cpp

           py::class_<MyClass, TrampolineClass>("MyClass")
             ....

       Importantly, both the original and the trampoline class are now
       specified as an arguments (in arbitrary order) to the ``py::class_``
       template, and the ``alias<..>()`` call is gone. The new scheme has zero
       overhead in cases when Python doesn't override any functions of the
       underlying C++ class. `rev. 86d825
       <https://github.com/pybind/pybind11/commit/86d825>`_.

* Added ``eval`` and ``eval_file`` functions for evaluating expressions and
  statements from a string or file. `rev. 0d3fc3
  <https://github.com/pybind/pybind11/commit/0d3fc3>`_.

* pybind11 can now create types with a modifiable dictionary.
  `#437 <https://github.com/pybind/pybind11/pull/437>`_ and
  `#444 <https://github.com/pybind/pybind11/pull/444>`_.

* Support for translation of arbitrary C++ exceptions to Python counterparts.
  `#296 <https://github.com/pybind/pybind11/pull/296>`_ and
  `#273 <https://github.com/pybind/pybind11/pull/273>`_.

* Report full backtraces through mixed C++/Python code, better reporting for
  import errors, fixed GIL management in exception processing.
  `#537 <https://github.com/pybind/pybind11/pull/537>`_,
  `#494 <https://github.com/pybind/pybind11/pull/494>`_,
  `rev. e72d95 <https://github.com/pybind/pybind11/commit/e72d95>`_, and
  `rev. 099d6e <https://github.com/pybind/pybind11/commit/099d6e>`_.

* Support for bit-level operations, comparisons, and serialization of C++
  enumerations. `#503 <https://github.com/pybind/pybind11/pull/503>`_,
  `#508 <https://github.com/pybind/pybind11/pull/508>`_,
  `#380 <https://github.com/pybind/pybind11/pull/380>`_,
  `#309 <https://github.com/pybind/pybind11/pull/309>`_.
  `#311 <https://github.com/pybind/pybind11/pull/311>`_.

* The ``class_`` constructor now accepts its template arguments in any order.
  `#385 <https://github.com/pybind/pybind11/pull/385>`_.

* Attribute and item accessors now have a more complete interface which makes
  it possible to chain attributes as in
  ``obj.attr("a")[key].attr("b").attr("method")(1, 2, 3)``. `#425
  <https://github.com/pybind/pybind11/pull/425>`_.

* Major redesign of the default and conversion constructors in ``pytypes.h``.
  `#464 <https://github.com/pybind/pybind11/pull/464>`_.

* Added built-in support for ``std::shared_ptr`` holder type. It is no longer
  necessary to to include a declaration of the form
  ``PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)`` (though continuing to
  do so won't cause an error).
  `#454 <https://github.com/pybind/pybind11/pull/454>`_.

* New ``py::overload_cast`` casting operator to select among multiple possible
  overloads of a function. An example:

    .. code-block:: cpp

        py::class_<Pet>(m, "Pet")
            .def("set", py::overload_cast<int>(&Pet::set), "Set the pet's age")
            .def("set", py::overload_cast<const std::string &>(&Pet::set), "Set the pet's name");

  This feature only works on C++14-capable compilers.
  `#541 <https://github.com/pybind/pybind11/pull/541>`_.

* C++ types are automatically cast to Python types, e.g. when assigning
  them as an attribute. For instance, the following is now legal:

    .. code-block:: cpp

        py::module m = /* ... */
        m.attr("constant") = 123;

  (Previously, a ``py::cast`` call was necessary to avoid a compilation error.)
  `#551 <https://github.com/pybind/pybind11/pull/551>`_.

* Redesigned ``pytest``-based test suite. `#321 <https://github.com/pybind/pybind11/pull/321>`_.

* Instance tracking to detect reference leaks in test suite. `#324 <https://github.com/pybind/pybind11/pull/324>`_

* pybind11 can now distinguish between multiple different instances that are
  located at the same memory address, but which have different types.
  `#329 <https://github.com/pybind/pybind11/pull/329>`_.

* Improved logic in ``move`` return value policy.
  `#510 <https://github.com/pybind/pybind11/pull/510>`_,
  `#297 <https://github.com/pybind/pybind11/pull/297>`_.

* Generalized unpacking API to permit calling Python functions from C++ using
  notation such as ``foo(a1, a2, *args, "ka"_a=1, "kb"_a=2, **kwargs)``. `#372 <https://github.com/pybind/pybind11/pull/372>`_.

* ``py::print()`` function whose behavior matches that of the native Python
  ``print()`` function. `#372 <https://github.com/pybind/pybind11/pull/372>`_.

* Added ``py::dict`` keyword constructor:``auto d = dict("number"_a=42,
  "name"_a="World");``. `#372 <https://github.com/pybind/pybind11/pull/372>`_.

* Added ``py::str::format()`` method and ``_s`` literal: ``py::str s = "1 + 2
  = {}"_s.format(3);``. `#372 <https://github.com/pybind/pybind11/pull/372>`_.

* Added ``py::repr()`` function which is equivalent to Python's builtin
  ``repr()``. `#333 <https://github.com/pybind/pybind11/pull/333>`_.

* Improved construction and destruction logic for holder types. It is now
  possible to reference instances with smart pointer holder types without
  constructing the holder if desired. The ``PYBIND11_DECLARE_HOLDER_TYPE``
  macro now accepts an optional second parameter to indicate whether the holder
  type uses intrusive reference counting.
  `#533 <https://github.com/pybind/pybind11/pull/533>`_ and
  `#561 <https://github.com/pybind/pybind11/pull/561>`_.

* Mapping a stateless C++ function to Python and back is now "for free" (i.e.
  no extra indirections or argument conversion overheads). `rev. 954b79
  <https://github.com/pybind/pybind11/commit/954b79>`_.

* Bindings for ``std::valarray<T>``.
  `#545 <https://github.com/pybind/pybind11/pull/545>`_.

* Improved support for C++17 capable compilers.
  `#562 <https://github.com/pybind/pybind11/pull/562>`_.

* Bindings for ``std::optional<t>``.
  `#475 <https://github.com/pybind/pybind11/pull/475>`_,
  `#476 <https://github.com/pybind/pybind11/pull/476>`_,
  `#479 <https://github.com/pybind/pybind11/pull/479>`_,
  `#499 <https://github.com/pybind/pybind11/pull/499>`_, and
  `#501 <https://github.com/pybind/pybind11/pull/501>`_.

* ``stl_bind.h``: general improvements and support for ``std::map`` and
  ``std::unordered_map``.
  `#490 <https://github.com/pybind/pybind11/pull/490>`_,
  `#282 <https://github.com/pybind/pybind11/pull/282>`_,
  `#235 <https://github.com/pybind/pybind11/pull/235>`_.

* The ``std::tuple``, ``std::pair``, ``std::list``, and ``std::vector`` type
  casters now accept any Python sequence type as input. `rev. 107285
  <https://github.com/pybind/pybind11/commit/107285>`_.

* Improved CMake Python detection on multi-architecture Linux.
  `#532 <https://github.com/pybind/pybind11/pull/532>`_.

* Infrastructure to selectively disable or enable parts of the automatically
  generated docstrings. `#486 <https://github.com/pybind/pybind11/pull/486>`_.

* ``reference`` and ``reference_internal`` are now the default return value
  properties for static and non-static properties, respectively. `#473
  <https://github.com/pybind/pybind11/pull/473>`_. (the previous defaults
  were ``automatic``). `#473 <https://github.com/pybind/pybind11/pull/473>`_.

* Support for ``std::unique_ptr`` with non-default deleters or no deleter at
  all (``py::nodelete``). `#384 <https://github.com/pybind/pybind11/pull/384>`_.

* Deprecated ``handle::call()`` method. The new syntax to call Python
  functions is simply ``handle()``. It can also be invoked explicitly via
  ``handle::operator<X>()``, where ``X`` is an optional return value policy.

* Print more informative error messages when ``make_tuple()`` or ``cast()``
  fail. `#262 <https://github.com/pybind/pybind11/pull/262>`_.

* Creation of holder types for classes deriving from
  ``std::enable_shared_from_this<>`` now also works for ``const`` values.
  `#260 <https://github.com/pybind/pybind11/pull/260>`_.

* ``make_iterator()`` improvements for better compatibility with various
  types (now uses prefix increment operator); it now also accepts iterators
  with different begin/end types as long as they are equality comparable.
  `#247 <https://github.com/pybind/pybind11/pull/247>`_.

* ``arg()`` now accepts a wider range of argument types for default values.
  `#244 <https://github.com/pybind/pybind11/pull/244>`_.

* Support ``keep_alive`` where the nurse object may be ``None``. `#341
  <https://github.com/pybind/pybind11/pull/341>`_.

* Added constructors for ``str`` and ``bytes`` from zero-terminated char
  pointers, and from char pointers and length. Added constructors for ``str``
  from ``bytes`` and for ``bytes`` from ``str``, which will perform UTF-8
  decoding/encoding as required.

* Many other improvements of library internals without user-visible changes


1.8.1 (July 12, 2016)
----------------------
* Fixed a rare but potentially very severe issue when the garbage collector ran
  during pybind11 type creation.

1.8.0 (June 14, 2016)
----------------------
* Redesigned CMake build system which exports a convenient
  ``pybind11_add_module`` function to parent projects.
* ``std::vector<>`` type bindings analogous to Boost.Python's ``indexing_suite``
* Transparent conversion of sparse and dense Eigen matrices and vectors (``eigen.h``)
* Added an ``ExtraFlags`` template argument to the NumPy ``array_t<>`` wrapper
  to disable an enforced cast that may lose precision, e.g. to create overloads
  for different precisions and complex vs real-valued matrices.
* Prevent implicit conversion of floating point values to integral types in
  function arguments
* Fixed incorrect default return value policy for functions returning a shared
  pointer
* Don't allow registering a type via ``class_`` twice
* Don't allow casting a ``None`` value into a C++ lvalue reference
* Fixed a crash in ``enum_::operator==`` that was triggered by the ``help()`` command
* Improved detection of whether or not custom C++ types can be copy/move-constructed
* Extended ``str`` type to also work with ``bytes`` instances
* Added a ``"name"_a`` user defined string literal that is equivalent to ``py::arg("name")``.
* When specifying function arguments via ``py::arg``, the test that verifies
  the number of arguments now runs at compile time.
* Added ``[[noreturn]]`` attribute to ``pybind11_fail()`` to quench some
  compiler warnings
* List function arguments in exception text when the dispatch code cannot find
  a matching overload
* Added ``PYBIND11_OVERLOAD_NAME`` and ``PYBIND11_OVERLOAD_PURE_NAME`` macros which
  can be used to override virtual methods whose name differs in C++ and Python
  (e.g. ``__call__`` and ``operator()``)
* Various minor ``iterator`` and ``make_iterator()`` improvements
* Transparently support ``__bool__`` on Python 2.x and Python 3.x
* Fixed issue with destructor of unpickled object not being called
* Minor CMake build system improvements on Windows
* New ``pybind11::args`` and ``pybind11::kwargs`` types to create functions which
  take an arbitrary number of arguments and keyword arguments
* New syntax to call a Python function from C++ using ``*args`` and ``*kwargs``
* The functions ``def_property_*`` now correctly process docstring arguments (these
  formerly caused a segmentation fault)
* Many ``mkdoc.py`` improvements (enumerations, template arguments, ``DOC()``
  macro accepts more arguments)
* Cygwin support
* Documentation improvements (pickling support, ``keep_alive``, macro usage)

1.7 (April 30, 2016)
----------------------
* Added a new ``move`` return value policy that triggers C++11 move semantics.
  The automatic return value policy falls back to this case whenever a rvalue
  reference is encountered
* Significantly more general GIL state routines that are used instead of
  Python's troublesome ``PyGILState_Ensure`` and ``PyGILState_Release`` API
* Redesign of opaque types that drastically simplifies their usage
* Extended ability to pass values of type ``[const] void *``
* ``keep_alive`` fix: don't fail when there is no patient
* ``functional.h``: acquire the GIL before calling a Python function
* Added Python RAII type wrappers ``none`` and ``iterable``
* Added ``*args`` and ``*kwargs`` pass-through parameters to
  ``pybind11.get_include()`` function
* Iterator improvements and fixes
* Documentation on return value policies and opaque types improved

1.6 (April 30, 2016)
----------------------
* Skipped due to upload to PyPI gone wrong and inability to recover
  (https://github.com/pypa/packaging-problems/issues/74)

1.5 (April 21, 2016)
----------------------
* For polymorphic types, use RTTI to try to return the closest type registered with pybind11
* Pickling support for serializing and unserializing C++ instances to a byte stream in Python
* Added a convenience routine ``make_iterator()`` which turns a range indicated
  by a pair of C++ iterators into a iterable Python object
* Added ``len()`` and a variadic ``make_tuple()`` function
* Addressed a rare issue that could confuse the current virtual function
  dispatcher and another that could lead to crashes in multi-threaded
  applications
* Added a ``get_include()`` function to the Python module that returns the path
  of the directory containing the installed pybind11 header files
* Documentation improvements: import issues, symbol visibility, pickling, limitations
* Added casting support for ``std::reference_wrapper<>``

1.4 (April 7, 2016)
--------------------------
* Transparent type conversion for ``std::wstring`` and ``wchar_t``
* Allow passing ``nullptr``-valued strings
* Transparent passing of ``void *`` pointers using capsules
* Transparent support for returning values wrapped in ``std::unique_ptr<>``
* Improved docstring generation for compatibility with Sphinx
* Nicer debug error message when default parameter construction fails
* Support for "opaque" types that bypass the transparent conversion layer for STL containers
* Redesigned type casting interface to avoid ambiguities that could occasionally cause compiler errors
* Redesigned property implementation; fixes crashes due to an unfortunate default return value policy
* Anaconda package generation support

1.3 (March 8, 2016)
--------------------------

* Added support for the Intel C++ compiler (v15+)
* Added support for the STL unordered set/map data structures
* Added support for the STL linked list data structure
* NumPy-style broadcasting support in ``pybind11::vectorize``
* pybind11 now displays more verbose error messages when ``arg::operator=()`` fails
* pybind11 internal data structures now live in a version-dependent namespace to avoid ABI issues
* Many, many bugfixes involving corner cases and advanced usage

1.2 (February 7, 2016)
--------------------------

* Optional: efficient generation of function signatures at compile time using C++14
* Switched to a simpler and more general way of dealing with function default
  arguments. Unused keyword arguments in function calls are now detected and
  cause errors as expected
* New ``keep_alive`` call policy analogous to Boost.Python's ``with_custodian_and_ward``
* New ``pybind11::base<>`` attribute to indicate a subclass relationship
* Improved interface for RAII type wrappers in ``pytypes.h``
* Use RAII type wrappers consistently within pybind11 itself. This
  fixes various potential refcount leaks when exceptions occur
* Added new ``bytes`` RAII type wrapper (maps to ``string`` in Python 2.7)
* Made handle and related RAII classes const correct, using them more
  consistently everywhere now
* Got rid of the ugly ``__pybind11__`` attributes on the Python side---they are
  now stored in a C++ hash table that is not visible in Python
* Fixed refcount leaks involving NumPy arrays and bound functions
* Vastly improved handling of shared/smart pointers
* Removed an unnecessary copy operation in ``pybind11::vectorize``
* Fixed naming clashes when both pybind11 and NumPy headers are included
* Added conversions for additional exception types
* Documentation improvements (using multiple extension modules, smart pointers,
  other minor clarifications)
* unified infrastructure for parsing variadic arguments in ``class_`` and cpp_function
* Fixed license text (was: ZLIB, should have been: 3-clause BSD)
* Python 3.2 compatibility
* Fixed remaining issues when accessing types in another plugin module
* Added enum comparison and casting methods
* Improved SFINAE-based detection of whether types are copy-constructible
* Eliminated many warnings about unused variables and the use of ``offsetof()``
* Support for ``std::array<>`` conversions

1.1 (December 7, 2015)
--------------------------

* Documentation improvements (GIL, wrapping functions, casting, fixed many typos)
* Generalized conversion of integer types
* Improved support for casting function objects
* Improved support for ``std::shared_ptr<>`` conversions
* Initial support for ``std::set<>`` conversions
* Fixed type resolution issue for types defined in a separate plugin module
* CMake build system improvements
* Factored out generic functionality to non-templated code (smaller code size)
* Added a code size / compile time benchmark vs Boost.Python
* Added an appveyor CI script

1.0 (October 15, 2015)
------------------------
* Initial release
Miscellaneous
#############

.. _macro_notes:

General notes regarding convenience macros
==========================================

pybind11 provides a few convenience macros such as
:func:`PYBIND11_DECLARE_HOLDER_TYPE` and ``PYBIND11_OVERRIDE_*``. Since these
are "just" macros that are evaluated in the preprocessor (which has no concept
of types), they *will* get confused by commas in a template argument; for
example, consider:

.. code-block:: cpp

    PYBIND11_OVERRIDE(MyReturnType<T1, T2>, Class<T3, T4>, func)

The limitation of the C preprocessor interprets this as five arguments (with new
arguments beginning after each comma) rather than three.  To get around this,
there are two alternatives: you can use a type alias, or you can wrap the type
using the ``PYBIND11_TYPE`` macro:

.. code-block:: cpp

    // Version 1: using a type alias
    using ReturnType = MyReturnType<T1, T2>;
    using ClassType = Class<T3, T4>;
    PYBIND11_OVERRIDE(ReturnType, ClassType, func);

    // Version 2: using the PYBIND11_TYPE macro:
    PYBIND11_OVERRIDE(PYBIND11_TYPE(MyReturnType<T1, T2>),
                      PYBIND11_TYPE(Class<T3, T4>), func)

The ``PYBIND11_MAKE_OPAQUE`` macro does *not* require the above workarounds.

.. _gil:

Global Interpreter Lock (GIL)
=============================

When calling a C++ function from Python, the GIL is always held.
The classes :class:`gil_scoped_release` and :class:`gil_scoped_acquire` can be
used to acquire and release the global interpreter lock in the body of a C++
function call. In this way, long-running C++ code can be parallelized using
multiple Python threads. Taking :ref:`overriding_virtuals` as an example, this
could be realized as follows (important changes highlighted):

.. code-block:: cpp
    :emphasize-lines: 8,9,31,32

    class PyAnimal : public Animal {
    public:
        /* Inherit the constructors */
        using Animal::Animal;

        /* Trampoline (need one for each virtual function) */
        std::string go(int n_times) {
            /* Acquire GIL before calling Python code */
            py::gil_scoped_acquire acquire;

            PYBIND11_OVERRIDE_PURE(
                std::string, /* Return type */
                Animal,      /* Parent class */
                go,          /* Name of function */
                n_times      /* Argument(s) */
            );
        }
    };

    PYBIND11_MODULE(example, m) {
        py::class_<Animal, PyAnimal> animal(m, "Animal");
        animal
            .def(py::init<>())
            .def("go", &Animal::go);

        py::class_<Dog>(m, "Dog", animal)
            .def(py::init<>());

        m.def("call_go", [](Animal *animal) -> std::string {
            /* Release GIL before calling into (potentially long-running) C++ code */
            py::gil_scoped_release release;
            return call_go(animal);
        });
    }

The ``call_go`` wrapper can also be simplified using the `call_guard` policy
(see :ref:`call_policies`) which yields the same result:

.. code-block:: cpp

    m.def("call_go", &call_go, py::call_guard<py::gil_scoped_release>());


Binding sequence data types, iterators, the slicing protocol, etc.
==================================================================

Please refer to the supplemental example for details.

.. seealso::

    The file :file:`tests/test_sequences_and_iterators.cpp` contains a
    complete example that shows how to bind a sequence data type, including
    length queries (``__len__``), iterators (``__iter__``), the slicing
    protocol and other kinds of useful operations.


Partitioning code over multiple extension modules
=================================================

It's straightforward to split binding code over multiple extension modules,
while referencing types that are declared elsewhere. Everything "just" works
without any special precautions. One exception to this rule occurs when
extending a type declared in another extension module. Recall the basic example
from Section :ref:`inheritance`.

.. code-block:: cpp

    py::class_<Pet> pet(m, "Pet");
    pet.def(py::init<const std::string &>())
       .def_readwrite("name", &Pet::name);

    py::class_<Dog>(m, "Dog", pet /* <- specify parent */)
        .def(py::init<const std::string &>())
        .def("bark", &Dog::bark);

Suppose now that ``Pet`` bindings are defined in a module named ``basic``,
whereas the ``Dog`` bindings are defined somewhere else. The challenge is of
course that the variable ``pet`` is not available anymore though it is needed
to indicate the inheritance relationship to the constructor of ``class_<Dog>``.
However, it can be acquired as follows:

.. code-block:: cpp

    py::object pet = (py::object) py::module_::import("basic").attr("Pet");

    py::class_<Dog>(m, "Dog", pet)
        .def(py::init<const std::string &>())
        .def("bark", &Dog::bark);

Alternatively, you can specify the base class as a template parameter option to
``class_``, which performs an automated lookup of the corresponding Python
type. Like the above code, however, this also requires invoking the ``import``
function once to ensure that the pybind11 binding code of the module ``basic``
has been executed:

.. code-block:: cpp

    py::module_::import("basic");

    py::class_<Dog, Pet>(m, "Dog")
        .def(py::init<const std::string &>())
        .def("bark", &Dog::bark);

Naturally, both methods will fail when there are cyclic dependencies.

Note that pybind11 code compiled with hidden-by-default symbol visibility (e.g.
via the command line flag ``-fvisibility=hidden`` on GCC/Clang), which is
required for proper pybind11 functionality, can interfere with the ability to
access types defined in another extension module.  Working around this requires
manually exporting types that are accessed by multiple extension modules;
pybind11 provides a macro to do just this:

.. code-block:: cpp

    class PYBIND11_EXPORT Dog : public Animal {
        ...
    };

Note also that it is possible (although would rarely be required) to share arbitrary
C++ objects between extension modules at runtime. Internal library data is shared
between modules using capsule machinery [#f6]_ which can be also utilized for
storing, modifying and accessing user-defined data. Note that an extension module
will "see" other extensions' data if and only if they were built with the same
pybind11 version. Consider the following example:

.. code-block:: cpp

    auto data = reinterpret_cast<MyData *>(py::get_shared_data("mydata"));
    if (!data)
        data = static_cast<MyData *>(py::set_shared_data("mydata", new MyData(42)));

If the above snippet was used in several separately compiled extension modules,
the first one to be imported would create a ``MyData`` instance and associate
a ``"mydata"`` key with a pointer to it. Extensions that are imported later
would be then able to access the data behind the same pointer.

.. [#f6] https://docs.python.org/3/extending/extending.html#using-capsules

Module Destructors
==================

pybind11 does not provide an explicit mechanism to invoke cleanup code at
module destruction time. In rare cases where such functionality is required, it
is possible to emulate it using Python capsules or weak references with a
destruction callback.

.. code-block:: cpp

    auto cleanup_callback = []() {
        // perform cleanup here -- this function is called with the GIL held
    };

    m.add_object("_cleanup", py::capsule(cleanup_callback));

This approach has the potential downside that instances of classes exposed
within the module may still be alive when the cleanup callback is invoked
(whether this is acceptable will generally depend on the application).

Alternatively, the capsule may also be stashed within a type object, which
ensures that it not called before all instances of that type have been
collected:

.. code-block:: cpp

    auto cleanup_callback = []() { /* ... */ };
    m.attr("BaseClass").attr("_cleanup") = py::capsule(cleanup_callback);

Both approaches also expose a potentially dangerous ``_cleanup`` attribute in
Python, which may be undesirable from an API standpoint (a premature explicit
call from Python might lead to undefined behavior). Yet another approach that
avoids this issue involves weak reference with a cleanup callback:

.. code-block:: cpp

    // Register a callback function that is invoked when the BaseClass object is collected
    py::cpp_function cleanup_callback(
        [](py::handle weakref) {
            // perform cleanup here -- this function is called with the GIL held

            weakref.dec_ref(); // release weak reference
        }
    );

    // Create a weak reference with a cleanup callback and initially leak it
    (void) py::weakref(m.attr("BaseClass"), cleanup_callback).release();

.. note::

    PyPy does not garbage collect objects when the interpreter exits. An alternative
    approach (which also works on CPython) is to use the :py:mod:`atexit` module [#f7]_,
    for example:

    .. code-block:: cpp

        auto atexit = py::module_::import("atexit");
        atexit.attr("register")(py::cpp_function([]() {
            // perform cleanup here -- this function is called with the GIL held
        }));

    .. [#f7] https://docs.python.org/3/library/atexit.html


Generating documentation using Sphinx
=====================================

Sphinx [#f4]_ has the ability to inspect the signatures and documentation
strings in pybind11-based extension modules to automatically generate beautiful
documentation in a variety formats. The python_example repository [#f5]_ contains a
simple example repository which uses this approach.

There are two potential gotchas when using this approach: first, make sure that
the resulting strings do not contain any :kbd:`TAB` characters, which break the
docstring parsing routines. You may want to use C++11 raw string literals,
which are convenient for multi-line comments. Conveniently, any excess
indentation will be automatically be removed by Sphinx. However, for this to
work, it is important that all lines are indented consistently, i.e.:

.. code-block:: cpp

    // ok
    m.def("foo", &foo, R"mydelimiter(
        The foo function

        Parameters
        ----------
    )mydelimiter");

    // *not ok*
    m.def("foo", &foo, R"mydelimiter(The foo function

        Parameters
        ----------
    )mydelimiter");

By default, pybind11 automatically generates and prepends a signature to the docstring of a function
registered with ``module_::def()`` and ``class_::def()``. Sometimes this
behavior is not desirable, because you want to provide your own signature or remove
the docstring completely to exclude the function from the Sphinx documentation.
The class ``options`` allows you to selectively suppress auto-generated signatures:

.. code-block:: cpp

    PYBIND11_MODULE(example, m) {
        py::options options;
        options.disable_function_signatures();

        m.def("add", [](int a, int b) { return a + b; }, "A function which adds two numbers");
    }

Note that changes to the settings affect only function bindings created during the
lifetime of the ``options`` instance. When it goes out of scope at the end of the module's init function,
the default settings are restored to prevent unwanted side effects.

.. [#f4] http://www.sphinx-doc.org
.. [#f5] http://github.com/pybind/python_example

.. _avoiding-cpp-types-in-docstrings:

Avoiding C++ types in docstrings
================================

Docstrings are generated at the time of the declaration, e.g. when ``.def(...)`` is called.
At this point parameter and return types should be known to pybind11.
If a custom type is not exposed yet through a ``py::class_`` constructor or a custom type caster,
its C++ type name will be used instead to generate the signature in the docstring:

.. code-block:: text

     |  __init__(...)
     |      __init__(self: example.Foo, arg0: ns::Bar) -> None
                                              ^^^^^^^


This limitation can be circumvented by ensuring that C++ classes are registered with pybind11
before they are used as a parameter or return type of a function:

.. code-block:: cpp

    PYBIND11_MODULE(example, m) {

        auto pyFoo = py::class_<ns::Foo>(m, "Foo");
        auto pyBar = py::class_<ns::Bar>(m, "Bar");

        pyFoo.def(py::init<const ns::Bar&>());
        pyBar.def(py::init<const ns::Foo&>());
    }
Smart pointers
##############

std::unique_ptr
===============

Given a class ``Example`` with Python bindings, it's possible to return
instances wrapped in C++11 unique pointers, like so

.. code-block:: cpp

    std::unique_ptr<Example> create_example() { return std::unique_ptr<Example>(new Example()); }

.. code-block:: cpp

    m.def("create_example", &create_example);

In other words, there is nothing special that needs to be done. While returning
unique pointers in this way is allowed, it is *illegal* to use them as function
arguments. For instance, the following function signature cannot be processed
by pybind11.

.. code-block:: cpp

    void do_something_with_example(std::unique_ptr<Example> ex) { ... }

The above signature would imply that Python needs to give up ownership of an
object that is passed to this function, which is generally not possible (for
instance, the object might be referenced elsewhere).

std::shared_ptr
===============

The binding generator for classes, :class:`class_`, can be passed a template
type that denotes a special *holder* type that is used to manage references to
the object.  If no such holder type template argument is given, the default for
a type named ``Type`` is ``std::unique_ptr<Type>``, which means that the object
is deallocated when Python's reference count goes to zero.

It is possible to switch to other types of reference counting wrappers or smart
pointers, which is useful in codebases that rely on them. For instance, the
following snippet causes ``std::shared_ptr`` to be used instead.

.. code-block:: cpp

    py::class_<Example, std::shared_ptr<Example> /* <- holder type */> obj(m, "Example");

Note that any particular class can only be associated with a single holder type.

One potential stumbling block when using holder types is that they need to be
applied consistently. Can you guess what's broken about the following binding
code?

.. code-block:: cpp

    class Child { };

    class Parent {
    public:
       Parent() : child(std::make_shared<Child>()) { }
       Child *get_child() { return child.get(); }  /* Hint: ** DON'T DO THIS ** */
    private:
        std::shared_ptr<Child> child;
    };

    PYBIND11_MODULE(example, m) {
        py::class_<Child, std::shared_ptr<Child>>(m, "Child");

        py::class_<Parent, std::shared_ptr<Parent>>(m, "Parent")
           .def(py::init<>())
           .def("get_child", &Parent::get_child);
    }

The following Python code will cause undefined behavior (and likely a
segmentation fault).

.. code-block:: python

   from example import Parent
   print(Parent().get_child())

The problem is that ``Parent::get_child()`` returns a pointer to an instance of
``Child``, but the fact that this instance is already managed by
``std::shared_ptr<...>`` is lost when passing raw pointers. In this case,
pybind11 will create a second independent ``std::shared_ptr<...>`` that also
claims ownership of the pointer. In the end, the object will be freed **twice**
since these shared pointers have no way of knowing about each other.

There are two ways to resolve this issue:

1. For types that are managed by a smart pointer class, never use raw pointers
   in function arguments or return values. In other words: always consistently
   wrap pointers into their designated holder types (such as
   ``std::shared_ptr<...>``). In this case, the signature of ``get_child()``
   should be modified as follows:

.. code-block:: cpp

    std::shared_ptr<Child> get_child() { return child; }

2. Adjust the definition of ``Child`` by specifying
   ``std::enable_shared_from_this<T>`` (see cppreference_ for details) as a
   base class. This adds a small bit of information to ``Child`` that allows
   pybind11 to realize that there is already an existing
   ``std::shared_ptr<...>`` and communicate with it. In this case, the
   declaration of ``Child`` should look as follows:

.. _cppreference: http://en.cppreference.com/w/cpp/memory/enable_shared_from_this

.. code-block:: cpp

    class Child : public std::enable_shared_from_this<Child> { };

.. _smart_pointers:

Custom smart pointers
=====================

pybind11 supports ``std::unique_ptr`` and ``std::shared_ptr`` right out of the
box. For any other custom smart pointer, transparent conversions can be enabled
using a macro invocation similar to the following. It must be declared at the
top namespace level before any binding code:

.. code-block:: cpp

    PYBIND11_DECLARE_HOLDER_TYPE(T, SmartPtr<T>);

The first argument of :func:`PYBIND11_DECLARE_HOLDER_TYPE` should be a
placeholder name that is used as a template parameter of the second argument.
Thus, feel free to use any identifier, but use it consistently on both sides;
also, don't use the name of a type that already exists in your codebase.

The macro also accepts a third optional boolean parameter that is set to false
by default. Specify

.. code-block:: cpp

    PYBIND11_DECLARE_HOLDER_TYPE(T, SmartPtr<T>, true);

if ``SmartPtr<T>`` can always be initialized from a ``T*`` pointer without the
risk of inconsistencies (such as multiple independent ``SmartPtr`` instances
believing that they are the sole owner of the ``T*`` pointer). A common
situation where ``true`` should be passed is when the ``T`` instances use
*intrusive* reference counting.

Please take a look at the :ref:`macro_notes` before using this feature.

By default, pybind11 assumes that your custom smart pointer has a standard
interface, i.e. provides a ``.get()`` member function to access the underlying
raw pointer. If this is not the case, pybind11's ``holder_helper`` must be
specialized:

.. code-block:: cpp

    // Always needed for custom holder types
    PYBIND11_DECLARE_HOLDER_TYPE(T, SmartPtr<T>);

    // Only needed if the type's `.get()` goes by another name
    namespace pybind11 { namespace detail {
        template <typename T>
        struct holder_helper<SmartPtr<T>> { // <-- specialization
            static const T *get(const SmartPtr<T> &p) { return p.getPointer(); }
        };
    }}

The above specialization informs pybind11 that the custom ``SmartPtr`` class
provides ``.get()`` functionality via ``.getPointer()``.

.. seealso::

    The file :file:`tests/test_smart_ptr.cpp` contains a complete example
    that demonstrates how to work with custom reference-counting holder types
    in more detail.
Classes
#######

This section presents advanced binding code for classes and it is assumed
that you are already familiar with the basics from :doc:`/classes`.

.. _overriding_virtuals:

Overriding virtual functions in Python
======================================

Suppose that a C++ class or interface has a virtual function that we'd like to
to override from within Python (we'll focus on the class ``Animal``; ``Dog`` is
given as a specific example of how one would do this with traditional C++
code).

.. code-block:: cpp

    class Animal {
    public:
        virtual ~Animal() { }
        virtual std::string go(int n_times) = 0;
    };

    class Dog : public Animal {
    public:
        std::string go(int n_times) override {
            std::string result;
            for (int i=0; i<n_times; ++i)
                result += "woof! ";
            return result;
        }
    };

Let's also suppose that we are given a plain function which calls the
function ``go()`` on an arbitrary ``Animal`` instance.

.. code-block:: cpp

    std::string call_go(Animal *animal) {
        return animal->go(3);
    }

Normally, the binding code for these classes would look as follows:

.. code-block:: cpp

    PYBIND11_MODULE(example, m) {
        py::class_<Animal>(m, "Animal")
            .def("go", &Animal::go);

        py::class_<Dog, Animal>(m, "Dog")
            .def(py::init<>());

        m.def("call_go", &call_go);
    }

However, these bindings are impossible to extend: ``Animal`` is not
constructible, and we clearly require some kind of "trampoline" that
redirects virtual calls back to Python.

Defining a new type of ``Animal`` from within Python is possible but requires a
helper class that is defined as follows:

.. code-block:: cpp

    class PyAnimal : public Animal {
    public:
        /* Inherit the constructors */
        using Animal::Animal;

        /* Trampoline (need one for each virtual function) */
        std::string go(int n_times) override {
            PYBIND11_OVERRIDE_PURE(
                std::string, /* Return type */
                Animal,      /* Parent class */
                go,          /* Name of function in C++ (must match Python name) */
                n_times      /* Argument(s) */
            );
        }
    };

The macro :c:macro:`PYBIND11_OVERRIDE_PURE` should be used for pure virtual
functions, and :c:macro:`PYBIND11_OVERRIDE` should be used for functions which have
a default implementation.  There are also two alternate macros
:c:macro:`PYBIND11_OVERRIDE_PURE_NAME` and :c:macro:`PYBIND11_OVERRIDE_NAME` which
take a string-valued name argument between the *Parent class* and *Name of the
function* slots, which defines the name of function in Python. This is required
when the C++ and Python versions of the
function have different names, e.g.  ``operator()`` vs ``__call__``.

The binding code also needs a few minor adaptations (highlighted):

.. code-block:: cpp
    :emphasize-lines: 2,3

    PYBIND11_MODULE(example, m) {
        py::class_<Animal, PyAnimal /* <--- trampoline*/>(m, "Animal")
            .def(py::init<>())
            .def("go", &Animal::go);

        py::class_<Dog, Animal>(m, "Dog")
            .def(py::init<>());

        m.def("call_go", &call_go);
    }

Importantly, pybind11 is made aware of the trampoline helper class by
specifying it as an extra template argument to :class:`class_`. (This can also
be combined with other template arguments such as a custom holder type; the
order of template types does not matter).  Following this, we are able to
define a constructor as usual.

Bindings should be made against the actual class, not the trampoline helper class.

.. code-block:: cpp
    :emphasize-lines: 3

    py::class_<Animal, PyAnimal /* <--- trampoline*/>(m, "Animal");
        .def(py::init<>())
        .def("go", &PyAnimal::go); /* <--- THIS IS WRONG, use &Animal::go */

Note, however, that the above is sufficient for allowing python classes to
extend ``Animal``, but not ``Dog``: see :ref:`virtual_and_inheritance` for the
necessary steps required to providing proper overriding support for inherited
classes.

The Python session below shows how to override ``Animal::go`` and invoke it via
a virtual method call.

.. code-block:: pycon

    >>> from example import *
    >>> d = Dog()
    >>> call_go(d)
    u'woof! woof! woof! '
    >>> class Cat(Animal):
    ...     def go(self, n_times):
    ...             return "meow! " * n_times
    ...
    >>> c = Cat()
    >>> call_go(c)
    u'meow! meow! meow! '

If you are defining a custom constructor in a derived Python class, you *must*
ensure that you explicitly call the bound C++ constructor using ``__init__``,
*regardless* of whether it is a default constructor or not. Otherwise, the
memory for the C++ portion of the instance will be left uninitialized, which
will generally leave the C++ instance in an invalid state and cause undefined
behavior if the C++ instance is subsequently used.

.. versionchanged:: 2.6
   The default pybind11 metaclass will throw a ``TypeError`` when it detects
   that ``__init__`` was not called by a derived class.

Here is an example:

.. code-block:: python

    class Dachshund(Dog):
        def __init__(self, name):
            Dog.__init__(self) # Without this, a TypeError is raised.
            self.name = name
        def bark(self):
            return "yap!"

Note that a direct ``__init__`` constructor *should be called*, and ``super()``
should not be used. For simple cases of linear inheritance, ``super()``
may work, but once you begin mixing Python and C++ multiple inheritance,
things will fall apart due to differences between Python's MRO and C++'s
mechanisms.

Please take a look at the :ref:`macro_notes` before using this feature.

.. note::

    When the overridden type returns a reference or pointer to a type that
    pybind11 converts from Python (for example, numeric values, std::string,
    and other built-in value-converting types), there are some limitations to
    be aware of:

    - because in these cases there is no C++ variable to reference (the value
      is stored in the referenced Python variable), pybind11 provides one in
      the PYBIND11_OVERRIDE macros (when needed) with static storage duration.
      Note that this means that invoking the overridden method on *any*
      instance will change the referenced value stored in *all* instances of
      that type.

    - Attempts to modify a non-const reference will not have the desired
      effect: it will change only the static cache variable, but this change
      will not propagate to underlying Python instance, and the change will be
      replaced the next time the override is invoked.

.. warning::

    The :c:macro:`PYBIND11_OVERRIDE` and accompanying macros used to be called
    ``PYBIND11_OVERLOAD`` up until pybind11 v2.5.0, and :func:`get_override`
    used to be called ``get_overload``. This naming was corrected and the older
    macro and function names may soon be deprecated, in order to reduce
    confusion with overloaded functions and methods and ``py::overload_cast``
    (see :ref:`classes`).

.. seealso::

    The file :file:`tests/test_virtual_functions.cpp` contains a complete
    example that demonstrates how to override virtual functions using pybind11
    in more detail.

.. _virtual_and_inheritance:

Combining virtual functions and inheritance
===========================================

When combining virtual methods with inheritance, you need to be sure to provide
an override for each method for which you want to allow overrides from derived
python classes.  For example, suppose we extend the above ``Animal``/``Dog``
example as follows:

.. code-block:: cpp

    class Animal {
    public:
        virtual std::string go(int n_times) = 0;
        virtual std::string name() { return "unknown"; }
    };
    class Dog : public Animal {
    public:
        std::string go(int n_times) override {
            std::string result;
            for (int i=0; i<n_times; ++i)
                result += bark() + " ";
            return result;
        }
        virtual std::string bark() { return "woof!"; }
    };

then the trampoline class for ``Animal`` must, as described in the previous
section, override ``go()`` and ``name()``, but in order to allow python code to
inherit properly from ``Dog``, we also need a trampoline class for ``Dog`` that
overrides both the added ``bark()`` method *and* the ``go()`` and ``name()``
methods inherited from ``Animal`` (even though ``Dog`` doesn't directly
override the ``name()`` method):

.. code-block:: cpp

    class PyAnimal : public Animal {
    public:
        using Animal::Animal; // Inherit constructors
        std::string go(int n_times) override { PYBIND11_OVERRIDE_PURE(std::string, Animal, go, n_times); }
        std::string name() override { PYBIND11_OVERRIDE(std::string, Animal, name, ); }
    };
    class PyDog : public Dog {
    public:
        using Dog::Dog; // Inherit constructors
        std::string go(int n_times) override { PYBIND11_OVERRIDE(std::string, Dog, go, n_times); }
        std::string name() override { PYBIND11_OVERRIDE(std::string, Dog, name, ); }
        std::string bark() override { PYBIND11_OVERRIDE(std::string, Dog, bark, ); }
    };

.. note::

    Note the trailing commas in the ``PYBIND11_OVERIDE`` calls to ``name()``
    and ``bark()``. These are needed to portably implement a trampoline for a
    function that does not take any arguments. For functions that take
    a nonzero number of arguments, the trailing comma must be omitted.

A registered class derived from a pybind11-registered class with virtual
methods requires a similar trampoline class, *even if* it doesn't explicitly
declare or override any virtual methods itself:

.. code-block:: cpp

    class Husky : public Dog {};
    class PyHusky : public Husky {
    public:
        using Husky::Husky; // Inherit constructors
        std::string go(int n_times) override { PYBIND11_OVERRIDE_PURE(std::string, Husky, go, n_times); }
        std::string name() override { PYBIND11_OVERRIDE(std::string, Husky, name, ); }
        std::string bark() override { PYBIND11_OVERRIDE(std::string, Husky, bark, ); }
    };

There is, however, a technique that can be used to avoid this duplication
(which can be especially helpful for a base class with several virtual
methods).  The technique involves using template trampoline classes, as
follows:

.. code-block:: cpp

    template <class AnimalBase = Animal> class PyAnimal : public AnimalBase {
    public:
        using AnimalBase::AnimalBase; // Inherit constructors
        std::string go(int n_times) override { PYBIND11_OVERRIDE_PURE(std::string, AnimalBase, go, n_times); }
        std::string name() override { PYBIND11_OVERRIDE(std::string, AnimalBase, name, ); }
    };
    template <class DogBase = Dog> class PyDog : public PyAnimal<DogBase> {
    public:
        using PyAnimal<DogBase>::PyAnimal; // Inherit constructors
        // Override PyAnimal's pure virtual go() with a non-pure one:
        std::string go(int n_times) override { PYBIND11_OVERRIDE(std::string, DogBase, go, n_times); }
        std::string bark() override { PYBIND11_OVERRIDE(std::string, DogBase, bark, ); }
    };

This technique has the advantage of requiring just one trampoline method to be
declared per virtual method and pure virtual method override.  It does,
however, require the compiler to generate at least as many methods (and
possibly more, if both pure virtual and overridden pure virtual methods are
exposed, as above).

The classes are then registered with pybind11 using:

.. code-block:: cpp

    py::class_<Animal, PyAnimal<>> animal(m, "Animal");
    py::class_<Dog, Animal, PyDog<>> dog(m, "Dog");
    py::class_<Husky, Dog, PyDog<Husky>> husky(m, "Husky");
    // ... add animal, dog, husky definitions

Note that ``Husky`` did not require a dedicated trampoline template class at
all, since it neither declares any new virtual methods nor provides any pure
virtual method implementations.

With either the repeated-virtuals or templated trampoline methods in place, you
can now create a python class that inherits from ``Dog``:

.. code-block:: python

    class ShihTzu(Dog):
        def bark(self):
            return "yip!"

.. seealso::

    See the file :file:`tests/test_virtual_functions.cpp` for complete examples
    using both the duplication and templated trampoline approaches.

.. _extended_aliases:

Extended trampoline class functionality
=======================================

.. _extended_class_functionality_forced_trampoline:

Forced trampoline class initialisation
--------------------------------------
The trampoline classes described in the previous sections are, by default, only
initialized when needed.  More specifically, they are initialized when a python
class actually inherits from a registered type (instead of merely creating an
instance of the registered type), or when a registered constructor is only
valid for the trampoline class but not the registered class.  This is primarily
for performance reasons: when the trampoline class is not needed for anything
except virtual method dispatching, not initializing the trampoline class
improves performance by avoiding needing to do a run-time check to see if the
inheriting python instance has an overridden method.

Sometimes, however, it is useful to always initialize a trampoline class as an
intermediate class that does more than just handle virtual method dispatching.
For example, such a class might perform extra class initialization, extra
destruction operations, and might define new members and methods to enable a
more python-like interface to a class.

In order to tell pybind11 that it should *always* initialize the trampoline
class when creating new instances of a type, the class constructors should be
declared using ``py::init_alias<Args, ...>()`` instead of the usual
``py::init<Args, ...>()``.  This forces construction via the trampoline class,
ensuring member initialization and (eventual) destruction.

.. seealso::

    See the file :file:`tests/test_virtual_functions.cpp` for complete examples
    showing both normal and forced trampoline instantiation.

Different method signatures
---------------------------
The macro's introduced in :ref:`overriding_virtuals` cover most of the standard
use cases when exposing C++ classes to Python. Sometimes it is hard or unwieldy
to create a direct one-on-one mapping between the arguments and method return
type.

An example would be when the C++ signature contains output arguments using
references (See also :ref:`faq_reference_arguments`). Another way of solving
this is to use the method body of the trampoline class to do conversions to the
input and return of the Python method.

The main building block to do so is the :func:`get_override`, this function
allows retrieving a method implemented in Python from within the trampoline's
methods. Consider for example a C++ method which has the signature
``bool myMethod(int32_t& value)``, where the return indicates whether
something should be done with the ``value``. This can be made convenient on the
Python side by allowing the Python function to return ``None`` or an ``int``:

.. code-block:: cpp

    bool MyClass::myMethod(int32_t& value)
    {
        pybind11::gil_scoped_acquire gil;  // Acquire the GIL while in this scope.
        // Try to look up the overridden method on the Python side.
        pybind11::function override = pybind11::get_override(this, "myMethod");
        if (override) {  // method is found
            auto obj = override(value);  // Call the Python function.
            if (py::isinstance<py::int_>(obj)) {  // check if it returned a Python integer type
                value = obj.cast<int32_t>();  // Cast it and assign it to the value.
                return true;  // Return true; value should be used.
            } else {
                return false;  // Python returned none, return false.
            }
        }
        return false;  // Alternatively return MyClass::myMethod(value);
    }


.. _custom_constructors:

Custom constructors
===================

The syntax for binding constructors was previously introduced, but it only
works when a constructor of the appropriate arguments actually exists on the
C++ side.  To extend this to more general cases, pybind11 makes it possible
to bind factory functions as constructors. For example, suppose you have a
class like this:

.. code-block:: cpp

    class Example {
    private:
        Example(int); // private constructor
    public:
        // Factory function:
        static Example create(int a) { return Example(a); }
    };

    py::class_<Example>(m, "Example")
        .def(py::init(&Example::create));

While it is possible to create a straightforward binding of the static
``create`` method, it may sometimes be preferable to expose it as a constructor
on the Python side. This can be accomplished by calling ``.def(py::init(...))``
with the function reference returning the new instance passed as an argument.
It is also possible to use this approach to bind a function returning a new
instance by raw pointer or by the holder (e.g. ``std::unique_ptr``).

The following example shows the different approaches:

.. code-block:: cpp

    class Example {
    private:
        Example(int); // private constructor
    public:
        // Factory function - returned by value:
        static Example create(int a) { return Example(a); }

        // These constructors are publicly callable:
        Example(double);
        Example(int, int);
        Example(std::string);
    };

    py::class_<Example>(m, "Example")
        // Bind the factory function as a constructor:
        .def(py::init(&Example::create))
        // Bind a lambda function returning a pointer wrapped in a holder:
        .def(py::init([](std::string arg) {
            return std::unique_ptr<Example>(new Example(arg));
        }))
        // Return a raw pointer:
        .def(py::init([](int a, int b) { return new Example(a, b); }))
        // You can mix the above with regular C++ constructor bindings as well:
        .def(py::init<double>())
        ;

When the constructor is invoked from Python, pybind11 will call the factory
function and store the resulting C++ instance in the Python instance.

When combining factory functions constructors with :ref:`virtual function
trampolines <overriding_virtuals>` there are two approaches.  The first is to
add a constructor to the alias class that takes a base value by
rvalue-reference.  If such a constructor is available, it will be used to
construct an alias instance from the value returned by the factory function.
The second option is to provide two factory functions to ``py::init()``: the
first will be invoked when no alias class is required (i.e. when the class is
being used but not inherited from in Python), and the second will be invoked
when an alias is required.

You can also specify a single factory function that always returns an alias
instance: this will result in behaviour similar to ``py::init_alias<...>()``,
as described in the :ref:`extended trampoline class documentation
<extended_aliases>`.

The following example shows the different factory approaches for a class with
an alias:

.. code-block:: cpp

    #include <pybind11/factory.h>
    class Example {
    public:
        // ...
        virtual ~Example() = default;
    };
    class PyExample : public Example {
    public:
        using Example::Example;
        PyExample(Example &&base) : Example(std::move(base)) {}
    };
    py::class_<Example, PyExample>(m, "Example")
        // Returns an Example pointer.  If a PyExample is needed, the Example
        // instance will be moved via the extra constructor in PyExample, above.
        .def(py::init([]() { return new Example(); }))
        // Two callbacks:
        .def(py::init([]() { return new Example(); } /* no alias needed */,
                      []() { return new PyExample(); } /* alias needed */))
        // *Always* returns an alias instance (like py::init_alias<>())
        .def(py::init([]() { return new PyExample(); }))
        ;

Brace initialization
--------------------

``pybind11::init<>`` internally uses C++11 brace initialization to call the
constructor of the target class. This means that it can be used to bind
*implicit* constructors as well:

.. code-block:: cpp

    struct Aggregate {
        int a;
        std::string b;
    };

    py::class_<Aggregate>(m, "Aggregate")
        .def(py::init<int, const std::string &>());

.. note::

    Note that brace initialization preferentially invokes constructor overloads
    taking a ``std::initializer_list``. In the rare event that this causes an
    issue, you can work around it by using ``py::init(...)`` with a lambda
    function that constructs the new object as desired.

.. _classes_with_non_public_destructors:

Non-public destructors
======================

If a class has a private or protected destructor (as might e.g. be the case in
a singleton pattern), a compile error will occur when creating bindings via
pybind11. The underlying issue is that the ``std::unique_ptr`` holder type that
is responsible for managing the lifetime of instances will reference the
destructor even if no deallocations ever take place. In order to expose classes
with private or protected destructors, it is possible to override the holder
type via a holder type argument to ``class_``. Pybind11 provides a helper class
``py::nodelete`` that disables any destructor invocations. In this case, it is
crucial that instances are deallocated on the C++ side to avoid memory leaks.

.. code-block:: cpp

    /* ... definition ... */

    class MyClass {
    private:
        ~MyClass() { }
    };

    /* ... binding code ... */

    py::class_<MyClass, std::unique_ptr<MyClass, py::nodelete>>(m, "MyClass")
        .def(py::init<>())

.. _destructors_that_call_python:

Destructors that call Python
============================

If a Python function is invoked from a C++ destructor, an exception may be thrown
of type :class:`error_already_set`. If this error is thrown out of a class destructor,
``std::terminate()`` will be called, terminating the process. Class destructors
must catch all exceptions of type :class:`error_already_set` to discard the Python
exception using :func:`error_already_set::discard_as_unraisable`.

Every Python function should be treated as *possibly throwing*. When a Python generator
stops yielding items, Python will throw a ``StopIteration`` exception, which can pass
though C++ destructors if the generator's stack frame holds the last reference to C++
objects.

For more information, see :ref:`the documentation on exceptions <unraisable_exceptions>`.

.. code-block:: cpp

    class MyClass {
    public:
        ~MyClass() {
            try {
                py::print("Even printing is dangerous in a destructor");
                py::exec("raise ValueError('This is an unraisable exception')");
            } catch (py::error_already_set &e) {
                // error_context should be information about where/why the occurred,
                // e.g. use __func__ to get the name of the current function
                e.discard_as_unraisable(__func__);
            }
        }
    };

.. note::

    pybind11 does not support C++ destructors marked ``noexcept(false)``.

.. versionadded:: 2.6

.. _implicit_conversions:

Implicit conversions
====================

Suppose that instances of two types ``A`` and ``B`` are used in a project, and
that an ``A`` can easily be converted into an instance of type ``B`` (examples of this
could be a fixed and an arbitrary precision number type).

.. code-block:: cpp

    py::class_<A>(m, "A")
        /// ... members ...

    py::class_<B>(m, "B")
        .def(py::init<A>())
        /// ... members ...

    m.def("func",
        [](const B &) { /* .... */ }
    );

To invoke the function ``func`` using a variable ``a`` containing an ``A``
instance, we'd have to write ``func(B(a))`` in Python. On the other hand, C++
will automatically apply an implicit type conversion, which makes it possible
to directly write ``func(a)``.

In this situation (i.e. where ``B`` has a constructor that converts from
``A``), the following statement enables similar implicit conversions on the
Python side:

.. code-block:: cpp

    py::implicitly_convertible<A, B>();

.. note::

    Implicit conversions from ``A`` to ``B`` only work when ``B`` is a custom
    data type that is exposed to Python via pybind11.

    To prevent runaway recursion, implicit conversions are non-reentrant: an
    implicit conversion invoked as part of another implicit conversion of the
    same type (i.e. from ``A`` to ``B``) will fail.

.. _static_properties:

Static properties
=================

The section on :ref:`properties` discussed the creation of instance properties
that are implemented in terms of C++ getters and setters.

Static properties can also be created in a similar way to expose getters and
setters of static class attributes. Note that the implicit ``self`` argument
also exists in this case and is used to pass the Python ``type`` subclass
instance. This parameter will often not be needed by the C++ side, and the
following example illustrates how to instantiate a lambda getter function
that ignores it:

.. code-block:: cpp

    py::class_<Foo>(m, "Foo")
        .def_property_readonly_static("foo", [](py::object /* self */) { return Foo(); });

Operator overloading
====================

Suppose that we're given the following ``Vector2`` class with a vector addition
and scalar multiplication operation, all implemented using overloaded operators
in C++.

.. code-block:: cpp

    class Vector2 {
    public:
        Vector2(float x, float y) : x(x), y(y) { }

        Vector2 operator+(const Vector2 &v) const { return Vector2(x + v.x, y + v.y); }
        Vector2 operator*(float value) const { return Vector2(x * value, y * value); }
        Vector2& operator+=(const Vector2 &v) { x += v.x; y += v.y; return *this; }
        Vector2& operator*=(float v) { x *= v; y *= v; return *this; }

        friend Vector2 operator*(float f, const Vector2 &v) {
            return Vector2(f * v.x, f * v.y);
        }

        std::string toString() const {
            return "[" + std::to_string(x) + ", " + std::to_string(y) + "]";
        }
    private:
        float x, y;
    };

The following snippet shows how the above operators can be conveniently exposed
to Python.

.. code-block:: cpp

    #include <pybind11/operators.h>

    PYBIND11_MODULE(example, m) {
        py::class_<Vector2>(m, "Vector2")
            .def(py::init<float, float>())
            .def(py::self + py::self)
            .def(py::self += py::self)
            .def(py::self *= float())
            .def(float() * py::self)
            .def(py::self * float())
            .def(-py::self)
            .def("__repr__", &Vector2::toString);
    }

Note that a line like

.. code-block:: cpp

            .def(py::self * float())

is really just short hand notation for

.. code-block:: cpp

    .def("__mul__", [](const Vector2 &a, float b) {
        return a * b;
    }, py::is_operator())

This can be useful for exposing additional operators that don't exist on the
C++ side, or to perform other types of customization. The ``py::is_operator``
flag marker is needed to inform pybind11 that this is an operator, which
returns ``NotImplemented`` when invoked with incompatible arguments rather than
throwing a type error.

.. note::

    To use the more convenient ``py::self`` notation, the additional
    header file :file:`pybind11/operators.h` must be included.

.. seealso::

    The file :file:`tests/test_operator_overloading.cpp` contains a
    complete example that demonstrates how to work with overloaded operators in
    more detail.

.. _pickling:

Pickling support
================

Python's ``pickle`` module provides a powerful facility to serialize and
de-serialize a Python object graph into a binary data stream. To pickle and
unpickle C++ classes using pybind11, a ``py::pickle()`` definition must be
provided. Suppose the class in question has the following signature:

.. code-block:: cpp

    class Pickleable {
    public:
        Pickleable(const std::string &value) : m_value(value) { }
        const std::string &value() const { return m_value; }

        void setExtra(int extra) { m_extra = extra; }
        int extra() const { return m_extra; }
    private:
        std::string m_value;
        int m_extra = 0;
    };

Pickling support in Python is enabled by defining the ``__setstate__`` and
``__getstate__`` methods [#f3]_. For pybind11 classes, use ``py::pickle()``
to bind these two functions:

.. code-block:: cpp

    py::class_<Pickleable>(m, "Pickleable")
        .def(py::init<std::string>())
        .def("value", &Pickleable::value)
        .def("extra", &Pickleable::extra)
        .def("setExtra", &Pickleable::setExtra)
        .def(py::pickle(
            [](const Pickleable &p) { // __getstate__
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(p.value(), p.extra());
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state!");

                /* Create a new C++ instance */
                Pickleable p(t[0].cast<std::string>());

                /* Assign any additional state */
                p.setExtra(t[1].cast<int>());

                return p;
            }
        ));

The ``__setstate__`` part of the ``py::picke()`` definition follows the same
rules as the single-argument version of ``py::init()``. The return type can be
a value, pointer or holder type. See :ref:`custom_constructors` for details.

An instance can now be pickled as follows:

.. code-block:: python

    try:
        import cPickle as pickle  # Use cPickle on Python 2.7
    except ImportError:
        import pickle

    p = Pickleable("test_value")
    p.setExtra(15)
    data = pickle.dumps(p, 2)


.. note::
    Note that only the cPickle module is supported on Python 2.7.

    The second argument to ``dumps`` is also crucial: it selects the pickle
    protocol version 2, since the older version 1 is not supported. Newer
    versions are also fine—for instance, specify ``-1`` to always use the
    latest available version. Beware: failure to follow these instructions
    will cause important pybind11 memory allocation routines to be skipped
    during unpickling, which will likely lead to memory corruption and/or
    segmentation faults.

.. seealso::

    The file :file:`tests/test_pickling.cpp` contains a complete example
    that demonstrates how to pickle and unpickle types using pybind11 in more
    detail.

.. [#f3] http://docs.python.org/3/library/pickle.html#pickling-class-instances

Deepcopy support
================

Python normally uses references in assignments. Sometimes a real copy is needed
to prevent changing all copies. The ``copy`` module [#f5]_ provides these
capabilities.

On Python 3, a class with pickle support is automatically also (deep)copy
compatible. However, performance can be improved by adding custom
``__copy__`` and ``__deepcopy__`` methods. With Python 2.7, these custom methods
are mandatory for (deep)copy compatibility, because pybind11 only supports
cPickle.

For simple classes (deep)copy can be enabled by using the copy constructor,
which should look as follows:

.. code-block:: cpp

    py::class_<Copyable>(m, "Copyable")
        .def("__copy__",  [](const Copyable &self) {
            return Copyable(self);
        })
        .def("__deepcopy__", [](const Copyable &self, py::dict) {
            return Copyable(self);
        }, "memo"_a);

.. note::

    Dynamic attributes will not be copied in this example.

.. [#f5] https://docs.python.org/3/library/copy.html

Multiple Inheritance
====================

pybind11 can create bindings for types that derive from multiple base types
(aka. *multiple inheritance*). To do so, specify all bases in the template
arguments of the ``class_`` declaration:

.. code-block:: cpp

    py::class_<MyType, BaseType1, BaseType2, BaseType3>(m, "MyType")
       ...

The base types can be specified in arbitrary order, and they can even be
interspersed with alias types and holder types (discussed earlier in this
document)---pybind11 will automatically find out which is which. The only
requirement is that the first template argument is the type to be declared.

It is also permitted to inherit multiply from exported C++ classes in Python,
as well as inheriting from multiple Python and/or pybind11-exported classes.

There is one caveat regarding the implementation of this feature:

When only one base type is specified for a C++ type that actually has multiple
bases, pybind11 will assume that it does not participate in multiple
inheritance, which can lead to undefined behavior. In such cases, add the tag
``multiple_inheritance`` to the class constructor:

.. code-block:: cpp

    py::class_<MyType, BaseType2>(m, "MyType", py::multiple_inheritance());

The tag is redundant and does not need to be specified when multiple base types
are listed.

.. _module_local:

Module-local class bindings
===========================

When creating a binding for a class, pybind11 by default makes that binding
"global" across modules.  What this means is that a type defined in one module
can be returned from any module resulting in the same Python type.  For
example, this allows the following:

.. code-block:: cpp

    // In the module1.cpp binding code for module1:
    py::class_<Pet>(m, "Pet")
        .def(py::init<std::string>())
        .def_readonly("name", &Pet::name);

.. code-block:: cpp

    // In the module2.cpp binding code for module2:
    m.def("create_pet", [](std::string name) { return new Pet(name); });

.. code-block:: pycon

    >>> from module1 import Pet
    >>> from module2 import create_pet
    >>> pet1 = Pet("Kitty")
    >>> pet2 = create_pet("Doggy")
    >>> pet2.name()
    'Doggy'

When writing binding code for a library, this is usually desirable: this
allows, for example, splitting up a complex library into multiple Python
modules.

In some cases, however, this can cause conflicts.  For example, suppose two
unrelated modules make use of an external C++ library and each provide custom
bindings for one of that library's classes.  This will result in an error when
a Python program attempts to import both modules (directly or indirectly)
because of conflicting definitions on the external type:

.. code-block:: cpp

    // dogs.cpp

    // Binding for external library class:
    py::class<pets::Pet>(m, "Pet")
        .def("name", &pets::Pet::name);

    // Binding for local extension class:
    py::class<Dog, pets::Pet>(m, "Dog")
        .def(py::init<std::string>());

.. code-block:: cpp

    // cats.cpp, in a completely separate project from the above dogs.cpp.

    // Binding for external library class:
    py::class<pets::Pet>(m, "Pet")
        .def("get_name", &pets::Pet::name);

    // Binding for local extending class:
    py::class<Cat, pets::Pet>(m, "Cat")
        .def(py::init<std::string>());

.. code-block:: pycon

    >>> import cats
    >>> import dogs
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: generic_type: type "Pet" is already registered!

To get around this, you can tell pybind11 to keep the external class binding
localized to the module by passing the ``py::module_local()`` attribute into
the ``py::class_`` constructor:

.. code-block:: cpp

    // Pet binding in dogs.cpp:
    py::class<pets::Pet>(m, "Pet", py::module_local())
        .def("name", &pets::Pet::name);

.. code-block:: cpp

    // Pet binding in cats.cpp:
    py::class<pets::Pet>(m, "Pet", py::module_local())
        .def("get_name", &pets::Pet::name);

This makes the Python-side ``dogs.Pet`` and ``cats.Pet`` into distinct classes,
avoiding the conflict and allowing both modules to be loaded.  C++ code in the
``dogs`` module that casts or returns a ``Pet`` instance will result in a
``dogs.Pet`` Python instance, while C++ code in the ``cats`` module will result
in a ``cats.Pet`` Python instance.

This does come with two caveats, however: First, external modules cannot return
or cast a ``Pet`` instance to Python (unless they also provide their own local
bindings).  Second, from the Python point of view they are two distinct classes.

Note that the locality only applies in the C++ -> Python direction.  When
passing such a ``py::module_local`` type into a C++ function, the module-local
classes are still considered.  This means that if the following function is
added to any module (including but not limited to the ``cats`` and ``dogs``
modules above) it will be callable with either a ``dogs.Pet`` or ``cats.Pet``
argument:

.. code-block:: cpp

    m.def("pet_name", [](const pets::Pet &pet) { return pet.name(); });

For example, suppose the above function is added to each of ``cats.cpp``,
``dogs.cpp`` and ``frogs.cpp`` (where ``frogs.cpp`` is some other module that
does *not* bind ``Pets`` at all).

.. code-block:: pycon

    >>> import cats, dogs, frogs  # No error because of the added py::module_local()
    >>> mycat, mydog = cats.Cat("Fluffy"), dogs.Dog("Rover")
    >>> (cats.pet_name(mycat), dogs.pet_name(mydog))
    ('Fluffy', 'Rover')
    >>> (cats.pet_name(mydog), dogs.pet_name(mycat), frogs.pet_name(mycat))
    ('Rover', 'Fluffy', 'Fluffy')

It is possible to use ``py::module_local()`` registrations in one module even
if another module registers the same type globally: within the module with the
module-local definition, all C++ instances will be cast to the associated bound
Python type.  In other modules any such values are converted to the global
Python type created elsewhere.

.. note::

    STL bindings (as provided via the optional :file:`pybind11/stl_bind.h`
    header) apply ``py::module_local`` by default when the bound type might
    conflict with other modules; see :ref:`stl_bind` for details.

.. note::

    The localization of the bound types is actually tied to the shared object
    or binary generated by the compiler/linker.  For typical modules created
    with ``PYBIND11_MODULE()``, this distinction is not significant.  It is
    possible, however, when :ref:`embedding` to embed multiple modules in the
    same binary (see :ref:`embedding_modules`).  In such a case, the
    localization will apply across all embedded modules within the same binary.

.. seealso::

    The file :file:`tests/test_local_bindings.cpp` contains additional examples
    that demonstrate how ``py::module_local()`` works.

Binding protected member functions
==================================

It's normally not possible to expose ``protected`` member functions to Python:

.. code-block:: cpp

    class A {
    protected:
        int foo() const { return 42; }
    };

    py::class_<A>(m, "A")
        .def("foo", &A::foo); // error: 'foo' is a protected member of 'A'

On one hand, this is good because non-``public`` members aren't meant to be
accessed from the outside. But we may want to make use of ``protected``
functions in derived Python classes.

The following pattern makes this possible:

.. code-block:: cpp

    class A {
    protected:
        int foo() const { return 42; }
    };

    class Publicist : public A { // helper type for exposing protected functions
    public:
        using A::foo; // inherited with different access modifier
    };

    py::class_<A>(m, "A") // bind the primary class
        .def("foo", &Publicist::foo); // expose protected methods via the publicist

This works because ``&Publicist::foo`` is exactly the same function as
``&A::foo`` (same signature and address), just with a different access
modifier. The only purpose of the ``Publicist`` helper class is to make
the function name ``public``.

If the intent is to expose ``protected`` ``virtual`` functions which can be
overridden in Python, the publicist pattern can be combined with the previously
described trampoline:

.. code-block:: cpp

    class A {
    public:
        virtual ~A() = default;

    protected:
        virtual int foo() const { return 42; }
    };

    class Trampoline : public A {
    public:
        int foo() const override { PYBIND11_OVERRIDE(int, A, foo, ); }
    };

    class Publicist : public A {
    public:
        using A::foo;
    };

    py::class_<A, Trampoline>(m, "A") // <-- `Trampoline` here
        .def("foo", &Publicist::foo); // <-- `Publicist` here, not `Trampoline`!

.. note::

    MSVC 2015 has a compiler bug (fixed in version 2017) which
    requires a more explicit function binding in the form of
    ``.def("foo", static_cast<int (A::*)() const>(&Publicist::foo));``
    where ``int (A::*)() const`` is the type of ``A::foo``.

Binding final classes
=====================

Some classes may not be appropriate to inherit from. In C++11, classes can
use the ``final`` specifier to ensure that a class cannot be inherited from.
The ``py::is_final`` attribute can be used to ensure that Python classes
cannot inherit from a specified type. The underlying C++ type does not need
to be declared final.

.. code-block:: cpp

    class IsFinal final {};

    py::class_<IsFinal>(m, "IsFinal", py::is_final());

When you try to inherit from such a class in Python, you will now get this
error:

.. code-block:: pycon

    >>> class PyFinalChild(IsFinal):
    ...     pass
    TypeError: type 'IsFinal' is not an acceptable base type

.. note:: This attribute is currently ignored on PyPy

.. versionadded:: 2.6

Custom automatic downcasters
============================

As explained in :ref:`inheritance`, pybind11 comes with built-in
understanding of the dynamic type of polymorphic objects in C++; that
is, returning a Pet to Python produces a Python object that knows it's
wrapping a Dog, if Pet has virtual methods and pybind11 knows about
Dog and this Pet is in fact a Dog. Sometimes, you might want to
provide this automatic downcasting behavior when creating bindings for
a class hierarchy that does not use standard C++ polymorphism, such as
LLVM [#f4]_. As long as there's some way to determine at runtime
whether a downcast is safe, you can proceed by specializing the
``pybind11::polymorphic_type_hook`` template:

.. code-block:: cpp

    enum class PetKind { Cat, Dog, Zebra };
    struct Pet {   // Not polymorphic: has no virtual methods
        const PetKind kind;
        int age = 0;
      protected:
        Pet(PetKind _kind) : kind(_kind) {}
    };
    struct Dog : Pet {
        Dog() : Pet(PetKind::Dog) {}
        std::string sound = "woof!";
        std::string bark() const { return sound; }
    };

    namespace pybind11 {
        template<> struct polymorphic_type_hook<Pet> {
            static const void *get(const Pet *src, const std::type_info*& type) {
                // note that src may be nullptr
                if (src && src->kind == PetKind::Dog) {
                    type = &typeid(Dog);
                    return static_cast<const Dog*>(src);
                }
                return src;
            }
        };
    } // namespace pybind11

When pybind11 wants to convert a C++ pointer of type ``Base*`` to a
Python object, it calls ``polymorphic_type_hook<Base>::get()`` to
determine if a downcast is possible. The ``get()`` function should use
whatever runtime information is available to determine if its ``src``
parameter is in fact an instance of some class ``Derived`` that
inherits from ``Base``. If it finds such a ``Derived``, it sets ``type
= &typeid(Derived)`` and returns a pointer to the ``Derived`` object
that contains ``src``. Otherwise, it just returns ``src``, leaving
``type`` at its default value of nullptr. If you set ``type`` to a
type that pybind11 doesn't know about, no downcasting will occur, and
the original ``src`` pointer will be used with its static type
``Base*``.

It is critical that the returned pointer and ``type`` argument of
``get()`` agree with each other: if ``type`` is set to something
non-null, the returned pointer must point to the start of an object
whose type is ``type``. If the hierarchy being exposed uses only
single inheritance, a simple ``return src;`` will achieve this just
fine, but in the general case, you must cast ``src`` to the
appropriate derived-class pointer (e.g. using
``static_cast<Derived>(src)``) before allowing it to be returned as a
``void*``.

.. [#f4] https://llvm.org/docs/HowToSetUpLLVMStyleRTTI.html

.. note::

    pybind11's standard support for downcasting objects whose types
    have virtual methods is implemented using
    ``polymorphic_type_hook`` too, using the standard C++ ability to
    determine the most-derived type of a polymorphic object using
    ``typeid()`` and to cast a base pointer to that most-derived type
    (even if you don't know what it is) using ``dynamic_cast<void*>``.

.. seealso::

    The file :file:`tests/test_tagbased_polymorphic.cpp` contains a
    more complete example, including a demonstration of how to provide
    automatic downcasting for an entire class hierarchy without
    writing one get() function for each class.

Accessing the type object
=========================

You can get the type object from a C++ class that has already been registered using:

.. code-block:: python

    py::type T_py = py::type::of<T>();

You can directly use ``py::type::of(ob)`` to get the type object from any python
object, just like ``type(ob)`` in Python.

.. note::

    Other types, like ``py::type::of<int>()``, do not work, see :ref:`type-conversions`.

.. versionadded:: 2.6
Exceptions
##########

Built-in C++ to Python exception translation
============================================

When Python calls C++ code through pybind11, pybind11 provides a C++ exception handler
that will trap C++ exceptions, translate them to the corresponding Python exception,
and raise them so that Python code can handle them.

pybind11 defines translations for ``std::exception`` and its standard
subclasses, and several special exception classes that translate to specific
Python exceptions. Note that these are not actually Python exceptions, so they
cannot be examined using the Python C API. Instead, they are pure C++ objects
that pybind11 will translate the corresponding Python exception when they arrive
at its exception handler.

.. tabularcolumns:: |p{0.5\textwidth}|p{0.45\textwidth}|

+--------------------------------------+--------------------------------------+
|  Exception thrown by C++             |  Translated to Python exception type |
+======================================+======================================+
| :class:`std::exception`              | ``RuntimeError``                     |
+--------------------------------------+--------------------------------------+
| :class:`std::bad_alloc`              | ``MemoryError``                      |
+--------------------------------------+--------------------------------------+
| :class:`std::domain_error`           | ``ValueError``                       |
+--------------------------------------+--------------------------------------+
| :class:`std::invalid_argument`       | ``ValueError``                       |
+--------------------------------------+--------------------------------------+
| :class:`std::length_error`           | ``ValueError``                       |
+--------------------------------------+--------------------------------------+
| :class:`std::out_of_range`           | ``IndexError``                       |
+--------------------------------------+--------------------------------------+
| :class:`std::range_error`            | ``ValueError``                       |
+--------------------------------------+--------------------------------------+
| :class:`std::overflow_error`         | ``OverflowError``                    |
+--------------------------------------+--------------------------------------+
| :class:`pybind11::stop_iteration`    | ``StopIteration`` (used to implement |
|                                      | custom iterators)                    |
+--------------------------------------+--------------------------------------+
| :class:`pybind11::index_error`       | ``IndexError`` (used to indicate out |
|                                      | of bounds access in ``__getitem__``, |
|                                      | ``__setitem__``, etc.)               |
+--------------------------------------+--------------------------------------+
| :class:`pybind11::key_error`         | ``KeyError`` (used to indicate out   |
|                                      | of bounds access in ``__getitem__``, |
|                                      | ``__setitem__`` in dict-like         |
|                                      | objects, etc.)                       |
+--------------------------------------+--------------------------------------+
| :class:`pybind11::value_error`       | ``ValueError`` (used to indicate     |
|                                      | wrong value passed in                |
|                                      | ``container.remove(...)``)           |
+--------------------------------------+--------------------------------------+
| :class:`pybind11::type_error`        | ``TypeError``                        |
+--------------------------------------+--------------------------------------+
| :class:`pybind11::buffer_error`      | ``BufferError``                      |
+--------------------------------------+--------------------------------------+
| :class:`pybind11::import_error`      | ``import_error``                     |
+--------------------------------------+--------------------------------------+
| Any other exception                  | ``RuntimeError``                     |
+--------------------------------------+--------------------------------------+

Exception translation is not bidirectional. That is, *catching* the C++
exceptions defined above above will not trap exceptions that originate from
Python. For that, catch :class:`pybind11::error_already_set`. See :ref:`below
<handling_python_exceptions_cpp>` for further details.

There is also a special exception :class:`cast_error` that is thrown by
:func:`handle::call` when the input arguments cannot be converted to Python
objects.

Registering custom translators
==============================

If the default exception conversion policy described above is insufficient,
pybind11 also provides support for registering custom exception translators.
To register a simple exception conversion that translates a C++ exception into
a new Python exception using the C++ exception's ``what()`` method, a helper
function is available:

.. code-block:: cpp

    py::register_exception<CppExp>(module, "PyExp");

This call creates a Python exception class with the name ``PyExp`` in the given
module and automatically converts any encountered exceptions of type ``CppExp``
into Python exceptions of type ``PyExp``.

It is possible to specify base class for the exception using the third
parameter, a `handle`:

.. code-block:: cpp

    py::register_exception<CppExp>(module, "PyExp", PyExc_RuntimeError);

Then `PyExp` can be caught both as `PyExp` and `RuntimeError`.

The class objects of the built-in Python exceptions are listed in the Python
documentation on `Standard Exceptions <https://docs.python.org/3/c-api/exceptions.html#standard-exceptions>`_.
The default base class is `PyExc_Exception`.

When more advanced exception translation is needed, the function
``py::register_exception_translator(translator)`` can be used to register
functions that can translate arbitrary exception types (and which may include
additional logic to do so).  The function takes a stateless callable (e.g.  a
function pointer or a lambda function without captured variables) with the call
signature ``void(std::exception_ptr)``.

When a C++ exception is thrown, the registered exception translators are tried
in reverse order of registration (i.e. the last registered translator gets the
first shot at handling the exception).

Inside the translator, ``std::rethrow_exception`` should be used within
a try block to re-throw the exception.  One or more catch clauses to catch
the appropriate exceptions should then be used with each clause using
``PyErr_SetString`` to set a Python exception or ``ex(string)`` to set
the python exception to a custom exception type (see below).

To declare a custom Python exception type, declare a ``py::exception`` variable
and use this in the associated exception translator (note: it is often useful
to make this a static declaration when using it inside a lambda expression
without requiring capturing).

The following example demonstrates this for a hypothetical exception classes
``MyCustomException`` and ``OtherException``: the first is translated to a
custom python exception ``MyCustomError``, while the second is translated to a
standard python RuntimeError:

.. code-block:: cpp

    static py::exception<MyCustomException> exc(m, "MyCustomError");
    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const MyCustomException &e) {
            exc(e.what());
        } catch (const OtherException &e) {
            PyErr_SetString(PyExc_RuntimeError, e.what());
        }
    });

Multiple exceptions can be handled by a single translator, as shown in the
example above. If the exception is not caught by the current translator, the
previously registered one gets a chance.

If none of the registered exception translators is able to handle the
exception, it is handled by the default converter as described in the previous
section.

.. seealso::

    The file :file:`tests/test_exceptions.cpp` contains examples
    of various custom exception translators and custom exception types.

.. note::

    Call either ``PyErr_SetString`` or a custom exception's call
    operator (``exc(string)``) for every exception caught in a custom exception
    translator.  Failure to do so will cause Python to crash with ``SystemError:
    error return without exception set``.

    Exceptions that you do not plan to handle should simply not be caught, or
    may be explicitly (re-)thrown to delegate it to the other,
    previously-declared existing exception translators.

.. _handling_python_exceptions_cpp:

Handling exceptions from Python in C++
======================================

When C++ calls Python functions, such as in a callback function or when
manipulating Python objects, and Python raises an ``Exception``, pybind11
converts the Python exception into a C++ exception of type
:class:`pybind11::error_already_set` whose payload contains a C++ string textual
summary and the actual Python exception. ``error_already_set`` is used to
propagate Python exception back to Python (or possibly, handle them in C++).

.. tabularcolumns:: |p{0.5\textwidth}|p{0.45\textwidth}|

+--------------------------------------+--------------------------------------+
|  Exception raised in Python          |  Thrown as C++ exception type        |
+======================================+======================================+
| Any Python ``Exception``             | :class:`pybind11::error_already_set` |
+--------------------------------------+--------------------------------------+

For example:

.. code-block:: cpp

    try {
        // open("missing.txt", "r")
        auto file = py::module_::import("io").attr("open")("missing.txt", "r");
        auto text = file.attr("read")();
        file.attr("close")();
    } catch (py::error_already_set &e) {
        if (e.matches(PyExc_FileNotFoundError)) {
            py::print("missing.txt not found");
        } else if (e.matches(PyExc_PermissionError)) {
            py::print("missing.txt found but not accessible");
        } else {
            throw;
        }
    }

Note that C++ to Python exception translation does not apply here, since that is
a method for translating C++ exceptions to Python, not vice versa. The error raised
from Python is always ``error_already_set``.

This example illustrates this behavior:

.. code-block:: cpp

    try {
        py::eval("raise ValueError('The Ring')");
    } catch (py::value_error &boromir) {
        // Boromir never gets the ring
        assert(false);
    } catch (py::error_already_set &frodo) {
        // Frodo gets the ring
        py::print("I will take the ring");
    }

    try {
        // py::value_error is a request for pybind11 to raise a Python exception
        throw py::value_error("The ball");
    } catch (py::error_already_set &cat) {
        // cat won't catch the ball since
        // py::value_error is not a Python exception
        assert(false);
    } catch (py::value_error &dog) {
        // dog will catch the ball
        py::print("Run Spot run");
        throw;  // Throw it again (pybind11 will raise ValueError)
    }

Handling errors from the Python C API
=====================================

Where possible, use :ref:`pybind11 wrappers <wrappers>` instead of calling
the Python C API directly. When calling the Python C API directly, in
addition to manually managing reference counts, one must follow the pybind11
error protocol, which is outlined here.

After calling the Python C API, if Python returns an error,
``throw py::error_already_set();``, which allows pybind11 to deal with the
exception and pass it back to the Python interpreter. This includes calls to
the error setting functions such as ``PyErr_SetString``.

.. code-block:: cpp

    PyErr_SetString(PyExc_TypeError, "C API type error demo");
    throw py::error_already_set();

    // But it would be easier to simply...
    throw py::type_error("pybind11 wrapper type error");

Alternately, to ignore the error, call `PyErr_Clear
<https://docs.python.org/3/c-api/exceptions.html#c.PyErr_Clear>`_.

Any Python error must be thrown or cleared, or Python/pybind11 will be left in
an invalid state.

.. _unraisable_exceptions:

Handling unraisable exceptions
==============================

If a Python function invoked from a C++ destructor or any function marked
``noexcept(true)`` (collectively, "noexcept functions") throws an exception, there
is no way to propagate the exception, as such functions may not throw.
Should they throw or fail to catch any exceptions in their call graph,
the C++ runtime calls ``std::terminate()`` to abort immediately.

Similarly, Python exceptions raised in a class's ``__del__`` method do not
propagate, but are logged by Python as an unraisable error. In Python 3.8+, a
`system hook is triggered
<https://docs.python.org/3/library/sys.html#sys.unraisablehook>`_
and an auditing event is logged.

Any noexcept function should have a try-catch block that traps
class:`error_already_set` (or any other exception that can occur). Note that
pybind11 wrappers around Python exceptions such as
:class:`pybind11::value_error` are *not* Python exceptions; they are C++
exceptions that pybind11 catches and converts to Python exceptions. Noexcept
functions cannot propagate these exceptions either. A useful approach is to
convert them to Python exceptions and then ``discard_as_unraisable`` as shown
below.

.. code-block:: cpp

    void nonthrowing_func() noexcept(true) {
        try {
            // ...
        } catch (py::error_already_set &eas) {
            // Discard the Python error using Python APIs, using the C++ magic
            // variable __func__. Python already knows the type and value and of the
            // exception object.
            eas.discard_as_unraisable(__func__);
        } catch (const std::exception &e) {
            // Log and discard C++ exceptions.
            third_party::log(e);
        }
    }

.. versionadded:: 2.6
Functions
#########

Before proceeding with this section, make sure that you are already familiar
with the basics of binding functions and classes, as explained in :doc:`/basics`
and :doc:`/classes`. The following guide is applicable to both free and member
functions, i.e. *methods* in Python.

.. _return_value_policies:

Return value policies
=====================

Python and C++ use fundamentally different ways of managing the memory and
lifetime of objects managed by them. This can lead to issues when creating
bindings for functions that return a non-trivial type. Just by looking at the
type information, it is not clear whether Python should take charge of the
returned value and eventually free its resources, or if this is handled on the
C++ side. For this reason, pybind11 provides a several *return value policy*
annotations that can be passed to the :func:`module_::def` and
:func:`class_::def` functions. The default policy is
:enum:`return_value_policy::automatic`.

Return value policies are tricky, and it's very important to get them right.
Just to illustrate what can go wrong, consider the following simple example:

.. code-block:: cpp

    /* Function declaration */
    Data *get_data() { return _data; /* (pointer to a static data structure) */ }
    ...

    /* Binding code */
    m.def("get_data", &get_data); // <-- KABOOM, will cause crash when called from Python

What's going on here? When ``get_data()`` is called from Python, the return
value (a native C++ type) must be wrapped to turn it into a usable Python type.
In this case, the default return value policy (:enum:`return_value_policy::automatic`)
causes pybind11 to assume ownership of the static ``_data`` instance.

When Python's garbage collector eventually deletes the Python
wrapper, pybind11 will also attempt to delete the C++ instance (via ``operator
delete()``) due to the implied ownership. At this point, the entire application
will come crashing down, though errors could also be more subtle and involve
silent data corruption.

In the above example, the policy :enum:`return_value_policy::reference` should have
been specified so that the global data instance is only *referenced* without any
implied transfer of ownership, i.e.:

.. code-block:: cpp

    m.def("get_data", &get_data, return_value_policy::reference);

On the other hand, this is not the right policy for many other situations,
where ignoring ownership could lead to resource leaks.
As a developer using pybind11, it's important to be familiar with the different
return value policies, including which situation calls for which one of them.
The following table provides an overview of available policies:

.. tabularcolumns:: |p{0.5\textwidth}|p{0.45\textwidth}|

+--------------------------------------------------+----------------------------------------------------------------------------+
| Return value policy                              | Description                                                                |
+==================================================+============================================================================+
| :enum:`return_value_policy::take_ownership`      | Reference an existing object (i.e. do not create a new copy) and take      |
|                                                  | ownership. Python will call the destructor and delete operator when the    |
|                                                  | object's reference count reaches zero. Undefined behavior ensues when the  |
|                                                  | C++ side does the same, or when the data was not dynamically allocated.    |
+--------------------------------------------------+----------------------------------------------------------------------------+
| :enum:`return_value_policy::copy`                | Create a new copy of the returned object, which will be owned by Python.   |
|                                                  | This policy is comparably safe because the lifetimes of the two instances  |
|                                                  | are decoupled.                                                             |
+--------------------------------------------------+----------------------------------------------------------------------------+
| :enum:`return_value_policy::move`                | Use ``std::move`` to move the return value contents into a new instance    |
|                                                  | that will be owned by Python. This policy is comparably safe because the   |
|                                                  | lifetimes of the two instances (move source and destination) are decoupled.|
+--------------------------------------------------+----------------------------------------------------------------------------+
| :enum:`return_value_policy::reference`           | Reference an existing object, but do not take ownership. The C++ side is   |
|                                                  | responsible for managing the object's lifetime and deallocating it when    |
|                                                  | it is no longer used. Warning: undefined behavior will ensue when the C++  |
|                                                  | side deletes an object that is still referenced and used by Python.        |
+--------------------------------------------------+----------------------------------------------------------------------------+
| :enum:`return_value_policy::reference_internal`  | Indicates that the lifetime of the return value is tied to the lifetime    |
|                                                  | of a parent object, namely the implicit ``this``, or ``self`` argument of  |
|                                                  | the called method or property. Internally, this policy works just like     |
|                                                  | :enum:`return_value_policy::reference` but additionally applies a          |
|                                                  | ``keep_alive<0, 1>`` *call policy* (described in the next section) that    |
|                                                  | prevents the parent object from being garbage collected as long as the     |
|                                                  | return value is referenced by Python. This is the default policy for       |
|                                                  | property getters created via ``def_property``, ``def_readwrite``, etc.     |
+--------------------------------------------------+----------------------------------------------------------------------------+
| :enum:`return_value_policy::automatic`           | **Default policy.** This policy falls back to the policy                   |
|                                                  | :enum:`return_value_policy::take_ownership` when the return value is a     |
|                                                  | pointer. Otherwise, it uses :enum:`return_value_policy::move` or           |
|                                                  | :enum:`return_value_policy::copy` for rvalue and lvalue references,        |
|                                                  | respectively. See above for a description of what all of these different   |
|                                                  | policies do.                                                               |
+--------------------------------------------------+----------------------------------------------------------------------------+
| :enum:`return_value_policy::automatic_reference` | As above, but use policy :enum:`return_value_policy::reference` when the   |
|                                                  | return value is a pointer. This is the default conversion policy for       |
|                                                  | function arguments when calling Python functions manually from C++ code    |
|                                                  | (i.e. via handle::operator()). You probably won't need to use this.        |
+--------------------------------------------------+----------------------------------------------------------------------------+

Return value policies can also be applied to properties:

.. code-block:: cpp

    class_<MyClass>(m, "MyClass")
        .def_property("data", &MyClass::getData, &MyClass::setData,
                      py::return_value_policy::copy);

Technically, the code above applies the policy to both the getter and the
setter function, however, the setter doesn't really care about *return*
value policies which makes this a convenient terse syntax. Alternatively,
targeted arguments can be passed through the :class:`cpp_function` constructor:

.. code-block:: cpp

    class_<MyClass>(m, "MyClass")
        .def_property("data"
            py::cpp_function(&MyClass::getData, py::return_value_policy::copy),
            py::cpp_function(&MyClass::setData)
        );

.. warning::

    Code with invalid return value policies might access uninitialized memory or
    free data structures multiple times, which can lead to hard-to-debug
    non-determinism and segmentation faults, hence it is worth spending the
    time to understand all the different options in the table above.

.. note::

    One important aspect of the above policies is that they only apply to
    instances which pybind11 has *not* seen before, in which case the policy
    clarifies essential questions about the return value's lifetime and
    ownership.  When pybind11 knows the instance already (as identified by its
    type and address in memory), it will return the existing Python object
    wrapper rather than creating a new copy.

.. note::

    The next section on :ref:`call_policies` discusses *call policies* that can be
    specified *in addition* to a return value policy from the list above. Call
    policies indicate reference relationships that can involve both return values
    and parameters of functions.

.. note::

   As an alternative to elaborate call policies and lifetime management logic,
   consider using smart pointers (see the section on :ref:`smart_pointers` for
   details). Smart pointers can tell whether an object is still referenced from
   C++ or Python, which generally eliminates the kinds of inconsistencies that
   can lead to crashes or undefined behavior. For functions returning smart
   pointers, it is not necessary to specify a return value policy.

.. _call_policies:

Additional call policies
========================

In addition to the above return value policies, further *call policies* can be
specified to indicate dependencies between parameters or ensure a certain state
for the function call.

Keep alive
----------

In general, this policy is required when the C++ object is any kind of container
and another object is being added to the container. ``keep_alive<Nurse, Patient>``
indicates that the argument with index ``Patient`` should be kept alive at least
until the argument with index ``Nurse`` is freed by the garbage collector. Argument
indices start at one, while zero refers to the return value. For methods, index
``1`` refers to the implicit ``this`` pointer, while regular arguments begin at
index ``2``. Arbitrarily many call policies can be specified. When a ``Nurse``
with value ``None`` is detected at runtime, the call policy does nothing.

When the nurse is not a pybind11-registered type, the implementation internally
relies on the ability to create a *weak reference* to the nurse object. When
the nurse object is not a pybind11-registered type and does not support weak
references, an exception will be thrown.

Consider the following example: here, the binding code for a list append
operation ties the lifetime of the newly added element to the underlying
container:

.. code-block:: cpp

    py::class_<List>(m, "List")
        .def("append", &List::append, py::keep_alive<1, 2>());

For consistency, the argument indexing is identical for constructors. Index
``1`` still refers to the implicit ``this`` pointer, i.e. the object which is
being constructed. Index ``0`` refers to the return type which is presumed to
be ``void`` when a constructor is viewed like a function. The following example
ties the lifetime of the constructor element to the constructed object:

.. code-block:: cpp

    py::class_<Nurse>(m, "Nurse")
        .def(py::init<Patient &>(), py::keep_alive<1, 2>());

.. note::

    ``keep_alive`` is analogous to the ``with_custodian_and_ward`` (if Nurse,
    Patient != 0) and ``with_custodian_and_ward_postcall`` (if Nurse/Patient ==
    0) policies from Boost.Python.

Call guard
----------

The ``call_guard<T>`` policy allows any scope guard type ``T`` to be placed
around the function call. For example, this definition:

.. code-block:: cpp

    m.def("foo", foo, py::call_guard<T>());

is equivalent to the following pseudocode:

.. code-block:: cpp

    m.def("foo", [](args...) {
        T scope_guard;
        return foo(args...); // forwarded arguments
    });

The only requirement is that ``T`` is default-constructible, but otherwise any
scope guard will work. This is very useful in combination with `gil_scoped_release`.
See :ref:`gil`.

Multiple guards can also be specified as ``py::call_guard<T1, T2, T3...>``. The
constructor order is left to right and destruction happens in reverse.

.. seealso::

    The file :file:`tests/test_call_policies.cpp` contains a complete example
    that demonstrates using `keep_alive` and `call_guard` in more detail.

.. _python_objects_as_args:

Python objects as arguments
===========================

pybind11 exposes all major Python types using thin C++ wrapper classes. These
wrapper classes can also be used as parameters of functions in bindings, which
makes it possible to directly work with native Python types on the C++ side.
For instance, the following statement iterates over a Python ``dict``:

.. code-block:: cpp

    void print_dict(py::dict dict) {
        /* Easily interact with Python types */
        for (auto item : dict)
            std::cout << "key=" << std::string(py::str(item.first)) << ", "
                      << "value=" << std::string(py::str(item.second)) << std::endl;
    }

It can be exported:

.. code-block:: cpp

    m.def("print_dict", &print_dict);

And used in Python as usual:

.. code-block:: pycon

    >>> print_dict({'foo': 123, 'bar': 'hello'})
    key=foo, value=123
    key=bar, value=hello

For more information on using Python objects in C++, see :doc:`/advanced/pycpp/index`.

Accepting \*args and \*\*kwargs
===============================

Python provides a useful mechanism to define functions that accept arbitrary
numbers of arguments and keyword arguments:

.. code-block:: python

   def generic(*args, **kwargs):
       ...  # do something with args and kwargs

Such functions can also be created using pybind11:

.. code-block:: cpp

   void generic(py::args args, py::kwargs kwargs) {
       /// .. do something with args
       if (kwargs)
           /// .. do something with kwargs
   }

   /// Binding code
   m.def("generic", &generic);

The class ``py::args`` derives from ``py::tuple`` and ``py::kwargs`` derives
from ``py::dict``.

You may also use just one or the other, and may combine these with other
arguments as long as the ``py::args`` and ``py::kwargs`` arguments are the last
arguments accepted by the function.

Please refer to the other examples for details on how to iterate over these,
and on how to cast their entries into C++ objects. A demonstration is also
available in ``tests/test_kwargs_and_defaults.cpp``.

.. note::

    When combining \*args or \*\*kwargs with :ref:`keyword_args` you should
    *not* include ``py::arg`` tags for the ``py::args`` and ``py::kwargs``
    arguments.

Default arguments revisited
===========================

The section on :ref:`default_args` previously discussed basic usage of default
arguments using pybind11. One noteworthy aspect of their implementation is that
default arguments are converted to Python objects right at declaration time.
Consider the following example:

.. code-block:: cpp

    py::class_<MyClass>("MyClass")
        .def("myFunction", py::arg("arg") = SomeType(123));

In this case, pybind11 must already be set up to deal with values of the type
``SomeType`` (via a prior instantiation of ``py::class_<SomeType>``), or an
exception will be thrown.

Another aspect worth highlighting is that the "preview" of the default argument
in the function signature is generated using the object's ``__repr__`` method.
If not available, the signature may not be very helpful, e.g.:

.. code-block:: pycon

    FUNCTIONS
    ...
    |  myFunction(...)
    |      Signature : (MyClass, arg : SomeType = <SomeType object at 0x101b7b080>) -> NoneType
    ...

The first way of addressing this is by defining ``SomeType.__repr__``.
Alternatively, it is possible to specify the human-readable preview of the
default argument manually using the ``arg_v`` notation:

.. code-block:: cpp

    py::class_<MyClass>("MyClass")
        .def("myFunction", py::arg_v("arg", SomeType(123), "SomeType(123)"));

Sometimes it may be necessary to pass a null pointer value as a default
argument. In this case, remember to cast it to the underlying type in question,
like so:

.. code-block:: cpp

    py::class_<MyClass>("MyClass")
        .def("myFunction", py::arg("arg") = static_cast<SomeType *>(nullptr));

Keyword-only arguments
======================

Python 3 introduced keyword-only arguments by specifying an unnamed ``*``
argument in a function definition:

.. code-block:: python

    def f(a, *, b):  # a can be positional or via keyword; b must be via keyword
        pass

    f(a=1, b=2)  # good
    f(b=2, a=1)  # good
    f(1, b=2)    # good
    f(1, 2)      # TypeError: f() takes 1 positional argument but 2 were given

Pybind11 provides a ``py::kw_only`` object that allows you to implement
the same behaviour by specifying the object between positional and keyword-only
argument annotations when registering the function:

.. code-block:: cpp

    m.def("f", [](int a, int b) { /* ... */ },
          py::arg("a"), py::kw_only(), py::arg("b"));

Note that you currently cannot combine this with a ``py::args`` argument.  This
feature does *not* require Python 3 to work.

.. versionadded:: 2.6

Positional-only arguments
=========================

Python 3.8 introduced a new positional-only argument syntax, using ``/`` in the
function definition (note that this has been a convention for CPython
positional arguments, such as in ``pow()``, since Python 2). You can
do the same thing in any version of Python using ``py::pos_only()``:

.. code-block:: cpp

   m.def("f", [](int a, int b) { /* ... */ },
          py::arg("a"), py::pos_only(), py::arg("b"));

You now cannot give argument ``a`` by keyword. This can be combined with
keyword-only arguments, as well.

.. versionadded:: 2.6

.. _nonconverting_arguments:

Non-converting arguments
========================

Certain argument types may support conversion from one type to another.  Some
examples of conversions are:

* :ref:`implicit_conversions` declared using ``py::implicitly_convertible<A,B>()``
* Calling a method accepting a double with an integer argument
* Calling a ``std::complex<float>`` argument with a non-complex python type
  (for example, with a float).  (Requires the optional ``pybind11/complex.h``
  header).
* Calling a function taking an Eigen matrix reference with a numpy array of the
  wrong type or of an incompatible data layout.  (Requires the optional
  ``pybind11/eigen.h`` header).

This behaviour is sometimes undesirable: the binding code may prefer to raise
an error rather than convert the argument.  This behaviour can be obtained
through ``py::arg`` by calling the ``.noconvert()`` method of the ``py::arg``
object, such as:

.. code-block:: cpp

    m.def("floats_only", [](double f) { return 0.5 * f; }, py::arg("f").noconvert());
    m.def("floats_preferred", [](double f) { return 0.5 * f; }, py::arg("f"));

Attempting the call the second function (the one without ``.noconvert()``) with
an integer will succeed, but attempting to call the ``.noconvert()`` version
will fail with a ``TypeError``:

.. code-block:: pycon

    >>> floats_preferred(4)
    2.0
    >>> floats_only(4)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: floats_only(): incompatible function arguments. The following argument types are supported:
        1. (f: float) -> float

    Invoked with: 4

You may, of course, combine this with the :var:`_a` shorthand notation (see
:ref:`keyword_args`) and/or :ref:`default_args`.  It is also permitted to omit
the argument name by using the ``py::arg()`` constructor without an argument
name, i.e. by specifying ``py::arg().noconvert()``.

.. note::

    When specifying ``py::arg`` options it is necessary to provide the same
    number of options as the bound function has arguments.  Thus if you want to
    enable no-convert behaviour for just one of several arguments, you will
    need to specify a ``py::arg()`` annotation for each argument with the
    no-convert argument modified to ``py::arg().noconvert()``.

.. _none_arguments:

Allow/Prohibiting None arguments
================================

When a C++ type registered with :class:`py::class_` is passed as an argument to
a function taking the instance as pointer or shared holder (e.g. ``shared_ptr``
or a custom, copyable holder as described in :ref:`smart_pointers`), pybind
allows ``None`` to be passed from Python which results in calling the C++
function with ``nullptr`` (or an empty holder) for the argument.

To explicitly enable or disable this behaviour, using the
``.none`` method of the :class:`py::arg` object:

.. code-block:: cpp

    py::class_<Dog>(m, "Dog").def(py::init<>());
    py::class_<Cat>(m, "Cat").def(py::init<>());
    m.def("bark", [](Dog *dog) -> std::string {
        if (dog) return "woof!"; /* Called with a Dog instance */
        else return "(no dog)"; /* Called with None, dog == nullptr */
    }, py::arg("dog").none(true));
    m.def("meow", [](Cat *cat) -> std::string {
        // Can't be called with None argument
        return "meow";
    }, py::arg("cat").none(false));

With the above, the Python call ``bark(None)`` will return the string ``"(no
dog)"``, while attempting to call ``meow(None)`` will raise a ``TypeError``:

.. code-block:: pycon

    >>> from animals import Dog, Cat, bark, meow
    >>> bark(Dog())
    'woof!'
    >>> meow(Cat())
    'meow'
    >>> bark(None)
    '(no dog)'
    >>> meow(None)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: meow(): incompatible function arguments. The following argument types are supported:
        1. (cat: animals.Cat) -> str

    Invoked with: None

The default behaviour when the tag is unspecified is to allow ``None``.

.. note::

    Even when ``.none(true)`` is specified for an argument, ``None`` will be converted to a
    ``nullptr`` *only* for custom and :ref:`opaque <opaque>` types. Pointers to built-in types
    (``double *``, ``int *``, ...) and STL types (``std::vector<T> *``, ...; if ``pybind11/stl.h``
    is included) are copied when converted to C++ (see :doc:`/advanced/cast/overview`) and will
    not allow ``None`` as argument.  To pass optional argument of these copied types consider
    using ``std::optional<T>``

.. _overload_resolution:

Overload resolution order
=========================

When a function or method with multiple overloads is called from Python,
pybind11 determines which overload to call in two passes.  The first pass
attempts to call each overload without allowing argument conversion (as if
every argument had been specified as ``py::arg().noconvert()`` as described
above).

If no overload succeeds in the no-conversion first pass, a second pass is
attempted in which argument conversion is allowed (except where prohibited via
an explicit ``py::arg().noconvert()`` attribute in the function definition).

If the second pass also fails a ``TypeError`` is raised.

Within each pass, overloads are tried in the order they were registered with
pybind11. If the ``py::prepend()`` tag is added to the definition, a function
can be placed at the beginning of the overload sequence instead, allowing user
overloads to proceed built in functions.

What this means in practice is that pybind11 will prefer any overload that does
not require conversion of arguments to an overload that does, but otherwise
prefers earlier-defined overloads to later-defined ones.

.. note::

    pybind11 does *not* further prioritize based on the number/pattern of
    overloaded arguments.  That is, pybind11 does not prioritize a function
    requiring one conversion over one requiring three, but only prioritizes
    overloads requiring no conversion at all to overloads that require
    conversion of at least one argument.

.. versionadded:: 2.6

    The ``py::prepend()`` tag.
.. _embedding:

Embedding the interpreter
#########################

While pybind11 is mainly focused on extending Python using C++, it's also
possible to do the reverse: embed the Python interpreter into a C++ program.
All of the other documentation pages still apply here, so refer to them for
general pybind11 usage. This section will cover a few extra things required
for embedding.

Getting started
===============

A basic executable with an embedded interpreter can be created with just a few
lines of CMake and the ``pybind11::embed`` target, as shown below. For more
information, see :doc:`/compiling`.

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.4)
    project(example)

    find_package(pybind11 REQUIRED)  # or `add_subdirectory(pybind11)`

    add_executable(example main.cpp)
    target_link_libraries(example PRIVATE pybind11::embed)

The essential structure of the ``main.cpp`` file looks like this:

.. code-block:: cpp

    #include <pybind11/embed.h> // everything needed for embedding
    namespace py = pybind11;

    int main() {
        py::scoped_interpreter guard{}; // start the interpreter and keep it alive

        py::print("Hello, World!"); // use the Python API
    }

The interpreter must be initialized before using any Python API, which includes
all the functions and classes in pybind11. The RAII guard class `scoped_interpreter`
takes care of the interpreter lifetime. After the guard is destroyed, the interpreter
shuts down and clears its memory. No Python functions can be called after this.

Executing Python code
=====================

There are a few different ways to run Python code. One option is to use `eval`,
`exec` or `eval_file`, as explained in :ref:`eval`. Here is a quick example in
the context of an executable with an embedded interpreter:

.. code-block:: cpp

    #include <pybind11/embed.h>
    namespace py = pybind11;

    int main() {
        py::scoped_interpreter guard{};

        py::exec(R"(
            kwargs = dict(name="World", number=42)
            message = "Hello, {name}! The answer is {number}".format(**kwargs)
            print(message)
        )");
    }

Alternatively, similar results can be achieved using pybind11's API (see
:doc:`/advanced/pycpp/index` for more details).

.. code-block:: cpp

    #include <pybind11/embed.h>
    namespace py = pybind11;
    using namespace py::literals;

    int main() {
        py::scoped_interpreter guard{};

        auto kwargs = py::dict("name"_a="World", "number"_a=42);
        auto message = "Hello, {name}! The answer is {number}"_s.format(**kwargs);
        py::print(message);
    }

The two approaches can also be combined:

.. code-block:: cpp

    #include <pybind11/embed.h>
    #include <iostream>

    namespace py = pybind11;
    using namespace py::literals;

    int main() {
        py::scoped_interpreter guard{};

        auto locals = py::dict("name"_a="World", "number"_a=42);
        py::exec(R"(
            message = "Hello, {name}! The answer is {number}".format(**locals())
        )", py::globals(), locals);

        auto message = locals["message"].cast<std::string>();
        std::cout << message;
    }

Importing modules
=================

Python modules can be imported using `module_::import()`:

.. code-block:: cpp

    py::module_ sys = py::module_::import("sys");
    py::print(sys.attr("path"));

For convenience, the current working directory is included in ``sys.path`` when
embedding the interpreter. This makes it easy to import local Python files:

.. code-block:: python

    """calc.py located in the working directory"""

    def add(i, j):
        return i + j


.. code-block:: cpp

    py::module_ calc = py::module_::import("calc");
    py::object result = calc.attr("add")(1, 2);
    int n = result.cast<int>();
    assert(n == 3);

Modules can be reloaded using `module_::reload()` if the source is modified e.g.
by an external process. This can be useful in scenarios where the application
imports a user defined data processing script which needs to be updated after
changes by the user. Note that this function does not reload modules recursively.

.. _embedding_modules:

Adding embedded modules
=======================

Embedded binary modules can be added using the `PYBIND11_EMBEDDED_MODULE` macro.
Note that the definition must be placed at global scope. They can be imported
like any other module.

.. code-block:: cpp

    #include <pybind11/embed.h>
    namespace py = pybind11;

    PYBIND11_EMBEDDED_MODULE(fast_calc, m) {
        // `m` is a `py::module_` which is used to bind functions and classes
        m.def("add", [](int i, int j) {
            return i + j;
        });
    }

    int main() {
        py::scoped_interpreter guard{};

        auto fast_calc = py::module_::import("fast_calc");
        auto result = fast_calc.attr("add")(1, 2).cast<int>();
        assert(result == 3);
    }

Unlike extension modules where only a single binary module can be created, on
the embedded side an unlimited number of modules can be added using multiple
`PYBIND11_EMBEDDED_MODULE` definitions (as long as they have unique names).

These modules are added to Python's list of builtins, so they can also be
imported in pure Python files loaded by the interpreter. Everything interacts
naturally:

.. code-block:: python

    """py_module.py located in the working directory"""
    import cpp_module

    a = cpp_module.a
    b = a + 1


.. code-block:: cpp

    #include <pybind11/embed.h>
    namespace py = pybind11;

    PYBIND11_EMBEDDED_MODULE(cpp_module, m) {
        m.attr("a") = 1;
    }

    int main() {
        py::scoped_interpreter guard{};

        auto py_module = py::module_::import("py_module");

        auto locals = py::dict("fmt"_a="{} + {} = {}", **py_module.attr("__dict__"));
        assert(locals["a"].cast<int>() == 1);
        assert(locals["b"].cast<int>() == 2);

        py::exec(R"(
            c = a + b
            message = fmt.format(a, b, c)
        )", py::globals(), locals);

        assert(locals["c"].cast<int>() == 3);
        assert(locals["message"].cast<std::string>() == "1 + 2 = 3");
    }


Interpreter lifetime
====================

The Python interpreter shuts down when `scoped_interpreter` is destroyed. After
this, creating a new instance will restart the interpreter. Alternatively, the
`initialize_interpreter` / `finalize_interpreter` pair of functions can be used
to directly set the state at any time.

Modules created with pybind11 can be safely re-initialized after the interpreter
has been restarted. However, this may not apply to third-party extension modules.
The issue is that Python itself cannot completely unload extension modules and
there are several caveats with regard to interpreter restarting. In short, not
all memory may be freed, either due to Python reference cycles or user-created
global data. All the details can be found in the CPython documentation.

.. warning::

    Creating two concurrent `scoped_interpreter` guards is a fatal error. So is
    calling `initialize_interpreter` for a second time after the interpreter
    has already been initialized.

    Do not use the raw CPython API functions ``Py_Initialize`` and
    ``Py_Finalize`` as these do not properly handle the lifetime of
    pybind11's internal data.


Sub-interpreter support
=======================

Creating multiple copies of `scoped_interpreter` is not possible because it
represents the main Python interpreter. Sub-interpreters are something different
and they do permit the existence of multiple interpreters. This is an advanced
feature of the CPython API and should be handled with care. pybind11 does not
currently offer a C++ interface for sub-interpreters, so refer to the CPython
documentation for all the details regarding this feature.

We'll just mention a couple of caveats the sub-interpreters support in pybind11:

 1. Sub-interpreters will not receive independent copies of embedded modules.
    Instead, these are shared and modifications in one interpreter may be
    reflected in another.

 2. Managing multiple threads, multiple interpreters and the GIL can be
    challenging and there are several caveats here, even within the pure
    CPython API (please refer to the Python docs for details). As for
    pybind11, keep in mind that `gil_scoped_release` and `gil_scoped_acquire`
    do not take sub-interpreters into account.
.. _numpy:

NumPy
#####

Buffer protocol
===============

Python supports an extremely general and convenient approach for exchanging
data between plugin libraries. Types can expose a buffer view [#f2]_, which
provides fast direct access to the raw internal data representation. Suppose we
want to bind the following simplistic Matrix class:

.. code-block:: cpp

    class Matrix {
    public:
        Matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols) {
            m_data = new float[rows*cols];
        }
        float *data() { return m_data; }
        size_t rows() const { return m_rows; }
        size_t cols() const { return m_cols; }
    private:
        size_t m_rows, m_cols;
        float *m_data;
    };

The following binding code exposes the ``Matrix`` contents as a buffer object,
making it possible to cast Matrices into NumPy arrays. It is even possible to
completely avoid copy operations with Python expressions like
``np.array(matrix_instance, copy = False)``.

.. code-block:: cpp

    py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
       .def_buffer([](Matrix &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),                               /* Pointer to buffer */
                sizeof(float),                          /* Size of one scalar */
                py::format_descriptor<float>::format(), /* Python struct-style format descriptor */
                2,                                      /* Number of dimensions */
                { m.rows(), m.cols() },                 /* Buffer dimensions */
                { sizeof(float) * m.cols(),             /* Strides (in bytes) for each index */
                  sizeof(float) }
            );
        });

Supporting the buffer protocol in a new type involves specifying the special
``py::buffer_protocol()`` tag in the ``py::class_`` constructor and calling the
``def_buffer()`` method with a lambda function that creates a
``py::buffer_info`` description record on demand describing a given matrix
instance. The contents of ``py::buffer_info`` mirror the Python buffer protocol
specification.

.. code-block:: cpp

    struct buffer_info {
        void *ptr;
        py::ssize_t itemsize;
        std::string format;
        py::ssize_t ndim;
        std::vector<py::ssize_t> shape;
        std::vector<py::ssize_t> strides;
    };

To create a C++ function that can take a Python buffer object as an argument,
simply use the type ``py::buffer`` as one of its arguments. Buffers can exist
in a great variety of configurations, hence some safety checks are usually
necessary in the function body. Below, you can see a basic example on how to
define a custom constructor for the Eigen double precision matrix
(``Eigen::MatrixXd``) type, which supports initialization from compatible
buffer objects (e.g. a NumPy matrix).

.. code-block:: cpp

    /* Bind MatrixXd (or some other Eigen type) to Python */
    typedef Eigen::MatrixXd Matrix;

    typedef Matrix::Scalar Scalar;
    constexpr bool rowMajor = Matrix::Flags & Eigen::RowMajorBit;

    py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
        .def(py::init([](py::buffer b) {
            typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

            /* Request a buffer descriptor from Python */
            py::buffer_info info = b.request();

            /* Some sanity checks ... */
            if (info.format != py::format_descriptor<Scalar>::format())
                throw std::runtime_error("Incompatible format: expected a double array!");

            if (info.ndim != 2)
                throw std::runtime_error("Incompatible buffer dimension!");

            auto strides = Strides(
                info.strides[rowMajor ? 0 : 1] / (py::ssize_t)sizeof(Scalar),
                info.strides[rowMajor ? 1 : 0] / (py::ssize_t)sizeof(Scalar));

            auto map = Eigen::Map<Matrix, 0, Strides>(
                static_cast<Scalar *>(info.ptr), info.shape[0], info.shape[1], strides);

            return Matrix(map);
        }));

For reference, the ``def_buffer()`` call for this Eigen data type should look
as follows:

.. code-block:: cpp

    .def_buffer([](Matrix &m) -> py::buffer_info {
        return py::buffer_info(
            m.data(),                                /* Pointer to buffer */
            sizeof(Scalar),                          /* Size of one scalar */
            py::format_descriptor<Scalar>::format(), /* Python struct-style format descriptor */
            2,                                       /* Number of dimensions */
            { m.rows(), m.cols() },                  /* Buffer dimensions */
            { sizeof(Scalar) * (rowMajor ? m.cols() : 1),
              sizeof(Scalar) * (rowMajor ? 1 : m.rows()) }
                                                     /* Strides (in bytes) for each index */
        );
     })

For a much easier approach of binding Eigen types (although with some
limitations), refer to the section on :doc:`/advanced/cast/eigen`.

.. seealso::

    The file :file:`tests/test_buffers.cpp` contains a complete example
    that demonstrates using the buffer protocol with pybind11 in more detail.

.. [#f2] http://docs.python.org/3/c-api/buffer.html

Arrays
======

By exchanging ``py::buffer`` with ``py::array`` in the above snippet, we can
restrict the function so that it only accepts NumPy arrays (rather than any
type of Python object satisfying the buffer protocol).

In many situations, we want to define a function which only accepts a NumPy
array of a certain data type. This is possible via the ``py::array_t<T>``
template. For instance, the following function requires the argument to be a
NumPy array containing double precision values.

.. code-block:: cpp

    void f(py::array_t<double> array);

When it is invoked with a different type (e.g. an integer or a list of
integers), the binding code will attempt to cast the input into a NumPy array
of the requested type. This feature requires the :file:`pybind11/numpy.h`
header to be included. Note that :file:`pybind11/numpy.h` does not depend on
the NumPy headers, and thus can be used without declaring a build-time
dependency on NumPy; NumPy>=1.7.0 is a runtime dependency.

Data in NumPy arrays is not guaranteed to packed in a dense manner;
furthermore, entries can be separated by arbitrary column and row strides.
Sometimes, it can be useful to require a function to only accept dense arrays
using either the C (row-major) or Fortran (column-major) ordering. This can be
accomplished via a second template argument with values ``py::array::c_style``
or ``py::array::f_style``.

.. code-block:: cpp

    void f(py::array_t<double, py::array::c_style | py::array::forcecast> array);

The ``py::array::forcecast`` argument is the default value of the second
template parameter, and it ensures that non-conforming arguments are converted
into an array satisfying the specified requirements instead of trying the next
function overload.

Structured types
================

In order for ``py::array_t`` to work with structured (record) types, we first
need to register the memory layout of the type. This can be done via
``PYBIND11_NUMPY_DTYPE`` macro, called in the plugin definition code, which
expects the type followed by field names:

.. code-block:: cpp

    struct A {
        int x;
        double y;
    };

    struct B {
        int z;
        A a;
    };

    // ...
    PYBIND11_MODULE(test, m) {
        // ...

        PYBIND11_NUMPY_DTYPE(A, x, y);
        PYBIND11_NUMPY_DTYPE(B, z, a);
        /* now both A and B can be used as template arguments to py::array_t */
    }

The structure should consist of fundamental arithmetic types, ``std::complex``,
previously registered substructures, and arrays of any of the above. Both C++
arrays and ``std::array`` are supported. While there is a static assertion to
prevent many types of unsupported structures, it is still the user's
responsibility to use only "plain" structures that can be safely manipulated as
raw memory without violating invariants.

Vectorizing functions
=====================

Suppose we want to bind a function with the following signature to Python so
that it can process arbitrary NumPy array arguments (vectors, matrices, general
N-D arrays) in addition to its normal arguments:

.. code-block:: cpp

    double my_func(int x, float y, double z);

After including the ``pybind11/numpy.h`` header, this is extremely simple:

.. code-block:: cpp

    m.def("vectorized_func", py::vectorize(my_func));

Invoking the function like below causes 4 calls to be made to ``my_func`` with
each of the array elements. The significant advantage of this compared to
solutions like ``numpy.vectorize()`` is that the loop over the elements runs
entirely on the C++ side and can be crunched down into a tight, optimized loop
by the compiler. The result is returned as a NumPy array of type
``numpy.dtype.float64``.

.. code-block:: pycon

    >>> x = np.array([[1, 3],[5, 7]])
    >>> y = np.array([[2, 4],[6, 8]])
    >>> z = 3
    >>> result = vectorized_func(x, y, z)

The scalar argument ``z`` is transparently replicated 4 times.  The input
arrays ``x`` and ``y`` are automatically converted into the right types (they
are of type  ``numpy.dtype.int64`` but need to be ``numpy.dtype.int32`` and
``numpy.dtype.float32``, respectively).

.. note::

    Only arithmetic, complex, and POD types passed by value or by ``const &``
    reference are vectorized; all other arguments are passed through as-is.
    Functions taking rvalue reference arguments cannot be vectorized.

In cases where the computation is too complicated to be reduced to
``vectorize``, it will be necessary to create and access the buffer contents
manually. The following snippet contains a complete example that shows how this
works (the code is somewhat contrived, since it could have been done more
simply using ``vectorize``).

.. code-block:: cpp

    #include <pybind11/pybind11.h>
    #include <pybind11/numpy.h>

    namespace py = pybind11;

    py::array_t<double> add_arrays(py::array_t<double> input1, py::array_t<double> input2) {
        py::buffer_info buf1 = input1.request(), buf2 = input2.request();

        if (buf1.ndim != 1 || buf2.ndim != 1)
            throw std::runtime_error("Number of dimensions must be one");

        if (buf1.size != buf2.size)
            throw std::runtime_error("Input shapes must match");

        /* No pointer is passed, so NumPy will allocate the buffer */
        auto result = py::array_t<double>(buf1.size);

        py::buffer_info buf3 = result.request();

        double *ptr1 = static_cast<double *>(buf1.ptr);
        double *ptr2 = static_cast<double *>(buf2.ptr);
        double *ptr3 = static_cast<double *>(buf3.ptr);

        for (size_t idx = 0; idx < buf1.shape[0]; idx++)
            ptr3[idx] = ptr1[idx] + ptr2[idx];

        return result;
    }

    PYBIND11_MODULE(test, m) {
        m.def("add_arrays", &add_arrays, "Add two NumPy arrays");
    }

.. seealso::

    The file :file:`tests/test_numpy_vectorize.cpp` contains a complete
    example that demonstrates using :func:`vectorize` in more detail.

Direct access
=============

For performance reasons, particularly when dealing with very large arrays, it
is often desirable to directly access array elements without internal checking
of dimensions and bounds on every access when indices are known to be already
valid.  To avoid such checks, the ``array`` class and ``array_t<T>`` template
class offer an unchecked proxy object that can be used for this unchecked
access through the ``unchecked<N>`` and ``mutable_unchecked<N>`` methods,
where ``N`` gives the required dimensionality of the array:

.. code-block:: cpp

    m.def("sum_3d", [](py::array_t<double> x) {
        auto r = x.unchecked<3>(); // x must have ndim = 3; can be non-writeable
        double sum = 0;
        for (py::ssize_t i = 0; i < r.shape(0); i++)
            for (py::ssize_t j = 0; j < r.shape(1); j++)
                for (py::ssize_t k = 0; k < r.shape(2); k++)
                    sum += r(i, j, k);
        return sum;
    });
    m.def("increment_3d", [](py::array_t<double> x) {
        auto r = x.mutable_unchecked<3>(); // Will throw if ndim != 3 or flags.writeable is false
        for (py::ssize_t i = 0; i < r.shape(0); i++)
            for (py::ssize_t j = 0; j < r.shape(1); j++)
                for (py::ssize_t k = 0; k < r.shape(2); k++)
                    r(i, j, k) += 1.0;
    }, py::arg().noconvert());

To obtain the proxy from an ``array`` object, you must specify both the data
type and number of dimensions as template arguments, such as ``auto r =
myarray.mutable_unchecked<float, 2>()``.

If the number of dimensions is not known at compile time, you can omit the
dimensions template parameter (i.e. calling ``arr_t.unchecked()`` or
``arr.unchecked<T>()``.  This will give you a proxy object that works in the
same way, but results in less optimizable code and thus a small efficiency
loss in tight loops.

Note that the returned proxy object directly references the array's data, and
only reads its shape, strides, and writeable flag when constructed.  You must
take care to ensure that the referenced array is not destroyed or reshaped for
the duration of the returned object, typically by limiting the scope of the
returned instance.

The returned proxy object supports some of the same methods as ``py::array`` so
that it can be used as a drop-in replacement for some existing, index-checked
uses of ``py::array``:

- ``r.ndim()`` returns the number of dimensions

- ``r.data(1, 2, ...)`` and ``r.mutable_data(1, 2, ...)``` returns a pointer to
  the ``const T`` or ``T`` data, respectively, at the given indices.  The
  latter is only available to proxies obtained via ``a.mutable_unchecked()``.

- ``itemsize()`` returns the size of an item in bytes, i.e. ``sizeof(T)``.

- ``ndim()`` returns the number of dimensions.

- ``shape(n)`` returns the size of dimension ``n``

- ``size()`` returns the total number of elements (i.e. the product of the shapes).

- ``nbytes()`` returns the number of bytes used by the referenced elements
  (i.e. ``itemsize()`` times ``size()``).

.. seealso::

    The file :file:`tests/test_numpy_array.cpp` contains additional examples
    demonstrating the use of this feature.

Ellipsis
========

Python 3 provides a convenient ``...`` ellipsis notation that is often used to
slice multidimensional arrays. For instance, the following snippet extracts the
middle dimensions of a tensor with the first and last index set to zero.
In Python 2, the syntactic sugar ``...`` is not available, but the singleton
``Ellipsis`` (of type ``ellipsis``) can still be used directly.

.. code-block:: python

   a = # a NumPy array
   b = a[0, ..., 0]

The function ``py::ellipsis()`` function can be used to perform the same
operation on the C++ side:

.. code-block:: cpp

   py::array a = /* A NumPy array */;
   py::array b = a[py::make_tuple(0, py::ellipsis(), 0)];

.. versionchanged:: 2.6
   ``py::ellipsis()`` is now also avaliable in Python 2.

Memory view
===========

For a case when we simply want to provide a direct accessor to C/C++ buffer
without a concrete class object, we can return a ``memoryview`` object. Suppose
we wish to expose a ``memoryview`` for 2x4 uint8_t array, we can do the
following:

.. code-block:: cpp

    const uint8_t buffer[] = {
        0, 1, 2, 3,
        4, 5, 6, 7
    };
    m.def("get_memoryview2d", []() {
        return py::memoryview::from_buffer(
            buffer,                                    // buffer pointer
            { 2, 4 },                                  // shape (rows, cols)
            { sizeof(uint8_t) * 4, sizeof(uint8_t) }   // strides in bytes
        );
    })

This approach is meant for providing a ``memoryview`` for a C/C++ buffer not
managed by Python. The user is responsible for managing the lifetime of the
buffer. Using a ``memoryview`` created in this way after deleting the buffer in
C++ side results in undefined behavior.

We can also use ``memoryview::from_memory`` for a simple 1D contiguous buffer:

.. code-block:: cpp

    m.def("get_memoryview1d", []() {
        return py::memoryview::from_memory(
            buffer,               // buffer pointer
            sizeof(uint8_t) * 8   // buffer size
        );
    })

.. note::

    ``memoryview::from_memory`` is not available in Python 2.

.. versionchanged:: 2.6
    ``memoryview::from_memory`` added.
Utilities
#########

Using Python's print function in C++
====================================

The usual way to write output in C++ is using ``std::cout`` while in Python one
would use ``print``. Since these methods use different buffers, mixing them can
lead to output order issues. To resolve this, pybind11 modules can use the
:func:`py::print` function which writes to Python's ``sys.stdout`` for consistency.

Python's ``print`` function is replicated in the C++ API including optional
keyword arguments ``sep``, ``end``, ``file``, ``flush``. Everything works as
expected in Python:

.. code-block:: cpp

    py::print(1, 2.0, "three"); // 1 2.0 three
    py::print(1, 2.0, "three", "sep"_a="-"); // 1-2.0-three

    auto args = py::make_tuple("unpacked", true);
    py::print("->", *args, "end"_a="<-"); // -> unpacked True <-

.. _ostream_redirect:

Capturing standard output from ostream
======================================

Often, a library will use the streams ``std::cout`` and ``std::cerr`` to print,
but this does not play well with Python's standard ``sys.stdout`` and ``sys.stderr``
redirection. Replacing a library's printing with `py::print <print>` may not
be feasible. This can be fixed using a guard around the library function that
redirects output to the corresponding Python streams:

.. code-block:: cpp

    #include <pybind11/iostream.h>

    ...

    // Add a scoped redirect for your noisy code
    m.def("noisy_func", []() {
        py::scoped_ostream_redirect stream(
            std::cout,                               // std::ostream&
            py::module_::import("sys").attr("stdout") // Python output
        );
        call_noisy_func();
    });

This method respects flushes on the output streams and will flush if needed
when the scoped guard is destroyed. This allows the output to be redirected in
real time, such as to a Jupyter notebook. The two arguments, the C++ stream and
the Python output, are optional, and default to standard output if not given. An
extra type, `py::scoped_estream_redirect <scoped_estream_redirect>`, is identical
except for defaulting to ``std::cerr`` and ``sys.stderr``; this can be useful with
`py::call_guard`, which allows multiple items, but uses the default constructor:

.. code-block:: py

    // Alternative: Call single function using call guard
    m.def("noisy_func", &call_noisy_function,
          py::call_guard<py::scoped_ostream_redirect,
                         py::scoped_estream_redirect>());

The redirection can also be done in Python with the addition of a context
manager, using the `py::add_ostream_redirect() <add_ostream_redirect>` function:

.. code-block:: cpp

    py::add_ostream_redirect(m, "ostream_redirect");

The name in Python defaults to ``ostream_redirect`` if no name is passed.  This
creates the following context manager in Python:

.. code-block:: python

    with ostream_redirect(stdout=True, stderr=True):
        noisy_function()

It defaults to redirecting both streams, though you can use the keyword
arguments to disable one of the streams if needed.

.. note::

    The above methods will not redirect C-level output to file descriptors, such
    as ``fprintf``. For those cases, you'll need to redirect the file
    descriptors either directly in C or with Python's ``os.dup2`` function
    in an operating-system dependent way.

.. _eval:

Evaluating Python expressions from strings and files
====================================================

pybind11 provides the `eval`, `exec` and `eval_file` functions to evaluate
Python expressions and statements. The following example illustrates how they
can be used.

.. code-block:: cpp

    // At beginning of file
    #include <pybind11/eval.h>

    ...

    // Evaluate in scope of main module
    py::object scope = py::module_::import("__main__").attr("__dict__");

    // Evaluate an isolated expression
    int result = py::eval("my_variable + 10", scope).cast<int>();

    // Evaluate a sequence of statements
    py::exec(
        "print('Hello')\n"
        "print('world!');",
        scope);

    // Evaluate the statements in an separate Python file on disk
    py::eval_file("script.py", scope);

C++11 raw string literals are also supported and quite handy for this purpose.
The only requirement is that the first statement must be on a new line following
the raw string delimiter ``R"(``, ensuring all lines have common leading indent:

.. code-block:: cpp

    py::exec(R"(
        x = get_answer()
        if x == 42:
            print('Hello World!')
        else:
            print('Bye!')
        )", scope
    );

.. note::

    `eval` and `eval_file` accept a template parameter that describes how the
    string/file should be interpreted. Possible choices include ``eval_expr``
    (isolated expression), ``eval_single_statement`` (a single statement, return
    value is always ``none``), and ``eval_statements`` (sequence of statements,
    return value is always ``none``). `eval` defaults to  ``eval_expr``,
    `eval_file` defaults to ``eval_statements`` and `exec` is just a shortcut
    for ``eval<eval_statements>``.
Python types
############

.. _wrappers:

Available wrappers
==================

All major Python types are available as thin C++ wrapper classes. These
can also be used as function parameters -- see :ref:`python_objects_as_args`.

Available types include :class:`handle`, :class:`object`, :class:`bool_`,
:class:`int_`, :class:`float_`, :class:`str`, :class:`bytes`, :class:`tuple`,
:class:`list`, :class:`dict`, :class:`slice`, :class:`none`, :class:`capsule`,
:class:`iterable`, :class:`iterator`, :class:`function`, :class:`buffer`,
:class:`array`, and :class:`array_t`.

.. warning::

    Be sure to review the :ref:`pytypes_gotchas` before using this heavily in
    your C++ API.

.. _casting_back_and_forth:

Casting back and forth
======================

In this kind of mixed code, it is often necessary to convert arbitrary C++
types to Python, which can be done using :func:`py::cast`:

.. code-block:: cpp

    MyClass *cls = ..;
    py::object obj = py::cast(cls);

The reverse direction uses the following syntax:

.. code-block:: cpp

    py::object obj = ...;
    MyClass *cls = obj.cast<MyClass *>();

When conversion fails, both directions throw the exception :class:`cast_error`.

.. _python_libs:

Accessing Python libraries from C++
===================================

It is also possible to import objects defined in the Python standard
library or available in the current Python environment (``sys.path``) and work
with these in C++.

This example obtains a reference to the Python ``Decimal`` class.

.. code-block:: cpp

    // Equivalent to "from decimal import Decimal"
    py::object Decimal = py::module_::import("decimal").attr("Decimal");

.. code-block:: cpp

    // Try to import scipy
    py::object scipy = py::module_::import("scipy");
    return scipy.attr("__version__");


.. _calling_python_functions:

Calling Python functions
========================

It is also possible to call Python classes, functions and methods
via ``operator()``.

.. code-block:: cpp

    // Construct a Python object of class Decimal
    py::object pi = Decimal("3.14159");

.. code-block:: cpp

    // Use Python to make our directories
    py::object os = py::module_::import("os");
    py::object makedirs = os.attr("makedirs");
    makedirs("/tmp/path/to/somewhere");

One can convert the result obtained from Python to a pure C++ version
if a ``py::class_`` or type conversion is defined.

.. code-block:: cpp

    py::function f = <...>;
    py::object result_py = f(1234, "hello", some_instance);
    MyClass &result = result_py.cast<MyClass>();

.. _calling_python_methods:

Calling Python methods
========================

To call an object's method, one can again use ``.attr`` to obtain access to the
Python method.

.. code-block:: cpp

    // Calculate e^π in decimal
    py::object exp_pi = pi.attr("exp")();
    py::print(py::str(exp_pi));

In the example above ``pi.attr("exp")`` is a *bound method*: it will always call
the method for that same instance of the class. Alternately one can create an
*unbound method* via the Python class (instead of instance) and pass the ``self``
object explicitly, followed by other arguments.

.. code-block:: cpp

    py::object decimal_exp = Decimal.attr("exp");

    // Compute the e^n for n=0..4
    for (int n = 0; n < 5; n++) {
        py::print(decimal_exp(Decimal(n));
    }

Keyword arguments
=================

Keyword arguments are also supported. In Python, there is the usual call syntax:

.. code-block:: python

    def f(number, say, to):
        ...  # function code

    f(1234, say="hello", to=some_instance)  # keyword call in Python

In C++, the same call can be made using:

.. code-block:: cpp

    using namespace pybind11::literals; // to bring in the `_a` literal
    f(1234, "say"_a="hello", "to"_a=some_instance); // keyword call in C++

Unpacking arguments
===================

Unpacking of ``*args`` and ``**kwargs`` is also possible and can be mixed with
other arguments:

.. code-block:: cpp

    // * unpacking
    py::tuple args = py::make_tuple(1234, "hello", some_instance);
    f(*args);

    // ** unpacking
    py::dict kwargs = py::dict("number"_a=1234, "say"_a="hello", "to"_a=some_instance);
    f(**kwargs);

    // mixed keywords, * and ** unpacking
    py::tuple args = py::make_tuple(1234);
    py::dict kwargs = py::dict("to"_a=some_instance);
    f(*args, "say"_a="hello", **kwargs);

Generalized unpacking according to PEP448_ is also supported:

.. code-block:: cpp

    py::dict kwargs1 = py::dict("number"_a=1234);
    py::dict kwargs2 = py::dict("to"_a=some_instance);
    f(**kwargs1, "say"_a="hello", **kwargs2);

.. seealso::

    The file :file:`tests/test_pytypes.cpp` contains a complete
    example that demonstrates passing native Python types in more detail. The
    file :file:`tests/test_callbacks.cpp` presents a few examples of calling
    Python functions from C++, including keywords arguments and unpacking.

.. _PEP448: https://www.python.org/dev/peps/pep-0448/

.. _implicit_casting:

Implicit casting
================

When using the C++ interface for Python types, or calling Python functions,
objects of type :class:`object` are returned. It is possible to invoke implicit
conversions to subclasses like :class:`dict`. The same holds for the proxy objects
returned by ``operator[]`` or ``obj.attr()``.
Casting to subtypes improves code readability and allows values to be passed to
C++ functions that require a specific subtype rather than a generic :class:`object`.

.. code-block:: cpp

    #include <pybind11/numpy.h>
    using namespace pybind11::literals;

    py::module_ os = py::module_::import("os");
    py::module_ path = py::module_::import("os.path");  // like 'import os.path as path'
    py::module_ np = py::module_::import("numpy");  // like 'import numpy as np'

    py::str curdir_abs = path.attr("abspath")(path.attr("curdir"));
    py::print(py::str("Current directory: ") + curdir_abs);
    py::dict environ = os.attr("environ");
    py::print(environ["HOME"]);
    py::array_t<float> arr = np.attr("ones")(3, "dtype"_a="float32");
    py::print(py::repr(arr + py::int_(1)));

These implicit conversions are available for subclasses of :class:`object`; there
is no need to call ``obj.cast()`` explicitly as for custom classes, see
:ref:`casting_back_and_forth`.

.. note::
    If a trivial conversion via move constructor is not possible, both implicit and
    explicit casting (calling ``obj.cast()``) will attempt a "rich" conversion.
    For instance, ``py::list env = os.attr("environ");`` will succeed and is
    equivalent to the Python code ``env = list(os.environ)`` that produces a
    list of the dict keys.

..  TODO: Adapt text once PR #2349 has landed

Handling exceptions
===================

Python exceptions from wrapper classes will be thrown as a ``py::error_already_set``.
See :ref:`Handling exceptions from Python in C++
<handling_python_exceptions_cpp>` for more information on handling exceptions
raised when calling C++ wrapper classes.

.. _pytypes_gotchas:

Gotchas
=======

Default-Constructed Wrappers
----------------------------

When a wrapper type is default-constructed, it is **not** a valid Python object (i.e. it is not ``py::none()``). It is simply the same as
``PyObject*`` null pointer. To check for this, use
``static_cast<bool>(my_wrapper)``.

Assigning py::none() to wrappers
--------------------------------

You may be tempted to use types like ``py::str`` and ``py::dict`` in C++
signatures (either pure C++, or in bound signatures), and assign them default
values of ``py::none()``. However, in a best case scenario, it will fail fast
because ``None`` is not convertible to that type (e.g. ``py::dict``), or in a
worse case scenario, it will silently work but corrupt the types you want to
work with (e.g. ``py::str(py::none())`` will yield ``"None"`` in Python).
Python C++ interface
####################

pybind11 exposes Python types and functions using thin C++ wrappers, which
makes it possible to conveniently call Python code from C++ without resorting
to Python's C API.

.. toctree::
   :maxdepth: 2

   object
   numpy
   utilities
STL containers
##############

Automatic conversion
====================

When including the additional header file :file:`pybind11/stl.h`, conversions
between ``std::vector<>``/``std::deque<>``/``std::list<>``/``std::array<>``/``std::valarray<>``,
``std::set<>``/``std::unordered_set<>``, and
``std::map<>``/``std::unordered_map<>`` and the Python ``list``, ``set`` and
``dict`` data structures are automatically enabled. The types ``std::pair<>``
and ``std::tuple<>`` are already supported out of the box with just the core
:file:`pybind11/pybind11.h` header.

The major downside of these implicit conversions is that containers must be
converted (i.e. copied) on every Python->C++ and C++->Python transition, which
can have implications on the program semantics and performance. Please read the
next sections for more details and alternative approaches that avoid this.

.. note::

    Arbitrary nesting of any of these types is possible.

.. seealso::

    The file :file:`tests/test_stl.cpp` contains a complete
    example that demonstrates how to pass STL data types in more detail.

.. _cpp17_container_casters:

C++17 library containers
========================

The :file:`pybind11/stl.h` header also includes support for ``std::optional<>``
and ``std::variant<>``. These require a C++17 compiler and standard library.
In C++14 mode, ``std::experimental::optional<>`` is supported if available.

Various versions of these containers also exist for C++11 (e.g. in Boost).
pybind11 provides an easy way to specialize the ``type_caster`` for such
types:

.. code-block:: cpp

    // `boost::optional` as an example -- can be any `std::optional`-like container
    namespace pybind11 { namespace detail {
        template <typename T>
        struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>> {};
    }}

The above should be placed in a header file and included in all translation units
where automatic conversion is needed. Similarly, a specialization can be provided
for custom variant types:

.. code-block:: cpp

    // `boost::variant` as an example -- can be any `std::variant`-like container
    namespace pybind11 { namespace detail {
        template <typename... Ts>
        struct type_caster<boost::variant<Ts...>> : variant_caster<boost::variant<Ts...>> {};

        // Specifies the function used to visit the variant -- `apply_visitor` instead of `visit`
        template <>
        struct visit_helper<boost::variant> {
            template <typename... Args>
            static auto call(Args &&...args) -> decltype(boost::apply_visitor(args...)) {
                return boost::apply_visitor(args...);
            }
        };
    }} // namespace pybind11::detail

The ``visit_helper`` specialization is not required if your ``name::variant`` provides
a ``name::visit()`` function. For any other function name, the specialization must be
included to tell pybind11 how to visit the variant.

.. warning::

    When converting a ``variant`` type, pybind11 follows the same rules as when
    determining which function overload to call (:ref:`overload_resolution`), and
    so the same caveats hold. In particular, the order in which the ``variant``'s
    alternatives are listed is important, since pybind11 will try conversions in
    this order. This means that, for example, when converting ``variant<int, bool>``,
    the ``bool`` variant will never be selected, as any Python ``bool`` is already
    an ``int`` and is convertible to a C++ ``int``. Changing the order of alternatives
    (and using ``variant<bool, int>``, in this example) provides a solution.

.. note::

    pybind11 only supports the modern implementation of ``boost::variant``
    which makes use of variadic templates. This requires Boost 1.56 or newer.
    Additionally, on Windows, MSVC 2017 is required because ``boost::variant``
    falls back to the old non-variadic implementation on MSVC 2015.

.. _opaque:

Making opaque types
===================

pybind11 heavily relies on a template matching mechanism to convert parameters
and return values that are constructed from STL data types such as vectors,
linked lists, hash tables, etc. This even works in a recursive manner, for
instance to deal with lists of hash maps of pairs of elementary and custom
types, etc.

However, a fundamental limitation of this approach is that internal conversions
between Python and C++ types involve a copy operation that prevents
pass-by-reference semantics. What does this mean?

Suppose we bind the following function

.. code-block:: cpp

    void append_1(std::vector<int> &v) {
       v.push_back(1);
    }

and call it from Python, the following happens:

.. code-block:: pycon

   >>> v = [5, 6]
   >>> append_1(v)
   >>> print(v)
   [5, 6]

As you can see, when passing STL data structures by reference, modifications
are not propagated back the Python side. A similar situation arises when
exposing STL data structures using the ``def_readwrite`` or ``def_readonly``
functions:

.. code-block:: cpp

    /* ... definition ... */

    class MyClass {
        std::vector<int> contents;
    };

    /* ... binding code ... */

    py::class_<MyClass>(m, "MyClass")
        .def(py::init<>())
        .def_readwrite("contents", &MyClass::contents);

In this case, properties can be read and written in their entirety. However, an
``append`` operation involving such a list type has no effect:

.. code-block:: pycon

   >>> m = MyClass()
   >>> m.contents = [5, 6]
   >>> print(m.contents)
   [5, 6]
   >>> m.contents.append(7)
   >>> print(m.contents)
   [5, 6]

Finally, the involved copy operations can be costly when dealing with very
large lists. To deal with all of the above situations, pybind11 provides a
macro named ``PYBIND11_MAKE_OPAQUE(T)`` that disables the template-based
conversion machinery of types, thus rendering them *opaque*. The contents of
opaque objects are never inspected or extracted, hence they *can* be passed by
reference. For instance, to turn ``std::vector<int>`` into an opaque type, add
the declaration

.. code-block:: cpp

    PYBIND11_MAKE_OPAQUE(std::vector<int>);

before any binding code (e.g. invocations to ``class_::def()``, etc.). This
macro must be specified at the top level (and outside of any namespaces), since
it adds a template instantiation of ``type_caster``. If your binding code consists of
multiple compilation units, it must be present in every file (typically via a
common header) preceding any usage of ``std::vector<int>``. Opaque types must
also have a corresponding ``class_`` declaration to associate them with a name
in Python, and to define a set of available operations, e.g.:

.. code-block:: cpp

    py::class_<std::vector<int>>(m, "IntVector")
        .def(py::init<>())
        .def("clear", &std::vector<int>::clear)
        .def("pop_back", &std::vector<int>::pop_back)
        .def("__len__", [](const std::vector<int> &v) { return v.size(); })
        .def("__iter__", [](std::vector<int> &v) {
           return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */
        // ....

.. seealso::

    The file :file:`tests/test_opaque_types.cpp` contains a complete
    example that demonstrates how to create and expose opaque types using
    pybind11 in more detail.

.. _stl_bind:

Binding STL containers
======================

The ability to expose STL containers as native Python objects is a fairly
common request, hence pybind11 also provides an optional header file named
:file:`pybind11/stl_bind.h` that does exactly this. The mapped containers try
to match the behavior of their native Python counterparts as much as possible.

The following example showcases usage of :file:`pybind11/stl_bind.h`:

.. code-block:: cpp

    // Don't forget this
    #include <pybind11/stl_bind.h>

    PYBIND11_MAKE_OPAQUE(std::vector<int>);
    PYBIND11_MAKE_OPAQUE(std::map<std::string, double>);

    // ...

    // later in binding code:
    py::bind_vector<std::vector<int>>(m, "VectorInt");
    py::bind_map<std::map<std::string, double>>(m, "MapStringDouble");

When binding STL containers pybind11 considers the types of the container's
elements to decide whether the container should be confined to the local module
(via the :ref:`module_local` feature).  If the container element types are
anything other than already-bound custom types bound without
``py::module_local()`` the container binding will have ``py::module_local()``
applied.  This includes converting types such as numeric types, strings, Eigen
types; and types that have not yet been bound at the time of the stl container
binding.  This module-local binding is designed to avoid potential conflicts
between module bindings (for example, from two separate modules each attempting
to bind ``std::vector<int>`` as a python type).

It is possible to override this behavior to force a definition to be either
module-local or global.  To do so, you can pass the attributes
``py::module_local()`` (to make the binding module-local) or
``py::module_local(false)`` (to make the binding global) into the
``py::bind_vector`` or ``py::bind_map`` arguments:

.. code-block:: cpp

    py::bind_vector<std::vector<int>>(m, "VectorInt", py::module_local(false));

Note, however, that such a global binding would make it impossible to load this
module at the same time as any other pybind module that also attempts to bind
the same container type (``std::vector<int>`` in the above example).

See :ref:`module_local` for more details on module-local bindings.

.. seealso::

    The file :file:`tests/test_stl_binders.cpp` shows how to use the
    convenience STL container wrappers.
Eigen
#####

`Eigen <http://eigen.tuxfamily.org>`_ is C++ header-based library for dense and
sparse linear algebra. Due to its popularity and widespread adoption, pybind11
provides transparent conversion and limited mapping support between Eigen and
Scientific Python linear algebra data types.

To enable the built-in Eigen support you must include the optional header file
:file:`pybind11/eigen.h`.

Pass-by-value
=============

When binding a function with ordinary Eigen dense object arguments (for
example, ``Eigen::MatrixXd``), pybind11 will accept any input value that is
already (or convertible to) a ``numpy.ndarray`` with dimensions compatible with
the Eigen type, copy its values into a temporary Eigen variable of the
appropriate type, then call the function with this temporary variable.

Sparse matrices are similarly copied to or from
``scipy.sparse.csr_matrix``/``scipy.sparse.csc_matrix`` objects.

Pass-by-reference
=================

One major limitation of the above is that every data conversion implicitly
involves a copy, which can be both expensive (for large matrices) and disallows
binding functions that change their (Matrix) arguments.  Pybind11 allows you to
work around this by using Eigen's ``Eigen::Ref<MatrixType>`` class much as you
would when writing a function taking a generic type in Eigen itself (subject to
some limitations discussed below).

When calling a bound function accepting a ``Eigen::Ref<const MatrixType>``
type, pybind11 will attempt to avoid copying by using an ``Eigen::Map`` object
that maps into the source ``numpy.ndarray`` data: this requires both that the
data types are the same (e.g. ``dtype='float64'`` and ``MatrixType::Scalar`` is
``double``); and that the storage is layout compatible.  The latter limitation
is discussed in detail in the section below, and requires careful
consideration: by default, numpy matrices and Eigen matrices are *not* storage
compatible.

If the numpy matrix cannot be used as is (either because its types differ, e.g.
passing an array of integers to an Eigen parameter requiring doubles, or
because the storage is incompatible), pybind11 makes a temporary copy and
passes the copy instead.

When a bound function parameter is instead ``Eigen::Ref<MatrixType>`` (note the
lack of ``const``), pybind11 will only allow the function to be called if it
can be mapped *and* if the numpy array is writeable (that is
``a.flags.writeable`` is true).  Any access (including modification) made to
the passed variable will be transparently carried out directly on the
``numpy.ndarray``.

This means you can can write code such as the following and have it work as
expected:

.. code-block:: cpp

    void scale_by_2(Eigen::Ref<Eigen::VectorXd> v) {
        v *= 2;
    }

Note, however, that you will likely run into limitations due to numpy and
Eigen's difference default storage order for data; see the below section on
:ref:`storage_orders` for details on how to bind code that won't run into such
limitations.

.. note::

    Passing by reference is not supported for sparse types.

Returning values to Python
==========================

When returning an ordinary dense Eigen matrix type to numpy (e.g.
``Eigen::MatrixXd`` or ``Eigen::RowVectorXf``) pybind11 keeps the matrix and
returns a numpy array that directly references the Eigen matrix: no copy of the
data is performed.  The numpy array will have ``array.flags.owndata`` set to
``False`` to indicate that it does not own the data, and the lifetime of the
stored Eigen matrix will be tied to the returned ``array``.

If you bind a function with a non-reference, ``const`` return type (e.g.
``const Eigen::MatrixXd``), the same thing happens except that pybind11 also
sets the numpy array's ``writeable`` flag to false.

If you return an lvalue reference or pointer, the usual pybind11 rules apply,
as dictated by the binding function's return value policy (see the
documentation on :ref:`return_value_policies` for full details).  That means,
without an explicit return value policy, lvalue references will be copied and
pointers will be managed by pybind11.  In order to avoid copying, you should
explicitly specify an appropriate return value policy, as in the following
example:

.. code-block:: cpp

    class MyClass {
        Eigen::MatrixXd big_mat = Eigen::MatrixXd::Zero(10000, 10000);
    public:
        Eigen::MatrixXd &getMatrix() { return big_mat; }
        const Eigen::MatrixXd &viewMatrix() { return big_mat; }
    };

    // Later, in binding code:
    py::class_<MyClass>(m, "MyClass")
        .def(py::init<>())
        .def("copy_matrix", &MyClass::getMatrix) // Makes a copy!
        .def("get_matrix", &MyClass::getMatrix, py::return_value_policy::reference_internal)
        .def("view_matrix", &MyClass::viewMatrix, py::return_value_policy::reference_internal)
        ;

.. code-block:: python

    a = MyClass()
    m = a.get_matrix()   # flags.writeable = True,  flags.owndata = False
    v = a.view_matrix()  # flags.writeable = False, flags.owndata = False
    c = a.copy_matrix()  # flags.writeable = True,  flags.owndata = True
    # m[5,6] and v[5,6] refer to the same element, c[5,6] does not.

Note in this example that ``py::return_value_policy::reference_internal`` is
used to tie the life of the MyClass object to the life of the returned arrays.

You may also return an ``Eigen::Ref``, ``Eigen::Map`` or other map-like Eigen
object (for example, the return value of ``matrix.block()`` and related
methods) that map into a dense Eigen type.  When doing so, the default
behaviour of pybind11 is to simply reference the returned data: you must take
care to ensure that this data remains valid!  You may ask pybind11 to
explicitly *copy* such a return value by using the
``py::return_value_policy::copy`` policy when binding the function.  You may
also use ``py::return_value_policy::reference_internal`` or a
``py::keep_alive`` to ensure the data stays valid as long as the returned numpy
array does.

When returning such a reference of map, pybind11 additionally respects the
readonly-status of the returned value, marking the numpy array as non-writeable
if the reference or map was itself read-only.

.. note::

    Sparse types are always copied when returned.

.. _storage_orders:

Storage orders
==============

Passing arguments via ``Eigen::Ref`` has some limitations that you must be
aware of in order to effectively pass matrices by reference.  First and
foremost is that the default ``Eigen::Ref<MatrixType>`` class requires
contiguous storage along columns (for column-major types, the default in Eigen)
or rows if ``MatrixType`` is specifically an ``Eigen::RowMajor`` storage type.
The former, Eigen's default, is incompatible with ``numpy``'s default row-major
storage, and so you will not be able to pass numpy arrays to Eigen by reference
without making one of two changes.

(Note that this does not apply to vectors (or column or row matrices): for such
types the "row-major" and "column-major" distinction is meaningless).

The first approach is to change the use of ``Eigen::Ref<MatrixType>`` to the
more general ``Eigen::Ref<MatrixType, 0, Eigen::Stride<Eigen::Dynamic,
Eigen::Dynamic>>`` (or similar type with a fully dynamic stride type in the
third template argument).  Since this is a rather cumbersome type, pybind11
provides a ``py::EigenDRef<MatrixType>`` type alias for your convenience (along
with EigenDMap for the equivalent Map, and EigenDStride for just the stride
type).

This type allows Eigen to map into any arbitrary storage order.  This is not
the default in Eigen for performance reasons: contiguous storage allows
vectorization that cannot be done when storage is not known to be contiguous at
compile time.  The default ``Eigen::Ref`` stride type allows non-contiguous
storage along the outer dimension (that is, the rows of a column-major matrix
or columns of a row-major matrix), but not along the inner dimension.

This type, however, has the added benefit of also being able to map numpy array
slices.  For example, the following (contrived) example uses Eigen with a numpy
slice to multiply by 2 all coefficients that are both on even rows (0, 2, 4,
...) and in columns 2, 5, or 8:

.. code-block:: cpp

    m.def("scale", [](py::EigenDRef<Eigen::MatrixXd> m, double c) { m *= c; });

.. code-block:: python

    # a = np.array(...)
    scale_by_2(myarray[0::2, 2:9:3])

The second approach to avoid copying is more intrusive: rearranging the
underlying data types to not run into the non-contiguous storage problem in the
first place.  In particular, that means using matrices with ``Eigen::RowMajor``
storage, where appropriate, such as:

.. code-block:: cpp

    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    // Use RowMatrixXd instead of MatrixXd

Now bound functions accepting ``Eigen::Ref<RowMatrixXd>`` arguments will be
callable with numpy's (default) arrays without involving a copying.

You can, alternatively, change the storage order that numpy arrays use by
adding the ``order='F'`` option when creating an array:

.. code-block:: python

    myarray = np.array(source, order='F')

Such an object will be passable to a bound function accepting an
``Eigen::Ref<MatrixXd>`` (or similar column-major Eigen type).

One major caveat with this approach, however, is that it is not entirely as
easy as simply flipping all Eigen or numpy usage from one to the other: some
operations may alter the storage order of a numpy array.  For example, ``a2 =
array.transpose()`` results in ``a2`` being a view of ``array`` that references
the same data, but in the opposite storage order!

While this approach allows fully optimized vectorized calculations in Eigen, it
cannot be used with array slices, unlike the first approach.

When *returning* a matrix to Python (either a regular matrix, a reference via
``Eigen::Ref<>``, or a map/block into a matrix), no special storage
consideration is required: the created numpy array will have the required
stride that allows numpy to properly interpret the array, whatever its storage
order.

Failing rather than copying
===========================

The default behaviour when binding ``Eigen::Ref<const MatrixType>`` Eigen
references is to copy matrix values when passed a numpy array that does not
conform to the element type of ``MatrixType`` or does not have a compatible
stride layout.  If you want to explicitly avoid copying in such a case, you
should bind arguments using the ``py::arg().noconvert()`` annotation (as
described in the :ref:`nonconverting_arguments` documentation).

The following example shows an example of arguments that don't allow data
copying to take place:

.. code-block:: cpp

    // The method and function to be bound:
    class MyClass {
        // ...
        double some_method(const Eigen::Ref<const MatrixXd> &matrix) { /* ... */ }
    };
    float some_function(const Eigen::Ref<const MatrixXf> &big,
                        const Eigen::Ref<const MatrixXf> &small) {
        // ...
    }

    // The associated binding code:
    using namespace pybind11::literals; // for "arg"_a
    py::class_<MyClass>(m, "MyClass")
        // ... other class definitions
        .def("some_method", &MyClass::some_method, py::arg().noconvert());

    m.def("some_function", &some_function,
        "big"_a.noconvert(), // <- Don't allow copying for this arg
        "small"_a            // <- This one can be copied if needed
    );

With the above binding code, attempting to call the the ``some_method(m)``
method on a ``MyClass`` object, or attempting to call ``some_function(m, m2)``
will raise a ``RuntimeError`` rather than making a temporary copy of the array.
It will, however, allow the ``m2`` argument to be copied into a temporary if
necessary.

Note that explicitly specifying ``.noconvert()`` is not required for *mutable*
Eigen references (e.g. ``Eigen::Ref<MatrixXd>`` without ``const`` on the
``MatrixXd``): mutable references will never be called with a temporary copy.

Vectors versus column/row matrices
==================================

Eigen and numpy have fundamentally different notions of a vector.  In Eigen, a
vector is simply a matrix with the number of columns or rows set to 1 at
compile time (for a column vector or row vector, respectively).  NumPy, in
contrast, has comparable 2-dimensional 1xN and Nx1 arrays, but *also* has
1-dimensional arrays of size N.

When passing a 2-dimensional 1xN or Nx1 array to Eigen, the Eigen type must
have matching dimensions: That is, you cannot pass a 2-dimensional Nx1 numpy
array to an Eigen value expecting a row vector, or a 1xN numpy array as a
column vector argument.

On the other hand, pybind11 allows you to pass 1-dimensional arrays of length N
as Eigen parameters.  If the Eigen type can hold a column vector of length N it
will be passed as such a column vector.  If not, but the Eigen type constraints
will accept a row vector, it will be passed as a row vector.  (The column
vector takes precedence when both are supported, for example, when passing a
1D numpy array to a MatrixXd argument).  Note that the type need not be
explicitly a vector: it is permitted to pass a 1D numpy array of size 5 to an
Eigen ``Matrix<double, Dynamic, 5>``: you would end up with a 1x5 Eigen matrix.
Passing the same to an ``Eigen::MatrixXd`` would result in a 5x1 Eigen matrix.

When returning an Eigen vector to numpy, the conversion is ambiguous: a row
vector of length 4 could be returned as either a 1D array of length 4, or as a
2D array of size 1x4.  When encountering such a situation, pybind11 compromises
by considering the returned Eigen type: if it is a compile-time vector--that
is, the type has either the number of rows or columns set to 1 at compile
time--pybind11 converts to a 1D numpy array when returning the value.  For
instances that are a vector only at run-time (e.g. ``MatrixXd``,
``Matrix<float, Dynamic, 4>``), pybind11 returns the vector as a 2D array to
numpy.  If this isn't want you want, you can use ``array.reshape(...)`` to get
a view of the same data in the desired dimensions.

.. seealso::

    The file :file:`tests/test_eigen.cpp` contains a complete example that
    shows how to pass Eigen sparse and dense data types in more detail.
Functional
##########

The following features must be enabled by including :file:`pybind11/functional.h`.


Callbacks and passing anonymous functions
=========================================

The C++11 standard brought lambda functions and the generic polymorphic
function wrapper ``std::function<>`` to the C++ programming language, which
enable powerful new ways of working with functions. Lambda functions come in
two flavors: stateless lambda function resemble classic function pointers that
link to an anonymous piece of code, while stateful lambda functions
additionally depend on captured variables that are stored in an anonymous
*lambda closure object*.

Here is a simple example of a C++ function that takes an arbitrary function
(stateful or stateless) with signature ``int -> int`` as an argument and runs
it with the value 10.

.. code-block:: cpp

    int func_arg(const std::function<int(int)> &f) {
        return f(10);
    }

The example below is more involved: it takes a function of signature ``int -> int``
and returns another function of the same kind. The return value is a stateful
lambda function, which stores the value ``f`` in the capture object and adds 1 to
its return value upon execution.

.. code-block:: cpp

    std::function<int(int)> func_ret(const std::function<int(int)> &f) {
        return [f](int i) {
            return f(i) + 1;
        };
    }

This example demonstrates using python named parameters in C++ callbacks which
requires using ``py::cpp_function`` as a wrapper. Usage is similar to defining
methods of classes:

.. code-block:: cpp

    py::cpp_function func_cpp() {
        return py::cpp_function([](int i) { return i+1; },
           py::arg("number"));
    }

After including the extra header file :file:`pybind11/functional.h`, it is almost
trivial to generate binding code for all of these functions.

.. code-block:: cpp

    #include <pybind11/functional.h>

    PYBIND11_MODULE(example, m) {
        m.def("func_arg", &func_arg);
        m.def("func_ret", &func_ret);
        m.def("func_cpp", &func_cpp);
    }

The following interactive session shows how to call them from Python.

.. code-block:: pycon

    $ python
    >>> import example
    >>> def square(i):
    ...     return i * i
    ...
    >>> example.func_arg(square)
    100L
    >>> square_plus_1 = example.func_ret(square)
    >>> square_plus_1(4)
    17L
    >>> plus_1 = func_cpp()
    >>> plus_1(number=43)
    44L

.. warning::

    Keep in mind that passing a function from C++ to Python (or vice versa)
    will instantiate a piece of wrapper code that translates function
    invocations between the two languages. Naturally, this translation
    increases the computational cost of each function call somewhat. A
    problematic situation can arise when a function is copied back and forth
    between Python and C++ many times in a row, in which case the underlying
    wrappers will accumulate correspondingly. The resulting long sequence of
    C++ -> Python -> C++ -> ... roundtrips can significantly decrease
    performance.

    There is one exception: pybind11 detects case where a stateless function
    (i.e. a function pointer or a lambda function without captured variables)
    is passed as an argument to another C++ function exposed in Python. In this
    case, there is no overhead. Pybind11 will extract the underlying C++
    function pointer from the wrapped function to sidestep a potential C++ ->
    Python -> C++ roundtrip. This is demonstrated in :file:`tests/test_callbacks.cpp`.

.. note::

    This functionality is very useful when generating bindings for callbacks in
    C++ libraries (e.g. GUI libraries, asynchronous networking libraries, etc.).

    The file :file:`tests/test_callbacks.cpp` contains a complete example
    that demonstrates how to work with callbacks and anonymous functions in
    more detail.
Strings, bytes and Unicode conversions
######################################

.. note::

    This section discusses string handling in terms of Python 3 strings. For
    Python 2.7, replace all occurrences of ``str`` with ``unicode`` and
    ``bytes`` with ``str``.  Python 2.7 users may find it best to use ``from
    __future__ import unicode_literals`` to avoid unintentionally using ``str``
    instead of ``unicode``.

Passing Python strings to C++
=============================

When a Python ``str`` is passed from Python to a C++ function that accepts
``std::string`` or ``char *`` as arguments, pybind11 will encode the Python
string to UTF-8. All Python ``str`` can be encoded in UTF-8, so this operation
does not fail.

The C++ language is encoding agnostic. It is the responsibility of the
programmer to track encodings. It's often easiest to simply `use UTF-8
everywhere <http://utf8everywhere.org/>`_.

.. code-block:: c++

    m.def("utf8_test",
        [](const std::string &s) {
            cout << "utf-8 is icing on the cake.\n";
            cout << s;
        }
    );
    m.def("utf8_charptr",
        [](const char *s) {
            cout << "My favorite food is\n";
            cout << s;
        }
    );

.. code-block:: python

    >>> utf8_test('🎂')
    utf-8 is icing on the cake.
    🎂

    >>> utf8_charptr('🍕')
    My favorite food is
    🍕

.. note::

    Some terminal emulators do not support UTF-8 or emoji fonts and may not
    display the example above correctly.

The results are the same whether the C++ function accepts arguments by value or
reference, and whether or not ``const`` is used.

Passing bytes to C++
--------------------

A Python ``bytes`` object will be passed to C++ functions that accept
``std::string`` or ``char*`` *without* conversion.  On Python 3, in order to
make a function *only* accept ``bytes`` (and not ``str``), declare it as taking
a ``py::bytes`` argument.


Returning C++ strings to Python
===============================

When a C++ function returns a ``std::string`` or ``char*`` to a Python caller,
**pybind11 will assume that the string is valid UTF-8** and will decode it to a
native Python ``str``, using the same API as Python uses to perform
``bytes.decode('utf-8')``. If this implicit conversion fails, pybind11 will
raise a ``UnicodeDecodeError``.

.. code-block:: c++

    m.def("std_string_return",
        []() {
            return std::string("This string needs to be UTF-8 encoded");
        }
    );

.. code-block:: python

    >>> isinstance(example.std_string_return(), str)
    True


Because UTF-8 is inclusive of pure ASCII, there is never any issue with
returning a pure ASCII string to Python. If there is any possibility that the
string is not pure ASCII, it is necessary to ensure the encoding is valid
UTF-8.

.. warning::

    Implicit conversion assumes that a returned ``char *`` is null-terminated.
    If there is no null terminator a buffer overrun will occur.

Explicit conversions
--------------------

If some C++ code constructs a ``std::string`` that is not a UTF-8 string, one
can perform a explicit conversion and return a ``py::str`` object. Explicit
conversion has the same overhead as implicit conversion.

.. code-block:: c++

    // This uses the Python C API to convert Latin-1 to Unicode
    m.def("str_output",
        []() {
            std::string s = "Send your r\xe9sum\xe9 to Alice in HR"; // Latin-1
            py::str py_s = PyUnicode_DecodeLatin1(s.data(), s.length());
            return py_s;
        }
    );

.. code-block:: python

    >>> str_output()
    'Send your résumé to Alice in HR'

The `Python C API
<https://docs.python.org/3/c-api/unicode.html#built-in-codecs>`_ provides
several built-in codecs.


One could also use a third party encoding library such as libiconv to transcode
to UTF-8.

Return C++ strings without conversion
-------------------------------------

If the data in a C++ ``std::string`` does not represent text and should be
returned to Python as ``bytes``, then one can return the data as a
``py::bytes`` object.

.. code-block:: c++

    m.def("return_bytes",
        []() {
            std::string s("\xba\xd0\xba\xd0");  // Not valid UTF-8
            return py::bytes(s);  // Return the data without transcoding
        }
    );

.. code-block:: python

    >>> example.return_bytes()
    b'\xba\xd0\xba\xd0'


Note the asymmetry: pybind11 will convert ``bytes`` to ``std::string`` without
encoding, but cannot convert ``std::string`` back to ``bytes`` implicitly.

.. code-block:: c++

    m.def("asymmetry",
        [](std::string s) {  // Accepts str or bytes from Python
            return s;  // Looks harmless, but implicitly converts to str
        }
    );

.. code-block:: python

    >>> isinstance(example.asymmetry(b"have some bytes"), str)
    True

    >>> example.asymmetry(b"\xba\xd0\xba\xd0")  # invalid utf-8 as bytes
    UnicodeDecodeError: 'utf-8' codec can't decode byte 0xba in position 0: invalid start byte


Wide character strings
======================

When a Python ``str`` is passed to a C++ function expecting ``std::wstring``,
``wchar_t*``, ``std::u16string`` or ``std::u32string``, the ``str`` will be
encoded to UTF-16 or UTF-32 depending on how the C++ compiler implements each
type, in the platform's native endianness. When strings of these types are
returned, they are assumed to contain valid UTF-16 or UTF-32, and will be
decoded to Python ``str``.

.. code-block:: c++

    #define UNICODE
    #include <windows.h>

    m.def("set_window_text",
        [](HWND hwnd, std::wstring s) {
            // Call SetWindowText with null-terminated UTF-16 string
            ::SetWindowText(hwnd, s.c_str());
        }
    );
    m.def("get_window_text",
        [](HWND hwnd) {
            const int buffer_size = ::GetWindowTextLength(hwnd) + 1;
            auto buffer = std::make_unique< wchar_t[] >(buffer_size);

            ::GetWindowText(hwnd, buffer.data(), buffer_size);

            std::wstring text(buffer.get());

            // wstring will be converted to Python str
            return text;
        }
    );

.. warning::

    Wide character strings may not work as described on Python 2.7 or Python
    3.3 compiled with ``--enable-unicode=ucs2``.

Strings in multibyte encodings such as Shift-JIS must transcoded to a
UTF-8/16/32 before being returned to Python.


Character literals
==================

C++ functions that accept character literals as input will receive the first
character of a Python ``str`` as their input. If the string is longer than one
Unicode character, trailing characters will be ignored.

When a character literal is returned from C++ (such as a ``char`` or a
``wchar_t``), it will be converted to a ``str`` that represents the single
character.

.. code-block:: c++

    m.def("pass_char", [](char c) { return c; });
    m.def("pass_wchar", [](wchar_t w) { return w; });

.. code-block:: python

    >>> example.pass_char('A')
    'A'

While C++ will cast integers to character types (``char c = 0x65;``), pybind11
does not convert Python integers to characters implicitly. The Python function
``chr()`` can be used to convert integers to characters.

.. code-block:: python

    >>> example.pass_char(0x65)
    TypeError

    >>> example.pass_char(chr(0x65))
    'A'

If the desire is to work with an 8-bit integer, use ``int8_t`` or ``uint8_t``
as the argument type.

Grapheme clusters
-----------------

A single grapheme may be represented by two or more Unicode characters. For
example 'é' is usually represented as U+00E9 but can also be expressed as the
combining character sequence U+0065 U+0301 (that is, the letter 'e' followed by
a combining acute accent). The combining character will be lost if the
two-character sequence is passed as an argument, even though it renders as a
single grapheme.

.. code-block:: python

    >>> example.pass_wchar('é')
    'é'

    >>> combining_e_acute = 'e' + '\u0301'

    >>> combining_e_acute
    'é'

    >>> combining_e_acute == 'é'
    False

    >>> example.pass_wchar(combining_e_acute)
    'e'

Normalizing combining characters before passing the character literal to C++
may resolve *some* of these issues:

.. code-block:: python

    >>> example.pass_wchar(unicodedata.normalize('NFC', combining_e_acute))
    'é'

In some languages (Thai for example), there are `graphemes that cannot be
expressed as a single Unicode code point
<http://unicode.org/reports/tr29/#Grapheme_Cluster_Boundaries>`_, so there is
no way to capture them in a C++ character type.


C++17 string views
==================

C++17 string views are automatically supported when compiling in C++17 mode.
They follow the same rules for encoding and decoding as the corresponding STL
string type (for example, a ``std::u16string_view`` argument will be passed
UTF-16-encoded data, and a returned ``std::string_view`` will be decoded as
UTF-8).

References
==========

* `The Absolute Minimum Every Software Developer Absolutely, Positively Must Know About Unicode and Character Sets (No Excuses!) <https://www.joelonsoftware.com/2003/10/08/the-absolute-minimum-every-software-developer-absolutely-positively-must-know-about-unicode-and-character-sets-no-excuses/>`_
* `C++ - Using STL Strings at Win32 API Boundaries <https://msdn.microsoft.com/en-ca/magazine/mt238407.aspx>`_
Custom type casters
===================

In very rare cases, applications may require custom type casters that cannot be
expressed using the abstractions provided by pybind11, thus requiring raw
Python C API calls. This is fairly advanced usage and should only be pursued by
experts who are familiar with the intricacies of Python reference counting.

The following snippets demonstrate how this works for a very simple ``inty``
type that that should be convertible from Python types that provide a
``__int__(self)`` method.

.. code-block:: cpp

    struct inty { long long_value; };

    void print(inty s) {
        std::cout << s.long_value << std::endl;
    }

The following Python snippet demonstrates the intended usage from the Python side:

.. code-block:: python

    class A:
        def __int__(self):
            return 123

    from example import print
    print(A())

To register the necessary conversion routines, it is necessary to add an
instantiation of the ``pybind11::detail::type_caster<T>`` template.
Although this is an implementation detail, adding an instantiation of this
type is explicitly allowed.

.. code-block:: cpp

    namespace pybind11 { namespace detail {
        template <> struct type_caster<inty> {
        public:
            /**
             * This macro establishes the name 'inty' in
             * function signatures and declares a local variable
             * 'value' of type inty
             */
            PYBIND11_TYPE_CASTER(inty, _("inty"));

            /**
             * Conversion part 1 (Python->C++): convert a PyObject into a inty
             * instance or return false upon failure. The second argument
             * indicates whether implicit conversions should be applied.
             */
            bool load(handle src, bool) {
                /* Extract PyObject from handle */
                PyObject *source = src.ptr();
                /* Try converting into a Python integer value */
                PyObject *tmp = PyNumber_Long(source);
                if (!tmp)
                    return false;
                /* Now try to convert into a C++ int */
                value.long_value = PyLong_AsLong(tmp);
                Py_DECREF(tmp);
                /* Ensure return code was OK (to avoid out-of-range errors etc) */
                return !(value.long_value == -1 && !PyErr_Occurred());
            }

            /**
             * Conversion part 2 (C++ -> Python): convert an inty instance into
             * a Python object. The second and third arguments are used to
             * indicate the return value policy and parent object (for
             * ``return_value_policy::reference_internal``) and are generally
             * ignored by implicit casters.
             */
            static handle cast(inty src, return_value_policy /* policy */, handle /* parent */) {
                return PyLong_FromLong(src.long_value);
            }
        };
    }} // namespace pybind11::detail

.. note::

    A ``type_caster<T>`` defined with ``PYBIND11_TYPE_CASTER(T, ...)`` requires
    that ``T`` is default-constructible (``value`` is first default constructed
    and then ``load()`` assigns to it).

.. warning::

    When using custom type casters, it's important to declare them consistently
    in every compilation unit of the Python extension module. Otherwise,
    undefined behavior can ensue.
Chrono
======

When including the additional header file :file:`pybind11/chrono.h` conversions
from C++11 chrono datatypes to python datetime objects are automatically enabled.
This header also enables conversions of python floats (often from sources such
as ``time.monotonic()``, ``time.perf_counter()`` and ``time.process_time()``)
into durations.

An overview of clocks in C++11
------------------------------

A point of confusion when using these conversions is the differences between
clocks provided in C++11. There are three clock types defined by the C++11
standard and users can define their own if needed. Each of these clocks have
different properties and when converting to and from python will give different
results.

The first clock defined by the standard is ``std::chrono::system_clock``. This
clock measures the current date and time. However, this clock changes with to
updates to the operating system time. For example, if your time is synchronised
with a time server this clock will change. This makes this clock a poor choice
for timing purposes but good for measuring the wall time.

The second clock defined in the standard is ``std::chrono::steady_clock``.
This clock ticks at a steady rate and is never adjusted. This makes it excellent
for timing purposes, however the value in this clock does not correspond to the
current date and time. Often this clock will be the amount of time your system
has been on, although it does not have to be. This clock will never be the same
clock as the system clock as the system clock can change but steady clocks
cannot.

The third clock defined in the standard is ``std::chrono::high_resolution_clock``.
This clock is the clock that has the highest resolution out of the clocks in the
system. It is normally a typedef to either the system clock or the steady clock
but can be its own independent clock. This is important as when using these
conversions as the types you get in python for this clock might be different
depending on the system.
If it is a typedef of the system clock, python will get datetime objects, but if
it is a different clock they will be timedelta objects.

Provided conversions
--------------------

.. rubric:: C++ to Python

- ``std::chrono::system_clock::time_point`` → ``datetime.datetime``
    System clock times are converted to python datetime instances. They are
    in the local timezone, but do not have any timezone information attached
    to them (they are naive datetime objects).

- ``std::chrono::duration`` → ``datetime.timedelta``
    Durations are converted to timedeltas, any precision in the duration
    greater than microseconds is lost by rounding towards zero.

- ``std::chrono::[other_clocks]::time_point`` → ``datetime.timedelta``
    Any clock time that is not the system clock is converted to a time delta.
    This timedelta measures the time from the clocks epoch to now.

.. rubric:: Python to C++

- ``datetime.datetime`` or ``datetime.date`` or ``datetime.time`` → ``std::chrono::system_clock::time_point``
    Date/time objects are converted into system clock timepoints. Any
    timezone information is ignored and the type is treated as a naive
    object.

- ``datetime.timedelta`` → ``std::chrono::duration``
    Time delta are converted into durations with microsecond precision.

- ``datetime.timedelta`` → ``std::chrono::[other_clocks]::time_point``
    Time deltas that are converted into clock timepoints are treated as
    the amount of time from the start of the clocks epoch.

- ``float`` → ``std::chrono::duration``
    Floats that are passed to C++ as durations be interpreted as a number of
    seconds. These will be converted to the duration using ``duration_cast``
    from the float.

- ``float`` → ``std::chrono::[other_clocks]::time_point``
    Floats that are passed to C++ as time points will be interpreted as the
    number of seconds from the start of the clocks epoch.
Overview
########

.. rubric:: 1. Native type in C++, wrapper in Python

Exposing a custom C++ type using :class:`py::class_` was covered in detail
in the :doc:`/classes` section. There, the underlying data structure is
always the original C++ class while the :class:`py::class_` wrapper provides
a Python interface. Internally, when an object like this is sent from C++ to
Python, pybind11 will just add the outer wrapper layer over the native C++
object. Getting it back from Python is just a matter of peeling off the
wrapper.

.. rubric:: 2. Wrapper in C++, native type in Python

This is the exact opposite situation. Now, we have a type which is native to
Python, like a ``tuple`` or a ``list``. One way to get this data into C++ is
with the :class:`py::object` family of wrappers. These are explained in more
detail in the :doc:`/advanced/pycpp/object` section. We'll just give a quick
example here:

.. code-block:: cpp

    void print_list(py::list my_list) {
        for (auto item : my_list)
            std::cout << item << " ";
    }

.. code-block:: pycon

    >>> print_list([1, 2, 3])
    1 2 3

The Python ``list`` is not converted in any way -- it's just wrapped in a C++
:class:`py::list` class. At its core it's still a Python object. Copying a
:class:`py::list` will do the usual reference-counting like in Python.
Returning the object to Python will just remove the thin wrapper.

.. rubric:: 3. Converting between native C++ and Python types

In the previous two cases we had a native type in one language and a wrapper in
the other. Now, we have native types on both sides and we convert between them.

.. code-block:: cpp

    void print_vector(const std::vector<int> &v) {
        for (auto item : v)
            std::cout << item << "\n";
    }

.. code-block:: pycon

    >>> print_vector([1, 2, 3])
    1 2 3

In this case, pybind11 will construct a new ``std::vector<int>`` and copy each
element from the Python ``list``. The newly constructed object will be passed
to ``print_vector``. The same thing happens in the other direction: a new
``list`` is made to match the value returned from C++.

Lots of these conversions are supported out of the box, as shown in the table
below. They are very convenient, but keep in mind that these conversions are
fundamentally based on copying data. This is perfectly fine for small immutable
types but it may become quite expensive for large data structures. This can be
avoided by overriding the automatic conversion with a custom wrapper (i.e. the
above-mentioned approach 1). This requires some manual effort and more details
are available in the :ref:`opaque` section.

.. _conversion_table:

List of all builtin conversions
-------------------------------

The following basic data types are supported out of the box (some may require
an additional extension header to be included). To pass other data structures
as arguments and return values, refer to the section on binding :ref:`classes`.

+------------------------------------+---------------------------+-------------------------------+
|  Data type                         |  Description              | Header file                   |
+====================================+===========================+===============================+
| ``int8_t``, ``uint8_t``            | 8-bit integers            | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``int16_t``, ``uint16_t``          | 16-bit integers           | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``int32_t``, ``uint32_t``          | 32-bit integers           | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``int64_t``, ``uint64_t``          | 64-bit integers           | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``ssize_t``, ``size_t``            | Platform-dependent size   | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``float``, ``double``              | Floating point types      | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``bool``                           | Two-state Boolean type    | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``char``                           | Character literal         | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``char16_t``                       | UTF-16 character literal  | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``char32_t``                       | UTF-32 character literal  | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``wchar_t``                        | Wide character literal    | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``const char *``                   | UTF-8 string literal      | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``const char16_t *``               | UTF-16 string literal     | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``const char32_t *``               | UTF-32 string literal     | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``const wchar_t *``                | Wide string literal       | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``std::string``                    | STL dynamic UTF-8 string  | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``std::u16string``                 | STL dynamic UTF-16 string | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``std::u32string``                 | STL dynamic UTF-32 string | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``std::wstring``                   | STL dynamic wide string   | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``std::string_view``,              | STL C++17 string views    | :file:`pybind11/pybind11.h`   |
| ``std::u16string_view``, etc.      |                           |                               |
+------------------------------------+---------------------------+-------------------------------+
| ``std::pair<T1, T2>``              | Pair of two custom types  | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``std::tuple<...>``                | Arbitrary tuple of types  | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``std::reference_wrapper<...>``    | Reference type wrapper    | :file:`pybind11/pybind11.h`   |
+------------------------------------+---------------------------+-------------------------------+
| ``std::complex<T>``                | Complex numbers           | :file:`pybind11/complex.h`    |
+------------------------------------+---------------------------+-------------------------------+
| ``std::array<T, Size>``            | STL static array          | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::vector<T>``                 | STL dynamic array         | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::deque<T>``                  | STL double-ended queue    | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::valarray<T>``               | STL value array           | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::list<T>``                   | STL linked list           | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::map<T1, T2>``               | STL ordered map           | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::unordered_map<T1, T2>``     | STL unordered map         | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::set<T>``                    | STL ordered set           | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::unordered_set<T>``          | STL unordered set         | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::optional<T>``               | STL optional type (C++17) | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::experimental::optional<T>`` | STL optional type (exp.)  | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::variant<...>``              | Type-safe union (C++17)   | :file:`pybind11/stl.h`        |
+------------------------------------+---------------------------+-------------------------------+
| ``std::function<...>``             | STL polymorphic function  | :file:`pybind11/functional.h` |
+------------------------------------+---------------------------+-------------------------------+
| ``std::chrono::duration<...>``     | STL time duration         | :file:`pybind11/chrono.h`     |
+------------------------------------+---------------------------+-------------------------------+
| ``std::chrono::time_point<...>``   | STL date/time             | :file:`pybind11/chrono.h`     |
+------------------------------------+---------------------------+-------------------------------+
| ``Eigen::Matrix<...>``             | Eigen: dense matrix       | :file:`pybind11/eigen.h`      |
+------------------------------------+---------------------------+-------------------------------+
| ``Eigen::Map<...>``                | Eigen: mapped memory      | :file:`pybind11/eigen.h`      |
+------------------------------------+---------------------------+-------------------------------+
| ``Eigen::SparseMatrix<...>``       | Eigen: sparse matrix      | :file:`pybind11/eigen.h`      |
+------------------------------------+---------------------------+-------------------------------+
.. _type-conversions:

Type conversions
################

Apart from enabling cross-language function calls, a fundamental problem
that a binding tool like pybind11 must address is to provide access to
native Python types in C++ and vice versa. There are three fundamentally
different ways to do this—which approach is preferable for a particular type
depends on the situation at hand.

1. Use a native C++ type everywhere. In this case, the type must be wrapped
   using pybind11-generated bindings so that Python can interact with it.

2. Use a native Python type everywhere. It will need to be wrapped so that
   C++ functions can interact with it.

3. Use a native C++ type on the C++ side and a native Python type on the
   Python side. pybind11 refers to this as a *type conversion*.

   Type conversions are the most "natural" option in the sense that native
   (non-wrapped) types are used everywhere. The main downside is that a copy
   of the data must be made on every Python ↔ C++ transition: this is
   needed since the C++ and Python versions of the same type generally won't
   have the same memory layout.

   pybind11 can perform many kinds of conversions automatically. An overview
   is provided in the table ":ref:`conversion_table`".

The following subsections discuss the differences between these options in more
detail. The main focus in this section is on type conversions, which represent
the last case of the above list.

.. toctree::
   :maxdepth: 1

   overview
   strings
   stl
   functional
   chrono
   eigen
   custom
CMake helpers
-------------

Pybind11 can be used with ``add_subdirectory(extern/pybind11)``, or from an
install with ``find_package(pybind11 CONFIG)``. The interface provided in
either case is functionally identical.

.. cmake-module:: ../../tools/pybind11Config.cmake.in
