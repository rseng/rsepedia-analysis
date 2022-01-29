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
