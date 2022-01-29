The easiest way to build and install the ParaMonte library and its examples is via the 
[Linux/macOS Bash](https://github.com/cdslaborg/paramonte/blob/main/install.sh.usage.txt) and 
[Windows Batch](https://github.com/cdslaborg/paramonte/blob/main/install.bat.usage.txt) 
install scripts that exist at the root of the ParaMonte repository on GitHub. 
The scripts will automatically install the library and run the tests (if requested) and the examples provided.  

However, if there are still good reasons to build library directly via the cmake software 
and manually go through all the rest of the build, testing, and installation processes, 
then you can follow the guidelines below to predefine the various settings of the library and its tests.  

**Note** that the ParaMonte library build via cmake is currently exclusively available and tested on Linux, macOS, and Microsoft WSL platforms.  

1.  Ensure a [ParaMonte-compatible](CHANGES.md) version of the GNU or Intel compiler suite 
    is already installed on the system and the paths to the different components of the 
    compiler and wrappers already exist in the environmental `PATH` variable. 
    (**This step would be automatically done via the [`install.sh` script](https://github.com/cdslaborg/paramonte/blob/main/install.sh.usage.txt).**)  
1.  Ensure a ParaMonte-compatible cmake installation (`3.14.0` or newer) exists on your system 
    and the path to the cmake executable already exist in the environmental `PATH` variable. 
    (**This step would be automatically done via the [`install.sh` script](https://github.com/cdslaborg/paramonte/blob/main/install.sh.usage.txt).**)  
1.  Navigate to the root directory of the library (where the `.git` file exists).  
1.  Then, assuming that you are using the GNU compiler suite, you can configure the library's build minimally via,  
    ```bash  
    mkdir build && \
    cd build && \
    cmake -DPMCS=GNU ..
    ```  
    (**This step would be automatically done via the [`install.sh` script](https://github.com/cdslaborg/paramonte/blob/main/install.sh.usage.txt).**)  
1.  The above command will configure the library with the default build settings. 
    To change the default settings, the following macros can be also defined,  
    ```bash  
    mkdir build && \
    cd build && \
    cmake \
    --verbose=1 \
    -DPMCS="the compiler suite: GNU/INTEL" \
    -DPMLIB_NAME="you_must_follow_the_guidelines_here:https://www.cdslab.org/paramonte/notes/installation/readme/#naming-convention-used-for-paramonte-library-builds" \
    -DCMAKE_Fortran_COMPILER="/path/to/the/specific/fortran/compiler/executable/to/be/used/example:/usr/bin/gfortran-10" \
    -DMATLAB_ROOT_DIR="/if/building/for/matlab/provide/path/to/the/matlab/root/installation/directory/where/bin/dir/exists" \
    -DMPIEXEC_EXECUTABLE="/path/to/the/specific/mpi/launcher/to/be/used" \
    -DINTERFACE_LANGUAGE="the single programming language interface: c/cpp/fortran/matlab/python/... (default = c)" \
    -DMPI_ENABLED="whether building for MPI parallelization: true/false (default = false)" \
    -DCAFTYPE="The Coarray parallelization type: none/single/shared/distributed (default = none)" \
    -DBTYPE="library build: debug/testing/release (default = release)" \
    -DLTYPE="library type: static/shared (default = static)" \
    -DHEAP_ARRAY_ENABLED="whether to perform all memory allocations on the heap: true/false (default = false)" \
    -DCFI_ENABLED="whether C-Fortran interoperation must be enabled: true/false (default = true)" \
    -DOMP_ENABLED="whether OpenMP parallelization must be enabled: true/false (default = false)" \
    -DMPILIB_NAME="The MPI library name: impi/mpich/openmpi (impi stands for Intel MPI)" \
    -DSAMPLER_TEST_ENABLED="whether the library's sampler tests should be enabled: true/false (default = false)" \
    -DBASIC_TEST_ENABLED="whether the library's basic tests should be enabled: true/false (default = false)" \
    -DPERFPROF_ENABLED="whether the library should be build for performance profiling: true/false (default = false)" \
    -DCODECOV_ENABLED="whether the library should be build for code coverage: true/false (default = false)" \
    -DOS_IS_WSL="whether the current operating system is Microsoft WSL" \
    ..
    ```  
    (**All of the above configuration settings could be easily done via the [`install.sh` script](https://github.com/cdslaborg/paramonte/blob/main/install.sh.usage.txt).**)  
1.  Once you configure the library's build, try,  
    ```bash  
    make -j 3 && make install
    ```  
    Replace the number `3` with as many number of cores you wish to use to build the library. 
    (**This step would be automatically done via the [`install.sh` script](https://github.com/cdslaborg/paramonte/blob/main/install.sh.usage.txt).**)  
1.  The ParaMonte library should be now built and installed in the `lib` directory. 
    If you built the library for testing or code coverage purposes, then a `testParaMonte` 
    executable must also exist in the subdirectory `test/bin/`. The onus is on the developer/user (presumably you) 
    to ensure all runtime libraries and dependencies can be found by the library or the test executable or the examples.  
    (**All dependency checks and test runs would be automatically done via the [`install.sh` script](https://github.com/cdslaborg/paramonte/blob/main/install.sh.usage.txt).**)  
  
<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank"><img src="https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png" alt="ParaMonte: Plain Powerful Parallel Monte Carlo Library" /></a>
<br><br>
<a href="https://github.com/cdslaborg/paramonte#license" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://travis-ci.com/cdslaborg/paramonte" target="_blank"><img src="https://travis-ci.com/cdslaborg/paramonte.svg?branch=main&style=flat-square" alt="Build Status" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/serial/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-serial%20kernel%20:%2098.2%25-brightgreen?style=flat-square" alt="kernel code coverage - serial" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/mpi/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-MPI%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - MPI" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/caf/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-Coarray%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - Coarray" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=PyPI%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/badge/dynamic/json?style=flat-square&labelColor=grey&color=brightgreen&maxAge=86400&label=PyPI%20downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fparamonte" alt="PyPI - Downloads Total" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
<br><br>
<a href="#paramonte-plain-powerful-parallel-monte-carlo-library">Overview</a> | 
<a href="#installation">Installation</a> | 
<a href="#dependencies">Dependencies</a> | 
<a href="#parallelism">Parallelism</a> | 
<a href="#example-usage-instructions">Examples</a> |
<a href="#citing-paramonte">Acknowledgments</a> | 
<a href="#license">License</a> | 
<a href="#authors-and-contributors">Authors</a>  
</div>
  
  
ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions 
of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, 
and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), 
**user-friendliness** (of the library), **accessibility** (from multiple programming environments), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  
  
  
ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  

-   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray 
    communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations up to 16 digits of precision as requested by the user.  

-   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in 
    external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

  
Installation  
============  

The pre-built ready-to-use libraries are available on [the release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases). 
Each prebuilt ParaMonte library automatically ships with a full-fledged set of example codes and build scripts.  

Alternatively, you can build the library from the source in the [GitHub repository of the project](https://github.com/cdslaborg/paramonte). 
The ParaMonte library installation/build process is fully automated for all of the supported programming languages. 
Currently, the following compiler suites are supported **for builds from source**:  
  
| Compiler Suite                    | Linux | macOS | Windows (64bit) |  
|----------------------------------:|:-----:|:-----:|:---------------:|  
| GNU Compiler Collection > 8.4     |&check;|&check;| &cross;         |  
| Intel Parallel Studio > 19.1.1    |&check;|&check;| &check;         |  

For more information and quick-start in the programming language of your choice, visit the [ParaMonte library homepage](https://www.cdslab.org/paramonte).  

  
Dependencies  
============  

Beyond an optional MPI runtime library for parallel simulations, the ParaMonte kernel has **zero dependency** on external third-party libraries or packages.  

  
Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. 
To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. 
In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). 
These automatic checks and installations happen when you download and install or use the library on your system, for the first time. 
If the automatic installation is unsuccessful, you can also install the libraries manually on your system:  

+   On **Windows** and **Linux** operating systems, we highly recommend downloading and installing the 
    [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), 
    which is available to the public free of charge, also available in the latest release of the 
    ParaMonte library on the [GitHub release page](https://github.com/cdslaborg/paramonte/releases) 
    (For Windows, look for the executable file that ends with `.exe`. For Linux, look for the file 
    that ends with `.tgz`, like `l_mpi-rt_2018.2.199.tgz`).
+   On **macOS**, the Intel MPI library is not available. Therefore, we recommend installing either 
    [Open-MPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/) MPI runtime libraries 
    depending the prebuilt version of the ParaMonte library that you have downloaded or 
    the configuration with which you intend to build the library.  

For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

  
Example usage instructions  
==========================  

+   For complete organized up-to-date instructions, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

+   For a quick look into *language-specific* README.md instructions, visit:  
    +   **C**: [https://github.com/cdslaborg/paramonte/tree/main/src/interface/C](https://github.com/cdslaborg/paramonte/tree/main/src/interface/C)  
    +   **C++**: [https://github.com/cdslaborg/paramonte/tree/main/src/interface/C++](https://github.com/cdslaborg/paramonte/tree/main/src/interface/C++)  
    +   **Fortran**: [https://github.com/cdslaborg/paramonte/tree/main/src/interface/Fortran](https://github.com/cdslaborg/paramonte/tree/main/src/interface/Fortran)  
    +   **MATLAB**: [https://github.com/cdslaborg/paramonte/tree/main/src/interface/MATLAB](https://github.com/cdslaborg/paramonte/tree/main/src/interface/MATLAB)  
    +   **Python**: [https://github.com/cdslaborg/paramonte/tree/main/src/interface/Python](https://github.com/cdslaborg/paramonte/tree/main/src/interface/Python)  

  
Citing ParaMonte  
================  

The ParaMonte library is an honor-ware, the currency of which is acknowledgment and citations.  
  
If you use ParaMonte or any ideas from the software, please acknowledge it by citing the ParaMonte library's 
main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) 
to access the PDF version of these files free of charge.  

  
License  
=======  

[MIT License](https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md)  

**What does this license mean?**  

Essentially, all we are asking from the users or developers is to  

>   explicitly acknowledge the use of this library or any concepts or parts of it in their education, research, or software (free or commercial).  

This is a free software, so help us keep it freely available to the public by redistributing the library and contributing to it. 
If you have questions or concerns about the license, do not hesitate to contact us (shahmoradi@utexas.edu).  

  
Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, 
        data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, 
        high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, 
        biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, 
        protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+   [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)  
    +   physicist / cosmologist by training,  
    +   currently a UTA Physics member,  
    +   deep philosophical thinker,  
    +   contact: [Fatemeh.Bagheri@uta.edu](mailto:"Fatemeh.Bagheri@uta.edu")  

+   [Shashank Kumbhare](https://www.cdslab.org/people/#shashank-kumbhare)  
    +   physicist / Computational Data Scientist,  
    +   currently a UTA Physics member,  
    +   contact: [shashankkumbhare8@gmail.com](mailto:"shashankkumbhare8@gmail.com")  

+   [Joshua Osborne](https://www.cdslab.org/people/#joshua-alexander-osborne)  
    +   physicist / Computational Data Scientist by training,  
    +   currently a UTA Physics member,  
    +   contact: [joshuaalexanderosborne@gmail.com](mailto:"joshuaalexanderosborne@gmail.com")  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

Thank you for your interest in contributing to ParaMonte! Your help is very much appreciated. 
Below are some tips and guidelines to get started with **contributing to the ParaMonte kernel routines**.  

## Initial steps  

+   First, read the [general development guidelines](CONTRIBUTING.md). 
+   Take a look at the [issues](https://github.com/cdslaborg/paramonte/issues) page. 
    Make sure you find an open issue **about the kernel routines** and that you do not duplicate someone else's work.  
+   If your contribution does not exist as an issue, post a [new issue](https://github.com/cdslaborg/paramonte/issues/new/choose) discussing the changes you're proposing to implement, 
    whether bug fix(es) or enhancement/feature request(s), or give the rest of developers a heads up that you are going to start work on [an open issue](https://github.com/cdslaborg/paramonte/issues).  
+   [Fork the ParaMonte git project](https://help.github.com/articles/fork-a-repo/) to your private account.  

## Kernel development conventions  

Pay careful attention to the following conventions used in the development of the kernel routines.

### Preprocessor directives  

The kernel routines of ParaMonte heavily rely on the compiler preprocessing directives to implement multiple completely different parallelism paradigms (serial, MPI, Coarray) that are specifically tailored for the needs of multiple completely different programming languages (C, C++, Fortran, Julia, MATLAB, Python, R, ...) for different builds (`debug`, `testing`, `release`), for different operating systems (`Windows`, `Linux`, `macOS`, `WSL`), and many more.  

Despite significant efforts to minimize the use of preprocessor directive in the kernel routines, their frequent usage is unavoidable as they greatly minimize the cost of development. Simultaneously, however, the preprocessor directive can increase confusion and and lead to implicit bugs that are hard to detect.  

> **WARNING**  
> Any code section that requires some specific operating system, platform, compiler, programming language, compiler settings, or parallelism paradigms must be fenced with the appropriate preprocessor macros such that it is executed only in the appropriate circumstances. Failure to do so, will lead to a non-portable codebase. In general, there must exist an `#else` preprocessor section for any fenced code section that appropriately handles other possible scenarios, if the other scenarios are not possible, it is followed by an `#error \"error message\"` preprocessor directive.  

> **IMPORTANT**  
> The preprocessor macros must be always defined with uppercase characters. This is to easily distinguish macros from other variables in the regular sections of the code.  

> **IMPORTANT**  
> Define macros as propositions that evaluate to either TRUE or FALSE. For example, `OS_IS_WIDOWS` or `MPI_ENABLED`.  

> **TIP**  
> When developing or debugging a particular feature that is bound to a specific set of preprocessor macros, you can request the compiler to generate intermediate preprocessed files that are free from preprocessor directives and let you focus on the problem better. In doing so, keep in mind that every change you make to the intermediate codes has be carefully implemented in the original set of codes.  

#### Platform preprocessor directives  

If there is any section of the code that must be exclusively executed on a particular platform or Operating System (OS), it must be properly fenced by the appropriate preprocessor flag.

1.  The `OS_IS_WINDOWS` preprocessor macro must be defined for any code section that exclusively belongs to Windows operating systems.  
1.  The `OS_IS_DARWIN` preprocessor macro must be defined for any code section that exclusively belongs to macOS operating systems.  
1.  The `OS_IS_LINUX` preprocessor macro must be defined for any code section that exclusively belongs to Linux operating systems.  
1.  The `OS_IS_WSL` preprocessor macro must be defined for any code section that exclusively belongs to Microsoft Subsystem for Linux.  
    > **NOTE**  
    > The `OS_IS_WSL` preprocessor macro is mostly useful for defining code sections that are written for compatibility with Microsoft Subsystem for Linux Version 1 (**WSL1**). For example, WSL1 is known to have issues with passing internal procedures to external routines since [WSL has a non-executable stack](https://github.com/Microsoft/WSL/issues/3083). Although, [this issue](https://github.com/Microsoft/WSL/issues/286) appears to have been resolved in [WSL2](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux#WSL_2), as a developer, you should be mindful of all users on all platforms.  

Several issues also need special care when developing ParaMonte on Windows OS:  

+   On Windows systems, **file-locking** by a processor often leads to problems and program crashes in parallel simulations. Since, currently only Intel compilers are supported and tested with the ParaMonte library, a quick remedy is activate the shared file manipulation via the Intel compiler `SHARED` argument that can be passed to `open()` statements. For example,  
    ```fortran  
    open( unit      = self%SampleFile%unit            &
        , file      = self%SampleFile%Path%original   &
        , status    = self%SampleFile%status          &
        , iostat    = self%SampleFile%Err%stat        &
    #if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
        , SHARED                                      &
    #endif
        , position  = self%SampleFile%Position%value  )
    ```  
+   To properly build shared (DLL) libraries on Windows OS, every procedure name must be correctly extracted via Intel compiler `DLLEXPORT` directive. Currently, all names are manually extracted immediately below each procedure's interface. This is normally not an issue, as long the symbols are correctly exported when the procedure is being developed. However, this can quickly become a challenge if the developer uses a wrong symbol name for the procedure or completely forgets to export the symbol.  
    > **NOTE**  
    > A potential solution to this problem is to enable builds with CMAKE on Windows and [request CMAKE to automatically export all symbols](https://blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/). Indeed, extending the current cmake build system of the ParaMonte library to Windows would be a great contribution to the package, if you can accomplish it. Before doing so, check the issues page of the project to make sure you would not duplicate someone else's work.  

#### Compiler suite preprocessor directives  

Generally, one should avoid the use of code that bind the library to a particular compiler suite. Sometimes, however, this is unavoidable. Currently, the library supports builds with the **Intel** and **GNU** compiler suites.

1.  The `INTEL_COMPILER_ENABLED` preprocessor macro must be defined for any code section that requires the Intel compiler to recognize the syntax. This preprocessor flag is automatically defined by the ParaMonte build scripts and passed to the Intel preprocessor when the Intel compiler is used to build the library.  
1.  The `GNU_COMPILER_ENABLED` preprocessor macro must be defined for any code section that requires the GNU compiler to recognize the syntax. This preprocessor flag is automatically defined by the GNU preprocessor on all platforms.  

#### Parallelism preprocessor directives  

There are currently two preprocessor directives that determine the type of parallelism to which the code section sandwiched by the preprocessor macro belongs.  

1.  The `MPI_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to MPI-parallelized versions of the kernel routines.  
1.  The `CAF_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to CAF-parallelized versions of the kernel routines (Coarray parallelism).  
1.  The `OMP_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to OMP-parallelized versions of the kernel routines (OpenMP parallelism). As of Dec 2020, there is no section of the kernel routines uses OpenMP parallelism.  

For example, if a section is MPI-parallelized, it must be fenced with the following preprocessor directive,  

```fortran  
#if defined MPI_ENABLED
call mpi_initialized( isInitialized, ierrMPI )
#endif
```  

or, if a section must be executed in either MPI or Coarray parallelisms, it must be fenced via the following preprocessor directive,  

```fortran  
#if defined CAF_ENABLED || defined MPI_ENABLED
integer :: imageID
#endif
```  

#### Library type preprocessor directives  

Any code section that must be activated for a particular library build type (`static` vs. `shared`) must be fenced with appropriate preprocessing macros.  

1.  The `DLL_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to shared library builds. Although `DLL` is reminiscent of Windows shared library files, the ParaMonte build scripts define this preprocessor macro for all shared library builds irrespective of the platform (Windows, Linux, macOS, ...).  
    > **IMPORTANT**  
    > Use this preprocessor macro in combination with **`IFORT_ENABLED`** to enable DLL symbol export by the Intel **`ifort`** compiler on macOS and Windows operating systems.  

#### Library build preprocessor directives  

Occasionally, it is necessary to define sections of code that should run only in a particular library build. Currently, the ParaMonte library build scripts support three different library builds (`debug`, `testing`, and `release`). Any code section that does not belong to the release (production) build must be appropriately fenced with appropriate preprocessor flags.  

1.  The `DEBUG_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to the `debug` library builds. Use this macro to fence sections of code that help with the debugging of the code.
1.  The `TESTING_ENABLED` preprocessor macro must be defined for any code section that exclusively belongs to the `testing` library builds. Use this macro to fence sections of code that must be activated during the testing of the library.
    > **TIP**  
    > The `TESTING_ENABLED` preprocessor macro often appears together with `DEBUG_ENABLED` and `CODECOV_ENABLED`. The testing mode deactivates all compiler optimizations, but also prevents the addition of the debugging symbols to the object files. This leads to faster compilation of the source files which is often desired for testing purposes.  

#### Library interface preprocessor directives  

To make the ParaMonte an inter-operable cross-language library, it is necessary to tailor sections of code for the particular languages of interest. The kernel connection to all languages other than Fortran is provided via the `FCI_ENABLED` preprocessor macro in the ParaMonte library.  

1.  The `CFI_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible from any programming language other than Fortran, most importantly, the C programming language.  

In addition, some code sections might have to be exclusively executed when the library is invoked from a particular programming language. In such cases, the relevant code section must be properly fenced with the appropriate language-specific preprocessor macros.  

1.  The `C_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the C programming language.  
1.  The `CPP_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the C++ programming language.  
1.  The `FORTRAN_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the Fortran programming language.  
1.  The `JULIA_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the Julia programming language.  
1.  The `MATLAB_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the MATLAB programming language.  
1.  The `MATHEMATICA_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the Mathematica programming language.  
1.  The `PYTHON_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the Python programming language.  
1.  The `R_ENABLED` preprocessor macro must be defined for any section of the code that is meant to be accessible only when the library is built for the R programming language.  

### Coding style conventions  

The following coding style are enforced within the kernel files. 
If you do not follows these rules in your contribution, please provide an explanation and justification of the alternative approach that you have taken in your contribution.  

+   [Strict naming conventions are enforced](CONTRIBUTING.md/#general-coding-style-conventions) within the entire library, including the kernel routines.

+   Strict semantic compliance with the latest Fortran standard (2018).

+   Strict source-code compliance with the latest Fortran standard. However, there are important vital exceptions to this rule,  
    +   It often preferable and sometimes essential to continue the source code lines beyond the current maximum line length limit specified by the standard, which is 132 characters. This maximum line length limit restriction is planned to be lifted in the next release of Fortran 202X standard. Nevertheless, a line length of approximately 200 characters is about where you should seriously think of breaking and continuing the code on multiple lines, if needed.  
    +   Strict parallelism compliance with the Coarray Fortran parallelism is not required (via coarrays). In other words, you can implement parallelism that is provided by external libraries, the most prominent example of which is the Message Passing Interface (MPI) Library. However, in doing so, you will have properly fence the MPI section with the proper preprocessor macro as explained in the previous section.  
    +   The minimal and efficient use of preprocessor macros and directives is allowed and encouraged to avoid code duplication and improve productivity.  

+   All Coarray variables or objects (that are communicated between processes) must begin with the `co_` prefix. Example: `co_logFunc`.  
    **Why?** This helps easily identify object that are responsible for inter-process communications in the Coarray implementation of the library.  

+   All **m**odule **v**ariables or module objects must begin with `mv_`, unless there is a good justification for not doing so. Example: `mv_state`, `comv_logFuncState`.  
    Why? module variables are global, at least within the scope of the module, and at most, everywhere in the program. Global variables are extremely susceptible of creating extremely difficult hard-to-identify bugs in the code and are hard to follow. Adding the prefix `mv_` to all global variables, resolves at least part of the problem by explicitly declaring the variable as global to the developer via its name.  

+   All **m**odule **c**onstant variables that are supposed to be defined only once throughout the entire runtime must be prefixed with `mc_` instead of `mv_`. Example: `mc_imageCount = 3`.  
    **Why?** Adding the prefix `mc_` helps other developers to realize that the object is not supposed to change value at runtime beyond the first initialization. Note that this type of variable is different from a compile-time constant, since its value must be defined at runtime, which could be different from one run to another. However, once set, the value is expected to remain constant throughout the current simulation.  

+   All constants (parameters) are upper-case separated by underscore. Example: `FILE_EXT = ".txt"`.  

+   All module names must begin with an uppercase letter and more importantly, must end with the prefix `_mod`. Example: `ParaMCMC_mod`.  
    **Why?** This helps create multiple different entity names from a single base name. For example,  
    ```fortran  
    use String_mod, only: String_type
    type(String_type) :: string
    ```  

+   Each module must be stored in a separate source file and the name of the source file must be the same as the module name, unless there is a good justification to do it otherwise. This is a strict rule to maintain the integrity of the library source file and the build scripts.  

+   All submodule names must begin with an uppercase letter and more importantly, must end with the prefix `_smod`. Example: `Routines_smod`.  

+   If a module has submodules, the name of the submodule file should follow the module name after `@`. For example the following submodule,  
    ```fortran  
    submodule (Decoration_mod) Routines_smod
    ...
    end submodule Routines_smod
    ```  
    will have to be stored in a source file with the name `Decoration_mod@Routines_smod.f90`.  

+   All type and class names must begin with an uppercase letter and more importantly, must end with the prefix `_type`. Example: `ParaMCMC_type`.  
    **Why?** This helps create multiple different entity names from a single base name. For example,  
    ```fortran  
    use String_mod, only: String_type
    type(String_type) :: string
    ```  

## Final steps  

Once you have implemented your contributions,  

+   Do not forget to test your contributions by adding new kernel tests to the unit-testing framework of the library.  
+   Also, generate code coverage report to ensure your contributions do not lower the overall code coverage of the kernel routines.  
+   Follow the [generic contribution guidelines](CONTRIBUTING.md/#all-contributors) to submit and merge your contributions with the main branch of the library on GitHub. 
**The ParaMonte library is an honor-ware and its currency is acknowledgment and citations**.  
  
If you use ParaMonte, please acknowledge it by citing the ParaMonte library's main publications as listed here:  

### The ParaMonte Python library  

+   Amir Shahmoradi, Fatemeh Bagheri, Joshua Alexander Osborne (2020). 
    **Fast fully-reproducible streamlined serial/parallel Monte Carlo/MCMC simulations and visualizations via `ParaMonte::Python` library.**. 
    Journal of Open Source Software (JOSS), to be submitted, [**PDF link**](https://www.cdslab.org/pubs/2020_Shahmoradi_III.pdf).  
    **BibTeX citation entries:**  
    ```text  
    
    @article{2020arXiv201000724S,
           author = { {Shahmoradi}, Amir and {Bagheri}, Fatemeh and {Osborne}, Joshua Alexand
            er},
            title = "{Fast fully-reproducible serial/parallel Monte Carlo and MCMC simulations and visualizations via ParaMonte::Python library}",
          journal = {arXiv e-prints},
         keywords = {Computer Science - Mathematical Software, Astrophysics - Instrumentation and Methods for Astrophysics, Quantitative Biology - Quantitative Methods, Statistics - Machine Learning},
             year = 2020,
            month = oct,
              eid = {arXiv:2010.00724},
            pages = {arXiv:2010.00724},
    archivePrefix = {arXiv},
           eprint = {2010.00724},
     primaryClass = {cs.MS},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv201000724S},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
    
    ``` 

### The ParaMonte C/C++/Fortran library  

+   Amir Shahmoradi, Fatemeh Bagheri (2020). 
    **ParaMonte: A high-performance serial/parallel Monte Carlo simulation library for C, C++, Fortran**. 
    Journal of Open Source Software (JOSS), submitted, [**PDF link**](https://www.cdslab.org/pubs/2020_Shahmoradi_II.pdf).  
    **BibTeX citation entries:**  
    ```text  
    
    @article{2020arXiv200914229S,
           author = { {Shahmoradi}, Amir and {Bagheri}, Fatemeh},
            title = "{ParaMonte: A high-performance serial/parallel Monte Carlo simulation library for C, C++, Fortran}",
          journal = {arXiv e-prints},
         keywords = {Computer Science - Mathematical Software, Astrophysics - Instrumentation and Methods for Astrophysics, Physics - Data Analysis, Statistics and Probability, Quantitative Biology - Quantitative Methods, Statistics - Machine Learning},
             year = 2020,
            month = sep,
              eid = {arXiv:2009.14229},
            pages = {arXiv:2009.14229},
    archivePrefix = {arXiv},
           eprint = {2009.14229},
     primaryClass = {cs.MS},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200914229S},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
    
    ```  

### The ParaDRAM sampler  

+   Amir Shahmoradi, Fatemeh Bagheri (2020). 
    **ParaDRAM: A Cross-Language Toolbox for Parallel High-Performance Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo Simulations**. 
    Journal of Computer Methods in Applied Mechanics and Engineering (CMAME), submitted, [**PDF link**](https://www.cdslab.org/pubs/2020_Shahmoradi_I.pdf).  
    **BibTeX citation entries:**  
    ```text  
        
        @article{2020arXiv200809589S,
                   author = { {Shahmoradi}, Amir and {Bagheri}, Fatemeh},
                    title = "{ParaDRAM: A Cross-Language Toolbox for Parallel High-Performance Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo Simulations}",
                  journal = {arXiv e-prints},
                 keywords = {Computer Science - Computational Engineering, Finance, and Science, Astrophysics - Instrumentation and Methods for Astrophysics, Physics - Data Analysis, Statistics and Probability, Statistics - Computation, Statistics - Machine Learning},
                     year = 2020,
                    month = aug,
                      eid = {arXiv:2008.09589},
                    pages = {arXiv:2008.09589},
            archivePrefix = {arXiv},
                   eprint = {2008.09589},
             primaryClass = {cs.CE},
                   adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200809589S},
                  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }
        
    ```  
  
### The ParaMonte MatDRAM MATLAB library  

+   Shashank Kumbhare, Amir Shahmoradi (2020). 
    **MatDRAM: A pure-MATLAB Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo Sampler**. 
    Journal of Computer Physics Communications (CPC), submitted, [**PDF link**](https://www.cdslab.org/pubs/2020_Kumbhare_I.pdf).  
    **BibTeX citation entries:**  
    ```text  
        
        @article{2020arXiv201004190K,
                   author = { {Kumbhare}, Shashank and {Shahmoradi}, Amir},
                    title = "{MatDRAM: A pure-MATLAB Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo Sampler}",
                  journal = {arXiv e-prints},
                 keywords = {Physics - Data Analysis, Statistics and Probability, Astrophysics - Instrumentation and Methods for Astrophysics, Computer Science - Computational Engineering, Finance, and Science, Quantitative Biology - Quantitative Methods, Statistics - Applications},
                     year = 2020,
                    month = oct,
                      eid = {arXiv:2010.04190},
                    pages = {arXiv:2010.04190},
            archivePrefix = {arXiv},
                   eprint = {2010.04190},
             primaryClass = {physics.data-an},
                   adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv201004190K},
                  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }

    ```  
  
<br>
  
For more information, visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work).  

Thank you for your interest in contributing to ParaMonte! Your help is very much appreciated. 
Below are some tips and guidelines to get started. 

## External contributors  

Here is a checklist to help you get started contributing to ParaMonte and walk you through the process,  

+   Take a look at the [issues](https://github.com/cdslaborg/paramonte/issues) page. Make sure that you're not about to duplicate someone else's work.  
+   Post a [new issue](https://github.com/cdslaborg/paramonte/issues/new/choose) discussing the changes you're proposing to implement, 
    whether bug fix(es) or enhancement(s)/feature request(s), 
    or give the rest of developers a heads up that you are going to start work on [an open issue](https://github.com/cdslaborg/paramonte/issues).  
+   [Fork the ParaMonte git project](https://help.github.com/articles/fork-a-repo/) to your private account.  
+   By contributing to the ParaMonte project you are automatically guaranteeing that,  
    1.  The contribution was created in whole or in part by you and you have the right to submit it under the open source license indicated in the file; **or**  
    1.  The contribution is based upon previous work that, to the best of your knowledge, 
    is covered under an appropriate open source license and you have the right under that license to submit that work with modifications, 
    whether created in whole or in part by you, under the same open source license as indicated in the GitHub repository; **or**  
    1.  The contribution was provided directly to me by some other person who guaranteed one of the aforementioned criteria and you have not modified it.  
    1.  You understand and agree that this project and the contribution are public and that a record of the contribution 
        (including all personal information you submit with it, including my sign-off) is maintained indefinitely 
        and may be redistributed consistent with this project or the open source license(s) involved.  
+   Follow the guidelines for [all contributors](#all-contributors) listed below.  

## All contributors  

+   [Create a branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/) and make sure to include the issue number(s) in the branch name, for example: `Provide-binary-OpenMPI-#5` or `fix-issue-#5`.  
+   Make your changes and commit them to your local repository, following these guidelines:  
    +   Each commit should be a logically atomic, self-consistent, cohesive set of changes.  
    +   The code should compile and pass all tests after each commit.  
    +   The code should be legible and any non-obvious features commented appropriately.  
    +   All unit tests should be run locally and pass (see the language-specific guidelines on how to run tests and generate code coverage report). 
    +   Tests should be added for new features and significant new code, 
        steps should be taken to ensure that the total coverage remains the same or increases.
    +   The [commit message](https://robots.thoughtbot.com/5-useful-tips-for-a-better-commit-message) should follow [these guidelines](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html):  
        +   The first line is the directive phrase, starting with a capitalized imperative verb, and is no longer than 50 characters summarizing your commit.  
        +   The next line, if necessary, is blank.  
        +   The following lines are all wrapped at 72 characters and can include additional paragraphs, bulleted lists, etc.  
        +   Use [Github keywords](https://help.github.com/articles/closing-issues-via-commit-messages/#closing-an-issue-in-a-different-repository), where appropriate, to indicate the commit resolves an open issue.  
        +   Do your best to keep a [clean and coherent history](https://www.notion.so/reviewboard/Keeping-Commit-Histories-Clean-0f717c4e802c4a0ebd852cf9337ce5d2). 
            The commands `git add -p ...`, `git commit --amend` and `git rebase --interactive <root-ref>` 
            can be helpful to rework your commits into a cleaner, clearer state.  
+   Next, [open up a pull request](https://github.com/cdslaborg/paramonte/compare) against the appropriate base branch, [`main` (formerly `master`)](https://github.com/cdslaborg/paramonte/tree/main) of [cdslaborg/paramonte](https://github.com/cdslaborg/paramonte).  
    +   In the title, please include the text `issue-<#>`, where `<#>` is replaced by 
        the issue number of the feature request or bug report corresponding to this pull request (PR).  
    +   If the PR is a work in progress, please add `WIP: ...` to the title, and rename it deleting that text once the PR is ready to be merged.  
    +   If the PR is problematic for any reason please add `DO NOT MERGE` to the title, until it is either abandoned or fixed.  
+   Please be patient and responsive to requests and comments from the PaaraMonte core team members.  
    You may be asked to amend or otherwise alter commits, or push new commits to your branch.  

## Contributors with Write Access  

The ParaMonte core developers and collaborators with push access must wait at least 24 hours before self-approving 
pull requests so that someone else has the chance to review the proposed changes and provide a formal code review. 
Due to the currently small size of the ParaMonte development team, it is unrealistic to *require* a code review under 
all circumstances. This policy ensures that there is at least an opportunity for a formal code review by another developer.  

## The ParaMonte Branches  

The ParaMonte project on GitHub uses the [Github flow](https://guides.github.com/introduction/flow/) workflow. If you are not familiar with [Github flow](https://guides.github.com/introduction/flow/), [this video](https://www.youtube.com/watch?v=EwWZbyjDs9c&feature=youtu.be&list=PLg7s6cbtAD17uAwaZwiykDci_q3te3CTY) might be a good start. The gist of it is that the `main` branch is always deploy-able and deployed. The means at anytime, a new tagged release could be shipped using the `main` branch.

### The main branch  

The `main` branch should remain pristine and stable at all times. Any changes should be applied atomically and exclusively via pull requests. It should be assumed that users are using the code on this branch,
and great care should be taken to ensure its stability. Most bug fixes and incremental improvements will get merged into the `main` branch as soon as they are deemed ready for production.

## General coding style conventions  

The following coding style are enforced within all ParaMonte source files in any programming language. If you do not follows these rules in your contribution, please provide a minimal explanation and justification of the alternative approach that you have taken in your contribution. Most importantly, strict naming conventions are enforced within the entire ParaMonte library.  

+   All names and variables and statements in the library must be self-explanatory to the highest level possible such that minimal comments would be necessary to explain the code behavior.  
    > **WARNING**  
    > Avoid short vague names that hard to decipher for variables and other objects. In particular, avoid single letter variable names, like `x`, `y`, `i`, ... .  
+   [**camelCase**](https://en.wikipedia.org/wiki/Camel_case) writing style is enforced in the entire ParaMonte library (except for constants), like, `sampleSize`, `domainUpperLimitVec`, ... .  Sometimes this convention may be against the common convention used within a particular language, for example, Python or Fortran. However, this deviation from the common practice is needed to bring homogeneity to the ParaMonte library across all programming languages. There are two advantages with using the `camelCase` naming convention:  
    +   The `camelCase` style naturally distinguishes some programming languages' intrinsic entities (for example, Python and Fortran) from the ParaMonte developer’s.  
    +   The `camelCase` style allows extremely long multi-segment variable names within the 63 character limits of many of the programming languages supported in ParaMonte.  
        > **NOTE**  
        > It is understandable that occasionally the `camelCase` style may be hard to follow and enforce in isolated locations in the library. In such cases, it is a good idea to briefly explain the reason for the deviation from the syntax rules of the library where the entity is defined for the first time in the code.  
+   Functions / subroutines / procedures in any programming language always begin with a verb. Example: `getCovarianceMatrix()` or, `getCorMatFromCovMat()`.  
+   All static functions or methods of classes begin with a lowercase verb.
+   Logical functions always begin with `is`. Example: `isDigit()`.  
+   All variables begin with a lower-case character. 
    > **TIP**  
    > An exception to this rule is the ParaMonte kernel routines where non-scalar objects used to begin with an upper-case letter. This old however, is now abandoned in favor of making the first character of all variables lower-case. This is to bring consistency with all other interfaces to the ParaMonte library.  
+   All logical variables must be English propositions that evaluate to either `true` or `false`. Example: `inputFileHasPriority`.  
+   All constants (parameters) or variable that are supposed to not change at runtime must be written upper-case, separated by underscore. Example: `FILE_EXT = ".txt"`.  
    > **NOTE**  
    > Exceptions to this rule sometimes happen in isolated scenarios. In such cases, we recommend that you provide a minimal comment next to the first appearance of the entity explaining why the deviation from the syntax rules of the library was necessary.  
+   The name of any variable that represents a vector of values, that is also anticipated to always represent a vector, is normally suffixed with `Vec`, for example: `startPointVec`, ...
+   The name of any variable that represents a matrix of values, that is also anticipated to always represent a matrix, is normally suffixed with `Mat`, for example: `proposalStartCorMat`, ...
+   The name of any variable that represents a list of varying-size values is normally suffixed with `List`, like: `variableNameList`, ...
+   Separate function arguments with a single space on the left side of the argument. For example,  
    ```python  
    def getCorMatFromCovMat(self, covMat):
    ```  
# ParaMonte kernel C/C++/Fortran release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases).  

See also,  

+   [CHANGES.md](https://github.com/cdslaborg/paramonte/blob/main/src/interface/MATLAB/CHANGES.md) for the `ParaMonte::MATLAB` release notes.  
+   [CHANGES.md](https://github.com/cdslaborg/paramonte/blob/main/src/interface/Python/CHANGES.md) for the `ParaMonte::Python` release notes.  


## **Version 1.x.x**  

### **Version 1.5.1** -- January 1, 2021  

**Minor Enhancements**  

+   This is a minor enhancement release, but is a major step toward 
    further portability of the ParaMonte library. All ParaMonte kernel 
    library dependencies are now properly handled and recognized at runtime
    without such aggressive actions as permanently redefining the environmental
    path variables, most importantly, `PATH` and `LD_LIBRARY_PATH` on Linux/macOS.

+   The ParaMonte kernel routines are now capable of handling user-input 
    file paths that contain white-space (blank) or other exotic characters.

+   The Bash build script for the ParaMonte C/C++/Fortran examples can now handle
    file paths that contain white-space (blank) or other exotic characters.

+   The shared (dynamic) library file naming convention is now changed from
    using "dynamic" in the file name to using "shared" in the file name.

+   Improved error-handling messages in the examples' build and run scripts.

+   Typo-fixes in the documentation of the library.

**Compiler Compatibility**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 8.4     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 19.1.1    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.1.1.216 Build 20200306 / Intel(R) MPI Library 2019 Update 7 for Windows  
+   **Linux**: Intel Parallel Studio Version 19.1.1.217 20200306 / Intel(R) MPI Library for Linux OS, Version 2019 Update 7 Build 20200312
+   **Linux**: GNU 10.1.0 / Open-MPI 4.0.3  
+   **Linux**: GNU 10.1.0 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  
+   **macOS**: GNU 10.2.0 / MPICH 3.3.2  

### **Version 1.5.0** -- December 17, 2020

**Major Enhancements**  

+   This version introduces numerous performance and accuracy enhancements to the ParaMonte library.
+   The entire kernel library is now fully documented and verified with over 866 tests that cover close 
    to 100% of all lines and functions in the kernel.
+   New prebuilt libraries with GNU compilers and Open-MPI on Linux are added.
+   New flags are now added to the build scripts of the library that automate the process of code coverage generation.
+   The `testing` builds are now removed from the ParaMonte release page as this build is mostly useful for development purposes.  
+   The issue of Windows file locking, that led to the occasional crashes of the 
    ParaDRAM and ParaDISE simulations in `multiChain` parallelism mode, is now resolved.

**Minor Enhancements**  

+   Minor enhancements to the ParaMonte C/C++/Fortran example build scripts `build.sh` and `build.bat`.  
+   The default build settings are now limited to `heap` memory allocation with `dynamic` library builds
    for only serial and MPI parallelization for all languages. The `stack` memory allocation results in 
    a ~10% gain in the efficiency of the code. The benefits of stack-memory builds are marginal and are 
    often problematic, in particular, for usage with non-compiled languages. However, users can still 
    build the library with `stack` memory allocation by specifying the appropriate build flags with
    the `install.sh` on Unix or `install.bat` script on Windows systems. For further information, 
    see the installation guidelines on the ParaMonte documentation website.
+   All temporary array creations in debug mode are now resolved, 
    except when Intel compilers are used, in which case, the debug warning messages are silenced.

**Compiler Compatibility**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 8.4     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 19.1.1    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.1.1.216 Build 20200306 / Intel(R) MPI Library 2019 Update 7 for Windows  
+   **Linux**: Intel Parallel Studio Version 19.1.1.217 20200306 / Intel(R) MPI Library for Linux OS, Version 2019 Update 7 Build 20200312
+   **Linux**: GNU 10.1.0 / Open-MPI 4.0.3  
+   **Linux**: GNU 10.1.0 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  
+   **macOS**: GNU 10.2.0 / MPICH 3.3.2  

### **Version 1.4.1** -- November 15, 2020

**Enhancements**  

+   The ParaMonte C/C++/Fortran example build scripts `build.sh` and `build.bat` 
    now check for the existence of both Intel and GNU compilers in the appropriate 
    order that is automatically inferred at the compilation time. Also, the dependencies 
    on the MPI compiler wrappers is now removed as the MPI libraries are not required to 
    build the ParaMonte examples, even in cases of parallel ParaMonte example builds.

**Compiler support**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 10.1.0 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  

### **Version 1.4.0** -- October 29, 2020

**Enhancements**  

+   The IO debugging info of all ParaMonte samplers have been enhanced. 
    In cases of wrong syntax or syntax-breaking input values in the simulation 
    output files, the error messages are now more informative and point directly 
    to the exact location of of error in the input file.  

+   The Integrated Autocorrelation (IAC) for sample refinement in ParaDRAM 
    sampler of ParaMonte is now set to the average of all variables' IAC values 
    instead of the maximum IAC value. This will lead to less aggressive decorrelation 
    of the final sample, which means significantly larger final sample sizes, without 
    compromising the i.i.d. property of the final refined sample. This behavior can 
    be reversed back to the original by specifying "max" or "maximum" along with 
    the requested refinement method, `SampleRefinementMethod = "batchmeans max"` 
    or `SampleRefinementMethod = "BatchMeans-max"` (case-insensitive).

**Compiler support**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 9.1 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  

### **Version 1.3.0** -- October 3, 2020

**Enhancements**  

+   A new simulation specification `overwriteRequested` has 
    been added to all ParaMonte samplers. If TRUE and the 
    ParaMonte sampler detects an existing set of old simulation 
    output files in the output path of the current simulation with 
    the same names as the output file names of the current simulation, 
    then, the ParaMonte sampler will overwrite the existing simulation files.  

**Compiler support**  
  
| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 9.1 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  

### **Version 1.2.0** -- September 22, 2020

**Enhancements**  

+   The post-processing report in the output report file 
    of ParaDRAM simulation has been significantly improved:  
    +   The parallel simulation summary now also provides the 
        predicted strong-scaling speedup behavior of the parallel 
        ParaDRAM simulations in "single chain" parallelism mode. 
        This can help make wiser decisions regarding the the number 
        of processors for similar parallel simulations in the future.  

+   The ParaDRAM restart output file in ASCII mode now contains all 
    proposal updates, including the first user-specified proposal specs.

+   All parallel simulations now avoid the unnecessary creation of 
    temporary files by all processors for System and OS operations. 
    This is particularly important for large-scale parallel simulations.
    As a side effect, this will also potentially improve the runtime 
    performances of the simulation.  

+   Major enhancements has been made to the parallel simulation 
    performance analysis reported in the post-processing section 
    of the ParaDRAM simulation output `_report.txt` files.  

+   The ParaMonte C/C++ build processes are now separate from each other.  

**Bug fixes**  

+   minor typo fixes

**Compiler support**  

| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 9.1 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.216 20200306  
+   **macOS**: GNU 10.2.0 / Open-MPI 4.0.5  

### **Version 1.1.0** -- May 27, 2020  

**New features**  

+   The MATLAB interface to the ParaMonte library is now ready to use, in addition 
    to the existing C/C++/Fortran/Python Programming language interfaces to ParaMonte.  

**Enhancements**  

+   The simulation summary report in the output report file has been improved.
+   The error handling in the build Batch-scripts on Windows is now greatly improved.

**Bug fixes**  

+   The build Batch-script for the ParaMonte examples on Windows now properly builds and runs coarray applications in parallel.
+   The fully-deterministic restart functionality is now functional also when chainFileFormat="verbose" in ParaDRAM simulations.

**Compiler support**  

| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 8.3 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.166 20191121  
+   **macOS**: GNU 9.3 / Open-MPI 4.0  

### **Version 1.0.0** -- January 1, 2020 -- Initial release  

This is the first public release of the ParaMonte library.  

**New features**  

+   ParaDRAM sampler: **Para**llel **D**elayed-**R**ejection **A**daptive Metropolis-Hastings **M**arkov Chain Monte Carlo Sampler.  
+   ParaMonte Interface to the C/C++/Fortran/Python Programming languages.  
+   ParaMonte simulation-output visualization via the ParaMonte Python interface.  

**Compiler support**  

| Compiler Suite                    | Windows (64bit) | Linux | macOS |  
|----------------------------------:|:---------------:|:-----:|:-----:|  
| GNU Compiler Collection > 7.5     | ✅              | ✅    | ✅    |  
| Intel Parallel Studio > 18.0.0    | ✅              | ✅    | ✅    |  
| Microsoft C/C++ Compiler > 16.0.0 | ✅              | ❌    | ❌    |  

**Compiler / MPI library used for this binary release**  

+   **Windows**: Intel Parallel Studio Version 19.0.4.245 Build 20190417 / Intel(R) MPI Library 2019 Update 4 for Windows  
+   **Linux**: Intel Parallel Studio Version 18.0.2 20180210 / Intel(R) MPI Library for Linux OS, Version 2018 Update 2 Build 20180125  
+   **Linux**: GNU 8.3 / MPICH 3.2  
+   **macOS**: Intel Parallel Studio Version 19.1.0.166 20191121  
+   **macOS**: GNU 9.3 / Open-MPI 4.0  
## Note

For more examples, in particular, Jupyter Notebooks and MATLAB live scripts, visit the ParaMonte's examples repository at [https://github.com/cdslaborg/paramonte.svg?branch=main](https://github.com/cdslaborg/paramonte.svg?branch=main)  


  
<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank"><img src="https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png" alt="ParaMonte: Plain Powerful Parallel Monte Carlo Library" /></a>
<br><br>
<a href="https://github.com/cdslaborg/paramonte#license" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://travis-ci.com/cdslaborg/paramonte" target="_blank"><img src="https://travis-ci.com/cdslaborg/paramonte.svg?branch=main&style=flat-square" alt="Build Status" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/serial/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-serial%20kernel%20:%2098.2%25-brightgreen?style=flat-square" alt="kernel code coverage - serial" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/mpi/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-MPI%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - MPI" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/caf/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-Coarray%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - Coarray" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=PyPI%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/badge/dynamic/json?style=flat-square&labelColor=grey&color=brightgreen&maxAge=86400&label=PyPI%20downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fparamonte" alt="PyPI - Downloads Total" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
<br><br>
<a href="#paramonte-plain-powerful-parallel-monte-carlo-library">Overview</a> | 
<a href="#installation">Installation</a> | 
<a href="#dependencies">Dependencies</a> | 
<a href="#parallelism">Parallelism</a> | 
<a href="#example-usage-instructions">Examples</a> |
<a href="#citing-paramonte">Acknowledgments</a> | 
<a href="#license">License</a> | 
<a href="#authors-and-contributors">Authors</a>  
</div>
  
  
ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions 
of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, 
and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), 
**user-friendliness** (of the library), **accessibility** (from multiple programming environments), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  
  
  
ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  

-   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray 
    communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations up to 16 digits of precision as requested by the user.  

-   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in 
    external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

  
Installation  
============  

The pre-built ready-to-use libraries are available on [the release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases). 
Each prebuilt ParaMonte library automatically ships with a full-fledged set of example codes and build scripts.  

Alternatively, you can build the library from the source in the [GitHub repository of the project](https://github.com/cdslaborg/paramonte). 
The ParaMonte library installation/build process is fully automated for all of the supported programming languages. 
Currently, the following compiler suites are supported **for builds from source**:  
  
| Compiler Suite                    | Linux | macOS | Windows (64bit) |  
|----------------------------------:|:-----:|:-----:|:---------------:|  
| GNU Compiler Collection > 8.4     |&check;|&check;| &cross;         |  
| Intel Parallel Studio > 19.1.1    |&check;|&check;| &check;         |  

For more information and quick-start in the programming language of your choice, visit the [ParaMonte library homepage](https://www.cdslab.org/paramonte).  

  
Dependencies  
============  

Beyond an optional MPI runtime library for parallel simulations, the ParaMonte kernel has **zero dependency** on external third-party libraries or packages.  

  
Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. 
To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. 
In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). 
These automatic checks and installations happen when you download and install or use the library on your system, for the first time. 
If the automatic installation is unsuccessful, you can also install the libraries manually on your system:  

+   On **Windows** and **Linux** operating systems, we highly recommend downloading and installing the 
    [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), 
    which is available to the public free of charge, also available in the latest release of the 
    ParaMonte library on the [GitHub release page](https://github.com/cdslaborg/paramonte/releases) 
    (For Windows, look for the executable file that ends with `.exe`. For Linux, look for the file 
    that ends with `.tgz`, like `l_mpi-rt_2018.2.199.tgz`).
+   On **macOS**, the Intel MPI library is not available. Therefore, we recommend installing either 
    [Open-MPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/) MPI runtime libraries 
    depending the prebuilt version of the ParaMonte library that you have downloaded or 
    the configuration with which you intend to build the library.  

For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

  
Example usage instructions  
==========================  

+   For complete clear organized up-to-date instructions on the build process and the installation of the ParaMonte library, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

## Quick start  

+   Go to the [release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases).  
+   Decide on the parallelism paradigm that you want to use: serial / MPI.  
+   Decide on the Operating System (OS) on which you want to run the ParaMonte simulations: Windows / macOS / Linux.  
+   Learn about the naming convention used for the ParaMonte prebuilt libraries
    [here](https://www.cdslab.org/paramonte/notes/installation/readme/#naming-convention-used-for-paramonte-library-builds).  
+   Download the prebuilt ParaMonte library of your choice based on the decisions you have made in the above. 
    If you are not sure which prebuilt library is suitable for your needs, use the prebuilt library recommended 
    [here for Windows](https://www.cdslab.org/paramonte/notes/installation/windows/#using-the-prebuilt-paramonte-library), 
    or [here for Linux](https://www.cdslab.org/paramonte/notes/installation/linux/#using-the-prebuilt-paramonte-library), 
    or [here for macOS](https://www.cdslab.org/paramonte/notes/installation/macos/#using-the-prebuilt-paramonte-library).  
+   Each prebuilt library ships with a full-fledged set of example codes and build scripts. Uncompress the prebuilt library:  
    +   On **Windows**: Simply double-click on the zip-file and select **extract files** from the Windows Explorer menu.  
    +   On **macOS/Linux**: Open a Bash terminal and navigate to the folder containing the compressed library. 
        Use the following command to untar the compressed file,  
        ```  
        ls libparamonte*.tar.gz* | xargs -i tar xvzf {}
        ```  
        to extract all libparamonte tar files in the current directory.  

### Building and running ParaMonte simulations on Windows  

+   **Note**: Theoretically, you can use any C/C++ compiler on Windows to build and link your applications against the ParaMonte library. 
    However, the ParaMonte library example build scripts, as described below, currently only recognize the Microsoft and Intel C/C++ compilers.  

+   **Note**: Theoretically, you can use any Fortran compiler on Windows to build and link your applications against the ParaMonte library.  
    A few options currently exist regarding the choice of compilers and environment:  
    +   Use the Intel Visual C/C++ compiler or the Microsoft Visual C++ Compiler along with the Windows **Batch build scripts** of the ParaMonte library examples.  
    +   Use the **Bash build scripts** that are also supplied with each ParaMonte example on Windows to build and run simulations via 
        the GNU compilers on Windows available in MinGW or Cygwin Linux environments installed on a Windows system.  
    +   Use the GNU compiler installed on a [Microsoft Windows Subsystem for Linux (WSL)](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). 
        In this case, you will have to download the prebuilt ParaMonte library for the Linux environment as opposed to the Windows OS. 

+   **Install the Microsoft Visual Studio (>2017)**: You will need to have a recent Microsoft Visual Studio (MSVS) installed on your system. 
    The community edition of this software is available free of charge. When installing MSVS, 
    make sure to install all the C++ components and compiler of the Visual Studio.  

+   **Install the Intel Parallel Studio**: If you are a student/teacher/faculty/open-source-developer, you can also download and install, free of charge, 
    the most recent **Intel Parallel Studio** on your system which, by default, includes the Intel MPI library. You can follow the instructions given 
    [on this page](https://www.cdslab.org/recipes/programming/intel-parallel-studio-installation-windows/intel-parallel-studio-installation-windows) 
    to install the Intel Parallel Studio on your system.  

+   **Open the right command-line interface to build/run the ParaMonte example**: 
    If the ParaMonte library that you intend to use is built for 64-bit architecture, 
    then make sure you open a 64-bit instance of the command-line interface in either of the two cases below:  
    +   If you have installed Intel Parallel Studio, open an instance of the **command-line interface** 
        that comes with Intel Parallel Studio from the list of programs in the Windows start menu. 
        This is simply a Windows command prompt that has all the necessary Intel compiler variables and paths predefined in it.  
    +   Otherwise, if you do not have Intel Parallel Studio, open an instance of the **command-line interface** that comes with 
        the Microsoft Visual Studio from the list of programs in the Windows start menu. This is simply a Windows command prompt 
        that has all the necessary Microsoft compiler variables and paths predefined in it.  

+   **Build and run the ParaMonte example**:  
    +   To Build the example via the Intel Parallel Studio command-line interface,  
        ```  
        build.bat  
        ```  
    +   To build the example via the Microsoft Visual Studio C/C++ compiler,  
        ```  
        build.bat msvc  
        ```  
        where the passed argument `msvc` implies the use of Microsoft Visual C++ compiler for building the application.  
    The build script will automatically detect whether a parallel simulation has been built. By default, the name of 
    the output executable is `main.exe`. *Note that the build script will both build and run the executable*.    

+   **Run the ParaMonte example executable**:  
    +   For serial simulations, simply type the name of the output executable,  
        ```
        main.exe
        ```
    +   For parallel simulations, invoke the MPI launcher `mpiexec`,  
        ```
        mpiexec -n NUM_PROCESSES main.exe
        ```
        where `NUM_PROCESSES` represents the number of processes on which the simulation will run. 
        If you are using the Intel MPI library to run your ParaMonte application in parallel, we also recommend using the `-localonly` flag. 
        See [this page](https://www.cdslab.org/paramonte/notes/run/#running-the-manually-generated-executable-on-multiple-processors-on-windows) 
        for usage and utilities of this Intel MPI launcher flag.  

### Building and running ParaMonte simulations on macOS / Linux  

+   **Note**: Theoretically, you can use any C/C++ compiler on macOS/Linux to build and link your applications against the ParaMonte library. 
    However, the ParaMonte library example build scripts, as described below, currently only recognize the Intel and GNU C/C++ compilers.  

+   If you intend to run **serial** ParaMonte simulations, install either,  
    +   **the Intel C/C++ compiler (icc/icpc >2018)**, or,  
    +   **the GNU C/C++ compiler (gcc/g++ >7.0.0)**,  
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  

+   If you intend to run **MPI parallel** ParaMonte simulations, install either,  
    +   **the Intel Parallel Studio (>2018)** on Linux, or,  
    +   **the GNU Compiler Collection (>7.0.0) and the MPICH (>3.2) library** on Linux,  
    +   **the GNU Compiler Collection (>7.0.0) and the OpenMPI (>4.0) library** on macOS,  
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  
    Note that on **macOS**, only the latter option (the GNU compilers) is 
    available since the Intel MPI library does not support the macOS platform.  

+   Open a Bash terminal, change directory to the ParaMonte prebuilt library's directory, then build the executable via,  
    ```  
    build.sh  
    ```  
    The build script will automatically detect whether a parallel simulation has to be built. 
    By default, the name of the output executable is `main.exe`. The script will also generate a new Bash script named `run.sh`. 
    To run the generated example executable, type,  
    ```  
    ./run.sh
    ```  
    The script will automatically detect whether the simulation has to be run in parallel or serial. 
    If the simulation is parallel, you can also pass the number of cores on which you want to run the example via,  
    ```  
    ./run.sh --nproc NUM_PROCESSOR
    ```  
    or,  
    ```  
    ./run.sh -n NUM_PROCESSOR
    ```  
    where you will have to replace `NUM_PROCESSOR` with your desired number of processes.  

  
Citing ParaMonte  
================  

The ParaMonte library is an honor-ware, the currency of which is acknowledgment and citations.  
  
If you use ParaMonte or any ideas from the software, please acknowledge it by citing the ParaMonte library's 
main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) 
to access the PDF version of these files free of charge.  

  
License  
=======  

[MIT License](https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md)  

**What does this license mean?**  

Essentially, all we are asking from the users or developers is to  

>   explicitly acknowledge the use of this library or any concepts or parts of it in their education, research, or software (free or commercial).  

This is a free software, so help us keep it freely available to the public by redistributing the library and contributing to it. 
If you have questions or concerns about the license, do not hesitate to contact us (shahmoradi@utexas.edu).  

  
Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, 
        data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, 
        high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, 
        biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, 
        protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+   [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)  
    +   physicist / cosmologist by training,  
    +   currently a UTA Physics member,  
    +   deep philosophical thinker,  
    +   contact: [Fatemeh.Bagheri@uta.edu](mailto:"Fatemeh.Bagheri@uta.edu")  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

  
<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank"><img src="https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png" alt="ParaMonte: Plain Powerful Parallel Monte Carlo Library" /></a>
<br><br>
<a href="https://github.com/cdslaborg/paramonte#license" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://travis-ci.com/cdslaborg/paramonte" target="_blank"><img src="https://travis-ci.com/cdslaborg/paramonte.svg?branch=main&style=flat-square" alt="Build Status" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/serial/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-serial%20kernel%20:%2098.2%25-brightgreen?style=flat-square" alt="kernel code coverage - serial" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/mpi/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-MPI%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - MPI" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/caf/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-Coarray%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - Coarray" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=PyPI%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/badge/dynamic/json?style=flat-square&labelColor=grey&color=brightgreen&maxAge=86400&label=PyPI%20downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fparamonte" alt="PyPI - Downloads Total" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
<br><br>
<a href="#paramonte-plain-powerful-parallel-monte-carlo-library">Overview</a> | 
<a href="#installation">Installation</a> | 
<a href="#dependencies">Dependencies</a> | 
<a href="#parallelism">Parallelism</a> | 
<a href="#example-usage-instructions">Examples</a> |
<a href="#citing-paramonte">Acknowledgments</a> | 
<a href="#license">License</a> | 
<a href="#authors-and-contributors">Authors</a>  
</div>
  
  
ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions 
of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, 
and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), 
**user-friendliness** (of the library), **accessibility** (from multiple programming environments), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  
  
  
ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  

-   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray 
    communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations up to 16 digits of precision as requested by the user.  

-   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in 
    external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

  
Installation  
============  

The pre-built ready-to-use libraries are available on [the release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases). 
Each prebuilt ParaMonte library automatically ships with a full-fledged set of example codes and build scripts.  

Alternatively, you can build the library from the source in the [GitHub repository of the project](https://github.com/cdslaborg/paramonte). 
The ParaMonte library installation/build process is fully automated for all of the supported programming languages. 
Currently, the following compiler suites are supported **for builds from source**:  
  
| Compiler Suite                    | Linux | macOS | Windows (64bit) |  
|----------------------------------:|:-----:|:-----:|:---------------:|  
| GNU Compiler Collection > 8.4     |&check;|&check;| &cross;         |  
| Intel Parallel Studio > 19.1.1    |&check;|&check;| &check;         |  

For more information and quick-start in the programming language of your choice, visit the [ParaMonte library homepage](https://www.cdslab.org/paramonte).  

  
Dependencies  
============  

Beyond an optional MPI runtime library for parallel simulations, the ParaMonte kernel has **zero dependency** on external third-party libraries or packages.  

  
Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. 
To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. 
In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). 
These automatic checks and installations happen when you download and install or use the library on your system, for the first time. 
If the automatic installation is unsuccessful, you can also install the libraries manually on your system:  

+   On **Windows** and **Linux** operating systems, we highly recommend downloading and installing the 
    [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), 
    which is available to the public free of charge, also available in the latest release of the 
    ParaMonte library on the [GitHub release page](https://github.com/cdslaborg/paramonte/releases) 
    (For Windows, look for the executable file that ends with `.exe`. For Linux, look for the file 
    that ends with `.tgz`, like `l_mpi-rt_2018.2.199.tgz`).
+   On **macOS**, the Intel MPI library is not available. Therefore, we recommend installing either 
    [Open-MPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/) MPI runtime libraries 
    depending the prebuilt version of the ParaMonte library that you have downloaded or 
    the configuration with which you intend to build the library.  

For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

  
Example usage instructions  
==========================  

+   For complete clear organized up-to-date instructions on the build process and the installation of the ParaMonte library, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

## Quick start  

+   Go to the [release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases).  
+   Decide on the parallelism paradigm that you want to use: serial / MPI.  
+   Decide on the Operating System (OS) on which you want to run the ParaMonte simulations: Windows / macOS / Linux.  
+   Learn about the naming convention used for the ParaMonte prebuilt libraries
    [here](https://www.cdslab.org/paramonte/notes/installation/readme/#naming-convention-used-for-paramonte-library-builds).  
+   Download the prebuilt ParaMonte library of your choice based on the decisions you have made in the above. 
    If you are not sure which prebuilt library is suitable for your needs, use the prebuilt library recommended 
    [here for Windows](https://www.cdslab.org/paramonte/notes/installation/windows/#using-the-prebuilt-paramonte-library), 
    or [here for Linux](https://www.cdslab.org/paramonte/notes/installation/linux/#using-the-prebuilt-paramonte-library), 
    or [here for macOS](https://www.cdslab.org/paramonte/notes/installation/macos/#using-the-prebuilt-paramonte-library).  
+   Each prebuilt library ships with a full-fledged set of example codes and build scripts. Uncompress the prebuilt library:  
    +   On **Windows**: Simply double-click on the zip-file and select **extract files** from the Windows Explorer menu.  
    +   On **macOS/Linux**: Open a Bash terminal and navigate to the folder containing the compressed library. 
        Use the following command to untar the compressed file,  
        ```  
        ls libparamonte*.tar.gz* | xargs -i tar xvzf {}
        ```  
        to extract all libparamonte tar files in the current directory.  

### Building and running ParaMonte simulations on Windows  

+   **Note**: Theoretically, you can use any C/C++ compiler on Windows to build and link your applications against the ParaMonte library. 
    However, the ParaMonte library example build scripts, as described below, currently only recognize the Microsoft and Intel C/C++ compilers.  

+   **Note**: Theoretically, you can use any Fortran compiler on Windows to build and link your applications against the ParaMonte library.  
    A few options currently exist regarding the choice of compilers and environment:  
    +   Use the Intel Visual C/C++ compiler or the Microsoft Visual C++ Compiler along with the Windows **Batch build scripts** of the ParaMonte library examples.  
    +   Use the **Bash build scripts** that are also supplied with each ParaMonte example on Windows to build and run simulations via 
        the GNU compilers on Windows available in MinGW or Cygwin Linux environments installed on a Windows system.  
    +   Use the GNU compiler installed on a [Microsoft Windows Subsystem for Linux (WSL)](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). 
        In this case, you will have to download the prebuilt ParaMonte library for the Linux environment as opposed to the Windows OS. 

+   **Install the Microsoft Visual Studio (>2017)**: You will need to have a recent Microsoft Visual Studio (MSVS) installed on your system. 
    The community edition of this software is available free of charge. When installing MSVS, 
    make sure to install all the C++ components and compiler of the Visual Studio.  

+   **Install the Intel Parallel Studio**: If you are a student/teacher/faculty/open-source-developer, you can also download and install, free of charge, 
    the most recent **Intel Parallel Studio** on your system which, by default, includes the Intel MPI library. You can follow the instructions given 
    [on this page](https://www.cdslab.org/recipes/programming/intel-parallel-studio-installation-windows/intel-parallel-studio-installation-windows) 
    to install the Intel Parallel Studio on your system.  

+   **Open the right command-line interface to build/run the ParaMonte example**: 
    If the ParaMonte library that you intend to use is built for 64-bit architecture, 
    then make sure you open a 64-bit instance of the command-line interface in either of the two cases below:  
    +   If you have installed Intel Parallel Studio, open an instance of the **command-line interface** 
        that comes with Intel Parallel Studio from the list of programs in the Windows start menu. 
        This is simply a Windows command prompt that has all the necessary Intel compiler variables and paths predefined in it.  
    +   Otherwise, if you do not have Intel Parallel Studio, open an instance of the **command-line interface** that comes with 
        the Microsoft Visual Studio from the list of programs in the Windows start menu. This is simply a Windows command prompt 
        that has all the necessary Microsoft compiler variables and paths predefined in it.  

+   **Build and run the ParaMonte example**:  
    +   To Build the example via the Intel Parallel Studio command-line interface,  
        ```  
        build.bat  
        ```  
    +   To build the example via the Microsoft Visual Studio C/C++ compiler,  
        ```  
        build.bat msvc  
        ```  
        where the passed argument `msvc` implies the use of Microsoft Visual C++ compiler for building the application.  
    The build script will automatically detect whether a parallel simulation has been built. By default, the name of 
    the output executable is `main.exe`. *Note that the build script will both build and run the executable*.    

+   **Run the ParaMonte example executable**:  
    +   For serial simulations, simply type the name of the output executable,  
        ```
        main.exe
        ```
    +   For parallel simulations, invoke the MPI launcher `mpiexec`,  
        ```
        mpiexec -n NUM_PROCESSES main.exe
        ```
        where `NUM_PROCESSES` represents the number of processes on which the simulation will run. 
        If you are using the Intel MPI library to run your ParaMonte application in parallel, we also recommend using the `-localonly` flag. 
        See [this page](https://www.cdslab.org/paramonte/notes/run/#running-the-manually-generated-executable-on-multiple-processors-on-windows) 
        for usage and utilities of this Intel MPI launcher flag.  

### Building and running ParaMonte simulations on macOS / Linux  

+   **Note**: Theoretically, you can use any C/C++ compiler on macOS/Linux to build and link your applications against the ParaMonte library. 
    However, the ParaMonte library example build scripts, as described below, currently only recognize the Intel and GNU C/C++ compilers.  

+   If you intend to run **serial** ParaMonte simulations, install either,  
    +   **the Intel C/C++ compiler (icc/icpc >2018)**, or,  
    +   **the GNU C/C++ compiler (gcc/g++ >7.0.0)**,  
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  

+   If you intend to run **MPI parallel** ParaMonte simulations, install either,  
    +   **the Intel Parallel Studio (>2018)** on Linux, or,  
    +   **the GNU Compiler Collection (>7.0.0) and the MPICH (>3.2) library** on Linux,  
    +   **the GNU Compiler Collection (>7.0.0) and the OpenMPI (>4.0) library** on macOS,  
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  
    Note that on **macOS**, only the latter option (the GNU compilers) is 
    available since the Intel MPI library does not support the macOS platform.  

+   Open a Bash terminal, change directory to the ParaMonte prebuilt library's directory, then build the executable via,  
    ```  
    build.sh  
    ```  
    The build script will automatically detect whether a parallel simulation has to be built. 
    By default, the name of the output executable is `main.exe`. The script will also generate a new Bash script named `run.sh`. 
    To run the generated example executable, type,  
    ```  
    ./run.sh
    ```  
    The script will automatically detect whether the simulation has to be run in parallel or serial. 
    If the simulation is parallel, you can also pass the number of cores on which you want to run the example via,  
    ```  
    ./run.sh --nproc NUM_PROCESSOR
    ```  
    or,  
    ```  
    ./run.sh -n NUM_PROCESSOR
    ```  
    where you will have to replace `NUM_PROCESSOR` with your desired number of processes.  

  
Citing ParaMonte  
================  

The ParaMonte library is an honor-ware, the currency of which is acknowledgment and citations.  
  
If you use ParaMonte or any ideas from the software, please acknowledge it by citing the ParaMonte library's 
main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) 
to access the PDF version of these files free of charge.  

  
License  
=======  

[MIT License](https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md)  

**What does this license mean?**  

Essentially, all we are asking from the users or developers is to  

>   explicitly acknowledge the use of this library or any concepts or parts of it in their education, research, or software (free or commercial).  

This is a free software, so help us keep it freely available to the public by redistributing the library and contributing to it. 
If you have questions or concerns about the license, do not hesitate to contact us (shahmoradi@utexas.edu).  

  
Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, 
        data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, 
        high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, 
        biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, 
        protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

  
<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank"><img src="https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png" alt="ParaMonte: Plain Powerful Parallel Monte Carlo Library" /></a>
<br><br>
<a href="https://github.com/cdslaborg/paramonte#license" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://travis-ci.com/cdslaborg/paramonte" target="_blank"><img src="https://travis-ci.com/cdslaborg/paramonte.svg?branch=main&style=flat-square" alt="Build Status" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/serial/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-serial%20kernel%20:%2098.2%25-brightgreen?style=flat-square" alt="kernel code coverage - serial" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/mpi/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-MPI%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - MPI" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/caf/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-Coarray%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - Coarray" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=PyPI%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/badge/dynamic/json?style=flat-square&labelColor=grey&color=brightgreen&maxAge=86400&label=PyPI%20downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fparamonte" alt="PyPI - Downloads Total" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
<br><br>
<a href="#paramonte-plain-powerful-parallel-monte-carlo-library">Overview</a> | 
<a href="#installation">Installation</a> | 
<a href="#dependencies">Dependencies</a> | 
<a href="#parallelism">Parallelism</a> | 
<a href="#example-usage-instructions">Examples</a> |
<a href="#citing-paramonte">Acknowledgments</a> | 
<a href="#license">License</a> | 
<a href="#authors-and-contributors">Authors</a>  
</div>
  
  
ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions 
of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, 
and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), 
**user-friendliness** (of the library), **accessibility** (from multiple programming environments), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  
  
  
ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  

-   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray 
    communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations up to 16 digits of precision as requested by the user.  

-   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in 
    external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

  
Quick start  
===========  

For a quick start with some Jupyter Notebook examples, visit [this ParaMonte documentation page](https://www.cdslab.org/paramonte/notes/examples/python/jupyter/). 
The corresponding example source files (the `*.ipynb` files) can be downloaded from the [paramontex GitHub repository](https://github.com/cdslaborg/paramontex/tree/main/Python/Jupyter), 
which is a repository dedicated to the ParaMonte library examples.  

The following example code samples a 4-dimensional MultiVariate Normal (MNV) distribution via the ParaDRAM sampler in serial mode,  

```python  
import numpy as np
import paramonte as pm
def getLogFunc(point): return -0.5 * np.dot(point, point)
pmpd = pm.ParaDRAM()
pmpd.runSampler ( ndim = 4 # assume 4-dimensional objective function
                , getLogFunc = getLogFunc   # the objective function
                )
```  

To learn about the post-processing and visualization tools of the `ParaMonte::Python` library, visit [this this documentation page](https://www.cdslab.org/paramonte/notes/examples/python/jupyter/).  

  
Installation  
============  

The latest release of ParaMonte can be installed from PyPI using `pip`:  

    pip3 install --user --upgrade paramonte  

or,  

    pip install --user --upgrade paramonte  

Alternatively, you can build the library from the source in the GitHub repository of the project ([https://github.com/cdslaborg/paramonte](https://github.com/cdslaborg/paramonte)). 
For instructions, please visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

  
Dependencies  
============  

The Python interface of ParaMonte depends on a very few third-party libraries. 
These include `numpy`, `scipy`, `pandas`, `matplotlib`, and `seaborn`. 
The last two (plotting) libraries are only used for the post-processing of simulation 
results and are therefore not needed if you do not plan to use the post-processing 
features of the ParaMonte library. If you have a recent version of Anaconda Python 
distribution installed on your system, then all of the dependencies already 
exist and are automatically installed on your system.  

  
Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. 
To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. 
In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). 
These automatic checks and installations happen when you download and install or use the library on your system, for the first time. 
If the automatic installation is unsuccessful, you can also install the libraries manually on your system:  

+   On **Windows** and **Linux** operating systems, we highly recommend downloading and installing the 
    [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), 
    which is available to the public free of charge, also available in the latest release of the 
    ParaMonte library on the [GitHub release page](https://github.com/cdslaborg/paramonte/releases) 
    (For Windows, look for the executable file that ends with `.exe`. For Linux, look for the file 
    that ends with `.tgz`, like `l_mpi-rt_2018.2.199.tgz`).
+   On **macOS**, the Intel MPI library is not available. Therefore, we recommend installing either 
    [Open-MPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/) MPI runtime libraries 
    depending the prebuilt version of the ParaMonte library that you have downloaded or 
    the configuration with which you intend to build the library.  

For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

  
Example usage instructions  
==========================  

+   **Install a Python 3 distribution**, preferably, the Anaconda distribution of Python. 
    The Anaconda distribution of Python automatically ships with all of the 
    ParaMonte Python package dependencies when installed on your system.  

+   **Optionally install a compatible MPI library** (or let the ParaMonte library take care of the installation 
    when you import the package into your Python session for the first time). For parallel simulations (via MPI), 
    you will need an MPI library already installed on your system. If you choose to install the library by yourself, 
    we recommend the Intel MPI library which is available for free from the Intel website. On macOS, the OpenMPI 
    library can be used in place of the Intel MPI library which currently does not support macOS.  

+   **Running the ParaMonte simulations**  
    +   Open an Anaconda command-line interface or `jupyter` notebook.  
    +   Suppose your mathematical objective function is a multivariate Normal distribution as implemented in this   
        [logfunc.py](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/mvn/Python/logfunc.py) file.  
    +   For **serial** simulations, download this example generic serial 
        [main.py](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main.py) 
        Python main file and save it in the same folder containing the `logfunc.py` file that you downloaded in the above. 
        Then, simply type the name of the Python main script, `python main.py` on the Bash terminal or the Anaconda command line.  
    +   For **parallel** simulations, download this example generic parallel 
        [main_mpi.py](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main_mpi.py) Python main file and save it in the same folder containing the `logfunc.py` file that you downloaded in the above. Then, simply invoke the MPI launcher followed by the name of the Python main script on the Bash terminal, similar to the following,  
        +   on Windows (within the Anaconda command line or a terminal that recognizes both `mpiexec` and `python` software),  
            ```  
            mpiexec -localonly -n 3 python main_mpi.py
            ```  
            where the `-localonly` flag is needed only if you are using the Intel MPI runtime libraries (which is the default MPI library used to build the ParaMonte libraries on Windows).  
        +   on macOS or Linux (within a Bash terminal),  
            ```  
            mpiexec -n 3 python main_mpi.py
            ```  
        Here, the parallel simulations are performed on 3 processes. Change the number 3 to any number of processes you wish to use, 
        but do not go beyond the maximum number of physical processes available on your system, otherwise, it will only degrade 
        the performance of your parallel simulations. For example, if you are running the parallel simulation on a personal 
        quad-cores laptop, set the number of processes to either 3 or 4 at most.  
    +   Enjoy the unification of simplicity, efficiency, and parallelism in Monte Carlo simulations!  
    +   The ParaMonte library samplers are extremely versatile with many adjustable input parameters. 
        To learn about the many advanced features of the ParaMonte routines, visit: https://www.cdslab.org/paramonte  

  
Citing ParaMonte  
================  

The ParaMonte library is an honor-ware, the currency of which is acknowledgment and citations.  
  
If you use ParaMonte or any ideas from the software, please acknowledge it by citing the ParaMonte library's 
main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) 
to access the PDF version of these files free of charge.  

  
License  
=======  

[MIT License](https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md)  

**What does this license mean?**  

Essentially, all we are asking from the users or developers is to  

>   explicitly acknowledge the use of this library or any concepts or parts of it in their education, research, or software (free or commercial).  

This is a free software, so help us keep it freely available to the public by redistributing the library and contributing to it. 
If you have questions or concerns about the license, do not hesitate to contact us (shahmoradi@utexas.edu).  

  
Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, 
        data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, 
        high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, 
        biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, 
        protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+   [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)  
    +   physicist / cosmologist by training,  
    +   currently a UTA Physics member,  
    +   deep philosophical thinker,  
    +   contact: [Fatemeh.Bagheri@uta.edu](mailto:"Fatemeh.Bagheri@uta.edu")  

+   [Joshua Osborne](https://www.cdslab.org/people/#joshua-alexander-osborne)  
    +   physicist / Computational Data Scientist by training,  
    +   currently a UTA Physics member,  
    +   contact: [joshuaalexanderosborne@gmail.com](mailto:"joshuaalexanderosborne@gmail.com")  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
# ParaMonte Python release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases) or [the ParaMonte page on the Python Package Index](https://pypi.org/project/paramonte/).  

## **Version 2.x.x**  

### Version  2.5.2 -- January 8, 2021  

+   Minor enhancements to the `checkForUpdate()` function of paramonte module. 
    This function now relies on a more robust method of latest-version-checking.  

### Version  2.5.1 -- January 3, 2021  

+   Minor enhancements and bug fixes.

### Version  2.5.0 -- January 1, 2021  

**Major enhancements**  

+   This release is a a major step toward further portability 
    of the kernel routines of the `ParaMonte::Python` library. The kernel 
    library dependencies are now properly handled and recognized at runtime
    without such aggressive actions as permanently redefining the environmental
    path variables, most importantly, `PATH` and `LD_LIBRARY_PATH` on Linux/macOS.

+   The `ParaMonte::Python` library is now capable of recognizing the existing MPI
    libraries such as MPICH and OpenMPI on user's system and avoid further 
    installation of a new MPI library if it is deemed unnecessary.  

+   The ParaMonte kernel routines are now capable of handling user-input 
    file paths that contain white-space (blank) or other exotic characters.  

+   Enhancements to the `build()` function of the paramonte module.  

**Minor enhancements**  

+   Typo-fixes in the documentation of the library.

+   The `ParaMonte::Python` library packages for different Operating systems and 
    processor architecture are now separate from each other. This change was made
    to lower the overall size of `ParaMonte::Python` by only keeping the relevant 
    files in each packaging of the library. The current release contains three 
    separate packages for `ParaMonte::Python`,
    +   `libparamonte_python_windows_x64`,  
    +   `libparamonte_python_darwin_x64`,  
    +   `libparamonte_python_linux_x64`.  

### Version  2.4.0 -- December 17, 2020  

**Major enhancements**  

+   This update presents several major performance, accuracy, 
    and verification enhancements to the ParaMonte kernel routines,
    in particular, to the ParaDRAM sampler.

+   An extensive set of over 866 tests have been added 
    that test all aspects of the ParaMonte kernel library.

+   The issue of Windows file locking, that led to the occasional crashes of the 
    ParaDRAM and ParaDISE simulations in `multiChain` parallelism mode, is now resolved.

+   All overflow / underflow exceptions are now properly handled.  

+   The `ParaDRAM` class in `paramonte` is now also available 
    as `Paradram` and `paradram`, although the original label 
    will remain the default preferred method of ParaDRAM 
    object instantiation.  

**Minor enhancements**  

+   The interactive mode for displaying plots is now automatically on. 
    The plots are not automatically displayed in `ipython` sessions.  

### Version  2.3.1 -- November 15, 2020  

**Minor enhancements**  

+   Minor enhancements to the Kernel library 
    build scripts and dependencies management.

### Version  2.3.0 -- October 29, 2020  

**Enhancements**  

+   The IO debugging info of all ParaMonte samplers have been enhanced. 
    In cases of wrong syntax or syntax-breaking input values in the simulation 
    output files, the error messages are now more informative and point directly 
    to the exact location of of error in the input file.  

+   The Integrated Autocorrelation (IAC) for sample refinement in ParaDRAM 
    sampler of ParaMonte is now set to the average of all variables' IAC values 
    instead of the maximum IAC value. This will lead to less aggressive decorrelation 
    of the final sample, which means significantly larger final sample sizes, without 
    compromising the i.i.d. property of the final refined sample. This behavior can 
    be reversed back to the original by specifying "max" or "maximum" along with 
    the requested refinement method, `SampleRefinementMethod = "batchmeans-max"` 
    or `SampleRefinementMethod = "BatchMeans-maximum"` (case-insensitive).

### Version  2.2.5 -- October 20, 2020  

**Minor enhancements**  

+   Enhancements to the messages for incompatible architecture.  

### Version  2.2.4 -- October 18, 2020  

**Minor enhancements**  

+   The bug preventing the setting of 
    SpecDRAM specifications is now fixed.  

### Version  2.2.3 -- October 15, 2020  

**Minor enhancements**  

+   Further enhancements and corrections to the 
    `checkForUpdate()` function of `paramonte` module.  

### Version  2.2.2 -- October 14, 2020  

**Minor enhancements**  

+   The `checkForUpdate()` function of `paramonte` module now 
    works fine when a newer version of the software is available.  

### Version  2.2.1 -- October 11, 2020  

**Minor enhancements**  

+   Documentation enhancements and typo fixes.

### Version  2.2.0 -- October 11, 2020  

**Major enhancements**  

+   The `ParaMonte::Python` output file parsers method are now 
    capable of parsing simulation output contents directly from 
    the web. All that is needed, is to provide the link to the 
    web file as the input file name to the simulation output 
    file parser methods (e.g., `readSample()`, `readChain()`, 
    `readReport()`, `readProgress()`, `readRestart()`, ...).

**Minor enhancements**  

+   A bug in the y-axis variable names in the heatmap plot is now fixed.

### Version  2.1.1 -- October 9, 2020  

**Minor enhancements**  

+   A Linux bug in the installation of the MPI library is now fixed.

### Version  2.1.0 -- October 3, 2020  

**Minor enhancements**  

+   A new simulation specification `overwriteRequested` has 
    been added to all ParaMonte samplers. If `True` and the 
    ParaMonte sampler detects an existing set of old simulation 
    output files in the output path of the current simulation with 
    the same names as the output file names of the current simulation, 
    then, the ParaMonte sampler will overwrite the existing simulation files.  

### Version  2.0.9 -- October 2, 2020  

**Minor enhancements**  

+   Minor correction to the value of `__version__`, 
    now representing solely the version number.  

+   A simple example-usage Python script is now 
    added to the README.md file of the package.  

### Version  2.0.8 -- September 29, 2020  

**Minor enhancements**  

+   Enhanced error messages for situations when 
    the MPI library cannot be found on the system.  

### Version  2.0.7 -- September 26, 2020  

**Minor enhancements**  

+   The guidelines for the installation of the MPI 
    library on macOS and Linux have been improved.  

### Version  2.0.6 -- September 25, 2020  

**Minor enhancements**  

+   The explicit dependencies on `scipy`, `matplotlib`, and `seaborn` are 
    now removed from the PyPI setup file of the ParaMonte library as these 
    are only required for the post-processing and visualizations of the 
    simulation results. From now on, only `numpy` and `pandas` are the 
    minimally-required Python modules, and practically, only `numpy`.  

+   Two new functions `verifyDependencyVersion()` and `getDependencyVersion()`
    are now added to the library that can check for the existence of the 
    ParaMonte library's visualization dependencies and their required 
    minimum versions.  

+   The `seaborn` Python library has now decided to deprecate the `distplot()` 
    function. The corresponding visualization method in the ParaMonte library 
    has been now updated to a more appropriate name and underlying function.  

### Version  2.0.4 -- September 22, 2020  

**Minor enhancements**  

+   The output of the plotting functions is now stored as a list in 
    the `currentFig` temporary component of the visualization objects.
    This way, access to multiple individual objects on the active plot 
    is maintained instead of only the last object. Overall, this is a 
    minor change that will not cause any noticeable change in the 
    behavior of the library in almost in all use cases.

+   A minor bug regarding the input value for the `outputDelimiter` 
    attribute of the `spec` component of the `ParaMonteSampler()` class,  
    used in the `readTabular()` internal method, is now fixed.

### Version  2.0.3 -- September 11, 2020  

**Minor enhancements**  

+   Minor enhancement to `checkForUpdate()` method of 
    the `paramonte` module.  

### Version  2.0.2 -- September 11, 2020  

**Minor enhancements**  

+   Minor enhancement to `checkForUpdate()` method of 
    the `paramonte` module.  

### Version  2.0.1 -- September 10, 2020  

**Minor enhancements**  

+   LGPL3 LICENSE is now switched to MIT LICENSE.md file.  

+   A fix to the `brew` software installation now avoids 
    the seemingly-unavoidable crash.  

### Version  2.0.0 -- September 6, 2020  

**Major enhancements to the ParaMonte / ParaDRAM sampler interfaces**  

+   The entire ParaMonte Python interface library has been revamped.
    The new naming conventions, visualization, and computing tools 
    are significantly nicer to deal with and in some cases, orders 
    of magnitude faster than the previous major release.

+   The kernel density estimates and visualization tools are now on average 
    **100 times or more faster than the previous release of the library**.

+   Several new post-processing functionalities have now been added, such as
    the ability to seamlessly parse the contents of the output `*_report.txt`, 
    `*_restart.txt`, and `*_progress.txt` simulation files, in addition to the
    other output files (`*_sample.txt` and `*_chain.txt`) that could be parsed
    in the previous versions.

+   The new major release also includes 3D visualization tools, such as 3D 
    line, scatter, or line+scatter plots as well as fast 2D and 3D kernel 
    density estimate contour plotting tools.

**Minor enhancements**  

+   The simulation output files reading is now completely overhauled. In particular, 
    the output file reader methods are now capable of handling input file paths that 
    point to a directory. In such cases, it will search the input directory for files 
    matching the requested file name pattern. If no input file is provided to the file 
    reader methods, the current working directory will be search for the the potential 
    simulation files that match the requested pattern. 

+   The error-signaling behavior of the library now is very much controlled, that is, 
    upon code failure, it does not automatically shutdown the Python kernel in Jupyter 
    Notebooks. The library now simply throws an error message upon failing instead of 
    restarting the environment.  

+   The single value assignment to `spec.targetAcceptanceRate` component of a ParaDRAM 
    object is now properly handled. For example, the following code is valid as expected,  
    ```python  
    import paramonte as pm
    pmpd = pm.ParaDRAM()
    pmpd.spec.targetAcceptanceRate = 0.23 # this is now valid
    pmpd.spec.targetAcceptanceRate = [0.2, 0.3] # this is also valid, which limits the acceptance rate to the specified range
    ```  

+   The minimum required dependency versions are now raised to the following,  
    ```python  
    python_requires = ">=3.5"
    install_requires = [ "numpy>=1.18.0"
                       , "scipy>=1.4.0"
                       , "pandas>=1.0.0"
                       , "seaborn>=0.10.0"
                       , "matplotlib>=3.2.0"
                       ]
    ```  


## **Version 1.x.x**  

### Version  1.1.1 -- June 7, 2020  

**Minor enhancements**  

+   The `_ScatterLinePlot` dangling class is removed from the package.  

### Version  1.1.0 -- June 1, 2020  

+   Major enhancements to the ParaMonte kernel library.  
+   Major bug fixes in the ParaMonte Python library.  
+   The ParaMonte kernel and Python interface versions are now reposted separately as components of the paramonte module. 

### Version  1.0.12 -- April 6, 2020  

+   Minor enhancements and bug fixes to the kernel routines.

### Version  1.0.11 -- April 4, 2020  

+   Minor enhancements and bug fixes to the GridPlot.

### Version  1.0.10 -- March 28, 2020  

+   Minor bug fix.

### Version  1.0.9 -- March 27, 2020  

+   Minor enhancements.

### Version  1.0.8 -- March 27, 2020  

+   Minor enhancements.

### Version  1.0.7 -- March 26, 2020  

+   Minor corrections.

### Version  1.0.6 -- March 22, 2020  

+   Minor bug fixes.

### Version  1.0.5 -- March 21, 2020  

+   Minor bug fix.

### Version  1.0.4 -- March 20, 2020  

+   support for macOS (Darwin) added.

### Version  1.0.3 -- February 13, 2020  

+   Minor bug fixes.

### Version  1.0.2 -- February 13, 2020  

+   Minor bug fixes to the parallel routines.

### Version  1.0.1 -- February 13, 2020  

+   Minor bug fixes.

### Version  1.0.0 -- January 1, 2020 -- Initial release  

+   This is the first public release of the ParaMonte library.  

**New features**  
+   ParaDRAM sampler: **Para**llel **D**elayed-**R**ejection **A**daptive Metropolis-Hastings **M**arkov Chain Monte Carlo Sampler.  
+   ParaMonte Interface to the Python Programming languages.  
+   ParaMonte simulation-output visualization via the ParaMonte Python interface.  

  
<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank"><img src="https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png" alt="ParaMonte: Plain Powerful Parallel Monte Carlo Library" /></a>
<br><br>
<a href="https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>  
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a> 
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a> 
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a> 
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=pypi%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>  
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a> 
</div>
  

> MatDRAM is part of the ParaMonte library.  


ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions 
of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, 
and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), 
**user-friendliness** (of the library), **accessibility** (from multiple programming environments), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  
  
  
ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  

-   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray 
    communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations up to 16 digits of precision as requested by the user.  

-   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in 
    external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

  
Quick start  
===========  

For a quick start with some MATLAB Live Script examples, visit [this ParaMonte documentation page](https://www.cdslab.org/paramonte/notes/examples/matlab/mlx/). 
The corresponding example source files (the `*.mlx` files) can be downloaded from the [paramonte.svg?branch=main GitHub repository](https://github.com/cdslaborg/paramonte.svg?branch=main/tree/main/MATLAB/mlx), 
a repository dedicated to the ParaMonte library examples.  

The following example code samples a 4-dimensional MultiVariate Normal (MNV) distribution via the ParaDRAM sampler in serial mode,  

```matlab  
addpath(genpath("./"),"-begin") % change this path to the root directory of paramonte
pm = paramonte();
pmpd = pm.ParaDRAM();
getLogFunc = @(x) -0.5 * sum( x.^2 );
pmpd.runSampler ( 4 ... assume a 4-dimensional objective function
                , getLogFunc ...           the objective function
                );
```  

To learn about the post-processing and visualization tools of the `ParaMonte::MATLAB` library, visit [this this documentation page](https://www.cdslab.org/paramonte/notes/examples/matlab/mlx/).  

  
Installation  
============  

The latest release of the ParaMonte MatDRAM library can be downloaded from the release page of the library's repository on GitHub:  

[https://github.com/cdslaborg/paramonte/releases/latest/](https://github.com/cdslaborg/paramonte/releases/latest/)  

Alternatively, you can build the library from the source in the GitHub repository of the project: https://github.com/cdslaborg/paramonte  
For instructions, please visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

  
Dependencies  
============  

None.

  
Citing ParaMonte  
================  

The ParaMonte library is an honor-ware and its currency is acknowledgment and citations.  
  
As per the ParaMonte library license agreement terms, if you use any parts of this library for any purposes, 
kindly acknowledge the use of the ParaMonte library in your work (education/research/industry/development/...) 
by citing the ParaMonte library's main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) 
to access the PDF version of these files free of charge.  

  
License  
=======  

[MIT License](https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md)  

**What does this license mean?**  

Essentially, all we are asking from the users or developers is to  

>   explicitly acknowledge the use of this library or any concepts or parts of it in their education, research, or software (free or commercial).  

This is a free software, so help us keep it freely available to the public by redistributing the library and contributing to it. 
If you have questions or concerns about the license, do not hesitate to contact us (shahmoradi@utexas.edu).  

  
Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, 
        data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, 
        high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, 
        biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, 
        protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+   [Shashank Kumbhare](https://www.cdslab.org/people/#shashank-kumbhare)  
    +   physicist / Computational Data Scientist,  
    +   currently a UTA Physics member,  
    +   contact: [shashankkumbhare8@gmail.com](mailto:"shashankkumbhare8@gmail.com")  

  
Example usage instructions  
==========================  

+   **Install a MATLAB >2017b distribution**, preferably, the the latest MATLAB. 
    Note that ParaMonte MATLAB library have been tested only with MATLAB version 2018b and newer.  

+   **Running the ParaMonte simulations**  
    For complete up-to-date detailed instructions, visit: https://www.cdslab.org/paramonte/notes/run/matlab/
    +   Open the MATLAB software.  
    +   Suppose your mathematical objective function is a multivariate Normal distribution as implemented in this 
        [logfunc.m](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/mvn/MATLAB/logfunc.m) file.  
    +   Download this example generic  
        [main.m](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main.m) 
        MATLAB main file and save it in the same folder containing the `logfunc.m` file that you downloaded in the above. 
        Then, simply type the name of this MATLAB main script, `main` on the MATLAB command prompt.  
    +   Enjoy the unification of simplicity, efficiency, and comprehensive reporting in Monte Carlo simulations!  
    +   The ParaMonte library samplers are extremely versatile with many adjustable input parameters. 
        To learn about the many advanced features of the ParaMonte samplers, visit: https://www.cdslab.org/paramonte  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

  
<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank"><img src="https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png" alt="ParaMonte: Plain Powerful Parallel Monte Carlo Library" /></a>
<br><br>
<a href="https://github.com/cdslaborg/paramonte#license" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://travis-ci.com/cdslaborg/paramonte" target="_blank"><img src="https://travis-ci.com/cdslaborg/paramonte.svg?branch=main&style=flat-square" alt="Build Status" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/serial/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-serial%20kernel%20:%2098.2%25-brightgreen?style=flat-square" alt="kernel code coverage - serial" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/mpi/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-MPI%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - MPI" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/caf/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-Coarray%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - Coarray" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=PyPI%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/badge/dynamic/json?style=flat-square&labelColor=grey&color=brightgreen&maxAge=86400&label=PyPI%20downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fparamonte" alt="PyPI - Downloads Total" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
<br><br>
<a href="#paramonte-plain-powerful-parallel-monte-carlo-library">Overview</a> | 
<a href="#installation">Installation</a> | 
<a href="#dependencies">Dependencies</a> | 
<a href="#parallelism">Parallelism</a> | 
<a href="#example-usage-instructions">Examples</a> |
<a href="#citing-paramonte">Acknowledgments</a> | 
<a href="#license">License</a> | 
<a href="#authors-and-contributors">Authors</a>  
</div>
  
  
ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions 
of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, 
and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), 
**user-friendliness** (of the library), **accessibility** (from multiple programming environments), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  
  
  
ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  

-   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray 
    communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations up to 16 digits of precision as requested by the user.  

-   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in 
    external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

  
Quick start  
===========  

For a quick start with some MATLAB Live Script examples, visit [this ParaMonte documentation page](https://www.cdslab.org/paramonte/notes/examples/matlab/mlx/). 
The corresponding example source files (the `*.mlx` files) can be downloaded from the [paramontex GitHub repository](https://github.com/cdslaborg/paramontex/tree/main/MATLAB/mlx), 
a repository dedicated to the ParaMonte library examples.  

The following example code samples a 4-dimensional MultiVariate Normal (MNV) distribution via the ParaDRAM sampler in serial mode,  

```matlab  
addpath(genpath("./"),"-begin") % change this path to the root directory of paramonte
pm = paramonte();
pmpd = pm.ParaDRAM();
getLogFunc = @(x) -0.5 * sum( x.^2 );
pmpd.runSampler ( 4 ... assume a 4-dimensional objective function
                , getLogFunc ...           the objective function
                );
```  

To learn about the post-processing and visualization tools of the `ParaMonte::MATLAB` library, visit [this this documentation page](https://www.cdslab.org/paramonte/notes/examples/matlab/mlx/).  

  
Installation  
============  

+   **Windows**  
    The latest release of the ParaMonte MATLAB library can be downloaded from the release page of the library's repository on GitHub:  
    [https://github.com/cdslaborg/paramonte/releases/latest/](https://github.com/cdslaborg/paramonte/releases/latest/)  

+   **Linux**  
    The latest release of the ParaMonte MATLAB library can be downloaded from the release page of the library's repository on GitHub:  
    [https://github.com/cdslaborg/paramonte/releases/latest/](https://github.com/cdslaborg/paramonte/releases/latest/)  
    Alternatively, you can download the library to your local system directly by calling the `wget` Linux application from the command line,  
    ```bash  
    libname=libparamonte_matlab_linux_x64
    wget https://github.com/cdslaborg/paramonte/releases/latest/download/$libname.tar.gz
    tar xvzf $libname.tar.gz && cd $libname
    matlab # run matlab from the command line, then call the supplied "main" example script 
    ```  

+   **macOS (darwin)**  
    We **strongly advise you** to download the ParaMonte library for macOS via the following commands in a `bash` / `zsh` terminal, 
    instead of downloading the library directly from the GitHub release page,  
    ```bash  
    libname=libparamonte_matlab_darwin_x64
    curl -OL https://github.com/cdslaborg/paramonte/releases/latest/download/$libname.tar.gz
    tar xvzf $libname.tar.gz && cd $libname
    matlab # run matlab from the command line, then call the supplied "main" example script 
    ```  

  
Dependencies  
============  

The serial version of the ParaMonte MATLAB library kernel library has **NO external library dependencies**. 
The parallel version requires an MPI runtime library.  

  
Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. 
To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. 
In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). 
These automatic checks and installations happen when you download and install or use the library on your system, for the first time. 
If the automatic installation is unsuccessful, you can also install the libraries manually on your system:  

+   On **Windows** and **Linux** operating systems, we highly recommend downloading and installing the 
    [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), 
    which is available to the public free of charge, also available in the latest release of the 
    ParaMonte library on the [GitHub release page](https://github.com/cdslaborg/paramonte/releases) 
    (For Windows, look for the executable file that ends with `.exe`. For Linux, look for the file 
    that ends with `.tgz`, like `l_mpi-rt_2018.2.199.tgz`).
+   On **macOS**, the Intel MPI library is not available. Therefore, we recommend installing either 
    [Open-MPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/) MPI runtime libraries 
    depending the prebuilt version of the ParaMonte library that you have downloaded or 
    the configuration with which you intend to build the library.  

For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

  
Example usage instructions  
==========================  

+   **Install a MATLAB >2017b distribution**, preferably, the the latest MATLAB. 
    Note that ParaMonte MATLAB library have been tested only with MATLAB version 2018b and newer.  

+   **Optionally install a compatible MPI library** (or let the ParaMonte library take care of the MPI installation 
    when you call the library for the first time). For parallel simulations (via MPI), you will need an MPI library 
    already installed on your system. If you choose to install the library by yourself, we recommend the Intel MPI 
    library which is available for free from the Intel website or from [the ParaMonte GitHub release page](https://github.com/cdslaborg/paramonte/releases/tag/v1.5.1).  
    On macOS, the OpenMPI (or MPICH) MPI library can be used in place of the Intel MPI library which currently does not support macOS.  

+   **Calling the ParaMonte library for the first time**  
    +   Depending on your platform,  
        +   **Windows**  
            Nothing special needs to be done. You are all set! Follow the instructions below on how to call the ParaMonte library for the first time.  
        +   **Linux/macOS**  
            Open a `Bash` or `zsh` terminal and open MATLAB from the command line by calling its name,  
            ```bash  
            matlab
            ```  
            If `matlab` is not recognized on your command line as an application, seek help from 
            [this ParaMonte documentation page](https://www.cdslab.org/paramonte/notes/troubleshooting/bash-matlab-command-not-found/). 
    +   Once the MATLAB interactive environment opens, navigate to the root folder of the ParaMonte library (where the LICENSE file exists) 
        and call the ParaMonte library for the first time via the following commands (simply type the commands on the MATLAB command prompt),  
        ```matlab  
        addpath(genpath("./"),"-begin"); % add the ParaMonte library directories to MATLAB's list of search paths.
        pm = paramonte(); % instantiate an object of class paramonte.
        pm.verify(); % verify the integrity of the ParaMonte library on your system.
        ```  
        If needed, follow any extra instructions provided by the library on your MATLAB command prompt.  
        **If you do not intend to run simulations in parallel, you can say NO (`n`) to any MPI library installation requests via ParaMonte**.  
        If you do not intend to run simulations in parallel, answer YES (`y`) to any permission requests by 
        the ParaMonte library to install the MPI libraries on your system.  
+   **Running the ParaMonte simulations**  
    For complete up-to-date detailed instructions, visit: https://www.cdslab.org/paramonte/notes/run/matlab/
    +   Open the MATLAB software. On **Linux** and **macOS**, call the matlab executable from a Bash command line.  
    +   The ParaMonte library typically ships with example scripts. If you see a file named `main.m` at the root directory of 
        your ParaMonte library, then simply call this MATLAB script to run the example simulation provided with the library. 
        If no such file exists, or if you intend to do simulations in parallel, then follow the rest of the instructions below.  
    +   Suppose your mathematical objective function is a multivariate Normal distribution as implemented in this 
        [logfunc.m](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/mvn/MATLAB/logfunc.m) file.  
    +   For **serial** simulations, download this example generic serial 
        [main.m](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main.m) 
        MATLAB main file and save it in the same folder containing the `logfunc.m` file that you downloaded in the above. 
        Then, simply type the name of this MATLAB main script, `main` on the MATLAB command prompt.  
    +   For **parallel** simulations, download this example generic parallel 
        [main_mpi.m](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main_mpi.m) 
        MATLAB main file and save it in the same folder containing the `logfunc.m` file that you downloaded in the above. 
        Then, simply invoke the MPI launcher followed by the name of the MATLAB main script on a 
        MATLAB-aware MPI-aware Windows or Bash command prompt, similar to the following,  
        +   on **Windows** (preferably, on an Intel Parallel Studio command prompt or, the Microsoft Visual Studio's command prompt or, 
            some other command prompt that recognizes both `matlab` and the Intel's `mpiexec` software),  
            ```  
            mpiexec -localonly -n 3 matlab -batch main_mpi
            ```  
            where the `-localonly` flag is needed only if you are using the Intel MPI runtime libraries 
            (which is the default MPI library used to build the ParaMonte libraries on Windows).  
        +   on **Linux** or **macOS** (within a Bash terminal),  
            ```  
            mpiexec -n 3 matlab -batch main_mpi
            ```  
        Here, the parallel simulations are performed on 3 processes. Change the number 3 to any number of processes you wish to use, 
        but do not go beyond the maximum number of physical processes available on your system, otherwise, it will only degrade 
        the performance of your parallel simulations. For example, if you are running the parallel simulation on a personal 
        quad-cores laptop, set the number of processes to either 3 or 4 at most.  
    +   Enjoy the unification of simplicity, efficiency, and parallelism in Monte Carlo simulations!  
    +   The ParaMonte library samplers are extremely versatile with many adjustable input parameters. 
        To learn about the many advanced features of the ParaMonte samplers, visit: https://www.cdslab.org/paramonte  

  
Citing ParaMonte  
================  

The ParaMonte library is an honor-ware, the currency of which is acknowledgment and citations.  
  
If you use ParaMonte or any ideas from the software, please acknowledge it by citing the ParaMonte library's 
main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) 
to access the PDF version of these files free of charge.  

  
License  
=======  

[MIT License](https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md)  

**What does this license mean?**  

Essentially, all we are asking from the users or developers is to  

>   explicitly acknowledge the use of this library or any concepts or parts of it in their education, research, or software (free or commercial).  

This is a free software, so help us keep it freely available to the public by redistributing the library and contributing to it. 
If you have questions or concerns about the license, do not hesitate to contact us (shahmoradi@utexas.edu).  

  
Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, 
        data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, 
        high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, 
        biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, 
        protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+   [Shashank Kumbhare](https://www.cdslab.org/people/#shashank-kumbhare)  
    +   physicist / Computational Data Scientist,  
    +   currently a UTA Physics member,  
    +   contact: [shashankkumbhare8@gmail.com](mailto:"shashankkumbhare8@gmail.com")  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
# ParaMonte MATLAB release notes

This project follows [Semantic Versioning](https://semver.org/). 
To access the latest release of the package, visit [the ParaMonte GitHub repository release page](https://github.com/cdslaborg/paramonte/releases) or [the ParaMonte page on MathWorks FileExchange central package repository](https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte).  

## **Version 2.x.x**  

### Version  2.5.0 -- January 1, 2021  

**Major enhancements**  

+   This release is a a major step toward further portability 
    of the kernel routines of the `ParaMonte::MATLAB` library. The kernel 
    library dependencies are now properly handled and recognized at runtime
    without such aggressive actions as permanently redefining the environmental
    path variables, most importantly, `PATH` and `LD_LIBRARY_PATH` on Linux/macOS.

+   The `ParaMonte::MATLAB` library is now capable of recognizing the existing MPI
    libraries such as MPICH and OpenMPI on user's system and avoid further 
    installation of a new MPI library if it is deemed unnecessary.  

+   The ParaMonte kernel routines are now capable of handling user-input 
    file paths that contain white-space (blank) or other exotic characters.  

**Minor enhancements**  

+   As of this version, when the `ParaMonte::MATLAB` library is called in MATLAB 
    `-batch` mode (from the command line) for the first time, the library avoids 
    asking the user's response to the question of installing an MPI library if it 
    is missing on the user's system. This will prevent undesired crashes of the 
    simulations for the first time when the simulation is run from outside the 
    MATLAB session. However, the onus will be on the user to ensure an MPI 
    library exists on the system if they intend to run simulations in parallel.  

+   The `ParaMonte::MATLAB` library packages for different Operating systems and 
    processor architecture are now separate from each other. This change was made 
    to lower the overall size of `ParaMonte::MATLAB` by only keeping the relevant 
    files in each packaging of the library. The current release contains three 
    separate packages for `ParaMonte::MATLAB`,  
    +   `libparamonte_matlab_windows_x64`,  
    +   `libparamonte_matlab_darwin_x64`,  
    +   `libparamonte_matlab_linux_x64`.  

+   Typo-fixes in the documentation of the library.  

**MATLAB versions used for this release**  

+   **Windows**: `MATLAB 9.6.0.1072779 (R2019a)`  
+     **Linux**: `MATLAB 9.8.0.1323502 (R2020a)`  
+     **macOS**: `MATLAB 9.8.0.1323502 (R2020a)`  

**MATLAB version compatibility**  

This release has been tested with MATLAB 2018, 2019, and 2020. 
It should be also compatible with MATLAB 2017, but is not tested.
If you notice an incompatibility with any of the above MATLAB versions,
please report this issue to the developers for a resolution at:

https://github.com/cdslaborg/paramonte/issues

### Version  2.4.0 -- December 23, 2020  

+   This version of the library was internal to 
    the developers and not released to the public.

### Version  2.3.0 -- December 17, 2020  

**Major enhancements**  

+   This update presents several major performance, accuracy, 
    and verification enhancements to the ParaMonte kernel routines, 
    in particular, to the ParaDRAM sampler.  

+   An extensive set of over 866 tests have been added 
    that test all aspects of the ParaMonte kernel library.  

+   The issue of Windows file locking, that led to the occasional crashes of the 
    ParaDRAM and ParaDISE simulations in `multiChain` parallelism mode, is now resolved.  

+   The `ParaDRAM` class in `paramonte` is now also available 
    as `Paradram` and `paradram`, although the original label 
    will remain the default preferred method of ParaDRAM 
    object instantiation.  

### Version  2.2.1 -- November 15, 2020  

**Minor enhancements**  

+   Minor enhancements to the Kernel library 
    build scripts and dependencies management.  

+   More informative error messages are now printed 
    on MATLAB console if any error happens during the 
    ParaMonte library setup on macOS for the first time.  

### Version  2.2.0 -- October 29, 2020  

**Enhancements**  

+   The `cmake` software dependency installation failure now 
    does not nullify the installation of other dependencies.

+   The IO debugging info of all ParaMonte samplers have been enhanced. 
    In cases of wrong syntax or syntax-breaking input values in the simulation 
    output files, the error messages are now more informative and point directly 
    to the exact location of of error in the input file.  

+   The Integrated Autocorrelation (IAC) for sample refinement in ParaDRAM 
    sampler of ParaMonte is now set to the average of all variables' IAC values 
    instead of the maximum IAC value. This will lead to less aggressive decorrelation 
    of the final sample, which means significantly larger final sample sizes, without 
    compromising the i.i.d. property of the final refined sample. This behavior can 
    be reversed back to the original by specifying "max" or "maximum" along with 
    the requested refinement method, `SampleRefinementMethod = "batchmeans max"` 
    or `SampleRefinementMethod = "BatchMeans-max"` (case-insensitive).

### Version  2.1.3 -- October 15, 2020  

**Minor enhancements**  

+   Further minor enhancements to the behavior of the 
    `checkForUpdate()` method of the `paramonte` class.  

### Version  2.1.2 -- October 15, 2020  

**Minor enhancements**  

+   The `checkForUpdate()` method of the `paramonte` 
    class now functions as expected.  

### Version  2.1.1 -- October 9, 2020  

**Minor enhancements**  

+   A Linux bug in the installation of the MPI library is now fixed.

### Version  2.1.0 -- October 3, 2020  

**Minor enhancements**  

+   A new simulation specification `overwriteRequested` has 
    been added to all ParaMonte samplers. If `True` and the 
    ParaMonte sampler detects an existing set of old simulation 
    output files in the output path of the current simulation with 
    the same names as the output file names of the current simulation, 
    then, the ParaMonte sampler will overwrite the existing simulation files.  

### Version  2.0.1 -- September 26, 2020  

**Minor enhancements**  

+   The guidelines for the installation of the 
    MPI library on macOS have been improved.  

+   The minor bug in GridPlot class method `rotateAxesLabels()` that caused 
    the `readSample()` , `readChain()`, `readMarkovChain()` to crash upon adding 
    Grid plots is now fixed.  

+   The minor bug in the naming of the ParaMonte kernel library files on macOS 
    (Darwin) is now fixed.  

### Version  2.0.0 -- September 22, 2020  

**Major enhancements to the ParaMonte / ParaDRAM sampler interfaces**  

+   The entire ParaMonte MATLAB interface library has been revamped.
    The new naming conventions, visualization, and computing tools 
    are significantly nicer to deal with and in some cases, orders 
    of magnitude faster than the previous major release.

+   The simulation output files reading is now completely overhauled. In particular, 
    the output file reader methods are now capable of handling input file paths that 
    point to a directory. In such cases, it will search the input directory for files 
    matching the requested file name pattern. If no input file is provided to the file 
    reader methods, the current working directory will be search for the the potential 
    simulation files that match the requested pattern. 

+   Several new post-processing functionalities have now been added, such as
    the ability to seamlessly parse the contents of the output `*_report.txt`, 
    `*_restart.txt`, and `*_progress.txt` simulation files, in addition to the
    other output files (`*_sample.txt` and `*_chain.txt`) that could be parsed
    in the previous versions.

+   The newly-added `readRestart()` method is now added to the ParaDRAM sampler 
    class. User can now parse the contents of the output ASCII-format restart files. 
    This is particularly useful to visualize the dynamics of the ParaDRAM sampler class, 
    such as the evolution of the proposal distribution's location, shape, and covariance 
    matrix.  

+   The `GridPlot()` class now has two additional methods `setAxesLabels()` and 
    `setAxesLimits()` which can directly set the labels and limits of axes, hassle-free.

**Minor enhancements**  

+   The single value assignment to `spec.targetAcceptanceRate` component of a ParaDRAM object is now properly handles. 
    For example, the following code is valid as expected,  
    ```matlab  
    import paramonte as pm
    pmpd = pm.ParaDRAM()
    pmpd.spec.targetAcceptanceRate = 0.23 # this is now valid
    pmpd.spec.targetAcceptanceRate = [0.2, 0.3] # this is also valid, which limits the acceptance rate to the specified range
    ```  

+   The default background color in all plots is now `"white"`.  
+   The `rotateAxisLabels()` of the `GridPlot()` class is now renamed to `rotateAxesLabels()`.  

**Bug fixes**  

+   ParaDRAM `readMarkovChain()` no-output-option bug is now fixed. 
    When calling `readMarkovChain()`, user can now either provide the output variable or not.  

## **Version 1.x.x**  

### Version  1.1.0 -- June 5, 2020  

+   Enhancements and bug fixes to the kernel routines.  
+   Several major enhancements and bug fixes to the MATLAB kernel and interface routines.  
+   MatDRAM now supports fully-deterministic restart functionality.  

### Version  1.0.0 -- June 1, 2020 -- Initial release  

+   This is the first public release of the ParaMonte MATLAB library.  

**New features**  

+   ParaDRAM sampler: **Para**llel **D**elayed-**R**ejection **A**daptive Metropolis-Hastings **M**arkov Chain Monte Carlo Sampler.  
+   ParaMonte Interface to the MATLAB Programming language.  
+   ParaMonte simulation-output visualization via the ParaMonte MATLAB interface.  
export_fig
==========

A toolbox for exporting figures from MATLAB to standard image and document formats nicely.

### Overview
Exporting a figure from MATLAB the way you want it (hopefully the way it looks on screen), can be a real headache for the unitiated, thanks to all the settings that are required, and also due to some eccentricities (a.k.a. features and bugs) of functions such as `print`. The first goal of export_fig is to make transferring a plot from screen to document, just the way you expect (again, assuming that's as it appears on screen), a doddle.
  
The second goal is to make the output media suitable for publication, allowing you to publish your results in the full glory that you originally intended. This includes embedding fonts, setting image compression levels (including lossless), anti-aliasing, cropping, setting the colourspace, alpha-blending and getting the right resolution.

Perhaps the best way to demonstrate what export_fig can do is with some examples.

*Note: `export_fig` currently supports only figures created with the `figure` function, or GUIDE. Figures created using `uifigure` or AppDesigner are only partially supported. See issues [#287](https://github.com/altmany/export_fig/issues/287), [#261](https://github.com/altmany/export_fig/issues/261) for details.*
  
### Examples
**Visual accuracy** - MATLAB's exporting functions, namely `saveas` and `print`, change many visual properties of a figure, such as size, axes limits and ticks, and background colour, in unexpected and unintended ways. Export_fig aims to faithfully reproduce the figure as it appears on screen. For example:  
```Matlab
plot(cos(linspace(0, 7, 1000)));
set(gcf, 'Position', [100 100 150 150]);
saveas(gcf, 'test.png');
export_fig test2.png
```
generates the following:

| Figure: | test.png: | test2.png: |
|:-------:|:---------:|:----------:|
|![](https://farm6.staticflickr.com/5616/15589249291_16e485c29a_o_d.png)|![](https://farm4.staticflickr.com/3944/15406302850_4d2e1c7afa_o_d.png)|![](https://farm6.staticflickr.com/5607/15568225476_8ce9bd5f6b_o_d.png)|

Note that the size and background colour of test2.png (the output of export_fig) are the same as those of the on screen figure, in contrast to test.png. Of course, if you want the figure background to be white (or any other colour) in the exported file then you can set this prior to exporting using:
```Matlab
set(gcf, 'Color', 'w');
```
  
Notice also that export_fig crops and anti-aliases (smooths, for bitmaps only) the output by default. However, these options can be disabled; see the Tips section below for details.
  
**Resolution** - by default, export_fig exports bitmaps at screen resolution. However, you may wish to save them at a different resolution. You can do this using either of two options: `-m<val>`, where <val> is a positive real number, magnifies the figure by the factor <val> for export, e.g. `-m2` produces an image double the size (in pixels) of the on screen figure; `-r<val>`, again where <val> is a positive real number, specifies the output bitmap to have <val> pixels per inch, the dimensions of the figure (in inches) being those of the on screen figure. For example, using:  
```Matlab
export_fig test.png -m2.5
```
on the figure from the example above generates:

![](https://farm4.staticflickr.com/3937/15591910915_dc7040c477_o_d.png)

Sometimes you might have a figure with an image in. For example:
```Matlab
imshow(imread('cameraman.tif'))
hold on
plot(0:255, sin(linspace(0, 10, 256))*127+128);
set(gcf, 'Position', [100 100 150 150]);
```
generates this figure:

![](https://farm4.staticflickr.com/3942/15589249581_ff87a56a3f_o_d.png)
  
Here the image is displayed in the figure at resolution lower than its native resolution. However, you might want to export the figure at a resolution such that the image is output at its native (i.e. original) size (in pixels). Ordinarily this would require some non-trivial computation to work out what that resolution should be, but export_fig has an option to do this for you. Using:
```Matlab
export_fig test.png -native
```
produces:

![](https://farm6.staticflickr.com/5604/15589249591_da2b2652e4_o_d.png)

with the image being the size (in pixels) of the original image. Note that if you want an image to be a particular size, in pixels, in the output (other than its original size) then you can resize it to this size and use the `-native` option to achieve this.

All resolution options (`-m<val>`, `-q<val>` and `-native`) correctly set the resolution information in PNG and TIFF files, as if the image were the dimensions of the on screen figure.

**Shrinking dots & dashes** - when exporting figures with dashed or dotted lines using either the ZBuffer or OpenGL (default for bitmaps) renderers, the dots and dashes can appear much shorter, even non-existent, in the output file, especially if the lines are thick and/or the resolution is high. For example:  
```Matlab
plot(sin(linspace(0, 10, 1000)), 'b:', 'LineWidth', 4);
hold on
plot(cos(linspace(0, 7, 1000)), 'r--', 'LineWidth', 3);
grid on
export_fig test.png
```
generates:

![](https://farm4.staticflickr.com/3956/15592747732_f943d4aa0a_o_d.png)

This problem can be overcome by using the painters renderer. For example:
```Matlab
export_fig test.png -painters
```
used on the same figure generates:

![](https://farm4.staticflickr.com/3945/14971168504_77692f11f5_o_d.png)

Note that not only are the plot lines correct, but the grid lines are too.

**Transparency** - sometimes you might want a figure and axes' backgrounds to be transparent, so that you can see through them to a document (for example a presentation slide, with coloured or textured background) that the exported figure is placed in. To achieve this, first (optionally) set the axes' colour to 'none' prior to exporting, using:  
```Matlab
set(gca, 'Color', 'none'); % Sets axes background
```
    
then use export_fig's `-transparent` option when exporting:
```Matlab
export_fig test.png -transparent
```

This will make the background transparent in PDF, EPS and PNG outputs. You can additionally save fully alpha-blended semi-transparent patch objects to the PNG format. For example:

```Matlab
logo;
alpha(0.5);
```

generates a figure like this:

![](https://farm4.staticflickr.com/3933/15405290339_b08de33528_o_d.png)

If you then export this to PNG using the `-transparent` option you can then put the resulting image into, for example, a presentation slide with fancy, textured background, like so:

![](https://farm6.staticflickr.com/5599/15406302920_59beaefff1_o_d.png)

and the image blends seamlessly with the background.

**Image quality** - when publishing images of your results, you want them to look as good as possible. By default, when outputting to lossy file formats (PDF, EPS and JPEG), export_fig uses a high quality setting, i.e. low compression, for images, so little information is lost. This is in contrast to MATLAB's print and saveas functions, whose default quality settings are poor. For example:
```Matlab
A = im2double(imread('peppers.png'));
B = randn(ceil(size(A, 1)/6), ceil(size(A, 2)/6), 3) * 0.1;
B = cat(3, kron(B(:,:,1), ones(6)), kron(B(:,:,2), ones(6)), kron(B(:,:,3), ones(6)));
B = A + B(1:size(A, 1),1:size(A, 2),:);
imshow(B);
print -dpdf test.pdf
```
generates a PDF file, a sub-window of which looks (when zoomed in) like this:

![](https://farm6.staticflickr.com/5613/15405290309_881b2774d6_o_d.png)

while the command

```Matlab
export_fig test.pdf
```
on the same figure produces this:

![](https://farm4.staticflickr.com/3947/14971168174_687473133f_o_d.png)

While much better, the image still contains some compression artifacts (see the low level noise around the edge of the pepper). You may prefer to export with no artifacts at all, i.e. lossless compression. Alternatively, you might need a smaller file, and be willing to accept more compression. Either way, export_fig has an option that can suit your needs: `-q<val>`, where <val> is a number from 0-100, will set the level of lossy image compression (again in PDF, EPS and JPEG outputs only; other formats are lossless), from high compression (0) to low compression/high quality (100). If you want lossless compression in any of those formats then specify a <val> greater than 100. For example:
```Matlab
export_fig test.pdf -q101
```
again on the same figure, produces this:

![](https://farm6.staticflickr.com/5608/15405803908_934512c1fe_o_d.png)

Notice that all the noise has gone.

### Tips
**Anti-aliasing** - the anti-aliasing which export_fig applies to bitmap outputs by default makes the images look nice, but it can also blur images and increase exporting time and memory requirements, so you might not always want it. You can set the level of anti-aliasing by using the `-a<val>` option, where <val> is 1 (no anti-aliasing), 2, 3 (default) or 4 (maximum anti-aliasing).  
  
**Cropping** - by default, export_fig crops its output to minimize the amount of empty space around the figure. If you'd prefer the figure to be uncropped, and instead have the same appearance (in terms of border width) as the on screen figure, then use the `-nocrop` option.  
  
**Colourspace** - by default, export_fig generates files in the RGB [colourspace](https://en.wikipedia.org/wiki/Color_space). However, you can also export in greyscale or the CMYK colourspace, using the `-grey` (or `-gray`) and `-cmyk` options respectively. The CMYK option is useful for publishers who require documents in this colourspace, but the option is only supported for PDF, EPS and TIFF files.

**Specifying a target directory** - you can get export_fig to save output files to any directory (for which you have write permission), simply by specifying the full or relative path in the filename. For example:
```Matlab
export_fig ../subdir/fig.png;
export_fig('C:/Users/Me/Documents/figures/myfig', '-pdf', '-png');
```

**Variable file names** - often you might want to save a series of figures in a for loop, each with a different name. For this you can use the functional form of input arguments, i.e. `export_fig(arg1, arg2)`,  and construct the filename string in a variable. Here's an example of this:  
```Matlab
for a = 1:5
    plot(rand(5, 2));
    export_fig(sprintf('plot%d.png', a));
end
```
When using the functional form like this, be sure to put string variables in quotes:
```Matlab
export_fig(sprintf('plot%d', a), '-a1', '-pdf', '-png');
```

**Specifying the figure/axes** - if you have multiple figures open you can specify which figure to export using its handle:  
```Matlab
export_fig(figure_handle, filename);
```
Equally, if your figure contains several subplots then you can export just one of them by giving export_fig the handle to the relevant axes:
```Matlab
export_fig(axes_handle, filename);
```

**Multiple formats** - save time by exporting to multiple formats simultaneously. E.g.: 
```Matlab
export_fig filename -pdf -eps -png -jpg -tiff
```

**Other file formats** - if you'd like to save your figure to a bitmap format that is not supported by export_fig, e.g. animated GIF, PPM file or a frame in a movie, then you can use export_fig to output the image, and optionally an alpha-matte, to the workspace. E.g.:  
```Matlab
frame = export_fig;
```
or
```Matlab
[frame, alpha] = export_fig;
```
These variables can then be saved to other image formats using other functions, such as imwrite.

**Appending to a file** - you can use the `-append` option to append the figure to the end of an image/document, if it already exists. This is supported for PDF and TIFF files only. Note that if you wish to append a lot of figures consecutively to a PDF, it can be more efficient to save all the figures to PDF separately then append them all in one go at the end (e.g. using [append_pdfs](http://www.mathworks.com/matlabcentral/fileexchange/31215-appendpdfs)).  
  
**Output to clipboard** - you can use the `-clipboard` option to copy the specified figure or axes to the system clipboard, for easy paste into other documents (e.g., Word or PowerPoint). Note that the image is copied in bitmap (not vector) format.  
  
**Font size** - if you want to place an exported figure in a document with the font a particular size then you need to set the font to that size in the figure, and not resize the output of export_fig in the document. To avoid resizing, simply make sure that the on screen figure is the size you want the output to be in the document before exporting.  
  
**Renderers** - MATLAB has three renderers for displaying and exporting figures: painters, OpenGL and ZBuffer. The different renderers have different [features](http://www.mathworks.com/access/helpdesk/help/techdoc/creating_plots/f3-84337.html#f3-102410), so if you aren't happy with the result from one renderer try another. By default, vector formats (i.e. PDF and EPS outputs) use the painters renderer, while other formats use the OpenGL renderer. Non-default renderers can be selected by using one of these three export_fig input options: `-painters`, `-opengl`, `-zbuffer`:  
```Matlab
export_fig test.png -painters
```
  
**Artifacts** - sometimes the output that you get from export_fig is not what you expected. If an output file contains artifacts that aren't in the on screen figure then make sure that the renderer used for rendering the figure on screen is the same as that used for exporting. To set the renderer used to display the figure, use:  
```Matlab
set(figure_handle, 'Renderer', 'opengl');
```
After matching the two renderers, if the artifact appears in the on screen figure then you'll need to fix that before exporting. Alternatively you can try changing the renderer used by export_fig. Finally check that it isn't one of the known issues mentioned in the section below.

**Smoothed/interpolated images in output PDF** - if you produce a PDF using export_fig and images in the PDF look overly smoothed or interpolated, this is because the software you are using to view the PDF is smoothing or interpolating the image data. The image is not smoothed in the PDF file itself. If the software has an option to disable this feature, you should select it. Alternatively, use another PDF viewer that doesn't exhibit this problem.  
  
**Locating Ghostscript/pdftops** - You may find a dialogue box appears when using export_fig, asking you to locate either [Ghostscript](http://www.ghostscript.com) or [pdftops (part of the Xpdf package)](http://www.xpdfreader.com). These are separate applications which export_fig requires to perform certain functions. If such a dialogue appears it is because export_fig can't find the application automatically. This is because you either haven't installed it, or it isn't in the normal place. Make sure you install the applications correctly first. They can be downloaded from the following places:  
 1. Ghostscript:     [www.ghostscript.com](http://www.ghostscript.com)
 2. pdftops (install the Xpdf package): [www.xpdfreader.com](http://www.xpdfreader.com)

If you choose to install them in a non-default location then point export_fig
to this location using the dialogue box.

**Undefined function errors** - If you download and run export_fig and get an error similar to this:  
```
??? Undefined function or method 'print2array' for input arguments of type 'double'.
```
then you are missing one or more of the files that come in the export_fig package. Make sure that you click the "Get from GitHub" button at the top-right of the download [page](http://www.mathworks.co.uk/matlabcentral/fileexchange/23629-exportfig), then extract all the files in the zip file to the same directory. You should then have all the necessary files.
  
### Known issues
There are lots of problems with MATLAB's exporting functions, especially `print`. Export_fig is simply a glorified wrapper for MATLAB's `print` function, and doesn't solve all of its bugs (yet?). Some of the problems I know about are:
  
**Fonts** - when using the painters renderer, MATLAB can only export a small number of fonts, details of which can be found [here](http://www.mathworks.com/help/releases/R2014a/matlab/creating_plots/choosing-a-printer-driver.html#f3-96545). Export_fig attempts to correct font names in the resulting EPS file (up to a maximum of 11 different fonts in one figure), but this is not always guaranteed to work. In particular, the text positions will be affected. It also does not work for text blocks where the 'Interpreter' property is set to 'latex'.

Also, when using the painters renderer, ghostscript will sometimes throw an error such as `Error: /undefined in /findfont`. This suggests that ghostscript could not find a definition file for one of your fonts. One possible fix for this is to make sure the file `EXPORT_FIG_PATH/.ignore/gs_font_path.txt` exists and contains a list of paths to the folder(s) containing the necessary font definitions (make sure that they are TrueType definitions!), separated by a semicolon.

**RGB color data not yet supported in Painter's mode** - you will see this as a warning if you try to export a figure which contains patch objects whose face or vertex colors are specified as an RGB colour, rather than an index into the colormap, using the painters renderer (the default renderer for vector output). This problem can arise if you use `pcolor`, for example. This is a problem with MATLAB's painters renderer, which also affects `print`; there is currently no fix available in export_fig (other than to export to bitmap). The suggested workaround is to avoid colouring patches using RGB. First, try to use colours in the figure's colourmap (instructions [here](http://www.mathworks.co.uk/support/solutions/en/data/1-6OTPQE/)) - change the colourmap, if necessary. If you are using `pcolor`, try using [uimagesc](http://www.mathworks.com/matlabcentral/fileexchange/11368) (on the file exchange) instead.  

**Dashed contour lines appear solid** - when using the painters renderer, MATLAB cannot generate dashed lines using the `contour` function (either on screen or in exported PDF and EPS files). Details can be found [here](http://www.mathworks.com/support/solutions/en/data/1-14PPHB/?solution=1-14PPHB).  
  
**Text size** - when using the OpenGL or ZBuffer renderers, large text can be resized relative to the figure when exporting at non-screen-resolution (including using anti-alising at screen resolution). This is a feature of MATLAB's `print `function. In this case, try using the `-painters` option.  
  
**Lighting and transparency** - when using the painters renderer, transparency and lighting effects are not supported. Sorry, but this is an inherent feature of MATLAB's painters renderer. To find out more about the capabilities of each rendering method, see [here](http://www.mathworks.com/access/helpdesk/help/techdoc/creating_plots/f3-84337.html#f3-102410). You can still export transparent objects to vector format (SVG) using the excellent [plot2svg](http://www.mathworks.com/matlabcentral/fileexchange/7401) package, then convert this to PDF, for example using [Inkscape](http://inkscape.org/). However, it can't handle lighting.  
  
**Lines in patch objects** - when exporting patch objects to PDF using the painters renderer (default), sometimes the output can appear to have lines across the middle of rectangular patches; these lines are the colour of the background, as if there is a crack in the patch, allowing you to see through. This appears to be due to bugs in MATLAB's internal vector rendering code. These lines can often be removed from the PDF using software such as [InkScape](https://inkscape.org). Sometimes disabling anti-aliasing in the PDF-reader software can get rid of the lines ([discussion](https://github.com/altmany/export_fig/issues/44)).  
  
**Out of memory** - if you run into memory issues when using export_fig, some ways to get round this are:  
 1. Reduce the level of anti-aliasing.
 2. Reduce the size of the figure.
 3. Reduce the export resolution (dpi). 
 4. Change the renderer to painters or ZBuffer.  
  
**Errors** - the other common type of errors people get with export_fig are OpenGL errors. This isn't a fault of export_fig, but either a bug in MATLAB's `print`, or your graphics driver getting itself into a state. Always make sure your graphics driver is up-to-date. If it still doesn't work, try using the ZBuffer renderer.  
  
### Raising issues
If you think you have found a genuine error or issue with export_fig **that is not listed above**, first ensure that the figure looks correct on screen when rendered using the renderer that export_fig is set to use (e.g. if exporting to PDF or EPS, does the figure look correct on screen using the painters renderer, or if exporting to bitmap, does the figure look correct on screen using the OpenGL renderer?). If it looks wrong then the problem is there, and I cannot help (other than to suggest you try exporting using a different renderer).

Secondly, if exporting to bitmap, do try all the renderers (i.e. try the options `-opengl`, `-zbuffer` and `-painters` separately), to see if one of them does produce an acceptable output, and if so, use that.

If this still does not help, then ensure that you are using the latest version of export_fig, which is available [here](https://github.com/altmany/export_fig/archive/master.zip).  
 
If the figure looks correct on screen, but an error exists in the exported output (which cannot be solved using a different renderer) then please feel free to raise an [issue](https://github.com/altmany/export_fig/issues). Please be sure to include the .fig file, the export_fig command you use, the output you get, and a description of what you expected. I can't promise anything, but if it's easy to fix I may indeed do it. Often I will find that the error is due to a bug in MATLAB's `print` function, in which case I will suggest you submit it as a bug to TheMathWorks, and inform me of any fix they suggest. Also, if there's a feature you'd like that isn't supported please tell me what it is and I'll consider implementing it.

### And finally...

![](https://farm4.staticflickr.com/3956/15591911455_b9008bd77e_o_d.jpg)

If you've ever wondered what's going on in the logo on the export_fig download page (reproduced here), then this explanantion is for you. The logo is designed to demonstrate as many of export_fig's features as possible: 
 
Given a figure containing a translucent mesh (top right), export_fig can export to pdf (bottom centre), which allows the figure to be zoomed-in without losing quality (because it's a vector graphic), but isn't able to reproduce the translucency. Also, depending on the PDF viewer program, small gaps appear between the patches, which are seen here as thin white lines. 
 
By contrast, when exporting to png (top left), translucency is preserved (see how the graphic below shows through), and the figure is anti-aliased. However, zooming-in does not reveal more detail since png is a bitmap format. Also, lines appear less sharp than in the pdf output.


  
<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank"><img src="https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png" alt="ParaMonte: Plain Powerful Parallel Monte Carlo Library" /></a>
<br><br>
<a href="https://github.com/cdslaborg/paramonte#license" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://travis-ci.com/cdslaborg/paramonte" target="_blank"><img src="https://travis-ci.com/cdslaborg/paramonte.svg?branch=main&style=flat-square" alt="Build Status" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/serial/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-serial%20kernel%20:%2098.2%25-brightgreen?style=flat-square" alt="kernel code coverage - serial" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/mpi/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-MPI%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - MPI" /></a>
<a href="https://cdslaborg.github.io/paramonte-codecov/kernel/caf/" target="_blank"><img src="https://img.shields.io/badge/code%20coverage-Coarray%20kernel%20:%2098.3%25-brightgreen?style=flat-square" alt="kernel code coverage - Coarray" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=PyPI%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/badge/dynamic/json?style=flat-square&labelColor=grey&color=brightgreen&maxAge=86400&label=PyPI%20downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fparamonte" alt="PyPI - Downloads Total" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
<br><br>
<a href="#paramonte-plain-powerful-parallel-monte-carlo-library">Overview</a> | 
<a href="#installation">Installation</a> | 
<a href="#dependencies">Dependencies</a> | 
<a href="#parallelism">Parallelism</a> | 
<a href="#example-usage-instructions">Examples</a> |
<a href="#citing-paramonte">Acknowledgments</a> | 
<a href="#license">License</a> | 
<a href="#authors-and-contributors">Authors</a>  
</div>
  
  
ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions 
of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, 
and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), 
**user-friendliness** (of the library), **accessibility** (from multiple programming environments), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  
  
  
ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  

-   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray 
    communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations up to 16 digits of precision as requested by the user.  

-   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in 
    external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

  
Installation  
============  

The pre-built ready-to-use libraries are available on [the release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases). 
Each prebuilt ParaMonte library automatically ships with a full-fledged set of example codes and build scripts.  

Alternatively, you can build the library from the source in the [GitHub repository of the project](https://github.com/cdslaborg/paramonte). 
The ParaMonte library installation/build process is fully automated for all of the supported programming languages. 
Currently, the following compiler suites are supported **for builds from source**:  
  
| Compiler Suite                    | Linux | macOS | Windows (64bit) |  
|----------------------------------:|:-----:|:-----:|:---------------:|  
| GNU Compiler Collection > 8.4     |&check;|&check;| &cross;         |  
| Intel Parallel Studio > 19.1.1    |&check;|&check;| &check;         |  

For more information and quick-start in the programming language of your choice, visit the [ParaMonte library homepage](https://www.cdslab.org/paramonte).  

  
Dependencies  
============  

Beyond an optional MPI runtime library for parallel simulations, the ParaMonte kernel has **zero dependency** on external third-party libraries or packages.  

  
Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. 
To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. 
In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). 
These automatic checks and installations happen when you download and install or use the library on your system, for the first time. 
If the automatic installation is unsuccessful, you can also install the libraries manually on your system:  

+   On **Windows** and **Linux** operating systems, we highly recommend downloading and installing the 
    [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), 
    which is available to the public free of charge, also available in the latest release of the 
    ParaMonte library on the [GitHub release page](https://github.com/cdslaborg/paramonte/releases) 
    (For Windows, look for the executable file that ends with `.exe`. For Linux, look for the file 
    that ends with `.tgz`, like `l_mpi-rt_2018.2.199.tgz`).
+   On **macOS**, the Intel MPI library is not available. Therefore, we recommend installing either 
    [Open-MPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/) MPI runtime libraries 
    depending the prebuilt version of the ParaMonte library that you have downloaded or 
    the configuration with which you intend to build the library.  

For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

  
Example usage instructions  
==========================  

+   For complete clear organized up-to-date instructions on the build process and the installation 
    of the ParaMonte library, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

## Quick start  

+   Go to the [release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases),  
+   Decide on the parallelism paradigm that you want to use: serial / MPI 
    (the Coarray Fortran implementation is not available as a prebuilt dynamic library),  
+   Decide on the Operating System (OS) on which you want to run the ParaMonte simulations: Windows / macOS / Linux,  
+   Learn about the naming convention used for the ParaMonte prebuilt libraries [here](https://www.cdslab.org/paramonte/notes/installation/readme/#naming-convention-used-for-paramonte-library-builds),  
+   Download the prebuilt ParaMonte library of your choice based on the decisions you have made in the above. 
    If you are not sure which prebuilt library is suitable for your needs, use the prebuilt library recommended 
    [here for Windows](https://www.cdslab.org/paramonte/notes/installation/windows/#using-the-prebuilt-paramonte-library), or 
    [here for Linux](https://www.cdslab.org/paramonte/notes/installation/linux/#using-the-prebuilt-paramonte-library), or 
    [here for macOS](https://www.cdslab.org/paramonte/notes/installation/macos/#using-the-prebuilt-paramonte-library).  
+   Each prebuilt library ships with a full-fledged set of example codes and build scripts. Uncompress the prebuilt library:  
    +   On **Windows**: Simply double-click on the zip-file and select **extract files** from the Windows Explorer menu.  
    +   On **macOS/Linux**: Open a Bash terminal and navigate to the folder containing the compressed library. 
        Use the following command to untar the compressed file,  
        ```  
        ls libparamonte*.tar.gz | xargs -i tar xvzf {}
        ```  
        to extract all libparamonte tar files in the current directory.  

### Building and running ParaMonte simulations on Windows  

+   **Note**: Theoretically, you can use any Fortran compiler on Windows to build and link your applications against the ParaMonte library.  
    A few options currently exist regarding the choice of compilers and environment:  
    +   Use the Intel Fortran compiler along with the Windows **Batch build scripts** of the ParaMonte library examples, as described below.  
    +   Use the **Bash build scripts** that are also supplied with each ParaMonte example on Windows to build and run simulations via 
        the GNU compilers on Windows available in MinGW or Cygwin Linux environments installed on a Windows system.  
    +   Use the GNU Fortran compiler installed on a [Microsoft Windows Subsystem for Linux (WSL)](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). 
        In this case, you will have to download the prebuilt ParaMonte library for the Linux environment as opposed to the Windows OS. 

+   **Install the Microsoft Visual Studio (>2017)**: You will need to have a recent Microsoft Visual Studio (MSVS) installed on your system. 
    The community edition of this software is available free of charge. When installing MSVS, 
    make sure to install all the C++ components of the Visual Studio.  

+   **Install the Intel Parallel Studio**: If you are a student/teacher/faculty/open-source-developer, you can also download and install, free of charge, 
    the most recent **Intel Parallel Studio** on your system which, by default, includes the Intel MPI library. You can follow the instructions given 
    [on this page](https://www.cdslab.org/recipes/programming/intel-parallel-studio-installation-windows/intel-parallel-studio-installation-windows) 
    to install the Intel Parallel Studio on your system.  

+   **Open the right command-line interface to build/run the ParaMonte example**: 
    If the ParaMonte library that you intend to use is built for 64-bit architecture, 
    then make sure you open a 64-bit instance of the command-line interface in either of the two cases below:  
    +   If you have installed Intel Parallel Studio, open an instance of the **command-line interface** that comes 
        with the Intel Parallel Studio from the list of programs in the Windows start menu. This is simply a Windows 
        command prompt that has all the necessary Intel compiler variables and paths predefined in it.  
    +   Otherwise, if you do not have Intel Parallel Studio, open an instance of the **command-line interface** that comes 
        with the Microsoft Visual Studio from the list of programs in the Windows start menu. This is simply a Windows 
        command prompt that has all the necessary compiler variables and paths predefined in it.  

+   **Build and run the ParaMonte example**:  
    Build the example via the Intel Parallel Studio command-line interface,  
    ```  
    build.bat  
    ```  
    The build script will automatically detect whether a parallel simulation has been built. 
    By default, the name of the output executable is `main.exe`. 
    *Note that the build script will both build and run the executable*.  

+   **Run the ParaMonte example executable**:  
    +   For serial simulations, simply type the name of the output executable,  
        ```  
        main.exe
        ```  
    +   For parallel simulations, invoke the MPI launcher `mpiexec`,  
        ```  
        mpiexec -n NUM_PROCESSES main.exe
        ```  
        where `NUM_PROCESSES` represents the number of processes on which the simulation will run. 
        If you are using the Intel MPI library to run your ParaMonte application in parallel, we also recommend using the `-localonly` flag. 
        See [this page](https://www.cdslab.org/paramonte/notes/run/#running-the-manually-generated-executable-on-multiple-processors-on-windows) 
        for usage and utilities of this Intel MPI launcher flag.  

### Building and running ParaMonte simulations on macOS / Linux  

+   **Note**: Theoretically, you can use any Fortran compiler on macOS/Linux to build and link your applications against the ParaMonte library. 
    However, the ParaMonte library example build scripts, as described below, currently only recognize the Intel and GNU Fortran compilers or 
    any compilers with similar executable names.  

+   If you intend to run **serial** ParaMonte simulations, install either,  
    +   **the Intel Fortran compiler (ifort >2018)**, or,  
    +   **the GNU Fortran compiler (gfortran >7.0.0)**,  
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  

+   If you intend to run **MPI parallel** ParaMonte simulations, install either,  
    +   **the Intel Parallel Studio (>2018)** on Linux, or,  
    +   **the GNU Compiler Collection (>7.0.0) and the MPICH (>3.2) library** on Linux,  
    +   **the GNU Compiler Collection (>7.0.0) and the OpenMPI (>4.0) library** on macOS,  
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  
    Note that on **macOS**, only the latter option (the GNU compilers) 
    is available since the Intel MPI library does not support the macOS platform.  

+   If you intend to run Coarray parallel ParaMonte simulations, install either,  
    +   **the Intel Parallel Studio (>2018)**, or,  
    +   **the GNU Compiler Collection (>7.0.0) and OpenCoarrays (>2.8.0)**, 
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  
    Note that on **macOS**, only the latter option (the GNU compilers) 
    is available since the Intel MPI library does not support the macOS platform.  

+   Open a Bash terminal, change directory to the ParaMonte prebuilt library's directory, then build the executable via,  
    ```  
    build.sh  
    ```  
    The build script will automatically detect whether a parallel simulation has to be built. 
    By default, the name of the output executable is `main.exe`. The script will also generate 
    a new Bash script named `run.sh`. To run the generated example executable, type,  
    ```  
    ./run.sh
    ```  
    The script will automatically detect whether the simulation has to be run in parallel or serial. 
    If the simulation is parallel, you can also pass the number of cores on which you want to run the example via,  
    ```  
    ./run.sh --nproc NUM_PROCESSOR
    ```  
    or,  
    ```  
    ./run.sh -n NUM_PROCESSOR
    ```  
    where you will have to replace `NUM_PROCESSOR` with your desired number of processes.  

  
Citing ParaMonte  
================  

The ParaMonte library is an honor-ware, the currency of which is acknowledgment and citations.  
  
If you use ParaMonte or any ideas from the software, please acknowledge it by citing the ParaMonte library's 
main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) 
to access the PDF version of these files free of charge.  

  
License  
=======  

[MIT License](https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md)  

**What does this license mean?**  

Essentially, all we are asking from the users or developers is to  

>   explicitly acknowledge the use of this library or any concepts or parts of it in their education, research, or software (free or commercial).  

This is a free software, so help us keep it freely available to the public by redistributing the library and contributing to it. 
If you have questions or concerns about the license, do not hesitate to contact us (shahmoradi@utexas.edu).  

  
Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, 
        data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, 
        high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, 
        biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, 
        protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+   [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)  
    +   physicist / cosmologist by training,  
    +   currently a UTA Physics member,  
    +   deep philosophical thinker,  
    +   contact: [Fatemeh.Bagheri@uta.edu](mailto:"Fatemeh.Bagheri@uta.edu")  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
---
name: Reproducible bug report
about: Use this template to report a reproducible bug to be fixed. Modify the template
  as deemed necessary.
title: 'Defect: '
labels: ''
assignees: ''

---

**The programming language (e.g., all/C/C++/Fortran/Julia/MATLAB/Python/R/...)**: 

**The compiler/interprerter (e.g., none/GNU/Intel/IBM/NAG/Cray/PGI-NVIDIA/Anaconda/...)**: 

**The operating system (e.g., all/Windows/Linux/macOS/WSL1/WSL2/...)**: 

**The architecture (e.g., all/x86_x64/Intelx64/ARM/POWER9/...)**: 

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots, if applicable**
If applicable, add screenshots to help explain your problem.

**Additional context**
Add any other context about the problem here.
---
name: Enhancement request
about: Suggest implementation of new ideas, samplers, and postprocessing tools in
  this library
title: 'Enhancement: '
labels: ''
assignees: ''

---

**The programming language (e.g., all/C/C++/Fortran/Julia/MATLAB/Python/R/...)**: 

**The operating system (e.g., all/Windows/Linux/macOS/WSL1/WSL2/...)**: 

**The architecture (e.g., all/x86_x64/Intelx64/ARM/POWER9/...)**: 

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. It would be great to have this new algorithm [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered if any**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: General inquiry
about: Generic questions about this library that are not directly related to a reproducible
  defect or enhancement request
title: 'General: '
labels: ''
assignees: ''

---


---
title: 'ParaMonte: A high-performance serial/parallel Monte Carlo simulation library for C, C++, Fortran'
tags:
  - C
  - C++
  - Fortran
  - Monte Carlo
  - Markov Chain Monte Carlo
  - Uncertainty Quantification
  - Metropolis-Hastings
  - adaptive sampling
  - MCMC
  - DRAM
authors:
  - name: Amir Shahmoradi
    orcid: 0000-0002-9720-8937
    affiliation: "1, 2"
  - name: Fatemeh Bagheri
    affiliation: "1"
affiliations:
 - name: Department of Physics, The University of Texas, Arlington, TX
   index: 1
 - name: Data Science Program, The University of Texas, Arlington, TX
   index: 2
date: 28 September 2020
bibliography: paper.bib
---

# Summary

ParaMonte (which stands for Parallel Monte Carlo) is a serial and MPI/Coarray-parallelized library of Monte Carlo routines for sampling mathematical objective functions of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, machine learning, and scientific inference. The ParaMonte library has been developed with the design goal of unifying the **automation**, **accessibility**, **high-performance**, **scalability**, and **reproducibility** of Monte Carlo simulations. The current implementation of the library includes **ParaDRAM**, a **Para**llel **D**elayed-**R**ejection **A**daptive **M**etropolis Markov Chain Monte Carlo sampler, accessible from a wide range of programming languages including C, C++, Fortran, with a unified Application Programming Interface and simulation environment across all supported programming languages. The ParaMonte library is MIT-licensed and is permanently located and maintained at [https://github.com/cdslaborg/paramonte](https://github.com/cdslaborg/paramonte).  

# Statement of need 

Monte Carlo simulation techniques [@metropolis1949monte], in particular, Markov Chain Monte Carlo (MCMC) [@metropolis1953equation] are among the most popular methods of quantifying uncertainty in scientific inference problems. Extensive work has been done over the past decades to develop Monte Carlo simulation programming environments that aim to partially or fully automate the problem of uncertainty quantification via Markov Chain Monte Carlo simulations. Example open-source libraries in C/C++/Fortran include `MCSim` in C [@Bois:2009], `MCMCLib` and `QUESO` [@Alexander:2012] libraries in C++, and `mcmcf90` in Fortran [@Haario:2006]. These packages, however, mostly serve the users of one particular programming language environment. Some can perform only serial simulations while others are inherently parallelized. Furthermore, the majority of the existing packages have significant dependencies on other external libraries. Such dependencies can potentially make the build process of the packages extremely complex and arduous due to software version incompatibilities, a phenomenon that has become known as the *dependency-hell* among software developers.  

The ParaMonte library presented in this work aims to address the aforementioned problems by providing a standalone high-performance serial/parallel Monte Carlo simulation environment with the following principal design goals:  

+   **Full automation** of the library's build process and all Monte Carlo simulations to ensure the highest level of user-friendliness of the library and minimal time investment requirements for building the library as well as running and post-processing the Monte Carlo simulations  
+   **Interoperability** of the core of the library with as many programming languages as currently possible, including C, C++, Fortran, as well as MATLAB and Python via the `ParaMonte::MATLAB` [@2020arXiv201004190S] and `ParaMonte::Python` [@2020arXiv201000724S] libraries  
+   **High-Performance**, meticulously-low-level implementation of the library that guarantees the fastest-possible Monte Carlo simulations **without** compromising the reproducibility of the simulations or the extensive external reporting of the simulation progress and results  
+   **Parallelizability** of all simulations via both MPI and PGAS/Coarray communication paradigms while **requiring zero-parallel-coding efforts from the user**   
+   **Zero external-library dependencies** to ensure hassle-free library builds and Monte Carlo simulation runs  
+   **Fully-deterministic reproducibility** and **automatically-enabled restart functionality** for all ParaMonte simulations (up to 16 digits of decimal precision if requested by the user)  
+   **Comprehensive-reporting and post-processing** of each simulation and its results as well as their efficient compact storage in external files to ensure the reproducibility and comprehensibility of the simulation results at any time in the future  

# The origins of ParaMonte  

The ParaMonte library grew out of the need for a free user-friendly high-performance parallel software for stochastic optimization, sampling, and integration problems in scientific inference and Data Science applications. The project started in 2012 to aid the research goals of the primary author of the package in the field of High Energy Astrophysics and Bioinformatics [@Shahmoradi:2013; @ShahmoradiA:2013; @Shahmoradi:2014; @Shahmoradi:2015; @Shahmoradi:2019; @osborne2020multilevelapj]. It remained in private usage over the years until summer 2020, when version 1.0 of the library was released as an open-source project for public usage and contributions.  

Many contemporary research problems are computationally demanding and require the analysis of vast amounts of high-dimensional data by a community of non-computer-science domain researchers who might be neither familiar with details of software dependencies, build processes, and complexities nor even with the inner-workings of the various stochastic optimization and sampling techniques. As such, since its inception, the ParaMonte library has been built upon the two pillars of user-friendliness and high-performance.  

# The build process  

The ParaMonte library is permanently located on GitHub and is available to view at: [https://github.com/cdslaborg/paramonte](https://github.com/cdslaborg/paramonte). The build process of the library is fully automated. Extensive detailed instructions are also available on the [documentation website of the library](https://www.cdslab.org/paramonte/).  

For the convenience of users, each versioned release of the library's source code also includes prebuilt, ready-to-use, copies of the library for `x64` architecture on Windows, Linux, macOS, in all supported programming languages, including C, C++, Fortran. These prebuilt libraries automatically ship with the language-specific example codes and build scripts that fully automate the process of building and running the examples. Users can either adapt the example's source files and build scripts to their own needs or in more sophisticated scenarios, they can simply link their applications against the supplied prebuilt library.  

Where the prebuilt libraries cannot be used, users can readily call the Unix-Bash and Windows-Batch build-scripts that are provided with the source code of the library to fully automate the build process of the library. These build scripts have been developed to automate the installation of any missing components that may be required for the successful build of the library, including the `cmake` build software, the GNU C/C++/Fortran compilers, as well as the MPI/Coarray libraries. All of these tasks are performed with the explicit permission granted by the user. The ParaMonte build scripts are heavily inspired by the impressive `OpenCoarrays` open-source software [@Fanfarillo:2014] developed and maintained by the [Sourcery Institute](http://www.sourceryinstitute.org/).  

# The ParaMonte samplers  

The current implementation of the ParaMonte library includes the **Para**llel **D**elayed-**R**ejection **A**daptive **M**etropolis Markov Chain Monte Carlo (**`ParaDRAM`**) sampler [@ShahmoradiGS:2020; @ShahmoradiADS:2020; @ShahmoradiASCL:2020; @Kumbhare:2020], and several other samplers whose development is in progress as of writing this manuscript. The ParaDRAM algorithm is a variant of the DRAM algorithm of [@Haario:2006] and can be used in either serial or parallel mode.  

In brief, the ParaDRAM sampler continuously adapts the shape and scale of the proposal distribution throughout the simulation to increase the efficiency of the sampler. This is in contrast to traditional MCMC samplers where the proposal distribution remains fixed throughout the simulation. The ParaDRAM sampler provides a highly customizable MCMC simulation environment whose complete description goes beyond the scope and limits of this manuscript. All of these simulation specifications are, however, expensively explained and discussed on [the documentation website of the ParaMonte library](https://www.cdslab.org/paramonte/). The description of all of these specifications are also automatically provided in the output `*_report.txt` files of every simulation performed by the ParaMonte samplers.  

Several additional more advanced samplers are also currently under development and are scheduled for release in 2021.  

## Parallelism  

Two modes of parallelism are currently implemented for all ParaMonte samplers,  

+   The **Perfect Parallelism** (multi-Chain): In this mode, independent instances of a particular ParaMonte sampler of choice run concurrently. Once all simulations are complete, the sampler compares the output samples from all processors with each other to calculate various statistics and to ensure that no evidence for a lack of convergence to the target density exists in any of the output chains.  

+   The **Fork-Join Parallelism** (single-Chain): In this mode, a single processor is responsible for collecting and dispatching information, generated by all processors, to construct the full sample from the target density function.  

For each parallel simulation in the Fork-Join mode, the ParaMonte samplers automatically compute the speedup gained compared to the serial mode. The speedup for a wide range of the number of processors is also automatically computed and reported in the output `*_report.txt` files that are automatically generated for all simulations. The processor contributions to the construction of the final output sample are also reported along with all visited states in the output `*_chain.*` files. Such information is particularly useful for finding the optimal number of processors for a given problem at hand, by first running a short simulation to predict the optimal number of processors from the sampler's output information, followed by the full production run using the optimal number of processors. For a comprehensive description and algorithmic details see [@ShahmoradiGS:2020; @ShahmoradiADS:2020].  

![An illustration of the contributions of 512 Intel Xeon Phi 7250 processors to a ParaMonte-ParaDRAM simulation parallelized via the Fork-Join paradigm. The predicted best-fit Geometric distribution from the post-processing phase of the ParaDRAM simulation is shown by the black line. The data used in this figure is automatically generated for each parallel simulation performed via any of the ParaMonte samplers. \label{fig:procContribution512}](procContribution512.png)  

As we argue in [@ShahmoradiGS:2020; @ShahmoradiADS:2020] for the particular case of MCMC simulations, the distribution of processor contributions to the construction of a final sample from an objective function in the Fork-Join parallelism paradigm follows a Geometric curve. \autoref{fig:procContribution512} depicts the processor contributions to an example ParaDRAM simulation of a variant of Himmelblau's function. Superimposed on the distribution of processor contributions is the Geometric fit to the data. \autoref{fig:PredictedSpeedupActualSpeedup} illustrates an example strong-scaling behavior of the sampler and the predicted speedup by the sampler for a range of processor counts. 

The exact strong scaling behavior of Fork-Join parallel simulations is highly dependent on the sampling efficiency of the problem. For a fixed number of processors, the parallel efficiency of the samplers monotonically increases toward the maximum theoretical speedup with decreasing sampling efficiency. This theoretical maximum parallel speedup and the corresponding predicted optimal number of processors are also computed for each parallel simulation and reported in the final output `*_report.txt` files.

![A comparison of the actual strong scaling behavior of an example ParaMonte-ParaDRAM simulation from 1 to 1088 processors with the strong-scaling behavior predicted during the post-processing phases of the corresponding individual ParaDRAM simulations. The data used in this figure is automatically generated for each parallel simulation performed via any of the ParaMonte samplers. The black dashed line represents the perfect parallelism.\label{fig:PredictedSpeedupActualSpeedup}](PredictedSpeedupActualSpeedup.png)  

## Efficient compact storage of the output samples  

Efficient continuous external storage of the output of Monte Carlo simulations is essential for both post-processing of the results and the restart functionality of the simulations, should any interruptions happen at runtime. However, as the number of dimensions or the complexity of the target density increases, such external storage of the output can easily become a challenge and a bottleneck in the speed of an otherwise high-performance sampler. Given the currently-available computational technologies, input/ouput (IO) to external hard-drives can be 2-3 orders of magnitude slower than the Random Access Memory (RAM) storage.  

To alleviate the effects of such external-IO speed bottlenecks, the ParaMonte samplers have been devised to carefully store the resulting samples in a small *compact*, yet ASCII human-readable format in external output files. This **compact** storage, as opposed to the **verbose** (or in the case of the ParaDRAM sampler, the **Markov-chain**) storage format, leads to significant speedup of the simulations while requiring 4-100 times less external memory to store the output samples in the external output files. The exact amount of reduction in the external memory usage depends on the runtime sampling efficiency of the samplers. Additionally, the format of output files can be set by the user to `binary`, further reducing the memory foot-print of the simulations while increasing the simulation speed. The implementation details of this compact storage format are extensively discussed in [@ShahmoradiGS:2020; @ShahmoradiADS:2020].  

## Restart functionality  

Each ParaMonte sampler is automatically capable of restarting an existing interrupted simulation, whether in serial or parallel. All that is required is to rerun the interrupted simulation with the same prefix for the simulation output file names. The ParaMonte samplers automatically detect the presence of an incomplete simulation in the output files and restart the simulation from where it was left off.  

Furthermore, if the user sets the seed of the random number generator of the sampler before running the simulation, *the ParaMonte samplers are capable of regenerating the same output sample that would have been produced if the simulation had not been interrupted in the first place*. Such **fully-deterministic reproducibility into-the-future** is guaranteed with 16 digits of decimal precision for the results of any ParaMonte simulation, whether serial or in parallel. To our knowledge, this is a unique feature of the ParaMonte library that does not appear to exist in any of the contemporary libraries for Monte Carlo simulations.  

# The ParaMonte Application Programming Interface  

Special care has been taken to develop highly-similar (if not the same) Application Programming Interface (API) to the ParaMonte library samplers across all supported programming languages. Ensuring backward-compatibility has been also one of the high priorities in the development of the library and will remain a priority for any future development and expansion of the project.  

The project currently uses semantic versioning to readily inform the users about the backward-compatibility of each release. This semantic versioning, however, only applies to the ParaMonte samplers' API. Conversely, the API to the foundational procedures of the ParaMonte kernel tend to be dynamic and subject to changes without incrementing the library's major version. This is line with the primary purpose of the library as a Monte Carlo optimization and sampling toolbox.  

Nevertheless, the kernel documentation is well developed and accessible to all, including the end-users, should they decide to use any of the basic routines in the library beyond the samplers. The kernel API documentation is permanently available on GitHub at its dedicated repository: \url{https://cdslaborg.github.io/paramonte-kernel-doc/html/} and is accessible from [the ParaMonte's main documentation portal](https://www.cdslab.org/paramonte/).

# Tests and code coverage  

The entire kernel routines of the ParaMonte library are currently tested with over 866 individual tests. Out of these, 510 test units currently cover and verify the functionalities of the basic building-block routines of the library and 156 tests verify the semantics and different aspects of the ParaMonte samplers built on top of the basic routines.  

The code coverage report analyses of the ParaMonte kernel library can be generated for a wide range of input build configurations. For example, corresponding to each programming language (C, C++, Fortran, MATLAB, Python, R, …), a separate code coverage report for the kernel routines can be generated. Similarly, different library builds (debug, testing, release), library types (static, dynamic), and memory allocation schemes (stack or heap) can lead to different code coverage reports. This is primarily because the kernel routines heavily utilize compiler preprocessor directives to accommodate the support for multiple programming languages, Operating Systems and architectures, as well as multiple parallelism paradigms.  

Among all build configurations, however, the single most important specification that yields the largest differences in code coverage reports is the parallelism paradigm (serial, MPI, and Coarray parallelism). The code coverage reports for the serial, MPI, and Coarray modes are remarkably different since each parallelism paradigm activates a different set of codes in the ParaMonte kernel routines.  

Thus, with each tagged release of the library, the code coverage reports for the three parallelism paradigms are also generated and stored in a separate repository dedicated to the code coverage analyses of the ParaMonte library. This repository is permanently located on GitHub at: \url{https://www.cdslab.org/paramonte/notes/codecov/} and the reports are automatically and conveniently accessible from [the ParaMonte documentation website](https://www.cdslab.org/paramonte/notes/codecov/).

The current set of 866 unit tests in version 1.5.0 release of the ParaMonte library collectively cover 47565 out of 48384 lines of code ($\sim 98\%$ line coverage) and 4467 out of 4476 procedures ($\sim 100\%$ function coverage) in the library's kernel for the three serial, MPI, and Coarray parallelism paradigms.  

# Documentation and Repository  

Extensive documentation and examples in C, C++, Fortran (as well as other programming languages) are available on the documentation website of the library at: [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/). The ParaMonte library is MIT-licensed and is permanently located and maintained at [https://github.com/cdslaborg/paramonte](https://github.com/cdslaborg/paramonte).  

A high number of comments in a codebase are frequently indicative of an organized well-documented project and considered a sign of a helpful and disciplined development team. As of January 2021, comments approximately make up $31\%$ of the entire ParaMonte codebase. Currently, the Open Hub indexer service ranks the ParaMonte project's level of documentation as "impressive" and among the top $10\%$ of all projects in the same category on Open Hub. 


# Funding and Acknowledgements  

The development cost of ParaMonte is estimated over 1 million dollars by [Open Hub](https://www.openhub.net/p/paramonte/estimated_cost). Considering the presence of  multiple high-level dynamic languages such as MATLAB, Python, and R, which together comprise approximately $40\%$ of the repository, an estimate of $500-600$ thousand dollars seems more realist for the development costs of the library's kernel routines. The current state of the library has resulted from multiple years of continuous development predominantly by the core developers of the library as graduate student, postdoc, and faculty in their free time.  The project was also partially supported in its final stages of development by the Peter O'Donnell, Jr. Postdoctoral Fellowship awarded by the Oden Institute for Computational Engineering and Sciences at The University of Texas at Austin to the primary developer of the library, Amir Shahmoradi. We thank the Texas Advanced Computing Center for providing the supercomputer time for testing and development of this library.  

# References  
