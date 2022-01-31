Download, Compile and Link OpenBLAS Library with R
================
  
  
  
**Package website**: **<https://prdm0.github.io/ropenblas>**  
  
**Summary**

-   [General aspects](#general-aspects)
-   [Advantages of using ropenblas
    package](#advantages-of-using-ropenblas-package)
-   [Dependencies](#dependencies)
-   [Installation](#installation)
-   [Use](#use)
    -   [‘ropenblas’ function](#ropenblas-function)
    -   [‘last\_version\_r’ function](#last_version_r-function)
    -   [‘last\_version\_openblas’
        function](#last_version_openblas-function)
    -   [‘rcompiler’ function](#rcompiler-function)
    -   [‘link\_again’ function](#link_again-function)
    -   [‘rnews’ function](#rnews-function)

**ropenblas package**

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02769/status.svg)](https://doi.org/10.21105/joss.02769)
[![Travis build
status](https://travis-ci.com/prdm0/ropenblas.svg?branch=master)](https://travis-ci.com/prdm0/ropenblas)
[![last](https://www.r-pkg.org/badges/last-release/ropenblas)](https://CRAN.R-project.org/package=ropenblas)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/ropenblas)](https://CRAN.R-project.org/package=ropenblas)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/ropenblas)](https://CRAN.R-project.org/package=ropenblas)
[![month](https://cranlogs.r-pkg.org/badges/ropenblas)](https://CRAN.R-project.org/package=ropenblas)

<!-- ```{r, fig.align='center', echo=FALSE, out.width=219} -->
<!-- knitr::include_graphics("https://raw.githubusercontent.com/prdm0/ropenblas/master/logo.png") -->
<!-- ``` -->

<img src="logo.png" height="230" width="200" align="right" />

# General aspects

The [**ropenblas**](https://prdm0.github.io/ropenblas/) is a package
designed to facilitate the linking of the library
[**OpenBLAS**](https://www.openblas.net/) with the language
[**R**](https://www.r-project.org/). The package, which works only for
Linux systems, will automatically download the latest source code from
the [**OpenBLAS**](https://www.openblas.net/) library and compile the
code. The package will automatically bind the language
[**R**](https://www.r-project.org/), through the `ropenblas()` function,
to use the [**OpenBLAS**](https://www.openblas.net/) library. Everything
will be done automatically regardless of the Linux distribution you are
using.

You can also specify older versions of the
[**OpenBLAS**](https://www.openblas.net/) library. Automatically, if no
version is specified, the
[**ropenblas**](https://prdm0.github.io/ropenblas/) package will
consider the latest version of the library
[**OpenBLAS**](https://www.openblas.net/).

Considering using the [**OpenBLAS**](https://www.openblas.net/) library
rather than the [**BLAS**](http://www.netlib.org/blas/) may bring extra
optimizations for your code and improved computational performance for
your simulations, since [**OpenBLAS**](https://www.openblas.net/) is an
optimized implementation of the library
[**BLAS**](http://www.netlib.org/blas/).

Some of the reasons why it is convenient to link
[**R**](https://www.r-project.org/) language to the use of
[**BLAS**](http://www.netlib.org/blas/) optimized alternatives can be
found [**here**](https://csantill.github.io/RPerformanceWBLAS/). Several
other [**benchmarks**](https://en.wikipedia.org/wiki/Benchmarking) that
point to improved computing performance by considering the library
[**OpenBLAS**](https://www.openblas.net/) can be found on the internet.

The ropenblas package, by `rcompiler()` function is also useful if you
want to install different versions of the
[**R**](https://www.r-project.org/) language. The different versions,
specified by the user of the [**R**](https://www.r-project.org/)
language, will be compiled and will also be linked to the
[**OpenBLAS**](https://www.openblas.net/) library. If you want to switch
between compiled versions of the [**R**](https://www.r-project.org/)
language, no compilation is needed anymore. This allows you to avoid
having to get your hands dirty with tedious operating system settings,
regardless of your GNU/Linux distribution. Another great use of the
`rcompiler()` function is that you will not be dependent on updating
your GNU/Linux distribution repositories and you can always have the
latest version of the [**R**](https://www.r-project.org/) language.

The use of the ropenblas package will return warnings that help you
proceed with the use of the functions. If your internet is not working
or if any dependency on the operating system is not present, the package
will let you know.

# Advantages of using ropenblas package

Some advantages of using the
[**ropenblas**](https://prdm0.github.io/ropenblas/) library:

-   Everything is done within the [**R**](https://www.r-project.org/)
    language;

-   The procedure will be the same for any Linux distribution;

-   The [**OpenBLAS**](https://www.openblas.net/) library will be
    compiled and you will choose which build version to bind to
    [**R**](https://www.r-project.org/), regardless of your Linux
    distribution;

-   If your GNU/Linux distribution does not have updated versions of
    [**OpenBLAS**](https://www.openblas.net/), it matters little. The
    ropenblas package fetches the latest stable release of the
    [**OpenBLAS**](https://www.openblas.net/) library development
    account on GitHub;

-   You do not need to know Linux well. In some distributions, it may
    not be so simple for a less experienced user to compile and link the
    library to the [**OpenBLAS**](https://www.openblas.net/) library
    with the [**R**](https://www.r-project.org/) language;

-   It is much easier to direct a person to link
    [**OpenBLAS**](https://www.openblas.net/) with
    [**R**](https://www.r-project.org/) saying “run `ropenblas()` within
    [**R**](https://www.r-project.org/)” than asking that person to
    verify that an unoptimized version of
    [**BLAS**](http://www.netlib.org/blas/) installed on the system.
    Then you have to guide the removal of the unoptimized version of
    [**BLAS**](http://www.netlib.org/blas/) and guide it to the
    installation of the library
    [**OpenBLAS**](https://www.openblas.net/) through the most diverse
    procedures depending on the GNU/Linux distribution used;

<!-- - As stated earlier, the procedure works for any Linux and this includes Android. If your Android is capable of running privileged commands (ROOT) and if you have [**R**](https://www.r-project.org/) installed via Termux with the required dependencies, you can compile and link [**OpenBLAS**](https://www.openblas.net/) with [**R**](https://www.r-project.org/) using ropenblas; -->

-   With the `rcompiler()` function you can build any version of
    [**R**](https://www.r-project.org/) into your computer architecture,
    which includes the most stable version of the language;

-   What is the latest stable version of
    [**R**](https://www.r-project.org/)? Are you too lazy to go to the
    [**R**](https://www.r-project.org/) site? Run
    `ropenblas::last_version_r(major = NULL)`.

# Dependencies

You must install the following dependencies on your operating system
(Linux):

1 - **GNU Make**: GNU Make utility to maintain groups of programs; <br/>

2 - **GNU GCC Compiler (C and Fortran)**: The GNU Compiler Collection -
C and Fortran frontends.

Do not worry that you will be notified if any of these dependencies are
not installed.

# Installation

Installing the [**ropenblas**](https://prdm0.github.io/ropenblas/)
library is easy and will require you to have installed the **devtools**
package. This will allow you to install the
[**ropenblas**](https://prdm0.github.io/ropenblas/) package directly
from GitHub. To install, after installing the **devtools** package, do:

    remotes::install_github(repo = "prdm0/ropenblas", force = TRUE)

or

    install.packages("ropenblas")

<!-- **Note**: If you want to access the latest features of the [**ropenblas**](https://prdm0.github.io/ropenblas/) package, install it using the first procedure. -->

# Use

The [**ropenblas**](https://prdm0.github.io/ropenblas/) package
currently provides five functions: `ropenblas()`, `rcompiler()`,
`last_version_openblas()`, `last_version_r()`, `link_again()` and
`rnews()`. First of all, do:

    library(ropenblas)

## ‘ropenblas’ function

Installing, compiling, and linking the
[**OpenBLAS**](https://www.openblas.net/) version **0.3.13** library to
the [**R**](https://www.r-project.org/) language:

<img src="https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/ropenblas.svg" style="display: block; margin: auto;" />

<!-- **Notes**:  -->
<!--    - You do not have to in every section of [**R**](https://www.r-project.org/) make use of the `ropenblas()` function. Once the function is used, [**R**](https://www.r-project.org/) will always consider using the [**OpenBLAS**](https://www.openblas.net/) library in future sections. -->
<!--    - [**OpenBLAS**](https://www.openblas.net/) versions tested: 0.3.0, 0.3.1, 0.3.2, 0.3.3, 0.3.4, 0.3.5, 0.3.6, 0.3.7, 0.3.8, 0.3.9, 0.3.10, 0.3.11, 0.3.12 and 0.3.13. These are the values that will be passed to `x` in `ropenblas(x)`;  -->
<!--    - If `x = NULL`, the latest stable version of the [**OpenBLAS**](https://www.openblas.net/) library will be compiled and linked to [**R**](https://www.r-project.org/). -->
<!-- **Details** -->
<!-- Your linux operating system may already be configured to use the [**OpenBLAS**](https://www.openblas.net/) library. Therefore, most likely [**R**](https://www.r-project.org/) will already be linked to this library. To find out if the [**R**](https://www.r-project.org/) language is using the [**OpenBLAS**](https://www.openblas.net/) library, at [**R**](https://www.r-project.org/), do: -->
<!-- ``` -->
<!-- extSoftVersion()["BLAS"] -->
<!-- ``` -->
<!-- If [**R**](https://www.r-project.org/) is using the [**OpenBLAS**](https://www.openblas.net/) library, something like `/any_directory/libopenblas.so` should be returned. Therefore, there should be the name **openblas** in the **s**hared **o**bject returned (file extension **.so**). -->
<!-- If the `ropenblas()` function can identify that the [**R**](https://www.r-project.org/) language is using the version of [**OpenBLAS**](https://www.openblas.net/) you wish to configure, a warning message will be returned asking if you really would like to proceed with the configuration again. -->
<!-- The `ropenblas()` function will download the desired version of the library [**OpenBLAS**](https://www.openblas.net/), compile and install the library in the `/opt` directory of your operational system. If the directory does not exist, it will be created so that the installation can be completed. Subsequently, files from the version of [**BLAS**](http://www.netlib.org/blas/) used in [**R**](https://www.r-project.org/) will be symbolically linked to the shared object files of the library version [**OpenBLAS**](https://www.openblas.net/) compiled and installed in `/opt`. -->
<!-- You must be the operating system administrator to use this library. Therefore, do not attempt to use it without telling your system administrator. If you have the ROOT password, you will be responsible for everything you do on your operating system. -->
<!-- You will not necessarily have to run `ropenblas()` on every section of [**R**](https://www.r-project.org/). Almost always it will not be necessary. However, it may be that the [**R**](https://www.r-project.org/) is updated by the operating system (GNU/Linux). Thus, it may be that in this update the [**R**](https://www.r-project.org/) unlink with the [**OpenBLAS**](https://www.openblas.net/) library. Therefore, from time to time check using the command `extSoftVersion()["BLAS"]` if the link with [**OpenBLAS**](https://www.openblas.net/) is correct, otherwise run the command `ropenblas()` again. -->

## ‘last\_version\_r’ function

Given the higher version, the function will return the latest stable
version of the [**R**](https://www.r-project.org/) language. See the
following example:

![](https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/last_version_r.svg)<!-- -->

<!-- **Note**: If `major = NULL`, the function will consider the major release number. -->

## ‘last\_version\_openblas’ function

The `last_version_openblas()` function automatically searches
[**OpenBLAS**](https://www.openblas.net/) library versions in the
official [**GitHub**](https://github.com/xianyi/OpenBLAS) project.

![](https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/last_version_openblas.svg)<!-- -->

## ‘rcompiler’ function

This function is responsible for compiling a version of the
[**R**](https://www.r-project.org/) language. The `x` argument is the
version of [**R**](https://www.r-project.org/) that you want to compile.
For example, `x = "4.0.4"` will compile and link **R-4.0.4** version as
the major version on your system. By default (`x = NULL`) will be
compiled the latest stable version of the
[**R**](https://www.r-project.org/).

For example, to compile the latest stable version of the
[**R**](https://www.r-project.org/) language, do:

![](https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/rcompiler.svg)<!-- -->

Regardless of your GNU/Linux distribution and what version of
[**R**](https://www.r-project.org/) is in your repositories, you can
have the latest stable version of the
[**R**](https://www.r-project.org/) language compiled into your computer
architecture.

You can use the `rcompiler()` function to compile different versions of
[**R**](https://www.r-project.org/). For example, running
`rcompiler(x = "3.6.3")` and `rcompiler()` will install versions 3.6.3
and 4.0.0 on its GNU/Linux distribution, respectively. If you are in
version 4.0.0 of [**R**](https://www.r-project.org/) and run the code
`rcompiler(x = "3.6.3")` again, the function will identify the existence
of version 3.6.3 in the system and give you the option to use the
binaries that were built in a previous compilation. This avoids
unnecessarys compilations.

<!-- In addition to the `x` argument, the` rcompiler()` function has two other arguments that will allow you to change and pass new compilation flags. Are they: -->
<!-- 1. `with_blas`: This argument sets the `--with-blas` flag in the R language compilation process and must be passed as a string. Details on the use of this flag can be found [**here**](https://cran.r-project.org/doc/manuals/r-devel/R-admin.html). If `with_blas = NULL` (default), then it will be considered: -->
<!--    ``` -->
<!--    ./configure --prefix=/opt/R/version_r --enable-memory-profiling --enable-R-shlib --enable-threads=posix -->
<!--    --with-blas="-L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include -lpthread -lm" -->
<!--    ``` -->
<!--    Most likely, you will have little reason to change this argument. Unless you know what you're doing, consider `with_blas = NULL`. Do not change the installation directory, that is, always consider `--prefix = /opt/R/version_r`, where` version_r` is a valid version of [**R**](https://www.r-project.org/). For a list of valid versions of [**R**](https://www.r-project.org/), run the `last_version_r()`. -->
<!-- Installing [**R**](https://www.r-project.org/) in the `/opt/R/version_r` directory is important because some functions in the package require this. Both the [**R**](https://www.r-project.org/) language and the [**OpenBLAS**](https://www.openblas.net/) library will be installed in the `/opt` directory. If this directory does not exist in your GNU/Linux distribution, it will be created; -->
<!-- 2. `complementary_flags`: String (`complementary_flags = NULL` by default) for adding complementary flags in the [**R**](https://www.r-project.org/) language compilation process. Passing a string to `complementary_flags` will compile it in the form: -->
<!--    ``` -->
<!--     ./configure --with-blas="..." complementary_flags -->
<!--    ``` -->

## ‘link\_again’ function

The `link_again` function links again the
[**OpenBLAS**](https://www.openblas.net/) library with the
[**R**](https://www.r-project.org/) language, being useful to correct
problems of untying the [**OpenBLAS**](https://www.openblas.net/)
library that is common when the operating system is updated.

The function `link_again` be able to link again the
[**R**](https://www.r-project.org/) language with the
[**OpenBLAS**](https://www.openblas.net/) library. Thus, link\_again
will only make the relinkagem when in some previous section of
[**R**](https://www.r-project.org/) the ropenblas function has been used
for the initial binding of the [**R**](https://www.r-project.org/)
language with the [**OpenBLAS**](https://www.openblas.net/) library.

For example, to relink the [**OpenBLAS**](https://www.openblas.net/)
library with the [**R**](https://www.r-project.org/) language, do:

![](https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/link_again.svg)<!-- -->

If `restart_r = TRUE` (default), a new section of
[**R**](https://www.r-project.org/) is started after linking the
[**OpenBLAS**](https://www.openblas.net/) library.

In situations where there was a disconnection due to an update of the
operating system, the `ropenblas` function can be used to relink the
[**OpenBLAS**](https://www.openblas.net/) library with the
[**R**](https://www.r-project.org/) language, however, it will be
necessary to compile the [**OpenBLAS**](https://www.openblas.net/)
library again. If you are interested in recompiling the
[**OpenBLAS**](https://www.openblas.net/) library and linking with
[**R**](https://www.r-project.org/), use the `ropenblas` function. If
the interest is to take advantage of a previous compilation of the
[**OpenBLAS**](https://www.openblas.net/) library, the function
`link_again` may be useful.

## ‘rnews’ function

Returns the contents of the
[**NEWS.html**](https://cran.r-project.org/doc/manuals/r-release/NEWS.html)
file in the standard browser installed on the operating system. The
[**NEWS.html**](https://cran.r-project.org/doc/manuals/r-release/NEWS.html)
file contains the main changes from the recently released versions of
the [**R**](https://www.r-project.org/) language. The goal is to
facilitate the query by invoking it directly from the
[**R**](https://www.r-project.org/) command prompt. The rnews function
is analogous to the news function of the utils package. However, using
the news command in a terminal style bash shell is possible to receive a
message like:

    > news()
    starting httpd help server ... done
    Error in browseURL(url) : 'browser' must be a non-empty character string

If `pdf = FALSE` (default), the
[**NEWS.html**](https://cran.r-project.org/doc/manuals/r-release/NEWS.html)
file will open in the browser, otherwise
[**NEWS.pdf**](https://cran.r-project.org/doc/manuals/r-release/NEWS.pdf)
will be opened. If `dev = FALSE` (default), it will not show changes
made to the language development version. To see changes in the
development version, do `dev = TRUE`.
# ropenblas (development version)

* Code review.

# ropenblas 0.2.10

* Code review;

* The **rvest** package is now a dependency.

# ropenblas 0.2.9

* Code review;

* Bug fix: adding missing symbolic links in the else statement of the `compiler_r` function. `compiler_r` is an internal function used by the `rcompiler` function.

# ropenblas 0.2.8

* Code review;

* Bug fix in submission **0.2.7**. Function not imported in the **NAMESPACE** file.

# ropenblas 0.2.7

* Code review;

* Bugs have been fixed;

* Add `rnews` function. This function will allow the user to view the changes in the recent release versions of the R language. It is also possible to check the changes in the development version of R. It is an alternative to the `news` function.

# ropenblas 0.2.6 

* Code review;

* Bugs have been fixed.

# ropenblas 0.2.5

* Code review;

* Bugs have been fixed;

* It is no longer necessary to enter the ROOT password several times;

* The arguments `with_blas` and `complementary_flags` have been added which allows the modification of some compilation flags;

* The `version_openblas` argument of the `rcompiler()` function has been removed. The `rcompiler()` function now keeps the pre-configured version of OpenBLAS or links the newest version of the OpenBLAS library to R;

* The fs package was removed as a dependency. 

# ropenblas 0.2.4

* Code review;

* Identifying bugs;

* Implementing the `link_again()` function. This function will allow you to link the R language again with the OpenBLAS library without having to recompile the OpenBLAS library. The `link_again()` function will be useful in situations where the operating system decouples the library from the OpenBLAS library. This is common in situations where the BLAS library binary file is updated.

# ropenblas 0.2.3

* Code review;

* Improve code efficiency.

# ropenblas 0.2.2 

* The `last_version_openblas()` function will be added.

# ropenblas 0.2.1

* Code review;

* Highlighting some messages;

* In the tests performed, apparently an error will not break the installed version of R.


# ropenblas 0.2.0

* Function `last_version_r(major = 3L)` implemented.  Given the higher version, the function will return the latest stable version of the R language. Major release number of R language (eg. `1L`, `2L`, `3L`, ...);

* The `rcompiler()` function has been exported. This function is responsible for compiling a version of the R language;

* It now depends on the version of R >= 3.1.0;

* Compiled versions of R will be kept. When attempting to compile a previously compiled version, you will be given the option to only switch between versions of R without compiling;

* An alert will be issued when attempting to compile and link a version of R earlier than version 3.6.1. The user will have to respond three times that they are aware that they will not be able to use the ropenblas package to return to a newer version of R.

# ropenblas 0.1.0

* Maintenance release;

* General documentation improvements;

* Check if there is no internet connection;

* The `ropenblas()` function documentation describes the details of changes that are made to the system;

* The `ropenblas()` function is capable of suggesting a stable and newer version of the OpenBLAS library. The user can decide whether to use the latest version or the one of his choice. By default if no argument is passed to the `ropenblas()` function, the latest version of the OpenBLAS library will be considered.Contributing to the development of the ropenbals package
================
**Package website**: **<https://prdm0.github.io/ropenblas>**  
  
**Summary**

-   [1 General information](#general-information)
-   [2 Main contributions](#main-contributions)

# 1 General information

The **ropenblas** package introduces a simple way, through its
functions, to link the library [**OpenBLAS**](https://www.openblas.net/)
with the language [**R**](https:%20//www.r-project.org/). In addition,
with ropenblas it is possible to link different versions of R, making it
possible to easily switch between the different versions of R compiled.

The package currently works only on GNU/Linux systems and is independent
of the official repositories of the user’s GNU/Linux distribution. This
allows the R user on GNU/Linux systems who is the administrator of the
operating system - OS to link the latest versions of R and OpenBLAS as
well as switch between versions if they wish.

# 2 Main contributions

The main contributions to the **ropenblas** package are listed below:

1.  Revision of the code to detect possible inconsistent points that may
    cause the functions provided by the package to fail in some GNU /
    Linux distribution;
2.  Test the package on different GNU / Linux distributions, reporting
    possible inconsistencies at
    <https://github.com/prdm0/ropenblas/issues>;
3.  Improved English for package messages and documentation;
4.  Improved writing of job documentation;
5.  Allow Windows system users to take advantage of the `ropenblas()`
    function. One suggestion would be to link precompiled binaries from
    the OpenBLAS library and distribute them in a `inst/` directory.
## Test environments
* local R installation, R 4.0.4
* ubuntu 16.04 (on travis-ci), R 4.0.4
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
---
# Example from https://joss.readthedocs.io/en/latest/submitting.html
title: 'ropenblas: Download, Compile and Link OpenBLAS Library with R'
tags:
  - R
  - Compiling R and OpenBLAS
  - Link OpenBLAS
  - Switch between versions of R
  - Fast algebraic computing
authors:
  - name: Pedro Rafael Diniz Marinho
    orcid: 0000-0003-1591-8300
    affiliation: 1
affiliations:
 - name: Department of Statistics, Federal University of Paraíba, João Pessoa, Paraíba - PB, Brazil
   index: 1
citation_author: Marinho
date: "2021-03-30"
year: "2021"
bibliography: paper.bib
output: rticles::joss_article
header-includes:
  - \usepackage{float}
csl: apa.csl
journal: JOSS
---

# Summary



<!-- The `ropenblas` package aims to facilitate the day-to-day life of R programmers who want more performance on GNU/Linux systems, without removing the possibility that specific configurations are made, if they deem convenient, that is, more technical users will be able to pass other flags that will be considered in the compilation process. Through the package's `ropenblas()` and `rcompiler()` functions, the library user will be able to compile and link the R language in his GNU/Linux distribution with the OpenBLAS library, all within R and in a very simple fashion. All functions work without being influenced by the GNU/Linux distribution and are independent of their repositories, that is, it does not matter which GNU/Linux distribution is being used. Linking the OpenBLAS library to R will bring better computational performance to the language in the most diverse algebraic operations commonly used in areas such as statistics, data science, and machine learning. -->

The `ropenblas` package provides improved performance on GNU/Linux systems and allows for passing additional compilation flags for more technical users. Through the `ropenblas()` and `rcompiler()` functions, the user will be able to compile and link the GNU/Linux distribution R language with the OpenBLAS library, all within R and in a very simple fashion. It works for all GNU/Linux distributions. Linking the OpenBLAS library to R brings better computational performance to the language in the most diverse algebraic operations commonly used in areas such as statistics, data science, and machine learning.

<!-- # Statement of Need -->

<!-- The `ropenblas` package aims to allow algebraic computing of R to be performed using the OpenBLAS library, that is, it allows easy linking of R to the OpenBLAS library on GNU/Linux systems without depending on the distribution repositories. This will allow several researchers in the areas of statistics, data science, and machine learning to take advantage of more efficient performance in algebraic calculations, for example, multiplication, factorization, and matrix inversion. The `ropenblas` library, version 0.2.9 will also allow the R programmer to have, in his GNU/Linux distribution, several compiled versions of the R language giving the possibility to easily switch between these versions, this being an Open Source functionality that is only possible in some commercial IDEs of R. All this is done within the R language, minimizing the chance of less experienced users to break their operating system by running several instructions that are not perfectly understood. -->

<!-- The fact that the `ropenblas` package does not depend on the repositories of the GNU/Linux distribution will allow that in more stable distributions the R programmer will have at his disposal the most recent version of the OpenBLAS libraries and the R language. Everything is done safely, since the `ropenblas` package uses the stable versions of the official development repositories of the OpenBLAS library and the R programming language, respectively. Until the present version, the package has more than 11000 downloads, having an average of more than 1000 downloads in the month before the date of this submission, based on the official R language repositories. -->

# Introduction

<!-- The term "computational efficiency" is very common for those who program statistical methods, in which a large part of them involve algebraic operations that are often reproduced in computationally intensive simulations, such as Monte-Carlo simulations - MC and resampling methods, as is the case with bootstrap resampling. Statistics is just one example within so many other areas that need performance and uses the R language. -->

Computational efficiency is important for those who program statistical methods since they often involve algebraic operations reproduced in computationally intensive simulations, such as Monte-Carlo simulations and resampling methods, as is the case with bootstrap resampling. Statistics is just one example within other areas that need good performance and use the R language.

In addition to the adoption of good programming practices and the maximum, efficient and adequate use of available computational resources, such as code parallelization, through multicore parallelism procedures allowed by most current processors and operating systems, small adjustments and linkage of libraries can provide useful benefits.

<!-- The `ropenblas` package aims to provide useful and simple experiences to R [@R] programmers who develop their activities on GNU/Linux operating systems. These experiences consist of being able to link any version of the OpenBLAS [@openblas] library to the R language, as well as allowing the programmer to install and link various versions of R and make them available on their operating system as well as switch between these versions as they see fit. -->

The `ropenblas` package aims to be easy to use for R [@R] programmers who work on GNU/Linux operating systems. For example, a user can link any version of the OpenBLAS [@openblas] library to the R language and install and link various versions of R to make them available on their operating system as well as switch between these versions.

Linking the R language to the OpenBLAS library can bring several benefits to algebraic computing in R. OpenBLAS is an Open-Source implementation of the Basic Linear Algebra Subprograms - BLAS library that is often the first library option for algebraic computing to be linked in the installation of R on many GNU/Linux distributions. The OpenBLAS library is available at <https://github.com/xianyi/OpenBLAS> and adds optimized implementations of linear algebra kernels that can run optimized on various processor architectures. OpenBLAS is based on the GotoBLAS2 project code in version 1.13 [@gotoblas2], code available under the terms of the BSD license.

<!-- The functions of the `ropenblas` package can help the R language user to make this link without leaving the R promprt, as well as allowing the user to choose the version of OpenBLAS to be considered, the last stable version being considered by default . Everything is accomplished by functions without many arguments and without major complications to be used. These functions can be very comforting for users of R on GNU/Linux systems who do not feel safe to run several codes in a terminal that runs Shell Script codes for the language configuration. -->

<!-- It is very common to see users of R in distributions with repositories not prone to immediate updates of programs using several tutorials found in communities on the internet suggesting several lines of code and steps that can be potentially dangerous or easy to be misunderstood by many users of R language. Being able to compile, enable resources if necessary and switch between versions of R using simple functions to be used and without leaving the command pompt of the R language is attractive. -->

<!-- # General package information -->

The `ropenblas` is a package designed to facilitate the linking of the library OpenBLAS with the language R. The package, which works only for Linux systems, will automatically download the latest source code from the OpenBLAS library and compile the code. The package will automatically bind the language R, through the `ropenblas()` function, to use the OpenBLAS library. Everything will be done automatically regardless of the Linux distribution you are using. Enumerating some advantages of the package:

1.  Everything is done within the R language;
2.  The procedure (use of functions) will be the same for any Linux distribution;
3.  The OpenBLAS library will be compiled and you will choose which build version to bind to R, regardless of your Linux distribution;
4.  The package allows you to install R $\geq 3.1.0$, also allowing you to install one more version, in addition to allowing you to easily switch between those versions;
5.  The linked versions of R will continue to be recognized by their Integrated Development Environment - IDE and nothing will have to be adjusted in your GNU/Linux distribution after using any function of the package;
6.  Unnecessary builds will be avoided. Therefore, if you need to switch between compiled versions of the R language, the use of binaries compiled compiled previously will be suggested;
7.  If any errors occur, the functions of the package will not damage the previous installation of the language;
8.  If something better can be done or if a newer version of what you want to install (R or OpenBLAS) exists, the functions will automatically suggest that you should consider installing newer versions.

The `ropenblas` package is already available on the Comprehensive R Archive Network - CRAN, currently in version 0.2.9, and the project is maintained on GitHub at <https://github.com/prdm0/ropenblas> where contributors can find other details of the code, information, as well as being able to contribute with the development of the project. On the website, it is also possible to read the `NEWS.md` file with details of the versions and the focus of the current development. The site is deposited at <https://prdm0.github.io/ropenblas/>. Suggestions for improvements and bug reports can be sent via the link <https://github.com/prdm0/ropenblas/issues>. You can find out how to contribute to the package by accessing the `CONTRIBUTING.md` file at <https://github.com/prdm0/ropenblas/blob/master/CONTRIBUTING.md>.

\begin{figure}[H]

{\centering \includegraphics[width=0.3\linewidth]{logo} 

}

\caption{Computer library logo.}\label{fig:logo}
\end{figure}

<!-- You can also specify older versions of the OpenBLAS library. Automatically, if no version is specified, the `ropenblas` package will consider the latest version of the library OpenBLAS.  -->

<!-- Considering using the OpenBLAS library rather than the BLAS may bring extra optimizations for your code and improved computational performance for your simulations, since OpenBLAS is an optimized implementation of the library BLAS. -->

<!-- Some of the reasons why it is convenient to link R language to the use of BLAS optimized alternatives can be found here. Several other benchmarks that point to improved computing performance by considering the library OpenBLAS can be found on the internet. -->

<!-- The `ropenblas` package, by `rcompiler()` function is also useful if you want to install different versions of the R language. The different versions, specified by the user of the R language, will be compiled and will also be linked to the OpenBLAS library. If you want to switch between compiled versions of the R language, no compilation is needed anymore. This allows you to avoid having to get your hands dirty with tedious operating system settings, regardless of your GNU/Linux distribution. Another great use of the `rcompiler()` function is that you will not be dependent on updating your GNU/Linux distribution repositories and you can always have the latest version of the R language. -->

<!-- The use of the `ropenblas` package will return warnings that help you proceed with the use of the functions. If your internet is not working or if any dependency on the operating system is not present, the package will let you know. -->

<!-- # Dependencies -->

<!-- In addition to dependencies in the form of other packages deposited with CRAN, the `ropenblas` package depends on external dependencies that are normally installed or are easily installed on any GNU/Linux distribution. Are they: -->

<!--   1. **GNU Make**: GNU Make utility to maintain groups of programs;  -->

<!--   2. **GNU GCC Compiler (C and Fortran)**: The GNU Compiler Collection - C and Fortran frontends. -->

<!-- These programs that are described in `SystemRequirements` in the package's `DESCRIPTION` file are essential for compiling the OpenBLAS library and the R language. The functions of the `ropenblas` package are designed to identify the lack of these dependencies external to CRAN, informing the package user which dependencies are missing and suggesting that they should be installed. The other dependencies indexed to CRAN are described in `Imports` in the file  `DESCRIPTION`. These will be installed automatically. -->

<!-- Other warnings can also be suggested, such as, for example, a problem with the internet connection. All warnings are given very clearly so that the user has no doubts about the problem that may be occurring. -->

# Brief explanation

<!-- # Installation -->

<!-- The `ropenblas` package can be installed in two ways. The first is using the `install.packages()` function of the `utils` package which is available in any basic language installation and the second is using the `devtools` package which will allow the package to be installed directly from the development directory on GitHub. -->

<!-- All code kept in the master branch of the package project on GitHub can be installed, since there will only be codes that are working properly and ready to use. The two forms of installation follow: -->

<!--   1. `install.packages("ropenblas")`: for installing the package available at CRAN; -->

<!--   2. `devtools::install_github(repo = "prdm0/ropenblas, ref = "master", force = TRUE)`: for installing the package from the project development directory on GitHub. -->

<!-- # Exported functions and usage -->

The `ropenblas` library exports six functions for use which are the `rcompiler()`, `ropenblas()`, `last_version_r()`, `last_version_openblas()`, `link_again()` and `rnews()`. All of them are very simple to use and have few arguments that are sufficient to maintain the flexibility of use. Also, functions like `rcompiler()` and `ropenblas()` do not return content or data structures that are of any practical use. What these functions do is configure the GNU/Linux system to use R, configure different versions of the language, switch between versions, and link with the OpenBLAS library. It is also possible to obtain a summary of the versions of R and the OpenBLAS library that are available.

Table 1 below presents the benefit of considering an optimized version of BLAS. Computational costs are presented in the calculation of the singular decomposition in a rectangular matrix (`svd()` function) and in the calculation of the inverse of that same matrix (`solve()` function). Some repetitions (100 repetitions) of each of the respective functions were performed. The benchmark can be better observed through the violin plots shown in Figure 2.

| Functions  | Library  | Time (seconds) |
|------------|----------|:--------------:|
| `svd(x)`   | BLAS     |     6.520      |
| `svd(x)`   | OpenBLAS |     0.641      |
| `solve(x)` | BLAS     |     1.640      |
| `solve(x)` | OpenBLAS |     0.640      |

: Comparison of the computational costs of the `svd()` and `solve()` functions (average of 100 repetitions).

Through a benchmark it is possible to better understand the performance gain that can be achieved by linking the R language to the OpenBLAS library. Figure 2 presents the benchmarks in the form of a violin plot, in which 100 reproductions of the `svd(X)` expression were considered, in the form of the code above, with the R linked to the BLAS library and linked to the OpenBLAS library, respectively, on the same hardware. It was observed that the average time of execution of the routine `svd(X)` considering the OpenBLAS library was less than 10 times the time necessary to execute it in R linking to a non-optimized version of BLAS, being the average time of 0.64 and 6.52 seconds, respectively.

\begin{figure}[H]

{\centering \includegraphics[width=0.5\linewidth]{benchmark} 

}

\caption{Benchmarks of a decomposition of singular and inverse value of a matrix of dimension 1000 x 1000.}\label{fig:bench}
\end{figure}

<!-- ## 'link_again' function -->

<!-- The `link_again()` function links again the OpenBLAS library with the R language, being useful to correct problems of untying the OpenBLAS library that is common when the operating system is updated. The function can link again the R language with the OpenBLAS library. -->

<!-- Thus, `link_again()` will only make the linkage when in some previous section of R the `ropenblas()` function has been used for the initial binding of the R language with the OpenBLAS library. -->

<!-- The use of the function is quite simple, just by running the code `link_again()` since the function has no arguments. It will automatically detect if there was a link break that will be rebuilt again without the need for any compilation. From time to time, after a major update of the operating system, it may be convenient to run the `link_again()` function. Link breakage rarely occurs, but if it does, it can be resolved quickly. The following code and image exemplify a possible reconstruction of symbolic links using the `link_again()` function: -->

<!-- ```{r, eval=FALSE, prompt=TRUE} -->

<!-- link_again() -->

<!-- ``` -->

<!-- ```{r, echo=FALSE, message=FALSE, warning=FALSE, fig.align='center', fig.cap="If an unlinking of the OpenBLAS library occurs, the function will re-link the library.", out.width="100%"} -->

<!-- knitr::include_graphics("link_again_01.png") -->

<!-- ``` -->

<!-- Running the `link_again()` function in a situation where there is no need will not generate problems. The function will return the message that everything is linked correctly, according to the code and image that follows: -->

<!-- ```{r, echo=FALSE, message=FALSE, warning=FALSE, fig.align='center', fig.cap="If relinking is not required, the function will make it clear.", out.width="100%"} -->

<!-- knitr::include_graphics("link_again_02.png") -->

<!-- ``` -->

<!-- ## 'rnews function' -->

<!-- The `rnews()` function returns the contents of the `NEWS.html` file in the standard browser installed on the operating system. The `NEWS.html` file contains the main changes from the recently released versions of the R language. The goal is to facilitate the query by invoking it directly from the R command prompt. -->

<!-- The `rnews()` function is analogous to the news function of the `utils` package. However, using the news command in a terminal style bash shell is possible to receive a message like: -->

<!-- ```{r, eval=FALSE} -->

<!-- news() -->

<!-- ## starting httpd help server ... done -->

<!-- ## Error in browseURL(url): 'browser' must be a non-empty character string -->

<!-- ``` -->

<!-- If `pdf = FALSE` (default), the `NEWS.html` file will open in the browser, otherwise `NEWS.pdf` will be opened. If `dev = FALSE` (default), it will not show changes made to the language development version. To see changes in the development version, do `dev = TRUE`. -->

# Improvements

The package will continue to evolve and code reviews will always be carried out. In addition, contributions to the development of the package are always welcome, especially those that aim to allow the `rcompiler()` and `ropenblas()` functions to work on Windows systems. There is also an interest that the `ropenblas` package will allow the linking of the Intel Math Kernel Library - MKL, just as it is done with the OpenBLAS library. All of these are improvements that we would like to see in future versions of the package.

# References
---
title: Contributing to the development of the ropenbals package
author: |
  | **Package website**: **<https://prdm0.github.io/ropenblas>**
  |
  | **Summary**
output: 
  bookdown::github_document2:
    number_sections: true
    toc: true
    toc_depth: 2
toccolor: 'blue'
---

# General information

The **ropenblas** package introduces a simple way, through its functions, to link the library [**OpenBLAS**](https://www.openblas.net/) with the language [**R**](https: //www.r-project.org/). In addition, with ropenblas it is possible to link different versions of R, making it possible to easily switch between the different versions of R compiled.

The package currently works only on GNU/Linux systems and is independent of the official repositories of the user's GNU/Linux distribution. This allows the R user on GNU/Linux systems who is the administrator of the operating system - OS to link the latest versions of R and OpenBLAS as well as switch between versions if they wish.

# Main contributions

The main contributions to the **ropenblas** package are listed below:

  1. Revision of the code to detect possible inconsistent points that may cause the functions provided by the package to fail in some GNU / Linux distribution;
  2. Test the package on different GNU / Linux distributions, reporting possible inconsistencies at <https://github.com/prdm0/ropenblas/issues>;
  3. Improved English for package messages and documentation;
  4. Improved writing of job documentation;
  5. Allow Windows system users to take advantage of the `ropenblas()` function. One suggestion would be to link precompiled binaries from the OpenBLAS library and distribute them in a `inst/` directory.---
title: "Download, Compile and Link OpenBLAS Library with R"
author: |
  |
  |
  |
  | **Package website**: **<https://prdm0.github.io/ropenblas>**
  |
  | **Summary**
output: 
  bookdown::github_document2:
    number_sections: false
    toc: true
    toc_depth: 2
toccolor: 'blue'
editor_options: 
  markdown: 
    wrap: 72
---

**ropenblas package**

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02769/status.svg)](https://doi.org/10.21105/joss.02769)
[![Travis build status](https://travis-ci.com/prdm0/ropenblas.svg?branch=master)](https://travis-ci.com/prdm0/ropenblas)
[![last](https://www.r-pkg.org/badges/last-release/ropenblas)](https://CRAN.R-project.org/package=ropenblas)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/ropenblas)](https://CRAN.R-project.org/package=ropenblas)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/ropenblas)](https://CRAN.R-project.org/package=ropenblas)
[![month](https://cranlogs.r-pkg.org/badges/ropenblas)](https://CRAN.R-project.org/package=ropenblas)

<!-- ```{r, fig.align='center', echo=FALSE, out.width=219} -->
<!-- knitr::include_graphics("https://raw.githubusercontent.com/prdm0/ropenblas/master/logo.png") -->
<!-- ``` -->

<img src="logo.png" height="230" width="200" align="right" />

# General aspects

The [**ropenblas**](https://prdm0.github.io/ropenblas/) is a package
designed to facilitate the linking of the library
[**OpenBLAS**](https://www.openblas.net/) with the language
[**R**](https://www.r-project.org/). The package, which works only for
Linux systems, will automatically download the latest source code from
the [**OpenBLAS**](https://www.openblas.net/) library and compile the
code. The package will automatically bind the language
[**R**](https://www.r-project.org/), through the `ropenblas()` function,
to use the [**OpenBLAS**](https://www.openblas.net/) library. Everything
will be done automatically regardless of the Linux distribution you are
using.

You can also specify older versions of the
[**OpenBLAS**](https://www.openblas.net/) library. Automatically, if no
version is specified, the
[**ropenblas**](https://prdm0.github.io/ropenblas/) package will
consider the latest version of the library
[**OpenBLAS**](https://www.openblas.net/).

Considering using the [**OpenBLAS**](https://www.openblas.net/) library
rather than the [**BLAS**](http://www.netlib.org/blas/) may bring extra
optimizations for your code and improved computational performance for
your simulations, since [**OpenBLAS**](https://www.openblas.net/) is an
optimized implementation of the library
[**BLAS**](http://www.netlib.org/blas/).

Some of the reasons why it is convenient to link
[**R**](https://www.r-project.org/) language to the use of
[**BLAS**](http://www.netlib.org/blas/) optimized alternatives can be
found [**here**](https://csantill.github.io/RPerformanceWBLAS/). Several
other [**benchmarks**](https://en.wikipedia.org/wiki/Benchmarking) that
point to improved computing performance by considering the library
[**OpenBLAS**](https://www.openblas.net/) can be found on the internet.

The ropenblas package, by `rcompiler()` function is also useful if you
want to install different versions of the
[**R**](https://www.r-project.org/) language. The different versions,
specified by the user of the [**R**](https://www.r-project.org/)
language, will be compiled and will also be linked to the
[**OpenBLAS**](https://www.openblas.net/) library. If you want to switch
between compiled versions of the [**R**](https://www.r-project.org/)
language, no compilation is needed anymore. This allows you to avoid
having to get your hands dirty with tedious operating system settings,
regardless of your GNU/Linux distribution. Another great use of the
`rcompiler()` function is that you will not be dependent on updating
your GNU/Linux distribution repositories and you can always have the
latest version of the [**R**](https://www.r-project.org/) language.

The use of the ropenblas package will return warnings that help you
proceed with the use of the functions. If your internet is not working
or if any dependency on the operating system is not present, the package
will let you know.

# Advantages of using ropenblas package

Some advantages of using the
[**ropenblas**](https://prdm0.github.io/ropenblas/) library:

-   Everything is done within the [**R**](https://www.r-project.org/)
    language;

-   The procedure will be the same for any Linux distribution;

-   The [**OpenBLAS**](https://www.openblas.net/) library will be
    compiled and you will choose which build version to bind to
    [**R**](https://www.r-project.org/), regardless of your Linux
    distribution;

-   If your GNU/Linux distribution does not have updated versions of
    [**OpenBLAS**](https://www.openblas.net/), it matters little. The
    ropenblas package fetches the latest stable release of the
    [**OpenBLAS**](https://www.openblas.net/) library development
    account on GitHub;

-   You do not need to know Linux well. In some distributions, it may
    not be so simple for a less experienced user to compile and link the
    library to the [**OpenBLAS**](https://www.openblas.net/) library
    with the [**R**](https://www.r-project.org/) language;

-   It is much easier to direct a person to link
    [**OpenBLAS**](https://www.openblas.net/) with
    [**R**](https://www.r-project.org/) saying "run `ropenblas()` within
    [**R**](https://www.r-project.org/)" than asking that person to
    verify that an unoptimized version of
    [**BLAS**](http://www.netlib.org/blas/) installed on the system.
    Then you have to guide the removal of the unoptimized version of
    [**BLAS**](http://www.netlib.org/blas/) and guide it to the
    installation of the library
    [**OpenBLAS**](https://www.openblas.net/) through the most diverse
    procedures depending on the GNU/Linux distribution used;

<!-- - As stated earlier, the procedure works for any Linux and this includes Android. If your Android is capable of running privileged commands (ROOT) and if you have [**R**](https://www.r-project.org/) installed via Termux with the required dependencies, you can compile and link [**OpenBLAS**](https://www.openblas.net/) with [**R**](https://www.r-project.org/) using ropenblas; -->

-   With the `rcompiler()` function you can build any version of
    [**R**](https://www.r-project.org/) into your computer architecture,
    which includes the most stable version of the language;

-   What is the latest stable version of
    [**R**](https://www.r-project.org/)? Are you too lazy to go to the
    [**R**](https://www.r-project.org/) site? Run
    `ropenblas::last_version_r(major = NULL)`.

# Dependencies

You must install the following dependencies on your operating system
(Linux):

1 - **GNU Make**: GNU Make utility to maintain groups of programs; <br/>

2 - **GNU GCC Compiler (C and Fortran)**: The GNU Compiler Collection -
C and Fortran frontends.

Do not worry that you will be notified if any of these dependencies are
not installed.

# Installation

Installing the [**ropenblas**](https://prdm0.github.io/ropenblas/)
library is easy and will require you to have installed the **devtools**
package. This will allow you to install the
[**ropenblas**](https://prdm0.github.io/ropenblas/) package directly
from GitHub. To install, after installing the **devtools** package, do:

    remotes::install_github(repo = "prdm0/ropenblas", force = TRUE)

or

    install.packages("ropenblas")

<!-- **Note**: If you want to access the latest features of the [**ropenblas**](https://prdm0.github.io/ropenblas/) package, install it using the first procedure. -->

# Use

The [**ropenblas**](https://prdm0.github.io/ropenblas/) package
currently provides five functions: `ropenblas()`, `rcompiler()`,
`last_version_openblas()`, `last_version_r()`, `link_again()` and
`rnews()`. First of all, do:

    library(ropenblas)

## 'ropenblas' function

Installing, compiling, and linking the
[**OpenBLAS**](https://www.openblas.net/) version **0.3.13** library to
the [**R**](https://www.r-project.org/) language:

```{r, fig.align='center',  echo=FALSE}
knitr::include_graphics(path = "https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/ropenblas.svg")
```

<!-- **Notes**:  -->

<!--    - You do not have to in every section of [**R**](https://www.r-project.org/) make use of the `ropenblas()` function. Once the function is used, [**R**](https://www.r-project.org/) will always consider using the [**OpenBLAS**](https://www.openblas.net/) library in future sections. -->

<!--    - [**OpenBLAS**](https://www.openblas.net/) versions tested: 0.3.0, 0.3.1, 0.3.2, 0.3.3, 0.3.4, 0.3.5, 0.3.6, 0.3.7, 0.3.8, 0.3.9, 0.3.10, 0.3.11, 0.3.12 and 0.3.13. These are the values that will be passed to `x` in `ropenblas(x)`;  -->

<!--    - If `x = NULL`, the latest stable version of the [**OpenBLAS**](https://www.openblas.net/) library will be compiled and linked to [**R**](https://www.r-project.org/). -->

<!-- **Details** -->

<!-- Your linux operating system may already be configured to use the [**OpenBLAS**](https://www.openblas.net/) library. Therefore, most likely [**R**](https://www.r-project.org/) will already be linked to this library. To find out if the [**R**](https://www.r-project.org/) language is using the [**OpenBLAS**](https://www.openblas.net/) library, at [**R**](https://www.r-project.org/), do: -->

<!-- ``` -->

<!-- extSoftVersion()["BLAS"] -->

<!-- ``` -->

<!-- If [**R**](https://www.r-project.org/) is using the [**OpenBLAS**](https://www.openblas.net/) library, something like `/any_directory/libopenblas.so` should be returned. Therefore, there should be the name **openblas** in the **s**hared **o**bject returned (file extension **.so**). -->

<!-- If the `ropenblas()` function can identify that the [**R**](https://www.r-project.org/) language is using the version of [**OpenBLAS**](https://www.openblas.net/) you wish to configure, a warning message will be returned asking if you really would like to proceed with the configuration again. -->

<!-- The `ropenblas()` function will download the desired version of the library [**OpenBLAS**](https://www.openblas.net/), compile and install the library in the `/opt` directory of your operational system. If the directory does not exist, it will be created so that the installation can be completed. Subsequently, files from the version of [**BLAS**](http://www.netlib.org/blas/) used in [**R**](https://www.r-project.org/) will be symbolically linked to the shared object files of the library version [**OpenBLAS**](https://www.openblas.net/) compiled and installed in `/opt`. -->

<!-- You must be the operating system administrator to use this library. Therefore, do not attempt to use it without telling your system administrator. If you have the ROOT password, you will be responsible for everything you do on your operating system. -->

<!-- You will not necessarily have to run `ropenblas()` on every section of [**R**](https://www.r-project.org/). Almost always it will not be necessary. However, it may be that the [**R**](https://www.r-project.org/) is updated by the operating system (GNU/Linux). Thus, it may be that in this update the [**R**](https://www.r-project.org/) unlink with the [**OpenBLAS**](https://www.openblas.net/) library. Therefore, from time to time check using the command `extSoftVersion()["BLAS"]` if the link with [**OpenBLAS**](https://www.openblas.net/) is correct, otherwise run the command `ropenblas()` again. -->

## 'last_version_r' function

Given the higher version, the function will return the latest stable
version of the [**R**](https://www.r-project.org/) language. See the
following example:

```{r, echo=FALSE}
knitr::include_graphics(path = "https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/last_version_r.svg")
```

<!-- **Note**: If `major = NULL`, the function will consider the major release number. -->

## 'last_version_openblas' function

The `last_version_openblas()` function automatically searches
[**OpenBLAS**](https://www.openblas.net/) library versions in the
official [**GitHub**](https://github.com/xianyi/OpenBLAS) project.

```{r, echo=FALSE}
knitr::include_graphics(path = "https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/last_version_openblas.svg")
```

## 'rcompiler' function

This function is responsible for compiling a version of the
[**R**](https://www.r-project.org/) language. The `x` argument is the
version of [**R**](https://www.r-project.org/) that you want to compile.
For example, `x = "4.0.4"` will compile and link **R-4.0.4** version as
the major version on your system. By default (`x = NULL`) will be
compiled the latest stable version of the
[**R**](https://www.r-project.org/).

For example, to compile the latest stable version of the
[**R**](https://www.r-project.org/) language, do:

```{r, echo=FALSE}
knitr::include_graphics(path = "https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/rcompiler.svg")
```

Regardless of your GNU/Linux distribution and what version of
[**R**](https://www.r-project.org/) is in your repositories, you can
have the latest stable version of the
[**R**](https://www.r-project.org/) language compiled into your computer
architecture.

You can use the `rcompiler()` function to compile different versions of
[**R**](https://www.r-project.org/). For example, running
`rcompiler(x = "3.6.3")` and `rcompiler()` will install versions 3.6.3
and 4.0.0 on its GNU/Linux distribution, respectively. If you are in
version 4.0.0 of [**R**](https://www.r-project.org/) and run the code
`rcompiler(x = "3.6.3")` again, the function will identify the existence
of version 3.6.3 in the system and give you the option to use the
binaries that were built in a previous compilation. This avoids
unnecessarys compilations.

<!-- In addition to the `x` argument, the` rcompiler()` function has two other arguments that will allow you to change and pass new compilation flags. Are they: -->

<!-- 1. `with_blas`: This argument sets the `--with-blas` flag in the R language compilation process and must be passed as a string. Details on the use of this flag can be found [**here**](https://cran.r-project.org/doc/manuals/r-devel/R-admin.html). If `with_blas = NULL` (default), then it will be considered: -->

<!--    ``` -->

<!--    ./configure --prefix=/opt/R/version_r --enable-memory-profiling --enable-R-shlib --enable-threads=posix -->

<!--    --with-blas="-L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include -lpthread -lm" -->

<!--    ``` -->

<!--    Most likely, you will have little reason to change this argument. Unless you know what you're doing, consider `with_blas = NULL`. Do not change the installation directory, that is, always consider `--prefix = /opt/R/version_r`, where` version_r` is a valid version of [**R**](https://www.r-project.org/). For a list of valid versions of [**R**](https://www.r-project.org/), run the `last_version_r()`. -->

<!-- Installing [**R**](https://www.r-project.org/) in the `/opt/R/version_r` directory is important because some functions in the package require this. Both the [**R**](https://www.r-project.org/) language and the [**OpenBLAS**](https://www.openblas.net/) library will be installed in the `/opt` directory. If this directory does not exist in your GNU/Linux distribution, it will be created; -->

<!-- 2. `complementary_flags`: String (`complementary_flags = NULL` by default) for adding complementary flags in the [**R**](https://www.r-project.org/) language compilation process. Passing a string to `complementary_flags` will compile it in the form: -->

<!--    ``` -->

<!--     ./configure --with-blas="..." complementary_flags -->

<!--    ``` -->

## 'link_again' function

The `link_again` function links again the
[**OpenBLAS**](https://www.openblas.net/) library with the
[**R**](https://www.r-project.org/) language, being useful to correct
problems of untying the [**OpenBLAS**](https://www.openblas.net/)
library that is common when the operating system is updated.

The function `link_again` be able to link again the
[**R**](https://www.r-project.org/) language with the
[**OpenBLAS**](https://www.openblas.net/) library. Thus, link_again will
only make the relinkagem when in some previous section of
[**R**](https://www.r-project.org/) the ropenblas function has been used
for the initial binding of the [**R**](https://www.r-project.org/)
language with the [**OpenBLAS**](https://www.openblas.net/) library.

For example, to relink the [**OpenBLAS**](https://www.openblas.net/)
library with the [**R**](https://www.r-project.org/) language, do:

```{r, echo=FALSE}
knitr::include_graphics(path = "https://raw.githubusercontent.com/prdm0/ropenblas/bce920c05f8daa13fc6568e6797d4b841e9d43ce/img/asciicast/link_again.svg")
```

If `restart_r = TRUE` (default), a new section of
[**R**](https://www.r-project.org/) is started after linking the
[**OpenBLAS**](https://www.openblas.net/) library.

In situations where there was a disconnection due to an update of the
operating system, the `ropenblas` function can be used to relink the
[**OpenBLAS**](https://www.openblas.net/) library with the
[**R**](https://www.r-project.org/) language, however, it will be
necessary to compile the [**OpenBLAS**](https://www.openblas.net/)
library again. If you are interested in recompiling the
[**OpenBLAS**](https://www.openblas.net/) library and linking with
[**R**](https://www.r-project.org/), use the `ropenblas` function. If
the interest is to take advantage of a previous compilation of the
[**OpenBLAS**](https://www.openblas.net/) library, the function
`link_again` may be useful.

## 'rnews' function

Returns the contents of the
[**NEWS.html**](https://cran.r-project.org/doc/manuals/r-release/NEWS.html)
file in the standard browser installed on the operating system. The
[**NEWS.html**](https://cran.r-project.org/doc/manuals/r-release/NEWS.html)
file contains the main changes from the recently released versions of
the [**R**](https://www.r-project.org/) language. The goal is to
facilitate the query by invoking it directly from the
[**R**](https://www.r-project.org/) command prompt. The rnews function
is analogous to the news function of the utils package. However, using
the news command in a terminal style bash shell is possible to receive a
message like:

    > news()
    starting httpd help server ... done
    Error in browseURL(url) : 'browser' must be a non-empty character string

If `pdf = FALSE` (default), the
[**NEWS.html**](https://cran.r-project.org/doc/manuals/r-release/NEWS.html)
file will open in the browser, otherwise
[**NEWS.pdf**](https://cran.r-project.org/doc/manuals/r-release/NEWS.pdf)
will be opened. If `dev = FALSE` (default), it will not show changes
made to the language development version. To see changes in the
development version, do `dev = TRUE`.
---
# Example from https://joss.readthedocs.io/en/latest/submitting.html
title: 'ropenblas: Download, Compile and Link OpenBLAS Library with R'
tags:
  - R
  - Compiling R and OpenBLAS
  - Link OpenBLAS
  - Switch between versions of R
  - Fast algebraic computing
authors:
  - name: Pedro Rafael Diniz Marinho
    orcid: 0000-0003-1591-8300
    affiliation: 1
affiliations:
 - name: Department of Statistics, Federal University of Paraíba, João Pessoa, Paraíba - PB, Brazil
   index: 1
citation_author: Marinho
date: "`r Sys.Date()`"
year: "`r format(Sys.Date(), '%Y')`"
bibliography: paper.bib
output: rticles::joss_article
header-includes:
  - \usepackage{float}
csl: apa.csl
journal: JOSS
---

# Summary

```{r global_options, include=FALSE}
knitr::opts_chunk$set(fig.pos = 'H')
```

<!-- The `ropenblas` package aims to facilitate the day-to-day life of R programmers who want more performance on GNU/Linux systems, without removing the possibility that specific configurations are made, if they deem convenient, that is, more technical users will be able to pass other flags that will be considered in the compilation process. Through the package's `ropenblas()` and `rcompiler()` functions, the library user will be able to compile and link the R language in his GNU/Linux distribution with the OpenBLAS library, all within R and in a very simple fashion. All functions work without being influenced by the GNU/Linux distribution and are independent of their repositories, that is, it does not matter which GNU/Linux distribution is being used. Linking the OpenBLAS library to R will bring better computational performance to the language in the most diverse algebraic operations commonly used in areas such as statistics, data science, and machine learning. -->

The `ropenblas` package provides improved performance on GNU/Linux systems and allows for passing additional compilation flags for more technical users. Through the `ropenblas()` and `rcompiler()` functions, the user will be able to compile and link the GNU/Linux distribution R language with the OpenBLAS library, all within R and in a very simple fashion. It works for all GNU/Linux distributions. Linking the OpenBLAS library to R brings better computational performance to the language in the most diverse algebraic operations commonly used in areas such as statistics, data science, and machine learning.

<!-- # Statement of Need -->

<!-- The `ropenblas` package aims to allow algebraic computing of R to be performed using the OpenBLAS library, that is, it allows easy linking of R to the OpenBLAS library on GNU/Linux systems without depending on the distribution repositories. This will allow several researchers in the areas of statistics, data science, and machine learning to take advantage of more efficient performance in algebraic calculations, for example, multiplication, factorization, and matrix inversion. The `ropenblas` library, version 0.2.9 will also allow the R programmer to have, in his GNU/Linux distribution, several compiled versions of the R language giving the possibility to easily switch between these versions, this being an Open Source functionality that is only possible in some commercial IDEs of R. All this is done within the R language, minimizing the chance of less experienced users to break their operating system by running several instructions that are not perfectly understood. -->

<!-- The fact that the `ropenblas` package does not depend on the repositories of the GNU/Linux distribution will allow that in more stable distributions the R programmer will have at his disposal the most recent version of the OpenBLAS libraries and the R language. Everything is done safely, since the `ropenblas` package uses the stable versions of the official development repositories of the OpenBLAS library and the R programming language, respectively. Until the present version, the package has more than 11000 downloads, having an average of more than 1000 downloads in the month before the date of this submission, based on the official R language repositories. -->

# Introduction

<!-- The term "computational efficiency" is very common for those who program statistical methods, in which a large part of them involve algebraic operations that are often reproduced in computationally intensive simulations, such as Monte-Carlo simulations - MC and resampling methods, as is the case with bootstrap resampling. Statistics is just one example within so many other areas that need performance and uses the R language. -->

Computational efficiency is important for those who program statistical methods since they often involve algebraic operations reproduced in computationally intensive simulations, such as Monte-Carlo simulations and resampling methods, as is the case with bootstrap resampling. Statistics is just one example within other areas that need good performance and use the R language.

In addition to the adoption of good programming practices and the maximum, efficient and adequate use of available computational resources, such as code parallelization, through multicore parallelism procedures allowed by most current processors and operating systems, small adjustments and linkage of libraries can provide useful benefits.

<!-- The `ropenblas` package aims to provide useful and simple experiences to R [@R] programmers who develop their activities on GNU/Linux operating systems. These experiences consist of being able to link any version of the OpenBLAS [@openblas] library to the R language, as well as allowing the programmer to install and link various versions of R and make them available on their operating system as well as switch between these versions as they see fit. -->

The `ropenblas` package aims to be easy to use for R [@R] programmers who work on GNU/Linux operating systems. For example, a user can link any version of the OpenBLAS [@openblas] library to the R language and install and link various versions of R to make them available on their operating system as well as switch between these versions.

Linking the R language to the OpenBLAS library can bring several benefits to algebraic computing in R. OpenBLAS is an Open-Source implementation of the Basic Linear Algebra Subprograms - BLAS library that is often the first library option for algebraic computing to be linked in the installation of R on many GNU/Linux distributions. The OpenBLAS library is available at <https://github.com/xianyi/OpenBLAS> and adds optimized implementations of linear algebra kernels that can run optimized on various processor architectures. OpenBLAS is based on the GotoBLAS2 project code in version 1.13 [@gotoblas2], code available under the terms of the BSD license.

<!-- The functions of the `ropenblas` package can help the R language user to make this link without leaving the R promprt, as well as allowing the user to choose the version of OpenBLAS to be considered, the last stable version being considered by default . Everything is accomplished by functions without many arguments and without major complications to be used. These functions can be very comforting for users of R on GNU/Linux systems who do not feel safe to run several codes in a terminal that runs Shell Script codes for the language configuration. -->

<!-- It is very common to see users of R in distributions with repositories not prone to immediate updates of programs using several tutorials found in communities on the internet suggesting several lines of code and steps that can be potentially dangerous or easy to be misunderstood by many users of R language. Being able to compile, enable resources if necessary and switch between versions of R using simple functions to be used and without leaving the command pompt of the R language is attractive. -->

<!-- # General package information -->

The `ropenblas` is a package designed to facilitate the linking of the library OpenBLAS with the language R. The package, which works only for Linux systems, will automatically download the latest source code from the OpenBLAS library and compile the code. The package will automatically bind the language R, through the `ropenblas()` function, to use the OpenBLAS library. Everything will be done automatically regardless of the Linux distribution you are using. Enumerating some advantages of the package:

1.  Everything is done within the R language;
2.  The procedure (use of functions) will be the same for any Linux distribution;
3.  The OpenBLAS library will be compiled and you will choose which build version to bind to R, regardless of your Linux distribution;
4.  The package allows you to install R $\geq 3.1.0$, also allowing you to install one more version, in addition to allowing you to easily switch between those versions;
5.  The linked versions of R will continue to be recognized by their Integrated Development Environment - IDE and nothing will have to be adjusted in your GNU/Linux distribution after using any function of the package;
6.  Unnecessary builds will be avoided. Therefore, if you need to switch between compiled versions of the R language, the use of binaries compiled compiled previously will be suggested;
7.  If any errors occur, the functions of the package will not damage the previous installation of the language;
8.  If something better can be done or if a newer version of what you want to install (R or OpenBLAS) exists, the functions will automatically suggest that you should consider installing newer versions.

The `ropenblas` package is already available on the Comprehensive R Archive Network - CRAN, currently in version 0.2.9, and the project is maintained on GitHub at <https://github.com/prdm0/ropenblas> where contributors can find other details of the code, information, as well as being able to contribute with the development of the project. On the website, it is also possible to read the `NEWS.md` file with details of the versions and the focus of the current development. The site is deposited at <https://prdm0.github.io/ropenblas/>. Suggestions for improvements and bug reports can be sent via the link <https://github.com/prdm0/ropenblas/issues>. You can find out how to contribute to the package by accessing the `CONTRIBUTING.md` file at <https://github.com/prdm0/ropenblas/blob/master/CONTRIBUTING.md>.

```{r logo, echo = FALSE, message = FALSE, fig.cap = "Computer library logo.", fig.align = "center", out.width="30%"}
knitr::include_graphics(path = "logo.png") 
```

<!-- You can also specify older versions of the OpenBLAS library. Automatically, if no version is specified, the `ropenblas` package will consider the latest version of the library OpenBLAS.  -->

<!-- Considering using the OpenBLAS library rather than the BLAS may bring extra optimizations for your code and improved computational performance for your simulations, since OpenBLAS is an optimized implementation of the library BLAS. -->

<!-- Some of the reasons why it is convenient to link R language to the use of BLAS optimized alternatives can be found here. Several other benchmarks that point to improved computing performance by considering the library OpenBLAS can be found on the internet. -->

<!-- The `ropenblas` package, by `rcompiler()` function is also useful if you want to install different versions of the R language. The different versions, specified by the user of the R language, will be compiled and will also be linked to the OpenBLAS library. If you want to switch between compiled versions of the R language, no compilation is needed anymore. This allows you to avoid having to get your hands dirty with tedious operating system settings, regardless of your GNU/Linux distribution. Another great use of the `rcompiler()` function is that you will not be dependent on updating your GNU/Linux distribution repositories and you can always have the latest version of the R language. -->

<!-- The use of the `ropenblas` package will return warnings that help you proceed with the use of the functions. If your internet is not working or if any dependency on the operating system is not present, the package will let you know. -->

<!-- # Dependencies -->

<!-- In addition to dependencies in the form of other packages deposited with CRAN, the `ropenblas` package depends on external dependencies that are normally installed or are easily installed on any GNU/Linux distribution. Are they: -->

<!--   1. **GNU Make**: GNU Make utility to maintain groups of programs;  -->

<!--   2. **GNU GCC Compiler (C and Fortran)**: The GNU Compiler Collection - C and Fortran frontends. -->

<!-- These programs that are described in `SystemRequirements` in the package's `DESCRIPTION` file are essential for compiling the OpenBLAS library and the R language. The functions of the `ropenblas` package are designed to identify the lack of these dependencies external to CRAN, informing the package user which dependencies are missing and suggesting that they should be installed. The other dependencies indexed to CRAN are described in `Imports` in the file  `DESCRIPTION`. These will be installed automatically. -->

<!-- Other warnings can also be suggested, such as, for example, a problem with the internet connection. All warnings are given very clearly so that the user has no doubts about the problem that may be occurring. -->

# Brief explanation

<!-- # Installation -->

<!-- The `ropenblas` package can be installed in two ways. The first is using the `install.packages()` function of the `utils` package which is available in any basic language installation and the second is using the `devtools` package which will allow the package to be installed directly from the development directory on GitHub. -->

<!-- All code kept in the master branch of the package project on GitHub can be installed, since there will only be codes that are working properly and ready to use. The two forms of installation follow: -->

<!--   1. `install.packages("ropenblas")`: for installing the package available at CRAN; -->

<!--   2. `devtools::install_github(repo = "prdm0/ropenblas, ref = "master", force = TRUE)`: for installing the package from the project development directory on GitHub. -->

<!-- # Exported functions and usage -->

The `ropenblas` library exports six functions for use which are the `rcompiler()`, `ropenblas()`, `last_version_r()`, `last_version_openblas()`, `link_again()` and `rnews()`. All of them are very simple to use and have few arguments that are sufficient to maintain the flexibility of use. Also, functions like `rcompiler()` and `ropenblas()` do not return content or data structures that are of any practical use. What these functions do is configure the GNU/Linux system to use R, configure different versions of the language, switch between versions, and link with the OpenBLAS library. It is also possible to obtain a summary of the versions of R and the OpenBLAS library that are available.

Table 1 below presents the benefit of considering an optimized version of BLAS. Computational costs are presented in the calculation of the singular decomposition in a rectangular matrix (`svd()` function) and in the calculation of the inverse of that same matrix (`solve()` function). Some repetitions (100 repetitions) of each of the respective functions were performed. The benchmark can be better observed through the violin plots shown in Figure 2.

| Functions  | Library  | Time (seconds) |
|------------|----------|:--------------:|
| `svd(x)`   | BLAS     |     6.520      |
| `svd(x)`   | OpenBLAS |     0.641      |
| `solve(x)` | BLAS     |     1.640      |
| `solve(x)` | OpenBLAS |     0.640      |

: Comparison of the computational costs of the `svd()` and `solve()` functions (average of 100 repetitions).

Through a benchmark it is possible to better understand the performance gain that can be achieved by linking the R language to the OpenBLAS library. Figure 2 presents the benchmarks in the form of a violin plot, in which 100 reproductions of the `svd(X)` expression were considered, in the form of the code above, with the R linked to the BLAS library and linked to the OpenBLAS library, respectively, on the same hardware. It was observed that the average time of execution of the routine `svd(X)` considering the OpenBLAS library was less than 10 times the time necessary to execute it in R linking to a non-optimized version of BLAS, being the average time of 0.64 and 6.52 seconds, respectively.

```{r bench, echo = FALSE, fig.cap = "Benchmarks of a decomposition of singular and inverse value of a matrix of dimension 1000 x 1000.", fig.align = "center", out.width="50%"}
knitr::include_graphics(path = "benchmark.png") 
```

<!-- ## 'link_again' function -->

<!-- The `link_again()` function links again the OpenBLAS library with the R language, being useful to correct problems of untying the OpenBLAS library that is common when the operating system is updated. The function can link again the R language with the OpenBLAS library. -->

<!-- Thus, `link_again()` will only make the linkage when in some previous section of R the `ropenblas()` function has been used for the initial binding of the R language with the OpenBLAS library. -->

<!-- The use of the function is quite simple, just by running the code `link_again()` since the function has no arguments. It will automatically detect if there was a link break that will be rebuilt again without the need for any compilation. From time to time, after a major update of the operating system, it may be convenient to run the `link_again()` function. Link breakage rarely occurs, but if it does, it can be resolved quickly. The following code and image exemplify a possible reconstruction of symbolic links using the `link_again()` function: -->

<!-- ```{r, eval=FALSE, prompt=TRUE} -->

<!-- link_again() -->

<!-- ``` -->

<!-- ```{r, echo=FALSE, message=FALSE, warning=FALSE, fig.align='center', fig.cap="If an unlinking of the OpenBLAS library occurs, the function will re-link the library.", out.width="100%"} -->

<!-- knitr::include_graphics("link_again_01.png") -->

<!-- ``` -->

<!-- Running the `link_again()` function in a situation where there is no need will not generate problems. The function will return the message that everything is linked correctly, according to the code and image that follows: -->

<!-- ```{r, echo=FALSE, message=FALSE, warning=FALSE, fig.align='center', fig.cap="If relinking is not required, the function will make it clear.", out.width="100%"} -->

<!-- knitr::include_graphics("link_again_02.png") -->

<!-- ``` -->

<!-- ## 'rnews function' -->

<!-- The `rnews()` function returns the contents of the `NEWS.html` file in the standard browser installed on the operating system. The `NEWS.html` file contains the main changes from the recently released versions of the R language. The goal is to facilitate the query by invoking it directly from the R command prompt. -->

<!-- The `rnews()` function is analogous to the news function of the `utils` package. However, using the news command in a terminal style bash shell is possible to receive a message like: -->

<!-- ```{r, eval=FALSE} -->

<!-- news() -->

<!-- ## starting httpd help server ... done -->

<!-- ## Error in browseURL(url): 'browser' must be a non-empty character string -->

<!-- ``` -->

<!-- If `pdf = FALSE` (default), the `NEWS.html` file will open in the browser, otherwise `NEWS.pdf` will be opened. If `dev = FALSE` (default), it will not show changes made to the language development version. To see changes in the development version, do `dev = TRUE`. -->

# Improvements

The package will continue to evolve and code reviews will always be carried out. In addition, contributions to the development of the package are always welcome, especially those that aim to allow the `rcompiler()` and `ropenblas()` functions to work on Windows systems. There is also an interest that the `ropenblas` package will allow the linking of the Intel Math Kernel Library - MKL, just as it is done with the OpenBLAS library. All of these are improvements that we would like to see in future versions of the package.

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ropenblas.R
\name{link_again}
\alias{link_again}
\title{Linking the OpenBLAS library with \R again}
\usage{
link_again(restart_r = TRUE)
}
\arguments{
\item{restart_r}{If \code{TRUE} (default), a new section of \R is started after linking the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library.}
}
\description{
The \code{link_again} function links again the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library with the \R language, being useful to correct problems
of untying the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library that is common when the operating system is updated.
}
\details{
The function \code{link_again} be able to link again the \R language with the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library. Thus, link_again will only make the
relinkagem when in some previous section of \R the ropenblas function has been used for the initial binding of the \R language with the
\href{https://www.openblas.net/}{\strong{OpenBLAS}} library.

Relinking is useful in situations of updating the operating system. In some update, it is possible that the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library compiled
in the \code{/opt} directory is unlinked. In this scenario, when the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library has already been compiled using the ropenblas
function, the \code{link_again} function performs a new link without the need to recompile, thus making the process less time
consuming.
}
\note{
In situations where there was a disconnection due to an update of the operating system, the \code{ropenblas} function can be used to
re-link the OpenBLAS library with the \R language, however, it will be necessary to compile the
\href{https://www.openblas.net/}{\strong{OpenBLAS}} library again. If you are interested in recompiling the
\href{https://www.openblas.net/}{\strong{OpenBLAS}} library and linking with \R, use the \code{\link{ropenblas}} function. If the
interest is to take advantage of a previous compilation of the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library, the
function \code{link_again} may be useful.
}
\examples{
# link_again()
}
\seealso{
\code{\link{ropenblas}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ropenblas.R
\name{rnews}
\alias{rnews}
\title{R News file}
\usage{
rnews(pdf = FALSE, dev = FALSE)
}
\arguments{
\item{pdf}{If \code{FALSE} (default), the \href{https://cran.r-project.org/doc/manuals/r-release/NEWS.html}{NEWS.html} file will open in the browser, otherwise
\href{https://cran.r-project.org/doc/manuals/r-release/NEWS.pdf}{NEWS.pdf} will be opened.}

\item{dev}{If \code{FALSE} (default), it will not show changes made to the language development version.
To see changes in the development version, do \code{dev = TRUE}.}
}
\description{
Returns the contents of the \href{https://cran.r-project.org/doc/manuals/r-release/NEWS.html}{NEWS.html} file in the standard browser installed on the operating system.
}
\details{
The \href{https://cran.r-project.org/doc/manuals/r-release/NEWS.html}{NEWS.html} file contains the main changes from the recently released versions of the \R language.
The goal is to facilitate the query by invoking it directly from the \R command prompt. The \link{rnews} function is
analogous to the \link{news} function of the \strong{utils} package. However, using the \link{news} command in a terminal style
bash shell is possible to receive a message like:\preformatted{news()
starting httpd help server ... done
Error in browseURL(url) : 'browser' must be a non-empty character string
}

This is an error that may occur depending on the installation of \R. Always prefer the use of the news
function but if you need to, use the \link{rnews} function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ropenblas.R
\name{last_version_r}
\alias{last_version_r}
\title{\R language versions}
\usage{
last_version_r(major = NULL)
}
\arguments{
\item{major}{Major release number of \R language (e.g. \code{1L}, \code{2L}, \code{3L}, ...). If \code{major = NULL}, the function
will consider the major release number.}
}
\value{
A list of two named elements will be returned. Are they:
\enumerate{
\item \code{last_version}: Returns the latest stable version of the language given a major version (major version).
If \code{major = NULL}, the latest stable version of the language will be returned based on the set of all language versions.
\item \code{versions}: Character vector with all language versions based on a major version (higher version).
If \code{major = NULL}, \code{versions} will be a vector with the latest language versions.
\item \code{n}: Total number of versions of \R based on major version. If \code{major = NULL}, \code{versions}
will be a vector with the latest language versions.
}
}
\description{
\R language versions
}
\details{
This function automatically searches \R language versions in the official language repositories. That way,
doing \code{last_version_r(major = NULL)} you will always be well informed about which latest stable version the
\R language is in. You can also set the higher version and do a search on the versions of the \R language whose major
version was \code{1L} or \code{2L}, for example.
}
\examples{
# last_version_r(major = NULL)
}
\seealso{
\code{\link{ropenblas}}, \code{\link{rcompiler}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ropenblas.R
\name{ropenblas}
\alias{ropenblas}
\title{Download, Compile and Link OpenBLAS Library with \R}
\usage{
ropenblas(x = NULL, restart_r = TRUE)
}
\arguments{
\item{x}{\href{https://www.openblas.net/}{\strong{OpenBLAS}} library version to be considered. By default, \code{x = NULL}.}

\item{restart_r}{If \code{TRUE}, a new section of \R is started after compiling and linking the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library.}
}
\value{
Returns a warning message informing you if the procedure occurred correctly. You will also be able to receive information about
missing dependencies.
}
\description{
Link \R with an optimized version of the \href{http://www.netlib.org/blas/}{\strong{BLAS}} library (\href{https://www.openblas.net/}{\strong{OpenBLAS}}).
}
\details{
The \code{ropenblas()} function will only work on Linux systems. When calling the \code{ropenblas()}
function on Windows, no settings will be made. Only a warning message will be issued informing you that the
configuration can only be performed on Linux systems.

The function will automatically download the latest version of the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library. However, it is possible to
inform olds versions to the single argument of \code{ropenblas()}. The \code{ropenblas()} function downloads,
compiles and link \R to use of the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library. Everything is done very simply, just loading the library and
invok the function \code{ropenblas()}.

Considering using the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library rather than the \href{http://www.netlib.org/blas/}{\strong{BLAS}} may bring extra optimizations for your code and improved
computational performance for your simulations, since \href{https://www.openblas.net/}{\strong{OpenBLAS}} is an optimized implementation of the library \href{http://www.netlib.org/blas/}{\strong{BLAS}}.

You must install the following dependencies on your operating system (Linux):
\enumerate{
\item \strong{GNU Make};
\item \strong{GNU GCC Compiler (C and Fortran)}.
}
Your linux operating system may already be configured to use the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library. Therefore, most likely \R will already be linked to this library. To find out if the \R language is using the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library,
at \R, do:

\code{extSoftVersion()["BLAS"]}

If \R is using the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library, something like \code{/any_directory/libopenblas.so} should be returned. Therefore, there should be the name openblas in the \strong{s}hared \strong{o}bject returned (file extension \strong{.so}).

If the \code{ropenblas()} function can identify that the \R language is using the version of \href{https://www.openblas.net/}{\strong{OpenBLAS}} you wish to configure, a warning message will be returned asking if you really would like to proceed with the
configuration again.

The \code{ropenblas()} function will download the desired version of the library \href{https://www.openblas.net/}{\strong{OpenBLAS}}, compile and install the library in the \code{/opt} directory of your operational system. If the directory does not exist, it will
be created so that the installation can be completed. Subsequently, files from the version of \href{http://www.netlib.org/blas/}{\strong{BLAS}} used in \R will be symbolically linked to the shared object files of the library version \href{https://www.openblas.net/}{\strong{OpenBLAS}} compiled and installed in \code{/opt}.

You must be the operating system administrator to use this library. Therefore, do not attempt to use it without telling your system administrator. If you have the ROOT password, you will be responsible for everything you do on your operating system. Other details you may also find \href{https://prdm0.github.io/ropenblas/index.html}{\strong{here}}.
}
\note{
You do not have to in every section of \R make use of the \code{ropenblas()} function. Once the function is used, \R
will always consider using the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library in future sections.
}
\examples{
# ropenblas()
}
\seealso{
\code{\link{rcompiler}}, \code{\link{last_version_r}}
}
\author{
Pedro Rafael D. Marinho (e-mail: \email{pedro.rafael.marinho@gmail.com})
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ropenblas.R
\name{last_version_openblas}
\alias{last_version_openblas}
\title{OpenBLAS library versions}
\usage{
last_version_openblas()
}
\description{
OpenBLAS library versions
}
\details{
This function automatically searches \href{https://www.openblas.net/}{\strong{OpenBLAS}} library versions in the official \href{https://github.com/xianyi/OpenBLAS}{\strong{GitHub}} project.
\enumerate{
\item \code{last_version}: Returns the latest stable version of the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library.
\item \code{versions}: All stable versions of the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library.
\item \code{n}: Total number of versions.
}
}
\examples{
# last_version_openblas()
}
\seealso{
\code{\link{last_version_r}}, \code{\link{ropenblas}}, \code{\link{rcompiler}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ropenblas.R
\name{rcompiler}
\alias{rcompiler}
\title{Compile a version of \R on GNU/Linux systems}
\usage{
rcompiler(x = NULL, with_blas = NULL, complementary_flags = NULL)
}
\arguments{
\item{x}{Version of \R you want to compile. By default (\code{x = NULL}) will be compiled the latest stable version of the \R
language. For example, \code{x = "3.6.2"} will compile and link \strong{R-3.6.2} version  as the major version on your system.}

\item{with_blas}{String, \code{--with-blas = NULL} by default, with flags for \code{--with-blas} used in the \R compilation process.
Details on the use of this flag can be found \href{https://cran.r-project.org/doc/manuals/r-devel/R-admin.html}{\strong{here}}.}

\item{complementary_flags}{String, \code{complementary_flags = NULL} by default, for adding complementary flags in the \href{https://www.r-project.org/}{\strong{R}} language
compilation process.}
}
\value{
Returns a warning message informing you if the procedure occurred correctly. You will also be able to receive information about
missing dependencies.
}
\description{
This function is responsible for compiling a version of the \R language.
}
\details{
This function is responsible for compiling a version of the \href{https://www.r-project.org/}{\strong{R}} language. The \code{x} argument is the version of \href{https://www.r-project.org/}{\strong{R}} that you want to compile.
For example, \code{x = "4.0.0"} will compile and link \strong{R-4.0.0} version  as the major version on your system. By default (\code{x = NULL}) will be compiled the latest stable version of the \href{https://www.r-project.org/}{\strong{R}}.

For example, to compile the latest stable version of the \href{https://www.r-project.org/}{\strong{R}} language, do:\preformatted{ rcompiler()
}

Regardless of your GNU/Linux distribution and what version of \href{https://www.r-project.org/}{\strong{R}} is in your repositories, you can have the latest stable version of the \href{https://www.r-project.org/}{\strong{R}} language compiled
into your computer architecture.

You can use the \code{rcompiler()} function to compile different versions of \href{https://www.r-project.org/}{\strong{R}}. For example, running \code{rcompiler(x = "3.6.3")} and \code{rcompiler()} will install versions 3.6.3 and 4.0.0 on its GNU/Linux distribution,
respectively. If you are in version 4.0.0 of \href{https://www.r-project.org/}{\strong{R}} and run the code \code{rcompiler(x = "3.6.3")} again, the function will identify the existence of version 3.6.3 in the system and give you the option to use the binaries
that were built in a previous compilation. This avoids unnecessarys compilations.

In addition to the \code{x} argument, the\code{ rcompiler()} function has two other arguments that will allow you to change and pass new compilation flags. Are they:
\enumerate{
\item \code{with_blas}: This argument sets the \code{--with-blas} flag in the R language compilation process and must be passed as a string. Details on the use of this flag
can be found \href{https://cran.r-project.org/doc/manuals/r-devel/R-admin.html}{\strong{here}}. If \code{with_blas = NULL} (default), then it will be considered:\preformatted{./configure --prefix=/opt/R/version_r --enable-memory-profiling --enable-R-shlib
 --enable-threads=posix --with-blas="-L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include
 -lpthread -lm"
}

Most likely, you will have little reason to change this aprgument. Unless you know what you're doing, consider \code{with_blas = NULL}. Do not change the installation directory,
that is, always consider \verb{--prefix = /opt/R/version_r}, where\code{ version_r} is a valid version of \href{https://www.r-project.org/}{\strong{R}}. For a list of valid versions of
\href{https://www.r-project.org/}{\strong{R}}, run the \code{last_version_r()}. Installing \href{https://www.r-project.org/}{\strong{R}} in the \verb{/opt/R/version_r} directory is important because some
functions in the package require this. Both the \href{https://www.r-project.org/}{\strong{R}} language and the \href{https://www.openblas.net/}{\strong{OpenBLAS}} library will be installed in the \verb{/opt} directory.
If this directory does not exist in your GNU/Linux distribution, it will be created;
\item \code{complementary_flags}: String (\code{complementary_flags = NULL} by default) for adding complementary flags in the \href{https://www.r-project.org/}{\strong{R}} language compilation process.
Passing a string to \code{complementary_flags} will compile it in the form:\preformatted{./configure --with-blas="..." complementary_flags
}
}
}
\examples{
# rcompiler()
}
\seealso{
\code{\link{ropenblas}}, \code{\link{last_version_r}}
}
