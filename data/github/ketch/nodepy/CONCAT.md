---
title: '`NodePy`: A package for the analysis of numerical ODE solvers'
tags:
  - Python
  - numerical analysis
  - differential equations
  - Runge-Kutta method
  - linear multistep method
authors:
  - name: David I. Ketcheson^[Corresponding author.]
    orcid: 0000-0002-1212-126X
    affiliation: 1
  - name: Hendrik Ranocha
    orcid: 0000-0002-3456-2277
    affiliation: 1
  - name: Matteo Parsani
    orcid: 0000-0001-7300-1280
    affiliation: 1
  - name: Umair bin Waheed
    orcid: 0000-0002-5189-0694
    affiliation: 2
  - name: Yiannis Hadjimichael
    orcid: 0000-0003-3517-8557
    affiliation: 3
affiliations:
 - name: King Abdullah University of Science and Technology
   index: 1
 - name: King Fahd University of Petroleum & Minerals
   index: 2
 - name: Eötvös Loránd Tudományegyetem
   index: 3
date: 9 July 2020
bibliography: paper.bib
---

# Summary

Ordinary differential equations (ODEs) are used to model a vast range of physical
and other phenomena.  They also arise in the discretization of partial differential
equations.  In most cases, solutions of differential equations must be approximated
by numerical methods.  The study of the properties of numerical methods for
ODEs comprises an important and large body of knowledge.  `NodePy`
(available from https://github.com/ketch/nodepy, with documentation at https://nodepy.readthedocs.io/en/latest/)
is a software
package for designing and studying the properties of numerical ODE solvers.
For the most important classes of methods, `NodePy` can automatically assess
their stability, accuracy, and many other properties.
`NodePy` has also been used as a catalog of coefficients for time integration methods
in PDE solver codes.

# Statement of need

There are many software packages that *implement* ODE solvers with the purpose
of efficiently providing numerical solutions; in contrast, the purpose of
`NodePy` is to facilitate understanding of the properties of the solver algorithms
themselves.  In this sense, it is a sort of meta-software, consisting of
algorithms whose purpose is to compute properties of other algorithms.
It also serves as a reference, providing precise definitions of many of the
algorithms themselves.

`NodePy` is written entirely in Python and provides software implementations
of many of the theoretical ideas contained for instance in reference texts
on numerical analysis of ODEs [@hairer1993;@Hairer:ODEs2].  It also contains implementations of
many theoretical ideas from the numerical analysis literature.
The implementation focuses on the two most important classes of methods;
namely, Runge-Kutta and linear multistep methods, but includes some
more exotic classes.  `NodePy` provides a means for numerical analysts to
quickly and easily determine the properties of existing methods or of new
methods they may develop.

`NodePy` development has been motivated largely by research needs and
it has been used in a number of papers (including some written by non-developers;
e.g. @jin2019higher and @horvathembedded) and also as a teaching tool for
graduate-level numerical analysis courses.  It relies on both SymPy [@meurer2017sympy]
and NumPy [@oliphant2006guide;@walt2011numpy]
in order to provide either exact or floating-point results based on the
nature of the inputs provided.  It makes use of Matplotlib [@hunter2007matplotlib] for all
graphical output.

# Features

`NodePy` includes object-oriented representations of the following classes
of numerical methods:

 - Runge-Kutta methods
   - Explicit and implicit
   - Embedded pairs
   - Classes of low-storage methods
   - Dense output formulas
   - Perturbed/additive and downwind methods
 - Linear multistep methods
 - Two-step Runge-Kutta methods
 - Additive (IMEX) linear multistep methods

The framework is designed to include general linear methods and even more
exotic classes.  Any method within these classes can be generated simply by entering
its coefficients.  Coefficients for many methods are catalogued or can be
automatically generated within `NodePy`, including:

 - Dozens of specific Runge-Kutta methods and pairs
 - General extrapolation methods, of any order of accuracy, based on a variety
   of building-block schemes and optionally including an error estimator
 - Deferred correction methods
 - Optimal strong stability preserving (SSP) Runge-Kutta methods
 - Adams-Bashforth, Adams-Moulton, and BDF methods of any order
 - A number of other specialized families of methods

For all of these numerical schemes, `NodePy` provides methods and functions to compute many
of their properties -- too many to list here.  The theory on which most of these
properties are based is outlined in standard references
[@hairer1993;@Hairer:ODEs2].  Many other properties are based on
recent research; usually the method docstring includes a reference to the relevant
paper.  Implementations of the methods themselves are also included as a convenience,
though they are not the primary purpose and are not expected to be efficient since
they are coded in pure Python.  Additional intermediate objects, such as the
absolute stability function of a method, are given their own software representation
and corresponding methods.

Additional features are provided to facilitate the analysis and testing of
these numerical methods.  This includes a range of initial value problems
for testing, such as the stiff and non-stiff DETEST suites, and a few simple
PDE semi-discretizations.  Also included is a library for dealing with rooted
trees, which are a class of graphs that play a key role in the theory of Runge-Kutta
methods.

`NodePy` is documented primarily through the Python docstrings for each method,
most of which contain examples that also serve as tests ("doctests").
These tests are executed automatically for all new commits and pull requests
using the Travis continuous integration service.  Some higher-level
documentation is also available at https://nodepy.readthedocs.io/en/latest/,
but it is not intended to be comprehensive.

# Related research and software

We are not aware of other software packages with a similar purpose.
`NodePy` development has proceeded in close connection to the RK-Opt package
[@ketcheson2020RK-Opt] (https://github.com/ketch/RK-Opt).
Whereas `NodePy` is focused in the analysis of numerical methods, RK-Opt is focused
more on their design through the use of numerical optimization to search
for optimal coefficients tailored to specific desired properties.
A common workflow involves generating new methods with RK-Opt and then studying
their properties in more detail using `NodePy`.

Some of the research projects that have made use of `NodePy` (most of which have led
to its further development) include development of:

 - Strong stability preserving (SSP) Runge-Kutta methods
   [@2008_explicit_ssp;@2009_implicit_ssp;@2013_effective_order_ssp]
 - SSP general linear methods [@2011_tsrk;@2017_msrk]
 - Low-storage Runge-Kutta methods [@2010_LSRK]
 - Additive and downwind SSP Runge-Kutta methods [@2011_dwssp;@2018_perturbations]
 - High-order parallel extrapolation and deferred correction methods [@2014_hork]
 - SSP linear multistep methods [@2016_ssp_lmm_vss;@2018_sspalmm]
 - Dense output formulas for Runge-Kutta methods [@2017_dense]
 - Internal stability theory for Runge-Kutta methods [@2014_internal_stability]
 - Embedded pairs for Runge-Kutta methods [@horvathembedded;@conde2018embedded]

Additional recent applications include [@norton2015structure;@jin2019higher;@ranocha2019some;@2019_energyRRK;@nusslein2020positivity].
As an example of a completely different kind of use,
in the fluid dynamics code SpectralDNS, `NodePy` is used simply for convenience as a way to enable
usage of a range of ODE solvers; here `NodePy` is used only for
retrieving the coefficients and not for the time-stepping implementation.
This facilitated the work in [@ketcheson2020more], for instance.
As can be seen from this list, applications have mostly stemmed from the
work of the main developer's research group, but have recently begun to expand
beyond that.

# Acknowledgements

Much of the initial `NodePy` development was performed by D. Ketcheson while
he was supported by a DOE Computational Science Graduate Fellowship.  Development
has also been supported by funding from King Abdullah University of Science and Technology.
Additional minor contributions to the code have been provided by Mikael Mortensen, Alex Fikl,
Sidafa Conde, John Sellers, Kevin Siswandi, and Colin Macdonald.

# References
If you use NodePy in a published work, please cite it as follows:

    Ketcheson, D. I.,  NodePy software version <version number>,
    http://github.com/ketch/nodepy/.

Please insert the version number that you used.

# NodePy: A package for the analysis of numerical ODE solvers

[![Build Status](https://travis-ci.com/ketch/nodepy.png)](https://travis-ci.com/ketch/nodepy)
[![Coverage Status](https://coveralls.io/repos/github/ketch/nodepy/badge.svg?branch=master)](https://coveralls.io/github/ketch/nodepy?branch=master)
[![codecov.io](https://codecov.io/github/ketch/nodepy/coverage.svg?branch=master)](https://codecov.io/github/ketch/nodepy?branch=master)

[![](https://readthedocs.org/projects/nodepy/badge)](https://readthedocs.org/projects/nodepy/)
[![version status](https://pypip.in/v/nodepy/badge.png)](https://pypi.python.org/pypi/nodepy)
[![downloads](https://pypip.in/d/nodepy/badge.png)](https://pypi.python.org/pypi/nodepy)

[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-success.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4275157.svg)](https://doi.org/10.5281/zenodo.4275157)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02515/status.svg)](https://doi.org/10.21105/joss.02515)


# Installation
NodePy requires Python 3.5 or later.  To install with pip, do:

    pip install nodepy

This will automatically fetch dependencies also.  It will not fetch
optional dependencies, which include networkx, cvxpy and scipy (that are used
only in a few specialized routines and/or examples).  The optional dependencies
can be installed with `pip`.

# Overview

NodePy (Numerical ODEs in Python) is a Python package for designing, analyzing,
and testing numerical methods for initial value ODEs. Its development was
motivated by my own research in time integration methods for PDEs. I found that
I was frequently repeating tasks that could be automated and integrated.
Initially I developed a collection of MATLAB scripts, but this became unwieldy
due to the large number of files that were necessary and the more limited
capability for code reuse.

NodePy represents an object-oriented approach, in which the basic object is a
numerical ODE solver. The idea is to design a laboratory for such methods in
the same sense that MATLAB is a laboratory for matrices.

Documentation can be found online at

http://nodepy.readthedocs.org/en/latest/

To get started, you can also have a look at the `examples` folder,
beginning with an [introduction as Jupyter notebook](examples/Introduction%20to%20NodePy.ipynb).

The development version can be obtained from

http://github.com/ketch/nodepy

# Citation

If you use NodePy in a published work, please cite it as follows:

    Ketcheson, D. I.  NodePy software version <version number>,
    http://github.com/ketch/nodepy/.

Please insert the version number that you used.

# Support

If you encounter an error or need help, please [raise an issue](https://github.com/ketch/nodepy/issues).

# Contributing

Contributions of new features or other improvements are very welcome!  Please
[submit a pull request](https://github.com/ketch/nodepy/pulls) or contact the authors.

# License

NodePy is distributed under the terms of the [modified Berkeley Software
Distribution (BSD) license](LICENSE.txt).


# Funding

NodePy development has been supported by:

* A U.S. Dept. of Energy Computational Science Graduate Fellowship
* Grants from King Abdullah University of Science & Technology


The following people have contributed to Nodepy:

- Umair bin Waheed: Embedded pairs of extrapolation methods
- Matteo Parsani: Internal stability functions for Runge-Kutta methods
- Sidafa Conde: Pretty-printing of TSRKs

Thank you!
