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
.. contents::

.. _create_rkm:

Runge-Kutta methods
============================

A Runge-Kutta method is a one-step method that computes the next time
step solution as follows:

    \\begin{align*}
    y_i = & u^{n} + \\Delta t \\sum_{j=1}^{s} + a_{ij} f(y_j)) & (1\\le j \\le s) \\\\
    u^{n+1} = & u^{n} + \\Delta t \\sum_{j=1}^{s} b_j f(y_j).
    \\end{align*}

The simplest way to load a Runge-Kutta method is using the 
loadRKM function::

    >> from nodepy import runge_kutta_method as rk
    >> import numpy as np
    >> rk44=rk.loadRKM('RK44')
    Classical RK4

     0.000 |
     0.500 |  0.500
     0.500 |  0.000  0.500
     1.000 |  0.000  0.000  1.000
    _______|________________________________
           |  0.167  0.333  0.333  0.167

Many well-known methods are available through the loadRKM() function.
Additionally, several classes of methods are available through the
following functions:

  * Optimal strong stability preserving methods: SSPRK2(s), SSPRK3(s), SSPIRK2(), etc.
  * Integral deferred correction methods: 
    :mod:`DC(s) <nodepy.runge_kutta_method.DC>`

  * Extrapolation methods: 
    :mod:`extrap(s) <nodepy.runge_kutta_method.extrap>`
  * Runge-Kutta Chebyshev methods: RKC1(s), RKC2(s)

See the documentation of these functions for more details.

More generally, any Runge-Kutta method may be instantiated by providing 
its Butcher coefficients, $A$ and $b$::

    >> A=np.array([[0,0],[0.5,0]])
    >> b=np.array([0,1.])
    >> rk22=rk.RungeKuttaMethod(A,b)

Note that, because NumPy arrays are indexed from zero, the Butcher coefficient
$a_{21}$, for instance, corresponds to my_rk.a[1,0].
The abscissas $c$ are automatically set to the row sums of $A$ (this
implies that every stage has has stage order at least equal to 1).
Alternatively, a method may be specified in Shu-Osher form, by coefficient
arrays $\alpha,\beta$::

    >> rk22=rk.RungeKuttaMethod(alpha=alpha,beta=beta)

A separate subclass is provided for explicit Runge-Kutta methods: *ExplicitRungeKuttaMethod*.
If a method is explicit, it is important to instantiate it as an
ExplicitRungeKuttaMethod and not simply a RungeKuttaMethod, since the
latter class has significantly less functionality in NodePy.
Most significantly, time stepping is currently implemented for explicit methods, 
but not for implicit methods.

.. automodule:: nodepy.runge_kutta_method
   :noindex:

Accuracy
--------------------------------------------
The principal measure of accuracy of a Runge-Kutta method is its
*order of accuracy*.  By comparing the Runge-Kutta solution with
the Taylor series for the exact solution, it can be shown that the
local truncation error for small enough step size $h$ is approximately

Error `\approx Ch^p`,

where $C$ is a constant independent of $h$.
Thus the expected asymptotic rate of convergence for small
step sizes is $p$.  This error corresponds to the lowest-order terms
that do not match those of the exact solution Taylor series.
Typically, a higher order accurate method will provide greater
accuracy than a lower order method even for practical step sizes.

In order to compare two methods with the same order of accuracy
more detailed information about the accuracy of a method may be
obtained by considering the relative size of the constant $C$.
This can be measured in various ways and is referred to as the 
principal error norm.

For example::

    >> rk22.order()
    2
    >> rk44.order()
    4
    >> rk44.principal_error_norm()
    0.014504582343198208
    >> ssp104=rk.loadRKM('SSP104')
    >> ssp104.principal_error_norm()
    0.002211223747053554
    
Since the SSP(10,4) method has smaller principal error norm, we
expect that it will provide better accuracy than the classical 4-stage
Runge-Kutta method for a given step size.  Of course, the SSP method
has 10 stages, so it requires more work per step.  In order to
determine which method is more efficient, we need to compare the relative
accuracy for a fixed amount of work::

    >> rk.relative_accuracy_efficiency(rk44,ssp104)
    1.7161905294239843

This indicates that, for a desired level of error, the SSP(10,4)
method will require about 72\% more work.

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.principal_error_norm
   :noindex:

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.error_metrics
   :noindex:

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.stage_order
   :noindex:

Classical (linear) stability
--------------------------------------------
.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.stability_function
   :noindex:

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.plot_stability_region
   :noindex:

Nonlinear stability
--------------------------------------------
.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.absolute_monotonicity_radius
   :noindex:

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.circle_contractivity_radius
   :noindex:

Reducibility of Runge-Kutta methods
--------------------------------------------
Two kinds of reducibility (*DJ-reducibility* and *HS-reducibility*) have
been identified in the literature.  NodePy contains functions for detecting
both and transforming a reducible method to an equivalent irreducible method.
Of course, reducibility is dealt with relative to some numerical tolerance,
since the method coefficients are floating point numbers.

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod._dj_reducible_stages
   :noindex:

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.dj_reduce
   :noindex:

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod._hs_reducible_stages
   :noindex:

Composing Runge-Kutta methods
--------------------------------------------

Butcher has developed an elegant theory of the group structure of 
Runge-Kutta methods.  The Runge-Kutta methods form a group under the
operation of composition.  The multiplication operator has been 
overloaded so that multiplying two Runge-Kutta methods gives the
method corresponding to their composition, with equal timesteps.

It is also possible to compose methods with non-equal timesteps using the
compose() function.


Embedded Runge-Kutta Pairs
============================
.. autoclass:: nodepy.runge_kutta_method.ExplicitRungeKuttaPair
   :noindex:


Low-Storage Runge-Kutta methods
=================================
.. automodule:: nodepy.low_storage_rk
   :noindex:

2S/3S methods
--------------------
.. autoclass:: nodepy.low_storage_rk.TwoSRungeKuttaMethod
   :noindex:

2S/3S embedded pairs
--------------------
.. autoclass:: nodepy.low_storage_rk.TwoSRungeKuttaPair
   :noindex:

2R/3R methods
--------------------
.. autoclass:: nodepy.low_storage_rk.TwoRRungeKuttaMethod
   :noindex:


.. _changes:

What's new in version 1.0
=========================
- Improvements to low-storage RK methods (thanks to @ranocha)
- Updates for compatibility with newer versions of Sympy
- Embedded pairs for SSP methods (thanks to @ranocha)
- Some additional tests
- Many improvements to the documentation and examples, suggested by reviewer @fruzsinaagocs
- NodePy no longer officially supports Python 2.x, although virtually everything still works in Python 2.7.

What's new in version 0.9
=========================
- Implicit extrapolation Runge-Kutta methods
- Updates for compatibility with Sympy
- Test for algebraic stability of RK methods
- Conversion from Shu-Osher to Butcher form now also works for implicit methods
- Add coefficients of two additional RK4(3) embedded pairs

What's new in version 0.8
=========================

- Compute E-polynomial for RK methods
- Plot step size controller stability for explicit RK pairs
- Now possible to integrate complex solutions
- Many new specific ODE methods added
- Fixed some bugs in generation of trees and order conditions of very high order
- Added some new example notebooks
- Compatibility with Matplotlib 3

What's new in version 0.7
=========================
*Released November 29, 2016*

- Support for Python 3 (thanks to Github user @alexfikl)
- Dense output for Runge-Kutta methods.
- Removal of a circular dependency.

What's new in Version 0.6.1
===========================
*Released May 14, 2015*

- Two algorithms for computing optimal downwind perturbations of Runge-Kutta methods.  One relies on CVXPY for solving linear programs.
- The Numipedia project has been moved to its own repository: https://github.com/ketch/numipedia
- Many new doctests; >80% test coverage.
- Many improvements to the two-step RK module, including stability region plots for arbitrary methods.
- Three-step RK methods removed from master (because most of the module was not working).
- Pretty-printing of linear multistep methods.
- New methods:
  - Several very-high-order RK methods
  - Some singly diagonally-implicit RK methods
  - Nystrom and Milne-Simpson families of multistep methods
- load_ivp() works similarly to load_RKM() (returns a dictionary).
- Improved computation of maximum linearly stable step sizes for semi-discretizations.
- Many bug fixes.

What's new in Version 0.6
==========================
Version 0.6 is a relatively small update, which includes the following:

- Computation of optimal perturbations (splittings) of Runge-Kutta methods
- Additive linear multistep methods
- More accurate calculation of imaginary stability intervals
- Rational coefficients for more of the built-in RK methods
- Faster computation of stability polynomials
- More general deferred correction methods
- Fixed major bug in deferred correction method construction
- Continuous integration via Travis-CI
- Added information on citing nodepy
- Corrections to the documentation
- Updates for compatitibility with sympy 0.7.6
- Fixed bug involving non-existence of alphahat attribute
- minor bug fixes



What's new in Version 0.5
==========================
*Released: Nov. 4, 2013*

Version 0.5 is a relatively small update, which includes the following:

* More Runge-Kutta methods available in rk.loadRKM(), including the 8(7) Prince-Dormand pair
* Lots of functionality and improvements for studying internal stability of RK methods
* Shu-Osher arrays used to construct an RK method are now stored and (by default) used for timestepping
* Ability to compute the effective order of a RK method
* More accurate computation of stability region intervals
* Use exact arithmetic (sympy) in many more functions
* Generation of Fortran code for order conditions
* Refactoring of how embedded Runge-Kutta pairs and low-storage methods are represented internally
* Plotting functions return a figure handle
* Better pretty-printing of RK methods with exact coefficients
* Updates for compatibility with sympy 0.7.3
* Improved reducibility for RK pairs
* More initial value problems
* Several bug fixes
* Automated testing with Travis

What's new in Version 0.4
==========================
*Released: Aug. 28, 2012*

Version 0.4 of NodePy inclues numerous bug fixes and new features.
The most significant new feature is the use of exact arithmetic for
construction and analysis of many methods, using SymPy.  Because exact
arithmetic can be slow, NodePy automatically switches to floating point
arithmetic for some operations, such as numerical integration of initial value
problems.  If you find operations that seem excessively slow let me know.
You can always revert to floating-point representation of a method by
using method.__num__().

Other new features and fixes include:

    * Improvements to linear multistep methods:
        * Stability region plotting
        * Zero-stability
        * `A(\alpha)`-stability angles
    * Automatic selection of plotting region for stability region plots
    * Code base now hosted on Github (github.com/ketch/nodepy)
    * Documentation corrections
    * Use MathJax (instead of jsMath) in docs
    * Much greater docstring coverage
    * Many more examples in docs (can be run as doctests)
        * For example, 95 doctests covering 25 items in runge_kutta_method.py
    * Extrapolation methods based on GBS (midpoint method) -- thanks to Umair bin Waheed
    * Construction of simple linear finite difference matrices
    * Analysis of the potential for parallelism in Runge-Kutta methods
        * Number of sequentially-dependent stages
        * Plotting of stage dependency graph
    * Automatic reduction of reducible Runge-Kutta methods
    * A heuristic method for possibly-optimal splittings of Runge-Kutta methods
      into upwind/downwind parts
    * Fix bugs in computation of stability intervals
    * Fix bugs in stability region plotting
    * New examples in nodepy/examples/
    * Spectral difference matrices for linear advection -- thanks to Matteo Parsani


=================================================
Testing Methods: Solving Initial Value Problems
=================================================
In addition to directly analyzing solver properties, NodePy also facilitates
the testing of solvers through application to problems of interest.
Furthermore, NodePy includes routines for automatically running sets of 
tests to compare the practical convergence or efficiency of various solvers
for given problem(s).


.. contents::


.. _ivp:

Initial Value Problems
==============================

The principal objects in NodePy are ODE solvers.  The object
upon which a solver acts is an initial value problem.  Mathematically,
an initial value problem (IVP) consists of one or more ordinary 
differential equations and an initial condition:

    \\begin{align*}
    u'(t) & = F(u) & u(0) & = u_0.
    \\end{align*}


In NodePy, 
an initial value problem is an object with the following properties:

    * rhs(): The right-hand-side function; i.e. F where $u(t)'=F(u)$.
    * u0:  The initial condition.
    * T:   The (default) final time of solution.

Optionally an IVP may possess the following:
    * exact(): a function that takes one argument (t) and returns the exact solution (Should we make this a function of u0 as well?)
    * dt0: The default initial timestep when a variable step size integrator is used.
    * Any other problem-specific parameters.

The module ivp contains functions for loading a variety of initial
value problems.  For instance, the van der Pol oscillator problem
can be loaded as follows::

    >> from NodePy import ivp
    >> myivp = ivp.load_ivp('vdp')
    

Instantiation
-----------------

.. automethod:: nodepy.ivp.detest


Solving Initial Value Problems
====================================

Any ODE solver object in NodePy can be used to solve an initial value
problem simply by calling the solver with an initial value problem object
as argument::

    >> t,u = rk44(my_ivp)


.. automethod:: nodepy.ode_solver.ODESolver.__call__


Convergence Testing
--------------------------------

.. automethod:: nodepy.convergence.ctest

Performance testing with automatic step-size control
----------------------------------------------------------------

.. automethod:: nodepy.convergence.ptest

.. contents::

.. _create_lmm:

Linear Multistep methods
==================================
A linear multistep method computes the next solution value from the values
at several previous steps:

    `\alpha_k y_{n+k} + \alpha_{k-1} y_{n+k-1} + ... + \alpha_0 y_n
    = h ( \beta_k f_{n+k} + ... + \beta_0 f_n )`

Note that different conventions for numbering the coefficients exist;
the above form is used in NodePy.
Methods are automatically normalized so that $\\alpha_k=1$.

.. automodule:: nodepy.linear_multistep_method
   :noindex:

------------------------
Instantiation
------------------------
The follwing functions return linear multistep methods of some 
common types:

  * Adams-Bashforth methods: Adams_Bashforth(k)
  * Adams-Moulton methods: Adams_Moulton(k)
  * backward_difference_formula(k)
  * Optimal explicit SSP methods (elm_ssp2(k))

In each case, the argument $k$ specifies the number of steps in the method.
Note that it is possible to generate methods for arbitrary $k$, but currently
for large $k$ there are large errors in the coefficients due to roundoff errors.
This begins to be significant at 7 steps.  However, members of these families
with many steps do not have good properties.

More generally, a linear multistep method can be instantiated by specifying
its coefficients $\\alpha,\\beta$::

    >> from nodepy import linear_multistep_method as lmm
    >> my_lmm=lmm.LinearMultistepMethod(alpha,beta)


Adams-Bashforth Methods
------------------------

.. automethod:: nodepy.linear_multistep_method.Adams_Bashforth
   :noindex:

Adams-Moulton Methods
------------------------

.. automethod:: nodepy.linear_multistep_method.Adams_Moulton
   :noindex:

Backward-difference formulas
------------------------------

.. automethod:: nodepy.linear_multistep_method.backward_difference_formula
   :noindex:

Optimal Explicit SSP methods
------------------------------

.. automethod:: nodepy.linear_multistep_method.elm_ssp2
   :noindex:

------------------------
Stability
------------------------

Characteristic Polynomials
------------------------------
.. automethod:: nodepy.linear_multistep_method.LinearMultistepMethod.characteristic_polynomials
   :noindex:

Plotting The Stability Region
------------------------------
.. automethod:: nodepy.linear_multistep_method.LinearMultistepMethod.plot_stability_region
   :noindex:
==============================
Quick Start Guide
==============================

.. .. contents::

Obtaining NodePy
================
The recommended installation method is via pip.
The current development version of NodePy can be obtained via Git::
    
    git clone git://github.com/ketch/nodepy.git


NodePy Documentation
====================

NodePy documentation can be found at 
http://numerics.kaust.edu.sa/nodepy

The documentation is also included in the nodepy/doc directory, and can
be built from your local install, if you have Sphinx.

Examples
====================

NodePy comes with some canned examples that can be run to confirm
your installation and to demonstrate capabilities of NodePy.
These can be found in the directory :file:`nodepy/examples`.
================================
Classes of ODE solvers
================================

The basic object in NodePy is an ODE solver.  Several types of solvers are
supported, including linear multistep methods, Runge-Kutta methods, and
Two-step Runge-Kutta methods.  For each class, individual methods
may be instantiated by specifying their coefficients.  
Many convenience functions are also provided for loading common methods
or common families of methods.  The Runge-Kutta method class also supports
some special low-storage RK method classes, integral deferred correction 
methods, and extrapolation methods.

.. toctree::

    rkm
    lmm
    tsrkm
.. _create_tsrkm:

Two-step Runge-Kutta methods
======================================

Two-step Runge-Kutta methods are a class of multi-stage multistep methods
that use two steps and (potentially) several stages.

.. autoclass:: nodepy.twostep_runge_kutta_method.TwoStepRungeKuttaMethod


==============================
Analyzing Stability Properties
==============================

Plotting the region of absolute stability
=========================================

Region of absolute stability for the optimal SSP 10-stage, 4th order
Runge-Kutta method:

.. plot::
   :include-source:

   from nodepy.runge_kutta_method import *
   ssp104=loadRKM('SSP104')
   ssp104.plot_stability_region(bounds=[-15,1,-10,10])

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.plot_stability_region
   :noindex:

Region of absolute stability for the 3-step Adams-Moulton method:

.. plot::
   :include-source:

   from nodepy.linear_multistep_method import *
   am3=Adams_Moulton(3)
   am3.plot_stability_region()

.. automethod:: nodepy.linear_multistep_method.LinearMultistepMethod.plot_stability_region
   :noindex:



Plotting the order star
=========================================

Order star for the optimal SSP 10-stage, 4th order
Runge-Kutta method:

.. plot::
   :include-source:

   from nodepy.runge_kutta_method import *
   ssp104=loadRKM('SSP104')
   ssp104.plot_order_star()

.. automethod:: nodepy.runge_kutta_method.RungeKuttaMethod.plot_stability_region
   :noindex:
================================
Planned Future Development
================================

Development of NodePy is based on research needs.  The following 
is a list of capabilities or features that would naturally fit 
into the package but have not yet been implemented.

Multistep methods
---------------------------

    * Time-stepping for multistep methods
    * Selection of startup method
    * Variable step size multistep methods
    * Properties of particular multistep + startup method combinations
    * Adaptive step size and order selection

Runge-Kutta Methods
---------------------------

    * Time stepping for implicit methods
    * Interpolants (dense output)
    * Adaptive order (extrapolation and deferred correction)


PDEs
---------------------------

Many common semi-discretizations of PDEs will be implemented as
`ivp` objects.  Initially this will be implemented purely in Python and
limited to simple 1D PDEs (e.g. advection, diffusion), since 
time-stepping for multi-dimensional or nonlinear PDEs will be too
slow in Python.  Eventually we plan to support wrapped Fortran and C
semi-discretizations.

Miscellaneous
---------------------------
 * Additional classes of multi-stage, multistep methods.
 * Analysis of geometric integrators.
 * Partitioned rooted trees, additive and partitioned RK methods.
 * Unit tests.
 * An automatically-generated encyclopedia of solvers.
.. NodePy documentation master file, created by
   sphinx-quickstart on Mon Mar  9 20:46:17 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. contents::

================
Overview
================
NodePy (Numerical ODEs in Python) is a Python package for designing,
analyzing, and testing numerical methods for initial value ODEs.
Its development was motivated by my own research in time integration
methods for PDEs.  I found that I was frequently repeating tasks that
could be automated and integrated.  Initially I developed a collection
of MATLAB scripts, but this became unwieldy due to the large number
of files that were necessary and the more limited capability for
code reuse.  

NodePy represents an object-oriented approach, in which the basic
object is a numerical ODE solver.  The idea is to design a laboratory for
such methods in the same sense that MATLAB is a laboratory for matrices.
Some distinctive design goals are:

  * **Plug-and-play**: any method can be applied to any problem using the
    same syntax.  Also, properties of different kinds of methods are 
    available through the same syntax.  This makes it easy to compare 
    different methods.
  * **Abstract representations**: Generally, the most abstract
    (hence powerful) representaton of an object is used whenever
    possible.  Thus, order conditions are generated using products
    on rooted trees (or other recursions) rather than being hard-coded.
  * **Numerical representation**: The most precise representation possible
    is used for quantities such as coefficients: rational numbers (using
    SymPy's Rational class) when available, floating-point numbers otherwise.
    Where necessary, method properties are determined by numerical
    calculations, using appropriate tolerances.  Thus the "order" of
    a method with floating-point coefficients is determined by checking whether
    the order conditions are
    satisfied to within a small value (near machine-epsilon).
    For efficiency reasons, coefficients are always converted to floating-point
    for purposes of applying the method to a problem.

In general, user-friendliness of the interface and readability of
the code are prioritized over performance.

NodePy includes capabilities for applying the methods to solve
systems of ODEs.  This is mainly intended for testing and comparison;
for realistic problems of interest in most fields, time-stepping in
Python will be too slow.  One way around this is to wrap Fortran or
C functions representing the right-hand-side of the ODE, and we are
looking into this.

.. note:: The online documentation is not comprehensive.
          For more complete documentation, it is best to
          refer to the docstrings of specific functions and classes.


Dependencies
================================
  * Works with Python 3.5+
  * Requires: Numpy, Matplotlib, Sympy
  * SymPy (note: NodePy is now compatible with SymPy 0.7.1)
  * Optional: networkx (for some Runge-Kutta stage dependency graphing), 
    cvxpy (for finding optimal downwind perturbations),
    scipy

Classes of Numerical ODE Solvers
================================

NodePy includes classes for the following types of methods:
    - :ref:`Runge-Kutta Methods <create_rkm>`
        - Implicit
        - Explicit
        - Embedded pairs
        - Low-storage methods
        - Extrapolation methods
        - Integral deferred correction methods
        - Strong stability preserving methods
        - Runge-Kutta-Chebyshev methods
    - :ref:`Linear Multistep Methods <create_lmm>`
    - :ref:`Two-step Runge-Kutta Methods <create_tsrkm>`

Arbitrary methods in these classes can be instantiated by specifying
their coefficients.


Analysis of Methods
===================

NodePy includes functions for analyzing many properties, including:
    - Stability:
        - Absolute stability (e.g., plot the region of absolute stability)
        - Strong stability preservation
    - Accuracy
        - Order of accuracy
        - Error coefficients
        - Relative accuracy efficiency
        - Generation of Python and MATLAB code for order conditions


Testing Methods
======================

NodePy includes implementation of the actual time-stepping algorithms
for the various classes of methods.  A wide range of 
:ref:`initial value ODEs <ivp>` can be loaded, including the DETEST suite of problems.
Arbitrary initial value problems
can be instantiated and solved simply by calling a method with the
initial value problem as argument.  For methods with error estimates,
adaptive time-stepping can be used based on a specified error tolerance.
NodePy also includes automated functions for convergence testing.

In the future, NodePy may also support solving semi-discretizations
of initial boundary value PDEs.


NodePy Manual
==================================

.. toctree::
   :maxdepth: 2

   quickstart
   methods
   solving
   stability
   rooted_trees
   changes
   future
   about


Modules Reference
===================

.. toctree::

    modules/runge_kutta_method
    modules/rooted_trees
    modules/linear_multistep_method
    modules/low_storage_rk

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. bibliography:: zrefs.bib
   :all:

=======================
About NodePy
=======================

NodePy is an open-source, free project.  It is primarily developed by David
Ketcheson.  Several other people have made important contributions to NodePy,
including Robert Bradshaw, Matteo Parsani, Umair bin Waheed, Hendrik Ranocha,
and Yiannis Hadjimichael.

Citing
=======================

If you use NodePy in work that is published, please cite

  David I. Ketcheson, NodePy Software <version number>

Contributing
=======================

Contributions to the package are most welcome.  If you have 
used NodePy for research, chances are that others would find your
code useful.  Feel free to either e-mail a patch or fork the
repository and 
`issue a pull request on Github <https://github.com/ketch/nodepy/compare>`_.
If you find a bug or there's an improvement you'd like to see, you can
also `open an issue <https://github.com/ketch/nodepy/issues>`_.

License
=======================
NodePy is distributed under the terms of the modified Berkeley Software Distribution
(BSD) license.  The license is in the file nodepy/LICENSE.txt and
reprinted below.

See http://www.opensource.org/licenses/bsd-license.php for more details.

Copyright (c) 2008-2020 David I. Ketcheson.  All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright 
    notice, this list of conditions and the following disclaimer in the 
    documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Funding
==========

NodePy development has been supported by:

  * A U.S. Dept. of Energy Computational Science Graduate Fellowship
  * Grants from King Abdullah University of Science & Technology
============
Rooted Trees
============
.. plot::

   from nodepy.rooted_trees import *
   tree=RootedTree('{T^2{T{T}}{T}}')
   tree.plot()


.. autoclass:: nodepy.rooted_trees.RootedTree
   :noindex:

Plotting trees
==============

A single tree can be plotted using the plot method of the RootedTree class.
For convenience, the method plot_all_trees() plots the whole forest of rooted
trees of a given order.

.. automethod:: nodepy.rooted_trees.RootedTree.plot
   :noindex:

.. plot::
   :include-source:

   from nodepy.rooted_trees import *
   tree=RootedTree('{T^2{T{T}}{T}}')
   tree.plot()

.. plot::
   :include-source:

   from nodepy.rooted_trees import *
   plot_all_trees(5)


Functions on rooted trees
===========================
.. automethod:: nodepy.rooted_trees.RootedTree.order
   :noindex:

.. automethod:: nodepy.rooted_trees.RootedTree.density
   :noindex:

.. automethod:: nodepy.rooted_trees.RootedTree.symmetry
   :noindex:

Computing products on trees
===========================
.. automethod:: nodepy.rooted_trees.RootedTree.Gprod
   :noindex:

.. automethod:: nodepy.rooted_trees.RootedTree.lamda
   :noindex:

Using rooted trees to generate order conditions
===============================================
:mod:`runge_kutta_method`
===========================

.. automodule:: nodepy.runge_kutta_method
   :members:
   :undoc-members:
:mod:`linear_multistep_method`
==============================

.. automodule:: nodepy.linear_multistep_method
   :members:
   :undoc-members:
:mod:`low_storage_rk`
===========================

.. automodule:: nodepy.low_storage_rk
   :members:
   :undoc-members:
:mod:`rooted_trees`
===========================

.. automodule:: nodepy.rooted_trees
   :members:
   :undoc-members:
