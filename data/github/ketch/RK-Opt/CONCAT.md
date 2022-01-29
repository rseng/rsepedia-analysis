---
title: '`RK-Opt`: A package for the design of numerical ODE solvers'
tags:
  - Python
  - numerical analysis
  - differential equations
  - Runge-Kutta method
  - linear multistep method
  - strong stability preservation
  - absolute stability
authors:
  - name: David I. Ketcheson^[Corresponding author.]
    orcid: 0000-0002-1212-126X
    affiliation: 1
  - name: Matteo Parsani
    orcid: 0000-0001-7300-1280
    affiliation: 1
  - name: Zachary J. Grant
    orcid: 0000-0002-1293-4770
    affiliation: 2
  - name: Aron J. Ahmadia
    orcid: 0000-0002-2573-2481
    affiliation: 3
  - name: Hendrik Ranocha
    orcid: 0000-0002-3456-2277
    affiliation: 1
affiliations:
 - name: King Abdullah University of Science and Technology, Saudi Arabia
   index: 1
 - name: Oak Ridge National Laboratory, USA
   index: 2
 - name: Capital One, USA
   index: 3
date: 9 July 2020
bibliography: paper.bib
---

# Summary

Ordinary and partial differential equations (ODEs and PDEs) are used to model
many important phenomena.  In most cases, solutions of these models must be
approximated by numerical methods.  Most of the relevant algorithms fall within
a few classes of methods, with the properties of individual methods determined
by their coefficients.  The choice of appropriate coefficients in the design of
methods for specific applications is an important area of research.
`RK-Opt` is a software package for designing numerical ODE solvers with
coefficients optimally chosen to provide desired properties.
It is available from https://github.com/ketch/RK-Opt, with documentation
at http://numerics.kaust.edu.sa/RK-Opt/.
The primary focus of the package is on the design of Runge-Kutta methods, but
some routines for designing other classes of methods such as multistep
Runge-Kutta and general linear methods are also included.

# Statement of need

Over the last several decades, a great deal of work has gone into the design
of numerical ODE solvers.  Initially this work was aimed at developing general
purpose solvers, but over time the emphasis shifted increasingly toward
development of optimized methods for specific applications.  Different
accuracy, stability, performance, and other properties may be relevant or
essential depending on the nature of the equations to be solved.

An s-stage Runge-Kutta method has roughly s^2 coefficients (roughly s^2/2 for
explicit methods), which can be chosen so as to provide high accuracy,
stability, or other properties. Historically, most interest in Runge-Kutta
methods has focused on methods using the minimum number of stages for a given
order of accuracy. However, in the past few decades there has been increasing
recognition that using extra stages can be worthwhile in order to improve other
method properties. Some areas where this is particularly useful are in the
enhancement of linear and nonlinear stability properties, the reduction of
storage requirements, and the design of embedded pairs. Methods with dozens or
even hundreds of stages are not unheard of.

At the same time, most existing Runge-Kutta methods have been designed by hand,
by researchers laboriously solving the order conditions. When using extra
stages, the number of available parameters makes the selection of a
near-optimal choice by hand impossible, and one resorts to computational
optimization. This leads to a different paradigm of numerical method design, in
which we use sophisticated numerical (optimization) algorithms to design
sophisticated numerical (integration) algorithms. It can be expected that this
trend will accelerate in the future, and perhaps one day simple
manually-constructed algorithms will be the exception.

RK-Opt contains a set of tools for designing Runge-Kutta methods in this paradigm.
It provides code that can enforce desired properties and/or objective
functions.  The constraints and objective are then used within an optimization
framework, to determine coefficients of methods that best achieve the desired
goal.  Thus, `RK-Opt` is a sort of meta-software, consisting of algorithms whose
purpose is to create other algorithms.

Typically, the most obvious formulation of the corresponding
optimization problem is intractable.  Therefore, these problems
are reformulated in ways that make them amenable to available techniques.
These reformulations include, for instance, turning a nonconvex problem into
a sequence of convex problems or even linear programs.  The resulting algorithms
can often guarantee optimality of their output.  However, for the
general problem of determining Runge-Kutta coefficients, the nonconvex problem
must be attacked directly and optimality cannot be guaranteed.

`RK-Opt` is written entirely in MATLAB, and leverages the MATLAB Optimization
Toolbox as well as the Global Optimization Toolbox.
Its development has been motivated largely by research needs and
it has been used in a number of papers (see below).

# Features

`RK-Opt` includes the following subpackages.

## `polyopt`

This package computes optimal stability functions for Runge-Kutta methods.
Here *optimal* means that the stable step size is maximized for a given ODE
spectrum.  The corresponding optimization problem is intractable under a
direct implementation.  The package uses the algorithm developed in
[@2012_optimal_stability_polynomials], which relaxes the global optimization
problem by solving a sequence of convex subproblems.  Under certain technical
assumptions, the result is guaranteed to be the optimal solution of the
original problem.  `polyopt` relies on CVX [@cvx;@gb08] to solve the convex
subproblems.  This package is usually used as the first step in designing a
Runge-Kutta method.

## `RK-Coeff-Opt`

This package computes optimal Runge-Kutta coefficients based on a desired
set of constraints and an objective.  Available constraints include:

 - The number of stages and order of accuracy
 - The class of method (explicit, implicit, diagonally implicit, low-storage)
 - The coefficients of the stability polynomial (usually determined using `polyopt`)

Two objective functions are provided; methods can be optimized for the
strong stability preserving (SSP) coefficient or the principal error norm
(a measure of the leading-order truncation error coefficients).
In addition to standard Runge-Kutta methods, various classes of multistep
Runge-Kutta methods can also be optimized.

The optimization problem in question is highly nonconvex and the available
solvers may fail to find a solution, or may converge to a non-optimal solution.
For this reason, the implementation is based on solving many local optimization
problems in parallel from different random initial points, using MATLAB's Global
Optimization Toolbox.

The packages `dwrk-opt` and `low-storage` are specialized but less full-featured
versions of `RK-Coeff-Opt` that were developed for specific research projects
involving downwind Runge-Kutta methods and low-storage Runge-Kutta methods, respectively.

## `am_radius-opt`

Whereas the previous two subpackages are fairly general-purpose tools,
this package solves a very specific and discrete set of problems described in
[@2009_monotonicity].  Specifically, the provided routines determine the coefficients of
multistep methods (including classes of general linear methods) with the
largest possible SSP coefficient (also known
as radius of absolute monotonicity).  The corresponding optimization problem
had previously been attacked using brute force search, but this limited
its solvability to methods with very few steps.  In this package the
problem is reformulated as a sequence of linear programming problems,
enabling its efficient solution for methods with many steps.


# Related research and software

`RK-Opt` development has proceeded in close connection to the `NodePy` package (https://github.com/ketch/NodePy).
Whereas `RK-Opt` is focused on the design of numerical methods, `NodePy` is focused
more on their analysis.  A common workflow involves generating new methods with
`RK-Opt` and then studying their properties in more detail using `NodePy`.

Some of the research projects that have made use of `RK-Opt` include development of:

 - SSP Runge-Kutta methods
   [@2008_explicit_ssp;@2009_implicit_ssp;@gottlieb2015optimal;@higueras2019strong]
 - SSP linear multistep methods [@2009_monotonicity]
 - SSP general linear methods [@2011_tsrk;@2017_msrk]
 - SSP IMEX Runge-Kutta methods [@conde2017implicit]
 - Low-storage Runge-Kutta methods [@2010_LSRK;@higueras2019new]
 - Optimal Runge-Kutta stability polynomials [@2012_optimal_stability_polynomials]
 - Additive and downwind SSP Runge-Kutta methods [@2011_dwssp;@2018_perturbations]
 - Optimal Runge-Kutta methods for specific PDE semi-discretizations [@parsani-eccomas;@Parsani_finnish;@2013_sd_erk;@2014_ssp_rkdg]
 - Optimal Runge-Kutta methods for pseudo-time stepping [@vermeire2019optimal;@vermeire2020optimal]
 - Embedded pairs for Runge-Kutta methods [@conde2018embedded]
 - Runge-Kutta methods with high weak stage order [@2018_wso]
 - SSP multistage, multiderivative methods [@christlieb2016explicit;@grant2019strong;@reynoso2017strong]

As can be seen from this list, applications have mostly stemmed from the
work of the main developer's research group, but have since expanded
beyond that.

Because of the nature of `RK-Opt`, applications often involve writing some additional
code to impose special constraints, or simply using the existing code as a template.
A number of related optimization routines written for similar purposes in this
vein can be found at https://github.com/SSPmethods.

# Acknowledgements

Much of the initial `RK-Opt` development was performed by D. Ketcheson while
he was supported by a DOE Computational Science Graduate Fellowship and by
AFOSR grant number FA9550-06-1-0255.  Development
has also been supported by funding from King Abdullah University of Science and Technology.

# References

# RK-Opt: A Package for the Design of Numerical ODE solvers

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rk-opt.readthedocs.io/en/latest/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-success.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4138076.svg)](https://doi.org/10.5281/zenodo.4138076)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02514/status.svg)](https://doi.org/10.21105/joss.02514)

See the full documentation [here](https://rk-opt.readthedocs.io/en/latest/).

RK-Opt is a collection of MATLAB codes for designing optimized numerical ODE solvers.
The main emphasis is on Runge-Kutta methods, but some routines deal with other classes of methods.
It is primarily developed and used by the
[KAUST Numerical Mathematics Group](http://numerics.kaust.edu.sa).
It includes the following sub-packages:

 - **polyopt**: Find optimal stability polynomials of a given degree and order of
   accuracy for a specified spectrum.
 - **RK-coeff-opt**: Find optimal Runge-Kutta method coefficients, for a prescribed
   order of accuracy and number of stages.
 - **am_rad-opt**: Find stability functions with optimal radius of absolute monotonicity.
   Includes capabilities for both multistep and multistage methods.
 - **RKtools**: A collection of routines for analyzing or computing various
   properties of Runge-Kutta methods.  For a much more extensive package along these
   lines, see [NodePy](http://nodepy.readthedocs.io/en/latest/).

A common workflow for designing Runge-Kutta methods is to use **polyopt** to find an
appropriate stability function and then **RK-coeff-opt** to determine the Runge-Kutta
method coefficients.

To run the tests, execute the MATLAB script `test.m`. This requires a relatively recent
version of MATLAB (tested with R2018a and later) with the following toolboxes.
 - MATLAB Optimization Toolbox
 - MATLAB Global Optimization Toolbox
 - CVX (http://cvxr.com/cvx/)
 - MATLAB Parallel Computing Toolbox (optional; allows faster searching for optimal methods in RK-Coeff-Opt)

# Citing
If you use RK-Opt in published work, please cite the following paper:

Ketcheson et al., (2020). RK-Opt: A package for the design of numerical ODE solvers. Journal of Open Source Software, 5(54), 2514, https://doi.org/10.21105/joss.02514

You can use the following bibtex entry:

```
@article{Ketcheson2020,
  doi = {10.21105/joss.02514},
  url = {https://doi.org/10.21105/joss.02514},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {54},
  pages = {2514},
  author = {David I. Ketcheson and Matteo Parsani and Zachary J. Grant and Aron J. Ahmadia and Hendrik Ranocha},
  title = {`RK-Opt`: A package for the design of numerical ODE solvers},
  journal = {Journal of Open Source Software}
}
```

# Authors
The code is primarily developed and maintained by David Ketcheson.
The following people have also made important contributions to RK-Opt (listed alphabetically):

 - Aron Ahmadia: Co-developer of **polyopt** algorithm and routines.
 - Zachary Grant: Extension of order conditions to multistep RK with more than two stages and
    addition of order conditions for orders 9-11.
 - Matteo Parsani: Many improvements to **RK-coeff-opt** routines and organization.
 - Hendrik Ranocha: General improvements and updates, including updating the test routines.
Contributions to RK-Opt are most welcome!  If you have an idea for a new
feature, or simply see something that should be fixed/improved, you can:

  - [Raise an issue to let us know](https://github.com/ketch/RK-Opt/issues)
  - Better yet, fork the RK-Opt repository on Github, implement the change
    yourself, and [open a pull request](https://github.com/ketch/RK-Opt/pulls)

Ideally, contributions that include new functionality should also include
documentation and testing, but we can help with that part as needed.

If you're raising an issue to report a bug, please provide enough detail to
let us reproduce the bug.
This directory contains experimental research code for
investigating downwind-biased Runge-Kutta methods.
Use at your own risk!
This directory contains experimental research code related to
finding optimal low-storage Runge-Kutta methods.  Use at your
own risk!
