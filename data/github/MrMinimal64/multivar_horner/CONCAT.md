---
title: 'multivar_horner: A Python package for computing Horner factorisations of multivariate polynomials'
tags:
    - python
    - mathematics
    - polynomial
    - evaluation
    - multivariate
    - horner
    - factorisation
    - factorization


authors:
    - name: Jannik Michelfeit
      orcid: 0000-0002-1819-6975
      affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
    - name: Technische Universität Dresden, Dresden, Germany
      index: 1
    - name: Max Planck Institute of Molecular Cell Biology and Genetics, Dresden, Germany
      index: 2
date: 20 April 2020
bibliography: paper.bib


---

# Abstract

Many applications in the sciences require numerically stable and computationally efficient evaluation of multivariate polynomials.
Finding beneficial representations of polynomials, such as Horner factorisations, is therefore crucial.
`multivar_horner` [@github], the Python package presented here, is, as far as we are aware, the first open-source software for computing multivariate Horner factorisations.
This paper briefly outlines the functionality of the package and places it in context with respect to previous work in the field.
Benchmarks additionally demonstrate the advantages of the implementation and Horner factorisations in general.


# Introduction

Polynomials are a central concept in mathematics and find application in a wide range of fields [@prasolov2009polynomials; @boumova2002applications; @cools2002advances; @akritas1989elements; @Hecht1].
(Multivariate) polynomials have different possible mathematical representations and the beneficial properties of some representations are in great demand in many applications [@LeeFactorization2013; @leiserson2010efficient; @Hecht1].

The *Horner factorisation* is such a representation with beneficial properties.
Compared to the unfactorised representation of a multivariate polynomial, in the following called *canonical form*, this representation offers some important advantages.
Firstly, the Horner factorisation is more compact, in the sense that it requires fewer mathematical operations in order to evaluate the polynomial (cf. \autoref{fig:num_ops_growth}).
Consequently, evaluating a multivariate polynomial in Horner factorisation is faster and numerically more stable [@pena2000multivariate; @pena2000multivariate2; @greedyHorner] (cf. \autoref{fig:num_err_growth}).
These advantages come at the cost of an initial computational effort required to find the factorisation.

The `multivar_horner` Python package implements a multivariate Horner scheme ("Horner's method", "Horner's rule") [@horner1819xxi] to compute Horner factorisations of multivariate polynomials given in canonical form.
The package offers the functionality of representing multivariate polynomials of arbitrary degree in Horner factorisation as well as in canonical form.
Additionally it allows to compute the partial derivatives of a polynomial and to evaluate a polynomial at a given point.
Accordingly the package presented here is useful whenever (multivariate) polynomials have to be evaluated efficiently, the numerical error of the polynomial evaluation has to be small, or a compact representation of the polynomial is required.
This holds true for many applications applying numerical analysis.
One example use case where this package is already being employed are novel response surface methods [@michelfeitresponse] based on multivariate Newton interpolation [@Hecht1].


# Functionality

`multivar_horner` implements a multivariate Horner scheme using the greedy heuristic presented in [@greedyHorner].
In the following the key functionality of this package is outlined.
For more details on polynomials and Horner factorisations please refer to the literature, e.g. [@neumaier2001introduction].

A polynomial in canonical form is a sum of monomials.
For a univariate polynomial $f(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_d x^d$ (canonical form), the Horner factorisation is unique: $f(x) = a_0 + x ( a_1 + x( \cdots x (a_d) \cdots )$.
In the multivariate case, however, the factorisation is ambiguous, as there are multiple possible factors to factorise with.
The key functionality of `multivar_horner` is finding a good instance among the many possible Horner factorisations of a multivariate polynomial.

Let's consider the example multivariate polynomial given in canonical form by $p(x) = 5 + x_1^3 x_2 + 2 x_1^2 x_3 + 3 x_1 x_2 x_3$.
The polynomial $p$ is the sum of $4$ monomials, has dimensionality $3$ and can also be written as $p(x) = 5 x_1^0 x_2^0 x_3^0 + 1 x_1^3 x_2^1 x_3^0 + 2 x_1^2 x_2^0 x_3^1 + 3 x_1^1 x_2^1 x_3^1$.
The coefficients of the monomials are $5$, $1$, $2$ and $3$ respectively.

From this formulation it is straightforward to represent a multivariate polynomial with a single vector of coefficients and one exponent matrix.
Due to its simplicity and universality this kind of representation is used for defining polynomials as input.
It should be noted that this most trivial representation is computationally expensive to evaluate.

The number of additions of a polynomial remains constant irrespective of the polynomial factorisation, since it depends solely on the number of monomials and a factorisation does not influence the number of monomials.
This holds true only without taking common subexpression elimination into account.
Hence the number of additions is irrelevant for evaluating the quality of a factorisation.
In the following we accordingly only count the number of multiplications for a less biased comparison to other polynomial representations.
Note that each exponentiation is counted as (exponent - 1) operations.

# Usage

The following code snippet shows how to use `multivar_horner` to compute a Horner factorisation of $p$:


```python
from multivar_horner import HornerMultivarPolynomial
coefficients = [5.0, 1.0, 2.0, 3.0]
exponents = [[0, 0, 0], [3, 1, 0], [2, 0, 1], [1, 1, 1]]
p = HornerMultivarPolynomial(coefficients, exponents, rectify_input=True,
	compute_representation=True)
````

The factorisation computed by `multivar_horner` is $p(x) =  x_1 (x_1 (x_1 (1 x_2) + 2 x_3) + 3 x_2 x_3) + 5$ and requires 7 multiplications for every polynomial evaluation.
The human readable representation of the polynomial can be accessed with:

```python
print(p.representation)
# [#ops=7] p(x) = x_1 (x_1 (x_1 (1.0 x_2) + 2.0 x_3) + 3.0 x_2 x_3) + 5.0
````

It should be noted that the implemented factorisation procedure is coefficient-agnostic and hence does not, for example, optimise multiplications with $1$.
This design choice has been made in order to have the ability to change the coefficients of a computed polynomial representation a posteriori.

With the default settings a Horner factorisation is computed by recursively factorising with respect to the factor most commonly used in all monomials.
When no leaves of the resulting binary "Horner factorisation tree" can be factorised any more, a "recipe" for evaluating the polynomial is  compiled.
This recipe encodes all operations required to evaluate the polynomial in `numpy` arrays [@numpy].
This enables the use of functions just-in-time compiled by `numba` [@numba], which allow the polynomial evaluation to be computationally efficient.
The just-in-time compiled functions are always used, since a pure-Python polynomial evaluation would to some extent outweigh the benefits of Horner factorisation representations.

`multivar_horner` allows to evaluate the polynomial $p$ at a point $x$:

```python
x = [-2.0, 3.0, 1.0]
p_x = p.eval(x, rectify_input=True) # -29.0
````


# Degrees of multivariate polynomials



![The number of coefficients of fully occupied polynomials of different degrees in 3 dimensions.\label{fig:num_coeff_growth}](num_coeff_growth.png)

In contrast to the one-dimensional case, there are several concepts of degree for polynomials in multiple dimensions.
Following the notation of [@trefethen2017multivariate] the usual notion of degree of a polynomial, the *total degree*, is the maximal sum of exponents of all monomials.
This is equal to the maximal $l_1$-norm of all exponent vectors of the monomials.
Accordingly the *euclidean degree* is the maximal $l_2$-norm and the *maximal degree* is the maximal $l_{\infty}$-norm of all exponent vectors.
Please refer to [@trefethen2017multivariate] for precise definitions.

A polynomial is called *fully occupied* with respect to a certain degree if all possible monomials having a smaller or equal degree are present.
The occupancy of a polynomial can then be defined as the number of monomials present, relative to the fully-occupied polynomial of this degree.
A fully-occupied polynomial hence has an occupancy of $1$.


The number of coefficients (equal to the number of possible monomials) in multiple dimensions highly depends on the type of degree a polynomial has (cf. \autoref{fig:num_coeff_growth}).
This effect intensifies as the dimensionality grows.


# Benchmarks

![Numerical error of evaluating randomly-generated polynomials of varying sizes.\label{fig:num_err_growth}](../docs/_static/num_err_growth.png)


For benchmarking our method the following procedure is used:
In order to sample polynomials with uniformly random occupancy, the probability of monomials being present is picked randomly.
For a fixed *maximal* degree $n$ in $m$ dimensions there are $(n+1)^m$ possible exponent vectors corresponding to monomials.
Each of these monomials is included with the chosen probability.

Five polynomials were sampled randomly for each maximal degree up to 7 and dimensionality up to 7.
In order to compute the numerical error, each polynomial is evaluated at a point sampled uniformly from $[-1; 1]^m$ with the different methods.
The polynomial evaluation algorithms use 64-bit floating point numbers, whereas the ground truth is computed with 128-bit accuracy in order to avoid numerical errors in the ground truth value.
To obtain more representative results, the numerical error is averaged over 100 runs with uniformly-random coefficients each in the range $[-1; 1]$.
All errors are displayed relative to the ground truth.

![Numerical error of evaluating randomly generated polynomials in canonical form relative to the Horner factorisation.\label{fig:num_err_heatmap}](../docs/_static/num_err_heatmap.png)

Note that even though the original monomials are not actually present in a Horner factorisation, the number of coefficients is nonetheless identical to the number of coefficients of its canonical form.
With increasing size in terms of the number of included coefficients, the numerical error of both the canonical form and the Horner factorisation found by `multivar_horner` grow exponentially (cf. \autoref{fig:num_err_growth}).
However, in comparison to the canonical form, the Horner factorisation is more numerically stable, as visualised in \autoref{fig:num_err_heatmap}.
The numerical stability of Horner factorisations has theoretically been shown in [@greedyHorner; @pena2000multivariate; @pena2000multivariate2].

Even though the number of operations required for evaluating the polynomials grows exponentially with their size, irrespective of the considered representation, \autoref{fig:num_ops_growth} shows that the rate of growth is lower for the Horner factorisation.
As a result, the Horner factorisations are computationally easier to evaluate.

![Number of operations required to evaluate randomly generated polynomials.\label{fig:num_ops_growth}](../docs/_static/num_ops_growth.png)

# Related work

This package has been created due to the recent advances in multivariate polynomial interpolation [@Hecht1; @Hecht2].
High-dimensional interpolants of large degrees create the demand for evaluating multivariate polynomials in a computationally efficient and numerically stable way.
These advances enable modeling the behaviour of (physical) systems with polynomials.
Obtaining an analytical, multidimensional and nonlinear representation of a system opens up many possibilities.
With so-called "interpolation response surface methods" [@michelfeitresponse], for example, a system can be analysed and optimised.


The commercial software [`Maple`](https://www.maplesoft.com/support/help/Maple/view.aspx?path=convert%2Fhorner) offers the ability to compute multivariate Horner factorisations. However `multivar_horner` is, as far as we are aware, the first open-source implementation of a multivariate Horner scheme.
The `Wolfram Mathematica` framework supports [univariate Horner factorisations](https://reference.wolfram.com/language/ref/HornerForm.html).
The `Julia` package [StaticPolynomials](https://github.com/JuliaAlgebra/StaticPolynomials.jl) has a functionality similar to `multivar_horner`, but does not support computing Horner factorisations.


[`NumPy`](https://numpy.org/doc/stable/reference/routines.polynomials.polynomial.html) [@numpy] offers functionality to represent and manipulate polynomials of dimensionality up to 3.
`SymPy` offers the dedicated module [`sympy.polys`](https://docs.sympy.org/latest/modules/polys/index.html) for symbolically operating with polynomials.
As stated in the [documentation](https://mattpap.github.io/masters-thesis/html/src/algorithms.html#evaluation-of-polynomials), `SymPy` does not use Horner factorisations for polynomial evaluation in the multivariate case.
[`Sage`](https://doc.sagemath.org/html/en/reference/polynomial_rings/index.html) covers the algebraic side of polynomials.

`multivar_horner` has no functions to directly interoperate with other software packages.
The generality of the required input parameters (coefficients and exponents), however, still ensures the compatibility with other approaches.
It is, for example, easy to manipulate a polynomial with other libraries and then compute the Horner factorisation representation of the resulting output polynomial with `multivar_horner` afterwards, by simply transferring coefficients and exponents.
Some intermediate operations to convert the parameters into the required format might be necessary.


# Further reading

The documentation of the package is hosted on [readthedocs.io](https://multivar_horner.readthedocs.io/en/latest/).
Any bugs or feature requests can be filed on [GitHub](https://github.com/MrMinimal64/multivar_horner/issues) [@github].
The [contribution guidelines](https://github.com/MrMinimal64/multivar_horner/blob/master/CONTRIBUTING.rst) can be found there as well.

The underlying basic mathematical concepts are explained in numerical analysis textbooks like [@neumaier2001introduction].
The Horner scheme at the core of `multivar_horner` has been theoretically outlined in [@greedyHorner].

Instead of using a heuristic to choose the next factor, one could instead search over all possible Horner factorisations in order to arrive at a minimal factorisation.
The number of possible factorisations, however, increases exponentially with the degree and dimensionality of a polynomial (number of monomials).
One possibility to avoid computing each factorisation is to employ a version of A-star search [@hart1968formal] adapted for factorisation trees.
`multivar_horner` also implements this approach, which is similar to the branch-and-bound method suggested in [@kojima2008efficient, ch. 3.1].

[@carnicer1990evaluation] shows how factorisation trees can be used to evaluate multivariate polynomials and their derivatives.
In [@kuipers2013improving] Monte Carlo tree search has been used to find more performant factorisations than with greedy heuristics.
Other beneficial representations of polynomials are specified, for example, in [@LeeFactorization2013] and [@leiserson2010efficient].


# Acknowledgements

Thanks to Michael Hecht (Max Planck Institute of Molecular Cell Biology and Genetics) and Steve Schmerler (Helmholtz-Zentrum Dresden-Rossendorf) for valuable input enabling this publication.
I also thank the editor David P. Sanders (Universidad Nacional Autónoma de México) as well as the reviewers Henrik Barthels (RWTH Aachen University) and Sascha Timme (TU Berlin) for their helpful feedback.


# References
=======================
Contribution Guidelines
=======================

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs via `Github Issues`_.

If you are reporting a bug, please include:

* Your version of this package, python and Numba (if you use it)
* Any other details about your local setup that might be helpful in troubleshooting, e.g. operating system.
* Detailed steps to reproduce the bug.
* Detailed description of the bug (error log etc.).


Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "help wanted"
and not assigned to anyone is open to whoever wants to implement it - please
leave a comment to say you have started working on it, and open a pull request
as soon as you have something working, so that Travis starts building it.

Issues without "help wanted" generally already have some code ready in the
background (maybe it's not yet open source), but you can still contribute to
them by saying how you'd find the fix useful, linking to known prior art, or
other such help.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Probably for some features the documentation is missing or unclear. You can help with that!


Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue via `Github Issues`_.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement. Create multiple issues if necessary.
* Remember that this is a volunteer-driven project, and that contributions  are welcome :)


Get Started!
------------

Ready to contribute? Here's how to set up this package for local development.

*  Fork this repo on GitHub.
*  Clone your fork locally

* To make changes, create a branch for local development:

   .. code-block:: sh

       $ git checkout -b name-of-your-bugfix-or-feature



* Check out the instructions and notes in ``publish.py``
* Install ``tox`` and run the tests:

   .. code-block:: sh

       $ pip install tox
       $ tox

   The ``tox.ini`` file defines a large number of test environments, for
   different Python etc., plus for checking codestyle. During
   development of a feature/fix, you'll probably want to run just one plus the
   relevant codestyle:

   .. code-block:: sh

       $ tox -e codestyle


* Commit your changes and push your branch to GitHub:

   .. code-block:: sh

       $ git add .
       $ git commit -m "Your detailed description of your changes."
       $ git push origin name-of-your-bugfix-or-feature

* Submit a pull request through the GitHub website. This will trigger the Travis CI build which runs the tests against all supported versions of Python.



.. _Github Issues: https://github.com/MrMinimal64/multivar_horner/issues
Changelog
=========


TODOs

* integrate CD PyPI upload into GHA workflow. including authentication...
* run speed and numerical tests with the new C evaluation method!
* Improve tests
* compare poly.num_ops of different factorisations. tests?
* num_ops currently will be 0 when caching is used (no factorisation will be computed)


POSSIBLE IMPROVEMENTS:

MultivarPoly:

- also make use of the concept of 'recipes' for efficiently evaluating the polynomial
    skipping the most unnecessary operations
- add option to skip this optimisation

HornerMultivarPoly:

- optimise factor evaluation (save instructions, 'factor factorisation'):
    a monomial factor consists of scalar factors and in turn some monomial factors consist of other monomial factors

-> the result of evaluating a factor can be reused for evaluating other factors containing it

-> find the optimal 'factorisation' of the factors themselves

-> set the factorisation_idxs of each factor in total_degree to link the evaluation appropriately

idea:
    choose  'Goedel IDs' as the monomial factor ids
    then the id of a monomial is the product of the ids of its scalar factors
    find the highest possible divisor among all factor ids
    (corresponds to the 'largest' factor included in the monomial)
    this leads to a minimal factorisation for evaluating the monomial values quickly

- add option to skip this optimisation to save build time

- optimise space requirement:
    after building a factorisation tree for the factors themselves,
    then use its structure to cleverly reuse storage space
    -> use compiler construction theory: minimal assembler register assignment, 'graph coloring'...

- optimise 'copy recipe': avoid copy operations for accessing values of x
    problem: inserting x into the value array causes operations as well and
        complicates address assigment and recipe compilation

-  when the polynomial does not depend on all variables, build a wrapper to maintain the same "interface"
    but internally reduce the dimensionality, this reduced the size of the numpy arrays -> speed, storage benefit

- the evaluation of subtrees is independent and could theoretically be done in parallel
    probably not worth the effort. more reasonable to just evaluate multiple polynomials in parallel


3.0.0 (2021-12-04)
__________________

ATTENTION: major changes:

* introduced the default behavior of compiling the evaluation instructions in C code (C compiler required)
* the previous ``numpy+numba`` evaluation using "recipes" is the fallback option in case the C file could not be compiled
* as a consequence dropping ``numba`` as a required dependency
* added the "extra" ``numba`` to install on demand with: ``pip install multivar_horner[numba]``
* introduced custom polynomial hashing and comparison operations
* using hash to cache and reuse the instructions for evaluation (for both C and recipe instructions)
* introduced constructions argument ``store_c_instr`` (``HornerMultivarPolynomial``) to force the storage of evaluation code in C for later usage
* introduced constructions argument ``store_numpy_recipe`` (``HornerMultivarPolynomial``) to force the storage of the custom "recipe" data structure required for the evaluation using ``numpy`` and ``numba``
* introduced class ``HornerMultivarPolynomialOpt`` for optimal Horner Factorisations to separate code and simplify tests
* as a consequence dropped construction argument ``find_optimal`` of class ``HornerMultivarPolynomial``
* introduced constructions argument ``verbose`` to show the output of status print statements
* dropping official python3.6 support because ``numba`` did so (supporting Python3.7+)

internal:

* using poetry for dependency management
* using GitHub Actions for CI instead of travis


2.2.0 (2021-02-04)
__________________

ATTENTION: API changes:

* removed ``validate_input`` arguments. input will now always be validated (otherwise the numba jit compiled functions will fail with cryptic error messages)
* black code style
* pre-commit checks


2.1.1 (2020-10-01)
__________________

Post-JOSS paper review release:

* Changed the method of counting the amount of operations of the polynomial representations. Only the multiplications are being counted. Exponentiations count as (exponent-1) operations.
* the numerical tests compute the relative average error with an increased precision now


2.1.0 (2020-06-15)
__________________


ATTENTION: API changes:

* ``TypeError`` and ``ValueError`` are being raised instead of ``AssertionError`` in case of invalid input parameters with ``validate_input=True``
* added same parameters and behavior of ``rectify_input`` and ``validate_input`` in the ``.eval()`` function of polynomials


internal:

* Use ``np.asarray()`` instead of ``np.array()`` to avoid unnecessary copies
* more test cases for invalid input parameters



2.0.0 (2020-04-28)
__________________

* BUGFIX: factor evaluation optimisation caused errors in rare cases. this optimisation has been removed completely. every factor occurring in a factorisation is being evaluated independently now. this simplifies the factorisation process. the concept of "Goedel ID" (=unique encoding using prime numbers) is not required any more
* ATTENTION: changed polynomial degree class attribute names to comply with naming conventions of the scientific literature
* added __call__ method for evaluating a polynomial in a simplified notation ``v=p(x)``
* fixed dependencies to: ``numpy>=1.16``, ``numba>=0.48``
* clarified docstrings (using Google style)
* more verbose error messages during input verification
* split up ``requirements.txt`` (into basic dependencies and test dependencies)
* added sphinx documentation
* updated benchmark results

tests:

* added test for numerical stability
* added plotting features for evaluating the numerical stability
* added tests comparing functionality to 1D ``numpy`` polynomials
* added tests comparing functionality to naive polynomial evaluation
* added basic API functionality test

internal:

* added class ``AbstractPolynomial``
* added typing
* adjusted publishing routine
* testing multiple python versions
* using the specific tags of the supported python version for the build wheels
* removed ``example.py``


1.3.0 (2020-03-14)
__________________


* NEW FEATURE: changing coefficients on the fly with ``poly.change_coefficients(coeffs)``
* NEW DEPENDENCY: ``python3.6+`` (for using f'' format strings)
* the real valued coefficients are now included in the string representation of a factorised polynomial
* add contribution guidelines
* added instructions in readme, ``example.py``
* restructured the factorisation routine (simplified, clean up)
* extended tests


1.2.0 (2019-05-19)
__________________

* support of newer numpy versions (ndarray.max() not supported)
* added plotting routine (partly taken from tests)
* added plots in readme
* included latest insights into readme


1.1.0 (2019-02-27)
__________________

* added option `find_optimal` to find an optimal factorisation with A* search, explanation in readme
* optimized heuristic factorisation (more clean approach using just binary trees)
* dropped option `univariate_factors`
* added option `compute_representation` to compute the string representation of a factorisation only when required
* added option `keep_tree` to keep the factorisation tree when required
* clarification and expansion of readme and `example.py`
* explained usage of optional parameters `rectify_input=True` and `validate_input=True`
* explained usage of functions `get_gradient()` and `get_partial_derivative(i)`
* averaged runtime in speed tests


1.0.1 (2018-11-12)
__________________

* introducing option to only factor out single variables with the highest usage with the optional parameter ``univariate_factors=True``
* compute the number of operations needed by the horner factorisation by the length of its recipe (instead of traversing the full tree)
* instead of computing the value of scalar factors with exponent 1, just copy the values from the given x vector ("copy recipe")
* compile the initial value array at construction time


1.0.0 (2018-11-08)
__________________

* first stable release


0.0.1 (2018-10-05)
__________________

* birth of this package
===============
multivar_horner
===============


.. image:: https://travis-ci.org/jannikmi/multivar_horner.svg?branch=master
    :alt: CI status
    :target: https://travis-ci.org/jannikmi/multivar_horner

.. image:: https://readthedocs.org/projects/multivar_horner/badge/?version=latest
    :alt: documentation status
    :target: https://multivar_horner.readthedocs.io/en/latest/?badge=latest

.. image:: https://img.shields.io/pypi/wheel/multivar-horner.svg
    :target: https://pypi.python.org/pypi/multivar-horner

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit

.. image:: https://pepy.tech/badge/multivar-horner
    :alt: Total PyPI downloads
    :target: https://pepy.tech/project/multivar-horner

.. image:: https://img.shields.io/pypi/v/multivar_horner.svg
    :alt: latest version on PyPI
    :target: https://pypi.python.org/pypi/multivar-horner

.. image:: https://joss.theoj.org/papers/0b514c6894780f3cc81ed88c141631d4/status.svg
    :alt: JOSS status
    :target: https://joss.theoj.org/papers/0b514c6894780f3cc81ed88c141631d4

.. image:: https://zenodo.org/badge/155578190.svg
   :target: https://zenodo.org/badge/latestdoi/155578190

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black


``multivar_horner`` is a python package implementing a multivariate
`Horner scheme ("Horner's method", "Horner's rule") <https://en.wikipedia.org/wiki/Horner%27s_method>`__
for efficiently evaluating multivariate polynomials.


Quick Guide:

::


    pip install multivar_horner


For efficiency this package is compiling the instructions required for polynomial evaluation to C by default.
If you don't have a C compiler (``gcc`` or ``cc``) installed you also need to install numba for using an alternative method:

::


    pip install multivar_horner[numba]


Simple example:

.. code-block:: python

    import numpy as np
    from multivar_horner.multivar_horner import HornerMultivarPolynomial

    # input parameters defining the polynomial
    #   p(x) = 5.0 + 1.0 x_1^3 x_2^1 + 2.0 x_1^2 x_3^1 + 3.0 x_1^1 x_2^1 x_3^1
    coefficients = np.array([[5.0], [1.0], [2.0], [3.0]], dtype=np.float64)
    exponents = np.array([[0, 0, 0], [3, 1, 0], [2, 0, 1], [1, 1, 1]], dtype=np.uint32)

    # [#ops=7] p(x) = x_1 (x_1 (x_1 (1.0 x_2) + 2.0 x_3) + 3.0 x_2 x_3) + 5.0
    horner_polynomial = HornerMultivarPolynomial(coefficients, exponents)
    x = np.array([-2.0, 3.0, 1.0], dtype=np.float64)
    p_x = horner_polynomial(x)


For more refer to the `documentation <https://multivar_horner.readthedocs.io/en/latest/>`__.


Also see:
`GitHub <https://github.com/jannikmi/multivar_horner>`__,
`PyPI <https://pypi.python.org/pypi/multivar_horner/>`__,
`arXiv paper <https://arxiv.org/abs/2007.13152>`__
.. _contributing:

.. include:: ../CONTRIBUTING.rst
===============
Getting started
===============


Installation
------------

Installation with pip:

::

    pip install multivar_horner



For efficiency this package is compiling the instructions required for polynomial evaluation to C by default.
If you don't have a C compiler (``gcc`` or ``cc``) installed you also need to install numba for using an alternative method:

::


    pip install multivar_horner[numba]



Basics
------

Let's consider this example multivariate polynomial:

:math:`p(x) = 5 + 1 x_1^3 x_2^1 + 2 x_1^2 x_3^1 + 3 x_1^1 x_2^1 x_3^1`


Which can also be written as:

:math:`p(x) = 5 x_1^0 x_2^0 x_3^0 + 1 x_1^3 x_2^1 x_3^0 + 2 x_1^2 x_2^0 x_3^1 + 3 x_1^1 x_2^1 x_3^1`


A polynomial is a sum of monomials.
Our example polynomial has :math:`M = 4` monomials and dimensionality :math:`N = 3`.

The coefficients of our example polynomial are: 5.0, 1.0, 2.0, 3.0

The exponent vectors of the corresponding monomials are:

* [0, 0, 0]
* [3, 1, 0]
* [2, 0, 1]
* [1, 1, 1]

To represent polynomials this package requires the coefficients and the exponent vectors as input.

This code shows how to compute the Horner factorisation of our example polynomial :math:`p`
and evaluating :math:`p` at a point :math:`x`:

.. code-block:: python

    import numpy as np
    from multivar_horner.multivar_horner import HornerMultivarPolynomial

    coefficients = np.array([[5.0], [1.0], [2.0], [3.0]], dtype=np.float64)  # shape: (M,1)
    exponents = np.array(
        [[0, 0, 0], [3, 1, 0], [2, 0, 1], [1, 1, 1]], dtype=np.uint32
    )  # shape: (M,N)
    p = HornerMultivarPolynomial(coefficients, exponents)

    x = np.array([-2.0, 3.0, 1.0], dtype=np.float64)  # shape: (1,N)
    p_x = p(x)  # -29.0



.. note::

    with the default settings the input is required to have these data types and shapes


With the class ``HornerMultivarPolynomial`` a polynomial can be represented in :ref:`Horner factorisation <horner_usage>`.

With the class ``HornerMultivarPolynomialOpt`` a polynomial can be represented in an :ref:`optimal Horner factorisation <optimal_usage>`.

With the class ``MultivarPolynomial`` a polynomial can be represented in :ref:`canonical form <canonical_usage>`.


All available features of this package are explained :ref:`HERE <usage>`.

The API documentation can be found :ref:`HERE <api>`.
.. _api:

=================
API documentation
=================


.. py:module:: multivar_horner


HornerMultivarPolynomial
------------------------
.. autoclass:: HornerMultivarPolynomial


MultivarPolynomial
------------------
.. autoclass:: MultivarPolynomial


HornerMultivarPolynomialOpt
---------------------------

.. autoclass:: HornerMultivarPolynomialOpt
.. include:: ../CHANGELOG.rst
===============
References
===============

.. bibliography:: refs.bib
   :style: unsrt

.. image:: https://travis-ci.org/jannikmi/multivar_horner.svg?branch=master
    :target: https://travis-ci.org/jannikmi/multivar_horner

.. image:: https://readthedocs.org/projects/multivar_horner/badge/?version=latest
    :alt: documentation status
    :target: https://multivar_horner.readthedocs.io/en/latest/?badge=latest

.. image:: https://img.shields.io/pypi/wheel/multivar-horner.svg
    :target: https://pypi.python.org/pypi/multivar-horner

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit

.. image:: https://pepy.tech/badge/multivar-horner
    :alt: Total PyPI downloads
    :target: https://pepy.tech/project/multivar-horner

.. image:: https://img.shields.io/pypi/v/multivar_horner.svg
    :alt: latest version on PyPI
    :target: https://pypi.python.org/pypi/multivar-horner

.. image:: https://joss.theoj.org/papers/0b514c6894780f3cc81ed88c141631d4/status.svg
    :alt: JOSS status
    :target: https://joss.theoj.org/papers/0b514c6894780f3cc81ed88c141631d4

.. image:: https://zenodo.org/badge/155578190.svg
   :target: https://zenodo.org/badge/latestdoi/155578190

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
.. _optimal:

=============================
Optimal Horner Factorisations
=============================

See :ref:`this code <optimal_usage>` for an example usage.


Instead of using a heuristic to choose the next factor one can allow a search over all possible (meaningful) factorisations in order to arrive at a minimal Horner factorisation.
The amount of possible factorisations however is increasing exponentially with the degree of a polynomial and its amount of monomials.
One possibility to avoid computing each factorisation is to employ a version of A*-search :cite:`hart1968formal` adapted for factorisation trees:
• Initialise a set of all meaningful possible first level Newton factorisations
• Rank all factorisation according to a lower bound (“heuristic”) of their lowest possible amount of operations
• Iteratively factorise the most promising factorisation and update the heuristic
• Stop when the most promising factorisation is fully factorised

This approach is guaranteed to yield a minimal Horner factorisation, but its performance highly depends on the heuristic in use: Irrelevant factorisations are only being ignored if the heuristic is not too optimistic in estimating the amount of operations. On the other hand the heuristic must be easy to compute, because it would otherwise be computationally cheaper to just try all different factorisations.
Even though it missing to cover exponentiations, the branch-and-bound method suggested in :cite:`kojima2008efficient` (ch. 3.1) is almost identical to this procedure.


Even with a good heuristic this method is only traceable for small polynomials because of its increased resource requirements.
Since experiments show that factorisations obtained by choosing one factorisation according to a heuristic have the same or only a slightly higher amount of included operations:cite:`kojima2008efficient` (ch. 7), the computational effort of this approach is not justifiable in most cases.
A use case however is to compute and store a minimal representation of a polynomial in advance if possible.

**NOTES:**

* for the small polynomial examples in the current tests, the results were identical (in terms of #ops) with the approach of just using the default heuristic = trying one factorisation (further analysis needed)!
* in some cases this approach currently is trying all possible factorisations, because the heuristic in use is too optimistic (= brute force, further analysis and improvements needed)
* this requires MUCH more computational resources than just trying one factorisation (the number of possible factorisations is growing exponentially with the size of the polynomial!).
* there are possibly many optimal Horner factorisations of a multivariate polynomial. one could easily adapt this approach to find all optimal Horner factorisations
* even an optimal Horner factorisation must not be the globally minimal representation of a polynomial. there are possibly better types of factorisations and techniques: e.g. "algebraic factorisation", "common subexpression elimination"
* there are multiple possible concepts of optimality (or minimality) of a polynomial


============================
multivar_horner
============================

.. include:: ./badges.rst



a python package implementing a multivariate `Horner scheme ("Horner's method", "Horner's rule") <https://en.wikipedia.org/wiki/Horner%27s_method>`__:cite:`horner1819xxi` for efficiently evaluating multivariate polynomials.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

    Getting Started <0_getting_started>
    Usage <1_usage>
    About <2_about>
    Optimal Horner Factorisations <3_optimal_factorisation>
    API  <4_api>
    Contributing <5_contributing>
    Changelog <6_changelog>
    References <7_references>

.. reference file must be loaded last for the bibtex citations to work (-> alphanumeric "last" name)!



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
=====
About
=====


.. include:: ./badges.rst

``multivar_horner`` is a python package implementing a multivariate `Horner scheme ("Horner's method", "Horner's rule") <https://en.wikipedia.org/wiki/Horner%27s_method>`__:cite:`horner1819xxi`  for efficiently evaluating multivariate polynomials.

For an explanation of multivariate Horner factorisations and the terminology used here refer to e.g. `Greedy algorithms for optimizing multivariate Horner schemes <https://dl.acm.org/doi/pdf/10.1145/980175.980179>`__ :cite:`ceberio2004greedy`

A given input polynomial in canonical form (or normal form) is being factorised according to the greedy heuristic described in :cite:`ceberio2004greedy` with some additional computational tweaks.
The resulting Horner factorisation requires less operations for evaluation and is being computed by growing a "Horner Factorisation Tree".
When the polynomial is fully factorized (= all leaves cannot be factorised any more), a computational "recipe" for evaluating the polynomial is being compiled.
This "recipe" (stored internally as ``numpy`` arrays) enables computationally efficient evaluation.
``Numba`` just in time compiled functions operating on the ``numpy`` arrays make this fast.
All factors appearing in the factorisation are being evaluated only once (reusing computed values).

**Pros:**

 * computationally efficient representation of a multivariate polynomial in the sense of space and time complexity of the evaluation
 * less roundoff errors :cite:`pena2000multivariate,pena2000multivariate2`
 * lower error propagation, because of fewer operations :cite:`ceberio2004greedy`


**Cons:**

 * increased initial computational requirements and memory to find and then store the factorisation


The impact of computing Horner factorisations has been evaluated in the :ref:`benchmarks below <benchmarks>`.

With this package it is also possible to represent polynomials in :ref:`canonical form <canonical_usage>` and to search for an :ref:`optimal Horner factorisation <optimal_usage>`.


Also see:
`GitHub <https://github.com/jannikmi/multivar_horner>`__,
`PyPI <https://pypi.python.org/pypi/multivar_horner/>`__,
`arXiv paper <https://arxiv.org/abs/2007.13152>`__



Dependencies
------------

``python>=3.6``,
``numpy>=1.16``,
``numba>=0.48``



License
-------

``multivar_horner`` is distributed under the terms of the MIT license
(see `LICENSE <https://github.com/jannikmi/multivar_horner/blob/master/LICENSE>`__).



.. _benchmarks:

Benchmarks
----------

To obtain meaningful results the benchmarks presented here use polynomials sampled randomly with the following procedure:
In order to draw polynomials with uniformly random occupancy, the probability of monomials being present is picked randomly.
For a fixed maximal degree n in m dimensions there are (n+1)^m possible exponent vectors corresponding to monomials.
Each of these monomials is being activated with the chosen probability.


Refer to :cite:`trefethen2017multivariate` for an exact definition of the maximal degree.

For each maximal degree up to 7 and until dimensionality 7, 5 polynomials were drawn randomly.
Note that even though the original monomials are not actually present in a Horner factorisation, the amount of coefficients however is identical to the amount of coefficients of its canonical form.


Even though the amount of operations required for evaluating the polynomials grow exponentially with their size irrespective of the representation, the rate of growth is lower for the Horner factorisation.


.. figure:: _static/num_ops_growth.png

    amount of operations required to evaluate randomly generated polynomials.



Due to this, the bigger the polynomial the more compact the Horner factorisation representation is relative to the canonical form.
As a result the Horner factorisations are computationally easier to evaluate.


Numerical error
^^^^^^^^^^^^^^^

In order to compute the numerical error, each polynomial has been evaluated at a point chosen uniformly random from $[-1; 1]^m$ with the different methods.
The polynomial evaluation algorithms use 64-bit floating point numbers, whereas the ground truth has been computed with 128-bit accuracy in order to avoid numerical errors in the ground truth value.
To receive more representative results, the obtained numerical error is being averaged over 100 tries with uniformly random coefficients each in the range $[-1; 1]$,
All errors are displayed as (averaged) absolute values.


With increasing size in terms of the amount of included coefficients the numerical error of both the canonical form and the Horner factorisation found by ``multivar_horner`` grow exponentially.


.. figure:: _static/num_err_growth.png

    numerical error of evaluating randomly generated polynomials of varying sizes.


In comparison to the canonical form however the Horner factorisation is much more numerically stable.
This has also been visualised in the following figure:


.. figure:: _static/num_err_heatmap.png

    numerical error of evaluating randomly generated polynomials in canonical form relative to the Horner factorisation.


.. note::

    if you require an even higher numerical stability you can set ``FLOAT_DTYPE = numpy.float128``
    or ``FLOAT_DTYPE = numpy.longfloat`` in ``global_settings.py``.
    Then however the jit compilation has to be removed in ``helper_fcts_numba.py`` (``Numba`` does not support float128).




Speed tests
^^^^^^^^^^^

The following speed benchmarks have been performed on a 2017 MacBook Pro: 4x2,8 GHz Intel Core i7 CPU, 16 GB 2133 MHz LPDDR3 RAM, macOS 10.13 High Sierra.
The software versions in use were: ``multivar_horner 2.0.0``, ``python 3.8.2``, ``numpy 1.18.1`` and ``numba 0.48.0``
Both evaluation algorithms (with and without Horner factorisation) make use of ``Numba`` just in time compiled functions.



::

    Speed test:
    testing 100 evenly distributed random polynomials
    average timings per polynomial:

     parameters   |  setup time (s)                         |  eval time (s)                       |  # operations                        | lucrative after
    dim | max_deg | canonical  | horner     | delta         | canonical  | horner     | delta      | canonical  | horner     | delta      |    # evals
    ================================================================================================================================================================
    1   | 1       | 4.90e-05   | 2.33e-04   | 3.8 x more    | 8.96e-06   | 1.28e-05   | 0.4 x more | 3          | 1          | 2.0 x less | -
    1   | 2       | 5.24e-05   | 1.95e-04   | 2.7 x more    | 3.42e-06   | 6.01e-06   | 0.8 x more | 4          | 2          | 1.0 x less | -
    1   | 3       | 5.07e-05   | 2.31e-04   | 3.6 x more    | 3.48e-06   | 5.86e-06   | 0.7 x more | 6          | 3          | 1.0 x less | -
    1   | 4       | 5.04e-05   | 2.65e-04   | 4.3 x more    | 3.59e-06   | 5.62e-06   | 0.6 x more | 7          | 4          | 0.8 x less | -
    1   | 5       | 5.08e-05   | 3.04e-04   | 5.0 x more    | 3.49e-06   | 8.47e-06   | 1.4 x more | 8          | 6          | 0.3 x less | -
    1   | 6       | 4.81e-05   | 4.65e-04   | 8.7 x more    | 3.54e-06   | 6.72e-06   | 0.9 x more | 10         | 7          | 0.4 x less | -
    1   | 7       | 5.39e-05   | 4.00e-04   | 6.4 x more    | 3.95e-06   | 6.49e-06   | 0.6 x more | 12         | 8          | 0.5 x less | -
    1   | 8       | 5.19e-05   | 3.83e-04   | 6.4 x more    | 5.63e-06   | 6.16e-06   | 0.1 x more | 12         | 8          | 0.5 x less | -
    1   | 9       | 4.88e-05   | 4.42e-04   | 8.0 x more    | 3.73e-06   | 6.05e-06   | 0.6 x more | 14         | 10         | 0.4 x less | -
    1   | 10      | 4.89e-05   | 5.41e-04   | 10 x more     | 3.80e-06   | 7.11e-06   | 0.9 x more | 15         | 10         | 0.5 x less | -

    2   | 1       | 8.34e-05   | 3.11e-04   | 2.7 x more    | 3.85e-06   | 6.09e-06   | 0.6 x more | 11         | 3          | 2.7 x less | -
    2   | 2       | 4.96e-05   | 7.05e-04   | 13 x more     | 3.80e-06   | 5.82e-06   | 0.5 x more | 26         | 10         | 1.6 x less | -
    2   | 3       | 5.20e-05   | 9.75e-04   | 18 x more     | 4.50e-06   | 6.70e-06   | 0.5 x more | 38         | 16         | 1.4 x less | -
    2   | 4       | 5.93e-05   | 1.44e-03   | 23 x more     | 5.53e-06   | 7.12e-06   | 0.3 x more | 63         | 27         | 1.3 x less | -
    2   | 5       | 5.26e-05   | 2.25e-03   | 42 x more     | 6.49e-06   | 6.46e-06   | -0.0 x more | 91         | 39         | 1.3 x less | 59828
    2   | 6       | 5.31e-05   | 2.90e-03   | 54 x more     | 7.65e-06   | 6.55e-06   | 0.2 x less | 127        | 54         | 1.4 x less | 2595
    2   | 7       | 5.72e-05   | 3.76e-03   | 65 x more     | 9.02e-06   | 6.03e-06   | 0.5 x less | 164        | 70         | 1.3 x less | 1238
    2   | 8       | 5.32e-05   | 4.39e-03   | 81 x more     | 9.71e-06   | 6.06e-06   | 0.6 x less | 198        | 84         | 1.4 x less | 1186
    2   | 9       | 5.27e-05   | 5.04e-03   | 95 x more     | 1.08e-05   | 7.25e-06   | 0.5 x less | 230        | 99         | 1.3 x less | 1418
    2   | 10      | 5.47e-05   | 6.74e-03   | 122 x more    | 1.36e-05   | 6.46e-06   | 1.1 x less | 310        | 132        | 1.3 x less | 935

    3   | 1       | 4.96e-05   | 5.69e-04   | 10 x more     | 3.70e-06   | 6.18e-06   | 0.7 x more | 26         | 7          | 2.7 x less | -
    3   | 2       | 5.34e-05   | 2.02e-03   | 37 x more     | 5.43e-06   | 6.70e-06   | 0.2 x more | 97         | 28         | 2.5 x less | -
    3   | 3       | 5.42e-05   | 4.47e-03   | 82 x more     | 8.88e-06   | 6.13e-06   | 0.4 x less | 222        | 68         | 2.3 x less | 1605
    3   | 4       | 5.59e-05   | 8.40e-03   | 149 x more    | 1.44e-05   | 6.92e-06   | 1.1 x less | 434        | 133        | 2.3 x less | 1115
    3   | 5       | 5.73e-05   | 1.35e-02   | 236 x more    | 2.10e-05   | 1.36e-05   | 0.5 x less | 685        | 211        | 2.2 x less | 1809
    3   | 6       | 7.70e-05   | 2.32e-02   | 300 x more    | 3.72e-05   | 8.75e-06   | 3.3 x less | 1159       | 355        | 2.3 x less | 811
    3   | 7       | 6.86e-05   | 3.46e-02   | 504 x more    | 5.71e-05   | 8.90e-06   | 5.4 x less | 1787       | 543        | 2.3 x less | 717
    3   | 8       | 7.07e-05   | 4.64e-02   | 655 x more    | 6.97e-05   | 9.97e-06   | 6.0 x less | 2402       | 730        | 2.3 x less | 775
    3   | 9       | 8.34e-05   | 6.90e-02   | 826 x more    | 1.05e-04   | 1.15e-05   | 8.2 x less | 3613       | 1084       | 2.3 x less | 736
    3   | 10      | 9.21e-05   | 9.54e-02   | 1034 x more   | 1.42e-04   | 1.35e-05   | 9.5 x less | 4988       | 1485       | 2.4 x less | 742

    4   | 1       | 5.45e-05   | 1.25e-03   | 22 x more     | 4.94e-06   | 6.49e-06   | 0.3 x more | 67         | 14         | 3.8 x less | -
    4   | 2       | 5.83e-05   | 7.20e-03   | 122 x more    | 1.19e-05   | 7.65e-06   | 0.6 x less | 390        | 91         | 3.3 x less | 1673
    4   | 3       | 6.57e-05   | 2.35e-02   | 357 x more    | 3.39e-05   | 7.93e-06   | 3.3 x less | 1295       | 303        | 3.3 x less | 903
    4   | 4       | 7.22e-05   | 4.96e-02   | 686 x more    | 6.68e-05   | 1.02e-05   | 5.6 x less | 2753       | 653        | 3.2 x less | 874
    4   | 5       | 9.85e-05   | 1.17e-01   | 1186 x more   | 1.56e-04   | 1.74e-05   | 8.0 x less | 6588       | 1535       | 3.3 x less | 843
    4   | 6       | 1.40e-04   | 1.98e-01   | 1416 x more   | 2.66e-04   | 1.96e-05   | 13 x less  | 11036      | 2582       | 3.3 x less | 802
    4   | 7       | 1.77e-04   | 3.27e-01   | 1846 x more   | 4.29e-04   | 2.93e-05   | 14 x less  | 18271      | 4276       | 3.3 x less | 820
    4   | 8       | 2.77e-04   | 5.97e-01   | 2153 x more   | 8.33e-04   | 4.72e-05   | 17 x less  | 33518      | 7736       | 3.3 x less | 760
    4   | 9       | 3.82e-04   | 8.90e-01   | 2330 x more   | 1.16e-03   | 6.35e-05   | 17 x less  | 47086      | 10944      | 3.3 x less | 812
    4   | 10      | 5.44e-04   | 1.30e+00   | 2388 x more   | 1.80e-03   | 8.80e-05   | 20 x less  | 73109      | 16873      | 3.3 x less | 758




Related work
------------

This package has been created due to the recent advances in multivariate polynomial interpolation :cite:`Hecht1,Hecht2`.
High dimensional interpolants of large degrees create the demand for evaluating multivariate polynomials computationally efficient and numerically stable.

:cite:`carnicer1990evaluation` shows how factorisation trees can be used to evaluate multivariate polynomials and their derivatives.

In :cite:`kuipers2013improving` Monte Carlo tree search has been used to find more performant factorisations than with greedy heuristics.

Other representations of polynomials are being presented, among others, in :cite:`LeeFactorization2013,leiserson2010efficient`.




Contact
--------


Tell me if and how your are using this package. This encourages me to develop and test it further.

Most certainly there is stuff I missed, things I could have optimized even further or explained more clearly, etc.
I would be really glad to get some feedback.

If you encounter any bugs, have suggestions etc. do not hesitate to **open an Issue** or **add a Pull Requests** on Git.
Please refer to the :ref:`contribution guidelines <contributing>`


Acknowledgements
----------------

Thanks to:

`Steve <https://github.com/elcorto>`__ for valuable feedback and writing tests.
.. _usage:

=====
Usage
=====

.. note::

   Also check out the :ref:`API documentation <api>` or the `code <https://github.com/MrMinimal64/multivar_horner>`__.


Let's look at the example multivariate polynomial:

:math:`p(x) = 5 + 1 x_1^3 x_2^1 + 2 x_1^2 x_3^1 + 3 x_1^1 x_2^1 x_3^1`


Which can also be written as:

:math:`p(x) = 5 x_1^0 x_2^0 x_3^0 + 1 x_1^3 x_2^1 x_3^0 + 2 x_1^2 x_2^0 x_3^1 + 3 x_1^1 x_2^1 x_3^1`

A polynomial is a sum of monomials.
Our example polynomial has :math:`M = 4` monomials and dimensionality :math:`N = 3`.

The coefficients of our example polynomial are: 5.0, 1.0, 2.0, 3.0

The exponent vectors of the corresponding monomials are:

* [0, 0, 0]
* [3, 1, 0]
* [2, 0, 1]
* [1, 1, 1]

To represent polynomials this package requires the coefficients and the exponent vectors as input.


.. code-block:: python

    import numpy as np

    coefficients = np.array(
        [[5.0], [1.0], [2.0], [3.0]], dtype=np.float64
    )  # numpy (M,1) ndarray
    exponents = np.array(
        [[0, 0, 0], [3, 1, 0], [2, 0, 1], [1, 1, 1]], dtype=np.uint32
    )  # numpy (M,N) ndarray


.. note::

    by default the Numba jit compiled functions require these data types and shapes



.. _horner_usage:

Horner factorisation
---------------------


to create a representation of the multivariate polynomial :math:`p` in Horner factorisation:

.. code-block:: python

    from multivar_horner.multivar_horner import HornerMultivarPolynomial

    horner_polynomial = HornerMultivarPolynomial(coefficients, exponents)


the found factorisation is :math:`p(x) = x_1^1 (x_1^1 (x_1^1 (1.0 x_2^1) + 2.0 x_3^1) + 3.0 x_2^1 x_3^1) + 5.0`.


pass ``rectify_input=True`` to automatically try converting the input to the required ``numpy`` data structures and types


.. code-block:: python


    coefficients = [5.0, 1.0, 2.0, 3.0]
    exponents = [[0, 0, 0], [3, 1, 0], [2, 0, 1], [1, 1, 1]]
    horner_polynomial = HornerMultivarPolynomial(
        coefficients, exponents, rectify_input=True
    )



pass ``keep_tree=True`` during construction of a Horner factorised polynomial,
when its factorisation tree should be kept after the factorisation process
(e.g. to be able to compute string representations of the polynomials later on)


.. code-block:: python

    horner_polynomial = HornerMultivarPolynomial(coefficients, exponents, keep_tree=True)



.. note::

    for increased efficiency the default for both options is ``False``


.. _canonical_usage:

canonical form
--------------

it is possible to represent the polynomial without any factorisation (refered to as 'canonical form' or 'normal form'):

.. code-block:: python

    from multivar_horner.multivar_horner import MultivarPolynomial

    polynomial = MultivarPolynomial(coefficients, exponents)


use this if ...

* the Horner factorisation takes too long
* the polynomial is going to be evaluated only a few times
* fast polynomial evaluation is not required or
* the numerical stability of the evaluation is not important


.. note::

    in the case of unfactorised polynomials many unnecessary operations are being done
    (internally uses naive numpy matrix operations)




string representation
---------------------


in order to compile a string representation of a polynomial pass ``compute_representation=True`` during construction

.. note::

    the number in square brackets indicates the number of multiplications required
    to evaluate the polynomial.

.. note::

    exponentiations are counted as exponent - 1 operations, e.g. x^3 <-> 2 operations

.. code-block:: python

    polynomial = MultivarPolynomial(coefficients, exponents)
    print(polynomial)  # [#ops=10] p(x)


    polynomial = MultivarPolynomial(coefficients, exponents, compute_representation=True)
    print(polynomial)
    # [#ops=10] p(x) = 5.0 x_1^0 x_2^0 x_3^0 + 1.0 x_1^3 x_2^1 x_3^0 + 2.0 x_1^2 x_2^0 x_3^1 + 3.0 x_1^1 x_2^1 x_3^1

    horner_polynomial = HornerMultivarPolynomial(
        coefficients, exponents, compute_representation=True
    )
    print(horner_polynomial.representation)
    # [#ops=7] p(x) = x_1 (x_1 (x_1 (1.0 x_2) + 2.0 x_3) + 3.0 x_2 x_3) + 5.0


the formatting of the string representation can be changed with the parameters ``coeff_fmt_str`` and ``factor_fmt_str``:

.. code-block:: python

    polynomial = MultivarPolynomial(
        coefficients,
        exponents,
        compute_representation=True,
        coeff_fmt_str="{:1.1e}",
        factor_fmt_str="(x{dim} ** {exp})",
    )


the string representation can be computed after construction as well.


.. note::

    for ``HornerMultivarPolynomial``: ``keep_tree=True`` is required at construction time


.. code-block:: python

    polynomial.compute_string_representation(
        coeff_fmt_str="{:1.1e}", factor_fmt_str="(x{dim} ** {exp})"
    )
    print(polynomial)
    # [#ops=10] p(x) = 5.0e+00 (x1 ** 0) (x2 ** 0) (x3 ** 0) + 1.0e+00 (x1 ** 3) (x2 ** 1) (x3 ** 0)
    #                   + 2.0e+00 (x1 ** 2) (x2 ** 0) (x3 ** 1) + 3.0e+00 (x1 ** 1) (x2 ** 1) (x3 ** 1)



change the coefficients of a polynomial
---------------------------------------

in order to access the polynomial string representation with the updated coefficients pass ``compute_representation=True``
with ``in_place=False`` a new polygon object is being generated


.. note::

    the string representation of a polynomial in Horner factorisation depends on the factorisation tree.
    the polynomial object must hence have keep_tree=True


.. code-block:: python

    new_coefficients = [
        7.0,
        2.0,
        0.5,
        0.75,
    ]  # must not be a ndarray, but the length must still fit
    new_polynomial = horner_polynomial.change_coefficients(
        new_coefficients,
        rectify_input=True,
        compute_representation=True,
        in_place=False,
    )



.. _optimal_usage:

optimal Horner factorisations
-----------------------------


use the class ``HornerMultivarPolynomialOpt`` for the construction of the polynomial
to trigger an adapted A* search to find the optimal factorisation.

See :ref:`this chapter <optimal>` for further information.


.. note::

    time and memory consumption is MUCH higher!

.. code-block:: python

    from multivar_horner import HornerMultivarPolynomialOpt

    horner_polynomial_optimal = HornerMultivarPolynomialOpt(
        coefficients,
        exponents,
        compute_representation=True,
        rectify_input=True,
    )




Caching
-------------------

by default the instructions required for evaluating a Horner factorised polynomial will be cached either as ``.c`` file or ``.pickle`` file in the case of ``numpy+numba`` evaluation.

One can explicitly force the compilation of the instructions in the required format:

.. code-block:: python

    horner_polynomial = HornerMultivarPolynomial(
        coefficients, exponents, store_c_instr=True, store_numpy_recipe=True
    )


If you construct a Horner polynomial with the same properties (= exponents) these cached instructions will be used for evaluation and a factorisation won't be computed again.
Note that as a consequence you won't be able to access the factorisation tree and string representation in these cases.

the cached files are being stored in ``<path/to/env/>multivar_horner/multivar_horner/__pychache__/``

.. code-block:: python

   horner_polynomial.c_file
   horner_polynomial.c_file_compiled
   horner_polynomial.recipe_file


you can read the content of the cached C instructions:

.. code-block:: python

   instr = horner_polynomial.get_c_instructions()
   print(instr)


you can also export the whole polynomial class (including the string representation etc.):

.. code-block:: python

    path = "file_name.pickle"
    polynomial.export_pickle(path=path)


to load again:

.. code-block:: python

    from multivar_horner import load_pickle

    polynomial = load_pickle(path)



evaluating a polynomial
-----------------------

in order to evaluate a polynomial at a point ``x``:


.. code-block:: python

    # define a query point and evaluate the polynomial
    x = np.array([-2.0, 3.0, 1.0], dtype=np.float64)  # numpy ndarray with shape [N]
    p_x = polynomial(x)  # -29.0


or


.. code-block:: python

    p_x = polynomial.eval(x)  # -29.0


or

.. code-block:: python

    x = [-2.0, 3.0, 1.0]
    p_x = polynomial.eval(x, rectify_input=True)  # -29.0


As during construction of a polynomial instance, pass ``rectify_input=True`` to automatically try converting the input to the required ``numpy`` data structure.


.. note::

    the default for both options is ``False`` for increased speed

.. note::

    the dtypes are fixed due to the just in time compiled ``Numba`` functions


computing the partial derivative of a polynomial
------------------------------------------------


.. note::

    BETA: untested feature


.. note::

    partial derivatives will be instances of the same parent class



.. note::

    all given additional arguments will be passed to the constructor of the derivative polynomial


.. note::

    dimension counting starts with 1 -> the first dimension is #1!


.. code-block:: python

    deriv_2 = polynomial.get_partial_derivative(2, compute_representation=True)
    # p(x) = x_1 (x_1^2 (1.0) + 3.0 x_3)




computing the gradient of a polynomial
------------------------------------------------

.. note::

    BETA: untested feature



.. note::

    all given additional arguments will be passed to the constructor of the derivative polynomials



.. code-block:: python

    grad = polynomial.get_gradient(compute_representation=True)
    # grad = [
    #     p(x) = x_1 (x_1 (3.0 x_2) + 4.0 x_3) + 3.0 x_2 x_3,
    #     p(x) = x_1 (x_1^2 (1.0) + 3.0 x_3),
    #     p(x) = x_1 (x_1 (2.0) + 3.0 x_2)
    # ]
