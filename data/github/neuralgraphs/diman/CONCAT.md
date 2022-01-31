---
title: 'diman: A Clojure Package for Dimensional Analysis'
tags:
  - Clojure
  - neuroscience
  - modeling
  - simulation
  - quantitative analysis
  - scientific computing
authors:
  - name: Lungsi Sharma^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-1607-0164
    affiliation: 1
affiliations:
 - name: Ronin Institute
   index: 1

date: 30 July 2021
bibliography: paper.bib
---

# Summary

`diman` (**dim**ensional **an**alysis) is a Clojure based scientific software package with the ability to: create dimensional formulas, create dimensional equations, check dimensional homogeneity (consistency), and derive dimensionless products.

`diman` provides functions for each step of the analytic process for checking dimensional homogeneity or deriving dimensionless products; the repetitive operations (computational) are hidden. Users can write compound functions that perform a desired process. Thus, not only is the computational labor saved, but also introspection of the analysis is possible; the analyst is able to go through the steps of dimensional analysis.

# Statement of need

Explaining the mechanism of a phenomenon is often the goal of experiments. As most mechanistic description is expressible in terms of some measurable quantity, its value is a function of other measurable quantities; the function represents the relationship among the quantities, which provides a mechanistic explanation. For example, $F = ma = m \frac{dv}{dt}$ where the measurable value of force $F$ is a function of the measurable quantities: mass, $m$; velocity, $v$; and time, $t$.

Some or all the independent variables of the parent (first or original) functions have dimensions. Since most of the functions are unknown, and hence conceptual, the researcher deals with many candidates for independent variable, whose considerations are based on experimental results. Although the mathematical expression of the function is unknown, knowledge of the relationship among the measurable quantities is profitable not only in putting together the series of experimental results to explain the mechanism, but also testing the hypothesis presented by the function.

If possible, it is beneficial to use the transformed parent function, where all the independent variables are dimensionless. Dimensionless products are magnitudes that contains information on the dimensional quantities that it is a product of. Therefore, not only are points in a graph of dimensionless products experimentally determinable, but also are more informative than dimensional graphs. Reducing the number of independent variables to a smaller collection of dimensionless products can assist in understanding the mechanism of the phenomenon [@Langhaar:1951; @Sharma:2021].

Numerous software packages have been developed to deal with dimensions in some shape or form [@Preussner:2018; @Sharma:2021]. Most incorporate the ability to tag quantities with units, however, few are capable of doing consistency checks and fewer still deal with dimensionless products let alone, deriving dimensionless products.

Frink and F sharp are two active languages that incorporate units of measurement. Frink is a calculating tool and programming language designed to make physical calculations as it tracks units of measure throughout calculations [@Eliasen:2004]. F# is a general purpose language which allows annotating floating point and integer values with statically-typed unit metadata. F#'s units of measure is based on a prototype Meta Language (ML) like language that has dimensional type; the language supports type polymorphism as well as dimension polymorphism [@Kennedy:1996].

In F# the notion of dimension is based on ideas from type theory; therefore, dimensions satisfy algebraic properties of Abelian groups (commutative groups). The concepts of dimensional type and dimensional invariance can be mathematically extended to apply Buckingham's Theorem (also called Pi theorem) [@Kennedy:1996]&mdash;this theorem is the basis for deriving the set of dimensionless products [@Ngwua:2020].

`diman` is designed with an emphasis on **analysis**; the  application of the algebraic theory of dimensionally homogeneous functions [@Langhaar:1951]. Unlike F#'s units of measure which performs dimensional checks during calculations, checks for dimensional homogeneity in diman is manual. Also, the possibility of applying Buckingham's theorem in F# is done by mathematical abstraction [@Kennedy:2010] while diman provides functions to derive the complete set of dimensionless products of a given equation. Therefore, diman is easily put into practice for providing complete information of experimental results of physical systems in a compact form and also transform hypothetical function that describes the physical system [@Sharma:2021].

# Design and implementation

Based on the International System of Units, `diman` uses the seven base (or elementary) dimensions: [M], [L], [T], [A], [K], [mol] and [cd] for the quantities: mass, length, time, electric current, thermodynamic temperature, amount of substance and luminous intensity respectively [@BIPM:2020]. They are defined in `base_dimensions`. Furthermore, some well-known dimensions derived from the `base_dimensions` are defined in `standard_formula`; a dimensional formula for respective quantity is its dimension.

## Consistency checking

This is done by the predicate `consistent?`. There are some preliminary steps before invoking the predicate. Consider the given function $E = \frac{1}{2}mv^2$.

We define the variables
```
=> (def variables [{:symbol "E", :quantity "energy"}
                   {:symbol "m", :quantity "mass"}
                   {:symbol "v", :quantity "velocity"}])
```
then the equation
```
=> (def equation {:lhs "E^(1)", :rhs "0.5*m^(1)*v^(2)"})
```
Finally, the predicate `consistent?` is used to check if the equation is dimensionally homogenous.
```
=> (consistent? variables equation)
true
```

## Derivation of set of dimensionless products

Imagine that the study of a system results in a hypothesis such that some measurable dimensionless product is a homogeneous function $f$ of the independent variables $P$, $Q$, $R$, $S$, $T$, $U$ and $V$. Also, assume that the independent variables have dimensions such that
```
=> (def dimensional_formulae_of_all_independent_variables
        [{:quantity "term-p", :dimension "[M^(2)*L^(1)]"}
         {:quantity "term-q", :dimension "[M^(-1)*T^(1)]"}
         {:quantity "term-r", :dimension "[M^(3)*L^(-1)]"}
         {:quantity "term-s", :dimension "[T^(3)]"}
         {:quantity "term-t", :dimension "[L^(2)*T^(1)]"}
         {:quantity "term-u", :dimension "[M^(-2)*L^(1)*T^(-1)]"}
         {:quantity "term-v", :dimension "[M^(1)*L^(2)*T^(2)]"}]) 
```
Supposing the independent variables of the parent function $f$ are not already defined in `standard_formula`, inject the dimensions of the independent variables into the `standard_formula` for the present read–eval–print loop session by
```
=> (update-sformula dimensional_formulae_of_all_independent_variables)
```
Thus, `diman` now contains dimensions of the independent variables of $f$. Hence, the independent variables can be defined as
```
=> (def independent_variables
        [{:symbol "P", :quantity "term-p"}
         {:symbol "Q", :quantity "term-q"}
         {:symbol "R", :quantity "term-r"}
         {:symbol "S", :quantity "term-s"}
         {:symbol "T", :quantity "term-t"}
         {:symbol "U", :quantity "term-u"}
         {:symbol "V", :quantity "term-v"}]) 
```
The theory of dimensionless products [@Ngwua:2020] tells us that the derivation of dimensionless products can be broken down into four steps: generate the dimensional matrix, solve the homogeneous equation, determine the solution matrix and get the set of dimensionless products. Compounding the first three steps into one code block we get,
```
=> (def solution_matrix
        (get-solution-matrix
            (solve (get-augmented-matrix
						(generate-dimmat independent_variables)))))
```
This is the solution matrix for a complete set of dimensionless products.
```
=> (view-matrix solution_matrix)
[1 0 0 0 -11N 5N 8N]
[0 1 0 0 9N -4N -7N]
[0 0 1 0 -9N 5N 7N]
[0 0 0 1 15N -6N -12N]
Size -> 4 x 7
```
The set of dimensionless products can be obtained from the solution matrix by using the function `get-dimensionless-products`. Thus
```
=> (println (get-dimensionless-products solution_matrix independent_variables))
  [{:symbol "pi0", :expression "P^(1)*T^(-11)*U^(5)*V^(8)"}
  {:symbol "pi1", :expression "Q^(1)*T^(9)*U^(-4)*V^(-7)"}
  {:symbol "pi2", :expression "R^(1)*T^(-9)*U^(5)*V^(7)"}
  {:symbol "pi3", :expression "S^(1)*T^(15)*U^(-6)*V^(-12)"}]
```
or
$$
\pi_0 = PT^{-11}U^5V^8, \quad \pi_1 = QT^9U^{-4}V^{-7}, \quad \pi_2 = RT^{-9}U^5V^7, \quad \pi_3 = ST^{15}U^{-6}V^{-12}
$$
Therefore, function $f$ is transformed into some function $f_1$ whose independent variables are the dimensionless products; $\pi_0$, $\pi_1$, $\pi_2$, and $\pi_3$&mdash;$\pi$ is the conventional notation for any dimensionless product. Thus, the number of variables is reduced from 7 to 4.

# Conclusion

`diman` is a Clojure library with no other dependencies. It has its own linear algebra submodule which provides all the necessary operations. Internally, the numerical data type is Clojure's *ratio*; a ratio between integers rather than floats [@Clojure:2020]. This avoids truncation and rounding errors. Since dimensional analysis does not often involve very large matrices, the hit on computational performance due to using the *ratio* number type is practically insignificant. `diman` supplies all the necessary functions for dimensional homogeneity operations and the derivation of dimensionless products; thus making the analysis steps transparent.


# Acknowledgements

The project received no funding.

# References
# Change Log
All notable changes to this project will be documented in this file. This change log follows the conventions of [keepachangelog.com](http://keepachangelog.com/).

## [Unreleased]
### Changed
- Add a new arity to `make-widget-async` to provide a different widget shape.

## [0.1.1] - 2019-04-26
### Changed
- Documentation on how to make the widgets.

### Removed
- `make-widget-sync` - we're all async, all the time.

### Fixed
- Fixed widget maker to keep working when daylight savings switches over.

## 0.1.0 - 2019-04-26
### Added
- Files from the new template.
- Widget maker public API - `make-widget-sync`.

[Unreleased]: https://github.com/your-name/diman/compare/0.1.1...HEAD
[0.1.1]: https://github.com/your-name/diman/compare/0.1.0...0.1.1
![diman logo](./resources/images/logo/diman.png)
# diman

[![Clojars Project](https://img.shields.io/clojars/v/com.neuralgraphs/diman.svg)](https://clojars.org/com.neuralgraphs/diman)

A Clojure library for applying dimensional analysis.

[Motivation for diman.](./doc/ProjectPlan.pdf)

## Current Features

- Create dimensional formulae.
- Create dimensional equations.
- Implement principle of dimensional homogeneity; Perform consistency checks.
- Derive dimensionless products.

## Usage

The easiest way to get all the built-in functions is to be in the default namespace `(in-ns 'diman.default)`. Then, `(println default-functions)` to list all the available functions. 

- Tutorial: Generate dimensional formulae and perform consistency checking; [AsciiDoc](./doc/tutorial1.adoc), [reStructuredText for Bitbucket](./doc/tutorial1.rst)
- Tutorial: Derive dimensionless products; [AsciiDoc](./doc/tutorial2.adoc), [reStructuredText for Bitbucket](./doc/tutorial2.rst)
- Example: Journal Bearing; [AsciiDoc](./doc/tutorial3.adoc), [reStructuredText for Bitbucket](./doc/tutorial3.rst)
- Rationale for the seven base dimensions; [AsciiDoc](./doc/rationale1.adoc), [reStructuredText for Bitbucket](./doc/rationale1.rst)
- Rationale for implementing the steps for deriving a complete set of dimensionless products; [AsciiDoc](./doc/rationale2.adoc), [reStructuredText for Bitbucket](./doc/rationale2.rst)
- [Source code documentation](https://cljdoc.org/d/com.neuralgraphs/diman)

### For Clojure Noobs

Since [Leiningen](https://leiningen.org/) is one of the easiest way to use Clojure, I recommend using Leiningen to run diman. Once Leiningen is installed you can use diman in two ways; by cloning this diman repo and starting up a repl ([Read-eval-print loop](https://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop)) inside the cloned directory `cd ~/diman`, and by making diman as a dependency to your clojure project.

#### 1. Running diman by cloning the repo

Once you have cloned the repository do `cd ~/path/to/diman`, then
```
lein repl
```

By default you should already be in the **default namespace**, that is, there is no need to `(in-ns 'diman.default)`. To list all the essential functions (and therefore all the functions for the tutorials) do `(println default-functions)`.

#### 2. Running diman as a dependency (recommended)

Assuming you already have a clojure project or you can create one with the command `lein new <project-name>`, then diman can be added as a dependency in the `project.clj` file by doing
```
...
:dependencies [[org.clojure/clojure "1.10.3"]
               [com.neuralgraphs/diman "x.y.z"]]
...
```

To go through the tutorials, startup a repl (`lein repl`) inside the created project (`cd /path/to/<project-name>`) load the diman libraries as follows
```
(require '[diman.dimensions :refer [base_dimensions standard_formula update-sformula]]
         '[diman.formula :refer [formula-term formula-eqn-side formula-eqn-side-manifold]]
         '[diman.analyze :refer [dimnames consistent?]]
         '[diman.buckingham [dimensional-matrix :refer [generate-dimmat]]
                            [homogeneous-equation :refer [get-augmented-matrix solve get-solution-matrix]]
                            [dimensionless-product :refer [get-dimensionless-products get-pi-expression]]]
         '[diman [core :refer [view-matrix]]]
         '[diman.linalg.matfun [rank :refer [rank]]])
```

These are all the diman libraries essential for dimensional analysis (you may copy-paste the above inside the repl).

To load specific diman libraries in specific namespace under the new project do
```
(ns <project-name>.<your-namespace>
  (:require [diman.analyze :refer [dimnames consistent?]]))
```

## Contributing to diman

Thank you for your interest in contributing to diman.
Please refer to the [guidelines](./doc/zCONTRIBUTING.adoc) on how to contribute.

## Publications

* Sharma, L., (2022). diman: A Clojure Package for Dimensional Analysis. Journal of Open Source Software, 7(69), 3735, [https://doi.org/10.21105/joss.03735](https://doi.org/10.21105/joss.03735)
    - [![DOI](https://joss.theoj.org/papers/10.21105/joss.03735/status.svg)](https://doi.org/10.21105/joss.03735)


## License

Copyright © 2021 Lungsi Ngwua

Distributed under BSD 3-Clause "New" or "Revised" License.
=========================================
Rationale for the seven base dimensions
=========================================

diman is based on the seven base dimensions for the quantities; mass, length, time, electric current, thermodynamic temperature, amount of substance and luminous intensity. Therefore, one may think of these seven names as the dimension names. But, they are quantities not dimensions. So,

What is a dimension?
====================

This is best defined by McNish as

    It is simply a tag we attach to a quantity in an equation expressing some physical law... [D]imensions are only symbols of an elementary algebra, involving neither addition nor subtraction. To ask what are the true or natural dimensions of a quantity makes no more sense than to ask what is the true or natural word for goldfish. [McNish1957]_


Dimension vs Quantity Units
---------------------------

Since dimensions are represented by symbols, they follow the rule of algebra. However, units (eg. meter is a unit for the quantity length) are not bound by the rules of an algebra. This is an important concept for dimensional analysis.

Due to this reason we can have a minimum number of dimension to create a dimensional system such that we can relate arbitrarily chosen chosen into one coherent system. But, limiting the units to certain "absolute" units from which other units are derived would be problem specific.

If one could just assign a magnitude to a single quantity and consider this as the "absolute" unit we might be able to derive a system of units (consistent around this absolute unit) to measure all other quantities by arbitrarily assigning values to several physical constants. But, in practice measuring the units for many of the quantities in this system will require difficult experiments that cannot be performed accurately. McNish says,

    The magnitudes of quantities are not determined by nature, but by the units we have arbitrarily chosen for our measurement system. [McNish1957]_


Why seven elemental dimensions is a good start?
===============================================

McNish's argument for it is

    [F]or the simplest mechanical quantities we need four dimensions to avoid ambiguities, five, if we include solid angle; that if we set one of these quantities equal to unity we can get along with four. But having one quantity equal to unity we cannot set another equal to unity without creating an ambiguity. Thus, unity itself becomes like a dimension, so again we may say we have five. Heat and electromagnetism add at least two more required dimensions. So I might venture to say that we should have seven elemental dimensions, at least, but I do not know, because I do not comprehend all of physics. One, of course, may get along with fewer dimensions if he will tolerate some ambiguities. [McNish1957]_


The seven base dimensions for diman
===================================

+---------------------------+----------------------+
| Quantity                  | Dimension (Notation) |
+===========================+======================+
| mass                      | \[M\]                |
+---------------------------+----------------------+
| length                    | \[L\]                | 
+---------------------------+----------------------+
| time                      | \[T\]                |
+---------------------------+----------------------+
| electric current          | \[A\]                |
+---------------------------+----------------------+
| thermodynamic temperature | \[K\]                |
+---------------------------+----------------------+
| amount of substance       | \[mol\]              |
+---------------------------+----------------------+
| luminous intensity        | \[cd\]               |
+---------------------------+----------------------+

The seven base dimensions for diman is based on the International System of Units, SI units [BIPM2020]_.

    A part of the secret of analysis is the art of using notation well. - Leibniz


This is implemented in diman as

::

    => (pprint base_dimensions)
    [{:quantity "mass", :dimension "[M]"}
    {:quantity "length", :dimension "[L]"}
    {:quantity "time", :dimension "[T]"}
    {:quantity "electric current", :dimension "[A]"}
    {:quantity "thermodynamic temperature", :dimension "[K]"}
    {:quantity "luminous intensity", :dimension "[cd]"}
    {:quantity "amount of substance", :dimension "[mol]"}]


Derived Dimensions
------------------

Some well known dimensions derived from the ``base_dimensions`` are placed in the ``standard_formula``.

::

    => (pprint standard_formula)
    [{:quantity "volume", :dimension "[L^(3)]"}
    {:quantity "frequency", :dimension "[T^(-1)]"}
    {:quantity "velocity", :dimension "[L^(1)*T^(-1)]"}
    {:quantity "acceleration", :dimension "[L^(1)*T^(-2)]"}
    {:quantity "force", :dimension "[M^(1)*L^(1)*T^(-2)]"}
    {:quantity "mass density", :dimension "[M^(1)*L^(-3)]"}
    {:quantity "energy", :dimension "[M^(1)*L^(2)*T^(-2)]"}
    {:quantity "work", :dimension "[M^(1)*L^(2)*T^(-2)]"}
    {:quantity "amount of heat", :dimension "[M^(1)*L^(2)*T^(-2)]"}
    {:quantity "pressure", :dimension "[M^(1)*L^(-1)*T^(-2)]"}
    {:quantity "stress", :dimension "[M^(1)*L^(-1)*T^(-2)]"}
    {:quantity "catalytic activity", :dimension "[mol^(1)*T^(-1)]"}
    {:quantity "charge", :dimension "[A^(1)*T^(1)]"}
    {:quantity "capacitance", :dimension "[M^(-1)*L^(-2)*T^(4)*A^(2)]"}
    {:quantity "inductance", :dimension "[M^(1)*L^(2)*T^(-2)*A^(-2)]"}
    {:quantity "resistance", :dimension "[M^(1)*L^(2)*T^(-3)*A^(-2)]"}
    {:quantity "conductance", :dimension "[M^(-1)*L^(-2)*T^(3)*A^(2)]"}
    {:quantity "magnetic flux density", :dimension "[M^(1)*T^(-2)*A^(-1)]"}
    {:quantity "electromotive force", :dimension "[M^(1)*L^(2)*T^(-3)*A^(-1)]"}
    {:quantity "power", :dimension "[M^(1)*L^(2)*T^(-3)]"}
    {:quantity "magnetic flux", :dimension "[M^(1)*L^(2)*T^(-2)*A^(-1)]"}]

Notice that the derived dimensions are in a sense the *dimensional formula for respective quantity*. Hence, the name ``standard_formula``.



References
==========

.. [McNish1957] McNish, A. G. (1957, April 1). Dimensions units and standards. *Physics Today*, 10(4), 19. `10.1063/1.3060330 <https://doi.org/10.1063/1.3060330>`_

.. [BIPM2020] BIPM (2020). *Base unit definitions*. Retrieved from the `Base units page <https://www.bipm.org/en/measurement-units/base-units.html>`_.
==========================================================================================
Rationale for implementing the steps for deriving a complete set of dimensionless products
==========================================================================================

The rationale for the derivation of a complete set dimensionless product is based on Buckingham's theorem [Ngwua2021]_. Langhaar provides an excellent guide on the systematic calculation [Langhaar1951]_. There are some minor difference in implementing the derivation steps in diman. Aside from the programmatic implementation for generating the dimensional matrix and returning the dimensionless products in its expression form the uniqueness in diman is in the steps for solving the homogeneous equation.

The typical approach to solving the system of homogeneous equation is followed by a substitution step which would then result in the solution matrix of the complete set of dimensionless products. In diman the substitution step is incorporated into the solving process. This is due to the combination of constructing an augmented matrix and using a modified Gaussian elimation process that returns a reduced row-echelon form of the augmented matrix. However, although this results in all the ingredients for the solution matrix the reduced row enchelon form is not the solution matrix. An additional step is required to obtain the solution matrix from the reduced row-echelon form of the augmented matrix.

Constructing the Dimensional Matrix
===================================

Create a matrix such that its entries are the exponents of the base dimensions. A row stands for a base dimension and a column stands for a quantity in |f|. Therefore, the size of the dimensional matrix is such that the number of rows is the number of base dimensions needed to derive all the quantities in |f| and the number of columns is the number of quantities in |f|, the number of independent variable for |f|.

The Dimensional Matrix and the System of Homogeneous Equations
--------------------------------------------------------------

The row entries of a dimensional matrix may be considered the coefficients of a homogeneous equation such that the unknown variables of the equation represent exponents of respective quantities, the independent variables of |f|. The resulting product of all the quantities with their respective exponents (the unknown variables of the homogeneous equation) is a dimensionless product.

Since each row of a dimensional matrix may be representative of a homogeneous equation, it follows that a dimensional matrix will represent a system of homogeneous equation.

Solving the Homogeneous Equation
================================

In diman the first step towards solving the system of homogeneous equations is to construct the augmented matrix.

Constructing the Augmented Matrix
---------------------------------

The term augmented matrix to most people will invoke the matrix |A_b| for the linear system |Ax_b|. Since our system is a system of homogeneous equation our augmented matrix is not the same as |A_b|. It is an augmentation of the dimensional matrix such that blocks of the matrix is rearranged.

How are the two blocks rearranged?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A arbitrary dimensional matrix |m| times |n| is partitioned into two blocks which are then shifted from tail-block to head-block and head-block to tail-block. Therefore, the dimensional and the augmented matrix have the same size.

The basic idea is that the augmented matrix has at its head section (initial columns of the matrix) the m by m block taken from the tail of the dimensional matrix. The remaining part of the augmented matrix is the rest of the head of the dimensional matrix multiplied by the scalar -1. Thus, the signs for every entry of the m by m block taken from the tail remain unchanged while the sign of the entries of the remaining block of the dimensional matrix is reversed when they become part of the augmented matrix.

Below is the pseudocode for constructing the augmented matrix.

::

    START
        D <- dimensional matrix
        m <- number of rows of D
        n <- number of columns of D
        r <- rank of D
        p <- n - r
        Dsub_LHS <- D(rows 1 to m, columns 1 to p) size is m x p
        Dsub_RHS <- D(rows 1 to m, columns p+1 to n) size is m x r
        Augmented matrix <- [Dsub_RHS | (-1)*Dsub_LHS] size is m x (p + r)
    STOP

Using Modified Gaussian Elimination Method to Solve the System of Homogeneous Equations
---------------------------------------------------------------------------------------

A modified Gaussian elimation method is used such that the |m| times |m| block at the head of the augmented matrix becomes an identity matrix. The diagonal of the identity block matrix represents the left hand side of the solution for the last |m| |k_i|'s (there are |n| |k_i|'s). In other words, the returned matrix is the **reduced row-echelon form** of the augmented matrix. But, this matrix is not the solution matrix, the solution matrix for complete set of dimensionless products.

Getting the Solution Matrix
===========================

The solution matrix for a complete set of dimensionless products is determined from the solution of the augmented matrix such that the number of rows of the solution matrix is the number of dimensionless products in the set and the number of columns is the number of quantities in |f|, the number of independent variable for |f|. Below shows the pseudocode for getting the solution matrix.

::

    START
        AugS <- solved augmented matrix
        A <- strip zero rows from AugS
        m <- number of rows of A
        n <- number of columns of A
        Asub <- A(rows 1 to m, columns m+1 to n) size is m x (n - m)
        Asub_tr <- Asub transpose size is (n - m) x m
        I <- identity matrix of size (n - m) x (n - m)
        Solution matrix <- [ I | Asub_tr ] size is (n - m) x n
    STOP



References
==========

.. [Langhaar1951] Langhaar, H. L. (1951). Systematic Calculation of Dimensionless Products. In *Dimensional Analysis and Theory of Models* (pp. 29-46). John Wiley & Sons, Inc.
.. [Ngwua2021] Ngwua, L. (2021, July 4). *Theory of Dimensionless Products*. NeuralGraphs. `<https://www.neuralgraphs.com/lectures/diman/lectp8.html>`_.


.. |f| image:: ../resources/math/f.gif

.. |m| image:: ../resources/math/small_m.gif

.. |n| image:: ../resources/math/small_n.gif

.. |A_b| image:: ../resources/math/augmented_Ab.gif

.. |Ax_b| image:: ../resources/math/Ax_b.gif

.. |k_i| image:: ../resources/math/unknown_ks.gif
=========================
Example: Journal Bearing
=========================

.. image:: ../resources/images/journal_bearing.png
   :width: 350px
   :align: center

If one were interested in studying the frictional coefficient |f| of the bearing, then we must consider the variables/parameters that may influence it.

* bearing length, |L|
* bearing diameter, |D|
* bearing load, |P|

  - The load on the bearing is represented in terms of the average bearing pressure |P_eq_W_div_LD| where |W| is the actual load of bearing.

* rotating speed, |N|

  - Assume that the resulting rotating speed of the bearing is the constant average speed *N*.

* viscosity of lubricating oil, |mu|

  - This is the viscosity at equillibrium temperature &mdash; the bearing rotating at an average of |N| produces heat which is conducts and convects.

* clearance between bearing and journal, |C|
* bearing moment, |M|

  - Load applied to the shaft passing through the bearing results in bearing moment.

We can therefore start our study with the assumption that the frictional coefficient of the bearing is a function of the above seven variables resulting in the function value

.. image:: ../resources/math/tutorial3_func_value_of_LDPNmuCM.gif
   :align: center

But, |L|, |D| and |C| have the same dimensions. Then, |LbyD| and |CbyD| are dimensionless. Therefore, if we temporarily disregard the variables |L| and |C|, then we can reduce seven variables to five.

Hence, the tentative function for the derivation is such that its value is given by

.. image:: ../resources/math/tutorial3_func_value_of_DPNmuM.gif
   :align: center

Since,

+-----------------+---------------+----------------+------------+
| quantity symbol | quantity name | unit (say, SI) | dimensions |
+=================+===============+================+============+
| |P|	          | pressure      | |Pa|           | |MbyLT2|   |
+-----------------+---------------+----------------+------------+
| |M|             | moment        | |Nm|           | |ML2byT2|  |
+-----------------+---------------+----------------+------------+
| |D|	          | diameter      | |m|            | |L|        |
+-----------------+---------------+----------------+------------+
| |mu|            | viscosity     | |Pas|          | |MbyLT|    |
+-----------------+---------------+----------------+------------+
| |N|	          | speed         | |1_divby_s|    | |1byT|     |
+-----------------+---------------+----------------+------------+

the dimensional system for the problem is MLT-system.

We can now proceed with the steps (four) for deriving the dimensionless products.

1. Generate Dimensional Formula for All the Terms
=================================================

The derivation of the dimensionless products will be based on the reduced |f| where
the parent function depends on the independent five variables. The generation of
the dimensional matrix is follows some preceding setup steps.

1.1. Setup for Generation
-------------------------

1.1.1. Definitions setup
~~~~~~~~~~~~~~~~~~~~~~~~

Since our problem uses MLT dimensional system

::

    (def varpars [{:symbol "x", :quantity "mass"}
                  {:symbol "y", :quantity "length"}
                  {:symbol "t", :quantity "time"}])

1.1.2. Expressions and equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We express the variables for the unknown function |f| as

::

    (def manifold_eqn [{:name "term-P", :eqn {:term1 "x^(1)*y^(-1)*t^(-2)"}}
                       {:name "term-M", :eqn {:term1 "x^(1)*y^(2)*t^(-2)"}}
                       {:name "term-D", :eqn {:term1 "y^(1)"}}
                       {:name "term-mu", :eqn {:term1 "x^(1)*y^(-1)*t^(-1)"}}
                       {:name "term-N", :eqn {:term1 "t^(-1)"}}])

1.2. Getting the Dimensional Formula
------------------------------------

The dimensional formula all the terms are

::

    => (pprint (formula-eqn-side-manifold varpars manifold_eqn))
      [{:quantity "term-P", :dimension "[M^(1)*T^(-2)*L^(-1)]"}
       {:quantity "term-M", :dimension "[M^(1)*T^(-2)*L^(2)]"}
       {:quantity "term-D", :dimension "[L^(1)]"}
       {:quantity "term-mu", :dimension "[M^(1)*T^(-1)*L^(-1)]"}
       {:quantity "term-N", :dimension "[T^(-1)]"}]

1.3 Standardize All the Generated Dimensional Formula
-----------------------------------------------------

We add the above dimensional formulae into the `standard_formula` 

::

    => (update-sformula (formula-eqn-side-manifold varpars manifold_eqn))
      [{:quantity "volume", :dimension "[L^(3)]"}
       {:quantity "frequency", :dimension "[T^(-1)]"}
       {:quantity "velocity", :dimension "[L^(1)*T^(-1)]"}
       {:quantity "acceleration", :dimension "[L^(1)*T^(-2)]"}
       {:quantity "force", :dimension "[M^(1)*L^(1)*T^(-2)]"}
       ...
       {:quantity "term-N", :dimension "[T^(-1)]"}
       {:quantity "term-mu", :dimension "[M^(1)*T^(-1)*L^(-1)]"}
       {:quantity "term-D", :dimension "[L^(1)]"}
       {:quantity "term-M", :dimension "[M^(1)*T^(-2)*L^(2)]"}
       {:quantity "term-P", :dimension "[M^(1)*T^(-2)*L^(-1)]"}]

1.4. Definitions setup for dimensional matrix
---------------------------------------------

::

    (def varpars2 [{:symbol "P", :quantity "term-P"}
                   {:symbol "M", :quantity "term-M"}
                   {:symbol "D", :quantity "term-D"}
                   {:symbol "mu", :quantity "term-mu"}
                   {:symbol "N", :quantity "term-N"}])

2. Generate Dimensional Matrix
==============================

::

    => (view-matrix (generate-dimmat varpars2))
      [-1N 2N 1N -1N 0]
      [-2N -2N 0 -1N -1N]
      [1N 1N 0 1N 0]
      Size -> 3 x 5

3. Get the Homogeneous equation of the Dimensional Matrix
=========================================================

3.1. Get the augmented matrix of the dimensional matrix
-------------------------------------------------------

::

    => (view-matrix (get-augmented-matrix (generate-dimmat varpars2)))
      [1N -1N 0 1N -2N]
      [0 -1N -1N 2N 2N]
      [0 1N 0 -1N -1N]
      Size -> 3 x 5

3.2. Solve the augmented matrix
-------------------------------

::

    => (view-matrix (solve (get-augmented-matrix (generate-dimmat varpars2))))
      [1N 0N 0N 0N -3N]
      [0 1N 0N -1N -1N]
      [0 0N 1N -1N -1N]
      Size -> 3 x 5

3.3. Get the solution matrix
----------------------------

::

    => (view-matrix (get-solution-matrix (solve (get-augmented-matrix (generate-dimmat varpars2)))))
      [1 0 0N -1N -1N]
      [0 1 -3N -1N -1N]
      Size -> 2 x 5

This is a 2 &times; 5 matrix. Therefore, two dimensionless products will be derived.

We can put all these individual steps involving matrix into one coding step such that it returns the solution matrix.

::

    => (def solution_matrix (get-solution-matrix
                                (solve
                                    (get-augmented-matrix
                                        (generate-dimmat varpars2)))))
    => (view-matrix solution_matrix)
      [1 0 0N -1N -1N]
      [0 1 -3N -1N -1N]
      Size -> 2 x 5

4. Get Dimensionless Products
=============================

::

    => (def all_dimless (get-dimensionless-products solution_matrix varpars2))

    => (pprint all_dimless)
      [{:symbol "pi0", :expression "P^(1)*mu^(-1)*N^(-1)"}
       {:symbol "pi1", :expression "M^(1)*D^(-3)*mu^(-1)*N^(-1)"}]

    => (get-pi-expression all_dimless "pi0")
      "P^(1)*mu^(-1)*N^(-1)"

Note that these two dimensionless products are derived from the tentative function |f| where we temporarily disregarded |LbyD| and |CbyD|.

But, |LbyD| and |CbyD| are dimensionless. Therefore, the number of products in the complete set of dimensionless products is four.
Hence, the frictional coefficient has the value

.. image:: ../resources/math/tutorial3_func_value_of_pis.gif
   :align: center

   

.. |f| image:: ../resources/math/f.gif

.. |L| image:: ../resources/math/L.gif

.. |D| image:: ../resources/math/D.gif

.. |P| image:: ../resources/math/P.gif

.. |W| image:: ../resources/math/W.gif

.. |N| image:: ../resources/math/N.gif

.. |C| image:: ../resources/math/C.gif

.. |M| image:: ../resources/math/M.gif

.. |mu| image:: ../resources/math/mu.gif

.. |P_eq_W_div_LD| image:: ../resources/math/P_eq_W_div_LD.gif

.. |Pa| image:: ../resources/math/Pascal.gif

.. |Nm| image:: ../resources/math/NewtonMeter.gif

.. |m| image:: ../resources/math/Meter.gif

.. |Pas| image:: ../resources/math/PascalSecond.gif

.. |1_divby_s| image:: ../resources/math/1overSecond.gif

.. |MbyLT2| image:: ../resources/math/MbyLT2.gif

.. |ML2byT2| image:: ../resources/math/ML2byT2.gif

.. |MbyLT| image:: ../resources/math/MbyLT.gif

.. |1byT| image:: ../resources/math/1byT.gif

.. |LbyD| image:: ../resources/math/LbyD.gif

.. |CbyD| image:: ../resources/math/CbyD.gif
===============================
Deriving dimensionless products
===============================

For this introductory tutorial imagine the equation derived from experimental observations is

.. image:: ../resources/math/tutorial2_unknown_f.gif
   :align: center

**What are its dimensionless products?**

The four steps for deriving the dimensionless products are as follows.

1. Generate Dimensional Formula for All the Terms (usually right hand side of equation)
=======================================================================================

Let us define

|term_p|, |term_q|, |term_r|, |term_s|

Replacing |p|, |q|, |r| and |s| for the terms in the main equation we get

.. image:: ../resources/math/tutorial2_unknown_f_reduced.gif
   :align: center

Broadly, the task is to

* generate dimensional formula for each term
* insert all the generated dimensional formula into ``standard_formula``

1.1. Setup for Generation
-------------------------

1.1.1. Definitions setup
~~~~~~~~~~~~~~~~~~~~~~~~

Define all the symbols in the parent mathematical expression that is associated with a dimension.

::

    (def varpars [{:symbol "x", :quantity "mass"}
                  {:symbol "y", :quantity "length"}
                  {:symbol "t", :quantity "time"}])

1.1.2. Expressions and equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Considering each term of the parent equation as some individual equation we define each of them as follows

::

    (def p_equation {:lhs "p^(1)", :rhs {:term1 "x^(2)*y^(-1)*t^(1)"}})
    (def q_equation {:lhs "q^(1)", :rhs {:term1 "x^(1)*y^(6)*t^(20)"}})
    (def r_equation {:lhs "r^(1)", :rhs {:term1 "x^(3)*y^(-3)*t^(-3)"}})
    (def s_equation {:lhs "s^(1)", :rhs {:term1 "x^(4)*t^(8)"}})


and then stacking the equations into a vector

::

    (def manifold_eqn [{:name "term-p", :eqn (:rhs p_equation)}
                       {:name "term-q", :eqn (:rhs q_equation)}
                       {:name "term-r", :eqn (:rhs r_equation)}
                       {:name "term-s", :eqn (:rhs s_equation)}])

The two steps are equivalent to

::

    (def manifold_eqn [{:name "term-p", :eqn {:term1 "x^(2)*y^(-1)*t^(1)"}}
                       {:name "term-q", :eqn {:term1 "x^(1)*y^(6)*t^(20)"}}
                       {:name "term-r", :eqn {:term1 "x^(3)*y^(-3)*t^(-3)"}}
                       {:name "term-s", :eqn {:term1 "x^(4)*t^(8)"}}])

However, the two step approach is recommended because it affords the user the flexibility to actually see the generation of individual dimensional formula and hence introspecting them.

1.2. Getting the Dimensional Formula
------------------------------------

The dimensional formula for one side of the expression (often right hand side) for every equation in the vector of all the equations defined earlier can be generated using the ``formula-eqn-side-manifold`` function.

::

    => (pprint (formula-eqn-side-manifold varpars manifold_eqn))
    [{:quantity "term-p", :dimension "[L^(-1)*M^(2)*T^(1)]"}
    {:quantity "term-q", :dimension "[L^(6)*M^(1)*T^(20)]"}
    {:quantity "term-r", :dimension "[T^(-3)*L^(-3)*M^(3)]"}
    {:quantity "term-s", :dimension "[M^(4)*T^(8)]"}]

1.3 Standardize All the Generated Dimensional Formula
-----------------------------------------------------

All the dimensional formula generated from each equation in the vector of equations is added to the ``standard_formula`` with 

::

    => (update-sformula (formula-eqn-side-manifold varpars manifold_eqn))
    [{:quantity "volume", :dimension "[L^(3)]"}
    {:quantity "frequency", :dimension "[T^(-1)]"}
    {:quantity "velocity", :dimension "[L^(1)*T^(-1)]"}
    {:quantity "acceleration", :dimension "[L^(1)*T^(-2)]"}
    {:quantity "force", :dimension "[M^(1)*L^(1)*T^(-2)]"}
    ...
    {:quantity "term-s", :dimension "[M^(4)*T^(8)]"}
    {:quantity "term-r", :dimension "[T^(-3)*L^(-3)*M^(3)]"}
    {:quantity "term-q", :dimension "[L^(6)*M^(1)*T^(20)]"}
    {:quantity "term-p", :dimension "[L^(-1)*M^(2)*T^(1)]"}]

1.4. Definitions setup for the reduced form of the parent expression
--------------------------------------------------------------------

Since all the dimensional formula of |p|, |q|, |r| and |s|, representing all the terms in the main equation are now part of the ``standard_formula``, we can now define all the symbols in the reduced form of the parent mathematical expression

.. image:: ../resources/math/tutorial2_unknown_f_reduced.gif
   :align: center

The definition will be such that each term symbol has the dimension name as defined in the preceeding step (and hence incorporated into the ``standard_formula``. For instance, since the term |p| (|term_p|) was named ``"term-p"`` in

::

    => (pprint manifold_eqn)
    [{:name "term-p", :eqn {:term1 "x^(2)*y^(-1)*t^(1)"}}
    {:name "term-q", :eqn {:term1 "x^(1)*y^(6)*t^(20)"}}
    {:name "term-r", :eqn {:term1 "x^(3)*y^(-3)*t^(-3)"}}
    {:name "term-s", :eqn {:term1 "x^(4)*t^(8)"}}]

we will have ``{:symbol "p", :dimension "term-p"}``. Therefore, we define

::

    (def varpars2 [{:symbol "p", :quantity "term-p"}
                   {:symbol "q", :quantity "term-q"}
                   {:symbol "r", :quantity "term-r"}
                   {:symbol "s", :quantity "term-s"}])

2. Generate Dimensional Matrix
==============================

The dimensional matrix of the parent equation is generated with the help of the ``generate-dimmat`` function.

::

    => (view-matrix (generate-dimmat varpars2))
    [1N 20N -3N 8N]
    [2N 1N 3N 4N]
    [-1N 6N -3N 0]
    Size -> 3 x 4

This is a 3 times 4 dimensional matrix.

3. Get the Homogeneous equation of the Dimensional Matrix
=========================================================

3.1. Get the augmented matrix of the dimensional matrix
-------------------------------------------------------


::

    => (view-matrix (get-augmented-matrix (generate-dimmat varpars2)))
    [-3N 8N -1N -20N]
    [3N 4N -2N -1N]
    [-3N 0 1N -6N]
    Size -> 3 x 4

3.2. Solve the augmented matrix
-------------------------------

::

    => (view-matrix (solve (get-augmented-matrix (generate-dimmat varpars2))))
    [1N 0N -1/3 2N]
    [0N 1N -1/4 -7/4]
    [0N 0N 0N 0N]
    Size -> 3 x 4

3.3. Get the solution matrix
----------------------------

::

    => (view-matrix (get-solution-matrix (solve (get-augmented-matrix (generate-dimmat varpars2)))))
    [1 0 -1/3 -1/4]
    [0 1 2N -7/4]
    Size -> 2 x 4

This is a 2 times 4 matrix. Therefore, there will be two dimensionless products.

We can put all these individual steps involving matrix into one coding step such that it returns the solution matrix.

::

    => (def solution_matrix (get-solution-matrix
                                (solve
                                    (get-augmented-matrix
                                        (generate-dimmat varpars2)))))
    => (view-matrix solution_matrix)
    [1 0 -1/3 -1/4]
    [0 1 2N -7/4]
    Size -> 2 x 4

4. Get Dimensionless Products
=============================

The dimensionless products are generated with the help of the ``get-dimensionless-products`` function.

::

    => (pprint (get-dimensionless-products solution_matrix varpars2))
    [{:symbol "pi0", :expression "p^(1)*r^(-1/3)*s^(-1/4)"}
    {:symbol "pi1", :expression "q^(1)*r^(2)*s^(-7/4)"}]

Since, &pi; is the conventional symbol for dimensionless products to get the |pi_i| th one use the ``get-pi-expression`` function. For example, for |pi0|

::

    => (def all-dimless (get-dimensionless-products solution_matrix varpars2))
    => (get-pi-expression all-dimless "pi0")
    "p^(1)*r^(-1/3)*s^(-1/4)"

    
.. |p| image:: ../resources/math/small_p.gif

.. |q| image:: ../resources/math/small_q.gif

.. |r| image:: ../resources/math/small_r.gif

.. |s| image:: ../resources/math/small_s.gif

.. |term_p| image:: ../resources/math/tutorial2_term_p.gif

.. |term_q| image:: ../resources/math/tutorial2_term_q.gif

.. |term_r| image:: ../resources/math/tutorial2_term_r.gif

.. |term_s| image:: ../resources/math/tutorial2_term_s.gif

.. |pi_i| image:: ../resources/math/pi_i.gif

.. |pi0| image:: ../resources/math/pi0.gif
==============================================================
Generate dimensional formulae and perform consistency checking
==============================================================

This introduction will teach how to

* generate dimensional formulae
* generate consistency checks
* standardize correct equation

How to Generate Dimensional Formulae
====================================

1. The Problem
--------------

Imagine you derived the equation below based on the experimental observations

.. image:: ../resources/math/tutorial1_unknown_f_generate_dimmat.gif
   :align: center

You want to know if this derived equation is correct. Using **diman** you can perform a *preliminary check with consistency analysis*. But before you can check for dimensional consistency you need to set it up for the analysis.

2. Setting Up for Generation
----------------------------

2.1. Definitions setup
~~~~~~~~~~~~~~~~~~~~~~

Define all the symbols in the mathematical expression that is associated with a dimension.

::

    (def varpars [{:symbol "x", :quantity "length"}
                  {:symbol "v", :quantity "velocity"}
                  {:symbol "t", :quantity "time"}
                  {:symbol "a", :quantity "acceleration"}])

2.2. Expressions and equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, define the equation in terms of its left and right hand side of the expression. If a side has more than one term they are expressed as a map with appropriate key ``:termi`` for respective term.

::

    (def lhs "x^(1)")
    (def rhs {:term1 "x^(1)",
              :term2 "v^(2)",
              :term3 "t^(1)",
              :term4 "0.5*a^(1)*t^(2)"})
    (def eqn {:lhs lhs, :rhs rhs})

3. Getting the Dimensional Formula
----------------------------------

The equation defined above is used for deriving the dimensional formula.

3.1. Sub-formula of the dimensional formula for one side of the equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the side of the equation that contains multiple terms the *sub-formula* is the dimensional formula for one of the terms. Notice that the sub-formula **IS** the dimensional formula for the expression if there is just one term.

Thus for the right hand side of the given equation (which was defined in the previous section)

::

    => rhs
    {:term1 "x^(1)", :term2 "v^(2)", :term3 "t^(1)", :term4 "0.5*a^(1)*t^(2)"}

the dimensional formula for ``:term4`` is given by

::

    => (formula-term varpars (:term4 rhs))
    "[T^(0)*L^(1)]"


3.2. Dimensional formula for one side of the equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dimensional formula for one side of the expression regardless of the number of terms in it can be generated using the ``formula-eqn-side`` function.

Dimensional formula for the ``rhs`` expression is given by

::

    => (formula-eqn-side varpars rhs)
    "[L^(1)] + [T^(-2)*L^(2)] + [T^(1)] + [T^(0)*L^(1)]"

3.3. Introspecting the dimensional formula
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3.3.1. Represent sub-formula of an expression term as dimension names
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

::

    => (dimnames (formula-term varpars (:term4 rhs)))
    "length^(1)"

3.3.2. Represent dimensional formula of an equation side as dimension names
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

::

    => (dimnames (formula-eqn-side varpars rhs))
    "length^(1) + time^(-2)*length^(2) + time^(1) + length^(1)"

How to do Consistency Checks
============================

1. The Problem
--------------

Consider the equation

.. image:: ../resources/math/tutorial1_unknown_f_consistency_check.gif
   :align: center

2. Setting Up for Checking
--------------------------

2.1. Definitions setup
~~~~~~~~~~~~~~~~~~~~~~

Define all the symbols in the mathematical expression that is associated with a dimension.

::

    (def varpars [{:symbol "x", :quantity "length"}
                  {:symbol "v", :quantity "velocity"}
                  {:symbol "t", :quantity "time"}
                  {:symbol "a", :quantity "acceleration"}])

2.2. Expressions and equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, define the equation in terms of its left and right hand side of the expression. If a side has more than one term they are expressed as a map with appropriate key ``:termi`` for respective term.

::

    (def lhs "x^(1)")
    (def rhs {:term1 "x^(1)",
              :term2 "v^(1)*t^(1)",
              :term3 "0.5*a^(1)*t^(2)"})
    (def eqn {:lhs lhs, :rhs rhs})

3. Perform Consistency Check
----------------------------

If the correctness of an equation is in doubt checking for dimensional consistency is a useful preliminary step.

To perform consistency check based on dimensional analysis in diman <sup>(c)</sup> you use the predicate function ``consistent?``. Thus, for the given problem

::

    => (consistent? varpars eqn)
    true

However, dimensionally consistent equation **does not guarantee** correct equation.

4. Consistency of multiple equations
------------------------------------

Let us consider the case of a problem where one derives multiple expressions thought to be potential candidates for representing the problem.

.. image:: ../resources/math/tutorial1_e_m2v2.gif
   :align: center
.. image:: ../resources/math/tutorial1_e_half_mv2.gif
   :align: center
.. image:: ../resources/math/tutorial1_e_ma.gif
   :align: center
.. image:: ../resources/math/tutorial1_e_3by16_mv2.gif
   :align: center
.. image:: ../resources/math/tutorial1_e_half_mv2_plus_ma.gif
   :align: center

the question is, **which of these equations are correct?** To tackle this question let us first look at the answer for *which of these equations are dimensionally correct?* In other words, let us perform dimensional consistency checks on each expression.

Thus

+---------------------------------------------------------------+--------------------------------------------------------------------------------------+
| Equation                                                      | Setup                                                                                |
+===============================================================+======================================================================================+
| .. image:: ../resources/math/tutorial1_e_m2v2.gif             | ``(def eqn1 {:lhs "e^(1)", :rhs "m^(2)*v^(2)"})``                                    |
+---------------------------------------------------------------+--------------------------------------------------------------------------------------+
| .. image:: ../resources/math/tutorial1_e_half_mv2.gif         | ``(def eqn2 {:lhs "e^(1)", :rhs "0.5*m^(1)*v^(2)"})``                                |
+---------------------------------------------------------------+--------------------------------------------------------------------------------------+
| .. image:: ../resources/math/tutorial1_e_ma.gif               | ``(def eqn3 {:lhs "e^(1)", :rhs "m^(1)*a^(1)"})``                                    |
+---------------------------------------------------------------+--------------------------------------------------------------------------------------+
| .. image:: ../resources/math/tutorial1_e_3by16_mv2.gif        | ``(def eqn4 {:lhs "e^(1)", :rhs "0.1875*m^(1)*v^(2)"})``                             |
+---------------------------------------------------------------+--------------------------------------------------------------------------------------+
| .. image:: ../resources/math/tutorial1_e_half_mv2_plus_ma.gif | ``(def eqn5 {:lhs "e^(1)", :rhs {:term1 "0.5*m^(1)*v^(2)", :term2 "m^(1)*a^(1)"}})`` |
+---------------------------------------------------------------+--------------------------------------------------------------------------------------+

and define the variables/parameters as

::

    (def varpars [{:symbol "e", :quantity "energy"}
                  {:symbol "m", :quantity "mass"}
                  {:symbol "v", :quantity "velocity"}
                  {:symbol "a", :quantity "acceleration"}])

Then

::

    => (consistent? varpars eqn1)
    false
    => (consistent? varpars eqn2)
    true
    => (consistent? varpars eqn3)
    false
    => (consistent? varpars eqn4)
    true
    => (consistent? varpars eqn5)
    false

which suggests |e_half_mv2| and |e_3by16_mv2| to be dimensionally consistent.

But both equations can't be correct, illustrating the point that

    a dimensionally consistent equation does not guarantee correct equation

How to Standardize the Correct Equation
=======================================

From the previous example of notice that kinetic ``"e"`` is not defined in the ``standard_formula``

::

    => (pprint standard_formula)
    [{:quantity "volume", :dimension "[L^(3)]"}
    {:quantity "velocity", :dimension "[L^(1)*T^(-1)]"}
    {:quantity "acceleration", :dimension "[L^(1)*T^(-2)]"}
    {:quantity "force", :dimension "[M^(1)*L^(1)*T^(-2)]"}
    {:quantity "mass density", :dimension "[M^(1)*L^(-3)]"}]

Since we already know that the kinetic energy is in Joules and |1Joule| whose dimensional formula is ``"[M^(1)*L^(2)*T(-2)]"`` this can be added to the ``standard_formula`` as

::

    => (update-sformula [{:quantity "energy", :dimension "[M^(1)*L^(2)*T(-2)]"}])
    [{:quantity "volume", :dimension "[M^(0)*L^(3)*T^(0)]"}
    {:quantity "velocity", :dimension "[M^(0)*L^(1)*T^(-1)]"}
    {:quantity "acceleration", :dimension "[M^(0)*L^(1)*T^(-2)]"}
    {:quantity "force", :dimension "[M^(1)*L^(1)*T^(-2)]"}
    {:quantity "mass density", :dimension "[M^(1)*L^(-3)*T^(0)]"}
    {:quantity "energy", :dimension "[M^(1)*L^(2)*T(-2)]"}]

Now since ``"energy"`` is one of the ``:quantity`` in the ``standard_formula``, we can now add the symbol ``"e"`` in our definition as follows

::

    => (def varpars (conj varpars {:symbol "e", :quantity "energy"}))
    => (pprint varpars)
    [{:symbol "m", :quantity "mass"}
    {:symbol "v", :quantity "velocity"}
    {:symbol "a", :quantity "acceleration"}
    {:symbol "e", :quantity "energy"}]


.. |e_half_mv2| image:: ../resources/math/tutorial1_e_half_mv2.gif

.. |e_3by16_mv2| image:: ../resources/math/tutorial1_e_3by16_mv2.gif

.. |1Joule| image:: ../resources/math/Joule.gif
