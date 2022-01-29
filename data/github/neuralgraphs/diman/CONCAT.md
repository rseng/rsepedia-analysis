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
