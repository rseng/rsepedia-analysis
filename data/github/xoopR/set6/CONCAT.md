set6
================

<img src="man/figures/logo.png" align="right" alt="" width="120" />

[![R CMD Check via
{tic}](https://github.com/xoopR/set6/workflows/R%20CMD%20Check%20via%20%7Btic%7D/badge.svg?branch=master)](https://github.com/xoopR/set6/actions)
[![codecov](https://app.codecov.io/gh/xoopR/set6/branch/master/graph/badge.svg)](https://app.codecov.io/gh/xoopR/set6)
[![CodeFactor](https://www.codefactor.io/repository/github/xoopr/set6/badge)](https://www.codefactor.io/repository/github/xoopr/set6)
[![dependencies](https://tinyverse.netlify.com/badge/set6)](https://CRAN.R-project.org/package=set6)

[![Repo
Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle
Badge](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://img.shields.io/badge/lifecycle-stable-brightgreen)

[![CRAN Status
Badge](https://www.r-pkg.org/badges/version-ago/set6)](https://cran.r-project.org/package=set6)
[![CRAN
Checks](https://cranchecks.info/badges/worst/set6)](https://cran.r-project.org/web/checks/check_results_set6.html)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/set6)](https://cran.r-project.org/package=set6)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02598/status.svg)](https://doi.org/10.21105/joss.02598)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Gitter
chat](https://badges.gitter.im/xoopR/set6.png)](https://gitter.im/xoopR/set6)

## What is set6?

`set6` is an R6 upgrade to the `sets` package in R that includes:

-   Multi-dimensional sets
-   Tuples
-   Finite and infinite intervals
-   Fuzzy sets and tuples
-   Set operations including union, intersect, (asymmetric and
    symmetric) difference, and product
-   Symbolic representation of infinite sets including common special
    sets such as the Reals and Integers
-   ConditionalSets for defining sets according to logical conditions

## Installation

The current CRAN release can be installed with

``` r
install.packages("set6")
```

Or for the latest stable build

``` r
remotes::install_github("xoopR/set6")
```

## Main Features

### A Clear Inheritance Structure

``` r
# Sets require elements to be unique and order doesn't matter
Set$new(1, 2, 1) == Set$new(1, 2)
#> [1] TRUE
Set$new(1, 2) == Set$new(2, 1)
#> [1] TRUE

# But tuples can enforce these restrictions
Tuple$new(1, 2, 1) != Tuple$new(1, 2)
#> [1] TRUE
Tuple$new(1, 2) != Tuple$new(2, 1)
#> [1] TRUE

# Fuzzy sets and tuples extend sets further
f = FuzzySet$new(1, 0, 2, 0.6, 3, 1)
f$inclusion(1)
#> [1] "Not Included"
f$inclusion(2)
#> [1] "Partially Included"
f$inclusion(3)
#> [1] "Fully Included"

# Symbolic intervals provide a clean way to represent infinite sets
Interval$new()
#> [-∞,+∞]
# Different closure types and classes are possible
Interval$new(1, 7, type = "(]") # half-open
#> (1,7]
Interval$new(1, 7, class = "integer") == Set$new(1:7)
#> [1] TRUE

# And SpecialSets inheriting from these intervals
Reals$new()
#> ℝ
PosRationals$new()
#> ℚ+
```

### Set operations

``` r
# Union
Set$new(1, 4, "a", "b") + Set$new(5)
#> {1, 4,...,a, b}
Interval$new(1, 5) + FuzzyTuple$new(1, 0.6)
#> [1,5]

# Power
Set$new(1:5)^2
#> {1, 2,...,4, 5}^2
# A symbolic representation is also possible
setpower(Set$new(1:5), power = 2, simplify = FALSE)
#> {1, 2,...,4, 5}^2
Reals$new()^5
#> ℝ^5

# Product
Set$new(1,2) * Set$new(5, 6)
#> {1, 2} × {5, 6}
Interval$new(1,5) * Tuple$new(3)
#> [1,5] × (3)

# Intersection
Set$new(1:5) & Set$new(4:10)
#> {4, 5}
ConditionalSet$new(function(x) x == 0) & Set$new(-2:2)
#> {0}
Interval$new(1, 10) & Set$new(5:6)
#> {5, 6}

# Difference
Interval$new(1, 10) - Set$new(5)
#> [1,5) ∪ (5,10]
Set$new(1:5) - Set$new(2:3)
#> {1, 4, 5}
```

### Containedness and Comparators

``` r
Interval$new(1, 10)$contains(5)
#> [1] TRUE
# check multiple elements
Interval$new(1, 10)$contains(8:12)
#> [1]  TRUE  TRUE  TRUE FALSE FALSE
# only return TRUE if all are TRUE
Interval$new(1, 10)$contains(8:12, all = TRUE)
#> [1] FALSE
# decide whether open bounds should be included
Interval$new(1, 10, type = "()")$contains(10, bound = TRUE)
#> [1] TRUE
Interval$new(1, 10, type = "()")$contains(10, bound = TRUE)
#> [1] TRUE

Interval$new(1, 5, class = "numeric")$equals(Set$new(1:5))
#> [1] FALSE
Interval$new(1, 5, class = "integer")$equals(Set$new(1:5))
#> [1] TRUE

Set$new(1) == FuzzySet$new(1, 1)
#> [1] TRUE

# proper subsets
Set$new(1:3)$isSubset(Set$new(1), proper = TRUE)
#> [1] TRUE
Set$new(1) < Set$new(1:3)
#> [1] TRUE

# (non-proper) subsets
Set$new(1:3)$isSubset(Set$new(1:3), proper = FALSE)
#> [1] TRUE
Set$new(1:3) <= Set$new(1:3)
#> [1] TRUE

# multi-dimensional checks
x = PosReals$new()^2
x$contains(list(Tuple$new(1, 1), Tuple$new(-2, 3)))
#> [1]  TRUE FALSE
```

## Usage

The primary use-cases of `set6` are:

1.  **Upgrading sets** Extend the R `sets` package to R6, which allows
    for generalised `Set` objects with a clear inheritance structure. As
    well as adding features including symbolic representation of
    infinite sets, and cartesian products.
2.  **Defining parameter interfaces** All objects inheriting from the
    `Set` parent class include methods `equals` and `contains`, which
    are used to check if two sets are equal or if elements lie in the
    given set. This makes `set6` ideal for parameter interfaces in which
    a range of values (possibly multi-dimensional or of mixed types)
    need to be defined.

## Short-term development plans

Whilst the `set6` API is stable, it is considered ‘maturing’, and
therefore whilst there are no plans for major updates, these may still
occur. There are a few features and refactoring we plan on implementing
before we consider the package to be in its first complete version.
These mainly include

-   Finalising all methods and fields - some are missing or possibly
    inaccurate for some wrappers. For example the cardinality of
    `ComplementSet`s is imprecise at the moment.
-   We are considering adding a `simplify` method to wrappers to reduce
    classes inheriting from `SetWrapper` to simpler sets. This allows
    users to perform operations with `simplify = FALSE` and then to
    change their mind.
-   There are known bottlenecks that need to be fixed to massively
    improve speed and efficiency.
-   Adding more tutorials to make the interface easier for beginners,
    especially people new to R6

At a later stage we may consider adding Venn diagrams for visualisation
of sets and intervals, but this is very low priority.

## Similar Packages

-   [sets](https://CRAN.R-project.org/package=sets) - The **sets**
    package uses S3 to define some symbolic representaton of
    mathematical sets, tuple, intervals, and fuzzy variants. However the
    symbolic representation is not consistent throughout the package,
    does not allow for clear inspection of set/interval elements, and
    there is no support for multi-dimensional sets.

-   [BaseSet](https://github.com/ropensci/BaseSet) - The **BaseSet**
    package focuses on storing and analysing sets in a ‘tidy’ way, with
    more options for data storage in long and wide formats. The primary
    usage is neat and efficient inspection of finite sets, there is no
    support for infinite sets or symbolic representation.

## Contributing

As `set6` is in its early stages, contributions are very welcome. If you
have any ideas for good features please open an issue or create a pull
request. Otherwise bug reports are very appreciated if you stumble
across any broken code, these can be posted to the [issue
tracker](https://github.com/xoopR/set6/issues). For updates on `set6`
follow/star this repo.

## Citing set6

If you use set6, please cite our [JOSS
article](https://doi.org/10.21105/joss.02598):

@Article{set6, title = {set6: R6 Mathematical Sets Interface}, author =
{Raphael Sonabend and Franz J. Kiraly}, journal = {Journal of Open
Source Software}, year = {2020}, month = {nov}, doi =
{10.21105/joss.02598}, url =
{<https://joss.theoj.org/papers/10.21105/joss.02598>}, }
# set6 0.2.4.9000

* Fixes bug that set all comparisons to empty sets as `TRUE` when using `$equals`

# set6 0.2.4

* Patch for R devel
* Imported ooplah

# set6 0.2.3

* Containedness checks for 'n' dimensional sets no longer require same length vectors if power is
"n"

# set6 0.2.2
* Fixed bug preventing `Logicals` from being deep cloned

# set6 0.2.1

* Bugfix in setcomplement (#65)
* Impossible intervals containing only one elements with type not equal to `[]` are now equal to the empty set
* Default `ConditionalSet` `condition` argument now `function(x) TRUE`
* Print method for `ConditionalSet` now omits RHS if only `"TRUE"`

# set6 0.2.0

* UniversalSet renamed Universal, old class will be removed in v0.4.0.
* LogicalSet renamed Logicals, old class will be removed in v0.4.0.
* `Complex` now inherits from `Set`, incorrect methods for `isSubset, equals` have been removed.
* Add `Multiset` for sets with non-unique elements but no ordering
* Small speed improvements in `Tuple` and `FuzzyTuple`
* For consistency most methods now return a `list` unless single elements requested
* Printing of `ConditionalSet` is fixed
* `Rationals` and child-classes now error on calls to `contains, isSubset, equals` as any prior results were likely wrong/misleading
* Removed erroneous complex boundaries in `Interval` class

# set6 0.1.8

* Patch for R-devel

# set6 0.1.7

* Critical patch

# set6 0.1.6

* Bugfix in set operation cleaner
* Bugfix causing `Interval$contains` to return `TRUE` for tuples
* Bugfix in union sets incorrectly unwrapping products
* Added variable length `ExponentSet`s

# set6 0.1.5

* Added `LogicalSet`, the set of $\{TRUE, FALSE\}$
* Added `as.Set.numeric` and `as.Tuple.numeric`

# set6 0.1.4

* Speed performance improvements for `$contains` method for `Interval` and `Set`. `Rcpp` now used for `Interval`.
* Now for any `Interval` if not bounded above and `upper` is $Inf$ then `max = .Machine$double.xmax`, analogously for `lower`.
* Default universe of `Interval` is now `ExtendedReals`
* Added default `as.Set` and `as.Interval` S3 methods

# set6 0.1.3

* Added assertion for testing if a set is countably finite
* Slight speed improvements to operations - still require a lot of work
* Fixed bug in `UnionSet` cardinality calculation
* Fixed bug in `UniversalSet` countability

# set6 0.1.2

### Patches
- Updated documentation to be compatible with roxygen2
- Fixed bug in typed Complex sets
- Added universe assertion check to `Set` constructor
- Bug fix in `setunion` causing some intervals not to be combined correctly
- `Interval$isSubset` now compares sets using `max` and `min` instead of `upper` and `lower`
- Calculation of `min` and `max` in `Interval` now uses `1e-15` instead of `.Machine$double.xmin`
- `$elements` now always returns a `list`

### Added classes, methods, and functions
- Add `$add` public method to sets, which mutates sets by adding given elements, and coercing to the typed-set class if appropriate
- Add `$remove` public method to sets, which mutates sets by removing given elements.
- Add assertion for checking if elements contained in a set, `test/check/assertContains`.
- Add assertion for checking if sets are subsets of another, `test/check/assertSubset`.

# set6 0.1.1

### Patches
- `absComplement` method is now deprecated, instead use `setcomplement` and omit the `y` argument
- Fixed error in `contains` default caused by `%inset%`
- Improved printing of `SpecialSet`s when `zero == TRUE`
- Added `UniversalSet` for the set containing all elements
- Changed default `universe` of sets to `UniversalSet`
- Coercions now error instead of producing a message when they fail
- On construction, `Set`s no longer guess the set class, instead an extra `class` argument is added to give a set the `typed` property
- The internal `Set` structure is slightly changed so that set elements are now stored in lists by default, which is only changed if the set is `typed`
- Added `element` argument to `Set` constructor, which takes a `list`. This is more efficient if passing lists of lists or lists of multiple types, and in line with the `FuzzySet` constructor
- Improved printing of `ConditionalSet`s
- Updated `powerset` to always return a `Set` of `Set`s (even if input is `Tuple`)
- Fixed bug in `Properties` causing an error if cardinality was too large
- Updated documentation
- Reduced `Set` constructor bottleneck by adding 'typed' sets
- Changed `use_unicode` default to `l10n_info()$UTF-8`

# set6 0.1.0

- `set6` upgrades the `sets` package to R6. Many forms of mathematical sets are implemented, including (countably finite) sets, tuples, intervals (countably infinite or uncountable), and fuzzy variants. Wrappers extend functionality by allowing symbolic representations of complex operations on sets, including unions, (cartesian) products, exponentiation, and differences (asymmetric and symmetric).
- See [the website](https://xoopR.github.io/set6/) for more details and the project readme
- See [getting started vignette](https://xoopR.github.io/set6/articles/set6.html) for a short tutorial and introduction
- `set6` is currently 'maturing', so whilst no major updates are planned they may happen. Constant minor updates should be expected.
## Comments

Early submission due to bug in downstream dependency.

## Test Environments

On GitHub actions:

* windows-latest (R release)
* macOS-latest (R release)
* ubuntu-latest (R release)
* ubuntu-latest (R oldrel)
* ubuntu-latest (R devel)

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
All downstream dependencies OK.
---
title: 'set6: R6 Mathematical Sets Interface'
tags:
  - R
  - statistics
  - sets
  - object-oriented
authors:
  - name: Raphael Sonabend
    orcid: 0000-0001-9225-4654
    affiliation: 1
  - name: Franz J. Kiraly
    affiliation: 1
affiliations:
 - name: Department of Statistical Science, University College London
   index: 1
date: 23 March 2020
bibliography: paper.bib
---

# Summary

`set6` makes use of the R6 object-oriented paradigm to introduce classes for important mathematical objects, including sets, tuples, and intervals (finite and infinite). Until now, the `R` [@packageR] programming language has traditionally supported mathematical sets in one of two ways: 1. via the five set operation functions: `union`, `intersect`, `setdiff`, `setequal`, `is.element`; and 2. via the `sets` [@pkgsets] package. `set6` uses `R6` [@packageR6] and has a clear class interface with minimal dependencies, which makes it the perfect dependency for any package that requires mathematical data types as R6 objects. Making use of design patterns [@Gamma1996], such as wrappers and compositors, `set6` allows for symbolic representation of sets to ensure maximum efficiency, and to provide neat and clear print methods.

`set6` is currently being used in `distr6` [@packagedistr6], which is an object-oriented probability distributions interface, that makes use of `set6` for distribution and parameter support. Additional uses of `set6` include representing infinite sets, and constructing assertions.

The speed and efficiency of ``R6`` and `Rcpp` [@pkgrcpp] allows `set6` to be a scalable and efficient interface. A focus on symbolic representation and neat printing methods means that `set6` can accurately and clearly represent complicated sets. `set6` has the ambitious long-term goal of being the only dependency package needed for object-oriented interfaces in `R` that require clear symbolic representations of mathematical sets.

Related software in `R` includes `sets` [@pkgsets], which uses the S3 and S4 object-oriented paradigms to implement mathematical sets. Whilst `sets` and `set6` have similar implemented features, as both offer multiple forms of sets (intervals, fuzzy, etc.) and with symbolic representation, the `sets` package does not offer lazy evaluation in set operations, which leads to significant overheads for large or possibly-infinite sets. Similar packages in other computing languages include: i) `portion` [@pkgportion] in Python, which only supports intervals without generalisation to other mathematical sets; ii) `IntervalSets.jl` [@pkgintervalsets] in Julia, which is again limited to intervals though with good features for infix operators and inspection; iii) a variety of packages in the `JuliaIntervals` suite (https://juliapackages.com/u/juliaintervals), which primarily focus on rigorous arithmetic implementation; and iv) `LazySets.jl` [@pkglazysets] and `DomainSets.jl` [@pkgdomainsets] which both faciliate symbolic representation of sets and infinite intervals in Julia.

The example below demonstrates construction of a set and interval, comparisons of these, the set complement operator, and printing of the final result. Note the example does not use unicode printing but that this is possible within the package.

```R
> a <- Set$new(1, 2, 3)
> a$print()
{1, 2, 3}
> a$contains(c(1, "a"))
[1] TRUE FALSE
> b <- Interval$new(1, Inf, class = "numeric")
> b$isSubset(a)
TRUE
> c <- b - a
> c$print()
(1,2) U (2,3) U (3,+Inf] 
```



# Key Design Principles

1. **Maximum user-control over set operations** - Users can select how operations act on sets, including a choice of  associativity, lazy evaluation, and unicode printing.
2. **Minimal dependencies** - `set6` has the goal of being a key dependency to any object-oriented package requiring representation of mathematical sets, for example for representing function inputs and supports. Therefore `set6` is itself dependent on only three packages.
3. **Inspectability and reactive user interface** - `set6` prioritises symbolic representation and lazy evaluation to allow for the package to be scalable and to fit into any other package. However it is ensured that this does not detract from a clear user interface, i.e. at all times it should be clear what an object contains both externally (how it prints) and internally (inspection methods). `set6` allows sets to be queried in many different ways, including calling the elements in the set (if finite), finding the bounds of the set (if numeric), and listing properties and traits.

```R
# Symbolic representation allows neat representation of infinite sets
# and intervals whilst making the elements clear by inspection

> I <- Interval$new(10, Inf, type = "[)")
> I$print()
[10,Inf)
> I$contains(9:11)
[1] FALSE TRUE TRUE
 
> n <- Naturals$new()
> n$print()
N0
# binary operators also available
> c(pi, 2) %inset% n
[1] FALSE TRUE
 
```

4. **Combination of lazy and greedy evaluation** - By default, 'multiplying' operations such as products and powersets are evaluated lazily, whereas 'adding' operations such as unions and differences are evaluated greedily. These prevent system crashes from trying to evaluate sets of very large cardinality. In all cases, the user can override defaults. Symbolic representation is used with lazy evaluation so that large sets can be printed neatly without the individual elements ever being evaluated.

```R
# Lazy evaluation allows very large sets to be queried without
# explicit evaluation

> s <- Set$new(1, 2, 3)
> p <- powerset(powerset(powerset(s)))
> p$print()
2^2^2^{1, 2, 3}
> p$properties$cardinality
[1] 1.157921e+77
> p$contains(Set$new(Set$new(Set$new(1), Set$new(1, 2, 3))))
[1] TRUE
```



# Key Use-Cases

1. **Constructing and querying mathematical sets** - Many mathematical Set-like objects can be constructed, including sets, tuples, intervals, and fuzzy variants. Sets can contain objects of any `R` type that have a valid `as.character` coercion method that creates a unique symbolic representation of the class; this is required as internally checks are performed symbolically on character representations.
2. **Containedness checks** - Public methods allow all objects inheriting from `Set` to check if elements are contained within them. This provides a powerful mechanism for use with parameter or distribution supports for other packages as it can be viewed as a 'type check', i.e. checks if a value fits within a specified mathematical type. A C++ implementation of these checks in Rcpp [@pkgrcpp] means that the computations are incredibly quick for sets and intervals of any size.
3. **Representation of infinite sets** - Symbolic representation and lazy evaluation allows infinite (or very large) sets and intervals to be constructed. This also allows operations such as powerset to be used without crashing the system.
4. **Comparison of, possibly infinite, sets** - Two `Set` objects can be compared to check if they are equal or (proper) sub/supersets. Infix operators allow quick and neat comparison.
5. **Creation of composite sets from simpler classes** - Common set operations, such as unions and complements are implemented, as well as products and exponents. These make use of S3 dispatch to allow quick calculation of composite sets.  In line with design principle 4), lazy and greedy evaluation with symbolic representation allow for composite sets to be created, inspected, and printed, without ever needing to be evaluated themselves.

# Software Availability

``set6`` is available on [GitHub](https://github.com/xoopR/set6) and [CRAN](https://CRAN.R-project.org/package=set6). It can either be installed from GitHub using the `devtools` [@packagedevtools] library or directly from CRAN with `install.packages`. The package uses the MIT open-source licence. Contributions, issues, feature requests, and general feedback can all be found and provided on the project [GitHub](https://github.com/xoopR/set6). Full tutorials and further details are available on the [project website](https://xoopR.github.io/set6/).

# References
