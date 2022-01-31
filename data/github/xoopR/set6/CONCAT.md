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
---
title: "set6"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(set6)
set.seed(42)
```

<img src="man/figures/logo.png" align="right" alt="" width="120" />

[![R CMD Check via {tic}](https://github.com/xoopR/set6/workflows/R%20CMD%20Check%20via%20%7Btic%7D/badge.svg?branch=master)](https://github.com/xoopR/set6/actions)
[![codecov](https://app.codecov.io/gh/xoopR/set6/branch/master/graph/badge.svg)](https://app.codecov.io/gh/xoopR/set6)
[![CodeFactor](https://www.codefactor.io/repository/github/xoopr/set6/badge)](https://www.codefactor.io/repository/github/xoopr/set6)
[![dependencies](https://tinyverse.netlify.com/badge/set6)](https://CRAN.R-project.org/package=set6)

[![Repo Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle Badge](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://img.shields.io/badge/lifecycle-stable-brightgreen)

[![CRAN Status Badge](https://www.r-pkg.org/badges/version-ago/set6)](https://cran.r-project.org/package=set6)
[![CRAN Checks](https://cranchecks.info/badges/worst/set6)](https://cran.r-project.org/web/checks/check_results_set6.html)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/set6)](https://cran.r-project.org/package=set6)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02598/status.svg)](https://doi.org/10.21105/joss.02598)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Gitter chat](https://badges.gitter.im/xoopR/set6.png)](https://gitter.im/xoopR/set6)

## What is set6?

`set6` is an R6 upgrade to the `sets` package in R that includes:

* Multi-dimensional sets
* Tuples
* Finite and infinite intervals
* Fuzzy sets and tuples
* Set operations including union, intersect, (asymmetric and symmetric) difference, and product
* Symbolic representation of infinite sets including common special sets such as the Reals and Integers
* ConditionalSets for defining sets according to logical conditions

## Installation

The current CRAN release can be installed with
```{r,eval=FALSE}
install.packages("set6")
```
Or for the latest stable build

```{r,eval=FALSE}
remotes::install_github("xoopR/set6")
```

## Main Features

### A Clear Inheritance Structure

```{r}
# Sets require elements to be unique and order doesn't matter
Set$new(1, 2, 1) == Set$new(1, 2)
Set$new(1, 2) == Set$new(2, 1)

# But tuples can enforce these restrictions
Tuple$new(1, 2, 1) != Tuple$new(1, 2)
Tuple$new(1, 2) != Tuple$new(2, 1)

# Fuzzy sets and tuples extend sets further
f = FuzzySet$new(1, 0, 2, 0.6, 3, 1)
f$inclusion(1)
f$inclusion(2)
f$inclusion(3)

# Symbolic intervals provide a clean way to represent infinite sets
Interval$new()
# Different closure types and classes are possible
Interval$new(1, 7, type = "(]") # half-open
Interval$new(1, 7, class = "integer") == Set$new(1:7)

# And SpecialSets inheriting from these intervals
Reals$new()
PosRationals$new()
```

### Set operations
```{r}
# Union
Set$new(1, 4, "a", "b") + Set$new(5)
Interval$new(1, 5) + FuzzyTuple$new(1, 0.6)

# Power
Set$new(1:5)^2
# A symbolic representation is also possible
setpower(Set$new(1:5), power = 2, simplify = FALSE)
Reals$new()^5

# Product
Set$new(1,2) * Set$new(5, 6)
Interval$new(1,5) * Tuple$new(3)

# Intersection
Set$new(1:5) & Set$new(4:10)
ConditionalSet$new(function(x) x == 0) & Set$new(-2:2)
Interval$new(1, 10) & Set$new(5:6)

# Difference
Interval$new(1, 10) - Set$new(5)
Set$new(1:5) - Set$new(2:3)
```

### Containedness and Comparators
```{r}
Interval$new(1, 10)$contains(5)
# check multiple elements
Interval$new(1, 10)$contains(8:12)
# only return TRUE if all are TRUE
Interval$new(1, 10)$contains(8:12, all = TRUE)
# decide whether open bounds should be included
Interval$new(1, 10, type = "()")$contains(10, bound = TRUE)
Interval$new(1, 10, type = "()")$contains(10, bound = TRUE)

Interval$new(1, 5, class = "numeric")$equals(Set$new(1:5))
Interval$new(1, 5, class = "integer")$equals(Set$new(1:5))

Set$new(1) == FuzzySet$new(1, 1)

# proper subsets
Set$new(1:3)$isSubset(Set$new(1), proper = TRUE)
Set$new(1) < Set$new(1:3)

# (non-proper) subsets
Set$new(1:3)$isSubset(Set$new(1:3), proper = FALSE)
Set$new(1:3) <= Set$new(1:3)

# multi-dimensional checks
x = PosReals$new()^2
x$contains(list(Tuple$new(1, 1), Tuple$new(-2, 3)))
```

## Usage

The primary use-cases of `set6` are:

1. **Upgrading sets** Extend the R `sets` package to R6, which allows for generalised `Set` objects with a clear inheritance structure. As well as adding features including symbolic representation of infinite sets, and cartesian products.
2. **Defining parameter interfaces** All objects inheriting from the `Set` parent class include methods `equals` and `contains`, which are used to check if two sets are equal or if elements lie in the given set. This makes `set6` ideal for parameter interfaces in which a range of values (possibly multi-dimensional or of mixed types) need to be defined.

## Short-term development plans

Whilst the `set6` API is stable, it is considered 'maturing', and therefore whilst there are no plans for major updates, these may still occur. There are a few features and refactoring we plan on implementing before we consider the package to be in its first complete version. These mainly include

* Finalising all methods and fields - some are missing or possibly inaccurate for some wrappers. For example the cardinality of `ComplementSet`s is imprecise at the moment.
* We are considering adding a `simplify` method to wrappers to reduce classes inheriting from `SetWrapper` to simpler sets. This allows users to perform operations with `simplify = FALSE` and then to change their mind.
* There are known bottlenecks that need to be fixed to massively improve speed and efficiency.
* Adding more tutorials to make the interface easier for beginners, especially people new to R6

At a later stage we may consider adding Venn diagrams for visualisation of sets and intervals, but this is very low priority.

## Similar Packages

* [sets](https://CRAN.R-project.org/package=sets) - The **sets** package uses S3 to define some symbolic representaton of mathematical sets,
tuple, intervals, and fuzzy variants. However the symbolic representation is not consistent throughout
the package, does not allow for clear inspection of set/interval elements, and there is no support for
multi-dimensional sets.

* [BaseSet](https://github.com/ropensci/BaseSet) - The **BaseSet** package focuses on storing and analysing
sets in a 'tidy' way, with more options for data storage in long and wide formats. The primary usage is
neat and efficient inspection of finite sets, there is no support for infinite sets or symbolic
representation.

## Contributing

As `set6` is in its early stages, contributions are very welcome. If you have any ideas for good features please open an issue or create a pull request. Otherwise bug reports are very appreciated if you stumble across any broken code, these can be posted to the [issue tracker](https://github.com/xoopR/set6/issues). For updates on `set6` follow/star this repo.

## Citing set6

If you use set6, please cite our [JOSS article](https://doi.org/10.21105/joss.02598):

@Article{set6,
  title = {set6: R6 Mathematical Sets Interface},
  author = {Raphael Sonabend and Franz J. Kiraly},
  journal = {Journal of Open Source Software},
  year = {2020},
  month = {nov},
  doi = {10.21105/joss.02598},
  url = {https://joss.theoj.org/papers/10.21105/joss.02598},
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_setpower.R
\name{setpower}
\alias{setpower}
\alias{^.Set}
\title{Power of a Set}
\usage{
setpower(x, power, simplify = FALSE, nest = FALSE)

\method{^}{Set}(x, power)
}
\arguments{
\item{x}{Set}

\item{power}{power to raise set to, if \code{"n"} then a variable dimension set is created,
see examples.`}

\item{simplify}{logical, if \code{TRUE} returns the result in its simplest (unwrapped) form, usually
a \code{Set}, otherwise an \code{ExponentSet}.}

\item{nest}{logical, if \code{FALSE} (default) returns the n-ary cartesian product, otherwise returns
the cartesian product applied n times.
\link{Set}s. See details and examples.}
}
\value{
An R6 object of class \code{Set} or \code{ExponentSet} inheriting from \code{ProductSet}.
}
\description{
A convenience wrapper for the n-ary cartesian product of a \code{Set} by itself,
possibly multiple times.
}
\details{
See the details of \link{setproduct} for a longer discussion on the use of the \code{nest}
argument, in particular with regards to n-ary cartesian products vs. 'standard' cartesian
products.
}
\examples{
# Power of a Set
setpower(Set$new(1, 2), 3, simplify = FALSE)
setpower(Set$new(1, 2), 3, simplify = TRUE)
Set$new(1, 2)^3

# Power of an interval
Interval$new(2, 5)^5
Reals$new()^3

# Use tuples for contains
(PosNaturals$new()^3)$contains(Tuple$new(1, 2, 3))

# Power of ConditionalSet is meaningless
ConditionalSet$new(function(x) TRUE)^2

# Power of FuzzySet
FuzzySet$new(1, 0.1, 2, 0.5)^2

# Variable length
x <- Interval$new(0, 1)^"n"
x$contains(Tuple$new(0))
x$contains(Tuple$new(0, 1))
x$contains(Tuple$new(0, 1, 0, 0, 1, 1, 0))
x$contains(list(Tuple$new(0, 2), Tuple$new(1, 1)))

}
\seealso{
Other operators: 
\code{\link{powerset}()},
\code{\link{setcomplement}()},
\code{\link{setintersect}()},
\code{\link{setproduct}()},
\code{\link{setsymdiff}()},
\code{\link{setunion}()}
}
\concept{operators}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testSet}
\alias{testSet}
\alias{checkSet}
\alias{assertSet}
\title{assert/check/test/Set}
\usage{
testSet(object, errormsg = "This is not an R6 Set object")

checkSet(object, errormsg = "This is not an R6 Set object")

assertSet(object, errormsg = "This is not an R6 Set object")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{Set}.
}
\examples{
testSet(Set$new(2, 3))
testSet(list(Set$new(2), Set$new(3)))
testSet(Tuple$new(2, 3))
testSet(Interval$new())
testSet(FuzzySet$new(2, 0.1))
testSet(FuzzyTuple$new(2, 0.1))
testSet(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testSetList}
\alias{testSetList}
\alias{checkSetList}
\alias{assertSetList}
\title{assert/check/test/SetList}
\usage{
testSetList(object, errormsg = "One or more items in the list are not Sets")

checkSetList(object, errormsg = "One or more items in the list are not Sets")

assertSetList(object, errormsg = "One or more items in the list are not Sets")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{SetList}.
}
\examples{
testSetList(Set$new(2, 3))
testSetList(list(Set$new(2), Set$new(3)))
testSetList(Tuple$new(2, 3))
testSetList(Interval$new())
testSetList(FuzzySet$new(2, 0.1))
testSetList(FuzzyTuple$new(2, 0.1))
testSetList(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set6-package.R
\docType{package}
\name{set6-package}
\alias{set6}
\alias{set6-package}
\title{set6: R6 Mathematical Sets Interface}
\description{
set6 upgrades the \code{{sets}} package to R6. Many forms of mathematical sets are implemented,
including (countably finite) sets, tuples, intervals (countably infinite or uncountable),
and fuzzy variants. Wrappers extend functionality by allowing symbolic representations of
complex operations on sets, including unions, (cartesian) products, exponentiation, and
differences (asymmetric and symmetric).
}
\details{
The main features of set6 are:

\itemize{
\item Object-oriented programming, which allows a clear inheritance structure for Sets, Intervals,
Tuples, and other variants.
\item Set operations and wrappers for both explicit and symbolic representations for algebra of sets.
\item Methods for assertions and comparison checks, including subsets, equality, and containedness.
}

To learn more about set6, start with the set6 vignette:

\code{vignette("set6", "set6")}

And for more advanced usage see the complete tutorials at

\href{https://github.com/xoopR/set6}{https://github.com/xoopR/set6}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://xoopR.github.io/set6/}
  \item \url{https://github.com/xoopR/set6}
  \item Report bugs at \url{https://github.com/xoopR/set6/issues}
}

}
\author{
\strong{Maintainer}: Raphael Sonabend \email{raphaelsonabend@gmail.com} (\href{https://orcid.org/0000-0001-9225-4654}{ORCID})

Authors:
\itemize{
  \item Franz Kiraly \email{f.kiraly@ucl.ac.uk}
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{PosReals}
\alias{PosReals}
\title{Set of Positive Real Numbers}
\description{
The mathematical set of positive real numbers,
defined as the union of the set of positive rationals and positive irrationals. i.e.
\deqn{I^+ \cup Q^+}{I+ U Q+}
where \eqn{I^+}{I+} is the set of positive irrationals and \eqn{Q^+}{Q+} is the set of positive rationals.
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{\link[set6:Reals]{set6::Reals}} -> \code{PosReals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PosReals$new()}}
\item \href{#method-clone}{\code{PosReals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{PosReals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PosReals$new(zero = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{zero}}{logical. If TRUE, zero is included in the set.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{PosReals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PosReals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testFuzzyTuple}
\alias{testFuzzyTuple}
\alias{checkFuzzyTuple}
\alias{assertFuzzyTuple}
\title{assert/check/test/FuzzyTuple}
\usage{
testFuzzyTuple(object, errormsg = "This is not an R6 FuzzyTuple object")

checkFuzzyTuple(object, errormsg = "This is not an R6 FuzzyTuple object")

assertFuzzyTuple(object, errormsg = "This is not an R6 FuzzyTuple object")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{FuzzyTuple}.
}
\examples{
testFuzzyTuple(Set$new(2, 3))
testFuzzyTuple(list(Set$new(2), Set$new(3)))
testFuzzyTuple(Tuple$new(2, 3))
testFuzzyTuple(Interval$new())
testFuzzyTuple(FuzzySet$new(2, 0.1))
testFuzzyTuple(FuzzyTuple$new(2, 0.1))
testFuzzyTuple(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_UniversalSet.R
\name{UniversalSet}
\alias{UniversalSet}
\title{Mathematical Universal Set}
\description{
The \code{UniversalSet} is defined as the \link{Set} containing all possible elements.
}
\details{
The Universal set is the default universe to all sets, and is the largest possible set.
The Universal set contains every single possible element. We denote the Universal set with \code{V}
instead of \code{U} to avoid confusion with the union symbol. The Universal set cardinality is set to
\code{Inf} where we assume \code{Inf} is greater than any \code{Aleph} or \code{Beth} numbers. The Universal set is
also responsible for a few set paradoxes, to resolve these we use the following results:

Let \eqn{V} be the universal set, \eqn{S} be any non-universal set, and \eqn{0} the empty set, then

\deqn{V \cup S = V}{V or S = V}
\deqn{V \cap S = S}{V and S = S}
\deqn{S - V = 0}
\deqn{V^n = V}
\deqn{P(V) = V}
}
\examples{
u <- UniversalSet$new()
print(u)
u$contains(c(1, letters, TRUE, Set$new()), all = TRUE)

## ------------------------------------------------
## Method `UniversalSet$equals`
## ------------------------------------------------

# Equals
Set$new(1,2)$equals(Set$new(5,6))
Set$new(1,2)$equals(Interval$new(1,2))
Set$new(1,2) == Interval$new(1,2, class = "integer")

# Not equal
!Set$new(1,2)$equals(Set$new(1,2))
Set$new(1,2) != Set$new(1,5)

## ------------------------------------------------
## Method `UniversalSet$isSubset`
## ------------------------------------------------

Set$new(1,2,3)$isSubset(Set$new(1,2), proper = TRUE)
Set$new(1,2) < Set$new(1,2,3) # proper subset

c(Set$new(1,2,3), Set$new(1)) < Set$new(1,2,3) # not proper
Set$new(1,2,3) <= Set$new(1,2,3) # proper

## ------------------------------------------------
## Method `UniversalSet$contains`
## ------------------------------------------------

s = Set$new(1:5)

# Simplest case
s$contains(4)
8 \%inset\% s

# Test if multiple elements lie in the set
s$contains(4:6, all = FALSE)
s$contains(4:6, all = TRUE)

# Check if a tuple lies in a Set of higher dimension
s2 = s * s
s2$contains(Tuple$new(2,1))
c(Tuple$new(2,1), Tuple$new(1,7), 2) \%inset\% s2
}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{UniversalSet}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{UniversalSet$new()}}
\item \href{#method-equals}{\code{UniversalSet$equals()}}
\item \href{#method-isSubset}{\code{UniversalSet$isSubset()}}
\item \href{#method-contains}{\code{UniversalSet$contains()}}
\item \href{#method-strprint}{\code{UniversalSet$strprint()}}
\item \href{#method-clone}{\code{UniversalSet$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{UniversalSet} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UniversalSet$new()}\if{html}{\out{</div>}}
}

\subsection{Details}{
The Universal set is the set containing every possible element.
}

\subsection{Returns}{
A new \code{UniversalSet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UniversalSet$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Equals
Set$new(1,2)$equals(Set$new(5,6))
Set$new(1,2)$equals(Interval$new(1,2))
Set$new(1,2) == Interval$new(1,2, class = "integer")

# Not equal
!Set$new(1,2)$equals(Set$new(1,2))
Set$new(1,2) != Set$new(1,5)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UniversalSet$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

When calling \verb{$isSubset} on objects inheriting from \link{Interval}, the method treats the interval as if
it is a \link{Set}, i.e. ordering and class are ignored. Use \verb{$isSubinterval} to test if one interval
is a subinterval of another.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}

Every \code{Set} is a subset of a \code{UniversalSet}. No \code{Set} is a super set of a \code{UniversalSet},
and only a \code{UniversalSet} is not a proper subset of a \code{UniversalSet}.
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Set$new(1,2,3)$isSubset(Set$new(1,2), proper = TRUE)
Set$new(1,2) < Set$new(1,2,3) # proper subset

c(Set$new(1,2,3), Set$new(1)) < Set$new(1,2,3) # not proper
Set$new(1,2,3) <= Set$new(1,2,3) # proper
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests to see if \code{x} is contained in the Set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UniversalSet$contains(x, all = FALSE, bound = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{ignored.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\code{x} can be of any type, including a Set itself. \code{x} should be a tuple if
checking to see if it lies within a set of dimension greater than one. To test for multiple \code{x}
at the same time, then provide these as a list.

If using the method directly, and not via one of the operators then the additional boolean
arguments \code{all} and \code{bound}. If \code{all = TRUE} then returns \code{TRUE} if all \code{x} are contained in the \code{Set}, otherwise
returns a vector of logicals. For \link{Interval}s, \code{bound} is used to specify if elements lying on the
(possibly open) boundary of the interval are considered contained (\code{bound = TRUE}) or not (\code{bound = FALSE}).
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all elements of \code{x} are contained in the \code{Set}, otherwise
\code{FALSE.} If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

The infix operator \verb{\%inset\%} is available to test if \code{x} is an element in the \code{Set},
see examples.

Every element is contained within the Universal set.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{s = Set$new(1:5)

# Simplest case
s$contains(4)
8 \%inset\% s

# Test if multiple elements lie in the set
s$contains(4:6, all = FALSE)
s$contains(4:6, all = TRUE)

# Check if a tuple lies in a Set of higher dimension
s2 = s * s
s2$contains(Tuple$new(2,1))
c(Tuple$new(2,1), Tuple$new(1,7), 2) \%inset\% s2
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UniversalSet$strprint(n = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UniversalSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{Rationals}
\alias{Rationals}
\title{Set of Rational Numbers}
\description{
The mathematical set of rational numbers,
defined as the set of numbers that can be written as a fraction of two integers. i.e.
\deqn{\\{\frac{p}{q} \ : \ p,q \ \in \ Z, \ q \ne 0 \\}}{p/q : p,q \epsilon Z, q != 0}
where \eqn{Z} is the set of integers.
}
\details{
The \verb{$contains} method does not work for the set of Rationals as it is notoriously
difficult/impossible to find an algorithm for determining if any given number is rational or not.
Furthermore, computers must truncate all irrational numbers to rational numbers.
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{Rationals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Rationals$new()}}
\item \href{#method-contains}{\code{Rationals$contains()}}
\item \href{#method-isSubset}{\code{Rationals$isSubset()}}
\item \href{#method-equals}{\code{Rationals$equals()}}
\item \href{#method-clone}{\code{Rationals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Rationals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Rationals$new(lower = -Inf, upper = Inf, type = "()")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lower}}{numeric. Where to start the set. Advised to ignore, used by child-classes.}

\item{\code{upper}}{numeric. Where to end the set. Advised to ignore, used by child-classes.}

\item{\code{type}}{character Set closure type. Advised to ignore, used by child-classes.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Rationals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Method not possible for Rationals.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Rationals$contains(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Method not possible for Rationals.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Rationals$isSubset(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Method not possible for Rationals.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Rationals$equals(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Rationals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testClosed}
\alias{testClosed}
\alias{checkClosed}
\alias{assertClosed}
\title{assert/check/test/Closed}
\usage{
testClosed(object, errormsg = "This is not a closed set")

checkClosed(object, errormsg = "This is not a closed set")

assertClosed(object, errormsg = "This is not a closed set")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is closed.
}
\examples{
testClosed(Interval$new(1, 10, type = "[]"))
testClosed(Interval$new(1, 10, type = "(]"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testClosedBelow}
\alias{testClosedBelow}
\alias{checkClosedBelow}
\alias{assertClosedBelow}
\title{assert/check/test/ClosedBelow}
\usage{
testClosedBelow(object, errormsg = "This is not a set closed below")

checkClosedBelow(object, errormsg = "This is not a set closed below")

assertClosedBelow(object, errormsg = "This is not a set closed below")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is closedbelow.
}
\examples{
testClosedBelow(Interval$new(1, 10, type = "[]"))
testClosedBelow(Interval$new(1, 10, type = "(]"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Universal.R
\name{Universal}
\alias{Universal}
\title{Mathematical Universal Set}
\description{
The \code{Universal} is defined as the \link{Set} containing all possible elements.
}
\details{
The Universal set is the default universe to all sets, and is the largest possible set.
The Universal set contains every single possible element. We denote the Universal set with \code{V}
instead of \code{U} to avoid confusion with the union symbol. The Universal set cardinality is set to
\code{Inf} where we assume \code{Inf} is greater than any \code{Aleph} or \code{Beth} numbers. The Universal set is
also responsible for a few set paradoxes, to resolve these we use the following results:

Let \eqn{V} be the universal set, \eqn{S} be any non-universal set, and \eqn{0} the empty set, then

\deqn{V \cup S = V}{V or S = V}
\deqn{V \cap S = S}{V and S = S}
\deqn{S - V = 0}
\deqn{V^n = V}
\deqn{P(V) = V}
}
\examples{
u <- Universal$new()
print(u)
u$contains(c(1, letters, TRUE, Set$new()), all = TRUE)

## ------------------------------------------------
## Method `Universal$equals`
## ------------------------------------------------

# Equals
Set$new(1,2)$equals(Set$new(5,6))
Set$new(1,2)$equals(Interval$new(1,2))
Set$new(1,2) == Interval$new(1,2, class = "integer")

# Not equal
!Set$new(1,2)$equals(Set$new(1,2))
Set$new(1,2) != Set$new(1,5)

## ------------------------------------------------
## Method `Universal$isSubset`
## ------------------------------------------------

Set$new(1,2,3)$isSubset(Set$new(1,2), proper = TRUE)
Set$new(1,2) < Set$new(1,2,3) # proper subset

c(Set$new(1,2,3), Set$new(1)) < Set$new(1,2,3) # not proper
Set$new(1,2,3) <= Set$new(1,2,3) # proper

## ------------------------------------------------
## Method `Universal$contains`
## ------------------------------------------------

s = Set$new(1:5)

# Simplest case
s$contains(4)
8 \%inset\% s

# Test if multiple elements lie in the set
s$contains(4:6, all = FALSE)
s$contains(4:6, all = TRUE)

# Check if a tuple lies in a Set of higher dimension
s2 = s * s
s2$contains(Tuple$new(2,1))
c(Tuple$new(2,1), Tuple$new(1,7), 2) \%inset\% s2
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}}
}
\concept{special sets}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{Universal}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Universal$new()}}
\item \href{#method-equals}{\code{Universal$equals()}}
\item \href{#method-isSubset}{\code{Universal$isSubset()}}
\item \href{#method-contains}{\code{Universal$contains()}}
\item \href{#method-strprint}{\code{Universal$strprint()}}
\item \href{#method-clone}{\code{Universal$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Universal} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Universal$new()}\if{html}{\out{</div>}}
}

\subsection{Details}{
The Universal set is the set containing every possible element.
}

\subsection{Returns}{
A new \code{Universal} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Universal$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Equals
Set$new(1,2)$equals(Set$new(5,6))
Set$new(1,2)$equals(Interval$new(1,2))
Set$new(1,2) == Interval$new(1,2, class = "integer")

# Not equal
!Set$new(1,2)$equals(Set$new(1,2))
Set$new(1,2) != Set$new(1,5)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Universal$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

When calling \verb{$isSubset} on objects inheriting from \link{Interval}, the method treats the interval as if
it is a \link{Set}, i.e. ordering and class are ignored. Use \verb{$isSubinterval} to test if one interval
is a subinterval of another.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}

Every \code{Set} is a subset of a \code{Universal}. No \code{Set} is a super set of a \code{Universal},
and only a \code{Universal} is not a proper subset of a \code{Universal}.
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Set$new(1,2,3)$isSubset(Set$new(1,2), proper = TRUE)
Set$new(1,2) < Set$new(1,2,3) # proper subset

c(Set$new(1,2,3), Set$new(1)) < Set$new(1,2,3) # not proper
Set$new(1,2,3) <= Set$new(1,2,3) # proper
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests to see if \code{x} is contained in the Set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Universal$contains(x, all = FALSE, bound = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{ignored.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\code{x} can be of any type, including a Set itself. \code{x} should be a tuple if
checking to see if it lies within a set of dimension greater than one. To test for multiple \code{x}
at the same time, then provide these as a list.

If using the method directly, and not via one of the operators then the additional boolean
arguments \code{all} and \code{bound}. If \code{all = TRUE} then returns \code{TRUE} if all \code{x} are contained in the \code{Set}, otherwise
returns a vector of logicals. For \link{Interval}s, \code{bound} is used to specify if elements lying on the
(possibly open) boundary of the interval are considered contained (\code{bound = TRUE}) or not (\code{bound = FALSE}).
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all elements of \code{x} are contained in the \code{Set}, otherwise
\code{FALSE.} If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

The infix operator \verb{\%inset\%} is available to test if \code{x} is an element in the \code{Set},
see examples.

Every element is contained within the Universal set.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{s = Set$new(1:5)

# Simplest case
s$contains(4)
8 \%inset\% s

# Test if multiple elements lie in the set
s$contains(4:6, all = FALSE)
s$contains(4:6, all = TRUE)

# Check if a tuple lies in a Set of higher dimension
s2 = s * s
s2$contains(Tuple$new(2,1))
c(Tuple$new(2,1), Tuple$new(1,7), 2) \%inset\% s2
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Universal$strprint(n = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Universal$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asInterval.R
\name{as.Interval}
\alias{as.Interval}
\alias{as.Interval.Set}
\alias{as.Interval.Interval}
\alias{as.Interval.list}
\alias{as.Interval.data.frame}
\alias{as.Interval.matrix}
\alias{as.Interval.numeric}
\alias{as.Interval.ConditionalSet}
\title{Coercion to R6 \code{Interval}}
\usage{
as.Interval(object)

\method{as.Interval}{Set}(object)

\method{as.Interval}{Interval}(object)

\method{as.Interval}{list}(object)

\method{as.Interval}{data.frame}(object)

\method{as.Interval}{matrix}(object)

\method{as.Interval}{numeric}(object)

\method{as.Interval}{ConditionalSet}(object)
}
\arguments{
\item{object}{object to coerce}
}
\description{
Coerces object to an R6 \link{Interval}.
}
\details{
\itemize{
\item \code{as.Interval.list/as.Interval.data.frame} - Assumes the \code{list}/\code{data.frame} has
named items/columns: \verb{lower, upper, type, class}.
\item \code{as.Interval.numeric} - If the \code{numeric} vector is a continuous interval with no breaks then
coerces to an \link{Interval} with: \verb{lower = min(object), upper = max(object), class = "integer"}.
Ordering is ignored.
\item \code{as.Interval.matrix} - Tries coercion via \link{as.Interval.numeric} on the first column of the
matrix.
\item \code{as.Interval.Set} - First tries coercion via \link{as.Interval.numeric}, if possible wraps result
in a \link{Set}.
\item \code{as.Interval.FuzzySet} - Tries coercion via \link{as.Interval.Set} on the support of the \link{FuzzySet}.
}
}
\seealso{
\link{Interval}

Other coercions: 
\code{\link{as.FuzzySet}()},
\code{\link{as.Set}()}
}
\concept{coercions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{PosNaturals}
\alias{PosNaturals}
\title{Set of Positive Natural Numbers}
\description{
The mathematical set of positive natural numbers, defined as the positive counting numbers. i.e.
\deqn{\\{1, 2, 3,...\\}}{0, 1, 2, 3,...}
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{\link[set6:Naturals]{set6::Naturals}} -> \code{PosNaturals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PosNaturals$new()}}
\item \href{#method-clone}{\code{PosNaturals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{PosNaturals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PosNaturals$new()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new \code{PosNaturals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PosNaturals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetWrapper_ComplementSet.R
\name{ComplementSet}
\alias{ComplementSet}
\title{Set of Complements}
\description{
ComplementSet class for symbolic complement of mathematical sets.
}
\details{
The purpose of this class is to provide a symbolic representation for the complement of sets that
cannot be represented in a simpler class. Whilst this is not an abstract class, it is not recommended to construct
this class directly but via the set operation methods.
}
\seealso{
Set operations: \link{setunion}, \link{setproduct}, \link{setpower}, \link{setcomplement}, \link{setsymdiff},  \link{powerset}, \link{setintersect}

Other wrappers: 
\code{\link{ExponentSet}},
\code{\link{PowersetSet}},
\code{\link{ProductSet}},
\code{\link{UnionSet}}
}
\concept{wrappers}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:SetWrapper]{set6::SetWrapper}} -> \code{ComplementSet}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{elements}}{Returns the elements in the object.}

\item{\code{length}}{Returns the number of elements in the object.}

\item{\code{addedSet}}{For the \code{ComplementSet} wrapper, \code{X-Y}, returns the set \code{X}.}

\item{\code{subtractedSet}}{For the \code{ComplementSet} wrapper, \code{X-Y}, returns the set \code{Y}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{ComplementSet$new()}}
\item \href{#method-strprint}{\code{ComplementSet$strprint()}}
\item \href{#method-contains}{\code{ComplementSet$contains()}}
\item \href{#method-clone}{\code{ComplementSet$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="equals">}\href{../../set6/html/SetWrapper.html#method-equals}{\code{set6::SetWrapper$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="isSubset">}\href{../../set6/html/SetWrapper.html#method-isSubset}{\code{set6::SetWrapper$isSubset()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{ComplementSet} object. It is not recommended to construct this class directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ComplementSet$new(addset, subtractset, lower = NULL, upper = NULL, type = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{addset}}{\link{Set} to be subtracted from.}

\item{\code{subtractset}}{\link{Set} to subtract.}

\item{\code{lower}}{lower bound of new object.}

\item{\code{upper}}{upper bound of new object.}

\item{\code{type}}{closure type of new object.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{ComplementSet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ComplementSet$strprint(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests if elements \code{x} are contained in \code{self}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ComplementSet$contains(x, all = FALSE, bound = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{logical}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all == TRUE} then returns \code{TRUE} if all \code{x} are contained in \code{self}, otherwise \code{FALSE}.
If \code{all == FALSE} returns a vector of logicals corresponding to the length of \code{x}, representing
if each is contained in \code{self}. If \code{bound == TRUE} then an element is contained in \code{self} if it
is on or within the (possibly-open) bounds of \code{self}, otherwise \code{TRUE} only if the element is within
\code{self} or the bounds are closed.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ComplementSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{NegRationals}
\alias{NegRationals}
\title{Set of Negative Rational Numbers}
\description{
The mathematical set of negative rational numbers,
defined as the set of numbers that can be written as a fraction of two integers and are non-positive. i.e.
\deqn{\\{\frac{p}{q} \ : \ p,q \ \in \ Z, \ p/q \le 0, \ q \ne 0\\}}{p/q : p,q \epsilon Z, p/q \le 0, q != 0}
where \eqn{Z} is the set of integers.
}
\details{
The \verb{$contains} method does not work for the set of Rationals as it is notoriously
difficult/impossible to find an algorithm for determining if any given number is rational or not.
Furthermore, computers must truncate all irrational numbers to rational numbers.
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{\link[set6:Rationals]{set6::Rationals}} -> \code{NegRationals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{NegRationals$new()}}
\item \href{#method-clone}{\code{NegRationals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Rationals" data-id="contains">}\href{../../set6/html/Rationals.html#method-contains}{\code{set6::Rationals$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Rationals" data-id="equals">}\href{../../set6/html/Rationals.html#method-equals}{\code{set6::Rationals$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Rationals" data-id="isSubset">}\href{../../set6/html/Rationals.html#method-isSubset}{\code{set6::Rationals$isSubset()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{NegRationals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NegRationals$new(zero = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{zero}}{logical. If TRUE, zero is included in the set.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{NegRationals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NegRationals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set.R
\name{Set}
\alias{Set}
\title{Mathematical Set}
\description{
A general Set object for mathematical sets. This also serves as the parent class to
intervals, tuples, and fuzzy variants.
}
\details{
Mathematical sets can loosely be thought of as a collection of objects of any kind. The Set class
is used for sets of finite elements, for infinite sets use \link{Interval}. These can be
expanded for fuzzy logic by using \link{FuzzySet}s. Elements in a set cannot be duplicated and ordering
of elements does not matter, \link{Tuple}s can be used if duplicates or ordering are required.
}
\examples{
# Set of integers
Set$new(1:5)

# Set of multiple types
Set$new("a", 5, Set$new(1))

# Each Set has properties and traits
s <- Set$new(1, 2, 3)
s$traits
s$properties

# Elements cannot be duplicated
Set$new(2, 2) == Set$new(2)

# Ordering does not matter
Set$new(1, 2) == Set$new(2, 1)

## ------------------------------------------------
## Method `Set$contains`
## ------------------------------------------------

s = Set$new(elements = 1:5)

# Simplest case
s$contains(4)
8 \%inset\% s

# Test if multiple elements lie in the set
s$contains(4:6, all = FALSE)
s$contains(4:6, all = TRUE)

# Check if a tuple lies in a Set of higher dimension
s2 = s * s
s2$contains(Tuple$new(2,1))
c(Tuple$new(2,1), Tuple$new(1,7), 2) \%inset\% s2

## ------------------------------------------------
## Method `Set$equals`
## ------------------------------------------------

# Equals
Set$new(1,2)$equals(Set$new(5,6))
Set$new(1,2)$equals(Interval$new(1,2))
Set$new(1,2) == Interval$new(1,2, class = "integer")

# Not equal
!Set$new(1,2)$equals(Set$new(1,2))
Set$new(1,2) != Set$new(1,5)

## ------------------------------------------------
## Method `Set$isSubset`
## ------------------------------------------------

Set$new(1,2,3)$isSubset(Set$new(1,2), proper = TRUE)
Set$new(1,2) < Set$new(1,2,3) # proper subset

c(Set$new(1,2,3), Set$new(1)) < Set$new(1,2,3) # not proper
Set$new(1,2,3) <= Set$new(1,2,3) # proper

## ------------------------------------------------
## Method `Set$add`
## ------------------------------------------------

Set$new(1,2)$add(3)$print()
Set$new(1,2,universe = Interval$new(1,3))$add(3)$print()
\dontrun{
# errors as 4 is not in [1,3]
Set$new(1,2,universe = Interval$new(1,3))$add(4)$print()
}
# coerced to complex
Set$new(0+1i, 2i, class = "complex")$add(4)$print()

# setunion vs. add
Set$new(1,2)$add(Interval$new(5,6))$print()
Set$new(1,2) + Interval$new(5,6)

## ------------------------------------------------
## Method `Set$remove`
## ------------------------------------------------

Set$new(1,2,3)$remove(1,2)$print()
Set$new(1,Set$new(1),2)$remove(Set$new(1))$print()
Interval$new(1,5)$remove(5)$print()
Interval$new(1,5)$remove(4)$print()

# setcomplement vs. remove
Set$new(1,2,3)$remove(Interval$new(5,7))$print()
Set$new(1,2,3) - Interval$new(5,7)

## ------------------------------------------------
## Method `Set$multiplicity`
## ------------------------------------------------

Set$new(1, 1, 2)$multiplicity()
Set$new(1, 1, 2)$multiplicity(1)
Set$new(1, 1, 2)$multiplicity(list(1, 2))
Tuple$new(1, 1, 2)$multiplicity(1)
Tuple$new(1, 1, 2)$multiplicity(2)
}
\seealso{
Other sets: 
\code{\link{ConditionalSet}},
\code{\link{FuzzyMultiset}},
\code{\link{FuzzySet}},
\code{\link{FuzzyTuple}},
\code{\link{Interval}},
\code{\link{Multiset}},
\code{\link{Tuple}}
}
\concept{sets}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{properties}}{Returns an object of class \code{Properties}, which lists the properties of the Set. Set
properties include:
\itemize{
\item \code{empty} - is the Set empty or does it contain elements?
\item \code{singleton} - is the Set a singleton? i.e. Does it contain only one element?
\item \code{cardinality} - number of elements in the Set
\item \code{countability} - One of: countably finite, countably infinite, uncountable
\item \code{closure} - One of: closed, open, half-open
}}

\item{\code{traits}}{List the traits of the Set. Set traits include:
\itemize{
\item \code{crisp} - is the Set crisp or fuzzy?
}}

\item{\code{type}}{Returns the type of the Set. One of: (), (], [), [], \{\}}

\item{\code{max}}{If the Set consists of numerics only then returns the maximum element in the Set. For open
or half-open sets, then the maximum is defined by
\deqn{upper - 1e-15}}

\item{\code{min}}{If the Set consists of numerics only then returns the minimum element in the Set. For open
or half-open sets, then the minimum is defined by
\deqn{lower + 1e-15}}

\item{\code{upper}}{If the Set consists of numerics only then returns the upper bound of the Set.}

\item{\code{lower}}{If the Set consists of numerics only then returns the lower bound of the Set.}

\item{\code{class}}{If all elements in the Set are the same class then returns that class, otherwise "ANY".}

\item{\code{elements}}{If the Set is finite then returns all elements in the Set as a \code{list}, otherwise "NA".}

\item{\code{universe}}{Returns the universe of the Set, i.e. the set of values that can be added to the Set.}

\item{\code{range}}{If the Set consists of numerics only then returns the range of the Set defined by
\deqn{upper - lower}}

\item{\code{length}}{If the Set is finite then returns the number of elements in the Set, otherwise Inf. See
the cardinality property for the type of infinity.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Set$new()}}
\item \href{#method-print}{\code{Set$print()}}
\item \href{#method-strprint}{\code{Set$strprint()}}
\item \href{#method-summary}{\code{Set$summary()}}
\item \href{#method-contains}{\code{Set$contains()}}
\item \href{#method-equals}{\code{Set$equals()}}
\item \href{#method-isSubset}{\code{Set$isSubset()}}
\item \href{#method-add}{\code{Set$add()}}
\item \href{#method-remove}{\code{Set$remove()}}
\item \href{#method-multiplicity}{\code{Set$multiplicity()}}
\item \href{#method-clone}{\code{Set$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Set} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$new(..., universe = Universal$new(), elements = NULL, class = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{\code{ANY} Elements can be of any class except \code{list}, as long as there is a unique
\code{as.character} coercion method available.}

\item{\code{universe}}{Set. Universe that the Set lives in, i.e. elements that could be added to
the Set. Default is \link{Universal}.}

\item{\code{elements}}{list. Alternative constructor that may be more efficient if passing objects
of multiple classes.}

\item{\code{class}}{character. Optional string naming a class that if supplied gives the set the
\code{typed} property. All elements will be coerced to this class and therefore there must be
a coercion method to this class available.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Set} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Prints a symbolic representation of the \code{Set}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$print(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
The function \code{\link[=useUnicode]{useUnicode()}} can be used to determine if unicode should be used when
printing the \code{Set}. Internally \code{print} first calls \code{strprint} to create a printable representation
of the Set.
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$strprint(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-summary"></a>}}
\if{latex}{\out{\hypertarget{method-summary}{}}}
\subsection{Method \code{summary()}}{
Summarises the \code{Set}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$summary(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
The function \code{\link[=useUnicode]{useUnicode()}} can be used to determine if unicode should be used when
printing the \code{Set}. Summarised details include the \code{Set} class, properties, and traits.
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests to see if \code{x} is contained in the Set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$contains(x, all = FALSE, bound = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{ignored, added for consistency.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\code{x} can be of any type, including a Set itself. \code{x} should be a tuple if
checking to see if it lies within a set of dimension greater than one. To test for multiple \code{x}
at the same time, then provide these as a list.

If \code{all = TRUE} then returns \code{TRUE} if all \code{x} are contained in the \code{Set}, otherwise
returns a vector of logicals.
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all elements of \code{x} are contained in the \code{Set}, otherwise
\code{FALSE.} If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

The infix operator \verb{\%inset\%} is available to test if \code{x} is an element in the \code{Set},
see examples.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{s = Set$new(elements = 1:5)

# Simplest case
s$contains(4)
8 \%inset\% s

# Test if multiple elements lie in the set
s$contains(4:6, all = FALSE)
s$contains(4:6, all = TRUE)

# Check if a tuple lies in a Set of higher dimension
s2 = s * s
s2$contains(Tuple$new(2,1))
c(Tuple$new(2,1), Tuple$new(1,7), 2) \%inset\% s2
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Two sets are equal if they contain the same elements. Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Equals
Set$new(1,2)$equals(Set$new(5,6))
Set$new(1,2)$equals(Interval$new(1,2))
Set$new(1,2) == Interval$new(1,2, class = "integer")

# Not equal
!Set$new(1,2)$equals(Set$new(1,2))
Set$new(1,2) != Set$new(1,5)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Set$new(1,2,3)$isSubset(Set$new(1,2), proper = TRUE)
Set$new(1,2) < Set$new(1,2,3) # proper subset

c(Set$new(1,2,3), Set$new(1)) < Set$new(1,2,3) # not proper
Set$new(1,2,3) <= Set$new(1,2,3) # proper
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add"></a>}}
\if{latex}{\out{\hypertarget{method-add}{}}}
\subsection{Method \code{add()}}{
Add elements to a set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$add(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{elements to add}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\verb{$add} is a wrapper around the \code{setunion} method with \code{setunion(self, Set$new(...))}.
Note a key difference is that any elements passed to \code{...} are first converted to a \code{Set}, this
important difference is illustrated in the examples by adding an \link{Interval} to a \code{Set}.

Additionally, \verb{$add} first coerces \code{...} to \verb{$class} if \code{self} is a typed-set (i.e. \verb{$class != "ANY"}),
and \verb{$add} checks if elements in \code{...} live in the universe of \code{self}.
}

\subsection{Returns}{
An object inheriting from \link{Set}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Set$new(1,2)$add(3)$print()
Set$new(1,2,universe = Interval$new(1,3))$add(3)$print()
\dontrun{
# errors as 4 is not in [1,3]
Set$new(1,2,universe = Interval$new(1,3))$add(4)$print()
}
# coerced to complex
Set$new(0+1i, 2i, class = "complex")$add(4)$print()

# setunion vs. add
Set$new(1,2)$add(Interval$new(5,6))$print()
Set$new(1,2) + Interval$new(5,6)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-remove"></a>}}
\if{latex}{\out{\hypertarget{method-remove}{}}}
\subsection{Method \code{remove()}}{
Remove elements from a set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$remove(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{elements to remove}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\verb{$remove} is a wrapper around the \code{setcomplement} method with
\code{setcomplement(self, Set$new(...))}. Note a key difference is that any elements passed to \code{...}
are first converted to a \code{Set}, this important difference is illustrated in the examples by
removing an \link{Interval} from a \code{Set}.
}

\subsection{Returns}{
If the complement cannot be simplified to a \code{Set} then a \link{ComplementSet} is returned
otherwise an object inheriting from \link{Set} is returned.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Set$new(1,2,3)$remove(1,2)$print()
Set$new(1,Set$new(1),2)$remove(Set$new(1))$print()
Interval$new(1,5)$remove(5)$print()
Interval$new(1,5)$remove(4)$print()

# setcomplement vs. remove
Set$new(1,2,3)$remove(Interval$new(5,7))$print()
Set$new(1,2,3) - Interval$new(5,7)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-multiplicity"></a>}}
\if{latex}{\out{\hypertarget{method-multiplicity}{}}}
\subsection{Method \code{multiplicity()}}{
Returns the number of times an element appears in a set,
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$multiplicity(element = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{element}}{element or list of elements in the \code{set}, if \code{NULL} returns multiplicity of all elements}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Value, or list of values, in R+.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Set$new(1, 1, 2)$multiplicity()
Set$new(1, 1, 2)$multiplicity(1)
Set$new(1, 1, 2)$multiplicity(list(1, 2))
Tuple$new(1, 1, 2)$multiplicity(1)
Tuple$new(1, 1, 2)$multiplicity(2)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Set$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testFinite}
\alias{testFinite}
\alias{checkFinite}
\alias{assertFinite}
\title{assert/check/test/Finite}
\usage{
testFinite(object, errormsg = "This is not finite")

checkFinite(object, errormsg = "This is not finite")

assertFinite(object, errormsg = "This is not finite")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is finite.
}
\examples{
testFinite(Interval$new(1, 10, class = "integer"))
testFinite(Interval$new(1, 10, class = "numeric"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetWrapper_ExponentSet.R
\name{ExponentSet}
\alias{ExponentSet}
\title{Set of Exponentiations}
\description{
ExponentSet class for symbolic exponentiation of mathematical sets.
}
\details{
The purpose of this class is to provide a symbolic representation for the exponentiation of sets that
cannot be represented in a simpler class. Whilst this is not an abstract class, it is not recommended to construct
this class directly but via the set operation methods.
}
\seealso{
Set operations: \link{setunion}, \link{setproduct}, \link{setpower}, \link{setcomplement}, \link{setsymdiff},  \link{powerset}, \link{setintersect}

Other wrappers: 
\code{\link{ComplementSet}},
\code{\link{PowersetSet}},
\code{\link{ProductSet}},
\code{\link{UnionSet}}
}
\concept{wrappers}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:SetWrapper]{set6::SetWrapper}} -> \code{\link[set6:ProductSet]{set6::ProductSet}} -> \code{ExponentSet}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{power}}{Returns the power that the wrapped set is raised to.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{ExponentSet$new()}}
\item \href{#method-strprint}{\code{ExponentSet$strprint()}}
\item \href{#method-contains}{\code{ExponentSet$contains()}}
\item \href{#method-clone}{\code{ExponentSet$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="equals">}\href{../../set6/html/SetWrapper.html#method-equals}{\code{set6::SetWrapper$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="isSubset">}\href{../../set6/html/SetWrapper.html#method-isSubset}{\code{set6::SetWrapper$isSubset()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{ExponentSet} object. It is not recommended to construct this
class directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ExponentSet$new(set, power)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{set}}{\link{Set} to wrap.}

\item{\code{power}}{numeric. Power to raise Set to.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{ExponentSet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ExponentSet$strprint(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests if elements \code{x} are contained in \code{self}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ExponentSet$contains(x, all = FALSE, bound = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{logical}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all == TRUE} then returns \code{TRUE} if all \code{x} are contained in \code{self},
otherwise \code{FALSE}. If \code{all == FALSE} returns a vector of logicals corresponding to the
length of \code{x}, representing if each is contained in \code{self}. If \code{bound == TRUE} then an
element is contained in \code{self} if it is on or within the (possibly-open) bounds of \code{self},
otherwise \code{TRUE} only if the element is within \code{self} or the bounds are closed.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ExponentSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{PosRationals}
\alias{PosRationals}
\title{Set of Positive Rational Numbers}
\description{
The mathematical set of positive rational numbers,
defined as the set of numbers that can be written as a fraction of two integers and are non-negative. i.e.
\deqn{\\{\frac{p}{q} \ : \ p,q \ \in \ Z, \ p/q \ge 0, \ q \ne 0\\}}{p/q : p,q \epsilon Z, p/q \ge 0, q != 0}
where \eqn{Z} is the set of integers.
}
\details{
The \verb{$contains} method does not work for the set of Rationals as it is notoriously
difficult/impossible to find an algorithm for determining if any given number is rational or not.
Furthermore, computers must truncate all irrational numbers to rational numbers.
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{\link[set6:Rationals]{set6::Rationals}} -> \code{PosRationals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PosRationals$new()}}
\item \href{#method-clone}{\code{PosRationals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Rationals" data-id="contains">}\href{../../set6/html/Rationals.html#method-contains}{\code{set6::Rationals$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Rationals" data-id="equals">}\href{../../set6/html/Rationals.html#method-equals}{\code{set6::Rationals$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Rationals" data-id="isSubset">}\href{../../set6/html/Rationals.html#method-isSubset}{\code{set6::Rationals$isSubset()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{PosRationals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PosRationals$new(zero = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{zero}}{logical. If TRUE, zero is included in the set.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{PosRationals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PosRationals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetWrapper_ProductSet.R
\name{ProductSet}
\alias{ProductSet}
\title{Set of Products}
\description{
ProductSet class for symbolic product of mathematical sets.
}
\details{
The purpose of this class is to provide a symbolic representation for the product of sets that
cannot be represented in a simpler class. Whilst this is not an abstract class, it is not recommended to construct
this class directly but via the set operation methods.
}
\seealso{
Set operations: \link{setunion}, \link{setproduct}, \link{setpower}, \link{setcomplement}, \link{setsymdiff},  \link{powerset}, \link{setintersect}

Other wrappers: 
\code{\link{ComplementSet}},
\code{\link{ExponentSet}},
\code{\link{PowersetSet}},
\code{\link{UnionSet}}
}
\concept{wrappers}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:SetWrapper]{set6::SetWrapper}} -> \code{ProductSet}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{length}}{Returns the number of elements in the object.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{ProductSet$new()}}
\item \href{#method-strprint}{\code{ProductSet$strprint()}}
\item \href{#method-contains}{\code{ProductSet$contains()}}
\item \href{#method-clone}{\code{ProductSet$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="equals">}\href{../../set6/html/SetWrapper.html#method-equals}{\code{set6::SetWrapper$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="isSubset">}\href{../../set6/html/SetWrapper.html#method-isSubset}{\code{set6::SetWrapper$isSubset()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{ProductSet} object. It is not recommended to construct this class directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ProductSet$new(
  setlist,
  lower = NULL,
  upper = NULL,
  type = NULL,
  cardinality = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{setlist}}{\code{list} of \link{Set}s to wrap.}

\item{\code{lower}}{lower bound of new object.}

\item{\code{upper}}{upper bound of new object.}

\item{\code{type}}{closure type of new object.}

\item{\code{cardinality}}{Either an integer, "Aleph0", or a beth number. If \code{NULL} then calculated automatically (recommended).}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{ProductSet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ProductSet$strprint(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests if elements \code{x} are contained in \code{self}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ProductSet$contains(x, all = FALSE, bound = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{logical}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all == TRUE} then returns \code{TRUE} if all \code{x} are contained in \code{self}, otherwise \code{FALSE}.
If \code{all == FALSE} returns a vector of logicals corresponding to the length of \code{x}, representing
if each is contained in \code{self}. If \code{bound == TRUE} then an element is contained in \code{self} if it
is on or within the (possibly-open) bounds of \code{self}, otherwise \code{TRUE} only if the element is within
\code{self} or the bounds are closed.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ProductSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_setunion.R
\name{setunion}
\alias{setunion}
\alias{+.Set}
\alias{|.Set}
\title{Union of Sets}
\usage{
setunion(..., simplify = TRUE)

\method{+}{Set}(x, y)

\method{|}{Set}(x, y)
}
\arguments{
\item{...}{\link{Set}s}

\item{simplify}{logical, if \code{TRUE} (default) returns the result in its simplest (unwrapped) form, usually a \code{Set},
otherwise a \code{UnionSet}.}

\item{x, y}{\link{Set}}
}
\value{
An object inheriting from \link{Set} containing the union of supplied sets.
}
\description{
Returns the union of objects inheriting from class \link{Set}.
}
\details{
The union of \eqn{N} sets, \eqn{X1, ..., XN}, is defined as the set of elements that exist
in one or more sets,
\deqn{U = \{x : x \epsilon X1 \quad or \quad x \epsilon X2 \quad or \quad...\quad or \quad x \epsilon XN\}}{U = {x : x \epsilon X1 or x \epsilon X2 or ... or x \epsilon XN}}

The union of multiple \link{ConditionalSet}s is given by combining their defining functions by an
'or', \code{|}, operator. See examples.

The union of fuzzy and crisp sets first coerces fuzzy sets to crisp sets by finding their support.
}
\examples{
# union of Sets

Set$new(-2:4) + Set$new(2:5)
setunion(Set$new(1, 4, "a"), Set$new("a", 6))
Set$new(1, 2) + Set$new("a", 1i) + Set$new(9)

# union of intervals

Interval$new(1, 10) + Interval$new(5, 15) + Interval$new(20, 30)
Interval$new(1, 2, type = "()") + Interval$new(2, 3, type = "(]")
Interval$new(1, 5, class = "integer") +
  Interval$new(2, 7, class = "integer")

# union of mixed types

Set$new(1:10) + Interval$new(5, 15)
Set$new(1:10) + Interval$new(5, 15, class = "integer")
Set$new(5, 7) | Tuple$new(6, 8, 7)

# union of FuzzySet
FuzzySet$new(1, 0.1, 2, 0.5) + Set$new(2:5)

# union of conditional sets

ConditionalSet$new(function(x, y) x >= y) +
  ConditionalSet$new(function(x, y) x == y) +
  ConditionalSet$new(function(x) x == 2)

# union of special sets
PosReals$new() + NegReals$new()
Set$new(-Inf, Inf) + Reals$new()
}
\seealso{
Other operators: 
\code{\link{powerset}()},
\code{\link{setcomplement}()},
\code{\link{setintersect}()},
\code{\link{setpower}()},
\code{\link{setproduct}()},
\code{\link{setsymdiff}()}
}
\concept{operators}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asFuzzySet.R
\name{as.FuzzySet}
\alias{as.FuzzySet}
\alias{as.FuzzySet.numeric}
\alias{as.FuzzySet.list}
\alias{as.FuzzySet.matrix}
\alias{as.FuzzySet.data.frame}
\alias{as.FuzzySet.Set}
\alias{as.FuzzySet.FuzzySet}
\alias{as.FuzzySet.Interval}
\alias{as.FuzzySet.ConditionalSet}
\alias{as.FuzzyTuple}
\alias{as.FuzzyTuple.numeric}
\alias{as.FuzzyTuple.list}
\alias{as.FuzzyTuple.matrix}
\alias{as.FuzzyTuple.data.frame}
\alias{as.FuzzyTuple.Set}
\alias{as.FuzzyTuple.FuzzySet}
\alias{as.FuzzyTuple.Interval}
\alias{as.FuzzyTuple.ConditionalSet}
\alias{as.FuzzyMultiset}
\alias{as.FuzzyMultiset.numeric}
\alias{as.FuzzyMultiset.list}
\alias{as.FuzzyMultiset.matrix}
\alias{as.FuzzyMultiset.data.frame}
\alias{as.FuzzyMultiset.Set}
\alias{as.FuzzyMultiset.FuzzySet}
\alias{as.FuzzyMultiset.Interval}
\alias{as.FuzzyMultiset.ConditionalSet}
\title{Coercion to R6 \code{FuzzySet}/\code{FuzzyTuple}}
\usage{
as.FuzzySet(object)

\method{as.FuzzySet}{numeric}(object)

\method{as.FuzzySet}{list}(object)

\method{as.FuzzySet}{matrix}(object)

\method{as.FuzzySet}{data.frame}(object)

\method{as.FuzzySet}{Set}(object)

\method{as.FuzzySet}{FuzzySet}(object)

\method{as.FuzzySet}{Interval}(object)

\method{as.FuzzySet}{ConditionalSet}(object)

as.FuzzyTuple(object)

\method{as.FuzzyTuple}{numeric}(object)

\method{as.FuzzyTuple}{list}(object)

\method{as.FuzzyTuple}{matrix}(object)

\method{as.FuzzyTuple}{data.frame}(object)

\method{as.FuzzyTuple}{Set}(object)

\method{as.FuzzyTuple}{FuzzySet}(object)

\method{as.FuzzyTuple}{Interval}(object)

\method{as.FuzzyTuple}{ConditionalSet}(object)

as.FuzzyMultiset(object)

\method{as.FuzzyMultiset}{numeric}(object)

\method{as.FuzzyMultiset}{list}(object)

\method{as.FuzzyMultiset}{matrix}(object)

\method{as.FuzzyMultiset}{data.frame}(object)

\method{as.FuzzyMultiset}{Set}(object)

\method{as.FuzzyMultiset}{FuzzySet}(object)

\method{as.FuzzyMultiset}{Interval}(object)

\method{as.FuzzyMultiset}{ConditionalSet}(object)
}
\arguments{
\item{object}{object to coerce}
}
\description{
Coerces object to an R6 \link{FuzzySet}/\link{FuzzyTuple}
}
\details{
\itemize{
\item \code{as.FuzzySet.list} - Assumes \code{list} has two items, named \code{elements} and \code{membership},
and that they are ordered to be corresponding.
\item \code{as.FuzzySet.matrix} - Assumes first column corresponds to \code{elements} and second column
corresponds to their respective \code{membership}.
\item \code{as.FuzzySet.data.frame} - First checks to see if one column is called \code{elements} and the
other is called \code{membership}. If not then uses \code{as.FuzzySet.matrix}.
\item \code{as.FuzzySet.Set} - Creates a \link{FuzzySet} by assuming \link{Set} elements all have \code{membership}
equal to \eqn{1}.
\item \code{as.FuzzySet.Interval} - First tries coercion via \link{as.Set.Interval} then uses
\link{as.FuzzySet.Set}.
}
}
\seealso{
\link{FuzzySet} \link{FuzzyTuple}

Other coercions: 
\code{\link{as.Interval}()},
\code{\link{as.Set}()}
}
\concept{coercions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testClosedAbove}
\alias{testClosedAbove}
\alias{checkClosedAbove}
\alias{assertClosedAbove}
\title{assert/check/test/ClosedAbove}
\usage{
testClosedAbove(object, errormsg = "This is not a set closed above")

checkClosedAbove(object, errormsg = "This is not a set closed above")

assertClosedAbove(object, errormsg = "This is not a set closed above")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is closedabove.
}
\examples{
testClosedAbove(Interval$new(1, 10, type = "[]"))
testClosedAbove(Interval$new(1, 10, type = "[)"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testEmpty}
\alias{testEmpty}
\alias{checkEmpty}
\alias{assertEmpty}
\title{assert/check/test/Empty}
\usage{
testEmpty(object, errormsg = "This is not an empty set")

checkEmpty(object, errormsg = "This is not an empty set")

assertEmpty(object, errormsg = "This is not an empty set")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is empty.
}
\examples{
testEmpty(Set$new())
testEmpty(Set$new(1))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_FuzzySet.R
\name{FuzzySet}
\alias{FuzzySet}
\title{Mathematical Fuzzy Set}
\description{
A general FuzzySet object for mathematical fuzzy sets, inheriting from \code{Set}.
}
\details{
Fuzzy sets generalise standard mathematical sets to allow for fuzzy relationships. Whereas a
standard, or crisp, set assumes that an element is either in a set or not, a fuzzy set allows
an element to be in a set to a particular degree, known as the membership function, which
quantifies the inclusion of an element by a number in [0, 1]. Thus a (crisp) set is a
fuzzy set where all elements have a membership equal to \eqn{1}. Similarly to \link{Set}s, elements
must be unique and the ordering does not matter, to establish order and non-unique elements,
\link{FuzzyTuple}s can be used.
}
\examples{
# Different constructors
FuzzySet$new(1, 0.5, 2, 1, 3, 0)
FuzzySet$new(elements = 1:3, membership = c(0.5, 1, 0))

# Crisp sets are a special case FuzzySet
# Note membership defaults to full membership
FuzzySet$new(elements = 1:5) == Set$new(1:5)

f <- FuzzySet$new(1, 0.2, 2, 1, 3, 0)
f$membership()
f$alphaCut(0.3)
f$core()
f$inclusion(0)
f$membership(0)
f$membership(1)

## ------------------------------------------------
## Method `FuzzySet$membership`
## ------------------------------------------------

f = FuzzySet$new(1, 0.1, 2, 0.5, 3, 1)
f$membership()
f$membership(2)
f$membership(list(1, 2))

## ------------------------------------------------
## Method `FuzzySet$alphaCut`
## ------------------------------------------------

f = FuzzySet$new(1, 0.1, 2, 0.5, 3, 1)
# Alpha-cut
f$alphaCut(0.5)

# Strong alpha-cut
f$alphaCut(0.5, strong = TRUE)

# Create a set from the alpha-cut
f$alphaCut(0.5, create = TRUE)

## ------------------------------------------------
## Method `FuzzySet$support`
## ------------------------------------------------

f = FuzzySet$new(0.1, 0, 1, 0.1, 2, 0.5, 3, 1)
f$support()
f$support(TRUE)

## ------------------------------------------------
## Method `FuzzySet$core`
## ------------------------------------------------

f = FuzzySet$new(0.1, 0, 1, 0.1, 2, 0.5, 3, 1)
f$core()
f$core(TRUE)

## ------------------------------------------------
## Method `FuzzySet$inclusion`
## ------------------------------------------------

f = FuzzySet$new(0.1, 0, 1, 0.1, 2, 0.5, 3, 1)
f$inclusion(0.1)
f$inclusion(1)
f$inclusion(3)
}
\seealso{
Other sets: 
\code{\link{ConditionalSet}},
\code{\link{FuzzyMultiset}},
\code{\link{FuzzyTuple}},
\code{\link{Interval}},
\code{\link{Multiset}},
\code{\link{Set}},
\code{\link{Tuple}}
}
\concept{sets}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{FuzzySet}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{FuzzySet$new()}}
\item \href{#method-strprint}{\code{FuzzySet$strprint()}}
\item \href{#method-membership}{\code{FuzzySet$membership()}}
\item \href{#method-alphaCut}{\code{FuzzySet$alphaCut()}}
\item \href{#method-support}{\code{FuzzySet$support()}}
\item \href{#method-core}{\code{FuzzySet$core()}}
\item \href{#method-inclusion}{\code{FuzzySet$inclusion()}}
\item \href{#method-equals}{\code{FuzzySet$equals()}}
\item \href{#method-isSubset}{\code{FuzzySet$isSubset()}}
\item \href{#method-clone}{\code{FuzzySet$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="contains">}\href{../../set6/html/Set.html#method-contains}{\code{set6::Set$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{FuzzySet} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$new(
  ...,
  elements = NULL,
  membership = rep(1, length(elements)),
  class = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Alternating elements and membership, see details.}

\item{\code{elements}}{Elements in the set, see details.}

\item{\code{membership}}{Corresponding membership of the elements, see details.}

\item{\code{class}}{Optional string naming a class that if supplied gives the set the \code{typed} property.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\code{FuzzySet}s can be constructed in one of two ways, either by supplying the elements and their
membership in alternate order, or by providing a list of elements to \code{elements} and a list of
respective memberships to \code{membership}, see examples. If the \code{class} argument is non-\code{NULL},
then all elements will be coerced to the given class in construction, and if elements of a
different class are added these will either be rejected or coerced.
}

\subsection{Returns}{
A new \code{FuzzySet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$strprint(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-membership"></a>}}
\if{latex}{\out{\hypertarget{method-membership}{}}}
\subsection{Method \code{membership()}}{
Returns the membership, i.e. value in [0, 1], of either the given element(s)
or all elements in the fuzzy set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$membership(element = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{element}}{element or list of element in the \code{set}, if \code{NULL} returns membership of all elements}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
For \code{FuzzySet}s this is straightforward and returns the membership of the given element(s),
however in \code{FuzzyTuple}s and \code{FuzzyMultiset}s when an element may be duplicated, the function returns the membership of
all instances of the element.
}

\subsection{Returns}{
Value, or list of values, in [0, 1].
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{f = FuzzySet$new(1, 0.1, 2, 0.5, 3, 1)
f$membership()
f$membership(2)
f$membership(list(1, 2))
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-alphaCut"></a>}}
\if{latex}{\out{\hypertarget{method-alphaCut}{}}}
\subsection{Method \code{alphaCut()}}{
The alpha-cut of a fuzzy set is defined as the set
\deqn{A_\alpha = \{x \epsilon F | m \ge \alpha\}}{A_\alpha = {x \epsilon F | m \ge \alpha}}
where \eqn{x} is an element in the fuzzy set, \eqn{F}, and \eqn{m} is the corresponding membership.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$alphaCut(alpha, strong = FALSE, create = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{alpha}}{numeric in [0, 1] to determine which elements to return}

\item{\code{strong}}{logical, if \code{FALSE} (default) then includes elements greater than or equal to alpha, otherwise only strictly greater than}

\item{\code{create}}{logical, if \code{FALSE} (default) returns the elements in the alpha cut, otherwise returns a crisp set of the elements}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Elements in \link{FuzzySet} or a \link{Set} of the elements.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{f = FuzzySet$new(1, 0.1, 2, 0.5, 3, 1)
# Alpha-cut
f$alphaCut(0.5)

# Strong alpha-cut
f$alphaCut(0.5, strong = TRUE)

# Create a set from the alpha-cut
f$alphaCut(0.5, create = TRUE)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-support"></a>}}
\if{latex}{\out{\hypertarget{method-support}{}}}
\subsection{Method \code{support()}}{
The support of a fuzzy set is defined as the set of elements whose membership
is greater than zero, or the strong alpha-cut with \eqn{\alpha = 0},
\deqn{A_\alpha = \{x \epsilon F | m > 0\}}{A_\alpha = {x \epsilon F | m > 0}}
where \eqn{x} is an element in the fuzzy set, \eqn{F}, and \eqn{m} is the corresponding
membership.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$support(create = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{create}}{logical, if \code{FALSE} (default) returns the support elements, otherwise returns a \link{Set} of the support elements}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Support elements in fuzzy set or a \link{Set} of the support elements.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{f = FuzzySet$new(0.1, 0, 1, 0.1, 2, 0.5, 3, 1)
f$support()
f$support(TRUE)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-core"></a>}}
\if{latex}{\out{\hypertarget{method-core}{}}}
\subsection{Method \code{core()}}{
The core of a fuzzy set is defined as the set of elements whose membership is equal to one,
or the alpha-cut with \eqn{\alpha = 1},
\deqn{A_\alpha = \{x \epsilon F \ : \ m \ge 1\}}{A_\alpha = {x \epsilon F : m \ge 1}}
where \eqn{x} is an element in the fuzzy set, \eqn{F}, and \eqn{m} is the corresponding membership.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$core(create = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{create}}{logical, if \code{FALSE} (default) returns the core elements, otherwise returns a \link{Set} of the core elements}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Core elements in \link{FuzzySet} or a \link{Set} of the core elements.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{f = FuzzySet$new(0.1, 0, 1, 0.1, 2, 0.5, 3, 1)
f$core()
f$core(TRUE)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-inclusion"></a>}}
\if{latex}{\out{\hypertarget{method-inclusion}{}}}
\subsection{Method \code{inclusion()}}{
An element in a fuzzy set, with corresponding membership \eqn{m}, is:
\itemize{
\item Included - If \eqn{m = 1}
\item Partially Included - If \eqn{0 < m < 1}
\item Not Included - If \eqn{m = 0}
}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$inclusion(element)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{element}}{element or list of elements in fuzzy set for which to get the inclusion level}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
For \link{FuzzySet}s this is straightforward and returns the inclusion level of the given element(s),
however in \link{FuzzyTuple}s and \link{FuzzyMultiset}s when an element may be duplicated, the function returns the inclusion level of
all instances of the element.
}

\subsection{Returns}{
One of: "Included", "Partially Included", "Not Included"
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{f = FuzzySet$new(0.1, 0, 1, 0.1, 2, 0.5, 3, 1)
f$inclusion(0.1)
f$inclusion(1)
f$inclusion(3)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Two fuzzy sets are equal if they contain the same elements with the same memberships.
Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzySet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetWrapper_PowersetSet.R
\name{PowersetSet}
\alias{PowersetSet}
\title{Set of Powersets}
\description{
PowersetSet class for symbolic powerset of mathematical sets.
}
\details{
The purpose of this class is to provide a symbolic representation for the powerset of sets that
cannot be represented in a simpler class. Whilst this is not an abstract class, it is not recommended to construct
this class directly but via the set operation methods.
}
\seealso{
Set operations: \link{setunion}, \link{setproduct}, \link{setpower}, \link{setcomplement}, \link{setsymdiff},  \link{powerset}, \link{setintersect}

Other wrappers: 
\code{\link{ComplementSet}},
\code{\link{ExponentSet}},
\code{\link{ProductSet}},
\code{\link{UnionSet}}
}
\concept{wrappers}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:SetWrapper]{set6::SetWrapper}} -> \code{\link[set6:ProductSet]{set6::ProductSet}} -> \code{PowersetSet}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PowersetSet$new()}}
\item \href{#method-strprint}{\code{PowersetSet$strprint()}}
\item \href{#method-contains}{\code{PowersetSet$contains()}}
\item \href{#method-isSubset}{\code{PowersetSet$isSubset()}}
\item \href{#method-clone}{\code{PowersetSet$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="equals">}\href{../../set6/html/SetWrapper.html#method-equals}{\code{set6::SetWrapper$equals()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{PowersetSet} object. It is not recommended to construct this class directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PowersetSet$new(set)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{set}}{\link{Set} to wrap.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{PowersetSet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PowersetSet$strprint(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests if elements \code{x} are contained in \code{self}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PowersetSet$contains(x, all = FALSE, bound = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{logical}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all == TRUE} then returns \code{TRUE} if all \code{x} are contained in \code{self}, otherwise \code{FALSE}.
If \code{all == FALSE} returns a vector of logicals corresponding to the length of \code{x}, representing
if each is contained in \code{self}. If \code{bound == TRUE} then an element is contained in \code{self} if it
is on or within the (possibly-open) bounds of \code{self}, otherwise \code{TRUE} only if the element is within
\code{self} or the bounds are closed.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Tests if \code{x} is a (proper) subset of \code{self}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PowersetSet$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all == TRUE} then returns \code{TRUE} if all \code{x} are (proper) subsets of \code{self}, otherwise \code{FALSE}.
If \code{all == FALSE} returns a vector of logicals corresponding to the length of \code{x}, representing
if each is a (proper) subset of \code{self}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PowersetSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{Reals}
\alias{Reals}
\title{Set of Real Numbers}
\description{
The mathematical set of real numbers,
defined as the union of the set of rationals and irrationals. i.e.
\deqn{I \cup Q}{I U Q}
where \eqn{I} is the set of irrationals and \eqn{Q} is the set of rationals.
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{Reals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Reals$new()}}
\item \href{#method-clone}{\code{Reals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Reals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Reals$new(lower = -Inf, upper = Inf, type = "()")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lower}}{numeric. Where to start the set. Advised to ignore, used by child-classes.}

\item{\code{upper}}{numeric. Where to end the set. Advised to ignore, used by child-classes.}

\item{\code{type}}{character Set closure type. Advised to ignore, used by child-classes.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Reals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Reals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operators.R
\name{equals}
\alias{equals}
\alias{==.Set}
\alias{!=.Set}
\title{equals Operator}
\usage{
\method{==}{Set}(x, y)

\method{!=}{Set}(x, y)
}
\arguments{
\item{x, y}{\link{Set}}
}
\description{
Operator for \verb{$equals} methods. See \link{Set}\verb{$equals} for full details.
Operators can be used for:
\tabular{lll}{
\strong{Name} \tab \strong{Description} \tab \strong{Operator} \cr
Equal \tab \code{x} equals \code{y} \tab \code{==} \cr
Not Equal \tab \code{x} does not equal \code{y} \tab \code{!=} \cr

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operators.R
\name{isSubset}
\alias{isSubset}
\alias{<.Set}
\alias{<=.Set}
\alias{>.Set}
\alias{>=.Set}
\title{isSubset Operator}
\usage{
\method{<}{Set}(x, y)

\method{<=}{Set}(x, y)

\method{>}{Set}(x, y)

\method{>=}{Set}(x, y)
}
\arguments{
\item{x, y}{\link{Set}}
}
\description{
Operator for \verb{$isSubset} methods. See \link{Set}\verb{$isSubset} for full details.
Operators can be used for:
\tabular{lll}{
\strong{Name} \tab \strong{Description} \tab \strong{Operator} \cr
Subset \tab \code{x} is a subset of \code{y} \tab \code{x <= y} \cr
Proper Subset \tab \code{x} is a proper subset of \code{y} \tab \code{x < y} \cr
Superset \tab \code{x} is a superset of \code{y} \tab \code{x >= y} \cr
Proper Superset \tab \code{x} is a proper superset of \code{y} \tab \code{x > y} \cr
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval.R
\name{Interval}
\alias{Interval}
\title{Mathematical Finite or Infinite Interval}
\description{
A general Interval object for mathematical intervals, inheriting from \link{Set}. Intervals
may be open, closed, or half-open; as well as bounded above, below, or not at all.
}
\details{
The Interval class can be used for finite or infinite intervals, but often Sets will be preferred for
integer intervals over a finite continuous range.
}
\examples{
# Set of Reals
Interval$new()

# Set of Integers
Interval$new(class = "integer")

# Half-open interval
i <- Interval$new(1, 10, "(]")
i$contains(c(1, 10))
i$contains(c(1, 10), bound = TRUE)

# Equivalent Set and Interval
Set$new(1:5) == Interval$new(1, 5, class = "integer")

# SpecialSets can provide more efficient implementation
Interval$new() == ExtendedReals$new()
Interval$new(class = "integer", type = "()") == Integers$new()

## ------------------------------------------------
## Method `Interval$equals`
## ------------------------------------------------

Interval$new(1,5) == Interval$new(1,5)
Interval$new(1,5, class = "integer") != Interval$new(1,5,class="numeric")

## ------------------------------------------------
## Method `Interval$contains`
## ------------------------------------------------

s = Set$new(1:5)

# Simplest case
s$contains(4)
8 \%inset\% s

# Test if multiple elements lie in the set
s$contains(4:6, all = FALSE)
s$contains(4:6, all = TRUE)

# Check if a tuple lies in a Set of higher dimension
s2 = s * s
s2$contains(Tuple$new(2,1))
c(Tuple$new(2,1), Tuple$new(1,7), 2) \%inset\% s2

## ------------------------------------------------
## Method `Interval$isSubset`
## ------------------------------------------------

Interval$new(1,3) < Interval$new(1,5)
Set$new(1,3) < Interval$new(0,5)

## ------------------------------------------------
## Method `Interval$isSubinterval`
## ------------------------------------------------

Interval$new(1,3)$isSubset(Set$new(1,2)) # TRUE
Interval$new(1,3)$isSubset(Set$new(2, 1)) # TRUE
Interval$new(1,3, class = "integer")$isSubinterval(Set$new(1, 2)) # TRUE
Interval$new(1,3)$isSubinterval(Set$new(1, 2)) # FALSE
Interval$new(1,3)$isSubinterval(Set$new(2, 1)) # FALSE

Reals$new()$isSubset(Integers$new()) # TRUE
Reals$new()$isSubinterval(Integers$new()) # FALSE
}
\seealso{
Other sets: 
\code{\link{ConditionalSet}},
\code{\link{FuzzyMultiset}},
\code{\link{FuzzySet}},
\code{\link{FuzzyTuple}},
\code{\link{Multiset}},
\code{\link{Set}},
\code{\link{Tuple}}
}
\concept{sets}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{Interval}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{length}}{If the \code{Interval} is countably finite then returns the number of elements in the \code{Interval},
otherwise \code{Inf}. See the cardinality property for the type of infinity.}

\item{\code{elements}}{If the \code{Interval} is finite then returns all elements in the \code{Interval}, otherwise \code{NA}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Interval$new()}}
\item \href{#method-strprint}{\code{Interval$strprint()}}
\item \href{#method-equals}{\code{Interval$equals()}}
\item \href{#method-contains}{\code{Interval$contains()}}
\item \href{#method-isSubset}{\code{Interval$isSubset()}}
\item \href{#method-isSubinterval}{\code{Interval$isSubinterval()}}
\item \href{#method-clone}{\code{Interval$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Interval} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Interval$new(
  lower = -Inf,
  upper = Inf,
  type = c("[]", "(]", "[)", "()"),
  class = "numeric",
  universe = ExtendedReals$new()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lower}}{numeric. Lower limit of the interval.}

\item{\code{upper}}{numeric. Upper limit of the interval.}

\item{\code{type}}{character. One of: '()', '(]', '[)', '[]', which specifies if interval is open, left-open, right-open, or closed.}

\item{\code{class}}{character. One of: 'numeric', 'integer', which specifies if interval is over the Reals or Integers.}

\item{\code{universe}}{Set. Universe that the interval lives in, default \link{Reals}.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\code{Interval}s are constructed by specifying the \code{Interval} limits, the boundary type,
the class, and the possible universe. The \code{universe} differs from \code{class} as it is primarily used
for the \link{setcomplement} method. Whereas \code{class} specifies if the interval takes integers or
numerics, the \code{universe} specifies what range the interval could take.
}

\subsection{Returns}{
A new \code{Interval} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Interval$strprint(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{ignored, added for consistency.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Interval$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Two \code{Interval}s are equal if they have the same: class, type, and bounds.
Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Interval$new(1,5) == Interval$new(1,5)
Interval$new(1,5, class = "integer") != Interval$new(1,5,class="numeric")
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests to see if \code{x} is contained in the Set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Interval$contains(x, all = FALSE, bound = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{logical.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\code{x} can be of any type, including a Set itself. \code{x} should be a tuple if
checking to see if it lies within a set of dimension greater than one. To test for multiple \code{x}
at the same time, then provide these as a list.

If \code{all = TRUE} then returns \code{TRUE} if all \code{x} are contained in the \code{Set}, otherwise
returns a vector of logicals. For \link{Interval}s, \code{bound} is used to specify if elements lying on the
(possibly open) boundary of the interval are considered contained (\code{bound = TRUE}) or not (\code{bound = FALSE}).
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all elements of \code{x} are contained in the \code{Set}, otherwise
\code{FALSE.} If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

The infix operator \verb{\%inset\%} is available to test if \code{x} is an element in the \code{Set},
see examples.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{s = Set$new(1:5)

# Simplest case
s$contains(4)
8 \%inset\% s

# Test if multiple elements lie in the set
s$contains(4:6, all = FALSE)
s$contains(4:6, all = TRUE)

# Check if a tuple lies in a Set of higher dimension
s2 = s * s
s2$contains(Tuple$new(2,1))
c(Tuple$new(2,1), Tuple$new(1,7), 2) \%inset\% s2
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Interval$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

When calling \code{isSubset} on objects inheriting from \code{Interval}, the method treats the interval as if
it is a \link{Set}, i.e. ordering and class are ignored. Use \code{isSubinterval} to test if one interval
is a subinterval of another.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Interval$new(1,3) < Interval$new(1,5)
Set$new(1,3) < Interval$new(0,5)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubinterval"></a>}}
\if{latex}{\out{\hypertarget{method-isSubinterval}{}}}
\subsection{Method \code{isSubinterval()}}{
Test if one interval is a (proper) subinterval of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Interval$isSubinterval(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\code{Set} or \code{list}}

\item{\code{proper}}{If \code{TRUE} then tests if \code{x} is a proper subinterval (i.e. subinterval and not equal to)
of \code{self}, otherwise \code{FALSE} tests if \code{x} is a (non-proper) subinterval.}

\item{\code{all}}{If \code{TRUE} then returns \code{TRUE} if all \code{x} are subintervals, otherwise returns a vector of logicals.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If \code{x} is a \link{Set} then will be coerced to an \link{Interval} if possible. \verb{$isSubinterval} differs
from \verb{$isSubset} in that ordering and class are respected in \verb{$isSubinterval}. See examples for
a clearer illustration of the difference.
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Interval$new(1,3)$isSubset(Set$new(1,2)) # TRUE
Interval$new(1,3)$isSubset(Set$new(2, 1)) # TRUE
Interval$new(1,3, class = "integer")$isSubinterval(Set$new(1, 2)) # TRUE
Interval$new(1,3)$isSubinterval(Set$new(1, 2)) # FALSE
Interval$new(1,3)$isSubinterval(Set$new(2, 1)) # FALSE

Reals$new()$isSubset(Integers$new()) # TRUE
Reals$new()$isSubinterval(Integers$new()) # FALSE
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Interval$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testSubset}
\alias{testSubset}
\alias{checkSubset}
\alias{assertSubset}
\title{assert/check/test/Subset}
\usage{
testSubset(
  object,
  sets,
  proper = FALSE,
  errormsg = "sets are not subsets of the object"
)

checkSubset(
  object,
  sets,
  proper = FALSE,
  errormsg = "sets are not subsets of the object"
)

assertSubset(
  object,
  sets,
  proper = FALSE,
  errormsg = "sets are not subsets of the object"
)
}
\arguments{
\item{object}{object to test}

\item{sets}{sets to check}

\item{proper}{logical. If TRUE tests for proper subsets.}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if given sets are subsets of a set.
}
\examples{
testSubset(Set$new(1,2,3), Set$new(1,2))
testSubset(Set$new(1,2,3), Set$new(3,4))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/useUnicode.R
\name{useUnicode}
\alias{useUnicode}
\title{Get/Set Unicode Printing Method}
\usage{
useUnicode(use)
}
\arguments{
\item{use}{logical, if \code{TRUE} unicode will be used in printing, otherwise simpler character strings. If missing the current setting is returned.}
}
\description{
Change whether unicode symbols should be used when printing sets.
}
\details{
Using unicode symbols makes the printing of sets and properties 'prettier', however may not
work on all machines or versions of \code{R}. Therefore this function is used to decide whether
unicode representations should be used, or standard alpha-numeric and special characters.

By default \code{set6} starts with unicode printing turned on.
}
\examples{
current <- useUnicode()
useUnicode(TRUE)
useUnicode()
useUnicode(current)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operators.R
\name{contains}
\alias{contains}
\alias{\%inset\%}
\title{contains Operator}
\usage{
x \%inset\% y
}
\arguments{
\item{x, y}{\link{Set}}
}
\description{
Operator for \verb{$contains} methods. See \link{Set}\verb{$contains} for full details.
Operators can be used for:
\tabular{lll}{
\strong{Name} \tab \strong{Description} \tab \strong{Operator} \cr
Contains \tab \code{x} contains \code{y} \tab \code{y $inset$ x} \cr

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{Naturals}
\alias{Naturals}
\title{Set of Natural Numbers}
\description{
The mathematical set of natural numbers, defined as the counting numbers. i.e.
\deqn{\\{0, 1, 2,...\\}}{0, 1, 2,...}
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{Naturals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Naturals$new()}}
\item \href{#method-clone}{\code{Naturals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Naturals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Naturals$new(lower = 0)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lower}}{numeric. Where to start the set. Advised to ignore, used by child-classes.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Naturals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Naturals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_setsymdiff.R
\name{setsymdiff}
\alias{setsymdiff}
\alias{\%-\%}
\title{Symmetric Difference of Two Sets}
\usage{
setsymdiff(x, y, simplify = TRUE)

x \%-\% y
}
\arguments{
\item{x, y}{Set}

\item{simplify}{logical, if \code{TRUE} (default) returns the result in its simplest form, usually a \code{Set} or
\link{UnionSet}, otherwise a \code{ComplementSet}.}
}
\value{
An object inheriting from \code{Set} containing the symmetric difference of elements in both \code{x} and \code{y}.
}
\description{
Returns the symmetric difference of two objects inheriting from class \code{Set}.
}
\details{
The symmetric difference, aka disjunctive union, of two sets, \eqn{X, Y}, is defined as the set
of elements that exist in set \eqn{X} or in \eqn{Y} but not both,
\deqn{\{z : (z \epsilon X \cup z \epsilon Y) \\ \cap \\ \neg(z \epsilon X \cap z \epsilon Y)\}}{\{z : (z \epsilon X or z \epsilon Y) and !(z \epsilon X and z \epsilon Y)\}}

The symmetric difference can also be expressed as the union of two sets minus the intersection.
Therefore \code{setsymdiff} is written as a thin wrapper over these operations, so for two sets, \verb{A,B}: \cr
\code{A \%-\% B = (A | B) - (A & B)}.

The symmetric difference of fuzzy and crisp sets first coerces fuzzy sets to crisp sets by finding their support.
}
\examples{
# symmetrical difference compared to union and intersection
Set$new(1, 2, 3) \%-\% Set$new(3, 4)
(Set$new(1, 2, 3) | Set$new(3, 4)) - (Set$new(1, 2, 3) & Set$new(3, 4))

# ConditionalSets demonstrate the internal logic
ConditionalSet$new(function(x) x > 0) \%-\%
  ConditionalSet$new(function(y) y == 0)
}
\seealso{
Other operators: 
\code{\link{powerset}()},
\code{\link{setcomplement}()},
\code{\link{setintersect}()},
\code{\link{setpower}()},
\code{\link{setproduct}()},
\code{\link{setunion}()}
}
\concept{operators}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testMultiset}
\alias{testMultiset}
\alias{checkMultiset}
\alias{assertMultiset}
\title{assert/check/test/Multiset}
\usage{
testMultiset(object, errormsg = "This is not an R6 Multiset object")

checkMultiset(object, errormsg = "This is not an R6 Multiset object")

assertMultiset(object, errormsg = "This is not an R6 Multiset object")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{Multiset}.
}
\examples{
testMultiset(Set$new(2, 3))
testMultiset(list(Set$new(2), Set$new(3)))
testMultiset(Tuple$new(2, 3))
testMultiset(Interval$new())
testMultiset(FuzzySet$new(2, 0.1))
testMultiset(FuzzyTuple$new(2, 0.1))
testMultiset(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testCrisp}
\alias{testCrisp}
\alias{checkCrisp}
\alias{assertCrisp}
\title{assert/check/test/Crisp}
\usage{
testCrisp(object, errormsg = "This is not crisp.")

checkCrisp(object, errormsg = "This is not crisp.")

assertCrisp(object, errormsg = "This is not crisp.")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is crisp.
}
\examples{
testCrisp(Set$new(1))
testCrisp(FuzzySet$new(1, 0.5))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listSpecialSets.R
\name{listSpecialSets}
\alias{listSpecialSets}
\title{Lists Implemented R6 Special Sets}
\usage{
listSpecialSets(simplify = FALSE)
}
\arguments{
\item{simplify}{logical. If \code{FALSE} (default) returns data.frame of set name and symbol,
otherwise set names as characters.}
}
\value{
Either a list of characters (if \code{simplify} is \code{TRUE}) or a \code{data.frame} of \code{SpecialSet}s
and their traits.
}
\description{
Lists special sets that can be used in Set.
}
\examples{
listSpecialSets()
listSpecialSets(TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{ExtendedReals}
\alias{ExtendedReals}
\title{Set of Extended Real Numbers}
\description{
The mathematical set of extended real numbers,
defined as the union of the set of reals with \eqn{\pm\infty}{±\infty}. i.e.
\deqn{R \cup \\{-\infty, \infty\\}}{R U {-\infty, \infty}}
where \eqn{R} is the set of reals.
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{\link[set6:Reals]{set6::Reals}} -> \code{ExtendedReals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{ExtendedReals$new()}}
\item \href{#method-clone}{\code{ExtendedReals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{ExtendedReals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ExtendedReals$new()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new \code{ExtendedReals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ExtendedReals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testCountablyFinite}
\alias{testCountablyFinite}
\alias{checkCountablyFinite}
\alias{assertCountablyFinite}
\title{assert/check/test/CountablyFinite}
\usage{
testCountablyFinite(object, errormsg = "This is not a countably finite set")

checkCountablyFinite(object, errormsg = "This is not a countably finite set")

assertCountablyFinite(object, errormsg = "This is not a countably finite set")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is countablyfinite.
}
\examples{
testCountablyFinite(Set$new(1,2,3))
testCountablyFinite(Interval$new(1,10))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_setcomplement.R
\name{setcomplement}
\alias{setcomplement}
\alias{setcomplement.Set}
\alias{setcomplement.Interval}
\alias{setcomplement.FuzzySet}
\alias{setcomplement.ConditionalSet}
\alias{setcomplement.Reals}
\alias{setcomplement.Rationals}
\alias{setcomplement.Integers}
\alias{setcomplement.ComplementSet}
\alias{-.Set}
\title{Complement of Two Sets}
\usage{
setcomplement(x, y, simplify = TRUE)

\method{setcomplement}{Set}(x, y, simplify = TRUE)

\method{setcomplement}{Interval}(x, y, simplify = TRUE)

\method{setcomplement}{FuzzySet}(x, y, simplify = TRUE)

\method{setcomplement}{ConditionalSet}(x, y, simplify = TRUE)

\method{setcomplement}{Reals}(x, y, simplify = TRUE)

\method{setcomplement}{Rationals}(x, y, simplify = TRUE)

\method{setcomplement}{Integers}(x, y, simplify = TRUE)

\method{setcomplement}{ComplementSet}(x, y, simplify = TRUE)

\method{-}{Set}(x, y)
}
\arguments{
\item{x, y}{Set}

\item{simplify}{logical, if \code{TRUE} (default) returns the result in its simplest form, usually a
\code{Set} or \link{UnionSet}, otherwise a \code{ComplementSet}.}
}
\value{
An object inheriting from \code{Set} containing the set difference of elements in \code{x} and \code{y}.
}
\description{
Returns the set difference of two objects inheriting from class \code{Set}. If \code{y} is missing
then the complement of \code{x} from its universe is returned.
}
\details{
The difference of two sets, \eqn{X, Y}, is defined as the set of elements that exist
in set \eqn{X} and not \eqn{Y},
\deqn{X-Y = \{z : z \epsilon X \quad and \quad \neg(z \epsilon Y)\}}{X-Y = {z : z \epsilon X and !(z \epsilon Y)}}

The set difference of two \link{ConditionalSet}s is defined by combining their defining functions by a negated
'and', \verb{!&}, operator. See examples.

The complement of fuzzy and crisp sets first coerces fuzzy sets to crisp sets by finding their support.
}
\examples{
# absolute complement
setcomplement(Set$new(1, 2, 3, universe = Reals$new()))
setcomplement(Set$new(1, 2, universe = Set$new(1, 2, 3, 4, 5)))

# complement of two sets

Set$new(-2:4) - Set$new(2:5)
setcomplement(Set$new(1, 4, "a"), Set$new("a", 6))

# complement of two intervals

Interval$new(1, 10) - Interval$new(5, 15)
Interval$new(1, 10) - Interval$new(-15, 15)
Interval$new(1, 10) - Interval$new(-1, 2)

# complement of mixed set types

Set$new(1:10) - Interval$new(5, 15)
Set$new(5, 7) - Tuple$new(6, 8, 7)

# FuzzySet-Set returns a FuzzySet
FuzzySet$new(1, 0.1, 2, 0.5) - Set$new(2:5)
# Set-FuzzySet returns a Set
Set$new(2:5) - FuzzySet$new(1, 0.1, 2, 0.5)

# complement of conditional sets

ConditionalSet$new(function(x, y, simplify = TRUE) x >= y) -
  ConditionalSet$new(function(x, y, simplify = TRUE) x == y)

# complement of special sets
Reals$new() - NegReals$new()
Rationals$new() - PosRationals$new()
Integers$new() - PosIntegers$new()
}
\seealso{
Other operators: 
\code{\link{powerset}()},
\code{\link{setintersect}()},
\code{\link{setpower}()},
\code{\link{setproduct}()},
\code{\link{setsymdiff}()},
\code{\link{setunion}()}
}
\concept{operators}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testContains}
\alias{testContains}
\alias{checkContains}
\alias{assertContains}
\title{assert/check/test/Contains}
\usage{
testContains(
  object,
  elements,
  errormsg = "elements are not contained in the set"
)

checkContains(
  object,
  elements,
  errormsg = "elements are not contained in the set"
)

assertContains(
  object,
  elements,
  errormsg = "elements are not contained in the set"
)
}
\arguments{
\item{object}{object to test}

\item{elements}{elements to check}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if given elements are contained in a set.
}
\examples{
testContains(Set$new(1,2,3), c(1,2))
testContains(Set$new(1,2,3), c(3,4))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_setintersect.R
\name{setintersect}
\alias{setintersect}
\alias{setintersect.Interval}
\alias{setintersect.ConditionalSet}
\alias{setintersect.UnionSet}
\alias{setintersect.ComplementSet}
\alias{setintersect.ProductSet}
\alias{&.Set}
\title{Intersection of Two Sets}
\usage{
setintersect(x, y)

\method{setintersect}{Interval}(x, y)

\method{setintersect}{ConditionalSet}(x, y)

\method{setintersect}{UnionSet}(x, y)

\method{setintersect}{ComplementSet}(x, y)

\method{setintersect}{ProductSet}(x, y)

\method{&}{Set}(x, y)
}
\arguments{
\item{x, y}{Set}
}
\value{
A \code{Set} consisting of elements in both \code{x} and \code{y}.
}
\description{
Returns the intersection of two objects inheriting from class \code{Set}.
}
\details{
The intersection of two sets, \eqn{X, Y}, is defined as the set of elements that exist
in both sets,
\deqn{X \cap Y = \{z : z \epsilon X \quad and \quad z \epsilon Y\}}{{z : z \epsilon X and z \epsilon Y}}
In the case where no elements are common to either set, then the empty set is returned.

The intersection of two \link{ConditionalSet}s is defined by combining their defining functions by an
'and', \code{&}, operator. See examples.

The intersection of fuzzy and crisp sets first coerces fuzzy sets to crisp sets by finding their support.
}
\examples{
# intersection of two sets

Set$new(-2:4) & Set$new(2:5)
setintersect(Set$new(1, 4, "a"), Set$new("a", 6))
Set$new(1:4) & Set$new(5:7)

# intersection of two intervals

Interval$new(1, 10) & Interval$new(5, 15)
Interval$new(1, 2) & Interval$new(2, 3)
Interval$new(1, 5, class = "integer") &
  Interval$new(2, 7, class = "integer")

# intersection of mixed set types

Set$new(1:10) & Interval$new(5, 15)
Set$new(5, 7) & Tuple$new(6, 8, 7)

# Ignores membership of FuzzySet

FuzzySet$new(1, 0.1, 2, 0.5) & Set$new(2:5)

# intersection of conditional sets

ConditionalSet$new(function(x, y) x >= y) &
  ConditionalSet$new(function(x, y) x == y)
ConditionalSet$new(function(x) x == 2) &
  ConditionalSet$new(function(y) y == 3)

# But be careful not to make an empty set

ConditionalSet$new(function(x) x == 2) &
  ConditionalSet$new(function(x) x == 3)
}
\seealso{
Other operators: 
\code{\link{powerset}()},
\code{\link{setcomplement}()},
\code{\link{setpower}()},
\code{\link{setproduct}()},
\code{\link{setsymdiff}()},
\code{\link{setunion}()}
}
\concept{operators}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set6.news.R
\name{set6News}
\alias{set6News}
\title{Show set6 NEWS.md File}
\usage{
set6News()
}
\value{
NEWS.md in viewer.
}
\description{
Displays the contents of the NEWS.md file for viewing set6
release information.
}
\examples{
set6News()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_setproduct.R
\name{setproduct}
\alias{setproduct}
\alias{*.Set}
\title{Cartesian Product of Sets}
\usage{
setproduct(..., simplify = FALSE, nest = FALSE)

\method{*}{Set}(x, y)
}
\arguments{
\item{...}{\link{Set}s}

\item{simplify}{logical, if \code{TRUE} returns the result in its simplest (unwrapped) form, usually a \code{Set}
otherwise a \code{ProductSet}.}

\item{nest}{logical, if \code{FALSE} (default) then will treat any \link{ProductSet}s passed to \code{...} as unwrapped
\link{Set}s. See details and examples.}

\item{x, y}{\link{Set}}
}
\value{
Either an object of class \code{ProductSet} or an unwrapped object inheriting from \code{Set}.
}
\description{
Returns the cartesian product of objects inheriting from class \code{Set}.
}
\details{
The cartesian product of multiple sets, the 'n-ary Cartesian product', is often
implemented in programming languages as being identical to the cartesian product of two sets applied recursively.
However, for sets \eqn{X, Y, Z},
\deqn{XYZ \ne (XY)Z}{X × Y × Z != (X × Y) × Z}
This is accommodated with the \code{nest} argument. If \code{nest == TRUE} then \eqn{X*Y*Z == (X × Y) × Z}, i.e. the cartesian
product for two sets is applied recursively. If \code{nest == FALSE} then \eqn{X*Y*Z == (X × Y × Z)} and
the n-ary cartesian product is computed. As it appears the latter (n-ary product) is more common, \code{nest = FALSE}
is the default. The N-ary cartesian product of \eqn{N} sets, \eqn{X1,...,XN}, is defined as
\deqn{X1 × ... × XN = \\{(x1,...,xN) : x1 \epsilon X1 \cap ... \cap xN \epsilon XN\\}}{X1 × ... × XN = {(x1,...,xN) : x1 \epsilon X1 and ... and xN \epsilon xN}}
where \eqn{(x1,...,xN)} is a tuple.

The product of fuzzy and crisp sets first coerces fuzzy sets to crisp sets by finding their support.
}
\examples{
# difference between nesting
Set$new(1, 2) * Set$new(2, 3) * Set$new(4, 5)
setproduct(Set$new(1, 2) * Set$new(2, 3), Set$new(4, 5), nest = FALSE) # same as above
setproduct(Set$new(1, 2) * Set$new(2, 3), Set$new(4, 5), nest = TRUE)
unnest_set <- setproduct(Set$new(1, 2) * Set$new(2, 3), Set$new(4, 5), nest = FALSE)
nest_set <- setproduct(Set$new(1, 2) * Set$new(2, 3), Set$new(4, 5), nest = TRUE)
# note the difference when using contains
unnest_set$contains(Tuple$new(1, 3, 5))
nest_set$contains(Tuple$new(Tuple$new(1, 3), 5))

# product of two sets
Set$new(-2:4) * Set$new(2:5)
setproduct(Set$new(1, 4, "a"), Set$new("a", 6))
setproduct(Set$new(1, 4, "a"), Set$new("a", 6), simplify = TRUE)

# product of two intervals
Interval$new(1, 10) * Interval$new(5, 15)
Interval$new(1, 2, type = "()") * Interval$new(2, 3, type = "(]")
Interval$new(1, 5, class = "integer") *
  Interval$new(2, 7, class = "integer")

# product of mixed set types
Set$new(1:10) * Interval$new(5, 15)
Set$new(5, 7) * Tuple$new(6, 8, 7)
FuzzySet$new(1, 0.1) * Set$new(2)

# product of FuzzySet
FuzzySet$new(1, 0.1, 2, 0.5) * Set$new(2:5)

# product of conditional sets
ConditionalSet$new(function(x, y) x >= y) *
  ConditionalSet$new(function(x, y) x == y)

# product of special sets
PosReals$new() * NegReals$new()
}
\seealso{
Other operators: 
\code{\link{powerset}()},
\code{\link{setcomplement}()},
\code{\link{setintersect}()},
\code{\link{setpower}()},
\code{\link{setsymdiff}()},
\code{\link{setunion}()}
}
\concept{operators}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{PosIntegers}
\alias{PosIntegers}
\title{Set of Positive Integers}
\description{
The mathematical set of positive integers, defined as the set of positive whole numbers. i.e.
\deqn{\\{0, 1, 2, 3,...\\}}{0, 1, 2, 3,...}
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{\link[set6:Integers]{set6::Integers}} -> \code{PosIntegers}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PosIntegers$new()}}
\item \href{#method-clone}{\code{PosIntegers$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{PosIntegers} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PosIntegers$new(zero = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{zero}}{logical. If TRUE, zero is included in the set.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{PosIntegers} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PosIntegers$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operation_powerset.R
\name{powerset}
\alias{powerset}
\title{Calculate a Set's Powerset}
\usage{
powerset(x, simplify = FALSE)
}
\arguments{
\item{x}{\link{Set}}

\item{simplify}{logical, if \code{TRUE} then tries to simplify the result to a \code{Set} otherwise
creates an object of class \link{PowersetSet}.}
}
\value{
\link{Set}
}
\description{
Calculates and returns the powerset of a Set.
}
\details{
A powerset of a set, S, is defined as the set of all subsets of S, including S itself
and the empty set.
}
\examples{
# simplify = FALSE is default
powerset(Set$new(1, 2))
powerset(Set$new(1, 2), simplify = TRUE)

# powerset of intervals
powerset(Interval$new())

# powerset of powersets
powerset(powerset(Reals$new()))
powerset(powerset(Reals$new()))$properties$cardinality
}
\seealso{
Other operators: 
\code{\link{setcomplement}()},
\code{\link{setintersect}()},
\code{\link{setpower}()},
\code{\link{setproduct}()},
\code{\link{setsymdiff}()},
\code{\link{setunion}()}
}
\concept{operators}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testTuple}
\alias{testTuple}
\alias{checkTuple}
\alias{assertTuple}
\title{assert/check/test/Tuple}
\usage{
testTuple(object, errormsg = "This is not an R6 Tuple object")

checkTuple(object, errormsg = "This is not an R6 Tuple object")

assertTuple(object, errormsg = "This is not an R6 Tuple object")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{Tuple}.
}
\examples{
testTuple(Set$new(2, 3))
testTuple(list(Set$new(2), Set$new(3)))
testTuple(Tuple$new(2, 3))
testTuple(Interval$new())
testTuple(FuzzySet$new(2, 0.1))
testTuple(FuzzyTuple$new(2, 0.1))
testTuple(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Multiset.R
\name{Multiset}
\alias{Multiset}
\title{Mathematical Multiset}
\description{
A general Multiset object for mathematical Multisets, inheriting from \code{Set}.
}
\details{
Multisets are generalisations of sets that allow an element to be repeated. They can be thought
of as \link{Tuple}s without ordering.
}
\examples{
# Multiset of integers
Multiset$new(1:5)

# Multiset of multiple types
Multiset$new("a", 5, Set$new(1), Multiset$new(2))

# Each Multiset has properties and traits
t <- Multiset$new(1, 2, 3)
t$traits
t$properties

# Elements can be duplicated
Multiset$new(2, 2) != Multiset$new(2)

# Ordering does not matter
Multiset$new(1, 2) == Multiset$new(2, 1)

## ------------------------------------------------
## Method `Multiset$equals`
## ------------------------------------------------

Multiset$new(1,2) ==  Multiset$new(1,2)
Multiset$new(1,2) != Multiset$new(1,2)
Multiset$new(1,1) != Set$new(1,1)

## ------------------------------------------------
## Method `Multiset$isSubset`
## ------------------------------------------------

Multiset$new(1,2,3) < Multiset$new(1,2,3,4)
Multiset$new(1,3,2) < Multiset$new(1,2,3,4)
Multiset$new(1,3,2,4) <= Multiset$new(1,2,3,4)
}
\seealso{
Other sets: 
\code{\link{ConditionalSet}},
\code{\link{FuzzyMultiset}},
\code{\link{FuzzySet}},
\code{\link{FuzzyTuple}},
\code{\link{Interval}},
\code{\link{Set}},
\code{\link{Tuple}}
}
\concept{sets}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{Multiset}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-equals}{\code{Multiset$equals()}}
\item \href{#method-isSubset}{\code{Multiset$isSubset()}}
\item \href{#method-clone}{\code{Multiset$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="contains">}\href{../../set6/html/Set.html#method-contains}{\code{set6::Set$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="initialize">}\href{../../set6/html/Set.html#method-initialize}{\code{set6::Set$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="strprint">}\href{../../set6/html/Set.html#method-strprint}{\code{set6::Set$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Multiset$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
An object is equal to a Multiset if it contains all the same elements, and in the same order.
Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Multiset$new(1,2) ==  Multiset$new(1,2)
Multiset$new(1,2) != Multiset$new(1,2)
Multiset$new(1,1) != Set$new(1,1)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Multiset$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

When calling \verb{$isSubset} on objects inheriting from \link{Interval}, the method treats the interval as if
it is a \link{Set}, i.e. ordering and class are ignored. Use \verb{$isSubinterval} to test if one interval
is a subinterval of another.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}

An object is a (proper) subset of a Multiset if it contains all (some) of the same elements,
and in the same order.
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Multiset$new(1,2,3) < Multiset$new(1,2,3,4)
Multiset$new(1,3,2) < Multiset$new(1,2,3,4)
Multiset$new(1,3,2,4) <= Multiset$new(1,2,3,4)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Multiset$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testInterval}
\alias{testInterval}
\alias{checkInterval}
\alias{assertInterval}
\title{assert/check/test/Interval}
\usage{
testInterval(object, errormsg = "This is not an R6 Interval object")

checkInterval(object, errormsg = "This is not an R6 Interval object")

assertInterval(object, errormsg = "This is not an R6 Interval object")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{Interval}.
}
\examples{
testInterval(Set$new(2, 3))
testInterval(list(Set$new(2), Set$new(3)))
testInterval(Tuple$new(2, 3))
testInterval(Interval$new())
testInterval(FuzzySet$new(2, 0.1))
testInterval(FuzzyTuple$new(2, 0.1))
testInterval(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testFuzzyMultiset}
\alias{testFuzzyMultiset}
\alias{checkFuzzyMultiset}
\alias{assertFuzzyMultiset}
\title{assert/check/test/FuzzyMultiset}
\usage{
testFuzzyMultiset(object, errormsg = "This is not an R6 FuzzyMultiset object")

checkFuzzyMultiset(object, errormsg = "This is not an R6 FuzzyMultiset object")

assertFuzzyMultiset(
  object,
  errormsg = "This is not an R6 FuzzyMultiset object"
)
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{FuzzyMultiset}.
}
\examples{
testFuzzyMultiset(Set$new(2, 3))
testFuzzyMultiset(list(Set$new(2), Set$new(3)))
testFuzzyMultiset(Tuple$new(2, 3))
testFuzzyMultiset(Interval$new())
testFuzzyMultiset(FuzzySet$new(2, 0.1))
testFuzzyMultiset(FuzzyTuple$new(2, 0.1))
testFuzzyMultiset(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Properties.R
\name{Properties}
\alias{Properties}
\title{Set Properties Class}
\description{
Used to store the properties of a \link{Set}. Though this is not an abstract class,
it should never be constructed outside of the \link{Set} constructor.
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{closure}}{Returns the closure of the \code{Set}. One of "open", "half-open", or "closed."}

\item{\code{countability}}{Returns the countability of the \code{Set}. One of "countably finite", "countably infinite", or "uncountable".}

\item{\code{cardinality}}{Returns the cardinality of the \code{Set}. Either an integer if the \code{Set} is countably finite, Aleph0 if countably infinite, or a Beth number.}

\item{\code{empty}}{Returns if the \code{Set} is empty or not. \code{TRUE} if the Set cardinality is \code{0}, \code{FALSE} otherwise.}

\item{\code{singleton}}{Returns if the \code{Set} is a singleton or not. \code{TRUE} if the Set cardinality is \code{1}, \code{FALSE} otherwise.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Properties$new()}}
\item \href{#method-print}{\code{Properties$print()}}
\item \href{#method-strprint}{\code{Properties$strprint()}}
\item \href{#method-clone}{\code{Properties$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Creates a new \code{Properties} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Properties$new(closure = character(0), cardinality = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{closure}}{One of "open", "half-open", or "closed."}

\item{\code{cardinality}}{If non-\code{NULL} then either an integer, "Aleph0", or a Beth number.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Properties} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Prints the \code{Properties} list.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Properties$print()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Prints \code{Properties} list to console.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the \code{Properties}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Properties$strprint()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{list} of properties.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Properties$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testFuzzySet}
\alias{testFuzzySet}
\alias{checkFuzzySet}
\alias{assertFuzzySet}
\title{assert/check/test/FuzzySet}
\usage{
testFuzzySet(object, errormsg = "This is not an R6 FuzzySet object")

checkFuzzySet(object, errormsg = "This is not an R6 FuzzySet object")

assertFuzzySet(object, errormsg = "This is not an R6 FuzzySet object")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{FuzzySet}.
}
\examples{
testFuzzySet(Set$new(2, 3))
testFuzzySet(list(Set$new(2), Set$new(3)))
testFuzzySet(Tuple$new(2, 3))
testFuzzySet(Interval$new())
testFuzzySet(FuzzySet$new(2, 0.1))
testFuzzySet(FuzzyTuple$new(2, 0.1))
testFuzzySet(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{Integers}
\alias{Integers}
\title{Set of Integers}
\description{
The mathematical set of integers, defined as the set of whole numbers. i.e.
\deqn{\\{...,-3, -2, -1, 0, 1, 2, 3,...\\}}{...,-3, -2, -1, 0, 1, 2, 3,...}
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{Integers}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Integers$new()}}
\item \href{#method-clone}{\code{Integers$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Integers} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Integers$new(lower = -Inf, upper = Inf, type = "()")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lower}}{numeric. Where to start the set. Advised to ignore, used by child-classes.}

\item{\code{upper}}{numeric. Where to end the set. Advised to ignore, used by child-classes.}

\item{\code{type}}{character Set closure type. Advised to ignore, used by child-classes.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Integers} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Integers$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testConditionalSet}
\alias{testConditionalSet}
\alias{checkConditionalSet}
\alias{assertConditionalSet}
\title{assert/check/test/ConditionalSet}
\usage{
testConditionalSet(
  object,
  errormsg = "This is not an R6 ConditionalSet object"
)

checkConditionalSet(
  object,
  errormsg = "This is not an R6 ConditionalSet object"
)

assertConditionalSet(
  object,
  errormsg = "This is not an R6 ConditionalSet object"
)
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is an R6 \code{ConditionalSet}.
}
\examples{
testConditionalSet(Set$new(2, 3))
testConditionalSet(list(Set$new(2), Set$new(3)))
testConditionalSet(Tuple$new(2, 3))
testConditionalSet(Interval$new())
testConditionalSet(FuzzySet$new(2, 0.1))
testConditionalSet(FuzzyTuple$new(2, 0.1))
testConditionalSet(ConditionalSet$new(function(x) x == 0))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asSet.R
\name{as.Set}
\alias{as.Set}
\alias{as.Set.default}
\alias{as.Set.numeric}
\alias{as.Set.list}
\alias{as.Set.matrix}
\alias{as.Set.data.frame}
\alias{as.Set.Set}
\alias{as.Set.FuzzySet}
\alias{as.Set.Interval}
\alias{as.Set.ConditionalSet}
\alias{as.Tuple}
\alias{as.Tuple.default}
\alias{as.Tuple.numeric}
\alias{as.Tuple.list}
\alias{as.Tuple.matrix}
\alias{as.Tuple.data.frame}
\alias{as.Tuple.FuzzySet}
\alias{as.Tuple.Set}
\alias{as.Tuple.Interval}
\alias{as.Tuple.ConditionalSet}
\alias{as.Multiset}
\alias{as.Multiset.default}
\alias{as.Multiset.numeric}
\alias{as.Multiset.list}
\alias{as.Multiset.matrix}
\alias{as.Multiset.data.frame}
\alias{as.Multiset.FuzzySet}
\alias{as.Multiset.Set}
\alias{as.Multiset.Interval}
\alias{as.Multiset.ConditionalSet}
\title{Coercion to R6 \code{Set}/\code{Tuple}}
\usage{
as.Set(object)

\method{as.Set}{default}(object)

\method{as.Set}{numeric}(object)

\method{as.Set}{list}(object)

\method{as.Set}{matrix}(object)

\method{as.Set}{data.frame}(object)

\method{as.Set}{Set}(object)

\method{as.Set}{FuzzySet}(object)

\method{as.Set}{Interval}(object)

\method{as.Set}{ConditionalSet}(object)

as.Tuple(object)

\method{as.Tuple}{default}(object)

\method{as.Tuple}{numeric}(object)

\method{as.Tuple}{list}(object)

\method{as.Tuple}{matrix}(object)

\method{as.Tuple}{data.frame}(object)

\method{as.Tuple}{FuzzySet}(object)

\method{as.Tuple}{Set}(object)

\method{as.Tuple}{Interval}(object)

\method{as.Tuple}{ConditionalSet}(object)

as.Multiset(object)

\method{as.Multiset}{default}(object)

\method{as.Multiset}{numeric}(object)

\method{as.Multiset}{list}(object)

\method{as.Multiset}{matrix}(object)

\method{as.Multiset}{data.frame}(object)

\method{as.Multiset}{FuzzySet}(object)

\method{as.Multiset}{Set}(object)

\method{as.Multiset}{Interval}(object)

\method{as.Multiset}{ConditionalSet}(object)
}
\arguments{
\item{object}{object to coerce}
}
\description{
Coerces object to an R6 \link{Set}/\link{Tuple}
}
\details{
\itemize{
\item \code{as.Set.default} - Creates a \link{Set} using the object as the elements.
\item \code{as.Set.list} - Creates a \link{Set} for each element in \code{list}.
\item \code{as.Set.matrix/as.Set.data.frame} - Creates a \link{Set} for each column in \code{matrix/data.frame}.
\item \code{as.Set.FuzzySet} - Creates a \link{Set} from the support of the \link{FuzzySet}.
\item \code{as.Set.Interval} - If the interval has finite cardinality then creates a \link{Set} from the
\link{Interval} elements.
}
}
\seealso{
\link{Set} \link{Tuple}

Other coercions: 
\code{\link{as.FuzzySet}()},
\code{\link{as.Interval}()}
}
\concept{coercions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{testFuzzy}
\alias{testFuzzy}
\alias{checkFuzzy}
\alias{assertFuzzy}
\title{assert/check/test/Fuzzy}
\usage{
testFuzzy(object, errormsg = "This is not fuzzy.")

checkFuzzy(object, errormsg = "This is not fuzzy.")

assertFuzzy(object, errormsg = "This is not fuzzy.")
}
\arguments{
\item{object}{object to test}

\item{errormsg}{error message to overwrite default if check fails}
}
\value{
If check passes then \code{assert} returns \code{object} invisibly and \code{test}/\code{check}
return \code{TRUE}. If check fails, \code{assert} stops code with error, \code{check} returns
an error message as string, and \code{test} returns \code{FALSE}.
}
\description{
Validation checks to test if a given object is fuzzy.
}
\examples{
testFuzzy(FuzzySet$new(1, 0.5))
testFuzzy(Set$new(1))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_FuzzySet_FuzzyTuple.R
\name{FuzzyTuple}
\alias{FuzzyTuple}
\title{Mathematical Fuzzy Tuple}
\description{
A general FuzzyTuple object for mathematical fuzzy tuples, inheriting from \code{FuzzySet}.
}
\details{
Fuzzy tuples generalise standard mathematical tuples to allow for fuzzy relationships. Whereas a
standard, or crisp, tuple assumes that an element is either in a tuple or not, a fuzzy tuple allows
an element to be in a tuple to a particular degree, known as the membership function, which
quantifies the inclusion of an element by a number in [0, 1]. Thus a (crisp) tuple is a
fuzzy tuple where all elements have a membership equal to \eqn{1}. Similarly to \link{Tuple}s, elements
do not need to be unique and the ordering does matter, \link{FuzzySet}s are special cases where the ordering
does not matter and elements must be unique.
}
\examples{
# Different constructors
FuzzyTuple$new(1, 0.5, 2, 1, 3, 0)
FuzzyTuple$new(elements = 1:3, membership = c(0.5, 1, 0))

# Crisp sets are a special case FuzzyTuple
# Note membership defaults to full membership
FuzzyTuple$new(elements = 1:5) == Tuple$new(1:5)

f <- FuzzyTuple$new(1, 0.2, 2, 1, 3, 0)
f$membership()
f$alphaCut(0.3)
f$core()
f$inclusion(0)
f$membership(0)
f$membership(1)

# Elements can be duplicated, and with different memberships,
#  although this is not necessarily sensible.
FuzzyTuple$new(1, 0.1, 1, 1)

# More important is ordering.
FuzzyTuple$new(1, 0.1, 2, 0.2) != FuzzyTuple$new(2, 0.2, 1, 0.1)
FuzzySet$new(1, 0.1, 2, 0.2) == FuzzySet$new(2, 0.2, 1, 0.1)
}
\seealso{
Other sets: 
\code{\link{ConditionalSet}},
\code{\link{FuzzyMultiset}},
\code{\link{FuzzySet}},
\code{\link{Interval}},
\code{\link{Multiset}},
\code{\link{Set}},
\code{\link{Tuple}}
}
\concept{sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:FuzzySet]{set6::FuzzySet}} -> \code{FuzzyTuple}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-equals}{\code{FuzzyTuple$equals()}}
\item \href{#method-isSubset}{\code{FuzzyTuple$isSubset()}}
\item \href{#method-alphaCut}{\code{FuzzyTuple$alphaCut()}}
\item \href{#method-clone}{\code{FuzzyTuple$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="contains">}\href{../../set6/html/Set.html#method-contains}{\code{set6::Set$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="core">}\href{../../set6/html/FuzzySet.html#method-core}{\code{set6::FuzzySet$core()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="inclusion">}\href{../../set6/html/FuzzySet.html#method-inclusion}{\code{set6::FuzzySet$inclusion()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="initialize">}\href{../../set6/html/FuzzySet.html#method-initialize}{\code{set6::FuzzySet$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="membership">}\href{../../set6/html/FuzzySet.html#method-membership}{\code{set6::FuzzySet$membership()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="strprint">}\href{../../set6/html/FuzzySet.html#method-strprint}{\code{set6::FuzzySet$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="support">}\href{../../set6/html/FuzzySet.html#method-support}{\code{set6::FuzzySet$support()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzyTuple$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Two fuzzy sets are equal if they contain the same elements with the same memberships and
in the same order. Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzyTuple$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-alphaCut"></a>}}
\if{latex}{\out{\hypertarget{method-alphaCut}{}}}
\subsection{Method \code{alphaCut()}}{
The alpha-cut of a fuzzy set is defined as the set
\deqn{A_\alpha = \{x \epsilon F | m \ge \alpha\}}{A_\alpha = {x \epsilon F | m \ge \alpha}}
where \eqn{x} is an element in the fuzzy set, \eqn{F}, and \eqn{m} is the corresponding membership.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzyTuple$alphaCut(alpha, strong = FALSE, create = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{alpha}}{numeric in [0, 1] to determine which elements to return}

\item{\code{strong}}{logical, if \code{FALSE} (default) then includes elements greater than or equal to alpha, otherwise only strictly greater than}

\item{\code{create}}{logical, if \code{FALSE} (default) returns the elements in the alpha cut, otherwise returns a crisp set of the elements}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Elements in \link{FuzzyTuple} or a \link{Set} of the elements.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzyTuple$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{SpecialSet}
\alias{SpecialSet}
\title{Abstract Class for Special Sets}
\description{
The 'special sets' are the group of sets that are commonly used in mathematics
and are thus given their own names.
}
\details{
This is an abstract class and should not be constructed directly. Use \link{listSpecialSets}
to see the list of implemented special sets.
}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{SpecialSet}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{SpecialSet$new()}}
\item \href{#method-strprint}{\code{SpecialSet$strprint()}}
\item \href{#method-clone}{\code{SpecialSet$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\code{SpecialSet} is an abstract class, the constructor cannot be used directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SpecialSet$new(lower = -Inf, upper = Inf, type = "()", class = "numeric")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lower}}{defines the lower bound of the interval.}

\item{\code{upper}}{defines the upper bound of the interval.}

\item{\code{type}}{defines the interval closure type.}

\item{\code{class}}{defines the interval class.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SpecialSet$strprint(n = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{ignored, added for consistency.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SpecialSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Tuple.R
\name{Tuple}
\alias{Tuple}
\title{Mathematical Tuple}
\description{
A general Tuple object for mathematical tuples, inheriting from \code{Set}.
}
\details{
Tuples are similar to sets, except that they drop the constraint for elements to be unique, and
ordering in a tuple does matter. Tuples are useful for methods including \verb{$contains} that may
require non-unique elements. They are also the return type of the product of sets. See examples.
}
\examples{
# Tuple of integers
Tuple$new(1:5)

# Tuple of multiple types
Tuple$new("a", 5, Set$new(1), Tuple$new(2))

# Each Tuple has properties and traits
t <- Tuple$new(1, 2, 3)
t$traits
t$properties

# Elements can be duplicated
Tuple$new(2, 2) != Tuple$new(2)

# Ordering does matter
Tuple$new(1, 2) != Tuple$new(2, 1)

## ------------------------------------------------
## Method `Tuple$equals`
## ------------------------------------------------

Tuple$new(1,2) ==  Tuple$new(1,2)
Tuple$new(1,2) != Tuple$new(1,2)
Tuple$new(1,1) != Set$new(1,1)

## ------------------------------------------------
## Method `Tuple$isSubset`
## ------------------------------------------------

Tuple$new(1,2,3) < Tuple$new(1,2,3,4)
Tuple$new(1,3,2) < Tuple$new(1,2,3,4)
}
\seealso{
Other sets: 
\code{\link{ConditionalSet}},
\code{\link{FuzzyMultiset}},
\code{\link{FuzzySet}},
\code{\link{FuzzyTuple}},
\code{\link{Interval}},
\code{\link{Multiset}},
\code{\link{Set}}
}
\concept{sets}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{Tuple}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-equals}{\code{Tuple$equals()}}
\item \href{#method-isSubset}{\code{Tuple$isSubset()}}
\item \href{#method-clone}{\code{Tuple$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="contains">}\href{../../set6/html/Set.html#method-contains}{\code{set6::Set$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="initialize">}\href{../../set6/html/Set.html#method-initialize}{\code{set6::Set$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="strprint">}\href{../../set6/html/Set.html#method-strprint}{\code{set6::Set$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tuple$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
An object is equal to a Tuple if it contains all the same elements, and in the same order.
Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Tuple$new(1,2) ==  Tuple$new(1,2)
Tuple$new(1,2) != Tuple$new(1,2)
Tuple$new(1,1) != Set$new(1,1)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tuple$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

When calling \verb{$isSubset} on objects inheriting from \link{Interval}, the method treats the interval as if
it is a \link{Set}, i.e. ordering and class are ignored. Use \verb{$isSubinterval} to test if one interval
is a subinterval of another.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}

An object is a (proper) subset of a Tuple if it contains all (some) of the same elements,
and in the same order.
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Tuple$new(1,2,3) < Tuple$new(1,2,3,4)
Tuple$new(1,3,2) < Tuple$new(1,2,3,4)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Tuple$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Complex.R
\name{Complex}
\alias{Complex}
\title{Set of Complex Numbers}
\description{
The mathematical set of complex numbers,
defined as the the set of reals with possibly imaginary components. i.e.
\deqn{\\{a + bi \\ : \\ a,b \in R\\}}{{a + bi : a,b \epsilon R}}
where \eqn{R} is the set of reals.
}
\details{
There is no inherent ordering in the set of complex numbers, hence only the \code{contains}
method is implemented here.
}
\examples{

## ------------------------------------------------
## Method `Complex$equals`
## ------------------------------------------------

# Equals
Set$new(1,2)$equals(Set$new(5,6))
Set$new(1,2)$equals(Interval$new(1,2))
Set$new(1,2) == Interval$new(1,2, class = "integer")

# Not equal
!Set$new(1,2)$equals(Set$new(1,2))
Set$new(1,2) != Set$new(1,5)

## ------------------------------------------------
## Method `Complex$isSubset`
## ------------------------------------------------

Set$new(1,2,3)$isSubset(Set$new(1,2), proper = TRUE)
Set$new(1,2) < Set$new(1,2,3) # proper subset

c(Set$new(1,2,3), Set$new(1)) < Set$new(1,2,3) # not proper
Set$new(1,2,3) <= Set$new(1,2,3) # proper
}
\seealso{
Other special sets: 
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{Complex}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Complex$new()}}
\item \href{#method-contains}{\code{Complex$contains()}}
\item \href{#method-equals}{\code{Complex$equals()}}
\item \href{#method-isSubset}{\code{Complex$isSubset()}}
\item \href{#method-strprint}{\code{Complex$strprint()}}
\item \href{#method-clone}{\code{Complex$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Complex} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Complex$new()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new \code{Complex} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests to see if \code{x} is contained in the Set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Complex$contains(x, all = FALSE, bound = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{logical.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\code{x} can be of any type, including a Set itself. \code{x} should be a tuple if
checking to see if it lies within a set of dimension greater than one. To test for multiple \code{x}
at the same time, then provide these as a list.

If \code{all = TRUE} then returns \code{TRUE} if all \code{x} are contained in the \code{Set}, otherwise
returns a vector of logicals. For \link{Interval}s, \code{bound} is used to specify if elements lying on the
(possibly open) boundary of the interval are considered contained (\code{bound = TRUE}) or not (\code{bound = FALSE}).
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all elements of \code{x} are contained in the \code{Set}, otherwise
\code{FALSE.} If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

The infix operator \verb{\%inset\%} is available to test if \code{x} is an element in the \code{Set},
see examples.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Complex$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Equals
Set$new(1,2)$equals(Set$new(5,6))
Set$new(1,2)$equals(Interval$new(1,2))
Set$new(1,2) == Interval$new(1,2, class = "integer")

# Not equal
!Set$new(1,2)$equals(Set$new(1,2))
Set$new(1,2) != Set$new(1,5)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Complex$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

When calling \verb{$isSubset} on objects inheriting from \link{Interval}, the method treats the interval as if
it is a \link{Set}, i.e. ordering and class are ignored. Use \verb{$isSubinterval} to test if one interval
is a subinterval of another.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}

Every \code{Set} is a subset of a \code{Universal}. No \code{Set} is a super set of a \code{Universal},
and only a \code{Universal} is not a proper subset of a \code{Universal}.
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{Set$new(1,2,3)$isSubset(Set$new(1,2), proper = TRUE)
Set$new(1,2) < Set$new(1,2,3) # proper subset

c(Set$new(1,2,3), Set$new(1)) < Set$new(1,2,3) # not proper
Set$new(1,2,3) <= Set$new(1,2,3) # proper
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Complex$strprint(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Complex$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetWrapper.R
\name{SetWrapper}
\alias{SetWrapper}
\title{Abstract SetWrapper Class}
\description{
This class should not be constructed directly. Parent class to \code{SetWrapper}s.
}
\details{
Wrappers in set6 are utilised to facilitate lazy evaluation and symbolic representation.
Each operation has an associated wrapper that will be returned if \code{simplify = FALSE} or if the
result would be too complex to return as a simple \link{Set}. Wrappers have an identical interface
to \link{Set}. Their primary advantage lies in a neat representation of any set composition (the
result of an operation) and the ability to query the set contents without ever directly
evaluating the set elements.
}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{SetWrapper}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{wrappedSets}}{Returns the list of \code{Set}s that are wrapped in the given wrapper.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{SetWrapper$new()}}
\item \href{#method-equals}{\code{SetWrapper$equals()}}
\item \href{#method-isSubset}{\code{SetWrapper$isSubset()}}
\item \href{#method-clone}{\code{SetWrapper$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="contains">}\href{../../set6/html/Set.html#method-contains}{\code{set6::Set$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="strprint">}\href{../../set6/html/Set.html#method-strprint}{\code{set6::Set$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{SetWrapper} object. It is not recommended to construct this class directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SetWrapper$new(
  setlist,
  lower = NULL,
  upper = NULL,
  type = NULL,
  class = NULL,
  cardinality
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{setlist}}{List of \link{Set}s to wrap.}

\item{\code{lower}}{\link{Set}. Lower bound of wrapper.}

\item{\code{upper}}{\link{Set}. Upper bound of wrapper.}

\item{\code{type}}{character. Closure type of wrapper.}

\item{\code{class}}{character. Ignored.}

\item{\code{cardinality}}{character or integer. Cardinality of wrapper.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{SetWrapper} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if \code{x} is equal to \code{self}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SetWrapper$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all == TRUE} then returns \code{TRUE} if all \code{x} are equal to \code{self}, otherwise \code{FALSE}.
If \code{all == FALSE} returns a vector of logicals corresponding to the length of \code{x}, representing
if each is equal to \code{self}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Tests if \code{x} is a (proper) subset of \code{self}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SetWrapper$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all == TRUE} then returns \code{TRUE} if all \code{x} are (proper) subsets of \code{self}, otherwise \code{FALSE}.
If \code{all == FALSE} returns a vector of logicals corresponding to the length of \code{x}, representing
if each is a (proper) subset of \code{self}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SetWrapper$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{NegIntegers}
\alias{NegIntegers}
\title{Set of Negative Integers}
\description{
The mathematical set of negative integers, defined as the set of negative whole numbers. i.e.
\deqn{\\{...,-3, -2, -1, 0\\}}{...,-3, -2, -1, 0}
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{\link[set6:Integers]{set6::Integers}} -> \code{NegIntegers}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{NegIntegers$new()}}
\item \href{#method-clone}{\code{NegIntegers$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{NegIntegers} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NegIntegers$new(zero = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{zero}}{logical. If TRUE, zero is included in the set.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{NegIntegers} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NegIntegers$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_ConditionalSet.R
\name{ConditionalSet}
\alias{ConditionalSet}
\title{Mathematical Set of Conditions}
\description{
A mathematical set defined by one or more logical conditions.
}
\details{
Conditional sets are a useful tool for symbolically defining possibly infinite sets. They can be combined
using standard 'and', \code{&}, and 'or', \code{|}, operators.
}
\examples{
# Set of Positive Naturals
s <- ConditionalSet$new(function(x) TRUE, argclass = list(x = PosNaturals$new()))

## ------------------------------------------------
## Method `ConditionalSet$contains`
## ------------------------------------------------

# Set of positives
s = ConditionalSet$new(function(x) x > 0)
s$contains(list(1,-1))

# Set via equality
s = ConditionalSet$new(function(x, y) x + y == 2)
s$contains(list(Set$new(2, 0), Set$new(0, 2)))

# Tuples are recommended when using contains as they allow non-unique elements
s = ConditionalSet$new(function(x, y) x + y == 4)
\dontrun{
s$contains(Set$new(2, 2)) # Errors as Set$new(2,2) == Set$new(2)
}

# Set of Positive Naturals
s = ConditionalSet$new(function(x) TRUE, argclass = list(x = PosNaturals$new()))
s$contains(list(-2, 2))
}
\seealso{
Other sets: 
\code{\link{FuzzyMultiset}},
\code{\link{FuzzySet}},
\code{\link{FuzzyTuple}},
\code{\link{Interval}},
\code{\link{Multiset}},
\code{\link{Set}},
\code{\link{Tuple}}
}
\concept{sets}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{ConditionalSet}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{condition}}{Returns the condition defining the ConditionalSet.}

\item{\code{class}}{Returns \code{argclass}, see \verb{$new}.}

\item{\code{elements}}{Returns \code{NA}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{ConditionalSet$new()}}
\item \href{#method-contains}{\code{ConditionalSet$contains()}}
\item \href{#method-equals}{\code{ConditionalSet$equals()}}
\item \href{#method-strprint}{\code{ConditionalSet$strprint()}}
\item \href{#method-summary}{\code{ConditionalSet$summary()}}
\item \href{#method-isSubset}{\code{ConditionalSet$isSubset()}}
\item \href{#method-clone}{\code{ConditionalSet$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{ConditionalSet} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionalSet$new(condition = function(x) TRUE, argclass = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{condition}}{function. Defines the set, see details.}

\item{\code{argclass}}{list. Optional list of sets that the function arguments live in, see details.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
The \code{condition} should be given as a function that when evaluated returns
either \code{TRUE} or \code{FALSE}. Further constraints can be given by providing the universe of the
function arguments as \link{Set}s, if these are not given then \link{Universal} is assumed.
See examples. Defaults construct the Universal set.
}

\subsection{Returns}{
A new \code{ConditionalSet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests to see if \code{x} is contained in the Set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionalSet$contains(x, all = FALSE, bound = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{ignored, added for consistency.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
\code{x} can be of any type, including a Set itself. \code{x} should be a tuple if
checking to see if it lies within a set of dimension greater than one. To test for multiple \code{x}
at the same time, then provide these as a list.

If \code{all = TRUE} then returns \code{TRUE} if all \code{x} are contained in the \code{Set}, otherwise
returns a vector of logicals.

An element is contained in a \code{ConditionalSet} if it returns \code{TRUE} as an argument in the defining function.
For sets that are defined with a function that takes multiple arguments, a \link{Tuple} should be
passed to \code{x}.
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all elements of \code{x} are contained in the \code{Set}, otherwise
\code{FALSE.} If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.

The infix operator \verb{\%inset\%} is available to test if \code{x} is an element in the \code{Set},
see examples.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Set of positives
s = ConditionalSet$new(function(x) x > 0)
s$contains(list(1,-1))

# Set via equality
s = ConditionalSet$new(function(x, y) x + y == 2)
s$contains(list(Set$new(2, 0), Set$new(0, 2)))

# Tuples are recommended when using contains as they allow non-unique elements
s = ConditionalSet$new(function(x, y) x + y == 4)
\dontrun{
s$contains(Set$new(2, 2)) # Errors as Set$new(2,2) == Set$new(2)
}

# Set of Positive Naturals
s = ConditionalSet$new(function(x) TRUE, argclass = list(x = PosNaturals$new()))
s$contains(list(-2, 2))
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionalSet$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Two sets are equal if they contain the same elements. Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionalSet$strprint(n = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{ignored, added for consistency.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-summary"></a>}}
\if{latex}{\out{\hypertarget{method-summary}{}}}
\subsection{Method \code{summary()}}{
See \code{strprint}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionalSet$summary(n = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{ignored, added for consistency.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Currently undefined for \code{ConditionalSet}s.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionalSet$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{ignored, added for consistency.}

\item{\code{proper}}{ignored, added for consistency.}

\item{\code{all}}{ignored, added for consistency.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionalSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Logicals.R
\name{Logicals}
\alias{Logicals}
\title{Set of Logicals}
\description{
The \code{Logicals} is defined as the \link{Set} containing the elements \code{TRUE} and \code{FALSE}.
}
\examples{
l <- Logicals$new()
print(l)
l$contains(list(TRUE, 1, FALSE))
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{NegReals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{Logicals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Logicals$new()}}
\item \href{#method-clone}{\code{Logicals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="contains">}\href{../../set6/html/Set.html#method-contains}{\code{set6::Set$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="equals">}\href{../../set6/html/Set.html#method-equals}{\code{set6::Set$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="isSubset">}\href{../../set6/html/Set.html#method-isSubset}{\code{set6::Set$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="strprint">}\href{../../set6/html/Set.html#method-strprint}{\code{set6::Set$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Logicals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Logicals$new()}\if{html}{\out{</div>}}
}

\subsection{Details}{
The Logical set is the set containing \code{TRUE} and \code{FALSE}.
}

\subsection{Returns}{
A new \code{Logicals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Logicals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set6-deprecated.R
\name{set6-deprecated}
\alias{set6-deprecated}
\title{Deprecated set6 Functions and Classes}
\description{
The functions listed below are deprecated and will be defunct in
the near future. When possible, alternative functions with similar
functionality are also mentioned. Help pages for deprecated functions are
available at \code{help("-deprecated")}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_FuzzySet_FuzzyMultiset.R
\name{FuzzyMultiset}
\alias{FuzzyMultiset}
\title{Mathematical Fuzzy Multiset}
\description{
A general FuzzyMultiset object for mathematical fuzzy multisets, inheriting from \code{FuzzySet}.
}
\details{
Fuzzy multisets generalise standard mathematical multisets to allow for fuzzy relationships. Whereas a
standard, or crisp, multiset assumes that an element is either in a multiset or not, a fuzzy multiset allows
an element to be in a multiset to a particular degree, known as the membership function, which
quantifies the inclusion of an element by a number in [0, 1]. Thus a (crisp) multiset is a
fuzzy multiset where all elements have a membership equal to \eqn{1}. Similarly to \link{Multiset}s, elements
do not need to be unique.
}
\examples{
# Different constructors
FuzzyMultiset$new(1, 0.5, 2, 1, 3, 0)
FuzzyMultiset$new(elements = 1:3, membership = c(0.5, 1, 0))

# Crisp sets are a special case FuzzyMultiset
# Note membership defaults to full membership
FuzzyMultiset$new(elements = 1:5) == Multiset$new(1:5)

f <- FuzzyMultiset$new(1, 0.2, 2, 1, 3, 0)
f$membership()
f$alphaCut(0.3)
f$core()
f$inclusion(0)
f$membership(0)
f$membership(1)

# Elements can be duplicated, and with different memberships,
#  although this is not necessarily sensible.
FuzzyMultiset$new(1, 0.1, 1, 1)

# Like FuzzySets, ordering does not matter.
FuzzyMultiset$new(1, 0.1, 2, 0.2) == FuzzyMultiset$new(2, 0.2, 1, 0.1)
}
\seealso{
Other sets: 
\code{\link{ConditionalSet}},
\code{\link{FuzzySet}},
\code{\link{FuzzyTuple}},
\code{\link{Interval}},
\code{\link{Multiset}},
\code{\link{Set}},
\code{\link{Tuple}}
}
\concept{sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:FuzzySet]{set6::FuzzySet}} -> \code{FuzzyMultiset}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-equals}{\code{FuzzyMultiset$equals()}}
\item \href{#method-isSubset}{\code{FuzzyMultiset$isSubset()}}
\item \href{#method-alphaCut}{\code{FuzzyMultiset$alphaCut()}}
\item \href{#method-clone}{\code{FuzzyMultiset$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="contains">}\href{../../set6/html/Set.html#method-contains}{\code{set6::Set$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="core">}\href{../../set6/html/FuzzySet.html#method-core}{\code{set6::FuzzySet$core()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="inclusion">}\href{../../set6/html/FuzzySet.html#method-inclusion}{\code{set6::FuzzySet$inclusion()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="initialize">}\href{../../set6/html/FuzzySet.html#method-initialize}{\code{set6::FuzzySet$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="membership">}\href{../../set6/html/FuzzySet.html#method-membership}{\code{set6::FuzzySet$membership()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="strprint">}\href{../../set6/html/FuzzySet.html#method-strprint}{\code{set6::FuzzySet$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="FuzzySet" data-id="support">}\href{../../set6/html/FuzzySet.html#method-support}{\code{set6::FuzzySet$support()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-equals"></a>}}
\if{latex}{\out{\hypertarget{method-equals}{}}}
\subsection{Method \code{equals()}}{
Tests if two sets are equal.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzyMultiset$equals(x, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{\link{Set} or vector of \link{Set}s.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Two fuzzy sets are equal if they contain the same elements with the same memberships and
in the same order. Infix operators can be used for:
\tabular{ll}{
Equal \tab \code{==} \cr
Not equal \tab \code{!=} \cr
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are equal to the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isSubset"></a>}}
\if{latex}{\out{\hypertarget{method-isSubset}{}}}
\subsection{Method \code{isSubset()}}{
Test if one set is a (proper) subset of another
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzyMultiset$isSubset(x, proper = FALSE, all = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{proper}}{logical. If \code{TRUE} tests for proper subsets.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
If using the method directly, and not via one of the operators then the additional boolean
argument \code{proper} can be used to specify testing of subsets or proper subsets. A Set is a proper
subset of another if it is fully contained by the other Set (i.e. not equal to) whereas a Set is a
(non-proper) subset if it is fully contained by, or equal to, the other Set.

Infix operators can be used for:
\tabular{ll}{
Subset \tab \code{<} \cr
Proper Subset \tab \code{<=} \cr
Superset \tab \code{>} \cr
Proper Superset \tab \code{>=}
}
}

\subsection{Returns}{
If \code{all} is \code{TRUE} then returns \code{TRUE} if all \code{x} are subsets of the Set, otherwise
\code{FALSE}. If \code{all} is \code{FALSE} then returns a vector of logicals corresponding to each individual
element of \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-alphaCut"></a>}}
\if{latex}{\out{\hypertarget{method-alphaCut}{}}}
\subsection{Method \code{alphaCut()}}{
The alpha-cut of a fuzzy set is defined as the set
\deqn{A_\alpha = \{x \epsilon F | m \ge \alpha\}}{A_\alpha = {x \epsilon F | m \ge \alpha}}
where \eqn{x} is an element in the fuzzy set, \eqn{F}, and \eqn{m} is the corresponding membership.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzyMultiset$alphaCut(alpha, strong = FALSE, create = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{alpha}}{numeric in [0, 1] to determine which elements to return}

\item{\code{strong}}{logical, if \code{FALSE} (default) then includes elements greater than or equal to alpha, otherwise only strictly greater than}

\item{\code{create}}{logical, if \code{FALSE} (default) returns the elements in the alpha cut, otherwise returns a crisp set of the elements}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Elements in \link{FuzzyMultiset} or a \link{Set} of the elements.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FuzzyMultiset$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_Interval_SpecialSet.R
\name{NegReals}
\alias{NegReals}
\title{Set of Negative Real Numbers}
\description{
The mathematical set of negative real numbers,
defined as the union of the set of negative rationals and negative irrationals. i.e.
\deqn{I^- \cup Q^-}{I- U Q-}
where \eqn{I^-}{I-} is the set of negative irrationals and \eqn{Q^-}{Q-} is the set of negative rationals.
}
\seealso{
Other special sets: 
\code{\link{Complex}},
\code{\link{ExtendedReals}},
\code{\link{Integers}},
\code{\link{Logicals}},
\code{\link{Naturals}},
\code{\link{NegIntegers}},
\code{\link{NegRationals}},
\code{\link{PosIntegers}},
\code{\link{PosNaturals}},
\code{\link{PosRationals}},
\code{\link{PosReals}},
\code{\link{Rationals}},
\code{\link{Reals}},
\code{\link{Universal}}
}
\concept{special sets}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:Interval]{set6::Interval}} -> \code{\link[set6:SpecialSet]{set6::SpecialSet}} -> \code{\link[set6:Reals]{set6::Reals}} -> \code{NegReals}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{NegReals$new()}}
\item \href{#method-clone}{\code{NegReals$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="contains">}\href{../../set6/html/Interval.html#method-contains}{\code{set6::Interval$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="equals">}\href{../../set6/html/Interval.html#method-equals}{\code{set6::Interval$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubinterval">}\href{../../set6/html/Interval.html#method-isSubinterval}{\code{set6::Interval$isSubinterval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Interval" data-id="isSubset">}\href{../../set6/html/Interval.html#method-isSubset}{\code{set6::Interval$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SpecialSet" data-id="strprint">}\href{../../set6/html/SpecialSet.html#method-strprint}{\code{set6::SpecialSet$strprint()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{NegReals} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NegReals$new(zero = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{zero}}{logical. If TRUE, zero is included in the set.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{NegReals} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NegReals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SetWrapper_UnionSet.R
\name{UnionSet}
\alias{UnionSet}
\title{Set of Unions}
\description{
UnionSet class for symbolic union of mathematical sets.
}
\details{
The purpose of this class is to provide a symbolic representation for the union of sets that
cannot be represented in a simpler class. Whilst this is not an abstract class, it is not recommended to construct
this class directly but via the set operation methods.
}
\seealso{
Set operations: \link{setunion}, \link{setproduct}, \link{setpower}, \link{setcomplement}, \link{setsymdiff},  \link{powerset}, \link{setintersect}

Other wrappers: 
\code{\link{ComplementSet}},
\code{\link{ExponentSet}},
\code{\link{PowersetSet}},
\code{\link{ProductSet}}
}
\concept{wrappers}
\section{Super classes}{
\code{\link[set6:Set]{set6::Set}} -> \code{\link[set6:SetWrapper]{set6::SetWrapper}} -> \code{UnionSet}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{elements}}{Returns the elements in the object.}

\item{\code{length}}{Returns the number of elements in the object.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{UnionSet$new()}}
\item \href{#method-strprint}{\code{UnionSet$strprint()}}
\item \href{#method-contains}{\code{UnionSet$contains()}}
\item \href{#method-clone}{\code{UnionSet$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="equals">}\href{../../set6/html/SetWrapper.html#method-equals}{\code{set6::SetWrapper$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="SetWrapper" data-id="isSubset">}\href{../../set6/html/SetWrapper.html#method-isSubset}{\code{set6::SetWrapper$isSubset()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{UnionSet} object. It is not recommended to construct this class directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnionSet$new(setlist, lower = NULL, upper = NULL, type = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{setlist}}{\code{list} of \link{Set}s to wrap.}

\item{\code{lower}}{lower bound of new object.}

\item{\code{upper}}{upper bound of new object.}

\item{\code{type}}{closure type of new object.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{UnionSet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-strprint"></a>}}
\if{latex}{\out{\hypertarget{method-strprint}{}}}
\subsection{Method \code{strprint()}}{
Creates a printable representation of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnionSet$strprint(n = 2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{numeric. Number of elements to display on either side of ellipsis when printing.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A character string representing the object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-contains"></a>}}
\if{latex}{\out{\hypertarget{method-contains}{}}}
\subsection{Method \code{contains()}}{
Tests if elements \code{x} are contained in \code{self}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnionSet$contains(x, all = FALSE, bound = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{any. Object or vector of objects to test.}

\item{\code{all}}{logical. If \code{FALSE} tests each \code{x} separately. Otherwise returns \code{TRUE} only if all \code{x} pass test.}

\item{\code{bound}}{logical.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{all == TRUE} then returns \code{TRUE} if all \code{x} are contained in \code{self}, otherwise \code{FALSE}.
If \code{all == FALSE} returns a vector of logicals corresponding to the length of \code{x}, representing
if each is contained in \code{self}. If \code{bound == TRUE} then an element is contained in \code{self} if it
is on or within the (possibly-open) bounds of \code{self}, otherwise \code{TRUE} only if the element is within
\code{self} or the bounds are closed.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnionSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Set_LogicalSet.R
\name{LogicalSet}
\alias{LogicalSet}
\title{Set of Logicals}
\description{
The \code{LogicalSet} is defined as the \link{Set} containing the elements \code{TRUE} and \code{FALSE}.
}
\examples{
l <- LogicalSet$new()
print(l)
l$contains(list(TRUE, 1, FALSE))
}
\section{Super class}{
\code{\link[set6:Set]{set6::Set}} -> \code{LogicalSet}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{LogicalSet$new()}}
\item \href{#method-clone}{\code{LogicalSet$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="add">}\href{../../set6/html/Set.html#method-add}{\code{set6::Set$add()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="contains">}\href{../../set6/html/Set.html#method-contains}{\code{set6::Set$contains()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="equals">}\href{../../set6/html/Set.html#method-equals}{\code{set6::Set$equals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="isSubset">}\href{../../set6/html/Set.html#method-isSubset}{\code{set6::Set$isSubset()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="multiplicity">}\href{../../set6/html/Set.html#method-multiplicity}{\code{set6::Set$multiplicity()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="print">}\href{../../set6/html/Set.html#method-print}{\code{set6::Set$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="remove">}\href{../../set6/html/Set.html#method-remove}{\code{set6::Set$remove()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="strprint">}\href{../../set6/html/Set.html#method-strprint}{\code{set6::Set$strprint()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="set6" data-topic="Set" data-id="summary">}\href{../../set6/html/Set.html#method-summary}{\code{set6::Set$summary()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{LogicalSet} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LogicalSet$new()}\if{html}{\out{</div>}}
}

\subsection{Details}{
The Logical set is the set containing \code{TRUE} and \code{FALSE}.
}

\subsection{Returns}{
A new \code{LogicalSet} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LogicalSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
