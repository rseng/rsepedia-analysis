
<!-- README.md is generated from README.Rmd. Please edit only README.Rmd! -->

# `caracas`: Computer algebra in R

<!-- badges: start -->

[![R build
status](https://github.com/r-cas/caracas/workflows/R-CMD-check/badge.svg)](https://github.com/r-cas/caracas/actions)
[![codecov.io](https://codecov.io/gh/r-cas/caracas/branch/master/graphs/badge.svg)](https://codecov.io/gh/r-cas/caracas?branch=master)
<!-- badges: end -->

## Installation

`caracas` is available on CRAN and can be installed as usual:

    install.packages('caracas')

Please ensure that you have SymPy installed, or else install it:

``` r
if (!caracas::has_sympy()) {
  caracas::install_sympy() 
}
```

To build and install from Github with vignettes run this command from
within `R` (please install `remotes` first if not already installed):

    # install.packages('remotes')
    remotes::install_github("r-cas/caracas", build_vignettes = TRUE)

You can also install the package without vignettes if needed as follows:

    remotes::install_github("r-cas/caracas")

## Configuring the Python environment

The `caracas` package uses the
[`reticulate`](https://github.com/rstudio/reticulate) package (to run
Python code). Thus, if you wish to configure your Python environment,
you need to 1) load `reticulate`, 2) configure the Python environment,
and 3) load `caracas`. The Python environment can be configured as
described
[here](https://rstudio.github.io/reticulate/articles/versions.html).
Again, this need to be done *before* loading `caracas`.

## Development site

See <https://github.com/r-cas/caracas>.

## Online documentation

See <https://r-cas.github.io/caracas/>.

## Origin of name

The name “caracas” is intended to mean “(inter)face to computer algebra
system(s)” - notice that “cara” is Spanish (Castellano to be precise)
for “face”.

## Code of conduct

Please note that the `caracas` project is released with a Contributor
Code of Conduct (available in `CODE_OF_CONDUCT.md`). By contributing to
this project, you agree to abide by its terms.

## Brief introduction

``` r
library(caracas)
```

``` r
x <- symbol('x')
eq <- 2*x^2 - x
eq
#> [caracas]:    2    
#>            2⋅x  - x
as.character(eq)
#> [1] "2*x^2 - x"
as_expr(eq)
#> expression(2 * x^2 - x)
tex(eq)
#> [1] "2 x^{2} - x"
```

``` r
solve_sys(eq, x)
#> Solution 1:
#>   x =  0 
#> Solution 2:
#>   x =  1/2
der(eq, x)
#> [caracas]: 4⋅x - 1
subs(eq, x, "y")
#> [caracas]:    2    
#>            2⋅y  - y
```

``` r
A <- matrix(c("x", 2, 0, "2*x"), 2, 2)
B <- as_sym(A)
B
#> [caracas]: ⎡x   0 ⎤
#>            ⎢      ⎥
#>            ⎣2  2⋅x⎦
Binv <- inv(B)
Binv
#> [caracas]: ⎡ 1      ⎤
#>            ⎢ ─    0 ⎥
#>            ⎢ x      ⎥
#>            ⎢        ⎥
#>            ⎢-1    1 ⎥
#>            ⎢───  ───⎥
#>            ⎢  2  2⋅x⎥
#>            ⎣ x      ⎦
tex(Binv)
#> [1] "\\left[\\begin{matrix}\\frac{1}{x} & 0\\\\- \\frac{1}{x^{2}} & \\frac{1}{2 x}\\end{matrix}\\right]"
```

``` r
eigenval(Binv)
#> [[1]]
#> [[1]]$eigval
#> [caracas]: 1
#>            ─
#>            x
#> 
#> [[1]]$eigmult
#> [1] 1
#> 
#> 
#> [[2]]
#> [[2]]$eigval
#> [caracas]:  1 
#>            ───
#>            2⋅x
#> 
#> [[2]]$eigmult
#> [1] 1
```

Please find more examples in the other vignettes available at
<https://r-cas.github.io/caracas/>.

## Contribute, issues, and support

Please use the issue tracker at
<https://github.com/r-cas/caracas/issues> if you want to notify us of an
issue or need support. If you want to contribute, please either create
an issue or make a pull request.
# caracas 1.1.2

* `sympy_func(x, fun)` first tries calling `fun` on `x`; and if it does not exist it tries from the global namespace
* New function: `mat_pow()` for raising a matrix to a power (not component-wise)
* New function: `expand_func()` added
* Bug with `Ops` (functions) fixed


# caracas 1.1.1

* Journal of Open Source Software submission

# caracas 1.1.0

* Global symbol assignment by `def_sym()` (#18)
* Linear algebra: New `do_la()` function with convinience functions like `eigenval()`, `eigenvec()`, `QRdecomposition()`; new vignette demonstrating these
* Assumptions being made available, see e.g. `symbol()` and `ask()`
* Arbitrary precision arithmetic: `N()` function and vignette on 
  "Arbitrary precision arithmetic"
* Rename `eigen_val()`/`eigen_vec()` to `eigenval()`/`eigenvec()`
* More clear naming convention: R has expressions and caracas has symbols; 
  in this connection `as_r()` was renamed to `as_expr()` and 
  `as_symbol()` to `as_sym()`. Also, `as_sym()` changed argument from `declare_variables` to `declare_symbols`.
* Changed internals such that `der()`, `der2()` and `solve_sys()` now takes multiple variables with `list()` (or as a vector symbol) instead of `c()`; see also `matrify()` and `listify()`
* Added `diag_()` and `matrix_()` (postfix `_` to avoid name clashes)
* `sumf()` renamed to `sum_()` and `prodf()` to `prod_()` (postfix `_` to avoid name clashes)
* `intf()` renamed to `int()` and `limf()` to `lim()` (because there are no name clashes with base R)
* Call SymPy functions directy with `sympy_func()`
* Added `taylor()` and `drop_remainder()`
* Minor bugs fixed

# caracas 1.0.1

* Require Python 3

# caracas 1.0.0

* An entire new interface for using SymPy, including symbols, symbolic 
  matrices, solving equations, limits and lots of other functionality.

# caracas 0.0.1

* Initial release
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our community a harassment-free experience for everyone, regardless of age, body size, visible or invisible disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming, diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes, and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of acceptable behavior and will take appropriate and fair corrective action in response to any behavior that they deem inappropriate, threatening, offensive, or harmful.

Community leaders have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, and will communicate reasons for moderation decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when an individual is officially representing the community in public spaces. Examples of representing our community include using an official e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported to the community leaders responsible for enforcement at the corresponding package maintainer. All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing clarity around the nature of the violation and an explanation of why the behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of actions.

**Consequence**: A warning with consequences for continued behavior. No interaction with the people involved, including unsolicited interaction with those enforcing the Code of Conduct, for a specified period of time. This includes avoiding interactions in community spaces as well as external channels like social media. Violating these terms may lead to a temporary or permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public communication with the community for a specified period of time. No public or private interaction with the people involved, including unsolicited interaction with those enforcing the Code of Conduct, is allowed during this period. Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community standards, including sustained inappropriate behavior,  harassment of an individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the project community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 2.0,
available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://www.contributor-covenant.org/translations.
---
title: 'caracas: Computer algebra in R'
authors:
- affiliation: 1
  name: Mikkel Meyer Andersen
  orcid: 0000-0002-0234-0266
- affiliation: 1
  name: Søren Højsgaard
  orcid: 0000-0002-3269-9552
date: "22 June 2021"
bibliography: paper.bib
tags:
- cas
- mathematics
- symbolic mathematics
- statistics
- tex
- latex
affiliations:
- index: 1
  name: Department of Mathematical Sciences, Aalborg University, Denmark
---

# Summary

`caracas` is an `R` [@R] package that enables a 
computer algebra system (CAS) within `R` via the open source Python 
CAS `SymPy` [@sympy], which is made possible via `reticulate` [@reticulate]. 
`caracas` is published at The Comprehensive R Archive Network (CRAN) [@R] at <https://cran.r-project.org/package=caracas>, its source is available at <https://github.com/r-cas/caracas> and the documentation is available at <https://r-cas.github.io/caracas/>.

Much work went into integrating `caracas` into `R` such that `caracas` behaves much like 
other `R` libraries and objects. 

`caracas` contains a number of vignettes demonstrating both basic functionality like solving equations 
as well as more advanced tasks like finding the concentration and covariance matrix in a dynamic linear model. 

Compared to other CAS `R` packages like `Ryacas` [@Andersen2019] based on `yacas` [@yacas;@Pinkus2002], 
`caracas` is more feature complete, for example with respect to solving equations.

# Statement of Need

From a statistician's perspective, `R` is excellent for data handling,
graphics, for model fitting and statistical inference and as a
programming environment. However, `R` largely lacks the ability to
perform symbolic computations. That is, `R` only supports to a small
extent the step from posing a problem (for example a model) in
mathematical terms over symbolic manipulations of the model and further onto a stage where a model can be combined with data. The `caracas` provides capabilities for these steps directly in `R`. Topics that can be handled in `caracas` include:

* Sums, 
* limits, 
* integration, 
* differentiation, 
* symbolic matrices and vectors,
* simplification of mathematical expressions and
* outputting in TeX format.

Several (commercial) systems are available for such tasks (and many more), e.g. Mathematica [@Mathematica] and Maple [@Maple]. However, we will argue that there is a virtue in being able to handle such tasks directly from within `R` using the familiar `R` syntax. Moreover, it is an integrated part of the design of `caracas` that it is straightforward to coerce a mathematical object into an `R` expression which can, e.g., be evaluated numerically. 


# Acknowledgements

We would like to thank the R Consortium for financial support for
creating the `caracas` package ([link to details on the funded project](https://www.r-consortium.org/projects/awarded-projects/2019-group-2#Symbolic+mathematics+in+R+with+SymPy)) and to users for pinpointing points
that can be improved in `caracas`.

# References
# Issues with derivatives, matrices and lists/vectors


```r
load_all()
```

```
## Loading caracas
```

```
## 
## Attaching package: 'testthat'
```

```
## The following object is masked from 'package:devtools':
## 
##     test_file
```

## Helper functions


```r
do_unbracket <- function(x){
  gsub("\\[(.+?)\\]", "\\1", x)
}

do_bracket <- function(x){
  paste0("[", x, "]")
}

do_split_rows <- function(x){
  xx3 <- strsplit(x, "\\],")[[1]]
  for (i in 1:(length(xx3)-1))
    xx3[i] <- paste(xx3[i], "]")
  xx3
  xx3 <- lapply(xx3, function(o) gsub("[[:space:]]*", "", o))
  xx3
  xx4 <- gsub("\\[+(.+?)\\]+", "\\1", xx3)
  ## split by ,
  out <- strsplit(xx4, ",")
}

do_comma <- function(x){
  if (is.list(x))
    lapply(x, paste, collapse=", ")
  else
    paste0(x, collapse=", ")
}
```

## Michaelis menten


```r
N <- 3
y <- as_sym(matrix(paste0("y", 1:N)))
x <- as_sym(matrix(paste0("x", 1:N)))
b <- as_sym(paste0("b", 1:2))

num <- b[1] * x
den <- b[2] + x
```

## Different representations



```r
## mu a matrix; a convention?
mu1 <- num / den
mu1
```

```
## [caracas]: ⎡ b₁⋅x₁    b₁⋅x₂    b₁⋅x₃ ⎤
##            ⎢───────  ───────  ───────⎥
##            ⎣b₂ + x₁  b₂ + x₂  b₂ + x₃⎦ᵀ
```

```r
symbol_is_matrix(mu1)
```

```
## [1] TRUE
```

```r
mu2 <-mu1 %>% remove_mat_prefix %>% do_unbracket %>% as_sym
mu2
```

```
## [caracas]: ⎡ b₁⋅x₁    b₁⋅x₂    b₁⋅x₃ ⎤
##            ⎢───────, ───────, ───────⎥
##            ⎣b₂ + x₁  b₂ + x₂  b₂ + x₃⎦
```

```r
symbol_is_matrix(mu2)
```

```
## [1] FALSE
```

```r
mu3 <- as_sym(do_unbracket(mu2))
mu3
```

```
## [caracas]: ⎛ b₁⋅x₁    b₁⋅x₂    b₁⋅x₃ ⎞
##            ⎜───────, ───────, ───────⎟
##            ⎝b₂ + x₁  b₂ + x₂  b₂ + x₃⎠
```

```r
symbol_is_matrix(mu3)
```

```
## [1] FALSE
```

```r
mu1$pyobj
```

```
## Matrix([[b1*x1/(b2 + x1)], [b1*x2/(b2 + x2)], [b1*x3/(b2 + x3)]])
```

```r
mu2$pyobj
```

```
## [b1*x1/(b2 + x1), b1*x2/(b2 + x2), b1*x3/(b2 + x3)]
```

```r
mu3$pyobj
```

```
## (b1*x1/(b2 + x1), b1*x2/(b2 + x2), b1*x3/(b2 + x3))
```

### gradients


```r
g1 <- der(mu1, b)
g2 <- der(mu2, b)
g3 <- der(mu3, b)
g1
```

```
## [caracas]: ⎡           ⎡ -b₁⋅x₁   ⎤⎤
##            ⎢           ⎢──────────⎥⎥
##            ⎢⎡   x₁  ⎤  ⎢         2⎥⎥
##            ⎢⎢───────⎥  ⎢(b₂ + x₁) ⎥⎥
##            ⎢⎢b₂ + x₁⎥  ⎢          ⎥⎥
##            ⎢⎢       ⎥  ⎢ -b₁⋅x₂   ⎥⎥
##            ⎢⎢   x₂  ⎥  ⎢──────────⎥⎥
##            ⎢⎢───────⎥  ⎢         2⎥⎥
##            ⎢⎢b₂ + x₂⎥  ⎢(b₂ + x₂) ⎥⎥
##            ⎢⎢       ⎥  ⎢          ⎥⎥
##            ⎢⎢   x₃  ⎥  ⎢ -b₁⋅x₃   ⎥⎥
##            ⎢⎢───────⎥  ⎢──────────⎥⎥
##            ⎢⎣b₂ + x₃⎦  ⎢         2⎥⎥
##            ⎣           ⎣(b₂ + x₃) ⎦⎦
```

```r
g2
```

```
## [caracas]: ⎡    x₁          x₂          x₃    ⎤
##            ⎢ ───────     ───────     ───────  ⎥
##            ⎢ b₂ + x₁     b₂ + x₂     b₂ + x₃  ⎥
##            ⎢                                  ⎥
##            ⎢ -b₁⋅x₁      -b₁⋅x₂      -b₁⋅x₃   ⎥
##            ⎢──────────  ──────────  ──────────⎥
##            ⎢         2           2           2⎥
##            ⎣(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦
```

```r
g3
```

```
## [caracas]: ⎡    x₁          x₂          x₃    ⎤
##            ⎢ ───────     ───────     ───────  ⎥
##            ⎢ b₂ + x₁     b₂ + x₂     b₂ + x₃  ⎥
##            ⎢                                  ⎥
##            ⎢ -b₁⋅x₁      -b₁⋅x₂      -b₁⋅x₃   ⎥
##            ⎢──────────  ──────────  ──────────⎥
##            ⎢         2           2           2⎥
##            ⎣(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦
```

```r
symbol_is_matrix(g1)
```

```
## [1] FALSE
```

```r
symbol_is_matrix(g2)
```

```
## [1] FALSE
```

```r
symbol_is_matrix(g3)
```

```
## [1] FALSE
```

```r
g1a <- g1$pyobj %>% do_unbracket %>% as_sym
g2a <- g2$pyobj %>% do_unbracket %>% as_sym # Dimension lost
g3a <- g3$pyobj %>% do_unbracket %>% as_sym # Dimension lost
g1a
```

```
## [caracas]: ⎡    x₁          x₂          x₃    ⎤
##            ⎢ ───────     ───────     ───────  ⎥
##            ⎢ b₂ + x₁     b₂ + x₂     b₂ + x₃  ⎥
##            ⎢                                  ⎥
##            ⎢ -b₁⋅x₁      -b₁⋅x₂      -b₁⋅x₃   ⎥
##            ⎢──────────  ──────────  ──────────⎥
##            ⎢         2           2           2⎥
##            ⎣(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦
```

```r
g2a
```

```
## [caracas]: ⎡   x₁       x₂       x₃     -b₁⋅x₁      -b₁⋅x₂      -b₁⋅x₃   ⎤
##            ⎢───────, ───────, ───────, ──────────, ──────────, ──────────⎥
##            ⎢b₂ + x₁  b₂ + x₂  b₂ + x₃           2           2           2⎥
##            ⎣                           (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦
```

```r
g3a
```

```
## [caracas]: ⎡   x₁       x₂       x₃     -b₁⋅x₁      -b₁⋅x₂      -b₁⋅x₃   ⎤
##            ⎢───────, ───────, ───────, ──────────, ──────────, ──────────⎥
##            ⎢b₂ + x₁  b₂ + x₂  b₂ + x₃           2           2           2⎥
##            ⎣                           (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦
```

```r
symbol_is_matrix(g1a)
```

```
## [1] TRUE
```

```r
symbol_is_matrix(g2a)
```

```
## [1] FALSE
```

```r
symbol_is_matrix(g3a)
```

```
## [1] FALSE
```

## Jacobi matrix


```r
h1 <- der2(mu1, b)
h2 <- der2(mu2, b)
h3 <- der2(mu3, b)
h1
```

```
## [caracas]: ⎡              ⎡   -x₁    ⎤⎤
##            ⎢              ⎢──────────⎥⎥
##            ⎢              ⎢         2⎥⎥
##            ⎢              ⎢(b₂ + x₁) ⎥⎥
##            ⎢              ⎢          ⎥⎥
##            ⎢    ⎡0⎤       ⎢   -x₂    ⎥⎥
##            ⎢    ⎢ ⎥       ⎢──────────⎥⎥
##            ⎢    ⎢0⎥       ⎢         2⎥⎥
##            ⎢    ⎢ ⎥       ⎢(b₂ + x₂) ⎥⎥
##            ⎢    ⎣0⎦       ⎢          ⎥⎥
##            ⎢              ⎢   -x₃    ⎥⎥
##            ⎢              ⎢──────────⎥⎥
##            ⎢              ⎢         2⎥⎥
##            ⎢              ⎣(b₂ + x₃) ⎦⎥
##            ⎢                          ⎥
##            ⎢⎡   -x₁    ⎤  ⎡ 2⋅b₁⋅x₁  ⎤⎥
##            ⎢⎢──────────⎥  ⎢──────────⎥⎥
##            ⎢⎢         2⎥  ⎢         3⎥⎥
##            ⎢⎢(b₂ + x₁) ⎥  ⎢(b₂ + x₁) ⎥⎥
##            ⎢⎢          ⎥  ⎢          ⎥⎥
##            ⎢⎢   -x₂    ⎥  ⎢ 2⋅b₁⋅x₂  ⎥⎥
##            ⎢⎢──────────⎥  ⎢──────────⎥⎥
##            ⎢⎢         2⎥  ⎢         3⎥⎥
##            ⎢⎢(b₂ + x₂) ⎥  ⎢(b₂ + x₂) ⎥⎥
##            ⎢⎢          ⎥  ⎢          ⎥⎥
##            ⎢⎢   -x₃    ⎥  ⎢ 2⋅b₁⋅x₃  ⎥⎥
##            ⎢⎢──────────⎥  ⎢──────────⎥⎥
##            ⎢⎢         2⎥  ⎢         3⎥⎥
##            ⎣⎣(b₂ + x₃) ⎦  ⎣(b₂ + x₃) ⎦⎦
```

```r
h2
```

```
## [caracas]: ⎡                                      ⎡   -x₁         -x₂         -x₃    ⎤⎤
##            ⎢⎡    0           0           0     ⎤  ⎢──────────  ──────────  ──────────⎥⎥
##            ⎢⎢                                  ⎥  ⎢         2           2           2⎥⎥
##            ⎢⎢   -x₁         -x₂         -x₃    ⎥  ⎢(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎥⎥
##            ⎢⎢──────────  ──────────  ──────────⎥  ⎢                                  ⎥⎥
##            ⎢⎢         2           2           2⎥  ⎢ 2⋅b₁⋅x₁     2⋅b₁⋅x₂     2⋅b₁⋅x₃  ⎥⎥
##            ⎢⎣(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦  ⎢──────────  ──────────  ──────────⎥⎥
##            ⎢                                      ⎢         3           3           3⎥⎥
##            ⎣                                      ⎣(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦⎦
```

```r
h3
```

```
## [caracas]: ⎡                                      ⎡   -x₁         -x₂         -x₃    ⎤⎤
##            ⎢⎡    0           0           0     ⎤  ⎢──────────  ──────────  ──────────⎥⎥
##            ⎢⎢                                  ⎥  ⎢         2           2           2⎥⎥
##            ⎢⎢   -x₁         -x₂         -x₃    ⎥  ⎢(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎥⎥
##            ⎢⎢──────────  ──────────  ──────────⎥  ⎢                                  ⎥⎥
##            ⎢⎢         2           2           2⎥  ⎢ 2⋅b₁⋅x₁     2⋅b₁⋅x₂     2⋅b₁⋅x₃  ⎥⎥
##            ⎢⎣(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦  ⎢──────────  ──────────  ──────────⎥⎥
##            ⎢                                      ⎢         3           3           3⎥⎥
##            ⎣                                      ⎣(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦⎦
```

```r
symbol_is_matrix(h1)
```

```
## [1] FALSE
```

```r
symbol_is_matrix(h2)
```

```
## [1] FALSE
```

```r
symbol_is_matrix(33)
```

```
## [1] FALSE
```

```r
h1a <- h1$pyobj %>% do_unbracket %>% do_unbracket %>% as_sym
h2a <- h2$pyobj %>% do_unbracket %>% do_unbracket %>% as_sym # Dimension lost
h3a <- h3$pyobj %>% do_unbracket %>% do_unbracket %>% as_sym # Dimension lost
h1a
```

```
## [caracas]: ⎡                                       -x₁         -x₂         -x₃    ⎤
##            ⎢    0           0           0       ──────────  ──────────  ──────────⎥
##            ⎢                                             2           2           2⎥
##            ⎢                                    (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎥
##            ⎢                                                                      ⎥
##            ⎢   -x₁         -x₂         -x₃       2⋅b₁⋅x₁     2⋅b₁⋅x₂     2⋅b₁⋅x₃  ⎥
##            ⎢──────────  ──────────  ──────────  ──────────  ──────────  ──────────⎥
##            ⎢         2           2           2           3           3           3⎥
##            ⎣(b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃)   (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦
```

```r
h2a
```

```
## [caracas]: ⎡            -x₁         -x₂         -x₃         -x₁         -x₂         -x₃  
##            ⎢0, 0, 0, ──────────, ──────────, ──────────, ──────────, ──────────, ────────
##            ⎢                  2           2           2           2           2          
##            ⎣         (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃)   (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃
##            
##                 2⋅b₁⋅x₁     2⋅b₁⋅x₂     2⋅b₁⋅x₃  ⎤
##            ──, ──────────, ──────────, ──────────⎥
##             2           3           3           3⎥
##            )   (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦
```

```r
h3a
```

```
## [caracas]: ⎡            -x₁         -x₂         -x₃         -x₁         -x₂         -x₃  
##            ⎢0, 0, 0, ──────────, ──────────, ──────────, ──────────, ──────────, ────────
##            ⎢                  2           2           2           2           2          
##            ⎣         (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃)   (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃
##            
##                 2⋅b₁⋅x₁     2⋅b₁⋅x₂     2⋅b₁⋅x₃  ⎤
##            ──, ──────────, ──────────, ──────────⎥
##             2           3           3           3⎥
##            )   (b₂ + x₁)   (b₂ + x₂)   (b₂ + x₃) ⎦
```

```r
symbol_is_matrix(h1a)
```

```
## [1] TRUE
```

```r
symbol_is_matrix(h2a)
```

```
## [1] FALSE
```

```r
symbol_is_matrix(h3a)
```

```
## [1] FALSE
```

