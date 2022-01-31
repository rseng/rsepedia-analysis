
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

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit only README.Rmd! -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)
```

# `caracas`: Computer algebra in R

<!-- badges: start -->
[![R build status](https://github.com/r-cas/caracas/workflows/R-CMD-check/badge.svg)](https://github.com/r-cas/caracas/actions) 
[![codecov.io](https://codecov.io/gh/r-cas/caracas/branch/master/graphs/badge.svg)](https://codecov.io/gh/r-cas/caracas?branch=master)
<!-- badges: end -->
  
## Installation

`caracas` is available on CRAN and can be installed as usual:

```
install.packages('caracas')
```

Please ensure that you have SymPy installed, or else install it:

```{r, eval = FALSE}
if (!caracas::has_sympy()) {
  caracas::install_sympy() 
}
```

To build and install from Github with vignettes run this command from within `R` (please install `remotes` first if not already installed):

```
# install.packages('remotes')
remotes::install_github("r-cas/caracas", build_vignettes = TRUE)
```

You can also install the package without vignettes if needed as follows:

```
remotes::install_github("r-cas/caracas")
```

## Configuring the Python environment

The `caracas` package uses the [`reticulate`](https://github.com/rstudio/reticulate) package (to run Python code). Thus, if you wish to configure your Python environment, you need to 1) load `reticulate`, 2) configure the Python environment, and 3) load `caracas`. The Python environment can be configured as described [here](https://rstudio.github.io/reticulate/articles/versions.html). Again, this need to be done *before* loading `caracas`.

## Development site

See <https://github.com/r-cas/caracas>.

## Online documentation

See <https://r-cas.github.io/caracas/>.

## Origin of name

The name "caracas" is intended to mean "(inter)face to computer algebra system(s)" - notice that "cara" is Spanish (Castellano to be precise) for "face".

## Code of conduct

Please note that the `caracas` project is released with a Contributor Code of Conduct (available in `CODE_OF_CONDUCT.md`). By contributing to this project, you agree to abide by its terms.

## Brief introduction

```{r, message=FALSE}
library(caracas)
```

```{r, include = FALSE}
if (!has_sympy()) {
  # SymPy not available, so the chunks shall not be evaluated
  knitr::opts_chunk$set(eval = FALSE)
}
```

```{r}
x <- symbol('x')
eq <- 2*x^2 - x
eq
as.character(eq)
as_expr(eq)
tex(eq)
```

```{r}
solve_sys(eq, x)
der(eq, x)
subs(eq, x, "y")
```

```{r}
A <- matrix(c("x", 2, 0, "2*x"), 2, 2)
B <- as_sym(A)
B
Binv <- inv(B)
Binv
tex(Binv)
```

```{r}
eigenval(Binv)
```

Please find more examples in the other vignettes available at <https://r-cas.github.io/caracas/>.

## Contribute, issues, and support

Please use the issue tracker at <https://github.com/r-cas/caracas/issues> 
if you want to notify us of an issue or need support.
If you want to contribute, please either create an issue or make a pull request.
---
title: "Reference"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Reference}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, message=FALSE}
library(caracas)
```

```{r, include = FALSE}
inline_code <- function(x) {
  x
}

if (!has_sympy()) {
  # SymPy not available, so the chunks shall not be evaluated
  knitr::opts_chunk$set(eval = FALSE)
  
  inline_code <- function(x) {
    deparse(substitute(x))
  }
}
```


## Quick start

```{r}
x <- symbol('x')
as.character(x)
x
as_expr(x)
```


```{r}
2*x
y <- symbol('y')
sqrt(3*x^y)
```

```{r}
z <- cos(x)^2 + sin(x)^2
z
simplify(z)
tex(z)
```

```{r}
z <- cos(x)*cos(y) - sin(x)*sin(y)
z
simplify(z)
z <- cos(x + y)
z
expand(z)
expand_trig(z)
```

```{r}
x <- symbol('x')
y <- symbol('y')
z <- log(x*y)
z
expand_log(z)
```


### Sums

```{r}
x <- symbol("x")
sum_(1/x, "x", 1, 10)
sum_(1/x, x, 1, 10)
s <- sum_(1/x, "x", 1, 10)
as_expr(s)
sum(1/(1:10))
n <- symbol("n")
simplify(sum_(x, x, 1, n))
```


### Products

```{r}
x <- symbol("x")
p <- prod_(1/x, "x", 1, 10)
p
as_expr(p)
prod(1/(1:10))
n <- symbol("n")
prod_(x, x, 1, n)
```


### Integrals

```{r}
x <- symbol("x")

int(1/x, x, 1, 10)
i1 <- int(1/x, x, 1, 10, doit = FALSE)
i1
tex(i1)
doit(i1)
int(1/x, x)
i1 <- int(1/x, x, doit = FALSE)
i1
tex(i1)
doit(i1)
```

### Limits


```{r}
x <- symbol("x")
lim(sin(x)/x, "x", 0)
lim(1/x, "x", 0, dir = '+')
lim(1/x, "x", 0, dir = '-')
```

We can also postpone evaluation:

```{r}
x <- symbol("x")
lim(sin(x)/x, "x", 0)
lim(sin(x)/x, x, 0)
```

```{r}
res <- lim(sin(x)/x, "x", 0, doit = FALSE)
res
as.character(res)
tex(res)
doit(res)
as_expr(res)
```

### Derivatives

Note that the function is called `d()` and not `deriv()`.

```{r}
x <- symbol("x")
y <- symbol("y")
f <- 3*x^2 + x*y^2
f
as_expr(f)
der(f, "x")
der(f, x)
der(f, c("x", "y"))
der(f, list(x, y))
f1 <- der(f, list(x, y))
f1
as.character(f1)
as_expr(f1)
eval(as_expr(f1), list(x = 1, y = 2))
der(f1, list(x, y))
f2 <- der2(f, list(x, y))
f2
as_expr(f2)
eval(as_expr(f2), list(x = 1, y = 2))
```


```{r}
x <- symbol("x")
y <- symbol("y")
f <- eval_to_symbol("[3*x**2 + x*y**2, 2*x, 5*y]")
f
der(f, list(x, y))
```

### Taylor expansion

```{r}
def_sym(x)
f <- cos(x)
ft_with_O <- taylor(f, x0 = 0, n = 4+1)
ft_with_O
ft_with_O %>% drop_remainder() %>% as_expr()
```



## Linear algebra

```{r}
A <- matrix(c("x", 0, 0, "2*x"), 2, 2)
A
B <- as_sym(A)
B
2*B
B*B # Component-wise / Hadamard product
dim(B)
sqrt(B)
log(B)
sum(B)
B %*% t(B)
diag(B)
cbind(B, B)
rbind(B, B)
```

```{r}
det(B)
QRdecomposition(B)
```

```{r}
A <- matrix(c("a", 0, 0, 0, "a", "a", "a", 0, 0), 3, 3)
B <- as_sym(A)
eigenval(B)
eigenvec(B)
eigen(eval(as_expr(B), list(a = 2)))
```

```{r}
B
diag(B)
diag(B) <- "b"
B
diag(B)[-2] <- "a"
B
```


## Solve

* Linear system of equations: `inv()` / `solve_lin()`
* Non-linear system of equations: `solve_sys()`

Below find an example with maximising the multinomial likelihood.

```{r}
p <- as_sym(paste0("p", 1:3))
y <- as_sym(paste0("y", 1:3))
a <- as_sym("a")
l <- sum(y*log(p))
l
L <- -l + a*(sum(p) - 1)
L
tex(L)
g <- der(L, list(p, a))
g
sol <- solve_sys(g, list(p, a))
sol
sol[[1L]]$p1
tex(sol[[1L]]$p1)
```

## Substitution

```{r}
x <- symbol('x')
eq <- 2*x^2 - x
eq
subs(eq, x, "y")
```

```{r}
p <- as_sym(paste0("p", 1:3))
y <- as_sym(paste0("y", 1:3))
a <- as_sym("a")
l <- sum(y*log(p))
L <- -l + a*(sum(p) - 1)
g <- der(L, c(a, p))
sols <- solve_sys(g, list(a, p))
sol <- sols[[1L]]
sol
H <- der2(L, list(p, a))
H
H_sol <- subs_lst(H, sol)
H_sol
```


## Subsetting

Note that all vectors in `caracas` are column vectors.

```{r}
A <- matrix(c("a", 0, 0, 0, "a", "a", "a", 0, 0), 3, 3)
B <- as_sym(A)
B[, 2]
B[, -2]
B[1, ]
B[1, , drop = FALSE] # Note this is a 1x3 matrix
B[, 2] <- "x"
B
```


## Using `SymPy` directly

```{r}
sympy <- get_sympy()
```

```{r}
sympy$diff("2*a*x", "x")
sympy$solve("x**2 - 1", "x")
```

## Assumptions

Below we give a brief example of assumptions.
First consider the Cholesky decomposition of a matrix:

```{r}
A <- matrix(c("x+1", 1, 1, 1), 2, 2) %>% as_sym()
A
```

```{r, error=TRUE}
do_la(A, "cholesky")
```

This fails as `A` is not positive (semi-)definite.

To ensure this, we need to impose restrictions on `x`. 
This is done by defining a symbol with an assumption about positivity:

```{r}
y <- symbol("y", positive = TRUE)
```

We continue and define `B`, where it is important that 
`declare_symbols = FALSE` or else a new `y` will automatically 
be defined by `caracas` overwriting the above definition:

```{r}
B <- as_sym("[[y + 1, 1], [1, 1]]", declare_symbols = FALSE)
B
do_la(B, "cholesky")
```

It is possible to ask for properties (see <https://docs.sympy.org/latest/modules/assumptions/ask.html>):

```{r}
ask(y, "positive")
ask(B, "hermitian")
ask(A, "hermitian")
```


## Output

```{r}
# Multinomial likelihood
p <- as_sym(paste0("p", 1:3))
y <- as_sym(paste0("y", 1:3))
a <- as_sym("a")
l <- sum(y*log(p))
L <- -l + a*(sum(p) - 1)
L
print(L, ascii = TRUE)
g <- der(L, list(p, a))
sol <- solve_sys(g, list(p, a))
sol
print(sol, simplify = FALSE)
```

```{r}
as.character(g)
as_character_matrix(g)
```



### Options

The following options are available:

* `caracas.print.prettyascii`
* `caracas.print.ascii`
* `caracas.print.rowvec`
* `caracas.print.sol.simplify`

```{r}
sol
L
options(caracas.print.prettyascii = TRUE) 
sol
L
options(caracas.print.prettyascii = NULL) # reset to default (FALSE)
```

```{r}
sol
L
options(caracas.print.ascii = TRUE) 
sol
L
options(caracas.print.ascii = NULL) # reset to default (FALSE)
```

```{r}
p
options(caracas.print.rowvec = FALSE)
p
options(caracas.print.rowvec = NULL) # reset to default (TRUE)
```

```{r}
sol
options(caracas.print.sol.simplify = FALSE)
sol
options(caracas.print.sol.simplify = NULL) # reset to default (TRUE)
```



---
title: "Linear algebra in `caracas`"
author: Mikkel Meyer Andersen and Søren Højsgaard
date: "`r date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 3
    number_sections: true
vignette: >
  %\VignetteIndexEntry{Linear algebra in `caracas`}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

```{r, message=FALSE, echo=FALSE}
library(caracas)
##packageVersion("caracas")
```

```{r, include = FALSE}
inline_code <- function(x) {
  x
}

if (!has_sympy()) {
  # SymPy not available, so the chunks shall not be evaluated
  knitr::opts_chunk$set(eval = FALSE)
  
  inline_code <- function(x) {
    deparse(substitute(x))
  }
}
```

This vignette is based on `caracas` version `r
packageVersion("caracas")`. `caracas` is avavailable on CRAN at
[https://cran.r-project.org/package=caracas] and on github at
[https://github.com/r-cas/caracas].

# Elementary matrix operations

## Creating matrices / vectors

We now show different ways to create a symbolic matrix:

```{r}
A <- matrix(c("a", "b", "0", "1"), 2, 2) %>% as_sym()
A
A <- matrix_(c("a", "b", "0", "1"), 2, 2) # note the '_' postfix
A
A <- as_sym("[[a, 0], [b, 1]]")
A

A2 <- matrix_(c("a", "b", "c", "1"), 2, 2)
A2

B <- matrix_(c("a", "b", "0", 
              "c", "c", "a"), 2, 3)
B

b <- matrix_(c("b1", "b2"), nrow = 2)

D <- diag_(c("a", "b")) # note the '_' postfix
D
```

Note that a vector is a matrix in which one of the dimensions is one. 


## Matrix-matrix sum and product

```{r}
A + A2
A %*% B
```

## Hadamard (elementwise) product

```{r}
A * A2
```

## Vector operations

```{r}
x <- as_sym(paste0("x", 1:3))
x
x + x
1 / x
x / x
```


## Reciprocal matrix

```{r}
reciprocal_matrix(A2)
reciprocal_matrix(A2, num = "2*a")
```


## Matrix inverse; solve system of linear equations

Solve $Ax=b$ for $x$:

```{r}
inv(A)
x <- solve_lin(A, b)
x
A %*% x ## Sanity check
```

## Generalized (Penrose-Moore) inverse; solve system of linear equations [TBW]

```{r}
M <- as_sym("[[1, 2, 3], [4, 5, 6]]")
pinv(M)
B <- as_sym("[[7], [8]]") 
B
z <- do_la(M, "pinv_solve", B)
print(z, rowvec = FALSE) # Do not print column vectors as transposed row vectors
```


# More special linear algebra functionality

Below we present convenient functions for performing linear algebra
operations.  If you need a function in SymPy for which we have not
supplied a convinience function (see
<https://docs.sympy.org/latest/modules/matrices/matrices.html>), you
can still call it with the `do_la()` (short for "do linear algebra") function presented at the end of
this section.


## QR decomposition

```{r}
A <- matrix(c("a", "0", "0", "1"), 2, 2) %>% as_sym()
A
qr_res <- QRdecomposition(A)
qr_res$Q
qr_res$R
```

## Eigenvalues and eigenvectors

```{r}
eigenval(A)
```

```{r}
evec <- eigenvec(A)
evec
evec1 <- evec[[1]]$eigvec
evec1
simplify(evec1)

lapply(evec, function(l) simplify(l$eigvec))
```


## Inverse, Penrose-Moore pseudo inverse

```{r}
inv(A)
pinv(cbind(A, A)) # pseudo inverse
```


## Additional functionality for linear algebra

`do_la` short for "do linear algebra"


```{r}
args(do_la)
```

The above functions can be called:

```{r}
do_la(A, "QRdecomposition") # == QRdecomposition(A)
do_la(A, "inv")             # == inv(A)
do_la(A, "eigenvec")        # == eigenvec(A)
do_la(A, "eigenvals")       # == eigenval(A)
```

### Characteristic polynomium

```{r}
cp <- do_la(A, "charpoly")
cp
as_expr(cp)
```

### Rank

```{r}
do_la(A, "rank")
```

### Cofactor

```{r}
A <- matrix(c("a", "b", "0", "1"), 2, 2) %>% as_sym()
A
do_la(A, "cofactor", 0, 1)
do_la(A, "cofactor_matrix")
```

### Echelon form

```{r}
do_la(cbind(A, A), "echelon_form")
```

### Cholesky factorisation

```{r}
B <- as_sym("[[9, 3*I], [-3*I, 5]]")
B
do_la(B, "cholesky")
```

### Gram Schmidt

```{r}
B <- t(as_sym("[[ 2, 3, 5 ], [3, 6, 2], [8, 3, 6]]"))
do_la(B, "GramSchmidt")
```


### Reduced row-echelon form (rref)

```{r}
B <- t(as_sym("[[ 2, 3, 5 ], [4, 6, 10], [8, 3, 6] ]"))
B
B_rref <- do_la(B, "rref")
B_rref
```

### Column space, row space and null space

```{r}
B <- matrix(c(1, 3, 0, -2, -6, 0, 3, 9, 6), nrow = 3) %>% as_sym()
B
columnspace(B)
rowspace(B)
x <- nullspace(B)
x
rref(B)
B %*% x[[1]]
```

### Singular values, svd

```{r}
B <- t(as_sym("[[ 2, 3, 5 ], [4, 6, 10], [8, 3, 6], [8, 3, 6] ]"))
B
singular_values(B)
```


<!-- ```{r} -->
<!-- ##svd(B) FIXME -->
<!-- # # FIXME -->
<!-- # B <- t(as_sym("[[ 2, 3, 5 ], [4, 6, 10], [8, 3, 6], [8, 3, 6] ]")) -->
<!-- # B -->
<!-- #  -->
<!-- # #do_la(B %*% t(B), "left_eigenvects") -->
<!-- #  -->
<!-- # U <- eigenvec(B %*% t(B)) %>% lapply(function(v) v$eigvec) %>% do.call(cbind, .) -->
<!-- # V <- eigenvec(t(B) %*% B) %>% lapply(function(v) v$eigvec) %>% do.call(cbind, .) -->
<!-- #  -->
<!-- # d <- singular_values(B)[seq_len(ncol(U))] -->
<!-- # D <- diag(ncol(U)) %>% as_sym() -->
<!-- # for (i in seq_along(d)) { -->
<!-- #   D[i, i] <- as.character(d[[i]]) -->
<!-- # } -->
<!-- #  -->
<!-- # d <- singular_values(B) %>% do.call(cbind, .) -->
<!-- # D <- as_diag(d) -->
<!-- #  -->
<!-- #  -->
<!-- # r_svd <- as_expr(B) %>% svd() -->
<!-- #  -->
<!-- # r_svd$d -->
<!-- # d %>% as_expr() -->
<!-- #  -->
<!-- # r_svd$u -->
<!-- # U %>% as_expr() -->
<!-- #  -->
<!-- # r_svd$v -->
<!-- # V %>% as_expr() -->
<!-- #  -->
<!-- # r_svd$u %*% diag(r_svd$d) %*% t(r_svd$v) -->
<!-- # U %*% D %*% t(V) %>% as_expr() -->
<!-- ``` -->

---
title: "Introduction to 'caracas'"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to 'caracas'}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, message=FALSE}
library(caracas)
```

```{r, include = FALSE}
inline_code <- function(x) {
  x
}

if (!has_sympy()) {
  # SymPy not available, so the chunks shall not be evaluated
  knitr::opts_chunk$set(eval = FALSE)
  
  inline_code <- function(x) {
    deparse(substitute(x))
  }
}
```

## Quick start

```{r}
x <- symbol('x')
eq <- 2*x^2 - x
eq
as.character(eq)
as_expr(eq)
tex(eq)
```

$$`r inline_code(tex(eq))`$$

```{r}
solve_sys(eq, x)
der(eq, x)
subs(eq, x, "y")
```

## Linear algebra

```{r}
A <- matrix(c("x", 2, 0, "2*x"), 2, 2)
B <- as_sym(A)
B
Binv <- inv(B) # or solve_lin(B)
Binv
tex(Binv)
```

$$`r inline_code(tex(Binv))`$$

```{r}
eigenval(Binv)
eigenvec(Binv)
```


## More examples

Please find more examples in the other vignettes.

---
title: "Using the 'SymPy' object directly"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using the 'SymPy' object directly}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, message=FALSE}
library(caracas)
```

```{r, include = FALSE}
inline_code <- function(x) {
  x
}

if (!has_sympy()) {
  # SymPy not available, so the chunks shall not be evaluated
  knitr::opts_chunk$set(eval = FALSE)
  
  inline_code <- function(x) {
    deparse(substitute(x))
  }
}
```

## Using `SymPy` directly

First we get the `SymPy` object:

```{r}
sympy <- get_sympy()
```

```{r}
sympy$diff("2*a*x", "x")
sympy$solve("x**2 - 1", "x")
```

## Elaborate example

How can we minimise the amount of material used to produce a cylindric
tin can that contains 1 litre.  The cylinder has diameter $d$ and
height $h$. The question is therefore: What is $d$ and $h$?

We introduce the variables `d` (diameter) and `h` (height):

```{r}
d <- sympy$symbols('d')
h <- sympy$symbols('h')
```

The problem is a constrained optimisation problem, and we solve it by
a Lagrange multiplier, and therefore we introduce `lam` (the Lagrange
multiplier):

```{r}
lam <- sympy$symbols('lam')
```

We now set up the problem:

```{r}
area_str <- "Pi/2 * d**2 + Pi * h * d"
vol_str <- "Pi/4 * d**2 * h"
lap_str <- paste0("(", area_str, ") - lam*((", vol_str, ") - 1)")
lap <- sympy$parsing$sympy_parser$parse_expr(
  lap_str,
  local_dict = list('d' = d, 'h' = h, 'lam' = lam))
```

We can now find the gradient:

```{r}
grad <- sympy$derive_by_array(lap, list(d, h, lam))
grad
```

And find the critical points:

```{r}
sol <- sympy$solve(grad, list(d, h, lam), dict = TRUE)
sol
```

We take the one with the real solution:

```{r}
sol[[1]]
```

We now have a short helper function to help getting appropriate `R`
expressions (such a function will be included in later versions of
this package):

```{r}
to_r <- function(x) {
  x <- as.character(x)
  x <- gsub("Pi", "pi", x, fixed = TRUE)
  x <- gsub("**", "^", x, fixed = TRUE)
  x <- parse(text = x)
  return(x)
}

sol_d <- to_r(sol[[1]]$d)
sol_d
eval(sol_d)
sol_h <- to_r(sol[[1]]$h)
sol_h
eval(sol_h)
```

(It is left as an exercise to the reader to show that the critical point indeed is a minimum.)

## Simple example with assumptions

```{r}
x <- sympy$symbols('x')
x$assumptions0
x <- sympy$symbols('x', positive = TRUE)
x$assumptions0
eq <- sympy$parsing$sympy_parser$parse_expr("x**2 - 1",
                                            local_dict = list('x' = x))
sympy$solve(eq, x, dict = TRUE)
```

## Another example with assumptions

```{r}
x <- sympy$symbols('x', positive = TRUE)
eq <- sympy$parsing$sympy_parser$parse_expr("x**3/3 - x",
                                            local_dict = list('x' = x))
eq
grad <- sympy$derive_by_array(eq, x)
grad
sympy$solve(grad, x, dict = TRUE)
```

---
title: "Concentration and covariance matrix in an autoregressive model and in a dynamic linear model"
author: Mikkel Meyer Andersen and Søren Højsgaard
date: "`r date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 3
    number_sections: true
vignette: >
  %\VignetteIndexEntry{Concentration and covariance matrix in an autoregressive model and in a dynamic linear model}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, message=FALSE}
library(caracas)
```

```{r, include = FALSE}
inline_code <- function(x) {
  x
}

if (!has_sympy()) {
  # SymPy not available, so the chunks shall not be evaluated
  knitr::opts_chunk$set(eval = FALSE)
  
  inline_code <- function(x) {
    deparse(substitute(x))
  }
}
```


## Autoregressive model ($AR(1)$)

```{r ar1, echo=F}
N <- 3
L1 <- diag(4)
L1[cbind(1 + (1:N), 1:N)] <- "-a"
L1 <- as_sym(L1)
```

```{r, echo=F}
e <- as_sym(paste0("e", 0:3))
x <- as_sym(paste0("x", 0:3))
u <- as_sym(paste0("u", 1:3))
y <- as_sym(paste0("y", 1:3))
eu <- rbind(e, u)
xy <- rbind(x, y)
```


Consider this model: 
$$
x_i = a x_{i-1} + e_i, \quad i=1, \dots, 3
$$ 
and $x_0=e_0$. All terms $e_0, \dots, e_3$ are independent and $N(0,v^2)$ distributed. 
Let $e=(e_0, \dots, e_3)$ and $x=(x_0, \dots x_3)$. Hence $e \sim N(0, v^2 I)$.
Isolating error terms gives
$$
e= `r inline_code(tex(e))` = `r inline_code(tex(L1))` `r inline_code(tex(x))` = L_1 x
$$

```{r}
<<ar1>>
```

Since $\mathbf{Var}(e)=v^2 I$ we have 
$\mathbf{Var}(e)=v^2 I=L \mathbf{Var}(x) L'$ so the covariance matrix of $x$ is $V_1=\mathbf{Var}(x) = v^2 L^- (L^-)'$ 
 while the concentration matrix (the inverse covariances matrix) is $K=v^{-2}L' L$.

```{r}
def_sym(v2)
L1inv <- inv(L1)
V1 <- v2 * L1inv %*% t(L1inv)
K1 <- (t(L1) %*% L1) / v2
```

```{r, echo=F, results="asis"}
cat(
  "\\begin{align} 
    K_1 &= ", tex(K1), " \\\\ 
   V_1 &= ", tex(V1), " 
  \\end{align}", sep = "")
```


## Dynamic linear model

```{r L2, echo=F}
N <- 3
L2 <- diag("1", 1 + 2*N)
L2[cbind(1 + (1:N), 1:N)] <- "-a"
L2[cbind(1 + N + (1:N), 1 + 1:N)] <- "-b"
L2 <- as_sym(L2)
```

Augment the $AR(1)$ process above with $y_i=b x_i + u_i$ for
$i=1,2,3$.  Suppose $u_i\sim N(0, w^2)$ and all $u_i$ are independent
and inpendent of $e$.
Then
$(e,u)$ can be expressed in terms of $(x,y)$ as
$$
(e,u) = `r inline_code(tex(eu))` = `r inline_code(tex(L2))` `r inline_code(tex(xy))` = L_2 (x,y)
$$
where
```{r}
<<L2>>
```


```{r}
Veu <- diag(1, 7)
diag(Veu)[1:4] <- "v2"
diag(Veu)[5:7] <- "w2"
Veu
Veu <- as_sym(Veu)
Veu
L2inv <- inv(L2) 
V2 <- L2inv %*% Veu %*% t(L2inv) 
K2 <- t(L2) %*% inv(Veu) %*% L2
```

```{r, results="asis", echo=F}
cat(
  "\\begin{align} K_2 &= ", tex(K2), " \\\\ 
                  V_2 &= ", tex(V2), " \\end{align}", sep = "")
```
---
title: "Arbitrary precision arithmetic"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Arbitrary precision arithmetic}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, message=FALSE}
library(caracas)
```

```{r, include = FALSE}
inline_code <- function(x) {
  x
}

if (!has_sympy()) {
  # SymPy not available, so the chunks shall not be evaluated
  knitr::opts_chunk$set(eval = FALSE)
  
  inline_code <- function(x) {
    deparse(substitute(x))
  }
}
```

## An example

```{r}
n_2 <- as_sym("2")
n_pi <- as_sym("pi", declare_symbols = FALSE) # pi is a constant in SymPy
x <- sqrt(n_2) * n_pi
x
N(x)
N(x, 5)
N(x, 50)
as.character(N(x, 50))
```

## Another example

As `SymPy` requires [`mpmath`](https://mpmath.org/) this can be used directly, for example like this (example due to Hans W Borchers, thanks!):

```{r}
mpmath <- reticulate::import('mpmath')
mpmath$mp$dps <- 30
z <- mpmath$mpc(real = '1.0', imag = '63.453')
zeta_z <- mpmath$zeta(z)
zeta_z
as.character(zeta_z$real)
as.character(zeta_z$imag)
```

---
title: "Fastest route through the forest"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Fastest route through the forest}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, message=FALSE}
library(caracas)
```

```{r, include = FALSE}
inline_code <- function(x) {
  x
}

if (!has_sympy()) {
  # SymPy not available, so the chunks shall not be evaluated
  knitr::opts_chunk$set(eval = FALSE)
  
  inline_code <- function(x) {
    deparse(substitute(x))
  }
}
```

A [problem](https://ing.dk/artikel/taenkeboks-hurtigste-vej-mellem-a-b-1150-m-lang-232871) was posted by the Danish newspaper, *Ingeniøren*, and it goes like this:

You are in the middle of a dense forest located at $A$. 
You need to get to $C$ in the fastest way possible, and you can only change direction once. 
You can walk directly via $AB$ to the dedicated walking path $BC$ where you can walk fast, 
you can take the direct path through the forest ($AC$) where you have to walk slower, 
or cross through the forest to the dedicated walking path ($AD$ and then $DC$).

```{r, echo = FALSE, fig.width=8, fig.height=5}
old_mar <- par("mar")
par(mar = c(0, 0, 0, 0))
A <- 300
kat <- sqrt(1000^2 - A^2)

pA <- c(0, A)
pB <- c(0, 0)
pC <- c(kat, 0)

D <- 400
pD <- c(D, 0)

pad <- 40
plot(0, type = 'n', axes = FALSE, ann = FALSE, 
     xlim = c(-pad, 1000+pad), 
     ylim = c(-pad, A+pad))

polygon(c(pA[1]-pad, pB[1]-pad, pC[1]+pad, pC[1]+pad), 
        c(pA[2]+pad, pB[2]-pad, pC[2]-pad, pA[2]+pad), border = NA, col = "green")

points(pA[1], pA[2], type = "p"); text(pA[1], pA[2], labels = "A", pos = 2)
points(pB[1], pB[2], type = "p"); text(pB[1], pB[2], labels = "B", pos = 2)
points(pC[1], pC[2], type = "p"); text(pC[1], pC[2], labels = "C", pos = 4)
points(pD[1], pD[2], type = "p"); text(pD[1], pD[2], labels = "D", pos = 1)

lines(c(pA[1], pB[1]), c(pA[2], pB[2]), type = "l", col = "black", lty = 2)
lines(c(pA[1], pC[1]), c(pA[2], pC[2]), type = "l", col = "black", lty = 2)
lines(c(pB[1], pC[1]), c(pB[2], pC[2]), type = "l", col = "black", lwd = 3)
lines(c(pA[1], pD[1]), c(pA[2], pD[2]), type = "l", col = "black", lty = 2)

legend(1000-200, A, 
       legend = c("5 m/s", "2 m/s"), 
       lty = c(1, 2), 
       lwd = c(3, 1))
text(-pad/3, A/2, labels = "|AB| = 300 m", pos = 3, srt = 90)
text(kat/2, A/1.7, labels = "|AC| = 1,000 m", pos = 3)

par(mar = old_mar)
```

## Information given

* Distances
  * $|AB| =$ 300 m
  * $|AC| =$ 1,000 m
* Velocities
  * $v_{AB} = 2$ m/s
  * $v_{AD} = 2$ m/s
  * $v_{AC} = 2$ m/s
  * $v_{BC} = 5$ m/s

### Length of line segments

We parameterise with $k = |BD|$, the distance between $B$ and $D$. 
That is, how much to walk on fast walking path before crossing into 
the forest.

Formulating using `caracas`:

```{r}
AB <- as_sym('300')
AB
AC <- as_sym('1000')
AC
BC <- sqrt(AC^2 - AB^2)
BC
k <- symbol('k') # |BD|
DC <- BC - k
AD <- sqrt(AB^2 + k^2)
AD
```

So for a distance of $|AD|$, you travel by 5 m/s, and then for a distance of 
$`r inline_code(tex(DC))`$ you travel by 2 m/s. 
Thus it takes $`r inline_code(tex(AD/2))`$ to travel $AD$ and $`r inline_code(tex(DC/5))`$ to travel $DC$.

The question is: What is the fastest way to get from $A$ to $C$?

First, the total duration of the route is:

```{r}
l <- AD/2 + DC/5
l
```

```{r, fig.width=7, fig.height=5}
lfun <- as_expr(l)
lfun
ks <- seq(0, as_expr(AC), length.out = 100)
ls <- eval(lfun, list(k = ks))
plot(ks, ls, type = "l", xlab = "k", ylab = "Time A to C")
```

It looks like a minimum around $k = 150$.

We find the analytical solution by first finding critical points:

```{r}
dl <- der(l, k)
dl
crit_points <- solve_sys(dl, k)
crit_points
best_k <- crit_points[[1]]$k
best_k
```

The type of the critical point is found by considering the Hessian:

```{r}
eval(as_expr(der(dl, k)), list(k = as_expr(best_k)))
```

Thus the critical point is indeed a minimum as suggested by the plot.

The fastest route is thus obtained for 
\[
  k = `r inline_code(tex(best_k))` \approx `r inline_code(round(as_expr(best_k), 2))` .
\]
It has a length of (in meters)
```{r}
DC_best <- BC - best_k
AD_best <- sqrt(AB^2 + best_k^2)
AD_best
best_route <- AD_best + DC_best
best_route
as_expr(best_route)
```
\[
  `r inline_code(tex(best_route))` \approx `r inline_code(round(as_expr(best_route), 2))`
\]
and takes (in seconds)
```{r}
best_l <- subs(l, "k", best_k)
best_l
as_expr(best_l)
```
\[
  `r inline_code(tex(best_l))` \approx `r inline_code(round(as_expr(best_l), 2))`
\]

The best route can be illustrated, too:


```{r, echo = FALSE, fig.width=8, fig.height=5}
old_mar <- par("mar")
par(mar = c(0, 0, 0, 0))
A <- 300
kat <- sqrt(1000^2 - A^2)

pA <- c(0, A)
pB <- c(0, 0)
pC <- c(kat, 0)

D <- as_expr(best_k)
pD <- c(D, 0)

pad <- 40
plot(0, type = 'n', axes = FALSE, ann = FALSE, 
     xlim = c(-pad, 1000+pad), 
     ylim = c(-pad, A+pad))

polygon(c(pA[1]-pad, pB[1]-pad, pC[1]+pad, pC[1]+pad), 
        c(pA[2]+pad, pB[2]-pad, pC[2]-pad, pA[2]+pad), border = NA, col = "green")

points(pA[1], pA[2], type = "p"); text(pA[1], pA[2], labels = "A", pos = 2)
points(pB[1], pB[2], type = "p"); text(pB[1], pB[2], labels = "B", pos = 2)
points(pC[1], pC[2], type = "p"); text(pC[1], pC[2], labels = "C", pos = 4)
points(pD[1], pD[2], type = "p"); text(pD[1], pD[2], labels = "D", pos = 1)

lines(c(pA[1], pB[1]), c(pA[2], pB[2]), type = "l", col = "black", lty = 2)
lines(c(pA[1], pC[1]), c(pA[2], pC[2]), type = "l", col = "black", lty = 2)
lines(c(pB[1], pC[1]), c(pB[2], pC[2]), type = "l", col = "black", lwd = 3)
lines(c(pA[1], pD[1]), c(pA[2], pD[2]), type = "l", col = "black", lty = 2)

lines(c(pA[1], pD[1]), c(pA[2], pD[2]), type = "l", col = "red", lwd = 2.5)
lines(c(pD[1], pC[1]), c(pD[2], pC[2]), type = "l", col = "red", lwd = 2.5)

legend(1000-200, A, 
       legend = c("5 m/s", "2 m/s"), 
       lty = c(1, 2), 
       lwd = c(3, 1))
text(-pad/3, A/2, labels = "|AB| = 300 m", pos = 3, srt = 90)
text(kat/2, A/1.7, labels = "|AC| = 1,000 m", pos = 3)

par(mar = old_mar)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{expand_func}
\alias{expand_func}
\title{Expand a function expression}
\usage{
expand_func(x)
}
\arguments{
\item{x}{A \code{caracas_symbol}}
}
\description{
Expand a function expression
}
\concept{simplify}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_sym.R
\name{as_sym}
\alias{as_sym}
\title{Convert object to symbol}
\usage{
as_sym(x, declare_symbols = TRUE)
}
\arguments{
\item{x}{R object to convert to a symbol}

\item{declare_symbols}{declare detected symbols automatically}
}
\description{
Variables are detected as a
character followed by a number of either:
character, number or underscore.
}
\details{
Default is to declare used variables. Alternatively, the user
must declare them first, e.g. by \code{\link[=symbol]{symbol()}}.

Note that matrices can be defined by specifying a Python matrix,
see below in examples.
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  A <- matrix(c("x", 0, 0, "2*x"), 2, 2)
  A
  B <- as_sym(A)
  B
  2*B
  dim(B)
  sqrt(B)
  D <- as_sym("[[1, 4, 5], [-5, 8, 9]]")
  D
}

}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{sympy_func}
\alias{sympy_func}
\title{Call a SymPy function directly on x}
\usage{
sympy_func(x, fun, ...)
}
\arguments{
\item{x}{Object to call \code{fun} on}

\item{fun}{Function to call}

\item{\dots}{Passed on to \code{fun}}
}
\description{
Call a SymPy function directly on x
}
\examples{
if (has_sympy()) {
  def_sym(x, a)
  p <- (x-a)^4
  p
  q <- p \%>\% sympy_func("expand")
  q
  q \%>\% sympy_func("factor")
  
  def_sym(x, y, z)
  expr <- x*y + x - 3 + 2*x^2 - z*x^2 + x^3
  expr
  expr \%>\% sympy_func("collect", x) 
  
  x <- symbol("x")
  y <- gamma(x+3)
  sympy_func(y, "expand_func")
  expand_func(y)
}
 
}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/numeric_eval.R
\name{N}
\alias{N}
\title{Numerical evaluation}
\usage{
N(x, digits = 15)
}
\arguments{
\item{x}{caracas object}

\item{digits}{Number of digits}
}
\description{
Numerical evaluation
}
\examples{
if (has_sympy()) {
  n_2 <- as_sym("2")
  n_pi <- as_sym("pi", declare_symbols = FALSE)
  x <- sqrt(n_2) * n_pi
  x
  N(x)
  N(x, 5)
  N(x, 50)
  as.character(N(x, 50))
}

}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ops-math.R
\name{matrix-products}
\alias{matrix-products}
\alias{\%*\%}
\alias{\%*\%.caracas_symbol}
\title{Matrix multiplication}
\usage{
x \%*\% y

\method{\%*\%}{caracas_symbol}(x, y)
}
\arguments{
\item{x}{Object \code{x}}

\item{y}{Object \code{y}}
}
\description{
Matrix multiplication

Matrix multiplication
}
\seealso{
\code{\link[base:matmult]{base::\%*\%()}}

\code{\link[base:matmult]{base::\%*\%()}}
}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{as_diag}
\alias{as_diag}
\title{Construct diagonal matrix from vector}
\usage{
as_diag(x)
}
\arguments{
\item{x}{Matrix with 1 row or 1 column that is the
diagonal in a new diagonal matrix}
}
\description{
Construct diagonal matrix from vector
}
\examples{
if (has_sympy()) {
  d <- as_sym(c("a", "b", "c"))
  D <- as_diag(d)
  D
}

}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{eval_to_symbol}
\alias{eval_to_symbol}
\title{Create a symbol from a string}
\usage{
eval_to_symbol(x)
}
\arguments{
\item{x}{String to evaluate}
}
\value{
A \code{caracas_symbol}
}
\description{
Create a symbol from a string
}
\examples{
if (has_sympy()) {
   x <- symbol('x')
   (1+1)*x^2
   lim(sin(x)/x, "x", 0)
}

}
\concept{lowlevel}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{der}
\alias{der}
\title{Symbolic differentiation of an expression}
\usage{
der(expr, vars, simplify = TRUE)
}
\arguments{
\item{expr}{A \code{caracas_symbol}}

\item{vars}{variables to take derivate with respect to}

\item{simplify}{Simplify result}
}
\description{
Symbolic differentiation of an expression
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  y <- symbol("y")
  f <- 3*x^2 + x*y^2
  der(f, x)
  g <- der(f, list(x, y))
  g
  dim(g)
  G <- matrify(g)
  G
  dim(G)
  
  h <- der(g, list(x, y))
  h
  dim(h)
  as.character(h)
  H <- matrify(h)
  H
  dim(H)
  
  g \%>\% 
    der(list(x, y)) \%>\% 
    der(list(x, y)) \%>\% 
    der(list(x, y))
}

}
\concept{calculus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{matrify}
\alias{matrify}
\title{Creates matrix from array symbol}
\usage{
matrify(x)
}
\arguments{
\item{x}{Array symbol to convert to matrix}
}
\description{
Creates matrix from array symbol
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  y <- symbol("y")
  f <- 3*x^2 + x*y^2
  h <- der2(f, list(x, y))
  h
  dim(h)
  H <- matrify(h)
  H
  dim(H)
}

}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{tex}
\alias{tex}
\title{Export object to TeX}
\usage{
tex(x)
}
\arguments{
\item{x}{A \code{caracas_symbol}}
}
\description{
Export object to TeX
}
\concept{output}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{diag.caracas_symbol}
\alias{diag.caracas_symbol}
\title{Matrix diagonal}
\usage{
\method{diag}{caracas_symbol}(x, ...)
}
\arguments{
\item{x}{Object \code{x}}

\item{\dots}{Not used}
}
\description{
Matrix diagonal
}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{expand_log}
\alias{expand_log}
\title{Expand a logarithmic expression}
\usage{
expand_log(x)
}
\arguments{
\item{x}{A \code{caracas_symbol}}
}
\description{
Note that \code{force} as described at
\url{https://docs.sympy.org/latest/tutorial/simplification.html#expand-log} is used
meaning that some assumptions are taken.
}
\examples{
if (has_sympy()) {
  x <- symbol('x')
  y <- symbol('y')
  z <- log(x*y)
  z
  expand_log(z)
}

}
\concept{simplify}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{unbracket}
\alias{unbracket}
\title{Remove inner-most dimension}
\usage{
unbracket(x)
}
\arguments{
\item{x}{Array symbol to collapse dimension from}
}
\description{
Remove inner-most dimension
}
\examples{
if (has_sympy()) {
  x <- as_sym("[[[x1/(b2 + x1)], 
                 [x2/(b2 + x2)], 
                 [x3/(b2 + x3)]], 
                [[-b1*x1/(b2 + x1)^2], 
                 [-b1*x2/(b2 + x2)^2], 
                 [-b1*x3/(b2 + x3)^2]]]")
  x
  unbracket(x)
  
  x <- as_sym("Matrix([[b1*x1/(b2 + x1)], [b1*x2/(b2 + x2)], [b1*x3/(b2 + x3)]])")
  
}

}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{symbol}
\alias{symbol}
\title{Create a symbol}
\usage{
symbol(x, ...)
}
\arguments{
\item{x}{Name to turn into symbol}

\item{\dots}{Assumptions like \code{positive = TRUE}}
}
\value{
A \code{caracas_symbol}
}
\description{
Find available assumptions at
\url{https://docs.sympy.org/latest/modules/core.html#module-sympy.core.assumptions}.
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  2*x
  
  x <- symbol("x", positive = TRUE)
  ask(x, "positive")
}

}
\seealso{
\code{\link[=as_sym]{as_sym()}}
}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg-advanced.R
\name{do_la}
\alias{do_la}
\title{Do linear algebra operation}
\usage{
do_la(x, slot, ...)
}
\arguments{
\item{x}{A matrix for which a property is requested}

\item{slot}{The property requested}

\item{...}{Auxillary arguments}
}
\value{
Returns the requested property of a matrix.
}
\description{
Do linear algebra operation
}
\examples{
if (has_sympy()) {
  A <- matrix(c("a", "0", "0", "1"), 2, 2) \%>\% as_sym()
  
  do_la(A, "QR")
  QRdecomposition(A)
  
  do_la(A, "eigenval")
  eigenval(A)
  
  do_la(A, "eigenvec")
  eigenvec(A)
  
  do_la(A, "inv")
  inv(A)
  
  do_la(A, "echelon_form")
  do_la(A, "rank")
  
  do_la(A, "det") # Determinant
  det(A)
}

}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{subs}
\alias{subs}
\title{Substitute symbol for value}
\usage{
subs(s, x, v)
}
\arguments{
\item{s}{Expression}

\item{x}{Name of symbol (character)}

\item{v}{Value for \code{x}}
}
\description{
Substitute symbol for value
}
\examples{
if (has_sympy()) {
   x <- symbol('x')
   e <- 2*x^2
   e
   subs(e, "x", "2")
   y <- as_sym("2")
   subs(e, "x", y)
}

}
\seealso{
\code{\link[=subs_vec]{subs_vec()}}, \code{\link[=subs_lst]{subs_lst()}}
}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{sum_}
\alias{sum_}
\title{Sum of a function}
\usage{
sum_(f, var, lower, upper, doit = TRUE)
}
\arguments{
\item{f}{Function to take sum of}

\item{var}{Variable to take sum for (either string or \code{caracas_symbol})}

\item{lower}{Lower limit}

\item{upper}{Upper limit}

\item{doit}{Evaluate the sum immediately (or later with \code{\link[=doit]{doit()}})}
}
\description{
Sum of a function
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  s <- sum_(1/x, "x", 1, 10)
  as_expr(s)
  sum(1/(1:10))
  n <- symbol("n")
  simplify(sum_(x, x, 1, n))
}

}
\concept{calculus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{subs_lst}
\alias{subs_lst}
\title{Substitute symbol for of value given by a list}
\usage{
subs_lst(s, x)
}
\arguments{
\item{s}{Expression}

\item{x}{Named list of values}
}
\description{
Useful for substituting solutions into expressions.
}
\examples{
if (has_sympy()) {
     p <- as_sym(paste0("p", 1:3))
     y <- as_sym(paste0("y", 1:3))
     a <- as_sym("a")
     l <- sum(y*log(p))
     L <- -l + a*(sum(p) - 1)
     g <- der(L, c(a, p))
     sols <- solve_sys(g, c(a, p))
     sol <- sols[[1L]]
     sol
     H <- der2(L, c(p, a))
     H
     H_sol <- subs_lst(H, sol)
     H_sol
}

}
\seealso{
\code{\link[=subs]{subs()}}, \code{\link[=subs_vec]{subs_vec()}}
}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{sympy_version}
\alias{sympy_version}
\title{Get 'SymPy' version}
\usage{
sympy_version()
}
\value{
The version of the 'SymPy' available
}
\description{
Get 'SymPy' version
}
\examples{
if (has_sympy()) {
  sympy_version()
}

}
\concept{sympy}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{as_character_matrix}
\alias{as_character_matrix}
\title{Get matrix as character matrix}
\usage{
as_character_matrix(x)
}
\arguments{
\item{x}{caracas symbol}
}
\description{
Get matrix as character matrix
}
\examples{
if (has_sympy()) {
  s  <- as_sym("[[r1, r2, r3], [u1, u2, u3]]")
  s2 <- apply(as_character_matrix(s), 2, function(x) (paste("1/(", x, ")")))
  as_sym(s2)
}

}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve.R
\name{solve_sys}
\alias{solve_sys}
\title{Solves a system of non-linear equations}
\usage{
solve_sys(lhs, rhs, vars)
}
\arguments{
\item{lhs}{Equation (or equations as row vector/1xn matrix)}

\item{rhs}{Equation (or equations as row vector/1xn matrix)}

\item{vars}{vector of variable names or symbols}
}
\value{
A list with solutions (with class \code{caracas_solve_sys_sol}
for compact printing), each element containing a named
list of the variables' values.
}
\description{
If called as \code{solve_sys(lhs, vars)}
the roots are found.
If called as \code{solve_sys(lhs, rhs, vars)}
the solutions to \code{lhs = rhs} for \code{vars} are found.
}
\examples{
if (has_sympy()) {
  x <- symbol('x')
  exp1 <- 2*x + 2
  exp2 <- x
  solve_sys(cbind(exp1), cbind(exp2), x)
  
  x <- symbol("x")
  y <- symbol("y")
  lhs <- cbind(3*x*y - y, x)
  rhs <- cbind(-5*x, y+4)
  sol <- solve_sys(lhs, rhs, list(x, y))
  sol
}

}
\concept{solve}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{drop_remainder}
\alias{drop_remainder}
\title{Remove remainder term}
\usage{
drop_remainder(x)
}
\arguments{
\item{x}{Expression to remove remainder term from}
}
\description{
Remove remainder term
}
\examples{
if (has_sympy()) {
  def_sym(x)
  f <- cos(x)
  ft_with_O <- taylor(f, x0 = 0, n = 4+1)
  ft_with_O
  ft_with_O \%>\% drop_remainder() \%>\% as_expr()
}

}
\seealso{
\code{\link[=taylor]{taylor()}}
}
\concept{calculus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{mat_pow}
\alias{mat_pow}
\title{Matrix power}
\usage{
mat_pow(x, pow = "1")
}
\arguments{
\item{x}{A \code{caracas_symbol}, a matrix.}

\item{pow}{Power to raise matrix \code{x} to}
}
\description{
Matrix power
}
\examples{
if (has_sympy()) {
  M <- matrix_(c("1", "a", "a", 1), 2, 2)
  M
  mat_pow(M, 1/2)
}

}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{t.caracas_symbol}
\alias{t.caracas_symbol}
\title{Transpose of matrix}
\usage{
\method{t}{caracas_symbol}(x)
}
\arguments{
\item{x}{If \code{caracas_symbol} treat as such, else
call \code{\link[base:t]{base::t()}}.}
}
\description{
Transpose of matrix
}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{taylor}
\alias{taylor}
\title{Taylor expansion}
\usage{
taylor(f, x0 = 0, n = 6)
}
\arguments{
\item{f}{Function to be expanded}

\item{x0}{Point to expand around}

\item{n}{Order of remainder term}
}
\description{
Taylor expansion
}
\examples{
if (has_sympy()) {
  def_sym(x)
  f <- cos(x)
  ft_with_O <- taylor(f, x0 = 0, n = 4+1)
  ft_with_O
  ft_with_O \%>\% drop_remainder() \%>\% as_expr()
}

}
\seealso{
\code{\link[=drop_remainder]{drop_remainder()}}
}
\concept{calculus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg-advanced.R
\name{linalg}
\alias{linalg}
\alias{columnspace}
\alias{nullspace}
\alias{rowspace}
\alias{singular_values}
\alias{inv}
\alias{eigenval}
\alias{eigenvec}
\alias{GramSchmidt}
\alias{pinv}
\alias{rref}
\alias{QRdecomposition}
\alias{det}
\title{Do linear algebra operation}
\usage{
columnspace(x)

nullspace(x)

rowspace(x)

singular_values(x)

inv(x)

eigenval(x)

eigenvec(x)

GramSchmidt(x)

pinv(x)

rref(x)

QRdecomposition(x)

det(x, ...)
}
\arguments{
\item{x}{A matrix for which a property is requested}

\item{...}{Auxillary arguments}
}
\value{
Returns the requested property of a matrix.
}
\description{
Performs various linear algebra operations like finding the inverse,
the QR decomposition, the eigenvectors and the eigenvalues.
}
\examples{
if (has_sympy()) {
  A <- matrix(c("a", "0", "0", "1"), 2, 2) \%>\% as_sym()
  
  QRdecomposition(A)
  eigenval(A)
  eigenvec(A)
  inv(A)
  det(A)
  
  
  A <- matrix(c("a", "b", "c", "d"), 2, 2) \%>\% as_sym()
  evec <- eigenvec(A)
  evec
  evec1 <- evec[[1]]$eigvec
  evec1
  simplify(evec1)
  
  lapply(evec, function(l) simplify(l$eigvec))

  A <- as_sym("[[1, 2, 3], [4, 5, 6]]")
  pinv(A)
}

}
\seealso{
\code{\link[=do_la]{do_la()}}
}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{diag<-.caracas_symbol}
\alias{diag<-.caracas_symbol}
\title{Replace diagonal}
\usage{
\method{diag}{caracas_symbol}(x) <- value
}
\arguments{
\item{x}{A \code{caracas_symbol}.}

\item{value}{Replacement value}
}
\description{
Replace diagonal
}
\examples{
if (has_sympy()) {
  A <- matrix(c("a", 0, 0, 0, "a", "a", "a", 0, 0), 3, 3)
  B <- as_sym(A)
  B
  diag(B)
  diag(B) <- "b"
  B
  diag(B)
}

}
\concept{vectors}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{subs_vec}
\alias{subs_vec}
\title{Substitute af vector of symbols for a vector of values}
\usage{
subs_vec(s, x, v)
}
\arguments{
\item{s}{Expression}

\item{x}{Names of symbol (vector)}

\item{v}{Values for \code{x} (vector)}
}
\description{
Substitute af vector of symbols for a vector of values
}
\examples{
if (has_sympy()) {
   x <- as_sym(paste0('x', 1:3))
   e <- 2*x^2
   e
   subs_vec(e, x, 1:3)
   subs_vec(e, x, x^2)
}

}
\seealso{
\code{\link[=subs]{subs()}}, \code{\link[=subs_lst]{subs_lst()}}
}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{matrix_}
\alias{matrix_}
\title{Symbolic matrix}
\usage{
matrix_(..., declare_symbols = TRUE)
}
\arguments{
\item{\dots}{Passed on to \code{\link[=matrix]{matrix()}}}

\item{declare_symbols}{Passed on to \code{as_sym()} when constructing symbolic matrix}
}
\description{
Symbolic matrix
}
\examples{
if (has_sympy()) {
  matrix_(1:9, nrow = 3)
  matrix_("a", 2, 2)
}

}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{install_sympy}
\alias{install_sympy}
\title{Install 'SymPy'}
\usage{
install_sympy(method = "auto", conda = "auto")
}
\arguments{
\item{method}{Installation method.
By default, "auto" automatically finds a method that will work
in the local environment.
Change the default to force a specific installation method.
Note that the "virtualenv" method is not available on Windows.}

\item{conda}{Path to conda executable (or "auto" to find conda
using the PATH and other conventional install locations).}
}
\value{
None
}
\description{
Install the 'SymPy' Python package into a
virtual environment or Conda environment.
}
\concept{sympy}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{der2}
\alias{der2}
\title{Symbolic differentiation of second order of an expression}
\usage{
der2(expr, vars, simplify = TRUE)
}
\arguments{
\item{expr}{A \code{caracas_symbol}}

\item{vars}{variables to take derivate with respect to}

\item{simplify}{Simplify result}
}
\description{
Symbolic differentiation of second order of an expression
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  y <- symbol("y")
  f <- 3*x^2 + x*y^2
  der2(f, x)
  h <- der2(f, list(x, y))
  h
  dim(h)
  H <- matrify(h)
  H
  dim(H)
}

}
\concept{calculus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{expand_trig}
\alias{expand_trig}
\title{Expand a trigonometric expression}
\usage{
expand_trig(x)
}
\arguments{
\item{x}{A \code{caracas_symbol}}
}
\description{
Expand a trigonometric expression
}
\concept{simplify}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/def_sym.R
\name{def_sym}
\alias{def_sym}
\title{Define caracas symbols in global environment}
\usage{
def_sym(..., charvec = NULL, warn = FALSE, env = parent.frame())
}
\arguments{
\item{...}{Names for new symbols, also supports non-standard evaluation}

\item{charvec}{Take each element in this character vector and define as caracas symbols}

\item{warn}{Warn if existing variable names are overwritten}

\item{env}{Environment to assign variable in}
}
\value{
Names of declared variables (invisibly)
}
\description{
Define caracas symbols in global environment
}
\examples{
if (has_sympy()) {
  ls()
  def_sym(n1, n2, n3)
  ls()
  def_sym("x1", "x2", "x3")
  ls()
  def_sym("x1", "x2", "x3", warn = TRUE)
  ls()
  def_sym(i, j, charvec = c("x", "y"))
  ls()
}

}
\seealso{
\code{\link[=symbol]{symbol()}}, \code{\link[=as_sym]{as_sym()}}
}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{as.character.caracas_symbol}
\alias{as.character.caracas_symbol}
\title{Convert symbol to character}
\usage{
\method{as.character}{caracas_symbol}(x, replace_I = TRUE, ...)
}
\arguments{
\item{x}{A \code{caracas_symbol}}

\item{replace_I}{Replace constant I (can both be identity and imaginary unit)}

\item{\dots}{not used}
}
\description{
Convert symbol to character
}
\concept{output}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ops-math.R
\name{Ops.caracas_symbol}
\alias{Ops.caracas_symbol}
\title{Math operators}
\usage{
\method{Ops}{caracas_symbol}(e1, e2)
}
\arguments{
\item{e1}{A \code{caracas_symbol}.}

\item{e2}{A \code{caracas_symbol}.}
}
\description{
Math operators
}
\concept{simple_algebra}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{listify}
\alias{listify}
\title{Convert object to list of elements}
\usage{
listify(x)
}
\arguments{
\item{x}{Object}
}
\description{
Convert object to list of elements
}
\examples{
if (has_sympy()) {
  x <- as_sym("Matrix([[b1*x1/(b2 + x1)], [b1*x2/(b2 + x2)], [b1*x3/(b2 + x3)]])")
  listify(x)
  
  xT <- t(x)
  listify(xT)
}

}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{tuplify}
\alias{tuplify}
\title{Convert object to tuple}
\usage{
tuplify(x)
}
\arguments{
\item{x}{Object}
}
\description{
Convert object to tuple
}
\examples{
if (has_sympy()) {
  x <- as_sym("Matrix([[b1*x1/(b2 + x1)], [b1*x2/(b2 + x2)], [b1*x3/(b2 + x3)]])")
  tuplify(x)
}

}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vector.R
\name{sum.caracas_symbol}
\alias{sum.caracas_symbol}
\title{Summation}
\usage{
\method{sum}{caracas_symbol}(..., na.rm = FALSE)
}
\arguments{
\item{\dots}{Elements to sum}

\item{na.rm}{Not used}
}
\description{
Summation
}
\concept{vectors}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ops-math.R
\name{Math.caracas_symbol}
\alias{Math.caracas_symbol}
\title{Math functions}
\usage{
\method{Math}{caracas_symbol}(x, ...)
}
\arguments{
\item{x}{\code{caracas_symbol}.}

\item{\dots}{further arguments passed to methods}
}
\description{
If \code{x} is a matrix, the function is applied component-wise.
}
\concept{simple_algebra}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset.R
\name{[.caracas_symbol}
\alias{[.caracas_symbol}
\title{Extract or replace parts of an object}
\usage{
\method{[}{caracas_symbol}(x, i, j, ..., drop = TRUE)
}
\arguments{
\item{x}{A \code{caracas_symbol}.}

\item{i}{row indices specifying elements to extract or replace}

\item{j}{column indices specifying elements to extract or replace}

\item{\dots}{Not used}

\item{drop}{Simplify dimensions of resulting object}
}
\description{
Extract or replace parts of an object
}
\examples{
if (has_sympy()) {
  A <- matrix(c("a", 0, 0, 0, "a", "a", "a", 0, 0), 3, 3)
  B <- as_sym(A)
  B[1:2, ]
  B[, 2]
  B[2, , drop = FALSE]
}

}
\concept{vectors}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{has_sympy}
\alias{has_sympy}
\title{Check if 'SymPy' is available}
\usage{
has_sympy()
}
\value{
\code{TRUE} if 'SymPy' is available, else \code{FALSE}
}
\description{
Check if 'SymPy' is available
}
\examples{
has_sympy()

}
\concept{sympy}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{prod_}
\alias{prod_}
\title{Product of a function}
\usage{
prod_(f, var, lower, upper, doit = TRUE)
}
\arguments{
\item{f}{Function to take product of}

\item{var}{Variable to take product for (either string or \code{caracas_symbol})}

\item{lower}{Lower limit}

\item{upper}{Upper limit}

\item{doit}{Evaluate the product immediately (or later with \code{\link[=doit]{doit()}})}
}
\description{
Product of a function
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  p <- prod_(1/x, "x", 1, 10)
  p
  as_expr(p)
  prod(1/(1:10))
  n <- symbol("n")
  prod_(x, x, 1, n)
}

}
\concept{calculus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{simplify}
\alias{simplify}
\title{Simplify expression}
\usage{
simplify(x)
}
\arguments{
\item{x}{A \code{caracas_symbol}}
}
\description{
Simplify expression
}
\concept{simplify}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{get_sympy}
\alias{get_sympy}
\title{Access 'SymPy' directly}
\usage{
get_sympy()
}
\value{
The 'SymPy' object with direct access to the library.
}
\description{
Get the 'SymPy' object.
Note that it gives you extra responsibilities
when you choose to access the 'SymPy' object directly.
}
\examples{
if (has_sympy()) {
  sympy <- get_sympy()
  sympy$solve("x**2-1", "x")
}

}
\concept{sympy}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe}
\arguments{
\item{lhs, rhs}{specify what lhs and rhs are}
}
\description{
Pipe operator
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{vec}
\alias{vec}
\title{Stacks matrix to vector}
\usage{
vec(x)
}
\arguments{
\item{x}{Matrix}
}
\description{
Stacks matrix to vector
}
\examples{
if (has_sympy()) {
  A <- as_sym(matrix(1:9, 3))
  vec(A)
}

}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset.R
\name{[<-.caracas_symbol}
\alias{[<-.caracas_symbol}
\title{Extract or replace parts of an object}
\usage{
\method{[}{caracas_symbol}(x, i, j, ...) <- value
}
\arguments{
\item{x}{A \code{caracas_symbol}.}

\item{i}{row indices specifying elements to extract or replace}

\item{j}{column indices specifying elements to extract or replace}

\item{\dots}{Not used}

\item{value}{Replacement value}
}
\description{
Extract or replace parts of an object
}
\examples{
if (has_sympy()) {
  A <- matrix(c("a", 0, 0, 0, "a", "a", "a", 0, 0), 3, 3)
  B <- as_sym(A)
  B[, 2] <- "x"
  B
}

}
\concept{vectors}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{print.caracas_solve_sys_sol}
\alias{print.caracas_solve_sys_sol}
\title{Print solution}
\usage{
\method{print}{caracas_solve_sys_sol}(
  x,
  simplify = getOption("caracas.print.sol.simplify", default = TRUE),
  ...
)
}
\arguments{
\item{x}{A \code{caracas_symbol}}

\item{simplify}{Print solution in a simple format}

\item{\dots}{Passed to \code{\link[=print.caracas_symbol]{print.caracas_symbol()}}}
}
\description{
Print solution
}
\concept{output}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{reciprocal_matrix}
\alias{reciprocal_matrix}
\title{Elementwise reciprocal matrix}
\usage{
reciprocal_matrix(x, numerator = 1)
}
\arguments{
\item{x}{Object \code{x}}

\item{numerator}{The numerator in the result.}
}
\description{
Elementwise reciprocal matrix
}
\examples{
if (has_sympy()) {
  s <- as_sym("[[r1, r2, r3], [u1, u2, u3]]")
  reciprocal_matrix(s, numerator = 7)
}

}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{get_py}
\alias{get_py}
\title{Access 'py' object}
\usage{
get_py()
}
\value{
The 'py' object with direct access to the library.
}
\description{
Get the 'py' object.
Note that it gives you extra responsibilities
when you choose to access the 'py' object directly.
}
\examples{
if (has_sympy()) {
  py <- get_py()
}

}
\concept{sympy}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{diag-set}
\alias{diag-set}
\alias{diag<-}
\title{Replace matrix diagonal}
\usage{
diag(x) <- value
}
\arguments{
\item{x}{Object \code{x}}

\item{value}{Replacement value}
}
\description{
Replace matrix diagonal
}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve.R
\name{solve_lin}
\alias{solve_lin}
\title{Solve a linear system of equations}
\usage{
solve_lin(A, b)
}
\arguments{
\item{A}{matrix}

\item{b}{vector}
}
\description{
Find \code{x} in \code{Ax = b}. If \code{b} not supplied,
the inverse of \code{A} is returned.
}
\concept{solve}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{int}
\alias{int}
\title{Integrate a function}
\usage{
int(f, var, lower, upper, doit = TRUE)
}
\arguments{
\item{f}{Function to integrate}

\item{var}{Variable to integrate with respect to (either string or \code{caracas_symbol})}

\item{lower}{Lower limit}

\item{upper}{Upper limit}

\item{doit}{Evaluate the integral immediately (or later with \code{\link[=doit]{doit()}})}
}
\description{
If no limits are provided, the
indefinite integral is calculated.
Otherwise, if both limits are provided,
the definite integral is calculated.
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  
  int(1/x, x, 1, 10)
  int(1/x, x, 1, 10, doit = FALSE)
  int(1/x, x)
  int(1/x, x, doit = FALSE)
  int(exp(-x^2/2), x, -Inf, Inf)
  int(exp(-x^2/2), x, -Inf, Inf, doit = FALSE)
}

}
\concept{calculus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{dim.caracas_symbol}
\alias{dim.caracas_symbol}
\title{Dimensions of a caracas symbol}
\usage{
\method{dim}{caracas_symbol}(x)
}
\arguments{
\item{x}{caracas symbol}
}
\description{
Dimensions of a caracas symbol
}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sh-additions.R
\name{all_vars}
\alias{all_vars}
\title{Get all variables used in a caracas symbol.}
\usage{
all_vars(sym)
}
\arguments{
\item{sym}{caracas symbol}
}
\description{
Get all variables used in a caracas symbol.
}
\examples{

all_vars(as_sym(paste0("v", 1:4)))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{diag_}
\alias{diag_}
\title{Symbolic diagonal matrix}
\usage{
diag_(x, n = 1L, declare_symbols = TRUE, ...)
}
\arguments{
\item{x}{Character vector with diagonal}

\item{n}{Number of times \code{x} should be repeated}

\item{declare_symbols}{Passed on to \code{as_sym()} when constructing symbolic matrix}

\item{\dots}{Passed on to \code{rep(x, n, ...)}}
}
\description{
Symbolic diagonal matrix
}
\examples{
if (has_sympy()) {
  diag_(c("a", "b", "c"))
  diag_("a", 2)
}

}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assumptions.R
\name{ask}
\alias{ask}
\title{Ask for a symbol's property}
\usage{
ask(x, property)
}
\arguments{
\item{x}{symbol}

\item{property}{property, e.g. 'positive'}
}
\description{
Ask for a symbol's property
}
\examples{
if (has_sympy()) {
  x <- symbol("x", positive = TRUE)
  ask(x, "positive")
}

}
\concept{assumptions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sh-additions.R
\name{eval_sym}
\alias{eval_sym}
\title{Evaluate caracas symbol.}
\usage{
eval_sym(sym, list_args)
}
\arguments{
\item{sym}{caracas symbol}

\item{list_args}{List of form name=value}
}
\description{
Evaluate caracas symbol after inserting values.
}
\examples{

x <- as_sym(rep(1:5))
x[2] <- "eps"
x
eval_sym(x, list(eps=10))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sh-additions.R
\name{generic-matrices}
\alias{generic-matrices}
\alias{n_vector}
\alias{nxn_matrix_diag}
\alias{nxn_matrix_symmetric}
\alias{nxm_matrix}
\alias{nxn_matrix}
\title{Generate generic vectors and matrices}
\usage{
n_vector(n, entry = "a")

nxn_matrix_diag(n, entry = "a")

nxn_matrix_symmetric(n, entry = "a")

nxm_matrix(n, m, entry = "a")

nxn_matrix(n, entry = "a", symmetric = FALSE)
}
\arguments{
\item{n, m}{Dimensions.}

\item{entry}{The symbolic name of each entry.}

\item{symmetric}{Logical}
}
\description{
Generate generic vectors and matrices
}
\examples{

n_vector(5, "b")
nxm_matrix(3, 2, "a")
nxn_matrix(5, "s", sym=TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sh-additions.R
\name{as_sym2}
\alias{as_sym2}
\title{Coerce (R) expression into caracas symbol.}
\usage{
as_sym2(expr, list_args = NULL)
}
\arguments{
\item{expr}{Expression (R)}

\item{list_args}{List of form name=value}
}
\description{
Coerce (R) expression into caracas symbol; possibly after substituting values into expression
}
\examples{

expr <- expression(1-a^d)

as_sym2(expr)
as_sym2(expr, list(d=8))
as_sym2("1-a^d")
as_sym2("1-a^d", list(d=8))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lin-alg.R
\name{diag}
\alias{diag}
\title{Matrix diagonal}
\usage{
diag(x, ...)
}
\arguments{
\item{x}{Object \code{x}}

\item{\dots}{Passed on}
}
\description{
Matrix diagonal
}
\concept{linalg}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{doit}
\alias{doit}
\title{Perform calculations setup previously}
\usage{
doit(x)
}
\arguments{
\item{x}{A \code{caracas_symbol}}
}
\description{
Perform calculations setup previously
}
\examples{
if (has_sympy()) {
   x <- symbol('x')
   res <- lim(sin(x)/x, "x", 0, doit = FALSE)
   res 
   doit(res)
}

}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{expand}
\alias{expand}
\title{Expand expression}
\usage{
expand(x)
}
\arguments{
\item{x}{A \code{caracas_symbol}}
}
\description{
Expand expression
}
\concept{simplify}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{print.caracas_symbol}
\alias{print.caracas_symbol}
\title{Print symbol}
\usage{
\method{print}{caracas_symbol}(
  x,
  caracas_prefix = TRUE,
  prettyascii = getOption("caracas.print.prettyascii", default = FALSE),
  ascii = getOption("caracas.print.ascii", default = FALSE),
  rowvec = getOption("caracas.print.rowvec", default = TRUE),
  ...
)
}
\arguments{
\item{x}{A \code{caracas_symbol}}

\item{caracas_prefix}{Print 'caracas' prefix}

\item{prettyascii}{\code{TRUE} to print in pretty ASCII format rather than in utf8}

\item{ascii}{\code{TRUE} to print in ASCII format rather than in utf8}

\item{rowvec}{\code{FALSE} to print column vectors as is}

\item{\dots}{not used}
}
\description{
Print symbol
}
\concept{output}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{lim}
\alias{lim}
\title{Limit of a function}
\usage{
lim(f, var, val, dir = NULL, doit = TRUE)
}
\arguments{
\item{f}{Function to take limit of}

\item{var}{Variable to take limit for (either string or \code{caracas_symbol})}

\item{val}{Value for \code{var} to approach}

\item{dir}{Direction from where \code{var} should approach \code{val}: \code{'+'} or \code{'-'}}

\item{doit}{Evaluate the limit immediately (or later with \code{\link[=doit]{doit()}})}
}
\description{
Limit of a function
}
\examples{
if (has_sympy()) {
  x <- symbol("x")
  lim(sin(x)/x, "x", 0)
  lim(1/x, "x", 0, dir = '+')
  lim(1/x, "x", 0, dir = '-')
}

}
\concept{calculus}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_python.R
\name{as_expr}
\alias{as_expr}
\title{Convert caracas object to R}
\usage{
as_expr(x, first_doit = TRUE)
}
\arguments{
\item{x}{caracas_symbol}

\item{first_doit}{Try \code{doit()} first}
}
\description{
Potentially calls \code{\link[=doit]{doit()}}.
}
\concept{caracas_symbol}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symbol.R
\name{fraction_parts}
\alias{fraction_parts}
\title{Get numerator and denominator of a fraction}
\usage{
fraction_parts(x)
}
\arguments{
\item{x}{Fraction}
}
\description{
Get numerator and denominator of a fraction
}
\examples{
if (has_sympy()) {
     x <- as_sym("a/b")
     frac <- fraction_parts(x)
     frac
     frac$numerator
     frac$denominator
 }

}
\concept{caracas_symbol}
