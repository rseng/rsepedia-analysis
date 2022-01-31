# autotest <a href='https://docs.ropensci.org/autotest'><img src='man/figures/autotest.png' align="right" height=210 width=182></a>

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R build
status](https://github.com/ropensci-review-tools/autotest/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-review-tools/autotest/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci-review-tools/autotest/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci-review-tools/autotest)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->

Automatic mutation testing of R packages. Mutation in the sense of
mutating inputs (parameters) to function calls. `autotest` primarily
works by scraping documented examples for all functions, and mutating
the parameters input to those functions.

## Installation

The easiest way to install this package is via the associated
[`r-universe`](https://ropensci-review-tools.r-universe.dev/ui#builds).
As shown there, simply enable the universe with

``` r
options(repos = c(
    ropenscireviewtools = "https://ropensci-review-tools.r-universe.dev",
    CRAN = "https://cloud.r-project.org"))
```

And then install the usual way with,

``` r
install.packages("autotest")
```

Alternatively, the package can be installed by running one of the
following lines:

``` r
# install.packages("remotes")
remotes::install_git("https://git.sr.ht/~mpadge/autotest")
remotes::install_bitbucket("mpadge/autotest")
remotes::install_gitlab("mpadge/autotest")
remotes::install_github("ropensci-review-tools/autotest")
```

The package can then be loaded the usual way:

``` r
library (autotest)
```

## Usage

The simply way to use the package is

``` r
x <- autotest_package ("<package>")
```

The main argument to the [`autotest_package()`
function](https://docs.ropensci.org/autotest/reference/autotest_package.html)
can either be the name of an installed package, or a path to a local
directory containing the source for a package. The result is a
`data.frame` of errors, warnings, and other diagnostic messages issued
during package `auotest`-ing. The function has an additional parameter,
`functions`, to restrict tests to specified functions only.

By default,
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html)
returns a list of all tests applied to a package without actually
running them. To implement those tests, set the parameter `test` to
`TRUE`. Results are only returned for tests in which functions do not
behave as expected, whether through triggering errors, warnings, or
other behaviour as described below. The ideal behaviour of
`autotest_package()` is to return nothing (or strictly, `NULL`),
indicating that all tests passed successfully. See the [main package
vignette](https://docs.ropensci.org/autotest/articles/autotest.html) for
an introductory tour of the package.

## What is tested?

The package includes a function which lists all tests currently
implemented.

``` r
autotest_types ()
#> # A tibble: 27 × 8
#>    type  test_name  fn_name parameter parameter_type operation   content   test 
#>    <chr> <chr>      <chr>   <chr>     <chr>          <chr>       <chr>     <lgl>
#>  1 dummy rect_as_o… <NA>    <NA>      rectangular    Convert on… "check f… TRUE 
#>  2 dummy rect_comp… <NA>    <NA>      rectangular    Convert on… "expect … TRUE 
#>  3 dummy rect_comp… <NA>    <NA>      rectangular    Convert on… "expect … TRUE 
#>  4 dummy rect_comp… <NA>    <NA>      rectangular    Convert on… "expect … TRUE 
#>  5 dummy extend_re… <NA>    <NA>      rectangular    Extend exi… "(Should… TRUE 
#>  6 dummy replace_r… <NA>    <NA>      rectangular    Replace cl… "(Should… TRUE 
#>  7 dummy vector_to… <NA>    <NA>      vector         Convert ve… "(Should… TRUE 
#>  8 dummy vector_cu… <NA>    <NA>      vector         Custom cla… "(Should… TRUE 
#>  9 dummy double_is… <NA>    <NA>      numeric        Check whet… "int par… TRUE 
#> 10 dummy trivial_n… <NA>    <NA>      numeric        Add trivia… "(Should… TRUE 
#> # … with 17 more rows
```

That functions returns a [`tibble`](https://tibble.tidyverse.org)
describing 27 unique tests. The default behaviour of
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html)
with `test = FALSE` uses these test types to identify which tests will
be applied to each parameter and function. The table returned from
[`autotest_types()`](https://docs.ropensci.org/autotest/reference/autotest_types.html)
can be used to selectively switch tests off by setting values in the
`test` column to `FALSE`, as demonstrated below.

## How Does It Work?

The package works by scraping documented examples from all `.Rd` help
files, and using those to identify the types of all parameters to all
functions. Usage therefore first requires that the usage of all
parameters be demonstrated in example code.

As described above, tests can also be selectively applied to particular
functions through the parameters `functions`, used to nominate functions
to include in tests, or `exclude`, used to nominate functions to exclude
from tests. The following code illustrates.

``` r
x <- autotest_package (package = "stats", functions = "var", test = FALSE)
#> 
#> ── autotesting stats ──
#> 
#> ✔ [1 / 6]: var
#> ✔ [2 / 6]: cor
#> ✔ [3 / 6]: cor
#> ✔ [4 / 6]: cov
#> ✔ [5 / 6]: cov
#> ✔ [6 / 6]: cor
print (x)
#> # A tibble: 185 × 9
#>    type    test_name  fn_name parameter parameter_type operation  content  test 
#>    <chr>   <chr>      <chr>   <chr>     <chr>          <chr>      <chr>    <lgl>
#>  1 warning par_is_de… var     use       <NA>           Check tha… Example… TRUE 
#>  2 warning par_is_de… cov     y         <NA>           Check tha… Example… TRUE 
#>  3 dummy   trivial_n… var     x         numeric        Add trivi… (Should… TRUE 
#>  4 dummy   vector_cu… var     x         vector         Custom cl… (Should… TRUE 
#>  5 dummy   vector_to… var     x         vector         Convert v… (Should… TRUE 
#>  6 dummy   negate_lo… var     na.rm     single logical Negate de… (Functi… TRUE 
#>  7 dummy   subst_int… var     na.rm     single logical Substitut… (Functi… TRUE 
#>  8 dummy   subst_cha… var     na.rm     single logical Substitut… should … TRUE 
#>  9 dummy   single_pa… var     na.rm     single logical Length 2 … Should … TRUE 
#> 10 dummy   return_su… var     (return … (return objec… Check tha… <NA>     TRUE 
#> # … with 175 more rows, and 1 more variable: yaml_hash <chr>
```

Testing the `var` function also tests `cor` and `cov`, because these are
all documented within a single `.Rd` help file. Typing `?var` shows that
the help topic is `cor`, and that the examples include the three
functions, `var`, `cor`, and `cov`. That result details the 185 tests
which would be applied to the `var` function from the `stats` package.
These 185 tests yield the following results when actually applied:

``` r
y <- autotest_package (package = "stats", functions = "var", test = TRUE)
#> ── autotesting stats ──
#> 
#> ✔ [1 / 6]: var
#> ✔ [2 / 6]: cor
#> ✔ [3 / 6]: cor
#> ✔ [4 / 6]: cov
#> ✔ [5 / 6]: cov
#> ✔ [6 / 6]: cor
print (y)
#> # A tibble: 23 × 9
#>    type       test_name fn_name parameter parameter_type operation content test 
#>    <chr>      <chr>     <chr>   <chr>     <chr>          <chr>     <chr>   <lgl>
#>  1 warning    par_is_d… var     use       <NA>           Check th… "Examp… TRUE 
#>  2 warning    par_is_d… cov     y         <NA>           Check th… "Examp… TRUE 
#>  3 diagnostic vector_t… var     x         vector         Convert … "Funct… TRUE 
#>  4 diagnostic vector_t… var     x         vector         Convert … "Funct… TRUE 
#>  5 diagnostic vector_t… var     y         vector         Convert … "Funct… TRUE 
#>  6 diagnostic single_c… cor     use       single charac… upper-ca… "is ca… TRUE 
#>  7 diagnostic single_c… cor     method    single charac… upper-ca… "is ca… TRUE 
#>  8 diagnostic vector_c… cor     x         vector         Custom c… "Funct… TRUE 
#>  9 diagnostic single_c… cor     method    single charac… upper-ca… "is ca… TRUE 
#> 10 diagnostic single_c… cor     use       single charac… upper-ca… "is ca… TRUE 
#> # … with 13 more rows, and 1 more variable: yaml_hash <chr>
```

And only 23 of the original 185 tests produced unexpected behaviour.
There were in fact only 4 kinds of tests which produced these 23
results:

``` r
unique (y$operation)
#> [1] "Check that parameter usage is demonstrated"
#> [2] "Convert vector input to list-columns"      
#> [3] "upper-case character parameter"            
#> [4] "Custom class definitions for vector input"
```

One of these involves conversion of a vector to a list-column
representation (via `I(as.list(<vec>))`). Relatively few packages accept
this kind of input, even though doing so is relatively straightforward.
The following lines demonstrate how these tests can be switched off when
`autotest`-ing a package. The `autotest_types()` function, used above to
extract information on all types of tests, also accepts a single
argument listing the `test_name` entries of any tests which are to be
switched off.

``` r
types <- autotest_types (notest = "vector_to_list_col")
y <- autotest_package (package = "stats", functions = "var",
                       test = TRUE, test_data = types)
#> ── autotesting stats ──
#> 
#> ✔ [1 / 6]: var
#> ✔ [2 / 6]: cor
#> ✔ [3 / 6]: cor
#> ✔ [4 / 6]: cov
#> ✔ [5 / 6]: cov
#> ✔ [6 / 6]: cor
print (y)
#> # A tibble: 28 × 9
#>    type       test_name fn_name parameter parameter_type operation content test 
#>    <chr>      <chr>     <chr>   <chr>     <chr>          <chr>     <chr>   <lgl>
#>  1 warning    par_is_d… var     use       <NA>           Check th… Exampl… TRUE 
#>  2 warning    par_is_d… cov     y         <NA>           Check th… Exampl… TRUE 
#>  3 diagnostic single_c… cor     use       single charac… upper-ca… is cas… TRUE 
#>  4 diagnostic single_c… cor     method    single charac… upper-ca… is cas… TRUE 
#>  5 diagnostic vector_c… cor     x         vector         Custom c… Functi… TRUE 
#>  6 diagnostic single_c… cor     method    single charac… upper-ca… is cas… TRUE 
#>  7 diagnostic single_c… cor     use       single charac… upper-ca… is cas… TRUE 
#>  8 diagnostic single_c… cor     use       single charac… upper-ca… is cas… TRUE 
#>  9 diagnostic single_c… cor     method    single charac… upper-ca… is cas… TRUE 
#> 10 diagnostic vector_c… cov     x         vector         Custom c… Functi… TRUE 
#> # … with 18 more rows, and 1 more variable: yaml_hash <chr>
```

Those tests are still returned from `autotest_package()`, but with
`test = FALSE` to indicate they were not run, and a `type` of “no_test”
rather than the previous “diagnostic”.

## Can `autotest` automatically create tests in my `tests` directory?

Not yet, but that should be possible soon. In the meantime, there are
[`testthat`](https://testthat.r-lib.org) expectations, listed in the
[main package
functions](https://docs.ropensci.org/autotest/reference/index.html),
which enable `autotest` to be used in a package’s test suite.

## Prior work

1.  The
    [`great-expectations`](https://github.com/great-expectations/great_expectations)
    framework for python, described in [this medium
    article](https://medium.com/@expectgreatdata/down-with-pipeline-debt-introducing-great-expectations-862ddc46782a).
2.  [`QuickCheck`](https://hackage.haskell.org/package/QuickCheck) for
    Haskell
3.  [`mutate`](https://github.com/mbj/mutant) for ruby
4.  [`mutant`](https://github.com/ropensci/mutant) for mutation of R
    code itself

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
---
title: "autotest"
output:
  md_document:
    variant: gfm

  rmarkdown::html_vignette:
    self_contained: no
---

# autotest <a href='https://docs.ropensci.org/autotest'><img src='man/figures/autotest.png' align="right" height=210 width=182></a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = TRUE,
  message = TRUE,
  width = 120,
  comment = "#>",
  fig.retina = 2,
  fig.path = "README-"
)
```

<!-- badges: start -->

[![R build
status](https://github.com/ropensci-review-tools/autotest/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-review-tools/autotest/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci-review-tools/autotest/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci-review-tools/autotest)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->


Automatic mutation testing of R packages. Mutation in the sense of mutating
inputs (parameters) to function calls. `autotest` primarily works by scraping
documented examples for all functions, and mutating the parameters input to
those functions.


## Installation

The easiest way to install this package is via the associated
[`r-universe`](https://ropensci-review-tools.r-universe.dev/ui#builds). As
shown there, simply enable the universe with

```{r options, eval = FALSE}
options(repos = c(
    ropenscireviewtools = "https://ropensci-review-tools.r-universe.dev",
    CRAN = "https://cloud.r-project.org"))
```

And then install the usual way with,

```{r install, eval = FALSE}
install.packages("autotest")
```

Alternatively, the package can be installed by running one of the following
lines:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_git("https://git.sr.ht/~mpadge/autotest")
remotes::install_bitbucket("mpadge/autotest")
remotes::install_gitlab("mpadge/autotest")
remotes::install_github("ropensci-review-tools/autotest")
```

The package can then be loaded the usual way:
```{r load-fakey, eval = FALSE}
library (autotest)
```
```{r load, echo = FALSE, message = FALSE}
devtools::load_all (".", export_all = FALSE)
```

## Usage

The simply way to use the package is

```{r autotest-example, eval = FALSE}
x <- autotest_package ("<package>")
```

The main argument to the [`autotest_package()`
function](https://docs.ropensci.org/autotest/reference/autotest_package.html)
can either be the name of an installed package, or a path to a local directory
containing the source for a package. The result is a `data.frame` of errors,
warnings, and other diagnostic messages issued during package `auotest`-ing.
The function has an additional parameter, `functions`, to restrict tests to
specified functions only.

By default,
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html)
returns a list of all tests applied to a package without actually running them.
To implement those tests, set the parameter `test` to `TRUE`. Results are only returned for tests 
in which functions do not behave as expected, whether through triggering
errors, warnings, or other behaviour as described below. The ideal behaviour of
`autotest_package()` is to return nothing (or strictly, `NULL`), indicating
that all tests passed successfully. See the [main package
vignette](https://docs.ropensci.org/autotest/articles/autotest.html) for an
introductory tour of the package.

## What is tested?

The package includes a function which lists all tests currently implemented.

```{r autotest_types}
autotest_types ()
```

That functions returns a [`tibble`](https://tibble.tidyverse.org) describing
`r nrow(autotest_types())` unique tests. The default behaviour of
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html)
with `test = FALSE` uses these test types to identify which tests will be
applied to each parameter and function. The table returned from
[`autotest_types()`](https://docs.ropensci.org/autotest/reference/autotest_types.html)
can be used to selectively switch tests off by setting values in the `test`
column to `FALSE`, as demonstrated below.

## How Does It Work?

The package works by scraping documented examples from all `.Rd` help files,
and using those to identify the types of all parameters to all functions. Usage
therefore first requires that the usage of all parameters be demonstrated in
example code.

As described above, tests can also be selectively applied to particular
functions through the parameters `functions`, used to nominate functions to
include in tests, or `exclude`, used to nominate functions to exclude from
tests. The following code illustrates.

```{r stats-var-no-test, fig.show = "hide"}
x <- autotest_package (package = "stats", functions = "var", test = FALSE)
print (x)
```

Testing the `var` function also tests `cor` and `cov`, because these are all
documented within a single `.Rd` help file. Typing `?var` shows that the help
topic is `cor`, and that the examples include the three functions, `var`,
`cor`, and `cov`. That result details the `r nrow (x)` tests which would be
applied to the `var` function from the `stats` package. These `r nrow (x)`
tests yield the following results when actually applied:

```{r stats-var-test}
y <- autotest_package (package = "stats", functions = "var", test = TRUE)
print (y)
```

And only `r nrow (y)` of the original `r nrow (x)` tests produced unexpected
behaviour. There were in fact only `r length (unique (y$operation))` kinds of
tests which produced these `r nrow (y)` results:


```{r unique-operations}
unique (y$operation)
```

One of these involves conversion of a vector to a list-column representation
(via `I(as.list(<vec>))`). Relatively few packages accept this kind of input,
even though doing so is relatively straightforward. The following lines
demonstrate how these tests can be switched off when `autotest`-ing a package.
The `autotest_types()` function, used above to extract information on all types
of tests, also accepts a single argument listing the `test_name` entries of any
tests which are to be switched off.

```{r stats-var-test-switch}
types <- autotest_types (notest = "vector_to_list_col")
y <- autotest_package (package = "stats", functions = "var",
                       test = TRUE, test_data = types)
print (y)
```

Those tests are still returned from `autotest_package()`, but with `test =
FALSE` to indicate they were not run, and a `type` of "no_test" rather than the
previous "diagnostic".


## Can `autotest` automatically create tests in my `tests` directory?

Not yet, but that should be possible soon. In the meantime, there are
[`testthat`](https://testthat.r-lib.org) expectations, listed in the [main
package
functions](https://docs.ropensci.org/autotest/reference/index.html),
which enable `autotest` to be used in a package's test suite.


## Prior work

1. The
   [`great-expectations`](https://github.com/great-expectations/great_expectations)
   framework for python, described in [this medium
   article](https://medium.com/@expectgreatdata/down-with-pipeline-debt-introducing-great-expectations-862ddc46782a).
2. [`QuickCheck`](https://hackage.haskell.org/package/QuickCheck) for Haskell
3. [`mutate`](https://github.com/mbj/mutant) for ruby
4. [`mutant`](https://github.com/ropensci/mutant) for mutation of R code itself

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

# hex sticker script

Images from
https://www.onlinewebfonts.com/icon/

Modify images in gimp by:
1. Layer -> Transparent -> Add Alpha Channel
2. (Fuzzy) Select bits to delete + select -> delete all background
3. Export as png with no background

Can then convert to eps if desired, but direct loading png below also works the
same. All black colours were adjusted to gray 0.8

```{r}
library (ggplot2)
# trace outline of hexagon from centre bottom point in anti-clockwise direction
s3 <- sqrt (3) / 2
border <- data.frame (x = 1 + c (rep (-s3, 2), 0, rep (s3, 2), 0, -s3),
                      y = 1 + c (0.5, -0.5, -1, -0.5, 0.5, 1, 0.5))
asp <- diff (range (border$x)) / diff (range (border$y)) # aspect ratio for image

f <- "clipboard.svg"
d <- data.frame (x = 1, y = 1, image = f)

col_fore <- "#444444"
col_back <- "#888888"
col_fill <- "#11DCDC"

hex <- ggplot() +
    geom_polygon (aes_ (x = ~x, y = ~y), data = border,
                 size = 6, fill = col_fill, color = col_back) +
    ggimage::geom_image (aes_ (x = ~x, y = ~y, image = ~image), d,
                         size = 0.65, asp = 1) +
    scale_x_continuous (expand = c (0.02, 0.02)) +
    scale_y_continuous (expand = c (0.02, 0.02))
print (hex)
```
```{r}
add_one_lab <- function (hex, lab_dat, aes, fs) {

    hex <- hex + ggplot2::geom_text (dat = lab_dat,
                                     mapping = aes,
                                     size = fs,
                                     colour = col_back,
                                     family = 'SF Alien Encounters', 
                                     fontface = 1,
                                     nudge_y = -0.02,
                                     nudge_x = 0.02)
    hex <- hex + ggplot2::geom_text (dat = lab_dat,
                                     mapping = aes,
                                     size = fs,
                                     colour = col_fore,
                                     fontface = 1,
                                     family = 'SF Alien Encounters')
    return (hex)
}

lab_dat <- data.frame (x = 1 - 0.0001,
                       y = 1.3 + 0.0001,
                       lab = 'auto')
aes <- ggplot2::aes (x, y, label = lab)
fs <- 48 # font size
hex <- add_one_lab (hex, lab_dat, aes, fs)

lab_dat <- data.frame (x = 1 - 0.0001,
                       y = 0.7 + 0.0001,
                       lab = 'test')
aes <- ggplot2::aes (x, y, label = lab)
fs <- 48 # font size
hex <- add_one_lab (hex, lab_dat, aes, fs)


th <- theme_minimal ()
th$panel.background <- element_rect (fill = "transparent", size = 0)
th$line <- element_blank ()
th$axis.text <- element_blank ()
th$axis.title <- element_blank ()
th$plot.margin <- margin (rep (unit (0, 'null'), 4))
#th$plot.margin <- margin (rep (unit (-0.5, 'line'), 4))
th$legend.position <- 'none'
th$axis.ticks.length <- unit (0, 'null')

hex <- hex + th
print (hex)
```



from https://databasefaq.com/index.php/answer/182192/r-fonts-ggplot2-eps-error-using-arial-in-eps-figure-with-extrafont-package
```{r}
fname <- "autotest.svg"
library(showtext)
## add the Arial font
font_add("SF Alien Encounters",
         regular = "SFAlienEncounters.ttf",
         bold = "SFAlienEncountersSolid.ttf",
         italic = "SFAlienEncounters-Italic.ttf",
         bolditalic = "SFAlienEncountersSolid-Ital.ttf")

setEPS()
postscript(fname)
showtext_begin() ## call this function after opening a device
hex + theme_minimal (base_family = "SF Alien Encounters") +
    theme (axis.line = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank())
dev.off()
ggsave (fname, hex)
ggsave ("autotest.png", hex)
```
---
title: "How to use autotest"
author: 
  - "Mark Padgham"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        number_sections: false
        theme: flatly
vignette: >
  %\VignetteIndexEntry{How to use autotest}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = TRUE,
  message = TRUE,
  width = 120,
  comment = "#>",
  fig.retina = 2,
  fig.path = "README-"
)
```

This vignette demonstrates the easiest way to use `autotest`, which is to apply
it continuously through the entire process of package development. The best way
to understand the process is to obtain a local copy of the vignette itself from
[this
link](https://github.com/ropensci-review-tools/autotest/blob/master/vignettes/autotest.Rmd),
and step through the code. We begin by constructing a simple package in the
local
[`tempdir()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/tempfile.html).

<details closed>

<summary> <span title="Click to Expand"> Package Construction </span> </summary>

To create a package in one simple line, we use
[`usethis::create_package()`](https://usethis.r-lib.org/reference/create_package.html),
and name our package `"demo"`.

```{r create_package}
path <- file.path (tempdir (), "demo")
usethis::create_package (path, check_name = FALSE, open = FALSE)
```

The structure looks like this:

```{r dir_tree}
fs::dir_tree (path)
```

</details><br>

Having constructed a minimal package structure, we can then insert some code in
the `R/` directory, including initial [`roxygen2`](https://roxygen2.r-lib.org)
documentation lines, and use the [`roxygenise()`
function](https://roxygen2.r-lib.org/reference/roxygenize.html) to create the
corresponding `man` files.

```{r first-fn}
code <- c ("#' my_function",
           "#'",
           "#' @param x An input",
           "#' @return Something else",
           "#' @export",
           "my_function <- function (x) {",
           "  return (x + 1)",
           "}")
writeLines (code, file.path (path, "R", "myfn.R"))
roxygen2::roxygenise (path)
```

Our package now looks like this:

```{r dir_tree2}
fs::dir_tree (path)
```

We can already apply `autotest` to that package to see what happens, first
ensuring that we've loaded the package ready to use.

```{r autotest1-fakey, eval = FALSE, echo = TRUE}
library (autotest)
x0 <- autotest_package (path)
```
```{r autotest1, eval = TRUE, echo = FALSE}
devtools::load_all (".", export_all = FALSE)
x0 <- autotest_package (path)
```

We use the [`DT` package](https://rstudio.github.io/DT) to display the results
here.

```{r}
DT::datatable (x0, options = list (dom = "t")) # display table only
```

## Adding examples to our code

That tells us straight away that we need to add an example to our function
documentation. So let's do that by inserting extra lines into the `code`
defined above, and see what happens.

```{r autotest-FALSE}
code <- c (code [1:4],
           "#' @examples",
           "#' y <- my_function (x = 1)",
           code [5:length (code)])
writeLines (code, file.path (path, "R", "myfn.R"))
roxygen2::roxygenise (path)
x1 <- autotest_package (path)
DT::datatable (x1, options = list (dom = "t"))
```

The first thing to notice is the first column, which has `test_type = "dummy"`
for all rows. The [`autotest_package()`
function](https://docs.ropensci.org/autotest/reference/autotest_package.html)
has a parameter `test` with a default value of `FALSE`, so that the default
call demonstrated above does not actually implement the tests, rather it
returns an object listing all tests that would be performed with actually doing
so. Applying the tests by setting `test = TRUE` gives the following result.

```{r autotest-TRUE}
x2 <- autotest_package (path, test = TRUE)
DT::datatable (x2, options = list (dom = "t"))
```

Of the `r nrow(x1)` tests which were performed, only `r nrow(x2)` yielded
unexpected behaviour. The first indicates that the parameter `x` has only been
used as an integer, yet was not specified as such. The second states that the
parameter `x` is "assumed to be a single numeric". `autotest` does its best to
figure out what types of inputs are expected for each parameter, and with the
example only demonstrating `x = 1`, assumes that `x` is always expected to be
a single value. We can resolve the first of these by replacing `x = 1` with `x
= 1.` to clearly indicate that it is not an integer, and the second by
asserting that `length(x) == 1`, as follows:

```{r assert-length}
code <- c ("#' my_function",
           "#'",
           "#' @param x An input",
           "#' @return Something else",
           "#' @examples",
           "#' y <- my_function (x = 1.)",
           "#' @export",
           "my_function <- function (x) {",
           "  if (length(x) > 1) {",
           "    warning(\"only the first value of x will be used\")",
           "    x <- x [1]",
           "  }",
           "  return (x + 1)",
           "}")
writeLines (code, file.path (path, "R", "myfn.R"))
roxygen2::roxygenise (path)
```

This is then sufficient to pass all `autotest` tests and so return `NULL`.

```{r autotest-TRUE2}
autotest_package (path, test = TRUE)
```

## Integer input

Note that `autotest` distinguishes integer and non-integer types by their
[`storage.mode`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mode.html)
of `"integer"` and `"double"`, and not by their respective classes of
`"integer"` and `"numeric"`, because `"numeric"` is ambiguous in R, and
`is.numeric(1L)` is `TRUE`, even though `storage.mode(1L)` is `"integer"`, and
not `"numeric"`. Replacing `x = 1` with `x = 1.` explicitly identifies that
parameter as a `"double"` parameter, and allowed the preceding tests to pass.
Note what happens if we instead specify that parameter as an integer (`x =
1L`).

```{r int-input}
code [6] <- gsub ("1\\.", "1L", code [6])
writeLines (code, file.path (path, "R", "myfn.R"))
roxygen2::roxygenise (path)
x3 <- autotest_package (path, test = TRUE)
DT::datatable (x3, options = list (dom = "t"))
```

That then generates two additional messages, the second of which reflects an
expectation that parameters assumed to be integer-valued should assert that,
for example by converting with `as.integer()`. The following suffices to remove
that message.

```{r use-as-int}
code <- c (code [1:12],
           "  if (is.numeric (x))",
           "    x <- as.integer (x)",
           code [13:length (code)])
```

The remaining message concerns integer ranges. For any parameters which
`autotest` identifies as single integers, routines will try a full range of
values between `+/- .Machine$integer.max`, to ensure that all values are
appropriately handled. Many routines may sensibly allow unrestricted ranges,
while many others may not implement explicit control over permissible ranges,
yet may error on, for example, unexpectedly large positive or negative values.
The content of the diagnostic message indicates one way to resolve this issue,
which is simply by describing the input as `"unrestricted"`.

```{r unrestricted}
code [3] <- gsub ("An input", "An unrestricted input", code [3])
writeLines (code, file.path (path, "R", "myfn.R"))
roxygen2::roxygenise (path)
autotest_package (path, test = TRUE)
```

An alternative, and frequently better way, is to ensure and document specific
control over permissible ranges, as in the following revision of our function.

```{r input-range}
code <- c ("#' my_function",
           "#'",
           "#' @param x An input between 0 and 10",
           "#' @return Something else",
           "#' @examples",
           "#' y <- my_function (x = 1L)",
           "#' @export",
           "my_function <- function (x) {",
           "  if (length(x) > 1) {",
           "    warning(\"only the first value of x will be used\")",
           "    x <- x [1]",
           "  }",
           "  if (is.numeric (x))",
           "    x <- as.integer (x)",
           "  if (x < 0 | x > 10) {",
           "    stop (\"x must be between 0 and 10\")",
           "  }",
           "  return (x + 1L)",
           "}")
writeLines (code, file.path (path, "R", "myfn.R"))
roxygen2::roxygenise (path)
autotest_package (path, test = TRUE)
```



## Vector input

The initial test results above suggested that the input was *assumed* to be of
length one. Let us now revert our function to its original format which
accepted vectors of length > 1, and include an example demonstrating such
input.

```{r vector-input}
code <- c ("#' my_function",
           "#'",
           "#' @param x An input",
           "#' @return Something else",
           "#' @examples",
           "#' y <- my_function (x = 1)",
           "#' y <- my_function (x = 1:2)",
           "#' @export",
           "my_function <- function (x) {",
           "  if (is.numeric (x)) {",
           "    x <- as.integer (x)",
           "  }",
           "  return (x + 1L)",
           "}")
writeLines (code, file.path (path, "R", "myfn.R"))
roxygen2::roxygenise (path)
```

Note that the first example no longer has `x = 1L`. This is because vector
inputs are identified as `integer` by examining all individual values, and
presuming `integer` representations for any parameters for which all values are
whole numbers, regardless of `storage.mode`.

```{r autotest-TRUE3}
x4 <- autotest_package (path, test = TRUE)
DT::datatable (x4, options = list (dom = "t"))
```

### List-column conversion

The above result reflects one of the standard tests, which is to determine
whether list-column formats are appropriately processed. List-columns commonly
arise when using (either directly or indirectly), the [`tidyr::nest()`
function](https://tidyr.tidyverse.org/reference/nest.html), or equivalently in
base R with the [`I` or `AsIs`
function](https://stat.ethz.ch/R-manual/R-devel/library/base/html/AsIs.html).
They look like this:

```{r list-col-demo}
dat <- data.frame (x = 1:3, y = 4:6)
dat$x <- I (as.list (dat$x)) # base R
dat <- tidyr::nest (dat, y = y)
print (dat)
```

The use of packages like [`tidyr`](https://tidyr.tidyverse.org) and
[`purrr`](https://purrr.tidyverse.org) quite often leads to
[`tibble`](https://tibble.tidyverse.org)-class inputs which contain
list-columns. Any functions which fail to identify and appropriately respond to
such inputs may generate unexpected errors, and this `autotest` is intended to
enforce appropriate handling of these kinds of inputs. The following lines
demonstrate the kinds of results that can arise without such checks.

```{r mtcars-error, error = TRUE}
m <- mtcars
head (m, n = 2L)
m$mpg <- I (as.list (m$mpg))
head (m, n = 2L) # looks exaxtly the same
cor (m)
```

In contrast, many functions either assume inputs to be lists, and convert when
not, or implicitly `unlist`. Either way, such functions may respond entirely
consistently regardless of the presence of list-columns, like this:

```{r mtcars-okay}
m$mpg <- paste0 ("a", m$mpg)
class (m$mpg)
```

The list-column `autotest` is intended to enforce consistent behaviour in
response to list-column inputs. One way to identify list-column formats is to
check the value of `class(unclass(.))` of each column. The `unclass` function
is necessary to first remove any additional class attributes, such as `I` in
`dat$x` above. A modified version of our function which identifies and responds
to list-column inputs might look like this:

```{r list-col-input}
code <- c ("#' my_function",
           "#'",
           "#' @param x An input",
           "#' @return Something else",
           "#' @examples",
           "#' y <- my_function (x = 1)",
           "#' y <- my_function (x = 1:2)",
           "#' @export",
           "my_function <- function (x) {",
           "  if (methods::is (unclass (x), \"list\")) {",
           "    x <- unlist (x)",
           "  }",
           "  if (is.numeric (x)) {",
           "    x <- as.integer (x)",
           "  }",
           "  return (x + 1L)",
           "}")
writeLines (code, file.path (path, "R", "myfn.R"))
roxygen2::roxygenise (path)
```

That change once again leads to clean `autotest` results:

```{r autotest-TRUE4}
autotest_package (path, test = TRUE)
```

Of course simply attempting to `unlist` a complex list-column may be dangerous,
and it may be preferable to issue some kind of message or warning, or even
either simply remove any list-columns entirely or generate an error. Replacing
the above, potentially dangerous, line, `x <- unlist (x)` with a simple
`stop("list-columns are not allowed")` will also produce clean `autotest`
results.



## Return results and documentation

Functions which return complicated results, such as objects with specific
classes, need to document those class types, and `autotest` compares return
objects with documentation to ensure that this is done. The following code
constructs a new function to demonstrate some of the ways `autotest` inspects
return objects, demonstrating a vector input (`length(x) > 1`) in the example
to avoid messages regarding length checks an integer ranges.

```{r return-val}
code <- c ("#' my_function3",
           "#'",
           "#' @param x An input",
           "#' @examples",
           "#' y <- my_function3 (x = 1:2)",
           "#' @export",
           "my_function3 <- function (x) {",
           "  return (datasets::iris)",
           "}")
writeLines (code, file.path (path, "R", "myfn3.R"))
roxygen2::roxygenise (path) # need to update docs with seed param
x5 <- autotest_package (path, test = TRUE)
DT::datatable (x5, options = list (dom = "t"))
```

Several new diagnostic messages are then issued regarding the description of
the returned value. Let's insert a description to see the effect.

```{r return-val-2}
code <- c (code [1:3],
           "#' @return The iris data set as dataframe",
           code [4:length (code)])
writeLines (code, file.path (path, "R", "myfn3.R"))
roxygen2::roxygenise (path) # need to update docs with seed param
x6 <- autotest_package (path, test = TRUE)
DT::datatable (x6, options = list (dom = "t"))
```

That result still contains a couple of diagnostic messages, but it is now
pretty clear what we need to do, which is to be precise with our specification
of the class of return object. The following then suffices to once again
generate clean `autotest` results.

```{r iris-update}
code [4] <- "#' @return The iris data set as data.frame"
writeLines (code, file.path (path, "R", "myfn3.R"))
roxygen2::roxygenise (path) # need to update docs with seed param
autotest_package (path, test = TRUE)
```

### Documentation of input parameters

Similar checks are performed on the documentation of input parameters, as
demonstrated by the following modified version of the preceding function.

```{r input-checks}
code <- c ("#' my_function3",
           "#'",
           "#' @param x An input",
           "#' @return The iris data set as data.frame",
           "#' @examples",
           "#' y <- my_function3 (x = datasets::iris)",
           "#' @export",
           "my_function3 <- function (x) {",
           "  return (x)",
           "}")
writeLines (code, file.path (path, "R", "myfn3.R"))
roxygen2::roxygenise (path) # need to update docs with seed param
x7 <- autotest_package (path, test = TRUE)
DT::datatable (x7, options = list (dom = "t"))
```

This warning again indicates precisely how it can be rectified, for example by
replacing the third line with

```{r input-fix, eval = FALSE}
code [3] <- "#' @param x An input which can be a data.frame"
```



## General Procedure

The demonstrations above hopefully suffice to indicate the general procedure
which `autotest` attempts to make as simple as possible. This procedure
consists of the following single point:

- From the moment you develop your first function, and every single time you
  modify your code, do whatever steps are necessary to ensure
  `autotest_package()` returns `NULL`.

This vignette has only demonstrated a few of the tests included in the package,
but as long as you use `autotest` throughout the entire process of package
development, any additional diagnostic messages should include sufficient
information for you to be able to restructure your code to avoid them.
---
title: "Control of tests"
author: 
  - "Mark Padgham"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        number_sections: false
        theme: flatly
vignette: >
  %\VignetteIndexEntry{Control of tests}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = TRUE,
  message = TRUE,
  width = 120,
  comment = "#>",
  fig.retina = 2,
  fig.path = "README-"
)
```

The [first
vignette](https://docs.ropensci.org/autotest/articles/autotest.html)
demonstrates the process of applying `autotest` at all stages of package
development. This vignette provides additional information for those applying
`autotest` to already developed packages, in particular through describing how
tests can be selectively applied to a package. By default the
[`autotest_package()`
function](https://docs.ropensci.org/autotest/reference/autotest_package.html)
tests an entire package, but testing can also be restricted to specified
functions only. This vignette will demonstrate application to a few functions
from the [`stats`
package](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html),
starting by loading the package.

```{r pkg-load-fakey, eval = FALSE, echo = TRUE}
library (autotest)
```
```{r pkg-lod, eval = TRUE, echo = FALSE, message = FALSE}
devtools::load_all (".", export_all = FALSE)
```


## 1. `.Rd` files, example code, and the autotest workflow

To understand what `autotest` does, it is first necessary to understand a bit
about the structure of documentation files for R package, which are contained
in files called `".Rd"` files. Tests are constructed by parsing individual
`.Rd` documentation files to extract the example code, identifying parameters
passed to the specified functions, and mutating those parameters.

The general procedure can be illustrated by examining a specific function, for
which we now choose the [`cor`
function](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html),
because of its relative simplicity. The following lines extract the
documentation for the `cor` function, a `.html` version of which can be seen by
clicking on the link above. Note that that web page reveals the name of the
`.Rd` file to be "cor" (in the upper left corner), meaning that the name of the
`.Rd` file is `"cor.Rd"`. The following lines extract that content, first by
loading the entire `.Rd` database for the [`stats`
package](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html).

```{r cor-rd}
rd <- tools::Rd_db (package = "stats")
cor_rd <- rd [[grep ("^cor\\.Rd", names (rd))]]
```

The database itself is a list, with each entry holding the contents of one
`.Rd` file in an object of class `r class(cor_rd)`, which is essentially
a large nested
list of components corresponding to the various `.Rd` tags such as
`\arguments`, `\details`, and `\value`. An internal function from the [`tools`
package](https://stat.ethz.ch/R-manual/R-devel/library/tools/html/00Index.html)
can be used to extract individual components (using the `:::` notation to
access internal functions). For example, a single `.Rd` file often describes
the functionality of several functions, each of which is identified by
specifying the function name as an `"alias"`. The aliases for the `"cor.Rd"`
file are:

```{r cor-aliases}
tools:::.Rd_get_metadata (cor_rd, "alias")
```

This one file thus contains documentation for those four functions. Example
code can be extracted with a dedicated function from the [`tools`
package](https://stat.ethz.ch/R-manual/R-devel/library/tools/html/Rd2HTML.html):

```{r Rd2ex-dummy, echo = TRUE, eval = FALSE}
tools::Rd2ex (cor_rd)
```
```{r Rd2ex, echo = FALSE, eval = TRUE}
tools::Rd2ex (cor_rd)
```

This is the entire content of the `\examples` portion of `"cor.Rd"`, as can be
confirmed by comparing with the [online
version](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html).

## 2. Internal structure of the `autotest` workflow

For each `.Rd` file in a package, `autotest` tests the code given in the
example section according to the following general steps:

1. Extract example lines from the `.Rd` file, as demonstrated above;
2. Identify all function `aliases` described by that file;
3. Identify all points at which those functions are called;
4. Identify all objects passed to those values, including values, classes,
   attributes, and other properties.
5. Identify any other parameters not explicitly passed in example code, but
   defined via default value;
6. Mutate the values of all parameters according to the kinds of test described
   in
   [`autotest_types()`](https://docs.ropensci.org/autotest/reference/autotest_types.html).

Calling `autotest_package(..., test = FALSE)` implements the first 5 of those
6 steps, and returns data on all possible mutations of each parameter, while
setting `test = TRUE` actually passes the mutated parameters to the specified
functions, and returns reports on any unexpected behaviour. 


## 3. `autotest`-ing the `stats::cov` function

The preceding sections describe how `autotest` actually works, while the
present section demonstrates how the package is typically used in practice. As
demonstrated in the
[`README`](https://docs.ropensci.org/autotest/), information on
all tests implemented within the package can be obtained by calling the
[`autotest_types()`
function](https://docs.ropensci.org/autotest/reference/autotest_types.html).
The main function for testing package is
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html).
The nominated package can be either an installed package, or the full path to
a local directory containing a package's source code. By default all `.Rd`
files of a package are tested, with restriction to specified functions possible
either by nominating functions to exclude from testing (via the `exclude`
parameter), or functions to include (via the `functions` parameter). The
`functions` parameter is intended to enable testing only of specified
functions, while the `exclude` parameter is intended to enable testing of all
functions except those specified with this parameter. Specifying values for
both of these parameters is not generally recommended.

### 3.1 Listing tests without conducting them

The following demonstrates the results of `autotest`-ing the [`cor`
function](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html) of
the `stats` package, noting that the default call uses `test = FALSE`, and so
returns details of all tests without actually implementing them (and for this
reason we name the object `xf` for "false"):

```{r at-cor-notest-dummy, echo = TRUE, eval = FALSE}
xf <- autotest_package (package = "stats",
                        functions = "cor")
```
```{r at-cor-notest, fig.show = "hide", collapse = TRUE, echo = FALSE}
xf <- autotest_package (package = "stats",
                        functions = "cor")
```

That reveals that the example code from the `.Rd` file has four main points at
which the aliases defined in that file are tested: once each for `var` and
`cor`, and twice for `cov`. The result looks like this:

```{r print-x-dummy, echo = TRUE, eval = FALSE}
print (xf)
```
```{r print-x, echo = FALSE, collapse = TRUE}
print (xf)
```

The object returned from `autotest_package()` is a simple
[`tibble`](https://tibble.tidyverse.org), with each row detailing one test
which would be applied to the each of the listed functions and parameters.
Because no tests were conducted, all tests will generally have a `type` of
`"dummy"`. In this case, however, we see the following:

```{r test-types}
table (xf$type)
```

In addition to the `r length (which (xf$type == "dummy"))` dummy tests, the
function also returns `r length (which (xf$type == "warning"))` warnings, the
corresponding rows of which are:

```{r test-warnings}
xf [xf$type != "dummy", c ("fn_name", "parameter", "operation", "content")]
```

Although the `auotest` package is primarily intended to apply mutation tests to
all parameters of all functions of a package, doing so requires identifying
parameter types and classes through parsing example code. Any parameters of
a function which are neither demonstrated within example code, nor given
default values can not be tested, because it is not possible to determine their
expected types. The above result reveals that neither the `use` parameter of
the `var` function, nor the `y` parameter of `cov`, are demonstrated in example
code, triggering a warning that these parameter are unable to be tested.

### 3.2 Conducting tests

The `r length (which (xf$type == "dummy"))` tests listed above with `type ==
"dummy"` can then be applied to all nominated functions and parameters by
calling the same function with `test = TRUE`. Doing so yields the following
results (as an object names `xt` for "true"):

```{r at-cor-test-dummy, echo = TRUE, eval = FALSE}
xt <- autotest_package (package = "stats",
                        functions = "cor",
                        test = TRUE)
```
```{r at-cor-test, fig.show = "hide", collapse = TRUE, echo = FALSE}
xt <- autotest_package (package = "stats",
                        functions = "cor",
                        test = TRUE)
```
```{r print-xt-dummy, echo = TRUE, eval = FALSE}
print (xt)
```
```{r print-xt, echo = FALSE, collapse = TRUE}
print (xt)
```

And the `r length (which (xf$type == "dummy"))` tests yielded `r nrow (xt)`
unexpected responses. The best way to understand these results is to examine
the object in detail, typically through `edit(xt)`, or equivalently in RStudio,
clicking on the listed object. The different types of tests which produced
unexpected responses were:

```{r xt-operations}
table (xt$operation)
```

Two of those reflect the previous results regarding parameters unable to be
tested, while the remainder come from only two types of tests. Information on
the precise results is contained in the `content` column, although in this case
it is fairly straightforward to see that the operation "upper case character
parameter" arises because the `use` and `method` parameters of the `cor` and
`cov` functions are case-dependent, and are only accepted in lower case form.
The other operation is the conversion of vectors to list-column format, as
described in the [first
vignette](https://docs.ropensci.org/autotest/articles/autotest.html).


### 3.4 Controlling which tests are conducted

The `test` parameter of the [`autotest_package()`
function](https://docs.ropensci.org/autotest/reference/autotest_package.html)
can be used to control whether all tests are conducted or not. Finer-level
control over tests can be achieved by specifying the `test_data` parameter.
This parameter must be an object of class `autotest_package`, as returned by
either the
[`autotest_types()`](https://docs.ropensci.org/autotest/reference/autotest_types.html)]
or
[`autotest_package()](https://docs.ropensci.org/autotest/reference/autotest_package.html)
functions. The former of these is the function which specifies all unique
tests, and so returns a relatively small `tibble` of `r nrow(autotest_types())`
rows. The following lines demonstrate how to switch off the list-column test
for all functions and parameters:

```{r test_data1-dummy, echo = TRUE, eval = FALSE}
types <- autotest_types()
types$test [grep ("list_col", types$test_name)] <- FALSE
xt2 <- autotest_package (package = "stats",
                         functions = "cor",
                         test = TRUE,
                         test_data = types)
```
```{r test_data1, echo = FALSE, eval = TRUE}
types <- autotest_types()
types$test [grep ("list_col", types$test_name)] <- FALSE
xt2 <- autotest_package (package = "stats",
                         functions = "cor",
                         test = TRUE,
                         test_data = types)
```
```{r print-xt2-dummy, echo = TRUE, eval = FALSE}
print (xt2)
```
```{r print-xt2, echo = FALSE, eval = TRUE}
print (xt2)
```


The result now has four rows with `test == FALSE`, and `type == "no_test"`,
indicating that these tests were not actually conducted. That also makes
apparent the role of these `test` flags. When initially calling 
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html)
with default `test = FALSE`, the result contains a `test` column in which all
values are `TRUE`. Although potentially perplexing at first, this value must be
understood in relation to the `type` column. A `type` of `"dummy"` indicates
that a test has not been conducted, in which case `test = TRUE` is a control
flag used to determine what would be conducted if these data were submitted as
the `test_data` parameter. For all `type` values other than `"dummy"`, the
`test` column specifies whether or not each test was actually conducted.

The preceding example showed how the results of 
[`autotest_types()`](https://docs.ropensci.org/autotest/reference/autotest_types.html)
can be used to control which tests are implemented for an entire package.
Finer-scale control can be achieved by modifying individual rows of the full
table returned by
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html).
The following code demonstrates by showing how list-column tests can be
switched off only for
particular functions, starting again with the `xf` data of dummy tests
generated above.

```{r test_data2-fakey, echo = TRUE, eval = FALSE}
xf <- autotest_package (package = "stats",
                        functions = "cor")
xf$test [grepl ("list_col", xf$test_name) & xf$fn_name == "var"] <- FALSE
xt3 <- autotest_package (package = "stats",
                         functions = "cor",
                         test = TRUE,
                         test_data = xf)
```
```{r test_data2, echo = FALSE, eval = TRUE}
xf <- autotest_package (package = "stats",
                        functions = "cor")
xf$test [grepl ("list_col", xf$test_name) & xf$fn_name == "var"] <- FALSE
xt3 <- autotest_package (package = "stats",
                         functions = "cor",
                         test = TRUE,
                         test_data = xf)
```
```{r print-xt3-dummy, echo = TRUE, eval = FALSE}
print (xt3)
```
```{r print-xt3, echo = FALSE, eval = TRUE}
print (xt3)
```

These procedures illustrate the three successively finer levels of control over
tests, by switching them off for:

1. Entire packages;
2. Specified functions only; or
3. Specific parameters of particular functions only.

## 4. `autotest`-ing your package

`autotest` can be very easily incorporated in your package's `tests/` directory
via to simple [`testthat`](https://testthat.r-lib.org) expectations:

- `expect_autotest_no_testdata`, which will expect `autotest_package` to work
  on your package with default values including no additional `test_data`
  specifying tests which are not to be run; or
- `expect_autotest_testdata`, to be used when specific tests are switched off.

Using these requires adding `autotest` to the `Suggests` list of your package's
`DESCRIPTION` file, along with `testthat`. Note that the use of testing
frameworks other than [`testthat`](https://testthat.r-lib.org) is possible
through writing custom expectations for the output of
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html),
but that is not considered here.

To use these expectations, you must first decide which, if any, tests you judge
to be not applicable to your package, and switch them off following the
procedure described above (that is, at the package level through modifying the
`test` flag of the object returned from
[`autotest_types()`](https://docs.ropensci.org/autotest/reference/autotest_types.html),
or at finer function- or parameter-levels by modifying equivalent values in the
object returned from [`autotest_package(..., test
= FALSE)`](https://docs.ropensci.org/autotest/reference/autotest_package.html).
These objects must then be passed as the `test_data` parameter to
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html).
If you consider all tests to be applicable, then
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html)
can be called without specifying this parameter.

If you switch tests off via a `test_data` parameter, then the
`expect_autotest` expectation requires you to append an additional column to
the `test_data` object called `"note"` (case-insensitive), and include a note
for each row which has `test = FALSE` explaining why those tests have been
switched off. Lines in your test directory should look something like this:

```{r testthat-demo-testdata, collapse = TRUE, echo = TRUE, eval = FALSE}
library (testthat) # as called in your test suite
# For example, to switch off vector-to-list-column tests:
test_data <- autotest_types (notest = "vector_to_list_col")
test_data$note <- ""
test_data$note [test_data$test == "vector_to_list_col"] <-
    "These tests are not applicable because ..."
expect_success (expect_autotest_testdata (test_data))
```
This procedure of requiring an additional `"note"` column ensures that your own
test suite will explicitly include all explanations of why you deem particular
tests not to be applicable to your package. 

In contrast, the following expectation should be used when `autotest_package()`
passes with all tests are implemented, in which case no parameters need be
passed to the expectation, and tests will confirm that no warnings or errors
are generated.

```{r testthat-demo-notestdata, collapse = TRUE, echo = TRUE, eval = FALSE}
expect_success (expect_autotest_no_testdata ())
```

### 4.2 Finer control over testing expectations

The two expectations shown above call the 
[`autotest_package()`](https://docs.ropensci.org/autotest/reference/autotest_package.html)
function internally, and assert that the results follow the expected pattern.
There are also three additional [`testthat`](https://testthat.r-lib.org)
expectations which can be applied to pre-generated `autotest` objects, to allow for finer control over testing expectations. These are:

- `expect_autotest_no_err` to expect no errors in results from `autotest_package()`;
- `expect_autotest_no_warn` to expect no warnings; and
- `expect_autotest_notes` to expect tests which have been switched off to have
  an additional `"note"` column explaining why.

These tests are demonstrated in one of the testing files used in this package,
which the following lines recreate to demonstrate the general process. The
first two expectations are that an object be free from both warnings and
errors. The tests implemented here are applied to the 
[`stats::cov()`
function](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html),
which actually triggers warnings because two parameters do not have their usage
demonstrated in the example code. The tests therefore `expect_failure()`, when
they generally should `expect_success()` throughout.

```{r testthat-no-err-warn-fakey, collapse = TRUE, echo = TRUE, eval = FALSE}
library (testthat) # as called in your test suite
# For example, to switch off vector-to-list-column tests:
test_data <- autotest_types (notest = "vector_to_list_col")
x <- autotest_package (package = "stats",
                       functions = "cov",
                       test = TRUE,
                       test_data = test_data)
       
expect_success (expect_autotest_no_err (x))
expect_failure (expect_autotest_no_warn (x)) # should expect_success!!
```
```{r testthat-no-err-warn, collapse = TRUE, echo = FALSE, eval = TRUE}
library (testthat) # as called in your test suite
# switch off vector-to-list-column tests:
test_data <- autotest_types (notest = "vector_to_list_col")
x <- autotest_package (package = "stats",
                       functions = "cov",
                       test = TRUE,
                       test_data = test_data)
       
expect_success (expect_autotest_no_err (x))
expect_failure (expect_autotest_no_warn (x)) # should expect_success!!
```

The test files then affirms that simply passing the object, `x`, which has
tests flagged as `type == "no_test"`,  yet without explaining why in an
additional `"note"` column, should cause `expect_autotest()` to fail. The
following line, removing the logical `testthat` expectation, demonstrates:

```{r testthat-demo-fail-fakey, echo = TRUE, eval = FALSE}
expect_autotest_notes (x)
```
```{r testthat-demo-fail, error = TRUE, echo = FALSE, eval = TRUE}
expect_autotest_notes (x)
```

As demonstrated above, these `expect_autotest_...` calls should always be
wrapped in a direct [`testhat`](https://testthat.r-lib.org) expectation of
[`expect_success()`](https://testthat.r-lib.org/reference/expect_success.html).
To achieve success in that case, we need to append an additional `"note"`
column containing explanations of why each test has been switched off:

```{r testthat-demo-success}
x$note <- ""
x [grep ("vector_to_list", x$test_name), "note"] <-
  "these tests have been switched off because ..."

expect_success (expect_autotest_notes (x))
```

In general, using `autotest` in a package's test suite should be as simple as
adding `autotest` to `Suggests`, and wrapping either
`expect_autotest_no_testdata` or `expect_autotest_testdata` in an
`expect_success` call.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testthat-expections.R
\name{expect_autotest_no_err}
\alias{expect_autotest_no_err}
\title{expect_autotest_no_err}
\usage{
expect_autotest_no_err(object)
}
\arguments{
\item{object}{An \code{autotest} object to be tested}
}
\value{
(invisibly) The same object
}
\description{
Expect \code{autotest_package()} to be clear of errors
}
\seealso{
Other expectations: 
\code{\link{expect_autotest_no_testdata}()},
\code{\link{expect_autotest_no_warn}()},
\code{\link{expect_autotest_notes}()},
\code{\link{expect_autotest_testdata}()}
}
\concept{expectations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autotest-class.R
\name{autotest_obj}
\alias{autotest_obj}
\title{autotest_obj class definition}
\usage{
autotest_obj(
  package = NA_character_,
  package_loc = NULL,
  test_name = NA_character_,
  fn_name = NA_character_,
  parameters = list(),
  parameter_types = NA_character_,
  class = NULL,
  classes = NULL,
  env = new.env(),
  test = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{package}{Name of package for which object is to be constructed.}

\item{package_loc}{Location of package on local system (for source packages
only)}

\item{test_name}{Name of test (use \link{autotest_types} to get all test
names).}

\item{fn_name}{Name of function to be tested.}

\item{parameters}{Names of all parameters for that function.}

\item{parameter_types}{Types of input parameters.}

\item{class}{Class of an individual parameter.}

\item{classes}{Classes of all parameters.}

\item{env}{Environment in which tests are to be run.}

\item{test}{If \code{FALSE}, return only descriptions of tests which would be run
with \code{test = TRUE}, without actually running them.}

\item{quiet}{If \code{FALSE}, issue progress and other messages during testing of
object.}
}
\description{
This function exists only to provide the class definitions for test objects,
and is not intended to be called directly.
}
\concept{class}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testthat-expections.R
\name{expect_autotest_no_warn}
\alias{expect_autotest_no_warn}
\title{expect_autotest_no_warn}
\usage{
expect_autotest_no_warn(object)
}
\arguments{
\item{object}{An \code{autotest} object to be tested}
}
\value{
(invisibly) The same object
}
\description{
Expect \code{autotest_package()} to be clear of warnings
}
\seealso{
Other expectations: 
\code{\link{expect_autotest_no_err}()},
\code{\link{expect_autotest_no_testdata}()},
\code{\link{expect_autotest_notes}()},
\code{\link{expect_autotest_testdata}()}
}
\concept{expectations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse-yaml.R
\name{at_yaml_template}
\alias{at_yaml_template}
\title{at_yaml_template}
\usage{
at_yaml_template(loc = tempdir())
}
\arguments{
\item{loc}{Location to generate template file. Append with filename and
'.yaml' suffix to overwrite default name of 'autotest.yaml', otherwise this
parameter will be used to specify directory only.}
}
\description{
Generate a 'yaml' template for an 'autotest'.
}
\seealso{
Other yaml: 
\code{\link{autotest_yaml}()},
\code{\link{examples_to_yaml}()}
}
\concept{yaml}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testthat-expections.R
\name{expect_autotest_notes}
\alias{expect_autotest_notes}
\title{expect_autotest_notes}
\usage{
expect_autotest_notes(object)
}
\arguments{
\item{object}{An \code{autotest} object to be tested}
}
\description{
Expect \code{test_data} param of \code{autotest_package} to have additional \code{note}
column explaining why tests have been switched off.
}
\seealso{
Other expectations: 
\code{\link{expect_autotest_no_err}()},
\code{\link{expect_autotest_no_testdata}()},
\code{\link{expect_autotest_no_warn}()},
\code{\link{expect_autotest_testdata}()}
}
\concept{expectations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autotest-package.R
\docType{package}
\name{autotest-package}
\alias{autotest-package}
\alias{_PACKAGE}
\title{autotest: Automatic Package Testing}
\description{
Automatic testing of R packages via a simple YAML schema.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/autotest/}
  \item \url{https://github.com/ropensci-review-tools/autotest}
  \item Report bugs at \url{https://github.com/ropensci-review-tools/autotest/issues}
}

}
\author{
\strong{Maintainer}: Mark Padgham \email{mark.padgham@email.com}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autotest-functions.R
\name{autotest_types}
\alias{autotest_types}
\title{autotest_types}
\usage{
autotest_types(notest = NULL)
}
\arguments{
\item{notest}{Character string of names of tests which should be switched off
by setting the \code{test} column to \code{FALSE}. Run this function first without this
parameter to get all names, then re-run with this parameter switch specified
tests off.}
}
\value{
An \code{autotest} object with each row listing one unique type of test
which can be applied to every parameter (of the appropriate class) of each
function.
}
\description{
List all types of 'autotests' currently implemented.
}
\seealso{
Other main_functions: 
\code{\link{autotest_package}()}
}
\concept{main_functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/examples-to-yaml.R
\name{examples_to_yaml}
\alias{examples_to_yaml}
\title{examples_to_yaml}
\usage{
examples_to_yaml(
  package = NULL,
  functions = NULL,
  exclude = NULL,
  quiet = FALSE
)
}
\arguments{
\item{package}{Name of package, as either
\enumerate{
\item Path to local package source
\item Name of installed package
\item Full path to location of installed package if not on \link{.libPaths}, or
\item Default which presumes current directory is within package to be
tested.
}}

\item{functions}{If specified, names of functions from which examples are to
be obtained.}

\item{exclude}{Names of functions to exclude from 'yaml' template}

\item{quiet}{If 'FALSE', provide printed output on screen.}
}
\description{
Convert examples for a specified package, optionally restricted to one or
more specified functions, to a list of 'autotest' 'yaml' objects to use to
automatically test package.
}
\seealso{
Other yaml: 
\code{\link{at_yaml_template}()},
\code{\link{autotest_yaml}()}
}
\concept{yaml}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autotest-functions.R
\name{autotest_package}
\alias{autotest_package}
\title{autotest_package}
\usage{
autotest_package(
  package = ".",
  functions = NULL,
  exclude = NULL,
  test = FALSE,
  test_data = NULL,
  quiet = FALSE
)
}
\arguments{
\item{package}{Name of package, as either
\enumerate{
\item Path to local package source
\item Name of installed package
\item Full path to location of installed package if not on \link{.libPaths}, or
\item Default which presumes current directory is within package to be
tested.
}}

\item{functions}{Optional character vector containing names of functions of
nominated package to be included in 'autotesting'.}

\item{exclude}{Optional character vector containing names of any functions of
nominated package to be excluded from 'autotesting'.}

\item{test}{If \code{FALSE}, return only descriptions of tests which would be run
with \code{test = TRUE}, without actually running them.}

\item{test_data}{Result returned from calling either \link{autotest_types} or
\link{autotest_package} with \code{test = FALSE} that contains a list of all tests
which would be conducted. These tests have an additional flag, \code{test}, which
defaults to \code{TRUE}. Setting any tests to \code{FALSE} will avoid running them when
\code{test = TRUE}.}

\item{quiet}{If 'FALSE', provide printed output on screen.}
}
\value{
An \code{autotest_package} object which is derived from a \pkg{tibble}
\code{tbl_df} object. This has one row for each test, and the following nine
columns:
\enumerate{
\item \code{type} The type of result, either "dummy" for \code{test = FALSE}, or one
of "error", "warning", "diagnostic", or "message".
\item \code{test_name} Name of each test
\item \code{fn_name} Name of function being tested
\item \code{parameter} Name of parameter being tested
\item \code{parameter_type} Expected type of parameter as identified by
\code{autotest}.
\item \code{operation} Description of the test
\item \code{content} For \code{test = FALSE}, the expected behaviour of the test; for
\code{test = TRUE}, the observed discrepancy with that expected behaviour
\item \code{test} If \code{FALSE} (default), list all tests without implementing them,
otherwise implement all tests.
\item \verb{yaml_hash' A unique hash which may be be used to extract the }yaml`
specification of each test.
}
Some columns may contain NA values, as explained in the Note.
}
\description{
Automatically test an entire package by converting examples to \code{yaml} format
and submitting each to the \link{autotest_yaml} function.
}
\note{
Some columns may contain NA values, including:
\itemize{
\item \code{parameer} and \code{parameter_type}, for tests applied to entire
functions, such as tests of return values.
\item \code{test_name} for warnings or errors generated through "normal"
function calls generated directly from example code, in which case \code{type}
will be "warning" or "error", and \code{content} will contain the content of
the corresponding message.
}
}
\seealso{
Other main_functions: 
\code{\link{autotest_types}()}
}
\concept{main_functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testthat-expections.R
\name{expect_autotest_testdata}
\alias{expect_autotest_testdata}
\title{expect_autotest_testdata}
\usage{
expect_autotest_testdata(object)
}
\arguments{
\item{object}{An \code{autotest_package} object with a \code{test} column flagging
tests which are not to be run on the local package.}
}
\value{
(invisibly) The autotest object
}
\description{
Expect \code{autotest_package()} to be clear of errors with some tests switched
off, and to have note column explaining why those tests are not run.
}
\seealso{
Other expectations: 
\code{\link{expect_autotest_no_err}()},
\code{\link{expect_autotest_no_testdata}()},
\code{\link{expect_autotest_no_warn}()},
\code{\link{expect_autotest_notes}()}
}
\concept{expectations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autotest-functions.R
\name{autotest_yaml}
\alias{autotest_yaml}
\title{autotest_yaml}
\usage{
autotest_yaml(
  yaml = NULL,
  filename = NULL,
  test = TRUE,
  test_data = NULL,
  quiet = FALSE
)
}
\arguments{
\item{yaml}{A 'yaml' template as a character vector, either hand-coded or
potentially loaded via \link{readLines} function or similar. Should generally
be left at default of 'NULL', with template specified by 'filename'
parameter.}

\item{filename}{Name (potentially including path) of file containing 'yaml'
template. See \link{at_yaml_template} for details of template. Default uses
template generated by that function, and held in local './tests' directory.}

\item{test}{If \code{FALSE}, return only descriptions of tests which would be run
with \code{test = TRUE}, without actually running them.}

\item{test_data}{Result returned from calling either \link{autotest_types} or
\link{autotest_package} with \code{test = FALSE} that contains a list of all tests
which would be conducted. These tests have an additional flag, \code{test}, which
defaults to \code{TRUE}. Setting any tests to \code{FALSE} will avoid running them when
\code{test = TRUE}.}

\item{quiet}{If 'FALSE', provide printed output on screen.}
}
\value{
An \code{autotest_pkg} object, derived from a \pkg{tibble}, detailing
instances of unexpected behaviour for every parameter of every function.
}
\description{
Automatically test inputs to functions specified in a 'yaml' template.
}
\examples{
\dontrun{
yaml_list <- examples_to_yaml (package = "stats", functions = "reshape")
res <- autotest_yaml (yaml = yaml_list)
}
}
\seealso{
Other yaml: 
\code{\link{at_yaml_template}()},
\code{\link{examples_to_yaml}()}
}
\concept{yaml}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testthat-expections.R
\name{expect_autotest_no_testdata}
\alias{expect_autotest_no_testdata}
\title{expect_autotest_no_testdata}
\usage{
expect_autotest_no_testdata(object = NULL)
}
\arguments{
\item{object}{Not used here, but required for \code{testthat} expectations}
}
\value{
(invisibly) The autotest object
}
\description{
Expect \code{autotest_package()} to be clear of errors with no tests switched off
}
\seealso{
Other expectations: 
\code{\link{expect_autotest_no_err}()},
\code{\link{expect_autotest_no_warn}()},
\code{\link{expect_autotest_notes}()},
\code{\link{expect_autotest_testdata}()}
}
\concept{expectations}
