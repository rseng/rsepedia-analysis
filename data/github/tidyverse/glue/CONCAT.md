## Current state

Glue is currently stable, it is used and installed quite extensively and has relatively few issues reported.
Most issues reported tend to center around `glue_sql()`.

## Known outstanding issues

I think we need a way to optionally disable interpretation of quotes, which would fix both https://github.com/r-lib/cli/issues/370, and the long standing https://github.com/tidyverse/glue/issues/158

Possibly the loading and closing blank line behaviors could be tweaked, so that the issues brought up in https://github.com/tidyverse/glue/issues/234 are more natural.
It may be difficult to do this while maintaining backwards compatibility.

## Future directions

Should we include more built-in support for custom numeric formats, as suggested by https://github.com/tidyverse/glue/issues/87 and https://github.com/tidyverse/glue/issues/86?
For now I have avoided this, but it is a common pain point for glue

Potentially we could revisit the use of the custom glue class, which I gather also causes people grief in practice?

<!-- README.md is generated from README.Rmd. Please edit that file -->

# glue <a href='https://glue.tidyverse.org'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/glue)](https://cran.r-project.org/package=glue)
[![R-CMD-check](https://github.com/tidyverse/glue/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tidyverse/glue/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/tidyverse/glue/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/tidyverse/glue/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

## Overview

Glue offers interpreted string literals that are small, fast, and
dependency-free. Glue does this by embedding R expressions in curly
braces which are then evaluated and inserted into the argument string.

## Installation

<div class=".pkgdown-release">

``` r
# Install released version from CRAN
install.packages("glue")
```

</div>

<div class=".pkgdown-devel">

``` r
# Install development version from GitHub
devtools::install_github("tidyverse/glue")
```

</div>

## Usage

##### Variables can be passed directly into strings.

``` r
library(glue)
name <- "Fred"
glue('My name is {name}.')
#> My name is Fred.
```

Note that `glue::glue()` is also made available via
`stringr::str_glue()`. So if you’ve already attached stringr (or perhaps
the whole tidyverse), you can access `glue()` like so:

``` r
library(stringr) # or library(tidyverse)

stringr_fcn <- "`stringr::str_glue()`"
glue_fcn    <- "`glue::glue()`"

str_glue('{stringr_fcn} is essentially an alias for {glue_fcn}.')
#> `stringr::str_glue()` is essentially an alias for `glue::glue()`.
```

##### Long strings are broken by line and concatenated together.

``` r
library(glue)

name <- "Fred"
age <- 50
anniversary <- as.Date("1991-10-12")
glue('My name is {name},',
  ' my age next year is {age + 1},',
  ' my anniversary is {format(anniversary, "%A, %B %d, %Y")}.')
#> My name is Fred, my age next year is 51, my anniversary is Saturday, October 12, 1991.
```

##### Named arguments are used to assign temporary variables.

``` r
glue('My name is {name},',
  ' my age next year is {age + 1},',
  ' my anniversary is {format(anniversary, "%A, %B %d, %Y")}.',
  name = "Joe",
  age = 40,
  anniversary = as.Date("2001-10-12"))
#> My name is Joe, my age next year is 41, my anniversary is Friday, October 12, 2001.
```

##### `glue_data()` is useful with [magrittr](https://cran.r-project.org/package=magrittr) pipes.

``` r
`%>%` <- magrittr::`%>%`
head(mtcars) %>% glue_data("{rownames(.)} has {hp} hp")
#> Mazda RX4 has 110 hp
#> Mazda RX4 Wag has 110 hp
#> Datsun 710 has 93 hp
#> Hornet 4 Drive has 110 hp
#> Hornet Sportabout has 175 hp
#> Valiant has 105 hp
```

##### Or within dplyr pipelines

``` r
library(dplyr)
head(iris) %>%
  mutate(description = glue("This {Species} has a petal length of {Petal.Length}"))
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 1          5.1         3.5          1.4         0.2  setosa
#> 2          4.9         3.0          1.4         0.2  setosa
#> 3          4.7         3.2          1.3         0.2  setosa
#> 4          4.6         3.1          1.5         0.2  setosa
#> 5          5.0         3.6          1.4         0.2  setosa
#> 6          5.4         3.9          1.7         0.4  setosa
#>                             description
#> 1 This setosa has a petal length of 1.4
#> 2 This setosa has a petal length of 1.4
#> 3 This setosa has a petal length of 1.3
#> 4 This setosa has a petal length of 1.5
#> 5 This setosa has a petal length of 1.4
#> 6 This setosa has a petal length of 1.7
```

##### Leading whitespace and blank lines from the first and last lines are automatically trimmed.

This lets you indent the strings naturally in code.

``` r
glue("
    A formatted string
    Can have multiple lines
      with additional indention preserved
    ")
#> A formatted string
#> Can have multiple lines
#>   with additional indention preserved
```

##### An additional newline can be used if you want a leading or trailing newline.

``` r
glue("

  leading or trailing newlines can be added explicitly

  ")
#> 
#> leading or trailing newlines can be added explicitly
```

##### `\\` at the end of a line continues it without a new line.

``` r
glue("
    A formatted string \\
    can also be on a \\
    single line
    ")
#> A formatted string can also be on a single line
```

##### A literal brace is inserted by using doubled braces.

``` r
name <- "Fred"
glue("My name is {name}, not {{name}}.")
#> My name is Fred, not {name}.
```

##### Alternative delimiters can be specified with `.open` and `.close`.

``` r
one <- "1"
glue("The value of $e^{2\\pi i}$ is $<<one>>$.", .open = "<<", .close = ">>")
#> The value of $e^{2\pi i}$ is $1$.
```

##### All valid R code works in expressions, including braces and escaping.

Backslashes do need to be doubled just like in all R strings.

``` r
  `foo}\`` <- "foo"
glue("{
      {
        '}\\'' # { and } in comments, single quotes
        \"}\\\"\" # or double quotes are ignored
        `foo}\\`` # as are { in backticks
      }
  }")
#> foo
```

##### `glue_sql()` makes constructing SQL statements safe and easy

Use backticks to quote identifiers, normal strings and numbers are
quoted appropriately for your backend.

``` r
library(glue)

con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
colnames(iris) <- gsub("[.]", "_", tolower(colnames(iris)))
DBI::dbWriteTable(con, "iris", iris)
var <- "sepal_width"
tbl <- "iris"
num <- 2
val <- "setosa"
glue_sql("
  SELECT {`var`}
  FROM {`tbl`}
  WHERE {`tbl`}.sepal_length > {num}
    AND {`tbl`}.species = {val}
  ", .con = con)
#> <SQL> SELECT `sepal_width`
#> FROM `iris`
#> WHERE `iris`.sepal_length > 2
#>   AND `iris`.species = 'setosa'

# `glue_sql()` can be used in conjunction with parameterized queries using
# `DBI::dbBind()` to provide protection for SQL Injection attacks
 sql <- glue_sql("
    SELECT {`var`}
    FROM {`tbl`}
    WHERE {`tbl`}.sepal_length > ?
  ", .con = con)
query <- DBI::dbSendQuery(con, sql)
DBI::dbBind(query, list(num))
DBI::dbFetch(query, n = 4)
#>   sepal_width
#> 1         3.5
#> 2         3.0
#> 3         3.2
#> 4         3.1
DBI::dbClearResult(query)

# `glue_sql()` can be used to build up more complex queries with
# interchangeable sub queries. It returns `DBI::SQL()` objects which are
# properly protected from quoting.
sub_query <- glue_sql("
  SELECT *
  FROM {`tbl`}
  ", .con = con)

glue_sql("
  SELECT s.{`var`}
  FROM ({sub_query}) AS s
  ", .con = con)
#> <SQL> SELECT s.`sepal_width`
#> FROM (SELECT *
#> FROM `iris`) AS s

# If you want to input multiple values for use in SQL IN statements put `*`
# at the end of the value and the values will be collapsed and quoted appropriately.
glue_sql("SELECT * FROM {`tbl`} WHERE sepal_length IN ({vals*})",
  vals = 1, .con = con)
#> <SQL> SELECT * FROM `iris` WHERE sepal_length IN (1)

glue_sql("SELECT * FROM {`tbl`} WHERE sepal_length IN ({vals*})",
  vals = 1:5, .con = con)
#> <SQL> SELECT * FROM `iris` WHERE sepal_length IN (1, 2, 3, 4, 5)

glue_sql("SELECT * FROM {`tbl`} WHERE species IN ({vals*})",
  vals = "setosa", .con = con)
#> <SQL> SELECT * FROM `iris` WHERE species IN ('setosa')

glue_sql("SELECT * FROM {`tbl`} WHERE species IN ({vals*})",
  vals = c("setosa", "versicolor"), .con = con)
#> <SQL> SELECT * FROM `iris` WHERE species IN ('setosa', 'versicolor')
```

##### Optionally combine strings with `+`

``` r
x <- 1
y <- 3
glue("x + y") + " = {x + y}"
#> x + y = 4
```

# Other implementations

Some other implementations of string interpolation in R (although not
using identical syntax).

-   [stringr::str_interp](https://stringr.tidyverse.org/reference/str_interp.html)
-   [R.utils::gstring](https://cran.r-project.org/package=R.utils)
-   [rprintf](https://cran.r-project.org/package=rprintf)

String templating is closely related to string interpolation, although
not exactly the same concept. Some packages implementing string
templating in R include.

-   [whisker](https://cran.r-project.org/package=whisker)
-   [brew](https://cran.r-project.org/package=brew)

## Code of Conduct

Please note that the glue project is released with a [Contributor Code
of Conduct](https://glue.tidyverse.org/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.
# glue (development version)

# glue 1.6.1

* glue now registers its custom knitr engines in a way that is more robust to namespace-loading edge cases that can arise during package installation (#254).

# glue 1.6.0

* `glue()`, `glue_data()`, `glue_col()`, and `glue_data_col()` gain a new `.literal` argument, which controls how quotes and the comment character are treated when parsing the expression string (#235). This is mostly useful when using a custom transformer.

* Trailing whitespace-only lines don't interfere with indentation (#247).

# glue 1.5.1

* Jennifer Bryan is now the maintainer.

* The existing custom language engines for knitr, `glue` and `glue_sql`, are documented in a new vignette (#71). *Detail added after release: glue now sets up registration of these engines in `.onLoad()`.*

* `glue_col()` gives special treatment to styling functions from the crayon package, e.g. `glue_col("{blue foo}")` "just works" now, even if crayon is not attached (but is installed) (#241).

* Unterminated backticks trigger the same error as unterminated single or double quotes (#237).

* `glue_sql()` collapses zero-length `DBI::SQL` object into `DBI::SQL("NULL")` (#244 @shrektan).

# glue 1.5.0

## Breaking changes

* Long deprecated function `collapse()` has been removed (#213)

## New functions and arguments

* New `glue_sql_collapse()` function to collapse inputs and return a `DBI::SQL()` object (#103).

* `glue()` gains a new `.comment` argument, to control the comment character (#193).
* `glue()` gains a new `.null` argument, to control the value to replace `NULL` values with (#217, @echasnovski).

## Bugfixes and minor changes

* `sql_quote_transformer()` is now allows whitespace after the trailing `*` (#218).
* `compare_proxy.glue()` method defined so glue objects can be compared to strings in testthat 3e without errors (#212)
* `print.glue()` no longer prints an empty newline for 0 length inputs (#214)
* Unterminated comments in glue expression now throw an error (#227, @gaborcsardi)
* Unterminated quotes in glue expressions now throw an error (#226, @gaborcsardi)


# glue 1.4.2

* `glue_safe()` gives a slightly nicer error message
* The required version of R is now 3.2 (#189)
* `glue_sql()` now collapses `DBI::SQL()` elements correctly (#192 @shrektan)
* The internal `compare()` method gains a `...` argument, for compatibility with testthat 3.0.0

# glue 1.4.1

* Internal changes for compatibility with vctrs 0.3.0 (#187).
* `glue_sql()` now replaces missing values correctly when collapsing values (#185).
* `glue_sql()` now always preserves the type of the column even in the presence of missing values (#130)

# glue 1.4.0

* `.envir = NULL` is now supported and is equivalent to passing `.envir = emptyenv()` (#140)

* New `glue_safe()` and `glue_data_safe()` functions, safer versions of
  `glue()` that do not execute code, only look up values (using `get()`). These
  alternatives are useful for things like shiny applications where you do not
  have control of the input for your glue expressions. (#140)

* Fixed memory access issue and memory leaks found by valgrind.

# glue 1.3.2

* glue now implements vctrs methods. This ensures that vectors of glue
  strings are compatible with tidyverse packages like tidyr
  (r-lib/tidyselect#170, tidyverse/tidyr#773, @lionel-).

* Fix a LTO type mismatch warning (#146)

* `glue_sql()` now quotes lists of values appropriate to their type, rather
  than coercing all values to characters (#153)

* `glue_data()` now implicitly coerces `.x` to a list.

* `glue()` gains the `.trim` argument, like `glue_data()`.

* `single_quote()` `double_quote()` and `backtick()` all return `NA` for `NA`
  inputs (#135).

* Improve `trim()`'s handling of lines containing only indentation (#162, #163, @alandipert)

# glue 1.3.1

## Features

* `glue()` now has a `+` method to combine strings.
* `glue_sql()` now collapses zero-length vector into `DBI::SQL("NULL")` (#134 @shrektan).

## Bugfixes and minor changes

* `glue_sql()` now supports unquoting lists of Id objects.
* `glue_sql()` now quotes characters with NAs appropriately (#115).
* `glue_sql()` now quotes Dates appropriately (#98).
* A potential protection error reported by rchk was fixed.

# glue 1.3.0

## Breaking changes

* The `evaluate()` function has been removed. Changes elsewhere in glue made
  the implementation trivial so it was removed for the sake of clarity.
  Previous uses can be replaced by `eval(parse(text = text), envir)`.

* `collapse()` has been renamed to `glue_collapse()` to avoid namespace
  collisions with `dplyr::collapse()`.

## Features

* `compare.glue()` was added, to make it easier to use glue objects in
  `testthat::expect_equal()` statements.

* `glue_col()` and `glue_data_col()` functions added to display strings with
  color.

## Bugfixes and minor changes

* Glue now throws an informative error message when it cannot interpolate a
  function into a string (#114, @haleyjeppson & @ijlyttle).

* Glue now evaluates unnamed arguments lazily with `delayedAssign()`, so there
  is no performance cost if an argument is not used. (#83, @egnha).

* Fixed a bug where names in the assigned expression of an interpolation
  variable would conflict with the name of the variable itself (#89, @egnha).

* Do not drop the `glue` class when subsetting (#66).

* Fix `glue()` and `collapse()` always return UTF-8 encoded strings (#81, @dpprdan)

# glue 1.2.0

* The implementation has been tweaked to be slightly faster in most cases.

* `glue()` now has a `.transformer` argument, which allows you to use custom
  logic on how to evaluate the code within glue blocks. See
  `vignette("transformers")` for more details and example transformer
  functions.

* `glue()` now returns `NA` if any of the results are `NA` and `.na` is `NULL`.
  Otherwise `NA` values are replaced by the value of `.na`.

* `trim()` to use the trimming logic from glue is now exported.

* `glue_sql()` and `glue_data_sql()` functions added to make constructing SQL
  statements with glue safer and easier.

* `glue()` is now easier to use when used within helper functions such as
  `lapply`.

* Fix when last expression in `glue()` is NULL.

# glue 1.1.1

* Another fix for PROTECT / REPROTECT found by the rchk static analyzer.

# glue 1.1.0

* Fix for PROTECT errors when resizing output strings.

* `glue()` always returns 'UTF-8' strings, converting inputs if in other
encodings if needed.

* `to()` and `to_data()` have been removed.

* `glue()` and `glue_data()` can now take alternative delimiters to `{` and `}`.
This is useful if you are writing to a format that uses a lot of braces, such
as LaTeX. (#23)

* `collapse()` now returns 0 length output if given 0 length input (#28).

# glue 0.0.0.9000

* Fix `glue()` to admit `.` as an embedded expression in a string (#15, @egnha).

* Added a `NEWS.md` file to track changes to the package.
## revdepcheck results

We checked 463 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 1 new problem
 * We failed to check 0 packages

Issues with CRAN packages are summarised below.

### New problems

* iotables

  As of 2022-01-21, I see very similar failure for iotables on CRAN, on
  flavours r-release-windows-ix86+x86_64 and
  r-oldrel-windows-ix86+x86_64.
  I suspect that the iotables package makes live API calls in its vignettes and
  is subject to intermittent failure.
  I don't think this has anything to do with glue.

  checking re-building of vignette outputs ... WARNING

  Error(s) in re-building vignettes:
  --- re-building ‘environmental_impact.Rmd’ using rmarkdown
  ...
  Quitting from lines 105-106 (environmental_impact.Rmd) 
  Error: processing vignette 'environmental_impact.Rmd' failed with diagnostics:
  get_eurostat_raw fails with the id env_ac_ainah_r2
  --- failed re-building ‘environmental_impact.Rmd’
  ...
  SUMMARY: processing the following file failed:
    ‘environmental_impact.Rmd’
  Error: Vignette re-building failed.
# `.literal` treats quotes and `#` as regular characters

    Code
      glue("{'fo`o\"#}", .transformer = function(x, ...) x)
    Error <simpleError>
      Unterminated quote (')

# glue_col() can exploit the `.literal` argument

    Code
      glue_col("Colorless {green idea's} sleep furiously")
    Error <simpleError>
      Unterminated quote (')

---

    Code
      glue_col("Colorless {green idea\"s} sleep furiously")
    Error <simpleError>
      Unterminated quote (")

---

    Code
      glue_col("Colorless {green idea`s} sleep furiously")
    Error <simpleError>
      Unterminated quote (`)

---

    Code
      glue_col("Hey a URL: {blue https://example.com/#section}")
    Error <simpleError>
      Unterminated comment

# Getting help with glue

Thanks for using glue!
Before filing an issue, there are a few places to explore and pieces to put together to make the process as smooth as possible.

## Make a reprex

Start by making a minimal **repr**oducible **ex**ample using the  [reprex](https://reprex.tidyverse.org/) package. 
If you haven't heard of or used reprex before, you're in for a treat! 
Seriously, reprex will make all of your R-question-asking endeavors easier (which is a pretty insane ROI for the five to ten minutes it'll take you to learn what it's all about). 
For additional reprex pointers, check out the [Get help!](https://www.tidyverse.org/help/) section of the tidyverse site.

## Where to ask?

Armed with your reprex, the next step is to figure out [where to ask](https://www.tidyverse.org/help/#where-to-ask). 

*   If it's a question: start with [community.rstudio.com](https://community.rstudio.com/), and/or StackOverflow. There are more people there to answer questions.  

*   If it's a bug: you're in the right place, [file an issue](https://github.com/batpigandme/glue/issues/new).  
  
*   If you're not sure: let the community help you figure it out! 
    If your problem _is_ a bug or a feature request, you can easily return here and report it. 

Before opening a new issue, be sure to [search issues and pull requests](https://github.com/batpigandme/glue/issues) to make sure the bug hasn't been reported and/or already fixed in the development version. 
By default, the search will be pre-populated with `is:issue is:open`. 
You can [edit the qualifiers](https://help.github.com/articles/searching-issues-and-pull-requests/)  (e.g. `is:pr`, `is:closed`) as needed. 
For example, you'd simply remove `is:open` to search _all_ issues in the repo, open or closed.

## What happens next?

To be as efficient as possible, development of tidyverse packages tends to be very bursty, so you shouldn't worry if you don't get an immediate response.
Typically we don't look at a repo until a sufficient quantity of issues accumulates, then there’s a burst of intense activity as we focus our efforts. 
That makes development more efficient because it avoids expensive context switching between problems, at the cost of taking longer to get back to you. 
This process makes a good reprex particularly important because it might be multiple months between your initial report and when we start working on it. 
If we can’t reproduce the bug, we can’t fix it!
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
# Contributing to glue

This outlines how to propose a change to glue. 
For more detailed info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib). 

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("batpigandme/glue", fork = TRUE)`.

*   Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://roxygen2.r-lib.org/articles/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  

## Code of Conduct

Please note that the glue project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
---
name: Bug report or feature request
about: Describe a bug you've seen or make a case for a new feature
---

Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, ask on <https://stackoverflow.com/> or <https://community.rstudio.com/>.

Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](http://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

Brief description of the problem

```r
# insert reprex here
```
## revdepcheck results

We checked 463 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 1 new problem
 * We failed to check 0 packages

Issues with CRAN packages are summarised below.

### New problems

* iotables

  As of 2022-01-21, I see very similar failure for iotables on CRAN, on
  flavours r-release-windows-ix86+x86_64 and
  r-oldrel-windows-ix86+x86_64.
  I suspect that the iotables package makes live API calls in its vignettes and
  is subject to intermittent failure.
  I don't think this has anything to do with glue.

  checking re-building of vignette outputs ... WARNING

  Error(s) in re-building vignettes:
  --- re-building ‘environmental_impact.Rmd’ using rmarkdown
  ...
  Quitting from lines 105-106 (environmental_impact.Rmd) 
  Error: processing vignette 'environmental_impact.Rmd' failed with diagnostics:
  get_eurostat_raw fails with the id env_ac_ainah_r2
  --- failed re-building ‘environmental_impact.Rmd’
  ...
  SUMMARY: processing the following file failed:
    ‘environmental_impact.Rmd’
  Error: Vignette re-building failed.
# Revdeps

## New problems (1)

|package                          |version |error |warning |note |
|:--------------------------------|:-------|:-----|:-------|:----|
|[iotables](problems.md#iotables) |0.4.7   |      |__+1__  |     |

# iotables

<details>

* Version: 0.4.7
* GitHub: https://github.com/rOpenGov/iotables
* Source code: https://github.com/cran/iotables
* Date/Publication: 2021-12-22 17:30:02 UTC
* Number of recursive dependencies: 115

Run `cloud_details(, "iotables")` for more info

</details>

## Newly broken

*   checking re-building of vignette outputs ... WARNING
    ```
    Error(s) in re-building vignettes:
    --- re-building ‘environmental_impact.Rmd’ using rmarkdown
    
    Attaching package: 'dplyr'
    
    The following objects are masked from 'package:stats':
    
        filter, lag
    
    The following objects are masked from 'package:base':
    ...
    Columns and rows of CPA_L68A, CPA_T, CPA_U are all zeros and will be removed.
    Joining, by = c("prod_na", "CPA_A01", "CPA_A02", "CPA_A03", "CPA_B", "CPA_C10-12", "CPA_C13-15", "CPA_C16", "CPA_C17", "CPA_C18", "CPA_C19", "CPA_C20", "CPA_C21", "CPA_C22", "CPA_C23", "CPA_C24", "CPA_C25", "CPA_C26", "CPA_C27", "CPA_C28", "CPA_C29", "CPA_C30", "CPA_C31_32", "CPA_C33", "CPA_D", "CPA_E36", "CPA_E37-39", "CPA_F", "CPA_G45", "CPA_G46", "CPA_H49", "CPA_H50", "CPA_H51", "CPA_H52", "CPA_H53", "CPA_I", "CPA_J58", "CPA_J59_60", "CPA_J61", "CPA_J62_63", "CPA_K64", "CPA_K65", "CPA_K66", "CPA_L68B", "CPA_M69_70", "CPA_M71", "CPA_M72", "CPA_M73", "CPA_M74_75", "CPA_N77", "CPA_N78", "CPA_N79", "CPA_N80-82", "CPA_O", "CPA_P", "CPA_Q86", "CPA_Q87_88", "CPA_R90-92", "CPA_R93", "CPA_S94", "CPA_S95", "CPA_S96")
    Columns and rows of CPA_L68A are all zeros and will be removed.
    --- finished re-building ‘working_with_eurostat.Rmd’
    
    SUMMARY: processing the following file failed:
      ‘environmental_impact.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

# Check times

|   |package         |version | check_time|
|:--|:---------------|:-------|----------:|
|3  |dplyr           |0.7.0   |      142.5|
|4  |shinycssloaders |0.2.0   |         37|
|1  |datadogr        |0.1.0   |       30.7|
|2  |dbplyr          |1.0.0   |        7.3|


*Wow, no problems at all. :)*---
output:
  github_document:
    html_preview: false
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
library(glue)
```

# glue <a href='https://glue.tidyverse.org'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/glue)](https://cran.r-project.org/package=glue)
[![R-CMD-check](https://github.com/tidyverse/glue/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tidyverse/glue/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/tidyverse/glue/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/tidyverse/glue/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

## Overview

Glue offers interpreted string literals that are small, fast, and dependency-free.  Glue does this by embedding R expressions in curly braces which are then evaluated and inserted into the argument string. 

## Installation

::: .pkgdown-release
```{r, eval = FALSE}
# Install released version from CRAN
install.packages("glue")
```
:::

::: .pkgdown-devel
```{r, eval = FALSE}
# Install development version from GitHub
devtools::install_github("tidyverse/glue")
```
:::

## Usage

##### Variables can be passed directly into strings.
```{r}
library(glue)
name <- "Fred"
glue('My name is {name}.')
```

Note that `glue::glue()` is also made available via `stringr::str_glue()`.
So if you've already attached stringr (or perhaps the whole tidyverse), you can access `glue()` like so:

```{r eval = FALSE}
library(stringr) # or library(tidyverse)

stringr_fcn <- "`stringr::str_glue()`"
glue_fcn    <- "`glue::glue()`"

str_glue('{stringr_fcn} is essentially an alias for {glue_fcn}.')
#> `stringr::str_glue()` is essentially an alias for `glue::glue()`.
```

##### Long strings are broken by line and concatenated together.
```{r}
library(glue)

name <- "Fred"
age <- 50
anniversary <- as.Date("1991-10-12")
glue('My name is {name},',
  ' my age next year is {age + 1},',
  ' my anniversary is {format(anniversary, "%A, %B %d, %Y")}.')
```

##### Named arguments are used to assign temporary variables.
```{r}
glue('My name is {name},',
  ' my age next year is {age + 1},',
  ' my anniversary is {format(anniversary, "%A, %B %d, %Y")}.',
  name = "Joe",
  age = 40,
  anniversary = as.Date("2001-10-12"))
```

##### `glue_data()` is useful with [magrittr](https://cran.r-project.org/package=magrittr) pipes.
```{r}
`%>%` <- magrittr::`%>%`
head(mtcars) %>% glue_data("{rownames(.)} has {hp} hp")
```

##### Or within dplyr pipelines
```{r, message = FALSE}
library(dplyr)
head(iris) %>%
  mutate(description = glue("This {Species} has a petal length of {Petal.Length}"))
```

##### Leading whitespace and blank lines from the first and last lines are automatically trimmed.
This lets you indent the strings naturally in code.
```{r}
glue("
    A formatted string
    Can have multiple lines
      with additional indention preserved
    ")
```

##### An additional newline can be used if you want a leading or trailing newline.
```{r}
glue("

  leading or trailing newlines can be added explicitly

  ")
```

##### `\\` at the end of a line continues it without a new line.
```{r}
glue("
    A formatted string \\
    can also be on a \\
    single line
    ")
```

##### A literal brace is inserted by using doubled braces.
```{r}
name <- "Fred"
glue("My name is {name}, not {{name}}.")
```

##### Alternative delimiters can be specified with `.open` and `.close`.
```{r}
one <- "1"
glue("The value of $e^{2\\pi i}$ is $<<one>>$.", .open = "<<", .close = ">>")
```

##### All valid R code works in expressions, including braces and escaping.
Backslashes do need to be doubled just like in all R strings.
```{r}
  `foo}\`` <- "foo"
glue("{
      {
        '}\\'' # { and } in comments, single quotes
        \"}\\\"\" # or double quotes are ignored
        `foo}\\`` # as are { in backticks
      }
  }")
```

##### `glue_sql()` makes constructing SQL statements safe and easy
Use backticks to quote identifiers, normal strings and numbers are quoted
appropriately for your backend.

```{r}
library(glue)

con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
colnames(iris) <- gsub("[.]", "_", tolower(colnames(iris)))
DBI::dbWriteTable(con, "iris", iris)
var <- "sepal_width"
tbl <- "iris"
num <- 2
val <- "setosa"
glue_sql("
  SELECT {`var`}
  FROM {`tbl`}
  WHERE {`tbl`}.sepal_length > {num}
    AND {`tbl`}.species = {val}
  ", .con = con)

# `glue_sql()` can be used in conjunction with parameterized queries using
# `DBI::dbBind()` to provide protection for SQL Injection attacks
 sql <- glue_sql("
    SELECT {`var`}
    FROM {`tbl`}
    WHERE {`tbl`}.sepal_length > ?
  ", .con = con)
query <- DBI::dbSendQuery(con, sql)
DBI::dbBind(query, list(num))
DBI::dbFetch(query, n = 4)
DBI::dbClearResult(query)

# `glue_sql()` can be used to build up more complex queries with
# interchangeable sub queries. It returns `DBI::SQL()` objects which are
# properly protected from quoting.
sub_query <- glue_sql("
  SELECT *
  FROM {`tbl`}
  ", .con = con)

glue_sql("
  SELECT s.{`var`}
  FROM ({sub_query}) AS s
  ", .con = con)

# If you want to input multiple values for use in SQL IN statements put `*`
# at the end of the value and the values will be collapsed and quoted appropriately.
glue_sql("SELECT * FROM {`tbl`} WHERE sepal_length IN ({vals*})",
  vals = 1, .con = con)

glue_sql("SELECT * FROM {`tbl`} WHERE sepal_length IN ({vals*})",
  vals = 1:5, .con = con)

glue_sql("SELECT * FROM {`tbl`} WHERE species IN ({vals*})",
  vals = "setosa", .con = con)

glue_sql("SELECT * FROM {`tbl`} WHERE species IN ({vals*})",
  vals = c("setosa", "versicolor"), .con = con)
```

##### Optionally combine strings with `+`

```{r}
x <- 1
y <- 3
glue("x + y") + " = {x + y}"
```

# Other implementations

Some other implementations of string interpolation in R (although not using identical syntax).

- [stringr::str_interp](https://stringr.tidyverse.org/reference/str_interp.html)
- [R.utils::gstring](https://cran.r-project.org/package=R.utils)
- [rprintf](https://cran.r-project.org/package=rprintf)

String templating is closely related to string interpolation, although not
exactly the same concept. Some packages implementing string templating in R
include.

- [whisker](https://cran.r-project.org/package=whisker)
- [brew](https://cran.r-project.org/package=brew)

## Code of Conduct

Please note that the glue project is released with a [Contributor Code of Conduct](https://glue.tidyverse.org/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
---
title: "glue custom knitr language engines"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{glue custom knitr language engines}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(glue)
```

Glue provides a few [custom language engines](https://bookdown.org/yihui/rmarkdown-cookbook/custom-engine.html#custom-engine) for knitr, which allows you to use glue directly in knitr chunks.

## `glue` engine

The first engine is the `glue` engine, which evaluates the chunk contents as a glue template.

```{glue}
1 + 1 = {1 + 1}
```

Maybe the most useful use of the `glue` engine is to set the knitr option `results = 'asis'` and output markdown or HTML directly into the document.

````markdown
`r '' ````{glue, results = 'asis', echo = FALSE}
#### mtcars has **{nrow(mtcars)} rows** and _{ncol(mtcars)} columns_.
```
````

```{glue, results = 'asis', echo = FALSE}
#### mtcars has **{nrow(mtcars)} rows** and _{ncol(mtcars)} columns_.
```

If you want to pass additional arguments into the glue call, simply include them as chunk options.

````markdown
`r '' ````{glue, .open = "<<", .close = ">>", results = 'asis', echo = FALSE}
The **median waiting time** between eruptions is <<median(faithful$waiting)>>.
```
````

```{glue, .open = "<<", .close = ">>", results = 'asis', echo = FALSE}
The **median waiting time** between eruptions is <<median(faithful$waiting)>>.
```

## `glue_sql` engine

The second engine is `glue_sql`, which will use `glue::glue_sql()` to generate a SQL query and then run the query using the [sql engine](https://bookdown.org/yihui/rmarkdown/language-engines.html#sql).

First we create a new connection to an in-memory SQLite database, and write a new table to it.
```{r}
con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
mtcars$model <- rownames(mtcars)
DBI::dbWriteTable(con, "mtcars", mtcars)
```

Next define some variables we that we can use with glue to interpolate.
```{r}
var <- "mpg"
tbl <- "mtcars"
num <- 150
```

Then we can use `glue_sql` to construct and run a query using those variables into that database. *Note* you need to provide the connection object as a `connection` chunk option.

In this example there are two type of quotes. The first is a bare backtick, these are passed directly to the SQL engine unchanged.
The second is backticks inside of braces, which are specially interpreted to do the proper quoting for the given SQL engine by glue.
In this example we use the `sqlite` engine, which uses backticks for quoting, but you would use the same backticks inside brace syntax for postgreSQL, and `glue_sql()` would automatically use double quotes for quoting instead.

````markdown
`r '' ````{glue_sql, connection = con}
SELECT `model`, `hp`, {`var`}
FROM {`tbl`}
WHERE {`tbl`}.hp > {num}
```
````

```{glue_sql, connection = con}
SELECT `model`, `hp`, {`var`}
  FROM {`tbl`}
  WHERE {`tbl`}.hp > {num}
```
---
title: "Speed of glue"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Speed of glue}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{R.utils
    R.utils,
    forcats,
    microbenchmark,
    rprintf,
    stringr,
    ggplot2}
---

Glue is advertised as

> Fast, dependency free string literals

So what do we mean when we say that glue is fast? This does not mean glue is
the fastest thing to use in all cases, however for the features it provides we
can confidently say it is fast.

A good way to determine this is to compare it's speed of execution to some alternatives.

- `base::paste0()`, `base::sprintf()` - Functions in base R implemented in C
  that provide variable insertion (but not interpolation).
- `R.utils::gstring()`, `stringr::str_interp()` - Provides a similar interface
  as glue, but using `${}` to delimit blocks to interpolate.
- `pystr::pystr_format()`[^1], `rprintf::rprintf()` - Provide a interfaces similar
  to python string formatters with variable replacement, but not arbitrary
  interpolation.

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  eval = as.logical(Sys.getenv("EVAL_VIGNETTES", "FALSE")),
  cache = FALSE)
library(glue)
```

```{r setup2, include = FALSE}
plot_comparison <- function(x, ...) {
  library(ggplot2)
  library(microbenchmark)
  x$expr <- forcats::fct_reorder(x$expr, x$time)
  colors <- ifelse(levels(x$expr) == "glue", "orange", "grey")
  autoplot(x, ...) +
    theme(axis.text.y = element_text(color = colors)) +
      aes(fill = expr) + scale_fill_manual(values = colors, guide = FALSE)
}
```

## Simple concatenation

```{r, message = FALSE}
bar <- "baz"

simple <-
  microbenchmark::microbenchmark(
  glue = glue::glue("foo{bar}"),
  gstring = R.utils::gstring("foo${bar}"),
  paste0 = paste0("foo", bar),
  sprintf = sprintf("foo%s", bar),
  str_interp = stringr::str_interp("foo${bar}"),
  rprintf = rprintf::rprintf("foo$bar", bar = bar)
)

print(unit = "eps", order = "median", signif = 4, simple)

plot_comparison(simple)
```

While `glue()` is slower than `paste0`,`sprintf()` it is
twice as fast as `str_interp()` and `gstring()`, and on par with `rprintf()`.

Although `paste0()`, `sprintf()` don't do string interpolation and will likely always be
significantly faster than glue, glue was never meant to be a direct replacement
for them.

`rprintf()` does only variable interpolation, not arbitrary expressions, which
was one of the explicit goals of writing glue.

So glue is ~2x as fast as the two functions (`str_interp()`, `gstring()`), which do have
roughly equivalent functionality.

It also is still quite fast, with over 6000 evaluations per second on this machine.

## Vectorized performance

Taking advantage of glue's vectorization is the best way to avoid performance.
For instance the vectorized form of the previous benchmark is able to generate
100,000 strings in only 22ms with performance much closer to that of
`paste0()` and `sprintf()`. NB: `str_interp()` does not support
vectorization, and so was removed.

```{r, message = FALSE}
bar <- rep("bar", 1e5)

vectorized <-
  microbenchmark::microbenchmark(
  glue = glue::glue("foo{bar}"),
  gstring = R.utils::gstring("foo${bar}"),
  paste0 = paste0("foo", bar),
  sprintf = sprintf("foo%s", bar),
  rprintf = rprintf::rprintf("foo$bar", bar = bar)
)

print(unit = "ms", order = "median", signif = 4, vectorized)

plot_comparison(vectorized, log = FALSE)
```

[^1]: pystr is no longer available from CRAN due to failure to correct
installation errors and was therefore removed from further testing.
---
title: "Transformers"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Transformers}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Transformers allow you to apply functions to the glue input and output, before
and after evaluation. This allows you to write things like `glue_sql()`, which
automatically quotes variables for you or add a syntax for automatically
collapsing outputs.

The transformer functions simply take two arguments `text` and `envir`, where
`text` is the unparsed string inside the glue block and `envir` is the
execution environment. Most transformers will then call `eval(parse(text = text,
keep.source = FALSE), envir)` which parses and evaluates the code.

You can then supply the transformer function to glue with the `.transformer`
argument. In this way users can manipulate the text before parsing and
change the output after evaluation.

It is often useful to write a `glue()` wrapper function which supplies a
`.transformer` to `glue()` or `glue_data()` and potentially has additional
arguments. One important consideration when doing this is to include
`.envir = parent.frame()` in the wrapper to ensure the evaluation environment
is correct.

Some example implementations of potentially useful transformers follow. The
aim right now is not to include most of these custom functions within the
`glue` package. Rather, users are encouraged to create custom functions using
transformers to fit their individual needs.

```{r, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

```{r}
library(glue)
```

### collapse transformer

A transformer which automatically collapses any glue block ending with `*`.

```{r}
collapse_transformer <- function(regex = "[*]$", ...) {
  function(text, envir) {
    collapse <- grepl(regex, text)
    if (collapse) {
      text <- sub(regex, "", text)
    }
    res <- identity_transformer(text, envir)
    if (collapse) {
      glue_collapse(res, ...)  
    } else {
      res
    }
  }
}

glue("{1:5*}\n{letters[1:5]*}", .transformer = collapse_transformer(sep = ", "))

glue("{1:5*}\n{letters[1:5]*}", .transformer = collapse_transformer(sep = ", ", last = " and "))

x <- c("one", "two")
glue("{x}: {1:5*}", .transformer = collapse_transformer(sep = ", "))
```

### Shell quoting transformer

A transformer which automatically quotes variables for use in shell commands,
e.g. via `system()` or `system2()`.

```{r}
shell_transformer <- function(type = c("sh", "csh", "cmd", "cmd2")) {
  type <- match.arg(type)
  function(text, envir) {
    res <- eval(parse(text = text, keep.source = FALSE), envir)
    shQuote(res)
  }
}

glue_sh <- function(..., .envir = parent.frame(), .type = c("sh", "csh", "cmd", "cmd2")) {
  .type <- match.arg(.type)
  glue(..., .envir = .envir, .transformer = shell_transformer(.type))

}

filename <- "test"
writeLines(con = filename, "hello!")

command <- glue_sh("cat {filename}")
command
system(command)
```

```{r include = FALSE}
if (file.exists("test")) {
  unlink("test")
}
```

### emoji transformer

A transformer which converts the text to the equivalent emoji.

```{r, eval = require("emo")}
emoji_transformer <- function(text, envir) {
  if (grepl("[*]$", text)) {
    text <- sub("[*]$", "", text)
    glue_collapse(ji_find(text)$emoji)
  } else {
    ji(text)
  }
}

glue_ji <- function(..., .envir = parent.frame()) {
  glue(..., .open = ":", .close = ":", .envir = .envir, .transformer = emoji_transformer)
}
glue_ji("one :heart:")
glue_ji("many :heart*:")
```

### sprintf transformer

A transformer which allows succinct sprintf format strings.

```{r}
sprintf_transformer <- function(text, envir) {
  m <- regexpr(":.+$", text)
  if (m != -1) {
    format <- substring(regmatches(text, m), 2)
    regmatches(text, m) <- ""
    res <- eval(parse(text = text, keep.source = FALSE), envir)
    do.call(sprintf, list(glue("%{format}"), res))
  } else {
    eval(parse(text = text, keep.source = FALSE), envir)
  }
}

glue_fmt <- function(..., .envir = parent.frame()) {
  glue(..., .transformer = sprintf_transformer, .envir = .envir)
}
glue_fmt("π = {pi:.3f}")
```

### safely transformer

A transformer that acts like `purrr::safely()`, which returns a value instead of an error.

```{r}
safely_transformer <- function(otherwise = NA) {
  function(text, envir) {
    tryCatch(
      eval(parse(text = text, keep.source = FALSE), envir),
      error = function(e) if (is.language(otherwise)) eval(otherwise) else otherwise)
  }
}

glue_safely <- function(..., .otherwise = NA, .envir = parent.frame()) {
  glue(..., .transformer = safely_transformer(.otherwise), .envir = .envir)
}

# Default returns missing if there is an error
glue_safely("foo: {xyz}")

# Or an empty string
glue_safely("foo: {xyz}", .otherwise = "Error")

# Or output the error message in red
library(crayon)
glue_safely("foo: {xyz}", .otherwise = quote(glue("{red}Error: {conditionMessage(e)}{reset}")))
```

### "Variables and Values" transformer

A transformer that expands input of the form `{var_name=}` into `var_name = var_value`, i.e. a shorthand for exposing variable names with their values. This is inspired by an [f-strings feature coming in Python 3.8](https://docs.python.org/3.8/whatsnew/3.8.html#f-strings-now-support-for-quick-and-easy-debugging). It's actually more general: you can use it with an expression input such as `{expr=}`.

```{r}
vv_transformer <- function(text, envir) {
  regex <- "=$"
  if (!grepl(regex, text)) {
    return(identity_transformer(text, envir))
  }

  text <- sub(regex, "", text)
  res <- identity_transformer(text, envir)
  n <- length(res)
  res <- glue_collapse(res, sep = ", ")
  if (n > 1) {
    res <- c("[", res, "]")
  }
  glue_collapse(c(text, " = ", res))
}
```

```{r}
set.seed(1234)
description <- "some random"
numbers <- sample(100, 4)
average <- mean(numbers)
sum <- sum(numbers)

glue("For {description} {numbers=}, {average=}, {sum=}.", .transformer = vv_transformer)

a <- 3
b <- 5.6
glue("{a=}\n{b=}\n{a * 9 + b * 2=}", .transformer = vv_transformer)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/safe.R
\name{glue_safe}
\alias{glue_safe}
\alias{glue_data_safe}
\title{Safely interpolate strings}
\usage{
glue_safe(..., .envir = parent.frame())

glue_data_safe(.x, ..., .envir = parent.frame())
}
\arguments{
\item{...}{[\code{expressions}]\cr Unnamed arguments are taken to be expression
string(s) to format. Multiple inputs are concatenated together before formatting.
Named arguments are taken to be temporary variables available for substitution.}

\item{.envir}{[\code{environment}: \code{parent.frame()}]\cr Environment to evaluate each expression in. Expressions are
evaluated from left to right. If \code{.x} is an environment, the expressions are
evaluated in that environment and \code{.envir} is ignored. If \code{NULL} is passed, it is equivalent to \code{\link[=emptyenv]{emptyenv()}}.}

\item{.x}{[\code{listish}]\cr An environment, list, or data frame used to lookup values.}
}
\description{
\code{glue_safe()} and \code{glue_data_safe()} differ from \code{\link[=glue]{glue()}} and \code{\link[=glue_data]{glue_data()}}
in that the safe versions only look up symbols from an environment using
\code{\link[=get]{get()}}. They do not execute any R code. This makes them suitable for use
with untrusted input, such as inputs in a Shiny application, where using the
normal functions would allow an attacker to execute arbitrary code.
}
\examples{
"1 + 1" <- 5
# glue actually executes the code
glue("{1 + 1}")

# glue_safe just looks up the value
glue_safe("{1 + 1}")

rm("1 + 1")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformer.R
\name{identity_transformer}
\alias{identity_transformer}
\title{Parse and Evaluate R code}
\usage{
identity_transformer(text, envir)
}
\arguments{
\item{text}{Text (typically) R code to parse and evaluate.}

\item{envir}{environment to evaluate the code in}
}
\description{
This is a simple wrapper around \code{eval(parse())}, used as the default
transformer.
}
\seealso{
\code{vignette("transformers", "glue")} for documentation on creating
custom glue transformers and some common use cases.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glue.R, R/sql.R
\name{glue_collapse}
\alias{glue_collapse}
\alias{glue_sql_collapse}
\title{Collapse a character vector}
\usage{
glue_collapse(x, sep = "", width = Inf, last = "")

glue_sql_collapse(x, sep = "", width = Inf, last = "")
}
\arguments{
\item{x}{The character vector to collapse.}

\item{sep}{a character string to separate the terms.  Not
    \code{\link[base]{NA_character_}}.}

\item{width}{The maximum string width before truncating with \code{...}.}

\item{last}{String used to separate the last two items if \code{x} has at least
2 items.}
}
\description{
\code{glue_collapse()} collapses a character vector of any length into a length 1 vector.
\code{glue_sql_collapse()} does the same but returns a \verb{[DBI::SQL()]}
object rather than a glue object.
}
\examples{
glue_collapse(glue("{1:10}"))

# Wide values can be truncated
glue_collapse(glue("{1:10}"), width = 5)

glue_collapse(1:4, ", ", last = " and ")
#> 1, 2, 3 and 4
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glue.R
\name{trim}
\alias{trim}
\title{Trim a character vector}
\usage{
trim(x)
}
\arguments{
\item{x}{A character vector to trim.}
}
\description{
This trims a character vector according to the trimming rules used by glue.
These follow similar rules to \href{https://www.python.org/dev/peps/pep-0257/}{Python Docstrings},
with the following features.
\itemize{
\item Leading and trailing whitespace from the first and last lines is removed.
\item A uniform amount of indentation is stripped from the second line on, equal
to the minimum indentation of all non-blank lines after the first.
\item Lines can be continued across newlines by using \verb{\\\\}.
}
}
\examples{
glue("
    A formatted string
    Can have multiple lines
      with additional indention preserved
    ")

glue("
  \ntrailing or leading newlines can be added explicitly\n
  ")

glue("
    A formatted string \\\\
    can also be on a \\\\
    single line
    ")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glue.R
\name{glue-deprecated}
\alias{glue-deprecated}
\title{Deprecated Functions}
\description{
These functions are Deprecated in this release of glue, they will be removed
in a future version.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color.R
\name{glue_col}
\alias{glue_col}
\alias{glue_data_col}
\title{Construct strings with color}
\usage{
glue_col(..., .envir = parent.frame(), .na = "NA", .literal = FALSE)

glue_data_col(.x, ..., .envir = parent.frame(), .na = "NA", .literal = FALSE)
}
\arguments{
\item{...}{[\code{expressions}]\cr Unnamed arguments are taken to be expression
string(s) to format. Multiple inputs are concatenated together before formatting.
Named arguments are taken to be temporary variables available for substitution.}

\item{.envir}{[\code{environment}: \code{parent.frame()}]\cr Environment to evaluate each expression in. Expressions are
evaluated from left to right. If \code{.x} is an environment, the expressions are
evaluated in that environment and \code{.envir} is ignored. If \code{NULL} is passed, it is equivalent to \code{\link[=emptyenv]{emptyenv()}}.}

\item{.na}{[\code{character(1)}: \sQuote{NA}]\cr Value to replace \code{NA} values
with. If \code{NULL} missing values are propagated, that is an \code{NA} result will
cause \code{NA} output. Otherwise the value is replaced by the value of \code{.na}.}

\item{.literal}{[\code{boolean(1)}: \sQuote{FALSE}]\cr Whether to treat single or
double quotes, backticks, and comments as regular characters (vs. as
syntactic elements), when parsing the expression string. Setting \code{.literal = TRUE} probably only makes sense in combination with a custom
\code{.transformer}, as is the case with \code{glue_col()}. Regard this argument
(especially, its name) as experimental.}

\item{.x}{[\code{listish}]\cr An environment, list, or data frame used to lookup values.}
}
\description{
The \link[crayon:crayon]{crayon} package defines a number of functions used to
color terminal output. \code{glue_col()} and \code{glue_data_col()} functions provide
additional syntax to make using these functions in glue strings easier.

Using the following syntax will apply the function \code{\link[crayon:crayon]{crayon::blue()}} to the text 'foo bar'.\preformatted{\{blue foo bar\}
}

If you want an expression to be evaluated, simply place that in a normal brace
expression (these can be nested).\preformatted{\{blue 1 + 1 = \{1 + 1\}\}
}

If the text you want to color contains, e.g., an unpaired quote or a comment
character, specify \code{.literal = TRUE}.
}
\examples{
\dontshow{if (require(crayon)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library(crayon)

glue_col("{blue foo bar}")

glue_col("{blue 1 + 1 = {1 + 1}}")

glue_col("{blue 2 + 2 = {green {2 + 2}}}")

white_on_black <- bgBlack $ white
glue_col("{white_on_black
  Roses are {red {colors()[[552]]}},
  Violets are {blue {colors()[[26]]}},
  `glue_col()` can show \\\\
  {red c}{yellow o}{green l}{cyan o}{blue r}{magenta s}
  and {bold bold} and {underline underline} too!
}")

# this would error due to an unterminated quote, if we did not specify
# `.literal = TRUE`
glue_col("{yellow It's} happening!", .literal = TRUE)

# `.literal = TRUE` also prevents an error here due to the `#` comment
glue_col(
  "A URL: {magenta https://github.com/tidyverse/glue#readme}",
  .literal = TRUE
)

# `.literal = TRUE` does NOT prevent evaluation
x <- "world"
y <- "day"
glue_col("hello {x}! {green it's a new {y}!}", .literal = TRUE)
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glue.R
\name{glue}
\alias{glue}
\alias{glue_data}
\title{Format and interpolate a string}
\usage{
glue_data(
  .x,
  ...,
  .sep = "",
  .envir = parent.frame(),
  .open = "{",
  .close = "}",
  .na = "NA",
  .null = character(),
  .comment = "#",
  .literal = FALSE,
  .transformer = identity_transformer,
  .trim = TRUE
)

glue(
  ...,
  .sep = "",
  .envir = parent.frame(),
  .open = "{",
  .close = "}",
  .na = "NA",
  .null = character(),
  .comment = "#",
  .literal = FALSE,
  .transformer = identity_transformer,
  .trim = TRUE
)
}
\arguments{
\item{.x}{[\code{listish}]\cr An environment, list, or data frame used to lookup values.}

\item{...}{[\code{expressions}]\cr Unnamed arguments are taken to be expression
string(s) to format. Multiple inputs are concatenated together before formatting.
Named arguments are taken to be temporary variables available for substitution.}

\item{.sep}{[\code{character(1)}: \sQuote{""}]\cr Separator used to separate elements.}

\item{.envir}{[\code{environment}: \code{parent.frame()}]\cr Environment to evaluate each expression in. Expressions are
evaluated from left to right. If \code{.x} is an environment, the expressions are
evaluated in that environment and \code{.envir} is ignored. If \code{NULL} is passed, it is equivalent to \code{\link[=emptyenv]{emptyenv()}}.}

\item{.open}{[\code{character(1)}: \sQuote{\\\{}]\cr The opening delimiter. Doubling the
full delimiter escapes it.}

\item{.close}{[\code{character(1)}: \sQuote{\\\}}]\cr The closing delimiter. Doubling the
full delimiter escapes it.}

\item{.na}{[\code{character(1)}: \sQuote{NA}]\cr Value to replace \code{NA} values
with. If \code{NULL} missing values are propagated, that is an \code{NA} result will
cause \code{NA} output. Otherwise the value is replaced by the value of \code{.na}.}

\item{.null}{[\code{character(1)}: \sQuote{character()}]\cr Value to replace
NULL values with. If \code{character()} whole output is \code{character()}. If
\code{NULL} all NULL values are dropped (as in \code{paste0()}). Otherwise the
value is replaced by the value of \code{.null}.}

\item{.comment}{[\code{character(1)}: \sQuote{#}]\cr Value to use as the comment
character.}

\item{.literal}{[\code{boolean(1)}: \sQuote{FALSE}]\cr Whether to treat single or
double quotes, backticks, and comments as regular characters (vs. as
syntactic elements), when parsing the expression string. Setting \code{.literal = TRUE} probably only makes sense in combination with a custom
\code{.transformer}, as is the case with \code{glue_col()}. Regard this argument
(especially, its name) as experimental.}

\item{.transformer}{[\verb{function]}\cr A function taking three parameters \code{code}, \code{envir} and
\code{data} used to transform the output of each block before, during, or after
evaluation. For example transformers see \code{vignette("transformers")}.}

\item{.trim}{[\code{logical(1)}: \sQuote{TRUE}]\cr Whether to trim the input
template with \code{\link[=trim]{trim()}} or not.}
}
\description{
Expressions enclosed by braces will be evaluated as R code. Long strings are
broken by line and concatenated together. Leading whitespace and blank lines
from the first and last lines are automatically trimmed.
}
\examples{
name <- "Fred"
age <- 50
anniversary <- as.Date("1991-10-12")
glue('My name is {name},',
  'my age next year is {age + 1},',
  'my anniversary is {format(anniversary, "\%A, \%B \%d, \%Y")}.')

# single braces can be inserted by doubling them
glue("My name is {name}, not {{name}}.")

# Named arguments can be used to assign temporary variables.
glue('My name is {name},',
  ' my age next year is {age + 1},',
  ' my anniversary is {format(anniversary, "\%A, \%B \%d, \%Y")}.',
  name = "Joe",
  age = 40,
  anniversary = as.Date("2001-10-12"))

# `glue()` can also be used in user defined functions
intro <- function(name, profession, country){
  glue("My name is {name}, a {profession}, from {country}")
}
intro("Shelmith", "Senior Data Analyst", "Kenya")
intro("Cate", "Data Scientist", "Kenya")

# `glue_data()` is useful in magrittr pipes
if (require(magrittr)) {

mtcars \%>\% glue_data("{rownames(.)} has {hp} hp")

# Or within dplyr pipelines
if (require(dplyr)) {

head(iris) \%>\%
  mutate(description = glue("This {Species} has a petal length of {Petal.Length}"))

}}

# Alternative delimiters can also be used if needed
one <- "1"
glue("The value of $e^{2\\\\pi i}$ is $<<one>>$.", .open = "<<", .close = ">>")
}
\seealso{
\url{https://www.python.org/dev/peps/pep-0498/} and
\url{https://www.python.org/dev/peps/pep-0257/} upon which this is based.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sql.R
\name{glue_sql}
\alias{glue_sql}
\alias{glue_data_sql}
\title{Interpolate strings with SQL escaping}
\usage{
glue_sql(..., .con, .envir = parent.frame(), .na = DBI::SQL("NULL"))

glue_data_sql(.x, ..., .con, .envir = parent.frame(), .na = DBI::SQL("NULL"))
}
\arguments{
\item{...}{[\code{expressions}]\cr Unnamed arguments are taken to be expression
string(s) to format. Multiple inputs are concatenated together before formatting.
Named arguments are taken to be temporary variables available for substitution.}

\item{.con}{[\code{DBIConnection}]: A DBI connection object obtained from
\code{\link[DBI:dbConnect]{DBI::dbConnect()}}.}

\item{.envir}{[\code{environment}: \code{parent.frame()}]\cr Environment to evaluate each expression in. Expressions are
evaluated from left to right. If \code{.x} is an environment, the expressions are
evaluated in that environment and \code{.envir} is ignored. If \code{NULL} is passed, it is equivalent to \code{\link[=emptyenv]{emptyenv()}}.}

\item{.na}{[\code{character(1)}: \sQuote{NA}]\cr Value to replace \code{NA} values
with. If \code{NULL} missing values are propagated, that is an \code{NA} result will
cause \code{NA} output. Otherwise the value is replaced by the value of \code{.na}.}

\item{.x}{[\code{listish}]\cr An environment, list, or data frame used to lookup values.}
}
\value{
A \code{\link[DBI:SQL]{DBI::SQL()}} object with the given query.
}
\description{
SQL databases often have custom quotation syntax for identifiers and strings
which make writing SQL queries error prone and cumbersome to do. \code{glue_sql()} and
\code{glue_data_sql()} are analogs to \code{\link[=glue]{glue()}} and \code{\link[=glue_data]{glue_data()}} which handle the
SQL quoting. \code{glue_sql_collapse()} can be used to collapse \code{\link[DBI:SQL]{DBI::SQL()}} objects.
}
\details{
They automatically quote character results, quote identifiers if the glue
expression is surrounded by backticks '\verb{`}' and do not quote
non-characters such as numbers. If numeric data is stored in a character
column (which should be quoted) pass the data to \code{glue_sql()} as a
character.

Returning the result with \code{\link[DBI:SQL]{DBI::SQL()}} will suppress quoting if desired for
a given value.

Note \href{https://db.rstudio.com/best-practices/run-queries-safely#parameterized-queries}{parameterized queries}
are generally the safest and most efficient way to pass user defined
values in a query, however not every database driver supports them.

If you place a \code{*} at the end of a glue expression the values will be
collapsed with commas. This is useful for the \href{https://www.w3schools.com/sql/sql_in.asp}{SQL IN Operator}
for instance.
}
\examples{
\dontshow{if (requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
iris2 <- iris
colnames(iris2) <- gsub("[.]", "_", tolower(colnames(iris)))
DBI::dbWriteTable(con, "iris", iris2)
var <- "sepal_width"
tbl <- "iris"
num <- 2
val <- "setosa"
glue_sql("
  SELECT {`var`}
  FROM {`tbl`}
  WHERE {`tbl`}.sepal_length > {num}
    AND {`tbl`}.species = {val}
  ", .con = con)

# If sepal_length is store on the database as a character explicitly convert
# the data to character to quote appropriately.
glue_sql("
  SELECT {`var`}
  FROM {`tbl`}
  WHERE {`tbl`}.sepal_length > {as.character(num)}
    AND {`tbl`}.species = {val}
  ", .con = con)


# `glue_sql()` can be used in conjuction with parameterized queries using
# `DBI::dbBind()` to provide protection for SQL Injection attacks
 sql <- glue_sql("
    SELECT {`var`}
    FROM {`tbl`}
    WHERE {`tbl`}.sepal_length > ?
  ", .con = con)
query <- DBI::dbSendQuery(con, sql)
DBI::dbBind(query, list(num))
DBI::dbFetch(query, n = 4)
DBI::dbClearResult(query)

# `glue_sql()` can be used to build up more complex queries with
# interchangeable sub queries. It returns `DBI::SQL()` objects which are
# properly protected from quoting.
sub_query <- glue_sql("
  SELECT *
  FROM {`tbl`}
  ", .con = con)

glue_sql("
  SELECT s.{`var`}
  FROM ({sub_query}) AS s
  ", .con = con)

# If you want to input multiple values for use in SQL IN statements put `*`
# at the end of the value and the values will be collapsed and quoted appropriately.
glue_sql("SELECT * FROM {`tbl`} WHERE sepal_length IN ({vals*})",
  vals = 1, .con = con)

glue_sql("SELECT * FROM {`tbl`} WHERE sepal_length IN ({vals*})",
  vals = 1:5, .con = con)

glue_sql("SELECT * FROM {`tbl`} WHERE species IN ({vals*})",
  vals = "setosa", .con = con)

glue_sql("SELECT * FROM {`tbl`} WHERE species IN ({vals*})",
  vals = c("setosa", "versicolor"), .con = con)

# If you need to reference variables from multiple tables use `DBI::Id()`.
# Here we create a new table of nicknames, join the two tables together and
# select columns from both tables. Using `DBI::Id()` and the special
# `glue_sql()` syntax ensures all the table and column identifiers are quoted
# appropriately.

iris_db <- "iris"
nicknames_db <- "nicknames"

nicknames <- data.frame(
  species = c("setosa", "versicolor", "virginica"),
  nickname = c("Beachhead Iris", "Harlequin Blueflag", "Virginia Iris"),
  stringsAsFactors = FALSE
)

DBI::dbWriteTable(con, nicknames_db, nicknames)

cols <- list(
  DBI::Id(table = iris_db, column = "sepal_length"),
  DBI::Id(table = iris_db, column = "sepal_width"),
  DBI::Id(table = nicknames_db, column = "nickname")
)

iris_species <- DBI::Id(table = iris_db, column = "species")
nicknames_species <- DBI::Id(table = nicknames_db, column = "species")

query <- glue_sql("
  SELECT {`cols`*}
  FROM {`iris_db`}
  JOIN {`nicknames_db`}
  ON {`iris_species`}={`nicknames_species`}",
  .con = con
)
query

DBI::dbGetQuery(con, query, n = 5)

DBI::dbDisconnect(con)
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link[=glue_sql_collapse]{glue_sql_collapse()}} to collapse \code{\link[DBI:SQL]{DBI::SQL()}} objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glue.R
\name{as_glue}
\alias{as_glue}
\title{Coerce object to glue}
\usage{
as_glue(x, ...)
}
\arguments{
\item{x}{object to be coerced.}

\item{...}{further arguments passed to methods.}
}
\description{
Coerce object to glue
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quoting.R
\name{quoting}
\alias{quoting}
\alias{single_quote}
\alias{double_quote}
\alias{backtick}
\title{Quoting operators}
\usage{
single_quote(x)

double_quote(x)

backtick(x)
}
\arguments{
\item{x}{A character to quote.}
}
\description{
These functions make it easy to quote each individual element and are useful
in conjunction with \code{\link[=glue_collapse]{glue_collapse()}}.
}
\examples{
x <- 1:5
glue('Values of x: {glue_collapse(backtick(x), sep = ", ", last = " and ")}')
}
