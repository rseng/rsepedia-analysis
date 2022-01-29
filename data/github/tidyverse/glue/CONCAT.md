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


*Wow, no problems at all. :)*